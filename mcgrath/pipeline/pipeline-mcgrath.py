#################################################################
#################################################################
############### McGrath Microarray Analysis ################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
from ruffus import *
import sys, os
import pandas as pd
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
import Support3 as S
import Mcgrath as P

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####

##### 2. R Connection #####
r.source('pipeline/scripts/mcgrath.R')

#######################################################
#######################################################
########## S1. Prepare Dataset
#######################################################
#######################################################

#############################################
########## 1. Read IDAT
#############################################

@follows(mkdir('s1-samples.dir'))

@transform('rawdata/*/*.idat',
           regex(r'(.*)/(.*).idat'),
           add_inputs(r'\1/HumanHT-12_V4_0_R2_15002873_B.bgx'),
           r's1-samples.dir/\2.txt')

def readIdat(infiles, outfile):

    # Read expression
    r.read_idat(infiles[0], infiles[1], outfile)

#############################################
########## 2. Merge
#############################################

@follows(mkdir('s2-expression.dir'))

@merge([readIdat, 'rawdata/conditions_correct.csv'],
       's2-expression.dir/mcgrath-rawdata.txt')

def mergeData(infiles, outfile):

    # Read metadata
    sample_metadata_dataframe = pd.read_csv(infiles.pop())

    # Initialize results
    results = []

    # Loop through infiles
    for infile in infiles:

        # Read dataframe
        sample_dataframe = pd.read_table(infile)

        # Add sample name
        sample_dataframe['sample_name'] = os.path.basename(infile)[:-len(".txt")]

        # Append
        results.append(sample_dataframe)

    # Rename
    rename_dict = {rowData['IDATfile'][:-len('.idat')]: rowData['Samples'] for index, rowData in sample_metadata_dataframe.iterrows()}

    # Concatenate and cast
    expression_dataframe = pd.concat(results).pivot_table(index='Symbol', columns='sample_name', values='expression').rename(columns=rename_dict)

    # Write
    expression_dataframe.to_csv(outfile, sep='\t')

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################


##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')
