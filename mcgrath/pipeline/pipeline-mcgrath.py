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
import sys, os, itertools
import pandas as pd
from functools import reduce
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
# Id columns
id_cols = ['Probe_Id', 'Array_Address_Id', 'Status', 'Symbol']
sample_metadata = 'rawdata/conditions_correct.csv'

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

@transform('rawdata/7196780027/*.idat',
           regex(r'(.*)/(.*).idat'),
           add_inputs(r'\1/HumanHT-12_V4_0_R2_15002873_B.bgx'),
           r's1-expression.dir/\2.txt')

def readIdat(infiles, outfile):

    # Read expression
    print(infiles, outfile)
#     r.read_idat(infiles[0], infiles[1], outfile)

#############################################
########## 2. Merge
#############################################

@follows(mkdir('s2-expression.dir'))

@merge([readIdat, sample_metadata],
       's2-expression.dir/mcgrath-rawdata.txt')

def mergeData(infiles, outfile):

    # Read metadata
    sample_metadata_dataframe = pd.read_csv(infiles.pop())

    # Initialize results
    dataframes = []

    # Loop through infiles 
    for infile in infiles:

        # Read dataframe
        sample_dataframe = pd.read_table(infile)[['Symbol', 'expression']]

        # Add sample
        sample_dataframe['sample'] = os.path.basename(infile)[:-len(".txt")]

        # Append
        dataframes.append(sample_dataframe)

    # Rename
    rename_dict = {rowData['IDATfile'][:-len('.idat')]: rowData['Samples'] for index, rowData in sample_metadata_dataframe.iterrows()}

    # Merge
    # merged_dataframe = reduce(lambda x, y: pd.merge(x, y, on='Symbol'), dataframes).rename(columns=rename_dict)[id_cols+list(rename_dict.values())]
    merged_dataframe = pd.concat(dataframes).pivot_table(index='Symbol', columns='sample', values='expression').rename(columns=rename_dict)[list(rename_dict.values())]

    # Write
    merged_dataframe.to_csv(outfile, sep='\t')

#############################################
########## 3. Normalize
#############################################

@transform(mergeData,
           suffix('rawdata.txt'),
           'quantile.txt')

def normalizeData(infile, outfile):

    # Read metadata
    r.quantile_normalization(infile, outfile)

#######################################################
#######################################################
########## S2. Differential Expression
#######################################################
#######################################################

#############################################
########## 1. Run limma
#############################################

@follows(mkdir('s3-differential_expression.dir'))

@subdivide(normalizeData,
           formatter(),
           add_inputs(sample_metadata),
           's3-differential_expression.dir/*-limma.txt',
           's3-differential_expression.dir/')

def runLimma(infiles, outfiles, outfileRoot):

    # Get comparisons
    comparisons = (['Control', 'Wounded'], ['Control', 'Normal'], ['Normal', 'Wounded'])

    # Loop through comparisons
    for comparison in comparisons:

        # Outfile
        outfile = '{outfileRoot}{comparison[0]}_vs_{comparison[1]}-limma.txt'.format(**locals())

        # Run Limma
        r.run_limma(expression_file=infiles[0], metadata_file=infiles[1], outfile=outfile, comparison=comparison)


##################################################
##################################################
########## Run pipeline
##################################################
##################################################
# with open('pipeline/pipeline.png', 'wb') as openfile:
#       pipeline_printout_graph(openfile, output_format='png')
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')
