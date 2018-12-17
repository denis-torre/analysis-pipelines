#################################################################
#################################################################
############### Mamalis RNA-seq Analysis ################
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
import Mamalis as P

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
counts_file = 'rawdata/counts/counts.txt'
metadata_file = 'rawdata/metadata/sample_metadata_processed.xlsx'

##### 2. R Connection #####
r.source('pipeline/scripts/mamalis.R')

#######################################################
#######################################################
########## S1. Quantify Expression
#######################################################
#######################################################

#############################################
########## 1. Align reads
#############################################

@follows(mkdir('s1-kallisto.dir'))

@collate('rawdata/fastq/*.fq.gz',
		 regex(r'.*/(.*)_..fq.gz'),
         r's1-kallisto.dir/\1')

def alignReads(infiles, outfile):

	# Align
	print('Doing {}...'.format(outfile))
	try:
		os.system('kallisto quant -t 8 -i rawdata/kallisto/Homo_sapiens_index.idx -o {outfile} {infiles[0]} {infiles[1]}'.format(**locals()))
	except:
		print('Error doing {}.'.format(outfile))

#############################################
########## 2. Merge data
#############################################

#############################################
########## 3. Split data
#############################################

@follows(mkdir('s2-expression.dir'))

@subdivide(counts_file,
		   formatter(),
		   add_inputs(metadata_file),
		   's2-expression.dir/*-counts.txt',
           's2-expression.dir/')

def splitData(infiles, outfiles, outfileRoot):

	# Read counts
	count_dataframe = pd.read_table(infiles[0], index_col='gene_symbol')

	# Read sample metadata
	metadata_dataframe = pd.read_excel(infiles[1])

	# Get groups
	sample_groups = metadata_dataframe.groupby('cell_line')['sample'].apply(list).to_dict()

	# Loop through groups
	for group_name, group_samples in sample_groups.items():

		# Get outfile
		outfile = '{outfileRoot}{group_name}-counts.txt'.format(**locals())

		# Write subset
		count_dataframe[group_samples].to_csv(outfile, sep='\t')

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
pipeline_run([sys.argv[-1]], multiprocess=4, verbose=1)
print('Done!')
