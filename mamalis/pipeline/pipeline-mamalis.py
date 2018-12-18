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
import sys, os, json, glob
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
########## 2. Merge JSON
#############################################

@merge('s1-kallisto.dir/*/run_info.json',
       's2-expression.dir/run_info.txt')

def mergeInfo(infiles, outfile):

	# Initialize result dict
	results = {}

	# Loop through infiles
	for infile in infiles:

		# Get sample name
		sample_name = infile.split('/')[1]

		# Add results
		with open(infile) as openfile:
			results[sample_name] = json.load(openfile)

	# Convert to dataframe
	result_dataframe = pd.DataFrame(results).T.rename_axis('sample_name')

	# Write
	result_dataframe.to_csv(outfile, sep='\t')

#############################################
########## 3. Merge counts
#############################################

@merge(['s1-kallisto.dir/*/abundance.tsv',
	   'rawdata/ensembl/mart_export.txt'],
	   's2-expression.dir/kallisto-counts.txt')

def mergeExpression(infiles, outfile):

	# Split infiles
	annotation_file = infiles.pop()
	abundance_files = infiles

	# Read data
	dataframes = []
	for infile in abundance_files:
		sample_dataframe = pd.read_table(infile)
		sample_dataframe['sample'] = infile.split('/')[1].split('_S')[0]
		dataframes.append(sample_dataframe)
	concatenated_dataframe = pd.concat(dataframes)

	# Read annotation
	annotation_dataframe = pd.read_table(annotation_file)

	# Cast data
	expression_dataframe = concatenated_dataframe.pivot_table(index='target_id', columns='sample', values='est_counts')

	# Merge
	merged_dataframe = expression_dataframe.merge(annotation_dataframe, left_index=True, right_on='Transcript stable ID version', how='inner')

	# Group counts to gene level
	result_dataframe = merged_dataframe.drop(['Transcript stable ID version', 'Gene type'], axis=1).groupby('Gene name').sum().astype(int).rename_axis('gene_symbol')

	# Write
	result_dataframe.to_csv(outfile, sep='\t')

#############################################
########## 4. Split data
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
