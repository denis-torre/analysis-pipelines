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
import sys
import os
import json
import glob
import urllib.request
import pandas as pd
import numpy as np
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
# import Support3 as S
# import Mamalis as P
sys.path.append('../../jupyter-notebook/biojupies-plugins/library/analysis_tools/enrichr/')
import enrichr
sys.path.append('../../jupyter-notebook/biojupies-plugins/library/core_scripts/shared/')
import shared

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
counts_file = 'rawdata/counts/counts.txt'
metadata_file = 'rawdata/metadata/sample_metadata_processed.xlsx'
metadata_file_txt = 'rawdata/metadata/sample_metadata_processed.txt'

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
       's1-kallisto.dir/run_info.txt')

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
########## S3. Differential Expression
#######################################################
#######################################################

#############################################
########## 1. RL Treatment
#############################################

@follows(mkdir('s3-differential_expression.dir'))

@transform(splitData,
		   regex(r'.*/(.*)-counts.txt'),
		   add_inputs(metadata_file_txt),
		   r's3-differential_expression.dir/\1-control_vs_RL-limma.txt')

def limmaRl(infiles, outfile):

	r.limma_rl(np.array(infiles), outfile)


#######################################################
#######################################################
########## S4. Enrichr
#######################################################
#######################################################

#############################################
########## 1. Enrichr
#############################################

@follows(mkdir('s4-enrichr.dir'))

@transform(limmaRl,
           regex(r'.*/(.*)-limma.txt'),
		   r's4-enrichr.dir/\1-listids.json')

def runEnrichr(infile, outfile):

	# Print
	print('Doing {}...'.format(outfile))
	control, perturbation = os.path.basename(infile).split('-')[1].split('_vs_')
	cell_line = os.path.basename(infile).split('_')[0]

	# Read dataframe
	limma_dataframe = pd.read_table(infile).set_index('gene_symbol').sort_values('t', ascending=False)

	# Get N
	n = 500

	# Get genesets
	listids = {
		'up': enrichr.submit_enrichr_geneset(limma_dataframe.index[:n].tolist(), 'Genes upregulated in {perturbation} condition (vs {cell_line} {control})'.format(**locals())),
		'down': enrichr.submit_enrichr_geneset(limma_dataframe.index[-n:].tolist(), 'Genes downregulated in {perturbation} condition (vs {cell_line} {control})'.format(**locals()))
	}

	# Write
	with open(outfile, 'w') as openfile:
		json.dump(listids, openfile, indent=4)

#######################################################
#######################################################
########## S5. Enrichment Results
#######################################################
#######################################################

#############################################
########## 1. Get results
#############################################

@follows(mkdir('s5-enrichment_results.dir'))

@transform(runEnrichr,
           regex(r'.*/(.*)-listids.json'),
           r's5-enrichment_results.dir/\1-enrichment.txt')

def getEnrichmentResults(infile, outfile):

	# Read IDs
	print('Doing {}...'.format(outfile))	
	with open(infile) as openfile:
		enrichr_ids = json.load(openfile)

	# Libraries
	libraries = ['GO_Biological_Process_2018', 'GO_Molecular_Function_2018', 'Reactome_2016', 'KEGG_2016','WikiPathways_2016', 'Human_Phenotype_Ontology', 'MGI_Mammalian_Phenotype_2017']

	# Results
	results = []

	# Loop
	for direction in ['up', 'down']:

		# Get results
		enrichment_results = shared.get_enrichr_results(enrichr_ids[direction]['userListId'], gene_set_libraries={x: x for x in libraries})
		
		# Add direction
		enrichment_results['direction'] = direction

		# Append
		results.append(enrichment_results)
		
	# Concatenate
	result_dataframe = pd.concat(results)

	# Write
	result_dataframe.to_csv(outfile, sep='\t', index=False)

#######################################################
#######################################################
########## S6. Excel Conversion
#######################################################
#######################################################

#############################################
########## 1. Convert
#############################################

@follows(mkdir('s6-excel.dir'))

@transform((limmaRl, getEnrichmentResults),
		   regex(r'.*/(.*).txt'),
		   r's6-excel.dir/\1.xls')

def convertExcel(infile, outfile):
	dataframe = pd.read_table(infile)
	dataframe.to_excel(outfile, index=False)

#######################################################
#######################################################
########## S7. ChEA
#######################################################
#######################################################

#############################################
########## 1. Run ChEA
#############################################

@follows(mkdir('s7-chea.dir'), convertExcel)

@transform(limmaRl,
           regex(r'.*/(.*)-limma.txt'),
           r's7-chea.dir/\1-chea.txt')

def runChea(infile, outfile):

	# Print
	print('Doing {}...'.format(outfile))
	control, perturbation = os.path.basename(infile).split('-')[1].split('_vs_')
	cell_line = os.path.basename(infile).split('_')[0]

	# Read dataframe
	limma_dataframe = pd.read_table(infile).set_index('gene_symbol').sort_values('t', ascending=False)

	# Get N
	n = 500

	# Get gene sets
	genesets = {
		'up': limma_dataframe.index[:n],
		'down': limma_dataframe.index[-n:]
	}

	# Define results
	results = []

	# Loop
	for key, value in genesets.items():

		# Get URL
		chea_url = 'https://amp.pharm.mssm.edu/chea3/api/enrich/'+','.join(value)

		# Get results
		chea_results =json.loads(urllib.request.urlopen(chea_url).read())
		chea_dataframe = pd.DataFrame(chea_results['Integrated--meanRank']).drop(['Library', 'Query Name'], axis=1)
		chea_dataframe['direction'] = key
		
		# Append
		results.append(chea_dataframe)

	# Concatenate
	result_dataframe = pd.concat(results)

	# Write
	result_dataframe.to_csv(outfile, sep='\t', index=False)

##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')
