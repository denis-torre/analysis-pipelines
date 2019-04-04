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
#%%
from ruffus import *
import sys
import os
import json
import glob
import requests
import urllib.request
import pandas as pd
import numpy as np

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
#%%

##### 2. R Connection #####
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()
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

@follows(alignReads)

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

@follows(alignReads)

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

#######################################################
#######################################################
########## S3. Differential Expression
#######################################################
#######################################################

#############################################
########## 1. RL Treatment
#############################################

@follows(mkdir('s3-differential_expression.dir'))

@merge([mergeExpression, metadata_file_txt],
	   's3-differential_expression.dir/timepoint-all-limma.txt')

def limmaRL(infiles, outfile):

	r.limma_rl(infiles, outfile)

#############################################
########## 2. Timepoint
#############################################

@subdivide(mergeExpression,
		   formatter(),
		   add_inputs(metadata_file_txt),
           's3-differential_expression.dir/timepoint-*',
           's3-differential_expression.dir/timepoint-')

def limmaTimepoint(infiles, outfile, outfileRoot):

	# Read metadata
	sample_metadata_dataframe = pd.read_table(infiles[1], index_col='sample')

	# Loop through timepoints
	for timepoint in sample_metadata_dataframe['timepoint'].unique():

		# Get outfile
		outfile = '{outfileRoot}{timepoint}-limma.txt'.format(**locals())

		# Timepoint
		r.limma_timepoint(np.array(infiles), outfile)

#######################################################
#######################################################
########## S4. Enrichr
#######################################################
#######################################################

#############################################
########## 1. Enrichr
#############################################

@follows(mkdir('s4-enrichr.dir'), limmaRL, limmaTimepoint)

@transform('s3-differential_expression.dir/*-limma.txt',
           regex(r'.*/(.*)-limma.txt'),
		   r's4-enrichr.dir/\1-listids.json')

def runEnrichr(infile, outfile):

	# Print
	print('Doing {}...'.format(outfile))
	timepoint = infile.split('-')[2]

	# Read dataframe
	limma_dataframe = pd.read_table(infile).set_index('gene_symbol').sort_values('t', ascending=False)

	# Get N
	n = 500

	# Get genesets
	listids = {
		'up': enrichr.submit_enrichr_geneset(limma_dataframe.index[:n].tolist(), 'Genes upregulated in RL at {timepoint} (all fibroblasts)'.format(**locals())),
		'down': enrichr.submit_enrichr_geneset(limma_dataframe.index[-n:].tolist(), 'Genes downregulated in RL at {timepoint} (all fibroblasts)'.format(**locals()))
	}

	# Write
	with open(outfile, 'w') as openfile:
		json.dump(listids, openfile, indent=4)

#############################################
########## 2. Get results
#############################################

# @follows(mkdir('s5-enrichment_results.dir'))

# @transform(runEnrichr,
#            regex(r'.*/(.*)-listids.json'),
#            r's5-enrichment_results.dir/\1-enrichment.txt')

@transform(runEnrichr,
           suffix('listids.json'),
           'enrichment.txt')

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
########## S6. Timepoint Analysis
#######################################################
#######################################################

#############################################
########## 1. Merge T
#############################################

@follows(mkdir('s6-timepoint.dir'))

@merge(limmaTimepoint,
	   's6-timepoint.dir/timepoint-table.txt')

def mergeT(infiles, outfile):

	# Initialize results
	results = []

	# Loop through infiles
	for infile in infiles:

		# Read data
		limma_dataframe = pd.read_table(infile)

		# Add timepoint
		limma_dataframe['timepoint'] = infile.split('-')[2]

		# Append
		results.append(limma_dataframe)

	# Result
	result_dataframe = pd.concat(results).pivot(index='gene_symbol', columns='timepoint', values='t').dropna()

	# Write
	result_dataframe.to_csv(outfile, sep='\t')

#############################################
########## 2. Cluster
#############################################

@transform(mergeT,
		   suffix('table.txt'),
		   'hclust.rda')

def clusterGenes(infile, outfile):

	# Cluster
	r.cluster_genes(infile, outfile)

#############################################
########## 3. Groups
#############################################

@transform(clusterGenes,
		   suffix('hclust.rda'),
		   'groups.txt')

def cutTree(infile, outfile):

	# Cut tree
	r.cut_tree(infile, outfile)

#############################################
########## 4. Line Plots
#############################################

@transform(mergeT,
		   suffix('table.txt'),
		   add_inputs(cutTree),
		   'lineplots.pdf')

def linePlots(infiles, outfile):

	# Cut tree
	r.line_plots(np.array(infiles), outfile)

#############################################
########## 5. NbClust
#############################################

@transform(clusterGenes,
           suffix('hclust.rda'),
		   'nbclust.rda')

def nbClust(infile, outfile):

	# Nbclust
	r.nbclust(infile, outfile)

#############################################
########## 6. Enrichment
#############################################

@transform(cutTree,
		   suffix('.txt'),
		   '-listids.txt')

def timepointEnrichr(infile, outfile):

	# Read groups
	group_dataframe = pd.read_table(infile).reset_index()

	# Group
	group_dict = group_dataframe.groupby('k15')['index'].apply(list).to_dict()

	# List IDs
	listids = {key: enrichr.submit_enrichr_geneset(value, 'Genes in cluster {key}'.format(**locals())) for key, value in group_dict.items()}

	# Write
	with open(outfile, 'w') as openfile:
		json.dump(listids, openfile, indent=4)

#############################################
########## 7. Enrichment Results
#############################################

@transform(timepointEnrichr,
           suffix('listids.txt'),
		   'enrichment.txt')

def timepointEnrichmentResults(infile, outfile):

	# Read
	with open(infile) as openfile:
		listids = json.load(openfile)

	# Libraries
	libraries = ['GO_Biological_Process_2018', 'GO_Molecular_Function_2018', 'GO_Cellular_Component_2018',
              'Reactome_2016', 'KEGG_2016', 'WikiPathways_2019_Human', 'ARCHS4_TFs_Coexp', 'ARCHS4_Kinases_Coexp']

	# Results
	results = []

	# Loop
	for key, value in listids.items():

		# Get results
		enrichment_results = shared.get_enrichr_results(value['userListId'], gene_set_libraries={x: x for x in libraries})
		
		# Add direction
		enrichment_results['group'] = key

		# Append
		results.append(enrichment_results)
		
	# Concatenate
	result_dataframe = pd.concat(results)

	# Write
	result_dataframe.to_csv(outfile, sep='\t', index=False)

#######################################################
#######################################################
########## S5. L1000FWD
#######################################################
#######################################################

#############################################
########## 1. Submit
#############################################

@follows(mkdir('s7-l1000fwd.dir'))

@transform((limmaRL, limmaTimepoint),
           regex(r'.*/(.*)-limma.txt'),
           r's7-l1000fwd.dir/\1-listids.json')

def runL1000FWD(infile, outfile):

   # Print
    print('Doing {}...'.format(outfile))

    # Read dataframe
    limma_dataframe = pd.read_table(infile).sort_values('t', ascending=False).set_index('gene_symbol')

    # Get N
    n = 500

    # Prepare payload
    payload = {
        'up_genes': limma_dataframe.index[:n].tolist(),
        'down_genes': limma_dataframe.index[-n:].tolist()
    }

    # Get results
    response = requests.post('http://amp.pharm.mssm.edu/L1000FWD/sig_search', json=payload)
    with open(outfile, 'w') as openfile:
        json.dump(response.json(), openfile, indent=4)

#############################################
########## 2. Get results
#############################################

# @follows(mkdir('s8-l1000fwd_results.dir'))

# @transform(runL1000FWD,
#            regex(r'.*/(.*)-listids.json'),
#            r's8-l1000fwd_results.dir/\1-l1000fwd.txt')

@transform(runL1000FWD,
           suffix('listids.json'),
           'l1000fwd.txt')

def getL1000FWDResults(infile, outfile):

   # Print
    print('Doing {}...'.format(outfile))

# Read ID
    with open(infile) as openfile:
        result_id = json.load(openfile)['result_id']

    # Perform request
    response = requests.get('http://amp.pharm.mssm.edu/L1000FWD/result/topn/' + result_id)

    # Get results
    results = []

    # Loop
    for direction, signatures in response.json().items():
        
        # Add metadata
        for signature in signatures:
            signature.update(requests.get('http://amp.pharm.mssm.edu/L1000FWD/sig/{sig_id}'.format(**signature)).json())
        
        # Append
        dataframe = pd.DataFrame(signatures)
        dataframe['direction'] = direction
        results.append(dataframe)

    # Concatenate
    result_dataframe = pd.concat(results).drop(['up_genes', 'down_genes', 'combined_genes'], axis=1)

    # Write
    result_dataframe.to_csv(outfile, sep='\t', index=False)

##################################################
##################################################
########## Run pipeline
##################################################
##################################################
if sys.argv[-1] == 'plot':
	with open('pipeline/pipeline.png', 'wb') as openfile:
		pipeline_printout_graph(openfile, output_format='png')
elif sys.argv[-1] == 'force':
	pipeline_run(forcedtorun_tasks=[sys.argv[-2]], multiprocess=1, verbose=1)
else:
	pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')


#%%
