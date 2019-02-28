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
import sys
import os
import json
import requests
import pandas as pd
import numpy as np
from functools import reduce
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
import Support3 as S
import Mcgrath as P
sys.path.append('../../jupyter-notebook/biojupies-plugins/library/analysis_tools/enrichr/')
import enrichr
sys.path.append('../../jupyter-notebook/biojupies-plugins/library/core_scripts/shared/')
import shared
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

@follows(mkdir('s1-expression.dir'))

@merge(['rawdata/7196780027/HumanHT-12_V4_0_R2_15002873_B.bgx', 'rawdata/7196780027/*.idat'],
        's1-expression.dir/mcgrath-rawdata.rda')

def readIdat(infiles, outfile):

    # Read background
    bgx_file = infiles.pop(0)
    idat_files = np.array(infiles)

    # Read expression
    r.read_idat(idat_files, bgx_file, outfile)

#############################################
########## 2. Normalize
#############################################

@transform(readIdat,
           suffix('rawdata.rda'),
           'normalized.rda')

def normalizeData(infile, outfile):

    # Normalize data
    r.normalize_data(infile, outfile)

#############################################
########## 3. Extract
#############################################

@transform(normalizeData,
           suffix('.rda'),
           '.txt')

def extractNormalizedData(infile, outfile):

    # Extract data
    r.extract_data(infile, outfile)

#######################################################
#######################################################
########## S2. Differential Expression
#######################################################
#######################################################

#############################################
########## 1. Run limma
#############################################

@follows(mkdir('s2-differential_expression.dir'))

@subdivide(normalizeData,
           formatter(),
           add_inputs(sample_metadata),
           's2-differential_expression.dir/*-limma.txt',
           's2-differential_expression.dir/')

def runLimma(infiles, outfiles, outfileRoot):

    # Get comparisons
    comparisons = (['Control', 'Wounded'], ['Control', 'Normal'], ['Normal', 'Wounded'])

    # Loop through comparisons
    for comparison in comparisons:

        # Outfile
        outfile = '{outfileRoot}{comparison[0]}_vs_{comparison[1]}-limma.txt'.format(**locals())

        # Run Limma
        r.run_limma(expression_file=infiles[0], metadata_file=infiles[1], outfile=outfile, comparison=comparison)


#######################################################
#######################################################
########## S3. Enrichr
#######################################################
#######################################################

#############################################
########## 1. Enrichr
#############################################

@follows(mkdir('s3-enrichr.dir'))

@transform(runLimma,
           regex(r'.*/(.*)-limma.txt'),
           r's3-enrichr.dir/\1-listids.json')

def runEnrichr(infile, outfile):

    # Print
    print('Doing {}...'.format(outfile))
    control, perturbation = os.path.basename(infile).split('-')[0].split('_vs_')

    # Read dataframe
    limma_dataframe = pd.read_table(infile).dropna().sort_values(['Symbol', 'probe_variance'], ascending=False).drop_duplicates('Symbol')
    limma_dataframe = limma_dataframe.sort_values('t', ascending=False).set_index('Symbol')

    # Get N
    n = 500

    # Get genesets
    listids = {
            'up': enrichr.submit_enrichr_geneset(limma_dataframe.index[:n].tolist(), 'Genes upregulated in {perturbation} condition (vs {control})'.format(**locals())),
            'down': enrichr.submit_enrichr_geneset(limma_dataframe.index[-n:].tolist(), 'Genes downregulated in {perturbation} condition (vs {control})'.format(**locals()))
    }

    # Write
    with open(outfile, 'w') as openfile:
            json.dump(listids, openfile, indent=4)

#######################################################
#######################################################
########## S4. Enrichment Results
#######################################################
#######################################################

#############################################
########## 1. Get results
#############################################

@follows(mkdir('s4-enrichment_results.dir'))

@transform(runEnrichr,
           regex(r'.*/(.*)-listids.json'),
           r's4-enrichment_results.dir/\1-enrichment.txt')

def getEnrichmentResults(infile, outfile):

	# Read IDs
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
########## S5. L1000FWD
#######################################################
#######################################################

#############################################
########## 1. Submit
#############################################

@follows(mkdir('s5-l1000fwd.dir'))

@transform(runLimma,
           regex(r'.*/(.*)-limma.txt'),
           r's5-l1000fwd.dir/\1-listids.json')

def runL1000FWD(infile, outfile):

   # Print
    print('Doing {}...'.format(outfile))
    control, perturbation = os.path.basename(infile).split('-')[0].split('_vs_')

    # Read dataframe
    limma_dataframe = pd.read_table(infile).dropna().sort_values(['Symbol', 'probe_variance'], ascending=False).drop_duplicates('Symbol')
    limma_dataframe = limma_dataframe.sort_values('t', ascending=False).set_index('Symbol')

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

#######################################################
#######################################################
########## S6. L1000FWD Results
#######################################################
#######################################################

#############################################
########## 1. Get results
#############################################

@follows(mkdir('s6-l1000fwd_results.dir'))

@transform(runL1000FWD,
           regex(r'.*/(.*)-listids.json'),
           r's6-l1000fwd_results.dir/\1-l1000fwd.txt')

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
# with open('pipeline/pipeline.png', 'wb') as openfile:
#       pipeline_printout_graph(openfile, output_format='png')
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')
