#################################################################
#################################################################
############### Baron RNA-seq Analysis ################
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
import sys, os, json
import pandas as pd
import numpy as np
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()

##### 2. Custom modules #####
# Pipeline running
sys.path.append('pipeline/scripts')
import Baron as P

# Alignment
sys.path.append('../alignment/pipeline')
import align
sys.path.append('../../jupyter-notebook/biojupies-plugins/library/analysis_tools/enrichr/')
import enrichr
sys.path.append('../../jupyter-notebook/biojupies-plugins/library/core_scripts/shared/')
import shared


#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####

##### 2. R Connection #####
r.source('pipeline/scripts/baron.R')

##### 3. Functions #####
# Rename samples
def rename_samples(colnames, names={'P4': 'P1', 'P5': 'P2', 'P6': 'P3', 'P7': 'DN'}):
	for key, value in names.items():
		colnames = [x.replace(key, value) for x in colnames]
	return colnames

#######################################################
#######################################################
########## S1. Align
#######################################################
#######################################################

#############################################
########## 1. STAR
#############################################

@follows(mkdir('s1-alignment.dir/star'))

@transform('rawdata.dir/*/*.fastq.gz',
           regex(r'rawdata.dir/Sample (.*)/.*.fastq.gz'),
           r's1-alignment.dir/star/\1/')

def runStar(infile, outfile):

    # Replace
	infile = infile.replace(' ', '\ ')

    # Directory
	if not os.path.exists(outfile):
		os.makedirs(outfile)

    # Align
	print('Doing {}...'.format(outfile))
	align.align(infile, outfile, method='STAR', genome='Mus_musculus.GRCm38')


#############################################
########## 2. kallisto
#############################################

@follows(mkdir('s1-alignment.dir/kallisto'))

@transform('rawdata.dir/*/*.fastq.gz',
           regex(r'rawdata.dir/Sample (.*)/.*.fastq.gz'),
           r's1-alignment.dir/kallisto/\1/')

def runKallisto(infile, outfile):

	# Replace
	infile = infile.replace(' ', '\ ')

	# Align
	print('Doing {}...'.format(outfile))
	align.align(infile, outfile, method='kallisto', genome='Mus_musculus.GRCm38')


#######################################################
#######################################################
########## S2. Merge
#######################################################
#######################################################

#############################################
########## 1. STAR
#############################################

@follows(mkdir('s2-expression.dir'))

@merge('s1-alignment.dir/star/*/ReadsPerGene.out.tab',
       's2-expression.dir/baron-star_counts.txt')

def mergeStar(infiles, outfile):
    
	# Merge
	count_dataframe = align.merge(infiles, method='STAR', annotation='GRCm38.p6', filter_biotype=['protein_coding'])
	count_dataframe.columns = rename_samples(count_dataframe.columns)

	# Write
	count_dataframe.to_csv(outfile, sep='\t')

#############################################
########## 2. kallisto
#############################################

@follows(mergeStar)

@merge('s1-alignment.dir/kallisto/*/abundance.tsv',
       's2-expression.dir/baron-kallisto_counts.txt')

def mergeKallisto(infiles, outfile):

	# Merge
	count_dataframe = align.merge(infiles, method='kallisto', annotation='GRCm38.p6', filter_biotype=['protein_coding'])
	count_dataframe.columns = rename_samples(count_dataframe.columns)

	# Write
	count_dataframe.to_csv(outfile, sep='\t')

#######################################################
#######################################################
########## S3. Differential Expression
#######################################################
#######################################################

#############################################
########## 1. limma
#############################################

@follows(mkdir('s3-differential_expression.dir'))

@subdivide(mergeStar,
		   formatter(),
		   's3-differential_expression.dir/baron-*',
		   's3-differential_expression.dir/baron-')

def runLimma(infile, outfiles, outfileRoot):

	# Loop through comparisons
	for comparison in [('P1', 'P3')]:

		# Get outfile
		outfile = '{outfileRoot}{comparison[0]}_vs_{comparison[1]}-limma.txt'.format(**locals())

		# Fix comparison
		comparison = np.array(['group'+x for x in comparison])

		# Run limma
		r.run_limma(infile, outfile, np.array(comparison))

#######################################################
#######################################################
########## S4. Enrichr
#######################################################
#######################################################

#############################################
########## 1. Enrichr
#############################################

@follows(mkdir('s4-enrichr.dir'))

@transform(runLimma,
           regex(r'.*/(.*)-limma.txt'),
		   r's4-enrichr.dir/\1-listids.json')

def runEnrichr(infile, outfile):

	# Print
	print('Doing {}...'.format(outfile))
	control, perturbation = os.path.basename(infile).split('-')[1].split('_vs_')

	# Read dataframe
	limma_dataframe = pd.read_table(infile, index_col='gene_symbol').sort_values('t', ascending=False)

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
	with open(infile) as openfile:
		enrichr_ids = json.load(openfile)

	# Libraries
	libraries = ['GO_Biological_Process_2018']

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

##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')
