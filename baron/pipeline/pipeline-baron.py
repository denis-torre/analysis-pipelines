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
import sys, os
import pandas as pd
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()

##### 2. Custom modules #####
# Pipeline running
sys.path.append('pipeline/scripts')
import Baron as P

# Alignment
sys.path.append('../alignment/pipeline')
import align

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
