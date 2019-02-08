#################################################################
#################################################################
############### RNA-seq Alignment ################
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
import sys, glob, os
import pandas as pd

##### 2. Custom modules #####


#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####

##### 2. R Connection #####

#######################################################
#######################################################
########## S1. Kallisto
#######################################################
#######################################################

#############################################
########## 1. Build Index
#############################################

@follows(mkdir('s1-kallisto.dir'))

@transform('ensembl.dir/*/*.cdna.all.fa.gz',
           regex(r'.*/(.*).fa.gz'),
		   r's1-kallisto.dir/\1.idx')

def buildKallistoIndex(infile, outfile):

	# Print
	print('Doing '+outfile+'...')

	# Execute
	os.system('kallisto index -i {outfile} {infile}'.format(**locals()))

#######################################################
#######################################################
########## S2. Star
#######################################################
#######################################################

#############################################
########## 1. Build Index
#############################################

@follows(mkdir('s2-star.dir'))

@transform('ensembl.dir/*/*.dna.primary_assembly.fa',
           regex(r'(.*)/(.*).dna.primary_assembly.fa'),
		   add_inputs(r'\1/*.gtf'),
		   r's2-star.dir/\2/genomeDir')

def buildStarIndex(infiles, outfile):

	# Create directory
	if not os.path.exists(outfile):
		os.makedirs(outfile)
	outdir = os.path.dirname(outfile)

	# Statement
	statement = '''
		STAR \
			--runMode genomeGenerate \
			--genomeDir {outfile} \
			--genomeFastaFiles {infiles[0]} \
			--sjdbGTFfile {infiles[1]} \
			--outFileNamePrefix {outdir}/ \
			--runThreadN 3 \
			--genomeSAsparseD 2 \
			--genomeSAindexNbases 13
	'''.format(**locals())

	os.system(statement)

##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')
