#################################################################
#################################################################
############### Geneshot Benchmarking ################
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
from ruffus.combinatorics import *
import sys, glob, os
import numpy as np
import pandas as pd
from multiprocessing.pool import ThreadPool
from sklearn.metrics import roc_auc_score
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
import Support3 as S
import GeneshotBenchmark as P

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####

##### 2. R Connection #####
r.source('pipeline/scripts/geneshot-benchmark.R')

##### 3. Functions #####
# Read GMT file
def readGMT(infile):
    gmt = {}
    with open(infile) as openfile:
        for line in openfile.read().split('\n'):
            split_line = line.strip().split('\t')
            gmt[split_line[0]] = [x.split(',')[0] for x in split_line[2:]]
    return gmt

#######################################################
#######################################################
########## S1. Prepare Data
#######################################################
#######################################################

#############################################
########## 1. Fraction
#############################################

def normalizeJobs():
    for method in ['fraction']:
        yield ['feather.dir/list_off_co.feather', 's1-normalized.dir/{}.feather'.format(method)]

@follows(mkdir('s1-normalized.dir'))

@files(normalizeJobs)

def normalizeCounts(infile, outfile):

    # Read data
    enrichr_dataframe = pd.read_feather(infile).set_index('gene_symbol')

    # Get method
    method = os.path.basename(outfile).split('.')[0]

    # Normalize
    if method == 'fraction':
        enrichr_dataframe = enrichr_dataframe.apply(lambda x: x/max(x)).reset_index()

    # Save
    enrichr_dataframe.to_feather(outfile)

#############################################
########## 2. Binarize genesets
#############################################

@follows(mkdir('s2-libraries_binary.dir'))

@transform(glob.glob('libraries.dir/*'),
           regex(r'.*/(.*).txt'),
           add_inputs('feather.dir/list_off_co.feather'),
           r's2-libraries_binary.dir/\1-binary.feather')

def binarizeGenesets(infiles, outfile):

    # Read GMT
    gmt = readGMT(infiles[0])

    # Read score dataframe
    score_dataframe = pd.read_feather(infiles[1]).set_index('gene_symbol')

    # Process gmt
    gene_dataframe = pd.DataFrame([{'geneset': geneset, 'gene_symbol': gene, 'value': 1} for geneset, genes in gmt.items() for gene in genes])

    # Pivot
    binary_dataframe = gene_dataframe.pivot_table(index='gene_symbol', columns='geneset', values='value', fill_value=0)

    # Add other genes
    binary_dataframe = binary_dataframe.reindex(score_dataframe.index, fill_value=0).reset_index()
    
    # Write
    binary_dataframe.to_feather(outfile)

#############################################
########## 3. Get average scores
#############################################

@follows(mkdir('s3-average_scores.dir'))

@product(normalizeCounts,
         formatter(r'.*/(.*).feather'),
         binarizeGenesets,
         formatter(r'.*/(.*)-binary.feather'),
         's3-average_scores.dir/{1[0][0]}-{1[1][0]}.feather')

def getAverageScores(infiles, outfile):

    # Get scores
    score_dataframe = pd.read_feather(infiles[0]).set_index('gene_symbol')

    # Get binary dataframe
    binary_dataframe = pd.read_feather(infiles[1]).set_index('gene_symbol')

    # Matrix product
    product_dataframe = (score_dataframe.dot(binary_dataframe)/binary_dataframe.sum()).reset_index()

    # Write
    product_dataframe.to_feather(outfile)

#############################################
########## 4. Get AUC scores
#############################################

@follows(mkdir('s4-auc_scores.dir'))

@transform(getAverageScores,
           regex(r'.*/(.*)-(.*).feather'),
           add_inputs(r's2-libraries_binary.dir/\2-binary.feather'),
           r's4-auc_scores.dir/\1-\2-scores.txt')

def getAucScores(infiles, outfile):

    # Get average
    average_score_dataframe = pd.read_feather(infiles[0]).set_index('gene_symbol')

    # Get binary dataframe
    binary_dataframe = pd.read_feather(infiles[1]).set_index('gene_symbol')

    # Initialize dictionary
    auc = {}

    # Loop through genes
    for index, gene_symbol in enumerate(binary_dataframe.index):
        
        # Get AUC
        binary = binary_dataframe.iloc[index]
        if any(binary):
            auc[gene_symbol] = roc_auc_score(binary, average_score_dataframe.iloc[index])
        else:
            auc[gene_symbol] = np.nan
            
    # Merge
    auc_dataframe = pd.Series(auc, name='auc').to_frame().rename_axis('gene_symbol')

    # Write
    auc_dataframe.to_csv(outfile, sep='\t')

#############################################
########## 5. Merge AUC scores
#############################################

@follows(mkdir('s5-merged_auc.dir'))

@merge(getAucScores,
       's5-merged_auc.dir/merged_auc.txt')

def mergeAucScores(infiles, outfile):

    # Initialize result
    result_dataframe = pd.DataFrame()

    # Loop through infiles
    for infile in infiles:
        
        # Get metadata
        normalization, library, scores = os.path.basename(infile).split('-')
        
        # Read AUC
        auc_dataframe = pd.read_table(infile)
        
        # Add metadata
        auc_dataframe['normalization'] = normalization
        auc_dataframe['library'] = library
        
        # Append
        result_dataframe = result_dataframe.append(auc_dataframe)
        
    # Write
    result_dataframe.to_csv(outfile, sep='\t', index=False)


##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=5, verbose=1)
print('Done!')
