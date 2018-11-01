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
import scipy.stats as ss
import seaborn as sns
import matplotlib.pyplot as plt
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
libraries = glob.glob('libraries.dir/*.txt')

##### 2. R Connection #####
r.source('pipeline/scripts/geneshot-benchmark.R')

##### 3. Functions #####
# Read GMT file
def readGMT(infile):
    gmt = {}
    with open(infile) as openfile:
        for line in openfile.read().split('\n'):
            split_line = line.strip().split('\t')
            if split_line[0]:
                gmt[split_line[0]] = [x.split(',')[0] for x in split_line[2:]]
    return gmt

# Reverse GMT
def reverseGMT(gmt):
    gmt_list = [{'term': term, 'gene_symbol': gene, 'value': 1} for term, genes in gmt.items() for gene in genes]
    gmt_dataframe = pd.DataFrame(gmt_list).pivot_table(index='gene_symbol', columns='term', values='value', fill_value=0)
    reverse_gmt = {gene_symbol: [key for key, value in rowData.items() if value] for gene_symbol, rowData in gmt_dataframe.iterrows()}
    return reverse_gmt

# Z-score
def zscoreDF(df): return ((df.T - df.T.mean())/df.T.std()).T

#######################################################
#######################################################
########## S1. Prepare Data
#######################################################
#######################################################

#############################################
########## 1. Fraction
#############################################

def normalizeJobs():
    for method in []:#['raw', 'fraction', 'zscore', 'random']:
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
        normalized_dataframe = enrichr_dataframe.apply(lambda x: x/max(x), axis=1)
    elif method == 'raw':
        normalized_dataframe = enrichr_dataframe
    elif method == 'zscore':
        normalized_dataframe = zscoreDF(enrichr_dataframe)
    elif method == 'zscore_nodiag':
        enrichr_dataframe = enrichr_dataframe.astype(float)
        np.fill_diagonal(enrichr_dataframe.values, np.nan)
        normalized_dataframe = zscoreDF(enrichr_dataframe)
    elif method == 'random':
        normalized_dataframe = enrichr_dataframe.T.apply(lambda x: np.random.rand(len(x))).T
    else:
        raise ValueError

    # Save
    normalized_dataframe.reset_index().to_feather(outfile)

#############################################
########## 2. Get average scores
#############################################

@follows(mkdir('s2-average_scores.dir'), normalizeCounts)

@product(glob.glob('s1-normalized.dir/*.feather'),
         formatter(r'.*/(.*).feather'),
         libraries,
         formatter(r'.*/(.*).txt'),
         's2-average_scores.dir/{1[0][0]}-{1[1][0]}.feather')

def getAverageScores(infiles, outfile):

    # Print
    print('Doing {}...'.format(outfile))

    # Read scores
    score_dataframe = pd.read_feather(infiles[0]).set_index('gene_symbol').astype(float)
    np.fill_diagonal(score_dataframe.values, np.nan)

    # Read GMT
    gmt = readGMT(infiles[1])

    # Initialize results
    results = {}

    # Loop through data
    for term, term_genes in gmt.items():
        results[term] = score_dataframe.reindex(term_genes, axis=1).mean(axis=1)

    # Get dataframe
    results_dataframe = pd.DataFrame(results).rename_axis('gene_symbol').reset_index()

    # Save
    results_dataframe.to_feather(outfile)

#############################################
########## 3. Get AUC scores
#############################################

@follows(mkdir('s3-auc_scores.dir'))

@transform(getAverageScores,
           regex(r'.*/(.*)-(.*).feather'),
           add_inputs(r'libraries.dir/\2.txt'),
           r's3-auc_scores.dir/\1-\2-library_scores.txt')

def getAucScores(infiles, outfile):

    # Print
    print('Doing {}...'.format(outfile))

    # Get average
    average_score_dataframe = pd.read_feather(infiles[0]).set_index('gene_symbol').dropna(axis=1)

    # Get binary dataframe
    gmt = readGMT(infiles[1])

    # Initialize results
    auc = {}

    # Loop through columns
    for column, colData in average_score_dataframe.items():
        
        # Binarize terms
        term_binary = [x in gmt[column] for x in average_score_dataframe.index]
        
        # Store results
        if any(term_binary):
            auc[column] = roc_auc_score(term_binary, colData)
        else:
            auc[column] = np.nan
            
    # Merge
    auc_dataframe = pd.Series(auc, name='auc').to_frame().rename_axis('term_name')

    # Write
    auc_dataframe.to_csv(outfile, sep='\t')

#############################################
########## 4. Get AUC scores
#############################################

@follows(mkdir('s3-auc_scores.dir'))

@transform(getAverageScores,
           regex(r'.*/(.*)-(.*).feather'),
           add_inputs(r'libraries.dir/\2.txt'),
           r's3-auc_scores.dir/\1-\2-gene_scores.txt')

def getGeneAucScores(infiles, outfile):

    # Print
    print('Doing {}...'.format(outfile))

    # Get average
    average_score_dataframe = pd.read_feather(infiles[0]).set_index('gene_symbol').dropna(axis=1)

    # Get binary dataframe
    gmt = readGMT(infiles[1])
    reverse_gmt = reverseGMT(gmt)

    # Initialize results
    auc = {}

    # Loop through columns
    for index, rowData in average_score_dataframe.iterrows():

        # Binarize genes
        gene_binary = [x in reverse_gmt.get(index, []) for x in average_score_dataframe.columns]

        # Store results
        if sum(gene_binary) > 10:
            auc[index] = roc_auc_score(gene_binary, rowData)
        else:
            auc[index] = np.nan

    # Merge
    auc_dataframe = pd.Series(auc, name='auc').to_frame().rename_axis('gene_symbol')

    # Write
    auc_dataframe.to_csv(outfile, sep='\t')

#############################################
########## 5. Merge AUC scores
#############################################

@follows(mkdir('s4-merged_auc.dir'), getAucScores, getGeneAucScores)

# @merge([getAucScores, glob.glob('libraries.dir/*.txt')],
    #    's4-merged_auc.dir/merged_auc.txt')
@collate('s3-auc_scores.dir/*_scores.txt',
         regex(r'.*/.*-(.*)_scores.txt'),
         r's4-merged_auc.dir/\1_auc.txt',
         libraries)

def mergeAucScores(infiles, outfile, libraries):

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

    # Get term sizes
    term_sizes = {term: len(term_genes) for gmt_file in libraries for term, term_genes in readGMT(gmt_file).items()}

    # Add term sizes
    if 'library' in outfile:
        result_dataframe = result_dataframe.merge(pd.Series(term_sizes).rename('nr_genes').to_frame(), left_on='term_name', right_index=True)
        
    # Write
    result_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 6. Plot AUC scores
#############################################

@follows(mkdir('s5-auc_plots.dir'))

@transform(mergeAucScores,
           regex(r'.*/(.*).txt'),
           r's5-auc_plots.dir/\1-.png')

def plotAucScores(infile, outfile):

    # Plot
    na=False
    for plot_type in ['density']: #'violin'
        # for na in [True, False]:
        r.plot_auc(infile, outfile.replace('.png', (plot_type+('_withna' if na else '')+'.png')), plot_type=plot_type, na=na)

##################################################
##################################################
########## Run pipeline
##################################################
##################################################
if 'plot' in sys.argv[-1]:
    pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
else:
    pipeline_run([sys.argv[-1]], multiprocess=4, verbose=1)
print('Done!')
