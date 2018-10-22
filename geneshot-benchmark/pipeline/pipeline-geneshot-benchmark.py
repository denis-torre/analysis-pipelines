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
            gmt[split_line[0]] = split_line[2:]
    return gmt

#######################################################
#######################################################
########## S1. RIF Processing
#######################################################
#######################################################

#############################################
########## 1. Convert to feather
#############################################

@follows(mkdir('s1-feather.dir'))

@transform(glob.glob('rawdata.dir/*.tsv'),
           regex(r'.*/(.*).tsv'),
           r's1-feather.dir/\1.feather')

def toFeather(infile, outfile):

    # Read data
    print('Doing {}...'.format(infile))
    data = pd.read_table(infile)

    # Fix dataframe
    if 'autorif' in infile or 'list' in infile:
        data = data.rename(columns={'Unnamed: 0': 'gene_symbol'})
    elif 'generif' in infile:
        data.columns = data.columns[1:].tolist()+['empty']
        data = data.drop('empty', axis=1).rename_axis('gene_symbol').reset_index()

    # Write 
    data.to_feather(outfile)

#############################################
########## 2. Get network edges
#############################################


@follows(mkdir('s2-networks.dir'))

@transform(glob.glob('s1-feather.dir/*_overlap.feather'),
           regex(r'.*/(.*).feather'),
           r's2-networks.dir/\1-edges.txt')

def getNetworks(infile, outfile):

    # Print
    print('Doing {}...'.format(outfile))

    # Read dataframe
    rif_dataframe = pd.read_feather(infile).set_index('gene_symbol')

    # Get nodes
    count_dataframe = pd.Series(index=rif_dataframe.index, data=np.diag(rif_dataframe), name='total_publications').to_frame()

    # Fill diagonal
    np.fill_diagonal(rif_dataframe.values, 0)

    # Filter edges
    edge_dataframe = pd.melt(rif_dataframe.reset_index(), id_vars='gene_symbol').rename(columns={'gene_symbol': 'source', 'variable': 'target', 'value': 'publications'}).query('publications > 15')

    # Get pairs
    edge_dataframe['pair'] = [str(set([rowData['source'], rowData['target']])) for index, rowData in edge_dataframe.iterrows()]

    # Drop duplicates
    edge_dataframe = edge_dataframe.drop_duplicates('pair').drop('pair', axis=1)

    # Write
    edge_dataframe.to_csv(outfile, sep='\t', index=False)
    count_dataframe.to_csv(outfile.replace('edges', 'nodes'), sep='\t')

#######################################################
#######################################################
########## S2. Prepare Data
#######################################################
#######################################################

#############################################
########## 1. Fraction
#############################################

@follows(mkdir('s3-normalized.dir'))

@transform('s1-feather.dir/list_off_co.feather',
           regex(r'.*/(.*).feather'),
           r's3-normalized.dir/\1-fraction.feather')

def normalizeFraction(infile, outfile):

    # Read data
    enrichr_dataframe = pd.read_feather(infile).set_index('gene_symbol')

    # Apply fraction
    enrichr_dataframe = enrichr_dataframe.apply(lambda x: x/max(x)).reset_index()

    # Save
    enrichr_dataframe.to_feather(outfile)

#############################################
########## 2. Binarize genesets
#############################################

@follows(mkdir('s4-libraries_binary.dir'))

@transform(glob.glob('libraries.dir/*'),
           regex(r'.*/(.*).txt'),
           add_inputs('s1-feather.dir/list_off_co.feather'),
           r's4-libraries_binary.dir/\1-binary.feather')

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
########## 1. Get ranks
#############################################

@follows(mkdir('s5-ranks.dir'))

@product(normalizeFraction,
         formatter(r'.*/.*-(.*).feather'),
         binarizeGenesets,
         formatter(r'.*/(.*)-binary.feather'),
         's5-ranks.dir/{1[0][0]}-{1[1][0]}.feather')

def getRanks(infiles, outfile):

    # Get scores
    score_dataframe = pd.read_feather(infiles[0]).set_index('gene_symbol')

    # Get binary dataframe
    binary_dataframe = pd.read_feather(infiles[1]).set_index('gene_symbol').T

    # Initialize dict
    average_score = {x: {} for x in binary_dataframe.index}

    # Loop through libraries
    for library, genes_binary in binary_dataframe.iterrows():
        for gene_symbol, gene_score in score_dataframe.iterrows():
            average_score[library][gene_symbol] = genes_binary.mul(gene_score).replace(0, np.nan).dropna().mean()

    # Get average score
    average_score_dataframe = pd.DataFrame(average_score).reset_index()

    # Write
    average_score_dataframe.to_feather(outfile)

##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=3, verbose=1)
print('Done!')
