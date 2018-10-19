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
import sys, glob
import numpy as np
import pandas as pd
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

#######################################################
#######################################################
########## S1. Process Data
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
pipeline_run([sys.argv[-1]], multiprocess=3, verbose=1)
print('Done!')
