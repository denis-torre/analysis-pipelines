#################################################################
#################################################################
############### Archana RNA-seq Analysis ################
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
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
import Support3 as S
import Archana as P

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
ensembl_file = 'rawdata/ensembl/ensembl_transcripts.txt'
metadata_file = 'rawdata/metadata/sample_metadata.txt'

##### 2. R Connection #####
r.source('pipeline/scripts/archana.R')

#######################################################
#######################################################
########## S1. Humble servant of the all knowing god of jupies
#######################################################
#######################################################

#############################################
########## 1. Merge Data
#############################################

@follows(mkdir('s1-expression.dir/merged'))

@collate('rawdata/featureCounts/*.txt',
         regex(r'.*/....(.*).txt'),
         r's1-expression.dir/merged/\1.txt')

def mergeData(infiles, outfile):

    # Results
    results = []

    # Loop through infiles
    for infile in infiles:
        counts = pd.read_table(infile, skiprows=1)
        counts = counts.rename(columns={counts.columns[-1]: 'counts'})
        counts['sample'] = os.path.basename(infile).split('.')[0]
        results.append(counts)

    # Melt
    count_dataframe = pd.concat(results).pivot_table(index='Geneid', columns='sample', values='counts', aggfunc='sum')

    # Write
    count_dataframe.to_csv(outfile, sep='\t')

#############################################
########## 2. Collapse Symbols
#############################################

@follows(mkdir('s1-expression.dir/ensembl'))

@transform(mergeData,
           regex(r'(.*)/merged/(.*).txt'),
           add_inputs(ensembl_file),
           r'\1/ensembl/\2-ensembl.txt')

def collapseGenes(infiles, outfile):

    # Read gene expression
    print(infiles)
    count_dataframe = pd.read_table(infiles[0])

    # Ensembl column
    ensembl_col = 'Transcript stable ID' if 'transcriptID' in outfile else 'Gene stable ID'

    # Read genes
    ensembl_genes = pd.read_table(infiles[1]).rename(columns={ensembl_col: 'Geneid'})[['Geneid', 'Gene name']]

    # Merge
    result_dataframe = count_dataframe.merge(ensembl_genes, on='Geneid').drop('Geneid', axis=1).groupby('Gene name').sum().rename_axis('gene_symbol')

    # Write
    result_dataframe.to_csv(outfile, sep='\t')

#############################################
########## 3. Filter genes
#############################################

@follows(mkdir('s1-expression.dir/filtered'))

def filterJobs():
    for value in [0, .5, 1, 3, 5, 10]:
        for percentage in [0, .1, .25]:
            infile = 's1-expression.dir/ensembl/exon.transcriptID-ensembl.txt'
            outfile = 's1-expression.dir/filtered/counts-{value}-{percentage}.txt'.format(**locals())
            yield [infile, outfile, value, percentage]

@files(filterJobs)

def filterGenes(infile, outfile, value, percentage):

    # Read dataframe
    count_dataframe = pd.read_table(infile, index_col='gene_symbol')

    # Get CPM
    cpm_dataframe = (count_dataframe/count_dataframe.sum())*10**6

    # Get genes
    genes_bool = ((cpm_dataframe > value).mean(axis=1) > percentage)
    genes_filtered = genes_bool.index[genes_bool]

    # Filter
    count_dataframe_filtered = count_dataframe.reindex(genes_filtered)

    # Write
    count_dataframe_filtered.to_csv(outfile, sep='\t')

#######################################################
#######################################################
########## S2. Differential Expression
#######################################################
#######################################################

#############################################
########## 1. Limma
#############################################

@follows(mkdir('s2-differential_expression.dir'))

@subdivide(filterGenes,
           formatter(),
           add_inputs(metadata_file),
           's2-differential_expression.dir/{basename[0]}-*.txt',
           's2-differential_expression.dir/{basename[0]}-')

def runLimma(infiles, outfiles, outfileRoot):

    # Loop through groups
    for group in ['EoE', 'PPI_REE']:

        # Get outfile
        outfile = '{outfileRoot}{group}.txt'.format(**locals())

        # Run limma
        r.run_limma(infiles=list(infiles), group=group, outfile=outfile)

#############################################
########## 2. Merge results
#############################################

@merge(runLimma,
       's2-differential_expression.dir/limma-merged.txt')

def mergeLimma(infiles, outfile):

    # Initialize results
    results = []

    # Loop through infiles
    for infile in infiles:
        limma_dataframe = pd.read_table(infile)
        limma_dataframe['group'] = os.path.basename(infile)[:-len('.txt')]
        results.append(limma_dataframe)

    # Concatenate
    result_dataframe = pd.concat(results)

    # Write
    result_dataframe.to_csv(outfile, sep='\t', index=False)

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
