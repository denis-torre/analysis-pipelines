#################################################################
#################################################################
############### Dubois RNA-seq Analysis ################
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
import sys, glob, os, json
import pandas as pd
import numpy as np
from rpy2.robjects import r, pandas2ri
# pandas2ri.activate()

##### 2. Custom modules #####
# Pipeline running
# sys.path.append('/Users/denis/Documents/Projects/scripts')
# sys.path.append('pipeline/scripts')
# sys.path.append('/Users/denis/Documents/Projects/jupyter-notebook/biojupies-plugins/library/core_scripts/normalize')
# import normalize as N
# sys.path.append('/Users/denis/Documents/Projects/jupyter-notebook/biojupies-plugins/library/core_scripts/signature')
# import signature as S
# import Support3 as S
# import Dubois as P
sys.path.append('../../jupyter-notebook/biojupies-plugins/library/analysis_tools/enrichr/')
import enrichr

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####

##### 2. R Connection #####
r.source('pipeline/scripts/dubois.R')

#######################################################
#######################################################
########## S1. Alignment
#######################################################
#######################################################

#############################################
########## 1. Concatenate Reads
#############################################

@collate(glob.glob('rawdata/fastq/*.fastq.gz'),
		 regex(r'rawdata/fastq/(.*)_(.*)_(.*)_001.fastq.gz'),
		 r'rawdata/fastq_concat/\1_\3_001.fastq.gz')

def concatenateReads(infiles, outfile):

	# Concatenate
	print('Doing {}...'.format(outfile))
	infiles = sorted(infiles)
	os.system('cat {infiles[0]} {infiles[1]} > {outfile}'.format(**locals()))

#############################################
########## 2. Kallisto
#############################################

@follows(mkdir('s1-alignment.dir/kallisto'))

@collate(glob.glob('rawdata/fastq_concat/*.fastq.gz'),
		 regex(r'rawdata/fastq_concat/(.*)_R._001.fastq.gz'),
         r's1-alignment.dir/kallisto/\1')

def runKallisto(infiles, outfile):

	# Align
	print('Doing {}...'.format(outfile))
	try:
		os.system('kallisto quant -t 8 -i rawdata/kallisto/Homo_sapiens_index.idx -o {outfile} {infiles[0]} {infiles[1]}'.format(**locals()))
	except:
		print('Error doing {}.'.format(outfile))

#############################################
########## 3. Star
#############################################

@follows(mkdir('s1-alignment.dir/star'))

@collate('rawdata/fastq/*.fastq.gz',
		 regex(r'rawdata/fastq/(.*)_L00._R._001.fastq.gz'),
         r's1-alignment.dir/star/\1/')

def runStar(infiles, outfile):

	# Directory
	if not os.path.exists(outfile):
		os.makedirs(outfile)

	# Align
	print('Doing {}...'.format(outfile))
	try:
		command = ''' STAR \
			--genomeDir /Users/denis/Data/star/human/genomeDir  \
			--outFileNamePrefix {outfile}  \
			--readFilesIn {infiles[0]},{infiles[2]} {infiles[1]},{infiles[3]}  \
			--readFilesCommand gunzip -c \
			--quantMode GeneCounts \
			--limitBAMsortRAM 10000000000  \
			--limitIObufferSize 50000000 \
			--outSAMstrandField intronMotif \
			--outFilterIntronMotifs RemoveNoncanonical  \
			--outSAMtype BAM SortedByCoordinate  \
			--outReadsUnmapped Fastx \
			--runThreadN 5
		'''.format(**locals())
		os.system(command)
	except:
		print('Error doing {}.'.format(outfile))

#######################################################
#######################################################
########## S2. Merge Expression
#######################################################
#######################################################

#############################################
########## 1. Merge Kallisto Stats
#############################################

@follows(mkdir('s2-expression.dir/kallisto'))

@merge(glob.glob('s1-alignment.dir/kallisto/*/run_info.json'),
       's2-expression.dir/kallisto/kallisto-run_info.txt')

def mergeKallistoInfo(infiles, outfile):

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
########## 2. Merge Kallisto Expression
#############################################

@merge((glob.glob('s1-alignment.dir/kallisto/*/abundance.tsv'),
	   'rawdata/ensembl/mart_export_transcript.txt'),
	   's2-expression.dir/kallisto/kallisto-counts.txt')

def mergeKallistoExpression(infiles, outfile):

	# Split infiles
	abundance_files, annotation_file = infiles

	# Read data
	dataframes = []
	for infile in abundance_files:
		if infile.split('/')[1].split('_')[0] in ['1', '2', '3', '4']:
			print("Skipping sample {}...".format(infile))
		else:
			sample_dataframe = pd.read_table(infile)
			sample_dataframe['sample'] = 'sample_'+infile.split('/')[-2].split('_S')[0]
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

#############################################
########## 3. Merge STAR Expression
#############################################

@follows(mkdir('s2-expression.dir/star'))

@merge((glob.glob('s1-alignment.dir/star/*/ReadsPerGene.out.tab'), 'rawdata/ensembl/mart_export_gene.txt'),
	   's2-expression.dir/star/star-counts.txt')

def mergeStarExpression(infiles, outfile):

	# Split infiles
	abundance_files, annotation_file = infiles

	# Read data
	dataframes = []
	for infile in abundance_files:
		# if infile.split('/')[1].split('_')[0] in ['1', '2', '3', '4']:
			# print("Skipping sample {}...".format(infile))
		# else:
		sample_dataframe = pd.read_table(infile, header=None, names=['gene_id', 'counts_unstranded', 'counts_strand1', 'counts_strand2'])
		sample_dataframe['sample'] = 'sample_'+infile.split('/')[-2].split('_S')[0]
		dataframes.append(sample_dataframe)
	concatenated_dataframe = pd.concat(dataframes)

	# Read annotation
	annotation_dataframe = pd.read_table(annotation_file)

	# Cast data
	expression_dataframe = concatenated_dataframe.pivot_table(index='gene_id', columns='sample', values='counts_unstranded')

	# Merge
	merged_dataframe = expression_dataframe.merge(annotation_dataframe, left_index=True, right_on='Gene stable ID', how='inner')

	# Group counts to gene level
	result_dataframe = merged_dataframe.drop(['Gene stable ID', 'Gene type'], axis=1).groupby('Gene name').sum().astype(int).rename_axis('gene_symbol')

	# Write
	result_dataframe.to_csv(outfile, sep='\t')

#############################################
########## 4. Filter
#############################################

@follows(mkdir('s2-expression.dir/kallisto/filtered'))
@follows(mkdir('s2-expression.dir/star/filtered'))

@subdivide(glob.glob('s2-expression.dir/kallisto/*-counts.txt'),
		   regex(r'(.*)/(.*).txt'),
		   r'\1/filtered/\2-*.txt',
           r'\1/filtered/\2')

def filterGenes(infile, outfile, outfileRoot):

	# Read dataframe
	count_dataframe = pd.read_table(infile, index_col='gene_symbol')

	# Get CPM
	cpm_dataframe = (count_dataframe/count_dataframe.sum())*10**6

	# Loop
	for value in [0, .5, 1, 3, 5, 10]:
		for percentage in [0, .1, .25]:

			# Get genes
			genes_bool = ((cpm_dataframe >= value).mean(axis=1) >= percentage)
			genes_filtered = genes_bool.index[genes_bool]

			# Filter
			count_dataframe_filtered = count_dataframe.reindex(genes_filtered)

			# Outfile
			outfile = '{outfileRoot}-{value}cpm-{percentage}pct.txt'.format(**locals())
			print('Doing {}...'.format(outfile))

			# Write
			count_dataframe_filtered.to_csv(outfile, sep='\t')

#############################################
########## 5. Median Expression
#############################################

@transform('s2-expression.dir/kallisto/*-counts.txt',
		   suffix('-counts.txt'),
		   add_inputs('rawdata/metadata/dubois-metadata.txt'),
		   '-median_group_logcpm.txt')

def getMedianLogcpm(infiles, outfile):

	# Get counts
	count_dataframe = pd.read_table(infiles[0], index_col='gene_symbol')

	# Get logCPM
	logcpm_dataframe = np.log10((count_dataframe/count_dataframe.sum())*10**6+1)
	round(logcpm_dataframe, ndigits=3).to_csv('dubois-logcpm.txt', sep='\t')

	# Melt
	melted_dataframe = pd.melt(logcpm_dataframe.reset_index(), id_vars='gene_symbol', var_name='Sample', value_name='logCPM')

	# Get metadata
	metadata_dataframe = pd.read_table(infiles[1])

	# Merge
	merged_dataframe = melted_dataframe.merge(metadata_dataframe, on='Sample')

	# Group
	median_dataframe = merged_dataframe.groupby(['gene_symbol', 'Condition'])['logCPM'].median().rename('logCPM').reset_index()

	# Cast
	cast_dataframe = median_dataframe.pivot(index='gene_symbol', columns='Condition', values='logCPM')
	cast_dataframe = round(cast_dataframe, ndigits=3)

	# Write
	cast_dataframe.to_csv(outfile, sep='\t')

#######################################################
#######################################################
########## S3. PCA
#######################################################
#######################################################

#############################################
########## 1. Collapse Expression
#############################################

# @follows(mkdir('s3-mean_expression.dir'))

# @transform(normalizeData,
# 		   regex(r'.*/dubois-(.*).txt'),
# 		   add_inputs(processMetadata),
#            r's3-mean_expression.dir/dubiois-grouped_\1.json')

# def collapseExpression(infiles, outfile):

# 	# Get expression
# 	expression_dataframe = pd.read_table(infiles[0])

# 	# Read metadata
# 	metadata_dataframe = pd.read_table(infiles[1])

# 	# Merge
# 	merged_dataframe = pd.melt(expression_dataframe, id_vars='gene_symbol').merge(metadata_dataframe, left_on='variable', right_on='Sample').rename(columns={'variable': 'sample', 'value': 'logCPM'}).drop('Diff', axis=1)

# 	# Group
# 	grouped_dataframe = merged_dataframe.groupby(['gene_symbol', 'Condition'])['logCPM'].apply(list).rename('logCPM').to_frame()

# 	# Initialize dictionary
# 	expression_dict = {x: {} for x in expression_dataframe['gene_symbol']}

# 	# Loop
# 	for index, rowData in grouped_dataframe.iterrows():
# 		expression_dict[index[0]][index[1]] = rowData['logCPM']

# 	# Write
# 	with open(outfile, 'w') as openfile:
# 		openfile.write(json.dumps(expression_dict))

#######################################################
#######################################################
########## S4. limma
#######################################################
#######################################################

#############################################
########## 1. Get Differential Expression
#############################################

def differentialJobs():
	metadata_series = pd.read_table('rawdata/metadata/dubois-metadata.txt').groupby('Condition')['Sample'].apply(tuple)
	control_samples = metadata_series.pop('Control')
	for group_name, perturbed_samples in metadata_series.items():
		group_name = group_name.replace(' ', '')
		for aligner in ['kallisto']:
			yield ('s2-expression.dir/{aligner}/nonbatch/{aligner}-counts.txt'.format(**locals()), 's3-limma.dir/{aligner}/{group_name}-{aligner}-limma.txt'.format(**locals()), group_name, control_samples, perturbed_samples)

@follows(mkdir('s3-limma.dir/star'))
@follows(mkdir('s3-limma.dir/kallisto/nonbatch'))

@files(differentialJobs)

def runDifferentialExpression(infile, outfile, group_name, control_samples, perturbed_samples):

	# Print
	print('Doing {}...'.format(outfile))

	# Read metadata
	count_dataframe = pd.read_table(infile, index_col='gene_symbol').rename_axis('index')

	# Run analysis
	signature_dataframe = S.limma({'rawdata': count_dataframe}, control_samples, perturbed_samples)

	# Add group name
	signature_dataframe['group_name'] = group_name

	# Write
	signature_dataframe.to_csv(outfile, sep='\t')

#############################################
########## 2. Batch differential expression
#############################################

def differentialJobs2():
	metadata_series = pd.read_table('rawdata/metadata/dubois-metadata.txt').groupby('Condition')['Sample'].apply(tuple)
	control_samples = metadata_series.pop('Control')
	for group_name, perturbed_samples in metadata_series.items():
		group_name = group_name.replace(' ', '')
		for aligner in ['kallisto']:
			yield ('s2-expression.dir/{aligner}/{aligner}-counts.txt'.format(**locals()), 's3-limma.dir/{aligner}/{group_name}-{aligner}-batch-limma.txt'.format(**locals()), group_name, control_samples, perturbed_samples)

@follows(mkdir('s3-limma.dir/kallisto/batch'))

@subdivide('s2-expression.dir/kallisto/kallisto-counts.txt',
           regex(r'.*/(.*)-counts.txt'),
           add_inputs('rawdata/metadata/dubois-metadata.txt', 'rawdata/metadata/dubois-comparisons.txt'),
           r's3-limma.dir/\1/.*.txt',
           r's3-limma.dir/\1/batch/\1-')

def runBatchDifferentialExpression(infiles, outfiles, outfileRoot):

	# Fix infiles
	infiles = list(infiles)

	# Comparison file
	comparison_dataframe = pd.read_table(infiles.pop())

	# Loop
	for index, groups in comparison_dataframe.iterrows():

		# Outfile
		outfile = '{outfileRoot}{control_condition}_vs_{perturbation_condition}-batch.txt'.format(**locals(), **groups)

		# Skip if exists
		if not os.path.exists(outfile):

			# Print
			print('Doing {}...'.format(outfile))

			# Run limma
			r.run_limma(
				count_file=infiles[0],
				metadata_file=infiles[1],
				control_condition=groups['control_condition'],
				perturbation_condition=groups['perturbation_condition'],
				outfile=outfile
			)

#############################################
########## 3. Merge limma
#############################################

@merge('s3-limma.dir/kallisto/batch/*.txt',
	   's3-limma.dir/kallisto/limma-batch-merged.txt')

def mergeLimma(infiles, outfile):

	# Data
	data = []

	# Loop
	for infile in infiles:
		
		# Groups
		control, perturbation = os.path.basename(infile).split('-')[1].split('_vs_')
		
		# Read
		dataframe = pd.read_table(infile)
		dataframe['control_condition'] = control
		dataframe['perturbation_condition'] = perturbation
		
		# Append
		data.append(dataframe)

	# Convert
	dataframe_concat = pd.concat(data).drop(['t', 'B'], axis=1)

	# Round
	for col in ['logFC', 'AveExpr']:
		dataframe_concat[col] = [round(x, ndigits=3) for x in dataframe_concat[col]]

	# Format
	for col in ['P.Value', 'adj.P.Val']:
		dataframe_concat[col] = ['{:.2e}'.format(x) for x in dataframe_concat[col]]

	# Write
	dataframe_concat.to_csv(outfile, index=False, sep='\t')

#######################################################
#######################################################
########## S4. Enrichr
#######################################################
#######################################################

#############################################
########## 1. Get gene sets
#############################################

@follows(mkdir('s4-enrichr.dir/listids'))

@transform('s3-limma.dir/*/batch/*.txt',
		   regex(r'.*/(.*)-batch.txt'),
		   r's4-enrichr.dir/listids/\1-listids.json')

def getEnrichrIds(infile, outfile):

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

#############################################
########## 2. Add Enrichr Links
#############################################

@merge(getEnrichrIds,
	   's4-enrichr.dir/enrichr-results.xls')

def mergeIds(infiles, outfile):

	# Data
	data = []

	# Loop
	for infile in infiles:
		
		# Groups
		control, perturbation = os.path.basename(infile).split('-')[1].split('_vs_')
		
		# Read
		with open(infile) as openfile:
			ids = json.load(openfile)
			
		# Append
		data.append({
			'control': control,
			'perturbation': perturbation,
			'enrichr_up': 'http://amp.pharm.mssm.edu/Enrichr/enrich?dataset='+ids['up']['shortId'],
			'enrichr_down': 'http://amp.pharm.mssm.edu/Enrichr/enrich?dataset='+ids['down']['shortId']
		})

	# Convert
	dataframe = pd.DataFrame(data)[['control', 'perturbation', 'enrichr_up', 'enrichr_down']]

	# write
	dataframe.to_excel(outfile, index=False)


# #############################################
# ########## 2. Table JSON
# #############################################

# @transform(mergeDifferentialExpression,
# 		   suffix('.txt'),
# 		   '_table.json')

# def toTableJson(infile, outfile):

# 	# Read data
# 	limma_dataframe = pd.read_table(infile).drop(['t', 'B'], axis=1).set_index('group_name').rename(columns={'P.Value': 'pvalue', 'adj.P.Val': 'fdr'})

# 	# Round
# 	limma_dataframe['AveExpr'] = [round(x, ndigits=2) for x in limma_dataframe['AveExpr']]
# 	limma_dataframe['logFC'] = [round(x, ndigits=2) for x in limma_dataframe['logFC']]
# 	limma_dataframe['pvalue'] = ['{:.2e}'.format(x) for x in limma_dataframe['pvalue']]
# 	limma_dataframe['fdr'] = ['{:.2e}'.format(x) for x in limma_dataframe['fdr']]

# 	# Convert to dict
# 	limma_dict = {x: limma_dataframe.loc[x].to_dict(orient='records') for x in limma_dataframe.index.unique()}

# 	# Write
# 	with open(outfile, 'w') as openfile:
# 		openfile.write(json.dumps(limma_dict))

# #############################################
# ########## 3. Plot JSON
# #############################################

# @transform(mergeDifferentialExpression,
# 		   suffix('.txt'),
# 		   '_plot.json')

# def toPlotJson(infile, outfile):

# 	# Read data
# 	limma_dataframe = pd.read_table(infile).drop(['t', 'B'], axis=1).set_index('group_name')

# 	# Significance
# 	limma_dataframe['significant'] = [x < 0.05 for x in limma_dataframe['adj.P.Val']]

# 	# Initialize results
# 	results = {}

# 	# Loop through significance
# 	for condition in limma_dataframe.index.unique():

# 		# Get subset
# 		dataframe_subset = limma_dataframe.loc[condition]

# 		# Convert to dict
# 		results[condition] = {x: dataframe_subset[dataframe_subset['significant'] == x].to_dict(orient='list') for x in [True, False]}

# 	# Write
# 	with open(outfile, 'w') as openfile:
# 		openfile.write(json.dumps(results))

#######################################################
#######################################################
########## S. STAR
#######################################################
#######################################################

#############################################
########## 1. Align Reads
#############################################


##################################################
##################################################
########## Run pipeline
##################################################
##################################################
with open('pipeline/pipeline.png', 'wb') as openfile:
      pipeline_printout_graph(openfile, output_format='png')
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')
