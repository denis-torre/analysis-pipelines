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
# from rpy2.robjects import r, pandas2ri
# pandas2ri.activate()

##### 2. Custom modules #####
# Pipeline running
# sys.path.append('/Users/denis/Documents/Projects/scripts')
# sys.path.append('pipeline/scripts')
sys.path.append('/Users/denis/Documents/Projects/jupyter-notebook/biojupies-plugins/library/core_scripts/normalize')
import normalize as N
sys.path.append('/Users/denis/Documents/Projects/jupyter-notebook/biojupies-plugins/library/core_scripts/signature')
import signature as S
# import Support3 as S
# import Dubois as P

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####

##### 2. R Connection #####
# r.source('pipeline/scripts/dubois.R')

#######################################################
#######################################################
########## S1. Align Files
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
########## 2. Align Reads
#############################################

@follows(mkdir('s1-kallisto.dir'))

@collate(glob.glob('rawdata/fastq_concat/*.fastq.gz'),
		 regex(r'rawdata/fastq_concat/(.*)_R._001.fastq.gz'),
         r's1-kallisto.dir/\1')

def alignReads(infiles, outfile):

	# Align
	print('Doing {}...'.format(outfile))
	try:
		os.system('kallisto quant -t 8 -i rawdata/kallisto/Homo_sapiens_index.idx -o {outfile} {infiles[0]} {infiles[1]}'.format(**locals()))
	except:
		print('Error doing {}.'.format(outfile))

#######################################################
#######################################################
########## S2. Merge Results
#######################################################
#######################################################

#############################################
########## 1. Merge Stats
#############################################

@follows(mkdir('s2-expression.dir'))

@merge(glob.glob('kallisto/*/run_info.json'),
       's2-expression.dir/run_info.txt')

def mergeInfo(infiles, outfile):

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
########## 2. Merge Expression
#############################################

def mergeJobs():

	# Get input files
	infiles = glob.glob('s1-kallisto.dir/*/abundance.tsv')

	# Loop through metrics
	for metric in ['est_counts']:

		# Get outfile
		outfile = 's2-expression.dir/dubois-{}.txt'.format(metric)

		# Yield
		yield [(infiles, 'rawdata/ensembl/mart_export.txt'), outfile]

@files(mergeJobs)

def mergeExpression(infiles, outfile):

	# Split infiles
	abundance_files, annotation_file = infiles

	# Read data
	dataframes = []
	for infile in abundance_files:
		if infile.split('/')[1].split('_')[0] in ['1', '2', '3', '4']:
			print("Skipping sample {}...".format(infile))
		else:
			sample_dataframe = pd.read_table(infile)
			sample_dataframe['sample'] = 'sample_'+infile.split('/')[1].split('_S')[0]
			dataframes.append(sample_dataframe)
	concatenated_dataframe = pd.concat(dataframes)

	# Read annotation
	annotation_dataframe = pd.read_table(annotation_file)

	# Get metric type
	metric = outfile.split('-')[-1][:-len('.txt')]

	# Cast data
	expression_dataframe = concatenated_dataframe.pivot_table(index='target_id', columns='sample', values=metric)

	# Merge
	merged_dataframe = expression_dataframe.merge(annotation_dataframe, left_index=True, right_on='Transcript stable ID version', how='inner')

	# Group counts to gene level
	if metric == 'est_counts':
		result_dataframe = merged_dataframe.drop(['Transcript stable ID version', 'Gene type'], axis=1).groupby('Gene name').sum().astype(int).rename_axis('gene_symbol')

	# Write
	result_dataframe.to_csv(outfile, sep='\t')

#############################################
########## 3. Process Metadata
#############################################

@files('rawdata/metadata/Dubois Lab Experimental Conditions.xlsx',
	   's2-expression.dir/dubois-metadata.txt')

def processMetadata(infile, outfile):

	# Read metadata
	print('Processing metadata...')
	metadata_dataframe = pd.read_excel(infile).iloc[:,:4].fillna('')

	# Sample labels
	metadata_dataframe['Sample'] = ['sample_'+str(x).replace('.', '_') for x in metadata_dataframe['Sample']]
	metadata_dataframe = metadata_dataframe.query('Diff not in ["hEB85", "hEB89", "hB90", "hEB91"]')

	# Old rowdata
	rowData_prev = []

	# Fill values
	for index, rowData in metadata_dataframe.iterrows():
		
		# Loop through columns
		for col in ['Condition', 'Description']:
			if not rowData[col]:
				rowData[col] = rowData_prev[col]
				
		# Set new rowdata
		rowData_prev = rowData.copy()

	# Write
	metadata_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 4. Normalize Expression
#############################################

def normalizeJobs():
	for normalization in ['logCPM']:
		yield ('s2-expression.dir/dubois-est_counts.txt', 's2-expression.dir/dubois-{}.txt'.format(normalization), normalization)

@follows(mergeExpression)

@files(normalizeJobs)

def normalizeData(infile, outfile, normalization):

	# Get expression
	count_dataframe = pd.read_table(infile, index_col='gene_symbol')

	# Normalize
	normalized_dataframe = getattr(N, normalization)({'rawdata': count_dataframe})

	# Write
	normalized_dataframe.to_csv(outfile, sep='\t')

#######################################################
#######################################################
########## S3. Mean Expression
#######################################################
#######################################################

#############################################
########## 1. Collapse Expression
#############################################

@follows(mkdir('s3-mean_expression.dir'))

@transform(normalizeData,
		   regex(r'.*/dubois-(.*).txt'),
		   add_inputs(processMetadata),
           r's3-mean_expression.dir/dubiois-grouped_\1.json')

def collapseExpression(infiles, outfile):

	# Get expression
	expression_dataframe = pd.read_table(infiles[0])

	# Read metadata
	metadata_dataframe = pd.read_table(infiles[1])

	# Merge
	merged_dataframe = pd.melt(expression_dataframe, id_vars='gene_symbol').merge(metadata_dataframe, left_on='variable', right_on='Sample').rename(columns={'variable': 'sample', 'value': 'logCPM'}).drop('Diff', axis=1)

	# Group
	grouped_dataframe = merged_dataframe.groupby(['gene_symbol', 'Condition'])['logCPM'].apply(list).rename('logCPM').to_frame()

	# Initialize dictionary
	expression_dict = {x: {} for x in expression_dataframe['gene_symbol']}

	# Loop
	for index, rowData in grouped_dataframe.iterrows():
		expression_dict[index[0]][index[1]] = rowData['logCPM']

	# Write
	with open(outfile, 'w') as openfile:
		openfile.write(json.dumps(expression_dict))

#############################################
########## 2. Get Genes
#############################################

@files(mergeExpression,
	   's3-mean_expression.dir/genes.json')

def getGenes(infile, outfile):

	# Get expression
	expression_dataframe = pd.read_table(infile[0])

	# Get results
	gene_dict = [{'gene_symbol': x} for x in expression_dataframe['gene_symbol']]

	# Write
	with open(outfile, 'w') as openfile:
		openfile.write(json.dumps(gene_dict))

#############################################
########## 3. Get Conditions
#############################################

@files(processMetadata,
	   's3-mean_expression.dir/conditions.json')

def getConditions(infile, outfile):

	# Get expression
	metadata_dataframe = pd.read_table(infile)

	# Get results
	gene_dict = [{'condition': x} for x in metadata_dataframe['Condition'].unique()]

	# Write
	with open(outfile, 'w') as openfile:
		openfile.write(json.dumps(gene_dict))

#######################################################
#######################################################
########## S4. limma
#######################################################
#######################################################

#############################################
########## 1. Get Differential Expression
#############################################

def differentialJobs():
	metadata_series = pd.read_table('s2-expression.dir/dubois-metadata.txt').groupby('Condition')['Sample'].apply(tuple)
	control_samples = metadata_series.pop('Control')
	for group_name, perturbed_samples in metadata_series.items():
		group_name = group_name.replace(' ', '')
		yield ('s2-expression.dir/dubois-est_counts.txt', 's4-limma.dir/{}-limma.txt'.format(group_name), group_name, control_samples, perturbed_samples)

@follows(mkdir('s4-limma.dir'))
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

#######################################################
#######################################################
########## S5. Merge results
#######################################################
#######################################################

#############################################
########## 1. Merge
#############################################

@follows(mkdir('s5-differential_expression.dir'))


@merge(runDifferentialExpression,
	   's5-differential_expression.dir/limma-merged.txt')

def mergeDifferentialExpression(infiles, outfile):

	# Concatenate
	merged_dataframe = pd.concat([pd.read_table(x) for x in infiles])

	# Write
	merged_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Table JSON
#############################################

@transform(mergeDifferentialExpression,
		   suffix('.txt'),
		   '_table.json')

def toTableJson(infile, outfile):

	# Read data
	limma_dataframe = pd.read_table(infile).drop(['t', 'B'], axis=1).set_index('group_name').rename(columns={'P.Value': 'pvalue', 'adj.P.Val': 'fdr'})

	# Round
	limma_dataframe['AveExpr'] = [round(x, ndigits=2) for x in limma_dataframe['AveExpr']]
	limma_dataframe['logFC'] = [round(x, ndigits=2) for x in limma_dataframe['logFC']]
	limma_dataframe['pvalue'] = ['{:.2e}'.format(x) for x in limma_dataframe['pvalue']]
	limma_dataframe['fdr'] = ['{:.2e}'.format(x) for x in limma_dataframe['fdr']]

	# Convert to dict
	limma_dict = {x: limma_dataframe.loc[x].to_dict(orient='records') for x in limma_dataframe.index.unique()}

	# Write
	with open(outfile, 'w') as openfile:
		openfile.write(json.dumps(limma_dict))

#############################################
########## 3. Plot JSON
#############################################

@transform(mergeDifferentialExpression,
		   suffix('.txt'),
		   '_plot.json')

def toPlotJson(infile, outfile):

	# Read data
	limma_dataframe = pd.read_table(infile).drop(['t', 'B'], axis=1).set_index('group_name')

	# Significance
	limma_dataframe['significant'] = [x < 0.05 for x in limma_dataframe['adj.P.Val']]

	# Initialize results
	results = {}

	# Loop through significance
	for condition in limma_dataframe.index.unique():

		# Get subset
		dataframe_subset = limma_dataframe.loc[condition]

		# Convert to dict
		results[condition] = {x: dataframe_subset[dataframe_subset['significant'] == x].to_dict(orient='list') for x in [True, False]}

	# Write
	with open(outfile, 'w') as openfile:
		openfile.write(json.dumps(results))

##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=3, verbose=1)
print('Done!')
