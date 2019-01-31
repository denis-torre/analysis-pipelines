#################################################################
#################################################################
############### Dubois RNA-seq Analysis - Python Support ############
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
import sys

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
import Support3 as S 

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####

#######################################################
#######################################################
########## S1. Limma
#######################################################
#######################################################

#############################################
########## 1. Limma
#############################################

def make_design_matrix(expression_dataframe, group_A, group_B, data):

	# Sample names
	group_A = [x.replace(':', '.').replace('-', '.') for x in group_A]
	group_B = [x.replace(':', '.').replace('-', '.') for x in group_B]
	expression_dataframe.columns = [x.replace(':', '.').replace('-', '.') for x in expression_dataframe.columns]

	# Collapse duplicate genes
	expression_dataframe = expression_dataframe.reset_index().groupby('index').sum()

	# Get expression dataframe
	if data == 'subset':
		expression_dataframe = expression_dataframe[group_A+group_B]

	# Create design dataframe
	sample_dict = {'A': group_A, 'B': group_B}
	design_dataframe = pd.DataFrame([{'index': x, 'A': int(x in group_A), 'B': int(x in group_B)} for x in expression_dataframe.columns]).set_index('index')

	# Return
	return {'expression': expression_dataframe, 'design': design_dataframe}

#############################################
########## 2. limma
#############################################

def limma(dataset, group_A, group_B, data='subset'):

	# Get design
	processed_data = make_design_matrix(dataset['rawdata'].copy(), group_A, group_B, data)

	# Add
	return pandas2ri.ri2py(r.limma(pandas2ri.py2ri(processed_data['expression']), pandas2ri.py2ri(processed_data['design']))).sort_values('logFC', ascending=False).set_index('gene_symbol')

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################