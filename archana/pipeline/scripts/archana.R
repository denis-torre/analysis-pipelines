#################################################################
#################################################################
############### Archana RNA-seq Analysis - R Support #################
#################################################################
#################################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. General support #####
source('/Users/denis/Documents/Projects/scripts/Support.R')
require(limma)
require(edgeR)

##### 2. Other libraries #####

#######################################################
#######################################################
########## S1. Limma
#######################################################
#######################################################

#############################################
########## 1. Limma
#############################################

run_limma <- function(infiles, group, outfile) {

    # Read data
    count_dataframe <- read.table(infiles[[1]], sep='\t', header=TRUE, row.names='gene_symbol')
    metadata_dataframe <- read.table(infiles[[2]], header=TRUE, row.names='sample')

    # Filter
    metadata_dataframe_subset <- metadata_dataframe[metadata_dataframe$group==group,]
    count_dataframe_subset <- count_dataframe[,rownames(metadata_dataframe_subset)]

    # Design
    patient_nr <- factor(metadata_dataframe_subset$patient_nr)
    treatment <- factor(metadata_dataframe_subset$treatment)
    design <- model.matrix(~patient_nr+treatment)

	# Create DGEList object
	dge <- DGEList(counts=count_dataframe_subset)

	# Calculate normalization factors
	dge <- calcNormFactors(dge)

	# Run VOOM
	v <- voom(dge, plot=FALSE)

	# Fit linear model
	fit <- lmFit(v, design)

	# Run DE
	fit <- eBayes(fit)

	# Get results
	limma_dataframe <- topTable(fit, adjust='BH', number=nrow(count_dataframe_subset), coef='treatmentPre_PPI')
	limma_dataframe <- cbind(gene_symbol=rownames(limma_dataframe), limma_dataframe)

    # Write
    write.table(limma_dataframe, outfile, sep='\t', quote=FALSE, row.names=FALSE)

}

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################