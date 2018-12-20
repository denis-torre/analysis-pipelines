#################################################################
#################################################################
############### McGrath Microarray Analysis - R Support #################
#################################################################
#################################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. General support #####
source('/Users/denis/Documents/Projects/scripts/Support.R')
library(limma)
library(Biobase)

##### 2. Other libraries #####

#######################################################
#######################################################
########## S1. IDAT Management
#######################################################
#######################################################

#############################################
########## 1. Read IDAT
#############################################

read_idat <- function(idat_file, bgx_file, outfile) {
    
    # Read Idat
    idat_data <- read.idat(idat_file, bgx_file)

    # Prepare dataframe
    idat_dataframe <- idat_data$genes
    idat_dataframe$expression <- idat_data$E

    # Write
    write.table(idat_dataframe, outfile, quote=FALSE, sep='\t')

}

#############################################
########## 2. Quantile Normalization
#############################################

quantile_normalization <- function(infile, outfile) {

    # Read expression
    expression_dataframe <- read.table(infile, sep='\t', header=TRUE, row.names='Symbol')

    # Normalize
    normalized_dataframe <- limma::normalizeQuantiles(expression_dataframe)

    # Write
    write.table(normalized_dataframe, outfile, quote=FALSE, sep='\t')

}

#############################################
########## 3. Run limma
#############################################

run_limma <- function(expression_file, metadata_file, outfile, comparison) {

    # Read data
    expression_dataframe <- read.table(expression_file)

    # Read metadata
    metadata_dataframe <- read.csv(metadata_file)

    # Create eset
    eset <- new("ExpressionSet", exprs=as.matrix(expression_dataframe))

    # Create design
    design <- model.matrix(~ 0 + metadata_dataframe$Condition)
    colnames(design) <- levels(factor(metadata_dataframe$Condition))

    # Fit
    fit <- lmFit(eset, design)

    # Make contrasts
    cmd <- paste("cont.matrix <- makeContrasts(", paste(comparison, collapse='-'), ", levels = design)", sep = '"')
    eval(parse(text = cmd))
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2)

    # Get table
    limma_dataframe <- topTable(fit2, adjust="BH", number=nrow(expression_dataframe))
    limma_dataframe <- cbind(gene_symbol=rownames(limma_dataframe), limma_dataframe)

    # Show results
    write.table(limma_dataframe, outfile, quote=FALSE, sep='\t', row.names=FALSE)
}




#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################