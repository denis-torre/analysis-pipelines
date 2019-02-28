#################################################################
#################################################################
############### Mamalis RNA-seq Analysis - R Support #################
#################################################################
#################################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. General support #####
# source('/Users/denis/Documents/Projects/scripts/Support.R')

##### 2. Other libraries #####

#######################################################
#######################################################
########## S1. 
#######################################################
#######################################################

#############################################
########## 1. 
#############################################

limma_rl <- function(infiles, outfile) {

    # Load
    require(edgeR)

    # Read data
    count_dataframe <- read.table(infiles[[1]], header=TRUE, row.names='gene_symbol')

    # Read metadata
    metadata_dataframe <- read.table(infiles[[2]], sep='\t', header=TRUE, row.names='sample')[colnames(count_dataframe),]
    metadata_dataframe$jsqcm <- as.character(metadata_dataframe$jsqcm)
    metadata_dataframe$timepoint <- make.names(metadata_dataframe$timepoint)

    # Design
    if (length(unique(metadata_dataframe$jsqcm)) > 1) {
        design <- model.matrix(~treatment+jsqcm+timepoint+0, data=metadata_dataframe);
    } else {
        design <- model.matrix(~treatment+timepoint+0, data=metadata_dataframe);
    }

    # Create DGEList object
    dge <- DGEList(counts=count_dataframe)

    # Filter genes
    keep <- filterByExpr(dge, design)
    dge <- dge[keep,]

    # Calculate normalization factors
    dge <- calcNormFactors(dge)

    # Run VOOM
    v <- voom(dge, plot=FALSE)

    # Fit linear model
    fit <- lmFit(v, design)

    # Make contrast matrix
    cont.matrix <- makeContrasts(de=treatmentRL-treatmentcontrol, levels=design)

    # Fit
    fit2 <- contrasts.fit(fit, cont.matrix)

    # Run DE
    fit2 <- eBayes(fit2)

    # Get results
    limma_dataframe <- topTable(fit2, adjust='BH', number=nrow(count_dataframe))
    limma_dataframe <- cbind(gene_symbol=rownames(limma_dataframe), limma_dataframe)

    # Write
    write.table(limma_dataframe, outfile, quote=FALSE, sep='\t')



}

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################