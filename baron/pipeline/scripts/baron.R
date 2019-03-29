#################################################################
#################################################################
############### Baron RNA-seq Analysis - R Support #################
#################################################################
#################################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. General support #####

##### 2. Other libraries #####

#######################################################
#######################################################
########## S1. 
#######################################################
#######################################################

#############################################
########## 1. 
#############################################

run_limma <- function(infile, outfile, comparison) {
    # Load
    require(edgeR)

    # Read counts
    count_dataframe <- read.table(infile, header=TRUE, row.names='gene_symbol')
    count_dataframe$P3.1 <- NULL

    # Get design
    group <- sapply(colnames(count_dataframe), function(x) strsplit(x, '.', fixed=TRUE)[[1]][1])
    design <- model.matrix(~group+0)

    # Create DGEList object
    dge <- DGEList(counts=count_dataframe)
    keep <- filterByExpr(dge, design)

    # Calculate normalization factors
    dge <- calcNormFactors(dge[keep,])

    # Run VOOM
    v <- voom(dge, plot=FALSE)

    # Fit linear model
    fit <- lmFit(v, design)

    # Make contrast matrix
    cmd <- paste("cont.matrix <- makeContrasts(", paste(rev(comparison), collapse='-'), ", levels = design)", sep = '"')
    eval(parse(text = cmd))

    # Fit
    fit2 <- contrasts.fit(fit, cont.matrix)

    # Run DE
    fit2 <- eBayes(fit2)

    # Get results
    limma_dataframe <- topTable(fit2, adjust='BH', number=nrow(count_dataframe))
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