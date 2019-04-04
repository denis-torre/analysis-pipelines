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
options(stringsAsFactors=FALSE)

##### 2. Other libraries #####

#######################################################
#######################################################
########## S1. Differential Gene Expression
#######################################################
#######################################################

#############################################
########## 1. RL vs Normal
#############################################

limma_rl <- function(infiles, outfile) {

    # Load
    require(edgeR)

    # Read data
    count_dataframe <- read.table(infiles[[1]], header=TRUE, row.names='gene_symbol')

    # Read metadata
    metadata_dataframe <- read.table(infiles[[2]], sep='\t', header=TRUE, row.names='sample')

    # Filter
    metadata_dataframe <- metadata_dataframe[metadata_dataframe$cell_type == 'normal_fibroblast',]
    count_dataframe <- count_dataframe[,rownames(metadata_dataframe)]

    # Design
    design <- model.matrix(~treatment+cell_line+timepoint+intensity+0, data=metadata_dataframe)

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

#############################################
########## 2. Timepoint DGE
#############################################

limma_timepoint <- function(infiles, outfile) {

    # Load
    require(edgeR)

    # Get timepoints
    timepoint <- strsplit(outfile, '-')[[1]][3]

    # Read data
    count_dataframe <- read.table(infiles[[1]], header=TRUE, row.names='gene_symbol')

    # Read metadata
    metadata_dataframe <- read.table(infiles[[2]], sep='\t', header=TRUE, row.names='sample')

    # Filter
    metadata_dataframe <- metadata_dataframe[(metadata_dataframe$cell_type == 'normal_fibroblast') & (metadata_dataframe$timepoint == timepoint),]
    count_dataframe <- count_dataframe[,rownames(metadata_dataframe)]

    # Design
    design <- model.matrix(~treatment+cell_line+intensity+0, data=metadata_dataframe)

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
    write.table(limma_dataframe, outfile, quote=FALSE, sep='\t', row.names=FALSE)

}

#######################################################
#######################################################
########## S2. Timepoint Analysis
#######################################################
#######################################################

#############################################
########## 1. Cluster Genes
#############################################

cluster_genes <- function(infile, outfile) {

    # Read data
    t_dataframe <- read.table(infile, sep='\t', header=TRUE, row.names='gene_symbol', check.names=FALSE)

    # Drop NAs
    # t_dataframe <- na.omit(t_dataframe)

    # Get variance
    gene_var <- apply(t_dataframe, 1, var)

    # Filter
    top_genes <- names(gene_var)[gene_var > median(gene_var)]
    t_dataframe_filtered <- t_dataframe[top_genes,]

    # Scale
    scaled_data <- t(scale(t(t_dataframe_filtered)))
    
    # Get distance
    d <- dist(t_dataframe_filtered, method="euclidean")
    # d <- dist(scaled_data, method="euclidean")
    # d <- dist(t_dataframe_filtered, method="euclidean")
    # d <- as.dist(1-cor(t(t_dataframe_filtered)))
    # d <- as.dist(1-cor(t(t_dataframe_filtered)))

    # Cluster
    tree <- hclust(d, method="ward.D2")

    # Write
    save(t_dataframe_filtered, scaled_data, d, tree, file=outfile)

}

#############################################
########## 2. Cut Tree
#############################################

cut_tree <- function(infile, outfile) {

    # Load
    load(infile)

    # Cluster
    group_matrix <- as.data.frame(sapply(seq(1, 15), function(x) cutree(tree, k=x)))

    # Rename
    colnames(group_matrix) <- sapply(seq(1, length(group_matrix)),  function(x) paste0('k', x))

    # Write
    write.table(group_matrix, file=outfile, quote=FALSE, sep='\t')

}

#############################################
########## 3. Line Plots
#############################################

line_plots <- function(infiles, outfile) {

    print(infiles)

    # Load
    require(reshape2)
    require(ggplot2)

    # Read data
    t_dataframe <- read.table(infiles[1], sep='\t', header=TRUE, row.names='gene_symbol', check.names=FALSE)

    # Read groups
    group_dataframe <- read.table(infiles[2])

    # Pdf
    pdf(outfile)

    # Loop through groups
    # for (k in colnames(group_dataframe)) {
    k = 'k15'

    # Merge
    merged_dataframe <- merge(t_dataframe, group_dataframe, by='row.names')
    colnames(merged_dataframe)[1] <- 'gene_symbol'
    
    # Rename
    colnames(merged_dataframe)[colnames(merged_dataframe)==k] <- "group"
    
    # Subset
    plot_dataframe <- melt(merged_dataframe[,c('gene_symbol', '0h', '4h', '12h', '24h', 'group')], id.vars = c('gene_symbol', 'group'), variable.name = 'timepoint', value.name = 't')
    
    # Labeller
    group_sizes <- table(merged_dataframe$group)
    label_group <- function(group) {
        paste0('Cluster ', group, ' (', group_sizes[as.numeric(as.character(group))],' genes)')
    }

    # Centroid plots
    centroid_dataframe <- as.data.frame(t(sapply(unique(merged_dataframe$group), function(x) colMeans(merged_dataframe[merged_dataframe$group == x,c('0h', '4h', '12h', '24h')]))))
    centroid_plot_dataframe <- melt(cbind(group=rownames(centroid_dataframe), gene_symbol=rownames(centroid_dataframe), centroid_dataframe), id.vars=c('gene_symbol', 'group'), variable.name = 'timepoint', value.name = 't')

    # Group order
    plot_dataframe$group <- factor(plot_dataframe$group, levels=sort(unique(plot_dataframe$group)))
    centroid_plot_dataframe$group <- factor(centroid_plot_dataframe$group, levels=sort(unique(centroid_plot_dataframe$group)))
    
    # Plot
    gp <- ggplot(plot_dataframe, aes(x=timepoint, y=t, group=gene_symbol, color=as.factor(group))) +
        geom_line(alpha=0.1) + 
        geom_line(data=centroid_plot_dataframe, color='black') + 
        facet_wrap(~group, ncol=3, labeller=labeller(group=label_group)) +
        theme_bw() + 
        guides(color=FALSE) +
        labs(x='Timepoint', y='t-statistic (Control vs RL)', title=paste0('Clusters of genes co-varying across timepoints (',gsub('k', 'k = ',k),')')) +
        theme(plot.title = element_text(hjust = 0.5))
    print(gp)
        
    # }

    # Close device
    dev.off()


}

#############################################
########## 4. NbClust
#############################################

nbclust <- function(infile, outfile) {

    # Load
    load(infile)

    # Cluster
    nbc <- NbClust::NbClust(t_dataframe_filtered, method="ward.D2")

    # Save
    save(nbc, t_dataframe_filtered, file=outfile)
}

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################
