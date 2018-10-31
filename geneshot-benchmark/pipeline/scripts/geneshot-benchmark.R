#################################################################
#################################################################
############### Geneshot Benchmarking - R Support #################
#################################################################
#################################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. General support #####
source('/Users/denis/Documents/Projects/scripts/Support.R')
require(ggplot2)

##### 2. Other libraries #####

#######################################################
#######################################################
########## S1. Plot AUC
#######################################################
#######################################################

#############################################
########## 1. Plot AUC Distribution
#############################################

plot_auc <- function(infile, outfile, plot_type, na) {

    # Print
    print(paste0('Plotting ', outfile, '...'))

	# Get AUC scores
    auc_dataframe <- read.csv2(infile, sep='\t')

    # Filter
    auc_dataframe <- auc_dataframe[auc_dataframe$normalization %in% c('correlation_v2', 'zscore', 'generif_overlap_zscore', 'autorif_overlap_zscore', 'random'),]

    # Convert variables
    auc_dataframe$Similarity <- as.character(auc_dataframe$normalization)
    auc_dataframe$auc <- as.numeric(as.character(auc_dataframe$auc))

    # Replace
    auc_dataframe[auc_dataframe == 'correlation_v2'] <- 'ARCHS4 Coexpression'
    auc_dataframe[auc_dataframe == 'zscore'] <- 'Enrichr Co-occurrence'
    auc_dataframe[auc_dataframe == 'generif_overlap_zscore'] <- 'GeneRIF'
    auc_dataframe[auc_dataframe == 'autorif_overlap_zscore'] <- 'AutoRIF'
    auc_dataframe[auc_dataframe == 'random'] <- 'Random'

    # NA
    if (na) {
        auc_dataframe[is.na(auc_dataframe)] <- 0
    }

    # Plot
    if (plot_type == 'violin') {
        gp <- ggplot(auc_dataframe, aes(x=Similarity, y=auc, fill=Similarity)) + geom_violin() + facet_wrap(library~., scales='free_y', ncol=5) + theme_minimal() + scale_y_continuous(lim=c(0,1)) + xlab('Similiarity') + ylab('AUC') + theme(axis.text.x=element_blank())
    } else if (plot_type == 'density') {
        gp <- ggplot(auc_dataframe, aes(x=auc, color=Similarity)) + geom_density() + facet_wrap(library~., scales='free_y', ncol=5) + theme_minimal() + scale_x_continuous(lim=c(0,1)) + xlab('AUC') + ylab('Density')
    }

    # Write
    ggsave(outfile, gp, width=8, height=4, scale=1.5)
}


#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################