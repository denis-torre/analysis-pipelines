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

plot_auc <- function(infile, outfile) {

	# Get AUC scores
    auc_dataframe <- read.csv2(infile, sep='\t')

    # Filter
    auc_dataframe <- auc_dataframe[auc_dataframe$normalization %in% c('correlation', 'zscore', 'generif_overlap_zscore', 'autorif_overlap_zscore'),]

    # Convert to numeric
    auc_dataframe$auc <- as.numeric(as.character(auc_dataframe$auc))

    # Plot
    gp <- ggplot(auc_dataframe, aes(x=auc, color=normalization)) + geom_density() + facet_wrap(library~., scales='free_y', ncol=3) + theme_minimal() + scale_x_continuous(lim=c(0,1)) + xlab('AUC') + ylab('Density')

    # Write
    ggsave(outfile, gp, width=6, height=4, scale=1.5)
}


#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################