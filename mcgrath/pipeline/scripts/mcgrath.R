#################################################################
#################################################################
############### McGrath Microarray Analysis - R Support #################
#################################################################
#################################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. General support #####
# source('/Users/denis/Documents/Projects/scripts/Support.R')
library(limma)

##### 2. Other libraries #####

#######################################################
#######################################################
########## S1. IDAT Management
#######################################################
#######################################################

#############################################
########## 1. Read IDAT
#############################################

read_idat <- function(idat_files, bgx_file, outfile) {

    # Read
    library(limma)
    
    # Read data
    data <- read.idat(idatfiles = idat_files, bgxfile = bgx_file)

    # Detection
    data$other$Detection <- detectionPValues(data)

    # Save
    save(data, file=outfile)

}

#############################################
########## 2. Normalize data
#############################################

normalize_data <- function(infile, outfile) {

    # Load
    load(infile)

    # Normalize
    data_normalized <- neqc(data)

    # Add variance and average
    data_normalized$genes$probe_variance <- apply(data_normalized$E, 1, var)
    data_normalized$genes$probe_mean <- apply(data_normalized$E, 1, mean)

    # Save
    save(data_normalized, file=outfile)

}

#############################################
########## 3. Normalize data
#############################################

extract_data <- function(infile, outfile) {

    # Load
    load(infile)

    # Get expression
    expression_dataframe <- data_normalized$E

    # Fix colnames
    colnames(expression_dataframe) <- sapply(colnames(expression_dataframe), function(x) paste0(strsplit(x, '/')[[1]][3], '.idat'))

    # Merge data
    merged_dataframe <- merge(data_normalized$genes, expression_dataframe, by.x='Array_Address_Id', by.y='row.names')

    # Write
    write.table(merged_dataframe, file = outfile, quote=FALSE, sep='\t')

}

#############################################
########## 3. Run limma
#############################################

run_limma <- function(expression_file, metadata_file, outfile, comparison) {

    # Load
    load(expression_file)

    # Read metadata
    metadata_dataframe <- read.csv(metadata_file, stringsAsFactors = FALSE)
    rownames(metadata_dataframe) <- sapply(metadata_dataframe$IDATfile, function(x) paste0('rawdata/7196780027/', strsplit(x, '.', fixed = TRUE)[[1]][1]))

    # Design
    design <- model.matrix(~Condition+0, data=metadata_dataframe)

    # Fit linear model
    fit <- lmFit(data_normalized, design)

    # Make contrast matrix
    eval(parse(text=paste0('cont.matrix <- makeContrasts(de=Condition',comparison[2],'-Condition',comparison[1],', levels=design)')))

    # Fit
    fit2 <- contrasts.fit(fit, cont.matrix)

    # Run DE
    fit2 <- eBayes(fit2)

    # Get results
    limma_dataframe <- topTable(fit2, adjust='BH', number=nrow(data_normalized$E))

    # Write
    write.table(limma_dataframe, file = outfile, quote=FALSE, sep='\t')

}




#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################