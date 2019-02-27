require(limma)

# Read metadata
setwd('/Users/denis/Documents/Projects/analysis-pipelines/mcgrath')
metadata_dataframe <- read.csv('rawdata/conditions_correct.csv')
metadata_dataframe

# Read data
setwd('rawdata/7196780027/')
idat_files <- Sys.glob('*.idat')
bgx_file <- 'HumanHT-12_V4_0_R2_15002873_B.bgx'
data <- read.idat(idatfiles = idat_files, bgxfile = bgx_file)

# Stats
data$other$Detection <- detectionPValues(data)
propexpr(data)

# Normalize
data_normalized <- neqc(data)

# Plot
sums <- sort(colMeans(data$E))
par(mar=c(5,10,4,2))
barplot(sums, horiz = TRUE, las=1, names.arg = sapply(names(sums), function(x) strsplit(x, '/')[[1]][3]))
