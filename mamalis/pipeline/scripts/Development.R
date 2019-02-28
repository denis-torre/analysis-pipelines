setwd('/Users/maayanlab/Documents/Projects/analysis-pipelines/mamalis')

infiles <- c('s2-expression.dir/kallisto/kallisto-counts.txt', 'rawdata/metadata/dubois-metadata.txt')

require(edgeR)

# Read metadata
metadata_dataframe <- read.table(infiles[2], sep='\t', header=TRUE, row.names='Sample')
metadata_dataframe$Condition <- make.names(metadata_dataframe$Condition)

# Read data
count_dataframe <- read.table(infiles[1], header=TRUE, row.names='gene_symbol')[,rownames(metadata_dataframe)]

# Get groups
groups <- make.names(c(' + FA + delta, 7d', ' + delta, 7d phased'))

# Design
Diff <- factor(metadata_dataframe$Diff)
Condition <- factor(metadata_dataframe$Condition)
design <- model.matrix(~Condition+Diff+0)

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
eval(parse(text=paste0('cont.matrix <- makeContrasts(de=Condition',groups[1],'-Condition',groups[2],', levels=design)')))

# Fit
fit2 <- contrasts.fit(fit, cont.matrix)

# Run DE
fit2 <- eBayes(fit2)

# Get results
limma_dataframe <- topTable(fit2, adjust='BH', number=nrow(count_dataframe))
limma_dataframe <- cbind(gene_symbol=rownames(limma_dataframe), limma_dataframe)

# Write
write.table(limma_dataframe, outfile, quote=FALSE, sep='\t')


