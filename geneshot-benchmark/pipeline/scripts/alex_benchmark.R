setwd('/Users/denis/Documents/Projects/analysis-pipelines/geneshot-benchmark')

# Read correlation
correlation_dataframe <- read.table('~/Downloads/correlation.tsv')
rownames(correlation_dataframe) <- correlation_dataframe[,1]
correlation_dataframe[,1] <- NULL
colnames(correlation_dataframe) <- rownames(correlation_dataframe)
correlation_dataframe[1:5,1:5]

# fraction
correlation_dataframe <- read.table('~/Desktop/fraction.txt')
correlation_dataframe <- correlation_dataframe[-1,]

# Change name
correlation_matrix <- sapply(correlation_dataframe, function(x){as.numeric(as.character(x))})
rownames(correlation_matrix) <- colnames(correlation_matrix)
correlation_matrix <- as.matrix(correlation_matrix)
correlation_matrix[1:5,1:5]

# Set NA
genes = rownames(correlation_matrix)      #expt=expression matrix, cor.score

#####
# Read GMT
gmt = readLines('libraries.dir/ARCHS4_Kinases_Coexp.txt')
genelists = list()                        #each element in list-> gene
for(ll in gmt){
  sp = unlist(strsplit(ll, "\t"))         #split strings of each line (into tabbed words)
  genesetname = sp[1]                     #1st string is the current geneset's name
  v = c()                                 #create an empty vector v, ends up being all of the genes in this gene set
  
  for(i in 3:length(sp)){                 #actual gene symbol starts at 3rd string(name of set, empty space for description, 1st gene)
    sp1 = unlist(strsplit(sp[i], ","))    #can read commas in a geneset names
    v = c(v, sp1[1])                      #dismiss other numbers only save gene symbols to v
  }
  
  ig = intersect(genes, v)                #unique genes of ref gmt & tissue corr expression matrix
  if(length(ig) > 1){                     #need at least 2 genes @ each geneset 
    genelists[[length(genelists)+1]] = ig #holds phenotypes & genes in list format, each set is a vector of its genes
    names(genelists)[length(genelists)] = gsub("_$", "", genesetname) #throwout underscore in genesets
  }
}

sgenes = rev(sort(table(unlist(genelists))))
sgenes = names(sgenes)[sgenes > 4]

reversegmt = list()

for(i in 1:length(genelists)){
  lname = names(genelists)[i]
  
  for(gene in genelists[[i]]){
    reversegmt[[gene]] = c(reversegmt[[gene]], lname)
  }
}
#####

#STEP II. Making a new expression matrix of only genes that are both in tissue-specific files & the reference gmt files
correlation_matrix[correlation_matrix > 0.99] = NA # remove self-self interaction
mean_matrix = matrix(0, nrow(correlation_matrix), length(genelists))
rownames(mean_matrix) = rownames(correlation_matrix)
colnames(mean_matrix) = names(genelists)

#STEP III. Very important-fast/efficient step for R,taking a subset of genes based on gene set & avging the corr all at once
for(p in names(genelists)){
  mean_matrix[, p] = rowMeans(correlation_matrix[,genelists[[p]]], na.rm=T)
}

mean_matrix[1:5,1:5]


########

library("pracma")
# define a function, scales to a range
scale_vector <- function(x, start, end) (x - min(x)) / max(x - min(x)) * (end - start) + start


aucgeneset = c() #gene function AUC vector
for (current_set in colnames(mean_matrix)) {
  setprob = mean_matrix[,current_set]
  v1 = rownames(mean_matrix) %in% genelists[[current_set]] 
  
  oo = rev(order(setprob))
  ov1 = v1[oo]
  
  cumulative = cumsum(ov1)
  scaled_y = scale_vector(cumulative, 0, 1)
  scaled_x= scale_vector(1:length(cumulative), 0, 1)
  
  aucgeneset[current_set] = trapz(scaled_x, scaled_y)
}

write.table(aucgeneset, '~/Desktop/auc_R.txt', sep='\t')
