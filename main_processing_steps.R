
# read the condition (tumor) TCGA-THCA dataset [TSS - 01] 
rna_expr_data_C <- read.table("./Data Extracted - Counts/TCGA-THCA/TCGA-THCA_rna_expr_data_C.txt", header=TRUE, sep="", check.names = F) 
# View(rna_expr_data_C)

rna_gene_data_C <- read.table("./Data Extracted - Counts/TCGA-THCA/TCGA-THCA_rna_genes_C.txt", header=FALSE, sep="", check.names = F) 
# View(rna_gene_data_C)

rna_patients_data_C <- read.table("./Data Extracted - Counts/TCGA-THCA/TCGA-THCA_rna_patients_C.txt", header=FALSE, sep="", check.names = F) 
# View(rna_patients_data_C)

# read the normal (natural) TCGA-THCA dataset [TSS - 11]
rna_expr_data_N <- read.table("./Data Extracted - Counts/TCGA-THCA/TCGA-THCA_rna_expr_data_N.txt", header=TRUE, sep="", check.names = F) 
# View(rna_expr_data_N)

rna_gene_data_N <- read.table("./Data Extracted - Counts/TCGA-THCA/TCGA-THCA_rna_genes_N.txt", header=FALSE, sep="", check.names = F) 
# View(rna_gene_data_N)

rna_patients_data_N <- read.table("./Data Extracted - Counts/TCGA-THCA/TCGA-THCA_rna_patients_N.txt", header=FALSE, sep="", check.names = F) 
# View(rna_patients_data_N)

# to be use (?) - for extra information, i think.. now its useless
clinical_data <- read.csv2("./Data Extracted - Counts/TCGA-THCA/TCGA-THCA_clinical_data.txt", header = FALSE, sep=",")
colnames(clinical_data) <- clinical_data[1,]; clinical_data <- clinical_data[-1,] 
# View(clinical_data)


# 1. Data - Intersection through the partecipants -------------------------

library(dplyr)
common_patients <- merge(rna_patients_data_C, rna_patients_data_N, by.x="V1", by.y="V1")

# subtract useless column between the two expression data.frames - check for multiple patterns and use the bool mask created to filter through the columns
# install.packages("stringr"), if you haven't the library
library(stringr)

mask_tuned_patientsC <- str_detect(as.vector(colnames(rna_expr_data_C)), paste(as.vector(unlist(common_patients)), collapse = '|'))
rna_expr_data_C <- rna_expr_data_C[, mask_tuned_patientsC]

mask_tuned_patientsN <- str_detect(as.vector(colnames(rna_expr_data_N)), paste(as.vector(unlist(common_patients)), collapse = '|')) # useless part
rna_expr_data_N <- rna_expr_data_N[, mask_tuned_patientsN] # useless part

# cleaning the rows with ANY zero values through each of the two data.frames
rna_expr_data_C <- rna_expr_data_C[apply(rna_expr_data_C, 1, function(x) all(x!=0)),]
rna_expr_data_N <- rna_expr_data_N[apply(rna_expr_data_N, 1, function(x) all(x!=0)),]

# intersect the rows (genes) through the two data.frames
rna_expr_data_C <- subset(rna_expr_data_C, rownames(rna_expr_data_C) %in% rownames(rna_expr_data_N))
rna_expr_data_N <- subset(rna_expr_data_N, rownames(rna_expr_data_N) %in% rownames(rna_expr_data_C))

# this is only a demonstration... we haven't duplicated partecipants, remove the dashes to try it
# dummy_split <- as.data.frame(t(apply(common_patients, 1, function(k){
#   sym <- unlist(strsplit(k, split = "-"))
#   return(c(sym[1], sym[2], sym[3]))
# }))); colnames(dummy_split) <- c("A", "B", "C"); dummy_split[duplicated(dummy_split$C),]


# 2. Differentially Expressed Genes (DEGs) --------------------------------

# It's important to handle Raw Counts of the number of mRNA producted by n genes, install DESeq2

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")

library(DESeq2)

# sort the column names for each dataset before merging the columns
rna_expr_data_C <- rna_expr_data_C[,order(colnames(rna_expr_data_C))]
rna_expr_data_N <- rna_expr_data_N[,order(colnames(rna_expr_data_N))]

# change a little bit the colnames between these two datasets
addInfoData <- function(dt, info){
  colInfo <- c()

  for(name in colnames(dt)){
    colInfo <- c(colInfo, paste(name, info, sep=""))  
  }
  
  colnames(dt) <- colInfo
  return(dt)
}

rna_expr_data_C <- addInfoData(rna_expr_data_C, "_tumor")
rna_expr_data_N <- addInfoData(rna_expr_data_N, "_normal")

# merges!
full_dt <- cbind(rna_expr_data_N, rna_expr_data_C) # sorted and aggregated! :)
# View(full_dt)

# create a dataset with additional information
PartecipantsCondition <- data.frame(partecipants = colnames(full_dt),
                                    condition = NA) 
PartecipantsCondition <- as.data.frame(t(apply(PartecipantsCondition, 1, function(x){
    if(grepl("tumor", x["partecipants"]))
      return(c(x["partecipants"], condition = "tumor"))
    else
      return(c(x["partecipants"], condition = "normal"))
})));


# creating the Object of our merged matrix to be "DESQued"
DESQ_dt <- DESeqDataSetFromMatrix(countData = full_dt,
                                  colData = PartecipantsCondition,
                                  design= ~ condition)

# remove the genes where there are very few reads (lowly expressed) and we reduce at the same time the memory size of the actual dataset, we increase also the future transformation and testing functions
# Here we perform a minimal pre-filtering to keep only rows that have at least 10 reads total. Note that more strict filtering to increase power is automatically applied via independent filtering on the mean of normalized counts within the results function.
keep <- rowSums(counts(DESQ_dt)) >= 10 # this could be a first threshold for us... more strict less genes
DESQ_dt <- DESQ_dt[keep,]

# main 
DESQ_dtOutput <- DESeq(DESQ_dt)

# get the results, obtaining the correction FDR through different tests, its important to note the adjpvalue!
res <- results(DESQ_dtOutput, alpha = 0.05, lfcThreshold = 1.2) # this applies a threshold of alpha = 0.5, to make correction of different tests applied, we apply some corrections optimizing the number of genes which will have an adjusted p value below a given FDR cutoff, alpha = 0.05 and not 0.1
res

# better to show the final summary of all
summary(res) 
# LFC = Log Fold Change
# LFC > 0 means the number of genes up-regulated and viceversa (down-regulated < 0)

# plotting the differentially expressed gene ordered from the lower at the larger value of adjustedPvalue to see the genes properly.
resOrdered <- res[order(res$padj), ]

# versus conditions
resultsNames(DESQ_dtOutput)

### more extra info: if a gene is higher expressed to be a normal condition this method gives the intensity of this direction, the logic explained is also considered in the opposite meaning
### FLC > 0 is more expressed in NORMAL condition
### FLC < 0 is more expressed in TUMOR condition
### this application tells the intensity of this direction!!

# the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down. 
plotMA(res, ylim = c(-8, 8))


# Result Volcano plot ------------------------------------------------------------

topT <- as.data.frame(res)

#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~p~value)))
with(subset(topT, padj<=0.05 & log2FoldChange>=1.2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))

with(subset(topT, padj<=0.05 & log2FoldChange<= -1.2), points(log2FoldChange, -log10(padj), pch=20, col="green", cex=0.5))

#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-1.2, col="black", lty=4, lwd=2.0)
abline(v=1.2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$padj<=0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)

# 3. Co-expression networks -----------------------------------------------

# extract the filtered DEGs
res_save <- subset(resOrdered, padj <= 0.05)
final_dt <- as.data.frame(res_save)

# extract the genes 
genes_passed <- rownames(final_dt)

# filter the two datasets expressed in tumor and normal 
rna_expr_data_C <- rna_expr_data_C[genes_passed, ]
rna_expr_data_N <- rna_expr_data_N[genes_passed, ]

# log-transform FPKM data using log2(x+1) of each count
rna_expr_data_C <- log2(rna_expr_data_C+1)
rna_expr_data_N <- log2(rna_expr_data_N+1)

# Extra Part - Find the good threshold / Emulate Fiscon lecture ------------------------------------

library(igraph)
createNet <- function(dt, dt2, x){
  # create the correlation datasets for plotting the network for each graph
  co_net_corr_dataC <- cor(t(dt), method = "pearson")
  co_net_corr_dataN <- cor(t(dt2), method = "pearson")

  tsh <- x
  co_net_corr_dataC <- ifelse(co_net_corr_dataC <= abs(tsh) & co_net_corr_dataC >= -abs(tsh), 1, 0)
  co_net_corr_dataN <- ifelse(co_net_corr_dataN <= abs(tsh) & co_net_corr_dataN >= -abs(tsh), 1, 0)

  gC <- graph_from_adjacency_matrix(co_net_corr_dataC, diag = FALSE)
  gN <- graph_from_adjacency_matrix(co_net_corr_dataN, diag = FALSE)
  
  avg_dens <- (edge_density(gN, loop = FALSE) + edge_density(gC, loop = FALSE))/2
  return(avg_dens)
}

possibletsh <- seq(0, 1, length.out = 20) # the behaviour is symmetric
densities <- unlist(lapply(possibletsh, function(x){
  return(createNet(rna_expr_data_C, rna_expr_data_N, x))
}))

plot(possibletsh, densities, col = "blue", type = "l", lwd = 3, xlab = "Pearson Correlation", ylab = "Fraction of nodes (Density)", main = "Optimal Threshold Between the two graphs")
abline(h = 0.8, lty=2) # abline(h = densities, lty=2)
abline(v = 0.495, lty=2) # abline(v = possibletsh, lty = 3)
points(x = 0.495, y = 0.8, pch = 20, col = "red", cex = 1.5) # this is our preferable thresholding

# Ended the Extra Part ----------------------------------------------------

# create the correlation datasets for plotting the network for each graph
co_net_corr_dataC <- cor(t(rna_expr_data_C), method = "pearson")
co_net_corr_dataN <- cor(t(rna_expr_data_N), method = "pearson")

# binary masks
tsh <- 0.495
co_net_corrBinary_dataC <- ifelse(co_net_corr_dataC <= abs(tsh) & co_net_corr_dataC >= -abs(tsh), 1, 0)
co_net_corrBinary_dataN <- ifelse(co_net_corr_dataN <= abs(tsh) & co_net_corr_dataN >= -abs(tsh), 1, 0)

# create the graph
gC <- graph_from_adjacency_matrix(co_net_corrBinary_dataC, diag = FALSE)
plot(gC, vertex.size=5, edge.curverd=.1, arrow.size=.1, vertex.color = "red", main = "Co-expression network in TUMOR",
     arrow.width=.1, edge.arrow.size=.1, layout= layout.kamada.kawai, vertex.label = NA) # , vertex.label.dist = .8 and .y, vertex.label.cex=1

gN <- graph_from_adjacency_matrix( co_net_corrBinary_dataN, diag = FALSE)
plot(gN, vertex.size=5, edge.curverd=.1, arrow.size=.1, vertex.color = "green", main = "Co-expression network in NORMAL",
     arrow.width=.1, edge.arrow.size=.1, layout= layout.kamada.kawai, vertex.label = NA) # , vertex.label.dist = .8 and .y, vertex.label.cex=1

# degree distribution of the graphs
par(mfrow=c(1,2))
hist(degree.distribution(gC), main = "Degree distribution in TUMOR condition", col = "red", xlab = "Degree Distribution")
hist(degree.distribution(gN), main = "Degree distribution in NORMAL condition", col = "green", xlab = "Degree Distribution")

# extract the 5% of HUBS, in their conditions
hubs_C <- sort(degree(gC, v = V(gC), mode = "all"), decreasing = TRUE) # normalized TRUE
hubs_C <- hubs_C[1:floor(0.05 * length(hubs_C))] 
hubs_N <- sort(degree(gN, v = V(gN), mode = "all"), decreasing = TRUE) # normalized TRUE
hubs_N <- hubs_N[1:floor(0.05 * length(hubs_N))] 

# find the common hubs
namesHUBS_C <- names(hubs_C)
namesHUBS_N <- names(hubs_N)

intersect(namesHUBS_C, namesHUBS_N) # Common HUBS!

# 4. Differential Co-expressed Network ------------------------------------




# 5. OPTIONAL TASKS -------------------------------------------------------

# the points provided by the professor
# adding the type of edge (red or blue) in the context of the pearson correlation > 0 or viceversa
