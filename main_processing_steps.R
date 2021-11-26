
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
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~adj.pvalue)))
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

# create the correlation datasets for plotting the network for each graph
co_net_corr_dataC <- cor(t(rna_expr_data_C), method = "pearson")
co_net_corr_dataN <- cor(t(rna_expr_data_N), method = "pearson")

# plot the distribution of the correlations to pick naive the good 
distroRho_pearsC <- co_net_corr_dataC[upper.tri(co_net_corr_dataC)]
distroRho_pearsN <- co_net_corr_dataN[upper.tri(co_net_corr_dataN)]

par(mfrow = c(1, 2))
hist(distroRho_pearsC, main = "Distribution of Correlations in TUMOR genes", col = "red", xlab = "Rho Correlation", breaks = 50)
hist(distroRho_pearsN, main = "Distribution of Correlations in NORMAL gene", col = "green", xlab = "Rho Correlation", breaks = 50)

# these are probably the candidates of having a good threshold
quantile(abs(distroRho_pearsC))
quantile(abs(distroRho_pearsN))

## Approximately is good to choose 0.55 instead 0.7!
# Extra Part - Find the good threshold / Emulate Fiscon lecture ------------------------------------

library(igraph)

fractionNodes <- function(graph){
  v_graph <- length(V(graph))
  component <- components(graph)
  ind <- which(component$membership == which.max(component$csize))
  LCC <- induced_subgraph(graph , V(graph)[ind])
  v_LCC <- length(V(LCC))
  
  frac_node_LCC <- v_LCC / v_graph
  return(frac_node_LCC)
}

OptimalThresholding <- function(dt, dt2, x){
  # create the correlation datasets for plotting the network for each graph
  co_net_corr_dataC <- cor(t(dt), method = "spearman")
  co_net_corr_dataN <- cor(t(dt2), method = "spearman")
  
  tsh <- x
  co_net_corr_dataC <- ifelse(co_net_corr_dataC <= -abs(tsh) | co_net_corr_dataC >= abs(tsh), 1, 0)
  co_net_corr_dataN <- ifelse(co_net_corr_dataN <= -abs(tsh) | co_net_corr_dataN >= abs(tsh), 1, 0)
  
  gC <- graph_from_adjacency_matrix(co_net_corr_dataC, diag = FALSE)
  gN <- graph_from_adjacency_matrix(co_net_corr_dataN, diag = FALSE)
  
  fracNodes_C <- fractionNodes(gC)
  fracNodes_N <- fractionNodes(gN)
  
  return(mean(fracNodes_C, fracNodes_N))
}

possibletsh <- seq(0, 1, length.out = 20) # the behaviour is symmetric
densities <- unlist(lapply(possibletsh, function(x){
  return(OptimalThresholding(rna_expr_data_C, rna_expr_data_N, x))
}))

plot(possibletsh, densities, col = "blue", type = "l", lwd = 3, xlab = "Threshold of Correlation", ylab = "Fraction of nodes", main = "Hard Threshold Choice")
rect(xleft = 0.5, ybottom = 0.0, xright = 0.65, ytop = 1.0, density = 5, border = "red", lty = 2, lwd = 1)

# we consider a naive threshold to consider a well connected network like a scale-free with few hubs, its however a hard thresholding in this case
abline(h = 0.97, lty=2)
abline(v = 0.55, lty=2)
points(x = 0.55, y = 0.97, pch = 20, col = "red", cex = 1.5) # this is our preferable naive-hard thresholding

# Ended the Extra Part ----------------------------------------------------

# binary masks
library(igraph)
tsh <- 0.55 # 0.7
co_net_corrBinary_dataC <- ifelse(co_net_corr_dataC <= -abs(tsh) | co_net_corr_dataC >= abs(tsh), 1, 0)
co_net_corrBinary_dataN <- ifelse(co_net_corr_dataN <= -abs(tsh) | co_net_corr_dataN >= abs(tsh), 1, 0)

par(mfrow = c(1, 2))
hist(distroRho_pearsC, main = "Distribution of Correlations in TUMOR genes", col = "red", xlab = "Rho Correlation", breaks = 50)
abline(v = 0.55, lty=2, col = "white")
points(x = 0.55, y = 0.0, pch = 20, col = "yellow", cex = 1.5) # this is our preferable naive-hard thresholding
hist(distroRho_pearsN, main = "Distribution of Correlations in NORMAL gene", col = "green", xlab = "Rho Correlation", breaks = 50)
abline(v = 0.55, lty=2, col = "white")
points(x = 0.55, y = 0.0, pch = 20, col = "yellow", cex = 1.5) # this is our preferable naive-hard thresholding

# create the graph
par(mfrow=c(1,1))
gC <- graph_from_adjacency_matrix(co_net_corrBinary_dataC, diag = FALSE)
plot(gC, vertex.size=5, edge.curverd=.1, arrow.size=.1, vertex.color = "red", main = "Co-expression network in TUMOR",
     arrow.width=.1, edge.arrow.size=.1, layout= layout.kamada.kawai, vertex.label = NA) # , vertex.label.dist = .8 and .y, vertex.label.cex=1


gN <- graph_from_adjacency_matrix( co_net_corrBinary_dataN, diag = FALSE)
plot(gN, vertex.size=5, edge.curverd=.1, arrow.size=.1, vertex.color = "green", main = "Co-expression network in NORMAL",
     arrow.width=.1, edge.arrow.size=.1, layout= layout.kamada.kawai, vertex.label = NA) # , vertex.label.dist = .8 and .y, vertex.label.cex=1

# degree distribution of the graphs, extract the 5% of HUBS, in their conditions
par(mfrow=c(1,2))
dgC <- degree(gC)
dgN <- degree(gN)
hubs_C <- sort(degree(gC, v = V(gC), mode = "all"), decreasing = TRUE) # normalized TRUE
hubs_C <- hubs_C[1:floor(0.05 * length(hubs_C))] 
hubs_N <- sort(degree(gN, v = V(gN), mode = "all"), decreasing = TRUE) # normalized TRUE
hubs_N <- hubs_N[1:floor(0.05 * length(hubs_N))] 

# plot the distribution highliting the 5% of hubs
hist(dgC[dgC != 0], main = "Degree distribution in TUMOR condition", col = "red", xlab = "Degree Distribution", breaks = 50)
hist(hubs_C, add = T, col = "orchid")
abline(v = tail(hubs_C, n=1), lty = 2)
hist(dgN[dgN != 0], main = "Degree distribution in NORMAL condition", col = "green", xlab = "Degree Distribution", breaks = 50)
hist(hubs_N, add = T, col = "orchid")
abline(v = tail(hubs_N, n=1), lty = 2)

# find the hubs and do something interesting
namesHUBS_C <- names(hubs_C)
namesHUBS_N <- names(hubs_N)

namesHUBS_C_df <- data.frame(gene = namesHUBS_C, degree = hubs_C, row.names = 1:length(namesHUBS_C)); namesHUBS_C_df
namesHUBS_N_df <- data.frame(gene = namesHUBS_N, degree = hubs_N, row.names = 1:length(namesHUBS_N)); namesHUBS_N_df

library(ggplot2)
library(cowplot)

histcol <- c("degree_N" = "green2", "degree" = "red")
hub_gene_C_plot <- ggplot(data = namesHUBS_C_df, aes(x = reorder(gene, degree), y = degree)) +
                      geom_col(aes(x = reorder(gene, degree), y = degree, fill = "degree")) +
                      theme_light() + 
                      theme(text = element_text(size=10),
                            axis.text.x = element_text(angle = 45,hjust = 1),
                            axis.title.x = element_text(face="bold"),
                            axis.title.y = element_text(face="bold"),
                            plot.margin=unit(c(t = 1.5, r = 1, b = 1.5, l = 1.2), "cm"),
                            plot.title = element_text(face = "bold")) +
                      labs( y = "Degree Index", x = "Gene") +
                      scale_fill_manual(values= histcol) +
                      ggtitle("Hub Gene Identification in Tumor Condition") +
                      guides(fill=FALSE)

hub_gene_N_plot <- ggplot(data = namesHUBS_N_df, aes(x = reorder(gene, degree), y = degree)) +
                      geom_col(aes(x = reorder(gene, degree), y = degree, fill = "degree_N")) +
                      theme_light() + 
                      theme(text = element_text(size=10),
                            axis.text.x = element_text(angle = 45,hjust = 1),
                            axis.title.x = element_text(face="bold"),
                            axis.title.y = element_text(face="bold"),
                            plot.margin=unit(c(t = 1.5, r = 1, b = 1.5, l = 1.2), "cm"),
                            plot.title = element_text(face = "bold")) +
                      labs( y = "Degree Index", x = "Gene") +
                      scale_fill_manual(values = histcol) +
                      ggtitle("Hub Gene Identification in Normal Condition") +
                      guides(fill=FALSE)

plot_grid(hub_gene_C_plot, hub_gene_N_plot, labels = NULL)

intersect(namesHUBS_C, namesHUBS_N) # Common HUBS!

## ENSG00000182580 <=> EPHB3
## A receptor gene involved on bring much signals and get two-way communication answers with other gene...
## Check the closeness to confirm!
sort(closeness(gC), decreasing = T)[1] # EUREKA! :D

## It's more considered as an increasing receptor respect the other cells...
## In fact...
ifelse(res["ENSG00000182580", "log2FoldChange"] >= 1.2, "UP-REGULATED", "DOWN-REGULATED") # DOUBLE EUREKA! :D

### Other Extra-Stats for HUB Genes -----------------------------------------------------

# # Up-regulated 
upgenes <- final_dt[final_dt$log2FoldChange >= 1.2, ]
hub_upgenes <- na.omit(upgenes[namesHUBS_C, ])

# # Down-regulated 
downgenes <- final_dt[final_dt$log2FoldChange <= 1.2, ]
hub_downgenes <- na.omit(downgenes[namesHUBS_C, ])

# # Up-regulated 
upgenes_N <- final_dt[final_dt$log2FoldChange >= 1.2, ]
hub_upgenes_N <- na.omit(upgenes_N[namesHUBS_N, ])

# # Down-regulated 
downgenes_N <- final_dt[final_dt$log2FoldChange <= 1.2, ]
hub_downgenes_N <- na.omit(downgenes_N[namesHUBS_N, ])

# 4. Differential Co-expressed Network ------------------------------------

par(mfrow=c(1,1)) #eventually reset the plot visualization

## Find the Optimal Z-treshold --------------------

OptimalThresholdingZ <- function(dt, dt2, x){
  # create the correlation datasets for plotting the network for each graph
  co_net_corr_dataC <- cor(t(dt), method = "pearson")
  co_net_corr_dataN <- cor(t(dt2), method = "pearson")
  
  # Application Z-Fisher Transform
  ZcoDataC <- log((1+co_net_corr_dataC)/(1-co_net_corr_dataC))/2
  ZcoDataN <- log((1+co_net_corr_dataN)/(1-co_net_corr_dataN))/2
  
  # Applying z-scores
  ZcoData <- (ZcoDataC-ZcoDataN)/sqrt((1/(nrow(rna_expr_data_C)-3)) + (1/(nrow(rna_expr_data_N)-3)))
  
  # Applying Z-tsh
  tshZ <- x
  ZcoData <- ifelse(ZcoData <= -abs(tshZ) | ZcoData >= abs(tshZ), 1, 0)
  
  gZ <- graph_from_adjacency_matrix(ZcoData, diag = FALSE)
  
  fracNodes_C <- fractionNodes(gZ)
  
  return(fracNodes_C)
}

possibletshZ <- seq(1, 15, by = 1) # the behaviour is symmetric
densitiesZ <- unlist(lapply(possibletshZ, function(x){
  return(OptimalThresholdingZ(rna_expr_data_C, rna_expr_data_N, x))
}))

plot(possibletshZ, densitiesZ, col = "blue", type = "l", lwd = 3, xlab = "Z-Threshold", ylab = "Fraction of nodes", main = "Z-Threshold Choice")
rect(xleft = 8.5, ybottom = 0.0, xright = 12, ytop = 1.0, density = 5, border = "red", lty = 2, lwd = 1)

# we consider a naive threshold to consider a well connected network like a scale-free with few hubs, its however a hard Zthresholding in this case
abline(h = 0.9965, lty=2)
abline(v = 11, lty=2)
points(x = 11, y = 0.9965, pch = 20, col = "red", cex = 1.5) # this is our preferable naive-hard Z-thresholding

## ------------------------------------------------

# create the correlation datasets for plotting the network for each graph
co_net_corr_dataC <- cor(t(rna_expr_data_C), method = "pearson")
co_net_corr_dataN <- cor(t(rna_expr_data_N), method = "pearson")

# Application Z-Fisher Transform
ZcoDataC <- log((1+co_net_corr_dataC)/(1-co_net_corr_dataC))/2
ZcoDataN <- log((1+co_net_corr_dataN)/(1-co_net_corr_dataN))/2

# Applying z-scores
ZcoData <- (ZcoDataC-ZcoDataN)/sqrt((1/(nrow(rna_expr_data_C)-3)) + (1/(nrow(rna_expr_data_N)-3)))

# Z-Threshold it
tshZ <- 11 # 11 is a general good Z-threshold 
ZcoData <- ifelse(ZcoData <= -abs(tshZ) | ZcoData >= abs(tshZ), 1, 0)

# Get the graph and plot it
gZcoData <- graph_from_adjacency_matrix( ZcoData, diag = FALSE)
plot(gZcoData, vertex.size=5, edge.curverd=.1, arrow.size=.1, vertex.color = "green", main = "Differential Co-expression network in TUMOR vs NORMAL",
     arrow.width=.1, edge.arrow.size=.1, layout= layout.auto, vertex.label = NA)

# dgZcoData degree distribution
dgZcoData <- degree(gZcoData)
hist(dgZcoData[dgZcoData != 0], main = "Degree distribution in the Differential Co-expression", col = "red", xlab = "Degree Distribution", breaks = 50)

# extract the 5% of HUBS, in their conditions
hubs_Z <- sort(degree(gZcoData, v = V(gZcoData), mode = "all"), decreasing = TRUE) # normalized TRUE
hubs_Z <- hubs_Z[1:floor(0.05 * length(hubs_Z))] 

# Comparing hubs in TUMORS and Z-Tum vs. Z-Norm
namesHUBS_Z <- names(hubs_Z)
namesHUBS_C <- names(hubs_C)
namesHUBS_N <- names(hubs_N)

hubs_commonZC <- intersect(namesHUBS_C, namesHUBS_Z) # Common HUBS!
print(hubs_commonZC)

# generalize, not always used... is to highlight the subnetwork of each hub! You can try changing the hubs_commonZC to see the results
for(hub_comm in hubs_commonZC){
  # save the genes connected for each hub in common
  subnodesHub <- names(na.omit(ZcoData[hub_comm, ZcoData[hub_comm, ] == 1]))
  
  subgraph_one_hub <- induced.subgraph(graph = gZcoData, vids = c(subnodesHub, hub_comm))
  
  V(subgraph_one_hub)$color <- ifelse(V(subgraph_one_hub)$name == hub_comm, "red", "green")
  plot(subgraph_one_hub, vertex.size=5, edge.curverd=.1, arrow.size=.1, vertex.color = V(subgraph_one_hub)$color, main = paste("Subgraph of Hub", hub_comm, "in Differential Co-expression", sep = " "),
       arrow.width=.1, edge.arrow.size=.1, layout= layout_on_sphere, vertex.label = NA)
}


# 5. OPTIONAL TASKS -------------------------------------------------------

# solution of the optional points 

# 1.

library(WGCNA)

# ------------- patients with Cancer ------------------
# Choose a set of soft-thresholding powers 
powers = c(c(1:10), seq(from = 12, to=20, by=2)) 

# Call the WGCNA function 
sft = pickSoftThreshold(t(rna_expr_data_C), powerVector = powers, verbose = 5,blockSize=8504) 


par(mfrow = c(1,2))

# Analysis of the scale-free fit index for various soft-thresholding powers (??)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = "Scale independence \n (patients with Cancer)");
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=0.9,col="red")
abline(h=0.85,col="red") # this line corresponds to using an R^2 cut-off of 0.85


#Analysis of the mean connectivity for various soft-thresholding powers  
plot(sft$fitIndices[-(1:2),1], sft$fitIndices[-(1:2),5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = "Mean connectivity \n (patients with Cancer)") 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9, col="red")


#it seems that 14 could be a power soft threshold to ensure a scale-free network, so:
soft_th = 14

#calculate the corresponding adjacency matrix:
adjacency=abs(cor(t(rna_expr_data_C),use="p"))^soft_th

# network connectivities:
k=as.vector(apply(adjacency,2,sum))

# The following histogram shows the frequency distribution of the connectivity;
# we can see a large number of low connected genes, and a small number 
# of highly connected genes:
hist(k, main="Histogram of Connectivity \n distribution when ??=14 \n (patients with Cancer)")

# The following log-log plot shows an R2 (the scale-free topology index)
# of 0.95, which means that this is a scale-free network.
scaleFreePlot(k, main="Check scale free topology \n (for patients with Cancer)\n")

# ------------- patients withOUT Cancer ------------------
# Choose a set of soft-thresholding powers 
powers = c(c(1:10), seq(from = 12, to=20, by=2)) 

# Call the WGCNA function 
sft = pickSoftThreshold(t(rna_expr_data_N), powerVector = powers, verbose = 5,blockSize=8504) 



# Analysis of the scale-free fit index for various soft-thresholding powers (??)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = "Scale independence \n (patients withOUT Cancer)");
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=0.9,col="red")
abline(h=0.85,col="red") # this line corresponds to using an R^2 cut-off of 0.85


#Analysis of the mean connectivity for various soft-thresholding powers  
plot(sft$fitIndices[-(1:2),1], sft$fitIndices[-(1:2),5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = "Mean connectivity \n (patients withOUT Cancer)") 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9, col="red")



#it seems that 12 could be a power soft threshold to ensure a scale-free network, so:
soft_th = 12

#calculate the corresponding adjacency matrix:
adjacency=abs(cor(t(rna_expr_data_N),use="p"))^soft_th

# network connectivity:
k=as.vector(apply(adjacency,2,sum))

# The following histogram shows the frequency distribution of the connectivity;
# we can see a large number of low connected genes, and a small number 
# of highly connected genes:
hist(k, main="Histogram of Connectivity \n distribution when ??=12 \n (patients withOUT Cancer)")

# The following log-log plot shows an R2 (the scale-free topology index)
# of 0.94, which means that this is a scale-free network.
scaleFreePlot(k, main="Check scale free topology \n (for patients withOUT Cancer)\n")


# 2.

# Check the overlapping between the 5% of the nodes with highest
# CI values and the "Degree"-based hubs

comparisonCIwithDegreeHUBS <- function(graph, TypegHUBS) {
  
  # calculates the centralities
  CI_index_between <- sort(betweenness(graph), decreasing = T) # Betweenness centrality
  CI_index_closeness <- sort(closeness(graph), decreasing = T) # Closeness centrality
  CI_index_eigen <- sort(eigen_centrality(graph)$vector, decreasing = T) # Eigen centrality
  
  # extract the hubs
  hubs_CI_between <- CI_index_between[1:floor(0.05 * length(CI_index_between))] # Find the Hubs (top 5%)
  namesHUBS_CI_between <- names(hubs_CI_between) # and their names
  
  hubs_CI_closeness <- CI_index_closeness[1:floor(0.05 * length(CI_index_closeness))] # Find the Hubs (top 5%)
  namesHUBS_CI_closeness <- names(hubs_CI_between) # and their names
  
  hubs_CI_eigen <- CI_index_eigen[1:floor(0.05 * length(CI_index_eigen))] # Find the Hubs (top 5%)
  namesHUBS_CI_eigen <- names(hubs_CI_eigen) # and their names
  
  # put all in a list
  localCentralities <- list(
    betweeness = intersect(TypegHUBS, namesHUBS_CI_between),
    closeness = intersect(TypegHUBS, namesHUBS_CI_closeness),
    eigen = intersect(TypegHUBS, namesHUBS_CI_eigen)
  )
  
  return(localCentralities) #Let's see which centralityIndex-based hubs are also "Degree"-based hubs
}

plotCommonHUBS <- function(matrix_type, graph_type, LocalCentralities_type){ # to be complete!!
  
  for(hubs_type in names(LocalCentralities_type)){
    hub_common <- LocalCentralities_type[[hubs_type]]
    
    V(graph_type)$color <- ifelse(V(graph_type)$name %in% hub_common, "red", "green")
    
    plot(graph_type, vertex.size=5, edge.curverd=.1, arrow.size=.1, vertex.color = V(graph_type)$color, main = paste(hubs_type, "=", length(hub_common), "hubs in common", sep = " "),
         arrow.width=.1, edge.arrow.size=.1, layout= layout_on_sphere, vertex.label = NA)
  }
  
}

LocalCentralities_C <- comparisonCIwithDegreeHUBS(gC, namesHUBS_C); LocalCentralities_C # Centralities
plotCommonHUBS(co_net_corrBinary_dataC, gC, LocalCentralities_C) # Plot the common hubs!
LocalCentralities_N <- comparisonCIwithDegreeHUBS(gN, namesHUBS_N); LocalCentralities_N # Centralities
plotCommonHUBS(co_net_corrBinary_dataN, gN, LocalCentralities_N) # Plot the common hubs!

####### 3. Extra part in which we use two different correlation method

## Print the values to check which value of correlation is better
par(mfrow=c(1,2))
hist(rna_expr_data_C[upper.tri(rna_expr_data_C)], main = "Distribution of the Tumor data", col = "red", xlab = "Value", breaks = 50)
rug(rna_expr_data_C[upper.tri(rna_expr_data_C)])
hist(rna_expr_data_N[upper.tri(rna_expr_data_N)], main = "Distribution of the Normal data", col = "green", xlab = "Value", breaks = 50)
rug(rna_expr_data_N[upper.tri(rna_expr_data_N)])

par(mfrow=c(2,2))
method_correlation <- c("pearson", "spearman")
for (i in  method_correlation) {
  # create the correlation datasets for plotting the network for each graph
  co_net_corr_dataC <- cor(t(rna_expr_data_C), method = i)
  co_net_corr_dataN <- cor(t(rna_expr_data_N), method = i)
  
  # binary masks
  tsh <- 0.55
  co_net_corrBinary_dataC <- ifelse(co_net_corr_dataC <= -abs(tsh) | co_net_corr_dataC >= abs(tsh), 1, 0)
  co_net_corrBinary_dataN <- ifelse(co_net_corr_dataN <= -abs(tsh) | co_net_corr_dataN >= abs(tsh), 1, 0)
  
  # create the graph
  gC <- graph_from_adjacency_matrix(co_net_corrBinary_dataC, diag = FALSE)
  plot(gN, vertex.size=5, edge.curverd=.1, arrow.size=.1, vertex.color = "RED", main = "Co-expression network TUMOR",
       arrow.width=.1, edge.arrow.size=.1, layout= layout.kamada.kawai, vertex.label = NA)
  
  
  gN <- graph_from_adjacency_matrix( co_net_corrBinary_dataN, diag = FALSE)
  plot(gN, vertex.size=5, edge.curverd=.1, arrow.size=.1, vertex.color = "green", main = "Co-expression network NORMAL",
       arrow.width=.1, edge.arrow.size=.1, layout= layout.kamada.kawai, vertex.label = NA) 
  
  # degree distribution of the graphs
  dgC <- degree(gC)
  dgN <- degree(gN)
  hist(dgC[dgC != 0], main = paste(i, "Corr. TUMOR", sep = " ") , col = "red", xlab = "Degree Distribution")
  hist(dgN[dgN != 0], main = paste(i, "Corr. NORMAL", sep = " "), col = "green", xlab = "Degree Distribution")
  
  # extract the 5% of HUBS, in their conditions
  hubs_C <- sort(degree(gC, v = V(gC), mode = "all"), decreasing = TRUE) # normalized TRUE
  hubs_C <- hubs_C[1:floor(0.05 * length(hubs_C))] 
  hubs_N <- sort(degree(gN, v = V(gN), mode = "all"), decreasing = TRUE) # normalized TRUE
  hubs_N <- hubs_N[1:floor(0.05 * length(hubs_N))] 
  
  nam <- paste("hubs_C_", i, sep = "")
  assign(nam, names(hubs_C))
  mamt <- paste("hubs_N_", i, sep = "")
  assign(mamt, names(hubs_N))
}

intersect(hubs_N_pearson, hubs_N_spearman) # 10 hub genes in common
intersect(hubs_C_pearson, hubs_C_spearman) # 19 hub genes in common

# TO DO: further works (eventually)
# 1. adding the type of edge (red or blue) in the context of the pearson correlation > 0 or viceversa