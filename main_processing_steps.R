
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


# 3. Co-expression networks -----------------------------------------------


# 4. Differential Co-expressed Network ------------------------------------


# 5. OPTIONAL TASKS -------------------------------------------------------


