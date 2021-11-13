library(WGCNA)

# ------------- patients with Cancer ------------------
# Choose a set of soft-thresholding powers 
powers = c(c(1:10), seq(from = 12, to=20, by=2)) 

# Call the WGCNA function 
sft = pickSoftThreshold(t(rna_expr_data_C), powerVector = powers, verbose = 5,blockSize=8504) 


par(mfrow = c(1,2))

# Analysis of the scale-free fit index for various soft-thresholding powers (β)
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
hist(k, main="Histogram of Connectivity \n distribution when β=14 \n (patients with Cancer)")

# The following log-log plot shows an R2 (the scale-free topology index)
# of 0.95, which means that this is a scale-free network.
scaleFreePlot(k, main="Check scale free topology \n (for patients with Cancer)\n")

# ------------- patients withOUT Cancer ------------------
# Choose a set of soft-thresholding powers 
powers = c(c(1:10), seq(from = 12, to=20, by=2)) 

# Call the WGCNA function 
sft = pickSoftThreshold(t(rna_expr_data_N), powerVector = powers, verbose = 5,blockSize=8504) 



# Analysis of the scale-free fit index for various soft-thresholding powers (β)
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

# network connectivities:
k=as.vector(apply(adjacency,2,sum))

# The following histogram shows the frequency distribution of the connectivity;
# we can see a large number of low connected genes, and a small number 
# of highly connected genes:
hist(k, main="Histogram of Connectivity \n distribution when β=12 \n (patients withOUT Cancer)")

# The following log-log plot shows an R2 (the scale-free topology index)
# of 0.94, which means that this is a scale-free network.
scaleFreePlot(k, main="Check scale free topology \n (for patients withOUT Cancer)\n")

