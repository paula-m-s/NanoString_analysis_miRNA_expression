# PIPELINE INFO ----
########################################################################################################################.
########################################################################################################################.
###
###   Nanostring data analysis - miRNA
###  
###   Authors: Paula Morales-Sanchez
###   Pipeline support: 
###     http://www.bioconductor.org/packages/release/bioc/vignettes/NanoTube/inst/doc/NanoTube.html 
###     https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8138885/
###       https://github.com/bhattacharya-a-bt/CBCS_normalization/
###   Date started: April 2023
###   Last date modified: August 2023
###   R version 
###   Script: Raw data and QC from NanoString miRNA RCC files. Versions: NS_H_miR_v3a and NS_H_miR_v3b
###
###   Public data source:  
###   Platform: NanoString
########################################################################################################################.
########################################################################################################################.


# 1. Sample cleaning and probe filtering ----

## 1.1. Cleaning pData dataframe ----

# For subsequent pipelines. Keep only the samples that passed all requirements
pData_f <- pData[pData$count_flag!="Flag" & pData$Sample_Title !="FVPTC-E518K-1-BIS" & pData$Sample_Diagnosis!="normal" & pData$Sample_Names != "FVPTC_WT_3" & pData$Sample_Names != "MNG_WT_3" & pData$Sample_Names != "MNG_D1_4" & pData$Sample_Names != "FTC_D8_1" & pData$Sample_Names != "PDTC_D1_5" & pData$Sample_Names != "FVPTC_D1_8" & pData$Sample_Names != "PDTC_D1_4" , ]
pData_f <- pData[pData$count_flag!="Flag" & pData$Sample_Title !="FVPTC-E518K-1-BIS" & pData$Sample_Diagnosis!="normal", ]
#pData_f <- pData[pData$count_flag!="Flag" & pData$Sample_Title!="FVPTC-E518K-1-BIS", ]
#pData_f <- pData[pData$count_flag!="Flag" & pData$Sample_Title!="FVPTC-E518K-1-BIS" & pData$Sample_Diagnosis!="normal" & pData$Sample_Diagnosis!="normal"  & pData$Sample_Diagnosis!="cPTC" & pData$Tumor_Type_PDTC!="TC_DGCR8_mutated" & pData$Tumor_Type_PDTC!="MNG_DGCR8_mutated" & pData$Sample_Diagnosis!="MNG", ]


## 1.2. Filtering probes ----

# Keep only endogenous for DE analysis
fData_f<- fData[fData$Class == "Endogenous",] # remove housekeeping (differentially expressed between; high variability across samples)
rownames(raw_expression) <- raw_expression$Gene
raw_f <- raw_expression[fData_f$Gene,pData_f$SampleID]


# Remove calibration factors attached to the name. Drop probes with different sensitivities between vr.a and vr.b rlf
version_a <- readRcc(files.RCC[1])$Code_Summary[, 2] # Check pData and get one of version a. Retrive names in column 2 of code summary
version_b <- readRcc(files.RCC[40])$Code_Summary[, 2] # Check pData and get one of version b. Retrive names in column 2 of code summary
mirna_diff_sen <- gsub("\\|.*", "", version_a[!(version_a %in% version_b)]) 
raw_f <- raw_f[!(rownames(raw_f) %in% mirna_diff_sen),]
table(!(rownames(raw_f) %in% mirna_diff_sen))


# Filter less expressed mirnas: at least 10 cpm in at least 70% of the samples
table(selectGenes(raw_f, min.count = 10, N = 0.7))
keep <- (selectGenes(raw_f, min.count = 10, N = 0.7))

## 1.3. Final count miRNA dataframe ----  
raw_f <- raw_f[keep,]
colnames(raw_f) <- pData_f$Sample_Names

# Filter mirna classes dataframe
fData_f = fData_f[rownames(raw_f),]
rownames(pData_f) <- pData_f$Sample_Names

# Order dataframes by group
pData_f <- pData_f[order(pData_f$Tumor_Type_Mut),]
raw_f <- raw_f[,pData_f$Sample_Names]

raw_f <- as.data.frame(raw_f)

pseudo_counts <- as.matrix(log2(raw_f+ 1))
head(pseudo_counts)


df_raw = pseudo_counts - apply(pseudo_counts,1,median)
df_raw <- melt(df_raw, id = rownames(df_raw))
names(df_raw)[1:2]<- c("id", "sample")
df_raw$method <- rep("Raw counts", nrow(df_raw))  
head(df_raw)

pal <- brewer.pal(8, "Dark2")

### Scatter Plot of Median and Variance per Sample

# Calculate median per sample
median_per_sample <- apply(pseudo_counts, 2, median)

# Calculate variance per sample
variance_per_sample <- apply(pseudo_counts, 2, var)

# Create a new dataframe with median and variance per sample
df_summary <- data.frame(median_per_sample, variance_per_sample)

# Add sample names as a column
df_summary$sample <- pData_f$Tumor_Type_Mut

# Convert 'sample' column to factor
df_summary$sample <- as.factor(df_summary$sample)

# Plot scatter plot using ggplot2
per_med_var_raw <- ggplot(df_summary, aes(x = median_per_sample, y = variance_per_sample, color = sample)) +
  geom_point(size = 4) +
  xlab("Per-Sample Median") +
  ylab("Per-Sample Variance") +
  labs(title = "Scatter Plot of Median and Variance per Sample \n Raw pseudocounts ") +
  theme_minimal() + geom_text_repel(aes(label = rownames(df_summary)), size = 3) + 
  theme(plot.title = element_text(hjust = 0.5))+
  #scale_shape_manual(values = c(15:19,7:9)) + 
  #scale_color_manual(values = pal) +  theme_bw(base_size = 18)  + 
  theme(axis.text.x = element_text(size = 25, colour = "black"), axis.text.y = element_text(size = 25, colour = "black"),axis.title.y = element_text(size = 18, colour = "black"), legend.text = element_text(size = 20, colour = "black") )
per_med_var_raw

#PCA
resPCA <- mixOmics::pca(t(pseudo_counts), center = T, scale = T, ncomp = 10)
plot(resPCA)
colores <- brewer.pal(3, "Dark2")

# Add legend

PCA <- prcomp(t(pseudo_counts), scale. = T, center = T)
pca.var <- PCA$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.var.per

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Group = pData_f$Tumor_Type_Mut, Version = pData_f$Version)
#dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Group = pData_f$Tumor_Type_Mut, Version = paste(pData_f$Version,pData_f$Source_Year, sep ="_"))
pal = brewer.pal(12, "Dark2")

pca_raw <- ggplot(dataGG, aes(PC1, PC2, label =rownames(dataGG))) +
  geom_point(aes(shape = Group, colour = Version), size=5) +
  ggtitle("PCA raw pseudo counts") +
  xlab(paste0("PC1, VarExp: ", pca.var.per[1], "%")) +
  ylab(paste0("PC2, VarExp: ", pca.var.per[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  #coord_fixed(ratio = 1) +
  scale_shape_manual(values = c(12:18)) + 
  #scale_color_manual(values = pal) +  
  theme_bw(base_size = 18)  + 
  theme(axis.text.x = element_text(size = 25, colour = "black"), axis.text.y = element_text(size = 25, colour = "black"),axis.title.y = element_text(size = 18, colour = "black"), legend.text = element_text(size = 20, colour = "black") ) + geom_text_repel(aes(label = rownames(dataGG)), size = 3)  #+  theme(aspect.ratio=3/3)
pca_raw

# 2. Normalization ----
## 2.0 Generate edgeR objects: ----  
group <- as.factor(pData_f$Tumor_Type_Mut)
miRNA_matrix <- DGEList(counts=raw_f, group = group)

## 2.1. Total counts ----
pseudo_TC <- log2(cpm(miRNA_matrix) + 1)

#PCA
resPCA <- mixOmics::pca(t(pseudo_TC), center = T, scale = T, ncomp = 10)
plot(resPCA)
colores <- brewer.pal(3, "Dark2")
plotIndiv(resPCA, group = pData_f$Version, col.per.group = colores[1:2], pch = c(1,2), legend = T, legend.title = "Group", ind.names = T, title = "PCA pseudo counts") 
# Add legend

PCA <- prcomp(t(pseudo_TC), scale. = T, center = T)
pca.var <- PCA$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Group = pData_f$Experiment_Run_date, Version = pData_f$Version) 
#dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Group = pData_f$Tumor_Type_Mut, Version = paste(pData_f$Version,pData_f$Source_Year, sep ="_"))
pal = brewer.pal(8, "Dark2")

pca_TC <- ggplot(dataGG, aes(PC1, PC2, label =rownames(dataGG))) +
  geom_point(aes(shape = Group, colour = Version), size=5) +
  ggtitle("PCA TC") +
  xlab(paste0("PC1, VarExp: ", pca.var.per[1], "%")) +
  ylab(paste0("PC2, VarExp: ", pca.var.per[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  #coord_fixed(ratio = 1) +
  scale_shape_manual(values = c(15,19,17,19,18))+ 
  #scale_color_manual(values = pal) +  
  theme_bw(base_size = 18)  + 
  theme(axis.text.x = element_text(size = 25, colour = "black"), axis.text.y = element_text(size = 25, colour = "black"),axis.title.y = element_text(size = 18, colour = "black"), legend.text = element_text(size = 20, colour = "black") ) + geom_text_repel(aes(label = rownames(dataGG)), size = 3)# +  theme(aspect.ratio=3/3)
pca_TC

# Scatter Plot of Median and Variance per Sample

# Calculate median per sample
median_per_sample <- apply(pseudo_TC, 2, median)

# Calculate variance per sample
variance_per_sample <- apply(pseudo_TC, 2, var)

# Create a new dataframe with median and variance per sample
df_summary <- data.frame(median_per_sample, variance_per_sample)

# Add sample names as a column
df_summary$sample <- pData_f$Tumor_Type_PDTC

# Convert 'sample' column to factor
df_summary$sample <- as.factor(df_summary$sample)

# Plot scatter plot using ggplot2
per_med_var_TC <- ggplot(df_summary, aes(x = median_per_sample, y = variance_per_sample, color = sample)) +
  geom_point(size = 3) +
  xlab("Per-Sample Median") +
  ylab("Per-Sample Variance") +
  labs(title = "Scatter Plot of Median and Variance per Sample \n log2 counts per million ") +
  theme_minimal() + geom_text_repel(aes(label = rownames(df_summary)), size = 3) + 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_shape_manual(values = c(15:19,7:9)) + 
  scale_color_manual(values = pal) +  theme_bw(base_size = 18)  + 
  theme(axis.text.x = element_text(size = 25, colour = "black"), axis.text.y = element_text(size = 25, colour = "black"),axis.title.y = element_text(size = 18, colour = "black"), legend.text = element_text(size = 20, colour = "black") )




df_TC = pseudo_TC - apply(pseudo_TC,1,median)
df_TC <- melt(df_TC, id = rownames(df_TC))
names(df_TC)[1:2] <- c ("id", "sample")
df_TC$method <- rep("Edger_TC", nrow(df_TC))

## 2.2. Upper quartile ----
miRNA_matrix_uq <- calcNormFactors(miRNA_matrix, method = "upperquartile")
miRNA_matrix_uq$samples

test_normcount <- sweep(miRNA_matrix_uq$counts, 2,
                        miRNA_matrix_uq$samples$lib.size*miRNA_matrix_uq$samples$norm.factors / 10^6,
                        "/")
range(as.vector(test_normcount - cpm(miRNA_matrix_uq)))

pseudo_UQ <- cpm(miRNA_matrix_uq,log = T,prior.count = 1)


#PCA
resPCA <- mixOmics::pca(t(pseudo_UQ), center = T, scale = T, ncomp = 10)
plot(resPCA)
colores <- brewer.pal(3, "Dark2")
plotIndiv(resPCA, group = pData_f$Version, col.per.group = colores[1:2], pch = c(1,2), legend = T, legend.title = "Group", ind.names = T, title = "PCA pseudo counts") 
# Add legend

PCA <- prcomp(t(pseudo_UQ), scale. = T, center = T)
pca.var <- PCA$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Group = pData_f$Tumor_Type_PDTC, Version = pData_f$Tumor_Type_PDTC) 
#dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Group = pData_f$Tumor_Type_Mut, Version = paste(pData_f$Version,pData_f$Source_Year, sep ="_"))
pal = brewer.pal(8, "Dark2")

pca_UQ <- ggplot(dataGG, aes(PC1, PC2, label =rownames(dataGG))) +
  geom_point(aes(shape = Group, colour = Version), size=5) +
  ggtitle("PCA UQ") +
  xlab(paste0("PC1, VarExp: ", pca.var.per[1], "%")) +
  ylab(paste0("PC2, VarExp: ", pca.var.per[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  #coord_fixed(ratio = 1) +
  scale_shape_manual(values =  c(15,19,17,19,18,20,21)) + 
  #scale_color_manual(values = pal) +  
  theme_bw(base_size = 18)  + 
  theme(axis.text.x = element_text(size = 25, colour = "black"), axis.text.y = element_text(size = 25, colour = "black"),axis.title.y = element_text(size = 18, colour = "black"), legend.text = element_text(size = 20, colour = "black") ) + geom_text_repel(aes(label = rownames(dataGG)), size = 3) #+  theme(aspect.ratio=3/3)
pca_UQ

# Scatter Plot of Median and Variance per Sample

# Calculate median per sample
median_per_sample <- apply(pseudo_UQ, 2, median)

# Calculate variance per sample
variance_per_sample <- apply(pseudo_UQ, 2, var)

# Create a new dataframe with median and variance per sample
df_summary <- data.frame(median_per_sample, variance_per_sample)

# Add sample names as a column
df_summary$sample <- pData_f$Group2

# Convert 'sample' column to factor
df_summary$sample <- as.factor(df_summary$sample)

# Plot scatter plot using ggplot2
per_med_var_UQ <- ggplot(df_summary, aes(x = median_per_sample, y = variance_per_sample, color = sample)) +
  geom_point(size = 3) +
  xlab("Per-Sample Median") +
  ylab("Per-Sample Variance") +
  labs(title = "Scatter Plot of Median and Variance per Sample \n EdgeR UQ counts") +
  theme_minimal() + geom_text_repel(aes(label = rownames(df_summary)), size = 3) + 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_shape_manual(values = c(15:19,7:9)) + 
  scale_color_manual(values = pal) +  theme_bw(base_size = 18)  + 
  theme(axis.text.x = element_text(size = 25, colour = "black"), axis.text.y = element_text(size = 25, colour = "black"),axis.title.y = element_text(size = 18, colour = "black"), legend.text = element_text(size = 20, colour = "black") )



df_UQ = pseudo_UQ - apply(pseudo_UQ,1,median)
df_UQ <- melt(df_UQ, id = rownames(df_UQ))
names(df_UQ)[1:2] <- c ("id", "sample")
df_UQ$method <- rep("Edger_UQ", nrow(df_UQ))

## 2.3. TMM ----
miRNA_matrix_tmm <- calcNormFactors(miRNA_matrix, method = "TMM")
miRNA_matrix_tmm$samples


pseudo_TMM <- log2(cpm(miRNA_matrix_tmm) + 1)

#PCA
resPCA <- mixOmics::pca(t(pseudo_TMM), center = T, scale = T, ncomp = 10)
plot(resPCA)
colores <- brewer.pal(3, "Dark2")
plotIndiv(resPCA, group = pData_f$Version, col.per.group = colores[1:2], pch = c(1,2), legend = T, legend.title = "Group", ind.names = T, title = "PCA pseudo counts") 
# Add legend

PCA <- prcomp(t(pseudo_TMM), scale. = T, center = T)
pca.var <- PCA$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Group = pData_f$Tumor_Type_Mut, Version = pData_f$Version) 
#dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Group = pData_f$Tumor_Type_Mut, Version = paste(pData_f$Version,pData_f$Source_Year, sep ="_"))
pal = brewer.pal(8, "Dark2")

pca_TMM <- ggplot(dataGG, aes(PC1, PC2, label =rownames(dataGG))) +
  geom_point(aes(shape = Group, colour = Version), size=5) +
  ggtitle("PCA TMM") +
  xlab(paste0("PC1, VarExp: ", pca.var.per[1], "%")) +
  ylab(paste0("PC2, VarExp: ", pca.var.per[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  #coord_fixed(ratio = 1) +
  scale_shape_manual(values =  c(15:19,1,8)) + 
  #scale_color_manual(values = pal) +  
  theme_bw(base_size = 18)  + 
  theme(axis.text.x = element_text(size = 25, colour = "black"), axis.text.y = element_text(size = 25, colour = "black"),axis.title.y = element_text(size = 18, colour = "black"), legend.text = element_text(size = 20, colour = "black") ) + geom_text_repel(aes(label = rownames(dataGG)), size = 3) #+   theme(aspect.ratio=3/3)
pca_TMM

# Scatter Plot of Median and Variance per Sample

# Calculate median per sample
median_per_sample <- apply(pseudo_TMM, 2, median)

# Calculate variance per sample
variance_per_sample <- apply(pseudo_TMM, 2, var)

# Create a new dataframe with median and variance per sample
df_summary <- data.frame(median_per_sample, variance_per_sample)

# Add sample names as a column
df_summary$sample <- pData_f$Tumor_Type_Mut

# Convert 'sample' column to factor
df_summary$sample <- as.factor(df_summary$sample)

# Plot scatter plot using ggplot2
per_med_var_TMM <- ggplot(df_summary, aes(x = median_per_sample, y = variance_per_sample, color = sample)) +
  geom_point(size = 3) +
  xlab("Per-Sample Median") +
  ylab("Per-Sample Variance") +
  labs(title = "Scatter Plot of Median and Variance per Sample \n EdgeR TMM counts") +
  theme_minimal() + geom_text_repel(aes(label = rownames(df_summary)), size = 3) + 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_shape_manual(values = c(15:19,7:9)) + 
  scale_color_manual(values = pal) +  theme_bw(base_size = 18)  + 
  theme(axis.text.x = element_text(size = 25, colour = "black"), axis.text.y = element_text(size = 25, colour = "black"),axis.title.y = element_text(size = 18, colour = "black"), legend.text = element_text(size = 20, colour = "black") )


df_TMM = pseudo_TMM - apply(pseudo_TMM,1,median)
df_TMM <- melt(df_TMM, id = rownames(df_TMM))
names(df_TMM)[1:2] <- c ("id", "sample")
df_TMM$method <- rep("Edger_TMM", nrow(df_TMM))


## 2.4. RLE ----
miRNA_matrix_rle <- calcNormFactors(miRNA_matrix, method = "RLE")
miRNA_matrix_rle$samples

pseudo_RLE <- log2(cpm(miRNA_matrix_rle) + 1)

#PCA
resPCA <- mixOmics::pca(t(pseudo_RLE), center = T, scale = T, ncomp = 10)
plot(resPCA)
colores <- brewer.pal(3, "Dark2")
plotIndiv(resPCA, group = pData_f$Version, col.per.group = colores[1:2], pch = c(1,2), legend = T, legend.title = "Group", ind.names = T, title = "PCA pseudo counts") 

# Add legend
PCA <- prcomp(t(pseudo_RLE), scale. = T, center = T)
pca.var <- PCA$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)


dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Group = pData_f$Group2, Version = pData_f$Group2) 
#dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Group = pData_f$Tumor_Type_Mut, Version = paste(pData_f$Version,pData_f$Source_Year, sep ="_"))
pal = brewer.pal(8, "Dark2")

pca_RLE <- ggplot(dataGG, aes(PC1, PC2, label =rownames(dataGG))) +
  geom_point(aes(shape = Group, colour = Version), size=5) +
  ggtitle("PCA RLE") +
  xlab(paste0("PC1, VarExp: ", pca.var.per[1], "%")) +
  ylab(paste0("PC2, VarExp: ", pca.var.per[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  #coord_fixed(ratio = 1) +
  scale_shape_manual(values = c(15:19,7:9)) + 
  #scale_color_manual(values = pal) +  
  theme_bw(base_size = 18)  + 
  theme(axis.text.x = element_text(size = 25, colour = "black"), axis.text.y = element_text(size = 25, colour = "black"),axis.title.y = element_text(size = 18, colour = "black"), legend.text = element_text(size = 20, colour = "black") ) + geom_text_repel(aes(label = rownames(dataGG)), size = 3) #+   theme(aspect.ratio=3/3)
pca_RLE


# Scatter Plot of Median and Variance per Sample
# Calculate median per sample
median_per_sample <- apply(pseudo_RLE, 2, median)

# Calculate variance per sample
variance_per_sample <- apply(pseudo_RLE, 2, var)

# Create a new dataframe with median and variance per sample
df_summary <- data.frame(median_per_sample, variance_per_sample)

# Add sample names as a column
df_summary$sample <- pData_f$Group2

# Convert 'sample' column to factor
df_summary$sample <- as.factor(df_summary$sample)

# Plot scatter plot using ggplot2
per_med_var_RLE <- ggplot(df_summary, aes(x = median_per_sample, y = variance_per_sample, color = sample)) +
  geom_point(size = 3) +
  xlab("Per-Sample Median") +
  ylab("Per-Sample Variance") +
  labs(title = "Scatter Plot of Median and Variance per Sample \n EdgeR RLE counts") +
  theme_minimal() + geom_text_repel(aes(label = rownames(df_summary)), size = 3) + 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_shape_manual(values = c(15:19,7:9)) + 
  scale_color_manual(values = pal) +  theme_bw(base_size = 18)  + 
  theme(axis.text.x = element_text(size = 25, colour = "black"), axis.text.y = element_text(size = 25, colour = "black"),axis.title.y = element_text(size = 18, colour = "black"), legend.text = element_text(size = 20, colour = "black") )


df_RLE = pseudo_RLE - apply(pseudo_RLE,1,median)
df_RLE <- melt(df_RLE, id = rownames(df_RLE))
names(df_RLE)[1:2] <- c ("id", "sample")
df_RLE$method <- rep("Edger_RLE", nrow(df_RLE))



ggpubr::ggarrange(pca_raw  +  theme_bw(base_size = 12)  + theme(legend.position = "none") , pca_TC +  theme_bw(base_size = 12) + theme(legend.position = "none") , pca_UQ  +  theme_bw(base_size = 12) + theme(legend.position = "none"), pca_TMM +  theme_bw(base_size = 12) + theme(legend.position = "none"), pca_RLE +  theme_bw(base_size = 12) + theme(legend.position = "none"), ncol = 3, nrow = 2, common.legend = T, legend = "right")


df_allnorm <- rbind(df_raw, df_TC, df_UQ, df_TMM, df_RLE)
df_allnorm$method <- factor(df_allnorm$method,
                            levels = c("Raw counts", "Edger_TC", 
                                       "Edger_UQ", "Edger_TMM","Edger_RLE"))

p <- ggplot(data=df_allnorm, aes(x=sample, y=value, fill=method))
p <- p + geom_boxplot()  
p <- p + theme_bw()
p <- p + ggtitle("RLE plot")
p <- p + facet_grid(. ~ method) 
p <- p + ylab('Median deviation of log expression') + xlab("")
p1 <- p + theme(title = element_text(size=10), axis.text.x = element_blank(), 
               axis.ticks.x = element_blank())
print(p1)


p <- ggplot(data=df_allnorm, aes(x=value, colour=sample))
p <- p + geom_density()  
p <- p + theme_bw()
p <- p + ggtitle("Density of Median deviation of log expression")
p <- p + facet_grid(. ~ method) 
p <- p + ylab(expression(log[2] ~ (normalized ~ count + 1))) + xlab("")
p2 <- p + theme(title = element_text(size=10), legend.position = "none")
print(p2)

pdf("RLEplot.pdf", width = 15, height = 4)

p1
p2

dev.off()


# PCA Cowplot
pc <- prcomp(t(pseudo_RLE), scale. = T, center = T)
pca.var <- pc$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
df <- cbind(pc$x[,1:2], pData_f$Sample_Diagnosis, pData_f$Tumor_Type_Mut) %>% as.data.frame()
#df <- cbind(pc$x[,1:2], pData_f$Version, pData_f$Tumor_Type_PDTC) %>% as.data.frame()
df$PC1 <- as.numeric(df$PC1) / (pc$sdev[1] * sqrt(nrow(pData_f))) # scales the principal component values
df$PC2 <- as.numeric(df$PC2) / (pc$sdev[2] * sqrt(nrow(pData_f))) # scales the principal component values
df$V3 <- as.factor(df$V3)
df$V4 <- as.factor(df$V4)

# plot
#p3 <- ggplot(df, aes(PC1, PC2, colour = V4)) +
#  geom_point(size = 5, aes(shape = V3,colour = V4)) +
  #scale_shape_manual(values = c(15:19,7:9)) +
#  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))),
#               data = df[df$V3 == "1" | df$V3 == "2",], size = 1) +
#  xlab(paste0("PC1, VarExp: ", pca.var.per[1], "%")) +
#  ylab(paste0("PC2, VarExp: ", pca.var.per[2], "%")) +  ggtitle("PCA pseudo TMM") +
#  theme(plot.title = element_text(hjust = 0.5)) +  theme_bw(base_size = 12)  + 
#  theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"),axis.title.y = element_text(size = 12, colour = "black"), legend.text = element_text(size = 12, colour = "black") ) 
p3 <- ggplot(df, aes(PC1, PC2, label = rownames(df))) +
  geom_point(size = 5, aes(shape = V3, colour = V4)) +
  xlab(paste0("PC1, VarExp: ", pca.var.per[1], "%")) +
  ylab(paste0("PC2, VarExp: ", pca.var.per[2], "%")) +
  ggtitle("PCA combat-seq pseudo TMM") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black"),
    legend.text = element_text(size = 12, colour = "black")
  ) +
  geom_text_repel(aes(label = rownames(df)), size = 3) #+ scale_shape_manual(values =  c(15:18)) 

# Create the main plot without the legend
main_plot <- p3 +
  theme(legend.position = "none")  # Remove legend from the main plot

# Create a new plot containing only the legend
legend_plot <- get_legend(p3)

legend_plot <- cowplot::plot_grid(legend_plot)

# Print or display the final arranged plot
print(legend_plot)

library(cowplot)
# Add density curves to y and x axis
xdens <- 
  axis_canvas(main_plot, axis = "x") + 
  geom_density(data = df, aes(x = PC1, fill = V4, colour = V4), alpha = 0.3)
ydens <-
  axis_canvas(main_plot, axis = "y", coord_flip = TRUE) + 
  geom_density(data = df, aes(x = PC2, fill = V4, colour = V4), alpha = 0.3) +
  coord_flip()
p_rle <- main_plot %>%
  insert_xaxis_grob(xdens, grid::unit(1, "in"), position = "top") %>%
  insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
  ggdraw()

# Step 1: Call the pdf command to start the plot
pdf(file = "cowplot_pca_ggarrange.pdf",   # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 12) # The height of the plot in inches

ggpubr::ggarrange(p_raw, p_tc, p_uq, p_tmm, p_rle, legend_plot, ncol = 3, nrow = 2)

dev.off()


# 3. Batch correction ----

## 3.1. Variance fractions from each model fit ----

version <- ((pData_f$Version))
source <- pData_f$Source
groups <- ((pData_f$Tumor_Type_Mut))
exp_date <- pData_f$Experiment_Run_date
batch <- pData_f$Batch
lib_size <- miRNA_matrix_rle[["samples"]][["lib.size"]]
individual <- as.character(c(1:38))

# Specify variables to consider
# Age is continuous so we model it as a fixed effect
# Individual and Tissue are both categorical, so we model them as random effects
info <- data.frame(colnames(pseudo_TMM),  source, groups, version, exp_date, individual)
rownames(info) <- colnames(pseudo_TMM)

# Specify variables to consider
# continuous variables,  we model it as a fixed effect
# categorical variables, we model them as random effects
form <- ~ +(1|source) + (1|groups) + (1|version) + (1|exp_date) 

# Step 1: fit linear mixed model on gene expression
# If categorical variables are specified, a linear mixed model is used
# If all variables are modeled as continuous, a linear model is used
# each entry in results is a regression model fit on a single gene
# Step 2: extract variance fractions from each model fit
# for each gene, returns fraction of variation attributable to each variable 
# Interpretation: the variance explained by each variable
# after correction for all other variables
library(variancePartition)
varPart <- fitExtractVarPartModel(pseudo_TMM, form, info)

vp <- sortCols(varPart)
plotPercentBars(vp[1:20, ])

# violin plot of contribution of each variable to total variance
plotVarPart( sortCols( varPart ), main = "Fraction of variation attributable to each variable: TMM + combat" )


pdf("varplot_combat.pdf", onefile = T, width = 15, height = 8)
ggpubr::ggarrange(b, c,  ncol = 2, nrow = 2)
dev.off()


# Compute Canonical Correlation Analysis (CCA) between all pairs of variables returns absolute correlation value
form <- ~ (source) + (groups) + (version) + (exp_date) 
# returns absolute correlation value
C <- canCorPairs(form, info)
# Plot correlation matrix
# between all pairs of variables
plotCorrMatrix(C, col = c(brewer.pal("BuPu", n = 8)))
corrplot(C, method="color",  order="hclust",   col=brewer.pal(n=8, name="BuPu"),  tl.col="black")



## 3.2. Corrections (Combat-seq, Combat, RemoveBatchEffect) ----

#### 3.2.1. Combat-seq ----
pData_f$Experiment_Run_date <- ifelse(is.na(pData_f$Experiment_Run_date), "-", pData_f$Experiment_Run_date)
batches <- as.numeric(factor(paste(pData_f$Experiment_Run_date, pData_f$Source, sep = "_")))
batches <- as.numeric(factor(pData_f$Version))
groups <- as.numeric(factor(pData_f$Tumor_Type_Mut))
corrected_data = ComBat_seq(counts = as.matrix(raw_f), batch = batches, group = groups)

group <- factor(pData_f$Tumor_Type_Mut)
miRNA_matrix <- DGEList(counts=corrected_data, group = group)
miRNA_matrix_tmm.combat <- calcNormFactors(miRNA_matrix, method = "TMM")
miRNA_matrix_tmm.combat$samples

pseudo_TMM.combat <- log2(cpm(miRNA_matrix_tmm.combat) + 1)


pc <- prcomp(t(pseudo_TMM), scale. = T, center = T)
pca.var <- pc$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
df <- cbind(pc$x[,1:2], paste(pData_f$Version, pData_f$Experiment_Run_date,sep = "_"), pData_f$Tumor_Type_Mut) %>% as.data.frame()
df$PC1 <- as.numeric(df$PC1) / (pc$sdev[1] * sqrt(nrow(pData_f)))
df$PC2 <- as.numeric(df$PC2) / (pc$sdev[2] * sqrt(nrow(pData_f)))
df$V3 <- as.factor(df$V3)
df$V4 <- as.factor(df$V4)

# plot
p1 <- ggplot(df, aes(PC1, PC2, colour = V4)) +
  geom_point(size = 5, aes(shape = V3,colour = V4)) +
  scale_shape_manual(values = c(15,15,19,19,19)) +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))),
               data = df[df$V3 == "1" | df$V3 == "2",], size = 1) +
  xlab(paste0("PC1, VarExp: ", pca.var.per[1], "%")) +
  ylab(paste0("PC2, VarExp: ", pca.var.per[2], "%")) +  ggtitle("source" ) +
  xlab(paste0("PC1, VarExp: ", pca.var.per[1], "%")) +
  ylab(paste0("PC2, VarExp: ", pca.var.per[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +  theme_bw(base_size = 10)  + 
  theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"),axis.title.y = element_text(size = 12, colour = "black"), legend.text = element_text(size = 12, colour = "black") ) 

# Create the main plot without the legend
main_plot <- p1 +
  theme(legend.position = "none")  # Remove legend from the main plot

# Create a new plot containing only the legend
legend_plot <- get_legend(p1)

legend_plot <- cowplot::plot_grid(legend_plot)

# Print or display the final arranged plot
print(legend_plot)

# Add density curves to y and x axis
xdens <- 
  axis_canvas(main_plot, axis = "x") + 
  geom_density(data = df, aes(x = PC1, fill = V4, colour = V4), alpha = 0.3)
ydens <-
  axis_canvas(main_plot, axis = "y", coord_flip = TRUE) + 
  geom_density(data = df, aes(x = PC2, fill = V4, colour = V4), alpha = 0.3) +
  coord_flip()
raw_p <- main_plot %>%
  insert_xaxis_grob(xdens, grid::unit(1, "in"), position = "top") %>%
  insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
  ggdraw()

pdf(file = "cowplot_pca_combat_ggarrange.pdf",   # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 12) # The height of the plot in inches
ggpubr::ggarrange(raw_p, version.exp.source_p, source.version_p, version.exp_p, exp.source_p, legend_plot,  ncol = 3, nrow = 2)
dev.off()
ggpubr::ggarrange(raw_p, exp_p, version_p, source_p, legend_plot,  ncol = 3, nrow = 2)


varPart <- fitExtractVarPartModel(pseudo_TMM.combat, form, info)

# violin plot of contribution of each variable to total variance
source_v <- plotVarPart( sortCols( varPart ))

ggpubr::ggarrange(raw_v + theme_bw(base_size = 10) + theme(legend.position = "none") +  ggtitle("no combat-seq" ), version.exp.source_v+ theme_bw(base_size = 10) + theme(legend.position = "none") +  ggtitle("version + experiment date + source" ), source.version_v+ theme_bw(base_size = 10) + theme(legend.position = "none")+  ggtitle("source + version" ), version.exp_v+ theme_bw(base_size = 10) + theme(legend.position = "none")+  ggtitle("version + experiment date" ), exp.source_v+ theme_bw(base_size = 10) + theme(legend.position = "none") +  ggtitle("experiment date + source" ), exp_v+ theme_bw(base_size = 10) + theme(legend.position = "none") +  ggtitle("experiment run date" ), version_v+ theme_bw(base_size = 10) + theme(legend.position = "none") +  ggtitle("version" ), source_v+ theme_bw(base_size = 10) + theme(legend.position = "none") +  ggtitle("source" ),  ncol = 4, nrow = 2)


### 3.2.2. Combat ----

batches <- as.numeric(factor(paste(pData_f$Source_Year, pData_f$Experiment_Run_date, sep = "_")))
batches <- as.numeric(factor(pData_f$Source_Year))
groups <- as.numeric(factor(pData_f$Tumor_Type_Mut))

pseudo_TMM.combat <- (ComBat((pseudo_TMM), batch = batches, 
                          mod = model.matrix(~0+groups), par.prior = F, prior.plots = F))

pc <- prcomp(t(pseudo_TMM.combat), scale. = T, center = T)
pca.var <- pc$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
df <- cbind(pc$x[,1:2], paste(pData_f$Version, pData_f$Experiment_Run_date,sep = "_"), pData_f$Tumor_Type_Mut) %>% as.data.frame()
df$PC1 <- as.numeric(df$PC1) / (pc$sdev[1] * sqrt(nrow(pData_f)))
df$PC2 <- as.numeric(df$PC2) / (pc$sdev[2] * sqrt(nrow(pData_f)))
df$V3 <- as.factor(df$V3)
df$V4 <- as.factor(df$V4)

# plot
p1 <- ggplot(df, aes(PC1, PC2, colour = V4)) +
  geom_point(size = 5, aes(shape = V3,colour = V4)) +
  scale_shape_manual(values = c(15,15,19,19,19)) +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))),
               data = df[df$V3 == "1" | df$V3 == "2",], size = 1) +
  xlab(paste0("PC1, VarExp: ", pca.var.per[1], "%")) +
  ylab(paste0("PC2, VarExp: ", pca.var.per[2], "%")) +  ggtitle("source" ) +
  xlab(paste0("PC1, VarExp: ", pca.var.per[1], "%")) +
  ylab(paste0("PC2, VarExp: ", pca.var.per[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +  theme_bw(base_size = 10)  + 
  theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"),axis.title.y = element_text(size = 12, colour = "black"), legend.text = element_text(size = 12, colour = "black") ) 

# Create the main plot without the legend
main_plot <- p1 +
  theme(legend.position = "none")  # Remove legend from the main plot

# Create a new plot containing only the legend
legend_plot <- get_legend(p1)

legend_plot <- cowplot::plot_grid(legend_plot)

# Print or display the final arranged plot
print(legend_plot)

# Add density curves to y and x axis
xdens <- 
  axis_canvas(main_plot, axis = "x") + 
  geom_density(data = df, aes(x = PC1, fill = V4, colour = V4), alpha = 0.3)
ydens <-
  axis_canvas(main_plot, axis = "y", coord_flip = TRUE) + 
  geom_density(data = df, aes(x = PC2, fill = V4, colour = V4), alpha = 0.3) +
  coord_flip()
source_p<- main_plot %>%
  insert_xaxis_grob(xdens, grid::unit(1, "in"), position = "top") %>%
  insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
  ggdraw()

ggpubr::ggarrange(raw_p, version.exp.source_p, source.version_p, version.exp_p, exp.source_p, legend_plot,  ncol = 3, nrow = 2)

ggpubr::ggarrange(raw_p, exp_p, version_p, source_p, legend_plot,  ncol = 3, nrow = 2)


varPart <- fitExtractVarPartModel(pseudo_TMM.combat, form, info)

### 3.2.3. Limma - removeBatchEffect ----

pseudo_TMM.limma <- removeBatchEffect(pseudo_TMM, batch = source,batch2 = paste(pData_f$Version, pData_f$Experiment_Run_date, sep = "_"), design = design)

pc <- prcomp(t(pseudo_TMM.combat), scale. = T, center = T)
pca.var <- pc$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
df <- cbind(pc$x[,1:2], paste(pData_f$Version, pData_f$Experiment_Run_date,sep = "_"), pData_f$Tumor_Type_Mut) %>% as.data.frame()
df$PC1 <- as.numeric(df$PC1) / (pc$sdev[1] * sqrt(nrow(pData_f)))
df$PC2 <- as.numeric(df$PC2) / (pc$sdev[2] * sqrt(nrow(pData_f)))
df$V3 <- as.factor(df$V3)
df$V4 <- as.factor(df$V4)

# plot
p1 <- ggplot(df, aes(PC1, PC2, colour = V4)) +
  geom_point(size = 5, aes(shape = V3,colour = V4)) +
  scale_shape_manual(values = c(15,15,19,19,19)) +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))),
               data = df[df$V3 == "1" | df$V3 == "2",], size = 1) +
  xlab(paste0("PC1, VarExp: ", pca.var.per[1], "%")) +
  ylab(paste0("PC2, VarExp: ", pca.var.per[2], "%")) +  ggtitle("source" ) +
  xlab(paste0("PC1, VarExp: ", pca.var.per[1], "%")) +
  ylab(paste0("PC2, VarExp: ", pca.var.per[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +  theme_bw(base_size = 10)  + 
  theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"),axis.title.y = element_text(size = 12, colour = "black"), legend.text = element_text(size = 12, colour = "black") ) 

# Create the main plot without the legend
main_plot <- p1 +
  theme(legend.position = "none")  # Remove legend from the main plot

# Create a new plot containing only the legend
legend_plot <- get_legend(p1)

legend_plot <- cowplot::plot_grid(legend_plot)

# Print or display the final arranged plot
print(legend_plot)

# Add density curves to y and x axis
xdens <- 
  axis_canvas(main_plot, axis = "x") + 
  geom_density(data = df, aes(x = PC1, fill = V4, colour = V4), alpha = 0.3)
ydens <-
  axis_canvas(main_plot, axis = "y", coord_flip = TRUE) + 
  geom_density(data = df, aes(x = PC2, fill = V4, colour = V4), alpha = 0.3) +
  coord_flip()
source_p <- main_plot %>%
  insert_xaxis_grob(xdens, grid::unit(1, "in"), position = "top") %>%
  insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
  ggdraw()

ggpubr::ggarrange(raw_p, version.exp.source_p, source.version_p, version.exp_p, exp.source_p, legend_plot,  ncol = 3, nrow = 2)

ggpubr::ggarrange(raw_p, exp_p, version_p, source_p, legend_plot,  ncol = 3, nrow = 2)


varPart <- fitExtractVarPartModel(pseudo_TMM.combat, form, info)

# violin plot of contribution of each variable to total variance
source_v <- plotVarPart( sortCols( varPart ))

ggpubr::ggarrange(raw_v + theme_bw(base_size = 10) + theme(legend.position = "none") +  ggtitle("no combat-seq" ), version.exp.source_v+ theme_bw(base_size = 10) + theme(legend.position = "none") +  ggtitle("version + experiment date + source" ), source.version_v+ theme_bw(base_size = 10) + theme(legend.position = "none")+  ggtitle("source + version" ), version.exp_v+ theme_bw(base_size = 10) + theme(legend.position = "none")+  ggtitle("version + experiment date" ), exp.source_v+ theme_bw(base_size = 10) + theme(legend.position = "none") +  ggtitle("experiment date + source" ), exp_v+ theme_bw(base_size = 10) + theme(legend.position = "none") +  ggtitle("experiment run date" ), version_v+ theme_bw(base_size = 10) + theme(legend.position = "none") +  ggtitle("version" ), source_v+ theme_bw(base_size = 10) + theme(legend.position = "none") +  ggtitle("source" ),  ncol = 4, nrow = 2)
