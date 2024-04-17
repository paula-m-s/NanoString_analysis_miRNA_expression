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
###     https://github.com/bhattacharya-a-bt/CBCS_normalization/
###   Date started: April 2023
###   Last date modified: April 2024
###   R version 
###   Script: Raw data and QC from NanoString miRNA RCC files. Versions: NS_H_miR_v3a and NS_H_miR_v3b
###
###   Platform: NanoString. RStudio
########################################################################################################################.
########################################################################################################################.

# 0. PIPELINE PREPARATION ----
# Install packages
#install.packages("remotes")
#remotes::install_github("calebclass/NanoTube")

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#devtools::install_github("Nanostring-Biostats/NanoStringNCTools", force = TRUE, ref = "master")

#install.packages("xlsx")
#BiocManager::install("sva")
# (!) Check the libraries that you would need in advance and install them

# Load packages

library(NanoStringNCTools) 
library(NanoStringQCPro)
library(RColorBrewer)
library(NanoTube)
library(NanoStringDiff)
library(pheatmap)
library(xlsx)
library(ComplexHeatmap)
library(viridis)
library(RUVSeq)
library(ggrepel)
library(plotly)
library(FactoMineR)
library(factoextra)
library(edgeR)
library(ruv)
library(rgl)
library(sva)
library(EnhancedVolcano)
library(EnvStats)
library(DESeq2)
library(MASS)
library(reshape2)

# Set working directory
dir<- "C:/path/to/your/files"
setwd(dir)


# 1. Import data, standarization, raw data extraction and quality control. ----

## 1.1. Load corrected phenodata for all samples ----
library(readr)
metadata <-  read_csv("C:/path/to/your/metadata.csv")
#View(metadata)

## 1.2. Load rcc files ----
nano_data <-  c(file.path(dir,"/folder/batch/two"), file.path(dir,"/folder/batch/two"))
sample_info <- file.path(dir,"/path/to/your/metadata.csv")

files.RCC <- list.files(nano_data, pattern = 'RCC', full.names = TRUE) # read in RCC files

## 1.3. Raw data extraction (raw counts, pData and fData) ----

# IMPORTANT TO NOTE: Need for probe homogenization. The two sets of samples are made with different versions of the NanoString for miRNA. The only differences between the V3a and V3b constructs are the ligation-dependent calibration factors: all probes and barcodes (miRNA accessions) are identical. As this is the case, with proposals to analyze the quality of all samples in full, the calibration factors that appear in the RCC files will be removed. The code below serves that purpose. Subsequently, both sets of samples are grouped for analysis.
# Raw data: raw counts (Endogenous, SpikeIn, Housekeeping, Ligation, Positive and Negative control probes), pData (phenodata table that contains all the information refered to the samples and the QC information), fData (Probe names, Accessions, Code Class)

### 1.3.1. Raw expression dataframe----
## Note: raw counts (Endogenous, SpikeIn, Housekeeping, Ligation, Positive and Negative control probes)

# First we make sure that the number of rows in 'Code_Summary' (other can be chosen) for each RCC file and store them in a vector
num_rows <- sapply(files.RCC, function(file) {
  nrow(readRcc(file)$Code_Summary)
})

# Check if all values in the vector are the same
all_same_ng <- all(num_rows == num_rows[1])
if (all_same_ng) {
  cat("All RCC files have the same number of rows in 'Code_Summary'.\n")
} else {
  cat("RCC files have different numbers of rows in 'Code_Summary'.\n")
}


#Retrieve the number of rows
ng <- nrow(readRcc(files.RCC[1])$Code_Summary)

#Retrieve the number of columns
ncol <- length(files.RCC)

#Empty raw counts data frame
raw_expression <- as.data.frame(matrix(nrow = ng, ncol = ncol + 2))
colnames(raw_expression)[1:2] <- c('Gene', 'Class')
raw_expression[, 1:2] <- readRcc(files.RCC[1])$Code_Summary[, c(2, 1)]
raw_expression$Class <- gsub("Endogenous1", "Endogenous", raw_expression$Class)
raw_expression$Gene <- gsub("\\|.*", "", raw_expression$Gene) #remove ligation-dependent calibration factors 

for (i in 1:length(files.RCC)){
  
  print(i)
  rcc <- readRcc(files.RCC[i])
  
  raw_expression[, i+2] <- as.numeric(rcc[["Code_Summary"]]$Count)
  colnames(raw_expression)[i+2] <- strsplit(files.RCC[i], '/')[[1]][10] #(!) Check the length of the string
}

### 1.3.2. pData dataframe ----
# Note: phenodata table that contains all the information refered to the samples and the QC information

#Empty pData data frame
pData <- as.data.frame(matrix(nrow = ncol, ncol = 15))
colnames(pData) <- c( 'SampleID', 'Owner', 'Comments', 
                      'Date', 'GeneRLF', 'SystemAPF', 'imagingQC_flag',
                      'imagingQC_num','bindingDensityQC_flag',
                      'bindingDensityQC_num', 'limitOfDetectionQC',
                      "posE", "posF", 'limitOfDetection',
                      'positiveLinearityQC')

# For each rcc file add the counts and QC information to the dataframes previously generated
for (i in 1:length(files.RCC)){
  
  print(i)
  rcc <- readRcc(files.RCC[i])
  
  pData[i, 1:6] <- as.vector(rcc$Sample_Attributes)
  pData[i, "SampleID"] <- strsplit(files.RCC[i], '/')[[1]][10] #(!) Check the length of the string
  pData$imagingQC_num[i] <- imagingQC(rcc)$fovRatio
  pData$imagingQC_flag[i] <- imagingQC(rcc)$flag
  pData$bindingDensityQC_flag[i] <- bindingDensityQC(rcc, .1, 2.25)$flag
  pData$bindingDensityQC_num[i] <- bindingDensityQC(rcc, .1, 2.25)$bd
  pData$limitOfDetectionQC[i] <- limitOfDetectionQC(rcc)$flag
  pData$limitOfDetection[i] <- limitOfDetectionQC(rcc)$lod
  pData$posE[i] <- limitOfDetectionQC(rcc)$pose
  pData$posF[i] <- limitOfDetectionQC(rcc)$posf
  pData$positiveLinearityQC[i] <- positiveLinQC(rcc)
}

# Add more QC information to the pData dataframe: positive scale factor, 'positive r squared', 'Positive probes QC flag', "Housekeeping scale factor", "Housekeeping flag" 

qc_analysis_2018 <- pos_HK_QC(nano_data[1])
qc_analysis_2022 <- pos_HK_QC(nano_data[2])
qc_analysis <- rbind(qc_analysis_2018, qc_analysis_2022)
pData <- merge(pData, qc_analysis, by.x = "SampleID", by.y = "SampleID")
rm(list = c("qc_analysis", "qc_analysis_2022", "qc_analysis_2018"))


# Add more QC information to the pData dataframe: Overal assay efficiency
pData$meanPos <- log2(colMeans(raw_expression[(raw_expression$Class == "Positive"), 3:46])) # means positive 
meanPosMean <- (mean(pData$meanPos)) # mean positive all samples
pData$meanPosRatio <- pData$meanPos/meanPosMean 

if (mean(pData$meanPosRatio) > 3 | mean(pData$meanPosRatio) < 1/3) {
  pData$Overal.assay.efficiency <- rep("Flag", length(pData$meanPosRatio))
} else {
  pData$Overal.assay.efficiency <- rep("No flag", length(pData$meanPosRatio))
}

table(pData$Overal.assay.efficiency)


#Review the results from the positive and negative controls. Positive controls with low counts or negative controls with counts significantly above background can trigger flags and should be checked to see if they indicate more serious issues with the data. 
for (i in 1:length(raw_expression[,3:46])){ #keep in mind the structure of your data
  print(i)
  pData$num.below.bg.end[i] <- colSums((raw_expression[raw_expression$Class == "Endogenous",3:46][i]) < pData$limitOfDetection[i])
  pData$num.below.bg.hk[i] <- colSums((raw_expression[raw_expression$Class == "Housekeeping",3:46][i]) < pData$limitOfDetection[i])
  pData$num.above.bg.neg[i] <- colSums((raw_expression[raw_expression$Class == "Negative",3:46][i]) > pData$limitOfDetection[i])
}

# Number of probes under 'n' counts
pData$end.blw.10.counts = apply((raw_expression[raw_expression$Class == "Endogenous",3:46]), 2, function(x) sum(x <= 10))
pData$hk.blw.20.counts = apply((raw_expression[raw_expression$Class == "Housekeeping",3:46]), 2, function(x) sum(x <= 20))
pData$neg.blw.10.counts = apply((raw_expression[raw_expression$Class == "Negative",3:46]), 2, function(x) sum(x <= 10))
pData$pos.blw.10.counts = apply((raw_expression[raw_expression$Class == "Positive",3:46]), 2, function(x) sum(x <= 10))
pData$end.blw.1.counts = apply((raw_expression[raw_expression$Class == "Endogenous",3:46]), 2, function(x) sum(x <= 1))
pData$hk.blw.1.counts = apply((raw_expression[raw_expression$Class == "Housekeeping",3:46]), 2, function(x) sum(x <= 1))
pData$neg.blw.1.counts = apply((raw_expression[raw_expression$Class == "Negative",3:46]), 2, function(x) sum(x <= 1))
pData$pos.blw.1.counts = apply((raw_expression[raw_expression$Class == "Positive",3:46]), 2, function(x) sum(x <= 1))
pData$spike.blw.10.counts = apply((raw_expression[raw_expression$Class == "SpikeIn",3:46]), 2, function(x) sum(x <= 10))

# Add count flags to the pData
pData$count_flag <- ifelse(as.numeric(pData$end.blw.10.counts) > 500, "Flag", "No flag")
table(pData$count_flag)


# Merge the sample information with the pData which contains the QC sample information
pData <- merge(pData, metadata, by.x = "SampleID", by.y = "RCC")


### 1.3.3. fData dataframe ---- 
#Dataframe and name column cleaning
fData <- readRcc(files.RCC[1])$Code_Summary[, c(2, 1, 3)]
colnames(fData)[1:3] <- c('Gene', 'Class', 'Accession')
fData$Class <- gsub("Endogenous1", "Endogenous", fData$Class)
fData$Gene <- gsub("\\|.*", "", fData$Gene)


### 1.3.4. Write results ----
# create folder if it doesn't exist
sub_dir_exists <- "/name/of/your/final/resutls/path/"
if (!file.exists(sub_dir_exists)){
  dir.create(file.path(dir, sub_dir_exists))
  message("\n Creating folder \n")
}else{message("\n directory already exists \n")
}
#write.csv(raw_expression  ,  "C:/name/of/your/final/resutls/path/raw_expression.csv", row.names = T)
#write.csv(pData,  "C:/name/of/your/final/resutls/path/pData.csv", row.names = T)
#write.csv(fData,  "C:/name/of/your/final/resutls/path/fData.csv.csv", row.names = T)

## 1.4. QC Plots ----

### 1.4.1. Counts per sample----
pData <- pData[order(pData$Tumor_Type_Mut),]
rownames(raw_expression) <- raw_expression$Gene
raw_expression2 <- raw_expression[,pData$SampleID]


end <-  raw_expression2[raw_expression$Class == "Endogenous",] #Change as needed for Endogenous, Housekeeping, SpikeIn...
colnames(end) <- pData$Sample_Names
data <- as.data.frame(colSums(end))
data$patients <- rownames(data)
names(data)[1] <- "Counts"
#data <- data[!grepl(pattern = "Removed", x = data$patients),]
data <- data[order(data$Counts),]
data$patients <- factor(data$patients, levels = data$patients[order(data$Counts)])
g_end <- ggplot(data, aes(x = Counts, y = patients)) +
  geom_point(shape = 21, fill = "#F8766D", size = 3) +
  theme_bw() +
  labs(x = "Counts Endogenous", y = "")

spike <- raw_expression2[raw_expression$Class == "SpikeIn",] #Change as needed for Endogenous, Housekeeping, SpikeIn...
colnames(spike) <- pData$Sample_Names
data <- as.data.frame(colSums(spike))
data$patients <- rownames(data)
names(data)[1] <- "Counts"
#data <- data[!grepl(pattern = "Removed", x = data$patients),]
data <- data[order(data$Counts),]
data$patients <- factor(data$patients, levels = data$patients[order(data$Counts)])
g_spike <- ggplot(data, aes(x = Counts, y = patients)) +
  geom_point(shape = 21, fill = "#E76BF3", size = 3) +
  theme_bw() +
  labs(x = "Counts SpikeIn", y = "")

hk <- raw_expression2[raw_expression$Class == "Housekeeping",] #Change as needed for Endogenous, Housekeeping, SpikeIn...
colnames(hk) <- pData$Sample_Names
data <- as.data.frame(colSums(hk))
data$patients <- rownames(data)
names(data)[1] <- "Counts"
#data <- data[!grepl(pattern = "Removed", x = data$patients),]
data <- data[order(data$Counts),]
data$patients <- factor(data$patients, levels = data$patients[order(data$Counts)])
g_hk <- ggplot(data, aes(x = Counts, y = patients)) +
  geom_point(shape = 21, fill = "#A3A500", size = 3) +
  theme_bw() +
  labs(x = "Counts Housekeeping", y = "")

posi <- raw_expression2[raw_expression$Class == "Positive",] #Change as needed for Endogenous, Housekeeping, SpikeIn...
colnames(posi) <- pData$Sample_Names
data <- as.data.frame(colSums(posi))
data$patients <- rownames(data)
names(data)[1] <- "Counts"
#data <- data[!grepl(pattern = "Removed", x = data$patients),]
data <- data[order(data$Counts),]
data$patients <- factor(data$patients, levels = data$patients[order(data$Counts)])
g_posi <- ggplot(data, aes(x = Counts, y = patients)) +
  geom_point(shape = 21, fill = "#00BF7D", size = 3) +
  theme_bw() +
  labs( x = "Counts Positive", y = "")

neg <- raw_expression2[raw_expression$Class == "Negative",] #Change as needed for Endogenous, Housekeeping, SpikeIn...
colnames(neg) <- pData$Sample_Names
data <- as.data.frame(colSums(neg))
data$patients <- rownames(data)
names(data)[1] <- "Counts"
#data <- data[!grepl(pattern = "Removed", x = data$patients),]
data <- data[order(data$Counts),]
data$patients <- factor(data$patients, levels = data$patients[order(data$Counts)])
g_neg <- ggplot(data, aes(x = Counts, y = patients)) +
  geom_point(shape = 21, fill = "#00B0F6", size = 3) +
  theme_bw() +
  labs(x = "Counts Negative", y = "")


ggarrange(g_end + theme(axis.text.y = element_text(size=7.5)), g_hk+ theme(axis.text.y = element_text(size=7.5)), g_posi+ theme(axis.text.y = element_text(size=7.5)), g_neg+ theme(axis.text.y = element_text(size=7.5)), g_spike+ theme(axis.text.y = element_text(size=7.5)), nrow = 1, ncol = 5)

data <- as.data.frame(colSums(end))
data$patients <- rownames(data)
names(data)[1] <- "Counts"
data <- data[!grepl(pattern = "Removed", x = data$patients),]
data <- data[order(data$Endogenous),]

# List of data frames to merge
dfs_list <- list(data, spike, posi, neg, hk)

# Merge data frames in the list based on the "ID" column
merged_df <- Reduce(function(x, y) merge(x, y, by = "patients"), dfs_list)
merged_df <- merged_df[order(merged_df$Endogenous),]

# Example: Assuming your data frame is named 'tidy_df'
tidy_df <- merged_df %>%
  gather(key = "variable", value = "value", -patients) %>%
  arrange(variable, desc(value))

# Your ggplot code
ggplot(tidy_df, aes(x = value, y = patients)) +
  geom_point(shape = 21, fill = "purple", size = 3) +
  theme_bw() +
  scale_x_reordered() +
  labs(title = "Counts by Patients", x = "Counts", y = "Patients") + 
  facet_grid(.~variable,  scales = "free_x")

### 2.2.2.  Overal.assay.efficiency----

a <- ggplot(pData, aes(x = Sample_Names, y = (meanPosRatio), color = Overal.assay.efficiency)) +
  geom_point( show.legend = FALSE) +
  theme_bw() +
  theme(#axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.text.x=element_text(size=8), 
    axis.text.x = element_text(size=15)) +
  ylab("mean pos ctrl counts/overall mean of pos ctrl counts") + xlab("") +
  ggtitle("Overal assay efficiency \n") +
  geom_hline(yintercept=c(1/3,3), color="red")+
  scale_color_manual(values=c("black", "red"), labels=c("No Flag", "Flag"), name = "Flag or no")  + coord_flip()
print(a)

### 2.2.3. Binding density QC ----

b <- ggplot(pData, aes(x = Sample_Names, y = bindingDensityQC_num, color = bindingDensityQC_flag)) + 
  geom_point( show.legend = FALSE) +
  theme_bw() +
  theme(#axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.text.x=element_text(size=8),
    axis.text.x = element_text(size=15)) +
  scale_y_continuous(expand = c(0,0)) +
  expand_limits(y = c(0,3)) +
  ylab("Binding density") +
  xlab("") +
  ggtitle("Binding Density QC Plot\nFlagged if < 0.05 or > 2.25") +
  geom_hline(yintercept=0.05, color="red") +
  geom_hline(yintercept=2.25, color="red") +
  scale_color_manual(values=c("black", "red"), labels=c("No Flag", "Flag"), name = "Flag or no") + coord_flip()
print(b)

### 2.2.4. Imaging QC ----

c<- ggplot(pData, aes(x = Sample_Names, y = (imagingQC_num*100), color = imagingQC_flag)) +
  geom_point( show.legend = FALSE) +
  theme_bw() +
  theme(#axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.text.x=element_text(size=8),
    axis.text.x = element_text(size=15)) +
  scale_y_continuous(expand = c(0,0)) +
  expand_limits(y = c(0,1.05 * max(pData$imagingQC_num*100))) +
  ylab("percentage of FOV counted") + xlab("") +
  ggtitle("FOV \nProblematic if < 75%") +
  geom_hline(yintercept=75, color="red")+
  scale_color_manual(values=c("black", "red"), labels=c("No Flag", "Flag"), name = "Flag or no") + coord_flip()
print(c)

library(ggpubr)
ggarrange(a+theme(axis.text.y=element_blank()),b+theme(axis.text.y=element_blank())+coord_flip(),c+theme(axis.text.y=element_blank())+coord_flip(), ncol = 3, nrow = 1)
ggarrange(a,b,c, ncol = 3, nrow = 1)
