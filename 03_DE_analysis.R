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
###   Script: Differential expression analysis_NanoString miRNA
###
###   Public data source:  
###   Platform: NanoString
########################################################################################################################.
########################################################################################################################.


# 1. Put the data into a DGEList object ----
group <- as.factor(pData_f$Tumor_Type_Mut)
miRNA_matrix <- DGEList(counts=raw_f, group = group)

miRNA_matrix_tmm <- calcNormFactors(miRNA_matrix, method = "TMM")
miRNA_matrix_tmm$samples


pseudo_TMM <- log2(cpm(miRNA_matrix_tmm) + 1)

# 2. Design matrix and address batch effects ----
batches <- as.factor(pData_f$Batch)
group <- as.factor(pData_f$Tumor_Type_Mut)
version <- gsub("NS_H_miR_", "", as.factor(pData_f$Version))
source <- as.factor(pData_f$Source_Year)
exp_date <- pData_f$Experiment_Run_date <- ifelse(is.na(pData_f$Experiment_Run_date), "-", pData_f$Experiment_Run_date)
exp_date <- gsub("^_", "_", exp_date)

design <- model.matrix(~0+group+batches+source, data=miRNA_matrix_tmm$samples) #add here to address batch effects
colnames(design) <- gsub("^group","",colnames(design))
design


# 3. Estimate data dispersion GLM approach ----
y <- estimateDisp(miRNA_matrix_tmm, design, robust = T)


# 4. Contrast matrix ----

# Create contrasts for each combination
mycontrast <- c("MNG_DGCR8_mutated-MNG_wildtype", 
                "MNG_DGCR8_mutated-TC_wildtype",
                "MNG_DICER1_mutated-MNG_wildtype", 
                "MNG_DICER1_mutated-TC_wildtype",
                "TC_DGCR8_mutated-TC_wildtype", 
                "TC_DICER1_mutated-TC_wildtype",
                "TC_DGCR8_mutated-MNG_wildtype", 
                "TC_DICER1_mutated-MNG_wildtype",
                "TC_DICER1_mutated-MNG_DICER1_mutated",
                "TC_DGCR8_mutated-MNG_DGCR8_mutated")

contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list
                                  (mycontrast),levels=list(design))))


# 5. Differential expression (DE) ----

# Fit the quasi-likelihood negative binomial generalized log-linear model 
fit <- glmQLFit(y, design)

# Create an empty list to store the DE gene results
de_genes_list <- list()


# Perform glmQLFTest for each contrast and save the results in the list
for (i in 1:ncol(contrast.matrix)) {
  contrast <- contrast.matrix[, i]
  qlfTest <- glmQLFTest(fit, contrast = contrast)
  
  # Save the results in the list with the contrast name
  contrast_name <- colnames(contrast.matrix)[i]
  de_genes <- topTags(qlfTest, n = Inf)$table
  de_genes <- de_genes[order(de_genes$FDR),]
  de_genes$names <- rownames(de_genes)
  de_genes_list[[contrast_name]] <- de_genes
}


# (!) NOTE (https://support.bioconductor.org/p/79149/): estimateGLMRobustDisp is distinct from robust=TRUE in estimateDisp. The former uses observation weights to reduce the impact of outlier counts within each gene. The latter doesn't protect each gene from outlier observations, but instead protects the empirical Bayes shrinkage from genes with outlier dispersions (that are caused by outlier counts). This difference in behaviour has some practical consequences. For example, a gene that is DE but has a couple of outlier counts in one group will (hopefully) still be detected as DE with estimateGLMRobustDisp. This is because those outlier counts should be downweighted, thus avoiding increases to the dispersion estimate. In contrast, with estimateDisp, the dispersion gets inflated by the outliers so there won't be enough power to detect DE - however, any deleterious effect of the inflated dispersion on EB shrinkage of all other genes is prevented with robust=TRUE. Using estimateDisp as it's faster (no need to iteratively compute weights), but there can be some benefit from robustifying against outlier observations - read the NAR:https://academic.oup.com/nar/article/42/11/e91/1427925 
# glmQLFit estimates QL dispersions, it doesn't re-estimate the NB dispersions; keep in mind that they are two separate sets of values. Technically, you don't need to set robust=TRUE in estimateDisp because that only affects the tagwise NB dispersions that are never used in the QL framework - only the trended NB dispersions are used. That said, robustifying doesn't do any harm so you might as well do it if you plan to use the tagwise NB dispersions for diagnostics later, e.g., with plotBCV. Note that the EB shrinkage is repeated within glmQLFit using the QL dispersions. Here, robust=TRUE in glmQLFit is analogous to that in estimateDisp, as it protects the EB shrinkage from outlier dispersions (QL, not NB). Finally, abundance.trend=TRUE is already the default, and just ensures that EB shrinkage of the QL dispersions is performed towards a mean-dependent trend based on the QL dispersions. This trend is distinct from that generated by estimateDisp, which concerns the NB dispersions.



# Access the results for each contrast as independent data frames
for (i in 1:length(de_genes_list)) {
  contrast_name <- names(de_genes_list)[i]
  de_genes <- as.data.frame(de_genes_list[[i]])
  assign(paste0("de_genes_", contrast_name), de_genes, envir = .GlobalEnv)
}


## 5.1 Save DE results ----
path_dir <- "Resultados/DE/All_Comp_Exp-Source_robust/"
if (!file.exists(path_dir)){
  dir.create(file.path(dir, path_dir))
  message("\n Creating folder \n")
}else{message("\n directory already exists \n")
}


# Get the names of all objects in your work environment that start with "DE"
object_names <- c(ls(pattern = "^de_genes_"))

# Save each object as CSV file to the specified path
for (nam in object_names) {
  arch <- paste(nam, ".csv", sep = "")
  path_archive <- file.path(path_dir, arch)
  write.csv(eval(as.symbol(nam)), path_archive, row.names = T)
}


# 6. DE exploration ----

## 6.1. Check number of significant miRNAs ----

# Create an empty dataframe to store the results
signif_res <- data.frame(Method = character(), P = numeric())

# Iterate over each element in the list. Distribution of all p-values
for (i in seq_along(de_genes_list)) {
  # element
  element <- de_genes_list[[i]]
  # name
  nam <- names(de_genes_list)[i]
  #Create the data frame using the current element and the corresponding name
  df <- data.frame(
    Method = rep(nam, each = nrow(element)),
    P = c(element$PValue)
  )
  # Join the current dataframe to the final result
  signif_res <- rbind(signif_res, df)
}


# Create an empty dataframe to store the results
signif_res <- data.frame(Method = character(), P = numeric())

# Iterate over each element in the list. Distribution of FDR in abs(log2FC) above 2
for (i in seq_along(de_genes_list)) {
  # element
  element <- de_genes_list[[i]]
  # names
  nam <- names(de_genes_list)[i]
  # Filter the data that meets the conditions: FDR < 0.05 and absolute logFC >= 2
  dat_filt <- element[element$FDR < 0.05 & abs(element$logFC) >= 2, ]
  # Create the data frame using the filtered data and the corresponding name
  df <- data.frame(
    Method = rep(nam, each = nrow(dat_filt)),
    P = c(dat_filt$FDR)
  )
  
  # Join the current dataframe to the final result
  signif_res <- rbind(signif_res, df)
}



# Filtrar el dataframe
#signif_res <- subset(signif_res, (signif_res$Method %in% methods_filt)) #subset list
signif_res$Method <- gsub("_", " ", signif_res$Method)
signif_res$Method <- gsub("mutated", "", signif_res$Method)
signif_res$Method <- gsub("-", "vs ", signif_res$Method)
signif_res$Method <- gsub("  ", " ", signif_res$Method)

pdf("pvalue_plots2.pdf", height = 8, width = 10)
ggplot(data = signif_res,
       aes(x = P, color = Method, fill = Method)) +
  geom_histogram(alpha = 0.2) +
  facet_wrap(~ Method, ncol = 2) +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 12),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(color = "grey", fill = NA, size = 0.1),
        legend.position = 'none',
        axis.ticks = element_line(size = 0.2)) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0.01, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0.001, linetype = "dashed", color = "seagreen") +
  ylab('Distribution of miRNAs') +
  xlab("p-value") +
  ggtitle("") +
  scale_x_continuous(trans = "log10", 
                     breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)))

dev.off()


## 6.2. Heatmap highly variable genes (option 1) ----

# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(pseudo_TMM, 1, var)
head(var_genes)

# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:50]
head(select_var)

# Subset logcounts matrix
highly_variable_lcpm <- pseudo_TMM[select_var,]
dim(highly_variable_lcpm)

# Annotation matrix
annotation_Group <- data.frame(Diagnosis = pData_f$Sample_Diagnosis, Source = pData_f$Source_Year, Group_Mutation_type = pData_f$Tumor_Type_Mut, Gender = pData_f$Sex, Development = pData_f$Development_stage, MAPK_Mutated = pData_f$MAPK_mutated)
rownames(annotation_Group) <- (pData_f$Sample_Names)
colnames(highly_variable_lcpm) <- (pData_f$Sample_Names)
colnames(annotation_Group) <- gsub("_", " ", colnames(annotation_Group))
annotation_Group$`Group Mutation type` <- gsub("_", " ", annotation_Group$`Group Mutation type`)

# Initialize an empty list to store the colors
annotation_colors <- list()
for (col in colnames(annotation_Group)) {
  levels <- unique(annotation_Group[[col]])
  colors <- brewer.pal(8, name = "Set2")
  names(colors) <- levels
  annotation_colors[[col]] <- colors[annotation_Group[[col]]]
}


# Z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
data_subset_norm <- t(apply(highly_variable_lcpm, 1, cal_z_score))


col <-  c("cadetblue","white", c(brewer.pal(5, "Purples")))

pdf("heatmap_50.pdf", width = 18, height = 10)
pheatmap(data_subset_norm, cluster_rows = T, cluster_cols = T, color = col, annotation_col = annotation_Group,annotation_colors = annotation_colors, clustering_method = "ward.D2", main = "", show_colnames = T, show_rownames = T)
pheatmap(highly_variable_lcpm, cluster_rows = T, cluster_cols = T, color = brewer.pal(8, "BuPu"), annotation_col = annotation_Group, annotation_colors = annotation_colors, clustering_method = "ward.D2", show_colnames = T, show_rownames = T, main = "")
dev.off()


## 6.3. Heatmap highly variable genes (option 2)----

hcopt <- function(d, HC=NULL, method = "ward.D2", members = NULL){
  require("cba")
  if ( is.null(HC) ) {
    HC <- hclust(d,method=method,members=members)
  }
  #optimal leaf ordering
  ORD <- cba::order.optimal(d,merge=HC$merge)
  HC$merge <- ORD$merge
  HC$order <- ORD$order
  HC
}

## use Eculidean distance for columns/samples
## use ward as agglomeration rule
hc01.col <- hcopt(dist(t(highly_variable_lcpm)),method="ward.D2")

## use 1-correlation as distance for for rows/genes
## use ward as agglomeration rule
hc01.row <- hcopt(as.dist(1-cor(t(highly_variable_lcpm))),method="ward.D2")


pheatmap(highly_variable_lcpm,
         color= brewer.pal(8, "BuPu"),
         annotation_col = annotation_Group,
         annotation_colors = annotation_colors,
         cluster_rows=hc01.row,
         cluster_cols=hc01.col,
         show_rownames = T,
         show_colnames = T,
         scale = "row")


## 6.4. Consensus clustering ----
library(ConsensusClusterPlus)

## next, run consensus clustering
CCout <- ConsensusClusterPlus(as.matrix(pseudo_TMM),maxK=10,reps=1000,pItem=0.8,pFeature=1,
                              innerLinkage="ward.D2", finalLinkage="ward.D2",
                              title="cc",clusterAlg="hc",distance="euclidean",seed=1262118388.71279,
                              plot=NULL, writeTable = T )

#To determine the number of clusters based on consensus clustering, we use the consensus distribution, and the proportion change in the area under the consensus CDF (see Section 3.3.1 of ConsensusClusterPlus documentation). The selected K will correspond to the number of clusters where the CDF levels off and the corresponding alpha(K) gets close to zero.

## specify the number of clusters
nc <-4

## remake heatmap, include both subtype and cluster assignments for visual comparison
annotation_Group <- data.frame(Diagnosis = pData_f$Sample_Diagnosis, Source = pData_f$Source_Year, Group_Mutation_type = pData_f$Tumor_Type_Mut, Gender = pData_f$Sex, Development = pData_f$Development_stage, MAPK_Mutated = pData_f$MAPK_mutated)
rownames(annotation_Group) <- (pData_f$Sample_Names)


annot1 <- data.frame(annotation_Group, cluster=as.factor(CCout[[nc]]$consensusClass))
annotation_Group <- annot1
annotation_Group <- annotation_Group[order(annotation_Group$cluster),]
colnames(annotation_Group) <- gsub("_", " ", colnames(annotation_Group))
annotation_Group$`Group Mutation type` <- gsub("_", " ", annotation_Group$`Group Mutation type`)
annotation_colors <- list()
for (col in colnames(annotation_Group)) {
  levels <- unique(annotation_Group[[col]])
  colors <- c(brewer.pal(10,"Set2"), brewer.pal(10,"Set1"))
  names(colors) <- levels
  annotation_colors[[col]] <- unlist(colors[annotation_Group[[col]]])
}


consensus_matrix <- (CCout[[nc]]$ml)
rownames(consensus_matrix) <- names(CCout[[nc]][["consensusClass"]])
colnames(consensus_matrix) <- names(CCout[[nc]][["consensusClass"]])
consensus_matrix <- consensus_matrix[rownames(annotation_Group),rownames(annotation_Group)]

# Define the callback function
callback <- function(hc, mat) {
  sv <- svd(t(mat))$v[, 1]
  dend <- reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

# Create the hclust object
hc <- hclust(dist(consensus_matrix))

# Use pheatmap with the callback function
pheatmap(
  consensus_matrix,
  color = colorRampPalette(c("white", "#542788"))(256),
  annotation_col = annotation_Group,
  annotation_colors = annotation_colors,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  border_color = NA,
  clustering_callback = function(...) callback(hc, consensus_matrix)
)



## 6.5. Heatmap of the sample-to-sample distances ----
#A heatmap of this distance matrix gives us an overview over similarities and dissimilarities between samples. We have to provide a hierarchical clustering hc to the heatmap function based on the sample distances, or else the heatmap function would calculate a clustering based on the distances between the rows/columns of the distance matrix.  Poisson dissimilarity, this measure of dissimilarity also takes the variance structure of counts into consideration when calculating the distances between samples. The PoissonDistance function takes the original count matrix (not normalized) with samples as rows instead of columns
library("PoiClaClu")
poisd <- PoissonDistance(t(raw_f))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- pData_f$Sample_Names
colnames(samplePoisDistMatrix) <- pData_f$Sample_Names
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=brewer.pal(9,"BuPu"), fontsize_col = 9.5, fontsize_row = 9.5)


## 6.6. Volcano plots ----


# List of data frame names
vol_plots <- c(ls(pattern = "^de_genes_"))[-1]

# Loop through the list
for (i in seq_along(vol_plots)) {
  # Get the dataframe
  df <- (mget(vol_plots[i])[[1]][,c(1,4,5)])
  colnames(df) <- c("LogFC", "Pvalue", "FDR")
  
  top10 <- df %>% 
    arrange(desc(abs(LogFC)), FDR) %>% 
    filter(abs(LogFC) > 1.5, FDR < 0.05) %>% 
    slice(1:20)
  
  # Generate the keyvals vector
  keyvals <- ifelse(
    df$LogFC > 2 & df$FDR < 0.05 , '#751094',
    ifelse(df$FDR < 0.05 & df$LogFC > 1, '#ca87de',
           ifelse(df$FDR < 0.05 & df$LogFC < -2, '#007863',
                  ifelse(df$FDR < 0.05 & df$LogFC < -1, '#67ab9f',
                         ifelse(df$Pvalue < 0.05 & df$LogFC > 1, '#e7b0f7',
                                ifelse(df$Pvalue < 0.05 & df$LogFC < -1, '#99f0e1',
                                       ifelse(df$Pvalue < 0.05, 'grey50', 'grey20')))))))
  
  keyvals[is.na(keyvals)] <- 'grey20'
  names(keyvals)[keyvals == '#ca87de'] <- 'Up: FDR < 0.05, logFC > 1'
  names(keyvals)[keyvals == 'grey20'] <- 'NS'
  names(keyvals)[keyvals == 'grey50'] <- 'p-value < 0.05'
  names(keyvals)[keyvals == '#67ab9f'] <- 'Down: FDR < 0.05, logFC < -1'
  names(keyvals)[keyvals == '#751094'] <- 'Up: FDR < 0.05, logFC > 2'
  names(keyvals)[keyvals == '#007863'] <- 'Down: FDR < 0.05, logFC < -2 '
  names(keyvals)[keyvals == '#e7b0f7'] <- 'Up: p-value < 0.05, logFC > 1'
  names(keyvals)[keyvals == '#99f0e1'] <- 'Down: p-value < 0.05, logFC < -1'
  
  # Additional coloring for FDR < 0.05
  
  
  names_vol <- vol_plots[i]
  names_vol <- sub("de_genes_", "", names_vol)
  names_vol2 <- gsub("_", " ", names_vol)
  
  # Generate the volcano plot with only top expressed genes labeled
  volcano_plot <- EnhancedVolcano(df,
                                  lab = rownames(df),
                                  selectLab = rownames(top10),
                                  x = 'LogFC',
                                  y = 'Pvalue',
                                  xlab = bquote(~Log[2]~ 'fold change'),
                                  ylim = c(0, -log10(min(df$Pvalue))),
                                  title = names_vol2, subtitle = "",
                                  pCutoff = 0.05,
                                  FCcutoff = 1,
                                  colCustom = keyvals,
                                  legendPosition = 'right',
                                  gridlines.major = FALSE,
                                  gridlines.minor = FALSE,
                                  cutoffLineCol = 'black',
                                  colAlpha = 1, pointSize = 2, drawConnectors = TRUE,
                                  arrowheads = FALSE,
                                  widthConnectors = 0.5,
                                  labSize = 2.5,
                                  colConnectors = 'black'
  ) + 
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(size = 16, colour = "black"),
          axis.text.y = element_text(size = 16, colour = "black"),
          axis.title.y = element_text(size = 16, colour = "black"),
          axis.line = element_line(linetype = "dashed", colour = "black")) +
    labs(color = "")
  
  # Save the plot as a PDF file
  pdf(paste0(path_dir, "volcano_plot_", names_vol, ".pdf"), width = 10, height = 10)
  print(volcano_plot)
  dev.off()
}

