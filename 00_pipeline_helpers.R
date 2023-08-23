RUV_total <- function(raw,pData,fData,k,hkgenes = NULL,exclude = NULL){
  
  ### INPUT: raw - p x n raw expressions with p genes and n samples
  ###        pData - phenotype metadata across samples
  ###        fData - feature metadata across genes
  ###        k - number of dimensions of unwanted variation estimated
  ###        exclude - vector of gene names to exclude
  
  library(RUVSeq)
  library(DESeq2)
  library(limma)
  library(matrixStats)
  
  if (!is.null(hkgenes)){
    
    fData(set)$Class[rownames(set) %in% hkgenes] = 'Housekeeping'
    
  }
  
  fData = fData[rownames(raw),]
  int = intersect(rownames(raw),rownames(fData))
  fData = fData[int,]
  raw = raw[int,]
  
  set <- newSeqExpressionSet(as.matrix(round(raw)),
                             phenoData=pData,
                             featureData=fData)
  
  cIdx <- rownames(set)[fData(set)$Class == "Housekeeping"]
  cIdx = cIdx[!(cIdx %in% exclude)]
  x <- as.factor(pData$Group)
  set <- betweenLaneNormalization(set, which="upper")
  set <- RUVg(set, cIdx, k=k)
  dds <- DESeqDataSetFromMatrix(counts(set),colData=pData(set),design=~1)
  rowData(dds) <- fData
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds)
  cts <- counts(dds, normalized=TRUE)
  disp <- pmax((rowVars(cts) - rowMeans(cts)),0)/rowMeans(cts)^2
  mcols(dds)$dispGeneEst <- disp
  dds <- estimateDispersionsFit(dds, fitType="mean")
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  mat <- assay(vsd)
  covars <- as.matrix(colData(dds)[,grep("W",colnames(colData(dds))),drop=FALSE])
  mat <- removeBatchEffect(mat, covariates=covars)
  assay(vsd) <- mat
  return(list(set = set,vsd = vsd))
  
  
  
}

  # QC Information. Return to NanoString manual for more information.
  
# Imaging QC: f fields of view(FOVs) the Digital Analyzer or Sprintwas able to capture. At least 75% of FOVs should be successfully counted to obtain robust data. An issue associated with the instrumentation
  #This percentage should be super high, close to 100%. If there isn't that usually means the cartridge needs to be rescanned. If it is < 3 weeks from the scan, rescanning should be okay. Nanostring folks call < 75% a failure.
  
# Binding Density QC: The ideal range 0.1 - 2.25 spots per square micron has been established for assays run on an nCounter MAX or FLEX system and 0.1 - 1.8 spots per square micron on the nCounter SPRINT system. Measurements outside of these ideal ranges will be flagged, but should be checked to see how much they deviate from the ideal range. If they are only slightly outside of range, they do not indicate a problem in the data and can be bypassed.
  #Nanostring folks call a sample problematic if the binding density is < 0.05 or > 2.25.
  
# Positive Control QC Parameters:
  # Overall assay efficiency. nSolver raises a warning flag when the geometric mean of positive controls is more than three-fold different from the geometric mean of all samples.
  # Assay linearity. Decreasing linear counts are expected from POS_A to POS_E (POS_F is considered below the limit of detection).
  # Limit of detection (LOD). It is expected that counts for POS_E will be higher than background, which is represented by the mean of the negative controls plus two standard deviations (for most assays) or simply the mean of the negative controls (for miRNA assays).
# NanoString recommends positive scaling factors between 0.3 and 3, and R-squared values greater than 0.95.

# Housekeeping normalization scale factors. CodeSet Content Normalization factor. Check which samples are outside the range of 0.1-10

# Ligation QC: the ligation-negative controls will yield counts in the negative control range (background level) and the ligation-positive control counts will be significantly higher, increasing from LIG_POS_C to LIG_POS_A.

imagingQC <- function(rcc){
  
  
  #### INPUT: rcc - input from rcc
  #### OUTPUT: flag for imaging quality
  
  fovRatio = as.numeric(rcc$Lane_Attributes[3]) / as.numeric(rcc$Lane_Attributes[2])
  
  if (!(fovRatio > .75)) {
    flag <- "Flag"
  } else {
    flag <- "No flag"
  }
  
  
  return(list(fovRatio = fovRatio, flag = flag))
  
}



pos_QC <- function(nano_data, i){
  require(NanoTube)
  dat_qc <- processNanostringData(nano_data, output.format = "list",
                                  includeQC = TRUE)
  bd_qc  <- as.data.frame(dat_qc$qc)
  
  posQC <- positiveQC(dat_qc)
  scale_factor <- posQC$tab[i,2]
  r_squared <- posQC$tab[i,3]
  if ( scale_factor > 3 | scale_factor < 0.3 | r_squared < 0.95) {
    flag <- "Flag"
  } else {
    flag <- "No flag"
  }
  
  cat("\n" , "\n" )
  return(list(pos_scale_factor = scale_factor, pos_r_squared = r_squared, flag = flag))
  
}


pos_HK_QC <- function(nano_data) {
  require(NanoTube)
  dat_qc <- processNanostringData(nano_data, output.format = "list", includeQC = TRUE)
  posQC <- positiveQC(dat_qc)
  SampleID <- posQC$tab[,1]
  scale_factor <- posQC$tab[,2]
  r_squared <- posQC$tab[,3]
  hk_scale_factor <-dat_qc$hk.scalefactors
  PosQC_flag <- ifelse(scale_factor > 3 | scale_factor < 0.3 | r_squared < 0.95, "Flag", "No flag")
  HK_flag <- ifelse(dat_qc$hk.scalefactors <= 0.1 | dat_qc$hk.scalefactors >= 10, "Flag", "No flag")
  
  result <- data.frame(SampleID = SampleID,
                       pos_scale_factor = scale_factor,
                       pos_r_squared = r_squared,
                       PosQC_flag = PosQC_flag,
                       hk_scale_factor = hk_scale_factor,
                       HK_flag = HK_flag)
}



bindingDensityQC <- function(rcc, low = 0.05, high = 2.25) {
  
  # extract binding density
  bd <- as.numeric(rcc[["Lane_Attributes"]][["BindingDensity"]])
  
  # check if bd is within specified limits
  if (bd < low | bd > high) {
    flag <- "Flag"
  } else {
    flag <- "No flag"
  }
  
  return(list(bd = bd, flag = flag))
}


limitOfDetectionQC <- function(rcc,numSD = 0){
  
  #### INPUT: rcc - input from rcc
  ####         numSD - number of standard deviations to calibrate the LOD
  #### OUTPUT: flag for limit of detection
  
  counts = rcc$Code_Summary
  posE = as.numeric(counts$Count[counts$Name == 'POS_E'])
  posF = as.numeric(counts$Count[counts$Name == 'POS_F'])
  negControls = as.numeric(counts$Count[grepl('NEG',counts$Name)])
  if(!(posE > mean(negControls) + numSD*sd(negControls))) {flag <- "Flag"}
  if (posE > mean(negControls) + numSD*sd(negControls)) {flag <- 'No flag'}
  return(list(flag = flag, lod = mean(negControls) + numSD*sd(negControls), pose = posE, posf = posF))
}



positiveLinQC <- function(rcc) {
  
  # extract counts for positive controls
  counts <- rcc$Code_Summary
  pos_names <- grep("^POS_", counts$Name)
  posControls <- as.numeric(counts$Count[pos_names])
  
  # expected counts for positive controls
  known <- c(128, 128/4, 128/16, 128/64, 128/256, 128/(256*4))
  
  # check linearity
  r2 <- summary(lm(sort(posControls)~sort(known)))$r.squared
  if (is.na(r2) || r2 < 0.95) {
    return('Flag')
  } else {
    return('No flag')
  }
}


#' Assess background expression. (!)(!) THIS CODE HAS BEEN MODIFIED FROM THE NANOTUBE PACKAGE. IN THIS CASE, IT ADMITS ARRAYS AND CAN PERFORM BOTH T.TEST AND THRESHOLD, RETURNING A LIST IN WHICH ALL THE DATA IS FOUND. IT IS NECESSARY TO PROVIDE AN FDATA DATAFRAME WITH THE FOLLOWING DATA: NAME OF THE GENE (THE SAME AS FOUND IN THE COUNTING MATRIX, AS fdat$Class) AND CODE/CLASS OF THE GENE (POSITIVE, NEGATIVE, ENDOGENOUS, SPIKE-IN, HOUSEKEEPING; AS fdat$Gene)
#' Original code: https://rdrr.io/github/calebclass/NanoTube/src/R/remove_background.R
#' Compare endogenous gene expression data against negative control genes and
#' remove data for genes that fail the comparison. This step is
#' conducted within processNanostringData, when normalization is set to 
#' "nCounter" in Nanotube package. 
#' 
#' @export
#'
#' @param dat Positive control-scaled NanoString data
#' @param fdat FDATA DATAFRAME WITH THE FOLLOWING DATA: NAME OF THE GENE (THE SAME AS FOUND IN THE COUNTING MATRIX, AS fdat$Class) AND CODE/CLASS OF THE GENE (POSITIVE, NEGATIVE, ENDOGENOUS, SPIKE-IN, HOUSEKEEPING; AS fdat$Gene)
#' @param mode Either "threshold" (default) or "t.test". If "threshold", 
#' requires proportionReq of samples to have expression numSD standard 
#' deviations among the mean of negative control genes. If "t.test", each gene 
#' will be compared with all negative control genes in a one-sided two-sample 
#' t-test. 
#' @param numSD Number of standard deviations above mean of negative control 
#' genes to used as background threshold for each sample: 
#' mean(negative_controls) + numSD * sd(negative_controls). 
#' Required if mode == "threshold" or subtract == TRUE
#' @param proportionReq Required proportion of sample expressions exceeding the
#' sample background threshold to include gene in further analysis. Required if
#' mode == "threshold" or subtract == TRUE
#' @param pval p-value (from one-sided t-test) threshold to declare gene 
#' expression above background expression level. Genes with p-values above 
#' this level are removed from further analysis. Required if mode == "t.test"
#' @param subtract Should calculated background levels be subtracted from 
#' reported expressions? If TRUE, will subtract mean+numSD*sd of the negative 
#' controls from the endogenous genes, and then set negative values to zero 
#' (default FALSE).
#' 
#' @return List with dataframes cleared with t.test or threshold and the dataframes in which de substraction is based
#' Expression levels are updated for all genes if subtract == TRUE.
#' 
#' @examples
#'
#' # Remove endogenous genes that fail to reject the null hypothesis
#' # in a one-sided t test against negative control genes with p < 0.05.
#' dat <- remove_background(dat, fdat,  mode = "t.test", pval = 0.05)
#' 
#' # Remove endogenous genes where fewer than 25% of samples have an expression
#' # 2 standard deviations above the average negative control gene. Also, 
#' # subtract this background level (mean + 2*sd) from endogenous genes.
#' dat <- remove_background(dat, mode = "threshold", 
#'                          numSD = 2, proportionReq = 0.25, subtract = TRUE)
#'Both threshold and t.test                          
#'dat <- remove_background(dat, fData, mode = c("threshold", "t.test"), numSD = 1, proportionReq = 0.5, pval = 0.05,
#'                          subtract = F)                          

remove_background <- function(dat, fdat,
                              mode = c("threshold", "t.test"),
                              numSD, proportionReq, pval,
                              subtract = FALSE) {
  
  # Check for negative control genes in the data set
  if (sum(fdat$Class == "Negative") == 0) {
    stop("Cannot conduct filtering with negative controls: No negative control genes found in input")
  }
  
  negative.mean <- apply(dat[grep("negative", fdat$Class, ignore.case = TRUE), ], 2, mean)
  negative.sd <- apply(dat[grep("negative", fdat$Class, ignore.case = TRUE), ], 2, sd)
  
  bg.stats <- NULL
  gene.stats <- NULL
  dat.ttest <- NULL
  
  if ("threshold" %in% mode | subtract) {
    bg.threshold <- negative.mean + numSD * negative.sd
    
    bg.stats <- data.frame(Mean.Neg = negative.mean,
                           Max.Neg = apply(dat[grep("negative", fdat$Class, ignore.case = TRUE), ], 2, max),
                           sd.Neg = negative.sd,
                           background = bg.threshold,
                           num.less.bg = colSums(t(t(dat[grep("endogenous", fdat$Class, ignore.case = TRUE), ]) < bg.threshold)),
                           frc.less.bg = colMeans(t(t(dat[grep("endogenous", fdat$Class, ignore.case = TRUE), ]) < bg.threshold)))
  }
  
  if ("t.test" %in% mode) {
    gene.stats <- data.frame(row.names = fdat$Gene,
                             t.stat = rep(NA, length(fdat$Gene)),
                             p.val = rep(NA, length(fdat$Gene)),
                             pass = rep(NA, length(fdat$Gene)))  # Initialize gene.stats with correct number of rows
    
    rowsKeep <- ifelse(grepl("endogenous", fdat$Class, ignore.case = TRUE), yes = FALSE, no = TRUE)
    background <- unlist(dat[grep("negative", fdat$Class, ignore.case = TRUE), ])
    
    for (i in grep("endogenous", fdat$Class, ignore.case = TRUE)) {
      ttest <- t.test(x = dat[i, ], y = background, var.equal = FALSE, alternative = "greater")
      if (ttest$p.value < pval) {
        rowsKeep[i] <- TRUE
      }
      gene.stats$t.stat[i] <- as.numeric(ttest$statistic)
      gene.stats$p.val[i] <- ttest$p.value
    }
    
    gene.stats$pass <- rowsKeep
    dat.ttest <- dat[rowsKeep, ]
  }
  
  if (subtract) {
    for (i in seq_len(ncol(dat))) {
      dat[, i] <- dat[, i] - bg.threshold[i]
      dat[dat < 0] <- 0
    }
  }
  
  rowsKeep <- which(!grepl("endogenous", fdat$Class, ignore.case = TRUE) | rowMeans(t(t(dat) > bg.threshold)) >= proportionReq)
  dat.threshold <- dat[rowsKeep, ]
  
  return(list(dat.threshold = dat.threshold, bg.stats = bg.stats, dat.ttest = dat.ttest, gene.stats = gene.stats))
}


selectGenes <- function(counts, min.count=10, N=0.90){
  
  lib.size <- colSums(counts)
  MedianLibSize <- median(lib.size)
  CPM.Cutoff <- min.count / MedianLibSize*1e6
  CPM <- edgeR::cpm(counts,lib.size=lib.size)
  
  min.samples <- round(N * ncol(counts))
  
  f1 <- genefilter::kOverA(min.samples, CPM.Cutoff)
  flist <- genefilter::filterfun(f1)
  keep <- genefilter::genefilter(CPM, flist)
  
  ## the same as:
  #keep <- apply(CPM, 1, function(x, n = min.samples){
  #  t = sum(x >= CPM.Cutoff) >= n
  #  t
  #})
  
  return(keep)
}


