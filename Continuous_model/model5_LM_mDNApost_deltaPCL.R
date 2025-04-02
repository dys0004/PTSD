# EWAS of treatment response
# Original source: https://github.com/PGC-PTSD-EWAS/PGC-PTSD-Longitudinal-Analysis
# Formula: mDNA_post ~ mDNA_pre + ∆PCL + ∆Age + ∆Cell Composition + PCs + Sex
################################################################################

# Load libraries
library(lme4)         # Mixed models (not used in this script, but loaded)
library(lmerTest)     # P-values for mixed models
library(data.table)   # Fast reading of large tables
library(tidyverse)    # Tidy data wrangling

# Set working directory
setwd("/your_folder/")

# Load pre-formatted wide phenotype file (each row = individual)
pheno <- read.csv("pheno_wideformat.csv")
pheno <- pheno[, !(names(pheno) == "X")]  # Drop index column if present

# Load normalized beta values (CpG x SampleID)
beta.norm <- fread("Beta_combat.csv", data.table = F)
rownames(beta.norm) <- beta.norm$V1
beta.norm <- beta.norm[, -1]

############################################################
# Define covariates and compute difference variables
############################################################

# Age difference (∆Age = post - pre)
ageVar <- "AgeDiff"
ageDiff <- pheno$agePost - pheno$agePre
pheno[[ageVar]] <- ageDiff

# Cell composition difference (∆CellType = post - pre)
cellTypes <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
for (cellType in cellTypes) {
  pre_col <- paste0(cellType, "_pre")
  post_col <- paste0(cellType, "_post")
  diff_col <- paste0(cellType, ".Epic.diff")
  pheno[[diff_col]] <- pheno[[post_col]] - pheno[[pre_col]]
}
# Drop "Neu" from model
cellTypes <- c("CD8T.Epic.diff", "CD4T.Epic.diff", "NK.Epic.diff", "Bcell.Epic.diff", "Mono.Epic.diff")

# PTSD symptom change
ptsdVar <- "PTSD_Diff"
ptsdDiff <- pheno$PCL_Completion_post - pheno$PCL_Baseline_pre
pheno[[ptsdVar]] <- ptsdDiff

# Ancestry PCs
pcs <- c("Comp.1_pre", "Comp.2_pre")

# Sex: recode as binary (0 = Male, 1 = Female)
pheno$Sex <- ifelse(pheno$Sex == "Male", 0, 1)
pheno$Sex <- as.numeric(as.character(pheno$Sex))

# Study name (used in output filename)
studyID <- "your_folder"

# Ensure Sample IDs are treated as characters
pheno$SampleID_pre <- as.character(pheno$SampleID_pre)
pheno$SampleID_post <- as.character(pheno$SampleID_post)

############################################################
# Match sample IDs in phenotype and methylation matrix
############################################################

# Confirm pre/post sample IDs exist in beta matrix
all(pheno$SampleID_pre %in% colnames(beta.norm))
all(pheno$SampleID_post %in% colnames(beta.norm))

# If mismatch, show which IDs are missing
unique_pre <- setdiff(pheno$SampleID_pre, colnames(beta.norm))
unique_post <- setdiff(pheno$SampleID_post, colnames(beta.norm))
unique_pre_data <- pheno[pheno$SampleID_pre %in% unique_pre, ]
unique_post_data <- pheno[pheno$SampleID_post %in% unique_post, ]

print("Unique SampleID_pre values not in beta.norm:")
print(unique_pre)
print("Unique SampleID_post values not in beta.norm:")
print(unique_post)

# Keep only samples present in both beta and phenotype
pheno_mod <- pheno[
  pheno$SampleID_pre %in% colnames(beta.norm) & 
    pheno$SampleID_post %in% colnames(beta.norm),
]

# Final ID check
all(pheno_mod$SampleID_pre %in% colnames(beta.norm))
all(pheno_mod$SampleID_post %in% colnames(beta.norm))

############################################################
# Convert beta values to M-values
############################################################

range(beta.norm, na.rm = TRUE)

# Fix boundary values to avoid log(0)
if (min(beta.norm, na.rm = TRUE) == 0) {
  beta.norm[which(beta.norm == 0)] <- 0.0001
}
if (max(beta.norm, na.rm = TRUE) == 1) {
  beta.norm[which(beta.norm == 1)] <- 0.9999
}

# Logit transform: beta to M-values
beta.norm <- log2(beta.norm / (1 - beta.norm))
range(beta.norm, na.rm = TRUE)

############################################################
# Build model formula
############################################################

# Model: mDNA_post ~ mDNA_pre + ∆PCL + ∆Age + ∆CellTypes + PCs + Sex
formula <- as.formula(
  paste("outCpG ~", 
        paste(c("expCpG", ptsdVar, ageVar, paste(cellTypes, collapse = "+"), paste(pcs, collapse = "+"), "Sex"),
              collapse = " + "))
)
print(formula)

# Variable names used in results matrices
vars <- c("(Intercept)", "expCpG", ptsdVar, cellTypes, ageVar, pcs, "Sex")

############################################################
# Initialize result matrices
############################################################

resultsBeta <- matrix(nrow = nrow(beta.norm), ncol = length(vars))
rownames(resultsBeta) <- rownames(beta.norm)
colnames(resultsBeta) <- vars

resultsSE <- resultsT <- resultsP <- resultsBeta
resultsDF <- matrix(nrow = nrow(beta.norm), ncol = 1)
rownames(resultsDF) <- rownames(beta.norm)
colnames(resultsDF) <- "df"

cpgs <- rownames(beta.norm)
errorProbes <- NULL

############################################################
# Fit linear model for each CpG site
############################################################

start <- proc.time()[3]  # Start timing

for (ii in 1:length(cpgs)) {
  tempPheno <- pheno_mod
  
  # Pull pre- and post-treatment M-values for current CpG
  outCpG <- t(beta.norm[cpgs[ii], pheno_mod$SampleID_post])
  expCpG <- t(beta.norm[cpgs[ii], pheno_mod$SampleID_pre])
  
  # Assign them to phenotype dataframe
  tempPheno$outCpG <- outCpG
  tempPheno$expCpG <- expCpG
  
  # Fit model
  fit <- lm(formula, data = tempPheno)
  res <- coef(summary(fit))
  
  # Store results
  resultsBeta[ii, ] <- res[, "Estimate"]
  resultsSE[ii, ] <- res[, "Std. Error"]
  resultsT[ii, ] <- res[, "t value"]
  resultsP[ii, ] <- res[, "Pr(>|t|)"]
  resultsDF[ii, ] <- fit$df.residual
  
  # Progress update
  if (ii %% 100 == 0) {
    print(ii)
  }
}

end <- proc.time()[3]
end - start  # Print time taken

############################################################
# Extract results 
############################################################

beta <- as.data.frame(resultsBeta)
se <- as.data.frame(resultsSE)
t <- as.data.frame(resultsT)
pval <- as.data.frame(resultsP)

final <- as.data.frame(cbind(beta[ptsdVar], se[ptsdVar], t[ptsdVar], pval[ptsdVar], resultsDF))
rownames(final) <- rownames(beta)
colnames(final) <- c("BETA", "SE", "t", "pval", "df")

# Save output
write.csv(final,
          file = paste0(studyID, "EWASresults.csv"),
          quote = FALSE, row.names = TRUE)
