# EWAS of treatment response: The code is in our GitHub repository 
# (https://github.com/PGC-PTSD-EWAS/PGC-PTSD-Longitudinal-Analysis). The name of the script is Model2.R.

################################################################################
# Linear model to evaluate CpGs associated with whether patients responded to treatment 
# (Response vs Nonresponse). This is a categorical model:
#   DNAmpost ~ DNAmpre + ResponseGroup + Age + Sex + ΔCell Comp + PC
# Phenotype data must be in wide format (each row = participant, columns = pre/post values)
################################################################################

# Load required libraries
library(lme4)
library(lmerTest)
library(data.table)
library(tidyverse)

# Set working directory
setwd("/yourfolder/")

# Load phenotype file (wide format: one row per participant, multiple columns for pre/post)
pheno <- read.csv("phenotype_wideformat.csv")
pheno <- pheno[, -1]  # Drop first column if it’s an index

############################################################
# Load preprocessed (Combat-normalized) beta values
# Beta matrix format: rows = CpGs, columns = Sample IDs
############################################################
beta.norm <- fread("combat_adjusted_beta.csv", data.table = F)
rownames(beta.norm) <- beta.norm$V1  # Assign CpG IDs to rownames
beta.norm <- beta.norm[, -1]         # Drop CpG ID column

############################################################
# Check if all SampleIDs in phenotype match beta matrix columns
############################################################
all(pheno$SampleID_pre %in% colnames(beta.norm))   # Should return TRUE
all(pheno$SampleID_post %in% colnames(beta.norm))  # Should return TRUE

# Debug mismatches if any
unique_pre <- setdiff(pheno$SampleID_pre, colnames(beta.norm))
unique_post <- setdiff(pheno$SampleID_post, colnames(beta.norm))
unique_pre_data <- pheno[pheno$SampleID_pre %in% unique_pre, ]
unique_post_data <- pheno[pheno$SampleID_post %in% unique_post, ]
print("Unique SampleID_pre values not in beta.norm:")
print(unique_pre)
print("Unique SampleID_post values not in beta.norm:")
print(unique_post)

# Filter phenotype data to keep only rows where pre/post SampleIDs are present in beta.norm
pheno_mod <- pheno[
  pheno$SampleID_pre %in% colnames(beta.norm) & 
    pheno$SampleID_post %in% colnames(beta.norm), 
]

# Confirm cleaned phenotype matches beta matrix
all(pheno_mod$SampleID_pre %in% colnames(beta.norm))
all(pheno_mod$SampleID_post %in% colnames(beta.norm))

############################################################
# Define variables and preprocess phenotype
############################################################

# Define binary outcome variable: 1 = Responder, 0 = Nonresponder
pheno_mod$Response_Group <- ifelse(pheno_mod$PCL_Completion_post <= 10, 1, 0)

# Encode Sex: Male = 0, Female = 1
pheno_mod$Sex <- ifelse(pheno_mod$Sex == "Male", 0, 1)
pheno_mod$Sex <- as.numeric(as.character(pheno_mod$Sex))

# Ensure age variable is numeric
pheno_mod$agePre <- as.numeric(pheno_mod$agePre)

# Calculate change in cell composition between pre/post
cellTypes <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
for (cellType in cellTypes) {
  pre_col <- paste0(cellType, "_pre")
  post_col <- paste0(cellType, "_post")
  diff_col <- paste0(cellType, ".Epic.diff")
  pheno_mod[[diff_col]] <- pheno_mod[[post_col]] - pheno_mod[[pre_col]]
}

# Final list of ΔCell Composition variables to include in the model
cellTypes <- c("CD8T.Epic.diff", "CD4T.Epic.diff", "NK.Epic.diff", "Bcell.Epic.diff", "Mono.Epic.diff")

# Define ancestry principal components
pcs <- c("Comp.1_pre", "Comp.2_pre")

# Define study ID (used in file output)
studyID <- "studynamehere"

# Ensure SampleID columns are character
pheno$SampleID_pre <- as.character(pheno$SampleID_pre)
pheno$SampleID_post <- as.character(pheno$SampleID_post)

############################################################
# Convert beta values to M-values for analysis
############################################################
range(beta.norm, na.rm = T)

# Adjust boundary values to avoid log(0) issues
if (min(beta.norm, na.rm = T) == 0) {
  beta.norm[which(beta.norm == 0)] <- 0.0001
}
if (max(beta.norm, na.rm = T) == 1) {
  beta.norm[which(beta.norm == 1)] <- 0.9999
}

# Convert to M-values
beta.norm <- log2(beta.norm / (1 - beta.norm))
range(beta.norm, na.rm = T)

############################################################
# Build model formula: DNAmpost ~ DNAmpre + covariates
############################################################
formula <- as.formula(
  paste("outCpG ~", 
        paste(c("expCpG", ptsdVar, "agePre", paste(cellTypes, collapse = "+"), 
                paste(pcs, collapse = "+"), "Sex"), collapse = " + "))
)

print(formula)

# Store variable names for results
vars <- c("(Intercept)", "expCpG", ptsdVar, "agePre", cellTypes, pcs, "Sex")

# Initialize result matrices
resultsBeta <- matrix(nrow = nrow(beta.norm), ncol = length(vars))
rownames(resultsBeta) <- rownames(beta.norm)
colnames(resultsBeta) <- vars

# Create matching structures for SEs, t-values, p-values
resultsSE <- resultsT <- resultsP <- resultsBeta
resultsDF <- matrix(nrow = nrow(beta.norm), ncol = 1)
rownames(resultsDF) <- rownames(beta.norm)
colnames(resultsDF) <- "df"

# Store CpG IDs
cpgs <- rownames(beta.norm)
errorProbes <- NULL

############################################################
# Run linear model for each CpG site
############################################################
start <- proc.time()[3]

for (ii in 1:length(cpgs)) {
  tempPheno <- pheno_mod
  outCpG <- t(beta.norm[cpgs[ii], pheno_mod$SampleID_post])  # post-treatment methylation
  expCpG <- t(beta.norm[cpgs[ii], pheno_mod$SampleID_pre])   # pre-treatment methylation
  
  fit <- lm(formula, data = tempPheno)
  res <- coef(summary(fit))
  
  resultsBeta[ii, ] <- res[, "Estimate"]
  resultsSE[ii, ] <- res[, "Std. Error"]
  resultsT[ii, ] <- res[, "t value"]
  resultsP[ii, ] <- res[, "Pr(>|t|)"]
  resultsDF[ii, ] <- fit$df.residual
  
  if (ii %% 100 == 0) { print(ii) }  # Progress tracker
}

end <- proc.time()[3]
end - start  # Execution time

############################################################
# Format and save results for Response_Group predictor
############################################################
beta <- as.data.frame(resultsBeta)
pval <- as.data.frame(resultsP)
se <- as.data.frame(resultsSE)
t <- as.data.frame(resultsT)

final <- as.data.frame(cbind(beta[ptsdVar], se[ptsdVar], t[ptsdVar], pval[ptsdVar], resultsDF))
rownames(final) <- rownames(beta)
colnames(final) <- c("BETA", "SE", "t", "pval", "df")

# Write results to file
write.csv(final, file = paste0(studyID, "results.csv"), quote = F, row.names = T)
