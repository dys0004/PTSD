# EWAS script
# Code adapted from GitHub repository: 
# https://github.com/PGC-PTSD-EWAS/PGC-PTSD-Longitudinal-Analysis
# Model 2: Continuous model
# Formula: ∆Beta ~ ∆PCL + ∆Cell Composition + Age + Sex + Ancestry PCs
################################################################################

# Load required libraries
library(lme4)
library(lmerTest)
library(data.table)
library(tidyverse)

############################################################
# Load phenotype and methylation data
############################################################

setwd("/Volumes/som-ts/GynOb/SmithLab/NewData/People/Dasia/RUSH/Step4_EWAS_DeltaBeta_DeltaPCL_DeltaCellComp_Age_Sex_RUSH/withPCs")

# Read phenotype data in wide format (one row per participant, pre/post data in columns)
pheno <- read.csv("wideformat.csv")
pheno <- pheno[, !(names(pheno) == "X")]  # Remove default index column if present

# Load normalized beta values (Combat-adjusted)
beta.norm <- fread("beta_combat.csv", data.table = F)
rownames(beta.norm) <- beta.norm$V1   # Assign CpG names to rownames
beta.norm <- beta.norm[, -1]          # Remove CpG ID column

############################################################
# Verify that phenotype sample IDs match methylation data
############################################################

all(pheno$SampleID_pre %in% colnames(beta.norm))   # Should return TRUE
all(pheno$SampleID_post %in% colnames(beta.norm))  # Should return TRUE

# Identify any mismatches
unique_pre <- setdiff(pheno$SampleID_pre, colnames(beta.norm))
unique_post <- setdiff(pheno$SampleID_post, colnames(beta.norm))

# Output rows with mismatched IDs
unique_pre_data <- pheno[pheno$SampleID_pre %in% unique_pre, ]
unique_post_data <- pheno[pheno$SampleID_post %in% unique_post, ]

print("Unique SampleID_pre values not in beta.norm:")
print(unique_pre)
print("Unique SampleID_post values not in beta.norm:")
print(unique_post)

# Filter phenotype to include only matched samples
pheno_mod <- pheno[
  pheno$SampleID_pre %in% colnames(beta.norm) & 
    pheno$SampleID_post %in% colnames(beta.norm), 
]

# Final check — all sample IDs should now match
all(pheno_mod$SampleID_pre %in% colnames(beta.norm))
all(pheno_mod$SampleID_post %in% colnames(beta.norm))

############################################################
# Prepare Variables for Model
############################################################

# Model: ∆Beta ~ ∆PCL + Age + Sex + ∆Cell Composition + PCs
studyID <- "Your_study_name"  # Set study label for file output

# Get pre- and post-treatment beta values
beta_pre <- beta.norm[, pheno_mod$SampleID_pre]
beta_post <- beta.norm[, pheno_mod$SampleID_post]

# Calculate delta methylation (post - pre) for each CpG
delta_beta <- beta_post - beta_pre

# Calculate PTSD score change (post - pre)
ptsdVar <- "PTSD_Diff"
ptsdDiff <- pheno_mod$PCL_Completion_post - pheno_mod$PCL_Baseline_pre
pheno_mod[[ptsdVar]] <- ptsdDiff

# Ensure Age is numeric
pheno_mod$agePre <- as.numeric(pheno_mod$agePre)

# Recode Sex: Male = 0, Female = 1
pheno_mod$Sex <- ifelse(pheno_mod$Sex == "Male", 0, 1)
pheno_mod$Sex <- as.numeric(as.character(pheno_mod$Sex))

# Calculate differences in cell type proportions (post - pre)
cellTypes <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
for (cellType in cellTypes) {
  pre_col <- paste0(cellType, "_pre")
  post_col <- paste0(cellType, "_post")
  diff_col <- paste0(cellType, ".Epic.diff")
  pheno_mod[[diff_col]] <- pheno_mod[[post_col]] - pheno_mod[[pre_col]]
}

# Update list of ∆Cell Composition covariates (Neu not included in final model)
cellTypes <- c("CD8T.Epic.diff", "CD4T.Epic.diff", "NK.Epic.diff", "Bcell.Epic.diff", "Mono.Epic.diff")

# Ancestry principal components
pcs <- c("Comp.1_pre", "Comp.2_pre")

#########################################################
# Build the model formula
#########################################################

# Dynamically create the model formula
formula <- as.formula(
  paste("delta_beta ~", paste(c(ptsdVar, "agePre", "Sex", paste(cellTypes, collapse = "+"), paste(pcs, collapse = "+")), collapse = " + "))
)

print(formula)

# List of variables used in the model
vars <- c("(Intercept)", ptsdVar, "agePre", "Sex", cellTypes, pcs)

#########################################################
# Initialize result matrices
#########################################################

# Create storage for beta coefficients, SEs, t-values, p-values, and degrees of freedom
resultsBeta <- matrix(NA, nrow = nrow(delta_beta), ncol = length(vars))
rownames(resultsBeta) <- rownames(delta_beta)
colnames(resultsBeta) <- vars

resultsSE <- resultsBeta
resultsT <- resultsBeta
resultsP <- resultsBeta

resultsDF <- matrix(NA, nrow = nrow(delta_beta), ncol = 1)
rownames(resultsDF) <- rownames(delta_beta)
colnames(resultsDF) <- "df"

cpgs <- rownames(delta_beta)  # Store CpG IDs
errorProbes <- NULL

#########################################################
# Loop through each CpG and fit model
#########################################################

start <- proc.time()[3]  # Start timing

for (i in 1:nrow(delta_beta)) {
  tempPheno <- pheno_mod
  
  # Assign delta_beta value for this CpG to each sample
  delta_values <- as.numeric(delta_beta[i, ])
  names(delta_values) <- colnames(delta_beta)
  tempPheno$delta_beta <- delta_values[match(tempPheno$SampleID_post, names(delta_values))]
  
  # Fit linear model
  fit <- lm(formula, data = tempPheno)
  res <- summary(fit)
  
  # Extract results
  resultsBeta[i, ] <- res$coefficients[, "Estimate"]
  resultsSE[i, ] <- res$coefficients[, "Std. Error"]
  resultsT[i, ] <- res$coefficients[, "t value"]
  resultsP[i, ] <- res$coefficients[, "Pr(>|t|)"]
  resultsDF[i, ] <- fit$df.residual
  
  # Progress update
  if (i %% 100 == 0) {
    print(paste("Processed CpG:", i))
  }
}

# End timing
end <- proc.time()[3]
print(paste("Time elapsed:", end - start, "seconds"))

#########################################################
# Save and export final results for PTSD_Diff effect
#########################################################

# Convert result matrices to data frames
beta <- as.data.frame(resultsBeta)
se <- as.data.frame(resultsSE)
t <- as.data.frame(resultsT)
p <- as.data.frame(resultsP)
df <- as.data.frame(resultsDF)

# Create final output with PTSD_Diff statistics only
final <- as.data.frame(cbind(beta[ptsdVar], se[ptsdVar], t[ptsdVar], p[ptsdVar], resultsDF))
rownames(final) <- rownames(beta)
colnames(final) <- c("BETA", "SE", "t", "pval", "df")

# Save results
write.csv(final, file = paste0(studyID, "resultsfile.csv"), quote = F, row.names = T)
