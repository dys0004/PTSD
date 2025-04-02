# EWAS of treatment response
# Original code adapted from: https://github.com/PGC-PTSD-EWAS/PGC-PTSD-Longitudinal-Analysis
# Model 2: Categorical model
# Formula: ∆Beta ~ Response_Group + age + sex + delta cell comp + Ancestry
################################################################################

# Load libraries
library(lme4)
library(lmerTest)
library(data.table)
library(tidyverse)

############################################################
# Load phenotype and methylation data
############################################################

# Set working directory to where the data is stored
setwd("/Volumes/som-ts/GynOb/SmithLab/NewData/People/Dasia/RUSH/Step4_EWAS_DeltaBetaResponseGroup_deltaCellComp_Age_Sex_RUSH/WithPCs")

# Read wide-format phenotype file (one row per participant, columns for pre/post)
pheno <- read.csv("pheno_wideformat.csv")
pheno <- pheno[, !(names(pheno) == "X")]  # Remove unwanted index column
head(pheno)

# Read Combat-normalized beta matrix (CpGs x Samples)
beta.norm <- fread("combat_beta_values.csv", data.table = F)
rownames(beta.norm) <- beta.norm$V1  # Assign CpG names
beta.norm <- beta.norm[, -1]         # Remove CpG ID column

############################################################
# Ensure phenotype IDs match beta sample IDs
############################################################

# Check that sample IDs in phenotype match methylation data
all(pheno$SampleID_pre %in% colnames(beta.norm))
all(pheno$SampleID_post %in% colnames(beta.norm))

# If mismatched, print and inspect them
unique_pre <- setdiff(pheno$SampleID_pre, colnames(beta.norm))
unique_post <- setdiff(pheno$SampleID_post, colnames(beta.norm))

unique_pre_data <- pheno[pheno$SampleID_pre %in% unique_pre, ]
unique_post_data <- pheno[pheno$SampleID_post %in% unique_post, ]

print("Unique SampleID_pre values not in beta.norm:")
print(unique_pre)

print("Unique SampleID_post values not in beta.norm:")
print(unique_post)

# Filter phenotype to only matched samples
pheno_mod <- pheno[
  pheno$SampleID_pre %in% colnames(beta.norm) & 
    pheno$SampleID_post %in% colnames(beta.norm),
]

# Final check after filtering
all(pheno_mod$SampleID_pre %in% colnames(beta.norm))
all(pheno_mod$SampleID_post %in% colnames(beta.norm))

############################################################
# Define variables and calculate delta beta
############################################################

studyID <- "your_study_name"  # Used for naming output file

# Extract pre and post DNAm matrices
beta_pre <- beta.norm[, pheno_mod$SampleID_pre]
beta_post <- beta.norm[, pheno_mod$SampleID_post]

# Calculate ∆Beta (post - pre) for each CpG
delta_beta <- beta_post - beta_pre

############################################################
# Sanity checks for ∆Beta calculation
############################################################

# Ensure CpGs are aligned
if (!all(rownames(beta_pre) == rownames(beta_post))) {
  stop("Mismatch in CpG site ordering!")
}

# Ensure sample columns are properly matched
pre_ids <- colnames(beta_pre)
post_ids <- colnames(beta_post)
expected_post_ids <- pheno_mod$SampleID_post[match(pre_ids, pheno_mod$SampleID_pre)]

if (!all(post_ids == expected_post_ids)) {
  stop("Mismatch in sample ID order!")
} else {
  print("Sample IDs correctly matched.")
}

# Validate ∆Beta calculation
expected_deltas <- beta_post - beta_pre
if (exists("delta_beta")) {
  discrepancy <- sum(abs(delta_beta - expected_deltas) > 1e-6, na.rm = TRUE)
  if (discrepancy == 0) {
    print("Beta differences correctly calculated.")
  } else {
    print(paste("Warning:", discrepancy, "mismatches found!"))
  }
} else {
  print("Delta beta matrix not found — calculating now...")
  delta_beta <- expected_deltas
}

############################################################
# Preprocess phenotype variables
############################################################

# Define outcome variable: binary treatment response
ptsdVar <- "Response_Group"

# Create binary variable: 1 = Responder (PCL ≤ 10), 0 = Nonresponder
pheno_mod$Response_Group <- ifelse(pheno_mod$PCL_Completion_post <= 10, 1, 0)

# Encode Sex: Male = 0, Female = 1
pheno_mod$Sex <- ifelse(pheno_mod$Sex == "Male", 0, 1)
pheno_mod$Sex <- as.numeric(as.character(pheno_mod$Sex))

# Ensure age is numeric
pheno_mod$agePre <- as.numeric(pheno_mod$agePre)

# Calculate delta cell proportions (post - pre) for immune cell types
cellTypes <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
for (cellType in cellTypes) {
  pre_col <- paste0(cellType, "_pre")
  post_col <- paste0(cellType, "_post")
  diff_col <- paste0(cellType, ".Epic.diff")
  pheno_mod[[diff_col]] <- pheno_mod[[post_col]] - pheno_mod[[pre_col]]
}

# Final cell types used (Neu removed)
cellTypes <- c("CD8T.Epic.diff", "CD4T.Epic.diff", "NK.Epic.diff", "Bcell.Epic.diff", "Mono.Epic.diff")

# Ancestry covariates
pcs <- c("Comp.1_pre", "Comp.2_pre")

############################################################
# Build model formula
############################################################

# Linear model: delta_beta ~ Response_Group + age + sex + ∆cellComp + PCs
formula <- as.formula(
  paste("delta_beta ~", 
        paste(c(ptsdVar, "agePre", "Sex", paste(cellTypes, collapse = "+"), paste(pcs, collapse = "+")), 
              collapse = " + "))
)

print(formula)

# List of variables for result matrices
vars <- c("(Intercept)", ptsdVar, "agePre", "Sex", cellTypes, pcs)

############################################################
# Initialize result matrices
############################################################

resultsBeta <- matrix(NA, nrow = nrow(delta_beta), ncol = length(vars))
rownames(resultsBeta) <- rownames(delta_beta)
colnames(resultsBeta) <- vars

resultsSE <- resultsT <- resultsP <- resultsBeta  # Same shape
resultsDF <- matrix(NA, nrow = nrow(delta_beta), ncol = 1)
rownames(resultsDF) <- rownames(delta_beta)
colnames(resultsDF) <- "df"

cpgs <- rownames(delta_beta)
errorProbes <- NULL

############################################################
# Run model across all CpGs
############################################################

start <- proc.time()[3]  # Start timing

for (i in 1:nrow(delta_beta)) {
  tempPheno <- pheno_mod
  
  # Assign CpG-specific delta_beta to phenotype
  delta_values <- as.numeric(delta_beta[i, ])
  names(delta_values) <- colnames(delta_beta)
  tempPheno$delta_beta <- delta_values[match(tempPheno$SampleID_post, names(delta_values))]
  
  # Fit linear model
  fit <- lm(formula, data = tempPheno)
  res <- summary(fit)
  
  # Store results
  resultsBeta[i, ] <- res$coefficients[, "Estimate"]
  resultsSE[i, ] <- res$coefficients[, "Std. Error"]
  resultsT[i, ] <- res$coefficients[, "t value"]
  resultsP[i, ] <- res$coefficients[, "Pr(>|t|)"]
  resultsDF[i, ] <- fit$df.residual
  
  # Progress update every 100 CpGs
  if (i %% 100 == 0) {
    print(paste("Processed CpG:", i))
  }
}

end <- proc.time()[3]
print(paste("Time elapsed:", end - start, "seconds"))

############################################################
# Extract final results for Response_Group
############################################################

# Convert results to data frames
beta <- as.data.frame(resultsBeta)
se <- as.data.frame(resultsSE)
t <- as.data.frame(resultsT)
p <- as.data.frame(resultsP)
df <- as.data.frame(resultsDF)

# Combine relevant results for Response_Group
final <- as.data.frame(cbind(beta[ptsdVar], se[ptsdVar], t[ptsdVar], p[ptsdVar], resultsDF))
rownames(final) <- rownames(beta)
colnames(final) <- c("BETA", "SE", "t", "p", "df")

# Export results
write.csv(final, file = "your_results.csv", quote = F, row.names = T)
