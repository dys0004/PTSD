# PTSD EWAS Analysis

This repository contains the pipeline and resources for performing Epigenome-Wide Association Studies (EWAS) on PTSD-related data. The codebase supports both **continuous** and **categorical** variable models. Most code adapted from (https://github.com/PGC-PTSD-EWAS/PGC-PTSD-Longitudinal-Analysis)

## Repository Structure

- **`Additional_files/`**  
  Contains supplementary files required for the analysis.

  - **Key Scripts**:
    - `Step1_Enmix_QC_pipeline.R`: Initial quality control using the Enmix framework.
    - `Step2_Calculate_cellcomposition.R`: Calculation of cell composition estimates for covariate adjustment.
    - `Step3.5_Calculate_ancestry_mPCs.R`: Calculation of ancestry-related principal components for covariate adjustment.
    - `Step3_Combat_normalization_RUSH.R`: Data normalization using the Combat method to minimize batch effects.

- **`Dichotomous_model/`**  
  EWAS models with **categorical** variables (e.g., treatment responders vs nonresponders).

  - `model2_LM_responseGroup.R`: Linear model evaluating treatment response (binary) using post-treatment DNAm.
  - `model4_LM_responseGroup_deltaBeta.R`: Linear model evaluating whether responders show greater DNAm change (∆Beta) from pre to post.

- **`Continuous_model/`**  
  EWAS models with **continuous** variables (e.g., PTSD severity scores or clinical change).

  - `model1_RE_PCL.R`: Random effects model testing PTSD symptom severity (PCL) and baseline DNAm.
  - `model3_LM_deltaPCL.R`: Linear model testing ∆PCL (symptom change) as a predictor of ∆Beta.
  - `model5_LM_mDNApost_deltaPCL.R`: Linear model testing whether symptom change predicts post-treatment DNAm, adjusting for pre-DNAm and other covariates.
 
 🔄 **Customize:** The PTSD variable in these scripts is set to **PCL**.
 You can replace it with your preferred variable (e.g., CAPS, symptom domain scores).

## How to Use the Code
Input file: Raw IDATS, phenotype file
1. Run the scripts sequentially:
   - `Step1_Enmix_QC_pipeline.R` → Quality control.
   - `Step2_Calculate_cellcomposition.R` → Estimate cell types.
   - `Step3.5_Calculate_ancestry_mPCs.R` → Calculate ancestry principal components.
   - `Step3_Combat_normalization_RUSH.R` → Normalize data for batch effects.

#### Variable Model
1. Navigate to the `Dichotomous_model/` directory OR `Continuous_model/`
2. Follow the provided scripts for analysis

### Generating Summary Statistics
- Use the `SummaryStats_phenotypefile.R` script to produce descriptive statistics for your phenotype data.
