# PTSD EWAS Analysis

This repository contains the pipeline and resources for performing Epigenome-Wide Association Studies (EWAS) on PTSD-related data. The codebase supports both **continuous** and **categorical** variable models. (https://github.com/PGC-PTSD-EWAS/PGC-PTSD-Longitudinal-Analysis)

## Repository Structure

- **`Additional_files/`**  
  Contains supplementary files required for the analysis.

  - **Key Scripts**:
    - `Step1_Enmix_QC_pipeline.R`: Initial quality control using the Enmix framework.
    - `Step2_Calculate_cellcomposition.R`: Calculation of cell composition estimates for covariate adjustment.
    - `Step3.5_Calculate_ancestry_mPCs.R`: Calculation of ancestry-related principal components for covariate adjustment.
    - `Step3_Combat_normalization_RUSH.R`: Data normalization using the Combat method to minimize batch effects.

- **`Dichotomous_model/`**  
  Scripts for conducting EWAS with categorical variables (e.g., PTSD case, Responder/NonResponder).
  
- **`Continuous_model/`**  
  Code and resources for performing EWAS with continuous variables (e.g., PTSD symptom severity, PCL Score).

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
