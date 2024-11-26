# PTSD EWAS Analysis

This repository contains the pipeline and resources for performing Epigenome-Wide Association Studies (EWAS) on PTSD-related data. The codebase supports both **continuous** and **categorical** variable models, enabling flexible and comprehensive analysis workflows.

## Repository Structure

- **`Additional_files/`**  
  Contains supplementary files required for the analysis, such as metadata, reference data, or supporting documents.

- **`Continuous_model/`**  
  Code and resources for performing EWAS with continuous variables (e.g., PTSD symptom severity).

  - **Key Scripts**:
    - `Step1_Enmix_QC_pipeline.R`: Initial quality control using the Enmix framework.
    - `Step2_Calculate_cellcomposition.R`: Calculation of cell composition estimates for downstream analysis.
    - `Step3.5_Calculate_ancestry_mPCs.R`: Calculation of ancestry-related principal components for covariate adjustment.
    - `Step3_Combat_normalization_RUSH.R`: Data normalization using the Combat method to mitigate batch effects.

- **`Dichotomous_model/`**  
  Scripts for conducting EWAS with categorical variables (e.g., PTSD case/control status).

## How to Use the Code
#### Continuous Variable Model
1. Navigate to the `Continuous_model/` directory.
2. Run the scripts sequentially:
   - `Step1_Enmix_QC_pipeline.R` → Quality control.
   - `Step2_Calculate_cellcomposition.R` → Estimate cell types.
   - `Step3.5_Calculate_ancestry_mPCs.R` → Calculate ancestry principal components.
   - `Step3_Combat_normalization_RUSH.R` → Normalize data for batch effects.

#### Categorical Variable Model
1. Navigate to the `Dichotomous_model/` directory.
2. Follow the provided scripts for analyzing categorical variables such as case-control comparisons.

### Generating Summary Statistics
- Use the `SummaryStats_phenotypefile.R` script to produce descriptive statistics for your data.
