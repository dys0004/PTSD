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


---

## Model Reference

| Type        | Formula                                                                 |
|-------------|-------------------------------------------------------------------------|
| Continuous  | ∆Beta ~ ∆PCL + ∆Cell Comp + Age + Sex + PCs                             |
| Continuous  | DNAm ~ PCL Score + Age + Sex + Cell Comp + PCs + (1 | Participant ID)   |
| Continuous  | mDNApost ~ mDNApre + ∆PCL + ∆Age + ∆Cell Composition + PCs + Sex        |
| Categorical | ∆Beta ~ ResponseGroup + Age + Sex + ∆Cell Comp + PCs                    |
| Categorical | DNAm_post ~ DNAm_pre + ResponseGroup + Age + Sex + ∆Cell Comp + PCs     |

> **Note:** All models use **PCL** (PTSD Checklist) as the PTSD variable.  
> You can modify this to use a different PTSD measure (e.g., CAPS) by editing the variable names in both the code and your phenotype file.


---

## Variable Definitions

| Variable        | Description                                                                 |
|-----------------|-----------------------------------------------------------------------------|
| `PCL`           | PTSD Checklist total score. Used as the primary PTSD symptom variable.      |
| `ResponseGroup` | Binary variable indicating treatment response: `1` = Responder (PCL ≤ 10), `0` = Nonresponder. |
| `PTSD_Diff`     | Change in PTSD symptoms: `PCL_post - PCL_pre`.                             |
| `Age` / `AgeDiff` | Age at time of sample collection, or age difference if using pre/post.     |
| `Sex`           | Encoded as `0` = Male, `1` = Female.                                        |
| `Cell Composition` | Proportions of immune cell types (CD8T, CD4T, NK, Bcell, Mono, Neu). Difference variables are computed as post - pre values. |
| `PCs`           | Ancestry principal components (e.g., `Comp.1_pre`, `Comp.2_pre`).           |
| `SampleID_pre` / `SampleID_post` | Sample identifiers for pre- and post-treatment timepoints.   |
| `mDNApre` / `mDNApost` | Methylation values before and after treatment (logit-transformed M-values). |

> You can customize variable names and thresholds (e.g., for `ResponseGroup` or `PCs`) in your phenotype file and the scripts.


