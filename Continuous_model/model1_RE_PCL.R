# EWAS script 
# (https://github.com/PGC-PTSD-EWAS/PGC-PTSD-Longitudinal-Analysis). The name of the script is Model1.R.

# Code was adjusted for this model 
# continuous model with random effects 
# mDNA ~ PCL Score (or CAPS) + Age + Sex, Cell Comp + Ancestry + (1 | Partcipant ID)

################################################################################
# Model 1: RE model to evaluate CpGs associated baseline and completion pcl scores 
################################################################################

## set working directory 
## Sample cohort : RUSH 
setwd("...")

# Load required libraries
library(lme4)         # for fitting linear mixed-effects models
library(lmerTest)     # to get p-values for fixed effects in lme4
library(data.table)   # fast data reading

# Load methylation beta values (CpG x Participants)
beta.norm<-fread("you_combat_adjusted_beta_values.csv",data.table=F)
rownames(beta.norm)<-beta.norm$V1
beta.norm<-beta.norm[,-1]

# Load phenotype data
pheno<-read.csv("you_phenotype_file.csv",stringsAsFactors=FALSE,header=TRUE, row.names = 1)


# Define key variables
ptsdVariable<-"Overall_PCL"    # Continuous PTSD variable 
covariates<-c("Age","Sex","CD8T","CD4T","NK","Bcell","Mono","Comp.1","Comp.2") # covariates
idVariable<-"Participant_ID"  # Random effect variable
studyID<-"RUSH"               # Study cohort ID

# Recode Sex: Male = 0, Female = 1
# check capitalization of male before you run (i.e is it MALE or male?)
pheno$Sex <-ifelse(pheno$Sex == "Male", 0, 1)
pheno$Sex <-as.numeric(as.character(pheno$Sex))


# Synchronize beta and phenotype data
beta.norm<-beta.norm[,names(beta.norm)%in%row.names(pheno)]
pheno<-pheno[row.names(pheno)%in%names(beta.norm),]
beta.norm<-beta.norm[,order(names(beta.norm))]
pheno<-pheno[order(row.names(pheno)),]
table(colnames(beta.norm)==rownames(pheno))   # Check alignment

range(pheno[,ptsdVariable])  # Check PTSD variable range

# Remove samples with missing data in key variables
naindex<-pheno[,c(ptsdVariable,covariates,idVariable)]
naindex<-complete.cases(naindex) 
beta.norm<-beta.norm[,naindex]
pheno<-pheno[naindex,]
table(rownames(pheno)==colnames(beta.norm))  # Confirm match

# Check beta value range (should be 0 < beta < 1)
range(beta.norm, na.rm=T)

# Replace problematic beta values
if(min(beta.norm, na.rm=T)==0){
  beta.norm[which(beta.norm==0)]<-0.0001
}
if(max(beta.norm, na.rm=T)==1.000000e+00){
  beta.norm[which(beta.norm==1.000000e+00)]<-0.9999
}

range(beta.norm, na.rm=T)
sum(is.na(beta.norm))  # Check number of NA values

# Convert beta values to M-values using logit transform
beta.norm<-log2(beta.norm/(1-beta.norm))
sum(beta.norm=="-Inf", na.rm=T)     # Check for log(0)
sum(is.na(beta.norm))               # Should match earlier NA count
range(beta.norm, na.rm=T)

# Prepare phenotype dataframe
pheno$meth<-NA
pheno[,idVariable]<-factor(pheno[,idVariable])  # Random effect as factor

# Preallocate results matrices for speed
resultsBeta<-matrix(nrow=nrow(beta.norm), ncol=2+length(covariates))
rownames(resultsBeta)<-rownames(beta.norm)
colnames(resultsBeta)<-c("(Intercept)", ptsdVariable, covariates)

resultsSE<-resultsT<-resultsP<-resultsDF<-resultsBeta  # Duplicate structure
errorProbes<-NULL
warningProbes<-NULL

# Construct formula dynamically
formula<-as.formula(paste("meth~",ptsdVariable, "+", 
                          paste(covariates, collapse="+"), 
                          "+(1|", idVariable, ")", sep=""))

print(formula)  # Confirm model formula

# Running analyses # updated March 13 2025

start<-proc.time()[3]  # Start timer

for(ii in 1:nrow(beta.norm)){
  
  # Assign CpG values to phenotype
  pheno$meth<-t(beta.norm[ii, rownames(pheno)])
  
  # Fit the mixed model
  fit<-try(lmer(formula, data=pheno), silent=F)
  
  if(class(fit)!="try-error"){
    
    res<-coef(summary(fit))
    
    resultsBeta[ii,]<-res[, "Estimate"]
    resultsSE[ii,]<-res[, "Std. Error"]
    resultsT[ii,]<-res[, "t value"]
    resultsP[ii,]<-res[, "Pr(>|t|)"]
    resultsDF[ii,]<-res[, "df"]
    
  }
  
  # Log any failed probes
  if(class(fit)=="try-error"){
    errorProbes<-append(errorProbes, rownames(beta.norm)[ii])
  }
  
  # Print progress every 10 CpGs
  if(ii%%10==0){print(ii)}
  
  # Reset meth column
  pheno$meth<-NA
}

end<-proc.time()[3]
end-start  # Total time

# Final output formatting
final<-as.data.frame(cbind(resultsBeta[,ptsdVariable],resultsSE[,ptsdVariable],
                           resultsT[,ptsdVariable],resultsP[,ptsdVariable],
                           resultsDF[,ptsdVariable]))
rownames(final)<-rownames(beta.norm)
colnames(final)<-c("BETA","SE","t","pval","df")

# Save results
write.csv(final,file = paste0(studyID,"_PCL_covar_AgeSexCellCompPCs_Mar13_25.csv"),quote = F,row.names = T)
