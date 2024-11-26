############ 3) ComBAT normalization to adjust for batch effects (chip and array position) ############
########################################################################################################
#' Clean
rm(list=ls())
gc()

library(impute)



# Load Noob Normalized Cross-Reactive Probes removed beta values
beta <- fread("RUSH_QC_beta_without_suffix_nov20_2024.csv", data.table = F) # cpgs x samples
row.names(beta)<-beta$V1
beta<-beta[,-1]

# Load the new phenotype file with QC info
phen <- read.csv("RUSH_Pheno_QC_withCell_types_PCLoverall_Nov20_24.csv", stringsAsFactors = FALSE, header = TRUE)
row.names(phen)<-phen$SampleID
#phen$X <- NULL # remove a column 

# Rename column
names(phen)[names(phen) == "Overall_PCL"] <- "PCL"
phen <- phen[, !(names(phen) %in% c("X"))]
# Remove the failed samples and sex mismatches
#phen <- subset(phen, failed == "FALSE" & SexMismatch == "FALSE")

# Methylation IDs (Row names of phenotype file and Column names of beta matrix) should match
beta <- beta[, which(names(beta) %in% row.names(phen))]
phen <- phen[row.names(phen) %in% names(beta), ]
beta <- beta[, order(names(beta))]
phen <- phen[order(row.names(phen)), ]
stopifnot(all(colnames(beta) == rownames(phen)))

## Define Variables for Combat step
## Define model.matrix, which includes the variables that you want to protect when adjusting for chip and position
## Generally variables that we use as covariates in the EWAS (sex, age, main phenotype -PTSD-, smoking) are included in the model.matrix
sex <- "Sex"         # Name of Sex Variable: Males coded as 0, females coded as 1
age <- "Age"         # Name of Age Variable
PTSD <- "PCL"       # Name of PTSD Variable: Cases coded as 1, controls coded as 0
chip <- "Sentrix_ID"
position <- "Sentrix_Position"

## You should not have NAs in model matrix, so we remove subjects with no phenotype info
print(paste0("Samples with no PTSD information = ", sum(is.na(phen[,PTSD]))))
print(paste0("Samples with no Sex information = ", sum(is.na(phen[,sex]))))
print(paste0("Samples with no Age information = ", sum(is.na(phen[,age]))))

naindex <- (!is.na(phen[,PTSD]) & !is.na(phen[, age]) & !is.na(phen[, sex]))
phen <- phen[naindex, ]
beta <- beta[, naindex]
stopifnot(all(colnames(beta) == rownames(phen)))

chip <- as.factor(phen[,chip])
position <- as.factor(phen[,position])
ptsd<-as.numeric(phen[,PTSD])
age<-as.numeric(phen[,age])
sex<-as.factor(phen[,sex])

moddata <- model.matrix(~ptsd+age+sex)

# Remaining samples
print(paste0("Remaining Samples = ", nrow(phen))) # N = 444

## ComBAT does not handle NAs in the methylation file,
## so we have to impute NAs in the methylation beta matrix
## Log transform to normalize the data (we'll reverse that later)
beta <- log((beta/(1-beta)))
beta.imputed <- impute.knn(as.matrix(beta))
beta.imputed <- beta.imputed$data

rm(beta)
gc()

# Run ComBat
SerialParam()
combat_beta <- ComBat(dat = beta.imputed, mod = moddata, batch = chip, BPPARAM = SerialParam())
combat_beta <- ComBat(dat = combat_beta, mod = moddata, batch = position, BPPARAM = SerialParam())

# Reverse Beta Values
reversbeta <- 1/(1+(1/exp(combat_beta)))

## We need to put NAs back (from the original matrix) to ComBAT adjusted beta matrix
## We don't want to use imputed beta values for missing data
norm <- fread("RUSH_QC_beta_without_suffix_nov20_2024.csv", data.table = F)
rownames(norm) <- norm$V1
norm<-norm[,-1]

beta <- as.data.frame(reversbeta)

norm <- norm[, names(norm) %in% names(beta)]
norm <- norm[, order(names(norm))]
beta <- beta[, order(names(beta))]
table(names(norm) == names(beta))

beta[is.na(norm)] <- NA
beta[beta >= 0.9999998] <- 0.9999999

write.csv(beta,file="RUSH_combat_CP_wcovar_age_pcl_sex_Nov20_24.csv",quote=FALSE,row.names=TRUE)
# This is the final beta values that you'll use in EWAS
