#!/usr/bin/env Rscript --vanilla

################################################################################
## Author: Dasia apdated code from Seyma's gitub https://github.com/PGC-PTSD-EWAS/EPICv2_pipeline
## Purpose: EPICv2 QC pipeline on part of the EVP

## Usage: Rscript --vanilla enmix_pipeline.R /my/input/path/ /my/output/path/ /my/otherfile/path/ my_project_name /my/sample/sheet
## See README file for usage, input, and output information.

################################################################################
## - Clean environ start.
rm(list = ls())
gc()
################################################################################
## If running line by line in Rstudio define the paths here.
setwd("/Volumes/som-ts/GynOb/SmithLab/NewData/People/Dasia/RUSH/")

INPUTPATH <- "/Volumes/som-ts/GynOb/SmithLab/NewData/Projects/PTSD/RawFiles/PTSDtreatment/"
OUTPUTPATH <- "/Volumes/som-ts/GynOb/SmithLab/NewData/People/Dasia/RUSH/Step1_Enmix_QC/Output/"
OTHERPATH <- "/Volumes/som-ts/GynOb/SmithLab/NewData/People/Dasia/RUSH/Step1_Enmix_QC/Other_path/"
PROJECTNAME <- "RUSH_QC_DYS_Nov142024"
SampleSheet <- "/Volumes/som-ts/GynOb/SmithLab/NewData/People/Dasia/RUSH/Step1_Enmix_QC/Input_files/RUSH_SampleSheet_DYS.csv"

#' See the example sample sheet (Samplesheet_Sample.csv). Must have columns: 
#' SampleID - Methylation ID as SentrixID_SentrixPosition
#' Sample_Name - Study specific sample name
#' CONTROL - Technical Control Samples, TRUE for technical controls, FALSE for real samples
#' Sex - MALE or FEMALE
#' Age
#' Pheno - Phenotype of Interest, 1 for cases and 0 for controls
################################################################################
## - Argument inputs
### - Setting global variables with args.
args = commandArgs(trailingOnly=T)
if (length(args)!=5) { stop("ERROR: Five args required. Exiting")
}else if (length(args)==5){
  INPUTPATH=args[1]
  OUTPUTPATH=args[2]
  OTHERPATH=args[3]
  PROJECTNAME=args[4]
  SampleSheet=args[5]
}
print("Input arguments:")
print(paste0("arg 1 INPUTPATH: ", INPUTPATH))
print(paste0("arg 2 OUTPUTPATH: ", OUTPUTPATH))
print(paste0("arg 3 OTHERPATH: ", OTHERPATH))
print(paste0("arg 4 PROJECTNAME: ", PROJECTNAME))
print(paste0("arg 5 SampleSheet: ", SampleSheet))
cat("\n")
################################################################################
# Saving the output
sink(paste0(OUTPUTPATH, PROJECTNAME, "_",Sys.Date(), "_Output.txt"), split = TRUE)
################################################################################
## - Load in libraries.
library(ENmix)
library(minfi)
library(ewastools)
library(impute)

library(tidyverse)

# Load EPICv2 annotation
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# Load in replicate functions.
source(paste0(OTHERPATH,"filter_functions_replicates.R"))

################################################################################
## - Functions.
file_list <- function(path){
  files <- list.files(path, recursive = TRUE, pattern = ".idat")
  files <- gsub(paste(c("_Grn.idat", "_Red.idat"), collapse = "|"), "", files)
  file_paths <- paste0(path, files)
  files <- gsub(".*/","",files)
  return(list(names = unique(files), Basename = unique(file_paths)))
}


changename = function(filename, newsuffix){
  if(file.exists(filename)){
    new = sub(".jpg$", "", filename)
    file.rename(filename, paste0(new,newsuffix))
  }else{
    print("FILE NOT FOUND:")
    print(filename)
  }
}


## Based on the ENmix sex prediction function, but set up to not require an rgDataSet object.
## When using Pidley's annotation, cutoff should be lower, try 0.5
sexprediction = function(mdat, PIDSLEYCSV, cutoff=2.0){
  
  chrXprobenames = PIDSLEYCSV[PIDSLEYCSV$CpG_chrm %in% c("chrX", "X"),]
  chrYprobenames = PIDSLEYCSV[PIDSLEYCSV$CpG_chrm %in% c("chrY", "Y"),]
  methMdat1 = assays(mdat)$Meth
  unmethMdat1 = assays(mdat)$Unmeth
  chrXrowsMeth = methMdat1[rownames(methMdat1) %in% chrXprobenames$IlmnID,]
  chrXrowsUnmeth = unmethMdat1[rownames(unmethMdat1) %in% chrXprobenames$IlmnID,]
  chrYrowsMeth = methMdat1[rownames(methMdat1) %in% chrYprobenames$IlmnID,]
  chrYrowsUnmeth = unmethMdat1[rownames(unmethMdat1) %in% chrYprobenames$IlmnID,]
  xCN = log2(chrXrowsMeth + chrXrowsUnmeth)
  yCN = log2(chrYrowsMeth + chrYrowsUnmeth)
  xMed = apply(xCN, 2, median, na.rm = TRUE)
  yMed = apply(yCN, 2, median, na.rm = TRUE)
  diff = xMed - yMed
  sex <- ifelse(diff > cutoff, "FEMALE", "MALE")
  return(sex)
  
}

## Negate 
`%ni%` <- Negate(`%in%`)

################################################################################
## - File reading.
print("Reading in phenotype CSV file, idat files.")

## - Read in:
pheno = read.csv(SampleSheet, 
                 header=T, 
                 sep=",", 
                 stringsAsFactors=F, 
                 fill=F)

print("Head and dimensions of phenotype file:")
head(pheno)
dim(pheno)

print("Printing file_paths:")
file_paths <- file_list(path = INPUTPATH)
file_paths <- as.data.frame(file_paths)
head(file_paths)

gc() # free up memory 

## Select idats just for our data
print("Printing paths and basenames:")
paths <- file_paths[which(file_paths$names %in% pheno$SampleID),]
head(paths)
stopifnot(all(grepl(paste(c(pheno$SampleID), collapse = "|"), paths$Basename))) # stop if all samples are not in your path
head(paths$Basename)
cat("\n")


## Reorder the columns make SampleID the first column in "pheno" dataframe so the next line of code can work DYS
#pheno <- pheno[, c("SampleID", setdiff(names(pheno), "SampleID"))]


## Merge with Pheno
pheno <- merge(pheno, paths, by = 1) 
head(pheno)

## Load idats
print("Reading in idats with minfi.")
RGset <- read.metharray.exp(targets = pheno, verbose = T, extended = T)

## Add EPIC v2 annotation
print("Idats read in:")
RGset@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38") 
RGset

print("Preparing input for ewastools.")
meth <- read_idats(paths$Basename, quiet = F)

metadata(RGset) <- as(meth, "list")
print("Input complete.")
cat("\n")

gc() # free up memory 

print("Saving RGset object.")
save(RGset, file = paste0(OUTPUTPATH, PROJECTNAME, "_RGset.RData"))


# load(paste0(OUTPUTPATH, PROJECTNAME, "_RGset.RData"))
print("Save complete.")
cat("\n")
################################################################################
## - Control plot generation
print("Generating control plots.")
plotCtrl(RGset)
jpglist = list.files(pattern = "\\.jpg$")
x = paste0(getwd(), "/")

if (x != OUTPUTPATH){
  print("Moving control plots from current working directory to command line argument output path.")
  sapply(jpglist, function(x) file.copy(from=x, to=OUTPUTPATH, copy.mode = T))
  sapply(jpglist, function(x) file.remove(from=x))
}
cat("\n")
################################################################################
## - ewastools control metric failures

print("Processing control metrics.")
ctrls = control_metrics(metadata(RGset))
failed <- sample_failure(ctrls)
names(failed) <- RGset@metadata[["meta"]][["sample_id"]]
failed <- as.data.frame(failed)
pheno <- merge(pheno, failed, by.x = "SampleID", by.y = "row.names")
fails = subset(pheno, failed == "TRUE")
print(paste0("Failed Samples: ", nrow(fails)))
fails[,c(1,2)]

write.csv(fails, "ewastools_control_metric_failures_fails.csv")

print("Control metric processing with ewastools finished.")
cat("\n")

## - EWAStools duplicate check
print("Checking duplicates")

snps <- getSnpBeta(RGset)
genotypes <- call_genotypes(snps, learn = FALSE)
donor_id <- enumerate_sample_donors(genotypes)
names(donor_id) <- RGset@metadata[["meta"]][["sample_id"]]
donor_id <- as.data.frame(donor_id)
pheno <- merge(pheno, donor_id, by.x = "SampleID", by.y = "row.names")

# List duplicates
pheno <- pheno %>%
  group_by(donor_id) %>%
  mutate(n = n()) %>%
  ungroup()

pheno[pheno$n > 1, c(1:2)]

print("EWAStools duplicate check finished.")
cat("\n")

################################################################################
## - ENmix QCinfo

## Notes: 
# 1. pvalue threshold at 0.05.
# 2. Missing/low quality data point threshold set to 0.05 
print("Running QCinfo of ENmix.")
setwd(OUTPUTPATH)
qc1 <- QCinfo(RGset, 
              detPtype = "negative",    # based on negtive internal control probes
              detPthre = 0.05,      # Set to 0.05.
              nbthre = 3,               # Number of bead threshold 
              samplethre = 0.05,         # % of low quality methylation data points across probes for each sample
              CpGthre = 0.05,            # % of low quality methylation data points across samples for each probe
              bisulthre = NULL,         # Threshold of bisulfite intensity
              outlier = TRUE, 
              distplot = TRUE)

save(qc1, file = "qc1.RData")
#load("qc1.RData")
bisulfthreshold = mean(qc1$bisul - (3*sd(qc1$bisul)))


changename("freqpolygon_beta_beforeQC.jpg", "_detP_beforeQC_ENmix.jpg")
changename("freqpolygon_beta_afterQC.jpg", "_detP_afterQC_ENmix.jpg")

################################################################################
## - ENmix primary data processing.

### Data preprocessing
# Background correction and dye bias correction
# Excluding low quality samples and probes

setwd(OUTPUTPATH)
print("Processing data with primary ENmix pipeline.")
mdat1 <- preprocessENmix(RGset, 
                         bgParaEst = "oob", 
                         dyeCorr = "RELIC", 
                         exQCsample = T, 
                         exQCcpg = T,
                         QCinfo = qc1, 
                         nCores = 8) # was 12 

## takes up a lot of memory need to clear everything out then continue code  with mdat1 object  #DYS
gc()             # This triggers garbage collection


#START RIGHT HERE #DYS

# Save the mdat1 object as an R data file #DYS
#save(mdat1, file = "mdat1.RData")

# Load the mdat1 object from the saved file #DYS
#load("mdat1.RData")


#########################################
#########################################


# Between-array normalization
qmdat1 <- norm.quantile(mdat1, method="quantile1")

#save(qmdat1, file = "qmdat1.RData") # saving R data file  just in case was computer runs out of memory DYS

# have to run then laod in different part becuase this computer (runs out of memory) DYS
#load("qmdat1.RData")

#rm(meth,RGset) # DYS # remove big ite
#gc() # DYS
#########################################
#########################################
#load(qmdat1.Rdata) # DYS  
#load("qc1.Rdata") #DYS

# Probe-type bias adjustment
betaMat1 <- rcp(qmdat1, qcscore = qc1)

# Filter low quality and outlier data points for each probe
# Remove rows and columns with too many missing values (>5%) will be removed if specify
beta1 <- qcfilter(betaMat1, qcscore = qc1, rmoutlier = F,
                  detPthre = 0.05, nbthre = 3,
                  rmcr = TRUE, rthre = 0.05, cthre = 0.05,
                  impute = F)

#save(beta1, file = "beta1.RData") # DYS
#rm(qmdat1)
#gc()

print("Saving RData file.")

# save(qc1, mdat1, qmdat1, betaMat1, beta1, file = paste0(OUTPUTPATH,PROJECTNAME,"_detP_0.05.RData")
save(qc1, mdat1, qmdat1, betaMat1, beta1, file = paste0(OUTPUTPATH,PROJECTNAME,"_preprocessENmix_Output.RData"))

################################################################################
print("Obtaining beta distribution plots.")
setwd(OUTPUTPATH)
jpeg(paste0(OUTPUTPATH,PROJECTNAME,"_BetaDist_ENmix.jpg"), height = 1600, width = 800)
par(mfrow=c(5,1))

#load("RUSH_QC_DYS_Nov142024_RGset.RData") #DYS

## Pre-QC plot
mraw <- getmeth(RGset)
beta <- getB(mraw)


png("Before QC.png",, width = 1200, height = 600)
multifreqpoly(beta, main = "Before QC", xlab = "Beta value", legend = F)
dev.off()
rm("RGset") #DYS
gc()
load("RUSH_QC_DYS_Nov142024_preprocessENmix_Output.RData")




#### START HERE 
#load("mdat1.RData")
## After preprocess
png("After_Preprocessing.png", width = 1200, height = 600)
multifreqpoly(getB(mdat1), main = "After Preprocessing", xlab = "Beta value", legend = F)
dev.off()
#rm(mdat1)


#load("qmdat1.RData")
## After between-array normalization
png("After_Between_Array_Normalization.png", width = 1200, height = 600)
multifreqpoly(getB(qmdat1), main = "After Between-Array Normalization", xlab = "Beta value", legend = F)
dev.off()
#rm(qmdat1)

#gc()
## After probe-type bias adjustment
png("After_Probe_Type_Bias_Adjustment.png", width = 1200, height = 600)
multifreqpoly(betaMat1, main = "After Probe-type Bias Adjustment", xlab = "Beta value", legend = F)
dev.off()


## After qc filter
#load("beta1.RData")
png("After_QC.png", width = 1200, height = 600)
multifreqpoly(beta1, main = "After QC", xlab = "Beta value", legend = F)
dev.off()

cat("\n")

################################################################################
## Outputting QCinfo information.
## Find the outlier sample
print("QCinfo processing finished. Outputting failed samples.")
out <- qc1$badsample
pheno[pheno$SampleID%in%out,c(1,2)]


bisulfite_fails = names(which(qc1$bisul < bisulfthreshold))
outlier_samples = qc1$outlier_sample

print("Samples failed from bisulfite failures:")
bifailures = pheno[pheno$SampleID %in% bisulfite_fails,][,1:2]
bifailures

write_csv(bifailures, file = "Samples_failed_from_bisulfite_failures.csv") #DYS

cat("\n")
print("Outlier samples:")
outlierfailures = pheno[pheno$SampleID %in% outlier_samples,][,1:2]
outlierfailures

write_csv(outlierfailures, file= "Outlier_samples.csv") #DYS

cat("\n")
print("QCinfo 'badsamples':")
enmixspecific_badsamples = pheno[pheno$SampleID %in% qc1$badsample,][,1:2]
enmixspecific_badsamples
write.csv(enmixspecific_badsamples, file = "Badsamples.csv") #DYS
cat("\n")

################################################################################
## Constructing flagged sample dataframe with pheno.csv input.
all_samples_for_flagged = pheno

all_samples_for_flagged = all_samples_for_flagged %>%
  mutate(ewastools_controls =
           case_when(all_samples_for_flagged$SampleID %in% fails$SampleID ~ "FAIL",
                     all_samples_for_flagged$SampleID %ni% fails$SampleID ~ "PASS"))

all_samples_for_flagged = all_samples_for_flagged %>%
  mutate(ENmix_bisulfite_flags =
           case_when(all_samples_for_flagged$SampleID %in% bifailures$SampleID ~ "FAIL",
                     all_samples_for_flagged$SampleID %ni% bifailures$SampleID ~ "PASS"))

all_samples_for_flagged = all_samples_for_flagged %>%
  mutate(ENmix_outlier_flags =
           case_when(all_samples_for_flagged$SampleID %in% outlierfailures$SampleID ~ "FAIL",
                     all_samples_for_flagged$SampleID %ni% outlierfailures$SampleID ~ "PASS"))

all_samples_for_flagged = all_samples_for_flagged %>%
  mutate(ENmix_flags =
           case_when(all_samples_for_flagged$SampleID %in% enmixspecific_badsamples$SampleID ~ "FAIL",
                     all_samples_for_flagged$SampleID %ni% enmixspecific_badsamples$SampleID ~ "PASS"))


control_samples = all_samples_for_flagged[all_samples_for_flagged$CONTROL == "TRUE",] # HERE

control_samples = control_samples %>%
  mutate(control_flagged_by_checks =
           case_when((control_samples$ewastools_controls == "FAIL") ~ TRUE,
                     (control_samples$ENmix_bisulfite_flags == "FAIL") ~ TRUE,
                     (control_samples$ENmix_outlier_flags == "FAIL") ~ TRUE,
                     (control_samples$ENmix_flags == "FAIL") ~ TRUE,
                     (control_samples$ewastools_controls == "PASS") ~ FALSE,
                     (control_samples$ENmix_bisulfite_flags == "PASS") ~ FALSE,
                     (control_samples$ENmix_outlier_flags == "PASS") ~ FALSE,
                     (control_samples$ENmix_flags == "PASS") ~ FALSE))
true_control_samples = control_samples[control_samples$control_flagged_by_checks == TRUE,]
false_control_samples = control_samples[control_samples$control_flagged_by_checks == FALSE,]
all_samples_for_flagged = all_samples_for_flagged %>%
  mutate(control_flagged =
           case_when(all_samples_for_flagged$SampleID %in% true_control_samples$SampleID ~ TRUE,
                     all_samples_for_flagged$SampleID %in% false_control_samples$SampleID ~ FALSE,
                     all_samples_for_flagged$SampleID %ni% control_samples$SampleID ~ NA))

################################################################################
## Cleaning to save memory
rm(qmdat1, mdat1)
gc()
################################################################################
## - Testing male female flags.
print("Sex prediction.")

epicv2manifestfilename = paste0(OTHERPATH, "pidsley2024.csv")
PIDSLEYCSV = data.table::fread(epicv2manifestfilename, header=T,
                               stringsAsFactors = F, sep=",", fill=T)


gc()
sexprediction_vector = sexprediction(mraw, PIDSLEYCSV, cutoff=0.5)
sexprediction_table = data.frame(sexprediction_vector)
sexprediction_table = tibble::rownames_to_column(sexprediction_table, "SampleID")
colnames(sexprediction_table) = c("SampleID", "ENmix_Predicted_Sex")

rm(mraw)
gc()
cat("\n")


all_samples_for_flagged = merge(all_samples_for_flagged, sexprediction_table, by="SampleID")
all_samples_for_flagged = all_samples_for_flagged %>%
  mutate(does_sex_match =
           case_when(toupper(all_samples_for_flagged$Sex) == all_samples_for_flagged$ENmix_Predicted_Sex ~ TRUE,
                     toupper(all_samples_for_flagged$Sex) != all_samples_for_flagged$ENmix_Predicted_Sex ~ FALSE))
################################################################################
print(paste0("Starting Samples = ", dim(RGset)[2])) #482
print(paste0("Starting Probes = ", dim(RGset)[1])) #1105209
gc()
#load("beta1.Rdata")
print(paste0("Remaining Samples = ", dim(beta1)[2])) #444
print(paste0("Remaining Probes = ", dim(beta1)[1])) #931782
print(paste0("Missing data = ", sum(is.na(beta1)))) #800143
cat("\n")
################################################################################
################################################################################
## - Filtering Illumina flagged, inaccurately mapped probes and cross reactive probes.

BETAFILE = beta1
print("Dimensions of beta file before removing chr0 probes, inaccurate probes, offtarget hitting probes, flagged probes, and managing replicates:")
dim(BETAFILE) # 931782    444  

######
# Pidlsey csv file:
print("EPICv2 manifest dimensions:")
dim(PIDSLEYCSV)
cat("\n")

# chr0 probe removal using Pidsley, no chr but remove NA
print("Obtaining chr0/unmapped probes from manifest and removing them from data:")
chr0probes = PIDSLEYCSV[is.na(PIDSLEYCSV$CpG_chrm),]$IlmnID
print(paste0("Number of chr0 probes: ",length(chr0probes)))

print(paste0("Number of probes currently in beta file: ", dim(BETAFILE)[1]))
BETAFILE = BETAFILE[rownames(BETAFILE) %ni% chr0probes,]
print(paste0("Number of probes in beta file following chr0 probe removal: ", dim(BETAFILE)[1])) #932414
cat("\n")

# inaccurately mapped probe removal using Illumina official csv file:
print("Reading in mapping inaccuracy file from Illumina - EPIC-8v2-0_A1-190MappingInaccuracies.csv.")
mapping_inaccuracy_filename = paste0(OTHERPATH,
                                     "EPIC-8v2-0_A1-190MappingInaccuracies.csv")
epicv2_mapping_inacc = read.table(mapping_inaccuracy_filename, sep=",", stringsAsFactors = F, header=T)
print('Dimensions of file:')
dim(epicv2_mapping_inacc) # 190
inaccurate_probes_from_illumina = epicv2_mapping_inacc$IlmnID

print("Number of inaccurate probes listed in file:")
length(inaccurate_probes_from_illumina) # 190
cat('\n')

print(paste0("Dimensions of beta file before inaccurate probe processing: ", dim(BETAFILE)[1]))
BETAFILE = BETAFILE[rownames(BETAFILE) %ni% inaccurate_probes_from_illumina,]
print(paste0("Dimensions of beta file post inaccurate probe processing: ", dim(BETAFILE)[1])) 
cat("\n")


# flagged probes using Illumina official csv file:
print("Reading in flagged probe file from Illlumina - EPIC-8v2-0_A1-FlaggedProbes.csv:")
flagged_probe_filename = paste0(OTHERPATH, 
                                "EPIC-8v2-0_A1-FlaggedProbes.csv")
epicv2_flagged_probes = read.table(flagged_probe_filename, sep=",", stringsAsFactors = F, header = T)
flagged_probes_from_illumina = epicv2_flagged_probes$IlmnID
print("Number of flagged probes listed in file:")
length(flagged_probes_from_illumina) # 50209

print(paste0("Dimensions of beta file before flagged probe processing: ", dim(BETAFILE)[1]))
BETAFILE = BETAFILE[rownames(BETAFILE) %ni% flagged_probes_from_illumina,]
print(paste0("Dimensions of beta file post flagged probe processing: ", dim(BETAFILE)[1])) #
cat("\n")

# offtarget probes using PIDSLEYCSV:
print("Removing offtargets using Pidsley 2024 csv file:")
print(paste0("Dimensions of beta file before offtarget processing: ",dim(BETAFILE)[1]))
offtargets = PIDSLEYCSV[PIDSLEYCSV$Num_offtargets != 0,]$IlmnID
BETAFILE = BETAFILE[rownames(BETAFILE) %ni% offtargets,]
print(paste0("Dimensions of beta file post offtarget processing: ", dim(BETAFILE)[1])) 
cat("\n")

# Non-mapping probes from PIDSLEY (Table S15):
print("Reading in Probes with no BLAT hits (Table S15):")
nonmapping_probe_filename = paste0(OTHERPATH,
                                   "NoMapping_Pidsley.txt")
nonmapping_probes = read.table(nonmapping_probe_filename, sep=" ", stringsAsFactors = F, header = F)
nonmapping_probes = nonmapping_probes$V1
print("Number of nonmapping probes listed in file:")
length(nonmapping_probes) # 24

print(paste0("Dimensions of beta file before nonmapping probe processing: ", dim(BETAFILE)[1]))
BETAFILE = BETAFILE[rownames(BETAFILE) %ni% nonmapping_probes,]
print(paste0("Dimensions of beta file post nonmapping probe processing: ", dim(BETAFILE)[1])) # 857773
cat("\n")

# Managing replicates based on PIDSLEYCSV file:
print("Processing replicates using Pidsley 2024 csv file:")
print(paste0("Dimensions of beta file before replicate processing: ", dim(BETAFILE)[1]))

PIDSLEYCSV = paste0(OTHERPATH, "pidsley2024.csv")
FLAGGEDCSV = flagged_probe_filename

BETAFILE = process_replicates(BETAFILE, PIDSLEYCSV, FLAGGEDCSV)
print(paste0("Dimensions of beta file post replicate processing: ", dim(BETAFILE)[1])) # 846179
cat("\n")

################################################################################)
################################################################################
## - Returning a beta matrix with suffixes and a beta matrix without suffixes.

BETAFILE = as.data.frame(BETAFILE)
print("Creating a beta matrix without epicv2 suffixes.")
BETAFILE_suffix_adj = tibble::rownames_to_column(BETAFILE, "probenames")
BETAFILE_suffix_adj$probenames = gsub("(.*)_.*", "\\1", BETAFILE_suffix_adj$probenames)

print("Are all probes in the beta matrix without suffixes unique?")
length(unique(sort(BETAFILE_suffix_adj$probenames))) == length(BETAFILE_suffix_adj$probenames)
cat("\n")
BETAFILE_suffix_adj = BETAFILE_suffix_adj %>% tibble::column_to_rownames(var = "probenames")
print("3x3 of suffix-altered beta matrix:")
BETAFILE_suffix_adj[1:3,1:3]


cat("\n")
print("Saving two beta value data frames - one with epicv2 suffixes and one without the suffixes.")
save(BETAFILE, file=paste0(OUTPUTPATH, PROJECTNAME, "_betas_with_suffix.RData"))
save(BETAFILE_suffix_adj, file=paste0(OUTPUTPATH, PROJECTNAME, "_betas_without_suffix.RData"))
print("Save complete.")
cat("\n")

print("Saving the FINAL beta value data frames without the suffixes and without failed samples and controls.")
all_samples_pass <-subset(all_samples_for_flagged, failed == FALSE & CONTROL == FALSE & ENmix_flags == "PASS")
BETAFILE_PASS <- BETAFILE_suffix_adj[,colnames(BETAFILE_suffix_adj) %in% all_samples_pass$SampleID]
print(paste0("Number of samples remaining after removing controls and failed samples: ", dim(all_samples_pass)[1])) #186
save(BETAFILE_PASS, file=paste0(OUTPUTPATH, PROJECTNAME, "_betas_without_suffix_passQC.RData"))

#write.csv(all_samples_pass,"SAMPLES_thatPass_QC_oct24.csv") #DYS


print("Saving pheno.csv updated with flags for ewastools control failures, ENmix internal check failures, and control issues, and sex prediction.")
write.csv(all_samples_for_flagged, file=paste0(OUTPUTPATH, PROJECTNAME, "_pheno_and_flagged_ENmix.csv"), row.names = F)
print("Save complete.")
cat("\n")

print("Saving pheno.csv that contains only samples pass QC")
write.csv(all_samples_pass, file=paste0(OUTPUTPATH, PROJECTNAME, "_pheno_ENmix_passQC.csv"), row.names = F)
print("Save complete.")
cat("\n")

################################################################################

print(paste0("Remaining Samples = ", dim(BETAFILE)[2])) 
print(paste0("Remaining Probes = ", dim(BETAFILE)[1])) 
print(paste0("Missing data = ", sum(is.na(BETAFILE))))

################################################################################
sessionInfo()

print("Script complete. Exiting.")

################################################################################
# Close sink
sink()
