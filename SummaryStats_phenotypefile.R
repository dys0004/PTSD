# Load necessary libraries
library(dplyr) # for data formating 
library(knitr) # for pretty tables 
library(lme4)       # linear mixed-effects models
library(lmerTest)   # test for linear mixed-effects models
library(kableExtra)
library(ggplot2)
library(tidyverse) 

# Read in data 
datum <- read.csv("/Volumes/som-ts/GynOb/SmithLab/NewData/People/Dasia/RUSH/Summary_Statistics/RUSH_Pheno_QC_withCell_types_PCLoverall_Nov20_24.csv")


## Demographics of the dataset 
# Number and percentage of males and females
# group by participant ID
MvsF <-datum %>%
  group_by(Participant_ID) %>%
  summarize(Sex = first(Sex)) %>%
  count(Sex) %>%
  mutate(Percentage = (n / sum(n)) * 100)
MvsF

library(dplyr)

# Summarize Age by Sex, ensuring only one row per individual
AgeSummary <- datum %>%
  group_by(Participant_ID) %>%          # Group by individual
  summarise(
    Sex = first(Sex),                   # Take the first value of Sex for each individual
    Age = first(Age)                    # Take the first value of Age for each individual
  ) %>%
  group_by(Sex) %>%                     # Group by Sex
  summarise(
    Mean_Age = mean(Age, na.rm = TRUE), # Calculate mean age
    SD_Age = sd(Age, na.rm = TRUE),     # Calculate standard deviation
    Count = n()                         # Count number of individuals
  )

# Print the Age summary by Sex
print(AgeSummary)


## Racial demographics with percentages
Race <-datum %>%
  group_by(Participant_ID) %>%
  summarize(Race = first(Race)) %>%
  count(Race) %>%
  mutate(Percentage = (n / sum(n)) * 100)
Race

# Count Overall_PCL <= 10 for post samples directly
PCL_count_below_10 <- datum %>%
  summarise(
    Count = sum(grepl("_Post", Sample_Name) & Overall_PCL <= 10)
  )
print(PCL_count_below_10)

# Calculate mean and standard deviation for PCL_Baseline using only "Pre" samples
pcl_baseline_summary <- datum %>%
  summarise(
    Mean_PCL_Baseline = mean(PCL_Baseline[grepl("_Pre", Sample_Name)], na.rm = TRUE),
    SD_PCL_Baseline = sd(PCL_Baseline[grepl("_Pre", Sample_Name)], na.rm = TRUE)
  )

# Print the result
print(pcl_baseline_summary)


# Calculate mean and standard deviation for PCL_Completion using only "Post" samples
pcl_completion_summary <- datum %>%
  summarise(
    Mean_PCL_Completion = mean(PCL_Completion[grepl("_Post", Sample_Name)], na.rm = TRUE),
    SD_PCL_Completion = sd(PCL_Completion[grepl("_Post", Sample_Name)], na.rm = TRUE)
  )

# Print the result
print(pcl_completion_summary)



#Chi-square to test between relationships between categorical variables 
#To determine if there is an association between the (Responder/Nonresponder) and (Male/Female) columns using a chi square test

# Create a column within the data frame with PCL_overall if equal to or less then 10 list as a reponder 
datum <- datum %>%
  mutate(Response = ifelse(Overall_PCL <= 10, "Responder", "Nonresponder"))


# Count the number of Responders and Nonresponders
response_counts <- datum %>%
  count(Response)


# Filter for Post samples and count responders and nonresponders
response_counts_post <- datum %>%
  filter(grepl("_Post", Sample_Name)) %>%
  mutate(Response = ifelse(Overall_PCL <= 10, "Responder", "Nonresponder")) %>%
  count(Response)

# Create the bar plot
ggplot(response_counts_post, aes(x = Response, y = n, fill = Response)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of Responders vs Nonresponders (Post Samples)",
       x = "Response Category",
       y = "Count") +
  theme_minimal()


# Filter for Post samples for chi-square test assumes each observation (row) is independent
post_data <- datum %>%
  filter(grepl("_Post", Sample_Name))

# Create the contingency table for Response and Sex
contingency_table <- table(post_data$Response, post_data$Sex)

# Print the contingency table
print(contingency_table)

# Perform Chi-square test
chi_square_result <- chisq.test(contingency_table)


