
#############################
# Pheno, geno and DNAm QC   #
#############################

library(tidyverse)  
library(knitr)
library(stringr)
library(readr)

mycolours <- c("#00AFBB", "#E7B800", "#FC4E07", "#0072B2", 
               "#009E73", "#F0E442", "#D55E00", "#CC79A7") 


GSpath="/Cluster_Filespace/Marioni_Group/GS"
GSocsa="/Local_Data/methylation/GS_20k/OSCA_files"
scracth=""

#############################
# Data description          #
#############################

# -   18,349 individuals
# -   831,733 DNAm probes
# -   3 Waves
# 
# -   Wave 1 = 5,087
# -   Wave 3 = 4,450
# -   Wave 4 = 8,876

# Phenotypes

# -   Body Mass Index (kg/m\^2)
# -   Waist / Hip Ratio
# -   Body Fat Composition (bio-impedance) (%)
# -   Glucose (mmol/l)
# -   HDL cholesterol (mmol/l)
# -   Total cholesterol (mmol/l)
# -   Height (control)




#############################
# Quality Control
#############################



# Subset to samples with covariate information; n = 18377
data <- readRDS(paste0(GSpath, "/GS_methylation/GS20k/GS20k_Targets.rds"))
# wave1 wave3 wave4 
# 5087  4450  8876 

grm.id <- read.table(paste0(scratch,"/unrelated/GS_GWAS.grm.id", header=F) # Fam
colnames(grm.id) <- c("FID", "IID")

oii <- read.table(paste0(GSpath, "/mvals-norm20k-18413-831733.oii"), header=F) # # OSCA ID, V1==V2
colnames(oii) <- c("FID_DNAm", "IID_DNAm", "V3", "V4", "V5")


withdrawn <- read.table(paste0(GSpath,"/withdrawals_v3.txt"), header=F) # Withdrawn
colnames(withdrawn) <- c("FID", "IID")

# Phenotypes
# GSpath/GS_dataset/clinical/body.csv. 
# They are labelled as bmi, whr and body_fat there. 

# GSpath/GS_dataset/clinical/biochemistry.csv. 
# These are glucose, HDL cholesterol and Total cholesterol. 

body <- read.table(paste0(GSpath"/GS_dataset/clinical/body.csv"), header=T, sep=",")
biochemistry <- read.table(paste0(GSpath,"/GS_dataset/clinical/biochemistry.csv"), header=T, sep=",")

# Merge
ds <- oii %>% 
  select(FID_DNAm, IID_DNAm) %>%
  inner_join(data, by=c("FID_DNAm"="Sample_Sentrix_ID")) %>%
  inner_join(grm.id, by=c("Sample_Name"="IID")) %>%
  left_join(body, by=c("Sample_Name"="id")) %>%
  left_join(biochemistry, by=c("Sample_Name"="ID")) %>% 
  filter(!Sample_Name %in% withdrawn$IID) %>%
  select(FID_DNAm, IID_DNAm, FID, IID=Sample_Name, Sample_Name, everything())
# 18349 with DNAm 

write.table(ds, file=paste0(scratch,"/covariates/oii_cov.txt"), col.name=T, row.name=F, quote=F)



#############################
### Remove related individuals
#############################

mkdir /Local_Scratch/Alesha/unrelated
gcta64 \
--bfile /Cluster_Filespace/Marioni_Group/GS/GS_GWAS/GS_GWAS \
--remove /Local_Scratch/Alesha/unrelated/withdrawals_v3.txt \
--make-grm \
--out /Local_Scratch/Alesha/unrelated/GS_GWAS
#   19989 samples, 561125 markers, 199790055 GRM elements

gcta64 --grm /Local_Scratch/Alesha/unrelated/GS_GWAS \
--grm-cutoff 0.05 \
--make-grm \
--out GS_GWAS_0.05
# After pruning the GRM, there are 8041 individuals (11948 individuals removed).

# However we are going to select only participants with DNAm first to make sure we dont loose more people than necessary
# Check unrelated doesnt change numbers when filtering for DNAm ind first

awk 'NR>1 {print $3, $4}' /Local_Scratch/Alesha/covariates/oii_cov.txt > /Local_Scratch/Alesha/covariates/oii_cov_FID_IID.txt
# 0.05 cut off for those with DNAm first
gcta64 --grm /Local_Scratch/Alesha/unrelated/GS_GWAS \
--keep /Local_Scratch/Alesha/covariates/oii_cov_FID_IID.txt \
--grm-cutoff 0.05 \
--make-grm \
--out /Local_Scratch/Alesha/unrelated/GS_GWAS_DNAm_0.05
# After pruning the GRM, there are 7758 individuals (10618 individuals removed).
# Note: Resulted in more individuals this way by ~ 300.


# Check unrelated makes sense
# 1. Set by waves and determine unrelated
awk 'NR>1 &&  $7=="wave4" {print $3, $4}' oii_cov.txt > wave4_FID.txt
awk 'NR>1 &&  $7=="wave3" {print $3, $4}' oii_cov.txt > wave3_FID.txt
awk 'NR>1 &&  $7=="wave1" {print $3, $4}' oii_cov.txt > wave1_FID.txt

gcta64 --grm /Local_Scratch/Alesha/unrelated/GS_GWAS \
--keep /Local_Scratch/Alesha/covariates/wave1_FID.txt \
--grm-cutoff 0.05  --make-grm  --out GS_GWAS_DNAm_0.05_wave1

gcta64 --grm /Local_Scratch/Alesha/unrelated/GS_GWAS \
--keep /Local_Scratch/Alesha/covariates/wave3_FID.txt \
--grm-cutoff 0.05  --make-grm  --out GS_GWAS_DNAm_0.05_wave3

gcta64 --grm /Local_Scratch/Alesha/unrelated/GS_GWAS \
--keep /Local_Scratch/Alesha/covariates/wave4_FID.txt \
--grm-cutoff 0.05  --make-grm  --out GS_GWAS_DNAm_0.05_wave4

#  7758 GS_GWAS_DNAm_0.05.grm.id
#  2710 GS_GWAS_DNAm_0.05_wave1.grm.id
#  4444 GS_GWAS_DNAm_0.05_wave3.grm.id
#  4267 GS_GWAS_DNAm_0.05_wave4.grm.id

# However when run all together
# wave1 wave3 wave4 
# 2035  4276  1435 
```

7,758 unrelated individuals

-   wave1 - 2,035
-   wave3 - 4,276
-   wave4 - 1,435

```{r Remove related from covariates}
# Check how many are unrelated from each wave individually
# 
# -   wave1 - 2710
# -   wave3 - 4444
# -   wave4 - 4267

# Related with DNAm
# wave1 wave3 wave4 
# 5080  4445  8851 

# Looks like almost everyone in wave4 was related or related to individuals in previous waves

ds <- read.table("/Cluster_Filespace/Marioni_Group/Alesha/covariates/oii_cov.txt", header=T)
unrelated <- read.table("/Cluster_Filespace/Marioni_Group/Alesha/unrelated/GS_GWAS_DNAm_0.05.grm.id", header=F)
colnames(unrelated) <- c("FID", "IID")

ds_unrelated <- ds %>%
  inner_join(unrelated, by=c("FID", "IID")) 
# [1] 7758   22 # Unrelated with DNAm and covariates

write.table(ds_unrelated, file="/Cluster_Filespace/Marioni_Group/Alesha/covariates/oii_cov_unrelated.txt", col.name=T, row.name=F, quote=F)
```

<br>
  
  
  
#############################
### Phenotypes
#############################

# For each trait, we adjusted the phenotype for age in each gender group and standardized the residuals by rank-based inverse normal transformation, which removed the age effect and potential difference in mean and variance between two gender groups. There was no adjustment for wave as there were minimal differences in the mean across waves.

# -   BMI
# -   Waist-to-hip ratio
# -   Body fat %
# -   Biochemical traits
#     -   Glucose
#     -   HDL cholesterol
#     -   Total cholesterol
# Other traits available include
# 
# -   Sodium
# -   Potassium
# -   Urea
# -   Creatinine
# -   Creat_mgdl

library(RNOmni) 
library(tidyverse)
# Unrelated
ds <- read.table("/Cluster_Filespace/Marioni_Group/Alesha/covariates/oii_cov_unrelated.txt", header=T)
# Related
# ds <- read.table("/Cluster_Filespace/Marioni_Group/Alesha/covariates/oii_cov.txt", header=T)


summary(ds)
dim(ds)

# Look at the data
library(psych) 
describe(ds %>% select(-FID_DNAm, -IID_DNAm, -IID, -FID, -Set, -Batch, -sex)) %>% 
  select(n, mean, median, sd, min, max, range)

# Plot each var
ds %>% str

t_ds_numeric_sex <- reshape2::melt(ds %>% select(FID_DNAm, IID_DNAm, FID, IID, sex, where(is.numeric), where(is.integer)), 
                                   id = c('FID_DNAm','IID_DNAm', 'FID', 'IID', "sex"))

t_ds_numeric_Set <- reshape2::melt(ds %>% select(FID_DNAm, IID_DNAm, FID, IID, Set, where(is.numeric), where(is.integer)), 
                                   id = c('FID_DNAm','IID_DNAm', 'FID', 'IID', "Set"))

# By sex
t_ds_numeric_sex %>% 
  filter(variable %in% c("age", "bmi", "whr", "body_fat", "Glucose", "HDL_cholesterol", "Total_cholesterol", "height", "waist", "hips")) %>%  
  ggplot(aes(x = value, fill=sex)) +
  geom_density(alpha=0.5)  +
  facet_wrap( ~ variable, ncol = 4,  scales = "free")


# By Set
t_ds_numeric_Set  %>% 
  filter(variable %in% c("age", "bmi", "whr", "body_fat", "Glucose", "HDL_cholesterol", "Total_cholesterol", "height", "waist", "hips")) %>% 
  ggplot(aes(x = value, fill=Set)) +
  geom_density(alpha=0.5)  +
  facet_wrap( ~ variable, ncol = 4,  scales = "free")



# Mean differences by set?
ds_Set <-   ds %>% group_by(Set) %>% # or Batch?
  summarise(sex_M = sum(sex == "M")/sum(!is.na(sex)) * 100,
            sex_F = sum(sex=="F")/sum(!is.na(sex)) * 100) %>%
  left_join(
    ds %>% group_by(Set) %>% # or Batch?
      select(-FID_DNAm, -IID_DNAm, -FID, -IID, -Batch) %>%
      summarise_if(is.numeric, list(mean = mean), na.rm = TRUE), 
    by=c("Set")) %>% as.data.frame()
# Minimal effect of set

ds_Set %>% reshape2::melt(id = c("Set")) %>% filter(Set=="wave1") %>% rename(wave1=value) %>% select(-Set) %>%
  left_join(
    ds_Set %>% reshape2::melt(id = c("Set")) %>% filter(Set=="wave3") %>% rename(wave3=value) %>% select(-Set),
    by=c("variable")) %>%
  left_join(
    ds_Set %>% reshape2::melt(id = c("Set")) %>% filter(Set=="wave4") %>% rename(wave4=value) %>% select(-Set),
    by=c("variable"))

# Large effect of Sex: Within sex, adjust for age
ds_Sex <- 
  ds %>% group_by(sex) %>% 
  select(-FID_DNAm, -IID_DNAm, -FID, -IID, -Batch, -Set) %>%
  summarise_if(is.numeric, list(mean = mean), na.rm = TRUE) %>% 
  as.data.frame()

ds_Sex %>% reshape2::melt(id = c("sex")) %>% filter(sex=="M") %>% rename(Male=value) %>% select(-sex) %>%
  left_join(
    ds_Sex %>% reshape2::melt(id = c("sex")) %>% filter(sex=="F") %>% rename(Female=value) %>% select(-sex),
    by=c("variable")) 




  
#############################
### Removing outliers
#############################

# Biochemical and anthropometric traits were trimmed for outliers (values that were ± 4 SDs from the mean). In addition BMI was trimmed for extreme values at \<17 and \>50 kg/m2.

library(RNOmni) 
library(tidyverse)
ds <- read.table("/Cluster_Filespace/Marioni_Group/Alesha/covariates/oii_cov_unrelated.txt", header=T)
# ds <- read.table("/Cluster_Filespace/Marioni_Group/Alesha/covariates/oii_cov.txt", header=T)

# 1. Remove outliers - as per Robs papers (4SD from the mean except for BMI (17-50))
ds_outliers <- ds %>% 
  mutate(
    across(
      .cols = c("height", "weight", "bmi", "waist", "hips", "whr", "body_fat", "Glucose", "HDL_cholesterol", "Total_cholesterol"),
      .fns = list(mean4sd_up = ~ mean(.x, na.rm=T) + 4*sqrt(var(.x, na.rm=T)),
                  mean4sd_low = ~ mean(.x, na.rm=T) - 4*sqrt(var(.x, na.rm=T))
      ))) 
ds_outliers %>% select(ends_with("mean4sd_up"), ends_with("mean4sd_low")) %>% unique()

ds_clean <- ds_outliers %>%
  mutate(bmi = case_when(bmi >= 17 & bmi <= 50  ~ bmi,    TRUE ~ as.numeric(NA)),
         height = case_when(height <= height_mean4sd_up & height >= height_mean4sd_low ~ height,    TRUE ~ as.numeric(NA)),
         weight = case_when(weight <= weight_mean4sd_up & weight >= weight_mean4sd_low ~ weight,    TRUE ~ as.numeric(NA)),
         waist = case_when(waist <= waist_mean4sd_up & waist >= waist_mean4sd_low ~ waist,    TRUE ~ as.numeric(NA)),
         hips = case_when(hips <= hips_mean4sd_up & hips >= hips_mean4sd_low ~ hips,    TRUE ~ as.numeric(NA)),
         whr = case_when(whr <= whr_mean4sd_up & whr >= whr_mean4sd_low ~ whr,    TRUE ~ as.numeric(NA)),
         body_fat = case_when(body_fat <= body_fat_mean4sd_up & body_fat >= body_fat_mean4sd_low ~ body_fat,    TRUE ~ as.numeric(NA)),
         Glucose = case_when(Glucose <= Glucose_mean4sd_up & Glucose >= Glucose_mean4sd_low ~ Glucose,    TRUE ~ as.numeric(NA)),
         HDL_cholesterol = case_when(HDL_cholesterol <= HDL_cholesterol_mean4sd_up & HDL_cholesterol >= HDL_cholesterol_mean4sd_low ~ HDL_cholesterol,    TRUE ~ as.numeric(NA)),
         Total_cholesterol = case_when(Total_cholesterol <= Total_cholesterol_mean4sd_up & Total_cholesterol >= Total_cholesterol_mean4sd_low ~ Total_cholesterol,    TRUE ~ as.numeric(NA)))



#############################
# Pheno Adjustment
#############################

### Adjust and normalise
# For each trait, we adjusted the phenotype for age in each gender group and standardized the residuals by RINT which removed the age effect and potential difference in mean and variance between two gender groups. There was no adjustment for wave as there were minimal differences in the mean across waves.

# Adjust each phenotype for age withine ach sex and rank normalise

ds_bmi <- ds_clean %>% 
  select(FID_DNAm, IID_DNAm, sex, age, bmi) %>%
  group_by(sex) %>% 
  filter(!is.na(bmi)) %>%
  mutate(bmi_RankNorm = residuals(lm(bmi ~ age, na.action = na.exclude)) %>% RankNorm)

ds_whr <- ds_clean %>% 
  select(FID_DNAm, IID_DNAm, sex, age, whr) %>%
  group_by(sex) %>% 
  filter(!is.na(whr)) %>%
  mutate(whr_RankNorm = residuals(lm(whr ~ age, na.action = na.exclude)) %>% RankNorm)

ds_body_fat <- ds_clean %>% 
  select(FID_DNAm, IID_DNAm, sex, age, body_fat) %>%
  group_by(sex) %>% 
  filter(!is.na(body_fat)) %>%
  mutate(body_fat_RankNorm = residuals(lm(body_fat ~ age, na.action = na.exclude)) %>% RankNorm)

ds_Glucose <- ds_clean %>% 
  select(FID_DNAm, IID_DNAm, sex, age, Glucose) %>%
  group_by(sex) %>% 
  filter(!is.na(Glucose)) %>%
  mutate(Glucose_RankNorm = residuals(lm(Glucose ~ age, na.action = na.exclude)) %>% RankNorm)

ds_HDL_cholesterol <- ds_clean %>% 
  select(FID_DNAm, IID_DNAm, sex, age, HDL_cholesterol) %>%
  group_by(sex) %>% 
  filter(!is.na(HDL_cholesterol)) %>%
  mutate(HDL_cholesterol_RankNorm = residuals(lm(HDL_cholesterol ~ age, na.action = na.exclude)) %>% RankNorm)

ds_Total_cholesterol <- ds_clean %>% 
  select(FID_DNAm, IID_DNAm, sex, age, Total_cholesterol) %>%
  group_by(sex) %>% 
  filter(!is.na(Total_cholesterol)) %>%
  mutate(Total_cholesterol_RankNorm = residuals(lm(Total_cholesterol ~ age, na.action = na.exclude)) %>% RankNorm)

ds_Height <- ds_clean %>% 
  select(FID_DNAm, IID_DNAm, sex, age, height) %>%
  group_by(sex) %>% 
  filter(!is.na(height)) %>%
  mutate(height_RankNorm = residuals(lm(height ~ age, na.action = na.exclude)) %>% RankNorm)

ds_weight <- ds_clean %>% 
  select(FID_DNAm, IID_DNAm, sex, age, weight) %>%
  group_by(sex) %>% 
  filter(!is.na(weight)) %>%
  mutate(weight_RankNorm = residuals(lm(weight ~ age, na.action = na.exclude)) %>% RankNorm)

ds_Waist <- ds_clean %>% 
  select(FID_DNAm, IID_DNAm, sex, age, waist) %>%
  group_by(sex) %>% 
  filter(!is.na(waist)) %>%
  mutate(waist_RankNorm = residuals(lm(waist ~ age, na.action = na.exclude)) %>% RankNorm)

ds_Hips <- ds_clean %>% 
  select(FID_DNAm, IID_DNAm, sex, age, hips) %>%
  group_by(sex) %>% 
  filter(!is.na(hips)) %>%
  mutate(hips_RankNorm = residuals(lm(hips ~ age, na.action = na.exclude)) %>% RankNorm)

# Merge back together
ds_pheno <- ds_clean %>% 
  left_join(ds_bmi, by=c("FID_DNAm", "IID_DNAm", "sex", "age", "bmi")) %>%
  left_join(ds_whr, by=c("FID_DNAm", "IID_DNAm", "sex", "age", "whr")) %>%
  left_join(ds_body_fat, by=c("FID_DNAm", "IID_DNAm", "sex", "age", "body_fat")) %>%
  left_join(ds_Glucose, by=c("FID_DNAm", "IID_DNAm", "sex", "age", "Glucose")) %>%
  left_join(ds_HDL_cholesterol, by=c("FID_DNAm", "IID_DNAm", "sex", "age", "HDL_cholesterol")) %>%
  left_join(ds_Total_cholesterol, by=c("FID_DNAm", "IID_DNAm", "sex", "age", "Total_cholesterol")) %>%
  left_join(ds_Height, by=c("FID_DNAm", "IID_DNAm", "sex", "age", "height")) %>%
  left_join(ds_Waist, by=c("FID_DNAm", "IID_DNAm", "sex", "age", "waist")) %>%
  left_join(ds_Hips, by=c("FID_DNAm", "IID_DNAm", "sex", "age", "hips")) %>%
  left_join(ds_weight, by=c("FID_DNAm", "IID_DNAm", "sex", "age", "weight")) %>%
  select(FID_DNAm, IID_DNAm, age, sex, Set, bmi, bmi_RankNorm, whr, whr_RankNorm, body_fat, body_fat_RankNorm, Glucose, Glucose_RankNorm,  HDL_cholesterol, HDL_cholesterol_RankNorm, Total_cholesterol, Total_cholesterol_RankNorm, height, height_RankNorm, weight, weight_RankNorm, waist, waist_RankNorm, hips, hips_RankNorm)

write.table(ds_pheno, file="/Local_Scratch/Alesha/DNAm/phenotypes/bodyfat_covariates_agesex_ranknorm_outliers4sd.txt", col.name=T, row.name=F, quote=F)


# Split for correlation analysis
ds_pheno <- read.table("/Local_Scratch/Alesha/DNAm/phenotypes/bodyfat_covariates_agesex_ranknorm_outliers4sd.txt", header=T)

# Plot again
t_ds_numeric_sex <- reshape2::melt(ds_pheno %>% select(FID_DNAm, IID_DNAm, sex, where(is.numeric), where(is.integer)), 
                                   id = c('FID_DNAm','IID_DNAm', "sex"))

t_ds_numeric_Set <- reshape2::melt(ds_pheno %>% select(FID_DNAm, IID_DNAm, Set, where(is.numeric), where(is.integer)), 
                                   id = c('FID_DNAm','IID_DNAm', "Set"))

# By sex
library(sjmisc)

t_ds_numeric_sex <- t_ds_numeric_sex %>% mutate(Method = case_when(grepl("RankNorm", variable) ~ "RINT",
                                                                   TRUE ~ "Raw"),
                                                Trait = 
                                                  case_when (variable %in% c("bmi", "bmi_RankNorm") ~ "bmi",
                                                             variable %in% c("body_fat", "body_fat_RankNorm") ~ "body_fat",
                                                             variable %in% c("whr", "whr_RankNorm") ~ "whr",
                                                             variable %in% c("Glucose", "Glucose_RankNorm") ~ "Glucose",
                                                             variable %in% c("HDL_cholesterol", "HDL_cholesterol_RankNorm") ~ "HDL_cholesterol",
                                                             variable %in% c("Total_cholesterol", "Total_cholesterol_RankNorm") ~ "Total_cholesterol"))

t_ds_numeric_sex$Trait <- recode_factor( t_ds_numeric_sex$Trait, 
                                         bmi  = "BMI", 
                                         body_fat = "Body Fat %", 
                                         whr = "Waist to Hip Ratio",
                                         Glucose = "Glucose",
                                         HDL_cholesterol = "HDL Cholesterol",
                                         Total_cholesterol = "Total Cholesterol")

p1 <- t_ds_numeric_sex %>% na.omit() %>%  ggplot(aes(x = value, fill=sex)) +
  geom_density(alpha=0.5)  +
  scale_fill_manual(values = c(mycolours[2], mycolours[1])) +
  facet_wrap( Method ~ Trait , scales = "free")


# By Set
p1 <- t_ds_numeric_Set  %>%  ggplot(aes(x = value, fill=Set)) +
  geom_density(alpha=0.5)  +
  facet_wrap( ~ variable, ncol = 4,  scales = "free")




awk '{print $1, $2, $3, $4}' oii_cov.txt > iid.txt







#############################
### Cell type proportions
#############################

# Includes technical effects and samplesheets


require(data.table)
require(tidyverse)

wave1 <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/stradl-samples-5087.rds")  # > dim(wave1) [1] 5087   15
wave3 <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave3-final/samplesheet.final.csv") # > dim(wave3) [1] 4450   13
wave4 <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave4/w4-samplesheet_v3.rds") # > dim(wave4) [1] 8877   12
wave4_cellcomp <- as.data.frame(fread("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave4/wave4_cellcomp.tsv")) # > dim(wave4_cellcomp) [1] 8877    7

comb <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds") %>% separate(col=Sample_Sentrix_ID, into=c('Sentrix_ID', 'Sentrix_Position'), sep='_', remove=FALSE)

# Merge all cell type information together
wave1_comb <- comb %>% 
  inner_join(wave1 %>% 
               mutate(Sentrix_ID = as.character(Sentrix_ID)) %>%
               select(-age, -sex, -father, -mother), 
             by=c("Sample_Name", "Sentrix_ID", "Sentrix_Position")) %>%
  select(Sample_Name, Sample_Sentrix_ID, Sentrix_ID, Sentrix_Position,
         age, sex, Set, Batch, CD8T,  CD4T,  NK, 
         Bcell , Mono,  Gran, plate_processing_batch)

wave3_comb <- comb %>%
  inner_join(wave3 %>%
               mutate(Sentrix_ID = as.character(Sentrix_ID)) %>%
               select(-age, -sex, -Batch), 
             by=c("Sample_Name", "Sample_Sentrix_ID", "Sentrix_ID", "Sentrix_Position")) %>%
  mutate(plate_processing_batch=NA) %>%
  select(Sample_Name, Sample_Sentrix_ID, Sentrix_ID, Sentrix_Position,
         age, sex, Set, Batch, CD8T,  CD4T,  NK, 
         Bcell , Mono,  Gran, plate_processing_batch)

wave4_comb <- comb %>%
  inner_join(wave4_cellcomp,
             by=c("Sample_Sentrix_ID" = "ID")) %>%
  mutate(plate_processing_batch=NA) %>%
  select(Sample_Name, Sample_Sentrix_ID, Sentrix_ID, Sentrix_Position,
         age, sex, Set, Batch, CD8T,  CD4T,  NK, 
         Bcell , Mono,  Gran, plate_processing_batch)

waves_comb <- rbind(wave1_comb, wave3_comb, wave4_comb)

write.table(waves_comb, file="/Local_Scratch/Alesha/DNAm/GS_20K_samplesheet.txt", col.name=T, row.name=F, quote=F)

ds_unrelated <- read.table("/Local_Scratch/Alesha/covariates/oii_cov_unrelated.txt", header=T)

waves_comb_unrelated <- waves_comb %>% filter(Sample_Sentrix_ID %in% ds_unrelated$IID_DNAm)
write.table(waves_comb_unrelated, file="/Local_Scratch/Alesha/DNAm/GS_20K_samplesheet_unrelated.txt", col.name=T, row.name=F, quote=F)



#############################
### **DNAm QC**
#############################

# -   GS20K 3 waves normalised together
# 
# -   We removed probes with almost invariable beta values across individuals (standard deviation \< 0.02).
# 
# -   We removed cross-hybridising/polymorphic probes based on publications by Chen *et al.* 2013 for 450k array data (<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3592906/>) or McCartney *et al.* 2016 for EPIC array data (<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4909830/>).
# 
# -   Adjust for covariates and cell types (Age, sex, wave, slide and cell type proportions)
# 
# -   Keep only unrelated individuals and exclude withdrawn individuals
# 
# -   781,379 DNAm probes


in_path=/Local_Data/methylation/GS_20k/OSCA_files
filename=mvals-norm20k-18413-831733
scratch=/Local_Scratch/Alesha

awk '{print $1, $2}' /Cluster_Filespace/Marioni_Group/Alesha/covariates/oii_cov.txt > /Cluster_Filespace/Marioni_Group/Alesha/covariates/oii_cov_related_FID_IID_DNAm.txt

osca_Linux --befile ${in_path}/${filename} \
--keep /Cluster_Filespace/Marioni_Group/Alesha/covariates/oii_cov_unrelated_FID_IID_DNAm.txt \
--exclude-probe ${scratch}/DNAm/CpGs_to_remove/exclude_epic_autosomes.txt \
--sd-min 0.02 \
--make-bod \
--out ${scratch}/DNAm/OSCA/${filename}_sd0.02_unrelated

# 0 of 781379 sites are Excluded with low variance (std < 0.020000) and 781379 sites left.
# 7758 individuals have been saved in the file /Local_Scratch/Alesha/DNAm/OSCA/mvals-norm20k-18413-831733_sd0.02_unrelated.oii .
# 781379 probes have been saved in the file /Local_Scratch/Alesha/DNAm/OSCA/mvals-norm20k-18413-831733_sd0.02_unrelated.opi .
# M-values of Methylation data for 781379 probes of 7758 individuals have been saved in the file /Local_Scratch/Alesha/DNAm/OSCA/mvals-norm20k-18413-831733_sd0.02_unrelated.bod.

# # related individuals
# <!-- 18413 individuals to be included from  /Local_Data/methylation/GS_20k/OSCA_files/mvals-norm20k-18413-831733.oii . -->
#   <!-- 18376 individuals are kept from /Cluster_Filespace/Marioni_Group/Alesha/covariates/oii_cov_related_FID_IID_DNAm.txt. -->
#   <!-- Reading probe information from /Local_Data/methylation/GS_20k/OSCA_files/mvals-norm20k-18413-831733.opi ... -->
#   <!-- WARNING: at least one gene id is missing. -->
#   <!-- 831733 probes to be included from  /Local_Data/methylation/GS_20k/OSCA_files/mvals-norm20k-18413-831733.opi . -->
#   <!-- Reading a list of probes from /Local_Scratch/Alesha/DNAm/CpGs_to_remove/exclude_epic_autosomes.txt . -->
#   <!-- 50354 probes are excluded from /Local_Scratch/Alesha/DNAm/CpGs_to_remove/exclude_epic_autosomes.txt and there are 781379 probe remaining. -->
#   <!-- Reserving 2047 MB memory for I/O buffer. -->

  
  
#############################
# Adjust probes for covariates
#############################

  # Attempt 1: Sentrix_Position, Set, Sex, age, CD8T, CD4T, NK, Bcell, Mono, Gran
  # Remove Sentrix_ID and Batch as too many levels
  
  # Sentrix_ID        Sentrix_Position Set          Batch      
  # 200784530081:   8   R02C01 : 989   wave1:2036   W3_21  : 186  
# 201057020031:   8   R06C01 : 986   wave3:4296   W3_18  : 185  
# 201057020068:   8   R04C01 : 977   wave4:1426   W3_22  : 184  
# 201057030041:   8   R08C01 : 975                W3_23  : 183  
# 201057030094:   8   R01C01 : 969                W3_25  : 183  
# 201057150040:   8   R05C01 : 962                W3_14  : 181  
# (Other)     :7710   (Other):1900                (Other):6656  

# Related individuals covariate list - age, sex, batch, cell types, smoking status, pack years
#  Qcov - "FID", "IID", "age", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "pack years"
# catcov -  "Sentrix_Position", "sex", "Batch", "smoking status"

awk -v OFS='\t' 'BEGIN { print "FID", "IID", "Sentrix_Position", "sex", "Set"} NR>1{print $2, $2,  $4, $6, $7}' /Local_Scratch/Alesha/DNAm/samplesheet/GS_20K_samplesheet_unrelated.txt >  /Local_Scratch/Alesha/DNAm/OSCA/GS_20K_covariates.txt

# Add batch - remove set? Are they colinear?
awk -v OFS='\t' 'BEGIN { print "FID", "IID", "Sentrix_Position", "sex", "Batch"} NR>1{print $2, $2,  $4, $6, $8}' /Local_Scratch/Alesha/DNAm/samplesheet/GS_20K_samplesheet_unrelated.txt >  /Local_Scratch/Alesha/DNAm/OSCA/GS_20K_covariates_batch.txt

awk -v OFS='\t' 'BEGIN { print "FID", "IID", "age", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"} NR>1{print $2, $2, $5, $9, $10, $11, $12, $13, $14}' /Local_Scratch/Alesha/DNAm/samplesheet/GS_20K_samplesheet_unrelated.txt > /Local_Scratch/Alesha/DNAm/OSCA/GS_20K_qcovariates.txt

# qcov and PCs
qcov <- read.table("/Local_Scratch/Alesha/DNAm/OSCA/GS_20K_qcovariates.txt", header=T)
path="/Local_Scratch/Alesha/DNAm/ORM/PC"
eigenvec <- read.table(paste0(path,"/mvals-norm20k-18413-831733_sd0.02_unrelated_covar-10PCS.eigenvec"))
colnames(eigenvec) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
qcov <- left_join(qcov, eigenvec, by=c("FID", "IID"))
write.table(qcov, file="/Local_Scratch/Alesha/DNAm/OSCA/GS_20K_qcovariates_PC.txt", col.name=T, row.name=F, quote=F)

# Smoking
library(tidyverse)
ever_smoke <- read.table("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/updated_smoking_jan_2019/ever_smoke.csv", header=T, sep=",")
pack_years <- read.table("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/updated_smoking_jan_2019/pack_years.csv", header=T, sep=",")
cov_batch <- read.table("/Local_Scratch/Alesha/DNAm/OSCA/GS_20K_covariates_batch.txt", header=T)
qcov <- read.table("/Local_Scratch/Alesha/DNAm/OSCA/GS_20K_qcovariates.txt", header=T)
data <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
ds_pheno <- read.table("/Local_Scratch/Alesha/DNAm/phenotypes/bodyfat_covariates_agesex_ranknorm_outliersremove.txt", header=T)

cov_batch_smoking <- cov_batch %>% 
  left_join(ever_smoke %>% 
              inner_join(data %>% select(Sample_Name, IID=Sample_Sentrix_ID), by="Sample_Name") %>% 
              select(-Sample_Name), by="IID") 
cov_batch_smoking %>% write.table(file="/Local_Scratch/Alesha/DNAm/OSCA/GS_20K_covariates_batch_smoking.txt", col.name=T, row.name=F, quote=F)

qcov_smoking <- qcov %>% 
  left_join(pack_years %>% 
              inner_join(data %>% select(Sample_Name, IID=Sample_Sentrix_ID), by="Sample_Name") %>% 
              select(-Sample_Name), by="IID") 

qcov_smoking %>% write.table(file="/Local_Scratch/Alesha/DNAm/OSCA/GS_20K_qcovariates_smoking.txt", col.name=T, row.name=F, quote=F)

cov_batch_smoking %>% drop_na() %>% select(FID, IID) %>% 
  inner_join(qcov_smoking %>% drop_na() %>% select(IID), by="IID") %>%
  write.table(file="/Local_Scratch/Alesha/DNAm/OSCA/GS_20K_keep_IID_smoking.txt", col.name=F, row.name=F, quote=F)




# Related individuals covariate list - age, sex, batch, cell types, smoking status, pack years
#  Qcov - "FID", "IID", "age", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "pack years"
# catcov -  "Sentrix_Position", "sex", "Batch", "smoking status"

library(tidyverse)
ever_smoke <- read.table("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/updated_smoking_jan_2019/ever_smoke.csv", header=T, sep=",")
pack_years <- read.table("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/updated_smoking_jan_2019/pack_years.csv", header=T, sep=",")
samplesheet <- read.table("/Local_Scratch/Alesha/DNAm/GS_20K_samplesheet_related.txt", header=T)
data <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
ds_pheno <- read.table("/Local_Scratch/Alesha/DNAm/phenotypes/bodyfat_covariates_agesex_ranknorm_outliers4sd_related.txt", header=T)

cov_batch_smoking <- samplesheet  %>%
  left_join(ever_smoke %>% 
              inner_join(data %>% select(Sample_Name, IID=Sample_Sentrix_ID), by="Sample_Name") %>% 
              select(-Sample_Name), by=c("Sample_Sentrix_ID"="IID") ) %>% 
  select(FID=Sample_Sentrix_ID, IID=Sample_Sentrix_ID, Sentrix_Position, sex,  Batch, ever_smoke)
cov_batch_smoking %>% write.table(file="/Local_Scratch/Alesha/DNAm/OSCA/GS_20K_covariates_batch_smoking_related.txt", col.name=T, row.name=F, quote=F) 

qcov_smoking <- samplesheet %>% 
  left_join(pack_years %>% 
              inner_join(data %>% select(Sample_Name, IID=Sample_Sentrix_ID), by="Sample_Name") %>% 
              select(-Sample_Name), by=c("Sample_Sentrix_ID"="IID") ) %>%
  select(FID=Sample_Sentrix_ID, IID=Sample_Sentrix_ID, age, CD8T, CD4T, NK, Bcell, Mono, Gran, pack_years)

qcov_smoking %>% write.table(file="/Local_Scratch/Alesha/DNAm/OSCA/GS_20K_qcovariates_smoking_related.txt", col.name=T, row.name=F, quote=F)


# Smoking and Weight
library(tidyverse)
ever_smoke <- read.table("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/updated_smoking_jan_2019/ever_smoke.csv", header=T, sep=",")
pack_years <- read.table("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/updated_smoking_jan_2019/pack_years.csv", header=T, sep=",")
cov_batch <- read.table("/Local_Scratch/Alesha/DNAm/OSCA/GS_20K_covariates_batch.txt", header=T)
qcov <- read.table("/Local_Scratch/Alesha/DNAm/OSCA/GS_20K_qcovariates.txt", header=T)
data <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
ds_pheno <- read.table("/Local_Scratch/Alesha/DNAm/phenotypes/bodyfat_covariates_agesex_ranknorm_outliersremove.txt", header=T)

qcov_smoking_bodyfat <- qcov %>% 
  left_join(pack_years %>% 
              inner_join(data %>% select(Sample_Name, IID=Sample_Sentrix_ID), by="Sample_Name") %>% 
              select(-Sample_Name), by="IID") %>%
  inner_join(ds_pheno %>% select(IID=IID_DNAm, body_fat), by="IID") 
qcov_smoking_bodyfat %>% write.table(file="/Local_Scratch/Alesha/DNAm/OSCA/GS_20K_qcovariates_smoking_bodyfat.txt", col.name=T, row.name=F, quote=F)

cov_batch_smoking %>% drop_na() %>% select(FID, IID) %>% 
  inner_join(qcov_smoking_bodyfat %>% drop_na() %>% select(IID), by="IID") %>%
  write.table(file="/Local_Scratch/Alesha/DNAm/OSCA/GS_20K_keep_IID_smoking_bodyfat.txt", col.name=F, row.name=F, quote=F)





# Adjust DNAm
filename=mvals-norm20k-18413-831733
scratch=/Local_Scratch/Alesha

osca_Linux --befile ${scratch}/DNAm/OSCA/${filename}_sd0.02_unrelated \
--qcovar ${scratch}/DNAm/OSCA/GS_20K_qcovariates.txt \
--covar ${scratch}/DNAm/OSCA/GS_20K_covariates.txt \
--adj-probe --make-bod \
--out ${scratch}/DNAm/OSCA/${filename}_sd0.02_unrelated_covar

osca_Linux --befile ${scratch}/DNAm/OSCA/${filename}_sd0.02_unrelated \
--qcovar ${scratch}/DNAm/OSCA/GS_20K_qcovariates.txt \
--covar ${scratch}/DNAm/OSCA/GS_20K_covariates_batch.txt \
--adj-probe --make-bod \
--out ${scratch}/DNAm/OSCA/${filename}_sd0.02_unrelated_covar_batch
# Note cant adjust for batch -  "Maybe there is multicollinearity in the covariates."
# Is this because I cant adjust for batch and set simultaneously? - removed set

# Adjust for standard covariates and PCs
osca_Linux --befile ${scratch}/DNAm/OSCA/${filename}_sd0.02_unrelated \
--qcovar ${scratch}/DNAm/OSCA/GS_20K_qcovariates_PC.txt \
--covar ${scratch}/DNAm/OSCA/GS_20K_covariates.txt \
--adj-probe --make-bod \
--out ${scratch}/DNAm/OSCA/${filename}_sd0.02_unrelated_covar_PC



# Std Cov + batch + smoking (pack years) 
# Note some individuals missing covariates so keep only those with complete information - 
osca_Linux --befile ${scratch}/DNAm/OSCA/${filename}_sd0.02_unrelated \
--qcovar ${scratch}/DNAm/OSCA/GS_20K_qcovariates_smoking.txt \
--covar ${scratch}/DNAm/OSCA/GS_20K_covariates_batch_smoking.txt \
--adj-probe --make-bod \
--thread-num 20 \
--keep ${scratch}/DNAm/OSCA/GS_20K_keep_IID_smoking.txt \
--out ${scratch}/DNAm/OSCA/${filename}_sd0.02_unrelated_covar_batch_eversmoke_packyears

# Std Cov + batch + smoking (pack years) + weight 
# Note some individuals missing covariates so keep only those with complete information - --keep extracts a subset of individuals.
osca_Linux --befile ${scratch}/DNAm/OSCA/${filename}_sd0.02_unrelated \
--qcovar ${scratch}/DNAm/OSCA/GS_20K_qcovariates_smoking_bodyfat.txt \
--covar ${scratch}/DNAm/OSCA/GS_20K_covariates_batch_smoking.txt \
--adj-probe --make-bod \
--thread-num 20 \
--keep ${scratch}/DNAm/OSCA/GS_20K_keep_IID_smoking_bodyfat.txt \
--out ${scratch}/DNAm/OSCA/${filename}_sd0.02_unrelated_covar_batch_eversmoke_packyears_bodyfat




# Related
# Std Cov + batch + smoking (pack years) 
# Note some individuals missing covariates so keep only those with complete information - 
osca_Linux --befile ${scratch}/DNAm/OSCA/${filename}_sd0.02_related \
--qcovar ${scratch}/DNAm/OSCA/GS_20K_qcovariates_smoking_related.txt \
--covar ${scratch}/DNAm/OSCA/GS_20K_covariates_batch_smoking_related.txt \
--adj-probe --make-bod \
--thread-num 20 \
--out ${scratch}/DNAm/OSCA/${filename}_sd0.02_related_covar_batch_eversmoke_packyears
# Non-missing quantitative covariates of 17831 individuals are included from /Local_Scratch/Alesha/DNAm/OSCA/GS_20K_qcovariates_smoking_related.txt.





#############################
# Checks
#############################

library(tidyverse)
ds_pheno <- read.table("/Local_Scratch/Alesha/DNAm/phenotypes/bodyfat_covariates_agesex_ranknorm_outliers4sd.txt", header=T)

# Confirm it is the same as..
# ds_pheno <- read.table("/Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds_pheno.txt", header=T) 

# And this too
# ds_pheno <- read.table("/Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds1v3.txt", header=T)

ds_pheno %>% 
  select(ends_with("RankNorm")) %>% 
  summarise(across(
    .cols = where(is.numeric), 
    .fns = list(Var=var), na.rm = TRUE, 
    .names = "{col}_{fn}"
  )) %>%
  reshape2::melt()

ds_pheno %>% 
  group_by(Set) %>% 
  select(ends_with("RankNorm")) %>% 
  summarise(across(
    .cols = where(is.numeric), 
    .fns = list(Var=var), na.rm = TRUE, 
    .names = "{col}_{fn}"
  ))%>% as.data.frame

# Check phenotypic variance matches expected & within wave
#                         variable     value
# 1               bmi_RankNorm_Var 0.9988650
# 2               whr_RankNorm_Var 0.9988578
# 3          body_fat_RankNorm_Var 0.9988474
# 4           Glucose_RankNorm_Var 0.9988308
# 5   HDL_cholesterol_RankNorm_Var 0.9988639
# 6 Total_cholesterol_RankNorm_Var 0.9988664

#     Set bmi_RankNorm_Var whr_RankNorm_Var body_fat_RankNorm_Var
# 1 wave1        1.0344391        1.0337696             0.9913359
# 2 wave3        0.9868144        0.9839824             1.0148148
# 3 wave4        0.9728171        0.9816899             0.9518543
#   Glucose_RankNorm_Var HDL_cholesterol_RankNorm_Var
# 1            1.0358339                    1.0357456
# 2            0.9623267                    0.9845110
# 3            1.0435387                    0.9873337
#   Total_cholesterol_RankNorm_Var
# 1                      1.0250887
# 2                      0.9760417
# 3                      1.0269719

# All phenotypic variances are 1






#############################
### ORM
#############################

# Make ORM adjusted for standard covariates
filename=mvals-norm20k-18413-831733_sd0.02_unrelated_covar
# Make ORM adjusted for standard covariates and PCs
# also make one adjusted for PCs
filename=mvals-norm20k-18413-831733_sd0.02_unrelated_covar_PC
# And adjusted for batch (removing set)
filename=mvals-norm20k-18413-831733_sd0.02_unrelated_covar_batch
# And adjusted for batch (removing set), ever_smoke, packyears
filename=mvals-norm20k-18413-831733_sd0.02_unrelated_covar_batch_eversmoke_packyears
# Adjusted for std cov, batch, eversmoke, packyears, bodyfat
filename=mvals-norm20k-18413-831733_sd0.02_unrelated_covar_batch_eversmoke_packyears_bodyfat

# Adjusted for std cov, batch, eversmoke, packyears, RELATED!!
# ${scratch}/DNAm/OSCA/${filename}_sd0.02_related_covar_batch_eversmoke_packyears -  Extract IID and match with GRM
filename=mvals-norm20k-18413-831733_sd0.02_related_covar_batch_eversmoke_packyears

ds_pheno <- read.table("/Local_Scratch/Alesha/DNAm/phenotypes/bodyfat_covariates_agesex_ranknorm_outliers4sd_related.txt", header=T)
GRM <- read.table("/Cluster_Filespace/Marioni_Group/Alesha/unrelated/GS_GWAS.grm.id", header=F)
ORM <- read.table("/Local_Scratch/Alesha/DNAm/OSCA/mvals-norm20k-18413-831733_sd0.02_related.oii", header=F)
data <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")

colnames(GRM) <- c("FID", "IID")
colnames(ORM) <- c("FID_DNAm", "IID_DNAm")

keep <- ds_pheno %>% select("FID_DNAm", "IID_DNAm") %>% 
  filter(IID_DNAm %in% ORM$IID_DNAm) %>%
  inner_join(data %>% select(Sample_Name, Sample_Sentrix_ID), by=c("IID_DNAm"="Sample_Sentrix_ID")) %>%
  inner_join(GRM, by=c("Sample_Name"="IID")) %>%
  select(FID_DNAm, IID_DNAm, FID, IID=Sample_Name)
write.table(keep, file="/Local_Scratch/Alesha/DNAm/OSCA/related_with_DNAm_DNA.txt", col.names=T, row.names=F, quote=F)


scratch=/Local_Scratch/Alesha
<!-- awk '{print $1, $2}' /Local_Scratch/Alesha/DNAm/OSCA/related_with_DNAm_DNA.txt > tmp -->
  
  osca_Linux  --befile ${scratch}/DNAm/OSCA/${filename} \
--make-orm --out ${scratch}/DNAm/ORM/${filename}

# Make by chr then merge together
for chr in {1..22}
do

osca_Linux  --befile ${scratch}/DNAm/OSCA/${filename} \
--chr ${chr} \
--make-orm --out ${scratch}/DNAm/ORM/${filename}_chr${chr}

done

# Make flist
echo ${scratch}/DNAm/ORM/${filename}_chr1 > ${scratch}/DNAm/ORM/myorm.flist 
for chr in {2..22}
do
echo ${scratch}/DNAm/ORM/${filename}_chr${chr} >> ${scratch}/DNAm/ORM/myorm.flist
done

# Merge orms
osca_Linux --multi-orm ${scratch}/DNAm/ORM/myorm.flist \
--make-orm \
--out ${scratch}/DNAm/ORM/${filename}



# Create sparse GRM for related individuals
scratch=/Local_Scratch/Alesha
awk '{print $3, $4}' /Local_Scratch/Alesha/DNAm/OSCA/related_with_DNAm_DNA.txt > tmp


gcta64 --grm /Cluster_Filespace/Marioni_Group/Alesha/unrelated/GS_GWAS \
--make-bK-sparse 0.05 \
--keep tmp \
--out /Cluster_Filespace/Marioni_Group/Alesha/unrelated/GS_GWAS_related_sp 


# Fix ID's for both sparse and normal GRM
# Rename GRM ID's to match DNAm ID's
# Save GRM ready for analysis
library(tidyverse)
ds <- read.table("/Cluster_Filespace/Marioni_Group/Alesha/covariates/oii_cov.txt", header=T)
id <- read.table("/Cluster_Filespace/Marioni_Group/Alesha/unrelated/GS_GWAS_related_sp.grm.id", header=F)

id %>%
  left_join(ds, by=c("V2"="IID")) %>%
  select(FID_DNAm, IID_DNAm) %>%
  write.table("/Cluster_Filespace/Marioni_Group/Alesha/unrelated/GS_GWAS_related.grm.id", col.names=F, row.names=F, quote=F)

cp /Cluster_Filespace/Marioni_Group/Alesha/unrelated/GS_GWAS.grm.N.bin /Cluster_Filespace/Marioni_Group/Alesha/unrelated/GS_GWAS_related.grm.N.bin
cp /Cluster_Filespace/Marioni_Group/Alesha/unrelated/GS_GWAS.grm.bin /Cluster_Filespace/Marioni_Group/Alesha/unrelated/GS_GWAS_related.grm.bin




# ORM for mQTLs from GoDMC and those not mQTLs
filename=mvals-norm20k-18413-831733_sd0.02_unrelated_covar_batch_eversmoke_packyears
scratch=/Local_Scratch/Alesha

osca_Linux  --befile ${scratch}/DNAm/OSCA/${filename} \
--extract-probe /Local_Scratch/Alesha/DNAm/GoDMC/cpg_1e-8.txt \
--make-orm --out ${scratch}/DNAm/ORM/${filename}_mQTLs
# 173099

osca_Linux  --befile ${scratch}/DNAm/OSCA/${filename} \
--exclude-probe /Local_Scratch/Alesha/DNAm/GoDMC/cpg_1e-8.txt \
--make-orm --out ${scratch}/DNAm/ORM/${filename}_nomQTLs --thread-num 20
# 608280



## Imputed GRM
gcta64 \
--bfile /Local_Scratch/Alesha/DNAm/ORM/GRM_GS_unrelated/imputed_HRS/GS20K_HRC.r1-1_nomono_I4_cpra 
--make-grm 
--threads 30 
--keep /Local_Scratch/Alesha/DNAm/ORM/GRM_GS_unrelated/imputed_HRS/unrelated_ID.txt 
--out /Local_Scratch/Alesha/DNAm/ORM/GRM_GS_unrelated/imputed_HRS/GS20K_HRC.r1-1_nomono_I4_cpra_unrelated 
# 24161581 SNPs to be included from BIM file. --- that seems like too many SNPs?

# Subset to HM3 SNPs
scp -P 8017 /Users/uqahatt2/Downloads/gera_hapmap3_m01_imprsq0.3_hwe1e6.bim v1ahatt4@localhost:/Local_Scratch/Alesha/DNAm/ORM/GRM_GS_unrelated/imputed_HRS/
  
  cd /Local_Scratch/Alesha/DNAm/ORM/GRM_GS_unrelated/imputed_HRS/
  
  library(tidyverse)
hm3 <- read.table("/Local_Scratch/Alesha/DNAm/ORM/GRM_GS_unrelated/imputed_HRS/gera_hapmap3_m01_imprsq0.3_hwe1e6.bim", header=F)
GS <- read.table("/Local_Scratch/Alesha/DNAm/ORM/GRM_GS_unrelated/imputed_HRS/GS20K_HRC.r1-1_nomono_I4_cpra.bim", header=F)

hm3_filter <- hm3 %>% left_join(GS, by=c("V1", "V4"))
hm3_filter %>% na.omit() %>% select(V2.y) %>% write.table("/Local_Scratch/Alesha/DNAm/ORM/GRM_GS_unrelated/imputed_HRS/hm3_snps.txt", col.names=F, row.names=F, quote=F)

gcta64 \
--bfile /Local_Scratch/Alesha/DNAm/ORM/GRM_GS_unrelated/imputed_HRS/GS20K_HRC.r1-1_nomono_I4_cpra \
--make-grm \
--threads 30 \
--extract hm3_snps.txt \
--keep /Local_Scratch/Alesha/DNAm/ORM/GRM_GS_unrelated/imputed_HRS/unrelated_ID.txt \
--out /Local_Scratch/Alesha/DNAm/ORM/GRM_GS_unrelated/imputed_HRS/GS20K_HRC.r1-1_nomono_I4_cpra_unrelated_hm3

# Check GRM from hm3 still all under 0.05

gcta64 \
--threads 30 \
--grm /Local_Scratch/Alesha/DNAm/ORM/GRM_GS_unrelated/imputed_HRS/GS20K_HRC.r1-1_nomono_I4_cpra_unrelated_hm3 \
--grm-cutoff 0.05 \
--make-grm \
--out /Local_Scratch/Alesha/DNAm/ORM/GRM_GS_unrelated/imputed_HRS/GS20K_HRC.r1-1_nomono_I4_cpra_unrelated_hm3_0.05



# Fix IDs to match DNAm IDs
cd /Local_Scratch/Alesha/DNAm/ORM/GRM_GS_unrelated/imputed_HRS
cp GS20K_HRC.r1-1_nomono_I4_cpra_unrelated.orm.id GS20K_HRC.r1-1_nomono_I4_cpra_unrelated_hm3.grm.id
cp GS20K_HRC.r1-1_nomono_I4_cpra_unrelated.orm.id GS20K_HRC.r1-1_nomono_I4_cpra_unrelated_hm3.orm.id
cp GS20K_HRC.r1-1_nomono_I4_cpra_unrelated_hm3.grm.bin GS20K_HRC.r1-1_nomono_I4_cpra_unrelated_hm3.orm.bin
cp GS20K_HRC.r1-1_nomono_I4_cpra_unrelated_hm3.grm.N.bin GS20K_HRC.r1-1_nomono_I4_cpra_unrelated_hm3.orm.N.bin
```

```{bash Rename ORM to GRM}

filename=mvals-norm20k-18413-831733_sd0.02_unrelated_covar_batch_eversmoke_packyears
filename=mvals-norm20k-18413-831733_sd0.02_unrelated_covar_batch_eversmoke_packyears_mQTLs
filename=mvals-norm20k-18413-831733_sd0.02_unrelated_covar_batch_eversmoke_packyears_nomQTLs

cp ${scratch}/DNAm/ORM/${filename}.orm.bin ${scratch}/DNAm/ORM/GRM/${filename}.grm.bin
cp ${scratch}/DNAm/ORM/${filename}.orm.id ${scratch}/DNAm/ORM/GRM/${filename}.grm.id
cp ${scratch}/DNAm/ORM/${filename}.orm.N.bin ${scratch}/DNAm/ORM/GRM/${filename}.grm.N.bin


```

```{r Smoking}
### Smoking

library(tidyverse)
ds_pheno <- read.table("/Local_Scratch/Alesha/DNAm/phenotypes/bodyfat_covariates_agesex_ranknorm_outliersremove.txt", header=T)
data <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")

ever_smoke <- read.table("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/ever_smoke.csv", header=T, sep=",")
pack_years <- read.table("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/pack_years.csv", header=T, sep=",")

data_smoke <- data %>% left_join(ever_smoke, by=c("Sample_Name")) %>%
  left_join(pack_years, by=c("Sample_Name")) %>%
  inner_join(ds_pheno, by=c("Sample_Sentrix_ID"="FID_DNAm"))

# 239 NA's when using both update_2019 and file in initial directory

# Regree variables against smoking

lm(bmi_RankNorm ~  ever_smoke + pack_years , data= data_smoke) %>% summary

Set               Batch             ever_smoke      pack_years     
Length:7758        Length:7758        Min.   :1.000   Min.   :  0.000  
Class :character   Class :character   1st Qu.:3.000   1st Qu.:  0.000  
Mode  :character   Mode  :character   Median :4.000   Median :  0.000  
Mean   :3.189   Mean   :  8.076  
3rd Qu.:4.000   3rd Qu.: 10.000  
Max.   :5.000   Max.   :133.000  
NA's   :239     NA's   :239     
```

