###################################
# Proportion of Variance Captured
###################################

# Format phenotypes, eval=FALSE}
library(tidyverse)
ds_pheno <- read.table("/Local_Scratch/Alesha/DNAm/phenotypes/bodyfat_covariates_agesex_ranknorm_outliers4sd.txt", header=T)

# Testing Sets
ds_pheno1 <- ds_pheno %>% filter(Set=="wave1") %>% 
  select(FID_DNAm, IID_DNAm, bmi_RankNorm1=bmi_RankNorm, whr_RankNorm1=whr_RankNorm, body_fat_RankNorm1=body_fat_RankNorm, Glucose_RankNorm1=Glucose_RankNorm,
         HDL_cholesterol_RankNorm1=HDL_cholesterol_RankNorm, Total_cholesterol_RankNorm1=Total_cholesterol_RankNorm, height_RankNorm1=height_RankNorm, weight_RankNorm1=weight_RankNorm, waist_RankNorm1=waist_RankNorm,  hips_RankNorm1=hips_RankNorm)

ds_pheno3 <- ds_pheno %>% filter(Set=="wave3") %>% 
  select(FID_DNAm, IID_DNAm, bmi_RankNorm3=bmi_RankNorm, whr_RankNorm3=whr_RankNorm, body_fat_RankNorm3=body_fat_RankNorm, Glucose_RankNorm3=Glucose_RankNorm,
         HDL_cholesterol_RankNorm3=HDL_cholesterol_RankNorm, Total_cholesterol_RankNorm3=Total_cholesterol_RankNorm, height_RankNorm3=height_RankNorm, weight_RankNorm3=weight_RankNorm, waist_RankNorm3=waist_RankNorm,  hips_RankNorm3=hips_RankNorm)

ds_pheno4 <- ds_pheno %>% filter(Set=="wave4") %>% 
  select(FID_DNAm, IID_DNAm, bmi_RankNorm4=bmi_RankNorm, whr_RankNorm4=whr_RankNorm, body_fat_RankNorm4=body_fat_RankNorm, Glucose_RankNorm4=Glucose_RankNorm,
         HDL_cholesterol_RankNorm4=HDL_cholesterol_RankNorm, Total_cholesterol_RankNorm4=Total_cholesterol_RankNorm, height_RankNorm4=height_RankNorm, weight_RankNorm4=weight_RankNorm, waist_RankNorm4=waist_RankNorm,  hips_RankNorm4=hips_RankNorm)

ds1v3 <- ds_pheno1 %>% full_join(ds_pheno3, by=c("FID_DNAm", "IID_DNAm")) %>%
  select(FID_DNAm, IID_DNAm, starts_with("bmi"), starts_with("whr"), starts_with("body"), 
         starts_with("Glucose"), starts_with("HDL"), starts_with("Total"), starts_with("height"),  starts_with("weight"),  starts_with("waist"),  starts_with("hips"))
write.table(ds1v3, file="/Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds1v3.txt", col.name=T, row.name=F, quote=F)

ds1v4 <- ds_pheno1 %>% full_join(ds_pheno4, by=c("FID_DNAm", "IID_DNAm")) %>%
  select(FID_DNAm, IID_DNAm, starts_with("bmi"), starts_with("whr"), starts_with("body"), 
         starts_with("Glucose"), starts_with("HDL"), starts_with("Total"), starts_with("height"),  starts_with("weight"),  starts_with("waist"),  starts_with("hips"))
write.table(ds1v4, file="/Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds1v4.txt", col.name=T, row.name=F, quote=F)

ds3v4 <- ds_pheno3 %>% full_join(ds_pheno4, by=c("FID_DNAm", "IID_DNAm")) %>%
  select(FID_DNAm, IID_DNAm, starts_with("bmi"), starts_with("whr"), starts_with("body"), 
         starts_with("Glucose"), starts_with("HDL"), starts_with("Total"), starts_with("height"),  starts_with("weight"),  starts_with("waist"),  starts_with("hips"))
write.table(ds3v4, file="/Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds3v4.txt", col.name=T, row.name=F, quote=F)


write.table(ds_pheno1, file="/Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds_wave1.txt", col.name=T, row.name=F, quote=F)
write.table(ds_pheno3, file="/Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds_wave3.txt", col.name=T, row.name=F, quote=F)
write.table(ds_pheno4, file="/Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds_wave4.txt", col.name=T, row.name=F, quote=F)



##################
# Split by Gender
#################

# Male
ds_phenoM <- ds_pheno %>% filter(sex=="M") %>% 
  select(FID_DNAm, IID_DNAm, bmi_RankNormM=bmi_RankNorm, whr_RankNormM=whr_RankNorm, body_fat_RankNormM=body_fat_RankNorm, Glucose_RankNormM=Glucose_RankNorm,
         HDL_cholesterol_RankNormM=HDL_cholesterol_RankNorm, Total_cholesterol_RankNormM=Total_cholesterol_RankNorm, height_RankNormM=height_RankNorm, weight_RankNormM=weight_RankNorm, waist_RankNormM=waist_RankNorm,  hips_RankNormM=hips_RankNorm)

ds_phenoF <- ds_pheno %>% filter(sex=="F") %>% 
  select(FID_DNAm, IID_DNAm, bmi_RankNormF=bmi_RankNorm, whr_RankNormF=whr_RankNorm, body_fat_RankNormF=body_fat_RankNorm, Glucose_RankNormF=Glucose_RankNorm,
         HDL_cholesterol_RankNormF=HDL_cholesterol_RankNorm, Total_cholesterol_RankNormF=Total_cholesterol_RankNorm, height_RankNormF=height_RankNorm, weight_RankNormF=weight_RankNorm, waist_RankNormF=waist_RankNorm,  hips_RankNormF=hips_RankNorm)

dsMvF <- ds_phenoM %>% full_join(ds_phenoF, by=c("FID_DNAm", "IID_DNAm")) %>%
  select(FID_DNAm, IID_DNAm, starts_with("bmi"), starts_with("whr"), starts_with("body"), 
         starts_with("Glucose"), starts_with("HDL"), starts_with("Total"), starts_with("height"), starts_with("weight"), starts_with("waist"),  starts_with("hips"))
write.table(dsMvF, file="/Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/dsMvF.txt", col.name=T, row.name=F, quote=F)

write.table(ds_phenoM, file="/Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds_pheno_Males.txt", col.name=T, row.name=F, quote=F)
write.table(ds_phenoF, file="/Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds_pheno_Females.txt", col.name=T, row.name=F, quote=F)



#################
# Analysis Sets
#################
ds_pheno %>% select(FID_DNAm, IID_DNAm, bmi_RankNorm, whr_RankNorm, body_fat_RankNorm ,  Glucose_RankNorm, HDL_cholesterol_RankNorm, Total_cholesterol_RankNorm, height_RankNorm, weight_RankNorm, waist_RankNorm, hips_RankNorm) %>%
  write.table(file="/Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds_pheno.txt", col.name=T, row.name=F, quote=F)

# Analysis Sets - related
ds_pheno <- read.table("/Local_Scratch/Alesha/DNAm/phenotypes/bodyfat_covariates_agesex_ranknorm_outliers4sd_related.txt", header=T)
ds_pheno %>% select(FID_DNAm, IID_DNAm, bmi_RankNorm, whr_RankNorm, body_fat_RankNorm ,  Glucose_RankNorm, HDL_cholesterol_RankNorm, Total_cholesterol_RankNorm, height_RankNorm, weight_RankNorm, waist_RankNorm, hips_RankNorm) %>%
  write.table(file="/Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds_pheno_related.txt", col.name=T, row.name=F, quote=F)



#################
# Descriptive tables
#################
library(tidyverse)
ds_pheno <- read.table("/Local_Scratch/Alesha/DNAm/phenotypes/bodyfat_covariates_agesex_ranknorm_outliers4sd.txt", header=T)
ORM <- read.table("/Local_Scratch/Alesha/DNAm/ORM/mvals-norm20k-18413-831733_sd0.02_unrelated_covar_batch_eversmoke_packyears.orm.id", header=F)
data <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")

ds <- ds_pheno %>% filter(IID_DNAm %in% ORM$V2) %>% left_join(data %>% select(Sample_Sentrix_ID, age), by=c("IID_DNAm"="Sample_Sentrix_ID"))





###############################
# Variance Component Analysis 
###############################

# We estimate te proportion of variance captured by DNAm probes and 
# compare this with the proportion of variance explained by genetic factors.

# Adjusting DNAm for age, sex, cell types, batch, smoking status and pack years

pheno=/Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds_pheno.txt
scratch=/Local_Scratch/Alesha

cov=covar_batch_eversmoke_packyears
filename=mvals-norm20k-18413-831733_sd0.02_unrelated_${cov}

cd ${scratch}/DNAm/phenotypes/greml.pheno/temp
mkdir ${scratch}/DNAm/GREML/prop_var_explained/${cov}

for trait in bmi whr body_fat Glucose HDL_cholesterol Total_cholesterol height weight waist hips
do

if [ "${trait}" == "bmi" ];              then  var=3; fi
if [ "${trait}" == "whr" ];              then  var=4; fi
if [ "${trait}" == "body_fat" ];         then  var=5; fi
if [ "${trait}" == "Glucose" ];          then  var=6; fi
if [ "${trait}" == "HDL_cholesterol" ];  then  var=7; fi
if [ "${trait}" == "Total_cholesterol" ]; then  var=8; fi
if [ "${trait}" == "height" ]; then  var=9; fi
if [ "${trait}" == "weight" ]; then  var=10; fi
if [ "${trait}" == "waist" ]; then  var=11; fi
if [ "${trait}" == "hips" ]; then  var=12; fi


awk -v "col1=${var}" '{print $1, $2, $col1}' ${pheno} > temp.${trait}.prop

osca_Linux --reml \
--orm ${scratch}/DNAm/ORM/${filename} \
--pheno temp.${trait}.prop \
--out ${scratch}/DNAm/GREML/prop_var_explained/${cov}/${trait}

done


pheno=/Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds_pheno.txt

scratch=/Local_Scratch/Alesha

mkdir ${scratch}/DNAm/GREML/prop_var_explained/GRM/
  
  cd ${scratch}/DNAm/phenotypes/greml.pheno/temp

cov=covar_batch_eversmoke_packyears
mkdir ${scratch}/DNAm/GREML/prop_var_explained/${cov}_GRM

for trait in bmi whr body_fat Glucose HDL_cholesterol Total_cholesterol height 
do

if [ "${trait}" == "bmi" ];              then  var=3; fi
if [ "${trait}" == "whr" ];              then  var=4; fi
if [ "${trait}" == "body_fat" ];         then  var=5; fi
if [ "${trait}" == "Glucose" ];          then  var=6; fi
if [ "${trait}" == "HDL_cholesterol" ];  then  var=7; fi
if [ "${trait}" == "Total_cholesterol" ]; then  var=8; fi
if [ "${trait}" == "height" ]; then  var=9; fi
if [ "${trait}" == "weight" ]; then  var=10; fi
if [ "${trait}" == "waist" ]; then  var=11; fi
if [ "${trait}" == "hips" ]; then  var=12; fi

awk -v "col1=${var}" '{print $1, $2, $col1}' ${pheno} > temp.${trait}.prop


gcta64 --reml --grm /Cluster_Filespace/Marioni_Group/Alesha/unrelated/GS_GWAS_unrelated \
--pheno temp.${trait}.prop \
--threads 10 \
--out ${scratch}/DNAm/GREML/prop_var_explained/${cov}_GRM/unrelated_${trait}


osca_Linux --reml --orm /Local_Scratch/Alesha/DNAm/ORM/GRM_GS_unrelated/imputed_HRS/GS20K_HRC.r1-1_nomono_I4_cpra_unrelated_hm3 \
--pheno temp.${trait}.prop \
--thread-num 10 \
--out ${scratch}/DNAm/GREML/prop_var_explained/${cov}_GRM/${trait}_imputed_hm3

done

###############
# Joint
###############

pheno=/Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds_pheno.txt
scratch=/Local_Scratch/Alesha

cov=covar_batch_eversmoke_packyears
filename=mvals-norm20k-18413-831733_sd0.02_unrelated_${cov}

mkdir ${scratch}/DNAm/GREML/prop_var_explained/GRM_ORM/
  mkdir ${scratch}/DNAm/GREML/prop_var_explained/GRM_ORM/${cov}

echo ${scratch}/DNAm/ORM/GRM/${filename} > ${scratch}/DNAm/GREML/prop_var_explained/GRM_ORM/${cov}/GRM_ORM_${cov}.txt
echo ${scratch}/DNAm/GREML/traits/GRM_ORM/GS_GWAS_unrelated >> ${scratch}/DNAm/GREML/prop_var_explained/GRM_ORM/${cov}/GRM_ORM_${cov}.txt


echo ${scratch}/DNAm/ORM/${filename} > ${scratch}/DNAm/GREML/prop_var_explained/GRM_ORM/${cov}/GRM_ORM_${cov}_imputed_hm3.txt
echo /Local_Scratch/Alesha/DNAm/ORM/GRM_GS_unrelated/imputed_HRS/GS20K_HRC.r1-1_nomono_I4_cpra_unrelated_hm3 >> ${scratch}/DNAm/GREML/prop_var_explained/GRM_ORM/${cov}/GRM_ORM_${cov}_imputed_hm3.txt



cd ${scratch}/DNAm/phenotypes/greml.pheno/temp

for trait in bmi whr body_fat Glucose HDL_cholesterol Total_cholesterol height  
do

if [ "${trait}" == "bmi" ];              then  var=3; fi
if [ "${trait}" == "whr" ];              then  var=4; fi
if [ "${trait}" == "body_fat" ];         then  var=5; fi
if [ "${trait}" == "Glucose" ];          then  var=6; fi
if [ "${trait}" == "HDL_cholesterol" ];  then  var=7; fi
if [ "${trait}" == "Total_cholesterol" ]; then  var=8; fi
if [ "${trait}" == "height" ]; then  var=9; fi

awk -v "col1=${var}" '{print $1, $2, $col1}' ${pheno} > temp.${trait}.joint

osca_Linux --reml \
--multi-orm ${scratch}/DNAm/GREML/prop_var_explained/GRM_ORM/${cov}/GRM_ORM_${cov}_imputed_hm3.txt \
--pheno temp.${trait}.joint \
--thread-num 30 \
--out ${scratch}/DNAm/GREML/prop_var_explained/GRM_ORM/${cov}/${trait}_imputed_hm3

done


path=${scratch}/DNAm/GREML/prop_var_explained/covar_batch_eversmoke_packyears
cd ${path}

awk  -F',' '{print FILENAME, "allsamples", "ORM",  $1, $2, $3, $4}' *.rsq  > ${path}/traits_imputed_hm3.hsq
cd ../${cov}_GRM
awk  -F',' '{print FILENAME, "allsamples", "GRM",  $1, $2, $3, $4}' *imputed_hm3.rsq  >> ${path}/traits_imputed_hm3.hsq

cd ../GRM_ORM/covar_batch_eversmoke_packyears
awk  -F',' '{print FILENAME, "allsamples", "Joint",  $1, $2, $3, $4}' *imputed_hm3.rsq  >> ${path}/traits_imputed_hm3.hsq



###############
# Output
#############
 
# Order is not importance. Same results are obtained either way e.g. GRM + ORM vs ORM + GRM Compared to literature both V(O) and V(G) are overestimated
library(tidyverse)
library(stringr)
library(readr)
library(xlsx)
traits <- read_table("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/OREML_propvar/covar_batch_eversmoke_packyears/traits_imputed_hm3.hsq") %>% 
  select(-X7) %>% 
  filter(Variance!="Variance" & !is.na(SE)) %>%
  mutate(Var = case_when(Source %in% c("V(O)","V(O1)") & ORM %in% c("ORM", "Joint") ~ "VO",
                         Source %in% c("V(O)", "V(O2)") & ORM %in% c("GRM", "Joint") ~ "Vg",
                         Source == "Vp" ~ "Vp",
                         Source == "V(e)" ~ "Ve",
                         Source %in% c("V(O)/Vp", "V(O1)/Vp") & ORM %in% c("ORM", "Joint") ~ "VO_Vp",
                         Source %in% c("V(O)/Vp", "V(O2)/Vp") & ORM %in% c("GRM", "Joint") ~ "Vg_Vp",
                         TRUE ~ Source)) %>%
  select(-Source, -allsamples) 
colnames(traits)[1] <- "Trait"
colnames(traits)[2] <- "Source"
traits$Trait <- str_extract(traits$Trait, ".*(?=\\.)")
traits$Trait <- recode_factor( traits$Trait, 
                               bmi  = "BMI", 
                               body_fat = "Body Fat %", 
                               whr = "Waist to Hip Ratio",
                               Glucose = "Glucose",
                               HDL_cholesterol = "HDL Cholesterol",
                               Total_cholesterol = "Total Cholesterol",
                               height = "Height",
                               weight = "Weight",
                               waist = "Waist Circumference",
                               hips = "Hip Circumference",
                               bmi_imputed_hm3  = "BMI", 
                               body_fat_imputed_hm3 = "Body Fat %", 
                               whr_imputed_hm3 = "Waist to Hip Ratio",
                               Glucose_imputed_hm3 = "Glucose",
                               HDL_cholesterol_imputed_hm3 = "HDL Cholesterol",
                               Total_cholesterol_imputed_hm3 = "Total Cholesterol",
                               height_imputed_hm3 = "Height")
# Recode to numeric
traits$Variance <- traits$Variance %>% as.numeric
traits$SE <- traits$SE %>% as.numeric

traits <- traits %>% mutate(Group = paste0(Source,"_", Var))

traits_wide <- traits %>%
  select(-Source, -Var) %>%
  nest(col_value = c(Variance, SE)) %>%
  spread(key = Group, value = col_value) %>%
  unnest(c(ORM_Ve, ORM_VO, ORM_VO_Vp, ORM_Vp,
           GRM_Ve, GRM_Vg, GRM_Vg_Vp, GRM_Vp,
           Joint_Ve, Joint_Vg, Joint_Vg_Vp, Joint_VO, Joint_VO_Vp, Joint_Vp),  .sep = '_') %>%
  select(Trait, 
         ORM_VO_Variance, ORM_VO_SE, ORM_Vp_Variance, ORM_Vp_SE, 
         ORM_VO_Vp_Variance, ORM_VO_Vp_SE, ORM_Ve_Variance, ORM_Ve_SE,
         GRM_Vg_Variance, GRM_Vg_SE, GRM_Vp_Variance, GRM_Vp_SE,
         GRM_Vg_Vp_Variance, GRM_Vg_Vp_SE, GRM_Ve_Variance, GRM_Ve_SE,
         Joint_VO_Variance, Joint_VO_SE, Joint_Vg_Variance, Joint_Vg_SE,
         Joint_Vp_Variance, Joint_Vp_SE, Joint_VO_Vp_Variance, Joint_VO_Vp_SE, 
         Joint_Vg_Vp_Variance, Joint_Vg_Vp_SE, Joint_Ve_Variance, Joint_Ve_SE)

traits_wide %>% 
  filter(!Trait %in% c("Weight", "Waist Circumference", "Hip Circumference")) %>%
  write.xlsx(file="/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/results/epigenetic_correlations.xlsx",
             sheetName="Proportion_variance", append=TRUE)


traits_wide_table <- traits %>%
  mutate(Var_source = paste0(Var,"_",Source)) %>%
  filter(Var %in% c("VO_Vp", "Vg_Vp")) %>%
  select(-Source, -Var) %>%
  nest(col_value = c(Variance, SE)) %>%
  spread(key = Var_source, value = col_value) %>%
  unnest(`Vg_Vp_GRM`, `Vg_Vp_Joint`, `VO_Vp_Joint`, `VO_Vp_ORM`, .sep = '_') %>%
  mutate(across(where(is.character), is.numeric)) %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  mutate(`VO/Vp` = paste0(VO_Vp_ORM_Variance, " (", VO_Vp_ORM_SE, ")"),
         `Vg/Vp` = paste0(Vg_Vp_GRM_Variance, " (", Vg_Vp_GRM_SE, ")"),
         `VO/Vp Joint` = paste0(VO_Vp_Joint_Variance, " (", VO_Vp_Joint_SE, ")"),
         `Vg/Vp Joint` = paste0(Vg_Vp_Joint_Variance, " (", Vg_Vp_Joint_SE, ")")) 

kable(traits_wide_table %>% 
        select(Trait, 
               `VO/Vp`, 
               `Vg/Vp`,
               `VO/Vp Joint`,
               `Vg/Vp Joint`), caption="Proportion of Variance captured") 



figure_ds <- traits %>% mutate(Var_source = paste0(Var,"_",Source)) %>%
  filter(Var %in% c("VO_Vp", "Vg_Vp")) %>%
  mutate(Var_comp  = case_when (Var_source == "VO_Vp_ORM" ~ "MRM",
                                Var_source == "Vg_Vp_GRM" ~ "GRM",
                                Var_source == "VO_Vp_Joint" ~ "MRM Joint",
                                Var_source == "Vg_Vp_Joint" ~ "GRM Joint" ) %>% as.factor)  %>%
  mutate(Source  = case_when (Var_source == "VO_Vp_ORM" ~ "Marginal",
                              Var_source == "Vg_Vp_GRM" ~ "Marginal",
                              Var_source == "VO_Vp_Joint" ~ "Joint",
                              Var_source == "Vg_Vp_Joint" ~ "Joint" ) %>% as.factor)
figure_ds$Source <- recode_factor( figure_ds$Source, 
                                   Marginal  = "Marginal", 
                                   Joint = "Joint")
figure_ds$Var_comp <- recode_factor( figure_ds$Var_comp, 
                                     MRM  = "MRM", 
                                     GRM = "GRM",
                                     `MRM Joint`="MRM Joint",
                                     `GRM Joint` ="GRM Joint")

figure_ds <- figure_ds %>% filter(Trait %in% c("BMI", "Body Fat %", "Waist to Hip Ratio", "Glucose", "HDL Cholesterol", "Total Cholesterol"))


figure_ds %>% ggplot(aes(fill=Var_comp, y=Variance, x=Source)) +
  geom_bar( stat="identity", width=0.8, position = position_dodge(width=0.9)) +
  geom_errorbar(aes(ymin=Variance-SE, ymax=Variance+SE), width=.2,
                position=position_dodge(0.9))  +
  facet_wrap(~Trait , nrow = 2) +
  theme_linedraw() +
  scale_fill_manual(values = mycolours, name = "Variance component") +
  labs(y="Proportion of variance captured", x="Model") 
# theme(legend.position="bottom")

