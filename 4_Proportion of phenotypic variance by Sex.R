
########################
### Males and Females
########################

# bash Proportion of variance for Males and Females}

scratch=/Local_Scratch/Alesha
cov=covar_batch_eversmoke_packyears
filename=mvals-norm20k-18413-831733_sd0.02_unrelated_${cov}

mkdir ${scratch}/DNAm/GREML/prop_var_explained/Male
mkdir ${scratch}/DNAm/GREML/prop_var_explained/Female

cd ${scratch}/DNAm/phenotypes/greml.pheno/temp


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
# Male
awk -v "col1=${var}" '{print $1, $2, $col1}' /Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds_pheno_Males.txt > temp.${trait}_Male.prop

osca_Linux --reml \
--orm ${scratch}/DNAm/ORM/${filename} \
--pheno temp.${trait}_Male.prop \
--out ${scratch}/DNAm/GREML/prop_var_explained/Male/${trait}_Male


# Female
awk -v "col1=${var}" '{print $1, $2, $col1}' /Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds_pheno_Females.txt > temp.${trait}_Female.prop

osca_Linux --reml \
--orm ${scratch}/DNAm/ORM/${filename} \
--pheno temp.${trait}_Female.prop \
--out ${scratch}/DNAm/GREML/prop_var_explained/Female/${trait}_Female

done


sex=Female
cd ${scratch}/DNAm/GREML/prop_var_explained/${sex}
rm  trait_${sex}.rsq
for trait in BMI whr body_fat Glucose HDL_cholesterol Total_cholesterol height weight waist hips
do
awk -v "trait=${trait}" -v sex=${sex} '{print trait, sex, $1, $2, $3}' ${trait}_${sex}.rsq >> trait_${sex}.rsq
done



#################
# Output 
################
library(stringr)
library(readr)
traits_F <- read_table("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/GREML_sex/trait_Female.rsq") %>% 
  mutate(Var = case_when(Source %in% c("V(O)") ~ "VO",
                         Source %in% c("Vp") ~ "Vp",
                         Source %in% c("rG") ~ "rG",
                         Source == "V(O)/Vp" ~ "VO/Vp")) %>% 
  select(-Source) %>%
  filter(!is.na(Var)) 

traits_M <- read_table("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/GREML_sex/trait_Male.rsq") %>% 
  mutate(Var = case_when(Source %in% c("V(O)") ~ "VO",
                         Source %in% c("Vp") ~ "Vp",
                         Source %in% c("rG") ~ "rG",
                         Source == "V(O)/Vp" ~ "VO/Vp")) %>% 
  select(-Source) %>%
  filter(!is.na(Var)) 
colnames(traits_F)[1] <- "Trait"
colnames(traits_F)[2] <- "Gender"
colnames(traits_M)[1] <- "Trait"
colnames(traits_M)[2] <- "Gender"

traits <- rbind(traits_F, traits_M)

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
                               hips = "Hip Circumference")

traits$Gender <- recode_factor( traits$Gender, 
                                Female  = "Female", 
                                Male = "Male")

# Recode to numeric
traits$Variance <- traits$Variance %>% as.numeric
traits$SE <- traits$SE %>% as.numeric

traits_wide <- traits %>%
  nest(col_value = c(Variance, SE)) %>%
  spread(key = Var, value = col_value) %>%
  unnest(c(VO, Vp, `VO/Vp`, Ve),  .sep = '_')  %>%
  select(Trait, Gender,
         VO_Variance, VO_SE, 
         Vp_Variance, Vp_SE, 
         Ve_Variance, Ve_SE,
         `VO/Vp_Variance`, `VO/Vp_SE`)
traits_wide %>% write.xlsx(file="/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/results/traits_proportion_of_variance_explained_waves.xlsx")


traits_wide_table <- traits %>%
  nest(col_value = c(Variance, SE)) %>%
  spread(key = Var, value = col_value) %>%
  unnest(c(`VO`, `Vp`, `VO/Vp`),  .sep = '_') %>%
  mutate(across(where(is.character), is.numeric)) %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  mutate(`VO` = paste0(`VO_Variance`, " (", `VO_SE`, ")"),
         `Vp` = paste0(`Vp_Variance`, " (", `Vp_SE`, ")"),
         `VO/Vp` = paste0(`VO/Vp_Variance`, " (", `VO/Vp_SE`, ")"))



figure_ds <- traits %>% 
  filter(Trait %in% c("BMI", "Body Fat %", "Waist to Hip Ratio", "Glucose", "HDL Cholesterol", "Total Cholesterol") &
           Var =="VO/Vp")

library(ggpubr)
ggarrange(
  figure_ds %>% ggplot(aes(fill=Gender, y=Variance, x=Gender)) +
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(ymin=Variance-SE, ymax=Variance+SE), width=.2,
                  position=position_dodge(.9))  +
    facet_wrap(~Trait, nrow=2) +
    theme_linedraw() +
    scale_fill_manual(values = mycolours, name = "Sex") +
    labs(y="Proportion of variance captured by DNAm", x="Sex") +
    theme(legend.position = "none"),
  figure_sex, ncol=2, widths=c(0.6, 0.4))


