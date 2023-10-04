#########################
### **Across waves**
##########################

# bash, Waves - Proportion of variance captured}

scratch=/Local_Scratch/Alesha
cov=covar_batch_eversmoke_packyears
filename=mvals-norm20k-18413-831733_sd0.02_unrelated_${cov}

cd ${scratch}/DNAm/phenotypes/greml.pheno/temp

mkdir ${scratch}/DNAm/GREML/prop_var_explained/waves/${cov}
wave=wave4

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

awk -v "col1=${var}" '{print $1, $2, $col1}' /Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds_${wave}.txt > temp.${trait}_${wave}.prop

osca_Linux --reml \
--orm ${scratch}/DNAm/ORM/${filename} \
--pheno temp.${trait}_${wave}.prop \ 
--keep /Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/${wave}.txt \ 
--out ${scratch}/DNAm/GREML/prop_var_explained/waves/${cov}/${trait}_${wave}
  
gcta64 --reml --grm /Cluster_Filespace/Marioni_Group/Alesha/unrelated/GS_GWAS_related \
--pheno temp.${trait}_${wave}.prop \
--threads 10 \
--out ${scratch}/DNAm/GREML/prop_var_explained/waves/${cov}/${trait}_${wave}_GRM

done



cd ${scratch}/GREML/prop_var_explained/waves/${cov}
for trait in bmi whr body_fat Glucose HDL_cholesterol Total_cholesterol height weight waist hips
do
for wave in wave1 wave3
do
awk -v "trait=${trait}" -v wave=${wave} '{print trait, wave, $1, $2, $3}' ${trait}_${wave}_1500.rsq >> trait_waves_1500.rsq
done
done




#################
# Ouptut
################

library(stringr)
library(readr)
traits <- read_table("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/OREML_propvar/covar_batch_eversmoke_packyears/trait_waves.rsq") %>% 
  filter(!is.na(SE) & SE!="SE" & Variance!="of" & !(Source %in% c("logL", "logL0", "LRT", "df", "Pval"))) %>%
  mutate(Var = case_when(Source %in% c("V(O)","V(G1)") ~ "VO",
                         Source == "Vp" ~ "Vp",
                         Source == "V(e)" ~ "Ve",
                         Source %in% c("V(O)/Vp", "V(G1)/Vp") ~ "VO_Vp")) %>%
  select(-Source) 
colnames(traits)[1] <- "Trait"
colnames(traits)[2] <- "Wave"

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

traits$Set <- recode_factor( traits$Wave, 
                             wave1  = "Set 1", 
                             wave3 = "Set 2", 
                             wave4 = "Set 3")
# Recode to numeric
traits$Variance <- traits$Variance %>% as.numeric
traits$SE <- traits$SE %>% as.numeric

traits_wide <- traits %>%
  nest(col_value = c(Variance, SE)) %>%
  spread(key = Var, value = col_value) %>%
  unnest(c(VO, Vp, VO_Vp, Ve),  .sep = '_')  %>%
  select(Trait, Wave,
         VO_Variance, VO_SE, 
         Vp_Variance, Vp_SE, 
         Ve_Variance, Ve_SE,
         `VO_Vp_Variance`, `VO_Vp_SE`)
traits_wide %>% write.xlsx(file="/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/results/traits_proportion_of_variance_explained_waves.xlsx")

traits_wide_table <- traits %>%
  nest(col_value = c(Variance, SE)) %>%
  spread(key = Var, value = col_value) %>%
  unnest(`VO`, `Vp`, `VO_Vp`, `Ve`,  .sep = '_') %>%
  mutate(across(where(is.character), is.numeric)) %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  mutate(`VO` = paste0(VO_Variance, " (", VO_SE, ")"),
         `Vp` = paste0(Vp_Variance, " (", Vp_SE, ")"),
         `VO/Vp` = paste0(VO_Vp_Variance, " (", VO_Vp_SE, ")"))

kable(traits_wide_table %>% 
        select(Trait, Wave, 
               `VO`, 
               `Vp`,
               `VO/Vp`), caption="Proportion of Variance captured by wave")

# Figure
figure_ds <- traits %>% 
  filter(Var %in% c("VO_Vp")) %>%
  mutate(Set = as.factor(Set))

figure_ds <- figure_ds %>% filter(Trait %in% c("BMI", "Body Fat %", "Waist to Hip Ratio", "Glucose", "HDL Cholesterol", "Total Cholesterol"))

figure_ds %>% ggplot(aes(fill=Set, y=Variance, x=Set)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=Variance-SE, ymax=Variance+SE), width=.2,
                position=position_dodge(.9))  +
  facet_wrap(~Trait, nrow=2) +
  theme_linedraw() +
  scale_fill_manual(values = mycolours, name = "Set") +
  labs(y="Proportion variance captured", x=" ") 
# theme(legend.position="bottom")




# Check numbers by batch
library(tidyverse)
cov_batch <- read.table("/Local_Scratch/Alesha/DNAm/OSCA/GS_20K_covariates_batch.txt", header=T)
ds_pheno <- read.table("/Local_Scratch/Alesha/DNAm/phenotypes/bodyfat_covariates_agesex_ranknorm_outliersremove.txt", header=T)

ds_batch <- ds_pheno %>% 
  left_join(cov_batch %>% select(IID, Batch), by=c("IID_DNAm"="IID"))  %>%
  group_by(Set, Batch) %>% tally()
ds_batch %>% write.table(file="/Local_Scratch/Alesha/DNAm/phenotypes/descriptive/batch_numbers.txt", col.name=T, row.name=F, quote=F)


# Read
library(tidyverse)
ds_batch <- read.table("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/phenotypes/batch_numbers.txt", header=T)

ds_batch$Set <- recode_factor(ds_batch$Set, 
                              wave1  = "Set 1", 
                              wave3 = "Set 2", 
                              wave4 = "Set 3")

ds_batch %>%
  ggplot(aes(fill=Set, y=n, x=Batch)) +
  geom_bar(position="dodge", stat="identity")  +
  theme_linedraw() +
  scale_fill_manual(values = mycolours, name = "Set") +
  labs(y="Samples per batch", x="Batch")  +
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(expand = c(0, 0))

  




#####################################################################
# EWAS for each trait to compare the marginal effects across waves
#####################################################################

pheno=/Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds_pheno.txt
scratch=/Local_Scratch/Alesha
cov=covar_batch_eversmoke_packyears
filename=mvals-norm20k-18413-831733_sd0.02_unrelated_${cov}


cd ${scratch}/DNAm/phenotypes/greml.pheno/temp
wave=wave4

mkdir ${scratch}/DNAm/GREML/EWAS

for trait in bmi whr body_fat Glucose HDL_cholesterol Total_cholesterol height
do

if [ "${trait}" == "bmi" ];              then  var=3; fi
if [ "${trait}" == "whr" ];              then  var=4; fi
if [ "${trait}" == "body_fat" ];         then  var=5; fi
if [ "${trait}" == "Glucose" ];          then  var=6; fi
if [ "${trait}" == "HDL_cholesterol" ];  then  var=7; fi
if [ "${trait}" == "Total_cholesterol" ]; then  var=8; fi
if [ "${trait}" == "height" ]; then  var=9; fi

awk -v "col1=${var}" '{print $1, $2, $col1}' /Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds_${wave}.txt > temp.${trait}_${wave}.EWAS

osca_Linux  \
--befile ${scratch}/DNAm/OSCA/${filename} \
--pheno temp.${trait}_${wave}.EWAS \
--linear \
--out ${scratch}/DNAm/GREML/EWAS/waves/${trait}_${wave}

done

# Filter out SNPs where p < 1e-4
scratch=/Local_Scratch/Alesha

for trait in bmi whr body_fat Glucose HDL_cholesterol Total_cholesterol height
do

wave=wave1
# awk -v "wave_var=${wave}" '{ if($8  < 0.001) { print $0, wave_var}}' ${scratch}/DNAm/GREML/EWAS/waves/${trait}_${wave}.linear > ${scratch}/DNAm/GREML/EWAS/waves/${trait}_0.001.linear
awk -v "wave_var=${wave}" '{{ print $0, wave_var}}' ${scratch}/DNAm/GREML/EWAS/waves/${trait}_${wave}.linear > ${scratch}/DNAm/GREML/EWAS/waves/${trait}.linear

for wave in wave3 wave4
do

awk -v "wave_var=${wave}" '{ { print $0, wave_var}}' ${scratch}/DNAm/GREML/EWAS/waves/${trait}_${wave}.linear >> ${scratch}/DNAm/GREML/EWAS/waves/${trait}.linear

done 
done





###########
# EWAS_r
###########

marginal_effects <- lapply(1:6, function(trait_n) {
  
  trait_list=c("bmi", "whr", "body_fat", "Glucose", "HDL_cholesterol", "Total_cholesterol", "height")
  trait=trait_list[trait_n]
  
  ds <- read.table(paste0("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/EWAS/",trait,".linear"), header=F)
  
  colnames(ds) <- c("Chr",	"Probe", "bp",	"Gene",	"Orientation",	"b",	"se",	"p",	"NMISS", "wave")
  ds <- ds %>% filter(Chr!="Chr") %>% mutate(b = as.numeric(b), se=as.numeric(se), p=as.numeric(p))
  
  ds_wide <- reshape(data=ds %>%  select(Chr, Probe, b, se, p, wave), 
                     idvar=c("Probe", "Chr"), v.names = c("b", "se", "p"), 
                     timevar = "wave", 
                     direction="wide") %>%
    na.exclude %>% mutate(trait=paste0(trait))
  ds_0.001_wide <- ds_wide %>% filter( p.wave1 < 0.001 & p.wave3 < 0.001 & p.wave4 < 0.001 )
  
  ds_0.001_wide_clean <- rbind(
    ds_0.001_wide %>% 
      select(Chr, Probe, trait, x=b.wave1, y=b.wave3) %>% 
      mutate(comparison="wave 1 vs 3"),
    ds_0.001_wide %>% 
      select(Chr, Probe, trait, x=b.wave1, y=b.wave4) %>% 
      mutate(comparison="wave 1 vs 4"),  
    ds_0.001_wide %>% 
      select(Chr, Probe, trait, x=b.wave3, y=b.wave4) %>% 
      mutate(comparison="wave 3 vs 4"))
  
  ds_0.001_wide
})

ds_marginal_effects <- do.call(rbind, marginal_effects)

# Signficiant in at least one wave..
# write.table(ds_marginal_effects, file="/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/EWAS/EWAS_0.001_atleast1wave.txt", col.names=T, row.names=F, quote=F)

# Signficiant in all waves..
# write.table(ds_marginal_effects, file="/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/EWAS/EWAS_0.001_allwave.txt", col.names=T, row.names=F, quote=F)


# Compare the difference between sig in at least one wave and sig in all waves
sig_one <- read.table("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/EWAS/EWAS_0.001_atleast1wave.txt", header=T) %>%
  filter(p.wave1 < 0.001 | p.wave3 < 0.001 | p.wave4 < 0.001)

sig_all <- read.table("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/EWAS/EWAS_0.001_allwave.txt", header=T) %>%
  filter(p.wave1 < 0.001 & p.wave3 < 0.001 & p.wave4 < 0.001)

sig_one %>% group_by(trait) %>% tally()
sig_all %>% group_by(trait) %>% tally()

# 0.001
# > sig_one %>% group_by(trait) %>% tally()
#   trait                 n
# 1 bmi               71136
# 2 body_fat          51521
# 3 Glucose            6210
# 4 HDL_cholesterol   56315
# 5 Total_cholesterol 10518
# 6 whr               42698

# > sig_all %>% group_by(trait) %>% tally()
#   trait                 n
# 1 bmi                2660
# 2 body_fat           1231
# 3 Glucose               8
# 4 HDL_cholesterol    1202
# 5 Total_cholesterol    13
# 6 whr                1012

sig_one$trait <- recode_factor(sig_one$trait, 
                               bmi  = "BMI", 
                               body_fat = "Body Fat %", 
                               whr = "Waist to Hip Ratio",
                               Glucose = "Glucose",
                               HDL_cholesterol = "HDL Cholesterol",
                               Total_cholesterol = "Total Cholesterol")

sig_all$trait <- recode_factor(sig_all$trait, 
                               bmi  = "BMI", 
                               body_fat = "Body Fat %", 
                               whr = "Waist to Hip Ratio",
                               Glucose = "Glucose",
                               HDL_cholesterol = "HDL Cholesterol",
                               Total_cholesterol = "Total Cholesterol")


# Clean
sig_one_clean <- rbind(
  sig_one %>% 
    select(Chr, Probe, trait, x=b.wave1, y=b.wave3) %>% 
    mutate(comparison="Set 1 v 2"),
  sig_one %>% 
    select(Chr, Probe, trait, x=b.wave1, y=b.wave4) %>% 
    mutate(comparison="Set 1 v 3"),  
  sig_one %>% 
    select(Chr, Probe, trait, x=b.wave3, y=b.wave4) %>% 
    mutate(comparison="Set 2 v 3"))

sig_all_clean <- rbind(
  sig_all %>% 
    select(Chr, Probe, trait, x=b.wave1, y=b.wave3) %>% 
    mutate(comparison="Set 1 v 2"),
  sig_all %>% 
    select(Chr, Probe, trait, x=b.wave1, y=b.wave4) %>% 
    mutate(comparison="Set 1 v 3"),  
  sig_all %>% 
    select(Chr, Probe, trait, x=b.wave3, y=b.wave4) %>% 
    mutate(comparison="Set 2 v 3"))

# Find correlations
sig_all_clean_pvals <- sig_all_clean %>% group_by(trait, comparison) %>%
  summarize(
    cor=cor.test(x, y)[["estimate"]],
    pvalue=cor.test(x,y)[["p.value"]],
    statistic=cor.test(x,y)[["statistic"]],
    parameter=cor.test(x,y)[["parameter"]]) %>%
  mutate(log.p = 2* pt(statistic,  df = parameter, lower.tail=FALSE, log.p=TRUE),
         edit_pvalue = case_when(pvalue==0 ~ 1e-100, 
                                 TRUE ~ pvalue),
         link = case_when(pvalue==0 ~ "<", 
                          TRUE ~ "=" ),) %>% as.data.frame %>%
  select(trait, comparison, cor, edit_pvalue, link)


# plot

library(viridis)
library(ggpubr)
sig_all_clean %>% 
  left_join(sig_all_clean_pvals, by=c("trait", "comparison")) %>%
  ggplot(aes(x=x, y=y)) +
  geom_point(aes(colour=trait), size=1) +
  # stat_cor(method = "pearson", size=3,   p.accuracy=1e-100) +
  facet_grid(rows=vars(trait), cols=vars(comparison)) +
  scale_color_viridis(discrete = TRUE, option = "D")+
  labs(y="", x="") + 
  theme_linedraw() +
  scale_fill_manual(values = mycolours) +
  labs(y="Probe association coefficient", x="Probe association coefficient")  +
  theme(legend.position="none") +
  geom_text(data = sig_all_clean_pvals, aes(x = -0.6, y = 1.5, 
                                            #                   label = paste0("R = ", round(cor,2), ",  p",link, scientific_10(signif(edit_pvalue, 2)))))
                                            
                                            label=paste("R=", round(cor,2), ", p", link, signif(edit_pvalue, 2))), 
            size=3.5)  

