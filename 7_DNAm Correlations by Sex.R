
###################################
## DNAm Correlation between sexes
##################################

# Bivariate REML

scratch=/Local_Scratch/Alesha
cov=covar_batch_eversmoke_packyears
filename=mvals-norm20k-18413-831733_sd0.02_unrelated_${cov}

mkdir ${scratch}/DNAm/phenotypes/greml.pheno/temp
cd ${scratch}/DNAm/phenotypes/greml.pheno/temp

for trait in bmi whr body_fat Glucose HDL_cholesterol Total_cholesterol height
do

if [ "${trait}" == "bmi" ];              then  var1=3 var2=4; fi
if [ "${trait}" == "whr" ];              then  var1=5 var2=6; fi
if [ "${trait}" == "body_fat" ];         then  var1=7 var2=8; fi
if [ "${trait}" == "Glucose" ];          then  var1=9 var2=10; fi
if [ "${trait}" == "HDL_cholesterol" ];  then  var1=11 var2=12; fi
if [ "${trait}" == "Total_cholesterol" ]; then  var1=13 var2=14; fi
if [ "${trait}" == "height" ];            then  var1=15 var2=16; fi
if [ "${trait}" == "weight" ];            then  var1=17 var2=18; fi
if [ "${trait}" == "waist" ];             then  var1=19 var2=20; fi
if [ "${trait}" == "hips" ];              then  var1=21 var2=22; fi

awk -v "col1=${var1}" -v "col2=${var2}" '{print $1, $2, $col1, $col2}' ${scratch}/DNAm/phenotypes/greml.pheno/dsMvF.txt > temp.${trait}.sex

gcta64 --reml-bivar 1 2 \
--grm ${scratch}/DNAm/ORM/GRM/${filename} \
--reml-bivar-lrt-rg 0  \
--pheno temp.${trait}.sex \
--reml-maxit 100 \
--threads 20 \
--out ${scratch}/DNAm/GREML/sex/${trait}_RankNorm_MvF_0

done

cd ${scratch}/DNAm/GREML/sex
for trait in bmi whr body_fat Glucose HDL_cholesterol Total_cholesterol height weight waist hips
do
awk -v "trait=${trait}"  '{print trait, $1, $2, $3}' ${trait}_RankNorm_MvF.hsq >> trait_MvF.rsq
done



############# 
# Output
#############

library(stringr)
library(readr)
traits <- read_table("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/GREML_sex/trait_MvF.rsq") %>% 
  mutate(Var = case_when(Source %in% c("V(G)/Vp_tr1") ~ "VO/Vp (Males)",
                         Source %in% c("V(G)/Vp_tr2") ~ "VO/Vp (Females)",
                         Source %in% c("rG") ~ "rG",
                         Source == "Pval" ~ "Pval"),
         SE = case_when(SE=="(one-tailed" ~ as.character(NA),
                        TRUE ~ SE)) %>%
  select(-Source) %>%
  filter(!is.na(Var)) 

colnames(traits)[1] <- "Trait"

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

# Recode to numeric
traits$Variance <- traits$Variance %>% as.numeric
traits$SE <- traits$SE %>% as.numeric

traits_wide_table <- traits %>%
  nest(col_value = c(Variance, SE)) %>%
  spread(key = Var, value = col_value) %>%
  unnest(c(`VO/Vp (Males)`, `VO/Vp (Females)`, rG, Pval),  .sep = '_') %>%
  mutate(across(where(is.character), is.numeric)) %>%
  mutate(across(where(is.numeric), round, 4)) %>%
  mutate(`VO/Vp (Males)` = paste0(`VO/Vp (Males)_Variance`, " (", `VO/Vp (Males)_SE`, ")"),
         `VO/Vp (Females)` = paste0(`VO/Vp (Females)_Variance`, " (", `VO/Vp (Females)_SE`, ")"),
         `rG` = paste0(rG_Variance, " (", rG_SE, ")"),
         `Pval` = Pval_Variance)

traits_wide_table %>% select(rG, Pval)

# Figure
figure_ds <- traits %>% filter(Var=="rG") %>% 
  left_join(traits %>% filter(Var=="Pval") %>% select(Trait, Pval=Variance), by=c("Trait")) %>%
  mutate(pval_annot = Pval %>% signif(1))

# options(scipen=999)
# ggarrange(
figure_sex <- figure_ds %>% filter(! Trait %in% c("Hip Circumference", "Waist Circumference", "Weight", "Height")) %>%
  mutate(sex="Male v Female") %>%
  ggplot(aes(y=fct_reorder(Trait, desc(Trait)), x=Variance, xmin=Variance-SE, xmax = Variance+SE, colour=Trait, label=pval_annot)) +
  geom_vline(xintercept = 1, colour="darkgrey") +
  geom_errorbar(position = position_dodge(width = 0.2), width=.1) +
  geom_point(position = position_dodge(width = 0.2), size=3) +
  # geom_text(aes(label=pval_annot),hjust=0.5, vjust=3, size=2,
  #       position = position_dodge(width = 1)) +
  theme_linedraw() +
  scale_colour_manual(values = mycolours, name = " ") +
  labs(y=" ", x="DNAm correlation") +
  theme(legend.position="bottom") +
  theme(legend.position = "none")+
  facet_wrap(~sex) +
  scale_y_discrete(position = "right")



#############################################################
# Pvalue for the likelihood test fixing the correlation at 1
############################################################

# Pheno ~ Probe:sex

# qq-plot of pvalues of interaction effects
path=/Local_Scratch/Alesha/DNAm/phenotypes/age_sex_interactions
scratch=/Local_Scratch/Alesha
cov=covar_batch_eversmoke_packyears
filename=mvals-norm20k-18413-831733_sd0.02_unrelated_${cov}


mkdir ${path}/DNAm
mkdir ${path}/pvalue_output

# Make by chr then merge together
for chr in {1..22}
do

osca_Linux  --befile ${scratch}/DNAm/OSCA/${filename} \
--chr ${chr} \
--make-efile --out ${path}/DNAm/DNAm_unrelated_chr${chr}.txt

done


# 1. extract DNAm 
# 2. Read into R along with ds_pheno..


library(tidyverse)
ds <- read.table("/Local_Scratch/Alesha/DNAm/phenotypes/bodyfat_covariates_agesex_ranknorm_outliers4sd_related.txt", header=T)

DNAm <- data.table::fread(paste0("/Local_Scratch/Alesha/DNAm/phenotypes/age_sex_interactions/DNAm/DNAm_unrelated_chr",CHR,".txt"), header=T) 

ds <- ds %>% mutate(sex=as.factor(sex))

DNAm_cov <- ds %>% select(FID=FID_DNAm, IID=IID_DNAm, trait=bmi_RankNorm, sex) %>%
  inner_join(DNAm, by=c("FID", "IID"))

# 3. Probelist to loop over..
probelist <- colnames(DNAm)[-c(1,2)] # 16874 probes in chr22
length(probelist)

# 4. LM

# Set up outut file
pval.matrix <- probelist %>% as.data.frame
colnames(pval.matrix) <- c("probeID")

pval.matrix$probe_Estimate <- NA
pval.matrix$probe_StdError <- NA
pval.matrix$probe_tvalue <- NA
pval.matrix$probe_P <- NA

pval.matrix$sex_Estimate <- NA
pval.matrix$sex_StdError <- NA
pval.matrix$sex_tvalue <- NA
pval.matrix$sex_P <- NA

pval.matrix$interaction_Estimate <- NA
pval.matrix$interaction_StdError <- NA
pval.matrix$interaction_tvalue <- NA
pval.matrix$interaction_P <- NA


# For each probe, create normalised residuals
for (i in 1:length(probelist)) {
  
  print(paste0(CHR, "   ", i))
  
  # Select probevalue
  probelist[i]
  this_probe <- probelist[i]
  this_ds <- DNAm_cov %>% select(FID, IID, trait, sex, DNAm=probelist[i])
  
  lm <- lm(trait ~ DNAm + sex + DNAm:sex, data=this_ds)
  
  # Add to normalised residual matrix
  
  # Probe
  pval.matrix[pval.matrix$probeID==this_probe,]$probe_Estimate <- summary(lm)$coefficients[2,1]
  pval.matrix[pval.matrix$probeID==this_probe,]$probe_StdError <- summary(lm)$coefficients[2,2]
  pval.matrix[pval.matrix$probeID==this_probe,]$probe_tvalue <- summary(lm)$coefficients[2,3]
  pval.matrix[pval.matrix$probeID==this_probe,]$probe_P <- summary(lm)$coefficients[2,4]
  
  # Sex
  pval.matrix[pval.matrix$probeID==this_probe,]$sex_Estimate <- summary(lm)$coefficients[3,1]
  pval.matrix[pval.matrix$probeID==this_probe,]$sex_StdError <- summary(lm)$coefficients[3,2]
  pval.matrix[pval.matrix$probeID==this_probe,]$sex_tvalue <- summary(lm)$coefficients[3,3]
  pval.matrix[pval.matrix$probeID==this_probe,]$sex_P <- summary(lm)$coefficients[3,4]
  
  # Interaction
  pval.matrix[pval.matrix$probeID==this_probe,]$interaction_Estimate <- summary(lm)$coefficients[4,1]
  pval.matrix[pval.matrix$probeID==this_probe,]$interaction_StdError <- summary(lm)$coefficients[4,2]
  pval.matrix[pval.matrix$probeID==this_probe,]$interaction_tvalue <- summary(lm)$coefficients[4,3]
  pval.matrix[pval.matrix$probeID==this_probe,]$interaction_P <- summary(lm)$coefficients[4,4]
  
  # Re-set all values at the end of each cycle
  this_probe <- NULL
  this_ds <- NULL
  lm <- NULL
  probe_Estimate = sex_Estimate = interaction_Estimate <- NULL
  probe_StdError = sex_StdError = interaction_StdError <- NULL
  probe_tvalue = sex_tvalue = interaction_tvalue <- NULL
  probe_P = sex_P = interaction_P <- NULL
  
}

write.table(pval.matrix, file=paste0("/Local_Scratch/Alesha/DNAm/phenotypes/age_sex_interactions/pvalue_output/interaction_pval_chr",CHR,".txt"), col.name=T, row.name=F, quote=F)




# }


# Combine all chr
DNAm <- data.table::fread(paste0("/Local_Scratch/Alesha/DNAm/phenotypes/age_sex_interactions/pvalue_output/interaction_pval_chr1.txt"), header=T)
DNAm$chr <- 1
for (i in 2:22) {
  DNAm_chr <- data.table::fread(paste0("/Local_Scratch/Alesha/DNAm/phenotypes/age_sex_interactions/pvalue_output/interaction_pval_chr",i,".txt"), header=T)
  DNAm_chr$chr <- i
  DNAm <- rbind(DNAm, DNAm_chr)
}
write.table(DNAm, file=paste0("/Local_Scratch/Alesha/DNAm/phenotypes/age_sex_interactions/interaction_pval.txt"), col.name=T, row.name=F, quote=F)

scp -P 8017 v1ahatt4@localhost:/Local_Scratch/Alesha/DNAm/phenotypes/age_sex_interactions/interaction_pval.txt  /Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/GREML_sex

# Reference is Male!






###########################
# Idenitify signficant probes
###########################



library(qqman)
DNAm <- data.table::fread(paste0("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/GREML_sex/interaction_pval.txt"), header=T)
qq(DNAm$probe_P)

DNAm %>% filter(probeID=="cg26651978")

pthreshold = 0.05/nrow(DNAm)

DNAm %>% filter(interaction_P < pthreshold)  # 8..
DNAm %>% filter(interaction_P < pthreshold) %>% select(probeID, chr, probe_Estimate, probe_StdError, probe_P, sex_Estimate, sex_StdError, sex_P, interaction_Estimate, interaction_StdError, interaction_P)


options(scipen=NA)


# Extract those DNAm values
path=/Local_Scratch/Alesha/DNAm/phenotypes/age_sex_interactions
scratch=/Local_Scratch/Alesha
cov=covar_batch_eversmoke_packyears
filename=mvals-norm20k-18413-831733_sd0.02_unrelated_${cov}

echo cg12269535 > ${path}/DNAm/sig_probes.txt
echo cg12054453 >> ${path}/DNAm/sig_probes.txt
echo cg18181703 >> ${path}/DNAm/sig_probes.txt
echo cg11047325 >> ${path}/DNAm/sig_probes.txt
echo cg19748455 >> ${path}/DNAm/sig_probes.txt
echo cg16936953 >> ${path}/DNAm/sig_probes.txt
echo cg00840791 >> ${path}/DNAm/sig_probes.txt
echo cg09349128 >> ${path}/DNAm/sig_probes.txt
echo cg15125798 >> ${path}/DNAm/sig_probes.txt
echo cg26651978 >> ${path}/DNAm/sig_probes.txt

osca_Linux  --befile ${scratch}/DNAm/OSCA/${filename} \
--extract-probe ${path}/DNAm/sig_probes.txt \
--make-efile --out ${path}/DNAm/DNAm_unrelated_sig_probes.txt


library(tidyverse)
ds <- read.table("/Local_Scratch/Alesha/DNAm/phenotypes/bodyfat_covariates_agesex_ranknorm_outliers4sd_related.txt", header=T)
age <- read.table("/Local_Scratch/Alesha/DNAm/GS_20K_samplesheet_related.txt", header=T)
DNAm <- data.table::fread(paste0("/Local_Scratch/Alesha/DNAm/phenotypes/age_sex_interactions/DNAm/DNAm_unrelated_sig_probes.txt"), header=T) 
ds <- ds %>% mutate(sex=as.factor(sex))
DNAm_cov <- ds %>% select(FID=FID_DNAm, IID=IID_DNAm, trait=bmi_RankNorm, sex) %>%
  inner_join(DNAm, by=c("FID", "IID")) %>%
  inner_join(age %>% select(age, Sample_Sentrix_ID), by=c("FID"="Sample_Sentrix_ID"))

write.table(DNAm_cov, file="/Local_Scratch/Alesha/DNAm/phenotypes/age_sex_interactions/DNAm/DNAm_unrelated_sig_probes_sex.txt", col.names=T, row.names=F, quote=F)



library(sjPlot)
library(sjmisc)
library(ggplot2)
library(tidyverse)

sig_probes <- read.table("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/GREML_sex/DNAm_unrelated_sig_probes_sex.txt", header=T) %>% 
  dplyr::rename(BMI = trait) %>% 
  na.omit() 


# Linear model in Males and Female
probelist <- colnames(sig_probes)[5:12] 
length(probelist)

# Set up outut file
pval.matrix <- probelist %>% as.data.frame
colnames(pval.matrix) <- c("probeID")
pval.matrix <- pval.matrix %>% mutate(M_Estimate=NA, M_StdError=NA, M_tvalue=NA, M_P=NA,
                                      F_Estimate=NA, F_StdError=NA, F_tvalue=NA, F_P=NA)

# For each probe, create normalised residuals
for (i in 1:length(probelist)) {
  
  # Select probevalue
  probelist[i]
  this_probe <- probelist[i]
  this_ds <- sig_probes %>% select(FID, IID, BMI, sex, DNAm=probelist[i])
  
  lm_M <- lm(BMI ~ DNAm, data=this_ds %>% filter(sex=="M"))
  lm_F <- lm(BMI ~ DNAm, data=this_ds %>% filter(sex=="F"))
  
  # Add to normalised residual matrix
  # Male
  pval.matrix[pval.matrix$probeID==this_probe,]$M_Estimate <- summary(lm_M)$coefficients[2,1]
  pval.matrix[pval.matrix$probeID==this_probe,]$M_StdError <- summary(lm_M)$coefficients[2,2]
  pval.matrix[pval.matrix$probeID==this_probe,]$M_tvalue <- summary(lm_M)$coefficients[2,3]
  pval.matrix[pval.matrix$probeID==this_probe,]$M_P <- summary(lm_M)$coefficients[2,4]
  
  # Female
  pval.matrix[pval.matrix$probeID==this_probe,]$F_Estimate <- summary(lm_F)$coefficients[2,1]
  pval.matrix[pval.matrix$probeID==this_probe,]$F_StdError <- summary(lm_F)$coefficients[2,2]
  pval.matrix[pval.matrix$probeID==this_probe,]$F_tvalue <- summary(lm_F)$coefficients[2,3]
  pval.matrix[pval.matrix$probeID==this_probe,]$F_P <- summary(lm_F)$coefficients[2,4]
  
  
  # Re-set all values at the end of each cycle
  this_probe <- NULL
  this_ds <- NULL
  lm <- NULL
}

pval.matrix #Estaimtes for Table2
# Done

plot_model(probe1, type = "pred", terms = c("cg12054453", "sex"))
plot_model(probe1, type = "int")

sig_probes_long <- 
  sig_probes %>% 
  gather(probeID, DNAm, cg12054453:cg09349128, factor_key=TRUE) %>%
  mutate(CHR = case_when(
    probeID %in% c("cg12269535") ~ 6,
    probeID %in% c("cg12054453", "cg18181703", "cg11047325", "cg19748455", "cg16936953") ~ 17,
    probeID %in% c("cg00840791") ~ 19,
    probeID %in% c("cg09349128") ~ 22),
    probe_POS = case_when(
      probeID=="cg12269535" ~ 43142014,
      probeID=="cg12054453" ~ 57915717, 
      probeID=="cg18181703" ~ 76354621, 
      probeID=="cg19748455" ~ 76274856, 
      probeID=="cg16936953" ~ 57915665, 
      probeID=="cg11047325" ~ 76354934, 
      probeID=="cg00840791" ~ 16453259,
      probeID=="cg09349128" ~ 50327986),
    order = CHR + probe_POS/max(probe_POS))

library(RNOmni)
sig_probes_long %>% group_by(probeID) %>% mutate(scale_DNAm = RankNorm(DNAm)/4) %>% ungroup %>% 
  mutate(agegrp = case_when(age<50 ~ 1,
                            TRUE ~ 2)) %>%
  ggplot(aes(x = DNAm, y = BMI, color = sex)) +
  # geom_point(alpha=0.1, colour="grey") +
  facet_wrap(~ fct_reorder(probeID, desc(-order)), nrow=2) +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c(mycolours[2], mycolours[1])) +
  theme_linedraw() +
  scale_x_continuous(limits=c(-1, 1), breaks=c(-1, 0, 1), name="DNAm") +
  guides(color=guide_legend(override.aes=list(fill=NA)))



# Checking co-methylation
sig_probes %>% 
  select(-FID, -IID, -trait, -sex) %>%
  as.matrix %>%
  cor %>%
  as.data.frame %>%
  rownames_to_column(var = 'var1') %>%
  gather(var2, value, -var1) %>% 
  filter(var1!=var)

cor(sig_probes$cg18181703, sig_probes$cg11047325)

cor(sig_probes$cg16936953, sig_probes$cg12054453)
cor(sig_probes$cg12054453, sig_probes$cg12269535)
cor(sig_probes$cg12054453, sig_probes$cg12269535)


sig_probes %>% 
  select(BMI, sex, cg18181703) %>% 
  ggplot(aes(x=BMI, y=cg18181703, colour=sex)) + 
  geom_smooth(method = "lm") 
sig_probes$sex2 <- ifelse(sig_probes$sex=="M" , 1, 0)

probes <- list("cg12269535", "cg16936953", "cg12054453", "cg19748455", "cg18181703", "cg11047325", "cg00840791", "cg09349128")
forms <- paste('BMI ~ ', probes)
lapply(forms, function(x) summary(lm(x, data =sig_probes %>% filter(sex=="F")))$coeff) 
lapply(forms, function(x) summary(lm(x, data =sig_probes %>% filter(sex=="M")))$coeff ) 
forms <- paste('BMI ~ ', probes,"*sex")
lapply(forms, function(x) summary(lm(x , data =sig_probes))$coeff ) 


# Investigating probes previously reported
# Wang et al. one CpG showed a statistically significant DNAm × sex interaction effect on BMI status transition across adolescence in the IoW cohort; at cg15125798, boys with higher DNAm were more likely to be in the transition group from normal weight to overweight or obese, but for girls, the association was in the opposite direction (p value for interaction effect: 2.3 × 10–3, Table 3).
# Mendelson et al. a significant sex interaction was demonstrated in the discovery cohorts for one unannotated CpG (cg26651978 on Chromosome 17q25.3; <3 kbp from the 3′ end of LGALS3BP [lectin galactoside-binding soluble 3-binding protein]), 


DNAm <- data.table::fread(paste0("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/GREML_sex/interaction_pval.txt"), header=T)
DNAm %>% filter(probeID %in% c("cg15125798", "cg26651978")) %>% select(probeID, interaction_P)
sig_probes <- read.table("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/GREML_sex/DNAm_unrelated_sig_probes_sex.txt", header=T) %>% 
  dplyr::rename(BMI = trait) %>% 
  na.omit()

library(RNOmni)
pvalue <- DNAm %>% filter(probeID %in% c("cg15125798", "cg26651978")) %>% select(probeID, interaction_P) %>% mutate(sex=)

rbind(sig_probes %>% select(BMI, Sex=sex, DNAm=cg15125798) %>% mutate(probeID="cg15125798"),
      sig_probes %>% select(BMI, Sex=sex, DNAm=cg26651978) %>% mutate(probeID="cg26651978")) %>% 
  group_by(probeID) %>% mutate(scale_DNAm = RankNorm(DNAm)/4) %>% ungroup %>% 
  ggplot(aes(x = DNAm, y = BMI, color = Sex)) +
  facet_wrap(~ probeID, ncol=2) +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c(mycolours[2], mycolours[1])) +
  theme_linedraw() +
  scale_x_continuous(limits=c(-1, 1), breaks=c(-1, 0, 1), name="DNAm") +
  labs(fill = "Sex") +
  geom_text(
    data    = pvalue,
    mapping = aes(x = -Inf, y = -Inf, label = paste0("Pval  = ",interaction_P %>% round(2))),
    colour="black", 
    hjust   = -0.1,
    vjust   = -1)
