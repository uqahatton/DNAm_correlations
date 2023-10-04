
###############################
## Correlation between traits
###############################

# Bivariate REML

# bmi_RankNorm $3
# whr_RankNorm $4
# body_fat_RankNorm $5
# Glucose_RankNorm $6
# HDL_cholesterol_RankNorm $7
# Total_cholesterol_RankNorm $8

scratch=/Local_Scratch/Alesha
cov=covar_batch_eversmoke_packyears
filename=mvals-norm20k-18413-831733_sd0.02_unrelated_${cov}
mkdir ${scratch}/DNAm/GREML/traits/${cov}


cd ${scratch}/DNAm/phenotypes/greml.pheno/temp

# for trait1 in bmi whr body_fat Glucose HDL_cholesterol
for trait1 in bmi whr Total_cholesterol
do

if [ "${trait1}" == "bmi" ];              then  var1=3; fi
if [ "${trait1}" == "whr" ];              then  var1=4; fi
if [ "${trait1}" == "body_fat" ];         then  var1=5; fi
if [ "${trait1}" == "Glucose" ];          then  var1=6; fi
if [ "${trait1}" == "HDL_cholesterol" ];  then  var1=7; fi
if [ "${trait1}" == "Total_cholesterol" ]; then  var1=8; fi


for trait2 in height

do

if [ "${trait2}" == "bmi" ];              then  var2=3; fi
if [ "${trait2}" == "whr" ];              then  var2=4; fi
if [ "${trait2}" == "body_fat" ];         then  var2=5; fi
if [ "${trait2}" == "Glucose" ];          then  var2=6; fi
if [ "${trait2}" == "HDL_cholesterol" ];  then  var2=7; fi
if [ "${trait2}" == "Total_cholesterol" ]; then  var2=8; fi
if [ "${trait2}" == "height" ]; then  var2=9; fi


trait_pair=${trait1}_${trait2}
awk -v "col1=${var1}" -v "col2=${var2}" '{print $1, $2, $col1, $col2}' /Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds_pheno.txt > temp3.${trait_pair}

gcta64 --reml-bivar 1 2 \
--grm ${scratch}/DNAm/ORM/GRM/${filename} \
--reml-bivar-lrt-rg 0 \
--thread-num 10 \
--pheno temp3.${trait_pair} \
--out ${scratch}/DNAm/GREML/traits/${cov}/${trait_pair}_0

done
done



cd ${scratch}/DNAm/GREML/traits/${cov}
for trait1 in bmi whr body_fat Glucose HDL_cholesterol Total_cholesterol height weight waist
do
for trait2 in bmi whr body_fat Glucose HDL_cholesterol Total_cholesterol height weight waist hips
do
awk -v "trait1=${trait1}" -v "trait2=${trait2}"  '{print trait1, trait2, $1, $2, $3}' ${trait1}_${trait2}.hsq >> trait_correlation.hsq
done
done


#####################
# Output
####################
  
library(stringr)
library(readr)
library(tidyverse)
library(forcats)

traits <- read_table("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/GREML_traits/covar_batch_eversmoke_packyears/trait_correlation.hsq") %>% 
  filter(!is.na(SE) & SE!="SE" & Variance!="of") %>%
  mutate(Var = case_when(Source %in% c("V(G)/Vp_tr1") ~ "V(O)/Vp1",
                         Source %in% c("V(G)/Vp_tr2") ~ "V(O)/Vp2",
                         Source == "rG" ~ "rG",
                         Source == "Pval" ~ "Pval")) %>%
  select(-Source) %>% filter(!is.na(Var)) 

colnames(traits)[1] <- "Trait1"
colnames(traits)[2] <- "Trait2"

# Recode to numeric
traits$Variance <- traits$Variance %>% as.numeric
traits$SE <- traits$SE %>% as.numeric

traits <- traits %>% mutate(Trait1a = ifelse(Trait1=="whr" & Trait2=="body_fat", "body_fat", Trait1),
                            Trait2a = ifelse(Trait1=="whr" & Trait2=="body_fat", "whr", Trait2)) %>%
  select(-Trait1, -Trait2, Variance , SE, Var, Trait1=Trait1a, Trait2= Trait2a ) 

traits_wide <- traits %>%
  nest(col_value = c(Variance, SE)) %>%
  spread(key = Var, value = col_value) %>%
  unnest(c(`V(O)/Vp1`, `V(O)/Vp2`, `rG`, `Pval`),  .sep = '_')


traits_wide$Trait1 <- recode_factor( traits_wide$Trait1, 
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
traits_wide$Trait2 <- recode_factor( traits_wide$Trait2, 
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

plot_ds <- traits_wide %>% filter(!Trait1 %in% c("Waist Circumference", "Weight") & !(Trait2 %in% c("Hip Circumference", "Weight", "Waist Circumference")) ) %>%
  mutate(Trait1 = fct_reorder(Trait1, desc(Trait1)),
         label = paste0(rG_Variance %>% signif(digits = 2), " (", rG_SE%>% signif(digits = 1), ")")) 

# label = paste0(rG_Variance %>% signif(digits = 2))) %>% select(Trait1, Trait2, rG_Variance, label)
plot_ds$Trait1 <- plot_ds$Trait1 %>% as.character
plot_ds$Trait2 <- plot_ds$Trait2 %>% as.character
plot_ds$rG_Variance <- plot_ds$rG_Variance %>% as.character

plot_ds <- plot_ds %>% select(Trait1, Trait2, rG_Variance, label)

plot_ds[nrow(plot_ds) + 1,] <- c("BMI", "BMI", NA, NA) %>% as.list
plot_ds[nrow(plot_ds) + 1,] <- c("Waist to Hip Ratio", "Waist to Hip Ratio", NA, NA) %>% as.list
plot_ds[nrow(plot_ds) + 1,] <- c("Body Fat %", "Body Fat %", NA, NA) %>% as.list
plot_ds[nrow(plot_ds) + 1,] <- c("Glucose", "Glucose",  NA, NA) %>% as.list
plot_ds[nrow(plot_ds) + 1,] <- c("HDL Cholesterol", "HDL Cholesterol",  NA, NA) %>% as.list
plot_ds[nrow(plot_ds) + 1,] <- c("Total Cholesterol", "Total Cholesterol", NA, NA) %>% as.list
# plot_ds[nrow(plot_ds) + 1,] <- c("Height", "Height",  NA, NA) %>% as.list

plot_ds$Trait1 <- plot_ds$Trait1 %>% as.factor
plot_ds$Trait2 <- plot_ds$Trait2 %>% as.factor



plot_ds$Trait1 <- recode_factor( plot_ds$Trait1, 
                                 "BMI"  = "BMI", 
                                 "Body Fat %" = "Body Fat %", 
                                 "Waist to Hip Ratio" = "Waist to Hip Ratio",
                                 "Glucose" = "Glucose",
                                 "HDL Cholesterol" = "HDL Cholesterol",
                                 "Total Cholesterol" = "Total Cholesterol",
                                 "Height" = "Height")
plot_ds$Trait2 <- recode_factor( plot_ds$Trait2, 
                                 "BMI"  = "BMI", 
                                 "Body Fat %" = "Body Fat %", 
                                 "Waist to Hip Ratio" = "Waist to Hip Ratio",
                                 "Glucose" = "Glucose",
                                 "HDL Cholesterol" = "HDL Cholesterol",
                                 "Total Cholesterol" = "Total Cholesterol",
                                 "Height" = "Height")
plot_ds$rG_Variance <- plot_ds$rG_Variance %>% as.numeric


# plot_ds1 <- plot_ds %>% left_join(traits_wide %>% select(Trait1, Trait2, Pval_Variance), by=c("Trait1", "Trait2"))
r_DNAm_plot <- plot_ds
rg_DNAm <- ggplot(plot_ds
                  %>% filter(Trait1!="Height" & Trait2!="Height" & Trait1!="Total Cholesterol" & Trait2!="BMI") ,
                  # %>% filter(Trait1!="Waist to Hip Ratio"& Trait2!="Body Fat %") ,
                  aes(x = Trait2, y = Trait1, fill = rG_Variance)) +
  # aes(x = fct_reorder(Trait1, desc(Trait1)), y = Trait2, fill = rG_Variance)) +
  geom_raster() +
  # geom_point(aes(size=Pval_Variance), shape=15) +
  geom_text(aes(label = label), size=3) +
  # scale_fill_distiller(palette = "Spectral", na.value = 'white',  limits=c(-0.7,1)) +
  theme_minimal() +
  theme(panel.grid = element_blank())+ 
  labs(fill = "rG (DNAm)", y="", x="") +
  scale_y_discrete(limits=rev) +
  scale_fill_gradient2(low = mycolours[1], mid = "beige", high = mycolours[3], na.value = 'white',  limits=c(-0.7,1), guide = "colourbar")  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(legend.title = element_text( size=12), legend.text=element_text(size=12))








##########################################
## Phenotypic correlation between traits
##########################################

library(tidyverse)
out="/Local_Scratch/Alesha/DNAm/phenotypes/descriptive"
ds_pheno <- read.table("/Local_Scratch/Alesha/DNAm/phenotypes/bodyfat_covariates_agesex_ranknorm_outliers4sd.txt", header=T)

traits <- c("bmi", "whr", "body_fat", "Glucose", "HDL_cholesterol", "Total_cholesterol", "height", "weight", "waist", "hips")

corr_with_RankNorm <- data.frame("Trait" = traits, 
                                 "Correlation with RankNorm" = c(
                                   cor(ds_pheno$bmi, ds_pheno$bmi_RankNorm, use = "complete.obs"),
                                   cor(ds_pheno$whr, ds_pheno$whr_RankNorm, use = "complete.obs"),
                                   cor(ds_pheno$body_fat, ds_pheno$body_fat_RankNorm, use = "complete.obs"),
                                   cor(ds_pheno$Glucose, ds_pheno$Glucose_RankNorm, use = "complete.obs"), 
                                   cor(ds_pheno$HDL_cholesterol, ds_pheno$HDL_cholesterol_RankNorm, use = "complete.obs"),
                                   cor(ds_pheno$Total_cholesterol, ds_pheno$Total_cholesterol_RankNorm, use = "complete.obs"),
                                   cor(ds_pheno$height, ds_pheno$height_RankNorm, use = "complete.obs"),
                                   cor(ds_pheno$weight, ds_pheno$weight_RankNorm, use = "complete.obs"),
                                   cor(ds_pheno$waist, ds_pheno$waist_RankNorm, use = "complete.obs"),
                                   cor(ds_pheno$hips, ds_pheno$hips_RankNorm, use = "complete.obs")))

write.table(corr_with_RankNorm, file=paste0(out,"/corr_with_RankNorm.txt"), row.names=F, col.names=T, quote=F)

# Between traits
library("Matrix")
library("psych")
pheno_cor <- ds_pheno %>% select(ends_with("RankNorm")) %>%
  rename_with(~str_remove(., '_RankNorm')) %>% corr.test

pheno_cor_tab <- pheno_cor$r %>% unlist %>% as.data.frame 
pheno_cor_se <- pheno_cor$se %>% unlist %>% as.data.frame 


write.table(pheno_cor_tab, file=paste0(out,"/pheno_correlations.txt"), col.names = T, row.names = T, quote = F)
write.table(pheno_cor_se, file=paste0(out,"/pheno_correlations_se.txt"), col.names = T, row.names = T, quote = F)



corr_with_RankNorm <- read.table("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/phenotypes/corr_with_RankNorm.txt", header=T) %>% mutate(Correlation.with.RankNorm=Correlation.with.RankNorm %>% signif(2))

corr_with_RankNorm$Trait <- recode_factor( corr_with_RankNorm$Trait, 
                                           bmi  = "BMI", 
                                           whr = "Waist to Hip Ratio",
                                           body_fat = "Body Fat %", 
                                           Glucose = "Glucose",
                                           HDL_cholesterol = "HDL Cholesterol",
                                           Total_cholesterol = "Total Cholesterol",
                                           height = "Height",
                                           weight = "Weight",
                                           waist = "Waist Circumference",
                                           hips = "Hip Circumference")

corr_with_RankNorm %>% select(Trait, Correlation.with.RankNorm) %>% kable(caption=" Correlation between each trait and rank normalised residuals")


pheno_correlations <- read.table("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/phenotypes/pheno_correlations.txt", header=T) 
names <- recode_factor(rownames(pheno_correlations) , 
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
colnames(pheno_correlations) = rownames(pheno_correlations) = names
pheno_correlations <- pheno_correlations[c(1, 3, 2, 4, 5), c(3, 2, 4, 5, 6)]

pheno_correlations_mat <- pheno_correlations %>% as.matrix %>% signif(2)
pheno_correlations_mat[lower.tri(pheno_correlations_mat)] <- NA

options(knitr.kable.NA = '')
pheno_correlations_mat %>% kable(caption = "Phenotypic correlation between traits")


pheno_correlations_se <- read.table("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/phenotypes/pheno_correlations_se.txt", header=T)
colnames(pheno_correlations_se) = rownames(pheno_correlations_se) = names
pheno_correlations_se_mat <- pheno_correlations_se[c(1, 3, 2, 4, 5), c(3, 2, 4, 5, 6)]
pheno_correlations_se_mat <- pheno_correlations_se_mat %>% as.matrix %>% signif(2)
pheno_correlations_se_mat[lower.tri(pheno_correlations_se_mat)] <- NA



library(reshape2)
pheno <- reshape2::melt(pheno_correlations_mat) %>% mutate(Trait1=Var2, Trait2=Var1)
pheno_se <- reshape2::melt(pheno_correlations_se_mat) %>% mutate(Trait1=Var2, Trait2=Var1)

pheno <- pheno %>% left_join(pheno_se %>% select(Trait1, Trait2, se=value), by=c("Trait1", "Trait2"))


plot_ds <- pheno %>% filter(
  # !Trait1 %in% c("Waist Circumference", "Weight", "Hip Circumference") & !(Trait2 %in% c("Hip Circumference", "Weight", "Waist Circumference")) & Trait1!=Trait2 & 
  !is.na(value) ) %>%
  mutate(Trait1 = fct_reorder(Trait1, desc(Trait1)),
         label = paste0(value %>% signif(digits = 2), " (", se %>% signif(digits=1), ")"))  %>% 
  select(Trait1, Trait2, value, label)

plot_ds$Trait1 <- plot_ds$Trait1 %>% as.character
plot_ds$Trait2 <- plot_ds$Trait2 %>% as.character
plot_ds$rG_Variance <- plot_ds$value %>% as.character



plot_ds[nrow(plot_ds) + 1,] <- c("BMI", "BMI", NA, NA, NA) %>% as.list
plot_ds[nrow(plot_ds) + 1,] <- c("Waist to Hip Ratio", "Waist to Hip Ratio", NA, NA, NA) %>% as.list
plot_ds[nrow(plot_ds) + 1,] <- c("Body Fat %", "Body Fat %", NA, NA, NA) %>% as.list
plot_ds[nrow(plot_ds) + 1,] <- c("Glucose", "Glucose",  NA, NA, NA) %>% as.list
plot_ds[nrow(plot_ds) + 1,] <- c("HDL Cholesterol", "HDL Cholesterol",  NA, NA, NA) %>% as.list
plot_ds[nrow(plot_ds) + 1,] <- c("Total Cholesterol", "Total Cholesterol", NA, NA, NA) %>% as.list
plot_ds[nrow(plot_ds) + 1,] <- c("Height", "Height",  NA, NA, NA) %>% as.list

plot_ds$Trait1 <- plot_ds$Trait1 %>% as.factor
plot_ds$Trait2 <- plot_ds$Trait2 %>% as.factor

# plot_ds <- plot_ds %>% mutate(Trait1a = ifelse(Trait1=="Waist to Hip Ratio" & Trait2=="Body Fat %" , "Body Fat %", as.character(Trait1)),
#                    Trait2a = ifelse(Trait1=="Waist to Hip Ratio" & Trait2=="Body Fat %" , "Waist to Hip Ratio",  as.character(Trait2)))
#   
plot_ds$Trait1 <- recode_factor( plot_ds$Trait1, 
                                 "BMI"  = "BMI", 
                                 "Body Fat %" = "Body Fat %", 
                                 "Waist to Hip Ratio" = "Waist to Hip Ratio",
                                 "Glucose" = "Glucose",
                                 "HDL Cholesterol" = "HDL Cholesterol",
                                 "Total Cholesterol" = "Total Cholesterol",
                                 "Height" = "Height")
plot_ds$Trait2 <- recode_factor( plot_ds$Trait2, 
                                 "BMI"  = "BMI", 
                                 "Body Fat %" = "Body Fat %", 
                                 "Waist to Hip Ratio" = "Waist to Hip Ratio",
                                 "Glucose" = "Glucose",
                                 "HDL Cholesterol" = "HDL Cholesterol",
                                 "Total Cholesterol" = "Total Cholesterol",
                                 "Height" = "Height")
plot_ds$value <- plot_ds$value %>% as.numeric


r_pheno_plot <- plot_ds  %>% filter(Trait1!="Height" & Trait2!="Height" & Trait1!="BMI" & Trait2!="Total Cholesterol") %>%
  rename(Trait1=Trait2, Trait2=Trait1) 



rg_pheno <- ggplot(plot_ds %>% filter(Trait1!="Height" & Trait2!="Height" & Trait1!="BMI" & Trait2!="Total Cholesterol"),
                   aes(x = Trait1, y = Trait2, fill = value)) +
  # aes(x = fct_reorder(Trait2, desc(Trait2)), y = Trait1, fill = value)) +
  geom_raster() +
  geom_text(aes(label = label), size=3) +
  scale_fill_gradient2(low = mycolours[1], mid = "beige", high = mycolours[3], na.value = 'white',  
                       limits=c(-0.7,1), guide = "colourbar")  +
  theme_minimal() +
  theme(panel.grid = element_blank())+ 
  labs(fill = "Correlation", y="", x="") +
  scale_y_discrete(limits=rev) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(legend.title = element_text( size=12), legend.text=element_text(size=12)) 





##########################################
## Genetic correlation between traits
##########################################

# cp GRM to dir
scratch=/Local_Scratch/Alesha/DNAm/GREML/traits/
  cp ${scratch}/GRM_ORM/GS_GWAS_unrelated.grm.id ${scratch}/GRM/GS_GWAS_unrelated.grm.id
cp ${scratch}/GRM_ORM/GS_GWAS_unrelated.grm.N.bin ${scratch}/GRM/GS_GWAS_unrelated.grm.N.bin
cp ${scratch}/GRM_ORM/GS_GWAS_unrelated.grm.bin ${scratch}/GRM/GS_GWAS_unrelated.grm.bin

# bmi_RankNorm $3
# whr_RankNorm $4
# body_fat_RankNorm $5
# Glucose_RankNorm $6
# HDL_cholesterol_RankNorm $7
# Total_cholesterol_RankNorm $8

scratch=/Local_Scratch/Alesha
cd ${scratch}/DNAm/phenotypes/greml.pheno/temp


#for trait1 in bmi whr body_fat Glucose HDL_cholesterol
#do

trait1=bmi
if [ "${trait1}" == "bmi" ];              then  var1=3; fi
if [ "${trait1}" == "whr" ];              then  var1=4; fi
if [ "${trait1}" == "body_fat" ];         then  var1=5; fi
if [ "${trait1}" == "Glucose" ];          then  var1=6; fi
if [ "${trait1}" == "HDL_cholesterol" ];  then  var1=7; fi
if [ "${trait1}" == "Total_cholesterol" ]; then  var1=8; fi
if [ "${trait1}" == "height" ]; then  var1=9; fi

#for trait2 in  height
for trait2 in  body_fat 
do

if [ "${trait2}" == "bmi" ];              then  var2=3; fi
if [ "${trait2}" == "whr" ];              then  var2=4; fi
if [ "${trait2}" == "body_fat" ];         then  var2=5; fi
if [ "${trait2}" == "Glucose" ];          then  var2=6; fi
if [ "${trait2}" == "HDL_cholesterol" ];  then  var2=7; fi
if [ "${trait2}" == "Total_cholesterol" ]; then  var2=8; fi
if [ "${trait2}" == "height" ]; then  var2=9; fi


trait_pair=${trait1}_${trait2}
awk -v "col1=${var1}" -v "col2=${var2}" '{print $1, $2, $col1, $col2}' /Local_Scratch/Alesha/DNAm/phenotypes/greml.pheno/ds_pheno.txt > temp.${trait_pair}_grm

# --grm ${scratch}/DNAm/GREML/traits/GRM/GS_GWAS_unrelated \ -->

gcta64 --reml-bivar 1 2 \
--grm /Local_Scratch/Alesha/DNAm/ORM/GRM_GS_unrelated/imputed_HRS/GS20K_HRC.r1-1_nomono_I4_cpra_unrelated_hm3 \
--reml-bivar-lrt-rg 0 \
--pheno temp.${trait_pair}_grm \
--reml-bivar-no-constrain \
--thread-num 10 \
--out ${scratch}/DNAm/GREML/traits/GRM/${trait_pair}_imputed_hm3_0

done
#done




cd ${scratch}/DNAm/GREML/traits/GRM
for trait1 in bmi whr body_fat Glucose HDL_cholesterol Total_cholesterol height
do
for trait2 in  whr body_fat Glucose HDL_cholesterol Total_cholesterol height 
do
awk -v "trait1=${trait1}" -v "trait2=${trait2}"  '{print trait1, trait2, $1, $2, $3}' ${trait1}_${trait2}_imputed_hm3.hsq >> trait_correlation_GRM_imputed_hm3.hsq
done
done




###########
# Output
###########

library(stringr)
library(readr)
library(tidyverse)

traits <- read_table("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/GREML_traits/covar_batch_eversmoke_packyears/trait_correlation_GRM_imputed_hm3.hsq") %>% 
  filter(!is.na(SE) & SE!="SE" & Variance!="of") %>%
  mutate(Var = case_when(Source %in% c("V(G)/Vp_tr1") ~ "V(O)/Vp1",
                         Source %in% c("V(G)/Vp_tr2") ~ "V(O)/Vp2",
                         Source == "rG" ~ "rG",
                         Source == "Pval" ~ "Pval")) %>%
  select(-Source) %>% filter(!is.na(Var)) 

colnames(traits)[1] <- "Trait1"
colnames(traits)[2] <- "Trait2"

traits <- traits %>% mutate(Trait1a = ifelse(Trait1=="whr" & Trait2=="body_fat", "body_fat", Trait1),
                            Trait2a = ifelse(Trait1=="whr" & Trait2=="body_fat", "whr", Trait2)) %>%
  select(-Trait1, -Trait2, Variance , SE, Var, Trait1=Trait1a, Trait2= Trait2a ) 


traits$Trait1 <- recode_factor( traits$Trait1, 
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
traits$Trait2 <- recode_factor( traits$Trait2, 
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

traits_wide <- traits %>%
  nest(col_value = c(Variance, SE)) %>%
  spread(key = Var, value = col_value) %>%
  unnest(c(`V(O)/Vp1`, `V(O)/Vp2`, `rG`, `Pval`),  .sep = '_')

# options(knitr.kable.NA = '')
# traits_wide %>% 
#   mutate(rG = rG_Variance %>% signif(digits = 2)) %>% 
#   select(Trait1, Trait2, rG) %>%
#   spread(key = Trait2, value = rG) %>%
#   kable(caption="Epigenetic Correlations") 


# Heatmap 
# library(RColorBrewer)
# library(ggplot2)
library(forcats)

plot_ds <- traits_wide %>% filter(!Trait1 %in% c("Waist Circumference", "Weight") & !(Trait2 %in% c("Hip Circumference", "Weight", "Waist Circumference")) ) %>%
  mutate(Trait1 = fct_reorder(Trait1, desc(Trait1)),
         label = paste0(rG_Variance %>% signif(digits = 2), " (", rG_SE %>% signif(digits = 1), ")"))
# label = paste0(rG_Variance %>% signif(digits = 2))) %>% select(Trait1, Trait2, rG_Variance, label)
plot_ds$Trait1 <- plot_ds$Trait1 %>% as.character
plot_ds$Trait2 <- plot_ds$Trait2 %>% as.character
plot_ds$rG_Variance <- plot_ds$rG_Variance %>% as.character


plot_ds <- plot_ds %>% select(Trait1, Trait2, rG_Variance, label)
plot_ds[nrow(plot_ds) + 1,] <- c("BMI", "BMI", NA, NA) %>% as.list
plot_ds[nrow(plot_ds) + 1,] <- c("Waist to Hip Ratio", "Waist to Hip Ratio", NA, NA) %>% as.list
plot_ds[nrow(plot_ds) + 1,] <- c("Body Fat %", "Body Fat %", NA, NA) %>% as.list
plot_ds[nrow(plot_ds) + 1,] <- c("Glucose", "Glucose",  NA, NA) %>% as.list
plot_ds[nrow(plot_ds) + 1,] <- c("HDL Cholesterol", "HDL Cholesterol",  NA, NA) %>% as.list
plot_ds[nrow(plot_ds) + 1,] <- c("Total Cholesterol", "Total Cholesterol", NA, NA) %>% as.list
plot_ds[nrow(plot_ds) + 1,] <- c("Height", "Height",  NA, NA) %>% as.list
plot_ds[nrow(plot_ds) + 1,] <- c("Total Cholesterol", "HDL Cholesterol",  NA, NA) %>% as.list


plot_ds$Trait1 <- plot_ds$Trait1 %>% as.factor
plot_ds$Trait2 <- plot_ds$Trait2 %>% as.factor
plot_ds$Trait1 <- recode_factor( plot_ds$Trait1, 
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
plot_ds$Trait2 <- recode_factor( plot_ds$Trait2, 
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
plot_ds$rG_Variance <- plot_ds$rG_Variance %>% as.numeric

plot_ds2 <- plot_ds %>% 
  filter( !(Trait1=="Waist to Hip Ratio" & Trait2 %in% c("BMI")) &
            !(Trait1=="Total Cholesterol" & Trait2 %in% c("BMI",  "Waist to Hip Ratio", "Body Fat %", "Glucose", "HDL Cholesterol","Height")) &
            !(Trait1=="Height" & Trait2 %in% c("BMI", "Waist to Hip Ratio")))


r_G_plot <- plot_ds2
rg_SNP <- ggplot(plot_ds2 %>% filter(Trait1!="Height" & Trait2!="Height" & Trait1!="Total Cholesterol" & Trait2!="" & !(Trait1=="BMI" & Trait2=="BMI")),
                 aes(x = Trait2, y = Trait1, fill = rG_Variance)) +
  geom_raster() +
  geom_text(aes(label = label), size=3) +
  # scale_fill_distiller(palette = "Spectral", na.value = 'white',  limits=c(-0.7,1)) +
  scale_fill_gradient2(low = mycolours[1], mid = "beige", high = mycolours[3], na.value = 'white',  limits=c(-0.7,1), guide = "colourbar")  +
  theme_minimal() +
  theme(panel.grid = element_blank())+ 
  labs(fill = "Correlation", y="", x="") +
  scale_y_discrete(limits=rev) +
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(legend.title = element_text( size=12), legend.text=element_text(size=12)) 








###################
# Combine plots
##################

s <- plot_spacer()
library(patchwork)
rg_DNAm + rg_SNP + rg_pheno + plot_layout(ncol = 3, nrow = 1)

# Add option to share legend

library(ggpubr)
ggarrange(
  rg_DNAm +
    theme(axis.title.x=element_blank(),
          # axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          # axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(legend.position="none") +
    ggtitle("DNAm"),
  rg_SNP +
    theme(axis.title.x=element_blank(),
          # axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(legend.position="none") +
    ggtitle("Genetic"),
  rg_pheno +
    theme(axis.title.x=element_blank(),
          # axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(legend.position="right") +
    ggtitle("Phenotypic"), 
  ncol=3, widths= c(1, 0.85, 1.0))



## Thesis figures

library(ggpubr)

ggarrange(
  rg_DNAm +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          # axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(legend.position="none") +
    ggtitle("DNAm"),
  rg_SNP +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          # axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(legend.position="none") +
    ggtitle("Genetic"),
  rg_pheno +
    theme(axis.title.x=element_blank(),
          # axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          # axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(legend.position="none") +
    ggtitle("Phenotypic"), 
  ncol=1, heights= c(1, 1, 1.5))

merge_heatmap <- rbind(
  r_DNAm_plot %>% select(Trait1, Trait2,  rG_Variance, label) %>% mutate(var="DNAm"),
  r_G_plot %>% select(Trait1, Trait2,  rG_Variance, label) %>% mutate(var="Genetic"),
  r_pheno_plot %>% select(Trait1, Trait2,  rG_Variance, label) %>% mutate(var="Phenotypic")) %>%
  filter(Trait2!="BMI" & Trait1!="Height" & Trait2!="Height" & Trait1!="Total Cholesterol" & Trait1!=Trait2) %>%
  mutate(rG_Variance = rG_Variance %>% as.numeric %>% round(2))


ggplot(merge_heatmap ,
       aes(x = Trait2, y = Trait1, fill = rG_Variance)) +
  geom_raster() +
  geom_text(aes(label =  rG_Variance), size=3) +
  theme_minimal() +
  theme(panel.grid = element_blank())+ 
  labs(fill = "Correlation", y="", x="") +
  scale_y_discrete(limits=rev) +
  scale_fill_gradient2(low = mycolours[1], mid = "beige", high = mycolours[3], na.value = 'white',  limits=c(-0.7,1), guide = "colourbar")  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  theme(axis.text.y = element_text(size=8)) +
  theme(legend.title = element_text( size=8), legend.text=element_text(size=7)) +
  facet_grid(~var) +
  theme(strip.text.x = element_text(size = 11, colour = "black")) +
  theme(panel.spacing = unit(0.1, "lines"))







