
###################################
## DNAm Correlation between waves
###################################

# Bivariate REML

# bmi - columns $3, $4
# whr - columns $5, $6
# body_fat - columns $7, $8
# Glucose - columns $9, $10
# HDL_cholesterol - columns $11, $12
# Total_cholesterol - columns $13, $14

scratch=/Local_Scratch/Alesha
cov=covar_batch_eversmoke_packyears
filename=mvals-norm20k-18413-831733_sd0.02_unrelated_${cov}

mkdir ${scratch}/DNAm/phenotypes/greml.pheno/temp
cd ${scratch}/DNAm/phenotypes/greml.pheno/temp

for trait in bmi whr body_fat Glucose HDL_cholesterol Total_cholesterol height weight waist hips
do

for ds in 1v3 1v4 3v4
do

if [ "${trait}" == "bmi" ];              then  var1=3 var2=4; fi
if [ "${trait}" == "whr" ];              then  var1=5 var2=6; fi
if [ "${trait}" == "body_fat" ];         then  var1=7 var2=8; fi
if [ "${trait}" == "Glucose" ];          then  var1=9 var2=10; fi
if [ "${trait}" == "HDL_cholesterol" ];  then  var1=11 var2=12; fi
if [ "${trait}" == "Total_cholesterol" ]; then  var1=13 var2=14; fi
if [ "${trait}" == "height" ]; then  var1=15 var2=16; fi


awk -v "col1=${var1}" -v "col2=${var2}" '{print $1, $2, $col1, $col2}' ${scratch}/DNAm/phenotypes/greml.pheno/ds${ds}.txt > temp.${trait}.waves
gcta64 --reml-bivar 1 2 \
--grm ${scratch}/DNAm/ORM/GRM/${filename} \
--reml-bivar-lrt-rg 0 \
--pheno temp.${trait}.waves \
--reml-maxit 200 \
--out ${scratch}/DNAm/GREML/waves/${trait}_RankNorm${ds}_0

done
done

cd ${scratch}/DNAm/GREML/waves
# for trait in bmi whr body_fat Glucose HDL_cholesterol Total_cholesterol height 
for trait in Total_cholesterol height 
do
for ds in 1v3 1v4 3v4
do
awk -v "trait=${trait}" -v ds=${ds} '{print trait, ds, $1, $2, $3}' ${trait}_RankNorm${ds}.hsq >> ${trait}_waves.rsq
done
done


####################
# Output
##################

library(stringr)
library(readr)
traits <- read_table("/Users/uqahatt2/Documents/2_project/4_GS/1_Correlations/data/GREML_waves/trait_waves.rsq") %>% 
  mutate(Var = case_when(Source %in% c("V(G)/Vp_tr1") ~ "VO/Vp1",
                         Source %in% c("V(G)/Vp_tr2")   ~ "VO/Vp2",
                         Source %in% c("rG") ~ "rG",
                         Source %in% c("Pval") ~ "Pval")) %>% 
  select(-Source) %>%
  filter(!is.na(Var)) 
colnames(traits)[1] <- "Trait"
colnames(traits)[2] <- "Waves"

traits$Trait <- recode_factor( traits$Trait, 
                               bmi  = "BMI", 
                               body_fat = "Body Fat %", 
                               whr = "Waist to Hip Ratio",
                               Glucose = "Glucose",
                               HDL_cholesterol = "HDL Cholesterol",
                               Total_cholesterol = "Total Cholesterol",
                               height = "Height")

# Recode to numeric
traits$Variance <- traits$Variance %>% as.numeric %>% round(3)
traits$SE <- traits$SE %>% as.numeric %>% round(3)

# Wide

traits_wide_table <- traits %>%
  mutate(var2 = paste0(Waves, " ", Var)) %>%
  select(-Waves, -Var) %>%
  nest(col_value = c(Variance, SE)) %>%
  spread(key = c(var2), value = col_value) %>%
  unnest(c(`1v3 VO/Vp1`, `1v3 VO/Vp2`, `1v3 rG`, `1v4 VO/Vp1`, `1v4 VO/Vp2`, `1v4 rG`, 
           `3v4 VO/Vp1`, `3v4 VO/Vp2`, `3v4 rG`),  .sep = '_') 

traits_wide_table %>% as.data.frame
# 
#   
#   mutate(`Vwave1/Vp` = case_when(Waves=="1v3" ~ `VO/Vp1_Variance`,
#                                  Waves=="1v4" ~ `VO/Vp1_Variance`),
#          `Vwave3/Vp` = case_when(Waves=="1v3" ~ `VO/Vp2_Variance`,
#                                  Waves=="3v4" ~ `VO/Vp1_Variance`),
#          `Vwave4/Vp` = case_when(Waves=="1v4" ~ `VO/Vp2_Variance`,
#                                  Waves=="3v4" ~ `VO/Vp1_Variance`))
#          
#          

traits %>% filter(Waves=="1v3" & Var=="rG") %>% 
  mutate(`Wave1 vs Wave3` = paste0(Variance, " (", SE, ")")) %>% select(Trait, `Wave1 vs Wave3`) %>%
  left_join(traits %>% filter(Waves=="1v4" & Var=="rG") %>% 
              mutate(`Wave1 vs Wave4` = paste0(Variance, " (", SE, ")")) %>% select(Trait, `Wave1 vs Wave4`), by="Trait") %>%
  left_join(traits %>% filter(Waves=="3v4" & Var=="rG") %>% 
              mutate(`Wave3 vs Wave4` = paste0(Variance, " (", SE, ")")) %>% select(Trait, `Wave3 vs Wave4`), by="Trait") %>%
  kable(caption="Epigenetic correlation between waves")


# Figure
# figure_ds <- traits %>% filter(Var=="rG") 
figure_ds <- traits %>% filter(Var=="rG") %>% 
  left_join(traits %>% filter(Var=="Pval") %>% select(Trait, Pval=Variance, Waves), by=c("Trait", "Waves")) %>%
  mutate(pval_annot = Pval %>% signif(1))


figure_ds$Set <- recode_factor(figure_ds$Waves, 
                               `1v3` = "Set 1 v Set 2", 
                               `1v4` = "Set 1 v Set 3", 
                               `3v4` = "Set 2 v Set 3")


figure_ds %>% filter(Trait!="Height") %>%
  ggplot(aes(y=fct_reorder(Trait, desc(Trait)), x=Variance, xmin=Variance-SE, xmax = Variance+SE, colour=Trait, label=pval_annot)) +
  geom_vline(xintercept = 1, colour="darkgrey") +
  geom_errorbar(position = position_dodge(width = 0.2), width=.1) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_text(aes(label=paste0("p = ",pval_annot), x=1.04), hjust=0, vjust=-1.5, size=2.5,
            position = position_dodge(width = 1), colour="black") +
  labs(x="rG (DNAm) with BMI", y="Trait") +
  facet_wrap(~Set) +
  theme_linedraw() +
  scale_colour_manual(values = mycolours, name = " ") +
  labs(y=" ", x="DNAm correlation") +
  theme(legend.position="bottom") +
  theme(legend.position = "none")

