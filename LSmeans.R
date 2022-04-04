rm(list=ls())

setwd ("F://Box Sync/BHASKAR UT AUSTIN/Research/WUE_20-21/Fresh/GWAS_WUE_v4/WUE_V4_Rfiles/")

library(tidyverse)
library(readxl)
library(ggpubr)
library(ggthemes)
library(lme4)
library(lsmeans)
library(emmeans)

##Large container experiment.

df3 <- read.csv("y1_y2_WUEV6.csv")#data file for the larger container experiment 
glimpse(df3)

df3$WUE <- as.numeric(df3$WUE)
glimpse(df3)

df4 <- df3 %>% select(Gene_Name,Gene_ID,Batch,mutant, mgDW,Water_loss, WUE, delta, SD, delta_LT)

### Mixed model, LS means and significance test for WUE. 

# At T-DNA LEVEL

# head(df4)
# df4$Gene_ID <- factor(df4$Gene_ID)
# levels(df4$Gene_ID)
 
mixed.lmer <- lmer(WUE~Gene_ID+(Gene_ID:mutant)+(1|Batch),  data = df4)#
summary(mixed.lmer)

emm_mutant_comp <- emmeans(mixed.lmer, ~Gene_ID)
emm_mutant_comp

contrast(emm_mutant_comp, "trt.vs.ctrl", ref = 27, adjust = "none")# ref = 27 corresponds to Columbia


#At GENE LEVEL

# df4$mutant <- factor(df4$mutant)
# 
# levels(df4$mutant)

mixed.lmer <- lmer(WUE~Gene_ID+ (Gene_ID:Genotype) +(1|RACK),  data = df_metaV1)#
summary(mixed.lmer)
emm_mutant_comp <- emmeans(mixed.lmer, ~Gene_ID)
emm_mutant_comp

contrast(emm_mutant_comp, "trt.vs.ctrl", ref = 71, adjust = "bonferroni")#ref = 71 corresponds to Columbia



##LS means and significance test for WUE measurements from small container experiment.

df_metaV1 <- read.csv("MetaV1 (1).csv") %>% #data file for smalle container experiment 
   select(Genotype, Gene_ID, WUE,mgDW, RACK) %>% 
   drop_na(WUE)
 
#   df_metaV1$Gene_ID <- factor(df_metaV1$Gene_ID)
#  
#  levels(df_metaV1$Gene_ID)
 
 #At GENE LEVEL

 mixed.lmer <- lmer(WUE~Gene_ID+ (Gene_ID:Genotype) +(1|RACK),  data = df_metaV1)#
 summary(mixed.lmer)
 emm_mutant_comp <- emmeans(mixed.lmer, ~Gene_ID)
 emm_mutant_comp
 
 contrast(emm_mutant_comp, "trt.vs.ctrl", ref = 71, adjust = "bonferroni")
 
 
 
 # At T_DNA LEVEL
 
 mixed.lmer <- lmer(mgDW~Genotype+(1|RACK), data = df_metaV1)#
 summary(mixed.lmer)
 emm_mutant_comp <- emmeans(mixed.lmer, ~Genotype)
 emm_mutant_comp
 
 contrast(emm_mutant_comp, "trt.vs.ctrl", ref = 1, adjust = "none")
 
 
 

