rm(list=ls())

setwd ("F://Box Sync/BHASKAR UT AUSTIN/Research/WUE_20-21/Fresh/GWAS_WUE_v4/WUE_V4_Rfiles/LT_jan22/")

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(readxl)
library(cowplot)
library(patchwork)

df3 <- read_excel("LT_jan.xlsx")

df176 <- df3 %>% 
  filter((ImageNo == "176") & mutant %in% c("Col","mef11-5"))

df176

my_comparisons <- list(c("Col", "mef11-5"))

my_comparisons

compare_means(LT~mutant, data = df176, method = "t.test")


p1 <-ggboxplot(df176, x = "mutant", y = "LT",
               color = "mutant", palette = "jco", add = "jitter")+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")+
  theme(legend.position = "NA")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_blank())
p1 

###pco5-2
df203 <- df3 %>% 
  filter((ImageNo == "203") & mutant %in% c("Col","pco5-1"))

df203

my_comparisons <- list(c("Col", "pco5-1"))

my_comparisons

compare_means(LT~mutant, data = df184, method = "t.test")


p2 <-ggboxplot(df203, x = "mutant", y = "LT",
               color = "mutant", palette = "jco", add = "jitter")+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")+
 
  theme(legend.position = "NA")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_blank())
p2 


##pco5-1
df184 <- df3 %>% 
  filter((ImageNo == "184") & mutant %in% c("Col","pco5-2"))

df184

my_comparisons <- list(c("Col", "pco5-2"))

my_comparisons

compare_means(LT~mutant, data = df184, method = "t.test")


p3 <-ggboxplot(df184, x = "mutant", y = "LT",
               color = "mutant", palette = "jco", add = "jitter")+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")+
  theme(legend.position = "NA")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_blank())
p3

#cpk23

df197 <- df3 %>% 
  filter((ImageNo == "197") & mutant %in% c("Col","cpk23"))

df197

my_comparisons <- list(c("Col", "cpk23"))

my_comparisons

compare_means(LT~mutant, data = df197, method = "t.test")


p4 <-ggboxplot(df197, x = "mutant", y = "LT",
               color = "mutant", palette = "jco", add = "jitter")+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")+
  theme(legend.position = "NA")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_blank())
p4 


#dur-1

df185 <- df3 %>% 
  filter((ImageNo == "185") & mutant %in% c("Col","dur-1"))

df185

my_comparisons <- list(c("Col", "dur-1"))

my_comparisons

compare_means(LT~mutant, data = df185, method = "t.test")


p5 <-ggboxplot(df185, x = "mutant", y = "LT",
               color = "mutant", palette = "jco", add = "jitter")+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")+
  theme(legend.position = "NA")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_blank())
p5 

##nsra-1

df193 <-  df3 %>% 
  filter((ImageNo == "193") & mutant %in% c("Col","nsra-1")) 

df193

my_comparisons <- list(c("Col", "nsra-1"))

my_comparisons

compare_means(LT~mutant, data = df193, method = "t.test")

p6 <-ggboxplot(df193, x = "mutant", y = "LT",
               color = "mutant", palette = "jco", add = "jitter")+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")+
  theme(legend.position = "NA")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_blank())

p6 

###nsra-2

df203 <-  df3 %>% 
  filter((ImageNo == "203") & mutant %in% c("Col","nsra-2"))

df203

my_comparisons <- list(c("Col", "nsra-2"))

my_comparisons

compare_means(LT~mutant, data = df203, method = "t.test")


p7 <-ggboxplot(df203, x = "mutant", y = "LT",
               color = "mutant", palette = "jco", add = "jitter")+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")+
  theme(legend.position = "NA")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_blank())
p7


LT_high <- plot_grid(p2, p3, p4, p5, p6,p7, p1,
                     ncol = 7)
LT_high

# p <- LT_high + 
#   plot_annotation(tag_levels = "A") &
#   theme(plot.tag = element_text(size = "15",  face = "bold"))


ggsave(LT_high, filename = "Fig3Snew_LT.tiff", width =9, height =4 , dpi = 300)


ggsave(LT_high, filename = "Fig3Snew_LT.svg", width =6, height =4 , dpi = 300)



ggsave(LT_high, filename = "Fig3Snew_LT.pdf", width = 11, height = 4)
