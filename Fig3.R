rm(list=ls())

setwd("F:/Box Sync/BHASKAR UT AUSTIN/Research/WUE_20-21/Fresh/GWAS_WUE_V6/Fig3/")

library(tidyverse)
library(viridis)
library(patchwork)
library(patchwork)


df3 <- read.csv("..//y1_y2_WUEV6.csv")
glimpse(df3)

df3$WUE <- as.numeric(df3$WUE)
glimpse(df3)

df4 <- df3 %>% select(Gene_ID,Genotype, Gene_Name, mutant, mgDW, Water_loss, WUE, 
                      delta,SD,delta_LT, Index, Pheno, Batch, Phe_status)%>% 
  filter(Pheno %in% c ("high_WUE", "Columbia"))

df4

df4$mutant <- factor(df4$mutant)

levels(df4$mutant)

#theme

WUE_theme <- theme_classic(base_size = 16) + 
  theme(axis.text = element_text(color = "black", size = 16 ),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1, face = "bold.italic"),
        axis.text.y = element_text(face = 'bold'),
        axis.title = element_text(color = "black", face = "bold"),
        axis.ticks = element_line(size = 1),
        axis.ticks.length.x = unit(-5, "pt"),
        axis.ticks.length.y = unit(5, "pt"),
        axis.ticks.x = element_line(color = 'black'),
        axis.line = element_line(size = 1),
        legend.position = "NA")


WUE_theme1 <- theme_classic(base_size = 16) + 
  theme(axis.text = element_text(color = "black",size = 16),
  #axis.text.x = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1, face = "bold.italic"),
  axis.text.y = element_text(face = 'bold'),
  axis.title = element_text(color = "black", face = "bold"),
  axis.ticks = element_line(size = 1),
  axis.ticks.length.x = unit(-5, "pt"),
  axis.ticks.length.y = unit(5, "pt"),
  axis.ticks.x = element_line(color = 'black'),
  axis.line = element_line(size = 1),
  legend.position = "NA")
                                                  
####WUE

p1<- df4 %>% 
  ggplot(aes(x=reorder(mutant, WUE, na.rm = TRUE), y = WUE))+
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, linetype = 1, size = 0.5, width = 0.5, coef = 1.5, alpha = 0.5)+
  geom_jitter(aes (fill = Gene_Name), alpha = 1, width = 0.2, size = 3, shape = 21)+
  scale_x_discrete(limits = c("Col", "pco5-1", "pco5-2", "cpk23", "dur-1", "nsra-1", "nsra-2", "mef11-5"),
                   labels = c("Col", "pco5-1", "pco5-2", "cpk23", "dur-1", "nsra-1", "nsra-2", "mef11-5"))+
  scale_y_continuous(breaks = seq(1, 2.2, 0.4),
                     limits = c(1, 2.2))+
  scale_fill_brewer(palette = 'Paired')+
  labs(y = "WUE(g/L)",
       x = "")+
  WUE_theme

p1

####biomass

p2 <- df4 %>% 
  mutate(mgDW = mgDW/1000) %>% 
  ggplot(aes(x=reorder(mutant, mgDW, na.rm = TRUE), y = mgDW))+
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, linetype = 1, size = 0.5, width = 0.5, coef = 1.5, alpha = 0.5)+
  geom_jitter(aes (fill = Gene_Name), alpha = 1, width = 0.2, size = 3, shape = 21)+
  scale_x_discrete(limits = c("Col", "pco5-1", "pco5-2", "cpk23", "dur-1", "nsra-1", "nsra-2", "mef11-5"),
                   labels = c("Col", "pco5-1", "pco5-2", "cpk23", "dur-1", "nsra-1", "nsra-2", "mef11-5"))+
  scale_y_continuous(limits = c(0, 0.25),
                     breaks = seq(0, 0.25,0.05))+
  #scale_fill_viridis_d(direction =1)+
  scale_fill_brewer(palette = "Paired")+
  labs(y = "Biomass(g)",
       x = "")+
  WUE_theme

p2

###Consumed Water

p3<- df4 %>% 
  ggplot(aes(x=reorder(mutant, Water_loss, na.rm = TRUE), y = Water_loss))+
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, linetype = 1, size = 0.5, width = 0.5, coef = 1.5, alpha = 0.5)+
  geom_jitter(aes (fill = Gene_Name), alpha = 1, width = 0.2, size = 3, shape = 21)+
  scale_x_discrete(limits = c("Col", "pco5-1", "pco5-2", "cpk23", "dur-1", "nsra-1", "nsra-2", "mef11-5"),
                   labels = c("Col", "pco5-1", "pco5-2", "cpk23", "dur-1", "nsra-1", "nsra-2", "mef11-5"))+
  scale_y_continuous(breaks = seq(0, 125, 25),
                     limits = c(0, 125))+
  scale_fill_brewer(palette = 'Paired')+
  labs(y = "Consumed water (mL)",
       x = "")+
  WUE_theme

p3


###delta 13C

p4<- df4 %>% 
  ggplot(aes(x=reorder(mutant, delta, na.rm = TRUE), y = delta))+
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, linetype = 1, size = 0.5, width = 0.5, coef = 1.5, alpha = 0.5)+
  geom_jitter(aes (fill = Gene_Name), alpha = 1, width = 0.2, size = 3, shape = 21)+
  scale_x_discrete(limits = c("Col", "pco5-1", "pco5-2", "cpk23", "dur-1", "nsra-1", "nsra-2", "mef11-5"),
                   labels = c("Col", "pco5-1", "pco5-2", "cpk23", "dur-1", "nsra-1", "nsra-2", "mef11-5"))+
  scale_y_continuous(breaks = seq(-33, -30, 1),
                     limits = c(-33, -30))+
  scale_fill_brewer(palette = 'Paired')+
  labs(y = expression(paste(delta^{13}, "C")),
       x = "")+
  WUE_theme

p4

##Stomatal Density

p5<- df4 %>%
  ggplot(aes(x=reorder(mutant, -SD, na.rm = TRUE), y = SD, fill = Gene_Name))+
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, linetype = 1, size = 0.5, width = 0.5, coef = 1.5, alpha = 0.5)+
  geom_jitter(aes (fill = Gene_Name), alpha = 1, width = 0.2, size = 3, shape = 21)+
  scale_x_discrete(limits = c("Col", "pco5-1", "pco5-2", "cpk23", "dur-1", "nsra-1", "nsra-2", "mef11-5"),
                   labels = c("Col", "pco5-1", "pco5-2", "cpk23", "dur-1", "nsra-1", "nsra-2", "mef11-5"))+
  scale_y_continuous(breaks = seq(100, 400, 50),
                     limits = c(100, 400))+
  #scale_fill_viridis_d(direction =1)+
  scale_fill_brewer(palette = "Paired")+
  labs(y = "Stomatal density",
       x = "")+
  WUE_theme1
p5


### leaf surface temperatures
p6<- df4 %>% 
  ggplot(aes(x= reorder(mutant, -delta_LT), y = delta_LT,  fill = Gene_Name ))+
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, linetype = 1, size = 0.5, width = 0.5, coef = 1.5, alpha = 0.5)+
  geom_jitter(aes (fill = Gene_Name), alpha = 1, width = 0.2, size = 3, shape = 21)+
  scale_x_discrete(limits = c("Col", "pco5-1", "pco5-2", "cpk23", "dur-1", "nsra-1", "nsra-2", "mef11-5"),
                   labels = c("Col", "pco5-1", "pco5-2", "cpk23", "dur-1", "nsra-1", "nsra-2", "mef11-5"))+
  geom_hline(yintercept = 0, color = "black", linetype = 1)+
  scale_y_continuous(breaks = seq(-1, 1, 0.5),
                     limits = c(-1, 1))+
  #scale_fill_viridis_d(direction =1)+
  scale_fill_brewer(palette = "Paired")+
  labs(y = expression(paste(Delta, "Leaf temp (°C)")),
       x = "")+
  WUE_theme1

p6


Fig2 <- p1/p2/p3/p4/p6

Fig2

ggsave(Fig2, filename = "Fig3_V6_R.png", dpi = 300, height = 13.5, width = 6)

ggsave(Fig2, filename = "Fig3_V6_R.svg", dpi = 300, height = 13.5, width = 6)


