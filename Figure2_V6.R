rm(list=ls())

setwd ("F://Box Sync/BHASKAR UT AUSTIN/Research/WUE_20-21/Fresh/GWAS_WUE_V6/Fig2/")


library(tidyverse)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(ggthemes)
library(ggplot2)
library(ggrepel)
library(patchwork)

df <- read.csv("..//y1_y2_WUEV6.csv")

df$Water_loss<-as.numeric(as.character(df$Water_loss))
df$WUE<-as.numeric(as.character(df$WUE))

Columbia<-as_tibble(df) %>% select(Index,Gene_ID,Water_loss,mgDW,WUE) %>% 
  filter(Gene_ID=="Columbia") %>% group_by(Gene_ID) %>% 
  summarise(across(where(is.numeric), 
                   list(Mean=~ mean(.x, na.rm = TRUE)))
  ) %>% 
  ungroup()

## log2 fold change at T-DNA level
Meandf<- as_tibble(df) %>% select(Index,mutant,Water_loss,mgDW,WUE,Phe_status) %>% 
  filter(!mutant=="Col") %>% group_by(mutant,Phe_status) %>% 
  summarise(across(where(is.numeric), 
                   list(Mean=~ mean(.x, na.rm = TRUE)))
  ) %>% 
  ungroup() %>% 
  mutate(Fold_Water_loss_Mean= 
           log2(Water_loss_Mean/Columbia$Water_loss_Mean[1]),
         Fold_mgDW_Mean= 
           log2(mgDW_Mean/Columbia$mgDW_Mean[1]),
         Fold_WUE_Mean= 
           log2(WUE_Mean/Columbia$WUE_Mean[1]),
         Delta_Water_loss_Mean=
           Water_loss_Mean-Columbia$Water_loss_Mean[1],
         Delta_mgDW_Mean= 
           mgDW_Mean-Columbia$mgDW_Mean[1],
         Delta_WUE_Mean= 
           WUE_Mean-Columbia$WUE_Mean[1]
         
  ) 

##log2 fold change at gene level

Meandf2<- as_tibble(df) %>% select(Index,Gene_Name,Water_loss,mgDW,WUE,Phe_status) %>% 
  filter(!Gene_Name=="Columbia") %>% group_by(Gene_Name,Phe_status) %>% 
  summarise(across(where(is.numeric), 
                   list(Mean=~ mean(.x, na.rm = TRUE)))
  ) %>% 
  ungroup() %>% 
  mutate(Fold_Water_loss_Mean= 
           log2(Water_loss_Mean/Columbia$Water_loss_Mean[1]),
         Fold_mgDW_Mean= 
           log2(mgDW_Mean/Columbia$mgDW_Mean[1]),
         Fold_WUE_Mean= 
           log2(WUE_Mean/Columbia$WUE_Mean[1]),
         Delta_Water_loss_Mean=
           Water_loss_Mean-Columbia$Water_loss_Mean[1],
         Delta_mgDW_Mean= 
           mgDW_Mean-Columbia$mgDW_Mean[1],
         Delta_WUE_Mean= 
           WUE_Mean-Columbia$WUE_Mean[1]
         
  ) 

## Plot log2fold changes at T-DNA level

p1 <- Meandf %>%   
  ggplot(aes(x=Fold_Water_loss_Mean,y=Fold_mgDW_Mean)) +
  geom_point(aes(fill=Phe_status),size=4,alpha=0.8,shape=21,color="black")+
  #geom_text(aes(label=mutant))+
  geom_vline(xintercept = 0,linetype=2,size=0.5, color = "red")+
  geom_hline(yintercept = 0,linetype=2,size=0.5, color = "red")+
  scale_x_continuous(limits = c(-3.2,2.2))+
  scale_y_continuous(limits = c(-3.2,2.2))+
  scale_fill_manual(values = c ( "blue",  "purple", "skyblue", "#E69F00", "darkgreen", "#999999"),
                    labels = c ("HDBDW", "HIBIW",  "HNBDW", "LDBDW", "LDBNW", "LNBNW"))+
  labs(x="LFC mutant-Col (Consumed water)",
       # y="\u0394  x-axis title")+
       y="LFC mutant-Col (Biomass)")+
  theme_classic(base_size = 16) + 
  theme(axis.title = element_text(size = 16, color = "black", face = "bold"),
        axis.text = element_text(size = 16, color = "black", face = "bold"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))+
  guides(fill= FALSE, shape = FALSE)

p1


####Correlation analysis

#1. Get mean values for all traits from large container experiment

glimpse(df)
df$Water_loss<-as.numeric(as.character(df$Water_loss))
df$WUE<-as.numeric(as.character(df$WUE))

head(df)

df1 <- df %>% select(mutant, mgDW, Water_loss, WUE,delta, SD, SI) %>% 
  group_by(mutant) %>% 
  summarise_all(funs(mean = mean(., na.rm = T),sd = sd (., na.rm = T),se=sd(., na.rm = T)/sqrt(n()) )) 


### Add new variale WUE_status and name them based on high, low and no change in WUE. 

df2 <- df1 %>% mutate(WUE_status = 
                        case_when(WUE_mean > 1.284426462~"high_WUE",
                                  WUE_mean < 1.284426462~ "low_WUE",
                                  WUE_mean == 1.284426462~ "nochange"))

write.csv(df2, "means_yogurt_V5.csv")

df3 <- read.csv("means_yogurt_V5.csv")

###regression analysis and plot for Biomass v/s water loss

model1<-lm(mgDW_mean ~ Water_loss_mean, data = df3)
summary(model1)

names(df3)

plot1<-df3 %>% 
  mutate(mgDW_mean = mgDW_mean/1000) %>% 
  ggplot(aes(Water_loss_mean, mgDW_mean, fill = Phe_status))+
  geom_smooth(aes(x= Water_loss_mean, y = mgDW_mean), method = lm, se = FALSE,
              color = "black", linetype = 1, fill = "black", size = 1,
              formula = y~x)+
  geom_point (shape = 21, size = 4, color = "black", alpha =0.8, show.legend = FALSE)+
  scale_fill_manual(values = c ("black", "blue",  "purple",  "skyblue", "#E69F00", "darkgreen", "#999999"),
                    labels = c ("Columbia","HDBDW", "HIBIW",  "HNBDW", "LDBDW", "LDBNW", "LNBNW"))+
     labs(x = "Consumed water (mL)",
       y= "Biomass (g)")+
  theme_classic(base_size = 16)+
  annotate("text", size = 6, x = 40, y = 0.2, label = "italic(R) ^ 2 == 0.89",parse = TRUE)+
  theme(axis.title = element_text(size = 16,color = "black", face = "bold"),
        axis.text = element_text(size = 16, color = "black", face = "bold"),
        legend.text = element_text(size = 30, color = "black"))
       
  
plot1

### regression analysis and plot for delta 13C v/s WUE

df3 <- read.csv("means_yogurt_V5.csv")

#df3 <- df3 %>% filter(mutant != "cyp707a3-3")
#glimpse(df3)

model2 <- lm(WUE_mean ~ delta_mean, data = df3)
summary(model2)

plot2<-ggplot(data = df3,aes(delta_mean, WUE_mean, fill = Phe_status))+
  geom_smooth(data = df3,aes(x= delta_mean, WUE_mean), method = lm, se = FALSE,
              color = "black", linetype = 1, fill = "lightgray", size = 1,
              formula = y~x)+
  geom_point (shape = 21,size = 4, alpha = 0.8, color = "black", show.legend = FALSE)+
  scale_fill_manual(values = c ("black", "blue",  "purple",  "skyblue", "#E69F00", "darkgreen", "#999999"),
                    labels = c ("Columbia","HDBDW", "HIBIW",  "HNBDW", "LDBDW", "LDBNW", "LNBNW"))+
  # scale_x_continuous(limits = c(-34.2, -30),
  #                    breaks = seq(-34.2, -30, 1))+
  #scale_shape_manual(values = c(22, 24, 21))+
  labs(x = expression(paste(delta^{13}, "C")),
       y = "WUE(g/L)",
       color = "Phe_status")+
  theme_classic(base_size = 16)+
  annotate("text", size = 6, x = -33.5, y = 2, label = "italic(R) ^ 2 == 0.546",parse = TRUE)+
  theme(axis.title = element_text(size = 16,color = "black", face = "bold"),
        axis.text = element_text(size = 16, color = "black", face = "bold"),
        legend.text = element_text(size = 16))
  
plot2

### regression analysis and plot for WUE v/s Stomatal density
df3 <- read.csv("means_yogurt_V5.csv")

df3 <- df3 %>%  
   drop_na(SD_mean)

model3<-lm(WUE_mean ~ SD_mean, data= df3)
summary(model3)


plot3<-ggplot(data = df3,aes(SD_mean, WUE_mean, fill = Phe_status))+
  geom_smooth(data = df3,aes(x= SD_mean, WUE_mean), method = lm, se = FALSE,
              color = "black", linetype = 1, fill = "lightgray", size = 1,
              formula = y~x)+
  geom_point (shape = 21,size = 4, alpha = 0.8, color = "black", show.legend = FALSE)+
  scale_fill_manual(values = c ("black", "blue",  "purple",  "skyblue", "#E69F00", "darkgreen", "#999999"),
                    labels = c ("Columbia","HDBDW", "HIBIW",  "HNBDW", "LDBDW", "LDBNW", "LNBNW"))+
  scale_y_continuous(limits = c(0.5, 2),
                     breaks = seq(0.5, 2, 0.5))+
  labs(x = "Stomatal density",
       y= "WUE (g/L)")+
  theme_classic(base_size = 16)+
  annotate("text", size = 6, x = 190, y = 2, label = "italic(R) ^ 2 == 0.142",parse = TRUE)+
  theme(axis.title = element_text(size = 16, color = "black", face = "bold"),
        axis.text = element_text(size = 16, color = "black", face = "bold"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))+
  guides(fill= FALSE, shape = FALSE)

plot3

##regression analysis and plot for WUE v/s mgDW_mean (Biomass)
df3 <- read.csv("means_yogurt_V5.csv")

WUEvsBiomass<-lm(WUE_mean ~ mgDW_mean,data = df3)
summary(WUEvsBiomass)

plot4<-df3 %>% 
  mutate(mgDW_mean = mgDW_mean/1000) %>% 
  ggplot(aes(mgDW_mean, WUE_mean, fill = Phe_status))+
  geom_smooth(aes(x= mgDW_mean, WUE_mean), method = lm, se = FALSE,
              color = "black", linetype = 1, fill = "lightgray", size = 1,
              formula = y~x)+
  geom_point (shape = 21,size = 4, alpha = 0.8, color = "black", show.legend = FALSE)+
  scale_fill_manual(values = c ("black", "blue",  "purple",  "skyblue", "#E69F00", "darkgreen", "#999999"),
                    labels = c ("Columbia","HDBDW", "HIBIW",  "HNBDW", "LDBDW", "LDBNW", "LNBNW"))+
  scale_x_continuous(limits = c(0, 0.22),
                     breaks = seq(0, 0.22, 0.05))+
  scale_y_continuous(limits = c(0.5, 2),
                     breaks = seq(0.5, 2, 0.5))+
  #scale_shape_manual(values = c(22, 24, 21))+
  labs(x = "Biomass(g)",
       y= "WUE (g/L)")+
  theme_classic(base_size = 16)+
  annotate("text", size = 6, x = .02, y = 2, label = "italic(R) ^ 2 == 0.628",parse = TRUE)+
  theme(axis.title = element_text(size = 16, color = "black", face = "bold"),
        axis.text = element_text(size = 16, color = "black", face = "bold"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))+
  guides(fill= FALSE, shape = FALSE)

plot4

##regression analysis and plot for WUE v/s water_loss)
df3 <- read.csv("means_yogurt_V5.csv")

WUEvswater_loss<-lm(WUE_mean ~ Water_loss_mean, data = df3)
summary(WUEvswater_loss)

plot5<-df3 %>% 
  ggplot(aes(Water_loss_mean, WUE_mean, fill = Phe_status))+
  geom_smooth(aes(x= Water_loss_mean, WUE_mean), method = lm, se = FALSE,
              color = "black", linetype = 1, fill = "lightgray", size = 1,
              formula = y~x)+
  geom_point (shape = 21,size = 4, alpha = 0.8, color = "black", show.legend = FALSE)+
  scale_fill_manual(values = c ("black", "blue",  "purple",  "skyblue", "#E69F00", "darkgreen", "#999999"),
                    labels = c ("Columbia","HDBDW", "HIBIW",  "HNBDW", "LDBDW", "LDBNW", "LNBNW"))+
  scale_x_continuous(limits = c(30, 110),
                     breaks = seq(30, 110, 20))+
  scale_y_continuous(limits = c(0.5, 2),
                     breaks = seq(0.5, 2, 0.5))+
  labs(x = "Consumed water (mL)",
       y= "WUE (g/L)")+
  theme_classic(base_size = 16)+
  annotate("text", size = 6, x = 40, y = 2, label = "italic(R) ^ 2 == 0.364",parse = TRUE)+
  theme(axis.title = element_text(size = 16, color = "black", face = "bold"),
        axis.text = element_text(size = 16, color = "black", face = "bold"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))+
  guides(fill= FALSE, shape = FALSE)

plot5

Fig2NEW <- (p1+plot1+plot2)/(plot3+plot4+plot5)

Fig2NEW+plot_annotation(tag_levels = "A")


ggsave("Fig2_V6.tiff", height = 9, width = 14, dpi = 300)

ggsave("Fig2_V6.svg", height = 9, width = 14, dpi = 300 )


##correlation analysis for supporting information

#regression analysis and plot for WUE v/s Stomatal Index (SI)

df3 <- read.csv("means_yogurt_V5.csv")

df3 <- df3 %>%  
  drop_na(SI_mean)

WUEvsSI <- lm(WUE_mean~SI_mean, data = df3)
summary(WUEvsSI)

plotSI<-ggplot(data = df3,aes(SI_mean, WUE_mean, fill = Phe_status))+
  geom_smooth(data = df3,aes(x= SI_mean, WUE_mean), method = lm, se = FALSE,
              color = "black", linetype = 1, fill = "lightgray", size = 1,
              formula = y~x)+
  geom_point (shape = 21,size = 4, alpha = 0.8, color = "black", show.legend = FALSE)+
  scale_fill_manual(values = c ("black", "blue",  "purple",  "skyblue", "#E69F00", "darkgreen", "#999999"),
                    labels = c ("Columbia","HDBDW", "HIBIW",  "HNBDW", "LDBDW", "LDBNW", "LNBNW"))+
  scale_x_continuous(limits = c(38, 53),
                     breaks = seq(35, 53, 3))+
  #scale_shape_manual(values = c(22, 24, 21))+
  labs(x = "Stomatal Index",
       y = "WUE(g/L)",
       color = "Phe_status")+
  theme_classic(base_size = 16)+
  annotate("text", size = 6, x = 40, y = 2, label = "italic(R) ^ 2 == 0.083",parse = TRUE)+
  theme(axis.title = element_text(size = 16,color = "black", face = "bold"),
        axis.text = element_text(size = 16, color = "black", face = "bold"),
        legend.text = element_text(size = 16))

plotSI


ggsave("WUEvsSI_V5.tiff", units="in", width=5, height=5, dpi=300, compression = 'lzw')


###Global outlook of correlation analysis for all the traits
df3 <- read.csv("means_yogurt_V5.csv")

df4 <- df3 %>% 
  select(mutant, Water_loss_mean, mgDW_mean, WUE_mean, delta_mean, SD_mean, SI_mean) %>% 
  filter(mutant != "dreb2a-4") %>% 
  filter (mutant != "nsra-2") %>% 
  drop_na(SD_mean)

df5 <- df4 %>% 
  select(Water_loss_mean, mgDW_mean, WUE_mean, delta_mean, SD_mean, SI_mean)

df5

mod <- lm(mgDW_mean~Water_loss_mean, dfa = df5)
mod

cor(df5, method = "pearson")

cor(Water_loss_mean,mgDW_mean)

cor.test(df, method = "pearson")

