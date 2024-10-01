install.packages("ggstance")

library("ggplot2")
library(Rmisc)
library(ggpubr)
library(tidyr)
library(ggstance)
library(nlme)
library(emmeans)
library(ggbreak)

setwd("C:/Users/sfindeisen/Desktop/Heathland/VeluweData/data fungi Veluwe")

pal1=c("#52854C","#0E958F", "#D16103","#7B2D26")

#### Data prep ####
fungi <- read.table("effect_fungi.txt", h=T)
fungi$Treatment <- factor(fungi$Treatment , levels=c("Control", "Soilfeed", "Biolit","Dolokal"))
fungi$Site <- factor(fungi$Site, levels = c("Dry","Wet"))
str(fungi)

fungi1 <- fungi %>%
  pivot_longer(
    c("AMF", "EMF", "Endo", "ErMF", "PP", "S", "U"),
    names_to="level",
    values_to="taxon"
    ) %>%
      mutate(
        level = factor(level,
                              levels=c("U", "S", "PP", "Endo","ErMF","EMF","AMF" ))
      )

#### Treatment plots - DOLOKAL

subset_data <- fungi1[fungi1$Treatment == "Dolokal", ]
sum_Do = summarySE(subset_data,
                     measurevar=c("taxon"),
                     groupvars=c("level", "Site"),
                     na.rm=TRUE)#in this last you could add more variables

sum_Dodry <- sum_Do[sum_Do$Site=="Dry",]
sum_DoWet <- sum_Do[sum_Do$Site=="Wet",]


dol_pointrange <- ggplot(data = sum_Dowet, aes(x = level, y = taxon)) + 
  geom_pointrange(data = sum_Dodry, aes(ymin = taxon - se, ymax = taxon + se, shape = Site), position = position_dodge(0.5), size = 1, color="#D16103") + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", fill = NA)
  ) + 
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size=12)
  ) +
  scale_x_discrete(labels=c("Endo"= "Endophytes", "PP"="Plant pathogens", "S"="Saprotrophs"  ,"U"='unassigned'))+
  coord_flip() + 
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black", linewidth = 0.5)


### Control ###

subset_data <- fungi1[fungi1$Treatment == "Control", ]
sum_CO = summarySE(subset_data,
                     measurevar=c("taxon"),
                     groupvars=c("level", "Site"),
                     na.rm=TRUE)#in this last you could add more variables


sum_COdry <- sum_CO[sum_CO$Site=="Dry",]
sum_CoWet <- sum_CO[sum_CO$Site=="Wet",]



con_pointrange <- ggplot(data = sum_AMns, aes(x = level, y = taxon, shape = Site)) + 
  scale_shape_manual(values = c("Wet" = 16, "Dry" = 17))+
  geom_pointrange(data = sum_AMns, aes(ymin = taxon - se, ymax = taxon + se, shape = Site), position = position_dodge(0.5), size = 1, color="#52854C") + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", fill = NA)
  ) + 
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size=12)
  ) +
  scale_x_discrete(labels=c("Endo"= "Endophytes", "PP"="Plant pathogens", "S"="Saprotrophs"  ,"U"='unassigned'))+
  coord_flip() + 
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black", linewidth = 0.5)


### Biolit ###

subset_data <- fungi1[fungi1$Treatment == "Biolit", ]
sum_Bi = summarySE(subset_data,
                     measurevar=c("taxon"),
                     groupvars=c("level", "Site"),
                     na.rm=TRUE)#in this last you could add more variables
sum_Bidry <- sum_Bi[sum_Bi$Site=="Dry",]
sum_Biwet <- sum_Bi[sum_Bi$Site=="Wet",]

bio_pointrange <- ggplot(data = sum_AMns, aes(x = level, y = taxon, shape = Site)) + 
  scale_shape_manual(values = c("Wet" = 16, "Dry" = 17))+
  geom_pointrange(data = sum_AMns, aes(ymin = taxon - se, ymax = taxon + se, shape = Site), position = position_dodge(0.5), size = 1, color="#0E958F") + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", fill = NA)
  ) + 
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size=12)
  ) +
  scale_x_discrete(labels=c("Endo"= "Endophytes", "PP"="Plant pathogens", "S"="Saprotrophs"  ,"U"='unassigned'))+
  coord_flip() + 
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black", linewidth = 0.5)


### Soilfeed ###


subset_data <- fungi1[fungi1$Treatment == "Soilfeed", ]
sum_So = summarySE(subset_data,
                     measurevar=c("taxon"),
                     groupvars=c("level", "Site"),
                     na.rm=TRUE)#in this last you could add more variables
sum_Sodry <- sum_So[sum_So$Site=="Dry",]
sum_Sowet <- sum_So[sum_So$Site=="Wet",]


soil_pointrange <- ggplot(data = sum_AMns, aes(x = level, y = taxon, shape = Site)) + 
  scale_shape_manual(values = c("Wet" = 16, "Dry" = 17))+
  geom_pointrange(data = sum_AMns, aes(ymin = taxon - se, ymax = taxon + se, shape = Site), position = position_dodge(0.5), size = 1, color="#7B2D26") + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", fill = NA)
  ) + 
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size=12)
  ) +
  scale_x_discrete(labels=c("Endo"= "Endophytes", "PP"="Plant pathogens", "S"="Saprotrophs"  ,"U"='unassigned'))+
  coord_flip() + 
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black", linewidth = 0.5)


## FINAL PLOT ##

ggarrange(con_pointrange,bio_pointrange, dol_pointrange,soil_pointrange, nrow = 2 ,ncol=2,common.legend = T)  

ylab = expression("Ratio of relative abundance by the mean of control plots")

soil_pointrange <- ggplot(data = sum_CoWet, aes(x = level, y = taxon)) +
  #geom_pointrange(data = sum_CoWet, aes(ymin = taxon - se, ymax = taxon + se), shape =16, size = 1.5, color = "#52854C") + 
  geom_pointrange(data = sum_DoWet, aes(ymin = taxon - se, ymax = taxon + se), shape=16, position = position_nudge(), size = 1, color = "#7B2D26") +  # Add another layer with different data
  geom_pointrange(data = sum_Sowet, aes(ymin = taxon - se, ymax = taxon + se), shape = 16, position = position_nudge(x= 0.4), size = 1,color = "#0E958F" )+
  geom_pointrange(data = sum_Biwet, aes(ymin = taxon - se, ymax = taxon + se), shape = 16,position = position_nudge(x= 0.2), size = 1, color = "#D16103") +  # Add another layer with different data
  # Add another layer with different data
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", fill = NA)
  ) + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size=12),
    axis.title = element_blank(),
  ) +
  scale_x_discrete(labels = c("Endo" = "Endophytes", "PP" = "Plant pathogens", "S" = "Saprotrophs", "U" = 'unassigned')) +
  scale_y_continuous(trans = "log2", limits = c(0.1, NA))+
  coord_flip() + 
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black", linewidth = 0.5)
 # geom_vline(aes(xintercept=level))

fungiwet <- soil_pointrange + scale_y_break(c(3,6), scales = 0.5)

soil_pointrangedry <- ggplot(data = sum_CO, aes(x = level, y = taxon)) +
  #geom_pointrange(data = sum_COdry, aes(ymin = taxon - se, ymax = taxon + se), shape =17, size = 1.5, color = "#52854C") + 
  geom_pointrange(data = sum_Dodry, aes(ymin = taxon - se, ymax = taxon + se)), shape=16, position = position_nudge(), size = 1, color = "#7B2D26") +  # Add another layer with different data
  geom_pointrange(data = sum_Bidry, aes(ymin = taxon - se, ymax = taxon + se)), shape = 16,position = position_nudge(x= 0.2), size = 1, color = "#D16103") +  # Add another layer with different data
  geom_pointrange(data = sum_Sodry, aes(ymin = taxon - se, ymax = taxon + se), shape = 16, position = position_nudge(x= 0.4), size = 1,color = "#0E958F" ) +  # Add another layer with different data
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", fill = NA)
  ) + 
  theme(
    axis.title  = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size=12)
  ) +
  scale_x_discrete(labels = c("Endo" = "Endophytes", "PP" = "Plant pathogens", "S" = "Saprotrophs", "U" = 'unassigned')) +
  coord_flip() + 
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black", linewidth = 0.5)
# geom_vline(aes(xintercept=level))

final <- ggarrange(soil_pointrangedry, soil_pointrange)

output_dir <- "OUTPUT"
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
file_path <- file.path(output_dir, "effectsfungiwet.png")
ggsave(file = file_path, plot = soil_pointrange, width = 6, height = 6, dpi = 300)  # Adjust width, height, and dpi as needed