library(pheatmap)
library("vegan")
library("rlang")
library("ggplot2")
library("ggpubr")
library(dplyr)
library(ggbreak)
library(tidyr)
library(Rmisc)
library(stringr)
library(devtools)
library(nlme)
library(emmeans)
library(pairwiseAdonis)
library(ggforce)


install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

setwd("C:/Users/sfindeisen/Desktop/Heathland/VeluweData/data vegetation Veluwe")

pal1=c("#52854C","#0E958F", "#D16103","#7B2D26")

### Read in data ####

veg <- read.table("vegetationdata_complete.txt", h=T)
veg$Treatment <- factor(veg$Treatment, levels = c("control", "soilfeed","biolit", "dolokal"), labels = c("Control", "Soilfeed", "Biolit", "Dolokal"))
veg$Site <- factor(veg$Site, levels = c("dry", "wet"), labels = c("Dry", "Wet"))
row.names(veg) = veg$Plot_Code
str(veg)


###################################
#####NMDS analysis (only 2020)#####
###################################
veg2020 = subset(veg, Year == "2020")
veg = veg2020[,6:90]
veg = veg[, colSums(veg) > 0]
veg <- veg %>%
  mutate_all(as.numeric)

veg2020 = cbind(veg2020[,1:5], veg)
rownames(veg2020)=veg2020$Plot

veg2020=veg2020[order(match(rownames(veg2020),rownames(soil))),]

### dry + wet apart
vegdry = subset(veg2020, Site=="Dry")
vegwet = subset(veg2020, Site == "Wet")

# Species names list

species = as.data.frame(colnames(veg2020))
species$`colnames(veg2020)` <- gsub("_", " ", species$`colnames(veg2020)`)

write.table(species, "Species_list.txt")

veg2020$Taraxacum_officinale[veg2020$Treatment=="Control"]

## Average number of species per treatment

species_columns <- colnames(vegwet[,6:57])

# Calculate the average number of species per treatment
data <- vegwet %>%
  mutate(Species_Count = rowSums(select(., all_of(species_columns)) > 0))

# Calculate the average number of species counts per treatment
average_species_per_treatment <- data %>%
  group_by(Treatment) %>%
  summarise(Average_Species_Count = mean(Species_Count))

### dry + wet combined
ordr<-metaMDS(sqrt(veg), distance = "bray", k = 2, binary=F)  # stress = 8.712131e-05
summary(ordr)
plot(ordr$points)
ordr$stress

#Plot using ggplot2
#Extract the axes scores
datascores=as.data.frame(scores(ordr)$sites) 


#Add/calculate spider diagram
scores <- cbind(as.data.frame(datascores), Treatment = veg2020$Site)
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Treatment, data = scores, FUN = mean)
seg <- merge(scores, setNames(centroids, c('Treatment','oNMDS1','oNMDS2')),
             by = 'Treatment', sort = FALSE)

#plot
plotNMDS1bac = ggplot(scores, aes(x = NMDS1, y = NMDS2, colour = Treatment)) +
  geom_segment(data = seg,
               mapping = aes(xend = oNMDS1, yend = oNMDS2), linewidth=.7) + # add spiders
  geom_point(data = centroids, size = 3) +                    # add centroids
  geom_point(size=1) +                                              
  #coord_fixed()+                                              
  theme_bw()+ 
  theme(legend.position="top",legend.background = element_rect(color = "black", linetype = "solid"),
        legend.text=element_text(size=12),legend.title= element_blank(),legend.direction='horizontal',
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))+
  scale_color_manual(values = c("#4E84C4", "#293352"), labels= c("Dry", "Wet"))
plotNMDS1bac




#Add/calculate spider diagram
scores <- cbind(as.data.frame(datascores), Treatment = veg2020$Treatment)
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Treatment, data = scores, FUN = mean)
seg <- merge(scores, setNames(centroids, c('Treatment','oNMDS1','oNMDS2')),
             by = 'Treatment', sort = FALSE)

#plot
plotNMDS2bac = ggplot(scores, aes(x = NMDS1, y = NMDS2, colour = Treatment)) +
  geom_segment(data = seg,
               mapping = aes(xend = oNMDS1, yend = oNMDS2), linewidth = .7) + # add spiders
  geom_point(data = centroids, size = 3) +                    # add centroids
  geom_point(size=1) +                                              
  #coord_fixed()+                                              
  theme_bw()+ 
  theme(legend.position="top",legend.background = element_rect(color = "black", linetype = "solid"),
        legend.text=element_text(size=12),legend.title= element_blank(),legend.direction='horizontal',
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))+
  scale_color_manual(values =pal1, labels=c('Control', 'Biolit', 'Dolokal', 'Soilfeed'))
plotNMDS2bac

### Check RL / heathland species abundances ###

mean(veg2020$Luzula_campestris[veg2020$Treatment=="Dolokal"], na.rm = T)




#vegdry = vegdry[-5,]# delete E2 for procrustes analysis 
veg = vegdry[,6:57]
rownames(vegdry)
str(veg)

ordrvegdry<-metaMDS(sqrt(veg), distance = "bray", k = 2, binary=F)  # stress = 0.1382159
summary(ordr)
plot(ordr$points)
ordr$stress

#Plot using ggplot2
#Extract the axes scores
datascores=as.data.frame(scores(ordr)$sites) 

#Add/calculate spider diagram
scores <- cbind(as.data.frame(datascores), Treatment = vegdry$Treatment)
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Treatment, data = scores, FUN = mean)
seg <- merge(scores, setNames(centroids, c('Treatment','oNMDS1','oNMDS2')),
             by = 'Treatment', sort = FALSE)

#plot
plotNMDSvegdry = ggplot(scores, aes(x = NMDS1, y = NMDS2, colour = Treatment)) +
  geom_segment(data = seg,
               mapping = aes(xend = oNMDS1, yend = oNMDS2), linetype = "dashed") + # add spiders
  geom_point(data = centroids, size = 3) +                    # add centroids
  geom_point(size=1) +                                              
  #coord_fixed()+                                              
  theme_bw()+ 
  theme(legend.background = element_rect(color = "black", linetype = "solid"),
        legend.text=element_text(size=12),legend.title= element_blank(),legend.position="none",
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))+
  scale_color_manual(values =pal1)+
  scale_y_continuous(n.breaks = 3, limits = c(-0.5,0.5))
plotNMDSvegdry

### Tests dry

testdry <- adonis2(sqrt(veg) ~ Treatment, data=vegdry, permutations=9999, method="bray")
testdry

pairdry <- pairwise.adonis(sqrt(veg), vegdry$Treatment, perm = 9999)


### wet
veg = vegwet[,6:57]

ordrvegwet<-metaMDS(sqrt(veg), distance = "bray", k = 2, binary=F)  # stress =  0.161113
summary(ordr)
plot(ordr$points)
ordr$stress

#Plot using ggplot2
#Extract the axes scores
datascores=as.data.frame(scores(ordr)$sites) 

#Add/calculate spider diagram
scores <- cbind(as.data.frame(datascores), Treatment = vegwet$Treatment)
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Treatment, data = scores, FUN = mean)
seg <- merge(scores, setNames(centroids, c('Treatment','oNMDS1','oNMDS2')),
             by = 'Treatment', sort = FALSE)

#plot
plotNMDSvegwet = ggplot(scores, aes(x = NMDS1, y = NMDS2, colour = Treatment)) +
  geom_segment(data = seg,
               mapping = aes(xend = oNMDS1, yend = oNMDS2), linetype = "dashed") + # add spiders
  geom_point(data = centroids, size = 3) +                    # add centroids
  geom_point(size=1) +                                              
  #coord_fixed()+                                              
  theme_bw()+ 
  theme(legend.position="none",legend.background = element_rect(color = "black", linetype = "solid"),
        legend.text=element_text(size=12),legend.title= element_blank(),legend.direction='horizontal',
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))+
  scale_color_manual(values =pal1)+
  scale_y_continuous(limits = c(-0.5,0.6), n.breaks = 4)+
  scale_x_continuous(limits = c(-0.5,1), n.breaks = 4)
plotNMDSvegwet

final <- ggarrange(plotNMDSvegdry,plotNMDSvegwet, align="hv")

### tests wet
testwet <- adonis2(sqrt(veg) ~ Treatment, data=vegwet, permutations=9999, method="bray")
testwet
pairwet <- pairwise.adonis(sqrt(veg), vegwet$Treatment, perm = 9999)


### save figure
output_dir <- "OUTPUT"
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
file_path <- file.path(output_dir, "nmdsveg.png")
ggsave(file = file_path, plot = final, width = 8, height = 6, dpi = 300)  # Adjust width, height, and dpi as needed



##### Vegetation Shannon-diversity #####
Vegdiv <- as.data.frame(diversity(veg2020[,c(6:57)]))   # shannon default
#write.table(Fungaldiversity,file="fungaldiversity.txt")
#Merge with Veluwemetadata
Veluwemetadata=merge(veg2020, Vegdiv,
                     by = 'row.names', all = TRUE)

colnames(Veluwemetadata)[59] = "Vegdiv"

Veluwemetadata$Block = rep(c("A","B","C","D","E","F","G","H","I","J"), each =4)

Vdry = subset(Veluwemetadata, Site=="Dry")
Vwet = subset(Veluwemetadata, Site=="Wet")

mean(Vdry$Vegdiv[Vdry$Treatment=="Dolokal"])

Vegdivdry <- ggplot(aes(y = Vegdiv, x = Treatment, fill = Treatment), data = Vdry) + geom_boxplot() + theme_bw() + xlab("")+ ylab("")+
  scale_fill_manual(values=pal1)+theme(legend.position = "none", axis.text.x = element_blank(), axis.title.y = element_blank(), plot.background = element_rect("transparent"))+
  scale_y_continuous(n.breaks = 3, limits = c(1.5, 2.6))

mean(Vwet$Vegdiv[Vwet$Treatment=="Dolokal"],na.rm=T )

Vegdivwet <- ggplot(aes(y = Vegdiv, x = Treatment, fill = Treatment), data = Vwet) + geom_boxplot() + theme_bw() + xlab("")+ ylab("")+
  scale_fill_manual(values=pal1)+theme(legend.position = "none", axis.text.x = element_blank(), axis.title.y = element_blank(), plot.background = element_rect("transparent"))+
  scale_y_continuous(n.breaks = 3, limits = c(NA,1.5))

str(Vdry)


lm1 <- lme(Vegdiv~Treatment, ~1|Plot, data=Vdry)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "Treatment", data=Vdry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")
 
lm1 <- lme(Vegdiv~Treatment, ~1|Plot, data=Vwet)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "Treatment", data=Vwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")


finalvegdry=  plotNMDSvegdry +
  inset_element(Vegdivdry, align_to="plot",left = 0.6, bottom = 0.6, right = 0.99, top = 0.99)

finalvegwet=  plotNMDSvegwet +
  inset_element(Vegdivwet, align_to="plot",left = 0.6, bottom = 0.6, right = 0.99, top = 0.99)


##### Ecological groups ####

plants <- read.table("normalized_veg.txt", h=T)
plantssummed <- aggregate(. ~ Group, data = plants, sum)
plantssummed = as.data.frame(t(plantssummed))
colnames(plantssummed)=plantssummed[1,]
plantssummed = plantssummed[-1,]

# save for manipulation in excel
write.table(plantssummed, "Lifstyles_veg.txt")

## START analysis here

effects_veg <- read.table("effects_veg.txt", h=T)
str(effects_veg)
effects_veg$Treatment <- factor(effects_veg$Treatment , levels=c("Control", "Soilfeed", "Biolit","Dolokal"))
effects_veg$Site <- factor(effects_veg$Site, levels = c("Dry","Wet"))
effects_veg$plot = row.names(effects_veg)


effects <- effects_veg %>%
  pivot_longer(
    c(  "Moss",  "Woody","Ruderal","Grassland","Buffered","Heathland"),
    names_to="level",
    values_to="taxon"
  ) %>%
  mutate(
    level = factor(level,
                   levels=c("Moss","Woody","Ruderal","Grassland","Buffered","Heathland"),
                   labels = c("Mosses","Woody","Ruderal","Eutrophic grassland","Nardus grassland","Heathland"))
  )

### Dolokal

subset_data <- effects[effects$Treatment == "Dolokal", ]
sum_Do = summarySE(subset_data,
                   measurevar=c("taxon"),
                   groupvars=c("level", "Site"),
                   na.rm=TRUE)#in this last you could add more variables

sum_Dodry <- sum_Do[sum_Do$Site=="Dry",]
sum_DoWet <- sum_Do[sum_Do$Site=="Wet",]

### Control

subset_data <- effects[effects$Treatment == "Control", ]
sum_Co = summarySE(subset_data,
                   measurevar=c("taxon"),
                   groupvars=c("level", "Site"),
                   na.rm=TRUE)#in this last you could add more variables

sum_Codry <- sum_Co[sum_Co$Site=="Dry",]
sum_CoWet <- sum_Co[sum_Co$Site=="Wet",]

### Soilfeed

subset_data <- effects[effects$Treatment == "Soilfeed", ]
sum_So = summarySE(subset_data,
                   measurevar=c("taxon"),
                   groupvars=c("level", "Site"),
                   na.rm=TRUE)#in this last you could add more variables

sum_Sodry <- sum_So[sum_So$Site=="Dry",]
sum_SoWet <- sum_So[sum_So$Site=="Wet",]


### Biolit

subset_data <- effects[effects$Treatment == "Biolit", ]
sum_Bi = summarySE(subset_data,
                   measurevar=c("taxon"),
                   groupvars=c("level", "Site"),
                   na.rm=TRUE)#in this last you could add more variables

sum_Bidry <- sum_Bi[sum_Bi$Site=="Dry",]
sum_Biwet <- sum_Bi[sum_Bi$Site=="Wet",]


#### Final plot

soil_pointrange <- ggplot(data = sum_CoWet, aes(x = level, y = taxon)) +
  #geom_pointrange(data = sum_CoWet, aes(ymin = taxon - se, ymax = taxon + se), shape =16, size = 1.5, color = "#52854C") + 
  geom_pointrange(data = sum_DoWet, aes(ymin = taxon - se, ymax = taxon + se), shape=16, position = position_nudge(), size = 1, color = "#7B2D26") +  # Add another layer with different data
  geom_pointrange(data = sum_SoWet, aes(ymin = taxon - se, ymax = taxon + se), shape = 16, position = position_nudge(x= 0.4), size = 1,color = "#0E958F" )+
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
  #scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  scale_y_continuous(n.breaks = 4, limits = c(0,3))+
  coord_flip() + 
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black", linewidth = 0.5)
# geom_vline(aes(xintercept=level))

soil_pointrangedry <- ggplot(data = sum_Co, aes(x = level, y = taxon)) +
  #geom_pointrange(data = sum_COdry, aes(ymin = taxon - se, ymax = taxon + se), shape =17, size = 1.5, color = "#52854C") + 
  geom_pointrange(data = sum_Dodry, aes(ymin = taxon - se, ymax = taxon + se), shape=16, position = position_nudge(), size = 1, color = "#7B2D26") +  # Add another layer with different data
  geom_pointrange(data = sum_Bidry, aes(ymin = taxon - se, ymax = taxon + se), shape = 16,position = position_nudge(x= 0.2), size = 1, color = "#D16103") +  # Add another layer with different data
  geom_pointrange(data = sum_Sodry, aes(ymin = taxon - se, ymax = taxon + se), shape = 16, position = position_nudge(x= 0.4), size = 1,color = "#0E958F" ) +  # Add another layer with different data
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", fill = NA)
  ) + 
  theme(
    axis.title  = element_blank(),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size = 12)
  ) +
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  coord_flip() + 
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black", linewidth = 0.5)
# geom_vline(aes(xintercept=level))

vegdry <- soil_pointrangedry + scale_y_break(c(3.5,6), scales = 0.5)


final <- ggarrange(vegdry, soil_pointrange)

output_dir <- "OUTPUT"
#if (!file.exists(output_dir)) {
 # dir.create(output_dir)
#}
file_path <- file.path(output_dir, "effectsvegdry.png")
ggsave(file = file_path, plot = last_plot(), width = 6, height = 6, dpi = 300)  # Adjust width, height, and dpi as needed


#### TESTS

RA <- read.table("Lifstyles_veg.txt", h=T)
RA$beh <- factor(RA$beh , levels=c("Control", "Soilfeed", "Biolit","Dolokal"))
str(RA)
RAdry = subset(RA, Droog.Vochtig=="Droog")
RAwet = subset(RA, Droog.Vochtig=="Vochtig")


### dry site
# Buffered
lm1 <- lme(Buffered~beh, ~1|Block, data=RAdry)
anova.lme(lm1)
fit2.emm.a <- emmeans(lm1, "beh", data=RAdry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(RAdry$Buffered[RAdry$beh=="Control"])
sd(RAdry$Buffered[RAdry$beh=="Control"])

#Heathland 
lm1 <- lme(Heathland~beh, ~1|Block, data=RAdry)
anova.lme(lm1)
fit2.emm.a <- emmeans(lm1, "beh", data=RAdry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(RAdry$Heathland[RAdry$beh=="Control"])
sd(RAdry$Heathland[RAdry$beh=="Control"])

# Grassland  

lm1 <- lme(Grassland~beh, ~1|Block, data=RAdry)
anova.lme(lm1)
fit2.emm.a <- emmeans(lm1, "beh", data=RAdry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(RAdry$Grassland[RAdry$beh=="Biolit"])
sd(RAdry$Grassland[RAdry$beh=="Control"])

# Ruderals

lm1 <- lme(Ruderal~beh, ~1|Block, data=RAdry)
anova.lme(lm1)
fit2.emm.a <- emmeans(lm1, "beh", data=RAdry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(RAdry$Ruderal[RAdry$beh=="Biolit"])
sd(RAdry$Ruderal[RAdry$beh=="Control"])

# woody
lm1 <- lme(Woody~beh, ~1|Block, data=RAdry)
anova.lme(lm1)
fit2.emm.a <- emmeans(lm1, "beh", data=RAdry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")


mean(RAdry$Woody[RAdry$beh=="Control"])
sd(RAdry$Woody[RAdry$beh=="Control"])

# mosses

lm1 <- lme(Moss~beh, ~1|Block, data=RAdry)
anova.lme(lm1)
fit2.emm.a <- emmeans(lm1, "beh", data=RAdry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(RAdry$Moss[RAdry$beh=="Control"])
sd(RAdry$Moss[RAdry$beh=="Control"])



# wet site 
# Buffered
lm1 <- lme(Buffered~beh, ~1|Block, data=RAwet)
anova.lme(lm1)
fit2.emm.a <- emmeans(lm1, "beh", data=RAwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(RAwet$Arbuscularmycorrhizal)
sd(RAwet$Arbuscularmycorrhizal)

#Heathland 
lm1 <- lme(Heathland~beh, ~1|Block, data=RAwet)
anova.lme(lm1)
fit2.emm.a <- emmeans(lm1, "beh", data=RAdry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(RAwet$Heathland[RAwet$beh=="Control"])
sd(RAwet$Heathland[RAwet$beh=="Control"])


# woody

lm1 <- lme(Woody~beh, ~1|Block, data=RAwet)
anova.lme(lm1)
fit2.emm.a <- emmeans(lm1, "beh", data=RAwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(RAwet$Woody[RAwet$beh=="Control"])
sd(RAwet$Woody[RAwet$beh=="Control"])

# mosses
lm1 <- lme(Moss~beh, ~1|Block, data=RAwet)
anova.lme(lm1)
fit2.emm.a <- emmeans(lm1, "beh", data=RAwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")


mean(RAwet$Moss[RAwet$beh=="Biolit"])
sd(RAwet$Moss)

# unassigned

lm1 <- lme(unassigned~beh, ~1|Block, data=RAwet)
anova.lme(lm1)
fit2.emm.a <- emmeans(lm1, "beh", data=RAwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(RAwet$unassigned)
sd(RAwet$unassigned)

#### DISTANCE MATRIX (BC) ####
# VEG DRY
Vdry$Vegdiv = NULL
str(Vdry)
vegdry[, sapply(vegdry, is.numeric)] <- lapply(vegdry[, sapply(vegdry, is.numeric)], sqrt)

 #calculate distance for control vs treatments.
  bc_t1 <- vegdist(vegdry[vegdry$Treatment %in% c("Control", "Soilfeed"),6:57], method = "bray")
  bc_t2 <- vegdist(vegdry[vegdry$Treatment %in% c("Control", "Biolit"), 6:57], method = "bray")
  bc_t3 <- vegdist(vegdry[vegdry$Treatment %in% c("Control", "Dolokal"),6:57], method = "bray")
  bc_t4<- vegdist(vegdry[vegdry$Treatment=="Control",6:57], method = "bray")
 
  
  
   distance1<- data.frame(
    BC = c(bc_t1, bc_t2, bc_t3),
    Treatment = rep(c("Soilfeed", "Biolit", "Dolokal"), each = length(bc_t1))
  )

 control <- data.frame(BC=bc_t4, Treatment = "Control")
 control$eucl = eucl_t4
 
 distance1$eucl= c(eucl_t1, eucl_t2, eucl_t3) #eucl in soil script
 
 distance1 = rbind(distance1,control)
 dup_rows <- duplicated(distance1$BC) | duplicated(distance1$BC, fromLast = TRUE)
 distance2 <- distance1[distance1$Treatment == "Control" | !dup_rows, ]
 
 distance2$Treatment = factor(distance2$Treatment, levels = c("Control", "Soilfeed", "Biolit", "Dolokal"))
 
vegabiodry <- ggplot(aes(y=BC, x=eucl, color = Treatment), data=distance2)+
       geom_point(aes(size=1, alpha=.5))+
       scale_color_manual(values=pal1)+theme_bw()+
  geom_smooth(data = subset(distance2, Treatment == "Biolit"), method = "lm", se = FALSE, aes(group = Treatment)) +  # Only for a specific treatment
  geom_smooth(data = subset(distance2, Treatment == "Dolokal"), method = "lm", se = FALSE, aes(group = Treatment)) +  # Only for a specific treatment 
       xlab("Euclidean distance to control (abiotics)")+ ylab("BC dissimilarity to control (vegetation)")+theme(legend.position = "")

vegabiodry

### Manteltest
# control
mantel(bc_t4, eucl_t4, method = "spearman")
#Soilfeed
mantel(bc_t1, eucl_t1, method = "spearman")
#Biolit
mantel(bc_t2, eucl_t2, method = "spearman")
#Dolokal
mantel(bc_t3, eucl_t3, method = "spearman")


# VEG Wet

vegwet[, sapply(vegwet, is.numeric)] <- lapply(vegwet[, sapply(vegwet, is.numeric)], sqrt)

bc_t1 <- vegdist(vegwet[vegwet$Treatment %in% c("Control", "Soilfeed"),6:57], method = "bray")
bc_t2 <- vegdist(vegwet[vegwet$Treatment %in% c("Control", "Biolit"), 6:57], method = "bray")
bc_t3 <- vegdist(vegwet[vegwet$Treatment %in% c("Control", "Dolokal"),6:57], method = "bray")
bc_t4<- vegdist(vegwet[vegwet$Treatment=="Control",6:57], method = "bray")

distance1<- data.frame(
  BC = c(bc_t1, bc_t2, bc_t3),
  Treatment = rep(c("Soilfeed", "Biolit", "Dolokal"), each = length(bc_t1))
)

control <- data.frame(BC=bc_t4, Treatment = "Control")
control$eucl = eucl_t4

distance1$eucl= c(eucl_t1, eucl_t2, eucl_t3) #eucl in soil script

distance1 = rbind(distance1,control)
dup_rows <- duplicated(distance1$BC) | duplicated(distance1$BC, fromLast = TRUE)
distance2 <- distance1[distance1$Treatment == "Control" | !dup_rows, ]

str(distance2)
distance2$Treatment = factor(distance2$Treatment, levels = c("Control", "Soilfeed", "Biolit", "Dolokal"))

vegabiowet <- ggplot(aes(y=BC, x=eucl, color = Treatment), data=distance2)+
  geom_point(aes(size=1, alpha=.5))+
  scale_color_manual(values=pal1)+theme_bw()+
 geom_smooth(data = subset(distance2, Treatment == "Biolit"), method = "lm", se = FALSE, aes(group = Treatment), linetype="dashed") +  # Only for a specific treatment
  geom_smooth(data = subset(distance2, Treatment == "Dolokal"), method = "lm", se = FALSE, aes(group = Treatment)) +  # Only for a specific treatment
  xlab("Euclidean distance to control (abiotics)")+ ylab("")+theme(legend.position = "")

vegabiowet

### Manteltest
# control
mantel(bc_t4, eucl_t4, method = "spearman")
#Soilfeed
mantel(bc_t1, eucl_t1, method = "spearman")
#Biolit
mantel(bc_t2, eucl_t2, method = "spearman")
#Dolokal
mantel(bc_t3, eucl_t3, method = "spearman")