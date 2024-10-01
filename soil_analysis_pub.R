library("ggplot2")
library(ggpubr)
library(DescTools)
library(FSA)
library(lme4)
library(nlme)
library(emmeans)
library(sidier)

setwd("C:/Users/sfindeisen/Desktop/Heathland/VeluweData/data soil Veluwe")

pal1=c("#52854C","#0E958F", "#D16103","#7B2D26")

### Read data ####
soil <- read.table("2020soilchemheathland.txt", h=T)
rownames(soil)= soil[,1]
soil$treatment <- factor(soil$treatment, levels = c("controle", "Soilfeed","biolit", "dolokal"), labels = c("Control", "Soilfeed", "Biolit", "Dolokal"))
soil$Block = as.factor(soil$Block)
soil$X.BS = as.numeric(soil$X.BS)
soil$CEC_meq.L = as.numeric(soil$CEC_meq.L)
soil$total_BC.Ca_K_Mg.= as.numeric(soil$total_BC.Ca_K_Mg.)

#### Wet and dry site apart ####
soildry = subset(soil, Dry_wet == "droog")
soilwet = subset(soil, Dry_wet == "nat")

# SOM
lm1 <- lme(OM~treatment, ~1|Block, data=soildry)
anova.lme(lm1)
summary(lm1)

mean(soildry$OM[soildry$treatment=="Control"])
sd(soildry$OM[soildry$treatment=="Control"])

fit2.emm.a <- emmeans(lm1, "treatment", data=soildry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

lm1 <- lme(OM~treatment, ~1|Block, data=soilwet)
anova.lme(lm1)
summary(lm1)

mean(soilwet$OM[soilwet$treatment=="Control"])
sd(soilwet$OM[soilwet$treatment=="Control"])

# POlsen
ylab=expression("P Olsen μmol kg soil" ^ -1)

p1 <- ggplot(aes(y = Olsen.P, x =treatment, fill = treatment), data = soil) + geom_boxplot() + theme_bw() + xlab("")+ ylab(ylab)+
  scale_fill_manual(values=pal1)+theme_bw()+
  facet_wrap(Dry_wet~., dir = 'h', labeller = labeller(Dry_wet = c(nat = "Wet", droog = "Dry")))+
  theme(legend.title=element_blank(),legend.position = "none")
p1

shapiro.test(soildry$Olsen.P)

lm1 <- lme(Olsen.P~treatment, ~1|Block, data=soildry)
anova.lme(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soildry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soildry$Olsen.P[soildry$treatment=="Control"])
sd(soildry$Olsen.P[soildry$treatment=="Control"])

anova.lme(lm1)
summary(lm1)

lm1 <- lme(Olsen.P~treatment, ~1|Block, data=soilwet)
anova.lme(lm1)
summary(lm1)
mean(soilwet$Olsen.P[soilwet$treatment=="Control"])
sd(soilwet$Olsen.P[soilwet$treatment=="Control"])
 
# NH4
shapiro.test(log(soil$NH4_NaCl))
ylab=expression(paste("log NH"[4], " μmol kg soil" ^-1))

p2 <- ggplot(aes(y = NH4_NaCl, x =treatment, fill = treatment), data = soil) + geom_boxplot() + theme_bw() + xlab("")+ ylab(ylab)+
  scale_fill_manual(values=pal1)+theme_bw()+
  facet_wrap(Dry_wet~., scales="free_y", dir = 'h', labeller = labeller(Dry_wet = c(nat = "Wet", droog = "Dry")))+
  theme(legend.title=element_blank(),legend.position = "none")
  
p2

lm1 <- lme(NH4_NaCl~treatment, ~1|Block, data=soildry)
anova.lme(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soildry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soildry$NH4_NaCl[soildry$treatment=="Control"])
sd(soildry$NH4_NaCl[soildry$treatment=="Control"])

lm1 <- lme(NH4_NaCl~treatment, ~1|Block, data=soilwet)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soilwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soilwet$NH4_NaCl[soilwet$treatment=="Control"])
sd(soilwet$NH4_NaCl[soilwet$treatment=="Control"])

# NO3
shapiro.test(soil$NO3_NaCl)
ylab=expression(paste("NO"[3] ^-"", " μmol kg soil" ^-1))

p2 <- ggplot(aes(y = NO3_NaCl, x =treatment, fill = treatment), data = soil) + geom_boxplot() + theme_bw() + xlab("")+ ylab(ylab)+
  scale_fill_manual(values=pal1)+theme_bw()+
  facet_wrap(Dry_wet~., scales="free_y", dir = 'h', labeller = labeller(Dry_wet = c(nat = "Wet", droog = "Dry")))+
  theme(legend.title=element_blank(),legend.position = "none")

p2

lm1 <- lme(NO3_NaCl~treatment, ~1|Block, data=soildry)
anova.lme(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soildry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soildry$NO3_NaCl[soildry$treatment=="Control"])
sd(soildry$NO3_NaCl[soildry$treatment=="Control"])

lm1 <- lme(NO3_NaCl~treatment, ~1|Block, data=soilwet)
anova.lme(lm1)
summary(lm1)

mean(soilwet$NO3_NaCl[soilwet$treatment=="Control"])
sd(soilwet$NO3_NaCl[soilwet$treatment=="Control"])

# BS
shapiro.test(soil$X.BS)
ylab=expression("Base saturation (%)")

p3 <- ggplot(aes(y =X.BS, x =treatment, fill = treatment), data = soil) + geom_boxplot(alpha=0.5) + theme_bw() + xlab("")+ ylab(ylab)+
 geom_point(aes(y =X.BS, x =treatment, fill = treatment), data = soil)+ scale_fill_manual(values=pal1)+theme_bw()+
  facet_wrap( ~ Dry_wet, dir = 'h', labeller = labeller(Dry_wet = c(nat = "Wet", droog = "Dry")))+
  theme(legend.title=element_blank(),legend.position = "none")
p3

lm1 <- lme(X.BS~treatment, ~1|Block, data=soildry)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soildry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soildry$X.BS[soildry$treatment=="Dolokal"])
sd(soildry$X.BS[soildry$treatment=="Control"])

lm1 <- lme(X.BS~treatment, ~1|Block, data=soilwet)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soilwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soilwet$X.BS[soilwet$treatment=="Dolokal"])
sd(soilwet$X.BS[soilwet$treatment=="Control"])


# pH NaCl
shapiro.test(soil$pH.NaCl)
ylab=expression("pH (NaCl)")

p4 <- ggplot(aes(y = pH.NaCl, x =treatment, fill = treatment), data = soil) + geom_boxplot() + theme_bw() + xlab("")+ ylab(ylab)+
  scale_fill_manual(values=pal1)+theme_bw()+
  facet_wrap(Dry_wet~., dir = 'h', labeller = labeller(Dry_wet = c(nat = "Wet", droog = "Dry")))+
  theme(legend.title=element_blank(),legend.position = "none")
       
p4

final = ggarrange(p3, p4, nrow = 2)


lm1 <- lme(pH.NaCl~treatment, ~1|Block, data=soildry)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soildry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soildry$pH.NaCl[soildry$treatment=="Dolokal"])
sd(soildry$pH.NaCl[soildry$treatment=="Control"])


lm1 <- lme(pH.NaCl~treatment, ~1|Block, data=soilwet)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soilwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soilwet$pH.NaCl[soilwet$treatment=="Dolokal"])
sd(soilwet$pH.NaCl[soilwet$treatment=="Control"])


ptogether = ggarrange(p3,p4, ncol = 1, nrow = 2)
output_dir <- "OUTPUT"
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
file_path <- file.path(output_dir, "soilabiotics.png")
ggsave(file = file_path, plot = final, width = 6, height = 6, dpi = 300)  # Adjust width, height, and dpi as needed



#Al:Ca
ylab=expression("Al/Ca ratio")
shapiro.test(log(soil$Al_Ca_NaCl))

p4 <- ggplot(aes(y = Al_Ca_NaCl, x =treatment, fill = treatment), data = soil) + geom_boxplot() + theme_bw() + xlab("")+ ylab(ylab)+
  scale_fill_manual(values=pal1)+theme_bw()+
  facet_wrap(Dry_wet~., dir = 'h')+
  theme(legend.title=element_blank(),legend.position = "none")
 # geom_hline(yintercept = 2, linetype = "dashed", colour = "black", linewidth = 0.5)

p4

lm1 <- lme(Al_Ca_NaCl~treatment, ~1|Block, data=soildry)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soildry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soildry$Al_Ca_NaCl[soildry$treatment=="Control"])
sd(soildry$Al_Ca_NaCl[soildry$treatment=="Control"])

lm1 <- lme(Al_Ca_NaCl~treatment, ~1|Block, data=soilwet)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soilwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soilwet$Al_Ca_NaCl[soilwet$treatment=="Control"])
sd(soilwet$Al_Ca_NaCl[soilwet$treatment=="Control"])

# C:N
shapiro.test(log(soil$C.N.Ratio))
ylab=expression("C/N ratio")

p5 <- ggplot(aes(y = C.N.Ratio, x =treatment, fill = treatment), data = soil) + geom_boxplot() + theme_bw() + xlab("")+ ylab(ylab)+
  scale_fill_manual(values=pal1)+theme_bw()+
  facet_wrap(Dry_wet~., dir = 'h', labeller = labeller(Dry_wet = c(nat = "Wet", droog = "Dry")))+
  theme(legend.title=element_blank(),legend.position = "none")
p5

lm1 <- lme(C.N.Ratio~treatment, ~1|Block, data=soildry)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soildry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soildry$C.N.Ratio[soildry$treatment=="Control"])
sd(soildry$C.N.Ratio[soildry$treatment=="Control"])

lm1 <- lme(C.N.Ratio~treatment, ~1|Block, data=soilwet)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soilwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soilwet$C.N.Ratio[soilwet$treatment=="Control"])
sd(soilwet$C.N.Ratio[soilwet$treatment=="Control"])

# Fe

lm1 <- lme(Fe_NaCl~treatment, ~1|Block, data=soildry)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soildry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soildry$Fe_NaCl[soildry$treatment=="Control"])
sd(soildry$Fe_NaCl[soildry$treatment=="Control"])

lm1 <- lme(Fe_NaCl~treatment, ~1|Block, data=soilwet)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soilwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soilwet$Fe_NaCl[soilwet$treatment=="Control"])
sd(soilwet$Fe_NaCl[soilwet$treatment=="Control"])

ylab=expression("Fe")

p5 <- ggplot(aes(y = Mg_NaCl, x =treatment, fill = treatment), data = soil) + geom_boxplot() + theme_bw() + xlab("")+ ylab("K")+
  scale_fill_manual(values=pal1)+theme_bw()+
  facet_wrap(Dry_wet~., dir = 'h', labeller = labeller(Dry_wet = c(nat = "Wet", droog = "Dry")))+
  theme(legend.title=element_blank(),legend.position = "none")
p5

# K
lm1 <- lme(K_NaCl~treatment, ~1|Block, data=soildry)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soildry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soildry$K_NaCl[soildry$treatment=="Control"])
sd(soildry$K_NaCl[soildry$treatment=="Control"])

lm1 <- lme(K_NaCl~treatment, ~1|Block, data=soilwet)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soilwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soilwet$K_NaCl[soilwet$treatment=="Control"])
sd(soilwet$K_NaCl[soilwet$treatment=="Control"])

# Mg
lm1 <- lme(Mg_NaCl~treatment, ~1|Block, data=soildry)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soildry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soildry$Mg_NaCl[soildry$treatment=="Control"])
sd(soildry$Mg_NaCl[soildry$treatment=="Control"])

lm1 <- lme(Mg_NaCl~treatment, ~1|Block, data=soilwet)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soilwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soilwet$Mg_NaCl[soilwet$treatment=="Control"])
sd(soilwet$Mg_NaCl[soilwet$treatment=="Control"])

# Mn

lm1 <- lme(Mn_NaCl~treatment, ~1|Block, data=soildry)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soildry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soildry$Mn_NaCl[soildry$treatment=="Control"])
sd(soildry$Mn_NaCl[soildry$treatment=="Control"])

lm1 <- lme(Mn_NaCl~treatment, ~1|Block, data=soilwet)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soilwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soilwet$Mn_NaCl[soilwet$treatment=="Control"])
sd(soilwet$Mn_NaCl[soilwet$treatment=="Control"])

# S
lm1 <- lme(S_NaCl~treatment, ~1|Block, data=soildry)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soildry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soildry$S_NaCl[soildry$treatment=="Control"])
sd(soildry$S_NaCl[soildry$treatment=="Control"])

lm1 <- lme(S_NaCl~treatment, ~1|Block, data=soilwet)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soilwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soilwet$S_NaCl[soilwet$treatment=="Control"])
sd(soilwet$S_NaCl[soilwet$treatment=="Control"])

# Si
lm1 <- lme(Si_NaCl~treatment, ~1|Block, data=soildry)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soildry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soildry$Si_NaCl[soildry$treatment=="Control"])
sd(soildry$Si_NaCl[soildry$treatment=="Control"])

lm1 <- lme(Si_NaCl~treatment, ~1|Block, data=soilwet)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soilwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soilwet$Si_NaCl[soilwet$treatment=="Control"])
sd(soilwet$Si_NaCl[soilwet$treatment=="Control"])

# Zn
lm1 <- lme(Zn_NaCl~treatment, ~1|Block, data=soildry)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soildry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soildry$Zn_NaCl[soildry$treatment=="Control"])
sd(soildry$Zn_NaCl[soildry$treatment=="Control"])

lm1 <- lme(Zn_NaCl~treatment, ~1|Block, data=soilwet)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "treatment", data=soilwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soilwet$Zn_NaCl[soilwet$treatment=="Control"])
sd(soilwet$Zn_NaCl[soilwet$treatment=="Control"])

# Al
p0 <- ggplot(aes(y = Al_NaCl, x =treatment, fill = treatment), data = soil) + geom_boxplot(alpha=0.5) + theme_bw() + xlab("")+ ylab("Organic matter content (%)")+
  scale_fill_manual(values=pal1)+theme_bw()+
  facet_wrap(Dry_wet~.,scales = "free_y", 
             dir = 'h', labeller = labeller(Dry_wet = c(nat = "Wet", droog = "Dry")))+
  theme(legend.title=element_blank(),legend.position = "none")
p0

lm1 <- lme(Al_NaCl~treatment, ~1|Block, data=soildry)
anova.lme(lm1)
summary(lm1)

mean(soilwet$Al_NaCl[soilwet$treatment=="Control"])
sd(soilwet$Al_NaCl[soilwet$treatment=="Control"])

fit2.emm.a <- emmeans(lm1, "treatment", data=soildry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

lm1 <- lme(Al_NaCl~treatment, ~1|Block, data=soilwet)
anova.lme(lm1)
summary(lm1)

fit2.emm.a <- emmeans(lm1, "treatment", data=soilwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soilwet$OM[soilwet$treatment=="Control"])
sd(soilwet$OM[soilwet$treatment=="Control"])

# Ca
p0 <- ggplot(aes(y = Ca_NaCl, x =treatment, fill = treatment), data = soil) + geom_boxplot(alpha=0.5) + theme_bw() + xlab("")+ ylab("Organic matter content (%)")+
  scale_fill_manual(values=pal1)+theme_bw()+
  facet_wrap(Dry_wet~.,scales = "free_y", 
             dir = 'h', labeller = labeller(Dry_wet = c(nat = "Wet", droog = "Dry")))+
  theme(legend.title=element_blank(),legend.position = "none")
p0

lm1 <- lme(Ca_NaCl~treatment, ~1|Block, data=soildry)
anova.lme(lm1)
summary(lm1)

mean(soildry$Ca_NaCl[soildry$treatment=="Control"])
sd(soildry$Ca_NaCl[soildry$treatment=="Control"])

fit2.emm.a <- emmeans(lm1, "treatment", data=soildry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

lm1 <- lme(Ca_NaCl~treatment, ~1|Block, data=soilwet)
anova.lme(lm1)
summary(lm1)

fit2.emm.a <- emmeans(lm1, "treatment", data=soilwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

mean(soilwet$Ca_NaCl[soilwet$treatment=="Control"])
sd(soilwet$Ca_NaCl[soilwet$treatment=="Control"])

# CEC
p0 <- ggplot(aes(y = CEC_meq.L, x =treatment, fill = treatment), data = soil) + geom_boxplot(alpha=0.5) + theme_bw() + xlab("")+ ylab("Organic matter content (%)")+
  scale_fill_manual(values=pal1)+theme_bw()+
  facet_wrap(Dry_wet~.,scales = "free_y", 
             dir = 'h', labeller = labeller(Dry_wet = c(nat = "Wet", droog = "Dry")))+
  theme(legend.title=element_blank(),legend.position = "none")
p0

lm1 <- lme(CEC_meq.L~treatment, ~1|Block, data=soildry)
anova.lme(lm1)
summary(lm1)

mean(soildry$CEC_meq.L[soildry$treatment=="Control"])
sd(soildry$CEC_meq.L[soildry$treatment=="Control"])

fit2.emm.a <- emmeans(lm1, "treatment", data=soilwet)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

lm1 <- lme(CEC_meq.L~treatment, ~1|Block, data=soilwet)
anova.lme(lm1)
summary(lm1)

mean(soilwet$CEC_meq.L[soilwet$treatment=="Control"])
sd(soilwet$CEC_meq.L[soilwet$treatment=="Control"])

#### EUCLIDIAN DISTANCE ####

# Dry subset
soildry = soildry[, -c(5:7)]
scalesoildry = as.data.frame(scale(soildry[,5:54]))
scalesoildry = cbind(scalesoildry, soildry[,1:4])
scalesoildry = scalesoildry[,-c(3:12,15,17,19,22,24,33:45)]# Only use available nutrient concentractions (NaCl extractions) 

scalesoildry = subset(scalesoildry , rownames(scalesoildry) %in% rownames(ASVdroog)) # only for fungi
  
  #calculate distance for control vs treatments.
  eucl_t1 <- dist(scalesoildry[scalesoildry$treatment %in% c("Control", "Soilfeed"),1:22], "euclidean")
  eucl_t2 <- dist(scalesoildry[scalesoildry$treatment %in% c("Control", "Biolit"),1:22], "euclidean")
  eucl_t3 <- dist(scalesoildry[scalesoildry$treatment %in% c("Control", "Dolokal"),1:22], "euclidean")
  eucl_t4 <- dist(scalesoildry[scalesoildry$treatment =="Control", 1:22], "euclidean")
  

#p1 <- ggplot(aes(y = euclidian_distance, x =treatment, fill = treatment), data = eucl_df) + geom_boxplot() + theme_bw() + xlab("")+ ylab("euclidian distance")+
 # scale_fill_manual(values=pal1)+theme_bw()


  # Wet subset
  soilwet = soilwet[, -c(5:7)]
  scalesoilwet = as.data.frame(scale(soilwet[,5:54]))
  scalesoilwet = cbind(scalesoilwet, soilwet[,1:4])
  scalesoilwet = scalesoilwet[,-c(3:12,15,17,19,22,24,33:45)]# Only use available nutrient concentractions (NaCl extractions) 
  
 
  
  #calculate distance for control vs treatments.
  eucl_t1 <- dist(scalesoilwet[scalesoilwet$treatment %in% c("Control", "Soilfeed"),1:22], "euclidean")
  eucl_t2 <- dist(scalesoilwet[scalesoilwet$treatment %in% c("Control", "Biolit"),1:22], "euclidean")
  eucl_t3 <- dist(scalesoilwet[scalesoilwet$treatment %in% c("Control", "Dolokal"),1:22], "euclidean")
  eucl_t4 <- dist(scalesoilwet[scalesoilwet$treatment =="Control", 1:22], "euclidean")
  