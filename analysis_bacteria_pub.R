###Multivariate analysis in VEGAN
install.packages('pheatmap')
library(pheatmap)
library("vegan")
library("rlang")
library("ggplot2")
library("ggpubr")
library(nlme)
#library("QsRutils")

pal1=c("#52854C","#0E958F", "#D16103","#7B2D26")

setwd("C:/Users/sfindeisen/Desktop/Heathland/VeluweData/data bacteria Veluwe")

#### SEQUENCE PREP ####
#Read ASV table
ASVtable_all <- read.table("ASVs_counts.tsv", header=TRUE)
colnames(ASVtable_all)=sub("_.*", "", colnames(ASVtable_all)) #replacement of _. in colnames

#Check the read nrs
readnr <- colSums(ASVtable_all)
readnr
#write.table(readnr,file="readnr.txt")

#######Data preparation########
#Read in tax data NOTE: some taxonomic levels have spaces in their names so can't be read into R easily
dat <- read.table("ASVs_taxonomy.tsv")
str(dat)

####FROM HERE ON, START REMOVING NOT BACTERIAL READS IN ASVTABLE, MITOCHONDRIA, etc.
#Data formatting => select all factors in the dataframe (sapply) and convert to characters (lapply)
i <- sapply(dat, is.factor)
dat[i] <- lapply(dat[i], as.character)
str(dat)
write.table(dat, file="taxformatted.txt")

#read the tax, then remove ASVs from the ASV table that were not matched to Bacteria earlier
tax1 <- read.table(as.matrix("taxformatted.txt"))
tax1=tax1[!is.na(tax1$domain),]
tax1=tax1[!tax1$family %in% c("Mitochondria"),]
tax1=tax1[!tax1$order %in% c("Chloroplast"),]
ASVtable <-  ASVtable_all[rownames(ASVtable_all) %in% rownames(tax1),] #(#4449 of 5629 left)

#Transpose matrix
ASVtransposed <- t(ASVtable) #transposes the matrix
length(ASVtransposed[1,]) #how many ASVs are there?
head(ASVtransposed)
#reads?

readnr3 <- colSums(t(ASVtransposed))
#write.table(readnr3,"readnr3.txt")
#readnrfinal <- read.table("readnr3.txt",h=T)

###Rarefaction => should be done per project seperately!
rarecurve(ASVtransposed, step = 500, xlab = "Sample Size", ylab = "Species", label = TRUE)

rarefied=rrarefy(ASVtransposed, sample= 17000) 
sum(colSums(rarefied)==0) ## how many are lost
ASVrarefiednozero =rarefied[,colSums(rarefied)>0]
length(ASVrarefiednozero[1,]) # how many are still there? from 4449 to 4431
rarecurve(ASVrarefiednozero, step = 500, xlab = "Sample Size", ylab = "Species", label = TRUE)

write.table(ASVrarefiednozero, file="rarefied.txt")

################################################################################################################################
######Analysis######
#Set session to source file location
####  Start here to load data
ASVmat=read.table(as.matrix("rarefied.txt"))

#ASVnrs
readnr_rarefied <- colSums(ASVmat)
readnr_Samples <- rowSums(ASVmat)
#write.table(readnr_rarefied,"readnr_rarefied.txt")
#write.table(readnr_Samples,"readnr_samples.txt")

##Load environmental data
env=read.table("env.txt", h=T)
rownames(env) <- paste("B", rownames(env), sep = "")
env=env[rownames(env) %in% rownames(ASVmat),]
ASVmat=ASVmat[rownames(ASVmat) %in% rownames(env),]
env=env[match(rownames(env),rownames(ASVmat)),]

#some samples that are double have been added to the env file. Among these, C22 seems very much off so removed from list

#Start here with clean files
ASV = read.table("ASVnormalized.txt",h=T)

##Load environmental data
env=read.table("VeluwemetadataBac.txt", h=T)
env$beh <- factor(env$beh , levels=c("Controle",  "Soilfeed","Biolit", "Dolokal"), labels= c("Control","Soilfeed","Biolit", "Dolokal"))
rownames(env)=env$Plotcode
env=env[rownames(env) %in% rownames(ASV),]
ASV=ASV[rownames(ASV) %in% rownames(env),]
env=env[match(rownames(env),rownames(ASV)),]

#Omit all the "empty" species (ASV) - 4185 left
ASV = ASV[, colSums(ASV) > 0]

####### DIVERSITY ##########
Bacdiv <- ggplot(aes(y = Diversity_prokaryotes, x = Droog.Vochtig, fill = beh), data = env) + geom_boxplot() + theme_bw() + xlab("")+ ylab("Shannon H diversity")+
  scale_fill_brewer(palette="Set2")+theme(legend.title=element_blank(),legend.text=element_text(size=13))+scale_x_discrete(labels=c("Dry","Wet"))
Bacdiv

#### Analysis for dry and wet site apart ####
envdroog = subset(env, Droog.Vochtig=="Droog")
envnat = subset(env, Droog.Vochtig=="Vochtig")

lm1 <- lme(Diversity_prokaryotes~beh, ~1|Plot, data=envdroog)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "beh", data=envdroog)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

lm1 <- lme(Diversity_prokaryotes~beh, ~1|Plot, data=envnat)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "beh", data=envnat)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

Bacdivdry <- ggplot(aes(y = Diversity_prokaryotes, x = beh, fill = beh), data = envdroog) + geom_boxplot() + theme_bw() + xlab("")+ ylab("Bacteriële diversiteit (Shannon H)")+
  scale_fill_manual(values=pal1)+theme(legend.position = "none", axis.text.x = element_blank(), axis.title.y = element_blank(), plot.background = element_rect("transparent"))+
  scale_y_continuous(breaks = c(5.0, 5.5, 6.0))
Bacdivdry

Bacdivwet <- ggplot(aes(y = Diversity_prokaryotes, x = beh, fill = beh), data = envnat) + geom_boxplot() + theme_bw() + xlab("")+ ylab("Bacteriële diversiteit (Shannon H)")+
  scale_fill_manual(values=pal1)+theme(legend.position="none", axis.title.y = element_blank(), axis.text.x = element_blank(), plot.background = element_rect("transparent"))+
  scale_y_continuous(limits =c(NA, 6), n.breaks = 3)
Bacdivwet

#Now first run fungal diversity script
ggarrange(Fungdiv1,Bacdiv1,align = "hv", common.legend=T)
ggarrange(Fungdiv2,Bacdiv2,align = "hv", common.legend=T)

Bacabund <- ggplot(aes(y = X16S_g, x = beh, fill = beh), data = envdroog) + geom_boxplot() + theme_bw() + xlab("")+ ylab("Bacteriële abundantie (aantal kopieën/g)")+
  scale_fill_manual(values=c("darkgreen","yellow","purple","red"))+theme(legend.title=element_blank(),legend.text=element_text(size=13))+scale_y_log10()
Bacabund 

Funabund <- ggplot(aes(y = FF_g, x = beh, fill = beh), data = envdroog) + geom_boxplot() + theme_bw() + xlab("")+ ylab("Schimmelabundantie (aantal kopieën/g)")+
  scale_fill_manual(values=c("darkgreen","yellow","purple","red"))+theme(legend.title=element_blank(),legend.text=element_text(size=13))+scale_y_log10()
Funabund

Funbacratio <- ggplot(aes(y = FF_16S_ratio, x = beh, fill = beh), data = envdroog) + geom_boxplot() + theme_bw() + xlab("")+ ylab("Schimmel:bacterie ratio")+
  scale_fill_manual(values=c("darkgreen","yellow","purple","red"))+theme(legend.title=element_blank(),legend.text=element_text(size=13))+ylim(0,0.8)
Funbacratio

ggarrange(Funabund,Bacabund,Funbacratio, align= "hv", ncol=3,common.legend=T)

#### Bacterial + Fungal abundance + FB ratio tests ####
lm1 <- lme(FF_g~beh, ~1|Plot, data=envnat)
anova.lme(lm1)
fit2.emm.a <- emmeans(lm1, "Treatment", data=Vdry)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

ylab = expression("Fungal:bacterial ratio")

Bacabund <- ggplot(aes(y = X16S_g, x = beh, fill = beh), data =env) + geom_boxplot() + theme_bw() + xlab("")+ ylab(ylab)+
  scale_fill_manual(values=pal1)+
  facet_wrap(Droog.Vochtig~., scales="fixed", dir = 'h', labeller = labeller(Droog.Vochtig = c(Vochtig = "Wet", Droog = "Dry")))+
  theme(legend.title=element_blank(),legend.position = "none")+scale_y_log10()
Bacabund 

mean(envdroog$FF_16S_ratio)

Funabund <- ggplot(aes(y = FF_g, x = beh, fill = beh), data =env) + geom_boxplot() + theme_bw() + xlab("")+ ylab(ylab)+
  scale_fill_manual(values=pal1)+
  facet_wrap(Droog.Vochtig~., scales="fixed", dir = 'h', labeller = labeller(Droog.Vochtig = c(Vochtig = "Wet", Droog = "Dry")))+
  theme(legend.title=element_blank(),legend.position = "none")+scale_y_log10()

Funbacratio <- ggplot(aes(y = FF_16S_ratio, x = beh, fill = beh), data =env) + geom_boxplot() + theme_bw() + xlab("")+ ylab(ylab)+
  scale_fill_manual(values=pal1)+
  facet_wrap(Droog.Vochtig~., scales="fixed", dir = 'h', labeller = labeller(Droog.Vochtig = c(Vochtig = "Wet", Droog = "Dry")))+
  theme(legend.title=element_blank(),legend.position = "none")+scale_y_log10()

ggarrange(Bacabund,Funabund,Funbacratio, align= "hv", ncol=2, nrow = 2)

######################################
##### NMDS for all data ###### 
ordr<-metaMDS(sqrt(ASV), distance = "bray", k = 2, binary=F)  # stress = 0.08045917
summary(ordr)
plot(ordr$points)
ordr$stress

#Plot using ggplot2
#Extract the axes scores
datascores=as.data.frame(scores(ordr)$sites) 

#Add/calculate spider diagram
scores <- cbind(as.data.frame(datascores), Treatment = env$Droog.Vochtig)
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Treatment, data = scores, FUN = mean)
seg <- merge(scores, setNames(centroids, c('Treatment','oNMDS1','oNMDS2')),
             by = 'Treatment', sort = FALSE)

#plot
plotNMDS1bac = ggplot(scores, aes(x = NMDS1, y = NMDS2, colour = Treatment)) +
  geom_segment(data = seg,
               mapping = aes(xend = oNMDS1, yend = oNMDS2), linewidth=.7) + # add spiders
  geom_point(data = centroids, size = 3) +                    # add centroids
  geom_point(size=1) +                                              
  coord_fixed()+                                              
  theme_bw()+ 
  theme(legend.position="top",legend.background = element_rect(color = "black", linetype = "solid"),
        legend.text=element_text(size=12),legend.title= element_blank(),legend.direction='horizontal',
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))+
  scale_color_manual(values = c("#4E84C4", "#293352"), labels= c("Dry", "Wet"))
plotNMDS1bac

#Add/calculate spider diagram
scores <- cbind(as.data.frame(datascores), Treatment = env$beh)
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Treatment, data = scores, FUN = mean)
seg <- merge(scores, setNames(centroids, c('Treatment','oNMDS1','oNMDS2')),
             by = 'Treatment', sort = FALSE)

#plot
plotNMDS2bac = ggplot(scores, aes(x = NMDS1, y = NMDS2, colour = Treatment)) +
  geom_segment(data = seg,
               mapping = aes(xend = oNMDS1, yend = oNMDS2), linewidth = .7) + # add spiders
  geom_point(data = centroids, size = 3) +                    # add centroids
  geom_point(size=1) +                                              
  coord_fixed()+                                              
  theme_bw()+ 
  theme(legend.position="top",legend.background = element_rect(color = "black", linetype = "solid"),
        legend.text=element_text(size=12),legend.title= element_blank(),legend.direction='horizontal',
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))+
  scale_color_manual(values =pal1, labels=c('Control', 'Biolit', 'Dolokal', 'Soilfeed'))
plotNMDS2bac

ggarrange(plotNMDS1bac,plotNMDS2bac, common.legend=F)

#Permanova
test1 <- adonis2(sqrt(ASV) ~ beh + Droog.Vochtig, data=env, permutations=999, method="bray")
test1

#######################################################################
####Dry and wet site apart ####

#Subset ASV data with env data from dry site
ASVdroog = subset(ASV , rownames(ASV) %in% rownames(envdroog))
###MATCH THE ORDER OF ENV DATA AND COMDATA = VERY IMPORTANT STEP!
ASVdroog <- ASVdroog[match(rownames(envdroog), rownames(ASVdroog)),]

#Subset ASV data with env data from wet
ASVnat = subset(ASV , rownames(ASV) %in% rownames(envnat))
###MATCH THE ORDER OF ENV DATA AND COMDATA = VERY IMPORTANT STEP!
ASVnat <- ASVnat[match(rownames(envnat), rownames(ASVnat)),]

####DRY SITE####
#NMDS
ordr<-metaMDS(sqrt(ASVdroog), distance = "bray", k = 2, binary=F)
summary(ordr)
plot(ordr$points)
ordr$stress    # stress = 0.1219795

#Plot using ggplot2
#Extract the axes scores
datascores=as.data.frame(scores(ordr)$sites) 

#Add/calculate spider diagram
scores <- cbind(as.data.frame(datascores), Treatment = envdroog$beh)
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Treatment, data = scores, FUN = mean)
seg <- merge(scores, setNames(centroids, c('Treatment','oNMDS1','oNMDS2')),
             byplotNMDS1bacdroog = 'Treatment', sort = FALSE)

#plot
plotNMDS1bacdroog= ggplot(scores, aes(x = NMDS1, y = NMDS2, colour = Treatment)) +
  geom_segment(data = seg,
               mapping = aes(xend = oNMDS1, yend = oNMDS2),linetype="dashed") + # add spiders
  geom_point(data = centroids, size = 3) +                    # add centroids
  geom_point(size=1, data=scores) +                                              
  #coord_fixed()+                                             
  theme_bw()+ 
  theme(legend.position = "none",axis.title = element_text(size = 12),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))+ scale_y_continuous(n.breaks = 3, limits = c(-0.5, 0.6))+
  scale_color_manual(values =pal1)

pal1=c("#52854C","#0E958F", "#D16103","#7B2D26")


finalbacdry=plotNMDS1bacdroog +
  inset_element(Bacdivdry, align_to="plot",left = 0.6, bottom = 0.6, right = 0.99, top = 0.99, clip = T)

final=ggarrange(finalbacdry,finalbacwet, finalfundry, finalfunwet,finalvegdry, finalvegwet, ncol = 2, nrow = 3)

file_path <- file.path(output_dir, "finalwithveg1.png")
ggsave(file = file_path, plot = final, width = 6, height = 8, dpi = 300)

#Permanova
test1 <- adonis2(sqrt(ASVdroog) ~ beh, data=envdroog, permutations=999, method="bray")
test1 #ns 

####WET SITE####
#NMDS 
ordr<-metaMDS(sqrt(ASVnat), distance = "bray", k = 2, binary=F)
summary(ordr)
plot(ordr$points)
ordr$stress #0.1270382

#Plot using ggplot2
#Extract the axes scores
datascores=as.data.frame(scores(ordr)$sites) 

#Add/calculate spider diagram
scores <- cbind(as.data.frame(datascores), Treatment = envnat$beh)
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Treatment, data = scores, FUN = mean)
seg <- merge(scores, setNames(centroids, c('Treatment','oNMDS1','oNMDS2')),
             by = 'Treatment', sort = FALSE)

#plot
plotNMDS2bacnat = ggplot(scores, aes(x = NMDS1, y = NMDS2, colour = Treatment)) +
  geom_segment(data = seg,
               mapping = aes(xend = oNMDS1, yend = oNMDS2), linetype="dashed") + # add spiders
  geom_point(data = centroids, size = 3) +                    # add centroids
  geom_point(size=1, data=scores) +                                              
  #coord_fixed()+                                              
  theme_bw()+ theme(legend.position = "none",axis.title = element_text(size = 12),plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))+
  scale_color_manual(values =pal1)+scale_y_continuous(n.breaks =4, limits = c(-0.6,0.8))

plotNMDS2bacnat

finalbacwet=  plotNMDS2bacnat +
  inset_element(Bacdivwet, align_to="plot",left = 0.6, bottom = 0.6, right = 0.99, top = 0.99)

final = ggarrange(plotNMDS1bacdroog,plotNMDS2bacnat, align="hv")

output_dir <- "OUTPUT"
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
file_path <- file.path(output_dir, "nmdsbacwet+dry.png")
ggsave(file = file_path, plot = final, width = 8, height = 6, dpi = 300)  # Adjust width, height, and dpi as needed


#Permanova
test1 <- adonis2(sqrt(ASVnat) ~ beh, data=envnat, permutations=9999, method="bray")
test1
pairwet <- pairwise.adonis(sqrt(ASVnat), envnat$beh, perm = 9999)


#### BC DISSIMILARITY DRY ####

ASVdroog[, sapply(ASVdroog, is.numeric)] <- lapply(ASVdroog[, sapply(ASVdroog, is.numeric)], sqrt)
ASVdroog$Treatment = envdroog$beh


#calculate distance for control vs treatments.
bc_t1 <- vegdist(ASVdroog[ASVdroog$Treatment %in% c("Control", "Soilfeed"),1:4200], method = "bray")
bc_t2 <- vegdist(ASVdroog[ASVdroog$Treatment %in% c("Control", "Biolit"), 1:4200], method = "bray")
bc_t3 <- vegdist(ASVdroog[ASVdroog$Treatment %in% c("Control", "Dolokal"),1:4200], method = "bray")
bc_t4<- vegdist(ASVdroog[ASVdroog$Treatment=="Control",1:4200], method = "bray")



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

bacabiodry <- ggplot(aes(y=BC, x=eucl, color = Treatment), data=distance2)+
  geom_point(aes(size=1, alpha=.5))+
  scale_color_manual(values=pal1)+theme_bw()+
  geom_smooth(data = subset(distance2, Treatment == "Dolokal"), method = "lm", se = FALSE, aes(group = Treatment)) +  # Only for a specific treatment
  xlab("")+ ylab("BC dissimilarity to control (prokaryotes)")+theme(legend.position = "")

bacabiodry

### Manteltest
# control
mantel(bc_t4, eucl_t4, method = "spearman")
#Soilfeed
mantel(bc_t1, eucl_t1, method = "spearman")
#Biolit
mantel(bc_t2, eucl_t2, method = "spearman")
#Dolokal
mantel(bc_t3, eucl_t3, method = "spearman")


#### BC DISSIMILARITY WET ####

ASVnat[, sapply(ASVnat, is.numeric)] <- lapply(ASVnat[, sapply(ASVnat, is.numeric)], sqrt)
ASVnat$Treatment = envnat$beh


bc_t1 <- vegdist(ASVnat[ASVnat$Treatment %in% c("Control", "Soilfeed"),1:4200], method = "bray")
bc_t2 <- vegdist(ASVnat[ASVnat$Treatment %in% c("Control", "Biolit"), 1:4200], method = "bray")
bc_t3 <- vegdist(ASVnat[ASVnat$Treatment %in% c("Control", "Dolokal"),1:4200], method = "bray")
bc_t4<- vegdist(ASVnat[ASVnat$Treatment=="Control",1:4200], method = "bray")

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

bacabiowet <- ggplot(aes(y=BC, x=eucl, color = Treatment), data=distance2)+
  geom_point(aes(size=1, alpha=.5))+
  scale_color_manual(values=pal1)+theme_bw()+
  # geom_smooth(data = subset(distance2, Treatment == "Biolit"), method = "lm", se = FALSE, aes(group = Treatment), linetype="dashed") +  # Only for a specific treatment
  geom_smooth(data = subset(distance2, Treatment == "Soilfeed"), method = "lm", se = FALSE, aes(group = Treatment)) + 
  geom_smooth(data = subset(distance2, Treatment == "Biolit"), method = "lm", se = FALSE, aes(group = Treatment)) + 
  geom_smooth(data = subset(distance2, Treatment == "Dolokal"), method = "lm", se = FALSE, aes(group = Treatment)) +  # Only for a specific treatment
  # Only for a specific treatment
  # Only for a specific treatment
  xlab("")+ ylab("")+theme(legend.position = "")

bacabiowet

### Manteltest
# control
mantel(bc_t4, eucl_t4, method = "spearman")
#Soilfeed
mantel(bc_t1, eucl_t1, method = "spearman")
#Biolit
mantel(bc_t2, eucl_t2, method = "spearman")
#Dolokal
mantel(bc_t3, eucl_t3, method = "spearman")