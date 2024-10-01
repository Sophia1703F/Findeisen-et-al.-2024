library("dplyr")
library("vegan")
library("rlang")
library("ggplot2")
library("RVAideMemoire")
library("ggpubr")
library("DescTools")
library(patchwork)
library(reshape2)
library(cowplot)
library(patchwork)
library(tidyr)


setwd("C:/Users/sfindeisen/Desktop/Heathland/VeluweData/data fungi Veluwe")

pal1=c("#52854C","#0E958F", "#D16103","#7B2D26","#4E84C4")

#### Sequence PREP ####
#Read AV table
ASVtable_all <- read.table("ASVs_counts.tsv", header=TRUE)
colnames(ASVtable_all)=sub("_.*", "", colnames(ASVtable_all))

#Check the read nrs
readnr <- colSums(ASVtable_all)
readnr
write.table(readnr,file="readnr.txt")

#Read in tax data
dat <- read.table("TaxonomyLifestyle.txt")
str(dat)

#Data formatting => select all factors in the dataframe (sapply) and convert to characters (lapply)
i <- sapply(dat, is.factor)
dat[i] <- lapply(dat[i], as.character)
str(dat)
write.table(dat, file="taxformatted.txt")

#read the ASV table, then remove ASVs that were not matched to UNITE earlier
tax1 <- read.table(as.matrix("taxformatted.txt"))
ASVtable <-  subset(ASVtable_all, rownames(ASVtable_all) %in% rownames(tax1)) 

#Number of ASVs
ASVtransposed <- t(ASVtable) 
length(ASVtransposed[1,])

#Number of reads
readnr3 <- colSums(t(ASVtransposed))
write.table(readnr3,"readnr3.txt")
readnrfinal <- read.table("readnr3.txt",h=T)

# WATCH OUT: sample FE2 has only 25 reads and should be removed
ASVtransposed <- ASVtransposed[!(rownames(ASVtransposed) %in% c("E01","E50", "E53", "E53.2", "E55", "E57.1", "H51.1", "H52.1", "H54.1", "H54.2", "H55", "H55.2", "L03.1", "L03","L06", "L09", "L07.1", "L12.1", "L12", "L68", "S47", "S49", "FE2", "FN01", "FN02", "FN03",  "FN04", "FN1", "FN2", "FN3", "FN4")),]

###Rarefaction => should be done per project seperately!
rarecurve(ASVtransposed, step = 500, xlab = "Sample Size", ylab = "Species", label = TRUE)
rarefied=rrarefy(ASVtransposed, sample= 10000)
sum(colSums(rarefied)==0) ## how many are lost: 367 
ASVrarefiednozero =rarefied[,colSums(rarefied)>0]
length(ASVrarefiednozero[1,]) # how many are still there? 7773 ASV's 
rarecurve(ASVrarefiednozero, step = 500, xlab = "Sample Size", ylab = "Species", label = TRUE)

write.table(ASVrarefiednozero, file="rarefied_10000.txt")

## clean files to read in later ##

ASVmat=read.table("rarefied_10000.txt")

#ASVnrs
readnr_rarefied <- colSums(ASVmat)
readnr_Samples <- rowSums(ASVmat)

##Load environmental data
env=read.table("envVeluweHeischraal.txt", h=T)
env$newnames=rownames(env)
env$newnames[1:40] <- paste("F", rownames(env[1:40,]), sep = "")
rownames(env)=env$newnames
test=setdiff(rownames(env), rownames(ASVmat)) #which names are in env that are not in ASVmat?
#subset veluwe project
env=env[env$Type=="Veluwe",]
env=env[rownames(env) %in% rownames(ASVmat),]
ASVmat=ASVmat[rownames(ASVmat) %in% rownames(env),]
env=env[order(match(rownames(env),rownames(ASVmat))),]

ASVmat =ASVmat[,colSums(ASVmat)>0]

#Normalize data because some plots < 10000 reads
ASV<- t(apply(ASVmat,1, function(x) x/sum(x, na.rm=TRUE)))#Normalize the data to sample (= relative gene abundance per sample)
rowSums(ASV)
write.table(ASV, "ASVnormalized.txt")

ASV = t(ASV)
ASV <- as.data.frame(ASV)
tax1=subset(tax1, rownames(tax1) %in% rownames(ASV))

write.table(tax1, "taxsubsetted.txt")

###########################################
#Start here with clean files
####NMDS analysis####
###########################################
#Veluwemetadata = subset(env, Type=="Veluwe")
#Remove row nr5
#Veluwemetadata <- Veluwemetadata[-c(5), ]
#write.table(Veluwemetadata, "Veluwemetadata.txt")

ASV = read.table("ASVnormalized.txt",h=T)
Veluwemetadata=read.table("Veluwemetadata.txt",h=T)
Veluwemetadata$Treatment_combi = paste(Veluwemetadata$Droog.Vochtig,Veluwemetadata$beh)
Veluwemetadata$beh  = factor(Veluwemetadata$beh, c("Controle",  "Soilfeed", "Biolit","Dolokal"), labels = c("Control", "Soilfeed", "Biolit","Dolokal" ))

#Subset env data with vegetation data (subsetting if needed)
ASV = subset(ASV , rownames(ASV) %in% rownames(Veluwemetadata))
###MATCH THE ORDER OF ENV DATA AND COMDATA = VERY IMPORTANT STEP!
ASV <- ASV[match(rownames(Veluwemetadata), rownames(ASV)),]
#Omit all the "empty" species
ASV = ASV[, colSums(ASV) > 0]

#Calculate fungal Shannon-diversity 
Fungaldiversity <- as.data.frame(diversity(ASV))   # shannon default
#write.table(Fungaldiversity,file="fungaldiversity.txt")
#Merge with Veluwemetadata
Veluwemetadata=merge(Veluwemetadata, Fungaldiversity,
      by = 'row.names', all = TRUE)
rownames(Veluwemetadata)=Veluwemetadata[,1]

#Change column name
colnames(Veluwemetadata)[14] ="Fungal_diversity"

Veluwemetadatadroog = subset(Veluwemetadata, Droog.Vochtig=="Droog")
rownames(Veluwemetadatadroog) = Veluwemetadatadroog$Row.names
Veluwemetadatanat = subset(Veluwemetadata, Droog.Vochtig=="Vochtig")
rownames(Veluwemetadatanat) = Veluwemetadatanat$Row.names

#Diversityplot
Fungdiv <- ggplot(aes(y = Fungal_diversity, x = Droog.Vochtig, fill = beh), data = Veluwemetadata) + geom_boxplot() + theme_bw() + xlab("")+ ylab("Schimmeldiversiteit (Shannon H)")+
  scale_fill_manual(values=c(pal1))+theme(legend.title=element_blank(),legend.text=element_text(size=13))
Fungdiv

#Dry and wet site seperated
Fungdivdry <- ggplot(aes(y = Fungal_diversity, x = beh, fill = beh), data = Veluwemetadatadroog) + geom_boxplot() + theme_bw() + xlab("")+ ylab("Fungal diversity (Shannon H)")+
  scale_fill_manual(values=c(pal1))+theme(legend.title=element_blank(),legend.position = "none", axis.title.y = element_blank(), axis.text.x = element_blank(), plot.background = element_rect("transparent"))
Fungdivdry

Fungdivwet <- ggplot(aes(y = Fungal_diversity, x = beh, fill = beh), data = Veluwemetadatanat) + geom_boxplot() + theme_bw() + xlab("")+ ylab("Fungal diversity (Shannon H)")+
  scale_fill_manual(values=c(pal1))+theme(legend.title=element_blank(),legend.position = "none", axis.title.y = element_blank(), axis.text.x = element_blank(), plot.background = element_rect("transparent"))
Fungdivwet

### Tests diversity
#Dry site
lm1 <- lme(Fungal_diversity~beh, ~1|Plot, data=Veluwemetadatadroog)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "beh", data=Veluwemetadatadroog)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

#Wet site
lm1 <- lme(Fungal_diversity~beh, ~1|Plot, data=Veluwemetadatanat)
anova.lme(lm1)
summary(lm1)
fit2.emm.a <- emmeans(lm1, "beh", data=Veluwemetadatanat)          #multiple comparisons with bonferoni adj.
pairs(fit2.emm.a, adjust="bonferroni")

###########################################
####NMDS for all data ####
ordr<-metaMDS(sqrt(ASV), distance = "bray", k = 2, binary=F) #stress = 0.111234
summary(ordr)
plot(ordr$points)
ordr$stress

#Plot using ggplot2
#Extract the axes scores
datascores=as.data.frame(scores(ordr)$sites) 

#Add/calculate spider diagram
scores <- cbind(as.data.frame(datascores), Treatment = Veluwemetadata$Droog.Vochtig)
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Treatment, data = scores, FUN = mean)
seg <- merge(scores, setNames(centroids, c('Treatment','oNMDS1','oNMDS2')),
             by = 'Treatment', sort = FALSE)

#plot
plotNMDS1a = ggplot(scores, aes(x = NMDS1, y = NMDS2, colour = Treatment)) +
  geom_segment(data = seg,
               mapping = aes(xend = oNMDS1, yend = oNMDS2), linewidth=.7) + # add spiders
  geom_point(data = centroids, size = 3) +                    # add centroids
  geom_point(size=1) +                                              
  coord_fixed()+                                              
  theme_bw()+ 
  theme(legend.position="top",legend.background = element_rect(color = "black", linetype = "solid"),
        legend.text=element_text(size=12),legend.title= element_blank(),legend.direction='horizontal',
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))+
   scale_color_manual(values=c("#4E84C4", "#293352"), labels=c('Dry', 'Wet'))

plotNMDS1a

#Add/calculate spider diagram
scores <- cbind(as.data.frame(datascores), Treatment = Veluwemetadata$beh)
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Treatment, data = scores, FUN = mean)
seg <- merge(scores, setNames(centroids, c('Treatment','oNMDS1','oNMDS2')),
             by = 'Treatment', sort = FALSE)

#plot
plotNMDS2a = ggplot(scores, aes(x = NMDS1, y = NMDS2, colour = Treatment)) +
  geom_segment(data = seg,
               mapping = aes(xend = oNMDS1, yend = oNMDS2), linewidth=.7) + # add spiders
  geom_point(data = centroids, size = 3) +                    # add centroids
  geom_point(size=1) +                                              
  coord_fixed()+                                              
  theme_bw()+ 
  theme(legend.position="top",legend.background = element_rect(color = "black", linetype = "solid"),
        legend.text=element_text(size=12),legend.title= element_blank(),legend.direction='horizontal',
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))+
  scale_color_manual(values =pal1, labels=c('Control', 'Biolit', 'Dolokal', 'Soilfeed'))
plotNMDS2a

ggarrange(plotNMDS1a,plotNMDS2a,align = "hv", common.legend=F)

#######################################################################
####Analyses of dry and wet site apart

#Subset ASV data with env data from dry site
ASVdroog = subset(ASV , rownames(ASV) %in% rownames(Veluwemetadatadroog))
###MATCH THE ORDER OF ENV DATA AND COMDATA = VERY IMPORTANT STEP!
ASVdroog <- ASVdroog[match(rownames(Veluwemetadatadroog), rownames(ASVdroog)),]

#Subset ASV data with env data from wet site
ASVnat = subset(ASV , rownames(ASV) %in% rownames(Veluwemetadatanat))
###MATCH THE ORDER OF ENV DATA AND COMDATA = VERY IMPORTANT STEP!
ASVnat <- ASVnat[match(rownames(Veluwemetadatanat), rownames(ASVnat)),]

####DRY####
#NMDS
ordrfundry<-metaMDS(sqrt(ASVdroog), distance = "bray", k = 2, binary=F) #0.1659256
summary(ordr)
plot(ordr$points)
ordr$stress

#Plot using ggplot2
#Extract the axes scores
datascores=as.data.frame(scores(ordr, display = "sites"))

#Add/calculate spider diagram
scores <- cbind(as.data.frame(datascores), Treatment = Veluwemetadatadroog$beh)
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Treatment, data = scores, FUN = mean)
seg <- merge(scores, setNames(centroids, c('Treatment','oNMDS1','oNMDS2')),
             by = 'Treatment', sort = FALSE)

#plot
plotNMDS1fundroog =ggplot(scores, aes(x = NMDS1, y = NMDS2, colour = Treatment)) +
  geom_segment(data = seg,
              mapping = aes(xend = oNMDS1, yend = oNMDS2), linetype="dashed") + # add spiders
 geom_point(data = scores, size = 1) +  
  geom_point(size=3, data = centroids)+ # add centroids
  #coord_fixed()+                                            
  theme_bw()+ 
  theme(legend.position="none",legend.background = element_rect(color = NA),
        legend.text=element_text(size=14),axis.title = element_text(size = 12), legend.title= element_blank(),legend.direction='horizontal',
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))+
  scale_color_manual(values =pal1,labels=c('Controle'='Control'))

plotNMDS1fundroog

finalfundry =plotNMDS1fundroog +
  inset_element(Fungdivdry, align_to="plot",left = 0.6, bottom = 0.6, right = 0.99, top = 0.99)


#Permanova
test1 <- adonis2(sqrt(ASVdroog) ~ beh, data=Veluwemetadatadroog, permutations=9999, method="bray")
test1

pairdry <- pairwise.adonis(sqrt(ASVdroog), Veluwemetadatadroog$beh, perm = 9999)


#### BC DISSIMILARIY ####
# Fungi DRY
ASVdroog$Treatment = Veluwemetadatadroog$beh
ASVdroog[, sapply(ASVdroog, is.numeric)] <- lapply(ASVdroog[, sapply(ASVdroog, is.numeric)], sqrt)
rownames(ASVdroog) <- gsub("^F", "", rownames(ASVdroog))


#calculate distance for control vs treatments.
bc_t1 <- vegdist(ASVdroog[ASVdroog$Treatment %in% c("Control", "Soilfeed"),1:884], method = "bray")
bc_t2 <- vegdist(ASVdroog[ASVdroog$Treatment %in% c("Control", "Biolit"), 1:884], method = "bray")
bc_t3 <- vegdist(ASVdroog[ASVdroog$Treatment %in% c("Control", "Dolokal"),1:884], method = "bray")
bc_t4<- vegdist(ASVdroog[ASVdroog$Treatment=="Control",1:884], method = "bray")

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

funabiodry <- ggplot(aes(y=BC, x=eucl, color = Treatment), data=distance2)+
  geom_point(aes(size=1, alpha=.5))+
  scale_color_manual(values=pal1)+theme_bw()+
  geom_smooth(data = subset(distance2, Treatment == "Soilfeed"), method = "lm", se = FALSE, aes(group = Treatment), linetype="dashed") +  # Only for a specific treatment
  geom_smooth(data = subset(distance2, Treatment == "Dolokal"), method = "lm", se = FALSE, aes(group = Treatment)) +  # Only for a specific treatment
  xlab("")+ ylab("BC dissimilarity to control (fungi)")+theme(legend.position = "")

funabiodry

### Manteltest
# control
mantel(bc_t4, eucl_t4, method = "spearman")
#Soilfeed
mantel(bc_t1, eucl_t1, method = "spearman")
#Biolit
mantel(bc_t2, eucl_t2, method = "spearman")
#Dolokal
mantel(bc_t3, eucl_t3, method = "spearman")


####WET site####
#NMDS  
ordrfunwet<-metaMDS(sqrt(ASVnat), distance = "bray", k = 2, binary=F) #0.1134634
summary(ordr)
plot(ordr$points)
ordr$stress

#Plot using ggplot2
#Extract the axes scores
datascores=as.data.frame(scores(ordr)$sites) 

#Add/calculate spider diagram
scores <- cbind(as.data.frame(datascores), Treatment = Veluwemetadatanat$beh)
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Treatment, data = scores, FUN = mean)
seg <- merge(scores, setNames(centroids, c('Treatment','oNMDS1','oNMDS2')),
             by = 'Treatment', sort = FALSE)

#plot
plotNMDS2funnat = ggplot(scores, aes(x = NMDS1, y = NMDS2, colour = Treatment)) +
  geom_segment(data = seg,
               mapping = aes(xend = oNMDS1, yend = oNMDS2),linetype="dashed") + # add spiders
  geom_point(data = centroids, size = 3) +                    # add centroids
  geom_point(size=1, data = scores) +                                              
  #coord_fixed()+                                              
  theme_bw()+
  theme(legend.position="none",legend.background = element_rect(color = NA),
        legend.text=element_text(size=14),axis.title = element_text(size = 12), legend.title= element_blank(),legend.direction='horizontal',
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))+ ylim(-0.5,0.8)+
  scale_color_manual(values =pal1)

plotNMDS2funnat

finalfunwet=  plotNMDS2funnat +
  inset_element(Fungdivwet, align_to="plot",left = 0.6, bottom = 0.6, right = 0.99, top = 0.99)

# final plot

output_dir <- "OUTPUT"
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
file_path <- file.path(output_dir, "combbacdry.png")
ggsave(file = file_path, plot = finalbacdry, width = 6, height = 6, dpi = 300)  # Adjust width, height, and dpi as needed


#Permanova
test1 <- adonis2(sqrt(ASVnat) ~ as.factor(Plot)+beh, data=Veluwemetadatanat, permutations=999, method="bray")
test1
pairwet <- pairwise.adonis(sqrt(ASVnat), Veluwemetadatanat$beh, perm = 9999)

#### BC Dissimilarity wet ####

# Fungi WET
ASVnat$Treatment = Veluwemetadatanat$beh
ASVnat[, sapply(ASVnat, is.numeric)] <- lapply(ASVnat[, sapply(ASVnat, is.numeric)], sqrt)

#calculate distance for control vs treatments.
bc_t1 <- vegdist(ASVnat[ASVnat$Treatment %in% c("Control", "Soilfeed"),1:884], method = "bray")
bc_t2 <- vegdist(ASVnat[ASVnat$Treatment %in% c("Control", "Biolit"), 1:884], method = "bray")
bc_t3 <- vegdist(ASVnat[ASVnat$Treatment %in% c("Control", "Dolokal"),1:884], method = "bray")
bc_t4<- vegdist(ASVnat[ASVnat$Treatment=="Control",1:884], method = "bray")

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

funabiowet <- ggplot(aes(y=BC, x=eucl, color = Treatment), data=distance2)+
  geom_point(aes(size=1, alpha=.5))+
  scale_color_manual(values=pal1)+theme_bw()+
  geom_smooth(data = subset(distance2, Treatment == "Dolokal"), method = "lm", se = FALSE, aes(group = Treatment), linetype="dashed") +  # Only for a specific treatment
  #geom_smooth(data = subset(distance2, Treatment == "Dolokal"), method = "lm", se = FALSE, aes(group = Treatment)) +  # Only for a specific treatment
  xlab("")+ ylab("")+theme(legend.position = "")

funabiowet

final <- ggarrange(bacabiodry, bacabiowet, funabiodry, funabiowet, vegabiodry, vegabiowet, ncol = 2, nrow = 3, labels = c("Dry", "Wet","",""), vjust = -2)+
  theme(plot.margin = margin(2,0.1,0.1,0.1, "cm")) 


file_path <- file.path(output_dir, "BCeuclFINAL.png")
ggsave(file = file_path, plot = final, width = 8, height = 10, dpi = 300)  # Adjust width, height, and dpi as needed


### Manteltest
# control
mantel(bc_t4, eucl_t4, method = "spearman")
#Soilfeed
mantel(bc_t1, eucl_t1, method = "spearman")
#Biolit
mantel(bc_t2, eucl_t2, method = "spearman")
#Dolokal
mantel(bc_t3, eucl_t3, method = "spearman")
