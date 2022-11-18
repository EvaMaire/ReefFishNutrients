#Supplemental figures
setwd("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based")
library(ggplot2)
library(tidyverse)
library(ggridges)
library(cowplot)
library(patchwork)

######################
## IMPORT DATA BASE ##
######################

load("dat.RData") 
#load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/finaldataFeb2022.RData") 

load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/nutyields_targetgroup.RData")
data <- nutyields_targetgroup %>% filter(targetgroup == 'target')
data <- data[-which(data$bm_kg_ha_sp==0),]

data$logbm <- log(data$bm_kg_ha_sp)
data$logzn <- log(data$Zinc_mg_m2)
data$logir <- log(data$Iron_mg_m2)
data$logca <- log(data$Calcium_mg_m2)
data$logom <- log(data$Omega_3_g_m2)
data$logva <- log(data$Vitamin_A_ug_m2)
data$logse <- log(data$Selenium_ug_m2)

#Add Nutrition Score
#Setting thresholds
#Vitamin A
VA_under5 <- 500*1/9 + 300*6/9 + 400*2/9 # 7–12 mo, 1–3y, 4–8y # 500*1/9 + 300*6/9 + 400*2/9

#Calcium
Ca_under5 <- 260*1/9 + 700*6/9 + 1000*2/9 # 7–12 mo, 1–3y, 4–8y

#Iron
Fe_under5 <-11*1/9 + 7*6/9 + 10*2/9 # 7–12 mo, 1–3y, 4–8y

#Total Omega 3
#Om_under5 <- 0.5*1/7 + 0.7*4/7 + 0.9*2/7 # 7–12 mo, 1–3y, 4–8y = > 0.5*1/9 + 0.7*6/9 + 0.9*2/9 # 7–12 mo, 1–3y, 4–8y
Om_under5 <- (100*1/9 + 100*2/9 + 250*4/9 + 250*2/9)/1000 # 7–11 mo, 1y, 2–3y, 4–17y

#Zinc
Zn_under5 <- 3*1/9 + 3*6/9 + 5*2/9 # 7–12 mo, 1–3y, 4–8y

#Selenium
Se_under5 <-  20*1/9 + 20*6/9 + 30*2/9 #0–6 mo, 7–12 mo, 1–3y, 4–8y

#Biomass
b500 <- log(501)
b1000 <- log(1001)

white_themejpg <-theme(axis.ticks=element_line(colour="black"),
                       axis.text=element_text(size=12,colour="black"),
                       axis.title=element_text(size=14),
                       panel.grid.minor=element_blank(),
                       panel.background=element_rect(fill="white",colour="black"),
                       plot.background=element_rect(fill="transparent",colour=NA),
                       legend.key = element_rect(fill = "white"),
                       plot.title = element_text(hjust = 0.5))

ca <- ggplot(data, aes(x=logbm, y = logca)) + 
  geom_segment(aes(x=min(data$logbm),y=min(data$logca),xend=max(data$logbm),yend=max(data$logca)))+
  geom_point(fill="dark grey",color="black",size=2,pch=21)+
  scale_x_continuous(name = "Log fish biomass")+
  scale_y_continuous(name = "Log nutrient yield")+
  ggtitle("Calcium")+
  white_themejpg 

ca

#Iron

ir <- ggplot(data, aes(x=logbm, y = logir)) + 
  geom_segment(aes(x=min(data$logbm),y=min(data$logir),xend=max(data$logbm),yend=max(data$logir)))+
  geom_point(fill="dark grey",color="black",size=2,pch=21)+
  scale_x_continuous(name = "Log fish biomass")+
  scale_y_continuous(name = " ")+
  ggtitle("Iron")+
  white_themejpg 

ir

#Zinc

zn <- ggplot(data, aes(x=logbm, y = logzn)) +  
  geom_segment(aes(x=min(data$logbm),y=min(data$logzn),xend=max(data$logbm),yend=max(data$logzn)))+
  geom_point(fill="dark grey",color="black",size=2,pch=21)+
  scale_x_continuous(name = "Log fish biomass")+
  scale_y_continuous(name = " ")+
  ggtitle("Zinc")+
  white_themejpg

zn

#Omega-3

om <- ggplot(data, aes(x=logbm, y = logom)) + 
  geom_point(fill="dark grey",color="black",size=2,pch=21)+
  scale_x_continuous(name = "Log fish biomass")+
  scale_y_continuous(name = "Log nutrient yield")+
  ggtitle("Omega-3")+
  white_themejpg

om

#Vitamin A

vit <- ggplot(data, aes(x=logbm, y = logva)) + 
  geom_point(fill="dark grey",color="black",size=2,pch=21)+
  scale_x_continuous(name = "Log fish biomass")+
  scale_y_continuous(name = " ")+
  ggtitle("Vitamin A")+
  white_themejpg

vit

#Selenium

se <- ggplot(data, aes(x=logbm, y = logse)) + 
  geom_point(fill="dark grey",color="black",size=2,pch=21)+
  scale_x_continuous(name = "Log fish biomass")+
  scale_y_continuous(name = " ")+
  ggtitle("Selenium")+
  white_themejpg

se

library(patchwork)
patch <- ca + ir + zn + plot_layout(nrow = 1) + plot_annotation(tag_levels = 'A')

jpeg("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/SI_figures/Nutrient Yields versus Fish Biomass.jpeg", res=300, width=6000, height=2000)
patch
graphics.off()

#Figure S3 - PCA on fish functional groups
#Run PCA on fish communities
library(cowplot)
library(dplyr)
library(tidyverse)
library(readr)
library(remotes)
library(RColorBrewer)
library(funrar)
library(factoextra)
library(FactoMineR)
require(gridExtra)
#remotes::install_github('jorvlan/raincloudplots')

setwd("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value")
#load("nutrient_US_DietsFeb2022.RData")
load("nutrient_US_Diets.RData")

##########################################################

comm <- nutrient_US_Diets %>% select(UniqueSite,Diets,biomass) %>%
  pivot_wider(names_from = Diets, values_from = biomass)

comm <- column_to_rownames(comm, var = "UniqueSite")
names(comm) <- c("herbivore/detritivore","invertivore (mobile)","invertivore (sessile)","piscivore","planktivore")

#Remove and clean
length(which(rowSums(comm)==0))
#comm <- comm[-which(rowSums(comm)==0),] # ONE TRANSECT WITH NO DATA 
comm[is.na(comm)] <- 0
comm <- make_relative(as.matrix(comm))
comm <- as.data.frame(comm)

#transformation => ARCSIN
commt <- asin(sqrt(comm))
comm_PCA <- prcomp(commt,scale = F) # Run PCA with no scaling

fviz_pca_var(comm_PCA, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)


#jpeg("eigenvalues.jpeg", res=300, width=2000, height=2000)
fviz_eig(comm_PCA) + labs(title ="", x = "PC", y = "% of explained variance")
#graphics.off()

d <- summary(comm_PCA)$importance #70% with the two first axes
d <- d[-1,]
eig <- tableGrob(round(d,2),theme = ttheme_minimal(base_size = 8))

comm_axes <- as.data.frame(comm_PCA$x[,1:2])
comm_axes$UniqueSite <- rownames(comm_axes)

#rename and save file
fishcomm_pca <- comm_axes
#save(fishcomm_pca,file="fishcomm_pca.RData")

#export figure

# Appendix 1 - Associations between environmental and benthic conditions of reefs through a Principal Component 
# Analysis and corresponding loadings. 

pca <- fviz_pca_var(comm_PCA, col.var = "cos2",
                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                    repel = TRUE,midpoint=20) + theme_minimal() +
  labs(title ="", x = "PC1", y = "PC2")

pca


summ <- as.data.frame(comm_PCA$rotation)
corr <- tableGrob(round(summ,2),theme = ttheme_minimal(base_size = 8))

# 2-panel Plot
grid.arrange(eig,pca,corr, widths = c(1,2),heights = c(1,1),
             layout_matrix = rbind(c(1, 2),
                                   c(3,2))) # Figure_S1 


jpeg("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/SI_figures/PCA_fishcomm.jpeg", res=300, width=3000, height=2000)
grid.arrange(eig,pca,corr, widths = c(1,2),heights = c(1,1),
             layout_matrix = rbind(c(1, 2),
                                   c(3,2))) # Figure_S1 
graphics.off()


#Figure S4
setwd("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based")
load('dat.RData')
load('/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/datasp.RData')

tim <- read.csv("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/Tim's data.csv",header=T,sep=",",dec=".")
siteTim <- unique(tim$UniqueSite)

total = dat$UniqueSite
length(which((total%in%siteTim)==F)) #1661 sites at species level
length(which((total%in%siteTim)==T)) #136 sites from Tim

goodsp <- total[which((total%in%siteTim)==F)]
length(goodsp)

#Import nutrient predictions
#Species level data
load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/nut_splevel_Sept2022.RData")
all = nut_splevel_Sept2022$UniqueSite
b = which(all%in%goodsp==T)

splevel = nut_splevel_Sept2022[b,]
splevel <- unique(splevel %>% select(FullSpeciesNut,Diets,Calcium_mu,Iron_mu,Zinc_mu,ca_rda,fe_rda,zn_rda))
splevel <- splevel[-which(is.na(splevel$ca_rda)==T),] #remove missing values = Anthias

#keep only species-level data (remove genus and family-level entries)
test <- vector()
for(k in 1:nrow(splevel)){
  test[k] <- length(unlist(str_split(splevel$FullSpeciesNut[k]," ")))
}

splevel <- splevel[-which(test==1),]
dim(splevel) 
length(unique(splevel$FullSpeciesNut)) #843 species

#Last check
sptot <- unique(splevel$FullSpeciesNut)

datfish <- unique(splevel)
summary(datfish)

#Add diets
#Complete missing values for Diets
nas <- datfish$Diets[which(datfish$Diets== -999)]
sp <- datfish$FullSpeciesNut[which(datfish$Diets== -999)] 
unique(sp)

datfish$Diets[which(datfish$FullSpeciesNut== "Decapterus macarellus")] <- "PK"
datfish$Diets[which(datfish$FullSpeciesNut== "Decapterus russelli")] <- "IM"

#nas <- datfish$Diets[which(datfish$Diets== -999)]
#sp <- datfish$FullSpeciesNut[which(datfish$Diets== -999)] 
#unique(sp)

#Family
#datfish$Diets[which(datfish$FullSpeciesNut== "Acanthuridae")] <- "HD"
#datfish$Diets[which(datfish$FullSpeciesNut== "Labridae")] <- "IM"
#datfish$Diets[which(datfish$FullSpeciesNut== "Monacanthidae")] <- "OM"
#datfish$Diets[which(datfish$FullSpeciesNut== "Serranidae")] <- "FC"
#datfish$Diets[which(datfish$FullSpeciesNut== "Holocentridae")] <- "PK"
#datfish$Diets[which(datfish$FullSpeciesNut== "Scaridae")] <- "OM"
#datfish$Diets[which(datfish$FullSpeciesNut== "Nemipteridae")] <- "IM"


#Genus
#x = datfish$Diets[which(datfish$Genus=="Acanthurus")]
#names(which.max(table(x)))
#table(x)

#datfish$Diets[which(datfish$FullSpeciesNut== "Acanthurus")] <-"HD"
#datfish$Diets[which(datfish$FullSpeciesNut== "Arothron")] <- "IS"
#datfish$Diets[which(datfish$FullSpeciesNut== "Cantherhines")] <- "OM"
#datfish$Diets[which(datfish$FullSpeciesNut== "Carangoides")] <-"FC"
#datfish$Diets[which(datfish$FullSpeciesNut== "Caranx" )] <-"FC"
#datfish$Diets[which(datfish$FullSpeciesNut== "Chaetodon")] <- "IS"
#datfish$Diets[which(datfish$FullSpeciesNut== "Epinephelus")] <- "FC"
#datfish$Diets[which(datfish$FullSpeciesNut== "Lethrinus")] <- "IM"
#datfish$Diets[which(datfish$FullSpeciesNut== "Lutjanus")] <-"IM"
#datfish$Diets[which(datfish$FullSpeciesNut== "Monotaxis")] <- "IM"
#datfish$Diets[which(datfish$FullSpeciesNut== "Mycteroperca")] <-"FC"
#datfish$Diets[which(datfish$FullSpeciesNut== "Naso")] <-"HM"
#datfish$Diets[which(datfish$FullSpeciesNut== "Scarus" )] <- "OM"

sp <- datfish$FullSpeciesNut[which(datfish$Diets== -999)] 
unique(sp)                   

datfish <- datfish %>% mutate(Diets = recode(Diets, "IM"="invertmobile",
                                             "HD" = "herbdetri",
                                             "PK" ="planktivore",
                                             "FC" = "piscivore",
                                             "IS" ="invertsessile",
                                             "OM" = "herbdetri",
                                             "HM"= "herbdetri"))
unique(datfish$Diets)

save(datfish,file="/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/datfish.RData")

head(datfish)
dim(datfish)

#VERSION 2
library(RColorBrewer)
pal <- brewer.pal(6, "PuBuGn")
#myPalette <- colorRampPalette(brewer.pal(6, "PuBuGn"))

fish.cols<-c('herbdetri'=pal[2], 
             'invertmobile'=pal[3], 
             'piscivore'=pal[4],
             'planktivore'=pal[5],
             'invertsessile'=pal[6])

fish.cols<-c('herbdetri'= "#3277B3", 
             'invertmobile'="#34B3AB", 
             'piscivore'= "#FACA21",
             'planktivore'="#E97979",
             'invertsessile'="#665C91")
#top 95% 
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

#% calcium
datcalcium <- datfish %>% select(Diets,ca_rda) 
ca <- ggplot(data=datcalcium,aes(x=Diets,y=ca_rda, fill = Diets, colour = Diets)) +
  #stat_summary(fun.data = quantiles_95, geom="boxplot",fill="dark grey")+
  #geom_hline(yintercept=mean(datcalcium$ca_rda),linetype="dashed",colour="black",size=.5)+
  geom_jitter(aes(color=Diets), size=0.5, alpha=1, width = .2) +
  geom_boxplot(outlier.shape = NA,colour="black",alpha=.6,width = .5) +
  #geom_text(data = stat,
  #aes(x=pos, y = pos90, label= paste0(round(mnuty, 0), '%')),  color="black", size=4,fontface = "bold")+
  scale_y_continuous(name="% RDA")+
  scale_x_discrete(name="",labels=c("herbivore/detritivore","invertivore (mobile)",
                                    "invertivore (sessile)","piscivore","planktivore"))+
  #geom_hline(yintercept=stat$mnuty,linetype="dashed",colour="black",size=.5)+
  coord_flip()+
  theme_cowplot()+guides(fill = FALSE, colour = FALSE) +
  scale_colour_manual(values = fish.cols)+
  scale_fill_manual(values=fish.cols) +ggtitle("calcium")

ca

datiron <- datfish %>% select(Diets,fe_rda) 
fe <- ggplot(data=datiron,aes(x=Diets,y=fe_rda, fill = Diets, colour = Diets)) +
  #geom_hline(yintercept=mean(datiron$fe_rda),linetype="dashed",colour="black",size=.5)+
  geom_jitter(aes(color=Diets), size=0.5, alpha=1, width = .2) +
  geom_boxplot(outlier.shape = NA,colour="black",alpha=.6,width = .5) +
  #geom_text(data = stat,
  #aes(x=pos, y = pos90, label= paste0(round(mnuty, 0), '%')),  color="black", size=4,fontface = "bold")+
  scale_y_continuous(name="% RDA")+
  scale_x_discrete(name="",labels=NULL)+
  #geom_hline(yintercept=stat$mnuty,linetype="dashed",colour="black",size=.5)+
  coord_flip()+
  theme_cowplot()+guides(fill = FALSE, colour = FALSE) +
  scale_colour_manual(values = fish.cols)+
  scale_fill_manual(values=fish.cols) +ggtitle("iron")

fe

datzinc <- datfish %>% select(Diets,zn_rda) 
zn <- ggplot(data=datzinc,aes(x=Diets,y=zn_rda, fill = Diets, colour = Diets)) +
  geom_jitter(aes(color=Diets), size=0.5, alpha=1, width = .2) +
  geom_boxplot(outlier.shape = NA,colour="black",alpha=.6,width = .5) +
  #geom_text(data = stat,
  #aes(x=pos, y = pos90, label= paste0(round(mnuty, 0), '%')),  color="black", size=4,fontface = "bold")+
  scale_y_continuous(name="% RDA")+
  scale_x_discrete(name="",labels=NULL)+
  #geom_hline(yintercept=stat$mnuty,linetype="dashed",colour="black",size=.5)+
  coord_flip()+
  theme_cowplot()+guides(fill = FALSE, colour = FALSE) +
  scale_colour_manual(values = fish.cols)+
  scale_fill_manual(values=fish.cols) +ggtitle("zinc")

zn

#plot1 <- plot_grid(gleg, gDiets, nrow=2, align='v', rel_heights=c(.1,1))
patchwork <- ( ca + fe + zn ) + plot_annotation(tag_levels = 'A')

jpeg("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/SI_figures/Figure_SX_percentRDA.jpeg", res=300, width=4000, height=2400)
patchwork
graphics.off()


#Version 2
#% calcium
datcalcium <- datfish %>% select(Diets,Calcium_mu) 
ca <- ggplot(data=datcalcium,aes(x=Diets,y=Calcium_mu, fill = Diets, colour = Diets)) +
  geom_hline(yintercept=mean(datcalcium$Calcium_mu),linetype="dashed",colour="black",size=1)+
  stat_summary(fun.data = quantiles_95,geom="boxplot",colour='black',alpha=.8)+
  geom_jitter(aes(color=Diets), size=0.5, alpha=1, width = .2) +
  #geom_jitter(aes(color=Diets), size=0.5, alpha=1, width = .2) +
  #geom_boxplot(outlier.shape = NA,colour="black",alpha=.6,width = .5) +
  #geom_hline(yintercept=Ca_under5/3,linetype="dashed",colour="dark grey",size=1)+
  #aes(x=pos, y = pos90, label= paste0(round(mnuty, 0), '%')),  color="black", size=4,fontface = "bold")+
  scale_y_continuous(name="Calcium concentration, mg",limits=c(0,150))+
  scale_x_discrete(name="",labels=c("herbivore/detritivore","invertivore (mobile)",
                                    "invertivore (sessile)","piscivore","planktivore"))+
  #geom_hline(yintercept=stat$mnuty,linetype="dashed",colour="black",size=.5)+
  coord_flip()+
  theme_cowplot()+guides(fill = FALSE, colour = FALSE) +
  scale_colour_manual(values = fish.cols)+
  scale_fill_manual(values=fish.cols) #+ggtitle("calcium")

ca

datiron <- datfish %>% select(Diets,Iron_mu) 
fe <- ggplot(data=datiron,aes(x=Diets,y=Iron_mu, fill = Diets, colour = Diets)) +
  geom_hline(yintercept=mean(datiron$Iron_mu),linetype="dashed",colour="black",size=1)+
  stat_summary(fun.data = quantiles_95,geom="boxplot",colour='black',alpha=.8)+
  geom_jitter(aes(color=Diets), size=0.5, alpha=1, width = .2) +
  #geom_boxplot(outlier.shape = NA,colour="black",alpha=.6,width = .5) +
  #geom_hline(yintercept=Fe_under5/3,linetype="dashed",colour="dark grey",size=1)+
  #aes(x=pos, y = pos90, label= paste0(round(mnuty, 0), '%')),  color="black", size=4,fontface = "bold")+
  scale_y_continuous(name="Iron concentration, mg",limits=c(0,1.2))+
  scale_x_discrete(name="",labels=NULL)+
  #geom_hline(yintercept=stat$mnuty,linetype="dashed",colour="black",size=.5)+
  coord_flip()+
  theme_cowplot()+guides(fill = FALSE, colour = FALSE) +
  scale_colour_manual(values = fish.cols)+
  scale_fill_manual(values=fish.cols) #+ggtitle("iron")

fe

datzinc <- datfish %>% select(Diets,Zinc_mu) 
zn <- ggplot(data=datzinc,aes(x=Diets,y=Zinc_mu, fill = Diets, colour = Diets)) +
  geom_hline(yintercept=mean(datzinc$Zinc_mu,na.rm=T),linetype="dashed",colour="black",size=1)+
  stat_summary(fun.data = quantiles_95,geom="boxplot",colour='black',alpha=.8)+
  geom_jitter(aes(color=Diets), size=0.5, alpha=1, width = .2) +
  #geom_boxplot(outlier.shape = NA,colour="black",alpha=.6,width = .5) +
  #geom_hline(yintercept=Zn_under5/3,linetype="dashed",colour="dark grey",size=1)+
  #aes(x=pos, y = pos90, label= paste0(round(mnuty, 0), '%')),  color="black", size=4,fontface = "bold")+
  scale_y_continuous(name="Zinc concentration, mg",limits=c(0,3))+
  scale_x_discrete(name="",labels=NULL)+
  #geom_hline(yintercept=stat$mnuty,linetype="dashed",colour="black",size=.5)+
  coord_flip()+
  theme_cowplot()+guides(fill = FALSE, colour = FALSE) +
  scale_colour_manual(values = fish.cols)+
  scale_fill_manual(values=fish.cols) #+ggtitle("zinc")

zn

#plot1 <- plot_grid(gleg, gDiets, nrow=2, align='v', rel_heights=c(.1,1))
patchwork <- ( ca + fe + zn ) + plot_annotation(tag_levels = 'A')

jpeg("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/SI_figures/Figure_SX_NutConc.jpeg", res=300, width=4500, height=1500)
patchwork
graphics.off()

#Figure S5
##########
library(ggrepel)
load('/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/nutrient_US_Diets.RData')
head(nutrient_US_Diets) 

nut_fg <- nutrient_US_Diets %>% select(UniqueSite,Diets,calcium,iron,zinc) %>% group_by(UniqueSite) %>%
  mutate(across(calcium:zinc, list(tot = ~sum(.x)))) %>% 
  ungroup() %>%
  mutate(propzinc = 100*zinc/zinc_tot,
         propiron = 100*iron/iron_tot,
         propcalcium = 100*calcium/calcium_tot)%>%#,
         #propomega_3 = 100*omega_3/omega_3_tot,
         #propvitaminA = 100*vitaminA/vitaminA_tot,
         #propselenium = 100*selenium/selenium_tot) 
  select(UniqueSite,Diets,propzinc:propcalcium) %>%
  ## turn to long format
  pivot_longer(propzinc:propcalcium, names_to = 'nutrient', values_to = 'nutprop') %>% group_by(Diets,nutrient) %>% summarise(se = funk::se(nutprop), nutprop = mean(nutprop)) %>% 
  mutate(lower = nutprop - 2*se, upper = nutprop + 2*se)

## average across nutrients for panel A
nut_fg2 <- nut_fg %>% group_by(Diets) %>% summarise(se = funk::se(nutprop), nutprop = mean(nutprop)) %>% 
  mutate(lower = nutprop - 2*se, upper = nutprop + 2*se)

d <- nut_fg %>% mutate(nutrientshort = recode(nutrient,
                                              "propcalcium"="calcium",
                                              "propiron"="iron",
                                              #"propomega_3"="omega-3",
                                              #"propselenium"="selenium",
                                              #"propvitaminA"="vit. A",
                                              "propzinc"="zinc" )) #creat new column with clean names

#devtools::install_github("slowkow/ggrepel")

pos <- position_jitter(width = 0.25, seed = 1)

g1 <- ggplot(d, aes(fct_reorder(Diets, nutprop),nutprop, fill=Diets, label=nutrientshort),col='black') + 
  #geom_point(size=3, pch=21, position = pos) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size=.5, pch=21,position = pos) +
  #geom_text(position = position_jitter(seed = 1))+
  geom_text_repel(position = pos,     max.iter = 3e3,segment.size = 0.5,
                  #point.padding = unit(0.2, "lines"),
                  max.overlaps = 20,
                  direction = "both",
                  # Add extra padding around each text label.
                  #box.padding = unit(.2, 'cm'),
                  # Color of the line segments.
                  segment.color = 'black',
                  # Width and transparency of the line segments.
                  segment.alpha= .4,
                  min.segment.length = 0.2)+
  coord_flip() +
  theme_cowplot() + theme(legend.position = 'none') +
  scale_colour_manual(values = fish.cols)+
  scale_fill_manual(values=fish.cols) +
  scale_y_continuous(limits=c(0,60),breaks=seq(0,60, by = 10), labels=seq(0, 60, by = 10)) +
  scale_x_discrete(name="",labels=rev(c("herbivore/detritivore","invertivore (mobile)",
                                        "piscivore","planktivore","invertivore (sessile)")))+
  labs(x = '', y = "Mean proportion of nutrient yield, %") 

g1

jpeg("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/SI_figures/Contri_nut_yields_FigureSX.jpeg", res=300, width=2600, height=2400)
g1
graphics.off()

########################################################################
#Figure S6 - Nutrient concentrations of fish families
setwd("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based")
load('dat.RData')

tim <- read.csv("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/Tim's data.csv",header=T,sep=",",dec=".")
siteTim <- unique(tim$UniqueSite)

total = dat$UniqueSite
length(which((total%in%siteTim)==F)) #1661 sites at species level
length(which((total%in%siteTim)==T)) #136 sites from Tim

goodsp <- total[which((total%in%siteTim)==F)]
length(goodsp)

#Import nutrient predictions
#Species level data
load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/nut_splevel_Sept2022.RData")
all = nut_splevel_Sept2022$UniqueSite
b = which(all%in%goodsp==T)

splevel = nut_splevel_Sept2022[b,]
splevel <- unique(splevel %>% select(FullSpeciesNut,Family,Calcium_mu,Iron_mu,Zinc_mu,ca_rda,fe_rda,zn_rda))
splevel <- splevel[-which(is.na(splevel$ca_rda)==T),] #remove missing values = Anthias

#keep only species-level data (remove genus and family-level entries)
test <- vector()
for(k in 1:nrow(splevel)){
  test[k] <- length(unlist(str_split(splevel$FullSpeciesNut[k]," ")))
}

splevel <- splevel[-which(test==1),]
dim(splevel) 
length(unique(splevel$FullSpeciesNut)) #843 species

#Last check
sptot <- unique(splevel$FullSpeciesNut)

datfish <- unique(splevel)
summary(datfish)

quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

#Calcium
ca <- ggplot(datfish,aes(x=Calcium_mu, y=Family)) +
  geom_vline(xintercept=mean(datfish$Calcium_mu),linetype="dashed",colour="black",size=.5)+
  stat_summary(fun.data = quantiles_95,geom="boxplot",colour='black',fill="light grey")+
  geom_jitter( size=0.5, alpha=1, width = .2) +
  #geom_boxplot(outlier.shape = NA,colour="black",alpha=.6,width = .5) +
  #geom_hline(yintercept=Zn_under5/3,linetype="dashed",colour="dark grey",size=1)+
  #aes(x=pos, y = pos90, label= paste0(round(mnuty, 0), '%')),  color="black", size=4,fontface = "bold")+
  scale_x_continuous(name="Calcium concentration, mg",limits=c(0,300))+
  scale_y_discrete(name="")+
  #geom_hline(yintercept=stat$mnuty,linetype="dashed",colour="black",size=.5)+
  #coord_flip()+
  theme_cowplot()+guides(fill = FALSE, colour = FALSE) +theme(legend.position = "none") 

ca 

zn <- ggplot(datfish,aes(x=Zinc_mu, y=Family)) +
  geom_vline(xintercept=mean(datfish$Zinc_mu),linetype="dashed",colour="black",size=.5)+
  stat_summary(fun.data = quantiles_95,geom="boxplot",colour='black',fill="light grey")+
  geom_jitter( size=0.5, alpha=1, width = .2) +
  #geom_boxplot(outlier.shape = NA,colour="black",alpha=.6,width = .5) +
  #geom_hline(yintercept=Zn_under5/3,linetype="dashed",colour="dark grey",size=1)+
  #aes(x=pos, y = pos90, label= paste0(round(mnuty, 0), '%')),  color="black", size=4,fontface = "bold")+
  scale_x_continuous(name="Zinc concentration, mg",limits=c(0,3.5))+
  scale_y_discrete(name="") +
  #geom_hline(yintercept=stat$mnuty,linetype="dashed",colour="black",size=.5)+
  #coord_flip()+
  theme_cowplot()+guides(fill = FALSE, colour = FALSE) +theme(legend.position = "none") 

zn

ir <- ggplot(datfish,aes(x=Iron_mu, y=Family)) +
  geom_vline(xintercept=mean(datfish$Iron_mu),linetype="dashed",colour="black",size=.5)+
  stat_summary(fun.data = quantiles_95,geom="boxplot",colour='black',fill="light grey")+
  geom_jitter( size=0.5, alpha=1, width = .2) +
  #geom_boxplot(outlier.shape = NA,colour="black",alpha=.6,width = .5) +
  #aes(x=pos, y = pos90, label= paste0(round(mnuty, 0), '%')),  color="black", size=4,fontface = "bold")+
  scale_x_continuous(name="Iron concentration, mg",limits=c(0,3))+
  scale_y_discrete(name="",labels=NULL) +
  #geom_hline(yintercept=stat$mnuty,linetype="dashed",colour="black",size=.5)+
  #coord_flip()+
  theme_cowplot()+guides(fill = FALSE, colour = FALSE) +theme(legend.position = "none") 

ir

patchwork <- ( ca + ir + zn ) + plot_annotation(tag_levels = 'A') +  plot_layout(ncol = 2)

jpeg("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/SI_figures/Figure_SX_NutConc_families.jpeg", res=300, width=4500, height=4500)
patchwork
graphics.off()


##############
#FIGURE S7 - Relative biomass of functional groups along a biomass gradient

#Import biomass of fish groups
setwd("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value")
load("nutrient_US_Diets.RData")
head(nutrient_US_Diets) 

#remove 1 outlier 
load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/dat.RData") 
finaldata <- dat %>% filter(targetgroup == 'target')
finaldata <- finaldata[-which(finaldata$biomasstarget==0),] #remove the 1 site with non target biomass
goodUS <- finaldata$UniqueSite

nut <- nutrient_US_Diets[-which((nutrient_US_Diets$UniqueSite%in%goodUS)==F),]

#Settings for plot

pal <- brewer.pal(6, "PuBuGn")
#myPalette <- colorRampPalette(brewer.pal(6, "PuBuGn"))

fish.cols<-c('herbdetri'=pal[2], 
             'invertmobile'=pal[3], 
             'piscivore'=pal[4],
             'planktivore'=pal[5],
             'invertsessile'=pal[6])

fish.cols<-c('herbdetri'= "#3277B3", 
             'invertmobile'="#34B3AB", 
             'piscivore'= "#FACA21",
             'planktivore'="#E97979",
             'invertsessile'="#665C91")
#Plot
white_themejpg <-theme(axis.ticks=element_line(colour="black"),
                       axis.text=element_text(size=12,colour="black"),
                       axis.title=element_text(size=14),
                       panel.grid.minor=element_blank(),
                       panel.background=element_rect(fill="white",colour="black"),
                       plot.background=element_rect(fill="transparent",colour=NA),
                       legend.key = element_rect(fill = "white"))

#Relationship with gravity
#% biomass
datbm <- nut %>% select(UniqueSite,Diets,biomass) %>%
  group_by(UniqueSite) %>%
  mutate(bm_tot = sum(biomass)) %>% #compute total bm for each site
  ungroup() %>%
  mutate(propbm = biomass/bm_tot) 

tomergeG <- finaldata[,c("UniqueSite","loggrav","logbm","Protection")]
datG <- merge(datbm,tomergeG,by="UniqueSite",all.x=T)

grav <- ggplot(datG, aes(x=loggrav, y=propbm, colour = Diets)) + 
  #geom_vline(xintercept=log(500),linetype="dashed",colour="black",size=0.5)+
  geom_point(aes(colour = Diets),size=1,alpha=0.4) +
  geom_smooth(method="glm", method.args = list(family = "quasipoisson"), formula = y ~ poly(x, 3), se = TRUE, fullrange = FALSE,
              size = 1, 
              aes(colour = Diets,fill=Diets))+
  scale_colour_manual(values = fish.cols,name="Fish functional groups",labels=c("herbivore/detritivore","invertivore (mobile)",
                                                                                "piscivore","planktivore","invertivore (sessile)"))+
  scale_fill_manual(values = fish.cols,name="Fish functional groups",labels=c("herbivore/detritivore","invertivore (mobile)",
                                                                                "piscivore","planktivore","invertivore (sessile)"))+
  scale_x_continuous(name = "LOG (Human gravity)",limits=c(0,10),breaks = c(0,2,5,7.5,10)) + 
  scale_y_continuous(name = "Relative biomass (%)", limits = c(0,1), 
                     breaks = c(0,.25,.50,.75,1),labels = c("0","25","50","75","100")) + 
  white_themejpg + theme(legend.position="none")#+ facet_wrap(~ Protection, ncol = 3) #+ theme(legend.position="none") +theme(plot.title = element_text(hjust = 0.5))


grav

biom <- ggplot(datG, aes(x=logbm, y=propbm, colour = Diets)) + 
  geom_vline(xintercept=log(500),linetype="dashed",colour="black",size=0.5)+
  geom_point(aes(colour = Diets),size=1,alpha=0.4) +
  geom_smooth(method="glm", method.args = list(family = "quasipoisson"), formula = y ~ poly(x, 3), se = TRUE, fullrange = FALSE,
              size = 1, 
              aes(colour = Diets,fill=Diets))+
  scale_colour_manual(values = fish.cols,name="Fish functional groups",labels=c("herbivore/detritivore","invertivore (mobile)",
                                                                                "piscivore","planktivore","invertivore (sessile)"))+
  scale_fill_manual(values = fish.cols,name="Fish functional groups",labels=c("herbivore/detritivore","invertivore (mobile)",
                                                                              "piscivore","planktivore","invertivore (sessile)"))+
  scale_x_continuous(name = "LOG (Standing fish biomass)",limits=c(.7,10.1),breaks = c(2,5,7.5,10)) + 
  scale_y_continuous(name = "", limits = c(0,1), 
                     breaks = c(0,.25,.50,.75,1),labels = c("0","25","50","75","100")) + 
  white_themejpg# + theme(legend.position="none")#+ facet_wrap(~ Protection, ncol = 3) #+ theme(legend.position="none") +theme(plot.title = element_text(hjust = 0.5))

biom

patchwork <- (grav + biom)  + plot_annotation(tag_levels = 'A')

#Final Plot
jpeg("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/SI_figures/Figure_SX_Relative_biomass_Grav_Biomass.jpeg", res=300, width=4200, height=2000)
patchwork
graphics.off()

#relationships with management
datbm$logbm <- log(datbm$bm_tot)

#Add protection
tomerge <- finaldata[,c("UniqueSite","Protection")]
dat <- merge(datbm,tomerge,by="UniqueSite",all.x=T)

#No-take
NT <- dat[which(dat$Protection=="UnfishedHigh"),]

bmNT <- ggplot(NT, aes(x=logbm, y=propbm, colour = Diets)) + 
  geom_vline(xintercept=log(500),linetype="dashed",colour="black",size=0.5)+
  geom_point(aes(colour = Diets),size=1,alpha=0.4) +
  geom_smooth(method="glm", method.args = list(family = "quasipoisson"), formula = y ~ poly(x, 3), se = TRUE, fullrange = FALSE,
              size = 1, 
              aes(colour = Diets,fill=Diets))+
  scale_fill_manual(values = fish.cols)+
  scale_colour_manual(values = fish.cols)+
  scale_x_continuous(name = "LOG (total fish biomass)",limits=c(.7,10.1),breaks = c(2,5,7.5,10)) + 
  scale_y_continuous(name = "Relative biomass (%)", limits = c(0,1.1), 
                     breaks = c(0,.25,.50,.75,1),labels = c("0","25","50","75","100")) + 
  white_themejpg+ theme(legend.position="none") + ggtitle("Marine reserves")+theme(plot.title = element_text(hjust = 0.5))


bmNT

#restricted
rest <- dat[which(dat$Protection=="Restricted"),]
#rest <- rest[-which(rest$logbm<1),]
bmrest <- ggplot(rest, aes(x=logbm, y=propbm, colour = Diets)) + 
  geom_vline(xintercept=log(500),linetype="dashed",colour="black",size=0.5)+
  geom_point(aes(colour = Diets),size=1,alpha=0.4) +
  geom_smooth(method="glm", method.args = list(family = "quasipoisson"), formula = y ~ poly(x, 3), se = TRUE, fullrange = FALSE,
              size = 1, 
              aes(colour = Diets,fill=Diets))+
  scale_fill_manual(values = fish.cols)+
  scale_colour_manual(values = fish.cols)+
  scale_x_continuous(name = "LOG (total fish biomass)",limits=c(.7,10.1),breaks = c(2,5,7.5,10)) + 
  scale_y_continuous(name = "Relative biomass (%)", limits = c(0,1.1), 
                     breaks = c(0,.25,.50,.75,1),labels = c("0","25","50","75","100")) + 
  white_themejpg+ theme(legend.position="none")+ggtitle("Restricted areas")+theme(plot.title = element_text(hjust = 0.5))

bmrest

#Fished
fished<- dat[which(dat$Protection=="Fished"),]
bmfished <- ggplot(fished, aes(x=logbm, y=propbm, colour = Diets)) + 
  geom_vline(xintercept=log(500),linetype="dashed",colour="black",size=0.5)+
  geom_point(aes(colour = Diets),size=1,alpha=0.4) +
  geom_smooth(method="glm", method.args = list(family = "quasipoisson"), formula = y ~ poly(x, 3), se = TRUE, fullrange = FALSE,
              size = 1, 
              aes(colour = Diets,fill=Diets))+
  scale_fill_manual(values = fish.cols)+
  scale_colour_manual(values = fish.cols)+
  scale_x_continuous(name = "LOG (total fish biomass)",limits=c(.7,10.1),breaks = c(2,5,7.5,10)) + 
  scale_y_continuous(name = "Relative biomass (%)", limits = c(0,1.1), 
                     breaks = c(0,.25,.50,.75,1),labels = c("0","25","50","75","100")) + 
  white_themejpg+ theme(legend.position="none")+ggtitle("Fished areas")+theme(plot.title = element_text(hjust = 0.5))

bmfished


pleg <- ggplot(NT, aes(x=logbm, y=propbm, colour = Diets)) + 
  geom_point(aes(colour = Diets),size=3,alpha=0.4) +
  geom_smooth(method="glm", method.args = list(family = "quasipoisson"), formula = y ~ poly(x, 3), se = TRUE, fullrange = TRUE,
              size = 1, 
              aes(colour = Diets))+
  scale_colour_manual(values = fish.cols,name="Fish functional groups",labels=c("herbivore/detritivore","invertivore (mobile)",
                                                                                "piscivore","planktivore","invertivore (sessile)"))+
  scale_x_continuous(name = "",limits=c(0,10.1)) + 
  scale_y_continuous(name = "Relative biomass (%)", limits = c(0,1), 
                     breaks = c(0,.25,.50,.75,1),labels = c("0","25","50","75","100")) + 
  white_themejpg

library(cowplot)
legend<- get_legend(
  pleg + 
    theme(legend.position = "bottom",legend.justification="center",legend.key.width = unit(3,"line"),
          legend.text=element_text(size=14),legend.title = element_text(size=16 ,face="bold"))
)


#Final Plot
jpeg("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/SI_figures/Figure_SX_Relative_biomass_3ManagementTypes.jpeg", res=300, width=5000, height=2400)
cowplot::plot_grid(bmfished,bmrest,bmNT,
                   NULL,legend,NULL,
                   nrow = 2, ncol = 3, rel_heights = c(1,.1), rel_widths = c(1,1,1),labels = c('A', 'B','C'),label_size=21,align = "h")
graphics.off()
###################################
#Figure S8
setwd("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based")

#Species level data
load('/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/speciesdata_serf.RData')
load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/nut_splevel_Sept2022.RData")

nutscore <- unique( nut_splevel_Sept2022 %>% select(FullSpecies,nutscore_3nutA,Diets))

spdat <- merge(speciesdata_serf,nutscore,by="FullSpecies",all.x=T)
dim(speciesdata_serf)
dim(spdat) 
length(unique(spdat$FullSpeciesNut)) #838 yeah!!!!

iucn <- read.csv("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/iucn_status_sp.csv",sep=";")
stat <- merge(spdat,iucn,by.x='FullSpeciesNut',by.y='species',all.x=T)

stat$IUCN_final[which(stat$IUCN_final=='Non Threatened')]<- "Not Threatened"
stat$IUCN_final[which(is.na(stat$IUCN_final)==T)]<- "No Status"
stat$IUCN_final <- as.factor(stat$IUCN_final)
stat$newnut <- stat$nutscore_3nutA/3

quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

#remove 2 duplicates species
stat2 <- unique(stat %>% select(FullSpeciesNut,newnut,IUCN_final))

#other panel
s8B <- ggplot(stat2,aes(x=newnut, y=IUCN_final)) +
  geom_vline(xintercept=mean(stat2$newnut),linetype="dashed",colour="black",size=.5)+
  stat_summary(fun.data = quantiles_95,geom="boxplot",colour='black',fill="light grey",alpha=.6)+
  annotate("text", y = 3, x = 27, size=7, label = "*", fontface = c("bold")) +
  #geom_jitter( size=0.5, alpha=1, width = .2) +
  #geom_boxplot(outlier.shape = NA,colour="black",alpha=.6,width = .5) +
  #geom_hline(yintercept=Zn_under5/3,linetype="dashed",colour="dark grey",size=1)+
  #aes(x=pos, y = pos90, label= paste0(round(mnuty, 0), '%')),  color="black", size=4,fontface = "bold")+
  scale_x_continuous(name="Micronutrient density score, %")+
  scale_y_discrete(name="")+
  #geom_hline(yintercept=stat$mnuty,linetype="dashed",colour="black",size=.5)+
  coord_flip()+
  theme_cowplot()+guides(fill = FALSE, colour = FALSE) 

s8B

#explore species
th <- unique(stat %>% select(FullSpeciesNut,Family, Diets:newnut) %>% filter(IUCN_final == 'Threatened'))

#same with functional distinctiveness
load('/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/functiodata_sp.RData')
head(functiodata_sp)

final <- functiodata_sp
final$newnut <- final$nutscore_3nutA/3

s8a <- ggplot(final, aes(meanFdist,newnut)) +
  geom_point(color="dark grey",size=2,alpha=.6)+
  stat_smooth(color="#2f3030",method = "gam", formula = y ~ s(x, bs = "cs", k = 3)
              ,se=T,fullrange = F,size=1,show.legend = F) +
  scale_x_continuous(name = "Functional distinctiveness",limits=c(0.3,0.65))+
  scale_y_continuous(name = "Micronutrient density score, %",limits=c(0,60))+
  #scale_fill_manual(name = "Protection",values=cols)+
  theme_cowplot()+guides(fill = FALSE, colour = FALSE) 

s8a

#by functional groups
load('/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/dat.RData')
data <- dat 

tim <- read.csv("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/Tim's data.csv",header=T,sep=",",dec=".")
siteTim <- unique(tim$UniqueSite)

total = data$UniqueSite
length(which((total%in%siteTim)==F)) #1662 sites at species level
length(which((total%in%siteTim)==T)) #136 sites from Tim

goodsp <- total[which((total%in%siteTim)==F)]
length(goodsp)

#Import nutrient predictions
#Species level data
all = nut_splevel_Sept2022$UniqueSite
b = which(all%in%goodsp==T)

splevel = nut_splevel_Sept2022[b,]
splevel <- unique(splevel %>% select(FullSpeciesNut,Family,Diets,ca_rda))
splevel <- splevel[-which(is.na(splevel$ca_rda)==T),] #remove missing values = Anthias

#keep only species-level data (remove genus and family-level entries)
test <- vector()
for(k in 1:nrow(splevel)){
  test[k] <- length(unlist(str_split(splevel$FullSpeciesNut[k]," ")))
}

splevel <- splevel[-which(test==1),]
#Complete missing data (-999)
splevel[which(splevel$Diets== -999),] #all Decapterus sp. ==> IM
splevel$Diets[which(splevel$Diets== -999)] <- "IM"

t <- merge(final,splevel)

t <- t %>% mutate(Diets = recode(Diets, "IM"="invertmobile",
                                             "HD" = "herbdetri",
                                             "PK" ="planktivore",
                                             "FC" = "piscivore",
                                             "IS" ="invertsessile",
                                             "OM" = "herbdetri",
                                             "HM"= "herbdetri"))

fish.cols<-c('herbdetri'= "#3277B3", 
             'invertmobile'="#34B3AB", 
             'piscivore'= "#FACA21",
             'planktivore'="#E97979",
             'invertsessile'="#665C91")

tt <- ggplot(t, aes(x=meanFdist,y=newnut,colour=Diets)) +
  geom_point(size=2,alpha=.6)+
  stat_smooth(color="#2f3030",method = "gam", formula = y ~ s(x, bs = "cs", k = 3)
              ,se=T,fullrange = F,size=1,show.legend = F) +
  scale_x_continuous(name = "Functional distinctiveness",limits=c(0.3,0.65))+
  scale_y_continuous(name = "Micronutrient density score, %",limits=c(0,60))+
  scale_colour_manual(values = fish.cols,name="Fish functional groups",labels=c("herbivore/detritivore","invertivore (mobile)",
                                                                                "piscivore","planktivore","invertivore (sessile)"))+
  theme_cowplot()

#version 2
ff <- unique(final %>% select(FullSpeciesNut,newnut,meanFdist))
ff <- ff %>% group_by(FullSpeciesNut) %>%
  summarise(nut = mean(newnut),
            dist = mean(meanFdist))

s8a2 <- ggplot(ff %>% filter(nut<42), aes(dist,nut)) +
  geom_point(color="dark grey",size=2,alpha=.6)+
  stat_smooth(color="#2f3030",method = "lm",se=T,fullrange = F,size=1,show.legend = F) +
  scale_x_continuous(name = "Functional distinctiveness",limits=c(0.3,0.65))+
  scale_y_continuous(name = "Micronutrient density score, %",limits=c(0,60))+
  #scale_fill_manual(name = "Protection",values=cols)+
  theme_cowplot()+guides(fill = FALSE, colour = FALSE) 

s8a2


res.aov <- aov(newnut ~ IUCN_final, data = stat)
par(mar=c(4.1,10,4,4))
plot(TukeyHSD(res.aov, conf.level = 0.95),las=1, col = "red")

patchwork <- (tt + s8B)  + plot_annotation(tag_levels = 'A')

jpeg("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/SI_figures/Figure_SX_Nutscore_IUCN.jpeg", res=300, width=4200, height=2000)
patchwork             
graphics.off()

#OLD SCRIPTS
#CORRELOGRAM
#Correlogram
require(corrgram)
require(corrplot)

d <- data[,c("calciumse100", "zincse100", "ironse100", "omegase100", "vitse100", "seleniumse100")]
#d <- data[,c("Calcium_mg_m2_US","Zinc_mg_m2_US","Iron_mg_m2_US","Omega_3_g_m2_US","Vitamin_A_ug_m2_US","Selenium_ug_m2_US")]
summary(d)
colnames(d)<-c("Calcium","Zinc","Iron","Omega-3","Vitamin A","Selenium")
ccor <- d
cor_mat <- cor(ccor)

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


# matrix of the p-value of the correlation
p.mat <- cor.mtest(ccor)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

raw <- corrplot(cor_mat, method="color", col=col(200),  
                type="upper",  
                addCoef.col = "black", # Add coefficient of correlation
                tl.col="black", tl.srt=45, tl.cex = 0.9, number.cex = .7, #Text label color and rotation
                # Combine with significance
                p.mat = p.mat, sig.level = 1, insig = "blank",  
                # hide correlation coefficient on the principal diagonal
                diag=FALSE 
)

jpeg("correlogram nutrient concentrations.jpeg", res=300, width=2400, height=2400)
corrplot(cor_mat, method="color", col=col(200),  
         type="upper",  
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex = 0.9, number.cex = .7, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 1, insig = "blank",  
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)
graphics.off() 

#Histograms nutrient concentrations
#GGPLOT
white_theme <-theme(axis.text.y=element_text(size=11.5, colour='black'),
                    axis.text.x=element_text(size=11.5, colour='black'),
                    axis.line = element_line(colour = "black"), 
                    axis.line.x.top = element_line(colour = "white"), 
                    axis.title = element_text(size=12, colour='black'),
                    panel.border = element_blank(),
                    strip.text = element_text(face="bold", colour='black', size=10),
                    legend.position ='none',
                    legend.title=element_blank(),
                    axis.ticks=element_line(colour="black"),
                    plot.margin = margin(0,1,0,1, "cm"))

cah <- ggplot(data, aes(calciumse100)) + 
  geom_histogram(alpha=0.7, binwidth=max(data$calciumse100)/30, fill = "grey", col='black') +
  geom_vline(xintercept=Ca_under5/3,linetype="dashed",colour="black",size=1)+
  scale_y_continuous(expand= c(0,0), name = "# reef sites") +
  scale_x_continuous(name= "Calcium, mg/100g")+
  white_theme

cah

znh <- ggplot(data, aes(zincse100)) + 
  geom_histogram(alpha=0.7, binwidth=max(data$zincse100)/30, fill = "grey", col='black') +
  geom_vline(xintercept=Zn_under5/3,linetype="dashed",colour="black",size=1)+
  scale_y_continuous(expand= c(0,0), name = "# reef sites") +
  scale_x_continuous(name= "Zinc, mg/100g")+
  white_theme

znh

irh <- ggplot(data, aes(ironse100)) + 
  geom_histogram(alpha=0.7, binwidth=max(data$ironse100)/30, fill = "grey", col='black') +
  geom_vline(xintercept=Fe_under5/3,linetype="dashed",colour="black",size=1)+
  scale_y_continuous(expand= c(0,0), name = "# reef sites") +
  scale_x_continuous(name= "Iron, mg/100g")+
  white_theme

irh

omh <- ggplot(data, aes(omegase100)) + 
  geom_histogram(alpha=0.7, binwidth=max(data$omegase100)/30, fill = "grey", col='black') +
  geom_vline(xintercept=Om_under5/3,linetype="dashed",colour="black",size=1)+
  scale_y_continuous(expand= c(0,0), name = "# reef sites") +
  scale_x_continuous(name= "Omega-3, g/100g")+
  white_theme

omh

seh <- ggplot(data, aes(seleniumse100)) + 
  geom_histogram(alpha=0.7, binwidth=max(data$seleniumse100)/30, fill = "grey", col='black') +
  geom_vline(xintercept=Se_under5/3,linetype="dashed",colour="black",size=1)+
  scale_y_continuous(expand= c(0,0), name = "# reef sites") +
  scale_x_continuous(name= expression(paste('Selenium, ', mu,'g/100g')))+
  white_theme

seh

vah <- ggplot(data, aes(vitse100)) + 
  geom_histogram(alpha=0.7, binwidth=max(data$vitse100)/30, fill = "grey", col='black') +
  geom_vline(xintercept=VA_under5/3,linetype="dashed",colour="black",size=1)+
  scale_y_continuous(expand= c(0,0), name = "# reef sites") +
  scale_x_continuous(name= expression(paste('Vitamin A, ', mu,'g/100g')))+
  white_theme

vah

patch <- cah + omh + irh +seh + znh + vah + plot_layout(nrow = 2, byrow = FALSE)

jpeg("ALL Hist.jpeg", res=300, width=6000, height=4000)
patch
graphics.off()

#Remove outlier
#data <- data[-which(data$UniqueSite=="7988"),]

##########################################################################
#figure S9
setwd("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value")
#load("nutrient_US_DietsFeb2022.RData")
load("nutrient_US_Diets.RData")
head(nutrient_US_Diets) 
nut_fg <- nutrient_US_Diets %>% select(UniqueSite,Diets,zinc:selenium) %>% group_by(UniqueSite) %>%
  mutate(across(zinc:selenium, list(tot = ~sum(.x)))) %>% 
  ungroup() %>%
  mutate(propzinc = 100*zinc/zinc_tot,
         propiron = 100*iron/iron_tot,
         propcalcium = 100*calcium/calcium_tot,
         propomega_3 = 100*omega_3/omega_3_tot,
         propvitaminA = 100*vitaminA/vitaminA_tot,
         propselenium = 100*selenium/selenium_tot) %>%
  select(UniqueSite,Diets,propzinc:propselenium) %>%
  ## turn to long format
  pivot_longer(propzinc:propselenium, names_to = 'nutrient', values_to = 'nutprop') %>% group_by(Diets,nutrient) %>% summarise(se = funk::se(nutprop), nutprop = mean(nutprop)) %>% 
  mutate(lower = nutprop - 2*se, upper = nutprop + 2*se)

## average across nutrients for panel A
nut_fg2 <- nut_fg %>% group_by(Diets) %>% summarise(se = funk::se(nutprop), nutprop = mean(nutprop)) %>% 
  mutate(lower = nutprop - 2*se, upper = nutprop + 2*se)

d <- nut_fg %>% mutate(nutrientshort = recode(nutrient,
                                              "propcalcium"="calcium",
                                              "propiron"="iron",
                                              "propomega_3"="omega-3",
                                              "propselenium"="selenium",
                                              "propvitaminA"="vit. A",
                                              "propzinc"="zinc" )) #creat new column with clean names

#devtools::install_github("slowkow/ggrepel")
library(ggrepel)

pos <- position_jitter(width = 0.25, seed = 1)

g1 <- ggplot(d, aes(fct_reorder(Diets, nutprop),nutprop, fill=Diets, label=nutrientshort),col='black') + 
  #geom_point(size=3, pch=21, position = pos) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size=.5, pch=21,position = pos) +
  #geom_text(position = position_jitter(seed = 1))+
  geom_text_repel(position = pos,     max.iter = 3e3,segment.size = 0.5,
                  #point.padding = unit(0.2, "lines"),
                  max.overlaps = 20,
                  direction = "both",
                  # Add extra padding around each text label.
                  #box.padding = unit(.2, 'cm'),
                  # Color of the line segments.
                  segment.color = 'black',
                  # Width and transparency of the line segments.
                  segment.alpha= .4,
                  min.segment.length = 0.2)+
  coord_flip() +
  theme_cowplot() + theme(legend.position = 'none') +
  scale_colour_manual(values = fish.cols)+
  scale_fill_manual(values=fish.cols) +
  scale_y_continuous(limits=c(0,60),breaks=seq(0,60, by = 10), labels=seq(0, 60, by = 10)) +
  scale_x_discrete(name="",labels=rev(c("herbivore/detritivore","invertivore (mobile)",
                                        "piscivore","planktivore","invertivore (sessile)")))+
  labs(x = '', y = "Mean proportion of nutrient yield, %") 

g1

jpeg("FigureS9.jpeg", res=300, width=2600, height=2400)
g1
graphics.off()



#FIGURE S8
#RELATIVE BIOMASS ALONG A BIOMASS GRADIENT
setwd("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value")
load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/nutrient_US_Diets.RData")
head(nutrient_US_Diets) 

load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/finaldataFeb2022.RData") 
data <- as.data.frame(finaldataFeb2022)

data$logbm <- log(data$bm_kg_ha_US)
data$logzn <- log(data$Zinc_mg_m2_US)
data$logir <- log(data$Iron_mg_m2_US)
data$logca <- log(data$Calcium_mg_m2_US)
data$logom <- log(data$Omega_3_g_m2_US)
data$logva <- log(data$Vitamin_A_ug_m2_US)
data$logse <- log(data$Selenium_ug_m2_US)

nutyields <- data %>% select(UniqueSite,logzn:logse)
save(nutyields,file='nutyields.RData')

fdata <- finaldata %>% select(UniqueSite,Larger)

#ADD region
load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/datMEOW.RData")
dd <- unique(datMEOW[,c("Larger","REALM")])
tomerge <- merge(fdata,dd,by="Larger")

nut <- merge(nutrient_US_Diets,tomerge,by="UniqueSite")

#FIGURE S8
#RELATIVE BIOMASS ALONG A BIOMASS GRADIENT
datbm <- nut %>% select(UniqueSite,Larger,REALM,Diets,biomass) %>%
  group_by(UniqueSite) %>%
  mutate(bm_tot = sum(biomass)) %>% #compute total bm for each site
  ungroup() %>%
  mutate(propbm = biomass/bm_tot) 

pal <- brewer.pal(6, "PuBuGn")
#myPalette <- colorRampPalette(brewer.pal(6, "PuBuGn"))

fish.cols<-c('herbdetri'=pal[2], 
             'invertmobile'=pal[3], 
             'piscivore'=pal[4],
             'planktivore'=pal[5],
             'invertsessile'=pal[6])
#Plot
white_themejpg <-theme(axis.ticks=element_line(colour="black"),
                       axis.text=element_text(size=12,colour="black"),
                       axis.title=element_text(size=14),
                       panel.grid.minor=element_blank(),
                       panel.background=element_rect(fill="white",colour="black"),
                       plot.background=element_rect(fill="transparent",colour=NA),
                       legend.key = element_rect(fill = "white"))

datbm$logbm <- log(datbm$bm_tot)

#Add protection
tomerge <- finaldata[,c("UniqueSite","Protection")]

dat <- merge(datbm,tomerge,by="UniqueSite",all.x=T)



ggplot(dat, aes(x=logbm, y=propbm, colour = Diets)) + 
  geom_vline(xintercept=log(500),linetype="dashed",colour="black",size=0.5)+
  geom_point(aes(colour = Diets),size=1,alpha=0.4) +
  geom_smooth(method="glm", method.args = list(family = "quasipoisson"), formula = y ~ poly(x, 3), se = TRUE, fullrange = FALSE,
              size = 1, 
              aes(colour = Diets,fill=Diets))+
  scale_fill_manual(values = fish.cols)+
  scale_colour_manual(values = fish.cols)+
  scale_x_continuous(name = "LOG (total fish biomass)",limits=c(.7,10.1),breaks = c(2,5,7.5,10)) + 
  scale_y_continuous(name = "Relative biomass (%)", limits = c(0,1.1), 
                     breaks = c(0,.25,.50,.75,1),labels = c("0","25","50","75","100")) + 
  white_themejpg+  facet_wrap(~ REALM)







