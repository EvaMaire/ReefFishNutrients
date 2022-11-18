setwd("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based")

######################
## LIBRARY PACKAGES ##
######################

library(data.table)
library(mapproj)
library(maptools)
library(plot3D)
library(pals)
#library(MuMIn)
library(bayestestR)
library(insight)
library(brms)
library(bayesplot)
library(performance)
library(parallel)
library(rstan)
library(readr)
library(tidyverse)
library(dplyr)
library(patchwork)
library(rethinking)


##############################
## IMPORT ALL FITTED MODELS ##
##############################

full_model <- readRDS("output/full_model.rds")
full_post <- as.data.frame(as.matrix(full_model)) %>%
  select('b_biomass':'b_MPAUnfishedHigh')
full_estimates <- data.frame(median=apply(full_post, 2, median))
full_estimates$abs_effect <- abs(full_estimates$median)
full_estimates <- full_estimates %>%
  arrange(desc(abs_effect))
full_post <- full_post[,order(match(colnames(full_post), rownames(full_estimates)))]

biomass_model <- readRDS("output/biomass_model.rds")
biomass_post <- as.data.frame(as.matrix(biomass_model)) %>%
  select('b_biomass')

#fishcompo_model <- readRDS("output/fishcompo_model.rds")
#fishcompo_post <- as.data.frame(as.matrix(fishcompo_model)) %>%
#select('b_fishcompo')
PC1_model <- readRDS("output/PC1_model.rds")
PC1_post <- as.data.frame(as.matrix(PC1_model)) %>%
  select('b_PC1')

PC2_model <- readRDS("output/PC2_model.rds")
PC2_post <- as.data.frame(as.matrix(PC2_model)) %>%
  select('b_PC2')

fishdiversity_model <- readRDS("output/fishdiversity_model.rds")
fishdiversity_post <- as.data.frame(as.matrix(fishdiversity_model)) %>%
  select('b_fishdiversity')

gravity_model <- readRDS("output/gravity_model.rds")
gravity_post <- as.data.frame(as.matrix(gravity_model)) %>%
  select('b_gravity')

depth_model <- readRDS("output/depth_model.rds")
depth_post <- as.data.frame(as.matrix(depth_model)) %>%
  select('b_depth>10m','b_depth0M4m')

NPP_model <- readRDS("output/NPP_model.rds")
NPP_post <- as.data.frame(as.matrix(NPP_model)) %>%
  select('b_NPP')

meantemp_model <- readRDS("output/meantemp_model.rds")
meantemp_post <- as.data.frame(as.matrix(meantemp_model)) %>%
  select('b_meantemp')

rangetemp_model <- readRDS("output/rangetemp_model.rds")
rangetemp_post <- as.data.frame(as.matrix(rangetemp_model)) %>%
  select('b_rangetemp')

HDI_model <- readRDS("output/HDI_model.rds")
HDI_post <- as.data.frame(as.matrix(HDI_model)) %>%
  select('b_HDI')

MPA_model <- readRDS("output/MPA_model.rds")
MPA_post <- as.data.frame(as.matrix(MPA_model)) %>%
  select('b_MPARestricted','b_MPAUnfishedHigh')

dhw_model <- readRDS("output/dhw_model.rds")
dhw_post <- as.data.frame(as.matrix(dhw_model)) %>%
  select('b_maxdhw')

bwsize_model <- readRDS("output/bwsize_model.rds")
bwsize_post <- as.data.frame(as.matrix(bwsize_model)) %>%
  select('b_bw_size')

voice_model <- readRDS("output/voice_model.rds")
voice_post <- as.data.frame(as.matrix(voice_model)) %>%
  select('b_voice')

#coralcover_model <- readRDS("output/coralcover_model.rds")
#coralcover_post <- as.data.frame(as.matrix(coralcover_model)) %>%
#select('bsp_nutrientscore_miCoralCover')

geomorphology_model <- readRDS("output/geomorphology_model.rds")
geomorphology_post <- as.data.frame(as.matrix(geomorphology_model)) %>%
  select('b_geomorphologyCrest','b_geomorphologyFlat','b_geomorphologyLagoon_Backreef')

wave_energy_model <- readRDS("output/wave_energy_model.rds")
wave_energy_post <- as.data.frame(as.matrix(wave_energy_model)) %>%
  select('b_wave_energy')

###############################################
## RECOMBINE DAG MODELS AND RENAME VARIABLES ##
###############################################

dag_output <- data.frame(biomass_post,PC1_post,PC2_post,#fishcompo_post,
                         fishdiversity_post,gravity_post,
                         depth_post,NPP_post,meantemp_post,rangetemp_post,
                         HDI_post,MPA_post,dhw_post,bwsize_post,voice_post,#coralcover_post,
                         geomorphology_post,wave_energy_post)

#names(full_post) <- gsub("b_", "", names(full_post))
names(dag_output) <- gsub("b_", "", names(dag_output))

dag_output <- dag_output %>%
  rename('Standing biomass'=biomass ,
         "PC1: herbivore/detritivore"=PC1,
         "PC2: piscivore (+) / invertivore (-)"=PC2,
         #'herb/detri'=fishcompo ,
         'Species richness'=fishdiversity  ,
         'Human gravity'=gravity ,
         'Deep reef: >10m'=depth.10m ,
         'Shallow reef: 0-4m'=depth0M4m ,      
         'Ocean productivity'=NPP  ,
         'SST (mean)'=meantemp ,
         'SST (range)'=rangetemp ,       
         'Restricted fishing'=MPARestricted ,
         'Marine reserve'=MPAUnfishedHigh  ,     
         'DHW'=maxdhw    ,
         'Mean fish size'=bw_size,
         'Voice and accountability'=voice ,
         #'coral cover'= bsp_nutrientscore_miCoralCover ,
         'Reef flat'=geomorphologyFlat ,
         'Lagoon'=geomorphologyLagoon_Backreef ,
         'Reef crest'=geomorphologyCrest,
         'Wave energy' = wave_energy)

dag_estimates <- data.frame(median=apply(dag_output, 2, median))
dag_estimates$abs_effect <- abs(dag_estimates$median)
dag_estimates <- dag_estimates %>%
  arrange(desc(abs_effect))

dag_output <- dag_output[,order(match(colnames(dag_output), rownames(dag_estimates)))]

p <- mcmc_intervals(dag_output) +
  #ggtitle("DAG-Based Models - Nutrient score calcium/iron/zinc") +
  xlim(range(dag_output))

tiff("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/Figure2.tiff", width=2000, height=2600, compression="lzw", res=300) 
p
graphics.off()

#Other option
post <- dag_output

test <- mcmc_intervals_data(dag_output, prob = 0.5,
  prob_outer = 0.9,
  point_est = c("median") )

##p_median <- apply(post,2,median)
#p_HPDI_1_95 <- hdi(post, prob = c(.5, .9)) #apply(post,2,HPDI,prob=0.9)
#p_HPDI_1_50 <- apply(post,2,HPDI,prob=0.5)
HPDI_1 <- data.frame(test)
HPDI_1 <- HPDI_1 %>% select(!c(outer_width:point_est))
names(HPDI_1)<-c("Var","ll0.05","ll0.25","median","ul0.75","ul0.95")
rownames(HPDI_1) <- HPDI_1$Var
#HPDI_1$Var <- rownames(HPDI_1)

HPDI_1$abs <- abs(HPDI_1$median)
HPDI_1$zero50 <- ifelse(HPDI_1$ll0.25<0 & HPDI_1$ul0.75 >0,1,0)
HPDI_1$zero95 <- ifelse(HPDI_1$ll0.05<0 & HPDI_1$ul0.95 >0,1,0)
HPDI_1$positive <- ifelse(HPDI_1$median >0,1,0)
HPDI_1$strong <- rep(0,nrow(HPDI_1))
HPDI_1$strong[which(HPDI_1$zero50 == 0 & HPDI_1$zero95 == 0)] <- 1
HPDI_1$light <- rep(0,nrow(HPDI_1))
HPDI_1$light[which(HPDI_1$zero50 == 0 & HPDI_1$zero95 == 1)] <- 1
HPDI_1$zero <- rep(0,nrow(HPDI_1))
HPDI_1$zero[which(HPDI_1$zero50 == 1 & HPDI_1$zero95 == 1)] <- 1
HPDI_1

#create matrix for reference levels (Slope, Fished Areas and average depth)

ref <- HPDI_1[c(1:3),]
rownames(ref) = c("Mid-depth reef (4-10m)","Slope","Fished areas")
ref$Var <- c("Mid-depth reef (4-10m)","Slope","Fished areas")
ref[,c(2:13)] <- 0

HPDI <- rbind(HPDI_1,ref)

roworder <- rev(c("PC1: herbivore/detritivore","Mean fish size","Species richness","Standing biomass","PC2: piscivore (+) / invertivore (-)", #Species composition
                  "Human gravity","Voice and accountability","HDI", #Social
                  "Shallow reef: 0-4m","Deep reef: >10m","Mid-depth reef (4-10m)","Reef flat","Reef crest","Lagoon","Slope", #local env.
                  "DHW","SST (range)","SST (mean)","Wave energy","Ocean productivity",#large scale env
                  "Marine reserve","Restricted fishing","Fished areas")) #Method

datFig <- HPDI[match(roworder,row.names(HPDI)),]
datFig$Var <- factor(datFig$Var, levels = datFig$Var)

#Set up color and size vectors
datFig$col <- rep("white",nrow(datFig)) 
datFig$col[which(datFig$zero == 1)] <- "white"
datFig$col[which(datFig$positive==1 & datFig$strong==1)] <- "#006666" #dark green ( strong positive)
datFig$col[which(datFig$positive==0 & datFig$strong==1)] <- "#862d59" #dark red (strong negative)
datFig$col[which(datFig$positive == 1 & datFig$light == 1)] <- "#00666680" #light green
datFig$col[which(datFig$positive == 0 & datFig$light == 1)] <- "#862d5980"

datFig$size <- rep(4,nrow(datFig))
datFig$size[which(datFig$strong==1)] <- 6

datFig$shape <- rep(21,nrow(datFig))
datFig$shape[c(1,9,13)] <- 22

datFig

#PLOT
white_themejpg <-theme(axis.ticks=element_line(colour="black"),
                       axis.text=element_text(size=15,colour="black"),
                       axis.title=element_text(size=18),
                       panel.grid.minor=element_blank(),
                       panel.background=element_rect(fill="white",colour="black"),
                       plot.background=element_rect(fill="transparent",colour=NA),
                       legend.key = element_rect(fill = "white"))


effect <- ggplot(datFig ,aes(x=Var, y=median)) + 
  #annotate("rect", ymin = min(datFig$ll0.05)+min(datFig$ll0.05)/6, ymax = max(datFig$ul0.95)+max(datFig$ll0.05)/6, xmin = 14.5, xmax = 16.5, alpha = .7, fill="#ffa31a10")+ #Management
  #annotate("rect", ymin = min(datFig$ll0.05)+min(datFig$ll0.05)/6, ymax = max(datFig$ul0.95)+max(datFig$ll0.05)/6, xmin = 9.5, xmax = 14.5, alpha = .7,  fill="#002b8050")+ #Social Drivers
  #annotate("rect", ymin = min(datFig$ll0.05)+min(datFig$ll0.05)/6, ymax = max(datFig$ul0.95)+max(datFig$ll0.05)/6, xmin = 2.5, xmax = 9.5, alpha = .7,   fill="#40bf8050")+ #Env
  #annotate("rect", ymin = min(datFig$ll0.05)+min(datFig$ll0.05)/6, ymax = max(datFig$ul0.95)+max(datFig$ll0.05)/6, xmin = 0.5, xmax = 2.5, alpha = .7,   fill="#52527a50")+ #Method
  geom_vline(xintercept=3.5, lwd=0.5, lty=1)+
  geom_vline(xintercept=15.5, lwd=0.5, lty=1)+
  geom_vline(xintercept=18.5, lwd=0.5, lty=1)+
  geom_hline(yintercept=0, lwd=0.5, lty=2)+
  geom_linerange(aes(x = Var,ymin = ll0.05, ymax = ul0.95),lwd=0.5)+
  geom_linerange(aes(x = Var,ymin = ll0.25,ymax = ul0.75),lwd=1.5)+
  geom_point(stat='identity', shape = datFig$shape,size=datFig$size, fill=datFig$col)  +
  scale_y_continuous(breaks=c(-10,-5,0,5,10,15),limits=c( (min(datFig$ll0.05)+min(datFig$ll0.05)/6), max(datFig$ul0.95)+max(datFig$ul0.95)/6),name="Standardised effect on micronutrient density")+
  scale_x_discrete(name="")+
  coord_flip()+
  white_themejpg

effect

tiff("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/Figure2V2.tiff", width=2500, height=2600, compression="lzw", res=300) 
effect
graphics.off()

#END
