######################
## LIBRARY PACKAGES ##
######################

library(bayestestR)
library(brms)
library(bayesplot)
library(parallel)
library(rstan)
library(tidyverse)
library(dplyr)
library(here)
library(patchwork)
library(rethinking)

here()

##############################
## IMPORT ALL FITTED MODELS ##
##############################

biomass_model <- readRDS("outputs/biomass_model.rds")
biomass_post <- as.data.frame(as.matrix(biomass_model)) %>%
  select('b_biomass')

PC1_model <- readRDS("outputs/PC1_model.rds")
PC1_post <- as.data.frame(as.matrix(PC1_model)) %>%
  select('b_PC1')

PC2_model <- readRDS("outputs/PC2_model.rds")
PC2_post <- as.data.frame(as.matrix(PC2_model)) %>%
  select('b_PC2')

fishdiversity_model <- readRDS("outputs/fishdiversity_model.rds")
fishdiversity_post <- as.data.frame(as.matrix(fishdiversity_model)) %>%
  select('b_fishdiversity')

gravity_model <- readRDS("outputs/gravity_model.rds")
gravity_post <- as.data.frame(as.matrix(gravity_model)) %>%
  select('b_gravity')

depth_model <- readRDS("outputs/depth_model.rds")
depth_post <- as.data.frame(as.matrix(depth_model)) %>%
  select('b_depth>10m','b_depth0M4m')

NPP_model <- readRDS("outputs/NPP_model.rds")
NPP_post <- as.data.frame(as.matrix(NPP_model)) %>%
  select('b_NPP')

meantemp_model <- readRDS("outputs/meantemp_model.rds")
meantemp_post <- as.data.frame(as.matrix(meantemp_model)) %>%
  select('b_meantemp')

HDI_model <- readRDS("outputs/HDI_model.rds")
HDI_post <- as.data.frame(as.matrix(HDI_model)) %>%
  select('b_HDI')

MPA_model <- readRDS("outputs/MPA_model.rds")
MPA_post <- as.data.frame(as.matrix(MPA_model)) %>%
  select('b_MPARestricted','b_MPAUnfishedHigh')

dhw_model <- readRDS("outputs/dhw_model.rds")
dhw_post <- as.data.frame(as.matrix(dhw_model)) %>%
  select('b_maxdhw')

bwsize_model <- readRDS("outputs/bwsize_model.rds")
bwsize_post <- as.data.frame(as.matrix(bwsize_model)) %>%
  select('b_bw_size')

voice_model <- readRDS("outputs/voice_model.rds")
voice_post <- as.data.frame(as.matrix(voice_model)) %>%
  select('b_voice')

geomorphology_model <- readRDS("outputs/geomorphology_model.rds")
geomorphology_post <- as.data.frame(as.matrix(geomorphology_model)) %>%
  select('b_geomorphologyCrest','b_geomorphologyFlat','b_geomorphologyLagoon_Backreef')

wave_energy_model <- readRDS("outputs/wave_energy_model.rds")
wave_energy_post <- as.data.frame(as.matrix(wave_energy_model)) %>%
  select('b_wave_energy')

###############################################
## RECOMBINE DAG MODELS AND RENAME VARIABLES ##
###############################################

dag_output <- data.frame(biomass_post,PC1_post,PC2_post,
                         fishdiversity_post,gravity_post,
                         depth_post,NPP_post,meantemp_post,
                         HDI_post,MPA_post,dhw_post,bwsize_post,voice_post,
                         geomorphology_post,wave_energy_post)

names(dag_output) <- gsub("b_", "", names(dag_output))

dag_output <- dag_output %>%
  rename('Fish biomass'=biomass ,
         "Trophic composition (PC1): HD (+) / PK (-)"=PC1,
         "Trophic composition (PC2): PS (+) / IM (-)"=PC2,
         'Species richness'=fishdiversity  ,
         'Human gravity'=gravity ,
         'Deep reef: >10m'=depth.10m ,
         'Shallow reef: 0-4m'=depth0M4m ,      
         'Ocean productivity'=NPP  ,
         'SST (mean)'=meantemp ,
         'Restricted fishing'=MPARestricted ,
         'Marine reserve'=MPAUnfishedHigh  ,     
         'DHW'=maxdhw    ,
         'Fish size'=bw_size,
         'Voice and accountability'=voice ,
         'Reef flat'=geomorphologyFlat ,
         'Lagoon'=geomorphologyLagoon_Backreef ,
         'Reef crest'=geomorphologyCrest,
         'Wave energy' = wave_energy)

dag_estimates <- data.frame(median=apply(dag_output, 2, median))
dag_estimates$abs_effect <- abs(dag_estimates$median)
dag_estimates <- dag_estimates %>%
  arrange(desc(abs_effect))

dag_output <- dag_output[,order(match(colnames(dag_output), rownames(dag_estimates)))]

post <- dag_output

test <- mcmc_intervals_data(dag_output, prob = 0.5,
                            prob_outer = 0.9,
                            point_est = c("median") )

HPDI_1 <- data.frame(test)
HPDI_1 <- HPDI_1 %>% select(!c(outer_width:point_est))
names(HPDI_1)<-c("Var","ll0.05","ll0.25","median","ul0.75","ul0.95")
rownames(HPDI_1) <- HPDI_1$Var

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

roworder <- rev(c("Trophic composition (PC1): HD (+) / PK (-)","Fish size","Species richness","Fish biomass","Trophic composition (PC2): PS (+) / IM (-)", #Species composition
                  "Human gravity","Voice and accountability","HDI", 
                  "Marine reserve","Restricted fishing","Fished areas",#management
                  "Shallow reef: 0-4m","Deep reef: >10m","Mid-depth reef (4-10m)","Reef flat","Reef crest","Lagoon","Slope", #local env.
                  "DHW",#"SST (range)",
                  "SST (mean)","Wave energy","Ocean productivity"#large scale env
)) #Method

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
baseline <- c(grep("Slope",datFig$Var),grep("Mid-depth",datFig$Var),grep("Fished",datFig$Var))
datFig$shape[baseline] <- 22

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
  geom_vline(xintercept=11.5, lwd=0.5, lty=1)+
  geom_vline(xintercept=17.5, lwd=0.5, lty=1)+
  geom_hline(yintercept=0, lwd=0.5, lty=2)+
  geom_linerange(aes(x = Var,ymin = ll0.05, ymax = ul0.95),lwd=0.5)+
  geom_linerange(aes(x = Var,ymin = ll0.25,ymax = ul0.75),lwd=1.5)+
  geom_point(stat='identity', shape = datFig$shape,size=datFig$size, fill=datFig$col)  +
  scale_y_continuous(breaks=c(-4,-3,-2,-1,0,1,2,3,4,5),limits=c( (min(datFig$ll0.05)+min(datFig$ll0.05)/6), max(datFig$ul0.95)+max(datFig$ul0.95)/6),name="Standardised effect on micronutrient density")+
  scale_x_discrete(name="")+
  coord_flip()+
  white_themejpg

effect

jpeg("figures/Figure2.jpeg", width=4000, height=2600, res=300) 
effect
graphics.off()

#END
