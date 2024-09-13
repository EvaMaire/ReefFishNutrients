#######################################################################################
#'  BAYSESIAN HIERARCHICAL MODELS TO ASSESS CAUSAL DRIVERS OF NUTRIENT DENSITY SCORE
#'  BASED ON DIRECTED ACYCLIC GRAPHS (DAG)
#'  THE CORRESPONDING DAG CAN BE FOUND AT http://dagitty.net/dags.html?id=g8iVuO
########################################################################################

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

rstan_options(auto_write = TRUE)

here()

###############################
##   IMPORT AND CLEAN DATA   ##
###############################

load('data/reefdata.RData')
dat <- reefdata %>% filter(nutrient_score>0) #remove site with nutrient_score=0 (zero target fish biomass)

#log-transform data
dat$human_gravity[which(is.na(dat$human_gravity))] <- 0 # gravity is assumed to be 0 when there is no population within 500km-buffer 
dat$loggrav <- log(dat$human_gravity+1)
dat$logNPP <- log(dat$meanNPP)
dat$logbm <- log(dat$bm_kg_ha+1)

#Replace NAs
dat$wave_energy[which(is.na(dat$wave_energy)==T)] <- mean(dat$wave_energy[which(is.na(dat$wave_energy)==F)])

#Factors
dat$depth<-as.factor(dat$depth)
dat$CensusMethod <- as.factor(dat$CensusMethod)
dat$MPA<-as.factor(dat$MPA)
dat$geomorphology<-as.factor(dat$geomorphology)

dat$depth=relevel(dat$depth,ref="4-10m")
dat$CensusMethod=relevel(dat$CensusMethod,ref="Standard belt transect")
dat$MPA=relevel(dat$MPA,ref="Fished")
dat$geomorphology=relevel(dat$geomorphology,ref="Slope")

#rename
bayes_model <- dat %>% rename(NPP = logNPP,
                              gravity = loggrav,
                              biomass = logbm,
                              voice = Voice_accountability)


summary(bayes_model)

#check missing data
sapply(bayes_model, function(x) paste(round(sum(is.na(x))/length(x),2)*100,"%",sep=""))

#standardise
bayes_model <- dplyr::mutate(bayes_model,
                             abs_latitude = abs(Site_Lat))

bayes_model <- bayes_model %>% dplyr::select(UniqueSite,ReefCluster,Larger,CensusMethod,
                                             abs_latitude,MPA,geomorphology,depth,nutrient_score,
                                             Total_sampling_area,NPP,HDId,fishdiversity,voice,wave_energy,biomass,
                                             maxdhw,meantemp,bw_size,gravity,PC1,PC2)


#standardised metric for species richness
bayes_model$fishdiversityst <- bayes_model$fishdiversity/(log(bayes_model$Total_sampling_area))

bayes_model[,c(10:ncol(bayes_model))] <- sapply(bayes_model[,c(10:ncol(bayes_model))], function(i){(i-mean(i, na.rm=T))/(2*sd(i, na.rm=T))} )
summary(bayes_model)

##################################
## PARALLEL SETTINGS FOR MODELS ##
##################################

ncores = detectCores()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###################
## BIOMASS MODEL ##
###################

biomass_model_formula <- 
  bf(nutrient_score ~ biomass +
       bw_size +
       PC1+PC2+ 
       meantemp + 
       CensusMethod + Total_sampling_area +
       (1 | Larger/ReefCluster),
     family=gaussian()) 

biomass_model <- brm(biomass_model_formula,
                     data=bayes_model,
                     chains=3, iter=2000, cores=ncores,
                     c(set_prior("normal(0,3)", class = "b"),
                       set_prior("normal(0,3)", class="Intercept")))

saveRDS(biomass_model, "outputs/biomass_model.rds")

##############################
## FISH COMPO MODEL 1 - PC1 ##
##############################

PC1_model_formula <- 
  bf(nutrient_score ~ PC1 +
       biomass +
       bw_size +
       meantemp +
       (1 | Larger/ReefCluster),
     family=gaussian()) 

PC1_model <- brm(PC1_model_formula,
                 data=bayes_model,
                 chains=3, iter=2000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(PC1_model, "outputs/PC1_model.rds")

##############################
## FISH COMPO MODEL 2 - PC2 ##
##############################

PC2_model_formula <- bf(nutrient_score ~ PC2 +
       biomass +
       bw_size +
       meantemp+
       (1 | Larger/ReefCluster),
     family=gaussian()) 

PC2_model <- brm(PC2_model_formula,
                 data=bayes_model,
                 chains=3, iter=2000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))


saveRDS(PC2_model, "outputs/PC2_model.rds")

#######################
##  DIVERSITY MODEL  ##
#######################

fishdiversity_model_formula <- bf(nutrient_score ~ fishdiversityst +
                                      CensusMethod +
                                      MPA + NPP + depth + maxdhw + geomorphology + meantemp +
                                      (1 | Larger/ReefCluster),
                                    family=gaussian()) 

fishdiversity_model <- brm(fishdiversity_model_formula,
                           data=bayes_model,
                           chains=3, iter=2000, cores=ncores,
                           c(set_prior("normal(0,3)", class = "b"),
                             set_prior("normal(0,3)", class="Intercept")))

saveRDS(fishdiversity_model, "outputs/fishdiversity_model.rds")

###################
## GRAVITY MODEL ##
###################

gravity_model_formula <- bf(nutrient_score ~ gravity + 
                              #MPA +
                              abs_latitude  + # ABSOLUTE VALUE FOR LINEARITY
                              (1 | Larger/ReefCluster),
                            family=gaussian())

gravity_model <- brm(gravity_model_formula,
                     data=bayes_model,
                     chains=3, iter=2000, cores=ncores,
                     c(set_prior("normal(0,3)", class = "b"),
                       set_prior("normal(0,3)", class="Intercept")))

saveRDS(gravity_model, "outputs/gravity_model.rds")

#################
## DEPTH MODEL ##
#################

depth_model_formula <- bf(nutrient_score ~ depth +
                            (1 | Larger/ReefCluster),
                          family=gaussian())

depth_model <- brm(depth_model_formula,
                   data=bayes_model,
                   chains=3, iter=2000, cores=ncores,
                   c(set_prior("normal(0,3)", class = "b"),
                     set_prior("normal(0,3)", class="Intercept")))

saveRDS(depth_model, "outputs/depth_model.rds")

###############
## NPP MODEL ##
###############

NPP_model_formula <- bf(nutrient_score ~ NPP + 
                          abs_latitude + wave_energy + 
                          (1 | Larger/ReefCluster),
                        family=gaussian())

NPP_model <- brm(NPP_model_formula,
                 data=bayes_model,
                 chains=3, iter=2000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(NPP_model, "outputs/NPP_model.rds")

################
## SST  MODEL ##
################

meantemp_model_formula <- bf(nutrient_score ~ meantemp + 
                               abs_latitude +
                               (1 | Larger/ReefCluster), 
                             family=gaussian())

meantemp_model <- brm(meantemp_model_formula,
                      data=bayes_model,
                      chains=3, iter=2000, cores=ncores,
                      c(set_prior("normal(0,3)", class = "b"),
                        set_prior("normal(0,3)", class="Intercept")))

saveRDS(meantemp_model, file= "outputs/meantemp_model.rds")

#################
## HDI   MODEL ##
#################

HDI_model_formula <- bf(nutrient_score ~ HDId + 
                          abs_latitude + # ABSOLUTE VALUE FOR LINEARITY
                          voice + 
                          (1 | Larger/ReefCluster),
                        family=gaussian())

HDI_model <- brm(HDI_model_formula,
                 data=bayes_model,
                 chains=3, iter=2000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(HDI_model, "outputs/HDId_model.rds")

###############
## MPA MODEL ##
###############

MPA_model_formula  <- bf(nutrient_score ~ MPA + 
                           abs_latitude +
                           (1 | Larger/ReefCluster),
                         family=gaussian())

MPA_model <- brm(MPA_model_formula,
                 data=bayes_model,
                 chains=3, iter=2000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(MPA_model, "outputs/MPA_model.rds")

###############
## DHW MODEL ##
###############

dhw_model_formula <- bf(nutrient_score ~ maxdhw +
                          (1 | Larger/ReefCluster),
                        family=gaussian())

dhw_model <- brm(dhw_model_formula,
                 data=bayes_model,
                 chains=3, iter=2000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(dhw_model, "outputs/dhw_model.rds")

###################
## BW SIZE MODEL ##
###################

bwsize_model_formula <- bf(nutrient_score ~ bw_size +
                             MPA + depth +
                             (1 | Larger/ReefCluster),
                           family=gaussian())

bwsize_model <- brm(bwsize_model_formula,
                    data=bayes_model,
                    chains=3, iter=2000, cores=ncores,
                    c(set_prior("normal(0,3)", class = "b"),
                      set_prior("normal(0,3)", class="Intercept")))

saveRDS(bwsize_model, "outputs/bwsize_model.rds")

#################
## VOICE MODEL ##
#################

voice_model_formula <- bf(nutrient_score ~ voice + 
                            abs_latitude +
                            (1 | Larger/ReefCluster),
                          family=gaussian())

voice_model <- brm(voice_model_formula,
                   data=bayes_model,
                   chains=3, iter=2000, cores=ncores,
                   c(set_prior("normal(0,3)", class = "b"),
                     set_prior("normal(0,3)", class="Intercept")))

saveRDS(voice_model, "outputs/voice_model.rds")

#######################
##   GEOMORPHO MODEL ##
#######################

geomorphology_model_formula <- bf(nutrient_score ~ geomorphology +
                                    depth +
                                    (1 | Larger/ReefCluster),
                                  family=gaussian())

geomorphology_model <- brm(geomorphology_model_formula,
                           data=bayes_model,
                           chains=3, iter=2000, cores=ncores,
                           c(set_prior("normal(0,3)", class = "b"),
                             set_prior("normal(0,3)", class="Intercept")))

saveRDS(geomorphology_model, "outputs/geomorphology_model.rds")

#######################
## WAVE ENERGY MODEL ##
#######################

wave_energy_model_formula <- bf(nutrient_score ~  wave_energy +
                                  depth + geomorphology +
                                  (1 | Larger/ReefCluster),
                                family=gaussian()) 

wave_energy_model <- brm(wave_energy_model_formula,
                         data=bayes_model,
                         chains=3, iter=2000, cores=ncores,
                         c(set_prior("normal(0,3)", class = "b"),
                           set_prior("normal(0,3)", class="Intercept")))

saveRDS(wave_energy_model, "outputs/wave_energy_model.rds")

############################
## CHECK DATA CONSISTENCY ##
############################

library(dagitty)

#rename variables to fit DAG: use only PC1 for fish composition
dag <- bayes_model %>% select(!c(PC2,fishdiversity)) %>% 
  rename(fishcompo = PC1,
         Latitude = abs_latitude, 
         dhw = maxdhw,HDI=HDId,
         Management_MPA = MPA,
         fishdiversity = fishdiversityst)

#Factors: must be converted as integers
dag$Management_MPA <- as.integer(dag$Management_MPA)
dag$depth <- as.integer(dag$depth)
dag$geomorphology <- as.integer(dag$geomorphology)

summary(dag)

# DOWNLOAD THE DAG 
# make sure your online DAG has the same variable names as your data; 
# make sure all unmeasured variables are labelled as unmeasured in the online DAG

DAG <- downloadGraph("dagitty.net/mW_c-oW") #september 19, 2023

dag <- dag[,names(dag) %in% names(DAG)]
str(dag)

length(names(dag))
length(names(DAG))

names(DAG)[!names(DAG) %in% names(dag)] #OK 

#evaluate the d-separation implications of the DAG
test <- localTests(DAG,dag) 

# turn results into data frame
test2 <- data.frame(test)

# create subset of data with a 0.5 correlation cut-off
testf <- subset(test2, estimate >= 0.5 | estimate <=-0.5)

# show final independencies that failed (NONE)
testf 

# END

