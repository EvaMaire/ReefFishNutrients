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

load("data/data_for_models.RData")

#remove 1 site with no target biomass (nutrient_score = 0)
hist(data_for_models$nutrient_score) 
bayes_model <- data_for_models %>% filter(nutrient_score >0) 
hist(bayes_model$nutrient_score) #normal distribution

#standardise all covariates
bayes_model[,c(9:ncol(bayes_model))] <- sapply(bayes_model[,c(9:ncol(bayes_model))], function(i){(i-mean(i, na.rm=T))/(2*sd(i, na.rm=T))} )
summary(bayes_model) #no missing values - OK

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
       meantemp+
       (1 | Larger/ReefCluster),
     family=gaussian()) 

biomass_model <- brm(biomass_model_formula,
                     data=bayes_model,
                     chains=3, iter=2000, cores=ncores,
                     c(set_prior("normal(0,3)", class = "b"),
                       set_prior("normal(0,3)", class="Intercept")))

saveRDS(biomass_model, "output/biomass_model.rds")

##############################
## FISH COMPO MODEL 1 - PC1 ##
##############################

PC1_model_formula <- 
  bf(nutrient_score ~ PC1 +
       fishdiversity+
       geomorphology +
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

PC2_model_formula <- 
  bf(nutrient_score ~ PC2 +
       fishdiversity+
       geomorphology +
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

fishdiversity_model_formula <-   bf(nutrient_score ~ fishdiversity +
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
                              HDI +
                              abs_latitude + # ABSOLUTE VALUE FOR LINEARITY
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
                          (1 | Larger/ReefCluster), 
                        family=gaussian())

meantemp_model <- brm(meantemp_model_formula,
                 data=bayes_model,
                 chains=3, iter=2000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(meantemp_model, file= "outputs/meantemp_model.rds")

######################
## SST RANGE  MODEL ##
######################

rangetemp_model_formula <- bf(nutrient_score ~ rangetemp +
                          (1 | Larger/ReefCluster), 
                        family=gaussian())

rangetemp_model <- brm(rangetemp_model_formula,
                 data=bayes_model,
                 chains=3, iter=2000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(rangetemp_model, file= "outputs/rangetemp_model.rds")

#################
## HDI   MODEL ##
#################

HDI_model_formula <- bf(nutrient_score ~ HDI + gravity +
                          abs_latitude + # ABSOLUTE VALUE FOR LINEARITY
                          (1 | Larger/ReefCluster),
                        family=gaussian())

HDI_model <- brm(HDI_model_formula,
                 data=bayes_model,
                 chains=3, iter=2000, cores=ncores,
                 c(set_prior("normal(0,3)", class = "b"),
                   set_prior("normal(0,3)", class="Intercept")))

saveRDS(HDI_model, "outputs/HDI_model.rds")

###############
## MPA MODEL ##
###############

MPA_model_formula  <- bf(nutrient_score ~ MPA +
                           HDI +
                           voice +
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
                            HDI +
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
dag <- bayes_model %>% select(!c(PC2)) %>% 
  rename(fishcompo = PC1,
         Latitude = abs_latitude, dhw = maxdhw)

#Factors: must be converted as integers
dag$MPA <- as.integer(dag$MPA)
dag$depth <- as.integer(dag$depth)
dag$geomorphology <- as.integer(dag$geomorphology)

summary(dag)

# DOWNLOAD THE DAG 
# make sure your online DAG has the same variable names as your data; 
# make sure all unmeasured variables are labelled as unmeasured in the online DAG

DAG <- downloadGraph("dagitty.net/mg8iVuO")

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

