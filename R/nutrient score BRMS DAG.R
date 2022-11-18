#######################################################################################
#'  BAYSESIAN HIERARCHICAL MODELS TO ASSESS CAUSAL DRIVERS OF NUTRIENT DENSITY SCORE
#'  BASED ON DIRECTED ACYCLIC GRAPHS (DAG)
#'  THE CORRESPONDING DAG CAN BE FOUND AT XXXXX
#'
#'     
#'
########################################################################################

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

#library(broom.mixed)
#library(janitor)

#library(hrbrthemes)

rstan_options(auto_write = TRUE)

# ADD PREDICTORS FOR IMPUTATIONS
# NOTE - IMPUTE ALL MISSING VALUES??
# CHECK CROSS VALIDATION FOR IMPUTATION WITH AND WITHOUT COVARIATES

###############################
##   IMPORT AND CLEAN DATA   ##
###############################

#load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/nutdata_outliers_removedFeb2022.RData") 
load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/finaldataFeb2022.RData") 
finaldata <- as.data.frame(finaldataFeb2022)
goodData <- finaldata %>% dplyr::select(UniqueSite:Larger,Protection:Total_sampling_area,Ocean_prod:HDI,
                                                       gravtotCluster:bm_kg_ha_US)

#load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/functiophylo.RData")

load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/nutscore_US_targetgroup.RData") 
d <- nutscore_US_targetgroup %>% filter(targetgroup=="target") %>% select(!biomass)
d[which(is.na(d$BWnut3A)==T),c(9:12)] <- 0

dd <- merge(goodData,d,by="UniqueSite")

#Add coral cover, voice 
last <- read.csv("/Volumes/EM2T/Lancaster/Nutrition/SERF/Raw SERF data for Eva/SERF DATA FOR GEORGIE_2015_11_18.csv",sep=";") %>% 
  dplyr::select(UniqueSite,Voice_accountability,CoralCover)

ddd <- merge(dd,last,by="UniqueSite")
#names(ddd)

load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/nutrient_US_Diets.RData")
bmdat <- nutrient_US_Diets %>% select(UniqueSite, Diets, biomass) %>% group_by(UniqueSite) %>%
  mutate(totbiomass = sum(biomass)) %>% dplyr::mutate(relbm = biomass/totbiomass) %>% select(UniqueSite,Diets,relbm) %>%
  pivot_wider(names_from = Diets, values_from = relbm)

#check relationships with target versus non target biomass and total biomass
load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/nutscore_US_targetgroup.RData")
target <- nutscore_US_targetgroup %>% select(UniqueSite:biomass) %>% 
  pivot_wider(names_from = targetgroup,values_from = biomass) %>%
  rename(biomasstarget = target,
         biomassnontarget = non_target)

dddd <- merge(ddd,target,by="UniqueSite",all.x=T)

#bmtotdat <- nutrient_US_Diets %>% select(UniqueSite, Diets, biomass) %>%
  #pivot_wider(names_from = Diets, values_from = biomass)

#explo <- merge(dd,bmtotdat,by="UniqueSite",all.x=T)
#explo$logbm <- log(explo$bm_kg_ha_US+1)
#explo$logherb <- log(explo$herbdetri+1)
#summary(explo) 
#explo <- explo %>% select(UniqueSite,logherb)

dat <- merge(dddd,bmdat,by="UniqueSite",all.x=T)
dat$logbm <- log(dat$bm_kg_ha_US+1)
dat$logbmtarget <- log(dat$biomasstarget+1)
dat$logbmnontarget <- log(dat$biomassnontarget+1)
summary(dat) 

#add temperature data
load("tempdata.RData")
tempdata <- tempdata %>% select(-dhw)
load('dhwdata.RData')
tdata <- merge(tempdata,dhwdata,by="ReefCluster")
dat <- merge(dat,tdata,by="ReefCluster")
summary(dat) 

#add biomass-weighted size of fish assemblages
load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/fishass_USlevel.RData")
tomerge <- fishass_USlevel %>% select(UniqueSite,meansize,meantl)
dat <- merge(dat,tomerge,by="UniqueSite")

#add wave energy
load("w_energy_data.RData")
dat <- merge(dat,w_energy_data,by="UniqueSite")

#Add Fish compo PCA
load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/fishcomm_pca.RData")
dat <- merge(dat,fishcomm_pca,by="UniqueSite")

#Add species richness
load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/newSR.RData")

load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/nutdatJuly2022.RData")
spr <- nutdatJuly2022 %>% select(UniqueSite,sp_richness)

spR <- merge(newSR,spr,by="UniqueSite") 
spR$sp_richness[which(is.na(spR$sp_richness)==T)] <- 0

plot(spR$sp_richness,spR$sp_richness_imp)

dat <- merge(dat,newSR,by="UniqueSite") 

#replace missing values
dat$gravtotCluster[which(is.na(dat$gravtotCluster))] <- 0
dat$loggrav <- log(dat$gravtotCluster+1)
dat$CoralCover[which(dat$CoralCover==-999)] <- NA

#create binomial response
dat$Increased_nutrient_score <- ifelse(dat$BWnut3A>median(dat$BWnut3A),1,0)

#save(dat,file='/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/dat.RData')

#rename
bayes_model <- dat %>% rename(MPA = Protection,
                              Latitude = Site_Lat,
                              depth = DepthCategory,
                              geomorphology = CleanHabitat,
                              NPP = Ocean_prod,
                              fishdiversity = sp_richness_imp,
                              gravity = loggrav,
                              meandhw = dhw,
                              maxdhw = mdhw,
                              #nutrient_score = BWnut3A,
                              bw_size= meansize,
                              fishcompo = herbdetri, #THIS SHOULD BE MODIFIED - USE PC1 AND PC2 INSTEAD
                              #Herbivore_Detritivore = logherb,
                              biomass = logbm,
                              voice = Voice_accountability) %>%
  mutate(nutrient_score = BWnut3A/3)

#Factors
bayes_model$depth<-as.factor(bayes_model$depth)
bayes_model$MPA<-as.factor(bayes_model$MPA)
bayes_model$geomorphology<-as.factor(bayes_model$geomorphology)

bayes_model$depth=relevel(bayes_model$depth,ref="4-10m")
bayes_model$MPA=relevel(bayes_model$MPA,ref="Fished")
bayes_model$geomorphology=relevel(bayes_model$geomorphology,ref="Slope")

summary(bayes_model)

#bayes_model$sst_range <- bayes_model$sst_max - bayes_model$sst_min

############################
## HOW MUCH MISSING DATA? ##
############################

sapply(bayes_model, function(x) paste(round(sum(is.na(x))/length(x),2)*100,"%",sep=""))

#standardise
bayes_model <- dplyr::mutate(bayes_model,
                             abs_latitude = abs(Latitude))

bayes_model <- bayes_model %>% select(UniqueSite,ReefCluster,Larger,abs_latitude,MPA,geomorphology,depth,nutrient_score,
                                      CoralCover,NPP,HDI,fishdiversity,voice,wave_energy,fishcompo,biomass,maxdhw,meantemp,rangetemp,bw_size,gravity,PC1,PC2)

bayes_model[,c(9:ncol(bayes_model))] <- sapply(bayes_model[,c(9:ncol(bayes_model))], function(i){(i-mean(i, na.rm=T))/(2*sd(i, na.rm=T))} )
summary(bayes_model)

#Replace NAs - 'mi' is not supported for family 'binomial(logit)'.
bayes_model$wave_energy[which(is.na(bayes_model$wave_energy)==T)] <- mean(bayes_model$wave_energy[which(is.na(bayes_model$wave_energy)==F)])
#bayes_model$CoralCover[which(is.na(bayes_model$CoralCover)==T)] <- mean(bayes_model$CoralCover[which(is.na(bayes_model$CoralCover)==F)])


hist(bayes_model$nutrient_score) #remove nutrient_score =0
bayes_model <- bayes_model[-which(bayes_model$nutrient_score==0),] #bayes_model[-which(bayes_model$nutrient_score==0),]
hist(bayes_model$nutrient_score)

#Correlogram
require(corrgram)
require(corrplot)

d <- bayes_model[,c(10:ncol(bayes_model))]
summary(d)
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

##################################
## PARALLEL SETTINGS FOR MODELS ##
##################################

ncores = detectCores()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

##########################
## INTERCEPT ONLY MODEL ##
##########################

int_only_formula <- bf(nutrient_score ~ 
                         (1 | Larger/ReefCluster),
                       family = gaussian())

int_only_model <- brm(int_only_formula,
                      data=bayes_model,
                      chains=3, iter=2000, cores=ncores,
                      refresh=500,
                      set_prior("normal(0,3)", class="Intercept"))

saveRDS(int_only_model, file="output/int_only_model.rds")
int_only_model <- readRDS("output/int_only_model.rds")

posterior_summary(int_only_model, pars = c("^b_", "^sd_", "sigma") )

library(brmstools)

forest(int_only_model, grouping = "Larger")

#test forest plot
#library(meta)
library(tidybayes)
library(ggridges)
library(glue)
library(stringr)
library(forcats)

study.draws <- spread_draws(int_only_model, r_Larger[Larger,], b_Intercept) %>% 
  mutate(b_Intercept = r_Larger + b_Intercept)

pooled.effect.draws <- spread_draws(int_only_model, b_Intercept) %>% 
  mutate(Larger = "Global average")

forest.data <- bind_rows(study.draws, 
                         pooled.effect.draws) %>% 
  ungroup() %>%
  mutate(Larger = str_replace_all(Larger, "[.]", " ")) %>% 
  mutate(Larger = reorder(Larger, b_Intercept))

forest.data.summary <- group_by(forest.data, Larger) %>% 
  mean_qi(b_Intercept)

ggplot(aes(b_Intercept, 
           relevel(Larger, "Global average", after = Inf)), data = forest.data) +
  
  # Add vertical lines for pooled effect and CI
  geom_vline(xintercept = fixef(int_only_model)[1, 1], color = "grey", size = 1) +
  geom_vline(xintercept = fixef(int_only_model)[1, 3:4], color = "grey", linetype = 2) +
  #geom_vline(xintercept = 0, color = "black", size = 1) +
  
  # Add densities
  geom_density_ridges(fill = "dark grey",  rel_min_height = 0.01,col = NA, scale = 1,alpha = 0.8) +
  #geom_pointintervalh(data = forest.data.summary, size = 1) +
  
  # Add text and labels
  geom_text(data = mutate_if(forest.data.summary,  is.numeric, round, 2),
            aes(label = glue("{b_Intercept} [{.lower}, {.upper}]"),  x = Inf), hjust = "inward") +
  labs(x = "Nutrient Density Score", y = element_blank()) +
  theme_minimal()

#######################################
#customised plot
dat <- study.draws %>% ungroup() %>%
  mutate(Larger = str_replace_all(Larger, "[.]", " ")) %>% 
  mutate(Larger = reorder(Larger, b_Intercept))

pint <- ggplot(aes(b_Intercept, Larger), data = dat) +
  
  # Add vertical lines for pooled effect and CI
  geom_vline(xintercept = fixef(int_only_model)[1, 1], color = "grey", size = .5,linetype=2) +
  #geom_vline(xintercept = fixef(int_only_model)[1, 3:4], color = "grey", linetype = 2) +
  #geom_vline(xintercept = 0, color = "black", size = 1) +
  
  # Add densities
  geom_density_ridges(fill = "dark grey",  rel_min_height = 0.01,col = NA, scale = 1,alpha = 0.8) +
  #geom_pointintervalh(data = forest.data.summary, size = 1) +
  
  # Add text and labels
  #geom_text(data = mutate_if(forest.data.summary,  is.numeric, round, 2),
            #aes(label = glue("{b_Intercept} [{.lower}, {.upper}]"),  x = Inf), hjust = "inward") +
  labs(x = "Nutrient Density Score", y = element_blank()) +
  theme_minimal()

pint

########################
## CAUSAL SALAD MODEL ##
########################

full_model_formula <- 
  bf(nutrient_score ~
       biomass +
       depth +
       geomorphology+
       fishdiversity +
       NPP +
       meantemp+
       rangetemp+
       maxdhw+
       bw_size+
       #mi(CoralCover) +
       wave_energy +
       HDI +
       gravity +
       voice +
       PC1+
       PC2+
       MPA +
       (1 | Larger/ReefCluster),
     family = gaussian()) #+ 
  
  #bf(CoralCover|mi() ~ 
       #(1|Larger/ReefCluster),
     #family=gaussian()) +

  #bf(wave_energy|mi() ~ 
      #(1|Larger/ReefCluster),
    #family=gaussian()) 

full_model <- brm(full_model_formula,
                  data=bayes_model,
                  chains=3, iter=2000, cores=ncores,
                  refresh=500,
                  c(set_prior("normal(0,3)", class = "b"),
                    set_prior("normal(0,3)", class="Intercept")))

saveRDS(full_model, "output/full_model.rds")
full_model <- readRDS("output/full_model.rds")

performance::r2_bayes(full_model) # 
full_post <- as.data.frame(as.matrix(full_model)) %>%
  select('b_biomass':'b_MPAUnfishedHigh')

full_estimates <- data.frame(median=apply(full_post, 2, median))
full_estimates$abs_effect <- abs(full_estimates$median)
full_estimates <- full_estimates %>%
  arrange(desc(abs_effect))

full_post <- full_post[,order(match(colnames(full_post), rownames(full_estimates)))]
mcmc_intervals(full_post)

pp_check(full_model)

#bayes.test <- diamonds.keep[sample(nrow(diamonds.keep), 20000), ]
pred.1 <- predict(full_model)#, newdata = bayes_model)

r.sq <- as.character(round(summary(lm(bayes_model$nutrient_score ~ pred.1[, 1]))$r.squared, 2))

lb1 <- paste("R^2 == ", r.sq)
ggplot()+  geom_point(aes(x = pred.1[,1], y = bayes_model$nutrient_score)) +
  geom_errorbarh(aes(x = pred.1[,1], y = bayes_model$nutrient_score, 
                     xmin = pred.1[,1] - pred.1[, 2], 
                     xmax = pred.1[,1] + pred.1[, 2])) + 
  geom_smooth(aes(x = pred.1[,1], y = bayes_model$nutrient_score), method = "lm", color = "red", lty = 2)+
  geom_text(aes(x=10, y=30, label = lb1, size = 8), parse=TRUE, show.legend = F) +
  xlab("Predicted") +
  ylab("Observed")

#compare intercept only model and full model
full.draws <- spread_draws(full_model, r_Larger[Larger,], b_Intercept) %>% 
  mutate(b_Intercept = r_Larger + b_Intercept)

pooled.effect.full.draws <- spread_draws(full_model, b_Intercept) %>% 
  mutate(Larger = "Global average")

forest.data.full <- bind_rows(full.draws, 
                         pooled.effect.full.draws) %>% 
  ungroup() %>%
  mutate(Larger = str_replace_all(Larger, "[.]", " ")) %>% 
  mutate(Larger = reorder(Larger, b_Intercept))

forest.data.full.summary <- group_by(forest.data.full, Larger) %>% 
  mean_qi(b_Intercept)

datfull <- full.draws %>% ungroup() %>%
  mutate(Larger = str_replace_all(Larger, "[.]", " ")) %>% 
  mutate(Larger = reorder(Larger, b_Intercept))

pfull <- ggplot(aes(b_Intercept, Larger), data = datfull) +
  
  # Add vertical lines for pooled effect and CI
  geom_vline(xintercept = fixef(full_model)[1, 1], color = "grey", size = .5,linetype=2) +
  #geom_vline(xintercept = fixef(full_model)[1, 3:4], color = "grey", linetype = 2) +
  #geom_vline(xintercept = 0, color = "black", size = 1) +
  
  # Add densities
  geom_density_ridges(fill = "dark grey",  rel_min_height = 0.01,col = NA, scale = 1,alpha = 0.8) +
  #geom_pointintervalh(data = forest.data.full.summary, size = 1) +
  
  # Add text and labels
  #geom_text(data = mutate_if(forest.data.full.summary,  is.numeric, round, 2),
  #aes(label = glue("{b_Intercept} [{.lower}, {.upper}]"),  x = Inf), hjust = "inward") +
  labs(x = "Nutrient Density Score", y = element_blank()) +
  ggtitle("Full model")+
  theme_minimal()

pfull

#rank int_only_model data to ease comparison
goodl <- levels(datfull$Larger)
dat2 <- study.draws %>% ungroup() %>%
  mutate(Larger = str_replace_all(Larger, "[.]", " ")) %>% 
  mutate(LargerOrd=fct_relevel(Larger,goodl)) 

pint2 <- ggplot(aes(b_Intercept, LargerOrd), data = dat2) +
  
  # Add vertical lines for pooled effect and CI
  geom_vline(xintercept = fixef(int_only_model)[1, 1], color = "grey", size = .5,linetype=2) +
  #geom_vline(xintercept = fixef(int_only_model)[1, 3:4], color = "grey", linetype = 2) +
  #geom_vline(xintercept = 0, color = "black", size = 1) +
  
  # Add densities
  geom_density_ridges(fill = "dark grey",  rel_min_height = 0.01,col = NA, scale = 1,alpha = 0.8) +
  #geom_pointintervalh(data = forest.data.summary, size = 1) +
  
  # Add text and labels
  #geom_text(data = mutate_if(forest.data.summary,  is.numeric, round, 2),
  #aes(label = glue("{b_Intercept} [{.lower}, {.upper}]"),  x = Inf), hjust = "inward") +
  labs(x = "Nutrient Density Score", y = element_blank()) +
  ggtitle("Null model")+
  theme_minimal()

pint2

(pfull + pint2)

jpeg("compare_full_versus_null_models_COUNTRIES_3nut.jpeg", res=300, width=4000, height=2000)
(pfull + pint2)
graphics.off()

#ADD region
load("datMEOW.RData")
tomerge <- unique(datMEOW[,c("Larger","REALM")])

datfull <- merge(datfull,tomerge,by="Larger")

pfullR <- ggplot(aes(b_Intercept, Larger), data = datfull) +
  
  # Add vertical lines for pooled effect and CI
  geom_vline(xintercept = fixef(full_model)[1, 1], color = "grey", size = .5,linetype=2) +
  #geom_vline(xintercept = fixef(full_model)[1, 3:4], color = "grey", linetype = 2) +
  #geom_vline(xintercept = 0, color = "black", size = 1) +
  
  # Add densities
  geom_density_ridges(fill = "dark grey",  rel_min_height = 0.01,col = NA, scale = 1,alpha = 0.8) +
  #geom_pointintervalh(data = forest.data.full.summary, size = 1) +
  
  # Add text and labels
  #geom_text(data = mutate_if(forest.data.full.summary,  is.numeric, round, 2),
  #aes(label = glue("{b_Intercept} [{.lower}, {.upper}]"),  x = Inf), hjust = "inward") +
  labs(x = "Nutrient Density Score", y = element_blank()) +
  scale_x_continuous(limits = c(20,100))+
  ggtitle("Full model")+
  theme_minimal()+
  facet_wrap(~REALM,scales = "free_y")

pfullR

###################
## BIOMASS MODEL ##
###################
#bw_size, PC1, PC2, meantemp

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
biomass_model <- readRDS("output/biomass_model.rds")

biomass_post <- as.data.frame(as.matrix(biomass_model)) %>%
  select('b_biomass')

##############################
## FISH COMPO MODEL 1 - PC1 ##
##############################
#fishdiversity , geomorphology

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

saveRDS(PC1_model, "output/PC1_model.rds")
PC1_model <- readRDS("output/PC1_model.rds")

PC1_post <- as.data.frame(as.matrix(PC1_model)) %>%
  select('b_PC1')

##############################
## FISH COMPO MODEL 2 - PC2 ##
##############################
#fishdiversity , geomorphology

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

saveRDS(PC2_model, "output/PC2_model.rds")
PC2_model <- readRDS("output/PC2_model.rds")

PC2_post <- as.data.frame(as.matrix(PC2_model)) %>%
  select('b_PC2')

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

saveRDS(fishdiversity_model, "output/fishdiversity_model.rds")
fishdiversity_model <- readRDS("output/fishdiversity_model.rds")

fishdiversity_post <- as.data.frame(as.matrix(fishdiversity_model)) %>%
  select('b_fishdiversity')

###################
## GRAVITY MODEL ##
###################
#HDI, Latitude
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

saveRDS(gravity_model, "output/gravity_model.rds")
gravity_model <- readRDS("output/gravity_model.rds")

gravity_post <- as.data.frame(as.matrix(gravity_model)) %>%
  select('b_gravity')

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

saveRDS(depth_model, "output/depth_model.rds")
depth_model <- readRDS("output/depth_model.rds")

depth_post <- as.data.frame(as.matrix(depth_model)) %>%
  select('b_depth>10m','b_depth0M4m')


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

saveRDS(NPP_model, "output/NPP_model.rds")
NPP_model <- readRDS("output/NPP_model.rds")

NPP_post <- as.data.frame(as.matrix(NPP_model)) %>%
  select('b_NPP')

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

saveRDS(meantemp_model, file= "output/meantemp_model.rds")
meantemp_model <- readRDS("output/meantemp_model.rds")

meantemp_post <- as.data.frame(as.matrix(meantemp_model)) %>%
  select('b_meantemp')

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

saveRDS(rangetemp_model, file= "output/rangetemp_model.rds")
rangetemp_model <- readRDS("output/rangetemp_model.rds")

rangetemp_post <- as.data.frame(as.matrix(rangetemp_model)) %>%
  select('b_rangetemp')

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

saveRDS(HDI_model, "output/HDI_model.rds")
HDI_model <- readRDS("output/HDI_model.rds")

HDI_post <- as.data.frame(as.matrix(HDI_model)) %>%
  select('b_HDI')

###############
## MPA MODEL ##
###############
#HDI, voice
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

saveRDS(MPA_model, "output/MPA_model.rds")
MPA_model <- readRDS("output/MPA_model.rds")

MPA_post <- as.data.frame(as.matrix(MPA_model)) %>%
  select('b_MPARestricted','b_MPAUnfishedHigh')


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

saveRDS(dhw_model, "output/dhw_model.rds")
dhw_model <- readRDS("output/dhw_model.rds")

dhw_post <- as.data.frame(as.matrix(dhw_model)) %>%
  select('b_maxdhw')

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

saveRDS(bwsize_model, "output/bwsize_model.rds")
bwsize_model <- readRDS("output/bwsize_model.rds")

bwsize_post <- as.data.frame(as.matrix(bwsize_model)) %>%
  select('b_bw_size')

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

saveRDS(voice_model, "output/voice_model.rds")
voice_model <- readRDS("output/voice_model.rds")

voice_post <- as.data.frame(as.matrix(voice_model)) %>%
  select('b_voice')


#######################
## CORAL COVER MODEL ##
#######################

#coralcover_model_formula <- bf(nutrient_score ~ mi(CoralCover) +
                           # (1 | Larger/ReefCluster),
                          #family=gaussian())+
  
  #bf(CoralCover|mi() ~ 
       #(1|Larger/ReefCluster),
     #family=gaussian())

#coralcover_model <- brm(coralcover_model_formula,
                   #data=bayes_model,
                   #chains=3, iter=2000, cores=ncores,
                  # c(set_prior("normal(0,3)", class = "b"),
                    # set_prior("normal(0,3)", class="Intercept")))

#saveRDS(coralcover_model, "output/coralcover_model.rds")
#coralcover_model <- readRDS("output/coralcover_model.rds")

#coralcover_post <- as.data.frame(as.matrix(coralcover_model)) %>%
  #select('bsp_nutrientscore_miCoralCover')

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

saveRDS(geomorphology_model, "output/geomorphology_model.rds")
geomorphology_model <- readRDS("output/geomorphology_model.rds")

geomorphology_post <- as.data.frame(as.matrix(geomorphology_model)) %>%
  select('b_geomorphologyCrest','b_geomorphologyFlat','b_geomorphologyLagoon_Backreef')

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

saveRDS(wave_energy_model, "output/wave_energy_model.rds")
wave_energy_model <- readRDS("output/wave_energy_model.rds")

wave_energy_post <- as.data.frame(as.matrix(wave_energy_model)) %>%
  select('b_wave_energy')

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

names(full_post) <- gsub("b_", "", names(full_post))
names(dag_output) <- gsub("b_", "", names(dag_output))

dag_output <- dag_output %>%
  rename('log biomass'=biomass ,
         "PC1: herbivore/detritivore"=PC1,
         "PC2: piscivore (+) / invertivore (-)"=PC2,
         #'herb/detri'=fishcompo ,
         'species richness'=fishdiversity  ,
         'gravity'=gravity ,
         'deep >10m'=depth.10m ,
         'shallow 0-4m'=depth0M4m ,      
         'NPP'=NPP  ,
         'mean SST'=meantemp ,
         'range SST'=rangetemp ,       
         'Restricted fishing'=MPARestricted ,
         'Reserve'=MPAUnfishedHigh  ,     
         'DHW'=dhw    ,
         'assemblage size'=bw_size,
         'voice and acc'=voice ,
         #'coral cover'= bsp_nutrientscore_miCoralCover ,
         'Flat'=geomorphologyFlat ,
         'Lagoon'=geomorphologyLagoon_Backreef ,
         'Crest'=geomorphologyCrest,
         'wave energy' = wave_energy)


full_post <- full_post %>%
  rename('log biomass'= biomass ,
         "PC1: herbivore/detritivore"= PC1 ,
         "PC2: piscivore (+) / invertivore (-)"= PC2 ,
         'species richness'= fishdiversity  ,
         'gravity'= gravity ,
         'deep >10m'='depth>10m' ,
         'shallow 0-4m'= depth0M4m ,      
         'NPP'= NPP  ,
         'mean SST'= meantemp ,
         'range SST'= rangetemp ,       
         'Restricted fishing'= MPARestricted ,
         'Reserve'= MPAUnfishedHigh  ,     
         'DHW'= dhw    ,
         "HDI" = HDI,
         'assemblage size'= bw_size,
         'voice and acc'= voice ,
         #'coral cover'= bsp_nutrientscore_miCoralCover ,
         'Flat'= geomorphologyFlat ,
         'Lagoon'= geomorphologyLagoon_Backreef ,
         'Crest'=geomorphologyCrest,
         'wave energy' = wave_energy)


dag_estimates <- data.frame(median=apply(dag_output, 2, median))
dag_estimates$abs_effect <- abs(dag_estimates$median)
dag_estimates <- dag_estimates %>%
  arrange(desc(abs_effect))

# MATCH DAG ORDER TO CAUSAL SALAD ORDER?
dag_output <- dag_output[,order(match(colnames(dag_output), colnames(full_post)))]

#########################
## SIMPLE FOREST PLOTS ##
#########################

color_scheme_set("darkgray")
mcmc_areas_ridges(full_post) +
  ggtitle("Causal Salad Model") +
  xlim(range(full_post, dag_output))

mcmc_areas_ridges(dag_output) +
  ggtitle("DAG-Based Models") +
  xlim(range(full_post, dag_output))

#########################
## SIMPLE FOREST PLOTS ##
#########################

color_scheme_set("darkgray")
mcmc_intervals(full_post) +
  ggtitle("Causal Salad Model") +
  xlim(range(full_post, dag_output))

mcmc_intervals(dag_output) +
  ggtitle("DAG-Based Models") +
  xlim(range(full_post, dag_output))

ggpubr::ggarrange(mcmc_intervals(full_post) +
                    ggtitle("Causal Salad Model") +
                    xlim(range(full_post, dag_output)),
                  
                  mcmc_intervals(dag_output) +
                    ggtitle("DAG-Based Models") +
                    xlim(range(full_post, dag_output)),
                  
                  ncol = 2)


####################################
## RANK DAG ORDER BY EFFECT SIZE? ##
####################################

dag_output <- dag_output[,order(match(colnames(dag_output), rownames(dag_estimates)))]
full_post <- full_post[,order(match(colnames(full_post), colnames(dag_output)))]

#########################
## SIMPLE FOREST PLOTS ##
#########################

p <- ggpubr::ggarrange(
  
  mcmc_intervals(dag_output) +
    ggtitle("DAG-Based Models - Nutrient score calcium/iron/zinc") +
    xlim(range(full_post, dag_output)),
  
  mcmc_intervals(full_post) +
    ggtitle("Causal Salad Model") +
    xlim(range(full_post, dag_output)),
  
  ncol = 2)

jpeg("compare_model_nutscore_cafezn.jpeg", res=300, width=4000, height=2000)
p
graphics.off()






####################
## MODEL CHECKING ##
####################

# NOT YET DONE

