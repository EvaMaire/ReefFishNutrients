#compute biodiversity levels

setwd("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based")
library(ggplot2)
library(tidyverse)
#library(mapdata)
library(plot3D)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(ggrepel)
library(gridExtra)
library(ggExtra)
library(cowplot)

load('dat.RData')
head(dat)

#remove 1 site which has no target biomass
dat <- dat[-which(dat$BWnut3A==0),]

load('/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/fd.RData')
#load('/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/pd.Rdata')

UniqueSite <- rownames(fd)
fd <- data.frame(fd,UniqueSite)

#UniqueSite <- rownames(pd)
#pd <- data.frame(pd,UniqueSite)

#biodiv <- merge(fd,pd,by="UniqueSite",all.x=T)

final <- merge(dat,fd,by="UniqueSite",all.x=T)

#plot(final$Phyl_Ric,final$BWnut3A)
#plot(final$Func_Ric,final$BWnut3A)

#plot(final$Phyl_Ent,final$BWnut3A)
#plot(final$Func_Ent,final$BWnut3A)

white_themejpg <-theme(axis.ticks=element_line(colour="black"),
                       axis.text=element_text(size=12,colour="black"),
                       axis.title=element_text(size=14),
                       panel.grid.minor=element_blank(),
                       panel.background=element_rect(fill="white",colour="black"),
                       plot.background=element_rect(fill="transparent",colour=NA),
                       legend.key = element_rect(fill = "white"),
                       plot.title = element_text(hjust = 0.5))

#topphylo <- quantile(final$Phyl_Ric,na.rm = T)[4]

#ggplot(final, aes(Phyl_Ric,BWnut3A)) +
  #geom_point(color="dark grey",size=2,pch=21,alpha=.6)+
  #geom_vline(xintercept=topphylo,linetype="dashed",colour="dark grey",size=0.5)+
  ##geom_hline(yintercept=Se_under5*0.333,linetype="dashed",colour="dark grey",size=0.5)+
  #stat_smooth(color="dark grey",method = "loess",se=T,fullrange = F,size=1,show.legend = F) +
  #scale_x_continuous(name = "Phylogenetic Richness",limits=c(0,25))+
  #scale_y_continuous(name = "Micronutrient density score of a 100g portion, %")+
  ##scale_fill_manual(name = "Protection",values=cols)+
  #ggtitle('Nutrient density score')+
  #white_themejpg

#topphyloE <- quantile(final$Phyl_Ent,na.rm = T)[4]
#ggplot(final, aes(Phyl_Ent,BWnut3A)) +
  #geom_point(color="dark grey",size=2,pch=21,alpha=.6)+
  #geom_vline(xintercept=topphyloE,linetype="dashed",colour="dark grey",size=0.5)+
  ##geom_hline(yintercept=Se_under5*0.333,linetype="dashed",colour="dark grey",size=0.5)+
  #stat_smooth(color="dark grey",method = "lm",se=T,fullrange = F,size=1,show.legend = F) +
  #scale_x_continuous(name = "Phylogenetic Entropy",limits=c(0,7))+
  #scale_y_continuous(name = "Micronutrient density score of a 100g portion, %")+
  ##scale_fill_manual(name = "Protection",values=cols)+
  #ggtitle('Nutrient density score')+
  #white_themejpg


final$newnut <- final$BWnut3A/3
topfunctioE <- quantile(final$Func_Ent,na.rm = T)[4]
topnut <- quantile(final$newnut,na.rm = T)[4]
pp <- ggplot(final, aes(Func_Ent,newnut)) +
  geom_point(color="dark grey",size=2,alpha=.6)+
  geom_vline(xintercept=topfunctioE,linetype="dashed",colour="dark grey",size=0.5)+
  geom_hline(yintercept=topnut,linetype="dashed",colour="dark grey",size=0.5)+
  stat_smooth(color="#2f3030",method = "lm",se=T,fullrange = F,size=1,show.legend = F) +
  scale_x_continuous(name = "Trait Diversity",limits=c(0,6))+
  scale_y_continuous(name = "Micronutrient density score of a 100g portion, %",limits=c(0,34))+
  #scale_fill_manual(name = "Protection",values=cols)+
  white_themejpg

pp

#per country
cdata <- final %>% dplyr::select(UniqueSite,Larger,newnut,Func_Ent) %>% filter(Func_Ent>0)

nations <- cdata %>% dplyr::select(UniqueSite,Larger) %>% group_by(Larger) %>% count(Larger)

goodnations <- nations %>% filter(n>10)

finalc <- cdata %>% filter(Larger %in% goodnations$Larger)  

country <- ggplot(finalc, aes(Func_Ent,newnut,fill=Larger,colour=Larger)) +
  geom_point(color="dark grey",size=2,alpha=.6)+
  geom_vline(xintercept=topfunctioE,linetype="dashed",colour="dark grey",size=0.5)+
  geom_hline(yintercept=topnut,linetype="dashed",colour="dark grey",size=0.5)+
  #stat_smooth(color="#2f3030",method = "lm",se=T,fullrange = F,size=1,show.legend = F) +
  scale_x_continuous(name = "Trait Diversity",limits=c(0,6))+
  scale_y_continuous(name = "Micronutrient density score of a 100g portion, %",limits=c(0,34))+
  scale_fill_discrete(name = "Nations")+
  stat_smooth(method = "lm",se=T,fullrange = F,size=1,show.legend = F) +
  white_themejpg

country

#Models
mod <- lm(newnut ~ Func_Ent, data=finalc)
summary(mod)

summary(mod)$coefficients[2,1]

#Mixed model
library(lme4)
m1 <- lmer(newnut ~ Func_Ent + (Func_Ent | Larger), finalc)
coef <- coef(m1)$Larger
coef %>% arrange(Func_Ent)
summary(m1)

# random intercepts model
lmod <- lm(newnut ~ Func_Ent , data = finalc)
            
fm0 <- lmer(newnut ~ Func_Ent + (1 | Larger), finalc,REML = TRUE)

# random intercepts and random slopes model    
fm1 <- lmer(newnut ~ Func_Ent + (Func_Ent | Larger), finalc,REML = TRUE)
coef %>% arrange(Func_Ent)
# likelihood ratio test between the two models
anova(fm0, fm1, method = "LRT")
AIC(fm0,fm1) # Keep fm1

#preds <- predictInterval(fm1, newdata = finalc, n.sims = 999)

#allpred <- cbind(finalc,preds)
allpred <- finalc
allpred$fit <- predict(fm1) 
allpred$lmfit <- predict(lmod)

high <- c("French Polynesia","Venezuela","Australia","Seychelles","Commonwealth of the Northern Mariana Islands")
low <- c("Belize","Marshall Islands","Netherlands Antilles","Solomon Islands","Reunion")

#finalc$col <- rep("dark grey",nrow(finalc))
#finalc$col[which((finalc$Larger %in% high)==T)] <- "blue"
#finalc$col[which((finalc$Larger %in% low)==T)] <- "orange"

allpred$type <- rep("average",nrow(allpred))
allpred$type[which((allpred$Larger %in% high)==T)] <- "high"
allpred$type[which((allpred$Larger %in% low)==T)] <- "low"

country.cols<-c('average'= "dark grey", 
             'low'="orange", 
             'high'= "dark blue")

tot <- c(high,low)

allpred$newcountry <- allpred$Larger
allpred$newcountry[which((allpred$Larger %in%tot)==F)] <-'other'

country.cols2 <-c("French Polynesia"="#08519c",
"Venezuela"= "#3182bd",
"Australia"="#6baed6",
"Seychelles"="#9ecae1",
"Commonwealth of the Northern Mariana Islands"="#c6dbef",
"Belize"="#993404",
"Marshall Islands"="#d95f0e",
"Netherlands Antilles"="#fe9929",
"Solomon Islands"="#fec44f",
"Reunion"="#fee391")

#add position for labels
pos <- allpred %>% dplyr::select(Larger,fit) %>% dplyr::group_by(Larger) %>% 
                 mutate(label_ypos=max(fit))

pos <- unique(pos[,c("Larger","label_ypos")])
pos$newcountry <- pos$Larger

plyr::revalue(pos$Larger, c("Netherlands Antilles" ="ANT",
                            "Commonwealth of the Northern Mariana Islands" = "MNP",
                            #"Comoro Islands"="Comoros",
                            #'Papua New Guinea' = 'PNG',
                            "Marshall Islands"="Marshall Isl."
                            #"Solomon Islands"="Solomon Isl.",
                            )) -> pos$Larger

xpos <- c(1.18,1.35,1.3,0.7,1.7,1.2,1.25,1.2,1.7,1.3)
ypos <- c(17,22.5,15.7,25.5,14.5,27,25.5,18,17.2,23.5)

p <- ggplot(allpred, aes(Func_Ent,newnut,fill=newcountry,colour=newcountry)) +
  geom_point(size=1,alpha=.2)+
  geom_vline(xintercept=topfunctioE,linetype="dashed",colour="dark grey",size=0.5)+
  geom_hline(yintercept=topnut,linetype="dashed",colour="dark grey",size=0.5)+
  #stat_smooth(color="#2f3030",method = "lm",se=T,fullrange = F,size=1,show.legend = F) +
  geom_line(aes(y=lmfit), color="#2f3030",size=1) +
  scale_x_continuous(name = "Trait Diversity",limits=c(0,6))+
  scale_y_continuous(name = "Micronutrient density score of a 100g portion, %",limits=c(0,34))+
  scale_fill_manual(name = "Nations",values=country.cols2)+
  scale_colour_manual(name = "Nations",values=country.cols2)+
  geom_line(data = allpred %>% filter(newcountry != 'other'), aes(y=fit,fill=newcountry,colour=newcountry), size=2) +
  geom_text(data = pos %>% filter(newcountry %in% tot),
            aes(y = ypos,x=xpos, label= Larger, color=newcountry), size=5,fontface = "bold")+
  white_themejpg + theme(legend.position = 'none') 

tiff("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/SI_figures/Fig3_Countries.tiff", width=2400, height=2400, compression="lzw", res=300) 
p
graphics.off()


#add threatened species
load('/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dat_iucn_withsharks.RData')
funct <- merge(final,dat_iucn_withsharks,by="UniqueSite")

load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/datMEOW.RData")
dd <- unique(datMEOW[,c("Larger","REALM")])
dd <- dd[-2,]

funct <- merge(funct,dd,by="Larger",all.x=T,all.y=F) 

datfinal <- funct %>% select(UniqueSite:Larger,Protection,REALM,targetgroup:BWnut6,
                             herbdetri:planktivore,logbm,Func_Ent,nb_sp_threatened)%>%
  filter(targetgroup =='target') 

datfinal$newnut <- datfinal$BWnut3A/3

#Quantile95
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

datfinal$nb_sp_threatened <- as.factor(datfinal$nb_sp_threatened) 
datfinal <- datfinal %>%  mutate(nb_sp_th = recode(nb_sp_threatened,
                                                   #'5' = '4',
                                                   '6' = '5',
                                                   '7' = '5'))

pbox1 <- ggplot(datfinal, aes(y=BWnut3A, x=nb_sp_th)) +
  stat_summary(fun.data = quantiles_95, geom="boxplot",fill="dark grey")+
  #annotate("text", x = 0.6, y = 100, size=6, label = "A", fontface = c("bold")) +
  #annotate("text", x = 4, y = 75, size=7, label = "*", fontface = c("bold")) +
  #annotate("text", x = 1, y = 75, size=7, label = "*", fontface = c("bold")) +
  labs(x="")+
  scale_fill_brewer(palette="Dark2") +
  scale_y_continuous("Micronutrient density score of a 100g portion, %") +
  scale_x_discrete("Nb of threatened species - IUCN",labels=c("0","1","2","3","4",">5"))+
  white_themejpg + theme(legend.position="none")

pbox1

boxp <- ggplot(datfinal, aes(y=newnut, x=nb_sp_th)) +
  geom_hline(yintercept=topnut,linetype="dashed",colour="dark grey",size=0.5)+
  stat_summary(fun.data = quantiles_95,geom="boxplot",fill="dark grey")+
  geom_jitter(color="black", size=0.3, alpha=0.6) +
  #stat_summary(fun="mean", geom = "crossbar",  width = 0.5)+
  #scale_fill_viridis(discrete=T,name="Thermal regimes",labels=c("Cold","Temperate","Subtropical","Tropical")) +
  scale_y_continuous("",limits=c(0,34))+
  scale_x_discrete("Nb of threatened species - IUCN",labels=c("0","1","2","3","4",expression("">=5)))+
  white_themejpg + theme(legend.position="none")

boxp 

library(patchwork)
patchV2 <- ( pp + boxp ) + plot_annotation(tag_levels = 'A')

tiff("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/Figure3.tiff", width=4800, height=2600, compression="lzw", res=300) 
patchV2
graphics.off()

#Tukey's test
res.aov <- aov(BWnut3A ~ nb_sp_th, data = datfinal)
#res.aov <- aov(NutScore_under5 ~ EnvTemp, data = s)
#summary(res.aov)
#TukeyHSD(res.aov)

#kruskal.test(Climate_V.index ~ EnvTemp, data = s)

par(mar=c(4.1,10,4,4))
plot(TukeyHSD(res.aov, conf.level = 0.95),las=1, col = "red")


#model
mod <- glm(newnut ~ Func_Ent, data=datfinal)
summary(mod)

summary(mod)$coefficients[2,1]


#try national averages
datcountry <- final %>% select(Larger,Func_Ent,newnut) %>% group_by(Larger) %>%
  summarise(across(Func_Ent:newnut, list(mean = ~mean(.x),
                                          min = ~min(.x),
                                         max = ~max(.x))))

load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/datMEOW.RData")
dd <- unique(datMEOW[,c("Larger","REALM")])
dd <- dd[-2,]

dat <- merge(datcountry,dd,by="Larger",all.x=T,all.y=F) 
dat <- dat %>% mutate(region = recode(REALM, 
                                      "Eastern Indo-Pacific"="Central Pacific",
                                      "Central Indo-Pacific" = "Indo-Pacific",
                                      "Tropical Atlantic" = "Western Atlantic",
                                      "Western Indo-Pacific" = "Indian Ocean"))

plyr::revalue(dat$Larger, c("British Indian Ocean Territory" = "BIOT",
                            "Federated States of Micronesia" = "Micronesia",
                            "Netherlands Antilles" ="ANT",
                            "Commonwealth of the Northern Mariana Islands" = "MNP",
                            "Comoro Islands"="Comoros",
                            'Papua New Guinea' = 'PNG',
                            "Marshall Islands"="Marshall Isl.",
                            "Solomon Islands"="Solomon Isl.",
                            "Cayman Islands"="Cayman Isl.")) -> dat$Larger

# Classic ggplot
p <- ggplot(dat,aes(x=Func_Ent_mean, y=newnut_mean, text=Larger)) +
  geom_point(alpha=0.7) +
  geom_errorbar(aes(ymin = newnut_min,ymax = newnut_max),color='dark grey') + 
  geom_errorbarh(aes(xmin = Func_Ent_min ,xmax = Func_Ent_max),color='dark grey')+
  geom_vline(xintercept=topfunctioE,linetype="dashed",colour="dark grey",size=0.5)+
  geom_hline(yintercept=topnut,linetype="dashed",colour="dark grey",size=0.5)+
  scale_x_continuous(name = "Trait Diversity",limits=c(0,6))+
  scale_y_continuous(name = "Micronutrient density score, %",limits=c(0,34))+
  #scale_fill_manual(name = "Protection",values=cols)+
  geom_text_repel(data = dat , aes(label=Larger), size = 6,
                  direction = "both",
                  max.overlaps = 100)+
                  # Add extra padding around each text label.
                  #box.padding = unit(.55, 'cm'))+
  white_themejpg +
  facet_wrap(~REALM)

p

jpeg("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/SI_figures/Potential_Fig4.jpeg", res=300, width=4000, height=4000)
p             
graphics.off()

#end

