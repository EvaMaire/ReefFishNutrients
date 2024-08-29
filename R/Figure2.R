#Figure 2 

library(ggplot2)
library(tidyverse)
library(cowplot)
library(patchwork)
library(RColorBrewer)
library(factoextra)

white_themejpg <-theme(axis.ticks=element_line(colour="black"),
                       axis.text=element_text(size=14,colour="black"),
                       axis.title=element_text(size=16),
                       panel.grid.minor=element_blank(),
                       panel.grid.major=element_blank(),
                       axis.line = element_line(color = 'black'),
                       panel.background=element_rect(fill="white",colour=NA),
                       #plot.background=element_rect(fill="transparent",colour=NA),
                       panel.border = element_blank(),
                       legend.key = element_rect(fill = "white"),
                       plot.title = element_text(hjust = 0.5))


#Species level data
load("data/speciesdata.RData")

ghistSP <-ggplot(speciesdata, aes(nutrient_score)) + 
  geom_histogram(alpha=0.7, binwidth=5, bins=13,center=2.5, fill = "grey", col='black') +
  #geom_vline(data = asm_adq, aes(xintercept = portion_adq), col=asm.cols[4],size = 1.2) +
  scale_y_continuous(expand= c(0,0), name = "Number of reef fish species",limits=c(0,250),breaks = c(0,50,100,150,200,250)) +
  scale_x_continuous(name= "Micronutrient density (%)",breaks = c(0,10,20,30,40,50,60,70))+
  #scale_x_continuous(expand= c(0,0))
  #breaks = seq(40, 160, by = 20),
  #expand= c(0.01,0.01),
  #sec.axis = dup_axis(name = "", 
  #)) +
  white_themejpg

ghistSP

#range values
load('data/reefdata.RData')

#remove site with zero target biomass where nutrient score = 0
dat <- reefdata %>% filter(nutrient_score>0)
dim(dat)

#create log biomass
dat$logbm <- log(dat$bm_kg_ha+1)

#create regular biomass bins
step <- seq(1.1,10.1,0.2)
length(step)

dat$stepbm <- rep(1.1,nrow(dat))

for(k in 1:nrow(dat)){
  
  for(s in 1:length(step)) {
    cut <- step[s]
    if( dat$logbm[k] > cut ) { dat$stepbm[k] <- step[s] }
  }
  
}# end of k

#summarise data
df <- dat %>% select(nutrient_score,stepbm) %>% group_by(stepbm) %>% dplyr::summarize(maxns = max(nutrient_score),minns = min(nutrient_score))
df$col <- rep("normal")
df$col[which(df$stepbm=="5.7")] <- "highlight"

pol.cols<-c('normal'='#a8a7a7','highlight'="#b0520b")

rangeplot <- ggplot()+
  geom_segment(data = df,aes(x = stepbm, xend = stepbm,y = minns, yend = maxns,colour=col), size = 5, alpha = 0.6) +
  #geom_segment(aes(x = Item, xend = Item, y = Median-1, yend = Median+1), size = 5, colour = "black") +
  geom_point(data = dat, aes(y=nutrient_score,x=logbm), alpha=.7,color="dark grey",size=1)+
  stat_smooth(data = dat, aes(y=nutrient_score,x=logbm),method = 'lm',colour='black',
              se=T,fullrange = F,linetype = "dashed",size=1,show.legend = F)+
  scale_colour_manual(values = pol.cols)+
  #scale_colour_manual(values = fish.cols)+
  scale_x_continuous("Log fish biomass")+
  scale_y_continuous(name="Micronutrient density (%)",limits=c(5,35),breaks=seq(5,35,5))+
  white_themejpg+theme(legend.position="none") + theme(plot.title = element_text(size = 15, face = "bold"))

rangeplot

#combine panels A & B
patch <- (rangeplot + ghistSP) + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 16,face = 'bold'))
patch

jpeg("figures/Figure2.jpeg", res=300, width=4800, height=2600)
patch
graphics.off()

#coefficients
#model
#mod <- glm(nutrient_score ~ logbm, data=dat)
#summary(mod)

#summary(mod)$coefficients[2,1]

#END
