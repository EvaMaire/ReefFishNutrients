#compute relationships between micronutrient density and biodiversity

library(ggplot2)
library(tidyverse)
library(cowplot)
library(patchwork)
library(RColorBrewer)
library(here)

here()

#customised theme
white_themejpg <-theme(axis.ticks=element_line(colour="black"),
                       axis.text=element_text(size=14,colour="black"),
                       axis.title=element_text(size=16),
                       panel.grid.minor=element_blank(),
                       panel.background=element_rect(fill="white",colour="black"),
                       plot.background=element_rect(fill="transparent",colour=NA),
                       legend.key = element_rect(fill = "white"),
                       plot.title = element_text(hjust = 0.5))

#import data
load('data/reefdata.RData')
dat <- reefdata %>% filter(nutrient_score>0) #remove site with nutrient_score=0 (zero target fish biomass)

topfunctioE <- quantile(dat$Func_Ent,na.rm = T)[4]
topnut <- quantile(dat$nutrient_score,na.rm = T)[4]

#plot micronutrient density versus trait diversity
traitd <- ggplot(dat, aes(Func_Ent,nutrient_score)) +
  geom_vline(xintercept=topfunctioE,linetype="dashed",colour="dark grey",size=0.5)+
  geom_hline(yintercept=topnut,linetype="dashed",colour="dark grey",size=0.5)+
  geom_point(color="dark grey",size=2,alpha=.6)+
  stat_smooth(color="#2f3030",method = "lm",se=T,fullrange = F,linewidth=1) +
  scale_x_continuous(name = "Trait Diversity",limits=c(0,6))+
  scale_y_continuous(name = "Micronutrient density, %",limits=c(0,34))+
  white_themejpg

traitd

#plot micronutrient density versus number of threatened species
dat$nb_sp_threatened <- as.factor(dat$nb_sp_threatened) 
dat <- dat %>%  mutate(nb_sp_th = recode(nb_sp_threatened,
                                                   #'5' = '4',
                                                   '6' = '5',
                                                   '7' = '5'))
#Quantile95
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

thsp <- ggplot(dat, aes(y=nutrient_score, x=nb_sp_th)) +
  geom_hline(yintercept=topnut,linetype="dashed",colour="dark grey",size=0.5)+
  stat_summary(fun.data = quantiles_95,geom="boxplot",fill="dark grey")+
  geom_jitter(color="black", size=0.3, alpha=0.6) +
  annotate("text", y = 1, x = 1, size=7, label = "*", fontface = c("bold")) +
  annotate("segment", x = 1.6, xend = 4.4, y = 3, yend = 3)+
  annotate("text", y = 1, x = 3, size=7, label = "*", fontface = c("bold")) +
  annotate("segment", x = 4.6, xend = 6.4, y = 3, yend = 3)+
  annotate("text", y = 1, x = 5.5, size=7, label = "*", fontface = c("bold")) +
  scale_y_continuous("",limits=c(0,34))+
  scale_x_discrete("Number of threatened species",labels=c("0","1","2","3","4",expression("">=5)))+
  white_themejpg + theme(legend.position="none")

thsp  

#panels A and B
figure4 <- ( traitd + thsp ) + plot_annotation(tag_levels = 'A')

jpeg("figures/Figure4.jpeg", width=4800, height=2600, res=300) 
figure4
graphics.off()

#end

