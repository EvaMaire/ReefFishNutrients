#libraries
library(here)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(ggrepel)
library(gridExtra)
library(ggExtra)
library(cowplot)

here()

load('dat.RData')
head(dat)
dat <- dat[-which(dat$BWnut3A==0),]

#customized theme 
thememap<-theme(axis.text=element_text(colour="black"),
                axis.ticks=element_line(colour="black"),
                panel.grid.minor=element_blank(),
                panel.grid.major=element_blank(),
                panel.background=element_rect(fill="white",colour="black"),
                plot.background=element_rect(fill="transparent",colour=NA),
                plot.margin = margin(-10,1,-10,1, "cm"))

map <- map_data("world")
dat$Site_Long2 <- ifelse(dat$Site_Long < -25, dat$Site_Long + 360, dat$Site_Long) # where d is your df
mapWorld <- map_data('world', wrap=c(-25,335), ylim=c(-55,75))

#myPalette <- colorRampPalette(brewer.pal(4, "Spectral"))
myPalette <- colorRampPalette(rev(brewer.pal(4, "BrBG")))

#Portion size
#dat$newnut <- dat$BWnut3A/3
quantile(dat$BWnut3A,probs = seq(0,1,.1))
q1 <- quantile(dat$BWnut3A,probs = seq(0,1,.1))[2]
q2 <- quantile(dat$BWnut3A,probs = seq(0,1,.1))[10]

nutscore <- ggplot() +
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group),fill = "grey", color = "grey") +
  coord_fixed(xlim = c(30, 320),  ylim = c(-30, 30), ratio = 1.3)+
  geom_point(data = dat %>% filter(BWnut3A<q1), aes(x=Site_Long2, y=Site_Lat, colour = BWnut3A), pch=16, size= 10, alpha=.6)+
  geom_point(data = dat %>% filter(BWnut3A<q2 & BWnut3A>=q1), aes(x=Site_Long2, y=Site_Lat, colour = BWnut3A), pch=16, size= 5, alpha=.6)+
  geom_point(data = dat %>% filter(BWnut3A>=q2), aes(x=Site_Long2, y=Site_Lat, colour = BWnut3A), pch=16, size=2,alpha=1)+
  scale_colour_gradientn(name = "Nutrient density score\nof a 100g portion, %\n", colours = rev(myPalette(100)),
                         limits=c(15,102),breaks=c(15,40,70,102),
                         labels=c(expression(paste("15 ", italic("(less nutritious)"))),
                                  "40","70",
                                  expression(paste("102 ", italic("(more nutritious)")))) ) +
  annotate("text", x = 70, y = -15, label = "Indian Ocean",size=5) +
  annotate("text", x = 150, y = 2, label = "Indo-Pacific",size=5) +
  annotate("text", x = 200, y = 12, label = "Central Pacific",size=5) +
  annotate("text", x = 315, y = 15, label = "Western\nAltantic",size=5) +
  #scale_size_area(breaks=quantile(dat$BWnut3A)[-5], max_size=4, name='')+
  scale_x_continuous("longitude",breaks=c(80,180,280),labels=c(-100,0,100))+
  scale_y_continuous("latitude",breaks=c(-20,-10,0,10,20))+
  thememap + guides(size = "none") + theme(legend.justification = "left",legend.key.height = unit(1,"cm"),legend.text.align = 0,
                                           legend.text=element_text(size=11),
                                           legend.title = element_text(size=12 ,face="bold")#,legend.margin=margin(l=-3))
  )

nutscore 

dat$ns100 <- dat$BWnut3A/3
nsq1 <- quantile(dat$ns100,probs = seq(0,1,.1))[2]
nsq2 <- quantile(dat$ns100,probs = seq(0,1,.1))[10]

nutscore2 <- ggplot() +
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group),fill = "grey", color = "grey") +
  coord_fixed(xlim = c(30, 320),  ylim = c(-30, 30), ratio = 1.3)+
  geom_point(data = dat %>% filter(ns100<nsq1), aes(x=Site_Long2, y=Site_Lat, colour = ns100), pch=16, size= 12, alpha=.6)+
  geom_point(data = dat %>% filter(ns100<nsq2 & ns100>=nsq1), aes(x=Site_Long2, y=Site_Lat, colour = ns100), pch=16, size= 7, alpha=.6)+
  geom_point(data = dat %>% filter(ns100>=nsq2), aes(x=Site_Long2, y=Site_Lat, colour = ns100), pch=16, size=4,alpha=1)+
  scale_colour_gradientn(name = "Micronutrient density\nscore of a 100g portion, %\n", colours = rev(myPalette(100)),
                         limits=c(5,35),breaks=c(5,10,15,20,25,30,35),
                         labels=c(expression(paste("5 ", italic("(less nutritious)"))),
                                  "10","15","20","25","30",
                                  expression(paste("35 ", italic("(more nutritious)")))) ) +
  annotate("text", x = 70, y = -15, label = "Indian Ocean",size=5) +
  annotate("text", x = 150, y = 2, label = "Indo-Pacific",size=5) +
  annotate("text", x = 200, y = 12, label = "Central Pacific",size=5) +
  annotate("text", x = 315, y = 15, label = "Western\nAltantic",size=5) +
  #scale_size_area(breaks=quantile(dat$BWnut3A)[-5], max_size=4, name='')+
  scale_x_continuous("longitude",breaks=c(80,180,280),labels=c(-100,0,100))+
  scale_y_continuous("latitude",breaks=c(-20,-10,0,10,20))+
  thememap + guides(size = "none") + theme(legend.justification = "left",legend.key.height = unit(1,"cm"),legend.text.align = 0,
                                           legend.text=element_text(size=11),
                                           legend.title = element_text(size=12 ,face="bold")#,legend.margin=margin(l=-3))
  )

nutscore2 

#histograms biomass and nutrient density score
load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/datMEOW.RData")
dd <- unique(datMEOW[,c("Larger","REALM")])
dd <- dd[-2,]
dat <- merge(dat,dd,by="Larger",all.x=T,all.y=F) 
dat <- dat %>% mutate(region = recode(REALM, 
                                       "Eastern Indo-Pacific"="Central Pacific",
                                       "Central Indo-Pacific" = "Indo-Pacific",
                                      "Tropical Atlantic" = "Western Atlantic",
                                      "Western Indo-Pacific" = "Indian Ocean"))

white_themejpg <-theme(axis.ticks=element_line(colour="black"),
                       axis.text=element_text(size=12,colour="black"),
                       axis.title=element_text(size=14),
                       panel.grid.minor=element_blank(),
                       panel.grid.major=element_blank(),
                       axis.line = element_line(color = 'black'),
                       panel.background=element_rect(fill="white",colour=NA),
                       #plot.background=element_rect(fill="transparent",colour=NA),
                       panel.border = element_blank(),
                       legend.key = element_rect(fill = "white"),
                       plot.title = element_text(hjust = 0.5))

reg.cols<-c('Central Pacific'='#ffa31a', 'Indo-Pacific'='#800a33', 'Western Atlantic'='#120b96', 'Indian Ocean'='#0f7280')

bm <- ggplot(dat, aes(logbm, fill = region, colour = region)) +
  geom_vline(xintercept=mean(dat$logbm),linetype="dashed",colour="dark grey",size=1)+
  geom_density(alpha = 0.3)+
  scale_fill_manual(values=reg.cols)+
  scale_color_manual(values=reg.cols)+
  scale_x_continuous("Log standing biomass")+
  scale_y_continuous("Density",breaks=NULL)+
  white_themejpg+theme(legend.position="none")

bm

ns <- ggplot(dat, aes(BWnut3A, fill = region, colour = region)) +
  geom_vline(xintercept=mean(dat$BWnut3A),linetype="dashed",colour="dark grey",size=1)+
  geom_density(alpha = .3)+
  scale_fill_manual(values=reg.cols)+
  scale_color_manual(values=reg.cols)+
  scale_x_continuous("Nutrient density score (%)",breaks=seq(20,110,10))+
  scale_y_continuous("Density",breaks=NULL)+
  white_themejpg+theme(legend.position="none")

ns

ns2 <- ggplot(dat, aes(ns100, fill = region, colour = region)) +
  geom_vline(xintercept=mean(dat$ns100),linetype="dashed",colour="dark grey",size=1)+
  geom_density(alpha = .3)+
  scale_fill_manual(values=reg.cols)+
  scale_color_manual(values=reg.cols)+
  scale_x_continuous("Micronutrient density score (%)",breaks=seq(5,35,5))+
  scale_y_continuous("",breaks=NULL)+
  white_themejpg+theme(legend.position="none")

ns2

#patch <-  (nutscore / (bm + ns)) + plot_layout(heights = c(1.5, 1))  + plot_annotation(tag_levels = 'a')

#legend
dat <- dat %>% mutate(regionOrder = recode(region,
                      "Indian Ocean"="aIndian Ocean",
                      "Indo-Pacific"= "bIndo-Pacific",
                      "Central Pacific"="cCentral Pacific",
                      "Western Atlantic"="dWestern Atlantic"))

level_order <- c("Indian Ocean","Indo-Pacific","Central Pacific","Western Atlantic")

dat$region <- factor(dat$region, levels=level_order)

reg.cols2 <-c('cCentral Pacific'='#ffa31a', 'bIndo-Pacific'='#800a33', 'dWestern Atlantic'='#120b96', 'aIndian Ocean'='#0f7280')

pleg <- ggplot(dat, aes(BWnut3A, fill = factor(region, levels=c("Indian Ocean","Indo-Pacific","Central Pacific","Western Atlantic")),
                        color= factor(region, levels=c("Indian Ocean","Indo-Pacific","Central Pacific","Western Atlantic")) )) +
  geom_vline(xintercept=mean(dat$BWnut3A),linetype="dashed",colour="dark grey",size=0.5)+
  geom_density(alpha = .3)+
  scale_fill_manual(values=reg.cols,name="",breaks=c("Indian Ocean","Indo-Pacific","Central Pacific","Western Atlantic"))+
  scale_color_manual(values=reg.cols,name="",breaks=c("Indian Ocean","Indo-Pacific","Central Pacific","Western Atlantic"))+
  scale_x_continuous("Nutrient density score (%)",breaks=seq(20,110,10))+
  scale_y_continuous("Density")+
  white_themejpg


#p <- pleg2 + guides(colour = guide_colourbar(order=1),
                    #size = guide_legend(order=2))
legend <- get_legend(
  pleg + 
    theme(legend.position = "top",legend.justification = "center",legend.direction = "horizontal",
          legend.text=element_text(size=14),legend.title = element_text(size=18 ,face="bold"),legend.key.width = unit(3,"line"))
)

#biplot
#Color scale
cols=c("Fished"="#d4d4d4","Restricted"="#595959","UnfishedHigh" = "#171717")

dat$alpha <- rep(0.9,nrow(dat))
dat$alpha[which(dat$Protection=="Fished")] <- .7

bmnut <- ggplot(dat ,aes(y=BWnut3A, x=logbm,fill=region) ) +
  geom_point(alpha=dat$alpha,shape = 21,color="dark grey",size=2)+
  stat_smooth(aes(fill=region,color=region),method = 'lm', #formula = y ~ splines::bs(x, 3),#"gam",formula = y ~ s(x, k = 3),
              se=T,fullrange = F,linetype = "dashed",size=1,show.legend = F) +
  #geom_segment(aes(x = 50, y = 0, xend = 50 , yend =100),linetype="dashed",size=0.5,color="darkgrey")+
  #geom_segment(aes(x = 0, y = 50, xend = 100 , yend =50),linetype="dashed",size=0.5,color="darkgrey")+
  scale_y_continuous(name = "Average nutrient density score, %",limits=c(15,102)) +
  scale_x_continuous(name = "Log fish biomass") +
  scale_fill_manual(values=reg.cols)+
  scale_color_manual(values=reg.cols)+
  #scale_color_brewer(palette="Greys")+
  white_themejpg + theme(legend.position = "right",legend.justification = "left",legend.key.height = unit(1,"cm"),legend.text.align = 0,
                      legend.text=element_text(size=12),
                      legend.title = element_text(size=14 ,face="bold"))+#,legend.margin=margin(l=-3))+
  guides(fill = guide_legend(override.aes = list(size=4))) 

bmnut

#Export Figure 2
p <- grid.arrange(nutscore2,
             legend,
             bm,ns2,
             ncol=2, nrow = 3, 
             layout_matrix = rbind(c(1,1), c(2,2),c(3,4)),
             widths = c(3, 3), heights = c(1.5,0.1,1))


# Add labels to the arranged plots
library(ggpubr)

final <- as_ggplot(p) + # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0, 0.5), y = c(1, 0.5, 0.5)) # Add labels
final

#tiff("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/Figure1.tiff", width=4800, height=2600, compression="lzw", res=300) 
#final
#graphics.off()


#rank countries
#change names
#library(plyr)
plyr::revalue(dat$Larger, c("British Indian Ocean Territory" = "BIOT",
  "Federated States of Micronesia" = "Micronesia",
  "Netherlands Antilles" ="ANT",
  "Commonwealth of the Northern Mariana Islands" = "MNP",
  "Comoro Islands"="Comoros",
  'Papua New Guinea' = 'PNG',
  "Marshall Islands"="Marshall Isl.",
  "Solomon Islands"="Solomon Isl.",
  "Cayman Islands"="Cayman Isl.")) -> dat$Larger

world_avg <- dat %>%
  summarize(meanns = mean(ns100, na.rm = TRUE)) %>%
  pull(meanns)

countrydat <- dat %>% select(Larger,region,ns100) %>% group_by(region,Larger) %>% 
  dplyr::summarize(meanns = mean(ns100)) 

reg.cols<-c('Central Pacific'='#ffa31a', 'Indo-Pacific'='#800a33', 'Western Atlantic'='#120b96', 'Indian Ocean'='#0f7280')

IO <- ggplot(dat %>% filter(region=='Indian Ocean'), aes(x=reorder(Larger, ns100, decreasing=F), y = ns100, color = ns100)) +
  coord_flip() +
  scale_y_continuous(limits = c(5, 35), expand = c(0.02, 0.02)) +
  #scale_color_manual(values=reg.cols)+
  scale_colour_gradientn(colours = rev(myPalette(100))) +
  labs(x = NULL, y = "") +
  geom_hline(aes(yintercept = world_avg), color = "gray70", size = 0.6,linetype = "dashed") +
  geom_jitter(size = 1, alpha = 0.6, width = 0.2)+
  geom_point(data = countrydat%>% filter(region=='Indian Ocean'), size = 2, aes(x= Larger, y =meanns),colour='#333333') +
  white_themejpg + theme(legend.position="none") +
  ggtitle('Indian Ocean')+ 
  theme(plot.title = element_text(size = 15, face = "bold"))

IO 

IP <- ggplot(dat %>% filter(region=='Indo-Pacific'), aes(x=reorder(Larger, ns100, decreasing=F), y = ns100, color = ns100)) +
  coord_flip() +
  scale_y_continuous(limits = c(5, 35), expand = c(0.02, 0.02)) +
  #scale_color_manual(values=reg.cols)+
  scale_colour_gradientn(colours = rev(myPalette(100))) +
  labs(x = NULL, y = "") +
  geom_hline(aes(yintercept = world_avg), color = "gray70", size = 0.6,linetype = "dashed") +
  geom_jitter(size = 1, alpha = 0.6, width = 0.2)+
  geom_point(data = countrydat%>% filter(region=='Indo-Pacific'), size = 2, aes(x= Larger, y =meanns),colour='#333333') +
  white_themejpg+theme(legend.position="none")+
  ggtitle('Indo-Pacific')+ 
  theme(plot.title = element_text(size = 15, face = "bold"))

IP 

CP <- ggplot(dat %>% filter(region=='Central Pacific'), aes(x=reorder(Larger, ns100, decreasing=F), y = ns100, color = ns100)) +
  coord_flip() +
  scale_y_continuous(limits = c(5, 35), expand = c(0.02, 0.02)) +
  #scale_color_manual(values=reg.cols)+
  scale_colour_gradientn(colours = rev(myPalette(100))) +
  labs(x = NULL, y = "") +
  geom_hline(aes(yintercept = world_avg), color = "gray70", size = 0.6,linetype = "dashed") +  
  geom_jitter(size = 1, alpha = 0.6, width = 0.2)+
  geom_point(data = countrydat %>% filter(region=='Central Pacific'), size = 2, aes(x= Larger, y =meanns),colour='#333333') +
  white_themejpg+theme(legend.position="none")+
  ggtitle('Central Pacific')+ 
  theme(plot.title = element_text(size = 15, face = "bold"))

CP 

WA <- ggplot(dat %>% filter(region=='Western Atlantic'), aes(x=reorder(Larger, ns100, decreasing=F), y = ns100, color = ns100)) +
  coord_flip() +
  scale_y_continuous(limits = c(5, 35), expand = c(0.02, 0.02)) +
  #scale_color_manual(values=reg.cols)+
  scale_colour_gradientn(colours = rev(myPalette(100))) +
  labs(x = NULL, y = "") +
  geom_hline(aes(yintercept = world_avg), color = "gray70", size = 0.6,linetype = "dashed") +
  geom_jitter(size = 1, alpha = 0.6, width = 0.2)+
  geom_point(data = countrydat%>% filter(region=='Western Atlantic'), size = 2, aes(x= Larger, y =meanns),colour='#333333') +
  white_themejpg+theme(legend.position="none")+
  ggtitle('Western Atlantic')+ 
  theme(plot.title = element_text(size = 15, face = "bold"))

WA 

(IO + IP + CP + WA)

#Export Figure 2
p2 <- grid.arrange(nutscore2,
                  IO, IP, CP, WA,
                  ncol=4, nrow = 2, 
                  layout_matrix = rbind(c(1,1,1,1), c(2,3,4,5)),
                  heights = c(1.4,1))


# Add labels to the arranged plots
library(ggpubr)

final <- as_ggplot(p2) + # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 15,
                  x = c(0, 0, 0.25 ,0.5, 0.75), y = c(1, 0.45, 0.45,0.45,0.45))  # Add labels

ff <- annotate_figure(final,
                  bottom= text_grob("Micronutrient density score (%)",  size = 14))

ff

tiff("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/Figure1_V2.tiff", width=4800, height=2600, compression="lzw", res=300) 
ff
graphics.off()


#jpeg("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/Figure1.jpeg", res=300, width=4000, height=4000)
#patch 
#graphics.off()

######################################
# SUPPLEMENTAL FIGURES
######################################

#Import nutrient predictions
#Species level data
load('/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/speciesdata_serf.RData')
load("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/nut_splevel_Sept2022.RData")

nutscore <- unique( nut_splevel_Sept2022 %>% select(FullSpecies,nutscore_3nutA))

spdat <- merge(speciesdata_serf,nutscore,by="FullSpecies",all.x=T)
dim(speciesdata_serf)
dim(spdat) #833 yeah!!!!

spdat$ns100 <- spdat$nutscore_3nutA/3

ghistSP <-ggplot(spdat, aes(ns100)) + 
  geom_histogram(alpha=0.7, binwidth=5, fill = "grey", col='black') +
  #geom_vline(data = asm_adq, aes(xintercept = portion_adq), col=asm.cols[4],size = 1.2) +
  scale_y_continuous(expand= c(0,0), name = "# coral reef fish species",limits=c(0,250),breaks = c(0,50,100,150,200,250)) +
  scale_x_continuous(name= "Micronutrient density score of a 100g portion, %",breaks = c(0,10,20,30,40,50,60,70))+
  #scale_x_continuous(expand= c(0,0))
  #breaks = seq(40, 160, by = 20),
  #expand= c(0.01,0.01),
  #sec.axis = dup_axis(name = "", 
  #)) +
  white_themejpg

ghistSP

#range values

head(dat)
dat$log10bm <- log10(dat$bm_kg_ha_US)
summary(dat$log10bm)

#step <- seq(0.3,4.4,0.1)
step <- seq(1.1,10.1,0.2)
length(step)

dat$stepbm <- rep(1.1,nrow(dat))

for(k in 1:nrow(dat)){
  
  for(s in 1:length(step)) {
    cut <- step[s]
    if( dat$logbm[k] > cut ) { dat$stepbm[k] <- step[s] }
  }
  
}# end of k


plot(dat$stepbm,dat$ns100)

df <- dat %>% select(ns100,stepbm) %>% group_by(stepbm) %>% dplyr::summarize(maxns = max(ns100),minns = min(ns100))
df$col <- rep("normal")
df$col[which(df$stepbm=="5.7")] <- "highlight"
  
pol.cols<-c('normal'='#a8a7a7','highlight'="#91130a")

rangeplot <- ggplot()+
  geom_segment(data = df,aes(x = stepbm, xend = stepbm,y = minns, yend = maxns,colour=col), size = 5, alpha = 0.6) +
  #geom_segment(aes(x = Item, xend = Item, y = Median-1, yend = Median+1), size = 5, colour = "black") +
  geom_point(data = dat, aes(y=ns100,x=logbm), alpha=.7,color="dark grey",size=1)+
  stat_smooth(data = dat, aes(y=ns100,x=logbm),method = 'lm',colour='black',
              se=T,fullrange = F,linetype = "dashed",size=1,show.legend = F)+
  scale_colour_manual(values = pol.cols)+
  #scale_colour_manual(values = fish.cols)+
  scale_x_continuous("Log fish biomass")+
  ylab("Average micronutrient density score, %")+
  white_themejpg+theme(legend.position="none")

rangeplot

patch <- (rangeplot + ghistSP) + plot_annotation(tag_levels = 'A')
patch

tiff("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/SI_figures/Figure S1X.tiff", width=4000, height=2400, compression="lzw", res=300) 
patch
graphics.off()

#END
