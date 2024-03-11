#libraries
library(here)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(ggrepel)
library(gridExtra)
library(ggExtra)
library(cowplot)

here()

load('data/reefdata.RData')
dat <- reefdata %>% filter(nutrient_score>0) #remove site with nutrient_score=0 (zero target fish biomass)

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

#Determine top and bottom quartiles
quantile(dat$nutrient_score,probs = seq(0,1,.1))
nsq1 <- quantile(dat$nutrient_score,probs = seq(0,1,.1))[2]
nsq2 <- quantile(dat$nutrient_score,probs = seq(0,1,.1))[10]

nutscoreR2 <- ggplot() +
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group),fill = "grey", color = "grey") +
  coord_fixed(xlim = c(30, 320),  ylim = c(-30, 30), ratio = 1.3)+
  #first highest values - top quartile
  geom_point(data = dat %>% filter(nutrient_score>=nsq2), aes(x=Site_Long2, y=Site_Lat, colour = nutrient_score), pch=16, size=12,alpha=1)+
  #intermediate values 
  geom_point(data = dat %>% filter(nutrient_score<nsq2 & nutrient_score>=nsq1), aes(x=Site_Long2, y=Site_Lat, colour = nutrient_score), pch=16, size= 7, alpha=.6)+
  #lowest values - bottom quartile
  geom_point(data = dat %>% filter(nutrient_score<nsq1), aes(x=Site_Long2, y=Site_Lat, colour = nutrient_score), pch=16, size= 3, alpha=.6)+
  scale_colour_gradientn(name = "Micronutrient density\nof a 100g portion, %\n", colours = rev(myPalette(100)),
                         limits=c(5,35),breaks=c(5,10,15,20,25,30,35),
                         labels=c(expression(paste("5 ", italic("(less nutritious)"))),
                                  "10","15","20","25","30",
                                  expression(paste("35 ", italic("(more nutritious)")))) ) +
  annotate("text", x = 70, y = -15, label = "Indian Ocean",size=5) +
  annotate("text", x = 150, y = 2, label = "Indo-Pacific",size=5) +
  annotate("text", x = 200, y = 12, label = "Central Pacific",size=5) +
  annotate("text", x = 315, y = 15, label = "Western\nAltantic",size=5) +
  #scale_size_area(breaks=quantile(dat$nutrient_score)[-5], max_size=4, name='')+
  scale_x_continuous("longitude",breaks=c(80,180,280),labels=c(-100,0,100))+
  scale_y_continuous("latitude",breaks=c(-20,-10,0,10,20))+
  thememap + guides(size = "none") + theme(legend.justification = "left",legend.key.height = unit(1.3,"cm"),legend.text.align = 0,
                                           legend.text=element_text(size=11),
                                           legend.title = element_text(size=12 ,face="bold")#,legend.margin=margin(l=-3))
  )

nutscoreR2 

#legend
level_order <- c("Indian Ocean","Indo-Pacific","Central Pacific","Western Atlantic")

dat$region <- factor(dat$region, levels=level_order)

reg.cols2 <-c('cCentral Pacific'='#ffa31a', 'bIndo-Pacific'='#800a33', 'dWestern Atlantic'='#120b96', 'aIndian Ocean'='#0f7280')

pleg <- ggplot(dat, aes(nutrient_score, fill = factor(region, levels=c("Indian Ocean","Indo-Pacific","Central Pacific","Western Atlantic")),
                        color= factor(region, levels=c("Indian Ocean","Indo-Pacific","Central Pacific","Western Atlantic")) )) +
  geom_vline(xintercept=mean(dat$nutrient_score),linetype="dashed",colour="dark grey",size=0.5)+
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

bmnut <- ggplot(dat ,aes(y=nutrient_score, x=logbm,fill=region) ) +
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

#regional boxplots with dots

IO <- ggplot(dat %>% filter(region=='Indian Ocean'), aes(x=reorder(Larger, nutrient_score, decreasing=F), y = nutrient_score, color = nutrient_score)) +
  coord_flip() +
  scale_y_continuous(limits = c(5, 35), expand = c(0.02, 0.02)) +
  #scale_color_manual(values=reg.cols)+
  scale_colour_gradientn(colours = rev(myPalette(100))) +
  labs(x = NULL, y = "") +
  geom_hline(aes(yintercept = world_avg), color = "gray70", size = 0.6, linetype = "dashed") +
  geom_boxplot(color="gray70") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(size = 1, alpha = 0.6, width = 0.2)+
  #geom_jitter(size = 1, alpha = 0.6, width = 0.2)+
  #geom_point(data = countrydat %>% filter(region=='Indian Ocean'), size = 2, aes(x= Larger, y =meanns),colour='#333333') +
  white_themejpg + theme(legend.position="none") +
  ggtitle('Indian Ocean')+ 
  theme(plot.title = element_text(size = 15, face = "bold"))

IO 

IP <- ggplot(dat %>% filter(region=='Indo-Pacific'), aes(x=reorder(Larger, nutrient_score, decreasing=F), y = nutrient_score, color = nutrient_score)) +
  coord_flip() +
  scale_y_continuous(limits = c(5, 35), expand = c(0.02, 0.02)) +
  #scale_color_manual(values=reg.cols)+
  scale_colour_gradientn(colours = rev(myPalette(100))) +
  labs(x = NULL, y = "") +
  geom_hline(aes(yintercept = world_avg), color = "gray70", size = 0.6,linetype = "dashed") +
  geom_boxplot(color="gray70") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(size = 1, alpha = 0.6, width = 0.2)+
  #geom_jitter(size = 1, alpha = 0.6, width = 0.2)+
  #geom_point(data = countrydat%>% filter(region=='Indo-Pacific'), size = 2, aes(x= Larger, y =meanns),colour='#333333') +
  white_themejpg+theme(legend.position="none")+
  ggtitle('Indo-Pacific')+ 
  theme(plot.title = element_text(size = 15, face = "bold"))

IP 

CP <- ggplot(dat %>% filter(region=='Central Pacific'), aes(x=reorder(Larger, nutrient_score, decreasing=F), y = nutrient_score, color = nutrient_score)) +
  coord_flip() +
  scale_y_continuous(limits = c(5, 35), expand = c(0.02, 0.02)) +
  #scale_color_manual(values=reg.cols)+
  scale_colour_gradientn(colours = rev(myPalette(100))) +
  labs(x = NULL, y = "") +
  geom_hline(aes(yintercept = world_avg), color = "gray70", size = 0.6,linetype = "dashed") + 
  geom_boxplot(color="gray70") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(size = 1, alpha = 0.6, width = 0.2)+
  #geom_jitter(size = 1, alpha = 0.6, width = 0.2)+
  #geom_point(data = countrydat %>% filter(region=='Central Pacific'), size = 2, aes(x= Larger, y =meanns),colour='#333333') +
  white_themejpg+theme(legend.position="none")+
  ggtitle('Central Pacific')+ 
  theme(plot.title = element_text(size = 15, face = "bold"))

CP 

WA <- ggplot(dat %>% filter(region=='Western Atlantic'), aes(x=reorder(Larger, nutrient_score, decreasing=F), y = nutrient_score, color = nutrient_score)) +
  coord_flip() +
  scale_y_continuous(limits = c(5, 35), expand = c(0.02, 0.02)) +
  #scale_color_manual(values=reg.cols)+
  scale_colour_gradientn(colours = rev(myPalette(100))) +
  labs(x = NULL, y = "") +
  geom_hline(aes(yintercept = world_avg), color = "gray70", size = 0.6,linetype = "dashed") +
  geom_boxplot(color="gray70") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(size = 1, alpha = 0.6, width = 0.2)+
  #geom_jitter(size = 1, alpha = 0.6, width = 0.2)+
  #geom_point(data = countrydat%>% filter(region=='Western Atlantic'), size = 2, aes(x= Larger, y =meanns),colour='#333333') +
  white_themejpg+theme(legend.position="none")+
  ggtitle('Western Atlantic')+ 
  theme(plot.title = element_text(size = 15, face = "bold"))

WA 

#reverse size scale
p2 <- grid.arrange(nutscoreR2,
                   IO, IP, CP, WA,
                   ncol=4, nrow = 2, 
                   layout_matrix = rbind(c(1,1,1,1), c(2,3,4,5)),
                   heights = c(1.4,1))


final <- as_ggplot(p2) + # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 15,
                  x = c(0, 0, 0.25 ,0.5, 0.75), y = c(1, 0.45, 0.45,0.45,0.45))  # Add labels

ff <- annotate_figure(final,
                      bottom= text_grob("Micronutrient density (%)",  size = 14))


#Try to embed a miniature distribution
globaldistri <- ggplot(dat, aes(nutrient_score)) +
  geom_vline(xintercept=world_avg,linetype="dashed",colour="dark grey",size=.5)+
  geom_density(alpha = .3,colour='dark grey',size=1)+
  scale_x_continuous(name="",breaks=seq(5,35,5))+
  scale_y_continuous("Density",breaks=NULL)+ 
  ggtitle("Global distribution")+
  white_themejpg+
  theme(legend.position="none",plot.title = element_text(size = 9, face = "bold"),
        axis.title=element_text(size=8),
        axis.text=element_text(size=8), 
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))

#globaldistri

#library(cowplot)

emb <- ggdraw(nutscoreR2) +
  draw_plot(
    {
      globaldistri
    },
    # The distance along a (0,1) x-axis to draw the left edge of the plot
    x = 0.585, 
    # The distance along a (0,1) y-axis to draw the bottom edge of the plot
    y = 0.15,
    # The width and height of the plot expressed as proportion of the entire ggdraw object
    width = 0.15, 
    height = 0.3)

p2 <- grid.arrange(emb,
                   IO, IP, CP, WA,
                   ncol=4, nrow = 2, 
                   layout_matrix = rbind(c(1,1,1,1), c(2,3,4,5)),
                   heights = c(1.4,1))


# Add labels to the arranged plots
final <- as_ggplot(p2) + # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 15,
                  x = c(0, 0, 0.25 ,0.5, 0.75), y = c(1, 0.45, 0.45,0.45,0.45))  # Add labels

ff <- annotate_figure(final,
                      bottom= text_grob("Micronutrient density (%)",  size = 14))

tiff("/Volumes/EM2T/Lancaster/Nutrition/SERF/nutrient_value/dag-based/Figure1_ReverseTestEmbedded.tiff", width=4800, height=2600, compression="lzw", res=300) 
ff
graphics.off()

#END
