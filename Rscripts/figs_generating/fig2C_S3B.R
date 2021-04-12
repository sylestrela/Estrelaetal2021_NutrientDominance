" R script for Figure 2C and Figure 2- figure supplement 3B in Estrela et al (2021) 
Nutrient dominance governs the assembly of microbial communities in mixed nutrient environments."

#.. Clear workspace
rm(list=ls())

#.. Load packages
library(data.table)
library(ggplot2)
library(dplyr)
library(plyr)
library(ggpubr)
library(reshape2)
library(tidyr)

#.. Assign path to data
#. usern = path to user directory
path.d<-paste0(usern, "/Estrelaetal_2021/data")
#.. Assign path to plots output
path.p<- paste0(usern,"/Estrelaetal_2021/Rplots")
#.. Assign path to data figures
path.ds<- paste0(usern,"/Estrelaetal_2021/data/data_figures/")

#.. Read data 
setwd(path.d)
dat<-read.csv('data_processed/data_mixedCS_obs_pred_fam_duos.csv')

#.. To plot predicted vs observed at the FAMILY-LEVEL
#. select data to process
dat.p<-select(dat,-od_duo, -od1_mean_mono,-od2_mean_mono,-ab1_mean_mono,-ab2_mean_mono,-CS1_mono,-CS2_mono,-ab1_corr,-ab2_corr)

#. calculate mean of predicted abundance and SE
mean.pred<-unique(select(dat.p,-ab_duo,-rep_duo)) %>% 
  dplyr::group_by(Inoculum, Family, CS_duo, CS_type_duo,CS1,CS2) %>%
  dplyr::summarize(ab_pred_mean=mean(ab_pred),
                   n=n(), 
                   sd_pred=sd(ab_pred), 
                   ab_pred_se=sd_pred/sqrt(n))

d.mean.pred<-as.data.frame(select(mean.pred,-n,-sd_pred))

#. calculate mean of observed abundance in pairs and SE
d.obs.mean = unique(select(dat.p,-rep_mono1,-rep_mono2,-ab_pred)) %>% 
  dplyr::group_by(Inoculum, Family, CS_duo, CS_type_duo,CS1,CS2) %>%
  dplyr::summarize(ab_obs_mean=mean(ab_duo),
                   n=n(), 
                   sd_obs=sd(ab_duo), 
                   ab_obs_se=sd_obs/sqrt(n))

d.obs.mean<-select(as.data.frame(d.obs.mean), -n,-sd_obs)
d.obs.mean[is.na(d.obs.mean)] <- 0 

#.. merge the 2 dataframes
dat.pred.obs<-merge(d.obs.mean,d.mean.pred)

#.. Colour scheme for plots
#. Family colors 
colEntbf<-'#516091'
colMoraxf<-'#74BEC1'
colPseuf<-'#E2F0CB'
colRhizf<-'#996666'
colOther<-'gray85'

#.. Plot Figure 2C
dp<-dat.pred.obs
#. select dominant families to plot
dps<-subset(dp, Family %in% c('Enterobacteriaceae', 'Moraxellaceae','Pseudomonadaceae', 'Rhizobiaceae'))
levels(factor(dps$CS2))
dps$cs.ord <- factor(dps$CS2, levels = c("fructose", "cellobiose", "ribose", "glycerol",
                                         "benzoate","fumarate","glycine","glutamine"))

p1<-ggplot(dps, aes(x=ab_pred_mean, y=ab_obs_mean, group=Family, colour=Family, shape=factor(Inoculum)))+
  geom_point(alpha=0.9,size=2.5) +
  geom_errorbar(aes(ymin=ab_obs_mean-ab_obs_se, ymax=ab_obs_mean+ab_obs_se))+
  geom_errorbarh(aes(xmin=ab_pred_mean-ab_pred_se, xmax=ab_pred_mean+ab_pred_se))+
  theme_minimal() +
  labs(x='predicted abundance', y='observed abundance')+
  geom_abline(linetype='dashed')+
  scale_colour_manual(values = c(colEntbf, colMoraxf,colPseuf,colRhizf))+
  guides(colour=guide_legend(title="Family"), shape=guide_legend(title='Inoculum'))+
  facet_grid(CS1~cs.ord)+
  theme(axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14),
        axis.title.y =element_text(size=14),
        axis.title.x = element_text(size = 14),
        panel.border = element_rect(color = "gray", fill = NA, size = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  scale_y_continuous(breaks=seq(0,1,1), labels = function(x) ifelse(x == 0, "0", x))+
  scale_x_continuous(breaks=seq(0,1,1), labels = function(x) ifelse(x == 0, "0", x))
p1
ggsave(path=path.p, "fig2C.pdf", width=14.5, height=3.7, onefile = TRUE)

#.. write data figure as csv file
ds<-dps
fname<-paste("fig2C_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)

#.. Plot Figure 2- figure supplement 3B (figure 2C in log-log scale)
p2<-ggplot(dps, aes(x=ab_pred_mean, y=ab_obs_mean, group=Family, colour=Family, shape=factor(Inoculum)))+
  geom_point(alpha=0.9,size=2) +
  geom_errorbar(aes(ymin=ab_obs_mean-ab_obs_se, ymax=ab_obs_mean+ab_obs_se))+
  geom_errorbarh(aes(xmin=ab_pred_mean-ab_pred_se, xmax=ab_pred_mean+ab_pred_se))+
  theme_minimal() +
  labs(x='predicted abundance', y='observed abundance')+
  geom_abline(linetype='dashed')+
  scale_colour_manual(values = c(colEntbf, colMoraxf,colPseuf,colRhizf))+
  guides(colour=guide_legend(title="Family"), shape=guide_legend(title='Inoculum'))+
  facet_grid(CS1~cs.ord)+
  theme(axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14),
        axis.title.y =element_text(size=14),
        axis.title.x = element_text(size = 14),
        panel.border = element_rect(color = "gray", fill = NA, size = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  scale_x_log10() +
  scale_y_log10()
p2
ggsave(path=path.p, "fig2S3B.pdf", width=21.75, height=5.7, onefile = TRUE)


