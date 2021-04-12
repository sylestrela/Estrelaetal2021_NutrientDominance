" R script for Figure 3D in Estrela et al (2021) 
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

##. 1) plot predicted vs observed at the FAMILY-LEVEL
#.. SELECT data to plot
dat.p<-dat
dat.p<-select(dat.p,-od_duo, -od1_mean_mono,-od2_mean_mono,-ab1_mean_mono,-ab2_mean_mono,-CS1_mono,-CS2_mono,-ab1_corr,-ab2_corr)

#. calculate mean of predicted abundance and SE
mean.pred<-unique(select(dat.p,-ab_duo,-rep_duo)) %>% 
  dplyr::group_by(Inoculum, Family, CS_duo, CS_type_duo,CS1,CS2) %>%
  dplyr::summarize(ab_pred_mean=mean(ab_pred),
                   n=n(), 
                   sd.pred=sd(ab_pred), 
                   ab_pred_se=sd.pred/sqrt(n))

d.mean.pred<-as.data.frame(select(mean.pred,-n,-sd.pred))
d.mean.pred[is.na(d.mean.pred)] <- 0 

#. calculate mean of observed abundance in pairs and SE
d.obs.mean<-unique(select(dat.p,-rep_mono1,-rep_mono2,-ab_pred)) %>% 
  dplyr::group_by(Inoculum, Family, CS_duo, CS_type_duo,CS1,CS2) %>%
  dplyr::summarize(ab_obs_mean=mean(ab_duo),
                   n=n(), 
                   sd.obs=sd(ab_duo), 
                   ab_obs_se=sd.obs/sqrt(n))

d.obs.mean<-select(d.obs.mean<-as.data.frame(d.obs.mean), -n,-sd.obs)
d.obs.mean[is.na(d.obs.mean)] <- 0 

##. merge 2 dataframes
dat.pred.obs<-merge(d.obs.mean,d.mean.pred)

##  ---------------------------------------------------------------------------
#. statistics- RMSE and Pearson's
library(Metrics)
dat<-dat.pred.obs

cor.test(dat$ab_obs_mean,dat$ab_pred_mean, method='pearson')
rmse(dat$ab_obs_mean, dat$ab_pred_mean)

dat.ss<-subset(dat, CS_type_duo=='SS')
fit.l.ss <- lm(ab_obs_mean~ab_pred_mean, dat.ss)
cor(dat.ss$ab_obs_mean,dat.ss$ab_pred_mean)
rmse(dat.ss$ab_obs_mean, dat.ss$ab_pred_mean)

dat.aa<-subset(dat, CS_type_duo=='AA')
fit.l.aa <- lm(ab_obs_mean~ab_pred_mean, dat.aa)
cor(dat.aa$ab_obs_mean,dat.aa$ab_pred_mean)
rmse(dat.aa$ab_obs_mean, dat.aa$ab_pred_mean)

dat.sa<-subset(dat, CS_type_duo=='SA')
fit.l.sa <- lm(ab_obs_mean~ab_pred_mean, dat.sa)
cor(dat.sa$ab_obs_mean,dat.sa$ab_pred_mean)
rmse(dat.sa$ab_obs_mean, dat.sa$ab_pred_mean)

#.. Plot Figure 3D
col.aa<-'tan1'
col.sa<-'gray'
col.ss<-'mediumpurple1'

dp<-dat.pred.obs
p1<-ggplot(dp, aes(x=ab_pred_mean, y=ab_obs_mean, group=CS_type_duo, colour=CS_type_duo,fill=CS_type_duo))+
  geom_point(alpha=0.8,size=2.5, shape=21) + 
  geom_errorbar(aes(ymin=ab_obs_mean-ab_obs_se, ymax=ab_obs_mean+ab_obs_se))+
  geom_errorbarh(aes(xmin=ab_pred_mean-ab_pred_se, xmax=ab_pred_mean+ab_pred_se))+
  theme_minimal() + 
  labs(x='predicted abundance' , y='observed abundance')+
  #     caption= paste0('Error bars represent mean +/- SE. For observed abundance, n =<4. For predicted abundance, n<=16.'))+
  geom_abline(linetype='dashed')+
  scale_colour_manual(values = c(col.aa, col.sa, col.ss))+
  scale_fill_manual(values = c( col.aa, col.sa, col.ss))+
  guides(fill=guide_legend(title="CS type"))+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14),
        axis.title.y =element_text(size=16),
        axis.title.x = element_text(size = 16),
        panel.border = element_rect(color = "gray", fill = NA, size = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))+
  scale_x_continuous( breaks=seq(0,1,1/2)) +
  scale_y_continuous( breaks=seq(0,1,1/2))+
  guides( colour = FALSE)
p1  
ggsave(path=path.p, "fig3D.pdf", width=5.2, height=4, onefile = TRUE)

#.. write data figure as csv file
ds<-dp
fname<-paste("fig3D_data.csv")
write.csv(ds, file.path(path.ds, fname), row.names=FALSE)

