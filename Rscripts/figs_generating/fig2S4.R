" R script for Figure 2- figure supplement 4 in Estrela et al (2021) 
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

dat.d.f<-read.csv('data_processed/data_mixedCS_obs_pred_fam_duos.csv')
dat.d.g<-read.csv('data_processed/data_mixedCS_obs_pred_gen_duos.csv')
dat.d.e<-read.csv('data_processed/data_mixedCS_obs_pred_esv_duos.csv')

##.. Predicted vs observed colored by focal CS for family, genus, and esv level

##  ---------------------------
#. 1)  FAMILY-LEVEL 
#.. select data to plot
dat.p<-dat.d.f
dat.p<-select(dat.p,-od_duo, -od1_mean_mono,-od2_mean_mono,-ab1_mean_mono,-ab2_mean_mono,-CS1_mono,-CS2_mono,-ab1_corr,-ab2_corr)

#. calculate mean of predicted and SE
mean.pred<-unique(select(dat.p,-ab_duo,-rep_duo)) %>% 
  dplyr::group_by(Inoculum, Family, CS_duo, CS_type_duo,CS1,CS2) %>%
  dplyr::summarize(ab_pred_mean=mean(ab_pred),
                   n=n(), sd.pred=sd(ab_pred), 
                   ab_pred_se=sd.pred/sqrt(n))

d.mean.pred<-as.data.frame(select(mean.pred,-n,-sd.pred))
d.mean.pred[is.na(d.mean.pred)] <- 0 

#. calculate mean of observed and SE
d.obs.mean = unique(select(dat.p,-rep_mono1,-rep_mono2,-ab_pred)) %>% 
  dplyr::group_by(Inoculum, Family, CS_duo, CS_type_duo,CS1,CS2) %>%
  dplyr::summarize(ab_obs_mean=mean(ab_duo),
                   n=n(), sd.obs=sd(ab_duo), 
                   ab_obs_se=sd.obs/sqrt(n))

d.obs.mean<-select(as.data.frame(d.obs.mean), -n,-sd.obs)
d.obs.mean[is.na(d.obs.mean)] <- 0 

##. merge 2 dataframes
dat.pred.obs<-merge(d.obs.mean,d.mean.pred)

##. Plot 
col.aa<-'tan1'
col.sa<-'gray'
col.ss<-'mediumpurple1'

dp<-dat.pred.obs
p.title<-'Family'
p.fam<-ggplot(dp, aes(x=ab_pred_mean, y=ab_obs_mean, group=CS1, colour=CS1,fill=CS1))+
  geom_point(alpha=0.5,size=2.5, shape=21) + 
  geom_errorbar(aes(ymin=ab_obs_mean-ab_obs_se, ymax=ab_obs_mean+ab_obs_se))+
  geom_errorbarh(aes(xmin=ab_pred_mean-ab_pred_se, xmax=ab_pred_mean+ab_pred_se))+
  theme_minimal() + 
  labs(x='predicted abundance' , y='observed abundance',title=p.title)+
  geom_abline(linetype='dashed')+
  scale_colour_manual(values = c( col.ss, col.aa))+
  scale_fill_manual(values = c( col.ss,col.aa))+
  guides(fill=guide_legend(title="Focal CS"))+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14),
        axis.title.y =element_text(size=16),
        axis.title.x = element_text(size = 16),
        plot.title = element_text(size = 18),
        panel.border = element_rect(color = "gray", fill = NA, size = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))+
  scale_x_continuous( breaks=seq(0,1,1), limits=c(0,1)) +
  scale_y_continuous( breaks=seq(0,1,1), limits=c(0,1))+
  guides( colour = FALSE)
p.fam

#.. write data figure as csv file
ds<-dp
fname<-paste("fig2S4_fam_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)

#. statistics- RMSE and Pearson's
library(Metrics)
fit.l <- lm(ab_obs_mean~ab_pred_mean, dp)
cor.test(dp$ab_obs_mean,dp$ab_pred_mean, method='pearson')
cor.fam<-cor(dp$ab_obs_mean,dp$ab_pred_mean)
rmse.fam<-rmse(dp$ab_obs_mean, dp$ab_pred_mean)

dat.glu<-subset(dp, CS1=='glucose')
dim(dat.glu)
fit.l.glu <- lm(ab_obs_mean~ab_pred_mean, dat.glu)
cor.test(dat.glu$ab_obs_mean,dat.glu$ab_pred_mean, method='pearson')
cor.fam.glu<-cor(dat.glu$ab_obs_mean,dat.glu$ab_pred_mean)
rmse.fam.glu<-rmse(dat.glu$ab_obs_mean, dat.glu$ab_pred_mean)

dat.succ<-subset(dp, CS1=='succinate')
dim(dat.succ)
fit.l.succ <- lm(ab_obs_mean~ab_pred_mean, dat.succ)
cor.test(dat.succ$ab_obs_mean,dat.succ$ab_pred_mean, method='pearson')
cor.fam.succ<-cor(dat.succ$ab_obs_mean,dat.succ$ab_pred_mean)
rmse.fam.succ<-rmse(dat.succ$ab_obs_mean, dat.succ$ab_pred_mean)

paste0('cor.fam:', round(cor.fam[1],3),'; rmse.fam:', round(rmse.fam[1],3))
paste0('cor.fam.glu:', round(cor.fam.glu[1],3),'; cor.fam.succ:', round(cor.fam.succ[1],3))
paste0('rmse.fam.glu:', round(rmse.fam.glu[1],3),'; rmse.fam.succ:', round(rmse.fam.succ[1],3))

##  ---------------------------
## 2) GENUS-LEVEL
#.. select data to plot
dat.p<-dat.d.g
dat.p<-select(dat.p,-od_duo, -od1_mean_mono,-od2_mean_mono,-ab1_mean_mono,-ab2_mean_mono,-CS1_mono,-CS2_mono,-ab1_corr,-ab2_corr)

#. calculate mean of predicted and SE
mean.pred<-unique(select(dat.p,-ab_duo,-rep_duo)) %>% 
  dplyr::group_by(Inoculum, Genus, CS_duo, CS_type_duo,CS1,CS2) %>%
  dplyr::summarize(ab_pred_mean=mean(ab_pred),
                   n=n(), sd.pred=sd(ab_pred), 
                   ab_pred_se=sd.pred/sqrt(n))

d.mean.pred<-as.data.frame(select(mean.pred,-n,-sd.pred))
d.mean.pred[is.na(d.mean.pred)] <- 0 

#. calculate mean of observed and SE
d.obs.mean = unique(select(dat.p,-rep_mono1,-rep_mono2,-ab_pred)) %>% 
  dplyr::group_by(Inoculum, Genus, CS_duo, CS_type_duo,CS1,CS2) %>%
  dplyr::summarize(ab_obs_mean=mean(ab_duo),
                   n=n(), sd.obs=sd(ab_duo), 
                   ab_obs_se=sd.obs/sqrt(n))

d.obs.mean<-select(as.data.frame(d.obs.mean), -n,-sd.obs)
d.obs.mean[is.na(d.obs.mean)] <- 0 

##. merge 2 dataframes
dat.pred.obs<-merge(d.obs.mean,d.mean.pred)

##. Plot
dp<-dat.pred.obs
p.title<-'Genus'
p.gen<-ggplot(dp, aes(x=ab_pred_mean, y=ab_obs_mean, group=CS1, colour=CS1,fill=CS1))+
  geom_point(alpha=0.5,size=2.5, shape=21) + 
  geom_errorbar(aes(ymin=ab_obs_mean-ab_obs_se, ymax=ab_obs_mean+ab_obs_se))+
  geom_errorbarh(aes(xmin=ab_pred_mean-ab_pred_se, xmax=ab_pred_mean+ab_pred_se))+
  theme_minimal() + 
  labs(x='predicted abundance' , y='observed abundance',title=p.title)+
  geom_abline(linetype='dashed')+
  scale_colour_manual(values = c( col.ss, col.aa))+
  scale_fill_manual(values = c( col.ss,col.aa))+
  guides(fill=guide_legend(title="Focal CS"))+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14),
        axis.title.y =element_text(size=16),
        axis.title.x = element_text(size = 16),
        plot.title= element_text(size = 18),
        panel.border = element_rect(color = "gray", fill = NA, size = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))+
  scale_x_continuous( breaks=seq(0,1,1), limits=c(0,1)) +
  scale_y_continuous( breaks=seq(0,1,1), limits=c(0,1))+
  guides( colour = FALSE)
p.gen

#.. write data figure as csv file
ds<-dp
fname<-paste("fig2S4_gen_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)

#. statistics-RMSE and Pearson's
fit.l <- lm(ab_obs_mean~ab_pred_mean, dp)
cor.test(dp$ab_obs_mean,dp$ab_pred_mean, method='pearson')
cor.gen<-cor(dp$ab_obs_mean,dp$ab_pred_mean)
rmse.gen<-rmse(dp$ab_obs_mean, dp$ab_pred_mean)

dat.glu<-subset(dp, CS1=='glucose')
dim(dat.glu)
fit.l.glu <- lm(ab_obs_mean~ab_pred_mean, dat.glu)
cor.test(dat.glu$ab_obs_mean,dat.glu$ab_pred_mean, method='pearson')
cor.gen.glu<-cor(dat.glu$ab_obs_mean,dat.glu$ab_pred_mean)
rmse.gen.glu<-rmse(dat.glu$ab_obs_mean, dat.glu$ab_pred_mean)

dat.succ<-subset(dp, CS1=='succinate')
dim(dat.succ)
fit.l.succ <- lm(ab_obs_mean~ab_pred_mean, dat.succ)
cor.test(dat.succ$ab_obs_mean,dat.succ$ab_pred_mean, method='pearson')
cor.gen.succ<-cor(dat.succ$ab_obs_mean,dat.succ$ab_pred_mean)
rmse.gen.succ<-rmse(dat.succ$ab_obs_mean, dat.succ$ab_pred_mean)

paste0('cor.gen:', round(cor.gen[1],3),'; rmse.gen:', round(rmse.gen[1],3))
paste0('cor.gen.glu:', round(cor.gen.glu[1],3),'; cor.gen.succ:', round(cor.gen.succ[1],3))
paste0('rmse.gen.glu:', round(rmse.gen.glu[1],3),'; rmse.gen.succ:', round(rmse.gen.succ[1],3))

## 3) ESV-LEVEL
#.. select data to plot
dat.p<-dat.d.e
dat.p<-select(dat.p,-od_duo, -od1_mean_mono,-od2_mean_mono,-ab1_mean_mono,-ab2_mean_mono,-CS1_mono,-CS2_mono,-ab1_corr,-ab2_corr)

#. calculate mean of predicted and SE
mean.pred<-unique(select(dat.p,-ab_duo,-rep_duo)) %>% 
  dplyr::group_by(Inoculum, ESV, CS_duo, CS_type_duo,CS1,CS2) %>%
  dplyr::summarize(ab_pred_mean=mean(ab_pred),
                   n=n(), sd.pred=sd(ab_pred), 
                   ab_pred_se=sd.pred/sqrt(n))

d.mean.pred<-as.data.frame(select(mean.pred,-n,-sd.pred))
d.mean.pred[is.na(d.mean.pred)] <- 0 

#. calculate mean of observed and SE
d.obs.mean = unique(select(dat.p,-rep_mono1,-rep_mono2,-ab_pred)) %>% 
  dplyr::group_by(Inoculum, ESV, CS_duo, CS_type_duo,CS1,CS2) %>%
  dplyr::summarize(ab_obs_mean=mean(ab_duo),
                   n=n(), sd.obs=sd(ab_duo), 
                   ab_obs_se=sd.obs/sqrt(n))

d.obs.mean<-select(as.data.frame(d.obs.mean), -n,-sd.obs)
d.obs.mean[is.na(d.obs.mean)] <- 0 

##. merge 2 dataframes
dat.pred.obs<-merge(d.obs.mean,d.mean.pred)

##. Plot
dp<-dat.pred.obs
p.title<-'ESV'
p.esv<-ggplot(dp, aes(x=ab_pred_mean, y=ab_obs_mean, group=CS1, colour=CS1,fill=CS1))+
  geom_point(alpha=0.5,size=2.5, shape=21) + 
  geom_errorbar(aes(ymin=ab_obs_mean-ab_obs_se, ymax=ab_obs_mean+ab_obs_se))+
  geom_errorbarh(aes(xmin=ab_pred_mean-ab_pred_se, xmax=ab_pred_mean+ab_pred_se))+
  theme_minimal() + 
  labs(x='predicted abundance' , y='observed abundance',title=p.title)+
  geom_abline(linetype='dashed')+
  scale_colour_manual(values = c( col.ss, col.aa))+
  scale_fill_manual(values = c( col.ss,col.aa))+
  guides(fill=guide_legend(title="Focal CS"))+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14),
        axis.title.y =element_text(size=16),
        axis.title.x = element_text(size = 16),
        plot.title = element_text(size = 18),
        panel.border = element_rect(color = "gray", fill = NA, size = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))+
  scale_x_continuous( breaks=seq(0,1,1), limits=c(0,1)) +
  scale_y_continuous( breaks=seq(0,1,1), limits=c(0,1))+
  guides( colour = FALSE)
p.esv

#.. write data figure as csv file
ds<-dp
fname<-paste("fig2S4_esv_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)

#. statistics- RMSE and Pearson's
fit.l <- lm(ab_obs_mean~ab_pred_mean, dp)
cor.test(dp$ab_obs_mean,dp$ab_pred_mean, method='pearson')
cor.esv<-cor(dp$ab_obs_mean,dp$ab_pred_mean)
rmse.esv<-rmse(dp$ab_obs_mean, dp$ab_pred_mean)

dat.glu<-subset(dp, CS1=='glucose')
dim(dat.glu)
fit.l.glu <- lm(ab_obs_mean~ab_pred_mean, dat.glu)
cor.test(dat.glu$ab_obs_mean,dat.glu$ab_pred_mean, method='pearson')
cor.esv.glu<-cor(dat.glu$ab_obs_mean,dat.glu$ab_pred_mean)
rmse.esv.glu<-rmse(dat.glu$ab_obs_mean, dat.glu$ab_pred_mean)

dat.succ<-subset(dp, CS1=='succinate')
dim(dat.succ)
fit.l.succ <- lm(ab_obs_mean~ab_pred_mean, dat.succ)
cor.test(dat.succ$ab_obs_mean,dat.succ$ab_pred_mean, method='pearson')
cor.esv.succ<-cor(dat.succ$ab_obs_mean,dat.succ$ab_pred_mean)
rmse.esv.succ<-rmse(dat.succ$ab_obs_mean, dat.succ$ab_pred_mean)

paste0('cor.esv:', round(cor.esv[1],3),'; rmse.esv:', round(rmse.esv[1],3))
paste0('cor.esv.glu:', round(cor.esv.glu[1],3),'; cor.esv.succ:', round(cor.esv.succ[1],3))
paste0('rmse.esv.glu:', round(rmse.esv.glu[1],3),'; rmse.esv.succ:', round(rmse.esv.succ[1],3))

ggarrange(p.fam,p.gen,p.esv, ncol=3, nrow=1, common.legend = TRUE, legend="right")
ggsave(path=path.p, "fig2S4.pdf", width=12, height=3.9, onefile = TRUE)







