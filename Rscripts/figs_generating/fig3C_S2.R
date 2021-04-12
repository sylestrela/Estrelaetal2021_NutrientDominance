" R script for Figure 3C and Figure 3- figure supplement 2 in Estrela et al (2021) 
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

dat.p<-select(dat,-od_duo,-od1_mean_mono,-od2_mean_mono,-CS1_mono,-CS2_mono,-ab1_corr,-ab2_corr,-ab_pred)

#. calculate mean of abundance in singles (observed)
mean.pred<-dat.p %>% 
  dplyr::group_by(Inoculum, Family, CS_duo, CS_type_duo) %>%
  dplyr::summarize(
    ab1_mean_mono=mean(ab1_mean_mono),
    ab2_mean_mono=mean(ab2_mean_mono))
d.mean.pred<-as.data.frame(mean.pred)

#. merge with observed abundance in pairs 
d.obs.mean <-unique(select(dat.p,-rep_mono1,-rep_mono2,-ab1_mean_mono,-ab2_mean_mono)) 
dat.pred.obs<-merge(d.obs.mean,d.mean.pred,by=c('Family','Inoculum', 'CS_duo','CS_type_duo'))

##-------------------------------------------------
## Outcome of interaction (Family-level)
##-------------------------------------------------

##. 1) get interaction (epsilon) and significance determined previously
d.epsi<-read.csv('data_processed/ttest_paired_outcome_family_all.csv')

#. combine observed abundances and epsilon 
dt<-merge(dat.pred.obs, select(d.epsi,-CS1_mono,-CS2_mono), by=c('Family','CS_duo','CS_type_duo'))

#. get no interaction group
dt$outcome<-with(dt, ifelse(t.outcome=='ns', 'no interaction',
                            'dominance'))
dt<-na.omit(dt)

##-------------------------------------------------
## identify which CS dominates for each Inocula (ie. lower dom_cs) 
##-------------------------------------------------

dt$cs1_type<- 'focal'
dt$cs2_type<- 'additional'

dt$dom_cs <- with(dt, ifelse(ab1_mean_mono>ab2_mean_mono & epsi_mean>0, 'CS1',
                             ifelse(ab1_mean_mono>ab2_mean_mono & epsi_mean<0, 'CS2',
                                    ifelse(ab1_mean_mono<ab2_mean_mono & epsi_mean>0, 'CS2',
                                           ifelse(ab1_mean_mono<ab2_mean_mono & epsi_mean<0, 'CS1','none')))))

dt1<-subset(dt,dom_cs=='none' )
dt1$outcome<-'no interaction'
dt2<-subset(dt,dom_cs!='none' )
#. 
no.interaction<-rbind(subset(dt, outcome=='no interaction'), dt1)
dominance<-subset(dt2, outcome=='dominance')

#. make sure that dominant CS in the 2 inocula is the same, if not then 'no interaction'.
dt<-dominance
dt$dom_cs_v<-with(dt, ifelse(dom_cs=='CS1', 0, 1))

dt<-dt %>%
  dplyr::group_by(Family,CS_duo) %>%
  dplyr::mutate(val=mean(dom_cs_v))

dt$outcome<-with(dt, ifelse(val %in% c(0,1), 'dominance',
                            'no interaction'))

dt<-as.data.frame(select(dt,-dom_cs_v,-val))
dominance<-as.data.frame(subset(dt, outcome=='dominance'))

#. update dataframe for no interaction
no.interaction<-rbind(no.interaction, subset(dt, outcome=='no interaction'))

#. test for synergy and antagonism. synergy and antagonism occur when f12> max(f1,f2) or f12<min(f1,f2)
dominance$dmax<-pmax(dominance$ab1_mean_mono, dominance$ab2_mean_mono)
dominance$dmin<-pmin(dominance$ab1_mean_mono, dominance$ab2_mean_mono)

my.t.test.p.val.greater <- function(...) { 
  obj<-try(t.test(..., alternative = "greater"), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
} 

my.t.test.p.val.less <- function(...) { 
  obj<-try(t.test(..., alternative='less'), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
} 

my.t.test.p.val.greater(numeric(0)) 

t.test2 = dominance %>%
  dplyr::group_by(Family, CS_duo) %>%
  dplyr::mutate(p.val.s=my.t.test.p.val.greater(abs(ab_duo), abs(dmax)),
                p.val.a=my.t.test.p.val.less(abs(ab_duo), abs(dmin)))
t.test2<-na.omit(t.test2)

#. get superdominance (synergy or antagonism group)
#. superdominance occurs when p-val< 0.05
t.test2$outcome<-with(t.test2, ifelse(p.val.s<0.05, 'synergy',
                                      ifelse(p.val.a<0.05 , 'antagonism',
                                             'dominance')))

dominance<-as.data.frame(t.test2)
dim(dominance)

##. combine all
names(no.interaction)
names(dominance)

#. create empty colums for no.interaction so that binds dominance and superdominance
no.interaction$p.val.s<-0
no.interaction$p.val.a<-0
dominance<-select(dominance,-dmax,-dmin)

d.all<-as.data.frame(rbind(no.interaction,dominance))
d.all<-select(d.all,-cs1_type,-cs2_type)

d.all$outcome2 <- with(d.all, ifelse(epsi_mean >0 & outcome=='dominance', "+ dominance", 
                                     ifelse( outcome=='synergy', "synergy", 
                                             ifelse(epsi_mean <0 & outcome=='dominance', "- dominance", 
                                                    ifelse(outcome=='antagonism', "antagonism", 
                                                           ifelse(outcome=='no interaction', "no interaction", 
                                                                  0))))))


##-------------------------------------------------
## Plot Figure 3C: outcome of dominance for SA
##-------------------------------------------------

dt<-subset(d.all, CS_type_duo=='SA')
dt$cs1_type<- with(dt, ifelse(CS1=='glucose', 'sugar' , 'acid'  ))
dt$cs2_type<- with(dt, ifelse(CS2 %in% c("glycerol","fructose", "cellobiose", "ribose" ), 'sugar', 'acid'))

dt$dom_cs_v <- with(dt, ifelse(dom_cs=='CS1' & cs1_type=='sugar', -abs(epsi_mean) , 
                               ifelse(dom_cs=='CS2' & cs2_type=='sugar', -abs(epsi_mean), 
                                      abs(epsi_mean) )))
dt$dom_cs_type <- with(dt, ifelse(dom_cs=='CS1', cs1_type , cs2_type  ))

#.. read data for the individual pairs (N=8) of one permutation test 
d.epsi.p<-read.csv('data_processed/epsi_paired_EPMR.csv')
d.epsi.p<-subset(d.epsi.p, CS_type_duo=='SA')
d.epsi.p<-select(d.epsi.p, -ab1_corr,-ab2_corr,-ab_pred,-ab_duo,-rep_mono1,-rep_mono2,-rep_duo,-rep_combo_inoc)
d.epsi.p$cs1_type<- with(d.epsi.p, ifelse(CS1_mono=='glucose', 'sugar' , 'acid'  ))
d.epsi.p$cs2_type<- with(d.epsi.p, ifelse(CS2_mono %in% c("glycerol","fructose", "cellobiose", "ribose" ), 'sugar', 'acid'))

d.epsi.p$dom_cs_v2 <- with(d.epsi.p, ifelse(dom_cs=='CS1' & cs1_type=='sugar', -abs(epsi) , 
                                            ifelse(dom_cs=='CS2' & cs2_type=='sugar', -abs(epsi), 
                                                   abs(epsi) )))

d.epsi.p$dom_cs_type <- with(d.epsi.p, ifelse(dom_cs=='CS1', cs1_type , cs2_type))

dt1<-subset(dt, Family %in% c('Enterobacteriaceae', 'Moraxellaceae','Pseudomonadaceae', 'Rhizobiaceae'))
#. define order of families to plot
dt1$Family<-factor(dt1$Family, levels = c('Enterobacteriaceae', 'Pseudomonadaceae','Moraxellaceae', 'Rhizobiaceae'))

#. define order of CS names to plot
dt1$CS_duo <- factor(dt1$CS_duo, levels = c("succinate_glycerol","succinate_ribose", "succinate_fructose", "succinate_cellobiose",
                                            "glucose_glutamine","glucose_glycine", "glucose_fumarate", "glucose_benzoate"))
d.epsi.p$CS_duo <- factor(d.epsi.p$CS_duo, levels = c("succinate_glycerol","succinate_ribose", "succinate_fructose", "succinate_cellobiose",
                                                      "glucose_glutamine","glucose_glycine", "glucose_fumarate", "glucose_benzoate"))

col.sd<-'mediumpurple1'
col.ss<-'mediumpurple4'
col.ad<-'tan1'
col.as<-'tan4'
col.ns<-'gray'

ggplot(dt1, aes(x = dom_cs_v, y = CS_duo))+
  geom_point(data=d.epsi.p, aes(x = dom_cs_v2, y = CS_duo, shape = factor(Inoculum),color = factor(dom_cs_type)),alpha=0.4)+
  geom_point(data=dt1, aes(color = factor(dom_cs_type):factor(outcome2)),alpha=0.8,size=2.5)+
  geom_errorbarh(data=dt1, aes(xmax = dom_cs_v + epsi_sd, xmin = dom_cs_v - epsi_sd,height=0,color = factor(dom_cs_type):factor(outcome2)))+
  theme_minimal()+
  facet_grid( Family~CS_type_duo)+
  geom_vline(aes(xintercept=0), linetype='dashed')+
  guides(colour=guide_legend(title="Dominant CS"),shape=guide_legend(title="Inoculum"))+
  scale_colour_manual(values = c(col.ad, col.ad, col.ad, col.ns,col.sd, col.sd,col.sd, col.ss,col.ns))+
  scale_shape_manual(values = c(5,6)) +
  labs(x=expression(delta), y='')+
  #   caption= paste0( 'For each CS in pair, d= (1- |xi(obs, mix)- mean(xi(obs, single)| )* |epsi|. CS with greater delta dominates. Only 4 most abundant families are shown.'))+
  theme(axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        axis.title.y =element_text(size=14),
        axis.title.x = element_text(size = 18),
        panel.border = element_rect(color = "gray", fill = NA, size = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
ggsave(path=path.p, "fig3C.pdf", width=6, height=8.5, onefile = TRUE)

#.. write data figure as csv file
ds<-unique(select(dt1, Family, CS_duo, CS_type_duo, dom_cs_type, outcome2, dom_cs_v,epsi_sd))
fname<-paste("fig3C_data.csv")
write.csv(ds, file.path(path.ds, fname), row.names=FALSE)

##-------------------------------------------------
## Plot Figure 3S2: outcome of dominance for AA and SS
##-------------------------------------------------
dt<-subset(d.all, CS_type_duo!='SA')
dt$cs1_type<- 'focal'
dt$cs2_type<- 'additional'

dt$dom_cs_v <- with(dt, ifelse(dom_cs=='CS1', -abs(epsi_mean) , abs(epsi_mean)))
dt$dom_cs_type <- with(dt, ifelse(dom_cs=='CS1', cs1_type , cs2_type  ))

#.. read data for individual pairs
d.epsi.p<-read.csv('data_processed/epsi_paired_EPMR.csv')
d.epsi.p<-subset(d.epsi.p, CS_type_duo!='SA')
d.epsi.p<-select(d.epsi.p, -ab1_corr,-ab2_corr,-ab_pred,-ab_duo,-rep_mono1,-rep_mono2,-rep_duo,-rep_combo_inoc)

d.epsi.p$dom_cs_v2 <- with(d.epsi.p, ifelse(dom_cs=='CS1', -abs(epsi) , 
                                            abs(epsi) ))

d.epsi.p$dom_cs_type <- with(d.epsi.p, ifelse(dom_cs=='CS1', as.character(cs1_type) , as.character(cs2_type)))

dp<-subset(dt, Family %in% c('Enterobacteriaceae', 'Moraxellaceae','Pseudomonadaceae', 'Rhizobiaceae'))
#. define order of CS names to plot
dp$CS_duo <- factor(dp$CS_duo, levels = c("glucose_fructose","glucose_cellobiose", "glucose_ribose", "glucose_glycerol",
                                          "succinate_benzoate","succinate_fumarate","succinate_glycine","succinate_glutamine"))
d.epsi.p$CS_duo <- factor(d.epsi.p$CS_duo, levels = c("glucose_fructose","glucose_cellobiose", "glucose_ribose", "glucose_glycerol",
                                                      "succinate_benzoate","succinate_fumarate","succinate_glycine","succinate_glutamine"))

ggplot(dp, aes(x = dom_cs_v, y = CS_duo))+
  geom_point(data=d.epsi.p, aes(x = dom_cs_v2, y = CS_duo, shape = factor(Inoculum),color = factor(dom_cs_type)),alpha=0.9,size=1.8)+
  geom_point(data=dp, aes(color = factor(dom_cs_type):factor(outcome2)),alpha=1,size=2.5)+
  geom_errorbarh(data=dp, aes(xmax = dom_cs_v + epsi_sd, xmin = dom_cs_v - epsi_sd,height=0,color = factor(dom_cs_type):factor(outcome2)))+
  theme_minimal()+
  facet_wrap(.~Family, scales='free_x', ncol=4)+
  geom_vline(aes(xintercept=0), linetype='dashed')+
  guides(colour=guide_legend(title="Dominant CS"),shape=guide_legend(title="Inoculum"))+
  scale_colour_manual(values = c( col.ad, col.ad,col.ad, col.as,col.ns,col.sd, col.ss,col.ns,col.ss))+
  scale_shape_manual(values = c(1,2)) +
  labs(x=expression(delta), y='')+
  scale_x_continuous(breaks = seq(-0.2, 0.2, by = 0.2))+
  #   caption= paste0( 'For each CS in pair, d= (1- |xi(obs, mix)- mean(xi(obs, single)| )* |epsi|. CS with greater delta dominates. Only 4 most abundant families are shown.'))+
  theme(axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        axis.title.y =element_text(size=14),
        axis.title.x = element_text(size = 18),
        panel.border = element_rect(color = "gray", fill = NA, size = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
ggsave(path=path.p, "fig3S2.pdf", width=11.3, height=4.2, onefile = TRUE)

#.. write data figure as csv file
ds<-select(dp, Family, Inoculum, CS_duo, CS_type_duo, dom_cs_type, outcome2, dom_cs_v, epsi_sd)
fname<-paste("fig3S2_data.csv")
write.csv(ds, file.path(path.ds, fname), row.names=FALSE)




