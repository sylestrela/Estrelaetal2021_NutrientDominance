" R script for Figure 3- figure supplement 1 in Estrela et al (2021)
Nutrient dominance governs the assembly of microbial communities in mixed nutrient environments.

Family-level outcome of interaction (no interaction, dominance, superdominance) 
determined using paired t-test significance level. Inocula are pooled together after pairing.
note: only pairs with n>=6 are considered, see companion script: stats_permutation_pairedttest_fam_ALL.R"

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

##. 1) get interaction (epsilon) and significance determined perviously
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

#. filter out replicates
d.all2 = d.all %>%
  dplyr::distinct(Family, CS_duo,.keep_all = TRUE) 
d.all2<-as.data.frame(d.all2)

dpo1<-as.data.frame(table(d.all2$outcome))
dpo2<-as.data.frame(table(d.all2$outcome2))

dpo1$frac<-round(dpo1$Freq/sum(dpo1$Freq)*100,1)
dpo1$Families<-'All'

dpo2$frac<-round(dpo2$Freq/sum(dpo2$Freq)*100,1)
dpo2$Families<-'All'

#. Select the top 4 dominant families only
dps<-subset(d.all2, Family %in% c('Enterobacteriaceae', 'Moraxellaceae','Pseudomonadaceae', 'Rhizobiaceae'))
dpo3<-as.data.frame(table(dps$outcome))
dpo4<-as.data.frame(table(dps$outcome2))

dpo3$frac<-round(dpo3$Freq/sum(dpo3$Freq)*100,1)
dpo4$frac<-round(dpo4$Freq/sum(dpo4$Freq)*100,1)

dpo3$Families<-'Dominant four'
dpo4$Families<-'Dominant four'

#. merge all and dominant families
dp1<-rbind(dpo2,dpo4)

#. define order of interactions for plot
dp1$type.ord <- factor(dp1$Var1, levels = c( 'synergy', '+ dominance', '- dominance','antagonism', 'no interaction'))

##. Plot Fig3S1A
dp1<-subset(dp1, Families=='All')                     
bp<-ggplot(dp1, aes(x='', y=frac, fill=type.ord))+geom_bar(stat="identity")
pie <- bp + coord_polar("y", start=0)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

#. Apply blank theme
library(scales)
p1<-pie + scale_fill_grey() +  blank_theme +
  theme(axis.text.x=element_blank()) +
  scale_fill_manual(values=c(  'steelblue4', "steelblue","indianred", "indianred4", 'gray'))+
  guides(fill=guide_legend(title=""))+
  theme(legend.text =element_text(size=18),
        legend.title =element_text(size=18))
p1
ggsave(path=path.p, "fig3S1A.pdf", width=10, height=3, onefile = TRUE)

#.. write data figure as csv file
ds<-select(dp1,-type.ord,-freq)
fname<-paste("fig3S1A_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)

##. Plot Fig3S1B
#.. by Carbon source type
dpo.cs1<-d.all2 %>%
  dplyr::group_by(CS_type_duo, outcome2) %>%
  dplyr::summarise(freq=table(outcome2))%>%
  dplyr::mutate(frac= freq/sum(freq)*100)
dpo.cs1$Families<-'All families'

dpo.cs1$type.ord <- factor(dpo.cs1$outcome2, levels = c( 'synergy', '+ dominance', '- dominance','antagonism', 'no interaction'))

p2<-ggplot(dpo.cs1, aes(x=CS_type_duo, y=frac, fill=type.ord))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c(  'steelblue4', "steelblue","indianred", "indianred4", 'gray'))+
  labs(y='Fraction',x='')+
  guides(fill=guide_legend(title=""))+
  theme_classic()+
  theme(axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        axis.title.y =element_text(size=18),
        legend.text =element_text(size=12),
        legend.title =element_text(size=12),
        axis.title.x = element_text(size = 16))
p2
ggsave(path=path.p, "fig3S1B.pdf", width=4.5, height=3.4, onefile = TRUE)

#.. write data figure as csv file
ds<-select(as.data.frame(dpo.cs1),-type.ord,-freq)
fname<-paste("fig3S1B_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)

##. Plot Fig3S1C
#. by Families (top 4 + others)
dpo.fam<-d.all2
levels(dpo.fam$Family)[!(levels(dpo.fam$Family) %in% c('Enterobacteriaceae', 'Moraxellaceae','Pseudomonadaceae', 'Rhizobiaceae'))] <- 'other'

dpo.fam<-dpo.fam %>%
  dplyr::group_by(Family, CS_type_duo, outcome2) %>%
  dplyr::summarise(freq=table(outcome2))%>%
  dplyr::mutate(frac= freq/sum(freq)*100)
as.data.frame(dpo.fam)

levels(factor(dpo.fam$outcome2))
dpo.fam$type.ord <- factor(dpo.fam$outcome2, levels = c( 'synergy', '+ dominance', '- dominance','antagonism', 'no interaction'))
dpo.fam$fam.ord <- factor(dpo.fam$Family, levels = c('Enterobacteriaceae','Pseudomonadaceae', 'Moraxellaceae', 'Rhizobiaceae', 'other'))

p3<-ggplot(dpo.fam, aes(x=CS_type_duo, y=frac, fill=type.ord))+
  geom_bar(stat="identity")+
  facet_grid(.~fam.ord)+
  scale_fill_manual(values=c( 'steelblue4', "steelblue","indianred", "indianred4", 'gray'))+
  labs(y='Fraction',x='')+
  guides(fill=guide_legend(title=""))+
  theme_classic()+
  theme(axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        axis.title.y =element_text(size=18),
        legend.text =element_text(size=12),
        legend.title =element_text(size=12),
        axis.title.x = element_text(size = 16))
p3
ggsave(path=path.p, "fig3S1C.pdf", width=11, height=3.5, onefile = TRUE)

#.. write data figure as csv file
ds<-select(as.data.frame(dpo.fam),-type.ord,-freq)
fname<-paste("fig3S1C_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)





