" R script for Figure 3- figure supplements 3 and 4 in Estrela et al (2021)
Nutrient dominance governs the assembly of microbial communities in mixed nutrient environments.

Genus-level outcome of interaction (no interaction, dominance, superdominance) 
determined using paired t-test significance level (inocula analysed separately).
note: only pairs with n>=3 are considered, see companion script:
stats_permutation_pairedttest_gen_ALL.R"

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
dat<-read.csv('data_processed/data_mixedCS_obs_pred_gen_duos.csv')

dat.p<-select(dat,-od_duo,-od1_mean_mono,-od2_mean_mono,-CS1_mono,-CS2_mono,-ab1_corr,-ab2_corr,-ab_pred)

#. calculate mean of abundance in singles (observed)
mean.pred<-dat.p %>% 
  dplyr::group_by(Inoculum, Genus, CS_duo, CS_type_duo) %>%
  dplyr::summarize(
    ab1_mean_mono=mean(ab1_mean_mono),
    ab2_mean_mono=mean(ab2_mean_mono))
d.mean.pred<-as.data.frame(mean.pred)

#. merge with observed abundance in pairs 
d.obs.mean <-unique(select(dat.p,-rep_mono1,-rep_mono2,-ab1_mean_mono,-ab2_mean_mono)) 
dat.pred.obs<-merge(d.obs.mean,d.mean.pred,by=c('Genus','Inoculum', 'CS_duo','CS_type_duo'))

##-------------------------------------------------
## Outcome of interaction (Genus-level)
##-------------------------------------------------

##. 1) get interaction (epsilon) and significance determined perviously
d.epsi<-read.csv('data_processed/ttest_paired_outcome_genus_all.csv')

#. combine observed abundances and epsilon 
dt<-merge(dat.pred.obs, select(d.epsi,-CS1_mono,-CS2_mono), by=c('Inoculum', 'Genus','CS_duo','CS_type_duo'))

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

no.interaction<-rbind(subset(dt, outcome=='no interaction'), dt1)
dominance<-subset(dt2, outcome=='dominance')

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
  dplyr::group_by(Genus, CS_duo) %>%
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
  dplyr::distinct(Inoculum, Genus, CS_duo,.keep_all = TRUE) 
d.all2<-as.data.frame(d.all2)

dpo1<-as.data.frame(table(d.all2$outcome))
dpo2<-as.data.frame(table(d.all2$outcome2))

dpo1$frac<-round(dpo1$Freq/sum(dpo1$Freq)*100,1)
dpo2$frac<-round(dpo2$Freq/sum(dpo2$Freq)*100,1)

##-------------------------------------------------
## Plot figure 3S3: outcome of interaction at genus-level
##-------------------------------------------------

##.. Plot 3S3A
dp1<-dpo2
#. define order of interactions for plot
levels(factor(dp1$Var1))
dp1$type.ord <- factor(dp1$Var1, levels = c( 'synergy', '+ dominance', '- dominance','antagonism', 'no interaction'))

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

# Apply blank theme
library(scales)
p1<-pie + scale_fill_grey() +  blank_theme +
  theme(axis.text.x=element_blank()) +
  scale_fill_manual(values=c( 'steelblue4', "steelblue","indianred", "indianred4", 'gray'))+
  guides(fill=guide_legend(title=""))+
  theme(legend.text =element_text(size=18),
        legend.title =element_text(size=18))
p1
ggsave(path=path.p, "fig3S3A.pdf", width=10, height=3, onefile = TRUE)

#.. write data figure as csv file
ds<-select(dp1,-type.ord,-Freq)
fname<-paste("fig3S3A_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)

##.. Plot 3S3B
#.. by Carbon source type
dpo.cs1<-d.all2 %>%
  dplyr::group_by(CS_type_duo, outcome2) %>%
  dplyr::summarise(freq=table(outcome2))%>%
  dplyr::mutate(frac= freq/sum(freq)*100)
dpo.cs1

dpo.cs1$type.ord <- factor(dpo.cs1$outcome2, levels = c( 'synergy', '+ dominance', '- dominance','antagonism', 'no interaction'))

p2<-ggplot(dpo.cs1, aes(x=CS_type_duo, y=frac, fill=type.ord))+
  geom_bar(stat="identity")+
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
p2
ggsave(path=path.p, "fig3S3B.pdf", width=4.5, height=3.4, onefile = TRUE)

#.. write data figure as csv file
ds<-select(as.data.frame(dpo.cs1),-type.ord)
fname<-paste("fig3S3B_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)

##.. Plot 3S3C
#. by Genus - select top 10 (confirm that spans top 4 Families
dp.gen<-d.all2 %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarise(sum=sum(ab_duo))
dp.gen<-as.data.frame(dp.gen)
dp.gen<-dp.gen[order(-dp.gen$sum),] 

##.. select top 10 genera
dp.gen.10<-dp.gen$Genus[c(1:10)]
dpo.gen<-d.all2
levels(dpo.gen$Genus)[!(levels(dpo.gen$Genus) %in% dp.gen.10)] <- 'other'

dpo.gen<-dpo.gen %>%
  dplyr::group_by(Genus, CS_type_duo, outcome2) %>%
  dplyr::summarise(freq=table(outcome2))%>%
  dplyr::mutate(frac= freq/sum(freq)*100)
dpo.gen

dpo.gen$type.ord <- factor(dpo.gen$outcome2, levels = c( 'synergy', '+ dominance', '- dominance','antagonism', 'no interaction'))
dpo.gen$gen.ord <- factor(dpo.gen$Genus, levels = c('Raoultella', 'Klebsiella', 'Enterobacteriaceae.8', 'Rosenbergiella', 'Pectobacterium', 'Enterobacter', 'Citrobacter', 'Pseudomonas', 'Acinetobacter', 'Shinella', 'other'))

p3<-ggplot(dpo.gen, aes(x=CS_type_duo, y=frac, fill=type.ord))+
  geom_bar(stat="identity")+
  facet_wrap(.~gen.ord, ncol=7)+
  scale_fill_manual(values=c(  'steelblue4', "steelblue","indianred", "indianred4", 'gray'))+
  labs(y='Fraction',x='')+
  guides(fill=guide_legend(title=""))+
  theme_classic()+
  theme(axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        strip.text.x = element_text(size=10),
        strip.text.y = element_text(size=12),
        axis.title.y =element_text(size=18),
        legend.text =element_text(size=12),
        legend.title =element_text(size=12),
        axis.title.x = element_text(size = 16))
p3
ggsave(path=path.p, "fig3S3C.pdf", width=13, height=6, onefile = TRUE)

#.. write data figure as csv file
ds<-select(as.data.frame(dpo.gen),-type.ord,-freq)
fname<-paste("fig3S3C_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)


##-------------------------------------------------
## Plot figure 3S4: outcome of dominance as density plot for SA
##-------------------------------------------------
dt<-subset(d.all2, CS_type_duo=='SA')
dt$cs1_type<- with(dt, ifelse(CS1=='glucose', 'sugar' , 'acid'  ))
dt$cs2_type<- with(dt, ifelse(CS2 %in% c("glycerol","fructose", "cellobiose", "ribose" ), 'sugar', 'acid'))

dt$dom_cs_v <- with(dt, ifelse(dom_cs=='CS1' & cs1_type=='sugar', -abs(epsi_mean) , 
                               ifelse(dom_cs=='CS2' & cs2_type=='sugar', -abs(epsi_mean), 
                                      abs(epsi_mean) )))

dt$dom_cs_type <- with(dt, ifelse(dom_cs=='CS1', cs1_type , cs2_type  ))

dt1<-dt
dp<-subset(dt1, Genus %in% dp.gen.10)
dp$gen.ord <- factor(dp$Genus, levels = c('Raoultella', 'Klebsiella', 'Enterobacteriaceae.8', 'Rosenbergiella', 'Pectobacterium', 'Enterobacter', 'Citrobacter', 'Pseudomonas', 'Acinetobacter', 'Shinella'))

#. define order of CS names to plot
dp$CS_duo <- factor(dp$CS_duo, levels = c("succinate_glycerol","succinate_ribose", "succinate_fructose", "succinate_cellobiose",
                                          "glucose_glutamine","glucose_glycine", "glucose_fumarate", "glucose_benzoate"))

# .. Colour scheme for plots
col.sd<-'mediumpurple1'
col.ss<-'mediumpurple4'
col.ad<-'tan1'
col.as<-'tan4'
col.ns<-'gray'

ggplot(dp, aes(x = dom_cs_v, y = CS_duo, color = factor(dom_cs_type):factor(outcome2), shape=factor(Inoculum)))+
  geom_point(alpha=1,size=3)+
  geom_errorbarh(aes(xmax = dom_cs_v + epsi_sd, xmin = dom_cs_v - epsi_sd,height=0))+
  theme_minimal()+
  facet_wrap( .~gen.ord, ncol=5)+
  geom_vline(aes(xintercept=0), linetype='dashed')+
  guides(colour=guide_legend(title="Dominant CS"),shape=guide_legend(title="Inoculum"))+
  scale_colour_manual(values = c(col.ad, col.ad, col.as, col.ns, col.as, col.sd,col.sd, col.ns,col.ss))+
  labs(x=expression(delta), y='')+
  theme(axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        strip.text.x = element_text(size=13),
        strip.text.y = element_text(size=12),
        axis.title.y =element_text(size=14),
        axis.title.x = element_text(size = 18),
        panel.border = element_rect(color = "gray", fill = NA, size = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
ggsave(path=path.p, "fig3S4.pdf", width=13, height=6, onefile = TRUE)

#.. write data figure as csv file
ds<-select(dp, Genus, Inoculum, CS_duo, CS_type_duo, dom_cs_type, outcome2, dom_cs_v, epsi_sd,gen.ord)
fname<-paste("fig3S4_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)


