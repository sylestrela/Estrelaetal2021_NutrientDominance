" R script used in Estrela et al (2021)
Nutrient dominance governs the assembly of microbial communities in mixed nutrient environments.

Script to calculate the INDIVIDUAL epsilons for one permutation pair 
and significance level for 4 dominant families in all CS 
(values to be used for plots in Fig3C and Fig. 3S2).
T-test: one-tailed paired t-test (run with alternative= 'greater' and then 'less')
No permutations, values reflect single pairing."

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

bs.rep<-1
paired.v<-TRUE
t.type<-"greater" #.. note: using greater or less does not matter here as we just want the mean difference

## ------------------------------
#.. FAMILY-LEVEL 
dat<-read.csv('data_processed/data_mixedCS_obs_pred_fam_duos.csv')
taxa.level<-'family'
dsp<-select(dat, Inoculum, Family, CS_type_duo, CS_duo, CS1_mono, CS2_mono, ab1_corr, ab2_corr, ab_pred, ab_duo,rep_mono1,rep_mono2,rep_duo)
dsp<-subset(dsp, Family %in% c('Enterobacteriaceae', 'Pseudomonadaceae', 'Rhizobiaceae', 'Moraxellaceae'))

##.. test whether Ent, Pseu, Morax and Rhiz are significantly MORE or LESS abundant than predicted by the null model in succinate-sugar pairs
dsp$rep_combo_inoc<-paste0(dsp$rep_mono1,dsp$rep_mono2, dsp$rep_duo,"_", dsp$Inoculum)
#. for all CS pairs except glycine (see script below for glycine)
ds<-subset(dsp, CS2_mono!='glycine')

#.. Permutation test
#.. create 8 unique singleCS-singleCS- pairedCS combinations keeping the identity of the inocula for each pair 
datalist = list()
set.seed(10)

for (i in 1:bs.rep) {
  
  ds1<-sample(unique(ds$rep_mono1), size=4, replace=FALSE)
  ds2<-sample(unique(ds$rep_mono2), size=4, replace=FALSE)
  ds3<-sample(unique(ds$rep_duo), size=4, replace=FALSE)
  
  rep.combo.1<-c(paste0(ds1[1], ds2[1], ds3[1],"_1"), 
                 paste0(ds1[2], ds2[2], ds3[2],"_1"),
                 paste0(ds1[3], ds2[3], ds3[3],"_1"),
                 paste0(ds1[4], ds2[4], ds3[4],"_1"))
  
  rep.combo.2<-c(paste0(ds1[1], ds2[1], ds3[1],"_2"), 
                 paste0(ds1[2], ds2[2], ds3[2],"_2"),
                 paste0(ds1[3], ds2[3], ds3[3],"_2"),
                 paste0(ds1[4], ds2[4], ds3[4],"_2"))
  
  pairing.1<-subset(ds, rep_combo_inoc %in% rep.combo.1)
  pairing.2<-subset(ds, rep_combo_inoc %in% rep.combo.2)
  pairing<-rbind(pairing.1,pairing.2)
  
}
dp1<-pairing

#  ---------------------------------------------------------------------------
##.. now test for GLYCINE only (there are 3 replicates per CS)
#.. create 3 unique singleCS-singleCS- pairedCS combinations 
ds<-subset(dsp, CS2_mono=='glycine')
datalist = list()

for (i in 1:bs.rep) {
  
  ds1<-sample(unique(ds$rep_mono1), size=3, replace=FALSE)
  ds2<-sample(unique(ds$rep_mono2), size=3, replace=FALSE)
  ds3<-sample(unique(ds$rep_duo), size=3, replace=FALSE)
  
  rep.combo.1<-c(paste0(ds1[1], ds2[1], ds3[1],"_1"), 
                 paste0(ds1[2], ds2[2], ds3[2],"_1"),
                 paste0(ds1[3], ds2[3], ds3[3],"_1"))
  
  rep.combo.2<-c(paste0(ds1[1], ds2[1], ds3[1],"_2"), 
                 paste0(ds1[2], ds2[2], ds3[2],"_2"),
                 paste0(ds1[3], ds2[3], ds3[3],"_2"))
  
  pairing.1<-subset(ds, rep_combo_inoc %in% rep.combo.1)
  pairing.2<-subset(ds, rep_combo_inoc %in% rep.combo.2)
  pairing<-rbind(pairing.1,pairing.2)
  
}

dp2<-pairing

#.. merge glycine and all other CS data
dp.all<-rbind(dp1,dp2)
dp.all$epsi<-dp.all$ab_duo-dp.all$ab_pred

##-------------------------------------------------
## identify which CS dominates for each Inocula (ie. lower dom_cs) 
##-------------------------------------------------
dt<-dp.all
dt$cs1_type<- 'focal'
dt$cs2_type<- 'additional'
dt$dom_cs <- with(dt, ifelse(ab1_corr>ab2_corr & epsi>0, 'CS1',
                             ifelse(ab1_corr>ab2_corr & epsi<0, 'CS2',
                                    ifelse(ab1_corr<ab2_corr & epsi>0, 'CS2',
                                           ifelse(ab1_corr<ab2_corr & epsi<0, 'CS1','none')))))

#.. write data as csv file for each paired t-test
#fname<-paste0("epsi_paired_EPMR.csv")
#write.csv(dt, file.path(path.d, fname), row.names=FALSE)

