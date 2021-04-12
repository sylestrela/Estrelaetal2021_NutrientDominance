" R script for analysis in Estrela et al (2021)
Nutrient dominance governs the assembly of microbial communities in mixed nutrient environments.

Script to determine whether the observed abundance in the CS mixture is 
significantly greater or lower than the predicted (additive) abundance.
Analysis for FAMILY-level.
T-test: one-tailed paired t-test (run first with alternative= 'greater' and then 'less') 
where the two inocula are pooled together but pairing is done within each inoculum. 
Extract t-statistic for each test.
Then perform M=1000 permutations of t-test to determine if effect is significant (95% CI). 

Glycine is analysed separately because it only has 3 replicates per CS (unlike 4 reps for other CS)."

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

##.. bootstrap and t-test params
bs.rep<-1000 #.. select number of Permutation tests
paired.v<-TRUE #.. select type of t-test
t.type<-"greater" #.. select alternative type ('greater' or 'less')

## ------------------------------
#.. FAMILY-LEVEL 
dat<-read.csv('data_processed/data_mixedCS_obs_pred_fam_duos.csv')
taxa.level<-'family'
dsp<-select(dat, Inoculum, Family, CS_type_duo, CS_duo, CS1_mono, CS2_mono, ab1_corr, ab2_corr, ab_pred, ab_duo,rep_mono1,rep_mono2,rep_duo)

#.. pre-processing: exclude families that are only present in one of the inocula (first get n inocula per family per CS, then subset)
dt<-unique(select(dsp, Family, CS_duo,Inoculum))
dt<-  dt %>%
  dplyr::group_by(Family,CS_duo) %>%
  dplyr::summarise(nt_inoc=n())
dt<-as.data.frame(dt)
dsp.s<-merge(dsp, dt)
dsp.s<-subset(dsp.s, nt_inoc==2)
dsp.s<-select(dsp.s,-nt_inoc)

#.. Test whether Family X is significantly MORE abundant than predicted by the null model (for each CS pair and inoculum)
dsp.s$rep_combo_inoc<-paste0(dsp.s$rep_mono1,dsp.s$rep_mono2, dsp.s$rep_duo,"_", dsp.s$Inoculum)
#. for all CS pairs except glycine (see script below for glycine)
ds<-subset(dsp.s, CS2_mono!='glycine')
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
  
  my.t.test.p.value <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
  } 
  
  my.t.test.mean <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$estimate) 
  } 
  
  my.t.test.t.stat <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$statistic) 
  } 
  
  my.t.test.altern <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$alternative) 
  }
  
  my.t.test.method <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$method) 
  }
  
  t.test = pairing %>%
    dplyr::group_by(Family, CS_duo) %>%
    dplyr::mutate( p.val=my.t.test.p.value(ab_duo, ab_pred),
                   mean_diff=my.t.test.mean(ab_duo, ab_pred),
                   t.stat=my.t.test.t.stat(ab_duo, ab_pred),
                   t.method=my.t.test.method(ab_duo, ab_pred),
                   t.altern=my.t.test.altern(ab_duo, ab_pred),
                   n=n(),
                   epsi_sd=sd(ab_duo- ab_pred),
                   epsi_se=epsi_sd/sqrt(n)
    )
  
  t.test<-as.data.frame(t.test)
  t.test<-select(t.test, Family, CS_type_duo, CS_duo, CS1_mono, CS2_mono, mean_diff, p.val, n, t.stat, t.altern, t.method,epsi_sd,epsi_se)
  t.test<-unique(t.test)
  
  #.index permutation number
  t.test$perm.repn<-i
  datalist[[i]] <- as.data.frame(t.test)
}

pairing.bs<-do.call("rbind", datalist)

dp1<-pairing.bs

#  ---------------------------------------------------------------------------
##.. now do test for GLYCINE only (this is because there are 3 replicates not 4)
#. create 6 unique singleCS-singleCS- pairedCS combinations 
ds<-subset(dsp.s, CS2_mono=='glycine')
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
  
  #.. paired t-test asking if means difference is different from 0
  
  my.t.test.p.value <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
  } 
  
  my.t.test.mean <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$estimate) 
  } 
  
  my.t.test.t.stat <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$statistic) 
  } 
  
  my.t.test.altern <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$alternative) 
  }
  
  my.t.test.method <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$method) 
  }
  
  t.test = pairing %>%
    dplyr::group_by(Family, CS_duo) %>%  
    dplyr::mutate( p.val=my.t.test.p.value(ab_duo, ab_pred),
                   mean_diff=my.t.test.mean(ab_duo, ab_pred),
                   t.stat=my.t.test.t.stat(ab_duo, ab_pred),
                   t.method=my.t.test.method(ab_duo, ab_pred),
                   t.altern=my.t.test.altern(ab_duo, ab_pred),
                   n=n(),
                   epsi_sd=sd(ab_duo- ab_pred),
                   epsi_se=epsi_sd/sqrt(n))
          
  t.test<-as.data.frame(t.test)
  t.test<-select(t.test, Family, CS_type_duo, CS_duo, CS1_mono, CS2_mono, mean_diff, p.val, n, t.stat, t.altern, t.method,epsi_sd,epsi_se)
  t.test<-unique(t.test)

  #.index permutation number
  t.test$perm.repn<-i
  datalist[[i]] <- as.data.frame(t.test)
}

pairing.bs<-do.call("rbind", datalist)
dp2<-pairing.bs

#.. merge glycine and all other CS data
dp.all.g<-rbind(dp1,dp2)

#  ---------------------------------------------------------------------------
##.. now do test for whether Family X is significantly LESS abundant than predicted by the null model (for each CS pair and inoculum)
t.type<-"less" #.. select alternative type ('greater' or 'less')

#. for all CS pairs except glycine (see script below for glycine)
ds<-subset(dsp.s, CS2_mono!='glycine')

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
  
  #.. paired t-test asking if means difference is different from 0
  #t.test(pairing$ab_duo, pairing$ab_pred, paired = paired.v, alternative = t.type)
  
  my.t.test.p.value <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
  } 
  
  my.t.test.mean <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$estimate) 
  } 
  
  my.t.test.t.stat <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$statistic) 
  } 
  
  my.t.test.altern <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$alternative) 
  }
  
  my.t.test.method <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$method) 
  }
  
  t.test = pairing %>%
    dplyr::group_by(Family, CS_duo) %>%
    dplyr::mutate( p.val=my.t.test.p.value(ab_duo, ab_pred),
                   mean_diff=my.t.test.mean(ab_duo, ab_pred),
                   t.stat=my.t.test.t.stat(ab_duo, ab_pred),
                   t.method=my.t.test.method(ab_duo, ab_pred),
                   t.altern=my.t.test.altern(ab_duo, ab_pred),
                   n=n(),
                   epsi_sd=sd(ab_duo- ab_pred),
                   epsi_se=epsi_sd/sqrt(n)
    )
  
  t.test<-as.data.frame(t.test)
  t.test<-select(t.test, Family, CS_type_duo, CS_duo, CS1_mono, CS2_mono, mean_diff, p.val, n, t.stat, t.altern, t.method,epsi_sd,epsi_se)
  t.test<-unique(t.test)
  
  #.index permutation number
  t.test$perm.repn<-i
  datalist[[i]] <- as.data.frame(t.test)
}

pairing.bs<-do.call("rbind", datalist)
dp1<-pairing.bs

#  ---------------------------------------------------------------------------
##.. now test for GLYCINE only 
#. create 6 unique singleCS-singleCS- pairedCS combinations 

ds<-subset(dsp.s, CS2_mono=='glycine')
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
  
  #.. paired t-test asking if means difference is different from 0
  
  my.t.test.p.value <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
  } 
  
  my.t.test.mean <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$estimate) 
  } 
  
  my.t.test.t.stat <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$statistic) 
  } 
  
  my.t.test.altern <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$alternative) 
  }
  
  my.t.test.method <- function(...) { 
    obj<-try(t.test(..., paired = paired.v, alternative = t.type), silent=TRUE) 
    if (is(obj, "try-error")) return(NA) else return(obj$method) 
  }
  
  t.test = pairing %>%
    dplyr::group_by(Family, CS_duo) %>%
    dplyr::mutate( p.val=my.t.test.p.value(ab_duo, ab_pred),
                   mean_diff=my.t.test.mean(ab_duo, ab_pred),
                   t.stat=my.t.test.t.stat(ab_duo, ab_pred),
                   t.method=my.t.test.method(ab_duo, ab_pred),
                   t.altern=my.t.test.altern(ab_duo, ab_pred),
                   n=n(),
                   epsi_sd=sd(ab_duo- ab_pred),
                   epsi_se=epsi_sd/sqrt(n)
    )
  
  t.test<-as.data.frame(t.test)
  t.test<-select(t.test, Family, CS_type_duo, CS_duo, CS1_mono, CS2_mono, mean_diff, p.val, n, t.stat, t.altern, t.method,epsi_sd,epsi_se)
  t.test<-unique(t.test)
  
  #.index permutation number
  t.test$perm.repn<-i
  datalist[[i]] <- as.data.frame(t.test)
}

pairing.bs<-do.call("rbind", datalist)
dp2<-pairing.bs

#.. merge glycine and all other CS data
dp.all.l<-rbind(dp1,dp2)

#.. combine output of the two t-tests
dp.all.lg<-rbind(dp.all.l,dp.all.g)

#  ---------------------------------------------------------------------------
#.. Combine the two one-sided paired t-tests (greater and less) 
#.. and plot epsilon for Family x CS pair
#.. NOTE: make sure to run both 'greater' and 'less' before running script below
#------

#.. only select pairs with n>6
d.lg<-subset(dp.all.lg, n>=6)
hist(d.lg$n)

#.. assign significance value for permutation test (depends on df: https://www.sjsu.edu/faculty/gerstman/StatPrimer/t-table.pdf)
sig.v<-0.95
t.c8<-1.895 # t-test t value value assuming that df=7 (i.e. N=8-1) and cf=0.95
t.c7<-1.943 # t-test t value value assuming that df=7 (i.e. N=7-1) and cf=0.95
t.c6<-2.015 # t-test t value value assuming that df=5 (i.e. N=6-1) and cf=0.95

d.lg$t.critic <- with(d.lg, ifelse(n ==8, t.c8,
                                   ifelse(n==7, t.c7,
                                          ifelse(n==6, t.c6, NA))))
#..
d.lg$t.score <- with(d.lg, ifelse(t.altern=='greater' & t.stat>t.critic, 'PASS',
                                  ifelse(t.altern=='greater' & t.stat<=t.critic, 'FAIL',
                                         ifelse(t.altern=='less' & t.stat<(-t.critic), 'PASS',
                                                ifelse(t.altern=='less' & t.stat>= (-t.critic), 'FAIL',NA)))))

d.lg<- d.lg %>% 
  dplyr::group_by(Family, CS_duo, CS_type_duo,CS1_mono,CS2_mono,t.altern) %>% 
  dplyr::summarize(
    npass=table(t.score)['PASS'],
    nfail=table(t.score)['FAIL'],
    epsi_mean=mean(mean_diff),
    epsi_sd=mean(epsi_sd))

d.lg<-as.data.frame(d.lg)
d.lg[is.na(d.lg)]<-0
d.lg$tt.successr<-d.lg$npass/bs.rep
d.lg<-select(d.lg, -npass,-nfail)

#..confirm significance by checking both sides of t-test
d.stat<-spread(d.lg, t.altern,tt.successr)

d.stat$greater[d.stat$greater>sig.v]<-'signif'
d.stat$greater[d.stat$greater<=sig.v]<-'ns'
d.stat$less[d.stat$less>sig.v]<-'signif'
d.stat$less[d.stat$less<=sig.v]<-'ns'

d.stat$t.outcome[d.stat$greater=='ns' & d.stat$less=='ns']<-'ns'
d.stat$t.outcome[d.stat$greater=='signif' | d.stat$less=='signif']<-'signif'
d.stat$t.outcome[d.stat$greater=='signif' & d.stat$less=='signif']<-'NA'

d.stat<-select(d.stat,-greater,-less)

#.. write data as csv file (outcome of all paired t-tests for each family in each CS pair)
dstat<-select(d.stat, -cs.ord)
#fname<-paste0("ttest_paired_outcome_",taxa.level, "_all.csv")
#write.csv(dstat, file.path(path.t, fname),row.names=FALSE)

