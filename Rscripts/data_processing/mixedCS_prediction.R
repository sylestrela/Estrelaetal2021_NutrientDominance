" R script for Estrela et al (2021)
Nutrient dominance governs the assembly of microbial communities in mixed nutrient environments.

1) calculate Frequencies in single CS weighted by OD;
2) calculate the predicted relative abundance in the pairs from abundance in singles 
assuming no interaction between CS (null additive model). 

Output: 3 csv files with observed and predicted data for each community at family/genus/esv level.
"

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
dat<-read.csv('data_master/data_16s_OD_singleCS_duoCS_processed.csv')

## ------------------------------
#.. FAMILY-LEVEL PREDICTION
##. 1) sum relative abundance of each family for each community
dp1<-select(dat,-Sequence,-Order,-Transfer,-ESV,-Genus)
sum_relab_fam<- dp1 %>%
  dplyr::group_by(Inoculum,Carbon_Source,Replicate,s_perc,cs_type,Experiment,nCS,Family,CS1,CS2) %>%
  dplyr::mutate(
    sum_relab=sum(Relative_Abundance))
dp2<-unique(select(as.data.frame(sum_relab_fam),-Relative_Abundance))

#. create exp_id based on inoculum, replicate and carbon source 
dp2$exp_id<-paste0(dp2$Inoculum,"_", dp2$Replicate,"_",dp2$Carbon_Source)
test<-select(dp2,exp_id, Family,sum_relab)
testd<-dcast(test, exp_id~Family, value.var='sum_relab')

#. if abundance is NA then add 0 (as it means that taxa is absent)
testd[is.na(testd)] <- 0
testm<-melt(testd)
names(testm)[names(testm) == "variable"] <- "Family"
names(testm)[names(testm) == "value"] <- "sum_relab"

#. merge 16s abundance and od data back together
test2<-unique(select(dp2,-sum_relab,-Family))
dp2m<-merge(test2, testm, by='exp_id')
dp2<-select(dp2m, -exp_id)

#. check it works (e.g. Moraxellaceae is absent in glycerol)
subset(dp2, Family=='Moraxellaceae' & Carbon_Source=='glycerol')

##. 2) get observed abundance in single CS
dp.s<-subset(dp2, Experiment=='CSS_Single')
dp.s<-select(dp.s,-CS2,-CS1)
d1<-dp.s
names(d1)[names(d1) == "Carbon_Source"] <- "CS_mono"
names(d1)[names(d1) == "Replicate"] <- "rep_mono"
names(d1)[names(d1) == "sum_relab"] <- "ab_mean_mono"
names(d1)[names(d1) == "od.mean"] <- "od_mean_mono"
d.obs.sm<-select(d1,-nCS,-cs_type,-s_perc,-Experiment)

##. 3) get observed abundance in mixed CS (PAIRS)
#. select pairs of CS
dp.d<-subset(dp2, nCS==2)
d2<-select(dp.d,- Experiment,-nCS,-s_perc)
names(d2)[names(d2) == "Carbon_Source"] <- "CS_duo"
names(d2)[names(d2) == "Replicate"] <- "rep_duo"
names(d2)[names(d2) == "sum_relab"] <- "ab_duo"
names(d2)[names(d2) == "od.mean"] <- "od_duo"
names(d2)[names(d2) == "cs_type"] <- "CS_type_duo"
d.obs.d<-d2

##. 4) merge observed single and observed mixture abundance
#. first need to separate cs1 and cs2 to add mono, and then merge
d.cs1<-subset(d.obs.sm, CS_mono %in% c('glucose', 'succinate'))
names(d.cs1)[names(d.cs1) == "CS_mono"] <- "CS1_mono"
names(d.cs1)[names(d.cs1) == "ab_mean_mono"] <- "ab1_mean_mono"
names(d.cs1)[names(d.cs1) == "od_mean_mono"] <- "od1_mean_mono"
names(d.cs1)[names(d.cs1) == "rep_mono"] <- "rep_mono1"

d.cs2<-subset(d.obs.sm, !(CS_mono %in% c('glucose', 'succinate')))
names(d.cs2)[names(d.cs2) == "CS_mono"] <- "CS2_mono"
names(d.cs2)[names(d.cs2) == "ab_mean_mono"] <- "ab2_mean_mono"
names(d.cs2)[names(d.cs2) == "od_mean_mono"] <- "od2_mean_mono"
names(d.cs2)[names(d.cs2) == "rep_mono"] <- "rep_mono2"

d.obs.all1<-left_join(d.cs1,d.obs.d, by=c('Inoculum','Family'))
d.obs.all2<-left_join(d.cs2, d.obs.all1, by=c('Inoculum','Family'))

d.obs.all3<-d.obs.all2[as.character(d.obs.all2$CS1)==as.character(d.obs.all2$CS1_mono),]
d.obs.all<-d.obs.all3[as.character(d.obs.all3$CS2)==as.character(d.obs.all3$CS2_mono),]

##. 5) calculate abundance corrected with OD for cs1 and cs2 
d.obs.all$ab1_corr<-d.obs.all$ab1_mean_mono*d.obs.all$od1_mean_mono/(d.obs.all$od1_mean_mono+d.obs.all$od2_mean_mono)
d.obs.all$ab2_corr<-d.obs.all$ab2_mean_mono*d.obs.all$od2_mean_mono/(d.obs.all$od1_mean_mono+d.obs.all$od2_mean_mono)

##. 6) calculate predicted abundance in mixture
d.obs.all$ab_pred<-(d.obs.all$ab1_corr +d.obs.all$ab2_corr)

#..remove taxa that have mono_obs, duo_obs all=0 for Inoculum x CS_duo (that is, we only want to keep taxa that are present in at least one sample for each inoc:CS to avoid artifacts)
d.all<- d.obs.all %>% 
  group_by(Inoculum, Family, CS_duo) %>%
  filter(ab1_mean_mono!=0 | ab2_mean_mono!=0 | ab_duo!=0)

## ..Write data table as csv file
#fname<-paste("data_mixedCS_obs_pred_fam_duos.csv")
#write.csv(as.data.frame(d.all),file.path(path.d, fname),row.names=FALSE)

## ------------------------------
#.. GENUS-LEVEL PREDICTION
#.. 1) sum relative abundance of each genus for each community
dp1<-select(dat,-Sequence,-Order,-Transfer,-ESV,-Family)
sum_relab_gen<- dp1 %>%
  dplyr::group_by(Inoculum,Carbon_Source,Replicate,s_perc,cs_type,Experiment,nCS,Genus,CS1,CS2) %>%
  dplyr::mutate(
    sum_relab=sum(Relative_Abundance))
dp2<-unique(select(as.data.frame(sum_relab_gen),-Relative_Abundance))

#. create exp_id based on inoculum, replicate and carbon source 
dp2$exp_id<-paste0(dp2$Inoculum,"_", dp2$Replicate,"_CS",dp2$Carbon_Source)
test<-select(dp2,exp_id, Genus,sum_relab)
testd<-dcast(test, exp_id~Genus, value.var='sum_relab')

#. if abundance is NA then add 0 (as it means that taxa is absent)
testd[is.na(testd)] <- 0
testm<-melt(testd)
names(testm)[names(testm) == "variable"] <- "Genus"
names(testm)[names(testm) == "value"] <- "sum_relab"

#. merge 16s abundance and od data back together
test2<-unique(select(dp2,-sum_relab,-Genus))
dp2m<-merge(test2, testm,by='exp_id')
dp2<-select(dp2m, -exp_id)

##.. 2) get observed abundance in single CS
dp.s<-subset(dp2, Experiment=='CSS_Single')
dp.s<-select(dp.s,-CS2,-CS1)
d1<-dp.s
names(d1)[names(d1) == "Carbon_Source"] <- "CS_mono"
names(d1)[names(d1) == "Replicate"] <- "rep_mono"
names(d1)[names(d1) == "sum_relab"] <- "ab_mean_mono"
names(d1)[names(d1) == "od.mean"] <- "od_mean_mono"
d.obs.sm<-select(d1,-nCS,-cs_type,-s_perc,-Experiment)

#.. 3) get observed abundance in mixed CS (PAIRS)
#.. select pairs of CS
dp.d<-subset(dp2, nCS==2)
d2<-select(dp.d,- Experiment,-nCS,-s_perc)
names(d2)[names(d2) == "Carbon_Source"] <- "CS_duo"
names(d2)[names(d2) == "Replicate"] <- "rep_duo"
names(d2)[names(d2) == "sum_relab"] <- "ab_duo"
names(d2)[names(d2) == "od.mean"] <- "od_duo"
names(d2)[names(d2) == "cs_type"] <- "CS_type_duo"
d.obs.d<-d2

##. 4) merge observed single and observed mixture abundance
#. first need to separate cs1 and cs2 to add mono, and then merge
d.cs1<-subset(d.obs.sm, CS_mono %in% c('glucose', 'succinate'))
names(d.cs1)[names(d.cs1) == "CS_mono"] <- "CS1_mono"
names(d.cs1)[names(d.cs1) == "ab_mean_mono"] <- "ab1_mean_mono"
names(d.cs1)[names(d.cs1) == "od_mean_mono"] <- "od1_mean_mono"
names(d.cs1)[names(d.cs1) == "rep_mono"] <- "rep_mono1"

d.cs2<-subset(d.obs.sm, !(CS_mono %in% c('glucose', 'succinate')))
names(d.cs2)[names(d.cs2) == "CS_mono"] <- "CS2_mono"
names(d.cs2)[names(d.cs2) == "ab_mean_mono"] <- "ab2_mean_mono"
names(d.cs2)[names(d.cs2) == "od_mean_mono"] <- "od2_mean_mono"
names(d.cs2)[names(d.cs2) == "rep_mono"] <- "rep_mono2"

d.obs.all1<-left_join(d.cs1,d.obs.d, by=c('Inoculum','Genus'))
d.obs.all2<-left_join(d.cs2, d.obs.all1, by=c('Inoculum','Genus'))

d.obs.all3<-d.obs.all2[as.character(d.obs.all2$CS1)==as.character(d.obs.all2$CS1_mono),]
d.obs.all<-d.obs.all3[as.character(d.obs.all3$CS2)==as.character(d.obs.all3$CS2_mono),]

##. 5) calculate abundance corrected for cs1 and cs2 
d.obs.all$ab1_corr<-d.obs.all$ab1_mean_mono*d.obs.all$od1_mean_mono/(d.obs.all$od1_mean_mono+d.obs.all$od2_mean_mono)
d.obs.all$ab2_corr<-d.obs.all$ab2_mean_mono*d.obs.all$od2_mean_mono/(d.obs.all$od1_mean_mono+d.obs.all$od2_mean_mono)

##. 6) calculate predicted abundance 
d.obs.all$ab_pred<-(d.obs.all$ab1_corr +d.obs.all$ab2_corr)

#..remove taxa that have mono_obs, duo_obs all=0 for Inoculum x CS_duo (that is, we only want to keep taxa that are present in at least one sample for each inoc: CS to avoid artifacts)
d.all<- d.obs.all %>% 
  group_by(Inoculum, Genus, CS_duo) %>%
  filter(ab1_mean_mono!=0 | ab2_mean_mono!=0 | ab_duo!=0)

#.. match Genus to Family
d.tax<-select(dat, Genus, Family)
d.tax<-d.tax %>% distinct()
d.all<-merge(d.tax,d.all, by='Genus')

## ..Write data table as csv file
#fname<-paste("data_mixedCS_obs_pred_gen_duos.csv")
#write.csv(d.all,file.path(path.d, fname))

## ------------------------------
#.. ESV-LEVEL PREDICTION
#.. 1) sum relative abundance of each genus for each community
dp1<-select(dat,-Sequence,-Order,-Transfer,-Genus,-Family)
sum_relab_gen<- dp1 %>%
  dplyr::group_by(Inoculum,Carbon_Source,Replicate,s_perc,cs_type,Experiment,nCS,ESV,CS1,CS2) %>%
  dplyr::mutate(
    sum_relab=sum(Relative_Abundance))
dp2<-as.data.frame(sum_relab_gen)
dp2<-unique(select(dp2,-Relative_Abundance))

#. create exp_id based on inoculum, replicate and carbon source 
dp2$exp_id<-paste0(dp2$Inoculum,"_", dp2$Replicate,"_CS",dp2$Carbon_Source)
test<-select(dp2,exp_id, ESV,sum_relab)
testd<-dcast(test, exp_id~ESV, value.var='sum_relab')

#. if abundance is NA then add 0 (as it means that taxa is absent)
testd[is.na(testd)] <- 0
testm<-melt(testd)
names(testm)[names(testm) == "variable"] <- "ESV"
names(testm)[names(testm) == "value"] <- "sum_relab"

#. merge 16s abundance and od data back together
test2<-unique(select(dp2,-sum_relab,-ESV))
dp2m<-merge(test2, testm, by='exp_id')
dp2<-select(dp2m, -exp_id)

##.. 2) get absorved abundance in single CS
dp.s<-subset(dp2, Experiment=='CSS_Single')
d1<-select(dp.s,-CS2, -CS1)
names(d1)[names(d1) == "Carbon_Source"] <- "CS_mono"
names(d1)[names(d1) == "Replicate"] <- "rep_mono"
names(d1)[names(d1) == "sum_relab"] <- "ab_mean_mono"
names(d1)[names(d1) == "od.mean"] <- "od_mean_mono"
d.obs.sm<-select(d1,-nCS,-cs_type,-s_perc,-Experiment)

#.. 3) get observed abundance in mixed CS (PAIRS, not trios)
#.. select pairs of CS
dp.d<-subset(dp2, nCS==2)
d2<-select(dp.d,- Experiment,-nCS,-s_perc)
names(d2)[names(d2) == "Carbon_Source"] <- "CS_duo"
names(d2)[names(d2) == "Replicate"] <- "rep_duo"
names(d2)[names(d2) == "sum_relab"] <- "ab_duo"
names(d2)[names(d2) == "od.mean"] <- "od_duo"
names(d2)[names(d2) == "cs_type"] <- "CS_type_duo"
d.obs.d<-d2

##. 4) merge observed single and observed mixture abundance
#. first need to separate cs1 and cs2 to add mono, and then merge
d.cs1<-subset(d.obs.sm, CS_mono %in% c('glucose', 'succinate'))
names(d.cs1)[names(d.cs1) == "CS_mono"] <- "CS1_mono"
names(d.cs1)[names(d.cs1) == "ab_mean_mono"] <- "ab1_mean_mono"
names(d.cs1)[names(d.cs1) == "od_mean_mono"] <- "od1_mean_mono"
names(d.cs1)[names(d.cs1) == "rep_mono"] <- "rep_mono1"

d.cs2<-subset(d.obs.sm, !(CS_mono %in% c('glucose', 'succinate')))
names(d.cs2)[names(d.cs2) == "CS_mono"] <- "CS2_mono"
names(d.cs2)[names(d.cs2) == "ab_mean_mono"] <- "ab2_mean_mono"
names(d.cs2)[names(d.cs2) == "od_mean_mono"] <- "od2_mean_mono"
names(d.cs2)[names(d.cs2) == "rep_mono"] <- "rep_mono2"

d.obs.all1<-left_join(d.cs1,d.obs.d, by=c('Inoculum','ESV'))
d.obs.all2<-left_join(d.cs2, d.obs.all1, by=c('Inoculum','ESV'))

d.obs.all3<-d.obs.all2[as.character(d.obs.all2$CS1)==as.character(d.obs.all2$CS1_mono),]
d.obs.all<-d.obs.all3[as.character(d.obs.all3$CS2)==as.character(d.obs.all3$CS2_mono),]

##. 5) calculate abundance corrected for cs1 and cs2 
d.obs.all$ab1_corr<-d.obs.all$ab1_mean_mono*d.obs.all$od1_mean_mono/(d.obs.all$od1_mean_mono+d.obs.all$od2_mean_mono)
d.obs.all$ab2_corr<-d.obs.all$ab2_mean_mono*d.obs.all$od2_mean_mono/(d.obs.all$od1_mean_mono+d.obs.all$od2_mean_mono)

##. 6) calculate predicted abundance 
d.obs.all$ab_pred<-(d.obs.all$ab1_corr +d.obs.all$ab2_corr)

#..remove taxa that have mono_obs, duo_obs all=0 for Inoculum x CS_duo (that is, we only want to keep taxa that are present in at least one sample for each inoc: CS to avoid artifacts)
d.all<- d.obs.all %>% 
  group_by(Inoculum, ESV, CS_duo) %>%
  filter(ab1_mean_mono!=0 | ab2_mean_mono!=0 | ab_duo!=0)

#.. match ESV to Genus and Family
d.tax<-select(dat, ESV, Genus, Family)
d.tax<-d.tax %>% distinct()
d.all<-merge(d.tax,d.all, by='ESV')

##.. Write data table as csv file
#fname<-paste("data_mixedCS_obs_pred_esv_duos.csv")
#write.csv(d.all,file.path(path.d, fname),row.names=FALSE)



