" R script for Figure 4C in Estrela et al (2021)
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
dat<-read.csv('data_processed/ravg_isolates.csv')

#.. get mean for each family-CS pair
summ<- dat %>%
  dplyr::group_by(Family,CS) %>%
  dplyr::summarise(
    ravg.mean = mean(ravg))
d.m<-as.data.frame(summ)
head(d.m)

#.. assign group to CS
d.m$CS_type[d.m$CS %in% c('glucose', 'cellobiose', 'fructose', 'ribose', 'glycerol')]<-'sugar'
d.m$CS_type[d.m$CS %in% c('succinate','fumarate','benzoate','glutamine','glycine')]<-'acid'

d.ent<-subset(d.m, Family=='Enterobacteriaceae')
d.other<-subset(d.m, Family!='Enterobacteriaceae')

d.ent$FamE<-d.ent$Family
d.ent$ravg.meanE<-d.ent$ravg.mean
d.ent<-select(d.ent,-Family,-ravg.mean)

d.other$FamO<-d.other$Family
d.other$ravg.meanO<-d.other$ravg.mean
d.other<-select(d.other,-Family,-ravg.mean)

dm<-merge(d.ent,d.other,by=c('CS','CS_type'))
dm$CS <- factor(dm$CS, levels = c("glucose", "fructose", 'cellobiose', 'ribose','glycerol', 'succinate', "fumarate", "benzoate", 'glutamine', 'glycine'))

#.. Calculate d_ravg(S) and d_ravg(A)
d.sug<-subset(dm, CS_type=='sugar')
d.oa<-subset(dm, CS_type=='acid')

d.sug$dravg.mean_sug<-d.sug$ravg.meanE-d.sug$ravg.meanO
d.sug$sugar<-d.sug$CS
d.sug<-select(d.sug,-CS, -CS_type,-ravg.meanE,-ravg.meanO)

d.oa$dravg.mean_oa<-d.oa$ravg.meanO-d.oa$ravg.meanE
d.oa$acid<-d.oa$CS
d.oa<-select(d.oa,-CS, -CS_type,-ravg.meanE,-ravg.meanO)

dm2<-merge(d.sug,d.oa, by=c('FamE', 'FamO'))

dm2$CSpair<-paste0(dm2$sugar,'_',dm2$acid)

##.. plot 
#. Family colors 
colEntbf<-'#516091'
colMoraxf<-'#74BEC1'
colPseuf<-'#E2F0CB'
colRhizf<-'#996666'

yl<-expression('q'[A])
xl<-expression('q'[S])

dp<-subset(dm2, (sugar=='glucose' & acid!='succinate') | (acid=='succinate' & sugar!='glucose') )
ggplot(dp, aes(x=dravg.mean_sug, y=dravg.mean_oa, group=factor(FamO), colour=factor(FamO) ))+
  geom_point(size=3, alpha=0.8)+
  geom_point(size=3, shape=1)+
  theme_minimal()+
  geom_abline(colour='gray', linetype='dashed')+
  geom_hline(yintercept=0, colour='gray')+
  geom_vline(xintercept=0, colour='gray')+
  scale_colour_manual(values=c( colMoraxf, colPseuf, colRhizf))+
  scale_fill_manual(values=c( colMoraxf, colPseuf, colRhizf))+
  theme(
    axis.text.x = element_text(size=18),
    axis.text.y = element_text(size=18),
    axis.title.y = element_text(size=20),
    axis.title.x =element_text(size=20),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())+
  labs(x=xl, y=yl)+
   theme(legend.title = element_blank(),
        legend.text=element_text(size=16))

ggsave(path=path.p, paste0("fig4C.pdf"), width=7.3, height=4.5, onefile = TRUE)

#.. write data figure as csv file
ds<-dp
fname<-paste("fig4C_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)

##.. calculate mean Qs and Qa
dp3<-select(dp, CSpair, FamE, FamO, dravg.mean_oa, dravg.mean_sug)
summ<- dp3 %>%
  dplyr::summarise(
    qs.mean = mean(dravg.mean_sug),
    qs.sd=sd(dravg.mean_sug),
    qa.mean = mean(dravg.mean_oa),
    qa.sd=sd((dravg.mean_oa)))
q.mean<-as.data.frame(summ)
head(q.mean)

#.. stats
t.test(dp3$dravg.mean_sug, dp3$dravg.mean_oa, paired=TRUE, alternative = "greater")




