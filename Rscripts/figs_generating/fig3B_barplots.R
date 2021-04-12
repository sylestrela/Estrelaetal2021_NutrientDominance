" R script for Figure 3B (barplots) in Estrela et al (2021) 
Nutrient dominance governs the assembly of microbial communities in mixed nutrient environments.
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

#.. Plot relative abundance
#. sum relative abundance of each family for each community
dp1<-select(dat,-Sequence,-Order,-Transfer,-ESV,-Genus)

sum_relab_fam<- dp1 %>%
  dplyr::group_by(Inoculum,Carbon_Source,Replicate,s_perc,cs_type,Experiment,nCS,Family,CS1,CS2) %>%
  dplyr::summarise(
    sum_relab=sum(Relative_Abundance))
dp1<-as.data.frame(sum_relab_fam)

#.. Select minimum abundance- taxa below threshold are 'other'
min.ab<-0.01
dfp.rare<-subset(dp1, sum_relab<min.ab)
dfp.rare$Family_type<-'Other'
dfp.com<-subset(dp1, sum_relab>=min.ab)
dfp.com$Family_type<-dfp.com$Family
dfp<-rbind(dfp.com,dfp.rare)

#.. Colour scheme for plots
#. Family colors 
colEntbf<-'#516091'
colMoraxf<-'#74BEC1'
colPseuf<-'#E2F0CB'
colRhizf<-'#996666'
colOther<-'gray85'

##.. Plot Figure 3B: barplots of succinate + fructose (representative examples)
dp1<-subset(dfp, Replicate==2 & Inoculum==2 & Carbon_Source %in% c("succinate_fructose",'succinate', 'fructose'))
#. define order of duo CS names to plot
dp1$cs.ord <- factor(dp1$Carbon_Source, levels = c('succinate', 'fructose', "succinate_fructose"))
p1<-ggplot(dp1, aes(x=cs.ord, y=sum_relab, fill=Family_type))+
  geom_bar(color="black", stat="identity", position="stack") + 
  labs(y= "Relative abundance", x="")+
  scale_fill_manual(values=c(colEntbf, colMoraxf, colPseuf, colRhizf, '#FECBA5', '#8EB8B6', '#FECB89','gray85'))+ 
  scale_y_continuous( breaks=seq(0,1,1)) +
  guides(fill=guide_legend(title="Family"))+
  theme_minimal()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=14),
        axis.title=element_text(size=16),
        strip.text.y = element_text(size = 11),
        strip.text.x = element_text(size = 8),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  theme(panel.spacing = unit(0, "lines"))
p1
ggsave(path=path.p, paste0("fig3B_barplots_succ_fruc.pdf"), width=6, height=4.5, onefile = TRUE)

#.. write data figure as csv file
ds<-select(dp1, Inoculum, Carbon_Source, Family, Replicate, sum_relab, Family_type, cs.ord)
fname<-paste("fig3B_barplot_succ_fruc_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)


##.. Figure 3B: barplots of plot glucose + glutamine (representative examples)
dp2_1<-subset(dfp, Carbon_Source %in% c('glucose', 'glutamine') & Inoculum==2 & Replicate==3)
dp2_2<-subset(dfp, Carbon_Source %in% c('glucose_glutamine') & Inoculum==2 & Replicate==2)
dp2<-rbind(dp2_1,dp2_2)
#. define order of duo CS names to plot
dp2$cs.ord <- factor(dp2$Carbon_Source, levels = c( 'glutamine','glucose', 'glucose_glutamine'))
dp2$fam.ord <- factor(dp2$Family_type, levels = c( 'Enterobacteriaceae','Moraxellaceae', 'Pseudomonadaceae','Rhizobiaceae','Comamonadaceae','Other'))

p2<-ggplot(dp2, aes(x=cs.ord, y=sum_relab, fill=fam.ord))+
  geom_bar(color="black", stat="identity", position="stack") + 
  labs(y= "Relative abundance", x="")+
  scale_fill_manual(values=c(colEntbf, colMoraxf, colPseuf,  colRhizf, colOther,colOther,colOther, '#8EB8B6', '#FECB89','gray85',colOther))+ 
  scale_y_continuous( breaks=seq(0,1,1)) +
  guides(fill=guide_legend(title="Family"))+
  theme_minimal()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=14),
        axis.title=element_text(size=16),
        strip.text.y = element_text(size = 11),
        strip.text.x = element_text(size = 8),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  theme(panel.spacing = unit(0, "lines"))
p2
ggsave(path=path.p, paste0("fig3B_barplots_glu_glut.pdf"), width=6, height=4.5, onefile = TRUE)

#.. write data figure as csv file
ds<-select(dp2, Inoculum, Carbon_Source, Family, Replicate, sum_relab, Family_type, cs.ord,fam.ord)
fname<-paste("fig3B_barplot_glu_glut_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)


