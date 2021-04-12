" R script for Figure 2- figure supplement 5 in Estrela et al (2021) 
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
dat<-read.csv('data_master/data_16s_OD_singleCS_duoCS_processed.csv')

d.od<-select(dat, Carbon_Source, Inoculum,Replicate,nCS,cs_type,od.mean)

#.. Select single carbon source
d.od<-unique(subset(d.od, nCS==1))

#. define order of CS names to plot
d.od$Carbon_Source <- factor(d.od$Carbon_Source, levels = c("glucose", "fructose", "cellobiose", "ribose",  "glycerol",
                                                            'succinate', 'benzoate','fumarate',"glycine", 'glutamine'))

ggplot(d.od, aes(x=Carbon_Source, y=od.mean,group=factor(Inoculum):factor(Replicate), colour=factor(Inoculum)))+
  geom_point(size=3, shape=1)+
  geom_point(size=3, shape=16, alpha=0.4)+
  ylim(0, max(d.od$od.mean)+0.01)+
  labs(x='', y='OD620')+
  theme_classic()+
  theme(
    axis.text.x = element_text(size=13, angle=45,hjust=1),
    axis.text.y = element_text(size=14),
    axis.title.y = element_text(size=18),
    axis.title.x =element_text(size=18),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linetype='dashed', size=0.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.text = element_text(size=13),
    legend.title = element_text(size=13),
    legend.position = 'top')+
  scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'), name='Inoculum')

ggsave(path=path.p, "fig2S5.pdf", width=8, height=4)

#.. write data figure as csv file
ds<-select(d.od,-nCS,-cs_type)
fname<-paste("fig2S5_data.csv")
write.csv(ds, file.path(path.ds, fname), row.names=FALSE)







