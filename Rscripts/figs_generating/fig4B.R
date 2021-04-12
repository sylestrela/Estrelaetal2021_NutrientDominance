" R script for Figure 4B in Estrela et al (2021)
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
library(Metrics)
library(operators)

#.. Assign path to data
#. usern = path to user directory
path.d<-paste0(usern, "/Estrelaetal_2021/data")
#.. Assign path to plots output
path.p<- paste0(usern,"/Estrelaetal_2021/Rplots")
#.. Assign path to data figures
path.ds<- paste0(usern,"/Estrelaetal_2021/data/data_figures/")

#.. Read CRM simulations data
path.crm <-paste0(path.d,"/data_MiCRM_Simulations")
setwd(path.crm)

family_pairs=fread('Family_Pairs.csv')
esv_pairs =fread('ESV_Pairs.csv')
error_table = fread('Error_table.csv')

family = family_pairs[q == 0.9 & q2 == 0.9 & S==200]
esv = esv_pairs[q == 0.9 & q2 == 0.9 & S==200]
#Remove runs that did not converge in at least one relevant environment
error_table = error_table[exp_id %in% esv$exp_id,]
error_table = error_table[SS>1e-5 | SA >1e-5 | AA >1e-5]
family =family[exp_id %!in% error_table$exp_id]
esv =esv[exp_id %!in% error_table$exp_id]

col.ss<-'mediumpurple1'
col.aa<-'tan1'
col.sa<-'gray'

p1 <- ggplot(family,aes(x=Predicted_Abundance,y=Abundance,col=Treatment)) +
  geom_point(alpha=0.9,size=2) + theme_minimal() +  
  labs(x='predicted abundance', y='observed abundance',col='') +
  scale_colour_manual(values = c( col.aa, col.sa, col.ss))+
  scale_fill_manual(values = c( col.aa, col.sa, col.ss))  +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14),
        axis.title.y =element_text(size=16),
        axis.title.x = element_text(size = 16),
        panel.border = element_rect(color = "gray", fill = NA, size = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14)) + geom_abline(linetype='dashed') +
  scale_y_continuous(limits = c(0,1),breaks=c(0,1)) + 
  scale_x_continuous(limits=c(0,1),breaks=c(0,1)) + ggtitle('Family-level')+
  theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(esv,aes(x=Predicted_Abundance,y=Abundance,col=Treatment)) +
  geom_point(alpha=0.9,size=2) + theme_minimal() +  
  labs(x='predicted abundance', y='observed abundance',col='') +
  scale_colour_manual(values = c(  col.aa, col.sa, col.ss))+
  scale_fill_manual(values = c(  col.aa, col.sa, col.ss))  +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14),
        axis.title.y =element_text(size=16),
        axis.title.x = element_text(size = 16),
        panel.border = element_rect(color = "gray", fill = NA, size = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14)) + 
  geom_abline(linetype='dashed') +
  scale_y_continuous(limits = c(0,1),breaks=c(0,1)) + 
  scale_x_continuous(limits=c(0,1),breaks=c(0,1)) + ggtitle('Species-level')+
  theme(plot.title = element_text(hjust = 0.5))

p1cd<-ggarrange(p2,p1, ncol=2, nrow=1, common.legend=TRUE)
p1cd
ggsave(path=path.p, "fig4B.pdf", width=7.1, height=4, onefile = TRUE)

#.. write data figure as csv file
esv$level<-'ESV'
family$level<-'family'
family$ESV<-NA
ds<-rbind(esv,family)
ds<-select(ds,-variable,-Replicate,-Abundance_1,-Abundance_2,-S)
fname<-paste("fig4B_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)





