" R script for Figure 4D in Estrela et al (2021)
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

dominance = fread('SA_Dominance_Family.csv')
error_table = fread('Error_table.csv')
dominance = dominance[S==200 & Family=='F1' & q !=0.9]
error_table = error_table[exp_id %in% dominance$exp_id,]
#Remove runs that did not converge in at least one relevant environment
dominance = dominance[exp_id %in% error_table[SA<1e-5,]$exp_id]
mean_dominance = dominance[,mean(Epsilon),by=c('q','q2')]
n = dominance[,length(Epsilon),by=c('q','q2')]

p1 <- ggplot(mean_dominance,aes(x=q,y=q2,width=0.1,fill= V1)) + geom_tile() + 
  scale_fill_gradient2(low="darkblue",mid='white', high="tan2",midpoint=0,limits=c(-0.1,0.1),breaks=c(-0.1,0,0.1))+ theme_pubr() +
  labs(x = expression(q[S]),y= expression(q[A]),fill = expression(mean(delta)))  +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  theme(legend.position = 'right')+
  theme_minimal() + 
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
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
        legend.title=element_text(size=14)) 
p1
ggsave(path=path.p, "fig4D.pdf", width=5, height=3.7, onefile = TRUE)

#.. write data figure as csv file
ds<-mean_dominance
setnames(ds, old = c('q','q2', 'V1'), new = c('qs','qa','delta'))
fname<-paste("fig4D_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)


