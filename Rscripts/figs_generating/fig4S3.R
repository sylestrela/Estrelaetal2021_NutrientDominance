" R script for Figure 4- figure supplement S3 in Estrela et al (2021)
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

family = fread('Family_Pairs.csv')
family = family[q!=0.9]
error_table = fread('Error_table.csv')
error_table = error_table[exp_id %in% family$exp_id,]
error_table = error_table[SS>1e-5 | SA >1e-5 | AA >1e-5]
family = family[exp_id %!in% error_table$exp_id]

cor_epsilon = family[,cor(Abundance,Predicted_Abundance),by=c('q','q2','Treatment')]
rmse_epsilon = family[,rmse(Abundance,Predicted_Abundance),by=c('q','q2','Treatment')]
rmse_epsilon = rmse_epsilon[q==q2]

ggplot(rmse_epsilon, aes(x=q,y=V1,col=Treatment)) +
  geom_point(size=3,shape=5, stroke=2) +
  scale_colour_manual(values=c('tan3','Grey','mediumpurple3')) +
  labs(x = 'specialization level (q)',y = 'RMSE') +
  theme_classic()+
  scale_x_continuous(breaks=seq(from = 0.05, to = 0.95, by = 0.1))+
  theme(axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=12),
        axis.title.y =element_text(size=18),
        axis.title.x = element_text(size = 18),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.text = element_text(size=11),
        legend.title = element_blank())
ggsave(path=path.p, "fig4S3.pdf", width=6.6, height=3.7)

#.. write data figure as csv file
ds<-select(rmse_epsilon, -q2)
names(ds)[names(ds)=='V1']<-'RMSE'
fname<-paste("fig4S3_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)



