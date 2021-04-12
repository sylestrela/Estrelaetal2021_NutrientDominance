" R script for Figure 4- figure supplement S2 in Estrela et al (2021)
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

seed_1_family = family_pairs[q == 0.9 & q2 == 0.9]
seed_1_esv = esv_pairs[q == 0.9 & q2 == 0.9]
#Remove runs that did not converge in at least one relevant environment
error_table = error_table[exp_id %in% seed_1_esv$exp_id,]
error_table = error_table[SS>1e-5 | SA >1e-5 | AA >1e-5]
seed_1_family =seed_1_family[exp_id %!in% error_table$exp_id]
seed_1_esv =seed_1_esv[exp_id %!in% error_table$exp_id]

seed_1_family$S = paste('n = ',seed_1_family$S)
seed_1_esv$S = paste('n = ',seed_1_esv$S)

p1 <- ggplot(seed_1_family,aes(x=Predicted_Abundance,y=Abundance,col=Treatment)) +
  geom_point(alpha=0.9,size=2) + theme_minimal() +  
  labs(x='predicted abundance', y='observed abundance') +
  scale_colour_manual(values = c( 'tan1', 'gray', 'mediumpurple1'))+
  scale_fill_manual(values = c( 'tan1', 'gray', 'mediumpurple1'))  +
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
  scale_x_continuous(limits=c(0,1),breaks=c(0,1)) + ggtitle('Family-level') +facet_wrap(~S,ncol=1)+coord_fixed()
p2 <- ggplot(seed_1_esv,aes(x=Predicted_Abundance,y=Abundance,col=Treatment)) +
  geom_point(alpha=0.9,size=2) + theme_minimal() +  
  labs(x='predicted abundance', y='observed abundance') +
  scale_colour_manual(values = c( 'tan1', 'gray', 'mediumpurple1'))+
  scale_fill_manual(values = c( 'tan1', 'gray', 'mediumpurple1'))  +
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
  scale_x_continuous(limits=c(0,1),breaks=c(0,1)) + ggtitle('Species-level') +facet_wrap(~S,ncol=1)  +coord_fixed()

p1cd<-ggarrange(p2,p1, ncol=2, nrow=1, labels=c('A','B'), common.legend=TRUE)
p1cd

ggsave(path=path.p, filename=paste0("fig4S2.pdf"), width=6, height=8)  

#.. write data figure as csv file
seed_1_esv$level<-'ESV'
seed_1_family$level<-'Family'
seed_1_family$ESV<-NA
ds<-rbind(seed_1_esv, seed_1_family)
ds<-select(ds,-Abundance_1,-Abundance_2,-Replicate,-variable)
fname<-paste("fig4S2_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)

