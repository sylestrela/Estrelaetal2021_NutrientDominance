" R script for Figure 4- figure supplement S5 in Estrela et al (2021)
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
library(stringr)

#.. Assign path to data
#. usern = path to user directory
path.d<-paste0(usern, "/Estrelaetal_2021/data")
#.. Assign path to plots output
path.p<- paste0(usern,"/Estrelaetal_2021/Rplots")
#.. Assign path to data figures
path.ds<- paste0(usern,"/Estrelaetal_2021/data/data_figures/")

#.. Read CRM simulations data
path.crm <-paste0(path.d,"/data_FBA_Simulations")
setwd(path.crm)

dat<-read.csv('FBA_oxygen.csv')

#. define order of CS names to plot
dat$CS <- factor(dat$CS, levels = c("glucose", "fructose", "cellobiose", "ribose",  "glycerol",
                                    'succinate', 'benzoate','fumarate',"glycine", 'glutamine'))

ggplot(dat, aes(x=CS, y=Oxygen))+
  geom_bar(stat="identity")+
  labs(y=expression('O'[2]*'/C'),x='')+
  theme_classic()+
  scale_y_continuous(breaks=c(0,0.5, 1), limit=c(0,1))+
  theme(axis.text.x = element_text(size=13, angle=45, hjust=1),
        axis.text.y = element_text(size=12),
        axis.title.y =element_text(size=18),
        axis.title.x = element_text(size = 18),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
ggsave(path=path.p, "fig4S5.pdf", width=8, height=3.5,onefile = TRUE)

#.. write data figure as csv file
ds<-select(dat,-Biomass,-Carbon_Source)
fname<-paste("fig4S5_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)






