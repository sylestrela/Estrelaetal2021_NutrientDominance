" R script for Figure 4- figure supplement S4 in Estrela et al (2021)
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

df = fread('C_matrices.csv')
df$Treatment = factor(df$Treatment,
                      levels=c('Unspecialised','Symmetric','Acid','Sugar'),
                      labels=c(expression(atop('Unspecialised',paste(q[A]==0,' and ',q[S]==0))),
                               expression(atop('Symmetric Specialisation',paste(q[A]==0.9,' and ',q[S]==0.9))),
                               expression(atop('Specialisation on A',paste(q[A]==0.9,' and ',q[S]==0))),
                               expression(atop('Specialisation on S',paste(q[A]==0.0,' and ',q[S]==0.9)))
                      ))
df$ESV = factor(df$ESV,levels=unique(df$ESV))
df$Resource = factor(df$Resource,levels=unique(df$Resource))

p1 <- ggplot(df) + geom_tile(aes(x=Resource,y=ESV,fill=c_ia)) +
  scale_fill_gradient(low="white",high="grey10") + 
  facet_wrap(~Treatment,ncol=2,nrow=2, labeller = label_parsed) + 
  labs(x = 'Resource',y = 'Species', fill = expression(c[i*alpha]))+
  theme_minimal()+  
  theme(axis.text = element_blank(),
        strip.text = element_text(size=15),
        axis.title = element_text(size=20),
        legend.title= element_text(size=20)) 
p1
ggsave(path=path.p, "fig4S4.pdf", width=7.8, height=8, onefile = TRUE)

#.. write data figure as csv file
ds<-as.data.frame(df)
fname<-paste("fig4S4_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)

