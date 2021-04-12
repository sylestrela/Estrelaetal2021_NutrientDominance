" R script for Figure 4- figure supplement 1 in Estrela et al (2021)
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

#.. Plot r.avg for each strain:CS
#. Family colors 
colEntbf<-'#516091'
colMoraxf<-'#74BEC1'
colPseuf<-'#E2F0CB'
colRhizf<-'#996666'

dp<-dat
xl<-""
yl<-expression('average growth rate (h'^-1* ')')

#. order CS and family for plotting
dp$CS <- factor(dp$CS, levels = c("glucose", "fructose", 'cellobiose', 'ribose','glycerol', 'succinate', "fumarate", "benzoate", 'glutamine', 'glycine'))
dp$Family<-factor(dp$Family, levels=c('Enterobacteriaceae','Pseudomonadaceae','Moraxellaceae', 'Rhizobiaceae'))

ymax<-max(dp$ravg)
ymin<-min(dp$ravg)

p<-ggplot(dp, aes(x=Family, y=ravg, fill=Family))+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  geom_jitter(width = 0.2, size=2, alpha = 0.8, aes(color=Family))+
  theme_classic()+
  facet_wrap(~CS, ncol=5, scales='free_y')+
  scale_y_continuous(limits = c(0, ymax+0.05))+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=18),
    axis.title.y = element_text(size=18),
    axis.title.x =element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linetype='dashed', size=0.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())+
  theme(strip.text.x = element_text(size = 18),
        strip.background = element_blank(),
        axis.ticks.x=element_blank())+
  scale_colour_manual(values=c(colEntbf, colPseuf, colMoraxf, colRhizf))+
  scale_fill_manual(values=c(colEntbf, colPseuf, colMoraxf, colRhizf))+
  labs(x=xl, y=yl)+
  theme(legend.title = element_blank(),
        legend.text=element_text(size=14),
        legend.spacing.x = unit(0.5, 'cm'),
        legend.position="bottom")
p

#. add significance
ps<-p+ 
  stat_compare_means(label = "p.signif", method = "t.test",label.y = 0.32,
                     ref.group = "Enterobacteriaceae") +   # Pairwise comparison against reference
  labs(x=xl, y=yl)
#       caption=paste0("Each dot is the mean for one isolate (n=3-4). r.avg = log2(OD_tf/ OD_ti)/ (tf-ti)) where tf= ",tf, 'h and ti= ',ti,'h. 
#               Significance (*) is measured by comparing ravg between Enterobacteriaceae (reference) and each other family. Paired t-test. '))
ps
ggsave(path=path.p, filename=paste0("fig4S1.pdf"), height=5.5, width=12,onefile = TRUE)  

#. write data figure as csv file
ds<-select(dp, -std)
fname<-paste("fig4S1_data.csv")
write.csv(ds, file.path(path.ds, fname), row.names=FALSE)


