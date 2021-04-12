" R script for Figure 2- figure supplement 2 in Estrela et al (2021) 
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

#--  To plot relative abundance
#.. sum relative abundance of each family for each community
dp1<-select(dat,-Sequence,-Order,-Transfer,-ESV,-Genus)

sum_relab_fam<- dp1 %>%
  dplyr::group_by(Inoculum,Carbon_Source,Replicate,s_perc,cs_type,Experiment,nCS,Family,CS1,CS2) %>%
  dplyr::summarise(
    sum_relab=sum(Relative_Abundance))
dp1<-as.data.frame(sum_relab_fam)

#.. SELECT minimum abundance for plot- taxa below threshold are 'other'
min.ab<-0.01
dfp.rare<-subset(dp1, sum_relab<min.ab)
dfp.rare$Family_type<-'Other'
dfp.com<-subset(dp1, sum_relab>=min.ab)
dfp.com$Family_type<-dfp.com$Family
dfp<-rbind(dfp.com,dfp.rare)

#.. Colour scheme for plots
#. Family colors 

colAlc<-'darkolivegreen1'
colBruc<- 'darkolivegreen2'
colChit<- 'darkolivegreen3'
colCom<- 'darkolivegreen4'
colEntb<- blues9[7]
colEntc<- 'skyblue1'
colFlav<-'skyblue2'
colLach<-'skyblue3'
colMorax<- "darkorchid4"
colOxal<- "darkorchid4"
colPseu<-"darkorchid2" 
colRhiz<-'tan1'
colSphing<-'tan2'
colXanth<-'tan3'
colOther<-'gray85'
colEntbf<-'#516091'
colMoraxf<-'#74BEC1'
colPseuf<-'#E2F0CB'
colRhizf<-'#996666'

#.. Make plot
#. select mixed carbon sources only
dp1<-subset(dfp, nCS ==2)

#. define order of duo CS names to plot
levels(factor(dp1$Carbon_Source))
dp1$cs.ord <- factor(dp1$Carbon_Source, levels = c("glucose_fructose","glucose_cellobiose", "glucose_ribose", "glucose_glycerol",
                                                   "glucose_benzoate","glucose_fumarate","glucose_glycine","glucose_glutamine",
                                                   "succinate_fructose","succinate_cellobiose", "succinate_ribose", "succinate_glycerol",
                                                   "succinate_benzoate","succinate_fumarate","succinate_glycine","succinate_glutamine"))

p1<-ggplot(dp1, aes(x=factor(Replicate), y=sum_relab, fill=Family_type))+
  geom_bar(color="black", stat="identity", position="stack") + 
  labs(y= "Relative abundance", x="Replicate community")+
  facet_grid(Inoculum~cs.ord,scales = "free")+
  scale_fill_manual(values=c(colBruc,colChit,colCom,colEntbf,colEntc,colFlav,colLach,colMoraxf,colOxal,colPseuf,colRhizf,colSphing,colXanth, colOther))+
  scale_y_continuous( breaks=seq(0,1,1)) +
  guides(fill=guide_legend(title="Family"))+
  theme_minimal()+
  theme(axis.text.x = element_text(size=14),
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
ggsave(path=path.p, paste0("fig2S2.pdf"), width=20, height=4, onefile = TRUE)

#.. write data figure as csv file
ds<-select(dp1, Inoculum, Carbon_Source, Family, Replicate, sum_relab, Family_type)
fname<-paste("fig2S2_data.csv")
write.csv(ds, file.path(path.ds, fname),row.names=FALSE)



