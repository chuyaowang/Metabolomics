library(openxlsx)
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(ggrepel)

featPerFile <- read.xlsx("data/20230307/cellsneg_features.xlsx")
featPerFile <- featPerFile %>% group_by(Study.File.ID) %>% summarise(feats=length(Molecular.Weight),avgArea=mean(Area))

cpdPerFile <- read.xlsx("data/20230307/cellsneg_compounds.xlsx")
cpdPerFile <- cpdPerFile %>% group_by(Study.File.ID) %>% summarise(cpds=length(Calc..MW),meanInt=mean(`Intensity.(Max.)`))

files <- read.xlsx("data/20230307/cellsneg_files.xlsx")
files <- left_join(files,featPerFile,by=c("ID"="Study.File.ID"))
files <- left_join(files,cpdPerFile,by=c("ID"="Study.File.ID"))

files.extraction.solvent <- files %>% filter(Cell.Count == 1.5e6)
files.cellcount <- files %>% filter(Cell.Count != 1.5e6)

ggplot(files.extraction.solvent,aes(fill=Solvent,y=cpds,x=Extraction.method,color=Name)) +
  geom_bar(position="dodge",just=.5,stat="identity",width=0.4) +
  theme_minimal() +
  scale_fill_viridis(discrete=TRUE,alpha = 1,begin=0.2,end=0.8,option="D")

ggplot(files.extraction.solvent,aes(x=cpds,y=meanInt,color=Solvent)) +
  geom_point() +
  geom_text_repel(aes(label = Extraction.method),hjust=0,vjust=1.4,size=4,show.legend = FALSE) +
  theme_minimal() +
  scale_color_viridis(option="A",discrete=TRUE,begin=0.2,end=0.8)

ggplot(files,aes(x=cpds,y=meanInt,color=factor(Cell.Count))) +
  geom_point()+
  geom_text_repel(aes(label = Extraction.method),hjust=0,vjust=1.4,size=4,show.legend = FALSE) +
  theme_minimal() +
  scale_color_viridis(option="turbo",discrete=TRUE,begin=0.2,end=0.8)
  
