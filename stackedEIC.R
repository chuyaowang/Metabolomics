library(openxlsx)
library(dplyr)
library(ggplot2)
library(viridis)
library(showtext)
library(cowplot)
showtext_auto()
font_add("Times", regular = "/System/Library/Fonts/Supplemental/Times New Roman.ttf",bold = "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf",italic = "/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf",bolditalic = "/System/Library/Fonts/Supplemental/Times New Roman Bold Italic.ttf")

equal_breaks <- function(n = 3, s = 0.05, r = 0,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    
    if ((min(x)+d)<0) {
      seqq <- seq(0, max(x)-d, length.out=n)
    }
    else {
      seqq <- seq(0,max(x)-d,length.out=n)
      # was seq(min(x)+d,max(x)-d,length.out=n)
    }
    
    if(seqq[2]-seqq[1] < 10^(-r)) seqq else round(seqq, r)
  }
}

# ph -----

datadir <- "data/20230308/ph_data.xlsx"
cpds <- c("Lac","Pyr","Glc","PEP","G6P/F6P","3PG")
concs <- c(0.1,0.5,1,2,3,3.5)

data <- lapply(seq_along(concs),function(i) {
  
  sheet <- read.xlsx(datadir,sheet = i,startRow = 5,colNames = FALSE)
  formated <- lapply(seq_along(cpds),function(j) {
    idxs <- c(1,2) + 2*(j-1)
    df <- sheet[,idxs] %>% 
      mutate(Compound = cpds[j],
             Concentration = concs[i])
    colnames(df)[1:2] <- c("Time","Intensity")
    return(df)
    
  }) %>% bind_rows
  return(formated)
}) %>% bind_rows

data$Compound <- factor(data$Compound, levels=c("Pyr","Glc","PEP","3PG","G6P/F6P","Lac"))

p1 <- ggplot(data,aes(x=Time,y=Intensity,color=factor(Concentration))) + 
  facet_wrap(~Compound,ncol=1,scales = "free",strip.position = "left") +
  geom_line(linewidth=0.2) +
  scale_color_viridis(name = "Concentration \n(%)",discrete = TRUE, alpha=1, begin=0.2, end=0.8, option="G") +
  xlab("Retention Time (min)") +
  ylab("Intensity") +
  theme_void() +
  theme(strip.placement = "outside",
        text = element_text(size = 6,family="Times"),
        strip.text = element_text(size = 4),
        plot.margin = margin(t=4,b=1,l=5),
        legend.position = "right",
        legend.justification = "top",
        legend.background = element_rect(fill = "white",color=NA),
        legend.direction = "vertical",
        legend.key.size = unit(10,"pt"),
        legend.key = element_rect(fill=NA,color=NA),
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 4),
        panel.grid.major = element_blank(),
        axis.line = element_line(linewidth = 0.1,color="black"),
        axis.text = element_text(size = 4),
        axis.text.y = element_text(margin=margin(r=1)),
        axis.text.x = element_text(margin=margin(t=1)),
        axis.title.x = element_text(size = 4),
        panel.spacing = unit(0.5,"lines")) +
  scale_y_continuous(breaks = equal_breaks(n=3,s=0.05), labels=scales::scientific_format(digits = 2), expand=c(0,0),limits=c(0,NA)) +
  scale_x_continuous(expand=c(0,0),limits=c(0,NA))

tiff(file = "~/Documents/work/paper/materials/figures/ph.tiff", width = 8.3, height = 7.5, units = "cm", res = 600)
p1
dev.off()

## salt -----

datadir <- "data/20230308/salt_data.xlsx"
cpds <- c("Lac","Pyr","Glc","PEP","G6P/F6P","3PG")
concs <- c(10,20,25)

data_salt <- lapply(seq_along(concs),function(i) {
  
  sheet <- read.xlsx(datadir,sheet = i,startRow = 5,colNames = FALSE)
  formated <- lapply(seq_along(cpds),function(j) {
    idxs <- c(1,2) + 2*(j-1)
    df <- sheet[,idxs] %>% 
      mutate(Compound = cpds[j],
             Concentration = concs[i])
    colnames(df)[1:2] <- c("Time","Intensity")
    return(df)
    
  }) %>% bind_rows
  return(formated)
}) %>% bind_rows

data_salt$Compound <- factor(data_salt$Compound, levels=c("Pyr","Glc","PEP","3PG","G6P/F6P","Lac"))

p2 <- ggplot(data_salt,aes(x=Time,y=Intensity,color=factor(Concentration, levels = c("10","25","20")))) + 
  facet_wrap(~Compound,ncol=1,scales = "free",strip.position = "left") +
  geom_line(linewidth=0.2) +
  scale_color_viridis(name = "Concentration \n(mM)", breaks = c("10","20","25"), discrete = TRUE, alpha=1, begin=0.2, end=0.8, option="G") +
  xlab("Retention Time (min)") +
  ylab("Intensity") +
  theme_void() +
  theme(strip.placement = "outside",
        text = element_text(size = 6,family="Times"),
        strip.text = element_text(size = 4),
        plot.margin = margin(t=4,b=1,l=5),
        legend.position = "right",
        legend.justification = "top",
        legend.background = element_rect(fill = "white",color=NA),
        legend.direction = "vertical",
        legend.key.size = unit(10,"pt"),
        legend.key = element_rect(fill=NA,color=NA),
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 4),
        panel.grid.major = element_blank(),
        axis.line = element_line(linewidth = 0.1,color="black"),
        axis.text = element_text(size = 4),
        axis.text.y = element_text(margin=margin(r=1)),
        axis.text.x = element_text(margin=margin(t=1)),
        axis.title.x = element_text(size = 4),
        panel.spacing = unit(0.5,"lines")) +
  scale_y_continuous(breaks = equal_breaks(n=3,s=0.05), labels=scales::scientific_format(digits = 2), expand = c(0,0),limits=c(0,NA)) +
  scale_x_continuous(expand=c(0,0),limits=c(0,NA))

tiff(file = "~/Documents/work/paper/materials/figures/salt.tiff", width = 8.3, height = 7.5, units = "cm", res = 600)
p2
dev.off()

## Column -----

datadir <- "data/20230308/column_data.xlsx"
cpds <- c("Lac","Pyr","Glc","PEP","G6P/F6P","3PG")
columns <- c("BEH Amide","C18")

data_column <- lapply(seq_along(columns),function(i) {
  
  sheet <- read.xlsx(datadir,sheet = i,startRow = 5,colNames = FALSE)
  formated <- lapply(seq_along(cpds),function(j) {
    idxs <- c(1,2) + 2*(j-1)
    df <- sheet[,idxs] %>% 
      mutate(Compound = cpds[j],
             Column = columns[i])
    colnames(df)[1:2] <- c("Time","Intensity")
    return(df)
    
  }) %>% bind_rows
  return(formated)
}) %>% bind_rows

data_column$Compound <- factor(data_column$Compound, levels=c("Pyr","Glc","PEP","3PG","G6P/F6P","Lac"))

p3 <- ggplot(data_column,aes(x=Time,y=Intensity,color=factor(Column, levels = c("C18", "BEH Amide")))) + 
  facet_wrap(~Compound,ncol=1,scales = "free",strip.position = "left") +
  geom_line(linewidth=0.2) +
  scale_color_viridis(name = "Column", labels = c("C18","BEH\nAmide"), discrete = TRUE, alpha=1, begin=0.2, end=0.8, option="G") +
  xlab("Retention Time (min)") +
  ylab("Intensity") +
  theme_void() +
  theme(strip.placement = "outside",
        text = element_text(size = 6,family="Times"),
        strip.text = element_text(size = 4),
        plot.margin = margin(t=4,b=1,l=5),
        legend.position = "right",
        legend.justification = "top",
        legend.background = element_rect(fill = "white",color=NA),
        legend.direction = "vertical",
        legend.key.size = unit(10,"pt"),
        legend.key = element_rect(fill=NA,color=NA),
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 4),
        panel.grid.major = element_blank(),
        axis.line = element_line(linewidth = 0.1,color="black"),
        axis.text = element_text(size = 4),
        axis.text.y = element_text(margin=margin(r=1)),
        axis.text.x = element_text(margin=margin(t=1)),
        axis.title.x = element_text(size = 4),
        panel.spacing = unit(0.5,"lines")) +
  scale_y_continuous(breaks = equal_breaks(n=3,s=0.05), labels=scales::scientific_format(digits = 2),expand=c(0,0),limits=c(0,NA)) +
  scale_x_continuous(expand=c(0,0),limits=c(0,NA))

tiff(file = "~/Documents/work/paper/materials/figures/column.tiff", width = 8.3, height = 7.5, units = "cm", res = 600)
p3
dev.off()

## Total -----

p_t1 <- plot_grid(p3,p2,p1,labels = c("A","B","C"),label_fontfamily = "Times",label_size = 6, label_fontface = "bold",nrow = 1)

tiff(file = "~/Documents/work/paper/materials/figures/lc_opt.tiff", width = 17.1, height = 6, units = "cm", res = 600)
p_t1
dev.off()