# Libraries
library(viridis)
library(showtext)
library(tiff)
library(openxlsx)
library(cowplot)
library(caret)
library(tidyverse)
library(dplyr)
library(gtools)

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

## Volt -----

datadir <- "./data/20230308/ms_data.xlsx"
data_volt <- read.xlsx(datadir,sheet=1) %>% 
  pivot_longer(.,2:6,names_to = "Volts")

min_max <- function(x) {
  process <- preProcess(as.data.frame(x),method=c("range"))
  out <- predict(process,as.data.frame(x)) %>% 
    unlist
  return(out)
}

data_volt <- data_volt %>% 
  group_by(Compound) %>% 
  mutate(value_norm = min_max(value)) %>% 
  ungroup

data_volt$Compound <- factor(data_volt$Compound, levels = c("Pyr",
                                                            "Glc",
                                                            "UDP-GlcNAc",
                                                            "PEP",
                                                            "3PG",
                                                            "S7P",
                                                            "G6P/F6P",
                                                            "R5P"))

p4 <- ggplot(data_volt, aes(x=Volts,y=value,group=Compound,color=Compound)) +
  geom_line(linewidth = 0.2) +
  geom_point(size=0.2) +
  geom_hline(yintercept = 1e6, linetype = "dashed", color = "red", linewidth = 0.1) +
  scale_color_viridis(discrete = TRUE, alpha=1, begin=0.8, end=0.2, option="G") +
  ylab("Intensity") +
  xlab("Voltage (V)") +
  theme_void() +
  theme(text = element_text(size = 6,family="Times"),
        plot.margin = margin(t=4,b=1,l=4),
        legend.position = "right",
        legend.justification = "top",
        legend.background = element_rect(fill = "white",color=NA),
        legend.direction = "vertical",
        legend.key.size = unit(6,"pt"),
        legend.key = element_rect(fill=NA,color=NA),
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 4),
        panel.grid.major = element_line(linewidth=0.1,color="grey"),
        axis.line = element_line(linewidth = 0.1,color="black"),
        axis.text = element_text(size = 4),
        axis.text.y = element_text(margin=margin(r=1)),
        axis.text.x = element_text(margin=margin(t=1)),
        axis.title.x = element_text(size = 4, margin = margin(t=2)),
        axis.title.y = element_text(size = 4, angle = 90, margin = margin(r=2))) +
  scale_y_continuous(breaks = sort(c(seq(0, max(data_volt$value), length.out=5), 1e6)),labels=scales::scientific_format(digits = 2),expand=c(0,0),limits=c(0,NA))

tiff(file = "~/Documents/work/paper/materials/figures/volt.tiff", width = 8.3, height = 9, units = "cm", res = 600)
p4
dev.off()

## Sheath gas -----

data_sheath <- read.xlsx(datadir,sheet=2) %>% 
  pivot_longer(.,2:5,names_to = "Pressure")

data_sheath$Compound <- factor(data_sheath$Compound, levels = c("Pyr",
                                                            "Glc",
                                                            "UDP-GlcNAc",
                                                            "PEP",
                                                            "3PG",
                                                            "S7P",
                                                            "G6P/F6P",
                                                            "R5P"))

p5 <- ggplot(data_sheath, aes(x=as.numeric(Pressure),y=value,group=Compound,color=Compound)) +
  geom_line(linewidth = 0.2) +
  geom_point(size=0.2) +
  geom_hline(yintercept = 1e6, linetype = "dashed", color = "red", linewidth = 0.1) +
  scale_color_viridis(discrete = TRUE, alpha=1, begin=0.8, end=0.2, option="G") +
  ylab("Intensity") +
  xlab("Sheath Gas (Arb)") +
  theme_void() +
  theme(text = element_text(size = 6,family="Times"),
        plot.margin = margin(t=4,b=1,l=4),
        legend.position = "right",
        legend.justification = "top",
        legend.background = element_rect(fill = "white",color=NA),
        legend.direction = "vertical",
        legend.key.size = unit(6,"pt"),
        legend.key = element_rect(fill=NA,color=NA),
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 4),
        panel.grid.major = element_line(linewidth=0.1,color="grey"),
        axis.line = element_line(linewidth = 0.1,color="black"),
        axis.text = element_text(size = 4),
        axis.text.y = element_text(margin=margin(r=1)),
        axis.text.x = element_text(margin=margin(t=1)),
        axis.title.x = element_text(size = 4, margin = margin(t=2)),
        axis.title.y = element_text(size = 4, angle = 90, margin = margin(r=2))) +
  scale_y_continuous(breaks = sort(c(seq(0, max(data_sheath$value), length.out=5), 1e6)),labels=scales::scientific_format(digits = 2),expand=c(0,0),limits=c(0,NA))

## Aux gas -----

data_aux <- read.xlsx(datadir,sheet=3) %>% 
  pivot_longer(.,2:4,names_to = "Pressure")

data_aux$Compound <- factor(data_aux$Compound, levels = c("Pyr",
                                                          "Glc",
                                                          "UDP-GlcNAc",
                                                          "PEP",
                                                          "3PG",
                                                          "S7P",
                                                          "G6P/F6P",
                                                          "R5P"))

p6 <- ggplot(data_aux, aes(x=factor(Pressure,levels=c("5","10","15")),y=value,group=Compound,color=Compound)) +
  geom_line(linewidth = 0.2) +
  geom_point(size=0.2) +
  geom_hline(yintercept = 1e6, linetype = "dashed", color = "red", linewidth = 0.1) +
  scale_color_viridis(discrete = TRUE, alpha=1, begin=0.8, end=0.2, option="G") +
  ylab("Intensity") +
  xlab("Auxiliary Gas (Arb)") +
  theme_void() +
  theme(text = element_text(size = 6,family="Times"),
        plot.margin = margin(t=4,b=1,l=4),
        legend.position = "right",
        legend.justification = "top",
        legend.background = element_rect(fill = "white",color=NA),
        legend.direction = "vertical",
        legend.key.size = unit(6,"pt"),
        legend.key = element_rect(fill=NA,color=NA),
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 4),
        panel.grid.major = element_line(linewidth=0.1,color="grey"),
        axis.line = element_line(linewidth = 0.1,color="black"),
        axis.text = element_text(size = 4),
        axis.text.y = element_text(margin=margin(r=1)),
        axis.text.x = element_text(margin=margin(t=1)),
        axis.title.x = element_text(size = 4, margin = margin(t=2)),
        axis.title.y = element_text(size = 4, angle = 90, margin = margin(r=2))) +
  scale_y_continuous(breaks = sort(c(seq(0, max(data_aux$value), length.out=5), 1e6)),labels=scales::scientific_format(digits = 2),expand=c(0,0),limits=c(0,NA))

## Total -----

p_t2 <- plot_grid(p5,p6,p4,nrow=1,labels="AUTO", label_fontfamily = "Times",label_size = 6, label_fontface = "bold")
tiff(file = "~/Documents/work/paper/materials/figures/ms_opt.tiff", width = 17.1, height = 6, units = "cm", res = 600)
p_t2
dev.off()

# source("stackedEIC.R")

# p_t <- plot_grid(p1,p2,p3,p4,p5,p6,nrow=2,labels = "AUTO", label_fontfamily = "Times",label_size = 6, label_fontface = "bold")
# tiff(file = "~/Documents/work/paper/materials/figures/opt.tiff", width = 17.1, height = 10, units = "cm", res = 600)
# p_t
# dev.off()