library(readr)
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(ggrepel)
library(showtext)

showtext_auto()
font_add("Times", regular = "/System/Library/Fonts/Supplemental/Times New Roman.ttf")

data <- read_tsv("~/R Projects/isocordb/0227_neg2/isocordata0227_neg2.tsv") %>% filter(metabolite %in% c("3PG","F6P","PEP","Pyr","GlcN6P","R5P","Glc"))
data$sample[which(data$sample == "Sample_ 1")] <- "F"
data$sample[which(data$sample == "Sample_ 2")] <- "B"
data <- data %>% mutate(isotopologue1 = paste("M",isotopologue,sep="+"))

p <- ggplot(data, aes(x=sample,y=area,fill=factor(isotopologue1))) +
  scale_y_continuous(labels = scales::percent_format(suffix=""),limits=c(0,NA),expand = c(0,0)) +
  geom_bar(position="fill",stat="identity",width=0.9) +
  facet_wrap(~metabolite,nrow=1) +
  ylab("Percentage (%)") +
  theme(panel.spacing = unit(0,"lines"),
        strip.text.x = element_text(margin = margin(b=2,t=2,l=0,r=0),size=4),
        strip.background = element_rect(colour = NA, fill=NA),
        panel.grid.major = element_line(colour = "grey", linewidth=0.1),
        panel.background = element_rect(fill="white"),
        axis.title.y = element_text(vjust=.5,size=4),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 4),
        axis.ticks = element_line(linewidth = 0.2),
        axis.line = element_line(linewidth = 0.1,color="black"),
        text = element_text(family="Times",size=4),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.position = "right",
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 4),
        legend.key.size = unit(8,"pt"),
        legend.key = element_rect(fill=NA),
        # aspect.ratio = 6.5,
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")) +
  scale_fill_viridis(name = "Mass \nIsotopologue", discrete=TRUE, option="mako",begin=0.9,end=0.3,alpha=1)

tiff(file = "~/Documents/work/paper/materials/figures/abundance.tiff", width = 8.3, height = 7, units = "cm", res = 600)
p
dev.off()
