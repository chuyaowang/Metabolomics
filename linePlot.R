library(openxlsx)
library(ggplot2)
library(tidyr)
library(viridis)
library(ggrepel)
library(gtools)
library(plyr)
library(dplyr)
library(showtext)

showtext_auto()
font_add("Times", regular = "/System/Library/Fonts/Supplemental/Times New Roman.ttf")
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
data <- read.xlsx("data/20230308/calibration.xlsx")
data$metabolite <- factor(data$metabolite,levels = c("Pyr",
                                                         "Glc",
                                                         "UDP-GlcNAc",
                                                         "PEP",
                                                         "3PG",
                                                         "S7P",
                                                         "GlcN6P",
                                                         "G6P",
                                                         "F6P",
                                                         "R5P"))
data <- data %>% 
  pivot_longer(.,colnames(data)[2:14],names_to = "conc",values_to = "area")
data$conc <- as.numeric(data$conc)
data$area <- as.numeric(data$area)
data$metabolite <- factor(data$metabolite)
data <- data %>% filter(!is.na(area))

lm_eqn <- function(data){
  m <- lm(area ~ conc,data);
  eq <- substitute(atop(italic(y) == a + b ~"\u2022" ~ italic(x),~~italic(r)^2~"="~r2), 
                   list(a = format(unname(coef(m)[1]), scientific = F, nsmall = 0, digits = 8),
                        b = format(unname(coef(m)[2]), scientific = F, nsmall =0, digits = 8),
                        r2 = format(summary(m)$r.squared, digits = 4, nsmall = 4)))
  as.character(as.expression(eq));
}
eq <- ddply(data,.(metabolite),lm_eqn)

equal_breaks <- function(n = 3, s = 0.05, r = 0,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    if(seq[2]-seq[1] < 10^(-r)) seq else round(seq, r)
  }
}

p <- ggplot(data, aes(x=conc, y=area, color=metabolite)) +
  geom_smooth(method=lm, se=TRUE, na.rm = TRUE, linewidth=.3, formula = y ~ x) +
  geom_point(size=.2,na.rm=TRUE) +
  geom_text(data=eq,aes(x=0, y=0, label=V1),size=1.5, vjust=-1, hjust= -0.05, parse = TRUE, inherit.aes=FALSE, family="Times") +
  facet_wrap(metabolite~.,scales = "free",ncol=3) +
  scale_color_viridis(option="mako",begin=0.8,end=0.2,discrete=T,alpha=0.8) +
  xlab("Concentration (ppm)") +
  ylab("Peak Area") +
  theme_void() +
  theme(legend.position = "none",
        plot.margin = margin(t=2,b=2,l=2,r=5),
        strip.text.x = element_text(margin = margin(b=2,t=2),size=6),
        strip.background = element_rect(color="black"),
        panel.grid.major = element_line(colour = "grey", linewidth=0.1),
        # panel.background = element_rect(fill="white"),
        axis.line = element_line(linewidth = 0.1,color="black"),
        axis.title.x = element_text(size = 6, margin = margin(t=2)),
        axis.title.y = element_text(size = 6, angle = 90, margin = margin(r=2)),
        axis.text = element_text(size = 4),
        axis.text.y = element_text(margin=margin(r=1)),
        axis.text.x = element_text(margin=margin(t=1)),
        axis.ticks = element_blank(),
        text = element_text(family="Times",size=6)) +
  scale_y_continuous(breaks = function(x){seq(0,max(x),length.out=5)},labels=scales::scientific_format(digits=2),limits=c(0,NA),expand=c(0,0))+
  scale_x_continuous(limits = c(0,NA),expand=c(0,0))

tiff(file = "~/Documents/work/paper/materials/figures/cali.tiff", width = 17.1, height = 10, units = "cm", res = 600)
p
dev.off()