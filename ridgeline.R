library(ggplot2)
library(ggridges)
library(openxlsx)
library(dplyr)
library(viridis)
library(hrbrthemes)
library(tidyr)
library(showtext)
library(tiff)

showtext_auto()

data <- read.xlsx("data/20230308/ninefigures.xlsx")
d1 <- data %>% 
  pivot_longer(.,colnames(data)[2:10],names_to="group",values_to = "height") %>%
  # filter(height>0) %>%
  filter(time<20) %>%
  group_by(group) %>% 
  mutate(heightnorm=height/max(height)) %>%
  ungroup
d1$group <- factor(d1$group, levels=c("Pyr",
                                          "Glc",
                                          "UDP-GlcNAc",
                                          "PEP",
                                          "3PG",
                                          "S7P",
                                          "GlcN6P",
                                          "G6P/F6P",
                                          "R5P"))
add <- seq(from=20,to=23,by=d1$time %>% diff %>% mean)
add.df <- data.frame(
  time = add,
  group = sample(d1$group, length(add), replace=TRUE),
  height = 0,
  heightnorm = 0
)
d1 <- rbind(d1,add.df)

maxs <- d1 %>% group_by(group) %>% summarize(x=which(height == max(height)), max_height=max(height),y=max(heightnorm)) %>% ungroup

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

p <- ggplot(d1,aes(x=time,y=group,group=group,height=heightnorm,color=group,fill=group))+
  geom_ridgeline(scale=1, size=.2) +
  scale_color_viridis(discrete = TRUE, alpha=1, begin=0.2, end=0.8, option="G") +
  scale_fill_viridis(discrete = TRUE, alpha=0) +
  theme_ridges(font_size=6,line_size=1,center_axis_labels = TRUE) +
  xlab("Retention Time (min)") +
  ylab("Intensity") +
  theme(text = element_text(size = 6,family="Times"),
        legend.position = "none",
        panel.grid.major = element_line(linewidth=0.3),
        axis.text = element_text(size = 4),
        axis.title.x = element_text(size = 6))

p1 <- ggplot(d1,aes(x=time,y=height,color=group)) +
  facet_wrap(~group,ncol=1,scales = "free",strip.position = "left") +
  geom_line(linewidth=0.3) +
  scale_color_viridis(discrete = TRUE, alpha=1, begin=0.8, end=0.2, option="G") +
  xlab("Retention Time (min)") +
  ylab("Intensity") +
  theme_void() +
  theme(strip.placement = "outside",
        text = element_text(size = 6,family="Times"),
        legend.position = "none",
        plot.margin = margin(t=2,b=1,l=2),
        axis.line = element_line(linewidth = 0.1,color="black"),
        axis.text = element_text(size = 4),
        axis.text.y = element_text(margin=margin(r=1)),
        axis.title.x = element_text(size = 4),
        strip.text = element_text(size = 4),
        panel.spacing = unit(0.5,"lines")) +
  scale_y_continuous(breaks = equal_breaks(n = 3, s = 0.05),labels=scales::scientific_format(digits = 2),expand=c(0,0),limits=c(0,NA)) +
  scale_x_continuous(expand=c(0,0))

tiff(file = "~/Documents/work/paper/materials/figures/EIC.tiff", width = 8.3, height = 9.5, units = "cm", res = 600)
p1
dev.off()
