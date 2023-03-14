library(ggplot2)
library(ggridges)
library(openxlsx)
library(dplyr)
library(viridis)
library(hrbrthemes)
library(tidyr)

data <- read.xlsx("data/20230308/ninefigures.xlsx")
d1 <- data %>% 
  pivot_longer(.,colnames(data)[2:10],names_to="group",values_to = "height") %>%
  # filter(height>0) %>%
  filter(time<20) %>%
  group_by(group) %>% 
  mutate(heightnorm=height/max(height)) %>%
  ungroup
d1$group <- factor(d1$group, levels=rev(c("Pyr",
                                          "Glc",
                                          "UDP-GlcNAc",
                                          "PEP",
                                          "3-PG",
                                          "S7P",
                                          "GlcN6P",
                                          "G6P/F6P",
                                          "R5P")))
add <- seq(from=20,to=23,by=d1$time %>% diff %>% mean)
add.df <- data.frame(
  time = add,
  group = sample(d1$group, length(add), replace=TRUE),
  height = 0,
  heightnorm = 0
)
d1 <- rbind(d1,add.df)

ggplot(d1,aes(x=time,y=group,group=group,height=heightnorm,fill=group))+
  geom_ridgeline(scale=1.3, size=1) +
  scale_fill_viridis(discrete = TRUE, alpha=0.7, begin=0.2, end=0.8, option="turbo") +
  theme_ridges(font_size=40,line_size=1,center_axis_labels = TRUE) +
  xlab("Retention Time (min)") +
  theme(legend.position = "none",axis.title.y = element_blank())
