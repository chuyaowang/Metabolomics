library(ggplot2)
library(ggpubr)
library(openxlsx)
library(extrafont)
library(dplyr)

data <- read.xlsx("./data/veggies.xlsx")

p <- ggplot(data, aes(x=type, y=value, fill=group)) + 
  geom_boxplot() +
  theme_bw() +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "t.test",
                     paired = FALSE,
                     hide.ns = FALSE,
                     vjust = .5
                     ) +
  facet_wrap(~type, scale="free") +
  theme(axis.text.x=element_blank()) +
  theme(axis.title.y=element_text(angle=0,vjust = 0.5)) +
  theme(text = element_text(size = 12)) +
  xlab("蔬菜种类") +
  ylab("比值") +
  labs(fill = "种类")

t_tests <- compare_means(value~group,data=data,method="t.test",group.by="type")

t_test2 <- data %>%
  group_by(type) %>%
  summarize(conf.int1 = t.test(value~group)$conf.int[1],
            conf.int2 = t.test(value~group)$conf.int[2])
t_tests <- left_join(t_tests,t_test2,by = "type")
