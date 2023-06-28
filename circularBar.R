library(readr)
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(ggrepel)

data <- read_tsv("~/R Projects/isocordb/0227_neg2/isocordata0227_neg2.tsv")

data <- data %>% 
  filter(sample=="Sample_ 2") %>% # sample 2 is bacteria body
  filter(metabolite %in% c("3PG","F6P","PEP","Pyr","GlcN6P","R5P","Glc")) %>%
  group_by(metabolite) %>% 
  mutate(
    fraction = 100*area / sum(area),
    ymax = cumsum(fraction),
    ymin = lag(ymax, default = 0),
    labelPosition = (ymax + ymin) / 2,
    label = paste0("M+",isotopologue, "\n", round(fraction,digits=2),"%")
  ) %>% 
  ungroup

# Calculate abundance
abundance <- data %>% 
  group_by(metabolite) %>% 
  summarize(abundance=sum(isotopologue*0.01*fraction)/max(isotopologue))

data <- data %>%
  filter(fraction != 0)
data$metabolite[data$metabolite=="F6P"] <- "F6P/G6P"

p <- ggplot(data, aes(ymax = ymax, ymin = ymin, xmax = 10, xmin = 7, fill = factor(isotopologue))) +
  facet_wrap(~`metabolite`) +
  geom_rect() +
  scale_fill_viridis(option="turbo",begin=0.2,end=0.8,discrete = TRUE,alpha=0.9) +
  coord_polar(theta = "y") +
  xlim(c(2, 10)) +
  geom_label_repel(x=10,aes(y = labelPosition, label = label, fill = factor(isotopologue)), size = 10, fontface = "bold") +
  theme_void() +
  theme(legend.position = "none",
        strip.text.x = element_text(margin = margin(b=6,t=6),size=32),
        strip.background = element_rect(colour = "black", fill=NA))


