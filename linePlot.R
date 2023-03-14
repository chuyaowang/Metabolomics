library(openxlsx)
library(ggplot2)
library(tidyr)
library(viridis)
library(ggrepel)
library(gtools)
library(plyr)
library(dplyr)

data <- read.xlsx("data/20230308/calibration.xlsx")
data <- data %>% 
  pivot_longer(.,colnames(data)[2:14],names_to = "conc",values_to = "area")
data$conc <- as.numeric(data$conc)
data$area <- as.numeric(data$area)
data$metabolite <- factor(data$metabolite)
data <- data %>% filter(!is.na(area))


lm_eqn <- function(data){
  m <- lm(area ~ conc,data);
  eq <- substitute(atop(italic(y) == a + b %.% italic(x),~~italic(r)^2~"="~r2), 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 4)))
  as.character(as.expression(eq));
}
eq <- ddply(data,.(metabolite),lm_eqn)

ggplot(data, aes(x=conc, y=area, color=metabolite)) +
  geom_smooth(method=lm, se=TRUE, na.rm = TRUE, linewidth=2, formula = y ~ x) +
  geom_point(size=2,na.rm=TRUE) +
  geom_text(data=eq,aes(x = 25, y=25,label=V1),size=8, vjust=-5, hjust= -0.1, parse = TRUE, inherit.aes=FALSE) +
  facet_wrap(metabolite~.,scales = "free_x") +
  scale_color_viridis(option="turbo",begin=0.2,end=0.8,discrete=T,alpha=0.8) +
  theme_minimal(base_size = 32) +
  xlab("Concentration (ppm)") +
  ylab("Peak Area") +
  theme(legend.position = "none",
        strip.text.x = element_text(margin = margin(b=6,t=6),size=32),
        strip.background = element_rect(colour = "black", fill=NA),
        axis.title.y = element_text(vjust=.5))
