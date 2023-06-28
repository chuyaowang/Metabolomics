library(openxlsx)
library(dplyr)
library(ggplot2)
library(showtext)
font_add("Songti", regular = "/System/Library/Fonts/Supplemental/Songti.ttc")

data <- "./data/20230516/Ile_abundance.xlsx"
data <- read.xlsx(data,sheet=5,startRow=5)

time <- c(4,8,12,24,36,48,60,72)
ab_jn <- data[8:1,2]
ab_sq <- data[16:9,2]
df_jn <- data.frame(time,ab_jn)
df_sq <- data.frame(time,ab_sq)

showtext_auto()
ggplot(df_sq[2:8,],aes(x=time,y=ab_sq)) +
  geom_point(size=1,na.rm=TRUE) + 
  stat_smooth(method=lm, se=TRUE, na.rm = TRUE, linewidth=1, formula = y ~ log(x)) +
  xlab("时间 (h)") +
  ylab("丰度") +
  ggtitle("上清液异亮氨酸丰度 vs. 时间") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + # Overrides theme_classic title position
  scale_x_continuous(breaks=seq(0,72,8)) +
  scale_y_continuous(breaks=seq(0,0.2,0.05))

ggplot(df_jn,aes(x=time,y=ab_jn)) +
  geom_point(size=1,na.rm=TRUE) +
  stat_smooth(method=lm, se=TRUE, na.rm = TRUE, linewidth=1, formula = y ~ log(x)) +
  xlab("时间 (h)") +
  ylab("丰度") +
  ggtitle("菌泥异亮氨酸丰度 vs. 时间") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + # Overrides theme_classic title position
  scale_x_continuous(breaks=seq(0,72,8)) +
  scale_y_continuous(breaks=seq(0,0.2,0.05))
