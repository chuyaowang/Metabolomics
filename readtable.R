library(dplyr)
library(openxlsx)
library(stringr)
library(ggplot2)

data <- read.xlsx("./Compounds.xlsx",sheet=1,na.strings = "")

areanames <- colnames(data)[str_detect(colnames(data),"Area")][2:14]
colnames(data)[9] <- "deltamasserror"
d2 <- data %>%
  dplyr::filter(mzCloud.Best.Match > 60 &
           deltamasserror<10 &
           deltamasserror>-10) %>%
  select(areanames)

totalnum <- apply(!is.na(d2),MARGIN=2,sum)

df <- as.data.frame(list(seq_along(totalnum),totalnum))
colnames(df) <- c("Files","num_detected")
# df <- df[2:14,]

ggplot(df,aes(x=Files, y=num_detected)) +
  geom_line() +
  geom_point()

