library(dplyr)
library(openxlsx)
library(stringr)

fls <- "./data/features.xlsx"

data <- read.xlsx(fls)

d2 <- data %>% select(Name) %>% unlist

compounds <- d2[!str_detect(d2,"F[:digit:]|Study")]
files <- d2[str_detect(d2,"F[:digit:]")] %>% unique

out <- matrix(rep(0,times = length(compounds)*length(files)),ncol=length(files)) %>%
  as.data.frame
rownames(out) <- compounds
colnames(out) <- files

d3 <- d2[!str_detect(d2,"Study")]

compoundIdx <- 0
for (i in 1:length(d3)) {
  if (!str_detect(d3[i],"F[:digit:]")) {
    compoundIdx <- compoundIdx + 1
  } else if (str_detect(d3[i],"F[:digit:]")) {
    fileIdx <- which(files == d3[i])
    out[compoundIdx,fileIdx] <- 1
  }
}

out <- mutate(out, total = apply(out,MARGIN=1,sum))
