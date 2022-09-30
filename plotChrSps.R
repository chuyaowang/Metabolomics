library(dplyr)
library(xcms)
library(stringr)
library(gtools)

register(SerialParam())

datadir <- "./data/20220927" # Edit data directory
fls_cent <- list.files(path = datadir, pattern = "cent",
                  all.files = TRUE, full.names = TRUE,
                  recursive = FALSE, ignore.case = FALSE,
                  include.dirs = FALSE, no.. = TRUE) %>%
  mixedsort(decreasing = T)

fls <- list.files(path = datadir, pattern = "ul.mzML",
                  all.files = TRUE, full.names = TRUE,
                  recursive = FALSE, ignore.case = FALSE,
                  include.dirs = FALSE, no.. = TRUE) %>%
  mixedsort(decreasing = T)

data <- readMSData(fls, mode = "onDisk")
data_cent <- readMSData(fls_cent, mode = "onDisk")

mzs <- c(326.11085,
         269.07681,
         313.10303,
         368.15646,
         340.11393)