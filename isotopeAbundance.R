library(dplyr)
library(xcms)
library(stringr)
library(gtools)
library(MSnbase)
library(purrr)

register(SerialParam())

datadir <- "./data/20220928" # Edit data directory
fls <- list.files(path = datadir, pattern = "^.*\\.mzML$",
                   all.files = TRUE, full.names = TRUE,
                   recursive = FALSE, ignore.case = FALSE,
                   include.dirs = FALSE, no.. = TRUE) %>%
  mixedsort(decreasing = T)

## Centroid the data and save to new file, only run once -----
# data_cent <- function(file) {
#   data <- readMSData(file,
#                      mode = "onDisk") %>%
#     pickPeaks()
# 
#   fls_new <- fileNames(data) %>% 
#     str_sub(end=-6) %>% 
#     paste("cent",sep = "_") %>%
#     paste("mzML",sep = ".")
#   writeMSData(data, file = fls_new)
# }
# 
# data_cent(file = fls)

## Read data -----
data <- readMSData(fls, mode = "onDisk")

## Edit parameters -----
all_labeled <- 79.04305 # mz ratio
mzs <- c(76.0393,77.04266,77.03634,78.04601,78.03969,79.04305) # mz ratio
label <- c("noLabel","C13","N15","2C13","C13_N15","2C13_N15")
ppm <- 5
rt_range <- c(1.2,1.4) # in minutes

## Getting intensities -----
chr <- data %>%
  filterMz(mz = all_labeled*c(1-ppm*1e-6,1+ppm*1e-6)) %>%
  filterRt(rt = 60*rt_range) %>%
  chromatogram(aggregationFun = "max")

int_mz <- lapply(mzs, function(x) {
  res <- list()
  for (f in 1:length(fileNames(data))) {
    sp <- data %>%
      filterFile(f) %>%
      filterAcquisitionNum(n = intensity(chr[1,f]) %>% 
                             which.max %>%
                             names %>%
                             str_sub(start = -4) %>% # may need to be changed
                             as.integer()  
      ) %>%
      filterMz(mz = x*c(1-ppm*1e-6,1+ppm*1e-6))
    
    vals <- c(mz(sp),intensity(sp)) %>% unlist
    
    if (sum(vals)==0) {
      res[[f]] <- c(f,x,0)
    } else {
    res[[f]] <- c(f,vals)
    }
  }
  
  return(res)
})

## Output -----
out <- lapply(seq_along(int_mz), function (x) {
  d <- int_mz[[x]] %>%
    as.data.frame() %>%
    t %>%
    as.data.frame()
  
  colnames(d) <- c("File",
                   paste("mz",x,sep = ""),
                   paste("Intensity",x,sep = ""))
  
  return(d)
}) %>%
  purrr::reduce(left_join, by = "File") %>%
  mutate(total = Intensity1+Intensity2+Intensity3+Intensity4+Intensity5+Intensity6) %>%
  mutate(Result = (0*Intensity1 + 
           1*(Intensity2+Intensity3) + 
           2*(Intensity4+Intensity5) + 
           3*Intensity6)/(3*total))

rsd <- sd(out$Result)/mean(out$Result)


