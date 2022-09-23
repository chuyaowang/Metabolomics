library(xcms)
library(CAMERA)
library(dplyr)

# Import data -----
datadir <- "./data/20220824"
file <- list.files(path = datadir, pattern = "^.*\\.mzML$",
                   all.files = TRUE, full.names = TRUE,
                   recursive = FALSE, ignore.case = FALSE,
                   include.dirs = FALSE, no.. = TRUE)
pd <- data.frame(file = basename(file),
                 sample = c("Sample1", "Sample2","Sample3"),
                 group = c("Blank","Sample","Sample"))
data <- readMSData(file, # filenames to load
                   # Stores data; stores sample info into @phenoData@data
                   # More info about each variable: file, injection. sample etc. can be entered in varMetadata
                   pdata = new("NAnnotatedDataFrame", pd), 
                   mode = "onDisk" # Only load metadata into memory and leave the m/z ratio vs intensity peaks on disk
) 

pos <- data %>%
  filterPolarity(1) %>%
  filterMsLevel(1)

neg <- data %>%
  filterPolarity(0) %>%
  filterMsLevel(1)

dataA <- list(pos,neg)
rm(data)
# Data summary -----
featureData <- fData(data)

polarity(data) %>% table
msLevel(data) %>% table

# Data processing -----
chr <- 

























xs <- xcmsSet(
  files=file, 
  method="centWaveWithPredictedIsotopeROIs",
  ppm=10,
  mslevel = 1,
  snthresh = 5,
  integrate = 1,
  peakwidth = c(2,30),
  prefilter = c(1,1e4),
  mzCenterFun = "wMeanApex3",
  mzdiff = -0.001
)

an <- xsAnnotate(xs, polarity="negative") # constructor; extracts peak table
an <- groupFWHM(an, perfwhm = 1) # group peaks by retention time
an <- findIsotopesWithValidation(object = an, ppm = 10,
                                 mzabs = 0.01, intval="intb",
                                 maxcharge = 3) # annotate isotopic peaks
## extract annotated peak table
peakTable <- getPeaklist(an) # extract peak list
