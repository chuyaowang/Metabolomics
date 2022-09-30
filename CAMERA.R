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

# Data summary -----
featureData <- fData(data)

polarity(data %>% filterPolarity(1)) %>% table
msLevel(data) %>% table

# Data processing -----
chr <- chromatogram(data, aggregationFun = "max")

dim(chr)

plot(chr %>% `[[`(1,1))
