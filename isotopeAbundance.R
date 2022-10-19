## libs and functions -----
library(dplyr)
library(xcms)
library(stringr)
library(gtools)
library(MSnbase)
library(purrr)

data_cent <- function(file) {
  data <- readMSData(file,
                     mode = "onDisk") %>%
    pickPeaks()
  
  fls_new <- fileNames(data) %>%
    str_sub(end=-6) %>%
    paste("cent",sep = "_") %>%
    paste("mzML",sep = ".")
  writeMSData(data, file = fls_new)
}

remove_original <- function() {
  fls <- list.files(path = datadir, pattern = "^.*\\.mzML$",
             all.files = TRUE, full.names = TRUE,
             recursive = FALSE, ignore.case = FALSE,
             include.dirs = FALSE, no.. = TRUE)  %>%
    mixedsort(decreasing = F)
  fls <- fls[!str_detect(fls,"cent")]
  file.remove(fls)
}

get_mz_cn <- function(C,N,H,O,pol) {
  # Get the mz ratios of a C and N isotope labeled compound
  # specify the number of c, n, h, o atoms
  # specify polarity for +1 charge or -1 charge
  
  m.c12 <- 12
  m.c13 <- 13.003355
  m.h <- 1.007825
  m.n14 <- 14.003074
  m.n15 <- 15.000109
  m.o16 <- 15.994915
  m.o17 <- 16.999131
  m.o18 <- 17.999159
  m.p <- 1.00727646677 # ex mass of proton
  
  b <- c(m.c12,m.c13,m.n14,m.n15,m.h,m.o16)
  
  a <- matrix(
    c(
      rep(0:C,each=N+1),
      rep(C:0,each=N+1),
      rep(0:N,times=C+1),
      rep(N:0,times=C+1),
      rep(H,times=(C+1)*(N+1)),
      rep(O,times=(C+1)*(N+1))
    ),
    ncol = length(b)
  )
  
  mz <- (a %*% b + pol*m.p) %>%
    matrix(nrow = N+1, ncol = C+1) %>%
    as.data.frame(
      row.names = paste("[15]N",N:0,sep = "")
    )
  
  colnames(mz) <- paste("[13]C",C:0,sep = "")
  
  return(mz)
}

get_ppm_range <- function(x,ppm) {
  return(x*c(1-ppm*1e-6,1+ppm*1e-6))
}

remove_zeros <- function(x) {
  x <- x[!x==0]
  return(x)
}

get_abundance <- function(datadir, ppm, rt_range, C, N, H, O, pol, parallel = FALSE, meanIntensity = TRUE, centroidData = TRUE, unlabeled = NA) {
  # Set up parallel processing
  if (parallel == FALSE) {
    register(SerialParam())
  } else if (parallel == TRUE) { # Not recommended, can have many errors
    if (.Platform$OS.type == "unix") {
      register(bpstart(MulticoreParam(6)))
    } else {
      register(bpstart(SnowParam()))
    }
  }
  
  # File names
  fls <- list.files(path = datadir, pattern = "^.*\\.mzML$",
                    all.files = TRUE, full.names = TRUE,
                    recursive = FALSE, ignore.case = FALSE,
                    include.dirs = FALSE, no.. = TRUE)  %>%
    mixedsort(decreasing = F)
  
  # Centroid data if not already centroided
  if (centroidData == TRUE) { # If using centroided data
    check <- (str_detect(fls,"cent") %>% sum) > 0 # check if centroided data already exist
    
    if (check == TRUE) { # if exist
      fls <- fls[str_detect(fls,"cent")] # read only the centroided data
    } else if (check == FALSE) { # if not exist
      print("Data not centroided yet. Computing centroids")
      data_cent(file = fls) # centroid data and save
      print("Data centroided and saved")
      # Get new file names
      fls <- list.files(path = datadir, pattern = "^.*\\cent.mzML$",
                        all.files = TRUE, full.names = TRUE,
                        recursive = FALSE, ignore.case = FALSE,
                        include.dirs = FALSE, no.. = TRUE)  %>%
        mixedsort(decreasing = F)
    }
  } else if (centroidData == FALSE) { # if not using centroided data
    fls <- fls[!str_detect(fls,"cent")] # get only uncentroided file names
  }

  # Read data
  data <- readMSData(fls, mode = "onDisk") %>%
    filterMsLevel(1)
  print("Data reading complete")
  
  # Get mz values
  mzs <- get_mz_cn(C=C,H=H,N=N,O=O,pol=pol) # mz ratio
  label <- sapply(colnames(mzs), function(x) {
    sapply(rownames(mzs), function(y) {
      paste(x,y,sep="_")
    })
  })
  mzs_vec <- as.vector(as.matrix(mzs))
  print(paste((C+1)*(N+1),"mz values computed"))
  
  # Get chromatogram
  chr <- data %>%
    filterMz(mz = get_ppm_range(x=max(mzs),ppm=ppm)) %>%
    filterRt(rt = 60*rt_range) %>%
    chromatogram(aggregationFun = "max", missing = 0)
  print("Chromatogram reading complete")
  
  print("Start getting intensities")
  # Get intensities
  mz_int <- lapply(seq_along(fls), function(f) {
    print(paste("Processing File",f))
    res <- list()
    
    # Get acquisition number(s)
    if (meanIntensity == TRUE) {
      ints <- intensity(chr[1,f])
      acNum <- ints[ints >= (0.5*max(ints))] %>%
        names %>%
        str_extract("(?<=S)[:digit:]+") %>%
        as.integer
    } else if (meanIntensity == FALSE) {
      acNum <- intensity(chr[1,f]) %>%
        which.max %>%
        names %>%
        str_extract("(?<=S)[:digit:]+") %>%
        as.integer
    }
    
    # Getting spectrum
    sp <- data %>%
      filterFile(f) %>%
      filterAcquisitionNum(n = acNum)
    
    for (i in seq_along(mzs_vec)) {
      x <- mzs_vec[i]
      
      # Matching mzs
      matched <- sp %>%
        filterMz(mz = get_ppm_range(x = x,ppm = ppm))
      
      # Obtain mz intensity pair
      vals <- c(mz(matched) %>% unlist %>% remove_zeros %>% mean,
                intensity(matched) %>% sapply(sum) %>% remove_zeros %>% mean)
      
      # Output result
      if (is.na(sum(vals))) {
        res[[i]] <- c(x,0)
      } else {
        res[[i]] <- vals
      }
    }
    
    return(res)
  })
  
  print("Done processing. Generating output...")
  # Generate output
  out <- sapply(mz_int, function(f) {
    d <- unlist(f)
  }) %>% 
    t %>%
    as.data.frame
  
  rownames(out) <- paste("File",seq_along(fls),sep = "_")
  colnames(out) <- paste(
    rep(label,each = 2),
    rep(c("mz","int"),times = length(label)),
    sep = "_"
  )
  
  if (!is.na(unlabeled) & unlabeled == "C") {
    print("C is unlabeled")
    out <- out %>%
      select(colnames(out)[str_detect(colnames(out),"C0")])
    C <- 0
    label <- label[str_detect(label,"C0")]
  } else if (!is.na(unlabeled) & unlabeled == "N") {
    print("N is unlabeled")
    out <- out %>%
      select(colnames(out)[str_detect(colnames(out),"N0")])
    N <- 0
    label <- label[str_detect(label,"N0")]
  }
  
  out <- out %>%
    mutate(total = out %>%
             select(colnames(out)[str_detect(colnames(out),"int")]) %>%
             apply(.,MARGIN = 1,sum)
           )
  
  if ((unlabeled == "N")|is.na(unlabeled)) {
    print("Computing [13]C abundance")
    out <- out %>%
      mutate(c13Abundance = sapply(seq_along(label), function(x) {
        a <- out %>%
          select(colnames(out)[str_detect(colnames(out),"int")])
        b <- rep(C:0,each = N+1)
        return(a[,x]*b[x])
      }) %>%
        apply(.,MARGIN = 1,sum) %>%
        `/`(C*total)
      ) %>%
      select(c13Abundance,everything())
    print(paste("RSD of [13]C:",rsd(out$c13Abundance)))
  }
  
  if ((unlabeled == "C")|is.na(unlabeled)) {
    print("Computing [15]N abundance")
    out <- out %>%
      mutate(n15Abundance = sapply(seq_along(label), function(x) {
        a <- out %>%
          select(colnames(out)[str_detect(colnames(out),"int")])
        b <- rep(N:0,times = C+1)
        return(a[,x]*b[x])
      }) %>%
        apply(.,MARGIN = 1,sum) %>%
        `/`(N*total)
      ) %>%
      select(n15Abundance,everything())
    print(paste("RSD of [15]N:",rsd(out$n15Abundance)))
  }
  
  if (is.na(unlabeled)) {
    print("Computing total isotope abundance")
    out <- out %>%
      mutate(allAbundance = sapply(seq_along(label), function(x) {
        a <- out %>%
          select(colnames(out)[str_detect(colnames(out),"int")])
        b1 <- rep(C:0,each = N+1)
        b2 <- rep(N:0,times = C+1)
        b <- b1+b2
        return(a[,x]*b[x])
      }) %>%
        apply(.,MARGIN = 1,sum) %>%
        `/` ((N+C)*total)
      ) %>%
      select(allAbundance,everything())
    print(paste("RSD of total isotope:",rsd(out$allAbundance)))
  }
  print("Output generated")
  return(out)
}

rsd <- function(x) {sd(x)/mean(x)} # relative standard deviation

retrieve_ints <- function(out) {
  ints <- out %>% select(colnames(out)[str_detect(colnames(out),"int")])
  return(ints)
}

## Edit parameters -----
out_mean_cent <- get_abundance(
  datadir = "./data/20221012", # Edit data directory
  ppm = 5,
  rt_range = c(2.3,3.0),
  C = 6,
  H = 14,
  N = 2,
  O = 2,
  pol = 1, # -1 for neg mode
  unlabeled = NA, # If "C" or "N" is unlabeled, or NA if C and N are both labeled
  centroidData = TRUE, # use centroided data or not; using centroided data takes the max intensity in a mz dimension cluster; using original data takes the sum intensity in a mz dimension cluster
  parallel = FALSE, # Use biocparallel or not; not recommended
  meanIntensity = TRUE # compute intensities from max spectra or mean of half-width spectrum
)

ints <- retrieve_ints(out = out)
## Utils -----
# calc <- get_mz_cn( # Compute mz values
#   C = 6,
#   H = 14,
#   N = 2,
#   O = 2,
#   pol = 1
# )

# remove_original() # Remove uncentroided original files to save disk space

# write.csv(out,"output.csv") # Write the output to a csv file