## libs and functions -----
library(dplyr)
library(xcms)
library(stringr)
library(gtools)
library(MSnbase)
library(purrr)
library(future.apply)

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

get_abundance <- function(datadir, ppm, rt_range, C, N, H, O, pol, parallel = FALSE, multiplier = 0.5, centroidData = TRUE, unlabeled = NA, backgroundRange = 1) {
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
    mixedsort(decreasing = T)
  
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

  # Read data -----
  data <- readMSData(fls, mode = "onDisk")
  print("Data reading complete")
  
  # Get mz values -----
  mzs <- get_mz_cn(C=C,H=H,N=N,O=O,pol=pol) # mz ratio
  label <- sapply(colnames(mzs), function(x) {
    sapply(rownames(mzs), function(y) {
      paste(x,y,sep="_")
    })
  })
  mzs_vec <- as.vector(as.matrix(mzs))
  print(paste((C+1)*(N+1),"mz values computed"))
  
  # Get chromatogram -----
  rt_range <- 60*rt_range
  chr <- data %>%
    filterMz(mz = get_ppm_range(x=max(mzs),ppm=ppm)) %>%
    filterRt(rt = rt_range) %>%
    chromatogram(aggregationFun = "sum", missing = 0)
  print("Chromatogram reading complete")

  # Get intensities -----
  mz_int <- future_lapply(seq_along(fls), function(f) {
    print(paste("Processing File",f))
    res <- list()
    
    # Get acquisition number(s) -----
    if (as.integer(multiplier) != 1) {
      ints <- intensity(chr[1,f])
      acNum <- ints[ints >= (multiplier*max(ints))] %>%
        names %>%
        str_extract("(?<=S)[:digit:]+") %>%
        as.integer
    } else if (as.integer(multiplier) == 1) {
      acNum <- intensity(chr[1,f]) %>%
        which.max %>%
        names %>%
        str_extract("(?<=S)[:digit:]+") %>%
        as.integer
    }
    
    # Getting analyte spectrum -----
    data_ana <- data %>%
      filterFile(f) %>%
      filterAcquisitionNum(n = acNum)
    
    # Getting control spectrum -----
    ctl_range <- data %>%
      filterFile(f) %>%
      fData %>%
      select(retentionTime) %>%
      unlist
    ctl_range <- c(min(ctl_range),max(ctl_range))
    ctl_range1 <- c(max(ctl_range[1],rt_range[1]-backgroundRange*60),rt_range[1])
    ctl_range2 <- c(rt_range[2],min(ctl_range[2],rt_range[2]+backgroundRange*60))
    
    sprintf("Ranges for background substraction: %.3f to %.3f seconds and %.3f to %.3f seconds",ctl_range1[1],ctl_range1[2],ctl_range2[1],ctl_range2[2])
    
    data_ctl1 <- data %>%
      filterFile(f) %>%
      filterRt(rt = ctl_range1)
    data_ctl2 <- data %>%
      filterFile(f) %>%
      filterRt(rt = ctl_range2)
    
    for (i in seq_along(mzs_vec)) {
      x <- mzs_vec[i]
      
      # Matching mzs
      matched <- data_ana %>%
        filterMz(mz = get_ppm_range(x = x,ppm = ppm))
      
      # Obtain mz intensity pair
      vals <- c(mz(matched) %>% unlist %>% remove_zeros %>% median,
                intensity(matched) %>% sapply(sum) %>% remove_zeros %>% median)
      
      # Background subtraction -----
      data_ctl1_med <- data_ctl1 %>%
        filterMz(mz = get_ppm_range(x=x,ppm=ppm)) %>%
        intensity %>%
        sapply(sum) %>%
        remove_zeros %>%
        median
      # data_ctl2_med <- data_ctl2 %>%
      #   filterMz(mz = get_ppm_range(x=x,ppm=ppm)) %>%
      #   intensity %>%
      #   sapply(sum) %>%
      #   remove_zeros %>%
      #   median
      data_ctl2_med <- data_ctl1_med
      
      if (is.na(data_ctl1_med)) {
        med <- data_ctl2_med
      } else if (is.na(data_ctl2_med)) {
        med <- data_ctl1_med
      } else {
        med <- median(c(data_ctl1_med,data_ctl2_med))
      }
      
      if (is.na(med)) {med <- 0}
      
      vals[2] <- vals[2] - med
      if (is.na(vals[2])) {
        vals[2] <- 0
      } else if (vals[2]<0) {
        vals[2] <- 0
      }
      
      # Output result -----
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
  
  # rownames(out) <- paste("File",seq_along(fls),sep = "_")
  rownames(out) <- fls %>%
    str_split("/") %>%
    as.data.frame %>%
    t %>%
    `[`(,ncol(.))
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
  ints <- out %>% select(colnames(out)[str_detect(colnames(out),"int")|str_detect(colnames(out),'Abundance')])
  return(ints)
}

## Edit parameters -----
plan(multisession, workers = 4)

out_arg_max_cent <- get_abundance(
  datadir = "./rdata/20221108/", # Edit data directory
  ppm = 5,
  rt_range = c(.4,2.5), # Retention time range for the ENTIRE PEAK
  C = 6,
  H = 14,
  N = 4,
  O = 2,
  pol = 1, # -1 for neg mode
  unlabeled = "C", # If "C" or "N" is unlabeled, or NA if C and N are both labeled
  centroidData = TRUE, # use centroided data or not; using centroided data takes the max intensity in a mz dimension peak; using original data takes the sum intensity in a mz dimension peak
  parallel = FALSE, # Use biocparallel or not; not recommended
  multiplier = 0.5, # Width around the max spectra to look for intensities; 1 is the max spectrum, 0.5 is half peak intensities spectrum, etc.
  backgroundRange = 1 # Minutes before and after the peak to be used for background subtraction
  )

# gly: 0.5 to 1.5
# urea: 0.8 to 2
# ints <- retrieve_ints(out = out_mean_cent)
## Utils -----
# calc <- get_mz_cn( # Compute mz values
#   C = 6,
#   H = 14,
#   N = 2,
#   O = 2,
#   pol = 1
# )

# remove_original() # Remove uncentroided original files to save disk space
# write.csv(out_stable,"./rdata/20221028/outputganansuan8zhen.csv") # Write the output to a csv file

## One-timers -----
# out_conc <- list.dirs("rdata/20221031/conc")[-1] %>%
#   mixedsort(decreasing = F) %>%
#   lapply(.,function(x){
#     a <- get_abundance(
#       datadir = x, # Edit data directory
#       ppm = 5,
#       rt_range = c(.5,1.5), # Retention time range for the ENTIRE PEAK
#       C = 2,
#       H = 5,
#       N = 1,
#       O = 2,
#       pol = 1, # -1 for neg mode
#       unlabeled = NA, # If "C" or "N" is unlabeled, or NA if C and N are both labeled
#       centroidData = TRUE, # use centroided data or not; using centroided data takes the max intensity in a mz dimension peak; using original data takes the sum intensity in a mz dimension peak
#       parallel = FALSE, # Use biocparallel or not; not recommended
#       multiplier = .5, # Width around the max spectra to look for intensities; 1 is the max spectrum, 0.5 is half peak intensities spectrum, etc.
#       backgroundRange = 2 # Minutes before and after the peak to be used for background subtraction
#     )
#     a <- a %>%
#       retrieve_ints() %>%
#       select(colnames(.)[str_detect(colnames(.),"int")],everything())
# 
#     write.csv(a,paste(x,"csv",sep="."))
#     return(a)
#   })
# 
# x <- rep(c(100,300,500,800,1000,2000),each = 3)
# y <- sapply(out_conc, function(x){
#   return(x %>% select(allAbundance) %>% unlist)
# }) %>% as.vector
# library(ggplot2)
# library(hrbrthemes)
# data <- as.data.frame(list(x,y))
# colnames(data) <- c("Conc","Abundance")
# 
# model <- lm(Abundance ~ log(Conc),data)
# data <- mutate(data,newy = predict(model,data))
# 
# ggplot(data, aes(x=Conc,y=Abundance)) +
#   geom_point(aes(x = Conc, y = Abundance, color = factor(Conc)), size = 1.5, show.legend = FALSE) +
#   ggtitle("C/N Abundance vs. Concentration") +
#   xlab("Concentration (ppb)") +
#   ylab("Abundance") +
#   scale_x_continuous(n.breaks = 10) +
#   # geom_label(aes(label=round(Abundance,digits=4), size = NULL), nudge_y = 0.006) +
#   geom_line(aes(x=Conc,y=newy),linetype = "dashed", color = "blue") +
#   theme_light()
# 
# ggplot(data,aes(x=factor(Conc),y=Abundance)) +
#   geom_boxplot()
# out_mean_uncent_urea %>% 
#   retrieve_ints() %>% 
#   select(colnames(.)[str_detect(colnames(.),"int")],everything()) %>%
#   write.csv("./rdata/20221024/UREA/urea.csv")

