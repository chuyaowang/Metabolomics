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

read_data <- function(fls, centroidData){
  # Check if fls contain centroided files
  check <- str_detect(fls,"cent") %>% sum %>% `>`(0)
  if (check & centroidData) {
    data <- readMSData(fls, mode = "onDisk")
  } else if (!check & centroidData) {
    data <- readMSData(fls, mode = "onDisk") %>%
      pickPeaks()
    fls_new <- fileNames(data) %>%
      str_sub(end=-6) %>%
      paste("cent",sep = "_") %>%
      paste("mzML",sep = ".")
    writeMSData(data, file = fls_new)
    data <- readMSData(fls_new %>% mixedsort(decreasing = T), mode = "onDisk")
  } else if (!centroidData) {
    data <- readMSData(fls, mode = "onDisk")
  }
  return(data)
}

get_abundance <- function(data, ppm, rt_range, bg_range, mzs, multiplier, background, unlabeled) {
  # Set up serial processing
  register(SerialParam())
  
  # Make mz vector and labels
  label <- sapply(colnames(mzs), function(x) {
    sapply(rownames(mzs), function(y) {
      paste(x,y,sep="_")
    })
  })
  mzs_vec <- as.vector(as.matrix(mzs))
  
  # Get intensities -----
  mz_int <- future_lapply(seq_along(fileNames(data)), function(f) {
    res <- list()
    
    # Filtering for peak range-----
    data_peak <- data %>%
      filterFile(f) %>%
      filterRt(rt = rt_range)
    
    # Filtering for background range-----
    if (background) {
      data_bg <- data %>%
        filterFile(f) %>%
        filterRt(rt = bg_range)
    }

    
    for (i in seq_along(mzs_vec)) {
      x <- mzs_vec[i]
      
      # Filtering for mz range and extract spectra
      data_peak_mz <- data_peak %>%
        filterMz(mz = get_ppm_range(x = x,ppm = ppm)) %>%
        spectra
      
      if (background) {
        data_bg_mz <- data_bg %>%
          filterMz(mz = get_ppm_range(x = x,ppm = ppm)) %>%
          spectra
      }
      
      # Get acquisition numbers
      if (x == max(mzs_vec)) {
        ints <- lapply(data_peak_mz,function(x){x@intensity}) %>% sapply(sum)
        acNum <- which(ints >= multiplier*max(ints))
      }
      
      # Get mz intensity pairs
      val.mz <- lapply(data_peak_mz, function(x) {x@mz}) %>% `[`(acNum) %>% unlist %>% remove_zeros %>% median
      val.intensity <- lapply(data_peak_mz, function(x) {x@intensity}) %>% `[`(acNum) %>% sapply(sum) %>% remove_zeros %>% median
      vals <- c(val.mz,val.intensity)

      # Background subtraction -----
      if (background) {
        bg.intensity <- lapply(data_bg_mz, function(x) {x@intensity}) %>%
          sapply(sum) %>%
          remove_zeros %>%
          median
        
        if (is.na(bg.intensity)) {bg.intensity <- 0}
        
        vals[2] <- vals[2] - bg.intensity
        if (is.na(vals[2])) {
          vals[2] <- 0
        } else if (vals[2]<0) {
          vals[2] <- 0
        }
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
  
  # Generate output
  out <- sapply(mz_int, function(f) {
    d <- unlist(f)
  }) %>% 
    t %>%
    as.data.frame
  
  rownames(out) <- fileNames(data) %>%
    basename %>%
    sub(pattern = "\\.mzX?ML", replacement = "")
  
  colnames(out) <- paste(
    rep(label,each = 2),
    rep(c("mz","int"),times = length(label)),
    sep = "_"
  )
  
  C <- ncol(mzs) - 1
  N <- nrow(mzs) - 1
  
  if (!is.na(unlabeled) & unlabeled == "C") {
    C <- 0
  } else if (!is.na(unlabeled) & unlabeled == "N") {
    N <- 0
  }
  
  out <- out %>%
    mutate(total = out %>%
             select(colnames(out)[str_detect(colnames(out),"int")]) %>%
             apply(.,MARGIN = 1,sum)
    )
  
  if ((unlabeled == "N")|is.na(unlabeled)) {
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
  }
  
  if ((unlabeled == "C")|is.na(unlabeled)) {
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
  }
  
  if (is.na(unlabeled)) {
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
  }
  return(out)
}

rsd <- function(x) {sd(x)/mean(x)} # relative standard deviation

retrieve_ints <- function(out) {
  ints <- out %>% select(colnames(out)[str_detect(colnames(out),"int")|str_detect(colnames(out),'Abundance')])
  return(ints)
}