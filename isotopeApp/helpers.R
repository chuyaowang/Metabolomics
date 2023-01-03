get_mz_cn <- function(C,N,H,O,S,c13,n15,pol) {
  # Get the mz ratios of a C and N isotope labeled compound
  # specify the number of c, n, h, o, s atoms
  # specify polarity for +1 charge or -1 charge
  
  m.c12 <- 12.000000
  m.c13 <- 13.003355
  m.h <- 1.007825
  m.n14 <- 14.003074
  m.n15 <- 15.000109
  m.o16 <- 15.994915
  m.o17 <- 16.999131
  m.o18 <- 17.999159
  m.s32 <- 31.972071
  m.p <- 1.00727646677 # ex mass of proton
  
  b <- c(m.c12,m.c13,m.n14,m.n15,m.h,m.o16,m.s32)
  
  a <- matrix(
    c(
      rep(0:C,each=N+1),
      rep(C:0,each=N+1),
      rep(0:N,times=C+1),
      rep(N:0,times=C+1),
      rep(H,times=(C+1)*(N+1)),
      rep(O,times=(C+1)*(N+1)),
      rep(S,times=(C+1)*(N+1))
    ),
    ncol = length(b)
  )
  
  mz <- (a %*% b + pol*m.p) %>%
    matrix(nrow = N+1, ncol = C+1) %>%
    as.data.frame(
      row.names = paste("[15]N",N:0,sep = "")
    )
  
  colnames(mz) <- paste("[13]C",C:0,sep = "")
  
  mz <- mz[(N-n15+1):nrow(mz),(C-c13+1):ncol(mz)]
  return(mz)
}

get_mz_ch <- function(C,N,H,O,S,c13,d,pol) {
  m.c12 <- 12.000000
  m.c13 <- 13.003355
  m.h <- 1.007825
  m.d <- 2.014102
  m.n14 <- 14.003074
  m.o16 <- 15.994915
  m.s32 <- 31.972071
  m.p <- 1.00727646677 # ex mass of proton
  
  b <- c(m.c12,m.c13,m.h,m.d,m.n14,m.o16,m.s32)
  
  a <- matrix(
    c(
      rep(0:C,each=H+1),
      rep(C:0,each=H+1),
      rep(0:H,times=C+1),
      rep(H:0,times=C+1),
      rep(N,times=(C+1)*(H+1)),
      rep(O,times=(C+1)*(H+1)),
      rep(S,times=(C+1)*(H+1))
    ),
    ncol = length(b)
  )
  
  mz <- (a %*% b + pol*m.p) %>%
    matrix(nrow = H+1, ncol = C+1) %>%
    as.data.frame(
      row.names = paste("D",H:0,sep = "")
    )
  
  colnames(mz) <- paste("[13]C",C:0,sep = "")
  
  mz <- mz[(H-d+1):nrow(mz),(C-c13+1):ncol(mz)]
  return(mz)
}

get_mz_nh <- function(C,N,H,O,S,n15,d,pol) {
  m.c12 <- 12.000000
  m.h <- 1.007825
  m.d <- 2.014102
  m.n14 <- 14.003074
  m.n15 <- 15.000109
  m.o16 <- 15.994915
  m.s32 <- 31.972071
  m.p <- 1.00727646677 # ex mass of proton
  
  b <- c(m.n14,m.n15,m.h,m.d,m.c12,m.o16,m.s32)
  
  a <- matrix(
    c(
      rep(0:N,each=H+1),
      rep(N:0,each=H+1),
      rep(0:H,times=N+1),
      rep(H:0,times=N+1),
      rep(C,times=(N+1)*(H+1)),
      rep(O,times=(N+1)*(H+1)),
      rep(S,times=(N+1)*(H+1))
    ),
    ncol = length(b)
  )
  
  mz <- (a %*% b + pol*m.p) %>%
    matrix(nrow = H+1, ncol = N+1) %>%
    as.data.frame(
      row.names = paste("D",H:0,sep = "")
    )
  
  colnames(mz) <- paste("[15]N",N:0,sep = "")
  
  mz <- mz[(H-d+1):nrow(mz),(N-n15+1):ncol(mz)]
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
  # Either [13]C or [15]N
  label1 <- colnames(mzs) %>% str_remove(.,"(?<!\\[)[:digit:]+(?!\\])") %>% unique
  # Either D or [15]N
  label2 <- rownames(mzs) %>% str_remove(.,"(?<!\\[)[:digit:]+(?!\\])") %>% unique
  
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
  
  # numLab1: number of [13]C or [15]N
  # numLab2: number of [15]N or D
  numLab1 <- ncol(mzs) - 1
  numLab2 <- nrow(mzs) - 1
  
  # unlabeled is [13]C or D or [15]N or NA
  if (!is.na(unlabeled)) {
    if (label1 == unlabeled) {
      numLab1 <- 0
    } else if (label2 == unlabeled) {
      numLab2 <- 0
    }
  }
  
  out <- out %>%
    mutate(total = out %>%
             select(colnames(out)[str_detect(colnames(out),"int")]) %>%
             apply(.,MARGIN = 1,sum)
    )
  
  if (numLab1 != 0) {
    out <- out %>%
      mutate(lab1Abundance = sapply(seq_along(label), function(x) {
        a <- out %>%
          select(colnames(out)[str_detect(colnames(out),"int")])
        b <- rep(numLab1:0,each = numLab2+1)
        return(a[,x]*b[x])
      }) %>%
        apply(.,MARGIN = 1,sum) %>%
        `/`(numLab1*total)
      ) %>%
      select(lab1Abundance,everything())
    colnames(out)[1] <- paste(label1,"Abundance") 
  }
  
  if (numLab2 != 0) {
    out <- out %>%
      mutate(lab2Abundance = sapply(seq_along(label), function(x) {
        a <- out %>%
          select(colnames(out)[str_detect(colnames(out),"int")])
        b <- rep(numLab2:0,times = numLab1+1)
        return(a[,x]*b[x])
      }) %>%
        apply(.,MARGIN = 1,sum) %>%
        `/`(numLab2*total)
      ) %>%
      select(lab2Abundance,everything())
    colnames(out)[1] <- paste(label2,"Abundance") 
  }
  
  if (numLab1 != 0 & numLab2 != 0) {
    out <- out %>%
      mutate(allAbundance = sapply(seq_along(label), function(x) {
        a <- out %>%
          select(colnames(out)[str_detect(colnames(out),"int")])
        b1 <- rep(numLab1:0,each = numLab2+1)
        b2 <- rep(numLab2:0,times = numLab1+1)
        b <- b1+b2
        return(a[,x]*b[x])
      }) %>%
        apply(.,MARGIN = 1,sum) %>%
        `/` ((numLab2+numLab1)*total)
      ) %>%
      select(allAbundance,everything())
    colnames(out)[1] <- paste("All","Abundance") 
  }
  return(out)
}

rsd <- function(x) {sd(x)/mean(x)} # relative standard deviation

retrieve_ints <- function(out) {
  ints <- out %>% select(colnames(out)[str_detect(colnames(out),"int")|str_detect(colnames(out),'Abundance')])
  return(ints)
}