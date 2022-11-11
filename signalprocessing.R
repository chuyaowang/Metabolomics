library(gsignal)
library(dplyr)
library(xcms)
library(MSnbase)
library(gtools)
library(ggplot2)

datadir <- "./rdata/20221031/stable/CN-GAN-1ppm-6.mzML"
data <- readMSData(datadir, mode = "onDisk")

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

plot.frequency.spectrum <- function(X.k, xlimits=c(0,length(X.k))) {
  plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))
  
  # TODO: why this scaling is necessary?
  plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2] 
  
  plot(plot.data, t="h", lwd=2, main="", 
       xlab="Frequency (Hz)", ylab="Strength", 
       xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
}

mzs <- get_mz_cn(C=2,H=5,N=1,O=2,pol=1) # mz ratio
label <- sapply(colnames(mzs), function(x) {
  sapply(rownames(mzs), function(y) {
    paste(x,y,sep="_")
  })
})
mzs_vec <- as.vector(as.matrix(mzs))

chr <- data %>%
  filterMz(mz = get_ppm_range(x=mzs_vec[1],ppm=5)) %>%
  filterRt(rt = 60*c(0.5,1.5)) %>%
  chromatogram(aggregationFun = "sum", missing=0)
plot(chr)

# rtime <- data %>%
#   filterMz(mz = get_ppm_range(x=mzs_vec[1],ppm=5)) %>%
#   filterRt(rt = 60*c(0.5,1.5)) %>%
#   fData %>%
#   select(retentionTime) %>%
#   `[`(which(.==max(.) | .==min(.)))

# ints <- data %>%
#   filterAcquisitionNum(as.integer(169)) %>%
#   intensity
# mzs <- data %>%
#   filterAcquisitionNum(as.integer(169)) %>%
#   mz
# vals <- data.frame(ints,mzs)
# colnames(vals) <- c("ints","mz")
# 
# ggplot(vals,aes(x=mz,y=ints)) +
#   geom_line()

ints <- intensity(chr[1,1])
fftint <- fft(ints)
plot.frequency.spectrum(fftint)
