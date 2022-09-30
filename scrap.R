# xs <- xcmsSet(
#   files=file, 
#   method="centWaveWithPredictedIsotopeROIs",
#   ppm=10,
#   mslevel = 1,
#   snthresh = 5,
#   integrate = 1,
#   peakwidth = c(2,30),
#   prefilter = c(1,1e4),
#   mzCenterFun = "wMeanApex3",
#   mzdiff = -0.001
# )
# 
# an <- xsAnnotate(xs, polarity="negative") # constructor; extracts peak table
# an <- groupFWHM(an, perfwhm = 1) # group peaks by retention time
# an <- findIsotopesWithValidation(object = an, ppm = 10,
#                                  mzabs = 0.01, intval="intb",
#                                  maxcharge = 3) # annotate isotopic peaks
# ## extract annotated peak table
# peakTable <- getPeaklist(an) # extract peak list