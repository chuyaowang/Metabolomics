library(openxlsx)
library(dplyr)

compounds <- read.xlsx("data/metabolites_neg.xlsx")

a <- rep("",times=length(compounds$Name))
for (cpd_idx in seq_along(compounds$Name)) {
  row <- compounds[cpd_idx,]
  mzs <- get_mz_ch(C=row$C,
                   N=row$N,
                   H=row$H,
                   S=row$S,
                   O=row$O,
                   P=row$P,
                   c13=row$C,
                   d=0,
                   pol=-1 ) %>% # -1 for neg mode
    unlist
  a[cpd_idx] <- mzs %>% round(.,digits=4) %>% as.character() %>% str_flatten(.,collapse = " ") 
}

compounds <- compounds %>% mutate(exmass_neg = a)
