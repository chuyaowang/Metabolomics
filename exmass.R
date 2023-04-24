library(openxlsx)
library(dplyr)
library(stringr)
source("isotopeApp/helpers.R")

compounds <- read.xlsx("data/metabolites_neg.xlsx",sheet=1,startRow = 1)
# compounds <- read.xlsx("data/metabolites_pos.xlsx",sheet=1,startRow = 1)
elements <- c("C","H","N","O","P","S")
count <- matrix(nrow=length(compounds$Formula),ncol = length(elements)) %>% as.data.frame
for (i in seq_along(elements)) {
  element_notexist <- str_extract(compounds$Formula,elements[i]) %>% is.na
  element_number <- str_extract(compounds$Formula,paste(elements[i],"\\d+",sep="")) %>%
    str_extract(.,"\\d+") %>%
    as.numeric()
  element_number[element_notexist] <- 0
  element_number[is.na(element_number)] <- 1
  
  count[,i] <- element_number
}
colnames(count) <- elements

count <- mutate(count, ion_formula = paste("C",C,"H",H-1,"N",N,"O",O,"P",P,"S",S,sep="")) # change to H+1 for positive mode
for (i in seq_along(count$ion_formula)) {
  vec <- str_extract_all(count$ion_formula[i],"\\D?\\d+") %>% unlist 
  vec[str_detect(vec,"^\\D?0$")] <- ""
  vec[str_detect(vec,"^\\D?1$")] <- str_remove(vec[str_detect(vec,"^\\D?1$")],"1")
  count$ion_formula[i] <- str_flatten(vec)
}
compounds <- bind_cols(compounds, count)

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
  a[cpd_idx] <- mzs[length(mzs)] %>% round(.,digits=4) %>% as.character() %>% str_flatten(.,collapse = " ") 
}

compounds <- compounds %>% mutate(exmass_neg = a) #change col name accordingly
write.csv(compounds,"outputs/exmass_neg.csv")
