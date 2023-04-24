library(dplyr)
library(ggplot2)
library(openxlsx)
library(KEGGREST)
library(stringr)
library(pathview)

. <- "data/20230413/msmls-lot-230-13.xlsx"
keggid <- read.xlsx(.,startRow=2)
# keggid <- keggid %>% 
#   mutate(KEGG=cpdidmap(CAS,"CAS Registry Number","KEGG")) # Converts CAS to KEGG

# keggquery <- keggid %>%
#   mutate(queryTerm = paste("compound:",KEGG,sep="")) %>%
#   select(queryTerm) %>%
#   filter(!str_detect(queryTerm,"NA")) %>%
#   apply(.,1,keggGet)
keggid <- keggid %>%
  mutate(BRITE = rep("",times=nrow(keggid)))
for (i in 554:nrow(keggid)) {
  id <- keggid[i,]$KEGG.ID_CSID %>% str_trim
  check <- !is.na(id)&!str_detect(id,"CSID")&!str_detect(id,"CID")
  if (!check) {
    next
  }
  keggquery <- keggGet(id)
  if ("BRITE" %in% names(keggquery[[1]])) {
    keggid[i,]$BRITE <- keggquery[[1]]$BRITE[3] %>% str_trim()
  }
}

p_mass <- 1.00727646677 # exact mass of a proton
keggid <- keggid %>%  
  mutate(name = keggquery %>%
           sapply(.,function(x){x %>% 
               `[[`(1) %>% # compound
               `[[`(2) %>% # name
               `[[`(1)})) %>% # first one
  mutate(exMass = keggquery %>%
           sapply(.,function(x){x %>%
               `[[`(1) %>% # compound
               `[[`(4)}) %>% # exact mass
           as.numeric) %>% # as number
  mutate(negMass = exMass - p_mass) %>% # Add proton
  mutate(posMass = exMass + p_mass) # Delete proton

write.csv(keggid,"keggid.csv")
