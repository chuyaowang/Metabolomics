library(dplyr)
library(ggplot2)
library(openxlsx)
library(KEGGREST)

. <- "./data/keggid.xlsx"
keggid <- read.xlsx(.)

p_mass <- 1.00727646677 # exact mass of a proton

keggquery <- keggid %>%
  mutate(queryTerm = paste("compound:",KEGG,sep="")) %>%
  select(queryTerm) %>%
  apply(.,1,keggGet)

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