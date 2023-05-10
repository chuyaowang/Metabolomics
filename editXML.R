library(XML)
library(xml2)
library(dplyr)
library(stringr)
library(openxlsx)
library(methods)
library(readr)

# Prepare compound list
info <- read.xlsx("data/20230413/msmls-lot-230-13.xlsx",startRow=2)
cpds <- read_csv("data/20230413/tca_sigma.csv",skip=1)

idx <- match(str_to_upper(cpds$`Sample Name`),info$PRIMARY_NAME)
cpds$`Sample Wt` <- rep("",nrow(cpds))
cpds$`Sample ID` <- info[idx,] %>% select("MOLECULAR_FORMULA") %>% unlist %>% str_remove(.,"\\([:alpha:]+\\)") %>% str_trim(.,"both") 
write_csv(cpds, "data/20230413/tca_sigma_addFormula.csv") # Add the first row back!

# Read data files
data <- xmlToDataFrame("data/20230413/sugar_individual.xml")

# Add info by compound name
idx <- match(str_to_upper(data$Compound),info$PRIMARY_NAME)
data$CasId <- info[idx,] %>% select(CAS.ID) %>% unlist
data$HMDBId <- info[idx,] %>% select("HMDB/YMDB.ID") %>% unlist
data$KEGGId <- info[idx,] %>% select("KEGG.ID_CSID") %>% unlist
data$PubChemId <- info[idx,] %>% select(PC_CID) %>% unlist
data$SmilesDescription <- info[idx,] %>% select(SMILES) %>% unlist
data$ChemSpiderId <- rep("",times=nrow(data))

# Convert to xml
# Create new xml file
newFile <- newXMLDoc(isHTML = FALSE)
# Create new node
array_node <- newXMLNode("ArrayOfCompoundInfo",doc=newFile)
# row data
row_data <- apply(data, 1, function(x) {
  z1 <- newXMLNode('CompoundInfo') # create a new node for each row
  addChildren(z1, lapply(names(x), function(y) newXMLNode(y, x[y])))
})

# add row data to table node
xmlParent(row_data) <- array_node

# save as xml file
saveXML(newFile, file = "data/20230413/sugar_individual_edit.xml")

# Now go to excel to finish editing
# First copy xml map from original file
# Then copy paste in the new info
