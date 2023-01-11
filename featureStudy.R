library(dplyr)
library(openxlsx)
library(stringr)
## Data handling -----
# Read data and clean up
data.pos <- read.xlsx("./data/features1.xlsx") 
data.neg <- read.xlsx("./data/features1_neg.xlsx")
data <- list(data.pos,data.neg) %>%
  lapply(.,function(x){
    x <- x %>%
      `[`(,!str_detect(colnames(.),"Blank")) %>%
      `[`(,!str_detect(colnames(.),"QC"))
    
    colnames(x)[str_detect(colnames(x),"Area")] <- colnames(x)[str_detect(colnames(x),"Area")] %>% str_remove(.,".raw.") %>% str_remove(.,"\\(F[:digit:]+\\)") %>% str_remove(.,"Sample_All_FS_") %>% str_replace(.,"\\..","\\_") %>% str_remove(.,"\\.") %>% str_remove(.,"_[:digit:]+")
    
    x <- x %>% 
      mutate(mode = colnames(x)[str_detect(colnames(x),"Area")] %>%
               str_sub(-3) %>%
               unique %>%
               rep(.,nrow(x)))
    
    colnames(x)[str_detect(colnames(x),"Area")] <- colnames(x)[str_detect(colnames(x),"Area")] %>% str_sub(.,end=-5)

    return(x)
  })

# Merge pos and neg data
common.colname <- data %>% lapply(.,colnames) %>% unlist %>% table %>% as.data.frame() %>% `[`(.$Freq>1,) %>% `[`(,1) %>% as.vector()

data <- data %>% lapply(.,function(x){
  x %>% select(all_of(common.colname))
}) %>%
  bind_rows
  
# Remove duplicates
duped <- data$Formula %>% table %>% as.data.frame() %>% `[`(.$Freq>1,) %>% `[`(,1) %>% as.vector()
for (i in duped) {
  idx <- which(data$Formula %in% i)
  newrow <- data[c(idx),] %>%
    select(colnames(data)[str_detect(colnames(data),"Area")]) %>%
    sapply(.,max)
  for (j in idx) {
    data[j,which(str_detect(colnames(data),"Area"))] <- newrow
  }
  data$Name[idx[-1]] <- "dup"
}
data <- data %>% filter((Name != "dup")|(is.na(Name)))

# Give NA names an ID
data$Name[which(is.na(data$Name))] <- paste("Unknown",seq_along(which(is.na(data$Name))),sep="_")

# Transpose and output
data<- data[which(!str_detect(data$Name,"Unknown")),]
areas.out <- data %>% select(colnames(data)[str_detect(colnames(data),"Area")]) %>% t %>% scale %>% as.data.frame
colnames(areas.out) <- data$Name
rownames(areas.out) <- rownames(areas.out) %>% str_sub(.,start=11)
areas.out <- areas.out %>%
  mutate(grp = rownames(areas.out) %>%
           str_remove(.,"\\-[:digit:]")) %>%
  select(grp, everything())
areas.out <- areas.out[2:nrow(areas.out),]
dummy.row <- c("768O",sapply(areas.out %>% filter(grp=="768O") %>% select(!grp),mean))
areas.out <- rbind(areas.out,dummy.row)
write.csv(areas.out,"knockout.csv")

## PCA -----
library(mixOmics)
library(viridis)
res.pca <- pca(areas.out[,-1],center=F,scale=F,ncomp=5)
plotIndiv(res.pca,
          ind.names = T,
          legend = T,
          ellipse = T,
          star = F,
          group = areas.out$grp,
          col.per.group = viridis(length(areas.out$grp %>% unique),begin=0.4,end=1),
          title = "Scores Plot",
          comp = c(1,3),
          style = "ggplot2",
          cex = 3,
          lwd = 0.1)

## PLSDA -----
library(mixOmics)
library(ggrepel)
data_noSH3 <- data[!str_detect(rownames(data),"SH3"),]
class_noSH3 <- as.character(class)[!str_detect(class,"SH3")] %>% as.factor

res.plsda <- splsda(data,class,keepX = c(35,35)) 
res.plsda.noSH3 <- splsda(data_noSH3,class_noSH3, ncomp = 3, keepX = c(35,35,35))

background <- background.predict(res.plsda.noSH3,comp.predicted = 2,dist = "max.dist",resolution = 200)
plotIndiv(res.plsda.noSH3, 
          ind.names = FALSE, # Whether to display sample names 
          legend=TRUE,
          ellipse = TRUE, # Draw confidence ellipse
          star = F, # Whether to draw lines from the centroid
          comp = c(1,2),
          title = 'sPLS-DA on Renal Cancer Cells',
          background = background)

plotVar(res.plsda.noSH3, 
        cutoff = 0.7,
        cex = 2,
        style = "ggplot2",
        var.names = TRUE)

plotLoadings(res.plsda.noSH3, 
             contrib = 'max', 
             method = 'mean',
             comp = 1,
             title = "Contribution on Comp 1",
             size.title = rel(1.5))

legend <- list(legend = levels(class_noSH3), # set of classes
            col = unique(color.mixo(class_noSH3)), # set of colours
            title = "Cell Type", # legend title
            cex = 0.7)
cim <- cim(res.plsda.noSH3, 
           row.sideColors = color.mixo(class_noSH3), 
           legend = legend,
           dist.method = c("correlation","correlation"),
           clust.method = c("ward","ward"))



## Data Transformation -----
# Standard scaling
d2 <- scale(data_noSH3) %>%
  as.data.frame


## ANOVA -----
res.anova <- lapply(d2,function(x){
  aov(x ~ class_noSH3, data = d2)
})
a <- sapply(res.anova,function(x){
  pval <- summary(x)[[1]][1,5]
}) %>%
  `>`(0)
