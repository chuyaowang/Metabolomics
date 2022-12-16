## Data handling -----
library(dplyr)
library(openxlsx)
library(stringr)

data <- read.xlsx("./data/features.xlsx") 
data <- data %>%
  select(Name, colnames(data) %>% 
           `[`(str_detect(.,"Norm..Area")) %>%
           `[`(str_detect(.,"Sample"))) 
data <- data %>%
  `[`(!is.na(data[,2]),)
rownames(data) <- data$Name
data <- data %>% 
  select(!Name)
colnames(data) <- colnames(data) %>%
  str_extract_all(.,"(?<=:.)[:graph:]+(?=_Sample)")
data <- t(data) %>%
  as.data.frame
class <- rownames(data) %>%
           str_extract(.,"(?<=)[:alnum:]+(?=-)")
class[3] <- "SH3"
class <- as.factor(class)

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
