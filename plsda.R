cat("\014")
rm(list = ls())

library(mixOmics)
data(srbct)

## sparse PLS-DA example -----
X <- srbct$gene # data matrix
Y <- srbct$class # class membership for each sample
summary(Y) ## class summary

# Run plsda with a subset of 100 variables, 50 for each component
# splsda performs feature selection and classification in one step
MyResult.splsda <- splsda(X, Y, keepX = c(50,50)) 
plotIndiv(MyResult.splsda)

# Select variables on component 1
selectVar(MyResult.splsda, comp=1)$name 

# Non sparse plsda without variable selection
# Non sparse algorithm uses too many parameters that are estimated with a high
# variance. Therefore, it is better to reduce the number of variables used by
# introducing sparsity -> improves clustering
MyResult.plsda <- plsda(X,Y) # 1 Run the method
plotIndiv(MyResult.plsda) 

# plot the correlation circle, selecting only the variables length>0.7
plotVar(MyResult.plsda, 
        cutoff = 0.7) # For non-sparse plsda, a cutoff can be selected to only
                      # show the most relevant variables (high loadings)

# More options for plotting
plotIndiv(MyResult.splsda, 
          ind.names = FALSE, # Whether to display sample names 
          legend=TRUE,
          ellipse = TRUE, # Draw confidence ellipse
          star = F, # Whether to draw lines from the centroid
          title = 'sPLS-DA on SRBCT',
          X.label = 'PLS-DA 1', 
          Y.label = 'PLS-DA 2')

plotVar(MyResult.splsda, 
        var.names=FALSE) # Whether to display variable names

## Other plots -----
# Background prediction
# Classifies each point in the coordinate system into one of the groups
# Size of white line depends on resolution
# background <- background.predict(MyResult.splsda, 
#                                  comp.predicted=2,
#                                  dist = "max.dist",
#                                  resolution = 200) 
# plotIndiv(MyResult.splsda, 
#           comp = 1:2, 
#           group = srbct$class,
#           ind.names = FALSE, 
#           title = "Maximum distance",
#           legend = TRUE,  
#           background = background)

# ROC curve
# area under the curve the more the merrier, each dot represents a threshold
# To find the best threshold, If sensitivity and specificity have the same 
# importance to you, one way of calculating the cut-off is choosing that value 
# that minimizes the Euclidean distance between your ROC curve and the upper left
# corner of your graph.
# Another way is using the value that maximizes (sensitivity + specificity - 1) 
# as a cut-off (Youden's index).

# However, ROC and AUC criteria may not be particularly insightful, or be in 
# agreement with the PLSDA performance, as the prediction threshold in PLS-DA is 
# based on specified distance. AUROC curves use a cutoff that maximises specificity
# and sensitivity rather than this distance and hence should be used a merely a 
# complementary tool.

auc.plsda <- auroc(MyResult.splsda)

## Variable selection -----
MyResult.splsda2 <- splsda(X,Y, ncomp=3, keepX=c(15,10,5))
selectVar(MyResult.splsda2, comp=3)$value

# Plot the highest contributing variables in comp 1 and color them with the group
# in which they have the highest abundance
plotLoadings(MyResult.splsda2, contrib = 'max', method = 'mean')

# A 3D plot
plotIndiv(MyResult.splsda2, style="3d")

## Tuning parameters -----
# Choose the optimal number of parameters
# Run a non-sparse plsda

MyResult.plsda2 <- plsda(X,Y, ncomp=10)

# Get predicative capability for each number of components
# k fold: partition sample into k partitions. Each time one partition is used as
# the validation set until all partitions are used up

set.seed(30) # for reproducibility in this vignette, otherwise increase nrepeat
MyPerf.plsda <- perf(MyResult.plsda2, 
                     validation = "Mfold", 
                     folds = 3, 
                     progressBar = FALSE, 
                     nrepeat = 10) # we suggest nrepeat = 50

plot(MyPerf.plsda, 
     col = color.mixo(5:7), 
     sd = TRUE, 
     legend.position = "horizontal")

# Many performance measures can be viewed
MyPerf.plsda

# Choose the optimal number of variables to keep
# Run 3 fold cross validation for 10 times with a maximum prediction distance
# max.dist: assigns the new sample to the class for which it has the highest
# probability value
# centroids.dist: assigns to the group with the shortest distance to centroid
# from the data point; more robust than max.dist; suited for moderately clustered
# classes
# mahalanobis.dist: similar to centroid distance but uses Mahalanobis distance.
# Mahalanobis distance takes into account the correlation between each of the 
# components, giving more weight to less correlated components; can give non
# linear boundaries.

# Balanced error rate (BER): calculates the error metric balancing for uneven sample
# sizes

list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX # to output the grid of values tested

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(X, 
                                 Y, 
                                 ncomp = 3, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, 
                                 dist = 'max.dist', 
                                 progressBar = FALSE,
                                 measure = "BER", 
                                 test.keepX = list.keepX,
                                 nrepeat = 10,  # we suggest nrepeat = 50
                                 cpus = 12)  

error <- tune.splsda.srbct$error.rate
# optimal number of components based on t-tests on the error rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp 
ncomp
# Optimal number of variables to keep
# This can differ from the results from perf since variable selection is done now
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

# Run splsda based on tuned results
MyResult.splsda.final <- splsda(X, 
                                Y, 
                                ncomp = ncomp, 
                                keepX = select.keepX)

plotIndiv(MyResult.splsda.final, 
          ind.names = FALSE, 
          legend=TRUE,
          ellipse = TRUE, 
          title="sPLS-DA - final result")
# Save a higher resolution plot
ggsave(filename = "test2.tiff",
       width = 5,
       height = 5,
       units = "in",
       dpi = 300)

## Performance assessment for the final model -----
perf.splsda <- perf(MyResult.splsda.final,
                    validation = "Mfold",
                    folds = 3,
                    nrepeat = 50,
                    progressBar = F,
                    dist = "max.dist",
                    auc = T,
                    cpus = 12)

plot(perf.splsda, 
     col = color.mixo(5), 
     sd = TRUE,
     legend.position = "horizontal")
perf.splsda$choice.ncomp

# Do the final plotting
plotIndiv(MyResult.splsda.final,
          comp = c(1,2),
          group = srbct$class, # colour by class label
          ind.names = FALSE, 
          ellipse = TRUE, # include 95% confidence ellipse
          legend = TRUE, 
          title = ' (a) sPLS-DA on SRBCT, comp 1 & 2')

plotIndiv(MyResult.splsda.final,
          comp = c(1,3),
          group = srbct$class, # colour by class label
          ind.names = FALSE, 
          ellipse = TRUE, # include 95% confidence ellipse
          legend = TRUE, 
          title = ' (b) sPLS-DA on SRBCT, comp 1 & 3')

# A Cluster Image Map: depicts gene expression levels and cluster similar genes
# together using Eucidean distance with a complete agglomeration method
# set the styling of the legend to be homogeneous with previous plots
legend=list(legend = levels(Y), # set of classes
            col = unique(color.mixo(Y)), # set of colours
            title = "Tumour Type", # legend title
            cex = 0.7) # legend size

# generate the CIM, using the legend and colouring rows by each sample's class
# Each group has a signature set of genes that are highly expressed
cim <- cim(MyResult.splsda.final, 
           row.sideColors = color.mixo(Y), 
           legend = legend)

# Stability plot
# The higher the bar for each variable the more frequently it is selected during
# cross validation as an important variable
# Less overall area means there exists multiple combinations that distinguish
# the clusters on that component
# plot the stability of each feature for the first three components, 'h' type refers to histogram
# Apparently this plot changes with the settings for cross validation... Hard to
# say how useful it is
par(mfrow=c(1,3))
plot(perf.splsda$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
plot(perf.splsda$features$stable[[3]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features',
     main = '(c) Comp 3', las =2)

# Correlation circle
# Overlay with sample plot to interpret: further to the edge and closer to a component
# means higher correlation between component and original variable
# form simplified gene names
var.name.short <- substr(srbct$gene.name[, 2], 1, 10) 

# generate correlation circle plot
plotVar(MyResult.splsda.final, 
        comp = c(1,2), 
        var.names = list(var.name.short), 
        cex = 3) 

