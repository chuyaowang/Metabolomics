cat("\014")
rm(list = ls())

library(readr)
library(dplyr)

data("iris")

# Get names and group name for convenience
samplename <- seq_along(iris$Species)
samplegroup <- iris$Species

# Get only the numerical part
sampleprocessed <- iris[,1:ncol(iris)-1]

# Optional horizontal normalization
sampleprocessed <- apply(sampleprocessed,1,function(x) {
    . <- norm(x,type="2")
    return(x/.)
  }) %>%
  t # apply merges the results by columns, so need transpose

## prcomp way -----

pca <- prcomp(as.matrix(sampleprocessed),
              retx = TRUE,
              center = TRUE,
              scale. = TRUE)

# Eigenvectors, all have unit length; linear combination coefficients to get
# the PC scores
eigvectors <- pca$rotation %>%
  as.data.frame

# Variances of each PC; also the eigenvalues for the eigenvectors
eigvals <- pca$sdev^2

# x is PC scores
# x is unstandardized
scores <- pca$x %>%
  as.data.frame
scores_std <- scale(scores) %>%
  as.data.frame

# Loadings are eigenvectors scaled by sqrt of eigenvalues;
# such that their magnitude equals to eigenvalues
# each loading equals to Pearson correlation between original data and PC scores
. <- colnames(eigvectors)
loadings <- as.matrix(eigvectors) %*% diag(sqrt(eigvals)) %>%
  as.data.frame
colnames(loadings) <- .

# Varimax rotation: rotate the principal axes to bring loadings closer to 0 or 1
# for better interpretation
# Usually done after deciding which PC to keep
. <- varimax(loadings %>%
               select(PC1:PC2) %>%
               as.matrix(),
             normalize = FALSE)
rLoadings <- .$loadings
rotmat <- .$rotmat

# Use rotmat to rotate the standardized PC scores
# Each rotated loading equals to Pearson correlation between original data and
# rotated PC scores (somehow only true when rotating standardized scores)
rScores <- as.matrix(select(scores_std,PC1:PC2)) %*% rotmat %>%
  as.data.frame
colnames(rScores) <- c("RC1","RC2")

## manual way -----
# samplescale <- scale(sampleprocessed,
#                      center = TRUE,
#                      scale = TRUE)
# eigs_manual <- eigen(cov(samplescale))
# eigs_manual_vectors <- eigs_manual$vectors
# eigs_manual_values <- eigs_manual$values

## Plot -----
# Scree plot: used to select PCs to keep
library(ggplot2)
library(hrbrthemes)

# scree <- as.data.frame(list(colnames(eigvectors),eigvals/sum(eigvals)))
# colnames(scree) <- c("PC","Percent") 
# ggplot(scree, aes(x = PC,y = Percent,group=1)) +
#   geom_line(color="grey") +
#   geom_point() +
#   ggtitle("Scree Plot") +
#   xlab("Principal Components") +
#   ylab("Percent Variance Explained") +
#   theme_modern_rc()

# PC1 PC2 biplot
library(ggbiplot)
ggbiplot(pca, 
         obs.scale = 1, 
         var.scale = 1, 
         labels = samplename,
         groups = samplegroup, 
         ellipse = TRUE, 
         circle = TRUE, 
         var.axes = T,
         varname.size = 1.5) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top') +
  theme_minimal()

# In correlation circle, variables with longer arrows are better represented for
# that PC
# The angle between the arrows represent correlation between the variables:
# small angle means high + correlation, 90 degrees means no correlation, 
# opposite means high - correlation
# Correlation circles are plotted using loadings, which represent correlation
# between original variables and the PCs. Hence, they have a radius of 1.

# The ellipses are confidence ellipses, 2D version of confidence intervals.
# Means that 95% percent sure of the circle contains the true population
# distribution (population from which the sample is drawn)