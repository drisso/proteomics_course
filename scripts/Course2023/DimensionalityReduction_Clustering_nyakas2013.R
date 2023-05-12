# library("devtools")
# install_github("sgibb/MALDIquantExamples")

# the main MALDIquant package
library("MALDIquant")
# the import/export routines for MALDIquant
library("MALDIquantForeign")
# example data
library("MALDIquantExamples")

# import the spectra
spectra_imp <- import(getPathNyakas2013(), verbose=FALSE)
# how many spectra we have?
length(spectra_imp)

# Now you can try to pre-process spectra!
# Attention: this data set comes from the analysis of intact protein, m/z range from 3000 to 20000
# plot a spectra for example to see the range of masses, before performed pre-processing on your data
# REMEMBER: Accessing the "x" element using double square bracket "[[]]"
plot(spectra_imp[[100]])

# which are the step to be performed?
# and the methods to be used for spectra of intact proteins?
spectra <- removeBaseline(spectra_imp, method="median")
spectra <- smoothIntensity (spectra, method = "MovingAverage",halfWindowSize = 5)
spectra <- calibrateIntensity(spectra, method="TIC")

#calculate mean spectrum of pre-processed spectra for some exploratory images 
meanSpectrum <- averageMassSpectra(spectra)
peaks_mean <- detectPeaks(meanSpectrum, SNR=3,
                          halfWindowSize=5)
plot(meanSpectrum, main="Mean Spectrum")
points(peaks_mean, col="red")

# find order of peak intensities
o <- order(intensity(peaks_mean), decreasing=TRUE)

# plot MSI slice for the highest one
plotMsiSlice(spectra, center=mass(peaks_mean)[o[1]], tolerance=0.5)

# plot MSI slice for the second and third one
plotMsiSlice(spectra, center=mass(peaks_mean)[o[2:3]], tolerance=0.5)

# how we can solve this Warning? Seems we are not able to see two peaks simultaneously
# REMEMBER: look at the help function, look at the Arguments "combine"
# Use different colors
plotMsiSlice(spectra, center=mass(peaks_mean)[o[1:2]], tolerance=0.5,
             combine=TRUE,
             colRamp=list(colorRamp(c("#000000", "#FF00FF")),
                          colorRamp(c("#000000", "#00FF00"))))


# the help function suggest to look "msiSlices" for details
slices <- msiSlices(spectra, center=mass(peaks_mean), tolerance=0.5)
# extract information about dimension, coordinates...
attributes(slices)
head(coordinates(spectra))
head(coordinates(spectra, adjust=TRUE))

# extract peaks from spectra
peaks <- detectPeaks(spectra, method="MAD",halfWindowSize=5,SNR=3)
peaks <- binPeaks(peaks)

# FINAL MATRIX
intMatrix <- intensityMatrix(peaks, spectra)
# dimension of the generated matrix
dim(intMatrix)

# Now we can control if we pre-processed spectra correctly.
# follow these lines of codes and see if you obtain the same image of Nyakas2013 in the pdf
km <- kmeans(intMatrix, centers=2)

coord <- coordinates(spectra, adjust=TRUE)
maxPixels <- apply(coord, MARGIN=2, FUN=max)
m <- matrix(NA, nrow=maxPixels["x"], ncol=maxPixels["y"])
m[coord] <- km$cluster

rgbCluster <- function(x) {
  col <- matrix(c(255, 0, 0,
                  0, 255, 0), nrow=2, byrow=TRUE)
  col[x, ]
}
plotMsiSlice(m, colRamp=rgbCluster, scale=FALSE)

# We performed it correctly!!

# Save the Final Matrix, we will need it in the next lessons
saveRDS(intMatrix, "intMatrix.RDS")



# 1. Dimensionality reduction
# Two main goals of dimensionality reduciton:
# - find structure in features 
# - aid in visualization

# Read the Final Matrix pre-processed
intMatrix <- readRDS("intMatrix.RDS")

# Simple PCA with prcomp function in R
# It can take a lot of time, so use "system.time" to return CPU (and other) times that the function "prcomp" used
system.time(pca <- prcomp(intMatrix, scale=TRUE))
# The PCA model produce additional diagnostic and output components:
# center - the column means used to center the data
# scale - the column sd used to scale the data
# rotation - the direction of the principal component vectors in terms of the original features/variables. This information somehow allows you to define new data in terms of the original principal components.
# x - the value of each observation in the original dataset projected to the principal components
head(pca$center)
head(pca$scale)
head(pca$rotation[,1:10])
head(pca$x[,1:10])
# Plot the variance explained by each Principal Component 
plot(pca, type='l')
# Extract the proportion of variance explained by the first 10 Principal Components
pca_summary <- summary(pca)
# Remember: the standard deviation in the square of the Variance
pca_summary$importance[,1:10]


# Visualize PCA
# Use "cbind" to combine the data frame of the PC components resulting from prcomp function and 
# the coordinates of the mouse kidney slice to visualized graphically the PCA results on the MsiSlice
?cbind
# cbind combine the two dataset by row, they must all have the same number or rows
df <- as.data.frame(cbind(pca$x, coord))
dim(pca$x)
dim(coord)
df$kmeans2 <- as.factor(km$cluster)
head(df)
# Look at the last three columns to see what we added to the pca data.frame
head(df[,2223:2225])

library(ggplot2)
# Plot the mouse kidney slice and colour by clustering resulting from kmeans
ggplot(df, aes(x, y, color=kmeans2)) +
  geom_point()
# Delete the gray background, improve image quality 
theme_set(theme_classic())
ggplot(df, aes(x, y, color=kmeans2)) +
  geom_point()

# Plot PCA components
ggplot(df, aes(PC1, PC2)) +
  geom_point()
# Color by k-mean clusters
ggplot(df, aes(PC1, PC2, color=kmeans2)) +
  geom_point()

# If we want to see the two graph simultaneously?
library(patchwork)
p<-ggplot(df, aes(PC1, PC2)) +
  geom_point()
p_pca<-ggplot(df, aes(PC1, PC2, color=kmeans2)) +
  geom_point()
p+p_pca

# Investigate the relationship between point and principal components
p1<-ggplot(df, aes(x, y, color=PC1)) +
  geom_point() +
  scale_color_viridis_c(option = "magma")

p2<-ggplot(df, aes(x, y, color=PC2)) +
  geom_point() +
  scale_color_viridis_c(option = "magma")

p3<-ggplot(df, aes(PC1, PC2, color=y)) +
  geom_point() +
  scale_color_viridis_c(option = "magma")

p4<-ggplot(df, aes(PC1, PC2, color=x)) +
  geom_point() +
  scale_color_viridis_c(option = "magma")

(p1+p2)/(p3+p4)

# We can also look at the distribution of the most impact variables over the sample
# Extract the most impact protein for the first principal component
idx <- names(sort(abs(pca$rotation[,1]), decreasing = TRUE))[1]
idx
# Extract intensity value of this protein over the slice sample
df$prot1 <- intMatrix[,idx]
# Plot the distribution of the protein along the PCA plot and the kidney mouse sample
p1 <- ggplot(df, aes(PC1, PC2, color=prot1)) +
  geom_point() +
  scale_color_viridis_c(option = "magma")
p2 <- ggplot(df, aes(x, y, color=prot1)) +
  geom_point() +
  scale_color_viridis_c(option = "magma")

# Extract the most impactful protein for the second principal component
idx <- names(sort(abs(pca$rotation[,2]), decreasing = TRUE))[1]
idx
# Extract intensity value of this protein over the slice sample
df$prot2 <- intMatrix[,idx]
# Plot the distribution of the protein along the PCA plot and the kidney mouse sample
p3 <- ggplot(df, aes(PC1, PC2, color=prot2)) +
  geom_point() +
  scale_color_viridis_c(option = "magma")

p4 <- ggplot(df, aes(x, y, color=prot2)) +
  geom_point() +
  scale_color_viridis_c(option = "magma")

(p1 + p2)/ (p3 + p4)


# Compare PCA with tsne and umap
## T-sne
library(Rtsne)
system.time(tsne1 <- Rtsne(intMatrix))
# extract the first two components
df$tsne1 <- tsne1$Y[,1]
df$tsne2 <- tsne1$Y[,2]
# plot
p_tsne<-ggplot(df, aes(tsne1, tsne2, color=kmeans2)) +
  geom_point()

# A lot of time...
# see the help function ?Rtsne to take a look of the Arguments of the function, 
# e.g. pca, pca_scale, perplexity
# We have already compute PCA, so use the first 50 PC as matrix and set pca=FALSE 
system.time(tsne <- Rtsne(pca$x[,1:50], pca=FALSE))
# extract the first two components
df$tsne1_pca <- tsne$Y[,1]
df$tsne2_pca <- tsne$Y[,2]
# plot
p_tsne_pca<-ggplot(df, aes(tsne1_pca, tsne2_pca, color=kmeans2)) +
  geom_point()

p_tsne+p_tsne_pca
# Why they are different?
# Look at the help function!
# pca:	logical; Whether an initial PCA step should be performed (default: TRUE)
# pca_scale:	logical; Should data be scaled before pca is applied? (default: FALSE)
# CONCLUSION: with the command Rtsne(intMatrix), the method performed a PCA without scaling
# So, How much change PCA without scaling?
pca_notScale <- prcomp(intMatrix, scale=FALSE)
df_pcaNotScale <- as.data.frame(cbind(pca_notScale$x, coord))
df_pcaNotScale$kmeans2 <- as.factor(km$cluster)
p_pcaNotScale<-ggplot(df_pcaNotScale, aes(PC1, PC2, color=kmeans2)) +
  geom_point()+xlab("PC1_notScale")+ylab("PC2_notScale")
p_pca+p_pcaNotScale

# PCA change consistently, tsne is less affected by outliers and high variance
(p_pca+p_pcaNotScale)/(p_tsne_pca+p_tsne)


## UMAP (similar to tsne)
library(umap)
## system.time(Umap1 <- umap(intMatrix))
system.time(Umap <- umap(pca$x[,1:50], pca=FALSE))
# extract the first two components
df$Umap1 <- Umap$layout[,1]
df$Umap2 <- Umap$layout[,2]
# plot
p_umap_pca<-ggplot(df, aes(Umap1, Umap2, color=kmeans2)) +
  geom_point()

# compare the results of the three methods
p_pca+p_tsne_pca+p_umap_pca


#RETURN to PCA!
# We learn how to plot results, but which original variables contribute most?
library(ggplot2)
library("factoextra")
# "factoextra" a package for quick Principal Component Analysis data visualization 

# The scree plot is used to visualize the importance of each principal component and 
# can be used to determine the number of principal components to retain. 
# The scree plot can be generated using the fviz_eig() function. 
fviz_eig(pca, addlabels = TRUE)

# Graph of the variables
# With the biplot, it is possible to visualize the similarities and dissimilarities between the samples,
# and further shows the impact of each attribute on each of the principal components.
fviz_pca_var(pca, col.var = "black")
# ok... too much... go to the help function! 
# There is a way to plot only the most informative proteins?
?fviz_pca_var
fviz_pca_var(pca, col.var = "black", select.var = list(contrib = 50))

# Graph of individuals
fviz_pca_ind(pca, geom="point")
# Change the size of points
fviz_pca_ind(pca, geom="point", pointsize = 4)
fviz_pca_ind(pca, geom="point", pointsize = 4, habillage=df$kmeans2)

# Contribution of each variable 
# The goal of this third visualization is to determine how much each variable is represented in a given component. 
# var$coord: coordinates of variables to create a scatter plot
# var$cos2: represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
# var$contrib: contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
var <- get_pca_var(pca)
var
# Such a quality of representation is called the Cos2 and corresponds to the square cosine, and it is computed using the fviz_cos2 function.
# A low value means that the variable is not perfectly represented by that component. 
# A high value, on the other hand, means a good representation of the variable on that component.
library("corrplot")
corrplot(var$cos2[1:10,1:10], is.corr=FALSE)
corrplot(var$cos2[1:50,1:50], is.corr=FALSE)
# PROBLEM: is impossible to plot all the variables, we want the top ones!
# The function fviz_contrib() creates a barplot of row/column contributions. 
fviz_contrib(pca, choice = "var", axes = 1:2)
# Select the top 50
fviz_contrib(pca, choice = "var", axes = 1:2, top=50)
?fviz_contrib
# A reference dashed line is also shown on the barplot. 
# This reference line corresponds to the expected value if the contribution where uniform.
# For a given dimension, any row/column with a contribution above the reference line 
# could be considered as important in contributing to the dimension.
1/dim(pca$x)[2]
1/50
# Ok.... we see the 50 impact protein in the PCA, but which are the most important in each component?
# Do you remember the distribution plot of the two main important protein over the kidney mouse slice?
# In another way we extract the same information
# Contributions of variables to PC1
graph1<-fviz_contrib(pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
graph2<-fviz_contrib(pca, choice = "var", axes = 2, top = 10)
graph1+graph2

# Rescale graph to have the same range on the y-axis
# Contributions of variables to PC1
graph1<-fviz_contrib(pca, choice = "var", axes = 1, top = 10, ylim=c(0,0.3))
# Contributions of variables to PC2
graph2<-fviz_contrib(pca, choice = "var", axes = 2, top = 10, ylim=c(0,0.3))
graph1+graph2


# Last... we left for a moment the correlation plot
# How we can extract more information with that graphical representation?
idx1 <- names(sort(abs(pca$rotation[,1]), decreasing = TRUE))[1:10]
idx1
idx2 <- names(sort(abs(pca$rotation[,2]), decreasing = TRUE))[1:10]
idx2
idx3 <- names(sort(abs(pca$rotation[,3]), decreasing = TRUE))[1:10]
idx3
idx<-c(idx1,idx2,idx3)
ind<-NULL
for(int in 1:length(idx)){
  ind[int]<-which(as.numeric(colnames(intMatrix))==as.numeric(idx)[int])
}
corrplot(var$cos2[ind,1:10], is.corr=FALSE)


## More information can be search here!
#https://rpkgs.datanovia.com/factoextra/


## The Bioconductor BiocSingular package implements several faster algorithms
# useful in high dimensional data set, we see how much time expensive is pca
#install.packages("BiocManager")
#BiocManager::install("BiocSingular")
#library(BiocSingular)
#system.time(rpca <- runPCA(intMatrix, rank=50, scale=TRUE, BSPARAM = RandomParam()))
#system.time(ipca <- runPCA(intMatrix, rank=50, scale=TRUE, BSPARAM = IrlbaParam()))




# 2. Clustering
# unsupervised learning: goal is to find structure in unlabeld data (i.e. find homogenous subgroups within a population)

# Because of curse of dimensionality, better to cluster in reduced space
dim(pca$x)
dim(intMatrix)
# simple k-means
km_rid <- kmeans(pca$x, centers=2)
km_rid
# Too much information, impossible to extract information
# Take a look to the first 10 components
# Cluster means:
km_rid$centers[,1:10]
df$kmeans2rid <- as.factor(km_rid$cluster)
# see clustering results of kmeans in reduced space
p_pca<-ggplot(df, aes(PC1, PC2, color=kmeans2rid)) +
  geom_point()
p_tsne<-ggplot(df, aes(tsne1_pca, tsne2_pca, color=kmeans2rid)) +
  geom_point()
p_pca+p_tsne

# Compare with k-means on full data
table(km$cluster, km_rid$cluster)

# Compare graphically k-means in reduced space vs. k-mean on original data matrix
p_Kpca<-ggplot(df, aes(x, y, color=kmeans2rid)) +
  geom_point()
p_Kor<-ggplot(df, aes(x, y, color=kmeans2)) +
  geom_point()
p_Kpca+p_Kor
# Add some information to the plot to make it more comprehensibly
p_Kpca<-ggplot(df, aes(x, y, color=kmeans2rid)) +
  geom_point(size=3)+xlab("X Coordinate")+ylab("Y Coordinate")+ggtitle("K-means on reduce space")+
  theme(plot.title = element_text(color="Black", size=20, face="bold"),
        axis.title.x = element_text(color="Black", size=15, face="bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=12),
        axis.title.y = element_text(color="Black", size=15, face="bold"),
        text = element_text(size = 15), panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray"),
        axis.line = element_line(colour = "black",size = 1), legend.position = "bottom",legend.box = "horizontal")+
  guides(color = guide_legend(nrow=1, byrow = TRUE, title="K-means"))
p_Kor<-ggplot(df, aes(x, y, color=kmeans2)) +
  geom_point(size=3)+xlab("X Coordinate")+ylab("Y Coordinate")+ggtitle("K-means on original matrix")+
  theme(plot.title = element_text(color="Black", size=20, face="bold"),
        axis.title.x = element_text(color="Black", size=15, face="bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=12),
        axis.title.y = element_text(color="Black", size=15, face="bold"),
        text = element_text(size = 15), panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray"),
        axis.line = element_line(colour = "black",size = 1), legend.position = "bottom",legend.box = "horizontal")+
  guides(color = guide_legend(nrow=1, byrow = TRUE, title="K-means"))
p_Kpca+p_Kor


# Handling random algorithms
# Compare different k-means runned on the same dataset
# How we can compare them?
# withinss:	Vector of within-cluster sum of squares, one component per cluster.
# tot.withinss:	Total within-cluster sum of squares, i.e. sum(withinss)
# Set up 2 x 3 plotting grid
par(mfrow = c(2, 3))
# Set seed: seed must be outside the cicle "for"
set.seed(1)
for(i in 1:6) {
  # Run kmeans() on x with three clusters and one start
  km.out <- kmeans(pca$x, centers = 2, nstart = 1)
  # Plot clusters
  plot(coord, col = km.out$cluster, 
       main = paste(km.out$withinss[1]," - ",km.out$withinss[2]), 
       xlab = "", ylab = "")
}
# Change seed at each iteration: seed must be inside the cicle "for"
# What can happened? 
# We saw Main title too long, retain only the first two decimal places 
for(i in 1:6) {
  set.seed(i)
  # Run kmeans() on x with three clusters and one start
  km.out <- kmeans(pca$x, centers = 2, nstart = 1)
  # Plot clusters
  plot(coord, col = km.out$cluster, 
       main = paste(round(km.out$withinss[1],2)," - ",round(km.out$withinss[2],2)), 
       xlab = "", ylab = "")
}


# How to select the right number of clusters?
# Initialize total within sum of squares error: wss
wss <- rep(0,15)
# For 1 to 15 cluster centers
for (i in 1:15) {
  km.out <- kmeans(pca$x, centers = i, nstart = 1)
  # Save total within sum of squares to wss variable
  wss[i] <- km.out$tot.withinss
}
# Plot total within sum of squares vs. number of clusters
plot(1:15, wss, type = "b", 
     xlab = "Number of Clusters", 
     ylab = "Within groups sum of squares")
# What is happening??
# Remember to set up the plotting grid to one, after you finished to use it
par(mfrow = c(1, 1))
plot(1:15, wss, type = "b", 
     xlab = "Number of Clusters", 
     ylab = "Within groups sum of squares")

# Set centers equal to the number of clusters corresponding to the elbow location
km_rid5 <- kmeans(pca$x, centers=5)
df$kmeans5rid <- as.factor(km_rid5$cluster)
ggplot(df, aes(x, y, color=kmeans5rid)) +
  geom_point()
# Better quality plot
ggplot(df, aes(x, y, color=kmeans5rid)) +
  geom_point(size=3)+xlab("X Coordinate")+ylab("Y Coordinate")+ggtitle("K-means on reduce space")+
  theme(plot.title = element_text(color="Black", size=20, face="bold"),
        axis.title.x = element_text(color="Black", size=15, face="bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=12),
        axis.title.y = element_text(color="Black", size=15, face="bold"),
        text = element_text(size = 15), panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray"),
        axis.line = element_line(colour = "black",size = 1), legend.position = "bottom",legend.box = "horizontal")+
  guides(color = guide_legend(nrow=1, byrow = TRUE, title="K-means"))


# Hierarchical clustering
# Two approaches. bottom up and top down. 
# We will focus on the function "hclust", which is bottom up process
# Assign each point to its own cluster
# Joint the two closest clusters/points into a new cluster
# Keep going until there is one cluster
# The way you calculate the distance between clusters is a parameter of the function
# First calculate the euclidean distance between all points using the dist() function
d <- dist(pca$x)
# Hierarchical Clustering function
hc <- hclust(d)
plot(as.dendrogram(hc))
# The cutree() function in R lets you split the hierarchical clusters 
# into set clusters by number (k) or by distance (height)
abline(h = 190, col = "red")
# Cut by height
cutree(hc, h = 190)
# How many spectra fall in each cluster? Very unbalanced! 
table(cutree(hc, h = 190))
# Cut by number of cluster
df$hclust5<- as.factor(cutree(hc, k=5))

p_pca_h<-ggplot(df, aes(PC1, PC2, color=hclust5)) +
  geom_point()
p_tsne_h<-ggplot(df, aes(tsne1_pca, tsne2_pca, color=hclust5)) +
  geom_point()
p_hclus<-ggplot(df, aes(x, y, color=hclust5)) +
  geom_point()

(p_pca+p_tsne+p_Kpca)/(p_pca_h+p_tsne_h+p_hclus)

# Clustering linkage
# 4 methods to measure distance between clusters
# complete: pairwise similarty between all observations in cluster 1 and 2, uses largest of similarities
# single: same as above but uses the smallest of similarities
# average: same as above but uses average of similarities
# centroid: finds centroid of cluster 1 and 2, uses similarity between two centroids

# Cluster using complete linkage: hclust.complete
hclust.complete <- hclust(d, method = "complete")
# Cluster using average linkage: hclust.average
hclust.average <- hclust(d, method = "average")
# Cluster using single linkage: hclust.single
hclust.single <- hclust(d, method = "single")

# Plot dendrogram of hclust.complete
plot(hclust.complete, main = "Complete")
# Plot dendrogram of hclust.average
plot(hclust.average, main = "Average")
# Plot dendrogram of hclust.single
plot(hclust.single, main = "Single")




# ## network-based clustering (walktrap method)
# 
# ## first compute the shared nearest neighbor graph
# library(bluster)
# graph <- makeSNNGraph(ipca$x)
# 
# ## then we can use several algorithms (here walktrap and louvain)
# cl <- igraph::cluster_walktrap(graph)
# cl2 <- igraph::cluster_louvain(graph)
# 
# df$cl_walktrap <- as.factor(cl$membership)
# df$cl_louvain <- as.factor(cl2$membership)
# 
# ggplot(df, aes(PC1, PC2, color=cl_louvain)) +
#   geom_point()
# 
# ggplot(df, aes(x, y, color=cl_louvain)) +
#   geom_point(size=3)
# 
# ggplot(df, aes(tsne1, tsne2, color=cl_louvain)) +
#   geom_point()
# 
# ## compare results
# library(mclust)
# adjustedRandIndex(df$cl_louvain, df$cl_walktrap)
# adjustedRandIndex(df$cl_louvain, df$kmeans5rid)
# adjustedRandIndex(df$cl_louvain, df$hclust5)
# table(df$hclust5,df$kmeans5rid)
