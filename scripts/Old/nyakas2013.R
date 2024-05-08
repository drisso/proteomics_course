# library("devtools")
# install_github("sgibb/MALDIquantExamples")

## the main MALDIquant package
library("MALDIquant")
## the import/export routines for MALDIquant
library("MALDIquantForeign")
## example data
library("MALDIquantExamples")

## import the spectra
spectra_imp <- import(getPathNyakas2013(), verbose=FALSE)
# spectra <- smoothIntensity(spectra_imp, method="SavitzkyGolay",
#                            halfWindowSize=5)
# baseline <- estimateBaseline(spectra[[100]], method="SNIP", iterations=50)
# plot(spectra_imp[[100]])
# lines(baseline, col="red", lwd=2)
# spectra <- removeBaseline(spectra, method="SNIP", iterations=100)
# plot(spectra[[50]])

spectra <- removeBaseline(spectra_imp, method="median")
spectra <- smoothIntensity (spectra, method = "MovingAverage",halfWindowSize = 5)
spectra <- calibrateIntensity(spectra, method="TIC")

meanSpectrum <- averageMassSpectra(spectra)
roi <- detectPeaks(meanSpectrum, SNR=3,
                   halfWindowSize=5)
plot(meanSpectrum, main="Mean Spectrum")
points(roi, col="red")
## find order of peak intensities
o <- order(intensity(roi), decreasing=TRUE)
## plot MSI slice for the highest one
plotMsiSlice(spectra, center=mass(roi)[o[1]], tolerance=0.5)
plotMsiSlice(spectra, center=mass(roi)[o[2:3]], tolerance=0.5)

plotMsiSlice(spectra, center=mass(roi)[o[1:2]], tolerance=0.5,
             combine=TRUE,
             colRamp=list(colorRamp(c("#000000", "#FF00FF")),
                          colorRamp(c("#000000", "#00FF00"))))

slices <- msiSlices(spectra, center=mass(roi), tolerance=0.5)
attributes(slices)
head(coordinates(spectra))
head(coordinates(spectra, adjust=TRUE))

peaks <- detectPeaks(spectra, method="MAD",halfWindowSize=5,SNR=3)
peaks <- binPeaks(peaks)

#FINAL MATRIX
intMatrix <- intensityMatrix(peaks, spectra)
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

## 1. Dimensionality reduction

## Simple PCA with prcomp
system.time(pca <- prcomp(intMatrix, scale=TRUE))
plot(pca, type='l')
pca_summary <- summary(pca)
pca_summary$importance[,1:10]

## The Bioconductor BiocSingular package implements several faster algorithms
library(BiocSingular)
system.time(rpca <- runPCA(intMatrix, rank=50, scale=TRUE, BSPARAM = RandomParam()))
system.time(ipca <- runPCA(intMatrix, rank=50, scale=TRUE, BSPARAM = IrlbaParam()))

## Visualize PCA
library(ggplot2)
theme_set(theme_classic())

df <- as.data.frame(cbind(ipca$x, coord))
df$kmeans2 <- as.factor(km$cluster)

ggplot(df, aes(x, y, color=kmeans2)) +
    geom_point()

ggplot(df, aes(PC1, PC2, color=kmeans2)) +
    geom_point()

ggplot(df, aes(x, y, color=PC1)) +
    geom_point() +
    scale_color_viridis_c(option = "magma")

ggplot(df, aes(x, y, color=PC2)) +
    geom_point() +
    scale_color_viridis_c(option = "magma")

ggplot(df, aes(PC1, PC2, color=y)) +
    geom_point() +
    scale_color_viridis_c(option = "magma")

ggplot(df, aes(PC1, PC2)) +
    geom_point()

idx <- names(sort(abs(pca$rotation[,1]), decreasing = TRUE))[1]

df$prot1 <- intMatrix[,idx]

p1 <- ggplot(df, aes(PC1, PC2, color=prot1)) +
    geom_point() +
    scale_color_viridis_c(option = "magma")

p2 <- ggplot(df, aes(x, y, color=prot1)) +
    geom_point() +
    scale_color_viridis_c(option = "magma")

idx <- names(sort(abs(pca$rotation[,2]), decreasing = TRUE))[1]

df$prot2 <- intMatrix[,idx]

p3 <- ggplot(df, aes(PC1, PC2, color=prot2)) +
    geom_point() +
    scale_color_viridis_c(option = "magma")

p4 <- ggplot(df, aes(x, y, color=prot2)) +
    geom_point() +
    scale_color_viridis_c(option = "magma")

library(patchwork)
p1 + p2 + p3 + p4

## T-sne
library(Rtsne)
## system.time(tsne1 <- Rtsne(intMatrix))
system.time(tsne <- Rtsne(ipca$x[,1:50], pca=FALSE))

df$tsne1 <- tsne$Y[,1]
df$tsne2 <- tsne$Y[,2]

t1 <- ggplot(df, aes(tsne1, tsne2, color=prot1)) +
    geom_point() +
    scale_color_viridis_c(option = "magma") +
    theme(legend.position = "none")

t2 <- ggplot(df, aes(tsne1, tsne2, color=prot2)) +
    geom_point() +
    scale_color_viridis_c(option = "magma") +
    theme(legend.position = "none")

p1 + t1 + p3 + t2

## UMAP
library(uwot)
system.time(um <- umap(ipca$x[,1:50]))

df$umap1 <- um[,1]
df$umap2 <- um[,2]

u1 <- ggplot(df, aes(umap1, umap2, color=prot1)) +
    geom_point() +
    scale_color_viridis_c(option = "magma") +
    theme(legend.position = "none")

u2 <- ggplot(df, aes(umap1, umap2, color=prot2)) +
    geom_point() +
    scale_color_viridis_c(option = "magma") +
    theme(legend.position = "none")

(p1 + t1 + u1) / (p3 + t2 + u2)


## 2. Clustering

## because of curse of dimensionality, better to cluster in reduced space

## simple k-means
km_rid <- kmeans(ipca$x, centers=2)

## compare with k-means on full data
table(km$cluster, km_rid$cluster)
df$kmeans2rid <- as.factor(km_rid$cluster)

ggplot(df, aes(PC1, PC2, color=kmeans2rid)) +
    geom_point()

ggplot(df, aes(x, y, color=kmeans2rid)) +
    geom_point()

ggplot(df, aes(tsne1, tsne2, color=kmeans2rid)) +
    geom_point()

## evaluate clustering: Silhouette
library(cluster)
d <- dist(ipca$x)
sil <- silhouette(km_rid$cluster, d)
plot(sil)

## use silhouette to choose number of clusters
ks <- 2:6
km_res <- lapply(ks, function(k) kmeans(ipca$x, centers=k))
sapply(km_res, function(x) mean(silhouette(x$cluster, d)[,3]))

df$kmeans4rid <- as.factor(km_res[[3]]$cluster)

ggplot(df, aes(PC1, PC2, color=kmeans4rid)) +
    geom_point()

ggplot(df, aes(x, y, color=kmeans4rid)) +
    geom_point()

ggplot(df, aes(tsne1, tsne2, color=kmeans4rid)) +
    geom_point()

## hierarchical clustering
hc <- hclust(d)
plot(as.dendrogram(hc))
df$hclust6<- as.factor(cutree(hc, k=6))

ggplot(df, aes(PC1, PC2, color=hclust6)) +
    geom_point()

ggplot(df, aes(x, y, color=hclust6)) +
    geom_point()

ggplot(df, aes(tsne1, tsne2, color=hclust6)) +
    geom_point()

## network-based clustering (walktrap method)

## first compute the shared nearest neighbor graph
library(bluster)
graph <- makeSNNGraph(ipca$x)

## then we can use several algorithms (here walktrap and louvain)
cl <- igraph::cluster_walktrap(graph)
cl2 <- igraph::cluster_louvain(graph)

df$cl_walktrap <- as.factor(cl$membership)
df$cl_louvain <- as.factor(cl2$membership)

ggplot(df, aes(PC1, PC2, color=cl_louvain)) +
    geom_point()

ggplot(df, aes(x, y, color=cl_louvain)) +
    geom_point(size=3)

ggplot(df, aes(tsne1, tsne2, color=cl_louvain)) +
    geom_point()

## compare results
library(mclust)
adjustedRandIndex(df$cl_louvain, df$cl_walktrap)
adjustedRandIndex(df$cl_louvain, df$kmeans4rid)
adjustedRandIndex(df$cl_louvain, df$hclust6)



