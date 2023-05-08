library(palmerpenguins)
library(ggplot2)
theme_set(theme_bw())

data("penguins")
ggplot(penguins, aes(x = flipper_length_mm, y = body_mass_g, color = species)) +
geom_point()

xy <- na.omit(penguins[,c("flipper_length_mm", "body_mass_g")])
xynorm <- scale(xy, scale = TRUE, center = TRUE)
e <- eigen(cov(xynorm))

plot(xynorm)
lines(xynorm[,1], e$vectors[2,1]/e$vectors[1,1] * xynorm[,1], col=2, lwd=3)
lines(xynorm[,1], e$vectors[2,2]/e$vectors[1,2] * xynorm[,1], col=4, lwd=3)

ggplot(as.data.frame(xynorm), aes(x = flipper_length_mm, y = body_mass_g)) +
    geom_point() +
    geom_abline(intercept = 0, slope = e$vectors[2,1]/e$vectors[1,1], color=palette()[2], size=1.5) +
    ggtitle("First PC")

ggplot(as.data.frame(xynorm), aes(x = flipper_length_mm, y = body_mass_g)) +
    geom_point() +
    geom_abline(intercept = 0, slope = e$vectors[2,1]/e$vectors[1,1], color=palette()[2], size=1.5) +
    geom_abline(intercept = 0, slope = e$vectors[2,2]/e$vectors[1,2], color=palette()[4], size=1.5) +
    ggtitle("First two PCs")

penguins <- na.omit(penguins)
pca <- prcomp(penguins[,c("flipper_length_mm", "body_mass_g")], scale. = TRUE)
df <- cbind(as.data.frame(pca$x), species=penguins$species)

ggplot(df, aes(x = PC1, y = PC2, color = species)) +
    geom_point() +
    geom_vline(xintercept=0, linetype = "dashed") +
    geom_hline(yintercept=0, linetype = "dashed")

ggplot(penguins, aes(x = flipper_length_mm, y = bill_length_mm, color = species)) +
    geom_point() +
    geom_point(data = penguins[1:2,], aes(x = flipper_length_mm, y = bill_length_mm), colour="black", size=4, pch=21, fill=NA, stroke=1)

penguins <- na.omit(penguins)
data <- penguins[, c("flipper_length_mm", "bill_length_mm")]
d <- dist(data)
hc <- hclust(d)
dend <- as.dendrogram(hc)

library(dendextend)
labels_colors(dend) <- palette()[-1][as.factor(penguins$species)]
plot(dend)

xy <- data.frame(x=c(5, -1, 1, -3), y=c(3, 1, -2, -2))
rownames(xy) <- c("A", "B", "C", "D")
xy$label <- factor(c("1", "1", "2", "2"))

ggplot(xy, aes(x, y, color=label)) +
    geom_point(size=4) +
    ggtitle("Iteration 1")

xy$label <- factor(c("1", "2", "2", "2"))

ggplot(xy, aes(x, y, color=label)) +
    geom_point(size=4) +
    ggtitle("Iteration 2")

k3 <- kmeans(data, centers=3)
penguins$k3 <- as.factor(k3$cluster)

ggplot(penguins, aes(x = flipper_length_mm, y = bill_length_mm, color = species, shape=k3)) +
    geom_point(size=3)

library(cluster)
png(filename="sil.png", width=900, height = 580)
plot(silhouette(k3$cluster, dist=d), col=palette()[sort(k3$cluster)])
dev.off()

xy <- data.frame(x=rnorm(10), y=rnorm(10))
xy$lab <- as.character(1:10)

ggplot(xy, aes(x, y)) +
    geom_point(size=5) +
    geom_point(data=xy[2,], color=palette()[2], size=6) +
    geom_point(data=xy[5,], color=palette()[3], size=6) +
    geom_point(data=xy[6,], color=palette()[4], size=6) +
    geom_point(data=xy[8,], color=palette()[4], size=6) +
    geom_text(aes(label=lab))

library(scran)
library(igraph)
g <- buildSNNGraph(t(xy), k=3)
plot(g, layout=as.matrix(xy[,1:2]))

ggplot(xy, aes(x, y)) +
    geom_point(size=5, color=palette()[7]) +
    geom_text(aes(label=lab))
