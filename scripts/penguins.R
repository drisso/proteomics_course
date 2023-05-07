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
