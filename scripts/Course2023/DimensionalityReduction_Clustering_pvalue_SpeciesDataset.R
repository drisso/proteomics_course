## ----packages------------------------------------------------------------
## the main MALDIquant package
library("MALDIquant")
## the import/export routines for MALDIquant
library("MALDIquantForeign")

#install.packages("devtools")
library(devtools)
#install_github("sgibb/MALDIquantExamples")
## example data
library("MALDIquantExamples")

## ----import--------------------------------------------------------------
spectra <- import(getPathSpecies(), verbose=FALSE)

## ----inspect------------------------------------------------------------------
any(sapply(spectra, isEmpty))
all(sapply(spectra, isRegular))

## ----trim----------------------------------------------------------------
spectra <- trim(spectra)

## ----seed------------------------------------------------
set.seed(123)

## ----plot----------------------------------------------------------------
plot(spectra[[1]])
idx <- sample(length(spectra), size=2)
plot(spectra[[idx[1]]])
plot(spectra[[idx[2]]])

## ----tI------------------------------------------------------------------
spectra <- transformIntensity(spectra, method="sqrt")

## ----Plot----------------------------------------------------------------
plot(spectra[[1]], type="b",
     xlim=c(2235.3, 2252.0), ylim=c(45, 100))
abline(h=72, col=4, lty=2)
plot(spectra[[1]], type="b",
     xlim=c(11220, 11250), ylim=c(24, 40))
abline(h=32, col=4, lty=2)

## ----sM------------------------------------------------------------------
spectra <- smoothIntensity(spectra, method="SavitzkyGolay",
                           halfWindowSize=10)

## ----eB------------------------------------------------------------------
## define iteration steps: 25, 50, ..., 100
iterations <- seq(from=25, to=100, by=25)
## define different colors for each step
col <- rainbow(length(iterations))

plot(spectra[[1]], xlim=c(2000, 12000))

## draw different baseline estimates
for (i in seq(along=iterations)) {
  baseline <- estimateBaseline(spectra[[1]], method="SNIP",
                               iterations=iterations[i])
  lines(baseline, col=col[i], lwd=2)
}

legend("topright", legend=iterations, col=col, lwd=1)

## ----rB------------------------------------------------------------------
spectra <- removeBaseline(spectra, method="SNIP",
                          iterations=25)
plot(spectra[[1]])

## ----cI------------------------------------------------------------------
spectra <- calibrateIntensity(spectra, method="TIC")

## ----aS------------------------------------------------------------------
spectra <- alignSpectra(spectra)

## ----metadata------------------------------------------------------------
metaData(spectra[[1]])$spot

## ----spots---------------------------------------------------------------
spots <- sapply(spectra, function(x)metaData(x)$spot)
species <- sapply(spectra, function(x)metaData(x)$sampleName)
head(spots)
head(species)

## ----average-------------------------------------------------------------
paste0(species,"-",spots)
avgSpectra <-
  averageMassSpectra(spectra, labels=paste0(species, spots))

## ----noise---------------------------------------------------------------
## define snrs steps: 1, 1.5, ... 2.5
snrs <- seq(from=1, to=2.5, by=0.5)
## define different colors for each step
col <- rainbow(length(snrs))

## estimate noise
noise <- estimateNoise(avgSpectra[[1]],
                       method="MAD")
#noise <- estimateNoise(avgSpectra[[1]],
#                       method="SuperSmoother")

plot(avgSpectra[[1]],
     xlim=c(6000, 16000), ylim=c(0, 0.0016))

for (i in seq(along=snrs)) {
  lines(noise[, "mass"],
        noise[, "intensity"]*snrs[i],
        col=col[i], lwd=2)
}
legend("topright", legend=snrs, col=col, lwd=1)

## ----pd------------------------------------------------------------------
peaks <- detectPeaks(avgSpectra, SNR=2, halfWindowSize=10)

## ----pdp-----------------------------------------------------------------
plot(avgSpectra[[1]], xlim=c(6000, 16000), ylim=c(0, 0.0016))
points(peaks[[1]], col="red", pch=4)

## ----bP------------------------------------------------------------------
peaks <- binPeaks(peaks)

## ----fP------------------------------------------------------------------
peaks <- filterPeaks(peaks, minFrequency=0.25)

## ----spots2--------------------------------------------------------------
spots <- sapply(avgSpectra, function(x)metaData(x)$spot)
species <- sapply(avgSpectra, function(x)metaData(x)$sampleName)
species <- factor(species) # convert to factor
table(species)
# (needed later in crossval)

## ----fm------------------------------------------------------------------
#FINAL MATRIX
featureMatrix <- intensityMatrix(peaks, avgSpectra)
rownames(featureMatrix) <- paste(species, spots, sep=".")
colnames(featureMatrix)
# If you want to reduce the number after the comma, use the function round in your original dataset
newcolnames<-as.numeric(colnames(featureMatrix))
newcolnames<-round(as.numeric(colnames(featureMatrix)),2)
newcolnames
as.character(newcolnames)
colnames(featureMatrix)<-as.character(newcolnames)

# 1. Dimensionality reduction

system.time(pca <- prcomp(featureMatrix, scale=TRUE))
# Extract the proportion of variance explained by the first 10 Principal Components
pca_summary <- summary(pca)
pca_summary$importance[,1:10]

library(ggplot2)
library("factoextra")
# "factoextra" a package for quick Principal Component Analysis data visualization 

# The scree plot is used to visualize the importance of each principal component and 
fviz_eig(pca, addlabels = TRUE)

# Graph of the variables
fviz_pca_var(pca, col.var = "black")
fviz_pca_var(pca, col.var = "black", select.var = list(contrib = 10))

# Graph of individuals
fviz_pca_ind(pca, geom="point")
# Change the size of points
fviz_pca_ind(pca, geom="point", pointsize = 4)
fviz_pca_ind(pca, geom="point", pointsize = 4, habillage=species)

# Contribution of each variable 
var <- get_pca_var(pca)
var
library("corrplot")
corrplot(var$cos2, is.corr=FALSE)
corrplot(var$cos2[1:50,1:10], is.corr=FALSE)

library(patchwork)
# The function fviz_contrib() creates a barplot of row/column contributions. 
fviz_contrib(pca, choice = "var", axes = 1:2)
fviz_contrib(pca, choice = "var", axes = 1:2, top=50)
# Contributions of variables to PC1
graph1<-fviz_contrib(pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
graph2<-fviz_contrib(pca, choice = "var", axes = 2, top = 10)
graph1+graph2

# All the information in the same plot
# Change the color by groups, add ellipses
fviz_pca_biplot(pca, label="var", habillage=species,
                addEllipses=TRUE, ellipse.level=0.95, select.var = list(cos2 = 10))
# Control automatically the color of individuals using the cos2
fviz_pca_biplot(pca, label ="var", col.ind="cos2", select.var = list(cos2 = 10)) +
  theme_minimal()
# Control variable colors using their contributions
fviz_pca_var(pca, col.var="contrib", select.var = list(cos2 = 10))+
  scale_color_gradient2(low="white", mid="blue",high="red", midpoint=96) +
  theme_minimal()
# Ops... pay attention to midpoint, must stay in the range of "contrib"
# Control variable colors using their contributions
fviz_pca_var(pca, col.var="contrib", select.var = list(cos2 = 10))+
  scale_color_gradient2(low="white", mid="blue",high="red", midpoint=0.403) +
  theme_minimal()


library(Rtsne)
# TSNE
system.time(tsne <- Rtsne(pca$x[,1:10], pca=FALSE,perplexity = 10))
df <- as.data.frame(cbind(pca$x))
# extract the first two components
df$tsne1_pca <- tsne$Y[,1]
df$tsne2_pca <- tsne$Y[,2]
df$species<-species
# plot
p_tsne_pca<-ggplot(df, aes(tsne1_pca, tsne2_pca, color=species)) +
  geom_point()+theme_classic()
p_pca<-ggplot(df, aes(PC1, PC2, color=species)) +
  geom_point()+theme_classic()

# compare tsne vs pca
p_pca+p_tsne_pca





# 2. Clustering

# Because of curse of dimensionality, better to cluster in reduced space
# simple k-means
table(species)
# We expect 4 clusters due to the fact we know that in the dataset there were 4 clusters
# But if we did not know the labels?
# Try to estimate the right number of clusters with the silhouette plot
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
# What we found? Exactly 4 cluster as the number of species

# K-means with 4 clusters
km_rid <- kmeans(pca$x, centers=4)
km_rid
# Cluster means of the first 10 components
km_rid$centers[,1:10]
df <- as.data.frame(pca$x)
df$kmeans_rid <- as.factor(km_rid$cluster)
# see clustering results of kmeans on dimensionality reduction space
p_pca<-ggplot(df, aes(PC1, PC2, color=kmeans_rid)) +
  geom_point()+theme_classic()
p_pca

# Compare k-means with the number of species
table(species, km_rid$cluster)


# Hierarchical clustering
# First calculate the euclidean distance between all points using the dist() function
d <- dist(pca$x)
# Hierarchical Clustering function
hc <- hclust(d)
plot(as.dendrogram(hc))
# set the heigth to obtain 4 cluster
abline(h = 30, col = "red")
# Cut by height
cutree(hc, h = 30)
# How many spectra fall in each cluster?
table(cutree(hc, h = 30))
# Cut by number of cluster
df$hclust4<- as.factor(cutree(hc, k=4))

#compare with the k-means clustering result
table(df$hclust4,df$kmeans_rid)
p_pca_h<-ggplot(df, aes(PC1, PC2, color=hclust4)) +
  geom_point()+theme_classic()
p_pca+p_pca_h

# Clustering linkage
# Cluster using complete linkage: hclust.complete
hclust.complete <- hclust(d, method = "complete")
# Cluster using average linkage: hclust.average
hclust.average <- hclust(d, method = "average")
# Cluster using single linkage: hclust.single
hclust.single <- hclust(d, method = "single")

# Plot dendrogram of hclust.complete
plot(hclust.complete, main = "Complete")
abline(h = 30, col = "red")
# Plot dendrogram of hclust.average
plot(hclust.average, main = "Average")
# The height of a dendogram can change
abline(h = 30, col = "red")
abline(h = 25, col = "blue")
# Plot dendrogram of hclust.single
plot(hclust.single, main = "Single")
abline(h = 30, col = "red")
abline(h = 25, col = "blue")
abline(h = 20.5, col = "green")

# We want to compare the clustering results that we will obtain
# with the different clustering linkage,
# asking always for 4 cluster.
# Looking at the plot,we can expected different results
df$hclust4_complete<- as.factor(cutree(hclust.complete, k=4))
df$hclust4_average<- as.factor(cutree(hclust.average, k=4))
df$hclust4_single<- as.factor(cutree(hclust.single, k=4))

table(df$hclust4_complete,df$hclust4_single)
table(df$hclust4_complete,df$hclust4_average)
table(df$hclust4_single,df$hclust4_average)
# who are the six samples clustered differently? Wrong assignment by hclust with average linkage
species[which(df$hclust4_complete==3 & df$hclust4_average==1)]
species[which(df$hclust4_single==3 & df$hclust4_average==1)]


table(df$hclust4,df$kmeans_rid)
p_pca_h_complete<-ggplot(df, aes(PC1, PC2, color=hclust4_complete)) +
  geom_point()+theme_classic()
p_pca_h_average<-ggplot(df, aes(PC1, PC2, color=hclust4_average)) +
  geom_point()+theme_classic()
p_pca_h_single<-ggplot(df, aes(PC1, PC2, color=hclust4_single)) +
  geom_point()+theme_classic()
p_pca+p_pca_h_complete+p_pca_h_average+p_pca_h_single


summary(featureMatrix[which(df$hclust4_average==1),])
summary(featureMatrix[which(df$hclust4_average==2),])
summary(featureMatrix[which(df$hclust4_average==3),])
summary(featureMatrix[which(df$hclust4_average==4),])




# 3. From p-value to penalized regression
library("tableone")
###TABELLE
##baseline characteristic
data_prot<-data.frame(featureMatrix)

colnames(data_prot)<-c(as.character(newcolnames))
summary(data_prot)
summary(data_prot[1:10])
data_prot$OUTCOME<-species

dim(data_prot)
#strata per OUTCOME
vars <-c(names(data_prot)[1:448])
strata<-c("OUTCOME")
nonorm<-c(names(data_prot)[1:448])
table1 <- CreateTableOne(vars=vars,data = data_prot,strata=strata,includeNA = T,addOverall = T)
table1 <- print(table1,nonnormal = nonorm,printToggle=FALSE,showAllLevels = T,formatOptions = list(big.mark = ","),explain=T)
View(table1)

# Change format
data_prot<-data_prot[1:448]*100000
data_prot$OUTCOME<-species
vars <-c(names(data_prot)[1:448])
strata<-c("OUTCOME")
nonorm<-c(names(data_prot)[1:448])
# calculate median and IQR, use in print the arguments: nonnormal = nonorm
table1 <- CreateTableOne(vars=vars,data = data_prot,strata=strata,includeNA = T,addOverall = T)
table1 <- print(table1,nonnormal = nonorm,printToggle=FALSE,showAllLevels = T,formatOptions = list(big.mark = ","),explain=T)
View(table1)
# calculate mean and SD, delete in print the arguments: nonnormal = nonorm
table1 <- CreateTableOne(vars=vars,data = data_prot,strata=strata,includeNA = T,addOverall = T)
table1 <- print(table1,printToggle=FALSE,showAllLevels = T,formatOptions = list(big.mark = ","),explain=T)
View(table1)


tableOne <- CreateTableOne(vars=vars,data = data_prot,strata=strata,includeNA = T,addOverall = T)
attributes(tableOne)
attributes(tableOne$ContTable)
attr(tableOne$ContTable, "pValues")
data_pvalue<-attr(tableOne$ContTable, "pValues")
p_value_norm<-data_pvalue$pNormal
p_value_NOnorm<-data_pvalue$pNonNormal
# Adjusted p_value: see the help function to see the possible adjustment methods include in the package.
#
?p.adjust
# The adjustment methods include the Bonferroni correction ("bonferroni") in which 
# the p-values are multiplied by the number of comparisons.
p_value_norm_adj_BON<-p.adjust(p_value_norm, method = "bonferroni", n = length(p_value_norm))
# The Holm method was designed to give strong control of the family-wise error rate
p_value_norm_adj_HOLM<-p.adjust(p_value_norm, method = "holm", n = length(p_value_norm))
# The Benjamini-Hochberg method control the false discovery rate, 
# the expected proportion of false discoveries amongst the rejected hypotheses. 
# The false discovery rate is a less stringent condition than the family-wise error rate, 
# so these methods are more powerful than the others.
p_value_norm_adj_BH<-p.adjust(p_value_norm, method = "BH", n = length(p_value_norm))
# See how much changes the number of p-value statistically significant
length(which(p_value_norm<0.05))
length(which(p_value_norm_adj_BON<0.05))
length(which(p_value_norm_adj_HOLM<0.05))
length(which(p_value_norm_adj_BH<0.05))

#...We can do the same analysis with non-normal test..
length(which(p_value_NOnorm<0.05))

# Create a table with all the p_value under different adjustment to see:
# - which variables are the most differentially expressed among group
# - how much change p-value under different adjustment methods
data_pvalue_all<-data.frame(p_value_norm,p_value_norm_adj_BON,p_value_norm_adj_BH,p_value_norm_adj_HOLM,p_value_NOnorm)
rownames(data_pvalue_all)<-names(data_prot)[1:448]
View(data_pvalue_all)
data_pvalue_all_order<-data_pvalue_all[order(data_pvalue_all$p_value_norm, decreasing = FALSE), ]
View(data_pvalue_all_order)


# When you have multiple variables in your logistic regression model, 
# it might be useful to find a reduced set of variables resulting to an optimal performing model.
# Penalized logistic regression imposes a penalty to the logistic model for having too many variables.
# This results in shrinking the coefficients of the less contributive variables toward zero. 
# This is also known as regularization.
# The most commonly used penalized regression include:
# - ridge regression: variables with minor contribution have their coefficients close to zero. 
#   However, all the variables are incorporated in the model. This is useful when all variables need to be incorporated in the model according to domain knowledge.
# - lasso regression: the coefficients of some less contributive variables are forced to be exactly zero. 
#   Only the most significant variables are kept in the final model.
# - elastic net regression: the combination of ridge and lasso regression. 
#   It shrinks some coefficients toward zero (like ridge regression) and set some coefficients to exactly zero (like lasso regression)

library(caret)
library(glmnet)
# Computing penalized logistic regression
# Data preparation
# The R function model.matrix() helps to create the matrix of predictors 
# and also automatically converts categorical predictors to appropriate dummy variables, 
# which is required for the glmnet() function.

# Dummy code categorical predictor variables
x <- model.matrix(~., data_prot)
x <- model.matrix(~., data_prot)[,-c(1,450:452)]
# outcome (class), must be transformed in numeric
y <- as.numeric(data_prot$OUTCOME)
# R functions
# We will use the R function glmnet() [glmnet package] 
# for computing penalized logistic regression.

#The simplified format is as follow:
#glmnet(x, y, family = "multinomial", alpha = 1, lambda = NULL)

# x: matrix of predictor variables
# y: the response or outcome variable, which is a binary variable.
# family: the response type. Use “binomial” for a binary outcome variable
# alpha: the elasticnet mixing parameter. Allowed values include:
#  “1”: for lasso regression
#  “0”: for ridge regression
# a value between 0 and 1 (say 0.3) for elastic net regression.
# lamba: a numeric value defining the amount of shrinkage. Should be specify by analyst.
# In penalized regression, you need to specify a constant lambda to adjust the amount of the coefficient shrinkage. 
# The best lambda for your data, can be defined as the lambda that minimize the cross-validation prediction error rate. 
# This can be determined automatically using the function cv.glmnet().

# ... To performed the LASSO model, which alpha value we must set?...
# See the help function!
?glmnet
# alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.
set.seed(123)
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "multinomial")
plot(cv.lasso)
cv.lasso$lambda.min
cv.lasso$lambda.1se
coef(cv.lasso, cv.lasso$lambda.min)
coef(cv.lasso, cv.lasso$lambda.1se)

# Final model with lambda.min
# Make prediction on test data
x_test <- model.matrix(~., data_prot)[,-c(1,450:452)]
y_multi_pred_class <- as.numeric(predict(cv.lasso, newx = x_test, type = "class", s = cv.lasso$lambda.min))
lasso_class<-y_multi_pred_class
table(lasso_class,species)


# Final model with lambda.1se
# Make prediction on test data
x_test <- model.matrix(~., data_prot)[,-c(1,450:452)]
y_multi_pred_class <- as.numeric(predict(cv.lasso, newx = x_test, type = "class", s = cv.lasso$lambda.min))
y_multi_pred_class_responseProb <- as.numeric(predict(cv.lasso, newx = x_test, type = "response", s = cv.lasso$lambda.min))
lasso_class<-y_multi_pred_class
table(lasso_class,species)
