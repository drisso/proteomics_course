## ----packages------------------------------------------------------------
## the main MALDIquant package
library("MALDIquant")
## the import/export routines for MALDIquant
library("MALDIquantForeign")

## example data
library("MALDIquantExamples")

## ----import--------------------------------------------------------------
## import the spectra
spectra <- import(getPathFiedler2009()["spectra"],
                  verbose=FALSE)

## import metadata
spectra.info <- read.table(getPathFiedler2009()["info"],
                           sep=",", header=TRUE)

## ----reduce--------------------------------------------------------------
isHeidelberg <- spectra.info$location == "heidelberg"

spectra <- spectra[isHeidelberg]
spectra.info <- spectra.info[isHeidelberg,]

## ----qc------------------------------------------------------------------
table(sapply(spectra, length))
any(sapply(spectra, isEmpty))
all(sapply(spectra, isRegular))

## ----trim----------------------------------------------------------------
## choose largest overlapping mass range for all spectra
spectra <- trim(spectra)

## ----plotseed, echo=FALSE------------------------------------------------
set.seed(123)

## ----plot----------------------------------------------------------------
idx <- sample(length(spectra), size=2)
plot(spectra[[idx[1]]])
plot(spectra[[idx[2]]])

## ----vs------------------------------------------------------------------
spectra <- transformIntensity(spectra, method="sqrt")

## ----sm------------------------------------------------------------------
spectra <- smoothIntensity(spectra, method="SavitzkyGolay",
                           halfWindowSize=20)

## ----be------------------------------------------------------------------
baseline <- estimateBaseline(spectra[[1]], method="SNIP",
                             iterations=150)
plot(spectra[[1]])
lines(baseline, col="red", lwd=2)

## ----bc------------------------------------------------------------------
spectra <- removeBaseline(spectra, method="SNIP",
                          iterations=150)
plot(spectra[[1]])

## ----cb------------------------------------------------------------------
spectra <- calibrateIntensity(spectra, method="TIC")

## ----pa------------------------------------------------------------------
spectra <- alignSpectra(spectra)

## ----avg-----------------------------------------------------------------
avgSpectra <-
  averageMassSpectra(spectra, labels=spectra.info$patientID)
avgSpectra.info <-
  spectra.info[!duplicated(spectra.info$patientID), ]

## ----noise---------------------------------------------------------------
noise <- estimateNoise(avgSpectra[[1]])
plot(avgSpectra[[1]], xlim=c(4000, 5000), ylim=c(0, 0.002))
lines(noise, col="red")                     # SNR == 1
lines(noise[, 1], 2*noise[, 2], col="blue") # SNR == 2

## ----pd------------------------------------------------------------------
peaks <- detectPeaks(avgSpectra, SNR=2, halfWindowSize=20)

## ----pdp-----------------------------------------------------------------
plot(avgSpectra[[1]], xlim=c(4000, 5000), ylim=c(0, 0.002))
points(peaks[[1]], col="red", pch=4)

## ----pb------------------------------------------------------------------
peaks <- binPeaks(peaks)

## ----pf------------------------------------------------------------------
peaks <- filterPeaks(peaks, minFrequency=c(0.5, 0.5),
                     labels=avgSpectra.info$health,
                     mergeWhitelists=TRUE)

## ----fm------------------------------------------------------------------
#FINAL MATRIX
featureMatrix <- intensityMatrix(peaks, avgSpectra)
rownames(featureMatrix) <- avgSpectra.info$patientID

## ----dda-----------------------------------------------------------------
library("sda")
Xtrain <- featureMatrix
Ytrain <- avgSpectra.info$health
ddar <- sda.ranking(Xtrain=featureMatrix, L=Ytrain, fdr=FALSE,
                    diagonal=TRUE)

## ----ddaresults, echo=FALSE, results="asis"------------------------------
#library(xtable)
#xtable(ddar[1:10, ], booktabs=TRUE)
View(ddar[1:10, ])

## ----hclust--------------------------------------------------------------
distanceMatrix <- dist(featureMatrix, method="euclidean")

hClust <- hclust(distanceMatrix, method="complete")

plot(hClust, hang=-1)

## ----hclustfs------------------------------------------------------------
top <- ddar[1:2, "idx"]

distanceMatrixTop <- dist(featureMatrix[, top],
                          method="euclidean")

hClustTop <- hclust(distanceMatrixTop, method="complete")

plot(hClustTop, hang=-1)

## ----cv------------------------------------------------------------------
library("crossval")
# create a prediction function for the cross validation
predfun.dda <- function(Xtrain, Ytrain, Xtest, Ytest,
                        negative) {
  dda.fit <- sda(Xtrain, Ytrain, diagonal=TRUE, verbose=FALSE)
  ynew <- predict(dda.fit, Xtest, verbose=FALSE)$class
  return(confusionMatrix(Ytest, ynew, negative=negative))
}

# set seed to get reproducible results
set.seed(1234)

cv.out <- crossval(predfun.dda,
                   X=featureMatrix[, top],
                   Y=avgSpectra.info$health,
                   K=10, B=20,
                   negative="control",
                   verbose=FALSE)
diagnosticErrors(cv.out$stat)


### pca
pca <- prcomp(featureMatrix, scale. = TRUE)
plot(pca$x)

library(ggplot2)
theme_set(theme_classic())

df <- as.data.frame(pca$x)
df$status <- avgSpectra.info$health

ggplot(df, aes(PC1, PC2, color=status)) +
    geom_point()

feature_df <- as.data.frame(featureMatrix)
feature_df <- cbind(feature_df, avgSpectra.info)

ggplot(feature_df, aes(x=health, y=`1011.95683040433`)) +
    geom_boxplot()

x <- feature_df$health
y <- feature_df$`1011.95683040433`
mtrue <- tapply(y, x, mean)
diff <- mtrue[1]-mtrue[2]

differences <- sapply(seq_len(10000), function(i) {
    newy <- sample(y)
    m <- tapply(newy, x, mean)
    m[1]-m[2]
})
hist(differences)
abline(v=diff, col=2, lwd=2)

dens <- density(differences, bw=4e-6)
plot(dens, type='h', col=as.numeric(abs(dens$x)>diff)+1, main="Differences")
abline(v=c(diff, -diff), col=2, lwd=2)

suby <- round(y*10^4, 2)[c(1:4, 37:40)]
subx <- x[c(1:4, 37:40)]

df <- data.frame(y = sample(suby), x = subx)
df
tmp <- tapply(df$y, df$x, mean)
tmp[1]-tmp[2]

t.test(y~x, var.equal=TRUE)
fit <- lm(y~x)
summary(fit)

x <- as.numeric(as.factor(x))-1
feature_df$y <- x
colnames(feature_df)[1:166] <- paste0("x", 1:166)
fit <- glm(y~x1+x2+x3+x4+x5, family = binomial(), data=feature_df)
summary(fit)
eta <- predict(fit)
expected <- exp(-eta)/(1+exp(-eta))
table(expected>.5, feature_df$health)

subset <- feature_df[,c(1:166, 172)]
fit <- glm(y~., family = binomial(), data=subset)
summary(fit)

library("crossval")
fun <- function(Xtrain, Ytrain, Xtest, Ytest, negative) {
    train <- cbind(Ytrain, Xtrain)
    names(train)[1] <- "y"
    fit <- glm(y~., family = binomial(), data=train)
    eta <- predict(fit, newdata = Xtest)
    expected <- exp(-eta)/(1+exp(-eta))
    ynew <- as.numeric(expected>.5)
    return(confusionMatrix(Ytest, ynew, negative=negative))
}

cv.out <- crossval(fun,
                   X=subset[,1:5],
                   Y=subset$y,
                   K=10, B=20,
                   negative="control",
                   verbose=FALSE)
diagnosticErrors(cv.out$stat)


library(glmnet)
Xtrain <- as.matrix(subset[,-ncol(subset)])
Xtrain <- scale(Xtrain)
cv <- cv.glmnet(x=Xtrain, y=subset$y, family=binomial())
plot(cv)
fit <- glmnet(x=Xtrain, y=subset$y, family=binomial())
plot(fit, xvar = "lambda")

fit <- glmnet(x=Xtrain, y=subset$y, family=binomial(), lambda=cv$lambda.1se)
coef(fit)
