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

## ----inspect------------------------------------------------------------------
any(sapply(spectra, isEmpty))
all(sapply(spectra, isRegular))

## ----trim----------------------------------------------------------------
## choose largest overlapping mass range for all spectra
spectra <- trim(spectra)

## ----seed------------------------------------------------
set.seed(123)

## ----plot----------------------------------------------------------------
idx <- sample(length(spectra), size=2)
plot(spectra[[idx[1]]])
plot(spectra[[idx[2]]])

## ----tI------------------------------------------------------------------
spectra <- transformIntensity(spectra, method="sqrt")

## ----sM------------------------------------------------------------------
spectra <- smoothIntensity(spectra, method="SavitzkyGolay",
                           halfWindowSize=20)

## ----eB------------------------------------------------------------------
baseline <- estimateBaseline(spectra[[1]], method="SNIP",
                             iterations=150)
plot(spectra[[1]])
lines(baseline, col="red", lwd=2)

## ----rB------------------------------------------------------------------
spectra <- removeBaseline(spectra, method="SNIP",
                          iterations=150)
plot(spectra[[1]])

## ----cI------------------------------------------------------------------
spectra <- calibrateIntensity(spectra, method="TIC")

## ----aS------------------------------------------------------------------
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

## ----dP------------------------------------------------------------------
peaks <- detectPeaks(avgSpectra, SNR=2, halfWindowSize=20)

## ----plot-----------------------------------------------------------------
plot(avgSpectra[[1]], xlim=c(4000, 5000), ylim=c(0, 0.002))
points(peaks[[1]], col="red", pch=4)

## ----bP------------------------------------------------------------------
peaks <- binPeaks(peaks)

## ----fP------------------------------------------------------------------
peaks <- filterPeaks(peaks, minFrequency=c(0.5, 0.5),
                     labels=avgSpectra.info$health,
                     mergeWhitelists=TRUE)

## ----fM------------------------------------------------------------------
#FINAL MATRIX
featureMatrix <- intensityMatrix(peaks, avgSpectra)
rownames(featureMatrix) <- avgSpectra.info$patientID
