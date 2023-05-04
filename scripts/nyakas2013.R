library("devtools")
install_github("sgibb/MALDIquantExamples")

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
