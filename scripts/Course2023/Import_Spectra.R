library("MALDIquant")
library("MALDIquantForeign")

pazienti<-c("I18-18417","290exvivo","I18-3561","I18-4057")
dir<-"C:/Users/Desktop/giulia/dati/"
n<-length(pazienti)
spectra_lista<-rep(0,n)

#for n different mass spectra, each one saved in a txt file
for(int in 1:length(pazienti)){

  filename <-paste (dir,pazienti[int], sep ="")
  spectra_lista[int]<-importTxt(filename, massRange=c(3000,20000))
  
}


#for an imzML file
filename <-paste (dir,"ID_10", sep ="")
spectra_imzml <- importImzMl(filename, massRange=c(3000,20000))
