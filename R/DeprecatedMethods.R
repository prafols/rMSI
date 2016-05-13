#Deprecated methods (starting with . to hide it from user)
#Deprectaed but may introduce some nice ideas.

#Convert MyImg Object to jpeg compressed mz-sliced image list (this could be a great trick to dramatically reduce dataset size yet keep relevant data)
.CompressMyImg2JpegSlices<-function(img_data, compression_quality = 0.9)
{
  #Generate all image slices
  dataCube<-list()
  for(mz_i in 1:length(img_data$mass))
  {
    img_slice<-matrix(rep(0, img_data$size["x"] * img_data$size["y"]), nrow = img_data$size["x"], ncol = img_data$size["y"])

    #Max in this slice
    max_spc<-max(unlist(lapply(img_data$data, function(x) { x$intensity[mz_i]  })))

    for(i in 1:length(img_data$data))
    {
      img_slice[img_data$data[[i]]$X, img_data$data[[i]]$Y]<- img_data$data[[i]]$intensity[mz_i]/max_spc
    }
    dataCube[[mz_i]]<-writeJPEG(img_slice, raw(), quality = compression_quality)
  }

  return(dataCube)
}

#Convert MyImg Object to MALDIquant MassSpectrum object
.ConvertMyImg2MassSpectrum<-function(img_data){
  dataCube<-list()

  #Generate img in MALDIquant format
  for(i in 1:length(img_data$data))
  {
    dataCube[[i]]<-createMassSpectrum(mass = img_data$mass, intensity = img_data$data[[i]]$intensity)
    dataCube[[i]]@metaData$imaging$pos["x"]<-img_data$data[[i]]$X
    dataCube[[i]]@metaData$imaging$pos["y"]<-img_data$data[[i]]$Y
    dataCube[[i]]@metaData$imaging$size["x"]<- img_data$size["x"]
    dataCube[[i]]@metaData$imaging$size["y"]<- img_data$size["y"]
  }
  return(dataCube)
}

#Convert MassSpectrum object to MyImg Object
.ConvertMassSpectrum2MyImg<-function(img_data)
{
  dataCube<-list()
  for(i in 1:length(img_data))
  {
    X<-img_data[[i]]@metaData$imaging$pos["x"]
    Y<-img_data[[i]]@metaData$imaging$pos["y"]
    spect<-img_data[[i]]@intensity
    spot_full<-list(X, Y, spect)
    names(spot_full)<-c("X","Y","intensity")
    dataCube[[length(dataCube)+1]]<-spot_full
  }

  mz_axis<-img_data[[1]]@mass
  x_size <- img_data[[1]]@metaData$imaging$size["x"]
  y_size <- img_data[[1]]@metaData$imaging$size["y"]
  dataCube<-list(mass = mz_axis, size = c(x_size, y_size), data = dataCube)

  return(dataCube)
}
