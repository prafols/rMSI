  # TODO: Add comment
#
# Author: sapista
###############################################################################

#Convert MyImg Object to jpeg compressed mz-sliced image list
CompressMyImg2JpegSlices<-function(img_data, compression_quality = 0.9)
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
ConvertMyImg2MassSpectrum<-function(img_data){
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
ConvertMassSpectrum2MyImg<-function(img_data)
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

##Stores a directory structure of raw_data
dataListing<-function(raw_data_folder, sample_spectra_name)
{
  topDataDir<-dir(raw_data_folder, full.names = F)
  if( sample_spectra_name %in% topDataDir )
  {
    #No sub-data dirs, spectra is directly in raw_data_folder
    rawDirs<-dir(raw_data_folder, full.names = T)
    names(rawDirs)<-topDataDir
  }
  else
  {
    #There exists sub-data dirs
    rawDirs<-unlist(lapply(topDataDir ,function(x) dir(file.path(raw_data_folder, x), full.names = T)))
    names(rawDirs)<-unlist(lapply(topDataDir ,function(x) dir(file.path(raw_data_folder, x), full.names = F)))
  }
  return(rawDirs)
}

#This is the top level function for RAW datacuve import
#RAW Data strcutre matrix (stored in HDD), nrow = npixels, ncol = length(mz_axis)
# | mz1 | mz2 | ... | mzN |
# +-----+-----+-----+-----+
# |     |     | ... |     |
# +-----+-----+-----+-----+
# |     |     | ... |     |
# +-----+-----+-----+-----+
# |     |     | ... |     |
# +-----+-----+-----+-----+
# |     |     | ... |     |
# +-----+-----+-----+-----+
#
#POS Data matrix (stored in RAM)
# | X  | Y  |
# +----+----+
# | x1 | y1 |
# +----+----+
# | x2 | y2 |
# +----+----+
# | .. | .. |
# +----+----+
# | xM | yM |
# +----+----+
#
#MZ_AXIS -> vector (stored in RAM)
importBrukerXmassImg<-function(raw_data_folder, xml_file, ff_data_folder, max_ff_file_size_MB = 50)
{
#1- Read XML input data to get a vector of spectrums
	spectraList<-readXmlSpectraList(xml_file)
  pb<-txtProgressBar(min = 0, max = length(spectraList), style = 3 )
	#Store directory structure to avoid wasting time in dir() fuctions
	rawDirs<-dataListing(raw_data_folder, spectraList[1])

#2- Read Spectrum acqu file and generate m/z vector
	mz_axis<-generateMzAxisFromAcquFile(file.path(rawDirs[spectraList[1]], "1", "1SRef","acqu"))

#3- Read Spectrum FID file and fill data structure if not zero intensity
  dataPos <- matrix(NA, ncol = 2, nrow = length(spectraList))
  colnames(dataPos)<-c("x","y")

  #Do not create to big files (> 50 MB by default)
  max_nrow <- floor((max_ff_file_size_MB*1024*1024)/(4*length(mz_axis)))
  max_nrow <- min(max_nrow, floor(.Machine$integer.max / (length(mz_axis))))

  dataCube<-list()
  for(i in 1:ceiling(length(spectraList)/ max_nrow))
  {
    nrows<-min(max_nrow, (length(spectraList) - sum(unlist(lapply(dataCube, nrow)))))
	  dataCube[[i]]<-ff::ff(vmode = "integer", dim = c(nrows, length(mz_axis)), filename = file.path(ff_data_folder, paste("ramdisk",i,".dat",sep = "")))
    names(dataCube)[i]<-paste("ramdisk",i,".dat",sep = "")
  }

  for(i in 1:length(spectraList))
  {
    setTxtProgressBar(pb, i)
    i_cube<-(1+((i-1) %/% max_nrow))
    dataCube[[i_cube]][(i - (i_cube -1) * max_nrow), ]<-readFidFile(file.path(rawDirs[spectraList[i]], "1", "1SRef","fid"), length(mz_axis))

    #Extract X Y Coords
    dataPos[i,"x"]<-as.integer(strsplit(strsplit(spectraList[i],"X")[[1]][2], "Y")[[1]][1])
    dataPos[i,"y"]<-as.integer(strsplit(strsplit(spectraList[i],"Y")[[1]][2], "X")[[1]][1])
  }
  close(pb)

#4- Calc offsets and subtract it
  x_offset<-min(dataPos[,"x"])
  y_offset<-min(dataPos[,"y"])
  for(i in 1:length(spectraList))
  {
    dataPos[i, "x"] <- dataPos[i, "x"] - x_offset + 1
    dataPos[i, "y"] <- dataPos[i, "y"] - y_offset + 1
  }

#5- Compute Motor coords range
  x_size<-max(dataPos[,"x"])
  y_size<-max(dataPos[,"y"])

#6- Map MALDI motor coords to image cords (1-pixels steps)
  #It is important to map MALDI motors coords to image coords.
  #Otherwise, null extra pixels may be added leading to bad reconstruction
  px_map <- matrix( 0, nrow = x_size, ncol = y_size)
  for(i in 1:nrow(dataPos))
  {
    xi <- dataPos[i, "x"]
    yi <- dataPos[i, "y"]
    px_map[xi, yi]<- i
  }

  colNull <- which( base::colSums(px_map) == 0)
  rowNull <- which( base::rowSums(px_map) == 0)
  remap<-FALSE
  if( length(colNull) > 0 && length(rowNull) > 0 )
  {
    px_map_ <- px_map[ -rowNull , -colNull ]
    remap<-TRUE
  }
  if( length(colNull) > 0 && length(rowNull) == 0 )
  {
    px_map_ <- px_map[ , -colNull ]
    remap<-TRUE
  }
  if( length(colNull) == 0 && length(rowNull) > 0 )
  {
    px_map_ <- px_map[ -rowNull , ]
    remap<-TRUE
  }

  if(remap)
  {
    for(ix in 1:nrow(px_map_))
    {
      for(iy in 1:ncol(px_map_))
      {
        if(px_map_[ix, iy] > 0)
        {
          dataPos[px_map_[ix, iy], "x"] <- ix
          dataPos[px_map_[ix, iy], "y"] <- iy
        }
      }
    }

    x_size <- max(dataPos[,"x"])
    y_size <- max(dataPos[,"y"])
  }

#7- Return dataCube, mz_axis, xsize, ysize as a list of elements
	return(list(mass = mz_axis, size = c(x = x_size, y = y_size), data = dataCube, pos=dataPos))
}

#Read XML input data
#Returns a vector of selected spectrums in a given ROI
readXmlSpectraList<-function(xml_file){
	xml_data<-xmlToList(xml_file)
	#spots<-xml_data$Class[[1]]["Spot"]
	spots<-c()
	for (i in 1:length(xml_data$Class))
	{
		if(!is.na(xml_data$Class[[i]]["Spot"]))
		{
			spots<-c(spots,xml_data$Class[[i]]["Spot"])
		}
	}
	names(spots)<-NULL
	spots
}

#Read Spectrum acqu file
#Returns a vector of all possible m/z values
generateMzAxisFromAcquFile<-function(acqu_file){
	acqu_con<-file(acqu_file, open="r")
	file_lines<-readLines(acqu_con)
	close(acqu_con)

	#Grab Mz axis generation data
	td<-as.numeric(strsplit(grep("##\\$TD=", file_lines, value=T),"=")[[1]][2])
	dw<-as.numeric(strsplit(grep("##\\$DW=", file_lines, value=T),"=")[[1]][2])
	delay<-as.numeric(strsplit(grep("##\\$DELAY=", file_lines, value=T),"=")[[1]][2])
	ml1<-as.numeric(strsplit(grep("##\\$ML1=", file_lines, value=T),"=")[[1]][2])
	ml2<-as.numeric(strsplit(grep("##\\$ML2=", file_lines, value=T),"=")[[1]][2])
	ml3<-as.numeric(strsplit(grep("##\\$ML3=", file_lines, value=T),"=")[[1]][2])

	tof<-delay + (0:(td-1))*dw
	#tof<-delay + (1:(td))*dw
	A<-ml3
	B<-sqrt(1e+12/ml1)
	C<-ml2 - tof
	mz<-((-B + sqrt((B * B) - (4 * A * C)))/(2 * A))^2

	###TODO: M'esta fallant la conversió a m/z, tinc algu d'error respecte els fitxers ascii de bruker...
	###Proposo seguir implementant la lectura RAW del FID file per poder representar espectres i ja seguire més enedevant amb aquest error d m/z
	###De moment ficu un eix mz per defecte no calculat:
# 	mz_axis_file<-"/home/sapista/R/workspace/test/rp_mz_axis.txt"
# 	mz<-read.table(mz_axis_file, header = F, sep = " ")[,1]

	#return mz axis
	mz
}

#Read FID intensities (Bruker Formated)
#Returns a vector of intensities
readFidFile<-function(fid_file, points)
{
	fid_con<-file(fid_file, open="rb")
	#specVect<-as.double(readBin(fid_con, integer(),n=points ,size = 4, signed=T, endian="little" ))
	specVect<-as.integer(readBin(fid_con, integer(),n=points ,size = 4, signed=T, endian="little" ))
	close(fid_con)
	specVect
}
