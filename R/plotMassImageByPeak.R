###Image re-construction from a compressed img and a given ion
buildImageByPeak<-function(img, mass.peak, tolerance=0.25, selectedPixels = NULL, NormCoefs = NULL, rotate = 0)
{
  #Grab only MassPeaks Lists
  if(typeof(img) != "list")
  {
  stop("Error in buildImageByPeak(), img is not a list")
  }

  img_width<-img$size["x"]
  img_height<-img$size["y"]

  #Auto fill selectedPixels if is empty
  if(is.null(selectedPixels))
  {
    selectedPixels <- 1:nrow(img$pos)
  }

  #Normalization coeficients, 1 if null
  if(is.null(NormCoefs))
  {
    NormCoefs<-rep(1, nrow(img$pos))
  }

  l1 <-which.min( abs( img$mass - ( mass.peak - tolerance ) ) )
  l2 <-which.min( abs( img$mass - ( mass.peak + tolerance ) ) )
#   cat(paste("Sel mass:",mass.peak,"\n"))
#   cat(paste("Sel tol:",tolerance,"\n"))
#   cat(paste("l1:",l1,"\n"))
#   cat(paste("l2:",l2,"\n"))

  l1 <- l1[1] #Store only the first element
  l2 <- l2[length(l2)] #store only the last element
  mass.peak <- round(mean(c(img$mass[l2] , img$mass[l1])), digits = 4)
  tolerance <- round(0.5*(img$mass[l2] - img$mass[l1]), digits = 4)

#   cat(paste("mz1:",img$mass[l1],"\n"))
#   cat(paste("mz2:",img$mass[l2],"\n"))
#   cat(paste("Real Mass:", mass.peak,"\n"))
#   cat(paste("Real Tol:", tolerance,"\n"))

  z_indexes <- l1:l2
  if(length(z_indexes) == 1)
  {
    z_indexes <- c(z_indexes,z_indexes)
  }
#   cat("z_indexes:")
#   print(z_indexes)


  #Pixel array loaded by HDD reading chunks
  zplots<-matrix(0, nrow=img_width, ncol=img_height) #Now I'm using a zero instead of NA to display a completely black background
  xy_i<-1
  for(icube in 1:length(img$data))
  {
    #iSelInCube <-  which((selectedPixels - (icube -1) * nrow(img$data[[1]])) %in% 1:nrow(img$data[[icube]]))
    #iSelInCube <- iSelInCube - (icube -1) * nrow(img$data[[1]]) #Substract the offset to get cube coordinate space


    SelinCubeCoords<-(selectedPixels - (icube -1) * nrow(img$data[[1]]))
    SelinCubeCoords<-SelinCubeCoords[which(SelinCubeCoords >= 1 & SelinCubeCoords <= nrow(img$data[[icube]]))]

    iSelInCube <-  which((1:nrow(img$data[[icube]])) %in% SelinCubeCoords )
#     cat("selectedPixels:\n")
#     print(selectedPixels)
#     cat("iSelInCube:\n")
#     print(iSelInCube)
#     cat("in cube:\n")
#     print(icube)
#     cat("-----------\n")
    if(length(iSelInCube) > 0)
    {
      if(length(iSelInCube) > 1)
      {
        #valpixel<-apply(img$data[[icube]][iSelInCube, z_indexes], 1, max)
        valpixel<-apply(img$data[[icube]][iSelInCube, z_indexes], 1, mean)
      }
      else
      {
        #valpixel<-max(img$data[[icube]][iSelInCube, z_indexes])
        valpixel<-mean(img$data[[icube]][iSelInCube, z_indexes])
      }
      for(i in 1:length(valpixel))
      {
        x<-img$pos[ selectedPixels[xy_i], "x" ]
        y<- 1 + img_height - img$pos[ selectedPixels[xy_i], "y" ] #Flip image verticaly to match ploting
        zplots[x,y]<- valpixel[i] * NormCoefs[selectedPixels[xy_i]] #Draw pixel apply normalization
        xy_i<- xy_i+1
      }
    }
  }

  #Apply rotation
  if(rotate >= 90 && rotate < 180)
  {
    zplots <- apply(zplots, 1, rev)
  }
  if(rotate >= 180 && rotate < 270)
  {
    zplots <- apply(apply(zplots, 1, rev), 1, rev)
  }
  if(rotate >= 270 && rotate < 360)
  {
    zplots <- t(apply(zplots, 2, rev))
  }

  #Return zplots matrix and some metadata in a list
  list(img.data=zplots, mass.peak = mass.peak, tolerance = tolerance)
}

#Image resize
DoubleImgResolution = function(im) {
  # initial width/height
  w.in <- nrow(im)
  h.in <- ncol(im)
  w.out <- 2*w.in - 1
  h.out <- 2*h.in - 1

  # Create empty matrix
  im.out <- matrix(0, nrow =w.out, ncol=h.out )

  #The original part
  im.out[seq(1,w.out,2), seq(1,h.out,2)]<-im

  #Filling the gaps by rows
  xi<-which((1:w.out)%%2 != 0 )
  yi<-which((1:h.out)%%2 == 0 )
  im.out[xi, yi] <- 0.5*im.out[ xi, yi -1] + 0.5*im.out[ xi , yi + 1]

  #Filling the gaps by cols
  xi<-which((1:w.out)%%2 == 0 )
  yi<-which((1:h.out)%%2 != 0 )
  im.out[xi, yi] <- 0.5*im.out[ xi - 1, yi] + 0.5*im.out[ xi + 1, yi]

  #Fill the corss gaps
  xi<-which((1:w.out)%%2 == 0 )
  yi<-which((1:h.out)%%2 == 0 )
  im.out[xi, yi] <- 0.125*im.out[ xi - 1, yi -1] + 0.125*im.out[ xi , yi -1] + 0.125*im.out[ xi + 1, yi -1] + 0.125*im.out[ xi - 1, yi ] + 0.125*im.out[ xi + 1, yi ] + 0.125*im.out[ xi - 1, yi +1] + 0.125*im.out[ xi , yi +1] + 0.125*im.out[ xi + 1, yi +1]

  return(im.out)
}


###Plot an image in matrix format
plotMassImage<-function(in_img, useColors=TRUE, smoothing=0.1, screenWidth = 100, screenHeight = 100 , XResLevel = "x1")
{
  #Comute data range auto as maximum spike
  data_range<-max(in_img[[1]] , na.rm=TRUE) ####Remove NA's to compute max , na.rm=TRUE)

  #Double the imatge resolution artificialy to improve visualitzation
  if(XResLevel == "x2")
  {
    in_img[[1]]<-DoubleImgResolution(in_img[[1]])
  }
  if(XResLevel == "x4")
  {
    in_img[[1]]<-DoubleImgResolution(in_img[[1]])
    in_img[[1]]<-DoubleImgResolution(in_img[[1]])
  }

  #Create the image in a regular 3D Grid
  img_width<-nrow(in_img[[1]])
  img_height<-ncol(in_img[[1]])
  img <- list(x=seq(from=1,to=img_width ,by=1), y=seq(from=1,to=img_height,by=1), z=in_img[[1]])

  #Create an image to ilustrate de control palette of 100 points aprox
  d<-seq(from=0, to=data_range, by=data_range/100)
  legend_data<-matrix(data=d, nrow= ceiling(length(d)*0.1), ncol=length(d), byrow = T)
  imgLegend<- list(x= 0.05* max(d) * seq(0,nrow(legend_data)-1, by = 1)/(nrow(legend_data)-1), y=d, z=legend_data)

  #Create a custom color palette
  if(useColors)
  {
  rcolor<-rev(rainbow(n=81, start=0.85, end=0.7))
	bcolor<-colorRampPalette(c(rgb(0,0,0,0), rcolor[1]))(20)
	colors<-c(bcolor,rcolor)
  }
  else
  {
    colors<-gray(seq(from=0, to=1, by=0.01))
  }

  #Use a black background
  par(bg = "black", fg =  "white", col.lab="white", col.axis = "white", col.main = "white", col.sub = "white")

  #Layout 2 images in one plot
  layout(matrix(c(2,1), 1, 2, byrow = TRUE), widths=c(7,1))

  #Plot Color Legend
  par(mar = c(2, 0, 2, 2))
  image(x=imgLegend,col=colors, useRaster=TRUE, xaxt="n", yaxt="n", asp = 1, bty = "n")
  axis(side=4, cex.axis = 0.7, pos = max(imgLegend$x))

  #Plot the image
  par(mar=c(2,1,2,0), mgp = c(2,0.2,0))

  #Smoothing
  img<-image.smooth(x=img, theta = smoothing)
  image(x=img,col=colors, useRaster=TRUE, xaxt="n", yaxt="n", cex.axis = 0.7, asp=1, bty = "n")
  Xresdivider<-switch(XResLevel, "x1" = 1, "x2" = 2, "x4" = 4)
  xAxis<- seq(0, (img_width + 1), by = ceiling(img_width/20))
  xLabels <- xAxis/Xresdivider
  yAxis<- seq(0, (img_height + 1), by = ceiling(img_height/20))
  yLabels <- yAxis/Xresdivider
  axis(side=2, tck = -0.015, cex.axis = 0.7, pos = 0, at = yAxis, labels = yLabels) #Y axes
  axis(side=1, tck = -0.015, cex.axis = 0.7, pos = 0, at = xAxis, labels = xLabels) #X axes
  title(paste("m/z",in_img[["mass.peak"]] , "±", in_img[["tolerance"]]), cex.main = 1.0, bty = "n")
}

###Hard Thresholding an image
thresholdImage<-function(img, Threshold)
{
  #Grab only MassPeaks Lists
  if(class(img[[1]]) != "matrix")
  {
    stop("Error in thresholdImage(), img is not a matrix")
  }

  #Thresholding
  img_width<-nrow(img[[1]])
  img_height<-ncol(img[[1]])
  out_img<-matrix(nrow=img_width, ncol=img_height)
  for(i in 1:img_width)
  {
    for(j in 1:img_height)
    {
      if(img[["img.data"]][i,j]>= Threshold)
      {
        out_img[i,j]<-1
      }
      else
      {
        out_img[i,j]<-0
      }
    }
  }

  #Return zplots matrix and some metadata in a list
  list(img.data=out_img, mass.peak = img[["mass.peak"]], tolerance = img[["tolerance"]])
}

###Image ploting by a peak
plotMassImageByPeak<-function(img, mass.peak, tolerance=0.25, useColors=TRUE, smoothFactor = 0.1, selectedPixels = NULL, NormalizationCoefs = NULL, rotation = 0, ...)
{
  image<-buildImageByPeak(img=img, mass.peak= mass.peak, tolerance=tolerance, selectedPixels, NormCoefs =  NormalizationCoefs, rotate = rotation)
  plotMassImage(in_img=image, useColors=useColors, smoothing = smoothFactor, ...)
  return(list(selMz = image$mass.peak, selTol = image$tolerance))
}

###Auto Thresholding image using K-Mean clustering over image
autoThresholdImage<-function(img)
{
  rowCount<- nrow(img[[1]]) * ncol(img[[1]])
  colCount<-1

  kmat<-matrix(nrow=rowCount, ncol=colCount)
  for(c in 1:ncol(img[[1]]))
  {
    for(r in 1:nrow(img[[1]]))
    {
      #Col1 is Ion mass
      kmat[r + (c - 1) * nrow(img[[1]]),1]<-img[[1]][r,c]
    }
  }

  #Clustering with  clara()
  clustData<-clara(x=kmat, k=2)
  print(clustData[["clusinfo"]])

  binaryImg<-matrix(nrow=nrow(img[[1]]), ncol=ncol(img[[1]]))
  for(c in 1:ncol(img[[1]]))
  {
    for(r in 1:nrow(img[[1]]))
    {
      binaryImg[r,c]<-clustData[["clustering"]][r + (c - 1) * nrow(img[[1]])]
    }
  }

  #Return binaryImgts matrix and some metadata in a list
  list(img.data=binaryImg, mass.peak = img[["mass.peak"]], tolerance = img[["tolerance"]])
}

#PlotZObject, I put this here cause i don't know where to put it....
plotImageZ<-function(Z, main_title = "", XResLevel = "x1")
{
  Z_hr<-Z

  if(XResLevel == "x2")
  {
    Z_hr<-DoubleImgResolution(Z_hr)
  }

  if(XResLevel == "x4")
  {
    Z_hr<-DoubleImgResolution(Z_hr)
    Z_hr<-DoubleImgResolution(Z_hr)
  }

  rcolor<-rev(rainbow(n=81, start=0.85, end=0.7))
  bcolor<-colorRampPalette(c(rgb(0,0,0,0), rcolor[1]))(20)
  mYcolors<-c(bcolor,rcolor)

  #Use a black background
  par(bg = "black", fg =  "white", col.lab="white", col.axis = "white", col.main = "white", col.sub = "white")

  #Layout 2 images in one plot
  layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(7,1))

  #Plot the image
  par(mar=c(2,1,2,0), mgp = c(2,0.2,0))
  image(1:nrow(Z_hr), 1:ncol(Z_hr), Z_hr, ylab = "", xlab = "", col = mYcolors, asp=1, cex.axis = 0.7)
  title(main_title)

  #Plot Color Legend
  d<-seq(min(Z, na.rm = T), max(Z, na.rm = T), by = (max(Z, na.rm = T) - min(Z, na.rm = T))/100)
  legend_data<-matrix(data=d, nrow= ceiling(length(d)*0.1), ncol=length(d), byrow = T)
  imgLegend<- list(x= 0.05* max(d) * seq(0,nrow(legend_data)-1, by = 1)/(nrow(legend_data)-1), y=d, z=legend_data)
  par(mar = c(2, 0, 2, 2))
  image(x=imgLegend,col=mYcolors, xaxt="n", yaxt="n", asp = 1, bty = "n")
  axis(side=4, cex.axis = 0.7, pos = max(imgLegend$x))
}
