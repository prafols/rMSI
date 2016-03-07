#Idem k anterior xo sense:
# - selectedPixels: Ho solucionare amb raster::crop raster::zoom i limatant el valor maxim de Z en plot del raster per cada layer RGB
# - rotate: Ho solucionare directament rotant objecte de raster en plotar, no em cal re-generar imatge per plotar
# - Retorna: una llista on
#         img.data correspondra al raster de la imatge,
#         mass.peak
#         tolerance
# El resultat d'aquest metode es guardara en la classe MSImagePlotWidget en una llista de capes R,G,B
.buildImageByPeak<-function(img, mass.peak, tolerance=0.25, NormCoefs = NULL)
{
  #Grab only MassPeaks Lists
  if(typeof(img) != "list")
  {
    stop("Error in rMSI::.buildImageByPeak(), img is not a list")
  }

  #Get the image slice
  img_slice<-builRasterImageFromMass( img, mass.peak, tolerance, "max", NormCoefs) #TODO implement more methods

  #Create the raster
  my_raster <- raster::raster( nrow = ncol(img_slice$pixels), ncol = nrow(img_slice$pixels), xmn= 0, xmx= nrow(img_slice$pixels), ymn= 0, ymx= ncol(img_slice$pixels))
  raster::values(my_raster) <- as.vector(img_slice$pixels)

  #Return zplots matrix and some metadata in a list
  list(raster = my_raster, mass = img_slice$Mass, tolerance = img_slice$Tolerance, cal_resolution = img$pixel_size_um)
}


#Create an empty raster object
.InitRGBEmptyRaster<-function( width, height)
{
  new_raster<-raster::raster( nrow = height, ncol = width, xmn= 0, xmx= width, ymn= 0, ymx= height)
  raster::values(new_raster) <- rep(0, width*height)
  return( list( raster = new_raster, mass = NULL, tolerance = NULL, cal_resolution = NULL))
}


#Normalize to 255 (8 bits per color channel)
.NormalizeTo255 <-function( m, scaling )
{
  maxN<- max(raster::values(m))
  if(maxN > 0)
  {
    raster::values(m) <- scaling*255 * raster::values(m) / (maxN*3) #TODO, check that 3 default scale
    raster::values(m)[raster::values(m)>255]<- 255 #Clipping
  }
  return(m)
}

#Idem k anterior plotMassImageRGB
#Com a param img_RGB es una llista amb atributs R,G,B on cada objecte (layer) es un objecte retornat per builImageByPeak
# XResLevel es ara un integer ja que delego la interpolacio a raster
# All RGB layer must have the same size
.BuildRGBImage <- function( imgR, imgG, imgB, XResLevel = 3 , light = 3)
{
  imgR <- .NormalizeTo255(imgR$raster, light)
  imgG <- .NormalizeTo255(imgG$raster, light)
  imgB <- .NormalizeTo255(imgB$raster, light)

  #Create an RGB image space
  RGB_raster <- raster::addLayer(imgR,imgG,imgB )
  interpolated_raster <- raster::raster( nrow= XResLevel*RGB_raster@nrows, ncol= XResLevel*RGB_raster@ncols, xmn= 0, xmx= RGB_raster@ncols, ymn= 0, ymx= RGB_raster@nrows)
  RGB_raster<-raster::resample(RGB_raster, interpolated_raster)
  raster::values(RGB_raster)[ raster::values(RGB_raster) < 0  ] <- 0 #Values below zero are dube interpolation artifacts, clip it to zero.

  return(RGB_raster)
}

#Remap a intensity vector to HSV coloring in a rainbow function
.ReMappingIntensity2HSV<-function(single_channel_raster, maxN = max(raster::values(single_channel_raster)), value_multiplier = 3)
{
  #Normalize to 1
  if(maxN > 0)
  {
    raster::values(single_channel_raster) <- raster::values(single_channel_raster) / maxN
  }

  #Remapping hue space
  hue_top <- 0.7
  hue_bottom <- 0.85
  hMapped <- (1 + hue_top - hue_bottom) * (-1*raster::values(single_channel_raster)  + 1) + hue_bottom
  over_one <- which(hMapped > 1)
  hMapped[ over_one ] <- hMapped[over_one] -1

  #Remapping value space
  vMapped <- raster::values(single_channel_raster) *value_multiplier
  vMapped[ vMapped > 1] <- 1

  #Creating HSV image
  rgbSpace <- col2rgb(hsv( h =  hMapped, s = rep(1, length(raster::values(single_channel_raster) )), v = vMapped))

  #Layering RGB raster
  R_raster<- single_channel_raster
  G_raster<- single_channel_raster
  B_raster<- single_channel_raster
  raster::values(R_raster) <-  rgbSpace[1, ] #Extract R
  raster::values(G_raster) <-  rgbSpace[2, ] #Extract G
  raster::values(B_raster) <-  rgbSpace[3, ] #Extract B

  #Create an RGB image space
  return (raster::addLayer(R_raster,G_raster,B_raster ))
}

#Build a RGB images raster using rainbow colors from only one raster layer
.BuildSingleIonRGBImage<-function( img,   XResLevel = 3 , global_intensity_scaling_factor = NULL, light = 3)
{
  #Create an RGB image space
  if(is.null(global_intensity_scaling_factor))
  {
    RGB_raster <- .ReMappingIntensity2HSV(img$raster, value_multiplier = light)
  }
  else
  {
    RGB_raster <- .ReMappingIntensity2HSV(img$raster, maxN = global_intensity_scaling_factor, value_multiplier = light)
  }
  interpolated_raster <- raster::raster( nrow= XResLevel*RGB_raster@nrows, ncol= XResLevel*RGB_raster@ncols, xmn= 0, xmx= RGB_raster@ncols, ymn= 0, ymx= RGB_raster@nrows)
  RGB_raster<-raster::resample(RGB_raster, interpolated_raster)
  raster::values(RGB_raster)[ raster::values(RGB_raster) < 0  ] <- 0 #Values below zero are due interpolation artifacts, clip it to zero.

  return(RGB_raster)
}

#roi_rectangle  is c(left, rigth, bottom, top)
.plotMassImageRGB <- function(rasterRGB, cal_um2pixels = 1,  rotation=0, flipV=F, flipH=F, display_axes=T, roi_rectangle = NULL, zoom = F)
{
  img_Xmax <- rasterRGB@extent@xmax
  img_Ymax <- rasterRGB@extent@ymax

  #Crop raster according the zoom window
  if( zoom && !is.null(roi_rectangle) )
  {
    rasterRGB<- raster::crop( rasterRGB, raster::extent( c(roi_rectangle[1] -1, roi_rectangle[2], rasterRGB@extent@ymax - roi_rectangle[4], rasterRGB@extent@ymax - roi_rectangle[3] + 1)))
  }

  #Setting my tricky par values...
  par( bg = "black", fg =  "white", col.lab="white", xaxt="n", yaxt="n", col.axis = "white", col.main = "white", col.sub = "white",
       cex.axis = 0.6, mar = c(1,1,1,1), mgp = c(2, 0.5, 0.5))

  #TODO la ROI no es mostra be amb zooms rotats! Xo oju k ara esta ben quadrat per no-zoom!

  #Apply rotation
  roi_rectangle_pre <- roi_rectangle
  if( rotation == 0 )
  {
    if(!is.null(roi_rectangle))
    {
      roi_rectangle[1] <- roi_rectangle_pre[1] - 1 #Left maps Left
      roi_rectangle[2] <- roi_rectangle_pre[2] #Rigth maps Right
      roi_rectangle[3] <- img_Ymax - roi_rectangle_pre[4] #Bottom  maps Top
      roi_rectangle[4] <- img_Ymax - roi_rectangle_pre[3] + 1 #Top maps Bottom
    }
  }

  if( rotation == 90 )
  {
    if(!is.null(roi_rectangle))
    {
      roi_rectangle[1] <- roi_rectangle_pre[3] - 1 #Left maps Bottom
      roi_rectangle[2] <- roi_rectangle_pre[4] #Rigth maps Top
      roi_rectangle[3] <- roi_rectangle_pre[1] - 1 #Bottom  maps Left
      roi_rectangle[4] <- roi_rectangle_pre[2] #Top maps Right
    }


    rasterRGB <- raster::flip(raster::t(rasterRGB), direction = "y") #90º rotation
  }
  if( rotation == 270 )
  {
    if(!is.null(roi_rectangle))
    {
      roi_rectangle[1] <- img_Ymax - roi_rectangle_pre[4] #Left maps top
      roi_rectangle[2] <- img_Ymax - roi_rectangle_pre[3] + 1 #Rigth maps bottom
      roi_rectangle[3] <- img_Xmax - roi_rectangle_pre[2] #Bottom  maps Right
      roi_rectangle[4] <- img_Xmax - roi_rectangle_pre[1] + 1 #Top maps Left
    }

    rasterRGB <- raster::flip(raster::t(rasterRGB), direction = "x") #270º rotation
  }
  if( rotation == 180 )
  {
    if(!is.null(roi_rectangle))
    {
      roi_rectangle[1] <- img_Xmax - roi_rectangle_pre[2] #Left maps Right
      roi_rectangle[2] <- img_Xmax - roi_rectangle_pre[1] + 1 #Rigth maps Left
      roi_rectangle[3] <- roi_rectangle_pre[3] -1 #Bottom  maps Top
      roi_rectangle[4] <- roi_rectangle_pre[4] #Top maps Bottom
    }

    rasterRGB <- raster::flip(raster::flip(rasterRGB, direction = "y"), direction = "x")#180º rotation
  }

  raster::plotRGB(rasterRGB, axes = T, asp = 1, interpolate = F ) ##TODO testing interpolation

  if(display_axes)
  {
    #Add calibrated axes
    xAxis<- seq(0, rasterRGB@extent@xmax, by = (rasterRGB@extent@xmax/10))
    xLabels <- sprintf( "%0.1f", xAxis * cal_um2pixels)
    yAxis<- seq(0, rasterRGB@extent@ymax, by = (rasterRGB@extent@ymax/10))
    yLabels <- sprintf( "%0.1f", yAxis * cal_um2pixels)
    par(xaxt = "l", yaxt = "l")
    axis(side=2, tck = -0.015, cex.axis = 0.7, pos = 0, at = yAxis, labels = yLabels, las = 1) #Y left axes
    axis(side=4, tck = -0.015, cex.axis = 0.7, pos = rasterRGB@extent@xmax, at = yAxis, labels = yLabels, las = 1) #Y right axes
    axis(side=1, tck = -0.015, cex.axis = 0.7, pos = 0, at = xAxis, labels = xLabels ) #X below axes
    axis(side=3, tck = -0.015, cex.axis = 0.7, pos = rasterRGB@extent@ymax, at = xAxis, labels = xLabels ) #X avobe axes
  }
  else
  {
    #Add calibrated scale bar
    Lp <- 0.05
    Hp <- 0.015
    WpTarget <- 0.15

    #Cal the most elegant nearest value
    legend_possible_values <- as.vector(sapply(10^(1:4), function(x){ x*(1:9) }))
    cal_length <- legend_possible_values[which.min(abs(legend_possible_values - WpTarget*rasterRGB@extent@xmax*cal_um2pixels))]
    Wp <- cal_length/(cal_um2pixels*rasterRGB@extent@xmax)
    yB <- Lp*rasterRGB@extent@ymax
    yT <- (Lp + Hp)*rasterRGB@extent@ymax

    if( rotation == 180  )
    {
      #Avoid overlapping scale and axes
      xL <- rasterRGB@extent@xmin + (Lp + Wp)*(rasterRGB@extent@xmax - rasterRGB@extent@xmin)
      xR <- rasterRGB@extent@xmin + Lp * (rasterRGB@extent@xmax - rasterRGB@extent@xmin)
    }
    else
    {
      xL <- rasterRGB@extent@xmax - (Lp + Wp)*(rasterRGB@extent@xmax - rasterRGB@extent@xmin)
      xR <- rasterRGB@extent@xmax - Lp * (rasterRGB@extent@xmax - rasterRGB@extent@xmin)
    }
    lines( c( xL, xL, xR, xR), c( yB, yT, yT, yB ), col = "white", lwd = 2 )
    text( x = (xL + 0.5*(xR - xL)), y = 1.3*yT, labels = sprintf("%0.0f um", cal_length), col = "white", cex = 0.8, adj = c(0.5,0))

    #Add coors system arrows
    raster_size <- c(rasterRGB@extent@xmax-rasterRGB@extent@xmin, rasterRGB@extent@ymax-rasterRGB@extent@ymin)
    arrow_length <- 0.25*min(raster_size)
    if( rotation == 0  )
    {
      P_0 <- c(rasterRGB@extent@xmin + 0.01*(raster_size[1]), rasterRGB@extent@ymax - 0.01*(raster_size[2]))
      P_X <- c(rasterRGB@extent@xmin + arrow_length, rasterRGB@extent@ymax - 0.01*(raster_size[2]))
      P_Y <- c(rasterRGB@extent@xmin + 0.01*(raster_size[1]), rasterRGB@extent@ymax - arrow_length)
      Txt_Adj <- c(0, 1)
    }
    if( rotation == 90 )
    {
      P_0 <- c(rasterRGB@extent@xmin + 0.01*(raster_size[1]), rasterRGB@extent@ymin + 0.01*(raster_size[2]))
      P_X <- c(rasterRGB@extent@xmin + 0.01*(raster_size[1]), rasterRGB@extent@ymin + arrow_length)
      P_Y <- c(rasterRGB@extent@xmin + arrow_length, rasterRGB@extent@ymin + 0.01*(raster_size[2]))
      Txt_Adj <- c(0, 0)
    }
    if( rotation == 270 )
    {
      P_0 <- c(rasterRGB@extent@xmax - 0.01*(raster_size[1]), rasterRGB@extent@ymax - 0.01*(raster_size[2]))
      P_X <- c(rasterRGB@extent@xmax - 0.01*(raster_size[1]), rasterRGB@extent@ymax - arrow_length)
      P_Y <- c(rasterRGB@extent@xmax - arrow_length, rasterRGB@extent@ymax - 0.01*(raster_size[2]))
      Txt_Adj <- c(1, 1)
    }
    if( rotation == 180 )
    {
      P_0 <- c(rasterRGB@extent@xmax - 0.01*(raster_size[1]), rasterRGB@extent@ymin + 0.01*(raster_size[2]))
      P_X <- c(rasterRGB@extent@xmax - arrow_length, rasterRGB@extent@ymin + 0.01*(raster_size[2]))
      P_Y <- c(rasterRGB@extent@xmax - 0.01*(raster_size[1]), rasterRGB@extent@ymin + arrow_length)
      Txt_Adj <- c(1, 0)
    }

    arrows(x0 = P_0[1], y0 = P_0[2], x1 = P_X[1], y1 = P_X[2], code = 2, lwd = 2, col = "white", length = 0.1)
    text ( x = P_X[1], y = P_X[2], labels = " X", col = "white", adj = Txt_Adj, cex = 0.8)
    arrows(x0 = P_0[1], y0 = P_0[2], x1 = P_Y[1], y1 = P_Y[2], code = 2, lwd = 2, col = "white", length = 0.1)
    text ( x = P_Y[1], y = P_Y[2], labels = " Y", col = "white", adj = Txt_Adj, cex = 0.8)
  }

  if(!is.null(roi_rectangle) && !zoom)
  {
    #roi_rectangle  is c(left, rigth, top, bottom)
    rect( xleft = roi_rectangle[1], xright = roi_rectangle[2], ytop = roi_rectangle[4], ybottom = roi_rectangle[3], border = "red", lwd = 2 )
  }
}

.plotIntensityScale<-function(img, color = NULL, light=3)
{
  max_int<-max(raster::values(img$raster))

  #Create the raster
  ncols_scale <- 10
  scale_raster <- raster::raster( nrow = 255, ncol = ncols_scale, xmn= 0, xmx= ncols_scale, ymn= 0, ymx= 255)
  raster::values(scale_raster) <- as.vector(matrix(seq(from=max_int, to=0, length.out = 255), nrow = ncols_scale, ncol = 255, byrow = T))

  #Check RGB channels
  if( is.null(color))
  {
    #Remap Color 2 rainbow Space  (24bits color space, 8 bits per channel 255 steps)
    RGB_raster<-.ReMappingIntensity2HSV(scale_raster, value_multiplier = light)
  }
  else
  {
    img_zero<-.InitRGBEmptyRaster( scale_raster@ncols, scale_raster@nrows )
    img_255 <-  .NormalizeTo255(scale_raster, light)
    if(color == "R")
    {
      RGB_raster <- raster::addLayer( img_255, img_zero$raster, img_zero$raster )
    }
    if( color == "G" )
    {
      RGB_raster <- raster::addLayer( img_zero$raster, img_255, img_zero$raster )
    }
    if( color == "B")
    {
      RGB_raster <- raster::addLayer( img_zero$raster, img_zero$raster, img_255 )
    }
  }


  #Setting my tricky par values...
  par( bg = "black", fg =  "white", col.lab="white", xaxt="n", yaxt="n", col.axis = "white", col.main = "white", col.sub = "white",
       cex.axis = 0.6, mar = c(1,0,1,2),  mgp = c(3, 0.5, 0.5))

  raster::plotRGB(RGB_raster, axes = T, asp = 1, interpolate = T  )

  #Add axes
  yAxis<- seq(0, RGB_raster@extent@ymax, length.out = 11)
  yLabels <- sprintf( "%0.1e", seq(0, max_int, length.out = 11))
  par(xaxt = "l", yaxt = "l")
  axis(side=2, tck = -0.015, cex.axis = 0.7, pos = 0, at = yAxis, labels = F, las = 1) #Y left axes
  if( max_int == 0 )
  {
    axis(side=4, tck = -0.015, cex.axis = 0.7, pos = RGB_raster@extent@xmax, at = yAxis, labels = F) #Y right axes
  }
  else
  {
    axis(side=4, tck = -0.015, cex.axis = 0.7, pos = RGB_raster@extent@xmax, at = yAxis, labels = yLabels, las = 1) #Y right axes
  }

  axis(side = 1, tck = -0.015, cex.axis = 0.7, labels = F, pos = 0, at = c(0,ncols_scale))
  axis(side = 3, tck = -0.015, cex.axis = 0.7, labels = F, pos = RGB_raster@extent@ymax, at = c(0,ncols_scale))

  #Add the main title
  if(max_int == 0)
  {
    mtext("Channel Disabled", side = 2, line = -1, cex = 0.8, adj = 0.5  )
  }
  else
  {
    mtext(sprintf("m/z: %0.3f+/-%0.2f Da", img$mass, img$tolerance), side = 2, line = -1, cex = 0.8, adj = 0.5  )
  }
}

#Combinar funcions anteior per fer aixo facil
#La idea es que si mass.peak i tolerance es passen com a vectors de 1 a 3 valors es fan imatges RGB
#' Plot a mass image of a up to 3 selected ions.
#'
#' @param img an rMSI data object.
#' @param mass.peak a vector of up to 3 elements containing the mass of ions to plot.
#' @param tolerance a vector of up to 3 elements containing the mass range to plot for each ion.
#' @param XResLevel the interpolation factor (default is 3).
#' @param NormalizationCoefs optional parameter with a matrix containing the scaling factors for each pixel.
#' @param rotation the rotation of image expressed in deg. Valid values are: 0, 90, 180 and 270.
#' @param show_axes if true axis will be plotted. Otherwise a um scale and xy arrows  are drawn.
#' @param scale_to_global_intensity scale the image intensity to fit into the global intensity (only for RGB).
#' @param vlight the lighting of the plotted image.
#' @param crop_area the region of image to show in a vector formated as c(left, right, bottom, top) if null the whole image is plot.
#' @param intensity_limit limit the all pixels intensities to a given value. If null no limiting is used.
#'
#' Plots a mass image. If only one mass ion are used a rainbow color code image is generated. If more ions are used each mass will
#' be encoded in an RGB color channel.
#'
#' @export
#'
plotMassImageByPeak<-function(img, mass.peak, tolerance=0.25, XResLevel = 3, NormalizationCoefs = NULL, rotation = 0, show_axes = F, scale_to_global_intensity = F, vlight = 3, crop_area = NULL, intensity_limit = NULL)
{
  numberOfChannels <- 1

  if(length(mass.peak) == 1 )
  {
    #Single Ion image
    im_sgn<-.buildImageByPeak(img, mass.peak, tolerance, NormalizationCoefs)

    #Apply limit intensity to raster object
    if(!is.null(intensity_limit) && !is.na(intensity_limit[1]) )
    {
      raster::values(im_sgn$raster)[ raster::values(im_sgn$raster) > intensity_limit ] <- intensity_limit[1]
    }

    if(scale_to_global_intensity)
    {
      if(class(img$mean) == "MassSpectrum")
      {
        #Handling old data format
        raster_RGB<-.BuildSingleIonRGBImage(im_sgn, XResLevel = XResLevel, global_intensity_scaling_factor = max(img$mean@intensity), light = vlight)
      }
      else
      {
        raster_RGB<-.BuildSingleIonRGBImage(im_sgn, XResLevel = XResLevel, global_intensity_scaling_factor = max(img$mean), light = vlight)
      }
    }
    else
    {
      raster_RGB<-.BuildSingleIonRGBImage(im_sgn, XResLevel = XResLevel, light = vlight)
    }
  }
  else
  {
    #Multiple Ions image
    if(length(tolerance) == 1)
    {
      tolerance[2] <- tolerance[1]
    }

    im_R<-.buildImageByPeak(img, mass.peak[1], tolerance[1], NormalizationCoefs)
    im_G<-.buildImageByPeak(img, mass.peak[2], tolerance[2], NormalizationCoefs)

    #Apply limit intensity to raster objects R and G
    if(!is.null(intensity_limit))
    {
      if(!is.na(intensity_limit[1]))
      {
        raster::values(im_R$raster)[ raster::values(im_R$raster) > intensity_limit ] <- intensity_limit[1]
      }
      if(!is.na(intensity_limit[2]))
      {
        raster::values(im_G$raster)[ raster::values(im_G$raster) > intensity_limit ] <- intensity_limit[2]
      }
    }

    numberOfChannels <- 2

    if( length(mass.peak) == 2 )
    {
      #Use Red and Green only
      im_B<-.InitRGBEmptyRaster( img$size["x"], img$size["y"] )
    }
    else
    {
      #Use RGB
      if(length(tolerance) == 2)
      {
        tolerance[3] <- tolerance[1]
      }
      im_B<-.buildImageByPeak(img, mass.peak[3], tolerance[3], NormalizationCoefs)

      #Apply limit intensity to raster object B
      if(!is.null(intensity_limit))
      {
        if(!is.na(intensity_limit[3]))
        {
          raster::values(im_B$raster)[ raster::values(im_B$raster) > intensity_limit ] <- intensity_limit[3]
        }
      }

      numberOfChannels <- 3
    }

    raster_RGB <-.BuildRGBImage( im_R, im_G, im_B,  XResLevel = XResLevel)
  }


  layout( matrix( (numberOfChannels+1):1, ncol = (1+numberOfChannels), nrow = 1, byrow = TRUE ), widths = c(7, rep(1, numberOfChannels)) )

  if(numberOfChannels == 1 )
  {
    .plotIntensityScale(im_sgn, light =  vlight)
  }
  else
  {
    if( numberOfChannels == 3 )
    {
      .plotIntensityScale(im_B, "B" )
    }
    .plotIntensityScale(im_G, "G" )
    .plotIntensityScale(im_R, "R" )
  }

  .plotMassImageRGB(raster_RGB, cal_um2pixels = img$pixel_size_um, rotation = rotation, display_axes = show_axes, roi_rectangle = crop_area, zoom = !is.null(crop_area))
}
