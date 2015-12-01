# Label Free FFT based alignment


####TODO: Improve the average spectra computation to reduce RAM usage in huge datasets

#TOP-METHOD: Align a whole dataset using multi-core capabilities and mean spectra as reference
# data -  a matrix of spectra to align
# iterations - number of iterations over aligned data
# returns - the aligned dataset
MatrixAligment<-function(data, ...)
{
  ref <- LabelFreeCreateRef( apply(data, 2, mean) , smoothing = F)
  data<-LabelFreeAlignDataSet(ref, data, ...)
  return(data)
}

#TOP-METHOD: Align a whole image using multi-core capabilities and mean spectra as reference
# raw - image data in ff custom data format
# returns - Nothing! The data is stored in ff object
FullDataSetAligment<-function(raw, ...)
{
  cat("Aligning Image...\n")
  pt<-proc.time()
  pb<-txtProgressBar(min = 0, max = length(raw$data), style = 3 )
  ref <- LabelFreeCreateRef( raw$mean@intensity , smoothing = F )
  CurrCube <- 1
  setTxtProgressBar(pb, CurrCube - 1)
  while( CurrCube <=  length(raw$data))
  {
    CurrCube <- AlignDataChunk(raw, ref, startingCube = CurrCube, ...)
    setTxtProgressBar(pb, CurrCube - 1)
  }
  gc()
  close(pb)
  pt<-proc.time() - pt
  cat(paste("Alignment time:",round(pt["elapsed"], digits = 1),"seconds\n"))
  cat(paste("Average time per spectrum:",round(pt["elapsed"]/length(raw$pos), digits = 1),"seconds\n" ))

  #Recompute mean spectrum
  cat("Calculating average spectrum...\n")
  pt<-proc.time()

  avgI<-apply(matrix(unlist(lapply(raw$data, function(x){ ff::ffrowapply(colSums(x[i1:i2,,drop=FALSE]), X=x, RETURN = TRUE, CFUN = "csum", FF_RETURN = FALSE) })), nrow = length(raw$data), byrow = T), 2, sum)
  avgI<-avgI/( sum(unlist(lapply(raw$data, nrow))) )

  meanSpc<-createMassSpectrum(mass =  raw$mass, intensity = avgI)
  meanSpc<-smoothIntensity(meanSpc, halfWindowSize=2)
  raw$mean<-meanSpc
  pt<-proc.time() - pt
  cat(paste("Average spectrun time:",round(pt["elapsed"], digits = 1),"seconds\n"))
  gc()
}

#Create Ref data structure to avoid unnecessary computation in ref spectrum
# ref - ref spectrum (for example the average spectrum)
# smoothing - if TRUE, oversamples the data to achive a better alignment but a 3 times slower computations
# returns - a list containing top and bottom windowed fft
LabelFreeCreateRef<-function(ref, smoothing = F)
{
  if(smoothing)
  {
    ref<- SmoothUpper(ref)
  }
  #refBottom <- Conj(FFT( TimeWindowLow(  ZeroPadding( ref, rev = F  ) )  ))
  #refTop    <- Conj(FFT( TimeWindowHigh( ZeroPadding( ref, rev = T  ) )  ))

  refBottom <- Conj(FFT( ZeroPadding( TimeWindowLow(ref) , rev = F )  ))
  refTop    <- Conj(FFT( ZeroPadding( TimeWindowHigh(ref), rev = T )  ))

  return(list( ref_fft_bot = refBottom, ref_fft_top = refTop, smooth =  smoothing))
}


#Aligns a Chunk of a dataset to fit the RAM machine limit. The Chunk consists in an intergar number of Cubes
# rawImg - The image file in custom ff data format
# refSpectrum - The reference to align (return value of LabelFreeCreateRef())
# startingCube - The starting ID of data cube to process
# subDataMemPercent - The maximum RAM mem usage in % units relative to machine installed memory
# returns - The next data cube ID to be processed
AlignDataChunk<-function( rawImg, refSpectrum, startingCube = 1, subDataMemPercent = 10, maxcubes = 5, ... )
{


  freeMem <-1024 * (as.numeric(system( "awk '/MemFree:/ {print $2}' /proc/meminfo", intern = T)) +
                    as.numeric(system( "awk '/Buffers:/ {print $2}' /proc/meminfo", intern = T)) +
                    as.numeric(system( "awk '/Cached:/ {print $2}' /proc/meminfo", intern = T))[1])  #Mult by 1024 because /proc/meminfo returns kB


  #This is only implemented for newest linux kernels!
  #freeMem <- 1024 * as.numeric(system( "awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = T)) #Mult by 1024 because /proc/meminfo returns kB

  #Matrix size = NumCols * NumRows * 4bytes (integer)
  NumCols <-length(rawImg$mass)
  NumRows <- round(round(freeMem * (subDataMemPercent/100)) / NumCols / 4)
  NumCubes <- ceiling( NumRows/ nrow(rawImg$data[[1]]))

  #Limit Num of cubes to the dataset size
  NumCubes <- min(NumCubes,  (length(rawImg$data) - startingCube + 1)) #Limit to remaning cubes
  NumCubes <- min(NumCubes, maxcubes)

  CubesId <- startingCube:(startingCube + NumCubes - 1)
  NumRows <- sum(sapply(CubesId, function(x) { nrow(rawImg$data[[x]]) }))

  #Extract data matrix
  m <- matrix(integer(), nrow = NumRows, ncol = NumCols)
  st_row <- 1
  for( i in CubesId)
  {
    m[st_row:(st_row + nrow(rawImg$data[[i]]) - 1),] <- rawImg$data[[i]][,]
    st_row <- st_row + nrow(rawImg$data[[i]])
  }
  gc()

  #Align the matrix
  m <- LabelFreeAlignDataSet(refSpectrum, m, ...) #TODO improve memory management inside this function
  gc()

  #Store the aligned data to rawImg
  st_row <- 1
  for( i in CubesId)
  {
    rawImg$data[[i]][,]<-m[st_row:(st_row + nrow(rawImg$data[[i]]) - 1),]
    st_row <- st_row + nrow(rawImg$data[[i]])
  }
  rm(m)
  gc()

  #Return the next Cube to process
  return( CubesId[length(CubesId)] + 1)
}



#Top-Level Multi Spectra aligment (threaded implemented)
# x - a spectrum to be aligned
# ref - a data structre created with LabelFreeCreateRef() function
# data - a matrix of spectra to align
# returns -  the aligned spectra to ref
LabelFreeAlignDataSet<-function(ref, data, iterations = 1, multithreading = T, ...)
{
  if( multithreading )
  {
    ###TODO: working here! makeForkCluster is better in ram usage
    ##clus <- makeCluster(parallel::detectCores())
    clus <- makeForkCluster(nnodes = parallel::detectCores())

    clusterExport(clus, "ref", envir = environment()) #Export aligned ref spectra to make it visible to all clusters

    #Export aligning functions to all custers
    clusterExport(clus, "LabelFreeAlign")
    clusterExport(clus, "TimeWindowLow")
    clusterExport(clus, "TimeWindowHigh")
    clusterExport(clus, "ZeroPadding")
    clusterExport(clus, "FourierBestCor")
    clusterExport(clus, "LinearScale")
    clusterExport(clus, "FourierLinearShift")
    clusterExport(clus, "FFT")
    clusterExport(clus, "IFFT")

    #Run cluster
    for(i in 1:iterations)
    {
      cat(paste("Align is running iteration ", i, " of ", iterations, "\n", sep = ""))
      data<-t(parApply(clus, data, 1, function(x) LabelFreeAlign(x, ref, ... )))
      gc()
    }
    stopCluster(clus)
  }
  else
  {
    #Use a single processor approach
    for(i in 1:iterations)
    {
      cat(paste("Align is running iteration ", i, " of ", iterations, "\n", sep = ""))
      data <- t(apply(data, 1, function(x) LabelFreeAlign(x, ref, ... )))
      gc()
    }
  }
  return(data)
}

# x - a spectrum to be aligned
# ref - a data structre created with LabelFreeCreateRef() function
# limit - shift limiting in ppm of spectra num of points
# returns -  the aligned spectra to ref
LabelFreeAlign<- function(x, ref, limit = 100)
{
  #0- Smoothing
  if(ref$smooth)
  {
    x<-SmoothUpper(x)
  }

  #1- Hanning Windowing
  x_Low <- TimeWindowLow(x)
  x_High <- TimeWindowHigh(x)

  #2- Zero-padding 2 improve performance
  x_High <- ZeroPadding(x_High, rev = T)
  x_Low <- ZeroPadding(x_Low, rev = F)

  #Get lags
  lagLow <- FourierBestCor(x_Low, ref$ref_fft_bot)
  lagHigh <- FourierBestCor(x_High, ref$ref_fft_top)

  #Limit lags
  lagMax <- (limit/1e6)*length(x)

  if( abs(lagLow) > lagMax) { lagLow = 0}
  if( abs(lagHigh) > lagMax) { lagHigh = 0}

  ##cat(paste("Lags Limited: max =", lagMax, "low =",lagLow,"high =", lagHigh,"\n"))

  #3- Aligment constansts mz'(n) = K*mz(n) + Sh
  x_ <- ZeroPadding(x, rev = F)

  ###This is the OLD constants only valid if Rl = 0 and Rh = N and shift is computed after scaling
  ##K  <- (length(x_) + lagHigh - lagLow)/length(x_)
  ##Sh <- lagLow

  Rl <- 0 #Currently I only use zero as the first reference
  Rh <- length(x_) * 0.9 #Currently I only use N as the first reference
  K  <- (Rh + lagHigh -Rl - lagLow)/(Rh - Rl) #New implementation supporting Other refs diferents than 0 and N. Extremes are Rl and Rh
  ##Sh <- (Rh*lagLow - Rl*lagHigh)/(Rh - Rl + lagHigh - lagLow) #If scling is performed after shift New implementation supporting Other refs diferents than 0 and N. Extremes are Rl and Rh
  Sh <- (Rh*lagLow - Rl*lagHigh)/(Rh - Rl) #If scaling is performed before shift. New implementation supporting Other refs diferents than 0 and N. Extremes are Rl and Rh

  ##cat(paste("K =",K,"Shift =",Sh,"\n"))

  #4- Apply the aligment
  x_ <- FourierLinerScaleShift( x_, scaling = K, shift = Sh) #New implementation with fft based interpolation
  #x_<- LinearScale(x_, K) ###TODO this is the slowest part, try changing it using a fft based time contraction/expansion
  #x_ <- FourierLinearShift(x_, Sh)

  #5- Remove the padded zeros to fit mass vector
  x_ <- x_[1:length(x)]

  #6- De-Smoothing (fit the rigth data size)
  if(ref$smooth)
  {
    x_ <- approx(x_, n = length(x_)/2)$y
  }

  return (x_)
}

#Time shift on a data vector using FFT
# X - FFT of Data vector to be shifted in X axis direction
# shift - The shift offset
FourierLinearShift<-function(x, shift)
{
  x_<-Mod(IFFT( FFT(x) * exp(-2i*pi*(1/length(x))*(1:length(x)) * shift) ))
  return(x_)
}

#Time compression/expansion
# x - Data vector to be scaled in X axis direction
# scaling - The scaling factor ( <1 compression / >1 expansion)
LinearScale<-function(x , scaling)
{
  x_ <- approx(x, n = length(x)*scaling)$y
  if(length(x_) < length(x))
  {
    x_ <- c(x_, rep(0, length(x) - length(x_)))
  }else
  {
    x_ <- x_[1:length(x)]
  }
  return(x_)
}

#Time scale and shift on a data vector using FFT
# X - FFT of Data vector to be shifted in X axis direction
# scaling - The scaling factor ( <1 compression / >1 expansion)
# shift - The shift offset
FourierLinerScaleShift<-function(x , scaling, shift)
{
  X <- fftw::FFT(x)

  #Shift... Sembla que funciona millor fent el shift despres de esclat,
      #el tema es que el shift introdueix  variacions complexes en espai de fourier, si s'introdueixen aquestes variacions a scale
      #es genere aliasing
  ###X <- X*exp(-2i*pi*(1/length(X))*(1:length(X)) * shift)

  #Scale...
  gap <- abs((scaling*length(X)) - length(X))
  X_ <- X #No scaling
  if( scaling < 1) #Time Contraction
  {
    X_ <- c( X[1: (length(X)/2 - floor(gap/2) )], X[( 1 + length(X)/2 + ceiling(gap/2) ):length(X)])
  }
  if(scaling > 1) #Time Expansion
  {
    X_ <- c( X[1:(length(X)/2)], rep(0, round(gap)), X[(1+length(X)/2):length(X)])
  }
  X_ <- X_ * length(X_)/length(X) #Magnitud adjusted to fit original aplitude

  #Shift...
  X_ <- X_*exp(-2i*pi*(1/length(X_))*(1:length(X_)) * shift)

  #Restore tof space...
  x_<- Mod( fftw::IFFT(X_) )
  if(length(x_) < length(x))
  {
    x_ <- c(x_, rep(0, ( length(x) - length(x_) )))
  }
  x_<-x_[1:length(x)]
  return(x_)
}

#Low-Windowing Hanning based
# x - Mass Spectra to apply window
# spectraSplit - the low part of spectra to keep (the resting points to unit will be removed)
# returns - The windowed data
TimeWindowLow<-function(x, spectraSplit = 0.6)
{
  x<-x[1:(spectraSplit*length(x))]
  hann = rev( 0.5 * (1 - cos( 2*pi*(1:(2*length(x)) ) /  ( 2*length(x)) ) )[1: length(x)])
  return(x*hann)
}

#High-Windowing Hanning based
# x - Mass Spectra to apply window
# spectraSplit - the high part of spectra to keep (the resting points to unit will be removed)
# returns - The windowed data
TimeWindowHigh<-function(x, spectraSplit = 0.6)
{
  x<-x[((1 - spectraSplit ) *length(x) + 1) : length(x)]
  hann = 0.5 * (1 - cos( 2*pi*(1:(2*length(x)) ) /  ( 2*length(x)) ) )[1: length(x)]
  return(x*hann)
}

#FFT Correlation
# x - vector to correlate with ref
# ref - Conj(FFT(ref))
# returns - the lag to shift x to fit ref with best correlation
FourierBestCor<-function(x, ref)
{
  X<-FFT(x)
  #Y<-FFT(ref)
  #comb <- X * Conj(Y)
  comb <- X * ref
  cor <- IFFT(comb)

  lag<-which.max(Mod(cor))
  if(lag > (length(cor)/2))
  {
    lag <- length(cor) - lag + 1
  }
  else
  {
    lag <- -(lag - 1)
  }

  return(lag)
}

#Zero-padding to fit 2^n length to speed-up fft
# x - vector to be padded with zeros
# rev - if is FALSE padd zeros at right side, else to left side
ZeroPadding<-function(x, rev = F)
{
  n <- ceiling( log2(length(x)) )
  padding <- (2^n) - length(x)
  if(rev)
  {
    x <- c(rep(0, padding), x)
  }
  else
  {
    x <- c(x, rep(0, padding))
  }
  return(x)
}

#Smooth a curve and interpolate it to double length
# x - spectrum to smooth
# returns - smoothed data with doubled length
SmoothUpper <- function(x)
{
  x_ <- approx(x, n = 2*length(x))$y
  x_ <- smooth(x_)
  return( x_ )
}
