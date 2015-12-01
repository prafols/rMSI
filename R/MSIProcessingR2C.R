#MSI Processing model based on eficient C++ procecing relaying on Rcpp package


#Main Processing function, but maybe it should be in a script run directly to use a global object pointing to image????
runMSIProcessing<-function()
{
  ncores<-parallel::detectCores()
  runMSIProcessingCpp( numberOfThreads = ncores, numberOfCubes = length(testData), rLoadCube = LoadCubeWrapper, rSaveCube = SaveCubeWrapper)
}



# Load a MSI ff cube from the image
# iCube: The cube ID to load
# Return: The numeric matrix containing intensities
LoadCubeWrapper<-function( iCube )
{
  return(testData[[iCube]])
}

#Save a MSI ff cube
# iCube: The cube ID to load
# cubeData: A numerix matrix contining the cube to be stored
# Return: Nothing
SaveCubeWrapper<-function( iCube, cubeData )
{
  testData[[iCube]] <<- cubeData #Yes it must point to the global var, otherwise a copy is created locally and no data is modified
}

##################################################### TESTING WITH GENERIC MATRICES #######################
# mz_length <- 10000
# cube_rows <- 100
# ncubes <- 200
# testData <- list()
# for(i in 1:ncubes)
# {
#   mat<-matrix(as.numeric(i), nrow = cube_rows, ncol = mz_length)
#   testData[[length(testData) + 1]] <- mat
# }

