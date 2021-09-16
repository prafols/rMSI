#' fuzzyCmeans
#'
#' @param intensityMatrix Intensity matrix where rows represent pixels and columns represent peaks.
#' @param numberOfClusters Number of clusters you want the data to be clusterized into. If ROIs is used, this parameter is effectless. 
#' @param ROIs Vector with a code value for each pixel in the intensity matrix.
#' @param m Fuzzifier. When m = 1, the algorithm behaves like a normal k-means, as m increases, the partitions will be more and more diffuse. Recommended value is 2. 
#' @param maxIterations Maximum number of iterations allowed during the execution.
#' @param epsilon Minimum step required between iterations to stop the algorithm. Values shoul be less than 1 and bigger than 0.                     
#'
#' @return The membership matrix, the centroids matrix and the hard partitioning clustering details.
#' @export
#'
#' @examples
fuzzyCmeans <- function(intensityMatrix, numberOfClusters = 2, ROIs = NULL, m = 2 , maxIterations = 100, epsilon = 1e-4, verbose = FALSE)
{
  if(m < 1)
  {
    writeLines("The parameter 'm' can't be less than 1, using m = 1 ...")
    m <- 1
  }
  
  if(!is.null(ROIs))
  {
    numericROIs <- rep(0, times = length(ROIs))
    cnt <- 0
    for(i in unique(ROIs))
    {
      numericROIs[which(ROIs==i)] <- cnt
      cnt <- cnt + 1
    }
    
    result <- C_fuzzyCmeansROIs(intensityMatrix, numericROIs, ncol(intensityMatrix), nrow(intensityMatrix),
                                      length(unique(numericROIs)), m, maxIterations, epsilon, verbose)
  }else
  {
    result <- C_fuzzyCmeansRandom(intensityMatrix, ncol(intensityMatrix), nrow(intensityMatrix),
                                  numberOfClusters, m, maxIterations, epsilon, verbose)
  }
  hardClustering <- apply(result[[2]],1,which.max)
  hardMembership <- rep(0, times = nrow(result[[2]]))
  for(i in 1:nrow(result[[2]]))
  {
    hardMembership[i] <- result[[2]][i, hardClustering[i]]  
  }
  result[[4]] <- data.frame(cluster = hardClustering, membership = hardMembership)
  result[[5]] <- apply(result[[2]],1,function(x, k) prod(x)*k^k, k = numberOfClusters)
  names(result) <- c("centroidMatrix", "membershipMatrix","objectiveFunction", "clustering", "heterogeneity")
  return(result)
}