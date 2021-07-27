/*************************************************************************
 *     rMSIproc - R package for MSI data processing
 *     Copyright (C) 2017 Pere Rafols Soler
 *     Copyright (C) 2021 Lluc Sementé Fernàndez
 * 
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 * 
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **************************************************************************/

#include "fuzzyCmeans.hpp"
using namespace Rcpp;

//************ Constructor with ROIs ************//
fuzzyCmeans::fuzzyCmeans(NumericMatrix r_peakMatrix, NumericVector r_ROIs, parametersFCM *r_parameters)
{
  // Parameters
  parameters = r_parameters;
  
  // ROI flag
  flagROIs = true;
  
  // Peak Matrix
  peakMatrix = new double*[parameters->numPixels];
  for(int i = 0; i < parameters->numPixels; i++)
  {
    peakMatrix[i] = new double[parameters->numPeaks];
    for(int j = 0; j < parameters->numPeaks; j++)
    {
      peakMatrix[i][j] = r_peakMatrix(i,j);
    }
  }
  
  // ROIs
  ROIs = new int[parameters->numPixels];
  for(int i = 0; i < parameters->numPixels; i++)
  {
    ROIs[i] = r_ROIs(i);
  }
  
  // Centroid Matrix
  centroidMatrix = new double*[parameters->numClusters];
  for(int i = 0; i < parameters->numClusters; i++)
  {
    centroidMatrix[i] = new double[parameters->numPeaks];
    for(int j = 0; j < parameters->numPeaks; j++)
    {
      centroidMatrix[i][j] = 0;
    }
  }
  
  // Membership Matrix
  membershipMatrix = new double*[parameters->numPixels];
  for(int i = 0; i < parameters->numPixels; i++)
  {
    membershipMatrix[i] = new double[parameters->numClusters];
    for(int j = 0; j < parameters->numClusters; j++)
    {
      membershipMatrix[i][j] = 0;
    }
  }
  
  // Distance Matrix
  distanceMatrix = new double*[parameters->numPixels];
  for(int i = 0; i < parameters->numPixels; i++)
  {
    distanceMatrix[i] = new double[parameters->numClusters];
    for(int j = 0; j < parameters->numClusters; j++)
    {
      distanceMatrix[i][j] = 0;
    }
  }
  
  // Objective function vector
  vectorJ = new double[parameters->maxIterations];
  for(int i = 0; i < (parameters->maxIterations); i++)
  {
    vectorJ[i] = 0;
  }
  
  //Initialization of Membership Matrix
  initializeMembershipMatrixWithROIs();
  
  //Preparing first Objective matrix computation
  updateCentroidMatrix();
  updateDistanceMatrix();
  computeObjectiveFunction();
  J = newJ;
  vectorJ[0] = newJ;
}



//************ Constructor without ROIs ************//
fuzzyCmeans::fuzzyCmeans(NumericMatrix r_peakMatrix, parametersFCM *r_parameters)
{
  // Parameters
  parameters = r_parameters;
  
  // ROI flag
  flagROIs = false;

  // Peak Matrix
  peakMatrix = new double*[parameters->numPixels];
  for(int i = 0; i < parameters->numPixels; i++)
  {
    peakMatrix[i] = new double[parameters->numPeaks];
    for(int j = 0; j < parameters->numPeaks; j++)
    {
      peakMatrix[i][j] = r_peakMatrix(i, j);
    }
  }
  
  // Centroid Matrix
  centroidMatrix = new double*[parameters->numClusters];
  for(int i = 0; i < parameters->numClusters; i++)
  {
    centroidMatrix[i] = new double[parameters->numPeaks];
    for(int j = 0; j < parameters->numPeaks; j++)
    {
      centroidMatrix[i][j] = 0;
    }
  }
  
  // Membership Matrix
  membershipMatrix = new double*[parameters->numPixels];
  for(int i = 0; i < parameters->numPixels; i++)
  {
    membershipMatrix[i] = new double[parameters->numClusters];
    for(int j = 0; j < parameters->numClusters; j++)
    {
      membershipMatrix[i][j] = 0;
    }
  }
  
  // Distance Matrix
  distanceMatrix = new double*[parameters->numPixels];
  for(int i = 0; i < parameters->numPixels; i++)
  {
    distanceMatrix[i] = new double[parameters->numClusters];
    for(int j = 0; j < parameters->numClusters; j++)
    {
      distanceMatrix[i][j] = 0;
    }
  }
  
  // Objective function vector
  vectorJ = new double[parameters->maxIterations];
  for(int i = 0; i < (parameters->maxIterations); i++)
  {
    vectorJ[i] = 0;
  }
  
  //Initialization of Membership Matrix
  initializeMembershipMatrixRandom();
  
  //Preparing first Objective matrix computation
  updateCentroidMatrix();
  updateDistanceMatrix();
  computeObjectiveFunction();
  J = newJ;
  vectorJ[0] = newJ;
}



//************ Destructor ************//
fuzzyCmeans::~fuzzyCmeans()
{
  // Peak Matrix
  for(int i = 0; i < parameters->numPixels; i++)
  {
    delete[] peakMatrix[i];
  }
  delete[] peakMatrix;
  
  
  // Centroid Matrix
  for(int i = 0; i < parameters->numClusters; i++)
  {
    delete[] centroidMatrix[i];
  }
  delete[] centroidMatrix;
  
  
  // Membership Matrix
  for(int i = 0; i < parameters->numPixels; i++)
  {
    delete[] membershipMatrix[i];
  }
  delete[] membershipMatrix;
  
  // Distance Matrix
  for(int i = 0; i < parameters->numPixels; i++)
  {
    delete[] distanceMatrix[i];
  }
  delete[] distanceMatrix;
  
  // ROIs
  if(flagROIs)
  {
    delete[] ROIs;
  }
  
  // vectorJ
  delete[] vectorJ;
}


//************ Initialize membership matrix with ROIs values ************//
void fuzzyCmeans::initializeMembershipMatrixWithROIs()
{
  for(int i = 0; i < parameters->numPixels; i++)
  {
    membershipMatrix[i][ROIs[i]] = 1;
  }
}


//************ Initialize membership matrix with random values ************//
void fuzzyCmeans::initializeMembershipMatrixRandom()
{
  int choosenCluster = 0;
  for(int i = 0; i < parameters->numPixels; i++)
  {
    choosenCluster = rand() % parameters->numClusters;
    membershipMatrix[i][choosenCluster] = 1;
  }
}


//************ Update membership matrix ************//
void fuzzyCmeans::updateMembershipMatrix()
{
  double totalDistance = 0;
  for(int i = 0; i < parameters->numClusters; i++)
  {
    for(int j = 0; j < parameters->numPixels; j++)
    {
      totalDistance = 0;
      for(int k = 0; k < parameters->numClusters; k++)
      {
        totalDistance += pow(distanceMatrix[j][i]/distanceMatrix[j][k], (2/((parameters->m)-1)));
      }
      membershipMatrix[j][i] = 1/totalDistance;
    }
  }
}



//************ Update centroid matrix ************//
void fuzzyCmeans::updateCentroidMatrix() //TODO: Optimize the loop nesting
{
  double numerator = 0;
  double denominator = 0;
  for(int i = 0; i < parameters->numClusters; i++)
  {
    denominator = 0;
    for(int k = 0; k < parameters->numPixels; k++)
    {
      denominator += pow(membershipMatrix[k][i], parameters->m);
    }
    
    for(int j = 0; j < parameters->numPeaks; j++)
    {
      numerator = 0;
      for(int k = 0; k < parameters->numPixels; k++)
      {
        numerator += pow(membershipMatrix[k][i], parameters->m)*peakMatrix[k][j];
      }
      centroidMatrix[i][j] = numerator/denominator;
    }
  }
}



//************ Update centroid matrix ************//
void fuzzyCmeans::updateDistanceMatrix() 
{
  // Distance Matrix
  distanceMatrix = new double*[parameters->numPixels];
  for(int i = 0; i < parameters->numPixels; i++)
  {
    distanceMatrix[i] = new double[parameters->numClusters];
    for(int j = 0; j < parameters->numClusters; j++)
    {
      distanceMatrix[i][j] = distancePixelToCluster(i, j);
    }
  }
}



//************ Distance between pixel and cluster center ************//
double fuzzyCmeans::distancePixelToCluster(int imagePixelID, int clusterCentroidID)
{
  double distance = 0;
  for(int j = 0; j < parameters->numPeaks; j++)
  {
    distance += pow((peakMatrix[imagePixelID][j]-centroidMatrix[clusterCentroidID][j]), 2);
  }
  distance = sqrt(distance);
  return distance;
}



//************ Update centroid matrix ************//
void fuzzyCmeans::computeObjectiveFunction()
{
  newJ = 0;
  for(int i = 0; i < parameters->numClusters; i++)
  {
    for(int j = 0; j < parameters->numPixels; j++)
    {
      newJ +=  pow(distanceMatrix[j][i], 2)*pow(membershipMatrix[j][i], parameters->m); 
    }
  }
}


  
//************ Executioning function ************//
List fuzzyCmeans::run()
{
  for(int iteration = 0; iteration < parameters->maxIterations; iteration++)
  {
    if(parameters->verbose)
    {
      Rcout << "Iteration: " << iteration << " | ";
      Rcout << "Update centroid matrix ... ";
    }
    updateCentroidMatrix();
    
    if(parameters->verbose)
    {
      Rcout << "OK | Update distance matrix ... ";
    }
    updateDistanceMatrix();
    
    if(parameters->verbose)
    {
      Rcout << "OK | Update membership matrix ... ";
    }
    updateMembershipMatrix();
    
    if(parameters->verbose)
    {
      Rcout << "OK | Objective function ... OK | ";
    }
    computeObjectiveFunction();
    
    if(parameters->verbose)
    {
     Rprintf(" OFD : %f \n", double(J-newJ));
    }
    vectorJ[iteration] = newJ;
    if((double(J-newJ) <= parameters->epsilon)) //End condition is epsilon and more than 4 iterations to prevent early false convergence
    {
      if(parameters->verbose)
      {
        Rcout << "End condition reached: OFD less than epsilon \n";
      }
      NumericMatrix cM(parameters->numClusters, parameters->numPeaks);
      NumericMatrix mM(parameters->numPixels, parameters->numClusters);
      NumericVector vecJ(iteration+1);

      for(int i = 0; i < parameters->numClusters; i++)
      {
        for(int j = 0; j < parameters->numPeaks; j++)
        {
          cM(i,j) = centroidMatrix[i][j];
        }
        
        for(int j = 0; j < parameters->numPixels; j++)
        {
          mM(j,i) = membershipMatrix[j][i];
        }
      }
      
      for(int i = 0; i <= iteration; i++)
      {
        vecJ(i) = vectorJ[i];
      }
      
      List result(3);
      result(0) = cM;
      result(1) = mM;
      result(2) = vecJ;
      return result;
    }
    J = newJ;
  }
  
  //End condition is iteration 
  if(parameters->verbose)
  {
    Rcout << "End condition reached: Maximum nomber of iteration reached \n";
  }
  NumericMatrix cM(parameters->numClusters, parameters->numPeaks);
  NumericMatrix mM(parameters->numPixels, parameters->numClusters);
  NumericVector vecJ(parameters->maxIterations);
  
  for(int i = 0; i < parameters->numClusters; i++)
  {
    for(int j = 0; j < parameters->numPeaks; j++)
    {
      cM(i,j) = centroidMatrix[i][j];
    }
    
    for(int j = 0; j < parameters->numPixels; j++)
    {
      mM(j,i) = membershipMatrix[j][i];
    }

  }
  for(int i = 0; i < parameters->maxIterations; i++)
  {
    vecJ(i) = vectorJ[i];
  }
  
  List result(3);
  result(0) = cM;
  result(1) = mM;
  result(2) = vecJ;
  return result;

}



// [[Rcpp::export]]
Rcpp::List C_fuzzyCmeansROIs(NumericMatrix peakMatrix, NumericVector ROIs, int numPeaks, int numPixels,
                             int numClusters, int m, int maxIterations, double epsilon, bool verbose) 
{
  fuzzyCmeans::parametersFCM myParameters;
  myParameters.numPeaks = numPeaks;
  myParameters.numPixels = numPixels;
  myParameters.numClusters = numClusters;
  myParameters.m = m;
  myParameters.maxIterations = maxIterations;
  myParameters.epsilon = epsilon;
  myParameters.verbose = verbose;
  fuzzyCmeans myFuzzyCmeans(peakMatrix, ROIs, &myParameters);
  List result = myFuzzyCmeans.run();
  
  return result;
}

// [[Rcpp::export]]
Rcpp::List C_fuzzyCmeansRandom(NumericMatrix peakMatrix, int numPeaks, int numPixels,
                               int numClusters, int m, int maxIterations, double epsilon, bool verbose) 
{
  fuzzyCmeans::parametersFCM myParameters;
  myParameters.numPeaks = numPeaks;
  myParameters.numPixels = numPixels;
  myParameters.numClusters = numClusters;
  myParameters.m = m;
  myParameters.maxIterations = maxIterations;
  myParameters.epsilon = epsilon;
  myParameters.verbose = verbose;
  fuzzyCmeans myFuzzyCmeans(peakMatrix, &myParameters);
  List result = myFuzzyCmeans.run();
  return result;
}


