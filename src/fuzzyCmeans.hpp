#include <Rcpp.h>
#include <stdlib.h>
#include <cmath>
using namespace Rcpp;

class fuzzyCmeans
{
public:
  //Types
  typedef struct{
  int numPeaks;
  int numPixels;
  int numClusters;
  int m;
  int maxIterations;
  double epsilon;
  bool verbose;
  }parametersFCM;
  
  // Functions
  fuzzyCmeans(NumericMatrix peakMatrix, 
              NumericVector ROIs,
              parametersFCM *parameters);  //Constructor with ROIs
  
  fuzzyCmeans(NumericMatrix peakMatrix, 
              parametersFCM *parameters);  //Constructor without ROIS
  
  ~fuzzyCmeans();      //Destructor
  
  List run();           //Algorithm starting function
  
private:
  parametersFCM *parameters;
  int *ROIs;
  double **membershipMatrix;
  double **centroidMatrix;
  double **distanceMatrix;
  double **peakMatrix;
  double* vectorJ;
  double J = 0;   //Objective function value at each iteration
  double newJ = 0;
  bool flagROIs = false;
  
  // Functions
  void updateMembershipMatrix();
  void updateCentroidMatrix();
  void updateDistanceMatrix();
  void initializeMembershipMatrixWithROIs();
  void initializeMembershipMatrixRandom();
  void computeObjectiveFunction();
  double distancePixelToCluster(int imagePixelID, int clusterCentroidID);
};





