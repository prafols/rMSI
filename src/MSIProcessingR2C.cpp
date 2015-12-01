#include <Rcpp.h>
#include <boost/thread.hpp>

using namespace Rcpp;

struct pdata_t
{
  Rcpp::NumericMatrix *cubes; //Pointer to matrices where each processing data cube will be sotred for each thread
  boost::mutex mtx;
  int *iCube; //This vector will porvide which cube ID is processing each thread
  bool *bDataReady; //This vector will contain true when a worker thread completes a datacube processing
  
  //Condition variable to notify thread ends
  boost::condition_variable  life_end_cond;
  boost::mutex               life_end_mutex;
  bool                       life_end;
};


/*
 * Example of summing a matrix items 1
 */
void processingFunction( Rcpp::NumericMatrix* data  )
{
  for( int irow = 0; irow < data->nrow(); irow++)
  {
    for( int icol = 0; icol < data->ncol(); icol++ )
    {
      (*data)(irow, icol)++; //Simple example of summing one to all matrix 
    }
  }
  
}

/*
 * Threading processing function, any new user processing must be called by this method
 */
void processingThread( int threadSlot, pdata_t *data )
{
  //Call the processing function for this thread
  processingFunction( &(data->cubes[threadSlot]));

  //Save data to R Session and Store the new state of total processed cubes
  data->mtx.lock();
  data->bDataReady[threadSlot] = true;
  data->mtx.unlock();

  //Notify Main thread
  boost::unique_lock<boost::mutex> lock(data->life_end_mutex);
  if(!data->life_end) //Notify only if it haven been done already by another thread
  {
    data->life_end = true;
    data->life_end_cond.notify_all();
  }
}

//wait for end of thread must be implemented in a function to be sure that the lock will be auto-released at function end (scope lock)
void waitForSomeThreadEnd(pdata_t *data)
{
  boost::unique_lock<boost::mutex> lock(data->life_end_mutex);
  while (!data->life_end)
  {
    data->life_end_cond.wait(lock);
  }
  data->life_end =false; //Reset it for the next thread
}

/*
 * Exported function to R session
 * Parameters:
 *  - numberOfThreads: The maximum number of thread that will be created to process the dataset (it's a good idea to keep it as number of cores of cpu )
 *  - numberOfCubes: The total number of cubes to proced
 *  - rLoadCube: The R function used to load a cube from the ff object
 *  - rSaveCube: The R function used to save a cube from the ff object
 */

// [[Rcpp::export]]
void runMSIProcessingCpp( int numberOfThreads, int numberOfCubes, Rcpp::Function rLoadCube, Rcpp::Function rSaveCube )
{
  //Initialize processing data cube
  pdata_t dc;
  dc.life_end = false;
  dc.cubes = new Rcpp::NumericMatrix[numberOfThreads];
  dc.iCube = new int[numberOfThreads];
  dc.bDataReady = new bool[numberOfThreads];
  
  for( int i = 0; i < numberOfThreads; i++)
  {
    dc.iCube[i] = 0;
    dc.bDataReady[i] = false;
  }

  boost::thread *tworkers = new boost::thread[numberOfThreads];
  int iCube = 0;
  bool end_of_program = false;
  while( !end_of_program ) 
  {
    //Load data and start working threads
    for(int iThread = 0; iThread < numberOfThreads; iThread++)
    {
      if(dc.iCube[iThread] == 0 && iCube < numberOfCubes) //No cube assigned then no thread running in this slot
      {
        iCube++;
        dc.iCube[iThread] = iCube;
        dc.cubes[iThread] = rLoadCube( iCube ); //R indexes cubes starting with 1 not zero
        tworkers[iThread] = boost::thread(processingThread, iThread, &dc);
      }
    }
    
    //Wait for thread ends
    waitForSomeThreadEnd(&dc);

    //Save data to RSession and free thread slots
    dc.mtx.lock(); //Any thread locked to this mutex is actually waiting to save data
    for(int iThread = 0; iThread < numberOfThreads; iThread++)
    {
      if( dc.bDataReady[iThread] )
      {
        rSaveCube( dc.iCube[iThread], dc.cubes[iThread] );
        dc.iCube[iThread] = 0; //Mark thread as stopped
        dc.bDataReady[iThread] = false; //Reset data ready state;
        tworkers[iThread].join();
      }
    }
    
    //Check end condition
    if( iCube == numberOfCubes )
    {
      end_of_program = true;
      for(int iThread = 0; iThread < numberOfThreads; iThread++)
      {
        end_of_program &= (dc.iCube[iThread] == 0);
      }
    }
    dc.mtx.unlock();
  }
  
  delete[] dc.cubes;
  delete[] dc.iCube;
  delete[] dc.bDataReady;
  delete[] tworkers;
}
