
#include "AlpsConfig.h"

#include <iostream>

#include "CoinError.hpp"
#include "CoinTime.hpp"

#ifdef COIN_HAS_MPI
#  include "AlpsKnowledgeBrokerMPI.h"
#else
#  include "AlpsKnowledgeBrokerSerial.h"
#endif

#include "headers/BiqModel.h"
#include "headers/BiqSolution.h"
#include "headers/BiqUtil.h"

int main(int argc, char * argv[])
{
    srand(time(NULL));

    int nVar = 6;
    int nBits = 0;
    
    double Q[49] = {0.0,  -0.25,  -0.25,   0.0,  -0.25,   0.0,   0.0, 
                   -0.25,  0.0,    0.0,    0.0,  -0.25,  -0.25,  0.0, 
                   -0.25,  0.0,    0.0,   -0.25, -0.25,  -0.25,  0.0, 
                    0.0,   0.0,   -0.25,   0.0,   0.0,   -0.25,  0.25,
                   -0.25, -0.25,  -0.25,   0.0,   0.0,   -0.25,  0.0, 
                    0.0,  -0.25,  -0.25,  -0.25, -0.25,   0.0,   0.25,
                    0.0,   0.0,    0.0,    0.25,  0.0,    0.25,  6.0};
   

    std::vector<Sparse> As;
    double *a = NULL;
    std::vector<Sparse> Bs;
    double *b = NULL;
   /*
   Sparse Bs_1 = { 
                     //some quad terms 
                    BiqSparseTriple(0, 0, 1),
                    BiqSparseTriple(1, 0, 2),
                    BiqSparseTriple(2, 0, 3),
                    BiqSparseTriple(3, 0, 4),  
                    BiqSparseTriple(2, 2, 1), 
                     // some linear terms 
                    BiqSparseTriple(nVar, 0, 1) 
                   };
    std::vector<Sparse> Bs = {Bs_1};
    double b[1] = {3.14};
    */    
    //double *pQ = &Q

    // we are only constructing a model object for now
    BiqModel model(nVar, Q, 1, As, a, Bs, b);
    //model.CreateSubProblem();
    return 0;
}