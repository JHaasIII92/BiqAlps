
#include "AlpsConfig.h"

#include <iostream>
#include <string>
#include <fstream>

#include "CoinError.hpp"
#include "CoinTime.hpp"

#ifdef COIN_HAS_MPI
#  include "AlpsKnowledgeBrokerMPI.h"
#else
#  include "AlpsKnowledgeBrokerSerial.h"
#endif

#include "headers/BiqModel.h"

int main(int argc, char * argv[])
{
    srand(16);

    int iBlasThreads = 1;
    openblas_set_num_threads(iBlasThreads);

    BiqModel model;
    #ifdef COIN_HAS_MPI
          AlpsKnowledgeBrokerMPI broker(argc, argv, model);
    #else
         AlpsKnowledgeBrokerSerial broker(argc, argv,  model);
    #endif

    broker.registerClass(AlpsKnowledgeTypeModel,new BiqModel());
    broker.registerClass(AlpsKnowledgeTypeSolution, new BiqSolution(&model));
    broker.registerClass(AlpsKnowledgeTypeNode, new BiqTreeNode(&model));
    const bool bSolutionProvided = model.BiqPar()->entry(BiqParams::bSolutionProvided);
    const double dSolutionValue = model.BiqPar()->entry(BiqParams::dSolutionValue);
    if (bSolutionProvided) {
        // if user provided a solution then add it to the knowledge
        // need to figure out how to recive a solution vector
    }

    broker.search(&model); 

    return 0;

}
