
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
          std::printf("MPI is enabled\n");
    #else
         AlpsKnowledgeBrokerSerial broker(argc, argv,  model);
    #endif

        broker.registerClass(AlpsKnowledgeTypeModel,new BiqModel());
        //std::printf("new BiqModel registered\n");
        broker.registerClass(AlpsKnowledgeTypeSolution,
                                 new BiqSolution(&model));
        //std::printf("new BiqSolution registered\n");
        broker.registerClass(AlpsKnowledgeTypeNode, new BiqTreeNode(&model));
        //std::printf("*************************** Begin search *************************** \n");
        broker.search(&model); 
    return 0;

}
