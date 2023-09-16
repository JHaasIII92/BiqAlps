
#include "AlpsConfig.h"

#include <iostream>

#include "CoinError.hpp"
#include "CoinTime.hpp"

#ifdef COIN_HAS_MPI
#  include "AlpsKnowledgeBrokerMPI.h"
#else
#  include "AlpsKnowledgeBrokerSerial.h"
#endif

#include "BiqModel.h"
#include "BiqSolution.h"

int main(int argc, char * argv[])
{
    try{
        BiqModel model;
        
#ifdef COIN_HAS_MPI
        AlpsKnowledgeBrokerMPI broker(argc, argv, model);
#else
        AlpsKnowledgeBrokerSerial broker(argc, argv, model);
#endif

        // 2: Register model, solution, and tree node
        broker.registerClass(AlpsKnowledgeTypeModel, new BiqModel());
        broker.registerClass(AlpsKnowledgeTypeSolution,
                             new BiqSolution(&model));
        broker.registerClass(AlpsKnowledgeTypeNode, new BiqTreeNode(&model));

        // 4: Solve the problem
        broker.search(&model);
    }

    catch(CoinError& er) {
        std::cerr << "ERROR:" << er.message() << std::endl
                  << " from function " << er.methodName() << std::endl
                  << " from class " << er.className() << std::endl;
    }
    catch(...) {
        std::cerr << "Something went wrong!" << std::endl;
    }

    return 0;
}