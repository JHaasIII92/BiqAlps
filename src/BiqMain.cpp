
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


void readModel(std::string strFileName, 
                int & nVar, bool & max_problem,
                double * &Q, Sparse & Qs,
                std::vector<Sparse> & As, double * &a,
                std::vector<Sparse> & Bs, double * &b);

int main(int argc, char * argv[])
{
    srand(time(NULL));

    int iBlasThreads = 1;
    openblas_set_num_threads(iBlasThreads);

    //std::string strDataFileName = std::string strDataFileName = argv[1];
    std::string strDataFileName = "/workspaces/BiqAlps/MLT-BiqCrunch_2.0/problems/max-cut/examples/g05_60.4.bc";
    //std::string strDataFileName = "/workspaces/BiqAlps/MLT-BiqCrunch_2.0/problems/generic/examples/sonetgr17.bc";
    //std::string strDataFileName = "/workspaces/BiqAlps/MLT-BiqCrunch_2.0/problems/generic/examples/randprob_prod.bc";
    //std::string strDataFileName = "/workspaces/BiqAlps/MLT-BiqCrunch_2.0/problems/generic/examples/randprob.bc";
    //std::string strDataFileName = "/workspaces/BiqAlps/MLT-BiqCrunch_2.0/problems/generic/examples/randprob_square.bc";

    // data for the problem
    int nVar = 0;
    bool max_problem = true;
    double *Q;
    Sparse Qs;
    std::vector<Sparse> As;
    double *a;
    std::vector<Sparse> Bs;
    double *b;
    
    readModel(strDataFileName, nVar, max_problem, Q, Qs, As, a, Bs, b);
    //std::printf("size of enq con %d",Bs.size());
    //print_vector(b, Bs.size());
    //PrintSparseMatrix(Bs.at(2));
    //PrintMatrix(Q, nVar+1,nVar+1);
    //exit(1);
    BiqModel model(nVar, max_problem, Q, Qs, As, a, Bs, b);
    #ifdef COIN_HAS_MPI
          AlpsKnowledgeBrokerMPI broker(argc, argv, model);
    #else
         AlpsKnowledgeBrokerSerial broker(argc, argv, model);
    #endif

        broker.registerClass(AlpsKnowledgeTypeModel,new BiqModel());
        broker.registerClass(AlpsKnowledgeTypeSolution,
                                 new BiqSolution(&model));
        broker.registerClass(AlpsKnowledgeTypeNode, new BiqTreeNode(&model));
        broker.search(&model); 


    return 0;

}

void readModel(std::string strFileName, 
                int & nVar, bool & max_problem,
                double * &Q, Sparse & Qs,
                std::vector<Sparse> & As, double * &a,
                std::vector<Sparse> & Bs, double * &b)
{
    int N = 0;
    int nCon = 0; 
    int nEqCon = 0;
    int nInEqCon = 0;
    int nBlocks= 0;
    int pos_a = 0;
    int pos_b = 0;
    int pos_rhs = 0;
    double *rhs;
    
    bool bInIneq = false;
    bool bInGeq = false;

    // data for the matrix read in
    int numMat;
    int numPrevMat;
    int numBlock; 
    int i_index;
    int j_index;
    int int_max_prob;

    double dBcVal;
    
    std::vector<double> tmpMatrix0;
    std::string strLine;



    

    std::ifstream myFile(strFileName);
    std::printf("Reading in: %s \n", strFileName.c_str());
    // Read out each comment for now
    while (getline (myFile, strLine) &&
            (strLine.at(0) == ';'   ||
             strLine.at(0) == '\\'  ||
             strLine.at(0) == '#'   ||
             strLine.at(0) == '*' 
            )
    ) 
    {
        std::printf("%s \n", strLine.c_str());
    }

    // read in problem sense
    std::sscanf(strLine.c_str(), "%d", &int_max_prob);
    max_problem = int_max_prob == 1;
    getline(myFile, strLine);
    // read in number of constraints 
    std::sscanf(strLine.c_str(), "%d", &nCon);
    // number of blocks
    getline(myFile, strLine);
    std::sscanf(strLine.c_str(), "%d", &nBlocks);
    // problem size and number of ineq
    getline(myFile, strLine);
    if(nBlocks == 1)
    {
        std::sscanf(strLine.c_str(), "%d", &N);
        nInEqCon = 0;
    }
    else
    {
        std::sscanf(strLine.c_str(), "%d, -%d", &N, &nInEqCon); 
    }
    nVar = N -1;
    nEqCon = nCon - nInEqCon;
    Q = new double[N*N];
    
    // the next line will have the RHS if any cons
        // add space to stor the data
        Bs.resize(nEqCon);
        As.resize(nInEqCon);
        b = new double[nEqCon];
        a = new double[nInEqCon];
        rhs = new double[nCon];

        // get the data 
        //getline (myFile, strLine);

        // read it into 
        for(int i = 0; i < nCon; ++i)
        {
            //std::fscanf(myFile, "%lf", &(rhs[i]));
            myFile >> rhs[i];
        }
        //myFile >> strLine;

    getline (myFile, strLine);
    // getline (myFile, strLine);
    // next we are going to read in matrix data
    tmpMatrix0.resize(N*N);

    zeroOutMatrix(tmpMatrix0);

    std::sscanf(strLine.c_str(), "%d %d %d %d %lf \n", &numMat, &numBlock, &i_index, &j_index, &dBcVal);
    --i_index;
    --j_index;

        if(numBlock == 2)
        {
            bInIneq = true;
            if(dBcVal < 0.0)
            {
                bInGeq = true;
            }
        }
        else
        {
            if(numMat == 0 && !max_problem)
            {
                dBcVal = -dBcVal;
            }
            // case in quad or
            // i_index < N - 1 && j_index < N - 1
 if(i_index < N-1 && j_index < N-1 && i_index != j_index)
            {
                tmpMatrix0.at(i_index + j_index*N) += 0.25 * dBcVal;
                tmpMatrix0.at(j_index + i_index*N) += 0.25 * dBcVal;
                // add some to the linear part
                tmpMatrix0.at(i_index + N*(N-1)) += 0.25 * dBcVal;
                tmpMatrix0.at(N-1 + N*i_index) += 0.25 * dBcVal;
                tmpMatrix0.at(N-1 + j_index*N) += 0.25 * dBcVal;
                tmpMatrix0.at(j_index + N*(N-1)) += 0.25 * dBcVal;
                // add some to the const
                tmpMatrix0.at(N-1 + N*(N-1)) += 0.5 * dBcVal;
            }
            // case we are on the diag of the quad
            // zero it out and move it to the const
            else if(i_index < N-1 && j_index < N-1 && i_index == j_index)
            {
                tmpMatrix0.at(i_index + j_index*N) = 0; 
                tmpMatrix0.at(N-1 + N*(N-1)) += 0.5 * dBcVal;
                // add some to the linear part
                tmpMatrix0.at(i_index + N*(N-1)) += 0.25 * dBcVal;
                tmpMatrix0.at(N-1 + N*i_index) += 0.25 * dBcVal;
            }
            // case we are linear 
            else if((i_index < N-1 && j_index == N-1) ||
                    (i_index == N-1 && j_index < N-1) )
            {
                tmpMatrix0.at(i_index + j_index*N) += 0.5 * dBcVal;
                tmpMatrix0.at(j_index + i_index*N) += 0.5 * dBcVal;

                // add some to the const
                tmpMatrix0.at(N-1 + N*(N-1)) += dBcVal;
            }
            // case in const 
            // we are in const when i_index == N - 1  && j_index == N - 1
            else
            {
                tmpMatrix0.at(j_index + i_index*N) += dBcVal;
            }
        }

    numPrevMat = numMat;

    // read until end of line
    bool bReading = true;
    
    while(bReading)
    {
        if(getline(myFile, strLine))
        {
            std::sscanf(strLine.c_str(), "%d %d %d %d %lf \n", &numMat, &numBlock, &i_index, &j_index, &dBcVal);
        }
        else
        {
            bReading = false;
        }

        if(numMat != numPrevMat || !bReading)
        {
            if(bInGeq)
            {
                for(auto &it: tmpMatrix0)
                {
                    it = -it;
                }
            }

            //printOutMatrix(tmpMatrix0);
            // now make a sparse matrix
            if(numPrevMat == 0)
            {
                for(int i = 0; i < N*N; ++i)
                {
                    dBcVal = tmpMatrix0.at(i);
                    Q[i] = dBcVal;
                }

                FillSparseMatrix(Qs, tmpMatrix0);
                //PrintSparseMatrix(Qs);
            }
            else if(!bInIneq) 
            {
                
                Sparse sparseCon;
                FillSparseMatrix(sparseCon, tmpMatrix0);
                Bs.at(pos_b) = sparseCon;
                b[pos_b] = rhs[pos_rhs];

                //PrintSparseMatrix(Bs.at(pos_b));

                ++pos_b;
                ++pos_rhs; 
            }
            else 
            {
                Sparse sparseCon;
                FillSparseMatrix(sparseCon, tmpMatrix0);
                As.at(pos_a) = sparseCon;
                if(bInGeq)
                {
                    a[pos_a] = - rhs[pos_rhs];
                }
                else
                {
                    a[pos_a] = rhs[pos_rhs];
                }
                ++pos_a;
                ++pos_rhs; 
            }

            // clear out for next matrix
            if(bReading)
            {
                zeroOutMatrix(tmpMatrix0);
                bInIneq = false;
                bInGeq = false;
            }
            else
            {
                break;
            }
        }

        if(numBlock == 2)
        {
            bInIneq = true;
            if(dBcVal < 0.0)
            {
                bInGeq = true;
            }
        }
        else
        {
            --i_index;
            --j_index;
            if(numMat == 0 && !max_problem)
            {
                dBcVal = -dBcVal;
            }
            // case in quad or
            // i_index < N - 1 && j_index < N - 1
            if(i_index < N-1 && j_index < N-1 && i_index != j_index)
            {
                tmpMatrix0.at(i_index + j_index*N) += 0.25 * dBcVal;
                tmpMatrix0.at(j_index + i_index*N) += 0.25 * dBcVal;
                // add some to the linear part
                tmpMatrix0.at(i_index + N*(N-1)) += 0.25 * dBcVal;
                tmpMatrix0.at(N-1 + N*i_index) += 0.25 * dBcVal;
                tmpMatrix0.at(N-1 + j_index*N) += 0.25 * dBcVal;
                tmpMatrix0.at(j_index + N*(N-1)) += 0.25 * dBcVal;
                // add some to the const
                tmpMatrix0.at(N-1 + N*(N-1)) += 0.5 * dBcVal;
            }
            // case we are on the diag of the quad
            // zero it out and move it to the const
            else if(i_index < N-1 && j_index < N-1 && i_index == j_index)
            {
                tmpMatrix0.at(i_index + j_index*N) = 0; 
                tmpMatrix0.at(N-1 + N*(N-1)) += 0.5 * dBcVal;
                // add some to the linear part
                tmpMatrix0.at(i_index + N*(N-1)) += 0.25 * dBcVal;
                tmpMatrix0.at(N-1 + N*i_index) += 0.25 * dBcVal;
            }
            // case we are linear 
            else if((i_index < N-1 && j_index == N-1) ||
                    (i_index == N-1 && j_index < N-1) )
            {
                tmpMatrix0.at(i_index + j_index*N) += 0.5 * dBcVal;
                tmpMatrix0.at(j_index + i_index*N) += 0.5 * dBcVal;

                // add some to the const
                tmpMatrix0.at(N-1 + N*(N-1)) += dBcVal;
            }
            // case in const 
            // we are in const when i_index == N - 1  && j_index == N - 1
            else
            {
                tmpMatrix0.at(j_index + i_index*N) += dBcVal;
            }
        }
        numPrevMat = numMat;
    }
    myFile.close();
}