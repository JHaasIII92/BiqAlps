#include "headers/BiqModel.h"
#include "AlpsKnowledgeBroker.h"
#include "AlpsKnowledge.h"
#include "headers/BiqSolution.h"
#include "headers/BiqParams.h"
#include "headers/BiqMessage.h"
#include <cstring>
#include <ctime>

void BiqModel::InitModel() 
{
    // get the params needed for initilization 
    const bool bProdCons = BiqPar_->entry(BiqParams::bAddProductConstraints);
    const bool bSolutionProvided = BiqPar_->entry(BiqParams::bSolutionProvided);
    const int nCuts = BiqPar_->entry(BiqParams::nCuts);
    const int MaxNineqAdded = BiqPar_->entry(BiqParams::MaxNineqAdded);

    ISUPPZ_ = new int[2*N_];
    X_ = new double[N_*N_];
    W_ = new double[N_];
    Z_ = new double[N_*N_];

    
    /* DSYEVR optimal workspace query */
    M_ = 0;
    nDWORK_ = 26 * N_;
    nIWORK_ =  10 * N_;
    DWORK_ = new double[nDWORK_];
    IWORK_ = new int[nIWORK_]; 
    
    
    mA_ = As_.size();

    mB_original_ = Bs_.size();
    int numEqualityConstraints = mB_original_ + N_;
    if(bProdCons) numEqualityConstraints += iEqualityIsLinear_.size()*nVar_;

    a_sub_ = new double[mA_];
    b_sub_ = new double[numEqualityConstraints];
    Q_sub_ = new double[N_*N_];
    vOffset_.resize(nVar_, 0);
    b_   = new double[numEqualityConstraints];
    std::copy(b_original_, b_original_ + mB_original_, b_);
    /* L-BFGS-B Data */
    int nMax = numEqualityConstraints + mA_ + MaxNineqAdded;
    int wa_length = (2 * mmax + 5) * nMax + 11 * mmax * mmax + 8 * mmax;
    int iwa_length = 3 * nMax;

    f_ = 0.0;
    g_ = new double[nMax];
    y_ = new double[nMax];
    binf_ = new double[nMax];
    bsup_ = new double[nMax];
    nbd_ = new int[nMax];

    wa_ = new double[wa_length];
    iwa_ = new int[iwa_length];

    
    
    AddDiagCons();
    if(bProdCons) AddProdCons();
    AllocSubCons();
    SetConSparseSize();
    // Set cut data 
    Cuts_.resize(MaxNineqAdded); // or start small and add space when needed
    container.reserve(nCuts);
    nIneq_ = 0;

    // we do not have past bests
    bHasPastBest_ = false;

    // set some memory for the sTemp_sub_ matrix
    sTmp_sub_.resize(N_*N_);
    pdTmp_sub_ = new double[N_*N_];
    pdTmpLinear_ = new double[nVar_];

    viSolution_1_.resize(N_);
    viSolution_2_.resize(N_);

    vdFracSol_.resize(nVar_);

    vdGainsZero_.resize(nVar_);
    vdGainsOne_.resize(nVar_);
    vbGreedyIndex_.resize(nVar_);

#ifdef DEBUG
    testEncodeDecode();
#endif
}


/** Read in Alps and Biq parameters. */
void BiqModel::readParameters(const int argnum, const char * const * arglist) 
{
        AlpsPar_->readFromArglist(argnum, arglist);
        BiqPar_->readFromArglist(argnum, arglist);
        BiqPar_->print();
}

void BiqModel::readInstance(const char* strDataFile)
{
    N_ = 0;
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
    bool bLinear = true;

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



    

    std::ifstream dataFile(strDataFile);
    std::printf("Reading in: %s \n", strDataFile);
    // Read out each comment for now
    while (getline (dataFile, strLine) &&
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
    max_problem_ = int_max_prob == 1;
    getline(dataFile, strLine);
    // read in number of constraints 
    std::sscanf(strLine.c_str(), "%d", &nCon);
    // number of blocks
    getline(dataFile, strLine);
    std::sscanf(strLine.c_str(), "%d", &nBlocks);
    // problem size and number of ineq
    getline(dataFile, strLine);
    if(nBlocks == 1)
    {
        std::sscanf(strLine.c_str(), "%d", &N_);
        nInEqCon = 0;
    }
    else
    {
        std::sscanf(strLine.c_str(), "%d, -%d", &N_, &nInEqCon); 
    }
    nVar_ = N_ - 1;
    nEqCon = nCon - nInEqCon;
    Q_ = new double[N_*N_];
    
    // the next line will have the RHS if any cons
    // add space to stor the data
    Bs_.resize(nEqCon);
    As_.resize(nInEqCon);
    b_original_ = new double[nEqCon];
    a_ = new double[nInEqCon];
    rhs = new double[nCon];
    
    // read all data to rhs split into a_ b_original_ later
    for(int i = 0; i < nCon; ++i)
    {
        dataFile >> rhs[i];           
    }
    if(nCon > 0) getline (dataFile, strLine);
    getline (dataFile, strLine);
    // getline (dataFile, strLine);
    // next we are going to read in matrix data
    tmpMatrix0.resize(N_*N_);

    zeroOutMatrix(tmpMatrix0);

    std::sscanf(strLine.c_str(), "%d %d %d %d %lf \n", &numMat, &numBlock, &i_index, &j_index, &dBcVal);
    --i_index;
    --j_index;
    //std::printf("%d %d %d %d %lf \n", numMat, numBlock, i_index, j_index, dBcVal); 
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
            if(numMat == 0 && !max_problem_)
            {
                dBcVal = -dBcVal;
            }

            if (i_index < N_-1 && j_index < N_-1)
            {
                bLinear = false;
            }
            
            // case in quad or
            // i_index < N - 1 && j_index < N - 1
            if(i_index < N_-1 && j_index < N_-1 && i_index != j_index)
            {
                tmpMatrix0.at(i_index + j_index*N_) += 0.25 * dBcVal;
                tmpMatrix0.at(j_index + i_index*N_) += 0.25 * dBcVal;
                // add some to the linear part
                tmpMatrix0.at(i_index + N_*(N_-1)) += 0.25 * dBcVal;
                tmpMatrix0.at(N_-1 + N_*i_index) += 0.25 * dBcVal;
                tmpMatrix0.at(N_-1 + j_index*N_) += 0.25 * dBcVal;
                tmpMatrix0.at(j_index + N_*(N_-1)) += 0.25 * dBcVal;
                // add some to the const
                tmpMatrix0.at(N_-1 + N_*(N_-1)) += 0.5 * dBcVal;
            }
            // case we are on the diag of the quad
            // zero it out and move it to the const
            else if(i_index < N_-1 && j_index < N_-1 && i_index == j_index)
            {
                tmpMatrix0.at(i_index + j_index*N_) = 0; 
                tmpMatrix0.at(N_-1 + N_*(N_-1)) += 0.5 * dBcVal;
                // add some to the linear part
                tmpMatrix0.at(i_index + N_*(N_-1)) += 0.25 * dBcVal;
                tmpMatrix0.at(N_-1 + N_*i_index) += 0.25 * dBcVal;
            }
            // case we are linear 
            else if((i_index < N_-1 && j_index == N_-1) ||
                    (i_index == N_-1 && j_index < N_-1) )
            {
                tmpMatrix0.at(i_index + j_index*N_) += 0.5 * dBcVal;
                tmpMatrix0.at(j_index + i_index*N_) += 0.5 * dBcVal;

                // add some to the const
                tmpMatrix0.at(N_-1 + N_*(N_-1)) += dBcVal;
            }
            // case in const 
            // we are in const when i_index == N_ - 1  && j_index == N_ - 1
            else
            {
                tmpMatrix0.at(j_index + i_index*N_) += dBcVal;
            }
        }

    numPrevMat = numMat;

    // read until end of line
    bool bReading = true;
    
    while(bReading)
    {
        if(getline(dataFile, strLine))
        {
            std::sscanf(strLine.c_str(), "%d %d %d %d %lf \n", &numMat, &numBlock, &i_index, &j_index, &dBcVal);
            //std::printf("%d %d %d %d %lf \n", numMat, numBlock, i_index, j_index, dBcVal);
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
            //PrintMatrix(tmpMatrix0);
            //std::printf("\n\n");
            // now make a sparse matrix
            if(numPrevMat == 0)
            {
                for(int i = 0; i < N_*N_; ++i)
                {
                    Q_[i] = tmpMatrix0.at(i);
                }

                FillSparseMatrix(Qs_, tmpMatrix0);
            }
            else if(!bInIneq) 
            {
                Sparse sparseCon;
                FillSparseMatrix(sparseCon, tmpMatrix0);

                // If linear, add product constraints
                if(bLinear)
                {
                    iEqualityIsLinear_.push_back(pos_b);
                }

                Bs_.at(pos_b) = sparseCon;
                b_original_[pos_b] = rhs[pos_rhs];
                
                ++pos_b;
                ++pos_rhs; 
            }
            else 
            {
                Sparse sparseCon;
                FillSparseMatrix(sparseCon, tmpMatrix0);
                As_.at(pos_a) = sparseCon;
                if(bInGeq)
                {
                    a_[pos_a] = - rhs[pos_rhs];
                }
                else
                {
                    a_[pos_a] = rhs[pos_rhs];
                }
                if(bLinear)
                {
                    //std::printf("pre push back\n");
                    iInequalityIsLinear_.push_back(pos_b);
                    //std::printf("post push back\n");
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
                bLinear = true;
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
            //std::printf("%d %d %d %d %lf \n", numMat, numBlock, i_index, j_index, dBcVal); 
                        if(numMat == 0 && !max_problem_)
            {
                dBcVal = -dBcVal;
            }

            if (i_index < N_-1 && j_index < N_-1)
            {
                bLinear = false;
            }

            // case in quad or
            // i_index < N_ - 1 && j_index < N_ - 1
            if(i_index < N_-1 && j_index < N_-1 && i_index != j_index)
            {
                tmpMatrix0.at(i_index + j_index*N_) += 0.25 * dBcVal;
                tmpMatrix0.at(j_index + i_index*N_) += 0.25 * dBcVal;
                // add some to the linear part
                tmpMatrix0.at(i_index + N_*(N_-1)) += 0.25 * dBcVal;
                tmpMatrix0.at(N_-1 + N_*i_index) += 0.25 * dBcVal;
                tmpMatrix0.at(N_-1 + j_index*N_) += 0.25 * dBcVal;
                tmpMatrix0.at(j_index + N_*(N_-1)) += 0.25 * dBcVal;
                // add some to the const
                tmpMatrix0.at(N_-1 + N_*(N_-1)) += 0.5 * dBcVal;
            }
            // case we are on the diag of the quad
            // zero it out and move it to the const
            else if(i_index < N_-1 && j_index < N_-1 && i_index == j_index)
            {
                tmpMatrix0.at(i_index + j_index*N_) = 0; 
                tmpMatrix0.at(N_-1 + N_*(N_-1)) += 0.5 * dBcVal;
                // add some to the linear part
                tmpMatrix0.at(i_index + N_*(N_-1)) += 0.25 * dBcVal;
                tmpMatrix0.at(N_-1 + N_*i_index) += 0.25 * dBcVal;
            }
            // case we are linear 
            else if((i_index < N_-1 && j_index == N_-1) ||
                    (i_index == N_-1 && j_index < N_-1) )
            {
                tmpMatrix0.at(i_index + j_index*N_) += 0.5 * dBcVal;
                tmpMatrix0.at(j_index + i_index*N_) += 0.5 * dBcVal;

                // add some to the const
                tmpMatrix0.at(N_-1 + N_*(N_-1)) += dBcVal;
            }
            // case in const 
            // we are in const when i_index == N_ - 1  && j_index == N_ - 1
            else
            {
                tmpMatrix0.at(j_index + i_index*N_) += dBcVal;
            }
        }
        numPrevMat = numMat;
    }
    dataFile.close();
    InitModel();
    return;
}


void BiqModel::InitEmptyModel()
{
    //std::printf("InitEmptyModel\n");
    BiqPar_ = new BiqParams;
    handler_ = new CoinMessageHandler();
    handler_->setLogLevel(2);
    defaultHandler_ = true;
    messages_ = BiqMessage();
    ISUPPZ_ = nullptr;
    X_ = nullptr;
    W_ = nullptr;
    Z_ = nullptr;
    DWORK_ = nullptr;
    IWORK_ = nullptr;
    Q_ = nullptr;
    a_ = nullptr;
    b_ = nullptr;
    b_sub_ = nullptr;
    a_sub_ = nullptr;
    Q_sub_ = nullptr;
    g_ = nullptr;
    y_ = nullptr; 
    binf_ = nullptr;
    bsup_ = nullptr; 
    nbd_ = nullptr; 
    wa_ = nullptr; 
    iwa_ = nullptr;
    pdTmpLinear_ = nullptr;
    pdTmp_sub_ = nullptr;
}

BiqModel::~BiqModel()
{
    FREE_DATA(ISUPPZ_); 
    FREE_DATA(X_);
    FREE_DATA(W_);
    FREE_DATA(Z_);
    FREE_DATA(DWORK_);
    FREE_DATA(IWORK_);
    FREE_DATA(b_sub_);
    FREE_DATA(a_sub_);
    FREE_DATA(Q_sub_);
    FREE_DATA(g_); 
    FREE_DATA(y_);
    FREE_DATA(binf_);
    FREE_DATA(bsup_);
    FREE_DATA(nbd_);
    FREE_DATA(wa_);
    FREE_DATA(iwa_);
    FREE_DATA(pdTmpLinear_);
    FREE_DATA(pdTmp_sub_);   
    FREE_DATA(Q_); 
    FREE_DATA(a_); 
    FREE_DATA(b_);  
    if ( handler_ != 0)
    {
        delete handler_;
        handler_ = 0;
    }
}

#ifdef DEBUG
void BiqModel::testEncodeDecode() {
    // Encode the original model
    AlpsEncoded* encoded = new AlpsEncoded();
    encode(encoded);

    // Create a new BiqModel instance and decode the encoded data into it
    BiqModel decodedModel = BiqModel(); // init empty model
    decodedModel.decodeToSelf(*encoded);


    // All the params used in bound from original 
    const bool bAddCuts = BiqPar_->entry(BiqParams::bAddCuts);

    const double dMinAlpha = BiqPar_->entry(BiqParams::dMinAlpha);
    const double dMinTol = BiqPar_->entry(BiqParams::dMinTol);
    const double dScaleAlpha = BiqPar_->entry(BiqParams::dScaleAlpha);
    const double dScaleTol = BiqPar_->entry(BiqParams::dScaleTol);

    const int MinNiter = BiqPar_->entry(BiqParams::nMinBoundingIter);
    const int maxNAiter = BiqPar_->entry(BiqParams::nMaxAlphaIter);
    const int nHeurRuns = BiqPar_->entry(BiqParams::nGoemanRuns);
    const int nMaxIter = BiqPar_->entry(BiqParams::nMaxBoundingIter);
    const int nMinAdded = BiqPar_->entry(BiqParams::nMinCuts);
    const int nMaxBFGSIter = BiqPar_->entry(BiqParams::nMaxBFGSIter);

    double dAlpha = BiqPar_->entry(BiqParams::dInitAlpha);
    double dTol = BiqPar_->entry(BiqParams::dInitTol);

    // All the params used in bound from decoded
    const bool bAddCuts_decodedModel = decodedModel.BiqPar_->entry(BiqParams::bAddCuts);

    const double dMinAlpha_decodedModel = decodedModel.BiqPar_->entry(BiqParams::dMinAlpha);
    const double dMinTol_decodedModel = decodedModel.BiqPar_->entry(BiqParams::dMinTol);
    const double dScaleAlpha_decodedModel = decodedModel.BiqPar_->entry(BiqParams::dScaleAlpha);
    const double dScaleTol_decodedModel = decodedModel.BiqPar_->entry(BiqParams::dScaleTol);

    const int MinNiter_decodedModel = decodedModel.BiqPar_->entry(BiqParams::nMinBoundingIter);
    const int maxNAiter_decodedModel = decodedModel.BiqPar_->entry(BiqParams::nMaxAlphaIter);
    const int nHeurRuns_decodedModel = decodedModel.BiqPar_->entry(BiqParams::nGoemanRuns);
    const int nMaxIter_decodedModel = decodedModel.BiqPar_->entry(BiqParams::nMaxBoundingIter);
    const int nMinAdded_decodedModel = decodedModel.BiqPar_->entry(BiqParams::nMinCuts);
    const int nMaxBFGSIter_decodedModel = decodedModel.BiqPar_->entry(BiqParams::nMaxBFGSIter);

    double dAlpha_decodedModel = decodedModel.BiqPar_->entry(BiqParams::dInitAlpha);
    double dTol_decodedModel = decodedModel.BiqPar_->entry(BiqParams::dInitTol);
    
    // Now check each param with assert
    std::printf("Checking Params ...\n");
    assert(bAddCuts == bAddCuts_decodedModel);
    assert(dMinAlpha == dMinAlpha_decodedModel);
    assert(dMinTol == dMinTol_decodedModel);
    assert(dScaleAlpha == dScaleAlpha_decodedModel);
    assert(dScaleTol == dScaleTol_decodedModel);
    assert(MinNiter == MinNiter_decodedModel);
    assert(maxNAiter == maxNAiter_decodedModel);
    assert(nHeurRuns == nHeurRuns_decodedModel);
    assert(nMaxIter == nMaxIter_decodedModel);
    assert(nMinAdded == nMinAdded_decodedModel);
    assert(nMaxBFGSIter == nMaxBFGSIter_decodedModel);
    assert(dAlpha == dAlpha_decodedModel);
    assert(dTol == dTol_decodedModel);
    std::printf("Params Match!\n:");


    // Now check that the decoded model has the same values as the original model
    //assert(originalModel == decodedModel);

    std::printf("Checking Model Data ...\n");
    assert(max_problem_ == decodedModel.max_problem_);
    assert(N_ == decodedModel.N_);
    assert(nVar_ == decodedModel.nVar_);
    assert(mB_original_ == decodedModel.mB_original_);
    assert(mA_ == decodedModel.mA_);
    assert(mB_ == decodedModel.mB_);

    bool bQmatch = true;
    bool bAmatch = true;
    bool bBmatch = true;

    for(int i = 0; i<N_*N_; ++i)
    {
        if(Q_[i] != decodedModel.Q_[i])
        {
            bQmatch = false;
        }
    }
    std::printf("bQmatch = %d\n", bQmatch);
    assert(bQmatch);

    for(int i = 0; i<mA_; ++i)
    {
        if(a_[i] != decodedModel.a_[i])
        {
            bAmatch = false;
        }
    }
    std::printf("bAmatch = %d\n", bAmatch);
    assert(bAmatch);

    for(int i = 0; i<mB_; ++i)
    {
        if(b_[i] != decodedModel.b_[i])
        {
            bBmatch = false;
        }
    }
    std::printf("bBmatch = %d\n", bBmatch);
    assert(bBmatch);
    
    // check sparse
    BiqSparseTriple bstThis;
    BiqSparseTriple bstDecoded;
    assert(Bs_.size() == decodedModel.Bs_.size());
    for(size_t i = 0; i < Bs_.size(); ++i)
    {
        assert(Bs_.at(i).size() == decodedModel.Bs_.at(i).size());
        for(size_t j = 0; j < Bs_.at(i).size(); ++j)
        {
            bstThis = (Bs_.at(i)).at(j);
            bstDecoded = (decodedModel.Bs_.at(i)).at(j);
            assert(bstThis.dVal_ == bstDecoded.dVal_);
            assert(bstThis.i_ == bstDecoded.i_);
            assert(bstThis.j_ == bstDecoded.j_);
    
        }

    }

    assert(As_.size() == decodedModel.As_.size());
    for(size_t i = 0; i < As_.size(); ++i)
    {
        assert(As_.at(i).size() == decodedModel.As_.at(i).size());
        for(size_t j = 0; j < As_.at(i).size(); ++j)
        {
            bstThis = (As_.at(i)).at(j);
            bstDecoded = (decodedModel.As_.at(i)).at(j);
            assert(bstThis.dVal_ == bstDecoded.dVal_);
            assert(bstThis.i_ == bstDecoded.i_);
            assert(bstThis.j_ == bstDecoded.j_);
    
        }

    }

    std::cout << "Encode/decode test passed!\n";
}
#endif

AlpsReturnStatus BiqModel::encode(AlpsEncoded * encoded) const
{
    AlpsReturnStatus status = AlpsReturnStatusOk;

    //------------------------------------------------------
    // Encode Alps part.
    //------------------------------------------------------

    status = AlpsModel::encode(encoded);
 
    //------------------------------------------------------
    // Pack up params.
    //------------------------------------------------------
    BiqPar_->pack(*encoded);
    //
    //------------------------------------------------------
    // Encode Biq part.
    //------------------------------------------------------
    // copy data from Sparse std:vectors
    int pos;
    int sizeQs, sizeAs, sizeBs;
    int NN = N_*N_;
    int* i_Qs;
    int* j_Qs;
    int* i_As;
    int* j_As;
    int* i_Bs;
    int *j_Bs;
    int* indexAs;
    int* indexBs;
    double* dVal_Qs;
    double* dVal_As;
    double* dVal_Bs;


    sizeQs = Qs_.size();
    i_Qs = new int[sizeQs];
    j_Qs = new int[sizeQs];
    dVal_Qs = new double[sizeQs];
    pos = 0;
    for(auto &it: Qs_)
    {
        i_Qs[pos] = it.i_;
        j_Qs[pos] = it.j_;
        dVal_Qs[pos] = it.dVal_;
        ++pos;
    }

    // need to figure out the size of As_
    sizeAs = 0;
    for(auto& it: As_)
    {
        sizeAs += it.size();
    }
    // if we have anything in As_
    if(sizeAs > 0)
    {
        // allocate memory needed;
        indexAs = new int[mA_];
        i_As = new int[sizeAs];
        j_As = new int[sizeAs];
        dVal_As = new double[sizeAs];

        //loop over sparse data structure
        pos = 0;
        for(int i = 0; i < mA_; ++i)
        {
            // set index to size maybe use better name
            indexAs[i] = As_.at(i).size();
            // loop over the vecotr of BiqTriples at i
            for(auto& it: As_.at(i))
            {
                i_As[pos] = it.i_;
                j_As[pos] = it.j_;
                dVal_As[pos] = it.dVal_;
                ++pos;
            }
        }
    }
    else
    {
        indexAs = 0;
        i_As = 0;
        j_As = 0;
        dVal_As = 0;
    }


    // need to figure out the size of Bs_
    sizeBs = 0;
    for(auto& it: Bs_)
    {
        sizeBs += it.size();
    }
    // if we have anything in Bs_
    if(sizeBs > 0)
    {
        // allocate memory needed;
        indexBs = new int[mB_];
        i_Bs = new int[sizeBs];
        j_Bs = new int[sizeBs];
        dVal_Bs = new double[sizeBs];

        // loop over sparse data structure
        pos = 0;
        for(int i = 0; i < mB_; ++i)
        {
            // set index to sie maybe use better name
            indexBs[i] = Bs_.at(i).size();
            // loop over the vecotr of BiqTriples at i
            for(auto& it: Bs_.at(i))
            {
                i_Bs[pos] = it.i_;
                j_Bs[pos] = it.j_;
                dVal_Bs[pos] = it.dVal_;
                ++pos;
            }
        }
    }
    else
    {
        indexBs = 0;
        i_Bs = 0;
        j_Bs = 0;
        dVal_Bs = 0;
    }

    // write the easy stuff Int, Double, Bool ...
    encoded->writeRep(max_problem_);
    encoded->writeRep(N_);
    encoded->writeRep(nVar_);
    encoded->writeRep(mB_original_); // needed for computing objective value
    encoded->writeRep(mA_);
    encoded->writeRep(mB_);
    encoded->writeRep(iPastBest_);
    encoded->writeRep(bHasPastBest_);

    // write the arrays
    encoded->writeRep(Q_,NN);
    encoded->writeRep(a_,mA_);
    encoded->writeRep(b_,mB_);

    // finally the data structures like std::vectors
    encoded->writeRep(sizeQs);
    encoded->writeRep(i_Qs,sizeQs);
    encoded->writeRep(j_Qs,sizeQs);
    encoded->writeRep(dVal_Qs,sizeQs);

    encoded->writeRep(sizeAs);
    encoded->writeRep(indexAs,mA_);
    encoded->writeRep(i_As,sizeAs);
    encoded->writeRep(j_As,sizeAs);
    encoded->writeRep(dVal_As,sizeAs);

    encoded->writeRep(sizeBs);
    encoded->writeRep(indexBs,mB_);
    encoded->writeRep(i_Bs,sizeBs);
    encoded->writeRep(j_Bs,sizeBs);
    encoded->writeRep(dVal_Bs,sizeBs);

    return status;
}

AlpsReturnStatus BiqModel::decodeToSelf(AlpsEncoded & encoded)
{
    AlpsReturnStatus status = AlpsReturnStatusOk;

    //------------------------------------------------------
    // Dencode Alps part.
    //------------------------------------------------------
    status = AlpsModel::decodeToSelf(encoded);

    //------------------------------------------------------
    // Dencode Params
    //------------------------------------------------------
    BiqPar_->unpack(encoded);
    //------------------------------------------------------
    // Decode Biq part.
    //------------------------------------------------------
    
    int pos, sizeTriple;
    int sizeQs, sizeAs, sizeBs;
    int* i_Qs = nullptr;
    int* j_Qs = nullptr;
    int* i_As = nullptr;
    int* j_As = nullptr;
    int* i_Bs = nullptr;
    int* j_Bs = nullptr;
    int* indexAs = nullptr;
    int* indexBs = nullptr;
    //int i_tmp, j_tmp;
    double* dVal_Qs = nullptr;
    double* dVal_As = nullptr;
    double* dVal_Bs = nullptr;
    double* Q_copy = nullptr;
    double* a_copy = nullptr;
    double* b_copy = nullptr;
    //double dVal_tmp;
    int NN;
    //size_t i;
    int i;
    
    // Read the data in the same order it was written
    encoded.readRep(max_problem_);
    encoded.readRep(N_);
    encoded.readRep(nVar_);
    encoded.readRep(mB_original_);
    encoded.readRep(mA_);
    encoded.readRep(mB_);
    encoded.readRep(iPastBest_);
    encoded.readRep(bHasPastBest_);

    // allocate for readRep
    NN = N_*N_;
    Q_ = new double[NN];
    a_ = new double[mA_];
    b_ = new double[mB_];

    encoded.readRep(Q_copy,NN);
    encoded.readRep(a_copy,mA_);
    encoded.readRep(b_copy,mB_);

    
    // copy Q_ 
    if(NN > 0)
    {
        for(i = 0; i < NN; ++i)
        {
            Q_[i] = Q_copy[i];
        }
    }
    // copy a_
    if(mA_ > 0) 
    {
        for(i = 0; i < mA_; ++i)
        {
            a_[i] = a_copy[i];
        }
    }
    // copy b_
    if(mB_ > 0)
    {
        for(i = 0; i < mB_; ++i)
        {
            b_[i] = b_copy[i];
        }
    }

    encoded.readRep(sizeQs);
    encoded.readRep(i_Qs,sizeQs);
    encoded.readRep(j_Qs,sizeQs);
    encoded.readRep(dVal_Qs,sizeQs);

    encoded.readRep(sizeAs);
    encoded.readRep(indexAs,mA_);
    encoded.readRep(i_As,sizeAs);
    encoded.readRep(j_As,sizeAs);
    encoded.readRep(dVal_As,sizeAs);

    encoded.readRep(sizeBs);
    encoded.readRep(indexBs,mB_);
    encoded.readRep(i_Bs,sizeBs);
    encoded.readRep(j_Bs,sizeBs);
    encoded.readRep(dVal_Bs,sizeBs);

    Qs_.resize(sizeQs);
    pos = 0;
    for(auto &it: Qs_)
    {
        it.i_ = i_Qs[pos];
        it.j_ = j_Qs[pos];
        it.dVal_ = dVal_Qs[pos];
        ++pos;
    }

    // now that data is read build As_ and Bs_
    As_.resize(mA_);
    pos = 0;
    for(int i = 0; i < mA_; ++i)
    {
        //
        sizeTriple = indexAs[i];
        As_.at(i).resize(sizeTriple);
        for(auto& it: As_.at(i))
        {
            it.i_ = i_As[pos];
            it.j_ = j_As[pos];
            it.dVal_ = dVal_As[pos];
            ++pos;
        }
    }

    Bs_.resize(mB_);
    pos = 0;
    for(int i = 0; i < mB_; ++i)
    {
        //
        sizeTriple = indexBs[i];
        Bs_.at(i).resize(sizeTriple);
        for(auto& it: Bs_.at(i))
        {
            it.i_ = i_Bs[pos];
            it.j_ = j_Bs[pos];
            it.dVal_ = dVal_Bs[pos];
            ++pos;
        }
    }

    // For now this is init() with out the data set above
    // get the params needed for initilization 
    const bool bProdCons = BiqPar_->entry(BiqParams::bAddProductConstraints);
    const int nCuts = BiqPar_->entry(BiqParams::nCuts);
    const int MaxNineqAdded = BiqPar_->entry(BiqParams::MaxNineqAdded);

    ISUPPZ_ = new int[2*N_];
    X_ = new double[N_*N_];
    W_ = new double[N_];
    Z_ = new double[N_*N_];

    /* DSYEVR optimal workspace query */
    M_ = 0;
    nDWORK_ = 26 * N_;
    nIWORK_ =  10 * N_;
    DWORK_ = new double[nDWORK_];
    IWORK_ = new int[nIWORK_]; 
    
    a_sub_ = new double[mA_];
    b_sub_ = new double[mB_];
    Q_sub_ = new double[N_*N_];
    vOffset_.resize(nVar_, 0);

    /* L-BFGS-B Data */
    int nMax = mB_ + mA_ + MaxNineqAdded;
    int wa_length = (2 * mmax + 5) * nMax + 11 * mmax * mmax + 8 * mmax;
    int iwa_length = 3 * nMax;

    f_ = 0.0;
    g_ = new double[nMax];
    y_ = new double[nMax];
    binf_ = new double[nMax];
    bsup_ = new double[nMax];
    nbd_ = new int[nMax];

    wa_ = new double[wa_length];
    iwa_ = new int[iwa_length];

    /* NK: Should we add diagonal and product constraints here? */
    //AddDiagCons();
    //if(bProdCons) AddProdCons();

    AllocSubCons();
    SetConSparseSize();
    // Set cut data 
    Cuts_.resize(MaxNineqAdded); // or start small and add space when needed
    container.reserve(nCuts);

    nIneq_ = 0;

    // set some memory for the sTemp_sub_ matrix
    sTmp_sub_.resize(N_*N_);
    pdTmp_sub_ = new double[N_*N_];
    pdTmpLinear_ = new double[nVar_];

    viSolution_1_.resize(N_);
    viSolution_2_.resize(N_);

    vdFracSol_.resize(nVar_);

    vdGainsZero_.resize(nVar_);
    vdGainsOne_.resize(nVar_);
    vbGreedyIndex_.resize(nVar_);

    return status;
}

AlpsKnowledge * BiqModel::decode(AlpsEncoded & encoded) const {
    std::cerr << "Not implemented!" << std::endl;
    throw std::exception();
}


void BiqModel::AddDiagCons()
{
    // resize ?
    // then loop 
    std::printf("---------------------- Adding Diag Cons to Bs_\n");

    // size of Bs_ pre diag cons use to shift entires of b_
    int nCons = Bs_.size();

                                 //  i, j, dVal
    Sparse bstTemp = {BiqSparseTriple(0,0,1.0)};
    for(int i = 0; i < N_; ++i)
    {
        // fill bstTemp 
        bstTemp.at(0).i_ = i;
        bstTemp.at(0).j_ = i;
        // them push bstTemp to Bs_
        Bs_.push_back(bstTemp);
        b_[i + nCons] = 1.0;
        //PrintSparseMatrix(bstTemp);
    }

    mB_ = Bs_.size();
}

/// @brief 
///
///
void BiqModel::AddProdCons()
{
    Sparse bstTemp;
    int jj;
    int posProdCon = Bs_.size();
    double dTemp;
    double dLinear = 0;

    std::printf("---------------------- Adding Prod Cons to Bs_\n");
    //std::printf("mB_ = %d \t Size(Bs_) = %d \t Number of Diag Cons = %d \t Equality Cons = %d\n", mB_, Bs_.size(), N_, iEqualityIsLinear_.size());

    // k = positon of equality con
    for(size_t k = 0; k < iEqualityIsLinear_.size(); ++k)
    {
        // r = variable index
        for(int r = 0; r < nVar_; ++r)
        {

            //std::printf("pos of the con %d\n",iEqualityIsLinear_.at(k));
            for(auto & it: Bs_.at(iEqualityIsLinear_.at(k)))
            {
                jj = it.j_;
                dTemp = it.dVal_;
                // dTemp xj xr
                // if linear 
                if(jj < nVar_)
                {
                    dTemp *= 2;
                    // on the diagonal
                    if(jj == r)
                    {
                        bstTemp.push_back(BiqSparseTriple(jj,jj,dTemp));
                    }
                    else
                    {
                        bstTemp.push_back(BiqSparseTriple(r,jj,0.5*dTemp));
                    }
                }
                // case we are in const
                else
                {
                    dLinear = dTemp;
                    
                }

            }
            //std::printf("Const = %f     RHS = %f \n", dLinear,  b_[k]);
            dLinear -= b_[k];
            if(dLinear != 0.0)
            {
                bstTemp.push_back(BiqSparseTriple(nVar_, r,0.5*dLinear));
            }
            
            // add data
            Bs_.push_back(bstTemp);
            b_[posProdCon] = 0;
            
            /*
            std::printf("\n og con");
            PrintSparseMatrix(Bs_.at(iEqualityIsLinear_.at(k)));
            std::printf("\n prod con");
            PrintSparseMatrix(Bs_.at(posProdCon));
            std::printf("\n");
            */
            // reset tmp data
            bstTemp.clear();
            dLinear = 0;
            // 
            posProdCon++;
        }
    }
    mB_ = Bs_.size();
}
/// @brief use this method
/// to set the size of each sub
/// con sparse matrix
void BiqModel::SetConSparseSize()
{

    for(int i = 0; i < mB_; ++i)
    {
        Bs_sub_.at(i).resize(Bs_.at(i).size());
    }

    
    for(int i = 0; i < mA_; ++i)
    {
        As_sub_.at(i).resize(As_.at(i).size());
    }

}

/// @brief
/// Need to set the size of the array of cons
/// Then need to set size for con
void BiqModel::AllocSubCons()
{
    Bs_sub_.resize(mB_);
    mB_sub_ = Bs_sub_.size();
    As_sub_.resize(mA_);
    mA_sub_ = As_sub_.size();
}

/// @brief 
/// @return 
double BiqModel::SDPbound(std::vector<BiqVarStatus> vbiqVarStatus , bool bRoot)
{


    clock_t start = clock();
    // reset cuts
    Map_.clear();
    nIneq_ = 0;
    //
    int i;
    int iStatus;
    int nbit = 0;
    int nFuncEvals = 0;
    int nAdded = 0;
    int nSubtracted = 0;
    int nbitalpha = 0;
    int nMaxFuncEvals = 0;    
    
    double dRetBound = 0.0;
    double dPrev_f = MAXFLOAT;
    double  bestVal;
    double dGap;
    
    bool prune = false;
    bool bGiveUp = false;
    bool bStopSDPBound = false;

    const bool bPrintTable = false;

    char cExitReason[100];
    // param
    //if (bRoot) BiqPar_->print();

    const bool bAddCuts = BiqPar_->entry(BiqParams::bAddCuts);
    const bool  bSolutionProvided = BiqPar_->entry(BiqParams::bSolutionProvided);
    const double dMinAlpha = BiqPar_->entry(BiqParams::dMinAlpha);
    const double dMinTol = BiqPar_->entry(BiqParams::dMinTol);
    const double dScaleAlpha = BiqPar_->entry(BiqParams::dScaleAlpha);
    const double dScaleTol = BiqPar_->entry(BiqParams::dScaleTol);

    const int MinNiter = BiqPar_->entry(BiqParams::nMinBoundingIter);
    const int maxNAiter = BiqPar_->entry(BiqParams::nMaxAlphaIter);
    const int nGoemanRuns = BiqPar_->entry(BiqParams::nGoemanRuns);
    const int nMaxIter = BiqPar_->entry(BiqParams::nMaxBoundingIter);
    const int nMinAdded = BiqPar_->entry(BiqParams::nMinCuts);
    int nMaxBFGSIter = BiqPar_->entry(BiqParams::nMaxBFGSIter);
    const int MaxNineqAdded = BiqPar_->entry(BiqParams::MaxNineqAdded);


    double dAlpha = BiqPar_->entry(BiqParams::dInitAlpha);
    double dTol = BiqPar_->entry(BiqParams::dInitTol);
    // Biq.par

    //int len_y = mA_ + mB_ + nIneq_;
    int len_y = mA_ + mB_ + MaxNineqAdded;
    int nHeurRuns = nGoemanRuns*nVar_;
    
    int iTreeDepth = GetOffset(vbiqVarStatus);

    if(bSolutionProvided)
    {
        addProvidedSol();
    }

    // 
    if(nVar_sub_ == 0)
    {
        dRetBound = Q_sub_[0];
        //std::printf("On a leaf bound is: %f\n", dRetBound);
        return dRetBound;
    }

    for(i = 0; i < len_y; ++i)
    {
        y_[i] = 0.0;
    }

    sim(dAlpha);
    for(i = 0; i < nMaxIter; ++i)
    {

        // update iteration countrer
        ++nbitalpha;

        dPrev_f = f_;
        // call BFGS solver
        iStatus = CallLBFGSB(dAlpha, dTol, nbit);
        // add nbits to the sum of function evaluations
        nFuncEvals += nbit;
        if(nbit > nMaxFuncEvals)
        {
            nMaxFuncEvals = nbit;
        }
        dRetBound = (max_problem_) ? floor(f_) : ceil(-f_);
        
        prune = pruneTest(dRetBound);

        if(!prune)
        {
            GWheuristic(nHeurRuns, vbiqVarStatus);
            prune = pruneTest(dRetBound);
        }


        bestVal = getBestVal();

        dGap = fabs(bestVal - dRetBound);

        if( i + 1 > MinNiter) 
        {
            bGiveUp = (fabs(f_ - dPrev_f) < 1.0) &&
                      (dAlpha  < 1e-4)           && 
                      (dGap > 2.0);
        }

        if(!prune && !bGiveUp && bAddCuts)
        {
            UpdateInequalities(nAdded, nSubtracted);
        }
        else
        {
            nAdded = 0;
            nSubtracted = 0;
        }

        bStopSDPBound = (dAlpha == dMinAlpha) && nAdded == 0;
    

        // check if we are done
        if(prune || bGiveUp || bStopSDPBound || iStatus == -1 || nbit >= nMaxBFGSIter)
        {
            clock_t end = clock();
            double duration = double(end - start) / CLOCKS_PER_SEC;
            // Build the messsage to go to output
            if(prune) std::strcpy(cExitReason,"0"); 
            else if(bGiveUp) std::strcpy(cExitReason,"1");
            else if(bStopSDPBound) std::strcpy(cExitReason,"2");
            else if(iStatus == -1) std::strcpy(cExitReason,"3");
            else if(nbit >= nMaxBFGSIter) std::strcpy(cExitReason,"4");
            else std::strcpy(cExitReason,"5");

            // BIQ_BOUND_DATA > "%f\t%d\t%d\t%d\t%s\t%f\t%f\t%d\t%d\t%d\t%f"   
            handler_->message(BIQ_BOUND_DATA, messages_)
                << duration
                << i+1
                << iTreeDepth
                << broker()->getNumNodesProcessed()
                << cExitReason
                << dRetBound
                << dGap
                << nFuncEvals
                << nMaxFuncEvals
                << nIneq_
                << dAlpha;
                
            break;
        }
        // update parameters
        
        if(nAdded < nMinAdded || nbitalpha > maxNAiter)
        {
            nbitalpha = 0;
            dAlpha *= dScaleAlpha;
            dAlpha = (dAlpha < dMinAlpha) ? dMinAlpha : dAlpha;
            dTol *= dScaleTol;
            dTol = (dTol < dMinTol) ? dMinTol : dTol;
        }

        
    }

    dRetBound = (max_problem_) ? f_ : -f_;

    return dRetBound;
}

void BiqModel::PrintBoundingTable(int iter, int nBit, int nAdded, int nSubtracted, double dAlpha, double dTol, double dMinAllIneq, double dGap /*double dTime*/)
{
    int len_y = mA_ + mB_ + nIneq_;  
    int incy = 1;
    double dYnorm = dnrm2_(len_y, y_, incy);

    if(iter == 1)
    {
        std::printf("===============================================================================================\n");
        std::printf("%4s  %6s  %8s  %5s  %5s  %4s  %5s  %5s  %5s  %6s  %6s  %6s  %6s\n", 
                    "Iter", 
                    "Time", 
                    //"gap", 
                    "bound",
                    "alpha", 
                    "tol", 
                    "nbit", 
                    "Enorm", 
                    "Inorm", 
                    "ynorm", 
                    "minCut", 
                    "NCut", 
                    "NSub", 
                    "NAdd");
        std::printf("===============================================================================================\n");
    }

    std::printf("%4d  %6.1f  %8.2f  %5.0e  %5.0e  %4d  %5.0e  %5.0e  %5.0e  %6.0e  %6d   -%-5d  +%-5d\n", 
                iter, 
                0.0, /* TODO time*/
                //dGap, 
                f_, 
                dAlpha, 
                dTol, 
                nBit, 
                gradEnorm_,
                gradInorm_,
                dYnorm, /* TODO ynorm*/
                dMinAllIneq, 
                nIneq_, 
                nSubtracted,
                nAdded
    );

}
int BiqModel::CallLBFGSB(double dAlpha, double dTol, int &nbit)
{
    const int nMaxBFGSIter = BiqPar_->entry(BiqParams::nMaxBFGSIter);

    int retStatus = 1;
    int mem = mmax;
    int i;
    int N = nVar_sub_ + 1;
    bool bStopBFGS = false;
    bool bPrune;
    double factr;
    double pgtol;
    double dMinTemp;
    double dBound;
    char task[60];
    int iprint = -10;
    char csave[60];
    int lsave[4];
    int isave[44];
    double dsave[29];
    int mA_sub = As_sub_.size();
    int mB_sub = Bs_sub_.size();
    // compute the number of variables in BFGS
    int len_y = mA_sub + mB_sub + nIneq_;
    // set task = "START" and the rest of it to empty char
    strcpy(task, "START");
    for(i = 5; i < 60; ++i) task[i] = ' ';

    
    // Initialize y
    int tmpPos = 0;
    for(auto it = Cuts_.begin(); it < Cuts_.begin() + nIneq_; ++it)
    {
        y_[mA_sub + mB_sub + tmpPos++] = it->y_;
    }


    // compute the bounds for y
    //// equalities
    for(i = 0; i < mB_sub; ++i)
    {
        nbd_[i] = 0; 
    }
    //// inequalities
    for(i = mB_sub; i < len_y; ++i)
    {
        nbd_[i] = 1; 
        binf_[i] = 0.0;
    }

    
    nbit = 0; 
    // BFGS main loop
    while (!bStopBFGS)
    {
        factr = 0.0;
        pgtol = 0.0;
        // this calls the main L-BFGS-B function
        setulb_(len_y, mem, y_, binf_, bsup_, nbd_, f_, g_,
                factr, pgtol, wa_, iwa_, task, iprint,
                csave, lsave, isave, dsave);
    
        if(strncmp(task, "FG", 2) == 0)
        {
           // L-BFGS-B requesting new f and g
           sim(dAlpha);
           ++nbit; 
        }
        else if(strncmp(task, "NEW_X", 5) == 0)
        {
            // L-BFGS-B found a new x
            // test if we should stop

            // Compute rhe Infinity-norm of gradE = g[g[0:mB_-1]
            gradEnorm_ = 0.0;
            for(i = 0; i < mB_sub; ++i)
            {
                gradEnorm_ = (gradEnorm_ < fabs(g_[i])) ? fabs(g_[i]) : gradEnorm_;
            }
            gradEnorm_ /= scaleEq;

            // Compute Infinity-norm of gradI = min(g[mB_sub:nLBFGSB_Vars-1], 0.0)
            gradInorm_ = 0.0;
            for(i = mB_sub; i < len_y; ++i)
            {
                dMinTemp = (g_[i] > 0.0) ? 0.0 : g_[i];
                gradInorm_ = (gradInorm_ < fabs(dMinTemp)) ? fabs(dMinTemp) : gradInorm_;
            }
            gradInorm_ /= scaleIneq;

            dBound = (max_problem_) ? f_  : - f_;
            bPrune = pruneTest(dBound);
            
            //std::printf("Grad E: %20.16f \t Grad I: %20.16f\n",gradEnorm_, gradInorm_);
            if(bPrune
               || (gradEnorm_ < dTol && gradInorm_ < dTol)
               || nbit >= nMaxBFGSIter)
               {
                /*
                std::printf("stopping ... prune: %d \t E norm %f tol %f \t E norm %f tol %f \tbits %d \n",
                             bPrune,gradEnorm_,dTol,gradInorm_,dTol,nbit);
                */
                strcpy(task, "STOP");
                for(i = 4; i < 60; ++i) task[i] = ' ';
               }
        }
        else if(strncmp(task, "STOP", 4) == 0)
        {
            bStopBFGS = true;
            retStatus = 0; 
        }
        else
        {
            bStopBFGS = true;
        }
    }

    // Scale: X = X / alpha
    double dAlphaInv = 1.0 / dAlpha;
    int N2 = N*N;
    int incx = 1; 
    dscal_(N2, dAlphaInv, X_, incx);
    // Update Z so that X = Z*Z'
    double scb = 10000.0 * dAlpha;
    double sca = 0.01/ sqrt(scb);
    int NM = N * M_;
    int incz = 1;
    dscal_(NM, sca, Z_, incz); 
    //std::printf("status => %s\n", task);
    return retStatus;
}

/// @brief 
/// @param alpha 
void BiqModel::sim(double alpha)
{
    int N  = nVar_sub_ + 1;
    int N2 = N*N;
    double dAlphaInv = 1.0 / alpha;
    // prepare to copy Q to X_
    // memory spacing for the two arrays
    int INCQ = 1; 
    int INCX = 1;
    // copy Q to X 
    dcopy_(N2, Q_sub_, INCQ, X_, INCX);

    // call the A and B methods transposed
    B(TRANSP, alpha);
    A(TRANSP, alpha);
    ProjSDP();


    f_ = f_ * 1000.0 * (dAlphaInv / 1000.0); // TODO see what happens 

    // call the A and B methods not transposed
    B(NOTRANSP, alpha);
    A(NOTRANSP, alpha);

    f_ += alpha * (0.5 * static_cast<double>(N2));

    //print_vector(g_,mA_+mB_+nIneq_);
    //print_symmetric_matrix(X_,N);
}
 
/// @brief 
///
///
///
void BiqModel::ProjSDP()
{
    double bornesup;
    double borneinf = 1e-8;

    char UPLO = 'L';
    int N = nVar_sub_ + 1;
    int LDX = N;

    // Operator norm for X to get a upper bound of the spectral radius of X
    // |X|_2 \leq \sqrt{ |X|_1 |X|_inf }  (Holder's inequality)
    //          = |X|_1 = |X|_inf  (since X is symmetric)
    //
    // Frobenius norm is also an upper bound on the spectral radius of X:
    //      |X|_2 <= |X|_F
    char NORM = 'I'; // in
    double norminf = dlansy_(NORM, UPLO, N, X_, LDX, DWORK_);
    NORM = 'F';
    double normfro = dlansy_(NORM, UPLO, N, X_, LDX, DWORK_);

    // bornesup = min(norminf, normfro)
    bornesup = (norminf < normfro) ? norminf : normfro;

    // Ensure that borneinf <= bornesup.
    if (bornesup < borneinf) {
        bornesup = 2.0 * borneinf;
    }
    
    //printf("\nX = \n"); print_symmetric_matrix(X, N);

    /* Compute the positive eigenvalues and associated eigenvectors of X.
     *
     * The M columns of Z will contain the orthonormal eigenvectors of the 
     * matrix X corresponding to the positive eigenvalues, the i-th column 
     * of Z holding the eigenvector associated with W[i
     */
    char JOBZ = 'V';    // eigenvalues and eigenvectors
    char RANGE = 'V';
    double VL = borneinf;
    double VU = bornesup;
    int IL = 0;
    int IU = 0;
    double ABSTOL = 1e-8; // TODO could be a param 
    int LDZ = N;
    int INFO;

    dsyevr_(JOBZ, RANGE, UPLO, N, X_, LDX, VL, VU, IL, IU, ABSTOL,
            M_, W_, Z_, LDZ, ISUPPZ_, DWORK_, nDWORK_, IWORK_, nIWORK_, INFO);

   //std::printf("M = %d\n", M_);
   // Check if the eigensolver failed (i.e., if INFO != 0)
    if (INFO) {
        VL = 0.0;
        ABSTOL = 0.0;
        dsyevr_(JOBZ, RANGE, UPLO, N, X_, LDX, VL, VU, IL, IU, ABSTOL,
            M_, W_, Z_, LDZ, ISUPPZ_, DWORK_, nDWORK_, IWORK_, nIWORK_, INFO);
        if (INFO) {
            fprintf(stderr, 
                    "Error: eigenvalue computation failed (INFO = %d)\n", 
                    INFO);
            exit(1);
        }
    }

    // Compute ff = 0.5*||X_+||^2 = 0.5*||W||^2
    int INCW = 1;
    double normW = dnrm2_(M_, W_, INCW);
    f_ = 0.5*(normW*normW);
    // Compute Z = Z*Diag(W)^{1/2}
    int INCZ = 1;
    double temp;
    for (int j = 0; j < M_; j++) {
        // Scale jth column of Z by sqrt(W[j])
        temp = sqrt(W_[j]);
        dscal_(N, temp, Z_ + j*N, INCZ);
    }

    char TRANS = 'N';
    double ALPHA = 1.0;
    double BETA = 0.0;

    /* X = ALPHA*Z*Z^T + BETA*X = Z*Z^T
     * Only writes to lower-triangular part of X.
     * When M = 0, we will obtain X = 0.
     */
    dsyrk_(UPLO, TRANS, N, M_, ALPHA, Z_, LDZ, BETA, X_, LDX);

    /*
    print_symmetric_matrix(X_, N_);
    */
}

/// @brief 
/// Constraint:  A(X) <= a
///
/// Evaluates:
///      if mode == TRANSP:
///          X = X + scaleIneq * A^*(y_)
///
///      if mode == NOTRANSP:
///          ff = ff + scaleIneq * <a,y_>
///          g_ = scaleIneq*( a + A(X/alpha) )
///
/// @param mode
/// @param y_ 
/// @param g_ 
/// @param alpha 
void BiqModel::A(int mode, double alpha)
{

    double dAlphaInv = 1.0 / alpha;
    double dTemp;
    int mA_sub = As_sub_.size();
    int mB_sub = Bs_sub_.size();
    int ineqCon;
    int N = nVar_sub_ + 1;

    if(mode == TRANSP)
    {
        // models inequalities 
        for(ineqCon = 0; ineqCon < mA_sub; ++ineqCon) // try a range based loop here too
        {
            dTemp = scaleIneq * y_[mB_sub + ineqCon];

            for(auto& it : As_sub_[ineqCon])
            {
                if(it.i_ == it.j_)
                {
                    X_[it.i_ + it.j_ * N] -= dTemp * it.dVal_;
                }
                else
                {
                    X_[it.i_ + it.j_ * N] -= dTemp * it.dVal_;
                    X_[it.j_ + it.i_ * N] -= dTemp * it.dVal_;
                } 
            }
        }
        // cut inequalities
        ineqCon = 0;
        for(auto it = Cuts_.begin(); it < Cuts_.begin() + nIneq_; ++it)
        {
            dTemp = 0.5 * scaleIneq * y_[mB_sub + mA_sub + ineqCon];
            
            switch (it->type_)
            {

            case ONE:
            {
                X_[it->i_ + it->j_ * N] += dTemp;
                X_[it->i_ + it->k_ * N] += dTemp;
                X_[it->j_ + it->k_ * N] += dTemp;
            } break;
            
            case TWO:
            {
                X_[it->i_ + it->j_ * N] += dTemp;
                X_[it->i_ + it->k_ * N] -= dTemp;
                X_[it->j_ + it->k_ * N] -= dTemp;
            } break;

            case THREE:
            {
                X_[it->i_ + it->j_ * N] -= dTemp;
                X_[it->i_ + it->k_ * N] += dTemp;
                X_[it->j_ + it->k_ * N] -= dTemp;
            } break;

            default: // FOUR
            {
                X_[it->i_ + it->j_ * N] -= dTemp;
                X_[it->i_ + it->k_ * N] -= dTemp;
                X_[it->j_ + it->k_ * N] += dTemp;
            } break;
            }
            ineqCon++;
        }

    }
    else
    {
        for(ineqCon = 0; ineqCon < mA_sub; ++ineqCon)
        {

            f_ += scaleIneq * a_sub_[ineqCon] * y_[mB_sub + ineqCon];

            g_[mB_sub + ineqCon] = 0.0;

            for(auto& it : As_sub_[ineqCon])
            {
                if(it.i_ == it.j_)
                {
                    g_[mB_sub + ineqCon] -= it.dVal_ * X_[it.i_ + it.j_ * N];
                }
                else
                {
                    g_[mB_sub + ineqCon] -= 2.0 * it.dVal_ * X_[it.i_ + it.j_ * N];
                }
            }
            g_[mB_sub + ineqCon] *= dAlphaInv;
            g_[mB_sub + ineqCon] += a_sub_[ineqCon];
            g_[mB_sub + ineqCon] *= scaleIneq;
        }

        // cut inequalities
        ineqCon = 0;
        
        for(auto it = Cuts_.begin(); it < Cuts_.begin() + nIneq_; ++it)
        {
            f_ += scaleIneq * y_[mB_sub + mA_sub + ineqCon];

            
            switch (it->type_)
            {

            case ONE:
            {
                dTemp = X_[it->i_ + it->j_ * N] + X_[it->i_ + it->k_ * N] + X_[it->j_ + it->k_ * N];
            } break;
            
            case TWO:
            {
                dTemp = X_[it->i_ + it->j_ * N] - X_[it->i_ + it->k_ * N] - X_[it->j_ + it->k_ * N];
            } break;

            case THREE:
            {
                dTemp = -X_[it->i_ + it->j_ * N] + X_[it->i_ + it->k_ * N] - X_[it->j_ + it->k_ * N];
            } break;

            default: // FOUR
            {
                dTemp = -X_[it->i_ + it->j_ * N] - X_[it->i_ + it->k_ * N] + X_[it->j_ + it->k_ * N];
            } break;
            }

            g_[mB_sub + mA_sub + ineqCon] = (dTemp * dAlphaInv + 1.0) * scaleIneq;
            
            ineqCon++;
        }
    }
}

/// @brief 
/// Constraint:  B(X) = b
///
/// Evaluates:
///      if mode == TRANSP:
///          X = X - scaleEq * B^*(y_)
///
///      if mode == NOTRANSP:
///          ff = ff + scaleEq * <b,y_>
///          g_ = scaleEq*( b - B(X/alpha) )
///
/// @param mode
/// @param y_ 
/// @param g_ 
/// @param alpha 
void BiqModel::B(int mode, double alpha)
{

    double dAlphaInv = 1.0 / alpha;
    double dTemp;
    int N = nVar_sub_ + 1;
    int mB_sub = Bs_sub_.size();
    //std::vector<BiqTriInequality>::iterator it = Bs_.begin();

    if(mode == TRANSP)
    {
        for(int eqCon = 0; eqCon < mB_sub; ++eqCon)
        {
            dTemp = scaleEq * y_[eqCon];
            for(auto& it : Bs_sub_[eqCon])
            {
                if(it.i_ == it.j_)
                {
                    X_[it.i_ + it.j_ * N] -= dTemp * it.dVal_;
                }
                else
                {
                    X_[it.i_ + it.j_ * N] -= dTemp * it.dVal_; // why not * 2.0 and skip the next line 
                    X_[it.j_ + it.i_ * N] -= dTemp * it.dVal_;
                }
            }
        }
    }
    else
    {
        std::vector<Sparse>::iterator itBs = Bs_sub_.begin();
        
        for(int eqCon = 0; eqCon < mB_sub; ++eqCon)   
        {
            f_ +=  scaleEq * b_sub_[eqCon] * y_[eqCon];

            g_[eqCon] = 0.0;

            for(auto &it : *itBs)
            {
                if(it.dVal_ == 0) continue;

                if(it.i_ == it.j_)
                {
                    g_[eqCon] -= it.dVal_ * X_[it.i_ + it.j_ * N];
                }
                else
                {
                    g_[eqCon] -= 2.0 * it.dVal_ * X_[it.i_ + it.j_ * N];
                }
            }
            g_[eqCon] *= dAlphaInv;
            g_[eqCon] += b_sub_[eqCon];
            g_[eqCon] *= scaleEq;

            

            ++itBs;
        }
    }

}

/// @brief 
void BiqModel::CreateSubProblem(std::vector<BiqVarStatus> vbiqVarStatus)
{
    int nFixed = GetOffset(vbiqVarStatus);

    nVar_sub_ = nVar_ - nFixed;

    for(auto it = Bs_sub_.begin(); it < Bs_sub_.end(); ++it)
    {
        it->clear();
    }
    for(auto it = As_sub_.begin(); it < As_sub_.end(); ++it)
    {
        it->clear();
    }

    BuildConstraints(mB_, b_, Bs_, b_sub_, Bs_sub_,
                        vbiqVarStatus, nFixed); 


    BuildConstraints(mA_, a_, As_, a_sub_, As_sub_,
                        vbiqVarStatus, nFixed); 


    BuildObjective(vbiqVarStatus, nFixed);

}


void BiqModel::BuildObjective(std::vector<BiqVarStatus> vbiqVarStatus, int nFixed)
{
    int i, j;
    double dConstant;
    //
    GetSubMatrix(vbiqVarStatus, nFixed);
    dConstant = GetConstant(Qs_, vbiqVarStatus);
    //
    GetLinear(Qs_, vbiqVarStatus, nFixed);

    // Set constant
    //Q_sub_[nVar_sub_ + nVar_sub_*(nVar_sub_ + 1)] = Q_[nVar_ + nVar_*(nVar_ + 1)]; // original 
    Q_sub_[nVar_sub_ + nVar_sub_*(nVar_sub_ + 1)] = dConstant; //  use GetConstant since the sol is shifted to {-1, 1}?
    
    // build Q_sub_
    for(i = 0; i < nVar_sub_; ++i)
    {
        for(j = 0; j < nVar_sub_; ++j)
        {
            Q_sub_[i + j * (nVar_sub_ + 1)] = pdTmp_sub_[i + j * nVar_sub_]; 
            /*
            // If we start with zeros on diag it should remain that way
            if(i != j)
            {
                Q_sub_[i + j * (nVar_sub_ + 1)] = pdTmp_sub_[i + j * nVar_sub_];
            }
            else
            {
                Q_sub_[i + j * (nVar_sub_ + 1)] = 0.0; 
            }
            */ 
        }
    }

    // lay down the linear terms
    for(i = 0; i < nVar_sub_; ++i)
    {
        Q_sub_[i + nVar_sub_*(nVar_sub_ + 1)] = pdTmpLinear_[i];
        Q_sub_[nVar_sub_ + i*(nVar_sub_ + 1)] = pdTmpLinear_[i];
    }
}

/// @brief 
/// @param RHSsource 
/// @param sMatSource 
/// @param RHSdest 
/// @param sMatdest 
void BiqModel::BuildConstraints(int nRows,
                                double *RHSsource, std::vector<Sparse> sMatArraySource,
                                double *RHSdest,   std::vector<Sparse> &sMatArraydest,
                                std::vector<BiqVarStatus> vbiqVarStatus, int nFixed)
{
    int nnzAdded;
    int i;
    double dCoefMatNorm;
    double dScaleFactor;
    double dConstant;

    Sparse::iterator itSource;
    
    for(int k = 0; k < nRows; ++k)
    {
        dCoefMatNorm = GetSubMatrixSparse(sMatArraySource.at(k), vbiqVarStatus, nnzAdded, nFixed);
        dConstant = GetConstant(sMatArraySource.at(k), vbiqVarStatus);
        GetLinear(sMatArraySource.at(k), vbiqVarStatus, nFixed);

        //
        //dScaleFactor = 1.0;
        dScaleFactor = 1.0 / (1.0 + fabs(RHSsource[k]) + dCoefMatNorm); 
        RHSdest[k] =  RHSsource[k];
        // Quad Part
        for(i = 0; i < nnzAdded; ++i)
        {
            if(fabs(sTmp_sub_.at(i).dVal_) > 0.0)
            {
                sMatArraydest[k].push_back(BiqSparseTriple(sTmp_sub_.at(i).i_, sTmp_sub_.at(i).j_, sTmp_sub_.at(i).dVal_ * dScaleFactor));
            }
        } 

        // Linear Part
        for(int j = 0; j < nVar_sub_; ++j)
        {
            if(fabs(pdTmpLinear_[j]) > 0.0)
            {
                sMatArraydest[k].push_back(BiqSparseTriple(nVar_sub_, j, pdTmpLinear_[j] * dScaleFactor));
            }
        }

        // new constant value on left-side
        if(fabs(dConstant) > 0.0)
        {
            sMatArraydest[k].push_back(BiqSparseTriple(nVar_sub_, nVar_sub_, dConstant * dScaleFactor));
        }

        RHSdest[k] *= dScaleFactor;
    }
}

/// @brief 
/// @param vbiqVarStatus 
/// @param nFixed 
void BiqModel::GetSubMatrix(std::vector<BiqVarStatus> vbiqVarStatus, int nFixed)
{
    int i, ii, jj;

    // init the submatrix zero out Q_sub_
    for(i = 0; i < nVar_sub_*nVar_sub_; ++i)
    {
        pdTmp_sub_[i] = 0.0;
    }

    for(auto & it: Qs_)
    {
        
        if(vbiqVarStatus[it.i_] == BiqVarFree && vbiqVarStatus[it.j_] == BiqVarFree && it.i_ < nVar_ && it.j_ <nVar_)
        {
            ii = it.i_ - vOffset_.at(it.i_);
            jj = it.j_ - vOffset_.at(it.j_);
            pdTmp_sub_[ii + jj * nVar_sub_] = it.dVal_;
            pdTmp_sub_[jj + ii * nVar_sub_] = it.dVal_;
        }
    }
}

/// @brief 
/// @param sSourceMat As_ or Bs_
/// @param vbiqVarStatus this is until the BiqNodeDesc
///         being able to pass the real :
/// @return 
double BiqModel::GetSubMatrixSparse(Sparse sSourceMat, 
                                    std::vector<BiqVarStatus> vbiqVarStatus, 
                                    int &nnzAdded, int nFixed)
{
    double dRetNorm = 0.0;
    double dVal;
    int ii, jj = 0;
    bool bBothFree;
    
    nnzAdded = 0;
    // fill in sTmp_sub_
    Sparse::iterator itDest =  sTmp_sub_.begin(); // TODO how to not be a pointer 
    for(auto& itSource : sSourceMat)
    {
        bBothFree =  vbiqVarStatus[itSource.i_] == BiqVarFree && vbiqVarStatus[itSource.j_] == BiqVarFree;
        if(bBothFree && itSource.i_ < nVar_ && itSource.j_ < nVar_)
        {
            ii = itSource.i_ - vOffset_.at(itSource.i_);
            jj = itSource.j_ - vOffset_.at(itSource.j_);
            dVal = itSource.dVal_;
            if(jj <= ii)
            {
                itDest->i_ = ii;
                itDest->j_ = jj;
            }
            else
            {
                itDest->i_ = jj;
                itDest->j_ = ii;
            } 
            itDest->dVal_ = dVal;
            dRetNorm += dVal*dVal;
            if(ii != jj)
            {
                dRetNorm += dVal*dVal;
            }
            // incress the itDest
            ++itDest;
            ++nnzAdded;
        }
    } 
    dRetNorm = sqrt(dRetNorm);

    return dRetNorm;
}

/// @brief loop over sMat and add up {-1, 1} const value from {0,1} sol
/// @param sMat 
/// @param solution 
/// @param piFixed 
/// @return 
double BiqModel::GetConstant(Sparse &sMat, std::vector<BiqVarStatus> vbiqVarStatus)
{
    double dRetConst = 0.0;
    double dTmp;
    double dSolVal_i;
    double dSolVal_j;
    for(auto &it : sMat)
    {

        // if in the quad part
        if(it.i_ < nVar_ && it.j_ < nVar_)
        {
            dSolVal_i = static_cast<double>(vbiqVarStatus.at(it.i_));
            dSolVal_j = static_cast<double>(vbiqVarStatus.at(it.j_));
            if(vbiqVarStatus.at(it.i_) != BiqVarFree && vbiqVarStatus.at(it.j_) != BiqVarFree)
            {

                dTmp = (2*dSolVal_i - 1.0)*(2*dSolVal_j - 1.0)*it.dVal_;
                dRetConst += dTmp;
                if(it.i_ != it.j_)
                {
                    dRetConst += dTmp;
                }
            }
        }
        // else if the corner const val
        else if(it.i_ == nVar_ && it.j_ == nVar_)
        {
            dRetConst += it.dVal_;
        }
        // else if in the linear column
        else if(it.j_ == nVar_)
        {
            if(vbiqVarStatus[it.i_] != BiqVarFree)
            {
                dSolVal_i = static_cast<double>(vbiqVarStatus.at(it.i_));
                dRetConst += (2*dSolVal_i- 1.0)*2*it.dVal_;
            }
        }
        // else if in the linear row
        else if(it.i_ == nVar_)
        {
            if(vbiqVarStatus[it.j_] != BiqVarFree)
            {
                dSolVal_j = static_cast<double>(vbiqVarStatus.at(it.j_));
                dRetConst += (2*dSolVal_j- 1.0)*2*it.dVal_;
            }
        }
    }

    return dRetConst;
}

void BiqModel::GetLinear(Sparse &sSource, std::vector<BiqVarStatus> vbiqVarStatus, int nFixed)
{
    double dVal = 0.0;
    double dSolVal_i;
    double dSolVal_j;
    int ii, jj = 0;
    bool bIfree, bJfree;
    
    for(int i = 0; i < nVar_sub_; ++i)
    {
        pdTmpLinear_[i] = 0.0;
    }

    for(auto &itSource : sSource)
    {
        ii = itSource.i_;
        jj = itSource.j_;
        dVal = itSource.dVal_;
        

        bIfree = vbiqVarStatus[ii] == BiqVarFree;
        bJfree = vbiqVarStatus[jj] == BiqVarFree;  
        // if we are already linear 
        if(ii < nVar_ && bIfree && jj == nVar_)
        {   
            pdTmpLinear_[ii - vOffset_.at(ii)] += dVal;
        }
        else if(jj < nVar_ && bJfree && ii == nVar_)
        {
            pdTmpLinear_[jj - vOffset_.at(jj)] += dVal;
        }
        // if in quad part
        else if(ii < nVar_ && jj < nVar_)
        {
            if(bIfree && !bJfree)
            {
                dSolVal_j = static_cast<double>(vbiqVarStatus.at(jj));
                pdTmpLinear_[ii - vOffset_.at(ii)] += (2 * dSolVal_j - 1.0)*dVal;
            }
            else if(!bIfree && bJfree)
            {
                dSolVal_i = static_cast<double>(vbiqVarStatus.at(ii));
                pdTmpLinear_[jj - vOffset_.at(jj)] += (2 * dSolVal_i - 1.0)*dVal;
            }
        }

    }

}

int BiqModel::GetOffset(std::vector<BiqVarStatus> vbiqVarStatus)
{
    int nFixedRet = 0;
    
    for(size_t i = 0; i < vOffset_.size(); ++i)
    {
        if(vbiqVarStatus.at(i) != BiqVarFree)
        {
            vOffset_.at(i) = 0;
            ++nFixedRet;  
        }
        else
        {
            vOffset_.at(i) = nFixedRet;
        }
    }

    return nFixedRet;
}

double BiqModel::UpdateInequalities(int &nAdded, int &nSubtracted)
{
    double dRetVal = 0.0;
    int yIndex;
    nSubtracted = 0;
    nAdded = 0;
    BiqTriInequality btiTemp;
    BiqTriTuple bttTemp;
    TriMap::iterator itMap;
    TriCuts::iterator itIneq;
    TriCuts::iterator itNextIneq;
    Heap_ = TriHeap(std::less<BiqTriInequality>(), std::move(container));

    

    // Fill the Heap_ with most violated cuts
    dRetVal = GetViolatedCuts();

    // remove cuts from Cuts_ array
    itNextIneq = Cuts_.begin();
    yIndex = mB_ + mA_;
    int count = 0 ;
    for(itIneq = Cuts_.begin(); itIneq < Cuts_.begin()+nIneq_; ++itIneq) 
    {
        itIneq->value_ = EvalInequalities(itIneq->type_, itIneq->i_, itIneq->j_, itIneq->k_);
        
        itIneq->y_ = y_[yIndex]; yIndex++; count++;
        // if multiplier is small and inequality is large
        if(itIneq->y_ < 1e-8 && itIneq->value_ > gradInorm_/10.0)
        {
            ++nSubtracted;

            bttTemp =  BiqTriTuple(itIneq->type_,
                                itIneq->i_, 
                                itIneq->j_, 
                                itIneq->k_);

            itMap = Map_.find(bttTemp);
            itMap->second = false;
            //std::printf("-  ");
        }
        else
        {
            itNextIneq->i_ = itIneq->i_;
            itNextIneq->j_ = itIneq->j_;
            itNextIneq->k_ = itIneq->k_;
            itNextIneq->type_ = itIneq->type_;
            itNextIneq->value_ = itIneq->value_;
            itNextIneq->y_ = itIneq->y_;
            ++itNextIneq;
             //std::printf("   ");
        }
    /*
        std::printf("%4d  t: %4d i: %4d j: %4d k: %4d \t value: %8f \t y: %8f\n",
                            count, itIneq->type_, itIneq->i_, itIneq->j_, itIneq->k_,
                            itIneq->value_, itIneq->y_);
    */
    }

    /*
    count = 0;
    for(auto it = Cuts_.begin(); it < Cuts_.begin() + nIneq_; ++it)
    {
        ++count;
        std::printf("%4d  t: %4d i: %4d j: %4d k: %4d \t value: %8f \t y: %8f\n",
                            count, it->type_, it->i_, it->j_, it->k_,
                            it->value_, it->y_);
    }
    */
    while (!Heap_.empty())
    {
        // get least violated cut
        btiTemp = Heap_.top();
        Heap_.pop();
        if(itNextIneq != Cuts_.end())
        {
            // update the map so we know this cut has been included
            bttTemp =  BiqTriTuple(btiTemp.type_, 
                                btiTemp.i_, 
                                btiTemp.j_, 
                                btiTemp.k_
            );
            /*
            std::printf("type: %d\t i: %d\t j: %d\t k: %d\t val: %f\n",
                        btiTemp.type_, btiTemp.i_, btiTemp.j_, btiTemp.k_,btiTemp.y_);
            */
                        
            itMap = Map_.find(bttTemp);
            if(itMap == Map_.end())
            {
                Map_.insert({bttTemp, true});
                // add it to the cut array
                itNextIneq->i_ = btiTemp.i_;
                itNextIneq->j_ = btiTemp.j_;
                itNextIneq->k_ = btiTemp.k_;
                itNextIneq->type_ = btiTemp.type_;
                itNextIneq->value_ = btiTemp.value_;
                itNextIneq->y_ = btiTemp.y_;
                ++itNextIneq;
                ++nAdded;
            }
            else if(itMap->second == false)
            {
                // back in the cuts array
                itMap->second = true;
                // add it to the cut array
                itNextIneq->i_ = btiTemp.i_;
                itNextIneq->j_ = btiTemp.j_;
                itNextIneq->k_ = btiTemp.k_;
                itNextIneq->type_ = btiTemp.type_;
                itNextIneq->value_ = btiTemp.value_;
                itNextIneq->y_ = btiTemp.y_;
                ++itNextIneq;
                ++nAdded;
            }
            else
            {
                /* Do not add the cut it is already in Cuts_*/
            }
        }
        
    }
    
    // TODO need to update the dual variable on each cut
    // update the number of inequalites
    nIneq_ = nIneq_ - nSubtracted + nAdded;
    return dRetVal;
}

/// @brief 
/// @return 
double BiqModel::GetViolatedCuts()
{ 
    double dRetMinIneq = INFINITY;
    double dTestIneq;
    int i, j, k, type = 0;
    double y;
    // params
    const double dGapCuts = BiqPar_->entry(BiqParams::dGapCuts);
    const int nCuts = BiqPar_->entry(BiqParams::nCuts);
    BiqTriInequality btiTemp;
    BiqTriTuple bttTemp;
    TriMap::iterator itMap;
    // loop through all the inequalities
    for(type = 1; type <= 4; ++type)
    {
        for(i = 0; i < (nVar_sub_ + 1); ++i)
        {
            for(j = 0; j < i; ++j)
            {
                for(k = 0; k < j; ++k)
                {
                    bttTemp = BiqTriTuple(type,i,j,k);
                    itMap = Map_.find(bttTemp);

                    // check if the key is in the map
                    if(itMap != Map_.end() && itMap->second)
                    {
                        continue;
                    }
                    
      
                    dTestIneq = EvalInequalities(static_cast<TriType>(type), i, j, k);
                    // keep track of the minimum ineq
                    dRetMinIneq = (dTestIneq < dRetMinIneq) ? dTestIneq : dRetMinIneq;

                    y = 0.0; // initialize the dual variable to be zero for a new cut
                    btiTemp = BiqTriInequality(
                                    static_cast<TriType>(type),
                                    i,j,k,
                                    dTestIneq,
                                    y);

                    if(static_cast<int>(Heap_.size()) < nCuts)
                    {
                        // (1) add cuts until full
                        // check if the cut has been added .. to the map
                        if(dTestIneq < dGapCuts)
                        {
                            Heap_.push(btiTemp);
                            /*
                            std::printf("%d \t %d \t %d \t %d \t %20.16f\n",
                                 type, i, j, k, dTestIneq);
                            */
                                 
                        }
                            

                    }
                    else if(btiTemp < Heap_.top())
                    {
                        // (2) the list is full, so remove the least violated 
                        // and then add in this new cut!
                        Heap_.pop();
                        Heap_.push(btiTemp);
                    }
                }
            }
        }
    }

    return dRetMinIneq;
}

inline double BiqModel::EvalInequalities(TriType triType, int ii, int jj, int kk)
{
    double dRetIneqVal = 0.0;

    // compute the inequality
    switch (triType)
    {
        case ONE:
        {
            dRetIneqVal = X_[ii+jj*(nVar_sub_ + 1)] + X_[ii+kk*(nVar_sub_ + 1)] + X_[jj+kk*(nVar_sub_ + 1)] + 1.0;
            //std::printf("%d, %d, %d \t dRetIneqVal = %f = %f + %f + %f + 1.0\n",ii, jj, kk, dRetIneqVal,X_[ii+jj*(nVar_sub_ + 1)] , X_[ii+kk*(nVar_sub_ + 1)] , X_[jj+kk*(nVar_sub_ + 1)]);
        }
        break;
        
        case TWO:
        {
            dRetIneqVal = X_[ii+jj*(nVar_sub_ + 1)] - X_[ii+kk*(nVar_sub_ + 1)] - X_[jj+kk*(nVar_sub_ + 1)] + 1.0;
        }
        break; 

        case THREE:
        {
            dRetIneqVal = -X_[ii+jj*(nVar_sub_ + 1)] + X_[ii+kk*(nVar_sub_ + 1)] - X_[jj+kk*(nVar_sub_ + 1)] + 1.0;
        }
        break;

        case FOUR:
        {
            dRetIneqVal = -X_[ii+jj*(nVar_sub_ + 1)] - X_[ii+kk*(nVar_sub_ + 1)] + X_[jj+kk*(nVar_sub_ + 1)] + 1.0;
        }
        break;

        default: // error
            break;
    }



    return dRetIneqVal;
}

double BiqModel::EvalSolution(std::vector<int> & solution)
{
    double dRetSol = 0.0;
    int i, j;
    double dVal;
    
    for(auto &it: Qs_)
    {
        i = it.i_;
        j = it.j_;
        dVal = it.dVal_;
        
        //printf("solution[%d] = %d \t solution[%d] = %d\n", i, solution.at(i), j, solution.at(j));
        if(i < nVar_ && j < nVar_)
        {
            dRetSol += 2*dVal*(2*solution.at(i) - 1.0)*(2*solution.at(j) - 1.0);
        } 
        else if(i == nVar_ && j < nVar_)
        {
            dRetSol += 2*dVal*(2*solution.at(j) - 1.0);
        }
        else if(i < nVar_ && j == nVar_)
        {
            dRetSol += 2*dVal*(2*solution.at(i) - 1.0);
        }
        else if(i == nVar_ && j == nVar_)
        {
            dRetSol += dVal;
        }
    }

    if(!max_problem_)
    {
        dRetSol = -dRetSol;
    }

    return dRetSol;
}


bool BiqModel::isFeasibleSolution(std::vector<int> solution)
{

    int i;
    double dSum, dTmp;
    bool bRet = true;
    // convert to {1, -1}
    for(auto &it: solution)
    {
        it = 2*it - 1;
    } 

    solution.at(nVar_) = 1; 

    // loop over equality cons of master problem;
    for(int eqCon = 0; eqCon < mB_original_; ++eqCon)
    {
        dSum  = 0;
        for(auto &it : Bs_.at(eqCon))
        {
            dTmp = it.dVal_ * solution.at(it.i_) * solution.at(it.j_);
            dSum += dTmp;
            if(it.i_ != it.j_)
            {
                dSum += dTmp;
            }
        }

        // check that the sum matches the rhs
        if(dSum != b_[eqCon])
        {
            //std::printf("dSum = %f \t b_[%d] = %f\n", dSum, i, b_[i]);
            bRet = false;
            break;
        }
    }

    // loop over inequality cons of master problem;
    i = 0;
    for(auto &itCons : As_)
    {
        dSum = 0.0;
        for(auto &it : itCons)
        {
            dTmp = it.dVal_ * solution.at(it.i_) * solution.at(it.j_);
            dSum += dTmp;
            if(it.i_ != it.j_)
            {
                dSum += dTmp;
            }
        }

        // check that the sum is grater than rhs
        if(dSum > a_[i])
        {
            //std::printf("dSum = %f \t a_[%d] = %f\n", dSum, i, a_[i]);
            bRet = false;
            break;
        }
        ++i; 
    }

    return bRet;
}

double BiqModel::primalHeuristic()
{

    assert(viSolution_1_.size() == viSolution_2_.size());

    const int nRuns = BiqPar_->entry(BiqParams::nPrimalRuns);

    double dRet = 0.0;
    double gamma;
    double dRand;
    double dTempEval;

    bool bStop = false;

    for(int i = 0; i < nRuns; ++i)
    {
        if(bStop) break;
        for(gamma = 0.0; gamma < 1.0; gamma += 0.01)
        {
            if(bStop) break;
            for(size_t j = 0; j < viSolution_1_.size(); ++j)
            {
                dRand = (static_cast<double>(rand())/ 2147483647.0);
                if(dRand > gamma)
                {
                    viSolution_1_.at(j) = 1;
                    viSolution_2_.at(j) = 0;
                }
                else
                {
                    viSolution_1_.at(j) = 0;
                    viSolution_2_.at(j) = 1;
                }
            }
            if(isFeasibleSolution(viSolution_1_))
            {
                //std::printf("BiqModel::primalHeuristic() Hit one!!\n");
                dTempEval = EvalSolution(viSolution_1_);
                UpdateSol(dTempEval, viSolution_1_);
                bStop = true;
            }
            if(!bStop && isFeasibleSolution(viSolution_2_))
            {
                //std::printf("BiqModel::primalHeuristic() Hit one!!\n");
                dTempEval = EvalSolution(viSolution_2_);
                UpdateSol(dTempEval, viSolution_2_);
                bStop = true;
            }
        }
    }
    
    return dRet;
}


void BiqModel::KCheuristic(std::vector<BiqVarStatus> vbiqVarStatus)
{
    // A Gready Heuristic used to find
    // Solutions to K Cluster Problems

    double dTempEval;

    double MaxFracSol;
    int posFix; 

    int K_Clusters = b_[0]/2;
    int nFixedToOne = 0;
    int K_diff;


    // get fractional solution for fixing
    std::vector<double> fracSol = GetFractionalSolution(vbiqVarStatus);

    // initilize a solution vector
    for(auto i = 0; i < nVar_; ++i)
    {
        if(vbiqVarStatus.at(i) == BiqVarFixedToOne)
        {
            viSolution_1_.at(i) = 1;
            ++nFixedToOne;
        }
        else
        {
            viSolution_1_.at(i) = 0;
        }
    }
    viSolution_1_.at(nVar_) = 0; 
    K_diff = K_Clusters - nFixedToOne;

    // for each remaining unfixed
    for(int k = 0; k < K_diff; ++k)
    {
        MaxFracSol = -1.0;
        posFix = -1;
        for(int i = 0; i < nVar_; ++i)
        {
            if( vbiqVarStatus.at(i) == BiqVarFree &&
                viSolution_1_.at(i) == 0
                )
                //printf("fracSol = %f  status = %d    solVal = %d\n", fracSol.at(i),vbiqVarStatus.at(i), viSolution_1_.at(i));
            if(
                vbiqVarStatus.at(i) == BiqVarFree &&
                fracSol.at(i) > MaxFracSol &&
                viSolution_1_.at(i) == 0
               )
               {
                    posFix = i;
                    MaxFracSol = fracSol.at(i);
               }
        }
        //printf("MaxFracSol = %f\n", MaxFracSol);
        assert(posFix >= 0);
        viSolution_1_.at(posFix) = 1;
    }
    
    
    //
#ifdef DEBUG
   int s = 0;
    for(auto &it: viSolution_1_)
    {
        s += it;
    }
    //printf("s = %d\n",s);
    assert(s == K_Clusters);
    assert(isFeasibleSolution(viSolution_1_));
#endif


    // should assert the we have k fixed
    if(isFeasibleSolution(viSolution_1_))
    {
        dTempEval = EvalSolution(viSolution_1_);

        //=if(dTempEval >= 255) printf("KC dTempEval = %f\n", dTempEval);
        UpdateSol(dTempEval, viSolution_1_);
    }

}

void BiqModel::NieveUpdateHeuristic(std::vector<BiqVarStatus> vbiqVarStatus)
{
    double dTempEval;

    // Bail if there is no known solutions
    if(broker()->hasKnowledge(AlpsKnowledgeTypeSolution) == false) return;
    // from the broker get the data stored from update solution
    std::pair< AlpsKnowledge *, double > bestSolPair = broker()->getBestKnowledge(AlpsKnowledgeTypeSolution);
    // cast the AlpsKnowledge pointer to a BiqSolution pointer
    BiqSolution* pBestSol = static_cast<BiqSolution*>(bestSolPair.first);
    
    viSolution_1_ = pBestSol->getSol();
    //printf("viSolution_1_.size():  %d     vbiqVarStatus.size():  %d\n", viSolution_1_, vbiqVarStatus);
    //assert(viSolution_1_.size() == vbiqVarStatus.size());
    // fix variables according to var status
    for(size_t i = 0; i < nVar_; ++i)
    {
        if(vbiqVarStatus.at(i) == BiqVarFixedToOne)
        {
            viSolution_1_.at(i) = 1;
        } 
        else if(vbiqVarStatus.at(i) == BiqVarFixedToZero)
        {
            viSolution_1_.at(i) = 0;
        }
        else
        {
            // Do Nothing
        }
    }
    // if feasible try to update solution

    // should assert the we have k fixed
    if(isFeasibleSolution(viSolution_1_))
    {
        dTempEval = EvalSolution(viSolution_1_);

        printf("NieveUpdateHeuristic dTempEval = %f\n", dTempEval);
        if(UpdateSol(dTempEval, viSolution_1_)) printf("::NieveUpdateHeuristic got a best!!!");
    }

    return;
}

double BiqModel::GWheuristic(int nPlanes, std::vector<BiqVarStatus> vbiqVarStatus)
{
    // update... lets keep track of the best solution
    // and only at the end attempt to add it....
    // if no cons then always feasible might save time!
    // Check the isfeasible somthing is wrong.. 

    double dRetBest = -1e+9;
    double dPlaneNorm;
    double dPlaneInvNorm;
    double sca;
    double dTempEval;
    double dBestVal;

    int index;
    int subN = nVar_sub_+1;

    bool bHasBest;
    const bool bOneOptSearch = BiqPar_->entry(BiqParams::bOneOptSearch);

    std::vector<double> hyperPlane;

    hyperPlane.resize(M_);

    viSolution_1_.at(nVar_) = 1;
    viSolution_2_.at(nVar_) = 1;
        
    dBestVal = getBestVal();
    bHasBest = hasBestVal();

    for(int k = 0; k < nPlanes; ++k)
    {
        // create a random hyperplane
        dPlaneNorm = 0.0;
        for(auto &it: hyperPlane)
        {
            it = 1.0 + (int) 100.0 * rand() / ((double) RAND_MAX + 1);
            //it = 1.0 + 100.0 * static_cast<double>(rand())/(static_cast<double>(RAND_MAX) + 1.0);
            dPlaneNorm += it*it;
        }
        dPlaneNorm = sqrt(dPlaneNorm);
        dPlaneInvNorm = 1.0/dPlaneNorm;

        for(auto &it: hyperPlane)
        {
            it *= dPlaneInvNorm;
        }

       // Compute cuts temp_x1 and temp_x2 generated by the hyperplane
       index = 0;
        for(int i = 0; i < nVar_; ++i)
        {
            switch (vbiqVarStatus.at(i))
            {
            case BiqVarFixedToOne:
            {
                viSolution_1_.at(i) = 1;
                viSolution_2_.at(i) = 1;
            }break;

            case BiqVarFixedToZero:
            {
                viSolution_1_.at(i) = 0;
                viSolution_2_.at(i) = 0;
            }break;

            case BiqVarFree:
            {
                sca = 0.0;
                for(int j = 0; j < M_; ++j)
                {
                    sca += hyperPlane.at(j)*Z_[j * subN + index];
                    /*
                    std::printf("%f \t %f \t %f \t %d \t %d \t %d\n", 
                                sca, hyperPlane.at(j),  Z_[j * subN + index],
                                j , subN , index);
                                */
                }
                if(sca < 0)
                {
                    viSolution_1_.at(i) = 0;
                    viSolution_2_.at(i) = 1;
                }
                else
                {
                    viSolution_1_.at(i) = 1;
                    viSolution_2_.at(i) = 0;  
                }

                ++index;
            }break;


            default:
                break;
            }

            

        }

        // now we have solutions see how they compare to the best ..  
        if(isFeasibleSolution(viSolution_1_))
        {
            dTempEval = EvalSolution(viSolution_1_);
            dBestVal = dTempEval;
            //printf("GW dTempEval = %f\n", dTempEval);
            if(UpdateSol(dTempEval, viSolution_1_) && bOneOptSearch)
            {
                OneOptLocalSearch(viSolution_1_);
            }
            
        }

        if(isFeasibleSolution(viSolution_2_))
        {
            dTempEval = EvalSolution(viSolution_2_);
            dBestVal = dTempEval;
            //printf("GW dTempEval = %f\n", dTempEval);
            if(UpdateSol(dTempEval, viSolution_2_) && bOneOptSearch)
            {
                OneOptLocalSearch(viSolution_2_);
            }
            
        }

      

    }

    //std::printf("best ... %f\n", dBestVal);


    return dRetBest;
}

bool BiqModel::UpdateSol(double dVal, std::vector<int> solution)
{

    double bestVal = 0;
    bool bForceAdd = false;
    bool bRet = false; 
    if(hasBestVal())
    {
        bestVal = getBestVal();
    }
    else
    {
        bForceAdd = true;
    }
    
    
    //printf("Val => %f  Beta => %f test => %d\n",dVal, bestVal, dVal < bestVal);
    BiqSolution* biqSol = new BiqSolution( this, solution, dVal);
    if(max_problem_ && (dVal > bestVal || bForceAdd))
    {
        
        broker()->addKnowledge(AlpsKnowledgeTypeSolution, biqSol, -dVal);
        bRet = true;
        //printf("Beta updated => %f\n",-dVal);
    }
    else if(!max_problem_ && (dVal < bestVal || bForceAdd))
    {     
        
        broker()->addKnowledge(AlpsKnowledgeTypeSolution, biqSol, dVal);
        bRet = true;
        //printf("Val => %f  Beta => %f Force? => %d \n",dVal, bestVal, bForceAdd);
    }   
    else
    {
        /* do nothing*/
    }

    return bRet;
}

std::vector<double> BiqModel::GetFractionalSolution(std::vector<BiqVarStatus> vbiqVarStatus)
{
    int ii, jj;
    double Xij;

    // Get Last row of X_ 
    for(int i = 0; i < nVar_; ++i)
    {
        if(vbiqVarStatus.at(i) == BiqVarFree)
        {
            ii = nVar_sub_;
            jj = i - vOffset_.at(i);
            Xij = X_[ii + jj * (nVar_sub_ + 1)];
            vdFracSol_.at(i) = (Xij + 1)/2;
        }
        else
        {
            vdFracSol_.at(i) = static_cast<double>(vbiqVarStatus.at(i));
        }
    }

    return vdFracSol_;
}

// This method wraps the brokers best with the desc best
// This is used since there is a delay in MPI sharing knowlage and it helps with prune
int BiqModel::getBestVal()
{
    int retBestVal;
    int brokerBestVal = static_cast<int>(broker()->getIncumbentValue());
    int pastBestVal = getPastBest();

    bool bHasPastBest = hasPastBest();
    bool bHasBrokerBest = broker()->hasKnowledge(AlpsKnowledgeTypeSolution);

    if(max_problem_)
    {
        brokerBestVal = -brokerBestVal;
    }
    

    

    if(max_problem_ && bHasPastBest && bHasBrokerBest && pastBestVal < brokerBestVal)
    {
        
        retBestVal = pastBestVal;
    } 
    else if(bHasPastBest && !bHasBrokerBest)
    {
        retBestVal = pastBestVal;
    }
    else
    {
        retBestVal = brokerBestVal;
    }

    

    //printf("bHasPastBest: %d  bHasBrokerBest: %d | brokerBestVal: %d   <>   pastBestVal: %d   <>   retBestVal: %d\n",bHasPastBest, bHasBrokerBest, brokerBestVal ,pastBestVal, retBestVal);

    return retBestVal;
}

// do we have a best??
// This method wraps brokers best with desc best. This helps with MPI
bool  BiqModel::hasBestVal()
{
    bool retHasBest = hasPastBest() || broker()->hasKnowledge(AlpsKnowledgeTypeSolution);
    return retHasBest;
}

bool BiqModel::pruneTest(double dBound)
{
    if(hasBestVal() == false)
    {
        // try primal 
        primalHeuristic();
        /*
        std::vector<int> solution;
        solution.resize(getNVar());
        BiqSolution* biqSol = new BiqSolution( this, solution, -1074);
        broker()->addKnowledge(AlpsKnowledgeTypeSolution, biqSol, -1074);
        */
        return false;
    }

    bool bRet = false;
    int bestVal = getBestVal();
    
    if((max_problem_ && static_cast<int>(floor(dBound)) <= bestVal) || 
       (!max_problem_ && static_cast<int>(ceil(dBound))  >= bestVal))
    {
        //std::printf("%f %d prune",dBound, bestVal);
        bRet = true;
    }
    //std::printf("bRet = %d \n",bRet);
    
    return bRet;
}

void BiqModel::addProvidedSol()
{
    // add the soluttion that is provided by the .par file
    double dSolutionValue = BiqPar()->entry(BiqParams::dSolutionValue);
    for(auto &it: viSolution_1_)
    {
        it = 0;
    }

    BiqSolution* biqSol = new BiqSolution( this, viSolution_1_, dSolutionValue);
    broker()->addKnowledge(AlpsKnowledgeTypeSolution, biqSol, dSolutionValue);

}


void BiqModel::GreedyUQBO()
{

    const bool  bSolutionProvided = BiqPar()->entry(BiqParams::bSolutionProvided);
    const double dSolutionValue = BiqPar()->entry(BiqParams::dSolutionValue);

    double dGreedySolVal;
    double dGreedyGain;
    double dGainTmp;
    double dGainSumTmp;
    double dBestGain;

    int indexBest;
    int bestVal;
    int tmpPos;
    int tmpIndex;
    


    GreedyGainsHeap ggHeap;
    GreedyGainsMap ggMap;
    GreedyGains ggBest = GreedyGains(); 
    
    std::vector<int> IndexAdded; 

    // print out
    //printf("Randomized Greedy Algo: \n");


        viSolution_1_.at(nVar_) = 1;
        // fill the compliment with 0 .. nVar_ 
        for(int i = 0; i < nVar_; ++i)
        {               
            ggMap.insert({i, GreedyGains(i,0,0)});
        }

        // pick a random first var
        indexBest = rand() % nVar_;
        bestVal =  rand() % 2;
        // place best in solution
        viSolution_1_.at(indexBest) = bestVal;
        IndexAdded.push_back(indexBest);
        // remove best from the map
        ggMap.erase(indexBest);

        for(int i = 0; i < nVar_ - 1; ++i)
        {

            dBestGain =-MAXFLOAT;
            indexBest = -1;
            bestVal = -1;
            // compute the gains 
            for(auto &it: ggMap)
            {   
                it.second.gainOne_ = 0;
                it.second.gainZero_ = 0;
                for(auto &j: IndexAdded)
                {
                    
                    it.second.gainOne_ += (2*viSolution_1_[j] - 1)*Q_[it.second.index_ + j*N_];
                    it.second.gainZero_ += (2*viSolution_1_[j] - 1)*Q_[it.second.index_ + j*N_];
                }

                it.second.gainOne_ *= 2;
                it.second.gainZero_ *= -2;


                it.second.gainOne_ +=  Q_[it.second.index_ + it.second.index_ *N_] + 2*Q_[it.second.index_ + nVar_*N_];
                it.second.gainZero_ +=  Q_[it.second.index_ + it.second.index_ *N_] - 2*Q_[it.second.index_ + nVar_*N_];
                //printf("Gain 0: %f \t Gain 1: %f \n",it.second.gainOne_,it.second.gainZero_);
                if(it.second.gainOne_ > dBestGain)
                {
                    dBestGain = it.second.gainOne_;
                    indexBest = it.second.index_;
                    bestVal = 1;
                }
                else if(it.second.gainZero_ > dBestGain)
                {
                    dBestGain = it.second.gainZero_;
                    indexBest = it.second.index_;
                    bestVal = 0;
                }
                else
                {
                    // do nothing ..
                }
            }

            assert(dBestGain != -INFINITY || dBestGain != INFINITY);
            assert(indexBest != -1);
            assert(bestVal != -1);
            
            //printf("Index: %d \t Val: %d \t Gain: %f \n",indexBest,bestVal,dBestGain);
            //(bestVal == 1) ? dGreedyGain = ggBest.gainOne_ : dGreedyGain = ggBest.gainZero_;
            //printf("index: %d \t %d Val: %f\n",indexBest,bestVal,dGreedyGain);
        
            // place best in solution
            viSolution_1_.at(indexBest) = bestVal;
            IndexAdded.push_back(indexBest);
            // remove best from the map
            ggMap.erase(indexBest);
        }

        
        dGreedySolVal = EvalSolution(viSolution_1_);
        //printf("nVar: %d Greed Heurisstic found solution with value: %f\t",nVar_,dGreedySolVal);
        //for(auto &it: viSolution_1_) printf("%d ", it);
    
        if(bSolutionProvided)
        {
            printf("Best Known: %f \t Gap: %f\n",dSolutionValue, abs(dSolutionValue - dGreedySolVal)/abs(dSolutionValue));
        }

        if(UpdateSol(dGreedySolVal, viSolution_1_))
        {
            printf("Greedy found new best: %f\n",dGreedySolVal);
        }
        OneOptLocalSearch(viSolution_1_);

}


void BiqModel::OneOptLocalSearch(std::vector<int>& vbiqSol)
{

    std::vector<double> vGain; 

    double dTempSum;
    double dBestGain; 
    double dBegVal;
    double dEndVal;
    double dTempMult;

    int iBestIndex;
    int iPrevBest;
    
    bool bHasBetter;

    // Initialize gains vector
    vGain.resize(nVar_);

    iPrevBest = iBestIndex;
    bHasBetter = true;

    iBestIndex = -1;
    dBestGain = -MAXFLOAT;
    for(int k = 0; k < nVar_; ++k)
    {
        dTempSum = 0.0;

        // g_k = -2y_k y'Q_k + Qkk - 2yk
        for(int i = 0; i <= nVar_; ++ i)
        {
            dTempSum += (2*vbiqSol.at(i) - 1)*Q_[i + k*N_];
        }
        dTempMult =  (vbiqSol.at(k) == 1) ?  -4: 4;
        dTempSum *= dTempMult;
        dTempSum += 4*Q_[k + k*N_];

        vGain.at(k) = dTempSum;

        if(vGain.at(k) > dBestGain)
        {
            dBestGain = vGain.at(k);
            iBestIndex = k;
        }

    }

    // update the solution with the best gain
    if(dBestGain > 0)
    {
        bHasBetter = true;
        vbiqSol.at(iBestIndex) = (vbiqSol.at(iBestIndex) == 1) ? 0: 1;
    }
    else
    {
        bHasBetter = false;
    }


    // while we can still find a gain
    while (bHasBetter)
    {
        iPrevBest = iBestIndex;
        iBestIndex = -1;
        dBestGain = -MAXFLOAT;
        // update gains
        for(int k = 0; k < nVar_; ++k)
        {
            if(k == iPrevBest)
            {
                vGain.at(k) = 0;
            }
            else
            {
                vGain.at(k) -= 8*(2*vbiqSol.at(k) - 1)*(2*vbiqSol.at(iPrevBest) - 1)*Q_[k + iPrevBest*N_];

            }
            
            
            if(vGain.at(k) > dBestGain)
            {
                dBestGain = vGain.at(k);
                iBestIndex = k;
            }
        }

        if(dBestGain > 0)
        {
            bHasBetter = true;
            vbiqSol.at(iBestIndex) = (vbiqSol.at(iBestIndex) == 1) ? 0: 1;
            //printf("iBestIndex: %d   Best gain: %f   EndSol - BegSol: %f\n",iBestIndex, dBestGain, dEndVal - dBegVal);
        }
        else
        {
            bHasBetter = false;
            //printf("Local Min\n");
        }
          
    }

    // compute the starting sol value for comparison
    dEndVal = EvalSolution(vbiqSol);
    if(UpdateSol(dEndVal, vbiqSol))
    {
        printf("1-Opt found new best: %f\n",dEndVal);
    }
}

