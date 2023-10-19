#include "headers/BiqModel.h"


BiqModel::BiqModel(
                int nVar, double *Q, int max_problem,
                std::vector<Sparse> As, double *a,
                std::vector<Sparse> Bs, double *b
   ) 
    : nVar_(nVar),
      Q_(Q),
      max_problem_(max_problem),
      As_(As),
      Bs_(Bs)
{
    N_ = nVar_ + 1;
    ISUPPZ_ = new int[2*N_];
    X_ = new double[N_*N_];
    W_ = new double[N_];
    Z_ = new double[N_*N_];

    /* DSYEVR optimal workspace query */
    M_ = 0;
    int INFO = 0;
    double *DWORK_SIZE = new double[1];
    int *IWORK_SIZE = new int[1];
    char JOBZ = 'N';    // eigenvalues and eigenvectors
    char RANGE = 'A';
    char UPLO = 'U';
    double VL = 0.0;
    double VU = 0.0;
    int IL = 1;
    int IU = 1;
    double ABSTOL = 0.0;
    int SIZE_FLAG = -1;

    dsyevr_(
        JOBZ, RANGE, UPLO, N_, X_, N_, 
        VL, VU, IL, IU, ABSTOL, M_, 
        W_, Z_, N_, ISUPPZ_, DWORK_SIZE, SIZE_FLAG, 
        IWORK_SIZE, SIZE_FLAG, INFO);

    if(INFO == 0)
    {
        nDWORK_ = static_cast<int>(DWORK_SIZE[0]);
        DWORK_ = new double[nDWORK_];
        nIWORK_ = IWORK_SIZE[0];
        IWORK_ = new int[nIWORK_]; 
    }
    else
    {
        //TODO exit on error
    }

    //temp 
    mA_ = As_.size();
    mB_ = Bs_.size() + N_;
    a_ = new double[mA_];
    b_ = new double[mB_];
    a_sub_ = new double[mA_];
    b_sub_ = new double[mB_];
    Q_sub_ = new double[N_*N_];
    FillSparseMatrix(Qs_,Q_,N_);
    // Copy the b vector in to the first Bs_.size() entries
    std::copy(b, b + Bs_.size(), b_);
    /* L-BFGS-B Data */
    int nMax = mB_ + mA_ + MaxNineqAdded;
    int wa_length = (2 * mmax + 5) * nMax + 11 * mmax * mmax + 8 * mmax;
    int iwa_length = 3 * nMax;

    nIneq_ = 0;
    f_ = 0.0;
    Cuts_.resize(MaxNineqAdded);   // TODO use push or resize ondemand 
    List_.resize(500); // TODO replace with param.cuts
    g_ = new double[nMax];
    y_ = new double[nMax];
    binf_ = new double[nMax];
    bsup_ = new double[nMax];
    nbd_ = new int[nMax];

    wa_ = new double[wa_length];
    iwa_ = new int[iwa_length];
 
    AddDiagCons();
    AllocSubCons();

    // set some memory for the sTemp_sub_ matrix
    sTmp_sub_.resize(N_*N_);
    pdTmp_sub_ = new double[N_*N_];
    pdTmpLinear = new double[nVar_];
    delete[] IWORK_SIZE;
    delete[] DWORK_SIZE;

    CreateSubProblem();
}

BiqModel::~BiqModel()
{
    // TODO clean up and make sure this is safe
    delete[] ISUPPZ_; 
    delete[]  X_;
    delete[]  W_;
    delete[] Z_;
    delete[]  g_; 
    delete[]  y_ ;
    delete[]  binf_;
    delete[]  bsup_ ;
    delete[]  nbd_ ;
    delete[]  wa_ ;
    delete[] iwa_;
    delete[] DWORK_;
    delete[] IWORK_;
    delete[] a_;
    delete[] b_;
    delete[] b_sub_;
    //delete[] a_sub_;
    delete[] pdTmpLinear;
    delete[] pdTmp_sub_;
}

AlpsKnowledge * BiqModel::decode(AlpsEncoded & encode) const
{
    std::cerr << "Not implemented!" << std::endl;
    throw std::exception();
}

void BiqModel::AddDiagCons()
{
    // resize ?
    // then loop 

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
    }
}

/// @brief
/// Need to set the size of the array of cons
/// Then need to set size for con
void BiqModel::AllocSubCons()
{
    Bs_sub_.resize(100);
    As_sub_.resize(100);

}

/// @brief 
/// @return 
double BiqModel::SDPbound()
{
    int i;
    int iStatus;
    int nbit = 0;
    double dAlpha = 0.1;
    double dFixedVal;
    double dTol = 1e-6;
    int len_y = mA_ + mB_ + nIneq_;

    // get fixed value
    dFixedVal = GetFixedValue();

    // 

    for(i = 0; i < len_y; ++i)
    {
        y_[i] = 0.0;
    }
    
    sim(dAlpha);
    iStatus = CallLBFGSB(dAlpha, dFixedVal, dTol, nbit);

    return 0.0;
}

double BiqModel::GetFixedValue()
{
    // TODO

    return 0.0;
}

int BiqModel::CallLBFGSB(double dAlpha, double dFixedVal, double dTol, int &nbit)
{
    int retStatus = 1;
    int mem = mmax;
    int i;
    bool bStopBFGS = false;
    bool bPrune;
    double factr = 0.000005;
    double pgtol = 0.0;
    double dMinTemp;
    double dBound;
    char task[60];
    int iprint = 101;
    char csave[60];
    int lsave[4];
    int isave[44];
    double dsave[29];
    // compute the number of variables in BFGS
    int len_y = mA_ + mB_ + nIneq_;
    // set task = "START" and the rest of it to empty char
    strcpy(task, "START");
    for(i = 5; i < 60; ++i) task[i] = ' ';


    
    // Initialize y
    for(i = 0;  i < nIneq_; ++i)
    {
        y_[mB_ + mA_ + i] = Cuts_.at(i).y_;
    }

    // compute the bounds for y
    //// equalities
    for(i = 0; i < mB_; ++i)
    {
        nbd_[i] = 0; 
    }
    //// inequalities
    for(i = mB_; i < len_y; ++i)
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
            for(i = 0; i < mB_; ++i)
            {
                gradEnorm_ = (gradEnorm_ < fabs(g_[i])) ? fabs(g_[i]) : gradEnorm_;
            }
            gradEnorm_ /= scaleEq;

            // Compute Infinity-norm of gradI = min(g[mB_:nLBFGSB_Vars-1], 0.0)
            gradInorm_ = 0.0;
            for(i = mB_; i < len_y; ++i)
            {
                dMinTemp = (g_[i] > 0.0) ? 0.0 : g_[i];
                gradInorm_ = (gradInorm_ < fabs(dMinTemp)) ? fabs(dMinTemp) : gradInorm_;
            }
            gradInorm_ /= scaleIneq;

            dBound = (max_problem_) ? f_ - dFixedVal : dFixedVal - f_;
            bPrune = Prune(); 

            if(bPrune
               || (gradEnorm_ < dTol && gradInorm_ < dTol)
               || nbit >= MAXITER)
               {
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
    int N2 = N_*N_;
    int incx = 1; 
    dscal_(N2, dAlphaInv, X_, incx);
    // Update Z so that X = Z*Z'
    double scb = 10000.0 * dAlpha;
    double sca = 0.01/ sqrt(scb);
    int NM = N_ * M_;
    int incz = 1;
    dscal_(NM, sca, Z_, incz); 
    return retStatus;
}

/// @brief 
/// @param alpha 
void BiqModel::sim(double alpha)
{
    int N2 = N_*N_;
    double dAlphaInv = 1.0 / alpha;
    // prepare to copy Q to X_
    // memory spacing for the two arrays
    int INCQ = 1; 
    int INCX = 1;
    // copy Q to X 
    dcopy_(N2, Q_, INCQ, X_, INCX);

    // call the A and B methods transposed
    B(TRANSP, alpha);
    A(TRANSP, alpha);
    ProjSDP();

    f_ = f_ * 1000.0 * (dAlphaInv / 1000.0); // TODO see what happens 


    // call the A and B methods not transposed
    B(NOTRANSP, alpha);
    A(NOTRANSP, alpha);



    f_ += alpha * (0.5 * static_cast<double>(N2));

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
    int LDX = N_;

    // Operator norm for X to get a upper bound of the spectral radius of X
    // |X|_2 \leq \sqrt{ |X|_1 |X|_inf }  (Holder's inequality)
    //          = |X|_1 = |X|_inf  (since X is symmetric)
    //
    // Frobenius norm is also an upper bound on the spectral radius of X:
    //      |X|_2 <= |X|_F
    char NORM = 'I'; // in
    double norminf = dlansy_(NORM, UPLO, N_, X_, LDX, DWORK_);
    NORM = 'F';
    double normfro = dlansy_(NORM, UPLO, N_, X_, LDX, DWORK_);

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
    double ABSTOL = 0; // TODO could be a param 
    int LDZ = N_;
    int INFO;

    int M = 0;

    dsyevr_(JOBZ, RANGE, UPLO, N_, X_, LDX, VL, VU, IL, IU, ABSTOL,
            M, W_, Z_, LDZ, ISUPPZ_, DWORK_, nDWORK_, IWORK_, nIWORK_, INFO);

   
   // Check if the eigensolver failed (i.e., if INFO != 0)
    if (INFO) {
        VL = 0.0;
        ABSTOL = 0.0;
        dsyevr_(JOBZ, RANGE, UPLO, N_, X_, LDX, VL, VU, IL, IU, ABSTOL,
            M, W_, Z_, LDZ, ISUPPZ_, DWORK_, nDWORK_, IWORK_, nIWORK_, INFO);
        if (INFO) {
            fprintf(stderr, 
                    "Error: eigenvalue computation failed (INFO = %d)\n", 
                    INFO);
            exit(1);
        }
    }
    // Compute ff = 0.5*||X_+||^2 = 0.5*||W||^2
    int INCW = 1;
    double normW = dnrm2_(M, W_, INCW);
    f_ = 0.5*(normW*normW);
    // Compute Z = Z*Diag(W)^{1/2}
    int INCZ = 1;
    double temp;
    for (int j = 0; j < M; j++) {
        // Scale jth column of Z by sqrt(W[j])
        temp = sqrt(W_[j]);
        dscal_(N_, temp, Z_ + j*N_, INCZ);
    }

    char TRANS = 'N';
    double ALPHA = 1.0;
    double BETA = 0.0;

    /* X = ALPHA*Z*Z^T + BETA*X = Z*Z^T
     * Only writes to lower-triangular part of X.
     * When M = 0, we will obtain X = 0.
     */
    dsyrk_(UPLO, TRANS, N_, M, ALPHA, Z_, LDZ, BETA, X_, LDX);

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
    int ineqCon;

    if(mode == TRANSP)
    {
        // models inequalities 
        for(ineqCon = 0; ineqCon < mA_; ++ineqCon) // try a range based loop here too
        {
            dTemp = scaleIneq * y_[mB_ + ineqCon];

            for(auto& it : As_[ineqCon])
            {
                if(it.i_ == it.j_)
                {
                    X_[it.i_ + it.j_ * N_] -= dTemp * it.dVal_;
                }
                else
                {
                    X_[it.i_ + it.j_ * N_] -= dTemp * it.dVal_;
                    X_[it.j_ + it.i_ * N_] -= dTemp * it.dVal_;
                } 
            }
        }
        // cut inequalities
        ineqCon = 0;
        for(auto& it : Cuts_)
        {
            dTemp = 0.5 * scaleIneq * y_[mB_ + mA_ + ineqCon++];
            
            switch (it.type_)
            {

            case ONE:
            {
                X_[it.i_ + it.j_ * N_] += dTemp;
                X_[it.i_ + it.k_ * N_] += dTemp;
                X_[it.j_ + it.k_ * N_] += dTemp;
            } break;
            
            case TWO:
            {
                X_[it.i_ + it.j_ * N_] += dTemp;
                X_[it.i_ + it.k_ * N_] -= dTemp;
                X_[it.j_ + it.k_ * N_] -= dTemp;
            } break;

            case THREE:
            {
                X_[it.i_ + it.j_ * N_] -= dTemp;
                X_[it.i_ + it.k_ * N_] += dTemp;
                X_[it.j_ + it.k_ * N_] -= dTemp;
            } break;

            default: // FOUR
            {
                X_[it.i_ + it.j_ * N_] -= dTemp;
                X_[it.i_ + it.k_ * N_] -= dTemp;
                X_[it.j_ + it.k_ * N_] += dTemp;
            } break;
            }
        }
    }
    else
    {
        for(ineqCon = 0; ineqCon < mA_; ++ineqCon)
        {

            f_ += scaleIneq * a_[ineqCon] * y_[mB_ + ineqCon];

            g_[mB_ + ineqCon] = 0.0;

            for(auto& it : As_[ineqCon])
            {
                if(it.i_ = it.j_)
                {
                    g_[mB_ + ineqCon] -= it.dVal_ * X_[it.i_ + it.j_ * N_];
                }
                else
                {
                    g_[mB_ + ineqCon] -= 2.0 * it.dVal_ * X_[it.i_ + it.j_ * N_];
                }
            }
            g_[mB_ + ineqCon] *= dAlphaInv;
            g_[mB_ + ineqCon] += a_[ineqCon];
            g_[mB_ + ineqCon] *= scaleIneq;
        }

        // cut inequalities
        ineqCon = 0;
        for(auto& it : Cuts_)
        {
            f_ += scaleIneq * y_[mB_ + mA_ + ineqCon];
            
            switch (it.type_)
            {

            case ONE:
            {
                X_[it.i_ + it.j_ * N_] += dTemp;
                X_[it.i_ + it.k_ * N_] += dTemp;
                X_[it.j_ + it.k_ * N_] += dTemp;
            } break;
            
            case TWO:
            {
                X_[it.i_ + it.j_ * N_] += dTemp;
                X_[it.i_ + it.k_ * N_] -= dTemp;
                X_[it.j_ + it.k_ * N_] -= dTemp;
            } break;

            case THREE:
            {
                X_[it.i_ + it.j_ * N_] -= dTemp;
                X_[it.i_ + it.k_ * N_] += dTemp;
                X_[it.j_ + it.k_ * N_] -= dTemp;
            } break;

            default: // FOUR
            {
                X_[it.i_ + it.j_ * N_] -= dTemp;
                X_[it.i_ + it.k_ * N_] -= dTemp;
                X_[it.j_ + it.k_ * N_] += dTemp;
            } break;
            }

            g_[mB_ + mA_ + ineqCon] = (dTemp * dAlphaInv + 1.0) * scaleIneq;

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

    //std::vector<BiqTriInequality>::iterator it = Bs_.begin();

    if(mode == TRANSP)
    {
        for(int eqCon = 0; eqCon < mB_; ++eqCon)
        {
            dTemp = scaleEq * y_[eqCon];
            for(auto& it : Bs_[eqCon])
            {
                if(it.i_ == it.j_)
                {
                    X_[it.i_ + it.j_ * N_] -= dTemp * it.dVal_;
                }
                else
                {
                    X_[it.i_ + it.j_ * N_] -= dTemp * it.dVal_; // why not * 2.0 and skip the next line 
                    X_[it.j_ + it.i_ * N_] -= dTemp * it.dVal_;
                }
            }
        }
    }
    else
    {
        for(int eqCon = 0; eqCon < mB_; ++eqCon)   
        {
            f_ +=  scaleEq * b_[eqCon] * y_[eqCon];

            g_[eqCon] = 0.0;
            for(auto& it: Bs_[eqCon])
            {
                if(it.i_ == it.j_)
                {
                    g_[eqCon] -= it.dVal_ * X_[it.i_ + it.j_ * N_];
                }
                else
                {
                    g_[eqCon] -= 2.0 * it.dVal_ * X_[it.i_ + it.j_ * N_];
                }
            }
            g_[eqCon] *= dAlphaInv;
            g_[eqCon] += b_[eqCon];
            g_[eqCon] *= scaleEq;
        }
    }

}

/// @brief 
/// @return 
bool BiqModel::Prune()
{
    // TODO add the Prune method

    return false;
}

/// @brief 
void BiqModel::CreateSubProblem()
{
    // temp data for building
    BiqVarStatus* pbiqVarStatus = new BiqVarStatus[nVar_];
    int *piSol = new int[nVar_];

    std::fill_n(pbiqVarStatus, nVar_, BiqVarFree);
    std::fill_n(piSol, nVar_, 0);

/*
    pbiqVarStatus[0] = BiqVarFixedToZero;
    piSol[0] = 0;
    pbiqVarStatus[1] = BiqVarFixedToOne;
    piSol[1] = 1;
*/
    int nFixed = 0;

    BuildConstraints(mB_, b_, Bs_, b_sub_, Bs_sub_,
                        piSol, pbiqVarStatus, nFixed); 

    buildObjective(piSol, pbiqVarStatus, nFixed);
}


void BiqModel::buildObjective(int *piSol, BiqVarStatus *pbiqVarStatus, int nFixed)
{
    int i;
    int nSubVar = nVar_ - nFixed;
    double dConst;
    //
    GetSubMatrix(pbiqVarStatus);
    
    //
    GetLinear(Qs_, piSol, pbiqVarStatus);
    

    // add the constant
    print_symmetric_matrix(pdTmp_sub_, nSubVar);
    print_vector(pdTmpLinear, nSubVar);
}

void BiqModel::GetSubMatrix(BiqVarStatus *pbiqVarStatus)
{

    int nFixed;
    int nFree;
    int i, ii, jj;
    int piOffset[nVar_];
    // figure out off set
    nFixed = GetOffset(piOffset, pbiqVarStatus);

    nFree = nVar_ - nFixed;

    // init the submatrix zero out Q_sub_
    for(i = 0; i < nFree*nFree; ++i)
    {
        pdTmp_sub_[i] = 0.0;
    }

    for(auto & it: Qs_)
    {
        if(pbiqVarStatus[it.i_] == BiqVarFree && pbiqVarStatus[it.j_] == BiqVarFree)
        {
            ii = it.i_ - piOffset[it.i_];
            jj = it.j_ - piOffset[it.j_];
            pdTmp_sub_[ii + jj * nFree] = it.dVal_;
            pdTmp_sub_[jj + ii * nFree] = it.dVal_;
        }
    }
}


/// @brief 
/// @param RHSsource 
/// @param sMatSource 
/// @param RHSdest 
/// @param sMatdest 
void BiqModel::BuildConstraints(int nRows,
                                double *RHSsource, std::vector<Sparse> sMatArraySource,
                                double *RHSdest,   std::vector<Sparse> sMatArraydest,
                                int *piSol, BiqVarStatus *pbiqVarStatus, int nFixed)
{
    int nnzAdded;
    int i;
    double dCoefMatNorm;
    double dScaleFactor;
    double dConstant;
    int nVarSub = nVar_ - nFixed;
    Sparse::iterator itSource;
    
    for(int k = 0; k < nRows; ++k)
    {
        dCoefMatNorm = GetSubMatrixSparse(sMatArraySource.at(k), pbiqVarStatus, nnzAdded);
        dConstant = GetConstant(sMatArraySource.at(k), piSol, pbiqVarStatus);
        GetLinear(sMatArraySource.at(k), piSol, pbiqVarStatus);

        //
        dScaleFactor = 1.0;
        std::printf("About to print the size\n");
        std::printf("Size of dest matrix : %d\n", sMatArraydest.at(k).size());
        // set itterator for fill
        sMatArraydest.at(k).resize(100);
        itSource = (sMatArraydest.at(k)).begin();
        // Quad Part
        for(i = 0; i < nnzAdded; ++i)
        {
            if(fabs(sTmp_sub_.at(i).dVal_) > 0.0)
            {
                itSource->i_    = sTmp_sub_.at(i).i_;
                itSource->j_    = sTmp_sub_.at(i).j_;
                itSource->dVal_ = sTmp_sub_.at(i).dVal_ * dScaleFactor;
                ++itSource;
                std::printf("adding quad part %d\n",i);
            }
        } 

        // Linear Part
        for(int j = 0; j < nVarSub; ++j)
        {
            if(fabs(pdTmpLinear[j]) > 0.0)
            {
                itSource->i_    = nVarSub;
                itSource->j_    = j;
                itSource->dVal_ = pdTmpLinear[j] * dScaleFactor;
                ++itSource;
                std::printf("adding quad part %d\n",j);
            }
        }

        // new constant value on left-side
        if(fabs(dConstant) > 0.0)
        {
            itSource->i_    = nVarSub;
            itSource->j_    = nVarSub;
            itSource->dVal_ = dConstant * dScaleFactor;
            ++itSource;
            std::printf("adding quad part %d\n",i);
        }

        RHSdest[k] *= dScaleFactor;
        PrintSparseMatrix(sMatArraySource.at(k));
        std::printf("BiqModel::BuildConstraints  New Matrix %d\n", k);
        PrintSparseMatrix(sMatArraydest.at(k));
    }
}


/// @brief 
/// @param sSourceMat As_ or Bs_
/// @param pbiqVarStatus this is until the BiqNodeDesc
///         being able to pass the real one
/// @return 
double BiqModel::GetSubMatrixSparse(Sparse sSourceMat, 
                                    BiqVarStatus *pbiqVarStatus, 
                                    int &nnzAdded)
{
    double dRetNorm = 0.0;
    double dVal;
    int nFixed = 0;
    int iArrOffset[nVar_];
    int ii, jj, pos = 0;
    bool bBothFree;
    // figure out the offset
    nFixed = GetOffset(iArrOffset, pbiqVarStatus);
    // zero out nnzAdded to be safe
    nnzAdded = 0;
    // fill in sTmp_sub_
    Sparse::iterator itDest =  sTmp_sub_.begin(); // TODO how to not be a pointer 
    for(auto& itSource : sSourceMat)
    {
        bBothFree =  pbiqVarStatus[itSource.i_] == BiqVarFree && pbiqVarStatus[itSource.j_] == BiqVarFree;
        if(bBothFree && itSource.i_ < nVar_ && itSource.j_ < nVar_)
        {
            ii = itSource.i_ - iArrOffset[itSource.i_];
            jj = itSource.j_ - iArrOffset[itSource.j_];
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
            std::printf(" %f * %f + ", dVal, dVal);
            dRetNorm += dVal*dVal;
            if(ii != jj)
            {
                std::printf(" %f * %f + ", dVal, dVal);
                dRetNorm += dVal*dVal;
            }
            // incress the itDest
            ++itDest;
            ++nnzAdded;
        }
    }
    std::printf("\n");    
    dRetNorm = sqrt(dRetNorm);

    return dRetNorm;
}

/// @brief loop over sMat and add up {-1, 1} const value from {0,1} sol
/// @param sMat 
/// @param piSol 
/// @param piFixed 
/// @return 
double BiqModel::GetConstant(Sparse &sMat, int *piSol, BiqVarStatus *pbiqVarStatus)
{
    double dRetConst = 0.0;
    double dTmp;

    for(auto &it : sMat)
    {

        // if in the quad part
        if(it.i_ < nVar_ && it.j_ < nVar_)
        {
            if(pbiqVarStatus[it.i_] != BiqVarFree && pbiqVarStatus[it.j_] != BiqVarFree)
            {
                std::printf("BiqModel::GetConstant fixed in Quad %d, %d\n",it.i_,  it.j_);
                dTmp = (2*piSol[it.i_] - 1.0)*(2*piSol[it.j_] - 1.0)*it.dVal_;
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
            std::printf("BiqModel::GetConstant fixed in Const %d, %d\n",it.i_,  it.j_);
            dRetConst += it.dVal_;
        }
        // else if in the linear column
        else if(it.j_ == nVar_)
        {
            if(pbiqVarStatus[it.i_] != BiqVarFree)
            {
                std::printf("BiqModel::GetConstant fixed in linear column %d, %d\n",it.i_,  it.j_);
                dRetConst += (2*piSol[it.i_] - 1.0)*2*it.dVal_;
            }
        }
        // else if in the linear row
        else if(it.i_ == nVar_)
        {
            if(pbiqVarStatus[it.j_] != BiqVarFree)
            {
                std::printf("BiqModel::GetConstant fixed in linear row %d, %d\n",it.i_,  it.j_);
                dRetConst += (2*piSol[it.j_] - 1.0)*2*it.dVal_;
            }
        }
    }

    return dRetConst;
}

void BiqModel::GetLinear(Sparse &sSource, int *piSol, BiqVarStatus *pbiqVarStatus)
{
    double dSum = 0.0;
    double dVal = 0.0;
    int nFixed = 0;
    int nFree = 0;
    int iArrOffset[nVar_];
    int ii, jj, pos = 0;

    bool bIfree, bJfree;
    // figure out the offset
    nFixed = GetOffset(iArrOffset, pbiqVarStatus);

    nFree = nVar_ - nFixed;

    for(int i = 0; i < nFree; ++i)
    {
        pdTmpLinear[i] = 0.0;
    }

    for(auto &itSource : sSource)
    {
        ii = itSource.i_;
        jj = itSource.j_;
        dVal = itSource.dVal_;

        bIfree = pbiqVarStatus[ii] == BiqVarFree;
        bJfree = pbiqVarStatus[jj] == BiqVarFree;  
        // if we are already linear 
        if(ii < nVar_ && bIfree && jj == nVar_)
        {
            pdTmpLinear[ii - iArrOffset[ii]] += dVal;
        }
        else if(jj < nVar_ && bJfree && ii == nVar_)
        {
            pdTmpLinear[jj - iArrOffset[jj]] += dVal;
        }
        // if in quad part
        else if(ii < nVar_ && jj < nVar_)
        {
            if(bIfree && !bJfree)
            {
                pdTmpLinear[ii - iArrOffset[ii]] += (2 * piSol[jj] - 1.0)*dVal;
            }
            else if(!bIfree && bJfree)
            {
                pdTmpLinear[jj - iArrOffset[jj]] += (2 * piSol[ii] - 1.0)*dVal;
            }
        }

    }

}

int BiqModel::GetOffset(int *piOffset, BiqVarStatus *pbiqVarStatus)
{
    int nFixedRet = 0;
    
    for(int i = 0; i < nVar_; ++i)
    {
        if(pbiqVarStatus[i] != BiqVarFree)
        {
            piOffset[i] = 0;
            ++nFixedRet;  
        }
        else
        {
            piOffset[i] = nFixedRet;
        }
    }

    return nFixedRet;
}