#include "headers/BiqModel.h"
#include "AlpsKnowledgeBroker.h"
#include "AlpsKnowledge.h"
#include "headers/BiqSolution.h"
BiqModel::BiqModel(
                int nVar, double *Q, bool max_problem,
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
    vOffset_.resize(nVar_, 0);
    FillSparseMatrix(Qs_,Q_,N_);
    // Copy the b vector in to the first Bs_.size() entries
    std::copy(b, b + Bs_.size(), b_);
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
    
    AddDiagCons();
    AllocSubCons();

    // Set cut data 
    Cuts_.resize(MaxNineqAdded); // or start small and add space when needed
    nIneq_ = 0;

    // set some memory for the sTemp_sub_ matrix
    sTmp_sub_.resize(N_*N_);
    pdTmp_sub_ = new double[N_*N_];
    pdTmpLinear_ = new double[nVar_];
    delete[] IWORK_SIZE;
    delete[] DWORK_SIZE;

    

    
}

BiqModel::~BiqModel()
{
    std::printf("BiqModel::~BiqModel()\n");
    // TODO clean up and make sure this is safe
    if(ISUPPZ_)
    {
        delete[] ISUPPZ_;
        ISUPPZ_ = NULL;
    }
    if(X_ )
    {
        delete[] X_ ;
        X_  = NULL;
    }
    if(W_)
    {
        delete[] W_;
        W_ = NULL;
    }
    if(Z_)
    {
        delete[] Z_;
        Z_ = NULL;
    }
    if(DWORK_)
    {
        delete[] DWORK_;
        DWORK_ = NULL;
    }
    if(IWORK_)
    {
        delete[] IWORK_;
        IWORK_ = NULL;
    }
    if(a_)
    {
        delete[] a_;
        a_ = NULL;
    }
    if(b_)
    {
        delete[] b_;
        b_ = NULL;
    }
    if(a_sub_)
    {
        delete[] a_sub_;
        a_sub_ = NULL;
    }
    if(b_sub_)
    {
        delete[] b_sub_;
        b_sub_ = NULL;
    }
    if(Q_sub_)
    {
        delete[] Q_sub_;
        Q_sub_ = NULL;
    }
    if(g_)
    {
        delete[] g_;
        g_ = NULL;
    }
    if(y_)
    {
        delete[] y_;
        y_ = NULL;
    }
    if(binf_)
    {
        delete[] binf_;
        binf_ = NULL;
    }
    if(bsup_)
    {
        delete[] bsup_;
        bsup_ = NULL;
    }
    if(nbd_)
    {
        delete[] nbd_;
        nbd_ = NULL;
    }
    if(wa_)
    {
        delete[] wa_;
        wa_ = NULL;
    }
    if(iwa_)
    {
        delete[] iwa_;
        iwa_ = NULL;
    }
    if(pdTmp_sub_)
    {
        delete[] pdTmp_sub_;
        pdTmp_sub_ = NULL;
    }
    if(pdTmpLinear_)
    {
        delete[] pdTmpLinear_;
        pdTmpLinear_ = NULL;
    }
    std::printf("BiqModel::~BiqModel() done\n");
}

void BiqModel::freeData(int *& data)
{
    std::printf("BiqModel::freeData int\n");
    if(data)
    {
        delete[] data;
        data = NULL;
    }
  
}

void BiqModel::freeData(double *& data)
{
    std::printf("BiqModel::freeData double\n");
    if(data)
    {
        delete[] data;
        data = NULL;
    }
  
}

/// @brief 
/// @param encode 
/// @return 
AlpsKnowledge * BiqModel::decode(AlpsEncoded & encode) const
{
    std::printf("AlpsKnowledge * BiqModel::decode\n");
    std::cerr << "Not implemented!" << std::endl;
    throw std::exception();
}

/// @brief 
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
        //PrintSparseMatrix(bstTemp);
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
    int nAdded;
    int nSubtracted;
    int nHeurRuns = 5;
    double dMinAllIneq;
    double dAlpha = 10;
    double dTol = 0.08;
    double dRetBound = 0.0;
    int len_y = mA_ + mB_ + nIneq_;
    // param
    int nMinAdded = 50;
int nMaxIter = 30;
    double dScaleAlpha = 0.5;
    double dScaleTol = 0.93;

    // 
    if(nVar_sub_ == 0)
    {
        dRetBound = Q_sub_[0];
        std::printf("N = 1: returning bound = %.1f\n", dRetBound); // TODO fprintf output
        return dRetBound;
    }

    for(i = 0; i < len_y; ++i)
    {
        y_[i] = 0.0;
    }
    
    sim(dAlpha);
    for(int i = 0; i <= nMaxIter; ++i)
    {

        // Call BFGS solver
        iStatus = CallLBFGSB(dAlpha, dTol, nbit);

        // run heur 
        GWheuristic(nHeurRuns);
        dMinAllIneq = UpdateInequalities(nAdded, nSubtracted);
        PrintBoundingTable(i+1, nbit, nAdded, nSubtracted, dAlpha, dTol, dMinAllIneq);

        if(nAdded < nMinAdded)
        {
                dAlpha *= dScaleAlpha;
                dTol *= dScaleTol;
        }
    }
    
    return dRetBound;
}

void BiqModel::PrintBoundingTable(int iter, int nBit, int nAdded, int nSubtracted, double dAlpha, double dTol, double dMinAllIneq /*double dTime*/)
{
    if(iter == 1)
    {
        std::printf("======================================================================================\n");
        std::printf("%4s  %6s  %5s  %5s  %5s  %4s  %5s  %5s  %5s  %6s  %4s  %4s  %4s\n", 
                    "Iter", 
                    "Time", 
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
        std::printf("======================================================================================\n");
    }

    std::printf("%4d  %6.1f  %5.1f  %5.0e  %5.0e  %4d  %5.0e  %5.0e  %5.0e  %6.0e  %4d  -%-3d  +%-3d\n", 
                iter, 
                0.0, /* TODO time*/
                f_, 
                dAlpha, 
                dTol, 
                nBit, 
                gradEnorm_,
                gradInorm_,
                0.0, /* TODO ynorm*/
                dMinAllIneq, 
                nIneq_, 
                nSubtracted,
                nAdded
    );

}
int BiqModel::CallLBFGSB(double dAlpha, double dTol, int &nbit)
{
    int retStatus = 1;
    int mem = mmax;
    int i;
    int N = nVar_sub_ + 1;
    bool bStopBFGS = false;
    bool bPrune;
    double factr = 0.000005;
    double pgtol = 0.0;
    double dMinTemp;
    double dBound;
    char task[60];
    int iprint = -1;
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

           /*
           std::printf("f = %f\n",f_);
           std::printf("=====\n");
           print_vector(g_, len_y);
           exit(0);
           */
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

            
            //std::printf("GradENorm = %f\n", gradEnorm_);
            //exit(0);
            

            dBound = (max_problem_) ? f_  : - f_;
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
    int N2 = N*N;
    int incx = 1; 
    dscal_(N2, dAlphaInv, X_, incx);
    // Update Z so that X = Z*Z'
    double scb = 10000.0 * dAlpha;
    double sca = 0.01/ sqrt(scb);
    int NM = N * M_;
    int incz = 1;
    dscal_(NM, sca, Z_, incz); 
    //std::printf("LBFGSB bound = %f\n", dBound);
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

    //std::printf("After ProjSdp f_ = %f\n", f_);

    f_ = f_ * 1000.0 * (dAlphaInv / 1000.0); // TODO see what happens 

    /*
    print_symmetric_matrix(X_, N);
    exit(0);
    */

    // call the A and B methods not transposed
    B(NOTRANSP, alpha);
    A(NOTRANSP, alpha);

    f_ += alpha * (0.5 * static_cast<double>(N2));

    //print_vector(g_,mA_+mB_+nIneq_);
    //exit(0);
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
    double ABSTOL = 0; // TODO could be a param 
    int LDZ = N;
    int INFO;

    /*
    std::printf("Printing X_..\n");
    print_symmetric_matrix(X_,N);
    exit(0);
    */
    //std::printf("") /// print out each input // which dsyevr_ // is there any rng used??

    /*
    std::printf("N = %d\t VL = %f\t VU = %f\t nDWORK_ = %d\t nIWORK_ = %d\n",
                 N, VL, VU, nDWORK_, nIWORK_);
    exit(0);
    */

    dsyevr_(JOBZ, RANGE, UPLO, N, X_, LDX, VL, VU, IL, IU, ABSTOL,
            M_, W_, Z_, LDZ, ISUPPZ_, DWORK_, nDWORK_, IWORK_, nIWORK_, INFO);

   
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

    /*
    std::printf("N = %d \t M = %d\n", N, M_);
    print_vector(Z_, N*M_);
    exit(0);
    */

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

    /*
    std::printf("Eigen Values\n");
    print_vector(W_, M_);
    std::printf("=====\n");
    exit(0);
    */

    char TRANS = 'N';
    double ALPHA = 1.0;
    double BETA = 0.0;

    /* X = ALPHA*Z*Z^T + BETA*X = Z*Z^T
     * Only writes to lower-triangular part of X.
     * When M = 0, we will obtain X = 0.
     */
    dsyrk_(UPLO, TRANS, N, M_, ALPHA, Z_, LDZ, BETA, X_, LDX);

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
    int N = nVar_sub_ + 1;

    if(mode == TRANSP)
    {
        // models inequalities 
        for(ineqCon = 0; ineqCon < mA_; ++ineqCon) // try a range based loop here too
        {
            dTemp = scaleIneq * y_[mB_ + ineqCon];

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
            dTemp = 0.5 * scaleIneq * y_[mB_ + mA_ + ineqCon++];
            
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
        }
    }
    else
    {
        for(ineqCon = 0; ineqCon < mA_; ++ineqCon)
        {

            f_ += scaleIneq * a_sub_[ineqCon] * y_[mB_ + ineqCon];

            g_[mB_ + ineqCon] = 0.0;

            for(auto& it : As_sub_[ineqCon])
            {
                if(it.i_ = it.j_)
                {
                    g_[mB_ + ineqCon] -= it.dVal_ * X_[it.i_ + it.j_ * N];
                }
                else
                {
                    g_[mB_ + ineqCon] -= 2.0 * it.dVal_ * X_[it.i_ + it.j_ * N];
                }
            }
            g_[mB_ + ineqCon] *= dAlphaInv;
            g_[mB_ + ineqCon] += a_sub_[ineqCon];
            g_[mB_ + ineqCon] *= scaleIneq;
        }

        // cut inequalities
        ineqCon = 0;
        
        for(auto it = Cuts_.begin(); it < Cuts_.begin() + nIneq_; ++it)
        {
            
            f_ += scaleIneq * y_[mB_ + mA_ + ineqCon];

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
            /*
            std::printf("i = %d  j = %d  k = %d type = %d\t  evaluated = %f\n", 
                        it->i_, it->j_, it->k_, static_cast<int>(it->type_), dTemp);
            */
            g_[mB_ + mA_ + ineqCon] = (dTemp * dAlphaInv + 1.0) * scaleIneq;
            ineqCon++;
        }
        //if(nIneq_>0) exit(0);
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
    //std::vector<BiqTriInequality>::iterator it = Bs_.begin();

    if(mode == TRANSP)
    {
        for(int eqCon = 0; eqCon < mB_; ++eqCon)
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
        
        for(int eqCon = 0; eqCon < mB_; ++eqCon)   
        {
            f_ +=  scaleEq * b_sub_[eqCon] * y_[eqCon];

            g_[eqCon] = 0.0;
            //PrintSparseMatrix(*itBs);
            //std::printf("-----------\n");
            for(auto &it : *itBs)
            {
                if(it.dVal_ == 0) continue;

                /*
                std::printf("X = %f \t dval = %f \n",
                            X_[it.i_ + it.j_ * N], it.dVal_);
                */
                
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
/// @return 
bool BiqModel::Prune()
{
    // TODO add the Prune method

    return false;
}


/// @brief 
void BiqModel::CreateSubProblem(std::vector<BiqVarStatus> vbiqVarStatus)
{
    int nFixed = GetOffset(vbiqVarStatus);

    nVar_sub_ = nVar_ - nFixed;

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

            }
        } 

        // Linear Part
        for(int j = 0; j < nVar_sub_; ++j)
        {
            if(fabs(pdTmpLinear_[j]) > 0.0)
            {
                itSource->i_    = nVar_sub_;
                itSource->j_    = j;
                itSource->dVal_ = pdTmpLinear_[j] * dScaleFactor;
                ++itSource;
            }
        }

        // new constant value on left-side
        if(fabs(dConstant) > 0.0)
        {
            itSource->i_    = nVar_sub_;
            itSource->j_    = nVar_sub_;
            itSource->dVal_ = dConstant * dScaleFactor;
            ++itSource;
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
    int piOffset[nVar_];

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
    int ii, jj, pos = 0;
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

    for(auto &it : sMat)
    {

        // if in the quad part
        if(it.i_ < nVar_ && it.j_ < nVar_)
        {
            if(vbiqVarStatus.at(it.i_) != BiqVarFree && vbiqVarStatus.at(it.j_) != BiqVarFree)
            {
                dTmp = (2*vbiqVarStatus.at(it.i_) - 1.0)*(2*vbiqVarStatus.at(it.j_) - 1.0)*it.dVal_;
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
                dRetConst += (2*vbiqVarStatus.at(it.i_)- 1.0)*2*it.dVal_;
            }
        }
        // else if in the linear row
        else if(it.i_ == nVar_)
        {
            if(vbiqVarStatus[it.j_] != BiqVarFree)
            {
                dRetConst += (2*vbiqVarStatus.at(it.j_)- 1.0)*2*it.dVal_;
            }
        }
    }

    return dRetConst;
}

void BiqModel::GetLinear(Sparse &sSource, std::vector<BiqVarStatus> vbiqVarStatus, int nFixed)
{
    double dSum = 0.0;
    double dVal = 0.0;
    int ii, jj, pos = 0;
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
                pdTmpLinear_[ii - vOffset_.at(ii)] += (2 * vbiqVarStatus.at(jj) - 1.0)*dVal;
            }
            else if(!bIfree && bJfree)
            {
                pdTmpLinear_[jj - vOffset_.at(jj)] += (2 * vbiqVarStatus.at(ii) - 1.0)*dVal;
            }
        }

    }

}

int BiqModel::GetOffset(std::vector<BiqVarStatus> vbiqVarStatus)
{
    int nFixedRet = 0;
    
    for(int i = 0; i < vOffset_.size(); ++i)
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
    double dIneqVal;
    int yIndex;
    nSubtracted = 0;
    nAdded = 0;
    BiqTriInequality btiTemp;
    BiqTriTuple bttTemp;
    TriMap::iterator itMap;
    TriCuts::iterator itIneq;
    TriCuts::iterator itNextIneq;
    

    // Fill the Heap_ with most violated cuts
    dRetVal = GetViolatedCuts();

    // remove cuts from Cuts_ array
    itNextIneq = Cuts_.begin();
    yIndex = mB_ + mA_;
    for(itIneq = Cuts_.begin(); itIneq < Cuts_.begin()+nIneq_; ++itIneq) 
    {
        dIneqVal = EvalInequalities(itIneq->type_, itIneq->i_, itIneq->j_, itIneq->k_);

        itIneq->y_ = y_[yIndex]; yIndex++;
        // if multiplier is small and inequality is large
        if(itIneq->y_ < 1e-8 && dIneqVal > gradInorm_/10.0)
        {
            ++nSubtracted;
        }
        else
        {
            itNextIneq = itIneq;
            ++itNextIneq;
        }
    }

    while (!Heap_.empty())
    {
        // get least violated cut
        btiTemp = Heap_.top();
        //std::printf("Cut[%d, %d, %d] \t value = %f \t type = %d\n", btiTemp.i_, btiTemp.j_, btiTemp.k_, btiTemp.value_, btiTemp.type_);
        Heap_.pop();
        if(itNextIneq != Cuts_.end())
        {
            // update the map so we know this cut has been included
            bttTemp =  BiqTriTuple(btiTemp.type_, 
                                btiTemp.i_, 
                                btiTemp.j_, 
                                btiTemp.k_
            );
            
            itMap = Map_.find(bttTemp);
            if(itMap == Map_.end())
            {
                Map_.insert({bttTemp, true});
                // add it to the cut array
                itNextIneq->i_ = btiTemp.i_;
                itNextIneq->j_ = btiTemp.j_;
                itNextIneq->k_ = btiTemp.k_;
                itNextIneq->type_ = btiTemp.type_;
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
    //std::printf("Added: %d \t Subtracted: %d \t nIneq = %d\n",nAdded, nSubtracted, nIneq_);
    return dRetVal;
}

/// @brief 
/// @return 
double BiqModel::GetViolatedCuts()
{ 
    double dRetMinIneq = INFINITY;
    double dTestIneq;
    int i, j, k, type, count = 0;
    // TODO make following data a param
    double dGapCuts = -5e-2;
    int nCuts = 500;
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
                

                    btiTemp = BiqTriInequality(
                                    static_cast<TriType>(type),
                                    i,j,k,
                                    dTestIneq,
                                    0.0);

                    if(Heap_.size() < nCuts)
                    {
                        // (1) add cuts until full
                        // check if the cut has been added .. to the map
                        if(dTestIneq < dGapCuts)
                            Heap_.push(btiTemp);

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

double BiqModel::EvalInequalities(TriType triType, int ii, int jj, int kk)
{
    double dRetIneqVal = 0.0;

    // compute the inequality
    switch (triType)
    {
        case ONE:
        {
            dRetIneqVal = X_[ii+jj*(nVar_sub_ + 1)] + X_[ii+kk*(nVar_sub_ + 1)] + X_[jj+kk*(nVar_sub_ + 1)] + 1.0;
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

double BiqModel::EvalSolution(std::vector<int> solution)
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

    return dRetSol;
}


bool BiqModel::isFeasibleSolution(std::vector<int> solution)
{

    int m, i, j, k;
    double dSum, dTmp;
    bool bRet = true;
    // convert to {1, -1}
    for(auto &it: solution)
    {
        it = 2*it - 1;
    } 

    solution.at(nVar_) = 1; 

    // loop over equality cons of master problem;
    i = 0;
    for(auto &itCons : Bs_)
    {
        dSum  = 0;
        for(auto &it : itCons)
        {
            dTmp = it.dVal_ * solution.at(it.i_) * solution.at(it.j_);
            dSum += dTmp;
            if(it.i_ != it.j_)
            {
                dSum += dTmp;
            }
        }

        // check that the sum matches the rhs
        if(dSum != b_[i])
        {
            bRet = false;
            break;
        }
        ++i; 
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
            bRet = false;
            break;
        }
        ++i; 
    }




    return bRet;
}

double BiqModel::primalHeuristic(std::vector<int> solution)
{
    double dRet = 0.0;
    double gamma;
    double dRand;
    bool bTest;
    for(gamma = 0.0; gamma < 1.0; gamma += 0.01)
    {
        for(auto &it: solution)
        {
            dRand = ((double) rand() / (double) RAND_MAX);
            if(dRand > gamma)
            {
                it = 1;
            }
            else
            {
                it = 0;
            }
        }
        bTest = isFeasibleSolution(solution);
        std::printf("Do we have a solution?? %d\n", bTest);
    }

    return dRet;
}

double BiqModel::GWheuristic(int nPlanes)
{
    double dRetBest = -1e+9;
    double dPlaneNorm;
    double dPlaneInvNorm;
    double sca;
    double dTmpVal;

    int index;
    int bestVal = 0;

    std::vector<double> hyperPlane;
    std::vector<int> solution_1;
    std::vector<int> solution_2;

    


    hyperPlane.resize(M_);
    solution_1.resize(nVar_+1);
    solution_2.resize(nVar_+1);

    solution_1.at(nVar_) = 1;
    solution_2.at(nVar_) = 1;

    // temp have a solutin vector
    std::vector<BiqVarStatus> vbiqVarStatus;
    vbiqVarStatus.resize(nVar_, BiqVarFree);
    bool bSol1Feasible;
    bool bSol2Feasible;
    

    for(int k = 0; k < nPlanes; ++k)
    {
        // create a random hyperplane
        dPlaneNorm = 0.0;
        for(auto &it: hyperPlane)
        {
            it = 1.0 + 100.0 * static_cast<double>(rand())/(static_cast<double>(RAND_MAX) + 1.0)
;            dPlaneNorm += it*it;
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
                solution_1.at(i) = 1;
                solution_2.at(i) = 1;
            }break;

            case BiqVarFixedToZero:
            {
                solution_1.at(i) = 0;
                solution_2.at(i) = 0;
            }break;

            case BiqVarFree:
            {
                sca = 0.0;
                for(int j = 0; j < M_; ++j)
                {
                    sca += hyperPlane.at(j)*Z_[j * nVar_sub_ + index];
                }
                if(sca < 0)
                {
                    solution_1.at(i) = 0;
                    solution_2.at(i) = 1;
                }
                else
                {
                    solution_1.at(i) = 1;
                    solution_2.at(i) = 0;  
                }

                ++index;
            }break;


            default:
                break;
            }

        }
        UpdateSol(solution_1);
        UpdateSol(solution_2);
    }


    return dRetBest;
}

void BiqModel::UpdateSol(std::vector<int> solution)
{
    bool bIsfeasible;
    // check if feasible
    bIsfeasible = isFeasibleSolution(solution);
    if(!bIsfeasible) return;
    

    int bestVal = 0;
    double dTmpVal;

    // The quality of a solution is the negative of the objective value
    //  since Alps consideres sols with lower quality values better.
    //bestVal = static_cast<int>(broker()->getIncumbentValue());
    if(max_problem_)
    {
        bestVal = -bestVal;
    }

    dTmpVal = EvalSolution(solution);
    //BiqSolution* biqSol = new BiqSolution( this, solution, dTmpVal);
    if(max_problem_ && dTmpVal > bestVal)
    {
        //broker()->addKnowledge(AlpsKnowledgeTypeSolution, biqSol, -dTmpVal);
    }
    else if(!max_problem_ && dTmpVal < bestVal)
    {     
        //broker()->addKnowledge(AlpsKnowledgeTypeSolution, biqSol, dTmpVal);
    }   
    else
    {
        /* do nothing*/
    }
}

void BiqModel::InitEmptyModel()
{
    ISUPPZ_ = 0;
    X_ = 0;
    W_ = 0;
    Z_ = 0;
    DWORK_ = 0;
    IWORK_ = 0;
    a_ = 0;
    b_ = 0;
    b_sub_ = 0;
    a_sub_ = 0;
    Q_sub_ = 0;
    g_ = 0;
    y_ = 0; 
    binf_ = 0;
    bsup_ = 0; 
    nbd_ = 0; 
    wa_ = 0; 
    iwa_ = 0;
    pdTmpLinear_ = 0;
    pdTmp_sub_ = 0;
}