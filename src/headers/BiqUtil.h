#ifndef BiqUtil_h_
#define BiqUtil_h_

#include <vector>
#include <list>
#include <map>
#include <tuple>
#include <queue>
#include <unordered_map>
#include <string>
#define TRANSP 1
#define NOTRANSP 0
#define scaleEq 1.0
#define scaleIneq 1.0 / (1.0 + sqrt(3.0)) // scales the inequality constraints (and dual vars)
#define MaxNineqAdded 100000
#define mmax 10
#define MAXITER 10000
#define INFVAL 2147483648

#define FREE_DATA(data)\
    if(data)\
    {\
        delete[] data;\
        data = NULL;\
    }

enum BiqVarStatus {
    BiqVarFree = -1,
    BiqVarFixedToOne = 1,
    BiqVarFixedToZero = 0
};

enum TriType
{
  ONE = 1,
  TWO,
  THREE,
  FOUR
};


class BiqModel;

class BiqSparseTriple
{
public:
    int i_;
    int j_;
    double dVal_;

    BiqSparseTriple(int i, int j, double dVal)
        : i_(i),
          j_(j),
          dVal_(dVal)
    {
    }
    BiqSparseTriple()
        : i_(0),
          j_(0),
          dVal_(0.0)
    {
    }
};

class BiqTriInequality
{
public: 
  TriType type_;
  int i_, j_, k_;
  double value_;
  double y_; // more discriptive name would be good

  BiqTriInequality(TriType type, int i, int j, int k, double value, double y)
    : type_(type),
      i_(i),
      j_(j),
      k_(k),
      value_(value),
      y_(y)
      {
      }
  BiqTriInequality()
     :type_(ONE),
      i_(-1),
      j_(-1),
      k_(-1),
      value_(0),
      y_(0)
      {
      }

      friend bool operator<(BiqTriInequality const &l, BiqTriInequality const &r);

};





using Sparse = std::vector<BiqSparseTriple>;
using BiqTriTuple = std::tuple<int,int,int,int>;
using TriCuts = std::vector<BiqTriInequality>;


/*
// unordered_map test does not seem like it is worth using
struct BiqTriTupleHash {
    std::size_t operator()(const BiqTriTuple& tuple) const {
        int key;
        int spacing = 100;
        key = std::get<0>(tuple) + spacing + std::get<1>(tuple) + (spacing*2) + std::get<2>(tuple) + (spacing*3) + std::get<3>(tuple);
        return std::hash<int>()(key);
    }
};

using TriMap  = std::unordered_map<BiqTriTuple, bool, BiqTriTupleHash>;
*/

using TriMap  = std::map<BiqTriTuple, bool>;
using TriHeap = std::priority_queue<BiqTriInequality, std::vector<BiqTriInequality>, std::less<BiqTriInequality>>;



int I(BiqTriTuple btt);
int J(BiqTriTuple btt);
int K(BiqTriTuple btt);
int TRI_TYPE(BiqTriTuple btt);

void FillSparseMatrix(Sparse& sMat, std::vector<double> data);

void PrintSparseMatrix(std::vector<BiqSparseTriple> bstVec);

void PrintMatrix(double * pdMat, int nRow, int nCol);

void print_vector(double *vec, int N);

void print_symmetric_matrix(double *Mat, int N);

void PrintMatrix(std::vector<double> mat);

void zeroOutMatrix(std::vector<double> &mat);





extern "C" int openblas_get_num_threads();
extern "C" void openblas_set_num_threads(int n);
extern "C" double dnrm2_(
        int &n, double *x, int &incx);

extern "C" void dscal_(
        int &n, double &da, double *x, int &incx);

extern "C" void dcopy_(
        int &n, double *dx, int &incx, double *dy, int &incy);

extern "C" void dsyevr_(
        char &JOBZ, char &RANGE, char &UPLO, int &N, double *A, int &LDA, 
        double &VL, double &VU, int &IL, int &IU, double &ABSTOL, int &M, 
        double *W, double *Z, int &LDZ, int *ISUPPZ, double *WORK, int &LWORK, 
        int *IWORK, int &LIWORK, int &INFO);

extern "C" void dsyrk_(
        char &UPLO, char &TRANS, int &N, int &K, double &ALPHA, double *A, 
        int &LDA, double &BETA, double *C, int &LDC);

extern "C" double dlansy_(
        char &NORM, char &UPLO, int &N, double *A, int &LDA, double *WORK);

extern "C" void setulb_( 
        int &n, // number of variables
        int &m, // size of the limited memory
        double *x, // current iterate (length n)
        double *l, // lower bounds (length n)
        double *u, // upper bounds (length n)
        int *nbd, // indicates which vars are bounded
        // nbd[i] = 0 : x[i] unbounded
        // nbd[i] = 1 : x[i] has only a lower bound
        // nbd[i] = 2 : x[i] has both lower and upper bounds
        // nbd[i] = 3 : x[i] has only an upper bound
        double &f, // value of function at the point x
        double *g, // value of the gradient at the point x (length n)
        double &factr, // termination tolerance
        // factr = 1.d+12 for low accuracy;
        //         1.d+7  for moderate accuracy;
        //         1.d+1  for extremely high accuracy.
        double &pgtol, // projected gradient tolerance
        // (suppress this termination test by pgtol = 0)
        double *wa, // workspace (length (2mmax + 4)nmax + 12mmax^2 + 12mmax)
        int *iwa, // workspace (length 3nmax)
        char *task, // character array (length 60)
        // 'START'    : when starting
        // 'FG'        : user must evaluate function f and gradient g
        // 'NEW_X'    : user can decide whether to continue or stop
        // 'CONV'    : termination test has been satisfied
        // 'ABNO'    : abnormal termination
        // 'ERROR'    : error in input parameters
        // 'STOP'    : set by user to stop L-BFGS-B
        int &iprint, // set level of output
        // iprint<0    no output is generated;
        // iprint=0    print only one line at the last iteration;
        // 0<iprint<99 print also f and |proj g| every iprint iterations;
        // iprint=99   print details of every iteration except n-vectors;
        // iprint=100  print also the changes of active set and final x;
        // iprint>100  print details of every iteration including x and g;
        char *save, // character array (length 60)
        int *lsave, // logicial array (length 4)
        int *isave, // integer array (length 44)
        double *dsave // double array (length 29)
        );

#endif