#ifndef BiqModel_h_
#define BiqModel_h_

#include "AlpsKnowledge.h"
#include "AlpsModel.h"
#include "BiqTreeNode.h"
#include "BiqNodeDesc.h"
#include "BiqSolution.h"
#include "BiqUtil.h"
#include "math.h"
#include<iostream>
#include <utility>
#include <list>
#include <map>
#include <unordered_map>
#include <tuple>
#include <queue>


class BiqModel: public AlpsModel
{
private:
    /* data */
    double *Q_;                           // Dense Objective matrix   Q_[i + j*nVar_] to get the i,j
    Sparse Qs_;
    int N_;
    int nVar_;                            // Number of variables

    // array of array
    std::vector<Sparse> As_;     // Inequality constraint matrix   
    double *a_;                            // Right-hand-side inequality constraints
    int mA_;                               // Number of inequality constraints 

    std::vector<Sparse> Bs_;     // Equality constraint matrix
    double *b_;                            // Right-hand-side equality constraints
    int mB_;                               // Number of equality constraints
                                 
    bool max_problem_;                      // Objective sense 

    /* sub model data*/
    double *X_;
    std::vector<int> vOffset_;

    /* data */
    double *Q_sub_;                           // Dense Objective matrix   Q_[i + j*nVar_] to get the i,j
    int nVar_sub_;                            // Number of variables

    // array of Sparse
    std::vector<Sparse> As_sub_;     // Inequality constraint matrix   
    double *a_sub_;                            // Right-hand-side inequality constraints
    int mA_sub_;                               // Number of inequality constraints 

    std::vector<Sparse> Bs_sub_;     // Equality constraint matrix
    double *b_sub_;                            // Right-hand-side equality constraints
    int mB_sub_;          

    Sparse sTmp_sub_;
    double* pdTmp_sub_;
    double *pdTmpLinear_;


    /* Triangle Inequalitie variables */
    TriCuts Cuts_;
    TriMap  Map_;
    TriHeap Heap_;
    
    std::vector<BiqTriInequality> container;


    int nIneq_;      // rename to nCuts_? 

    /* L-BFGS-B Data */
    double gradInorm_;
    double gradEnorm_;
    double *g_;
    double f_;
    double *y_;
    double *binf_;
    double *bsup_;
    int *nbd_;
    double *wa_;
    int *iwa_;

    /* dsyevr Data */
    int nEigenVal_;
    double *W_;
    double *Z_;
    int *ISUPPZ_;
    double *DWORK_; 
    int *IWORK_;
    int nDWORK_;
    int nIWORK_;
    int M_;

    std::vector<int> viSolution_1_;
    std::vector<int> viSolution_2_;

    std::vector<double> vdFracSol_;

    BiqModel(const BiqModel&);
    BiqModel& operator=(const BiqModel&);
    
public:
    BiqModel(){InitEmptyModel();};
    BiqModel(int nVar, bool max_problem,
             double *Q, Sparse Qs,
             std::vector<Sparse> As, double *a,
             std::vector<Sparse> Bs, double *b);
    ~BiqModel();
  /** Create the root node. Default: do nothing */
  virtual AlpsTreeNode * createRoot(){
     return (new BiqTreeNode(this));
    }

    inline int getNVar() const { return nVar_;}
    inline void setNVar( int nVar) { nVar_ = nVar;}
    double GetObjective(){return floor(f_);}

    inline bool isMax() const { return max_problem_;}
    virtual AlpsKnowledge * decode(AlpsEncoded & encode) const;

    double SDPbound(std::vector<BiqVarStatus> vbiqVarStatus, bool bRoot);

    void CreateSubProblem(std::vector<BiqVarStatus> vbiqVarStatus);

    double primalHeuristic();

    double GWheuristic(int nPlanes,  std::vector<BiqVarStatus> vbiqVarStatus);
    
    
    std::vector<double> GetFractionalSolution(std::vector<BiqVarStatus> vbiqVarStatus);
private:

    void AddDiagCons();

    void AllocSubCons();

    int CallLBFGSB(double dAlpha, double dTol, int &nbit);

    void ProjSDP();

    void sim(double alpha);

    void A(int mode, double alpha);

    void B(int mode, double alpha);
     
    void BuildConstraints(int nRows,
                          double *pdRHSsource, std::vector<Sparse> sMatSource,
                          double *pdRHSdest,   std::vector<Sparse> &sMatdest,
                          std::vector<BiqVarStatus> vbiqVarStatus, int nFixed);

    double GetSubMatrixSparse(Sparse sSourceMat, std::vector<BiqVarStatus> vbiqVarStatus, int & nnzAdded, int nFixed);

    void BuildObjective(std::vector<BiqVarStatus> vbiqVarStatus, int nFixed);

    void GetSubMatrix(std::vector<BiqVarStatus> vbiqVarStatus, int nFixed);

    double GetConstant(Sparse &sMat, std::vector<BiqVarStatus> vbiqVarStatus); 

    void GetLinear(Sparse &sSource, std::vector<BiqVarStatus> vbiqVarStatus, int nFixed);

    int GetOffset(std::vector<BiqVarStatus> vbiqVarStatus);

    double UpdateInequalities(int &nAdded, int &nSubtracted);

    double GetViolatedCuts();

    double EvalInequalities(TriType triType, int ii, int jj, int kk);

    void PrintBoundingTable(int iter, int nBit, int nAdded, int nSubtracted, double dAlpha, double dTol, double dMinAllIneq, double dGap/*double dTime*/);

    bool isFeasibleSolution(std::vector<int> solution);

    double EvalSolution(std::vector<int> solution);

    void UpdateSol(double dVal, std::vector<int> solution);

    void freeData(int *& data);

    void freeData(double *& data);

    void InitEmptyModel();
    
    void SetConSparseSize();

    bool pruneTest(double bound);
};

#endif
