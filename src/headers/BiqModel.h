#ifndef BiqModel_h_
#define BiqModel_h_

#include "AlpsKnowledge.h"
#include "AlpsModel.h"
#include "BiqTreeNode.h"
#include "BiqNodeDesc.h"
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
                                 
    int max_problem_;                      // Objective sense 

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
    BiqModel(const BiqModel&);
    BiqModel& operator=(const BiqModel&);
    
public:
    BiqModel() {}
    BiqModel(int nVar, double *Q, int max_problem,
             std::vector<Sparse> As, double *a,
             std::vector<Sparse> Bs, double *b);
    ~BiqModel();
  /** Create the root node. Default: do nothing */
  virtual AlpsTreeNode * createRoot(){
     return (new BiqTreeNode(this));
    }

    inline int getNVar() const { return nVar_;}
    inline void setNVar( int nVar) { nVar_ = nVar;}

    virtual AlpsKnowledge * decode(AlpsEncoded & encode) const;

    double SDPbound();

    void CreateSubProblem();

private:

    void AddDiagCons();

    void AllocSubCons();

    int CallLBFGSB(double dAlpha, double dTol, int &nbit);

    void ProjSDP();

    void sim(double alpha);

    void A(int mode, double alpha);

    void B(int mode, double alpha);
    
    bool Prune();
    
    void BuildConstraints(int nRows,
                          double *pdRHSsource, std::vector<Sparse> sMatSource,
                          double *pdRHSdest,   std::vector<Sparse> &sMatdest,
                          int *piSol, std::vector<BiqVarStatus> vbiqVarStatus, int nFixed);

    double GetSubMatrixSparse(Sparse sSourceMat, std::vector<BiqVarStatus> vbiqVarStatus, int & nnzAdded, int nFixed);

    void BuildObjective(int *piSol, std::vector<BiqVarStatus> vbiqVarStatus, int nFixed);

    void GetSubMatrix(std::vector<BiqVarStatus> vbiqVarStatus, int nFixed);

    double GetConstant(Sparse &sMat, int *piSol, std::vector<BiqVarStatus> vbiqVarStatus); 

    void GetLinear(Sparse &sSource, int *piSol, std::vector<BiqVarStatus> vbiqVarStatus, int nFixed);

    int GetOffset(std::vector<BiqVarStatus> vbiqVarStatus);

    double UpdateInequalities(int &nAdded, int &nSubtracted);

    double GetViolatedCuts();

    double EvalInequalities(TriType triType, int ii, int jj, int kk);

    void PrintBoundingTable(int iter, int nBit, int nAdded, int nSubtracted, double dAlpha, double dTol, double dMinAllIneq/*double dTime*/);
};

#endif
