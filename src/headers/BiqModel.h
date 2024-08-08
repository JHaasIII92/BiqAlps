#ifndef BiqModel_h_
#define BiqModel_h_

#include "AlpsKnowledge.h"
#include "AlpsModel.h"
#include "CoinMessageHandler.hpp"
#include "BiqTreeNode.h"
#include "BiqNodeDesc.h"
#include "BiqSolution.h"
#include "BiqUtil.h"
#include "BiqParams.h"
#include "BiqMessage.h"
#include "math.h"
#include<iostream>
#include <utility>
#include <list>
#include <map>
#include <unordered_map>
#include <tuple>
#include <queue>
#include <string>
#include <fstream>

class BiqModel: public AlpsModel
{
private:
    /* data */
    int nVar_;                            // Number of variables
    bool max_problem_;                      // Objective sense  
    double *Q_;                           // Dense Objective matrix   Q_[i + j*nVar_] to get the i,j
    Sparse Qs_;
    int N_;
    
    // array of array
    std::vector<Sparse> As_;     // Inequality constraint matrix   
    double *a_;                            // Right-hand-side inequality constraints
    int mA_;                               // Number of inequality constraints 
    std::vector<int> iInequalityIsLinear_;

    std::vector<Sparse> Bs_;     // Equality constraint matrix
    double *b_;                            // Right-hand-side equality constraints
    double *b_original_;                            // Right-hand-side equality constraints
    int mB_original_;                               // Number of equality constraints from original data
    int mB_;                               // Number of equality constraints
    
    std::vector<int> iEqualityIsLinear_;

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

    // past knowlage
    int iPastBest_;
    bool bHasPastBest_;

    // storage space
    std::vector<int> viSolution_1_;
    std::vector<int> viSolution_2_;

    std::vector<double> vdFracSol_;

    BiqModel(const BiqModel&);
    BiqModel& operator=(const BiqModel&);


    BiqParams *BiqPar_;

    /// Message handler
    CoinMessageHandler * handler_;
    /// Biq messages
    CoinMessages messages_;
    
    bool defaultHandler_;
public:
    BiqModel(){InitEmptyModel();};
    ~BiqModel();
  /** Create the root node. Default: do nothing */
  virtual AlpsTreeNode * createRoot(){
     return (new BiqTreeNode(this));
    }

    void setPastBest(int iPastBest) {iPastBest_ = iPastBest;};
    int getPastBest() {return iPastBest_;};
 
    void setHasPastBest(bool bHasPastBest) {bHasPastBest_ = bHasPastBest;};
    bool hasPastBest() {return bHasPastBest_;};

    inline int getNVar() const { return nVar_;}
    inline void setNVar( int nVar) { nVar_ = nVar;}

    double GetObjective(){return floor(f_);}

    inline bool isMax() const { return max_problem_;}
    
    bool hasBestVal();
    int getBestVal();

    /// Get encode from AlpsModel.
    using AlpsKnowledge::encode;
    /// Encode this into the given AlpsEncoded object.
    virtual AlpsReturnStatus encode(AlpsEncoded * encoded) const;
    /// Decode the given AlpsEncoded object into this.
    virtual AlpsReturnStatus decodeToSelf(AlpsEncoded & encoded);
    virtual AlpsKnowledge * decode(AlpsEncoded & encoded) const;

    double SDPbound(std::vector<BiqVarStatus> vbiqVarStatus, bool bRoot);

    void CreateSubProblem(std::vector<BiqVarStatus> vbiqVarStatus);

    double primalHeuristic();

    double GWheuristic(int nPlanes,  std::vector<BiqVarStatus> vbiqVarStatus);
    
    int GetOffset(std::vector<BiqVarStatus> vbiqVarStatus);
    
    //double primalHeuristicKC();
    void KCheuristic(std::vector<BiqVarStatus> vbiqVarStatus);
    
    
    void NieveUpdateHeuristic(std::vector<BiqVarStatus> vbiqVarStatus);

    std::vector<double> GetFractionalSolution(std::vector<BiqVarStatus> vbiqVarStatus);

    virtual void readInstance(const char* dataFile);

    BiqParams *BiqPar() { return BiqPar_; }

    void readParameters(const int argnum, const char * const * arglist);

    /**@name Message handling */
    //@{
    /// Pass in Message handler (not deleted at end)
    void passInMessageHandler(CoinMessageHandler * handler)
        {
            if (defaultHandler_) {
                delete handler_;
                handler_ = NULL;
            }
            defaultHandler_ = false;
            handler_ = handler;
        }

    /// Set language
    void newLanguage(CoinMessages::Language language)
        { messages_ = BiqMessage(language); }
    void setLanguage(CoinMessages::Language language)
        { newLanguage(language); }
    /// Return handler
    CoinMessageHandler * messageHandler() const
        { return handler_; }
    /// Return messages
    CoinMessages messages()
        { return messages_; }
    /// Return pointer to messages
    CoinMessages * messagesPointer()
        { return &messages_; }
    //@}


private:

    void AddDiagCons();

    void AddProdCons();

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

    double UpdateInequalities(int &nAdded, int &nSubtracted);

    double GetViolatedCuts();

    inline double EvalInequalities(TriType triType, int ii, int jj, int kk);

    void PrintBoundingTable(int iter, int nBit, int nAdded, int nSubtracted, double dAlpha, double dTol, double dMinAllIneq, double dGap/*double dTime*/);

    bool isFeasibleSolution(std::vector<int> solution);

    double EvalSolution(std::vector<int> & solution);

    bool UpdateSol(double dVal, std::vector<int> solution);

    void freeData(int *& data);

    void freeData(double *& data);

    void InitEmptyModel();

    void InitModel();
    
    void SetConSparseSize();

    bool pruneTest(double bound);

    void testEncodeDecode();

    void addProvidedSol();
};

#endif
