#ifndef BiqSubProb_h_
#define BiqSubProb_h_

#include "BiqUtil.h";


class BiqSubProb
{
private:
    
    /* Problem Data */
    double *Q_;
    BiqSparseTriple * Qs_;
    int nVar_;

    /* L-BFGS-B Data */
    double gradInorm_;
    double granEnorm_;
    double *g_;
    double f_;
    double *y_;
    double *binf_;
    double *bsup_;
    int *nbd;
    double *wa;
    int *iwa;

    /* dsyevr Data */
    int nEigenVal_;
    double *W_;
    double *Z_;
    int *ISUPPZ_;
    double *DWORK_; 
    int *IWORK_;
    int nDWORK_;
    int nIWORK_;


public:
    BiqSubProb(int nVar, int nEigenVal, int nIWORK, int nDWORk);
    ~BiqSubProb();
};










#endif