#ifndef BiqModel_h_
#define BiqModel_h_

#include "AlpsKnowledge.h"
#include "AlpsModel.h"
#include "BiqTreeNode.h"

class BiqModel: public AlpsModel
{
private:
    /* data */

    // Q[i + j*nVar] to get i,j entry of the matrix
    //double *Q;
    //int nVar;

    
    BiqModel(const BiqModel&);
    BiqModel& operator=(const BiqModel&);
    
public:

  /** Create the root node. Default: do nothing */
  virtual AlpsTreeNode * createRoot(){
     return (new BiqTreeNode(this));
    }

    virtual AlpsKnowledge * decode(AlpsEncoded & encode) const;

};

#endif
