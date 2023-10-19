#ifndef BiqNodeDesc_h_
#define BiqNodeDesc_h_

#include "AlpsNodeDesc.h"


enum BiqVarStatus {
    BiqVarFree = 0,
    BiqVarFixedToOne,
    BiqVarFixedToZero
};

class BiqModel;

/*!
  This is extened from the base class AlpsNodeDesc.
  The bellow is from AlpsNodeDesc.h:
  
  This is an abstract base class for subproblem data to be stored in a tree
  node. An instance of this class is a member in #AlpsTreeNode.

  #AlpsTreeNode keeps data related to node's position within the tree and
  #AlpsNodeDesc holds data directly related to the corresponding
  subproblem. This design is prefered due to its simplicity and convenience of
  separating subproblem data from tree search related data.

 */

class BiqNodeDesc :  public AlpsNodeDesc 
{
private:
    /* data */
    BiqModel * model_;
    BiqVarStatus * varStatus_;


    // Disable copy constructor
    BiqNodeDesc(BiqNodeDesc const & other);
    BiqNodeDesc & operator=(BiqNodeDesc const & rhs);
public:

    BiqNodeDesc(BiqModel * model);
    BiqNodeDesc(BiqModel * model, BiqVarStatus *& st);

    virtual ~BiqNodeDesc();

    void setVarStatus(const int i, const BiqVarStatus status);
    BiqVarStatus getVarStatus(const int i) { return varStatus_[i];};
    BiqVarStatus const * getVarStati() const;

    virtual AlpsKnowledge * decode(AlpsEncoded & encode) const;

    BiqModel* model() { return model_; };
    BiqModel const * model() const { return model_; };


    void bound(int & iBound, double * pdSol);
};





#endif