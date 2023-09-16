#ifndef BiqNodeDesc_h_
#define BiqNodeDesc_h_

#include "AlpsNodeDesc.h"

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
public:

    BiqNodeDesc(BiqModel * model);
    virtual AlpsKnowledge * decode(AlpsEncoded & encode) const;

    BiqModel* model() { return model_; };
};





#endif