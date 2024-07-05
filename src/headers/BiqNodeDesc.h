#ifndef BiqNodeDesc_h_
#define BiqNodeDesc_h_

#include "AlpsNodeDesc.h"
#include "BiqModel.h"
#include "BiqUtil.h"
#include <iostream>



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
    std::vector<BiqVarStatus> varStatus_;
    double dQuality_;
    int iParentsBest_;
    int iBranchedOn_;
    bool bHasBest_;
    // Disable copy constructor
    BiqNodeDesc(BiqNodeDesc const & other);
    BiqNodeDesc & operator=(BiqNodeDesc const & rhs);
public:

    BiqNodeDesc(BiqModel * model);
    BiqNodeDesc(BiqModel * model, std::vector<BiqVarStatus> & st);

    void setQuality(double val){dQuality_ = val;};
    double getQuality(){return dQuality_;}
    
    void setBranchedOn(int val){iBranchedOn_ = val;};
    int getBranchedOn(){return iBranchedOn_;}

    void setParentsBest(int val){iParentsBest_ = val;};
    int getParentsBest(){return iParentsBest_;}

    void setHasBest(int val){bHasBest_ = val;};
    bool getHasBest(){return bHasBest_;}

    virtual ~BiqNodeDesc();

    void setVarStatus(const int i, const BiqVarStatus status);
    BiqVarStatus getVarStatus(const int i) { return varStatus_.at(i);};
    std::vector<BiqVarStatus> const getVarStati() const;

    /// Grab encode function from #AlpsKnowledge
    using AlpsKnowledge::encode;
    /// Pack this node description into the given #AlpsEncoded object.
    virtual AlpsReturnStatus encode(AlpsEncoded * encoded) const;
    /// Unpack fields from the given #AlpsEncoded object.
    virtual AlpsReturnStatus decodeToSelf(AlpsEncoded & encoded);
    virtual AlpsNodeDesc * decode(AlpsEncoded & encoded) const;

    BiqModel* model() { return model_; };
    BiqModel const * model() const { return model_; };


    void bound(int & iBound, std::vector<int> solution);
};





#endif