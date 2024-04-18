#ifndef BiqTreeNode_h_
#define BiqTreeNode_h_

#include <utility>

#include "AlpsTreeNode.h"


#include "AlpsKnowledgeBroker.h"
#include "AlpsKnowledge.h"

#include "BiqUtil.h"

class BiqModel;
class BiqNodeDesc;


//#############################################################################
/** BiqTree node is extended from AlpsTreeNode.
 
    This class holds one node of the search tree. Note that the generic search
    procedure doesn't know anything about the nodes in the
    tree other than their index, lower bound, etc. Other application-specific
    data is contained in derived classes, but is not needed for the basic
    operation of the search tree.*/
//#############################################################################

class BiqTreeNode : public AlpsTreeNode
{

private:

    int branchOn_;
    bool bCloseToZero_;
    
public:

   BiqTreeNode( BiqModel * model);
   BiqTreeNode( BiqNodeDesc *& desc);
   ~BiqTreeNode();
        /** The purpose of this function is be able to create the children of
        a node after branching. */
    /* FIXME: I think that we probably want the argument to be a diff'd
       description, but this is open to debate. Maybe we should
       overload this method and have a version that creates a diff'd
       description, as well as one that creates an explicit
       description. */
    virtual AlpsTreeNode* createNewTreeNode(AlpsNodeDesc*& desc) const;


        /**
     *  Perform the processing of the node. For branch and bound, this would
     *  mean performing the bounding operation. The minimum requirement for
     *  this method is that it set the status of the node.
     *  <ul>
     *  <li> active_ flag is set
     *  <li> the description is explicit
     *  <li> status_ is not pregnant nor fathomed
     *  </ul>
     *
     *  The status of the node is set to one of the following (and children may
     *  be internally stored, etc.):
     *  <ul>
     *  <li> [candidate:] the node is put back into the priority queue for
     *       later processing
     *  <li> [evaluated:] the node is put back into the priority queue for
     *       later branching
     *  <li> [pregnant:] the node has already decided how to branch and will
     *       give birth to the children when instructed with the \c branch()
     *       function which will be the next method invoked for the node.
     *  <li> [branched:] the node has tried to create new nodes by calling
     *       <code>createChildren<\code>.
     *  <li> [fathomed:] the node will not be further considered in this phase.
     *  </ul>
     *  \param isRoot  Indicate if this node is a root of a subtree.
     *  \param rampUp  Indicate if it is in ramp up period. Only useful
     *                 for parallel code.
     * Currently, the return code of this method is not used.
     */
    virtual int process(bool isRoot = false, bool rampUp = false);

        /** This method must be invoked on a \c pregnant node (which has all the
        information needed to create the children) and should create the
        children in the data structure of the node. The stati of the children
        can be any of the ones \c process() can return.  The third component
        is the \c priority.

        NOTE: This method may be almost empty if the descriptions of the
        children were created when the node became pregnant. In that case this
        routine is pretty much just copying data. */
    virtual std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > 
        branch();
    

    using AlpsTreeNode::encode;
    /// Encode this into the given AlpsEncoded object.
    virtual AlpsReturnStatus encode(AlpsEncoded * encoded) const;
    /// Decode the given AlpsEncoded object into this.
    virtual AlpsReturnStatus decodeToSelf(AlpsEncoded & encoded);
    virtual AlpsKnowledge * decode(AlpsEncoded & encoded) const;

    /*
    
    
    */
private:
    BiqTreeNode(BiqTreeNode const &);
    BiqTreeNode & operator=(BiqTreeNode const &);

    void SetBranchingVariable(std::vector<double> fracSol, std::vector<BiqVarStatus> varStatus);
};

#endif