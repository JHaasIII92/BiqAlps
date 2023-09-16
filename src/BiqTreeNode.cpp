

#include <iostream>
#include <utility>

//#include "CoinFloat.h"
#include "CoinUtility.hpp"

#include "AlpsKnowledgeBroker.h"
#include "AlpsKnowledge.h"

#include "BiqTreeNode.h"
#include "BiqNodeDesc.h"
#include "BiqSolution.h"
#include "BiqModel.h"

BiqTreeNode::BiqTreeNode(BiqModel * model)
{
    desc_ = new BiqNodeDesc(model); // make biq
}

BiqTreeNode::BiqTreeNode(BiqNodeDesc *& desc)
{
    desc_ = desc;
    desc = 0;
}

AlpsTreeNode* BiqTreeNode::createNewTreeNode(AlpsNodeDesc*& desc) const 
{
    // cast to BiqTreeNode
    BiqNodeDesc * d = dynamic_cast<BiqNodeDesc*>(desc);
    // call the copy constructor
    BiqTreeNode * node = new BiqTreeNode(d);
    // TODO figure out why desc = 0
    desc = 0;
    return node;
}



int BiqTreeNode::process(bool isRoot, bool rampUp)
{
    return 0;
}


std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > BiqTreeNode::branch()
{

    std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > newNodes;

    return newNodes;

}


AlpsKnowledge * BiqTreeNode::decode(AlpsEncoded & encoded) const{
    BiqTreeNode * nn = new BiqTreeNode(dynamic_cast<BiqNodeDesc*>(desc_)->model());
    return nn;
}