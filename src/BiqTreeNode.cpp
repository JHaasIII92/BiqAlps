

#include <iostream>
#include <utility>

//#include "CoinFloat.h"
#include "CoinUtility.hpp"


#include "headers/BiqTreeNode.h"
#include "headers/BiqNodeDesc.h"
#include "headers/BiqSolution.h"
#include "headers/BiqModel.h"


BiqTreeNode::BiqTreeNode(BiqModel * model)
{    
    std::cout << "BiqTreeNode(BiqModel * model) ..." << std::endl;
    desc_ = new BiqNodeDesc(model); // make biq
}

BiqTreeNode::BiqTreeNode(BiqNodeDesc *& desc)
{
    std::cout << "BiqTreeNode::BiqTreeNode(BiqNodeDesc *& desc) ..." << std::endl;
    std::cout << "BiqTreeNode::BiqTreeNode ... desc->varStatus_[0] = " << desc->getVarStatus(0) << " ... " << std::endl;
    desc_ = desc;
    desc = 0;
}

BiqTreeNode::~BiqTreeNode()
{

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
    setStatus(AlpsNodeStatusFathomed);
    return true;
}


std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > BiqTreeNode::branch()
{
    std::cout << "BiqTreeNode::branch index_ " << index_ << std::endl;

    BiqNodeDesc * desc = dynamic_cast<BiqNodeDesc *>(desc_);

    BiqModel * model = dynamic_cast<BiqModel *>(desc->model());

    std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > newNodes;

    int nVar = model->getNVar();
    const std::vector<BiqVarStatus> oldStati = desc->getVarStati();
    std::vector<BiqVarStatus> newStatiLeft = oldStati;
    std::vector<BiqVarStatus> newStatiRight = oldStati;

    newStatiRight.at(branchOn_) = BiqVarFixedToZero;
    AlpsNodeDesc *descRight = new BiqNodeDesc(model, newStatiRight);
    newNodes.push_back(CoinMakeTriple(descRight, AlpsNodeStatusCandidate, 0.0));

    newStatiLeft.at(branchOn_)  = BiqVarFixedToOne;
    AlpsNodeDesc *descLeft = new BiqNodeDesc(model, newStatiLeft);
    newNodes.push_back(CoinMakeTriple(descLeft, AlpsNodeStatusCandidate, 1.0)); 

    return newNodes;
}


AlpsKnowledge * BiqTreeNode::decode(AlpsEncoded & encoded) const{
    
    std::printf("AlpsKnowledge * BiqTreeNode::decode\n");
    std::cout << "BiqTreeNode::decode" << std::endl;
    BiqTreeNode * nn = new BiqTreeNode(dynamic_cast<BiqNodeDesc*>(desc_)->model());
    return nn;
}