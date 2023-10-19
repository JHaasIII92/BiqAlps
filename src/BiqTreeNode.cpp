

#include <iostream>
#include <utility>

//#include "CoinFloat.h"
#include "CoinUtility.hpp"

#include "AlpsKnowledgeBroker.h"
#include "AlpsKnowledge.h"

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
    std::cout << "BiqTreeNode::~BiqTreeNode() ..." << std::endl;
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
    std::cout << "BiqTreeNode::process index_ " << index_ << " ..." << std::endl;
    bool bFoundSolution = false;
    BiqNodeDesc * desc = dynamic_cast<BiqNodeDesc *>(desc_);
    const BiqModel * model = dynamic_cast<BiqModel *>(desc->model());
    // get the status from our desc
    const BiqVarStatus* bvStati = desc->getVarStati();
    int nVar = model->getNVar();
    // call the bounding method
    int iBound;
    double *pdSol = new double[nVar];
    desc->bound(iBound, pdSol);

    int objVal = 5;
    int * piSol = new int[nVar];
    std::fill_n(piSol, nVar, 0);
    BiqSolution* biqSol = new BiqSolution( dynamic_cast<BiqModel *>(desc->model()), piSol, objVal);
    broker()->addKnowledge(AlpsKnowledgeTypeSolution, biqSol, objVal);
    //  since Alps consideres sols with lower quality values better.
    double bestVal;
    bestVal = static_cast<double>(broker()->getIncumbentValue());

    if(iBound < bestVal)
    {
        // pregnent
        branchOn_ = 0;
        setStatus(AlpsNodeStatusPregnant);
    }
    else
    {
        setStatus(AlpsNodeStatusFathomed);
    }
    return bFoundSolution;
}


std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > BiqTreeNode::branch()
{
    std::cout << "BiqTreeNode::branch index_ " << index_ << std::endl;
    BiqNodeDesc * desc = dynamic_cast<BiqNodeDesc *>(desc_);
    BiqModel * model = dynamic_cast<BiqModel *>(desc->model());
    std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > newNodes;
    int nVar = model->getNVar();
    const BiqVarStatus * oldStati = desc->getVarStati();
    BiqVarStatus *newStatiLeft = new BiqVarStatus[nVar];
    BiqVarStatus *newStatiRight = new BiqVarStatus[nVar];
    std::copy(oldStati, oldStati+nVar, newStatiLeft);
    std::copy(oldStati, oldStati+nVar, newStatiRight);
    // First Triple
    // AlpsNodeDesc:
    // make a new desc
    // fix the branchOn_ variable of varStatus_ BiqVarFixedToZero 
    // AlpsNodeStatus: use AlpsNodeStatusCandidate so it can get to BiqTreeNode::process
    newStatiRight[branchOn_] = BiqVarFixedToZero;
    std::cout << " BiqTreeNode::branch ... newStatiRight[" << branchOn_ <<"] = " << newStatiRight[branchOn_] << " ..." << std::endl;
    AlpsNodeDesc *descRight = new BiqNodeDesc(model, newStatiRight);
        std::cout << " BiqTreeNode::branch ... newStatiRight[" << branchOn_ <<"] = " << dynamic_cast<BiqNodeDesc*>(descRight)->getVarStatus(branchOn_) << " ..." << std::endl;
    newNodes.push_back(CoinMakeTriple(descRight, AlpsNodeStatusCandidate, 1.0));

    // Second Triple
    // Same but use BiqVarFixedToOne
    newStatiLeft[branchOn_] = BiqVarFixedToOne;
    AlpsNodeDesc *descLeft = new BiqNodeDesc(model, newStatiLeft);
    newNodes.push_back(CoinMakeTriple(descLeft, AlpsNodeStatusCandidate, 1.0)); 

    return newNodes;
}


AlpsKnowledge * BiqTreeNode::decode(AlpsEncoded & encoded) const{
    
    std::cout << "BiqTreeNode::decode" << std::endl;
    BiqTreeNode * nn = new BiqTreeNode(dynamic_cast<BiqNodeDesc*>(desc_)->model());
    return nn;
}