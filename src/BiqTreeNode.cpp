

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
    desc_ = new BiqNodeDesc(model); // make biq
}

BiqTreeNode::BiqTreeNode(BiqNodeDesc *& desc)
{
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


/// @brief 
/// @param isRoot 
/// @param rampUp 
/// @return 
int BiqTreeNode::process(bool isRoot, bool rampUp)
{
    bool bFoundSolution = false;
    bool bmaxProblem;

    // get a pointer to out desc class
    BiqNodeDesc * desc = dynamic_cast<BiqNodeDesc *>(desc_);
    // BiqModel *model=dynamic_cast<BiqModel *>(broker()->getModel()); // maybe do this?
    // get a pointer to out model class
    BiqModel * model = dynamic_cast<BiqModel *>(desc->model());
    
    bmaxProblem = model->isMax();

    // get the best solution so far (for parallel, it is the incumbent)
    int bestVal = static_cast<int>(broker()->getIncumbentValue());

    // build sub prob
    const std::vector<BiqVarStatus> biqVarStatus = desc->getVarStati();

    model->CreateSubProblem(biqVarStatus);

    // check all variables ar fixed..
    int nFixed = 0;
    for(int i = 0; i < biqVarStatus.end()-biqVarStatus.begin();++i)
    {
        if(biqVarStatus.at(i) != BiqVarFree) ++nFixed;
    }
    //if(nFixed == model->getNVar)


    // call bounding 
    double valRelax = model->SDPbound();

    std::vector<double> vdFracSol = model->GetFractionalSolution(biqVarStatus);

    // get the best solution so far (for parallel, it is the incumbent)
    bestVal = static_cast<int>(broker()->getIncumbentValue());
    if(bmaxProblem)
    {
        bestVal = -bestVal;
    }

    if(
        ( bmaxProblem && static_cast<int>(floor(valRelax)) <= bestVal) || 
        (!bmaxProblem && static_cast<int>(ceil(valRelax))  >= bestVal)
      )
    {
        setStatus(AlpsNodeStatusFathomed);
    }
    else
    {
        // determine variable to branch on using vdFracSol
        SetBranchingVariable(vdFracSol, biqVarStatus);
        //setStatus(AlpsNodeStatusPregnant);
        setStatus(AlpsNodeStatusPregnant);
    }

    //setStatus(AlpsNodeStatusFathomed);

    return bFoundSolution;
}


std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > BiqTreeNode::branch()
{
    // get a pointer to out desc class
    BiqNodeDesc * desc = dynamic_cast<BiqNodeDesc *>(desc_);
    // get a pointer to out model class
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


void BiqTreeNode::SetBranchingVariable(std::vector<double> fracSol, std::vector<BiqVarStatus> varStatus)
{
    double dMaxVal = -INFINITY;
    branchOn_ = 0;

   for(int i = 0; i < fracSol.size(); ++i)
   {
        if(varStatus.at(i) == BiqVarFree && fabs(0.5 - fracSol.at(i)) > dMaxVal)
        {
            branchOn_ = i;
            dMaxVal = fabs(0.5 - fracSol.at(i));
        }
   }


}