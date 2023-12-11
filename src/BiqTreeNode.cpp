

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

    const std::vector<BiqVarStatus> biqVarStatus = desc->getVarStati();
    
    if(isRoot)
    {
        // try the brute force heuristic 
        model->primalHeuristic();
    }

    // get the best solution so far (for parallel, it is the incumbent)
    int bestVal = static_cast<int>(broker()->getIncumbentValue());

    // build sub prob
    model->CreateSubProblem(biqVarStatus);

    // check all variables ar fixed..
    int nFixed = 0;
    for(int i = 0; i < biqVarStatus.size();++i)
    {
        if(biqVarStatus.at(i) != BiqVarFree) ++nFixed;
    }
    //if(nFixed == model->getNVar)

    double valRelax = model->SDPbound();

    std::vector<double> vdFracSol = model->GetFractionalSolution(biqVarStatus);

    // call bounding 
    for(auto &it: biqVarStatus)
    {
        switch(it) {
            case BiqVarFixedToZero:
                std::printf("0   ");
                break;
            case BiqVarFixedToOne:
                std::printf("1   ");
                break;
            case BiqVarFree:
                std::printf("F   ");
                break;
            default:
                std::printf("E   ");
                break;
        }
    }
    std::printf("\n");

    // print vdFracSol vector in one line
    for(auto &it: vdFracSol)
    {
        std::printf("%3.1f ", it);
    }
    std::printf("\n");
    // get the best solution so far (for parallel, it is the incumbent)
    bestVal = static_cast<int>(broker()->getIncumbentValue());
    if(bmaxProblem)
    {
        bestVal = -bestVal;
    }
    std::printf("BestVal = %d\n",bestVal);

    if(
        ( bmaxProblem && static_cast<int>(floor(valRelax)) <= bestVal) || 
        (!bmaxProblem && static_cast<int>(ceil(valRelax))  >= bestVal) ||
        nFixed ==  model->getNVar()
      )
    {
        std::printf("This node is being fathomed\n");
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

    //std::printf("x[%d] = 0\n",branchOn_);
    newStatiRight.at(branchOn_) = BiqVarFixedToZero;
    AlpsNodeDesc *descRight = new BiqNodeDesc(model, newStatiRight);
    newNodes.push_back(CoinMakeTriple(descRight, AlpsNodeStatusCandidate, 0.0));

    //std::printf("x[%d] = 1\n",branchOn_);
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
    double dMinVal = INFINITY;
    branchOn_ = 0;

   for(int i = 0; i < fracSol.size(); ++i)
   {
        /*
        // Branch on the variable x[i] that has the least fractional value
        //std::printf("%d =>\t %f \t %f\n", i, fracSol.at(i), fabs(0.5 - fracSol.at(i)));
        if(varStatus.at(i) == BiqVarFree && fabs(0.5 - fracSol.at(i)) > dMaxVal)
        {
            branchOn_ = i;
            dMaxVal = fabs(0.5 - fracSol.at(i));
        }
        */

        // Branch on the variable x[i] that has the most fractional value
        if(varStatus.at(i) == BiqVarFree && fabs(0.5 - fracSol.at(i)) < dMinVal)
        {
            branchOn_ = i;
            dMinVal = fabs(0.5 - fracSol.at(i));
        }

   }
}