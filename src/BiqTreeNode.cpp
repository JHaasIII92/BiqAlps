

#include <iostream>
#include <utility>

//#include "CoinFloat.h"
#include "CoinUtility.hpp"


#include "headers/BiqTreeNode.h"
#include "headers/BiqNodeDesc.h"
#include "headers/BiqSolution.h"
#include "headers/BiqModel.h"
#include "headers/BiqParams.h"


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
    bool bHasBestVal = false;
    bool bHasPastBestVal;
    int bestVal;
    int pastBestVal;
    // get a pointer to out desc class
    BiqNodeDesc * desc = dynamic_cast<BiqNodeDesc *>(desc_);
    bHasPastBestVal = desc->getHasBest();
    // BiqModel *model=dynamic_cast<BiqModel *>(broker()->getModel()); // maybe do this?
    // get a pointer to out model class
    BiqModel * model = dynamic_cast<BiqModel *>(broker()->getModel());
    bmaxProblem = model->isMax();
    // get the 1, 0, free vector to pass to model
    const std::vector<BiqVarStatus> biqVarStatus = desc->getVarStati();
    
    // get the best solution so far (for parallel, it is the incumbent)
    bestVal = static_cast<int>(broker()->getIncumbentValue());
    bHasBestVal = broker()->hasKnowledge(AlpsKnowledgeTypeSolution);
    pastBestVal = desc->getParentsBest();
    // pass desc knowlage of best to model for bounding
    if(bHasPastBestVal)
    {
        model->setHasPastBest(true);
        model->setPastBest(pastBestVal);
    }

    // build sub prob
    model->CreateSubProblem(biqVarStatus);

    // check all variables ar fixed..
    int nFixed = 0;
    for(size_t i = 0; i < biqVarStatus.size();++i)
    {
        if(biqVarStatus.at(i) != BiqVarFree) ++nFixed;
    }

    double valRelax = model->SDPbound(biqVarStatus, isRoot);

    std::vector<double> vdFracSol = model->GetFractionalSolution(biqVarStatus);

    bestVal = model->getBestVal();
    bHasBestVal =  model->hasBestVal();

    if( bHasBestVal && (
        ( bmaxProblem && static_cast<int>(floor(valRelax)) <= bestVal) || 
        (!bmaxProblem && static_cast<int>(ceil(valRelax))  >= bestVal) || 
        nFixed ==  model->getNVar())
      )
    {
        //std::printf("This node is being fathomed \t valRelax: %f\t bestVal: %d  parent BV: %d\n",valRelax, bestVal, desc->getParentsBest());
        setStatus(AlpsNodeStatusFathomed);
    }
    else
    {
        //std::printf("This node is being Pregnant \t valRelax: %f\t bestVal: %d  parent BV: %d\n",valRelax, bestVal, desc->getParentsBest());
        //setStatus(AlpsNodeStatusPregnant);
        setStatus(AlpsNodeStatusPregnant);
    }

    //setStatus(AlpsNodeStatusFathomed);

    return bFoundSolution;
}


std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > BiqTreeNode::branch()
{
    double dQuality;
    // get a pointer to our desc class
    BiqNodeDesc * desc = dynamic_cast<BiqNodeDesc *>(desc_);
    // get a pointer to our model class
    BiqModel * model = dynamic_cast<BiqModel *>(broker()->getModel());
    
    const std::vector<BiqVarStatus> biqVarStatus = desc->getVarStati();
    std::vector<double> vdFracSol = model->GetFractionalSolution(biqVarStatus);
    // determine variable to branch on using vdFracSol
    SetBranchingVariable(vdFracSol, biqVarStatus);

    if(model->isMax()) 
    {
        dQuality = model->GetObjective(); // = -model->f_
    }
    else
    {
        dQuality = -model->GetObjective();
    }
    
    int bestVal = model->getBestVal();

    //
    bool bHasBrokerBestVal = broker()->hasKnowledge(AlpsKnowledgeTypeSolution);
    bool bHasPastBestVal = desc->getHasBest();
    bool bHasBestVal = bHasBrokerBestVal || bHasPastBestVal;




    // compare models best with past best
    std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > newNodes;

    const std::vector<BiqVarStatus> oldStati = desc->getVarStati();
    std::vector<BiqVarStatus> newStatiLeft = oldStati;
    std::vector<BiqVarStatus> newStatiRight = oldStati;

    //std::printf("BestVal: %d  Has Best: %d  x[%d] = 0\n",bestVal,bHasBestVal,branchOn_);

    newStatiRight.at(branchOn_) = BiqVarFixedToZero;
    BiqNodeDesc *descRight = new BiqNodeDesc(model, newStatiRight);
    descRight->setBroker(broker());
    descRight->setHasBest(bHasBestVal);
    descRight->setParentsBest(bestVal);
    descRight->setQuality(dQuality);
    descRight->setBranchedOn(branchOn_);
    newNodes.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc *>(descRight), AlpsNodeStatusCandidate, dQuality));

    //std::printf("x[%d] = 1\n",branchOn_);
    newStatiLeft.at(branchOn_)  = BiqVarFixedToOne;
    BiqNodeDesc *descLeft = new BiqNodeDesc(model, newStatiLeft);
    descLeft->setBroker(broker());
    descLeft->setHasBest(bHasBestVal);
    descLeft->setParentsBest(bestVal);
    descLeft->setQuality(dQuality);
    descLeft->setBranchedOn(branchOn_);
    newNodes.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc *>(descLeft), AlpsNodeStatusCandidate, dQuality)); 


    return newNodes;
}

void BiqTreeNode::SetBranchingVariable(std::vector<double> fracSol, std::vector<BiqVarStatus> varStatus)
{
    BiqModel * model = dynamic_cast<BiqModel *>(broker()->getModel());
    double dBestVal = -INFINITY;
    int branchingStrategy = model->BiqPar()->entry(BiqParams::branchingStrategy);

    branchOn_ = 0;

   // branching stratagy...
   // find the least fractional 
   // then flag for prioity ....
   for(size_t i = 0; i < fracSol.size(); ++i)
   {
        
        if(branchingStrategy == CLOSEST_TO_ONE)
        {
            // Branch on the variable x[i] that has the least fractional value
            dBestVal = 0;
            if(varStatus.at(i) == BiqVarFree && fracSol.at(i) > dBestVal)
            {
                branchOn_ = i;
                dBestVal = fracSol.at(i);
            }
        }
        else if(branchingStrategy == LEAST_FRACTIONAL)
        {
            // Branch on the variable x[i] that has the least fractional value
            dBestVal = 0.5;
            if(varStatus.at(i) == BiqVarFree && fabs(0.5 - fracSol.at(i)) > dBestVal)
            {
                branchOn_ = i;
                dBestVal = fabs(fracSol.at(i));
            }
        }
        else
        {
            // default MOST_FRACTIONAL
            // Branch on the variable x[i] that has the most fractional value
            dBestVal = 1;
            if(varStatus.at(i) == BiqVarFree && fabs(0.5 - fracSol.at(i)) < dBestVal)
            {
                //std::printf("frac sol => %f\n", fracSol.at(i));
                branchOn_ = i;
                dBestVal = fabs(fracSol.at(i));
            }
        }
   }
   //printf("fracSol.at(%d) = %f\n", branchOn_, fracSol.at(branchOn_));
}


AlpsReturnStatus BiqTreeNode::encode(AlpsEncoded * encoded) const {
  AlpsReturnStatus status;
  // Alps part
  status = AlpsTreeNode::encode(encoded);
  
  // Biq part
  encoded->writeRep(branchOn_);
  status = dynamic_cast<BiqNodeDesc*>(desc_)->encode(encoded);
  assert(status==AlpsReturnStatusOk);

  return status;
}

/// Decode the given AlpsEncoded object into this.
AlpsReturnStatus BiqTreeNode::decodeToSelf(AlpsEncoded & encoded) {
    AlpsReturnStatus status;
    status = AlpsTreeNode::decodeToSelf(encoded);
    //std::printf("BiqTreeNode::decodeToSelf\n");
    encoded.readRep(branchOn_);
    status = dynamic_cast<BiqNodeDesc*>(desc_)->decodeToSelf(encoded);
    assert(status==AlpsReturnStatusOk);
    return status;
}

AlpsKnowledge * BiqTreeNode::decode(AlpsEncoded & encoded) const{
    
    //std::printf("AlpsKnowledge * BiqTreeNode::decode \n");
    //std::printf("AlpsKnowledge * BiqTreeNode::decode desc_ pointer: %p\n", dynamic_cast<BiqNodeDesc*>(desc_));
    BiqTreeNode * nn = new BiqTreeNode(dynamic_cast<BiqNodeDesc*>(desc_)->model());
    nn->decodeToSelf(encoded);
    return nn;
}
