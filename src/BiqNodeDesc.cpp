
#include <iostream>
#include "headers/BiqNodeDesc.h"
#include "headers/BiqModel.h"



BiqNodeDesc::BiqNodeDesc(BiqModel * model) : 
AlpsNodeDesc(),
model_(model)
{
    // make a all free BiqVarStatus
    int nVar = model_->getNVar();
    std::cout << " BiqNodeDesc::BiqNodeDesc ... The pointer to varStatus " << varStatus_ << std::endl; 
    std::cout << " BiqNodeDesc::BiqNodeDesc ...  nVar from model = " << nVar << std::endl;
    varStatus_ = new BiqVarStatus[nVar];
    std::fill_n(varStatus_, nVar, BiqVarFree);
    std::cout << "BiqNodeDesc::BiqNodeDesc ... just filled in varStatus_" << std::endl;
}

BiqNodeDesc::BiqNodeDesc(BiqModel * model, BiqVarStatus *& st) :
AlpsNodeDesc(),
model_(model),
varStatus_(st)
{
    std::cout << "BiqNodeDesc::BiqNodeDesc(BiqModel * model, BiqVarStatus *& st)" << std::endl;
    std::cout << "BiqNodeDesc::BiqNodeDesc  ...  varStatus_[0] = " << varStatus_[0] << std::endl;
}

BiqNodeDesc::~BiqNodeDesc() 
{
    std::cout << " BiqNodeDesc::~BiqNodeDesc ... Freeing BiqNodeDesc ..." << std::endl;
    if(varStatus_)
    {
        delete[] varStatus_;
        varStatus_ = NULL;
    }
    std::cout << " BiqNodeDesc::~BiqNodeDesc ... Freed BiqNodeDesc ..." << std::endl;
}

BiqVarStatus const * BiqNodeDesc::getVarStati() const
{
    return varStatus_;
}

AlpsKnowledge * BiqNodeDesc::decode(AlpsEncoded & encode) const
{
    std::cerr << "BiqNodeDesc::decode ... Not implemented!" << std::endl;
    throw std::exception();
}


void BiqNodeDesc::bound(int & iBound, double * pdSol)
{
    // filling solution vector with 0.5 and setting bound to 0;
    int nVar = model_->getNVar();
    std::fill_n(pdSol, nVar, 0.5);
    // make a random bound between 1 and 10
    iBound = rand() % 10 + 1;
}


                            