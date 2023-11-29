#include "headers/BiqNodeDesc.h"




BiqNodeDesc::BiqNodeDesc(BiqModel * model) : 
AlpsNodeDesc(),
model_(model)
{
    std::cout << "BiqNodeDesc::BiqNodeDesc(BiqModel * model)" << std::endl;
    int nVar = model_->getNVar();
    varStatus_.resize(nVar, BiqVarFree);
}

BiqNodeDesc::BiqNodeDesc(BiqModel * model, std::vector<BiqVarStatus> & st) :
AlpsNodeDesc(),
model_(model),
varStatus_(st)
{
    std::cout << "BiqNodeDesc::BiqNodeDesc(BiqModel * model, std::vector<BiqVarStatus> & st) " << std::endl;
}

BiqNodeDesc::~BiqNodeDesc() 
{
    std::cout << " BiqNodeDesc::~BiqNodeDesc ... Freeing BiqNodeDesc ..." << std::endl;
}

std::vector<BiqVarStatus> const BiqNodeDesc::getVarStati() const
{
    return varStatus_;
}

AlpsKnowledge * BiqNodeDesc::decode(AlpsEncoded & encode) const
{
    std::cerr << "BiqNodeDesc::decode ... Not implemented!" << std::endl;
    throw std::exception();
}


void BiqNodeDesc::bound(int & iBound, std::vector<int> solution)
{
    std::cerr << "BiqNodeDesc::bound(int & iBound, std::vector<int> solution)" << std::endl;
}


                            