#include "headers/BiqNodeDesc.h"




BiqNodeDesc::BiqNodeDesc(BiqModel * model) : 
AlpsNodeDesc(),
model_(model),
dQuality_(0.0)
{
    int nVar = model_->getNVar();
    varStatus_.resize(nVar, BiqVarFree);
}

BiqNodeDesc::BiqNodeDesc(BiqModel * model, std::vector<BiqVarStatus> & st) :
AlpsNodeDesc(),
model_(model),
varStatus_(st),
dQuality_(0.0)
{
}

BiqNodeDesc::~BiqNodeDesc() 
{
}

std::vector<BiqVarStatus> const BiqNodeDesc::getVarStati() const
{
    return varStatus_;
}

void BiqNodeDesc::bound(int & iBound, std::vector<int> solution)
{
    std::cerr << "BiqNodeDesc::bound(int & iBound, std::vector<int> solution)" << std::endl;
}

/*
    MPI methods... 
    TODO: Encode object data using writeRep and Decode with readRep
*/

/// Pack this node description into the given #AlpsEncoded object.
AlpsReturnStatus BiqNodeDesc::encode(AlpsEncoded * encoded) const 
{
  return AlpsReturnStatusOk;
}

/// Unpack fields from the given #AlpsEncoded object.
AlpsReturnStatus BiqNodeDesc::decodeToSelf(AlpsEncoded & encoded) 
{
  return AlpsReturnStatusOk;
}

AlpsNodeDesc * BiqNodeDesc::decode(AlpsEncoded & encoded) const {
  std::cerr << "Not implemented!" << std::endl;
  throw std::exception();
}