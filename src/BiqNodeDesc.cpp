#include "headers/BiqNodeDesc.h"




BiqNodeDesc::BiqNodeDesc(BiqModel * model) : 
AlpsNodeDesc(),
model_(model),
dQuality_(0.0)
{
    int nVar = model_->getNVar();
    varStatus_.resize(nVar, BiqVarFree);
    bHasBest_ = false;
}

BiqNodeDesc::BiqNodeDesc(BiqModel * model, std::vector<BiqVarStatus> & st) :
AlpsNodeDesc(),
model_(model),
varStatus_(st),
dQuality_(0.0)
{
  bHasBest_ = false;
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
  //std::printf("BiqNodeDesc::encode\n");
  int *ipStatus;
  int sizeStatus;
  
  sizeStatus = varStatus_.size();
  //std::printf("encode stat size: %d\n", sizeStatus);
  if(sizeStatus > 0)
  {
    // give ipstatus the memory needed
    ipStatus = new int[sizeStatus];
    // fill ipstatus
    for(int i = 0; i < sizeStatus; ++i)
    {
      ipStatus[i] = static_cast<int>(varStatus_.at(i));
    }
  }
  else
  {
    ipStatus = 0;
  }
  // easy part
  encoded->writeRep(dQuality_);
  encoded->writeRep(iBranchedOn_);
  encoded->writeRep(iParentsBest_);
  encoded->writeRep(bHasBest_);
  // tricky
  encoded->writeRep(sizeStatus);
  encoded->writeRep(ipStatus,sizeStatus);
  return AlpsReturnStatusOk;
}


/// Unpack fields from the given #AlpsEncoded object.
AlpsReturnStatus BiqNodeDesc::decodeToSelf(AlpsEncoded & encoded) 
{
  int *ipStatus = NULL;
  int sizeStatus;

  encoded.readRep(dQuality_);
  encoded.readRep(iBranchedOn_);
  encoded.readRep(iParentsBest_);
  encoded.readRep(bHasBest_);
  encoded.readRep(sizeStatus);
  //std::printf("decodeToSelf size status: %d\n",sizeStatus);
  encoded.readRep(ipStatus,sizeStatus);
  //std::printf("encoded.readRep(ipStatus,sizeStatus);... \n");
  // build the vector by looping over ipStatus
  varStatus_.resize(sizeStatus);
  for(int i = 0; i < sizeStatus; ++i)
  {
    varStatus_.at(i) = static_cast<BiqVarStatus>(ipStatus[i]);
  }

  return AlpsReturnStatusOk;
}

AlpsNodeDesc * BiqNodeDesc::decode(AlpsEncoded & encoded) const {
  std::cerr << "Not implemented!" << std::endl;
  throw std::exception();
}