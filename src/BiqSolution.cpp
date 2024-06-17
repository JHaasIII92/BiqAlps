#include "headers/BiqSolution.h"
#include<iostream>

BiqSolution::BiqSolution(BiqModel * model, std::vector<int> sol, int val) : 
model_(model),
solution_(sol),
value_(val)
{ 

}

BiqSolution::BiqSolution(BiqModel * model) :
model_(model),
solution_(0),
value_(0)
{
 
}

BiqSolution::~BiqSolution()
{

}
AlpsReturnStatus decodeToSelf(AlpsEncoded & encoded)
{
    std::printf("AlpsKnowledge * BiqSolution::decode\n");
    std::cerr << "Not implemented!" << std::endl;
    throw std::exception();
}


void BiqSolution::print(std::ostream& os) const 
{

  os << "BiqAlps Solution Report: " << std::endl;
  os << "Value: " << value_<< std::endl;
  int i = 1;
  for(auto & it: solution_)
  {
    os << "x[" << i << "] = " << it << " ";
    if(i%5 == 0) os << std::endl;

    ++i;
  }

  os << std::endl;

}

AlpsReturnStatus BiqSolution::encode(AlpsEncoded * encoded) const  
{
  
  int* ipSolution;
  int sizeSol;

  sizeSol = solution_.size();
  if(sizeSol > 0)
  {
    ipSolution = new int[sizeSol];  
    for(int i = 0; i < sizeSol; ++i)
    {
      ipSolution[i] = solution_.at(i);
    }
  }

  encoded->writeRep(value_);

  encoded->writeRep(sizeSol);
  encoded->writeRep(ipSolution,sizeSol);

  return AlpsReturnStatusOk;
}

// Note: write and read sequence MUST same!
AlpsKnowledge * BiqSolution::decode(AlpsEncoded& encoded) const {
  
  int* ipSolution = NULL;
  int sizeSol = 0;

  std::vector<int> sol;
  int val;

  encoded.readRep(val);

  encoded.readRep(sizeSol);
  encoded.readRep(ipSolution, sizeSol);

  
  sol.resize(sizeSol);
  for(int i = 0; i < sizeSol; ++i)
  {
    sol.at(i) = ipSolution[i];
  }

  return new BiqSolution(model_, sol, val);
}