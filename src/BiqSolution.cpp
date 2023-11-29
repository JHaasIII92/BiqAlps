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
AlpsKnowledge * BiqSolution::decode(AlpsEncoded & encode) const
{
    std::printf("AlpsKnowledge * BiqSolution::decode\n");
    std::cerr << "Not implemented!" << std::endl;
    throw std::exception();
}