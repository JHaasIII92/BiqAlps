#include "headers/BiqSolution.h"
#include<iostream>

BiqSolution::BiqSolution(BiqModel * model, int *& sol, int val) : 
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
    std::cout << "Freeing BiqSolution ..." << std::endl;
    if(solution_)
    {
        delete[] solution_;
        solution_ = 0;
    }
    std::cout << "Freeing BiqSolution ..." << std::endl;
}
AlpsKnowledge * BiqSolution::decode(AlpsEncoded & encode) const
{
    std::cerr << "Not implemented!" << std::endl;
    throw std::exception();
}