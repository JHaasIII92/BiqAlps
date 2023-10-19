#ifndef BiqSolution_h_
#define BiqSolution_h_

#include "AlpsSolution.h"
#include "AlpsKnowledge.h"

class BiqModel;
class BiqSolution : public AlpsSolution
{
private:

    /* data */
    BiqModel * model_;
    int * solution_;
    int value_;

public:

    /* No pure virtuals */
    BiqSolution(BiqModel * model);
    BiqSolution(BiqModel * model, int *& sol, int val);
    ~BiqSolution();
    virtual AlpsKnowledge * decode(AlpsEncoded & encode) const;
    
private:
    BiqSolution(BiqSolution const &);
    BiqSolution & operator=(BiqSolution const &);
};

#endif