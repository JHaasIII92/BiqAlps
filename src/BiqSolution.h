#ifndef BiqSolution_h_
#define BiqSolution_h_

#include "AlpsSolution.h"
#include "AlpsKnowledge.h"

class BiqModel;
class BiqSolution : public AlpsSolution
{
    /* data */
public:

    /* No pure virtuals */
    BiqSolution(BiqModel * model) {};
    virtual AlpsKnowledge * decode(AlpsEncoded & encode) const;
};

#endif