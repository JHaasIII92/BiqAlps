#ifndef BiqSolution_h_
#define BiqSolution_h_

#include "AlpsSolution.h"
#include "AlpsKnowledge.h"
#include <vector>
class BiqModel;
class BiqSolution : public AlpsSolution
{
private:

    /* data */
    BiqModel * model_;
    std::vector<int> solution_;
    int value_;

public:

    /* No pure virtuals */
    BiqSolution(BiqModel * model);
    BiqSolution(BiqModel * model, std::vector<int> sol, int val);
    ~BiqSolution();
    virtual AlpsKnowledge * decode(AlpsEncoded & encode) const;

   /** Get the best solution value */
  double getObjValue() const { return value_; }

  virtual double getQuality() const { return getObjValue(); }

    
private:
    BiqSolution(BiqSolution const &);
    BiqSolution & operator=(BiqSolution const &);
};

#endif