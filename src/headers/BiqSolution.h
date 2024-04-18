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
    
      /** Print out the solution. TODO **/
    //virtual void print(std::ostream& os) const;

    using AlpsSolution::encode;
    /// Encode this into the given AlpsEncoded object.
    virtual AlpsReturnStatus encode(AlpsEncoded * encoded) const;
    /// Decode the given AlpsEncoded object into this.
    virtual AlpsReturnStatus decodeToSelf(AlpsEncoded & encoded) {
    std::cerr << "Not implemented!" << std::endl;
    throw std::exception();
    }
    /// Decode the given AlpsEncoded object into a new BiqSolution
    /// and return a pointer to it.
    virtual AlpsKnowledge * decode(AlpsEncoded & encoded) const;

   /** Get the best solution value */
  double getObjValue() const { return value_; }

  virtual double getQuality() const { return getObjValue(); }

    
private:
    BiqSolution(BiqSolution const &);
    BiqSolution & operator=(BiqSolution const &);
};

#endif