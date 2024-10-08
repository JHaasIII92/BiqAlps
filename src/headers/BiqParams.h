#ifndef BiqParams_h_
#define BiqParams_h_

#include "AlpsKnowledge.h"
#include "AlpsParameterBase.h"
#include "stdio.h"

// The three branching stratagies
#define MOST_FRACTIONAL 0
#define LEAST_FRACTIONAL 1
#define CLOSEST_TO_ONE 2

#define BIQ_MAX_CUT    0 
#define BIQ_K_CLUSTER  1
#define BIQ_GENERIC    2
//#############################################################################
//#############################################################################
//** Parameters used in Biq. */
class BiqParams : public AlpsParameterSet {
 public:
  /** Character parameters. All of these variable are used as booleans
      (ture = 1, false = 0). */
  enum boolParams{
    bAddProductConstraints,    
    bAddCuts,
    bScale,
    bSolutionProvided,
    bOneOptSearch,
    bRunWide,
    //
    endOfBoolParams
  };

  /** Integer paramters. */
  enum intParams{
      nCuts,
      nMinCuts,
      MaxNineqAdded,
      nMaxBoundingIter,
      nMaxBFGSIter,
      nMinBoundingIter,
      nMaxAlphaIter,
      nGoemanRuns,
      nPrimalRuns,
      branchingStrategy,
      ProblemType,
      nGreedyRootRuns,
      //
      endOfIntParams
  };

  /** Double parameters. */
  enum dblParams{
    dInitAlpha,
    dScaleAlpha,
    dMinAlpha,
    dScaleTol,
    dMinTol,
    dGapCuts,
    dInitTol,
    dSolutionValue,

    //
    endOfDblParams
  };

  /** String parameters. */
  enum strParams{
    strDummy,

    //
    endOfStrParams
  };

  /** There are no string array parameters. */
  enum strArrayParams{
    strArrayDummy,

    //
    endOfStrArrayParams
  };

 public:
  /**@name Constructors. */
  /*@{*/
  /** The default constructor creates a parameter set with from the template
      argument structure. The keyword list is created and the defaults are
      set. */
  BiqParams() :
    AlpsParameterSet(
    static_cast<int>(endOfBoolParams),
    static_cast<int>(endOfIntParams),
    static_cast<int>(endOfDblParams),
    static_cast<int>(endOfStrParams),
    static_cast<int>(endOfStrArrayParams)
    )
    {
      createKeywordList();
      setDefaultEntries();
    }
  /*@}*/

  /** Method for creating the list of keyword looked for in the parameter
      file. */
  virtual void createKeywordList();
  /** Method for setting the default values for the parameters. */
  virtual void setDefaultEntries();
  /*@}*/


  public:
  //===========================================================================
  /** For user application:
   *   Following code are do NOT need to change.
   *   The reason can not put following functions in base class
   *   <CODE> AlpsParameterSet </CODE> is that <CODE> boolParams </CODE>
   *   and <CODE> endOfBoolParams </CODE> etc., are NOT the same as those
   *   declared in base class.
   */
  //===========================================================================


  /**@name Query methods

     The members of the parameter set can be queried for using the overloaded
     entry() method. Using the example in the class
     documentation the user can get a parameter with the
     "<code>param.entry(USER_par::parameter_name)</code>" expression.
  */
  /*@{*/
  ///
  inline char
    entry(const boolParams key) const { return bpar_[key]; }
  ///
  inline int
    entry(const intParams key) const { return ipar_[key]; }
  ///
  inline double
    entry(const dblParams key) const { return dpar_[key]; }
  ///
  inline const std::string&
    entry(const strParams key) const { return spar_[key]; }
  ///
  inline const std::vector<std::string>&
    entry(const strArrayParams key) const { return sapar_[key]; }
  /*@}*/

  //---------------------------------------------------------------------------
  /// char* is true(1) or false(0), not used
  void setEntry(const boolParams key, const char * val) {
    bpar_[key] = atoi(val) ? true : false; }
  /// char is true(1) or false(0), not used
  void setEntry(const boolParams key, const char val) {
	  bpar_[key] = val ? true : false; }
  /// This method is the one that ever been used.
  void setEntry(const boolParams key, const bool val) {
    bpar_[key] = val; }
  ///
  void setEntry(const intParams key, const char * val) {
    ipar_[key] = atoi(val); }
  ///
  void setEntry(const intParams key, const int val) {
    ipar_[key] = val; }
  ///
  void setEntry(const dblParams key, const char * val) {
    dpar_[key] = atof(val); }
  ///
  void setEntry(const dblParams key, const double val) {
    dpar_[key] = val; }
  ///
  void setEntry(const strParams key, const char * val) {
    spar_[key] = val; }
  ///
  void setEntry(const strArrayParams key, const char *val) {
    sapar_[key].push_back(val); }

  //---------------------------------------------------------------------------

  void print() const {
      printf("Biq Parameters:\n");
      printf("Biq_bAddProductConstraints %d\n", entry(bAddProductConstraints));
      printf("Biq_bAddCuts %d\n", entry(bAddCuts));
      printf("Biq_bScale %d\n", entry(bScale));
      printf("Biq_bSolutionProvided %d\n", entry(bSolutionProvided));
      printf("Biq_nCuts %d\n", entry(nCuts));
      printf("Biq_MaxNineqAdded %d\n", entry(MaxNineqAdded));
      printf("Biq_nMaxBoundingIter %d\n", entry(nMaxBoundingIter));
      printf("Biq_nMinBoundingIter %d\n", entry(nMinBoundingIter));
      printf("Biq_nMinCuts %d\n", entry(nMinCuts));
      printf("Biq_nMaxBFGSIter %d\n", entry(nMaxBFGSIter));
      printf("Biq_nMaxAlphaIter %d\n", entry(nMaxAlphaIter));
      printf("Biq_nGoemanRuns %d\n", entry(nGoemanRuns));
      printf("Biq_nPrimalRuns %d\n", entry(nPrimalRuns));  
      printf("Biq_branchingStrategy %d\n", entry(branchingStrategy)); 
      printf("Biq_dInitAlpha %.2e\n", entry(dInitAlpha));
      printf("Biq_dScaleAlpha %.2f\n", entry(dScaleAlpha));
      printf("Biq_dMinAlpha %.2e\n", entry(dMinAlpha));
      printf("Biq_dInitTol %.2e\n", entry(dInitTol));
      printf("Biq_dScaleTol %.2f\n", entry(dScaleTol));
      printf("Biq_dMinTol %.2e\n", entry(dMinTol));
      printf("Biq_dGapCuts %.2e\n", entry(dGapCuts));
      printf("Biq_dSolutionValue %.2e\n", entry(dSolutionValue));  
  }

  /**@name Packing/unpacking methods */
  /*@{*/
  /** Pack the parameter set into the buffer (AlpsEncoded is used
      as buffer Here). */
  void pack(AlpsEncoded& buf) {
    buf.writeRep(bpar_, endOfBoolParams)
        .writeRep(ipar_, endOfIntParams)
        .writeRep(dpar_, endOfDblParams);
    for (int i = 0; i < endOfStrParams; ++i)
      buf.writeRep(spar_[i]);
    for (int i = 0; i < endOfStrArrayParams; ++i) {
      buf.writeRep(sapar_[i].size());
      for (size_t j = 0; j < sapar_[i].size(); ++j)
        buf.writeRep(sapar_[i][j]);
    }
  }
  /** Unpack the parameter set from the buffer. */
  void unpack(AlpsEncoded& buf) {
    int dummy;
    // No need to allocate the arrays, they are of fixed length
    dummy = static_cast<int>(endOfBoolParams);
    buf.readRep(bpar_, dummy, false);
    dummy = static_cast<int>(endOfIntParams);
    buf.readRep(ipar_, dummy, false);
    dummy = static_cast<int>(endOfDblParams);
    buf.readRep(dpar_, dummy, false);
    for (int i = 0; i < endOfStrParams; ++i)
      buf.readRep(spar_[i]);
    for (int i = 0; i < endOfStrArrayParams; ++i) {
      size_t str_size;
      buf.readRep(str_size);
      sapar_[i].reserve(str_size);
      for (size_t j = 0; j < str_size; ++j){
        //	sapar_[i].unchecked_push_back(std::string());
        sapar_[i].push_back(std::string());
        buf.readRep(sapar_[i].back());
      }
    }
  }
  /*@}*/

};

#endif
