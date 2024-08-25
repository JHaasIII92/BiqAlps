#include "AlpsParameterBase.h"
#include "headers/BiqParams.h"

using std::make_pair;


void
BiqParams::createKeywordList() {
    // Create the list of keywords for parameter file reading
    //-------------------------------------------------------------------------
    // CharPar

    //-------------------------------------------------------------------------
    // BoolPar
    keys_.push_back(make_pair(std::string("Biq_bAddProductConstraints"),
			      AlpsParameter(AlpsBoolPar, bAddProductConstraints)));
    keys_.push_back(make_pair(std::string("Biq_bAddCuts"),
			      AlpsParameter(AlpsBoolPar, bAddCuts)));
    keys_.push_back(make_pair(std::string("Biq_bScale"),
			      AlpsParameter(AlpsBoolPar, bScale)));
    keys_.push_back(make_pair(std::string("Biq_bSolutionProvided"),
                  AlpsParameter(AlpsBoolPar, bSolutionProvided)));
    keys_.push_back(make_pair(std::string("Biq_bOneOptSearch"),
                  AlpsParameter(AlpsBoolPar, bOneOptSearch)));
    keys_.push_back(make_pair(std::string("Biq_bRunWide"),
                  AlpsParameter(AlpsBoolPar, bRunWide)));                 
    //-------------------------------------------------------------------------
    // IntPar
    keys_.push_back(make_pair(std::string("Biq_nCuts"),
			      AlpsParameter(AlpsIntPar, nCuts)));
    keys_.push_back(make_pair(std::string("Biq_MaxNineqAdded"),
			      AlpsParameter(AlpsIntPar, MaxNineqAdded)));
    keys_.push_back(make_pair(std::string("Biq_nMinCuts"),
			      AlpsParameter(AlpsIntPar, nMinCuts)));
    keys_.push_back(make_pair(std::string("Biq_nMaxBoundingIter"),
			      AlpsParameter(AlpsIntPar, nMaxBoundingIter)));
    keys_.push_back(make_pair(std::string("Biq_nMinBoundingIter"),
			      AlpsParameter(AlpsIntPar, nMinBoundingIter)));
    keys_.push_back(make_pair(std::string("Biq_nMaxBFGSIter"),
			      AlpsParameter(AlpsIntPar, nMaxBFGSIter)));
    keys_.push_back(make_pair(std::string("Biq_nMaxAlphaIter"),
			      AlpsParameter(AlpsIntPar, nMaxAlphaIter)));
    keys_.push_back(make_pair(std::string("Biq_nGoemanRuns"),
			      AlpsParameter(AlpsIntPar, nGoemanRuns)));
    keys_.push_back(make_pair(std::string("Biq_nPrimalRuns"),
			      AlpsParameter(AlpsIntPar, nPrimalRuns))); 
    keys_.push_back(make_pair(std::string("Biq_branchingStrategy"),
			      AlpsParameter(AlpsIntPar, branchingStrategy)));
    keys_.push_back(make_pair(std::string("Biq_ProblemType"),
                  AlpsParameter(AlpsIntPar, ProblemType)));
    keys_.push_back(make_pair(std::string("Biq_nGreedyRootRuns"),
                  AlpsParameter(AlpsIntPar, nGreedyRootRuns)));
    //-------------------------------------------------------------------------
    // DoublePar
    keys_.push_back(make_pair(std::string("Biq_dInitAlpha"),
			      AlpsParameter(AlpsDoublePar, dInitAlpha)));
    keys_.push_back(make_pair(std::string("Biq_dScaleAlpha"),
			      AlpsParameter(AlpsDoublePar, dScaleAlpha)));
    keys_.push_back(make_pair(std::string("Biq_dMinAlpha"),
			      AlpsParameter(AlpsDoublePar, dMinAlpha)));
    keys_.push_back(make_pair(std::string("Biq_dInitTol"),
			      AlpsParameter(AlpsDoublePar, dInitTol)));
    keys_.push_back(make_pair(std::string("Biq_dScaleTol"),
			      AlpsParameter(AlpsDoublePar, dScaleTol)));
    keys_.push_back(make_pair(std::string("Biq_dMinTol"),
			      AlpsParameter(AlpsDoublePar, dMinTol)));
    keys_.push_back(make_pair(std::string("Biq_dGapCuts"),
			      AlpsParameter(AlpsDoublePar, dGapCuts)));
    keys_.push_back(make_pair(std::string("Biq_dSolutionValue"),
                    AlpsParameter(AlpsDoublePar, dSolutionValue)));
    //-------------------------------------------------------------------------
    // StringPar

}

//#############################################################################

void
BiqParams::setDefaultEntries() {
    //std::printf("BiqParams::setDefaultEntries()\n");
    //-------------------------------------------------------------------------
    // CharPar

    //-------------------------------------------------------------------------
    // BoolPar
    setEntry(bAddProductConstraints, true);
    setEntry(bAddCuts, true);
    setEntry(bScale, true);
    setEntry(bSolutionProvided, false);
    setEntry(bOneOptSearch, false);
    setEntry(bRunWide, false);

    //-------------------------------------------------------------------------
    // IntPar
    setEntry(nCuts, 500);
    setEntry(MaxNineqAdded, 10000);
    setEntry(nMaxBoundingIter, 100);
    setEntry(nMinBoundingIter, 12);
    setEntry(nMinCuts, 50);
    setEntry(nMaxBFGSIter, 2000);
    setEntry(nMaxAlphaIter, 50);
    setEntry(nGoemanRuns, 10);
    setEntry(nPrimalRuns, 1);
    setEntry(branchingStrategy, MOST_FRACTIONAL);
    setEntry(ProblemType, BIQ_MAX_CUT);
    setEntry(nGreedyRootRuns, 1);
    //-------------------------------------------------------------------------
    // DoublePar
    setEntry(dInitAlpha, 1e-1);
    setEntry(dScaleAlpha, 0.5);
    setEntry(dMinAlpha, 5e-5);
    setEntry(dInitTol, 1e-1);
    setEntry(dScaleTol, 0.95);
    setEntry(dMinTol, 1e-2);
    setEntry(dGapCuts, -5e-2);
    setEntry(dSolutionValue, 0.0);
    //-------------------------------------------------------------------------
    // StringPar

}