#include "AlpsParameterBase.h"
#include "headers/BiqParams.h"

using std::make_pair;


void
BiqParams::createKeywordList() {
    // Create the list of keywords for parameter file reading
    //-------------------------------------------------------------------------
    // CharPar

    //-------------------------------------------------------------------------
    // BoolArrayPar

    //-------------------------------------------------------------------------
    // IntPar
    keys_.push_back(make_pair(std::string("Biq_boundStatusInterval"),
			      AlpsParameter(AlpsIntPar, boundStatusInterval)));

    //-------------------------------------------------------------------------
    // DoublePar

    //-------------------------------------------------------------------------
    // StringPar

}

//#############################################################################

void
BiqParams::setDefaultEntries() {
    //-------------------------------------------------------------------------
    // CharPar

    //-------------------------------------------------------------------------
    // IntPar
    setEntry(boundStatusInterval, 2);
    //-------------------------------------------------------------------------
    // DoublePar

    //-------------------------------------------------------------------------
    // StringPar

}