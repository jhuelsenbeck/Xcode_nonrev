#include <iomanip>
#include <iostream>
#include "Alignment.h"
#include "ConditionalLikelihoods.h"
#include "Msg.h"


ConditionalLikelihoods::ConditionalLikelihoods(Alignment* a, int nc) {

    // initialize some instance variables
    numGammaCats   = nc;
    numSites       = a->getNumChar();
    clsSizeInBytes = numSites * 4 * numGammaCats * sizeof(double);
    
    // dynamically allocate the information for the conditional likelihoods
    int clsSize = numSites * 4 * numGammaCats;
    cls = new double[clsSize];
    for (int i=0; i<clsSize; i++)
        cls[i] = 0.0;
    
    // allocate the log scaler vector
    lnScaler = new double[numSites];
    lnScalerDp = new double[numSites];
    for (int i=0; i<numSites; i++)
        {
        lnScaler[i] = 0.0;
        lnScalerDp[i] = 0.0;
        }
}

ConditionalLikelihoods::~ConditionalLikelihoods(void) {

    delete [] cls;
    delete [] lnScaler;
    delete [] lnScalerDp;
}

ConditionalLikelihoods& ConditionalLikelihoods::operator=(ConditionalLikelihoods& c) {

    if (this != &c)
        {
        numGammaCats = c.numGammaCats;
        numSites = c.numSites;
        clsSizeInBytes = c.clsSizeInBytes;
        memcpy(cls, c.cls, clsSizeInBytes);
        memcpy(lnScaler, c.lnScaler, numSites*sizeof(double));
        memcpy(lnScalerDp, c.lnScalerDp, numSites*sizeof(double));
        }
    return *this;
}

void ConditionalLikelihoods::initializeTipConditonalLikelihoods(Alignment* a, std::string tName) {

    int idx = a->getTaxonIndex(tName);
    if (idx == -1)
        Msg::error("Couldn't find taxon " + tName + " in alignment");
    
    double* p = &cls[0];
    for (int c=0; c<numSites; c++)
        {
        int nucCode = a->getNucleotide(idx, c);
        int possibleNucs[4];
        a->getPossibleNucs (nucCode, possibleNucs);
        for (int k=0; k<numGammaCats; k++)
            {
            for (int s=0; s<4; s++)
                p[s] = (double)possibleNucs[s];
            p += 4;
            }
        }
}

void ConditionalLikelihoods::print(void) {

    double* p = &cls[0];
    for (int c=0; c<numSites; c++)
        {
        std::cout << std::fixed << std::setprecision(0);
        for (int s=0; s<4; s++)
            std::cout << p[s];
        std::cout << " ";
        p += 4 * numGammaCats;
        }
    std::cout << std::endl;
}
