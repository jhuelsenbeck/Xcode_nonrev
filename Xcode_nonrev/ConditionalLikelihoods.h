#ifndef ConditionalLikelihoods_H
#define ConditionalLikelihoods_H

class Alignment;


class ConditionalLikelihoods {

    public:
                                    ConditionalLikelihoods(void) = delete;
                                    ConditionalLikelihoods(Alignment* a, int nc);
                                   ~ConditionalLikelihoods(void);
        ConditionalLikelihoods&     operator=(ConditionalLikelihoods& c);
        double*                     getCondLike(void) { return cls; }
        double*                     getLnScaler(void) { return lnScaler; }
        double*                     getLnScalerDp(void) { return lnScalerDp; }
        void                        initializeTipConditonalLikelihoods(Alignment* a, std::string tName);
        void                        print(void);
        
    protected:
        double*                     cls;
        double*                     lnScaler;
        double*                     lnScalerDp;
        int                         numGammaCats;
        int                         numSites;
        int                         clsSizeInBytes;
};

#endif
