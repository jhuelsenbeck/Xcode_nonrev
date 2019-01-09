#ifndef ParameterExchangabilityRates_H
#define ParameterExchangabilityRates_H

#include <vector>
#include "Parameter.h"



class ParameterExchangabilityRates : public Parameter {

    public:
                                        ParameterExchangabilityRates(RandomVariable* rp, Model* mp, std::string nm, double tn, bool isTimeReversible);
                                        ParameterExchangabilityRates(ParameterExchangabilityRates& b);
        Parameter&                      operator=(Parameter& b);
        ParameterExchangabilityRates&   operator=(ParameterExchangabilityRates& b);
        double&                         operator[](int idx);
        void                            print(void);
        double                          update(void);
        double                          lnPriorProb(void);
        std::string                     getParmString(int n);
        std::string                     getParmHeader(int n);
        std::vector<double>&            getExchangabilityRates(void) { return f; }
        std::vector<double>&            getExchangabilityAlpha(void) { return a; }
        double                          getTuning(void) { return alpha0; }
        void                            setTuning(double x) { alpha0 = x; }
        int                             size(void) { return (int)f.size(); }

    protected:
        void                            clone(ParameterExchangabilityRates& b);
        double                          logit(double x);
        double                          ilogit(double x);
        double                          alpha0;
        std::vector<double>             a;
        std::vector<double>             f;
};

#endif
