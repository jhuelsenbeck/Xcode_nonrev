#ifndef ParameterBaseFrequencies_H
#define ParameterBaseFrequencies_H

#include <vector>
#include "Parameter.h"



class ParameterBaseFrequencies : public Parameter {

    public:
                                        ParameterBaseFrequencies(RandomVariable* rp, Model* mp, std::string nm, double tn);
                                        ParameterBaseFrequencies(ParameterBaseFrequencies& b);
        Parameter&                      operator=(Parameter& b);
        ParameterBaseFrequencies&       operator=(ParameterBaseFrequencies& b);
        double&                         operator[](int idx);
        void                            print(void);
        double                          update(void);
        double                          lnPriorProb(void);
        std::string                     getParmString(int n);
        std::string                     getParmHeader(int n);
        std::vector<double>             getBaseFrequencies(void) { return f; }
        std::vector<double>&            getBaseFrequencyAlpha(void) { return a; }
        double                          getTuning(void) { return alpha0; }
        void                            setTuning(double x) { alpha0 = x; }

    protected:
        void                            clone(ParameterBaseFrequencies& b);
        double                          alpha0;
        std::vector<double>             a;
        std::vector<double>             f;
};

#endif
