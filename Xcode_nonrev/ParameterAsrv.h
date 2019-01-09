#ifndef ParameterAsrv_H
#define ParameterAsrv_H

#include <string>
#include <vector>
#include "Parameter.h"



class ParameterAsrv : public Parameter {

    public:
                                ParameterAsrv(RandomVariable* rp, Model* mp, std::string nm, double tn, int n, double alp);
                                ParameterAsrv(ParameterAsrv& b);
        Parameter&              operator=(Parameter& b);
        ParameterAsrv&          operator=(ParameterAsrv& b);
        double&                 operator[](int idx);
        std::string             getParmHeader(int n);
        std::string             getParmString(int n);
        double                  lnPriorProb(void);
        void                    print(void);
        double                  update(void);
        double&                 getAsrv(void) { return alpha; }
        std::vector<double>&    getRates(void) { return rates; }
        double                  getTuning(void) { return tuning; }
        void                    setTuning(double x) { tuning = x; }

    protected:
        void                    clone(ParameterAsrv& b);
        double                  alpha;
        double                  alphaPrior;
        int                     numCategories;
        std::vector<double>     rates;
        double                  tuning;
};

#endif
