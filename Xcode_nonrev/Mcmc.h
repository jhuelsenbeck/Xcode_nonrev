#ifndef Mcmc_H
#define Mcmc_H

#include <fstream>
#include <string>
class Alignment;
class Model;
class RandomVariable;
class Settings;


class Mcmc {

    public:
                            Mcmc(void) = delete;
                            Mcmc(Settings* s, Model* m, RandomVariable* r, Alignment* a);
                           ~Mcmc(void);
        void                run(void);
        void                runPowerPosterior(void);

    protected:
        double              safeExponentiation(double lnX);
        void                sampleChain(int n, int len, double lnL, double lnP);
        void                samplePower(int sn, int n, int len, double lnL, double power);
        Model*              myModel;
        RandomVariable*     rv;
        Settings*           mySettings;
        Alignment*          myAlignment;
        std::ofstream       parmOut;
        std::ofstream       treeOut;
        std::ofstream       powerOut;
        std::string         parmFile;
        std::string         treeFile;
        std::string         powerFile;
        int                 burnin;
};

#endif
