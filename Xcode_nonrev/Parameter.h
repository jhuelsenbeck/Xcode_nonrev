#ifndef Parameter_H
#define Parameter_H

#include <string>
class Model;
class RandomVariable;

class Parameter {

    public:
                                Parameter(RandomVariable* rp, Model* mp, std::string nm);
        virtual                ~Parameter(void);
        virtual Parameter&      operator=(Parameter& b)=0;
        virtual void            print(void)=0;
        virtual double          update(void)=0;
        virtual double          lnPriorProb(void)=0;
        virtual std::string     getParmString(int n)=0;
        std::string             getName(void) { return name; }
        virtual std::string     getParmHeader(int n)=0;

    protected:
        RandomVariable*         rv;
        Model*                  modelPtr;
        std::string             name;
};

#endif
