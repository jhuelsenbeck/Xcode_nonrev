#include <iomanip>
#include <iostream>
#include "ParameterExchangabilityRates.h"
#include "RandomVariable.h"



ParameterExchangabilityRates::ParameterExchangabilityRates(RandomVariable* rp, Model* mp, std::string nm, double tn, bool isTimeReversible) : Parameter(rp, mp, nm) {

    int n = 6;
    if (isTimeReversible == false)
        n = 12;
    f.resize( n );
    a.resize( n );
    alpha0 = tn;
    for (int i=0; i<a.size(); i++)
        a[i] = 1.0;
    rv->dirichletRv(a, f);
    for (int i=0; i<f.size(); i++)
        f[i] = a[i] / (double)n; // expected value, for now
}

ParameterExchangabilityRates::ParameterExchangabilityRates(ParameterExchangabilityRates& b) : Parameter(b.rv, b.modelPtr, b.name) {

    f.resize( b.f.size() );
    a.resize( b.a.size() );
    clone(b);
}

Parameter& ParameterExchangabilityRates::operator=(Parameter& b) {

    if (this != &b)
        {
        ParameterExchangabilityRates* dc = dynamic_cast<ParameterExchangabilityRates*>(&b);
        clone(*dc);
        }
    return *this;
}

ParameterExchangabilityRates& ParameterExchangabilityRates::operator=(ParameterExchangabilityRates& b) {

    if (this != &b)
        {
        clone(b);
        }
    return *this;
}

double& ParameterExchangabilityRates::operator[](int idx) {

    if (idx >= f.size())
        {
        std::cout << "Exchangability index out of bounds";
        exit(0);
        }
    return f[idx];
}

void ParameterExchangabilityRates::clone(ParameterExchangabilityRates& b) {

    alpha0 = b.alpha0;
    for (int i=0; i<b.a.size(); i++)
        {
        a[i] = 1.0;
        f[i] = b.f[i];
        }
}

void ParameterExchangabilityRates::print(void) {

    std::cout << "Exchangability Rates = ";
    for (int i=0; i<f.size(); i++)
        std::cout << std::fixed << std::setprecision(4) << f[i] << " ";
    std::cout << '\n';
}

double ParameterExchangabilityRates::update(void) {
    
    int n = (int)a.size();
    
    std::vector<double> aForward(n);
    std::vector<double> aReverse(n);
    std::vector<double> oldFreqs(n);
    for (int i=0; i<n; i++)
        {
        oldFreqs[i] = f[i];
        aForward[i] = f[i] * alpha0;
        if (aForward[i] < 1E-100)
            return -(1E-310);
        }
    rv->dirichletRv(aForward, f);
    for (int i=0; i<n; i++)
        {
        aReverse[i] = f[i] * alpha0;
        if (aReverse[i] < 1E-100)
            return -(1E-310);
        }
    return rv->lnDirichletPdf(aReverse, oldFreqs) - rv->lnDirichletPdf(aForward, f);
}

double ParameterExchangabilityRates::lnPriorProb(void) {

    return rv->lnDirichletPdf(a, f);
}

std::string ParameterExchangabilityRates::getParmString(int n) {

    std::string tempStr = "";
    for (int i=0; i<f.size(); i++)
        {
        char temp[20];
        sprintf(temp, "%1.4lf\t", f[i]);
        tempStr += temp;
        }
    return tempStr;
}

std::string ParameterExchangabilityRates::getParmHeader(int n) {

    char temp[500];
    sprintf (temp, "R[AC]\tR[AG]\tR[AT]\tR[CG]\tR[CT]\tR[GT]\t");
    std::string tempString = temp;
    return tempString;
}