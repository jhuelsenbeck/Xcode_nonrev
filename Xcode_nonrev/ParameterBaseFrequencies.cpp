#include <cmath>
#include <iomanip>
#include <iostream>
#include "ParameterBaseFrequencies.h"
#include "RandomVariable.h"



ParameterBaseFrequencies::ParameterBaseFrequencies(RandomVariable* rp, Model* mp, std::string nm, double tn) : Parameter(rp, mp, nm) {

    f.resize( 4 );
    a.resize( 4 );
    alpha0 = tn;
    for (int i=0; i<4; i++)
        a[i] = 1.0;
    rv->dirichletRv(a, f);
    for (int i=0; i<4; i++)
        f[i] = a[i] / 4.0; // set to the expected value, for now
}

ParameterBaseFrequencies::ParameterBaseFrequencies(ParameterBaseFrequencies& b) : Parameter(b.rv, b.modelPtr, b.name) {

    f.resize( 4 );
    a.resize( 4 );
    clone(b);
}

Parameter& ParameterBaseFrequencies::operator=(Parameter& b) {

    if (this != &b)
        {
        ParameterBaseFrequencies* dc = dynamic_cast<ParameterBaseFrequencies*>(&b);
        clone(*dc);
        }
    return *this;
}

ParameterBaseFrequencies& ParameterBaseFrequencies::operator=(ParameterBaseFrequencies& b) {

    if (this != &b)
        {
        clone(b);
        }
    return *this;
}

double& ParameterBaseFrequencies::operator[](int idx) {

    if (idx >= 4)
        {
        std::cout << "Base frequency index out of bounds";
        exit(0);
        }
    return f[idx];
}

void ParameterBaseFrequencies::clone(ParameterBaseFrequencies& b) {

    alpha0 = b.alpha0;
    for (int i=0; i<4; i++)
        {
        a[i] = 1.0;
        f[i] = b.f[i];
        }
}

void ParameterBaseFrequencies::print(void) {

    std::cout << "Base Frequencies = ";
    for (int i=0; i<4; i++)
        std::cout << std::fixed << std::setprecision(4) << f[i] << " ";
    std::cout << '\n';
}

double ParameterBaseFrequencies::update(void) {

#   if 1
    
    // pick the element to update
    int k = (int)(rv->uniformRv() * 4);
    
    //std::cout << "(" << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << ") -> ";
    // change value
    std::vector<double> oldFreqs(2);
    std::vector<double> newFreqs(2);
    std::vector<double> aForward(2);
    std::vector<double> aReverse(2);
    oldFreqs[0] = f[k];
    oldFreqs[1] = 1.0 - f[k];
    for (int i=0; i<2; i++)
        {
        aForward[i] = oldFreqs[i] * alpha0;
        if (aForward[i] < 1E-100)
            return -(1E-310);
        }
    rv->dirichletRv(aForward, newFreqs);
    for (int i=0; i<2; i++)
        {
        aReverse[i] = newFreqs[i] * alpha0;
        if (aReverse[i] < 1E-100)
            return -(1E-310);
        }
    f[k] = newFreqs[0];
    double factor = newFreqs[1] / oldFreqs[1];
    for (int i=0; i<4; i++)
        {
        if (i != k)
            f[i] *= factor;
        }

    // renormalize, just in case
    double sum = 0.0;
    for (int i=0; i<4; i++)
        sum += f[i];
    for (int i=0; i<4; i++)
        f[i] /= sum;

    //std::cout << "(" << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << ")" << std::endl;

    // return the proposal ratio
    double lnProposalRatio  = rv->lnDirichletPdf(aReverse, oldFreqs) - rv->lnDirichletPdf(aForward, newFreqs); // Hastings Ratio
    lnProposalRatio += (4 - 2) * log(factor); // Jacobian
    return lnProposalRatio;

#   else

    std::vector<double> aForward(4);
    std::vector<double> aReverse(4);
    std::vector<double> oldFreqs(4);
    for (int i=0; i<4; i++)
        {
        oldFreqs[i] = f[i];
        aForward[i] = f[i] * alpha0;
        if (aForward[i] < 1E-100)
            return -(1E-310);
        }
    rv->dirichletRv(aForward, f);
    for (int i=0; i<4; i++)
        {
        aReverse[i] = f[i] * alpha0;
        if (aReverse[i] < 1E-100)
            return -(1E-310);
        }
    return rv->lnDirichletPdf(aReverse, oldFreqs) - rv->lnDirichletPdf(aForward, f);
    
#   endif
}

double ParameterBaseFrequencies::lnPriorProb(void) {

    return rv->lnDirichletPdf(a, f);
}

std::string ParameterBaseFrequencies::getParmString(int n) {

    std::string tempStr = "";
    for (int i=0; i<4; i++)
        {
        char temp[20];
        sprintf(temp, "%1.4lf\t", f[i]);
        tempStr += temp;
        }
    return tempStr;
}

std::string ParameterBaseFrequencies::getParmHeader(int n) {

    char temp[500];
    sprintf (temp, "Pi[A]\tPi[C]\tPi[G]\tPi[T]\t");
    std::string tempString = temp;
    return tempString;
}
