#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include "Model.h"
#include "Node.h"
#include "ParameterAsrv.h"
#include "TransitionProbabilities.h"
//#include "ParameterTree.h"



TransitionProbabilities::TransitionProbabilities(int s, int k, Model* m) {

    numStates = 4;
    numGammaCats = k;
    mySpace = s;
    myModel = m;
    tiSizeInBytes = numStates * numStates * numGammaCats * sizeof(double);
    
    // dynamically allocate a vector of matrices
#   if 0
    p = new double**[numGammaCats];
    for (int c=0; c<numGammaCats; c++)
        {
        p[c] = new double*[numStates];
        p[c][0] = new double[numStates * numStates];
        for (int i=1; i<numStates; i++)
            p[c][i] = p[c][i-1] + numStates;
        for (int i=0; i<numStates; i++)
            for (int j=0; j<numStates; j++)
                p[c][i][j] = 0.0;
        }
#   else
    p = new double**[numGammaCats];
    p[0] = new double*[numGammaCats * numStates];
    for (int c=1; c<numGammaCats; c++)
        p[c] = p[c-1] + numStates;
    p[0][0] = new double[numStates * numStates * numGammaCats];
    double* ptr = p[0][0];
    for (int c=0; c<numGammaCats; c++)
        {
        for (int i=0; i<numStates; i++)
            {
            p[c][i] = ptr;
            ptr += numStates;
            }
        }
    for (int c=0; c<numGammaCats; c++)
        {
        for (int i=0; i<numStates; i++)
            for (int j=0; j<numStates; j++)
                p[c][i][j] = 0.0;
        }
#   endif
}

TransitionProbabilities::~TransitionProbabilities(void) {

    delete [] p[0][0];
    delete [] p[0];
    delete [] p;
}

TransitionProbabilities& TransitionProbabilities::operator=(TransitionProbabilities& c) {

    if (this != &c)
        {
        memcpy(p[0][0], c.p[0][0], tiSizeInBytes);
        }
    return *this;
}

void TransitionProbabilities::print(void) {

    for (int i=0; i<numStates; i++)
        {
        for (int c=0; c<numGammaCats; c++)
            {
            for (int j=0; j<numStates; j++)
                {
                std::cout << std::fixed << std::setprecision(4);
                std::cout << p[c][i][j] << " ";
                }
            std::cout << " | ";
            }
        std::cout << std::endl;
        }
}

void TransitionProbabilities::tiProbs(double v) {

    Eigen::Matrix<std::complex<double>, 4, 1> ceigenvalue = myModel->getEigenValues(mySpace);
    std::complex<double>* cc_ijk = myModel->getCijk(mySpace);
    ParameterAsrv* asrv = myModel->getParameterAsrv(mySpace);
    
    for (int c=0; c<numGammaCats; c++)
        {
        double r = (*asrv)[c];

        std::complex<double> ceigValExp[4];
        for (int s=0; s<4; s++)
            ceigValExp[s] = exp(ceigenvalue[s] * v * r);

        std::complex<double>* ptr = cc_ijk;
        for (int i=0; i<4; i++)
            {
            for (int j=0; j<4; j++)
                {
                std::complex<double> sum = std::complex<double>(0.0, 0.0);
                for(int s=0; s<4; s++)
                    sum += (*ptr++) * ceigValExp[s];
                p[c][i][j] = (sum.real() < 0.0) ? 0.0 : sum.real();
            }
        }
    }
}

