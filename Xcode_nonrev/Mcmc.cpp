#include <cmath>
#include <iomanip>
#include <iostream>
#include "Alignment.h"
#include "Mcmc.h"
#include "Model.h"
#include "RandomVariable.h"
#include "Settings.h"

#include "ParameterTree.h"


Mcmc::Mcmc(Settings* s, Model* m, RandomVariable* r, Alignment* a) {

    mySettings = s;
    myModel = m;
    rv = r;
    myAlignment = a;
}

Mcmc::~Mcmc(void) {

}

void Mcmc::run(void) {

    parmFile = mySettings->getOutPutFileName() + ".p";
    treeFile = mySettings->getOutPutFileName() + ".t";
    parmOut.open( parmFile.c_str(), std::ios::out );
    treeOut.open( treeFile.c_str(), std::ios::out );

    burnin = mySettings->getBurninLength();
    int chainLength = mySettings->getChainLength();
    int printFreq   = mySettings->getPrintFrequency();
    int sampleFreq  = mySettings->getSampleFrequency();
    
    std::cout << "   * Running Markov chain Monte Carlo Analysis for " << chainLength << " generations" << std::endl;

    // get probabilities of initial state
    myModel->getParameterTree(0)->updateNodeFlags(true);
    myModel->getParameterTree(1)->updateNodeFlags(true);
    myModel->getParameterTree(0)->updateBranchFlags(true);
    myModel->getParameterTree(1)->updateBranchFlags(true);
    double curLnL = myModel->lnLikelihood(0);
    double curLnPrior = myModel->lnPrior(0);
    std::cout << "   * Initial log likelihood[1] = " << std::fixed << std::setprecision(3) << curLnL << std::endl;
    std::cout << "   * Initial log likelihood[2] = " << std::fixed << std::setprecision(3) << myModel->lnLikelihood(1) << std::endl;
    std::cout << "   * Initial log prior[1] = " << std::fixed << std::setprecision(3) << curLnPrior << std::endl;
    std::cout << "   * Initial log prior[2] = " << std::fixed << std::setprecision(3) << myModel->lnPrior(1) << std::endl;

    sampleChain(0, chainLength, curLnL, curLnPrior);
    
    auto t1 = std::chrono::high_resolution_clock::now();
    for (int n=1; n<=chainLength; n++)
        {
        // update parameter
        double lnProposalProb = myModel->update();
        
        // calculate the acceptance probability
        double newLnL = myModel->lnLikelihood(1);
        double newLnPrior = myModel->lnPrior(1);
        double lnR = (newLnL - curLnL) + (newLnPrior - curLnPrior) + lnProposalProb;
        double R = safeExponentiation(lnR);
        
        // print to screen (part 1)
        if (n % printFreq == 0)
            std::cout << "   * " << std::setw(5) << n << " -- " << std::fixed << std::setprecision(3) << curLnL << " -> " << newLnL;;
        
        // accept or reject the proposed state
        double u = rv->uniformRv();
        bool didAccept = false;
        if (u < R)
            {
            myModel->acceptMove();
            curLnL = newLnL;
            curLnPrior = newLnPrior;
            didAccept = true;
            }
        else
            {
            myModel->rejectMove();
            }
        
        // print to screen (part 2)
        if (n % printFreq == 0)
            {
            if (didAccept == true)
                std::cout << " (accepted " << myModel->getProposalName() << ")" << std::endl;
            else
                std::cout << " (rejected " << myModel->getProposalName() << ")" << std::endl;
            }

        // save the state of the chain to a file
        if (n % sampleFreq == 0 || n == chainLength)
            sampleChain(n, chainLength, curLnL, curLnPrior);
        
        }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t2 - t1;

    // result time
    std::cout << "   * MCMC took " << elapsed.count() << " seconds" << std::endl;

    myModel->summarizeProposals();

    parmOut.close();
    treeOut.close();
}

void Mcmc::runPowerPosterior(void) {

    // determine the powers that will be used
    int numStones = 127;
    std::vector<double> powers;
    double alpha = 0.3;
    double beta = 1.0;
    double intervalProb = (double)1.0 / numStones;
    powers.push_back(1.0);
    for (int i=numStones-1; i>0; i--)
        powers.push_back( rv->betaQuantile(alpha, beta, i * intervalProb) );
    powers.push_back(0.0);
    //for (int i=0; i<powers.size(); i++)
    //    std::cout << i+1 << " " << std::fixed << std::setprecision(20) << powers[i] << std::endl;

    // open the file for output of the power posterior information
    powerFile = mySettings->getOutPutFileName() + ".pp";
    powerOut.open( powerFile.c_str(), std::ios::out );

    // parameters of the chain
    int preburninLength = mySettings->getPreburninLength();
    int tuneLength      = mySettings->getTuneLength();
    int burninLength    = mySettings->getBurninLength();
    int sampleLength    = mySettings->getSampleLength();
    int printFreq       = mySettings->getPrintFrequency();
    int sampleFreq      = mySettings->getSampleFrequency();
    
    // get probabilities of initial state
    std::cout << "   * Running Power Posterior Analysis for " << preburninLength+tuneLength+burninLength+sampleLength << " generations for each of " << powers.size() << " stones" << std::endl;
    myModel->getParameterTree(0)->updateNodeFlags(true);
    myModel->getParameterTree(0)->updateBranchFlags(true);
    myModel->getParameterTree(1)->updateNodeFlags(true);
    myModel->getParameterTree(1)->updateBranchFlags(true);
    double curLnL = myModel->lnLikelihood(0);
    double curLnPrior = myModel->lnPrior(0);
    std::cout << "   * Initial log likelihood = " << std::fixed << std::setprecision(3) << curLnL << std::endl;
    //std::cout << "   * Initial log likelihood[2] = " << std::fixed << std::setprecision(3) << myModel->lnLikelihood(1) << std::endl;
    std::cout << "   * Initial log prior = " << std::fixed << std::setprecision(3) << curLnPrior << std::endl;
    //std::cout << "   * Initial log prior[2] = " << std::fixed << std::setprecision(3) << myModel->lnPrior(1) << std::endl;

    samplePower(0, sampleLength, sampleLength, curLnL, 0.0);
    
    auto t1 = std::chrono::high_resolution_clock::now();
    
    int sampleNum = 0;
    for (int pow=0; pow<powers.size(); pow++)
        {
        // set the power
        double power = powers[pow];
        
        // pre-burnin period
        for (int n=1; n<=preburninLength; n++)
            {
            // update parameter
            double lnProposalProb = myModel->update();
            
            // calculate the acceptance probability
            double newLnL = myModel->lnLikelihood(1);
            double newLnPrior = myModel->lnPrior(1);
            double lnR = power * (newLnL - curLnL) + (newLnPrior - curLnPrior) + lnProposalProb;
            double R = safeExponentiation(lnR);
            
            // print to screen (part 1)
            if (n % printFreq == 0)
                std::cout << "   * " << std::setw(5) << n << " -- " << std::fixed << std::setprecision(3) << curLnL << " -> " << newLnL << " -- " << power << " -- ";
            
            // accept or reject the proposed state
            double u = rv->uniformRv();
            bool didAccept = false;
            if (u < R)
                {
                myModel->acceptMove();
                curLnL = newLnL;
                curLnPrior = newLnPrior;
                didAccept = true;
                }
            else
                {
                myModel->rejectMove();
                }
            
            // print to screen (part 2)
            if (n % printFreq == 0)
                {
                std::cout << "pre-burnin period for stone " << pow+1 << "/" << powers.size();
                if (didAccept == true)
                    std::cout << " (accepted " << myModel->getProposalName() << ")" << std::endl;
                else
                    std::cout << " (rejected " << myModel->getProposalName() << ")" << std::endl;
                }
            }
            
        // tune period
        myModel->resetAcceptanceInformation();
        int numTunes = 0, tn = 0;
        while (numTunes < tuneLength)
            {
            tn++;
            
            // update parameter
            double lnProposalProb = myModel->update();
            
            // calculate the acceptance probability
            double newLnL = myModel->lnLikelihood(1);
            double newLnPrior = myModel->lnPrior(1);
            double lnR = power * (newLnL - curLnL) + (newLnPrior - curLnPrior) + lnProposalProb;
            double R = safeExponentiation(lnR);
            
            // print to screen (part 1)
            if (tn % printFreq == 0)
                std::cout << "   * " << std::setw(5) << tn << " -- " << std::fixed << std::setprecision(3) << curLnL << " -> " << newLnL << " -- " << power << " -- ";
            
            // accept or reject the proposed state
            double u = rv->uniformRv();
            bool didAccept = false;
            if (u < R)
                {
                myModel->acceptMove();
                curLnL = newLnL;
                curLnPrior = newLnPrior;
                didAccept = true;
                }
            else
                {
                myModel->rejectMove();
                }
            
            // print to screen (part 2)
            if (tn % printFreq == 0)
                {
                std::cout << "tuning period for stone " << pow+1 << "/" << powers.size() << " " << numTunes << " tunes";
                if (didAccept == true)
                    std::cout << " (accepted " << myModel->getProposalName() << ")" << std::endl;
                else
                    std::cout << " (rejected " << myModel->getProposalName() << ")" << std::endl;
                }
                
            // adjust parameters
            if (myModel->minimumNumberOfTries() >= 100)
                {
                myModel->adjustTuningParameters();
                myModel->resetAcceptanceInformation();
                numTunes++;
                }
            }

        // second burnin period
        for (int n=1; n<=burninLength; n++)
            {
            // update parameter
            double lnProposalProb = myModel->update();
            
            // calculate the acceptance probability
            double newLnL = myModel->lnLikelihood(1);
            double newLnPrior = myModel->lnPrior(1);
            double lnR = power * (newLnL - curLnL) + (newLnPrior - curLnPrior) + lnProposalProb;
            double R = safeExponentiation(lnR);
            
            // print to screen (part 1)
            if (n % printFreq == 0)
                std::cout << "   * " << std::setw(5) << n << " -- " << std::fixed << std::setprecision(3) << curLnL << " -> " << newLnL << " -- " << power << " -- ";
            
            // accept or reject the proposed state
            double u = rv->uniformRv();
            bool didAccept = false;
            if (u < R)
                {
                myModel->acceptMove();
                curLnL = newLnL;
                curLnPrior = newLnPrior;
                didAccept = true;
                }
            else
                {
                myModel->rejectMove();
                }
            
            // print to screen (part 2)
            if (n % printFreq == 0)
                {
                std::cout << "burnin period for stone " << pow+1 << "/" << powers.size();
                if (didAccept == true)
                    std::cout << " (accepted " << myModel->getProposalName() << ")" << std::endl;
                else
                    std::cout << " (rejected " << myModel->getProposalName() << ")" << std::endl;
                }
            }
            
        // sample chain
        for (int n=1; n<=sampleLength; n++)
            {
            // update parameter
            double lnProposalProb = myModel->update();
            
            // calculate the acceptance probability
            double newLnL = myModel->lnLikelihood(1);
            double newLnPrior = myModel->lnPrior(1);
            double lnR = power * (newLnL - curLnL) + (newLnPrior - curLnPrior) + lnProposalProb;
            double R = safeExponentiation(lnR);
            
            // print to screen (part 1)
            if (n % printFreq == 0)
                std::cout << "   * " << std::setw(5) << n << " -- " << std::fixed << std::setprecision(3) << curLnL << " -> " << newLnL << " -- " << power << " -- ";
            
            // accept or reject the proposed state
            double u = rv->uniformRv();
            bool didAccept = false;
            if (u < R)
                {
                myModel->acceptMove();
                curLnL = newLnL;
                curLnPrior = newLnPrior;
                didAccept = true;
                }
            else
                {
                myModel->rejectMove();
                }
            
            // print to screen (part 2)
            if (n % printFreq == 0)
                {
                std::cout << "sample period for stone " << pow+1 << "/" << powers.size();
                if (didAccept == true)
                    std::cout << " (accepted " << myModel->getProposalName() << ")" << std::endl;
                else
                    std::cout << " (rejected " << myModel->getProposalName() << ")" << std::endl;
                }

            // save the state of the chain to a file
            if ( n % sampleFreq == 0 || n == sampleLength )
                samplePower(++sampleNum, n, sampleLength, curLnL, power);
            }
            
        }
    
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t2 - t1;

    // result time
    std::cout << "   * MCMC took " << elapsed.count() << " seconds" << std::endl;

    powerOut.close();
}

double Mcmc::safeExponentiation(double lnX) {

    if (lnX < -300.0)
        return 0.0;
    else if (lnX > 0.0)
        return 1.0;
    return exp(lnX);
}

void Mcmc::sampleChain(int n, int len, double lnL, double lnP) {

    // parameter file
    if (n == 0)
        {
        char nucs[4] = { 'A', 'C', 'G', 'T' };
        parmOut << "Gen" << '\t';
        parmOut << "lnL" << '\t';
        parmOut << "lnP" << '\t';
        parmOut << "Shape" << '\t';
        for (int i=0; i<4; i++)
            {
            for (int j=0; j<4; j++)
                {
                if (i != j)
                    parmOut << "r[" << nucs[i] << "->" << nucs[j] << "]" << '\t';
                }
            }
        parmOut << "pi[A]" << '\t' << "pi[C]" << '\t' << "pi[G]" << '\t' << "pi[T]" << '\t';
        parmOut << "D" << '\t';
        parmOut << std::endl;
        }
    
    parmOut << n << '\t';
    parmOut << std::fixed << std::setprecision(4) << lnL << '\t';
    parmOut << std::fixed << std::setprecision(4) << lnP << '\t';
    parmOut << myModel->getGammaShape(1);
    parmOut << myModel->getSubstitutionRates(1);
    parmOut << myModel->getPiString(1);
    parmOut << myModel->getDStatistic(1) << '\t';
    parmOut << std::endl;

    // tree file
    if (n == 0)
        {
        std::vector<std::string> tnames = myAlignment->getTaxonNames();

        treeOut << "#NEXUS" << std::endl << std::endl;
        treeOut << "begin trees;" << std::endl;
        treeOut << "   translate" << std::endl;
        for (int i=0; i<tnames.size(); i++)
            {
            treeOut << "   " << i+1 << " " << tnames[i];
            if (i+1 != tnames.size())
                treeOut << ",";
            treeOut << std::endl;
            }
        treeOut << "   ;" << std::endl;
        }
    
    treeOut << "   tree sample_" << n << " = " << myModel->getTreeString(1) << ";" << std::endl;
    
    if (n == len)
        {
        treeOut << "end;" << std::endl << std::endl;
        }
}

void Mcmc::samplePower(int sn, int n, int len, double lnL, double power) {

    // parameter file
    if (n == 0)
        {
        powerOut << "num" << '\t';
        powerOut << "sample" << '\t';
        powerOut << "state" << '\t';
        powerOut << "power" << '\t';
        powerOut << "likelihood" << '\t';
        powerOut << std::endl;
        }
    else
        {
        powerOut << sn << '\t';
        powerOut << n << '\t';
        powerOut << std::fixed << std::setprecision(20) << power << '\t';
        powerOut << std::fixed << std::setprecision(4) << lnL << '\t';
        powerOut << std::endl;
        }
}
