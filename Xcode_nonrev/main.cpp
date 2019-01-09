#include <iostream>
#include "Alignment.h"
#include "Mcmc.h"
#include "Model.h"
#include "RandomVariable.h"
#include "Settings.h"

void printHeader(void);



int main(int argc, char* argv[]) {

    // print nifty header
    printHeader();
    
    // get the user settings
    Settings mySettings(argc, argv);

    // make the random number object
    RandomVariable rv;

    // read the data matrix
    std::cout << &mySettings << std::endl;
    Alignment myAlignment( mySettings.getInputFileName() );
    //myAlignment.print();

    // set up the phylogenetic model
    Model myModel(&mySettings, &myAlignment, &rv);
    
    // run the Markov chain Monte Carlo analysis
    Mcmc myMcmc(&mySettings, &myModel, &rv, &myAlignment);
    myMcmc.runPowerPosterior();

    return 0;
}

void printHeader(void) {

    std::cout << std::endl;
    std::cout << "   Kasey's First Phylogeny Program, version 1.0" << std::endl;
    std::cout << std::endl;
    std::cout << "   * John Huelsenbeck and Kasey Arzumanova" << std::endl;
    std::cout << "   * University of California, Berkeley" << std::endl;
    std::cout << "   * " << std::endl;
}
