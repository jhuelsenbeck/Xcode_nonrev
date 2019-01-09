#include <iostream>
#include <cmath>
#include <string>
#include "Settings.h"



Settings::Settings(int argc, char* argv[]) {

#   if 1
    /* set up fake command-line argument string */
    char *cmdString[33];
    cmdString[ 0] = (char*)"ck";
    cmdString[ 1] = (char*)"-i";
    //cmdString[ 2] = (char*)"/Users/johnh/Repositories/nonrev/data/bglobin/bglobin.in";
    cmdString[ 2] = (char*)"/Users/johnh/Repositories/nonrev/data/allen07/allen2007.in";
    cmdString[ 3] = (char*)"-o";
    //cmdString[ 4] = (char*)"/Users/johnh/Repositories/nonrev/data/output/bglobin.out";
    cmdString[ 4] = (char*)"/Users/johnh/Repositories/nonrev/data/output/allen2007.out";
    cmdString[ 5] = (char*)"-t";
    //cmdString[ 6] = (char*)"/Users/johnh/Repositories/nonrev/data/bglobin/bglobin_altnex.trees";
    cmdString[ 6] = (char*)"/Users/johnh/Repositories/nonrev/data/allen07/allen2007.trees";
    cmdString[ 7] = (char*)"-tr";
    cmdString[ 8] = (char*)"no";
    cmdString[ 9] = (char*)"-pb";
    cmdString[10] = (char*)"10000";
    cmdString[11] = (char*)"-tt";
    cmdString[12] = (char*)"10000";
    cmdString[13] = (char*)"-d";
    cmdString[14] = (char*)"10000";
    cmdString[15] = (char*)"-tt";
    cmdString[16] = (char*)"20";
    cmdString[17] = (char*)"-st";
    cmdString[18] = (char*)"10000";
    cmdString[19] = (char*)"-p";
    cmdString[20] = (char*)"100";
    cmdString[21] = (char*)"-s";
    cmdString[22] = (char*)"100";
    cmdString[23] = (char*)"-lambda";
    cmdString[24] = (char*)"12.2";
    cmdString[25] = (char*)"-g";
    cmdString[26] = (char*)"4";
    cmdString[27] = (char*)"-e";
    cmdString[28] = (char*)"2.0";
    argc = 29;
    argv = cmdString;
#   endif

    enum Mode { DATA_FILE, TREE_FILE, OUTPUT_FILE, CHAIN_LENGTH, PREBURNIN_LENGTH, TUNE_LENGTH, BURNIN_LENGTH, SAMPLE_LENGTH, PRINT_FREQ, SAMPLE_FREQ, BRLEN_LAMBDA, NUM_GAMMA_CATS, ASRV_LAMBDA, IS_REVERSIBLE, TUNE_T1, TUNE_T2, TUNE_T3, TUNE_T4, TUNE_T5, NONE };

    /* set default values for parameters */
    inputFileName          = "";
    treeFileName           = "";
    outPutFileName         = "";
    preburninLength        = 10000;
    tuneLength             = 10000;
    burninLength           = 100000;
    sampleLength           = 100000;
    chainLength            = 1000000;
    printFrequency         = 100;
    sampleFrequency        = 100;
    brlenLambda            = 10.0;
    numGammaCats           = 4;
    asrvLambda             = 2.0;
    isReversible           = true;
    tuningParm[0]          = 500.0;
    tuningParm[1]          = log(1.4);
    tuningParm[2]          = 300.0;
    tuningParm[3]          = 300.0;
    tuningParm[4]          = log(1.2);
    
    if (argc > 1)
        {
        if (argc % 2 == 0)
            {
            printUsage();
            }
            
        /* read the command-line arguments */
        int status = NONE;
        for (int i=1; i<argc; i++)
            {
            std::string cmd = argv[i];
            //std::cout << cmd << std::endl;
            if (status == NONE)
                {
                /* read the parameter specifier */
                if ( cmd == "-i" )
                    status = DATA_FILE;
                else if ( cmd == "-t" )
                    status = TREE_FILE;
                else if ( cmd == "-o" )
                    status = OUTPUT_FILE;
                else if ( cmd == "-l" )
                    status = CHAIN_LENGTH;
                else if ( cmd == "-pb" )
                    status = PREBURNIN_LENGTH;
                else if ( cmd == "-tt" )
                    status = TUNE_LENGTH;
                else if ( cmd == "-d" )
                    status = BURNIN_LENGTH;
                else if ( cmd == "-st" )
                    status = SAMPLE_LENGTH;
                else if ( cmd == "-p" )
                    status = PRINT_FREQ;
                else if ( cmd == "-s" )
                    status = SAMPLE_FREQ;
                else if ( cmd == "-tr" )
                    status = IS_REVERSIBLE;
                else if ( cmd == "-lambda" )
                    status = BRLEN_LAMBDA;
                else if ( cmd == "-g" )
                    status = NUM_GAMMA_CATS;
                else if ( cmd == "-e" )
                    status = ASRV_LAMBDA;
                else if ( cmd == "-t1" )
                    status = TUNE_T1;
                else if ( cmd == "-t2" )
                    status = TUNE_T2;
                else if ( cmd == "-t3" )
                    status = TUNE_T3;
                else if ( cmd == "-t4" )
                    status = TUNE_T4;
                else if ( cmd == "-t5" )
                    status = TUNE_T5;
                else
                    {
                    std::cerr << "Could not interpret option \"" << cmd << "\"." << std::endl;
                    exit(1);
                    }
                }
            else
                {
                /* read the parameter */
                if ( status == DATA_FILE )
                    inputFileName = argv[i];
                else if ( status == OUTPUT_FILE )
                    outPutFileName = argv[i];
                else if ( status == TREE_FILE )
                    treeFileName = argv[i];
                else if ( status == CHAIN_LENGTH )
                    chainLength = atoi(argv[i]);
                else if ( status == PREBURNIN_LENGTH )
                    preburninLength = atoi(argv[i]);
                else if ( status == TUNE_LENGTH )
                    tuneLength = atoi(argv[i]);
                else if ( status == BURNIN_LENGTH )
                    burninLength = atoi(argv[i]);
                else if ( status == SAMPLE_LENGTH )
                    sampleLength = atoi(argv[i]);
                else if ( status == PRINT_FREQ )
                    printFrequency = atoi(argv[i]);
                else if ( status == SAMPLE_FREQ )
                    sampleFrequency = atoi(argv[i]);
                else if ( status == IS_REVERSIBLE )
                    {
                    if (argv[i][0] == 'Y' || argv[i][0] == 'y')
                        isReversible = true;
                    else if (argv[i][0] == 'N' || argv[i][0] == 'n')
                        isReversible = false;
                    else
                        printUsage();
                    }
                else if ( status == BRLEN_LAMBDA )
                    brlenLambda = atof(argv[i]);
                else if ( status == NUM_GAMMA_CATS )
                    numGammaCats = atoi(argv[i]);
                else if ( status == ASRV_LAMBDA )
                    asrvLambda = atof(argv[i]);
                else if ( status == TUNE_T1 )
                    tuningParm[0] = atof(argv[i]);
                else if ( status == TUNE_T2 )
                    tuningParm[1] = atof(argv[i]);
                else if ( status == TUNE_T3 )
                    tuningParm[2] = atof(argv[i]);
                else if ( status == TUNE_T4 )
                    tuningParm[3] = atof(argv[i]);
                else if ( status == TUNE_T5 )
                    tuningParm[4] = atof(argv[i]);
                else
                    {
                    std::cerr << "Unknown status reading command line information" << std::endl;
                    exit(1);
                    }
                status = NONE;
                }
            }
        }
    else
        {
        printUsage();
        }
}

double Settings::getTuningParm(std::string parmNameStr) {

    if ( parmNameStr == "tree" )
        return tuningParm[0];
    else if ( parmNameStr == "asrv" )
        return tuningParm[1];
    else if ( parmNameStr == "freqs" )
        return tuningParm[2];
    else if ( parmNameStr == "rates" )
        return tuningParm[3];
    else if ( parmNameStr == "length" )
        return tuningParm[4];
    else
        {
        std::cerr << "ERROR: Unknown parameter " << parmNameStr << std::endl;
        exit(1);
        }
    return 0.0;
}

void Settings::printUsage(void) {

    std::cout << "Usage:" << std::endl;
    std::cout << "     -i : Input file name" << std::endl;
    std::cout << "     -t : Tree file name (for constraining the analysis to a fixed tree)" << std::endl;
    std::cout << "     -o : Output file name" << std::endl;
    std::cout << "     -l : Number of MCMC cycles" << std::endl;
    std::cout << "     -d : Number of MCMC cycles to discard as the \"burn-in\"" << std::endl;
    std::cout << "     -p : Print frequency" << std::endl;
    std::cout << "     -s : Sample frequency" << std::endl;
    std::cout << "    -tr : Is the model time reversible (yes or no)" << std::endl;
    std::cout << "-lambda : Exponential parameter for branch lengths" << std::endl;
    std::cout << "     -g : Number of gamma rate categories" << std::endl;
    std::cout << "     -e : Exponential parameter for shape parameter describing ASRV" << std::endl;
    std::cout << "    -t1 : MCMC tuning parameter for the tree topology parameter" << std::endl;
    std::cout << "    -t2 : MCMC tuning parameter for the gamma shape parameter" << std::endl;
    std::cout << "    -t3 : MCMC tuning parameter for the base frequencies" << std::endl;
    std::cout << "    -t4 : MCMC tuning parameter for the substitution rates" << std::endl;
    std::cout << "    -t5 : MCMC tuning parameter for the tree length parameter" << std::endl;
    std::cout << std::endl;
    std::cout << "Example:" << std::endl;
    std::cout << "   ./kaseysprg -i <input file> -o <output file>" << std::endl;
    exit(0);
}

