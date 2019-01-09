#ifndef Settings_H
#define Settings_H

#include <string>


class Settings {

    public:
                        Settings(int argc, char* argv[]);
        double          getAsrvLambda(void) { return asrvLambda; }
        double          getBrlenLambda(void) { return brlenLambda; }
        int             getBurninLength(void) { return burninLength; }
        int             getChainLength(void) { return chainLength; }
        std::string     getInputFileName(void) { return inputFileName; }
        bool            getIsReversible(void) { return isReversible; }
        int             getNumGammaCats(void) { return numGammaCats; }
        std::string     getOutPutFileName(void) { return outPutFileName; }
        int             getPreburninLength(void) { return preburninLength; }
        int             getPrintFrequency(void) { return printFrequency; }
        int             getSampleFrequency(void) { return sampleFrequency; }
        int             getSampleLength(void) { return sampleLength; }
        std::string     getTreeFileName(void) { return treeFileName; }
        int             getTuneLength(void) { return tuneLength; }
        void            setAsrvLambda(double x) { asrvLambda = x; }
        void            setInputFileName(std::string s) { inputFileName = s; }
        void            setIsReversible(bool tf) { isReversible = tf; }
        void            setNumGammaCats(int x) { numGammaCats = x; }
        void            setOutPutFileName(std::string s) { outPutFileName = s; }
        void            setTreeFileName(std::string s ) { treeFileName = s; }

    private:
        void            printUsage(void);
        double          brlenLambda;
        std::string     inputFileName;
        std::string     outPutFileName;
        std::string     treeFileName;
        double          asrvLambda;
        int             chainLength;
        int             preburninLength;
        int             tuneLength;
        int             burninLength;
        int             sampleLength;
        int             printFrequency;
        int             sampleFrequency;
        int             numGammaCats;
        bool            isReversible;
};

#endif
