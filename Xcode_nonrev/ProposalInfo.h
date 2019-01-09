#ifndef ProposalInfo_H
#define ProposalInfo_H

#include <string>

class ProposalInfo {

    public:
                    ProposalInfo(void) = delete;
                    ProposalInfo(double w, double p, std::string s);
        void        accept(void);
        double      getProbability(void) { return probability; }
        double      getWeight(void) { return weight; }
        std::string getName(void) { return name; }
        int         getNumAcceptances(void) { return numAcceptances; }
        int         getNumAttempts(void) { return numAttempts; }
        void        reject(void);
        void        reset(void);
        void        setProbability(double x) { probability = x; }
        void        setWeight(double x) { weight = x; }
        double      acceptanceFrequency(void);
    
    protected:
        double      weight;
        double      probability;
        std::string name;
        int         numAttempts;
        int         numAcceptances;
};

#endif
