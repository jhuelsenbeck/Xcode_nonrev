#include "ProposalInfo.h"


ProposalInfo::ProposalInfo(double w, double p, std::string s) {

    weight = w;
    probability = p;
    name = s;
    numAttempts = 0;
    numAcceptances = 0;
}

double ProposalInfo::acceptanceFrequency(void) {

    return (double) numAcceptances / numAttempts;
}

void ProposalInfo::accept(void) {

    numAttempts++;
    numAcceptances++;
}

void ProposalInfo::reject(void) {

    numAttempts++;
}

void ProposalInfo::reset(void) {

    numAttempts = 0;
    numAcceptances = 0;
}
