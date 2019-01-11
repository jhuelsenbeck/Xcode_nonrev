#ifndef Model_H
#define Model_H

#include <complex>
#include <string>
#include <vector>
#include <Eigen/Core>

#undef RUN_ON_PRIOR

typedef Eigen::Matrix<double, 4, 1> NucleotideVector_t;
typedef Eigen::Matrix<double, 4, 4> NucleotideSquareMatrix_t;

class Alignment;
class ConditionalLikelihoods;
class Parameter;
class ParameterTree;
class ParameterAsrv;
class ParameterBaseFrequencies;
class ParameterExchangabilityRates;
class ProposalInfo;
class RandomVariable;
class Settings;
class TransitionProbabilities;
class ParameterTree;


class Model {

    public:
                                                        Model(Settings* s, Alignment* a, RandomVariable* r);
                                                       ~Model(void);
        void                                            acceptMove(void);
        void                                            adjustTuningParameters(void);
        double                                          lnLikelihood(int whichSpace);
        double                                          lnPrior(int whichSpace);
        Eigen::Matrix<std::complex<double>, 4, 1>&      getEigenValues(int whichSpace) { return eigenValues[whichSpace]; }
        std::complex<double>*                           getCijk(int whichSpace) { return cc_ijk[whichSpace]; }
        double                                          getDStatistic(int whichSpace) { return dStatistic[whichSpace]; }
        std::string                                     getGammaShape(int whichSpace);
        std::vector<double>                             getOrderedBranchLengths(int whichSpace);
        ParameterBaseFrequencies*                       getParameterBaseFrequencies(int whichSpace);
        ParameterTree*                                  getParameterTree(int whichSpace);
        ParameterExchangabilityRates*                   getParameterExchangabilityRates(int whichSpace);
        ParameterAsrv*                                  getParameterAsrv(int whichSpace);
        std::vector<double>                             getPi(int whichSpace) { return pi[whichSpace]; }
        std::string                                     getProposalName(void);
        std::string                                     getSubstitutionRates(int whichSpace);
        std::string                                     getTreeString(int whichSpace);
        std::string                                     getPiString(int whichSpace);
        ParameterTree*                                  getTree(int whichSpace);
        int                                             minimumNumberOfTries(void);
        void                                            printTransitionProbabilities(int whichSpace);
        void                                            rejectMove(void);
        void                                            resetAcceptanceInformation(void);
        void                                            summarizeProposals(void);
        double                                          update(void);
        
    private:
        std::vector<double>                             calcStationaryFreq(NucleotideSquareMatrix_t& Q);
        double                                          updateBaseFrequencies(void);
        void                                            updateRateMatrix(int whichSpace);
        double                                          updateGammaShapeParm(void);
        double                                          updateSubstitutionRatesParm(void);
        double                                          updateBranchLengthParm(void);
        double                                          updateRootPosition(void);
        double                                          updateTreeParmSpr(void);
        double                                          updateTreeParmLocal(void);
        void                                            updateTiProbs(int whichSpace);
    
        Alignment*                                      alnPtr;
        Settings*                                       settingsPtr;
        RandomVariable*                                 rv;
    
        std::vector<Parameter*>                         parameters[2];
        Parameter*                                      updatedParameters[2];
        std::vector<ProposalInfo*>                      proposals;

        std::vector<ConditionalLikelihoods*>            cls[2];
        TransitionProbabilities**                       tis[2];
        std::vector<double>                             pi[2];
        Eigen::Matrix<std::complex<double>, 4, 1>       eigenValues[2];
        std::complex<double>*                           cc_ijk[2];
    
        double                                          dStatistic[2];
        double                                          gammaShapePrior;
        int                                             numSites;
        int                                             numGammaCats;
        int                                             lastProposal;
        int                                             numNodes;
        bool                                            isTimeReversible;
};

#endif
