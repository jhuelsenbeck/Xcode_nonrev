#include <cmath>
#include <iomanip>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include "Alignment.h"
#include "ConditionalLikelihoods.h"
#include "Model.h"
#include "Msg.h"
#include "Node.h"
#include "ParameterAsrv.h"
#include "ParameterBaseFrequencies.h"
#include "ParameterExchangabilityRates.h"
#include "ParameterTree.h"
#include "ProposalInfo.h"
#include "RandomVariable.h"
#include "Settings.h"
#include "TransitionProbabilities.h"

#undef USE_SSE

#if defined(USE_SSE)
#   include <xmmintrin.h>
#   include <immintrin.h>
#endif

Model::Model(Settings* s, Alignment* a, RandomVariable* r) {

    // remember the locations of some important objects
    alnPtr      = a;
    settingsPtr = s;
    rv          = r;
    
    // initialize some instance variables
    isTimeReversible = settingsPtr->getIsReversible();
    numSites = alnPtr->getNumChar();
    numGammaCats = settingsPtr->getNumGammaCats();
    numNodes = 2 * alnPtr->getNumTaxa() - 1;

    // initialize the conditional likelihoods
    std::cout << "   * Initializing conditional likelihoods" << std::endl;
    for (int i=0; i<numNodes; i++)
        {
        cls[0].push_back( new ConditionalLikelihoods(alnPtr, numGammaCats) );
        cls[1].push_back( new ConditionalLikelihoods(alnPtr, numGammaCats) );
        }
    for (int i=0; i<alnPtr->getNumTaxa(); i++)
        {
        cls[0][i]->initializeTipConditonalLikelihoods(alnPtr, alnPtr->getTaxonNames()[i]);
        cls[1][i]->initializeTipConditonalLikelihoods(alnPtr, alnPtr->getTaxonNames()[i]);
        }

    // instantiate the tree parameter of the model
    std::cout << "   * Initializing parameters" << std::endl;
    ParameterTree* t0 = NULL;
    if (settingsPtr->getTreeFileName() == "")
        t0 = new ParameterTree(rv, this, "Tree", alnPtr->getTaxonNames(), true, settingsPtr->getBrlenLambda(), 0.0, isTimeReversible);
    else
        t0 = new ParameterTree(rv, this, "Tree", alnPtr->getTaxonNames(), settingsPtr->getTreeFileName(), true, settingsPtr->getBrlenLambda(), 0.0, isTimeReversible);
    ParameterTree* t1 = new ParameterTree(*t0);
    parameters[0].push_back(t0);
    parameters[1].push_back(t1);
    numNodes = (int)t0->getNumNodes();
    if (numNodes != t0->getNumNodes())
        Msg::error("Tree is not the right size for the alignment!");

    // initialize the base frequencies, if time reversible
    if (isTimeReversible == true)
        {
        ParameterBaseFrequencies* f0 = new ParameterBaseFrequencies(rv, this, "Base Frequencies", 100.0);
        ParameterBaseFrequencies* f1 = new ParameterBaseFrequencies(*f0);
        parameters[0].push_back(f0);
        parameters[1].push_back(f1);
        }

    // initialize the substituton rate parameters
    std::cout << "   * Initializing rate matrix" << std::endl;
    cc_ijk[0] = new std::complex<double>[64];
    cc_ijk[1] = new std::complex<double>[64];
    ParameterExchangabilityRates* e0 = new ParameterExchangabilityRates(rv, this, "Exchangabilities", 1000.0, isTimeReversible);
    ParameterExchangabilityRates* e1 = new ParameterExchangabilityRates(*e0);
    parameters[0].push_back(e0);
    parameters[1].push_back(e1);
    pi[0].resize(4);
    pi[1].resize(4);
#   if 0
    for (int i=0; i<substitutionRates[0].size(); i++)
        substitutionRates[0][i] = 1.0;
    substitutionRates[1] = substitutionRates[0];
    for (int i=0; i<substitutionRates[0].size(); i++)
        std::cout << i << " -- " << substitutionRates[0][i] << " " << substitutionRates[1][i] << std::endl;
#   endif
    updateRateMatrix(0);
    updateRateMatrix(1);
    
    // initialize the gamma shape parameter
    gammaShapePrior = settingsPtr->getAsrvLambda();
    ParameterAsrv* s0 = new ParameterAsrv(rv, this, "Gamma Shape", std::log(1.1), numGammaCats, gammaShapePrior);
    ParameterAsrv* s1 = new ParameterAsrv(*s0);
    parameters[0].push_back(s0);
    parameters[1].push_back(s1);

    // initialize the transition probabilities
    std::cout << "   * Initializing transition probabilities" << std::endl;
    for (int i=0; i<2; i++)
        {
        // allocate the transition probabilities for each branch of the tree
        tis[i] = new TransitionProbabilities*[numNodes];
        for (int n=0; n<numNodes; n++)
            tis[i][n] = new TransitionProbabilities(i, numGammaCats, this);
            
        // initialize the transition probabilities
        std::vector<Node*> dps = t0->getDownPassSequence();
        for (int n=0; n<dps.size(); n++)
            {
            Node* p = dps[n];
            if (p->getAncestor() != NULL)
                {
                Branch* b = t0->findBranch(p, p->getAncestor());
                double v = b->getLength();
                TransitionProbabilities* tp = tis[i][p->getIndex()];
                tp->tiProbs(v);
                //tp->print();
                }
            }
        }
    
    // initialize the proposals
    std::cout << "   * Initializing proposal weights" << std::endl;
    proposals.push_back( new ProposalInfo( 5.0, 0.0, "gamma shape update") );
    proposals.push_back( new ProposalInfo(10.0, 0.0, "substitution rates update") );
    if (settingsPtr->getTreeFileName() == "")
        {
        proposals.push_back( new ProposalInfo(10.0, 0.0, "branch length update") );
        proposals.push_back( new ProposalInfo(10.0, 0.0, "tree and branch lengths update (SPR)") );
        proposals.push_back( new ProposalInfo(10.0, 0.0, "tree and branch lengths update (LOCAL)") );
        if (isTimeReversible == false)
            proposals.push_back( new ProposalInfo( 0.0, 0.0, "root position update") );
        }
    else
        {
        proposals.push_back( new ProposalInfo(75.0, 0.0, "branch length update") );
        proposals.push_back( new ProposalInfo( 0.0, 0.0, "tree and branch lengths update (SPR)") );
        proposals.push_back( new ProposalInfo( 0.0, 0.0, "tree and branch lengths update (LOCAL)") );
        if (isTimeReversible == false)
            proposals.push_back( new ProposalInfo( 5.0, 0.0, "root position update") );
        }
    if (isTimeReversible == true)
        proposals.push_back( new ProposalInfo( 5.0, 0.0, "base frequency update") );
    else
        proposals.push_back( new ProposalInfo( 0.0, 0.0, "base frequency update") );

    // rescale weights
    if (numGammaCats == 1)
        proposals[0]->setWeight(0.0);
    double sumWeight = 0.0;
    for (int i=0; i<proposals.size(); i++)
        sumWeight += proposals[i]->getWeight();
    for (int i=0; i<proposals.size(); i++)
        proposals[i]->setProbability( proposals[i]->getWeight() / sumWeight );
}

Model::~Model(void) {
    
    for (int i=0; i<2; i++)
        {
        for (int n=0; n<numNodes; n++)
            delete tis[i][n];
        delete [] tis[i];
        }
    for (int i=0; i<numNodes; i++)
        {
        delete cls[0][i];
        delete cls[1][i];
        }
    delete [] cc_ijk[0];
    delete [] cc_ijk[1];
    for (int i=0; i<proposals.size(); i++)
        delete proposals[i];
    for (int i=0; i<parameters[0].size(); i++)
        {
        delete parameters[0][i];
        delete parameters[1][i];
        }
}

void Model::acceptMove(void) {

    proposals[lastProposal]->accept();
    
    // copy the parameters
    (*updatedParameters[0]) = (*updatedParameters[1]);
    
    // update the transition probability information (TEMP: This can go later)
    /*for (int i=0; i<numNodes; i++)
        {
        (*cls[0][i]) = (*cls[1][i]);
        (*tis[0][i]) = (*tis[1][i]);
        }*/

    // update other information, including the tree if necessary
    bool properlyUpdated = false;
    
    ParameterBaseFrequencies* p1 = dynamic_cast<ParameterBaseFrequencies*>(updatedParameters[0]);
    if (p1 != NULL)
        {
        ParameterTree* t0 = getParameterTree(0);
        ParameterTree* t1 = getParameterTree(1);
        (*t0) = (*t1);
        
        eigenValues[0] = eigenValues[1];
        pi[0] = pi[1];
        dStatistic[0] = dStatistic[1];
        for (int i=0; i<64; i++)
            cc_ijk[0][i] = cc_ijk[1][i];
        properlyUpdated = true;
        }

    ParameterExchangabilityRates* p2 = dynamic_cast<ParameterExchangabilityRates*>(updatedParameters[0]);
    if (p2 != NULL)
        {
        ParameterTree* t0 = getParameterTree(0);
        ParameterTree* t1 = getParameterTree(1);
        (*t0) = (*t1);

        eigenValues[0] = eigenValues[1];
        pi[0] = pi[1];
        dStatistic[0] = dStatistic[1];
        for (int i=0; i<64; i++)
            cc_ijk[0][i] = cc_ijk[1][i];
        properlyUpdated = true;
        }

    ParameterAsrv* p3 = dynamic_cast<ParameterAsrv*>(updatedParameters[0]);
    if (p3 != NULL)
        {
        ParameterTree* t0 = getParameterTree(0);
        ParameterTree* t1 = getParameterTree(1);
        (*t0) = (*t1);

        properlyUpdated = true;
        }

    ParameterTree* p4 = dynamic_cast<ParameterTree*>(updatedParameters[0]);
    if (p4 != NULL)
        {
        // no need to copy the tree as it was already done, above
        
        properlyUpdated = true;
        }
    
    if (properlyUpdated == false)
        Msg::error("Problem squaring parameters in acceptMove");
}

void Model::adjustTuningParameters(void) {

    for (int i=0; i<proposals.size(); i++)
        {
        if (proposals[i]->getProbability() < 0.000001)
            continue;
            
        std::string proposalName = proposals[i]->getName();
        double acceptanceRate = proposals[i]->acceptanceFrequency();

        if (proposalName == "gamma shape update")
            {
            ParameterAsrv* r0 = getParameterAsrv(0);
            ParameterAsrv* r1 = getParameterAsrv(1);
            double p = 0.44;
            double oldTuning = r0->getTuning();
            double newTuning = oldTuning;
            if (acceptanceRate > p)
                newTuning *= (1.0 + ((acceptanceRate-p) / (1.0-p)) );
            else
                newTuning /= (2.0 - acceptanceRate / p);
            if (newTuning < log(1.01) )
                newTuning = log(1.01);
            else if (newTuning > log(1000.0))
                newTuning = log(1000.0);
            r0->setTuning(newTuning);
            r1->setTuning(newTuning);
            std::cout << "   * Adjusting gamma shape tuning parameter from " << oldTuning << " to " << newTuning << " (" << acceptanceRate << ")" << std::endl;
            }
        else if (proposalName == "branch length update")
            {
            ParameterTree* t0 = getParameterTree(0);
            ParameterTree* t1 = getParameterTree(1);
            double p = 0.44;
            double oldTuning = t0->getBranchTuning();
            double newTuning = oldTuning;
            if (acceptanceRate > p)
                newTuning *= (1.0 + ((acceptanceRate-p) / (1.0-p)) );
            else
                newTuning /= (2.0 - acceptanceRate / p);
            if (newTuning < log(1.01) )
                newTuning = log(1.01);
            else if (newTuning > log(1000.0))
                newTuning = log(1000.0);
            t0->setBranchTuning(newTuning);
            t1->setBranchTuning(newTuning);
            std::cout << "   * Adjusting branch length tuning parameter from " << oldTuning << " to " << newTuning << " (" << acceptanceRate << ")" << std::endl;
            }
        else if (proposalName == "substitution rates update")
            {
            ParameterExchangabilityRates* e0 = getParameterExchangabilityRates(0);
            ParameterExchangabilityRates* e1 = getParameterExchangabilityRates(1);
            double p = 0.22;
            double oldTuning = e0->getTuning();
            double newTuning = oldTuning;
            if (acceptanceRate > p)
                newTuning /= (1.0 + ((acceptanceRate-p) / (1.0-p)) );
            else
                newTuning *= (2.0 - acceptanceRate / p);
            if (newTuning < e0->getExchangabilityAlpha().size() )
                newTuning = e0->getExchangabilityAlpha().size();
            else if (newTuning > 10000.0)
                newTuning = 10000.0;
            e0->setTuning(newTuning);
            e1->setTuning(newTuning);
            std::cout << "   * Adjusting exchangability rates tuning parameter from " << oldTuning << " to " << newTuning << " (" << acceptanceRate << ")" << std::endl;
            }
        else if (proposalName == "base frequency update")
            {
            ParameterBaseFrequencies* b0 = getParameterBaseFrequencies(0);
            ParameterBaseFrequencies* b1 = getParameterBaseFrequencies(1);
            double p = 0.22;
            double oldTuning = b0->getTuning();
            double newTuning = oldTuning;
            if (acceptanceRate > p)
                newTuning /= (1.0 + ((acceptanceRate-p) / (1.0-p)) );
            else
                newTuning *= (2.0 - acceptanceRate / p);
            if (newTuning < 4.0 )
                newTuning = 4.0;
            else if (newTuning > 10000.0)
                newTuning = 10000.0;
            b0->setTuning(newTuning);
            b1->setTuning(newTuning);
            std::cout << "   * Adjusting base frequencies tuning parameter from " << oldTuning << " to " << newTuning << " (" << acceptanceRate << ")" << std::endl;
            }
        else if (proposalName == "root position update")
            {
            
            }
        else if (proposalName == "tree and branch lengths update (SPR)")
            {
            
            }
        else if (proposalName == "tree and branch lengths update (LOCAL)")
            {
            
            }
        }
}

/*
 * This function calculates the stationary frequencies of the rate matrix. The
 * rate matrix, Q, is the infinitesimal generator of the Markov chain. It is an
 * n X n matrix whose off-diagonal elements are q_ij >= 0 and whose diagonal elements
 * are specified such that each row sums to zero. The rate matrix is finite (has
 * a fixed number of states) and we assume that the input matrix is irreducible, as
 * is the usual case for substitution models. Because Q is irreducible and finite,
 * it has a stationary distribution, pi, which is a row vector of n probabilities.
 * The stationary probabilities can be calculated by solving the homogeneous system
 * of equations, pi*Q = 0, where 0 is a vector of zeros.
 *
 * We do the following to calculate the stationary frequencies.
 *
 * 1. We perform an LU decomposition of the transpose of the matrix Q.
 *
 *    Q' = LU
 *
 * 2. Now we set Ux = z (x will eventually hold the stationary probabilities).
 *    Because L is nonsingular, we have z = 0. We proceed to back substitute on
 *    Ux = z = 0. When u_nn = 0, we can put in any solution for x. Here, we put
 *    in x_n = 1. We then solve the other values of x through back substitution.
 *
 * 3. The solution obtained in 2 is not a probability vector. We normalize the
 *    vector such that the sum of the elements is 1.
 *
 * Note that the only time we need to use this function is when we don't
 * know the stationary frequencies of the rate matrix beforehand. For most
 * substitution models used in molecular evolution, the stationary frequencies
 * are built into the rate matrix itself. These models are time-reversible.
 * This function is useful for the non-reversible models.
 *
 * \brief Calculate stationary frequencies of non-reversible model.
 *
 * \see
 * Stewart, W. J. 1999. Numerical methods for computing stationary distributions of
 *    finite irreducible Markov chains. In "Advances in Computational
 *    Probability", W. Grassmann, ed. Kluwer Academic Publishers.
 *
 */
std::vector<double> Model::calcStationaryFreq(NucleotideSquareMatrix_t& Q) {
    
    Eigen::VectorXd f = Q.transpose().fullPivLu().kernel();
    f = f / f.sum();
    //std::cout << f << std::endl;
    
    std::vector<double> stationaryFrequencies(4);
    for (int i=0; i<4; i++)
        stationaryFrequencies[i] = f(i);
    return stationaryFrequencies;
}

std::string Model::getProposalName(void) {

    return proposals[lastProposal]->getName();
}

std::string Model::getGammaShape(int whichSpace) {

    ParameterAsrv* asrv = getParameterAsrv(whichSpace);
    char cStr[50];
    sprintf(cStr, "%1.4lf", asrv->getAsrv());
    std::string str = cStr;
    str += '\t';
    return str;
}

ParameterBaseFrequencies* Model::getParameterBaseFrequencies(int whichSpace) {

    for (Parameter* parm : parameters[whichSpace])
        {
        ParameterBaseFrequencies* p = dynamic_cast<ParameterBaseFrequencies*>(parm);
        if (p != NULL)
            return p;
        }
    return NULL;
}

ParameterTree* Model::getParameterTree(int whichSpace) {

    for (Parameter* parm : parameters[whichSpace])
        {
        ParameterTree* p = dynamic_cast<ParameterTree*>(parm);
        if (p != NULL)
            return p;
        }
    return NULL;
}

ParameterExchangabilityRates* Model::getParameterExchangabilityRates(int whichSpace) {

    for (Parameter* parm : parameters[whichSpace])
        {
        ParameterExchangabilityRates* p = dynamic_cast<ParameterExchangabilityRates*>(parm);
        if (p != NULL)
            return p;
        }
    return NULL;
}

ParameterAsrv* Model::getParameterAsrv(int whichSpace) {

    for (Parameter* parm : parameters[whichSpace])
        {
        ParameterAsrv* p = dynamic_cast<ParameterAsrv*>(parm);
        if (p != NULL)
            return p;
        }
    return NULL;
}

std::string Model::getPiString(int whichSpace) {

    char cStr[50];
    std::string str = "";
    for (int i=0; i<4; i++)
        {
        sprintf(cStr, "%1.4lf", pi[whichSpace][i]);
        str += cStr;
        str += '\t';
        }
    return str;
}

std::string Model::getSubstitutionRates(int whichSpace) {

    ParameterExchangabilityRates* sr = getParameterExchangabilityRates(whichSpace);
    char cStr[50];
    std::string str = "";
    for (int i=0; i<12; i++)
        {
        sprintf(cStr, "%1.4lf", (*sr)[i]);
        str += cStr;
        str += '\t';
        }
    return str;
}

ParameterTree* Model::getTree(int whichSpace) {

    for (int i=0; i<parameters[whichSpace].size(); i++)
        {
        Parameter* parm = parameters[whichSpace][i];
        ParameterTree* dc = dynamic_cast<ParameterTree*>(parm);
        if (dc != NULL)
            return dc;
        }
    return NULL;
}

std::string Model::getTreeString(int whichSpace) {

    ParameterTree* t = getParameterTree(whichSpace);
    return t->getNewick();
}

# if defined(USE_SSE)

double Model::lnLikelihood(int whichSpace) {

#   if defined(RUN_ON_PRIOR)
    return 0.0;
#   endif

    // get the tree
    ParameterTree* t = getParameterTree(whichSpace);
    
    // get the downpass sequence
    std::vector<Node*> dps = t->getDownPassSequence();
    
    // fill in the conditional likelihoods using the Felsenstein pruning algorithm
    int stride = 4 * numGammaCats;
    for (int n=0; n<dps.size(); n++)
        {
        Node* p = dps[n];
        if (p->getIsLeaf() == false && p->getUpdate() == true)
            {
            std::vector<Node*> des = p->getDescendants();
            if (des.size() != 2)
                Msg::error("There should be only two descendants");
            Node* ndeL = des[0];
            Node* ndeR = des[1];
            int lftIdx = ndeL->getIndex();
            int rhtIdx = ndeR->getIndex();
            ConditionalLikelihoods* clInfoL = cls[ndeL->getActiveCl()][lftIdx];
            ConditionalLikelihoods* clInfoR = cls[ndeR->getActiveCl()][rhtIdx];
            ConditionalLikelihoods* clInfoP = cls[p->getActiveCl()][p->getIndex()];

            double* clL = clInfoL->getCondLike();
            double* clR = clInfoR->getCondLike();
            double* clP = clInfoP->getCondLike();
            double*** pL = tis[ndeL->getActiveTi()][lftIdx]->getTiProbs();
            double*** pR = tis[ndeR->getActiveTi()][rhtIdx]->getTiProbs();
            for (int c=0; c<numSites; c++)
                {
                for (int k=0; k<numGammaCats; k++)
                    {
#                   if 1
                    __m256d clL_0123 = _mm256_load_pd(clL);
                    __m256d clR_0123 = _mm256_load_pd(clR);
                    
                    // A
                    __m256d a_tis_l_0123 = _mm256_load_pd(pL[k][0]);
                    __m256d a_tis_r_0123 = _mm256_load_pd(pR[k][0]);
                    __m256d a_l_p0123 = _mm256_mul_pd(clL_0123, a_tis_l_0123);
                    __m256d a_r_p0123 = _mm256_mul_pd(clR_0123, a_tis_r_0123);

                    // C
                    __m256d c_tis_l_0123 = _mm256_load_pd(pL[k][1]);
                    __m256d c_tis_r_0123 = _mm256_load_pd(pR[k][1]);
                    __m256d c_l_p0123 = _mm256_mul_pd(clL_0123, c_tis_l_0123);
                    __m256d c_r_p0123 = _mm256_mul_pd(clR_0123, c_tis_r_0123);

                    // G
                    __m256d g_tis_l_0123 = _mm256_load_pd(pL[k][2]);
                    __m256d g_tis_r_0123 = _mm256_load_pd(pR[k][2]);
                    __m256d g_l_p0123 = _mm256_mul_pd(clL_0123, g_tis_l_0123);
                    __m256d g_r_p0123 = _mm256_mul_pd(clR_0123, g_tis_r_0123);

                    // T
                    __m256d t_tis_l_0123 = _mm256_load_pd(pL[k][3]);
                    __m256d t_tis_r_0123 = _mm256_load_pd(pR[k][3]);
                    __m256d t_l_p0123 = _mm256_mul_pd(clL_0123, t_tis_l_0123);
                    __m256d t_r_p0123 = _mm256_mul_pd(clR_0123, t_tis_r_0123);

                    // sum left
                    __m256d sumL_AC = _mm256_hadd_pd(a_l_p0123, c_l_p0123);
                    __m256d sumL_GT = _mm256_hadd_pd(g_l_p0123, t_l_p0123);
                    __m256d blendL = _mm256_blend_pd(sumL_AC, sumL_GT, 0b1100);
                    __m256d permL = _mm256_permute2f128_pd(sumL_AC, sumL_GT, 0x21);
                    __m256d sumL = _mm256_add_pd(permL, blendL);

                    // sum right
                    __m256d sumR_AC = _mm256_hadd_pd(a_r_p0123, c_r_p0123);
                    __m256d sumR_GT = _mm256_hadd_pd(g_r_p0123, t_r_p0123);
                    __m256d blendR = _mm256_blend_pd(sumR_AC, sumR_GT, 0b1100);
                    __m256d permR = _mm256_permute2f128_pd(sumR_AC, sumR_GT, 0x21);
                    __m256d sumR = _mm256_add_pd(permR, blendR);

                    // product
                    __m256d res = _mm256_mul_pd(sumL, sumR);
                    _mm256_store_pd(clP, res);

#                   else

                    __m128d clL_01 = _mm_load_pd(clL);
                    __m128d clL_23 = _mm_load_pd(clL+2);
                    
                    __m128d clR_01 = _mm_load_pd(clR);
                    __m128d clR_23 = _mm_load_pd(clR+2);
                    
                    /* A */
                    __m128d a_tis_l_01 = _mm_load_pd(pL[k][0]  );
                    __m128d a_tis_l_23 = _mm_load_pd(pL[k][0]+2);
                    
                    __m128d a_tis_r_01 = _mm_load_pd(pR[k][0]  );
                    __m128d a_tis_r_23 = _mm_load_pd(pR[k][0]+2);
                    
                    __m128d a_l_p01 = _mm_mul_pd(clL_01,a_tis_l_01);
                    __m128d a_l_p23 = _mm_mul_pd(clL_23,a_tis_l_23);
                    
                    __m128d a_r_p01 = _mm_mul_pd(clR_01,a_tis_r_01);
                    __m128d a_r_p23 = _mm_mul_pd(clR_23,a_tis_r_23);
                    
                    __m128d a_sum_L = _mm_hadd_pd(a_l_p01,a_l_p23);
                    __m128d a_sum_R = _mm_hadd_pd(a_r_p01,a_r_p23);
                    
                    /* C */
                    __m128d c_tis_l_01 = _mm_load_pd(pL[k][1]  );
                    __m128d c_tis_l_23 = _mm_load_pd(pL[k][1]+2);
                    
                    __m128d c_tis_r_01 = _mm_load_pd(pR[k][1]  );
                    __m128d c_tis_r_23 = _mm_load_pd(pR[k][1]+2);
                    
                    __m128d c_l_p01 = _mm_mul_pd(clL_01,c_tis_l_01);
                    __m128d c_l_p23 = _mm_mul_pd(clL_23,c_tis_l_23);
                    
                    __m128d c_r_p01 = _mm_mul_pd(clR_01,c_tis_r_01);
                    __m128d c_r_p23 = _mm_mul_pd(clR_23,c_tis_r_23);
                    
                    __m128d c_sum_L = _mm_hadd_pd(c_l_p01,c_l_p23);
                    __m128d c_sum_R = _mm_hadd_pd(c_r_p01,c_r_p23);
                    
                    /* G */
                    __m128d g_tis_l_01 = _mm_load_pd(pL[k][2]  );
                    __m128d g_tis_l_23 = _mm_load_pd(pL[k][2]+2);
                    
                    __m128d g_tis_r_01 = _mm_load_pd(pR[k][2]  );
                    __m128d g_tis_r_23 = _mm_load_pd(pR[k][2]+2);
                    
                    __m128d g_l_p01 = _mm_mul_pd(clL_01,g_tis_l_01);
                    __m128d g_l_p23 = _mm_mul_pd(clL_23,g_tis_l_23);
                    
                    __m128d g_r_p01 = _mm_mul_pd(clR_01,g_tis_r_01);
                    __m128d g_r_p23 = _mm_mul_pd(clR_23,g_tis_r_23);
                    
                    __m128d g_sum_L = _mm_hadd_pd(g_l_p01,g_l_p23);
                    __m128d g_sum_R = _mm_hadd_pd(g_r_p01,g_r_p23);
                    
                    /* T */
                    __m128d t_tis_l_01 = _mm_load_pd(pL[k][3]  );
                    __m128d t_tis_l_23 = _mm_load_pd(pL[k][3]+2);
                    
                    __m128d t_tis_r_01 = _mm_load_pd(pR[k][3]  );
                    __m128d t_tis_r_23 = _mm_load_pd(pR[k][3]+2);
                    
                    __m128d t_l_p01 = _mm_mul_pd(clL_01,t_tis_l_01);
                    __m128d t_l_p23 = _mm_mul_pd(clL_23,t_tis_l_23);
                    
                    __m128d t_r_p01 = _mm_mul_pd(clR_01,t_tis_r_01);
                    __m128d t_r_p23 = _mm_mul_pd(clR_23,t_tis_r_23);
                    
                    __m128d t_sum_L = _mm_hadd_pd(t_l_p01,t_l_p23);
                    __m128d t_sum_R = _mm_hadd_pd(t_r_p01,t_r_p23);
                    
                    // now put it all together
                    __m128d ac_sum_L = _mm_hadd_pd(a_sum_L,c_sum_L);
                    __m128d ac_sum_R = _mm_hadd_pd(a_sum_R,c_sum_R);
                    __m128d gt_sum_L = _mm_hadd_pd(g_sum_L,t_sum_L);
                    __m128d gt_sum_R = _mm_hadd_pd(g_sum_R,t_sum_R);
                    
                    __m128d ac = _mm_mul_pd(ac_sum_L,ac_sum_R);
                    __m128d gt = _mm_mul_pd(gt_sum_L,gt_sum_R);
                    
                    _mm_store_pd(clP,ac);
                    _mm_store_pd(clP+2,gt);
#                   endif
                    
                    clP += 4;
                    clL += 4;
                    clR += 4;
                    }
                }
                
            // scale
            double* scL = clInfoL->getLnScalerDp();
            double* scR = clInfoR->getLnScalerDp();
            double* scP = clInfoP->getLnScalerDp();
            double* lnS = clInfoP->getLnScaler();
            clP = clInfoP->getCondLike();
            for (int c=0; c<numSites; c++)
                {
                double maxCl = 0.0;
                for (int i=0; i<stride; i++)
                    {
                    if (clP[i] > maxCl)
                        maxCl = clP[i];
                    }
                
#               if 1
                __m256d scaleFactor = _mm256_set1_pd(1.0 / maxCl);
                for (int k=0; k<numGammaCats; k++)
                    {
                    __m256d x = _mm256_load_pd(clP);
                    _mm256_store_pd(clP, _mm256_mul_pd(x, scaleFactor) );
                    clP += 4;
                    }
                lnS[c] = log(maxCl);
                scP[c] = scL[c] + scR[c] + lnS[c];
#               else
                double scaleFactor = 1.0 / maxCl;
                for (int i=0; i<stride; i++)
                    clP[i] *= scaleFactor;
                lnS[c] = log(maxCl);
                scP[c] = scL[c] + scR[c] + lnS[c];
                clP += stride;
#               endif
                }
                
            }
        }
    
    // calculate the likelihood, using the conditional likelihoods at the root of the tree
    double gammaWeight = 1.0 / numGammaCats;
    ConditionalLikelihoods* clInfoR = cls[t->getRoot()->getActiveCl()][t->getRoot()->getIndex()];
    double* clR = clInfoR->getCondLike();
    double* scR = clInfoR->getLnScalerDp();
    double lnL = 0.0;
    for (int c=0; c<numSites; c++)
        {
#       if 1
        __m256d likeV = _mm256_set1_pd(0.0);
        __m256d gammaWeightV = _mm256_set1_pd(gammaWeight);
        __m256d piV = _mm256_load_pd(&pi[whichSpace][0]);
        for (int k=0; k<numGammaCats; k++)
            {
            __m256d clVR = _mm256_load_pd(clR);
            __m256d prod = _mm256_mul_pd( gammaWeightV, _mm256_mul_pd(clVR, piV) );
            likeV = _mm256_add_pd(likeV, prod);
            clR += 4;
            }
        __m256d s = _mm256_hadd_pd(likeV, likeV);
        double like = ((double*)&s)[0] + ((double*)&s)[2];
        lnL += (log(like) + scR[c]);
#       else
        double like = 0.0;
        for (int k=0; k<numGammaCats; k++)
            {
            for (int i=0; i<4; i++)
                like += clR[i] * pi[whichSpace][i] * gammaWeight;
            clR += 4;
            }
        lnL += (log(like) + scR[c]);
#       endif
        }
    
    return lnL;
}

# else

double Model::lnLikelihood(int whichSpace) {

#   if defined(RUN_ON_PRIOR)
    return 0.0;
#   endif

    // get the tree
    ParameterTree* t = getParameterTree(whichSpace);
    
    // get the downpass sequence
    std::vector<Node*> dps = t->getDownPassSequence();
    
    // fill in the conditional likelihoods using the Felsenstein pruning algorithm
    for (int n=0; n<dps.size(); n++)
        {
        Node* p = dps[n];
        if (p->getIsLeaf() == false && p->getUpdate() == true)
            {
            std::vector<Node*> des = p->getDescendants();
            if (des.size() != 2)
                Msg::error("There should be only two descendants");
            Node* ndeL = des[0];
            Node* ndeR = des[1];
            int lftIdx = ndeL->getIndex();
            int rhtIdx = ndeR->getIndex();
            ConditionalLikelihoods* clInfoL = cls[ ndeL->getActiveCl() ][ lftIdx        ];
            ConditionalLikelihoods* clInfoR = cls[ ndeR->getActiveCl() ][ rhtIdx        ];
            ConditionalLikelihoods* clInfoP = cls[    p->getActiveCl() ][ p->getIndex() ];

            double* clL = clInfoL->getCondLike();
            double* clR = clInfoR->getCondLike();
            double* clP = clInfoP->getCondLike();
            double* scL = clInfoL->getLnScalerDp();
            double* scR = clInfoR->getLnScalerDp();
            double* scP = clInfoP->getLnScalerDp();
            double* lnS = clInfoP->getLnScaler();
            double*** pL = tis[ndeL->getActiveTi()][lftIdx]->getTiProbs();
            double*** pR = tis[ndeR->getActiveTi()][rhtIdx]->getTiProbs();
            for (int c=0; c<numSites; c++)
                {
                double* clPStart = clP;
                double maxCl = 0.0;
                for (int k=0; k<numGammaCats; k++)
                    {
                    for (int i=0; i<4; i++)
                        {
                        double sumL = 0.0, sumR = 0.0;
                        for (int j=0; j<4; j++)
                            {
                            sumL += pL[k][i][j] * clL[j];
                            sumR += pR[k][i][j] * clR[j];
                            }
                        clP[i] = sumL * sumR;
                        if (clP[i] > maxCl)
                            maxCl = clP[i];
                        }
                    clP += 4;
                    clL += 4;
                    clR += 4;
                    }
                double scaleFactor = 1.0 / maxCl;
                for (int i=0; i<numGammaCats*4; i++)
                    clPStart[i] *= scaleFactor;
                lnS[c] = log(maxCl);
                scP[c] = scL[c] + scR[c] + lnS[c];
                }
            }
        }
    
    // calculate the likelihood, using the conditional likelihoods at the root of the tree
    double gammaWeight = 1.0 / numGammaCats;
    ConditionalLikelihoods* clInfoR = cls[t->getRoot()->getActiveCl()][t->getRoot()->getIndex()];
    double* clR = clInfoR->getCondLike();
    double* scR = clInfoR->getLnScalerDp();
    double lnL = 0.0;
    for (int c=0; c<numSites; c++)
        {
        double like = 0.0;
        for (int k=0; k<numGammaCats; k++)
            {
            for (int i=0; i<4; i++)
                like += clR[i] * pi[whichSpace][i] * gammaWeight;
            clR += 4;
            }
        lnL += (log(like) + scR[c]);
        }
    
    return lnL;
}

# endif

double Model::lnPrior(int whichSpace) {

    double lnP = 0.0;
    
    // tree prior
    ParameterTree* t = getParameterTree(whichSpace);
    lnP += t->lnPriorProb();
    
    // gamma shape prior
    ParameterAsrv* asrv = getParameterAsrv(whichSpace);
    lnP += rv->lnExponentialPdf( gammaShapePrior, asrv->getAsrv() );
    
    // substitution rates
    ParameterExchangabilityRates* sr = getParameterExchangabilityRates(whichSpace);
    lnP += rv->lnDirichletPdf( sr->getExchangabilityAlpha(), sr->getExchangabilityRates() );
    
    // base frequencies
    if (isTimeReversible == true)
        {
        ParameterBaseFrequencies* bf = getParameterBaseFrequencies(whichSpace);
        lnP += rv->lnDirichletPdf( bf->getBaseFrequencyAlpha(), bf->getBaseFrequencies() );
        }
    
    return lnP;
}

int Model::minimumNumberOfTries(void) {

    int m = 0;
    bool initializedMin = false;
    for (int i=0; i<proposals.size(); i++)
        {
        if (proposals[i]->getProbability() > 0.0000001)
            {
            int n = proposals[i]->getNumAttempts();
            if (initializedMin == false)
                {
                m = n;
                initializedMin = true;
                }
            else
                {
                if (n < m)
                    m = n;
                }
            }
        }
    return m;
}

void Model::printTransitionProbabilities(int whichSpace) {

    ParameterTree* t = getParameterTree(whichSpace);
    int nNodes = (int)t->getNumNodes();
    
    for (int i=0; i<nNodes; i++)
        {
        TransitionProbabilities* tp = tis[whichSpace][i];
        std::cout << "Node " << i << std::endl;
        tp->print();
        }
}

void Model::summarizeProposals(void) {

    std::cout << std::endl;
    std::cout << "   * MCMC proposal summary" << std::endl;
    for (int i=0; i<proposals.size(); i++)
        {
        if (proposals[i]->getNumAttempts() > 0)
            {
            double percent = proposals[i]->acceptanceFrequency() * 100.0;
            std::cout << "   *    Accepted " << std::fixed << std::setprecision(1) << percent << "% of proposed moves for " << proposals[i]->getName() << std::endl;
            }
        }
}

void Model::rejectMove(void) {

    proposals[lastProposal]->reject();

    (*updatedParameters[1]) = (*updatedParameters[0]);

    /*for (int i=0; i<numNodes; i++)
        {
        (*cls[1][i]) = (*cls[0][i]);
        (*tis[1][i]) = (*tis[0][i]);
        }*/

    bool properlyUpdated = false;
    
    ParameterBaseFrequencies* p1 = dynamic_cast<ParameterBaseFrequencies*>(updatedParameters[0]);
    if (p1 != NULL)
        {
        ParameterTree* t0 = getParameterTree(0);
        ParameterTree* t1 = getParameterTree(1);
        (*t1) = (*t0);

        eigenValues[1] = eigenValues[0];
        pi[1] = pi[0];
        dStatistic[1] = dStatistic[0];
        for (int i=0; i<64; i++)
            cc_ijk[1][i] = cc_ijk[0][i];
        properlyUpdated = true;
        }

    ParameterExchangabilityRates* p2 = dynamic_cast<ParameterExchangabilityRates*>(updatedParameters[0]);
    if (p2 != NULL)
        {
        ParameterTree* t0 = getParameterTree(0);
        ParameterTree* t1 = getParameterTree(1);
        (*t1) = (*t0);

        eigenValues[1] = eigenValues[0];
        pi[1] = pi[0];
        dStatistic[1] = dStatistic[0];
        for (int i=0; i<64; i++)
            cc_ijk[1][i] = cc_ijk[0][i];
        properlyUpdated = true;
        }

    ParameterAsrv* p3 = dynamic_cast<ParameterAsrv*>(updatedParameters[0]);
    if (p3 != NULL)
        {
        ParameterTree* t0 = getParameterTree(0);
        ParameterTree* t1 = getParameterTree(1);
        (*t1) = (*t0);

        properlyUpdated = true;
        }

    ParameterTree* p4 = dynamic_cast<ParameterTree*>(updatedParameters[0]);
    if (p4 != NULL)
        {
        // no need to copy the tree because it was already done, above
        
        properlyUpdated = true;
        }
    
    if (properlyUpdated == false)
        Msg::error("Problem squaring parameters in rejectMove");

}

void Model::resetAcceptanceInformation(void) {

    for (int i=0; i<proposals.size(); i++)
        proposals[i]->reset();
}

double Model::update(void) {

    // pick a parameter to change (and remember it)
    double u = rv->uniformRv();
    double sum = 0.0;
    int whichParm = 0;
    for (int i=0; i<proposals.size(); i++)
        {
        sum += proposals[i]->getProbability();
        if (u < sum)
            {
            whichParm = i;
            break;
            }
        }
    lastProposal = whichParm;
    // change the parameter
    double lnProposalProb = 0.0;
    if (whichParm == 0)
        lnProposalProb = updateGammaShapeParm();
    else if (whichParm == 1)
        lnProposalProb = updateSubstitutionRatesParm();
    else if (whichParm == 2)
        lnProposalProb = updateBranchLengthParm();
    else if (whichParm == 3)
        lnProposalProb = updateTreeParmSpr();
    else if (whichParm == 4)
        lnProposalProb = updateTreeParmLocal();
    else if (whichParm == 5)
        lnProposalProb = updateRootPosition();
    else if (whichParm == 6)
        lnProposalProb = updateBaseFrequencies();

    return lnProposalProb;
}

double Model::updateBaseFrequencies(void) {

    ParameterTree* t = getParameterTree(1);
    t->updateNodeFlags(true);
    t->updateBranchFlags(true);
    t->flipAllActiveCls();
    t->flipAllActiveTis();
    ParameterBaseFrequencies* bf = getParameterBaseFrequencies(1);
    updatedParameters[1] = bf;
    updatedParameters[0] = getParameterBaseFrequencies(0);
    double lnProposalProb = bf->update();
    updateRateMatrix(1);
    updateTiProbs(1);
    return lnProposalProb;
}

double Model::updateGammaShapeParm(void) {

    ParameterTree* t = getParameterTree(1);
    t->updateNodeFlags(true);
    t->updateBranchFlags(true);
    t->flipAllActiveCls();
    t->flipAllActiveTis();
    ParameterAsrv* asrv = getParameterAsrv(1);
    updatedParameters[1] = asrv;
    updatedParameters[0] = getParameterAsrv(0);
    double lnProposalProb = asrv->update();
    updateTiProbs(1);
    return lnProposalProb;
}

double Model::updateSubstitutionRatesParm(void) {

    ParameterTree* t = getParameterTree(1);
    t->updateNodeFlags(true);
    t->updateBranchFlags(true);
    t->flipAllActiveCls();
    t->flipAllActiveTis();
    ParameterExchangabilityRates* sr = getParameterExchangabilityRates(1);
    updatedParameters[1] = sr;
    updatedParameters[0] = getParameterExchangabilityRates(0);
    double lnProposalProb = sr->update();
    updateRateMatrix(1);
    updateTiProbs(1);
    
    return lnProposalProb;
}

double Model::updateBranchLengthParm(void) {

    ParameterTree* t = getParameterTree(1);
    updatedParameters[1] = t;
    updatedParameters[0] = getParameterTree(0);
    double lnProposalProb = t->updateBrlen(isTimeReversible);
    //t->print("After branch length change");
    updateTiProbs(1);
    return lnProposalProb;
}

void Model::updateRateMatrix(int whichSpace) {

    //NucleotideSquareMatrix_t Q;
    NucleotideSquareMatrix_t Q;
    
    // substitution rates
    ParameterExchangabilityRates* sr = getParameterExchangabilityRates(whichSpace);
    ParameterBaseFrequencies* bf = NULL;

    // fill in off diagonal components of rate matrix
    if (isTimeReversible == false)
        {
        // non-reversible model
        int k = 0;
        for (int i=0; i<4; i++)
            {
            for (int j=0; j<4; j++)
                {
                if (i != j)
                    Q(i,j) = (*sr)[k++];
                }
            }
        }
    else
        {
        // gtr model
        bf = getParameterBaseFrequencies(whichSpace);
        int k = 0;
        for (int i=0; i<4; i++)
            {
            for (int j=i+1; j<4; j++)
                {
                Q(i,j) = (*sr)[k] * (*bf)[j];
                Q(j,i) = (*sr)[k] * (*bf)[i];
                k++;
                }
            }
        }
    
    // fill in the diagonal elements of the rate matrix
    for (int i=0; i<4; i++)
        {
        double sum = 0.0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += Q(i,j);
            }
        Q(i,i) = -sum;
        }
    
    // calculate the stationary frequencies
    if (isTimeReversible == false)
        pi[whichSpace] = calcStationaryFreq(Q);
    else
        pi[whichSpace] = bf->getBaseFrequencies();

    // rescale the rate matrix
    double averageRate = 0.0;
    for (int i=0; i<4; i++)
        averageRate += -pi[whichSpace][i] * Q(i,i);
    double scaleFactor = 1.0 / averageRate;
    Q *= scaleFactor;
    
    // calculate D statistic
    dStatistic[whichSpace] = 0.0;
    for (int i=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            dStatistic[whichSpace] += ( pi[whichSpace][i] * Q(i,j) - pi[whichSpace][j] * Q(j,i) ) * ( pi[whichSpace][i] * Q(i,j) - pi[whichSpace][j] * Q(j,i) );
            }
        }
    
    // calculate the Eigenvalues and Eigenvectors and do some precomputation
    Eigen::EigenSolver< NucleotideSquareMatrix_t > eigSolver;
    eigSolver.compute( Q, true );
    eigenValues[whichSpace] = eigSolver.eigenvalues();
    Eigen::Matrix<std::complex<double>, 4, 4> eigenVectors = eigSolver.eigenvectors();
    Eigen::Matrix<std::complex<double>, 4, 4> inverseEigenVectors = eigenVectors.inverse();

    // calculate cc_ijk
    std::complex<double>* pc = cc_ijk[whichSpace];
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            for (int k=0; k<4; k++)
                 *(pc++) = eigenVectors(i,k) * inverseEigenVectors(k,j);

    //std::cout << pi[whichSpace] << std::endl;
    //std::cout << Q << std::endl;
}

double Model::updateRootPosition(void) {

    ParameterTree* t = getParameterTree(1);
    updatedParameters[1] = t;
    updatedParameters[0] = getParameterTree(0);
    double lnProposalProb = t->updateRoot();
    updateTiProbs(1);
    return lnProposalProb;
}

double Model::updateTreeParmSpr(void) {

    ParameterTree* t = getParameterTree(1);
    updatedParameters[1] = t;
    updatedParameters[0] = getParameterTree(0);
    double lnProposalProb = t->updateSpr();
    updateTiProbs(1);
    return lnProposalProb;
}

double Model::updateTreeParmLocal(void) {

    ParameterTree* t = getParameterTree(1);
    updatedParameters[1] = t;
    updatedParameters[0] = getParameterTree(0);
    double lnProposalProb = t->updateLocal();
    updateTiProbs(1);
    return lnProposalProb;
}

void Model::updateTiProbs(int whichSpace) {

    ParameterTree* t = getParameterTree(whichSpace);

    std::vector<Node*> dps = t->getDownPassSequence();
    for (int n=0; n<dps.size(); n++)
        {
        Node* p = dps[n];
        if (p->getAncestor() != NULL)
            {
            Branch* b = t->findBranch( p, p->getAncestor() );
            if (b->getUpdate() == true)
                {
                //TransitionProbabilities* tp = tis[whichSpace][p->getIndex()];
                TransitionProbabilities* tp = tis[p->getActiveTi()][p->getIndex()];
                tp->tiProbs( b->getLength() );
                }
            //tp->print();
            }
        }

}
