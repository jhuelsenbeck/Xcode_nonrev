#ifndef ParameterTree_H
#define ParameterTree_H

#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "Branch.h"
#include "Parameter.h"
class Node;


class ParameterTree : public Parameter {

    public:
                                        ParameterTree(RandomVariable* rp, Model* mp, std::string nm, std::vector<std::string> tn, bool useExpPrior, double a, double b, bool tao);
                                        ParameterTree(RandomVariable* rp, Model* mp, std::string nm, std::vector<std::string> tn, std::string fn, bool useExpPrior, double a, double b, bool tao);
                                        ParameterTree(RandomVariable* rp, Model* mp, std::string nm);
                                        ParameterTree(ParameterTree& t);
                                       ~ParameterTree(void);
        Parameter&                      operator=(Parameter& t);
        ParameterTree&                  operator=(ParameterTree& t);
        std::set<Branch*>               getBranchesAround(Node* p);
        double                          getBranchTuning(void) { return branchTuning; }
        std::vector<Node*>              getDownPassSequence(void) { return downPassSequence; }
        int                             getNumBranches(void) { return (int)branches.size(); }
        size_t                          getNumNodes(void) { return nodes.size(); }
        size_t                          getNumTaxa(void);
        std::string                     getNewick(void);
        std::string                     getParmHeader(int n);
        std::string                     getParmString(int n);
        Node*                           getRoot(void) { return root; }
        double                          getTreeLength(void);
        Branch*                         findBranch(Node* e1, Node* e2);
        void                            flipAllActiveCls(void);
        void                            flipAllActiveTis(void);
        double                          lnPriorProb(void);
        void                            print(void);
        void                            print(std::string header);
        void                            printNodes(void);
        void                            setAllFlags(bool tf);
        void                            setBranchTuning(double x) { branchTuning = x; }
        double                          treeLength(void);
        double                          update(void);
        double                          updateBrlen(bool isTimeReversible);
        double                          updateLocal(void);
        double                          updateRoot(void);
        double                          updateSpr(void);
        void                            updateBranchFlags(bool tf);
        void                            updateNodeFlags(bool tf);
        void                            updateNodeFlags(std::vector<Node*> nds, bool tf);

    protected:
        Branch*                         addBranch(Node* e1, Node* e2);
        Node*                           addNode(void);
        Node*                           addNode(int idx);
        Node*                           addNode(Node* nodeToCopy);
        double                          branchLengthSum(void);
        void                            buildRandomTree(std::vector<std::string> names);
        void                            buildTreeFromNewickString(std::vector<std::string> tn, std::string fn);
        void                            clone(ParameterTree& b);
        Node*                           copyNode(std::map<Node*, Node*>& visited, Node* orig);
        Node*                           cutTree(Branch* b, Node* p, Node* pAnc, double& lnP, double lambda0);
        void                            initializeDownPassSequence(void);
        int                             getIndexForTaxon(std::string n);
        Node*                           getNodeForTaxon(std::string n);
        void                            listNodes(Node* p, Node* anc, size_t indent);
        std::vector<std::string>        parseNewickString(std::string ns);
        void                            passDown(Node* p, Node* anc);
        Branch*                         randomBranch(void);
        void                            randomContiguousBranches(std::vector<Branch*>& contiguousBranches, std::vector<Node*>& contiguousNodes);
        std::string                     readNewickTree(std::string fn);
        void                            removeBranch(Branch* b);
        void                            removeBranches(void);
        void                            removeNode(Node* p);
        void                            removeNodes(void);
        void                            writeTree(Node* p, std::stringstream& ss);

        Node*                           root;
        std::set<Node*>                 nodes;
        std::set<Branch*,CompBranch>    branches;
        std::vector<std::string>        taxonNames;
        std::vector<Node*>              downPassSequence;
        double                          lambda;
        bool                            exponentialPrior;
        bool                            treatBasalBranchAsOne;
        double                          branchTuning;
};

#endif
