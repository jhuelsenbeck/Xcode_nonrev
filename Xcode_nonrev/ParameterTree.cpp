#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include "Branch.h"
#include "BranchFactory.h"
#include "Model.h"
#include "Msg.h"
#include "Node.h"
#include "NodeFactory.h"
#include "NodePair.h"
#include "ParameterTree.h"
#include "RandomVariable.h"



ParameterTree::ParameterTree(RandomVariable* rp, Model* mp, std::string nm, std::vector<std::string> tn, bool useExpPrior, double p1, double p2, bool tao) : Parameter(rp, mp, nm) {

    exponentialPrior = useExpPrior;
    branchTuning = log(4.0);
    treatBasalBranchAsOne = tao;
    if (exponentialPrior == true)
        lambda = p1;
    else
        Msg::error("Non-exponential prior not yet implemented");
    buildRandomTree(tn);
}

ParameterTree::ParameterTree(RandomVariable* rp, Model* mp, std::string nm, std::vector<std::string> tn, std::string fn, bool useExpPrior, double p1, double p2, bool tao) : Parameter(rp, mp, nm) {

    exponentialPrior = useExpPrior;
    branchTuning = log(4.0);
    treatBasalBranchAsOne = tao;
    if (exponentialPrior == true)
        lambda = p1;
    else
        Msg::error("Non-exponential prior not yet implemented");
    buildTreeFromNewickString(tn, fn);
}

ParameterTree::ParameterTree(RandomVariable* rp, Model* mp, std::string nm) : Parameter(rp, mp, nm) {

}

ParameterTree::ParameterTree(ParameterTree& t) : Parameter(t.rv, t.modelPtr, t.name) {

    clone(t);
}

ParameterTree::~ParameterTree(void) {

    removeNodes();
    removeBranches();
}

Parameter& ParameterTree::operator=(Parameter& t) {

    if (this != &t)
        {
        ParameterTree* dc = dynamic_cast<ParameterTree*>(&t);
        clone(*dc);
        }
    return *this;
}

ParameterTree& ParameterTree::operator=(ParameterTree& t) {

    if (this != &t)
        {
        clone(t);
        }
    return *this;
}

Branch* ParameterTree::addBranch(Node* e1, Node* e2) {

    BranchFactory& bf = BranchFactory::branchFactoryInstance();
    Branch* b = bf.getBranch();
    b->setEnds(e1, e2);
    branches.insert(b);
    return b;
}

Node* ParameterTree::addNode(void) {

    NodeFactory& nf = NodeFactory::nodeFactoryInstance();
    Node* p = nf.getNode();
    nodes.insert(p);
    return p;
}

Node* ParameterTree::addNode(int idx) {

    Node* p = addNode();
    p->setIndex(idx);
    return p;
}

Node* ParameterTree::addNode(Node* nodeToCopy) {

    Node* p = addNode();
    p->setIndex( nodeToCopy->getIndex() ); // copy everything but pointers (neighbors & myBranch)
    p->setIsLeaf( nodeToCopy->getIsLeaf() );
    p->setName( nodeToCopy->getName() );
    p->setActiveCl( nodeToCopy->getActiveCl() );
    p->setActiveTi( nodeToCopy->getActiveTi() );
    return p;
}

double ParameterTree::branchLengthSum(void) {

    double sum = 0.0;
    for (Branch* b : branches)
        sum += b->getLength();
    return sum;
}

void ParameterTree::buildRandomTree(std::vector<std::string> names) {

    // copy the taxon names
    taxonNames = names;
    
    // find the branch factory
    BranchFactory& bf = BranchFactory::branchFactoryInstance();

    // allocate and initialize the tips
    std::vector<Node*> tipNodes;
    for (size_t i=0; i<names.size(); i++)
        {
        Node* p = addNode( getIndexForTaxon(names[i]) );
        p->setName(names[i]);
        p->setIsLeaf(true);
        tipNodes.push_back(p);
        }
    
    // initialize the two-species tree
    int nextTip = 0, nextInt = (int)taxonNames.size();;
    root = addNode(nextInt++);
    std::vector<Node*> availableNodes;
    availableNodes.push_back(root);
    for (size_t i=0; i<2; i++)
        {
        Node* p = tipNodes[nextTip++];
        root->addNeighbor(p);
        p->addNeighbor(root);
        p->setAncestor(root);
        
        availableNodes.push_back(p);
        
        Branch* b = bf.getBranch();
        b->setEnds(p, root);
        branches.insert(b);
        }
    
    // build the remaining portion of the tree
    for (size_t i=2; i<names.size(); i++)
        {
        // pick a node at random
        Node* n1 = availableNodes[(int)(rv->uniformRv()*availableNodes.size())];
        if (n1 != root)
            {
            Node* n2 = n1->getAncestor();
            Branch* b0 = findBranch( n1, n2 );

            n1->removeNeighbor(n2);
            n2->removeNeighbor(n1);
            Node* newTip = tipNodes[nextTip++];
            Node* newInt = addNode(nextInt++);
            newTip->addNeighbor(newInt);
            newInt->addNeighbor(newTip);
            n1->addNeighbor(newInt);
            newInt->addNeighbor(n1);
            n2->addNeighbor(newInt);
            newInt->addNeighbor(n2);
            newTip->setAncestor(newInt);
            newInt->setAncestor(n2);
            n1->setAncestor(newInt);

            branches.erase(b0);
            b0->setEnds(n1, newInt);
            Branch* b1 = bf.getBranch();
            Branch* b2 = bf.getBranch();
            b1->setEnds(n2, newInt);
            b2->setEnds(newTip, newInt);
            branches.insert(b0);
            branches.insert(b1);
            branches.insert(b2);
            }
        else
            {
            Node* newTip = tipNodes[nextTip++];
            Node* newInt = addNode(nextInt++);
            newTip->addNeighbor(newInt);
            newInt->addNeighbor(newTip);
            newInt->addNeighbor(n1);
            n1->addNeighbor(newInt);
            root = newInt;
            n1->setAncestor(newInt);
            newTip->setAncestor(newInt);
 
            Branch* b1 = bf.getBranch();
            Branch* b2 = bf.getBranch();
            b1->setEnds(n1, newInt);
            b2->setEnds(newTip, newInt);
            branches.insert(b1);
            branches.insert(b2);
           }
        }
    
    // initialize the branch length proportions
    for (Branch* b : branches)
        {
        if (exponentialPrior == true)
            {
            double v = rv->exponentialRv(lambda);
            b->setLength(v);
            }
        else
            Msg::error("Compound dirichlet prior not implemented");
        }
    
    // initialize the traversal sequence for the tree
    initializeDownPassSequence();
    
    // set sizes of taxon bipartitions
    for (Node* p : nodes)
        p->getBipartition().resize((int)taxonNames.size());

    //printNodes();
    print();
}

void ParameterTree::buildTreeFromNewickString(std::vector<std::string> tn, std::string fn) {

    // copy the taxon names
    taxonNames = tn;
    
    std::string newickStr = readNewickTree(fn);

    std::vector<std::string> newickTokens;
    newickTokens = parseNewickString(newickStr);
    
    Node* p = NULL;
    int intIdx = (int)tn.size();
    bool readingBrlen = false;
    for (int i=0; i<newickTokens.size(); i++)
        {
        std::string token = newickTokens[i];
        if (token == "(")
            {
            if (p == NULL)
                {
                p = addNode(intIdx++);
                root = p;
                }
            else
                {
                Node* q = addNode(intIdx++);
                p->addNeighbor(q);
                q->addNeighbor(p);
                q->setAncestor(p);
                Branch* b = addBranch(p, q);
                if (exponentialPrior == true)
                    b->setLength(rv->exponentialRv(lambda));
                else
                    Msg::error("Compound dirichlet prior not yet implemented");
                p = q;
                p->setMyBranch(b);
                }
            p->setIsLeaf(false);
            readingBrlen = false;
            }
        else if (token == ")" || token == ",")
            {
            if (p->getAncestor() != NULL)
                p = p->getAncestor();
            else
                Msg::error("Tried to move down tree");
            readingBrlen = false;
            }
        else if (token == ":")
            {
            readingBrlen = true;
            }
        else if (token == ";")
            {
            
            }
        else
            {
            if (readingBrlen == false)
                {
                Node* q = addNode();
                p->addNeighbor(q);
                q->addNeighbor(p);
                q->setAncestor(p);
                Branch* b = addBranch(p, q);
                if (exponentialPrior == true)
                    b->setLength(rv->exponentialRv(lambda));
                else
                    Msg::error("Compound dirichlet prior not yet implemented");
                p = q;
                int idx = getIndexForTaxon(token);
                p->setIndex(idx);
                p->setName(token);
                p->setIsLeaf(true);
                p->setMyBranch(b);
                }
            else
                {
                Branch* b = findBranch(p, p->getAncestor());
                if (b != NULL)
                    {
                    double x = std::stod(token);
                    b->setLength(x);
                    }
                }
            readingBrlen = false;
            }
        }
    
    // randomly root the tree
    if (root->numNeighbors() > 2)
        {
        Branch* b = randomBranch();
        double v = b->getLength();
        Node* n1 = b->getEnd1();
        Node* n2 = b->getEnd2();
        n1->removeNeighbor(n2);
        n2->removeNeighbor(n1);
        removeBranch(b);
        
        Node* p = addNode(intIdx++);
        root = p;
        p->addNeighbor(n1);
        p->addNeighbor(n2);
        n1->addNeighbor(p);
        n2->addNeighbor(p);
        Branch* b1 = addBranch(n1, p);
        Branch* b2 = addBranch(n2, p);
        double v1 = v * rv->uniformRv();
        double v2 = v-v1;
        if (v1 < 0.0)
            v1 = -v1;
        if (v2 < 0.0)
            v2 = -v2;
        b1->setLength(v1);
        b2->setLength(v2);
        
        }

    initializeDownPassSequence();

    // set sizes of taxon bipartitions
    for (Node* p : nodes)
        p->getBipartition().resize((int)taxonNames.size());
}

void ParameterTree::clone(ParameterTree& t) {

    // copy prior info
    exponentialPrior = t.exponentialPrior;
    branchTuning = t.branchTuning;
    treatBasalBranchAsOne = t.treatBasalBranchAsOne;
    lambda = t.lambda;
    
    // copy the taxon names
    taxonNames = t.taxonNames;
    
    // remove any nodes we have in this tree
    removeNodes();
    
    // make a deep copy of the nodes and neighbors
    std::map<Node*, Node*> visited;
    root = copyNode(visited, t.root);
    
    // copy the downpass sequence
    downPassSequence.clear();
    for (size_t i=0; i<t.downPassSequence.size(); i++)
        {
        std::map<Node*, Node*>::iterator it = visited.find(t.downPassSequence[i]) ;
        if ( it != visited.end() )
            downPassSequence.push_back( it->second );
        else
            Msg::error("Couldn't find node in down pass sequence");
        }
    
    // copy the branches
    removeBranches();
    for (int n=0; n<t.downPassSequence.size(); n++)
        {
        Node* pOrig = t.downPassSequence[n];
        Node* pOrigAnc = pOrig->getAncestor();
        if (pOrigAnc != NULL)
            {
            Node* p = NULL;
            std::map<Node*, Node*>::iterator it = visited.find( pOrig ) ;
            if ( it != visited.end() )
                p = it->second;
            Node* pAnc = p->getAncestor();
            if (p == NULL || pAnc == NULL)
                Msg::error("p or pAnc during copy");
            Branch* bOrig = t.findBranch(pOrig, pOrigAnc);
            Branch* b = addBranch(p, pAnc);
            b->setLength( bOrig->getLength() );
            p->setMyBranch(b);
            }
        }

    // set sizes of taxon bipartitions
    for (Node* p : nodes)
        p->getBipartition().resize((int)taxonNames.size());
}

Node* ParameterTree::copyNode(std::map<Node*, Node*>& visited, Node* orig) {

    if (orig == NULL)
        return NULL;
    
    std::map<Node*, Node*>::iterator it = visited.find(orig) ;
    if ( it != visited.end() )
        {
        return it->second;
        }

    Node* copy = addNode(orig);
    visited.insert( std::make_pair(orig, copy) );

    std::set<Node*> oNeighbors = orig->getNeighbors();
    for (std::set<Node*>::iterator it = oNeighbors.begin(); it != oNeighbors.end(); it++)
        {
        copy->addNeighbor( copyNode(visited, *it) );
        }
    copy->setAncestor( copyNode(visited, orig->getAncestor()) );
    return copy;
}

Node* ParameterTree::cutTree(Branch* b, Node* p, Node* pAnc, double& lnP, double lambda0) {

    std::set<Branch*> affectedBranches = getBranchesAround(pAnc);
    affectedBranches.erase(b);
    if (affectedBranches.size() == 1)
        {
        std::set<Node*> nbs = pAnc->getNeighborsCopy();
        nbs.erase(p);
        Node* n = *nbs.begin();
        if (n == NULL || nbs.size() > 1)
            Msg::error("Couldn't find node to reconnect in SPR move");
        
        n->removeNeighbor(pAnc);
        n->setAncestor(NULL);
        root = n;
        pAnc->removeNeighbor(n);
        pAnc->setAncestor(NULL);

        std::set<Branch*>::iterator it = affectedBranches.begin();
        double vOld = (*it)->getLength();
        removeBranch(*it);
        
        lnP += -log(lambda0) - lambda0 * vOld;
        
        return pAnc;
        }
    else if (affectedBranches.size() == 2)
        {
        Node* pAncAnc = pAnc->getAncestor();
        if (pAncAnc == NULL)
            Msg::error("Couldn't find node to reconnect in SPR move");
        std::set<Node*> nbs = pAnc->getNeighborsCopy();
        nbs.erase(p);
        nbs.erase(pAncAnc);
        Node* n = *nbs.begin();
        if (n == NULL)
            Msg::error("Couldn't find node to reconnect in SPR move");
            
        Branch* b1 = findBranch(n, pAnc);
        Branch* b2 = findBranch(pAnc, pAncAnc);
        if (b1 == NULL || b2 == NULL)
            Msg::error("Couldn't find branches to merge in SPR move");
        
        n->removeNeighbor(pAnc);
        n->addNeighbor(pAncAnc);
        n->setAncestor(pAncAnc);
        pAncAnc->removeNeighbor(pAnc);
        pAncAnc->addNeighbor(n);
        pAnc->removeNeighbor(n);
        pAnc->removeNeighbor(pAncAnc);
        pAnc->setAncestor(NULL);
        
        double sumV = b1->getLength() + b2->getLength();
        removeBranch(b1);
        removeBranch(b2);
        Branch* nb = addBranch(n, pAncAnc);
        nb->setLength(sumV);
        n->setMyBranch(nb);
        
        lnP += -log(sumV);

        return pAnc;
        }
    else
        Msg::error("Expecting one or two branches in SPR move");

    return NULL;
}

Branch* ParameterTree::findBranch(Node* e1, Node* e2) {

    if (e1 > e2)
        {
        Node* temp = e1;
        e1 = e2;
        e2 = temp;
        }

    for (Branch* b : branches)
        {
        if ( b->getEnd1() == e1 )
            {
            if ( b->getEnd2() == e2 )
                {
                return b;
                }
            }
        }
    
    return NULL;
}

void ParameterTree::flipAllActiveCls(void) {

    for (Node* n : nodes)
        n->flipActiveCl();
}

void ParameterTree::flipAllActiveTis(void) {

    for (Node* n : nodes)
        n->flipActiveTi();
}

std::set<Branch*> ParameterTree::getBranchesAround(Node* p) {

    std::set<Branch*> br;
    for (Branch* b : branches)
        {
        if (b->getEnd1() == p || b->getEnd2() == p)
            br.insert(b);
        }
    return br;
}

int ParameterTree::getIndexForTaxon(std::string n) {

    for (int i=0; i<taxonNames.size(); i++)
        {
        if (n == taxonNames[i])
            return i;
        }
    Msg::error("Cannot find taxon with name \"" + n + "\"");
    return 0;
}

std::string ParameterTree::getNewick(void) {

    std::stringstream ss;
    if (root->getIsLeaf() == true)
        {
        Node* oldRoot = root;
        std::vector<Node*> nbs = root->getNeighborsAsVector();
        if (nbs.size() > 1)
            Msg::error("Expecting only a single neighbor at the root of the tree");
        Node* newRoot = nbs[0];
        root->setAncestor(newRoot);
        oldRoot->setAncestor(newRoot);
        newRoot->setAncestor(NULL);
        root = newRoot;

        writeTree(root, ss);

        newRoot->setAncestor(oldRoot);
        oldRoot->setAncestor(NULL);
        root = oldRoot;
        }
    else
        {
        writeTree(root, ss);
        }
    std::string newick = ss.str();
    return newick;
}

Node* ParameterTree::getNodeForTaxon(std::string n) {

    for (std::set<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++)
        {
        if ((*it)->getName() == n)
            return (*it);
        }
    return NULL;
}

size_t ParameterTree::getNumTaxa(void) {

    size_t x = 0;
    for (std::set<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++)
        {
        if ( (*it)->getIsLeaf() == true )
            x++;
        }
    return x;
}

double ParameterTree::getTreeLength(void) {

    double sum = 0.0;
    for (Branch* b : branches)
        sum += b->getLength();
    return sum;
}

void ParameterTree::initializeDownPassSequence(void) {

    downPassSequence.clear();
    
    for (Node* n : nodes)
        n->setAncestor(NULL);
    passDown(root, root);
    if (root != NULL)
        root->setAncestor(NULL);
}

void ParameterTree::listNodes(Node* p, Node* anc, size_t indent) {

    if (p != NULL)
        {
        std::set<Node*> neighbors = p->getNeighbors();
        
        for (int i=0; i<indent; i++)
            std::cout << " ";
        std::cout << p->getIndex() << " ( ";
        for (Node* n : neighbors)
            {
            if (n == p->getAncestor())
                std::cout << "a.";
            std::cout << n->getIndex() << " ";
            }
        std::cout << ") ";
        if (p->getMyBranch() != NULL)
            std::cout << std::fixed << std::setprecision(6) << p->getMyBranch()->getLength();
        else
            std::cout << "NULL";
        if (p->getIsLeaf() == true)
            std::cout << " (" << p->getName() << ")";
            
        
        /* Branch* b = p->getMyBranch();
        std::cout << "     [" << p << "  ";
        if (b != NULL)
            {
            std::cout << b->getEnd1() << "  " << b->getEnd2() << " : " << b->getEnd1()->getIndex() << " - " << b->getEnd2()->getIndex() << " : " << b->getLength();
            }
        std::cout << "]";*/
        
        std::cout << "  [" << p->getUpdate() << " " << p->getActiveCl() << " ";
        if (p->getAncestor() != NULL)
            std::cout << findBranch(p, p->getAncestor())->getUpdate() << " " << p->getActiveTi();
        std::cout << "]";

        std::cout << std::endl;

        for (std::set<Node*>::iterator it = neighbors.begin(); it != neighbors.end(); it++)
            {
            if ( (*it) != anc )
                listNodes( (*it), p, indent+3 );
            }
        }
}

double ParameterTree::lnPriorProb(void) {

    if (exponentialPrior == true)
        {
        double nBranches = (double)getNumBranches();
        double T = treeLength();
        double lnP = nBranches * log(lambda) - lambda * T;
        return lnP;
        }
    else
        {
        Msg::error("Compound dirichlet prior not yet implemented");
        return 0.0;
        }
}

std::vector<std::string> ParameterTree::parseNewickString(std::string ns) {

    std::vector<std::string> tks;
    for (int i=0; i<ns.size(); i++)
        {
        char c = ns[i];
        if (c == '(' || c == ')' || c == ',' || c == ':' || c == ';')
            {
            std::string tempStr;
            tempStr = c;
            tks.push_back(tempStr);
            }
        else
            {
            int j = i;
            std::string tempStr = "";
            while ( !(c == '(' || c == ')' || c == ',' || c == ':' || c == ';') )
                {
                tempStr += c;
                j++;
                c = ns[j];
                }
            i = j-1;
            tks.push_back(tempStr);
            }
        }
#   if 0
    std::cout << "The Newick string, broken into its parts:" << std::endl;
    for (int i=0; i<tks.size(); i++)
        std::cout << "   tks[" << i << "] = \"" << tks[i] << "\"" << std::endl;
#   endif
    
    return tks;
}

void ParameterTree::passDown(Node* p, Node* anc) {

    if (p != NULL)
        {
        std::set<Node*> neighbors = p->getNeighbors();
        for (std::set<Node*>::iterator it = neighbors.begin(); it != neighbors.end(); it++)
            {
            if ( (*it) != anc )
                passDown( (*it), p );
            }
        p->setMyBranch( findBranch(p, anc) );
        p->setAncestor(anc);
        downPassSequence.push_back(p);
        }
}

void ParameterTree::print(void) {

    listNodes(root, root, 3);
}

void ParameterTree::print(std::string header) {

    std::cout << header << std::endl;
    print();
}

void ParameterTree::printNodes(void) {

    for (std::set<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++)
        (*it)->print();
    std::cout << "Down Pass Sequence: ";
    for (size_t i=0; i<downPassSequence.size(); i++)
        std::cout << downPassSequence[i]->getIndex() << " ";
    std::cout << std::endl;
}

Branch* ParameterTree::randomBranch(void) {

    size_t whichBranch = (int)(rv->uniformRv()*branches.size()), k = 0;
    for (Branch* b : branches)
        {
        if (k == whichBranch)
            return b;
        k++;
        }
    return NULL;
}

void ParameterTree::randomContiguousBranches(std::vector<Branch*>& contiguousBranches, std::vector<Node*>& contiguousNodes) {

    // pick a random internal branch
    Branch* b1 = NULL;
    do
        {
        b1 = randomBranch();
        } while (b1->isTip() == true);
    
    Node* p2 = b1->getEnd1();
    Node* p3 = b1->getEnd2();
    
    // pick branch from p2
    std::set<Branch*> s2 = getBranchesAround(p2);
    s2.erase(b1);
    int whichBranch = (int)(rv->uniformRv() * s2.size()), k = 0;
    Branch* b2 = NULL;
    for (std::set<Branch*>::iterator it = s2.begin(); it != s2.end(); it++)
        {
        if (k == whichBranch)
            {
            b2 = (*it);
            break;
            }
        k++;
        }
    Node* p1 = b2->getEnd1();
    if (p1 == p2)
        p1 = b2->getEnd2();
    
    // pick branch from p3
    std::set<Branch*> s3 = getBranchesAround(p3);
    s3.erase(b1);
    whichBranch = (int)(rv->uniformRv() * s3.size());
    k = 0;
    Branch* b3 = NULL;
    for (std::set<Branch*>::iterator it = s3.begin(); it != s3.end(); it++)
        {
        if (k == whichBranch)
            {
            b3 = (*it);
            break;
            }
        k++;
        }
    Node* p4 = b3->getEnd1();
    if (p4 == p3)
        p4 = b3->getEnd2();
    
    // add the nodes and branches to the vectors in a sensible order
    // p1 .. b2 .. p2 .. b1 .. p3 .. b3 .. p4
    contiguousNodes.push_back(p1);
    contiguousNodes.push_back(p2);
    contiguousNodes.push_back(p3);
    contiguousNodes.push_back(p4);
    contiguousBranches.push_back(b2);
    contiguousBranches.push_back(b1);
    contiguousBranches.push_back(b3);
}

std::string ParameterTree::readNewickTree(std::string fn) {

    std::ifstream treeStream(fn.c_str());
    if (!treeStream)
        Msg::error("Cannot open tree file \"" + fn + "\"");

    std::string ns = "";
    
    if (getline(treeStream, ns).good() == false)
        {
        std::cout << ns << std::endl;
        Msg::error("Failed to read tree file \"" + fn + "\"");
        }
    
    treeStream.close();
    
    return ns;
}

void ParameterTree::removeBranch(Branch* b) {

    branches.erase(b);
    BranchFactory& bf = BranchFactory::branchFactoryInstance();
    bf.returnBranchToPool(b);
}

void ParameterTree::removeBranches(void) {

    BranchFactory& bf = BranchFactory::branchFactoryInstance();
    for (Branch* b : branches)
        bf.returnBranchToPool(b);
    branches.clear();
}

void ParameterTree::removeNode(Node* p) {

    NodeFactory& nf = NodeFactory::nodeFactoryInstance();
    nf.returnNodeToPool(p);
    nodes.erase(p);
}

void ParameterTree::removeNodes(void) {

    NodeFactory& nf = NodeFactory::nodeFactoryInstance();
    for (Node* p : nodes)
        nf.returnNodeToPool(p);
    nodes.clear();
}

double ParameterTree::treeLength(void) {

    double sum = 0.0;
    for (Branch* b : branches)
        sum += b->getLength();
    return sum;
}

double ParameterTree::update(void) {

    return 0.0;
}

double ParameterTree::updateBrlen(bool isTimeReversible) {

    updateNodeFlags(false);
    updateBranchFlags(false);
    
    Branch* b = randomBranch();
    //std::cout << "updating branch length with ends " << b->getEnd1()->getIndex() << " " << b->getEnd2()->getIndex() << std::endl;
    
    bool treatRootBranchesAsOne = false;
    if (isTimeReversible == true)
        {
        if (b->getEnd1() == root || b->getEnd2() == root)
            treatRootBranchesAsOne = true;
        }
    
    double lambda0 = branchTuning;
    double newBrlen, curBrlen;
    if (treatRootBranchesAsOne == false)
        {
        curBrlen = b->getLength();
        b->setUpdate(true);
        b->getDescendantNode()->flipActiveTi();
        newBrlen = curBrlen * exp( lambda0*(rv->uniformRv()-0.5) );
        b->setLength( newBrlen );
        Node* n = b->getAncestralNode();
        while (n != NULL)
            {
            n->setUpdate(true);
            n->flipActiveCl();
            n = n->getAncestor();
            }
        if (newBrlen > MAX_LENGTH || newBrlen < MIN_LENGTH)
            return -1000000000000.0;
        }
    else
        {
        std::vector<Node*> ndes = root->getNeighborsAsVector();
        if (ndes.size() != 2)
            Msg::error("Expecting two branches at root");
        Branch* b1 = findBranch(root, ndes[0]);
        Branch* b2 = findBranch(root, ndes[1]);
        b1->setUpdate(true);
        b2->setUpdate(true);
        b1->getDescendantNode()->flipActiveTi();
        b2->getDescendantNode()->flipActiveTi();
        curBrlen = b1->getLength() + b2->getLength();
        newBrlen = curBrlen * exp( lambda0*(rv->uniformRv()-0.5) );
        b1->setLength(newBrlen * 0.5);
        b2->setLength(newBrlen * 0.5);
        root->setUpdate(true);
        root->flipActiveCl();
        if (newBrlen > MAX_LENGTH || newBrlen < MIN_LENGTH)
            return -1000000000000.0;
        }
    
    return log(newBrlen / curBrlen);
}

double ParameterTree::updateLocal(void) {

    updateNodeFlags(false);

    // get the branches and nodes for the rearrangement
    std::vector<Node*> contiguousNodes;
    std::vector<Branch*> contiguousBranches;
    randomContiguousBranches(contiguousBranches, contiguousNodes);
    
    // expand or contract the backbone branches
    double tuning = log(4.0);
    double randomFactor = exp( tuning * (rv->uniformRv()-0.5) );
    double oldM = 0.0, newM = 0.0;
    for (int i=0; i<contiguousBranches.size(); i++)
        {
        double vOld = contiguousBranches[i]->getLength();
        oldM += vOld;
        double vNew = vOld * randomFactor;
        newM += vNew;
        contiguousBranches[i]->setLength( vNew );
        }
    
    // name nodes and branches
    Node* p1 = contiguousNodes[0];
    Node* p2 = contiguousNodes[1];
    Node* p3 = contiguousNodes[2];
    Node* p4 = contiguousNodes[3];
    Branch* b1 = contiguousBranches[0];
    Branch* b2 = contiguousBranches[1];
    Branch* b3 = contiguousBranches[2];
    double v1 = b1->getLength();
    double v2 = b2->getLength();
    double v3 = b3->getLength();

    double newV1 = 0.0, newV2 = 0.0, newV3 = 0.0;
    // choose the node to slide
    if (rv->uniformRv() < 0.5)
        {
        // move p2
        if (rv->uniformRv() < (v1+v2)/newM)
            {
            // no topology change
            double u = rv->uniformRv() * (v1+v2);
            b1->setLength(u);
            b2->setLength(v1+v2-u);
            newV1 = u;
            newV2 = v1+v2-u;
            }
        else
            {
            // topology change
            removeBranch(b1);
            removeBranch(b2);
            removeBranch(b3);
            p1->removeNeighbor(p2);
            p2->removeNeighbor(p1);
            p3->removeNeighbor(p4);
            p4->removeNeighbor(p3);
            p1->addNeighbor(p3);
            p2->addNeighbor(p4);
            p3->addNeighbor(p1);
            p4->addNeighbor(p2);
            b1 = addBranch(p1, p3);
            b2 = addBranch(p3, p2);
            b3 = addBranch(p2, p4);
            b1->setLength(v1+v2);
            double u = rv->uniformRv() * v3;
            b2->setLength(u);
            b3->setLength(v3-u);
            newV1 = v1+v2;
            newV2 = u;
            newV3 = v3-u;
            }
        }
    else
        {
        // move p3
        if (rv->uniformRv() < (v2+v3)/newM)
            {
            // no topology change
            double u = rv->uniformRv() * (v2+v3);
            b2->setLength(u);
            b3->setLength(v2+v3-u);
            newV2 = u;
            newV3 = v2+v3-u;
            }
        else
            {
            // topology change
            removeBranch(b1);
            removeBranch(b2);
            removeBranch(b3);
            p1->removeNeighbor(p2);
            p2->removeNeighbor(p1);
            p3->removeNeighbor(p4);
            p4->removeNeighbor(p3);
            p1->addNeighbor(p3);
            p2->addNeighbor(p4);
            p3->addNeighbor(p1);
            p4->addNeighbor(p2);
            b1 = addBranch(p1, p3);
            b2 = addBranch(p3, p2);
            b3 = addBranch(p2, p4);
            b3->setLength(v2+v3);
            double u = rv->uniformRv() * v1;
            b1->setLength(u);
            b2->setLength(v1-u);
            newV1 = u;
            newV2 = v1-u;
            newV3 = v2+v3;
            }
        }
    
    initializeDownPassSequence();
    
    for (Node* n : contiguousNodes)
        {
        while (n != NULL)
            {
            n->setUpdate(true);
            n = n->getAncestor();
            }
        }
    
    if (v1 > MAX_LENGTH || v2 > MAX_LENGTH || v3 > MAX_LENGTH)
        return -1000000000.0;
    if (v1 < MIN_LENGTH || v2 < MIN_LENGTH || v3 < MIN_LENGTH)
        return -1000000000.0;

    return 3.0 * log( newM / oldM );
}

double ParameterTree::updateRoot(void) {

    // remove the current root node
    std::vector<Node*> nbs = root->getNeighborsAsVector();
    if (nbs.size() != 2)
        Msg::error("Expecting two descendants of root node");
    Branch* b1 = findBranch(nbs[0], root);
    Branch* b2 = findBranch(nbs[1], root);
    if (b1 == NULL || b2 == NULL)
        Msg::error("Couldn't find branches at root");
    double v1 = b1->getLength();
    double v2 = b2->getLength();
    removeBranch(b1);
    removeBranch(b2);
    nbs[0]->removeNeighbor(root);
    root->removeNeighbor(nbs[0]);
    nbs[1]->removeNeighbor(root);
    root->removeNeighbor(nbs[1]);
    nbs[0]->addNeighbor(nbs[1]);
    nbs[1]->addNeighbor(nbs[0]);
    Branch* b = addBranch(nbs[0], nbs[1]);
    b->setLength(v1 + v2);
    double lnReverseProb = log(1.0 / (v1+v2));
    
    // choose new branch at random for new root
    Branch* rb = randomBranch();
    Node* e1 = rb->getEnd1();
    Node* e2 = rb->getEnd2();
    double v = rb->getLength();
    v1 = v * rv->uniformRv();
    v2 = v - v1;
    removeBranch(rb);
    e1->removeNeighbor(e2);
    e2->removeNeighbor(e1);
    e1->addNeighbor(root);
    root->addNeighbor(e1);
    e2->addNeighbor(root);
    root->addNeighbor(e2);
    b = addBranch(e1, root);
    b->setLength(v1);
    b = addBranch(e2, root);
    b->setLength(v2);
    double lnForwardProb = log(1.0 / v);

    initializeDownPassSequence();
    
    // TEMP: This could be better
    updateNodeFlags(true);
    updateBranchFlags(true);
    flipAllActiveCls();
    flipAllActiveTis();
    
    return lnReverseProb - lnForwardProb;
}

double ParameterTree::updateSpr(void) {

    // parameter of exponential for length of branch to root
    double lambda0 = lambda;

    double lnP = 0.0;

    // pick a branch at random
    Branch* b = randomBranch();
    
    // get the ancestor and descendant nodes of the branch
    Node* p = b->getDescendantNode();
    Node* pAnc = b->getAncestralNode();
    //std::cout << "Cutting at branch " << p->getIndex() << "-" << pAnc->getIndex() << std::endl;
    
    if (p->getAncestor() != pAnc || p == NULL || pAnc == NULL)
        Msg::error("Problem in SPR move (1)");
 
    // remove subtree defined by p and pAnc
    Node* subRoot = cutTree(b, p, pAnc, lnP, lambda0);
    initializeDownPassSequence();
    
    //std::cout << "subtree(root)" << std::endl;
    //listNodes(root, root, 3);
    //std::cout << "subtree(subRoot)" << std::endl;
    //listNodes(subRoot, subRoot, 3);

    // pick a new branch at random
    p = downPassSequence[(int)(rv->uniformRv()*downPassSequence.size())];
    pAnc = p->getAncestor();
    
    // reconnect
    if (pAnc == NULL)
        {
        subRoot->addNeighbor(p);
        subRoot->setAncestor(NULL);
        p->addNeighbor(subRoot);
        p->setAncestor(subRoot);
        root = subRoot;

        Branch* nb = addBranch(p, subRoot);
        p->setMyBranch(nb);
        double vNew = rv->exponentialRv(lambda0);
        nb->setLength(vNew);

        lnP -= log(lambda0) - lambda0 * vNew;
        }
    else
        {
        b = findBranch(p, pAnc);
        double oldLength = b->getLength();
        removeBranch(b);
        
        if (p == NULL || pAnc == NULL || b == NULL)
            Msg::error("Problem in SPR move (2)");
        pAnc->removeNeighbor(p);
        pAnc->addNeighbor(subRoot);
        subRoot->addNeighbor(p);
        subRoot->addNeighbor(pAnc);
        subRoot->setAncestor(pAnc);
        p->removeNeighbor(pAnc);
        p->addNeighbor(subRoot);
        p->setAncestor(subRoot);
        
        Branch* b1 = addBranch(p, subRoot);
        Branch* b2 = addBranch(subRoot, pAnc);
        p->setMyBranch(b1);
        subRoot->setMyBranch(b2);
        double v1 = rv->uniformRv() * oldLength;
        double v2 = oldLength - v1;
        b1->setLength(v1);
        b2->setLength(v2);

        lnP -= -log(oldLength);
        }
    
    // get the new downpass sequence
    initializeDownPassSequence();

    // TEMP: This could be better
    updateNodeFlags(true);
    updateBranchFlags(true);

    // check branch lengths
    for (std::set<Branch*>::iterator it = branches.begin(); it != branches.end(); it++)
        {
        if ( (*it)->getLength() > MAX_LENGTH || (*it)->getLength() < MIN_LENGTH)
            return -1000000000.0;
        }

    return lnP;
}

void ParameterTree::updateBranchFlags(bool tf) {

    for (Branch* b : branches)
        b->setUpdate(tf);
}

void ParameterTree::updateNodeFlags(bool tf) {

    for (Node* n : nodes)
        n->setUpdate(tf);
}

void ParameterTree::updateNodeFlags(std::vector<Node*> nds, bool tf) {

    for (Node* n : nds)
        n->setUpdate(tf);
}

void ParameterTree::writeTree(Node* p, std::stringstream& ss) {

    if (p != NULL)
        {
        Branch* b = findBranch(p, p->getAncestor());
        if (p->getIsLeaf() == true)
            {
            ss << p->getIndex()+1;
            ss << ":" << b->getLength();
            }
        else
            {
            ss << "(";
            }
        std::vector<Node*> myDescendants = p->getDescendants();
        for (int i=0; i<(int)myDescendants.size(); i++)
            {
            writeTree(myDescendants[i], ss);
            if ( (i + 1) != (int)myDescendants.size() )
                ss << ",";
            }
        if (p->getIsLeaf() == false)
            {
            ss << ")";
            if (b != NULL)
                ss << ":" << b->getLength();
            }
        }
}

std::string ParameterTree::getParmString(int n) {

    return "";
}

std::string ParameterTree::getParmHeader(int n) {

    return "";
}


