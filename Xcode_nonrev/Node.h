#ifndef Node_H
#define Node_H

#include <set>
#include <string>
#include <vector>
#include "TaxonBipartition.h"
class Branch;



class Node {

    public:
                            Node(void);
        void                addNeighbor(Node* p) { neighbors.insert(p); }
        void                clean(void);
        int                 degree(void) { return (int)neighbors.size(); }
        void                flipActiveCl(void);
        void                flipActiveTi(void);
        Node*               getAncestor(void) { return ancestor; }
        TaxonBipartition&   getBipartition(void) { return bipartition; }
        std::vector<Node*>  getDescendants(void);
        Branch*             getMyBranch(void) { return myBranch; }
        std::set<Node*>&    getNeighbors(void) { return neighbors; }
        std::set<Node*>     getNeighborsCopy(void) { return neighbors; }
        std::vector<Node*>  getNeighborsAsVector(void);
        int                 getIndex(void) { return index; }
        bool                getIsLeaf(void) { return isLeaf; }
        std::string         getName(void) { return name; }
        bool                getUpdate(void) { return update; }
        int                 getActiveCl(void) { return activeCl; }
        int                 getActiveTi(void) { return activeTi; }
        bool                isDescendant(Node* p);
        size_t              numNeighbors(void) { return neighbors.size(); }
        void                print(void);
        void                removeNeighbor(Node* p) { neighbors.erase(p); }
        void                removeNeighbors(void) { neighbors.clear(); }
        void                setAncestor(Node* p) { ancestor = p; }
        void                setIndex(int x) { index = x; }
        void                setIsLeaf(bool tf) { isLeaf = tf; }
        void                setMyBranch(Branch* b) { myBranch = b; }
        void                setName(std::string s) { name = s; }
        void                setUpdate(bool tf) { update = tf; }
        void                setActiveCl(int x) { activeCl = x; }
        void                setActiveTi(int x) { activeTi = x; }

    protected:
        std::set<Node*>     neighbors;
        Node*               ancestor;
        int                 index;
        bool                isLeaf;
        std::string         name;
        Branch*             myBranch;
        bool                update;
        int                 activeCl;
        int                 activeTi;
        TaxonBipartition    bipartition;
};

#endif
