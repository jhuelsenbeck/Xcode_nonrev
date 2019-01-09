#ifndef NodePair_H
#define NodePair_H

#include "Node.h"

class NodePair {

    public:
                                NodePair(Node* f, Node* s);
        bool                    operator <(const NodePair& rhs) const;
        friend std::ostream&    operator<<(std::ostream& output, const NodePair& p)
                                    {
                                    output << "(" << p.first->getIndex() << ", " << p.second->getIndex() << ")";
                                    return output;
                                    }
        Node*                   getFirstNode(void) { return first; }
        Node*                   getSecondNode(void) { return second; }
        Node*                   getFirstNode(void) const { return first; }
        Node*                   getSecondNode(void) const { return second; }
        void                    print(void) const;

    protected:
                                NodePair(void) { }
        Node*                   first;
        Node*                   second;
};

#endif
