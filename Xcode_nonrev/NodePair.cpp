#include <iostream>
#include "Node.h"
#include "NodePair.h"



NodePair::NodePair(Node* f, Node* s) {

    first  = f;
    second = s;
}

bool NodePair::operator <(const NodePair& rhs) const {

    if (first < rhs.first)
        return true;
    else if (first == rhs.first)
        {
        if (second < rhs.second)
            return true;
        }
    return false;
}

void NodePair::print(void) const {

    std::cout << first->getIndex() << " " << second->getIndex() << std::endl;
}
