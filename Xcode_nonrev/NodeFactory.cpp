#include "Node.h"
#include "NodeFactory.h"





NodeFactory::NodeFactory(void) {

}

NodeFactory::~NodeFactory(void) {

    for (std::set<Node*>::iterator nde=allocatedNodes.begin(); nde != allocatedNodes.end(); nde++)
        delete (*nde);
}

void NodeFactory::drainNodePool(void) {

    for (std::vector<Node*>::iterator nde=nodePool.begin(); nde != nodePool.end(); nde++)
        {
        allocatedNodes.erase( *nde );
        delete (*nde);
        }
}

Node* NodeFactory::getNode(void) {

    if ( nodePool.empty() == true )
        {
        /* If the node pool is empty, we allocate a new node and return it. We
           do not need to add it to the node pool. */
        Node* nde = new Node;
        allocatedNodes.insert( nde );
        return nde;
        }
    
    // Return a node from the node pool, remembering to remove it from the pool.
    Node* nde = nodePool.back();
    nodePool.pop_back();
    return nde;
}

void NodeFactory::returnNodeToPool(Node* nde) {

    nde->clean();
    nodePool.push_back( nde );
}
