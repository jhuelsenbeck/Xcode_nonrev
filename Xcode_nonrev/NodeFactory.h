#ifndef NodeFactory_H
#define NodeFactory_H

#include <set>
#include <vector>
class Node;

class NodeFactory {

    public:
        static NodeFactory&     nodeFactoryInstance(void)
                                    {
                                    static NodeFactory singleNodeFactory;
                                    return singleNodeFactory;
                                    }
        void                    drainNodePool(void);
        Node*                   getNode(void);
        void                    returnNodeToPool(Node* nde);

    private:
                                NodeFactory(void);
                                NodeFactory(const NodeFactory&);
                                NodeFactory& operator=(const NodeFactory&);
                               ~NodeFactory(void);
        std::vector<Node*>      nodePool;
        std::set<Node*>         allocatedNodes;
};

#endif
