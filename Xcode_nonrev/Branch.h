#ifndef Branch_H
#define Branch_H

#define MIN_LENGTH 0.000001
#define MAX_LENGTH 10.0

class Node;


class Branch {

    public:
                    Branch(void);
                    Branch(Node* e1, Node* e2);
                    Branch(Node* e1, Node* e2, double len);
        Branch&     operator=(const Branch& a);
        bool        operator==(const Branch& a) const;
        bool        operator<(const Branch& a) const;
        void        clean(void);
        Node*       getAncestralNode(void);
        Node*       getDescendantNode(void);
        Node*       getEnd1(void) const { return end1; }
        Node*       getEnd2(void) const { return end2; }
        double      getLength(void) { return length; }
        bool        getUpdate(void) { return update; }
        bool        isTip(void);
        void        print(void);
        void        setEnds(Node* e1, Node* e2);
        void        setLength(double x);
        void        setUpdate(bool tf) { update = tf; }

    protected:
        double      length;
        Node*       end1;
        Node*       end2;
        bool        update;
};

class CompBranch {

    public:
                     bool   operator()(Branch* b1, Branch* b2) const;
};

#endif
