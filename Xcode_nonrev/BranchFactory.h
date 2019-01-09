#ifndef BranchFactory_H
#define BranchFactory_H

#include <set>
#include <vector>
class Branch;

class BranchFactory {

    public:
        static BranchFactory&   branchFactoryInstance(void)
                                    {
                                    static BranchFactory singleBranchFactory;
                                    return singleBranchFactory;
                                    }
        void                    drainBranchPool(void);
        Branch*                 getBranch(void);
        void                    returnBranchToPool(Branch* nde);

    private:
                                BranchFactory(void);
                                BranchFactory(const BranchFactory&);
                                BranchFactory& operator=(const BranchFactory&);
                               ~BranchFactory(void);
        std::vector<Branch*>    branchPool;
        std::set<Branch*>       allocatedBranches;
};

#endif
