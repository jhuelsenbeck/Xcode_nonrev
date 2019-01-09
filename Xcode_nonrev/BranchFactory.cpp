#include "Branch.h"
#include "BranchFactory.h"



BranchFactory::BranchFactory(void) {

}

BranchFactory::~BranchFactory(void) {

    for (std::set<Branch*>::iterator b=allocatedBranches.begin(); b != allocatedBranches.end(); b++)
        delete (*b);
}

void BranchFactory::drainBranchPool(void) {

    for (std::vector<Branch*>::iterator b=branchPool.begin(); b != branchPool.end(); b++)
        {
        allocatedBranches.erase( *b );
        delete (*b);
        }
}

Branch* BranchFactory::getBranch(void) {

    if ( branchPool.empty() == true )
        {
        /* If the branch pool is empty, we allocate a new branch and return it. We
           do not need to add it to the branch pool. */
        Branch* b = new Branch;
        allocatedBranches.insert( b );
        return b;
        }
    
    // Return a branch from the branch pool, remembering to remove it from the pool.
    Branch* b = branchPool.back();
    branchPool.pop_back();
    return b;
}

void BranchFactory::returnBranchToPool(Branch* b) {

    b->clean();
    branchPool.push_back( b );
}
