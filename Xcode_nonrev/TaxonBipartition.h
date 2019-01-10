#ifndef TaxonBipartition_H
#define TaxonBipartition_H

#include <vector>

class TaxonBipartition {

    public:
                            TaxonBipartition(void) = delete;
                            TaxonBipartition(int s);
        bool                operator<(const TaxonBipartition& b) const;
        void                flip(void);
        int                 getSize(void) { return size; }
        std::vector<bool>&  getBipartition(void) { return bits; }
        void                setBit(int idx) { bits[idx] = true; }
        void                setAllBits(bool tf);
    
    protected:
        int                 size;
        std::vector<bool>   bits;
};

struct BipartitionComparator {

    bool operator()(const TaxonBipartition& left, const TaxonBipartition& right) const {
    
        return (left < right);
    }
};

#endif
