#ifndef TaxonBipartition_H
#define TaxonBipartition_H

#include <vector>

class TaxonBipartition {

    public:
                            TaxonBipartition(void);
                            TaxonBipartition(int s);
                            TaxonBipartition(const TaxonBipartition& b);
        bool                operator<(const TaxonBipartition& b) const;
        bool                operator==(const TaxonBipartition& b) const;
        void                operator|=(const TaxonBipartition& b);
        bool                operator[](int idx);
        void                flip(void);
        void                resize(int s);
        int                 size(void) { return _size; }
        std::vector<bool>&  getBipartition(void) { return bits; }
        void                set(int idx) { bits[idx] = true; }
        void                setAll(bool tf);
    
    protected:
        int                 _size;
        std::vector<bool>   bits;
};

struct BipartitionComparator {

    bool operator()(const TaxonBipartition& left, const TaxonBipartition& right) const {
    
        return (left < right);
    }
};

#endif
