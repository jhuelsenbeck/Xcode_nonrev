#ifndef Alignment_H
#define Alignment_H

#include <string>
#include <vector>


class Alignment {

    public:
                                    Alignment(void) = delete;
                                    Alignment(std::string fileName);
                                   ~Alignment(void);
        void                        compress(void);
        bool                        getIsExcluded(size_t i) { return isExcluded[i]; }
        int                         getNumTaxa(void) { return numTaxa; }
        int                         getNumChar(void) { return numChar; }
        void                        getPossibleNucs (int nucCode, int* nuc);
        int                         getNucleotide(size_t i, size_t j);
        size_t                      getNumSubsets(void);
        int                         getTaxonIndex(std::string ns);
        std::vector<std::string>    getTaxonNames(void);
        size_t                      getPartitionId(size_t i) { return partitionId[i]; }
        std::string                 getTaxonName(int i);
        void                        listTaxa(void);
        void                        print(void);
        void                        uncompress(void);

    protected:
        void                        interpretString(std::string s, bool* v, int n);
        bool                        isNumber(std::string s);
        int                         nucID(char nuc);
        std::vector<std::string>    taxonNames;
        int**                       matrix;
        size_t*                     partitionId;
        bool*                       isExcluded;
        int                         numTaxa;
        int                         numChar;
        int                         numSitePatterns;
};

#endif
