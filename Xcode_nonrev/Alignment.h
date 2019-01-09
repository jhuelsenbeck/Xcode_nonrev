#ifndef Alignment_H
#define Alignment_H

#include <string>
#include <vector>


class Alignment {

    public:
                                    Alignment(void) = delete;
                                    Alignment(std::string fileName);
                                   ~Alignment(void);
        int                         getNumTaxa(void) { return numTaxa; }
        int                         getNumChar(void) { return numChar; }
        void                        getPossibleNucs (int nucCode, int* nuc);
        int                         getNucleotide(size_t i, size_t j);
        int                         getTaxonIndex(std::string ns);
        std::vector<std::string>    getTaxonNames(void);
        std::string                 getTaxonName(int i);
        void                        listTaxa(void);
        void                        print(void);

    protected:
        int                         nucID(char nuc);
        std::vector<std::string>    taxonNames;
        int**                       matrix;
        int                         numTaxa;
        int                         numChar;
};

#endif
