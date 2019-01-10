#include "TaxonBipartition.h"



TaxonBipartition::TaxonBipartition(int s) {

    size = s;
    bits.resize(size);
    for (int i=0; i<size; i++)
        bits[i] = false;
}

bool TaxonBipartition::operator<(const TaxonBipartition& b) const {

    // assumes this and b are of the same size
    for (int i=0; i<size; i++)
        {
        if (bits[i] != b.bits[i])
            return b.bits[i];
        }
    return false;
}

void TaxonBipartition::flip(void) {

    for (int i=0; i<size; i++)
        {
        if (bits[i] == true)
            bits[i] = false;
        else
            bits[i] = true;
        }
}

void TaxonBipartition::setAllBits(bool tf) {

    for (int i=0; i<size; i++)
        bits[i] = tf;
}
