#include "TaxonBipartition.h"



TaxonBipartition::TaxonBipartition(void) {

    _size = 0;
}

TaxonBipartition::TaxonBipartition(int s) {

    _size = s;
    bits.resize(_size);
    for (int i=0; i<_size; i++)
        bits[i] = false;
}

TaxonBipartition::TaxonBipartition(const TaxonBipartition& b) {

    _size = b._size;
    bits.resize(_size);
    for (int i=0; i<_size; i++)
        bits[i] = b.bits[i];
}

bool TaxonBipartition::operator<(const TaxonBipartition& b) const {

    // assumes this and b are of the same size
    for (int i=0; i<_size; i++)
        {
        if (bits[i] != b.bits[i])
            return b.bits[i];
        }
    return false;
}

bool TaxonBipartition::operator==(const TaxonBipartition& b) const {

    if (_size != b._size)
        return false;
    for (int i=0; i<_size; i++)
        {
        if (bits[i] != b.bits[i])
            return false;
        }
    return true;
}

void TaxonBipartition::operator|=(const TaxonBipartition& b) {

    // assumes this and b are of the same size
    for (int i=0; i<_size; i++)
        {
        if (b.bits[i] == true)
            bits[i] = true;
        }
}

bool TaxonBipartition::operator[](int idx) {

    return bits[idx];
}

void TaxonBipartition::flip(void) {

    for (int i=0; i<_size; i++)
        {
        if (bits[i] == true)
            bits[i] = false;
        else
            bits[i] = true;
        }
}

void TaxonBipartition::resize(int s) {

    _size = s;
    bits.resize(_size);
    for (int i=0; i<_size; i++)
        bits[i] = false;
}

void TaxonBipartition::setAll(bool tf) {

    for (int i=0; i<_size; i++)
        bits[i] = tf;
}
