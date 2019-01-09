#include "Parameter.h"



Parameter::Parameter(RandomVariable* rp, Model* mp, std::string nm) {

    rv          = rp;
    modelPtr    = mp;
    name        = nm;
}

Parameter::~Parameter(void) {

}
