#ifndef __MYSV_HH__
#define __MYSV_HH__
#include "TROOT.h"
#include "TObject.h"
using namespace std;
#include <vector>

class mySV : public TObject {
public:

    mySV() {
        ;
    }

    ~mySV() {
        ;
    }

    float pt, eta, phi, mass, met;
    
    ClassDef(mySV, 1)
};
#endif
