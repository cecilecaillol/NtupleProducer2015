#ifndef __MYTRIGGEROBJECT_HH__
#define __MYTRIGGEROBJECT_HH__
#include "TROOT.h"
#include "TObject.h"
using namespace std;
#include <vector>

class myTriggerObject : public TObject {
public:

    myTriggerObject() {
        ;
    }

    ~myTriggerObject() {
        ;
    }

    string path;
    bool isLastFilter;
    bool isL3;
    float phi, pt, eta, px, py, pz;

    ClassDef(myTriggerObject, 1)
};
#endif
