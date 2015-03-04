#ifndef __MYGENOBJECT_HH__
#define __MYGENOBJECT_HH__
#include "TROOT.h"
#include "TObject.h"
using namespace std;
#include <vector>

class myGenObject : public TObject {
public:

    myGenObject() {
        ;
    }

    ~myGenObject() {
        ;
    }

    //General



    float E, Energy, px, py, pz, pt, eta,  phi, charge, z, mass, et;
    float mod_pt, mod_eta, mod_phi, mod_charge, mod_z, mod_mass;
    float Gmod_pt, Gmod_eta, Gmod_phi, Gmod_charge, Gmod_z, Gmod_mass;

    int pdgId, status, mod_pdgId, mod_status, Gmod_pdgId, Gmod_status;
    int gen_index, decay_mode;

    ClassDef(myGenObject, 1)
};
#endif
