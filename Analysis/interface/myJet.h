#ifndef __MYJET_HH__
#define __MYJET_HH__
#include "TROOT.h"
#include "TObject.h"
using namespace std;
#include <vector>

class myJet : public TObject {
public:

    myJet() {
        ;
    }

    ~myJet() {
        ;
    }

    float pt, eta, px, py, phi, E, pz, z, mass, Energy;
    float vtxMass, vtxNtracks, vtx3DVal, vtx3DSig;
    int partonFlavour, gen_index;
    bool jetId_Loose, jetId_Medium, jetId_Tight;
    float CSV;
    bool puJetId;
    float puJetIdraw;

    ClassDef(myJet, 1)
};
#endif
