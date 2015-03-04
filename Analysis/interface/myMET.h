#ifndef __MYMET_HH__
#define __MYMET_HH__
#include "TROOT.h"
#include "TObject.h"
using namespace std;
#include <vector>

class myMET : public TObject {
public:

    myMET() {
        ;
    }

    ~myMET() {
        ;
    }

    float pt, eta, px, py, phi, et;
    float genMET, pt_JetEnUp, pt_JetEnDown, pt_MuonEnUp, pt_MuonEnDown, pt_ElectronEnUp, pt_ElectronEnDown, pt_UnclusteredEnUp, pt_UnclusteredEnDown, pt_JetResUp, pt_JetResDown;
    float phi_JetEnUp, phi_JetEnDown, phi_MuonEnUp, phi_MuonEnDown, phi_ElectronEnUp, phi_ElectronEnDown, phi_UnclusteredEnUp, phi_UnclusteredEnDown, phi_JetResUp, phi_JetResDown;
    
    ClassDef(myMET, 1)
};
#endif
