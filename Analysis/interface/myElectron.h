#ifndef __MYELECTRON_HH__
#define __MYELECTRON_HH__
#include "TROOT.h"
#include "TObject.h"
using namespace std;
#include <vector>

class myElectron : public TObject {
public:

    myElectron() {
        ;
    }

    ~myElectron() {
        ;
    }

    float pt, eta, px, py, phi, charge, E, et, pz, z, mass, Energy, mt, eta_SC, phi_SC, dxy, dz, ecalEnergy, full5x5, absiso, reliso;
    float pfIsoCharged, pfIsoNeutral, pfIsoGamma, pfIsoPU;
     int gen_index, numLostHitEleInner;
    float deltaPhiSuperClusterTrackAtVtx, deltaEtaSuperClusterTrackAtVtx, sigmaIetaIeta, sigmaEtaEta;
    float ecalIso, hcalIso, caloIso, hcalOverEcal, SIP, ooEmooP;
    bool passConversionVeto;
    float rawE_SC, preshowerE_SC;
    float MVAtrigID, MVAnontrigID;
    bool cutID_loose, cutID_medium, cutID_tight, cutID_veto, MVAID_nontrig_Loose, MVAID_nontrig_Tight;

    ClassDef(myElectron, 1)
};
#endif
