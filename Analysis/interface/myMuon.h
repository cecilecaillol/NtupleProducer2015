#ifndef __MYMUON_HH__
#define __MYMUON_HH__
#include "TROOT.h"
#include "TObject.h"
using namespace std;
#include <vector>

class myMuon : public TObject {
public:

    myMuon() {
        ;
    }

    ~myMuon() {
        ;
    }

    float pt, eta, px, py, phi, charge, E, et, pz, z, mass, Energy, mt, absiso, reliso;
    float pfIsoAll, pfIsoCharged, pfIsoNeutral, pfIsoGamma, pfIsoPU;
     int gen_index;
     float dB, d0, emfraction;
    float DepositR03Ecal, DepositR03Hcal, DepositR03TrackerOfficial;
    bool GlobalMuonPromptTight;
    bool isLooseMuon, isTightMuon, isSoftMuon, isHighPtMuon;
    bool TMOneStationLoose;
    bool TM2DCompatibilityLoose;
    bool isGlobalMuon, isTrackerMuon, isStandAloneMuon, isPFMuon;
    int numberOfValidMuonHits, numberOfHits,numMatchStation;
   int normalizedChi2_innTrk, numberOfValidMuonHits_innTrk,  numberOfHits_innTrk;
   float normalizedChi2;
   int trkLayerMeasure , intrkLayerMeasure,  intrkLayerpixel;
    float dxy, dz ;
    float segmentCompatibility, validFraction, trkKink, chi2LocalPosition;

    ClassDef(myMuon, 1)
};
#endif
