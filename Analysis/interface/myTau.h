#ifndef __MYTAU_HH__
#define __MYTAU_HH__
#include "TROOT.h"
#include "TObject.h"
using namespace std;
#include <vector>

class myTau : public TObject {
public:

    myTau() {
        ;
    }

    ~myTau() {
        ;
    }

    float pt, eta, px, py, phi, charge, E, et, pz, z, mass, dz_Ver_match, Energy, mt, jetMass, eta_SC, phi_SC;
     int gen_index;
    float dxy, dz, zImpact;
    bool isFirstVtx;
    int decayMode;
    float leadChargedParticlePt;
    float trackRefPt;
    int numChargedParticlesSignalCone, numNeutralHadronsSignalCone, numPhotonsSignalCone, numParticlesSignalCone, signalPiZeroCandidates;
    int numChargedParticlesIsoCone, numNeutralHadronsIsoCone, numPhotonsIsoCone, numParticlesIsoCone;
    float ptSumChargedParticlesIsoCone, ptSumPhotonsIsoCone;
    bool discriminationByDecayModeFinding;
    bool discriminationByElectronMVA5VLoose;
    bool discriminationByElectronMVA5Loose;
    bool discriminationByElectronMVA5Medium;
    bool discriminationByMuonLoose3;
    bool discriminationByMuonTight3;
    bool byLooseCombinedIsolationDeltaBetaCorr3Hits;
    bool byMediumCombinedIsolationDeltaBetaCorr3Hits;
    bool byTightCombinedIsolationDeltaBetaCorr3Hits;
    bool byVLooseIsolationMVA3oldDMwLT, byLooseIsolationMVA3oldDMwLT, byMediumIsolationMVA3oldDMwLT, byTightIsolationMVA3oldDMwLT, byVTightIsolationMVA3oldDMwLT;
    float byRawCombinedIsolationDeltaBetaCorr3Hits;
    bool discriminationByDecayModeFindingNewDMs;
    bool discriminationByDecayModeFindingOldDMs;
    float IsolationChargedIsoPtSum;
    float IsolationNeutralIsoPtSum;
    float IsolationPUcorrPtSum;

    ClassDef(myTau, 1)
};
#endif
