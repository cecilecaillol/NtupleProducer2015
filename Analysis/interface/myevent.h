#ifndef __MYEVENT_HH__
#define __MYEVENT_HH__
#include "TObject.h"
using namespace std;
#include <vector>
#include <map>
#include <utility>
#include "myVertex.h"
#include "myTriggerObject.h"
#include "myGenObject.h"
#include "myMuon.h"
#include "myElectron.h"
#include "myTau.h"
#include "myJet.h"
#include "myMET.h"
#include "mySV.h"

class myevent : public TObject {
public:

    myevent() {
        ;
    }

    ~myevent() {
        ;
    }

    vector<myGenObject> RecGenParticles;
    vector<myTriggerObject> TriggerObj;
    vector<myJet> PreselectedJets;
    vector<myElectron> PreSelectedElectrons;
    vector<myMuon> PreSelectedMuons;
    vector<myTau> PreSelectedTaus;
    vector<myTau> LooseTaus;
    vector<myMET> RecPFMet;
    vector<myVertex> Vertex;

    vector<myMET> PairMet_etau;
    vector<myMET> PairMet_mutau;
    vector<myMET> PairMet_tautau;
    vector<myMET> PairMet_emu;
    vector<mySV> SV_etau;
    vector<mySV> SV_mutau;
    vector<mySV> SV_tautau;
    vector<mySV> SV_emu;
    vector<myMET> PairRecoilMet_etau;
    vector<myMET> PairRecoilMet_mutau;
    vector<myMET> PairRecoilMet_tautau;
    vector<myMET> PairRecoilMet_emu;
    vector<float> PairMet_etau_sigMatrix_00;
    vector<float> PairMet_etau_sigMatrix_10;
    vector<float> PairMet_etau_sigMatrix_01;
    vector<float> PairMet_etau_sigMatrix_11;
    vector<float> PairMet_mutau_sigMatrix_00;
    vector<float> PairMet_mutau_sigMatrix_10;
    vector<float> PairMet_mutau_sigMatrix_01;
    vector<float> PairMet_mutau_sigMatrix_11;
    vector<float> PairMet_tautau_sigMatrix_00;
    vector<float> PairMet_tautau_sigMatrix_10;
    vector<float> PairMet_tautau_sigMatrix_01;
    vector<float> PairMet_tautau_sigMatrix_11;
    vector<float> PairMet_emu_sigMatrix_00;
    vector<float> PairMet_emu_sigMatrix_10;
    vector<float> PairMet_emu_sigMatrix_01;
    vector<float> PairMet_emu_sigMatrix_11;

    map<string, int> HLT;
    map<string, int> METfilters;

    unsigned int runNumber;
    unsigned int eventNumber;
    unsigned int lumiNumber;
    float PUInfo;
    float PUInfo_true;
    int PUInfo_Bunch0;
    float Rho;
    int NUP;

private:

    ClassDef(myevent, 1)
};
#endif
