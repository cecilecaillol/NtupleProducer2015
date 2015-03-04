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
    vector<myJet> RecPFJetsAK5;
    vector<myElectron> PreSelectedElectrons;
    vector<myMuon> PreSelectedMuons;
    vector<myTau> PreSelectedHPSTaus;
    vector<myMET> RecPFMet;
    vector<myMET> RecMVAMet;
    vector<myVertex> Vertex;

    map<string, int> HLT;
    map<string, int> METfilters;

    unsigned int runNumber;
    unsigned int eventNumber;
    unsigned int lumiNumber;
    float PUInfo;
    float PUInfo_true;
    int PUInfo_Bunch0;
    float RhoCorr;
    float RhoCenNeutral;
    float RhoCenCharged;
    float RhoCenNeutralTight;
    float Rho;

private:

    ClassDef(myevent, 1)
};
#endif
