#include "NtupleProducer2015/Analysis/interface/NtupleProducer.h"

// ------------ method called once each job just before starting event loop  ------------

void
NtupleProducer::beginJob() {
    hOutputFile = new TFile(fOutputFileName.c_str(), "RECREATE");
    t = new TTree("t", "tree");
    m = new myevent;
    t->Branch("myevent", "myevent", &m, 256000, 1);

}

// ------------ method called once each job just after ending the event loop  ------------

void
NtupleProducer::endJob() {

    hOutputFile->Write();
    hOutputFile->Close();
}

NtupleProducer::NtupleProducer(const edm::ParameterSet& iConfig) {

    edm::InputTag theMuonLabel("slimmedMuons");
    edm::InputTag theElectronLabel("slimmedElectrons");
    edm::InputTag theTauLabel("slimmedTaus");
    edm::InputTag theJetLabel("slimmedJets");
    edm::InputTag theVertexLabel("offlineSlimmedPrimaryVertices");
    edm::InputTag thePULabel("addPileupInfo");
    edm::InputTag theMetLabel("slimmedMETs");
    edm::InputTag theBitsLabel("TriggerResults","","HLT");
    edm::InputTag theMETFiltersLabel("TriggerResults","","PAT");
    edm::InputTag theObjectsLabel("selectedPatTrigger");
    edm::InputTag theGenLabel("prunedGenParticles");
    //edm::InputTag thePATElectronLabel("slimmedElectrons","","PAT");
    //edm::InputTag theMVAIdLabel("mvaTrigV050nsCSA14","","addMVAid");

    muonCollToken = consumes<pat::MuonCollection>(theMuonLabel);
    electronCollToken = consumes<pat::ElectronCollection>(theElectronLabel);
    tauCollToken = consumes<pat::TauCollection>(theTauLabel);
    jetCollToken = consumes<pat::JetCollection>(theJetLabel);
    pileupToken=consumes<edm::View<PileupSummaryInfo>>(thePULabel);
    vtxToken=consumes<reco::VertexCollection>(theVertexLabel);
    metToken=consumes<pat::METCollection>(theMetLabel);
    objectsToken=consumes<pat::TriggerObjectStandAloneCollection>(theObjectsLabel);
    bitsToken=consumes<edm::TriggerResults>(theBitsLabel);
    metFiltersToken=consumes<edm::TriggerResults>(theMETFiltersLabel);
    prunedGenToken=consumes<edm::View<reco::GenParticle>>(theGenLabel);

///////////////////////////////////////////////////////////////////////////
/*    std::vector<std::string> myManualCatWeigths;
    myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/CSA14/TrigIDMVA_50ns_EB_BDT.weights.xml");
    myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/CSA14/TrigIDMVA_50ns_EE_BDT.weights.xml");
    
    vector<string> myManualCatWeigthsTrig;
    string the_path;
    for (unsigned i  = 0 ; i < myManualCatWeigths.size() ; i++){
        the_path = edm::FileInPath ( myManualCatWeigths[i] ).fullPath();
        myManualCatWeigthsTrig.push_back(the_path);
    }
    
    myMVATrig = new EGammaMvaEleEstimatorCSA14();
    myMVATrig->initialize("BDT", EGammaMvaEleEstimatorCSA14::kTrig, true, myManualCatWeigthsTrig);*/
////////////////////////////////////////////////////////////////////////////

    fOutputFileName = iConfig.getUntrackedParameter<string > ("HistOutFile");

    Include_Vertex = iConfig.getParameter<bool>("Include_Vertex");
    Include_Electron = iConfig.getParameter<bool>("Include_Electron");
    Include_Muon = iConfig.getParameter<bool>("Include_Muon");
    Include_HPSTau = iConfig.getParameter<bool>("Include_HPSTau");
    Include_GenParticles = iConfig.getParameter<bool>("Include_GenParticles");
    Include_Jet = iConfig.getParameter<bool>("Include_Jet");
    Include_MET = iConfig.getParameter<bool>("Include_MET");
    Include_HLT = iConfig.getParameter<bool>("Include_HLT");
    IsMC = iConfig.getParameter<bool>("Is_MC");
    IsEmbedded = iConfig.getParameter<bool>("Is_Embedded");
    IsMT = iConfig.getParameter<bool>("Is_MT");
    IsET = iConfig.getParameter<bool>("Is_ET");
    IsEM = iConfig.getParameter<bool>("Is_EM");
    IsTT = iConfig.getParameter<bool>("Is_TT");

}

NtupleProducer::~NtupleProducer() {

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}

void
NtupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {


    if (Include_Vertex) DoVertexAnalysis(iEvent);
    if (Include_Electron) DoElectronAnalysis(iEvent, iSetup);
    if (Include_Muon) DoMuonAnalysis(iEvent, iSetup);
    if (Include_HPSTau) DoHPSTauAnalysis(iEvent, iSetup);
    if (Include_GenParticles) DoGenParticlesAnalysis(iEvent);
    if (Include_Jet) DoJetAnalysis(iEvent);
    if (Include_MET) DoMetAnalysis(iEvent);
    if (Include_HLT) DoHLTAnalysis(iEvent);

    t->Fill();

}//analyze



//define this as a plug-in
//DEFINE_FWK_MODULE(NtupleProducer);


//DEFINE_SEAL_MODULE();
DEFINE_FWK_MODULE(NtupleProducer);
