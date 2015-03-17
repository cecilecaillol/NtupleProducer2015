#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"

#include <memory>

#include <memory>
#include <string>
#include <iostream>
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Utilities/interface/typelookup.h"
#include "Math/GenVector/VectorUtil.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/BTauReco/interface/IsolatedTauTagInfo.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"//for muon namespace
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/ClusterShapeFwd.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/BTauReco/interface/TrackCountingTagInfo.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronIsoCollection.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/ElectronTkIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaEcalIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "RecoEcal/EgammaClusterProducers/interface/Multi5x5ClusterProducer.h"
#include "PhysicsTools/SelectorUtils/interface/SimpleCutBasedElectronIDSelectionFunctor.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include <CLHEP/Vector/LorentzVector.h>
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <vector>
#include <iostream>
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/GenVector/VectorUtil.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronIsoCollection.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/ElectronTkIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaEcalIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "RecoEcal/EgammaClusterProducers/interface/Multi5x5ClusterProducer.h"
#include "PhysicsTools/SelectorUtils/interface/SimpleCutBasedElectronIDSelectionFunctor.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include <CLHEP/Vector/LorentzVector.h>
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/BaseTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauTagInfo.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "RecoTauTag/TauTagTools/interface/PFTauElementsOperators.h"
#include "RecoTauTag/RecoTau/interface/TauDiscriminationProducerBase.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Type1MET/interface/JetCorrExtractorT.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "NtupleProducer2015/Analysis/interface/mySV.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace std;
using namespace edm;
using namespace pat;


class SVmassProducer : public edm::EDProducer {
public:
  explicit SVmassProducer(const edm::ParameterSet&);
  ~SVmassProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run&, edm::EventSetup const&);
  virtual void endRun(edm::Run&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  // ----------member data ---------------------------

  // ----------member data ---------------------------

  edm::InputTag lepton1Src_;
  edm::InputTag lepton2Src_;
  unsigned int lepton1Index_;
  unsigned int lepton2Index_;
  unsigned int TYPE_;
  bool MC_;
  string NAME_;

};

SVmassProducer::SVmassProducer(const edm::ParameterSet& iConfig):
lepton1Src_(iConfig.getParameter<edm::InputTag>("lepton1Src" )),
lepton2Src_(iConfig.getParameter<edm::InputTag>("lepton2Src" )),
lepton1Index_(iConfig.getParameter<unsigned int>("lepton1Index" )),
lepton2Index_(iConfig.getParameter<unsigned int>("lepton2Index" )),
TYPE_(iConfig.getParameter<unsigned int>("TYPE" )),
MC_(iConfig.getParameter<bool>("MC")),
NAME_(iConfig.getParameter<string>("NAME" ))
{
  produces<vector<float>>(NAME_).setBranchAlias(NAME_);
  //produces<mySV>(NAME_).setBranchAlias(NAME_);
}


SVmassProducer::~SVmassProducer()
{
}

void
SVmassProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  double mass=-1;
  char buffer [50];
  reco::METCovMatrix covMVAMet;

  edm::Handle <pat::METCollection>  mvamet;
  if (MC_){
     if (TYPE_==1)
         sprintf (buffer, "eMuRecoil%ix%i",lepton1Index_,lepton2Index_);
     if (TYPE_==2)
         sprintf (buffer, "eTauRecoil%ix%i",lepton1Index_,lepton2Index_);
     if (TYPE_==3)
         sprintf (buffer, "muTauRecoil%ix%i",lepton1Index_,lepton2Index_);
     if (TYPE_==4)
         sprintf (buffer, "tauTauRecoil%ix%i",lepton1Index_,lepton2Index_);
  }
  else{
     if (TYPE_==1)
         sprintf (buffer, "eMuMet%ix%i",lepton1Index_,lepton2Index_);
     if (TYPE_==2)
         sprintf (buffer, "eTauMet%ix%i",lepton1Index_,lepton2Index_);
     if (TYPE_==3)
         sprintf (buffer, "muTauMet%ix%i",lepton1Index_,lepton2Index_);
     if (TYPE_==4)
         sprintf (buffer, "tauTauMet%ix%i",lepton1Index_,lepton2Index_);
  }
  iEvent.getByLabel(buffer,"N", mvamet);
  if (mvamet.isValid()){
     const pat::MET mvaMETpf =  (*mvamet)[0];
     covMVAMet = mvaMETpf.getSignificanceMatrix();
     double measuredMETx = mvaMETpf.px();
     double measuredMETy = mvaMETpf.py();
     TMatrixD covMET(2, 2);
     covMET[0][0] = covMVAMet[0][0];
     covMET[1][0] = covMVAMet[1][0];
     covMET[0][1] = covMVAMet[0][1];
     covMET[1][1] = covMVAMet[1][1];
     
     if (TYPE_==1){ //EM
        edm::Handle<edm::View<pat::Electron> > ele1;
        iEvent.getByLabel(lepton1Src_,ele1);
        edm::Handle<edm::View<pat::Muon> > muon2;
        iEvent.getByLabel(lepton2Src_,muon2);

	if (ele1->size()>0 && muon2->size()>0){
           const pat::Electron & myEle1 = ele1->at(0);
           const pat::Muon & myMuon2 = muon2->at(0);
           if (myEle1.pt()>20 && myMuon2.pt()>20 && deltaR(myEle1.eta(),myEle1.phi(),myMuon2.eta(),myMuon2.phi())>0.5){
              std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
              measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToElecDecay, myEle1.pt(),myEle1.eta(),myEle1.phi(),myEle1.mass()));
              measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, myMuon2.pt(),myMuon2.eta(),myMuon2.phi(),myMuon2.mass()));
              unsigned verbosity = 2;
              SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMETx, measuredMETy, covMET, verbosity);
              algo.addLogM(false);
              algo.integrateMarkovChain();
              mass = algo.getMass(); 
           }
	}
     }

     if (TYPE_==2){ //ET
        edm::Handle<edm::View<pat::Electron> > ele1;
        iEvent.getByLabel(lepton1Src_,ele1);
        edm::Handle<edm::View<pat::Tau> > tau2;
        iEvent.getByLabel(lepton2Src_,tau2);
	if (ele1->size()>0 && tau2->size()>0){
	   const pat::Electron & myEle1 = ele1->at(0);
	   const pat::Tau & myTau2 = tau2->at(0);
           if (myEle1.pt()>24 && myTau2.pt()>20 && deltaR(myEle1.eta(),myEle1.phi(),myTau2.eta(),myTau2.phi())>0.5){
              std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
              measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToElecDecay, myEle1.pt(),myEle1.eta(),myEle1.phi(),myEle1.mass()));
              measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, myTau2.pt(),myTau2.eta(),myTau2.phi(),myTau2.mass()));
              unsigned verbosity = 2;
              SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMETx, measuredMETy, covMET, verbosity);
              algo.addLogM(false);
              algo.integrateVEGAS();
              mass = algo.getMass();
           }
	}
     }

     if (TYPE_==3){ //MT
        edm::Handle<edm::View<pat::Muon> > muon1;
        iEvent.getByLabel(lepton1Src_,muon1);
        edm::Handle<edm::View<pat::Tau> > tau2;
        iEvent.getByLabel(lepton2Src_,tau2);
	if (muon1->size()>0 && tau2->size()>0){
	   const pat::Muon & myMuon1 = muon1->at(0);
	   const pat::Tau & myTau2 = tau2->at(0);
           if (myMuon1.pt()>20 && myTau2.pt()>20 && deltaR(myMuon1.eta(),myMuon1.phi(),myTau2.eta(),myTau2.phi())>0.5){
              std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
              measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, myMuon1.pt(),myMuon1.eta(),myMuon1.phi(),myMuon1.mass()));
              measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, myTau2.pt(),myTau2.eta(),myTau2.phi(),myTau2.mass()));
              unsigned verbosity = 2;
              SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMETx, measuredMETy, covMET, verbosity);
              algo.addLogM(false);
              algo.integrateVEGAS();
              mass = algo.getMass();
           }
	}
     }

     if (TYPE_==4){ //TT
        edm::Handle<edm::View<pat::Tau> > tau1;
        iEvent.getByLabel(lepton1Src_,tau1);
        edm::Handle<edm::View<pat::Tau> > tau2;
        iEvent.getByLabel(lepton2Src_,tau2);
	if (tau1->size()>0 && tau2->size()>0){
	   const pat::Tau & myTau1 = tau1->at(0);
	   const pat::Tau & myTau2 = tau2->at(0);
           if (myTau1.pt()>40 && myTau2.pt()>40 && deltaR(myTau1.eta(),myTau1.phi(),myTau2.eta(),myTau2.phi())>0.5){
              std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
              measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, myTau1.pt(),myTau1.eta(),myTau1.phi(),myTau1.mass()));
              measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, myTau2.pt(),myTau2.eta(),myTau2.phi(),myTau2.mass()));
              unsigned verbosity = 2;
              SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMETx, measuredMETy, covMET, verbosity);
              algo.addLogM(false);
              algo.integrateVEGAS();
              mass = algo.getMass();	
           }
	}
     }
  }

  //mySV* res=new mySV();
  //res->mass=mass;
  std::vector<float> * storedSV = new std::vector<float>();
  storedSV->push_back(mass);

  std::auto_ptr<std::vector<float> > eptr(storedSV);
  iEvent.put(eptr,NAME_);
}

// ------------ method called once each job just before starting event loop  ------------
void
SVmassProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
SVmassProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void
SVmassProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
SVmassProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
SVmassProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
SVmassProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SVmassProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SVmassProducer);
