// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// needed by ntuple tau producer
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


typedef math::XYZTLorentzVector LorentzVector;
using namespace std;
using namespace edm;
using namespace pat;

class ScaleTauProducer : public edm::EDProducer {
public:
  explicit ScaleTauProducer(const edm::ParameterSet&);
  ~ScaleTauProducer();

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

  edm::InputTag tauSrc_;
  unsigned int INDEX_;
  string NAME_;

};

ScaleTauProducer::ScaleTauProducer(const edm::ParameterSet& iConfig):
tauSrc_(iConfig.getParameter<edm::InputTag>("tauSrc" )),
NAME_(iConfig.getParameter<string>("NAME" ))
{
  produces<vector<pat::Tau>>(NAME_).setBranchAlias(NAME_);
}


ScaleTauProducer::~ScaleTauProducer()
{
}


void
ScaleTauProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // get tau collection
/*  edm::Handle<edm::View<pat::Tau> > taus;
  iEvent.getByLabel(tauSrc_,taus);

  std::vector<pat::Tau> * storedTaus = new std::vector<pat::Tau>();

  const TauCollection &tault = *(taus.product());
  pat::TauCollection::const_iterator itault = tault.begin();
  pat::TauCollection::const_iterator jtault = tault.end();
  int i=0;
  for (; itault != jtault; ++itault, i++) {
*/
  std::vector<pat::Tau> * storedTaus = new std::vector<pat::Tau>();
  Handle<pat::TauCollection> tausHandleLT;
  iEvent.getByLabel(tauSrc_, tausHandleLT);
  const TauCollection &tault = *(tausHandleLT.product());
  pat::TauCollection::const_iterator itault = tault.begin();
  pat::TauCollection::const_iterator jtault = tault.end();
  for (int i=0; itault != jtault; ++itault, i++) {
  //for (int i=0; i<taus.size(); ++i){
      pat::Tau aTau( (*tausHandleLT)[i] );

      double shiftP = 1.;
      double shiftMass = 1.;
//std::cout<<"ici "<<aTau.genJet()<<std::endl;
      if ( aTau.genJet() && deltaR(aTau.p4(), aTau.genJet()->p4()) < 0.5 && aTau.genJet()->pt() > 8. ) {
	//Olivier TES
	if((aTau.signalPFChargedHadrCands()).size()==1 && (aTau.signalPFGammaCands()).size()>0){
	  shiftP = 1.01; //New correction for Winter2013
	  shiftMass = 1.01;
	}
	else if((aTau.signalPFChargedHadrCands()).size()==1 && (aTau.signalPFGammaCands()).size()==0){
	  shiftP = 1.01; //New correction for Winter2013
	  shiftMass = 1.;
	}
	else if((aTau.signalPFChargedHadrCands()).size()==3) {
	  shiftP = 1.01; //New correction for Winter2013
	  shiftMass = 1.01;
	}
      }
      
      double pxS = aTau.px()*shiftP;
      double pyS = aTau.py()*shiftP;
      double pzS = aTau.pz()*shiftP;
      double massS = aTau.mass()*shiftMass;
      double enS = TMath::Sqrt(pxS*pxS + pyS*pyS + pzS*pzS + massS*massS);
      math::XYZTLorentzVectorD p4S( pxS, pyS, pzS, enS );
      aTau.setP4(p4S);
      storedTaus->push_back(aTau);
  }
  // add the taus to the event output
  std::auto_ptr<std::vector<pat::Tau> > eptr(storedTaus);
  iEvent.put(eptr,NAME_);
}

// ------------ method called once each job just before starting event loop  ------------
void
ScaleTauProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
ScaleTauProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void
ScaleTauProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
ScaleTauProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
ScaleTauProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
ScaleTauProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ScaleTauProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ScaleTauProducer);
