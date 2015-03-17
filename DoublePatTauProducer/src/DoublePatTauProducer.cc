#include <memory>

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

typedef math::XYZTLorentzVector LorentzVector;
using namespace std;
using namespace edm;
using namespace pat;

class DoublePatTauProducer : public edm::EDProducer {
public:
  explicit DoublePatTauProducer(const edm::ParameterSet&);
  ~DoublePatTauProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run&, edm::EventSetup const&);
  virtual void endRun(edm::Run&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  edm::InputTag tauSrc1_;
  edm::InputTag tauSrc2_;
  string NAME_;

};

DoublePatTauProducer::DoublePatTauProducer(const edm::ParameterSet& iConfig):
tauSrc1_(iConfig.getParameter<edm::InputTag>("tauSrc1" )),
tauSrc2_(iConfig.getParameter<edm::InputTag>("tauSrc2" )),
NAME_(iConfig.getParameter<string>("NAME" ))
{
  produces<vector<pat::Tau>>(NAME_).setBranchAlias(NAME_);
}


DoublePatTauProducer::~DoublePatTauProducer()
{
}

void
DoublePatTauProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::vector<pat::Tau> * storedTaus = new std::vector<pat::Tau>();

  Handle<pat::TauCollection> taus1Handle;
  iEvent.getByLabel(tauSrc1_, taus1Handle);
  const TauCollection &tau1 = *(taus1Handle.product());
  pat::TauCollection::const_iterator itau1 = tau1.begin();
  pat::TauCollection::const_iterator jtau1 = tau1.end();
  for (int i=0; itau1 != jtau1; ++itau1, i++) {
    pat::Tau aTau( (*taus1Handle)[i] );
    storedTaus->push_back(aTau);
  }

  Handle<pat::TauCollection> taus2Handle;
  iEvent.getByLabel(tauSrc2_, taus2Handle);
  const TauCollection &tau2 = *(taus2Handle.product());
  pat::TauCollection::const_iterator itau2 = tau2.begin();
  pat::TauCollection::const_iterator jtau2 = tau2.end();
  for (int i=0; itau2 != jtau2; ++itau2, i++) {
    pat::Tau aTau( (*taus2Handle)[i] );
    storedTaus->push_back(aTau);
  }

  std::auto_ptr<std::vector<pat::Tau> > eptr(storedTaus);
  iEvent.put(eptr,NAME_);
}

// ------------ method called once each job just before starting event loop  ------------
void
DoublePatTauProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
DoublePatTauProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void
DoublePatTauProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
DoublePatTauProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
DoublePatTauProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
DoublePatTauProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DoublePatTauProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DoublePatTauProducer);
