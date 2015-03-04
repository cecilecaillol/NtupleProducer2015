#include <FWCore/Framework/interface/Event.h>

#include "NtupleProducer2015/Analysis/interface/NtupleProducer.h"

void NtupleProducer::DoMetAnalysis(const edm::Event& iEvent) {
    (m->RecPFMet).clear();
    edm::Handle<pat::METCollection> mets;
    iEvent.getByToken(metToken, mets);
    const pat::MET &met = mets->front();
    myMET mymet;
    mymet.pt=met.pt();
    mymet.phi=met.phi();
    mymet.et=met.sumEt();
    mymet.pt_JetEnUp=met.shiftedPt(pat::MET::JetEnUp);
    mymet.pt_JetEnDown=met.shiftedPt(pat::MET::JetEnDown);
    mymet.phi_JetEnUp=met.shiftedPhi(pat::MET::JetEnUp);
    mymet.phi_JetEnDown=met.shiftedPhi(pat::MET::JetEnDown);
    mymet.pt_ElectronEnUp=met.shiftedPt(pat::MET::ElectronEnUp);
    mymet.pt_ElectronEnDown=met.shiftedPt(pat::MET::ElectronEnDown);
    mymet.phi_ElectronEnUp=met.shiftedPhi(pat::MET::ElectronEnUp);
    mymet.phi_ElectronEnDown=met.shiftedPhi(pat::MET::ElectronEnDown);
    mymet.pt_MuonEnUp=met.shiftedPt(pat::MET::MuonEnUp);
    mymet.pt_MuonEnDown=met.shiftedPt(pat::MET::MuonEnDown);
    mymet.phi_MuonEnUp=met.shiftedPhi(pat::MET::MuonEnUp);
    mymet.phi_MuonEnDown=met.shiftedPhi(pat::MET::MuonEnDown);
    mymet.pt_JetResUp=met.shiftedPt(pat::MET::JetResUp);
    mymet.pt_JetResDown=met.shiftedPt(pat::MET::JetResDown);
    mymet.phi_JetResUp=met.shiftedPhi(pat::MET::JetResUp);
    mymet.phi_JetResDown=met.shiftedPhi(pat::MET::JetResDown);
    mymet.pt_UnclusteredEnUp=met.shiftedPt(pat::MET::UnclusteredEnUp);
    mymet.pt_UnclusteredEnDown=met.shiftedPt(pat::MET::UnclusteredEnDown);
    mymet.phi_UnclusteredEnUp=met.shiftedPhi(pat::MET::UnclusteredEnUp);
    mymet.phi_UnclusteredEnDown=met.shiftedPhi(pat::MET::UnclusteredEnDown);
    mymet.genMET=met.genMET()->pt();
    (m->RecPFMet).push_back(mymet);

    (m->METfilters).clear();
    edm::Handle<edm::TriggerResults> trigFilter;
    iEvent.getByToken(metFiltersToken, trigFilter);

    const edm::TriggerNames &names = iEvent.triggerNames(*trigFilter);
    for (unsigned int i = 0, n = trigFilter->size(); i < n; ++i) {
        (m->METfilters)[names.triggerName(i)] = trigFilter->accept(i);
    }


/*/////////////////// Other variables in Run 1 ntuplizer /////////////
 
MVA MET

*/
}
