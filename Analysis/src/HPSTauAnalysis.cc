#include "NtupleProducer2015/Analysis/interface/NtupleProducer.h"

void NtupleProducer::DoHPSTauAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
     (m->PreSelectedTaus).clear();
     (m->LooseTaus).clear();
 
     edm::Handle<reco::VertexCollection> vertices;
     iEvent.getByToken(vtxToken, vertices);
     if (vertices->empty()) return; // skip the event if no PV found
     const reco::Vertex &PV = vertices->front();
    
     edm::Handle<vector<pat::Tau>> loosetaus;
     iEvent.getByToken(tauCollToken, loosetaus);

     Handle<pat::TauCollection> taus;
     iEvent.getByLabel("skimmedCorrPatTaus","skimmedCorrPatTaus", taus);

     int ii=0;
     for (auto it = taus->cbegin(); it != taus->cend(); ++it) {
	myTau tau;
	tau.gen_index=ii;
        tau.pt = it->pt(); 
        tau.eta = it->eta();
        tau.phi = it->phi();
        tau.mass = it->mass();
        tau.mt = it->mt();
        tau.et = it->et();
        tau.px = it->px();
        tau.py = it->py();
        tau.pz = it->pz();
        tau.z = it->vz();
        tau.E = it->p();
        tau.Energy = it->energy();
        tau.charge = it->charge();

        tau.byLooseCombinedIsolationDeltaBetaCorr3Hits = it->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5 ? true : false;
        tau.byMediumCombinedIsolationDeltaBetaCorr3Hits = it->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits") > 0.5 ? true : false;
        tau.byTightCombinedIsolationDeltaBetaCorr3Hits = it->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits") > 0.5 ? true : false;
        tau.byRawCombinedIsolationDeltaBetaCorr3Hits = it->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
	tau.byVLooseIsolationMVA3oldDMwLT=it->tauID("byVLooseIsolationMVA3oldDMwLT");
        tau.byLooseIsolationMVA3oldDMwLT=it->tauID("byLooseIsolationMVA3oldDMwLT");
        tau.byMediumIsolationMVA3oldDMwLT=it->tauID("byMediumIsolationMVA3oldDMwLT");
        tau.byTightIsolationMVA3oldDMwLT=it->tauID("byTightIsolationMVA3oldDMwLT");
        tau.byVTightIsolationMVA3oldDMwLT=it->tauID("byVTightIsolationMVA3oldDMwLT");
        tau.discriminationByElectronMVA5Loose = it->tauID("againstElectronLooseMVA5") > 0.5 ? true : false;
        tau.discriminationByElectronMVA5Medium = it->tauID("againstElectronMediumMVA5") > 0.5 ? true : false;
        tau.discriminationByElectronMVA5VLoose = it->tauID("againstElectronVLooseMVA5") > 0.5 ? true : false;
        tau.discriminationByMuonLoose3 = it->tauID("againstMuonLoose3") > 0.5 ? true : false;
        tau.discriminationByMuonTight3 = it->tauID("againstMuonTight3") > 0.5 ? true : false;
        tau.IsolationChargedIsoPtSum = it->tauID("chargedIsoPtSum");
        tau.IsolationNeutralIsoPtSum = it->tauID("neutralIsoPtSum");
        tau.IsolationPUcorrPtSum = it->tauID("puCorrPtSum");
        //tau.discriminationByDecayModeFindingOldDMs = it->tauID("decayModeFindingOldDMs") > 0.5 ? true : false;
        tau.discriminationByDecayModeFindingNewDMs = it->tauID("decayModeFindingNewDMs") > 0.5 ? true : false;
        tau.discriminationByDecayModeFinding = it->tauID("decayModeFinding") > 0.5 ? true : false;

        tau.numChargedParticlesSignalCone = it->signalPFChargedHadrCands().size();
        tau.numNeutralHadronsSignalCone = it->signalPFNeutrHadrCands().size();
        tau.numPhotonsSignalCone = it->signalPFGammaCands().size();
        tau.numParticlesSignalCone = it->signalPFCands().size();
        tau.numChargedParticlesIsoCone = it->isolationPFChargedHadrCands().size();
        tau.numNeutralHadronsIsoCone = it->isolationPFNeutrHadrCands().size();
        tau.numPhotonsIsoCone = it->isolationPFGammaCands().size();
        tau.numParticlesIsoCone = it->isolationPFCands().size();

	tau.leadChargedParticlePt=it->leadCand()->pt();
        tau.trackRefPt = (it->leadPFChargedHadrCand().isNonnull() ? it->leadPFChargedHadrCand()->pt() : 0.);

        tau.decayMode = it->decayMode();
        tau.dxy=it->dxy();
        tau.dz=(PV.z()-it->vertex().z()) - ((PV.x()-it->vertex().x())*it->p4().x()+(PV.y()-it->vertex().y())*it->p4().y())/ it->p4().pt() *  it->p4().z()/ it->p4().pt();
        tau.zImpact=it->vertex().z() + 130./tan(it->theta());
        tau.isFirstVtx=(PV.z()==it->vertex().z());
	ii++;
	(m->PreSelectedTaus).push_back(tau);
    }
     for (auto it = loosetaus->cbegin(); it != loosetaus->cend(); ++it) {
        myTau tau;
        tau.pt = it->pt();
        tau.eta = it->eta();
        tau.phi = it->phi();
        tau.mass = it->mass();
        tau.mt = it->mt();
        tau.et = it->et();
        tau.px = it->px();
        tau.py = it->py();
        tau.pz = it->pz();
        tau.z = it->vz();
        tau.E = it->p();
        tau.Energy = it->energy();
        tau.charge = it->charge();
        tau.byLooseCombinedIsolationDeltaBetaCorr3Hits = it->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5 ? true : false;
        tau.byMediumCombinedIsolationDeltaBetaCorr3Hits = it->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits") > 0.5 ? true : false;
        tau.byTightCombinedIsolationDeltaBetaCorr3Hits = it->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits") > 0.5 ? true : false;
        tau.byRawCombinedIsolationDeltaBetaCorr3Hits = it->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
        tau.byVLooseIsolationMVA3oldDMwLT=it->tauID("byVLooseIsolationMVA3oldDMwLT");
        tau.byLooseIsolationMVA3oldDMwLT=it->tauID("byLooseIsolationMVA3oldDMwLT");
        tau.byMediumIsolationMVA3oldDMwLT=it->tauID("byMediumIsolationMVA3oldDMwLT");
        tau.byTightIsolationMVA3oldDMwLT=it->tauID("byTightIsolationMVA3oldDMwLT");
        tau.byVTightIsolationMVA3oldDMwLT=it->tauID("byVTightIsolationMVA3oldDMwLT");
        tau.discriminationByElectronMVA5Loose = it->tauID("againstElectronLooseMVA5") > 0.5 ? true : false;
        tau.discriminationByElectronMVA5Medium = it->tauID("againstElectronMediumMVA5") > 0.5 ? true : false;
        tau.discriminationByElectronMVA5VLoose = it->tauID("againstElectronVLooseMVA5") > 0.5 ? true : false;
        tau.discriminationByMuonLoose3 = it->tauID("againstMuonLoose3") > 0.5 ? true : false;
        tau.discriminationByMuonTight3 = it->tauID("againstMuonTight3") > 0.5 ? true : false;
        tau.IsolationChargedIsoPtSum = it->tauID("chargedIsoPtSum");
        tau.IsolationNeutralIsoPtSum = it->tauID("neutralIsoPtSum");
        tau.IsolationPUcorrPtSum = it->tauID("puCorrPtSum");
        tau.discriminationByDecayModeFindingNewDMs = it->tauID("decayModeFindingNewDMs") > 0.5 ? true : false;
        tau.discriminationByDecayModeFinding = it->tauID("decayModeFinding") > 0.5 ? true : false;
        tau.numChargedParticlesSignalCone = it->signalChargedHadrCands().size();
        tau.numNeutralHadronsSignalCone = it->signalNeutrHadrCands().size();
        tau.numPhotonsSignalCone = it->signalGammaCands().size();
        tau.numParticlesSignalCone = it->signalCands().size();
        tau.numChargedParticlesIsoCone = it->isolationChargedHadrCands().size();
        tau.numNeutralHadronsIsoCone = it->isolationNeutrHadrCands().size();
        tau.numPhotonsIsoCone = it->isolationGammaCands().size();
        tau.numParticlesIsoCone = it->isolationCands().size();
        tau.leadChargedParticlePt=it->leadCand()->pt();
        tau.trackRefPt = (it->leadChargedHadrCand().isNonnull() ? it->leadChargedHadrCand()->pt() : 0.);
        tau.decayMode = it->decayMode();
	tau.dxy=it->dxy();
        tau.dz=(PV.z()-it->vertex().z()) - ((PV.x()-it->vertex().x())*it->p4().x()+(PV.y()-it->vertex().y())*it->p4().y())/ it->p4().pt() *  it->p4().z()/ it->p4().pt();
	tau.zImpact=it->vertex().z() + 130./tan(it->theta());
	tau.isFirstVtx=(PV.z()==it->vertex().z());
        (m->LooseTaus).push_back(tau);
    }
/* //////////// Other variables in Run 1 ntuplizer //////////////////
 
        tauu.jetPt = itault->pfJetRef().get()->pt();
        tauu.jetEta = itault->pfJetRef().get()->eta();
        tauu.jetPhi = itault->pfJetRef().get()->phi();
        tauu.jetMass = itault->pfJetRef().get()->mass();

        tauu.leadTrackD0 = (itault->leadTrack().isNonnull() ? (itault->leadTrack())->dxy(PrVx->front().position()) : 0.);

	tauu.signalPiZeroCandidates = itault->signalPiZeroCandidates().size();
        tauu.ptSumChargedParticlesIsoCone = itault->isolationPFChargedHadrCandsPtSum();
        tauu.ptSumPhotonsIsoCone = itault->isolationPFGammaCandsEtSum();

*/
}
