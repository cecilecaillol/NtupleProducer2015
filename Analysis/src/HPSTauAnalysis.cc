#include "NtupleProducer2015/Analysis/interface/NtupleProducer.h"

void NtupleProducer::DoHPSTauAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
     (m->PreSelectedHPSTaus).clear();
 
     edm::Handle<reco::VertexCollection> vertices;
     iEvent.getByToken(vtxToken, vertices);
     if (vertices->empty()) return; // skip the event if no PV found
     const reco::Vertex &PV = vertices->front();
    
     edm::Handle<vector<pat::Tau>> taus;
     iEvent.getByToken(tauCollToken, taus);
     for (auto it = taus->cbegin(); it != taus->cend(); ++it) {
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
	//tau.leadChargedParticlePdgID=it->leadCand().pdgId();//FIXME
	//tau.dxy_Sig()
        tau.mva_e_pi = (it->leadPFChargedHadrCand().isNonnull() ? it->leadPFChargedHadrCand()->mva_e_pi() : 0.);
        tau.mva_pi_mu = (it->leadPFChargedHadrCand().isNonnull() ? it->leadPFChargedHadrCand()->mva_pi_mu() : 0.);
        tau.mva_e_mu = (it->leadPFChargedHadrCand().isNonnull() ? it->leadPFChargedHadrCand()->mva_e_mu() : 0.);
        tau.hcalEnergy = (it->leadPFChargedHadrCand().isNonnull() ? it->leadPFChargedHadrCand()->hcalEnergy() : 0.);
        tau.ecalEnergy = (it->leadPFChargedHadrCand().isNonnull() ? it->leadPFChargedHadrCand()->ecalEnergy() : 0.);
        tau.trackRefPt = (it->leadPFChargedHadrCand().isNonnull() ? it->leadPFChargedHadrCand()->pt() : 0.);

        tau.decayMode = it->decayMode();

	(m->PreSelectedHPSTaus).push_back(tau);
    }

/* //////////// Other variables in Run 1 ntuplizer //////////////////
 
        tauu.jetPt = itault->pfJetRef().get()->pt();
        tauu.jetEta = itault->pfJetRef().get()->eta();
        tauu.jetPhi = itault->pfJetRef().get()->phi();
        tauu.jetMass = itault->pfJetRef().get()->mass();

        tauu.leadTrackD0 = (itault->leadTrack().isNonnull() ? (itault->leadTrack())->dxy(PrVx->front().position()) : 0.);

        if( PrVx->size()!=0 && (itault->leadPFChargedHadrCand()).isNonnull() ){
          if( (itault->leadPFChargedHadrCand()->trackRef()).isNonnull() ){
	    tauu.dxy_PV = itault->leadPFChargedHadrCand()->trackRef()->dxy( (*PrVx)[0].position() );
	    tauu.dz_PV = itault->leadPFChargedHadrCand()->trackRef()->dz( (*PrVx)[0].position() );
          }
        }

	tauu.signalPiZeroCandidates = itault->signalPiZeroCandidates().size();
        tauu.ptSumChargedParticlesIsoCone = itault->isolationPFChargedHadrCandsPtSum();
        tauu.ptSumPhotonsIsoCone = itault->isolationPFGammaCandsEtSum();

        tauu.discriminationByDecayModeFinding = itault->tauID("decayModeFinding") > 0.5 ? true : false;

*/
}
