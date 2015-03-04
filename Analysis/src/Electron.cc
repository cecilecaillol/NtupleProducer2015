#include "NtupleProducer2015/Analysis/interface/NtupleProducer.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

void NtupleProducer::DoElectronAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
   (m->PreSelectedElectrons).clear();

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken, vertices);
    if (vertices->empty()) return; // skip the event if no PV found
    const reco::Vertex &PV = vertices->front();

   edm::Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronCollToken, electrons);
   for (const pat::Electron &el : *electrons) {
      myElectron elo;
        //elo.gen_index=index;
        elo.pt = el.pt();
        elo.eta = el.eta();
        elo.phi = el.phi();
        elo.px = el.px();
        elo.py = el.py();
        elo.pz = el.pz();
        elo.z = el.vz();
        elo.E = el.p();
        elo.et = el.et();
        elo.mass = el.mass();
        elo.mt = el.mt();
        elo.Energy = el.energy();
        elo.charge = el.charge();
	elo.eta_SC=el.superCluster()->eta() ;
        elo.phi_SC=el.superCluster()->phi();
        elo.rawE_SC  = el.superCluster()->rawEnergy(); 
        elo.preshowerE_SC = el.superCluster()->preshowerEnergy();
        elo.hcalOverEcal=el.hcalOverEcal();
        elo.deltaPhiSuperClusterTrackAtVtx = el.deltaPhiSuperClusterTrackAtVtx();
        elo.deltaEtaSuperClusterTrackAtVtx = el.deltaEtaSuperClusterTrackAtVtx();
        elo.sigmaIetaIeta = el.sigmaIetaIeta();
        elo.sigmaEtaEta = el.sigmaEtaEta();
        elo.ecalIso = el.ecalIso();
        elo.hcalIso = el.hcalIso();
        //elo.trackIso = el.trackIso();
        elo.caloIso = el.caloIso();
        elo.ecalEnergy = el.ecalEnergy();
        
  	GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();
	elo.pfIsoCharged = ( pfIso.sumChargedHadronPt );
	elo.pfIsoNeutral = ( pfIso.sumNeutralHadronEt );
	elo.pfIsoGamma = ( pfIso.sumPhotonEt );
	elo.pfIsoPU = ( pfIso.sumPUPt );
	elo.absiso = pfIso.sumChargedHadronPt + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
        elo.reliso=elo.absiso/el.pt();
        //elo.numLostHitEleInner = el.gsfTrack().trackerExpectedHitsInner().numberOfLostHits(); 
        elo.full5x5=el.full5x5_sigmaIetaIeta();
        elo.passConversionVeto=el.passConversionVeto();
        elo.dz=el.gsfTrack()->dz(PV.position());
        elo.dxy=el.gsfTrack()->dxy(PV.position());
	if( el.ecalEnergy() == 0 ){
           elo.ooEmooP = 1e30;
        }
	else if( !std::isfinite(el.ecalEnergy())){
           elo.ooEmooP = 1e30;
        }
	else{
           elo.ooEmooP = fabs(1.0/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy() );
        }
	//elo.MVAtrigID=myMVATrig->mvaValue(el,false); // doesn't work in CMSSW_7_3_0_pre3
        (m->PreSelectedElectrons).push_back(elo);
   }


/*  Other variables 
	MVA ID
        elo.numValidHitEle = (iElectron->gsfTrack().isNonnull() ? iElectron->gsfTrack()->numberOfValidHits() : 0.);
        elo.IP3D = iElectron->dB(pat::Electron::PV3D);
	vertex fit probability

*/
}
