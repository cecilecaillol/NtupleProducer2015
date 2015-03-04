#include "NtupleProducer2015/Analysis/interface/NtupleProducer.h"

void NtupleProducer::DoMuonAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    (m->PreSelectedMuons).clear();

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken, vertices);
    if (vertices->empty()) return; // skip the event if no PV found
    const reco::Vertex &PV = vertices->front();

     edm::Handle<vector<pat::Muon>> muons;
     iEvent.getByToken(muonCollToken, muons);
     for (auto it = muons->cbegin(); it != muons->cend(); ++it) {
	myMuon muo;
	muo.pt=it->pt();
        muo.px=it->px();
        muo.py=it->py();
        muo.pz=it->pz();
        muo.phi=it->phi();
        muo.eta=it->eta();
        muo.E = it->p();
        muo.Energy = it->energy();
        muo.mass = it->mass();
        muo.mt = it->mt();
        muo.et = it->et();
        muo.charge = it->charge();
	muo.isGlobalMuon = it->isGlobalMuon();
        muo.z = it->vz();

	muo.segmentCompatibility=it->segmentCompatibility();
	//muo.validFraction=it->innerTrack()->validFraction();
	muo.trkKink=it->combinedQuality().trkKink;
	muo.chi2LocalPosition=it->combinedQuality().chi2LocalPosition;
	//muo.normalizedChi2=std::min(it->globalTrack()->normalizedChi2(), 10.99);
	//muo.numberOfValidMuonHits=(std::min(it->globalTrack()->hitPattern().numberOfValidMuonHits(), (int)39));
	muo.numMatchStation=(std::min(it->numberOfMatchedStations(), (int)5));
	//muo.numberOfValidPixelHits=(std::min(it->innerTrack()->hitPattern().numberOfValidPixelHits(), (int)5));//FIXME
	//muo.intrkLayerMeasure=(std::min(it->innerTrack()->hitPattern().trackerLayersWithMeasurement(), (int)14));
	if (it->isIsolationValid())
	{
	   reco::MuonPFIsolation pfR04 = it->pfIsolationR04();
	   muo.pfIsoCharged=(std::min(pfR04.sumChargedHadronPt, (float)9.9));
	   muo.pfIsoNeutral=(std::min(pfR04.sumNeutralHadronEt, (float)9.9));
	   muo.pfIsoGamma=(std::min(pfR04.sumPhotonEt, (float)9.9));
	   muo.pfIsoPU=(std::min(pfR04.sumPUPt, (float)9.9));
	   muo.absiso= pfR04.sumChargedHadronPt + std::max(0., pfR04.sumNeutralHadronEt+pfR04.sumPhotonEt-0.5*pfR04.sumPUPt);
	   muo.reliso=muo.absiso/it->pt();
	}
	muo.isLooseMuon=it->isLooseMuon();
	muo.isTightMuon=it->isTightMuon(PV);
	muo.isSoftMuon = it->isSoftMuon();
	muo.isHighPtMuon = it-> isHighPtMuon();
        muo.isTrackerMuon = it->isTrackerMuon();
        muo.isStandAloneMuon = it->isStandAloneMuon();
	muo.isPFMuon = it->isPFMuon();
	muo.dz=it->muonBestTrack()->dz(PV.position());
	muo.dxy=it->muonBestTrack()->dxy(PV.position());
        (m->PreSelectedMuons).push_back(muo);
     }


/*  ///////////// Other variables in Run 1 ntuplizer ///////////////////
        
	double DepInTracker = amuon->isolationR03().sumPt;
        double DepInEcal = amuon->isolationR03().emEt;
        double DepInHcal = amuon->isolationR03().hadEt;
	muo.gen_index=qq;
        muo.dB = amuon->dB();
	muo.IP3D = amuon->dB(pat::Muon::PV3D);
        muo.DepositR03Ecal = DepInEcal;
        muo.DepositR03Hcal = DepInHcal;
        muo.DepositR03TrackerOfficial = DepInTracker;
        muo.GlobalMuonPromptTight = muon::isGoodMuon(*amuon, muon::GlobalMuonPromptTight);
        muo.TMOneStationLoose = muon::isGoodMuon(*amuon, muon::TMOneStationLoose);
        muo.TM2DCompatibilityLoose = muon::isGoodMuon(*amuon, muon::TM2DCompatibilityLoose);
        muo.trkLayerMeasure =(amuon->track().isNonnull() ? amuon->track()->hitPattern().trackerLayersWithMeasurement():0);
        muo.intrkLayerpixel =(amuon->innerTrack().isNonnull() ? amuon->innerTrack()->hitPattern().numberOfValidPixelHits():0);
        muo.dxy_in =(amuon->innerTrack().isNonnull() ? amuon->innerTrack()->dxy(primVertex.position()) : 0.);
        muo.dZ_in =(amuon->innerTrack().isNonnull() ? amuon->innerTrack()->dz(primVertex.position()) : 0.);
	muo.normalizedChi2_innTrk = (amuon->innerTrack().isNonnull() ? amuon->innerTrack()->normalizedChi2() : 0);
        muo.numberOfValidMuonHits_innTrk = (amuon->innerTrack().isNonnull() ? amuon->innerTrack()->hitPattern().numberOfValidMuonHits() : 0);
        muo.numberOfHits_innTrk = (amuon->innerTrack().isNonnull() ? amuon->innerTrack()->hitPattern().numberOfHits() : 0);
        muo.pfIsoAll = (*isoAllMuMap)[muRef];
*/
}

