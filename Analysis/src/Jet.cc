#include <PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h>
#include <PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h>
#include "NtupleProducer2015/Analysis/interface/NtupleProducer.h"

void NtupleProducer::DoJetAnalysis(const edm::Event& iEvent) {

     (m->RecPFJetsAK5).clear();

     edm::Handle<vector<pat::Jet>> jets;
     iEvent.getByToken(jetCollToken, jets);
     for (auto it = jets->cbegin(); it != jets->cend(); ++it) {
	myJet myjet;
        myjet.pt = it->pt();
        myjet.eta = it->eta();
        myjet.phi = it->phi();
        myjet.px = it->px();
        myjet.py = it->py();
        myjet.pz = it->pz();
        myjet.z = it->vz();
        myjet.E = it->p();
        myjet.Energy = it->energy();
        myjet.mass = it->mass();
        myjet.mt = it->mt();
        myjet.et = it->et();

	myjet.CSV = it->bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
        myjet.partonFlavour=it->partonFlavour();
 	myjet.puJetId = it->userFloat("pileupJetId:fullDiscriminant");
	myjet.vtxMass = it->userFloat("vtxMass");
	myjet.vtxNtracks = it->userFloat("vtxNtracks");
	myjet.vtx3DVal = it->userFloat("vtx3DVal");
	myjet.vtx3DSig = it->userFloat("vtx3DSig");

        (m->RecPFJetsAK5).push_back(myjet);
     }

/* //////////////// Other variables in Run1 ntuplizer ////////////////////////
 
        PFJetIDSelectionFunctor jetIDLoose(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE);
        PFJetIDSelectionFunctor jetIDTight(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT);

        myjet.jetId_loose = jetIDLoose(*jet, ret_loose);
        myjet.jetId_tight = jetIDTight(*jet, ret_tight);

        myjet.bDiscriminatiors_JP = jet->bDiscriminator("jetProbabilityBJetTags");
        myjet.bDiscriminatiors_TCHPT = jet->bDiscriminator("trackCountingHighPurBJetTags");

	myjet.puJetIdLoose= PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kLoose );
	myjet.puJetIdMedium = PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kMedium );
	myjet.puJetIdTight = PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kTight );

+ Add jet variations?

*/
}

