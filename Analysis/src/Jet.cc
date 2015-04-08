#include <PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h>
#include <PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h>
#include "NtupleProducer2015/Analysis/interface/NtupleProducer.h"

void NtupleProducer::DoJetAnalysis(const edm::Event& iEvent) {

     (m->PreselectedJets).clear();

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
        myjet.E = it->p();
        myjet.Energy = it->energy();
        myjet.mass = it->mass();

	myjet.CSV = it->bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
        myjet.partonFlavour=it->partonFlavour();
 	myjet.puJetIdraw = it->userFloat("pileupJetId:fullDiscriminant");
	myjet.puJetId=false;
	if (fabs(it->eta())>=0 && fabs(it->eta())<2.5 && myjet.puJetIdraw>-0.63) myjet.puJetId=true;
	else if (fabs(it->eta())>=2.5 && fabs(it->eta())<2.75 && myjet.puJetIdraw>-0.60) myjet.puJetId=true;
	else if (fabs(it->eta())>=2.75 && fabs(it->eta())<3.0 && myjet.puJetIdraw>-0.55) myjet.puJetId=true;
	else if (fabs(it->eta())>=3.0 && fabs(it->eta())<5.2 && myjet.puJetIdraw>-0.45) myjet.puJetId=true;
	myjet.vtxMass = it->userFloat("vtxMass");
	myjet.vtxNtracks = it->userFloat("vtxNtracks");
	myjet.vtx3DVal = it->userFloat("vtx3DVal");
	myjet.vtx3DSig = it->userFloat("vtx3DSig");

	float chf = it->chargedHadronEnergy()/it->energy();
        float nhf = it->neutralHadronEnergy()/it->energy();
        float phf = it->neutralEmEnergy()/it->energy();
        float elf = it->chargedEmEnergy()/it->energy();
        float chm = it->chargedHadronMultiplicity();
        float npr = it->chargedMultiplicity() + it->neutralMultiplicity();

	myjet.jetId_Loose=(npr>1 and phf<0.99 and nhf<0.99) and (it->eta()>2.4 or (elf<0.99 and chf>0 and chm>0));
	myjet.jetId_Medium=(npr>1 and phf<0.95 and nhf<0.95) and (it->eta()>2.4 or (elf<0.99 and chf>0 and chm>0));
        myjet.jetId_Tight=(npr>1 and phf<0.90 and nhf<0.90) and (it->eta()>2.4 or (elf<0.99 and chf>0 and chm>0));

        (m->PreselectedJets).push_back(myjet);
     }

/* //////////////// Other variables in Run1 ntuplizer ////////////////////////
 
+ Add jet variations?

*/
}

