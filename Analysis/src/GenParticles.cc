#include "NtupleProducer2015/Analysis/interface/NtupleProducer.h"

void NtupleProducer::DoGenParticlesAnalysis(const edm::Event& iEvent) {
    (m->RecGenParticles).clear();
     edm::Handle<edm::View<reco::GenParticle>> genPart;
     iEvent.getByToken(prunedGenToken, genPart);

     for(size_t i=0; i<genPart->size();i++){
	myGenObject gen;
	const Candidate * particle = &(*genPart)[i];
        gen.pdgId = particle->pdgId();
        gen.status = particle->status();
        gen.pt = particle->pt();
        gen.eta = particle->eta();
        gen.phi = particle->phi();
        gen.charge = particle->charge();
        gen.z = particle->vz();
        gen.mass = particle->mass();
        gen.px = particle->px();
        gen.px = particle->px();
        gen.px = particle->px();
        gen.E = particle->p();
        gen.Energy = particle->energy();

        gen.mod_pdgId = (particle->mother() != NULL ? particle->mother()->pdgId() : -1000);
        gen.mod_status = (particle->mother() != NULL ? particle->mother()->status() : -1000);
        gen.mod_pt = (particle->mother() != NULL ? particle->mother()->pt() : -1000);
        gen.mod_eta = (particle->mother() != NULL ? particle->mother()->eta() : -1000);
        gen.mod_phi = (particle->mother() != NULL ? particle->mother()->phi() : -1000);
        gen.mod_charge = (particle->mother() != NULL ? particle->mother()->charge() : -1000);
        gen.mod_z = (particle->mother() != NULL ? particle->mother()->vz() : -1000);
        gen.mod_mass = (particle->mother() != NULL ? particle->mother()->mass() : -1000);

        gen.Gmod_pdgId = ((particle->mother() != NULL && particle->mother()->mother() != NULL) ? particle->mother()->mother()->pdgId() : -1000);
        gen.Gmod_status = ((particle->mother() != NULL && particle->mother()->mother() != NULL) ? particle->mother()->mother()->status() : -1000);
        gen.Gmod_pt = ((particle->mother() != NULL && particle->mother()->mother() != NULL) ? particle->mother()->mother()->pt() : -1000);
        gen.Gmod_eta = ((particle->mother() != NULL && particle->mother()->mother() != NULL) ? particle->mother()->mother()->eta() : -1000);
        gen.Gmod_phi = ((particle->mother() != NULL && particle->mother()->mother() != NULL) ? particle->mother()->mother()->phi() : -1000);
        gen.Gmod_charge = ((particle->mother() != NULL && particle->mother()->mother() != NULL) ? particle->mother()->mother()->charge() : -1000);
        gen.Gmod_z = ((particle->mother() != NULL && particle->mother()->mother() != NULL) ? particle->mother()->mother()->vz() : -1000);
        gen.Gmod_mass = ((particle->mother() != NULL && particle->mother()->mother() != NULL) ? particle->mother()->mother()->mass() : -1000);

	(m->RecGenParticles).push_back(gen);
    }
/* //////////// Other variables in Run 1 ntuplizer /////////////////
 

	    // storing visible tau 4-vectors
	    if(abs(gen->pdgId())==15)
	      {
		const reco::GenParticleRefVector& mRefs = gen->daughterRefVector();
		reco::Particle::LorentzVector invisibleP4( 0.0, 0.0, 0.0, 0.0 );
		int decayMode = -1; // 0= hadronic, 1=electron, 2=muon
		unsigned int nNu = 0;
		for(reco::GenParticleRefVector::const_iterator imr = mRefs.begin(); imr!= mRefs.end(); ++imr) {
		  if(abs((*imr)->pdgId())==11) decayMode = 1;
		  if(abs((*imr)->pdgId())==13) decayMode = 2;
		  if(abs((*imr)->pdgId())==16 || abs((*imr)->pdgId())==14 || abs((*imr)->pdgId())==12)
		    {
		      invisibleP4+=(*imr)->p4();
		      nNu++;
		    }
		}
		if(nNu==1) decayMode = 0;
		reco::Particle::LorentzVector visibleP4 = gen->p4()-invisibleP4;
		myGenobject myVisTau;
		myVisTau.pdgId=gen->pdgId();
		myVisTau.status = gen->status();
		myVisTau.gen_index = index;
		myVisTau.pt = visibleP4.pt();
		myVisTau.eta = visibleP4.eta();
		myVisTau.phi = visibleP4.phi();
		myVisTau.charge = gen->charge();
		myVisTau.mass = visibleP4.M();
		myVisTau.decay_mode = decayMode;
		(m->RecGenTauVisible).push_back(myVisTau);
	      }
	    index++;
        }//loop over gens
*/
}

