#include "NtupleProducer2015/Analysis/interface/NtupleProducer.h"

void NtupleProducer::DoPairWiseMetAnalysis(const edm::Event& iEvent) {

    using namespace std;
    using namespace reco;
    using namespace edm;
    using namespace pat;

    (m->PairMet_etau).clear();
    (m->PairMet_mutau).clear();
    (m->PairMet_tautau).clear();
    (m->PairMet_emu).clear();
    (m->PairRecoilMet_etau).clear();
    (m->PairRecoilMet_mutau).clear();
    (m->PairRecoilMet_tautau).clear();
    (m->PairRecoilMet_emu).clear();
    (m->PairMet_etau_sigMatrix_00).clear();
    (m->PairMet_etau_sigMatrix_01).clear();
    (m->PairMet_etau_sigMatrix_10).clear();
    (m->PairMet_etau_sigMatrix_11).clear();
    (m->PairMet_mutau_sigMatrix_00).clear();
    (m->PairMet_mutau_sigMatrix_01).clear();
    (m->PairMet_mutau_sigMatrix_10).clear();
    (m->PairMet_mutau_sigMatrix_11).clear();
    (m->PairMet_tautau_sigMatrix_00).clear();
    (m->PairMet_tautau_sigMatrix_01).clear();
    (m->PairMet_tautau_sigMatrix_10).clear();
    (m->PairMet_tautau_sigMatrix_11).clear();
    (m->PairMet_emu_sigMatrix_00).clear();
    (m->PairMet_emu_sigMatrix_01).clear();
    (m->PairMet_emu_sigMatrix_10).clear();
    (m->PairMet_emu_sigMatrix_11).clear();

    char buffer [50];
  if (IsET){
    for (int i = 0; i < 10; ++i){
        for (int j = 0; j < 10; ++j){
           edm::Handle <std::vector<reco::PFMET> >  mvamet;
	   sprintf (buffer, "eTauMet%ix%i",i,j);
           iEvent.getByLabel(buffer, mvamet);
           const reco::PFMET mvaMETpf =  (*mvamet)[0];
	   myMET mymet;
           mymet.pt = mvaMETpf.pt();
           mymet.px = mvaMETpf.px();
           mymet.py = mvaMETpf.py();
           mymet.phi = mvaMETpf.phi();
           mymet.et = mvaMETpf.et();
           (m->PairMet_etau).push_back(mymet);
	   reco::METCovMatrix covMVAMet;
	   covMVAMet = mvaMETpf.getSignificanceMatrix();
	   (m->PairMet_etau_sigMatrix_00).push_back(covMVAMet[0][0]);
           (m->PairMet_etau_sigMatrix_01).push_back(covMVAMet[0][1]);
           (m->PairMet_etau_sigMatrix_10).push_back(covMVAMet[1][0]);
           (m->PairMet_etau_sigMatrix_11).push_back(covMVAMet[1][1]);
	}
    }
  }
  if (IsMT){
    for (int i = 0; i < 10; ++i){
        for (int j = 0; j < 10; ++j){
           edm::Handle <std::vector<reco::PFMET> >  mvamet;
           sprintf (buffer, "muTauMet%ix%i",i,j);
           iEvent.getByLabel(buffer, mvamet);
           const reco::PFMET mvaMETpf =  (*mvamet)[0];
           myMET mymet;
           mymet.pt = mvaMETpf.pt();
           mymet.px = mvaMETpf.px();
           mymet.py = mvaMETpf.py();
           mymet.phi = mvaMETpf.phi();
           mymet.et = mvaMETpf.et();
           (m->PairMet_mutau).push_back(mymet);
           reco::METCovMatrix covMVAMet;
           covMVAMet = mvaMETpf.getSignificanceMatrix();
           (m->PairMet_mutau_sigMatrix_00).push_back(covMVAMet[0][0]);
           (m->PairMet_mutau_sigMatrix_01).push_back(covMVAMet[0][1]);
           (m->PairMet_mutau_sigMatrix_10).push_back(covMVAMet[1][0]);
           (m->PairMet_mutau_sigMatrix_11).push_back(covMVAMet[1][1]);
        }
    }
  }
  if (IsTT){
    for (int i = 0; i < 9; ++i){
        for (int j = i+1; j < 10; ++j){
           edm::Handle <std::vector<reco::PFMET> >  mvamet;
           sprintf (buffer, "tauTauMet%ix%i",i,j);
           iEvent.getByLabel(buffer, mvamet);
           const reco::PFMET mvaMETpf =  (*mvamet)[0];
           myMET mymet;
           mymet.pt = mvaMETpf.pt();
           mymet.px = mvaMETpf.px();
           mymet.py = mvaMETpf.py();
           mymet.phi = mvaMETpf.phi();
           mymet.et = mvaMETpf.et();
           (m->PairMet_tautau).push_back(mymet);
           reco::METCovMatrix covMVAMet;
           covMVAMet = mvaMETpf.getSignificanceMatrix();
           (m->PairMet_tautau_sigMatrix_00).push_back(covMVAMet[0][0]);
           (m->PairMet_tautau_sigMatrix_01).push_back(covMVAMet[0][1]);
           (m->PairMet_tautau_sigMatrix_10).push_back(covMVAMet[1][0]);
           (m->PairMet_tautau_sigMatrix_11).push_back(covMVAMet[1][1]);
        }
    }
  }
  if (IsEM){
    for (int i = 0; i < 10; ++i){
        for (int j = 0; j < 10; ++j){
           edm::Handle <std::vector<reco::PFMET> >  mvamet;
           sprintf (buffer, "eMuMet%ix%i",i,j);
           iEvent.getByLabel(buffer, mvamet);
           const reco::PFMET mvaMETpf =  (*mvamet)[0];
           myMET mymet;
           mymet.pt = mvaMETpf.pt();
           mymet.px = mvaMETpf.px();
           mymet.py = mvaMETpf.py();
           mymet.phi = mvaMETpf.phi();
           mymet.et = mvaMETpf.et();
           (m->PairMet_emu).push_back(mymet);
           reco::METCovMatrix covMVAMet;
           covMVAMet = mvaMETpf.getSignificanceMatrix();
           (m->PairMet_emu_sigMatrix_00).push_back(covMVAMet[0][0]);
           (m->PairMet_emu_sigMatrix_01).push_back(covMVAMet[0][1]);
           (m->PairMet_emu_sigMatrix_10).push_back(covMVAMet[1][0]);
           (m->PairMet_emu_sigMatrix_11).push_back(covMVAMet[1][1]);
        }
    }
  }
  if (IsET){
    for (int i = 0; i < 10; ++i){
        for (int j = 0; j < 10; ++j){
           edm::Handle <pat::METCollection>  mvamet;
           sprintf (buffer, "eTauRecoil%ix%i",i,j);
           iEvent.getByLabel(buffer,"N", mvamet);
	   if (mvamet.isValid()){
              const pat::MET mvaMETpf =  (*mvamet)[0];
              myMET mymet;
              mymet.pt = mvaMETpf.pt();
              mymet.px = mvaMETpf.px();
              mymet.py = mvaMETpf.py();
              mymet.phi = mvaMETpf.phi();
              mymet.et = mvaMETpf.et();
              (m->PairRecoilMet_etau).push_back(mymet);
	   }
        }
    }
  }
  if (IsMT){
    for (int i = 0; i < 10; ++i){
        for (int j = 0; j < 10; ++j){
           edm::Handle <pat::METCollection>  mvamet;
           sprintf (buffer, "muTauRecoil%ix%i",i,j);
           iEvent.getByLabel(buffer,"N", mvamet);
           if (mvamet.isValid()){
              const pat::MET mvaMETpf =  (*mvamet)[0];
              myMET mymet;
              mymet.pt = mvaMETpf.pt();
              mymet.px = mvaMETpf.px();
              mymet.py = mvaMETpf.py();
              mymet.phi = mvaMETpf.phi();
              mymet.et = mvaMETpf.et();
              (m->PairRecoilMet_mutau).push_back(mymet);
	   }
        }
    }
  }
  if (IsTT){
    for (int i = 0; i < 9; ++i){
        for (int j = i+1; j < 10; ++j){
           edm::Handle <pat::METCollection>  mvamet;
           sprintf (buffer, "tauTauRecoil%ix%i",i,j);
           iEvent.getByLabel(buffer,"N", mvamet);
           if (mvamet.isValid()){
              const pat::MET mvaMETpf =  (*mvamet)[0];
              myMET mymet;
              mymet.pt = mvaMETpf.pt();
              mymet.px = mvaMETpf.px();
              mymet.py = mvaMETpf.py();
              mymet.phi = mvaMETpf.phi();
              mymet.et = mvaMETpf.et();
              (m->PairRecoilMet_tautau).push_back(mymet);
	   }
        }
    }
  }
  if (IsEM){
    for (int i = 0; i < 10; ++i){
        for (int j = 0; j < 10; ++j){
           edm::Handle <pat::METCollection>  mvamet;
           sprintf (buffer, "eMuRecoil%ix%i",i,j);
           iEvent.getByLabel(buffer,"N", mvamet);
           if (mvamet.isValid()){
              const pat::MET mvaMETpf =  (*mvamet)[0];
              myMET mymet;
              mymet.pt = mvaMETpf.pt();
              mymet.px = mvaMETpf.px();
              mymet.py = mvaMETpf.py();
              mymet.phi = mvaMETpf.phi();
              mymet.et = mvaMETpf.et();
              (m->PairRecoilMet_emu).push_back(mymet);
	   }
        }
    }
  }
}

