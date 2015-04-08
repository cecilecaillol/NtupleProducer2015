#include "NtupleProducer2015/Analysis/interface/NtupleProducer.h"

void NtupleProducer::DoSVAnalysis(const edm::Event& iEvent) {

  using namespace std;
  using namespace reco;
  using namespace edm;
  using namespace pat;

  (m->SV_etau).clear();
  (m->SV_mutau).clear();
  (m->SV_emu).clear();
  (m->SV_tautau).clear();

  char buffer [50];
  if (IsEM){
    for (int i = 0; i < 10; ++i){
        for (int j = 0; j < 10; ++j){
           edm::Handle <std::vector<float> >  sv;
           sprintf (buffer, "eMuSV%ix%i",i,j);
           edm::InputTag theSVLabel(buffer,buffer,"NtupleProducer");
	   iEvent.getByLabel(theSVLabel,  sv);
           mySV SVmass;
	   SVmass.mass =  (*sv)[0];
           SVmass.pt =  (*sv)[1];
           SVmass.eta =  (*sv)[2];
           SVmass.phi =  (*sv)[3];
           SVmass.met =  (*sv)[4];
	   (m->SV_emu).push_back(SVmass);
	}
    }
  }
  if (IsET){
    for (int i = 0; i < 10; ++i){
        for (int j = 0; j < 10; ++j){
           edm::Handle <std::vector<float> >  sv;
           sprintf (buffer, "eTauSV%ix%i",i,j);
           edm::InputTag theSVLabel(buffer,buffer,"NtupleProducer");
           iEvent.getByLabel(theSVLabel,  sv);
           mySV SVmass;
	   SVmass.mass =  (*sv)[0];
           SVmass.pt =  (*sv)[1];
           SVmass.eta =  (*sv)[2];
           SVmass.phi =  (*sv)[3];
           SVmass.met =  (*sv)[4];
           (m->SV_etau).push_back(SVmass);
        }
    }
  }
  if (IsMT){
    for (int i = 0; i < 10; ++i){
        for (int j = 0; j < 10; ++j){
           edm::Handle <std::vector<float> >  sv;
           sprintf (buffer, "muTauSV%ix%i",i,j);
           edm::InputTag theSVLabel(buffer,buffer,"NtupleProducer");
           iEvent.getByLabel(theSVLabel,  sv);
           mySV SVmass;
	   SVmass.mass =  (*sv)[0];
           SVmass.pt =  (*sv)[1];
           SVmass.eta =  (*sv)[2];
           SVmass.phi =  (*sv)[3];
           SVmass.met =  (*sv)[4];
           (m->SV_mutau).push_back(SVmass);
        }
    }
  }
  if (IsTT){
    for (int i = 0; i < 9; ++i){
        for (int j = i+1; j < 10; ++j){
           edm::Handle <std::vector<float> >  sv;
           sprintf (buffer, "tauTauSV%ix%i",i,j);
           edm::InputTag theSVLabel(buffer,buffer,"NtupleProducer");
           iEvent.getByLabel(theSVLabel,  sv);
           mySV SVmass;
	   SVmass.mass =  (*sv)[0];
           SVmass.pt =  (*sv)[1];
           SVmass.eta =  (*sv)[2];
           SVmass.phi =  (*sv)[3];
           SVmass.met =  (*sv)[4];
           (m->SV_tautau).push_back(SVmass);
        }
    }
  }
}


