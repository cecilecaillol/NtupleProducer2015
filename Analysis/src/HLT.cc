#include <FWCore/Framework/interface/Event.h>
#include "NtupleProducer2015/Analysis/interface/NtupleProducer.h"

void NtupleProducer::DoHLTAnalysis(const edm::Event& iEvent) {

    (m->HLT).clear();
    (m->TriggerObj).clear();
    edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;

    iEvent.getByToken(bitsToken, triggerBits);
    iEvent.getByToken(objectsToken, triggerObjects);

    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
	(m->HLT)[names.triggerName(i)] = triggerBits->accept(i);
    }
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
        obj.unpackPathNames(names);
        //std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
        //std::cout << "\t   Collection: " << obj.collection() << std::endl;
        //std::cout << "\t   Type IDs:   ";
        //for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
        //std::cout << "\t   Filters:    ";
        //for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
        std::vector<std::string> pathNamesAll  = obj.pathNames(false);
        std::vector<std::string> pathNamesLast = obj.pathNames(true);
        //std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
        for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
	    myTriggerObject trig;
            bool isBoth = obj.hasPathName( pathNamesAll[h], true, true ); 
            bool isL3   = obj.hasPathName( pathNamesAll[h], false, true ); 
            bool isLF   = obj.hasPathName( pathNamesAll[h], true, false ); 
	    trig.pt=obj.pt();
	    trig.eta=obj.eta();
	    trig.phi=obj.phi();
            trig.px=obj.px();
            trig.py=obj.py();
            trig.pz=obj.pz();
	    trig.path=pathNamesAll[h];
	    trig.isLastFilter= isBoth || isLF;
	    trig.isL3= isBoth || isL3;
	    (m->TriggerObj).push_back(trig);
        }
    }
}
