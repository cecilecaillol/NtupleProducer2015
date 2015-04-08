#include "NtupleProducer2015/Analysis/interface/NtupleProducer.h"

void NtupleProducer::DoVertexAnalysis(const edm::Event& iEvent) {
     (m->Vertex).clear();
     (m->runNumber) = 0;
     (m->eventNumber) = 0;
     (m->lumiNumber) = 0;
     m->PUInfo = 0;
     m->PUInfo_true = 0;
     m->PUInfo_Bunch0 = 0;

     m->runNumber = iEvent.id().run();
     m->eventNumber = iEvent.id().event();
     m->lumiNumber = iEvent.id().luminosityBlock();
     Handle<edm::View<PileupSummaryInfo> > pileupHandle;
     iEvent.getByToken(pileupToken, pileupHandle);
     for( auto & puInfoElement : *pileupHandle){
        if( puInfoElement.getBunchCrossing() == 0 ){
           m->PUInfo = puInfoElement.getPU_NumInteractions();
           m->PUInfo_true= puInfoElement.getTrueNumInteractions();
        }
     }

     edm::Handle<LHEEventProduct > LHEHandle;
     const LHEEventProduct* LHE = 0;
     iEvent.getByLabel("source",LHEHandle);
     if(LHEHandle.isValid())
 	LHE = LHEHandle.product();
     m->NUP=(LHE->hepeup()).NUP;

     edm::Handle<reco::VertexCollection> vertices;
     iEvent.getByToken(vtxToken, vertices);
     if (vertices->empty()) return; // skip the event if no PV found
     //const reco::Vertex &PV = vertices->front();
     VertexCollection::const_iterator firstGoodVertex = vertices->end();
     int firstGoodVertexIdx = 0;
     for (VertexCollection::const_iterator vtx = vertices->begin();
           vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
        myVertex vo;
        vo.Num_Vertex = vertices->size();
        vo.ndof=vtx->ndof();
        vo.normalizedChi2=vtx->chi2();
        vo.px=vtx->x();
        vo.py=vtx->y();
	vo.pz=vtx->z();
        vo.position_Rho=vtx->position().Rho();
        vo.tracksSize = vtx->tracksSize();
	if (firstGoodVertexIdx==0){
           m->Rho=vtx->position().Rho();
	}
        (m->Vertex).push_back(vo);
        /*if ( !(vtx->chi2()==0 && vtx->ndof()==0) && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0) {
              firstGoodVertex = vtx;
              break;
        }*/
     }
     //isFake, rho
     /*if ( firstGoodVertex==vertices->end() )
         return; // skip event if there are no good PVs
     pvNTracks_ = firstGoodVertex->nTracks();
*/

//isFake, rho, bunchcrossing0, embedding weight, isValid, puweight

}

