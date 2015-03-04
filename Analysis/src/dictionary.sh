#!/bin/csh

#eval `scramv1 runtime -csh`
rootcint -f eventdict.cc -c -I${PWD}/../../.. \
            NtupleProducer2015/Analysis/interface/myVertex.h \
            NtupleProducer2015/Analysis/interface/myTriggerObject.h \
            NtupleProducer2015/Analysis/interface/myGenObject.h \
            NtupleProducer2015/Analysis/interface/myJet.h \
            NtupleProducer2015/Analysis/interface/myTau.h \
            NtupleProducer2015/Analysis/interface/myElectron.h \
            NtupleProducer2015/Analysis/interface/myMuon.h \
            NtupleProducer2015/Analysis/interface/myMET.h \
            NtupleProducer2015/Analysis/interface/myevent.h \
            NtupleProducer2015/Analysis/interface/LinkDef.h

