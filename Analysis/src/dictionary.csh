#!/bin/csh

#eval `scramv1 runtime -csh`
rootcint -f eventdict.cc -c -I${PWD}/../../.. \
            NtupleProducer2014/Analysis/interface/myobject.h \
            NtupleProducer2014/Analysis/interface/myevent.h \
            NtupleProducer2014/Analysis/interface/LinkDef.h
