NtupleProducer2014
==================
Directory with all packages

```shell
mkdir NtupleProducer_help
cd NtupleProducer_help/
cmsrel CMSSW_5_3_14
cd CMSSW_5_3_14/src/
cmsenv
```

MVA MET
```shell
git cms-merge-topic cms-analysis-tools:5_3_14-updateSelectorUtils
git cms-merge-topic -u TaiSakuma:53X-met-131120-01
git-cms-merge-topic -u cms-met:53X-MVaNoPuMET-20131217-01
```

EGamma tools
```shell
cvs co -r V00-00-09 EgammaAnalysis/ElectronTools
#cvs co -r jakob19April2013_2012ID EgammaAnalysis/ElectronTools
cvs co -r V09-00-01 RecoEgamma/EgammaTools
cd EgammaAnalysis/ElectronTools/data
cat download.url | xargs wget
cd ../../..
```

CMG tools
```shell
wget --no-check-certificate https://jez.web.cern.ch/jez/CMGTools.tgz
tar xzvf CMGTools.tgz
```

Ntuple producer
```shell
git clone https://github.com/cecilecaillol/NtupleProducer2014.git
```

Remove some conflicting packages
```shell
rm -rf DataFormats/TauReco
rm -rf RecoTauTag/RecoTau
rm -rf RecoTauTag/Configuration
rm -rf RecoTauTag/ImpactParameter
rm -rf RecoTauTag/TauTagTools
rm -rf PhysicsTools/PatAlgos
rm -rf DataFormats/PatCandidates
rm -rf GeneratorInterface
```

Directory with soft muon packages
```shell
cd ../../..
mkdir NtupleLLR
cd NtupleLLR
cmsrel CMSSW_5_3_14
cd CMSSW_5_3_14/src/
cmsenv
cvs co -r b5_3_X_cvMEtCorr_2013Feb22 JetMETCorrections/Type1MET
cvs co -r b5_3_X_cvMEtCorr_2013Feb26 CommonTools/ParticleFlow
cvs co -r caloMEtForSoftLep_26Apr_v2 -d LLRAnalysis UserCode/LLRAnalysis
```

Directory with tau ID packages
```shell
cd ../../..
mkdir NtupleProducer_main
cd NtupleProducer_main
cmsrel CMSSW_5_3_14
cd CMSSW_5_3_14/src/
cmsenv
```

New tau ID
```shell
git cms-merge-topic -u cms-tau-pog:CMSSW_5_3_X_boostedTaus_2013Dec17
cp -rf ../../../NtupleProducer_help/CMSSW_5_3_14/src/CMGTools .
cp -rf ../../../NtupleProducer_help/CMSSW_5_3_14/src/EgammaAnalysis/ .
cp -rf ../../../NtupleProducer_help/CMSSW_5_3_14/src/RecoBTag/ .
cp -rf ../../../NtupleProducer_help/CMSSW_5_3_14/src/RecoEgamma/ .
cp -rf ../../../NtupleProducer_help/CMSSW_5_3_14/src/RecoMET/ .
cp -rf ../../../NtupleProducer_help/CMSSW_5_3_14/src/DataFormats/ .
cp -rf ../../../NtupleProducer_help/CMSSW_5_3_14/src/RecoJets/ .
cp -rf ../../../NtupleProducer_help/CMSSW_5_3_14/src/PhysicsTools/ .
cp -rf ../../../NtupleProducer_help/CMSSW_5_3_14/src/JetMETCorrections/ .
cp -r ../../../NtupleProducer_help/CMSSW_5_3_14/src/NtupleProducer2014/ .
cp ../../../NtupleLLR/CMSSW_5_3_14/src/JetMETCorrections/Type1MET/plugins/CaloTowerMETcorrInputProducer.* JetMETCorrections/Type1MET/plugins/.
cp -rf ../../../NtupleLLR/CMSSW_5_3_14/src/LLRAnalysis .
git cms-addpkg DataFormats/HcalDetId
git cms-addpkg PhysicsTools/PatUtils
```

Add this line to JetMETCorrections/Type1MET/BuildFile.xml: 
```shell
<use   name="DataFormats/HcalDetId"/>
```

MVA MET configuration
```shell
cp -rf /afs/cern.ch/user/c/ccaillol/public/MVAMET_python/*.py RecoMET/METPUSubtraction/python/.
```

Correct TauMETAlgo.cc
```shell
vi JetMETCorrections/Type1MET/src/TauMETAlgo.cc 
:%s/PFCandidateRefVector/vector<reco::PFCandidatePtr>/g
:wq
```

PileUp reweighting
```shell
git cms-addpkg PhysicsTools/Utilities
```

Modify some files (RSS problem, duplicated weights, ...)
```shell
cp -f /afs/cern.ch/user/c/ccaillol/public/NtupleProducer/mvaPFMET*.py RecoMET/METPUSubtraction/python/.
cp -f /afs/cern.ch/user/v/veelken/public/forCecile3/JetMETCorrections/METPUSubtraction/data/*.db RecoMET/METPUSubtraction/data/.
cp -f /afs/cern.ch/user/c/ccaillol/public/NtupleProducer/puJetIDAlgo_cff.py CMGTools/External/python/.
cp -f /afs/cern.ch/user/v/veelken/public/forCecile3/PFJetIDSelectionFunctor.h PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h
```

Remove cout in PhysicsTools/Utilities/interface/LumiReWeighting.h.

```shell
scram b -j 20
```

