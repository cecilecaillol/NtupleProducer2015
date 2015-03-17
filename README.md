cmsrel CMSSW_7_2_0
cd CMSSW_7_2_0/src
git cms-addpkg FWCore/Version
git cms-addpkg PhysicsTools/PatAlgos
git-cms-merge-topic -u cms-met:72X-13TeV-Training-30Jan15
git cms-merge-topic HuguesBrun:trigElecIdInCommonIsoSelection720
git clone https://github.com/rfriese/RecoMET-METPUSubtraction.git EgammaAnalysis/ElectronTools/data
git clone https://github.com/cecilecaillol/NtupleProducer2015.git
scram b -j 20
cmsRun runTT_cfg.py


