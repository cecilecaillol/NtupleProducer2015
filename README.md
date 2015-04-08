cmsrel CMSSW_7_2_3

cd CMSSW_7_2_3/src

git cms-addpkg FWCore/Version

git cms-addpkg PhysicsTools/PatAlgos

git-cms-merge-topic -u cms-met:72X-MetSig-150311

git-cms-merge-topic -u cms-met:72X-mvaMETForMiniAOD

cd RecoMET/METPUSubtraction/ ; git clone https://github.com/rfriese/RecoMET-METPUSubtraction data -b 72X-13TeV-Phys14_25_V4-26Mar15 

git clone https://github.com/veelken/SVfit_standalone TauAnalysis/SVfitStandalone

# EDIT THE RecoMET/METPUSubtraction/python/mvaPFMET_cff.py @ LINE 75
#inputFileNames = cms.PSet(
#        U     = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmet_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
#        DPhi  = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
#        CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
#        CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root')
#    ),
echo "Please edit mvaPFMET_cff.py file with correct weights!"

git clone https://github.com/cecilecaillol/NtupleProducer2015.git

scram b -j 20

cmsRun runTT_cfg.py


