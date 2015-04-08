import FWCore.ParameterSet.Config as cms

process = cms.Process("NtupleProducer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring(
#'/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root',
'/store/mc/Phys14DR/VBF_HToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/147B369C-9F77-E411-B99D-00266CF9B184.root',
)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )
#process.MessageLogger.cerr.FwkReport.reportEvery=cms.untracked.int32(10)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

TRIGGERLIST = [#"HLT_*", #["HLT_Mu17_Mu8_v*", "HLT_Mu17_TkMu8_v*"] # to run on DATA/MC 2012 # "HLT_*" is a empty path
    "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1",
    "HLT_IsoMu17_eta2p1_v1",
    "HLT_IsoMu17_eta2p1_MediumIsoPFTau40_Trk1_eta2p1_Reg_v1",
    "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_v1",
    "HLT_IsoMu24_eta2p1_IterTrk01_v1",
    "HLT_IsoMu24_eta2p1_IterTrk02_v1",
    "HLT_IsoMu24_eta2p1_IterTrk02_LooseIsoPFTau20_v1",
    "HLT_Ele22_eta2p1_WP85_Gsf_LooseIsoPFTau20_v1",
    "HLT_Ele32_eta2p1_WP85_Gsf_v1",
    "HLT_Ele32_eta2p1_WP85_Gsf_LooseIsoPFTau20_v1",
    "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v1",
    "HLT_IsoMu16_eta2p1_CaloMET30_LooseIsoPFTau50_Trk30_eta2p1_v1",
    "HLT_IsoMu16_eta2p1_CaloMET30_v1",
    "HLT_Mu16_eta2p1_CaloMET30_v1",
    "HLT_LooseIsoPFTau50_Trk30_eta2p1_v1",
    "HLT_DoubleIsoMu17_eta2p1_v1",
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1",
    "HLT_Ele27_eta2p1_WP85_Gsf_LooseIsoPFTau20_v1",
    "HLT_Ele27_eta2p1_WP85_Gsf_v1"]


isMC=True
isEmbedded=False
isTT=True
isET=True
isMT=True
isEM=True
doPairMet=True
doSV=True
HLTProcessName="HLT"

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck", 
			ignoreTotal = cms.untracked.int32(1) # default is one 
		)

process.ntupler = cms.EDAnalyzer('NtupleProducer',
                                    #OutPut Filles (NTuples)
                                    HistOutFile=cms.untracked.string('output_Ntuples.root'),

                                    # Include or Exclude the objects
                                    Include_HPSTau=cms.bool(True),
                                    Include_Muon=cms.bool(True),
                                    Include_Electron=cms.bool(True),
                                    Include_Jet=cms.bool(True),
                                    Include_MET=cms.bool( True),
                                    Include_GenParticles=cms.bool(True),
                                    Include_HLT=cms.bool(True),
                                    Include_Vertex=cms.bool(True),
				    Include_PairMet=cms.bool(doPairMet),
				    Include_SV=cms.bool(doSV),
                                    Is_MC=cms.bool(isMC),
                                    Is_Embedded=cms.bool(isEmbedded),
                                    Is_TT=cms.bool(isTT),
                                    Is_EM=cms.bool(isEM),
                                    Is_MT=cms.bool(isMT),
                                    Is_ET=cms.bool(isET),

)

import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt

process.hltFilter = hlt.hltHighLevel.clone(
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    HLTPaths = TRIGGERLIST,
    andOr = cms.bool(True), # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = cms.bool(False) #if True: throws exception if a trigger path is invalid  
)

process.nEventsTotal = cms.EDProducer("EventCountProducer") 
process.nEventsPassTrigger = cms.EDProducer("EventCountProducer")

process.skimmedPatElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string("abs(eta) < 2.5")
)

process.skimmedPatMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string("abs(eta) < 2.5")
)

process.skimmedPatTaus = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("slimmedTaus"),
    cut = cms.string('pt > 19 && abs(eta) < 2.5 &&  tauID("decayModeFinding") > 0.5 && tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 4.0')
)

process.skimmedCorrPatTaus = cms.EDProducer('ScaleTauProducer',
    tauSrc=cms.InputTag('skimmedPatTaus'),
    NAME=cms.string('skimmedCorrPatTaus')
)


singlePatLeptons = cms.Sequence()
 
for eINDEX in range(10):
  eModuleName = "selectedPatElectrons%i" % (eINDEX)
  eModule = cms.EDProducer('SinglePatElectronProducer' ,
    electronSrc =cms.InputTag('skimmedPatElectrons'),
    INDEX = cms.uint32(eINDEX),
    NAME=cms.string(eModuleName)
    )
  setattr(process, eModuleName, eModule)
  if (isET or isEM):
     singlePatLeptons += eModule
 
for mINDEX in range(10):
  mModuleName = "selectedPatMuons%i" % (mINDEX)
  mModule = cms.EDProducer('SinglePatMuonProducer' ,
    muonSrc =cms.InputTag('skimmedPatMuons'),
    INDEX = cms.uint32(mINDEX),
    NAME=cms.string(mModuleName)
    )
  setattr(process, mModuleName, mModule)
  if (isMT or isEM):
     singlePatLeptons += mModule
 
for tINDEX in range(10):
  tModuleName = "selectedPatTaus%i" % (tINDEX)
  tModule = cms.EDProducer('SinglePatTauProducer' ,
    tauSrc =cms.InputTag('skimmedCorrPatTaus:skimmedCorrPatTaus:NtupleProducer'),
    INDEX = cms.uint32(tINDEX),
    NAME=cms.string(tModuleName)
    )
  setattr(process, tModuleName, tModule)
  if (isTT or isET or isMT):
     singlePatLeptons += tModule

for INDEX1 in range(9):
   for INDEX2 in range(INDEX1+1,10):
      tModuleName = "doublePatTaus%ix%i" % (INDEX1,INDEX2)
      tModule = cms.EDProducer('DoublePatTauProducer' ,
        tauSrc1 =cms.InputTag("selectedPatTaus%i:selectedPatTaus%i:NtupleProducer" % (INDEX1,INDEX1)),
        tauSrc2 =cms.InputTag("selectedPatTaus%i:selectedPatTaus%i:NtupleProducer" % (INDEX2,INDEX2)),
        NAME=cms.string(tModuleName)
      )
      setattr(process, tModuleName, tModule)
      if (isTT):
         singlePatLeptons += tModule


process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.load("RecoJets.JetProducers.ak4PFJets_cfi")
process.ak4PFJets.src = cms.InputTag("packedPFCandidates")

from JetMETCorrections.Configuration.DefaultJEC_cff import ak4PFJetsL1FastL2L3

process.load("RecoMET.METPUSubtraction.mvaPFMET_cff")
process.pfMVAMEt.srcPFCandidates = cms.InputTag("packedPFCandidates")
process.pfMVAMEt.srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices")

process.puJetIdForPFMVAMEt.jec =  cms.string('AK4PF')
process.puJetIdForPFMVAMEt.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
process.puJetIdForPFMVAMEt.rho = cms.InputTag("fixedGridRhoFastjetAll")

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

#if isMC:
#    process.GlobalTag.globaltag = 'PHYS14_25_V1::All' #MC in PHYS14
#else :
#    process.GlobalTag.globaltag = 'GR_70_V2_AN1::All'

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True)
)

pairWiseMvaMETs = cms.Sequence()

###########
# eTau METs
###########
 
for eINDEX in range(10):
  for tINDEX in range(10):
    metModuleName = "eTauMet%ix%i" % (eINDEX,tINDEX)
    eModuleName = "selectedPatElectrons%i:selectedPatElectrons%i:NtupleProducer" % (eINDEX,eINDEX)
    tModuleName = "selectedPatTaus%i:selectedPatTaus%i:NtupleProducer" % (tINDEX,tINDEX)
    metModule = process.pfMVAMEt.clone(
      srcLeptons = cms.VInputTag(cms.InputTag(eModuleName),cms.InputTag(tModuleName)),
      )
    setattr(process, metModuleName, metModule)
    metModule.minNumLeptons = cms.int32(2)
    if (isET):
       pairWiseMvaMETs += metModule
 
###########
# muTau METs
###########

for mINDEX in range(10):
  for tINDEX in range(10):
    metModuleName = "muTauMet%ix%i" % (mINDEX,tINDEX)
    mModuleName = "selectedPatMuons%i:selectedPatMuons%i:NtupleProducer" % (mINDEX,mINDEX)
    tModuleName = "selectedPatTaus%i:selectedPatTaus%i:NtupleProducer" % (tINDEX,tINDEX)
    metModule = process.pfMVAMEt.clone(
      srcLeptons = cms.VInputTag(cms.InputTag(mModuleName),cms.InputTag(tModuleName)),
      )
    setattr(process, metModuleName, metModule)
    metModule.minNumLeptons = cms.int32(2)
    if (isMT):
       pairWiseMvaMETs += metModule

###########
# tauTau METs
###########

for mINDEX in range(9):
  for tINDEX in range(mINDEX+1,10):
    metModuleName = "tauTauMet%ix%i" % (mINDEX,tINDEX)
    mModuleName = "selectedPatTaus%i:selectedPatTaus%i:NtupleProducer" % (mINDEX,mINDEX)
    tModuleName = "selectedPatTaus%i:selectedPatTaus%i:NtupleProducer" % (tINDEX,tINDEX)
    metModule = process.pfMVAMEt.clone(
      srcLeptons = cms.VInputTag(cms.InputTag(mModuleName),cms.InputTag(tModuleName)),
      )
    setattr(process, metModuleName, metModule)
    metModule.minNumLeptons = cms.int32(2)
    if (isTT):
       pairWiseMvaMETs += metModule

###########
# eMu METs
###########

for mINDEX in range(10):
  for tINDEX in range(10):
    metModuleName = "eMuMet%ix%i" % (mINDEX,tINDEX)
    mModuleName = "selectedPatElectrons%i:selectedPatElectrons%i:NtupleProducer" % (mINDEX,mINDEX)
    tModuleName = "selectedPatMuons%i:selectedPatMuons%i:NtupleProducer" % (tINDEX,tINDEX)
    metModule = process.pfMVAMEt.clone(
      srcLeptons = cms.VInputTag(cms.InputTag(mModuleName),cms.InputTag(tModuleName)),
      )
    setattr(process, metModuleName, metModule)
    metModule.minNumLeptons = cms.int32(2)
    if (isEM):
       pairWiseMvaMETs += metModule

pairWiseRecoil = cms.Sequence()

process.metRecoilCorPFMet = cms.EDProducer("MEtRecoilCorrectorProducer",
        genParticleTag = cms.InputTag("prunedGenParticles"),
        jetTag = cms.InputTag("slimmedJets"),
        metTag = cms.InputTag("pfMet"),
        electronTag = cms.InputTag("slimmedElectrons"),
        muonTag = cms.InputTag("slimmedMuons"),
        tauTag = cms.InputTag("slimmedTaus"),
        inputFileNamezmm42X = cms.FileInPath("NtupleProducer2015/Analysis/data/RecoilCorrector_v7/recoilfits/recoilfit_zmm53XRR_2012_njet.root"),
        inputFileNamedatamm = cms.FileInPath("NtupleProducer2015/Analysis/data/RecoilCorrector_v7/recoilfits/recoilfit_datamm53XRR_2012_njet.root"),
        inputFileNamewjets = cms.FileInPath("NtupleProducer2015/Analysis/data/RecoilCorrector_v7/recoilfits/recoilfit_wjets53X_20pv_njet.root"),
        inputFileNamezjets = cms.FileInPath("NtupleProducer2015/Analysis/data/RecoilCorrector_v7/recoilfits/recoilfit_zmm53X_2012_njet.root"),
        inputFileNamehiggs = cms.FileInPath("NtupleProducer2015/Analysis/data/RecoilCorrector_v7/recoilfits/recoilfit_higgs53X_20pv_njet.root"),
        numOfSigmas = cms.double(1.0),
        minJetPt = cms.double(30.0),
        verbose = cms.bool(False),
        isMC = cms.bool(isMC),
        idxTau = cms.int32(-1),
    )

##############   
# eTau Recoil
##############

for eINDEX in range(10):
  for tINDEX in range(10):
    recoilModuleName = "eTauRecoil%ix%i" % (eINDEX,tINDEX)
    eModuleName = "selectedPatElectrons%i:selectedPatElectrons%i:NtupleProducer" % (eINDEX,eINDEX)
    tModuleName = "selectedPatTaus%i:selectedPatTaus%i:NtupleProducer" % (tINDEX,tINDEX)
    recoilModule = process.metRecoilCorPFMet.clone(
        metTag = cms.InputTag("eTauMet%ix%i" % (eINDEX,tINDEX)),
        electronTag = cms.InputTag(eModuleName),
        muonTag = cms.InputTag(""),
        tauTag = cms.InputTag(tModuleName)
      )
    setattr(process, recoilModuleName, recoilModule)
    if (isET):
       pairWiseRecoil += recoilModule

############### 
# muTau Recoil
###############

for eINDEX in range(10):
  for tINDEX in range(10):
    recoilModuleName = "muTauRecoil%ix%i" % (eINDEX,tINDEX)
    eModuleName = "selectedPatMuons%i:selectedPatMuons%i:NtupleProducer" % (eINDEX,eINDEX)
    tModuleName = "selectedPatTaus%i:selectedPatTaus%i:NtupleProducer" % (tINDEX,tINDEX)
    recoilModule = process.metRecoilCorPFMet.clone(
        metTag = cms.InputTag("muTauMet%ix%i" % (eINDEX,tINDEX)),
        electronTag = cms.InputTag(""),
        muonTag = cms.InputTag(eModuleName),
        tauTag = cms.InputTag(tModuleName)
      )
    setattr(process, recoilModuleName, recoilModule)
    if (isMT):
       pairWiseRecoil += recoilModule

################   
# tauTau Recoil
################

for eINDEX in range(9):
  for tINDEX in range(eINDEX+1,10):
    recoilModuleName = "tauTauRecoil%ix%i" % (eINDEX,tINDEX)
    dModuleName = "doublePatTaus%ix%i:doublePatTaus%ix%i:NtupleProducer" % (eINDEX,tINDEX,eINDEX,tINDEX)
    recoilModule = process.metRecoilCorPFMet.clone(
        metTag = cms.InputTag("tauTauMet%ix%i" % (eINDEX,tINDEX)),
        electronTag = cms.InputTag(""),
        muonTag = cms.InputTag(""),
        tauTag = cms.InputTag(dModuleName)
      )
    setattr(process, recoilModuleName, recoilModule)
    if (isTT):
       pairWiseRecoil += recoilModule

#############   
# eMu Recoil
#############

for eINDEX in range(10):
  for tINDEX in range(10):
    recoilModuleName = "eMuRecoil%ix%i" % (eINDEX,tINDEX)
    eModuleName = "selectedPatElectrons%i:selectedPatElectrons%i:NtupleProducer" % (eINDEX,eINDEX)
    tModuleName = "selectedPatMuons%i:selectedPatMuons%i:NtupleProducer" % (tINDEX,tINDEX)
    recoilModule = process.metRecoilCorPFMet.clone(
        metTag = cms.InputTag("eMuMet%ix%i" % (eINDEX,tINDEX)),
        electronTag = cms.InputTag(eModuleName),
        muonTag = cms.InputTag(tModuleName),
        tauTag = cms.InputTag("")
      )
    setattr(process, recoilModuleName, recoilModule)
    if (isEM):
       pairWiseRecoil += recoilModule


pairWiseSV = cms.Sequence()

#############   
# eMu SV
#############

for eINDEX in range(10):
  for mINDEX in range(10):
    SVModuleName = "eMuSV%ix%i" % (eINDEX,mINDEX)
    eModuleName = "selectedPatElectrons%i:selectedPatElectrons%i:NtupleProducer" % (eINDEX,eINDEX)
    mModuleName = "selectedPatMuons%i:selectedPatMuons%i:NtupleProducer" % (mINDEX,mINDEX)
    SVModule = cms.EDProducer('SVmassProducer',
	lepton1Src=cms.InputTag(eModuleName),
 	lepton2Src=cms.InputTag(mModuleName),
        lepton1Index=cms.uint32(eINDEX),
        lepton2Index=cms.uint32(mINDEX),
	TYPE=cms.uint32(1),
        MC=cms.bool(isMC),
	NAME=cms.string(SVModuleName)
    )
    setattr(process, SVModuleName, SVModule)
    if (isEM):
       pairWiseSV += SVModule


for eINDEX in range(10):
  for tINDEX in range(10):
    SVModuleName = "eTauSV%ix%i" % (eINDEX,tINDEX)
    eModuleName = "selectedPatElectrons%i:selectedPatElectrons%i:NtupleProducer" % (eINDEX,eINDEX)
    tModuleName = "selectedPatTaus%i:selectedPatTaus%i:NtupleProducer" % (tINDEX,tINDEX)
    SVModule = cms.EDProducer('SVmassProducer',
        lepton1Src=cms.InputTag(eModuleName),
        lepton2Src=cms.InputTag(tModuleName),
        lepton1Index=cms.uint32(eINDEX),
        lepton2Index=cms.uint32(tINDEX),
        TYPE=cms.uint32(2),
        MC=cms.bool(isMC),
        NAME=cms.string(SVModuleName)
    )
    setattr(process, SVModuleName, SVModule)
    if (isET):
       pairWiseSV += SVModule

for mINDEX in range(10):
  for tINDEX in range(10):
    SVModuleName = "muTauSV%ix%i" % (mINDEX,tINDEX)
    eModuleName = "selectedPatMuons%i:selectedPatMuons%i:NtupleProducer" % (mINDEX,mINDEX)
    mModuleName = "selectedPatTaus%i:selectedPatTaus%i:NtupleProducer" % (tINDEX,tINDEX)
    SVModule = cms.EDProducer('SVmassProducer',
        lepton1Src=cms.InputTag(eModuleName),
        lepton2Src=cms.InputTag(mModuleName),
        lepton1Index=cms.uint32(mINDEX),
        lepton2Index=cms.uint32(tINDEX),
        MC=cms.bool(isMC),
        TYPE=cms.uint32(3),
        NAME=cms.string(SVModuleName)
    )
    setattr(process, SVModuleName, SVModule)
    if (isMT):
       pairWiseSV += SVModule

for eINDEX in range(9):
  for mINDEX in range(eINDEX+1,10):
    SVModuleName = "tauTauSV%ix%i" % (eINDEX,mINDEX)
    eModuleName = "selectedPatTaus%i:selectedPatTaus%i:NtupleProducer" % (eINDEX,eINDEX)
    mModuleName = "selectedPatTaus%i:selectedPatTaus%i:NtupleProducer" % (mINDEX,mINDEX)
    SVModule = cms.EDProducer('SVmassProducer',
        lepton1Src=cms.InputTag(eModuleName),
        lepton2Src=cms.InputTag(mModuleName),
	lepton1Index=cms.uint32(eINDEX),
        lepton2Index=cms.uint32(mINDEX),
	MC=cms.bool(isMC),
        TYPE=cms.uint32(4),
        NAME=cms.string(SVModuleName)
    )
    setattr(process, SVModuleName, SVModule)
    if (isTT):
       pairWiseSV += SVModule

if doSV and doPairMet:
   process.p = cms.Path(process.nEventsTotal*process.hltFilter*process.nEventsPassTrigger*process.skimmedPatElectrons*process.skimmedPatMuons*process.skimmedPatTaus*process.skimmedCorrPatTaus*singlePatLeptons*pairWiseMvaMETs*pairWiseRecoil*pairWiseSV*process.ntupler)
if not doSV and doPairMet:
   process.p = cms.Path(process.nEventsTotal*process.hltFilter*process.nEventsPassTrigger*process.skimmedPatElectrons*process.skimmedPatMuons*process.skimmedPatTaus*process.skimmedCorrPatTaus*singlePatLeptons*pairWiseMvaMETs*pairWiseRecoil*process.ntupler)
if not doSV and not doPairMet:
   process.p = cms.Path(process.nEventsTotal*process.hltFilter*process.nEventsPassTrigger*process.skimmedPatElectrons*process.skimmedPatMuons*process.skimmedPatTaus*process.skimmedCorrPatTaus*process.ntupler)


