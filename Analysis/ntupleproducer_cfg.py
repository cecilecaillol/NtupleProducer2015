import FWCore.ParameterSet.Config as cms

process = cms.Process("NtupleProducer")

process.out = cms.OutputModule("PoolOutputModule",
                               fileName=cms.untracked.string('PATLayer1_Output.fromAOD_full.root'),
                               SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('p')),
                               outputCommands=cms.untracked.vstring('drop *')
                               )
process.out.outputCommands.append('keep *_metNoHFresidualCorrected*_*_*') #  caloMEtNoHF corrected for residaul corrections between data and MC
process.out.outputCommands.append('keep *_genTauMatchedCaloJet*_*_*') # calo-deposits around generated taus
process.out.outputCommands.append('keep *_correctedCaloMEtNoHF_*_*') # caloMEtNoHF corrected for presence of calo-deposits around generated taus
process.out.outputCommands.append('keep *_uncorrectedL1ETM_*_*') # l1extraParticles_MET converted into caloMet data format
process.out.outputCommands.append('keep *_correctedL1ETM_*_*') # l1extraParticles_MET converted into caloMet data format and corrected for presence of calo-deposits around generated taus
process.out.outputCommands.append('keep *_l1extraParticles_MET_*') # L1ETM = MEt at L1-trigger
process.out.outputCommands.append('keep *_metNoHF_*_*') # caloMET w/o HF
process.out.outputCommands.append('keep *_generator_minVisPtFilter_*')

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("RecoEcal.EgammaClusterProducers.ecalClusteringSequence_cff")
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("PhysicsTools/PatAlgos/patSequences_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi")
process.load("NtupleProducer2014.SinglePatElectronProducer.singlepatelectronproducer_cfi")
process.load("NtupleProducer2014.SinglePatMuonProducer.singlepatmuonproducer_cfi")
process.load("NtupleProducer2014.SinglePatTauProducer.singlepattauproducer_cfi")
process.load("NtupleProducer2014.ScaleTauProducer.scaletauproducer_cfi")
process.load("NtupleProducer2014.DoublePatTauProducer.doublepattauproducer_cfi")
process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesHZZElectronIdentificationV06_cfi")
process.load("NtupleProducer2014.Analysis.New_hTozzTo4leptonsPFIsolationProducer_cff")
process.load("NtupleProducer2014.Analysis.JES_Uncertainty_FR_cff")
process.MessageLogger = cms.Service("MessageLogger")
process.load("RecoLocalCalo/EcalRecAlgos/EcalSeverityLevelESProducer_cfi")
process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi') #new location of the code
process.load("CMGTools.External.pujetidsequence_cff") #load PU JetID sequence
process.load("EventFilter.HcalRawToDigi.hcallasereventfilter2012_cff") #laser correction
process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff") #MET correction
#process.load("RecoMET.METPUSubtraction.mvaPFMET_cff")
process.load("RecoMET.METPUSubtraction.mvaPFMET_db_cfi")
process.load("RecoMET.METPUSubtraction.mvaPFMET_leptons_cff") #new location of the code
process.load("RecoMET.METPUSubtraction.mvaPFMET_leptons_ditau_cff")
process.load("RecoMET.METPUSubtraction.mvaPFMET_leptons_mutau_cff")
process.load("RecoMET.METPUSubtraction.mvaPFMET_leptons_etau_cff")
process.load("RecoMET.METPUSubtraction.mvaPFMET_leptons_emu_cff")
process.load("LLRAnalysis.TauTauStudies.sumCaloTowersInEtaSlices_cfi")
process.load("LLRAnalysis.TauTauStudies.sumCaloTowersInEtaSlices_cfi")
process.load('RecoBTag.Configuration.RecoBTag_cff')  #to run btagging again for embedded samples
process.load('RecoJets.JetAssociationProducers.ak5JTA_cff') #to run btagging again for embedded samples

#################################################   Samples and GlobalTag   ############################
#Source File
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
#'file:/afs/cern.ch/work/a/abdollah/Samples/VBF_HToTauTau_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_169CDF2F-A4ED-E111-8318-0017A477002C.root'
#'file:/localgrid_mnt/localgrid/abdollah/Samples/ZH_HToTauTau_M-125_lepdecay_8TeV-pythia6-tauola_PU_S10_START53_V7C-v1_SYNC.root'
    #'file:/afs/cern.ch/work/c/ccaillol/DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_PU_S7_START52_V9_v2.root'
    #'file:/tmp/abdollah/Run2012A_DoubleMu_AOD_13Jul2012-v1_F87222D4-D6CF-E111-B6ED-0026189438BC.root'
    #'file:/afs/cern.ch/user/a/abdollah/public/For_Olivier/pickevents_1_1_3Gd_SUSYBB_120GeV_285Events.root'
    'root://xrootd2.ihepa.ufl.edu//store/results/tau-pflow/DoubleMuParked/StoreResults-Run2012C_22Jan2013_v1_PFembedded_trans1_tau122_ptmu1_18mu2_8_v1-5ef1c0fd428eb740081f19333520fdc8/DoubleMuParked/USER/StoreResults-Run2012C_22Jan2013_v1_PFembedded_trans1_tau122_ptmu1_18mu2_8_v1-5ef1c0fd428eb740081f19333520fdc8/0000/00400AD0-31E8-E211-9652-0023AEFDE8D8.root'
    )
                            )
process.maxEvents = cms.untracked.PSet(
                                       input=cms.untracked.int32(20)
                                       )

#isMC = True      
isMC = False     
isEmbedded = True
#isEmbedded = False 

isTT = True
isMT = True
isET = True
isEM = True

process.myProducerLabel = cms.EDProducer('NtupleProducer')
#################################

if isMC:
    HLTProcessName = 'HLT'
    process.GlobalTag.globaltag = cms.string('START53_V23::All') #updated JEC

else:
    HLTProcessName = 'HLT'
    process.GlobalTag.globaltag = cms.string('FT_53_V21_AN4::All') 
    process.patPFMet.addGenMET = cms.bool(False)

#################################################   EDANALYZER   ##################################
process.myanalysis = cms.EDAnalyzer("NtupleProducer",

                                    #OutPut Filles (NTuples)
                                    HistOutFile=cms.untracked.string('output_Ntuples.root'),

                                    # Include or Exclude the objects
                                    Include_HPSTau=cms.bool(True),
                                    Include_Muon=cms.bool(True),
                                    Include_Electron=cms.bool(True),
                                    Include_Jet=cms.bool(True),
                                    Include_JetCorrection=cms.bool(False),
                                    Include_MET=cms.bool( True),
                                    Include_L1=cms.bool(False),#new to have L1 trigger objects
                                    Include_MET_Uncertaity=cms.bool(False),
                                    Include_GenPartiles=cms.bool(True),
                                    Include_HLT=cms.bool(True),
                                    Include_Vertex=cms.bool(True),
				    Include_PairWiseMet=cms.bool(True),
                                    Is_MC=cms.bool(isMC),
                                    Is_Embedded=cms.bool(isEmbedded),
				    Is_TT=cms.bool(isTT),
                                    Is_EM=cms.bool(isEM),
                                    Is_MT=cms.bool(isMT),
                                    Is_ET=cms.bool(isET),
                                    
                                    #Storing only certain trigger strings
                                    filterTriggerResults=cms.bool(True),

                                    #Vertices and Tracks
                                    vertices=cms.InputTag('offlinePrimaryVertices'),
                                    tracks=cms.InputTag("generalTracks"),

                                    #MET
                                    met=cms.InputTag("met"),
                                    PFmet=cms.InputTag("pfMet"),
                                    tcmet=cms.InputTag("tcMet"),
                                    Type1CorMET=cms.InputTag("patType1CorrectedPFMet"),
                                    MVAmet=cms.InputTag("pfMEtMVA"),
				    MVAmet_ditau=cms.InputTag("pfMEtMVAditau"),
                                    MVAmet_mutau=cms.InputTag("pfMEtMVAmutau"),
                                    MVAmet_etau=cms.InputTag("pfMEtMVAetau"),
                                    MVAmet_emu=cms.InputTag("pfMEtMVAemu"),
				    RecoilCorrectedPFMET=cms.InputTag("metRecoilCorPFMet","N"),
				    RecoilCorrectedMVAMetditau=cms.InputTag("metRecoilCorMVAditau","N"),
				    RecoilCorrectedMVAMetmutau=cms.InputTag("metRecoilCorMVAmutau","N"),
				    RecoilCorrectedMVAMetetau=cms.InputTag("metRecoilCorMVAetau","N"),
				    RecoilCorrectedMVAMetemu=cms.InputTag("metRecoilCorMVAemu","N"),
			
                                    #Jets
                                    bjets=cms.InputTag("trackCountingHighPurBJetTags"),
                                    PFAK5=cms.InputTag("selectedPatJets"),
                                    rhoJetsLabel=cms.InputTag("kt6PFJets", "rho"),
				    rhoCenChargedPU=cms.InputTag("kt6PFJetsCentralChargedPileUp", "rho"),
				    rhoCenNeutral=cms.InputTag("kt6PFJetsCentralNeutral", "rho"),
				    rhoCenNeutralTight=cms.InputTag("kt6PFJetsCentralNeutralTight", "rho"),

                                    #Leptons
                                    preselectedelectrons=cms.InputTag("skimmedPatElectrons"),
                                    preselectedmuons=cms.InputTag("skimmedPatMuons"),
                                    preselectedHPSTaus=cms.InputTag("skimmedPatTausLT"),
                                    vertexCollectionForLeptonIP=cms.InputTag("selectPrimaryVertex"),
                                    tauPtCut=cms.double(15.0),

                                    # CIC electron Identification
                                    eleID_VeryLooseTag=cms.InputTag("eidVeryLoose"),
                                    eleID_LooseTag=cms.InputTag("eidLoose"),
                                    eleID_MediumTag=cms.InputTag("eidMedium"),
                                    eleID_TightTag=cms.InputTag("eidTight"),

                                    #Trigger
                                    srcTriggerResults=cms.InputTag("TriggerResults", "", HLTProcessName),

                                    # MC Information
                                    PileUpInfo=cms.InputTag("addPileupInfo"),
                                    genParticlesInfo=cms.InputTag("genParticles"),

                                    #Trigger and TriggerMatching
                                    triggerEvent=cms.InputTag("patTriggerEvent"),
                                    tauMatch_Mu17Tau20=cms.string('tauTriggerMu17Tau20'),
                                    muMatch_Mu17Tau20=cms.string('muTriggerMu17Tau20'),
                                    tauMatch_Mu18Tau25=cms.string('tauTriggerMu18Tau25'),
                                    muMatch_Mu18Tau25=cms.string('muTriggerMu18Tau25'),
                                    tauMatch_Ele20Tau20=cms.string('tauTriggerEle20Tau20'),
                                    eleMatch_Ele20Tau20=cms.string('eleTriggerEle20Tau20'),
                                    tauMatch_Ditau35=cms.string('tauTriggerDitau35'),
                                    tauMatch_Ditau30Jet30=cms.string('tauTriggerDitau30Jet30'),
                                    jetMatch_Ditau30Jet30=cms.string('jetTriggerDitau30Jet30'),
                                    muMatch_Mu24=cms.string('muTriggerMu24'),
				    muMatch_SoftMu=cms.string('muTriggerSoftMu'),
                                    tauMatch_SoftMu=cms.string('tauTriggerSoftMu'),
				    muMatch_Ele8Mu17=cms.string('muTriggerEle8Mu17'),
                                    muMatch_Ele17Mu8=cms.string('muTriggerEle17Mu8'),
                                    eleMatch_Ele17Mu8=cms.string('eleTriggerEle17Mu8'),
                                    muMatch_Mu17Mu8=cms.string('muTriggerMu17Mu8'),
                                    eleMatch_Ele17Ele8=cms.string('eleTriggerEle17Ele8'),
                                    eleMatch_Ele8Mu17=cms.string('eleTriggerEle8Mu17'),
                                    eleMatch_SoftEle=cms.string('eleTriggerSoftEle'),
                                    tauMatch_SoftEle=cms.string('tauTriggerSoftEle'),

                                    puJetIdFlag=cms.InputTag("puJetMva","fullId"),
				    rhoProducer = cms.InputTag('kt6PFJets','rho') #new
                                    )
#################################################   PFIsolation  ################################
from CommonTools.ParticleFlow.pfParticleSelection_cff import *

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso
process.eleIsoSequence = setupPFElectronIso(process, 'skimmedPatElectrons')

process.muPFIsoDepositCharged.src=cms.InputTag("skimmedPatMuons")
process.muPFIsoDepositChargedAll.src=cms.InputTag("skimmedPatMuons")
process.muPFIsoDepositGamma.src=cms.InputTag("skimmedPatMuons")
process.muPFIsoDepositNeutral.src=cms.InputTag("skimmedPatMuons")
process.muPFIsoDepositPU.src=cms.InputTag("skimmedPatMuons")
process.muPFIsoValueCharged03.src=cms.InputTag("skimmedPatMuons")
process.muPFIsoValueChargedAll03.src=cms.InputTag("skimmedPatMuons")
process.muPFIsoValueGamma03.src=cms.InputTag("skimmedPatMuons")
process.muPFIsoValueNeutral03.src=cms.InputTag("skimmedPatMuons")
process.muPFIsoValueGammaHighThreshold03.src=cms.InputTag("skimmedPatMuons")
process.muPFIsoValueNeutralHighThreshold03.src=cms.InputTag("skimmedPatMuons")
process.muPFIsoValuePU03.src=cms.InputTag("skimmedPatMuons")
process.muPFIsoValueCharged04.src=cms.InputTag("skimmedPatMuons")
process.muPFIsoValueChargedAll04.src=cms.InputTag("skimmedPatMuons")
process.muPFIsoValueGamma04.src=cms.InputTag("skimmedPatMuons")
process.muPFIsoValueNeutral04.src=cms.InputTag("skimmedPatMuons")
process.muPFIsoValueGammaHighThreshold04.src=cms.InputTag("skimmedPatMuons")
process.muPFIsoValueNeutralHighThreshold04.src=cms.InputTag("skimmedPatMuons")
process.muPFIsoValuePU04.src=cms.InputTag("skimmedPatMuons")

# cone vetos as used in Htautau
process.elPFIsoValueChargedAll04NoPFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')
process.elPFIsoValueChargedAll04PFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')
process.elPFIsoValueChargedAll03NoPFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')
process.elPFIsoValueChargedAll03PFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')

process.elPFIsoValueGamma04NoPFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.08)','EcalEndcaps:ConeVeto(0.08)')
process.elPFIsoValueGamma04PFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.08)','EcalEndcaps:ConeVeto(0.08)')

#################################################   scrapingVeto + hcal laser veto  ################################
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                    applyfilter=cms.untracked.bool(True),
                                    debugOn=cms.untracked.bool(False),
                                    numtrack=cms.untracked.uint32(10),
                                    thresh=cms.untracked.double(0.25)
                                    )
if not isMC:
    process.dataFilter = cms.Sequence(process.hcallLaserEvent2012Filter+ process.scrapingVeto )
else:
    process.dataFilter = cms.Sequence()

#################################################   Good Primary Vertex ################################
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.selectPrimaryVertex = cms.EDFilter(
                                    "PrimaryVertexObjectFilter",
                                    filterParams=pvSelector.clone(minNdof=cms.double(4.0), maxZ=cms.double(24.0), maxRho=cms.double(2.0)),
                                    src=cms.InputTag('offlinePrimaryVertices')
                                    )

process.GoodVertexFilter = cms.EDFilter("VertexSelector",
                                           src = cms.InputTag("offlinePrimaryVertices"),
                                           cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
                                           filter = cms.bool(True)   
                                        )

#################################################   PAT APPLICATIONS   ##################################

#Removing MC Matching
from PhysicsTools.PatAlgos.tools.coreTools import *
#if not (isMC):
if (not isMC  and  not isEmbedded):
   removeMCMatching(process, ['All'])

##Switch to ak5PFJets (L2 and L3 Corrections are included)
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(process,
                    cms.InputTag('ak5PFJets'),
                    doJTA=True,
                    doBTagging=True,
                    jetCorrLabel=('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])),
                    doType1MET=False,
                    genJetCollection=cms.InputTag("ak5GenJets"),
                    doJetID=True,
                    jetIdLabel="ak5"
                    )
if not isMC:
    process.patJetCorrFactors.levels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual')

#$ Electron Id
from RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi import *
process.mvaID = cms.Sequence(process.mvaTrigV0 + process.mvaNonTrigV0 + process.mvaTrigNoIPV0)
#process.mvaID = cms.Sequence(process.mvaNonTrigV0 + process.mvaTrigNoIPV0)
process.ElectronIDs = cms.Sequence(process.simpleEleIdSequence)
process.patElectrons.electronIDSources = cms.PSet(
#    mvaTrigV0=cms.InputTag("mvaTrigV0"),
    mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0"),
    mvaTrigNoIPV0 = cms.InputTag("mvaTrigNoIPV0"),
    simpleEleId95relIso = cms.InputTag("simpleEleId95relIso"),
    simpleEleId90relIso = cms.InputTag("simpleEleId90relIso"),
    simpleEleId85relIso = cms.InputTag("simpleEleId85relIso"),
    simpleEleId80relIso = cms.InputTag("simpleEleId80relIso"),
    simpleEleId70relIso = cms.InputTag("simpleEleId70relIso"),
    simpleEleId60relIso = cms.InputTag("simpleEleId60relIso")
    )
#change PV source for electrons and muons
process.patElectrons.pvSrc = cms.InputTag("selectPrimaryVertex")
process.patMuons.pvSrc = cms.InputTag("selectPrimaryVertex")

from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

process.pfPileUp.Enable = cms.bool(True)

process.elPFIsoValueChargedAll04PFId.deposits[0].vetos= cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')
process.elPFIsoValueGamma04PFId.deposits[0].vetos= cms.vstring('EcalBarrel:ConeVeto(0.08)','EcalEndcaps:ConeVeto(0.08)')
process.elPFIsoValueChargedAll04NoPFId.deposits[0].vetos= cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')
process.elPFIsoValueGamma04NoPFId.deposits[0].vetos= cms.vstring('EcalBarrel:ConeVeto(0.08)','EcalEndcaps:ConeVeto(0.08)')

process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
#singleMu
process.muTriggerMu24 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                               src=cms.InputTag("skimmedPatMuons"),
                                                               matched=cms.InputTag("patTrigger"),
                                                               matchedCuts=cms.string('path("HLT_IsoMu24_eta2p1_v*")'),
                                                               maxDPtRel=cms.double(0.5),
                                                               maxDeltaR=cms.double(0.5),
                                                               resolveAmbiguities=cms.bool(True),
                                                               resolveByMatchQuality=cms.bool(True)
                                                                                                      )
#MuTau
process.tauTriggerMu18Tau25 = process.muTriggerMu24.clone(
                                                src=cms.InputTag('selectedPatTaus'),
                                                matchedCuts=cms.string('path("HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk1_eta2p1_v*") || path("HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v*")')
                                                )

process.muTriggerMu18Tau25 = process.muTriggerMu24.clone(
                                                src=cms.InputTag("skimmedPatMuons"),
                                                matchedCuts=cms.string('path("HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk1_eta2p1_v*") || path("HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v*")')
                                                )

process.tauTriggerMu17Tau20 = process.muTriggerMu24.clone(
                                                src=cms.InputTag("skimmedPatTausLT"),
                                                matchedCuts=cms.string('path("HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v*") || path("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v*")')
                                               )

process.muTriggerMu17Tau20 = process.muTriggerMu24.clone(
                                                src=cms.InputTag("skimmedPatMuons"),
                                                matchedCuts=cms.string('path("HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v*") || path("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v*")')
                                                )

#EleTau
process.tauTriggerEle20Tau20 = process.muTriggerMu24.clone(
                                                src=cms.InputTag("skimmedPatTausLT"),
                                                matchedCuts=cms.string('path("HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*") || path("HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_v*")')
                                                )

process.eleTriggerEle20Tau20 = process.muTriggerMu24.clone(
                                                src=cms.InputTag("skimmedPatElectrons"),
                                                matchedCuts=cms.string('path("HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*") || path("HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_v*")')
                                                )

#Double lepton
process.eleTriggerEle17Ele8 = process.muTriggerMu24.clone(
                                                src=cms.InputTag("skimmedPatElectrons"),
                                                matchedCuts=cms.string('path("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")')
                                                )

process.muTriggerMu17Mu8 = process.muTriggerMu24.clone(
                                                src=cms.InputTag("skimmedPatMuons"),
                                                matchedCuts=cms.string('path("HLT_Mu17_Mu8_v*") || path("HLT_Mu17_TkMu8_v*")')
                                                )

#EMu
process.eleTriggerEle17Mu8 = process.muTriggerMu24.clone(
                                                src=cms.InputTag("skimmedPatElectrons"),
                                                matchedCuts=cms.string('path("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")')
                                                )

process.eleTriggerEle8Mu17 = process.muTriggerMu24.clone(
                                                src=cms.InputTag("skimmedPatElectrons"),
                                                matchedCuts=cms.string('path("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")')
                                                )

process.muTriggerEle17Mu8 = process.muTriggerMu24.clone(
                                                src=cms.InputTag("skimmedPatMuons"),
                                                matchedCuts=cms.string('path("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")')
                                                )

process.muTriggerEle8Mu17 = process.muTriggerMu24.clone(
                                                src=cms.InputTag("skimmedPatMuons"),
                                                matchedCuts=cms.string('path("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")')
                                                )

#TauTau
process.tauTriggerDitau35 = process.muTriggerMu24.clone(
                                                src=cms.InputTag("skimmedPatTausTT"),
                                                matchedCuts=cms.string('path("HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1_v*") || path("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_v*")'),
                                                )

process.tauTriggerDitau30Jet30 = process.muTriggerMu24.clone(
                                                src=cms.InputTag("skimmedPatTausTT"),
                                                matchedCuts=cms.string('path("HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30_v*") || path("HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30_v*") || path("HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30_v*")')
                                                )

process.jetTriggerDitau30Jet30 = process.muTriggerMu24.clone(
                                                               src=cms.InputTag("selectedPatJets"),
                                                               matchedCuts=cms.string('path("HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30_v*") || path("HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30_v*") || path("HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30_v*")')
                                                                                                      )

#SoftMu

process.tauTriggerSoftMu = process.muTriggerMu24.clone(
                                                src=cms.InputTag("skimmedPatTausLT"),
                                                matchedCuts=cms.string('path("HLT_IsoMu8_eta2p1_LooseIsoPFTau20_L1ETM26_v*")')
                                                )

process.muTriggerSoftMu = process.muTriggerMu24.clone(
                                                src=cms.InputTag("skimmedPatMuons"),
                                                matchedCuts=cms.string('path("HLT_IsoMu8_eta2p1_LooseIsoPFTau20_L1ETM26_v*")')
                                                )

#SoftEle

process.tauTriggerSoftEle = process.muTriggerMu24.clone(
                                                src=cms.InputTag("skimmedPatTausLT"),
                                                matchedCuts=cms.string('path("HLT_Ele13_eta2p1_WP90Rho_LooseIsoPFTau20_L1ETM36_v*")')
                                                )

process.eleTriggerSoftEle = process.muTriggerMu24.clone(
                                                src=cms.InputTag("skimmedPatElectrons"),
                                                matchedCuts=cms.string('path("HLT_Ele13_eta2p1_WP90Rho_LooseIsoPFTau20_L1ETM36_v*")')
                                                )

#Vertexing
process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")

###--------------------------------------------------------------!!!!!!!!!!!!!!!!!!
#process.load("JetMETCorrections/Type1MET/pfMETsysShiftCorrections_cfi")
#from PhysicsTools.PatUtils.tools.metUncertaintyTools import runMEtUncertainties
#if isMC == False:
#    runMEtUncertainties(process,
#                        electronCollection = cms.InputTag('cleanPatElectrons'),
#                        photonCollection = '',
#                        muonCollection = 'selectedPatMuons',
#                        tauCollection = 'selectedPatTaus',
#                        jetCollection = cms.InputTag('selectedPatJets'),
#                        jetCorrLabel = 'L2L3Residual',
#                        doSmearJets = False,
#                        makeType1corrPFMEt = True,
#                        makeType1p2corrPFMEt = False,
#                        makePFMEtByMVA = False,
#                        makeNoPileUpPFMEt = False,
#                        doApplyType0corr = False,
#                        sysShiftCorrParameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data,
#                        doApplySysShiftCorr = False,
#                        addToPatDefaultSequence = False,
#                        )
#
#
#else:
#    runMEtUncertainties(process,
#                        electronCollection = cms.InputTag('cleanPatElectrons'),
#                        photonCollection = '',
#                        muonCollection = 'selectedPatMuons',
#                        tauCollection = 'selectedPatTaus',
#                        jetCollection = cms.InputTag('selectedPatJets'),
#                        jetCorrLabel = 'L3Absolute',
#                        doSmearJets = True,
#                        makeType1corrPFMEt = True,
#                        makeType1p2corrPFMEt = False,
#                        makePFMEtByMVA = False,
#                        makeNoPileUpPFMEt = False,
#                        doApplyType0corr = False,
#                        sysShiftCorrParameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc, 
#                        doApplySysShiftCorr = False,
#                        addToPatDefaultSequence = False,
#                        )
#
#process.patJetsNotOverlappingWithLeptonsForMEtUncertainty.checkOverlaps.taus.preselection = cms.string('pt > 20 && abs(eta) < 2.3 && leadPFChargedHadrCand.isNonnull() && leadPFChargedHadrCand.pt() > 5 && leadPFChargedHadrCand.mva_e_pi() < 0.6 && tauID("decayModeFinding") > 0.5 && tauID("againstMuonTight") > 0.5 && tauID("againstElectronLoose") > 0.5')
#process.patJetsNotOverlappingWithLeptonsForMEtUncertainty.checkOverlaps.muons.preselection = cms.string('pt > 25 && abs(eta) < 2.4')

#Making a cleanPatJet similar to smeraedJet
process.cleanPatJets.checkOverlaps.taus.src=cms.InputTag("selectedPatTaus")
process.cleanPatJets.checkOverlaps.taus.preselection= cms.string('pt > 20 && abs(eta) < 2.3 && leadPFChargedHadrCand.isNonnull() && leadPFChargedHadrCand.pt() > 5 && leadPFChargedHadrCand.mva_e_pi() < 0.6 && tauID("decayModeFinding") > 0.5 && tauID("againstMuonTight") > 0.5 && tauID("againstElectronLoose") > 0.5')
process.cleanPatJets.checkOverlaps.muons.src= cms.InputTag("selectedPatMuons")
process.cleanPatJets.checkOverlaps.muons.preselection = cms.string('pt > 25 && abs(eta) < 2.4')
process.cleanPatJets.checkOverlaps.tkIsoElectrons.requireNoOverlaps= cms.bool(True)
process.cleanPatJets.checkOverlaps.photons.requireNoOverlaps= cms.bool(True)

###--------------------------------------------------------------
## MVA MET configuration

from RecoMET.METPUSubtraction.mvaPFMET_leptons_cfi import isotaus 
process.isotaus.discriminators = cms.VPSet(
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),       selectionCut=cms.double(0.5)),
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits"), selectionCut=cms.double(0.5)),
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseElectronRejection"), selectionCut=cms.double(0.5)),
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection"),     selectionCut=cms.double(0.5))
            )
# One can overwrite the parameters in the config files for MVA MET
from RecoMET.METPUSubtraction.mvaPFMET_leptons_ditau_cfi import isotausditau
process.isotausditau.discriminators = cms.VPSet(
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),       selectionCut=cms.double(0.5)),
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr"),           selectionCut=cms.double(0.5)),
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseElectronRejection"), selectionCut=cms.double(0.5)),
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection"),     selectionCut=cms.double(0.5))
            )

process.load('JetMETCorrections.Configuration.JetCorrectionProducers_cff')
from RecoMET.METPUSubtraction.mvaPFMET_cff import calibratedAK5PFJetsForPFMEtMVA
process.load("RecoMET.METPUSubtraction.mvaPFMET_cff")
process.calibratedAK5PFJetsForPFMEtMVA = calibratedAK5PFJetsForPFMEtMVA.clone()
if isMC:
        process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3")
else:
        process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3Residual")
        
process.metRecoilCorPFMet = cms.EDProducer("MEtRecoilCorrectorProducer",
        genParticleTag = cms.InputTag("genParticles"),
        jetTag = cms.InputTag("selectedPatJets"),
        metTag = cms.InputTag("pfMet"),
        electronTag = cms.InputTag("isoelectrons"),
        muonTag = cms.InputTag("isomuons"),
        tauTag = cms.InputTag("isotaus"),
        inputFileNamezmm42X = cms.FileInPath("NtupleProducer2014/Analysis/data/RecoilCorrector_v7/recoilfits/recoilfit_zmm53X_2012_njet.root"),
        inputFileNamedatamm = cms.FileInPath("NtupleProducer2014/Analysis/data/RecoilCorrector_v7/recoilfits/recoilfit_datamm53X_2012_njet.root"),
        inputFileNamewjets = cms.FileInPath("NtupleProducer2014/Analysis/data/RecoilCorrector_v7/recoilfits/recoilfit_wjets53X_20pv_njet.root"),
        inputFileNamezjets = cms.FileInPath("NtupleProducer2014/Analysis/data/RecoilCorrector_v7/recoilfits/recoilfit_zmm53X_2012_njet.root"),
        inputFileNamehiggs = cms.FileInPath("NtupleProducer2014/Analysis/data/RecoilCorrector_v7/recoilfits/recoilfit_higgs53X_20pv_njet.root"),
        numOfSigmas = cms.double(1.0),
        minJetPt = cms.double(30.0),
        verbose = cms.bool(False),
        isMC = cms.bool(isMC),
        idxTau = cms.int32(-1),
    )

process.metRecoilCorMVAditau= process.metRecoilCorPFMet.clone(
metTag = cms.InputTag("pfMEtMVAditau")
)
process.metRecoilCorMVAmutau= process.metRecoilCorPFMet.clone(
metTag = cms.InputTag("pfMEtMVAmutau")
)
process.metRecoilCorMVAetau= process.metRecoilCorPFMet.clone(
metTag = cms.InputTag("pfMEtMVAetau")
)
process.metRecoilCorMVAemu= process.metRecoilCorPFMet.clone(
metTag = cms.InputTag("pfMEtMVAemu")
)

process.METRecoilCorrections = cms.Sequence(
                                         process.metRecoilCorPFMet
                                         * process.metRecoilCorMVAditau
                                         * process.metRecoilCorMVAmutau
                                         * process.metRecoilCorMVAetau
					 * process.metRecoilCorMVAemu
                                         )

process.skimmedPatElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("patElectrons"),
    cut = cms.string("pt > 9 && abs(eta) < 2.5")
)

process.skimmedPatMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("pt > 9 && abs(eta) < 2.5")
)

process.skimmedPatTausTT = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("patTaus"),
    cut = cms.string('pt > 40 && abs(eta) < 2.5 &&  tauID("decayModeFindingOldDMs") > 0.5 && (tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 4.0 || tauID("byVLooseIsolationMVA3oldDMwLT") > 0.5)')
)

process.skimmedPatTausLT = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("patTaus"),
    cut = cms.string('pt > 19 && abs(eta) < 2.5 &&  tauID("decayModeFindingOldDMs") > 0.5 && (tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 4.0 || tauID("byVLooseIsolationMVA3oldDMwLT") > 0.5)')
)

process.skimmedCorrPatTausLT = cms.EDProducer('ScaleTauProducer',
    tauSrc=cms.InputTag('skimmedPatTausLT'),
    NAME=cms.string('skimmedCorrPatTausLT')
)

process.skimmedCorrPatTausTT = cms.EDProducer('ScaleTauProducer',
    tauSrc=cms.InputTag('skimmedPatTausTT'),
    NAME=cms.string('skimmedCorrPatTausTT')
)

#skimLeptons=cms.Sequence()
skimLeptons=(process.skimmedPatElectrons * process.skimmedPatMuons * process.skimmedPatTausTT * process.skimmedPatTausLT * process.skimmedCorrPatTausLT * process.skimmedCorrPatTausTT)

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
  tModuleName = "selectedPatTausTT%i" % (tINDEX)
  tModule = cms.EDProducer('SinglePatTauProducer' ,
    tauSrc =cms.InputTag('skimmedCorrPatTausTT:skimmedCorrPatTausTT:NtupleProducer'),
    INDEX = cms.uint32(tINDEX),
    NAME=cms.string(tModuleName)
    )
  setattr(process, tModuleName, tModule)
  if (isTT):
     singlePatLeptons += tModule

for INDEX1 in range(9):
   for INDEX2 in range(INDEX1+1,10):
      tModuleName = "doublePatTausTT%ix%i" % (INDEX1,INDEX2)
      tModule = cms.EDProducer('DoublePatTauProducer' ,
        tauSrc1 =cms.InputTag("selectedPatTausTT%i:selectedPatTausTT%i:NtupleProducer" % (INDEX1,INDEX1)),
        tauSrc2 =cms.InputTag("selectedPatTausTT%i:selectedPatTausTT%i:NtupleProducer" % (INDEX2,INDEX2)),
        NAME=cms.string(tModuleName)
      )
      setattr(process, tModuleName, tModule)
      if (isTT):
         singlePatLeptons += tModule

for tINDEX in range(10):
  tModuleName = "selectedPatTausLT%i" % (tINDEX)
  tModule = cms.EDProducer('SinglePatTauProducer' ,
    tauSrc =cms.InputTag('skimmedCorrPatTausLT:skimmedCorrPatTausLT:NtupleProducer'),
    INDEX = cms.uint32(tINDEX),
    NAME=cms.string(tModuleName)
    )
  setattr(process, tModuleName, tModule)
  if (isMT or isET):
     singlePatLeptons += tModule

################################
# create the pair-wise mva mets
################################
 
pairWiseMvaMETs = cms.Sequence()

###########
# eTau METs
###########
 
for eINDEX in range(10):
  for tINDEX in range(10):
    metModuleName = "eTauMet%ix%i" % (eINDEX,tINDEX)
    eModuleName = "selectedPatElectrons%i:selectedPatElectrons%i:NtupleProducer" % (eINDEX,eINDEX)
    tModuleName = "selectedPatTausLT%i:selectedPatTausLT%i:NtupleProducer" % (tINDEX,tINDEX)
    metModule = process.pfMEtMVA.clone(
      corrector = cms.string('ak5PFL1FastL2L3'),
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
    tModuleName = "selectedPatTausLT%i:selectedPatTausLT%i:NtupleProducer" % (tINDEX,tINDEX)
    metModule = process.pfMEtMVA.clone(
      corrector = cms.string('ak5PFL1FastL2L3'),
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
    mModuleName = "selectedPatTausTT%i:selectedPatTausTT%i:NtupleProducer" % (mINDEX,mINDEX)
    tModuleName = "selectedPatTausTT%i:selectedPatTausTT%i:NtupleProducer" % (tINDEX,tINDEX)
    metModule = process.pfMEtMVA.clone(
      corrector = cms.string('ak5PFL1FastL2L3'),
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
    metModule = process.pfMEtMVA.clone(
      corrector = cms.string('ak5PFL1FastL2L3'),
      srcLeptons = cms.VInputTag(cms.InputTag(mModuleName),cms.InputTag(tModuleName)),
      )
    setattr(process, metModuleName, metModule)
    metModule.minNumLeptons = cms.int32(2)
    if (isEM):
       pairWiseMvaMETs += metModule

process.selectedPatJets.cut = cms.string('abs(eta) < 5.0')

##################################
# create the pair-wise recoil mets
##################################

pairWiseRecoil = cms.Sequence()

process.metRecoilCorPFMet = cms.EDProducer("MEtRecoilCorrectorProducer",
        genParticleTag = cms.InputTag("genParticles"),
        jetTag = cms.InputTag("selectedPatJets"),
        metTag = cms.InputTag("pfMet"),
        electronTag = cms.InputTag("isoelectrons"),
        muonTag = cms.InputTag("isomuons"),
        tauTag = cms.InputTag("isotaus"),
        inputFileNamezmm42X = cms.FileInPath("NtupleProducer2014/Analysis/data/RecoilCorrector_v7/recoilfits/recoilfit_zmm53XRR_2012_njet.root"),
        inputFileNamedatamm = cms.FileInPath("NtupleProducer2014/Analysis/data/RecoilCorrector_v7/recoilfits/recoilfit_datamm53XRR_2012_njet.root"),
        inputFileNamewjets = cms.FileInPath("NtupleProducer2014/Analysis/data/RecoilCorrector_v7/recoilfits/recoilfit_wjets53X_20pv_njet.root"),
        inputFileNamezjets = cms.FileInPath("NtupleProducer2014/Analysis/data/RecoilCorrector_v7/recoilfits/recoilfit_zmm53X_2012_njet.root"),
        inputFileNamehiggs = cms.FileInPath("NtupleProducer2014/Analysis/data/RecoilCorrector_v7/recoilfits/recoilfit_higgs53X_20pv_njet.root"),
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
    tModuleName = "selectedPatTausLT%i:selectedPatTausLT%i:NtupleProducer" % (tINDEX,tINDEX)
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
    tModuleName = "selectedPatTausLT%i:selectedPatTausLT%i:NtupleProducer" % (tINDEX,tINDEX)
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
    eModuleName = "selectedPatTausTT%i:selectedPatTausTT%i:NtupleProducer" % (eINDEX,eINDEX)
    tModuleName = "selectedPatTausTT%i:selectedPatTausTT%i:NtupleProducer" % (tINDEX,tINDEX)
    dModuleName = "doublePatTausTT%ix%i:doublePatTausTT%ix%i:NtupleProducer" % (eINDEX,tINDEX,eINDEX,tINDEX)
    recoilModule = process.metRecoilCorPFMet.clone(
        metTag = cms.InputTag("tauTauMet%ix%i" % (eINDEX,tINDEX)),
        electronTag = cms.InputTag(""),
        muonTag = cms.InputTag(""),
        tauTag = cms.InputTag(dModuleName)
	#verbose=cms.bool(True)
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


###----------------------------------------------------------------
## Calo MET
#
#process.produceCaloMEtNoHF = cms.Sequence(process.sumCaloTowersInEtaSlicesNoHF)
#
#if isMC:
#    process.metNoHFresidualCorrected.residualCorrLabel = cms.string("ak5CaloResidual")
#    process.metNoHFresidualCorrected.extraCorrFactor = cms.double(1.05)
#    process.metNoHFresidualCorrected.isMC = cms.bool(True)
#    process.produceCaloMEtNoHF += process.metNoHFresidualCorrected
#    process.metNoHFresidualCorrectedUp = process.metNoHFresidualCorrected.clone(
#        extraCorrFactor = cms.double(1.10)
#    )
#    process.produceCaloMEtNoHF += process.metNoHFresidualCorrectedUp
#    process.metNoHFresidualCorrectedDown = process.metNoHFresidualCorrected.clone(
#        extraCorrFactor = cms.double(1.0)
#    )
#    process.produceCaloMEtNoHF += process.metNoHFresidualCorrectedDown
#else:
#    process.metNoHFresidualCorrected.residualCorrLabel = cms.string("")
#    process.metNoHFresidualCorrected.extraCorrFactor = cms.double(1.0)
#    process.metNoHFresidualCorrected.isMC = cms.bool(False)
#    process.produceCaloMEtNoHF += process.metNoHFresidualCorrected
#
### Add caloJets matched to genTaus (for embedded only)
#process.load("LLRAnalysis.Utilities.genTauMatchedCaloJet_cff")
#
### Correct caloMetNoHF and L1ETM for calo-deposits around genTaus (for embedded only)
#process.load("LLRAnalysis.TauTauStudies.calibrateCaloMETandL1ETMforEmbedded_cff")

if isEmbedded:
    process.ak5JetTracksAssociatorAtVertex.tracks = cms.InputTag("tmfTracks")  #(for embedded only)
    process.jetTracksAssociatorAtVertex.tracks = cms.InputTag("tmfTracks")


makemet=cms.Sequence()
if isMC:
   makemet= singlePatLeptons * pairWiseMvaMETs * pairWiseRecoil
else:
   makemet= singlePatLeptons * pairWiseMvaMETs

#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck", 
#	ignoreTotal = cms.untracked.int32(1) # default is one 
#)

#################################################   PATH of CONFIG FILE ##################################
if isEmbedded:
   process.p = cms.Path (
                      process.dataFilter *
                      process.selectPrimaryVertex *
                      (process.jetTracksAssociatorAtVertex+process.ak5JetTracksAssociatorAtVertex +process.btagging )*
                      process.PFTau * 
                      process.recoTauClassicHPSSequence *
                      process.mvaID   *
                      process.ElectronIDs  *
                      process.inclusiveVertexing *
                      process.pfParticleSelectionSequence +
                      process.type0PFMEtCorrection + 
                      process.patPFMETtype0Corr + 
                      process.patDefaultSequence +
                      #skimLeptons*
                      process.eleIsoSequence*
                      process.pfMuonIsolationSequence                     *
                      #process.metUncertaintySequence*
                      process.producePatPFMETCorrections *
                      process.puJetIdSqeuence *
                      process.pfMEtMVAsequence *
                      process.pfMEtMVAsequenceditau *
                      process.pfMEtMVAsequencemutau *
                      process.pfMEtMVAsequenceetau *
                      process.pfMEtMVAsequenceemu *
                      process.METRecoilCorrections *                     
		      makemet*
                      #process.makeTauMatchedCaloJets+process.calibrateCaloMETandL1ETMforEmbedded+
                      #process.produceCaloMEtNoHF *
                      process.myanalysis 
                      )
else:
   process.p = cms.Path (
                      process.dataFilter *
                      process.selectPrimaryVertex *
                      process.PFTau *
                      process.recoTauClassicHPSSequence *
                      process.mvaID   *
                      process.ElectronIDs  *
                      process.inclusiveVertexing *
                      process.pfParticleSelectionSequence +
                      process.type0PFMEtCorrection + 
                      process.patPFMETtype0Corr + 
                      process.patDefaultSequence +
#                      skimLeptons*
                      process.pfMuonIsolationSequence*
                      process.eleIsoSequence*
                      #process.metUncertaintySequence*
                      process.producePatPFMETCorrections *
                      process.puJetIdSqeuence *
                      process.pfMEtMVAsequence *
                      process.pfMEtMVAsequenceditau *
                      process.pfMEtMVAsequencemutau *
                      process.pfMEtMVAsequenceetau *
                      process.pfMEtMVAsequenceemu *
                      process.METRecoilCorrections *
		      makemet*
                      process.myanalysis 
                      )

#################################################   NEEDED FOR TRIGGER MATCHING   #######################

process.patDefaultSequence*=process.skimmedPatTausLT
process.patDefaultSequence*=process.skimmedPatTausTT
process.patDefaultSequence*=process.skimmedPatMuons
process.patDefaultSequence*=process.skimmedPatElectrons
process.patDefaultSequence*=process.skimmedCorrPatTausLT
process.patDefaultSequence*=process.skimmedCorrPatTausTT

#process.seq=cms.Sequence(process.patDefaultSequence+skimLeptons)
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger(process)
removeCleaningFromTriggerMatching(process)
#process.patTriggerEvent.patTriggerMatches = []
#switchOnTrigger(process)
switchOnTriggerMatching(process, triggerMatchers=['eleTriggerEle17Ele8','muTriggerMu17Mu8','muTriggerEle17Mu8','eleTriggerEle17Mu8','muTriggerEle8Mu17','eleTriggerEle8Mu17','tauTriggerDitau35','tauTriggerDitau30Jet30','jetTriggerDitau30Jet30','eleTriggerEle20Tau20','tauTriggerEle20Tau20','muTriggerMu17Tau20','tauTriggerMu17Tau20','muTriggerMu18Tau25','tauTriggerMu18Tau25','muTriggerMu24','tauTriggerSoftMu','tauTriggerSoftEle','muTriggerSoftMu','eleTriggerSoftEle'])
#process.patTriggerEvent.processName = "HLT"
#process.patTrigger.processName = "HLT"
