import FWCore.ParameterSet.Config as cms

process = cms.Process("NtupleProducer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring(
'/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root',
)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
#process.MessageLogger.cerr.FwkReport.reportEvery=cms.untracked.int32(10)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

isMC=True
isEmbedded=False
isTT=True
isET=False
isMT=False
isEM=False
HLTProcessName="HLT"

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
                                    Is_MC=cms.bool(isMC),
                                    Is_Embedded=cms.bool(isEmbedded),
                                    Is_TT=cms.bool(isTT),
                                    Is_EM=cms.bool(isEM),
                                    Is_MT=cms.bool(isMT),
                                    Is_ET=cms.bool(isET),

)


process.p = cms.Path(process.ntupler)
