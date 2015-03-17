import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:myfile.root'
	'file:/afs/cern.ch/work/c/ccaillol/DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_PU_S7_START52_V9_v2.root'
    )
)

process.myProducerLabel = cms.EDProducer('SinglePatElectronProducer'
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)


process.p = cms.Path(process.myProducerLabel)

process.e = cms.EndPath(process.out)
