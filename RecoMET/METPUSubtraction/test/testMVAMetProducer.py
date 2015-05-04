import FWCore.ParameterSet.Config as cms

process = cms.Process("analysis")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10)
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionProducers_cff')
process.load('RecoMET.METPUSubtraction.mvaPFMET_cff')

# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00CC714A-F86B-E411-B99A-0025904B5FB8.root'
        ),
                            skipEvents = cms.untracked.uint32(0)                        
)

process.output = cms.OutputModule("PoolOutputModule",
                                  outputCommands = cms.untracked.vstring('keep *'),
                                  fileName = cms.untracked.string("MVaTest.root")
)       

process.pfMVAMEt.inputFileNames = cms.PSet(
    U     = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmet_7_2_X_MINIAOD_BX50PU40_Jan2015.root'),
    DPhi  = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_2_X_MINIAOD_BX50PU40_Jan2015.root'),
    CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_2_X_MINIAOD_BX50PU40_Jan2015.root'),
    CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_2_X_MINIAOD_BX50PU40_Jan2015.root')
    )

process.ana      = cms.Sequence(process.pfMVAMEtSequence)
process.p        = cms.Path(process.ana)
process.outpath  = cms.EndPath(process.output)
