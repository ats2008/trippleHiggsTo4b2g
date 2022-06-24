import FWCore.ParameterSet.Config as cms
process = cms.Process("Validation")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15')
process.load("PhysicsTools.PatAlgos.slimming.slimmedAddPileupInfo_cfi")
process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
"file:/afs/cern.ch/work/m/mstamenk/public/HHH/HHH4b2y_RunIISummer20UL17/RunIISummer20UL17MINIAODSIM_71.root",
"file:/afs/cern.ch/work/m/mstamenk/public/HHH/HHH4b2y_RunIISummer20UL17/RunIISummer20UL17MINIAODSIM_73.root",
"file:/afs/cern.ch/work/m/mstamenk/public/HHH/HHH4b2y_RunIISummer20UL17/RunIISummer20UL17MINIAODSIM_70.root",
                    ),
                duplicateCheckMode=cms.untracked.string("noDuplicateCheck"),

                            )
process.bRegression = cms.EDAnalyzer(
    "genNtuple",
    label_pileUp = cms.untracked.InputTag("slimmedAddPileupInfo"),
    #GenRunInf  = cms.InputTag("generator", "", "GEN"),
    #GenInf     = cms.InputTag("generator", "", "GEN"),
    pruned     = cms.InputTag("prunedGenParticles"),
    photons = cms.InputTag("slimmedPhotons"),
    recJet = cms.InputTag("slimmedJets"),
    genJet = cms.InputTag("slimmedGenJets"),
    #genParticles = cms.InputTag('genParticles'),
    )
process.TFileService = cms.Service("TFileService", fileName = cms.string('pu_trees_VH.root'))
process.p = cms.Path(
     process.bRegression
    )
process.schedule = cms.Schedule(process.p)

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))
