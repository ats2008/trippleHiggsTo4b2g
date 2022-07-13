import FWCore.ParameterSet.Config as cms

Ntuple = cms.EDAnalyzer(
        "genNtuple",
        #GenRunInf  = cms.InputTag("generator", "", "GEN"),
        #GenInf     = cms.InputTag("generator", "", "GEN"),
        #genParticles = cms.InputTag('genParticles'),
        pruned     = cms.InputTag("prunedGenParticles"),
        label_pileUp = cms.untracked.InputTag("slimmedAddPileupInfo"),
        photons = cms.InputTag("slimmedPhotons"),
        doGenJets = cms.bool(False),
        doRecoJets = cms.bool(False),
        doRecoPhotons = cms.bool(False),
        recJet = cms.InputTag("slimmedJets"),
        genJet = cms.InputTag("slimmedGenJets"),
    )

TFileService = cms.Service("TFileService",
                      fileName = cms.string('ntuple.root')
                      )

def switchOnRecoPhotons(process,switch=True):
    process.Ntuple.doRecoPhotons = switch
def switchOnRecoJets(process,switch=True):
    process.Ntuple.doRecoJets = switch
def switchOnGenJets(process,switch=True):
    process.Ntuple.doGenJets = switch

def getDefaultProcess(maxEvents=-1):
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
        
    process.TFileService = TFileService
    process.Ntuple = Ntuple
    process.p = cms.Path( process.Ntuple )
    
    process.schedule = cms.Schedule(process.p)
    process.MessageLogger.cerr.FwkReport.reportEvery = 100
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(maxEvents))
    
    
    return process

def addBjetRegression(process):
    
    ########## b Jet Energy Corrections ################
    # Modules from NanoAOD
    from PhysicsTools.NanoAOD.jets_cff import bJetVars,bjetNN
    #### update modules parameters
    # create variables for regression
    
    process.bJetVars = bJetVars.clone(
        src = cms.InputTag("slimmedJets"),
    )
    
    # add variables to the jet collection
    # (using only bJetVars, which are variables relevant to regression)
    process.updatedJetsWithUserData = cms.EDProducer("PATJetUserDataEmbedder",
        src = cms.InputTag("slimmedJets"),
        userFloats = cms.PSet(
             leadTrackPt = cms.InputTag("bJetVars:leadTrackPt"),
             leptonPtRel = cms.InputTag("bJetVars:leptonPtRel"),
             leptonPtRatio = cms.InputTag("bJetVars:leptonPtRatio"),
             leptonPtRelInv = cms.InputTag("bJetVars:leptonPtRelInv"),
             leptonPtRelv0 = cms.InputTag("bJetVars:leptonPtRelv0"),
             leptonPtRatiov0 = cms.InputTag("bJetVars:leptonPtRatiov0"),
             leptonPtRelInvv0 = cms.InputTag("bJetVars:leptonPtRelInvv0"),
             leptonDeltaR = cms.InputTag("bJetVars:leptonDeltaR"),
             leptonPt = cms.InputTag("bJetVars:leptonPt"),
             vtxPt = cms.InputTag("bJetVars:vtxPt"),
             vtxMass = cms.InputTag("bJetVars:vtxMass"),
             vtx3dL = cms.InputTag("bJetVars:vtx3dL"),
             vtx3deL = cms.InputTag("bJetVars:vtx3deL"),
             ptD = cms.InputTag("bJetVars:ptD"),
             genPtwNu = cms.InputTag("bJetVars:genPtwNu"),
    #        pileupJetIdDiscriminant = cms.InputTag("pileupJetIdUpdated:fullDiscriminant"),
             ),
        userInts = cms.PSet(
            vtxNtrk = cms.InputTag("bJetVars:vtxNtrk"),
            leptonPdgId = cms.InputTag("bJetVars:leptonPdgId"),
    #       pileupJetId = cms.InputTag("pileupJetIdUpdated:fullId"),
        ),
        )
    # obtain the b jet regression corrections
    process.bjetNN = bjetNN.clone(
        src = cms.InputTag("updatedJetsWithUserData"),
    )
    
    # add the b jet regression corrections to the jet collection
    process.slimmedJetsbReg = cms.EDProducer("PATJetUserDataEmbedder",
         src = cms.InputTag("updatedJetsWithUserData"),
         userFloats = cms.PSet(
             bJetRegCorr = cms.InputTag("bjetNN:corr"),
             bJetRegRes = cms.InputTag("bjetNN:res"),
    #         pileupJetIdDiscriminant = cms.InputTag("pileupJetIdUpdated:fullDiscriminant"),
             ),
    #     userInts = cms.PSet(
    #         pileupJetId = cms.InputTag("pileupJetIdUpdated:fullId"),
    #         ),
    )
    
    #process.BJetRegression = cms.Task()
    #process.BJetRegression.add(process.bJetVars)
    #process.BJetRegression.add(process.updatedJetsWithUserData)
    #process.BJetRegression.add(process.bjetNN)
    #process.BJetRegression.add(process.slimmedJetsbReg)
    
    process.bJetRegressionSequence=cms.Sequence(     process.bJetVars* process.updatedJetsWithUserData*  process.bjetNN*  process.slimmedJetsbReg )
    process.Ntuple.recJet = cms.InputTag("slimmedJetsbReg")
    process.p.insert(0,process.bJetRegressionSequence)

    return process

def getMCConfig(maxEvents=-1):
    process=getDefaultProcess(maxEvents)
    switchOnGenJets(process,True)
    switchOnRecoJets(process,True)
    switchOnRecoPhotons(process,True)
    return process

def getGenOnlyConfig(maxEvents):
    process=getDefaultProcess(maxEvents)
    switchOnGenJets(process,True)
    switchOnRecoJets(process,False)
    switchOnRecoPhotons(process,False)
    return process

