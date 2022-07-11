import FWCore.ParameterSet.Config as cms

import trippleHiggsTo4b2g.GenValidation.customize_ntuplizer as customizer

process = customizer.getGenOnlyConfig(100)
customizer.addBjetRegression(process)

process.source.fileNames = cms.untracked.vstring(
            "file:/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/fastSimStudy/CMSSW_10_2_22/eventGen/MiniAODSIM/condor/results/mc/v1/c3_0_c4_0_HHHto4b2gamma/c3_0_c4_0_0.root",
            )
