#flake8: noqa

'''
>>--------------------<<
CEP NTuple Producer
>>--------------------<<

Goal:
Produce CEP ntuple.

Usage:
cmsRun DijetsAnalyzer_cfg.py

Example:
cmsRun DijetsAnalyzer_cfg.py

Authors: Brazilian/PPS Group
'''

# Loading Modules and Python Libraries.
import FWCore.ParameterSet.Config as cms
import os, sys
import atexit

class config: pass
config.debug = False # Debugger
config.Ebeam = 6500. # Beam Energy
config.jettag = "ak5PFJets" #ak5PFJets or ak5PFJetsCHS

# Setting Code Options.
print("")
print("#########################")
print("Running NTuple ProducerMC")
print("#########################")
print("")

print("")
print(">> Running: CEP")
print("")
fileout = 'out.root'

# CMSSW Default Configuration.
process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Openning file
process.source = cms.Source("PoolSource",
      duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
      fileNames = cms.untracked.vstring(
          'file:/storage5/dmf/TestSamples/YR/CEP.root'
     )
   )

# Check Inputs
print("")
print(">>> Input Options:")
print("Energy beam: %.2f" % config.Ebeam)
print("Debug: %s" % config.debug)
print("Output file: %s" % fileout)
print("")

# Adding Good Primary Vertex.
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )


# EDAnalyzer Parameters.
process.DijetsAnalyzer = cms.EDAnalyzer('DijetsAnalyzer',
			       ParticleFlowTag = cms.InputTag("particleFlow"),
                               VertexTag = cms.InputTag("goodOfflinePrimaryVertices"),
                               pTPFThresholdCharged = cms.double(0.1),
                               energyPFThresholdBar = cms.double(1.5),
                               energyPFThresholdEnd = cms.double(3.5),
                               energyPFThresholdHF = cms.double(4.0),
                               GenParticleTag = cms.InputTag("genParticles"),
                               EBeam = cms.double(config.Ebeam),
                               debug = cms.bool(config.debug),
                               JetTag = cms.InputTag(config.jettag)
)

# TFileService.
process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string(fileout)
                                   )


process.p = cms.Path(process.goodOfflinePrimaryVertices*process.DijetsAnalyzer)

