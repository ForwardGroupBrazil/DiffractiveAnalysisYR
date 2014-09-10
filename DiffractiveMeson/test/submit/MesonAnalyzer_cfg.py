#flake8: noqa

'''
>>--------------------<<
MESON Jpsi NTuple Producer
>>--------------------<<

Goal:
Produce MesonJpsi ntuple.

Usage:
cmsRun MesonAnalyzer_cfg.py

Example:
cmsRun MesonAnalyzer_cfg.py Run=mesonJpsi

Authors: Brazilian/PPS Group
'''

# Loading Modules and Python Libraries.
import FWCore.ParameterSet.Config as cms
import os, sys
import atexit
from FWCore.ParameterSet.VarParsing import VarParsing

# Setting Input Parameters from Line Command.
options = VarParsing ('analysis')
options.register('Run','mesonJpsi',VarParsing.multiplicity.singleton, VarParsing.varType.string,"Option to Run:  mesonJpsi")
options.parseArguments()

# Some Variables.
mesonJpsi = False

class config: pass
config.debug = False # Debugger
config.Ebeam = 4000. # Beam Energy


# Setting Code Options.
print("")
print("#########################")
print("Running NTuple ProducerMC")
print("#########################")
print("")

if options.Run == "mesonJpsi":
  print("")
  print(">> Running: mesonJpsi")
  print("")
  mesonJpsi = True
  fileout = 'out.root'


else:
  print("")
  print("")
  raise RuntimeError, "Unknown option. EXIT! YOU NEED TO SETUP WITH ONE OF THE CORRECT OPTIONS."
  print("")

# CMSSW Default Configuration.
process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Openning files mesonJpsi.
if mesonJpsi:
   process.source = cms.Source("PoolSource",
      duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
      fileNames = cms.untracked.vstring(
          'file:JpsiToMuMuSDMinus.root'
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
process.mesonAnalyzer = cms.EDAnalyzer('MesonAnalyzer',
			       ParticleFlowTag = cms.InputTag("particleFlow"),
                               VertexTag = cms.InputTag("goodOfflinePrimaryVertices"),
                               pTPFThresholdCharged = cms.double(0.1),
                               energyPFThresholdBar = cms.double(1.5),
                               energyPFThresholdEnd = cms.double(3.5),
                               energyPFThresholdHF = cms.double(4.0),
                               GenParticleTag = cms.InputTag("genParticles"),
                               EBeam = cms.double(config.Ebeam),
                               debug = cms.bool(config.debug),
                               muonTag = cms.InputTag("muons")
)

# TFileService.
process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string(fileout)
                                   )


process.p = cms.Path(process.goodOfflinePrimaryVertices*process.mesonAnalyzer)

