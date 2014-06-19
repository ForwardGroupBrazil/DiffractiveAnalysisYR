#flake8: noqa

'''
>>--------------------<<
BOSON Z/W NTuple Producer
>>--------------------<<

Goal:
Produce BosonZ ntuple.

Usage:
cmsRun BosonAnalyzer_cfg.py

Example:
cmsRun BosonAnalyzer_cfg.py Run=bosonZ

Optional arguments:
Run = bosonW or bosonZ

Authors: Brazilian/PPS Group
'''

# Loading Modules and Python Libraries.
import FWCore.ParameterSet.Config as cms
import os, sys
import atexit
from FWCore.ParameterSet.VarParsing import VarParsing

# Setting Input Parameters from Line Command.
options = VarParsing ('analysis')
options.register('Run','bosonZ',VarParsing.multiplicity.singleton, VarParsing.varType.string,"Option to Run: bosonW or bosonZ")
options.parseArguments()

# Some Variables.
bosonZ = False
bosonW = False

class config: pass
config.debug = False # Debugger
config.Ebeam = 6500. # Beam Energy


# Setting Code Options.
print("")
print("#########################")
print("Running NTuple ProducerMC")
print("#########################")
print("")

if options.Run == "bosonZ":
  print("")
  print(">> Running: BosonZ")
  print("")
  bosonZ = True
  fileout = 'out.root'

elif options.Run == "bosonW":
  print("")
  print(">> Running: BosonW")
  print("")
  bosonW = True
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

# Openning files BosonZ or BosonW.
if bosonZ:
   process.source = cms.Source("PoolSource",
      duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
      fileNames = cms.untracked.vstring(
          'file:/storage1/dmf/TestSamples/YR/bosonzRECO.root'
     )
   )

if bosonW:
   process.source = cms.Source("PoolSource",
      duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
      fileNames = cms.untracked.vstring(
          'file:/storage1/dmf/TestSamples/YR/bosonwRECO.root'
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
process.BosonAnalyzer = cms.EDAnalyzer('BosonAnalyzer',
			       ParticleFlowTag = cms.InputTag("particleFlow"),
                               VertexTag = cms.InputTag("goodOfflinePrimaryVertices"),
                               pTPFThresholdCharged = cms.double(0.1),
                               energyPFThresholdBar = cms.double(1.5),
                               energyPFThresholdEnd = cms.double(3.5),
                               energyPFThresholdHF = cms.double(4.0),
                               GenParticleTag = cms.InputTag("genParticles"),
                               EBeam = cms.double(config.Ebeam),
                               debug = cms.bool(config.debug),
                               electronTag = cms.InputTag("gsfElectrons"),
                               muonTag = cms.InputTag("muons"),
                               metTag = cms.InputTag("pfMet")
)

# TFileService.
process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string(fileout)
                                   )


process.p = cms.Path(process.goodOfflinePrimaryVertices*process.BosonAnalyzer)

