#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''

# Instructions: macro root parameters
# void BosonRPTotem(string inputfile, string outputfile,double XSmcJpsi, double lumi)

#
#  plus
#

# Muon
os.system("root -l -q MesonRPTotem.C\'(\"/afs/cern.ch/work/e/eliza/private/TOTEM_JPSI/MonteCarlo_13TeV/CMSSW_6_2_5/src/DiffractiveAnalysisYR/DiffractiveMeson/test/samples/JpsitoMuMu_plus/JpsitoMuMu_plus.root\",\"histo_JpsitoMuMuPlus.root\",1.7e+4,100.0)\'")

os.system("root -l -q MesonRPTotem.C\'(\"/afs/cern.ch/work/e/eliza/private/TOTEM_JPSI/MonteCarlo_13TeV/CMSSW_6_2_5/src/DiffractiveAnalysisYR/DiffractiveMeson/test/samples/JpsitoMuMu_plus/JpsitoMuMu_plus.root\",\"histo_JpsitoMuMuPlus10pb.root\",1.7e+4,10.0)\'")  

#
# minus
#

# Muon
os.system("root -l -q MesonRPTotem.C\'(\"/afs/cern.ch/work/e/eliza/private/TOTEM_JPSI/MonteCarlo_13TeV/CMSSW_6_2_5/src/DiffractiveAnalysisYR/DiffractiveMeson/test/samples/JpsitoMuMu_minus/JpsitoMuMu_minus.root\",\"histo_JpsitoMuMuMinus10pb.root\",1.5e+4,10.0)\'")


print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
