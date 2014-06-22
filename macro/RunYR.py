#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''

# Instructions: macro root parameters
# void BosonRPTotem(string inputfile, string outputfile,double XSmcW, double XSmcZ, double lumi)

#
# DY to Z, plus
#

# Muon
os.system("root -l BosonRPTotem.C\'(\"/storage1/dmf/YellowReport/CMSSW_6_2_5/src/DiffractiveAnalysisYR/DiffractiveBoson/test/samples/ZtoMuMu_plus/ZtoMuMu_plus.root\",\"histo_ZtoMuMuPlus.root\",1200.0,300.0,100.0)\'")

# Electron
os.system("root -l BosonRPTotem.C\'(\"/storage1/dmf/YellowReport/CMSSW_6_2_5/src/DiffractiveAnalysisYR/DiffractiveBoson/test/samples/ZtoEE_plus/ZtoEE_plus.root\",\"histo_ZtoEEPlus.root\",1200.0,150.0,100.0)\'")

#
# DY to Z, minus
#

# Muon
os.system("root -l BosonRPTotem.C\'(\"/storage1/dmf/YellowReport/CMSSW_6_2_5/src/DiffractiveAnalysisYR/DiffractiveBoson/test/samples/ZtoMuMu_minus/ZtoMuMu_minus.root\",\"histo_ZtoMuMuMinus.root\",1200.0,300.0,100.0)\'")

# Electron
os.system("root -l BosonRPTotem.C\'(\"/storage1/dmf/YellowReport/CMSSW_6_2_5/src/DiffractiveAnalysisYR/DiffractiveBoson/test/samples/ZtoEE_minus/ZtoEE_minus.root\",\"histo_ZtoEEMinus.root\",1200.0,150.0,100.0)\'")


#
# W, plus
#

# Muon
os.system("root -l BosonRPTotem.C\'(\"/storage1/dmf/YellowReport/CMSSW_6_2_5/src/DiffractiveAnalysisYR/DiffractiveBoson/test/samples/WtoMuNu_plus/WtoMuNu_plus.root\",\"histo_WtoMuNuPlus.root\",1200.0,300.0,100.0)\'")

# Electron
os.system("root -l BosonRPTotem.C\'(\"/storage1/dmf/YellowReport/CMSSW_6_2_5/src/DiffractiveAnalysisYR/DiffractiveBoson/test/samples/WtoENu_plus/WtoENu_plus.root\",\"histo_WtoENuPlus.root\",1200.0,150.0,100.0)\'")

#
# W, minus
#

# Muon
os.system("root -l BosonRPTotem.C\'(\"/storage1/dmf/YellowReport/CMSSW_6_2_5/src/DiffractiveAnalysisYR/DiffractiveBoson/test/samples/WtoMuNu_minus/WtoMuNu_minus.root\",\"histo_WtoMuNuMinus.root\",1200.0,300.0,100.0)\'")

# Electron
os.system("root -l BosonRPTotem.C\'(\"/storage1/dmf/YellowReport/CMSSW_6_2_5/src/DiffractiveAnalysisYR/DiffractiveBoson/test/samples/WtoENu_minus/WtoENu_minus.root\",\"histo_WtoENuMinus.root\",1200.0,150.0,100.0)\'")

print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
