#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''

# Instructions: macro root parameters
# void CEPRPTotem(string inputfile, string outputfile,double XSmcW, double XSmcZ, double lumi)

#
# EXHUME
#

os.system("root -l -q CEPRPTotem.C\'(\"out.root\",\"histo_CEP.root\",20.0,10.0)\'")

print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
