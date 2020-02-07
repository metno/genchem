#!/usr/bin/env python3
Usage="""

  reads  boxChem output files such as AeroEmChem19a.csv
  and extracts concentrations for one pollutant.

  Usage:

     getConcs OH  AeroEmChem19a.csv

"""
import pandas as pd
import sys

if len(sys.argv) == 1:
  poll='O3'
  ifile='boxEmChem19a.csv'
else:
  poll  = sys.argv[1]
  ifile = sys.argv[2]

x=pd.read_csv(ifile)
concs=x[poll].values   # 1st cell will be unit, e.g. ppb
unit= concs[0]
vals = [ float(v) for v in concs[1:] ]

ofile='ResConcs_' + ifile.replace(".csv","_%s_%s.txt"%(poll, unit))
print(ofile)
print('%s MAXMIN %s %s: %f %f' % ( ifile, poll, unit, min(vals), max(vals) ))
with open(ofile,'w') as f:
  f.write( '%s\n' %poll )
  for v in vals:
    f.write('%f\n' % v )


