#!/usr/bin/env python3
import argparse
import pandas as pd
import sys
Usage="""

  reads  boxChem output files such as AeroEmChem19a.csv
  and extracts concentrations for one pollutant.

  Usage:

     getConcs -v OH  -i AeroEmChem19a.csv

"""
#----------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-i', required=True,
                     help='input file, e.g. -i Ouput/test1.csv')
parser.add_argument('-v',required=True,
    help='species wanted , e.g. -v OH')
parser.print_usage = parser.print_help
args = parser.parse_args()

print(args, len(args.i))
#----------------------------------------------------------
poll= args.v
ifile= args.i


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


