#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

Usage = '''
  From Output, type e.g.

    ../scripts/boxplots.py  -v OH -i Output/EmChemXX.csv  [--png] [--dbg]
    ../scripts/boxplots.py  -v OH NO2  -i Output/EmChemXX.csv  [--png] [--dbg]

  Can have two or more inputs, e.g.

    ../scripts/boxplots.py  -v O3 -i Output/EmChem09.csv Output/MCM_v3.3.csv  --png

   'ALL' or 'DEF' is special either ALL or a default list, e.g.

    ../scripts/boxplots.py  -v ALL  -i Output/EmChemXX.csv  [--png]

    ALL or DEF also trigger --png

   --png  saves plots as png (not on screen)

   --dbg  extra debug info

'''

def rdboxdata(fname):
  """ Extracts data from BoxChem output.
      Scans 2 lines for headers and units, then grabs data
  """
  f=open(fname,'r')
  headers = f.readline().rstrip().split(',')
  units   = f.readline().rstrip().split(',')
  f.close()
  data=np.genfromtxt(fname,delimiter=',',skip_header=2)
  return headers, units, data

def find_index(str,strlist):
  for n in range(len(strlist)):
    if str == strlist[n]: return n 
  return -1

#----------------------------------------------------------
# Some colours & a nice set of variables for comparisons:
# ALSO if needed:
#http://stackoverflow.com/questions/4805048/how-to-get-different-colored-lines-for-different-plots-in-a-single-figure
cols = ['blue', 'green', 'orange', 'magenta', 'black', 'red', # 0-5
    'greenyellow', 'cornflowerblue', 'darkviolet', 'forestgreen']
def_variabs = ['O3', 'OH', 'HO2', 'CO', 'NO2', 'NO3', 'NO', 'OD', 'CH3O2', 
    'H2O2', 'HNO3', 'C2H4', 'C5H8', 'OP', 'HCHO', 'HONO', 'PAN', 'PANs_SUM', 'RO2s_SUM']
linestyle = '- -. -- :'.split()
#----------------------------------------------------------

#parser = argparse.ArgumentParser(usage=Usage)
parser = argparse.ArgumentParser()
parser.add_argument('-i',nargs='+', required=True,
                     help='input file(s), e.g. -i Ouput/test1.csv')
parser.add_argument('-v','--vars',nargs='+',required=True,
    help='species wanted , e.g. -v OH or -v O3 NO2 or -v ALL or DEF')
parser.add_argument('-x','--outlabel', 
                     help='base label for plot output')
parser.add_argument('-p','--png', 
                     help='produce png file(s)',action='store_true')
parser.add_argument('-s','--statsname', 
                     help='name for stats output')
parser.add_argument('-t','--txt', 
                     help='produce txt file(s)',action='store_true')
parser.add_argument('-u','--unitless', 
                     help='no units in plotfile name',action='store_true')
parser.add_argument('--dbg', help='extra debug info',action='store_true')
args = parser.parse_args()

print(args, len(args.i))
#----------------------------------------------------------

# First, collect all data, headers, etc. (confusing mix
# of dict and list. Needs clean.

vars={}; h={}; u={}; d={}; flabel=[]
several_files = len(args.i) > 1

modelcol=dict()
for nf, fname in enumerate(args.i):
  if args.dbg: print('FILE ', fname) # , args.i)
  h[nf],u[nf],d[nf] = rdboxdata( fname )
  #if several_files:
  flabel.append( os.path.basename(fname).replace('.csv','') )

  if args.vars[0] == 'ALL':
    vars[nf] =  h[nf][1:]
  elif args.vars[0] == 'DEF':
    vars[nf] = def_variabs
  else:
    vars[nf] = args.vars

  print(nf, 'VARS', flabel[nf], fname,vars[nf][0], vars[nf][-1] )

# ------ pre-screen values to get global max across files
gmax = dict()   # max value across all nf files
gunits = dict()
for v in vars[0]: # O3, OH, etc
  gmax[v] = 0.0
  gunits[v] = ''
  for n in range(nf+1):
    ind  = find_index(v,h[n]) # returns -1 if not found
    if ind>0:
      vals=d[n][:,ind]
      gmax[v] = max( gmax[v], max(vals) )
      gunits[v] = u[n][ind]
  print('Global max ', v, gmax[v], gunits[v] )

# ---- start label for output files:
olabel = '_'.join(flabel) # + '.txt'
if args.outlabel: olabel = args.outlabel
if args.dbg:print('olabel: '+olabel)
assert len(olabel)<100, \
  "\n\nERROR: olabel too long!!! Use --outlabel option\n\n"+olabel

# ---------- to output mean, min, max:
statsname = 'Stats_' + olabel + '.txt'
if args.statsname: statsname=args.statsname
if args.dbg:print('statsname: '+statsname)

stats = open(statsname,'w')
stats.write('%-20s %12s %12s %12s %12s %12s\n' % ('Label',
            'v','unit','mean','min','max' ) )

# ---------- loop over species (from 1st file):
for v in vars[0]:
   plot_ok = gmax[v] > 0.0   # True
   if not plot_ok: break # See below
   unit=gunits[v]
   stats.write('\n')
   nls=0; nco=0
   for n in range(nf+1):

      #if not plot_ok: break # See below

      ind  = find_index(v,h[n]) # returns -1 if not found
      maxv = -999.0
      nls += 1; nco += 1

      if ind > 0:
        vals=d[n][:,ind]
        maxv = np.max(vals)   # Still in input units
      print('PLOT OK ', args.i[n], ind, maxv)
      stats.write('%-20s %12s %12s %12.3e %12.3e %12.3e\n' % (flabel[n], 
                    v, unit, np.mean(vals), np.min(vals), np.max(vals) ) )

     # We decide units on 1st file and whether to make plot or not.

      if n== 0:
         txtUnit = unit
         factor = 1.0
         if unit == 'ppb' and gmax[v] < 0.001:
            factor = 2.5e10
            txtUnit = 'molec/cm3'
         elif unit == 'ppb' and gmax[v] < 0.1:
            factor = 1.0e3
            txtUnit = 'ppt'
         txtUnit2 = txtUnit.replace('/','-') # Needed for filenames
         print('NEW UNIT ', v, 'MAX = ', gmax[v], maxv, 
                 ' in: ', unit, ' out ',  txtUnit, maxv*factor)
         if args.dbg : print('VALS: ', factor*vals)

         if factor*gmax[v] < 1.0e-6: 
           print(n, v, 'MAX = ', maxv, unit, ' TOO LOW. BREAK ', args.i[n]  )
           break
         else:

           fig=plt.figure()  # Start new fig or each variable
           # prefix for output file names if several_files
           #dirname= os.path.normpath(os.path.dirname(fname))  # Avoids empty, gets .
           #prefix = dirname + '/Comp_%s_%s' % ( txtUnit2, v ) 
           #if args.unitless: 
           #  prefix = dirname + '/Comp_%s' % v 
           #print("PREFIX ", prefix )

      if args.txt : 
           pfile='Out_%s_%s_%s.txt' % ( flabel[n], txtUnit2, v )
           with open(pfile,'w') as f:
              for val in vals:
                vf=val * factor # needed to stop error in next?
                f.write('%12.3f\n' % vf )

      if several_files:
         linelab=flabel[n]
      else:
         linelab=v

      vals = vals * factor 
      print(n, v, 'MAX = ', maxv, unit, gmax[v], fname, olabel, plot_ok )
      if nls == len(linestyle): nls = 0
      if nco == len(cols)     : nco = 0
      if maxv>0.0: plt.plot(vals,label=linelab,ls=linestyle[nls],color=cols[nco])
      #   plt.ylim(bottom=0.0)

      # Needed to prevent plt using strange constant offsets (e.g. 1701 -> 1)
      # See http://stackoverflow.com/questions/9303728/matplotlib-yaxis-range-
      #  display-using-absolute-values-rather-than-offset-values

   if plot_ok: 
      ax = plt.gca()
      ax.ticklabel_format(useOffset=False)
      plt.ylim(bottom=0.0)
      plt.tight_layout()
      ftitle = fname.replace('.csv','')
      if several_files:
         ftitle = 'Comp'  # files names would be too long!
      plt.title(ftitle + ':' + v + ' (%s)' % txtUnit )
      plt.xlabel('time in h')
      plt.ylabel('concentration in ' + txtUnit)
      if several_files:
          plt.legend(loc='upper center')
      if args.dbg:print('OK HERE', plot_ok, olabel, args.png , args.vars[0])
  
  # If the user wants "ALL" or "DEF"  species or adds "p" we make pngs directly

      if args.png  or args.vars[0] == 'ALL' or args.vars[0] == 'DEF' : # png wanted
          if args.dbg:print('OK HERE TOO', plot_ok, olabel)

        # savefig didn't like . in MCM_v3.3, and a simple replace can also
        # change . as directory. So, long-winded, we use:
          #old=os.path.basename(fname) # e.g. EmChem09 or MCMv3.3
          #new=old.replace('.','pt')
          nlabel=olabel.replace('.','pt')
          fname= 'Plot_' + nlabel + '_' +  v + '.png'  # test!
          print('FNAME P ', olabel, nlabel, fname )
          plt.savefig('%s' % fname, bbox_inches='tight')
      else:
          plt.show()
      plt.close(fig)  
