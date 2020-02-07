#!/usr/bin/env python3
# -*- coding: utf-8 -*-
Usage="""\
Usage:

  from e.g. box/tmp_work directory,

  ../scripts/emep_setup.py XXXX  [-d]

  - where XXXX is e.g. EmChem19a or CRI-R5-emep or ...

  -d is for debug mode. Waits for user input.

Available XXX are: """  # Will add choices later

# Sets up EMEP chemical mechanisms, and does quick test with boxChem
import argparse
import os
import sys
import subprocess
import time

# start script from  e.g. .../esx/src or  .../box/src
cdir=os.environ['PWD']      # e.g. .../esx/src or  .../box/src
chemscripts='../../chem/scripts'  # ecosx/chem/scripts directory (do.GenChem)

# some useful collections
# Comments: Aqueous_EmChem16x is still standard
#  omments: Aqueous_EmChem16z uses Fgas(SO2) and still need checking
common= ' Aqueous_EmChem16x Aero2017nx ShipNOx FFireInert SeaSalt DustExtended Ash PM_WoodFFuelInert BVOC_SQT_NV' # Most compact but typical EMEP
dcommon=' Aqueous_EmChem16x Aero2017nx ShipNOx FFireInert SeaSalt Dust BVOC_SQT_NV' # More compact dust
ocommon=' ShipNOx FFireInert SeaSalt DustExtended Pollen' # Added pollen and extended dust
fcommon=' Aqueous_EmChem16x Aero2017nx Ash ShipNOx FFireInert SeaSalt DustExtended Pollen BVOC_SQT_NV'
ncommon=' Aqueous_EmChem16x Aero2017nx Ash ShipNOx FFireInert16z SeaSalt DustExtended16z Pollen'

# Common emep schemes:
# (feel free to define your own combinations, but call them something else!)
cmdx=dict()
#cmdx['EmChem19']  ='-b EmChem19  -e APINENE_BVOC_emis BVOC_EmChem19 VBS_EmChem16z'+common
# Nov 2019, use simpler Dust rather than DustExtended
#cmdx['EmChem19a'] ='-b EmChem19a -e APINENE_BVOC_emis BVOC_EmChem19 VBS_EmChem19'+common
# order from rv4.33



#NOT WORKING or SUPPORTED: cmdx['EmChem09']  ='-b EmChem09  -e PM_VBS_acp2012 APINENE_BVOC_emis'+common
#OLD cmdx['EmChem19']  ='-b EmChem19  -e PM_VBS_EmChem19 BVOC_IsoMT1_emis'+common
cmdx['EmChem19a'] ='-b EmChem19a -e PM_VBS_EmChem19 BVOC_IsoMT1_emis'+common
cmdx['EmChem19p'] = cmdx['EmChem19a'] + ' Pollen'
cmdx['CB6r2']     ='-b CB6r2emep -e PM_Hodzic_CB6'+common # BoxAero # CB6 has rcbio in base
#cmdx['TestEm'] ='-b EmChem19a -e FFireInert SeaSalt' # TEST
cmdx['CRI-R5-emep'] ='-b CRI-R5-emep -e '+common
available_emeps= list(cmdx.keys())

Usage += '  '.join(available_emeps)

parser = argparse.ArgumentParser(usage=Usage)

if len(sys.argv) < 2:
    parser.print_usage()
    sys.exit(1)

parser.add_argument('chem_mech', help='Base chemical mechanisms (one  needed), e.g. EmChem19',choices=available_emeps)
parser.add_argument('-d','--dbg',
    help='debug mode. Waits for user input at key steps', action='store_true')
args = parser.parse_args()
print(args)

chem=args.chem_mech
dbg=args.dbg

# Prepare
cmds=['pwd', 'make clean', 'rm boxChem', 'rm CM_*f90 CM_*inc CMX_*inc *.o *.mod', 'rm build/*' ]
for cmd in cmds:
  print('Runs:', cmd)
  subprocess.call(cmd,shell=True)  # needs string, not args list, wit shell=True


gchem='../../chem/scripts/do.GenChem'

txt= gchem + ' ' + cmdx[chem] 
if dbg: txt += ' -d'
args=txt.split()
for a in args: print('Arg ', a)
if dbg:  input('Press enter key to continue...')

########### RUNS GENCHEM #########################
subprocess.call(args)
########### END GENCHEM #########################

assert os.path.isfile('CM_ChemDims_mod.f90'),'Failed GENCHEM!'

print('FINISHED GENCHEM', chem)
print(' Making....', chem)
time.sleep(3)
subprocess.call('make')
print(' MADE...', chem)

## Create directory to store CM and emissplit files: use these for EMEP model
# All files and  emissplit_run directory created by do.GenChem 

odir='ZCM_' + chem
splitdir=odir+'/emissplit_run'
os.makedirs(odir,exist_ok=True)
os.makedirs(splitdir,exist_ok=True)
print('COPY CM FILES ', odir)
subprocess.call('wc CM*inc',shell=True)
subprocess.call('cp CM_* %s' % odir,shell=True)
subprocess.call('cp CMX_* %s' % odir,shell=True)
subprocess.call('mv emissplit_run/emissplit*  %s' % splitdir,shell=True )

## Run example
os.makedirs('Output',exist_ok=True) # just in case ;-)

