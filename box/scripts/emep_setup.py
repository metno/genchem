#!/usr/bin/env python3
# -*- coding: utf-8 -*-
Usage="""\
Usage:

  from e.g. box/tmp_work directory,

  ../scripts/emep_setup.py XXXX  [-d]

  - where XXXX is e.g. EmChem19a-vbs or CRIv2R5Em-H3 or (use emep_setup.py -h to see list)

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
#           Aqueous_EmChem16z uses Fgas(SO2) and still need checking
common= ' Aqueous_EmChem16x Aero2017nx ShipNOx PM_FFireInert SeaSalt DustExtended Ash PM_WoodFFuelInert EC_ageing' # Most typical EMEP
# with ResNonRes PM, EC aging moved in PM_ResNonResInert
commonRNR= ' Aqueous_EmChem16x Aero2017nx ShipNOx PM_FFireInert SeaSalt DustExtended Ash PM_ResNonResInert' 
dcommon=' Aqueous_EmChem16x Aero2017nx ShipNOx PM_FFireInert SeaSalt Dust BVOC_SQT_NV' # More compact dust
ocommon=' ShipNOx PM_FFireInert SeaSalt DustExtended Pollen' # Added pollen and extended dust
fcommon=' Aqueous_EmChem16x Aero2017nx Ash ShipNOx PM_FFireInert SeaSalt DustExtended Pollen BVOC_SQT_NV'
ncommon=' Aqueous_EmChem16x Aero2017nx Ash ShipNOx FFireInert16z SeaSalt DustExtended16z Pollen'

# Many schemes use common + same BVOC options. We add
# IsoMT1 adds code for BVOC_SQT_NV, isoprene and apinene as surrogate for all monoterpenes
# Mainly emissions in these common packages, since chem reactions and species
# differ betweem EmChems, CRIs and CB6 schemes
common_IsoMT1  = common + ' BVOC_SQT_NV BVOC_IsoMT1_emis'
common_RNR     = commonRNR + ' BVOC_SQT_NV BVOC_IsoMT1_emis'
common_IsoMT2  = common + ' BVOC_SQT_NV BVOC_IsoMT2_emis'
common_IsoMT3  = common + ' BVOC_SQT_NV BVOC_IsoMT3_emis'

# Mix n' match system - can combine various gas and particle packages
# Common emep combinations given here.
# (feel free to define your own combinations, but call them something else!)
# postfixes -vbs.. indicates VBS SOA scheme similar to the "standard" EMEP scheme
#           -H..   indicates SOA scheme loosely based on Hodzic et al., 2016
#                   doi:10.5194/acp-16-7917-2016
#           -M19'  indicates SOA scheme used in McFiggans et al., 2019

cmdx=dict()

cmdx['EmChem19a-vbs']  ='-b EmChem19a -e PM_VBS_EmChem19 '+ common_IsoMT1
cmdx['EmChem19c-vbs']  ='-b EmChem19c -e PM_VBS_EmChem19 '+ common_IsoMT1 # mimics EmChem19a but with EmChem19c base mechanism
cmdx['EmChem19p-vbs']  = cmdx['EmChem19a-vbs'] + ' Pollen'
# To preserve older notation, we allow two simple aliases:
cmdx['EmChem19a']  = cmdx['EmChem19a-vbs']
cmdx['EmChem19c']  = cmdx['EmChem19c-vbs']
cmdx['EmChem19p']  = cmdx['EmChem19p-vbs']
cmdx['EmChem19cAsh7'] = '-b EmChem19c -e EmAsh ' #+ ' ShipNOx PM_FFireInert SeaSalt  DustExtended PM_WoodFFuelInert EC_ageing'
cmdx['EmergencyAsh'] = '-b Emergency -e EmAsh '
cmdx['EmergencyRadiation'] = '-b Emergency -e Radiation '
cmdx['Emergency'] = '-b Emergency '
# schemes which uses Res, nonRes split instead of wood/ffuel
cmdx['EmChem19r']  ='-b EmChem19a -e PM_VBS_EmChem19 '+ common_RNR
cmdx['EmChem19rp']  ='-b EmChem19a -e PM_VBS_EmChem19 '+ common_RNR + ' Pollen'

cmdx['EmChem19rc']  ='-b EmChem19c -e PM_VBS_EmChem19 ' + common_RNR
cmdx['EmChem19rcp'] ='-b EmChem19c -e PM_VBS_EmChem19 ' + common_RNR + ' Pollen'

cmdx['EmChem19c-vbs3'] ='-b EmChem19c -e BVOC_ExtraMTs PM_VBS_EmChem19 PM_VBS_ExtraMTs'+common_IsoMT3 
cmdx['EmChem19c-H']    ='-b EmChem19c -e PM_Hodzic_EmChem19'+common_IsoMT1
#
#cmdx['EmChem19X-vbs']  ='-b EmChem19X -e PM_VBS_EmChem19 '+ common_IsoMT1
#cmdx['TestEm'] ='-b EmChem19a -e FFireInert SeaSalt' # minimal gas/particle
# CB6s:
cmdx['CB6r2Em-vbs']   ='-b CB6r2Em -e PM_VBS_CB6r2Em'+common_IsoMT1
cmdx['CB6r2Em-H']     ='-b CB6r2Em -e PM_Hodzic_CB6' +common_IsoMT1
# CRIs:
cmdx['CRIv2R5Em-vbs'] ='-b CRIv2R5Em -e PM_VBS_EmChem19'+common_IsoMT1
cmdx['CRIv2R5Em-M19'] ='-b CRIv2R5Em -e PM_JPAC_MT3 PM_Hodzic_Aromatics BVOC_XTERP_CRI'+common_IsoMT3
cmdx['CRIv2R5Em-H1'] ='-b CRIv2R5Em -e PM_Hodzic_EmChem19'+common_IsoMT1
cmdx['CRIv2R5Em-H2'] ='-b CRIv2R5Em -e PM_Hodzic_EmChem19 PM_Hodzic_BPINENE'+common_IsoMT2
cmdx['CRIv2R5Em-H3'] ='-b CRIv2R5Em -e PM_Hodzic_EmChem19 PM_Hodzic_BPINENE PM_Hodzic_XTERP BVOC_XTERP_CRI'+common_IsoMT3
available_emeps= list(cmdx.keys())

Usage += '  '.join(available_emeps)

parser = argparse.ArgumentParser(usage=Usage)

if len(sys.argv) < 2:
    parser.print_usage()
    sys.exit(1)

parser.add_argument('chem_mech', help='Base chemical mechanisms (one  needed), e.g. EmChem19a-vbs',choices=available_emeps)
parser.add_argument('-g','--gnfr',
    help='use GNFR sectors. Only EmChem19a type so far', action='store_true')
parser.add_argument('-d','--dbg',
    help='debug mode. Waits for user input at key steps', action='store_true')
args = parser.parse_args()
print('sys.args:', sys.argv )
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
if dbg:      txt += ' -d'
if args.gnfr:
    txt += ' -g gnfr' # gnfr 
else:
    txt += ' -g snap' # SNAP was originally used in boxChem 
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
subprocess.call('cp emissplit_run/emissplit*  %s' % splitdir,shell=True )

## Run example
os.makedirs('Output',exist_ok=True) # just in case ;-)

