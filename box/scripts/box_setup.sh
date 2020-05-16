#/usr/bin/env bash
Usage="
Usage:

   box_setup.sh -h 

    or

   box_setup.sh  workdir_name (e.g. tmp_work)

(This script just saves time for the first setup. For more info you can
  browse the GenChem documentation
  with e.g. firefox doc/UsersGuide/html/index.html)"

#echo "NARGS $#"
if [[ $# ==  0  ]]; then  # we have arguments ($# is number of args)
  echo "No arguments $Usage"
  exit 0
elif [[ $1  ==  "-h" ]]; then  # we have arguments ($# is number of args)
  printf "%s\n" "help: $Usage"
  exit 0
else
   #echo ARG NOW $1
   workdir=$1   # e.g. tmp_boxwork
fi
echo Creates workdir: $workdir

# Sets up fortan code and scripts in temporary work-dir
# Run from this location just once to setup initial code and links.
# Thereafter see XXXX

mkdir -p $workdir
cd $workdir
cp    ../src/*f90             . # Copy src f90 files from box/src to here
cp -i ../src/config_box.nml   . # ask before potential overwrite
cp -i ../src/Makefile         . # ask before potential overwrite
cp  ../scripts/do.testChems .
cp  ../scripts/emep_setup.py .
# do.GenChem  is not usually needed as it is better to use do.testChems, but
# sometimes it helps as a first test
cp ../../chem/scripts/do.GenChem . 
