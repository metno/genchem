# Script to build a makefile for a desired driver (e.g. esx_tester.f90)

#     -p ".;./DO3SE/src;./MAFOR/src;./ESX_MEGAN/src" \

python "$(dirname $0)/gen-makefile.py" \
     -i iso_fortran_env \
     -f90 "gfortran -pedantic -Wall -fbounds-check -fdefault-real-8 -finit-real=nan -ffree-line-length-none" \
     -b build \
     "$@"
