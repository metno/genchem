# Scripts

## Makefile creation scripts

  * do.Makefile - a helper script to run gen-makefile.py. Here one can set e.g.
    the gfortran flags, and if needed sub-directories to search

  * do.testChems - harminoses usage of BoxAero with base-mechanisms, to make
    testing more comparable



## Plotting scripts

### boxplots.py

  NEW - merged elements of original compare_2schemes and  extract_poll. 

  Reads one or more .csv files, plotting individually or as comparison plots

  From box/src, type e.g.

    ../../scripts/boxplots.py -v O3 -i mcm.csv   [--png]

    ../../scripts/boxplots.py -v O3 NO2 OD -i file1.csv file2.csv file3.csv  [--png]

  One can plot all of a defined list, or just all species:

    ../../scripts/boxplots.py -v DEF -i file1.csv file2.csv file3.csv

    ../../scripts/boxplots.py -v ALL -i file1.csv file2.csv file3.csv



## Misc scripts

### getboxconcs.py

  Reads data for a specific pollutant, and stores in ResConcs file::

     ../../scripts/getboxconcs.py O3 file.csv

  would produce ResConcs\_file_O3_ppb.txt

