To compile and run a cross-section fortran file using gfortran

#############################################################################

gfortran will then create an executable from your code. This executable will be called a.out on unix systems and a.exe on Windows. Of course, you might consider giving it a more appropriate name: you can specify the name you want with the -o option: 

>> gfortran myfile.f -o program.exe

then run by

>> ./program.exe



