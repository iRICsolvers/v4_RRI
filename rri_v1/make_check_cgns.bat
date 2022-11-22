call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2019

del *.obj *.mod

ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ../1.4.2.3/iric.f90
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ../check_cgns/check_cgns.f90
ifort /MD /Qopenmp /O3 *.obj ../lib/iriclib.lib /o  ./check_cgns.exe


del *.obj 
del *.mod 



