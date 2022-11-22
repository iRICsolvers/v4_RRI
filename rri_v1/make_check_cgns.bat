call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2019

del *.obj *.mod
ifort /c /O3 ../check_cgns/check_cgns.f90
ifort /O3 *.obj ../lib/iriclib.lib /o  ./check_cgns.exe


del *.obj 
del *.mod 



