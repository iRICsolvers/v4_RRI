call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2019

del *.obj *.mod
ifort /MD /c /Qopenmp /O3 /nostandard-realloc-lhs ./src/RRI/iric.f90
ifort /MD /c /Qopenmp /O3 /nostandard-realloc-lhs ./src/demAdjust2/RRI_IRIC_Mod2.f90
ifort /MD /c /Qopenmp /O3 /nostandard-realloc-lhs ./src/demAdjust2/demAdjust2.f90
ifort /MD /Qopenmp /O3 *.obj ./lib/iriclib.lib /o  ./install/demAdjust2.exe
editbin /stack:40000000 ./install/demAdjust2.exe

del *.obj 
del *.mod 



