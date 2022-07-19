call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2019

del *.obj *.mod
ifort /c /O3 ../demAdjust2/RRI_IRIC_Mod2.f90  -I ../include
ifort /c /O3 ../demAdjust2/demAdjust2.f90  -I ../include
ifort /O3 *.obj ../lib/cgnsdll_x64_ifort.lib ../lib/iriclib_x64_ifort.lib /o  ./demAdjust2.exe
editbin /stack:40000000 ./demAdjust2.exe

del *.obj 
del *.mod 



