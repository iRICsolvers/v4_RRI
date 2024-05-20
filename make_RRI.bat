call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2019

del *.obj *.mod
ifort /MD /c /Qopenmp /O3 /nostandard-realloc-lhs ./src/RRI/iric.f90
ifort /MD /c /Qopenmp /O3 /nostandard-realloc-lhs ./src/RRI/RRI_Mod.f90
ifort /MD /c /Qopenmp /O3  /nostandard-realloc-lhs ./src/RRI/RRI_Mod2.f90
ifort /MD /c /Qopenmp /O3  /nostandard-realloc-lhs ./src/RRI/RRI_Mod_Dam.f90
ifort /MD /c /Qopenmp /O3  /nostandard-realloc-lhs ./src/RRI/RRI_Mod_Tecout.f90
ifort /MD /c /Qopenmp /O3  /nostandard-realloc-lhs ./src/RRI/RRI_IRIC_Mod.f90
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ./src/RRI/RRI.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ./src/RRI/RRI_Bound.f90
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ./src/RRI/RRI_Dam.f90
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ./src/RRI/RRI_Div.f90
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ./src/RRI/RRI_DT_Check.f90
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ./src/RRI/RRI_Evp.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ./src/RRI/RRI_GW.f90
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ./src/RRI/RRI_Infilt.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ./src/RRI/RRI_Read.f90
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ./src/RRI/RRI_Riv.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ./src/RRI/RRI_RivSlo.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ./src/RRI/RRI_section.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ./src/RRI/RRI_Slope.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ./src/RRI/RRI_Sub.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ./src/RRI/RRI_Tecout.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ./src/RRI/RRI_TSAS.f90

ifort /MD /Qopenmp /O3 *.obj ./lib/iriclib.lib /o  ./install/rri.exe
editbin /stack:40000000  ./install/rri.exe

del *.obj 
del *.mod 



