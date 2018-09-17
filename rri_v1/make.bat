del *.obj *.mod
ifort /MD /c /Qopenmp /O3 /nostandard-realloc-lhs ../1.4.2.3/RRI_Mod.f90
ifort /MD /c /Qopenmp /O3  /nostandard-realloc-lhs ../1.4.2.3/RRI_Mod2.f90
ifort /MD /c /Qopenmp /O3  /nostandard-realloc-lhs ../1.4.2.3/RRI_Mod_Dam.f90
ifort /MD /c /Qopenmp /O3  /nostandard-realloc-lhs ../1.4.2.3/RRI_Mod_Tecout.f90
ifort /MD /c /Qopenmp /O3  /nostandard-realloc-lhs ../1.4.2.3/RRI_IRIC_Mod.f90 -I ../include
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ../1.4.2.3/RRI.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ../1.4.2.3/RRI_Bound.f90
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ../1.4.2.3/RRI_Dam.f90
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ../1.4.2.3/RRI_Div.f90
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ../1.4.2.3/RRI_DT_Check.f90
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ../1.4.2.3/RRI_Evp.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ../1.4.2.3/RRI_GW.f90
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ../1.4.2.3/RRI_Infilt.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ../1.4.2.3/RRI_Read.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ../1.4.2.3/RRI_Riv.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ../1.4.2.3/RRI_RivSlo.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ../1.4.2.3/RRI_section.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ../1.4.2.3/RRI_Slope.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ../1.4.2.3/RRI_Sub.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ../1.4.2.3/RRI_Tecout.f90 
ifort /MD /c  /Qopenmp /O3   /nostandard-realloc-lhs ../1.4.2.3/RRI_TSAS.f90

ifort /MD /Qopenmp /O3 *.obj ../lib/cgnsdll_x64_ifort.lib ../lib/iriclib_x64_ifort.lib /o  ./rri.exe
editbin /stack:40000000  ./rri.exe

del *.obj 
del *.mod 



