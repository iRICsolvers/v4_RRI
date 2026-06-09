call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2022

del *.obj *.mod
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/iric.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_Mod.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_Mod2.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RSR_Mod.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_Mod_Dam.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_Mod_Tecout.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_IRIC_Mod.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_Bound.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_Dam.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_Div.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_DT_Check.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_Evp.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_GW.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_Infilt.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_Read.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_Riv.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_RivSlo.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_section.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_Slope.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_Sub.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_Sediment.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_Sed2.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_Tecout.f90
ifx /MT /c /heap-arrays /Qopenmp /O3 /nostandard-realloc-lhs /traceback ./src/RRI/RRI_TSAS.f90

ifx /Qopenmp-link:static /libs:static /MT /Qopenmp /O3 *.obj ./lib/iriclib.lib /o  ./rri.exe
editbin /stack:40000000  ./rri.exe

del *.obj 
del *.mod 



REM Copy iriclib.dll next to the EXE if we have it locally
if exist ".\lib\iriclib.dll" copy /Y ".\lib\iriclib.dll" "."
