@echo off
setlocal

call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
if errorlevel 1 exit /b 1

call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2022
if errorlevel 1 exit /b 1

del /q *.obj *.mod 2>nul

ifx /MT /c /O3 /nostandard-realloc-lhs ./src/RRI/iric.f90
if errorlevel 1 goto :error
ifx /MT /c /O3 /nostandard-realloc-lhs ./src/RRI_GridGenerator/RRI_GridGenerator_IRIC.f90
if errorlevel 1 goto :error
ifx /MT /c /O3 /nostandard-realloc-lhs ./src/RRI_GridGenerator/RRI_GridGenerator.f90
if errorlevel 1 goto :error
ifx /MT /O3 *.obj ./lib/iriclib.lib /o ./gridcreator/RRI_demAdjust2/RRI_demAdjust2.exe
if errorlevel 1 goto :error

editbin /stack:40000000 ./gridcreator/RRI_demAdjust2/RRI_demAdjust2.exe
if errorlevel 1 goto :error

del /q *.obj *.mod 2>nul
exit /b 0

:error
del /q *.obj *.mod 2>nul
del /q ./gridcreator/RRI_demAdjust2/RRI_demAdjust2.exe 2>nul
exit /b 1
