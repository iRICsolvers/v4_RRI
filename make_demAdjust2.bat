@echo off
setlocal

call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
if errorlevel 1 exit /b 1

call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64
if errorlevel 1 exit /b 1

del /q *.obj *.mod 2>nul

ifx /MD /c /Qopenmp /O3 /nostandard-realloc-lhs ./src/RRI/iric.f90
if errorlevel 1 goto :error
ifx /MD /c /Qopenmp /O3 /nostandard-realloc-lhs ./src/demAdjust2/RRI_IRIC_Mod2.f90
if errorlevel 1 goto :error
ifx /MD /c /Qopenmp /O3 /nostandard-realloc-lhs ./src/demAdjust2/demAdjust2.f90
if errorlevel 1 goto :error
ifx /MD /Qopenmp /O3 *.obj ./lib/iriclib.lib /o ./install/demAdjust2_new.exe
if errorlevel 1 goto :error

editbin /stack:40000000 ./install/demAdjust2_new.exe
if errorlevel 1 goto :error

move /y .\install\demAdjust2_new.exe .\install\demAdjust2.exe >nul
del /q *.obj *.mod 2>nul
exit /b 0

:error
del /q *.obj *.mod 2>nul
del /q .\install\demAdjust2_new.exe 2>nul
exit /b 1
