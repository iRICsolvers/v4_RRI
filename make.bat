@echo off

call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2019

rem ----------------------------------------------------------------------
rem ifort compile
rem ----------------------------------------------------------------------
ifort .\src\iric.f90  /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\global_variables.f90 /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\read_cgns.f90 /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\set_initial.f90 /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\calc_main.f90 /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\write_cgns.f90 /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\main.f90 /Qopenmp /nostandard-realloc-lhs /MD /c
ifort *.obj .\lib\iriclib.lib -o solver.exe
del *.obj
del *.mod

copy /y .\solver.exe .\INSTALL
