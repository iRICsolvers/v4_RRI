@echo off
setlocal ENABLEDELAYEDEXPANSION

REM 1) Load Intel oneAPI env
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2022

REM 2) Make sure the correct redist path is at the FRONT of PATH
if defined INTEL_DEV_REDIST (
  set "OMP_REDIST=%INTEL_DEV_REDIST%\redist\intel64\compiler"
  if exist "%OMP_REDIST%\libiomp5md.dll" (
    set "PATH=%OMP_REDIST%;%PATH%"
  )
)

REM 3) (Optional) warn if a local libiomp5md.dll exists in app dir
pushd %~dp0
if exist libiomp5md.dll (
  echo [WARNING] A local libiomp5md.dll exists in the app folder and may cause conflicts.
  echo          Consider renaming it temporarily:  ren libiomp5md.dll libiomp5md.dll.bak
)
echo Using libiomp5md.dll from:
where libiomp5md.dll

REM 4) Run the app
.\rri.exe %*
set rc=%ERRORLEVEL%

popd
endlocal & exit /b %rc%
