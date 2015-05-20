@echo off

setlocal

set command=%1

if "%command%" == "all" (
  del .\bin\*.exe
  gfortran -o .\bin\compile.exe .\src\compile_win.f90
) else if "%command%" == "clean" (
  del .\bin\*.exe
) else (
  del .\bin\*.exe
  gfortran -o .\bin\compile.exe .\src\compile_win.f90
)