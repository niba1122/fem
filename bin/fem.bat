@echo off

setlocal

set current_dir=%~dp0
set command=%1
set solver=%2

if "%command%" == "compile" (
  cd %current_dir%
  cd ..\
  call ".\bin\compile.exe" %solver%
) else if "%command%" == "exec" (
  cd %current_dir%
  cd ..\solvers\%solver%\
pwd
  call ".\%solver%.exe"
)
