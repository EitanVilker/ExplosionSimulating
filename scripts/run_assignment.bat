@echo off

set assignment=%1

if "%~1"=="" goto :error

set test=%2

if "%test%"=="" set test=1

cmake --build ./build --config Release --target %assignment%
if not %ERRORLEVEL% == 0 goto :error

build\proj\%assignment%\Release\%assignment%.exe  -test %test%
if not %ERRORLEVEL% == 0 goto :error

goto :endofscript

:error
echo =========================
echo There was an error!
echo =========================

:endofscript
