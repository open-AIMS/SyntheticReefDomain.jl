@echo off
rem Starts the synthetic reef builder (a Pluto notebook) in your web browser.
rem Needs Julia 1.11 or newer: https://julialang.org/install
rem The first start installs the packages, which takes 10-15 minutes.
cd /d "%~dp0"
where julia >nul 2>nul || (echo Julia was not found. Install it from https://julialang.org/install and try again. & pause & exit /b 1)
julia --project=notebooks notebooks\run.jl
pause
