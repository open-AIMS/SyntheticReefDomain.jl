@echo off
rem Starts the reef builder + C~scape (R) export notebook in your web browser.
rem Needs Julia 1.11 or newer: https://julialang.org/install
rem To build the R files directly, also install R and the cscape package
rem (https://github.com/open-AIMS/c-scape); without them the notebook writes the inputs
rem and a script to run later in R.
cd /d "%~dp0"
where julia >nul 2>nul || (echo Julia was not found. Install it from https://julialang.org/install and try again. & pause & exit /b 1)
julia --project=notebooks notebooks\run_cscape.jl
pause
