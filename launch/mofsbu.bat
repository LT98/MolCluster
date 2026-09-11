@echo off
rem Start the mofsbu viewer on Windows.  Double-click this file.
rem
rem   mofsbu.bat                  CPU, finds the registry on its own
rem   mofsbu.bat --gpu            use CUDA for ml_go runs
rem   mofsbu.bat --gpu --workers 8
rem   mofsbu.bat --db scratch.db --port 8001
rem
rem Same reasoning as mofsbu.sh: the interpreter is used by FULL PATH rather than
rem through `conda activate`, which on Windows needs either an Anaconda Prompt or a
rem `conda init`-ed shell, and a double-click gives you neither.
setlocal EnableDelayedExpansion

set "HERE=%~dp0"
for %%I in ("%HERE%..") do set "REPO=%%~fI"

rem ── options ────────────────────────────────────────────────────────────────
set "DEVICE="
set "WORKERS="
set "PORT=8000"
set "OPENFLAG=--open"
set "PASS="

:parse
if "%~1"=="" goto parsed
if /I "%~1"=="--gpu"        ( set "DEVICE=cuda" & shift & goto parse )
if /I "%~1"=="--cpu"        ( set "DEVICE=cpu"  & shift & goto parse )
if /I "%~1"=="--device"     ( set "DEVICE=%~2"  & shift & shift & goto parse )
if /I "%~1"=="--workers"    ( set "WORKERS=%~2" & shift & shift & goto parse )
if /I "%~1"=="--port"       ( set "PORT=%~2"    & shift & shift & goto parse )
if /I "%~1"=="--no-browser" ( set "OPENFLAG="   & shift & goto parse )
set "PASS=!PASS! %~1"
shift
goto parse
:parsed

rem ── which python? ──────────────────────────────────────────────────────────
rem Probed, not assumed: environment.yml names the env `mofsbu`, but an env is whatever
rem it was actually created as.  The imports are the test.
set "PY="
if defined MOFSBU_PYTHON (
  call :usable "%MOFSBU_PYTHON%" && set "PY=%MOFSBU_PYTHON%"
)

if not defined PY (
  for %%R in ("%USERPROFILE%\miniforge3" "%USERPROFILE%\mambaforge" "%USERPROFILE%\miniconda3" "%USERPROFILE%\anaconda3" "%LOCALAPPDATA%\miniforge3" "%LOCALAPPDATA%\Continuum\anaconda3" "C:\ProgramData\miniforge3" "C:\ProgramData\Anaconda3") do (
    if not defined PY (
      for %%N in ("%MOFSBU_ENV%" "mofsbu") do (
        if not defined PY if not "%%~N"=="" (
          call :usable "%%~R\envs\%%~N\python.exe" && set "PY=%%~R\envs\%%~N\python.exe"
        )
      )
    )
  )
)

rem Any environment that can actually do the job — the case that matters when the env
rem was created under a different name.
if not defined PY (
  for %%R in ("%USERPROFILE%\miniforge3" "%USERPROFILE%\mambaforge" "%USERPROFILE%\miniconda3" "%USERPROFILE%\anaconda3" "%LOCALAPPDATA%\miniforge3" "C:\ProgramData\miniforge3" "C:\ProgramData\Anaconda3") do (
    if exist "%%~R\envs" (
      for /D %%E in ("%%~R\envs\*") do (
        if not defined PY call :usable "%%~E\python.exe" && set "PY=%%~E\python.exe"
      )
    )
  )
)

if not defined PY (
  where python >nul 2>&1 && (
    for /F "delims=" %%P in ('where python') do (
      if not defined PY call :usable "%%P" && set "PY=%%P"
    )
  )
)

if not defined PY (
  echo.
  echo mofsbu: could not find a Python with the dependencies installed.
  echo.
  echo Looked for a conda environment containing fastapi, uvicorn and rdkit.
  echo Create one once, from an Anaconda Prompt in the repo folder:
  echo.
  echo     conda env create -f "%REPO%\environment.yml"
  echo.
  echo If the environment exists under a different name, say which:
  echo.
  echo     set MOFSBU_ENV=the-name
  echo.
  echo.
  pause
  exit /b 1
)

rem ── go ─────────────────────────────────────────────────────────────────────
set "ARGS=--port %PORT%"
if defined OPENFLAG set "ARGS=%ARGS% %OPENFLAG%"
if defined DEVICE   set "ARGS=%ARGS% --device %DEVICE%"
if defined WORKERS  set "ARGS=%ARGS% --workers %WORKERS%"

echo mofsbu: %PY%
cd /d "%REPO%"
set PYTHONUNBUFFERED=1
"%PY%" scripts\viewer.py %ARGS% %PASS%
rem The window must not vanish on an error, taking the message with it.
if errorlevel 1 pause
exit /b %errorlevel%

:usable
if not exist "%~1" exit /b 1
"%~1" -c "import fastapi, uvicorn, rdkit" >nul 2>&1
exit /b %errorlevel%
