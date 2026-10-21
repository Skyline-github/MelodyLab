@echo off
setlocal
title MelodyLab - Install Dependencies
cd /d "%~dp0.."

echo ============================================
echo  MelodyLab dependency installer
echo ============================================

set PYEXE=python
where python >nul 2>nul
if errorlevel 1 (
    where py >nul 2>nul
    if errorlevel 1 (
        echo [ERROR] Python not found. Install Python 3.10+ from
        echo         https://www.python.org/downloads/
        echo         IMPORTANT: check "Add python.exe to PATH" in installer.
        pause
        exit /b 1
    )
    set PYEXE=py
)

echo Using interpreter:
%PYEXE% --version
if errorlevel 1 (
    echo [ERROR] Failed to run Python.
    pause
    exit /b 1
)

echo.
echo Installing backend dependencies (requires internet)...
cd backend
%PYEXE% -m pip install -r requirements.txt
if errorlevel 1 (
    echo.
    echo [ERROR] pip install failed. Check your network and retry,
    echo         or run manually:  python -m pip install -r requirements.txt
    cd /d "%~dp0.."
    pause
    exit /b 1
)

echo.
echo ============================================
echo  DONE. Now double-click start_windows.bat
echo ============================================
cd /d "%~dp0.."
pause
