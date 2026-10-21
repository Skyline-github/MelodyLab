@echo off
setlocal
title MelodyLab Server
cd /d "%~dp0.."

set PYEXE=python
where python >nul 2>nul
if errorlevel 1 (
    where py >nul 2>nul
    if errorlevel 1 (
        echo [ERROR] Python not found. Run install_windows.bat first,
        echo         or install Python 3.10+ from python.org.
        pause
        exit /b 1
    )
    set PYEXE=py
)

echo ============================================
echo  MelodyLab starting on http://127.0.0.1:8300
echo  Browser will open in a few seconds...
echo  Keep this window OPEN. Press Ctrl+C to stop.
echo ============================================
echo.

start "" cmd /c "timeout /t 8 /nobreak >nul & start http://127.0.0.1:8300"

cd backend
%PYEXE% run.py --port 8300

echo.
echo [INFO] Server exited. If this was unexpected, read the error above.
pause
