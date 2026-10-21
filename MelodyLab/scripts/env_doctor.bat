@echo off
setlocal
title MelodyLab Environment Doctor
cd /d "%~dp0.."

set PYEXE=python
where python >nul 2>nul
if errorlevel 1 (
    where py >nul 2>nul
    if errorlevel 1 (
        echo [ERROR] Python not found. Install Python 3.10+ first.
        pause
        exit /b 1
    )
    set PYEXE=py
)

echo Scanning local conda/python environments ...
echo.
%PYEXE% scripts\doctor.py %*

echo.
pause
