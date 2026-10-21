@echo off
title MelodyLab Server (port 8300)
cd /d D:\Files\TEST\MelodyLab_win64\MelodyLab\backend

echo ==============================================
echo  MelodyLab starting on http://127.0.0.1:8300
echo  Keep this window OPEN while using the app.
echo  Close this window to stop the server.
echo ==============================================

rem open browser after 3 s (server boots faster than that)
start "" cmd /c "timeout /t 3 >nul & start http://127.0.0.1:8300"

D:\Anaconda\envs\tav\python.exe run.py --port 8300

echo.
echo Server stopped.
pause
