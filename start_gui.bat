@echo off
rem Starts the TuxNet-R GUI in WSL and opens it in your browser.
set "HERE=%~dp0"
set "HERE=%HERE:\=/%"
set "HERE=%HERE:~0,-1%"
for /f "tokens=1,* delims=:" %%a in ("%HERE%") do set "DRIVE=%%a" & set "REST=%%b"
for %%L in (a b c d e f g h i j k l m n o p q r s t u v w x y z) do if /i "%DRIVE%"=="%%L" set "DRIVE=%%L"
start "StRym server (close this window to stop)" wsl -d Ubuntu -u root -e bash "/mnt/%DRIVE%%REST%/start_gui.sh"
timeout /t 8 >nul
start http://localhost:3838
