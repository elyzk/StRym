@echo off
setlocal EnableExtensions
rem Starts the StRym GUI in WSL and opens it in the browser. STRYM_PORT defaults to 3838.
if not defined STRYM_PORT set "STRYM_PORT=3838"
powershell.exe -NoProfile -Command "$p = 0; if (-not [int]::TryParse($env:STRYM_PORT, [ref]$p) -or $p -lt 1 -or $p -gt 65535) { exit 2 }"
if errorlevel 1 (
  echo Invalid STRYM_PORT: %STRYM_PORT%. Use a number from 1 to 65535.
  exit /b 2
)
set "HERE=%~dp0"
set "HERE=%HERE:\=/%"
set "HERE=%HERE:~0,-1%"
for /f "tokens=1,* delims=:" %%a in ("%HERE%") do set "DRIVE=%%a" & set "REST=%%b"
for %%L in (a b c d e f g h i j k l m n o p q r s t u v w x y z) do if /i "%DRIVE%"=="%%L" set "DRIVE=%%L"
start "StRym server (close this window to stop)" wsl -d Ubuntu -u root -e env STRYM_PORT=%STRYM_PORT% bash "/mnt/%DRIVE%%REST%/start_gui.sh"
timeout /t 8 >nul
start "" "http://localhost:%STRYM_PORT%"
