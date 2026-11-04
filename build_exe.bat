@echo off
setlocal
cd /d "%~dp0"

echo ============================================
echo Mica - ZO1 Junction GUI EXE Builder
echo ============================================
echo.

where py >nul 2>&1
if %errorlevel%==0 (
    set PY=py
) else (
    set PY=python
)

%PY% -m pip install --upgrade pip
if errorlevel 1 goto :error
%PY% -m pip install -r requirements.txt
if errorlevel 1 goto :error
%PY% -m PyInstaller --clean --noconfirm Mica_ZO1_Junction_GUI.spec
if errorlevel 1 goto :error

echo.
echo Build completed.
echo EXE: dist\Mica_ZO1_Junction_GUI.exe
echo.
pause
exit /b 0

:error
echo.
echo Build failed. Check the error messages above.
pause
exit /b 1
