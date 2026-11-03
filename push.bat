@echo off
echo ========================================================
echo Pushing Nigeria Malaria PDE Project to GitHub: IsaacMoon1
echo ========================================================
echo.

git branch -M main
git remote remove origin 2>nul
git remote add origin https://github.com/IsaacMoon1/malaria-pde-nigeria.git

echo Remote origin configured:
git remote -v
echo.
echo Pushing to GitHub (you may be prompted to log in)...
git push -u origin main

if %ERRORLEVEL% EQU 0 (
    echo.
    echo ========================================================
    echo SUCCESS! Repository pushed to:
    echo https://github.com/IsaacMoon1/malaria-pde-nigeria
    echo ========================================================
) else (
    echo.
    echo [!] If the push failed because the repository does not exist yet:
    echo 1. Go to https://github.com/new
    echo 2. Create a new repository named: malaria-pde-nigeria
    echo 3. Run this script again: push.bat
)
pause
