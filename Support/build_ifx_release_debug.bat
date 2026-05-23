@echo off
setlocal EnableExtensions EnableDelayedExpansion

rem Build both Release and Debug EMsoftOO trees with Intel ifx + NMake.
rem Usage:
rem   build_ifx_release_debug.bat [SDK_ROOT] [REPO_ROOT]
rem Example:
rem   build_ifx_release_debug.bat C:\Users\westraadt.1\Repos\EMsoftSuperbuild\EMsoftOO_SDK C:\Users\westraadt.1\Repos\EMsoftOO

set "SCRIPT_DIR=%~dp0"
set "DEFAULT_REPO_ROOT=%SCRIPT_DIR%.."
for %%I in ("%DEFAULT_REPO_ROOT%") do set "DEFAULT_REPO_ROOT=%%~fI"

set "SDK_ROOT=C:\Users\westraadt.1\Repos\EMsoftSuperbuild\EMsoftOO_SDK"
set "REPO_ROOT=%DEFAULT_REPO_ROOT%"

if not "%~1"=="" set "SDK_ROOT=%~1"
if not "%~2"=="" set "REPO_ROOT=%~2"

for %%I in ("%SDK_ROOT%") do set "SDK_ROOT=%%~fI"
for %%I in ("%REPO_ROOT%") do set "REPO_ROOT=%%~fI"

if not exist "%REPO_ROOT%\CMakeLists.txt" (
  echo [ERROR] Repo root is invalid: "%REPO_ROOT%"
  echo [ERROR] Could not find CMakeLists.txt
  exit /b 1
)

if not exist "%SDK_ROOT%\EMsoftOO_SDK.cmake" (
  echo [ERROR] SDK root is invalid: "%SDK_ROOT%"
  echo [ERROR] Could not find EMsoftOO_SDK.cmake
  exit /b 1
)

set "ONEAPI_SETVARS=C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
if not exist "%ONEAPI_SETVARS%" (
  echo [ERROR] Intel oneAPI setvars script not found:
  echo         "%ONEAPI_SETVARS%"
  exit /b 1
)

set "IFX_EXE=C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\ifx.exe"
if not exist "%IFX_EXE%" (
  echo [ERROR] ifx compiler not found:
  echo         "%IFX_EXE%"
  exit /b 1
)

set "SDK_ROOT_CMAKE=%SDK_ROOT:\=/%"
set "IFX_EXE_CMAKE=%IFX_EXE:\=/%"

echo [INFO] Repo root : "%REPO_ROOT%"
echo [INFO] SDK root  : "%SDK_ROOT%"
echo [INFO] Loading Intel+VS toolchain...
call "%ONEAPI_SETVARS%" intel64 vs2022 >nul
if errorlevel 1 (
  echo [ERROR] Failed to initialize Intel/VS toolchain.
  exit /b 1
)

call :configure_and_build Release build-ifx-release
if errorlevel 1 exit /b 1

call :configure_and_build Debug build-ifx-debug
if errorlevel 1 exit /b 1

echo [OK] Build complete.
echo [OK] Release binaries: "%REPO_ROOT%\build-ifx-release\Bin"
echo [OK] Debug binaries  : "%REPO_ROOT%\build-ifx-debug\Bin"
exit /b 0

:configure_and_build
set "CFG=%~1"
set "BUILD_DIR=%~2"
set "BUILD_PATH=%REPO_ROOT%\%BUILD_DIR%"

echo [INFO] ===============================================================
echo [INFO] Configuring %CFG%: "%BUILD_PATH%"

cmake -S "%REPO_ROOT%" -B "%BUILD_PATH%" -G "NMake Makefiles" ^
  -DBUILD_SHARED_LIBS=ON ^
  -DCMAKE_BUILD_TYPE=%CFG% ^
  -DCMAKE_Fortran_COMPILER="%IFX_EXE_CMAKE%" ^
  -DEMsoftOO_SDK="%SDK_ROOT_CMAKE%" ^
  -DEMsoftOO_ENABLE_TESTING=OFF ^
  -DJSONFORTRAN_INSTALL="%SDK_ROOT_CMAKE%/jsonfortran-4.2.1-%CFG%" ^
  -DJSONFORTRAN_DIR="%SDK_ROOT_CMAKE%/jsonfortran-4.2.1-%CFG%/lib/cmake/jsonfortran-intelllvm-4.2.1" ^
  -Djsonfortran-intelllvm_DIR="%SDK_ROOT_CMAKE%/jsonfortran-4.2.1-%CFG%/lib/cmake/jsonfortran-intelllvm-4.2.1" ^
  -DNLopt_DIR="%SDK_ROOT_CMAKE%/nlopt-2.7.0-%CFG%/lib/cmake/nlopt" ^
  -DBSPLINEFORTRAN_INSTALL="%SDK_ROOT_CMAKE%/BSPLINEFORTRAN-7.4.0-%CFG%" ^
  -DBSPLINEFORTRAN_DIR="%SDK_ROOT_CMAKE%/BSPLINEFORTRAN-7.4.0-%CFG%/lib" ^
  -DHDF5_DIR="%SDK_ROOT_CMAKE%/hdf5-1.12.2-%CFG%/cmake" ^
  -DCMAKE_C_FLAGS="/FS" ^
  -DCMAKE_CXX_FLAGS="/FS"
if errorlevel 1 (
  echo [ERROR] CMake configure failed for %CFG%.
  exit /b 1
)

echo [INFO] Building %CFG%: "%BUILD_PATH%"
cmake --build "%BUILD_PATH%"
if errorlevel 1 (
  echo [ERROR] Build failed for %CFG%.
  exit /b 1
)

exit /b 0
