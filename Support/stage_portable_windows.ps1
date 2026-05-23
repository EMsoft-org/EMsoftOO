[CmdletBinding()]
param(
  [string]$RepoRoot = 'C:\Users\westraadt.1\Repos\EMsoftOO',
  [string]$BuildRoot = 'C:\Users\westraadt.1\EMSOFT\EMsoftOOBuild\Release',
  [string]$StageRoot = 'C:\Users\westraadt.1\EMSOFT\EMsoftOO-portable',
  [string]$SDKRoot = 'C:\Users\westraadt.1\EMSOFT\EMsoftOO_SDK',
  [string]$XtalFolder = '',
  [string[]]$InputData = @(),
  [string]$ZipPath = '',
  [switch]$Clean
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

function Resolve-ExistingPath {
  param(
    [Parameter(Mandatory = $true)]
    [string]$Path,
    [Parameter(Mandatory = $true)]
    [string]$Label
  )

  if (-not (Test-Path -LiteralPath $Path)) {
    throw "$Label not found: $Path"
  }

  return (Resolve-Path -LiteralPath $Path).Path
}

function Ensure-Directory {
  param(
    [Parameter(Mandatory = $true)]
    [string]$Path
  )

  if (-not (Test-Path -LiteralPath $Path)) {
    New-Item -ItemType Directory -Path $Path | Out-Null
  }
}

function Copy-DirectoryContents {
  param(
    [Parameter(Mandatory = $true)]
    [string]$Source,
    [Parameter(Mandatory = $true)]
    [string]$Destination
  )

  Ensure-Directory -Path $Destination
  Copy-Item -LiteralPath (Join-Path $Source '*') -Destination $Destination -Recurse -Force
}

function Copy-OptionalFile {
  param(
    [Parameter(Mandatory = $true)]
    [string]$Source,
    [Parameter(Mandatory = $true)]
    [string]$DestinationDirectory
  )

  if (Test-Path -LiteralPath $Source) {
    Copy-Item -LiteralPath $Source -Destination $DestinationDirectory -Force
    return $true
  }

  return $false
}

function Remove-StageRoot {
  param(
    [Parameter(Mandatory = $true)]
    [string]$Path
  )

  $resolved = Resolve-Path -LiteralPath $Path
  $target = $resolved.Path

  if ($target -match '^[A-Za-z]:\\?$') {
    throw "Refusing to remove a drive root: $target"
  }

  Remove-Item -LiteralPath $target -Recurse -Force
}

$RepoRoot = Resolve-ExistingPath -Path $RepoRoot -Label 'RepoRoot'
$BuildRoot = Resolve-ExistingPath -Path $BuildRoot -Label 'BuildRoot'
$BuildBin = Resolve-ExistingPath -Path (Join-Path $BuildRoot 'Bin') -Label 'BuildRoot\Bin'
$SDKRoot = Resolve-ExistingPath -Path $SDKRoot -Label 'SDKRoot'

if (Test-Path -LiteralPath $StageRoot) {
  if (-not $Clean) {
    throw "StageRoot already exists: $StageRoot. Re-run with -Clean to replace it."
  }

  Remove-StageRoot -Path $StageRoot
}

Ensure-Directory -Path $StageRoot

$StageBin = Join-Path $StageRoot 'bin'
$StageTemplates = Join-Path $StageRoot 'NamelistTemplates'
$StageResources = Join-Path $StageRoot 'resources'
$StageOpenCL = Join-Path $StageRoot 'opencl'
$StageXtal = Join-Path $StageRoot 'XtalFolder'
$StageInputData = Join-Path $StageRoot 'InputData'
$StageDocs = Join-Path $StageRoot 'Documentation'

Ensure-Directory -Path $StageBin
Ensure-Directory -Path $StageTemplates
Ensure-Directory -Path $StageResources
Ensure-Directory -Path $StageOpenCL
Ensure-Directory -Path $StageXtal
Ensure-Directory -Path $StageInputData
Ensure-Directory -Path $StageDocs

# Copy the complete runtime output rather than trying to enumerate binaries.
Copy-DirectoryContents -Source $BuildBin -Destination $StageBin
Copy-DirectoryContents -Source (Join-Path $RepoRoot 'NamelistTemplates') -Destination $StageTemplates
Copy-DirectoryContents -Source (Join-Path $RepoRoot 'resources') -Destination $StageResources
Copy-DirectoryContents -Source (Join-Path $RepoRoot 'opencl') -Destination $StageOpenCL

Copy-Item -LiteralPath (Join-Path $RepoRoot 'README.md') -Destination $StageDocs -Force
Copy-Item -LiteralPath (Join-Path $RepoRoot 'License.txt') -Destination $StageDocs -Force

if ($XtalFolder -ne '') {
  $ResolvedXtalFolder = Resolve-ExistingPath -Path $XtalFolder -Label 'XtalFolder'
  Copy-DirectoryContents -Source $ResolvedXtalFolder -Destination $StageXtal
}

foreach ($InputPath in $InputData) {
  $ResolvedInput = Resolve-ExistingPath -Path $InputPath -Label 'InputData entry'
  $LeafName = Split-Path -Path $ResolvedInput -Leaf

  if (Test-Path -LiteralPath $ResolvedInput -PathType Container) {
    Copy-DirectoryContents -Source $ResolvedInput -Destination (Join-Path $StageInputData $LeafName)
  }
  else {
    Copy-Item -LiteralPath $ResolvedInput -Destination $StageInputData -Force
  }
}

# These DLLs are not guaranteed to be staged by the current CMake install rules.
$RuntimeDlls = @(
  (Join-Path $SDKRoot 'nlopt-2.7.0-Release\bin\nlopt.dll'),
  (Join-Path $SDKRoot 'bcls-0.1-Release\bin\bcls.dll'),
  (Join-Path $SDKRoot 'tbb-2020.1-win\tbb\bin\intel64\vc14\tbb.dll'),
  (Join-Path $SDKRoot 'tbb-2020.1-win\tbb\bin\intel64\vc14\tbbmalloc.dll')
)

$CopiedOptionalDlls = @()
foreach ($DllPath in $RuntimeDlls) {
  if (Copy-OptionalFile -Source $DllPath -DestinationDirectory $StageBin) {
    $CopiedOptionalDlls += (Split-Path -Path $DllPath -Leaf)
  }
}

$PortableReadme = @"
EMsoftOO Portable Runtime
=========================

Folder layout:
- bin\ contains the executables and runtime DLLs.
- NamelistTemplates\, resources\, and opencl\ are required at runtime.
- XtalFolder\ is the default location for crystal files.
- InputData\ is an optional place for large HDF5/master-pattern inputs.

Setup on the target machine:
1. Unzip this folder to a final location, for example C:\Tools\EMsoftOO-portable.
2. Open a PowerShell window.
3. Set EMSOFTPATHNAME to the unzip root:
   `$env:EMSOFTPATHNAME = 'C:\Tools\EMsoftOO-portable\'
4. Run:
   .\bin\EMsoftinit.exe
5. Edit:
   %USERPROFILE%\.config\EMsoft\EMsoftConfig.json
6. Set at minimum:
   - EMsoftpathname = C:/Tools/EMsoftOO-portable/
   - EMXtalFolderpathname = C:/Tools/EMsoftOO-portable/XtalFolder/
   - EMdatapathname = a writable output folder
   - EMtmppathname = a writable temp folder

Notes:
- EMEBSDFull and other OpenCL programs still require a working GPU OpenCL driver on the target machine.
- Oxford .h5oina files that depend on an HDF5 LZF plugin are a separate case; this package does not add an HDF5 plugin directory.
"@

Set-Content -LiteralPath (Join-Path $StageRoot 'README-portable.txt') -Value $PortableReadme -Encoding ascii

if ($ZipPath -ne '') {
  $ZipDirectory = Split-Path -Path $ZipPath -Parent
  if ($ZipDirectory -ne '') {
    Ensure-Directory -Path $ZipDirectory
  }

  if (Test-Path -LiteralPath $ZipPath) {
    Remove-Item -LiteralPath $ZipPath -Force
  }

  Compress-Archive -Path (Join-Path $StageRoot '*') -DestinationPath $ZipPath -Force
}

Write-Host "Portable runtime staged at: $StageRoot"
Write-Host "Build bin copied from:      $BuildBin"
if ($CopiedOptionalDlls.Count -gt 0) {
  Write-Host "Extra SDK DLLs copied:      $($CopiedOptionalDlls -join ', ')"
}
if ($ZipPath -ne '') {
  Write-Host "Zip created at:             $ZipPath"
}
