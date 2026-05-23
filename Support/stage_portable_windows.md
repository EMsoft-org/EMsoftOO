# Portable Export on Windows

This is the manual procedure for making a portable EMsoftOO runtime that can be moved to another Windows machine.

## 1. Create a staging folder

Create a folder such as:

```text
C:\Users\westraadt.1\Desktop\EMsoftOO-portable
```

Use this layout:

```text
EMsoftOO-portable/
  bin/
  NamelistTemplates/
  resources/
  opencl/
  XtalFolder/
  InputData/
```

## 2. Copy the runtime files

Copy all files from:

```text
C:\Users\westraadt.1\EMSOFT\EMsoftOOBuild\Release\Bin\
```

into:

```text
EMsoftOO-portable\bin\
```

Do not select individual executables. Copy the complete contents of `Bin`.

## 3. Copy the repo data folders

Copy these folders from the repo:

```text
C:\Users\westraadt.1\Repos\EMsoftOO\NamelistTemplates\
C:\Users\westraadt.1\Repos\EMsoftOO\resources\
C:\Users\westraadt.1\Repos\EMsoftOO\opencl\
```

into:

```text
EMsoftOO-portable\NamelistTemplates\
EMsoftOO-portable\resources\
EMsoftOO-portable\opencl\
```

Copy the full contents of each folder.

Important:

- `resources\ShapeFiles\` must be included
- `opencl\*.cl` files must be included

## 4. Copy extra DLLs from the SDK

If these DLLs are not already present in `bin`, copy them into:

```text
EMsoftOO-portable\bin\
```

from:

```text
C:\Users\westraadt.1\EMSOFT\EMsoftOO_SDK\nlopt-2.7.0-Release\bin\nlopt.dll
C:\Users\westraadt.1\EMSOFT\EMsoftOO_SDK\bcls-0.1-Release\bin\bcls.dll
C:\Users\westraadt.1\EMSOFT\EMsoftOO_SDK\tbb-2020.1-win\tbb\bin\intel64\vc14\tbb.dll
C:\Users\westraadt.1\EMSOFT\EMsoftOO_SDK\tbb-2020.1-win\tbb\bin\intel64\vc14\tbbmalloc.dll
```

You do not need to copy the full `EMsoftOO_SDK`.

## 5. Copy your input data

Add your own files:

- Put `.xtal` files in `XtalFolder\`
- Put HDF5 files, master patterns, and other large inputs in `InputData\`
- You can also keep input files elsewhere, but then the config and namelists must point to those paths

## 6. Zip the folder

Zip the full `EMsoftOO-portable` folder and move it to the target machine.

## 7. Set it up on the target machine

Unzip it to a permanent location, for example:

```text
C:\Tools\EMsoftOO-portable
```

Open PowerShell and set:

```powershell
$env:EMSOFTPATHNAME = 'C:\Tools\EMsoftOO-portable\'
```

Then run:

```powershell
C:\Tools\EMsoftOO-portable\bin\EMsoftinit.exe
```

This creates:

```text
%USERPROFILE%\.config\EMsoft\EMsoftConfig.json
```

Edit that file and set at least:

- `EMsoftpathname = C:/Tools/EMsoftOO-portable/`
- `EMXtalFolderpathname = C:/Tools/EMsoftOO-portable/XtalFolder/`
- `EMdatapathname =` a writable output folder
- `EMtmppathname =` a writable temp folder
- `Release = Yes`
- `Develop = No`

Make sure path values end with `/` or `\`.

## 8. Notes

- The target machine still needs a working GPU OpenCL driver if you want to run `EMEBSDFull` or other OpenCL programs.
- The portable folder must contain `NamelistTemplates`, `resources`, and `opencl` next to `bin`. Do not place them inside `bin`.
- If a run references additional input files, those must also be copied to the target machine.
