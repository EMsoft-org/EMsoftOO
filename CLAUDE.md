# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

EMsoftOO (version 6.0) is a complete rewrite of the EMsoft electron microscopy simulation package in object-oriented Fortran 2018. It simulates diffraction patterns and performs indexing for various electron/X-ray microscopy modalities (EBSD, ECP, TKD, TEM, Laue, etc.). The code is BSD-3 licensed and developed at Carnegie Mellon University (Marc De Graef Research Group).

> Machine-specific notes (where the SDK lives, the local build directory, the Metal toolchain, and how builds are driven on this machine) are kept in `CLAUDE.local.md`, which is gitignored. This file holds only the project-wide, shareable guidance.

## Build System

EMsoftOO uses CMake and requires a pre-built SDK (`EMsoftOO_SDK`) from the [EMsoftSuperbuild](https://github.com/EMsoft-org/EMsoftSuperbuild) repository (use the `DevelopOO` branch). It also requires the `EMsoftData` repository cloned at the same directory level as `EMsoftOO`.

### Standard build

```bash
mkdir -p EMsoftOOBuild/Release
cd EMsoftOOBuild/Release
cmake -DCMAKE_BUILD_TYPE=Release -DEMsoftOO_SDK=/path/to/EMsoftOO_SDK ../../EMsoftOO
make -j

# Debug build (slower but better error messages)
mkdir -p ../Debug
cd ../Debug
cmake -DCMAKE_BUILD_TYPE=Debug -DEMsoftOO_SDK=/path/to/EMsoftOO_SDK ../../EMsoftOO
make -j
```

Binaries are placed in `EMsoftOOBuild/Release/Bin/`.

### Selective compilation

Individual modalities can be toggled via CMake options. From the build directory:

```bash
ccmake ../../EMsoftOO
```

Key options (toggle ON/OFF):
- `EMsoftOO_ENABLE_SEM` — EBSD, ECP, TKD, and related SEM programs
- `EMsoftOO_ENABLE_DictionaryIndexing` — Dictionary indexing programs
- `EMsoftOO_ENABLE_TEM` — TEM programs
- `EMsoftOO_ENABLE_Utilities` — Utility programs
- `EMsoftOO_ENABLE_HDF5_SUPPORT` — HDF5 I/O (required for most programs)
- `EMsoftOO_ENABLE_OpenCL_SUPPORT` — GPU acceleration via OpenCL
- `EMsoftOO_ENABLE_XRay`, `EMsoftOO_ENABLE_GBs`, `EMsoftOO_ENABLE_QC`, etc.

### Running tests

Tests are disabled by default. Enable with:

```bash
cmake -DEMsoftOO_ENABLE_TESTING=ON ...
make -j
ctest
```

Unit test sources are in `Source/Test/` (e.g., `MODRotationsTest.f90`, `MODQuaternionsTest.f90`).

## Running Programs

All executables take a single namelist (`.nml`) file as the argument:

```bash
EMEBSDmaster EMEBSDmaster.nml
```

Template `.nml` files for every program live in `NamelistTemplates/`. Copy and modify them to set up a run. The `EMsoft_T` constructor reads the nml file path from the command line.

## Code Architecture

### Two-layer program structure

Every simulation program has two files:

1. **Thin driver** in `Source/<Modality>/EM<ProgramName>.f90` — instantiates a few objects and calls the compute method. These are ~50 lines.
2. **Implementation module** in `Source/EMsoftOOLib/program_mods/mod_<ProgramName>.f90` — contains the `<ProgramName>_T` class with the full computation logic, namelist I/O, and HDF5 output.

Example: `Source/SEM/EMEBSDmaster.f90` creates `EBSDmaster_T` and calls `MP%EBSDmaster(...)`, which is implemented in `Source/EMsoftOOLib/program_mods/mod_EBSDmaster.f90`.

### Core library: `Source/EMsoftOOLib/`

The main library (`EMsoftOOLib`) is compiled first and linked by all executables. Key foundational modules (in dependency order):

- `mod_kinds.f90` — kind parameters: `sgl` (real32), `dbl` (real64), `irg` (int32), `ill` (int64), `ish` (int16)
- `mod_global.f90` — global constants, `fnlen` string length parameter
- `mod_memory.f90` — `memory_T` class: wraps `allocate`/`deallocate` with bounds checking and initialization; use `mem%alloc(arr, shape, 'name')` / `mem%dealloc(arr, 'name')` instead of raw Fortran allocate
- `mod_io.f90` — `IO_T` class: all terminal I/O goes through this; replaces direct `write`/`read` calls
- `mod_EMsoft.in.f90` (generated as `mod_EMsoft.f90`) — `EMsoft_T` class: program entry point, config management, command-line argument handling. Every program starts with `EMsoft = EMsoft_T(progname, progdesc, tpl = (/...,0/))`
- `mod_HDFsupport.f90` — `HDF_T` class: all HDF5 file I/O
- `mod_HDFnames.f90` — `HDFnames_T` class: standardizes HDF5 dataset/group name strings via `stringconstants`
- `mod_rotations.f90` — `Rotations_T` class: rotation representations (eu, om, ax, ro, qu, ho, cu, st, rv) and conversions between them
- `mod_quaternions.f90` — quaternion arithmetic
- `mod_so3.f90` — SO(3) sampling, Rodrigues fundamental zones
- `mod_Lambert.f90` — Lambert/cubochoric sphere-plane mappings
- `mod_crystallography.f90` — `Crystal_T` class: crystal structure, lattice parameters
- `mod_symmetry.f90` — space group / point group symmetry operations
- `mod_diffraction.f90` — electron diffraction structure factors
- `mod_patterns.f90` — detector/pattern generation utilities
- `mod_DIsupport.f90` — dictionary indexing support routines
- `mod_JSONsupport.f90` — JSON config file reading (EMsoftConfig.json)

### Modality source directories under `Source/`

| Directory | Contents |
|-----------|----------|
| `SEM/` | EBSD, ECP, TKD, ECCI, ISE, HREBSD programs |
| `TEM/` | CBED, LACBED, HOLZ, defect, STEM-DCI programs |
| `DictionaryIndexing/` | EBSD/ECP/TKD dictionary indexing (EMDI) |
| `SphInx/` | Spherical indexing programs |
| `GBs/` | Grain boundary programs |
| `QC/` | Quasicrystal programs |
| `XRay/` | Laue diffraction, XRD programs |
| `Utilities/` | Orientation tools, format converters, crystal tools |
| `CTEMbook/` | Programs from the CTEM textbook |
| `Demag/` | Demagnetization programs |
| `Shapes/` | Shape amplitude programs |
| `OM/` | Optical microscopy programs |
| `Test/` | Unit test sources |

### String constants

`stringconstants` module (generated from `Source/EMsoftOOLib/stringconstants.in.f90`) provides `SC_*` constants for all HDF5 group/dataset name strings. Use these instead of hardcoded strings.

### Generated files

Two source files are generated by CMake at configure time (do not edit directly — edit the `.in.f90` templates):
- `mod_EMsoft.f90` from `mod_EMsoft.in.f90`
- `mod_platformsupport.f90` from `mod_platformsupport.in.f90`

### GPU programs (OpenCL / Apple Metal)

GPU-accelerated programs (e.g., `EMMCGPU`, `EMDI`) use the separate `EMOpenCLLib` library in `Source/EMOpenCLLib/`. That library exposes a single GPU abstraction — `mod_GPUsupport` / `GPU_T` — with two interchangeable implementations selected at build time by `EMsoftOO_ENABLE_Metal_SUPPORT`: the OpenCL backend (`mod_GPUsupport.f90`) and the Apple Metal backend (`mod_GPUsupport_metal.f90` + the metal-cpp C shim in `metal/emtl_shim.{h,cpp}`). On macOS, Metal is the default. See `MetalMigrationPlan.md`.

Two programs were renamed once the GPU path stopped being OpenCL-specific (the old names are still built as aliases): `EMMCOpenCL` → **`EMMCGPU`** (Monte Carlo; `Source/SEM/EMMCGPU.f90`) and `EMOpenCLinfo` → **`EMGPUinfo`** (GPU/device info; `Source/Utilities/EMGPUinfo.f90`). Note: `EMMCGPU` still writes `'EMMCOpenCL.f90'` as the program identity into its Monte Carlo HDF5 output (`MCfileProgName` in the driver) so that existing files and the downstream master-pattern readers (`mod_HDFFileInfo`, EBSD/ECP/TKD master) stay compatible — do not change that embedded string.

## Fortran Coding Conventions

- All modules are named `mod_<name>` and define a class `<Name>_T` with private data and public methods
- Procedure bindings use the pattern `procedure, pass(self) :: method_` (trailing underscore for private name), exposed publicly via a generic interface without the underscore
- Use `mem%alloc` / `mem%dealloc` (from `mod_memory`) for array management instead of `allocate`/`deallocate`
- Use `call Message%printMessage(...)` or `call Message%printError(...)` (from `mod_io`) for all console output
- Always `use mod_kinds` and `use mod_global` in every program unit; declare `IMPLICIT NONE`
- All new source files must include the BSD-3 license header block at the top
