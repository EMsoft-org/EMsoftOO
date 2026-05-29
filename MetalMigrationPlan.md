# Migrating EMsoftOO GPU code from OpenCL to Apple Metal

**Status:** Phase 0 in progress (on branch `feature/metal-backend`)
**Author:** drafted with Claude Code, 2026-05-29
**Decision summary:** Add a **Metal backend alongside** the existing OpenCL backend (do not
remove OpenCL — it is still needed for Windows/Linux and NVIDIA/AMD hardware). Metal
becomes the default GPU backend on Apple Silicon. Metal kernels are precompiled to
`.metallib` at **build time**. metal-cpp is **vendored in-tree**.

---

## Progress log

- **2026-05-29 — branch `feature/metal-backend` created** off `develop`.
- **metal-cpp vendored** at `ExternalProjects/metal-cpp/` (Apple release **macOS 15 / iOS 18**,
  Apache-2.0; 95 headers incl. Foundation/Metal/MetalFX/QuartzCore + `SingleHeader/`).
  This keeps the Metal build self-contained (no SDK/network dependency at configure time).
- **Phase 0 seam landed (OpenCL-only, no behavior change):**
  - Added GPU-operation wrapper methods to `OpenCL_T` in
    `Source/EMOpenCLLib/mod_CLsupport.f90`: `build_program`, `get_kernel`,
    `release_program`, `create_buffer`, `write_buffer`, `read_buffer`,
    `set_kernel_arg`, `enqueue_kernel`, `finish`, `release_buffer`, `release_kernel`,
    `release_context_queue`. `init_PDCCQ` now caches `context`/`command_queue`/device
    list inside the object so callers no longer pass them on every call.
  - **Migrated `mod_MCOpenCL.f90` (`MCOpenCL_`)** — the Monte Carlo PoC module — off the
    raw clfortran API and onto these wrappers (build/kernel, buffers, ~40 kernel args
    across the standard/foil/Ivol branches, dispatch, finish, reads, releases). No raw
    `clCreate*/clEnqueue*/clSetKernelArg/clRelease*` calls remain in that routine.
- **Built and verified** — compiles cleanly against the installed EMsoftOO_SDK (build dir
  `../EMsoftOObuild/MTL`). `EMMCOpenCL` was run on the same `EMMCOpenCL.nml` with the
  pre-refactor (`develop`) binary and the Phase 0 binary; `h5diff -v` on the `/EMData` group
  reports **0 differences** across every dataset (`accum_e`, `accum_z`, `accumSP`,
  `multiplier`, `numEbins`, `numzbins`, `totnum_el`) and the `HDF_FileVersion` attribute.
  The MC output is bit-identical, confirming the wrapper seam is behavior-preserving. (Only
  the timestamp in the metadata groups differs, as expected.)
- **2026-05-29 — `mod_DI.f90` migrated** onto the wrappers (commit 19551a5 covered the MC
  milestone; DI follows). Converted the shared `InnerProdGPU` helper (buffer/5 args/dispatch
  with explicit 16×16 local size/finish/read/release), the build+buffer blocks in all three
  GPU drivers (`DIdriver`, `OSMDIdriver`, `DIRAMdriver`), the `cl_expt`/`cl_dict` host writes,
  and every release (incl. `OSMDIdriver`'s queue/context via `release_context_queue`). No
  active raw `cl*` verb calls remain (only a pre-existing commented-out block). The CPU
  `sgemm` path (`DIRAMCPUdriver`) is untouched. **Needs build + `h5diff` verification** of an
  EMDI run vs a `develop` binary, same protocol as MC.
- **Remaining Phase 0 work:** `mod_EBSDFull.f90`, `mod_SEMCLwrappers.f90` (and confirm
  `mod_HROSM.f90`, which currently has no direct cl calls).

---

## 1. Motivation

Apple deprecated the OpenCL framework with macOS 10.14 and has signalled it will not be
maintained. On Apple Silicon (M-series) it currently runs only through a compatibility
layer, capped at OpenCL 1.2, with degrading reliability and performance, and is a
candidate for removal in a future macOS release. EMsoftOO's GPU-accelerated programs
(Monte Carlo, dictionary indexing, full-physics EBSD, multibeam scattering) should be
ported to Metal — Apple's native, actively developed GPU API — to stay viable on M3/M4
and later machines.

This is **future-proofing**, not an emergency fix: the OpenCL path still functions on
M3/M4 today. That gives us room to migrate incrementally and validate each step against
the existing OpenCL/CPU results.

---

## 2. Current GPU architecture (inventory)

### 2.1 Kernels — `opencl/` (6 files, 8 kernels, all single precision)

Verified: **no `double` precision anywhere** in the kernels. This matters because Apple
GPUs do **not** support fp64 in Metal — a port would have been blocked if any kernel
relied on doubles. All kernels use `float` / `float2`, so they are Metal-compatible as-is.

| File | Kernel(s) | Nature | Parallelism | Port difficulty |
|------|-----------|--------|-------------|-----------------|
| `EMMC.cl` | `MC` | Monte Carlo BSE electron trajectories (LFSR113 RNG + Lambert projection per thread) | Embarrassingly parallel, **no** barriers / local memory | Easy |
| `EMMCfoil.cl` | `MC` | MC for foil geometry (adds transmission outputs) | same | Easy |
| `EMMCxyz.cl` | `MCxyz` | MC with xyz direction-cosine output | same | Easy |
| `DictIndx.cl` | `InnerProd` | Tiled GEMM (BLOCK_SIZE×BLOCK_SIZE `__local` tiles + `barrier`) for dictionary dot products | Classic blocked matmul | Easy/medium — or replace (see §5.2) |
| `DictIndx.cl` | `ParamEstm` | Concentration-parameter estimation | Per-element | Medium |
| `MBmoduleOpenCL.cl` | `ScatMat`, `CalcLgh`, `CalcLghMaster` | Multibeam complex (`float2`) scattering-matrix propagation (~80 complex ops) | Matrix-heavy | Hardest |

Thread indexing uses `get_global_id(0/1)` (2-D grids). `InnerProd` is the only kernel
using `get_group_id` / `get_local_id` / `__local` / `barrier(CLK_LOCAL_MEM_FENCE)`.

### 2.2 Host side — `Source/EMOpenCLLib/`

- **`mod_CLsupport.f90`** (`OpenCL_T` class, ~1220 lines): wraps the `clfortran` Fortran
  bindings for platform/device enumeration, context + command-queue creation, source-file
  reading, memory estimation, and error checking.
- **Program modules** (`program_mods/`): `mod_MCOpenCL.f90`, `mod_DI.f90`,
  `mod_EBSDFull.f90`, `mod_SEMCLwrappers.f90`, `mod_HROSM.f90`.
- **Key coupling problem:** the program modules call the raw OpenCL C API **directly** —
  `clCreateBuffer`, `clCreateProgramWithSource`, `clBuildProgram`, `clCreateKernel`,
  `clSetKernelArg`, `clEnqueueWriteBuffer`, `clEnqueueNDRangeKernel`, `clEnqueueReadBuffer`
  — not through `OpenCL_T`. The GPU API surface is therefore **not** behind a single
  abstraction. This is the largest structural obstacle to a second backend.
- Kernels are loaded as **`.cl` source at runtime** (`read_source_file` → copied to
  `Bin/opencl/`) and JIT-compiled with `clBuildProgram`.

### 2.3 Build wiring

- `Source/Source.cmake`: `option(EMsoftOO_ENABLE_OpenCL_SUPPORT "Enable OpenCL support" ON)`
  gates `add_subdirectory(Source/EMOpenCLLib)`.
- `Source/EMOpenCLLib/CMakeLists.txt`: builds `EMOpenCLLib`, links `${OpenCL_LIBRARY}` +
  `clfortran`, includes `${CLFortran_INSTALL}/include`.
- `opencl/SourceList.cmake`: copies the 6 `.cl` files to `${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/opencl/`
  and installs them to `<install>/opencl/`.
- Per-modality `CMakeLists.txt` (DictionaryIndexing, OM, Demag, CTEMbook, …) gate targets
  on `EMsoftOO_ENABLE_OpenCL_SUPPORT` and link `EMOpenCLLib`.

---

## 3. Why Metal is not a drop-in replacement

1. **No Fortran bindings.** Metal's public API is Objective-C / Swift, or Apple's
   header-only C++ wrapper **metal-cpp**. There is no `metalfortran` equivalent to
   `clfortran`. Fortran must therefore call Metal through a **C-ABI shim** bound via
   `ISO_C_BINDING`.
2. **Different kernel language.** `.cl` (OpenCL C) → `.metal` (Metal Shading Language, a
   C++14 dialect). Translation is mostly mechanical (see §5.1), with one real semantic
   change: scalar kernel arguments.
3. **Different argument model.** OpenCL sets each scalar argument individually with
   `clSetKernelArg`. MSL prefers scalars **bundled into a `constant` struct buffer** bound
   to one buffer index. Each ported kernel gets a small parameter `struct` shared between
   the `.metal` file and the Fortran caller (via a matching derived type / C struct).
4. **Coexistence is required.** OpenCL stays for non-Apple platforms. We are *adding* a
   backend, which is the strongest argument for a clean abstraction so the kernels become
   the only duplicated artifact.

---

## 4. Target architecture

```
 Program modules (mod_MCOpenCL, mod_DI, mod_EBSDFull, …)
        │  call ONLY the backend-neutral interface
        ▼
 mod_GPUsupport.f90        ← new: abstract GPU device/buffer/kernel/dispatch API
        │
   ┌────┴───────────────┐
   ▼                    ▼
 mod_CLsupport.f90    mod_MTLsupport.f90      ← OpenCL impl (existing) / Metal impl (new)
   │ clfortran          │ ISO_C_BINDING
   ▼                    ▼
 OpenCL runtime       metal_shim  (metal-cpp, C++ → C ABI)
                         │
                         ▼
                       Metal.framework  +  *.metallib (built at build time)
```

### 4.1 Backend abstraction (`mod_GPUsupport.f90`, new)

Define a `GPU_T` class (or a thin façade) whose public methods cover exactly the API
surface the program modules use today:

- device/queue lifecycle: `init`, `select_device`, `destroy`
- program/kernel: `build_library`, `get_kernel(name)`
- memory: `create_buffer(size, access)`, `write_buffer`, `read_buffer`, `release_buffer`
- launch: `set_args(...)`, `dispatch(global, local)`, `finish`
- diagnostics: `error_check`, `device_info`, `memory_estimate`

The backend is selected at **compile time** by preprocessor flag (`EMSOFT_METAL` vs
`EMSOFT_OPENCL`), since a given binary targets one platform. `OpenCL_T` and `MTLsupport`
become the two concrete implementations; the existing `OpenCL_T` API is the template for
the interface, minimizing churn.

> **Prerequisite refactor (Phase 0):** move the raw `clCreate*/clEnqueue*/clSetKernelArg`
> calls *out of* the program modules and *into* the OpenCL implementation behind this
> interface. This has standalone value (cleaner code, single error-handling path) and is
> mandatory before a second backend can be slotted in cleanly.

### 4.2 Metal shim (`Source/EMOpenCLLib/metal/metal_shim.{hpp,cpp}`, new)

A small C++ translation unit using **metal-cpp** (header-only, single compiler — no
Objective-C runtime juggling) exposing a flat C ABI (~15–20 functions), e.g.:

```c
void*  emtl_create_device(void);
void*  emtl_create_queue(void* dev);
void*  emtl_load_library(void* dev, const char* metallib_path);   // build-time .metallib
void*  emtl_get_pipeline(void* lib, const char* fn_name);
void*  emtl_create_buffer(void* dev, size_t nbytes, int access);
void   emtl_write_buffer(void* buf, const void* src, size_t nbytes);
void   emtl_read_buffer (void* buf, void* dst, size_t nbytes);
void*  emtl_begin_dispatch(void* queue, void* pipeline);
void   emtl_set_buffer(void* enc, int index, void* buf);
void   emtl_dispatch(void* enc, /*grid*/uint64_t gx,uint64_t gy, /*tg*/uint64_t lx,uint64_t ly);
void   emtl_commit_wait(void* enc);
void   emtl_release(void* obj);
```

Fortran binds these with `ISO_C_BINDING` (`type(c_ptr)` handles, `bind(C)` interfaces) in
`mod_MTLsupport.f90`. On Apple Silicon's **unified memory**, buffers can be created with
`MTLResourceStorageModeShared`, so `write/read_buffer` are plain `memcpy`s with no PCIe
transfer — a natural performance win over discrete-GPU OpenCL.

### 4.3 Metal kernels (`metal/*.metal`, new)

One `.metal` file per current `.cl` file. See §5 for the translation playbook.

---

## 5. Kernel translation

### 5.1 Mechanical OpenCL C → MSL mapping

| OpenCL C | Metal Shading Language |
|----------|------------------------|
| `__kernel void K(...)` | `kernel void K(..., uint2 gid [[thread_position_in_grid]])` |
| `__global float* p` (arg n) | `device float* p [[buffer(n)]]` |
| `const float E` scalar (arg n) | field of a `constant Params& p [[buffer(n)]]` struct (see §3.3) |
| `__local float t[N][N]` | `threadgroup float t[N][N]` |
| `barrier(CLK_LOCAL_MEM_FENCE)` | `threadgroup_barrier(mem_flags::mem_threadgroup)` |
| `get_global_id(0/1)` | `gid.x` / `gid.y` |
| `get_local_id(0/1)` | `uint2 lid [[thread_position_in_threadgroup]]` |
| `get_group_id(0/1)` | `uint2 grp [[threadgroup_position_in_grid]]` |
| `float2` complex | `float2` (define `cmul`/`cadd` helpers; no built-in complex mul) |
| `sqrt/sin/cos/exp/...` | same names under `metal::` (usually unqualified is fine) |
| `rsqrt`, `mad`, `clamp` | available in `metal::` |

`float2` complex multiply must be written explicitly in both languages anyway; provide a
shared helper header for the multibeam kernels.

### 5.2 Per-kernel notes

- **MC / MCxyz / MCfoil (`EMMC*.cl`)** — easiest. Pure per-thread compute, 2-D grid, no
  shared memory. Port LFSR113 RNG and Lambert projection verbatim; bundle the ~12 scalar
  args (`E, count, z, rho, A, num_el, sig, omega, steps, …`) into one `Params` struct.
  **Recommended first proof-of-concept.**
- **`InnerProd` (`DictIndx.cl`)** — tiled GEMM with threadgroup tiles + barriers; translates
  directly. **But** consider replacing it on Apple with **`MPSMatrixMultiplication`**
  (Metal Performance Shaders) or even **Accelerate `cblas_sgemm`** (Apple AMX). The DI code
  *already has a CPU `sgemm` path* (`mod_DI.f90`), so on unified-memory M-series a library
  call may beat a hand-written kernel with far less code to maintain. Worth benchmarking
  before committing to a hand-translated kernel.
- **`ParamEstm` (`DictIndx.cl`)** — per-element with its own RNG/lookup; straightforward
  port after MC establishes the RNG pattern.
- **`ScatMat` / `CalcLgh` / `CalcLghMaster` (`MBmoduleOpenCL.cl`)** — hardest: heavy
  complex (`float2`) matrix work. Do last. Lean on **MPS** for the matmul-dominated parts;
  hand-translate the rest. Define complex helpers in a shared `.metal` header.

---

## 6. Build system changes (build-time `.metallib`)

1. **Option & detection** — in `Source/Source.cmake`:
   ```cmake
   option(EMsoftOO_ENABLE_Metal_SUPPORT "Enable Apple Metal GPU backend" ${APPLE})
   ```
   On Apple, default Metal ON. Guard so Metal is only ever enabled on `APPLE`.

2. **Frameworks** — in `EMOpenCLLib/CMakeLists.txt` (or a renamed `EMGPULib`):
   ```cmake
   if(EMsoftOO_ENABLE_Metal_SUPPORT)
     find_library(METAL_LIB Metal)
     find_library(FOUNDATION_LIB Foundation)
     find_library(QUARTZ_LIB QuartzCore)
     # metal-cpp headers (vendored under ExternalProjects/ or fetched)
     target_include_directories(EMOpenCLLib PRIVATE ${METAL_CPP_INCLUDE})
     target_link_libraries(EMOpenCLLib ${METAL_LIB} ${FOUNDATION_LIB} ${QUARTZ_LIB})
     set_source_files_properties(metal/metal_shim.cpp PROPERTIES
       COMPILE_FLAGS "-std=c++17 -fno-objc-arc")
   endif()
   ```

3. **Compile `.metal` → `.metallib` at build time** — a custom command per kernel file:
   ```cmake
   # for each kernel.metal:
   add_custom_command(OUTPUT ${OUT}/metal/EMMC.metallib
     COMMAND xcrun -sdk macosx metal   -c ${SRC}/metal/EMMC.metal -o ${TMP}/EMMC.air
     COMMAND xcrun -sdk macosx metallib   ${TMP}/EMMC.air -o ${OUT}/metal/EMMC.metallib
     DEPENDS ${SRC}/metal/EMMC.metal)
   ```
   Aggregate the `.metallib` outputs into a custom target that `EMOpenCLLib` depends on.
   Install them to `<install>/metal/` (parallel to the current `opencl/` install rule).
   The Metal `mod_MTLsupport` loads `metal/<name>.metallib` at runtime via
   `emtl_load_library` (instead of the OpenCL JIT-from-source flow).

4. **Per-modality gating** — change the many
   `if((EMsoftOO_ENABLE_HDF5_SUPPORT) AND (EMsoftOO_ENABLE_OpenCL_SUPPORT))` guards to
   `... AND (EMsoftOO_ENABLE_OpenCL_SUPPORT OR EMsoftOO_ENABLE_Metal_SUPPORT)`. Consider a
   single derived variable `EMsoftOO_ENABLE_GPU` set when either is on, to avoid editing
   each file's condition repeatedly.

5. **Library naming** — optional but recommended: the GPU library can be renamed from
   `EMOpenCLLib` to a backend-neutral `EMGPULib` once it hosts both backends. Defer to a
   later cleanup to keep the first diffs small; the export targets/install rules touch many
   places.

---

## 7. Validation strategy

For each ported kernel, add a numerical regression check comparing **Metal vs OpenCL vs
CPU** on a fixed input/seed:

- **MC kernels:** fixed RNG seeds → bit-for-bit or tight-tolerance match of the
  accumulated Lambert histograms / energy / depth arrays.
- **InnerProd / DI:** compare dot-product matrices against the existing CPU `sgemm` path
  (already present) to a tight relative tolerance.
- **Multibeam:** compare scattering-matrix / Lgh outputs against the OpenCL result on a
  small reference structure.

Run the existing programs end-to-end on a known `.nml` + reference dataset (from the
`EMsoftData` repo) and diff the HDF5 outputs. Because all kernels are single precision,
expect small floating-point differences from reordered reductions — use relative
tolerances, not exact equality, for the reductions; exact match where ordering is fixed.

---

## 8. Phasing (recommended order)

| Phase | Work | Risk |
|-------|------|------|
| **0** | Refactor: route program modules' raw `clCreate*/clEnqueue*` calls through a backend-neutral `mod_GPUsupport` interface; OpenCL stays the only implementation. Validate no behavior change. | Low, high value |
| **1** | Metal plumbing + MC proof-of-concept: `metal_shim` (metal-cpp), `mod_MTLsupport`, `EMMC.metal`, CMake `.metallib` rule, wire `mod_MCOpenCL` through the abstraction, validate. | Medium |
| **2** | DI: port `InnerProd`/`ParamEstm`, or replace `InnerProd` with MPS/Accelerate; benchmark. | Medium |
| **3** | Multibeam `ScatMat`/`CalcLgh`/`CalcLghMaster` (complex matrices, MPS-assisted). | High |
| **4** | Cleanup: rename `EMOpenCLLib`→`EMGPULib`, consolidate CMake GPU gating, docs. | Low |

Each phase is independently shippable and leaves OpenCL working throughout.

---

## 9. Risks & open questions

- **Dual-backend maintenance.** Two GPU code paths forever. Mitigated by the §4.1
  abstraction so only the `.metal`/`.cl` kernel pairs are duplicated, not host logic.
- **fp64:** none used today — but any *future* kernel must stay single precision on Metal,
  or fall back to CPU on Apple. Document this constraint for contributors.
- **metal-cpp sourcing:** vendor the header pack under `ExternalProjects/` or add to the
  EMsoftSuperbuild SDK (`DevelopOO` branch) so the build is self-contained. **Open question:
  add to the Superbuild SDK or vendor in-tree?**
- **Toolchain dependency:** build-time `.metallib` requires Xcode command-line tools
  (`xcrun metal`/`metallib`) on the build machine. Acceptable for a macOS build; document it.
- **MPS vs hand-written kernels:** decide per kernel by benchmark (DI especially). MPS
  reduces code but adds a framework dependency and less control over numerics.
- **Reductions / ordering:** GPU reduction order differs from CPU/OpenCL; validation must
  use tolerances. Identify any kernel relying on a specific accumulation order.

---

## 10. First concrete deliverable (when implementation starts)

Phase 0 + Phase 1 together: a backend abstraction, a working Metal path for the Monte
Carlo kernel end-to-end (shim → `mod_MTLsupport` → `EMMC.metal` → `.metallib` build rule),
and a validation run of `EMMCOpenCL` comparing Metal vs OpenCL output on a reference
`.nml`. That proves the entire toolchain (Fortran↔C↔metal-cpp↔Metal, build-time
`.metallib`, unified-memory buffers) on the lowest-risk kernel before tackling DI and
multibeam.
