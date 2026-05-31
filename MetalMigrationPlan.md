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
  `sgemm` path (`DIRAMCPUdriver`) is untouched.
- **`mod_DI` verified** — EMDI run, develop vs Phase 0 binaries on the same fixed Euler-angle
  dictionary, `h5diff` on `/Scan 1`: `TopDotProductList`, `TopMatchIndices`, `EulerAngles`,
  `CI`, `Phi/Phi1/Phi2`, `KAM`, `OSM` are all **0 differences** — the direct `InnerProd` GPU
  output and every orientation result are bit-identical, confirming the refactor is
  behavior-preserving. Four datasets still differ (`DictionaryEulerAngles`, `ISM`, `ISMap`,
  `IndexingSuccessRate`); these are *not* downstream of the (identical) dot products and trace
  to pre-existing non-determinism in the FZ reduction (`mod_so3`/`mod_sampleRFZ`, which carry
  unrelated uncommitted edits) — a separate known issue, not caused by Phase 0.
- **2026-05-29 — `mod_EBSDFull.f90` migrated** onto the wrappers. Its GPU section is the
  embedded Monte Carlo simulation in `ComputeFullEBSDPatterns_` (`EMMC.cl`, the same 14-arg
  `MC` kernel as `mod_MCOpenCL`'s full mode — already verified bit-identical). Converted the
  build, 5 buffers, seed write, 14 kernel args, dispatch/finish/4 reads, and all releases
  (`release_context_queue` for queue+context). The CPU dynamical-pattern code
  (`ComputeFullDynamicalPatterns`/`CalcLghSM`) is untouched. No active raw `cl*` verbs remain.
  **Pending build + verification** (EMEBSDFull); low risk as the kernel/args are identical to
  the verified MC path.
- **2026-05-29 — wrappers gained an optional `quiet` flag.** `mod_CLsupport`'s GPU-op
  wrappers now take `quiet` (via a small `checkq_` helper): when `.TRUE.` they skip the fatal
  `error_check_` (and, for `build_program`, the build-log print). Default is unchanged
  (fatal checking), so the already-verified `mod_MCOpenCL`/`mod_DI`/`mod_EBSDFull` call paths
  are byte-for-byte unaffected.
- **2026-05-29 — `mod_SEMCLwrappers.f90` migrated** (the C-callable `EMsoftCgetMCOpenCL`, an
  embedded MC sim using the same `MC` kernel/14 args). Routed through the wrappers with
  `quiet=.TRUE.` throughout to preserve its original "ignore CL errors, defer to the caller"
  semantics — important because it is invoked from external host programs and `error_check_`
  would otherwise `stop` the process. No active raw `cl*` verbs remain.
- **`mod_HROSM.f90` confirmed** — no GPU code (only matched the scan via `MCOpenCL`
  HDF group-name constants); nothing to migrate.
- **Phase 0 active scope complete.** All compiled GPU consumers
  (`mod_MCOpenCL`, `mod_DI`, `mod_EBSDFull`, `mod_SEMCLwrappers`) route through `OpenCL_T`'s
  wrappers; `mod_HROSM` has no GPU code.
- **Deferred / TODO:** `mod_DIPCA.f90` is **not compiled** (absent from every CMakeLists — a
  dead WIP PCA-DI variant) and still contains raw `clfortran` calls. It must be migrated *if
  it is ever re-enabled*; left untouched now since it cannot be built or verified.

### Phase 1 — Metal backend (in progress, on `feature/metal-backend`)

- **`Source/EMOpenCLLib/metal/EMMC.metal`** — MSL translation of `opencl/EMMC.cl` (the MC
  kernel). Behaviour-preserving: LFSR113 RNG, Lambert projection and MC physics reproduced
  exactly; all single precision (Apple GPUs have no fp64). Buffer indices match the OpenCL
  `clSetKernelArg` indices so the *same* host code drives either backend. Can be compiled
  standalone for an early check: `xcrun -sdk macosx metal -c EMMC.metal -o EMMC.air`.
- **`Source/EMOpenCLLib/metal/emtl_shim.h` + `emtl_shim.cpp`** — C-ABI shim over metal-cpp
  (the vendored headers). Key design choices:
  - Handles are `int64_t` (= metal-cpp pointers cast), so the Fortran side keeps using
    `integer(c_intptr_t)` handles exactly as for OpenCL.
  - Buffers use `MTLResourceStorageModeShared` (unified memory) — write/read are plain
    `memcpy` to/from `buffer->contents()`.
  - **Arg discrimination:** `emtl_set_arg(index, ptr, size)` checks a live-buffer registry —
    if `*ptr` is a registered buffer handle it binds with `setBuffer`, else `setBytes`. This
    lets the host call one uniform `set_kernel_arg(index, size, ptr)` for both buffers and
    scalars, exactly like OpenCL — so **no program-module changes** are needed.
  - **Set-once / dispatch-many:** the OpenCL host sets args once then enqueues in a loop;
    Metal binds per-encoder, so the shim *caches* arg bindings per pipeline and re-applies
    them inside `emtl_enqueue` (which builds a fresh command buffer/encoder, commits without
    waiting). `emtl_finish` waits on the last commit — mirroring `clEnqueueNDRangeKernel` +
    `clFinish`. Dispatch supports auto threadgroup (`dispatchThreads`) and explicit local
    size (`dispatchThreadgroups`, for the tiled DI kernel in Phase 2).
- **Backend-selection decision (chosen): CMake source-swap.** The Metal backend will live in
  `Source/EMOpenCLLib/mod_CLsupport_metal.f90`, which defines the **same** module name
  (`mod_CLsupport`) and type (`OpenCL_T`) with the **same public method surface** as the
  OpenCL `mod_CLsupport.f90`, but backed by `emtl_shim`. CMake compiles **one** of the two
  files depending on `EMsoftOO_ENABLE_Metal_SUPPORT`. Result: **zero changes** to the program
  modules (`mod_MCOpenCL` etc. keep `use mod_CLsupport` / `type(OpenCL_T)`), no Fortran
  preprocessor needed. The `OpenCL_*` names remain under Metal until the Phase 4 rename to
  `GPU_T`/`EMGPULib`. (Considered and rejected for the PoC: `#ifdef` preprocessor — needs
  `-cpp`; and an immediate `GPU_T` rename — re-touches verified Phase 0 code.)
- **Still to do this phase:** `mod_CLsupport_metal.f90` (drop-in `OpenCL_T` over `emtl_shim`,
  mirroring every public method incl. the `quiet` flag; `read_source_file` derives the
  `.metallib` path, `build_program` loads it, `get_kernel` builds the pipeline,
  `enqueue_kernel` maps global/local sizes); CMake (Metal option, compile `emtl_shim.cpp` as
  C++17 with the vendored metal-cpp includes, link `-framework Metal -framework Foundation
  -framework QuartzCore`, build `EMMC.metal`→`EMMC.metallib` at build time and install to
  `Bin/metal/`, swap the Fortran backend source); then run `EMMCOpenCL` on Metal and
  `h5diff` against `reference.h5`.
- **Standalone checks (DONE — both compile clean):**
  `xcrun -sdk macosx metal -c .../EMMC.metal` validated the kernel; `clang++ -std=c++17
  -I ExternalProjects/metal-cpp -c .../emtl_shim.cpp` validated the shim against metal-cpp.
- **2026-05-29 — Phase 1 implementation landed (CMake source-swap backend):**
  - `Source/EMOpenCLLib/mod_CLsupport_metal.f90` — drop-in `module mod_CLsupport` /
    `type OpenCL_T` over `emtl_shim`. Mirrors every public method of the OpenCL backend
    (incl. the `quiet` flag). `read_source_file` maps the requested `*.cl` name to the
    prebuilt `*.metallib` (resolved via `OpenCLpathname`); `build_program` loads the
    library; `get_kernel` builds the compute pipeline; `enqueue_kernel` maps global/local
    sizes; buffers/args/finish/release go through the shim. Handles stay
    `integer(c_intptr_t)`.
  - `Source/EMOpenCLLib/clfortran_metal_stub.f90` — a tiny `module clfortran` providing just
    `CL_MEM_READ_WRITE/WRITE_ONLY/READ_ONLY` (+ `CL_TRUE/FALSE`), compiled only in a Metal
    build so the program modules' `use clfortran` resolves without the real library. (Audit
    confirmed those are the only active clfortran symbols the program modules reference.)
  - **CMake:** `EMsoftOO_ENABLE_Metal_SUPPORT` option (OFF by default, Apple-only);
    `Source.cmake` widens the EMOpenCLLib gate to `OpenCL OR Metal`; `EMOpenCLLib/CMakeLists.txt`
    swaps the backend source set, compiles `emtl_shim.cpp` as C++17 with the vendored
    metal-cpp include, links `-framework Metal/Foundation/QuartzCore` + `c++`, and builds each
    `*.metal`→`*.metallib` at build time into `Bin/opencl/` (next to the `.cl` files, so the
    runtime path resolution is reused). Only `EMMC` is enabled in the kernel list (DI/MB are
    Phase 2/3).
  - **No program-module changes** — `mod_MCOpenCL` etc. are untouched; the source-swap makes
    `mod_CLsupport`/`OpenCL_T` resolve to the Metal implementation.
- **2026-05-29 — Metal MC PoC BUILT and VALIDATED.** `cmake -DEMsoftOO_ENABLE_Metal_SUPPORT=ON`
  + `make` built clean (MSL→metallib, C++ shim, Fortran backend, stub, mixed C++/Fortran link
  all OK). `EMMCOpenCL` ran on Metal; `h5diff /EMData` vs the OpenCL `reference.h5`:
  - `accumSP` (the energy-summed master pattern that feeds EBSD): **0 differences** (bit-identical).
  - `accum_e`: 4126 bins differ, but only **4** by more than 1 count and **none** by more than 2.
  - `accum_z`: 2733 differ, **25** by more than 1, **none** by more than 2.
  Every difference is essentially a single electron (±1) crossing a histogram bin edge; no
  systematic shift. This is inherent chaotic-MC floating-point divergence: the kernel's
  transcendentals (`pow`/`log`/`acos`/`sin`/…) and FMA contraction differ ~1 ULP between the
  Metal and OpenCL compilers, and over ~300 random-walk steps that occasionally nudges an
  electron across a bin. Bit-identity is unattainable for a chaotic MC across two different
  GPU compilers (an OpenCL run on a different GPU vendor would jitter the same way); `accumSP`
  averaging confirms the physics is equivalent. **The Metal backend is proven end-to-end.**
- **Runtime path fix:** the dev `OpenCLpathname` resolves to the source-tree `opencl/` folder
  (where the `.cl` files are version-controlled), so CMake now also copies each built
  `*.metallib` there (in addition to `Bin/opencl/`); `opencl/*.metallib` is gitignored.
- **Phase 1 complete** (committed `666086a`). Remaining: Phase 2 (DI), Phase 3 (multibeam),
  Phase 4 (rename to `GPU_T`/`EMGPULib`, metallib install rules, default Metal ON on Apple).
  Validation note for Phase 2+: bit-identity only where the computation is linear (DI dot
  products); tolerance/statistical comparison for chaotic or reduction-order-sensitive output.
- **Footnote (observed):** on Apple Silicon the Metal MC ran ~2x faster than the OpenCL MC.
  Expected — Apple OpenCL is a deprecated non-native 1.2 compatibility layer, while Metal is
  native with a modern compiler, hardware-tuned scheduling, true unified memory (our
  `StorageModeShared` buffers need no host<->device copy), and build-time `.metallib` (no JIT).

### Phase 2 — dictionary indexing (in progress)

- **`Source/EMOpenCLLib/metal/DictIndx.metal`** — MSL port of the `InnerProd` tiled GEMM
  (BLOCK_SIZE=16). `get_group_id`/`get_local_id`/`get_global_id` →
  `[[threadgroup_position_in_grid]]`/`[[thread_position_in_threadgroup]]`/`[[thread_position_in_grid]]`;
  `get_global_size(0)` → `[[threads_per_grid]].x`; `__local`→`threadgroup`;
  `barrier(CLK_LOCAL_MEM_FENCE)`→`threadgroup_barrier(mem_flags::mem_threadgroup)`. Arg indices
  match (0 expt,1 dict,2 Wexp,3 Wdict,4 result). The (16,16) local size routes through
  `emtl_enqueue`'s `dispatchThreadgroups` path. `ParamEstm` omitted (unused by any compiled module).
- **CMake:** `DictIndx` added to `EMsoftOO_METAL_KERNELS` → built to `DictIndx.metallib`.
- **No other changes:** `mod_DI` was wrapper-routed in Phase 0, so everything else flows through
  the existing Metal backend (read_source_file→DictIndx.metallib, get_kernel('InnerProd'),
  expt/dict/result auto-detected as buffers, Wexp/Wdict as scalar bytes, (16,16) dispatch).
- **2026-05-29 — Phase 2 VALIDATED.** Rebuilt Metal ON, ran `EMDI`, `h5diff /Scan 1` vs the
  OpenCL reference: `TopDotProductList`, `TopMatchIndices`, `EulerAngles`, `CI`, `Phi/Phi1/Phi2`,
  `KAM`, `OSM` are all **0 differences** — the Metal `InnerProd` dot products are **bit-identical**
  to OpenCL (the fixed-order tiled accumulation reproduced exactly; no FMA divergence, so
  `-fno-fast-math` is unnecessary). The only residuals — `DictionaryEulerAngles` (131),
  `ISM` (121), `ISMap` (36), `IndexingSuccessRate` (1) — are the same datasets that differed in
  the Phase 0 OpenCL-vs-OpenCL DI run: pre-existing FZ-sampling non-determinism
  (`mod_so3`/`mod_sampleRFZ`), not the Metal backend. **Phase 2 complete.**
- **Runtime path:** as with MC, the built `DictIndx.metallib` is auto-copied next to the `.cl`
  files (CMake), so no manual copy is needed after a rebuild.

### Phase 3 — remaining active GPU kernels (scoping correction + translation)

- **Scoping finding: the multibeam kernels are dead code.** `MBmoduleOpenCL.cl`
  (`ScatMat`/`CalcLgh`/`CalcLghMaster`) is **not loaded by any module** — the only `MBmodule`
  reference in `Source/` is the CMake comment + the `opencl/SourceList.cmake` copy rule. The
  master-pattern programs (`EBSDmaster`/`ECPmaster`/`TKDmaster`) compute the scattering matrix
  on the **CPU** via `mod_gvectors::CalcLgh_` (LAPACK `ZGEEV`). So the original "Phase 3 =
  multibeam" is a no-op; `MBmoduleOpenCL.cl` is dead (like `mod_DIPCA`) and is intentionally
  **not** built to a metallib.
- **The real remaining active GPU kernels were the MC variants** used by `mod_MCOpenCL`'s
  `foil` and `Ivol` modes (Phase 1 only covered the default MC mode):
  - **`metal/EMMCfoil.metal`** — `MC` (15 args) for foil geometry: tracks electrons through a
    slab of `thickness`, accumulates transmitted electrons (z >= thickness) in the southern
    hemisphere (`LamxSH`/`LamySH`). Uses LFSR113 + Lambert. Port of `EMMCfoil.cl`.
  - **`metal/EMMCxyz.metal`** — `MCxyz` (13 args) for the Ivol (interaction-volume) mode:
    outputs (x,y,z) exit positions (`Lamx`/`Lamy`/`Lamz`), no Lambert. Port of `EMMCxyz.cl`.
  Both share `EMMC.metal`'s validated MC structure (only geometry/output differ) and
  `mod_MCOpenCL` already routes them through the wrappers, so no host changes are needed.
- **CMake:** `EMMCfoil` and `EMMCxyz` added to `EMsoftOO_METAL_KERNELS`.
- **Coverage:** with these, **every live OpenCL kernel** (`EMMC`/`EMMCfoil`/`EMMCxyz`/`InnerProd`)
  is now ported to Metal; the only untranslated `.cl` are dead code (`MBmoduleOpenCL.cl`,
  `DictIndx.cl`'s `ParamEstm`).
- **To verify:** rebuild Metal ON; run `EMMCOpenCL` with `mode='Ivol'` and `mode='foil'`
  namelists, compare `accumSP`/outputs to OpenCL (expect chaotic-MC ±1 jitter like the default
  MC mode, `accumSP` bit-identical).

### Phase 4 — finalization

- **DONE — metallib install rules.** `opencl/SourceList.cmake` now installs the built
  `*.metallib` (from `Bin/opencl/`) to `<install>/opencl/` when Metal is enabled, via
  `INSTALL(DIRECTORY ... FILES_MATCHING PATTERN "*.metallib")` — so `make install`/packaged
  builds resolve them, with no hard-coded kernel list.
- **DONE — default Metal ON for Apple, with OpenCL fallback** (user decision). `Source.cmake`
  defaults `EMsoftOO_ENABLE_Metal_SUPPORT` to ON on `APPLE` (OFF/forced-off elsewhere). At
  configure it probes for the Metal compiler (`xcrun -sdk macosx -f metal`); if absent (e.g. a
  Command-Line-Tools-only install, or Xcode 16+ without the Metal Toolchain component) it
  **automatically falls back to the OpenCL backend** — forces `EMsoftOO_ENABLE_Metal_SUPPORT`
  OFF (re-enabling OpenCL if it had been turned off) and emits a `WARNING` explaining how to
  enable Metal later. So a clean default build succeeds on any Mac: Metal where the toolchain
  exists, OpenCL otherwise. (When Metal is ON, the EMOpenCLLib source-swap selects the Metal
  backend and OpenCL/clfortran are not linked.)
- **DONE (partial) — backend-neutral type/module rename.** Renamed the abstraction type
  `OpenCL_T` → `GPU_T` and the module `mod_CLsupport` → `mod_GPUsupport` (files
  `mod_CLsupport.f90`/`mod_CLsupport_metal.f90` → `mod_GPUsupport.f90`/`mod_GPUsupport_metal.f90`),
  across both backends and all consumers (`mod_MCOpenCL`, `mod_DI`, `mod_EBSDFull`,
  `mod_SEMCLwrappers`, the dead `mod_DIPCA`, `Utilities/EMOpenCLinfo`, the clfortran stub, and
  the CMake source paths). The rename used a word-boundary-safe substitution so `MCOpenCL_T`
  (the MC program class) was left intact. **By user decision the library/folder name
  `EMOpenCLLib` was kept** (renaming it to `EMGPULib` is pure churn across every modality
  CMakeLists + export targets, with no functional benefit). Local handle variables stay named
  `CL`. The OpenCL backend file still uses `clfortran`; only the public abstraction names changed.

## Build-directory hygiene (Metal vs OpenCL)

`clfortran_metal_stub.f90` defines `module clfortran` (constants only) and is compiled **only**
in a Metal build, producing a small `clfortran.mod` in that build tree. If the *same* build
directory is later reconfigured for OpenCL (`-DEMsoftOO_ENABLE_Metal_SUPPORT=OFF`), the OpenCL
backend's `use clfortran` can pick up that stale stub `.mod` (which lacks `CL_SUCCESS`,
`clGetPlatformIDs`, …) instead of the real `clfortran`, giving "Symbol … has no IMPLICIT type"
errors in `mod_GPUsupport.f90`. **Use separate build directories for Metal and OpenCL** (a clean
OpenCL tree never compiles the stub and uses the real `clfortran`). Non-Apple OpenCL builds are
unaffected. A more invasive alternative (not taken) is to drop the stub and put the real
`clfortran` headers on the Metal include path (constants only, no OpenCL link) — that removes
the footgun but couples Metal builds to `CLFortran_INSTALL` and needs Metal re-validation.

## Migration complete

All phases done on `feature/metal-backend`: every live OpenCL kernel ported to Metal, a
backend-neutral `GPU_T`/`mod_GPUsupport` abstraction selected by CMake source-swap, Metal
default-on for Apple with automatic OpenCL fallback, install rules, README + plan docs.
Runtime-verified by the user (MC default/foil/Ivol, DI, EBSDFull, SEMCLwrappers): no issues.
Only `EMOpenCLLib` library name and the dead `MBmoduleOpenCL.cl`/`ParamEstm`/`mod_DIPCA`
remain as-is (deliberately).

## Migration outcome

The OpenCL→Metal migration is functionally complete on `feature/metal-backend`. Every live
OpenCL kernel (`EMMC`/`EMMCfoil`/`EMMCxyz` MC family, `InnerProd` for DI) is ported to MSL and
runs through a drop-in Metal backend selected by a CMake source-swap; Metal is the default GPU
backend on Apple. Validated: MC (default mode) and DI dot products. Remaining before declaring
production-ready: runtime-verify the foil/Ivol MC modes and the `mod_EBSDFull`/`mod_SEMCLwrappers`
Metal paths (all reuse already-validated kernels, so low risk). Deferred polish: the
`GPU_T`/`EMGPULib` rename.

## Post-migration: program renames (2026-05-31)

Now that the GPU compute is backend-neutral, two user-facing programs were renamed away from
the OpenCL-specific names (both old names remain installed as build aliases, so existing
scripts keep working):

- **`EMOpenCLinfo` → `EMGPUinfo`** (`Source/Utilities/EMGPUinfo.f90`). The Metal backend's
  `print_platform_info` was also fleshed out: instead of the previous one-liner it now
  enumerates every Metal device (`MTL::CopyAllDevices` via new `emtl_device_*` shim calls) and
  reports name, unified-memory/low-power/headless/removable attributes, location, recommended
  working-set size, max buffer length, max threadgroup memory, and max threads per
  threadgroup — the Metal analogues of the OpenCL per-device report.
- **`EMMCOpenCL` → `EMMCGPU`** (`Source/SEM/EMMCGPU.f90`). **Important:** the program identity
  embedded in the Monte Carlo HDF5 output is deliberately kept as `'EMMCOpenCL.f90'`
  (`MCfileProgName` in the driver) because downstream readers (`mod_HDFFileInfo`, the
  EBSD/ECP/TKD master programs) match MC files on that exact string and existing files carry
  it. Only the executable/CLI name and the console header changed.

The implementation module/type names (`mod_MCOpenCL`/`MCOpenCL_T`) and the `EMOpenCLLib`
library name are intentionally left unchanged.

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
