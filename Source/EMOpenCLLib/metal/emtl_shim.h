/*
 * emtl_shim.h
 *
 * C-ABI shim over Apple's metal-cpp, so Fortran (mod_MTLsupport) can drive the
 * Metal compute backend through ISO_C_BINDING.  Metal has no Fortran bindings
 * and its API is C++/Obj-C only; this thin C layer is the bridge.
 *
 * Design notes (see MetalMigrationPlan.md, Phase 1):
 *  - All object handles are passed as int64_t (emtl_handle) so the Fortran side
 *    can keep using integer(c_intptr_t) handles exactly as it does for OpenCL.
 *    Internally these are the metal-cpp pointers reinterpret_cast to int64_t.
 *  - Buffers use MTLResourceStorageModeShared (Apple unified memory): write/read
 *    are plain memcpy to/from buffer->contents(); no staging copies.
 *  - The OpenCL host sets kernel args ONCE and then enqueues the kernel many
 *    times in a loop.  Metal binds args per command-encoder, so this shim CACHES
 *    each arg binding (per pipeline) at emtl_set_arg time and (re)applies them
 *    inside emtl_enqueue, which builds a fresh command buffer + encoder each call.
 *    emtl_enqueue commits without waiting; emtl_finish waits on the last commit —
 *    mirroring clEnqueueNDRangeKernel (non-blocking) + clFinish.
 *  - emtl_set_arg distinguishes a buffer argument from raw scalar bytes by
 *    consulting the live-buffer registry: if the value at `ptr` (when its size
 *    matches a handle) is a registered buffer handle it is bound with setBuffer,
 *    otherwise the bytes are bound with setBytes.  This lets the host call a
 *    single uniform set_kernel_arg(index, size, ptr) for both, exactly as OpenCL.
 *
 * Error model: handle-returning calls return 0 on failure; all calls record a
 * message retrievable via emtl_last_error().  mod_MTLsupport mirrors OpenCL_T's
 * error_check using that.
 */
#ifndef EMTL_SHIM_H
#define EMTL_SHIM_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* opaque handle; 0 == null/failure.  intptr_t so it maps 1:1 onto the Fortran
 * integer(c_intptr_t) handles the OpenCL backend already uses (no conversions). */
typedef intptr_t emtl_handle;

/* ---- device / queue lifecycle ---- */
emtl_handle emtl_create_device(void);
emtl_handle emtl_create_queue(emtl_handle device);

/* ---- program (library) / pipeline ----
 * emtl_load_library loads a prebuilt .metallib from disk (built at build time
 * via xcrun metal/metallib).  emtl_get_pipeline builds a compute pipeline state
 * for the named kernel function. */
emtl_handle emtl_load_library(emtl_handle device, const char* metallib_path);
emtl_handle emtl_get_pipeline(emtl_handle device, emtl_handle library, const char* fn_name);

/* ---- buffers (unified memory, StorageModeShared) ----
 * access is accepted for API symmetry with OpenCL flags but ignored. */
emtl_handle emtl_create_buffer(emtl_handle device, size_t nbytes, int access);
void        emtl_write_buffer(emtl_handle buffer, const void* src, size_t nbytes);
void        emtl_read_buffer (emtl_handle buffer, void* dst, size_t nbytes);

/* ---- argument binding (cached per pipeline) ----
 * size==sizeof(emtl_handle) and *ptr a registered buffer  -> setBuffer
 * otherwise                                               -> setBytes */
void emtl_set_arg(emtl_handle pipeline, int index, const void* ptr, size_t size);
/* clear cached bindings for a pipeline (optional; bindings are overwritten by
 * index otherwise) */
void emtl_clear_args(emtl_handle pipeline);

/* ---- dispatch ----
 * Builds a command buffer + compute encoder, applies the pipeline's cached
 * args, dispatches a (gx,gy,gz) grid, commits (no wait).
 * If (lx,ly,lz) are all > 0 they are used as the threadgroup size
 * (dispatchThreadgroups, grid must be a multiple); otherwise a threadgroup
 * size is chosen automatically (dispatchThreads, non-uniform allowed). */
void emtl_enqueue(emtl_handle queue, emtl_handle pipeline,
                  uint64_t gx, uint64_t gy, uint64_t gz,
                  uint64_t lx, uint64_t ly, uint64_t lz);

/* block until the last committed command buffer completes */
void emtl_finish(void);

/* ---- release ---- */
void emtl_release(emtl_handle h);

/* ---- error reporting ----
 * returns 0 if no error since last clear, nonzero otherwise; copies up to
 * buflen-1 chars of the last message into buf (NUL-terminated). */
int  emtl_last_error(char* buf, int buflen);

/* ---- device enumeration / properties (informational; drives EMGPUinfo) ----
 * On macOS every Metal device is enumerated via MTL::CopyAllDevices(); idx is
 * 0-based in [0, emtl_device_count()).  These mirror, as closely as Metal
 * allows, the per-device fields the OpenCL backend reports.  Sizes are in bytes;
 * the Fortran side converts to GB/MB/KB for display.  Out-of-range idx yields 0
 * (or an empty name). */
int      emtl_device_count(void);
/* copies up to buflen-1 chars of the device name into buf (NUL-terminated);
 * returns the number of chars copied. */
int      emtl_device_name(int idx, char* buf, int buflen);
uint64_t emtl_device_recommended_working_set(int idx);  /* recommendedMaxWorkingSetSize */
uint64_t emtl_device_max_buffer_length(int idx);         /* maxBufferLength */
uint64_t emtl_device_max_threadgroup_memory(int idx);    /* maxThreadgroupMemoryLength */
uint64_t emtl_device_current_allocated(int idx);         /* currentAllocatedSize */
uint64_t emtl_device_registry_id(int idx);               /* registryID */
/* maxThreadsPerThreadgroup as a 3-tuple (each component written if non-NULL) */
void     emtl_device_max_threads_per_threadgroup(int idx,
                                                 uint64_t* x, uint64_t* y, uint64_t* z);
/* packed boolean attributes:
 *   bit0 = hasUnifiedMemory, bit1 = lowPower, bit2 = headless, bit3 = removable */
int      emtl_device_flags(int idx);
/* DeviceLocation enum (0 built-in, 1 slot, 2 external, else unspecified) and its
 * location number */
int      emtl_device_location(int idx);
uint64_t emtl_device_location_number(int idx);

#ifdef __cplusplus
}
#endif

#endif /* EMTL_SHIM_H */
