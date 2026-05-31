/*
 * emtl_shim.cpp
 *
 * C-ABI shim over Apple metal-cpp implementing emtl_shim.h.  Compiled as C++17.
 * This is the single translation unit that pulls in the metal-cpp private
 * implementation (the NS_PRIVATE_IMPLEMENTATION / MTL_PRIVATE_IMPLEMENTATION
 * macros must be defined in exactly one .cpp in the whole library).
 *
 * See emtl_shim.h and MetalMigrationPlan.md (Phase 1) for the design.
 */

#define NS_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION

#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>

#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <unordered_map>

#include "emtl_shim.h"

namespace {

// ---- handle <-> pointer helpers ------------------------------------------
inline emtl_handle toHandle(const void* p) { return (emtl_handle)(intptr_t)p; }
template <class T> inline T* fromHandle(emtl_handle h) { return reinterpret_cast<T*>((intptr_t)h); }

// ---- shim state (single-threaded host usage) -----------------------------
std::string g_lastError;

// handles that are MTL::Buffer*; used by emtl_set_arg to pick setBuffer vs setBytes
std::unordered_set<emtl_handle> g_liveBuffers;

struct Arg {
    bool isBuffer = false;
    emtl_handle buf = 0;
    std::vector<uint8_t> bytes;
};
// per-pipeline cached argument bindings, keyed by argument index
std::unordered_map<emtl_handle, std::map<int, Arg>> g_args;

// last committed command buffer (retained until emtl_finish)
MTL::CommandBuffer* g_lastCmd = nullptr;

void setError(const std::string& s) { g_lastError = s; }
void clearError() { g_lastError.clear(); }

} // namespace

extern "C" {

emtl_handle emtl_create_device(void)
{
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    clearError();
    MTL::Device* dev = MTL::CreateSystemDefaultDevice();
    if (dev == nullptr) setError("emtl_create_device: no Metal device available");
    emtl_handle h = toHandle(dev);
    pool->release();
    return h;
}

emtl_handle emtl_create_queue(emtl_handle device)
{
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    clearError();
    MTL::Device* dev = fromHandle<MTL::Device>(device);
    MTL::CommandQueue* q = dev ? dev->newCommandQueue() : nullptr;
    if (q == nullptr) setError("emtl_create_queue: newCommandQueue failed");
    emtl_handle h = toHandle(q);
    pool->release();
    return h;
}

emtl_handle emtl_load_library(emtl_handle device, const char* metallib_path)
{
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    clearError();
    MTL::Device* dev = fromHandle<MTL::Device>(device);
    emtl_handle h = 0;
    if (dev != nullptr) {
        NS::String* path = NS::String::string(metallib_path, NS::UTF8StringEncoding);
        NS::Error* err = nullptr;
        MTL::Library* lib = dev->newLibrary(path, &err);
        if (lib == nullptr) {
            std::string msg = "emtl_load_library: failed to load ";
            msg += metallib_path ? metallib_path : "(null)";
            if (err && err->localizedDescription())
                msg += std::string(": ") + err->localizedDescription()->utf8String();
            setError(msg);
        }
        h = toHandle(lib);
    } else {
        setError("emtl_load_library: null device");
    }
    pool->release();
    return h;
}

emtl_handle emtl_get_pipeline(emtl_handle device, emtl_handle library, const char* fn_name)
{
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    clearError();
    MTL::Device*  dev = fromHandle<MTL::Device>(device);
    MTL::Library* lib = fromHandle<MTL::Library>(library);
    emtl_handle h = 0;
    if (dev && lib) {
        NS::String* name = NS::String::string(fn_name, NS::UTF8StringEncoding);
        MTL::Function* fn = lib->newFunction(name);
        if (fn == nullptr) {
            std::string msg = "emtl_get_pipeline: function not found: ";
            msg += fn_name ? fn_name : "(null)";
            setError(msg);
        } else {
            NS::Error* err = nullptr;
            MTL::ComputePipelineState* pso = dev->newComputePipelineState(fn, &err);
            if (pso == nullptr) {
                std::string msg = "emtl_get_pipeline: newComputePipelineState failed";
                if (err && err->localizedDescription())
                    msg += std::string(": ") + err->localizedDescription()->utf8String();
                setError(msg);
            }
            h = toHandle(pso);
            fn->release();
        }
    } else {
        setError("emtl_get_pipeline: null device or library");
    }
    pool->release();
    return h;
}

emtl_handle emtl_create_buffer(emtl_handle device, size_t nbytes, int /*access*/)
{
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    clearError();
    MTL::Device* dev = fromHandle<MTL::Device>(device);
    MTL::Buffer* buf = dev ? dev->newBuffer((NS::UInteger)nbytes, MTL::ResourceStorageModeShared)
                           : nullptr;
    if (buf == nullptr) setError("emtl_create_buffer: newBuffer failed");
    emtl_handle h = toHandle(buf);
    if (h != 0) g_liveBuffers.insert(h);
    pool->release();
    return h;
}

void emtl_write_buffer(emtl_handle buffer, const void* src, size_t nbytes)
{
    MTL::Buffer* buf = fromHandle<MTL::Buffer>(buffer);
    if (buf && src) std::memcpy(buf->contents(), src, nbytes);
    else setError("emtl_write_buffer: null buffer or source");
}

void emtl_read_buffer(emtl_handle buffer, void* dst, size_t nbytes)
{
    MTL::Buffer* buf = fromHandle<MTL::Buffer>(buffer);
    if (buf && dst) std::memcpy(dst, buf->contents(), nbytes);
    else setError("emtl_read_buffer: null buffer or destination");
}

void emtl_set_arg(emtl_handle pipeline, int index, const void* ptr, size_t size)
{
    Arg a;
    emtl_handle cand = 0;
    if (size == sizeof(emtl_handle) && ptr != nullptr)
        cand = *reinterpret_cast<const emtl_handle*>(ptr);
    if (cand != 0 && g_liveBuffers.count(cand) != 0) {
        a.isBuffer = true;
        a.buf = cand;
    } else {
        a.isBuffer = false;
        a.bytes.assign(reinterpret_cast<const uint8_t*>(ptr),
                       reinterpret_cast<const uint8_t*>(ptr) + size);
    }
    g_args[pipeline][index] = std::move(a);
}

void emtl_clear_args(emtl_handle pipeline)
{
    g_args.erase(pipeline);
}

void emtl_enqueue(emtl_handle queue, emtl_handle pipeline,
                  uint64_t gx, uint64_t gy, uint64_t gz,
                  uint64_t lx, uint64_t ly, uint64_t lz)
{
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    clearError();
    MTL::CommandQueue*         q   = fromHandle<MTL::CommandQueue>(queue);
    MTL::ComputePipelineState* pso = fromHandle<MTL::ComputePipelineState>(pipeline);
    if (!q || !pso) { setError("emtl_enqueue: null queue or pipeline"); pool->release(); return; }

    MTL::CommandBuffer*        cmd = q->commandBuffer();
    MTL::ComputeCommandEncoder* enc = cmd->computeCommandEncoder();
    enc->setComputePipelineState(pso);

    auto it = g_args.find(pipeline);
    if (it != g_args.end()) {
        for (auto& kv : it->second) {
            int idx = kv.first;
            Arg& a = kv.second;
            if (a.isBuffer)
                enc->setBuffer(fromHandle<MTL::Buffer>(a.buf), 0, (NS::UInteger)idx);
            else
                enc->setBytes(a.bytes.data(), (NS::UInteger)a.bytes.size(), (NS::UInteger)idx);
        }
    }

    if (gy == 0) gy = 1;
    if (gz == 0) gz = 1;
    MTL::Size grid = MTL::Size::Make(gx, gy, gz);

    if (lx > 0 && ly > 0 && lz > 0) {
        // explicit threadgroup size (e.g. tiled GEMM): grid must be a multiple
        MTL::Size tg      = MTL::Size::Make(lx, ly, lz);
        MTL::Size ngroups = MTL::Size::Make(gx / lx, gy / ly, gz / lz);
        enc->dispatchThreadgroups(ngroups, tg);
    } else {
        // auto threadgroup; non-uniform grids allowed via dispatchThreads
        NS::UInteger w    = pso->threadExecutionWidth();
        NS::UInteger maxt = pso->maxTotalThreadsPerThreadgroup();
        NS::UInteger tgx  = (w < gx) ? w : (NS::UInteger)gx;
        if (tgx < 1) tgx = 1;
        NS::UInteger tgy  = (tgx > 0) ? (maxt / tgx) : 1;
        if (tgy < 1) tgy = 1;
        if (tgy > gy) tgy = (NS::UInteger)gy;
        MTL::Size tg = MTL::Size::Make(tgx, tgy, 1);
        enc->dispatchThreads(grid, tg);
    }

    enc->endEncoding();
    cmd->commit();

    if (g_lastCmd) g_lastCmd->release();
    g_lastCmd = cmd->retain();   // keep alive until emtl_finish
    pool->release();
}

void emtl_finish(void)
{
    if (g_lastCmd) {
        g_lastCmd->waitUntilCompleted();
        g_lastCmd->release();
        g_lastCmd = nullptr;
    }
}

void emtl_release(emtl_handle h)
{
    if (h == 0) return;
    g_liveBuffers.erase(h);
    g_args.erase(h);  // harmless if h is not a pipeline
    NS::Object* obj = fromHandle<NS::Object>(h);
    if (obj) obj->release();
}

int emtl_last_error(char* buf, int buflen)
{
    if (buf && buflen > 0) {
        int n = (int)g_lastError.size();
        if (n > buflen - 1) n = buflen - 1;
        std::memcpy(buf, g_lastError.c_str(), (size_t)n);
        buf[n] = '\0';
    }
    return g_lastError.empty() ? 0 : 1;
}

// ---- device enumeration / properties (informational; drives EMGPUinfo) ------
//
// MTL::CopyAllDevices() (macOS) returns every Metal device.  We cache the array
// on first use; emtl_device(idx) returns the idx-th MTL::Device* (or nullptr).
namespace {

NS::Array* g_allDevices = nullptr;

NS::Array* allDevices()
{
    if (g_allDevices == nullptr) {
        g_allDevices = MTL::CopyAllDevices();   // retained; never released (process-lifetime)
        if (g_allDevices == nullptr || g_allDevices->count() == 0) {
            // headless / older systems: fall back to the system default device
            MTL::Device* def = MTL::CreateSystemDefaultDevice();
            if (def != nullptr) {
                const NS::Object* objs[1] = { def };
                g_allDevices = NS::Array::array(objs, 1);
                g_allDevices->retain();
            }
        }
    }
    return g_allDevices;
}

MTL::Device* deviceAt(int idx)
{
    NS::Array* arr = allDevices();
    if (arr == nullptr || idx < 0 || (NS::UInteger)idx >= arr->count()) return nullptr;
    return arr->object<MTL::Device>((NS::UInteger)idx);
}

} // namespace

int emtl_device_count(void)
{
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    NS::Array* arr = allDevices();
    int n = arr ? (int)arr->count() : 0;
    pool->release();
    return n;
}

int emtl_device_name(int idx, char* buf, int buflen)
{
    if (buf == nullptr || buflen <= 0) return 0;
    buf[0] = '\0';
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    MTL::Device* dev = deviceAt(idx);
    int n = 0;
    if (dev && dev->name()) {
        const char* nm = dev->name()->utf8String();
        if (nm) {
            n = (int)std::strlen(nm);
            if (n > buflen - 1) n = buflen - 1;
            std::memcpy(buf, nm, (size_t)n);
            buf[n] = '\0';
        }
    }
    pool->release();
    return n;
}

uint64_t emtl_device_recommended_working_set(int idx)
{
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    MTL::Device* dev = deviceAt(idx);
    uint64_t v = dev ? (uint64_t)dev->recommendedMaxWorkingSetSize() : 0;
    pool->release();
    return v;
}

uint64_t emtl_device_max_buffer_length(int idx)
{
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    MTL::Device* dev = deviceAt(idx);
    uint64_t v = dev ? (uint64_t)dev->maxBufferLength() : 0;
    pool->release();
    return v;
}

uint64_t emtl_device_max_threadgroup_memory(int idx)
{
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    MTL::Device* dev = deviceAt(idx);
    uint64_t v = dev ? (uint64_t)dev->maxThreadgroupMemoryLength() : 0;
    pool->release();
    return v;
}

uint64_t emtl_device_current_allocated(int idx)
{
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    MTL::Device* dev = deviceAt(idx);
    uint64_t v = dev ? (uint64_t)dev->currentAllocatedSize() : 0;
    pool->release();
    return v;
}

uint64_t emtl_device_registry_id(int idx)
{
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    MTL::Device* dev = deviceAt(idx);
    uint64_t v = dev ? (uint64_t)dev->registryID() : 0;
    pool->release();
    return v;
}

void emtl_device_max_threads_per_threadgroup(int idx, uint64_t* x, uint64_t* y, uint64_t* z)
{
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    MTL::Device* dev = deviceAt(idx);
    MTL::Size s = dev ? dev->maxThreadsPerThreadgroup() : MTL::Size::Make(0, 0, 0);
    if (x) *x = (uint64_t)s.width;
    if (y) *y = (uint64_t)s.height;
    if (z) *z = (uint64_t)s.depth;
    pool->release();
}

int emtl_device_flags(int idx)
{
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    MTL::Device* dev = deviceAt(idx);
    int f = 0;
    if (dev) {
        if (dev->hasUnifiedMemory()) f |= 1;
        if (dev->lowPower())         f |= 2;
        if (dev->headless())         f |= 4;
        if (dev->removable())        f |= 8;
    }
    pool->release();
    return f;
}

int emtl_device_location(int idx)
{
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    MTL::Device* dev = deviceAt(idx);
    int loc = dev ? (int)dev->location() : -1;
    pool->release();
    return loc;
}

uint64_t emtl_device_location_number(int idx)
{
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();
    MTL::Device* dev = deviceAt(idx);
    uint64_t v = dev ? (uint64_t)dev->locationNumber() : 0;
    pool->release();
    return v;
}

} // extern "C"
