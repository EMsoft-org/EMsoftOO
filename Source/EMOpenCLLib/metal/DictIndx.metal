//
// DictIndx.metal
//
// Metal Shading Language translation of the InnerProd kernel from
// opencl/DictIndx.cl — the tiled (block) matrix-product used by dictionary
// indexing to compute experimental x dictionary dot products.  Direct,
// behaviour-preserving port of the BLOCK_SIZE=16 tiled GEMM.
//
// Argument buffer indices match the OpenCL clSetKernelArg indices used by
// InnerProdGPU in mod_DI.f90:
//   0 expt   1 dict   2 Wexp   3 Wdict   4 result
//
// The host dispatches a 2-D grid (Ne, Nd) with a (16,16) threadgroup, so this
// uses dispatchThreadgroups (see emtl_shim emtl_enqueue with a nonzero local
// size).  threads_per_grid gives the OpenCL get_global_size().
//
// Note: only InnerProd is ported; ParamEstm (also in DictIndx.cl) is not used by
// any compiled module, so it is omitted.  If reproducibility vs the OpenCL dot
// products needs tightening, build this kernel with -fno-fast-math (disables FMA
// contraction).  See MetalMigrationPlan.md (Phase 2).
//

#include <metal_stdlib>
using namespace metal;

#define BLOCK_SIZE 16

kernel void InnerProd(device   float* expt   [[buffer(0)]],
                      device   float* dict   [[buffer(1)]],
                      constant int&   Wexp   [[buffer(2)]],
                      constant int&   Wdict  [[buffer(3)]],
                      device   float* result [[buffer(4)]],
                      uint2 gid [[thread_position_in_grid]],
                      uint2 lid [[thread_position_in_threadgroup]],
                      uint2 grp [[threadgroup_position_in_grid]],
                      uint2 tpg [[threads_per_grid]])
{
    // Block index
    int bx = int(grp.x);
    int by = int(grp.y);

    // Thread index inside the block
    int tx = int(lid.x);
    int ty = int(lid.y);

    int aBegin = Wexp * BLOCK_SIZE * by;
    int aEnd   = aBegin + Wexp - 1;
    int aStep  = BLOCK_SIZE;

    int bBegin = BLOCK_SIZE * bx;
    int bstep  = BLOCK_SIZE * Wdict;
    float Csub = 0.0f;

    threadgroup float As[BLOCK_SIZE][BLOCK_SIZE];
    threadgroup float Bs[BLOCK_SIZE][BLOCK_SIZE];

    for (int a = aBegin, b = bBegin; a <= aEnd; a += aStep, b += bstep){

        As[ty][tx] = expt[a + Wexp * ty + tx];
        Bs[ty][tx] = dict[b + Wdict * ty + tx];

        threadgroup_barrier(mem_flags::mem_threadgroup);

        for (int k = 0; k < BLOCK_SIZE; ++k){
            Csub += As[ty][k] * Bs[k][tx];
        }

        threadgroup_barrier(mem_flags::mem_threadgroup);
    }

    result[int(gid.y) * int(tpg.x) + int(gid.x)] = Csub;
}
