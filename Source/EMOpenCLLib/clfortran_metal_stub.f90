! ###################################################################
! Copyright (c) 2013-2025, Marc De Graef Research Group/Carnegie Mellon University
! All rights reserved. (BSD-3 license — see other EMsoftOO source headers.)
! ###################################################################

module clfortran
  !! author: Claude Code
  !! version: 1.0
  !! date: 05/29/26
  !!
  !! Minimal stub of the `clfortran` module for the Apple Metal build
  !! (EMsoftOO_ENABLE_Metal_SUPPORT). When the Metal backend is selected the real
  !! clfortran library is not linked, but the GPU program modules still
  !! `use clfortran` to obtain a couple of memory-flag constants that they pass to
  !! GPU_T%create_buffer. The Metal backend ignores the flag value (unified
  !! memory), so only the symbols' existence and integer kind matter here.
  !!
  !! This file is compiled ONLY in a Metal build, in place of the real clfortran
  !! `.mod`. See MetalMigrationPlan.md (Phase 1) and mod_GPUsupport_metal.f90.

use ISO_C_BINDING

IMPLICIT NONE
  public

! cl_mem_flags (cl_bitfield, 64-bit) — standard OpenCL values, unused by Metal
  integer(c_int64_t), parameter :: CL_MEM_READ_WRITE = 1_c_int64_t
  integer(c_int64_t), parameter :: CL_MEM_WRITE_ONLY = 2_c_int64_t
  integer(c_int64_t), parameter :: CL_MEM_READ_ONLY  = 4_c_int64_t

! cl_bool — provided for completeness (only referenced by disabled/commented code)
  integer(c_int32_t), parameter :: CL_FALSE = 0_c_int32_t
  integer(c_int32_t), parameter :: CL_TRUE  = 1_c_int32_t

end module clfortran
