! ###################################################################
! Copyright (c) 2013-2025, Marc De Graef Research Group/Carnegie Mellon University
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are
! permitted provided that the following conditions are met:
!
!     - Redistributions of source code must retain the above copyright notice, this list
!        of conditions and the following disclaimer.
!     - Redistributions in binary form must reproduce the above copyright notice, this
!        list of conditions and the following disclaimer in the documentation and/or
!        other materials provided with the distribution.
!     - Neither the names of Marc De Graef, Carnegie Mellon University nor the names
!        of its contributors may be used to endorse or promote products derived from
!        this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
! OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
! SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
! TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
! BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
! DAMAGE.
! ###################################################################

module mod_GPUsupport
  !! author: MDG / Claude Code
  !! version: 1.0
  !! date: 05/29/26
  !!
  !! Apple Metal backend for the GPU abstraction (Phase 1 of the OpenCL->Metal
  !! migration; see MetalMigrationPlan.md).  This module is a DROP-IN replacement
  !! for the OpenCL mod_GPUsupport.f90: it defines the same module name
  !! (mod_GPUsupport) and the same type (GPU_T) with the identical public
  !! method surface, but every GPU operation is routed through the metal-cpp C
  !! shim (emtl_shim.{h,cpp}).  CMake compiles EITHER this file OR the OpenCL
  !! mod_GPUsupport.f90 depending on EMsoftOO_ENABLE_Metal_SUPPORT, so the program
  !! modules (mod_MCOpenCL, mod_DI, ...) are unchanged.
  !!
  !! Handles (device/queue/library/pipeline/buffer) are integer(c_intptr_t),
  !! exactly as in the OpenCL backend.  Internally these are metal-cpp pointers.
  !! Kernel "programs" are Metal libraries (.metallib, built at build time);
  !! "kernels" are compute pipeline states.  read_source_file derives the
  !! .metallib path from the requested .cl name and build_program loads it.

use mod_EMsoft
use mod_global
use ISO_C_BINDING

IMPLICIT NONE
  private

!--------------------------------------------------------------------------
! C-ABI interfaces to emtl_shim (metal-cpp).  Handles are c_intptr_t.
!--------------------------------------------------------------------------
  interface
    function emtl_create_device() bind(C, name='emtl_create_device') result(h)
      import :: c_intptr_t
      integer(c_intptr_t) :: h
    end function emtl_create_device

    function emtl_create_queue(device) bind(C, name='emtl_create_queue') result(h)
      import :: c_intptr_t
      integer(c_intptr_t), value :: device
      integer(c_intptr_t) :: h
    end function emtl_create_queue

    function emtl_load_library(device, path) bind(C, name='emtl_load_library') result(h)
      import :: c_intptr_t, c_char
      integer(c_intptr_t), value :: device
      character(kind=c_char), dimension(*) :: path
      integer(c_intptr_t) :: h
    end function emtl_load_library

    function emtl_get_pipeline(device, library, fn) bind(C, name='emtl_get_pipeline') result(h)
      import :: c_intptr_t, c_char
      integer(c_intptr_t), value :: device, library
      character(kind=c_char), dimension(*) :: fn
      integer(c_intptr_t) :: h
    end function emtl_get_pipeline

    function emtl_create_buffer(device, nbytes, access) bind(C, name='emtl_create_buffer') result(h)
      import :: c_intptr_t, c_size_t, c_int
      integer(c_intptr_t), value :: device
      integer(c_size_t), value   :: nbytes
      integer(c_int), value      :: access
      integer(c_intptr_t)        :: h
    end function emtl_create_buffer

    subroutine emtl_write_buffer(buffer, src, nbytes) bind(C, name='emtl_write_buffer')
      import :: c_intptr_t, c_ptr, c_size_t
      integer(c_intptr_t), value :: buffer
      type(c_ptr), value         :: src
      integer(c_size_t), value   :: nbytes
    end subroutine emtl_write_buffer

    subroutine emtl_read_buffer(buffer, dst, nbytes) bind(C, name='emtl_read_buffer')
      import :: c_intptr_t, c_ptr, c_size_t
      integer(c_intptr_t), value :: buffer
      type(c_ptr), value         :: dst
      integer(c_size_t), value   :: nbytes
    end subroutine emtl_read_buffer

    subroutine emtl_set_arg(pipeline, indx, ptr, sz) bind(C, name='emtl_set_arg')
      import :: c_intptr_t, c_int, c_ptr, c_size_t
      integer(c_intptr_t), value :: pipeline
      integer(c_int), value      :: indx
      type(c_ptr), value         :: ptr
      integer(c_size_t), value   :: sz
    end subroutine emtl_set_arg

    subroutine emtl_clear_args(pipeline) bind(C, name='emtl_clear_args')
      import :: c_intptr_t
      integer(c_intptr_t), value :: pipeline
    end subroutine emtl_clear_args

    subroutine emtl_enqueue(queue, pipeline, gx, gy, gz, lx, ly, lz) bind(C, name='emtl_enqueue')
      import :: c_intptr_t, c_int64_t
      integer(c_intptr_t), value :: queue, pipeline
      integer(c_int64_t), value  :: gx, gy, gz, lx, ly, lz
    end subroutine emtl_enqueue

    subroutine emtl_finish() bind(C, name='emtl_finish')
    end subroutine emtl_finish

    subroutine emtl_release(h) bind(C, name='emtl_release')
      import :: c_intptr_t
      integer(c_intptr_t), value :: h
    end subroutine emtl_release

    function emtl_last_error(buf, buflen) bind(C, name='emtl_last_error') result(s)
      import :: c_char, c_int
      character(kind=c_char), dimension(*) :: buf
      integer(c_int), value :: buflen
      integer(c_int) :: s
    end function emtl_last_error

    ! ---- device enumeration / properties (informational; drives EMGPUinfo) ----
    function emtl_device_count() bind(C, name='emtl_device_count') result(n)
      import :: c_int
      integer(c_int) :: n
    end function emtl_device_count

    function emtl_device_name(idx, buf, buflen) bind(C, name='emtl_device_name') result(n)
      import :: c_int, c_char
      integer(c_int), value :: idx
      character(kind=c_char), dimension(*) :: buf
      integer(c_int), value :: buflen
      integer(c_int) :: n
    end function emtl_device_name

    function emtl_device_recommended_working_set(idx) &
             bind(C, name='emtl_device_recommended_working_set') result(v)
      import :: c_int, c_int64_t
      integer(c_int), value :: idx
      integer(c_int64_t) :: v
    end function emtl_device_recommended_working_set

    function emtl_device_max_buffer_length(idx) &
             bind(C, name='emtl_device_max_buffer_length') result(v)
      import :: c_int, c_int64_t
      integer(c_int), value :: idx
      integer(c_int64_t) :: v
    end function emtl_device_max_buffer_length

    function emtl_device_max_threadgroup_memory(idx) &
             bind(C, name='emtl_device_max_threadgroup_memory') result(v)
      import :: c_int, c_int64_t
      integer(c_int), value :: idx
      integer(c_int64_t) :: v
    end function emtl_device_max_threadgroup_memory

    function emtl_device_current_allocated(idx) &
             bind(C, name='emtl_device_current_allocated') result(v)
      import :: c_int, c_int64_t
      integer(c_int), value :: idx
      integer(c_int64_t) :: v
    end function emtl_device_current_allocated

    function emtl_device_registry_id(idx) bind(C, name='emtl_device_registry_id') result(v)
      import :: c_int, c_int64_t
      integer(c_int), value :: idx
      integer(c_int64_t) :: v
    end function emtl_device_registry_id

    subroutine emtl_device_max_threads_per_threadgroup(idx, x, y, z) &
               bind(C, name='emtl_device_max_threads_per_threadgroup')
      import :: c_int, c_int64_t
      integer(c_int), value :: idx
      integer(c_int64_t) :: x, y, z
    end subroutine emtl_device_max_threads_per_threadgroup

    function emtl_device_flags(idx) bind(C, name='emtl_device_flags') result(f)
      import :: c_int
      integer(c_int), value :: idx
      integer(c_int) :: f
    end function emtl_device_flags

    function emtl_device_location(idx) bind(C, name='emtl_device_location') result(loc)
      import :: c_int
      integer(c_int), value :: idx
      integer(c_int) :: loc
    end function emtl_device_location

    function emtl_device_location_number(idx) &
             bind(C, name='emtl_device_location_number') result(v)
      import :: c_int, c_int64_t
      integer(c_int), value :: idx
      integer(c_int64_t) :: v
    end function emtl_device_location_number
  end interface

!--------------------------------------------------------------------------
  type, public :: GPU_T
    private
      integer(c_intptr_t) :: device        = 0
      integer(c_intptr_t) :: queue         = 0
      integer(c_intptr_t) :: context       = 0   ! alias of device (for API symmetry)
      integer(kind=irg)   :: numdev        = 0
      integer(kind=irg)   :: seldev        = 1

    contains
      private
        procedure, pass(self) :: error_check_
        procedure, pass(self) :: query_platform_info_
        procedure, pass(self) :: print_platform_info_
        procedure, pass(self) :: read_source_file_
        procedure, pass(self) :: read_source_file_wrapper_
        procedure, pass(self) :: init_PDCCQ_
        procedure, pass(self) :: init_multiPDCCQ_
        procedure, pass(self) :: DI_memory_estimate_
        procedure, pass(self) :: build_program_
        procedure, pass(self) :: get_kernel_
        procedure, pass(self) :: release_program_
        procedure, pass(self) :: create_buffer_
        procedure, pass(self) :: write_buffer_
        procedure, pass(self) :: read_buffer_
        procedure, pass(self) :: set_kernel_arg_
        procedure, pass(self) :: enqueue_kernel_
        procedure, pass(self) :: finish_
        procedure, pass(self) :: release_buffer_
        procedure, pass(self) :: release_kernel_
        procedure, pass(self) :: release_context_queue_

        generic, public :: error_check => error_check_
        generic, public :: query_platform_info => query_platform_info_
        generic, public :: print_platform_info => print_platform_info_
        generic, public :: read_source_file => read_source_file_
        generic, public :: read_source_file_wrapper => read_source_file_wrapper_
        generic, public :: init_PDCCQ => init_PDCCQ_, init_multiPDCCQ_
        generic, public :: DI_memory_estimate => DI_memory_estimate_
        generic, public :: build_program => build_program_
        generic, public :: get_kernel => get_kernel_
        generic, public :: release_program => release_program_
        generic, public :: create_buffer => create_buffer_
        generic, public :: write_buffer => write_buffer_
        generic, public :: read_buffer => read_buffer_
        generic, public :: set_kernel_arg => set_kernel_arg_
        generic, public :: enqueue_kernel => enqueue_kernel_
        generic, public :: finish => finish_
        generic, public :: release_buffer => release_buffer_
        generic, public :: release_kernel => release_kernel_
        generic, public :: release_context_queue => release_context_queue_

  end type GPU_T

  interface GPU_T
    module procedure CL_constructor
  end interface GPU_T

contains

!--------------------------------------------------------------------------
type(GPU_T) function CL_constructor( verb, skipCPU ) result(CL)
!DEC$ ATTRIBUTES DLLEXPORT :: CL_constructor
  !! Metal constructor: create the system default device.

IMPLICIT NONE

logical, INTENT(IN), OPTIONAL :: verb
logical, INTENT(IN), OPTIONAL :: skipCPU

CL%device  = emtl_create_device()
CL%context = CL%device
CL%numdev  = 1
CL%seldev  = 1

end function CL_constructor

!--------------------------------------------------------------------------
subroutine CL_destructor( CL )
!DEC$ ATTRIBUTES DLLEXPORT :: CL_destructor
  !! release the queue and device

IMPLICIT NONE

type(GPU_T), INTENT(INOUT) :: CL

if (CL%queue.ne.0)  call emtl_release(CL%queue)
if (CL%device.ne.0) call emtl_release(CL%device)
CL%queue  = 0
CL%device = 0

end subroutine CL_destructor

!--------------------------------------------------------------------------
! internal helper: surface the last Metal shim error (unless quiet)
!--------------------------------------------------------------------------
recursive subroutine mtl_check_(self, routine, quiet)

use mod_io

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)        :: self
character(*), INTENT(IN)              :: routine
logical, INTENT(IN), OPTIONAL         :: quiet

type(IO_T)                            :: Message
character(len=512, kind=c_char)       :: ebuf
integer(c_int)                        :: status, i
character(len=512)                    :: emsg

if (present(quiet)) then
  if (quiet.eqv..TRUE.) return
end if

ebuf = ''
status = emtl_last_error(ebuf, 512_c_int)
if (status.ne.0) then
  emsg = ''
  do i = 1, 512
    if (ebuf(i:i).eq.C_NULL_CHAR) exit
    emsg(i:i) = ebuf(i:i)
  end do
  call Message%printError(trim(routine), trim(emsg))
end if

end subroutine mtl_check_

!--------------------------------------------------------------------------
recursive subroutine error_check_(self, routine, ierr, nonfatal)
!DEC$ ATTRIBUTES DLLEXPORT :: error_check_
  !! API-compatible error check.  For Metal there is no OpenCL error code, so a
  !! nonzero ierr supplied by the caller is reported; pending shim errors are
  !! also surfaced.

use mod_io

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)          :: self
character(*), INTENT(IN)                :: routine
integer(kind=c_int32_t), INTENT(IN)     :: ierr
logical, INTENT(IN), OPTIONAL           :: nonfatal

type(IO_T)                              :: Message
integer(kind=irg)                       :: iout(1)

if (ierr.ne.0) then
  if (present(nonfatal)) then
    if (nonfatal.eqv..TRUE.) then
      iout(1) = ierr
      call Message%WriteValue('mod_GPUsupport(Metal):'//trim(routine)//' non-fatal error code ', iout, 1)
      return
    end if
  end if
  iout(1) = ierr
  call Message%WriteValue('mod_GPUsupport(Metal):'//trim(routine)//' error code ', iout, 1)
end if
call mtl_check_(self, routine, nonfatal)

end subroutine error_check_

!--------------------------------------------------------------------------
recursive subroutine query_platform_info_(self, p_id, verbose, skCPU)
!DEC$ ATTRIBUTES DLLEXPORT :: query_platform_info_
  !! Metal: single system-default device; minimal informational output.

use mod_io

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)   :: self
integer(kind=irg), INTENT(IN)    :: p_id
logical, INTENT(IN), OPTIONAL    :: verbose
logical, INTENT(IN), OPTIONAL    :: skCPU

type(IO_T)                       :: Message

if (present(verbose)) then
  if (verbose.eqv..TRUE.) call Message%printMessage(' GPU backend: Apple Metal (system default device)')
end if

end subroutine query_platform_info_

!--------------------------------------------------------------------------
recursive subroutine print_platform_info_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: print_platform_info_
  !! author: MDG / Claude Code
  !!
  !! Enumerate every Metal device (MTL::CopyAllDevices via the shim) and report
  !! the per-device properties that most closely correspond to the OpenCL
  !! backend's output (working-set/global memory, max buffer/allocation size,
  !! threadgroup memory, max threads per threadgroup), plus Metal-specific
  !! attributes (unified memory, low-power, headless, removable, location).

use mod_io
use ISO_C_BINDING

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)            :: self
type(IO_T)                                :: Message

integer(c_int)                            :: ndev, i, nc, flags, loc
character(len=256, kind=c_char)           :: cname
character(fnlen)                          :: dname, line, attr, locstr
integer(c_int64_t)                        :: ws, mbl, tgm, calloc, regid, locnum
integer(c_int64_t)                        :: tx, ty, tz
integer(kind=irg)                         :: io_int(4)

ndev = emtl_device_count()

call Message%printMessage(' ')
call Message%printMessage('GPU backend: Apple Metal')
call Message%printMessage('------------------------')
io_int(1) = int(ndev)
call Message%WriteValue('Number of Metal devices: ', io_int, 1, "(I2)")
call Message%printMessage('------------------------')

if (ndev.lt.1) then
  call Message%printMessage( &
     (/ 'No Metal devices were found; this means that EMsoftOO programs with GPU      ', &
        'functionality will not work properly.  Please check your Metal configuration.' /) )
  return
end if

do i = 0, ndev-1
! device name
  cname = ''
  nc = emtl_device_name(i, cname, 256_c_int)
  dname = ''
  call c2f_string_(cname, dname)

  io_int(1) = int(i+1)
  call Message%WriteValue('Device #', io_int, 1, "(I2)")
  call pv_(Message, 'Name:', trim(dname))

! boolean attributes
  flags = emtl_device_flags(i)
  attr = ''
  if (iand(flags, 1).ne.0) attr = trim(attr)//' unified-memory'
  if (iand(flags, 2).ne.0) attr = trim(attr)//' low-power'
  if (iand(flags, 4).ne.0) attr = trim(attr)//' headless'
  if (iand(flags, 8).ne.0) attr = trim(attr)//' removable'
  if (len_trim(attr).eq.0) attr = ' (none)'
  call pv_(Message, 'Attributes:', trim(adjustl(attr)))

! location
  loc    = emtl_device_location(i)
  locnum = emtl_device_location_number(i)
  select case (loc)
    case (0)
      locstr = 'built-in'
    case (1)
      write (locstr,'(A,I0,A)') 'slot (', locnum, ')'
    case (2)
      write (locstr,'(A,I0,A)') 'external (', locnum, ')'
    case default
      locstr = 'unspecified'
  end select
  call pv_(Message, 'Location:', trim(locstr))

! registry id
  regid = emtl_device_registry_id(i)
  write (line,'(I0)') regid
  call pv_(Message, 'Registry ID:', trim(adjustl(line)))

! recommended max working set size (analogous to OpenCL global memory size)
  ws = emtl_device_recommended_working_set(i)
  write (line,'(F12.2,A)') real(ws,dbl)/1024.0_dbl/1024.0_dbl/1024.0_dbl, ' GB'
  call pv_(Message, 'Recommended working set:', trim(adjustl(line)))

! max buffer length (analogous to OpenCL max allocatable memory size)
  mbl = emtl_device_max_buffer_length(i)
  write (line,'(F14.1,A)') real(mbl,dbl)/1024.0_dbl/1024.0_dbl, ' MB'
  call pv_(Message, 'Max buffer length:', trim(adjustl(line)))

! current allocated size
  calloc = emtl_device_current_allocated(i)
  write (line,'(F14.1,A)') real(calloc,dbl)/1024.0_dbl/1024.0_dbl, ' MB'
  call pv_(Message, 'Currently allocated:', trim(adjustl(line)))

! max threadgroup memory (analogous to OpenCL local memory size)
  tgm = emtl_device_max_threadgroup_memory(i)
  write (line,'(I0,A)') tgm/1024_c_int64_t, ' KB'
  call pv_(Message, 'Max threadgroup memory:', trim(adjustl(line)))

! max threads per threadgroup (analogous to OpenCL max work item sizes, 3D)
  call emtl_device_max_threads_per_threadgroup(i, tx, ty, tz)
  write (line,'(I0,A,I0,A,I0)') tx, ' x ', ty, ' x ', tz
  call pv_(Message, 'Max threads/threadgroup:', trim(adjustl(line)))

  call Message%printMessage('------------------------')
end do

call Message%printMessage( &
  (/ '                                                                  ', &
     ' Notes: Apple GPUs use unified memory, so the recommended working ', &
     ' set is the suggested resident budget rather than a hard limit.   ', &
     ' Max buffer length is the largest single allocation; max          ', &
     ' threadgroup memory is the per-threadgroup (local) memory.        ' /) )

end subroutine print_platform_info_

!--------------------------------------------------------------------------
! print a left-justified "  <label> <value>" line with the value column
! aligned across all device properties.
!--------------------------------------------------------------------------
recursive subroutine pv_(Message, label, val)

use mod_io

IMPLICIT NONE

type(IO_T), INTENT(INOUT)    :: Message
character(*), INTENT(IN)     :: label
character(*), INTENT(IN)     :: val

character(fnlen)             :: line
character(len=26)            :: lab

! assigning to a fixed-length variable left-justifies and right-pads the label,
! so the value column lines up (the Aw edit descriptor would right-justify)
lab = label
write (line,'(2X,A,A)') lab, trim(val)
call Message%printMessage(trim(line))

end subroutine pv_

!--------------------------------------------------------------------------
! copy a NUL-terminated c_char buffer into a Fortran character variable
!--------------------------------------------------------------------------
recursive subroutine c2f_string_(cstr, fstr)

IMPLICIT NONE

character(len=*, kind=c_char), INTENT(IN) :: cstr
character(len=*), INTENT(OUT)             :: fstr

integer(kind=irg)                         :: i

fstr = ''
do i = 1, min(len(cstr), len(fstr))
  if (cstr(i:i).eq.C_NULL_CHAR) exit
  fstr(i:i) = cstr(i:i)
end do

end subroutine c2f_string_

!--------------------------------------------------------------------------
recursive subroutine DI_memory_estimate_(self, Nr, Nd, Ne, pl, gpu)
!DEC$ ATTRIBUTES DLLEXPORT :: DI_memory_estimate_
  !! Metal stub (unified memory): no separate device-memory estimate performed.

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)   :: self
integer(kind=8), INTENT(IN)      :: Nr
integer(kind=8), INTENT(IN)      :: Nd
integer(kind=8), INTENT(IN)      :: Ne
integer(kind=4), INTENT(IN)      :: pl
integer(kind=4), INTENT(IN)      :: gpu

! Apple GPUs use unified memory; no separate allocation estimate is needed here.

end subroutine DI_memory_estimate_

!--------------------------------------------------------------------------
! map a requested OpenCL '.cl' source name to the prebuilt '.metallib' path
! (installed alongside the .cl files in the OpenCLpathname folder).
!--------------------------------------------------------------------------
recursive subroutine metallib_path_(EMsoft, sourcefile, fullpath)

use mod_io

IMPLICIT NONE

type(EMsoft_T), intent(INOUT)   :: EMsoft
character(fnlen), INTENT(IN)    :: sourcefile
character(fnlen), INTENT(OUT)   :: fullpath

type(IO_T)                      :: Message
character(fnlen)                :: base, mlib
integer(kind=irg)              :: i, islash, idot
logical                         :: fexist

! strip any directory component -> base name
islash = 0
do i = 1, len_trim(sourcefile)
  if (sourcefile(i:i).eq.'/' .or. sourcefile(i:i).eq.'\') then ! '
    islash = i
  end if
end do
base = trim(sourcefile(islash+1:len_trim(sourcefile)))

! replace trailing '.cl' with '.metallib'
idot = index(base, '.cl', back=.TRUE.)
if (idot.gt.0) then
  mlib = trim(base(1:idot-1))//'.metallib'
else
  mlib = trim(base)//'.metallib'
end if

! resolve under the OpenCLpathname folder (where CMake also installs the metallibs)
fullpath = EMsoft%generateFilePath('OpenCLpathname', mlib)
inquire(file=trim(fullpath), exist=fexist)
if (.not.fexist) then
  call Message%printError('mod_GPUsupport(Metal):metallib_path', &
                          'Metal library not found: '//trim(fullpath))
end if

end subroutine metallib_path_

!--------------------------------------------------------------------------
recursive subroutine read_source_file_(self, EMsoft, sourcefile, csource, slength)
!DEC$ ATTRIBUTES DLLEXPORT :: read_source_file_
  !! Metal: instead of reading OpenCL source, resolve the matching prebuilt
  !! .metallib path and stash it in csource for build_program to load.

IMPLICIT NONE

integer, parameter                                    :: source_length = 50000
class(GPU_T), INTENT(IN)                           :: self
type(EMsoft_T), intent(INOUT)                         :: EMsoft
character(fnlen), INTENT(IN)                          :: sourcefile
character(len=source_length, KIND=c_char), INTENT(OUT):: csource
integer(c_size_t), INTENT(OUT)                        :: slength

character(fnlen)                                      :: fullpath

call metallib_path_(EMsoft, sourcefile, fullpath)
csource = ''
csource(1:len_trim(fullpath)) = trim(fullpath)
slength = len_trim(fullpath)

end subroutine read_source_file_

!--------------------------------------------------------------------------
recursive subroutine read_source_file_wrapper_(self, sourcefile, csource, slength)
!DEC$ ATTRIBUTES DLLEXPORT :: read_source_file_wrapper_
  !! Metal: treat sourcefile as a literal path; map '.cl' -> '.metallib' and
  !! stash it in csource.

IMPLICIT NONE

integer, parameter                                    :: source_length = 50000
class(GPU_T), INTENT(IN)                           :: self
character(fnlen), INTENT(IN)                          :: sourcefile
character(len=source_length, KIND=c_char), INTENT(OUT):: csource
integer(c_size_t), INTENT(OUT)                        :: slength

character(fnlen)                                      :: mlib
integer(kind=irg)                                     :: idot

idot = index(sourcefile, '.cl', back=.TRUE.)
if (idot.gt.0) then
  mlib = trim(sourcefile(1:idot-1))//'.metallib'
else
  mlib = trim(sourcefile)//'.metallib'
end if
csource = ''
csource(1:len_trim(mlib)) = trim(mlib)
slength = len_trim(mlib)

end subroutine read_source_file_wrapper_

!--------------------------------------------------------------------------
recursive subroutine init_PDCCQ_(self, platform, nump, selnump, device, numd, selnumd, devinfo, &
                                  context, command_queue)
!DEC$ ATTRIBUTES DLLEXPORT :: init_PDCCQ_
  !! Metal: create the device (if needed) and a command queue; report a single
  !! device through the OpenCL-style out arguments.

use mod_io

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)           :: self
integer(c_intptr_t), allocatable, target :: platform(:)
integer(kind=irg), INTENT(OUT)           :: nump
integer(kind=irg), INTENT(IN)            :: selnump
integer(c_intptr_t), allocatable, target :: device(:)
integer(kind=irg), INTENT(OUT)           :: numd
integer(kind=irg), INTENT(IN)            :: selnumd
character(fnlen), INTENT(OUT)            :: devinfo
integer(c_intptr_t), target              :: context
integer(c_intptr_t), target              :: command_queue

if (self%device.eq.0) self%device = emtl_create_device()
call mtl_check_(self, 'init_PDCCQ:create_device')
self%queue   = emtl_create_queue(self%device)
call mtl_check_(self, 'init_PDCCQ:create_queue')
self%context = self%device
self%numdev  = 1
self%seldev  = 1

nump = 1
numd = 1
if (allocated(platform)) deallocate(platform)
if (allocated(device))   deallocate(device)
allocate(platform(1), device(1))
platform(1)   = self%device
device(1)     = self%device
context       = self%device
command_queue = self%queue
devinfo       = 'Apple Metal (system default device)'

end subroutine init_PDCCQ_

!--------------------------------------------------------------------------
recursive subroutine init_multiPDCCQ_(self, platform, nump, selnump, device, numd, usenumd, &
                                       selnumd, devinfo, context, command_queue)
!DEC$ ATTRIBUTES DLLEXPORT :: init_multiPDCCQ_
  !! Metal: only the system default device is exposed (single device).

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)            :: self
integer(c_intptr_t), allocatable, target  :: platform(:)
integer(kind=irg), INTENT(OUT)            :: nump
integer(kind=irg), INTENT(IN)             :: selnump
integer(c_intptr_t), allocatable, target  :: device(:)
integer(kind=irg), INTENT(OUT)            :: numd
integer(kind=irg), INTENT(INOUT)          :: usenumd
integer(kind=irg), INTENT(IN)             :: selnumd(usenumd)
character(fnlen), allocatable, INTENT(OUT):: devinfo(:)
integer(c_intptr_t), allocatable, target  :: context(:)
integer(c_intptr_t), allocatable, target  :: command_queue(:)

if (self%device.eq.0) self%device = emtl_create_device()
self%queue   = emtl_create_queue(self%device)
self%context = self%device
self%numdev  = 1
self%seldev  = 1

nump    = 1
numd    = 1
usenumd = 1
if (allocated(platform))      deallocate(platform)
if (allocated(device))        deallocate(device)
if (allocated(context))       deallocate(context)
if (allocated(command_queue)) deallocate(command_queue)
if (allocated(devinfo))       deallocate(devinfo)
allocate(platform(1), device(1), context(1), command_queue(1), devinfo(1))
platform(1)      = self%device
device(1)        = self%device
context(1)       = self%device
command_queue(1) = self%queue
devinfo(1)       = 'Apple Metal (system default device)'

end subroutine init_multiPDCCQ_

!--------------------------------------------------------------------------
recursive function build_program_(self, csource, slength, quiet) result(prog)
!DEC$ ATTRIBUTES DLLEXPORT :: build_program_
  !! Metal: csource holds the .metallib path (set by read_source_file); load it.

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)                   :: self
character(len=*, kind=c_char), target, INTENT(IN):: csource
integer(c_size_t), target, INTENT(IN)            :: slength
logical, INTENT(IN), OPTIONAL                    :: quiet
integer(c_intptr_t)                              :: prog

character(len=:, kind=c_char), allocatable       :: cpath
integer(kind=irg)                                :: n

n = int(slength, irg)
allocate(character(len=n+1, kind=c_char) :: cpath)
cpath = csource(1:n)//C_NULL_CHAR
prog = emtl_load_library(self%device, cpath)
call mtl_check_(self, 'build_program:load_library', quiet)
deallocate(cpath)

end function build_program_

!--------------------------------------------------------------------------
recursive function get_kernel_(self, prog, kernelname, quiet) result(kernel)
!DEC$ ATTRIBUTES DLLEXPORT :: get_kernel_
  !! Metal: build a compute pipeline state for the named function.

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)                       :: self
integer(c_intptr_t), INTENT(IN)                      :: prog
character(len=*), INTENT(IN)                         :: kernelname
logical, INTENT(IN), OPTIONAL                        :: quiet
integer(c_intptr_t)                                  :: kernel

character(len=len_trim(kernelname)+1, kind=c_char)   :: cname

cname = trim(kernelname)//C_NULL_CHAR
kernel = emtl_get_pipeline(self%device, prog, cname)
call mtl_check_(self, 'get_kernel:'//trim(kernelname), quiet)

end function get_kernel_

!--------------------------------------------------------------------------
recursive subroutine release_program_(self, prog, quiet)
!DEC$ ATTRIBUTES DLLEXPORT :: release_program_

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)          :: self
integer(c_intptr_t), INTENT(IN)         :: prog
logical, INTENT(IN), OPTIONAL           :: quiet

if (prog.ne.0) call emtl_release(prog)

end subroutine release_program_

!--------------------------------------------------------------------------
recursive function create_buffer_(self, flags, nbytes, label, quiet) result(buf)
!DEC$ ATTRIBUTES DLLEXPORT :: create_buffer_

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)          :: self
integer(c_int64_t), INTENT(IN)          :: flags
integer(c_size_t), INTENT(IN)           :: nbytes
character(len=*), INTENT(IN)            :: label
logical, INTENT(IN), OPTIONAL           :: quiet
integer(c_intptr_t)                     :: buf

buf = emtl_create_buffer(self%device, nbytes, 0_c_int)
call mtl_check_(self, 'create_buffer:'//trim(label), quiet)

end function create_buffer_

!--------------------------------------------------------------------------
recursive subroutine write_buffer_(self, buf, hostptr, nbytes, label, quiet)
!DEC$ ATTRIBUTES DLLEXPORT :: write_buffer_

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)          :: self
integer(c_intptr_t), INTENT(IN)         :: buf
type(c_ptr), INTENT(IN)                 :: hostptr
integer(c_size_t), INTENT(IN)           :: nbytes
character(len=*), INTENT(IN)            :: label
logical, INTENT(IN), OPTIONAL           :: quiet

call emtl_write_buffer(buf, hostptr, nbytes)
call mtl_check_(self, 'write_buffer:'//trim(label), quiet)

end subroutine write_buffer_

!--------------------------------------------------------------------------
recursive subroutine read_buffer_(self, buf, hostptr, nbytes, label, quiet)
!DEC$ ATTRIBUTES DLLEXPORT :: read_buffer_

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)          :: self
integer(c_intptr_t), INTENT(IN)         :: buf
type(c_ptr), INTENT(IN)                 :: hostptr
integer(c_size_t), INTENT(IN)           :: nbytes
character(len=*), INTENT(IN)            :: label
logical, INTENT(IN), OPTIONAL           :: quiet

call emtl_read_buffer(buf, hostptr, nbytes)
call mtl_check_(self, 'read_buffer:'//trim(label), quiet)

end subroutine read_buffer_

!--------------------------------------------------------------------------
recursive subroutine set_kernel_arg_(self, kernel, argindex, argsize, argptr, label, quiet)
!DEC$ ATTRIBUTES DLLEXPORT :: set_kernel_arg_
  !! Metal: cache the binding for `kernel` (a pipeline) at index argindex.  The
  !! shim decides setBuffer vs setBytes from the live-buffer registry.

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)          :: self
integer(c_intptr_t), INTENT(IN)         :: kernel
integer(c_int32_t), INTENT(IN)          :: argindex
integer(c_size_t), INTENT(IN)           :: argsize
type(c_ptr), INTENT(IN)                 :: argptr
character(len=*), INTENT(IN)            :: label
logical, INTENT(IN), OPTIONAL           :: quiet

call emtl_set_arg(kernel, int(argindex, c_int), argptr, argsize)
call mtl_check_(self, 'set_kernel_arg:'//trim(label), quiet)

end subroutine set_kernel_arg_

!--------------------------------------------------------------------------
recursive subroutine enqueue_kernel_(self, kernel, globalsize, label, localsize, quiet)
!DEC$ ATTRIBUTES DLLEXPORT :: enqueue_kernel_
  !! Metal: dispatch the pipeline.  Work dimension inferred from globalsize.
  !! If localsize is present it is used as the threadgroup size.

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)             :: self
integer(c_intptr_t), INTENT(IN)            :: kernel
integer(c_int64_t), INTENT(IN)             :: globalsize(:)
character(len=*), INTENT(IN)               :: label
integer(c_int64_t), INTENT(IN), OPTIONAL   :: localsize(:)
logical, INTENT(IN), OPTIONAL              :: quiet

integer(c_int64_t)                         :: gx, gy, gz, lx, ly, lz
integer(kind=irg)                          :: nd

nd = size(globalsize)
gx = globalsize(1)
gy = 1_c_int64_t
gz = 1_c_int64_t
if (nd.ge.2) gy = globalsize(2)
if (nd.ge.3) gz = globalsize(3)

lx = 0_c_int64_t
ly = 0_c_int64_t
lz = 0_c_int64_t
if (present(localsize)) then
  nd = size(localsize)
  lx = localsize(1)
  if (nd.ge.2) ly = localsize(2)
  if (nd.ge.3) lz = localsize(3)
  if (ly.eq.0_c_int64_t) ly = 1_c_int64_t
  if (lz.eq.0_c_int64_t) lz = 1_c_int64_t
end if

call emtl_enqueue(self%queue, kernel, gx, gy, gz, lx, ly, lz)
call mtl_check_(self, 'enqueue_kernel:'//trim(label), quiet)

end subroutine enqueue_kernel_

!--------------------------------------------------------------------------
recursive subroutine finish_(self, quiet)
!DEC$ ATTRIBUTES DLLEXPORT :: finish_

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)          :: self
logical, INTENT(IN), OPTIONAL           :: quiet

call emtl_finish()
call mtl_check_(self, 'finish', quiet)

end subroutine finish_

!--------------------------------------------------------------------------
recursive subroutine release_buffer_(self, buf, quiet)
!DEC$ ATTRIBUTES DLLEXPORT :: release_buffer_

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)          :: self
integer(c_intptr_t), INTENT(IN)         :: buf
logical, INTENT(IN), OPTIONAL           :: quiet

if (buf.ne.0) call emtl_release(buf)

end subroutine release_buffer_

!--------------------------------------------------------------------------
recursive subroutine release_kernel_(self, kernel, quiet)
!DEC$ ATTRIBUTES DLLEXPORT :: release_kernel_

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)          :: self
integer(c_intptr_t), INTENT(IN)         :: kernel
logical, INTENT(IN), OPTIONAL           :: quiet

if (kernel.ne.0) call emtl_release(kernel)

end subroutine release_kernel_

!--------------------------------------------------------------------------
recursive subroutine release_context_queue_(self, quiet)
!DEC$ ATTRIBUTES DLLEXPORT :: release_context_queue_
  !! release the command queue (the device is released by the destructor).

IMPLICIT NONE

class(GPU_T), INTENT(INOUT)          :: self
logical, INTENT(IN), OPTIONAL           :: quiet

if (self%queue.ne.0) then
  call emtl_release(self%queue)
  self%queue = 0
end if

end subroutine release_context_queue_

end module mod_GPUsupport
