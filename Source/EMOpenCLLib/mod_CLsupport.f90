! ###################################################################
! Copyright (c) 2013-2020, Marc De Graef Research Group/Carnegie Mellon University
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
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
! USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! ###################################################################
module mod_CLsupport
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/12/20
  !!
  !! OpenCL module; this module is based on the following code, but modified 
  !! substantially and turned into an OpenCL_T class :
  !!--------------------------------------------------------------------------
  !!--------------------------------------------------------------------------
  !! original Copyright information (clfortran's query_platforms_devices.f90)
  !! -----------------------------------------------------------------------------
  !!
  !! Copyright (C) 2013-2014 Company for Advanced Supercomputing Solutions LTD
  !! Bosmat 2a St.
  !! Shoham
  !! Israel 60850
  !! http://www.cass-hpc.com
  !!
  !! Author: Mordechai Butrashvily <support@cass-hpc.com>
  !!
  !! -----------------------------------------------------------------------------
  !!
  !! This program is free software: you can redistribute it and/or modify
  !! it under the terms of the GNU Lesser General Public License as published by
  !! the Free Software Foundation, either version 3 of the License, or
  !! (at your option) any later version.
  !!
  !! This program is distributed in the hope that it will be useful,
  !! but WITHOUT ANY WARRANTY; without even the implied warranty of
  !! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  !! GNU Lesser General Public License for more details.
  !!
  !! You should have received a copy of the GNU Lesser General Public License
  !! along with this program.  If not, see <http://www.gnu.org/licenses/>.
  !!
  !! -----------------------------------------------------------------------------

use mod_EMsoft
use clfortran
use ISO_C_BINDING
use mod_global

IMPLICIT NONE
  private 

  character(45)   ::  errorStrings(68) = (/ &
        'CL_DEVICE_NOT_FOUND                          ', &  ! = -1
        'CL_DEVICE_NOT_AVAILABLE                      ', &  ! = -2
        'CL_COMPILER_NOT_AVAILABLE                    ', &  ! = -3
        'CL_MEM_OBJECT_ALLOCATION_FAILURE             ', &  ! = -4
        'CL_OUT_OF_RESOURCES                          ', &  ! = -5
        'CL_OUT_OF_HOST_MEMORY                        ', &  ! = -6
        'CL_PROFILING_INFO_NOT_AVAILABLE              ', &  ! = -7
        'CL_MEM_COPY_OVERLAP                          ', &  ! = -8
        'CL_IMAGE_FORMAT_MISMATCH                     ', &  ! = -9
        'CL_IMAGE_FORMAT_NOT_SUPPORTED                ', &  ! = -10
        'CL_BUILD_PROGRAM_FAILURE                     ', &  ! = -11
        'CL_MAP_FAILURE                               ', &  ! = -12
        'CL_MISALIGNED_SUB_BUFFER_OFFSET              ', &  ! = -13
        'CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST ', &  ! = -14
        'CL_COMPILE_PROGRAM_FAILURE                   ', &  ! = -15
        'CL_LINKER_NOT_AVAILABLE                      ', &  ! = -16
        'CL_LINK_PROGRAM_FAILURE                      ', &  ! = -17
        'CL_DEVICE_PARTITION_FAILED                   ', &  ! = -18
        'CL_KERNEL_ARG_INFO_NOT_AVAILABLE             ', &  ! = -19
        'UNDEFINED                                    ', &  ! = -20
        'UNDEFINED                                    ', &  ! = -21
        'UNDEFINED                                    ', &  ! = -22
        'UNDEFINED                                    ', &  ! = -23
        'UNDEFINED                                    ', &  ! = -24
        'UNDEFINED                                    ', &  ! = -25
        'UNDEFINED                                    ', &  ! = -26
        'UNDEFINED                                    ', &  ! = -27
        'UNDEFINED                                    ', &  ! = -28
        'UNDEFINED                                    ', &  ! = -29
        'CL_INVALID_VALUE                             ', &  ! = -30
        'CL_INVALID_DEVICE_TYPE                       ', &  ! = -31
        'CL_INVALID_PLATFORM                          ', &  ! = -32
        'CL_INVALID_DEVICE                            ', &  ! = -33
        'CL_INVALID_CONTEXT                           ', &  ! = -34
        'CL_INVALID_QUEUE_PROPERTIES                  ', &  ! = -35
        'CL_INVALID_COMMAND_QUEUE                     ', &  ! = -36
        'CL_INVALID_HOST_PTR                          ', &  ! = -37
        'CL_INVALID_MEM_OBJECT                        ', &  ! = -38
        'CL_INVALID_IMAGE_FORMAT_DESCRIPTOR           ', &  ! = -39
        'CL_INVALID_IMAGE_SIZE                        ', &  ! = -40
        'CL_INVALID_SAMPLER                           ', &  ! = -41
        'CL_INVALID_BINARY                            ', &  ! = -42
        'CL_INVALID_BUILD_OPTION                      ', &  ! = -43
        'CL_INVALID_PROGRAM                           ', &  ! = -44
        'CL_INVALID_PROGRAM_EXECUTABLE                ', &  ! = -45
        'CL_INVALID_KERNEL_NAME                       ', &  ! = -46
        'CL_INVALID_KERNEL_DEFINITION                 ', &  ! = -47
        'CL_INVALID_KERNEL                            ', &  ! = -48
        'CL_INVALID_ARG_INDEX                         ', &  ! = -49
        'CL_INVALID_ARG_VALUE                         ', &  ! = -50
        'CL_INVALID_ARG_SIZE                          ', &  ! = -51
        'UNDEFINED                                    ', &  ! = -52
        'CL_INVALID_WORK_DIMENSION                    ', &  ! = -53
        'CL_INVALID_WORK_GROUP_SIZE                   ', &  ! = -54
        'CL_INVALID_WORK_ITEM_SIZE                    ', &  ! = -55
        'CL_INVALID_GLOBAL_OFFSET                     ', &  ! = -56
        'CL_INVALID_EVENT_WAIT_LIST                   ', &  ! = -57
        'CL_INVALID_EVENT                             ', &  ! = -58
        'CL_INVALID_OPERATION                         ', &  ! = -59
        'CL_INVALID_GL_OBJECT                         ', &  ! = -60
        'CL_INVALID_BUFFER_SIZE                       ', &  ! = -61
        'CL_INVALID_MIP_LEVEL                         ', &  ! = -62
        'CL_INVALID_GLOBAL_WORK_SIZE                  ', &  ! = -63
        'CL_INVALID_PROPERTY                          ', &  ! = -64
        'CL_INVALID_IMAGE_DESCRIPTOR                  ', &  ! = -65
        'CL_INVALID_COMPILER_OPTIONS                  ', &  ! = -66
        'CL_INVALID_LINKER_OPTIONS                    ', &  ! = -67
        'CL_INVALID_DEVICE_PARTITION_COUNT            ' /)  ! = -68


  type,public :: OpenCL_T 
    private 
    ! platform variables 
      character(fnlen), allocatable             :: p_profile(:)
      character(fnlen), allocatable             :: p_version(:)
      character(fnlen), allocatable             :: p_name(:)
      character(fnlen), allocatable             :: p_vendor(:)
      character(fnlen), allocatable             :: p_extensions(:)
      integer(c_intptr_t), allocatable          :: p_ids(:)
! CPU information      
      integer(c_intptr_t), allocatable          :: d_CPUids(:,:)
      integer(c_size_t), allocatable            :: d_CPUmwgs(:,:) 
      integer(c_size_t), allocatable            :: d_CPUmwis(:,:,:) 
      integer(c_size_t), allocatable            :: d_CPUmaxalloc(:,:)
      integer(c_int32_t), allocatable           :: d_CPUcu(:,:)
      integer(c_int64_t), allocatable           :: d_CPUgms(:,:) 
      integer(c_int64_t), allocatable           :: d_CPUmmas(:,:)
      character(fnlen), allocatable             :: d_CPUname(:,:)
! GPU information 
      integer(c_intptr_t), allocatable          :: d_GPUids(:,:)
      integer(c_size_t), allocatable            :: d_GPUmwgs(:,:) 
      integer(c_size_t), allocatable            :: d_GPUmwis(:,:,:) 
      integer(c_size_t), allocatable            :: d_GPUmaxalloc(:,:)
      integer(c_int32_t), allocatable           :: d_GPUcu(:,:)
      integer(c_int64_t), allocatable           :: d_GPUgms(:,:) 
      integer(c_int64_t), allocatable           :: d_GPUmmas(:,:)
      character(fnlen), allocatable             :: d_GPUname(:,:)
! other parameters
      integer(kind=irg)                         :: num_platforms
      integer(kind=irg)                         :: maxCPUdev
      integer(kind=irg)                         :: maxGPUdev
      integer(c_int32_t),allocatable            :: num_CPUdevices(:)
      integer(c_int32_t),allocatable            :: num_GPUdevices(:)
      logical,allocatable                       :: noCPUdevices(:)
      logical,allocatable                       :: noGPUdevices(:) 

    contains
      private
        procedure, pass(self) :: error_check_
        procedure, pass(self) :: query_platform_info_
        procedure, pass(self) :: print_platform_info_
        procedure, pass(self) :: CLread_source_file_
        procedure, pass(self) :: CLread_source_file_wrapper_
        procedure, pass(self) :: CLinit_PDCCQ_
        procedure, pass(self) :: CLinit_multiPDCCQ_
        final :: CL_destructor 

        generic, public :: error_check => error_check_
        generic, public :: query_platform_info => query_platform_info_
        generic, public :: print_platform_info => print_platform_info_
        generic, public :: CLread_source_file => CLread_source_file_
        generic, public :: CLread_source_file_wrapper => CLread_source_file_wrapper_
        generic, public :: CLinit_PDCCQ => CLinit_PDCCQ_, CLinit_multiPDCCQ_

  end type OpenCL_T

!DEC$ ATTRIBUTES DLLEXPORT :: error_check_
!DEC$ ATTRIBUTES DLLEXPORT :: query_platform_info_
!DEC$ ATTRIBUTES DLLEXPORT :: print_platform_info_
!DEC$ ATTRIBUTES DLLEXPORT :: CLread_source_file_
!DEC$ ATTRIBUTES DLLEXPORT :: CLread_source_file_wrapper_
!DEC$ ATTRIBUTES DLLEXPORT :: CLinit_PDCCQ_

  ! the constructor routine for this class 
  interface OpenCL_T
    module procedure CL_constructor
  end interface OpenCL_T

contains

!--------------------------------------------------------------------------
type(OpenCL_T) function CL_constructor( ) result(CL)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/12/20
  !!
  !! constructor for the OpenCL Class
  !!
  !! This constructor queries the OpenCL platforms and devices, and fills in 
  !! all the corresponding arrays  
  
use mod_io 

IMPLICIT NONE

type(IO_T)                      :: Message
integer(c_int32_t)              :: err
integer(c_int)                  :: nplatforms
integer(c_size_t)               :: zero_size = 0
integer(c_size_t)               :: temp_size
integer(kind=irg)               :: i
integer(c_intptr_t)             :: platform_id
integer(c_int32_t)              :: num_devices
integer(c_intptr_t), allocatable, target :: platform_ids(:)


! Get the number of platforms, prior to allocating arrays.
  err = clGetPlatformIDs(0, C_NULL_PTR, nplatforms)
  if (err /= CL_SUCCESS) call Message%printError('clGetPlatformIDs: ','Error quering platforms')
  CL%num_platforms = nplatforms

  if (CL%num_platforms.gt.0) then 
    allocate(CL%p_profile(CL%num_platforms))
    allocate(CL%p_version(CL%num_platforms))
    allocate(CL%p_name(CL%num_platforms))
    allocate(CL%p_vendor(CL%num_platforms))
    allocate(CL%p_extensions(CL%num_platforms))
    allocate(CL%p_ids(CL%num_platforms))
    allocate(CL%num_CPUdevices(CL%num_platforms))
    allocate(CL%num_GPUdevices(CL%num_platforms))
    allocate(CL%noCPUdevices(CL%num_platforms))
    allocate(CL%noGPUdevices(CL%num_platforms))
    CL%noCPUdevices = .FALSE.
    CL%noGPUdevices = .FALSE.
    allocate(platform_ids(CL%num_platforms))

  ! Get platforms IDs.
    err = clGetPlatformIDs(nplatforms, C_LOC(platform_ids), nplatforms)
    if (err /= CL_SUCCESS) call Message%printError('clGetPlatformIDs: ','Error quering platforms')
    CL%p_ids(:) = platform_ids(:)

  ! for each platform, get the number of devices so we can allocate the d_*PUids array 
    CL%maxCPUdev = 0
    CL%maxGPUdev = 0
    do i=1, CL%num_platforms  
      platform_id = platform_ids(i)

  ! device_type = CL_DEVICE_TYPE_CPU
      err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_CPU, 0, C_NULL_PTR, num_devices)
      CL%maxCPUdev = maxval( (/ CL%maxCPUdev, num_devices /) )

  ! device_type = CL_DEVICE_TYPE_GPU
      err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 0, C_NULL_PTR, num_devices)
      CL%maxGPUdev = maxval( (/ CL%maxGPUdev, num_devices /) )
    end do 
    allocate(CL%d_CPUids(CL%num_platforms, CL%maxCPUdev))
    allocate(CL%d_CPUmwgs(CL%num_platforms, CL%maxCPUdev))
    allocate(CL%d_CPUmwis(CL%num_platforms, CL%maxCPUdev,3))
    allocate(CL%d_CPUmaxalloc(CL%num_platforms, CL%maxCPUdev))
    allocate(CL%d_CPUcu(CL%num_platforms, CL%maxCPUdev))
    allocate(CL%d_CPUgms(CL%num_platforms, CL%maxCPUdev))
    allocate(CL%d_CPUmmas(CL%num_platforms, CL%maxCPUdev))
    allocate(CL%d_CPUname(CL%num_platforms, CL%maxCPUdev))
    allocate(CL%d_GPUids(CL%num_platforms, CL%maxGPUdev))
    allocate(CL%d_GPUmwgs(CL%num_platforms, CL%maxGPUdev))
    allocate(CL%d_GPUmwis(CL%num_platforms, CL%maxGPUdev,3))
    allocate(CL%d_GPUmaxalloc(CL%num_platforms, CL%maxGPUdev))
    allocate(CL%d_GPUcu(CL%num_platforms, CL%maxGPUdev))
    allocate(CL%d_GPUgms(CL%num_platforms, CL%maxGPUdev)) 
    allocate(CL%d_GPUmmas(CL%num_platforms, CL%maxGPUdev))
    allocate(CL%d_GPUname(CL%num_platforms, CL%maxGPUdev))

  ! get all relevant information for each platform
    do i=1, CL%num_platforms  
      call CL%query_platform_info_(i)
    end do
  else
  ! the number of platforms is 0 which means that OpenCL is either absent or incorrectly set up
    call Message%printMessage( &
       (/ 'No OpenCL platforms were found; this means that EMsoft programs with OpenCL   ', &
          'functionality will not work properly.  Please check your OpenCL configuration.', &
          '    ----> EMOpenCLinfo: No OpenCL functionality detected on this system       ' /) )
end if

end function CL_constructor

! destructor ... 
!--------------------------------------------------------------------------
subroutine CL_destructor( CL )
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/12/20
  !!
  !! destructor for the OpenCL Class
  !!
IMPLICIT NONE

type(OpenCL_T),INTENT(INOUT)  :: CL 

  deallocate(CL%p_profile)
  deallocate(CL%p_version)
  deallocate(CL%p_name)
  deallocate(CL%p_vendor)
  deallocate(CL%p_extensions)
  deallocate(CL%p_ids)
  deallocate(CL%num_CPUdevices)
  deallocate(CL%num_GPUdevices)
  deallocate(CL%noCPUdevices)
  deallocate(CL%noGPUdevices)
  deallocate(CL%d_CPUids)
  deallocate(CL%d_GPUids)
  deallocate(CL%d_CPUmwgs)
  deallocate(CL%d_CPUmwis)
  deallocate(CL%d_CPUmaxalloc)
  deallocate(CL%d_CPUcu)
  deallocate(CL%d_CPUgms)
  deallocate(CL%d_CPUmmas)
  deallocate(CL%d_CPUname)

end subroutine CL_destructor
! -----------------------------------------------------------------------------
recursive subroutine query_platform_info_(self, p_id)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/12/20
  !!
  !! extract information about OpenCL platforms and devices 
  !! (based on (clfortran's query_platforms_devices.f90)

use ISO_C_BINDING
use mod_global 

IMPLICIT NONE

class(OpenCL_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)  :: p_id 

! Input variable.
integer(c_intptr_t)            :: platform_id

! Helper variables to work with OpenCL API.
integer(c_int32_t)             :: err
integer(c_size_t)              :: zero_size = 0
integer(c_size_t)              :: temp_size
! For quering devices.
integer(c_int64_t)             :: device_type
integer(c_int32_t)             :: num_devices
integer(c_int)                 :: i, j
integer(c_intptr_t), allocatable, target :: device_ids(:)

! String arrays for holding platform details.
character, allocatable, target :: platform_profile(:)
character, allocatable, target :: platform_version(:)
character, allocatable, target :: platform_name(:)
character, allocatable, target :: platform_vendor(:)
character, allocatable, target :: platform_extensions(:)

! String array for holding device name.
character, allocatable, target :: device_name(:)
! Maximum compute units for device.
integer(c_size_t), target      :: device_mwgs, device_mwis(3), device_maxalloc
integer(c_int32_t), target     :: device_cu

integer(c_int64_t), target     :: device_gms, device_mmas

platform_id = self%p_ids(p_id)

! Profile.
err = clGetPlatformInfo(platform_id, CL_PLATFORM_PROFILE, zero_size, C_NULL_PTR, temp_size)
call error_check_(self, 'CLquery_platform_info:clGetPlatformInfo',err)
allocate(platform_profile(temp_size))
err = clGetPlatformInfo(platform_id, CL_PLATFORM_PROFILE, temp_size, C_LOC(platform_profile), temp_size)
call error_check_(self, 'CLquery_platform_info:clGetPlatformInfo',err)
self%p_profile(p_id) = trim(cv_a2s(platform_profile))
deallocate(platform_profile)

! Version.
err = clGetPlatformInfo(platform_id, CL_PLATFORM_VERSION, zero_size, C_NULL_PTR, temp_size)
call error_check_(self, 'CLquery_platform_info:clGetPlatformInfo',err)
allocate(platform_version(temp_size))
err = clGetPlatformInfo(platform_id, CL_PLATFORM_VERSION, temp_size, C_LOC(platform_version), temp_size)
call error_check_(self, 'CLquery_platform_info:clGetPlatformInfo',err)
self%p_version(p_id) = trim(cv_a2s(platform_version))
deallocate(platform_version)

! Name.
err = clGetPlatformInfo(platform_id, CL_PLATFORM_NAME, zero_size, C_NULL_PTR, temp_size)
call error_check_(self, 'CLquery_platform_info:clGetPlatformInfo',err)
allocate(platform_name(temp_size))
err = clGetPlatformInfo(platform_id, CL_PLATFORM_NAME, temp_size, C_LOC(platform_name), temp_size)
call error_check_(self, 'CLquery_platform_info:clGetPlatformInfo',err)
self%p_name(p_id) = trim(cv_a2s(platform_name))
deallocate(platform_name)

! Vendor.
err = clGetPlatformInfo(platform_id, CL_PLATFORM_VENDOR, zero_size, C_NULL_PTR, temp_size)
call error_check_(self, 'CLquery_platform_info:clGetPlatformInfo',err)
allocate(platform_vendor(temp_size))
err = clGetPlatformInfo(platform_id, CL_PLATFORM_VENDOR, temp_size, C_LOC(platform_vendor), temp_size)
call error_check_(self, 'CLquery_platform_info:clGetPlatformInfo',err)
self%p_vendor(p_id) = trim(cv_a2s(platform_vendor))
deallocate(platform_vendor)

! Extensions.
err = clGetPlatformInfo(platform_id, CL_PLATFORM_EXTENSIONS, zero_size, C_NULL_PTR, temp_size)
call error_check_(self, 'CLquery_platform_info:clGetPlatformInfo',err)
allocate(platform_extensions(temp_size))
err = clGetPlatformInfo(platform_id, CL_PLATFORM_EXTENSIONS, temp_size, C_LOC(platform_extensions), temp_size)
call error_check_(self, 'CLquery_platform_info:clGetPlatformInfo',err)
self%p_extensions(p_id) = trim(cv_a2s(platform_extensions))
deallocate(platform_extensions)

!
! Get device information for this platform.
!
! Get CPU device count.

! device_type = CL_DEVICE_TYPE_CPU
err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_CPU, 0, C_NULL_PTR, num_devices)
call error_check_(self, 'CLquery_platform_info:clGetDeviceIDs',err,.TRUE.)

if (err /= CL_SUCCESS .or. num_devices < 1) then
  self%noCPUdevices(p_id) = .TRUE. 
else
  self%num_CPUdevices(p_id) = num_devices 
  allocate(device_ids(num_devices))

! Get device IDs.
  err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_CPU, num_devices, C_LOC(device_ids), num_devices)
  call error_check_(self, 'CLquery_platform_info:clGetDeviceIDs',err)
  self%d_CPUids(p_id,1:num_devices) = device_ids

! Loop over devices and print information.
  do i = 1, num_devices
! Maximum compute units.
    temp_size = 8
    err = clGetDeviceInfo(device_ids(i), CL_DEVICE_MAX_COMPUTE_UNITS, temp_size, C_LOC(device_cu), temp_size)
    call error_check_(self, 'CLquery_platform_info:clGetDeviceInfo',err)
    self%d_CPUcu(p_id, i) = device_cu 

    temp_size = 8
    err = clGetDeviceInfo(device_ids(i), CL_DEVICE_GLOBAL_MEM_SIZE, temp_size, C_LOC(device_gms), temp_size)
    call error_check_(self, 'CLquery_platform_info:clGetDeviceInfo',err)
    device_gms = device_gms/1024/1024/1024
    self%d_CPUgms(p_id, i) = device_gms

! CL_DEVICE_MAX_WORK_GROUP_SIZE
    temp_size = 8 
    err = clGetDeviceInfo(device_ids(i), CL_DEVICE_MAX_WORK_GROUP_SIZE, temp_size, C_LOC(device_mwgs), temp_size)
    call error_check_(self, 'CLquery_platform_info:clGetDeviceInfo',err)
    self%d_CPUmwgs(p_id, i) = device_mwgs

! CL_DEVICE_MAX_WORK_ITEM_SIZES
    temp_size = 8 * 3
    err = clGetDeviceInfo(device_ids(i), CL_DEVICE_MAX_WORK_ITEM_SIZES, temp_size, C_LOC(device_mwis), temp_size)
    call error_check_(self, 'CLquery_platform_info:clGetDeviceInfo',err)
    self%d_CPUmwis(p_id, i, 1:3) = device_mwis

! Name.
    temp_size = 4
    err = clGetDeviceInfo(device_ids(i), CL_DEVICE_NAME, zero_size, C_NULL_PTR, temp_size)
    call error_check_(self, 'CLquery_platform_info:clGetDeviceInfo',err)
    allocate(device_name(temp_size))
    err = clGetDeviceInfo(device_ids(i), CL_DEVICE_NAME, temp_size, C_LOC(device_name), temp_size)
    call error_check_(self, 'CLquery_platform_info:clGetDeviceInfo',err)
    self%d_CPUname(p_id, i) = cv_a2s(device_name)
    deallocate(device_name)
  end do
end if

! Get GPU device count.
! device_type = CL_DEVICE_TYPE_GPU
err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 0, C_NULL_PTR, num_devices)
call error_check_(self, 'CLquery_platform_info:clGetDeviceIDs',err,.TRUE.)

if (err /= CL_SUCCESS .or. num_devices < 1) then
  self%noGPUdevices(p_id) = .TRUE. 
else
  self%num_GPUdevices(p_id) = num_devices 

! Allocate an array to hold device handles.
  if (allocated(device_ids)) deallocate(device_ids)
  allocate(device_ids(num_devices))

! Get device IDs.
  err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, num_devices, C_LOC(device_ids), num_devices)
  call error_check_(self, 'CLquery_platform_info:clGetDeviceIDs',err)
  self%d_GPUids(p_id,1:num_devices) = device_ids


! Loop over devices and print information.
  do i = 1, num_devices
! Maximum compute units.
    temp_size = 8
    err = clGetDeviceInfo(device_ids(i), CL_DEVICE_MAX_COMPUTE_UNITS, temp_size, C_LOC(device_cu), temp_size)
    call error_check_(self, 'CLquery_platform_info:clGetDeviceInfo',err)
    self%d_GPUcu(p_id, i) = device_cu 

    temp_size = 8
    err = clGetDeviceInfo(device_ids(i), CL_DEVICE_GLOBAL_MEM_SIZE, temp_size, C_LOC(device_gms), temp_size)
    call error_check_(self, 'CLquery_platform_info:clGetDeviceInfo',err)
    device_gms = device_gms/1024/1024/1024
    self%d_GPUgms(p_id, i) = device_gms 

! CL_DEVICE_MAX_WORK_GROUP_SIZE
    temp_size = 8 
    err = clGetDeviceInfo(device_ids(i), CL_DEVICE_MAX_WORK_GROUP_SIZE, temp_size, C_LOC(device_mwgs), temp_size)
    call error_check_(self, 'CLquery_platform_info:clGetDeviceInfo',err)
    self%d_GPUmwgs(p_id, i) = device_mwgs 


! CL_DEVICE_MAX_WORK_ITEM_SIZES
    temp_size = 8 * 3
    err = clGetDeviceInfo(device_ids(i), CL_DEVICE_MAX_WORK_ITEM_SIZES, temp_size, C_LOC(device_mwis), temp_size)
    call error_check_(self, 'CLquery_platform_info:clGetDeviceInfo',err)
    self%d_GPUmwis(p_id, i, 1:3) = device_mwis(1:3) 

! Name.
    temp_size = 4
    err = clGetDeviceInfo(device_ids(i), CL_DEVICE_NAME, zero_size, C_NULL_PTR, temp_size)
    call error_check_(self, 'CLquery_platform_info:clGetDeviceInfo',err)
    allocate(device_name(temp_size))
    err = clGetDeviceInfo(device_ids(i), CL_DEVICE_NAME, temp_size, C_LOC(device_name), temp_size)
    call error_check_(self, 'CLquery_platform_info:clGetDeviceInfo',err)
    self%d_GPUname(p_id, i) = cv_a2s(device_name) 
    deallocate(device_name) 

! CL_DEVICE_MAX_MEM_ALLOC_SIZE
    temp_size = 8 
    err = clGetDeviceInfo(device_ids(i), CL_DEVICE_MAX_MEM_ALLOC_SIZE, temp_size, C_LOC(device_maxalloc), temp_size)
    call error_check_(self, 'CLquery_platform_info:clGetDeviceInfo',err)
    device_maxalloc = device_maxalloc/1024/1024
    self%d_GPUmaxalloc(p_id, i) = device_maxalloc

  end do
end if

end subroutine query_platform_info_

! -----------------------------------------------------------------------------
function cv_a2s(inp) result(outp)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/12/20
  !!
  !! (private) auxiliary conversion routine 

IMPLICIT NONE 

character, intent(in),pointer       :: inp(:)

character(fnlen)                    :: outp 
integer(kind=irg)                   :: sz(1), i 

  sz = shape(inp)
  outp = ''

  do i= 1, sz(1)
    outp(i:i) = inp(i)
  end do

end function cv_a2s

! -----------------------------------------------------------------------------
recursive subroutine print_platform_info_(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/12/20
  !!
  !! display information about an OpenCL device (based on (clfortran's query_platforms_devices.f90)

use ISO_C_BINDING
use mod_io
use mod_global 

IMPLICIT NONE

class(OpenCL_T),INTENT(IN)     :: self 

type(IO_T)                     :: Message
integer(kind=irg)              :: io_int(8), i, j

io_int(1) = self%num_platforms
call Message%WriteValue('Number of Platforms: ',io_int,1,"(I2)") 

call Message%printMessage('------------------------')

! Loop over platforms and print information.
do i = 1, self%num_platforms
! Iterate over platforms and get number of devices.
  io_int(1) = i
  call Message%WriteValue('Platform: ', io_int, 1, "(I2/)")

! Profile.
  call Message%printMessage('Profile:      '//self%p_profile(i), frm="(' ',A)")
! Version.
  call Message%printMessage('Version:      '//self%p_version(i), frm="(' ',A)")
! Name.
  call Message%printMessage('Name:         '//self%p_name(i), frm="(' ',A)")
! Vendor.
  call Message%printMessage('Vendor:       '//self%p_vendor(i), frm="(' ',A)")
! Extensions.
  call Message%printMessage('Extensions :  '//self%p_extensions(i), frm="(' ',A/)")
!
! Print CPU device information for this platform.
  if (self%noCPUdevices(i).eqv..TRUE.) then
    call Message%printMessage( 'No CPU devices found on this platform', frm="(/A)")
  else
    io_int(1) = self%num_CPUdevices(i)
    call Message%WriteValue(' # CPU Devices: ', io_int, 1, "(I2)")
    call Message%printMessage('')

! Loop over devices and print information.
    do j = 1, self%num_CPUdevices(i)
      io_int(1:7) = (/ j, self%d_CPUcu(i,j), int(self%d_CPUmwgs(i,j)), int(self%d_CPUmwis(i,j,1)), &
        int(self%d_CPUmwis(i,j,2)),int(self%d_CPUmwis(i,j,3)),int(self%d_CPUgms(i,j))  /)
      call Message%WriteValue('', io_int, 7, &
            "(' Device (#',I2,', CU/MWGS/MWIS/GMS: ',I4,'/',I4,'/',I4,',',I4,',',I4,'/',I3,$)") 
      call Message%printMessage(') - '//trim(self%d_CPUname(i,j)), frm="(A)" )
    end do
  end if

! Print GPU device information for this platform.

  if (self%noGPUdevices(i).eqv..TRUE.) then 
    call Message%printMessage('No GPU devices found on this platform ', frm="(/A)")
  else
    io_int(1) = self%num_GPUdevices(i)
    call Message%WriteValue(' # GPU Devices: ', io_int, 1, "(I2)") 
    call Message%printMessage('')

! Loop over devices and print information.
    do j = 1, self%num_GPUdevices(i)
      io_int = (/ j, self%d_GPUcu(i,j), int(self%d_GPUmwgs(i,j)), int(self%d_GPUmwis(i,j,1)), &
        int(self%d_GPUmwis(i,j,2)),int(self%d_GPUmwis(i,j,3)),int(self%d_GPUgms(i,j)), int(self%d_GPUmaxalloc(i,j))  /)
      call Message%WriteValue('', io_int, 8, &
            "(' Device (#',I2,', CU/MWGS/MWIS/GMS/MAS: ',I4,'/',I4,'/',I4,',',I4,',',I4,'/',I3,',',I4,$)") 
      call Message%printMessage(') - '//trim(self%d_GPUname(i,j)), frm="(A)" )
    end do

    call Message%printMessage('------------------------')
  end if

end do 

call Message%printMessage( &
  (/ '                                            ', &
     ' CU = Compute Units;                        ', &
     ' MWGS = Maximum Work Group Size;            ', &
     ' MWIS = Maximum Work Item Sizes (3D);       ', &
     ' GMS = Global Memory Size (Gb);             ', &
     ' MAS = Maximum Allocatable Memory Size (Mb) '/) )

end subroutine print_platform_info_

!--------------------------------------------------------------------------
recursive subroutine CLread_source_file_(self, EMsoft, sourcefile, csource, slength)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/13/20
  !!
  !! read an OpenCL source file and return the source properly formatted

use mod_EMsoft
use mod_global
use mod_io
use stringconstants
use ISO_C_BINDING

IMPLICIT NONE

integer, parameter                      :: source_length = 50000

class(OpenCL_T),INTENT(IN)              :: self 
type(EMsoft_T),intent(INOUT)            :: EMsoft
character(fnlen), INTENT(IN)            :: sourcefile
character(len=source_length, KIND=c_char),INTENT(OUT) :: csource
integer(c_size_t),INTENT(OUT)           :: slength

type(IO_T)                              :: Message

character(len=source_length),target     :: source
character(fnlen)                        :: fname, clpath, clpath2, tcf
integer(kind=irg)                       :: irec, ierr, ipos, i, j
logical                                 :: fexist
character(1)                            :: EMsoftnativedelimiter
integer(kind=irg)                       :: idx
character(3)                            :: develop

! find the cl file in the main opencl folder or the private folder if the Develop mode equals Yes...
clpath = trim(EMsoft%getConfigParameter('OpenCLpathname') )

! then, determine whether or not the user is working in develop mode by checking for the 
! Develop keyword in the EMsoftconfig.json file... Regular users will only have a single
! opencl folder, but developers have two, so we need to make sure we check both
! locations.  The second location is the private folder...
develop = EMsoft%getConfigParameter('EMdevelop')
clpath2 = ''
if (develop.eq.'Yes') then
  ipos = index(clpath,'Public')
  do i=1,ipos-1
    clpath2(i:i) = clpath(i:i)
  end do
  tcf = 'Private/opencl/'
  do i=ipos,ipos+15
    j = i-ipos+1
    clpath2(i:i) = tcf(j:j)
  end do
end if

if (trim(EMsoft%getConfigParameter('EMsoftplatform')).eq.SC_Windows) then
  EMsoftnativedelimiter = ':'
  idx = 2
else
  EMsoftnativedelimiter = '/'
  idx = 1
end if

if (sourcefile(idx:idx).ne.EMsoftnativedelimiter) then
  fname = trim(clpath)//trim(sourcefile)
else
  fname = trim(sourcefile)
endif
fname = EMsoft%generateFilePath('OpenCLpathname',fname)
inquire(file=trim(fname),exist=fexist)
if (.not.fexist) then 
  if (develop.eq.'Yes') then
   fname = trim(clpath2)//trim(sourcefile)
   fname = EMsoft%generateFilePath('OpenCLpathname',fname)
   inquire(file=trim(fname),exist=fexist)
   if (.not.fexist) then 
     call Message%printError('CLread_source_file','opencl source file '//trim(sourcefile)// &
       ' not found in either opencl folder.'//trim(fname))
   end if
  else
   call Message%printError('CLread_source_file','opencl source  file '//trim(fname)//' not found')
  end if
 end if

! read the source file from the opencl folder
open(unit = dataunit, file = trim(fname), access='direct', status = 'old', &
     action = 'read', iostat = ierr, recl = 1)
if (ierr /= 0) call Message%printError("CLread_source_file: ",'Cannot open file '//fname)

source = ''
irec = 1
do
  read(unit = dataunit, rec = irec, iostat = ierr) source(irec:irec)
  if (ierr /= 0) exit
  if(irec == source_length) call Message%printError("CLread_source_file: ",'Error: CL source file is too big')
  irec = irec + 1
end do
close(unit=dataunit)

csource = trim(source)
csource(irec:irec) = C_NULL_CHAR
slength = irec

end subroutine CLread_source_file_

!--------------------------------------------------------------------------
recursive subroutine CLread_source_file_wrapper_(self, sourcefile, csource, slength)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/13/20
  !!
  !! read an OpenCL source file and return the source properly formatted

use mod_global
use mod_io
use ISO_C_BINDING

IMPLICIT NONE

integer, parameter                      :: source_length = 50000

class(OpenCL_T), INTENT(IN)             :: self
character(fnlen), INTENT(IN)            :: sourcefile
character(len=source_length, KIND=c_char),INTENT(OUT) :: csource
integer(c_size_t),INTENT(OUT)           :: slength

type(IO_T)                              :: Message
character(len=source_length),target     :: source
character(fnlen)                        :: fname, clpath, clpath2, tcf
integer(kind=irg)                       :: irec, ierr, ipos, i, j
logical                                 :: develop, fexist


! find the cl file in the main opencl folder or the private folder if the Develop mode equals Yes...
fname = trim(sourcefile)

inquire(file=trim(fname),exist=fexist)
if (.not.fexist) then 
   call Message%printError('CLread_source_file','opencl source  file '//trim(fname)//' not found')
end if

! read the source file from the opencl folder
open(unit = dataunit, file = trim(fname), access='direct', status = 'old', &
     action = 'read', iostat = ierr, recl = 1)
if (ierr /= 0) call Message%printError("CLread_source_file: ",'Cannot open file '//fname)

source = ''
irec = 1
do
  read(unit = dataunit, rec = irec, iostat = ierr) source(irec:irec)
  if (ierr /= 0) exit
  if(irec == source_length) call Message%printError("CLread_source_file: ",'Error: CL source file is too big')
  irec = irec + 1
end do
close(unit=dataunit)

csource = trim(source)
csource(irec:irec) = C_NULL_CHAR
slength = irec

end subroutine CLread_source_file_wrapper_

!--------------------------------------------------------------------------
recursive subroutine CLinit_PDCCQ_(self, platform, nump, selnump, device, numd, selnumd, devinfo, &
                                  context, command_queue)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/13/20
  !!
  !! initialize a CL platform, device, context, and command queue

use ISO_C_BINDING
use mod_io
use mod_global

IMPLICIT NONE

class(OpenCL_T),INTENT(INOUT)            :: self
integer(c_intptr_t),allocatable, target  :: platform(:)
 !! platform
integer(kind=irg), INTENT(OUT)           :: nump
 !! nump number of platforms available
integer(kind=irg), INTENT(IN)            :: selnump
 !! selnump selected platform number
integer(c_intptr_t),allocatable, target  :: device(:)
 !! device
integer(kind=irg), INTENT(OUT)           :: numd
 !! numd number of devices available in selected platform 
integer(kind=irg), INTENT(IN)            :: selnumd
 !! selnumd selected device number
character(fnlen),INTENT(OUT)             :: devinfo
 !! devinfo
integer(c_intptr_t),target               :: context
 !! context
integer(c_intptr_t),target               :: command_queue
 !! command_queue

type(IO_T)                               :: Message
integer(c_int32_t)                       :: ierr
integer(c_size_t)                        :: cnuminfo
character(fnlen),target                  :: info 
integer(c_intptr_t),target               :: ctx_props(3)
integer(c_int64_t)                       :: cmd_queue_props

! get the platform ID
ierr = clGetPlatformIDs(0, C_NULL_PTR, nump)
call error_check_(self, 'CLinit_PDCCQ:clGetPlatformIDs',ierr)
allocate(platform(nump))
ierr = clGetPlatformIDs(nump, C_LOC(platform), nump)
call error_check_(self, 'CLinit_PDCCQ:clGetPlatformIDs',ierr)

if (selnump.gt.nump) then
  call Message%printError("CLinit_PDCCQ","non-existing platform id requested")
end if

! get the device ID
ierr =  clGetDeviceIDs(platform(selnump), CL_DEVICE_TYPE_GPU, 0, C_NULL_PTR, numd)
call error_check_(self, 'CLinit_PDCCQ:clGetDeviceIDs',ierr)
allocate(device(numd))
ierr =  clGetDeviceIDs(platform(selnump), CL_DEVICE_TYPE_GPU, numd, C_LOC(device), numd)
call error_check_(self, 'CLinit_PDCCQ:clGetDeviceIDs',ierr)

if (selnumd.gt.numd) then
  call Message%printError("CLinit_PDCCQ","non-existing device id requested")
end if

! get the device name and return it as devinfo
ierr = clGetDeviceInfo(device(selnumd), CL_DEVICE_NAME, sizeof(info), C_LOC(info), cnuminfo)
call error_check_(self, 'CLinit_PDCCQ:clGetDeviceInfo',ierr)

if (cnuminfo.gt.fnlen) then 
  call Message%WriteValue("CLinit_PDCCQ","device info string truncated to 132 characters")
  devinfo = trim(info(1:132))
else
  devinfo = trim(info(1:cnuminfo))
end if

! create the context and the command queue
ctx_props(1) = CL_CONTEXT_PLATFORM
ctx_props(2) = platform(selnump)
ctx_props(3) = 0
context = clCreateContext(C_LOC(ctx_props), numd, C_LOC(device),C_NULL_FUNPTR, C_NULL_PTR, ierr)
call error_check_(self, 'CLinit_PDCCQ:clCreateContext',ierr)

cmd_queue_props = CL_QUEUE_PROFILING_ENABLE
command_queue = clCreateCommandQueue(context, device(selnumd), cmd_queue_props, ierr)
call error_check_(self, 'CLinit_PDCCQ:clCreateCommandQueue',ierr)

end subroutine CLinit_PDCCQ_

!--------------------------------------------------------------------------
recursive subroutine CLinit_multiPDCCQ_(self, platform, nump, selnump, device, numd, usenumd, &
                                       selnumd, devinfo, context, command_queue)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/13/20
  !!
  !! initialize a CL platform, multiple devices, a context, and command queues for multi-GPU runs
  !!
  !! this routine still needs to be tested 

use ISO_C_BINDING
use mod_io
use mod_global

IMPLICIT NONE

class(OpenCL_T),INTENT(INOUT)            :: self
integer(c_intptr_t),allocatable, target  :: platform(:)
 !! platform
integer(kind=irg), INTENT(OUT)           :: nump
 !! nump number of platforms available
integer(kind=irg), INTENT(IN)            :: selnump
 !! selnump selected platform number
integer(c_intptr_t),allocatable, target  :: device(:)
 !! device
integer(kind=irg), INTENT(OUT)           :: numd
 !! numd number of devices available in selected platform 
integer(kind=irg), INTENT(INOUT)         :: usenumd
 !! number of devices to be used
integer(kind=irg), INTENT(IN)            :: selnumd(usenumd)
 !! array of selected device numbers
character(fnlen),allocatable,INTENT(OUT) :: devinfo(:)
 !! devinfo
integer(c_intptr_t),allocatable,target   :: context(:)
 !! context
integer(c_intptr_t),allocatable,target   :: command_queue(:)
 !! command_queue

type(IO_T)                               :: Message
integer(c_int32_t)                       :: ierr
integer(c_size_t)                        :: cnuminfo
character(fnlen),target                  :: info 
integer(c_intptr_t),target               :: ctx_props(3)
integer(c_int64_t)                       :: cmd_queue_props
integer(kind=irg)                        :: i

! get the platform ID
ierr = clGetPlatformIDs(0, C_NULL_PTR, nump)
call error_check_(self, 'CLinit_PDCCQ:clGetPlatformIDs',ierr)
allocate(platform(nump))
ierr = clGetPlatformIDs(nump, C_LOC(platform), nump)
call error_check_(self, 'CLinit_PDCCQ:clGetPlatformIDs',ierr)

if (selnump.gt.nump) then
  call Message%printError("CLinit_PDCCQ","non-existing platform id requested")
end if

! get the device ID
ierr =  clGetDeviceIDs(platform(selnump), CL_DEVICE_TYPE_GPU, 0, C_NULL_PTR, numd)
call error_check_(self, 'CLinit_PDCCQ:clGetDeviceIDs',ierr)
allocate(device(numd))
ierr =  clGetDeviceIDs(platform(selnump), CL_DEVICE_TYPE_GPU, numd, C_LOC(device), numd)
call error_check_(self, 'CLinit_PDCCQ:clGetDeviceIDs',ierr)

if(usenumd .gt. numd) then
  call Message%printMessage('')
  call Message%printMessage('Number of devices requested is greater than number of available devices')
  call Message%printMessage('setting number of devices to maximum number of available devices')
  call Message%printMessage('')
  usenumd = numd
end if

do i=1,usenumd
  if (selnumd(i).gt.numd) then
    call Message%printError("CLinit_PDCCQ","non-existing device id requested")
  end if
end do

allocate(context(usenumd), command_queue(usenumd), devinfo(usenumd))

! get the device name and return it as devinfo
do i=1,usenumd
  ierr = clGetDeviceInfo(device(selnumd(i)), CL_DEVICE_NAME, sizeof(info), C_LOC(info), cnuminfo)
  call error_check_(self, 'CLinit_PDCCQ:clGetDeviceInfo',ierr)

  if (cnuminfo.gt.fnlen) then 
    call Message%WriteValue("CLinit_PDCCQ","device info string truncated")
    devinfo(i) = trim(info(1:fnlen))
  else
    devinfo(i) = trim(info(1:cnuminfo))
  end if
end do

! create the context and the command queue
ctx_props(1) = CL_CONTEXT_PLATFORM
ctx_props(2) = platform(selnump)
ctx_props(3) = 0

cmd_queue_props = CL_QUEUE_PROFILING_ENABLE
do i=1,usenumd

  context(i) = clCreateContext(C_LOC(ctx_props), 1, C_LOC(device(selnumd(i))),C_NULL_FUNPTR, C_NULL_PTR, ierr)
  call error_check_(self, 'CLinit_PDCCQ:clCreateContext',ierr)

  command_queue(i) = clCreateCommandQueue(context(i), device(selnumd(i)), cmd_queue_props, ierr)
  call error_check_(self, 'CLinit_PDCCQ:clCreateCommandQueue',ierr)
end do
end subroutine CLinit_multiPDCCQ_

!--------------------------------------------------------------------------
recursive subroutine error_check_(self, routine, ierr, nonfatal)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/07/20
  !!
  !! checks whether or not there was a CL 1.2 error and returns the error message

use ISO_C_BINDING
use mod_io
use mod_global

IMPLICIT NONE

class(OpenCL_T), INTENT(INOUT)          :: self 
character(*),INTENT(IN)                 :: routine
integer(kind=c_int32_t),INTENT(IN)      :: ierr
logical,INTENT(IN),OPTIONAL             :: nonfatal

type(IO_T)                              :: Message
character(fnlen)                        :: estr
integer(kind=irg)                       :: iout(1)

if (ierr.ne.0) then
  if ( (abs(ierr).le.68) .and. (abs(ierr).gt.0) ) then 
    estr = errorStrings(-ierr)
  else
    iout(1) = ierr
    call Message%WriteValue('Unknown CL error code : ', iout, 1)
  end if

  if (present(nonfatal)) then
    if (nonfatal.eqv..TRUE.) then
      call Message%printMessage('error_check', ' Non-fatal error: '//trim(estr) )
    end if
  else
    call Message%printError(trim(routine),trim(estr))
  end if

end if

end subroutine error_check_


end module mod_CLsupport
