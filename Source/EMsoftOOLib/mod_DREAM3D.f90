! ###################################################################
! Copyright (c) 2013-2023, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_DREAM3D
  !! author: MDG 
  !! version: 1.0 
  !! date: 05/28/21
  !!
  !! module with DREAM3D support functions (mostly reading .dream3d files)

use mod_kinds
use mod_global
use iso_fortran_env, only: int64
use mod_quaternions

IMPLICIT NONE 

type microstructure 

  type(Quaternion3DArray_T)         :: Quaternions
  integer(kind=irg),allocatable     :: FeatureIDs(:,:,:) 
  integer(kind=int64),allocatable   :: dimensions(:)
  real(kind=sgl),allocatable        :: origin(:)
  real(kind=sgl),allocatable        :: gridspacing(:)
  real(kind=sgl)                    :: samplescalefactor
  integer(kind=irg)                 :: numvoxels
end type microstructure 

contains

!--------------------------------------------------------------------------
subroutine ReadDREAM3Dfile(EMsoft, dname, microstr, EApath, FIDpath) 
!DEC$ ATTRIBUTES DLLEXPORT :: ReadDream3Dfile
!! author: MDG 
!! version: 1.0 
!! date: 05/28/21
!!
!! load a microstructure from an HDF5 file 
!! this assumes that the HDF5 interface has been initialized by the calling program 

use mod_EMsoft
use mod_IO 
use HDF5
use mod_HDFsupport 
use mod_rotations

IMPLICIT NONE 

type(EMsoft_T),INTENT(INOUT)          :: EMsoft 
character(fnlen),INTENT(IN)           :: dname 
type(microstructure),INTENT(INOUT)    :: microstr 
character(fnlen),INTENT(IN)           :: EApath(10) 
character(fnlen),INTENT(IN)           :: FIDpath(10) 

type(HDF_T)                           :: HDF 
type(IO_T)                            :: Message 
type(e_T)                             :: eu 
type(q_T)                             :: qq 
type(Quaternion_T)                    :: quat 

character(fnlen)                      :: fname, groupname, dataset  
logical                               :: f_exists, readonly 
integer(kind=irg)                     :: hdferr, iq, ix, iy, iz, d3(3) 
integer(HSIZE_T)                      :: dims(1), dims3(3), dims4(4)
integer(kind=int64),allocatable       :: dimensions(:)
real(kind=sgl),allocatable            :: origin(:)
real(kind=sgl),allocatable            :: gridspacing(:)
integer(kind=irg),allocatable         :: FeatureIDs(:,:,:,:) 
real(kind=sgl),allocatable            :: EulerAngles(:,:,:,:) 

HDF = HDF_T()
fname = trim(EMsoft%generateFilePath('EMdatapathname',dname))

! make sure the file actually exists
inquire(file=trim(fname), exist=f_exists)
if (.not.f_exists) then
  call Message%printError('ReadDREAM3Dfile','DREAM3D file '//trim(fname)//' does not exist')
end if

! open the file in read-only mode
readonly = .TRUE.
hdferr =  HDF%openFile(fname, readonly)

! open the DataContainers group
groupname = 'DataContainers'
  hdferr = HDF%openGroup(groupname)

groupname = EApath(2)
  hdferr = HDF%openGroup(groupname)

!--------------------------------------------------------------------------
! in that group we have the _SIMPL_GEOMETRY group with the array dimensions 
groupname = '_SIMPL_GEOMETRY'
  hdferr = HDF%openGroup(groupname)

dataset = 'DIMENSIONS'
  call HDF%readDatasetInteger64Array(dataset, dims, hdferr, dimensions)
  allocate(microstr%dimensions(dims(1)))
  microstr%dimensions = dimensions
  d3(1:3) = dimensions(1:3)

dataset = 'ORIGIN'
  call HDF%readDatasetFloatArray(dataset, dims, hdferr, origin)
  allocate(microstr%origin(dims(1)))
  microstr%origin = origin

dataset = 'SPACING'
  call HDF%readDatasetFloatArray(dataset, dims, hdferr, gridspacing)
  allocate(microstr%gridspacing(dims(1)))
  microstr%gridspacing = gridspacing

call HDF%pop()
!--------------------------------------------------------------------------

! next we get the EulerAngles and the FeatureIDs arrays from where ever they are located ... 
groupname = EApath(3)
  hdferr = HDF%openGroup(groupname)

dataset = trim(EApath(4))
  call HDF%readDatasetFloatArray(dataset, dims4, hdferr, EulerAngles)

! convert them to a QuaternionArray_T
  microstr%numvoxels = product(dimensions)
  microstr%Quaternions = Quaternion3DArray_T( d3, s = 'd' )
  do iz=1,dimensions(3)
    do iy=1,dimensions(2)
      do ix=1,dimensions(1)
        eu = e_T( edinp = dble(EulerAngles(1:3,ix,iy,iz)) )
        qq = eu%eq()
        quat = Quaternion_T( qd = qq%q_copyd() )
        call microstr%Quaternions%insertQuatin3DArray( (/ ix, iy, iz /), quat)
      end do
    end do
  end do

dataset = trim(FIDpath(4))
  call HDF%readDatasetIntegerArray(dataset, dims4, hdferr, FeatureIDs)
  allocate(microstr%FeatureIDs(dims4(2),dims4(3),dims4(4)))
  microstr%FeatureIDs = FeatureIDs(1,:,:,:)

! clean up some large arrays
  deallocate(EulerAngles, FeatureIDs, dimensions, origin, gridspacing)

! that's it, so we close the file 
call HDF%popall()

call Message%printMessage(' Microstructure data read from DREAM.3D file '//trim(fname))

end subroutine ReadDREAM3Dfile

!--------------------------------------------------------------------------
recursive subroutine getRayOrientations(microstr, rslist, nsp, orlist) 
!DEC$ ATTRIBUTES DLLEXPORT :: getRayOrientations
!! author: MDG 
!! version: 1.0 
!! date: 07/15/21
!!
!! extract orientations from the microstructure along a ray 

type(microstructure),INTENT(INOUT)    :: microstr 
integer(kind=irg),INTENT(IN)          :: nsp
real(kind=dbl),INTENT(IN)             :: rslist(3, nsp)
type(QuaternionArray_T),INTENT(INOUT) :: orlist

integer(kind=irg)                     :: i

! declare a standard QuaternionArray structure
orlist = QuaternionArray_T( n = nsp, s = 'd' )

! get the orientation at each ray position from the 3D quaternion array
do i=1,nsp 
  call orlist%insertQuatinArray( i, microstr%Quaternions%getQuatfrom3DArray( nint(rslist(:,i)) ) )
end do 

! expansion for multi-phase microstructures:  get Feature and Phase IDs as well
end subroutine getRayOrientations

end module mod_DREAM3D
