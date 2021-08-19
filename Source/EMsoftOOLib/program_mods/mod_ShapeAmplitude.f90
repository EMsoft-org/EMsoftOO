! ###################################################################
! Copyright (c) 2013-2021, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_ShapeAmplitude
  !! author: MDG 
  !! version: 1.0 
  !! date: 08/13/21
  !!
  !! class definition for the EMShapeAmplitude program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMShapeAmplitude program
type, public :: ShapeAmplitudeNameListType
  real(kind=dbl)        :: polyEdgeL
  real(kind=dbl)        :: geometry(10)
  real(kind=dbl)        :: dxyz
  real(kind=dbl)        :: dk
  real(kind=dbl)        :: isovalue
  integer(kind=irg)     :: dims(3)
  integer(kind=irg)     :: nthreads
  character(80)         :: shapetype
  logical               :: shapeIntensity
  character(fnlen)      :: polyhedronFilename 
  character(fnlen)      :: STLFilename 
  character(80)         :: STLheader
  character(fnlen)      :: shampFilename 
end type ShapeAmplitudeNameListType

! class definition
type, public :: ShapeAmplitude_T
private 
  character(fnlen)                  :: nmldeffile = 'EMShapeAmplitude.nml'
  type(ShapeAmplitudeNameListType)  :: nml 
  real(kind=sgl)                    :: volume

contains
private 
  procedure, pass(self) :: get_polyEdgeL_
  procedure, pass(self) :: get_geometry_
  procedure, pass(self) :: get_volume_
  procedure, pass(self) :: get_dxyz_
  procedure, pass(self) :: get_dk_
  procedure, pass(self) :: get_isovalue_
  procedure, pass(self) :: get_dims_
  procedure, pass(self) :: get_nthreads_
  procedure, pass(self) :: get_shapetype_
  procedure, pass(self) :: get_polyhedronFilename_
  procedure, pass(self) :: get_STLFilename_
  procedure, pass(self) :: get_STLheader_
  procedure, pass(self) :: get_shampFilename_
  procedure, pass(self) :: get_shapeIntensity_
  procedure, pass(self) :: set_shapeIntensity_
  procedure, pass(self) :: set_polyEdgeL_
  procedure, pass(self) :: set_volume_
  procedure, pass(self) :: set_geometry_
  procedure, pass(self) :: set_dxyz_
  procedure, pass(self) :: set_dk_
  procedure, pass(self) :: set_isovalue_
  procedure, pass(self) :: set_dims_
  procedure, pass(self) :: set_nthreads_
  procedure, pass(self) :: set_shapetype_
  procedure, pass(self) :: set_polyhedronFilename_
  procedure, pass(self) :: set_STLFilename_
  procedure, pass(self) :: set_STLheader_
  procedure, pass(self) :: set_shampFilename_
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: ShapeAmplitude_
  procedure, pass(self) :: readShapeAmplitude_
  procedure, pass(self) :: ellipsoid_shapefunction_
  procedure, pass(self) :: ellipsoid_shapeamplitude_
  procedure, pass(self) :: prism_shapefunction_
  procedure, pass(self) :: prism_shapeamplitude_
  procedure, pass(self) :: ellipticcylinder_shapefunction_
  procedure, pass(self) :: ellipticcylinder_shapeamplitude_
  procedure, pass(self) :: torus_shapefunction_
  ! procedure, pass(self) :: torus_shapeamplitude_
  ! procedure, pass(self) :: readShapeFunction_

  generic, public :: get_polyEdgeL => get_polyEdgeL_
  generic, public :: get_geometry => get_geometry_
  generic, public :: get_volume => get_volume_
  generic, public :: get_dxyz => get_dxyz_
  generic, public :: get_dk => get_dk_
  generic, public :: get_isovalue => get_isovalue_
  generic, public :: get_dims => get_dims_
  generic, public :: get_nthreads => get_nthreads_
  generic, public :: get_shapetype => get_shapetype_
  generic, public :: get_polyhedronFilename => get_polyhedronFilename_
  generic, public :: get_STLFilename => get_STLFilename_
  generic, public :: get_STLheader => get_STLheader_
  generic, public :: get_shampFilename => get_shampFilename_
  generic, public :: get_shapeIntensity => get_shapeIntensity_
  generic, public :: set_shapeIntensity => set_shapeIntensity_
  generic, public :: set_polyEdgeL => set_polyEdgeL_
  generic, public :: set_volume => set_volume_
  generic, public :: set_geometry => set_geometry_
  generic, public :: set_dxyz => set_dxyz_
  generic, public :: set_dk => set_dk_
  generic, public :: set_isovalue => set_isovalue_
  generic, public :: set_dims => set_dims_
  generic, public :: set_nthreads => set_nthreads_
  generic, public :: set_shapetype => set_shapetype_
  generic, public :: set_polyhedronFilename => set_polyhedronFilename_
  generic, public :: set_STLFilename => set_STLFilename_
  generic, public :: set_STLheader => set_STLheader_
  generic, public :: set_shampFilename => set_shampFilename_
  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: ShapeAmplitude => ShapeAmplitude_
  generic, public :: readShapeAmplitude => readShapeAmplitude_
  ! generic, public :: readShapeFunction => readShapeFunction_

end type ShapeAmplitude_T

! the constructor routine for this class 
interface ShapeAmplitude_T
  module procedure ShapeAmplitude_constructor
end interface ShapeAmplitude_T

contains

!--------------------------------------------------------------------------
type(ShapeAmplitude_T) function ShapeAmplitude_constructor( nmlfile ) result(ShapeAmplitude)
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! constructor for the ShapeAmplitude_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

if (present(nmlfile)) call ShapeAmplitude%readNameList(nmlfile)

end function ShapeAmplitude_constructor

!--------------------------------------------------------------------------
subroutine ShapeAmplitude_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! destructor for the ShapeAmplitude_T Class
 
IMPLICIT NONE

type(ShapeAmplitude_T), INTENT(INOUT)  :: self 

call reportDestructor('ShapeAmplitude_T')

end subroutine ShapeAmplitude_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! read the namelist from an nml file for the ShapeAmplitude_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

real(kind=dbl)        :: polyEdgeL
real(kind=dbl)        :: geometry(10)
real(kind=dbl)        :: dxyz
real(kind=dbl)        :: dk
real(kind=dbl)        :: isovalue
integer(kind=irg)     :: dims(3)
integer(kind=irg)     :: nthreads
logical               :: shapeIntensity
character(80)         :: shapetype
character(fnlen)      :: polyhedronFilename 
character(fnlen)      :: STLFilename 
character(80)         :: STLheader
character(fnlen)      :: shampFilename 

namelist / shampdata / polyEdgeL, geometry, dxyz, dk, isovalue, dims, nthreads, shapetype, polyhedronFilename, &
                       STLFilename, STLheader, shampFilename, shapeIntensity 

! set the input parameters to default values
shapetype = 'polyhedron'
polyhedronFilename = 'undefined'
polyEdgeL = 1.D0
geometry = 0.D0
dxyz = 1.D0
dk = 0.1D0
dims = (/ 128, 128, 128 /)
shapeIntensity = .FALSE.
nthreads = 1
STLFilename = 'undefined'
isovalue = 0.D0
STLheader = '                                                                                '
shampFilename = 'undefined'

! read the name list, depending on the class type
if (.not.skipread) then
! read the namelist file
  open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
  read(UNIT=dataunit,NML=shampdata)
  close(UNIT=dataunit,STATUS='keep')

! check for required entries
  if (trim(shampFilename).eq.'undefined') then
    call Message%printError('EMShapeAmplitude:',' shampFilename is undefined in '//nmlfile)
  end if
end if

self%nml%shapetype = shapetype
self%nml%polyhedronFilename = polyhedronFilename
self%nml%polyEdgeL = polyEdgeL
self%nml%geometry = geometry
self%nml%dxyz = dxyz
self%nml%dk = dk
self%nml%dims = dims
self%nml%shapeIntensity = shapeIntensity
self%nml%nthreads = nthreads
self%nml%STLFilename = STLFilename
self%nml%isovalue = isovalue
self%nml%STLheader = STLheader
self%nml%shampFilename = shampFilename

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! pass the namelist for the ShapeAmplitude_T Class to the calling program

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)          :: self
type(ShapeAmplitudeNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(ShapeAmplitude_T), INTENT(INOUT)  :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 2, n_real = 4
integer(kind=irg)                       :: hdferr,  io_int(n_int), ii
real(kind=dbl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( emnl => self%nml )

! create the group for this namelist
groupname = trim(HDFnames%get_NMLlist())
hdferr = HDF%createGroup(groupname)

! write all the single integers
ii = 0 
if (emnl%shapeIntensity.eqv..TRUE.) ii = 1
io_int = (/ emnl%nthreads, ii /)
intlist(1) = 'nthreads'
intlist(1) = 'shapeIntensity'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the double reals
io_real = (/ emnl%polyEdgeL, emnl%dk, emnl%dxyz, emnl%isovalue /)
reallist(1) = 'polyEdgeL'
reallist(2) = 'dk'
reallist(3) = 'dxyz'
reallist(4) = 'isovalue'
call HDF%writeNMLdbles(io_real, reallist, n_real)

! vectors
dataset = 'dims'
hdferr = HDF%writeDatasetIntegerArray(dataset, emnl%dims, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList_: unable to create dims dataset', hdferr)

dataset = 'geometry'
hdferr = HDF%writeDatasetDoubleArray(dataset, emnl%geometry, 10)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList_: unable to create geometry dataset', hdferr)

! write all the strings
dataset = 'shapetype'
line2(1) = trim(emnl%shapetype)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList_: unable to create shapetype dataset', hdferr)

dataset = 'polyhedronFilename'
line2(1) = trim(emnl%polyhedronFilename)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList_: unable to create polyhedronFilename dataset', hdferr)

dataset = 'STLFilename'
line2(1) = trim(emnl%STLFilename)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList_: unable to create STLFilename dataset', hdferr)

dataset = 'STLheader'
line2(1) = trim(emnl%STLheader)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList_: unable to create STLheader dataset', hdferr)

dataset = 'shampFilename'
line2(1) = trim(emnl%shampFilename)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList_: unable to create shampFilename dataset', hdferr)

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
function get_volume_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_volume_
!! author: MDG 
!! version: 1.0 
!! date: 08/19/21
!!
!! get volume from the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
real(kind=sgl)                             :: out

out = self%volume

end function get_volume_

!--------------------------------------------------------------------------
subroutine set_volume_(self, v) 
!DEC$ ATTRIBUTES DLLEXPORT :: set_volume_
!! author: MDG 
!! version: 1.0 
!! date: 08/19/21
!!
!! set volume from the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
real(kind=sgl),INTENT(IN)                  :: v

self%volume = v

end subroutine set_volume_

!--------------------------------------------------------------------------
function get_polyEdgeL_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_polyEdgeL_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get polyEdgeL from the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
real(kind=dbl)                             :: out

out = self%nml%polyEdgeL

end function get_polyEdgeL_

!--------------------------------------------------------------------------
subroutine set_polyEdgeL_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_polyEdgeL_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set polyEdgeL in the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)                 :: inp

self%nml%polyEdgeL = inp

end subroutine set_polyEdgeL_

!--------------------------------------------------------------------------
function get_geometry_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_geometry_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get geometry from the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
real(kind=dbl)                             :: out(10)

out = self%nml%geometry

end function get_geometry_

!--------------------------------------------------------------------------
subroutine set_geometry_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_geometry_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set geometry in the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)                 :: inp(10)

self%nml%geometry = inp

end subroutine set_geometry_

!--------------------------------------------------------------------------
function get_dxyz_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dxyz_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get dxyz from the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
real(kind=dbl)                             :: out

out = self%nml%dxyz

end function get_dxyz_

!--------------------------------------------------------------------------
subroutine set_dxyz_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_dxyz_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set dxyz in the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)                 :: inp

self%nml%dxyz = inp

end subroutine set_dxyz_

!--------------------------------------------------------------------------
function get_dk_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dk_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get dk from the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
real(kind=dbl)                             :: out

out = self%nml%dk

end function get_dk_

!--------------------------------------------------------------------------
subroutine set_dk_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_dk_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set dk in the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)                 :: inp

self%nml%dk = inp

end subroutine set_dk_

!--------------------------------------------------------------------------
function get_isovalue_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_isovalue_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get isovalue from the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
real(kind=dbl)                             :: out

out = self%nml%isovalue

end function get_isovalue_

!--------------------------------------------------------------------------
subroutine set_isovalue_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_isovalue_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set isovalue in the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)                 :: inp

self%nml%isovalue = inp

end subroutine set_isovalue_

!--------------------------------------------------------------------------
function get_dims_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dims_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get dims from the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
integer(kind=irg)                          :: out(3)

out = self%nml%dims

end function get_dims_

!--------------------------------------------------------------------------
subroutine set_dims_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_dims_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set dims in the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)              :: inp(3)

self%nml%dims = inp

end subroutine set_dims_

!--------------------------------------------------------------------------
function get_nthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nthreads_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get nthreads from the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
integer(kind=irg)                          :: out

out = self%nml%nthreads

end function get_nthreads_

!--------------------------------------------------------------------------
subroutine set_nthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nthreads_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set nthreads in the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)              :: inp

self%nml%nthreads = inp

end subroutine set_nthreads_

!--------------------------------------------------------------------------
function get_shapetype_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_shapetype_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get shapetype from the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
character(10)                              :: out

out = self%nml%shapetype

end function get_shapetype_

!--------------------------------------------------------------------------
subroutine set_shapetype_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_shapetype_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set shapetype in the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
character(10), INTENT(IN)                  :: inp

self%nml%shapetype = inp

end subroutine set_shapetype_

!--------------------------------------------------------------------------
function get_polyhedronFilename_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_polyhedronFilename_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get polyhedronFilename from the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
character(fnlen)                           :: out

out = self%nml%polyhedronFilename

end function get_polyhedronFilename_

!--------------------------------------------------------------------------
subroutine set_polyhedronFilename_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_polyhedronFilename_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set polyhedronFilename in the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)               :: inp

self%nml%polyhedronFilename = inp

end subroutine set_polyhedronFilename_

!--------------------------------------------------------------------------
function get_STLFilename_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_STLFilename_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get STLFilename from the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
character(fnlen)                           :: out

out = self%nml%STLFilename

end function get_STLFilename_

!--------------------------------------------------------------------------
subroutine set_STLFilename_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_STLFilename_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set STLFilename in the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)               :: inp

self%nml%STLFilename = inp

end subroutine set_STLFilename_

!--------------------------------------------------------------------------
function get_STLheader_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_STLheader_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get STLheader from the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
character(80)                              :: out

out = self%nml%STLheader

end function get_STLheader_

!--------------------------------------------------------------------------
subroutine set_STLheader_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_STLheader_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set STLheader in the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
character(80), INTENT(IN)                  :: inp

self%nml%STLheader = inp

end subroutine set_STLheader_

!--------------------------------------------------------------------------
function get_shampFilename_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_shampFilename_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get shampFilename from the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
character(fnlen)                           :: out

out = self%nml%shampFilename

end function get_shampFilename_

!--------------------------------------------------------------------------
subroutine set_shampFilename_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_shampFilename_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set shampFilename in the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)               :: inp

self%nml%shampFilename = inp

end subroutine set_shampFilename_

!--------------------------------------------------------------------------
function get_shapeIntensity_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_shapeIntensity_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get shapeIntensity from the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
logical                                    :: out

out = self%nml%shapeIntensity

end function get_shapeIntensity_

!--------------------------------------------------------------------------
subroutine set_shapeIntensity_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_shapeIntensity_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set shapeIntensity in the ShapeAmplitude_T class

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)     :: self
logical, INTENT(IN)                        :: inp

self%nml%shapeIntensity = inp

end subroutine set_shapeIntensity_

!--------------------------------------------------------------------------
subroutine ellipsoid_shapefunction_(self, shapearray, dims, dxyz, geom)
!DEC$ ATTRIBUTES DLLEXPORT :: ellipsoid_shapefunction_
!! author: MDG 
!! version: 1.0 
!! date: 08/18/21
!!
!! determine the shape function of an ellipsoid

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)    :: self
integer(kind=irg),INTENT(IN)              :: dims(3)
real(kind=sgl),INTENT(INOUT)              :: shapearray( -dims(1):dims(1),-dims(2):dims(2),-dims(3):dims(3) )
real(kind=sgl),INTENT(IN)                 :: dxyz
real(kind=sgl),INTENT(IN)                 :: geom(3) 

integer(kind=irg)                         :: i, j, k 
real(kind=sgl)                            :: x, y, z, a, b, c 

shapearray = 0.0

a = 1.0/geom(1)
b = 1.0/geom(2)
c = 1.0/geom(3)

do i=-dims(1),dims(1)
  x = (i*dxyz*a)**2
  do j=-dims(2),dims(2)
    y = (j*dxyz*b)**2
    do k=-dims(3),dims(3)
      z = (k*dxyz*c)**2
      if (x+y+z.le.1.0) shapearray(i,j,k) = 1.0
    end do 
  end do 
end do

self%volume = 4.0*sngl(cPi)*product(geom)/3.0

write (*,*) maxval(shapearray), sum(shapearray), 4.0*cPi*product(geom)*dxyz**3/3.0


end subroutine ellipsoid_shapefunction_

!--------------------------------------------------------------------------
subroutine ellipticcylinder_shapefunction_(self, shapearray, dims, dxyz, geom)
!DEC$ ATTRIBUTES DLLEXPORT :: ellipticcylinder_shapefunction_
!! author: MDG 
!! version: 1.0 
!! date: 08/18/21
!!
!! determine the shape function of an (elliptic) cylinder 

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)    :: self
integer(kind=irg),INTENT(IN)              :: dims(3)
real(kind=sgl),INTENT(INOUT)              :: shapearray( -dims(1):dims(1),-dims(2):dims(2),-dims(3):dims(3) )
real(kind=sgl),INTENT(IN)                 :: dxyz
real(kind=sgl),INTENT(IN)                 :: geom(3) 

integer(kind=irg)                         :: i, j, k 
real(kind=sgl)                            :: x, y, z, a, b

shapearray = 0.0

a = 1.0/geom(1)
b = 1.0/geom(2)

do i=-dims(1),dims(1)
  x = (i*dxyz*a)**2
  do j=-dims(2),dims(2)
    y = (j*dxyz*b)**2
    do k=-dims(3),dims(3)
      z = abs(k*dxyz)
      if ( (x+y.le.1.0).and.(z.lt.geom(3))) shapearray(i,j,k) = 1.0
    end do 
  end do 
end do

end subroutine ellipticcylinder_shapefunction_

!--------------------------------------------------------------------------
subroutine prism_shapefunction_(self, shapearray, dims, dxyz, geom)
!DEC$ ATTRIBUTES DLLEXPORT :: prism_shapefunction_
!! author: MDG 
!! version: 1.0 
!! date: 08/18/21
!!
!! determine the shape function of a prism

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)    :: self
integer(kind=irg),INTENT(IN)              :: dims(3)
real(kind=sgl),INTENT(INOUT)              :: shapearray( -dims(1):dims(1),-dims(2):dims(2),-dims(3):dims(3) )
real(kind=sgl),INTENT(IN)                 :: dxyz
real(kind=sgl),INTENT(IN)                 :: geom(3) 

integer(kind=irg)                         :: i, j, k 
real(kind=sgl)                            :: x, y, z

shapearray = 0.0

do i=-dims(1),dims(1)
  x = abs(i*dxyz)
  do j=-dims(2),dims(2)
    y = abs(j*dxyz)
    do k=-dims(3),dims(3)
      z = abs(k*dxyz)
      if ( (x.le.geom(1)).and.(y.le.geom(2)).and.(z.le.geom(3)) ) shapearray(i,j,k) = 1.0
    end do 
  end do 
end do

end subroutine prism_shapefunction_

!--------------------------------------------------------------------------
subroutine torus_shapefunction_(self, shapearray, dims, dxyz, geom)
!DEC$ ATTRIBUTES DLLEXPORT :: torus_shapefunction_
!! author: MDG 
!! version: 1.0 
!! date: 08/18/21
!!
!! determine the shape function of a torus  (x^2+y^2+z^2-(a^2+b^2))^2-4a^2(b^2-z^2)

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)    :: self
integer(kind=irg),INTENT(IN)              :: dims(3)
real(kind=sgl),INTENT(INOUT)              :: shapearray( -dims(1):dims(1),-dims(2):dims(2),-dims(3):dims(3) )
real(kind=sgl),INTENT(IN)                 :: dxyz
real(kind=sgl),INTENT(IN)                 :: geom(2) 

integer(kind=irg)                         :: i, j, k 
real(kind=sgl)                            :: x, y, z, a, b, d, ab

shapearray = 0.0

a = geom(1)**2
b = geom(2)**2
ab = a+b 

do i=-dims(1),dims(1)
  x = (i*dxyz)**2
  do j=-dims(2),dims(2)
    y = (j*dxyz)**2
    do k=-dims(3),dims(3)
      z = (k*dxyz)**2
      d = (x+y+z-ab)**2 - 4.0*a*(b-z)
      if ( d.ge.0.0 ) shapearray(i,j,k) = 1.0
    end do 
  end do 
end do

end subroutine torus_shapefunction_

!--------------------------------------------------------------------------
subroutine ellipsoid_shapeamplitude_(self, shamp, dims, dk, geom)
!DEC$ ATTRIBUTES DLLEXPORT :: ellipsoid_shapeamplitude_
!! author: MDG 
!! version: 1.0 
!! date: 08/18/21
!!
!! determine the shape amplitude of an ellipsoid

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)    :: self
integer(kind=irg),INTENT(IN)              :: dims(3)
complex(kind=dbl),INTENT(INOUT)           :: shamp( -dims(1):dims(1)-1,-dims(2):dims(2)-1,-dims(3):dims(3)-1 )
real(kind=dbl),INTENT(IN)                 :: dk
real(kind=sgl),INTENT(IN)                 :: geom(3) 

integer(kind=irg)                         :: i, j, k 
real(kind=dbl)                            :: kvec(3), q
real(kind=dbl)                            :: pre

pre = 4.D0*cPi*dble(product(geom))

do i=-dims(1),dims(1)-1
  do j=-dims(2),dims(2)-1
    do k=-dims(3),dims(3)-1
      kvec = dk * (/ geom(1)*dble(i), geom(2)*dble(j), geom(3)*dble(k) /)  ! this is the Fourier space frequency vector k
      q = NORM2(kvec)
      if (q.eq.0.0) then 
        shamp(0,0,0) = cmplx(pre/3.D0)
      else
        shamp(i,j,k) = cmplx(pre * (sin(q)/q**3 - cos(q)/q**2)) 
      end if 
    end do 
  end do 
end do

end subroutine ellipsoid_shapeamplitude_

!--------------------------------------------------------------------------
function sinc( x ) result(s)
!! author: MDG 
!! version: 1.0 
!! date: 08/19/21
!!
!! sinc(x) = sin(x)/x

IMPLICIT NONE

real(kind=dbl),INTENT(IN)       :: x
real(kind=dbl)                  :: s 

if (x.eq.0.D0) then 
  s = 1.D0 
else 
  s = sin(x)/x 
end if 

end function sinc

!--------------------------------------------------------------------------
subroutine prism_shapeamplitude_(self, shamp, dims, dk, geom)
!DEC$ ATTRIBUTES DLLEXPORT :: prism_shapeamplitude_
!! author: MDG 
!! version: 1.0 
!! date: 08/18/21
!!
!! determine the shape amplitude of a rectangular prism

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)    :: self
integer(kind=irg),INTENT(IN)              :: dims(3)
complex(kind=dbl),INTENT(INOUT)           :: shamp( -dims(1):dims(1)-1,-dims(2):dims(2)-1,-dims(3):dims(3)-1 )
real(kind=dbl),INTENT(IN)                 :: dk
real(kind=sgl),INTENT(IN)                 :: geom(3) 

integer(kind=irg)                         :: i, j, k 
real(kind=dbl)                            :: kvec(3), q
real(kind=dbl)                            :: pre

pre = 8.D0*dble(product(geom))

do i=-dims(1),dims(1)-1
  do j=-dims(2),dims(2)-1
    do k=-dims(3),dims(3)-1
      kvec = dk * (/ geom(1)*dble(i), geom(2)*dble(j), geom(3)*dble(k) /)  ! this is the Fourier space frequency vector k
      q = NORM2(kvec)
      if (q.eq.0.0) then 
        shamp(0,0,0) = cmplx(pre)
      else
        shamp(i,j,k) = cmplx(pre * sinc(kvec(1)) * sinc(kvec(2)) * sinc(kvec(3)) ) 
      end if 
    end do 
  end do 
end do

end subroutine prism_shapeamplitude_

!--------------------------------------------------------------------------
subroutine ellipticcylinder_shapeamplitude_(self, shamp, dims, dk, geom)
!DEC$ ATTRIBUTES DLLEXPORT :: ellipticcylinder_shapeamplitude_
!! author: MDG 
!! version: 1.0 
!! date: 08/18/21
!!
!! determine the shape amplitude of a rectangular prism

use mod_math

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)    :: self
integer(kind=irg),INTENT(IN)              :: dims(3)
complex(kind=dbl),INTENT(INOUT)           :: shamp( -dims(1):dims(1)-1,-dims(2):dims(2)-1,-dims(3):dims(3)-1 )
real(kind=dbl),INTENT(IN)                 :: dk
real(kind=sgl),INTENT(IN)                 :: geom(3) 

integer(kind=irg)                         :: i, j, k 
real(kind=dbl)                            :: qq(2), q
real(kind=dbl)                            :: pre, beta

beta = geom(2)/geom(1) 
pre = 4.D0* cPi * geom(1) * geom(3)

do i=-dims(1),dims(1)-1
  do j=-dims(2),dims(2)-1
    do k=-dims(3),dims(3)-1
      qq = dk * (/ dble(i), beta*dble(j)/)  ! this is the Fourier space frequency vector k
      q = NORM2(qq)
      if ((q.eq.0.0).and.(k.eq.0)) then 
        shamp(0,0,0) = cmplx( 2.D0*cPi*product(dble(geom)) )
      else
        shamp(i,j,k) = cmplx(pre * bessj(1, geom(1)*q) * sinc(geom(3)*dble(k)) / q ) 
      end if 
    end do 
  end do 
end do

end subroutine ellipticcylinder_shapeamplitude_

!--------------------------------------------------------------------------
subroutine ShapeAmplitude_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: ShapeAmplitude_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! perform the computations

use mod_EMsoft
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_timing
use mod_PGA3D 
use mod_PGA3Dsupport
use mod_polyhedra
use stringconstants
use mod_STL
use mod_MCA 
use mod_IO

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)  :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

type(STL_T)                             :: STL 
type(MCA_T)                             :: MCA 
type(polyhedron_T)                      :: shape     
type(PGA3D_T)                           :: mv
type(HDF_T)                             :: HDF 
type(IO_T)                              :: Message
type(Timing_T)                          :: timer

character(fnlen,kind=c_char)            :: line2(1)
character(fnlen)                        :: polyname, stlname, sname, datafile, dataset, datagroupname, attributename, &
                                           HDF_FileVersion 
character(80)                           :: header
character(11)                           :: dstr
character(15)                           :: tstrb
character(15)                           :: tstre
real(kind=dbl)                          :: edgeL, dk 
real(kind=sgl)                          :: iso, mi, ma, dxyz 
real(kind=sgl),allocatable              :: sf(:,:,:)
complex(kind=dbl),allocatable           :: shamp(:,:,:)
real(kind=sgl),allocatable              :: shampreal(:,:,:), shint(:,:,:)
integer(kind=irg)                       :: i, j, k, nthr, dims(3), ntriangles, hdferr
logical                                 :: analytical

associate(emnl => self%nml)
dk = emnl%dk
nthr = emnl%nthreads
dims = emnl%dims
dxyz = sngl(emnl%dxyz)

if (trim(emnl%shapetype).eq.'polyhedron') then 
  call Message%printMessage(' Initializing 3D Projective Geometric Algebra module')
  call PGA3D_initialize()
  edgeL = emnl%polyEdgeL 
  ! check for / in the filename !!!
  if (index(trim(emnl%polyhedronFilename),'/').ne.0) then 
    polyname = ''
    shape = polyhedron_T( polyname, edgeL, shapefile=emnl%polyhedronFilename )
  else 
    polyname = trim(emnl%polyhedronFilename)
    shape = polyhedron_T( polyname, edgeL )
  end if 
  call shape%polyhedron_info()
  call self%set_volume_(sngl(shape%get_volume()))

! first generate the shape function 
  call Message%printMessage(' Generating the Shape Function')
  allocate(sf(-dims(1):dims(1),-dims(2):dims(2),-dims(3):dims(3)))
  call shape%polyhedron_shapefunction(sf, dims, emnl%dxyz)

! then the shape amplitude using the Komrska approach
  call Message%printMessage(' Generating the Shape Amplitude')
  allocate(shamp(-dims(1):dims(1)-1,-dims(2):dims(2)-1,-dims(3):dims(3)-1))
  call shape%polyhedron_shapeamplitude(shamp, dims, dk, nthr)

else   ! the shape is not a polyhedron 

! first generate the shape function 
  allocate(sf(-dims(1):dims(1),-dims(2):dims(2),-dims(3):dims(3)))
  allocate(shamp(-dims(1):dims(1)-1,-dims(2):dims(2)-1,-dims(3):dims(3)-1))

  select case (trim(emnl%shapetype))
    case('sphere') 
      call Message%printMessage(' Generating the Sphere Shape Function')
      call self%ellipsoid_shapefunction_( sf, dims, dxyz, sngl((/ emnl%geometry(1), emnl%geometry(1), emnl%geometry(1) /) ))
      call Message%printMessage(' Generating the Sphere Shape Amplitude')
      call self%ellipsoid_shapeamplitude_( shamp, dims, dk, sngl((/ emnl%geometry(1), emnl%geometry(1), emnl%geometry(1) /) ))
      analytical = .TRUE.
    case('ellipsoid')
      call Message%printMessage(' Generating the Ellipsoid Shape Function')
      call self%ellipsoid_shapefunction_( sf, dims, dxyz, sngl(emnl%geometry(1:3) ))
      call Message%printMessage(' Generating the Ellipsoid Shape Amplitude')
      call self%ellipsoid_shapeamplitude_( shamp, dims, dk, sngl(emnl%geometry(1:3) ))
      analytical = .TRUE.
    case('prism')
      call Message%printMessage(' Generating the Prism Shape Function')
      call self%prism_shapefunction_( sf, dims, dxyz, sngl(emnl%geometry(1:3) ))
      call Message%printMessage(' Generating the Prism Shape Amplitude')
      call self%prism_shapeamplitude_( shamp, dims, dk, sngl(emnl%geometry(1:3) ))
      analytical = .TRUE.
    case('cylinder')
      call Message%printMessage(' Generating the Cylinder Shape Function')
      call self%ellipticcylinder_shapefunction_( sf, dims, dxyz, sngl((/ emnl%geometry(1), emnl%geometry(1), emnl%geometry(2) /) ))
      call Message%printMessage(' Generating the Cylinder  Shape Amplitude')
      call self%ellipticcylinder_shapeamplitude_( shamp, dims, dk, sngl((/ emnl%geometry(1), emnl%geometry(1), emnl%geometry(2)/)))
      analytical = .TRUE.
    case('ellipticcylinder')
      call Message%printMessage(' Generating the Elliptic Cylinder Shape Function')
      call self%ellipticcylinder_shapefunction_( sf, dims, dxyz, sngl(emnl%geometry(1:3) ))
      call Message%printMessage(' Generating the Elliptic Cylinder  Shape Amplitude')
      call self%ellipticcylinder_shapeamplitude_( shamp, dims, dk, sngl(emnl%geometry(1:3) ))
      analytical = .TRUE.
    case('torus')
      call Message%printMessage(' Generating the Torus Shape Function')
      call self%torus_shapefunction_( sf, dims, dxyz, sngl(emnl%geometry(1:2) ))
      call Message%printMessage(' Generating the Torus Shape Amplitude')
      ! call self%torus_shapeamplitude_( sf, shamp, dims, dk )
      analytical = .FALSE.
    case default 
      call Message%printError(' EMShapeAmplitude', 'Unknown shapetype: '//trim(emnl%shapetype))
  end select




end if 

! get the abs value for the STL file if needed
  if (trim(emnl%STLFilename).ne.'undefined') then
    allocate(shampreal(2*dims(1)+1,2*dims(2)+1,2*dims(3)+1))
    do i=1,2*dims(1)
        do j=1,2*dims(2)
            do k=1,2*dims(3)
                shampreal(i,j,k) = abs(shamp(i-dims(1)-1,j-dims(2)-1,k-dims(3)-1))
            end do 
        end do 
    end do 
  end if 

! shape intensity needed ? 
if (emnl%shapeIntensity.eqv..TRUE.) then 
  call Message%printMessage(' Computing the Shape Intensity')
  allocate(shint(2*dims(1)+1,2*dims(2)+1,2*dims(3)+1))
  shint = 0.D0
  do i=1,2*dims(1)
      do j=1,2*dims(2)
          do k=1,2*dims(3)
              shint(i,j,k) = abs(shamp(i-dims(1)-1,j-dims(2)-1,k-dims(3)-1))**2
          end do 
      end do 
  end do 
end if 

!====================================
!====================================
!====================================
! generate an HDF5 file with all the necessary arrays ... 
call openFortranHDFInterface()
HDF = HDF_T()

timer = Timing_T()
tstrb = timer%getTimeString()
dstr = timer%getDateString()

! Create a new file using the default properties.
datafile = EMsoft%generateFilePath('EMdatapathname', emnl%shampFilename)

hdferr =  HDF%createFile(datafile)
if (hdferr.ne.0) call HDF%error_check('HDF_createFile ', hdferr)

dataset = SC_Manufacturer
line2(1) = 'EMsoftOO'
line2(1) = cstringify(line2(1))
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)

! write the EMheader to the file
datagroupname = trim(HDFnames%get_ProgramData()) 
call HDF%writeEMheader(EMsoft,dstr, tstrb, tstre, progname, datagroupname)

! create a namelist group to write all the namelist files into
hdferr = HDF%createGroup(HDFnames%get_NMLfiles())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup NMLfiles', hdferr)

call HDF%pop()

! create a NMLparameters group to write all the namelist entries into
hdferr = HDF%createGroup(HDFnames%get_NMLparameters())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup NMLparameters', hdferr)

! and leave this group
call HDF%pop()

! then the remainder of the data in a EMData group
hdferr = HDF%createGroup(HDFnames%get_EMData())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup EMData', hdferr)

! create the sub group and add a HDF_FileVersion attribute to it
hdferr = HDF%createGroup(datagroupname)
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup ShapeAmplitude', hdferr)
HDF_FileVersion = '4.1'
attributename = SC_HDFFileVersion
hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)
 
dataset = SC_ShapeFunction
hdferr = HDF%writeDatasetFloatArray(dataset, sf, 2*dims(1)+1, 2*dims(2)+1, 2*dims(3)+1)

dataset = SC_ShapeAmplitudeReal
hdferr = HDF%writeDatasetDoubleArray(dataset, real(shamp), 2*dims(1), 2*dims(2), 2*dims(3))

dataset = SC_ShapeAmplitudeImaginary
hdferr = HDF%writeDatasetDoubleArray(dataset, aimag(shamp), 2*dims(1), 2*dims(2), 2*dims(3))

if (emnl%shapeIntensity.eqv..TRUE.) then 
  dataset = SC_ShapeIntensity
  hdferr = HDF%writeDatasetFloatArray(dataset, shint, 2*dims(1)+1, 2*dims(2)+1, 2*dims(3)+1)
end if 

call HDF%pop(.TRUE.)
call closeFortranHDFInterface()

!====================================
!====================================
!====================================
! do we need to create STL files?
if (trim(emnl%STLFilename).ne.'undefined') then
  header = emnl%STLheader
  stlname = EMsoft%generateFilePath('EMdatapathname',trim(emnl%STLFilename))
  MCA = MCA_T()

! first the shape function 
  sname = trim(stlname)//'_shapefunction.stl'  
  call Message%printMessage(' Generating '//trim(sname))
  if (emnl%isovalue.ne.0.D0) then   ! normalize the array to range [0,1]
    iso = emnl%isovalue 
  else
    iso = 0.50
  end if 
  call MCA%doMCA( sf, 2*dims+1, sngl(emnl%dxyz), iso )

  ntriangles = MCA%getNtriangles()
  STL = STL_T(sname, header, ntriangles, MCAlist=MCA%getMCAptr()) 

! then the shape amplitude 
  MCA = MCA_T()
  sname = trim(stlname)//'_shapeamplitude.stl'  
  call Message%printMessage(' Generating '//trim(sname))
  if (emnl%isovalue.ne.0.D0) then   ! normalize the array to range [0,1]
    iso = emnl%isovalue 
    mi = minval(shampreal)
    ma = maxval(shampreal)
    shampreal = (shampreal - mi) / (ma-mi)
  else
    iso = 0.10 * sngl(self%get_volume())
  end if 

write (*,*) 'shamp range ', minval(shampreal), maxval(shampreal), iso

  call MCA%doMCA( shampreal, 2*dims+1, sngl(emnl%dk), iso )

  ntriangles = MCA%getNtriangles()
  STL = STL_T(sname, header, ntriangles, MCAlist=MCA%getMCAptr()) 

! finally the shape intensity if needed 
  if (emnl%shapeIntensity.eqv..TRUE.) then 
    MCA = MCA_T()
    sname = trim(stlname)//'_shapeintensity.stl'  
    call Message%printMessage(' Generating '//trim(sname))
    if (emnl%isovalue.ne.0.D0) then   ! normalize the array to range [0,1]
      iso = emnl%isovalue 
      mi = minval(shint)
      ma = maxval(shint)
      shint = (shint - mi) / (ma-mi)
    else
      iso = 0.01 * sngl(self%get_volume())**2
    end if 
    call MCA%doMCA( shint, 2*dims+1, sngl(emnl%dk), iso )

    ntriangles = MCA%getNtriangles()
    STL = STL_T(sname, header, ntriangles, MCAlist=MCA%getMCAptr()) 
  end if 
end if 

end associate 

end subroutine ShapeAmplitude_

!--------------------------------------------------------------------------
subroutine readShapeAmplitude_(self, shfname, shamp) ! , shampnml)
!DEC$ ATTRIBUTES DLLEXPORT :: readShapeAmplitude_
!! author: MDG 
!! version: 1.0 
!! date: 08/16/21
!!
!! read an existing shape amplitude from an HDF5 file

use mod_EMsoft
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_io
use mod_HDFnames 
use stringconstants
use mod_EMsoft

IMPLICIT NONE

class(ShapeAmplitude_T), INTENT(INOUT)          :: self        
character(fnlen),INTENT(IN)                     :: shfname 
complex(kind=dbl),INTENT(INOUT),allocatable     :: shamp(:,:,:)
! type(ShapeAmplitudeNameListType),INTENT(INOUT)  :: shampnml 

type(HDF_T)                                     :: HDF
type(IO_T)                                      :: Message
type(HDFnames_T)                                :: HDFnames
type(EMsoft_T)                                  :: EMsoft
character(fnlen)                                :: datafile, dataset, tmpnmlname, p
integer(kind=irg)                               :: hdferr, nlines, i, d, io_int(3) 
integer(HSIZE_T)                                :: sz(1), dims(3)
logical                                         :: fexists
real(kind=dbl),allocatable                      :: shampr(:,:,:), shampi(:,:,:)
character(fnlen, KIND=c_char),allocatable,TARGET:: stringarray(:)

p = ''
EMsoft = EMsoft_T(p, p, silent=.TRUE.)

call openFortranHDFInterface()
HDF = HDF_T()

! Open an existing file (make sure it exists)
datafile = EMsoft%generateFilePath('EMdatapathname', trim(shfname))
inquire(file=trim(datafile),exist=fexists)

if (.not.fexists) then 
  call Message%printMessage(' looking for file '//trim(datafile))
  call Message%printError('readShapeAmplitude_','Shape amplitude file does not exist')
end if 

hdferr =  HDF%openFile(datafile, readonly=.TRUE.)
if (hdferr.ne.0) call HDF%error_check('HDF_openFile ', hdferr)

HDFnames = HDFnames_T() 
call HDFnames%set_ProgramData(SC_ShapeAmplitude) 
call HDFnames%set_NMLfilename(SC_SHAMPNML)

hdferr = HDF%openGroup(HDFnames%get_EMData())
if (hdferr.ne.0) call HDF%error_check('HDF_openGroup EMData', hdferr)

hdferr = HDF%openGroup(HDFnames%get_ProgramData())
if (hdferr.ne.0) call HDF%error_check('HDF_openGroup EMProgramData', hdferr)

dataset = SC_ShapeAmplitudeReal
call HDF%readDatasetDoubleArray(dataset, dims, hdferr, shampr)

dataset = SC_ShapeAmplitudeImaginary
call HDF%readDatasetDoubleArray(dataset, dims, hdferr, shampi)

io_int = dims 
call Message%writeValue(' Found shape amplitude array of size ',io_int,3)

! put the origin at the point (1,1,1)
d = dims(1)/2
shampr = cshift(shampr,d,1)
shampr = cshift(shampr,d,2)
shampr = cshift(shampr,d,3)
shampi = cshift(shampi,d,1)
shampi = cshift(shampi,d,2)
shampi = cshift(shampi,d,3)

allocate(shamp(dims(1),dims(2),dims(3)))
shamp = cmplx(shampr,shampi)
deallocate(shampr, shampi)

write (*,*) 'Shape Amplitude value in point 1,1,1 : ', shamp(1,1,1)

call HDF%pop()
call HDF%pop()

!!!! the following code is correct but for some reason the gfortran compiler
!!!! gets confused about module dependencies between this module and mod_demag

! ! next we get the namelist from the file
! hdferr = HDF%openGroup(HDFnames%get_NMLfiles())
! dataset = trim(HDFnames%get_NMLfilename())
! call HDF%readdatasetstringarray(dataset, nlines, hdferr, stringarray)
! sz = shape(stringarray)
! tmpnmlname = trim(EMsoft%generateFilePath('EMtmppathname'))//'tmp.nml'
! open(unit=65,file=trim(tmpnmlname),status='unknown',form='formatted')
! do i=1,sz(1)
!   write (65,"(A)") trim(stringarray(i))
! end do
! close(unit=65,status='keep')
! call self%readNameList(tmpnmlname)
! shampnml = self%nml 

! ! delete the tmp file
! open(unit=65,file=trim(tmpnmlname),status='unknown',form='formatted')
! close(unit=65,status='delete')

call HDF%pop(.TRUE.)
call closeFortranHDFInterface()

call Message%printMessage(' Read shape amplitude from file '//trim(datafile))

end subroutine readShapeAmplitude_

! !--------------------------------------------------------------------------
! subroutine readShapeFunction_(self, shfname, shfunc, shampnml)
! !DEC$ ATTRIBUTES DLLEXPORT :: readShapeFunction_
! !! author: MDG 
! !! version: 1.0 
! !! date: 08/16/21
! !!
! !! read an existing shape function from an HDF5 file

! use mod_EMsoft
! use HDF5
! use mod_HDFsupport
! use mod_HDFnames


! end subroutine readShapeFunction_

end module mod_ShapeAmplitude