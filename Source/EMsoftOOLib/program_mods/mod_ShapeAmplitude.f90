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
  character(10)         :: shapetype
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

contains
private 
  procedure, pass(self) :: get_polyEdgeL_
  procedure, pass(self) :: get_geometry_
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

  generic, public :: get_polyEdgeL => get_polyEdgeL_
  generic, public :: get_geometry => get_geometry_
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

call ShapeAmplitude%readNameList(nmlfile)

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
character(10)         :: shapetype
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
subroutine ShapeAmplitude_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: ShapeAmplitude_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! perform the computations

use mod_EMsoft
use mod_PGA3D 
use mod_PGA3Dsupport
use mod_polyhedra
use mod_STL
use mod_MCA 
use mod_IO

IMPLICIT NONE 

class(ShapeAmplitude_T), INTENT(INOUT)  :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 

type(STL_T)                             :: STL 
type(MCA_T)                             :: MCA 
type(polyhedron_T)                      :: shape     
type(PGA3D_T)                           :: mv
type(IO_T)                              :: Message

character(fnlen)                        :: polyname, stlname, sname 
character(80)                           :: header
real(kind=dbl)                          :: edgeL, dk 
real(kind=sgl)                          :: iso, mi, ma 
real(kind=sgl),allocatable              :: sf(:,:,:)
complex(kind=dbl),allocatable           :: shamp(:,:,:)
real(kind=sgl),allocatable              :: shampreal(:,:,:), shint(:,:,:)
integer(kind=irg)                       :: i, j, k, nthr, dims(3), ntriangles

associate(emnl => self%nml)

if (trim(emnl%shapetype).eq.'polyhedron') then 
  call Message%printMessage(' Initializing 3D Projective Geometric Algebra module')
  call PGA3D_initialize()
  ! check for / in the filename !!!
  polyname = trim(emnl%polyhedronFilename)
  edgeL = emnl%polyEdgeL 
  shape = polyhedron_T( polyname, edgeL )
  dims = emnl%dims

! first generate the shape function 
  call Message%printMessage(' Generating the Shape Function')
  allocate(sf(-dims(1):dims(1),-dims(2):dims(2),-dims(3):dims(3)))
  call shape%polyhedron_shapefunction(sf, dims, emnl%dxyz)

! then the shape amplitude using the Komrska approach
  call Message%printMessage(' Generating the Shape Amplitude')
  allocate(shamp(-dims(1):dims(1)-1,-dims(2):dims(2)-1,-dims(3):dims(3)-1))
  dk = emnl%dk
  nthr = emnl%nthreads
  call shape%polyhedron_shapeamplitude(shamp, dims, dk, nthr)
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
  if (emnl%isovalue.ne.0.D0) then   ! normalize the array to range [0,1]
    mi = minval(shampreal)
    ma = maxval(shampreal)
    shampreal = (shampreal - mi) / (ma-mi)
  end if 

! shape intensity needed ? 
  if (emnl%shapeIntensity.eqv..TRUE.) then 
    call Message%printMessage(' Computing the Shape Intensity')
    allocate(shint(2*dims(1)+1,2*dims(2)+1,2*dims(3)+1))
    do i=1,2*dims(1)
        do j=1,2*dims(2)
            do k=1,2*dims(3)
                shint(i,j,k) = abs(shamp(i-dims(1)-1,j-dims(2)-1,k-dims(3)-1))**2
            end do 
        end do 
    end do 
  end if 
  if (emnl%isovalue.ne.0.D0) then   ! normalize the array to range [0,1]
    mi = minval(shint)
    ma = maxval(shint)
    shint = (shint - mi) / (ma-mi)
  end if 

else 

end if 

! generate an HDF5 file with all the necessary arrays ... 


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
  call MCA%doMCA( sf, 2*dims+1, sngl(emnl%dk), iso )

  ntriangles = MCA%getNtriangles()
  STL = STL_T(sname, header, ntriangles, MCAlist=MCA%getMCAptr()) 

! then the shape amplitude 
  MCA = MCA_T()
  sname = trim(stlname)//'_shapeamplitude.stl'  
  call Message%printMessage(' Generating '//trim(sname))
  if (emnl%isovalue.ne.0.D0) then   ! normalize the array to range [0,1]
    iso = emnl%isovalue 
  else
    iso = 0.10 * sngl(shape%get_volume())
  end if 
  call MCA%doMCA( shampreal, 2*dims+1, sngl(emnl%dxyz), iso )

  ntriangles = MCA%getNtriangles()
  STL = STL_T(sname, header, ntriangles, MCAlist=MCA%getMCAptr()) 

! finally the shape intensity if needed 
  if (emnl%shapeIntensity.eqv..TRUE.) then 
    MCA = MCA_T()
    sname = trim(stlname)//'_shapeintensity.stl'  
    call Message%printMessage(' Generating '//trim(sname))
    if (emnl%isovalue.ne.0.D0) then   ! normalize the array to range [0,1]
      iso = emnl%isovalue 
    else
      iso = 0.01 * sngl(shape%get_volume())**2
    end if 
    call MCA%doMCA( shint, 2*dims+1, sngl(emnl%dk), iso )

    ntriangles = MCA%getNtriangles()
    STL = STL_T(sname, header, ntriangles, MCAlist=MCA%getMCAptr()) 
  end if 
end if 

end associate 

end subroutine ShapeAmplitude_



end module mod_ShapeAmplitude