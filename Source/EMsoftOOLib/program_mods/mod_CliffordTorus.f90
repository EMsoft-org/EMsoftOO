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

module mod_CliffordTorus
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/28/22
  !!
  !! class definition for the EMCliffordTorus program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMCliffordTorus program
type, public :: CliffordTorusNameListType
  integer(kind=irg)       :: reducetoRFZ 
  integer(kind=irg)       :: symmetrize
  integer(kind=irg)       :: shownegativeq0
  integer(kind=irg)       :: n
  integer(kind=irg)       :: pgnum
  character(fnlen)        :: anglefile
  character(fnlen)        :: sqtfile
  character(fnlen)        :: zpfile 
end type CliffordTorusNameListType

! class definition
type, public :: CliffordTorus_T
private 
  character(fnlen)       :: nmldeffile = 'EMCliffordTorus.nml'
  type(CliffordTorusNameListType)  :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: CliffordTorus_
  procedure, pass(self) :: setreducetoRFZ_
  procedure, pass(self) :: getreducetoRFZ_
  procedure, pass(self) :: setsymmetrize_
  procedure, pass(self) :: getsymmetrize_
  procedure, pass(self) :: setshownegativeq0_
  procedure, pass(self) :: getshownegativeq0_
  procedure, pass(self) :: setn_
  procedure, pass(self) :: getn_
  procedure, pass(self) :: setpgnum_
  procedure, pass(self) :: getpgnum_
  procedure, pass(self) :: setanglefile_
  procedure, pass(self) :: getanglefile_
  procedure, pass(self) :: setsqtfile_
  procedure, pass(self) :: getsqtfile_
  procedure, pass(self) :: setzpfile_
  procedure, pass(self) :: getzpfile_
  procedure, pass(self) :: createZonePlate_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: CliffordTorus => CliffordTorus_
  generic, public :: setreducetoRFZ => setreducetoRFZ_
  generic, public :: getreducetoRFZ => getreducetoRFZ_
  generic, public :: setsymmetrize => setsymmetrize_
  generic, public :: getsymmetrize => getsymmetrize_
  generic, public :: setshownegativeq0 => setshownegativeq0_
  generic, public :: getshownegativeq0 => getshownegativeq0_
  generic, public :: setpgnum => setpgnum_
  generic, public :: getpgnum => getpgnum_
  generic, public :: setn => setn_
  generic, public :: getn => getn_
  generic, public :: setanglefile => setanglefile_
  generic, public :: getanglefile => getanglefile_
  generic, public :: setsqtfile => setsqtfile_
  generic, public :: getsqtfile => getsqtfile_
  generic, public :: setzpfile => setzpfile_
  generic, public :: getzpfile => getzpfile_
  generic, public :: createZonePlate => createZonePlate_

end type CliffordTorus_T

! the constructor routine for this class 
interface CliffordTorus_T
  module procedure CliffordTorus_constructor
end interface CliffordTorus_T

contains

!--------------------------------------------------------------------------
type(CliffordTorus_T) function CliffordTorus_constructor( nmlfile ) result(CliffordTorus)
!! author: MDG 
!! version: 1.0 
!! date: 12/28/22
!!
!! constructor for the CliffordTorus_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

if (present(nmlfile)) call CliffordTorus%readNameList(nmlfile)

end function CliffordTorus_constructor

!--------------------------------------------------------------------------
subroutine CliffordTorus_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 12/28/22
!!
!! destructor for the CliffordTorus_T Class
 
IMPLICIT NONE

type(CliffordTorus_T), INTENT(INOUT)  :: self 

call reportDestructor('CliffordTorus_T')

end subroutine CliffordTorus_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 12/28/22
!!
!! read the namelist from an nml file for the CliffordTorus_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(CliffordTorus_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg)       :: reducetoRFZ 
integer(kind=irg)       :: symmetrize
integer(kind=irg)       :: shownegativeq0
integer(kind=irg)       :: n
integer(kind=irg)       :: pgnum
character(fnlen)        :: anglefile
character(fnlen)        :: sqtfile
character(fnlen)        :: zpfile 

namelist  / CliffordTorus / reducetoRFZ, symmetrize, shownegativeq0, n, pgnum, anglefile, sqtfile, zpfile

! set the input parameters to default values
anglefile = 'undefined' 
reducetoRFZ = 1
symmetrize = 0 
shownegativeq0 = 0
n = 500
pgnum = 32
sqtfile = 'undefined'
zpfile = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=CliffordTorus)
    close(UNIT=dataunit,STATUS='keep')

! check for required entries
    if (trim(anglefile).eq.'undefined') then
        call Message%printError('readNameList:',' anglefile file name is undefined in '//nmlfile)
    end if

    if ( (trim(sqtfile).eq.'undefined').and.(trim(zpfile).eq.'undefined') ) then
        call Message%printError('readNameList:',' at least one of sqtfile and zpfile need to be defined in '//nmlfile)
    end if
 end if

self%nml%anglefile = anglefile
self%nml%reducetoRFZ = reducetoRFZ
self%nml%symmetrize = symmetrize
self%nml%shownegativeq0 = shownegativeq0
self%nml%n = n
self%nml%pgnum = pgnum
self%nml%sqtfile = sqtfile
self%nml%zpfile = zpfile

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 12/28/22
!!
!! pass the namelist for the CliffordTorus_T Class to the calling program

IMPLICIT NONE 

class(CliffordTorus_T), INTENT(INOUT)          :: self
type(CliffordTorusNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine setreducetoRFZ_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setreducetoRFZ_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! set reducetoRFZ in the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%reducetoRFZ = inp

end subroutine setreducetoRFZ_

!--------------------------------------------------------------------------
function getreducetoRFZ_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getreducetoRFZ_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! get reducetoRFZ from the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%reducetoRFZ

end function getreducetoRFZ_

!--------------------------------------------------------------------------
subroutine setsymmetrize_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setsymmetrize_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! set symmetrize in the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%symmetrize = inp

end subroutine setsymmetrize_

!--------------------------------------------------------------------------
function getsymmetrize_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getsymmetrize_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! get symmetrize from the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%symmetrize

end function getsymmetrize_

!--------------------------------------------------------------------------
subroutine setshownegativeq0_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setshownegativeq0_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! set shownegativeq0 in the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%shownegativeq0 = inp

end subroutine setshownegativeq0_

!--------------------------------------------------------------------------
function getshownegativeq0_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getshownegativeq0_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! get shownegativeq0 from the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%shownegativeq0

end function getshownegativeq0_

!--------------------------------------------------------------------------
subroutine setn_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setn_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! set n in the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%n = inp

end subroutine setn_

!--------------------------------------------------------------------------
function getn_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getn_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! get n from the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%n

end function getn_

!--------------------------------------------------------------------------
subroutine setpgnum_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setpgnum_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! set pgnum in the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%pgnum = inp

end subroutine setpgnum_

!--------------------------------------------------------------------------
function getpgnum_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getpgnum_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! get pgnum from the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%pgnum

end function getpgnum_

!--------------------------------------------------------------------------
subroutine setanglefile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setanglefile_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! set anglefile in the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%anglefile = trim(inp)

end subroutine setanglefile_

!--------------------------------------------------------------------------
function getanglefile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getanglefile_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! get anglefile from the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%anglefile)

end function getanglefile_

!--------------------------------------------------------------------------
subroutine setsqtfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setsqtfile_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! set sqtfile in the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%sqtfile = trim(inp)

end subroutine setsqtfile_

!--------------------------------------------------------------------------
function getsqtfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getsqtfile_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! get sqtfile from the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%sqtfile)

end function getsqtfile_

!--------------------------------------------------------------------------
subroutine setzpfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setzpfile_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! set zpfile in the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%zpfile = trim(inp)

end subroutine setzpfile_

!--------------------------------------------------------------------------
function getzpfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getzpfile_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! get zpfile from the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%zpfile)

end function getzpfile_

!--------------------------------------------------------------------------
recursive subroutine createZonePlate_(self, EMsoft, SO, listmode)
!DEC$ ATTRIBUTES DLLEXPORT :: createZonePlate_
  !! author: MDG
  !! version: 1.0
  !! date: 12/21/22
  !!
  !! Create an orientation "zone plate", as described in DOI 10.1109/CVPR52688.2022.00811

use mod_quaternions
use mod_rotations
use mod_EMsoft
use h5im
use h5lt
use mod_image
use ISO_C_BINDING
use mod_io
use mod_so3

use, intrinsic :: iso_fortran_env

IMPLICIT NONE

class(CliffordTorus_T),INTENT(INOUT)   :: self
type(EMsoft_T),INTENT(INOUT)    :: EMsoft
type(so3_T),INTENT(INOUT)       :: SO
character(2),INTENT(IN)         :: listmode

type(q_T)                       :: q 
type(IO_T)                      :: Message

integer(kind=irg)               :: i, j, k, cnt, num, w, nn, offset, ixx, iyy, px, py, VZcnt
type(FZpointd), pointer         :: FZtmp, FZviz, VZlist, VZtmp
real(kind=dbl),parameter        :: s2 = 1.D0/sqrt(2.D0), r(4) = (/ s2, 0.D0, s2, 0.D0 /)
real(kind=dbl),allocatable      :: qu(:,:), z1(:), z2(:), h(:,:), h2(:,:), xx(:,:), yy(:,:), l(:), g(:,:)
real(kind=dbl)                  :: x(4), d, d1, d2, dx, dy, ss, zx, ee, kk 
character(fnlen)                :: vizname, fname

! declare variables for use in object oriented image module
integer                         :: iostat
character(len=128)              :: iomsg
logical                         :: isInteger
type(image_t)                   :: im
integer(int8)                   :: i8 (3,4)
integer(int8), allocatable      :: TIFF_image(:,:)

! output image size
nn = self%nml%n
num = 2*nn+1
w = 5
offset = w+(num-1)/2
ss = 0.25D0 * dble(self%nml%n) / 500.D0    ! initial plots were made on a 1001x1001 grid
kk = 40.D0 * dble(self%nml%n) / 500.D0

! get the list of orientations 
select case(listmode)
  case('FZ')
    FZtmp => SO%getListHead('FZ' )
    cnt = SO%getListCount('FZ')
  case('CM')
    FZtmp => SO%getListHead('CM' )
    cnt = SO%getListCount('CM')
  case('CO')
    FZtmp => SO%getListHead('CO' )
    cnt = SO%getListCount('CO')
  case('FB')
    FZtmp => SO%getListHead('FB' )
    cnt = SO%getListCount('FB')
  case('SF')
    FZtmp => SO%getListHead('SF' )
    cnt = SO%getListCount('SF')
  case('MA')
    FZtmp => SO%getListHead('MA' )
    cnt = SO%getListCount('MA')
  case('UN')
    FZtmp => SO%getListHead('UN' )
    cnt = SO%getListCount('UN')
  case default
    FZtmp => SO%getListHead('FZ' )
    cnt = SO%getListCount('FZ')
end select

! qu will hold all the orientation quaternions projected onto the Clifford torus
if (self%nml%shownegativeq0.eq.1) then 
  allocate(qu(4,2*cnt), z1(2*cnt), z2(2*cnt))
else
  allocate(qu(4,cnt), z1(cnt), z2(cnt))
end if 

allocate(VZlist)
VZtmp => VZlist
nullify(VZtmp%next)

do i=1,cnt
  q = FZtmp%qu
  x = q%q_copyd()
  d = sqrt(x(1)**2+x(2)**2)
  if (d.eq.0.D0) then
    d1 = 1.D0 
  else 
    d1 = 1.D0/d
  end if 
  d = sqrt(x(3)**2+x(4)**2)
  if (d.eq.0.D0) then
    d2 = 1.D0 
  else 
    d2 = 1.D0/d
  end if 
  qu(1:4,i) = (/ x(1)*d1, x(2)*d1, x(3)*d2, x(4)*d2 /) * s2 
! store the Clifford torus projections into a linked list for visualization purposes
  VZtmp%qu = q_T( qdinp = qu(1:4,i) )
  allocate(VZtmp%next)
  VZtmp => VZtmp%next
  if (self%nml%shownegativeq0.eq.1) then 
    qu(1:4,i+cnt) = - qu(1:4, i)
    VZtmp%qu = q_T( qdinp = qu(1:4,i) )
    allocate(VZtmp%next)
    VZtmp => VZtmp%next
  end if 
  nullify(VZtmp%next)
  FZtmp => FZtmp%next 
end do 

if (self%nml%shownegativeq0.eq.1) then 
  VZcnt = 2*cnt 
else
  VZcnt = cnt 
end if 

! compute the arc-tangent coordinates by projecting the Clifford torus onto a square
do i=1,VZcnt 
  z1(i) = atan2(qu(2,i), qu(1,i))
  z2(i) = atan2(qu(4,i), qu(3,i))
end do 

! rescale the arctangent coordinates to the output grid
z1 = z1*dble(nn)/cPi
z2 = z2*dble(nn)/cPi

! allocate the main padded output arrays
allocate( h(num+2*w,num+2*w), h2(num+2*w,num+2*w) )
allocate( xx(2*w+1, 2*w+1), yy(2*w+1, 2*w+1), l(2*w+1), g(2*w+1, 2*w+1) )
do i=1,2*w+1
  l(i) = dble(i-1-w)
end do 
do i=1,2*w+1 
  yy(i,:) = l(:)
  xx(:,i) = l(:)
end do 

! and fill the h array to obtain the zone plate
h = 0.D0
do i=1,VZcnt 
  ixx = int(z1(i))
  iyy = int(z2(i))
  dx = z1(i) - dble(ixx)
  dy = z2(i) - dble(iyy)
  zx = 0.5D0 * (1.D0 + cos( kk * acos( sum( qu(1:4,i) * r(1:4) ))**2 ))
  do j=1,2*w+1
    px = offset+ixx-w-1
    do k=1,2*w+1
      py = offset+iyy-w-1 
      ee = exp(-( (xx(j,k)-dx)**2 + (yy(j,k)-dy)**2 ) * ss )
      h(px+j,py+k) = h(px+j,py+k) + zx * ee
      h2(px+j,py+k) = h2(px+j,py+k) + ee
    end do 
  end do 
end do 

! and finally prepare for tiff output
if (trim(self%nml%zpfile).ne.'undefined') then
  h = h-minval(h)
  h = h / maxval(h)
  h = h*255.D0

  allocate(TIFF_image(num,num))
  do j=1,num
   do i=1,num
    TIFF_image(i,j) = int(h(w+i,w+j))
   end do
  end do

  fname = EMsoft%generateFilePath('EMdatapathname',self%nml%zpfile)

  ! set up the image_t structure
  im = image_t(TIFF_image)
  if(im%empty()) call Message%printMessage("createZonePlate_","failed to convert array to image")

  ! create the file
  call im%write(trim(fname), iostat, iomsg) ! format automatically detected from extension
  if(0.ne.iostat) then
    call Message%printMessage("failed to write image to file : "//iomsg)
  else
    call Message%printMessage('orientation zone plate written to '//trim(fname))
  end if
  deallocate(TIFF_image)
end if 

if (trim(self%nml%sqtfile).ne.'undefined') then
  h2 = h2-minval(h2)
  h2 = h2 / maxval(h2)
  h2 = h2*255.D0

  allocate(TIFF_image(num,num))
  do j=1,num
   do i=1,num
    TIFF_image(i,j) = int(h2(w+i,w+j))
   end do
  end do

  fname = EMsoft%generateFilePath('EMdatapathname',self%nml%sqtfile)

  ! set up the image_t structure
  im = image_t(TIFF_image)
  if(im%empty()) call Message%printMessage("createZonePlate_","failed to convert array to image")

  ! create the file
  call im%write(trim(fname), iostat, iomsg) ! format automatically detected from extension
  if(0.ne.iostat) then
    call Message%printMessage("failed to write image to file : "//iomsg)
  else
    call Message%printMessage('orientation zone plate written to '//trim(fname))
  end if
  deallocate(TIFF_image)
end if

end subroutine createZonePlate_

!--------------------------------------------------------------------------
subroutine CliffordTorus_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: CliffordTorus_
!! author: MDG 
!! version: 1.0 
!! date: 12/28/22
!!
!! perform the computations

use mod_EMsoft
use mod_so3
use mod_io

IMPLICIT NONE 

class(CliffordTorus_T), INTENT(INOUT)   :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 

type(so3_T)                             :: SO
type(IO_T)                              :: Message 

character(fnlen)                        :: oname 

! initialize the SO3 class
SO = so3_T( self%nml%pgnum, zerolist='FZ' )

! read all the orientations from the anglefile
oname = EMsoft%generateFilePath('EMdatapathname',self%nml%anglefile)
call Message%printMessage(' Reading orientations from file '//trim(oname))
call SO%getOrientationsfromFile( oname )

! we have the list, so what do we need to do with it before computing the Clifford Torus representation ?

! to be written [12/28/22]


! generate the output images using the Clifford Torus sprojection
call Message%printMessage(' Generating square torus/zone plate representation(s)')
call self%createZonePlate_(EMsoft, SO, listmode='FZ' )

end subroutine CliffordTorus_

end module mod_CliffordTorus