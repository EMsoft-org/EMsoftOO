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
  integer(kind=irg)       :: overlayRFZ 
  integer(kind=irg)       :: symmetrize
  integer(kind=irg)       :: shownegativeq0
  integer(kind=irg)       :: logarithmic
  integer(kind=irg)       :: doRiesz
  integer(kind=irg)       :: n
  integer(kind=irg)       :: pgnum
  character(fnlen)        :: hdffile
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
  procedure, pass(self) :: setoverlayRFZ_
  procedure, pass(self) :: getoverlayRFZ_
  procedure, pass(self) :: setreducetoRFZ_
  procedure, pass(self) :: getreducetoRFZ_
  procedure, pass(self) :: setsymmetrize_
  procedure, pass(self) :: getsymmetrize_
  procedure, pass(self) :: setshownegativeq0_
  procedure, pass(self) :: getshownegativeq0_
  procedure, pass(self) :: setlogarithmic_
  procedure, pass(self) :: getlogarithmic_
  procedure, pass(self) :: setdoRiesz_
  procedure, pass(self) :: getdoRiesz_
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
  procedure, pass(self) :: projectqtoCT_
  procedure, pass(self) :: convertqtoSquareTorus_
  procedure, pass(self) :: overlayRFZ_
  procedure, pass(self) :: makeSquareTorus_
  procedure, pass(self) :: createZonePlate_
  procedure, pass(self) :: generateTIFFfile_
  procedure, pass(self) :: computeRieszEnergy_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: CliffordTorus => CliffordTorus_
  generic, public :: setoverlayRFZ => setoverlayRFZ_
  generic, public :: getoverlayRFZ => getoverlayRFZ_
  generic, public :: setreducetoRFZ => setreducetoRFZ_
  generic, public :: getreducetoRFZ => getreducetoRFZ_
  generic, public :: setsymmetrize => setsymmetrize_
  generic, public :: getsymmetrize => getsymmetrize_
  generic, public :: setshownegativeq0 => setshownegativeq0_
  generic, public :: getshownegativeq0 => getshownegativeq0_
  generic, public :: setlogarithmic => setlogarithmic_
  generic, public :: getlogarithmic => getlogarithmic_
  generic, public :: setdoRiesz => setdoRiesz_
  generic, public :: getdoRiesz => getdoRiesz_
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
  generic, public :: projectqtoCT => projectqtoCT_
  generic, public :: convertqtoSquareTorus => convertqtoSquareTorus_
  generic, public :: createZonePlate => createZonePlate_
  generic, public :: overlayRFZ => overlayRFZ_
  generic, public :: makeSquareTorus => makeSquareTorus_
  generic, public :: generateTIFFfile => generateTIFFfile_
  generic, public :: computeRieszEnergy => computeRieszEnergy_ 

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
integer(kind=irg)       :: overlayRFZ 
integer(kind=irg)       :: symmetrize
integer(kind=irg)       :: shownegativeq0
integer(kind=irg)       :: logarithmic
integer(kind=irg)       :: doRiesz
integer(kind=irg)       :: n
integer(kind=irg)       :: pgnum
character(fnlen)        :: hdffile
character(fnlen)        :: anglefile
character(fnlen)        :: sqtfile
character(fnlen)        :: zpfile 

namelist  / CliffordTorus / reducetoRFZ, symmetrize, shownegativeq0, n, pgnum, anglefile, sqtfile, &
                            zpfile, doRiesz, overlayRFZ, logarithmic, hdffile

! set the input parameters to default values
anglefile = 'undefined' 
reducetoRFZ = 1
overlayRFZ = 0
symmetrize = 0 
shownegativeq0 = 0
logarithmic = 0
doRiesz = 0
n = 500
pgnum = 32
hdffile = 'undefined'
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
self%nml%hdffile = hdffile
self%nml%reducetoRFZ = reducetoRFZ
self%nml%overlayRFZ = overlayRFZ
self%nml%symmetrize = symmetrize
self%nml%shownegativeq0 = shownegativeq0
self%nml%logarithmic =logarithmic  
self%nml%doRiesz =doRiesz 
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
subroutine setoverlayRFZ_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setoverlayRFZ_
!! author: MDG
!! version: 1.0
!! date: 01/02/23
!!
!! set overlayRFZ in the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%overlayRFZ = inp

end subroutine setoverlayRFZ_

!--------------------------------------------------------------------------
function getoverlayRFZ_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getoverlayRFZ_
!! author: MDG
!! version: 1.0
!! date: 01/02/23
!!
!! get overlayRFZ from the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%overlayRFZ

end function getoverlayRFZ_

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
subroutine setlogarithmic_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setlogarithmic_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! set logarithmic in the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%logarithmic = inp

end subroutine setlogarithmic_

!--------------------------------------------------------------------------
function getlogarithmic_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getlogarithmic_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! get logarithmic from the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%logarithmic

end function getlogarithmic_

!--------------------------------------------------------------------------
subroutine setdoRiesz_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdoRiesz_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! set doRiesz in the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)             :: inp

self%nml%doRiesz = inp

end subroutine setdoRiesz_

!----------------doRiesz--------------------------------------------
function getdoRiesz_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdoRiesz_
!! author: MDG
!! version: 1.0
!! date: 12/28/22
!!
!! get doRiesz from the CliffordTorus_T class

IMPLICIT NONE

class(CliffordTorus_T), INTENT(INOUT)     :: self
integer(kind=irg)                         :: out

out = self%nml%doRiesz

end function getdoRiesz_

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
recursive subroutine computeRieszEnergy_(self, SO, listmode)
!DEC$ ATTRIBUTES DLLEXPORT :: computeRieszEnergy_
  !! author: MDG
  !! version: 1.0
  !! date: 12/31/22
  !!
  !! Compute the Riesz energies for this orientation set and compare them to the 
  !! theoretical optimal values

use mod_quaternions
use mod_rotations
use mod_io
use mod_so3
use omp_lib

IMPLICIT NONE

class(CliffordTorus_T),INTENT(INOUT)    :: self
type(so3_T),INTENT(INOUT)               :: SO
character(2),INTENT(IN)                 :: listmode

type(IO_T)                              :: Message
type(q_T)                               :: q

type(FZpointd), pointer                 :: FZtmp, FZviz, VZlist, VZtmp
integer(kind=irg)                       :: cnt, nthreads, i, j, TID, io_int(1)
real(kind=dbl)                          :: rieszoptimal(3), rieszsum(3), tsum(3), eu(3), y, r, thr
real(kind=dbl),allocatable              :: qar(:,:), qar2(:,:)
integer(kind=irg),allocatable           :: done(:)

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

if (self%nml%shownegativeq0.eq.1) then 
  allocate(qar(4,2*cnt))
else
  allocate(qar(4,cnt))
end if 

do i=1,cnt
  q = FZtmp%qu
  qar(1:4,i) = q%q_copyd()
  if (self%nml%shownegativeq0.eq.1) then 
    qar(1:4,i+cnt) = - qar(1:4, i)
  end if 
  FZtmp => FZtmp%next
end do 

if (self%nml%shownegativeq0.eq.1) cnt = 2*cnt

allocate(qar2(4,cnt), done(cnt))
qar2 = qar

! optimal values for the Riesz energies
rieszoptimal = (/ 4.D0*dble(cnt)**2/3.D0/cPi, dble(cnt)**2*0.5D0, dble(cnt)**2*dlog(dble(cnt))/3.D0/cPi /)

tsum = 0.D0
thr = 1.0
nthreads = OMP_GET_MAX_THREADS() 
call OMP_SET_NUM_THREADS(nthreads)
io_int(1) = nthreads
call Message%WriteValue(' number of threads : ', io_int, 1)
rieszsum = 0.D0
done = 0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tsum,j,y,TID)

TID = OMP_GET_THREAD_NUM()

! compute the Riesz energy
!$OMP DO SCHEDULE(DYNAMIC)
do i=1,cnt
  do j=i+1,cnt
    y = dsqrt( sum( (qar(1:4,i) - qar2(1:4,j))**2 ) )
    if (y.ne.0.D0) then 
      y = 1.D0 / dsqrt( sum( (qar(1:4,i) - qar2(1:4,j))**2 ) )
      tsum(1) = tsum(1) + y
      tsum(2) = tsum(2) + y**2
      tsum(3) = tsum(3) + y**3
    end if
  end do
  done(i)=1
  if (TID.eq.0) then
     r = 100.0*float(sum(done))/float(cnt)
     if (r.gt.thr) then 
       thr = thr+1.0
       io_int(1) = r
       call Message%WriteValue('percent completed : ', io_int, 1) 
     end if
  end if
  end do
!$OMP CRITICAL
  rieszsum = rieszsum + tsum
!$OMP END CRITICAL

!$OMP END PARALLEL

eu = rieszsum/rieszoptimal

write (*,"('Riesz optimal : ',3(F25.3,A1))") rieszoptimal(1),' ',rieszoptimal(2),' ',rieszoptimal(3)
write (*,"('Riesz energy  : ',3(F25.3,A1))") rieszsum(1),' ',rieszsum(2),' ',rieszsum(3)
write (*,"('Riesz ratio   : ',3(F25.20,A1))") eu(1),' ',eu(2),' ',eu(3)

end subroutine computeRieszEnergy_

!--------------------------------------------------------------------------
recursive subroutine overlayRFZ_(self, SO, num, XYZ, TIFF_image)
!DEC$ ATTRIBUTES DLLEXPORT :: overlayRFZ_
  !! author: MDG
  !! version: 1.0
  !! date: 01/02/23
  !!
  !! overlay the projection of the RFZ onto the zone plate image

use mod_so3
use mod_rotations
use mod_povray
use ISO_C_BINDING

use, intrinsic :: iso_fortran_env

IMPLICIT NONE

class(CliffordTorus_T),INTENT(INOUT)    :: self
type(so3_T),INTENT(INOUT)               :: SO
character(4),INTENT(IN)                 :: XYZ
integer(kind=irg),INTENT(IN)            :: num
integer(int8),INTENT(INOUT)             :: TIFF_image(num,num)

type(PoVRay_T)                          :: POV
type(r_T)                               :: ro1, ro2, ro
type(q_T)                               :: qu, q 

integer(kind=irg)                       :: dims(3), FZtype, FZorder, ns, nt, i, j, intXY(2), offset
real(kind=dbl)                          :: d, dx, tpi, hpi, aux(4), xx, rod(3), XY(2), scl
integer(kind=irg),allocatable           :: s_edge(:,:), t_edge(:,:)
real(kind=dbl),allocatable              :: cpos(:,:)
logical                                 :: twostep

call setRotationPrecision('Double')
tpi = 2.D0 * cPi
hpi = 0.5D0 * cPi
scl = dble((num-1)/2)/cPi
offset = (num-1)/2 + 1

call SO%getFZtypeandorder(FZtype, FZorder)

! there is no need to instantiate the PoVRay_T class; we only need to call 
! one of its methods which does not affect any parameters in the class 
! we need to get the coordinates and connectivity of the vertices of the RFZ

if (FZtype.eq.2) then
  if (FZorder.eq.6) then
    twostep = .TRUE.
    dims = (/ 24, 24, 12 /)
    allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
    call POV%getpos_FZ622(dims, cpos, s_edge, t_edge, ns, d, nt)
  end if
  if (FZorder.eq.4) then
    twostep = .TRUE.
    dims = (/ 16, 16, 8 /)
    allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
    call POV%getpos_FZ422(dims, cpos, s_edge, t_edge, ns, d, nt)
  end if
  if (FZorder.eq.3) then
    twostep = .TRUE.
    dims = (/ 12, 12, 6 /)
    allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
    call POV%getpos_FZ32(dims, cpos, s_edge, t_edge, ns, d, nt)
  end if
  if (FZorder.eq.2) then
    twostep = .TRUE.
    dims = (/ 8, 8, 4 /)
    allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
    call POV%getpos_FZ222(dims, cpos, s_edge, t_edge, ns, d, nt)
  end if
end if

if (FZtype.eq.3) then
! rotational group 23
  twostep = .FALSE.
  dims = (/ 6, 12, 1 /)
  allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
  call POV%getpos_FZ23(dims, cpos, s_edge, t_edge, ns, d, nt)
end if

if (FZtype.eq.4) then
! rotational group 432
  twostep = .TRUE.
  dims = (/ 24, 12, 24 /)
  allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
  call POV%getpos_FZ432(dims, cpos, s_edge, t_edge, ns, d, nt)
end if

! overlay the outline onto the TIFF_image at byte level 255

if (FZtype.ne.1) then 
  ! create the square edges first
   dx = 1.D0/dble(ns)
   do i=1,dims(2)
    ro1 = r_T( rdinp = (/ cpos(1:3,s_edge(1,i)), d/) )
    ro2 = r_T( rdinp = (/ cpos(1:3,s_edge(2,i)), d/) )
    do j=1,ns+1
      aux = d*ro1%r_copyd() + d*(ro2%r_copyd() - ro1%r_copyd()) * j * dx
      xx = dsqrt( sum (aux(1:3)**2) )
      ro = r_T( rdinp = (/ aux(1:3)/xx, xx /) )
  ! project this point onto the Clifford Torus
      q = ro%rq()
      qu = self%projectqtoCT_( q, XYZ )
  ! convert to Square Torus coordinates
      XY = self%convertqtoSquareTorus_( qu, 'XZ_Y')
  ! and draw this point on the zone plate 
      intXY = nint( XY * scl )
      TIFF_image( offset+intXY(1), offset+intXY(2) ) = -1_int8
    end do
   end do

   if (twostep) then
     dx = 1.D0/dble(nt)
     do i=1,dims(3)
      ro1 = r_T( rdinp = (/ cpos(1:3,t_edge(1,i)), d/) )
      ro2 = r_T( rdinp = (/ cpos(1:3,t_edge(2,i)), d/) )
      do j=1,nt+1
        aux = d*ro1%r_copyd() + d*(ro2%r_copyd() - ro1%r_copyd()) * j * dx
        xx = dsqrt( sum (aux(1:3)**2) )
        ro = r_T( rdinp = (/ aux(1:3)/xx, xx /) )
  ! project this point onto the Clifford Torus
        q = ro%rq()
        qu = self%projectqtoCT_( q, XYZ )
  ! convert to Square Torus coordinates
        XY = self%convertqtoSquareTorus_( qu, 'XZ_Y' )
  ! and draw this point on the zone plate 
        intXY = nint( XY * scl )
        TIFF_image( offset+intXY(1), offset+intXY(2) ) = -1_int8
      end do
     end do
    end if
else 

end if 

end subroutine overlayRFZ_

!--------------------------------------------------------------------------
recursive function projectqtoCT_(self, q, XYZ) result(qout)
!DEC$ ATTRIBUTES DLLEXPORT :: projectqtoCT_
  !! author: MDG
  !! version: 1.0
  !! date: 01/02/23
  !!
  !! Project a unit quaternion onto the Clifford Torus 

use mod_rotations

class(CliffordTorus_T),INTENT(INOUT)    :: self
type(q_T),INTENT(INOUT)                 :: q  
character(4),INTENT(IN)                 :: XYZ
type(q_T)                               :: qout

real(kind=dbl),parameter                :: s2 = 1.D0/sqrt(2.D0)
real(kind=dbl)                          :: x(4), d, d1, d2

call setRotationPrecision('Double')

select case(XYZ)
  case('XZ_Y')
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
    qout = q_T( qdinp = (/ x(1)*d1, x(2)*d1, x(3)*d2, x(4)*d2 /) * s2 )
  case('YX_Z')
    x = q%q_copyd()
    d = sqrt(x(1)**2+x(3)**2)
    if (d.eq.0.D0) then
      d1 = 1.D0 
    else 
      d1 = 1.D0/d
    end if 
    d = sqrt(x(2)**2+x(4)**2)
    if (d.eq.0.D0) then
      d2 = 1.D0 
    else 
      d2 = 1.D0/d
    end if 
    qout = q_T( qdinp = (/ x(1)*d1, x(3)*d1, x(4)*d2, x(2)*d2 /) * s2 )
  case('ZY_X')
    x = q%q_copyd()
    d = sqrt(x(1)**2+x(4)**2)
    if (d.eq.0.D0) then
      d1 = 1.D0 
    else 
      d1 = 1.D0/d
    end if 
    d = sqrt(x(2)**2+x(3)**2)
    if (d.eq.0.D0) then
      d2 = 1.D0 
    else 
      d2 = 1.D0/d
    end if 
    qout = q_T( qdinp = (/ x(1)*d1, x(4)*d1, x(2)*d2, x(3)*d2 /) * s2 )
end select 

end function projectqtoCT_

!--------------------------------------------------------------------------
recursive function convertqtoSquareTorus_(self, q, XYZ) result(XY)
!DEC$ ATTRIBUTES DLLEXPORT :: convertqtoSquareTorus_
  !! author: MDG
  !! version: 1.0
  !! date: 01/02/23
  !!
  !! Project a quaternion from the Clifford Torus to the 2-D Square Torus

use mod_rotations

class(CliffordTorus_T),INTENT(INOUT)    :: self
type(q_T),INTENT(INOUT)                 :: q  
character(4),INTENT(IN)                 :: XYZ
real(kind=dbl)                          :: XY(2)

real(kind=dbl)                          :: qu(4)

qu = q%q_copyd()

! which projection do we need ?
select case(XYZ)
  case('XZ_Y')
    XY(1) = atan2(qu(2), qu(1))
    XY(2) = atan2(qu(4), qu(3))
  case('YX_Z')
    XY(1) = atan2(qu(3), qu(1))
    XY(2) = atan2(qu(2), qu(4))
  case('ZY_X')
    XY(1) = atan2(qu(4), qu(1))
    XY(2) = atan2(qu(3), qu(2))
end select

end function convertqtoSquareTorus_

!--------------------------------------------------------------------------
recursive subroutine makeSquareTorus_(self, num, w, cnt, xx, yy, offset, qu, h, h2) !, XYZ) 
!DEC$ ATTRIBUTES DLLEXPORT :: makeSquareTorus_
  !! author: MDG
  !! version: 1.0
  !! date: 01/02/23
  !!
  !! generate the zone plate and square torus intensity arrays

use mod_rotations
use mod_io
use omp_lib

class(CliffordTorus_T),INTENT(INOUT)    :: self
integer(kind=irg),INTENT(IN)            :: num
integer(kind=irg),INTENT(IN)            :: w
integer(kind=irg),INTENT(IN)            :: cnt
real(kind=dbl),INTENT(IN)               :: xx(2*w+1, 2*w+1)
real(kind=dbl),INTENT(IN)               :: yy(2*w+1, 2*w+1)
integer(kind=irg),INTENT(IN)            :: offset
real(kind=dbl),INTENT(IN)               :: qu(5,cnt)
real(kind=dbl),INTENT(INOUT)            :: h(num+2*w,num+2*w)
real(kind=dbl),INTENT(INOUT)            :: h2(num+2*w,num+2*w)
! character(4),INTENT(IN)                 :: XYZ

type(q_T)                               :: q
type(IO_T)                              :: Message 

integer(kind=irg)                       :: i, j, k, nn, ixx, iyy, px, py, io_int(2), nthreads, TID  
real(kind=dbl)                          :: ss, kk, dx, dy, zx, ee, XY(2), z1(cnt), z2(cnt)
real(kind=dbl),parameter                :: s2 = 1.D0/sqrt(2.D0), r(4) = (/ s2, 0.D0, s2, 0.D0 /)
real(kind=dbl),allocatable              :: hlocal(:,:), h2local(:,:)

nn = self%nml%n
ss = 0.25D0 * dble(self%nml%n) / 500.D0    ! initial plots were made on a 1001x1001 grid
kk = 40.D0 * dble(self%nml%n) / 500.D0

! compute the arc-tangent coordinates by projecting the Clifford torus onto a square
do i=1,cnt 
  q = q_T( qdinp=qu(1:4,i) )
  XY = self%convertqtoSquareTorus_( q, 'XZ_Y')
  z1(i) = XY(1)
  z2(i) = XY(2)
end do 

! rescale the arctangent coordinates to the output grid
z1 = z1*dble(nn)/cPi
z2 = z2*dble(nn)/cPi

call Message%printMessage('  - adding orientations to zone plate ')
! and fill the h arrays to obtain the zone plate; we'll use parallel threads to do this...
nthreads = OMP_GET_MAX_THREADS()
call OMP_SET_NUM_THREADS(nthreads)
io_int(1) = nthreads
call Message%WriteValue(' number of threads : ', io_int, 1)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(hlocal, h2local, i, j, k, ixx, iyy, dx, dy, zx, px, py, ee, io_int, TID)

TID = OMP_GET_THREAD_NUM()

allocate( hlocal(num+2*w,num+2*w), h2local(num+2*w,num+2*w) )
hlocal = 0.D0
h2local = 0.D0

!$OMP DO SCHEDULE(DYNAMIC)
do i=1,cnt 
  ixx = int(z1(i))
  iyy = int(z2(i))
  dx = z1(i) - dble(ixx)
  dy = z2(i) - dble(iyy)
  zx = 0.5D0 * (1.D0 + cos( kk * acos( sum( qu(1:4,i) * r(1:4) ))**2 ))
  px = offset+ixx-w-1
  py = offset+iyy-w-1 
  if ((px.gt.0).and.(px.lt.num+2*w).and.(py.gt.0).and.(py.lt.num+2*w)) then
    do j=1,2*w+1
      do k=1,2*w+1
        ee = exp(-( (xx(j,k)-dx)**2 + (yy(j,k)-dy)**2 ) * ss )
        hlocal(px+j,py+k) = hlocal(px+j,py+k) + qu(5,i) * zx * ee
        h2local(px+j,py+k) = h2local(px+j,py+k) + qu(5,i) * ee
      end do 
    end do 
  end if 
  if (mod(i,1000000).eq.0) then 
    io_int(1) = i
    io_int(2) = TID
    call Message%WriteValue('  - current orientation #, TID ', io_int, 2)
  end if 
end do 

!$OMP CRITICAL
h = h + hlocal
h2 = h2 + h2local
!$OMP END CRITICAL

!$OMP END PARALLEL

end subroutine makeSquareTorus_

!--------------------------------------------------------------------------
recursive subroutine generateTIFFfile_(self, EMsoft, SO, num, w, h, mode, XYZ) 
!DEC$ ATTRIBUTES DLLEXPORT :: generateTIFFfile_
  !! author: MDG
  !! version: 1.0
  !! date: 01/02/23
  !!
  !! generate the zone plate and square torus TIFF files

use mod_EMsoft
use mod_image
use ISO_C_BINDING
use mod_io
use mod_so3

use, intrinsic :: iso_fortran_env

class(CliffordTorus_T),INTENT(INOUT)    :: self
type(EMsoft_T),INTENT(INOUT)            :: EMsoft
type(so3_T),INTENT(INOUT)               :: SO
integer(kind=irg),INTENT(IN)            :: num
integer(kind=irg),INTENT(IN)            :: w
real(kind=dbl),INTENT(INOUT)            :: h(num+2*w,num+2*w)
character(*),INTENT(IN)                 :: mode
character(4),INTENT(IN)                 :: XYZ

type(IO_T)                              :: Message

integer(kind=irg)                       :: i, j 
character(fnlen)                        :: str, fname

! declare variables for use in object oriented image module
integer                                 :: iostat
character(len=128)                      :: iomsg
logical                                 :: isInteger
type(image_t)                           :: im
integer(int8)                           :: i8 (3,4)
integer(int8), allocatable              :: TIFF_image(:,:)

if (trim(mode).eq.'SQT') then 
  str = trim(self%nml%sqtfile)//XYZ//'_SQT.tiff'
else
  str = trim(self%nml%zpfile)//XYZ//'_ZP.tiff'
end if 

h = h - minval(h)
h = h / maxval(h)
h = h * 255.D0

allocate(TIFF_image(num,num))
do j=1,num
 do i=1,num
  TIFF_image(i,num+1-j) = int(h(w+i,w+j))
 end do
end do

if (self%nml%overlayRFZ.eq.1) then 
  call self%overlayRFZ( SO, num, XYZ, TIFF_image )
end if 

fname = EMsoft%generateFilePath('EMdatapathname',str)

! set up the image_t structure
im = image_t(TIFF_image)
if(im%empty()) call Message%printMessage("createZonePlate_","failed to convert array to image")

! create the file
call im%write(trim(fname), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage("failed to write image to file : "//iomsg)
else
  call Message%printMessage(' - orientation zone plate written to '//trim(fname))
end if
deallocate(TIFF_image)

end subroutine generateTIFFfile_

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
use HDF5
use mod_HDFsupport
use mod_io
use mod_so3

IMPLICIT NONE

class(CliffordTorus_T),INTENT(INOUT)   :: self
type(EMsoft_T),INTENT(INOUT)    :: EMsoft
type(so3_T),INTENT(INOUT)       :: SO
character(2),INTENT(IN)         :: listmode

type(q_T)                       :: q 
type(IO_T)                      :: Message
type(HDF_T)                     :: HDF

integer(kind=irg)               :: i, j, k, cnt, num, w, nn, offset, io_int(1), hdferr
type(FZpointd), pointer         :: FZtmp, FZviz
real(kind=dbl),allocatable      :: qu(:,:), h(:,:), h2(:,:), xx(:,:), yy(:,:), l(:), g(:,:)
real(kind=dbl)                  :: x(4), d, d1, d2, dx, dy, ss, zx, ee, kk, logoffset 
character(fnlen)                :: vizname, fname, dataset, groupname 

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

! save the head pointer
FZviz => FZtmp

! qu will hold all the orientation quaternions projected onto the Clifford torus
! the fifth entry holds the weight factor
allocate(qu(5,cnt))

! do we need to export the raw data to an HDF5 file ?
if (trim(self%nml%hdffile).ne.'undefined') then 
  fname = EMsoft%generateFilePath('EMdatapathname',self%nml%hdffile)
  call openFortranHDFInterface()

  HDF = HDF_T()
  hdferr = HDF%createFile( fname )

  groupname = 'EMCliffordTorus'
  hdferr = HDF%createGroup( groupname )
end if

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

! first the standard (X,Z_Y) projection
call Message%printMessage('  - projecting quaternions onto Clifford Torus (X,Z_Y)')
do i=1,cnt
  q = self%projectqtoCT_( FZtmp%qu, 'XZ_Y' ) 
  qu(1:4,i) = q%q_copyd()
  qu(5,i) = FZtmp%weight
  FZtmp => FZtmp%next 
end do 
call Message%printMessage('  - projecting Clifford Torus onto Square Torus (X,Z_Y)')
call self%makeSquareTorus(num, w, cnt, xx, yy, offset, qu, h, h2)!, 'XZ_Y')

! do we need to export the raw data to an HDF5 file ?
if (trim(self%nml%hdffile).ne.'undefined') then 
  dataset = 'SquareTorusXZ_Y'
  hdferr = HDF%writeDatasetDoubleArray( dataset, h2, num+2*w, num+2*w )

  dataset = 'ZonePlateXZ_Y'
  hdferr = HDF%writeDatasetDoubleArray( dataset, h, num+2*w, num+2*w )
end if 

! prepare tiff output
if (trim(self%nml%zpfile).ne.'undefined') then
  call self%generateTIFFfile_(EMsoft, SO, num, w, h, 'ZP', 'XZ_Y')
end if 

if (trim(self%nml%sqtfile).ne.'undefined') then
! logarithmic intensity scaling ?
  if (self%nml%logarithmic.eq.1) then 
    logoffset = (maxval(h2)-minval(h2))*0.2D0
    h2 = log10(h2+logoffset)
  end if 
  call self%generateTIFFfile_(EMsoft, SO, num, w, h2, 'SQT', 'XZ_Y')
end if

! then the (Y,X_Z) projection
FZtmp => FZviz
call Message%printMessage('  - projecting quaternions onto Clifford Torus (Y,X_Z)')
do i=1,cnt
  q = self%projectqtoCT_( FZtmp%qu, 'YX_Z' ) 
  qu(1:4,i) = q%q_copyd()
  FZtmp => FZtmp%next 
end do 
call Message%printMessage('  - projecting Clifford Torus onto Square Torus (Y,X_Z)')
h = 0.D0 
h2 = 0.D0
call self%makeSquareTorus(num, w, cnt, xx, yy, offset, qu, h, h2) ! , 'YX_Z')

! do we need to export the raw data to an HDF5 file ?
if (trim(self%nml%hdffile).ne.'undefined') then 
  dataset = 'SquareTorusYX_Z'
  hdferr = HDF%writeDatasetDoubleArray( dataset, h2, num+2*w, num+2*w )

  dataset = 'ZonePlateYX_Z'
  hdferr = HDF%writeDatasetDoubleArray( dataset, h, num+2*w, num+2*w )
end if 

! prepare tiff output
if (trim(self%nml%zpfile).ne.'undefined') then
  call self%generateTIFFfile_(EMsoft, SO, num, w, h, 'ZP', 'YX_Z')
end if 

if (trim(self%nml%sqtfile).ne.'undefined') then
! logarithmic intensity scaling ?
  if (self%nml%logarithmic.eq.1) then 
    logoffset = (maxval(h2)-minval(h2))*0.2D0
    h2 = log10(h2+logoffset)
  end if 
  call self%generateTIFFfile_(EMsoft, SO, num, w, h2, 'SQT', 'YX_Z')
end if

! and finally the (Z,Y_X) projection
FZtmp => FZviz
call Message%printMessage('  - projecting quaternions onto Clifford Torus (Z,Y_X)')
do i=1,cnt
  q = self%projectqtoCT_( FZtmp%qu, 'ZY_X' ) 
  qu(1:4,i) = q%q_copyd()
  FZtmp => FZtmp%next 
end do 
call Message%printMessage('  - projecting Clifford Torus onto Square Torus (Z,Y_X)')
h = 0.D0 
h2 = 0.D0
call self%makeSquareTorus(num, w, cnt, xx, yy, offset, qu, h, h2) ! , 'ZY_X')

! do we need to export the raw data to an HDF5 file ?
if (trim(self%nml%hdffile).ne.'undefined') then 
  dataset = 'SquareTorusZY_X'
  hdferr = HDF%writeDatasetDoubleArray( dataset, h2, num+2*w, num+2*w )

  dataset = 'ZonePlateZY_X'
  hdferr = HDF%writeDatasetDoubleArray( dataset, h, num+2*w, num+2*w )
end if 

! prepare tiff output
if (trim(self%nml%zpfile).ne.'undefined') then
  call self%generateTIFFfile_(EMsoft, SO, num, w, h, 'ZP', 'ZY_X')
end if 

if (trim(self%nml%sqtfile).ne.'undefined') then
! logarithmic intensity scaling ?
  if (self%nml%logarithmic.eq.1) then 
    logoffset = (maxval(h2)-minval(h2))*0.2D0
    h2 = log10(h2+logoffset)
  end if 
  call self%generateTIFFfile_(EMsoft, SO, num, w, h2, 'SQT', 'ZY_X')
end if

if (trim(self%nml%hdffile).ne.'undefined') then 
  call HDF%pop( .TRUE. )
  call closeFortranHDFInterface()
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
use mod_rotations
use mod_quaternions
use mod_io

IMPLICIT NONE 

class(CliffordTorus_T), INTENT(INOUT)   :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 

type(so3_T)                             :: SO
type(IO_T)                              :: Message 
type(QuaternionArray_T)                 :: Pm, dummy 
type(Quaternion_T)                      :: qm, qus 
type(q_T)                               :: qq

character(fnlen)                        :: oname 
integer(kind=irg)                       :: i, k, num, FZcnt, oldFZcnt, io_int(1)
real(kind=dbl)                          :: xx(4), qqd(4)
type(FZpointd),pointer                  :: FZtail, FZtmp, FZhead
logical                                 :: weights

! initialize the SO3 class
SO = so3_T( self%nml%pgnum, zerolist='FZ' )

! read all the orientations from the anglefile
oname = EMsoft%generateFilePath('EMdatapathname',self%nml%anglefile)
call Message%printMessage(' Reading orientations from file '//trim(oname))
call SO%getOrientationsfromFile( oname )

! are we using weighted orientations ?  This could happen with .wxt files that are derived
! from programs like POPLA that extract orientations from an ODG based on pole figures
weights = SO%getuseweights()

! we have the list, so what do we need to do with it before computing the Clifford Torus representation ?

! do we need to generate the -q orientations as well ?
! it is up to the user to make sure that they are not already in the input file !!!
if (self%nml%shownegativeq0.eq.1) then 
  call Message%printMessage(' generating -q orientations ')
! get a pointer to the end of the current linked list
  FZtmp => SO%getListHead('FZ')
  FZtail => SO%getListHead('FZ')
  FZcnt = SO%getListCount('FZ')
  do i=1,FZcnt
    FZtail => FZtail%next
  end do
! loop over all the orientations in the list and generate -q for each one; append it to the end of the list
  do i=1,FZcnt
    qq = FZtmp%qu
    qqd = -qq%q_copyd()
    qq = q_T( qdinp = qqd )
    allocate(FZtail%next)
    FZtail%qu = qq
    FZtail%weight = FZtmp%weight
    FZtail => FZtail%next
    nullify(FZtail%next)
    FZtmp => FZtmp%next
  end do
  FZcnt = 2*FZcnt 
  call SO%setFZcnt(FZcnt, 'FZ')
    io_int(1) = FZcnt
    call Message%WriteValue(' FZcnt after adding -q quaternions = ', io_int, 1)
end if

! do we need to reduce to the RFZ ? 
if (self%nml%reducetoRFZ.eq.1) then 
! get the symmetry operator quaternions for the point group
  call dummy%QSym_Init(self%nml%pgnum, Pm)
  num = Pm%getQnumber()
  io_int(1) = num
  call Message%WriteValue(' Number of symmetry operators ', io_int, 1)
! note that the following step does not propagate any orientation weights... that would require a 
! major rewrite of a lot of the code in mod_so3.f90 ...
  call SO%ReducelisttoRFZ(Pm)
else
! do we need to generate all the equivalent orientations instead ?
  if (self%nml%symmetrize.eq.1) then 
    call Message%printMessage(' applying crystal symmetry to the orientation set')

  ! get a pointer to the end of the current linked list
    FZhead => SO%getListHead('FZ')
    FZtail => SO%getListHead('FZ')
    FZcnt = SO%getListCount('FZ')
    oldFZcnt = FZcnt
    do i=1,FZcnt
      FZtail => FZtail%next
    end do

  ! get the symmetry operator quaternions for the point group
    call dummy%QSym_Init(self%nml%pgnum, Pm)
    num = Pm%getQnumber()
    io_int(1) = num
    call Message%WriteValue(' Number of symmetry operators ', io_int, 1)

  ! loop over all current orientations and generate the equivalent ones
  ! keep in mind that the identity operator is always the first one in the list so we skip it
    do k=2,num 
      FZtmp => FZhead
      qm = Pm%getQuatfromArray(k)
      do i=1,oldFZcnt
        xx = FZtmp%qu%q_copyd() 
        qus = qm * Quaternion_T( qd = xx )
        qq = q_T( qdinp = qus%get_quatd() )
        allocate(FZtail%next)
        FZtail%qu = qq
        FZtail%weight = FZtmp%weight
        FZtail => FZtail%next
        nullify(FZtail%next)       
        FZtmp => FZtmp%next
      end do
    end do
    FZcnt = num*FZcnt 
    call SO%setFZcnt(FZcnt, 'FZ')
    io_int(1) = FZcnt
    call Message%WriteValue(' FZcnt after symmetrization = ', io_int, 1)
  end if 
end if 

! generate the output images using the Clifford Torus sprojection
call Message%printMessage(' Generating square torus/zone plate representation(s)')
call self%createZonePlate_(EMsoft, SO, listmode='FZ' )

! compute Riesz energies ?
if (self%nml%doRiesz.eq.1) then 
  call self%computeRieszEnergy_(SO, listmode='FZ')
end if

end subroutine CliffordTorus_

end module mod_CliffordTorus