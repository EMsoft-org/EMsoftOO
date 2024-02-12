! ###################################################################
! Copyright (c) 2013-2024, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_CBED
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/01/24
  !!
  !! class definition for the EMCBED program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMCBED program
type, public :: CBEDNameListType
  integer(kind=irg)   :: maxHOLZ
  integer(kind=irg)   :: npix
  integer(kind=irg)   :: numthick
  integer(kind=irg)   :: nthreads
  integer(kind=irg)   :: k(3)
  integer(kind=irg)   :: fn(3)
  real(kind=sgl)      :: voltage
  real(kind=sgl)      :: camlen
  real(kind=sgl)      :: klaue(2)
  real(kind=sgl)      :: dmin
  real(kind=sgl)      :: convergence
  real(kind=sgl)      :: startthick
  real(kind=sgl)      :: thickinc
  character(fnlen)    :: xtalname
  character(fnlen)    :: outname
end type CBEDNameListType

! class definition
type, public :: CBED_T
private 
  character(fnlen)       :: nmldeffile = 'EMCBED.nml'
  type(CBEDNameListType)  :: nml 

contains
private 
  procedure, pass(self) :: setmaxHOLZ_
  procedure, pass(self) :: getmaxHOLZ_
  procedure, pass(self) :: setnpix_
  procedure, pass(self) :: getnpix_
  procedure, pass(self) :: setnumthick_
  procedure, pass(self) :: getnumthick_
  procedure, pass(self) :: setnthreads_
  procedure, pass(self) :: getnthreads_
  procedure, pass(self) :: setk_
  procedure, pass(self) :: getk_
  procedure, pass(self) :: setfn_
  procedure, pass(self) :: getfn_
  procedure, pass(self) :: setvoltage_
  procedure, pass(self) :: getvoltage_
  procedure, pass(self) :: setcamlen_
  procedure, pass(self) :: getcamlen_
  procedure, pass(self) :: setklaue_
  procedure, pass(self) :: getklaue_
  procedure, pass(self) :: setdmin_
  procedure, pass(self) :: getdmin_
  procedure, pass(self) :: setconvergence_
  procedure, pass(self) :: getconvergence_
  procedure, pass(self) :: setstartthick_
  procedure, pass(self) :: getstartthick_
  procedure, pass(self) :: setthickinc_
  procedure, pass(self) :: getthickinc_
  procedure, pass(self) :: setxtalname_
  procedure, pass(self) :: getxtalname_
  procedure, pass(self) :: setoutname_
  procedure, pass(self) :: getoutname_
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: CBED_

  generic, public :: setmaxHOLZ => setmaxHOLZ_
  generic, public :: getmaxHOLZ => getmaxHOLZ_
  generic, public :: setnpix => setnpix_
  generic, public :: getnpix => getnpix_
  generic, public :: setnumthick => setnumthick_
  generic, public :: getnumthick => getnumthick_
  generic, public :: setnthreads => setnthreads_
  generic, public :: getnthreads => getnthreads_
  generic, public :: setk => setk_
  generic, public :: getk => getk_
  generic, public :: setfn => setfn_
  generic, public :: getfn => getfn_
  generic, public :: setvoltage => setvoltage_
  generic, public :: getvoltage => getvoltage_
  generic, public :: setcamlen => setcamlen_
  generic, public :: getcamlen => getcamlen_
  generic, public :: setklaue => setklaue_
  generic, public :: getklaue => getklaue_
  generic, public :: setdmin => setdmin_
  generic, public :: getdmin => getdmin_
  generic, public :: setconvergence => setconvergence_
  generic, public :: getconvergence => getconvergence_
  generic, public :: setstartthick => setstartthick_
  generic, public :: getstartthick => getstartthick_
  generic, public :: setthickinc => setthickinc_
  generic, public :: getthickinc => getthickinc_
  generic, public :: setxtalname => setxtalname_
  generic, public :: getxtalname => getxtalname_
  generic, public :: setoutname => setoutname_
  generic, public :: getoutname => getoutname_
  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: CBED => CBED_

end type CBED_T

! the constructor routine for this class 
interface CBED_T
  module procedure CBED_constructor
end interface CBED_T

contains

!--------------------------------------------------------------------------
type(CBED_T) function CBED_constructor( nmlfile ) result(CBED)
!! author: MDG 
!! version: 1.0 
!! date: 02/01/24
!!
!! constructor for the CBED_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call CBED%readNameList(nmlfile)

end function CBED_constructor

!--------------------------------------------------------------------------
subroutine CBED_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 02/01/24
!!
!! destructor for the CBED_T Class
 
IMPLICIT NONE

type(CBED_T), INTENT(INOUT)  :: self 

call reportDestructor('CBED_T')

end subroutine CBED_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/01/24
!!
!! read the namelist from an nml file for the CBED_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(CBED_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg)       :: k(3)
integer(kind=irg)       :: fn(3)
integer(kind=irg)       :: maxHOLZ
integer(kind=irg)       :: numthick
integer(kind=irg)       :: npix
integer(kind=irg)       :: nthreads
real(kind=sgl)          :: voltage
real(kind=sgl)          :: camlen
real(kind=sgl)          :: klaue(2)
real(kind=sgl)          :: dmin
real(kind=sgl)          :: convergence
real(kind=sgl)          :: startthick
real(kind=sgl)          :: thickinc
character(fnlen)        :: xtalname
character(fnlen)        :: outname

namelist /CBEDlist/ xtalname, voltage, k, fn, dmin, convergence, klaue, camlen, &
                    nthreads, startthick, thickinc, numthick, outname, npix, maxHOLZ

k = (/ 0, 0, 1 /)               ! beam direction [direction indices]
fn = (/ 0, 0, 1 /)              ! foil normal [direction indices]
maxHOLZ = 2                     ! maximum HOLZ layer index to be used for the output file; note that his number
                                ! does not affect the actual computations; it only determines which reflection 
                                ! families will end up in the output file
klaue = (/ 0.0, 0.0 /)          ! Laue center coordinates
numthick = 10                   ! number of increments
npix = 256                      ! output arrays will have size npix x npix
nthreads = 1                    ! number of computational threads
voltage = 200.0                 ! acceleration voltage [kV]
camlen = 1000.0                 ! camera length [mm]
dmin = 0.025                    ! smallest d-spacing to include in dynamical matrix [nm]
convergence = 25.0              ! beam convergence angle [mrad]
startthick = 10.0               ! starting thickness [nm]
thickinc = 10.0                 ! thickness increment
xtalname = 'undefined'          ! initial value to check that the keyword is present in the nml file
outname = 'undefined'           ! output filename

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=CBEDlist)
close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(outname).eq.'undefined') then
  call Message%printError('readNameList:',' output file name is undefined in '//nmlfile)
 end if

 if (trim(xtalname).eq.'undefined') then
  call Message%printError('readNameList:',' xtalname is undefined in '//nmlfile)
 end if
end if

self%nml%k = k
self%nml%fn = fn
self%nml%maxHOLZ = maxHOLZ
self%nml%numthick = numthick
self%nml%npix = npix
self%nml%nthreads = nthreads
self%nml%klaue = klaue
self%nml%voltage = voltage
self%nml%camlen = camlen
self%nml%dmin = dmin
self%nml%convergence = convergence
self%nml%startthick = startthick
self%nml%thickinc = thickinc
self%nml%xtalname = xtalname
self%nml%outname = outname

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/01/24
!!
!! pass the namelist for the CBED_T Class to the calling program

IMPLICIT NONE 

class(CBED_T), INTENT(INOUT)          :: self
type(CBEDNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/01/24
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)        :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 4, n_real = 6
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( enl => self%nml )

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! write all the single integers
io_int = (/ enl%maxHOLZ, enl%numthick, enl%npix, enl%nthreads /)
intlist(1) = 'maxHOLZ'
intlist(2) = 'numthick'
intlist(3) = 'npix'
intlist(4) = 'nthreads'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ enl%voltage, enl%camlen, enl%dmin, enl%convergence, enl%startthick, enl%thickinc /)
reallist(1) = 'voltage'
reallist(2) = 'camlen'
reallist(3) = 'dmin'
reallist(4) = 'convergence'
reallist(5) = 'startthick'
reallist(6) = 'thickinc'
call HDF%writeNMLreals(io_real, reallist, n_real)

! 3-vectors
dataset = 'k'
hdferr = HDF%writeDatasetIntegerArray(dataset, enl%k, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create k dataset', hdferr)

dataset = 'fn'
hdferr = HDF%writeDatasetIntegerArray(dataset, enl%fn, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create fn dataset', hdferr)

dataset = 'klaue'
hdferr = HDF%writeDatasetFloatArray(dataset, enl%klaue, 2)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create klaue dataset', hdferr)

! write all the strings
dataset = SC_xtalname
line2(1) = trim(enl%xtalname)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create xtalname dataset', hdferr)

dataset = SC_outname
line2(1) = trim(enl%outname)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create outname dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine setmaxHOLZ_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setmaxHOLZ_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! set maxHOLZ in the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%maxHOLZ = inp

end subroutine setmaxHOLZ_

!--------------------------------------------------------------------------
function getmaxHOLZ_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getmaxHOLZ_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! get maxHOLZ from the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%maxHOLZ

end function getmaxHOLZ_

!--------------------------------------------------------------------------
subroutine setnpix_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnpix_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! set npix in the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%npix = inp

end subroutine setnpix_

!--------------------------------------------------------------------------
function getnpix_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnpix_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! get npix from the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%npix

end function getnpix_

!--------------------------------------------------------------------------
subroutine setnumthick_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnumthick_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! set numthick in the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%numthick = inp

end subroutine setnumthick_

!--------------------------------------------------------------------------
function getnumthick_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnumthick_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! get numthick from the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%numthick

end function getnumthick_

!--------------------------------------------------------------------------
subroutine setnthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnthreads_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! set nthreads in the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%nthreads = inp

end subroutine setnthreads_

!--------------------------------------------------------------------------
function getnthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnthreads_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! get nthreads from the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%nthreads

end function getnthreads_

!--------------------------------------------------------------------------
subroutine setk_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setk_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! set k in the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp(3)

self%nml%k = inp

end subroutine setk_

!--------------------------------------------------------------------------
function getk_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getk_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! get k from the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out(3)

out = self%nml%k

end function getk_

!--------------------------------------------------------------------------
subroutine setfn_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setfn_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! set fn in the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp(3)

self%nml%fn = inp

end subroutine setfn_

!--------------------------------------------------------------------------
function getfn_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getfn_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! get fn from the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out(3)

out = self%nml%fn

end function getfn_

!--------------------------------------------------------------------------
subroutine setvoltage_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setvoltage_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! set voltage in the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%voltage = inp

end subroutine setvoltage_

!--------------------------------------------------------------------------
function getvoltage_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getvoltage_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! get voltage from the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%voltage

end function getvoltage_

!--------------------------------------------------------------------------
subroutine setcamlen_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setcamlen_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! set camlen in the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%camlen = inp

end subroutine setcamlen_

!--------------------------------------------------------------------------
function getcamlen_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getcamlen_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! get camlen from the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%camlen

end function getcamlen_

!--------------------------------------------------------------------------
subroutine setklaue_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setklaue_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! set klaue in the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp(2)

self%nml%klaue = inp

end subroutine setklaue_

!--------------------------------------------------------------------------
function getklaue_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getklaue_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! get klaue from the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out(2)

out = self%nml%klaue

end function getklaue_

!--------------------------------------------------------------------------
subroutine setdmin_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdmin_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! set dmin in the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%dmin = inp

end subroutine setdmin_

!--------------------------------------------------------------------------
function getdmin_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdmin_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! get dmin from the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%dmin

end function getdmin_

!--------------------------------------------------------------------------
subroutine setconvergence_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setconvergence_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! set convergence in the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%convergence = inp

end subroutine setconvergence_

!--------------------------------------------------------------------------
function getconvergence_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getconvergence_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! get convergence from the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%convergence

end function getconvergence_

!--------------------------------------------------------------------------
subroutine setstartthick_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setstartthick_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! set startthick in the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%startthick = inp

end subroutine setstartthick_

!--------------------------------------------------------------------------
function getstartthick_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getstartthick_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! get startthick from the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%startthick

end function getstartthick_

!--------------------------------------------------------------------------
subroutine setthickinc_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setthickinc_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! set thickinc in the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%thickinc = inp

end subroutine setthickinc_

!--------------------------------------------------------------------------
function getthickinc_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getthickinc_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! get thickinc from the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%thickinc

end function getthickinc_

!--------------------------------------------------------------------------
subroutine setxtalname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setxtalname_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! set xtalname in the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%xtalname = trim(inp)

end subroutine setxtalname_

!--------------------------------------------------------------------------
function getxtalname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getxtalname_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! get xtalname from the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%xtalname)

end function getxtalname_

!--------------------------------------------------------------------------
subroutine setoutname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setoutname_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! set outname in the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%outname = trim(inp)

end subroutine setoutname_

!--------------------------------------------------------------------------
function getoutname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getoutname_
!! author: MDG
!! version: 1.0
!! date: 02/01/24
!!
!! get outname from the CBED_T class

IMPLICIT NONE

class(CBED_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%outname)

end function getoutname_

!--------------------------------------------------------------------------
subroutine CBED_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: CBED_
!! author: MDG 
!! version: 1.0 
!! date: 02/01/24
!!
!! perform the computations

use mod_EMsoft
use mod_HDFnames
use mod_crystallography
use mod_diffraction
use mod_symmetry
use mod_symmetry2D
use mod_kvectors
use mod_gvectors
use mod_Lambert
use mod_math
use mod_io
use mod_initializers
use mod_HDFsupport
use HDF5
use mod_timing
use omp_lib
use mod_memory
use mod_HOLZ
use stringconstants
use ISO_C_BINDING
use mod_image

use, intrinsic :: iso_fortran_env


IMPLICIT NONE 

class(CBED_T), INTENT(INOUT)          :: self
type(EMsoft_T), INTENT(INOUT)         :: EMsoft
character(fnlen), INTENT(INOUT)       :: progname 
type(HDFnames_T), INTENT(INOUT)       :: HDFnames

real(kind=sgl)                        :: ktmax, io_real(3), galen, bragg, RR, gg(3), thetac, &
                                         sc, scmax, PX, frac, klaue(2), pxy(2), mi, ma
real(kind=dbl)                        :: WL, kinp(3), kv(3) 
integer(kind=irg)                     :: ijmax,ga(3),gb(3), cnt, skip, istat, dgn, badpoints, DynNbeamsLinked, maxthreads, &
                                         newcount,count_rate,count_max, io_int(6), ii, i, j, isym, ir, pgnum, ierr, DynNbeams, &
                                         npx, npy, numt, numk, npix, ik, ip, jp, iequiv(2,12), nequiv, it, nref, hdferr
character(3)                          :: method
type(kvectorlist),pointer             :: khead, ktmp
type(reflisttype),pointer             :: reflist, rltmpa, rltmpb, gtmp

type(Cell_T)                          :: cell
type(SpaceGroup_T)                    :: SG
type(HDF_T)                           :: HDF
type(Diffraction_T)                   :: Diff
type(Timing_T)                        :: timer
type(memory_T)                        :: mem
type(HOLZ_T)                          :: HOLZ
type(IO_T)                            :: Message
type(HOLZentries)                     :: HOLZdata 
type(gvectors_T)                      :: gvec
type(kvectors_T)                      :: kvec
type(symdata2D)                       :: TDPG

type(reflisttype),pointer             :: listrootw

character(fnlen)                      :: dataset, instring, outname, groupname, datagroupname, TIFF_filename
real(kind=sgl),parameter              :: xoff(0:5)=(/0.0,3.3125,0.0,3.3125,0.0,3.3125/),yoff(0:5)=(/6.0,6.0,3.0,3.0,0.0,0.0/)
real(kind=sgl),allocatable            :: disk(:,:,:), thick(:), slice(:,:)
integer(kind=irg),allocatable         :: diskoffset(:,:)
real(kind=sgl),allocatable            :: inten(:,:)
complex(kind=dbl),allocatable         :: DynMat(:,:)
logical                               :: usesym=.TRUE., bp, verbose, f_exists, g_exists, overwrite=.TRUE.
character(11)                         :: dstr
character(15)                         :: tstrb
character(15)                         :: tstre
character(2)                          :: str

! declare variables for use in object oriented image module
integer                               :: iostat
character(len=128)                    :: iomsg
logical                               :: isInteger
type(image_t)                         :: im
integer(int8)                         :: i8 (3,4)
integer(int8), allocatable            :: TIFF_image(:,:)

call openFortranHDFInterface()

associate( enl => self%nml, Dyn => Diff%Dyn )

! max number of OpenMP threads on this platform
maxthreads = omp_get_max_threads()

! memory handling class
mem = memory_T()

! start the timer
timer = Timing_T()
tstrb = timer%getTimeString()
dstr = timer%getDateString()

! initialize the HDF class
HDF = HDF_T() 

! crystallography section
verbose = .TRUE.

call cell%setFileName(enl%xtalname)
call Diff%setrlpmethod('WK')
Diff%Dyn%FN = enl%fn 
call Diff%setV(dble(enl%voltage))
! initialize the HOLZ geometry type
HOLZ = HOLZ_T()
HOLZ%maxHOLZ = self%getmaxHOLZ()
! initialize the gvectors list
gvec = gvectors_T()

! get the crystal structure data and initialize the LUTs under the assumption of 
! a single zone axis and consequently a set of HOLZ layers
call Initialize_Cell_HOLZ(cell, Diff, SG, HOLZ, HOLZdata, gvec, TDPG, EMsoft, & 
                          enl%dmin, enl%k, ga, gb, verbose, useHDF=HDF)

! force dynamical matrix routine to read new Bethe parameters from file
call Diff%SetBetheParameters(EMsoft,silent=.TRUE.)

! construct the list of all possible reflections
thetac = enl%convergence/1000.0
call Diff%setDynNbeamsLinked( gvec%get_nref() )
DynNbeamsLinked = gvec%get_nref()

! enter range of incident beam directions
bragg = Diff%CalcDiffAngle(cell, ga)*0.5
  
! convert to ktmax along ga
ktmax = 0.5*thetac/bragg

! compute number of pixels along diameter of central disk for given camera length
RR = 300.0/25.4   ! dots per millimeter for 300 dots per inch; legacy code from when the output was in PostScript
npx = int(RR*enl%camlen*thetac)
npy = npx
io_int(1) = 2.0*npx
call Message%WriteValue('Number of image pixels along diameter of central disk = ', io_int, 1, "(I4)")
call Message%printMessage(' ', "(A/)")
  
! get number of thicknesses for which to compute the CBED pattern
numt = enl%numthick
call mem%alloc(thick,(/ numt /),'thick', 0.0)
thick = enl%startthick + enl%thickinc* (/ (float(i),i=0,numt-1) /)

! if the Laue center is at the origin, then we can use symmetry groups to 
! speed up the simulation; otherwise we have to cover each incident wave
! vector separately.
if (maxval(abs(enl%klaue)).eq.0.0) then
  usesym=.TRUE.
else
  usesym=.FALSE.
end if

! determine all independent incident beam directions (use a linked list starting at khead)
! isym = 1
! determine the point group number
 j=0
 do i=1,32
  if (SGPG(i).le.SG%getSpaceGroupNumber()) j=i
 end do

dgn = SG%GetPatternSymmetry(enl%k,j,.TRUE.)
pgnum = j
isym = WPPG(dgn)
ijmax = float(npx)**2   ! truncation value for beam directions
kvec = kvectors_T()
kinp = dble(enl%k)
call cell%TransSpace( kinp, kv, 'd', 'r') 
call cell%NormVec( kv, 'r' )
call kvec%set_kinp( kv / Diff%getWaveLength() )
write (*,*) 'isym value = ', isym
call kvec%set_isym( isym )
call kvec%set_ktmax( dble(ktmax) )
call kvec%CalckvectorsSymmetry(cell,Diff,TDPG,dble(ga),npx,npy,ijmax,enl%klaue,.TRUE.)
numk = kvec%get_numk()
  
! allocate the disk variable which will hold the entire computed pattern
call mem%alloc(disk, (/ numt,enl%npix,enl%npix /), 'disk', 0.0)

WL = Diff%getWaveLength() 
sc = WL * enl%camlen * RR
PX = enl%npix/2
scmax = PX + npx
npix = enl%npix

! allocate the offset array
! to get this array, we need to do a mock initialization of the dynamical matrix in zone axis orientation
! point to the first beam direction
ktmp => kvec%get_ListHead()

call gvec%GetDynMatHOLZ(cell, Diff, HOLZ, HOLZdata, EMsoft, 'BLOCHBETHE', dble(ktmp%k), dble(ktmp%kt), .FALSE. )
DynNbeamsLinked = Diff%getDynNbeamsLinked()
call mem%alloc(diskoffset, (/ DynNbeamsLinked,3 /), 'diskoffset', 0)

! project every g-vector onto gx and gy to get the components
! and keep only the ones that will fall on the viewing region
rltmpa => gvec%get_ListHead()
rltmpa => rltmpa%next
do i=1,gvec%get_nref()
  gg(1:3)=rltmpa%hkl

  pxy =  sc * HOLZ%GetHOLZcoordinates(cell, HOLZdata, gg, (/ 0.0, 0.0, 0.0 /), sngl(WL))

  if ((abs(pxy(1)).lt.scmax).and.(abs(pxy(2)).lt.scmax)) then
    diskoffset(i,1) = 1
    diskoffset(i,2) = nint(pxy(1))
    diskoffset(i,3) = nint(pxy(2))
  else
    diskoffset(i,1) = 0
  end if
  rltmpa => rltmpa%next
end do


frac = 0.05
badpoints = 0

call Message%printMessage(' ', "(A/)")
io_int(1)=numk
call Message%WriteValue(' Starting computation for # beam directions = ', io_int, 1, "(I6)")

ktmp => kvec%get_ListHead()

do ik = 1, numk
  call gvec%GetDynMatHOLZ(cell, Diff, HOLZ, HOLZdata, EMsoft, 'BLOCHBETHE', dble(ktmp%k), dble(ktmp%kt), .FALSE. )

! allocate the intensity array to include both strong beams and weak beams (in that order)
  DynNbeams = Diff%getDynNbeams()
  DynNbeamsLinked = Diff%getDynNbeamsLinked()
  call mem%alloc(inten, (/ numt,DynNbeams+Diff%BetheParameters%nnw /), 'inten', 0.0)

! write (*,*) ik, Diff%BetheParameters%nns, Diff%BetheParameters%nnw, Diff%BetheParameters%sgcutoff

! solve the dynamical eigenvalue equation and return the intensities of ALL reflections,
! both strong and weak; the weak intensities should also be plotted at the correct locations....
! this uses a new version of the CalcBWint routine that implements both strong and weak beam intensities.
  call gvec%CalcBWint(cell,Diff,ktmp,Diff%BetheParameters%nns,Diff%BetheParameters%nnw,numt,thick,inten)

! if the Bethe truncation parameters are not set correctly, then there may 
! be bad intensity points for some beams.  If this happens, we set the 
! corresponding intensities to zero and raise a flag that will cause the
! program to print a warning message at the end.  We also count the number of 
! bad points.  
  bp = .FALSE. 
  if (minval(inten).lt.0.0) then
    where (inten.lt.0.0) inten = 0.0
    bp = .TRUE.
  end if
  if (maxval(inten).gt.1.0) then
    where (inten.gt.1.0) inten = 0.0
    bp = .TRUE.
  end if
  if (bp) badpoints = badpoints+1

! we combine the reflistindex and weakreflistindex into one to make the CBED drawing.
  Diff%BetheParameters%reflistindex = Diff%BetheParameters%reflistindex + Diff%BetheParameters%weakreflistindex

! and copy the strong intensities in the correct locations
   do i=1,DynNbeamsLinked
    if ( (diskoffset(i,1).eq.1).and.(Diff%BetheParameters%reflistindex(i).ne.0)) then
     if (usesym) then ! use 2D point group symmetry
      ip = diskoffset(i,2) - ktmp%i
      jp = diskoffset(i,3) + ktmp%j
      call Apply2DPGSymmetry(TDPG,ip,jp,isym,iequiv,nequiv)
      iequiv = iequiv+PX

      do ii=1,nequiv
! is this point inside the viewing square ?
        ip = iequiv(1,ii)
        jp = iequiv(2,ii)
        if (((ip.ge.1).and.(ip.le.npix)).and.((jp.ge.1).and.(jp.le.npix))) then
         if ( (Diff%BetheParameters%reflistindex(i).eq.1) ) then
          disk(1:numt,ip,jp) = disk(1:numt,ip,jp) + inten(1:numt,Diff%BetheParameters%reflistindex(i))
         else
! this is a placeholder; it is not technically correct to do this, but the result looks quite reasonable
! it works fine when there is no overlap between diffraction disks, but when there is, then the result will
! be incorrect.  This should be noted in the manual. 
          do it=1,numt
           disk(it,ip,jp) = maxval( (/ inten(it,Diff%BetheParameters%reflistindex(i)), disk(it,ip,jp) /) )
          end do
         end if
        end if
      end do
     else  ! do not use symmetry
       ip = PX + diskoffset(i,2) - ktmp%i
       jp = PX + diskoffset(i,3) + ktmp%j
       if (((ip.ge.1).and.(ip.le.npix)).and.((jp.ge.1).and.(jp.le.npix))) then
          disk(1:numt,ip,jp) = disk(1:numt,ip,jp) + inten(1:numt,Diff%BetheParameters%reflistindex(i))
        end if
     end if
    end if
   end do

  call mem%dealloc(inten, 'inten')

! select next beam direction
   if (ik.ne.numk) ktmp => ktmp%next

! update computation progress
   if (float(ik)/float(numk) .gt. frac) then
    io_int(1) = nint(100.0*frac) 
    call Message%WriteValue('       ', io_int, 1, "(1x,I3,' percent completed')") 
    frac = frac + 0.05
   end if  

end do 

TIFF_filename = 'CBEDpattern.tiff'
allocate(TIFF_image(enl%npix,enl%npix))

ma = maxval(log10(disk(numt,:,:)+1.0E-4))
mi = minval(log10(disk(numt,:,:)+1.0E-4))

TIFF_image = int(255 * (log10(disk(numt,1:enl%npix,1:enl%npix)+1.0E-4)-mi)/(ma-mi))

! set up the image_t structure
im = image_t(TIFF_image)
if(im%empty()) call Message%printMessage("EMCBED","failed to convert array to image")

! create the file
call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage("failed to write image to file : "//iomsg)
else
  call Message%printMessage('CBED pattern written to '//trim(TIFF_filename))
end if


! stop the timer 
tstre = timer%getTimeString()

! HDF output 
outname = trim(EMsoft%generateFilePath('EMdatapathname',enl%outname))
hdferr =  HDF%createFile(outname)

! write the EMheader to the file
datagroupname = trim(HDFnames%get_ProgramData())
call HDF%writeEMheader(EMsoft,dstr, tstrb, tstre, progname, datagroupname)

! open or create a namelist group to write all the namelist files into
groupname = SC_NMLfiles
  hdferr = HDF%createGroup(HDFnames%get_NMLfiles())

! read the text file and write the array to the file
dataset = trim(HDFnames%get_NMLfilename())
  hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%nmldeffile)

! leave this group
call HDF%pop()

! create a namelist group to write all the namelist files into
hdferr = HDF%createGroup(HDFnames%get_NMLparameters())

call self%writeHDFNameList(HDF, HDFnames)
call Diff%writeBetheparameterNameList(HDF)

! leave this group
call HDF%pop()

! then the remainder of the data in a EMData group
hdferr = HDF%createGroup(HDFnames%get_EMData())
hdferr = HDF%createGroup(HDFnames%get_ProgramData())

dataset = 'CBEDpatterns'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetFloatArray(dataset, disk, numt, npix, npix, overwrite)
  else
    hdferr = HDF%writeDatasetFloatArray(dataset, disk, numt, npix, npix)
  end if

call HDF%popall()


! and clean up
end associate

call mem%dealloc(thick, 'thick')
call mem%dealloc(disk, 'disk')
call mem%dealloc(diskoffset, 'diskofset')

call closeFortranHDFInterface()

end subroutine CBED_



end module mod_CBED