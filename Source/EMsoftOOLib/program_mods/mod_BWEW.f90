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

module mod_BWEW
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/19/24
  !!
  !! class definition for the EMBWEW program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMBWEW program
type, public :: BWEWNameListType
  integer(kind=irg)   :: k(3)
  integer(kind=irg)   :: npix
  integer(kind=irg)   :: numthick
  real(kind=sgl)      :: maxsamint
  real(kind=sgl)      :: startthick
  real(kind=sgl)      :: thickinc
  real(kind=sgl)      :: voltage
  real(kind=sgl)      :: klauec(2)
  character(fnlen)    :: xtalname 
  character(fnlen)    :: outname 
  character(fnlen)    :: tiffprefix
end type BWEWNameListType

! class definition
type, public :: BWEW_T
private 
  character(fnlen)       :: nmldeffile = 'EMBWEW.nml'
  type(BWEWNameListType)  :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: BWEW_
  procedure, pass(self) :: setk_
  procedure, pass(self) :: getk_  
  procedure, pass(self) :: setklauec_
  procedure, pass(self) :: getklauec_
  procedure, pass(self) :: setnpix_
  procedure, pass(self) :: getnpix_
  procedure, pass(self) :: setnumthick_
  procedure, pass(self) :: getnumthick_
  procedure, pass(self) :: setmaxsamint_
  procedure, pass(self) :: getmaxsamint_
  procedure, pass(self) :: setstartthick_
  procedure, pass(self) :: getstartthick_
  procedure, pass(self) :: setthickinc_
  procedure, pass(self) :: getthickinc_
  procedure, pass(self) :: setvoltage_
  procedure, pass(self) :: getvoltage_
  procedure, pass(self) :: setxtalname_
  procedure, pass(self) :: getxtalname_
  procedure, pass(self) :: setoutname_
  procedure, pass(self) :: getoutname_
  procedure, pass(self) :: settiffprefix_
  procedure, pass(self) :: gettiffprefix_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: BWEW => BWEW_
  generic, public :: setk => setk_
  generic, public :: getk => getk_
  generic, public :: setklauec => setklauec_
  generic, public :: getklauec => getklauec_
  generic, public :: setnpix => setnpix_
  generic, public :: getnpix => getnpix_
  generic, public :: setnumthick => setnumthick_
  generic, public :: getnumthick => getnumthick_
  generic, public :: setmaxsamint => setmaxsamint_
  generic, public :: getmaxsamint => getmaxsamint_
  generic, public :: setstartthick => setstartthick_
  generic, public :: getstartthick => getstartthick_
  generic, public :: setthickinc => setthickinc_
  generic, public :: getthickinc => getthickinc_
  generic, public :: setvoltage => setvoltage_
  generic, public :: getvoltage => getvoltage_
  generic, public :: setxtalname => setxtalname_
  generic, public :: getxtalname => getxtalname_
  generic, public :: setoutname => setoutname_
  generic, public :: getoutname => getoutname_
  generic, public :: settiffprefix => settiffprefix_
  generic, public :: gettiffprefix => gettiffprefix_
end type BWEW_T

! the constructor routine for this class 
interface BWEW_T
  module procedure BWEW_constructor
end interface BWEW_T

contains

!--------------------------------------------------------------------------
type(BWEW_T) function BWEW_constructor( nmlfile ) result(BWEW)
!! author: MDG 
!! version: 1.0 
!! date: 02/19/24
!!
!! constructor for the BWEW_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call BWEW%readNameList(nmlfile)

end function BWEW_constructor

!--------------------------------------------------------------------------
subroutine BWEW_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 02/19/24
!!
!! destructor for the BWEW_T Class
 
IMPLICIT NONE

type(BWEW_T), INTENT(INOUT)  :: self 

call reportDestructor('BWEW_T')

end subroutine BWEW_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/19/24
!!
!! read the namelist from an nml file for the BWEW_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(BWEW_T), INTENT(INOUT)         :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg)                    :: k(3)
integer(kind=irg)                    :: npix
integer(kind=irg)                    :: numthick
real(kind=sgl)                       :: klauec(2)
real(kind=sgl)                       :: maxsamint
real(kind=sgl)                       :: startthick
real(kind=sgl)                       :: thickinc
real(kind=sgl)                       :: voltage
character(fnlen)                     :: xtalname 
character(fnlen)                     :: outname 
character(fnlen)                     :: tiffprefix

namelist /BWEWlist/ xtalname, voltage, k, tiffprefix, startthick, thickinc, numthick, outname, &
                    npix, maxsamint, klauec

k = (/ 0, 0, 1 /)
klauec = (/ 0.0, 0.0 /)
npix = 256 
numthick = 10
maxsamint = 0.025
startthick = 10.0
thickinc = 10.0
voltage = 200.0
xtalname = 'undefined'
outname = 'undefined'
tiffprefix = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=BWEWlist)
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
self%nml%klauec = klauec
self%nml%npix = npix
self%nml%numthick = numthick
self%nml%maxsamint = maxsamint
self%nml%startthick = startthick
self%nml%thickinc = thickinc
self%nml%voltage = voltage
self%nml%xtalname = xtalname
self%nml%outname = outname
self%nml%tiffprefix = tiffprefix

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/19/24
!!
!! pass the namelist for the BWEW_T Class to the calling program

IMPLICIT NONE 

class(BWEW_T), INTENT(INOUT)          :: self
type(BWEWNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/19/24
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)            :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 2, n_real = 4
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( enl => self%nml )

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! write all the single integers
io_int = (/ enl%numthick, enl%npix /)
intlist(1) = 'numthick'
intlist(2) = 'npix'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ enl%voltage, enl%startthick, enl%thickinc, enl%maxsamint /)
reallist(1) = 'voltage'
reallist(2) = 'startthick'
reallist(3) = 'thickinc'
reallist(4) = 'maxsamint'
call HDF%writeNMLreals(io_real, reallist, n_real)

! 3-vector
dataset = 'k'
hdferr = HDF%writeDatasetIntegerArray(dataset, enl%k, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create k dataset', hdferr)

! 2-vector
dataset = 'klauec'
hdferr = HDF%writeDatasetFloatArray(dataset, enl%klauec, 2)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create klauec dataset', hdferr)

! write all the strings
dataset = SC_xtalname
line2(1) = trim(enl%xtalname)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create xtalname dataset', hdferr)

dataset = SC_outname
line2(1) = trim(enl%outname)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create outname dataset', hdferr)

dataset = SC_tiffprefix
line2(1) = trim(enl%tiffprefix)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create tiffprefix dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine setk_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setk_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! set k in the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)    :: inp(3)

self%nml%k = inp

end subroutine setk_

!--------------------------------------------------------------------------
function getk_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getk_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! get k from the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
integer(kind=irg)                :: out(3)

out = self%nml%k

end function getk_

!--------------------------------------------------------------------------
subroutine setnpix_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnpix_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! set npix in the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%npix = inp

end subroutine setnpix_

!--------------------------------------------------------------------------
function getnpix_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnpix_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! get npix from the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%npix

end function getnpix_

!--------------------------------------------------------------------------
subroutine setnumthick_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnumthick_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! set numthick in the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%numthick = inp

end subroutine setnumthick_

!--------------------------------------------------------------------------
function getnumthick_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnumthick_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! get numthick from the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%numthick

end function getnumthick_

!--------------------------------------------------------------------------
subroutine setmaxsamint_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setmaxsamint_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! set maxsamint in the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%maxsamint = inp

end subroutine setmaxsamint_

!--------------------------------------------------------------------------
function getmaxsamint_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getmaxsamint_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! get maxsamint from the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%maxsamint

end function getmaxsamint_

!--------------------------------------------------------------------------
subroutine setstartthick_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setstartthick_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! set startthick in the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%startthick = inp

end subroutine setstartthick_

!--------------------------------------------------------------------------
function getstartthick_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getstartthick_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! get startthick from the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%startthick

end function getstartthick_

!--------------------------------------------------------------------------
subroutine setklauec_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setklauec_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! set klauec in the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp(2)

self%nml%klauec = inp

end subroutine setklauec_

!--------------------------------------------------------------------------
function getklauec_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getklauec_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! get klauec from the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out(2)

out = self%nml%klauec

end function getklauec_

!--------------------------------------------------------------------------
subroutine setthickinc_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setthickinc_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! set thickinc in the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%thickinc = inp

end subroutine setthickinc_

!--------------------------------------------------------------------------
function getthickinc_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getthickinc_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! get thickinc from the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%thickinc

end function getthickinc_

!--------------------------------------------------------------------------
subroutine setvoltage_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setvoltage_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! set voltage in the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%voltage = inp

end subroutine setvoltage_

!--------------------------------------------------------------------------
function getvoltage_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getvoltage_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! get voltage from the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%voltage

end function getvoltage_

!--------------------------------------------------------------------------
subroutine setxtalname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setxtalname_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! set xtalname in the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%xtalname = trim(inp)

end subroutine setxtalname_

!--------------------------------------------------------------------------
function getxtalname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getxtalname_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! get xtalname from the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%xtalname)

end function getxtalname_

!--------------------------------------------------------------------------
subroutine setoutname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setoutname_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! set outname in the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%outname = trim(inp)

end subroutine setoutname_

!--------------------------------------------------------------------------
function getoutname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getoutname_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! get outname from the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%outname)

end function getoutname_

!--------------------------------------------------------------------------
subroutine settiffprefix_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: settiffprefix_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! set tiffprefix in the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%tiffprefix = trim(inp)

end subroutine settiffprefix_

!--------------------------------------------------------------------------
function gettiffprefix_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: gettiffprefix_
!! author: MDG
!! version: 1.0
!! date: 02/19/24
!!
!! get tiffprefix from the BWEW_T class

IMPLICIT NONE

class(BWEW_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%tiffprefix)

end function gettiffprefix_

!--------------------------------------------------------------------------
subroutine BWEW_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: BWEW_
!! author: MDG 
!! version: 1.0 
!! date: 02/19/24
!!
!! perform the computations

use mod_kinds
use mod_global
use mod_EMsoft
use mod_gvectors
use mod_HDFnames
use mod_crystallography
use mod_diffraction
use mod_initializers
use mod_symmetry
use mod_symmetry2D
use mod_kvectors
use mod_math
use mod_io
use mod_HDFsupport 
use mod_memory
use mod_timing
use mod_initializers
use mod_image
use HDF5
use mod_HDFsupport

use, intrinsic :: iso_fortran_env

IMPLICIT NONE 

class(BWEW_T), INTENT(INOUT)      :: self
type(EMsoft_T), INTENT(INOUT)     :: EMsoft
character(fnlen), INTENT(INOUT)   :: progname 
type(HDFnames_T), INTENT(INOUT)   :: HDFnames

type(Cell_T)                      :: Cell 
type(Spacegroup_T)                :: SG 
type(Diffraction_T)               :: Diff
type(DynType)                     :: Dyn
type(IO_T)                        :: Message
type(Memory_T)                    :: mem
type(Timing_T)                    :: timer
type(gnode)                       :: rlp
type(gvectors_T)                  :: gvec
type(kvectors_T)                  :: kvec
type(HDF_T)                       :: HDF
type(reflisttype),pointer         :: reflist, rltmpa, rl, firstw

real(kind=sgl)                    :: laL,kt,z0,thc,thb,hkl(3),ind(3),fn(3),thick,mi, ma, galen, &
                                     dom,glen,c(3),RR,gx(3),gy(3),gg(3),kstar(3),gad(3),gbd(3),frac,gal,gbl, &
                                     gmax,gamax,gbmax,gadl,gbdl, io_real(3), dmin, DynWV(3), lambda, Upz, xgpz
real(kind=dbl)                    :: arg,th,maxsamint,dummy(3),pixpos(3)
integer(kind=irg)                 :: g(3),ira,dpcnt,ppi,ijmax,ga(3),gb(3),k(3),fcnt,ccnt,count, &
                                     newcount,count_rate,count_max,nn,i,j,ig,ih,isym,ir,iorder,nt, &
                                     npix,npiy,istat,numt,ip,jp,ik,numk,mins, io_int(6), dgn, nns, nnw, pgnum
character(1)                      :: ans
character(2)                      :: srza
character(20)                     :: fname
character(11)                     :: dstr
character(15)                     :: tstrb
character(15)                     :: tstre
character(3)                      :: str
character(fnlen)                  :: TIFF_filename

integer(kind=irg),allocatable     :: IPIV(:)
complex(kind=dbl),allocatable     :: CGinv(:,:), Minp(:,:),diag(:),amp(:,:),pw(:), W(:), CG(:,:), alpha(:), DynMat(:,:)
complex(kind=dbl)                 :: czero = cmplx(0.0,0.0,dbl)
complex(kind=sgl),allocatable     :: ew(:,:,:)
real(kind=sgl),allocatable        :: inten(:,:)
real(kind=dbl),allocatable        :: kg(:,:)

! declare variables for use in object oriented image module
integer                           :: iostat
character(len=128)                :: iomsg
logical                           :: isInteger, verbose
type(image_t)                     :: im
integer(int8)                     :: i8 (3,4)
integer(int8), allocatable        :: TIFF_image(:,:)


call openFortranHDFInterface()

associate( enl => self%nml )

! memory handling class
mem = memory_T()

! start the timer
timer = Timing_T()
tstrb = timer%getTimeString()
dstr = timer%getDateString()

! initialize the HDF class
HDF = HDF_T() 

! maximum sampling interval
maxsamint = enl%maxsamint 
numt = enl%numthick
npix = enl%npix

! incident wave vector direction 
k = enl%k 

! crystallography section
verbose = .TRUE.

call cell%setFileName(enl%xtalname)
call Diff%setrlpmethod('WK')
! foil normal is assumed to be parallel to k 
Diff%Dyn%FN = dble(k) 
dmin = 0.005D0   
call Diff%setV(dble(enl%voltage))

! initialize the gvectors list
gvec = gvectors_T()

! get the crystal structure data and initialize the LUTs 
call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, dmin, verbose, useHDF=HDF)

! determine the point group number
j=0
do i=1,32
    if (SGPG(i).le.SG%getSpaceGroupNumber()) j=i
end do
isym = j

! use the new routine to get the whole pattern 2D symmetry group, since that
! is the one that determines the independent beam directions.
dgn = SG%GetPatternSymmetry(enl%k,j,.TRUE.)
pgnum = j
isym = WPPG(dgn)

! force dynamical matrix routine to read new Bethe parameters from file
call Diff%SetBetheParameters(EMsoft)
call Diff%setDynNbeamsLinked( gvec%get_nref() )

lambda = Diff%getWaveLength()

!determine the shortest reciprocal lattice points for this zone
call cell%ShortestG(SG,enl%k,ga,gb,isym)
io_int(1:3)=ga(1:3)
io_int(4:6)=gb(1:3)
call Message%WriteValue(' Reciprocal lattice vectors : ', io_int, 6,"('(',3I3,') and (',3I3,')')")
call Message%printMessage('  (the first lattice vector is horizontal in the CBED pattern)')
call Message%printMessage(' ')
galen = cell%CalcLength(float(ga), 'r')

! normal aborption factor
call Diff%CalcUcg(cell, (/ 0, 0, 0 /))
rlp = Diff%getrlp()
xgpz= rlp%xgp
Upz = rlp%Upmod

! compute real space vectors delimiting the unit cell
call cell%TransSpace(float(ga),gad,'r','d')
call cell%TransSpace(float(gb),gbd,'r','d')
call cell%NormVec(gad,'d')
call cell%Normvec(gbd,'d')
gal = cell%CalcLength(float(ga),'r')
gbl = cell%CalcLength(float(gb),'r')
gad = gad/gal
gbd = gbd/gbl
call Message%printMessage(' Direct space unit cell vectors :')
io_real(1:3) = gad(1:3)
call Message%WriteValue('   x-axis = ', io_real, 3, "(3F10.5)")
io_real(1:3) = gbd(1:3)
call Message%WriteValue('   y-axis = ', io_real, 3, "(3F10.5)")
! get number of pixels along x and y
! if we take a sampling with 0.025 nm per pixel as maximum,
! then the minimum number of sampling points is given by
mins = nint(1.0/gal/maxsamint)
io_int(1) = mins
call Message%WriteValue(' Minimum number of sampling pixels along x-axis :', io_int, 1, "(I4)")
! do not allow less than the minimum
if (npix.lt.mins) npix = mins
npiy = nint(npix * gal/gbl)
io_int(1:2) = (/ npix,npiy /)
call Message%WriteValue(' Number of image pixels set to : ', io_int, 2,"(2I6)")

! this sets the reciprocal range of reflections to be taken into account (in principle)
gamax = 1.0/(2.0/gal/float(npix))
gbmax = 1.0/(2.0/gbl/float(npiy))
gmax = maxval( (/ gamax, gbmax /) )
! suggested g-range:
io_real(1)=gmax
call Message%WriteValue(' For the given sampling, the maximum g value equals : ', io_real, 1,"(F10.4,' nm^-1'//)")

! compute incident wave vector in reciprocal space
call cell%TransSpace(float(k),kstar,'d','r')      ! transform incident direction to reciprocal space
call cell%NormVec(kstar,'r')                      ! normalize reciprocal beam vector
DynWV = kstar/lambda                              ! divide by wavelength and assign


! generate the dynamical matrix for Bloch wave computations (only one single wave vector !)
! we are in zone axis here so there is no tangential component...
gvec = gvectors_T()
FN = enl%k
write (*,*) ' Incident wave vector/foil normal : ', DynWV, FN
call gvec%Initialize_ReflectionList(cell, SG, Diff, FN, DynWV, dmin, verbose=.TRUE.)

! go through the reflist and reset all strong and weak links to null pointers and .FALSE. identifiers;
! at the same time, compute the excitation errors for this incident beam direction kk
rl => gvec%Get_ListHead()
rl => rl%next
do
  if (.not.associated(rl)) EXIT
  nullify(rl%nexts, rl%nextw)
  rl%strong = .FALSE.
  rl%weak = .FALSE.
  rl%sg = Diff%Calcsg(cell,float(rl%hkl(1:3)),DynWV,FN)
  rl => rl%next
end do

! determine strong and weak reflections
nullify(firstw)
nns = 0
nnw = 0
call gvec%Apply_BethePotentials(Diff, firstw, nns, nnw)
io_int(1:2) = (/ nns, nnw /)
call Message%WriteValue(' Total # strong/weak reflections after Bethe potentials : ',io_int,2)

! generate the dynamical matrix
nn = nns
call mem%alloc(DynMat, (/ nn, nn /), 'DynMat', initval=czero)
call gvec%GetDynMat(cell, Diff, firstw, DynMat, nn, nnw)

! allocate various arrays
mem = memory_T()
call mem%alloc(Minp, (/ nn, nn /), 'Minp', initval=czero)
call mem%alloc(W, (/ nn /), 'W', initval = czero)
call mem%alloc(CG, (/ nn, nn /), 'CG', initval = czero)
call mem%alloc(alpha, (/ nn /), 'alpha', initval = czero)
call mem%alloc(CGinv, (/ nn, nn /), 'CGinv', initval=czero)
call mem%alloc(diag, (/ nn /), 'diag', initval = czero)
call mem%alloc(IPIV, (/ nn /), 'IPIV', initval = 0)
call mem%alloc(amp, (/ numt, nn /), 'amp', initval = czero)
call mem%alloc(pw, (/ nn /), 'pw', initval = czero)
call mem%alloc(ew, (/ numt, npix, npiy /), 'ew', initval = cmplx(0.0,0.0))

! solve the eigenvalue problem
Minp = DynMat
call Diff%BWsolve(Minp, W, CG, CGinv, nn, IPIV)
call Message%printMessage(' Eigenvalues/vectors computed, beginning exit wave computation')

! the alpha coefficients are in the first column of the inverse matrix
W = cPi*W*lambda
do i=1,nn
 alpha(i) = CGinv(i,1)
end do
! loop over all thicknesses
do i=1,enl%numthick
 th = dble(enl%startthick+float(i-1) * enl%thickinc)
 diag(1:nn)=exp(-th*imag(W(1:nn)))*cmplx(cos(th*real(W(1:nn))),sin(th*real(W(1:nn))),dbl)*alpha(1:nn)
 do j=1,nn
  amp(i,j) = sum(CG(j,1:nn)*diag(1:nn))
 end do 
end do
! loop over all pixels and compute the wave function for each thickness
reflist => gvec%get_ListHead()
frac = 0.05
ik=0
numk = npix*npiy

! precompute the list of k0+g vectors for the plane waves
call cell%TransSpace(float(enl%k),DynWV,'d','c')
call cell%NormVec(DynWV,'c')
DynWV = DynWV/lambda
call cell%TransSpace(DynWV, DynWV, 'c', 'r')
write (*,*) 'incident wave in reciprocal space reference frame : ', DynWV, 1.0/lambda

call mem%alloc(kg, (/ 3, nn /), 'kg', initval=0.D0)
reflist => gvec%get_ListHead()
rltmpa => reflist%next
do i=1,nn 
  gg = dble(rltmpa%hkl)
  ! glen = cell%CalcLength(gg,'r')
  ! call cell%NormVec(gg,'r')
  ! kg(1:3,i) = DynWV + gg * glen
  kg(1:3,i) = dble(DynWV) + gg
  rltmpa => rltmpa%nexts
end do

gad = gad/float(npix)
gbd = gbd/float(npiy)
reflist => gvec%get_ListHead()
do ip=1,npix
 do jp=1,npiy
  ik=ik+1
  pixpos = dble(ip-1)*gad + dble(jp-1)*gbd
! precompute the plane waves for this location
  ! rltmpa => reflist%next
  do i=1,nn
   ! arg = 2.0*cPi*dot_product(DynWV+dble(rltmpa%hkl),pixpos)
   arg = 2.0*cPi*sum(kg(1:3,i)*pixpos(1:3))
   pw(i) = cmplx(cos(arg),sin(arg),dbl)
   ! rltmpa => rltmpa%nexts
  end do
! loop over all thicknesses
  do i=1,enl%numthick
! loop over all reflections to get the exit wave at the selected location
   ew(i,ip,jp) = sum(amp(i,1:nn)*pw(1:nn))
  end do
! update computation progress
  if (float(ik)/float(numk) .gt. frac) then
   io_int(1) = int(100.0*frac) 
   call Message%WriteValue(' ',io_int, 1,"(1x,I3,' percent completed')") 
   frac = frac + 0.05
  end if
 end do
end do
write (*,*) ' Range ew ', minval(abs(ew)**2), maxval(abs(ew)**2)

! get rid of useless variables
call mem%dealloc(W, 'W')
call mem%dealloc(CG, 'CG')
call mem%dealloc(alpha, 'alpha')
call mem%dealloc(CGinv, 'CGinv')
call mem%dealloc(Minp, 'Minp')
call mem%dealloc(diag, 'diag')
call mem%dealloc(IPIV, 'IPIV')
call mem%dealloc(amp, 'amp')
call mem%dealloc(pw, 'pw')
call mem%dealloc(kg, 'kg')

! and store the exit wave functions for all thicknesses in an HDF5 file;
! also save the intensities as tiff files...
call mem%alloc(inten, (/ npix, npiy /), 'inten')
do i=1,numt
  inten = abs(ew(i,:,:))**2
  mi = minval(inten)
  ma = maxval(inten)

write (*,*) mi, ma 

  TIFF_image = nint(255.0*(inten-mi)/(ma-mi))

  ! set up the image_t structure
  im = image_t(TIFF_image)
  if(im%empty()) call Message%printMessage("EMgetADP","failed to convert array to image")

  write (str,"(I3.3)") i
  TIFF_filename = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(enl%tiffprefix)//'_'//str//'.tiff'

  ! create the file
  call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
  if(0.ne.iostat) then
    call Message%printMessage(" failed to write image to file : "//iomsg)
  else
    call Message%printMessage(' systematic row image written to '//trim(TIFF_filename))
  end if 

end do

call mem%dealloc(ew, 'ew')
call mem%dealloc(inten, 'inten')

end associate

end subroutine BWEW_



end module mod_BWEW