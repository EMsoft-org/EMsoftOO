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

module mod_SRCBED
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/16/24
  !!
  !! class definition for the EMSRCBED program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMSRCBED program
type, public :: SRCBEDNameListType
  integer(kind=irg) :: nthreads
  integer(kind=irg) :: SRG(3)
  integer(kind=irg) :: SRK(3)
  integer(kind=irg) :: Grange
  integer(kind=irg) :: numrows
  real(kind=sgl)    :: ktonG(10)
  real(kind=sgl)    :: camlen
  real(kind=sgl)    :: voltage
  real(kind=sgl)    :: convergence
  real(kind=sgl)    :: thick
  character(fnlen)  :: outname
  character(fnlen)  :: xtalname
  character(fnlen)  :: tiffprefix
end type SRCBEDNameListType

! class definition
type, public :: SRCBED_T
private 
  character(fnlen)          :: nmldeffile = 'EMSRCBED.nml'
  type(SRCBEDNameListType)  :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: SRCBED_
  procedure, pass(self) :: setnthreads_
  procedure, pass(self) :: getnthreads_
  procedure, pass(self) :: setSRG_
  procedure, pass(self) :: getSRG_
  procedure, pass(self) :: setSRK_
  procedure, pass(self) :: getSRK_
  procedure, pass(self) :: setGrange_
  procedure, pass(self) :: getGrange_
  procedure, pass(self) :: setnumrows_
  procedure, pass(self) :: getnumrows_
  procedure, pass(self) :: setktonG_
  procedure, pass(self) :: getktonG_
  procedure, pass(self) :: setvoltage_
  procedure, pass(self) :: getvoltage_
  procedure, pass(self) :: setcamlen_
  procedure, pass(self) :: getcamlen_
  procedure, pass(self) :: setconvergence_
  procedure, pass(self) :: getconvergence_
  procedure, pass(self) :: setthick_
  procedure, pass(self) :: getthick_
  procedure, pass(self) :: setoutname_
  procedure, pass(self) :: getoutname_
  procedure, pass(self) :: setxtalname_
  procedure, pass(self) :: getxtalname_
  procedure, pass(self) :: settiffprefix_
  procedure, pass(self) :: gettiffprefix_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: SRCBED => SRCBED_
  generic, public :: setnthreads => setnthreads_
  generic, public :: getnthreads => getnthreads_
  generic, public :: setSRG => setSRG_
  generic, public :: getSRG => getSRG_
  generic, public :: setSRK => setSRK_
  generic, public :: getSRK => getSRK_
  generic, public :: setGrange => setGrange_
  generic, public :: getGrange => getGrange_
  generic, public :: setnumrows => setnumrows_
  generic, public :: getnumrows => getnumrows_
  generic, public :: setktonG => setktonG_
  generic, public :: getktonG => getktonG_
  generic, public :: setvoltage => setvoltage_
  generic, public :: getvoltage => getvoltage_
  generic, public :: setcamlen => setcamlen_
  generic, public :: getcamlen => getcamlen_
  generic, public :: setconvergence => setconvergence_
  generic, public :: getconvergence => getconvergence_
  generic, public :: setthick => setthick_
  generic, public :: getthick => getthick_
  generic, public :: setoutname => setoutname_
  generic, public :: getoutname => getoutname_
  generic, public :: setxtalname => setxtalname_
  generic, public :: getxtalname => getxtalname_
  generic, public :: settiffprefix => settiffprefix_
  generic, public :: gettiffprefix => gettiffprefix_
end type SRCBED_T

! the constructor routine for this class 
interface SRCBED_T
  module procedure SRCBED_constructor
end interface SRCBED_T

contains

!--------------------------------------------------------------------------
type(SRCBED_T) function SRCBED_constructor( nmlfile ) result(SRCBED)
!! author: MDG 
!! version: 1.0 
!! date: 02/16/24
!!
!! constructor for the SRCBED_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call SRCBED%readNameList(nmlfile)

end function SRCBED_constructor

!--------------------------------------------------------------------------
subroutine SRCBED_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 02/16/24
!!
!! destructor for the SRCBED_T Class
 
IMPLICIT NONE

type(SRCBED_T), INTENT(INOUT)  :: self 

call reportDestructor('SRCBED_T')

end subroutine SRCBED_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/16/24
!!
!! read the namelist from an nml file for the SRCBED_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(SRCBED_T), INTENT(INOUT)       :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg) :: nthreads
integer(kind=irg) :: SRG(3)
integer(kind=irg) :: SRK(3)
integer(kind=irg) :: Grange
integer(kind=irg) :: numrows
real(kind=sgl)    :: ktonG(10)
real(kind=sgl)    :: voltage
real(kind=sgl)    :: camlen
real(kind=sgl)    :: convergence
real(kind=sgl)    :: thick
character(fnlen)  :: outname
character(fnlen)  :: xtalname
character(fnlen)  :: tiffprefix

namelist /SRCBEDlist/ nthreads, SRG, SRK, Grange, numrows, ktonG, voltage, xtalname, convergence, &
                        thick, outname, tiffprefix, camlen

nthreads = 6
voltage = 200.0
camlen = 1000.0
SRG = (/ 1, 0, 0 /)
SRK = (/ 0, 0, 1 /)
Grange = 4
numrows = 10
ktonG = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
convergence = 5.0
thick = 100.0
xtalname = 'undefined'
outname = 'undefined'
tiffprefix = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=SRCBEDlist)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(xtalname).eq.'undefined') then
  call Message%printError('readNameList:',' xtalname file name is undefined in '//nmlfile)
 end if

 if (trim(outname).eq.'undefined') then
  call Message%printError('readNameList:',' outname file name is undefined in '//nmlfile)
 end if

 if (trim(tiffprefix).eq.'undefined') then
  call Message%printError('readNameList:',' tiffprefix file name is undefined in '//nmlfile)
 end if
end if

self%nml%nthreads = nthreads
self%nml%SRG = SRG
self%nml%SRK = SRK
self%nml%Grange = Grange
self%nml%numrows = numrows
self%nml%ktonG = ktonG
self%nml%voltage = voltage
self%nml%camlen = camlen
self%nml%convergence = convergence
self%nml%thick = thick
self%nml%outname = outname
self%nml%xtalname = xtalname
self%nml%tiffprefix = tiffprefix

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/16/24
!!
!! pass the namelist for the SRCBED_T Class to the calling program

IMPLICIT NONE 

class(SRCBED_T), INTENT(INOUT)          :: self
type(SRCBEDNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine setnthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnthreads_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! set nthreads in the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%nthreads = inp

end subroutine setnthreads_

!--------------------------------------------------------------------------
function getnthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnthreads_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! get nthreads from the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%nthreads

end function getnthreads_

!--------------------------------------------------------------------------
subroutine setSRG_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setSRG_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! set SRG in the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp(3)

self%nml%SRG = inp

end subroutine setSRG_

!--------------------------------------------------------------------------
function getSRG_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getSRG_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! get SRG from the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out(3)

out = self%nml%SRG

end function getSRG_

!--------------------------------------------------------------------------
subroutine setSRK_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setSRK_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! set SRK in the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp(3)

self%nml%SRK = inp

end subroutine setSRK_

!--------------------------------------------------------------------------
function getSRK_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getSRK_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! get SRK from the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out(3)

out = self%nml%SRK

end function getSRK_

!--------------------------------------------------------------------------
subroutine setGrange_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setGrange_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! set Grange in the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%Grange = inp

end subroutine setGrange_

!--------------------------------------------------------------------------
function getGrange_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getGrange_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! get Grange from the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%Grange

end function getGrange_

!--------------------------------------------------------------------------
subroutine setnumrows_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnumrows_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! set numrows in the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%numrows = inp

end subroutine setnumrows_

!--------------------------------------------------------------------------
function getnumrows_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnumrows_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! get numrows from the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%numrows

end function getnumrows_

!--------------------------------------------------------------------------
subroutine setktonG_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setktonG_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! set ktonG in the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)   :: self
real(kind=sgl), INTENT(IN)       :: inp(10)

self%nml%ktonG = inp

end subroutine setktonG_

!--------------------------------------------------------------------------
function getktonG_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getktonG_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! get ktonG from the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
real(kind=sgl),allocatable         :: out(:)

allocate( out(self%nml%numrows))
out(1:self%nml%numrows) = self%nml%ktonG(1:self%nml%numrows)

end function getktonG_

!--------------------------------------------------------------------------
subroutine setvoltage_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setvoltage_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! set voltage in the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%voltage = inp

end subroutine setvoltage_

!--------------------------------------------------------------------------
function getvoltage_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getvoltage_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! get voltage from the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%voltage

end function getvoltage_

!--------------------------------------------------------------------------
subroutine setcamlen_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setcamlen_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! set camlen in the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%camlen = inp

end subroutine setcamlen_

!--------------------------------------------------------------------------
function getcamlen_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getcamlen_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! get camlen from the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%camlen

end function getcamlen_

!--------------------------------------------------------------------------
subroutine setconvergence_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setconvergence_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! set convergence in the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%convergence = inp

end subroutine setconvergence_

!--------------------------------------------------------------------------
function getconvergence_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getconvergence_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! get convergence from the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%convergence

end function getconvergence_

!--------------------------------------------------------------------------
subroutine setthick_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setthick_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! set thick in the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%thick = inp

end subroutine setthick_

!--------------------------------------------------------------------------
function getthick_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getthick_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! get thick from the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%thick

end function getthick_

!--------------------------------------------------------------------------
subroutine setoutname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setoutname_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! set outname in the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%outname = trim(inp)

end subroutine setoutname_

!--------------------------------------------------------------------------
function getoutname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getoutname_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! get outname from the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%outname)

end function getoutname_

!--------------------------------------------------------------------------
subroutine setxtalname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setxtalname_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! set xtalname in the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%xtalname = trim(inp)

end subroutine setxtalname_

!--------------------------------------------------------------------------
function getxtalname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getxtalname_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! get xtalname from the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%xtalname)

end function getxtalname_

!--------------------------------------------------------------------------
subroutine settiffprefix_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: settiffprefix_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! set tiffprefix in the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%tiffprefix = trim(inp)

end subroutine settiffprefix_

!--------------------------------------------------------------------------
function gettiffprefix_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: gettiffprefix_
!! author: MDG
!! version: 1.0
!! date: 02/16/24
!!
!! get tiffprefix from the SRCBED_T class

IMPLICIT NONE

class(SRCBED_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%tiffprefix)

end function gettiffprefix_

!--------------------------------------------------------------------------
subroutine SRCBED_( self, EMsoft, progname )
!DEC$ ATTRIBUTES DLLEXPORT :: SRCBED_
!! author: MDG 
!! version: 1.0 
!! date: 02/16/24
!!
!! perform the computations

use mod_kinds
use mod_global
use mod_EMsoft
use mod_crystallography
use mod_diffraction
use mod_symmetry
use mod_math
use mod_io
use mod_HDFsupport 
use mod_memory
use mod_timing
use mod_initializers
use mod_image

use, intrinsic :: iso_fortran_env

IMPLICIT NONE

class(SRCBED_T),INTENT(INOUT)     :: self
type(EMsoft_T),INTENT(INOUT)      :: EMsoft 
character(fnlen),INTENT(IN)       :: progname

type(Cell_T)                      :: Cell 
type(Spacegroup_T)                :: SPG 
type(Diffraction_T)               :: Diff
type(DynType)                     :: Dyn
type(IO_T)                        :: Message
type(Memory_T)                    :: mem
type(Timing_T)                    :: timer
type(gnode)                       :: rlp

real(kind=dbl)                    :: laL,z0,alp,thc,thb,omega_c,omega_min,omega_max,dmin, mi, ma, camlen, Upz,&
                                     dom,glen,xgpz, io_real(2),sc,omega,exer,sl,thr,zmax,att,gc,gci, lambda
integer(kind=irg)                 :: g(3),ira,dpcnt,ppi,io_int(2),nn,izero,npix,i,j,numi,n,l,ll,np2,nps,nsl,&
                                     is,iq,npx,npy,ipos,istart,istop,rowmax,imo,k, wpix,hkl(3), ii, imanum
character(1)                      :: ans,c
complex(kind=dbl)                 :: czero = cmplx(0.D0,0.D0,dbl), cone = cmplx(1.D0,0.D0,dbl)
logical                           :: overlap,np
character(fnlen)                  :: TIFF_filename
character(3)                      :: tnum
real(kind=dbl),allocatable        :: row(:,:), kt(:)
real(kind=dbl),parameter          :: xoff(0:5)=(/0.0,3.3125,0.0,3.3125,0.0,3.3125/),yoff(0:5)=(/6.0,6.0,3.0,3.0,0.0,0.0/)
complex(kind=dbl),allocatable     :: SMz(:,:,:),disk(:,:,:),p(:),Sphiz(:,:), DHWMz(:,:)
complex(kind=dbl),allocatable     :: q(:,:),qin(:,:),qout(:,:),r(:,:), Az(:,:)

! declare variables for use in object oriented image module
integer                           :: iostat
character(len=128)                :: iomsg
logical                           :: isInteger
type(image_t)                     :: im
integer(int8)                     :: i8 (3,4)
integer(int8), allocatable        :: TIFF_image(:,:)

call openFortranHDFInterface()

associate( enl => self%nml )

! first get the crystal data and microscope voltage
dmin = 0.05D0   ! should be set in namelist file
call Diff%setV( dble(enl%voltage) )
call cell%setFileName( enl%xtalname )
call Initialize_Cell(cell, Diff, SPG, Diff%Dyn, EMsoft, sngl(dmin), noLUT=.TRUE. )
call cell%calcPositions(SPG, 'v')
call Diff%setrlpmethod('WK')
lambda = Diff%getWaveLength()

! normal aborption factor
call Diff%CalcUcg(cell, (/ 0, 0, 0 /))
rlp = Diff%getrlp()
xgpz= rlp%xgp
Upz = rlp%Upmod

! general parameters
g = enl%SRG
glen = cell%CalcLength(float(g),'r')
z0 = enl%thick
rowmax = enl%numrows
mem = memory_T()
kt = self%getktonG_()

! various angles in milli-radians
thb = Diff%CalcDiffAngle(cell, g)*0.5D0*1000.D0
thc = enl%convergence

io_real(1) = thb 
io_real(2) = thc
call Message%WriteValue(' Bragg and convergence angles [mrad] : ', io_real, 2)

if (thc.ge.thb) call Message%printMessage(' Warning: the diffraction disks appear to overlap...')

! camera length
camlen = enl%camlen
laL = lambda * camlen
io_real(1) = lambda
call Message%WriteValue(' wavelength [nm] = ', io_real, 1, "(F10.6)")
io_real(1) = camlen
call Message%WriteValue('  L         [mm] = ', io_real, 1, "(f10.2)")
io_real(1) = laL
call Message%WriteValue(' camera length lambda*L [mm nm] = ', io_real, 1, "(f10.5)")

! determine the number of pixels for a diffraction disk
gc = thc/lambda     ! radius of disk in nm^-1

! scale bar (sc is the conversion factor from nm-1 to inches)
sc = laL/25.4D0
gci = gc*sc                 ! disk radius in inches
! define the number of pixels per inch for the computation
ppi = 300

! for 600 dpi output we need this many pixels per disk diameter (odd number)
npix = int(ppi * 2.D0 * gci/1000.D0)
if (mod(npix,2).eq.0) npix=npix+1
io_int(2) = npix

! next we need the width of the systematic row image; each image row is npix+2 pixels tall
wpix = nint( 2.D0*thb/thc * (npix-1)/2.D0 ) * 2 * enl%Grange + npix  
! the last term manages to get the entire outer disks into the image
io_int(1) = wpix
call Message%WriteValue(' Dimensions of single output image : ', io_int, 2, "(I6,' x ',I6)")
 
! the number of beams is the same for all patterns 
ira = enl%Grange

! total number of beams
nn = 2*ira+1
izero = (nn-1)/2+1
imanum = 0

! absorption factor (amplitude)
att = exp(-cPi*z0/xgpz)
io_real(1) = att
call Message%WriteValue(' Absorption factor (amplitude) : ', io_real, 1)

! allocate dynamical variables 
call mem%alloc(SMz, (/ npix,nn,nn /), 'SMz', initval=czero)
call mem%alloc(Sphiz, (/ npix,nn /), 'Sphiz', initval=czero)
call mem%alloc(Az, (/ nn,nn /), 'Az', initval=czero)
call mem%alloc(p, (/ nn /), 'p', initval=czero)
call mem%alloc(DHWMz, (/ nn,nn /), 'DHWMz', initval=czero)
call mem%alloc(q, (/ nn,nn /), 'q', initval=czero)
call mem%alloc(qin, (/ nn,nn /), 'qin', initval=czero)
call mem%alloc(qout, (/ nn,nn /), 'qout', initval=czero)
call mem%alloc(r, (/ nn,nn /), 'r', initval=czero)
call mem%alloc(disk, (/ nn,npix,npix /), 'disk', initval=czero)
call mem%alloc(row, (/ wpix, npix /), 'row', initval = 0.D0)

allocate(TIFF_image(wpix, npix))

! compute the complex DHW matrix
do i=1,nn
 do j=1,nn
  hkl=(-ira+i-1)*g-(-ira+j-1)*g
  if (i.ne.j) then
   call Diff%CalcUcg(cell, hkl)
   rlp = Diff%getrlp()
   DHWMz(i,j) = cPi*cmplx(-aimag(rlp%qg),real(rlp%qg),dbl)
  else
   DHWMz(i,j) = cmplx(0.0,0.0,dbl)
  endif
 end do
end do
call Message%printMessage(' Darwin-Howie-Whelan matrix initialized')

!========================
! main loop over patterns
 do ii = 1, enl%numrows

! convert k_t to the alp and omega angles
  alp = -2.D0*kt(ii)*thb/1000.D0
  omega_c = cPi*0.5D0+alp
  omega_min = omega_c - thc/1000.D0
  omega_max = omega_c + thc/1000.D0

! step size in omega angle
  dom = (omega_max - omega_min)/float(npix-1)

! loop over all incident beam directions
  numi = 0
  do j=1,npix

! set omega angle
   omega = (omega_min+float(j-1)*dom) 
   do i=1,nn
    n = -ira+i-1
! exer = excitation error
    ! exer = -n*glen*cos(omega)-(1.D0-sqrt(1.D0-(n*lambda*glen*sin(omega))**2))/lambda
    exer = -n*glen*0.5D0*(2.D0*cos(omega)+n*lambda*glen)/(1.D0+n*lambda*glen*cos(omega))
    DHWMz(i,i)=cmplx(0.D0,2.D0*cPi*exer,dbl)
   end do

! compute the first slice scattering matrix 
   sl = z0/100.D0
   nsl = 100

! make sure slice thickness does not exceed 0.4 nm (somewhat arbitrary)
   do while (sl.gt.0.4D0)
     sl = sl/2.D0
     nsl = nsl*2
   end do

! use the MatrixExponential routine from the math module to compute the 
! starting scattering matrix. 
   call MatrixExponential(DHWMz, Az, dble(sl), 'Pade', nn)

   do i=1,nn
    SMz(j,i,1:nn) = Az(i,1:nn)
   end do

  end do

! this completes the first slice scattering matrix; next multiply
! it with itself to the desired thickness
! the izero column of SMz is the actual initial wavefunction at z=sl
  Sphiz(1:npix,1:nn)=SMz(1:npix,1:nn,izero)

! loop to maximum thickness
  do l=2,nsl
   do ll=1,npix
    p(1:nn)=cmplx(0.D0,0.D0,dbl)
    do i=1,nn
     p(i)=p(i)+sum(SMz(ll,i,1:nn)*Sphiz(ll,1:nn))
    end do
    Sphiz(ll,1:nn)=p(1:nn)
   end do
  end do

! fill the disk variable with amplitudes; fill vertical columns, then apply a mask
  do j=1,npix
   do i=1,nn
    disk(i,j,1:npix) = Sphiz(j,i)
   end do
  end do

! next we apply a circular mask (using fourfold symmetry)
  np2 = npix/2
  nps = np2**2
  do i=1,npix/2
   is = (i-np2)**2
   do j=1,npix/2
    iq = is+(j-np2)**2
    if (iq.gt.nps) then 
     do k=1,nn
      disk(k,i,j) = czero
      disk(k,npix+1-i,j) = czero
      disk(k,i,npix+1-j) = czero
      disk(k,npix+1-i,npix+1-j) = czero
     end do
    end if
   end do
  end do

  disk = disk * cmplx(att,0.D0) 

! the disks may overlap, but we are adding intensities, not amplitudes
! so the overlap does not matter.
! first the zero order beam
   row = 0.D0
   ipos = wpix/2+np2+1
   do i=1,npix
    row(ipos-i,1:npix) = abs(disk(izero,i,1:npix))**2
   end do

! negative disks
   do j=1,nn/2
     ipos = wpix/2 - np2 - 1 + int(glen*(-ira+j-1)*ppi*sc)
      do i=1,npix
       if (ipos+i.gt.0) row(ipos+i,1:npix) = row(ipos+i,1:npix) + abs(disk(j,i,1:npix))**2
      end do
   end do

! positive disks
   do j=nn/2+2,nn
    ipos = wpix/2 - np2 - 1 + int(glen*(-ira+j-1)*ppi*sc)
     do i=1,npix
      if (ipos+i.lt.wpix) row(ipos+i,1:npix) = row(ipos+i,1:npix) + abs(disk(j,i,1:npix))**2
     end do
   end do

! deallocate the disk variable
  where(row.lt.1.D-5) row = 1.D-5
  row = log10(1.D0/row)
  rowmax = maxval(row)
  row = row/rowmax 
  row = 1.D0-row

! write the systematic row image to a tiff file
  mi = minval(row)
  ma = maxval(row)

  TIFF_image = nint(255.0*(row-mi)/(ma-mi))

  ! set up the image_t structure
  im = image_t(TIFF_image)
  if(im%empty()) call Message%printMessage("EMgetADP","failed to convert array to image")

  write (tnum,"(I3.3)") ii
  TIFF_filename = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(enl%tiffprefix)//'_'//tnum//'.tiff'

  ! create the file
  call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
  if(0.ne.iostat) then
    call Message%printMessage(" failed to write image to file : "//iomsg)
  else
    call Message%printMessage(' systematic row image written to '//trim(TIFF_filename))
  end if 

end do  ! main loop over patterns

! deallocate all arrays
  call mem%dealloc(DHWMz, 'DHWMz')
  call mem%dealloc(Az, 'Az')
  call mem%dealloc(SMz, 'Smz')
  call mem%dealloc(Sphiz, 'Sphiz')
  call mem%dealloc(p, 'p')
  call mem%dealloc(q, 'q')
  call mem%dealloc(qin, 'qin')
  call mem%dealloc(qout, 'qout')
  call mem%dealloc(r, 'r')
  call mem%dealloc(row, 'row')
  call mem%dealloc(disk, 'disk')

end associate

end subroutine SRCBED_

end module mod_SRCBED