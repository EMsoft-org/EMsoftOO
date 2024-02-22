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

module mod_BW
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/21/24
  !!
  !! class definitions used by the EMTBSRBW program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMTBSRBW program
type, public :: BWNameListType
  integer(kind=irg)  :: g(3)
  integer(kind=irg)  :: k(3)
  integer(kind=irg)  :: f(3)
  integer(kind=irg)  :: numkt
  integer(kind=irg)  :: maxng
  real(kind=sgl)     :: ktmax
  real(kind=sgl)     :: voltage
  character(2)       :: TBSR
  character(fnlen)   :: xtalname
  character(fnlen)   :: outname
end type BWNameListType

! class definition for the two beam Bloch wave program
type, public :: TBSRBW_T
private 
  character(fnlen)        :: nmldeffile = 'EMTBSRBW.nml'
  type(BWNameListType)    :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: TBSRBW_
  procedure, pass(self) :: setg_
  procedure, pass(self) :: getg_
  procedure, pass(self) :: setk_
  procedure, pass(self) :: getk_
  procedure, pass(self) :: setf_
  procedure, pass(self) :: getf_
  procedure, pass(self) :: setnumkt_
  procedure, pass(self) :: getnumkt_
  procedure, pass(self) :: setmaxng_
  procedure, pass(self) :: getmaxng_
  procedure, pass(self) :: setktmax_
  procedure, pass(self) :: getktmax_
  procedure, pass(self) :: setvoltage_
  procedure, pass(self) :: getvoltage_
  procedure, pass(self) :: setxtalname_
  procedure, pass(self) :: getxtalname_
  procedure, pass(self) :: setoutname_
  procedure, pass(self) :: getoutname_
  procedure, pass(self) :: setTBSR_
  procedure, pass(self) :: getTBSR_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: TBSRBW => TBSRBW_
  generic, public :: setg => setg_
  generic, public :: getg => getg_
  generic, public :: setk => setk_
  generic, public :: getk => getk_
  generic, public :: setf => setf_
  generic, public :: getf => getf_
  generic, public :: setnumkt => setnumkt_
  generic, public :: getnumkt => getnumkt_
  generic, public :: setmaxng => setmaxng_
  generic, public :: getmaxng => getmaxng_
  generic, public :: setktmax => setktmax_
  generic, public :: getktmax => getktmax_
  generic, public :: setvoltage => setvoltage_
  generic, public :: getvoltage => getvoltage_
  generic, public :: setxtalname => setxtalname_
  generic, public :: getxtalname => getxtalname_
  generic, public :: setoutname => setoutname_
  generic, public :: getoutname => getoutname_
  generic, public :: setTBSR => setTBSR_
  generic, public :: getTBSR => getTBSR_
end type TBSRBW_T

! the constructor routine for this class 
interface TBSRBW_T
  module procedure TBSRBW_constructor
end interface TBSRBW_T

contains

!--------------------------------------------------------------------------
type(TBSRBW_T) function TBSRBW_constructor( nmlfile ) result(TBSRBW)
!! author: MDG 
!! version: 1.0 
!! date: 02/21/24
!!
!! constructor for the TBSRBW_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call TBSRBW%readNameList(nmlfile)

end function TBSRBW_constructor

!--------------------------------------------------------------------------
subroutine TBSRBW_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 02/21/24
!!
!! destructor for the TBSRBW_T Class
 
IMPLICIT NONE

type(TBSRBW_T), INTENT(INOUT)  :: self 

call reportDestructor('TBSRBW_T')

end subroutine TBSRBW_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/21/24
!!
!! read the namelist from an nml file for the TBSRBW_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(TBSRBW_T), INTENT(INOUT)       :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg)                    :: g(3)
integer(kind=irg)                    :: k(3)
integer(kind=irg)                    :: f(3)
integer(kind=irg)                    :: numkt
integer(kind=irg)                    :: maxng
real(kind=sgl)                       :: ktmax
real(kind=sgl)                       :: voltage
character(2)                         :: TBSR
character(fnlen)                     :: xtalname
character(fnlen)                     :: outname

namelist /TBSRBWlist/ g, k, f, numkt, maxng, ktmax, voltage, xtalname, outname, TBSR

! set the input parameters to default values 
voltage = 200.0
g = (/ 1, 0, 0 /)
k = (/ 0, 0, 1 /)
f = (/ 0, 0, 1 /)
ktmax = 1.0
maxng = 3
numkt = 256
TBSR = 'TB'
xtalname = 'undefined'
outname = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=TBSRBWlist)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(xtalname).eq.'undefined') then
  call Message%printError('readNameList:',' xtalname file name is undefined in '//nmlfile)
 end if

 if (trim(outname).eq.'undefined') then
  call Message%printError('readNameList:',' outname file name is undefined in '//nmlfile)
 end if
end if

self%nml%voltage = voltage
self%nml%g = g
self%nml%k = k
self%nml%f = f
self%nml%ktmax = ktmax
self%nml%numkt = numkt
self%nml%maxng = maxng
self%nml%xtalname = xtalname
self%nml%outname = outname
self%nml%TBSR = TBSR

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/21/24
!!
!! pass the namelist for the TBSRBW_T Class to the calling program

IMPLICIT NONE 

class(TBSRBW_T), INTENT(INOUT)          :: self
type(BWNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/21/24
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)          :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 2, n_real = 2
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( enl => self%nml )

! create the group for this namelist
groupname = trim(HDFnames%get_NMLlist())
hdferr = HDF%createGroup(groupname)

io_int = (/ enl%numkt, enl%maxng /)
intlist(1) = 'numkt'
intlist(2) = 'maxng'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single integers
io_real = (/ enl%ktmax, enl%voltage /)
reallist(1) = 'ktmax'
reallist(2) = 'voltage'
call HDF%writeNMLreals(io_real, reallist, n_real)

! vectors
dataset = SC_k
hdferr = HDF%writeDatasetIntegerArray(dataset, enl%k, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList_: unable to create k dataset', hdferr)

dataset = SC_fn
hdferr = HDF%writeDatasetIntegerArray(dataset, enl%f, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList_: unable to create fn dataset', hdferr)

dataset = 'g'
hdferr = HDF%writeDatasetIntegerArray(dataset, enl%g, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList_: unable to create g dataset', hdferr)

! write all the strings
dataset = SC_outname
line2(1) = enl%outname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList_: unable to create outname dataset', hdferr)

dataset = SC_xtalname
line2(1) = enl%xtalname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList_: unable to create xtalname dataset', hdferr)

dataset = 'TBSR'
line2(1) = enl%TBSR
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList_: unable to create TBSR dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine setg_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setg_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! set g in the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp(3)

self%nml%g = inp

end subroutine setg_

!--------------------------------------------------------------------------
function getg_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getg_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! get g from the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out(3)

out = self%nml%g

end function getg_

!--------------------------------------------------------------------------
subroutine setk_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setk_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! set k in the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp(3)

self%nml%k = inp

end subroutine setk_

!--------------------------------------------------------------------------
function getk_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getk_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! get k from the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out(3)

out = self%nml%k

end function getk_

!--------------------------------------------------------------------------
subroutine setf_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setf_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! set f in the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp(3)

self%nml%f = inp

end subroutine setf_

!--------------------------------------------------------------------------
function getf_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getf_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! get f from the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out(3)

out = self%nml%f

end function getf_

!--------------------------------------------------------------------------
subroutine setnumkt_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnumkt_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! set numkt in the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%numkt = inp

end subroutine setnumkt_

!--------------------------------------------------------------------------
function getnumkt_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnumkt_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! get numkt from the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%numkt

end function getnumkt_

!--------------------------------------------------------------------------
subroutine setmaxng_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setmaxng_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! set maxng in the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%maxng = inp

end subroutine setmaxng_

!--------------------------------------------------------------------------
function getmaxng_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getmaxng_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! get maxng from the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%maxng

end function getmaxng_

!--------------------------------------------------------------------------
subroutine setktmax_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setktmax_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! set ktmax in the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%ktmax = inp

end subroutine setktmax_

!--------------------------------------------------------------------------
function getktmax_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getktmax_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! get ktmax from the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%ktmax

end function getktmax_

!--------------------------------------------------------------------------
subroutine setvoltage_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setvoltage_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! set voltage in the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%voltage = inp

end subroutine setvoltage_

!--------------------------------------------------------------------------
function getvoltage_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getvoltage_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! get voltage from the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%voltage

end function getvoltage_

!--------------------------------------------------------------------------
subroutine setxtalname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setxtalname_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! set xtalname in the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%xtalname = trim(inp)

end subroutine setxtalname_

!--------------------------------------------------------------------------
function getxtalname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getxtalname_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! get xtalname from the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%xtalname)

end function getxtalname_

!--------------------------------------------------------------------------
subroutine setoutname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setoutname_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! set outname in the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%outname = trim(inp)

end subroutine setoutname_

!--------------------------------------------------------------------------
function getoutname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getoutname_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! get outname from the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%outname)

end function getoutname_

!--------------------------------------------------------------------------
subroutine setTBSR_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setTBSR_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! set TBSR in the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)    :: self
character(2), INTENT(IN)          :: inp

self%nml%TBSR = trim(inp)

end subroutine setTBSR_

!--------------------------------------------------------------------------
function getTBSR_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getTBSR_
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! get TBSR from the TBSRBW_T class

IMPLICIT NONE

class(TBSRBW_T), INTENT(INOUT)  :: self
character(2)                    :: out

out = trim(self%nml%TBSR)

end function getTBSR_

!--------------------------------------------------------------------------
subroutine TBSRBW_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: TBSRBW_
!! author: MDG 
!! version: 1.0 
!! date: 02/21/24
!!
!! perform the computations

use mod_EMsoft
use mod_HDFnames
use mod_kinds
use mod_global
use mod_gvectors
use mod_crystallography
use mod_diffraction
use mod_initializers
use mod_symmetry
use mod_kvectors
use mod_math
use mod_timing
use mod_io
use mod_HDFsupport 
use mod_memory
use mod_initializers
use HDF5
use mod_HDFsupport
use stringconstants

use, intrinsic :: iso_fortran_env
use ISO_C_BINDING

IMPLICIT NONE 

class(TBSRBW_T), INTENT(INOUT)      :: self
type(EMsoft_T), INTENT(INOUT)     :: EMsoft
character(fnlen), INTENT(INOUT)   :: progname 
type(HDFnames_T), INTENT(INOUT)   :: HDFnames

type(Cell_T)                      :: Cell 
type(Spacegroup_T)                :: SG 
type(Diffraction_T)               :: Diff
type(DynType)                     :: Dyn
type(IO_T)                        :: Message
type(Memory_T)                    :: mem
type(gnode)                       :: rlp
type(gvectors_T)                  :: gvec
type(kvectors_T)                  :: kvec
type(HDF_T)                       :: HDF
type(Timing_T)                    :: timer
type(reflisttype),pointer         :: reflist, rltmpa, rl, firstw

real(kind=sgl)                    :: Vmod,Vphase,Vpmod,Vpphase,pre,upzero,find(3), dmin,&
                                     kk,kt(3),kz,io_real(1),pre2,dkt,gg,s,ktmax, duration
real(kind=dbl)                    :: lambda
complex(kind=dbl),allocatable     :: M(:,:),alph(:),CGinv(:,:),Mcp(:,:),CG(:,:),W(:)
complex(kind=dbl),allocatable     :: alpha(:,:),CGarray(:,:,:),Warray(:,:)
real(kind=sgl),allocatable        :: kttb(:), kn(:)
complex(kind=dbl)                 :: czero = cmplx(0.0,0.0,dbl)
integer(kind=irg)                 :: ind(3),ivec(3),ik,izero,IPIV(2),io_int(2),i,j,nn,ns,g(3),k(3),fn(3),&
                                     ir1, ir2, hdferr
character(fnlen)                  :: oname, datagroupname, groupname, dataset 
logical                           :: verbose 
character(11)                     :: dstr
character(15)                     :: tstrb
character(15)                     :: tstre
character(2)                      :: str
character(fnlen,kind=c_char)      :: line2(1)

call openFortranHDFInterface()

associate( enl => self%nml )

! memory handling class
mem = memory_T()

! initialize the HDF class
HDF = HDF_T() 

! initialize the timing routines
timer = Timing_T()
tstrb = timer%getTimeString()

call timer%Time_tick(1)

! extract parameters from the namelist
g = enl%g
k = enl%k
fn = enl%f
ktmax = enl%ktmax
ns = enl%numkt 
oname  = EMsoft%generateFilePath('EMdatapathname',enl%outname)

! are we doing a two-beam or systematic row computation ?
if (enl%TBSR.eq.'TB') then 
  ir1 = 1
  ir2 = 2
  izero = 1
  nn = 2
else ! systematic row case
  ir1 = 1
  ir2 = 2*enl%maxng+1
  izero = enl%maxng+1
  nn = ir2
end if 

! allocate the larger arrays that will be used for the HDF5 output file 
call mem%alloc(kttb, (/ ns /), 'kttb', initval = 0.0)
call mem%alloc(kn, (/ ns /), 'kn', initval = 0.0)
call mem%alloc(Warray, (/ nn, ns /), 'Warray', initval = czero)
call mem%alloc(CGarray, (/ nn, nn, ns /), 'CGarray', initval = czero)
call mem%alloc(alpha, (/ nn, ns /), 'alpha', initval = czero)
! and the regular arrays
call mem%alloc(W, (/ nn, ns /), 'W', initval = czero)
call mem%alloc(M, (/ nn, nn, ns /), 'M', initval = czero)
call mem%alloc(Mcp, (/ nn, nn, ns /), 'Mcp', initval = czero)
call mem%alloc(CG, (/ nn, nn, ns /), 'CG', initval = czero)
call mem%alloc(CGinv, (/ nn, nn, ns /), 'CGinv', initval = czero)
call mem%alloc(alph, (/ nn, ns /), 'alph', initval = czero)

! crystallography section
verbose = .TRUE.

call cell%setFileName(enl%xtalname)
call Diff%setrlpmethod('WK')
! foil normal is assumed to be parallel to k 
Diff%Dyn%FN = dble(fn) 
call Diff%setV(dble(enl%voltage))

! initialize the gvectors list
gvec = gvectors_T()
dmin = 0.05
! get the crystal structure data and initialize the LUTs 
call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, dmin, verbose, noLUT=.TRUE., useHDF=HDF)
lambda = Diff%getWaveLength()

! pre converts from V to U
 pre = 2.0*sngl(cRestmass*cCharge/cPlanck**2)*1.0E-18

! scaling factor for excitation error (2*k_0)
 pre2 = 2.0/sngl(lambda)

! normal aborption potential Uprime_0
 ind = [0,0,0]
 call Diff%CalcUcg(cell, (/ 0, 0, 0 /))
 rlp = Diff%getrlp()
 Vmod = rlp%Vmod
 Vpmod = rlp%Vpmod
 Vphase = rlp%Vphase
 Vpphase = rlp%Vpphase
 upzero = pre*Vpmod

! determine the dynamical matrix M (all but the diagonal)
! i is the row index
 do i=ir1,ir2
  ind = (i-enl%maxng-1)*int(g) 
! j is the column index
  do j=ir1,ir2
   if (j.ne.i) then
    ivec = ind - (j-enl%maxng-1)*int(g)
! use Weickenmeier-Kohl scattering parameters and form factors
    call Diff%CalcUcg(cell, ivec)
    rlp = Diff%getrlp()
    M(i,j) = -cPi*cmplx(-aimag(rlp%qg),real(rlp%qg))*cmplx(0.0,1.0/cPi/lambda)
   end if
  end do
 end do

! next we iterate over all incident beam directions, and for
! each direction we complete and diagonalize the M-matrix 
!
 dkt = 2.0*ktmax/float(ns-1)
 io_real(1) = dkt
 call Message%WriteValue(' Beam tilt step size = ', io_real, 1, "(F8.4)")
 find = float(g)
 kk = cell%CalcLength(float(k),'r')
 gg = cell%CalcLength(find,'r')
 k = k/sngl(lambda)/kk
 kz = 1.0/sngl(lambda)

! loop over the beam directions
 do ik = 1,ns
  if (mod(ik,25).eq.0) then
   io_int(1) = ik
   io_int(2) = ns
   call Message%WriteValue(' ', io_int, 2,"('  -> completed column ',I4,' of ',I4)")
  endif

! rescale the wavevector and foil normal
  kt = k + dkt*(float(ik-ns/2)-0.5)*g
  kk = cell%CalcLength(kt,'r')
  kt = kt/sngl(lambda)/kk

! then complete the diagonal of the M matrix
! i is the row index
  do i=ir1,ir2
   ind = (i-enl%maxng-1)*int(g) 
! get the excitation error
   find = float(ind)
   if (i.eq.izero) then
    s = 0.0
   else
    s = Diff%Calcsg(cell,find,kt,float(fn) )
   endif
! and multiply with 2k_0 and store, along with Uprime_0
   M(i,i) = cmplx(pre2*s,upzero)
  end do
!
! next, compute the eigenvalues and eigenvectors
! using the LAPACK CGEEV, CGETRF, and CGETRI routines

! first, make a copy of M, since BWsolve destroys M
  Mcp = M

! then get the eigenvalues and eigenvectors
  call Diff%BWsolve(Mcp,W,CG,CGinv,nn,IPIV)
  CGarray(:,:,ik) = CG(:,:)

! the alpha coefficients are in the izero column of the inverse matrix
! the minus sign in W(i) stems from the fact that k_n is in the direction
! opposite to the foil normal
  kttb(ik) = dkt*(float(ik-ns/2)-0.5)
  kn(ik) = -sqrt(kz**2-(kttb(ik)*gg)**2)
  Warray(:,ik) = W(:)/cmplx(2.0*kn(ik),0.0)
  do i=1,2
   alpha(i, ik) = CGinv(i,izero)
  end do
end do

call timer%Time_tock(1)
duration = timer%getInterval(1)
call timer%makeTimeStamp()
dstr = timer%getDateString()
tstre = timer%getTimeString()

! store everything in an HDF5 file that can then be read by the EMBWshow program
! for visualizations
! Open an existing file or create a new file using the default properties.
 hdferr =  HDF%createFile(oname)

 ! write the EMheader to the file
 datagroupname = trim(HDFnames%get_ProgramData())
 call HDF%writeEMheader(EMsoft, dstr, tstrb, tstre, progname, datagroupname)

   ! create a namelist group to write all the namelist files into
 groupname = SC_NMLfiles
 hdferr = HDF%createGroup(groupname)

! read the text file and write the array to the file
 dataset = SC_TBSRBWNameList
 hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%nmldeffile)

 ! leave this group
 call HDF%pop()
  
 ! create a namelist group to write all the namelist files into
 hdferr = HDF%createGroup(HDFnames%get_NMLparameters())
 call self%writeHDFNameList(HDF, HDFnames)

 ! leave this group
 call HDF%pop()

! then the remainder of the data in a EMData group
 groupname = SC_EMData
 hdferr = HDF%createGroup(groupname)
 hdferr = HDF%createGroup(datagroupname)

! here we distinguish between two beam, systematic row, and (eventually) zone axis cases
 dataset = 'TBSR'
  line2(1) = enl%TBSR
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)

 dataset = SC_xtalname
  line2(1) = trim(enl%xtalname)
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)

 dataset = 'nn'
  hdferr = HDF%writeDatasetInteger(dataset, nn)

 dataset = 'ns'
  hdferr = HDF%writeDatasetInteger(dataset, ns)

 dataset = 'g'
  hdferr = HDF%writeDatasetIntegerArray(dataset, enl%g, 3)

 dataset = 'k'
  hdferr = HDF%writeDatasetIntegerArray(dataset, enl%k, 3)

 dataset = 'kz'
  hdferr = HDF%writeDatasetFloat(dataset, kz)

 dataset = 'kttb'
  hdferr = HDF%writeDatasetFloatArray(dataset, kttb, ns)

 dataset = 'kn'
  hdferr = HDF%writeDatasetFloatArray(dataset, kn, ns)

 dataset = 'W_R'
  hdferr = HDF%writeDatasetDoubleArray(dataset, Warray%re, nn, ns)

 dataset = 'W_I'
  hdferr = HDF%writeDatasetDoubleArray(dataset, Warray%im, nn, ns)

 dataset = 'CG_R'
  hdferr = HDF%writeDatasetDoubleArray(dataset, CGarray%re, nn, nn, ns)

 dataset = 'CG_I'
  hdferr = HDF%writeDatasetDoubleArray(dataset, CGarray%im, nn, nn, ns)

 dataset = 'alpha_R'
  hdferr = HDF%writeDatasetDoubleArray(dataset, alpha%re, nn, ns)

 dataset = 'alpha_I'
  hdferr = HDF%writeDatasetDoubleArray(dataset, alpha%im, nn, ns)

! leave this group and close the file
call HDF%popall()

! deallocate arrays
call mem%dealloc(CGarray, 'CGarray')
call mem%dealloc(Warray, 'Warray')
call mem%dealloc(kttb, 'kttb')
call mem%dealloc(kn, 'kn')
call mem%dealloc(alpha, 'alpha')
call mem%dealloc(alph, 'alph')
call mem%dealloc(CG, 'CG')
call mem%dealloc(CGinv, 'CGinv')
call mem%dealloc(W, 'W')
call mem%dealloc(M, 'M')
call mem%dealloc(Mcp, 'Mcp')

call closeFortranHDFInterface()

end associate 

end subroutine TBSRBW_

end module mod_BW