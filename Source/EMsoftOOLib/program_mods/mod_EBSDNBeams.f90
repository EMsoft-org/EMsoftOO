! ###################################################################
! Copyright (c) 2013-2026, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_EBSDNBeams
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/25
  !!
  !! class definition for the EMEBSDNBeams program

use mod_kinds
use mod_global
use mod_MPfiles
use mod_platformsupport

IMPLICIT NONE
private

type EBSDNBeamsNameListType
 real(kind=sgl)         :: dmin
 integer(kind=irg)      :: npx
 integer(kind=irg)      :: nthreads
 character(fnlen)       :: energyfile 
 character(fnlen)       :: tiffprefix
 character(fnlen)       :: BetheParametersFile
end type EBSDNBeamsNameListType

! class definition
type, public :: EBSDNBeams_T
  private
    character(fnlen)              :: nmldeffile = 'EMEBSDNBeams.nml'
    type(EBSDNBeamsNameListType)  :: nml

  contains
  private

    procedure, pass(self) :: readNameList_
    procedure, pass(self) :: getNameList_
    procedure, pass(self) :: writeHDFNameList_
    procedure, pass(self) :: EBSDNBeams_
    procedure, pass(self) :: setdmin_
    procedure, pass(self) :: getdmin_
    procedure, pass(self) :: setnpx_
    procedure, pass(self) :: getnpx_
    procedure, pass(self) :: setnthreads_
    procedure, pass(self) :: getnthreads_
    procedure, pass(self) :: setenergyfile_
    procedure, pass(self) :: getenergyfile_
    procedure, pass(self) :: settiffprefix_
    procedure, pass(self) :: gettiffprefix_
    procedure, pass(self) :: setBetheParametersFile_
    procedure, pass(self) :: getBetheParametersFile_
    
    generic, public :: getNameList => getNameList_
    generic, public :: readNameList => readNameList_
    generic, public :: writeHDFNameList => writeHDFNameList_
    generic, public :: EBSDNBeams => EBSDNBeams_
    generic, public :: setdmin => setdmin_
    generic, public :: getdmin => getdmin_
    generic, public :: setnpx => setnpx_
    generic, public :: getnpx => getnpx_
    generic, public :: setnthreads => setnthreads_
    generic, public :: getnthreads => getnthreads_
    generic, public :: setenergyfile => setenergyfile_
    generic, public :: getenergyfile => getenergyfile_
    generic, public :: settiffprefix => settiffprefix_
    generic, public :: gettiffprefix => gettiffprefix_
    generic, public :: setBetheParametersFile => setBetheParametersFile_
    generic, public :: getBetheParametersFile => getBetheParametersFile_
 end type EBSDNBeams_T

! the constructor routine for this class
interface EBSDNBeams_T
  module procedure EBSDNBeams_constructor
end interface EBSDNBeams_T

contains

!--------------------------------------------------------------------------
type(EBSDNBeams_T) function EBSDNBeams_constructor( nmlfile ) result(EBSDNBeams)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDNBeams_constructor
!! author: MDG
!! version: 1.0
!! date: 01/07/25
!!
!! constructor for the EBSDNBeams_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call EBSDNBeams%readNameList(nmlfile)

end function EBSDNBeams_constructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 01/07/25
!!
!! read the namelist from an nml file for the EBSDNBeams_T Class

use mod_io

IMPLICIT NONE

class(EBSDNBeams_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)                 :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)                 :: initonly
 !! fill in the default values only; do not read the file

type(IO_T)                                  :: Message
logical                                     :: skipread = .FALSE.

real(kind=sgl)                              :: dmin
integer(kind=irg)                           :: npx
integer(kind=irg)                           :: nthreads
character(fnlen)                            :: energyfile 
character(fnlen)                            :: tiffprefix
character(fnlen)                            :: BetheParametersFile

! define the IO namelist to facilitate passing variables to the program.
namelist /EBSDNBeams/ dmin, npx, nthreads, energyfile, BetheParametersFile, tiffprefix

! set the input parameters to default values (except for xtalname, which must be present)
dmin = 0.05                     ! smallest d-spacing [nm]
npx = 500                       ! Nx pixels (total = 2Nx+1)
nthreads = 1
energyfile = 'undefined'        ! default filename for z_0(E_e) data from EMMC Monte Carlo simulations
BetheParametersFile='BetheParameters.nml'
tiffprefix = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=EBSDNBeams)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(energyfile).eq.'undefined') then
  call Message%printError('readNameList:',' output (energy) file name is undefined in '//nmlfile)
 end if
end if

! if we get here, then all appears to be ok, and we need to fill in the nml fields
self%nml%npx = npx
self%nml%nthreads = nthreads
self%nml%dmin = dmin
self%nml%energyfile = energyfile
self%nml%tiffprefix = tiffprefix
self%nml%BetheParametersFile = BetheParametersFile

end subroutine readNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG
!! version: 1.0
!! date: 01/07/25
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants

use ISO_C_BINDING

IMPLICIT NONE

class(EBSDNBeams_T), INTENT(INOUT)      :: self
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 2, n_real = 1
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: reallist(n_real)
character(20)                           :: intlist(n_int)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( enl => self%nml )

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! write all the single integers
io_int = (/ enl%npx, enl%nthreads /)
intlist(1) = 'npx'
intlist(2) = 'nthreads'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ enl%dmin /)
reallist(1) = 'dmin'
call HDF%writeNMLreals(io_real, reallist, n_real)

dataset = SC_energyfile
line2(1) = trim(enl%energyfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create energyfile dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 01/07/25
!!
!! pass the namelist for the EBSDNBeams_T Class to the calling program

IMPLICIT NONE

class(EBSDNBeams_T), INTENT(INOUT)          :: self
type(EBSDNBeamsNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine setdmin_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdmin_
!! author: MDG
!! version: 1.0
!! date: 01/07/25
!!
!! set dmin in the EBSDNBeams_T class

IMPLICIT NONE

class(EBSDNBeams_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%dmin = inp

end subroutine setdmin_

!--------------------------------------------------------------------------
function getdmin_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdmin_
!! author: MDG
!! version: 1.0
!! date: 01/07/25
!!
!! get dmin from the EBSDNBeams_T class

IMPLICIT NONE

class(EBSDNBeams_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%dmin

end function getdmin_

!--------------------------------------------------------------------------
subroutine setnpx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnpx_
!! author: MDG
!! version: 1.0
!! date: 01/07/25
!!
!! set npx in the EBSDNBeams_T class

IMPLICIT NONE

class(EBSDNBeams_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%npx = inp

end subroutine setnpx_

!--------------------------------------------------------------------------
function getnpx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnpx_
!! author: MDG
!! version: 1.0
!! date: 01/07/25
!!
!! get npx from the EBSDNBeams_T class

IMPLICIT NONE

class(EBSDNBeams_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%npx

end function getnpx_

!--------------------------------------------------------------------------
subroutine setnthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnthreads_
!! author: MDG
!! version: 1.0
!! date: 01/07/25
!!
!! set nthreads in the EBSDNBeams_T class

IMPLICIT NONE

class(EBSDNBeams_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%nthreads = inp

end subroutine setnthreads_

!--------------------------------------------------------------------------
function getnthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnthreads_
!! author: MDG
!! version: 1.0
!! date: 01/07/25
!!
!! get nthreads from the EBSDNBeams_T class

IMPLICIT NONE

class(EBSDNBeams_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%nthreads

end function getnthreads_

!--------------------------------------------------------------------------
subroutine setenergyfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setenergyfile_
!! author: MDG
!! version: 1.0
!! date: 01/07/25
!!
!! set energyfile in the EBSDNBeams_T class

IMPLICIT NONE

class(EBSDNBeams_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%energyfile = trim(inp)

end subroutine setenergyfile_

!--------------------------------------------------------------------------
function getenergyfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getenergyfile_
!! author: MDG
!! version: 1.0
!! date: 01/07/25
!!
!! get energyfile from the EBSDNBeams_T class

IMPLICIT NONE

class(EBSDNBeams_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%energyfile)

end function getenergyfile_

!--------------------------------------------------------------------------
subroutine settiffprefix_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: settiffprefix_
!! author: MDG
!! version: 1.0
!! date: 01/08/25
!!
!! set tiffprefix in the EBSDNBeams_T class

IMPLICIT NONE

class(EBSDNBeams_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%tiffprefix = trim(inp)

end subroutine settiffprefix_

!--------------------------------------------------------------------------
function gettiffprefix_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: gettiffprefix_
!! author: MDG
!! version: 1.0
!! date: 01/08/25
!!
!! get tiffprefix from the EBSDNBeams_T class

IMPLICIT NONE

class(EBSDNBeams_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%tiffprefix)

end function gettiffprefix_

!--------------------------------------------------------------------------
subroutine setBetheParametersFile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setBetheParametersFile_
!! author: MDG
!! version: 1.0
!! date: 01/07/25
!!
!! set BetheParametersFile in the EBSDNBeams_T class

IMPLICIT NONE

class(EBSDNBeams_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%BetheParametersFile = trim(inp)

end subroutine setBetheParametersFile_

!--------------------------------------------------------------------------
function getBetheParametersFile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getBetheParametersFile_
!! author: MDG
!! version: 1.0
!! date: 01/07/25
!!
!! get BetheParametersFile from the EBSDNBeams_T class

IMPLICIT NONE

class(EBSDNBeams_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%BetheParametersFile)

end function getBetheParametersFile_

!--------------------------------------------------------------------------
subroutine EBSDNBeams_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDNBeams_
!! author: MDG
!! version: 1.0
!! date: 01/07/25
!!
!! compute an EBSD number-of-beams master pattern... could be useful for 4D-EBSD work
!!

use mod_EMsoft
use mod_initializers
use mod_symmetry
use mod_crystallography
use mod_gvectors
use mod_kvectors
use mod_io
use mod_math
use mod_diffraction
use mod_timing
use mod_Lambert
use HDF5
use mod_HDFsupport
use mod_HDFnames
use ISO_C_BINDING
use omp_lib
use mod_OMPsupport
use mod_memory
use mod_notifications
use stringconstants
use mod_MCfiles
use mod_image
use, intrinsic :: iso_fortran_env

IMPLICIT NONE

class(EBSDNBeams_T), INTENT(INOUT)  :: self
type(EMsoft_T), INTENT(INOUT)       :: EMsoft
character(fnlen),INTENT(IN)         :: progname
type(HDFnames_T),INTENT(INOUT)      :: HDFnames

type(Cell_T)            :: cell
type(DynType)           :: Dyn
type(Timing_T)          :: timer
type(IO_T)              :: Message
type(Lambert_T)         :: L
type(HDF_T)             :: HDF
type(SpaceGroup_T)      :: SG
type(Diffraction_T),save:: Diff
type(MCfile_T)          :: MCFT
type(MPfile_T)          :: MPFT
type(kvectors_T)        :: kvec
type(gvectors_T)        :: reflist
type(HDFnames_T)        :: saveHDFnames
type(memory_T)          :: mem, memth

real(kind=dbl)          :: ctmp(192,3), arg, Radius, xyz(3)
integer(HSIZE_T)        :: dims4(4), cnt4(4), offset4(4)
integer(HSIZE_T)        :: dims3(3), cnt3(3), offset3(3)
integer(kind=irg)       :: isym,i,j,ik,npy,ipx,ipy,ipz,debug,iE,izz, izzmax, iequiv(3,48), nequiv, num_el, MCnthreads, & ! counters
                           numk, timestart, timestop, numsites, nthreads, & ! number of independent incident beam directions
                           ir,nat(maxpasym),kk(3), skip, ijmax, one, NUMTHREADS, TID, SamplingType, &
                           numset,n,ix,iy,iz, io_int(6), nns, nnw, nref, Estart, sz(3), &
                           istat,gzero,ic,ip,ikk, totstrong, totweak, jh, ierr, nix, niy, nixp, niyp     ! counters
real(kind=dbl)          :: tpi,Znsq, kkl, DBWF, kin, delta, h, lambda, omtl, srt, dc(3), xy(2), edge, scl, tmp, dx, dxm, dy, dym, &
                           kkk(3), sxy(2) !
real(kind=sgl)          :: io_real(5), selE, kn, FN(3), tstop, nabsl, etotal, density, Ze, at_wt, bp(4)
real(kind=sgl),allocatable      :: EkeVs(:), auxNH(:,:,:), auxSH(:,:,:)  ! results
real(kind=sgl),allocatable      :: mLPNH(:,:,:), mLPSH(:,:,:), masterSPNH(:,:,:), masterSPSH(:,:,:), EBSDmaster(:,:,:)
integer(kind=irg),allocatable   :: accum_z(:,:,:,:), minbeam(:), maxbeam(:)
real(kind=dbl),allocatable      :: SGrecip(:,:,:), SGdirec(:,:,:)
complex(kind=dbl)               :: czero
logical                         :: usehex, switchmirror, verbose
character(fnlen)                :: xtalname

! Monte Carlo derived quantities
integer(kind=irg)               :: numEbins, nsx, nsy, hdferr, nlines, lastEnergy    ! variables used in MC energy file
character(fnlen)                :: oldprogname, groupname, energyfile, outname, datagroupname, attributename, &
                                   HDF_FileVersion, fname
character(8)                    :: MCscversion
character(11)                   :: dstr
character(15)                   :: tstrb
character(15)                   :: tstre
character(4)                    :: number
logical                         :: f_exists, readonly, overwrite=.TRUE., insert=.TRUE., stereog, g_exists, xtaldataread, FL, &
                                   doLegendre, isTKD = .FALSE.
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)
character(fnlen,kind=c_char)                     :: line2(1)

type(gnode),save                :: rlp
real(kind=dbl),allocatable      :: karray(:,:)
integer(kind=irg),allocatable   :: kij(:,:)
complex(kind=dbl),allocatable   :: DynMat(:,:)
character(fnlen)                :: dataset, instring, TIFF_filename
type(MCOpenCLNameListType)      :: mcnl
type(EBSDmasterNameListType)    :: mpnl
type(kvectorlist), pointer      :: ktmp
type(reflisttype), pointer      :: firstw

character(fnlen),ALLOCATABLE    :: MessageLines(:)
integer(kind=irg)               :: NumLines, info
character(fnlen)                :: SlackUsername, exectime
character(100)                  :: c
integer(kind=4)                 :: hnStat

! declare variables for use in object oriented image module
integer(kind=irg),allocatable   :: RGBmap(:,:,:)
integer                         :: iostat
character(len=128)              :: iomsg
logical                         :: isInteger, OPC, PUC
type(image_t)                   :: im
integer(int8), allocatable      :: TIFF_image(:,:)
integer                         :: dim2(2), Pm
integer(c_int32_t)              :: result

!$OMP THREADPRIVATE(Diff)

call openFortranHDFInterface()

! set the HDF group names for this program
HDF = HDF_T()
call MPFT%setModality('EBSD')

! simplify the notation a little
associate( emnl => self%nml, MPDT => MPFT%MPDT )

! initialize the timing routines
timer = Timing_T()
tstrb = timer%getTimeString()

! initialize the memory class 
mem = memory_T()
stereog = .TRUE.

tpi = 2.D0*cPi
czero = cmplx(0.D0,0.D0)

!=============================================
!=============================================
! ---------- read Monte Carlo .h5 output file and extract necessary parameters
! set the HDF group names for reading the MC input file
saveHDFnames = HDFnames  
call HDFnames%set_ProgramData(SC_MCOpenCL)
call HDFnames%set_NMLlist(SC_MCCLNameList)
call HDFnames%set_NMLfilename(SC_MCOpenCLNML)
fname = EMsoft%generateFilePath('EMdatapathname',trim(emnl%energyfile))
call MCFT%setFileName(fname)
call MCFT%readMCfile(HDF, HDFnames, getAccumz=.TRUE.)
mcnl = MCFT%getnml()
call MCFT%copyaccumz(accum_z)

! set the HDFnames to the correct strings for this program
HDFnames = saveHDFnames

nsx = (mcnl%numsx - 1)/2
nsy = nsx
etotal = float(mcnl%totnum_el)

io_int(1) = mcnl%totnum_el
call Message%WriteValue(' --> total number of BSE electrons in MC data set ', io_int, 1)
!=============================================
!=============================================
numEbins = MCFT%getnumEbins()

call mem%alloc( EkeVs, (/ numEbins /), 'EkeVs')
call mem%alloc( minbeam, (/ numEbins /), 'minbeam')
call mem%alloc( maxbeam, (/ numEbins /), 'maxbeam')

do i=1,numEbins
  EkeVs(i) = mcnl%Ehistmin + float(i-1)*mcnl%Ebinsize
end do

!=============================================
! should we create a new file or open an existing file?
!=============================================
energyfile = trim(EMsoft%generateFilePath('EMdatapathname',emnl%energyfile))
outname = trim(energyfile)

saveHDFnames = HDFnames  
call HDFnames%set_ProgramData(SC_EBSDmaster)
call HDFnames%set_NMLlist(SC_EBSDmasterNameList)
call HDFnames%set_NMLfilename(SC_EBSDmasterNML)
call HDFnames%set_Variable(SC_MCOpenCL)

call MPFT%setFileName(energyfile)
call MPFT%readMPfile(HDF, HDFnames, mpnl, getmasterSPNH=.TRUE.)
HDFnames = saveHDFnames

!=============================================
!=============================================
! crystallography section;
verbose = .TRUE.

call cell%setFileName(mcnl%xtalname)
call Diff%setrlpmethod('WK')

call Diff%setV(dble(mcnl%EkeV))
call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, emnl%dmin, verbose, useHDF=HDF)

! check the crystal system and setting; abort the program for trigonal with rhombohedral setting with
! an explanation for the user

if ((SG%getSpaceGroupXtalSystem().eq.5).and.(cell%getLatParm('b').eq.cell%getLatParm('c'))) then
    call Message%printMessage( (/ &
    '                                                                         ', &
    ' ========Program Aborted========                                         ', &
    ' The EBSD master pattern simulation for rhombohedral/trigonal structures ', &
    ' requires that the structure be described using the hexagonal reference  ', &
    ' frame.  Please re-enter the crystal structure in this setting and re-run', &
    ' the Monte Carlo calculation and this master pattern program.            '/) )
    stop
end if

! allocate and compute the Sgh loop-up table
 numset = cell%getNatomtype()
 call Diff%Initialize_SghLUT(cell, SG, emnl%dmin, numset, nat, verbose)

! determine the point group number
 j=0
 do i=1,32
  if (SGPG(i).le.SG%getSpaceGroupNumber()) j=i
 end do
 isym = j

! here is new code dealing with all the special cases (quite a few more compared to the
! Laue group case)...  isym is the point group number. Once the symmetry case has been
! fully determined (taking into account things like 31m and 3m1 an such), then the only places
! that symmetry is handled are the modified Calckvectors routine, and the filling of the modified
! Lambert projections after the dynamical simulation step.  We are also changing the name of the
! sr array (or srhex) to mLPNH and mLPSH (modified Lambert Projection Northern/Southern Hemisphere),
! and we change the output HDF5 file a little as well. We need to make sure that the EMEBSD program
! issues a warning when an old format HDF5 file is read.

! Here, we encode isym into a new number that describes the sampling scheme; the new schemes are
! described in detail in the EBSD manual pdf file.

SamplingType = PGSamplingType(isym)

! next, intercept the special cases (hexagonal vs. rhombohedral cases that require special treatment)
if ((SamplingType.eq.-1).or.(isym.eq.14).or.(isym.eq.26)) then
  SamplingType = SG%getHexvsRho(isym)
end if

! if the point group is trigonal or hexagonal, we need to switch usehex to .TRUE. so that
! the program will use the hexagonal sampling method
usehex = .FALSE.
if ((SG%getSpaceGroupXtalSystem().eq.4).or.(SG%getSpaceGroupXtalSystem().eq.5)) usehex = .TRUE.

! ---------- end of symmetry and crystallography section
!=============================================
!=============================================

!=============================================
!=============================================
! ---------- a couple of initializations
   npy = emnl%npx
   
   gzero = 1  ! index of incident beam
   debug = 0  ! no longer used
! ----------
!=============================================
!=============================================

!=============================================
!=============================================
! ---------- allocate memory for the master patterns
  call mem%alloc(mLPNH, (/ emnl%npx, npy, 1 /), 'mLPNH', startdims=(/-emnl%npx,-npy,1/))
  call mem%alloc(mLPSH, (/ emnl%npx, npy, 1 /), 'mLPSH', startdims=(/-emnl%npx,-npy,1/))
  call mem%alloc(masterSPNH, (/ emnl%npx, npy, 1 /), 'masterSPNH', startdims=(/-emnl%npx,-npy,1/))
  call mem%alloc(masterSPSH, (/ emnl%npx, npy, 1 /), 'masterSPSH', startdims=(/-emnl%npx,-npy,1/))
  call mem%alloc(RGBmap, (/ 3, 2*emnl%npx+1, 2*npy+1 /), 'RGBmap')
  call mem%alloc(EBSDmaster, (/ 2*emnl%npx+1, 2*npy+1, numEbins /), 'EBSDmaster')
  EBSDmaster = MPDT%masterSPNH
  deallocate( MPDT%masterSPNH )

! set various arrays to zero
   mLPNH = 0.0
   mLPSH = 0.0
   masterSPNH = 0.0
   masterSPSH = 0.0
   RGBmap = 0
! ---------- end allocate memory for the master patterns
!=============================================
!=============================================

! force dynamical matrix routine to read new Bethe parameters from file
! this will all be changed with the new version of the Bethe potentials
  call Diff%SetBetheParameters(EMsoft, .FALSE., emnl%BetheParametersFile)

!=============================================
! create or update the HDF5 output file
!=============================================
  HDF = HDF_T()

! Open the energy file
  outname = EMsoft%generateFilePath('EMdatapathname',trim(emnl%energyfile))
  hdferr =  HDF%openFile(outname)

! write the EMheader to the file
  datagroupname = trim(HDFnames%get_ProgramData())
  call HDF%writeEMheader(EMsoft,dstr, tstrb, tstre, progname, datagroupname)

! add the Duration field to the EMheader group
  hdferr = HDF%openGroup(HDFnames%get_EMheader())
  hdferr = HDF%openGroup(HDFnames%get_ProgramData())

dataset = SC_Duration
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  tstop = 0
  if (g_exists) then
    hdferr = HDF%writeDatasetFloat(dataset, tstop, overwrite)
  else
    hdferr = HDF%writeDatasetFloat(dataset, tstop)
  end if
  call HDF%pop()
  call HDF%pop()

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
  hdferr = HDF%openGroup(HDFnames%get_EMData())
  hdferr = HDF%createGroup(HDFnames%get_ProgramData())

! create the EBSDNBeams group and add a HDF_FileVersion attribute to it
  HDF_FileVersion = '4.0'
  HDF_FileVersion = cstringify(HDF_FileVersion)
  attributename = SC_HDFFileVersion
  call H5Aexists_f(HDF%getobjectID(),trim(attributename),g_exists, hdferr)
  if (.not.g_exists) then
    hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)
  end if

! =====================================================
! The following write commands constitute HDF_FileVersion = 4.0
! =====================================================
dataset = SC_xtalname
  allocate(stringarray(1))
  stringarray(1)= trim(mcnl%xtalname)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1, overwrite)
  else
    hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)
  end if

dataset = SC_BetheParameters
  bp(1) = Diff%getBetheParameter('c1')
  bp(2) = Diff%getBetheParameter('c2')
  bp(3) = Diff%getBetheParameter('c3')
  bp(4) = Diff%getBetheParameter('sg')
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetFloatArray(dataset, bp, 4, overwrite)
  else
    hdferr = HDF%writeDatasetFloatArray(dataset, bp, 4)
  end if


dataset = SC_numEbins
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then
      hdferr = HDF%writeDatasetInteger(dataset, numEbins, overwrite)
    else
      hdferr = HDF%writeDatasetInteger(dataset, numEbins)
    end if

dataset = SC_EkeVs
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then
      hdferr = HDF%writeDatasetFloatArray(dataset, EkeVs, numEbins, overwrite)
    else
      hdferr = HDF%writeDatasetFloatArray(dataset, EkeVs, numEbins)
    end if

! create the hyperslabs and write zeroes to them for now
dataset = SC_mLPNH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1 /)
  offset3 = (/ 0, 0, 0 /)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims3, offset3, cnt3, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims3, offset3, cnt3)
  end if

dataset = SC_mLPSH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1 /)
  offset3 = (/ 0, 0, 0 /)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims3, offset3, cnt3, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims3, offset3, cnt3)
  end if

dataset = SC_masterSPNH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1 /)
  offset3 = (/ 0, 0, 0 /)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPNH, dims3, offset3, cnt3, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPNH, dims3, offset3, cnt3)
  end if

dataset = SC_masterSPSH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1 /)
  offset3 = (/ 0, 0, 0 /)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPSH, dims3, offset3, cnt3, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPSH, dims3, offset3, cnt3)
  end if

! =====================================================
! end of HDF_FileVersion = 4.0 write statements
! =====================================================

  call HDF%popall()

!=============================================
!=============================================
Estart = numEbins

!=============================================
!=============================================
! ---------- from here on, we need to repeat the entire computation for each energy value
! so this is where we could in principle implement an OpenMP approach; alternatively,
! we could do the inner loop over the incident beam directions in OpenMP (probably simpler)

! we use two times, one (1) for each individual energy level, the other (2) for the overall time
call timer%Time_tick(2)
reflist = gvectors_T()

! instantiate the memory class for the OpenMP section
memth = memory_T( nt = emnl%nthreads, silent=.TRUE. )

energyloop: do iE=Estart,1,-1
   ! start the energy level timer
   call timer%Time_tick(1)

! print a message to indicate where we are in the computation
   io_int(1)=iE
   io_int(2)=Estart
   call Message%printMessage(' Starting computation for energy bin (in reverse order)', frm = "(/A)",advance="no")
   call Message%WriteValue(' ',io_int,2,"(I4,' of ',I4)",advance="no")
   io_real(1) = EkeVs(iE)
   call Message%WriteValue('; energy [keV] = ',io_real,1,"(F6.2/)")
   selE = EkeVs(iE)

! set the accelerating voltage
   call cell%setFileName(mcnl%xtalname)
   call Diff%setV(dble(EkeVs(iE)))
   call Diff%setrlpmethod('WK')
   if(iE .ne. Estart) then
    verbose = .FALSE.
    call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, emnl%dmin, useHDF=HDF, verbose=verbose)
   end if

!=============================================
! ---------- create the incident beam directions list
! determine all independent incident beam directions (use a linked list starting at khead)
! numk is the total number of k-vectors to be included in this computation;
! note that this needs to be redone for each energy, since the wave vector changes with energy
   kvec = kvectors_T()   ! initialize the wave vector list
   call kvec%set_kinp( (/ 0.D0, 0.D0, 1.D0 /) )
   call kvec%set_ktmax( 0.D0 )
   call kvec%set_SamplingType( SamplingType )

   call kvec%set_mapmode('RoscaLambert')
   if (usehex) then
     call kvec%Calckvectors(cell, SG, Diff, (/ 0.D0, 0.D0, 0.D0 /),emnl%npx,npy, ijmax,usehex)
   else
     call kvec%Calckvectors(cell, SG, Diff, (/ 0.D0, 0.D0, 0.D0 /),emnl%npx,npy, ijmax,usehex)
   end if
   numk = kvec%get_numk()
   io_int(1)=numk
   call Message%WriteValue('# independent beam directions to be considered = ', io_int, 1, "(I8)")

! are using a Hall space group with potentially different setting ?  If so, then we
! must transform the k-vectors to a different reference frame before using them
  if (SG%getuseHallSG().eqv..TRUE.) then 
! point to the first beam direction
write (*,*) 'using kvector transform rule'
    ktmp => kvec%get_ListHead()
    ktmp%k(1:3) = matmul(SG%HallSG%kvec_transform, ktmp%k(1:3))
    do ik=2,numk
      ktmp => ktmp%next
      if (ik.eq.55000) write (*,*) 'original k-210 ', ktmp%k(1:3)
      ktmp%k(1:3) = matmul(SG%HallSG%kvec_transform, ktmp%k(1:3))
      if (ik.eq.55000) write (*,*) 'transformed k-210 ', ktmp%k(1:3)
    end do
  end if

! convert part of the kvector linked list into arrays for OpenMP
  call mem%alloc(karray, (/4,numk/), 'karray')
  call mem%alloc(kij, (/3,numk/), 'kij')
! point to the first beam direction
  ktmp => kvec%get_ListHead()
  ik = 1
! and loop through the list, keeping k, kn, and i,j
  karray(1:3,ik) = ktmp%k(1:3)
  karray(4,ik) = ktmp%kn
  kij(1:3,ik) = (/ ktmp%i, ktmp%j, ktmp%hs /)
  do ik=2,numk
    ktmp => ktmp%next
    karray(1:3,ik) = ktmp%k(1:3)
    karray(4,ik) = ktmp%kn
    kij(1:3,ik) = (/ ktmp%i, ktmp%j, ktmp%hs /)
  end do
! and remove the linked list
  call kvec%Delete_kvectorlist()

  verbose = .FALSE.
  totstrong = 0
  totweak = 0

  scl = float(emnl%npx)

! ---------- end of "create the incident beam directions list"
!=============================================

! here's where we introduce the OpenMP calls, to speed up the overall calculations...

! set the number of OpenMP threads
  call OMP_setNThreads(emnl%nthreads)

! use OpenMP to run on multiple cores ...
!$OMP PARALLEL COPYIN(Diff) &
!$OMP& PRIVATE(ik,FN,TID,kn,ipx,ipy,ipz,ix,iequiv,nequiv,reflist,firstw) &
!$OMP& PRIVATE(kkk,nns,nnw,nref,io_int,L,ierr,sxy)

  NUMTHREADS = OMP_GET_NUM_THREADS()
  TID = OMP_GET_THREAD_NUM()

!$OMP DO SCHEDULE(DYNAMIC,100)
! ---------- and here we start the beam direction loop
   beamloop:do ik = 1,numk

!=============================================
! ---------- create the master reflection list for this beam direction
! Then we must determine the masterlist of reflections (also a linked list);
! This list basically samples a large reciprocal space volume; it does not
! distinguish between zero and higher order Laue zones, since that
! distinction becomes meaningless when we consider the complete
! reciprocal lattice.
     reflist = gvectors_T()
     kkk = karray(1:3,ik)
     FN = kkk

     call reflist%Initialize_ReflectionList(cell, SG, Diff, FN, sngl(kkk), self%nml%dmin, verbose)
     nref = reflist%get_nref()
! ---------- end of "create the master reflection list"
!=============================================

! determine number of strong and weak reflections
     nullify(firstw)
     nns = 0
     nnw = 0
     call reflist%Apply_BethePotentials(Diff, firstw, nns, nnw)
! and we are done... simply store this number nns in the mLPNH and mLPSH arrays as a float

! If we are using the Hall space groups, then we need to compute the indices directly
! from the unit wave vector since the conventional kij indices will be incorrect.
     if (SG%getuseHallSG().eqv..FALSE.) then
       ipx = kij(1,ik)
       ipy = kij(2,ik)
       ipz = kij(3,ik)
     else  ! we are using the Hall space group symbols so extract (ipx, ipy, ipz) from unit k
       ! first transform to the Cartesian frame 
       call cell%TransSpace(karray(1:3,ik), kkk, 'r', 'c') ! then normalize this vector
       call cell%NormVec(kkk,'c')   ! then transform the unit vector into square Lambert components 
       L = Lambert_T( xyzd = kkk )
       ierr = L%LambertSphereToSquare(sxy)
       ipx = nint(sxy(1)*scl)
       ipy = nint(sxy(2)*scl)
       ipz = -1
       if (kkk(3).ge.0.0) then 
        ipz = 1
       end if
     end if 
  !
     if (usehex) then
       call L%Apply3DPGSymmetry(cell,SG,ipx,ipy,ipz,self%nml%npx,iequiv,nequiv,usehex)
     else
       if ((SG%getSpaceGroupNumber().ge.195).and.(SG%getSpaceGroupNumber().le.230)) then
         call L%Apply3DPGSymmetry(cell,SG,ipx,ipy,ipz,self%nml%npx,iequiv,nequiv,cubictype=SamplingType)
       else
         call L%Apply3DPGSymmetry(cell,SG,ipx,ipy,ipz,self%nml%npx,iequiv,nequiv)
       end if
     end if

!$OMP CRITICAL
     do ix=1,nequiv
       if (iequiv(3,ix).eq.-1) mLPSH(iequiv(1,ix),iequiv(2,ix),1) = real(nns)
       if (iequiv(3,ix).eq.1) mLPNH(iequiv(1,ix),iequiv(2,ix),1) = real(nns)
     end do
!$OMP END CRITICAL
     if (mod(ik,5000).eq.0) then
       io_int(1) = ik
       io_int(2) = numk
       call Message%WriteValue('  completed beam direction ',io_int, 2, "(I8,' of ',I8)")
     end if

     call reflist%Delete_gvectorlist()
    end do beamloop
    
! end of OpenMP portion
!$OMP END PARALLEL

! deallocate arrays that will need to be re-allocated in the next cycle
  call mem%dealloc(karray, 'karray')
  call mem%dealloc(kij, 'kij')

  if (usehex) then
! and finally, we convert the hexagonally sampled array to a square Lambert projection which will be used
! for all EBSD pattern interpolations;  we need to do this for both the Northern and Southern hemispheres

! we begin by allocating auxiliary arrays to hold copies of the hexagonal data; the original arrays will
! then be overwritten with the newly interpolated data.
    call mem%alloc(auxNH, (/emnl%npx,npy,1/), 'auxNH', startdims=(/-emnl%npx,-npy,1/))
    call mem%alloc(auxSH, (/emnl%npx,npy,1/), 'auxSH', startdims=(/-emnl%npx,-npy,1/))
    auxNH = mLPNH
    auxSH = mLPSH

    edge = 1.D0 / dble(emnl%npx)
    scl = float(emnl%npx)
    do i=-emnl%npx,emnl%npx
      do j=-npy,npy
! determine the spherical direction for this point
        L = Lambert_T( xyd = (/ dble(i), dble(j) /) * edge )
        ierr = L%LambertSquareToSphere(dc)
! convert direction cosines to hexagonal Lambert projections
        L = Lambert_T( xyzd = dc )
        ierr = L%LambertSphereToHex(xy)
        xy = xy * scl
! interpolate intensity from the neighboring points
        if (ierr.eq.0) then
          nix = floor(xy(1))
          niy = floor(xy(2))
          nixp = nix+1
          niyp = niy+1
          if (nixp.gt.emnl%npx) nixp = nix
          if (niyp.gt.emnl%npx) niyp = niy
          dx = xy(1) - nix
          dy = xy(2) - niy
          dxm = 1.D0 - dx
          dym = 1.D0 - dy
          mLPNH(i,j,1) = auxNH(nix,niy,1)*dxm*dym + auxNH(nixp,niy,1)*dx*dym + &
                               auxNH(nix,niyp,1)*dxm*dy + auxNH(nixp,niyp,1)*dx*dy
          mLPSH(i,j,1) = auxSH(nix,niy,1)*dxm*dym + auxSH(nixp,niy,1)*dx*dym + &
                               auxSH(nix,niyp,1)*dxm*dy + auxSH(nixp,niyp,1)*dx*dy
        end if
      end do
    end do
    call mem%dealloc(auxNH, 'auxNH')
    call mem%dealloc(auxSH, 'auxSH')
  end if

! make sure that the outer pixel rim of the mLPSH patterns is identical to
! that of the mLPNH array.
  mLPSH(-emnl%npx,-emnl%npx:emnl%npx,1) = mLPNH(-emnl%npx,-emnl%npx:emnl%npx,1)
  mLPSH( emnl%npx,-emnl%npx:emnl%npx,1) = mLPNH( emnl%npx,-emnl%npx:emnl%npx,1)
  mLPSH(-emnl%npx:emnl%npx,-emnl%npx,1) = mLPNH(-emnl%npx:emnl%npx,-emnl%npx,1)
  mLPSH(-emnl%npx:emnl%npx, emnl%npx,1) = mLPNH(-emnl%npx:emnl%npx, emnl%npx,1)

! get stereographic projections
  Radius = 1.0
  do i=-emnl%npx,emnl%npx
    do j=-emnl%npx,emnl%npx
      L = Lambert_T( xyd = (/ dble(i), dble(j) /) / dble(emnl%npx) )
      ierr = L%StereoGraphicInverse( xyz, Radius )
      xyz = xyz/vecnorm(xyz)
      if (ierr.ne.0) then
        masterSPNH(i,j,1) = 0.0
        masterSPSH(i,j,1) = 0.0
      else
        masterSPNH(i,j,1) = InterpolateLambert(xyz, mLPNH, emnl%npx, 1 )
        masterSPSH(i,j,1) = InterpolateLambert(xyz, mLPSH, emnl%npx, 1 )
      end if
    end do
  end do

! combine the original master pattern with this new NBeams pattern as an RGB image
! and store it in a tiff image file; master patter nin red channel, NBeams in green and blue.
! they  are scaled to a common but arbitrary intensity scale 
  RGBmap(1,1:2*npy+1,1:2*npy+1)=EBSDmaster(1:2*npy-1,1:2*npy-1,iE) * &
                                maxval(mLPNH(:,:,1))/maxval(EBSDmaster(:,:,iE))
  RGBmap(2,1:2*npy+1,1:2*npy+1) = masterSPNH(-npy:npy,-npy:npy,1)
  RGBmap(3,1:2*npy+1,1:2*npy+1) = masterSPNH(-npy:npy,-npy:npy,1)
  minbeam(iE) = minval(mLPNH(:,:,1))
  maxbeam(iE) = maxval(mLPNH(:,:,1))

  RGBmap = RGBmap * 255.0 / maxval(RGBmap)

! and generate a color tiff file 
  fname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(emnl%tiffprefix)
  write (number, "(F4.1)") EkeVs(iE)
  TIFF_filename = trim(fname)//'_'//number//'keV.tiff'

! allocate memory for a color image; each pixel has 3 bytes (RGB)
  allocate(TIFF_image(3*(2*npy+1),2*npy+1))
  TIFF_image = reshape( RGBmap, (/ 3*(2*npy+1),2*npy+1 /) )

! set up the image_t structure
  im = image_t(TIFF_image)
  im%dims = (/ 2*npy+1,2*npy+1 /)
  im%samplesPerPixel = 3
  if(im%empty()) call Message%printMessage("EMdpmerge: failed to convert array to rgb image")

! create the file
  call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
  if(0.ne.iostat) then
    call Message%printMessage(" Failed to write image to file : "//iomsg)
  else
    call Message%printMessage(' IPF map written to '//trim(TIFF_filename),"(A)")
  end if

  deallocate(TIFF_image)

! since these computations can take a long time, here we store
! all the output at the end of each pass through the energyloop.

  io_int(1) = nint(float(totstrong)/float(numk))
  call Message%WriteValue(' -> Average number of strong reflections = ',io_int, 1, "(I5)")
  io_int(1) = nint(float(totweak)/float(numk))
  call Message%WriteValue(' -> Average number of weak reflections   = ',io_int, 1, "(I5)")

  call timer%makeTimeStamp()
  dstr = timer%getDateString()
  tstre = timer%getTimeString()

  datagroupname = trim(HDFnames%get_ProgramData())

! open the existing file using the default properties.
  hdferr =  HDF%openFile(outname)

! update the time string
  hdferr = HDF%openGroup(HDFnames%get_EMheader())
  hdferr = HDF%openGroup(datagroupname)

dataset = SC_StopTime
  call timer%Time_tock(1)
  tstop = timer%getInterval(1)
  call timer%Time_reset(1)
  line2(1) = dstr//', '//tstre
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)

  io_int(1) = tstop
  call Message%WriteValue(' Execution time [s]: ',io_int,1)

dataset = SC_Duration
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetFloat(dataset, tstop, overwrite)
  else
    hdferr = HDF%writeDatasetFloat(dataset, tstop)
  end if

  call HDF%pop()
  call HDF%pop()

  hdferr = HDF%openGroup(HDFnames%get_EMData())
  hdferr = HDF%openGroup(datagroupname)

dataset = 'minbeams'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetIntegerArray(dataset, minbeam, numEbins, overwrite)
  else
    hdferr = HDF%writeDatasetIntegerArray(dataset, minbeam, numEbins)
  end if

dataset = 'maxbeams'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetIntegerArray(dataset, maxbeam, numEbins, overwrite)
  else
    hdferr = HDF%writeDatasetIntegerArray(dataset, maxbeam, numEbins)
  end if

! add data to the hyperslab
dataset = SC_mLPNH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1 /)
  offset3 = (/ 0, 0, iE-1 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims3, offset3, cnt3, insert)

dataset = SC_mLPSH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1 /)
  offset3 = (/ 0, 0, iE-1 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims3, offset3, cnt3, insert)

dataset = SC_masterSPNH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1 /)
  offset3 = (/ 0, 0, iE-1 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPNH, dims3, offset3, cnt3, insert)

dataset = SC_masterSPSH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1 /)
  offset3 = (/ 0, 0, iE-1 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPSH, dims3, offset3, cnt3, insert)

  call HDF%popall()

 if (iE.ne.1) then
  call Message%printMessage(' Intermediate data stored in file '//trim(emnl%energyfile), frm = "(A/)")
 end if

 if (iE.eq.1) then
  call Message%printMessage(' Final data stored in file '//trim(emnl%energyfile), frm = "(A/)")
 end if

end do energyloop

call timer%Time_tock(2)
io_int(1) = timer%getInterval(2)
call Message%WriteValue(' Total execution time [s] ',io_int,1)

end associate

call closeFortranHDFInterface()

call Message%printMessage(' min/max number of strong beams vs. energy :')
do i=1,numEbins 
  write (*,*) EkeVs(i), minbeam(i), maxbeam(i)
end do 


call mem%dealloc(EkeVs, 'EkeVs')
call mem%dealloc(mLPNH, 'mLPNH')
call mem%dealloc(mLPSH, 'mLPSH')
call mem%dealloc(masterSPNH, 'masterSPNH')
call mem%dealloc(masterSPSH, 'masterSPSH')

end subroutine EBSDNBeams_

end module mod_EBSDNBeams
