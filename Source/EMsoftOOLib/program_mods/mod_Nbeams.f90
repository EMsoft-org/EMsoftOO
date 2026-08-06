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

module mod_Nbeams
  !! author: MDG
  !! version: 1.0
  !! date: 08/06/26
  !!
  !! class definition for the EMgetNbeams program

use mod_kinds
use mod_global
use mod_MPfiles
use mod_platformsupport

IMPLICIT NONE
private

type, public :: NbeamsNameListType
  integer(kind=irg)   :: npx
  integer(kind=irg)   :: nthreads
  real(kind=sgl)      :: dmin
  real(kind=sgl)      :: EkeV
  character(fnlen)    :: xtalname
  character(fnlen)    :: csvfilename
  character(fnlen)    :: BetheParametersFile
end type NBeamsNameListType

! class definition
type, public :: Nbeams_T
  private
    character(fnlen)          :: nmldeffile = 'EMgetNbeams.nml'
    type(NbeamsNameListType)  :: nml

  contains
  private

    procedure, pass(self) :: setnpx_
    procedure, pass(self) :: getnpx_
    procedure, pass(self) :: setnthreads_
    procedure, pass(self) :: getnthreads_
    procedure, pass(self) :: setdmin_
    procedure, pass(self) :: getdmin_
    procedure, pass(self) :: setEkeV_
    procedure, pass(self) :: getEkeV_
    procedure, pass(self) :: setxtalname_
    procedure, pass(self) :: getxtalname_
    procedure, pass(self) :: setcsvfilename_
    procedure, pass(self) :: getcsvfilename_
    procedure, pass(self) :: setBetheParametersFile_
    procedure, pass(self) :: getBetheParametersFile_
    procedure, pass(self) :: getNameList_
    procedure, pass(self) :: readNameList_
    procedure, pass(self) :: Nbeams_

    generic, public :: setnpx => setnpx_
    generic, public :: getnpx => getnpx_
    generic, public :: setnthreads => setnthreads_
    generic, public :: getnthreads => getnthreads_
    generic, public :: setdmin => setdmin_
    generic, public :: getdmin => getdmin_
    generic, public :: setEkeV => setEkeV_
    generic, public :: getEkeV => getEkeV_
    generic, public :: setxtalname => setxtalname_
    generic, public :: getxtalname => getxtalname_
    generic, public :: setcsvfilename => setcsvfilename_
    generic, public :: getcsvfilename => getcsvfilename_
    generic, public :: setBetheParametersFile => setBetheParametersFile_
    generic, public :: getBetheParametersFile => getBetheParametersFile_
    generic, public :: getNameList => getNameList_
    generic, public :: readNameList => readNameList_
    generic, public :: Nbeams => Nbeams_
 
 end type Nbeams_T

! the constructor routine for this class
interface Nbeams_T
  module procedure Nbeams_constructor
end interface Nbeams_T

contains

!--------------------------------------------------------------------------
type(Nbeams_T) function Nbeams_constructor( nmlfile ) result(Nbeams)
!DEC$ ATTRIBUTES DLLEXPORT :: Nbeams_constructor
!! author: MDG
!! version: 1.0
!! date: 08/06/26
!!
!! constructor for the Nbeams_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call Nbeams%readNameList(nmlfile)

end function Nbeams_constructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 08/06/26
!!
!! read the namelist from an nml file for the Nbeams_T Class

use mod_io

IMPLICIT NONE

class(Nbeams_T), INTENT(INOUT)     :: self
character(fnlen),INTENT(IN)        :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)        :: initonly
 !! fill in the default values only; do not read the file

type(IO_T)                         :: Message
logical                            :: skipread = .FALSE.

integer(kind=irg)   :: npx
integer(kind=irg)   :: nthreads
real(kind=sgl)      :: dmin
real(kind=sgl)      :: EkeV
character(fnlen)    :: xtalname
character(fnlen)    :: csvfilename
character(fnlen)    :: BetheParametersFile


! define the IO namelist to facilitate passing variables to the program.
namelist /getNbeams/ dmin,npx,nthreads,EkeV,xtalname,csvfilename,BetheParametersFile

! set the input parameters to default values (except for xtalname, which must be present)
npx = 500                       ! Nx pixels (total = 2Nx+1)
dmin = 0.05                     ! smallest d-spacing to include in dynamical matrix [nm]
nthreads = 1
xtalname = 'undefined'
EkeV = 30.0
csvfilename = 'undefined'
BetheParametersFile = 'BetheParameters.nml'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=getNbeams)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(xtalname).eq.'undefined') then
  call Message%printError('readNameList:',' xtalname file name is undefined in '//nmlfile)
 end if

 if (trim(csvfilename).eq.'undefined') then
  call Message%printError('readNameList:',' csvfilename file name is undefined in '//nmlfile)
 end if
end if

! if we get here, then all appears to be ok, and we need to fill in the nml fields
self%nml%npx = npx
self%nml%EkeV = EkeV
self%nml%nthreads = nthreads
self%nml%dmin = dmin
self%nml%BetheParametersFile = BetheParametersFile
self%nml%xtalname = xtalname
self%nml%csvfilename = csvfilename

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 08/06/26
!!
!! pass the namelist for the Nbeams_T Class to the calling program

IMPLICIT NONE

class(Nbeams_T), INTENT(INOUT)          :: self
type(NbeamsNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine setnpx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnpx_
!! author: MDG
!! version: 1.0
!! date: 8/6/2026
!!
!! set npx in the Nbeams_T class

IMPLICIT NONE

class(Nbeams_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%npx = inp

end subroutine setnpx_

!--------------------------------------------------------------------------
function getnpx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnpx_
!! author: MDG
!! version: 1.0
!! date: 8/6/2026
!!
!! get npx from the Nbeams_T class

IMPLICIT NONE

class(Nbeams_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%npx

end function getnpx_

!--------------------------------------------------------------------------
subroutine setnthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnthreads_
!! author: MDG
!! version: 1.0
!! date: 8/6/2026
!!
!! set nthreads in the Nbeams_T class

IMPLICIT NONE

class(Nbeams_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%nthreads = inp

end subroutine setnthreads_

!--------------------------------------------------------------------------
function getnthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnthreads_
!! author: MDG
!! version: 1.0
!! date: 8/6/2026
!!
!! get nthreads from the Nbeams_T class

IMPLICIT NONE

class(Nbeams_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%nthreads

end function getnthreads_

!--------------------------------------------------------------------------
subroutine setdmin_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdmin_
!! author: MDG
!! version: 1.0
!! date: 8/6/2026
!!
!! set dmin in the Nbeams_T class

IMPLICIT NONE

class(Nbeams_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%dmin = inp

end subroutine setdmin_

!--------------------------------------------------------------------------
function getdmin_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdmin_
!! author: MDG
!! version: 1.0
!! date: 8/6/2026
!!
!! get dmin from the Nbeams_T class

IMPLICIT NONE

class(Nbeams_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%dmin

end function getdmin_

!--------------------------------------------------------------------------
subroutine setEkeV_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setEkeV_
!! author: MDG
!! version: 1.0
!! date: 8/6/2026
!!
!! set EkeV in the Nbeams_T class

IMPLICIT NONE

class(Nbeams_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%EkeV = inp

end subroutine setEkeV_

!--------------------------------------------------------------------------
function getEkeV_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getEkeV_
!! author: MDG
!! version: 1.0
!! date: 8/6/2026
!!
!! get EkeV from the Nbeams_T class

IMPLICIT NONE

class(Nbeams_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%EkeV

end function getEkeV_

!--------------------------------------------------------------------------
subroutine setxtalname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setxtalname_
!! author: MDG
!! version: 1.0
!! date: 8/6/2026
!!
!! set xtalname in the Nbeams_T class

IMPLICIT NONE

class(Nbeams_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%xtalname = trim(inp)

end subroutine setxtalname_

!--------------------------------------------------------------------------
function getxtalname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getxtalname_
!! author: MDG
!! version: 1.0
!! date: 8/6/2026
!!
!! get xtalname from the Nbeams_T class

IMPLICIT NONE

class(Nbeams_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%xtalname)

end function getxtalname_

!--------------------------------------------------------------------------
subroutine setcsvfilename_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setcsvfilename_
!! author: MDG
!! version: 1.0
!! date: 8/6/2026
!!
!! set csvfilename in the Nbeams_T class

IMPLICIT NONE

class(Nbeams_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%csvfilename = trim(inp)

end subroutine setcsvfilename_

!--------------------------------------------------------------------------
function getcsvfilename_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getcsvfilename_
!! author: MDG
!! version: 1.0
!! date: 8/6/2026
!!
!! get csvfilename from the Nbeams_T class

IMPLICIT NONE

class(Nbeams_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%csvfilename)

end function getcsvfilename_

!--------------------------------------------------------------------------
subroutine setBetheParametersFile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setBetheParametersFile_
!! author: MDG
!! version: 1.0
!! date: 8/6/2026
!!
!! set BetheParametersFile in the Nbeams_T class

IMPLICIT NONE

class(Nbeams_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%BetheParametersFile = trim(inp)

end subroutine setBetheParametersFile_

!--------------------------------------------------------------------------
function getBetheParametersFile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getBetheParametersFile_
!! author: MDG
!! version: 1.0
!! date: 8/6/2026
!!
!! get BetheParametersFile from the Nbeams_T class

IMPLICIT NONE

class(Nbeams_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%BetheParametersFile)

end function getBetheParametersFile_

!--------------------------------------------------------------------------
subroutine Nbeams_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: Nbeams_
!! author: MDG
!! version: 1.0
!! date: 08/06/26
!!
!! compute the number of beams for an EBSD master pattern and output histogram as csv
!!

use mod_EMsoft
use mod_initializers
use mod_symmetry
use mod_crystallography
use mod_gvectors
use mod_kvectors
use mod_io
use HDF5
use mod_HDFsupport
use mod_math
use mod_diffraction
use mod_timing
use mod_memory
use mod_Lambert
use ISO_C_BINDING
use omp_lib
use mod_OMPsupport
use stringconstants

IMPLICIT NONE

class(Nbeams_T), INTENT(INOUT)     :: self
type(EMsoft_T), INTENT(INOUT)      :: EMsoft
character(fnlen),INTENT(IN)        :: progname

type(Cell_T)                       :: cell
type(DynType)                      :: Dyn
type(Timing_T)                     :: timer
type(IO_T)                         :: Message
type(Lambert_T)                    :: L
type(HDF_T)                        :: HDF
type(SpaceGroup_T)                 :: SG
type(Diffraction_T),save           :: Diff
type(kvectors_T)                   :: kvec
type(gvectors_T)                   :: reflist
type(memory_T)                     :: mem, memth

real(kind=dbl)                     :: ctmp(192,3), arg, Radius, xyz(3)
integer(kind=irg)                  :: isym,i,j,ik,npy,ipx,ipy,ipz,debug,iE,izz, izzmax, iequiv(3,48), nequiv, num_el, MCnthreads, & ! counters
                           numk, timestart, timestop, numsites, nthreads, & ! number of independent incident beam directions
                           ir,nat(maxpasym),kk(3), skip, ijmax, one, NUMTHREADS, TID, SamplingType, &
                           numset,n,ix,iy,iz, io_int(6), nns, nnw, nref, Estart, sz(3), ma, &
                           istat,gzero,ic,ip,ikk, totstrong, totweak, jh, ierr, nix, niy, nixp, niyp     ! counters
real(kind=dbl)          :: tpi,Znsq, kkl, DBWF, kin, delta, h, lambda, omtl, srt, dc(3), xy(2), edge, scl, tmp, dx, dxm, dy, dym, &
                           kkk(3), sxy(2) !
real(kind=sgl)          :: io_real(5), selE, kn, FN(3), tstop, nabsl, etotal, density, Ze, at_wt, bp(4)
complex(kind=dbl)               :: czero
logical                         :: usehex, switchmirror, verbose
character(fnlen)                :: xtalname

! Monte Carlo derived quantities
integer(kind=irg)               :: numEbins, nsx, nsy, hdferr, nlines, lastEnergy    ! variables used in MC energy file
character(fnlen)                :: oldprogname, groupname, energyfile, outname, datagroupname, attributename, &
                                   HDF_FileVersion, fname
logical                         :: f_exists, readonly, overwrite=.TRUE., insert=.TRUE., stereog, g_exists, xtaldataread, FL, &
                                   doLegendre, isTKD = .FALSE.

type(gnode),save                :: rlp
real(kind=dbl),allocatable      :: karray(:,:)
integer(kind=irg),allocatable   :: kij(:,:)
integer(kind=irg),allocatable   :: nnsarray(:), histo(:)
type(kvectorlist), pointer      :: ktmp
type(reflisttype), pointer      :: firstw


!$OMP THREADPRIVATE(Diff)

call openFortranHDFInterface()

! set the HDF group names for this program
HDF = HDF_T()

! simplify the notation a little
associate( emnl => self%nml )

tpi = 2.D0*cPi
czero = cmplx(0.D0,0.D0)

!=============================================
!=============================================
! crystallography section;
verbose = .TRUE.

call cell%setFileName(emnl%xtalname)
call Diff%setrlpmethod('WK')

call Diff%setV(dble(emnl%EkeV))
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

write (*,*) '========================'
write (*,*) 'isym = ',isym
write (*,*) 'SamplingType = ', SamplingType
write (*,*) 'usehex = ', usehex
write (*,*) 'SG%trigonal = ', SG%getSpaceGrouptrigonal()
write (*,*) '========================'

! ---------- end of symmetry and crystallography section
!=============================================
!=============================================

!=============================================
!=============================================
! ---------- a couple of initializations
npy = emnl%npx
gzero = 1  ! index of incident beam
!=============================================
!=============================================

! force dynamical matrix routine to read new Bethe parameters from file
! this will all be changed with the new version of the Bethe potentials
call Diff%SetBetheParameters(EMsoft, .FALSE., emnl%BetheParametersFile)

Estart = 1

!=============================================
!=============================================
! ---------- from here on, we need to repeat the entire computation for each energy value
! so this is where we could in principle implement an OpenMP approach; alternatively,
! we could do the inner loop over the incident beam directions in OpenMP (probably simpler)

! we use two times, one (1) for each individual energy level, the other (2) for the overall time
reflist = gvectors_T()

! instantiate the memory class for the OpenMP section
memth = memory_T( nt = emnl%nthreads, silent=.TRUE. )

iE=Estart  ! this used to be the energyloop
selE = emnl%EkeV

! set the accelerating voltage
call cell%setFileName(emnl%xtalname)
call Diff%setV(dble(emnl%EkeV))
call Diff%setrlpmethod('WK')

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
  call Message%printMessage(' Using kvector transform rule for Hall space groups')
  ktmp => kvec%get_ListHead()
  ktmp%k(1:3) = matmul(SG%HallSG%kvec_transform, ktmp%k(1:3))
  do ik=2,numk
    ktmp => ktmp%next
    ktmp%k(1:3) = matmul(SG%HallSG%kvec_transform, ktmp%k(1:3))
  end do
end if

! convert part of the kvector linked list into arrays for OpenMP
mem = memory_T()
call mem%alloc(karray, (/4,numk/), 'karray')
call mem%alloc(nnsarray, (/ numk /), 'nnsarray', initval = 0)
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

write (*,*) ' entering parallel section '

! use OpenMP to run on multiple cores ...
!$OMP PARALLEL COPYIN(Diff) &
!$OMP& PRIVATE(ik,TID,kn,ipx,ipy,ipz,ix,iequiv,nequiv,reflist,firstw) &
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

! determine strong and weak reflections
     nullify(firstw)
     nns = 0
     nnw = 0
     call reflist%Apply_BethePotentials(Diff, firstw, nns, nnw)

! and store the number of strong beams, which is the only program output
     nnsarray(ik) = nns

     if (mod(ik,5000).eq.0) then
       io_int(1) = ik
       io_int(2) = numk
       call Message%WriteValue('  completed beam direction ',io_int, 2, "(I8,' of ',I8)")
     end if

     call reflist%Delete_gvectorlist()

    end do beamloop

! end of OpenMP portion
!$OMP END PARALLEL

! reorganize the nnsarray into a histogram and store it in a csv file
ma = maxval(nnsarray)
call mem%alloc(histo, (/ ma /), 'histo', initval = 0)
do ik=1,numk
  histo(nnsarray(ik)) = histo(nnsarray(ik)) + 1
end do 

fname = EMsoft%generateFilePath('EMdatapathname', emnl%csvfilename)
open(dataunit, file=trim(fname), status='unknown', form='formatted')
write(dataunit, "(A)") 'Beams, Frequency'
do ik=1,ma
  write (dataunit,"(I5,',',I8)") ik, histo(ik)
end do 
close(dataunit,status='keep')
call Message%printMessage(' data stored in file '//trim(fname))


end associate

call closeFortranHDFInterface()

end subroutine Nbeams_

end module mod_Nbeams