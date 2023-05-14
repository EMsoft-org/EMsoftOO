! ###################################################################
! Copyright (c) 2013-2022, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_ECPmaster
  !! author: MDG
  !! version: 1.0
  !! date: 03/03/20
  !!
  !! class definition for the EMECPmaster program

use mod_kinds
use mod_global
use mod_MPfiles

IMPLICIT NONE

! class definition
type, public :: ECPmaster_T
private
  character(fnlen)             :: nmldeffile = 'EMECPmaster.nml'
  type(ECPmasterNameListType)  :: nml

contains
private
  procedure, pass(self) :: readNameList_
!  procedure, pass(self) :: writeHDFNameList_ replaced by routine in MPfiles.f90
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: ECPmaster_

  generic, public :: getNameList => getNameList_
  ! generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: ECPmaster => ECPmaster_

end type ECPmaster_T

! the constructor routine for this class
interface ECPmaster_T
  module procedure ECPmaster_constructor
end interface ECPmaster_T

contains

!--------------------------------------------------------------------------
type(ECPmaster_T) function ECPmaster_constructor( nmlfile ) result(ECPmaster)
!DEC$ ATTRIBUTES DLLEXPORT :: ECPmaster_constructor
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! constructor for the ECPmaster_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call ECPmaster%readNameList(nmlfile)

end function ECPmaster_constructor

!--------------------------------------------------------------------------
subroutine ECPmaster_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: ECPmaster_destructor
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! destructor for the ECPmaster_T Class

IMPLICIT NONE

type(ECPmaster_T), INTENT(INOUT)  :: self

call reportDestructor('ECPmaster_T')

end subroutine ECPmaster_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! read the namelist from an nml file for the ECPmaster_T Class

use mod_io
use mod_EMsoft

IMPLICIT NONE

class(ECPmaster_T), INTENT(INOUT)           :: self
character(fnlen),INTENT(IN)                 :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)                 :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                              :: EMsoft
type(IO_T)                                  :: Message
logical                                     :: skipread = .FALSE.

integer(kind=irg)                           :: npx
integer(kind=irg)                           :: nthreads
real(kind=sgl)                              :: dmin
character(3)                                :: Notify
character(fnlen)                            :: copyfromenergyfile
character(fnlen)                            :: h5copypath
character(fnlen)                            :: BetheParametersFile
character(fnlen)                            :: energyfile
logical                                     :: combinesites
logical                                     :: kinematical

! define the IO namelist to facilitate passing variables to the program.
namelist /ECPmastervars/ dmin, Notify, h5copypath, energyfile, npx, nthreads, copyfromenergyfile, combinesites, & 
                         kinematical, BetheParametersFile

! set the input parameters to default values (except for xtalname, which must be present)
nthreads = 1
dmin = 0.04                     ! smallest d-spacing to include in dynamical matrix [nm]
npx = 256
Notify = 'Off'
h5copypath = 'undefined'
energyfile = 'undefined'        ! default filename for z_0(E_e) data from EMMC Monte Carlo simulations
copyfromenergyfile = 'undefined'
BetheParametersFile = 'undefined'
kinematical = .FALSE.
combinesites = .FALSE.

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
  open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
  read(UNIT=dataunit,NML=ECPmastervars)
  close(UNIT=dataunit,STATUS='keep')

! check for required entries
  if (trim(energyfile).eq.'undefined') then
    call Message%printError('EMECPmaster:',' energy file name is undefined in '//nmlfile)
  end if
end if

! if we get here, then all appears to be ok, and we need to fill in the emnl fields
self%nml%npx = npx
self%nml%nthreads = nthreads
self%nml%dmin = dmin
self%nml%Notify = Notify
self%nml%h5copypath = h5copypath
self%nml%copyfromenergyfile = copyfromenergyfile
self%nml%BetheParametersFile = BetheParametersFile
self%nml%energyfile = energyfile
self%nml%combinesites = combinesites
self%nml%kinematical = kinematical

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! pass the namelist for the ECPmaster_T Class to the calling program

IMPLICIT NONE

class(ECPmaster_T), INTENT(INOUT)          :: self
type(ECPmasterNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine ECPmaster_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: ECPmaster_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! compute an ECP master pattern for a given energy

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
use mod_notifications
use stringconstants
use mod_MCfiles

IMPLICIT NONE

class(ECPmaster_T), INTENT(INOUT)       :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen),INTENT(IN)             :: progname

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
type(HDFnames_T)        :: HDFnames

real(kind=dbl)          :: frac
integer(kind=irg)       :: gzero, istat, tickstart
type(MCOpenCLNameListType) :: mcnl

integer(kind=irg)       :: numangle, numzbins, nx, ny, npy, totnum_el, numsites ! reading from MC file
real(kind=dbl)          :: EkeV, Ehistmin, Ebinsize, depthmax, depthstep, sig, omega  ! reading from MC file
integer(kind=irg), allocatable :: acc_z(:,:,:,:),accum_z(:,:,:,:) ! reading from MC file

integer(kind=irg)       :: io_int_sgl(1), io_int(6) ! integer output variable
real(kind=dbl)          :: io_real(5) ! real output variable

integer(kind=irg)       :: i, j, ik, kkk, isym, pgnum, SamplingType, nix, nixp, niy, niyp, hkl(3) ! variables for point group and Laue group
integer(kind=irg),parameter     :: LaueTest(11) = (/ 149, 151, 153, 156, 158, 160, 161, 164, 165, 166, 167 /)  ! space groups with 2 or mirror at 30 degrees
integer(kind=irg)       :: npyhex, ijmax, numk, skip ! parameters for calckvectors and calcwavelength subroutine

integer(kind=irg)       :: ga(3), gb(3) ! shortest reciprocal lattice vector for zone axis
real(kind=sgl), allocatable :: thick(:), mLPNH(:,:,:), mLPSH(:,:,:), svals(:), lambdaZ(:), klist(:,:), knlist(:),&
                               masterSPNH(:,:,:), masterSPSH(:,:,:), auxNH(:,:,:), auxSH(:,:,:)
real(kind=dbl)          :: intthick, dc(3), dx, dxm, dy, dym, edge, scl, xy(2), Radius, tpi
complex(kind=dbl),allocatable   :: Lgh(:,:),Sgh(:,:),Sghtmp(:,:,:)
complex(kind=dbl),allocatable   :: DynMat(:,:)
complex(kind=dbl)       :: czero

integer(kind=irg)       :: nt, nns, nnw, tots, totw ! thickness array and BetheParameters strong and weak beams
real(kind=sgl)          :: FN(3), kk(3), fnat, kn, tstop
integer(kind=irg)       :: numset, nref, ipx, ipy, ipz, iequiv(3,48), nequiv, ip, jp, izz, IE, iz, one,ierr, &
                           NUMTHREADS, totstrong, totweak, nat(maxpasym)
integer(kind=irg),allocatable   :: kij(:,:)
real(kind=dbl)          :: res(2), xyz(3), ind, nabsl

character(fnlen)        :: oldprogname, energyfile, outname
character(fnlen)        :: xtalname, groupname, datagroupname, HDF_FileVersion, attributename
character(8)            :: MCscversion
character(4)            :: MCmode
character(6)            :: projtype
character(11)           :: dstr
character(15)           :: tstrb
character(15)           :: tstre
logical                 :: f_exists, readonly, overwrite=.TRUE., insert=.TRUE., g_exists, xtaldataread, stereog
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)
character(fnlen,kind=c_char)                     :: line2(1)


logical                             :: verbose, usehex, switchmirror

type(gnode),save                    :: rlp
type(kvectorlist), pointer          :: ktmp ! linked list for incident wave vectors for master list
type(kvectorlist), pointer          :: kheadcone,ktmpcone ! linked list for incident wave vectors for individual pattern
real(kind=dbl),allocatable          :: ecpattern(:,:)
type(reflisttype),pointer           :: firstw,rltmp
integer(kind=irg)                   :: nthreads,TID,ix,hdferr,num_el,etotal, nlines,nsx,nsy,SelE
character(fnlen)                    :: dataset, instring, fname
character(fnlen)                    :: mode
integer(HSIZE_T)                    :: dims4(4), cnt4(4), offset4(4), dims3(3), cnt3(3), offset3(3)

character(fnlen),ALLOCATABLE        :: MessageLines(:)
integer(kind=irg)                   :: NumLines
character(fnlen)                    :: SlackUsername, exectime
character(100)                      :: c

!$OMP THREADPRIVATE(Diff)

call openFortranHDFInterface()
HDF = HDF_T()

! set the HDF group names for this program
HDFnames = HDFnames_T()

call MPFT%setModality('ECP')

! simplify the notation a little
associate( emnl => self%nml )

! initialize the timing routines
timer = Timing_T()
tstrb = timer%getTimeString()

! if copyfromenergyfile is different from 'undefined', then we need to
! copy all the Monte Carlo data from that file into a new file, which
! will then be read from and written to by the ComputeMasterPattern routine.
if (emnl%copyfromenergyfile.ne.'undefined') then
  call MCFT%copyMCdata(EMsoft, HDF, emnl%copyfromenergyfile, emnl%energyfile, emnl%h5copypath)
end if

stereog = .TRUE.

timer = Timing_T()
tstrb = timer%getTimeString()

tpi = 2.D0*cPi
czero = cmplx(0.D0,0.D0)
gzero = 1
frac = 0.05

!=============================================
!=============================================
! ---------- read Monte Carlo .h5 output file and extract necessary parameters
call HDFnames%set_ProgramData(SC_MCOpenCL)
call HDFnames%set_NMLlist(SC_MCCLNameList)
call HDFnames%set_NMLfilename(SC_MCOpenCLNML)
energyfile = ''
energyfile = EMsoft%generateFilePath('EMdatapathname',trim(emnl%energyfile))
outname = trim(energyfile)
call MCFT%setFileName(energyfile)
call MCFT%readMCfile(HDF, HDFnames, getAccumz=.TRUE.)
mcnl = MCFT%getnml()
call MCFT%copyaccumz(accum_z)

numzbins = MCFT%getnumzbins()
numangle = MCFT%getnumangles()

nsx = (mcnl%numsx - 1)/2
nsy = nsx
etotal = sum(accum_z(numangle,:,:,:))

io_int(1) = mcnl%totnum_el
call Message%WriteValue(' --> total number of BSE electrons in MC data set ', io_int, 1)

!=============================================
!=============================================
! crystallography section
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
    ' The ECP master pattern simulation for rhombohedral/trigonal structures  ', &
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
   ijmax = float(emnl%npx)**2   ! truncation value for beam directions
! ----------
!=============================================
!=============================================

!=============================================
!=============================================
! if the combinesites parameter is .TRUE., then we only need to
! allocate a dimension of 1 in the master pattern array since we are adding
! together the master patterns for all sites in the asymmetric unit.
  if (emnl%combinesites.eqv..TRUE.) then
    numsites = 1
  else
    numsites = numset
  end if

! ---------- allocate memory for the master patterns
  allocate(mLPNH(-emnl%npx:emnl%npx,-npy:npy,1:numsites),stat=istat)
  allocate(mLPSH(-emnl%npx:emnl%npx,-npy:npy,1:numsites),stat=istat)
  allocate(masterSPNH(-emnl%npx:emnl%npx,-npy:npy,1))
  allocate(masterSPSH(-emnl%npx:emnl%npx,-npy:npy,1))

! set various arrays to zero
  if (emnl%uniform.eqv..TRUE.) then
   mLPNH = 1.0
   mLPSH = 1.0
   masterSPNH = 1.0
   masterSPSH = 1.0
  else
   mLPNH = 0.0
   mLPSH = 0.0
   masterSPNH = 0.0
   masterSPSH = 0.0
  end if

! force dynamical matrix routine to read new Bethe parameters from file
  call Diff%SetBetheParameters(EMsoft, .FALSE.) ! , emnl%BetheParametersFile)

!=============================================
! create or update the HDF5 output file
!=============================================
call HDFnames%set_ProgramData(SC_ECPmaster)
call HDFnames%set_NMLlist(SC_ECPmasterNameList)
call HDFnames%set_NMLfilename(SC_ECPmasterNML)
! Open an existing file or create a new file using the default properties.
  if (trim(energyfile).eq.trim(outname)) then
    hdferr =  HDF%openFile(outname)
  else
    hdferr =  HDF%createFile(outname)
  end if

! write the EMheader to the file
  datagroupname = trim(HDFnames%get_ProgramData())
  call HDF%writeEMheader(EMsoft,dstr, tstrb, tstre, progname, datagroupname)

! open or create a namelist group to write all the namelist files into
  hdferr = HDF%createGroup(HDFnames%get_NMLfiles())

! read the text file and write the array to the file
  dataset = HDFnames%get_NMLfilename()
  hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%nmldeffile)

! leave this group
  call HDF%pop()

! create a namelist group to write all the namelist files into
  hdferr = HDF%createGroup(HDFnames%get_NMLparameters())

  call MPFT%writeHDFNameList(HDF, HDFnames, emnl)
  call Diff%writeBetheparameterNameList(HDF)

! leave this group
  call HDF%pop()

! then the remainder of the data in a EMData group
  hdferr = HDF%createGroup(HDFnames%get_EMData())

! create the EBSDmaster group and add a HDF_FileVersion attribbute to it
  hdferr = HDF%createGroup(HDFnames%get_ProgramData())
  HDF_FileVersion = '4.0'
  HDF_FileVersion = cstringify(HDF_FileVersion)
  attributename = SC_HDFFileVersion
  hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)

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

dataset = SC_numset
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetInteger(dataset, numset, overwrite)
  else
    hdferr = HDF%writeDatasetInteger(dataset, numset)
  end if

dataset = SC_EkeV
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetDouble(dataset, EkeV, overwrite)
  else
    hdferr = HDF%writeDatasetDouble(dataset, EkeV)
  end if

! create the hyperslabs and write zeroes to them for now
dataset = SC_mLPNH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, numsites /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, numsites /)
  offset3 = (/ 0, 0, 0 /)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims3, offset3, cnt3, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims3, offset3, cnt3)
  end if

dataset = SC_mLPSH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, numsites /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, numsites /)
  offset3 = (/ 0, 0, 0 /)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims3, offset3, cnt3, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims3, offset3, cnt3)
  end if

dataset = SC_masterSPNH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, 1 /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1 /)
  offset3 = (/ 0, 0, 0 /)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPNH, dims3, offset3, cnt3, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPNH, dims3, offset3, cnt3)
  end if

dataset = SC_masterSPSH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, 1 /)
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

! we use two times, one (1) for each individual energy level, the other (2) for the overall time
call timer%Time_tick(1)
reflist = gvectors_T()

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

! ! convert part of the kvector linked list into arrays for OpenMP
   izz = MCFT%getnumzbins()

   allocate(lambdaZ(1:izz),stat=istat)
   allocate(kij(3,numk),stat=istat)
   allocate(klist(3,numk),knlist(numk),stat=istat)

! ! point to the first beam direction
   ktmp => kvec%get_ListHead()
! ! and loop through the list, keeping k, kn, and i,j
   kij(1:3,1) = (/ ktmp%i, ktmp%j, ktmp%hs /)
   klist(1:3,1) = ktmp%k
   knlist(1) = ktmp%kn

   do i = 2,numk
       ktmp => ktmp%next
       kij(1:3,i) = (/ ktmp%i, ktmp%j, ktmp%hs /)
       klist(1:3,i) = ktmp%k
       knlist(i) = ktmp%kn
   end do

! ! and remove the linked list
   call kvec%Delete_kvectorlist()

   call Diff%setrlpmethod('WK')
  hkl=(/0,0,0/)
  call Diff%CalcUcg(cell,hkl)
  rlp = Diff%getrlp()
  nabsl = rlp%xgp

  do iz=1,izz
    lambdaZ(iz) = float(sum(accum_z(numangle,iz,:,:)))/float(etotal)
    lambdaZ(iz) = lambdaZ(iz) * exp(2.0*sngl(cPi)*(iz-1)*mcnl%depthstep/nabsl)
  end do

verbose = .FALSE.
totstrong = 0
totweak = 0

fnat = 1.0/float(sum(cell%getnumat()))
intthick = dble(mcnl%depthmax)

! here's where we introduce the OpenMP calls, to speed up the overall calculations...

! set the number of OpenMP threads
nthreads = emnl%nthreads
!call OMP_SET_NUM_THREADS(nthreads)
call OMP_setNThreads(nthreads)

! use OpenMP to run on multiple cores ...
!$OMP PARALLEL COPYIN(Diff) &
!$OMP& PRIVATE(DynMat,Sgh,sghtmp,Lgh,i,FN,TID,kn,ipx,ipy,ix,ip,iequiv,nequiv,reflist,firstw) &
!$OMP& PRIVATE(kk,nns,nnw,nref,nat,io_int,io_int_sgl,nthreads,svals)

NUMTHREADS = OMP_GET_NUM_THREADS()
TID = OMP_GET_THREAD_NUM()

allocate(svals(numset))

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
     kk = klist(1:3,ik)
     FN = kk

     ! 02-06-2020 CLément Lafond : using emnl% cause fortran error 157 on Windows, replaced by self%nml% while trying to fix it
     call reflist%Initialize_ReflectionList(cell, SG, Diff, FN, kk, self%nml%dmin, verbose)
     nref = reflist%get_nref()
! ---------- end of "create the master reflection list"
!=============================================

! determine strong and weak reflections
     nullify(firstw)
     nns = 0
     nnw = 0
     call reflist%Apply_BethePotentials(Diff, firstw, nns, nnw)

     if (self%nml%kinematical.eqv..FALSE.) then
! generate the dynamical matrix
       if (allocated(DynMat)) deallocate(DynMat)
       allocate(DynMat(nns,nns))
       call reflist%GetDynMat(cell, Diff, firstw, DynMat, nns, nnw)
       totstrong = totstrong + nns
       totweak = totweak + nnw
       if (ik.eq.1) write (*,*) ' maxval(DynMat) = ', maxval(abs(DynMat))
     else
! all reflections are strong, but they are not coupled to each other, only to the
! incident beam; all q_{g-g'} are zero except the ones with g'=0.  In addition, there
! is no anomalous absorption, only normal absorption.
       if (allocated(DynMat)) deallocate(DynMat)
       allocate(DynMat(nns,nns))
       call reflist%GetDynMatKin(cell, Diff, firstw, DynMat, nns)
       totstrong = totstrong + nns
       totweak = 0
     end if

! then we need to initialize the Sgh and Lgh arrays
     if (allocated(Sghtmp)) deallocate(Sghtmp)
     if (allocated(Lgh)) deallocate(Lgh)
     allocate(Sghtmp(nns,nns,numset),Lgh(nns,nns))
     Sghtmp = czero
     Lgh = czero
     call reflist%getSghfromLUT(Diff,nns,numset,Sghtmp)

! solve the dynamical eigenvalue equation for this beam direction
     kn = knlist(ik)
     call reflist%CalcLgh(DynMat,Lgh,intthick,dble(kn),nns,gzero,mcnl%depthstep,lambdaZ,izz)

! sum over the element-wise (Hadamard) product of the Lgh and Sgh arrays
     svals = 0.0
     do ix=1,numset
       svals(ix) = real(sum(Lgh(1:nns,1:nns)*Sghtmp(1:nns,1:nns,ix)))
     end do
     svals = svals * fnat

! and store the resulting svals values, applying point group symmetry where needed.
     ipx = kij(1,ik)
     ipy = kij(2,ik)
     ipz = kij(3,ik)
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
  if (self%nml%combinesites.eqv..FALSE.) then
     do ix=1,nequiv
       if (iequiv(3,ix).eq.-1) mLPSH(iequiv(1,ix),iequiv(2,ix),1:numset) = svals(1:numset)
       if (iequiv(3,ix).eq.1) mLPNH(iequiv(1,ix),iequiv(2,ix),1:numset) = svals(1:numset)
     end do
  else
     do ix=1,nequiv
       if (iequiv(3,ix).eq.-1) mLPSH(iequiv(1,ix),iequiv(2,ix),1) = sum(svals)
       if (iequiv(3,ix).eq.1) mLPNH(iequiv(1,ix),iequiv(2,ix),1) = sum(svals)
     end do
  end if
!$OMP END CRITICAL

    totw = totw + nnw
    tots = tots + nns

    deallocate(Lgh, Sghtmp)

    if (mod(ik,5000).eq.0) then
       io_int(1) = ik
       io_int(2) = numk
       call Message%WriteValue('  completed beam direction ',io_int, 2, "(I8,' of ',I8)")
    end if

     call reflist%Delete_gvectorlist()

  end do beamloop

!end of OpenMP portion
!$OMP END PARALLEL

io_int(1) = nint(float(tots)/float(numk))
call Message%WriteValue(' -> Average number of strong reflections = ',io_int, 1, "(I5)")
io_int(1) = nint(float(totw)/float(numk))
call Message%WriteValue(' -> Average number of weak reflections   = ',io_int, 1, "(I5)")

if (usehex) then
! and finally, we convert the hexagonally sampled array to a square Lambert projection which will be used
! for all ECP pattern interpolations;  we need to do this for both the Northern and Southern hemispheres

! we begin by allocating auxiliary arrays to hold copies of the hexagonal data; the original arrays will
! then be overwritten with the newly interpolated data.
    allocate(auxNH(-emnl%npx:emnl%npx,-npy:npy,1:numset),stat=istat)
    allocate(auxSH(-emnl%npx:emnl%npx,-npy:npy,1:numset),stat=istat)
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
          mLPNH(i,j,1:numset) = auxNH(nix,niy,1:numset)*dxm*dym + auxNH(nixp,niy,1:numset)*dx*dym + &
                                auxNH(nix,niyp,1:numset)*dxm*dy + auxNH(nixp,niyp,1:numset)*dx*dy
          mLPSH(i,j,1:numset) = auxSH(nix,niy,1:numset)*dxm*dym + auxSH(nixp,niy,1:numset)*dx*dym + &
                                auxSH(nix,niyp,1:numset)*dxm*dy + auxSH(nixp,niyp,1:numset)*dx*dy
        end if
      end do
    end do
    deallocate(auxNH, auxSH)
end if

! make sure that the outer pixel rim of the mLPSH patterns is identical to
! that of the mLPNH array.
mLPSH(-emnl%npx,-emnl%npx:emnl%npx,1:numset) = mLPNH(-emnl%npx,-emnl%npx:emnl%npx,1:numset)
mLPSH( emnl%npx,-emnl%npx:emnl%npx,1:numset) = mLPNH( emnl%npx,-emnl%npx:emnl%npx,1:numset)
mLPSH(-emnl%npx:emnl%npx,-emnl%npx,1:numset) = mLPNH(-emnl%npx:emnl%npx,-emnl%npx,1:numset)
mLPSH(-emnl%npx:emnl%npx, emnl%npx,1:numset) = mLPNH(-emnl%npx:emnl%npx, emnl%npx,1:numset)

! get stereographic projections (summed over the atomic positions)
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
      masterSPNH(i,j,1) = InterpolateLambert(xyz, mLPNH, emnl%npx, numset)
      masterSPSH(i,j,1) = InterpolateLambert(xyz, mLPSH, emnl%npx, numset)
    end if
  end do
end do

! open the existing file using the default properties.
  hdferr =  HDF%openFile(outname)

! update the time string
  hdferr = HDF%openGroup(HDFnames%get_EMheader())
  hdferr = HDF%openGroup(HDFnames%get_ProgramData())

dataset = SC_StopTime
  call timer%Time_tock(1)
  tstop = timer%getInterval(1)
  call timer%Time_reset(1)

  call timer%makeTimeStamp()
  dstr = timer%getDateString()
  tstre = timer%getTimeString()

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
  hdferr = HDF%openGroup(HDFnames%get_ProgramData())

! add data to the hyperslab
dataset = SC_mLPNH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, numsites /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, numsites /)
  offset3 = (/ 0, 0, 0 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims3, offset3, cnt3, insert)

dataset = SC_mLPSH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, numsites /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, numsites /)
  offset3 = (/ 0, 0, 0 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims3, offset3, cnt3, insert)

dataset = SC_masterSPNH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, numsites /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, numsites /)
  offset3 = (/ 0, 0, 0 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPNH, dims3, offset3, cnt3, insert)

dataset = SC_masterSPSH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, numsites /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, numsites /)
  offset3 = (/ 0, 0, 0 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPSH, dims3, offset3, cnt3, insert)

call HDF%popall()

call Message%printMessage(' Final data stored in file '//trim(emnl%energyfile), frm = "(A/)")

! if requested, we notify the user that this program has completed its run
if (trim(EMsoft%getConfigParameter('EMNotify')).ne.'Off') then
  if (trim(emnl%Notify).eq.'On') then
    NumLines = 3
    allocate(MessageLines(NumLines))

    call hostnm(c)

    MessageLines(1) = 'EMECPmaster program has ended successfully'
    MessageLines(2) = 'Master pattern data stored in '//trim(outname)
    write (exectime,"(F10.4)") timer%getInterval(2)
    MessageLines(3) = 'Total execution time [s]: '//trim(exectime)
    SlackUsername = 'EMsoft on '//trim(c)
    i = PostMessage(EMsoft, MessageLines, NumLines, SlackUsername)
  end if
end if

end associate

end subroutine ECPmaster_



end module mod_ECPmaster
