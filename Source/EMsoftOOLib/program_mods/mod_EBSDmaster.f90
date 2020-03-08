! ###################################################################
! Copyright (c) 2013-2020, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_EBSDmaster
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/05/20
  !!
  !! class definition for the EMEBSDmaster program

use mod_kinds
use mod_global
use mod_MPfiles

IMPLICIT NONE 
private 


! class definition
type, public :: EBSDmaster_T
  private 
    character(fnlen)              :: nmldeffile = 'EMEBSDmaster.nml'
    type(EBSDmasterNameListType)  :: nml 

  contains
  private 

    procedure, pass(self) :: readNameList_
    procedure, pass(self) :: getNameList_
    procedure, pass(self) :: EBSDmaster_

    generic, public :: getNameList => getNameList_
    generic, public :: readNameList => readNameList_
    generic, public :: EBSDmaster => EBSDmaster_

end type EBSDmaster_T

! the constructor routine for this class 
interface EBSDmaster_T
  module procedure EBSDmaster_constructor
end interface EBSDmaster_T

contains

!--------------------------------------------------------------------------
type(EBSDmaster_T) function EBSDmaster_constructor( nmlfile ) result(EBSDmaster)
!! author: MDG 
!! version: 1.0 
!! date: 02/05/20
!!
!! constructor for the EBSDmaster_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call EBSDmaster%readNameList(nmlfile)

end function EBSDmaster_constructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!! author: MDG 
!! version: 1.0 
!! date: 02/05/20
!!
!! read the namelist from an nml file for the EBSDmaster_T Class 

use mod_io 

IMPLICIT NONE 

class(EBSDmaster_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)                 :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)                 :: initonly
 !! fill in the default values only; do not read the file

type(IO_T)                                  :: Message       
logical                                     :: skipread = .FALSE.

integer(kind=irg) :: npx
integer(kind=irg) :: Esel
integer(kind=irg) :: nthreads
real(kind=sgl)    :: dmin
character(3)      :: Notify
character(fnlen)  :: copyfromenergyfile
character(fnlen)  :: energyfile
character(fnlen)  :: BetheParametersFile
character(fnlen)  :: h5copypath
logical           :: combinesites
logical           :: restart
logical           :: uniform
logical           :: kinematical

! define the IO namelist to facilitate passing variables to the program.
namelist /EBSDmastervars/ dmin,npx,nthreads,copyfromenergyfile,energyfile,Esel,restart,uniform,Notify, &
                          combinesites,h5copypath,BetheParametersFile,kinematical

! set the input parameters to default values (except for xtalname, which must be present)
npx = 500                       ! Nx pixels (total = 2Nx+1)
nthreads = 1
Esel = -1                       ! selected energy value for single energy run
dmin = 0.025                    ! smallest d-spacing to include in dynamical matrix [nm]
Notify = 'Off'
copyfromenergyfile = 'undefined'! default filename for z_0(E_e) data from a different Monte Carlo simulation
h5copypath = 'undefined'
energyfile = 'undefined'        ! default filename for z_0(E_e) data from EMMC Monte Carlo simulations
BetheParametersFile='BetheParameters.nml'
combinesites = .FALSE.          ! combine all atom sites into one BSE yield or not
restart = .FALSE.               ! when .TRUE. an existing file will be assumed 
uniform = .FALSE.               ! when .TRUE., the output master patterns will contain 1.0 everywhere
kinematical = .FALSE.           ! use the kinematical approximation if .TRUE.

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=EBSDmastervars)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(energyfile).eq.'undefined') then
  call Message%printError('readNameList:',' output (energy) file name is undefined in '//nmlfile)
 end if
end if

! if we get here, then all appears to be ok, and we need to fill in the nml fields
self%nml%npx = npx
self%nml%Esel = Esel
self%nml%nthreads = nthreads
self%nml%dmin = dmin
self%nml%copyfromenergyfile = copyfromenergyfile
self%nml%h5copypath = h5copypath
self%nml%energyfile = energyfile
self%nml%BetheParametersFile = BetheParametersFile
self%nml%Notify = Notify
self%nml%combinesites = combinesites
self%nml%restart = restart
self%nml%uniform = uniform
self%nml%kinematical = kinematical

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!! author: MDG 
!! version: 1.0 
!! date: 02/05/20
!!
!! pass the namelist for the EBSDmaster_T Class to the calling program

IMPLICIT NONE 

class(EBSDmaster_T), INTENT(INOUT)          :: self
type(EBSDmasterNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine EBSDmaster_(self, EMsoft, progname)
!! author: MDG 
!! version: 1.0 
!! date: 02/05/20
!!
!! compute an EBSD master pattern as a function of energy

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
use mod_notifications
use stringconstants
use mod_MCfiles

IMPLICIT NONE 

class(EBSDmaster_T), INTENT(INOUT)  :: self
type(EMsoft_T), INTENT(INOUT)       :: EMsoft
character(fnlen),INTENT(IN)         :: progname

type(Cell_T)            :: cell
type(DynType)           :: Dyn
type(Timing_T)          :: timer
type(IO_T)              :: Message
type(Lambert_T)         :: L
type(HDF_T)             :: HDF
type(SpaceGroup_T)      :: SG
type(Diffraction_T)     :: Diff
type(MCfile_T)          :: MCFT
type(MPfile_T)          :: MPFT
type(kvectors_T)        :: kvec
type(gvectors_T)        :: reflist
type(HDFnames_T)        :: HDFnames

real(kind=dbl)          :: ctmp(192,3), arg, Radius, xyz(3)
integer(HSIZE_T)        :: dims4(4), cnt4(4), offset4(4)
integer(HSIZE_T)        :: dims3(3), cnt3(3), offset3(3)
integer(kind=irg)       :: isym,i,j,ik,npy,ipx,ipy,ipz,debug,iE,izz, izzmax, iequiv(3,48), nequiv, num_el, MCnthreads, & ! counters
                           numk, timestart, timestop, numsites, nthreads, & ! number of independent incident beam directions
                           ir,nat(100),kk(3), skip, ijmax, one, NUMTHREADS, TID, SamplingType, &
                           numset,n,ix,iy,iz, io_int(6), nns, nnw, nref, Estart, &
                           istat,gzero,ic,ip,ikk, totstrong, totweak, jh, ierr, nix, niy, nixp, niyp     ! counters
real(kind=dbl)          :: tpi,Znsq, kkl, DBWF, kin, delta, h, lambda, omtl, srt, dc(3), xy(2), edge, scl, tmp, dx, dxm, dy, dym !
real(kind=sgl)          :: io_real(5), selE, kn, FN(3), kkk(3), tstop, bp(4), nabsl, etotal, density, Ze, at_wt
real(kind=sgl),allocatable      :: EkeVs(:), svals(:), auxNH(:,:,:,:), auxSH(:,:,:,:), Z2percent(:)  ! results
real(kind=sgl),allocatable      :: mLPNH(:,:,:,:), mLPSH(:,:,:,:), masterSPNH(:,:,:), masterSPSH(:,:,:)
real(kind=dbl),allocatable      :: LegendreArray(:), upd(:), diagonal(:)
integer(kind=irg),allocatable   :: accum_z(:,:,:,:)
complex(kind=dbl)               :: czero
complex(kind=dbl),allocatable   :: Lgh(:,:), Sgh(:,:,:)
logical                 :: usehex, switchmirror, verbose
character(fnlen)        :: xtalname

! Monte Carlo derived quantities
integer(kind=irg)       :: numEbins, nsx, nsy, hdferr, nlines, lastEnergy    ! variables used in MC energy file
integer(kind=irg),allocatable :: thick(:)
real(kind=sgl),allocatable :: lambdaE(:,:)
character(fnlen)        :: oldprogname, groupname, energyfile, outname, datagroupname, attributename, HDF_FileVersion, fname
character(8)            :: MCscversion
character(11)           :: dstr
character(15)           :: tstrb
character(15)           :: tstre
logical                 :: f_exists, readonly, overwrite=.TRUE., insert=.TRUE., stereog, g_exists, xtaldataread, FL, doLegendre
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)
character(fnlen,kind=c_char)                     :: line2(1)

type(gnode),save                :: rlp
real(kind=sgl),allocatable      :: karray(:,:)
integer(kind=irg),allocatable   :: kij(:,:)
complex(kind=dbl),allocatable   :: DynMat(:,:)
character(fnlen)                :: dataset, instring
type(MCOpenCLNameListType)      :: mcnl
type(kvectorlist), pointer      :: ktmp
type(reflisttype), pointer      :: firstw

character(fnlen),ALLOCATABLE    :: MessageLines(:)
integer(kind=irg)               :: NumLines, info
character(fnlen)                :: SlackUsername, exectime
character(100)                  :: c

!$OMP THREADPRIVATE(rlp) 

call openFortranHDFInterface()
HDF = HDF_T() 

! set the HDF group names for this program
HDFnames = HDFnames_T() 
call HDFnames%set_ProgramData(SC_EBSDmaster) 
call HDFnames%set_NMLlist(SC_EBSDmasterNameList) 
call HDFnames%set_NMLfilename(SC_EBSDmasterNML) 
call HDFnames%set_Variable(SC_MCOpenCL) 
call MPFT%setModality('EBSD')

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

! is the master pattern used for spherical indexing only ?  If so, then we need to modifiy the k-vector sampling
doLegendre = .FALSE.

!=============================================
!=============================================
! ---------- read Monte Carlo .h5 output file and extract necessary parameters
fname = EMsoft%generateFilePath('EMdatapathname',trim(emnl%energyfile))
call MCFT%setFileName(fname)
call MCFT%readMCfile(HDF, getAccumz=.TRUE.)
mcnl = MCFT%getnml()
call MCFT%copyaccumz(accum_z)


nsx = (mcnl%numsx - 1)/2
nsy = nsx
etotal = float(mcnl%totnum_el)

io_int(1) = mcnl%totnum_el
call Message%WriteValue(' --> total number of BSE electrons in MC data set ', io_int, 1)
!=============================================
!=============================================
numEbins = MCFT%getnumEbins() 

allocate(EkeVs(numEbins),thick(numEbins))

do i=1,numEbins
  EkeVs(i) = mcnl%Ehistmin + float(i-1)*mcnl%Ebinsize
end do

!=============================================
! should we create a new file or open an existing file?
!=============================================
  lastEnergy = -1
  energyfile = trim(EMsoft%generateFilePath('EMdatapathname',emnl%energyfile))
  outname = trim(energyfile)

if (emnl%restart.eqv..TRUE.) then
! in this case we need to check whether or not the file exists, then open
! it and read the value of the last energy level that was simulated and written
! to that file; if this level is different from the lowest energy level we 
! know that there is at least one more level to be simulated.  If it is equal,
! then we can abort the program here.

  inquire(file=trim(outname), exist=f_exists)
  if (.not.f_exists) then 
    call Message%printError('EBSDmaster','restart HDF5 file does not exist')
  end if
  
!=============================================
! open the existing HDF5 file 
!=============================================
  datagroupname = 'EBSDmaster'
  HDF = HDF_T()

! Create a new file using the default properties.
  readonly = .TRUE.
  hdferr =  HDF%openFile(outname, readonly)

! all we need to get from the file is the lastEnergy parameter
groupname = SC_EMData
  hdferr = HDF%openGroup(groupname)
  hdferr = HDF%openGroup(datagroupname)

dataset = SC_lastEnergy
  call HDF%readDatasetInteger(dataset, hdferr, lastEnergy)

  call HDF%pop(.TRUE.)
end if
!=============================================
!=============================================
! crystallography section; 
verbose = .TRUE.

call cell%setFileName(mcnl%xtalname)
call Diff%setrlpmethod('WK')

if (emnl%restart.eqv..TRUE.) then 
  call Diff%setV(dble(EkeVs(lastEnergy-1)))
  call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, emnl%dmin, verbose, useHDF=HDF)
else
  call Diff%setV(dble(mcnl%EkeV))
  call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, emnl%dmin, verbose, useHDF=HDF)
end if

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

! then calculate density, average atomic number and average atomic weight
call cell%calcDensity(Z2percent)
io_real(1:3) = cell%getDensity()
density = io_real(1)
at_wt = io_real(2)
Ze = io_real(3)
call Message%WriteValue('Density, avA, avZ = ',io_real,3,"(2f10.5,',',f10.5)")

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
! this is where we determine the value for the thickness integration limit for the CalcLgh3 routine...


! then, for each energy determine the 95% histogram thickness
izzmax = 0
do iE = 1,numEbins
 do ix=-nsx/10,nsx/10
  do iy=-nsy/10,nsy/10
   istat = sum(accum_z(iE,:,ix,iy))
   izz = 1
   do while (sum(accum_z(iE,1:izz,ix,iy)).lt.(0.99*istat)) 
    izz = izz+1
   end do
   if (izz.gt.izzmax) izzmax = izz
  end do
 end do
 thick(iE) = dble(izzmax) * mcnl%depthstep
end do

izz = nint(maxval(thick)/mcnl%depthstep)
allocate(lambdaE(1:numEbins,1:izz),stat=istat)
do iE=1,numEbins
 call Diff%setV(dble(Ekevs(iE)))
 call Diff%CalcUcg(cell,(/0,0,0/))
 rlp = Diff%getrlp()
 nabsl = rlp%xgp
 do iz=1,izz
  lambdaE(iE,iz) = float(sum(accum_z(iE,iz,-nsx/10:nsx/10,-nsy/10:nsy/10)))/etotal
  lambdaE(iE,iz) = lambdaE(iE,iz) * exp(2.0*sngl(cPi)*(iz-1)*mcnl%depthstep/nabsl)
 end do
end do

! ---------- end of 'read Monte Carlo output file and extract necessary parameters' section
!=============================================
!=============================================

!=============================================
!=============================================
! ---------- a couple of initializations
   npy = emnl%npx
   allocate(svals(numset),stat=istat)
   gzero = 1  ! index of incident beam
   debug = 0  ! no longer used
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
  allocate(mLPNH(-emnl%npx:emnl%npx,-npy:npy,1,1:numsites),stat=istat)
  allocate(mLPSH(-emnl%npx:emnl%npx,-npy:npy,1,1:numsites),stat=istat)
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
! ---------- end allocate memory for the master patterns
!=============================================
!=============================================

! force dynamical matrix routine to read new Bethe parameters from file
! this will all be changed with the new version of the Bethe potentials
  call Diff%SetBetheParameters(EMsoft, .TRUE., emnl%BetheParametersFile)

if (emnl%restart.eqv..FALSE.) then
!=============================================
! create or update the HDF5 output file
!=============================================
  HDF = HDF_T() 

! Open an existing file or create a new file using the default properties.
  if (trim(energyfile).eq.trim(outname)) then
    hdferr =  HDF%openFile(outname)
  else
    hdferr =  HDF%createFile(outname)
  end if

! write the EMheader to the file
  datagroupname = 'EBSDmaster'
  call HDF%writeEMheader(dstr, tstrb, tstre, progname, datagroupname)

! open or create a namelist group to write all the namelist files into
groupname = SC_NMLfiles
  hdferr = HDF%createGroup(groupname)

! read the text file and write the array to the file
dataset = SC_EBSDmasterNML
  hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%nmldeffile)

! leave this group
  call HDF%pop()

! create a namelist group to write all the namelist files into
groupname = SC_NMLparameters
  hdferr = HDF%createGroup(groupname)

  call MPFT%writeHDFNameList(HDF, HDFnames, emnl)
  call Diff%writeBetheparameterNameList(HDF)

! leave this group
  call HDF%pop()

! then the remainder of the data in a EMData group
groupname = SC_EMData
  hdferr = HDF%createGroup(groupname)

! create the EBSDmaster group and add a HDF_FileVersion attribbute to it 
  hdferr = HDF%createGroup(datagroupname)
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

dataset = SC_lastEnergy
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then 
    hdferr = HDF%writeDatasetInteger(dataset, lastEnergy, overwrite)
  else
    hdferr = HDF%writeDatasetInteger(dataset, lastEnergy)
  end if
  
dataset = 'Z2percent'  ! SC_Z2percent
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then 
    hdferr = HDF%writeDatasetFloatArray(dataset, Z2percent, cell%getNatomtype(), overwrite)
  else
    hdferr = HDF%writeDatasetFloatArray(dataset, Z2percent, cell%getNatomtype())
  end if

  if (emnl%Esel.eq.-1) then
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
  else
dataset = SC_numEbins
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then 
      hdferr = HDF%writeDatasetInteger(dataset, one, overwrite)
    else
      hdferr = HDF%writeDatasetInteger(dataset, one)
    end if
  
dataset = SC_selE
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then 
      hdferr = HDF%writeDatasetFloat(dataset, selE, overwrite)
    else
      hdferr = HDF%writeDatasetFloat(dataset, selE)
    end if
  end if  

! create the hyperslabs and write zeroes to them for now
dataset = SC_mLPNH
  dims4 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins, numsites /)
  cnt4 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1, numsites /)
  offset4 = (/ 0, 0, 0, 0 /)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then 
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims4, offset4, cnt4, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims4, offset4, cnt4)
  end if

dataset = SC_mLPSH
  dims4 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins, numsites /)
  cnt4 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1, numsites /)
  offset4 = (/ 0, 0, 0, 0 /)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then 
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims4, offset4, cnt4, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims4, offset4, cnt4)
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

  call HDF%pop(.TRUE.)
end if

!=============================================
!=============================================
! do we need to precompute the Legendre array for the new lattitudinal grid values?
if (doLegendre.eqv..TRUE.) then
  call Message%printMessage(' Computing Legendre lattitudinal grid values')
  allocate(diagonal(2*emnl%npx+1),upd(2*emnl%npx+1))
  diagonal = 0.D0
  upd = (/ (dble(i) / dsqrt(4.D0 * dble(i)**2 - 1.D0), i=1,2*emnl%npx+1) /)
  call dsterf(2*emnl%npx-1, diagonal, upd, info) 
! the eigenvalues are stored from smallest to largest and we need them in the opposite direction
  allocate(LegendreArray(2*emnl%npx+1))
  LegendreArray(1:2*emnl%npx+1) = diagonal(2*emnl%npx+1:1:-1)
! set the center eigenvalue to 0
  LegendreArray(emnl%npx+1) = 0.D0
  deallocate(diagonal, upd)
end if

!=============================================
!=============================================
! figure out what the start energy value is for the energyloop
if (lastEnergy.ne.-1) then
  Estart = lastEnergy-1
  if (Estart.eq.0) then 
    call Message%printMessage('All energy levels are present in the HDF5 file')
    call Message%printMessage('No further computations needed.')
    stop 'Program halted.'
  end if
else
  Estart = numEbins
end if

!=============================================
!=============================================
! ---------- from here on, we need to repeat the entire computation for each energy value
! so this is where we could in principle implement an OpenMP approach; alternatively, 
! we could do the inner loop over the incident beam directions in OpenMP (probably simpler)

! we use two times, one (1) for each individual energy level, the other (2) for the overall time
call timer%Time_tick(2)
reflist = gvectors_T()

energyloop: do iE=Estart,1,-1
 if (emnl%uniform.eqv..FALSE.) then
! is this a single-energy run ?
   if (emnl%Esel.ne.-1) then
     if (emnl%Esel.ne.iE) CYCLE energyloop
   end if
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
    verbose = .TRUE.
!   call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, emnl%dmin, useHDF=HDF, noLUT=.TRUE., verbose=verbose)
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

   if (doLegendre.eqv..FALSE.) then
    call kvec%set_mapmode('RoscaLambert')
    if (usehex) then
      call kvec%Calckvectors(cell, SG, Diff, (/ 0.D0, 0.D0, 0.D0 /),emnl%npx,npy, ijmax,usehex)
    else 
      call kvec%Calckvectors(cell, SG, Diff, (/ 0.D0, 0.D0, 0.D0 /),emnl%npx,npy, ijmax,usehex)
    end if
   else 
    call kvec%set_mapmode('RoscaLambertLegendre')
    if (usehex) then
      call kvec%Calckvectors(cell, SG, Diff, (/ 0.D0, 0.D0, 0.D0 /),emnl%npx,npy, ijmax,usehex, LegendreArray)
    else 
      call kvec%Calckvectors(cell, SG, Diff, (/ 0.D0, 0.D0, 0.D0 /),emnl%npx,npy, ijmax,usehex, LegendreArray)
    end if
   end if
   numk = kvec%get_numk()
   io_int(1)=numk
   call Message%WriteValue('# independent beam directions to be considered = ', io_int, 1, "(I8)")

! convert part of the kvector linked list into arrays for OpenMP
  allocate(karray(4,numk), kij(3,numk),stat=istat)
! point to the first beam direction
  ktmp => kvec%get_ListHead()
! and loop through the list, keeping k, kn, and i,j
  karray(1:3,1) = sngl(ktmp%k(1:3))
  karray(4,1) = sngl(ktmp%kn)
  kij(1:3,1) = (/ ktmp%i, ktmp%j, ktmp%hs /)
   do ik=2,numk
     ktmp => ktmp%next
     karray(1:3,ik) = sngl(ktmp%k(1:3))
     karray(4,ik) = sngl(ktmp%kn)
     kij(1:3,ik) = (/ ktmp%i, ktmp%j, ktmp%hs /)
   end do
! and remove the linked list
  call kvec%Delete_kvectorlist()

  verbose = .FALSE.
  totstrong = 0
  totweak = 0

! ---------- end of "create the incident beam directions list"
!=============================================

! here's where we introduce the OpenMP calls, to speed up the overall calculations...

! set the number of OpenMP threads 
  if (emnl%nthreads.eq.0) then 
    nthreads = OMP_GET_MAX_THREADS()
  else
    nthreads = emnl%nthreads
  end if
  call OMP_SET_NUM_THREADS(nthreads)
  io_int(1) = nthreads
  call Message%WriteValue(' Attempting to set number of threads to ',io_int, 1, frm = "(I4)")

! use OpenMP to run on multiple cores ... 
!$OMP PARALLEL COPYIN(rlp) &
!$OMP& PRIVATE(DynMat,Sgh,Lgh,ik,FN,TID,kn,ipx,ipy,ix,iequiv,nequiv,reflist,firstw) &
!$OMP& PRIVATE(kkk,nns,nnw,nref,svals,io_int)

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

     call reflist%Initialize_ReflectionList(cell, SG, Diff, FN, kkk, emnl%dmin, verbose)
     nref = reflist%get_nref()
! ---------- end of "create the master reflection list"
!=============================================

! determine strong and weak reflections
     nullify(firstw)
     nns = 0
     nnw = 0
     call reflist%Apply_BethePotentials(Diff, firstw, nns, nnw)

     if (emnl%kinematical.eqv..FALSE.) then 
! generate the dynamical matrix
       if (allocated(DynMat)) deallocate(DynMat)
       allocate(DynMat(nns,nns))
       call reflist%GetDynMat(cell, Diff, firstw, DynMat, nns, nnw)
       totstrong = totstrong + nns
       totweak = totweak + nnw
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
     if (allocated(Sgh)) deallocate(Sgh)
     if (allocated(Lgh)) deallocate(Lgh)
     allocate(Sgh(nns,nns,numset),Lgh(nns,nns))
     Sgh = czero
     Lgh = czero
     call reflist%getSghfromLUT(Diff,nns,numset,Sgh)


! solve the dynamical eigenvalue equation for this beam direction  
     kn = karray(4,ik)
     call reflist%CalcLgh(DynMat,Lgh,dble(thick(iE)),dble(kn),nns,gzero,mcnl%depthstep,lambdaE(iE,1:izzmax),izzmax)

! sum over the element-wise (Hadamard) product of the Lgh and Sgh arrays 
     svals = 0.0
     do ix=1,numset
       svals(ix) = real(sum(Lgh(1:nns,1:nns)*Sgh(1:nns,1:nns,ix)))
     end do
     svals = svals/float(sum(nat(1:numset)))


! and store the resulting svals values, applying point group symmetry where needed.
     ipx = kij(1,ik)
     ipy = kij(2,ik)
     ipz = kij(3,ik)
!
     if (usehex) then 
       call L%Apply3DPGSymmetry(cell,SG,ipx,ipy,ipz,emnl%npx,iequiv,nequiv,usehex)
     else
       if ((SG%getSpaceGroupNumber().ge.195).and.(SG%getSpaceGroupNumber().le.230)) then
         call L%Apply3DPGSymmetry(cell,SG,ipx,ipy,ipz,emnl%npx,iequiv,nequiv,cubictype=SamplingType)
       else
         call L%Apply3DPGSymmetry(cell,SG,ipx,ipy,ipz,emnl%npx,iequiv,nequiv)
       end if
     end if
!$OMP CRITICAL
  if (emnl%combinesites.eqv..FALSE.) then
     do ix=1,nequiv
       if (iequiv(3,ix).eq.-1) mLPSH(iequiv(1,ix),iequiv(2,ix),1,1:numset) = svals(1:numset)
       if (iequiv(3,ix).eq.1) mLPNH(iequiv(1,ix),iequiv(2,ix),1,1:numset) = svals(1:numset)
     end do
  else
     do ix=1,nequiv
       if (iequiv(3,ix).eq.-1) mLPSH(iequiv(1,ix),iequiv(2,ix),1,1) = sum(svals)
       if (iequiv(3,ix).eq.1) mLPNH(iequiv(1,ix),iequiv(2,ix),1,1) = sum(svals)
     end do
  end if
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
  deallocate(karray, kij)

  if (usehex) then
! and finally, we convert the hexagonally sampled array to a square Lambert projection which will be used 
! for all EBSD pattern interpolations;  we need to do this for both the Northern and Southern hemispheres

! we begin by allocating auxiliary arrays to hold copies of the hexagonal data; the original arrays will
! then be overwritten with the newly interpolated data.
    allocate(auxNH(-emnl%npx:emnl%npx,-npy:npy,1,1:numsites),stat=istat)
    allocate(auxSH(-emnl%npx:emnl%npx,-npy:npy,1,1:numsites),stat=istat)
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
          mLPNH(i,j,1,1:numsites) = auxNH(nix,niy,1,1:numsites)*dxm*dym + auxNH(nixp,niy,1,1:numsites)*dx*dym + &
                               auxNH(nix,niyp,1,1:numsites)*dxm*dy + auxNH(nixp,niyp,1,1:numsites)*dx*dy
          mLPSH(i,j,1,1:numsites) = auxSH(nix,niy,1,1:numsites)*dxm*dym + auxSH(nixp,niy,1,1:numsites)*dx*dym + &
                               auxSH(nix,niyp,1,1:numsites)*dxm*dy + auxSH(nixp,niyp,1,1:numsites)*dx*dy
        end if
      end do
    end do
    deallocate(auxNH, auxSH)
  end if

! make sure that the outer pixel rim of the mLPSH patterns is identical to
! that of the mLPNH array.
  mLPSH(-emnl%npx,-emnl%npx:emnl%npx,1,1:numsites) = mLPNH(-emnl%npx,-emnl%npx:emnl%npx,1,1:numsites)
  mLPSH( emnl%npx,-emnl%npx:emnl%npx,1,1:numsites) = mLPNH( emnl%npx,-emnl%npx:emnl%npx,1,1:numsites)
  mLPSH(-emnl%npx:emnl%npx,-emnl%npx,1,1:numsites) = mLPNH(-emnl%npx:emnl%npx,-emnl%npx,1,1:numsites)
  mLPSH(-emnl%npx:emnl%npx, emnl%npx,1,1:numsites) = mLPNH(-emnl%npx:emnl%npx, emnl%npx,1,1:numsites)


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
        masterSPNH(i,j,1) = InterpolateLambert(xyz, mLPNH, emnl%npx, numsites)
        masterSPSH(i,j,1) = InterpolateLambert(xyz, mLPSH, emnl%npx, numsites)
      end if
    end do
  end do

! since these computations can take a long time, here we store 
! all the output at the end of each pass through the energyloop.

  io_int(1) = nint(float(totstrong)/float(numk))
  call Message%WriteValue(' -> Average number of strong reflections = ',io_int, 1, "(I5)")
  io_int(1) = nint(float(totweak)/float(numk))
  call Message%WriteValue(' -> Average number of weak reflections   = ',io_int, 1, "(I5)")

  call timer%makeTimeStamp()
  dstr = timer%getDateString() 
  tstre = timer%getTimeString()
 end if  ! (emnl%uniform.eqv..FALSE.) 

  datagroupname = 'EBSDmaster'
  HDF = HDF_T() 

! open the existing file using the default properties.
  hdferr =  HDF%openFile(outname)

! update the time string
groupname = SC_EMheader
  hdferr = HDF%openGroup(groupname)
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
  if (iE.eq.numEbins) then 
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then     
      hdferr = HDF%writeDatasetFloat(dataset, tstop, overwrite)
    else
      hdferr = HDF%writeDatasetFloat(dataset, tstop)
    end if
  else
    hdferr = HDF%writeDatasetFloat(dataset, tstop, overwrite)
  end if

  call HDF%pop()
  call HDF%pop()
  
groupname = SC_EMData
  hdferr = HDF%openGroup(groupname)
  hdferr = HDF%openGroup(datagroupname)

! update the current energy level counter, so that the restart option will function
dataset = SC_lastEnergy
  hdferr = HDF%writeDataSetInteger(dataset, iE, overwrite)

! add data to the hyperslab
dataset = SC_mLPNH
  dims4 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins, numsites /)
  cnt4 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1, numsites /)
  offset4 = (/ 0, 0, iE-1, 0 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims4, offset4, cnt4, insert)

dataset = SC_mLPSH
  dims4 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins, numsites /)
  cnt4 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1, numsites /)
  offset4 = (/ 0, 0, iE-1, 0 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims4, offset4, cnt4, insert)

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

  call HDF%pop(.TRUE.)

 if ((emnl%Esel.eq.-1).and.(iE.ne.1)) then 
  call Message%printMessage(' Intermediate data stored in file '//trim(emnl%energyfile), frm = "(A/)")
 end if

 if ((emnl%Esel.eq.-1).and.(iE.eq.1)) then 
  call Message%printMessage(' Final data stored in file '//trim(emnl%energyfile), frm = "(A/)")
 end if

end do energyloop

call timer%Time_tock(2) 
io_int(1) = timer%getInterval(2)
call Message%WriteValue(' Total execution time [s] ',io_int,1)

! if requested, we notify the user that this program has completed its run
if (trim(EMsoft%getConfigParameter('EMNotify')).ne.'Off') then
  if (trim(emnl%Notify).eq.'On') then 
    NumLines = 3
    allocate(MessageLines(NumLines))

    call hostnm(c)
 
    MessageLines(1) = 'EMEBSDmaster program has ended successfully'
    MessageLines(2) = 'Master pattern data stored in '//trim(outname)
    write (exectime,"(F10.4)") timer%getInterval(2)  
    MessageLines(3) = 'Total execution time [s]: '//trim(exectime)
    SlackUsername = 'EMsoft on '//trim(c)
    i = PostMessage(EMsoft, MessageLines, NumLines, SlackUsername)
  end if
end if

end associate 

end subroutine EBSDmaster_


end module mod_EBSDmaster