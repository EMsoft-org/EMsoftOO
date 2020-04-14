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

module mod_ECP
  !! author: MDG 
  !! version: 1.0 
  !! date: 03/14/20
  !!
  !! class definition for the EMECP program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMECP program
type, public :: ECPNameListType
  integer(kind=irg)       :: nthreads
  integer(kind=irg)       :: npix
  real(kind=sgl)          :: conesemiangle
  character(1)            :: maskpattern
  character(fnlen)        :: xtalname
  character(fnlen)        :: energyfile
  character(fnlen)        :: masterfile
  character(fnlen)        :: datafile
  character(fnlen)        :: anglefile
  character(3)            :: eulerconvention
  real(kind=sgl)          :: gammavalue
  character(3)            :: outputformat
  real(kind=dbl)          :: sampletilt
  real(kind=sgl)          :: workingdistance
  real(kind=sgl)          :: Rin
  real(kind=sgl)          :: Rout
end type ECPNameListType

type ECPDetectorType
  real(kind=sgl),allocatable  :: rgx(:,:), rgy(:,:), rgz(:,:)  ! auxiliary detector arrays needed for interpolation
  integer(kind=irg)           :: npolar, nazimuth
  ! real(kind=sgl),allocatable  :: accum_e_detector(:,:,:)
  ! type(EBSDPixel),allocatable :: detector(:,:) 
end type ECPDetectorType

type, public :: IncidentListECP
  integer(kind=irg)               :: i, j, numk
  real(kind=dbl)                  :: k(3)
  type(IncidentListECP),pointer   :: next
end type IncidentListECP

! class definition
type, public :: ECP_T
private 
  character(fnlen)                :: nmldeffile = 'EMECP.nml'
  type(ECPNameListType),public    :: nml 
  type(ECPDetectorType)           :: det
  type(IncidentListECP),pointer   :: klist

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: get_ListHead_
  procedure, pass(self) :: ECP_
  procedure, pass(self) :: ECPGenerateDetector_
  procedure, pass(self) :: ECPGetWeightFactors_
  procedure, pass(self) :: GetVectorsCone_
  procedure, pass(self) :: get_numk_
  procedure, pass(self) :: set_numk_
  procedure, pass(self) :: CalcECPatternSingle_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: get_ListHead => get_ListHead_
  generic, public :: ECP => ECP_
  generic, public :: ECPGenerateDetector => ECPGenerateDetector_
  generic, public :: ECPGetWeightFactors => ECPGetWeightFactors_
  generic, public :: GetVectorsCone => GetVectorsCone_
  generic, public :: get_numk => get_numk_
  generic, public :: set_numk => set_numk_
  generic, public :: CalcECPatternSingle => CalcECPatternSingle_

end type ECP_T

! the constructor routine for this class 
interface ECP_T
  module procedure ECP_constructor
end interface ECP_T

contains

!--------------------------------------------------------------------------
type(ECP_T) function ECP_constructor( nmlfile ) result(ECP)
!! author: MDG 
!! version: 1.0 
!! date: 03/14/20
!!
!! constructor for the ECP_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 
integer(kind=irg)            :: istat 

if (present(nmlfile)) call ECP%readNameList(nmlfile)

allocate(ECP%klist,stat=istat)
nullify(ECP%klist%next)

end function ECP_constructor

!--------------------------------------------------------------------------
subroutine ECP_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 03/14/20
!!
!! destructor for the ECP_T Class
 
IMPLICIT NONE

type(ECP_T), INTENT(INOUT)  :: self 

call reportDestructor('ECP_T')

end subroutine ECP_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!! author: MDG 
!! version: 1.0 
!! date: 03/14/20
!!
!! read the namelist from an nml file for the ECP_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(ECP_T), INTENT(INOUT)      :: self
character(fnlen),INTENT(IN)      :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)      :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                   :: EMsoft 
type(IO_T)                       :: Message       
logical                          :: skipread = .FALSE.

! parameters for the standard ECP program
integer(kind=irg)       :: nthreads
integer(kind=irg)       :: npix
real(kind=sgl)          :: conesemiangle
character(fnlen)        :: xtalname
character(fnlen)        :: energyfile
character(fnlen)        :: masterfile
character(fnlen)        :: datafile
character(1)            :: maskpattern
character(fnlen)        :: anglefile
character(3)            :: eulerconvention
real(kind=sgl)          :: gammavalue
character(3)            :: outputformat
real(kind=dbl)          :: sampletilt
real(kind=sgl)          :: workingdistance
real(kind=sgl)          :: Rin
real(kind=sgl)          :: Rout


namelist /ECPlist/ xtalname, anglefile, nthreads, conesemiangle, npix, maskpattern, eulerconvention, Rin, Rout, &
                   energyfile, masterfile, datafile, gammavalue, outputformat, sampletilt, workingdistance

! set the common default parameters 
npix = 200                              ! number of pixels in final image (npix x npix)
nthreads = 1                            ! number of OpenMP threads
conesemiangle = 0.0                     ! beam convergence in mrad (either ktmax or thetac must be given)
xtalname = 'undefined'                  ! initial value to check that the keyword is present in the nml file
energyfile = 'undefined'
masterfile = 'undefined'
datafile = 'undefined'
maskpattern = 'y'
anglefile = 'undefined'
eulerconvention = 'hkl'
gammavalue = 1.0
outputformat = 'gui'
sampletilt = 0.D0
workingdistance = 13.0
Rin = 2.0
Rout = 6.0

! read the name list, depending on the class type 
if (.not.skipread) then
! read the namelist file
  open(UNIT=dataunit,FILE=trim(EMsoft%toNativePath(nmlfile)),DELIM='apostrophe',STATUS='old')
  read(UNIT=dataunit,NML=ECPlist)
  close(UNIT=dataunit,STATUS='keep')

! check for required entries
  if (trim(xtalname).eq.'undefined') then
    call Message%printError('EMECP:',' crystal file name is undefined in '//nmlfile)
  end if
end if 

self%nml%npix = npix
self%nml%nthreads = nthreads
self%nml%conesemiangle = conesemiangle
self%nml%datafile = datafile
self%nml%xtalname = xtalname
self%nml%energyfile = energyfile
self%nml%masterfile = masterfile
self%nml%maskpattern = maskpattern
self%nml%anglefile = anglefile
self%nml%eulerconvention = eulerconvention
self%nml%gammavalue = gammavalue
self%nml%outputformat = outputformat
self%nml%sampletilt = sampletilt
self%nml%workingdistance = workingdistance
self%nml%Rin = Rin
self%nml%Rout = Rout

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!! author: MDG 
!! version: 1.0 
!! date: 03/14/20
!!
!! pass the namelist for the ECP_T Class to the calling program

IMPLICIT NONE 

class(ECP_T), INTENT(INOUT)   :: self
type(ECPNameListType)         :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
function get_numk_(self) result(numk)
!! author: MDG 
!! version: 1.0 
!! date: 03/17/20
!!
!! get the number of incident wave vectors in the list 

IMPLICIT NONE 

class(ECP_T), INTENT(INOUT)   :: self
integer(kind=irg)             :: numk

numk = self%klist%numk

end function get_numk_

!--------------------------------------------------------------------------
subroutine set_numk_(self, numk)
!! author: MDG 
!! version: 1.0 
!! date: 03/17/20
!!
!! set the number of incident wave vectors in the list 

IMPLICIT NONE 

class(ECP_T), INTENT(INOUT)   :: self
integer(kind=irg),INTENT(OUT) :: numk

self%klist%numk = numk

end subroutine set_numk_

!--------------------------------------------------------------------------
recursive function get_ListHead_(self) result(klist)
!! author: MDG 
!! version: 1.0 
!! date: 02/12/20
!!
!! return a pointer ot the head of the list 

IMPLICIT NONE

class(ECP_T), INTENT(INOUT)     :: self
type(IncidentListECP), pointer  :: klist

klist => self%klist 

end function get_ListHead_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!! author: MDG 
!! version: 1.0 
!! date: 03/14/20
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 
use ISO_C_BINDING

IMPLICIT NONE

class(ECP_T), INTENT(INOUT)             :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames 

integer(kind=irg),parameter             :: n_int = 2, n_real = 6
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( ecpnl => self%nml )

! create the group for this namelist
groupname = trim(HDFnames%get_NMLlist())
hdferr = HDF%createGroup(groupname)

! write all the single integers
io_int = (/ ecpnl%nthreads, ecpnl%npix /)
intlist(1) = 'nthreads'
intlist(2) = 'npix'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ ecpnl%conesemiangle, sngl(ecpnl%sampletilt), &
             ecpnl%gammavalue, ecpnl%workingdistance, ecpnl%Rin, ecpnl%Rout /)
reallist(1) = 'conesemiangle'
reallist(2) = 'sampletilt'
reallist(3) = 'gammavalue'
reallist(4) = 'workingdistance'
reallist(5) = 'Rin'
reallist(6) = 'Rout'
call HDF%writeNMLreals(io_real, reallist, n_real)

! write all the strings

dataset = SC_energyfile
line2(1) = ecpnl%energyfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECPNameList: unable to create energyfile dataset', hdferr)

dataset = SC_datafile
line2(1) = ecpnl%datafile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECPNameList: unable to create datafile dataset', hdferr)

dataset = SC_xtalname
line2(1) = ecpnl%xtalname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECPNameList: unable to create xtalname dataset', hdferr)

dataset = SC_masterfile
line2(1) = ecpnl%masterfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECPNameList: unable to create masterfile dataset', hdferr)

dataset = SC_maskpattern
line2(1) = ecpnl%maskpattern
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECPNameList: unable to create maskpattern dataset', hdferr)

dataset = SC_anglefile
line2(1) = ecpnl%anglefile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECPNameList: unable to create anglefile dataset', hdferr)

dataset = SC_outputformat
line2(1) = trim(ecpnl%outputformat)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECPNameList: unable to create outputformat dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine ECP_(self, EMsoft, progname, nmldeffile)
!! author: MDG 
!! version: 1.0 
!! date: 03/14/20
!!
!! perform the computations

use mod_EMsoft
use mod_crystallography
use mod_symmetry
use mod_Lambert
! use mod_initializers
use mod_io
use mod_so3
use mod_timing
use mod_rotations
use mod_quaternions
use mod_MPfiles
use mod_MCfiles
use HDF5
use mod_HDFsupport
use mod_HDFnames
use ISO_C_BINDING
use omp_lib
! use distortion
use mod_filters 
use stringconstants

IMPLICIT NONE 

class(ECP_T), INTENT(INOUT)             :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 
character(fnlen), INTENT(INOUT)         :: nmldeffile

type(MCfile_T)                          :: MCFT
type(MPfile_T)                          :: MPFT
type(HDF_T)                             :: HDF
type(HDFnames_T)                        :: HDFnames
type(so3_T)                             :: SO
type(Timing_T)                          :: timer
type(IO_T)                              :: Message
type(Quaternion_T)                      :: quat, qq 
type(QuaternionArray_T)                 :: qAR

type(ECPmasterNameListType)             :: mpnl
type(MCOpenCLNameListType)              :: mcnl

type(IncidentListECP),pointer           :: ktmp
real(kind=dbl),allocatable              :: klist(:,:)
real(kind=sgl),allocatable              :: mLPNH(:,:) , mLPSH(:,:)


integer(kind=irg)                       :: npx,npy,numset,istat,val, numangles
integer(kind=irg),allocatable           :: ATOM_type(:)
real(kind=dbl)                          :: EkeV
real(kind=sgl)                          :: dmin, FN(3), tstop
real(kind=sgl),allocatable              :: mask(:,:), lx(:), ly(:)
integer(kind=irg)                       :: maskradius, io_int(1), hdferr
logical                                 :: verbose
real(kind=sgl),allocatable              :: anglewf(:)
integer(kind=irg)                       :: nsig, isig, isigp
integer(kind=irg)                       :: numk, nix, niy, nixp, niyp, i, j, ierr, &
                                           ipx, ipy, iang, idir
real(kind=dbl)                          :: scl, x, dx, dy, dxm, dym, wf
real(kind=dbl)                          :: dc(3), ixy(2), qu(4)
integer(kind=irg),allocatable           :: kij(:,:)
real(kind=sgl),allocatable              :: ECPpattern(:,:), ECPpatternintd(:,:), ECP_tmp(:,:,:)
integer(kind=irg),allocatable           :: ECPpatterninteger(:,:), ECPpatternad(:,:)
real(kind=sgl)                          :: time_start, time_end, ma, mi
character(len=1),allocatable            :: bpat(:,:), bpat_tmp(:,:,:)
integer(kind=irg)                       :: TID, nthreads
real(kind=dbl)                          :: dp, MCangle
logical                                 :: switchwfoff = .FALSE.

integer(HSIZE_T), dimension(1:3)        :: hdims, offset 
integer(HSIZE_T)                        :: dims3(3)
character(fnlen,kind=c_char)            :: line2(1)
character(fnlen)                        :: groupname, dataset, attributename, HDF_FileVersion, fname
character(11)                           :: dstr
character(15)                           :: tstrb
character(15)                           :: tstre
character(fnlen)                        :: datafile
logical                                 :: overwrite = .TRUE., insert = .TRUE., g_exists

! open the HDF interface
call openFortranHDFInterface()

! construct the HDF and HDFnames classes
HDF = HDF_T() 
HDFnames = HDFnames_T() 

! associate a couple of variable types
associate( enl => self%nml, ECPdetector => self%det, ECPMCdata => MCFT%MCDT )
call MPFT%setModality('ECP')

!=================================================================
! 1. read the angle array from file
!=================================================================
verbose = .TRUE.
call setRotationPrecision('d')

fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%anglefile))
call SO%nullifyList()
call SO%getOrientationsfromFile(fname)
numangles = SO%getListCount('FZ')
call SO%listtoQuaternionArray( qAR )
call SO%delete_FZlist()

! set the HDF group names for reading the MC input file 
call HDFnames%set_ProgramData(SC_MCOpenCL) 
call HDFnames%set_NMLlist(SC_MCCLNameList) 
call HDFnames%set_NMLfilename(SC_MCOpenCLNML)
! call HDFnames%get_AllNames()
!=================================================================
! 2. read the Monte Carlo data file (HDF format)
!=================================================================
fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%energyfile))
call MCFT%setFileName(fname)
call MCFT%readMCfile(HDF, HDFnames, getAccume=.TRUE.)
mcnl = MCFT%getnml()

! set the HDF group names for reading the MP input file 
call HDFnames%set_ProgramData(SC_ECPmaster) 
call HDFnames%set_NMLlist(SC_ECPMasterNameList) 
call HDFnames%set_NMLfilename(SC_ECPmasterNML) 
! call HDFnames%get_AllNames()
!=================================================================
! 3. read ECP master pattern file (HDF format); and sum to 2D arrays
!=================================================================
fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfile))
call MPFT%setFileName(fname)
call MPFT%readMPfile(HDF, HDFnames, mpnl, getmLPNH=.TRUE., getmLPSH=.TRUE.)
call MPFT%copysummLPNH(mLPNH)
call MPFT%copysummLPSH(mLPSH)

! reset the HDFnames to the ones needed by the EMECP program
call HDFnames%set_ProgramData(SC_ECP) 
call HDFnames%set_NMLlist(SC_ECPNameList) 
call HDFnames%set_NMLfilename(SC_EMECPNML) 
! call HDFnames%get_AllNames()

!=================================================================
! 4. generate the detector 
!=================================================================
call self%ECPGenerateDetector(verbose)

!=================================================================
! get the weight factors 
!=================================================================
nsig = nint((enl%conesemiangle) + abs(enl%sampletilt)) + 1
allocate(anglewf(1:nsig),stat=istat)

call Message%printMessage(' -> Calculating weight factors', frm = "(A)" )
call self%ECPGetWeightFactors(mcnl, MCFT, anglewf, nsig, verbose=.TRUE.)

!=================================================================
! check if there are enough angles in MC for detector geometry
!=================================================================
if (mcnl%sigend .lt. (abs(enl%sampletilt) + enl%conesemiangle)) then
  call Message%printMessage('Not enough angles in Monte carlo file...interpolation will be done without &
  appropriate weight factors',frm = "(A)")
  switchwfoff = .TRUE.
end if

if ((-mcnl%sigend .gt. (enl%conesemiangle - abs(enl%sampletilt))) .and. (switchwfoff .eqv. .FALSE.)) then
  call Message%printMessage('Not enough angles in Monte carlo file...interpolation will be done without &
  appropriate weight factors',frm = "(A)")
  switchwfoff = .TRUE.
end if

!=================================================================
! generate list of incident vectors
!=================================================================
numk = 0
call self%GetVectorsCone()
numk = self%get_numk()
allocate(kij(2,numk),klist(3,numk),stat=istat)

io_int(1) = numk
call Message%WriteValue('Number of beams for which interpolation will be done = ',io_int,1) 

ktmp => self%klist
! converting to array for OpenMP parallelization
do i = 1,numk
   klist(1:3,i) = ktmp%k(1:3)
   kij(1:2,i) = (/ktmp%i,ktmp%j/)
   ktmp => ktmp%next
end do

!=================================================================
! open the output file (only thread 0 can write to this file)
!================================================================
! initialize the timing routines
timer = Timing_T()
tstrb = timer%getTimeString()
dstr = timer%getDateString()
tstre = ''
call timer%Time_tick(1)

! Create a new file using the default properties.
datafile = trim(EMsoft%generateFilePath('EMdatapathname',enl%datafile))
hdferr =  HDF%createFile(datafile)

! write the EMheader to the file
groupname = trim(HDFnames%get_ProgramData())
call HDF%writeEMheader(EMsoft,dstr, tstrb, tstre, progname, groupname)

! create a namelist group to write all the namelist files into
groupname = trim(HDFnames%get_NMLfiles())
hdferr = HDF%createGroup(groupname)

groupname = trim(HDFnames%get_ProgramData())
hdferr = HDF%createGroup(groupname)

! read the text file and write the array to the file
dataset = trim(HDFnames%get_NMLfilename())
hdferr = HDF%writeDatasetTextFile(dataset, nmldeffile)

call HDF%pop()
call HDF%pop()

! create a NMLparameters group to write all the namelist entries into
groupname = trim(HDFnames%get_NMLparameters())
hdferr = HDF%createGroup(groupname)

call self%writeHDFNameList(HDF, HDFnames)

! and leave this group
call HDF%pop()

! then the remainder of the data in a EMData group
groupname = trim(HDFnames%get_EMData())
hdferr = HDF%createGroup(groupname)

! create the ECP group and add a HDF_FileVersion attribbute to it 
groupname = trim(HDFnames%get_ProgramData())
  hdferr = HDF%createGroup(groupname)

! before Feb. 19, 2019, an undetected error caused all patterns to be upside down in the Kikuchi bands only,
! not in the background intensity profile.  This was compensated by a pattern flip of all experimental 
! patterns in the dictionary indexing program, but when taking individual patterns from this program, they
! are actually upside down in all versions through HDF_FileVersion 4.0.  As of 4.1, the patterns are in the
! correct orientation.  This was detected by manually indexing a simulated pattern.

HDF_FileVersion = '4.0'
HDF_FileVersion = cstringify(HDF_FileVersion)
attributename = SC_HDFFileVersion
hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)

! =====================================================
! The following write commands constitute HDF_FileVersion = 4.0 and above
! =====================================================
! we need to write the image dimensions
dataset = SC_npix
hdferr = HDF%writeDatasetInteger(dataset, enl%npix) 

dataset = SC_numangledictionary
hdferr = HDF%writeDatasetInteger(dataset, numangles) 
! =====================================================
! end of HDF_FileVersion = 4.0 and above write statements
! =====================================================

!====================================
! new since Release 4.3: add a Manufacturer string (null terminated)
dataset = SC_Manufacturer
line2(1) = 'EMsoft'
line2(1) = cstringify(line2(1))
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
!====================================

! and we leave this group open for further data output ... 

! define the circular mask
allocate(mask(1:enl%npix, 1:enl%npix),stat=istat)
if (istat .ne. 0) then
   call Message%printError('ECP_','could not allocate mask array')
end if

mask = 1.0
if (enl%maskpattern .eq. 'y') then
! create the circular mask
  maskradius = (float(enl%npix)/2.0)**2
  allocate(lx(1:enl%npix), ly(1:enl%npix), stat=istat)
  lx = (/ (float(i),i=1,enl%npix) /) - float(enl%npix+1)/2.0
  ly = (/ (float(i),i=1,enl%npix) /) - float(enl%npix+1)/2.0
  do i= 1,enl%npix
    do j= 1,enl%npix
      if ((lx(i)**2+ly(j)**2).gt.maskradius) mask(i,j) = 0.0
    end do
  end do
  deallocate(lx, ly)
end if

! determine the scale factor for the Lambert interpolation; the square has
! an edge length of 2 x sqrt(pi/2)
scl = dble(mpnl%npx)

! allocate various arrays
allocate(ECPpattern(1:enl%npix, 1:enl%npix),&
         ECPpatternintd(1:enl%npix,1:enl%npix),&
         ECPpatterninteger(1:enl%npix,1:enl%npix),&
         ECPpatternad(1:enl%npix,1:enl%npix),stat=istat)

ECPpattern = 0.0
ECPpatternintd = 0.0
ECPpatterninteger = 0
ECPpatternad = 0

! create the pattern data set in the output file 
dataset = SC_ECpatterns

if (enl%outputformat .eq. 'bin') then
    allocate(bpat(1:enl%npix,1:enl%npix),bpat_tmp(1:enl%npix,1:enl%npix,1),stat=istat)
    if (istat .ne. 0) call Message%printError('ECpattern','cannot allocate bpat array')
    bpat = char(nint(255.0*ECPpattern))

! write dictionary pattern to h5 file
    offset = (/ 0, 0, 0 /)
    hdims = (/ enl%npix, enl%npix, numangles /)
    dims3 = (/ enl%npix, enl%npix, 1 /)
    hdferr = HDF%writeHyperslabCharArray(dataset, bpat_tmp, hdims, offset, dims3)
end if

if (enl%outputformat .eq. 'gui') then
    allocate(ECP_tmp(enl%npix, enl%npix, 1),stat=istat)
    offset = (/ 0, 0, 0 /)
    hdims = (/ enl%npix, enl%npix, numangles /)
    dims3 = (/ enl%npix, enl%npix, 1 /)
    hdferr = HDF%writeHyperslabFloatArray(dataset, ECP_tmp, hdims, offset, dims3)
end if

! set the number of OpenMP threads
io_int(1) = enl%nthreads
call Message%WriteValue(' Attempting to set number of threads to ',io_int,1,"(I4)")

call OMP_SET_NUM_THREADS(enl%nthreads)

! use OpenMP to run on multiple cores
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(TID,nthreads,dc,ixy,istat,nix,niy,nixp,niyp,dx,dy,dxm,dym,MCangle,isig,dp,isigp) &
!$OMP& PRIVATE(ipx,ipy,ECPpattern,bpat,ECPpatternintd,ma,mi,offset,hdims,dims3,hdferr,qq,idir,wf)

TID = OMP_GET_THREAD_NUM()
nthreads = OMP_GET_NUM_THREADS()

!$OMP DO SCHEDULE(DYNAMIC)
angleloop: do iang = 1,numangles
    qq = qAR%getQuatfromArray(iang)
    ! qu = qq%get_quatd()

    imageloop: do idir = 1,numk

! do the active coordinate transformation for this euler angle

        dc = klist(1:3,idir)
        dp = DOT_PRODUCT(dc(1:3),(/dsin(enl%sampletilt*dtor),0.D0,dcos(enl%sampletilt*dtor)/))        
        MCangle = acos(dp)*rtod
      
! find index closest to the list of MC runs we already have and interpolate the weight factor
        isig = int(MCangle) + 1
        if (isig .gt. nsig) isig = nsig

        isigp = isig + 1
        if (isigp .gt. nsig) isigp = nsig

        dx = MCangle - int(MCangle)
        dxm =  1.0 - dx
 
        wf = anglewf(isig) * dxm + anglewf(isigp) * dx
        
        dc = qq%quat_LP(dc)
        dc = dc/dsqrt(sum(dc*dc))

! convert these direction cosines to coordinates in the Rosca-Lambert projection
        call LambertgetInterpolation(dc, scl, enl%npix, enl%npix, nix, niy, nixp, niyp, dx, dy, dxm, dym)

! interpolate the intensity
        ipx = kij(1,idir)
        ipy = kij(2,idir)
        
! including the detector model with some sample tilt
        if (switchwfoff .eqv. .FALSE.) then
            if (dc(3) .ge. 0.0) then 
                ECPpattern(ipx,ipy) = wf * ( mLPNH(nix,niy) * dxm * dym + &
                                             mLPNH(nixp,niy) * dx * dym + &
                                             mLPNH(nix,niyp) * dxm * dy + &
                                             mLPNH(nixp,niyp) * dx * dy )
            else
                ECPpattern(ipx,ipy) =  wf * ( mLPSH(nix,niy) * dxm * dym + &
                                              mLPSH(nixp,niy) * dx * dym + &
                                              mLPSH(nix,niyp) * dxm * dy + &
                                              mLPSH(nixp,niyp) * dx * dy )
            end if
        else
            if (dc(3) .ge. 0.0) then 
                ECPpattern(ipx,ipy) =  mLPNH(nix,niy) * dxm * dym + &
                                       mLPNH(nixp,niy) * dx * dym + &
                                       mLPNH(nix,niyp) * dxm * dy + &
                                       mLPNH(nixp,niyp) * dx * dy 
            else
                ECPpattern(ipx,ipy) =  mLPSH(nix,niy) * dxm * dym + &
                                       mLPSH(nixp,niy) * dx * dym + &
                                       mLPSH(nix,niyp) * dxm * dy + &
                                       mLPSH(nixp,niyp) * dx * dy 
            end if
        end if
    end do imageloop
    
    !call BarrelDistortion(D,ECPpattern,enl%npix,enl%npix)
    !ma = maxval(ECPpattern)
    !mi = minval(ECPpattern)
    !ECPpatternintd = ((ECPpattern - mi)/ (ma-mi))
    !ECPpatterninteger = nint(ECPpatternintd*255.0)
    !ECPpatternad =  adhisteq(10,enl%npix,enl%npix,ECPpatterninteger)
    !ECPpattern = float(ECPpatternad)
! =====================================================
! The following write commands constitute HDF_FileVersion = 4.0
! =====================================================

    if (mod(iang,2500) .eq. 0) then
        io_int(1) = iang
        call Message%WriteValue(' completed pattern # ',io_int,1)
    end if
!$OMP CRITICAL
    if (enl%outputformat .eq. 'bin') then
        ma = maxval(ECPpattern)
        mi = minval(ECPpattern)
        ECPpatternintd = ((ECPpattern - mi)/ (ma-mi))
        if (enl%maskpattern.eq.'y')  ECPpatternintd = ECPpatternintd * mask 
        bpat = char(nint(255.0*ECPpatternintd))

! write dictionary pattern to h5 file
        bpat_tmp(:,:,1) = bpat
        offset = (/ 0, 0, iang-1 /)
        hdims = (/ enl%npix, enl%npix, numangles /)
        dims3 = (/ enl%npix, enl%npix, 1 /)
        hdferr = HDF%writeHyperslabCharArray(dataset, bpat_tmp, hdims, offset, dims3, insert)
 
    end if

    if (enl%outputformat .eq. 'gui') then
          if (enl%maskpattern.eq.'y')  ECPpattern = ECPpattern * mask
          ECP_tmp(:,:,1) = ECPpattern
          offset = (/ 0, 0, iang-1 /)
          hdims = (/ enl%npix, enl%npix, numangles /)
          dims3 = (/ enl%npix, enl%npix, 1 /)
          hdferr = HDF%writeHyperslabFloatArray(dataset, ECP_tmp, hdims, offset, dims3, insert)
    end if
!$OMP END CRITICAL

end do angleloop 
!$OMP END DO
!$OMP END PARALLEL

call HDF%pop()
call HDF%pop()

! and update the end time
call timer%makeTimeStamp()
tstre = timer%getTimeString()

groupname = trim(HDFnames%get_EMheader())
hdferr = HDF%openGroup(groupname)

groupname = trim(HDFnames%get_ProgramData())
hdferr = HDF%openGroup(groupname)

! stop time /EMheader/StopTime 'character'
dataset = SC_StopTime
call timer%Time_tock(1) 
tstop = timer%getInterval(1)
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

! close all groups and the file
call HDF%pop(.TRUE.)

! close the Fortran interface
call closeFortranHDFInterface()

call Message%printMessage(' -> Execution completed.',frm='(A)')

end associate 

end subroutine ECP_

!--------------------------------------------------------------------------
recursive subroutine ECPGenerateDetector_(self, verbose)
!! author: MDG 
!! version: 1.0 
!! date: 03/15/20
!!
!! generate the ECP detector arrays

use mod_io
use mod_quaternions
use mod_rotations

IMPLICIT NONE

class(ECP_T),INTENT(INOUT)     :: self
logical, INTENT(IN), OPTIONAL  :: verbose

type(IO_T)                     :: Message 
real(kind=sgl)                 :: thetain, thetaout, polar, azimuthal, delpolar, delazimuth
real(kind=sgl)                 :: io_real(2), om(3,3), sampletilt, dc(3)
integer(kind=irg)              :: iazimuth, ipolar, istat

associate(ecpnl => self%nml, det => self%det)

if (ecpnl%Rin .gt. ecpnl%Rout) then
    call Message%printError('ECPGenerateDetector','Inner radius of annular detector cannot be greater than outer radius')
end if

thetain = atan2(ecpnl%Rin,ecpnl%workingdistance)
thetaout = atan2(ecpnl%Rout,ecpnl%workingdistance)

sampletilt = ecpnl%sampletilt*cPi/180.0
om(1,:) = (/cos(sampletilt),0.0,sin(sampletilt)/)
om(2,:) = (/0.0,1.0,0.0/)
om(3,:) = (/-sin(sampletilt),0.0,cos(sampletilt)/)

if (present(verbose)) then
    if(verbose) then
       io_real(1) = thetain*rtod
       io_real(2) = thetaout*rtod
       call Message%WriteValue('Inner and outer polar angles for detector (in degrees) are ',io_real,2)
    end if
end if

det%npolar = nint((thetaout - thetain)*rtod) + 1
delpolar = (thetaout - thetain)/float(det%npolar-1)

det%nazimuth = 361
delazimuth = 2.0*cPi/float(det%nazimuth-1)

write (*,*) det%npolar, det%nazimuth 

allocate(det%rgx(det%npolar,det%nazimuth))
allocate(det%rgy(det%npolar,det%nazimuth))
allocate(det%rgz(det%npolar,det%nazimuth),stat=istat)
if (istat .ne. 0) then 
  write (*,*) 'istat = ', istat
  call Message%printError('ECPGenerateDetector','cannot allocate the rgx, rgy and rgz arrays')
end if

det%rgx = 0.0
det%rgy = 0.0
det%rgz = 0.0

! compute the direction cosines of the detector elements in the sample reference frame.
do ipolar = 1,det%npolar
    polar = thetain + float(ipolar-1)*delpolar

    do iazimuth = 1,det%nazimuth
         azimuthal = float(iazimuth-1)*delazimuth

         dc(1) = cos(azimuthal)*sin(polar)
         dc(2) = sin(azimuthal)*sin(polar)
         dc(3) = cos(polar)

         dc = matmul(om,dc)

         det%rgx(ipolar,iazimuth) = dc(1)
         det%rgy(ipolar,iazimuth) = dc(2)
         det%rgz(ipolar,iazimuth) = dc(3)
    end do
end do

if (present(verbose)) then
    if(verbose) then
        call Message%printMessage(' -> Finished generating detector',frm='(A)')
    end if
end if

end associate 

end subroutine ECPGenerateDetector_

!--------------------------------------------------------------------------
recursive subroutine ECPGetWeightFactors_(self, mcnl, MCFT, weightfact, nsig, verbose)
!! author: MDG 
!! version: 1.0 
!! date: 03/15/20
!!
!! generate the ECP weight factor array

use mod_io
use mod_quaternions
use mod_rotations
use mod_MCfiles
use mod_Lambert

IMPLICIT NONE

class(ECP_T),INTENT(INOUT)              :: self
type(MCOpenCLNameListType),INTENT(INOUT):: mcnl
type(MCfile_T),INTENT(INOUT)            :: MCFT
integer(kind=irg), INTENT(IN)           :: nsig
real(kind=sgl), INTENT(OUT)             :: weightfact(nsig)
logical, INTENT(IN), OPTIONAL           :: verbose

type(IO_T)                              :: Message
integer(kind=irg)                       :: isig, ipolar, iazimuth, istat
integer(kind=irg)                       :: nix, niy, nixp, niyp, isampletilt
real(kind=sgl)                          :: dx, dy, dxm, dym, acc_sum, samplenormal(3), dp
real(kind=sgl)                          :: dc(3), ixy(2), scl, deltheta, thetac, x, MCangle
integer(kind=irg),allocatable           :: acc(:,:,:)

associate( ecpnl => self%nml, det => self%det )

scl = mcnl%numsx
call MCFT%copyaccume(acc)

thetac = ecpnl%conesemiangle
deltheta = (thetac+abs(ecpnl%sampletilt))/float(nsig-1)

weightfact = 0.0

do isig = 1,nsig
    acc_sum = 0.0
    MCangle = (isig - 1)*deltheta
    isampletilt = nint((MCangle - mcnl%sigstart)/mcnl%sigstep)

    if (isampletilt .lt. 1) then
        isampletilt = abs(isampletilt) + 1
    else
        isampletilt = isampletilt + 1
    end if

    do ipolar = 1,det%npolar
        do iazimuth = 1,det%nazimuth
            dc(1:3) = (/det%rgx(ipolar,iazimuth),det%rgy(ipolar,iazimuth),det%rgz(ipolar,iazimuth)/)

! convert to Rosca-lambert projection
            call LambertgetInterpolation(dc, scl, mcnl%numsx, mcnl%numsx, nix, niy, nixp, niyp, dx, dy, dxm, dym)

            acc_sum = 0.25*(acc(isampletilt,nix,niy) * dxm * dym + &
                            acc(isampletilt,nixp,niy) * dx * dym + &
                            acc(isampletilt,nix,niyp) * dxm * dy + &
                            acc(isampletilt,nixp,niyp) * dx * dy)

            weightfact(isig) = weightfact(isig) + acc_sum

        end do
    end do
end do

weightfact(1:nsig) = weightfact(1:nsig)/weightfact(1)

if (present(verbose)) then
    if (verbose) call Message%printMessage(' -> Finished calculating the weight factors',frm='(A)')
end if

end associate

end subroutine ECPGetWeightFactors_

!--------------------------------------------------------------------------
recursive subroutine GetVectorsCone_(self)
!! author: MDG 
!! version: 1.0 
!! date: 03/15/20
!!
!! generate incident wave vectors inside the cone

IMPLICIT NONE 

class(ECP_T),INTENT(INOUT)          :: self

type(IncidentListECP),pointer       :: ktmp
real(kind=dbl)                      :: kk(3), thetacr, delta, ktmax
integer(kind=irg)                   :: imin, imax, jmin, jmax, numk
integer(kind=irg)                   :: ii, jj, istat

associate( ecpnl => self%nml, klist => self%klist )

numk = 0
kk = (/0.D0,0.D0,1.D0/)
thetacr = dtor*ecpnl%conesemiangle
ktmax = tan(thetacr)
delta = 2.0*ktmax/(float(ecpnl%npix)-1.0)

imin = 1
imax = ecpnl%npix
jmin = 1
jmax = ecpnl%npix

ktmp => self%get_ListHead()

do ii = imin, imax
    do jj = jmin, jmax
        ktmp%k(1:3) = (/-ktmax+delta*(ii-1),-ktmax+delta*(jj-1),0.D0/) + kk(1:3)
        ktmp%k = ktmp%k/sqrt(sum(ktmp%k**2))
        ktmp%i = ii
        ktmp%j = jj
        numk = numk + 1
        allocate(ktmp%next)
        ktmp => ktmp%next
        nullify(ktmp%next)
    end do
end do

end associate 

call self%set_numk(numk) 

end subroutine GetVectorsCone_

!--------------------------------------------------------------------------
recursive subroutine CalcECPatternSingle_(self, ipar, qu, anglewf, mLPNH, mLPSH, kij, klist, ECPattern)
!! author: MDG 
!! version: 1.0 
!! date: 03/15/20
!!
!! Calculate a single EC pattern for a given orientation 

use mod_io 
use mod_quaternions 
use mod_Lambert 

IMPLICIT NONE

class(ECP_T), INTENT(INOUT)                     :: self
integer(kind=irg), INTENT(IN)                   :: ipar(4)
type(Quaternion_T),INTENT(INOUT)                :: qu
real(kind=sgl),INTENT(IN)                       :: anglewf(ipar(1))
real(kind=sgl),INTENT(IN)                       :: mLPNH(-ipar(4):ipar(4),-ipar(4):ipar(4))
real(kind=sgl),INTENT(IN)                       :: mLPSH(-ipar(4):ipar(4),-ipar(4):ipar(4))
integer(kind=irg),INTENT(IN)                    :: kij(2,ipar(2))
real(kind=sgl),INTENT(IN)                       :: klist(3,ipar(2))
real(kind=sgl),INTENT(OUT)                      :: ECPattern(1:ipar(3),1:ipar(3))

integer(kind=irg)                               :: numk, idir, isig, isigp, nsig, istat
real(kind=dbl)                                  :: dc(3), dc2(3), dp, MCangle, scl, ixy(2)
real(kind=dbl)                                  :: wf, dx, dy, dxm, dym
integer(kind=irg)                               :: nix, niy, nixp, niyp, ipx, ipy

associate(ecpnl=>self%nml)

numk = ipar(2)
nsig = ipar(1)
scl = dble(ipar(3))

do idir = 1,numk

! do the active coordinate transformation for this euler angle
    dc = klist(1:3,idir)
    dp = DOT_PRODUCT(dc(1:3),(/dsin(ecpnl%sampletilt*dtor),0.D0,dcos(ecpnl%sampletilt*dtor)/))        
      

    MCangle = acos(dp)*rtod
! find index closest to the list of MC runs we already have and interpolate the weight factor
    isig = int(MCangle) + 1
    if (isig .gt. nsig) isig = nsig

    isigp = isig + 1
    if (isigp .gt. nsig) isigp = nsig

    dx = MCangle - int(MCangle)
    dxm =  1.0 - dx
 
    wf = anglewf(isig) * dxm + anglewf(isigp) * dx

    dc2 = qu%quat_Lp(dc)
    dc = dc2/sqrt(sum(dc2*dc2))

! convert these direction cosines to coordinates in the Rosca-Lambert projection
    call LambertgetInterpolation(dc, scl, ecpnl%npix, ecpnl%npix, nix, niy, nixp, niyp, dx, dy, dxm, dym)

! interpolate the intensity
    ipx = kij(1,idir)
    ipy = kij(2,idir)
        
! including the detector model with some sample tilt
    if (dc(3) .gt. 0.0) then 

        ECPattern(ipx,ipy) = wf * (mLPNH(nix,niy) * dxm * dym + &
                                   mLPNH(nixp,niy) * dx * dym + &
                                   mLPNH(nix,niyp) * dxm * dy + &
                                   mLPNH(nixp,niyp) * dx * dy)

    else

        ECPattern(ipx,ipy) = wf * (mLPSH(nix,niy) * dxm * dym + &
                                   mLPSH(nixp,niy) * dx * dym + &
                                   mLPSH(nix,niyp) * dxm * dy + &
                                   mLPSH(nixp,niyp) * dx * dy)

    end if

end do 

end associate 

end subroutine CalcECPatternSingle_




end module mod_ECP
