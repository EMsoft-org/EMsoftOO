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

module mod_EBSD
  !! author: MDG
  !! version: 1.0
  !! date: 02/17/20
  !!
  !! class definition for the EMEBSD program

use mod_kinds
use mod_global
use mod_quaternions
use mod_MPfiles

IMPLICIT NONE
private

! namelist for the EMEBSD program
type, public :: EBSDNameListType
  integer(kind=irg)       :: numsx
  integer(kind=irg)       :: numsy
  integer(kind=irg)       :: binning
  integer(kind=irg)       :: nthreads
  ! integer(kind=irg)       :: energyaverage    [removed option, 02/17/20, MDG]
  integer(kind=irg)       :: maskradius
  integer(kind=irg)       :: nregions
  real(kind=sgl)          :: L
  real(kind=sgl)          :: thetac
  real(kind=sgl)          :: delta
  real(kind=sgl)          :: omega
  real(kind=sgl)          :: xpc
  real(kind=sgl)          :: ypc
  real(kind=sgl)          :: energymin
  real(kind=sgl)          :: energymax
  real(kind=sgl)          :: gammavalue
  real(kind=sgl)          :: axisangle(4)
  real(kind=sgl)          :: alphaBD
  real(kind=sgl)          :: hipassw
  real(kind=dbl)          :: Ftensor(3,3)
  real(kind=dbl)          :: beamcurrent
  real(kind=dbl)          :: dwelltime
  character(1)            :: makedictionary
  character(1)            :: poisson
  character(1)            :: includebackground
  character(1)            :: applyDeformation
  character(1)            :: maskpattern
  character(3)            :: scalingmode
  character(3)            :: eulerconvention
  character(3)            :: outputformat
  character(1)            :: spatialaverage
  character(5)            :: bitdepth
  character(fnlen)        :: anglefile
  character(fnlen)        :: anglefiletype
  character(fnlen)        :: masterfile
  character(fnlen)        :: energyfile
  character(fnlen)        :: datafile
end type EBSDNameListType

! angles, patterns centers, and deformation tensor arrays
type EBSDAnglePCDefType
  real(kind=sgl),allocatable  :: pcs(:,:)
  real(kind=sgl),allocatable  :: deftensors(:,:,:)
  real(kind=dbl),allocatable  :: pcfield(:,:,:)
  real(kind=dbl),allocatable  :: deftensorfield(:,:,:,:)
end type EBSDAnglePCDefType

type EBSDPixel
  real(kind=sgl),allocatable  :: lambdaEZ(:,:)
  real(kind=dbl)              :: dc(3) ! direction cosine in sample frame
  real(kind=dbl)              :: cfactor
end type EBSDPixel

type EBSDDetectorType
  real(kind=sgl),allocatable  :: rgx(:,:), rgy(:,:), rgz(:,:)  ! auxiliary detector arrays needed for interpolation
  real(kind=sgl),allocatable  :: accum_e_detector(:,:,:)
  type(EBSDPixel),allocatable :: detector(:,:)
end type EBSDDetectorType


! class definition
type, public :: EBSD_T
private
  character(fnlen)                :: nmldeffile = 'EMEBSD.nml'
  type(EBSDNameListType), public  :: nml
  type(EBSDDetectorType), public  :: det

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: EBSD_
  procedure, pass(self) :: EBSDreadorpcdef_
  procedure, pass(self) :: GenerateDetector_
  procedure, pass(self) :: GeneratemyEBSDDetector_
  procedure, pass(self) :: ComputeEBSDPatterns_
  procedure, pass(self) :: ComputedeformedEBSDpatterns_
  procedure, pass(self) :: CalcEBSDPatternSingleFull_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: EBSD => EBSD_
  generic, public :: EBSDreadorpcdef => EBSDreadorpcdef_
  generic, public :: GenerateDetector => GenerateDetector_
  generic, public :: GeneratemyEBSDDetector => GeneratemyEBSDDetector_
  generic, public :: ComputeEBSDPatterns => ComputeEBSDPatterns_
  generic, public :: ComputedeformedEBSDpatterns => ComputedeformedEBSDpatterns_
  generic, public :: CalcEBSDPatternSingleFull => CalcEBSDPatternSingleFull_

end type EBSD_T

! the constructor routine for this class
interface EBSD_T
  module procedure EBSD_constructor
end interface EBSD_T

contains

!--------------------------------------------------------------------------
type(EBSD_T) function EBSD_constructor( nmlfile, isTKD ) result(EBSD)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSD_constructor
!! author: MDG
!! version: 1.0
!! date: 02/17/20
!!
!! constructor for the EBSD_T Class; reads the name list

IMPLICIT NONE

character(fnlen)    :: nmlfile
logical, OPTIONAL   :: isTKD

if (present(isTKD)) then
  call EBSD%readNameList(nmlfile, isTKD=isTKD)
else
  call EBSD%readNameList(nmlfile)
end if

end function EBSD_constructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly, isTKD)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 02/17/20
!!
!! read the namelist from an nml file for the EBSD_T Class

use mod_io

IMPLICIT NONE

class(EBSD_T), INTENT(INOUT)                :: self
character(fnlen),INTENT(IN)                 :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)                 :: initonly
 !! fill in the default values only; do not read the file
logical,OPTIONAL,INTENT(IN)                 :: isTKD

type(IO_T)                                  :: Message
logical                                     :: skipread = .FALSE.

integer(kind=irg)       :: numsx
integer(kind=irg)       :: numsy
integer(kind=irg)       :: binning
integer(kind=irg)       :: nthreads
integer(kind=irg)       :: maskradius
integer(kind=irg)       :: nregions
real(kind=sgl)          :: L
real(kind=sgl)          :: thetac
real(kind=sgl)          :: delta
real(kind=sgl)          :: xpc
real(kind=sgl)          :: ypc
real(kind=sgl)          :: omega
real(kind=sgl)          :: energymin
real(kind=sgl)          :: energymax
real(kind=sgl)          :: gammavalue
real(kind=sgl)          :: alphaBD
real(kind=sgl)          :: axisangle(4)
real(kind=sgl)          :: hipassw
real(kind=dbl)          :: Ftensor(3,3)
real(kind=dbl)          :: beamcurrent
real(kind=dbl)          :: dwelltime
character(1)            :: includebackground
character(1)            :: poisson
character(1)            :: makedictionary
character(1)            :: applyDeformation
character(1)            :: maskpattern
character(1)            :: spatialaverage
character(3)            :: scalingmode
character(3)            :: eulerconvention
character(3)            :: outputformat
character(5)            :: bitdepth
character(fnlen)        :: anglefile
character(fnlen)        :: anglefiletype
character(fnlen)        :: masterfile
character(fnlen)        :: energyfile  ! removed from template file 05/16/19 [MDG]
character(fnlen)        :: datafile

! define the IO namelist to facilitate passing variables to the program.
namelist  / EBSDdata / L, thetac, delta, numsx, numsy, xpc, ypc, anglefile, eulerconvention, masterfile, bitdepth, &
                       energyfile, datafile, beamcurrent, dwelltime, energymin, energymax, binning, gammavalue, alphaBD, &
                       scalingmode, axisangle, nthreads, outputformat, maskpattern, omega, spatialaverage, &
                       applyDeformation, Ftensor, includebackground, anglefiletype, makedictionary, hipassw, nregions, &
                       maskradius, poisson

! define the IO namelist to facilitate passing variables to the program.
! the two name lists are currently identical, but we allow for the fact that
! they might diverge in the future.
namelist  / TKDdata /  L, thetac, delta, numsx, numsy, xpc, ypc, anglefile, eulerconvention, masterfile, bitdepth, &
                       energyfile, datafile, beamcurrent, dwelltime, energymin, energymax, binning, gammavalue, alphaBD, &
                       scalingmode, axisangle, nthreads, outputformat, maskpattern, omega, spatialaverage, &
                       applyDeformation, Ftensor, includebackground, anglefiletype, makedictionary, hipassw, nregions, &
                       maskradius, poisson

! set the input parameters to default values (except for xtalname, which must be present)
numsx           = 0             ! [dimensionless]
numsy           = 0             ! [dimensionless]
binning         = 1             ! binning mode  (1, 2, 4, or 8)
L               = 20000.0       ! [microns]
nthreads        = 1             ! number of OpenMP threads
nregions        = 10            ! number of regions in adaptive histogram equalization
thetac          = 0.0           ! [degrees]
delta           = 25.0          ! [microns]
xpc             = 0.0           ! [pixels]
ypc             = 0.0           ! [pixels]
omega           = 0.0
energymin       = 15.0          ! minimum energy to consider
energymax       = 30.0          ! maximum energy to consider
gammavalue      = 1.0           ! gamma factor
alphaBD         = 0.0           ! transfer lens barrel distortion parameter
maskradius      = 240           ! mask radius
hipassw         = 0.05          ! hi-pass filter radius
axisangle       = (/0.0, 0.0, 1.0, 0.0/)        ! no additional axis angle rotation
Ftensor         = reshape( (/ 1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 1.D0 /), (/ 3,3 /) )
beamcurrent     = 14.513D0      ! beam current (actually emission current) in nano ampere
dwelltime       = 100.0D0       ! in microseconds
makedictionary  = 'y'
poisson         = 'n'           ! apply poisson noise ?
includebackground = 'y'         ! set to 'n' to remove realistic background intensity profile
applyDeformation = 'n'          ! should we apply a deformation tensor to the unit cell?
maskpattern     = 'n'           ! 'y' or 'n' to include a circular mask
scalingmode     = 'not'         ! intensity selector ('lin', 'gam', or 'not')
eulerconvention = 'tsl'         ! convention for the first Euler angle ['tsl' or 'hkl']
outputformat    = 'gui'         ! output format for 'bin' or 'gui' use
bitdepth        = '8bit'        ! format for output; '8char' for [0..255], '##int' for integers, 'float' for floats
! the '##int' notation stands for the actual bitdepth; all values are stored as 32bit integers, but they are scaled
! from the float values to a maximum that is given by the first two digits, which indicate the bit depth; so, valid
! values would be '10int' for a 10-bit integer scale, '16int' for a 16-bit integer scale, and so on.
anglefile       = 'undefined'   ! filename
anglefiletype   = 'orientations'! 'orientations' or 'orpcdef'
masterfile      = 'undefined'   ! filename
energyfile      = 'undefined'   ! name of file that contains energy histograms for all scintillator pixels (output from MC program)
datafile        = 'undefined'   ! output file name
spatialaverage  = 'n'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 if (present(isTKD)) then
   read(UNIT=dataunit,NML=TKDdata)
 else
   read(UNIT=dataunit,NML=EBSDdata)
 end if
 close(UNIT=dataunit,STATUS='keep')

! check for required entries

! we no longer require the energyfile parameter, but for backwards compatibility
! we still allow the user to include it (it doesn't do anything though)
! if (trim(energyfile).eq.'undefined') then
!  call FatalError('GetEBSDNameList:',' energy file name is undefined in '//nmlfile)
! end if

 if (trim(anglefile).eq.'undefined') then
  call Message%printError('readNameList:',' angle file name is undefined in '//nmlfile)
 end if

 if (trim(masterfile).eq.'undefined') then
  call Message%printError('readNameList:',' master pattern file name is undefined in '//nmlfile)
 end if

 if (trim(datafile).eq.'undefined') then
  call Message%printError('readNameList:',' output file name is undefined in '//nmlfile)
 end if

 if (numsx.eq.0) then
  call Message%printError('readNameList:',' pattern size numsx is zero '//nmlfile)
 end if

 if (numsx.eq.0) then
  call Message%printError('readNameList:',' pattern size numsy is zero '//nmlfile)
 end if
end if

! if we get here, then all appears to be ok, and we need to fill in the enl fields
self%nml%numsx = numsx
self%nml%numsy = numsy
self%nml%binning = binning
self%nml%nregions = nregions
self%nml%maskradius = maskradius
self%nml%L = L
self%nml%nthreads = nthreads
self%nml%thetac = thetac
self%nml%delta = delta
self%nml%xpc = xpc
self%nml%ypc = ypc
self%nml%energymin = energymin
self%nml%energymax = energymax
self%nml%gammavalue = gammavalue
self%nml%alphaBD = alphaBD
self%nml%hipassw = hipassw
self%nml%axisangle = axisangle
self%nml%Ftensor = Ftensor
self%nml%beamcurrent = beamcurrent
self%nml%dwelltime = dwelltime
self%nml%includebackground = includebackground
self%nml%makedictionary = makedictionary
self%nml%poisson = poisson
self%nml%applyDeformation = applyDeformation
self%nml%maskpattern = maskpattern
self%nml%scalingmode = scalingmode
self%nml%eulerconvention = eulerconvention
self%nml%outputformat = outputformat
self%nml%bitdepth = bitdepth
self%nml%anglefile = anglefile
self%nml%anglefiletype = anglefiletype
self%nml%masterfile = masterfile
! we require energyfile to be identical to masterfile, so the
! user definition, if any, in the namelist file is overwritten here...
self%nml%energyfile = masterfile       ! changed on 05/16/19 [MDG]
self%nml%datafile = datafile
self%nml%omega = omega
self%nml%spatialaverage = spatialaverage

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 02/17/20
!!
!! pass the namelist for the EBSD_T Class to the calling program

IMPLICIT NONE

class(EBSD_T), INTENT(INOUT)          :: self
type(EBSDNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames, isTKD)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG
!! version: 1.0
!! date: 02/17/20
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants

use ISO_C_BINDING

IMPLICIT NONE

class(EBSD_T), INTENT(INOUT)            :: self
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames
logical,OPTIONAL,INTENT(IN)             :: isTKD

integer(kind=irg),parameter             :: n_int = 6, n_real = 10
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
io_int = (/ enl%numsx, enl%numsy, enl%binning, enl%nthreads, enl%nregions, enl%maskradius /)
intlist(1) = 'numsx'
intlist(2) = 'numsy'
intlist(3) = 'binning'
intlist(4) = 'nthreads'
intlist(5) = 'nregions'
intlist(6) = 'maskradius'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ enl%L, enl%thetac, enl%delta, enl%xpc, enl%ypc, enl%energymin, enl%energymax, enl%gammavalue, &
             enl%alphaBD, enl%hipassw /)
reallist(1) = 'L'
reallist(2) = 'thetac'
reallist(3) = 'delta'
reallist(4) = 'xpc'
reallist(5) = 'ypc'
reallist(6) = 'energymin'
reallist(7) = 'energymax'
reallist(8) = 'gammavalue'
reallist(9) = 'alphaBD'
reallist(10)= 'hipassw'
call HDF%writeNMLreals(io_real, reallist, n_real)

! a 4-vector
dataset = SC_axisangle
hdferr = HDF%writeDatasetFloatArray(dataset, enl%axisangle, 4)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create axisangle dataset', hdferr)

! a 3x3 matrix
dataset = SC_Ftensor
hdferr = HDF%writeDatasetDoubleArray(dataset, enl%Ftensor, 3, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create Ftensor dataset', hdferr)

! a few doubles
dataset = SC_beamcurrent
hdferr = HDF%writeDatasetDouble(dataset, enl%beamcurrent)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create beamcurrent dataset', hdferr)

dataset = SC_dwelltime
hdferr = HDF%writeDatasetDouble(dataset, enl%dwelltime)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create dwelltime dataset', hdferr)

! write all the strings
dataset = SC_maskpattern
line2(1) = trim(enl%maskpattern)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create maskpattern dataset', hdferr)

dataset = SC_makedictionary
line2(1) = trim(enl%makedictionary)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create makedictionary dataset', hdferr)

dataset = SC_poisson
line2(1) = trim(enl%poisson)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create poisson dataset', hdferr)

dataset = SC_includebackground
line2(1) = trim(enl%includebackground)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create includebackground dataset', hdferr)

dataset = SC_applyDeformation
line2(1) = trim(enl%applyDeformation)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create applyDeformation dataset', hdferr)

dataset = SC_scalingmode
line2(1) = trim(enl%scalingmode)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create scalingmode dataset', hdferr)

dataset = SC_eulerconvention
line2(1) = trim(enl%eulerconvention)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create eulerconvention dataset', hdferr)

dataset = SC_outputformat
line2(1) = trim(enl%outputformat)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create outputformat dataset', hdferr)

dataset = SC_bitdepth
line2(1) = trim(enl%bitdepth)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create bitdepth dataset', hdferr)

dataset = SC_energyfile
line2(1) = trim(enl%energyfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create energyfile dataset', hdferr)

dataset = SC_masterfile
line2(1) = trim(enl%masterfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create masterfile dataset', hdferr)

dataset = SC_anglefile
line2(1) = trim(enl%anglefile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create anglefile dataset', hdferr)

dataset = SC_anglefiletype
line2(1) = trim(enl%anglefiletype)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create anglefiletype dataset', hdferr)

dataset = SC_datafile
line2(1) = trim(enl%datafile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create datafile dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine EBSD_(self, EMsoft, progname, HDFnames, TKD)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSD_
!! author: MDG
!! version: 1.0
!! date: 02/17/20
!!
!! perform the computations

use mod_EMsoft
use mod_so3
use mod_quaternions
use mod_MCfiles
use mod_MPfiles
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_io
use mod_rotations
use stringconstants

IMPLICIT NONE

class(EBSD_T), INTENT(INOUT)        :: self
type(EMsoft_T), INTENT(INOUT)       :: EMsoft
character(fnlen), INTENT(INOUT)     :: progname
type(HDFnames_T), INTENT(INOUT)     :: HDFnames
logical, OPTIONAL, INTENT(IN)       :: TKD

type(MCfile_T)                      :: MCFT
type(MPfile_T)                      :: MPFT
type(HDF_T)                         :: HDF
type(HDFnames_T)                    :: saveHDFnames
type(so3_T)                         :: SO
type(IO_T)                          :: Message
type(Quaternion_T)                  :: quat
type(QuaternionArray_T)             :: qAR

type(EBSDmasterNameListType)        :: mpnl
type(MCOpenCLNameListType)          :: mcnl
type(EBSDAnglePCDefType)            :: orpcdef

logical                             :: verbose, isTKD = .FALSE.
character(fnlen)                    :: fname, nmldeffile
integer(kind=irg)                   :: numangles, istat
type(FZpointd),pointer              :: FZtmp
type(r_T)                           :: rr

if (present(TKD)) then
  saveHDFnames = HDFnames
  isTKD = .TRUE.
  call MPFT%setModality('TKD')
else
  call MPFT%setModality('EBSD')
end if
nmldeffile = trim(EMsoft%nmldeffile)

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()

associate( enl => self%nml, EBSDdetector => self%det, EBSDMCdata => MCFT%MCDT )

! set the HDF group names for this program
HDFnames = HDFnames_T()

! 1. read the angle array from file
verbose = .TRUE.

call setRotationPrecision('d')

if (trim(enl%anglefiletype).eq.'orientations') then
  fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%anglefile))
  call SO%nullifyList()
  call SO%getOrientationsfromFile(fname)
  numangles = SO%getListCount('FZ')
  call SO%listtoQuaternionArray( qAR )
  call SO%delete_FZlist()
else if (trim(enl%anglefiletype).eq.'orpcdef') then
! this requires a conversion from the Euler angles in the file to quaternions
! plus storage of the pattern center and deformation tensor arrays
  call self%EBSDreadorpcdef(EMsoft, numangles, qAR, orpcdef, verbose)
else
  call Message%printError('EBSD','unknown anglefiletype')
end if

! 2. read the Monte Carlo data file (HDF format)
call HDFnames%set_ProgramData(SC_MCOpenCL)
call HDFnames%set_NMLlist(SC_MCCLNameList)
call HDFnames%set_NMLfilename(SC_MCOpenCLNML)
fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%energyfile))
call MCFT%setFileName(fname)
call MCFT%readMCfile(HDF, HDFnames, getAccume=.TRUE.)
mcnl = MCFT%getnml()

! 3. read EBSD master pattern file (HDF format)
if (isTKD.eqv..TRUE.) then
  call HDFnames%set_ProgramData(SC_TKDmaster)
  call HDFnames%set_NMLlist(SC_TKDmasterNameList)
  call HDFnames%set_NMLfilename(SC_TKDmasterNML)
else
  call HDFnames%set_ProgramData(SC_EBSDmaster)
  call HDFnames%set_NMLlist(SC_EBSDmasterNameList)
  call HDFnames%set_NMLfilename(SC_EBSDmasterNML)
end if
call HDFnames%set_Variable(SC_MCOpenCL)

fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfile))
call MPFT%setFileName(fname)
call MPFT%readMPfile(HDF, HDFnames, mpnl, getmLPNH=.TRUE., getmLPSH=.TRUE.)

if (isTKD.eqv..TRUE.) then
  HDFnames = saveHDFnames
else
  call HDFnames%set_ProgramData(SC_EBSD)
  call HDFnames%set_NMLlist(SC_EBSDNameList)
  call HDFnames%set_NMLfilename(SC_EBSD)
end if

! for a regular Euler angle file, we precompute the detector arrays here; for the 'orpcdef' mode
! we compute them later (for each pattern separately)
if (trim(enl%anglefiletype).eq.'orientations') then
  allocate(EBSDdetector%rgx(enl%numsx,enl%numsy), &
           EBSDdetector%rgy(enl%numsx,enl%numsy), &
           EBSDdetector%rgz(enl%numsx,enl%numsy), &
           EBSDdetector%accum_e_detector(EBSDMCdata%numEbins,enl%numsx,enl%numsy), stat=istat)
! 4. generate detector arrays
  call self%GenerateDetector(MCFT, verbose, isTKD)

  ! perform the pattern computations
  call self%ComputeEBSDPatterns(EMsoft, MCFT, MPFT, HDF, HDFnames, mpnl, numangles, qAR, progname, nmldeffile)
end if

if (trim(enl%anglefiletype).eq.'orpcdef') then
  call self%ComputedeformedEBSDPatterns(EMsoft, MCFT, MPFT, HDF, HDFnames, mpnl, numangles, qAR, orpcdef, progname, nmldeffile)
end if

end associate

call closeFortranHDFInterface()

end subroutine EBSD_

!--------------------------------------------------------------------------
recursive subroutine EBSDreadorpcdef_(self, EMsoft, numangles, qAR, orpcdef, verbose)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDreadorpcdef_
!! author: MDG
!! version: 1.0
!! date: 02/17/20
!!
!! read angles, pattern centers, and deformation tensors from a text file

use mod_EMsoft
use mod_io
use mod_rotations
use mod_quaternions

IMPLICIT NONE

class(EBSD_T), INTENT(INOUT)            :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
integer(kind=irg),INTENT(OUT)           :: numangles
type(QuaternionArray_T), INTENT(INOUT)  :: qAR
type(EBSDAnglePCDefType),INTENT(INOUT)  :: orpcdef
logical,INTENT(IN),OPTIONAL             :: verbose

type(IO_T)                              :: Message
type(e_T)                               :: e
type(q_T)                               :: q
type(Quaternion_T)                      :: qq

integer(kind=irg)                       :: io_int(1), i, istat
character(2)                            :: atype
real(kind=sgl)                          :: eulang(3)
character(fnlen)                        :: fname

associate( nml => self%nml )

!====================================
! get the angular information, either in Euler angles or in quaternions, from a text file
!====================================
! open the angle file
  fname = EMsoft%generateFilePath('EMdatapathname',nml%anglefile)
  open(unit=dataunit,file=trim(fname),status='old',action='read')

! get the type of angle first [ 'eu' or 'qu' ]
  read(dataunit,*) atype
  if (atype.ne.'eu') then
    call Message%printError("EBSDreadorpcdef","Other orientation formats to be implemented; only Euler for now")
  end if

! then the number of angles in the file
  read(dataunit,*) numangles

  if (present(verbose)) then
    io_int(1) = numangles
    call Message%WriteValue(' Number of angle entries = ',io_int,1)
  end if

! allocate the euler angle, pattern center, and deformation tensor arrays
  allocate(orpcdef%pcs(3,numangles),stat=istat)
  allocate(orpcdef%deftensors(3,3,numangles),stat=istat)
  qAR = QuaternionArray_T( n = numangles, s='d' )

! if istat.ne.0 then do some error handling ...
  do i=1,numangles
    read(dataunit,*) eulang(1:3), orpcdef%pcs(1:3,i), orpcdef%deftensors(1:3,1:3,i)
    if (nml%eulerconvention.eq.'hkl') eulang(1) = eulang(1) + 90.0
! insert the angles into the qAR quaternion array
    call e%e_setd( eulang*dtor )
    q = e%eq()
    qq = Quaternion_T( qd = q%q_copyd() )
    call qAR%insertQuatinArray( i, qq )
  end do
  close(unit=dataunit,status='keep')

! convert the euler angle triplets to quaternions
  if (present(verbose)) call Message%printMessage('  -> converting Euler angles to quaternions', frm = "(A/)")

  call Message%printMessage(' Completed reading Euler angles, pattern centers, and deformation tensors')

end associate

end subroutine EBSDreadorpcdef_

!--------------------------------------------------------------------------
recursive subroutine GenerateDetector_(self, MCFT, verbose, isTKD)
!DEC$ ATTRIBUTES DLLEXPORT :: GenerateDetector_
!! author: MDG
!! version: 1.0
!! date: 02/05/20
!!
!! generate the detector arrays for both EBSD and TKD modalities

use mod_io
use mod_Lambert
use mod_math
use mod_MCfiles

IMPLICIT NONE

class(EBSD_T), INTENT(INOUT)            :: self
type(MCfile_T), INTENT(INOUT)           :: MCFT
logical,INTENT(IN),OPTIONAL             :: verbose
logical,INTENT(IN),OPTIONAL             :: isTKD

type(IO_T)                              :: Message
type(Lambert_T)                         :: L

real(kind=sgl),allocatable              :: scin_x(:), scin_y(:), testarray(:,:)  ! scintillator coordinate arrays [microns]
real(kind=sgl)                          :: alp, ca, sa, cw, sw
real(kind=sgl)                          :: L2, Ls, Lc, calpha     ! distances
real(kind=sgl),allocatable              :: z(:,:)
integer(kind=irg)                       :: nix, niy, binx, biny , i, j, Emin, Emax, istat, k, ipx, ipy, nsx, nsy, elp  ! various parameters
real(kind=sgl)                          :: dc(3), scl, alpha, theta, g, pcvec(3), s, dp           ! direction cosine array
real(kind=sgl)                          :: sx, dx, dxm, dy, dym, rhos, x, bindx         ! various parameters
real(kind=sgl)                          :: ixy(2)
logical                                 :: TKD = .FALSE.

associate( enl => self%nml, mcnl => MCFT%nml, EBSDMCdata => MCFT%MCDT, EBSDdetector => self%det )

if (present(isTKD)) TKD=.TRUE.

!====================================
! ------ generate the detector arrays
!====================================
! This needs to be done only once for a given detector geometry
allocate(scin_x(enl%numsx),scin_y(enl%numsy),stat=istat)

! change to the detector point of view necessitates negating the x pattern center coordinate
scin_x = - ( -enl%xpc - ( 1.0 - enl%numsx ) * 0.5 - (/ (i-1, i=1,enl%numsx) /) ) * enl%delta
scin_y = ( enl%ypc - ( 1.0 - enl%numsy ) * 0.5 - (/ (i-1, i=1,enl%numsy) /) ) * enl%delta

! auxiliary angle to rotate between reference frames
alp = 0.5 * cPi - (mcnl%sig - enl%thetac) * dtor
ca = cos(alp)
sa = sin(alp)

cw = cos(mcnl%omega * dtor)
sw = sin(mcnl%omega * dtor)

! we will need to incorporate a series of possible distortions
! here as well, as described in Gert Nolze's paper; for now we
! just leave this place holder comment instead

! compute auxilliary interpolation arrays
! if (istat.ne.0) then ...

elp = enl%numsy + 1
L2 = enl%L * enl%L
do j=1,enl%numsx
  sx = L2 + scin_x(j) * scin_x(j)
  Ls = -sw * scin_x(j) + enl%L*cw
  Lc = cw * scin_x(j) + enl%L*sw
  do i=1,enl%numsy
   rhos = 1.0/sqrt(sx + scin_y(i)**2)
   EBSDdetector%rgx(j,elp-i) = (scin_y(i) * ca + sa * Ls) * rhos!Ls * rhos
   EBSDdetector%rgy(j,elp-i) = Lc * rhos!(scin_x(i) * cw + Lc * sw) * rhos
   EBSDdetector%rgz(j,elp-i) = (-sa * scin_y(i) + ca * Ls) * rhos!(-sw * scin_x(i) + Lc * cw) * rhos
  end do
end do
deallocate(scin_x, scin_y)

! normalize the direction cosines.
allocate(z(enl%numsx,enl%numsy))
  z = 1.0/sqrt(EBSDdetector%rgx*EBSDdetector%rgx+EBSDdetector%rgy*EBSDdetector%rgy+EBSDdetector%rgz*EBSDdetector%rgz)
  EBSDdetector%rgx = EBSDdetector%rgx*z
  EBSDdetector%rgy = EBSDdetector%rgy*z
  EBSDdetector%rgz = EBSDdetector%rgz*z
deallocate(z)
!====================================

!====================================
! ------ create the equivalent detector energy array
!====================================
! from the Monte Carlo energy data, we need to extract the relevant
! entries for the detector geometry defined above.  Once that is
! done, we can get rid of the larger energy array
!
! in the old version, we either computed the background model here, or
! we would load a background pattern from file.  In this version, we are
! using the background that was computed by the MC program, and has
! an energy histogram embedded in it, so we need to interpolate this
! histogram to the pixels of the scintillator.  In other words, we need
! to initialize a new accum_e array for the detector by interpolating
! from the Lambert projection of the MC results.
!
  nsx = (mcnl%numsx - 1)/2
  nsy = nsx
! determine the scale factor for the Lambert interpolation; the square has
! an edge length of 2 x sqrt(pi/2)
  scl = float(nsx) !  / LPs%sPio2  [removed on 09/01/15 by MDG for new Lambert routines]

! get the indices of the minimum and maximum energy
  Emin = nint((enl%energymin - mcnl%Ehistmin)/mcnl%Ebinsize) +1
  if (Emin.lt.1)  Emin=1
  if (Emin.gt.EBSDMCdata%numEbins)  Emin=EBSDMCdata%numEbins

  Emax = nint((enl%energymax - mcnl%Ehistmin)/mcnl%Ebinsize) +1
  if (Emax.lt.1)  Emax=1
  if (Emax.gt.EBSDMCdata%numEbins)  Emax=EBSDMCdata%numEbins

! correction of change in effective pixel area compared to equal-area Lambert projection
  alpha = atan(enl%delta/enl%L/sqrt(sngl(cPi)))
  ipx = enl%numsx/2 + nint(enl%xpc)
  ipy = enl%numsy/2 + nint(enl%ypc)
  if (TKD.eqv..TRUE.) then
    if ((abs(ipy).gt.enl%numsy).or.(abs(ipx).gt.enl%numsx)) then
      pcvec = (/enl%ypc*enl%delta*ca + enl%xpc*enl%delta*sa*sw + enl%L*cw*sa, &
               enl%L*sw - enl%xpc*enl%delta*cw,&
               enl%L*ca*cw + enl%xpc*enl%delta*ca*sw - enl%ypc*enl%delta*sa/)
      pcvec = pcvec/vecnorm(pcvec)
    else
      pcvec = (/ EBSDdetector%rgx(ipx,ipy), EBSDdetector%rgy(ipx,ipy), EBSDdetector%rgz(ipx,ipy) /)
    end if
  else
    pcvec = (/ EBSDdetector%rgx(ipx,ipy), EBSDdetector%rgy(ipx,ipy), EBSDdetector%rgz(ipx,ipy) /)
  end if
  calpha = cos(alpha)
  do i=1,enl%numsx
    do j=1,enl%numsy
! do the coordinate transformation for this detector pixel
       dc = (/ EBSDdetector%rgx(i,j),EBSDdetector%rgy(i,j),EBSDdetector%rgz(i,j) /)
! make sure the third one is positive; if not, switch all
       if (dc(3).lt.0.0) dc = -dc

! convert these direction cosines to coordinates in the Rosca-Lambert projection
        call L%setxyz( dc )
        istat = L%LambertSphereToSquare( ixy )
        ixy = ixy * scl
        x = ixy(1)
        ixy(1) = ixy(2)
        ixy(2) = -x
! four-point interpolation (bi-quadratic)
        nix = int(nsx+ixy(1))-nsx
        niy = int(nsy+ixy(2))-nsy
        dx = ixy(1)-nix
        dy = ixy(2)-niy
        dxm = 1.0-dx
        dym = 1.0-dy
! do the area correction for this detector pixel
        dp = dot_product(pcvec,dc)
        theta = acos(dp)
        if ((i.eq.ipx).and.(j.eq.ipy)) then
          g = 0.25
        else
          g = ((calpha*calpha + dp*dp - 1.0)**1.5)/(calpha**3) * 0.25
        end if
! interpolate the intensity
        do k=Emin,Emax
          s = EBSDMCdata%accum_e(k,nix,niy) * dxm * dym + &
              EBSDMCdata%accum_e(k,nix+1,niy) * dx * dym + &
              EBSDMCdata%accum_e(k,nix,niy+1) * dxm * dy + &
              EBSDMCdata%accum_e(k,nix+1,niy+1) * dx * dy
! EBSD intensities do not need to be flipped vertically, but TKD intensities apparently
! do need to be flipped... we need to look into this a bit more to make sure it is correct.
          if (TKD.eqv..TRUE.) then
            EBSDdetector%accum_e_detector(k,i,elp-j) = g * s
          else
            EBSDdetector%accum_e_detector(k,i,j) = g * s
          end if
        end do
    end do
  end do

if (present(verbose)) call Message%printMessage(' -> completed detector generation', frm = "(A)")

end associate

end subroutine GenerateDetector_

!--------------------------------------------------------------------------
subroutine ComputeEBSDPatterns_(self, EMsoft, MCFT, MPFT, HDF, HDFnames, mpnl, numangles, angles, progname, nmldeffile)
!DEC$ ATTRIBUTES DLLEXPORT :: ComputeEBSDPatterns_
!! author: MDG
!! version: 1.0
!! date: 02/17/20
!!
!! compute an energy-weighted EBSD pattern

use mod_EMsoft
use mod_symmetry
use mod_crystallography
use mod_io
use mod_diffraction
use mod_Lambert
use mod_quaternions
use mod_rotations
use HDF5
use mod_HDFsupport
use mod_HDFnames
use ISO_C_BINDING
use omp_lib
use mod_OMPsupport
use mod_timing
use stringconstants
use mod_math
use mod_filters
use mod_MCfiles
use mod_MPfiles

IMPLICIT NONE

class(EBSD_T), INTENT(INOUT)            :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
type(MCfile_T), INTENT(INOUT)           :: MCFT
type(MPfile_T), INTENT(INOUT)           :: MPFT
type(HDF_T),INTENT(INOUT)               :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames
type(EBSDmasterNameListType),INTENT(INOUT) :: mpnl
integer(kind=irg),INTENT(IN)            :: numangles
type(QuaternionArray_T), INTENT(IN)     :: angles
character(fnlen),INTENT(IN)             :: progname
character(fnlen),INTENT(IN)             :: nmldeffile

type(SpaceGroup_T)                      :: SG
type(IO_T)                              :: Message
type(q_T)                               :: qq, qq1, qq2, qq3
type(o_T)                               :: om
type(e_T)                               :: eu
type(Quaternion_T)                      :: quat
type(Timing_T)                          :: timer
type(Cell_T)                            :: cell

! all geometrical parameters and filenames
real(kind=dbl)                          :: prefactor, qz(3)

! allocatable arrays
real(kind=sgl),allocatable              :: EBSDpattern(:,:), binned(:,:)        ! array with EBSD patterns
real(kind=sgl),allocatable              :: z(:,:)               ! used to store the computed patterns before writing to disk
real(kind=sgl),allocatable              :: energywf(:), eulerangles(:,:)

! arrays for each OpenMP thread
real(kind=sgl),allocatable              :: tmLPNH(:,:,:) , tmLPSH(:,:,:)
real(kind=sgl),allocatable              :: trgx(:,:), trgy(:,:), trgz(:,:)          ! auxiliary detector arrays needed for interpolation
real(kind=sgl),allocatable              :: taccum(:,:,:)

! various items
integer(kind=irg)                       :: i, j, iang, jang, k, io_int(6), hdferr, L, correctsize, dim1, dim2          ! various counters
integer(kind=irg)                       :: istat, ipar(7), tick, tock, tickstart
integer(kind=irg)                       :: nix, niy, binx, biny, nixp, niyp, maxthreads,nextra,ninlastbatch,nlastremainder, npy     ! various parameters
integer(kind=irg)                       :: NUMTHREADS, TID   ! number of allocated threads, thread ID
integer(kind=irg)                       :: ninbatch, nbatches,nremainder,ibatch,nthreads,maskradius,nlastbatches, totnumbatches
integer(kind=irg),allocatable           :: istart(:,:), istop(:,:), patinbatch(:)

real(kind=sgl)                          :: bindx, ma, mi, tstart, tstop, io_real(3)
real(kind=dbl),parameter                :: nAmpere = 6.241D+18   ! Coulomb per second
integer(kind=irg),parameter             :: storemax = 20        ! number of EBSD patterns stored in one output block
integer(kind=irg)                       :: Emin, Emax      ! various parameters
real(kind=dbl)                          :: dc(3), scl, nel, emult           ! direction cosine array
real(kind=dbl)                          :: sx, dx, dxm, dy, dym, rhos, x         ! various parameters
real(kind=dbl)                          :: ixy(2), tmp

real(kind=sgl),allocatable              :: mask(:,:), lx(:), ly(:), masklin(:), binnedvec(:)
character(kind=c_char),allocatable      :: batchpatterns(:,:,:), bpat(:,:), threadbatchpatterns(:,:,:)
integer(kind=irg),allocatable           :: batchpatternsint(:,:,:), bpatint(:,:), threadbatchpatternsint(:,:,:)
real(kind=sgl),allocatable              :: batchpatterns32(:,:,:), threadbatchpatterns32(:,:,:), threadbatchpatterns32lin(:,:)
real(kind=sgl),allocatable              :: batchpatterns32lin(:,:)
integer(kind=irg),allocatable           :: acc_array(:,:)
real(kind=sgl),allocatable              :: master_arrayNH(:,:), master_arraySH(:,:), wf(:)
character(len=3)                        :: outputformat
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)

! parameter for random number generator
integer, parameter                      :: K4B=selected_int_kind(9)      ! used by ran function in math.f90
integer(K4B)                            :: idum

integer(HSIZE_T), dimension(1:3)        :: hdims, offset
integer(HSIZE_T), dimension(1:2)        :: hdims2, offset2
integer(HSIZE_T)                        :: dims2(2), dims3(3)
character(fnlen,kind=c_char)            :: line2(1)
character(fnlen)                        :: groupname, dataset, datagroupname, attributename, HDF_FileVersion
character(11)                           :: dstr
character(15)                           :: tstrb
character(15)                           :: tstre
character(10)                           :: char10
character(fnlen)                        :: datafile
logical                                 :: overwrite = .TRUE., insert = .TRUE., singlebatch
character(5)                            :: bitmode
integer(kind=irg)                       :: numbits
real(kind=sgl)                          :: bitrange

! new stuff: deformation tensor
real(kind=dbl)                          :: Umatrix(3,3), Fmatrix(3,3), Smatrix(3,3), quF(4), Fmatrix_inverse(3,3), &
                                           Gmatrix(3,3)
logical                                 :: includeFmatrix=.FALSE., noise, isTKD=.FALSE.

if (trim(HDFnames%get_ProgramData()).eq.trim(SC_TKD)) isTKD = .TRUE.

associate( enl => self%nml, mcnl => MCFT%nml, &
           EBSDMCdata => MCFT%MCDT, EBSDMPdata => MPFT%MPDT, EBSDdetector => self%det )

!====================================
! max number of OpenMP threads on this platform
maxthreads = omp_get_max_threads()

!====================================
! what is the output format?  GUI or BIN ?
outputformat = enl%outputformat

!====================================
! bit depth and format of output
call get_bit_parameters(enl%bitdepth, numbits, bitrange, bitmode)

if (enl%makedictionary.eq.'y') then
  bitmode = 'dict'
  call Message%printMessage('Program will work in dictionary generation mode')
end if

! define some energy-related parameters derived from MC input parameters
!====================================
noise = .FALSE.
if (enl%poisson.eq.'y') noise = .TRUE.

! make sure the requested energy range is within the range available from the Monte Carlo computation
if (enl%energymin.lt.mcnl%Ehistmin) enl%energymin = mcnl%Ehistmin
if (enl%energymax.gt.mcnl%EkeV) enl%energymax = mcnl%EkeV

! get the indices of the minimum and maximum energy
Emin = nint((enl%energymin - mcnl%Ehistmin)/mcnl%Ebinsize) + 1
if (Emin.lt.1)  Emin=1
if (Emin.gt.EBSDMCdata%numEbins)  Emin=EBSDMCdata%numEbins

Emax = nint((enl%energymax - mcnl%Ehistmin)/mcnl%Ebinsize) + 1
if (Emax.lt.1)  Emax=1
if (Emax.gt.EBSDMCdata%numEbins)  Emax=EBSDMCdata%numEbins

! modified by MDG, 03/26/18
!nel = sum(acc%accum_e_detector)
nel = float(mcnl%totnum_el) * float(mcnl%multiplier)
emult = nAmpere * 1e-9 / nel  ! multiplicative factor to convert MC data to an equivalent incident beam of 1 nanoCoulomb
io_real(1) = emult
call Message%WriteValue(' Multiplicative factor to generate 1 nC of incident electrons ', io_real, 1)
! intensity prefactor  (redefined by MDG, 3/23/18)
! prefactor = 0.25D0 * nAmpere * enl%beamcurrent * enl%dwelltime * 1.0D-15/ nel
prefactor = emult * enl%beamcurrent * enl%dwelltime * 1.0D-6
io_real(1) = prefactor
call Message%WriteValue(' Intensity scaling prefactor = ', io_real, 1)

allocate(energywf(Emin:Emax), wf(EBSDMCdata%numEbins),stat=istat)
energywf = 0.0
wf = 0.0

wf = sum(sum(EBSDdetector%accum_e_detector,3),2)
energywf(Emin:Emax) = wf(Emin:Emax)
energywf = energywf/sum(energywf)
deallocate(wf)

!====================================
! init a bunch of parameters
!====================================
! binned pattern array
  binx = enl%numsx/enl%binning
  biny = enl%numsy/enl%binning
  bindx = 1.0/float(enl%binning)**2
!====================================

! get the crystal structure data
call cell%getCrystalData(mcnl%xtalname, SG, EMsoft, useHDF=HDF)

!====================================
! ------ and open the output file (only thread 0 can write to this file)
!====================================
! we need to write the image dimensions, and also how many of those there are...

timer = Timing_T()
tstrb = timer%getTimeString()
dstr = timer%getDateString()

! Create a new file using the default properties.
datafile = EMsoft%generateFilePath('EMdatapathname', enl%datafile)

hdferr =  HDF%createFile(datafile)
if (hdferr.ne.0) call HDF%error_check('HDF_createFile ', hdferr)

!====================================
! new in Release 4.3: add a Manufacturer string (null terminated)
dataset = SC_Manufacturer
line2(1) = 'EMsoft'
line2(1) = cstringify(line2(1))
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
!====================================

! write the EMheader to the file
datagroupname = trim(HDFnames%get_ProgramData()) ! 'EBSD' or 'TKD'
call HDF%writeEMheader(EMsoft,dstr, tstrb, tstre, progname, datagroupname)

! add the CrystalData group at the top level of the file
call cell%addXtalDataGroup(SG, EMsoft, HDF)

! create a namelist group to write all the namelist files into
hdferr = HDF%createGroup(HDFnames%get_NMLfiles())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup NMLfiles', hdferr)

! read the text file and write the array to the file
if (isTKD.eqv..TRUE.) then
  dataset = SC_EMTKDNML
else
  dataset = SC_EMEBSDNML
end if
hdferr = HDF%writeDatasetTextFile(dataset, nmldeffile)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetTextFile EMEBSDNML/EMTKDNML', hdferr)

call HDF%pop()

! create a NMLparameters group to write all the namelist entries into
hdferr = HDF%createGroup(HDFnames%get_NMLparameters())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup NMLparameters', hdferr)

if (isTKD.eqv..TRUE.) then
  call self%writeHDFNameList_(HDF, HDFnames, isTKD)
else
  call self%writeHDFNameList_(HDF, HDFnames)
end if

! and leave this group
call HDF%pop()

! then the remainder of the data in a EMData group
hdferr = HDF%createGroup(HDFnames%get_EMData())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup EMData', hdferr)

! create the EBSD group and add a HDF_FileVersion attribute to it
hdferr = HDF%createGroup(datagroupname)
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup EBSD/TKD', hdferr)
! before Feb. 19, 2019, an undetected error caused all patterns to be upside down in the Kikuchi bands only,
! not in the background intensity profile.  This was compensated by a pattern flip of all experimental
! patterns in the dictionary indexing program, but when taking individual patterns from this program, they
! are actually upside down in all versions through HDF_FileVersion 4.0.  As of 4.1, the patterns are in the
! correct orientation.  This was detected by manually indexing a simulated pattern.
HDF_FileVersion = '4.1'
attributename = SC_HDFFileVersion
hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)

! =====================================================
! The following write commands constitute HDF_FileVersion = 4.0 and above
! =====================================================
dataset = SC_xtalname
allocate(stringarray(1))
stringarray(1)= trim(mcnl%xtalname)
hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetStringArray xtalname', hdferr)

dataset = SC_numangles
hdferr = HDF%writeDatasetInteger(dataset, numangles)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetInteger numangles', hdferr)

! and add the Euler angles to the output file
allocate(eulerangles(3,numangles))
do i=1,numangles
  quat = angles%getQuatfromArray(i)
  call qq%q_set( quat%get_quats() )
  eu = qq%qe()
  eulerangles(1:3,i) = eu%e_copy()
end do
dataset = SC_Eulerangles
hdferr = HDF%writeDatasetFloatArray(dataset, eulerangles/sngl(dtor), 3, numangles)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetFloatArray2D Eulerangles', hdferr)

! =====================================================
! end of HDF_FileVersion = 4.0 and above write statements
! =====================================================

! and we leave this group open for further data output from the main program loop ...

!====================================
! generate the deformation matrix and its polar decomposition
!====================================
includeFmatrix = .FALSE.
if (enl%applyDeformation.eq.'y') then
  includeFmatrix = .TRUE.
! importantly, we need to transpose the input deformation tensor for the
! computations to be performed correctly; that's because the namelist file
! has the tensor defined row by row instead of column-major
  Fmatrix = transpose(enl%Ftensor)

! the following lines are commented out for now; they are simply informational
! and not really essential

! ! perform the polar decomposition on the deformation tensor
!   call getPolarDecomposition(Fmatrix, Umatrix, Smatrix)
!   call Message%printMessage('')

! ! and convert the unitary matrix to a quaternion
!   call om%o_set(Umatrix)
!   quF = conjg(om%oq())

!   call Message%printMessage('Polar Decomposition')
!   call Message%printMessage('  --> Unitary Matrix')
!   call PrintMatrixd('U = ',Umatrix)
!   char10 = ''
!   call print_orientation(init_orientation(quF,'qu'),'qu',char10)
!   write(*,*) 'rotation angle = ',2.0*acos(quF(1))*180.D0/cPi
!   call Message%printMessage('  --> Stretch Matrix')
!   call PrintMatrixd('S = ',Smatrix)

! ! compute the effective lattice parameters for the given deformation, based on the
! ! undeformed unit cell
!   Gmatrix = matmul(matmul(transpose(Fmatrix),cell%dsm), transpose(matmul(transpose(Fmatrix),cell%dsm)) )
!   call Message%printMessage('Metric tensor for distorted cell:')
!   call PrintMatrixd('Gdis=',Gmatrix)
!   io_real(1:3) = (/ sqrt(Gmatrix(1,1)), sqrt(Gmatrix(2,2)), sqrt(Gmatrix(3,3)) /)
!   call WriteValue('(a, b, c) = ',io_real,3)

!   io_real(1:3) = (/ acos(Gmatrix(2,3)/sqrt(Gmatrix(2,2))/sqrt(Gmatrix(3,3)))*180.D0/cPi  , &
!                     acos(Gmatrix(1,3)/sqrt(Gmatrix(1,1))/sqrt(Gmatrix(3,3)))*180.D0/cPi  , &
!                     acos(Gmatrix(1,2)/sqrt(Gmatrix(1,1))/sqrt(Gmatrix(2,2)))*180.D0/cPi  /)
!   call WriteValue('(alpha, beta, gamma) = ',io_real,3)
!   call Message%printMessage('')

! invert the deformation tensor and pass it on to the pattern computation routine
  call mInvert(Fmatrix, Fmatrix_inverse, .FALSE.)
end if

!====================================
! ------ start the actual image computation loop
!====================================

!====================================
! to speed things up, we'll split the computation into batches of 1,024 patterns per thread; once those
! are computed, we leave the OpenMP part to write them to a file
!====================================


! and allocate space to store each batch; this requires some careful analysis
! since we are doing things in multiple threads
  nthreads = enl%nthreads
  if (nthreads.gt.maxthreads) then
    io_int(1) = maxthreads
    call Message%WriteValue('# threads requested is larger than available number; resetting to ',io_int, 1)
    nthreads = maxthreads
  end if
  ninbatch = 1024
  nlastbatches = 0
  singlebatch = .FALSE.
  if (numangles.ge.ninbatch*nthreads) then
    nbatches = numangles/(ninbatch*nthreads)
    nremainder = mod(numangles,ninbatch*nthreads)
    nextra = 0
    if (nremainder.gt.0) then
      singlebatch = .TRUE.
      nextra = 1
      ninlastbatch = nremainder/nthreads+1
      nlastremainder = nremainder - (nthreads-1)*ninlastbatch
    end if
  else
! if there are fewer patterns than ninbatch*nthreads we need to redefine ninbatch
    singlebatch = .TRUE.
    if (numangles.le.nthreads) then
      nthreads = 1
    end if
    nbatches = 0
    ninlastbatch = numangles/nthreads+1
    nlastremainder = numangles - (nthreads-1)*ninlastbatch
    nlastbatches = 1
    nextra = 0
    if (nlastremainder.gt.0) nextra = 1
end if
  if (nbatches.ne.0) then
    io_int(1) = numangles
    io_int(2) = ninbatch
    io_int(3) = nthreads
    io_int(4) = nbatches
    io_int(5) = nremainder
    call Message%WriteValue('  OpenMP loop variables : ',io_int,5,"(I10,' = ',I4,' * ',I2,' * ',I4,' + ',I6)")
  end if
  if ((ninlastbatch.ne.0).and.(nextra.ne.0)) then
    io_int(1) = numangles - nbatches * nthreads * ninbatch
    io_int(2) = ninlastbatch
    io_int(3) = nthreads-1
    io_int(4) = 1
    io_int(5) = nlastremainder
    call Message%WriteValue('  Remainder loop variables : ',io_int,5,"(I10,' = ',I4,' * ',I2,' * ',I4,' + ',I6)")
  end if

! allocate the istart and istop arrays for all the separate runs
  totnumbatches = nbatches + nextra
  allocate(istart(0:nthreads-1,totnumbatches))
  allocate(istop(0:nthreads-1,totnumbatches))
  allocate(patinbatch(totnumbatches))
  do i=1,nbatches
    do j=0,nthreads-1
      istart(j,i) = 1 + ninbatch * ( j + nthreads*(i-1) )
      istop(j,i)  = ninbatch * ( j+1 + nthreads*(i-1) )
    end do
  end do
  if (nextra.eq.1) then
    i = nbatches+1
    do j=0,nthreads-1
      istart(j,i) = nthreads*ninbatch*nbatches + 1 + ninlastbatch * j
      if (j.ne.nthreads-1) then
         istop(j,i) = nthreads*ninbatch*nbatches + ninlastbatch * ( j + 1 )
      else
         istop(j,i) = nthreads*ninbatch*nbatches + ninlastbatch * ( nthreads - 1 ) + nlastremainder
      end if
    end do
  end if
  patinbatch = sum(istop-istart,1) + nthreads

! and allocate the batchpatterns array for hyperslab writing [modified 8/25/17 for different output formats]
L = binx*biny
! make sure that correctsize is a multiple of 16; if not, make it so
if (mod(L,16) .ne. 0) then
    correctsize = 16*ceiling(float(L)/16.0)
else
    correctsize = L
end if

if (trim(bitmode).eq.'char') then
  allocate(batchpatterns(binx,biny,ninbatch*nthreads),stat=istat)
end if
if (trim(bitmode).eq.'int') then
  allocate(batchpatternsint(binx,biny,ninbatch*nthreads),stat=istat)
end if
if (trim(bitmode).eq.'float') then
  allocate(batchpatterns32(binx,biny,ninbatch*nthreads),stat=istat)
end if
if (trim(bitmode).eq.'dict') then
  allocate(batchpatterns32lin(correctsize,ninbatch*nthreads),stat=istat)
end if

!====================================
! here we also create a mask if necessary
  allocate(mask(binx,biny),masklin(binx*biny),stat=istat)
  mask = 1.0
  masklin = 1.0
  if (enl%maskpattern.eq.'y') then
! create the circular mask in a potentially rectangular array
    maskradius = (minval( (/ binx, biny /) ) / 2 )**2
    allocate(lx(binx), ly(biny), stat=istat)
    lx = (/ (float(i),i=1,binx) /) - float(binx/2)
    ly = (/ (float(i),i=1,biny) /) - float(biny/2)
    do i=1,binx
      do j=1,biny
        if ((lx(i)**2+ly(j)**2).gt.maskradius) mask(i,j) = 0.0
      end do
    end do
    deallocate(lx, ly)
    if (trim(bitmode).eq.'dict') then
      do j = 1,biny
        do i = 1,binx
          masklin((j-1)*binx+i) = mask(i,j)
        end do
      end do
    end if
  end if

!====================================
! determine the scale factor for the Lambert interpolation
scl = dble(mpnl%npx)

!====================================
! define the integer parameter list for the CalcEBSDPatternSingleFull call
ipar(1) = enl%binning
ipar(2) = enl%numsx
ipar(3) = enl%numsy
ipar(4) = mpnl%npx
ipar(5) = mpnl%npx
ipar(6) = EBSDMCdata%numEbins
ipar(7) = EBSDMCdata%numEbins

!====================================
! set the number of OpenMP threads
call OMP_setNThreads(nthreads)

call timer%Time_tick()

!====================================
!====================================
do ibatch=1,totnumbatches

! use OpenMP to run on multiple cores ...
!$OMP PARALLEL default(shared)  PRIVATE(TID,iang,i,j,istat,EBSDpattern,binned,idum,bpat,ma,mi,threadbatchpatterns,bpatint)&
!$OMP& PRIVATE(tmLPNH, tmLPSH, trgx, trgy, trgz, taccum, dims2, dims3, threadbatchpatternsint, threadbatchpatterns32)&
!$OMP& PRIVATE(binnedvec, threadbatchpatterns32lin)

  NUMTHREADS = OMP_GET_NUM_THREADS()
  TID = OMP_GET_THREAD_NUM()

! initialize the random number generator for the Poison noise
  if (noise.eqv..TRUE.) then
    idum = -1-TID
  else
    idum = 0_K4B
  end if

! each thread needs a private copy of the master and accum arrays; not having
! those can produce poor scaling...
  dims2 = shape(self%det%rgx)
  allocate(trgx(dims2(1),dims2(2)), trgy(dims2(1),dims2(2)), trgz(dims2(1),dims2(2)))
  dims3 = shape(self%det%accum_e_detector)
  allocate(taccum(dims3(1),dims3(2),dims3(3)))
  dims3 = shape(MPFT%MPDT%mLPNH)
  allocate(tmLPNH(dims3(1),dims3(2),dims3(3)), tmLPSH(dims3(1),dims3(2),dims3(3)))
! and copy the data in
  trgx = self%det%rgx
  trgy = self%det%rgy
  trgz = self%det%rgz
  taccum = self%det%accum_e_detector
  tmLPNH = MPFT%MPDT%mLPNH
  tmLPSH = MPFT%MPDT%mLPSH

! allocate the arrays that will hold the computed pattern
  allocate(binned(binx,biny),stat=istat)
  if (trim(bitmode).eq.'char') then
    allocate(bpat(binx,biny),stat=istat)
  end if
  if (trim(bitmode).eq.'int') then
    allocate(bpatint(binx,biny),stat=istat)
  end if
  if (trim(bitmode).eq.'dict') then
    allocate(bpatint(binx,biny),stat=istat)
    allocate(binnedvec(correctsize),stat=istat)
  end if

! this array requires some care in terms of its size parameters...
  if ((singlebatch.eqv..TRUE.).AND.(ibatch.eq.totnumbatches)) then
     if (TID.eq.nthreads-1) then
      if (trim(bitmode).eq.'char') then
        allocate(threadbatchpatterns(binx,biny,nlastremainder),stat=istat)
      end if
      if (trim(bitmode).eq.'int') then
        allocate(threadbatchpatternsint(binx,biny,nlastremainder),stat=istat)
      end if
      if (trim(bitmode).eq.'float') then
        allocate(threadbatchpatterns32(binx,biny,nlastremainder),stat=istat)
      end if
      if (trim(bitmode).eq.'dict') then
        allocate(threadbatchpatterns32lin(correctsize,nlastremainder),stat=istat)
      end if
    else
      if (trim(bitmode).eq.'char') then
        allocate(threadbatchpatterns(binx,biny,ninlastbatch),stat=istat)
      end if
      if (trim(bitmode).eq.'int') then
        allocate(threadbatchpatternsint(binx,biny,ninlastbatch),stat=istat)
      end if
      if (trim(bitmode).eq.'float') then
        allocate(threadbatchpatterns32(binx,biny,ninlastbatch),stat=istat)
      end if
      if (trim(bitmode).eq.'dict') then
        allocate(threadbatchpatterns32lin(correctsize,ninlastbatch),stat=istat)
      end if
    end if
  else
    if (trim(bitmode).eq.'char') then
      allocate(threadbatchpatterns(binx,biny,ninbatch),stat=istat)
    end if
    if (trim(bitmode).eq.'int') then
      allocate(threadbatchpatternsint(binx,biny,ninbatch),stat=istat)
    end if
    if (trim(bitmode).eq.'float') then
      allocate(threadbatchpatterns32(binx,biny,ninbatch),stat=istat)
    end if
    if (trim(bitmode).eq.'dict') then
      allocate(threadbatchpatterns32lin(correctsize,ninbatch),stat=istat)
    end if
  end if

  if (trim(bitmode).eq.'char') then
    threadbatchpatterns = ' '
  end if
  if (trim(bitmode).eq.'int') then
    threadbatchpatternsint = 0_irg
  end if
  if (trim(bitmode).eq.'float') then
    threadbatchpatterns32 = 0.0
  end if
  if (trim(bitmode).eq.'dict') then
    threadbatchpatterns32lin = 0.0
  end if

  do iang=istart(TID,ibatch),istop(TID,ibatch)
! convert the direction cosines to quaternions, include the
! sample quaternion orientation, and then back to direction cosines...
! then convert these individually to the correct EBSD pattern location
    binned = 0.0

    if (includeFmatrix.eqv..TRUE.) then
     if (self%nml%includebackground.eq.'y') then
      call self%CalcEBSDPatternSingleFull(ipar,angles%getQuatfromArray(iang),taccum,tmLPNH,tmLPSH,trgx,trgy,trgz,binned, &
                                          Emin,Emax,mask,prefactor,Fmatrix_inverse,applynoise=idum)
     else
      call self%CalcEBSDPatternSingleFull(ipar,angles%getQuatfromArray(iang),taccum,tmLPNH,tmLPSH,trgx,trgy,trgz,binned, &
                                          Emin,Emax,mask,prefactor,Fmatrix_inverse,removebackground='y',applynoise=idum)
     end if
    else
     if (self%nml%includebackground.eq.'y') then
      call self%CalcEBSDPatternSingleFull(ipar,angles%getQuatfromArray(iang),taccum,tmLPNH,tmLPSH,trgx,trgy,trgz,binned, &
                                          Emin,Emax,mask,prefactor,applynoise=idum)
     else
      call self%CalcEBSDPatternSingleFull(ipar,angles%getQuatfromArray(iang),taccum,tmLPNH,tmLPSH,trgx,trgy,trgz,binned, &
                                          Emin,Emax,mask,prefactor,removebackground='y',applynoise=idum)
     end if
    end if

    if (self%nml%scalingmode .eq. 'gam') then
        binned = binned**self%nml%gammavalue
    end if

    if (trim(bitmode).eq.'dict') then  ! pre-process the patterns for dictionary indexing
! this step includes adaptive histogram equalization, masking, and normalization

! adaptive histogram equalization
        ma = maxval(binned)
        mi = minval(binned)
        bpatint = nint(((binned - mi)/ (ma-mi))*255.0)
        binned =  float(adhisteq(self%nml%nregions,binx,biny,bpatint))

! linearize the array and apply the mask
        binnedvec = 0.0
        do j= 1,biny
          do i = 1,binx
            binnedvec((j-1)*binx+i) = binned(i,j)
          end do
        end do

! apply circular mask and normalize
        binnedvec(1:L) = binnedvec(1:L) * masklin(1:L)
        binnedvec(1:correctsize) = binnedvec(1:correctsize)/vecnorm(binnedvec(1:correctsize))

! store in array for hyperslab writing
        threadbatchpatterns32lin(1:correctsize, iang-istart(TID,ibatch)+1) = binnedvec

    else  ! don't make a dictionary so no preprocessing of patterns... just intensity scaling to the correct range
      if (trim(bitmode).eq.'char') then
        ma = maxval(binned)
        mi = minval(binned)
        binned = mask * ((binned - mi)/ (ma-mi))
        bpat = char(nint(bitrange*binned))

        threadbatchpatterns(1:binx,1:biny, iang-istart(TID,ibatch)+1) = bpat
      end if

      if (trim(bitmode).eq.'int') then
        ma = maxval(binned)
        mi = minval(binned)
        binned = mask * ((binned - mi)/ (ma-mi))
        bpatint = nint(bitrange*binned)

        threadbatchpatternsint(1:binx,1:biny, iang-istart(TID,ibatch)+1) = bpatint
      end if

      if (trim(bitmode).eq.'float') then
        threadbatchpatterns32(1:binx,1:biny, iang-istart(TID,ibatch)+1) = binned
      end if
    end if

  end do ! end of iang loop

! and now we write the threadbatchpatterns arrays to the main batch patterns array; we
! need to intercept the special case when there are remainder patterns!
!$OMP CRITICAL
  if ((singlebatch.eqv..TRUE.).AND.(ibatch.eq.totnumbatches)) then
    if (TID.eq.nthreads-1) then
      if (trim(bitmode).eq.'char') then
        batchpatterns(1:binx,1:biny,TID*ninlastbatch+1:TID*ninlastbatch+nlastremainder)=&
          threadbatchpatterns(1:binx,1:biny, 1:nlastremainder)
      end if

      if (trim(bitmode).eq.'int') then
        batchpatternsint(1:binx,1:biny,TID*ninlastbatch+1:TID*ninlastbatch+nlastremainder)=&
          threadbatchpatternsint(1:binx,1:biny, 1:nlastremainder)
      end if

      if (trim(bitmode).eq.'float') then
        batchpatterns32(1:binx,1:biny,TID*ninlastbatch+1:TID*ninlastbatch+nlastremainder)=&
          threadbatchpatterns32(1:binx,1:biny, 1:nlastremainder)
      end if

      if (trim(bitmode).eq.'dict') then
        batchpatterns32lin(1:correctsize,TID*ninlastbatch+1:TID*ninlastbatch+nlastremainder)=&
          threadbatchpatterns32lin(1:correctsize, 1:nlastremainder)
      end if

    else
      if (trim(bitmode).eq.'char') then
        batchpatterns(1:binx,1:biny, TID*ninlastbatch+1:(TID+1)*ninlastbatch) = &
          threadbatchpatterns(1:binx,1:biny, 1:ninlastbatch)
      end if

      if (trim(bitmode).eq.'int') then
        batchpatternsint(1:binx,1:biny, TID*ninlastbatch+1:(TID+1)*ninlastbatch) = &
          threadbatchpatternsint(1:binx,1:biny, 1:ninlastbatch)
      end if

      if (trim(bitmode).eq.'float') then
        batchpatterns32(1:binx,1:biny, TID*ninlastbatch+1:(TID+1)*ninlastbatch) = &
          threadbatchpatterns32(1:binx,1:biny, 1:ninlastbatch)
      end if

      if (trim(bitmode).eq.'dict') then
        batchpatterns32lin(1:correctsize,TID*ninlastbatch+1:(TID+1)*ninlastbatch)=&
          threadbatchpatterns32lin(1:correctsize, 1:ninlastbatch)
      end if
    end if
  else
    if (trim(bitmode).eq.'char') then
      batchpatterns(1:binx,1:biny, TID*ninbatch+1:(TID+1)*ninbatch) = threadbatchpatterns(1:binx,1:biny, 1:ninbatch)
    end if

    if (trim(bitmode).eq.'int') then
      batchpatternsint(1:binx,1:biny, TID*ninbatch+1:(TID+1)*ninbatch) = threadbatchpatternsint(1:binx,1:biny, 1:ninbatch)
    end if

    if (trim(bitmode).eq.'float') then
      batchpatterns32(1:binx,1:biny, TID*ninbatch+1:(TID+1)*ninbatch) = threadbatchpatterns32(1:binx,1:biny, 1:ninbatch)
    end if

    if (trim(bitmode).eq.'dict') then
      batchpatterns32lin(1:correctsize, TID*ninbatch+1:(TID+1)*ninbatch) = &
         threadbatchpatterns32lin(1:correctsize, 1:ninbatch)
    end if

  end if
!$OMP END CRITICAL

!$OMP END PARALLEL

! here we write all the entries in the batchpatterns array to the HDF file as a hyperslab
! =====================================================
! The following write commands constitute HDF_FileVersion = 4.0
! =====================================================
if (isTKD.eqv..TRUE.) then
  dataset = SC_TKDpatterns
else
  dataset = SC_EBSDpatterns
end if
 !if (outputformat.eq.'bin') then
   if (trim(bitmode).eq.'dict') then
     offset2 = (/ 0, (ibatch-1)*ninbatch*enl%nthreads /)
     hdims2 = (/ correctsize, numangles /)
     dims2 = (/ correctsize, patinbatch(ibatch) /)
     dim1 = patinbatch(ibatch)
   else
     offset = (/ 0, 0, (ibatch-1)*ninbatch*enl%nthreads /)
     hdims = (/ binx, biny, numangles /)
     dims3 = (/ binx, biny, patinbatch(ibatch) /)
     dim2 = patinbatch(ibatch)
   end if
   if (ibatch.eq.1) then
     if (trim(bitmode).eq.'char') then
       hdferr = HDF%writeHyperslabCharArray(dataset, batchpatterns(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabCharArray EBSDpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'int') then
       hdferr = HDF%writeHyperslabIntegerArray(dataset, batchpatternsint(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabIntegerArray EBSDpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'float') then
       hdferr = HDF%writeHyperslabFloatArray(dataset, batchpatterns32(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabFloatArray EBSDpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'dict') then
       hdferr = HDF%writeHyperslabFloatArray(dataset, batchpatterns32lin(1:correctsize,1:dim1), hdims2, offset2, &
                                              dims2)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabFloatArray EBSDpatterns', hdferr)
     end if
   else
     if (trim(bitmode).eq.'char') then
       hdferr = HDF%writeHyperslabCharArray(dataset, batchpatterns(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3, insert)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabCharArray EBSDpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'int') then
       hdferr = HDF%writeHyperslabIntegerArray(dataset, batchpatternsint(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3, insert)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabIntegerArray EBSDpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'float') then
       hdferr = HDF%writeHyperslabFloatArray(dataset, batchpatterns32(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3, insert)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabFloatArray EBSDpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'dict') then
       hdferr = HDF%writeHyperslabFloatArray(dataset, batchpatterns32lin(1:correctsize,1:dim1), hdims2, offset2, &
                                              dims2, insert)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabFloatArray EBSDpatterns', hdferr)
     end if
   end if
 !end if
! =====================================================
! end of HDF_FileVersion = 4.0 write statements
! =====================================================


 io_int(1) = ibatch
 io_int(2) = totnumbatches
 call Message%WriteValue('Completed cycle ',io_int,2,"(I4,' of ',I4)")

end do
!====================================
!====================================

call timer%Time_tock()
tstop = timer%getInterval()

io_int(1) = tstop
call Message%WriteValue('Execution time [system_clock()] = ',io_int,1,"(I8,' [s]')")

call HDF%pop()
call HDF%pop()

! and update the end time
timer = timing_T()
tstre = timer%getTimeString()

hdferr = HDF%openGroup(HDFnames%get_EMheader())
if (hdferr.ne.0) call HDF%error_check('HDF_openGroup EMheader', hdferr)

hdferr = HDF%openGroup(HDFnames%get_ProgramData())
if (hdferr.ne.0) call HDF%error_check('HDF_openGroup EBSD/TKD', hdferr)

! stop time /EMheader/StopTime 'character'
dataset = SC_StopTime
line2(1) = dstr//', '//tstre
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetStringArray StopTime', hdferr)

dataset = SC_Duration
hdferr = HDF%writeDatasetFloat(dataset, tstop)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetFloat Duration', hdferr)

! close the datafile
call HDF%pop(.TRUE.)


end associate

end subroutine ComputeEBSDPatterns_

!--------------------------------------------------------------------------
recursive subroutine CalcEBSDPatternSingleFull_(self,ipar,qq,accum,mLPNH,mLPSH,rgx,rgy,rgz,binned,Emin,Emax,mask, &
                                                prefactor, Fmatrix, removebackground, applynoise)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcEBSDPatternSingleFull_
 !! author: MDG
 !! version: 1.0
 !! date: 01/06/20
 !!
 !! compute a single EBSD pattern, used in many programs

use mod_Lambert
use mod_quaternions
use mod_rotations

IMPLICIT NONE

integer, parameter                              :: K4B=selected_int_kind(9)

class(EBSD_T), INTENT(INOUT)                    :: self
integer(kind=irg),INTENT(IN)                    :: ipar(7)
type(Quaternion_T),INTENT(IN)                   :: qq
real(kind=dbl),INTENT(IN)                       :: prefactor
integer(kind=irg),INTENT(IN)                    :: Emin, Emax
real(kind=sgl),INTENT(IN)                       :: accum(ipar(6),ipar(2),ipar(3))
real(kind=sgl),INTENT(IN)                       :: mLPNH(-ipar(4):ipar(4),-ipar(5):ipar(5),ipar(7))
real(kind=sgl),INTENT(IN)                       :: mLPSH(-ipar(4):ipar(4),-ipar(5):ipar(5),ipar(7))
real(kind=sgl),INTENT(IN)                       :: rgx(ipar(2),ipar(3))
real(kind=sgl),INTENT(IN)                       :: rgy(ipar(2),ipar(3))
real(kind=sgl),INTENT(IN)                       :: rgz(ipar(2),ipar(3))
real(kind=sgl),INTENT(OUT)                      :: binned(ipar(2)/ipar(1),ipar(3)/ipar(1))
real(kind=sgl),INTENT(IN)                       :: mask(ipar(2)/ipar(1),ipar(3)/ipar(1))
real(kind=dbl),INTENT(IN),optional              :: Fmatrix(3,3)
character(1),INTENT(IN),OPTIONAL                :: removebackground
integer(K4B),INTENT(INOUT),OPTIONAL             :: applynoise

real(kind=sgl),allocatable                      :: EBSDpattern(:,:)
real(kind=sgl),allocatable                      :: wf(:)
real(kind=sgl)                                  :: dc(3),ixy(2),scl,bindx, tmp
real(kind=sgl)                                  :: dx,dy,dxm,dym, x, y, z
integer(kind=irg)                               :: ii,jj,kk,istat
integer(kind=irg)                               :: nix,niy,nixp,niyp
logical                                         :: nobg, noise


! ipar(1) = ebsdnl%binning
! ipar(2) = ebsdnl%numsx
! ipar(3) = ebsdnl%numsy
! ipar(4) = ebsdnl%npx
! ipar(5) = ebsdnl%npy
! ipar(6) = ebsdnl%numEbins
! ipar(7) = ebsdnl%nE

nobg = .FALSE.
if (present(removebackground)) then
  if (removebackground.eq.'y') nobg = .TRUE.
end if

noise = .FALSE.
if (present(applynoise)) then
  if (applynoise.ne.0_K4B) noise = .TRUE.
end if

allocate(EBSDpattern(ipar(2),ipar(3)),stat=istat)

binned = 0.0
EBSDpattern = 0.0

scl = float(ipar(4))

do ii = 1,ipar(2)
    do jj = 1,ipar(3)
! get the pixel direction cosines from the pre-computed array
        dc = (/ rgx(ii,jj),rgy(ii,jj),rgz(ii,jj) /)
! apply the grain rotation
        dc = sngl( qq%quat_Lp( dble(dc) ) )
! apply the deformation if present
        if (present(Fmatrix)) then
          dc = matmul(sngl(Fmatrix), dc)
        end if
! and normalize the direction cosines (to remove any rounding errors)
        dc = dc/sqrt(sum(dc**2))

! convert these direction cosines to interpolation coordinates in the Rosca-Lambert projection
        call LambertgetInterpolation(dc, scl, ipar(4), ipar(5), nix, niy, nixp, niyp, dx, dy, dxm, dym)

! interpolate the intensity
        if (nobg.eqv..TRUE.) then
          if (dc(3) .ge. 0.0) then
            do kk = Emin, Emax
                EBSDpattern(ii,jj) = EBSDpattern(ii,jj) + ( mLPNH(nix,niy,kk) * dxm * dym + &
                                               mLPNH(nixp,niy,kk) * dx * dym + mLPNH(nix,niyp,kk) * dxm * dy + &
                                               mLPNH(nixp,niyp,kk) * dx * dy )

            end do
          else
            do kk = Emin, Emax
                EBSDpattern(ii,jj) = EBSDpattern(ii,jj) + ( mLPSH(nix,niy,kk) * dxm * dym + &
                                               mLPSH(nixp,niy,kk) * dx * dym + mLPSH(nix,niyp,kk) * dxm * dy + &
                                               mLPSH(nixp,niyp,kk) * dx * dy )

            end do

          end if
        else
          if (dc(3) .ge. 0.0) then
            do kk = Emin, Emax
                EBSDpattern(ii,jj) = EBSDpattern(ii,jj) + accum(kk,ii,jj) * ( mLPNH(nix,niy,kk) * dxm * dym + &
                                               mLPNH(nixp,niy,kk) * dx * dym + mLPNH(nix,niyp,kk) * dxm * dy + &
                                               mLPNH(nixp,niyp,kk) * dx * dy )

            end do
          else
            do kk = Emin, Emax
                EBSDpattern(ii,jj) = EBSDpattern(ii,jj) + accum(kk,ii,jj) * ( mLPSH(nix,niy,kk) * dxm * dym + &
                                               mLPSH(nixp,niy,kk) * dx * dym + mLPSH(nix,niyp,kk) * dxm * dy + &
                                               mLPSH(nixp,niyp,kk) * dx * dy )

            end do
          end if
        end if
    end do
end do

EBSDpattern = prefactor * EBSDpattern

! do we need to apply Poisson noise ?  (slow...)
! if (noise.eqv..TRUE.) then
!   EBSDpattern = applyPoissonNoise( EBSDpattern, ipar(2), ipar(3), applynoise )
! end if

! do we need to bin the pattern ?

! 17/07/2020 Clément Lafond, temporary fix to avoid NaN value, and negatives
! values when EBSD pattern size is large
do ii=1,ipar(2)1
    do jj=1,ipar(3)
        if (isnan(EBSDpattern(ii,jj)).or.EBSDpattern(ii,jj).lt.0.0) then
          EBSDpattern(ii,jj) = 0.0
        end if
    end do
end do

if (ipar(1) .ne. 1) then
    do ii=1,ipar(2),ipar(1)
        do jj=1,ipar(3),ipar(1)
            binned(ii/ipar(1)+1,jj/ipar(1)+1) = &
            sum(EBSDpattern(ii:ii+ipar(1)-1,jj:jj+ipar(1)-1))
        end do
    end do
! and divide by binning^2
!   binned = binned * bindx
else
    binned = EBSDpattern
end if

binned = binned * mask

end subroutine CalcEBSDPatternSingleFull_

!--------------------------------------------------------------------------
subroutine ComputedeformedEBSDPatterns_(self, EMsoft, MCFT, MPFT, HDF, HDFnames, &
                                        mpnl, numangles, angles, orpcdef, progname, nmldeffile)
!DEC$ ATTRIBUTES DLLEXPORT :: ComputedeformedEBSDPatterns_
  !! author: MDG
  !! version: 1.0
  !! date: 02/18/20
  !!
  !! compute a energy-weighted EBSD patterns with different pattern centers and deformation states

use mod_EMsoft
use mod_symmetry
use mod_crystallography
use mod_io
use mod_diffraction
use mod_Lambert
use mod_quaternions
use mod_rotations
use HDF5
use mod_HDFsupport
use mod_HDFnames
use ISO_C_BINDING
use omp_lib
use mod_OMPsupport
use mod_timing
use stringconstants
use mod_math
use mod_filters
use mod_MCfiles
use mod_MPfiles

IMPLICIT NONE

class(EBSD_T), INTENT(INOUT)            :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
type(MCfile_T), INTENT(INOUT)           :: MCFT
type(MPfile_T), INTENT(INOUT)           :: MPFT
type(HDF_T),INTENT(INOUT)               :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames
type(EBSDmasterNameListType),INTENT(INOUT) :: mpnl
integer(kind=irg),INTENT(IN)            :: numangles
type(QuaternionArray_T), INTENT(IN)     :: angles
type(EBSDAnglePCDefType),INTENT(IN)     :: orpcdef
character(fnlen),INTENT(IN)             :: progname
character(fnlen),INTENT(IN)             :: nmldeffile

type(SpaceGroup_T)                      :: SG
type(IO_T)                              :: Message
type(q_T)                               :: qq, qq1, qq2, qq3
type(o_T)                               :: om
type(e_T)                               :: eu
type(Quaternion_T)                      :: quat
type(Timing_T)                          :: timer
type(Cell_T)                            :: cell

! all geometrical parameters and filenames
real(kind=dbl)                          :: prefactor, qz(3)

! allocatable arrays
real(kind=sgl),allocatable              :: EBSDpattern(:,:), binned(:,:)        ! array with EBSD patterns
real(kind=sgl),allocatable              :: z(:,:)               ! used to store the computed patterns before writing to disk
real(kind=sgl),allocatable              :: energywf(:), eulerangles(:,:)

! arrays for each OpenMP thread
real(kind=sgl),allocatable              :: tmLPNH(:,:,:) , tmLPSH(:,:,:)
real(kind=sgl),allocatable              :: trgx(:,:), trgy(:,:), trgz(:,:)          ! auxiliary detector arrays needed for interpolation
real(kind=sgl),allocatable              :: taccum(:,:,:)

! various items
integer(kind=irg)                       :: i, j, iang, jang, k, io_int(6), hdferr, dim2          ! various counters
integer(kind=irg)                       :: istat, ipar(7), tick, tock, tickstart
integer(kind=irg)                       :: nix, niy, binx, biny, nixp, niyp, maxthreads,nextra,ninlastbatch,nlastremainder     ! various parameters
integer(kind=irg)                       :: NUMTHREADS, TID   ! number of allocated threads, thread ID
integer(kind=irg)                       :: ninbatch, nbatches,nremainder,ibatch,nthreads,maskradius,nlastbatches, totnumbatches
integer(kind=irg),allocatable           :: istart(:,:), istop(:,:), patinbatch(:)

real(kind=sgl)                          :: bindx, ma, mi, tstart, tstop, io_real(3)
real(kind=dbl),parameter                :: nAmpere = 6.241D+18   ! Coulomb per second
integer(kind=irg),parameter             :: storemax = 20        ! number of EBSD patterns stored in one output block
integer(kind=irg)                       :: Emin, Emax      ! various parameters
real(kind=dbl)                          :: dc(3), scl, nel, emult           ! direction cosine array
real(kind=dbl)                          :: sx, dx, dxm, dy, dym, rhos, x         ! various parameters
real(kind=dbl)                          :: ixy(2), tmp

real(kind=sgl),allocatable              :: mask(:,:), lx(:), ly(:)
character(kind=c_char),allocatable      :: batchpatterns(:,:,:), bpat(:,:), threadbatchpatterns(:,:,:)
integer(kind=irg),allocatable           :: batchpatternsint(:,:,:), bpatint(:,:), threadbatchpatternsint(:,:,:)
real(kind=sgl),allocatable              :: batchpatterns32(:,:,:), threadbatchpatterns32(:,:,:)
integer(kind=irg),allocatable           :: acc_array(:,:)
real(kind=sgl),allocatable              :: master_arrayNH(:,:), master_arraySH(:,:), wf(:)
character(len=3)                        :: outputformat
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)

! parameter for random number generator
integer, parameter                      :: K4B=selected_int_kind(9)      ! used by ran function in math.f90
integer(K4B)                            :: idum

integer(HSIZE_T), dimension(1:3)        :: hdims, offset
integer(HSIZE_T)                        :: dims2(2), dims3(3)
character(fnlen,kind=c_char)            :: line2(1)
character(fnlen)                        :: groupname, dataset, datagroupname, attributename, HDF_FileVersion
character(11)                           :: dstr
character(15)                           :: tstrb
character(15)                           :: tstre
character(10)                           :: char10
character(fnlen)                        :: datafile
logical                                 :: overwrite = .TRUE., insert = .TRUE., singlebatch
character(5)                            :: bitmode
integer(kind=irg)                       :: numbits
real(kind=sgl)                          :: bitrange

! new stuff: deformation tensor
real(kind=dbl)                          :: Umatrix(3,3), Fmatrix(3,3), Smatrix(3,3), quF(4), Fmatrix_inverse(3,3), &
                                           Gmatrix(3,3)
logical                                 :: includeFmatrix=.FALSE., isTKD=.FALSE.

if (trim(HDFnames%get_ProgramData()).eq.trim(SC_TKD)) isTKD = .TRUE.

associate( enl => self%nml, mcnl => MCFT%nml, &
           EBSDMCdata => MCFT%MCDT, EBSDMPdata => MPFT%MPDT, EBSDdetector => self%det )

!====================================
! max number of OpenMP threads on this platform
maxthreads = omp_get_max_threads()

!====================================
! what is the output format?  GUI or BIN ?
outputformat = enl%outputformat

!====================================
! bit depth and format of output
call get_bit_parameters(enl%bitdepth, numbits, bitrange, bitmode)

! define some energy-related parameters derived from MC input parameters
!====================================
! make sure the requested energy range is within the range available from the Monte Carlo computation
if (enl%energymin.lt.mcnl%Ehistmin) enl%energymin = mcnl%Ehistmin
if (enl%energymax.gt.mcnl%EkeV) enl%energymax = mcnl%EkeV

! get the indices of the minimum and maximum energy
Emin = nint((enl%energymin - mcnl%Ehistmin)/mcnl%Ebinsize) +1
if (Emin.lt.1)  Emin=1
if (Emin.gt.EBSDMCdata%numEbins)  Emin=EBSDMCdata%numEbins

Emax = nint((enl%energymax - mcnl%Ehistmin)/mcnl%Ebinsize) +1
if (Emax.lt.1)  Emax=1
if (Emax.gt.EBSDMCdata%numEbins)  Emax=EBSDMCdata%numEbins

! modified by MDG, 03/26/18
nel = float(mcnl%totnum_el) * float(mcnl%multiplier)
emult = nAmpere * 1e-9 / nel  ! multiplicative factor to convert MC data to an equivalent incident beam of 1 nanoCoulomb
io_real(1) = emult
call Message%WriteValue(' Multiplicative factor to generate 1 nC of incident electrons ', io_real, 1)

!====================================
! init a bunch of parameters
!====================================
! binned pattern array
  binx = enl%numsx/enl%binning
  biny = enl%numsy/enl%binning
  bindx = 1.0/float(enl%binning)**2

!====================================

! get the crystal structure data
call cell%getCrystalData(mcnl%xtalname, SG, EMsoft, useHDF=HDF)

!====================================
! ------ and open the output file for IDL visualization (only thread 0 can write to this file)
!====================================
! we need to write the image dimensions, and also how many of those there are...

timer = Timing_T()
tstrb = timer%getTimeString()
dstr = timer%getDateString()

! Create a new file using the default properties.
datafile = EMsoft%generateFilePath('EMdatapathname', enl%datafile)

hdferr =  HDF%createFile(datafile)
if (hdferr.ne.0) call HDF%error_check('HDF_createFile ', hdferr)

! write the EMheader to the file
datagroupname = trim(HDFnames%get_ProgramData()) ! 'EBSD' or 'TKD'
call HDF%writeEMheader(EMsoft,dstr, tstrb, tstre, progname, datagroupname)

! add the CrystalData group at the top level of the file
call cell%addXtalDataGroup(SG, EMsoft, HDF)

! create a namelist group to write all the namelist files into
hdferr = HDF%createGroup(HDFnames%get_NMLfiles())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup NMLfiles', hdferr)

! read the text file and write the array to the file
if (isTKD.eqv..TRUE.) then
  dataset = SC_EMTKDNML
else
  dataset = SC_EMEBSDNML
end if
hdferr = HDF%writeDatasetTextFile(dataset, nmldeffile)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetTextFile EMEBSDNML/EMTKDNML', hdferr)

call HDF%pop()

! create a NMLparameters group to write all the namelist entries into
hdferr = HDF%createGroup(HDFnames%get_NMLparameters())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup NMLparameters', hdferr)

if (isTKD.eqv..TRUE.) then
  call self%writeHDFNameList_(HDF, HDFnames, isTKD)
else
  call self%writeHDFNameList_(HDF, HDFnames)
end if

! and leave this group
call HDF%pop()

! then the remainder of the data in a EMData group
hdferr = HDF%createGroup(HDFnames%get_EMData())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup EMData', hdferr)

! create the EBSD group and add a HDF_FileVersion attribute to it
hdferr = HDF%createGroup(datagroupname)
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup EBSD/TKD', hdferr)
! before Feb. 19, 2019, an undetected error caused all patterns to be upside down in the Kikuchi bands only,
! not in the background intensity profile.  This was compensated by a pattern flip of all experimental
! patterns in the dictionary indexing program, but when taking individual patterns from this program, they
! are actually upside down in all versions through HDF_FileVersion 4.0.  As of 4.1, the patterns are in the
! correct orientation.  This was detected by manually indexing a simulated pattern.
HDF_FileVersion = '4.1'
attributename = SC_HDFFileVersion
hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)

dataset = SC_xtalname
allocate(stringarray(1))
stringarray(1)= trim(mcnl%xtalname)
hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetStringArray xtalname', hdferr)

dataset = SC_numangles
hdferr = HDF%writeDatasetInteger(dataset, numangles)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetInteger numangles', hdferr)

! and add the Euler angles to the output file
allocate(eulerangles(3,numangles))
do i=1,numangles
  quat = angles%getQuatfromArray(i)
  call qq%q_set( quat%get_quats() )
  eu = qq%qe()
  eulerangles(1:3,i) = eu%e_copy()
end do
dataset = SC_Eulerangles
hdferr = HDF%writeDatasetFloatArray(dataset, eulerangles/sngl(dtor), 3, numangles)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetFloatArray2D Eulerangles', hdferr)

! and we leave this group open for further data output from the main program loop ...

!====================================
! in this particular routine we always include the deformation tensors, no matter what the nml file says
!====================================
includeFmatrix = .TRUE.
!====================================

!====================================
! ------ start the actual image computation loop
!====================================

!====================================
! to speed things up, we'll split the computation into batches of 1,024 patterns per thread; once those
! are computed, we leave the OpenMP part to write them to a file
!====================================


! and allocate space to store each batch; this requires some careful analysis
! since we are doing things in multiple threads
  nthreads = enl%nthreads
  if (nthreads.gt.maxthreads) then
    io_int(1) = maxthreads
    call Message%WriteValue('# threads requested is larger than available number; resetting to ',io_int, 1)
    nthreads = maxthreads
  end if
  ninbatch = 1024
  nlastbatches = 0
  singlebatch = .FALSE.
  if (numangles.ge.ninbatch*nthreads) then
    nbatches = numangles/(ninbatch*nthreads)
    nremainder = mod(numangles,ninbatch*nthreads)
    nextra = 0
    if (nremainder.gt.0) then
      singlebatch = .TRUE.
      nextra = 1
      ninlastbatch = nremainder/nthreads+1
      nlastremainder = nremainder - (nthreads-1)*ninlastbatch
    end if
  else
! if there are fewer patterns than ninbatch*nthreads we need to redefine ninbatch
    singlebatch = .TRUE.
    if (numangles.le.nthreads) then
      nthreads = 1
    end if
    nbatches = 0
    ninlastbatch = numangles/nthreads+1
    nlastremainder = numangles - (nthreads-1)*ninlastbatch
    nlastbatches = 1
    nextra = 0
    if (nlastremainder.gt.0) nextra = 1
end if
  if (nbatches.ne.0) then
    io_int(1) = numangles
    io_int(2) = ninbatch
    io_int(3) = nthreads
    io_int(4) = nbatches
    io_int(5) = nremainder
    call Message%WriteValue('  OpenMP loop variables : ',io_int,5,"(I10,' = ',I4,' * ',I2,' * ',I4,' + ',I6)")
  end if
  if ((ninlastbatch.ne.0).and.(nextra.ne.0)) then
    io_int(1) = numangles - nbatches * nthreads * ninbatch
    io_int(2) = ninlastbatch
    io_int(3) = nthreads-1
    io_int(4) = 1
    io_int(5) = nlastremainder
    call Message%WriteValue('  Remainder loop variables : ',io_int,5,"(I10,' = ',I4,' * ',I2,' * ',I4,' + ',I6)")
  end if

! allocate the istart and istop arrays for all the separate runs
  totnumbatches = nbatches + nextra
  allocate(istart(0:nthreads-1,totnumbatches))
  allocate(istop(0:nthreads-1,totnumbatches))
  allocate(patinbatch(totnumbatches))
  do i=1,nbatches
    do j=0,nthreads-1
      istart(j,i) = 1 + ninbatch * ( j + nthreads*(i-1) )
      istop(j,i)  = ninbatch * ( j+1 + nthreads*(i-1) )
    end do
  end do
  if (nextra.eq.1) then
    i = nbatches+1
    do j=0,nthreads-1
      istart(j,i) = nthreads*ninbatch*nbatches + 1 + ninlastbatch * j
      if (j.ne.nthreads-1) then
         istop(j,i) = nthreads*ninbatch*nbatches + ninlastbatch * ( j + 1 )
      else
         istop(j,i) = nthreads*ninbatch*nbatches + ninlastbatch * ( nthreads - 1 ) + nlastremainder
      end if
    end do
  end if
  patinbatch = sum(istop-istart,1) + nthreads

! and allocate the batchpatterns array for hyperslab writing [modified 8/25/17 for different output formats]
if (trim(bitmode).eq.'char') then
  allocate(batchpatterns(binx,biny,ninbatch*nthreads),stat=istat)
end if
if (trim(bitmode).eq.'int') then
  allocate(batchpatternsint(binx,biny,ninbatch*nthreads),stat=istat)
end if
if (trim(bitmode).eq.'float') then
  allocate(batchpatterns32(binx,biny,ninbatch*nthreads),stat=istat)
end if

!====================================
! here we also create a mask if necessary
  allocate(mask(binx,biny),stat=istat)
  mask = 1.0
  if (enl%maskpattern.eq.'y') then
! create the circular mask in a potentially rectangular array
    maskradius = (minval( (/ binx, biny /) ) / 2 )**2
    allocate(lx(binx), ly(biny), stat=istat)
    lx = (/ (float(i),i=1,binx) /) - float(binx/2)
    ly = (/ (float(i),i=1,biny) /) - float(biny/2)
    do i=1,binx
      do j=1,biny
        if ((lx(i)**2+ly(j)**2).gt.maskradius) mask(i,j) = 0.0
      end do
    end do
    deallocate(lx, ly)
  end if

!====================================
! determine the scale factor for the Lambert interpolation
scl = dble(mpnl%npx)

!====================================
! define the integer parameter list for the CalcEBSDPatternSingleFull call
ipar(1) = enl%binning
ipar(2) = enl%numsx
ipar(3) = enl%numsy
ipar(4) = mpnl%npx
ipar(5) = mpnl%npx
ipar(6) = EBSDMCdata%numEbins
ipar(7) = EBSDMCdata%numEbins

!====================================
! set the number of OpenMP threads
call OMP_setNThreads(nthreads)

call timer%Time_tick()

!====================================
!====================================
do ibatch=1,totnumbatches

! use OpenMP to run on multiple cores ...
!$OMP PARALLEL default(shared)  PRIVATE(TID,iang,i,j,istat,EBSDpattern,binned,idum,bpat,ma,mi,threadbatchpatterns,bpatint)&
!$OMP& PRIVATE(tmLPNH, tmLPSH, trgx, trgy, trgz, taccum, threadbatchpatternsint, threadbatchpatterns32, prefactor)&
!$OMP& PRIVATE(Fmatrix_inverse, nel, Fmatrix)

  NUMTHREADS = OMP_GET_NUM_THREADS()
  TID = OMP_GET_THREAD_NUM()

! each thread needs a private copy of the master and accum arrays; not having
! those can produce poor scaling... in addition, they need to be recomputed for each pattern !
  allocate(trgx(enl%numsx,enl%numsy), trgy(enl%numsx,enl%numsy), trgz(enl%numsx,enl%numsy))
  allocate(taccum(EBSDMCdata%numEbins,enl%numsx,enl%numsy))
  allocate(tmLPNH(enl%numsx,enl%numsy,EBSDMCdata%numEbins), tmLPSH(enl%numsx,enl%numsy,EBSDMCdata%numEbins))
! and copy the data in
  tmLPNH = EBSDMPdata%mLPNH
  tmLPSH = EBSDMPdata%mLPSH

! allocate the arrays that will hold the computed pattern
  allocate(binned(binx,biny),stat=istat)
  if (trim(bitmode).eq.'char') then
    allocate(bpat(binx,biny),stat=istat)
  end if
  if (trim(bitmode).eq.'int') then
    allocate(bpatint(binx,biny),stat=istat)
  end if

! this array requires some care in terms of its size parameters...
  if ((singlebatch.eqv..TRUE.).AND.(ibatch.eq.totnumbatches)) then
     if (TID.eq.nthreads-1) then
      if (trim(bitmode).eq.'char') then
        allocate(threadbatchpatterns(binx,biny,nlastremainder),stat=istat)
      end if
      if (trim(bitmode).eq.'int') then
        allocate(threadbatchpatternsint(binx,biny,nlastremainder),stat=istat)
      end if
      if (trim(bitmode).eq.'float') then
        allocate(threadbatchpatterns32(binx,biny,nlastremainder),stat=istat)
      end if
    else
      if (trim(bitmode).eq.'char') then
        allocate(threadbatchpatterns(binx,biny,ninlastbatch),stat=istat)
      end if
      if (trim(bitmode).eq.'int') then
        allocate(threadbatchpatternsint(binx,biny,ninlastbatch),stat=istat)
      end if
      if (trim(bitmode).eq.'float') then
        allocate(threadbatchpatterns32(binx,biny,ninlastbatch),stat=istat)
      end if
    end if
  else
    if (trim(bitmode).eq.'char') then
      allocate(threadbatchpatterns(binx,biny,ninbatch),stat=istat)
    end if
    if (trim(bitmode).eq.'16bit') then
      allocate(threadbatchpatternsint(binx,biny,ninbatch),stat=istat)
    end if
    if (trim(bitmode).eq.'float') then
      allocate(threadbatchpatterns32(binx,biny,ninbatch),stat=istat)
    end if
  end if

  if (trim(enl%bitdepth).eq.'char') then
    threadbatchpatterns = ' '
  end if
  if (trim(enl%bitdepth).eq.'int') then
    threadbatchpatternsint = 0_irg
  end if
  if (trim(enl%bitdepth).eq.'float') then
    threadbatchpatterns32 = 0.0
  end if

  do iang=istart(TID,ibatch),istop(TID,ibatch)
! invert the transposed deformation tensor for this pattern
    Fmatrix = transpose(orpcdef%deftensors(1:3,1:3,iang))
    call mInvert(Fmatrix, Fmatrix_inverse, .FALSE.)

! for each pattern we need to compute the detector arrays
    if (enl%includebackground.eq.'y') then
      call self%GeneratemyEBSDDetector(MCFT, enl%numsx, enl%numsy, EBSDMCdata%numEbins, trgx, trgy, trgz, taccum, &
                                  orpcdef%pcs(1:3,iang),bg=.TRUE.)
! intensity prefactor
      prefactor = emult * enl%beamcurrent * enl%dwelltime * 1.0D-6
    else
      call self%GeneratemyEBSDDetector(MCFT, enl%numsx, enl%numsy, EBSDMCdata%numEbins, trgx, trgy, trgz, taccum, &
                                  orpcdef%pcs(1:3,iang),bg=.FALSE.)
! we pick a reasonable value here ...
      prefactor = 3.D0 * enl%beamcurrent * enl%dwelltime * 1.0D-6
    end if

    binned = 0.0

!   write (*,*) TID, nel, maxval(trgx), maxval(taccum), maxval(Fmatrix_inverse)

    if (includeFmatrix.eqv..TRUE.) then
     if (enl%includebackground.eq.'y') then
      call self%CalcEBSDPatternSingleFull(ipar,angles%getQuatfromArray(iang),taccum,tmLPNH,tmLPSH,trgx,trgy,trgz,binned, &
                                     Emin,Emax,mask,prefactor,Fmatrix_inverse)
     else
      call self%CalcEBSDPatternSingleFull(ipar,angles%getQuatfromArray(iang),taccum,tmLPNH,tmLPSH,trgx,trgy,trgz,binned, &
                                     Emin,Emax,mask,prefactor,Fmatrix_inverse,removebackground='y')
     end if
    else
     if (enl%includebackground.eq.'y') then
      call self%CalcEBSDPatternSingleFull(ipar,angles%getQuatfromArray(iang),taccum,tmLPNH,tmLPSH,trgx,trgy,trgz,binned, &
                                     Emin,Emax,mask,prefactor)
     else
      call self%CalcEBSDPatternSingleFull(ipar,angles%getQuatfromArray(iang),taccum,tmLPNH,tmLPSH,trgx,trgy,trgz,binned, &
                                     Emin,Emax,mask,prefactor,removebackground='y')
     end if
    end if

! from here on everything is the same as before...

    if (enl%scalingmode .eq. 'gam') then
        binned = binned**enl%gammavalue
    end if

    if (trim(bitmode).eq.'char') then
      ma = maxval(binned)
      mi = minval(binned)
      binned = mask * ((binned - mi)/ (ma-mi))
      bpat = char(nint(bitrange*binned))

      threadbatchpatterns(1:binx,1:biny, iang-istart(TID,ibatch)+1) = bpat
    end if

    if (trim(bitmode).eq.'int') then
      ma = maxval(binned)
      mi = minval(binned)
      binned = mask * ((binned - mi)/ (ma-mi))
      bpatint = nint(bitrange*binned)

      threadbatchpatternsint(1:binx,1:biny, iang-istart(TID,ibatch)+1) = bpatint
    end if

    if (trim(bitmode).eq.'float') then
      threadbatchpatterns32(1:binx,1:biny, iang-istart(TID,ibatch)+1) = binned
    end if

  end do ! end of iang loop

! and now we write the threadbatchpatterns arrays to the main batch patterns array; we
! need to intercept the special case when there are remainder patterns!
!$OMP CRITICAL
  if ((singlebatch.eqv..TRUE.).AND.(ibatch.eq.totnumbatches)) then
    if (TID.eq.nthreads-1) then
      if (trim(bitmode).eq.'char') then
        batchpatterns(1:binx,1:biny,TID*ninlastbatch+1:TID*ninlastbatch+nlastremainder)=&
          threadbatchpatterns(1:binx,1:biny, 1:nlastremainder)
      end if

      if (trim(bitmode).eq.'int') then
        batchpatternsint(1:binx,1:biny,TID*ninlastbatch+1:TID*ninlastbatch+nlastremainder)=&
          threadbatchpatternsint(1:binx,1:biny, 1:nlastremainder)
      end if

      if (trim(bitmode).eq.'float') then
        batchpatterns32(1:binx,1:biny,TID*ninlastbatch+1:TID*ninlastbatch+nlastremainder)=&
          threadbatchpatterns32(1:binx,1:biny, 1:nlastremainder)
      end if

    else
      if (trim(bitmode).eq.'char') then
        batchpatterns(1:binx,1:biny, TID*ninlastbatch+1:(TID+1)*ninlastbatch) = &
          threadbatchpatterns(1:binx,1:biny, 1:ninlastbatch)
      end if

      if (trim(bitmode).eq.'int') then
        batchpatternsint(1:binx,1:biny, TID*ninlastbatch+1:(TID+1)*ninlastbatch) = &
          threadbatchpatternsint(1:binx,1:biny, 1:ninlastbatch)
      end if

      if (trim(bitmode).eq.'float') then
        batchpatterns32(1:binx,1:biny, TID*ninlastbatch+1:(TID+1)*ninlastbatch) = &
          threadbatchpatterns32(1:binx,1:biny, 1:ninlastbatch)
      end if
    end if
  else
    if (trim(bitmode).eq.'char') then
      batchpatterns(1:binx,1:biny, TID*ninbatch+1:(TID+1)*ninbatch) = threadbatchpatterns(1:binx,1:biny, 1:ninbatch)
    end if

    if (trim(bitmode).eq.'int') then
      batchpatternsint(1:binx,1:biny, TID*ninbatch+1:(TID+1)*ninbatch) = threadbatchpatternsint(1:binx,1:biny, 1:ninbatch)
    end if

    if (trim(bitmode).eq.'float') then
      batchpatterns32(1:binx,1:biny, TID*ninbatch+1:(TID+1)*ninbatch) = threadbatchpatterns32(1:binx,1:biny, 1:ninbatch)
    end if
  end if
!$OMP END CRITICAL

!$OMP END PARALLEL

! here we write all the entries in the batchpatterns array to the HDF file as a hyperslab
if (isTKD.eqv..TRUE.) then
  dataset = SC_TKDpatterns
else
  dataset = SC_EBSDpatterns
end if
 !if (outputformat.eq.'bin') then
   offset = (/ 0, 0, (ibatch-1)*ninbatch*enl%nthreads /)
   hdims = (/ binx, biny, numangles /)
   dims3 = (/ binx, biny, patinbatch(ibatch) /)
   dim2 = patinbatch(ibatch)
   if (ibatch.eq.1) then
     if (trim(bitmode).eq.'char') then
       hdferr = HDF%writeHyperslabCharArray(dataset, batchpatterns(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabCharArray EBSDpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'int') then
       hdferr = HDF%writeHyperslabIntegerArray(dataset, batchpatternsint(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabIntegerArray EBSDpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'float') then
       hdferr = HDF%writeHyperslabFloatArray(dataset, batchpatterns32(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabFloatArray EBSDpatterns', hdferr)
     end if
   else
     if (trim(bitmode).eq.'char') then
       hdferr = HDF%writeHyperslabCharArray(dataset, batchpatterns(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3, insert)
       if (hdferr.ne.0) call  HDF%error_check('HDF_writeHyperslabCharArray EBSDpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'int') then
       hdferr = HDF%writeHyperslabIntegerArray(dataset, batchpatternsint(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3, insert)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabIntegerArray EBSDpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'float') then
       hdferr = HDF%writeHyperslabFloatArray(dataset, batchpatterns32(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3, insert)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabFloatArray EBSDpatterns', hdferr)
     end if
   end if
 !end if

 io_int(1) = ibatch
 io_int(2) = totnumbatches
 call Message%WriteValue('Completed cycle ',io_int,2,"(I4,' of ',I4)")

end do
!====================================
!====================================

call timer%Time_tock()
tstop = timer%getInterval()

io_int(1) = tstop
call Message%WriteValue('Execution time [system_clock()] = ',io_int,1,"(I8,' [s]')")

call HDF%pop()
call HDF%pop()

! and update the end time
timer = timing_T()
tstre = timer%getTimeString()

hdferr = HDF%openGroup(HDFnames%get_EMheader())
if (hdferr.ne.0) call HDF%error_check('HDF_openGroup EMheader', hdferr)

if (isTKD.eqv..TRUE.) then
  datagroupname = "TKD"
else
  datagroupname = "EBSD"
end if
hdferr = HDF%openGroup(datagroupname)
if (hdferr.ne.0) call HDF%error_check('HDF_openGroup EBSD/TKD', hdferr)

! stop time /EMheader/StopTime 'character'
dataset = SC_StopTime
line2(1) = dstr//', '//tstre
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetStringArray StopTime', hdferr)

dataset = SC_Duration
hdferr = HDF%writeDatasetFloat(dataset, tstop)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetFloat Duration', hdferr)

! close the datafile
call HDF%pop(.TRUE.)

end associate

end subroutine ComputedeformedEBSDPatterns_

!--------------------------------------------------------------------------
recursive subroutine GeneratemyEBSDDetector_(self, MCFT, nsx, nsy, numE, tgx, tgy, tgz, accum_e_detector, patcntr, bg)
!DEC$ ATTRIBUTES DLLEXPORT :: GeneratemyEBSDDetector_
!! author: MDG
!! version: 1.0
!! date: 02/05/20
!!
!! generate the detector arrays for the case where each pattern has a (slightly) different detector configuration

use mod_io
use mod_Lambert
use mod_MCfiles

IMPLICIT NONE

class(EBSD_T), INTENT(INOUT)            :: self
type(MCfile_T), INTENT(INOUT)           :: MCFT
integer(kind=irg),INTENT(IN)            :: nsx
integer(kind=irg),INTENT(IN)            :: nsy
integer(kind=irg),INTENT(IN)            :: numE
real(kind=sgl),INTENT(INOUT)            :: tgx(nsx,nsy)
!f2py intent(in,out) ::  tgx
real(kind=sgl),INTENT(INOUT)            :: tgy(nsx,nsy)
!f2py intent(in,out) ::  tgy
real(kind=sgl),INTENT(INOUT)            :: tgz(nsx,nsy)
!f2py intent(in,out) ::  tgz
real(kind=sgl),INTENT(INOUT)            :: accum_e_detector(numE,nsx,nsy)
!f2py intent(in,out) ::  accum_e_detector
real(kind=sgl),INTENT(IN)               :: patcntr(3)
logical,INTENT(IN),OPTIONAL             :: bg

type(Lambert_T)                         :: La

real(kind=sgl),allocatable              :: scin_x(:), scin_y(:), testarray(:,:)                 ! scintillator coordinate ararays [microns]
real(kind=sgl)                          :: alp, ca, sa, cw, sw
real(kind=sgl)                          :: L2, Ls, Lc, calpha     ! distances
real(kind=sgl),allocatable              :: z(:,:)
integer(kind=irg)                       :: nix, niy, binx, biny , i, j, Emin, Emax, istat, k, ipx, ipy, nx, ny, elp     ! various parameters
real(kind=sgl)                          :: dc(3), scl, alpha, theta, g, pcvec(3), s, dp           ! direction cosine array
real(kind=sgl)                          :: sx, dx, dxm, dy, dym, rhos, x, bindx, xpc, ypc, L         ! various parameters
real(kind=sgl)                          :: ixy(2)

associate( enl => self%nml, mcnl => MCFT%nml, EBSDMCdata => MCFT%MCDT )

!====================================
! ------ generate the detector arrays
!====================================
xpc = patcntr(1)
ypc = patcntr(2)
L = patcntr(3)

allocate(scin_x(nsx),scin_y(nsy),stat=istat)
! if (istat.ne.0) then ...
scin_x = - ( -xpc - ( 1.0 - nsx ) * 0.5 - (/ (i-1, i=1,nsx) /) ) * enl%delta
scin_y = ( ypc - ( 1.0 - nsy ) * 0.5 - (/ (i-1, i=1,nsy) /) ) * enl%delta

! auxiliary angle to rotate between reference frames
alp = 0.5 * cPi - (mcnl%sig - enl%thetac) * dtor
ca = cos(alp)
sa = sin(alp)

cw = cos(mcnl%omega * dtor)
sw = sin(mcnl%omega * dtor)

! we will need to incorporate a series of possible distortions
! here as well, as described in Gert nolze's paper; for now we
! just leave this place holder comment instead

! compute auxilliary interpolation arrays
! if (istat.ne.0) then ...

elp = nsy + 1
L2 = L * L
do j=1,nsx
  sx = L2 + scin_x(j) * scin_x(j)
  Ls = -sw * scin_x(j) + L*cw
  Lc = cw * scin_x(j) + L*sw
  do i=1,nsy
   rhos = 1.0/sqrt(sx + scin_y(i)**2)
   tgx(j,elp-i) = (scin_y(i) * ca + sa * Ls) * rhos!Ls * rhos
   tgy(j,elp-i) = Lc * rhos!(scin_x(i) * cw + Lc * sw) * rhos
   tgz(j,elp-i) = (-sa * scin_y(i) + ca * Ls) * rhos!(-sw * scin_x(i) + Lc * cw) * rhos
  end do
end do
deallocate(scin_x, scin_y)

! normalize the direction cosines.
allocate(z(enl%numsx,enl%numsy))
  z = 1.0/sqrt(tgx*tgx+tgy*tgy+tgz*tgz)
  tgx = tgx*z
  tgy = tgy*z
  tgz = tgz*z
deallocate(z)
!====================================

!====================================
! ------ create the equivalent detector energy array
!====================================
! from the Monte Carlo energy data, we need to extract the relevant
! entries for the detector geometry defined above.  Once that is
! done, we can get rid of the larger energy array
!
! in the old version, we either computed the background model here, or
! we would load a background pattern from file.  In this version, we are
! using the background that was computed by the MC program, and has
! an energy histogram embedded in it, so we need to interpolate this
! histogram to the pixels of the scintillator.  In other words, we need
! to initialize a new accum_e array for the detector by interpolating
! from the Lambert projection of the MC results.
!
nx = (mcnl%numsx - 1)/2
ny = nsx
if (present(bg)) then
 if (bg.eqv..TRUE.) then
! determine the scale factor for the Lambert interpolation; the square has
! an edge length of 2 x sqrt(pi/2)
  scl = float(nx) !  / LPs%sPio2  [removed on 09/01/15 by MDG for new Lambert routines]

! get the indices of the minimum and maximum energy
  Emin = nint((enl%energymin - mcnl%Ehistmin)/mcnl%Ebinsize) +1
  if (Emin.lt.1)  Emin=1
  if (Emin.gt.EBSDMCdata%numEbins)  Emin=EBSDMCdata%numEbins

  Emax = nint((enl%energymax - mcnl%Ehistmin)/mcnl%Ebinsize) +1
  if (Emax.lt.1)  Emax=1
  if (Emax.gt.EBSDMCdata%numEbins)  Emax=EBSDMCdata%numEbins

! correction of change in effective pixel area compared to equal-area Lambert projection
  alpha = atan(enl%delta/L/sqrt(sngl(cPi)))
  ipx = nsx/2 + nint(xpc)
  ipy = nsy/2 + nint(ypc)
  pcvec = (/ tgx(ipx,ipy), tgy(ipx,ipy), tgz(ipx,ipy) /)
  calpha = cos(alpha)
  do i=1,nsx
    do j=1,nsy
! do the coordinate transformation for this detector pixel
       dc = (/ tgx(i,j),tgy(i,j),tgz(i,j) /)
! make sure the third one is positive; if not, switch all
       if (dc(3).lt.0.0) dc = -dc
! convert these direction cosines to coordinates in the Rosca-Lambert projection
        call La%setxyz( dc )
        istat = La%LambertSphereToSquare( ixy )
        ixy = ixy * scl
        x = ixy(1)
        ixy(1) = ixy(2)
        ixy(2) = -x
! four-point interpolation (bi-quadratic)
        nix = int(nx+ixy(1))-nx
        niy = int(ny+ixy(2))-ny
        dx = ixy(1)-nix
        dy = ixy(2)-niy
        dxm = 1.0-dx
        dym = 1.0-dy
! do the area correction for this detector pixel
        dp = dot_product(pcvec,dc)
        theta = acos(dp)
        if ((i.eq.ipx).and.(j.eq.ipy)) then
          g = 0.25
        else
          g = ((calpha*calpha + dp*dp - 1.0)**1.5)/(calpha**3)
        end if
! interpolate the intensity
        do k=Emin,Emax
          s = EBSDMCdata%accum_e(k,nix,niy) * dxm * dym + &
              EBSDMCdata%accum_e(k,nix+1,niy) * dx * dym + &
              EBSDMCdata%accum_e(k,nix,niy+1) * dxm * dy + &
              EBSDMCdata%accum_e(k,nix+1,niy+1) * dx * dy
          accum_e_detector(k,i,j) = g * s
        end do
    end do
  end do
 else
   accum_e_detector = 1.0
 end if
end if

end associate

end subroutine GeneratemyEBSDDetector_

end module mod_EBSD
