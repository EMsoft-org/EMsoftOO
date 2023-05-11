! ###################################################################
! Copyright (c) 2023-2023, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_EBSDPCA
  !! author: MDG
  !! version: 1.0
  !! date: 05/08/23
  !!
  !! class definition for the EMEBSDPCA program

use mod_kinds
use mod_global
use mod_quaternions
use mod_MPfiles
use mod_EBSD

IMPLICIT NONE
private

! namelist for the EMEBSDPCA program
type, public :: EBSDPCANameListType
  integer(kind=irg)       :: numsx
  integer(kind=irg)       :: numsy
  integer(kind=irg)       :: binning
  integer(kind=irg)       :: nthreads
  ! integer(kind=irg)       :: energyaverage    [removed option, 05/08/23, MDG]
  integer(kind=irg)       :: maskradius
  integer(kind=irg)       :: nregions
  real(kind=sgl)          :: L
  real(kind=sgl)          :: thetac
  real(kind=sgl)          :: delta
  real(kind=sgl)          :: xpc
  real(kind=sgl)          :: ypc
  real(kind=sgl)          :: energymin
  real(kind=sgl)          :: energymax
  real(kind=sgl)          :: gammavalue
  real(kind=sgl)          :: axisangle(4)
  real(kind=sgl)          :: alphaBD
  real(kind=sgl)          :: hipassw
  real(kind=dbl)          :: beamcurrent
  real(kind=dbl)          :: dwelltime
  character(1)            :: poisson
  character(1)            :: includebackground
  character(1)            :: maskpattern
  character(3)            :: scalingmode
  character(3)            :: outputformat
  character(1)            :: spatialaverage
  character(5)            :: bitdepth
  character(fnlen)        :: anglefile
  character(fnlen)        :: masterfile
  character(fnlen)        :: energyfile
  character(fnlen)        :: datafile
end type EBSDPCANameListType

! class definition
type, public :: EBSDPCA_T
private
  character(fnlen)                    :: nmldeffile = 'EMEBSDPCA.nml'
  type(EBSDPCANameListType), public   :: nml
  type(EBSD_T), public                :: EBSD

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: EBSDPCA_
  procedure, pass(self) :: ComputeEBSDPCA_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: EBSDPCA => EBSDPCA_
  generic, public :: ComputeEBSDPCA => ComputeEBSDPCA_

end type EBSDPCA_T

! the constructor routine for this class
interface EBSDPCA_T
  module procedure EBSDPCA_constructor
end interface EBSDPCA_T

contains

!--------------------------------------------------------------------------
type(EBSDPCA_T) function EBSDPCA_constructor( nmlfile ) result(EBSDPCA)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDPCA_constructor
!! author: MDG
!! version: 1.0
!! date: 05/08/23
!!
!! constructor for the EBSDPCA_T Class; reads the name list

IMPLICIT NONE

character(fnlen)    :: nmlfile

call EBSDPCA%readNameList(nmlfile)

end function EBSDPCA_constructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 05/08/23
!!
!! read the namelist from an nml file for the EBSDPCA_T Class

use mod_io

IMPLICIT NONE

class(EBSDPCA_T), INTENT(INOUT)             :: self
character(fnlen),INTENT(IN)                 :: nmlfile
logical,OPTIONAL,INTENT(IN)                 :: initonly

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
real(kind=sgl)          :: energymin
real(kind=sgl)          :: energymax
real(kind=sgl)          :: gammavalue
real(kind=sgl)          :: alphaBD
real(kind=sgl)          :: axisangle(4)
real(kind=sgl)          :: hipassw
real(kind=dbl)          :: beamcurrent
real(kind=dbl)          :: dwelltime
character(1)            :: includebackground
character(1)            :: poisson
character(1)            :: maskpattern
character(1)            :: spatialaverage
character(3)            :: scalingmode
character(3)            :: eulerconvention
character(3)            :: outputformat
character(5)            :: bitdepth
character(fnlen)        :: anglefile
character(fnlen)        :: masterfile
character(fnlen)        :: energyfile  ! removed from template file 05/16/19 [MDG]
character(fnlen)        :: datafile

! define the IO namelist to facilitate passing variables to the program.
namelist  / EBSDPCAdata / L, thetac, delta, numsx, numsy, xpc, ypc, anglefile, eulerconvention, masterfile, bitdepth, &
                          energyfile, datafile, beamcurrent, dwelltime, energymin, energymax, binning, gammavalue, alphaBD, &
                          scalingmode, axisangle, nthreads, outputformat, maskpattern, spatialaverage, &
                          includebackground, hipassw, nregions, &
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
energymin       = 15.0          ! minimum energy to consider
energymax       = 30.0          ! maximum energy to consider
gammavalue      = 1.0           ! gamma factor
alphaBD         = 0.0           ! transfer lens barrel distortion parameter
maskradius      = 240           ! mask radius
hipassw         = 0.05          ! hi-pass filter radius
axisangle       = (/0.0, 0.0, 1.0, 0.0/)        ! no additional axis angle rotation
beamcurrent     = 14.513D0      ! beam current (actually emission current) in nano ampere
dwelltime       = 100.0D0       ! in microseconds
poisson         = 'n'           ! apply poisson noise ?
includebackground = 'y'         ! set to 'n' to remove realistic background intensity profile
maskpattern     = 'n'           ! 'y' or 'n' to include a circular mask
scalingmode     = 'not'         ! intensity selector ('lin', 'gam', or 'not')
eulerconvention = 'tsl'         ! convention for the first Euler angle ['tsl' or 'hkl']
outputformat    = 'gui'         ! output format for 'bin' or 'gui' use
bitdepth        = '8bit'        ! format for output; '8char' for [0..255], '##int' for integers, 'float' for floats
! the '##int' notation stands for the actual bitdepth; all values are stored as 32bit integers, but they are scaled
! from the float values to a maximum that is given by the first two digits, which indicate the bit depth; so, valid
! values would be '10int' for a 10-bit integer scale, '16int' for a 16-bit integer scale, and so on.
anglefile       = 'undefined'   ! filename
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
 read(UNIT=dataunit,NML=EBSDPCAdata)
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
self%nml%beamcurrent = beamcurrent
self%nml%dwelltime = dwelltime
self%nml%includebackground = includebackground
self%nml%poisson = poisson
self%nml%maskpattern = maskpattern
self%nml%scalingmode = scalingmode
self%nml%outputformat = outputformat
self%nml%bitdepth = bitdepth
self%nml%anglefile = anglefile
self%nml%masterfile = masterfile
! we require energyfile to be identical to masterfile, so the
! user definition, if any, in the namelist file is overwritten here...
self%nml%energyfile = masterfile       ! changed on 05/16/19 [MDG]
self%nml%datafile = datafile
self%nml%spatialaverage = spatialaverage

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 05/08/23
!!
!! pass the namelist for the EBSDPCA_T Class to the calling program

IMPLICIT NONE

class(EBSDPCA_T), INTENT(INOUT)          :: self
type(EBSDPCANameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames, isTKD)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG
!! version: 1.0
!! date: 05/08/23
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants

use ISO_C_BINDING

IMPLICIT NONE

class(EBSDPCA_T), INTENT(INOUT)            :: self
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

dataset = SC_poisson
line2(1) = trim(enl%poisson)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create poisson dataset', hdferr)

dataset = SC_includebackground
line2(1) = trim(enl%includebackground)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create includebackground dataset', hdferr)

dataset = SC_scalingmode
line2(1) = trim(enl%scalingmode)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create scalingmode dataset', hdferr)

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

dataset = SC_datafile
line2(1) = trim(enl%datafile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create datafile dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine EBSDPCA_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDPCA_
!! author: MDG
!! version: 1.0
!! date: 05/08/23
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
use mod_memory

IMPLICIT NONE

class(EBSDPCA_T), INTENT(INOUT)     :: self
type(EMsoft_T), INTENT(INOUT)       :: EMsoft
character(fnlen), INTENT(INOUT)     :: progname
type(HDFnames_T), INTENT(INOUT)     :: HDFnames

type(MCfile_T)                      :: MCFT
type(MPfile_T)                      :: MPFT
type(HDF_T)                         :: HDF
type(HDFnames_T)                    :: saveHDFnames
type(so3_T)                         :: SO
type(IO_T)                          :: Message
type(Quaternion_T)                  :: quat
type(QuaternionArray_T)             :: qAR
type(memory_T)                      :: mem

type(EBSDmasterNameListType)        :: mpnl
type(MCOpenCLNameListType)          :: mcnl

logical                             :: verbose
character(fnlen)                    :: fname, nmldeffile
integer(kind=irg)                   :: numangles, istat, io_int(1)
type(FZpointd),pointer              :: FZtmp
type(r_T)                           :: rr

call MPFT%setModality('EBSD')
nmldeffile = trim(EMsoft%nmldeffile)

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()

associate( enl => self%nml, EBSDdetector => self%EBSD%det, EBSDMCdata => MCFT%MCDT )

! set the HDF group names for this program
HDFnames = HDFnames_T()

! 1. read the angle array from file and randomize the order using a Fisher-Yates shuffle
verbose = .TRUE.

call setRotationPrecision('d')

fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%anglefile))
call SO%nullifyList()
call SO%getOrientationsfromFile(fname)
numangles = SO%getListCount('FZ')
call SO%listtoQuaternionArray( qAR )
call SO%delete_FZlist()
io_int(1) = numangles
call Message%WriteValue(' Number of orientations read from file: ', io_int, 1)

! randomize the order of the orientations (this is needed to remove any 
! "spatial" correlations in the set of orientations generated by cubochoric sampling)
call Message%printMessage(' Randomizing orientation order (Fisher-Yates + Mersenne)')
call SO%randomizeQuaternionArray( qAR, numangles )


! 2. read the Monte Carlo data file (HDF format)
call HDFnames%set_ProgramData(SC_MCOpenCL)
call HDFnames%set_NMLlist(SC_MCCLNameList)
call HDFnames%set_NMLfilename(SC_MCOpenCLNML)
fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%energyfile))
call MCFT%setFileName(fname)
call MCFT%readMCfile(HDF, HDFnames, getAccume=.TRUE.)
mcnl = MCFT%getnml()

! 3. read EBSD master pattern file (HDF format)
call HDFnames%set_ProgramData(SC_EBSDmaster)
call HDFnames%set_NMLlist(SC_EBSDmasterNameList)
call HDFnames%set_NMLfilename(SC_EBSDmasterNML)
call HDFnames%set_Variable(SC_MCOpenCL)

fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfile))
call MPFT%setFileName(fname)
call MPFT%readMPfile(HDF, HDFnames, mpnl, getmLPNH=.TRUE., getmLPSH=.TRUE.)

! we need to create an EMEBSD style output file with a few extra datasets, 
! so we set the EBSD HDFnames parameters
call HDFnames%set_ProgramData(SC_EBSD)
call HDFnames%set_NMLlist(SC_EBSDNameList)
call HDFnames%set_NMLfilename(SC_EBSD)

! for a regular Euler angle file, we precompute the detector arrays here; for the 'orpcdef' mode
! we compute them later (for each pattern separately)
mem = memory_T()
call mem%alloc(EBSDdetector%rgx, (/ enl%numsx,enl%numsy /), 'EBSDdetector%rgx' )
call mem%alloc(EBSDdetector%rgy, (/ enl%numsx,enl%numsy /), 'EBSDdetector%rgy' )
call mem%alloc(EBSDdetector%rgz, (/ enl%numsx,enl%numsy /), 'EBSDdetector%rgz' )
call mem%alloc(EBSDdetector%accum_e_detector, (/ EBSDMCdata%numEbins,enl%numsx,enl%numsy /), 'EBSDdetector%accum_e_detector' )

! 4. generate detector arrays
! first copy a few of the necessary parameters into the EBSD namelist
call self%EBSD%setnumsx( enl%numsx )
call self%EBSD%setnumsy( enl%numsy )
call self%EBSD%setxpc( enl%xpc )
call self%EBSD%setypc( enl%ypc )
call self%EBSD%setdelta( enl%delta )
call self%EBSD%setthetac( enl%thetac )
call self%EBSD%setL( enl%L )
call self%EBSD%setenergymin( enl%energymin )
call self%EBSD%setenergymax( enl%energymax )

call self%EBSD%GenerateDetector(MCFT, verbose)

  ! perform the pattern computations
call self%ComputeEBSDPCA(EMsoft, MCFT, MPFT, HDF, HDFnames, mpnl, mcnl, mem, numangles, qAR, progname, nmldeffile)

end associate

call closeFortranHDFInterface()

end subroutine EBSDPCA_

!--------------------------------------------------------------------------
subroutine ComputeEBSDPCA_(self, EMsoft, MCFT, MPFT, HDF, HDFnames, mpnl, mcnl, mem, numangles, angles, progname, nmldeffile)
!DEC$ ATTRIBUTES DLLEXPORT :: ComputeEBSDPCA_
!! author: MDG
!! version: 1.0
!! date: 05/08/23
!!
!! compute an energy-weighted EBSD Principal Component decomposition for use
!! in the EMDIPCA indexing program

use mod_EMsoft
use mod_symmetry
use mod_crystallography
use mod_io
use mod_diffraction
use mod_Lambert
use mod_quaternions
use mod_rotations
use mod_so3
use HDF5
use h5im
use h5lt
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
use mod_memory
use mod_image

use, intrinsic :: iso_fortran_env

IMPLICIT NONE

class(EBSDPCA_T), INTENT(INOUT)         :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
type(MCfile_T), INTENT(INOUT)           :: MCFT
type(MPfile_T), INTENT(INOUT)           :: MPFT
type(HDF_T),INTENT(INOUT)               :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames
type(EBSDmasterNameListType),INTENT(INOUT) :: mpnl
type(MCOpenCLNameListType),INTENT(INOUT):: mcnl
type(memory_T), INTENT(INOUT)           :: mem 
integer(kind=irg),INTENT(IN)            :: numangles
type(QuaternionArray_T), INTENT(INOUT)  :: angles
character(fnlen),INTENT(IN)             :: progname
character(fnlen),INTENT(IN)             :: nmldeffile


type(IO_T)                                          :: Message
type(memory_T)                                      :: memth 
type(Timing_T)                                      :: timer
type(e_T)                                           :: eu, eee
type(Quaternion_T)                                  :: qu, ququ
type(q_T)                                           :: quat, qqq
type(r_T)                                           :: ro
type(SpaceGroup_T)                                  :: SG
type(Cell_T)                                        :: cell

integer(kind=irg)                                   :: num,ierr,irec,istat, jpar(7), SGnum, nlines
integer(kind=irg),parameter                         :: iunit = 40
integer(kind=irg),parameter                         :: iunitexpt = 41
integer(kind=irg),parameter                         :: iunitdict = 42
character(fnlen)                                    :: info ! info about the GPU
real(kind=dbl),parameter                            :: nAmpere = 6.241D+18   ! Coulomb per second


integer(kind=irg)                                   :: Ne,Nd,L,totnumexpt,numdictsingle,numexptsingle,imght,imgwd,nnk,numE,&
                                                       recordsize, fratio, cratio, fratioE, cratioE, iii, itmpexpt, hdferr, &
                                                       nsig, numk, recordsize_correct, patsz, tickstart, tickstart2, tock, &
                                                       npy, sz(3), jjj
real(kind=sgl),allocatable                          :: mask(:,:)
real(kind=sgl),allocatable                          :: binned(:,:)
integer(kind=irg),allocatable                       :: acc_array(:,:), ppend(:), ppendE(:)
real(kind=sgl),allocatable                          :: mLPNH(:,:,:),mLPSH(:,:,:),accum_e_MC(:,:,:)
real(kind=sgl),allocatable                          :: patternintd(:,:), lp(:), cp(:), EBSDpat(:,:)
integer(kind=irg),allocatable                       :: patterninteger(:,:), patternad(:,:), EBSDpint(:,:), kij(:,:)
real(kind=dbl),allocatable                          :: rdata(:,:), fdata(:,:), rrdata(:,:), ffdata(:,:), ksqarray(:,:), & 
                                                       klist(:,:), dbinned(:,:), pcs(:,:), svals(:), eulers(:,:)
complex(kind=dbl),allocatable                       :: hpmask(:,:)
complex(C_DOUBLE_COMPLEX),allocatable               :: inp(:,:), outp(:,:)
real(kind=dbl)                                      :: w, Jres
integer(kind=irg)                                   :: dims(2)
character(11)                                       :: dstr
character(15)                                       :: tstrb
character(15)                                       :: tstre
character(3)                                        :: vendor
character(fnlen, KIND=c_char),allocatable,TARGET    :: stringarray(:)
character(fnlen)                                    :: groupname, dataset, fname, clname, ename, sourcefile, &
                                                       datagroupname, dictfile, attname, IPFmode
integer(hsize_t)                                    :: expwidth, expheight
integer(hsize_t),allocatable                        :: iPhase(:), iValid(:)
integer(c_size_t),target                            :: slength
integer(c_int)                                      :: numd, nump
type(C_PTR)                                         :: planf, HPplanf, HPplanb
integer(HSIZE_T)                                    :: dims2(2), offset2(2), dims3(3), offset3(3)

integer(kind=irg)                                   :: i,j,ii,jj,kk,ll,mm,pp,qq, cn, dn, totn, icnt
integer(kind=irg)                                   :: FZcnt, pgnum, io_int(4), ncubochoric, pc, ecpipar(4)
type(FZpointd),pointer                              :: FZlist, FZtmp
integer(kind=irg),allocatable                       :: indexlist(:),indexarray(:),indexmain(:,:),indextmp(:,:)
real(kind=sgl)                                      :: dmin,voltage,scl,ratio, mi, ma, ratioE, io_real(2), tstart, tmp, &
                                                       totnum_el, vlen, tstop, ttime
real(kind=dbl)                                      :: prefactor
character(fnlen)                                    :: xtalname, IPFmapfile
integer(kind=irg)                                   :: binx,biny,TID,nthreads,Emin,Emax, iiistart, iiiend, jjend
real(kind=sgl)                                      :: sx,dx,dxm,dy,dym,rhos,x,projweight, dp, mvres, nel, emult
real(kind=sgl)                                      :: dc(3),ixy(2),bindx, MCsig, WD, fpar1(1), fpar2(2)
real(kind=dbl)                                      :: sm, smsqr, mn, sdev
integer(kind=irg)                                   :: nix,niy,nixp,niyp
real(kind=sgl)                                      :: euler(3)
integer(kind=irg)                                   :: indx
integer(kind=irg)                                   :: correctsize
logical                                             :: f_exists, init, ROIselected, Clinked, cancelled, isTKD = .FALSE., &
                                                       isEBSD = .FALSE., isECP = .FALSE., switchwfoff, verbose, overwrite=.TRUE.

integer(kind=irg)                                   :: ipar(10)

character(fnlen),ALLOCATABLE                        :: MessageLines(:)
integer(kind=irg)                                   :: NumLines
character(fnlen)                                    :: TitleMessage, exectime
character(100)                                      :: c
character(1000)                                     :: charline
character(3)                                        :: stratt
character(fnlen)                                    :: progdesc, TIFF_filename, datafile, attributename, HDF_FileVersion
character(4)                                        :: nstr
character(fnlen,kind=c_char)                        :: line2(1)

! parameters for BLAS:dsyrk() and LAPACK:dsyevd() calls 
character(1)                                        :: UPLO 
character(1)                                        :: TRANS
character(1)                                        :: JOBZ
integer(kind=irg)                                   :: NNNN 
integer(kind=irg)                                   :: KKKK
real(kind=dbl)                                      :: alpha 
real(kind=dbl),allocatable                          :: dict(:,:)
real(kind=dbl),allocatable                          :: WEIG(:)
real(kind=dbl),allocatable                          :: WORK(:)
integer(kind=irg)                                   :: LDA 
integer(kind=irg)                                   :: LWORK
integer(kind=irg),allocatable                       :: IWORK(:)
integer(kind=irg)                                   :: LIWORK
real(kind=dbl)                                      :: BETA
real(kind=dbl),allocatable                          :: covmat(:,:)
integer(kind=irg)                                   :: LINFO
integer(kind=irg)                                   :: LDC
integer(kind=irg)                                   :: iar(3,2)

! declare variables for use in object oriented image module
integer                                             :: iostat
character(len=128)                                  :: iomsg
logical                                             :: isInteger
type(image_t)                                       :: im
integer(int8)                                       :: i8 (3,4)
integer(int8), allocatable                          :: TIFF_image(:,:)

call setRotationPrecision('d')

associate( dinl=>self%nml, MPDT=>MPFT%MPDT, MCDT=>MCFT%MCDT, det=>self%EBSD%det, enl=>self%EBSD%nml )

! initialize the multi-threaded memory allocation classes
memth = memory_T( nt = dinl%nthreads )

!=====================================================
! make sure the minimum energy is set smaller than the maximum
!=====================================================
if (dinl%energymin.gt.dinl%energymax) then
    call Message%printMessage(' Minimum energy is larger than maximum energy; please correct input file')
    stop
end if

!=====================================================
! get the indices of the minimum and maximum energy
!=====================================================
Emin = nint((dinl%energymin - mcnl%Ehistmin)/mcnl%Ebinsize) +1
if (Emin.lt.1)  Emin=1
if (Emin.gt.MCDT%numEbins)  Emin=MCDT%numEbins

Emax = nint((dinl%energymax - mcnl%Ehistmin)/mcnl%Ebinsize) + 1
if (Emax .lt. 1) Emax = 1
if (Emax .gt. MCDT%numEbins) Emax = MCDT%numEbins

sz = shape(MPDT%mLPNH)
numE = sz(3)

! intensity prefactor
nel = float(mcnl%totnum_el) * float(mcnl%multiplier)
emult = nAmpere * 1e-9 / nel  ! multiplicative factor to convert MC data to an equivalent incident beam of 1 nanoCoulomb
! intensity prefactor  (redefined by MDG, 3/23/18)
prefactor = emult * dinl%beamcurrent * dinl%dwelltime * 1.0D-6

npy = mpnl%npx
binx = enl%numsx 
biny = enl%numsy 
call mem%alloc(mLPNH, (/ mpnl%npx,npy,MCDT%numEbins /), 'mLPNH', startdims= (/ -mpnl%npx,-npy,1 /))
call mem%alloc(mLPSH, (/ mpnl%npx,npy,MCDT%numEbins /), 'mLPSH', startdims= (/ -mpnl%npx,-npy,1 /))
call mem%alloc(accum_e_MC, (/ MCDT%numEbins,dinl%numsx,dinl%numsy /), 'accum_e_MC')
accum_e_MC = det%accum_e_detector
mLPNH = MPDT%mLPNH
mLPSH = MPDT%mLPSH
deallocate(det%accum_e_detector, MPDT%mLPNH, MPDT%mLPSH)

FZcnt = numangles
L = dinl%numsx*dinl%numsy/dinl%binning**2
w = dinl%hipassw

!=========================================
! ALLOCATION AND INITIALIZATION OF ARRAYS
!=========================================
call Message%printMessage(' --> Allocating various arrays for covariance matrix computation')

call mem%alloc(dict, (/ L, FZcnt /), 'dict', initval = 0.D0)
call mem%alloc(covmat, (/ L, L /), 'covmat', initval = 0.D0)
call mem%alloc(mask, (/ binx,biny /), 'mask', initval = 1.0)
! call mem%alloc(rdata, (/ binx,biny /), 'rdata', initval = 0.D0) 
! call mem%alloc(fdata, (/ binx,biny /), 'fdata', initval = 0.D0)

if (dinl%maskpattern.eq.'y') then
  do ii = 1,biny
      do jj = 1,binx
          if((ii-biny/2)**2 + (jj-binx/2)**2 .ge. dinl%maskradius**2) then
              mask(jj,ii) = 0.0
          end if
      end do
  end do
end if

!=====================================================
! MAIN COMPUTATIONAL LOOP (finally...) somewhat similar to EMDI Driver routine !!!
!=====================================================

timer = Timing_T()
tstrb = timer%getTimeString()
dstr = timer%getDateString()

call timer%makeTimeStamp()
call timer%Time_tick(1)
call timer%Time_tick(2)

call OMP_setNThreads(dinl%nthreads)

! define the jpar array of integer parameters
jpar(1) = 1 ! dinl%binning
jpar(2) = dinl%numsx
jpar(3) = dinl%numsy
jpar(4) = mpnl%npx
jpar(5) = npy
jpar(6) = MCDT%numEbins
jpar(7) = numE

verbose = .FALSE.

! open a temporary file that will contain the original pre-processed dictionary patterns
! for re-use later on when we compute the PCA projected patterns.



! compute the entire array of dictionary patterns as a L * FZcnt array of doubles
dict = 0.D0



!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID,iii,jj,ll,mm,pp,ierr,io_int, vlen, tock, ttime) &
!$OMP& PRIVATE(binned, ma, mi, patternintd, patterninteger, patternad, qu, ro, quat, icnt) &
!$OMP& PRIVATE(eee, qqq, ququ, sm, smsqr, mn, sdev, dbinned)

TID = OMP_GET_THREAD_NUM()

! allocate the local arrays that are used by each thread
call memth%alloc(patterninteger, (/ binx,biny /), 'patterninteger', TID=TID, initval = 0)
call memth%alloc(patternad, (/ binx,biny /), 'patternad', TID=TID, initval = 0) 
call memth%alloc(patternintd, (/ binx,biny /), 'patternintd', TID=TID, initval = 0.0)
call memth%alloc(binned, (/ binx,biny /), 'binned', TID=TID, initval = 0.0)
call memth%alloc(dbinned, (/ binx,biny /), 'dbinned', TID=TID, initval = 0.D0)

if (TID.eq.0) then
  io_int(1) = OMP_GET_NUM_THREADS()
  call Message%WriteValue(' actual number of OpenMP threads  = ', io_int, 1)
end if


!$OMP DO SCHEDULE(DYNAMIC)
dictionaryloop: do ii = 1,FZcnt
   binned = 0.0
   qu = angles%getQuatfromArray( ii ) 
   call self%EBSD%CalcEBSDPatternSingleFull(jpar,qu,accum_e_MC,mLPNH,mLPSH,self%EBSD%det%rgx,&
                                           self%EBSD%det%rgy,self%EBSD%det%rgz,binned,Emin,Emax,mask,prefactor)

   if (dinl%scalingmode .eq. 'gam') then
     binned = binned**dinl%gammavalue
   end if

! hi pass filtering
   ! rdata = dble(binned)
   ! fdata = HiPassFilter(rdata,dims,w, init=.TRUE.)
   ! binned = sngl(fdata)

! adaptive histogram equalization
   ma = maxval(binned)
   mi = minval(binned)

   patternintd = ((binned - mi)/ (ma-mi))
   patterninteger = nint(patternintd*255.0)
   patternad =  adhisteq(dinl%nregions,binx,biny,patterninteger)
   dbinned = dble(patternad) * mask

! zero mean and unit standard deviation
   sm = sum(dbinned)
   smsqr = sum(dbinned**2)
   mn = sm / dble(L)
   sdev = sqrt( (smsqr - sm*sm/dble(L)) / dble(L-1) )
   if (sdev.eq.0.0) sdev = 1.0
   dbinned = (dbinned - mn) / sdev 
! and place it in the array
   dict(1:L,ii) = reshape( dbinned, (/ L /) )
! message every now and then
   if (mod(ii,1000).eq.0) then 
       io_int(1) = ii
       io_int(2) = FZcnt
       call Message%WriteValue(' completed pattern ',io_int,2,"(I6,' out of ',I8)")
   end if
end do dictionaryloop
!$OMP END DO

call memth%dealloc(binned, 'binned', TID=TID)
call memth%dealloc(dbinned, 'dbinned', TID=TID)
call memth%dealloc(patterninteger, 'patterninteger', TID=TID)
call memth%dealloc(patternad, 'patternad', TID=TID)
call memth%dealloc(patternintd, 'patternintd', TID=TID)

! and we end the parallel section here (all threads will synchronize).
!$OMP END PARALLEL

call Message%printMessage(' Computing covariance matrix')
! next we compute the covariance array using the BLAS dsyrk routine
LDA = L
LDC = L
NNNN = L
KKKK = FZcnt
ALPHA = 1.D0
BETA = 0.D0 
UPLO = 'U'
TRANS = 'N'
call dsyrk( UPLO, TRANS, NNNN, KKKK, ALPHA, dict, LDA, BETA, covmat, LDC)
covmat = covmat / dble(L-1)
call Message%printMessage('  ---> done')

! next is the SVD step
call Message%printMessage(' Performing SVD decomposition')
! first get the workspace information
JOBZ = 'V'
call mem%alloc(WEIG, (/ L /), 'WEIG')
call mem%alloc(WORK, (/ 1 /), 'WORK')
call mem%alloc(IWORK, (/ 1 /), 'IWORK')
LWORK = -1 
LIWORK = -1 
call dsyevd( JOBZ, UPLO, NNNN, covmat, LDA, WEIG, WORK, LWORK, IWORK, LIWORK, LINFO)
if (LINFO.eq.0) then 
  LWORK = WORK(1)
  LIWORK = IWORK(1)
  call mem%dealloc(WORK, 'WORK')
  call mem%dealloc(IWORK, 'IWORK')
  call mem%alloc(WORK, (/ LWORK /), 'WORK')
  call mem%alloc(IWORK, (/ LIWORK /), 'IWORK')
else 
  call Message%printError('LAPACK:dsyevd','Error determining LWORK/LIWORK parameters')
end if 
call dsyevd( JOBZ, UPLO, NNNN, covmat, LDA, WEIG, WORK, LWORK, IWORK, LIWORK, LINFO)
call mem%dealloc(WORK, 'WORK')
call mem%dealloc(IWORK, 'IWORK')
call Message%printMessage('  ---> done')

! reorder the eigenvalues and eigenvectors to be in decreasing order 
WEIG = WEIG(L:1:-1)
covmat = covmat(:,L:1:-1)

! singular values 
call mem%alloc(svals, (/ L /), 'svals')
svals = sqrt(abs(WEIG))

! Apply whitening by dividing the components by the singular values and scaling by sqrt(n_samples)
do i=1,L 
  covmat(:,i) = covmat(:,i) / svals(i)
end do 
covmat = covmat * sqrt(dble(FZcnt))
covmat = transpose(covmat)

! determine the principal components 
call Message%printMessage(' Projecting dictionary patterns')
call mem%alloc(pcs, (/ L, FZcnt /), 'pcs')
pcs = matmul( covmat, dict )
call Message%printMessage('  ---> done')

!========================
!======Output============
!========================
! This output file must be compatible with the standard dictionary file produced
! by the EMEBSD program, with the same data set names, only with additional 
! data sets for the PCA part.   That way, we should be able to use EMDI as before
! in static indexing mode...

! get the crystal structure data
call cell%getCrystalData(mcnl%xtalname, SG, EMsoft, useHDF=HDF)

! in the EBSDPCA.h5 output file, we store the usual stuff, and in terms of 
! data arrays:  pcs, covmat, svals, FZcnt, qAR 
! Create a new file using the default properties.
datafile = EMsoft%generateFilePath('EMdatapathname', dinl%datafile)

call Message%printMessage(' Creating HDF5 file '//trim(datafile))
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
datagroupname = trim(HDFnames%get_ProgramData()) ! 'EBSD' 
call HDF%writeEMheader(EMsoft,dstr, tstrb, tstre, progname, datagroupname)

! add the CrystalData group at the top level of the file
call cell%addXtalDataGroup(SG, EMsoft, HDF)

! create a namelist group to write all the namelist files into
hdferr = HDF%createGroup(HDFnames%get_NMLfiles())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup NMLfiles', hdferr)

! read the text file and write the array to the file
dataset = SC_EMEBSDNML
hdferr = HDF%writeDatasetTextFile(dataset, nmldeffile)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetTextFile EMEBSDNML', hdferr)

call HDF%pop()

! create a NMLparameters group to write all the namelist entries into
hdferr = HDF%createGroup(HDFnames%get_NMLparameters())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup NMLparameters', hdferr)

call self%writeHDFNameList_(HDF, HDFnames)

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
call mem%alloc(stringarray, (/ 1 /), 'stringarray')
stringarray(1)= trim(mcnl%xtalname)
hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetStringArray xtalname', hdferr)

dataset = SC_numangles
hdferr = HDF%writeDatasetInteger(dataset, numangles)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetInteger numangles', hdferr)

! and add the Euler angles to the output file
call mem%alloc(eulers, (/ 3,numangles /), 'eulerangles')
do i=1,numangles
  qu = angles%getQuatfromArray(i)
  call quat%q_setd( qu%get_quatd() )
  eu = quat%qe()
  eulers(1:3,i) = sngl( eu%e_copyd() )
end do
dataset = SC_Eulerangles
hdferr = HDF%writeDatasetFloatArray(dataset, sngl(eulers/dtor), 3, numangles)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetFloatArray2D Eulerangles', hdferr)

! =====================================================
! end of HDF_FileVersion = 4.0 and above write statements
! =====================================================

! normally we would store the dict array as EBSDpatterns in the output file;
! in this case, we store the PCA principal components of all the patterns in 
! an array of the same size and name, but we also store the covariance matrix
! and the singular values...
dataset = SC_EBSDpatterns
hdferr = HDF%writeDatasetFloatArray(dataset, sngl(pcs), L, FZcnt)
if (hdferr.ne.0) call HDF%error_check('writeDatasetDoubleArray pcs', hdferr)

dataset = 'CovarianceMatrixWhitened'
hdferr = HDF%writeDatasetFloatArray(dataset, sngl(covmat), L, L)
if (hdferr.ne.0) call HDF%error_check('writeDatasetDoubleArray covmat', hdferr)

dataset = 'SingularValues'
hdferr = HDF%writeDatasetFloatArray(dataset, sngl(svals), L)
if (hdferr.ne.0) call HDF%error_check('writeDatasetDoubleArray svals', hdferr)

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
call HDF%popall()


end associate 

end subroutine ComputeEBSDPCA_


end module mod_EBSDPCA
