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

module mod_EBSDkin
  !! author: MDG
  !! version: 1.0
  !! date: 11/18/21
  !!
  !! class definition for the EMEBSDkin program

use mod_kinds
use mod_global
use mod_quaternions
use mod_MPfiles

IMPLICIT NONE
private

! namelist for the EMEBSDkin program
type, public :: EBSDkinNameListType
  integer(kind=irg)       :: numsx
  integer(kind=irg)       :: numsy
  integer(kind=irg)       :: binning
  integer(kind=irg)       :: nthreads
  integer(kind=irg)       :: maskradius
  real(kind=sgl)          :: L
  real(kind=sgl)          :: sigma
  real(kind=sgl)          :: omega
  real(kind=sgl)          :: thetac
  real(kind=sgl)          :: delta
  real(kind=sgl)          :: xpc
  real(kind=sgl)          :: ypc
  real(kind=sgl)          :: gammavalue
  real(kind=sgl)          :: axisangle(4)
  real(kind=sgl)          :: alphaBD
  real(kind=dbl)          :: Ftensor(3,3)
  character(1)            :: makedictionary
  character(1)            :: applyDeformation
  character(1)            :: maskpattern
  character(3)            :: eulerconvention
  character(3)            :: outputformat
  character(5)            :: bitdepth
  character(fnlen)        :: anglefile
  character(fnlen)        :: anglefiletype
  character(fnlen)        :: masterfile
  character(fnlen)        :: energyfile
  character(fnlen)        :: datafile
end type EBSDkinNameListType

! angles, patterns centers, and deformation tensor arrays
type EBSDkinAnglePCDefType
  real(kind=sgl),allocatable  :: pcs(:,:)
  real(kind=sgl),allocatable  :: deftensors(:,:,:)
  real(kind=dbl),allocatable  :: pcfield(:,:,:)
  real(kind=dbl),allocatable  :: deftensorfield(:,:,:,:)
end type EBSDkinAnglePCDefType

type EBSDkinPixel
  real(kind=sgl),allocatable  :: lambdaEZ(:,:)
  real(kind=dbl)              :: dc(3) ! direction cosine in sample frame
  real(kind=dbl)              :: cfactor
end type EBSDkinPixel

type EBSDkinDetectorType
  real(kind=sgl),allocatable  :: rgx(:,:), rgy(:,:), rgz(:,:)  ! auxiliary detector arrays needed for interpolation
  type(EBSDkinPixel),allocatable :: detector(:,:)
end type EBSDkinDetectorType


! class definition
type, public :: EBSDkin_T
private
  character(fnlen)                    :: nmldeffile = 'EMEBSDkin.nml'
  character(fnlen)                    :: MPfile
  character(fnlen)                    :: xtalname
  integer(kind=irg)                   :: npx
  type(EBSDkinNameListType), public   :: nml
  type(EBSDkinDetectorType), public   :: det
  real(kind=sgl),allocatable          :: mLPNH(:,:)
  real(kind=sgl),allocatable          :: mLPSH(:,:)
  real(kind=sgl),allocatable          :: masterSPNH(:,:)
  real(kind=sgl),allocatable          :: masterSPSH(:,:)

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: readMPfile_
  procedure, pass(self) :: setFileName_
  procedure, pass(self) :: EBSDkin_
  procedure, pass(self) :: EBSDkinreadorpcdef_
  procedure, pass(self) :: GenerateDetector_
  procedure, pass(self) :: GeneratemyEBSDkinDetector_
  procedure, pass(self) :: ComputeEBSDkinPatterns_
  procedure, pass(self) :: ComputedeformedEBSDkinpatterns_
  procedure, pass(self) :: CalcEBSDkinPatternSingleFull_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: readMPfile => readMPfile_
  generic, public :: setFileName => setFileName_
  generic, public :: EBSDkin => EBSDkin_
  generic, public :: EBSDkinreadorpcdef => EBSDkinreadorpcdef_
  generic, public :: GenerateDetector => GenerateDetector_
  generic, public :: GeneratemyEBSDkinDetector => GeneratemyEBSDkinDetector_
  generic, public :: ComputeEBSDkinPatterns => ComputeEBSDkinPatterns_
  generic, public :: ComputedeformedEBSDkinpatterns => ComputedeformedEBSDkinpatterns_
  generic, public :: CalcEBSDkinPatternSingleFull => CalcEBSDkinPatternSingleFull_

end type EBSDkin_T

! the constructor routine for this class
interface EBSDkin_T
  module procedure EBSDkin_constructor
end interface EBSDkin_T

contains

!--------------------------------------------------------------------------
type(EBSDkin_T) function EBSDkin_constructor( nmlfile ) result(EBSDkin)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDkin_constructor
!! author: MDG
!! version: 1.0
!! date: 11/18/21
!!
!! constructor for the EBSDkin_T Class; reads the name list

IMPLICIT NONE

character(fnlen)    :: nmlfile

call EBSDkin%readNameList(nmlfile)

end function EBSDkin_constructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 11/18/21
!!
!! read the namelist from an nml file for the EBSDkin_T Class

use mod_io

IMPLICIT NONE

class(EBSDkin_T), INTENT(INOUT)             :: self
character(fnlen),INTENT(IN)                 :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)                 :: initonly
 !! fill in the default values only; do not read the file

type(IO_T)                                  :: Message
logical                                     :: skipread = .FALSE.

integer(kind=irg)       :: numsx
integer(kind=irg)       :: numsy
integer(kind=irg)       :: binning
integer(kind=irg)       :: nthreads
integer(kind=irg)       :: maskradius
real(kind=sgl)          :: L
real(kind=sgl)          :: thetac
real(kind=sgl)          :: sigma
real(kind=sgl)          :: omega
real(kind=sgl)          :: delta
real(kind=sgl)          :: xpc
real(kind=sgl)          :: ypc
real(kind=sgl)          :: alphaBD
real(kind=dbl)          :: Ftensor(3,3)
character(1)            :: makedictionary
character(1)            :: applyDeformation
character(1)            :: maskpattern
character(3)            :: eulerconvention
character(3)            :: outputformat
character(5)            :: bitdepth
character(fnlen)        :: anglefile
character(fnlen)        :: anglefiletype
character(fnlen)        :: masterfile
character(fnlen)        :: datafile

! define the IO namelist to facilitate passing variables to the program.
namelist  / EBSDkindata / L, thetac, delta, numsx, numsy, xpc, ypc, anglefile, eulerconvention, masterfile, bitdepth, &
                       datafile, binning, alphaBD, nthreads, outputformat, maskpattern, sigma, omega, &
                       applyDeformation, Ftensor, anglefiletype, makedictionary, maskradius


! set the input parameters to default values (except for xtalname, which must be present)
numsx           = 0             ! [dimensionless]
numsy           = 0             ! [dimensionless]
binning         = 1             ! binning mode  (1, 2, 4, or 8)
L               = 20000.0       ! [microns]
nthreads        = 1             ! number of OpenMP threads
thetac          = 0.0           ! [degrees]
sigma           = 70.0          ! [degrees]
omega           = 0.0          ! [degrees]
delta           = 25.0          ! [microns]
xpc             = 0.0           ! [pixels]
ypc             = 0.0           ! [pixels]
alphaBD         = 0.0           ! transfer lens barrel distortion parameter
maskradius      = 240           ! mask radius
Ftensor         = reshape( (/ 1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 1.D0 /), (/ 3,3 /) )
makedictionary  = 'n'
applyDeformation = 'n'          ! should we apply a deformation tensor to the unit cell?
maskpattern     = 'n'           ! 'y' or 'n' to include a circular mask
eulerconvention = 'tsl'         ! convention for the first Euler angle ['tsl' or 'hkl']
outputformat    = 'gui'         ! output format for 'bin' or 'gui' use
bitdepth        = '8bit'        ! format for output; '8char' for [0..255], '##int' for integers, 'float' for floats
! the '##int' notation stands for the actual bitdepth; all values are stored as 32bit integers, but they are scaled
! from the float values to a maximum that is given by the first two digits, which indicate the bit depth; so, valid
! values would be '10int' for a 10-bit integer scale, '16int' for a 16-bit integer scale, and so on.
anglefile       = 'undefined'   ! filename
anglefiletype   = 'orientations'! 'orientations' or 'orpcdef'
masterfile      = 'undefined'   ! filename
datafile        = 'undefined'   ! output file name

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=EBSDkindata)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries

! we no longer require the energyfile parameter, but for backwards compatibility
! we still allow the user to include it (it doesn't do anything though)
! if (trim(energyfile).eq.'undefined') then
!  call FatalError('GetEBSDkinNameList:',' energy file name is undefined in '//nmlfile)
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
self%nml%maskradius = maskradius
self%nml%L = L
self%nml%nthreads = nthreads
self%nml%thetac = thetac
self%nml%sigma = sigma
self%nml%omega = omega
self%nml%delta = delta
self%nml%xpc = xpc
self%nml%ypc = ypc
self%nml%alphaBD = alphaBD
self%nml%Ftensor = Ftensor
self%nml%makedictionary = makedictionary
self%nml%applyDeformation = applyDeformation
self%nml%maskpattern = maskpattern
self%nml%eulerconvention = eulerconvention
self%nml%outputformat = outputformat
self%nml%bitdepth = bitdepth
self%nml%anglefile = anglefile
self%nml%anglefiletype = anglefiletype
self%nml%masterfile = masterfile
self%nml%datafile = datafile

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 11/18/21
!!
!! pass the namelist for the EBSDkin_T Class to the calling program

IMPLICIT NONE

class(EBSDkin_T), INTENT(INOUT)          :: self
type(EBSDkinNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine setFileName_(self, fname) 
!DEC$ ATTRIBUTES DLLEXPORT :: setFileName_
!! author: MDG
!! version: 1.0
!! date: 1/19/22
!!
!! set the MPfile parameter

IMPLICIT NONE

class(EBSDkin_T), INTENT(INOUT)          :: self
character(fnlen), INTENT(IN)             :: fname

self%MPfile = trim(fname) 

end subroutine setFileName_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames, isTKD)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG
!! version: 1.0
!! date: 11/18/21
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants

use ISO_C_BINDING

IMPLICIT NONE

class(EBSDkin_T), INTENT(INOUT)         :: self
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames
logical,OPTIONAL,INTENT(IN)             :: isTKD

integer(kind=irg),parameter             :: n_int = 5, n_real = 8
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
io_int = (/ enl%numsx, enl%numsy, enl%binning, enl%nthreads, enl%maskradius /)
intlist(1) = 'numsx'
intlist(2) = 'numsy'
intlist(3) = 'binning'
intlist(4) = 'nthreads'
intlist(5) = 'maskradius'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ enl%L, enl%thetac, enl%delta, enl%xpc, enl%ypc, enl%alphaBD, enl%sigma, enl%omega /)
reallist(1) = 'L'
reallist(2) = 'thetac'
reallist(3) = 'delta'
reallist(4) = 'xpc'
reallist(5) = 'ypc'
reallist(6) = 'alphaBD'
reallist(7) = 'sigma'
reallist(8) = 'omega'
call HDF%writeNMLreals(io_real, reallist, n_real)

! a 3x3 matrix
dataset = SC_Ftensor
hdferr = HDF%writeDatasetDoubleArray(dataset, enl%Ftensor, 3, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create Ftensor dataset', hdferr)

! write all the strings
dataset = SC_maskpattern
line2(1) = trim(enl%maskpattern)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create maskpattern dataset', hdferr)

dataset = SC_makedictionary
line2(1) = trim(enl%makedictionary)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create makedictionary dataset', hdferr)

dataset = SC_applyDeformation
line2(1) = trim(enl%applyDeformation)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create applyDeformation dataset', hdferr)

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
recursive subroutine readMPfile_(self, HDF, HDFnames, mpnl, getmLPNH, getmLPSH, &
                                 getmasterSPNH, getmasterSPSH)
!DEC$ ATTRIBUTES DLLEXPORT :: readMPfile_
!! author: MDG
!! version: 1.0
!! date: 02/17/20
!!
!! read a kinematical EBSD Master Pattern file into the correct namelist and data structure

use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_io
use stringconstants
use ISO_C_BINDING

IMPLICIT NONE

class(EBSDkin_T), INTENT(INOUT)                  :: self
type(HDF_T), INTENT(INOUT)                       :: HDF
type(HDFnames_T), INTENT(INOUT)                  :: HDFnames
type(EBSDkinNameListType), INTENT(INOUT)         :: mpnl
logical,INTENT(IN),OPTIONAL                      :: getmLPNH
logical,INTENT(IN),OPTIONAL                      :: getmLPSH
logical,INTENT(IN),OPTIONAL                      :: getmasterSPNH
logical,INTENT(IN),OPTIONAL                      :: getmasterSPSH

type(HDFnames_T)                                 :: saveHDFnames
type(IO_T)                                       :: Message
character(fnlen)                                 :: infile, groupname, datagroupname, dataset
logical                                          :: stat, readonly, g_exists, f_exists, FL
integer(kind=irg)                                :: ii, nlines, istat, hdferr
integer(kind=irg),allocatable                    :: iarray(:)
real(kind=sgl),allocatable                       :: farray(:)
integer(HSIZE_T)                                 :: dims(1), dims2(2), dims3(3), offset3(3), dims4(4)
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)

! set the HDF group names to those for the kinematical master pattern file 
saveHDFnames = HDFnames 
call HDFnames%set_ProgramData(SC_EMkinematical)
call HDFnames%set_NMLlist(SC_EMkinematicalNameList)
call HDFnames%set_NMLfilename(SC_EMkinematicalNML)

associate( MPDT => self%nml )

! we assume that the calling program has opened the HDF interface
inquire(file=trim(self%MPfile), exist=f_exists)

if (.not.f_exists) then
  call Message%printError('readMPfile','Master Pattern input file does not exist')
end if

! is this a proper HDF5 file ?
call h5fis_hdf5_f(trim(self%MPfile), stat, hdferr)

if (stat.eqv..FALSE.) then ! the file exists, so let's open it an first make sure it is an EBSD dot product file
   call Message%printError('readMPfile','This is not a proper HDF5 file')
end if

! open the Master Pattern file
readonly = .TRUE.
hdferr =  HDF%openFile(self%MPfile, readonly)

! check whether or not the MC file was generated using DREAM.3D
! this is necessary so that the proper reading of fixed length vs. variable length strings will occur.
! this test sets a flag in side the HDFsupport module so that the proper reading routines will be employed

hdferr = HDF%openGroup(HDFnames%get_EMheader())
datagroupname = trim(HDFnames%get_ProgramData())
call H5Lexists_f(HDF%getobjectID(),trim(datagroupname),g_exists, hdferr)
if (.not.g_exists) then
  call Message%printError('readMPfile','This HDF file does not contain Master Pattern header data')
end if

hdferr = HDF%openGroup(HDFnames%get_ProgramData())
FL = .FALSE.
datagroupname = 'FixedLength'
FL = HDF%CheckFixedLengthflag(datagroupname)
if (FL.eqv..TRUE.) then
  call Message%printMessage('Input file was generated by a program using fixed length strings')
end if
call HDF%pop()
call HDF%pop()

!====================================
! make sure this is a Master Pattern file
!====================================
hdferr = HDF%openGroup(HDFnames%get_NMLfiles())
dataset = trim(HDFnames%get_NMLfilename())
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists.eqv..FALSE.) then
    call HDF%popall()
    call Message%printError('readMPfile','this is not a valid Kinematical Master Pattern file')
end if
call HDF%pop()

!====================================
! read all NMLparameters group datasets
!====================================
hdferr = HDF%openGroup(HDFnames%get_NMLparameters())
hdferr = HDF%openGroup(HDFnames%get_NMLlist())

dataset = SC_xtalname
    call HDF%readDatasetStringArray(dataset, nlines, hdferr, stringarray)
    self%xtalname = trim(stringarray(1))
    deallocate(stringarray)

dataset = 'nx'
    call HDF%readDatasetInteger(dataset, hdferr, self%npx)

! and close the NMLparameters group
call HDF%pop()
call HDF%pop()

!====================================
!====================================

! open the Master Pattern data group
hdferr = HDF%openGroup(HDFnames%get_EMData())
hdferr = HDF%openGroup(HDFnames%get_ProgramData())

if (present(getmLPNH)) then
  if (getmLPNH.eqv..TRUE.) then
    dataset = 'masterNH'
    call HDF%readDatasetFloatArray(dataset, dims2, hdferr, self%mLPNH)
  end if
end if

if (present(getmLPSH)) then
  if (getmLPSH.eqv..TRUE.) then
    dataset = 'masterSH'
    call HDF%readDatasetFloatArray(dataset, dims2, hdferr, self%mLPSH)
  end if
end if

if (present(getmasterSPNH)) then
  if (getmasterSPNH.eqv..TRUE.) then
    dataset = 'stereoNH'
    call HDF%readDatasetFloatArray(dataset, dims2, hdferr, self%masterSPNH)
  end if
end if

if (present(getmasterSPSH)) then
  if (getmasterSPSH.eqv..TRUE.) then
    dataset = 'stereoSH'
    call HDF%readDatasetFloatArray(dataset, dims2, hdferr, self%masterSPSH)
  end if
end if

! and close the HDF5 Master Pattern file
call HDF%popall()

call Message%printMessage(' --> Completed reading master pattern data from '//trim(self%MPfile), frm = "(A/)")

end associate

HDFnames = saveHDFnames

end subroutine readMPfile_

!--------------------------------------------------------------------------
subroutine EBSDkin_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDkin_
!! author: MDG
!! version: 1.0
!! date: 11/18/21
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

class(EBSDkin_T), INTENT(INOUT)     :: self
type(EMsoft_T), INTENT(INOUT)       :: EMsoft
character(fnlen), INTENT(INOUT)     :: progname
type(HDFnames_T), INTENT(INOUT)     :: HDFnames

type(HDF_T)                         :: HDF
type(HDFnames_T)                    :: saveHDFnames
type(so3_T)                         :: SO
type(IO_T)                          :: Message
type(Quaternion_T)                  :: quat
type(QuaternionArray_T)             :: qAR
type(memory_T)                      :: mem

type(MCOpenCLNameListType)          :: mcnl
type(EBSDkinAnglePCDefType)         :: orpcdef

logical                             :: verbose
character(fnlen)                    :: fname, nmldeffile
integer(kind=irg)                   :: numangles, istat
type(FZpointd),pointer              :: FZtmp
type(r_T)                           :: rr

nmldeffile = trim(EMsoft%nmldeffile)

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()

associate( enl => self%nml, EBSDkindetector => self%det )

! set the HDF group names for this program
! saveHDFnames = HDFnames

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
  write (*,*) ' Number of orientations read from file: ', numangles
else if (trim(enl%anglefiletype).eq.'orpcdef') then
! this requires a conversion from the Euler angles in the file to quaternions
! plus storage of the pattern center and deformation tensor arrays
  call self%EBSDkinreadorpcdef(EMsoft, numangles, qAR, orpcdef, verbose)
else
  call Message%printError('EBSDkin','unknown anglefiletype')
end if


! 2. read EBSDkin master pattern file (HDF format)
fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfile))
call self%setFileName(fname)
call self%readMPfile(HDF, HDFnames, enl, getmLPNH=.TRUE., getmLPSH=.TRUE.)

write (*,*) 'range of mLPNH: ', minval(self%mLPNH), maxval(self%mLPNH)

! for a regular Euler angle file, we precompute the detector arrays here; for the 'orpcdef' mode
! we compute them later (for each pattern separately)
if (trim(enl%anglefiletype).eq.'orientations') then
  mem = memory_T()
  call mem%alloc(EBSDkindetector%rgx, (/ enl%numsx,enl%numsy /), 'EBSDkindetector%rgx' )
  call mem%alloc(EBSDkindetector%rgy, (/ enl%numsx,enl%numsy /), 'EBSDkindetector%rgy' )
  call mem%alloc(EBSDkindetector%rgz, (/ enl%numsx,enl%numsy /), 'EBSDkindetector%rgz' )

! 3. generate detector arrays
  call self%GenerateDetector(verbose)

  ! perform the pattern computations
  call self%ComputeEBSDkinPatterns(EMsoft, HDF, HDFnames, mem, numangles, qAR, progname, nmldeffile)
end if

if (trim(enl%anglefiletype).eq.'orpcdef') then
  call self%ComputedeformedEBSDkinPatterns(EMsoft, HDF, HDFnames, numangles, qAR, orpcdef, progname, nmldeffile)
end if

end associate

call closeFortranHDFInterface()

end subroutine EBSDkin_

!--------------------------------------------------------------------------
recursive subroutine EBSDkinreadorpcdef_(self, EMsoft, numangles, qAR, orpcdef, verbose)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDkinreadorpcdef_
!! author: MDG
!! version: 1.0
!! date: 11/18/21
!!
!! read angles, pattern centers, and deformation tensors from a text file

use mod_EMsoft
use mod_io
use mod_rotations
use mod_quaternions

IMPLICIT NONE

class(EBSDkin_T), INTENT(INOUT)            :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
integer(kind=irg),INTENT(OUT)           :: numangles
type(QuaternionArray_T), INTENT(INOUT)  :: qAR
type(EBSDkinAnglePCDefType),INTENT(INOUT)  :: orpcdef
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
    call Message%printError("EBSDkinreadorpcdef","Other orientation formats to be implemented; only Euler for now")
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

end subroutine EBSDkinreadorpcdef_

!--------------------------------------------------------------------------
recursive subroutine GenerateDetector_(self, verbose)
!DEC$ ATTRIBUTES DLLEXPORT :: GenerateDetector_
!! author: MDG
!! version: 1.0
!! date: 11/18/21
!!
!! generate the detector array (no energy arrays needed)

use mod_io
use mod_Lambert
use mod_math

IMPLICIT NONE

class(EBSDkin_T), INTENT(INOUT)         :: self
logical,INTENT(IN),OPTIONAL             :: verbose

type(IO_T)                              :: Message
type(Lambert_T)                         :: L

real(kind=sgl),allocatable              :: scin_x(:), scin_y(:), testarray(:,:)  ! scintillator coordinate arrays [microns]
real(kind=sgl)                          :: alp, ca, sa, cw, sw, rhos, sx
real(kind=sgl)                          :: L2, Ls, Lc, calpha     ! distances
real(kind=sgl),allocatable              :: z(:,:)
integer(kind=irg)                       :: nix, niy, binx, biny , i, j,  istat, k, ipx, ipy, nsx, nsy, elp  ! various parameters

associate( enl => self%nml, EBSDkindetector => self%det )

!====================================
! ------ generate the detector arrays
!====================================
! This needs to be done only once for a given detector geometry
allocate(scin_x(enl%numsx),scin_y(enl%numsy),stat=istat)

! change to the detector point of view necessitates negating the x pattern center coordinate
scin_x = - ( -enl%xpc - ( 1.0 - enl%numsx ) * 0.5 - (/ (i-1, i=1,enl%numsx) /) ) * enl%delta
scin_y = ( enl%ypc - ( 1.0 - enl%numsy ) * 0.5 - (/ (i-1, i=1,enl%numsy) /) ) * enl%delta

! auxiliary angle to rotate between reference frames
alp = 0.5 * cPi - (enl%sigma - enl%thetac) * dtor
ca = cos(alp)
sa = sin(alp)

! for now, don't take omega angle into account
cw = 1.0 ! cos(mcnl%omega * dtor)
sw = 0.0 ! sin(mcnl%omega * dtor)

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
   EBSDkindetector%rgx(j,elp-i) = (scin_y(i) * ca + sa * Ls) * rhos!Ls * rhos
   EBSDkindetector%rgy(j,elp-i) = Lc * rhos!(scin_x(i) * cw + Lc * sw) * rhos
   EBSDkindetector%rgz(j,elp-i) = (-sa * scin_y(i) + ca * Ls) * rhos!(-sw * scin_x(i) + Lc * cw) * rhos
  end do
end do
deallocate(scin_x, scin_y)

! normalize the direction cosines.
allocate(z(enl%numsx,enl%numsy))
  z = 1.0/sqrt(EBSDkindetector%rgx*EBSDkindetector%rgx+EBSDkindetector%rgy*EBSDkindetector%rgy & 
              +EBSDkindetector%rgz*EBSDkindetector%rgz)
  EBSDkindetector%rgx = EBSDkindetector%rgx*z
  EBSDkindetector%rgy = EBSDkindetector%rgy*z
  EBSDkindetector%rgz = EBSDkindetector%rgz*z
deallocate(z)
!====================================

if (present(verbose)) then
  if (verbose.eqv..TRUE.) call Message%printMessage(' -> completed detector generation', frm = "(A)")
end if 

end associate

end subroutine GenerateDetector_

!--------------------------------------------------------------------------
subroutine ComputeEBSDkinPatterns_(self, EMsoft, HDF, HDFnames, mem, numangles, angles, progname, nmldeffile)
!DEC$ ATTRIBUTES DLLEXPORT :: ComputeEBSDkinPatterns_
!! author: MDG
!! version: 1.0
!! date: 11/18/21
!!
!! compute an energy-weighted EBSDkin pattern

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
use mod_memory

IMPLICIT NONE

class(EBSDkin_T), INTENT(INOUT)         :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
type(HDF_T),INTENT(INOUT)               :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames
type(memory_T), INTENT(INOUT)           :: mem 
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
type(memory_T)                          :: memth

! all geometrical parameters and filenames
real(kind=dbl)                          :: prefactor, qz(3)

! allocatable arrays
real(kind=sgl),allocatable              :: EBSDkinpattern(:,:), binned(:,:)        ! array with EBSDkin patterns
real(kind=sgl),allocatable              :: z(:,:)               ! used to store the computed patterns before writing to disk
real(kind=sgl),allocatable              :: energywf(:), eulerangles(:,:)

! arrays for each OpenMP thread
real(kind=sgl),allocatable              :: tmLPNH(:,:) , tmLPSH(:,:)
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
integer(kind=irg),parameter             :: storemax = 20        ! number of EBSDkin patterns stored in one output block
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

associate( enl => self%nml, EBSDkindetector => self%det )

!====================================
! max number of OpenMP threads on this platform
maxthreads = omp_get_max_threads()
call setRotationPrecision('d')

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

!====================================
! init a bunch of parameters
!====================================
! binned pattern array
  binx = enl%numsx/enl%binning
  biny = enl%numsy/enl%binning
  bindx = 1.0/float(enl%binning)**2
!====================================

! get the crystal structure data
call cell%getCrystalData(self%xtalname, SG, EMsoft, useHDF=HDF)

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
datagroupname = trim(HDFnames%get_ProgramData()) ! 'EBSDkin' or 'TKD'
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
  dataset = SC_EMEBSDkinNML
end if
hdferr = HDF%writeDatasetTextFile(dataset, nmldeffile)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetTextFile EMEBSDkinNML/EMTKDNML', hdferr)

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

! create the EBSDkin group and add a HDF_FileVersion attribute to it
hdferr = HDF%createGroup(datagroupname)
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup EBSDkin/TKD', hdferr)
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
stringarray(1)= trim(self%xtalname)
hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetStringArray xtalname', hdferr)

dataset = SC_numangles
hdferr = HDF%writeDatasetInteger(dataset, numangles)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetInteger numangles', hdferr)

! and add the Euler angles to the output file
call mem%alloc(eulerangles, (/ 3,numangles /), 'eulerangles')
do i=1,numangles
  quat = angles%getQuatfromArray(i)
  call qq%q_setd( quat%get_quatd() )
  eu = qq%qe()
  eulerangles(1:3,i) = sngl( eu%e_copyd() )
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
  call mem%alloc(istart, (/ nthreads-1,totnumbatches /), 'istart', startdims = (/ 0, 1/) )
  call mem%alloc(istop, (/ nthreads-1,totnumbatches/), 'istop', startdims = (/ 0, 1 /) )
  call mem%alloc(patinbatch, (/ totnumbatches /), 'patinbatch')
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
  allocate(batchpatterns(binx,biny,ninbatch*nthreads))
end if
if (trim(bitmode).eq.'int') then
  call mem%alloc(batchpatternsint, (/binx,biny,ninbatch*nthreads /), 'batchpatternsint')
end if
if (trim(bitmode).eq.'float') then
  call mem%alloc(batchpatterns32, (/binx,biny,ninbatch*nthreads /), 'batchpatterns32')
end if
if (trim(bitmode).eq.'dict') then
  call mem%alloc(batchpatterns32lin, (/correctsize,ninbatch*nthreads /), 'batchpatterns32lin')
end if

!====================================
! here we also create a mask if necessary
  call mem%alloc(mask, (/ binx,biny /), 'mask', 1.0) 
  call mem%alloc(masklin, (/ binx*biny /), 'masklin', 1.0)
  if (enl%maskpattern.eq.'y') then
! create the circular mask in a potentially rectangular array
    maskradius = (minval( (/ binx, biny /) ) / 2 )**2
    call mem%alloc(lx, (/ binx /), 'lx')
    call mem%alloc(ly, (/ biny /), 'ly')
    lx = (/ (float(i),i=1,binx) /) - float(binx/2)
    ly = (/ (float(i),i=1,biny) /) - float(biny/2)
    do i=1,binx
      do j=1,biny
        if ((lx(i)**2+ly(j)**2).gt.maskradius) mask(i,j) = 0.0
      end do
    end do
    call mem%dealloc(lx, 'lx')
    call mem%dealloc(ly, 'ly')
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
scl = dble(self%npx)

!====================================
! define the integer parameter list for the CalcEBSDkinPatternSingleFull call
ipar(1) = enl%binning
ipar(2) = enl%numsx
ipar(3) = enl%numsy
ipar(4) = self%npx
ipar(5) = self%npx
! ipar(6) = EBSDkinMCdata%numEbins
! ipar(7) = EBSDkinMCdata%numEbins

!====================================
! set the number of OpenMP threads
call OMP_setNThreads(nthreads)

memth = memory_T( nt = nthreads )

call timer%Time_tick()

write (*,*) ' prefactor : ', prefactor 
prefactor = 1.D0

!====================================
!====================================
do ibatch=1,totnumbatches

! use OpenMP to run on multiple cores ...
!$OMP PARALLEL default(shared)  PRIVATE(TID,iang,i,j,istat,EBSDkinpattern,binned,idum,bpat,ma,mi,threadbatchpatterns,bpatint)&
!$OMP& PRIVATE(tmLPNH, tmLPSH, trgx, trgy, trgz, dims2, dims3, threadbatchpatternsint, threadbatchpatterns32)&
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
  call memth%alloc(trgx, shape(self%det%rgx), 'trgx', TID=TID)
  call memth%alloc(trgy, shape(self%det%rgx), 'trgy', TID=TID)
  call memth%alloc(trgz, shape(self%det%rgx), 'trgz', TID=TID)
  call memth%alloc(tmLPNH, shape(self%mLPNH), 'tmLPNH', TID=TID)
  call memth%alloc(tmLPSH, shape(self%mLPNH), 'tmLPSH', TID=TID)

! and copy the data in
  trgx = self%det%rgx
  trgy = self%det%rgy
  trgz = self%det%rgz
  tmLPNH = self%mLPNH
  tmLPSH = self%mLPSH

! allocate the arrays that will hold the computed pattern
  call memth%alloc(binned, (/ binx,biny /), 'binned', TID=TID) 
  if (trim(bitmode).eq.'char') then
    allocate(bpat(binx,biny))
  end if
  if (trim(bitmode).eq.'int') then
    call memth%alloc(bpatint, (/ binx,biny /), 'bpatint', TID=TID)
  end if
  if (trim(bitmode).eq.'dict') then
    call memth%alloc(bpatint, (/ binx,biny /), 'bpatint', TID=TID) 
    call memth%alloc(binnedvec, (/ correctsize /), 'binnedvec', TID=TID)
  end if

! this array requires some care in terms of its size parameters...
  if ((singlebatch.eqv..TRUE.).AND.(ibatch.eq.totnumbatches)) then
     if (TID.eq.nthreads-1) then
      if (trim(bitmode).eq.'char') then
        allocate(threadbatchpatterns(binx,biny,nlastremainder))
      end if
      if (trim(bitmode).eq.'int') then
        call memth%alloc(threadbatchpatternsint, (/ binx,biny,nlastremainder /),'threadbatchpatternsint', TID=TID) 
      end if
      if (trim(bitmode).eq.'float') then
        call memth%alloc(threadbatchpatterns32, (/ binx,biny,nlastremainder /),'threadbatchpatterns32', TID=TID) 
      end if
      if (trim(bitmode).eq.'dict') then
        call memth%alloc(threadbatchpatterns32lin, (/ correctsize,nlastremainder /),'threadbatchpatterns32lin', TID=TID) 
      end if
    else
      if (trim(bitmode).eq.'char') then
        allocate(threadbatchpatterns(binx,biny,ninlastbatch)) 
      end if
      if (trim(bitmode).eq.'int') then
        call memth%alloc(threadbatchpatternsint, (/ binx,biny,ninlastbatch /), 'threadbatchpatternsint', TID=TID) 
      end if
      if (trim(bitmode).eq.'float') then
        call memth%alloc(threadbatchpatterns32, (/ binx,biny,ninlastbatch /),'threadbatchpatterns32', TID=TID) 
      end if
      if (trim(bitmode).eq.'dict') then
        call memth%alloc(threadbatchpatterns32lin, (/ correctsize,ninlastbatch /),'threadbatchpatterns32lin', TID=TID) 
      end if
    end if
  else
    if (trim(bitmode).eq.'char') then
      allocate(threadbatchpatterns(binx,biny,ninbatch)) 
    end if
    if (trim(bitmode).eq.'int') then
      call memth%alloc(threadbatchpatternsint, (/ binx,biny,ninbatch /),'threadbatchpatternsint', TID=TID) 
    end if
    if (trim(bitmode).eq.'float') then
      call memth%alloc(threadbatchpatterns32, (/ binx,biny,ninbatch /),'threadbatchpatterns32', TID=TID) 
    end if
    if (trim(bitmode).eq.'dict') then
      call memth%alloc(threadbatchpatterns32lin, (/ correctsize,ninbatch /),'threadbatchpatterns32lin', TID=TID) 
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
! then convert these individually to the correct EBSDkin pattern location
    binned = 0.0

    if (includeFmatrix.eqv..TRUE.) then
      call self%CalcEBSDkinPatternSingleFull(ipar,angles%getQuatfromArray(iang),tmLPNH,tmLPSH,trgx,trgy,trgz,binned, &
                                            mask,prefactor,Fmatrix_inverse)
    else
      call self%CalcEBSDkinPatternSingleFull(ipar,angles%getQuatfromArray(iang),tmLPNH,tmLPSH,trgx,trgy,trgz,binned, &
                                            mask,prefactor)
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
      write (*,*) 'float range : ',minval(binned), maxval(binned), minval(tmLPNH), maxval(tmLPNH)
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

if (trim(bitmode).eq.'char') then
  deallocate(threadbatchpatterns)
  deallocate(bpat)
end if
if (trim(bitmode).eq.'int') then
  call memth%dealloc(threadbatchpatternsint, 'threadbatchpatternsint', TID=TID) 
  call memth%dealloc(bpatint, 'bpatint', TID=TID)
end if
if (trim(bitmode).eq.'float') then
  call memth%dealloc(threadbatchpatterns32, 'threadbatchpatterns32', TID=TID) 
end if
if (trim(bitmode).eq.'dict') then
  call memth%dealloc(threadbatchpatterns32lin, 'threadbatchpatterns32lin', TID=TID) 
  call memth%dealloc(bpatint, 'bpatint', TID=TID) 
  call memth%dealloc(binnedvec, 'binnedvec', TID=TID)
end if
call memth%dealloc(trgx, 'trgx', TID=TID)
call memth%dealloc(trgy, 'trgy', TID=TID)
call memth%dealloc(trgz, 'trgz', TID=TID)
call memth%dealloc(tmLPNH, 'tmLPNH', TID=TID)
call memth%dealloc(tmLPSH, 'tmLPSH', TID=TID)
call memth%dealloc(binned, 'binned', TID=TID)

!$OMP END PARALLEL

! test for memory allocations in the threaded region 
! call memth%thread_memory_use()

! here we write all the entries in the batchpatterns array to the HDF file as a hyperslab
! =====================================================
! The following write commands constitute HDF_FileVersion = 4.0
! =====================================================
if (isTKD.eqv..TRUE.) then
  dataset = SC_TKDpatterns
else
  dataset = SC_EBSDkinpatterns
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
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabCharArray EBSDkinpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'int') then
       hdferr = HDF%writeHyperslabIntegerArray(dataset, batchpatternsint(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabIntegerArray EBSDkinpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'float') then
       hdferr = HDF%writeHyperslabFloatArray(dataset, batchpatterns32(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabFloatArray EBSDkinpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'dict') then
       hdferr = HDF%writeHyperslabFloatArray(dataset, batchpatterns32lin(1:correctsize,1:dim1), hdims2, offset2, &
                                              dims2)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabFloatArray EBSDkinpatterns', hdferr)
     end if
   else
     if (trim(bitmode).eq.'char') then
       hdferr = HDF%writeHyperslabCharArray(dataset, batchpatterns(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3, insert)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabCharArray EBSDkinpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'int') then
       hdferr = HDF%writeHyperslabIntegerArray(dataset, batchpatternsint(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3, insert)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabIntegerArray EBSDkinpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'float') then
       hdferr = HDF%writeHyperslabFloatArray(dataset, batchpatterns32(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3, insert)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabFloatArray EBSDkinpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'dict') then
       hdferr = HDF%writeHyperslabFloatArray(dataset, batchpatterns32lin(1:correctsize,1:dim1), hdims2, offset2, &
                                              dims2, insert)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabFloatArray EBSDkinpatterns', hdferr)
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
if (hdferr.ne.0) call HDF%error_check('HDF_openGroup EBSDkin/TKD', hdferr)

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

call mem%dealloc(EBSDkindetector%rgx, 'EBSDkindetector%rgx')
call mem%dealloc(EBSDkindetector%rgy, 'EBSDkindetector%rgy')
call mem%dealloc(EBSDkindetector%rgz, 'EBSDkindetector%rgz')
call mem%dealloc(stringarray, 'stringarray')
call mem%dealloc(eulerangles, 'eulerangles')
call mem%dealloc(istart, 'istart')
call mem%dealloc(istop, 'istop')
call mem%dealloc(patinbatch, 'patinbatch')
call mem%dealloc(batchpatterns32lin, 'batchpatterns32lin')
call mem%dealloc(mask, 'mask')
call mem%dealloc(masklin, 'masklin')

end associate

! test for memory deallocations 
! call mem%allocated_memory_use()

end subroutine ComputeEBSDkinPatterns_

!--------------------------------------------------------------------------
recursive subroutine CalcEBSDkinPatternSingleFull_(self,ipar,qq,mLPNH,mLPSH,rgx,rgy,rgz,binned,mask, &
                                                  prefactor, Fmatrix)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcEBSDkinPatternSingleFull_
 !! author: MDG
 !! version: 1.0
 !! date: 01/06/20
 !!
 !! compute a single EBSDkin pattern 

use mod_Lambert
use mod_quaternions
use mod_rotations

IMPLICIT NONE

integer, parameter                              :: K4B=selected_int_kind(9)

class(EBSDkin_T), INTENT(INOUT)                 :: self
integer(kind=irg),INTENT(IN)                    :: ipar(7)
type(Quaternion_T),INTENT(IN)                   :: qq
real(kind=dbl),INTENT(IN)                       :: prefactor
real(kind=sgl),INTENT(IN)                       :: mLPNH(-ipar(4):ipar(4),-ipar(5):ipar(5))
real(kind=sgl),INTENT(IN)                       :: mLPSH(-ipar(4):ipar(4),-ipar(5):ipar(5))
real(kind=sgl),INTENT(IN)                       :: rgx(ipar(2),ipar(3))
real(kind=sgl),INTENT(IN)                       :: rgy(ipar(2),ipar(3))
real(kind=sgl),INTENT(IN)                       :: rgz(ipar(2),ipar(3))
real(kind=sgl),INTENT(OUT)                      :: binned(ipar(2)/ipar(1),ipar(3)/ipar(1))
real(kind=sgl),INTENT(IN)                       :: mask(ipar(2)/ipar(1),ipar(3)/ipar(1))
real(kind=dbl),INTENT(IN),optional              :: Fmatrix(3,3)

real(kind=sgl),allocatable                      :: EBSDkinpattern(:,:)
real(kind=sgl),allocatable                      :: wf(:)
real(kind=sgl)                                  :: dc(3),ixy(2),scl,bindx, tmp
real(kind=sgl)                                  :: dx,dy,dxm,dym, x, y, z
integer(kind=irg)                               :: ii,jj,kk,istat
integer(kind=irg)                               :: nix,niy,nixp,niyp


! ipar(1) = ebsdnl%binning
! ipar(2) = ebsdnl%numsx
! ipar(3) = ebsdnl%numsy
! ipar(4) = ebsdnl%npx
! ipar(5) = ebsdnl%npy

write (*,*) 'ipar = ', ipar, prefactor  

allocate(EBSDkinpattern(ipar(2),ipar(3)),stat=istat)

binned = 0.0
EBSDkinpattern = 0.0

scl = float(ipar(4))

do ii = 1,ipar(2)
    do jj = 1,ipar(3)
      if ((ii.eq.1).and.(jj.eq.1)) write (*,*) 'inside ',minval(mLPNH),maxval(mLPNH)
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
        if (dc(3) .ge. 0.0) then
              EBSDkinpattern(ii,jj) = EBSDkinpattern(ii,jj) + ( mLPNH(nix,niy) * dxm * dym + &
                                             mLPNH(nixp,niy) * dx * dym + mLPNH(nix,niyp) * dxm * dy + &
                                             mLPNH(nixp,niyp) * dx * dy )

        else
              EBSDkinpattern(ii,jj) = EBSDkinpattern(ii,jj) + ( mLPSH(nix,niy) * dxm * dym + &
                                             mLPSH(nixp,niy) * dx * dym + mLPSH(nix,niyp) * dxm * dy + &
                                             mLPSH(nixp,niyp) * dx * dy )

        end if
    end do
end do

EBSDkinpattern = prefactor * EBSDkinpattern

! do we need to bin the pattern ?

! 17/07/2020 Clément Lafond, temporary fix to avoid NaN value, and negatives
! values when EBSDkin pattern size is large
do ii=1,ipar(2)
    do jj=1,ipar(3)
        if (isnan(EBSDkinpattern(ii,jj)).or.EBSDkinpattern(ii,jj).lt.0.0) then
          EBSDkinpattern(ii,jj) = 0.0
        end if
    end do
end do

if (ipar(1) .ne. 1) then
    do ii=1,ipar(2),ipar(1)
        do jj=1,ipar(3),ipar(1)
            binned(ii/ipar(1)+1,jj/ipar(1)+1) = &
            sum(EBSDkinpattern(ii:ii+ipar(1)-1,jj:jj+ipar(1)-1))
        end do
    end do
! and divide by binning^2
!   binned = binned * bindx
else
    binned = EBSDkinpattern
end if

binned = binned * mask

deallocate(EBSDkinpattern)

end subroutine CalcEBSDkinPatternSingleFull_

!--------------------------------------------------------------------------
subroutine ComputedeformedEBSDkinPatterns_(self, EMsoft,  HDF, HDFnames, &
                                          numangles, angles, orpcdef, progname, nmldeffile)
!DEC$ ATTRIBUTES DLLEXPORT :: ComputedeformedEBSDkinPatterns_
  !! author: MDG
  !! version: 1.0
  !! date: 02/18/20
  !!
  !! compute a energy-weighted EBSDkin patterns with different pattern centers and deformation states

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

class(EBSDkin_T), INTENT(INOUT)         :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
type(HDF_T),INTENT(INOUT)               :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames
integer(kind=irg),INTENT(IN)            :: numangles
type(QuaternionArray_T), INTENT(IN)     :: angles
type(EBSDkinAnglePCDefType),INTENT(IN)     :: orpcdef
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
real(kind=sgl),allocatable              :: EBSDkinpattern(:,:), binned(:,:)        ! array with EBSDkin patterns
real(kind=sgl),allocatable              :: z(:,:)               ! used to store the computed patterns before writing to disk
real(kind=sgl),allocatable              :: energywf(:), eulerangles(:,:)

! arrays for each OpenMP thread
real(kind=sgl),allocatable              :: tmLPNH(:,:) , tmLPSH(:,:)
real(kind=sgl),allocatable              :: trgx(:,:), trgy(:,:), trgz(:,:)          ! auxiliary detector arrays needed for interpolation

! various items
integer(kind=irg)                       :: i, j, iang, jang, k, io_int(6), hdferr, dim2          ! various counters
integer(kind=irg)                       :: istat, ipar(7), tick, tock, tickstart
integer(kind=irg)                       :: nix, niy, binx, biny, nixp, niyp, maxthreads,nextra,ninlastbatch,nlastremainder     ! various parameters
integer(kind=irg)                       :: NUMTHREADS, TID   ! number of allocated threads, thread ID
integer(kind=irg)                       :: ninbatch, nbatches,nremainder,ibatch,nthreads,maskradius,nlastbatches, totnumbatches
integer(kind=irg),allocatable           :: istart(:,:), istop(:,:), patinbatch(:)

real(kind=sgl)                          :: bindx, ma, mi, tstart, tstop, io_real(3)
real(kind=dbl),parameter                :: nAmpere = 6.241D+18   ! Coulomb per second
integer(kind=irg),parameter             :: storemax = 20        ! number of EBSDkin patterns stored in one output block
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

associate( enl => self%nml, EBSDkindetector => self%det )

!====================================
! max number of OpenMP threads on this platform
maxthreads = omp_get_max_threads()

!====================================
! what is the output format?  GUI or BIN ?
outputformat = enl%outputformat

!====================================
! bit depth and format of output
call get_bit_parameters(enl%bitdepth, numbits, bitrange, bitmode)

!====================================
! init a bunch of parameters
!====================================
! binned pattern array
  binx = enl%numsx/enl%binning
  biny = enl%numsy/enl%binning
  bindx = 1.0/float(enl%binning)**2

!====================================

! get the crystal structure data
call cell%getCrystalData(self%xtalname, SG, EMsoft, useHDF=HDF)

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
datagroupname = trim(HDFnames%get_ProgramData()) ! 'EBSDkin' or 'TKD'
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
  dataset = SC_EMEBSDkinNML
end if
hdferr = HDF%writeDatasetTextFile(dataset, nmldeffile)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetTextFile EMEBSDkinNML/EMTKDNML', hdferr)

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

! create the EBSDkin group and add a HDF_FileVersion attribute to it
hdferr = HDF%createGroup(datagroupname)
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup EBSDkin/TKD', hdferr)
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
stringarray(1)= trim(self%xtalname)
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
    if (nlastremainder.le.0) then 
      ninlastbatch = numangles/nthreads
      nlastremainder = numangles - (nthreads-1)*ninlastbatch
    end if
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
scl = dble(self%npx)

!====================================
! define the integer parameter list for the CalcEBSDkinPatternSingleFull call
ipar(1) = enl%binning
ipar(2) = enl%numsx
ipar(3) = enl%numsy
ipar(4) = self%npx
ipar(5) = self%npx
! ipar(6) = EBSDkinMCdata%numEbins
! ipar(7) = EBSDkinMCdata%numEbins

!====================================
! set the number of OpenMP threads
call OMP_setNThreads(nthreads)

call timer%Time_tick()

!====================================
!====================================
do ibatch=1,totnumbatches

! use OpenMP to run on multiple cores ...
!$OMP PARALLEL default(shared)  PRIVATE(TID,iang,i,j,istat,EBSDkinpattern,binned,idum,bpat,ma,mi,threadbatchpatterns,bpatint)&
!$OMP& PRIVATE(tmLPNH, tmLPSH, trgx, trgy, trgz, threadbatchpatternsint, threadbatchpatterns32, prefactor)&
!$OMP& PRIVATE(Fmatrix_inverse, nel, Fmatrix)

  NUMTHREADS = OMP_GET_NUM_THREADS()
  TID = OMP_GET_THREAD_NUM()

! each thread needs a private copy of the master and accum arrays; not having
! those can produce poor scaling... in addition, they need to be recomputed for each pattern !
  allocate(trgx(enl%numsx,enl%numsy), trgy(enl%numsx,enl%numsy), trgz(enl%numsx,enl%numsy))
  allocate(tmLPNH(enl%numsx,enl%numsy), tmLPSH(enl%numsx,enl%numsy))
! and copy the data in
  tmLPNH = self%mLPNH
  tmLPSH = self%mLPSH

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
    call self%GeneratemyEBSDkinDetector(enl%numsx, enl%numsy, trgx, trgy, trgz, &
                                        orpcdef%pcs(1:3,iang),bg=.TRUE.)
! we pick a reasonable value here ...
    prefactor = 1.D0 
    binned = 0.0

    if (includeFmatrix.eqv..TRUE.) then
      call self%CalcEBSDkinPatternSingleFull(ipar,angles%getQuatfromArray(iang),tmLPNH,tmLPSH,trgx,trgy,trgz,binned, &
                                     mask,prefactor,Fmatrix_inverse)
    else
      call self%CalcEBSDkinPatternSingleFull(ipar,angles%getQuatfromArray(iang),tmLPNH,tmLPSH,trgx,trgy,trgz,binned, &
                                     mask,prefactor)
    end if

! from here on everything is the same as before...

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
  dataset = SC_EBSDkinpatterns
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
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabCharArray EBSDkinpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'int') then
       hdferr = HDF%writeHyperslabIntegerArray(dataset, batchpatternsint(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabIntegerArray EBSDkinpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'float') then
       hdferr = HDF%writeHyperslabFloatArray(dataset, batchpatterns32(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabFloatArray EBSDkinpatterns', hdferr)
     end if
   else
     if (trim(bitmode).eq.'char') then
       hdferr = HDF%writeHyperslabCharArray(dataset, batchpatterns(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3, insert)
       if (hdferr.ne.0) call  HDF%error_check('HDF_writeHyperslabCharArray EBSDkinpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'int') then
       hdferr = HDF%writeHyperslabIntegerArray(dataset, batchpatternsint(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3, insert)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabIntegerArray EBSDkinpatterns', hdferr)
     end if
     if (trim(bitmode).eq.'float') then
       hdferr = HDF%writeHyperslabFloatArray(dataset, batchpatterns32(1:binx,1:biny,1:dim2), hdims, offset, &
                                              dims3, insert)
       if (hdferr.ne.0) call HDF%error_check('HDF_writeHyperslabFloatArray EBSDkinpatterns', hdferr)
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
  datagroupname = "EBSDkin"
end if
hdferr = HDF%openGroup(datagroupname)
if (hdferr.ne.0) call HDF%error_check('HDF_openGroup EBSDkin/TKD', hdferr)

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

end subroutine ComputedeformedEBSDkinPatterns_

!--------------------------------------------------------------------------
recursive subroutine GeneratemyEBSDkinDetector_(self, nsx, nsy, tgx, tgy, tgz, patcntr, bg)
!DEC$ ATTRIBUTES DLLEXPORT :: GeneratemyEBSDkinDetector_
!! author: MDG
!! version: 1.0
!! date: 02/05/20
!!
!! generate the detector arrays for the case where each pattern has a (slightly) different detector configuration

use mod_io
use mod_Lambert
use mod_MCfiles

IMPLICIT NONE

class(EBSDkin_T), INTENT(INOUT)         :: self
integer(kind=irg),INTENT(IN)            :: nsx
integer(kind=irg),INTENT(IN)            :: nsy
real(kind=sgl),INTENT(INOUT)            :: tgx(nsx,nsy)
!f2py intent(in,out) ::  tgx
real(kind=sgl),INTENT(INOUT)            :: tgy(nsx,nsy)
!f2py intent(in,out) ::  tgy
real(kind=sgl),INTENT(INOUT)            :: tgz(nsx,nsy)
!f2py intent(in,out) ::  tgz
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

associate( enl => self%nml )

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
alp = 0.5 * cPi - (enl%sigma - enl%thetac) * dtor
ca = cos(alp)
sa = sin(alp)

cw = cos(enl%omega * dtor)
sw = sin(enl%omega * dtor)

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


end associate

end subroutine GeneratemyEBSDkinDetector_

end module mod_EBSDkin
