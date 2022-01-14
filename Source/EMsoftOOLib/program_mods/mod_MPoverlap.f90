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

module mod_MPoverlap
  !! author: MDG
  !! version: 1.0
  !! date: 03/22/20
  !!
  !! class definition for the EMMPoverlap program

use mod_kinds
use mod_global

IMPLICIT NONE

! namelist for the EMMPoverlap program
type, public :: MPoverlapNameListType
  integer(kind=irg)       :: newpgnum
  integer(kind=irg)       :: PatternAxisA(3)
  integer(kind=irg)       :: HorizontalAxisA(3)
  real(kind=sgl)          :: tA(3)
  real(kind=sgl)          :: tB(3)
  real(kind=sgl)          :: tA2(3)
  real(kind=sgl)          :: tC(3)
  real(kind=sgl)          :: tA3(3)
  real(kind=sgl)          :: tD(3)
  real(kind=sgl)          :: gA(3)
  real(kind=sgl)          :: gB(3)
  real(kind=sgl)          :: gA2(3)
  real(kind=sgl)          :: gC(3)
  real(kind=sgl)          :: gA3(3)
  real(kind=sgl)          :: gD(3)
  real(kind=sgl)          :: fracB
  real(kind=sgl)          :: fracC
  real(kind=sgl)          :: fracD
  character(fnlen)        :: masterfileA
  character(fnlen)        :: masterfileB
  character(fnlen)        :: masterfileC
  character(fnlen)        :: masterfileD
  character(fnlen)        :: datafile
  character(fnlen)        :: h5copypath
  character(fnlen)        :: overlapmode
  character(fnlen)        :: modality
end type MPoverlapNameListType

! class definition
type, public :: MPoverlap_T
private
  character(fnlen)              :: nmldeffile = 'EMMPoverlap.nml'
  type(MPoverlapNameListType)   :: nml

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: MPoverlap_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: MPoverlap => MPoverlap_

end type MPoverlap_T

! the constructor routine for this class
interface MPoverlap_T
  module procedure MPoverlap_constructor
end interface MPoverlap_T

contains

!--------------------------------------------------------------------------
type(MPoverlap_T) function MPoverlap_constructor( nmlfile ) result(MPoverlap)
!DEC$ ATTRIBUTES DLLEXPORT :: MPoverlap_constructor
!! author: MDG
!! version: 1.0
!! date: 03/22/20
!!
!! constructor for the MPoverlap_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call MPoverlap%readNameList(nmlfile)

end function MPoverlap_constructor

!--------------------------------------------------------------------------
subroutine MPoverlap_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: MPoverlap_destructor
!! author: MDG
!! version: 1.0
!! date: 03/22/20
!!
!! destructor for the MPoverlap_T Class

IMPLICIT NONE

type(MPoverlap_T), INTENT(INOUT)  :: self

call reportDestructor('MPoverlap_T')

end subroutine MPoverlap_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 03/22/20
!!
!! read the namelist from an nml file for the MPoverlap_T Class

use mod_io
use mod_EMsoft

IMPLICIT NONE

class(MPoverlap_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft
type(IO_T)                           :: Message
logical                              :: skipread = .FALSE.

integer(kind=irg)       :: newpgnum
integer(kind=irg)       :: PatternAxisA(3)
integer(kind=irg)       :: HorizontalAxisA(3)
real(kind=sgl)          :: tA(3)
real(kind=sgl)          :: tB(3)
real(kind=sgl)          :: tA2(3)
real(kind=sgl)          :: tC(3)
real(kind=sgl)          :: tA3(3)
real(kind=sgl)          :: tD(3)
real(kind=sgl)          :: gA(3)
real(kind=sgl)          :: gB(3)
real(kind=sgl)          :: gA2(3)
real(kind=sgl)          :: gC(3)
real(kind=sgl)          :: gA3(3)
real(kind=sgl)          :: gD(3)
real(kind=sgl)          :: fracB
real(kind=sgl)          :: fracC
real(kind=sgl)          :: fracD
character(fnlen)        :: masterfileA
character(fnlen)        :: masterfileB
character(fnlen)        :: masterfileC
character(fnlen)        :: masterfileD
character(fnlen)        :: h5copypath
character(fnlen)        :: overlapmode
character(fnlen)        :: datafile
character(fnlen)        :: modality

! define the IO namelist to facilitate passing variables to the program.
namelist / MPoverlapdata / PatternAxisA, tA, tB, gA, gB, masterfileA, masterfileB, modality, &
                           datafile, HorizontalAxisA, overlapmode, newpgnum, tC, gC, tD, gD, fracB, &
                           fracC, fracD, masterfileC, masterfileD, gA2, gA3, tA2, tA3, h5copypath

! set the input parameters to default values (except for xtalname, which must be present)
newpgnum        = -1                            ! -1 means 'use the point group of phase A'
PatternAxisA    = (/ 0, 0, 1 /)                 ! center axis for output pattern
HorizontalAxisA = (/ 1, 0, 0 /)                 ! horizontal axis for output pattern
tA              = (/0.0, 0.0, 1.0/)             ! direction vector in crystal A
tB              = (/0.0, 0.0, 1.0/)             ! direction vector in crystal B
tA2             = (/0.0, 0.0, 1.0/)             ! direction vector in crystal A
tC              = (/0.0, 0.0, 1.0/)             ! direction vector in crystal C
tA3             = (/0.0, 0.0, 1.0/)             ! direction vector in crystal A
tD              = (/0.0, 0.0, 1.0/)             ! direction vector in crystal D
gA              = (/1.0, 0.0, 0.0/)             ! plane normal in crystal A
gB              = (/1.0, 0.0, 0.0/)             ! plane normal in crystal B
gA2             = (/1.0, 0.0, 0.0/)             ! plane normal in crystal A
gC              = (/1.0, 0.0, 0.0/)             ! plane normal in crystal C
gA3             = (/1.0, 0.0, 0.0/)             ! plane normal in crystal A
gD              = (/1.0, 0.0, 0.0/)             ! plane normal in crystal D
fracB           = 0.25                          ! volume fraction of phase B
fracC           = 0.25                          ! volume fraction of phase C
fracD           = 0.25                          ! volume fraction of phase D
masterfileA     = 'undefined'   ! filename
masterfileB     = 'undefined'   ! filename
masterfileC     = 'undefined'   ! filename
masterfileD     = 'undefined'   ! filename
h5copypath      = 'undefined'   ! filename
datafile        = 'undefined'   ! output file name
overlapmode     = 'series'      ! options are 'full' or 'series'
modality        = 'EBSD'        ! 'EBSD, 'ECP', or 'TKD'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=MPoverlapdata)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(masterfileA).eq.'undefined') then
  call Message%printError('readNameList:',' master pattern file name A is undefined in '//nmlfile)
 end if

 if (trim(masterfileB).eq.'undefined') then
  call Message%printError('readNameList:',' master pattern file name B is undefined in '//nmlfile)
 end if

 if (trim(datafile).eq.'undefined') then
  call Message%printError('readNameList:',' output file name is undefined in '//nmlfile)
 end if
end if

! if we get here, then all appears to be ok, and we need to fill in the emnl fields
self%nml%newpgnum = newpgnum
self%nml%PatternAxisA = PatternAxisA
self%nml%HorizontalAxisA = HorizontalAxisA
self%nml%tA = tA
self%nml%tB = tB
self%nml%tA2= tA2
self%nml%tC = tC
self%nml%tA3= tA3
self%nml%tD = tD
self%nml%gA = gA
self%nml%gB = gB
self%nml%gA2= gA2
self%nml%gC = gC
self%nml%gA3= gA3
self%nml%gD = gD
self%nml%fracB = fracB
self%nml%fracC = fracC
self%nml%fracD = fracD
self%nml%masterfileA = masterfileA
self%nml%masterfileB = masterfileB
self%nml%masterfileC = masterfileC
self%nml%masterfileD = masterfileD
self%nml%h5copypath = h5copypath
self%nml%datafile = datafile
self%nml%overlapmode = overlapmode
self%nml%modality = modality

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 03/22/20
!!
!! pass the namelist for the MPoverlap_T Class to the calling program

IMPLICIT NONE

class(MPoverlap_T), INTENT(INOUT)          :: self
type(MPoverlapNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG
!! version: 1.0
!! date: 03/22/20
!!
!! write namelist to HDF file

use mod_HDFsupport
use HDF5
use mod_HDFnames
use stringconstants

use ISO_C_BINDING

IMPLICIT NONE

class(MPoverlap_T), INTENT(INOUT)       :: self
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T),INTENT(INOUT)          :: HDFnames

integer(kind=irg),parameter             :: n_int = 1, n_real = 1
integer(kind=irg)                       :: hdferr, io_int(n_int)
character(20)                           :: intlist(n_int)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)
logical                                 :: g_exists, overwrite = .TRUE.

associate( emnl => self%nml )

! create the group for this namelist
groupname = trim(HDFnames%get_NMLlist()) !  SC_MPoverlapNameList
hdferr = HDF%createGroup(groupname)

! write all the single integers
io_int = (/ emnl%newpgnum /)
intlist(1) = 'newpgnum'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write a single real
dataset = 'fracB'
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetFloat(dataset, emnl%fracB, overwrite)
else
  hdferr = HDF%writeDatasetFloat(dataset, emnl%fracB)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create fracB dataset',hdferr)

! write a single real
dataset = 'fracC'
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetFloat(dataset, emnl%fracC, overwrite)
else
  hdferr = HDF%writeDatasetFloat(dataset, emnl%fracC)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create fracC dataset',hdferr)

! write a single real
dataset = 'fracD'
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetFloat(dataset, emnl%fracD, overwrite)
else
  hdferr = HDF%writeDatasetFloat(dataset, emnl%fracD)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create fracD dataset',hdferr)


! vectors
dataset = 'tA'
hdferr = HDF%writeDatasetFloatArray(dataset, emnl%tA, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create tA dataset',hdferr)

dataset = 'tA2'
hdferr = HDF%writeDatasetFloatArray(dataset, emnl%tA2, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create tA2 dataset',hdferr)

dataset = 'tA3'
hdferr = HDF%writeDatasetFloatArray(dataset, emnl%tA3, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create tA3 dataset',hdferr)

dataset = 'tB'
hdferr = HDF%writeDatasetFloatArray(dataset, emnl%tB, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create tB dataset',hdferr)

dataset = 'tC'
hdferr = HDF%writeDatasetFloatArray(dataset, emnl%tC, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create tC dataset',hdferr)

dataset = 'tD'
hdferr = HDF%writeDatasetFloatArray(dataset, emnl%tD, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create tD dataset',hdferr)


dataset = 'gA'
hdferr = HDF%writeDatasetFloatArray(dataset, emnl%gA, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create gA dataset',hdferr)

dataset = 'gA2'
hdferr = HDF%writeDatasetFloatArray(dataset, emnl%gA2, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create gA2 dataset',hdferr)

dataset = 'gA3'
hdferr = HDF%writeDatasetFloatArray(dataset, emnl%gA3, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create gA3 dataset',hdferr)

dataset = 'gB'
hdferr = HDF%writeDatasetFloatArray(dataset, emnl%gB, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create gB dataset',hdferr)

dataset = 'gC'
hdferr = HDF%writeDatasetFloatArray(dataset, emnl%gC, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create gC dataset',hdferr)

dataset = 'gD'
hdferr = HDF%writeDatasetFloatArray(dataset, emnl%gD, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create gD dataset',hdferr)

dataset = 'PatternAxisA'
hdferr = HDF%writeDatasetIntegerArray(dataset, emnl%PatternAxisA, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create PatternAxisA dataset',hdferr)

dataset = 'HorizontalAxisA'
hdferr = HDF%writeDatasetIntegerArray(dataset, emnl%HorizontalAxisA, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create HorizontalAxisA dataset',hdferr)

dataset = 'masterfileA'
line2(1) = emnl%masterfileA
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create masterfileA dataset',hdferr)

dataset = 'masterfileB'
line2(1) = emnl%masterfileB
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create masterfileB dataset',hdferr)

dataset = 'masterfileC'
line2(1) = emnl%masterfileC
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create masterfileC dataset',hdferr)

dataset = 'masterfileD'
line2(1) = emnl%masterfileD
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create masterfileD dataset',hdferr)

dataset = 'datafile'
line2(1) = emnl%datafile
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create datafile dataset',hdferr)


dataset = 'overlapmode'
line2(1) = emnl%overlapmode
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create overlapmode dataset',hdferr)


dataset = 'modality'
line2(1) = emnl%modality
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create overlapmode dataset',hdferr)

! and pop this group off the stack
call HDF%pop()


end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine MPoverlap_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: MPoverlap_
!! author: MDG
!! version: 1.0
!! date: 03/22/20
!!
!! perform the computations

use mod_EMsoft
use mod_HDFnames
use mod_crystallography
use mod_initializers
use mod_io
use mod_gvectors
use mod_symmetry
use mod_misc
use mod_rotations
use mod_quaternions
use mod_Lambert
use mod_HDFsupport
use HDF5
use mod_HDFnames
use mod_MCfiles
use mod_MPfiles
use mod_diffraction
use stringconstants

IMPLICIT NONE

class(MPoverlap_T), INTENT(INOUT)      :: self
type(EMsoft_T), INTENT(INOUT)          :: EMsoft
character(fnlen), INTENT(INOUT)        :: progname

type(MCOpenCLNameListType)             :: mcnlA, mcnlB, mcnlC, mcnlD
type(EBSDmasterNameListType)           :: mpnlEBSD
type(ECPmasterNameListType)            :: mpnlECP
type(TKDmasterNameListType)            :: mpnlTKD

type(IO_T)                             :: Message
type(HDF_T)                            :: HDF
type(HDFnames_T)                       :: HDFnames
type(cell_T)                           :: cellA, cellB, cellC, cellD
type(SpaceGroup_T)                     :: SGA, SGB, SGC, SGD
type(Lambert_T)                        :: L
type(Diffraction_T)                    :: Diff, DiffA, DiffB, DiffC, DiffD
type(MCfile_T)                         :: MCFTA, MCFTB, MCFTC, MCFTD
type(MPfile_T)                         :: MPFTA, MPFTB, MPFTC, MPFTD
type(OrientationRelation)              :: orelAB, orelAC, orelAD

integer(kind=irg)                      :: istat, ierr, i, j, ii, hdferr, npx, npy, numvariants, io_int(4)
integer(kind=irg),allocatable          :: sA(:), sB(:), sC(:), sD(:)
integer(kind=irg),parameter            :: numfrac = 21
integer(kind=irg)                      :: paxA(3)
integer(kind=irg)                      :: haxA(3)
real(kind=sgl)                         :: PP(3), HH(3), CCm(3)
logical                                :: verbose, iv, overwrite, isEBSD, isECP, isTKD
character(6)                           :: sqorheA, sqorheB, sqorheC, sqorheD
character(fnlen)                       :: outstr, datafile, xtalnameA, xtalnameB, xtalnameC, xtalnameD, nmldeffile
type(DynType),save                     :: DynA, DynB, DynC, DynD
type(gnode),save                       :: rlpA, rlpB, rlpC, rlpD
real(kind=sgl)                         :: dmin, voltage, TTAB(3,3), io_real(3), cA, cB, scl, fA(numfrac), om(3,3), &
                                          TTAC(3,3), TTAD(3,3), cC, cD, om2(3,3), om3(3,3), fother(numfrac), &
                                          fracA, TT(3,3)
real(kind=dbl)                         :: edge, xy(2), xyz(3), txyz(3), txy(2), txyz1(3), txyz2(3), txyz3(3), xyz1(3), xyz2(3), &
                                          xyz3(3), Radius, dc(3)
real(kind=sgl),allocatable             :: master(:,:,:), masterLC(:,:,:), masterSP(:,:,:), masterNH(:,:,:,:), &
                                          masterSH(:,:,:,:), ccA(:), ccB(:), ccC(:), ccD(:), SPNH(:,:,:), SPSH(:,:,:)
character(fnlen,kind=c_char)           :: line2(1)
integer(HSIZE_T)                       :: dims4(4), cnt4(4), offset4(4)
integer(HSIZE_T)                       :: dims3(3), cnt3(3), offset3(3)
character(fnlen)                       :: groupname, datagroupname, dataset, fname

associate( enl=>self%nml )

paxA = enl%PatternAxisA
haxA = enl%HorizontalAxisA

! determine how many variant phases there are:
numvariants = 1
if (trim(enl%masterfileC).ne.'undefined') numvariants = numvariants+1
if (trim(enl%masterfileD).ne.'undefined') numvariants = numvariants+1

isEBSD = .FALSE.
isTKD = .FALSE.
isECP = .FALSE.

select case(enl%modality)
  case('EBSD')
    isEBSD = .TRUE.
    call MPFTA%setModality('EBSD')
    call MPFTB%setModality('EBSD')
    call MPFTC%setModality('EBSD')
    call MPFTD%setModality('EBSD')
  case('TKD')
    isTKD = .TRUE.
    call MPFTA%setModality('TKD')
    call MPFTB%setModality('TKD')
    call MPFTC%setModality('TKD')
    call MPFTD%setModality('TKD')
  case('ECP')
    isECP = .TRUE.
    call MPFTA%setModality('ECP')
    call MPFTB%setModality('ECP')
    call MPFTC%setModality('ECP')
    call MPFTD%setModality('ECP')
  end select

nmldeffile = trim(EMsoft%nmldeffile)

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()

! set the HDF group names for the Monte Carlo program first
HDFnames = HDFnames_T()
call HDFnames%set_ProgramData(SC_MCOpenCL)
call HDFnames%set_NMLlist(SC_MCCLNameList)
call HDFnames%set_NMLfilename(SC_MCOpenCLNML)

! read the Monte Carlo data first for A and B
fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfileA))
call MCFTA%setFileName(fname)
call MCFTA%readMCfile(HDF, HDFnames)
mcnlA = MCFTA%getnml()
xtalnameA = trim(mcnlA%xtalname)
call cellA%setFileName(xtalnameA)
fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfileB))
call MCFTB%setFileName(fname)
call MCFTB%readMCfile(HDF, HDFnames)
mcnlB = MCFTB%getnml()
xtalnameB = trim(mcnlB%xtalname)
call cellB%setFileName(xtalnameB)

! then C and D if needed
if (numvariants.gt.1) then
  fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfileC))
  call MCFTC%setFileName(fname)
  call MCFTC%readMCfile(HDF, HDFnames)
  mcnlC = MCFTC%getnml()
  xtalnameC = trim(mcnlC%xtalname)
  call cellC%setFileName(xtalnameC)
  if (numvariants.gt.2) then
    fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfileD))
    call MCFTD%setFileName(fname)
    call MCFTD%readMCfile(HDF, HDFnames)
    mcnlD = MCFTD%getnml()
    xtalnameD = trim(mcnlD%xtalname)
    call cellD%setFileName(xtalnameD)
  end if
end if

! read master pattern files
if (isTKD.eqv..TRUE.) then
  call HDFnames%set_ProgramData(SC_TKDmaster)
  call HDFnames%set_NMLlist(SC_TKDmasterNameList)
  call HDFnames%set_NMLfilename(SC_TKDmasterNML)
end if
if (isEBSD.eqv..TRUE.) then
  call HDFnames%set_ProgramData(SC_EBSDmaster)
  call HDFnames%set_NMLlist(SC_EBSDmasterNameList)
  call HDFnames%set_NMLfilename(SC_EBSDmasterNML)
end if
if (isECP.eqv..TRUE.) then
  call HDFnames%set_ProgramData(SC_ECPmaster)
  call HDFnames%set_NMLlist(SC_ECPmasterNameList)
  call HDFnames%set_NMLfilename(SC_ECPmasterNML)
end if
call HDFnames%set_Variable(SC_MCOpenCL)

if (trim(enl%overlapmode).eq.'series') then
  fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfileA))
  call MPFTA%setFileName(fname)
  if (isEBSD) call MPFTA%readMPfile(HDF, HDFnames, mpnlEBSD, getmLPNH=.TRUE., getmLPSH=.TRUE.)
  if (isECP) call MPFTA%readMPfile(HDF, HDFnames, mpnlECP, getmLPNH=.TRUE., getmLPSH=.TRUE.)
  if (isTKD) call MPFTA%readMPfile(HDF, HDFnames, mpnlTKD, getmLPNH=.TRUE., getmLPSH=.TRUE.)
  fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfileB))
  call MPFTB%setFileName(fname)
  if (isEBSD) call MPFTB%readMPfile(HDF, HDFnames, mpnlEBSD, getmLPNH=.TRUE., getmLPSH=.TRUE.)
  if (isECP) call MPFTB%readMPfile(HDF, HDFnames, mpnlECP, getmLPNH=.TRUE., getmLPSH=.TRUE.)
  if (isTKD) call MPFTB%readMPfile(HDF, HDFnames, mpnlTKD, getmLPNH=.TRUE., getmLPSH=.TRUE.)
  allocate(sA(3), sB(3))
  sA = shape(MPFTA%MPDT%mLPNH)
  sB = shape(MPFTB%MPDT%mLPNH)
  write (*,*) maxval(MPFTA%MPDT%mLPNH), maxval(MPFTB%MPDT%mLPNH)
  if (numvariants.gt.1) then
    fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfileC))
    call MPFTC%setFileName(fname)
    if (isEBSD) call MPFTC%readMPfile(HDF, HDFnames, mpnlEBSD, getmLPNH=.TRUE., getmLPSH=.TRUE.)
    if (isECP) call MPFTC%readMPfile(HDF, HDFnames, mpnlECP, getmLPNH=.TRUE., getmLPSH=.TRUE.)
    if (isTKD) call MPFTC%readMPfile(HDF, HDFnames, mpnlTKD, getmLPNH=.TRUE., getmLPSH=.TRUE.)
    allocate(sC(3))
    sC = shape(MPFTC%MPDT%mLPNH)
    if (numvariants.gt.2) then
      fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfileD))
      call MPFTD%setFileName(fname)
      if (isEBSD) call MPFTD%readMPfile(HDF, HDFnames, mpnlEBSD, getmLPNH=.TRUE., getmLPSH=.TRUE.)
      if (isECP) call MPFTD%readMPfile(HDF, HDFnames, mpnlECP, getmLPNH=.TRUE., getmLPSH=.TRUE.)
      if (isTKD) call MPFTD%readMPfile(HDF, HDFnames, mpnlTKD, getmLPNH=.TRUE., getmLPSH=.TRUE.)
      allocate(sD(3))
      sD = shape(MPFTD%MPDT%mLPNH)
    end if
  end if
else
  ! we must keep the 4D master arrays in this case, since we want to create a correct new master pattern file
  fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfileA))
  call MPFTA%setFileName(fname)
  if (isEBSD) call MPFTA%readMPfile(HDF, HDFnames, mpnlEBSD, getmLPNH=.TRUE., getmLPSH=.TRUE., keep4=.TRUE.)
  if (isECP) call MPFTA%readMPfile(HDF, HDFnames, mpnlECP, getmLPNH=.TRUE., getmLPSH=.TRUE., keep4=.TRUE.)
  if (isTKD) call MPFTA%readMPfile(HDF, HDFnames, mpnlTKD, getmLPNH=.TRUE., getmLPSH=.TRUE., keep4=.TRUE.)
  fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfileB))
  call MPFTB%setFileName(fname)
  if (isEBSD) call MPFTB%readMPfile(HDF, HDFnames, mpnlEBSD, getmLPNH=.TRUE., getmLPSH=.TRUE., keep4=.TRUE.)
  if (isECP) call MPFTB%readMPfile(HDF, HDFnames, mpnlECP, getmLPNH=.TRUE., getmLPSH=.TRUE., keep4=.TRUE.)
  if (isTKD) call MPFTB%readMPfile(HDF, HDFnames, mpnlTKD, getmLPNH=.TRUE., getmLPSH=.TRUE., keep4=.TRUE.)
  allocate(sA(4), sB(4))
  sA = shape(MPFTA%MPDT%mLPNH4)
  sB = shape(MPFTB%MPDT%mLPNH4)
  if (numvariants.gt.1) then
    fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfileC))
    call MPFTC%setFileName(fname)
    if (isEBSD) call MPFTC%readMPfile(HDF, HDFnames, mpnlEBSD, getmLPNH=.TRUE., getmLPSH=.TRUE., keep4=.TRUE.)
    if (isECP) call MPFTC%readMPfile(HDF, HDFnames, mpnlECP, getmLPNH=.TRUE., getmLPSH=.TRUE., keep4=.TRUE.)
    if (isTKD) call MPFTC%readMPfile(HDF, HDFnames, mpnlTKD, getmLPNH=.TRUE., getmLPSH=.TRUE., keep4=.TRUE.)
    allocate(sC(4))
    sC = shape(MPFTC%MPDT%mLPNH4)
    if (numvariants.gt.2) then
      fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfileD))
      call MPFTD%setFileName(fname)
      if (isEBSD) call MPFTD%readMPfile(HDF, HDFnames, mpnlEBSD, getmLPNH=.TRUE., getmLPSH=.TRUE., keep4=.TRUE.)
      if (isECP) call MPFTD%readMPfile(HDF, HDFnames, mpnlECP, getmLPNH=.TRUE., getmLPSH=.TRUE., keep4=.TRUE.)
      if (isTKD) call MPFTD%readMPfile(HDF, HDFnames, mpnlTKD, getmLPNH=.TRUE., getmLPSH=.TRUE., keep4=.TRUE.)
      allocate(sD(4))
      sD = shape(MPFTD%MPDT%mLPNH4)
    end if
  end if
end if

! make sure that the master pattern arrays have the same dimensions
call Message%printMessage('array sizes : ')
if (trim(enl%overlapmode).eq.'series') then
  ii = 3
else
  ii = 4
end if
io_int(1:ii) = sA
call Message%WriteValue('A: ', io_int, ii)
io_int(1:ii) = sB
call Message%WriteValue('B: ', io_int, ii)
if (numvariants.gt.1) then
  io_int(1:ii) = sC
  call Message%WriteValue('C: ', io_int, ii)
  if (numvariants.gt.2) then
    io_int(1:ii) = sD
    call Message%WriteValue('D: ', io_int, ii)
  end if
end if

if (sum(sA(1:3)-sB(1:3)).ne.0) then
  call Message%printError('EMMPoverlap','master patterns A and B have different dimensions')
end if
if (numvariants.gt.1) then
  if (sum(sA(1:3)-sC(1:3)).ne.0) then
    call Message%printError('EMMPoverlap','master patterns A and C have different dimensions')
  end if
  if (numvariants.gt.2) then
    if (sum(sA(1:3)-sD(1:3)).ne.0) then
      call Message%printError('EMMPoverlap','master patterns A and D have different dimensions')
    end if
  end if
end if

! if the overlapmode is 'full', then we need to copy the phase A master pattern file in its entirety,
! and replace the master pattern arrays by overlap arrays; this can be useful for dictionary or spherical
! indexing runs when the majority of the patterns are overlap patterns between phases with a known OR.
! In that case we will also need to replace the crystal structure by a low symmetry one, possibly triclinic,
! since the symmetry of the overlap master pattern will generally be really low.
if (trim(enl%overlapmode).eq.'full') then
  if (enl%newpgnum.eq.-1) then
    call MPFTA%copyMPoverlapdata(EMsoft, HDF, HDFnames, enl%masterfileA, enl%datafile, enl%h5copypath)
  else
    call MPFTA%copyMPoverlapdata(EMsoft, HDF, HDFnames, enl%masterfileA, enl%datafile, enl%h5copypath, skipCrystalData=.TRUE.)
  end if
end if

!=============================
!=============================
! ok, we're in business... let's initialize the crystal structures so that we
! can compute the orientation relation matrix
dmin=0.05
voltage = 30000.0

call Diff%setV(dble(voltage))
call Initialize_Cell(cellA, DiffA, SGA, DynA, EMsoft, dmin, verbose)
call Initialize_Cell(cellB, DiffB, SGB, DynB, EMsoft, dmin, verbose)
if (numvariants.gt.1) then
  call Initialize_Cell(cellC, DiffC, SGC, DynC, EMsoft, dmin, verbose)
  if (numvariants.gt.2) then
    call Initialize_Cell(cellD, DiffD, SGD, DynD, EMsoft, dmin, verbose)
  end if
end if

! ensure that PatternAxisA and HorizontalAxisA are orthogonal vectors
! make sure that these two directions are orthogonal
if (abs(cellA%CalcDot(float(enl%PatternAxisA),float(enl%HorizontalAxisA),'d')).gt.1.e-6) then
  call Message%printError('MPoverlap','PatternAxisA and HorizontalAxisA must be orthogonal !!!')
end if

! check the OR for orthogonality
if (sum(enl%gA*enl%tA).ne.0.0) then
  call Message%printError('MPoverlap','gA and tA must be orthogonal !!!')
end if

if (sum(enl%gB*enl%tB).ne.0.0) then
  call Message%printError('MPoverlap','gB and tB must be orthogonal !!!')
end if
verbose = .TRUE.

orelAB%tA = enl%tA
orelAB%tB = enl%tB
orelAB%gA = enl%gA
orelAB%gB = enl%gB
TT = ComputeOR(orelAB, cellA, cellB, 'AB')
TTAB = transpose(matmul( cellA%getrsm(), matmul( TT, transpose(cellB%getdsm()) )))

! convert the vectors to the standard cartesian frame for structure A
call cellA%TransSpace(float(paxA), PP, 'd', 'c')
call cellA%TransSpace(float(haxA), HH, 'd', 'c')
call cellA%CalcCross(PP, HH, CCm, 'c', 'c', 0)

call cellA%NormVec(PP, 'c')
call cellA%NormVec(HH, 'c')
call cellA%NormVec(CCm, 'c')

! place these normalized vectors in the om array
om(1,1:3) = HH(1:3)
om(2,1:3) = CCm(1:3)
om(3,1:3) = PP(1:3)
om = transpose(om)


if (numvariants.gt.1) then
  if (sum(enl%gC*enl%tC).ne.0.0) then
    call Message%printError('MPoverlap','gC and tC must be orthogonal !!!')
  end if
  orelAC%tA = enl%tA
  orelAC%tB = enl%tC
  orelAC%gA = enl%gA
  orelAC%gB = enl%gC
  TT = ComputeOR(orelAC, cellA, cellC, 'AB')
  TTAC = transpose(matmul( cellA%getrsm(), matmul( TT, transpose(cellC%getdsm()) )))
  om2 = om

  if (numvariants.gt.2) then
    if (sum(enl%gD*enl%tD).ne.0.0) then
      call Message%printError('MPoverlap','gD and tD must be orthogonal !!!')
    end if
    orelAD%tA = enl%tA
    orelAD%tB = enl%tD
    orelAD%gA = enl%gA
    orelAD%gB = enl%gD
    TT = ComputeOR(orelAD, cellA, cellD, 'AB')
    TTAD = transpose(matmul( cellA%getrsm(), matmul( TT, transpose(cellD%getdsm()) )))
    om3 = om

  end if
end if

!=============================
!=============================
!=============================
!=============================
! next, if we are in 'series' overlapmode, allocate a new master array into which we'll write the superimposed patterns
! we discard all energies except for the highest energy

if (trim(enl%overlapmode).eq.'series') then
  npx = (sA(2)-1)/2
  npy = npx
  allocate(master(-npx:npx,-npy:npy,numfrac))

  fA = (/ (float(i)*0.05,i=0,numfrac-1) /)
  fother = (1.0 - fA) / float(numvariants)

  call Message%WriteValue('','Each master pattern has its own intensity range.',"(/A)")
  call Message%WriteValue('','This means that one pattern may dominate over another')
  call Message%WriteValue('','even when the volume fractions of A and B are equal. ',"(A/)")
  io_real(1) = maxval(MPFTA%MPDT%mLPNH)
  call Message%WriteValue('maximum intensity in master A: ',io_real, 1)
  io_real(1) = maxval(MPFTB%MPDT%mLPNH)
  call Message%WriteValue('maximum intensity in master B: ',io_real, 1)
  if (numvariants.gt.1) then
    io_real(1) = maxval(MPFTC%MPDT%mLPNH)
    call Message%WriteValue('maximum intensity in master C: ',io_real, 1)
    if (numvariants.gt.2) then
      io_real(1) = maxval(MPFTD%MPDT%mLPNH)
      call Message%WriteValue('maximum intensity in master D: ',io_real, 1)
    end if
  end if

!=============================
!=============================
  sA(1:2) = (sA(1:2)-1)/2
  sB(1:2) = (sB(1:2)-1)/2
  if (numvariants.gt.1) then
    sC(1:2) = (sC(1:2)-1)/2
  else if (numvariants.gt.2) then
          sD(1:2) = (sD(1:2)-1)/2
  end if

  edge = 1.0D0 / dble(npx)
  do i=-npx,npx
    do j=-npy,npy
! determine the spherical direction for this point
      L = Lambert_T( xyd = (/ dble(i), dble(j) /) * edge )
      ierr = L%LambertSquareToSphere(xyz)
! apply the overall pattern rotation with rotation matrix om
      xyz1 = matmul(om, xyz)
      call cellA%NormVec(xyz1, 'c')
! since A is already square Lambert, all we need to do is compute the
! beam orientation in crystal B, and sample the master pattern for that
! location.
      txyz1 = matmul(TTAB, xyz1)
! normalize these direction cosines (they are already in a cartesian reference frame!)
      txyz1 = txyz1/sqrt(sum(txyz1*txyz1))
! and interpolate the masterB pattern
      cA = InterpolateMaster(xyz1, MPFTA%MPDT%mLPNH, sA)
      cB = InterpolateMaster(txyz1, MPFTB%MPDT%mLPNH, sB)
      if (numvariants.gt.1) then
        xyz2 = matmul(om2, xyz)
        call cellA%NormVec(xyz2, 'c')
        txyz2 = matmul(TTAC, xyz2)
        txyz2 = txyz2/sqrt(sum(txyz2*txyz2))
        cC = InterpolateMaster(txyz2, MPFTC%MPDT%mLPNH, sC)
        if (numvariants.gt.2) then
          xyz3 = matmul(om3, xyz)
          call cellA%NormVec(xyz3, 'c')
          txyz3 = matmul(TTAD, xyz3)
          txyz3 = txyz3/sqrt(sum(txyz3*txyz3))
          cD = InterpolateMaster(txyz3, MPFTD%MPDT%mLPNH, sD)
        end if
      end if

      select case(numvariants)
        case(1)
           master(i,j,1:numfrac) = fA(1:numfrac)*cA + fother(1:numfrac)*cB
        case(2)
           master(i,j,1:numfrac) = fA(1:numfrac)*cA + fother(1:numfrac)*(cB+cC)
        case(3)
           master(i,j,1:numfrac) = fA(1:numfrac)*cA + fother(1:numfrac)*(cB+cC+cD)
      end select

    end do
  end do
  call Message%printMessage(' completed interpolation ')

  !=============================
  !=============================
  ! convert the square Lambert projection to a stereographic projection
  ! with the PatternAxis of structure A at the center
  allocate(masterSP(-npx:npx,-npy:npy,numfrac))
  Radius = 1.0
  do i=-npx,npx
    do j=-npy,npy
      L = Lambert_T( xyd = (/ dble(i), dble(j) /) / dble(npx) )
      ierr = L%StereoGraphicInverse( xyz, Radius )
      if (ierr.ne.0) then
        masterSP(i,j,1:numfrac) = 0.0
      else
        masterSP(i,j,1:numfrac) = OverlapInterpolateLambert(xyz, master, npx, numfrac)
      end if
    end do
  end do

  call Message%printMessage(' completed SP conversion')

  ! convert the square Lambert projection to a circular Lambert projection
  ! with the PatternAxis of structure A at the center
  allocate(masterLC(-npx:npx,-npy:npy,numfrac))
  Radius = 1.0
  do i=-npx,npx
    do j=-npy,npy
      xy = sqrt(2.0) * (/ float(i), float(j) /) / float(npx)
      if (sum(xy*xy).gt.2.0) then
        masterLC(i,j,1:numfrac) = 0.0
      else
        L = Lambert_T( xyd = xy )
        ierr = L%LambertInverse( xyz, Radius )
        masterLC(i,j,1:numfrac) = OverlapInterpolateLambert(xyz, master, npx, numfrac)
      end if
    end do
  end do

  call Message%printMessage(' completed LC conversion')

write (*,*) 'maxvals : ', maxval(master), maxval(masterLC), maxval(masterSP)

  ! finally, create simple HDF5 file with only the overlap master array in it
  ! Create a new file using the default properties.
  datafile = EMsoft%generateFilePath('EMdatapathname',trim(enl%datafile))
  hdferr =  HDF%createFile(datafile)

  ! create datasets
  dataset = 'MasterLambertSquare'
  hdferr = HDF%writeDatasetFloatArray(dataset, master, 2*npx+1, 2*npx+1, numfrac)

  dataset = 'MasterLambertCircle'
  hdferr = HDF%writeDatasetFloatArray(dataset, masterLC, 2*npx+1, 2*npx+1, numfrac)

  dataset = 'MasterStereographic'
  hdferr = HDF%writeDatasetFloatArray(dataset, masterSP, 2*npx+1, 2*npx+1, numfrac)

  call HDF%pop(.TRUE.)

  call Message%WriteValue('Output data stored in '//trim(datafile),'',"(//A/)")
else ! we're doing a full merge of master patterns...
! in this mode, we compute the complete master pattern array for a fracA volume fraction of pattern A,
! and store it in the output file which was generated earlier by an h5copy command.
  npx = (sA(2)-1)/2
  npy = npx
  allocate(masterNH(-npx:npx,-npy:npy,1:sA(3),1:sA(4)), masterSH(-npx:npx,-npy:npy,1:sA(3),1:sA(4)))
  masterNH = 0.0
  masterSH = 0.0
  allocate(ccA(sA(3)), ccB(sA(3)))
  if (numvariants.gt.1) then
    allocate(ccC(sA(3)))
    if (numvariants.gt.2) then
      allocate(ccD(sA(3)))
    end if
  end if

!=============================
!=============================
  fracA = 1.0 - enl%fracB
  sA(1:2) = (sA(1:2)-1)/2
  sB(1:2) = (sB(1:2)-1)/2
  if (numvariants.gt.1) then
    fracA = fracA - enl%fracC
    sC(1:2) = (sC(1:2)-1)/2
    if (numvariants.gt.2) then
      fracA = fracA - enl%fracD
      sD(1:2) = (sD(1:2)-1)/2
    end if
  end if

  edge = 1.0D0 / dble(npx)
  do i=-npx,npx
    do j=-npy,npy
! determine the spherical direction for this point
      xy = (/ dble(i), dble(j) /) * edge
      L = Lambert_T( xyd = (/ dble(i), dble(j) /) * edge )
      ierr = L%LambertSquareToSphere(xyz)
! apply the overall pattern rotation with rotation matrix om
      xyz1 = matmul(om, xyz)
      call cellA%NormVec(xyz1, 'c')
! since A is already square Lambert, all we need to do is compute the
! beam orientation in crystal B, and sample the master pattern for that
! location.
      txyz1 = matmul(TTAB, xyz1)
! normalize these direction cosines (they are already in a cartesian reference frame!)
      txyz1 = txyz1/sqrt(sum(txyz1*txyz1))
! and interpolate the master patterns
      ccA = sum(InterpolateCompleteMaster(xyz1, MPFTA%MPDT%mLPNH4, sA), 2)
      ccB = sum(InterpolateCompleteMaster(txyz1, MPFTB%MPDT%mLPNH4, sB), 2)
      if (numvariants.gt.1) then
        xyz2 = matmul(om2, xyz)
        call cellA%NormVec(xyz2, 'c')
        txyz2 = matmul(TTAC, xyz2)
        txyz2 = txyz2/sqrt(sum(txyz2*txyz2))
        ccC = sum(InterpolateCompleteMaster(txyz2, MPFTC%MPDT%mLPNH4, sC), 2)
        if (numvariants.gt.2) then
          xyz3 = matmul(om3, xyz)
          call cellA%NormVec(xyz3, 'c')
          txyz3 = matmul(TTAD, xyz3)
          txyz3 = txyz3/sqrt(sum(txyz3*txyz3))
          ccD = sum(InterpolateCompleteMaster(txyz3, MPFTD%MPDT%mLPNH4, sD), 2)
        end if
      end if

      select case(numvariants)
        case(1)
          masterNH(i,j,1:sA(3),1) = fracA*ccA(1:sA(3))+enl%fracB*ccB(1:sA(3))

        case(2)
          masterNH(i,j,1:sA(3),1) = fracA*ccA(1:sA(3))+enl%fracB*ccB(1:sA(3)) + &
                                                       enl%fracC*ccC(1:sA(3))

        case(3)
          masterNH(i,j,1:sA(3),1) = fracA*ccA(1:sA(3))+enl%fracB*ccB(1:sA(3)) + &
                                                       enl%fracC*ccC(1:sA(3)) + &
                                                       enl%fracD*ccD(1:sA(3))
      end select
      masterSH(-i,-j,1:sA(3),1) = masterNH(i,j,1:sA(3),1)
    end do
  end do
  call Message%printMessage(' completed interpolation ')

! recompute the stereographic projection arrays as well
  allocate(SPNH(-npx:npx,-npy:npy,1:sA(3)), SPSH(-npx:npx,-npy:npy,1:sA(3)))
  SPNH = 0.0
  SPSH = 0.0
  allocate(master(-npx:npx,-npy:npy,1:sA(3)))
  master = sum(masterNH,4)   ! sum over all the atom positions

  Radius = 1.0
  do i=-npx,npx
    do j=-npy,npy
      L = Lambert_T( xyd = (/ dble(i), dble(j) /) / dble(npx) )
      ierr = L%StereoGraphicInverse( xyz, Radius )
      if (ierr.ne.0) then
        SPNH(i,j,1:sA(3)) = 0.0
      else
        SPNH(i,j,1:sA(3)) = OverlapInterpolateLambert(xyz, master, npx, sA(3))
        SPSH(-i,-j,1:sA(3)) = SPNH(i,j,1:sA(3))
      end if
    end do
  end do

  deallocate(master)

  call Message%printMessage(' completed SP conversion')

! add the output arrays to the output file
  call HDFnames%set_ProgramData(SC_MPoverlap)
  call HDFnames%set_NMLlist(SC_MPoverlapNameList)
  call HDFnames%set_NMLfilename(SC_MPoverlapNML)

! Initialize FORTRAN HDF interface.
  overwrite = .TRUE.

! open the existing file using the default properties.
  datafile = EMsoft%generateFilePath('EMdatapathname',trim(enl%datafile))
  hdferr =  HDF%openFile(datafile)

! open or create a namelist group to write all the namelist files into
groupname = SC_NMLfiles
  hdferr = HDF%openGroup(groupname)

! read the text file and write the array to the file
dataset = trim(HDFnames%get_NMLfilename())
  hdferr = HDF%writeDatasetTextFile(dataset, nmldeffile)

! leave this group
  call HDF%pop()

! create a namelist group to write all the namelist files into
groupname = SC_NMLparameters
  hdferr = HDF%createGroup(groupname)

  call self%writeHDFNameList(HDF, HDFnames)

! leave this group
  call HDF%pop()

! if the CrystalData group was not written to the output file by h5copy,
! then we need to create it here with only the new SpaceGroupNumber
! and the PointGroupNumber as data sets (since the merged pattern may have a different
! symmetry than either/any of the member phases).
if (enl%newpgnum.ne.-1) then
  groupname = SC_CrystalData
    hdferr = HDF%createGroup(groupname)

! write the PointGroupNumber data set
  dataset = 'PointGroupNumber'
    hdferr = HDF%writeDataSetInteger(dataset, enl%newpgnum)

! write the SpaceGroupNumber data set
  dataset = SC_SpaceGroupNumber
    hdferr = HDF%writeDataSetInteger(dataset, SGPG(enl%newpgnum))

  call HDF%pop()
end if

groupname = SC_EMData
datagroupname = trim(HDFnames%get_ProgramData())
  hdferr = HDF%createGroup(groupname)
  hdferr = HDF%createGroup(datagroupname)

! add data to the hyperslab
dataset = SC_mLPNH
  dims4 = (/  2*npx+1, 2*npx+1, sA(3), sA(4) /)
  cnt4 = (/ 2*npx+1, 2*npx+1, sA(3), sA(4) /)
  offset4 = (/ 0, 0, 0, 0 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, masterNH, dims4, offset4, cnt4)

dataset = SC_mLPSH
  hdferr = HDF%writeHyperslabFloatArray(dataset, masterSH, dims4, offset4, cnt4)

  deallocate(masterNH, masterSH)

dataset = SC_masterSPNH
  dims3 = (/  2*npx+1, 2*npx+1, sA(3) /)
  cnt3 = (/ 2*npx+1, 2*npx+1, sA(3) /)
  offset3 = (/ 0, 0, 0 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, SPNH, dims3, offset3, cnt3)

dataset = SC_masterSPSH
  hdferr = HDF%writeHyperslabFloatArray(dataset, SPSH, dims3, offset3, cnt3)

  deallocate(SPNH, SPSH)

  call HDF%pop()
  call HDF%pop()

! and, at the top level of the file, add a string that states that this is a modified master pattern file
dataset = SC_READMEFIRST
  line2(1) = 'Caution: This master pattern file was generated by the EMMPoverlap program!  See the NMLparameters group. '
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)

  call HDF%pop(.TRUE.)

! and close the fortran hdf interface
  call closeFortranHDFInterface()

  call Message%printMessage(' merged patterns written to output file.')

end if

end associate

end subroutine MPoverlap_

!--------------------------------------------------------------------------
function InterpolateMaster(dc, mLPNH, s) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: InterpolateMaster
  !! author: MDG
  !! version: 1.0
  !! date: 03/23/20
  !!
  !! perform a Lambert Master Pattern interpolation

use mod_Lambert

IMPLICIT NONE

real(kind=dbl),INTENT(INOUT)            :: dc(3)
integer(kind=irg),INTENT(IN)            :: s(3)
real(kind=sgl),INTENT(IN)               :: mLPNH(-s(1):s(1),-s(2):s(2), 1:s(3))
real(kind=sgl)                          :: res

type(Lambert_T)                         :: L
integer(kind=irg)                       :: nix, niy, nixp, niyp, istat, npx, ierr
real(kind=sgl)                          :: xy(2), dx, dy, dxm, dym, scl, tmp

npx = s(1)
if (dc(3).lt.0.0) dc = -dc

! convert direction cosines to lambert projections
scl = float(npx)
L = Lambert_T( xyz = sngl(dc) )
ierr = L%LambertSphereToSquare(xy)
xy = xy * scl
res = 0.0

if (istat.eq.0) then
! interpolate intensity from the neighboring points
  nix = floor(xy(1))
  niy = floor(xy(2))
  nixp = nix+1
  niyp = niy+1
  if (nixp.gt.npx) nixp = nix
  if (niyp.gt.npx) niyp = niy
  dx = xy(1) - nix
  dy = xy(2) - niy
  dxm = 1.0 - dx
  dym = 1.0 - dy

  res = mLPNH(nix,niy,s(3))*dxm*dym + mLPNH(nixp,niy,s(3))*dx*dym + &
        mLPNH(nix,niyp,s(3))*dxm*dy + mLPNH(nixp,niyp,s(3))*dx*dy
end if

end function InterpolateMaster

!--------------------------------------------------------------------------
function OverlapInterpolateLambert(dc, master, npx, nf) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: OverlapInterpolateLambert
  !! author: MDG
  !! version: 1.0
  !! date: 03/23/20
  !!
  !! perform a Lambert Master Pattern interpolation

use mod_Lambert

IMPLICIT NONE

real(kind=dbl),INTENT(INOUT)            :: dc(3)
integer(kind=irg),INTENT(IN)            :: npx
real(kind=sgl),INTENT(IN)               :: master(-npx:npx,-npx:npx, 1:nf)
integer(kind=irg),INTENT(IN)            :: nf
real(kind=sgl)                          :: res(nf)

type(Lambert_T)                         :: L
integer(kind=irg)                       :: nix, niy, nixp, niyp, istat, ierr
real(kind=sgl)                          :: xy(2), dx, dy, dxm, dym, scl

scl = float(npx) !/ LPs%sPio2

if (dc(3).lt.0.0) dc = -dc

! convert direction cosines to lambert projections
L = Lambert_T( xyz = sngl(dc) )
ierr = L%LambertSphereToSquare(xy)
xy = xy * scl
res = 0.0

if (istat.eq.0) then
! interpolate intensity from the neighboring points
  nix = floor(xy(1))
  niy = floor(xy(2))
  nixp = nix+1
  niyp = niy+1
  if (nixp.gt.npx) nixp = nix
  if (niyp.gt.npx) niyp = niy
  dx = xy(1) - nix
  dy = xy(2) - niy
  dxm = 1.0 - dx
  dym = 1.0 - dy

  res(1:nf) = master(nix,niy,1:nf)*dxm*dym + master(nixp,niy,1:nf)*dx*dym + &
              master(nix,niyp,1:nf)*dxm*dy + master(nixp,niyp,1:nf)*dx*dy
end if

end function OverlapInterpolateLambert

!--------------------------------------------------------------------------
function InterpolateCompleteMaster(dc, mLPNH4, s) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: InterpolateCompleteMaster
  !! author: MDG
  !! version: 1.0
  !! date: 03/23/20
  !!
  !! perform a Lambert Master Pattern interpolation

use mod_Lambert

IMPLICIT NONE

real(kind=dbl),INTENT(INOUT)            :: dc(3)
integer(kind=irg),INTENT(IN)            :: s(4)
real(kind=sgl),INTENT(IN)               :: mLPNH4(-s(1):s(1),-s(2):s(2), 1:s(3), 1:s(4))
real(kind=sgl)                          :: res(s(3),s(4))

type(Lambert_T)                         :: L
integer(kind=irg)                       :: nix, niy, nixp, niyp, istat, npx, ierr
real(kind=sgl)                          :: xy(2), dx, dy, dxm, dym, scl, tmp

npx = s(1)
if (dc(3).lt.0.0) dc = -dc

! convert direction cosines to lambert projections
scl = float(npx) !/ LPs%sPio2
L = Lambert_T( xyz = sngl(dc) )
ierr = L%LambertSphereToSquare(xy)
xy = xy * scl
res = 0.0

if (istat.eq.0) then
! interpolate intensity from the neighboring points
  nix = floor(xy(1))
  niy = floor(xy(2))
  nixp = nix+1
  niyp = niy+1
  if (nixp.gt.npx) nixp = nix
  if (niyp.gt.npx) niyp = niy
  dx = xy(1) - nix
  dy = xy(2) - niy
  dxm = 1.0 - dx
  dym = 1.0 - dy

  res(1:s(3),1:s(4)) = mLPNH4(nix,niy ,1:s(3),1:s(4))*dxm*dym + mLPNH4(nixp, niy,1:s(3),1:s(4))*dx*dym+&
                       mLPNH4(nix,niyp,1:s(3),1:s(4))*dxm*dy  + mLPNH4(nixp,niyp,1:s(3),1:s(4))*dx*dy
end if

end function InterpolateCompleteMaster





end module mod_MPoverlap
