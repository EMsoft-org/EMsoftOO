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

module mod_DIfiles
  !! author: MDG
  !! version: 1.0
  !! date: 04/03/20
  !!
  !! Class definition for Dictionary Indexing programs.
  !! Contains everything dealing with dictionary indexing namelists.
  !! Also has the read routines for the HDF5 dot product file format.

use mod_kinds
use mod_global
use HDF5
use h5im
use h5lt
use mod_HDFsupport
use stringconstants

IMPLICIT NONE

! parent namelist for dictionary indexing programs; this is the minimum needed
! for EBSD indexing, but other modalities may need additional parameters that
! are defined via inherited classes
type, public :: DictionaryIndexingNameListType
  integer(kind=irg)  :: ncubochoric
  integer(kind=irg)  :: numexptsingle
  integer(kind=irg)  :: numdictsingle
  integer(kind=irg)  :: ipf_ht
  integer(kind=irg)  :: ipf_wd
  integer(kind=irg)  :: ROI(4)
  integer(kind=irg)  :: nnk
  integer(kind=irg)  :: nnav
  integer(kind=irg)  :: nosm
  integer(kind=irg)  :: nism
  integer(kind=irg)  :: maskradius
  integer(kind=irg)  :: numsx
  integer(kind=irg)  :: numsy
  integer(kind=irg)  :: exptnumsx
  integer(kind=irg)  :: exptnumsy
  integer(kind=irg)  :: binning
  integer(kind=irg)  :: nthreads
  integer(kind=irg)  :: devid
  integer(kind=irg)  :: usenumd
  integer(kind=irg)  :: multidevid(8)
  integer(kind=irg)  :: platid
  integer(kind=irg)  :: nregions
  integer(kind=irg)  :: nlines
  integer(kind=irg)  :: sw
  real(kind=sgl)     :: L
  real(kind=sgl)     :: thetac
  real(kind=sgl)     :: delta
  real(kind=sgl)     :: omega
  real(kind=sgl)     :: xpc
  real(kind=sgl)     :: ypc
  real(kind=sgl)     :: isangle
  real(kind=sgl)     :: energymin
  real(kind=sgl)     :: energymax
  real(kind=sgl)     :: gammavalue
  real(kind=sgl)     :: axisangle(4)
  real(kind=sgl)     :: beamcurrent
  real(kind=sgl)     :: dwelltime
  real(kind=sgl)     :: hipassw
  real(kind=sgl)     :: stepX
  real(kind=sgl)     :: stepY
  real(kind=sgl)     :: lambda
  logical            :: doNLPAR
  character(1)       :: maskpattern
  character(3)       :: scalingmode
  character(3)       :: Notify
  character(3)       :: similaritymetric
  character(1)       :: keeptmpfile
  character(1)       :: usetmpfile
  character(fnlen)   :: exptfile
  character(fnlen)   :: masterfile
  character(fnlen)   :: energyfile
  character(fnlen)   :: datafile
  character(fnlen)   :: tmpfile
  character(fnlen)   :: ctffile
  character(fnlen)   :: avctffile
  character(fnlen)   :: angfile
  character(fnlen)   :: eulerfile
  character(fnlen)   :: dictfile
  character(fnlen)   :: maskfile
  character(fnlen)   :: indexingmode
  character(fnlen)   :: refinementNMLfile
  character(fnlen)   :: inputtype
  character(fnlen)   :: HDFstrings(10)
  character(fnlen)   :: DIModality
  ! ECP parameters
  real(kind=sgl)     :: workingdistance
  real(kind=sgl)     :: Rin
  real(kind=sgl)     :: Rout
  real(kind=sgl)     :: conesemiangle
  real(kind=sgl)     :: sampletilt
  integer(kind=irg)  :: npix
end type DictionaryIndexingNameListType

type, public, extends(DictionaryIndexingNameListType) :: EBSDDINameListType
end type EBSDDINameListType

type, public, extends(DictionaryIndexingNameListType) :: ECPDINameListType
end type ECPDINameListType

type, public, extends(DictionaryIndexingNameListType) :: TKDDINameListType
end type TKDDINameListType


type, public :: DIdataType
  integer(kind=irg)             :: FZcnt
  integer(kind=irg)             :: Nexp
  integer(kind=irg)             :: pgnum
  real(kind=sgl)                :: MCsig
  integer(kind=sgl),allocatable :: ADP(:,:)
  real(kind=sgl),allocatable    :: AverageOrientations(:,:)
  real(kind=sgl),allocatable    :: CI(:)
  real(kind=sgl),allocatable    :: EulerAngles(:,:)
  real(kind=sgl),allocatable    :: DictionaryEulerAngles(:,:)
  real(kind=sgl),allocatable    :: Fit(:)
  real(kind=sgl),allocatable    :: IQ(:)
  real(kind=sgl),allocatable    :: KAM(:,:)
  real(kind=sgl),allocatable    :: OSM(:,:)
  integer(kind=irg),allocatable :: Phase(:)
  real(kind=sgl),allocatable    :: Phi1(:)
  real(kind=sgl),allocatable    :: Phi(:)
  real(kind=sgl),allocatable    :: Phi2(:)
  integer(kind=irg),allocatable :: SEMsignal(:)
  real(kind=sgl),allocatable    :: TopDotProductList(:,:)
  integer(kind=irg),allocatable :: TopMatchIndices(:,:)
  integer(kind=irg),allocatable :: Valid(:)
  real(kind=sgl),allocatable    :: XPosition(:)
  real(kind=sgl),allocatable    :: YPosition(:)
  real(kind=sgl),allocatable    :: RefinedEulerAngles(:,:)
  real(kind=sgl),allocatable    :: RefinedDotProducts(:)
  logical                       :: orthocomment
end type DIdataType

! class definition
type, public :: DIfile_T
private
  character(fnlen)                              :: DIfile
  type(DIdataType),public                       :: DIDT
  character(fnlen)                              :: Modality = 'unknown'
  type(DictionaryIndexingNameListType), public  :: nml

contains
private

  procedure, pass(self) :: getrefinementfilename_
  procedure, pass(self) :: getfilename_
  procedure, pass(self) :: setfilename_
  procedure, pass(self) :: getModality_
  procedure, pass(self) :: setModality_
  procedure, pass(self) :: readDIModality_
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: readDotProductFile_
  procedure, pass(self) :: h5_writeFile_
  procedure, pass(self) :: h5_writeInfo

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: getrefinementfilename => getrefinementfilename_
  generic, public :: getfilename => getfilename_
  generic, public :: setfilename => setfilename_
  generic, public :: getModality => getModality_
  generic, public :: setModality => setModality_
  generic, public :: readDIModality => readDIModality_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readDotProductFile => readDotProductFile_
  generic, public :: h5_writeFile => h5_writeFile_

end type DIfile_T

private :: h5_writePhaseGroup, h5_write2DImageFromVector, h5_writeCoordinateSystemGroup, &
           h5_writePatternCenterGroup

! the constructor routine for this class
interface DIfile_T
  module procedure DIfile_constructor
end interface DIfile_T

contains

!--------------------------------------------------------------------------
type(DIfile_T) function DIfile_constructor( nmlfile, fname, inRAM ) result(DIfile)
!DEC$ ATTRIBUTES DLLEXPORT :: DIfile_constructor
!! author: MDG
!! version: 1.0
!! date: 04/03/20
!!
!! constructor for the DIfile_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile
character(fnlen), OPTIONAL   :: fname
logical, OPTIONAL            :: inRAM 

if (present(nmlfile)) then 
  if (present(inRAM)) then 
    if (inRAM.eqv..TRUE.) call DIfile%readNameList(nmlfile, inRAM=.TRUE.)
  else
    call DIfile%readNameList(nmlfile)
  end if 
end if 
if (present(fname)) call DIfile%setfilename(fname)

end function DIfile_constructor

!--------------------------------------------------------------------------
subroutine DIfile_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: DIfile_destructor
!! author: MDG
!! version: 1.0
!! date: 04/03/20
!!
!! destructor for the DIfile_T Class

IMPLICIT NONE

type(DIfile_T), INTENT(INOUT)  :: self

call reportDestructor('DIfile_T')

end subroutine DIfile_destructor

!--------------------------------------------------------------------------
function getrefinementfilename_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getrefinementfilename_
!! author: MDG
!! version: 1.0
!! date: 11/07/21
!!
!! get refinement filename from the DIfile_T class

IMPLICIT NONE

class(DIfile_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%refinementNMLfile)

end function getrefinementfilename_

!--------------------------------------------------------------------------
function getfilename_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getfilename_
!! author: MDG
!! version: 1.0
!! date: 04/03/20
!!
!! get filename from the DIfile_T class

IMPLICIT NONE

class(DIfile_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%DIfile)

end function getfilename_

!--------------------------------------------------------------------------
subroutine setfilename_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setfilename_
!! author: MDG
!! version: 1.0
!! date: 04/03/20
!!
!! set filename in the DIfile_T class

IMPLICIT NONE

class(DIfile_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%DIfile = trim(inp)

end subroutine setfilename_

!--------------------------------------------------------------------------
function getModality_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getModality_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! get Modality from the DIfile_T class

IMPLICIT NONE

class(DIfile_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%Modality)

end function getModality_

!--------------------------------------------------------------------------
subroutine setModality_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setModality_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! set Modality in the DIfile_T class

use mod_io

IMPLICIT NONE

class(DIfile_T), INTENT(INOUT)     :: self
character(*), INTENT(IN)           :: inp

type(IO_T)                         :: Message
self%Modality = trim(inp)

! call Message%printMessage('DIFT%setModality: setting modality to '//trim(inp)//trim(self%Modality))

end subroutine setModality_

!--------------------------------------------------------------------------
subroutine readDIModality_(self, HDF, DIfile)
!DEC$ ATTRIBUTES DLLEXPORT :: readDIModality_
!! author: MDG
!! version: 1.0
!! date: 04/21/20
!!
!! read the DIModality parameter from the DI file and set it in the DIfile_T class

use HDF5
use mod_HDFsupport
use mod_io
use stringconstants
use ISO_C_BINDING

IMPLICIT NONE

class(DIfile_T), INTENT(INOUT)    :: self
type(HDF_T), INTENT(INOUT)        :: HDF
character(fnlen), INTENT(IN)      :: DIfile

type(IO_T)                        :: Message
character(fnlen)                  :: dataset
logical                           :: f_exists, g_exists, stat
integer(kind=irg)                 :: hdferr, nlines
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)

! we assume that DIfile contains the full path to the master pattern file
inquire(file=trim(DIfile), exist=f_exists)

if (.not.f_exists) then
  call Message%printError('readDIModality','Dot product file '//trim(DIfile)//' does not exist')
end if

! is this a proper HDF5 file ?
call h5fis_hdf5_f(trim(DIfile), stat, hdferr)
if (stat.eqv..FALSE.) then
  call Message%printError('readDIModality','This is not an HDF5 file.')
end if

! open the file
hdferr =  HDF%openFile(DIfile)

! check whether or not the DIModality data set exists at the top level
dataset = SC_DIModality
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
! read the DIModality string from the top level (only present for EMsoft 6.X versions)
  call HDF%readDatasetStringArray(dataset, nlines, hdferr, stringarray)
  if (hdferr.ne.0) call HDF%error_check('readDIModality: unable to read DIModality dataset', hdferr)
  self%nml%DIModality = trim(stringarray(1))
  self%Modality = trim(stringarray(1))
  deallocate(stringarray)
else
  call Message%printMessage(' readDIModality: this file does not contain a DIModality data set.',"(/A)")
  call Message%printMessage('  --> program will continue assuming this is an old dot product file',"(A/)")
end if

! close the file
call HDF%pop(.TRUE.)

end subroutine readDIModality_


!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly, inRAM)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! read the namelist from an nml file for the DIfile_T Class

use mod_io
use mod_EMsoft

IMPLICIT NONE

class(DIfile_T), INTENT(INOUT)      :: self
character(fnlen),INTENT(IN)         :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)         :: initonly
 !! fill in the default values only; do not read the file
logical,OPTIONAL,INTENT(IN)         :: inRAM

type(EMsoft_T)                      :: EMsoft
type(IO_T)                          :: Message
logical                             :: skipread = .FALSE.

integer(kind=irg)  :: numsx
integer(kind=irg)  :: numsy
integer(kind=irg)  :: exptnumsx
integer(kind=irg)  :: exptnumsy
integer(kind=irg)  :: ROI(4)
integer(kind=irg)  :: binning
integer(kind=irg)  :: devid
integer(kind=irg)  :: multidevid(8)
integer(kind=irg)  :: usenumd
integer(kind=irg)  :: platid
integer(kind=irg)  :: nregions
integer(kind=irg)  :: nlines
integer(kind=irg)  :: nthreads
integer(kind=irg)  :: ncubochoric
integer(kind=irg)  :: numexptsingle
integer(kind=irg)  :: numdictsingle
integer(kind=irg)  :: ipf_ht
integer(kind=irg)  :: ipf_wd
integer(kind=irg)  :: nnk
integer(kind=irg)  :: nnav
integer(kind=irg)  :: nosm
integer(kind=irg)  :: nism
integer(kind=irg)  :: maskradius
integer(kind=irg)  :: sw
real(kind=sgl)     :: L
real(kind=sgl)     :: thetac
real(kind=sgl)     :: delta
real(kind=sgl)     :: xpc
real(kind=sgl)     :: ypc
real(kind=sgl)     :: isangle
real(kind=sgl)     :: gammavalue
real(kind=sgl)     :: omega
real(kind=sgl)     :: stepX
real(kind=sgl)     :: stepY
real(kind=sgl)     :: energymin
real(kind=sgl)     :: energymax
real(kind=sgl)     :: beamcurrent
real(kind=sgl)     :: dwelltime
real(kind=sgl)     :: hipassw
real(kind=sgl)     :: lambda
logical            :: doNLPAR
character(1)       :: maskpattern
character(1)       :: keeptmpfile
character(1)       :: usetmpfile
character(3)       :: scalingmode
character(3)       :: Notify
character(3)       :: similaritymetric
character(fnlen)   :: dotproductfile
character(fnlen)   :: masterfile
character(fnlen)   :: tmpfile
character(fnlen)   :: datafile
character(fnlen)   :: ctffile
character(fnlen)   :: avctffile
character(fnlen)   :: angfile
character(fnlen)   :: eulerfile
character(fnlen)   :: inputtype
character(fnlen)   :: HDFstrings(10)
character(fnlen)   :: refinementNMLfile
character(fnlen)   :: exptfile
character(fnlen)   :: dictfile
character(fnlen)   :: maskfile
character(fnlen)   :: indexingmode
character(fnlen)   :: DIModality
! ECP parameters
real(kind=sgl)     :: workingdistance
real(kind=sgl)     :: Rin
real(kind=sgl)     :: Rout
real(kind=sgl)     :: conesemiangle
real(kind=sgl)     :: sampletilt
integer(kind=irg)  :: npix


! define the IO namelist to facilitate passing variables to the program.
namelist  / DIdata / thetac, delta, numsx, numsy, xpc, ypc, masterfile, devid, platid, inputtype, DIModality, &
                     beamcurrent, dwelltime, binning, gammavalue, energymin, nregions, nlines, maskfile, &
                     scalingmode, maskpattern, L, omega, nthreads, energymax, datafile, angfile, ctffile, &
                     ncubochoric, numexptsingle, numdictsingle, ipf_ht, ipf_wd, nnk, nnav, exptfile, maskradius, &
                     dictfile, indexingmode, hipassw, stepX, stepY, tmpfile, avctffile, nosm, eulerfile, Notify, &
                     HDFstrings, ROI, keeptmpfile, multidevid, usenumd, nism, isangle, refinementNMLfile, &
                     workingdistance, Rin, Rout, conesemiangle, sampletilt, npix, doNLPAR, sw, lambda, similaritymetric, &
                     exptnumsx, exptnumsy, usetmpfile

namelist  / DIRAMdata / thetac, delta, numsx, numsy, xpc, ypc, masterfile, devid, platid, inputtype, DIModality, &
                     beamcurrent, dwelltime, binning, gammavalue, energymin, nregions, nlines, maskfile, &
                     scalingmode, maskpattern, L, omega, nthreads, energymax, datafile, angfile, ctffile, &
                     ncubochoric, numexptsingle, numdictsingle, ipf_ht, ipf_wd, nnk, nnav, exptfile, maskradius, &
                     dictfile, indexingmode, hipassw, stepX, stepY, tmpfile, avctffile, nosm, eulerfile, Notify, &
                     HDFstrings, ROI, keeptmpfile, multidevid, usenumd, nism, isangle, refinementNMLfile, &
                     workingdistance, Rin, Rout, conesemiangle, sampletilt, npix, doNLPAR, sw, lambda

! set the input parameters to default values (except for xtalname, which must be present)
ncubochoric     = 50
numexptsingle   = 1024
numdictsingle   = 1024
platid          = 1
devid           = 1
usenumd         = 1
multidevid      = (/ 0, 0, 0, 0, 0, 0, 0, 0 /)
nregions        = 10
nlines          = 3
nnk             = 50
nnav            = 20
nosm            = 20
nism            = 5
exptfile        = 'undefined'
numsx           = 0             ! [dimensionless]
numsy           = 0             ! [dimensionless]
exptnumsx       = 0             ! [dimensionless]
exptnumsy       = 0             ! [dimensionless]
ROI             = (/ 0, 0, 0, 0 /)  ! Region of interest (/ x0, y0, w, h /)
maskradius      = 240
binning         = 1             ! binning mode  (1, 2, 4, or 8)
sw              = 3
L               = 20000.0       ! [microns]
thetac          = 0.0           ! [degrees]
delta           = 25.0          ! [microns]
xpc             = 0.0           ! [pixels]
ypc             = 0.0           ! [pixels]
gammavalue      = 1.0           ! gamma factor
isangle         = 1.5
beamcurrent     = 14.513        ! beam current (actually emission current) in nano ampere
dwelltime       = 100.0         ! in microseconds
hipassw         = 0.05          ! hi pass inverted Gaussian mask parameter
stepX           = 1.0           ! sampling step size along X
stepY           = 1.0           ! sampling step size along Y
lambda          = 0.375
doNLPAR         = .FALSE.
keeptmpfile     = 'n'
usetmpfile      = 'n'
maskpattern     = 'n'           ! 'y' or 'n' to include a circular mask
Notify          = 'Off'
scalingmode     = 'not'         ! intensity selector ('lin', 'gam', or 'not')
masterfile      = 'undefined'   ! filename
dotproductfile  = 'undefined'
energymin       = 10.0
energymax       = 20.0
ipf_ht          = 100
ipf_wd          = 100
nthreads        = 1
datafile        = 'undefined'
ctffile         = 'undefined'
avctffile       = 'undefined'
angfile         = 'undefined'
eulerfile       = 'undefined'
omega           = 0.0
tmpfile         = 'EMEBSDDict_tmp.data'
dictfile        = 'undefined'
maskfile        = 'undefined'
refinementNMLfile = 'undefined'
indexingmode    = 'dynamic'
similaritymetric= 'ndp'
inputtype       = 'Binary'    ! Binary, EMEBSD, TSLHDF, TSLup2, OxfordHDF, OxfordBinary, BrukerHDF
HDFstrings      = ''
DIModality      = 'EBSD'      ! EBSD, TKD, ECP, ...
! ECP
npix            = 256
conesemiangle   = 5.0
workingdistance = 13.0
Rin             = 2.0
Rout            = 6.0

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    if (present(inRAM)) then
      if (inRAM.eqv..TRUE.) read(UNIT=dataunit,NML=DIRAMdata)
    else
      read(UNIT=dataunit,NML=DIdata)
    end if 
    close(UNIT=dataunit,STATUS='keep')

    if (trim(indexingmode) .eq. 'static') then
        if (trim(dictfile) .eq. 'undefined') then
            call Message%printError('readNameList:',' dictionary file name is undefined in '//nmlfile)
        end if
    end if

! check for required entries
    if (trim(indexingmode) .eq. 'dynamic') then
        if (trim(masterfile).eq.'undefined') then
            call Message%printError('readNameList:',' master pattern file name is undefined in '//nmlfile)
        end if
    end if

    if (trim(exptfile).eq.'undefined') then
        call Message%printError('readNameList:',' experimental file name is undefined in '//nmlfile)
    end if

    if (exptnumsx.eq.0) then
        call Message%printError('readNameList:',' pattern size exptnumsx is zero in '//nmlfile)
    end if

    if (exptnumsy.eq.0) then
        call Message%printError('readNameList:',' pattern size exptnumsy is zero in '//nmlfile)
    end if

end if

! if we get here, then all appears to be ok, and we need to fill in the enl fields

self%nml%devid         = devid
self%nml%multidevid    = multidevid
self%nml%usenumd       = usenumd
self%nml%platid        = platid
self%nml%nregions      = nregions
self%nml%nlines        = nlines
self%nml%maskpattern   = maskpattern
self%nml%keeptmpfile   = keeptmpfile
self%nml%usetmpfile    = usetmpfile
self%nml%exptfile      = trim(exptfile)
self%nml%nnk           = nnk
self%nml%nnav          = nnav
self%nml%nosm          = nosm
self%nml%nism          = nism
self%nml%isangle       = isangle
self%nml%ipf_ht        = ipf_ht
self%nml%ipf_wd        = ipf_wd
self%nml%nthreads      = nthreads
self%nml%sw            = sw 
self%nml%lambda        = lambda 
self%nml%doNLPAR       = doNLPAR
self%nml%datafile      = trim(datafile)
self%nml%tmpfile       = trim(tmpfile)
self%nml%ctffile       = trim(ctffile)
self%nml%avctffile     = trim(avctffile)
self%nml%angfile       = trim(angfile)
self%nml%eulerfile     = trim(eulerfile)
self%nml%maskradius    = maskradius
self%nml%numdictsingle = numdictsingle
self%nml%numexptsingle = numexptsingle
self%nml%hipassw       = hipassw
self%nml%masterfile    = trim(masterfile)
self%nml%energyfile    = trim(masterfile)
self%nml%maskfile      = trim(maskfile)
self%nml%StepX         = stepX
self%nml%StepY         = stepY
self%nml%indexingmode  = trim(indexingmode)
self%nml%similaritymetric = similaritymetric
self%nml%Notify        = Notify
self%nml%inputtype     = inputtype
self%nml%HDFstrings    = HDFstrings
self%nml%L             = L
self%nml%numsx         = numsx
self%nml%numsy         = numsy
self%nml%exptnumsx     = exptnumsx
self%nml%exptnumsy     = exptnumsy
self%nml%ROI           = ROI
self%nml%binning       = binning
self%nml%thetac        = thetac
self%nml%delta         = delta
self%nml%xpc           = xpc
self%nml%ypc           = ypc
self%nml%gammavalue    = gammavalue
self%nml%beamcurrent   = beamcurrent
self%nml%dwelltime     = dwelltime
self%nml%scalingmode   = scalingmode
self%nml%ncubochoric   = ncubochoric
self%nml%omega         = omega
self%nml%energymin     = energymin
self%nml%energymax     = energymax
self%nml%DIModality    = DIModality
self%nml%dictfile      = trim(dictfile)
self%nml%refinementNMLfile = refinementNMLfile
! ECP
self%nml%npix            = npix
self%nml%conesemiangle   = conesemiangle
self%nml%workingdistance = workingdistance
self%nml%Rin             = Rin
self%nml%Rout            = Rout

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! pass the namelist for the DI_T Class to the calling program

IMPLICIT NONE

class(DIfile_T), INTENT(INOUT)        :: self
type(DictionaryIndexingNameListType)  :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames, emnl)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants
use mod_IO
use ISO_C_BINDING

IMPLICIT NONE

class(DIfile_T), INTENT(INOUT)                      :: self
type(HDF_T), INTENT(INOUT)                          :: HDF
type(HDFnames_T), INTENT(INOUT)                     :: HDFnames
class(DictionaryIndexingNameListType), INTENT(INOUT):: emnl

type(IO_T)                                          :: Message
integer(kind=irg)                                   :: n_int, n_real, NLPAR
integer(kind=irg)                                   :: hdferr
integer(kind=irg),allocatable                       :: io_int(:)
real(kind=sgl),allocatable                          :: io_real(:)
character(20),allocatable                           :: intlist(:), reallist(:)
character(fnlen)                                    :: dataset, sval(1),groupname, modality
character(fnlen,kind=c_char)                        :: line2(1), line10(10)
logical                                             :: g_exists, overwrite=.TRUE., isEBSD=.FALSE., &
                                                       isECP=.FALSE., isTKD=.FALSE., isEBSDSHT=.FALSE.
! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! in the future we might decide to use inherited classes for the different
! modalities but for now we don't; but we leave appropriate code in place...
modality = trim(self%getModality())

select case(trim(modality))
  case('EBSD')
    isEBSD = .TRUE.
    n_int = 22
    n_real = 20
    allocate( io_int(n_int), intlist(n_int), io_real(n_real), reallist(n_real) )
  case('ECP')
    isECP = .TRUE.
    n_int = 20
    n_real = 19
    allocate( io_int(n_int), intlist(n_int), io_real(n_real), reallist(n_real) )
  case('TKD')
    isTKD = .TRUE.
    n_int = 20
    n_real = 19
    allocate( io_int(n_int), intlist(n_int), io_real(n_real), reallist(n_real) )
  case default
    call Message%printError('writeHDFNameList', 'unknown name list type requested')
end select

NLPAR = 0
if (emnl%doNLPAR.eqv..TRUE.) NLPAR=1

! write all the single integers
io_int = (/ emnl%ncubochoric, emnl%numexptsingle, emnl%numdictsingle, emnl%ipf_ht, &
            emnl%ipf_wd, emnl%nnk, emnl%maskradius, emnl%numsx, emnl%numsy, emnl%binning, &
            emnl%nthreads, emnl%devid, emnl%platid, emnl%nregions, emnl%nnav, &
            emnl%nosm, emnl%nlines, emnl%usenumd, emnl%nism, emnl%npix, emnl%sw, NLPAR /)
intlist(1) = 'Ncubochoric'
intlist(2) = 'numexptsingle'
intlist(3) = 'numdictsingle'
intlist(4) = 'ipf_ht'
intlist(5) = 'ipf_wd '
intlist(6) = 'nnk'
intlist(7) = 'maskradius'
intlist(8) = 'numsx'
intlist(9) = 'numsy'
intlist(10) = 'binning'
intlist(11) = 'nthreads'
intlist(12) = 'devid'
intlist(13) = 'platid'
intlist(14) = 'nregions'
intlist(15) = 'nnav'
intlist(16) = 'nosm'
intlist(17) = 'nlines'
intlist(18) = 'usenumd'
intlist(19) = 'nism'
intlist(20) = 'npix'
intlist(21) = 'sw'
intlist(22) = 'NLPAR'
call HDF%writeNMLintegers(io_int, intlist, n_int)

io_real = (/ emnl%L, emnl%thetac, emnl%delta, emnl%omega, emnl%xpc, &
             emnl%ypc, emnl%energymin, emnl%energymax, emnl%gammavalue, emnl%StepX, &
             emnl%stepY, emnl%isangle, emnl%beamcurrent, emnl%dwelltime, emnl%hipassw, &
             emnl%workingdistance, emnl%conesemiangle, emnl%Rin, emnl%Rout, emnl%lambda /)
reallist(1) = 'L'
reallist(2) = 'thetac'
reallist(3) = 'delta'
reallist(4) = 'omega'
reallist(5) = 'xpc'
reallist(6) = 'ypc'
reallist(7) = 'energymin'
reallist(8) = 'energymax'
reallist(9) = 'gammavalue'
reallist(10) = 'stepX'
reallist(11) = 'stepY'
reallist(12) = 'isangle'
reallist(13) = 'beamcurrent'
reallist(14) = 'dwelltime'
reallist(15) = 'hipassw'
reallist(16) = 'workingdistance'
reallist(17) = 'conesemiangle'
reallist(18) = 'Rin'
reallist(19) = 'Rout'
reallist(20) = 'lambda'
call HDF%writeNMLreals(io_real, reallist, n_real)

! a 4-vector
dataset = SC_axisangle
hdferr = HDF%writeDatasetFloatArray(dataset, emnl%axisangle, 4)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create axisangle dataset', hdferr)

! an integer 4-vector
dataset = SC_ROI
hdferr = HDF%writeDatasetIntegerArray(dataset, emnl%ROI, 4)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create ROI dataset', hdferr)

! an integer 8-vector
dataset = 'multidevid'
hdferr = HDF%writeDatasetIntegerArray(dataset, emnl%multidevid, 8)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create multidevid dataset', hdferr)

! strings
dataset = SC_maskpattern
line2(1) = emnl%maskpattern
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create maskpattern dataset', hdferr)

dataset = SC_scalingmode
line2(1) = emnl%scalingmode
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create scalingmode dataset', hdferr)

dataset = SC_exptfile
line2(1) = emnl%exptfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create exptfile dataset', hdferr)

dataset = SC_masterfile
line2(1) = emnl%masterfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create masterfile dataset', hdferr)

dataset = SC_energyfile
line2(1) = emnl%energyfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create energyfile dataset', hdferr)

dataset = SC_datafile
line2(1) = emnl%datafile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create datafile dataset', hdferr)

dataset = SC_tmpfile
line2(1) = emnl%tmpfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create tmpfile dataset', hdferr)

dataset = SC_ctffile
line2(1) = emnl%ctffile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create ctffile dataset', hdferr)

dataset = SC_angfile
line2(1) = emnl%angfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create angfile dataset', hdferr)

dataset = SC_eulerfile
line2(1) = emnl%eulerfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create eulerfile dataset', hdferr)

dataset = SC_maskfile
line2(1) = emnl%maskfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create maskfile dataset', hdferr)

dataset = SC_inputtype
line2(1) = emnl%inputtype
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create inputtype dataset', hdferr)

dataset = SC_HDFstrings
line10 = emnl%HDFstrings
hdferr = HDF%writeDatasetStringArray(dataset, line10, 10)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create HDFstrings dataset', hdferr)

! pop this group off the stack
call HDF%pop()
call HDF%pop()

! and write the DIModality string at the top level (only present for EMsoft 6.X versions)
dataset = SC_DIModality
line2(1) = emnl%DIModality
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create DIModality dataset', hdferr)

! return to the correct group
hdferr = HDF%openGroup(HDFnames%get_NMLparameters())

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
recursive subroutine readDotProductFile_(self, EMsoft, HDF, HDFnames, dpfile, hdferr, getADP, getAverageOrientations, getCI, &
                                         getEulerAngles, getFit, getIQ, getKAM, getOSM, getPhase, getPhi1, &
                                         getPhi, getPhi2, getSEMsignal, getTopDotProductList, getTopMatchIndices, &
                                         getValid, getXPosition, getYPosition, getRefinedDotProducts, &
                                         getRefinedEulerAngles, getDictionaryEulerAngles)
!DEC$ ATTRIBUTES DLLEXPORT :: readDotProductFile_
!! author: MDG
!! version: 1.0
!! date: 04/02/20
!!
!! read a Dot Product File from Dictionary Indexing into the correct namelist and data structure

use mod_EMsoft
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_io
use ISO_C_BINDING

IMPLICIT NONE

class(DIfile_T),INTENT(INOUT)                       :: self
type(EMsoft_T), INTENT(INOUT)                       :: EMsoft
type(HDF_T),INTENT(INOUT)                           :: HDF
type(HDFnames_T),INTENT(INOUT)                      :: HDFnames
character(fnlen),INTENT(IN)                         :: dpfile
integer(kind=irg),INTENT(OUT)                       :: hdferr
logical,INTENT(IN),OPTIONAL                         :: getADP
logical,INTENT(IN),OPTIONAL                         :: getAverageOrientations
logical,INTENT(IN),OPTIONAL                         :: getCI
logical,INTENT(IN),OPTIONAL                         :: getEulerAngles
logical,INTENT(IN),OPTIONAL                         :: getDictionaryEulerAngles
logical,INTENT(IN),OPTIONAL                         :: getFit
logical,INTENT(IN),OPTIONAL                         :: getIQ
logical,INTENT(IN),OPTIONAL                         :: getKAM
logical,INTENT(IN),OPTIONAL                         :: getOSM
logical,INTENT(IN),OPTIONAL                         :: getPhase
logical,INTENT(IN),OPTIONAL                         :: getPhi1
logical,INTENT(IN),OPTIONAL                         :: getPhi
logical,INTENT(IN),OPTIONAL                         :: getPhi2
logical,INTENT(IN),OPTIONAL                         :: getSEMsignal
logical,INTENT(IN),OPTIONAL                         :: getTopDotProductList
logical,INTENT(IN),OPTIONAL                         :: getTopMatchIndices
logical,INTENT(IN),OPTIONAL                         :: getValid
logical,INTENT(IN),OPTIONAL                         :: getXPosition
logical,INTENT(IN),OPTIONAL                         :: getYPosition
logical,INTENT(IN),OPTIONAL                         :: getRefinedDotProducts
logical,INTENT(IN),OPTIONAL                         :: getRefinedEulerAngles

type(IO_T)                                          :: Message
type(HDFnames_T)                                    :: saveHDFnames

character(fnlen)                                    :: infile, groupname, dataset, tmpnmlname, Modality
logical                                             :: stat, readonly, g_exists, h_exists
character(3)                                        :: DIModality
integer(kind=irg)                                   :: ii, nlines, i
integer(kind=irg),allocatable                       :: iarray(:)
real(kind=sgl),allocatable                          :: farray(:)
integer(HSIZE_T)                                    :: dims(1), dims2(2), dims3(3), offset3(3), sz(1)
character(fnlen, KIND=c_char),allocatable,TARGET    :: stringarray(:)


associate(DIDT=>self%DIDT, ebsdnl=>self%nml)

! we assume that the calling program has opened the HDF interface,
! and that it passes the full path filename to this routine.

! is this a proper HDF5 file ?
call h5fis_hdf5_f(trim(dpfile), stat, hdferr)

!===================================================================================
!===============read dot product file===============================================
!===================================================================================

if (stat.eqv..FALSE.) then ! the file exists, so let's open it an first make sure it is an EBSD dot product file
   call Message%printError('readDotProductFile','This is not a proper HDF5 file')
end if

! open the dot product file
readonly = .TRUE.
hdferr =  HDF%openFile(dpfile, readonly)

! check the modality for this DI file...
! Starting with EMsoft 6.X, the DIModality string is present at the top level.
! If it is not present, then we have an older DI file, which different group names...
DIModality = 'new'
dataset = SC_DIModality
  call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists.eqv..FALSE.) then
    DIModality = 'old'
    Modality = 'EBSD'
    saveHDFnames = HDFnames
    call HDFnames%set_NMLfiles(SC_NMLfiles)
    call HDFnames%set_NMLfilename(SC_EBSDDictionaryIndexingNML)
    call HDFnames%set_NMLparameters(SC_NMLparameters)
    call HDFnames%set_NMLlist(SC_EBSDIndexingNameListType)
  else
! read the Modality parameter
    dataset = SC_DIModality
    call HDF%readDatasetStringArray(dataset, nlines, hdferr, stringarray)
    Modality = trim(stringarray(1))
    deallocate(stringarray)
    call Message%printMessage(' This file has the '//trim(Modality)//' modality.')
    ! select case(trim(Modality))
    !   case('EBSD')
    !   case('TKD')
    !   case('ECP')
    !   case default
    ! end select
  end if 

! make sure this is an EBSD dot product file
hdferr = HDF%openGroup(HDFnames%get_NMLfiles())

dataset = trim(HDFnames%get_NMLfilename())
  call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)

! dataset = SC_IndexEBSD
!   call H5Lexists_f(HDF%getObjectID(),trim(dataset),h_exists, hdferr)

!   if ((g_exists.eqv..FALSE.).and.(h_exists.eqv..FALSE.)) then
!       call Message%printError('readDotProductFile','this is not a dot product file')
!   end if

!   if (g_exists) then
!     call Message%printMessage(' --> EBSD dictionary indexing file found')
!   end if

!   if (h_exists) then
!     call Message%printMessage(' --> EBSD spherical indexing file found')
!   end if

call HDF%pop()

! set this value to -1 initially to trigger steps in the calling routine

DIDT%Nexp = -1

if (g_exists.eqv..TRUE.) then
!====================================
! read all NMLparameters group datasets by writing the NMLfiles string array to a
! temporary file and calling the regular nml read routine. This should replace a
! large number of individual read routines, so it simplifies the code ...
!====================================
    hdferr = HDF%openGroup(HDFnames%get_NMLfiles())
    dataset = trim(HDFnames%get_NMLfilename())
    call HDF%readdatasetstringarray(dataset, nlines, hdferr, stringarray)
    sz = shape(stringarray)
    tmpnmlname = trim(EMsoft%generateFilePath('EMtmppathname'))//'tmp.nml'
    open(unit=65,file=trim(tmpnmlname),status='unknown',form='formatted')
    write (65,"(A)") '&DIdata'
    do i=2,sz(1)
      write (65,"(A)") trim(stringarray(i))
    end do
    close(unit=65,status='keep')
    call self%readNameList(tmpnmlname)
! delete the tmp file
    open(unit=65,file=trim(tmpnmlname),status='unknown',form='formatted')
    close(unit=65,status='delete')
    call HDF%pop()
end if
!====================================
!====================================

! for .ctf files we also will need the sample tilt, which we can get from
! the 'Scan 1/EBSD/Header/Sample Tilt' data set

groupname = 'Scan 1'
    hdferr = HDF%openGroup(groupname)
if (trim(Modality).eq.'EBSD') groupname = SC_EBSD
if (trim(Modality).eq.'TKD') groupname = SC_TKD
if (trim(Modality).eq.'ECP') groupname = SC_ECP
    hdferr = HDF%openGroup(groupname)
groupname = 'Header'
    hdferr = HDF%openGroup(groupname)

dataset = 'Sample Tilt'
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetFloat(dataset, hdferr, DIDT%MCsig)
    else
      call Message%printMessage('  --> no Sample Tilt data set found ... continuing ... ')
    end if
    call HDF%pop()

! open the Scan 1/(Modality)/Data group; dictionary indexing files only have one "scan" in them...
groupname = SC_Data
    hdferr = HDF%openGroup(groupname)

dataset = 'Comment'
    DIDT%orthocomment = .FALSE. 
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
        DIDT%orthocomment = .TRUE. 
        call Message%printMessage('  --> This dot product file was modified by EMEBSDDIchangesetting !!! ')
    end if 

! integers
dataset = SC_FZcnt
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetInteger(dataset, hdferr, DIDT%FZcnt)
    else
      call Message%printMessage('  --> no FZcnt data set found ... continuing ... ')
    end if

dataset = SC_NumExptPatterns
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetInteger(dataset, hdferr, DIDT%Nexp)
    else
      call Message%printMessage('  --> no NumExptPatterns data set found ... continuing ... ')
    end if

dataset = SC_PointGroupNumber
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetInteger(dataset, hdferr, DIDT%pgnum)
    else
      call Message%printMessage('  --> no PointGroupNumber data set found ... continuing ... ')
    end if

! various optional arrays
if (present(getADP)) then
  if (getADP.eqv..TRUE.) then
   dataset = SC_AvDotProductMap
   call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
   if (g_exists.eqv..TRUE.) then
    if (sum(ebsdnl%ROI).ne.0) then 
      allocate(DIDT%ADP(ebsdnl%ROI(3), ebsdnl%ROI(4)))
      call HDF%read2DImage(dataset, DIDT%ADP, ebsdnl%ROI(3), ebsdnl%ROI(4))
    else
      allocate(DIDT%ADP(ebsdnl%ipf_wd, ebsdnl%ipf_ht))
      call HDF%read2DImage(dataset, DIDT%ADP, ebsdnl%ipf_wd, ebsdnl%ipf_ht)
    end if 
   else
        call Message%printMessage('  --> no AvDotProductMap data set found ... continuing ... ')
   end if
  end if
end if

if (present(getAverageOrientations)) then
  if (getAverageOrientations.eqv..TRUE.) then
    dataset = SC_AverageOrientations
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetFloatArray(dataset, dims2, hdferr, DIDT%AverageOrientations)
    else
      call Message%printMessage('  --> no AverageOrientations data set found ... continuing ... ')
    end if
  end if
end if

if (present(getCI)) then
  if (getCI.eqv..TRUE.) then
! if this is a Spherical Indexing file, then we should look for the 'Metric' data set,
! otherwise the CI data set
    dataset = SC_CI
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetFloatArray(dataset, dims, hdferr, DIDT%CI)
    else
      call Message%printMessage('  --> no CI data set found ... continuing ... ')
      dataset = 'Metric'
      call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
      if (g_exists.eqv..TRUE.) then
        call HDF%readDatasetFloatArray(dataset, dims, hdferr, DIDT%CI)
      else
        call Message%printMessage('  --> no Metric data set found ... continuing ... ')
      end if
    end if
  end if
end if

if (present(getEulerAngles)) then
  if (getEulerAngles.eqv..TRUE.) then
    dataset = SC_EulerAngles
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetFloatArray(dataset, dims2, hdferr, DIDT%EulerAngles)
   else
      call Message%printMessage('  --> no EulerAngles data set found ... continuing ... ')
   end if
  end if
end if

if (present(getDictionaryEulerAngles)) then
  if (getDictionaryEulerAngles.eqv..TRUE.) then
    dataset = 'DictionaryEulerAngles'
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetFloatArray(dataset, dims2, hdferr, DIDT%DictionaryEulerAngles)
    else
      call Message%printMessage('  --> no DictionaryEulerAngles data set found ... continuing ... ')
    end if
  end if
end if

if (present(getFit)) then
  if (getFit.eqv..TRUE.) then
    dataset = SC_Fit
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetFloatArray(dataset, dims, hdferr, DIDT%Fit)
    else
      call Message%printMessage('  --> no Fit data set found ... continuing ... ')
    end if
  end if
end if

if (present(getIQ)) then
  if (getIQ.eqv..TRUE.) then
    dataset = SC_IQ
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetFloatArray(dataset, dims, hdferr, DIDT%IQ)
    else
      call Message%printMessage('  --> no IQ data set found ... continuing ... ')
    end if
  end if
end if

if (present(getKAM)) then
  if (getKAM.eqv..TRUE.) then
    dataset = SC_KAM
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetFloatArray(dataset, dims2, hdferr, DIDT%KAM)
    else
      call Message%printMessage('  --> no KAM data set found ... continuing ... ')
    end if
  end if
end if

if (present(getOSM)) then
  if (getOSM.eqv..TRUE.) then
    dataset = SC_OSM
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetFloatArray(dataset, dims2, hdferr, DIDT%OSM)
    else
      call Message%printMessage('  --> no OSM data set found ... continuing ... ')
    end if
  end if
end if

if (present(getPhase)) then   ! this is a 1-byte integer, to be implemented
  if (getPhase.eqv..TRUE.) then
!   dataset = SC_Phase
!   call HDF%readDatasetIntegerArray(dataset, dims, hdferr, EBSDDIdata%Phase)
    call Message%printMessage('Phase','reading the Phase variable is not yet implemented')
  end if
end if

if (present(getPhi)) then
  if (getPhi.eqv..TRUE.) then
    dataset = SC_Phi
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetFloatArray(dataset, dims, hdferr, DIDT%Phi)
    else
      call Message%printMessage('  --> no Phi data set found ... continuing ... ')
    end if
  end if
end if

if (present(getPhi1)) then
  if (getPhi1.eqv..TRUE.) then
    dataset = SC_Phi1
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetFloatArray(dataset, dims, hdferr, DIDT%Phi1)
    else
      call Message%printMessage('  --> no Phi1 data set found ... continuing ... ')
    end if
  end if
end if

if (present(getPhi2)) then
  if (getPhi2.eqv..TRUE.) then
    dataset = SC_Phi2
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetFloatArray(dataset, dims, hdferr, DIDT%Phi2)
    else
      call Message%printMessage('  --> no Phi2 data set found ... continuing ... ')
    end if
  end if
end if

if (present(getSEMsignal)) then
  if (getSEMsignal.eqv..TRUE.) then
    dataset = SC_SEMsignal
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetIntegerArray(dataset, dims, hdferr, DIDT%SEMsignal)
    else
      call Message%printMessage('  --> no SEMsignal data set found ... continuing ... ')
    end if
  end if
end if

if (present(getTopDotProductList)) then
  if (getTopDotProductList.eqv..TRUE.) then
    dataset = SC_TopDotProductList
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetFloatArray(dataset, dims2, hdferr, DIDT%TopDotProductList)
    else
      call Message%printMessage('  --> no TopDotProductList data set found ... continuing ... ')
    end if
  end if
end if

if (present(getTopMatchIndices)) then
  if (getTopMatchIndices.eqv..TRUE.) then
    dataset = SC_TopMatchIndices
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetIntegerArray(dataset, dims2, hdferr, DIDT%TopMatchIndices)
    else
      call Message%printMessage('  --> no TopMatchIndices data set found ... continuing ... ')
    end if
  end if
end if

if (present(getValid)) then
  if (getValid.eqv..TRUE.) then
!   dataset = SC_Valid
!   call HDF%readDatasetIntegerArray(dataset, dims, hdferr, EBSDDIdata%Valid)
    call Message%printMessage('Valid','reading the Valid variable is not yet implemented')
  end if
end if

if (present(getXPosition)) then
  if (getXPosition.eqv..TRUE.) then
    dataset = SC_XPosition
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetFloatArray(dataset, dims, hdferr, DIDT%XPosition)
    else
      call Message%printMessage('  --> no XPosition data set found ... continuing ... ')
    end if
  end if
end if

if (present(getYPosition)) then
  if (getYPosition.eqv..TRUE.) then
    dataset = SC_YPosition
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetFloatArray(dataset, dims, hdferr, DIDT%YPosition)
    else
      call Message%printMessage('  --> no YPosition data set found ... continuing ... ')
    end if
  end if
end if

if (present(getRefinedDotProducts)) then
  if (getRefinedDotProducts.eqv..TRUE.) then
    dataset = SC_RefinedDotProducts
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then
      call HDF%readDatasetFloatArray(dataset, dims, hdferr, DIDT%RefinedDotProducts)
    else
      call Message%printMessage('readDotProductFile','There is no RefinedDotProducts data set in this file')
    end if
  end if
end if

if (present(getRefinedEulerAngles)) then
  if (getRefinedEulerAngles.eqv..TRUE.) then
    dataset = SC_RefinedEulerAngles
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then
      call HDF%readDatasetFloatArray(dataset, dims2, hdferr, DIDT%RefinedEulerAngles)
    else
      call Message%printMessage('readDotProductFile','There is no RefinedEulerAngles data set in this file')
    end if
  end if
end if

call HDF%pop()

groupname = SC_Header
    hdferr = HDF%openGroup(groupname)

dataset = SC_StepX
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then
      call HDF%readDatasetFloat(dataset, hdferr, ebsdnl%StepX)
    else
      call Message%printMessage('readDotProductFile','There is no StepX data set in this file')
    end if

dataset = SC_StepY
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then
      call HDF%readDatasetFloat(dataset, hdferr, ebsdnl%StepY)
    else
      call Message%printMessage('readDotProductFile','There is no StepY data set in this file')
    end if

! and close the HDF5 dot product file
call HDF%pop(.TRUE.)

end associate

if (DIModality.eq.'old') HDFnames = saveHDFnames

end subroutine readDotProductFile_


!--------------------------------------------------------------------------
subroutine h5_writeFile_(self, EMsoft, HDF, HDFnames, vendor, mcnl, xtalname, dstr, tstrb, tstre, ipar, resultmain, exptIQ, &
                         indexmain, dicteulerarray, dpmap, progname, nmldeffile, OSMmap)
!DEC$ ATTRIBUTES DLLEXPORT :: h5_writeFile_
!! author: MDG
!! version: 1.0
!! date: 04/03/20
!!
!! perform the computations

use mod_EMsoft
use mod_HDFnames
use mod_io
use mod_MCfiles
use mod_DIsupport

IMPLICIT NONE

class(DIfile_T), INTENT(INOUT)                      :: self
type(EMsoft_T), INTENT(INOUT)                       :: EMsoft
type(HDF_T), INTENT(INOUT)                          :: HDF
type(HDFnames_T), INTENT(INOUT)                     :: HDFnames
character(3),INTENT(IN)                             :: vendor   ! 'TSL' 'HKL' 'BRU'
! type(DictionaryIndexingNameListType),INTENT(INOUT)  :: ebsdnl
!f2py intent(in,out) ::  ebsdnl
type(MCOpenCLNameListType),INTENT(INOUT)            :: mcnl
character(fnlen),INTENT(IN)                         :: xtalname
character(11),INTENT(INOUT)                         :: dstr
!f2py intent(in,out) ::  dstr
character(15),INTENT(IN)                            :: tstrb
character(15),INTENT(IN)                            :: tstre
integer(kind=irg),INTENT(INOUT)                     :: ipar(10)
!f2py intent(in,out) ::  ipar
real(kind=sgl),INTENT(IN)                           :: resultmain(ipar(1),ipar(2))
real(kind=sgl),INTENT(IN)                           :: exptIQ(ipar(3))
integer(kind=irg),INTENT(IN)                        :: indexmain(ipar(1),ipar(2))
real(kind=sgl),INTENT(INOUT)                        :: dicteulerarray(3,ipar(4))
!f2py intent(in,out) ::  dicteulerarray
real(kind=sgl),INTENT(IN)                           :: dpmap(ipar(3))
character(fnlen),INTENT(IN)                         :: progname
character(fnlen),INTENT(IN)                         :: nmldeffile
real(kind=sgl),INTENT(OUT)                          :: OSMmap(ipar(7),ipar(8))

type(IO_T)                                          :: Message

character(fnlen, KIND=c_char),allocatable,TARGET    :: stringarray(:)
integer(kind=irg)                                   :: hdferr, filetype, i, j, ii, jj,indx, istat, ipar2(6), L
character(fnlen)                                    :: groupname, dataset, h5ebsdfile, savefile, Modality
logical                                             :: noindex, g_exists, overwrite=.TRUE.
real(kind=sgl)                                      :: eulerarray(3,ipar(4)), WD

real(kind=sgl),allocatable                          :: kam(:,:), ISMap(:)

real(kind=sgl),allocatable                          :: exptCI(:), eangle(:), eangles(:,:), results(:), avEuler(:,:), &
                                                       lresultmain(:,:), eulers(:,:)
integer(kind=1),allocatable                         :: iPhase(:), valid(:)
integer(kind=irg),allocatable                       :: SEMsignal(:), lindexmain(:,:)
real(kind=sgl)                                      :: isratio, io_real(1)

associate(ebsdnl=>self%nml)

! copy the dictionary euler angle array
eulerarray = dicteulerarray

!=====================================================
! write the output in the format of an h5ebsd file
!!!! THIS PART IS STILL UNDER DEVELOPMENT !!!!
! we use the TSL h5ebsd file as a template for now; this 
! can be extended later for other vendor formats
!=====================================================

if (vendor.ne.'TSL') then
  call Message%printMessage(' Only TSL h5ebsd file format is implemented in this version.')
  call Message%printMessage(' Program results will be saved in this format.')
end if

allocate(stringarray(1))

! Create a new file using the default properties.
  hdferr =  HDF%createFile(self%getfilename())
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error opening file', hdferr)

  filetype = 1  ! for vendor = 'TSL'; others not yet implemented
  call self%h5_writeInfo(EMsoft, HDF, HDFnames, filetype, dstr, tstrb, tstre, progname, ebsdnl, nmldeffile)

! here we start with the h5ebsd-specific stuff
  groupname = 'Scan 1'
  hdferr = HDF%createGroup(groupname)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error opening group Scan 1', hdferr)

! we need to name the following group according to the scattering modality 
Modality = self%getModality()
if (trim(Modality).eq.'EBSD') groupname = SC_EBSD
if (trim(Modality).eq.'TKD') groupname = SC_TKD
if (trim(Modality).eq.'ECP') groupname = SC_ECP
  hdferr = HDF%createGroup(groupname)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error opening group EBSD/TKD/ECP', hdferr)

!=====================================================
!=====================================================
groupname = SC_Data
  hdferr = HDF%createGroup(groupname)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error opening group Data', hdferr)

! there are 15 datasets in this Data group: CI, Fit, IQ, PRIAS Bottom Strip,
! PRIAS Center Square, PRIAS Top Strip, Pattern, Phase, Phi, Phi1, Phi2,
! SEM Signal, Valid, X Position, Y Position

!=====================================================
! CI Confidence Index: real(kind=sgl), one for each pattern... we take this
! to be the largest dot product
dataset = SC_CI
  allocate(exptCI(ipar(3)))
  exptCI = resultmain(1,1:ipar(3))
  hdferr = HDF%writeDatasetFloatArray(dataset, exptCI, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset CI', hdferr)

! we also insert a visual map of the Confidence Index, resampled on a rectangular array
dataset = SC_CIMap
  call h5_write2DImageFromVector(HDF, dataset, exptCI, ipar(3), ebsdnl)
  deallocate(exptCI)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset CImap', hdferr)

!=====================================================
! Fit
dataset = SC_Fit
  allocate(eangle(ipar(3)),stat=istat)
  eangle = 1.0
  hdferr = HDF%writeDatasetFloatArray(dataset, eangle, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Fit', hdferr)
  deallocate(eangle)

!=====================================================
! this option is disabled starting in version 4.3 [05/19/19, MDG]; code block can be deleted
! Averaged Orientation Map (using near-match list and quaternion logarithm averaging)
! define the ipar2 entries
!   allocate(avEuler(3,ipar(3)))
  ipar2(1) = ipar(6)
  ipar2(2) = ipar(5)
  ipar2(3) = ipar(3)
  ipar2(4) = ipar(1)
  ipar2(5) = ipar(2)
  ipar2(6) = ebsdnl%nnav
! ! get the avEuler array
!   eulerarray = eulerarray * sngl(cPi)/180.0
!   call EBSDgetAverageOrientations(ipar2, eulerarray, indexmain(1:ipar2(4),1:ipar2(5)), resultmain(1:ipar2(4),1:ipar2(5)), &
!                                   avEuler)
!   eulerarray = eulerarray * 180.0/sngl(cPi)

! ! and write it to the HDF file
! dataset = SC_AverageOrientations
!   hdferr = HDF_writeDatasetFloatArray2D(dataset, avEuler, 3, ipar(3))
!   if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset AverageOrientations')

! get the nearest neighbor Kernel Average Misorientation Map (KAM)
  allocate(eulers(3,ipar(3)))
  do i=1,ipar(3)
    eulers(1:3,i) = eulerarray(1:3,indexmain(1,i))
  end do
  eulers = eulers*dtor
dataset = SC_KAM
  if (sum(ebsdnl%ROI).ne.0) then
    allocate(kam(ebsdnl%ROI(3),ebsdnl%ROI(4)))
    call getKAMMap(ipar(3), eulers, ebsdnl%ROI(3), ebsdnl%ROI(4), ipar(6), kam)
    kam = kam*rtod
    hdferr = HDF%writeDatasetFloatArray(dataset, kam, ebsdnl%ROI(3), ebsdnl%ROI(4))
  else
    allocate(kam(ebsdnl%ipf_wd,ebsdnl%ipf_ht))
    call getKAMMap(ipar(3), eulers, ebsdnl%ipf_wd, ebsdnl%ipf_ht, ipar(6), kam)
    kam = kam*rtod
    hdferr = HDF%writeDatasetFloatArray(dataset, kam, ebsdnl%ipf_wd, ebsdnl%ipf_ht)
  end if
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset KAM', hdferr)
  deallocate(kam, eulers)

! get the Orientation Similarity Map (OSM); map is now returned to calling routine [MDG, 3/5/18]
dataset = SC_OSM
  if (sum(ebsdnl%ROI).ne.0) then
!   allocate(osm(ebsdnl%ROI(3),ebsdnl%ROI(4)))
    call getOrientationSimilarityMap( (/ipar(1), ipar(2)/), indexmain, ebsdnl%nosm, ebsdnl%ROI(3), ebsdnl%ROI(4), OSMmap)
    hdferr = HDF%writeDatasetFloatArray(dataset, OSMmap, ebsdnl%ROI(3), ebsdnl%ROI(4))
  else
!   allocate(osm(ebsdnl%ipf_wd,ebsdnl%ipf_ht))
    call getOrientationSimilarityMap( (/ipar(1), ipar(2)/), indexmain, ebsdnl%nosm, ebsdnl%ipf_wd, ebsdnl%ipf_ht, OSMmap)
    hdferr = HDF%writeDatasetFloatArray(dataset, OSMmap, ebsdnl%ipf_wd, ebsdnl%ipf_ht)
  end if
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset OSM', hdferr)

! we also insert a visual map of the Confidence Index, resampled on a rectangular array
! dataset = 'FitMap'
! call h5ebsd_write2DImageFromVector(dataset, totnumexpt, exptCI, ebsdnl)

!=====================================================
! IQ Image Quality; computed using the second moment of the pattern power spectrum
dataset = SC_IQ
  hdferr = HDF%writeDatasetFloatArray(dataset, exptIQ, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset IQ', hdferr)

! we also insert a visual map of the Image Quality, resampled on a rectangular array
dataset = SC_IQMap
  call h5_write2DImageFromVector(HDF, dataset, exptIQ, ipar(3), ebsdnl)

!=====================================================
! generate the indexing success map (ISM)
  eulerarray = eulerarray * dtor
  allocate(ISMap(ipar(7)*ipar(8)))
  call getIndexingSuccessMap(ipar, indexmain, eulerarray, ebsdnl%nism, ebsdnl%nnk, ebsdnl%nthreads, ISMap)
  eulerarray = eulerarray * rtod

dataset = SC_ISM
  hdferr = HDF%writeDatasetFloatArray(dataset, ISMap, ipar(7)*ipar(8))

dataset = SC_ISMap
  call h5_write2DImageFromVector(HDF, dataset, ISMap, ipar(7)*ipar(8), ebsdnl, binary=ebsdnl%isangle)
  j = 0
  do i=1,ipar(7)*ipar(8)
    if (ISMap(i).le.ebsdnl%isangle) j = j+1
  end do
  isratio = 100.0 * real(j) / real(ipar(7)*ipar(8))
  io_real(1) = isratio
  call Message%WriteValue(' Indexing Success Rate (%) : ',io_real,1)
  deallocate(ISMap)

dataset = SC_ISR
  hdferr = HDF%writeDatasetFloat(dataset, isratio)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset ISR', hdferr)


!=====================================================
! PRIAS Bottom Strip: to be implemented
!   call Message('h5ebsd_writeFile: writing of ->PRIAS Bottom Strip<- data not yet implemented.')

!=====================================================
! PRIAS Center Square: to be implemented
!   call Message('h5ebsd_writeFile: writing of ->PRIAS Center Strip<- data not yet implemented.')

!=====================================================
! PRIAS Top Strip: to be implemented
!   call Message('h5ebsd_writeFile: writing of ->PRIAS Top Strip<- data not yet implemented.')

!=====================================================
! Pattern: in principle, this is where the fitted patterns could be stored
! This will require re-computing them for the best match orientations; we
! could leave this as an option for the user, to be implemented.
!   call Message('h5ebsd_writeFile: writing of ->Pattern<- data not yet implemented.')

!=====================================================
! Phase: Phase identifier (all zero for now)
dataset = SC_Phase
  allocate(iPhase(ipar(3)),stat=istat)
  iPhase = 0
  hdferr = HDF%writeDatasetInteger1byteArray1D(dataset, iPhase, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Phase', hdferr)
  deallocate(iPhase)

!=====================================================
! SEM Signal: all 0 for now
dataset = SC_SEMSignal
  allocate(SEMsignal(ipar(3)),stat=istat)
  SEMsignal = 10000
  hdferr = HDF%writeDatasetIntegerArray(dataset, SEMsignal, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset SEM Signal', hdferr)
  deallocate(SEMsignal)

!=====================================================
! Valid : all 0 for now
dataset = SC_Valid
  allocate(valid(ipar(3)),stat=istat)
  valid = 0
  hdferr = HDF%writeDatasetInteger1byteArray1D(dataset, valid, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Valid', hdferr)
  deallocate(valid)

dataset = SC_EulerAngles
  call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
  allocate(eangles(3,ipar(3)),stat=istat)
  do ii = 1,ipar(3)
    indx = indexmain(1,ii)
    eangles(1:3,ii) = eulerarray(1:3,indx)
  end do
  eangles = eangles * dtor
  if (g_exists) then
     hdferr = HDF%writeDatasetFloatArray(dataset, eangles, 3, ipar(3), overwrite)
  else
     hdferr = HDF%writeDatasetFloatArray(dataset, eangles, 3, ipar(3))
  end if
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset EulerAngles', hdferr)

!=====================================================
! Euler angles: Phi
dataset = SC_Phi
  allocate(eangle(ipar(3)),stat=istat)
  do ii = 1,ipar(3)
    indx = indexmain(1,ii)
    eangle(ii) = eulerarray(2,indx)
  end do
  eangle = eangle * dtor
  hdferr = HDF%writeDatasetFloatArray(dataset, eangle, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Phi', hdferr)

!=====================================================
! Euler angles: Phi1
dataset = SC_Phi1
  do ii = 1,ipar(3)
    indx = indexmain(1,ii)
    eangle(ii) = eulerarray(1,indx)
  end do
  eangle = eangle * dtor
  hdferr = HDF%writeDatasetFloatArray(dataset, eangle, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Phi1', hdferr)

!=====================================================
! Euler angles: Phi2
dataset = SC_Phi2
  do ii = 1,ipar(3)
    indx = indexmain(1,ii)
    eangle(ii) = eulerarray(3,indx)
  end do
  eangle = eangle * dtor
  hdferr = HDF%writeDatasetFloatArray(dataset, eangle, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Phi2', hdferr)
  deallocate(eangle)

!=====================================================
! X Position: list of x positions for sampling points; requires knowledge of step size
! from Header
dataset = SC_XPos
  allocate(results(ipar(3)),stat=istat)
  if (sum(ebsdnl%ROI).eq.0) then
    do jj=1,ebsdnl%ipf_ht
      do ii=1,ebsdnl%ipf_wd
        results(ebsdnl%ipf_wd*(jj-1)+ii) = (ii-1)*ebsdnl%StepX
      end do
    end do
  else
    do jj=1,ebsdnl%ROI(4)
      do ii=1,ebsdnl%ROI(3)
        results(ebsdnl%ROI(3)*(jj-1)+ii) = (ii-1)*ebsdnl%StepX
      end do
    end do
  end if
  hdferr = HDF%writeDatasetFloatArray(dataset, results, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset X Position', hdferr)


!=====================================================
! Y Position: list of y positions for sampling points; requires knowledge of step size
! from Header
dataset = SC_YPos
  if (sum(ebsdnl%ROI).eq.0) then
    do jj=1,ebsdnl%ipf_ht
      do ii=1,ebsdnl%ipf_wd
        results(ebsdnl%ipf_wd*(jj-1)+ii) = (ii-1)*ebsdnl%StepY
      end do
    end do
  else
    do jj=1,ebsdnl%ROI(4)
      do ii=1,ebsdnl%ROI(3)
        results(ebsdnl%ROI(3)*(jj-1)+ii) = (ii-1)*ebsdnl%StepY
      end do
    end do
  end if
  hdferr = HDF%writeDatasetFloatArray(dataset, results, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Y position', hdferr)
  deallocate(results)

!=====================================================
!=====================================================
! this concludes the standard data sets in the Data group
! here, we have additional data sets based on results from the
! dictionary indexing program; these are not part of the standard
! TSL h5ebsd file format.
!=====================================================
! EBSD average dot product map
dataset = SC_AvDotProductMap
  call h5_write2DImageFromVector(HDF, dataset, dpmap, ipar(3), ebsdnl)

! number of samples in dictionary
dataset = SC_FZcnt
  hdferr = HDF%writeDatasetInteger(dataset, ipar(5))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset FZcnt', hdferr)

! point group number
dataset = SC_PointGroupNumber
  hdferr = HDF%writeDatasetInteger(dataset, ipar(6))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset PointGroupNumber', hdferr)

! Ncubochoric
dataset = SC_Ncubochoric
  hdferr = HDF%writeDatasetInteger(dataset, ebsdnl%ncubochoric)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Ncubochoric', hdferr)

! write the list of sampled Euler angles
dataset = SC_DictionaryEulerAngles
  hdferr = HDF%writeDatasetFloatArray(dataset, dicteulerarray, 3, ipar(4))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset EulerAngles', hdferr)

! number of experimental patterns
dataset = SC_NumExptPatterns
  hdferr = HDF%writeDatasetInteger(dataset, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset NumExptPatterns', hdferr)

! list of top nnk dot product values
dataset = SC_TopDotProductList
  hdferr = HDF%writeDatasetFloatArray(dataset, resultmain, ipar(1), ipar(2))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset TopDotProductList', hdferr)

! indices of top matches into Euler angle list
dataset = SC_TopMatchIndices
  hdferr = HDF%writeDatasetIntegerArray(dataset, indexmain, ipar(1), ipar(2))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset TopMatchIndices', hdferr)

! leave this group
  call HDF%pop()
!=====================================================
!=====================================================

!=====================================================
!=====================================================
! create the Header group
groupname = SC_Header
  hdferr = HDF%createGroup(groupname)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error creating group Header', hdferr)

! there are 15 datasets in this group: Camera Azimuth Angle, Camera Elevation Angle,
! Grid Type, Notes, Operator, Pattern Height, Pattern Width, Sample ID, Sample Tilt,
! Scan ID, Step X, Step Y, Working Distance, nColumns, nRows
! there are also 3 groups: Coordinate System, Pattern Center Calibration, and Phase

!=====================================================
! Camera Azimuthal Angle
dataset = SC_CameraAzimuthalAngle
  hdferr = HDF%writeDatasetFloat(dataset, 0.0)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Camera Azimuthal Angle', hdferr)

!=====================================================
! Camera Elevation Angle
dataset = SC_CameraElevationAngle
  hdferr = HDF%writeDatasetFloat(dataset, ebsdnl%thetac)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset TopMatchIndices', hdferr)

!=====================================================
! Coordinate System group
  call h5_writeCoordinateSystemGroup(EMsoft, HDF)

!=====================================================
! Grid Type
dataset = SC_GridType
  stringarray(1) = 'SqrGrid'
  hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Grid Type', hdferr)

!=====================================================
! Notes
dataset = SC_Notes
  stringarray(1) = ''
  hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Notes', hdferr)

!=====================================================
! Operator
dataset = SC_Operator
  stringarray(1) = trim(EMsoft%getConfigParameter('Username'))//' ['// &
                   trim(EMsoft%getConfigParameter('Useremail'))//']'
  hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Operator', hdferr)

!=====================================================
! Pattern Center Calibration group
  call h5_writePatternCenterGroup(HDF, ebsdnl%xpc, ebsdnl%ypc, ebsdnl%L, ebsdnl%delta, (/ebsdnl%numsx, ebsdnl%numsy/))

!=====================================================
! Pattern height
dataset = SC_PatternHeight
  hdferr = HDF%writeDatasetInteger(dataset, ebsdnl%numsx/ebsdnl%binning)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Pattern Height', hdferr)

!=====================================================
! Pattern width
dataset = SC_PatternWidth
  hdferr = HDF%writeDatasetInteger(dataset, ebsdnl%numsy/ebsdnl%binning)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Pattern Width', hdferr)

!=====================================================
! Phase group
groupname = SC_Phase
  hdferr = HDF%createGroup(groupname)
groupname = "1"

  call h5_writePhaseGroup(EMsoft, HDF, groupname, xtalname)

! close the Phase group
  call HDF%pop()

!=====================================================
! Sample ID
dataset = SC_SampleID
  stringarray(1) = ''
  hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

!=====================================================
! Sample Tilt
dataset = SC_SampleTilt
  hdferr = HDF%writeDatasetFloat(dataset, sngl(mcnl%sig))

!=====================================================
! Scan ID
dataset = SC_ScanID
  stringarray(1) = ''
  hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

!=====================================================
! Step X
dataset = SC_StepX
  hdferr = HDF%writeDatasetFloat(dataset, ebsdnl%StepX)

!=====================================================
! Step Y
dataset = SC_StepY
  hdferr = HDF%writeDatasetFloat(dataset, ebsdnl%StepY)

!=====================================================
! Working Distance
dataset = SC_WorkingDistance
  WD = 0.0
  hdferr = HDF%writeDatasetFloat(dataset, WD)

!=====================================================
! nColumns
dataset = SC_nColumns
  hdferr = HDF%writeDatasetInteger(dataset, ebsdnl%ipf_wd)

!=====================================================
! nRows
dataset = SC_nRows
  hdferr = HDF%writeDatasetInteger(dataset, ebsdnl%ipf_ht)

!=====================================================
!=====================================================

! once all these have been written, we simply pop all the way to the top and close the file
  call HDF%pop(.TRUE.)

end associate

end subroutine h5_writeFile_

!--------------------------------------------------------------------------
subroutine h5_writeInfo(self, EMsoft, HDF, HDFnames, filetype, dstr, tstrb, tstre, progname, ebsdnl, nmldeffile)
!DEC$ ATTRIBUTES DLLEXPORT :: h5_writeInfo
!! author: MDG
!! version: 1.0
!! date: 04/03/20
!!
!! write general information fields to the h5ebsd file, including EMsoft specific fields

use mod_EMsoft
use mod_HDFnames

IMPLICIT NONE

class(DIfile_T),INTENT(INOUT)                       :: self
type(EMsoft_T),INTENT(INOUT)                        :: EMsoft
type(HDF_T),INTENT(INOUT)                           :: HDF
type(HDFnames_T),INTENT(INOUT)                      :: HDFnames
integer(kind=irg),INTENT(IN)                        :: filetype
character(11),INTENT(INOUT)                         :: dstr
!f2py intent(in,out) ::  dstr
character(15),INTENT(IN)                            :: tstrb
character(15),INTENT(IN)                            :: tstre
character(fnlen),INTENT(IN)                         :: progname
type(DictionaryIndexingNameListType),INTENT(INOUT)  :: ebsdnl
!f2py intent(in,out) ::  ebsdnl
character(fnlen),INTENT(IN)                         :: nmldeffile

character(fnlen, KIND=c_char),allocatable,TARGET    :: stringarray(:)
character(fnlen)                                    :: groupname, dataset, nmlname, manufacturer
integer(kind=irg)                                   :: hdferr

! to be replaced by HDFnames code
if (filetype.eq.1) then ! EBSDDictionarIndexing file
  manufacturer = 'EMDI.f90'
  nmlname = 'DictionaryIndexingNML'
else
  manufacturer = ''
  nmlname = ''
end if

allocate(stringarray(1))

! set the Manufacturer and Version data sets
dataset = SC_Manufacturer
  stringarray(1)= trim(manufacturer)
  hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

dataset = SC_Version
  stringarray(1)= 'EMsoft '//EMsoft%getConfigParameter('EMsoftversion')
  hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

! add the EMsoft header group
! write the EMheader to the file
groupname = SC_h5EBSD
  call HDF%writeEMheader(EMsoft, dstr, tstrb, tstre, progname)

! create a namelist group to write all the namelist files into
groupname = SC_NMLfiles
  hdferr = HDF%createGroup(groupname)

! and write the nml file for this program to the HDF5 file
! read the text file and write the array to the file
  dataset = trim(nmlname)
  hdferr = HDF%writeDatasetTextFile(dataset, nmldeffile)

! leave this group
  call HDF%pop()

! create a namelist group to write all the namelist files into
groupname = SC_NMLparameters
  hdferr = HDF%createGroup(groupname)
  if (filetype.eq.1) then
    call self%writeHDFNameList(HDF, HDFnames, ebsdnl)
  end if

! leave this group
  call HDF%pop()

end subroutine h5_writeInfo

!--------------------------------------------------------------------------
subroutine h5_write2DImageFromVector(HDF, dataset, inpvec, nump, ebsdnl, binary)
!DEC$ ATTRIBUTES DLLEXPORT :: h5_write2DImageFromVector
!! author: MDG
!! version: 1.0
!! date: 04/03/20
!!
!! write a gray scale image to the HDF5 file starting from a 1D vector

use mod_io

IMPLICIT NONE

type(HDF_T),INTENT(INOUT)               :: HDF
character(fnlen),INTENT(IN)             :: dataset
integer(kind=irg),INTENT(IN)            :: nump
real(kind=sgl),INTENT(IN)               :: inpvec(nump)
type(DictionaryIndexingNameListType),INTENT(IN)  :: ebsdnl
real(kind=sgl),OPTIONAL,INTENT(IN)      :: binary

type(IO_T)                              :: Message
real(kind=sgl)                          :: mi, ma
integer(kind=irg)                       :: istat, ii, jj, hdferr
real(kind=sgl),allocatable              :: newvec(:)
integer(kind=irg),allocatable           :: image(:,:)
integer(HSIZE_T)                        :: width, height
logical                                 :: isbinary

isbinary = .FALSE.
if (present(binary)) isbinary=.TRUE.

allocate(newvec(nump),stat=istat)
if (istat.ne.0) call Message%printError('h5_write2DImageFromVector','Could not allocate array for copy of input image')

newvec = inpvec

if (sum(ebsdnl%ROI).ne.0) then
  width = ebsdnl%ROI(3)
  height = ebsdnl%ROI(4)
else
  width = ebsdnl%ipf_wd
  height = ebsdnl%ipf_ht
end if
allocate(image(width,height),stat=istat)
if (istat.ne.0) call Message%printError('h5_write2DImageFromVector','Could not allocate array for output image')

if (isbinary.eqv..TRUE.) then
  do jj = 1,height
    do ii = 1, width
      if (newvec((jj-1)*width+ii).gt.ebsdnl%isangle) then
        image(ii,jj) = 0
      else
        image(ii,jj) = 255
      end if
    end do
  end do
else
  mi = minval(newvec)
  newvec = newvec - mi
  ma = maxval(newvec)

  do jj = 1,height
    image(1:width,jj) = int(255.0*newvec((jj-1)*width+1:jj*width)/ma)
  end do
end if

call h5immake_image_8bit_f(HDF%getObjectID(),dataset,width,height,image,hdferr)
deallocate(image, newvec)

end subroutine h5_write2DImageFromVector

!--------------------------------------------------------------------------
subroutine h5_writeCoordinateSystemGroup(EMsoft, HDF)
!DEC$ ATTRIBUTES DLLEXPORT :: h5_writeCoordinateSystemGroup
!! author: MDG
!! version: 1.0
!! date: 02/17/20
!!
!! write information about the sample/detector coordinate frames

use mod_EMsoft
use mod_io
use h5im

IMPLICIT NONE

type(EMsoft_T),INTENT(INOUT)       :: EMsoft
type(HDF_T),INTENT(INOUT)          :: HDF

type(IO_T)                         :: Message
character(fnlen)                   :: groupname, dataset, fname, resourcepathname
integer(kind=irg)                  :: hdferr
integer(HSIZE_T),allocatable       :: EBSDview(:,:,:), schematic(:,:,:)
character(1),allocatable           :: chararr(:,:,:)
integer(kind=irg)                  :: dims(3), istat
integer(HSIZE_T)                   :: width, height

! create the Coordinate System group
groupname = 'Coordinate System'
hdferr = HDF%createGroup(groupname)

!=====================================================
! EBSD View Reference Frame
fname = trim(EMsoft%generateFilePath('Resourcepathname'))//'EBSDview.data'
open(unit=50,file=trim(fname),status='old',form='unformatted')
read(50) dims
allocate(EBSDview(dims(1),dims(2),dims(3)),chararr(dims(1),dims(2),dims(3)),stat=istat)
if (istat.ne.0) call Message%printError('h5_writeCoordinateSystemGroup','Could not allocate array for EBSD view output image')
read(50) chararr
close(unit=50,status='keep')
EBSDview = ichar(chararr)

dataset = 'EBSD View Reference Frame'
width = dims(2)
height = dims(3)
call h5immake_image_24bit_f(HDF%getObjectID(),dataset,width,height,'INTERLACE_PIXEL',int(EBSDview),hdferr)
deallocate(EBSDview,chararr)

!=====================================================
! Schematic 1
fname = trim(EMsoft%generateFilePath('Resourcepathname'))//'Schematic1.data'
open(unit=50,file=fname,status='old',form='unformatted')
read(50) dims
allocate(schematic(dims(1),dims(2),dims(3)),chararr(dims(1),dims(2),dims(3)),stat=istat)
if (istat.ne.0) call Message%printError('h5_writeCoordinateSystemGroup','Could not allocate array for Schematic output image')
read(50) chararr
close(unit=50,status='keep')
schematic = ichar(chararr)

dataset = 'Schematic 1'
width = dims(2)
height = dims(3)
call h5immake_image_24bit_f(HDF%getObjectID(),dataset,width,height,'INTERLACE_PIXEL',int(schematic),hdferr)

!=====================================================
! Schematic 2
fname = trim(EMsoft%generateFilePath('Resourcepathname'))//'Schematic2.data'
open(unit=50,file=fname,status='old',form='unformatted')
read(50) dims
read(50) chararr
close(unit=50,status='keep')
schematic = ichar(chararr)

dataset = 'Schematic 2'
width = dims(2)
height = dims(3)
call h5immake_image_24bit_f(HDF%getObjectID(),dataset,width,height,'INTERLACE_PIXEL',int(schematic),hdferr)

!=====================================================
! Schematic 3
fname = trim(EMsoft%generateFilePath('Resourcepathname'))//'Schematic3.data'
open(unit=50,file=fname,status='old',form='unformatted')
read(50) dims
read(50) chararr
close(unit=50,status='keep')
schematic = ichar(chararr)

dataset = 'Schematic 3'
width = dims(2)
height = dims(3)
call h5immake_image_24bit_f(HDF%getObjectID(),dataset,width,height,'INTERLACE_PIXEL',int(schematic),hdferr)

!=====================================================
! Schematic 4
fname = trim(EMsoft%generateFilePath('Resourcepathname'))//'Schematic4.data'
open(unit=50,file=fname,status='old',form='unformatted')
read(50) dims
read(50) chararr
close(unit=50,status='keep')
schematic = ichar(chararr)

dataset = 'Schematic 4'
width = dims(2)
height = dims(3)
call h5immake_image_24bit_f(HDF%getObjectID(),dataset,width,height,'INTERLACE_PIXEL',int(schematic),hdferr)

deallocate(schematic,chararr)
!=====================================================
! and finally the selected type
dataset = SC_ID
hdferr = HDF%writeDatasetInteger(dataset, 2)

call HDF%pop()

end subroutine h5_writeCoordinateSystemGroup

!--------------------------------------------------------------------------
subroutine h5_writePatternCenterGroup(HDF, xpc, ypc, L, delta, scdim)
!DEC$ ATTRIBUTES DLLEXPORT :: h5_writePatternCenterGroup
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! write the pattern center group

use mod_io

IMPLICIT NONE

type(HDF_T),INTENT(INOUT)            :: HDF
real(kind=sgl),INTENT(IN)            :: xpc      ! pattern center x [pixels]
real(kind=sgl),INTENT(IN)            :: ypc      ! pattern center y [pixels]
real(kind=sgl),INTENT(IN)            :: L        ! sample-scintillator distance [micron]
real(kind=sgl),INTENT(IN)            :: delta    ! scintillator pixel size [micron]
integer(kind=irg),INTENT(IN)         :: scdim(2) ! scintillator dimensions [pixels]

character(fnlen)                     :: groupname, dataset
integer(kind=irg)                    :: hdferr
real(kind=sgl)                       :: xstar, ystar, zstar

! create the Coordinate System group
groupname = 'Pattern Center Calibration'
hdferr = HDF%createGroup(groupname)

! we assume that we are writing a TSL file

! in EMsoft, the pattern center is measured in units of pixels from the
! center of the scintillator.  For TSL, the pattern center is measured
! from the bottom left of the scintillator (when looking towards it from the
! sample) and in units of the width of the scintillator.

xstar = ( float(scdim(1))*0.5 + xpc ) / float(scdim(1))
ystar = ( float(scdim(2))*0.5 + ypc ) / float(scdim(2))
zstar = L / ( delta * float(scdim(1)) )

dataset = SC_xstar
hdferr = HDF%writeDatasetFloat(dataset, xstar)

dataset = SC_ystar
hdferr = HDF%writeDatasetFloat(dataset, ystar)

dataset = SC_zstar
hdferr = HDF%writeDatasetFloat(dataset, zstar)

call HDF%pop()

end subroutine h5_writePatternCenterGroup

!--------------------------------------------------------------------------
subroutine h5_writePhaseGroup(EMsoft, HDF, groupname, xtalname)
!DEC$ ATTRIBUTES DLLEXPORT :: h5ebsd_writePhaseGroup
 !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! write the phase group, describing the crystal structure

use mod_EMsoft
use mod_io
use mod_crystallography
use mod_symmetry

IMPLICIT NONE

type(EMsoft_T),INTENT(INOUT)                      :: EMsoft
type(HDF_T),INTENT(INOUT)                         :: HDF
character(fnlen),intent(IN)                       :: groupname
character(fnlen),intent(IN)                       :: xtalname

type(HDF_T)                                       :: localHDF
type(cell_T)                                      :: cell
type(SpaceGroup_T)                                :: SG

character(fnlen)                                  :: dataset, grname, filename
integer(kind=irg)                                 :: istat, SGnum, hdferr
real(kind=dbl),allocatable                        :: cellparams(:)
integer(HSIZE_T)                                  :: dims(1)
logical                                           :: readonly, stat


integer(kind=irg)                                 :: i, pgnum
character(fnlen, KIND=c_char),allocatable,TARGET  :: stringarray(:)

! TSL point group labels [courtesy of S. Wright]
character(26),parameter       :: TSLpgname(32) = (/ "Triclinic (C1) [1]        ", "Triclinic (S2, Ci) [-1]   ",&
                        "Monoclinic b (C2)[2]      ", "Monoclinic b (C1h, Cs) [m]", "Monoclinic b (C2h) [2/m]  ",&
                        "Orthorhombic (D2) [222]   ", "Orthorhombic (C2v) [mm2]  ", "Orthorhombic (D2h) [mmm]  ",&
                        "Tetragonal (C4) [4]       ", "Tetragonal (S4) [-4]      ", "Tetragonal (C4h) [4/m]    ",&
                        "Tetragonal (D4) [422]     ", "Tetragonal (C4v) [4mm]    ", "Tetragonal (D2d) [-42m]   ",&
                        "Tetragonal (D4h) [4/mmm]  ", "Trigonal (C3) [3]         ", "Trigonal (S6, C3i) [-3]   ",&
                        "Trigonal (D3) [32]        ", "Trigonal (C3v) [3m]       ", "Trigonal (D3d) [-3m]      ",&
                        "Hexagonal (C6) [6]        ", "Hexagonal (C3h) [-6]      ", "Hexagonal (C6h) [6/m]     ",&
                        "Hexagonal (D6) [622]      ", "Hexagonal (C6v) [6mm]     ", "Hexagonal (D3h) [-6m2]    ",&
                        "Hexagonal (D6h) [6/mmm]   ", "Cubic (T) [23]            ", "Cubic (Th) [m3]           ",&
                        "Cubic (O) [432]           ", "Cubic (Td) [-43m]         ", "Cubic (Oh) [m3m]          " /)

! TSL old symmetry identifiers [courtesy of S. Wright]
integer(kind=irg),parameter    :: TSLoldID(32) = (/ 1,1,2,2,2,22,22,22,4,4,4,42,42,42,42,3,3, &
                                                  32,32,32,6,6,6,62,62,62,62,23,23,43,43,43 /)


! this routine first extracts information from the xtal file and then
! puts it in the right format for the h5ebsd file format.
! This is organized by phase, so each phase is a separate numbered
! subgroup; the subgroupname is passed in as groupname

! test to make sure the input file exists and is HDF5 format
cell = cell_T( )
call cell%setFileName(xtalname)
call cell%readDataHDF(SG, EMsoft, useHDF=localHDF)

! create the subgroup [now we are back in the original HDF5 file]
hdferr = HDF%createGroup(groupname)

! the following data sets need to be created: Formula, Info, Lattice Constant a,
! b, c, alpha, beta, gamma, Laue Group, MaterialName, NumberFamilies, Point Group,
! Symmetry, hkl Families.  These last ones are typically used by the EDAX/TSL
! software, so we do not necessarily have to fill them in here.

cellparams = cell%getLatParm()
sgnum = SG%getSpaceGroupNumber()

! lattice parameters [in Angstrom]
dataset = 'Lattice Constant a'
hdferr = HDF%writeDatasetFloat(dataset, sngl(cellparams(1))*10.0)
dataset = 'Lattice Constant b'
hdferr = HDF%writeDatasetFloat(dataset, sngl(cellparams(2))*10.0)
dataset = 'Lattice Constant c'
hdferr = HDF%writeDatasetFloat(dataset, sngl(cellparams(3))*10.0)
dataset = 'Lattice Constant alpha'
hdferr = HDF%writeDatasetFloat(dataset, sngl(cellparams(4)))
dataset = 'Lattice Constant beta'
hdferr = HDF%writeDatasetFloat(dataset, sngl(cellparams(5)))
dataset = 'Lattice Constant gamma'
hdferr = HDF%writeDatasetFloat(dataset, sngl(cellparams(6)))

allocate(stringarray(1))

! point group
pgnum = 0
do i=1,32
  if (SGPG(i).le.sgnum) pgnum = i
end do
dataset = 'Point Group'
stringarray(1)= trim(TSLpgname(pgnum))
hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

! Laue group
dataset = 'Laue Group'
stringarray(1)= trim(TSLpgname(PGrot(pgnum)))
hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

! Symmetry
dataset = SC_Symmetry
hdferr = HDF%writeDatasetInteger(dataset, TSLoldID(pgnum))

! various other strings

! Formula [extract this from the first part of xtalname]
dataset = SC_Formula
i = scan(trim(xtalname),'.')
stringarray(1) = xtalname(1:i-1)
hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

! Material name [same as Formula; this would require adding a field to the .xtal files]
dataset = SC_MaterialName
stringarray(1) = xtalname(1:i-1)
hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

! Info [empty string most of the time]
dataset = SC_Info
stringarray(1) = ''
hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

! hkl Families [this will require a bit of work !!!!!]
! this item uses the Compound data type; we will need to generate the
! families of unique planes, and compute structure factors ...

! in this version of the software [EMsoft 3.1], we leave these datasets empty
dataset = SC_NumberFamilies
i = 0
hdferr = HDF%writeDatasetInteger(dataset, i)
! call Message('h5ebsd_writePhaseGroup: writing of ->NumberFamilies<- data not yet implemented.')

dataset = 'hkl Families'
hdferr = HDF%writeDatasetInteger(dataset, i)
! call Message('h5ebsd_writePhaseGroup: writing of ->hkl Families<- data not yet implemented.')


! and leave this group
call HDF%pop()

end subroutine h5_writePhaseGroup

end module mod_DIfiles
