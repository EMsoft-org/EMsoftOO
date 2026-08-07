! ###################################################################
! Copyright (c) 2026-2026, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_icoMP
  !! author: MDG
  !! version: 1.0
  !! date: 08/07/26
  !!
  !! class definition for the resampling of a square Lambert master pattern onto
  !! a subdivided icosahedron (icosphere)
  !!
  !! This module is a Fortran 2018 port of the master_to_icosphere.py program; it
  !! reads an EBSD/ECP/TKD master pattern file (through mod_MPfiles), resamples the
  !! two square Lambert hemispheres onto the vertices of an icosphere, and writes
  !! the result as a numpy .npz file (through mod_npz).  Two representations are
  !! produced: a flat vertex list signal(C,V) for mesh convolutions, and a set of
  !! five rectangular charts chart_signal(C,5,H,W) with H = 2^L and W = 2^(L+1)
  !! for icosahedral CNNs.
  !!
  !! All index arrays that end up in the .npz file (faces, neighbours, charts,
  !! pole_index and the Laplacian row/column arrays) are ZERO-BASED, since they are
  !! meant to be consumed by python.  For the same reason most of the internal
  !! arrays carry a zero lower bound as well; this keeps the port a line-by-line
  !! match with the python original, which is what guarantees that the vertex
  !! numbering (and hence everything derived from it) is identical.

use mod_kinds
use mod_global

IMPLICIT NONE

private

  integer(kind=irg), parameter        :: icoMPmaxlevel = 8
   !! largest allowed subdivision level (655362 vertices)
  integer(kind=irg), parameter        :: icoMPmaxdegree = 8
   !! upper bound on the number of edges incident on a vertex (really 5 or 6)
  integer(kind=irg), parameter        :: icoMPmaxweights = 100
   !! size of the energy weight array in the namelist
  character(len=*), parameter         :: icoMPversion = '1.0'

! namelist for the EMicoMP program
type, public :: icoMPNameListType
  integer(kind=irg)  :: level
  integer(kind=irg)  :: energybin
  integer(kind=irg)  :: aarings
  integer(kind=irg)  :: aaazim
  real(kind=sgl)     :: energyweights(icoMPmaxweights)
  logical            :: antialias
  logical            :: doubleprecision
  logical            :: laplacian
  logical            :: verify
  logical            :: selftest
  character(fnlen)   :: masterfile
  character(fnlen)   :: npzfile
  character(fnlen)   :: modality
  character(fnlen)   :: energymode
  character(fnlen)   :: normalize
  character(fnlen)   :: interpolation
end type icoMPNameListType

type, public :: icoMP_T
  private
    character(fnlen)                :: nmldeffile = 'EMicoMP.nml'
    type(icoMPNameListType)         :: nml

! control parameters
    integer(kind=irg)               :: level = 5
     !! subdivision level of the icosphere
    character(fnlen)                :: energymode = 'sum'
     !! 'sum', 'all', 'weighted' or 'index'
    integer(kind=irg)               :: energybin = 0
     !! ZERO-BASED energy bin number, used when energymode = 'index'
    real(kind=dbl), allocatable     :: eweights(:)
     !! energy weight factors, used when energymode = 'weighted'
    character(fnlen)                :: normmode = 'none'
     !! 'none', 'zscore', 'minmax' or 'mean1'
    character(fnlen)                :: interpmode = 'auto'
     !! 'auto', 'bilinear' or 'bicubic'
    character(fnlen)                :: kernel = 'bilinear'
     !! the interpolation kernel actually used (resolved from interpmode)
    logical                         :: antialias = .TRUE.
     !! average over the vertex cell instead of point sampling
    integer(kind=irg)               :: aarings = 6
     !! number of quadrature rings for the cell average
    integer(kind=irg)               :: aaazim = 10
     !! number of azimuthal quadrature points for the cell average
    logical                         :: doubleprec = .FALSE.
     !! write float64 instead of float32 for the real output arrays
    logical                         :: dolaplacian = .FALSE.
     !! also compute and store the cotangent Laplacian

! icosphere geometry
    integer(kind=irg)               :: nverts = 0
    integer(kind=irg)               :: nfaces = 0
    integer(kind=irg)               :: hchart = 0
     !! number of chart rows, including the shared border; 2^level+1
    integer(kind=irg)               :: wchart = 0
     !! number of chart columns, including the shared border; 2^(level+1)+1
    real(kind=dbl), allocatable     :: verts(:,:)
     !! (0:nverts-1,3) unit vectors
    integer(kind=irg), allocatable  :: faces(:,:)
     !! (0:nfaces-1,3) zero-based vertex indices
    integer(kind=irg), allocatable  :: neigh(:,:)
     !! (0:nverts-1,6) zero-based neighbours in CCW order, -1 padded
    integer(kind=irg), allocatable  :: chrt(:,:,:)
     !! (0:4,0:hchart-1,0:wchart-1) chart index lattices
    integer(kind=irg), allocatable  :: counts(:)
     !! (0:level) number of vertices at each level

! master pattern
    integer(kind=irg)               :: npx = 0
     !! master pattern semi-edge length; arrays are (-npx:npx,-npx:npx,...)
    integer(kind=irg)               :: numEbins = 0
    real(kind=dbl), allocatable     :: mLPNH(:,:,:)
    real(kind=dbl), allocatable     :: mLPSH(:,:,:)
    real(kind=dbl), allocatable     :: keVs(:)
    logical                         :: haskeVs = .FALSE.
    logical                         :: hasmaster = .FALSE.
    character(fnlen)                :: masterfile = ''
    character(fnlen)                :: modality = 'EBSD'
    real(kind=dbl)                  :: eqmismatch = 0.D0

! energy-collapsed master pattern; the channel axis leads for fast interpolation
    integer(kind=irg)               :: nchannels = 0
    real(kind=dbl), allocatable     :: nhsel(:,:,:)
     !! (1:nchannels,-npx:npx,-npx:npx)
    real(kind=dbl), allocatable     :: shsel(:,:,:)
    real(kind=dbl), allocatable     :: energies(:)

! results
    real(kind=dbl), allocatable     :: signal(:,:)
     !! (1:nchannels,0:nverts-1)
    real(kind=dbl), allocatable     :: chartsignal(:,:,:,:)
     !! (1:nchannels,5,2^level,2^(level+1))
    integer(kind=irg), allocatable  :: poleindex(:)
    integer(kind=irg), allocatable  :: laprow(:)
    integer(kind=irg), allocatable  :: lapcol(:)
    real(kind=dbl), allocatable     :: lapval(:)
    real(kind=dbl), allocatable     :: vertexarea(:)

  contains
  private
    procedure, pass(self) :: readNameList_
    procedure, pass(self) :: getNameList_
    procedure, pass(self) :: icoMP_
    procedure, pass(self) :: setLevel_
    procedure, pass(self) :: getLevel_
    procedure, pass(self) :: setEnergyMode_
    procedure, pass(self) :: setNormalization_
    procedure, pass(self) :: setInterpolation_
    procedure, pass(self) :: setAntialias_
    procedure, pass(self) :: setDoublePrecision_
    procedure, pass(self) :: setLaplacian_
    procedure, pass(self) :: getNumVertices_
    procedure, pass(self) :: getNumFaces_
    procedure, pass(self) :: getKernel_
    procedure, pass(self) :: getEquatorMismatch_
    procedure, pass(self) :: buildIcosphere_
    procedure, pass(self) :: cotangentLaplacian_
    procedure, pass(self) :: readMasterPattern_
    procedure, pass(self) :: setMasterPattern_
    procedure, pass(self) :: equatorMismatch_
    procedure, pass(self) :: pixelDegrees_
    procedure, pass(self) :: selectEnergy_
    procedure, pass(self) :: chooseKernel_
    procedure, pass(self) :: sampleDirections_
    procedure, pass(self) :: sampleCellAverage_
    procedure, pass(self) :: normalizeSignal_
    procedure, pass(self) :: projectMaster_
    procedure, pass(self) :: metaJSON_
    procedure, pass(self) :: writeNPZ_
    procedure, pass(self) :: selfTest_
    procedure, pass(self) :: verifySampling_

    generic, public :: readNameList => readNameList_
    generic, public :: getNameList => getNameList_
    generic, public :: icoMP => icoMP_
    generic, public :: setLevel => setLevel_
    generic, public :: getLevel => getLevel_
    generic, public :: setEnergyMode => setEnergyMode_
    generic, public :: setNormalization => setNormalization_
    generic, public :: setInterpolation => setInterpolation_
    generic, public :: setAntialias => setAntialias_
    generic, public :: setDoublePrecision => setDoublePrecision_
    generic, public :: setLaplacian => setLaplacian_
    generic, public :: getNumVertices => getNumVertices_
    generic, public :: getNumFaces => getNumFaces_
    generic, public :: getKernel => getKernel_
    generic, public :: getEquatorMismatch => getEquatorMismatch_
    generic, public :: buildIcosphere => buildIcosphere_
    generic, public :: cotangentLaplacian => cotangentLaplacian_
    generic, public :: readMasterPattern => readMasterPattern_
    generic, public :: setMasterPattern => setMasterPattern_
    generic, public :: equatorMismatch => equatorMismatch_
    generic, public :: pixelDegrees => pixelDegrees_
    generic, public :: sampleDirections => sampleDirections_
    generic, public :: projectMaster => projectMaster_
    generic, public :: writeNPZ => writeNPZ_
    generic, public :: selfTest => selfTest_
    generic, public :: verifySampling => verifySampling_

end type icoMP_T

! the constructor routines for this class
interface icoMP_T
  module procedure icoMP_constructor
  module procedure icoMP_nml_constructor
end interface icoMP_T

contains

!--------------------------------------------------------------------------
type(icoMP_T) function icoMP_constructor( level ) result(icoMP)
!DEC$ ATTRIBUTES DLLEXPORT :: icoMP_constructor
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! constructor for the icoMP_T Class

IMPLICIT NONE

integer(kind=irg), INTENT(IN), OPTIONAL   :: level

if (present(level)) call icoMP%setLevel_(level)

end function icoMP_constructor

!--------------------------------------------------------------------------
type(icoMP_T) function icoMP_nml_constructor( nmlfile ) result(icoMP)
!DEC$ ATTRIBUTES DLLEXPORT :: icoMP_nml_constructor
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! constructor for the icoMP_T Class; reads the name list

IMPLICIT NONE

character(fnlen), INTENT(IN)   :: nmlfile

icoMP%nmldeffile = trim(nmlfile)
call icoMP%readNameList_(nmlfile)

end function icoMP_nml_constructor

!--------------------------------------------------------------------------
recursive subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! read the namelist from an nml file for the icoMP_T Class

use mod_io
use mod_EMsoft

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
character(fnlen), INTENT(IN)    :: nmlfile
 !! full path to namelist file
logical, OPTIONAL, INTENT(IN)   :: initonly
 !! fill in the default values only; do not read the file

type(IO_T)                      :: Message
logical                         :: skipread = .FALSE.

integer(kind=irg)               :: level, energybin, aarings, aaazim
real(kind=sgl)                  :: energyweights(icoMPmaxweights)
logical                         :: antialias, doubleprecision, laplacian, verify, selftest
character(fnlen)                :: masterfile, npzfile, modality, energymode, normalize, interpolation

! define the IO namelist to facilitate passing variables to the program.
namelist / icoMPdata / level, energybin, aarings, aaazim, energyweights, antialias, &
                       doubleprecision, laplacian, verify, selftest, masterfile, npzfile, &
                       modality, energymode, normalize, interpolation

! set the input parameters to default values
level = 5
energybin = 0
aarings = 6
aaazim = 10
energyweights = 0.0
antialias = .TRUE.
doubleprecision = .FALSE.
laplacian = .FALSE.
verify = .FALSE.
selftest = .FALSE.
masterfile = 'undefined'
npzfile = 'undefined'
modality = 'undefined'
energymode = 'sum'
normalize = 'none'
interpolation = 'auto'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=icoMPdata)
    close(UNIT=dataunit,STATUS='keep')

! check for required entries; the self test needs neither of them
    if (selftest.eqv..FALSE.) then
      if (trim(masterfile).eq.'undefined') then
        call Message%printError('readNameList:', ' master pattern file name is undefined in '//nmlfile)
      end if
      if (trim(npzfile).eq.'undefined') then
        call Message%printError('readNameList:', ' output npz file name is undefined in '//nmlfile)
      end if
    end if

! catch the misspelled keyword strings here rather than after the master pattern
! has been read; the energy bin number can only be checked once the number of
! energy bins is known, so that one is left to selectEnergy
    select case (trim(energymode))
      case ('sum', 'all', 'index', 'weighted')
      case default
        call Message%printError('readNameList:', ' unknown energymode '//trim(energymode)//' in '//nmlfile)
    end select

    select case (trim(modality))
      case ('undefined', 'EBSD', 'ECP', 'TKD')
      case default
        call Message%printError('readNameList:', ' unknown modality '//trim(modality)//' in '//nmlfile)
    end select
end if

! if we get here, then all appears to be ok, and we need to fill in the nml fields
self%nml%level = level
self%nml%energybin = energybin
self%nml%aarings = aarings
self%nml%aaazim = aaazim
self%nml%energyweights = energyweights
self%nml%antialias = antialias
self%nml%doubleprecision = doubleprecision
self%nml%laplacian = laplacian
self%nml%verify = verify
self%nml%selftest = selftest
self%nml%masterfile = masterfile
self%nml%npzfile = npzfile
self%nml%modality = modality
self%nml%energymode = energymode
self%nml%normalize = normalize
self%nml%interpolation = interpolation

end subroutine readNameList_

!--------------------------------------------------------------------------
recursive function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! pass the namelist for the icoMP_T Class to the calling program

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
type(icoMPNameListType)         :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine icoMP_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: icoMP_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! resample a master pattern onto an icosphere and write it as a .npz file
!!
!! This is the entry point used by the EMicoMP program; it transfers the namelist
!! parameters into the class, reads the master pattern, and runs the pipeline.

use mod_EMsoft
use mod_HDFsupport
use mod_io
use mod_memory

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
type(EMsoft_T), INTENT(INOUT)   :: EMsoft
character(fnlen), INTENT(IN)    :: progname

type(IO_T)                      :: Message
type(HDF_T)                     :: HDF
type(memory_T)                  :: mem
real(kind=dbl), allocatable     :: wts(:)
character(fnlen)                :: fname, charline
integer(kind=irg)               :: io_int(3), nw
real(kind=dbl)                  :: io_dbl(2)
logical                         :: ok

associate( nml => self%nml )

mem = memory_T()

! the geometry self test needs no master pattern at all, so it short circuits
if (nml%selftest.eqv..TRUE.) then
  call self%selfTest_(min(nml%level, 6))
  return
end if

!--------------------------------------------------------------------------
! transfer the namelist parameters into the class
!--------------------------------------------------------------------------
call self%setLevel_(nml%level)
call self%setNormalization_(trim(nml%normalize))
call self%setInterpolation_(trim(nml%interpolation))
call self%setAntialias_(nml%antialias, rings = nml%aarings, azim = nml%aaazim)
call self%setDoublePrecision_(nml%doubleprecision)
call self%setLaplacian_(nml%laplacian)

!--------------------------------------------------------------------------
! read the master pattern
!--------------------------------------------------------------------------
call openFortranHDFInterface()
HDF = HDF_T()

if (trim(nml%modality).eq.'undefined') then
  call self%readMasterPattern_(EMsoft, HDF, nml%masterfile)
else
  call self%readMasterPattern_(EMsoft, HDF, nml%masterfile, trim(nml%modality))
end if

call closeFortranHDFInterface()

call Message%printMessage(' Master pattern  : '//trim(self%masterfile))
call Message%printMessage(' Modality        : '//trim(self%modality))
io_int(1) = self%npx
io_int(2) = 2*self%npx+1
io_int(3) = self%numEbins
call Message%WriteValue(' npx, edge, bins : ', io_int, 3, "(I6,I8,I6)")
io_dbl(1) = self%eqmismatch
io_dbl(2) = self%pixelDegrees_()
call Message%WriteValue(' NH/SH mismatch, pixel size [deg] : ', io_dbl, 2, "(ES12.4,F10.5)")

!--------------------------------------------------------------------------
! the energy mode needs the number of energy bins, so it is set after the read
!--------------------------------------------------------------------------
select case (trim(nml%energymode))
  case ('sum', 'all')
    call self%setEnergyMode_(trim(nml%energymode))
  case ('index')
    call self%setEnergyMode_('index', bin = nml%energybin)
  case ('weighted')
    nw = self%numEbins
    if (nw.gt.icoMPmaxweights) then
      call Message%printError('icoMP', 'the master pattern has more energy bins than the weight array can hold')
    end if
    call mem%alloc(wts, (/ nw /), 'wts')
    wts = dble(nml%energyweights(1:nw))
    if (sum(dabs(wts)).eq.0.D0) then
      call Message%printError('icoMP', 'energymode is weighted but all energyweights entries are zero')
    end if
    call self%setEnergyMode_('weighted', weights = wts)
    call mem%dealloc(wts, 'wts')
  case default
    call Message%printError('icoMP', 'unknown energymode '//trim(nml%energymode))
end select

!--------------------------------------------------------------------------
! optionally check the sampling against the master pattern's own grid nodes
!--------------------------------------------------------------------------
if (nml%verify.eqv..TRUE.) then
  call Message%printMessage(' Verifying the interpolation ')
  ok = self%verifySampling_()
  if (ok.eqv..FALSE.) then
    call Message%printError('icoMP', 'the sampling verification failed')
  end if
  call Message%printMessage(' Verification passed ')
end if

!--------------------------------------------------------------------------
! and run the pipeline
!--------------------------------------------------------------------------
call self%projectMaster_()

io_int(1) = self%level
io_int(2) = self%nverts
io_int(3) = self%nfaces
call Message%WriteValue(' level, vertices, faces : ', io_int, 3, "(I4,2I10)")
write (charline,"(' chart_signal    : (',I3,', 5,',I5,',',I5,')')") &
      self%nchannels, self%hchart-1, self%wchart-1
call Message%printMessage(trim(charline))
call Message%printMessage(' Interpolation   : '//trim(self%kernel))
if (self%antialias.eqv..TRUE.) then
  write (charline,"(' Sampling        : cell average, ',I3,'x',I3,' equal-area disc quadrature')") &
        self%aarings, self%aaazim
else
  charline = ' Sampling        : point sample'
end if
call Message%printMessage(trim(charline))

fname = EMsoft%generateFilePath('EMdatapathname', trim(nml%npzfile))
call self%writeNPZ_(fname)

end associate

end subroutine icoMP_

!--------------------------------------------------------------------------
recursive subroutine setLevel_(self, level)
!DEC$ ATTRIBUTES DLLEXPORT :: setLevel_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! set the icosphere subdivision level

use mod_io

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
integer(kind=irg), INTENT(IN)   :: level

type(IO_T)                      :: Message
character(fnlen)                :: charline

if ((level.lt.0).or.(level.gt.icoMPmaxlevel)) then
  write (charline,"('requested level = ',I4)") level
  call Message%printError('setLevel', 'subdivision level must lie in the range 0-8; '//trim(charline))
end if

self%level = level

end subroutine setLevel_

!--------------------------------------------------------------------------
recursive function getLevel_(self) result(level)
!DEC$ ATTRIBUTES DLLEXPORT :: getLevel_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! get the icosphere subdivision level

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
integer(kind=irg)               :: level

level = self%level

end function getLevel_

!--------------------------------------------------------------------------
recursive subroutine setEnergyMode_(self, mode, bin, weights)
!DEC$ ATTRIBUTES DLLEXPORT :: setEnergyMode_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! set the way in which the energy axis of the master pattern is collapsed
!!
!! mode is one of 'sum' (add all energy bins), 'all' (keep all bins as separate
!! channels), 'weighted' (weighted sum, requires weights) or 'index' (a single
!! bin, requires bin).  NOTE that bin is ZERO-BASED, to match the python program.

use mod_io

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)             :: self
character(*), INTENT(IN)                  :: mode
integer(kind=irg), INTENT(IN), OPTIONAL   :: bin
real(kind=dbl), INTENT(IN), OPTIONAL      :: weights(:)

type(IO_T)                                :: Message

select case (trim(mode))
  case ('sum', 'all')
    self%energymode = trim(mode)
  case ('index')
    if (.not.present(bin)) then
      call Message%printError('setEnergyMode', 'energy mode index requires the bin argument')
    end if
    self%energymode = 'index'
    self%energybin = bin
  case ('weighted')
    if (.not.present(weights)) then
      call Message%printError('setEnergyMode', 'energy mode weighted requires the weights argument')
    end if
    self%energymode = 'weighted'
    if (allocated(self%eweights)) deallocate(self%eweights)
    allocate(self%eweights(size(weights)))
    self%eweights = weights
  case default
    call Message%printError('setEnergyMode', 'unknown energy mode '//trim(mode))
end select

end subroutine setEnergyMode_

!--------------------------------------------------------------------------
recursive subroutine setNormalization_(self, mode)
!DEC$ ATTRIBUTES DLLEXPORT :: setNormalization_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! set the per-channel normalization mode: 'none', 'zscore', 'minmax' or 'mean1'

use mod_io

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
character(*), INTENT(IN)        :: mode

type(IO_T)                      :: Message

select case (trim(mode))
  case ('none', 'zscore', 'minmax', 'mean1')
    self%normmode = trim(mode)
  case default
    call Message%printError('setNormalization', 'unknown normalization mode '//trim(mode))
end select

end subroutine setNormalization_

!--------------------------------------------------------------------------
recursive subroutine setInterpolation_(self, mode)
!DEC$ ATTRIBUTES DLLEXPORT :: setInterpolation_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! set the interpolation kernel: 'auto', 'bilinear' or 'bicubic'

use mod_io

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
character(*), INTENT(IN)        :: mode

type(IO_T)                      :: Message

select case (trim(mode))
  case ('auto', 'bilinear', 'bicubic')
    self%interpmode = trim(mode)
  case default
    call Message%printError('setInterpolation', 'unknown interpolation mode '//trim(mode))
end select

end subroutine setInterpolation_

!--------------------------------------------------------------------------
recursive subroutine setAntialias_(self, antialias, rings, azim)
!DEC$ ATTRIBUTES DLLEXPORT :: setAntialias_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! turn the cell average on or off and set the quadrature size

use mod_io

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)             :: self
logical, INTENT(IN)                       :: antialias
integer(kind=irg), INTENT(IN), OPTIONAL   :: rings
integer(kind=irg), INTENT(IN), OPTIONAL   :: azim

type(IO_T)                                :: Message

self%antialias = antialias

if (present(rings)) then
  if (rings.lt.1) call Message%printError('setAntialias', 'the number of quadrature rings must be positive')
  self%aarings = rings
end if

if (present(azim)) then
  if (azim.lt.1) call Message%printError('setAntialias', 'the number of azimuthal points must be positive')
  self%aaazim = azim
end if

end subroutine setAntialias_

!--------------------------------------------------------------------------
recursive subroutine setDoublePrecision_(self, dp)
!DEC$ ATTRIBUTES DLLEXPORT :: setDoublePrecision_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! write the real output arrays as float64 instead of float32

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
logical, INTENT(IN)             :: dp

self%doubleprec = dp

end subroutine setDoublePrecision_

!--------------------------------------------------------------------------
recursive subroutine setLaplacian_(self, dl)
!DEC$ ATTRIBUTES DLLEXPORT :: setLaplacian_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! also compute and store the cotangent Laplacian and the vertex areas

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
logical, INTENT(IN)             :: dl

self%dolaplacian = dl

end subroutine setLaplacian_

!--------------------------------------------------------------------------
recursive function getNumVertices_(self) result(nv)
!DEC$ ATTRIBUTES DLLEXPORT :: getNumVertices_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! number of icosphere vertices

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
integer(kind=irg)               :: nv

nv = self%nverts

end function getNumVertices_

!--------------------------------------------------------------------------
recursive function getNumFaces_(self) result(nf)
!DEC$ ATTRIBUTES DLLEXPORT :: getNumFaces_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! number of icosphere faces

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
integer(kind=irg)               :: nf

nf = self%nfaces

end function getNumFaces_

!--------------------------------------------------------------------------
recursive function getKernel_(self) result(k)
!DEC$ ATTRIBUTES DLLEXPORT :: getKernel_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! the interpolation kernel that was actually used

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
character(fnlen)                :: k

k = trim(self%kernel)

end function getKernel_

!--------------------------------------------------------------------------
recursive function getEquatorMismatch_(self) result(m)
!DEC$ ATTRIBUTES DLLEXPORT :: getEquatorMismatch_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! the largest NH/SH difference on the square border

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
real(kind=dbl)                  :: m

m = self%eqmismatch

end function getEquatorMismatch_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! icosphere construction
!
! Vertices are ordered N, upper ring U0..U4, lower ring L0..L4, S, so chart i is
! just the 2x3 lattice
!
!     N        U_i      L_(i-1)
!     U_(i+1)  L_i      S
!
! Every row, column and anti-diagonal step in it is a real icosahedron edge, and
! refinement doubles the lattice while keeping that true.
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine baseIcosahedron_(v)
!DEC$ ATTRIBUTES DLLEXPORT :: baseIcosahedron_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! the twelve vertices of the regular icosahedron

IMPLICIT NONE

real(kind=dbl), INTENT(INOUT)   :: v(0:,:)

real(kind=dbl)                  :: a, sa, ca, az, az2, d2r
integer(kind=irg)               :: k

d2r = cPi/180.D0
a = datan(2.D0)
sa = dsin(a)
ca = dcos(a)

v(0,:) = (/ 0.D0, 0.D0, 1.D0 /)

do k = 0, 4
  az = (72.D0 * dble(k)) * d2r
  az2 = az + 36.D0 * d2r
  v(1+k,:) = (/ sa*dcos(az),  sa*dsin(az),   ca /)
  v(6+k,:) = (/ sa*dcos(az2), sa*dsin(az2), -ca /)
end do

v(11,:) = (/ 0.D0, 0.D0, -1.D0 /)

end subroutine baseIcosahedron_

!--------------------------------------------------------------------------
recursive subroutine baseCharts_(g)
!DEC$ ATTRIBUTES DLLEXPORT :: baseCharts_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! the five 2x3 chart lattices of the base icosahedron
!!
!! note the use of modulo() rather than mod(); the python original relies on
!! -1 % 5 = 4, which is what modulo() does and mod() does not

IMPLICIT NONE

integer(kind=irg), INTENT(INOUT)  :: g(0:,0:,0:)

integer(kind=irg)                 :: i

do i = 0, 4
  g(i,0,0) = 0
  g(i,0,1) = 1 + modulo(i,5)
  g(i,0,2) = 6 + modulo(i-1,5)
  g(i,1,0) = 1 + modulo(i+1,5)
  g(i,1,1) = 6 + modulo(i,5)
  g(i,1,2) = 11
end do

end subroutine baseCharts_

!--------------------------------------------------------------------------
recursive subroutine getMidpoint_(i, j, v, nv, mkey, mval, mcnt, idx)
!DEC$ ATTRIBUTES DLLEXPORT :: getMidpoint_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! return the index of the normalized midpoint of the edge (i,j), creating it if
!! this is the first time the edge is encountered
!!
!! The lookup table is a bucket list indexed by the lower of the two vertex
!! numbers.  Since every pair passed to this routine is an edge of the current
!! mesh, and every vertex of an icosphere has degree 5 or 6, each bucket holds at
!! most six entries and a short linear scan is all that is needed.  New vertices
!! are appended in call order, which is what keeps the numbering identical to the
!! python original.

use mod_io

IMPLICIT NONE

integer(kind=irg), INTENT(IN)     :: i
integer(kind=irg), INTENT(IN)     :: j
real(kind=dbl), INTENT(INOUT)     :: v(0:,:)
integer(kind=irg), INTENT(INOUT)  :: nv
integer(kind=irg), INTENT(INOUT)  :: mkey(0:,:)
integer(kind=irg), INTENT(INOUT)  :: mval(0:,:)
integer(kind=irg), INTENT(INOUT)  :: mcnt(0:)
integer(kind=irg), INTENT(OUT)    :: idx

type(IO_T)                        :: Message
integer(kind=irg)                 :: imin, imax, p
real(kind=dbl)                    :: m(3)

imin = min(i,j)
imax = max(i,j)

do p = 1, mcnt(imin)
  if (mkey(imin,p).eq.imax) then
    idx = mval(imin,p)
    return
  end if
end do

if (mcnt(imin).ge.icoMPmaxdegree) then
  call Message%printError('getMidpoint', 'edge bucket overflow; this mesh is not an icosphere')
end if

m = v(i,:) + v(j,:)
m = m / dsqrt(sum(m*m))

v(nv,:) = m
idx = nv
nv = nv + 1

mcnt(imin) = mcnt(imin) + 1
mkey(imin,mcnt(imin)) = imax
mval(imin,mcnt(imin)) = idx

end subroutine getMidpoint_

!--------------------------------------------------------------------------
recursive subroutine buildIcosphere_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: buildIcosphere_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! build the icosphere for the requested subdivision level
!!
!! fills verts, faces, neigh, chrt and counts.  The charts include the border
!! shared with the next chart; slice (:,1:,:-1) in python terms for the disjoint
!! CNN tiling.  Level k is a prefix of level k+1, so pooling is a slice.

use mod_io
use mod_memory

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)     :: self

type(IO_T)                        :: Message
type(memory_T)                    :: mem
integer(kind=irg), allocatable    :: gold(:,:,:), gnew(:,:,:)
integer(kind=irg), allocatable    :: mkey(:,:), mval(:,:), mcnt(:)
integer(kind=irg)                 :: nvmax, nv, h, w, hn, wn, lev, i, r, c, idx, f
real(kind=dbl)                    :: q

mem = memory_T()

nvmax = 10 * 4**self%level + 2
self%nfaces = 20 * 4**self%level

! the vertex array is allocated at its final size right away; getMidpoint_ fills
! it in as new midpoints are created
call mem%alloc(self%verts, (/ nvmax-1, 3 /), 'verts', startdims = (/ 0, 1 /))
call mem%alloc(self%counts, (/ self%level /), 'counts', startdims = (/ 0 /))

call baseIcosahedron_(self%verts)
nv = 12
self%counts(0) = nv

h = 2
w = 3
call mem%alloc(gold, (/ 4, h-1, w-1 /), 'gold', startdims = (/ 0, 0, 0 /))
call baseCharts_(gold)

do lev = 1, self%level
  hn = 2*h - 1
  wn = 2*w - 1
  call mem%alloc(gnew, (/ 4, hn-1, wn-1 /), 'gnew', startdims = (/ 0, 0, 0 /))
  call mem%alloc(mkey, (/ nv-1, icoMPmaxdegree /), 'mkey', startdims = (/ 0, 1 /))
  call mem%alloc(mval, (/ nv-1, icoMPmaxdegree /), 'mval', startdims = (/ 0, 1 /))
  call mem%alloc(mcnt, (/ nv-1 /), 'mcnt', initval = 0, startdims = (/ 0 /))

! the loop order below must not be changed; it fixes the order in which new
! vertices are created and hence the entire vertex numbering
  do i = 0, 4
    do r = 0, h-1
      do c = 0, w-1
        gnew(i,2*r,2*c) = gold(i,r,c)
      end do
    end do
    do r = 0, h-1
      do c = 0, w-2
        call getMidpoint_(gold(i,r,c), gold(i,r,c+1), self%verts, nv, mkey, mval, mcnt, idx)
        gnew(i,2*r,2*c+1) = idx
      end do
    end do
    do r = 0, h-2
      do c = 0, w-1
        call getMidpoint_(gold(i,r,c), gold(i,r+1,c), self%verts, nv, mkey, mval, mcnt, idx)
        gnew(i,2*r+1,2*c) = idx
      end do
    end do
    do r = 0, h-2
      do c = 0, w-2
        call getMidpoint_(gold(i,r,c+1), gold(i,r+1,c), self%verts, nv, mkey, mval, mcnt, idx)
        gnew(i,2*r+1,2*c+1) = idx
      end do
    end do
  end do

  call mem%dealloc(mkey, 'mkey')
  call mem%dealloc(mval, 'mval')
  call mem%dealloc(mcnt, 'mcnt')

  call mem%dealloc(gold, 'gold')
  call mem%alloc(gold, (/ 4, hn-1, wn-1 /), 'gold', startdims = (/ 0, 0, 0 /))
  gold = gnew
  call mem%dealloc(gnew, 'gnew')

  h = hn
  w = wn
  self%counts(lev) = nv
end do

if (nv.ne.nvmax) then
  call Message%printError('buildIcosphere', 'unexpected number of vertices; the refinement went wrong')
end if

self%nverts = nv
self%hchart = h
self%wchart = w

! renormalize, exactly as the python original does
do i = 0, nv-1
  q = dsqrt(sum(self%verts(i,:)**2))
  self%verts(i,:) = self%verts(i,:) / q
end do

call mem%alloc(self%chrt, (/ 4, h-1, w-1 /), 'chrt', startdims = (/ 0, 0, 0 /))
self%chrt = gold
call mem%dealloc(gold, 'gold')

! split every cell along its anti-diagonal; the two triangle blocks of a chart
! are emitted one after the other, and the cells are visited in row-major order
call mem%alloc(self%faces, (/ self%nfaces-1, 3 /), 'faces', startdims = (/ 0, 1 /))
f = 0
do i = 0, 4
  do r = 0, h-2
    do c = 0, w-2
      self%faces(f,1) = self%chrt(i,r,c)
      self%faces(f,2) = self%chrt(i,r,c+1)
      self%faces(f,3) = self%chrt(i,r+1,c)
      f = f + 1
    end do
  end do
  do r = 0, h-2
    do c = 0, w-2
      self%faces(f,1) = self%chrt(i,r,c+1)
      self%faces(f,2) = self%chrt(i,r+1,c+1)
      self%faces(f,3) = self%chrt(i,r+1,c)
      f = f + 1
    end do
  end do
end do

if (f.ne.self%nfaces) then
  call Message%printError('buildIcosphere', 'unexpected number of faces')
end if

call buildNeighbours_(self%verts, self%nverts, self%faces, self%nfaces, self%neigh, mem)

end subroutine buildIcosphere_

!--------------------------------------------------------------------------
recursive subroutine buildNeighbours_(v, nv, faces, nf, neigh, mem)
!DEC$ ATTRIBUTES DLLEXPORT :: buildNeighbours_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! (nv,6) neighbour table in CCW order, -1 padded for the twelve pentagons
!!
!! The tangent frame is anchored on the lowest-numbered neighbour, so that the
!! cyclic order is deterministic, which is all a gauge-equivariant kernel needs.

use mod_io
use mod_memory

IMPLICIT NONE

real(kind=dbl), INTENT(IN)                    :: v(0:,:)
integer(kind=irg), INTENT(IN)                 :: nv
integer(kind=irg), INTENT(IN)                 :: faces(0:,:)
integer(kind=irg), INTENT(IN)                 :: nf
integer(kind=irg), INTENT(INOUT), allocatable :: neigh(:,:)
type(memory_T), INTENT(INOUT)                 :: mem

type(IO_T)                                    :: Message
integer(kind=irg), allocatable                :: adj(:,:), adjcnt(:)
integer(kind=irg)                             :: i, k, p, q, deg, ns(icoMPmaxdegree), tmp
real(kind=dbl)                                :: n(3), e1(3), e2(3), d(3), w0(3)
real(kind=dbl)                                :: ang(icoMPmaxdegree), tang, twopi

twopi = 2.D0 * cPi

call mem%alloc(adj, (/ nv-1, icoMPmaxdegree /), 'adj', startdims = (/ 0, 1 /))
call mem%alloc(adjcnt, (/ nv-1 /), 'adjcnt', initval = 0, startdims = (/ 0 /))

do i = 0, nf-1
  call addNeighbour_(adj, adjcnt, faces(i,1), faces(i,2))
  call addNeighbour_(adj, adjcnt, faces(i,1), faces(i,3))
  call addNeighbour_(adj, adjcnt, faces(i,2), faces(i,1))
  call addNeighbour_(adj, adjcnt, faces(i,2), faces(i,3))
  call addNeighbour_(adj, adjcnt, faces(i,3), faces(i,1))
  call addNeighbour_(adj, adjcnt, faces(i,3), faces(i,2))
end do

call mem%alloc(neigh, (/ nv-1, 6 /), 'neigh', initval = -1, startdims = (/ 0, 1 /))

do i = 0, nv-1
  deg = adjcnt(i)
  if (deg.gt.6) then
    call Message%printError('buildNeighbours', 'vertex with more than six neighbours')
  end if
  ns(1:deg) = adj(i,1:deg)
  n = v(i,:)

! tangent frame anchored on the lowest-indexed neighbour
  w0 = v(ns(1),:)
  e1 = w0 - n * dot_product(w0, n)
  e1 = e1 / dsqrt(sum(e1*e1))
  e2 = cross_(n, e1)

! the anchor lies along e1 by construction, so its azimuth is exactly zero and is
! set rather than computed.  Evaluating it would give atan2 of a quantity that is
! mathematically zero but numerically a tiny number of either sign, and the mod
! then sends it to either 0 or 2*pi, which rotates the whole list by one entry.
! The python original does compute it, and as a result starts the cycle on the
! anchor for only about half of the vertices; the cyclic order itself is the same
! either way, and it is all that a gauge-equivariant kernel uses.
  ang(1) = 0.D0
  do p = 2, deg
    d = v(ns(p),:) - n * dot_product(v(ns(p),:), n)
    ang(p) = modulo(datan2(dot_product(d, e2), dot_product(d, e1)), twopi)
  end do

! insertion sort on the azimuthal angle
  do p = 2, deg
    tang = ang(p)
    tmp = ns(p)
    q = p - 1
    do while (q.ge.1)
      if (ang(q).le.tang) exit
      ang(q+1) = ang(q)
      ns(q+1) = ns(q)
      q = q - 1
    end do
    ang(q+1) = tang
    ns(q+1) = tmp
  end do

  do k = 1, deg
    neigh(i,k) = ns(k)
  end do
end do

call mem%dealloc(adj, 'adj')
call mem%dealloc(adjcnt, 'adjcnt')

end subroutine buildNeighbours_

!--------------------------------------------------------------------------
recursive subroutine addNeighbour_(adj, adjcnt, i, j)
!DEC$ ATTRIBUTES DLLEXPORT :: addNeighbour_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! add j to the sorted neighbour list of i, skipping duplicates

use mod_io

IMPLICIT NONE

integer(kind=irg), INTENT(INOUT)  :: adj(0:,:)
integer(kind=irg), INTENT(INOUT)  :: adjcnt(0:)
integer(kind=irg), INTENT(IN)     :: i
integer(kind=irg), INTENT(IN)     :: j

type(IO_T)                        :: Message
integer(kind=irg)                 :: p

do p = 1, adjcnt(i)
  if (adj(i,p).eq.j) return
end do

if (adjcnt(i).ge.icoMPmaxdegree) then
  call Message%printError('addNeighbour', 'neighbour list overflow')
end if

p = adjcnt(i)
do while (p.ge.1)
  if (adj(i,p).lt.j) exit
  adj(i,p+1) = adj(i,p)
  p = p - 1
end do
adj(i,p+1) = j
adjcnt(i) = adjcnt(i) + 1

end subroutine addNeighbour_

!--------------------------------------------------------------------------
recursive function cross_(a, b) result(c)
!DEC$ ATTRIBUTES DLLEXPORT :: cross_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! vector cross product

IMPLICIT NONE

real(kind=dbl), INTENT(IN)   :: a(3)
real(kind=dbl), INTENT(IN)   :: b(3)
real(kind=dbl)               :: c(3)

c(1) = a(2)*b(3) - a(3)*b(2)
c(2) = a(3)*b(1) - a(1)*b(3)
c(3) = a(1)*b(2) - a(2)*b(1)

end function cross_

!--------------------------------------------------------------------------
recursive subroutine cotangentLaplacian_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: cotangentLaplacian_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! COO representation (rows, cols, vals) of L = D - W, plus barycentric areas
!!
!! Duplicate (row,col) entries are kept rather than summed, exactly as the python
!! original emits them; scipy.sparse sums them on assembly.

use mod_memory

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self

type(memory_T)                  :: mem
integer(kind=irg)               :: nnz, t, f, ia, ib, ic, base1, base2, i
integer(kind=irg)               :: pa(3), pb(3), pc(3)
real(kind=dbl)                  :: u(3), w(3), x(3), crs, cot, val
real(kind=dbl), allocatable     :: diag(:)

mem = memory_T()

nnz = 6*self%nfaces + self%nverts

call mem%alloc(self%laprow, (/ nnz-1 /), 'laprow', startdims = (/ 0 /))
call mem%alloc(self%lapcol, (/ nnz-1 /), 'lapcol', startdims = (/ 0 /))
call mem%alloc(self%lapval, (/ nnz-1 /), 'lapval', startdims = (/ 0 /))
call mem%alloc(self%vertexarea, (/ self%nverts-1 /), 'vertexarea', initval = 0.D0, startdims = (/ 0 /))
call mem%alloc(diag, (/ self%nverts-1 /), 'diag', initval = 0.D0, startdims = (/ 0 /))

! the three cyclic permutations (i,j,k), (j,k,i), (k,i,j) of the face vertices;
! in each one the first entry is the apex at which the angle is measured
pa = (/ 1, 2, 3 /)
pb = (/ 2, 3, 1 /)
pc = (/ 3, 1, 2 /)

do t = 1, 3
  base1 = (2*(t-1)) * self%nfaces
  base2 = (2*(t-1) + 1) * self%nfaces
  do f = 0, self%nfaces-1
    ia = self%faces(f,pa(t))
    ib = self%faces(f,pb(t))
    ic = self%faces(f,pc(t))

    u = self%verts(ib,:) - self%verts(ia,:)
    w = self%verts(ic,:) - self%verts(ia,:)
    x = cross_(u, w)
    crs = dsqrt(sum(x*x))
    cot = dot_product(u, w) / max(crs, 1.0D-30)
    val = -0.5D0 * cot

    self%laprow(base1+f) = ib
    self%lapcol(base1+f) = ic
    self%lapval(base1+f) = val

    self%laprow(base2+f) = ic
    self%lapcol(base2+f) = ib
    self%lapval(base2+f) = val

    self%vertexarea(ia) = self%vertexarea(ia) + crs / 6.D0
  end do
end do

do i = 0, 6*self%nfaces-1
  diag(self%laprow(i)) = diag(self%laprow(i)) - self%lapval(i)
end do

do i = 0, self%nverts-1
  self%laprow(6*self%nfaces+i) = i
  self%lapcol(6*self%nfaces+i) = i
  self%lapval(6*self%nfaces+i) = diag(i)
end do

call mem%dealloc(diag, 'diag')

end subroutine cotangentLaplacian_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! master pattern input and energy selection
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine readMasterPattern_(self, EMsoft, HDF, mpfile, modality)
!DEC$ ATTRIBUTES DLLEXPORT :: readMasterPattern_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! read the two square Lambert hemispheres from an EBSD, ECP or TKD master file
!!
!! When modality is absent the file is interrogated for it.  The caller must have
!! opened the Fortran HDF interface before calling this routine.

use mod_EMsoft
use mod_HDFsupport
use mod_HDFnames
use mod_MPfiles
use mod_io
use stringconstants

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)             :: self
type(EMsoft_T), INTENT(INOUT)             :: EMsoft
type(HDF_T), INTENT(INOUT)                :: HDF
character(fnlen), INTENT(IN)              :: mpfile
character(*), INTENT(IN), OPTIONAL        :: modality

type(MPfile_T)                            :: MPFT
type(HDFnames_T)                          :: HDFnames
type(IO_T)                                :: Message
type(EBSDmasterNameListType)              :: ebsdnl
type(ECPmasterNameListType)               :: ecpnl
type(TKDmasterNameListType)               :: tkdnl
real(kind=sgl), allocatable               :: nh(:,:,:), sh(:,:,:), ev(:)
character(fnlen)                          :: fname, dgname
integer(kind=irg)                         :: nx

MPFT = MPfile_T()
HDFnames = HDFnames_T()

fname = EMsoft%generateFilePath('EMdatapathname', trim(mpfile))
call MPFT%setFileName(fname)

if (present(modality)) then
  call MPFT%setModality(trim(modality))
else
  call MPFT%determineModality(HDF, fname)
end if
self%modality = trim(MPFT%getModality())

select case (trim(self%modality))
  case ('EBSD')
    call HDFnames%set_ProgramData(SC_EBSDmaster)
    call HDFnames%set_NMLlist(SC_EBSDmasterNameList)
    call HDFnames%set_NMLfilename(SC_EBSDmasterNML)
    call HDFnames%set_Variable(SC_MCOpenCL)
    call MPFT%readMPfile(HDF, HDFnames, ebsdnl, getmLPNH = .TRUE., getmLPSH = .TRUE.)
    nx = ebsdnl%npx
    dgname = SC_EBSDmaster
  case ('TKD')
    call HDFnames%set_ProgramData(SC_TKDmaster)
    call HDFnames%set_NMLlist(SC_TKDmasterNameList)
    call HDFnames%set_NMLfilename(SC_TKDmasterNML)
    call HDFnames%set_Variable(SC_MCOpenCL)
    call MPFT%readMPfile(HDF, HDFnames, tkdnl, getmLPNH = .TRUE., getmLPSH = .TRUE.)
    nx = tkdnl%npx
    dgname = SC_TKDmaster
  case ('ECP')
    call HDFnames%set_ProgramData(SC_ECPmaster)
    call HDFnames%set_NMLlist(SC_ECPmasterNameList)
    call HDFnames%set_NMLfilename(SC_ECPmasterNML)
    call HDFnames%set_Variable(SC_MCOpenCL)
    call MPFT%readMPfile(HDF, HDFnames, ecpnl, getmLPNH = .TRUE., getmLPSH = .TRUE.)
    nx = ecpnl%npx
    dgname = SC_ECPmaster
  case default
    call Message%printError('readMasterPattern', 'unsupported master pattern modality '//trim(self%modality))
end select

call MPFT%copymLPNH(nh)
call MPFT%copymLPSH(sh)

! the energy list is not read through readMPfile, because that routine looks for
! a dataset named keVs whereas the master pattern programs write EkeVs
call readEnergies_(HDF, fname, dgname, ev)

if (allocated(ev)) then
  call self%setMasterPattern_(nx, nh, sh, ev)
else
  call self%setMasterPattern_(nx, nh, sh)
end if

self%masterfile = trim(fname)

deallocate(nh, sh)
if (allocated(ev)) deallocate(ev)

end subroutine readMasterPattern_

!--------------------------------------------------------------------------
recursive subroutine readEnergies_(HDF, fname, datagroupname, ev)
!DEC$ ATTRIBUTES DLLEXPORT :: readEnergies_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! read the energy bin centers from a master pattern file, if they are there
!!
!! The master pattern programs store this array as EkeVs; some other file types
!! call it keVs, so both names are tried.  ev is left unallocated when neither is
!! present, which is a normal outcome for single-energy modalities.

use HDF5
use mod_HDFsupport
use stringconstants

IMPLICIT NONE

type(HDF_T), INTENT(INOUT)                  :: HDF
character(fnlen), INTENT(IN)                :: fname
character(fnlen), INTENT(IN)                :: datagroupname
real(kind=sgl), allocatable, INTENT(INOUT)  :: ev(:)

character(fnlen)                            :: groupname, dataset
logical                                     :: g_exists
integer(kind=irg)                           :: hdferr
integer(HSIZE_T)                            :: dims(1)

if (allocated(ev)) deallocate(ev)

hdferr = HDF%openFile(fname, readonly = .TRUE.)

groupname = SC_EMData
hdferr = HDF%openGroup(groupname)
groupname = trim(datagroupname)
hdferr = HDF%openGroup(groupname)

dataset = SC_EkeVs
call H5Lexists_f(HDF%getobjectID(), trim(dataset), g_exists, hdferr)
if (.not.g_exists) then
  dataset = SC_keVs
  call H5Lexists_f(HDF%getobjectID(), trim(dataset), g_exists, hdferr)
end if

if (g_exists) call HDF%readDatasetFloatArray(dataset, dims, hdferr, ev)

call HDF%popall()

end subroutine readEnergies_

!--------------------------------------------------------------------------
recursive subroutine setMasterPattern_(self, npx, nh, sh, keVs)
!DEC$ ATTRIBUTES DLLEXPORT :: setMasterPattern_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! set the master pattern arrays directly, bypassing the HDF5 reader
!!
!! nh and sh must be dimensioned (-npx:npx,-npx:npx,numEbins).  This entry point
!! exists so that the module can be exercised without an HDF5 file.

use mod_io
use mod_memory

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)             :: self
integer(kind=irg), INTENT(IN)             :: npx
real(kind=sgl), INTENT(IN)                :: nh(-npx:,-npx:,:)
real(kind=sgl), INTENT(IN)                :: sh(-npx:,-npx:,:)
real(kind=sgl), INTENT(IN), OPTIONAL      :: keVs(:)

type(IO_T)                                :: Message
type(memory_T)                            :: mem
integer(kind=irg)                         :: s(3)

mem = memory_T()

s = shape(nh)
if ((s(1).ne.2*npx+1).or.(s(2).ne.2*npx+1)) then
  call Message%printError('setMasterPattern', 'master pattern arrays must be square and of size (2*npx+1)^2')
end if
if (any(shape(sh).ne.s)) then
  call Message%printError('setMasterPattern', 'the two hemispheres have different shapes')
end if

self%npx = npx
self%numEbins = s(3)

call mem%alloc(self%mLPNH, (/ npx, npx, s(3) /), 'mLPNH', startdims = (/ -npx, -npx, 1 /))
call mem%alloc(self%mLPSH, (/ npx, npx, s(3) /), 'mLPSH', startdims = (/ -npx, -npx, 1 /))

self%mLPNH = dble(nh)
self%mLPSH = dble(sh)

if (present(keVs)) then
  if (size(keVs).eq.s(3)) then
    call mem%alloc(self%keVs, (/ s(3) /), 'keVs')
    self%keVs = dble(keVs)
    self%haskeVs = .TRUE.
  else
    self%haskeVs = .FALSE.
  end if
else
  self%haskeVs = .FALSE.
end if

self%hasmaster = .TRUE.
self%eqmismatch = self%equatorMismatch_()

end subroutine setMasterPattern_

!--------------------------------------------------------------------------
recursive function equatorMismatch_(self) result(m)
!DEC$ ATTRIBUTES DLLEXPORT :: equatorMismatch_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! largest |NH - SH| on the square border; this should be about zero

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
real(kind=dbl)                  :: m

integer(kind=irg)               :: n

n = self%npx

m = maxval(dabs(self%mLPNH(:,-n,:) - self%mLPSH(:,-n,:)))
m = max(m, maxval(dabs(self%mLPNH(:, n,:) - self%mLPSH(:, n,:))))
m = max(m, maxval(dabs(self%mLPNH(-n,:,:) - self%mLPSH(-n,:,:))))
m = max(m, maxval(dabs(self%mLPNH( n,:,:) - self%mLPSH( n,:,:))))

end function equatorMismatch_

!--------------------------------------------------------------------------
recursive function pixelDegrees_(self) result(p)
!DEC$ ATTRIBUTES DLLEXPORT :: pixelDegrees_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! angular size of one master pattern pixel, in degrees

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
real(kind=dbl)                  :: p

p = dsqrt(2.D0*cPi/dble(2*self%npx+1)**2) * 180.D0/cPi

end function pixelDegrees_

!--------------------------------------------------------------------------
recursive subroutine selectEnergy_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: selectEnergy_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! collapse the energy axis and store the result with the channel axis leading
!!
!! The channel axis is moved to the front so that the innermost loop of the
!! interpolation runs over contiguous memory.

use, intrinsic :: ieee_arithmetic
use mod_io
use mod_memory

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self

type(IO_T)                      :: Message
type(memory_T)                  :: mem
integer(kind=irg)               :: n, c, ix, iy, ib
real(kind=dbl)                  :: sw

mem = memory_T()
n = self%npx

select case (trim(self%energymode))
  case ('all')
    self%nchannels = self%numEbins
  case default
    self%nchannels = 1
end select

call mem%alloc(self%nhsel, (/ self%nchannels, n, n /), 'nhsel', startdims = (/ 1, -n, -n /))
call mem%alloc(self%shsel, (/ self%nchannels, n, n /), 'shsel', startdims = (/ 1, -n, -n /))
call mem%alloc(self%energies, (/ self%nchannels /), 'energies')

select case (trim(self%energymode))
  case ('all')
    do iy = -n, n
      do ix = -n, n
        do c = 1, self%nchannels
          self%nhsel(c,ix,iy) = self%mLPNH(ix,iy,c)
          self%shsel(c,ix,iy) = self%mLPSH(ix,iy,c)
        end do
      end do
    end do
    if (self%haskeVs.eqv..TRUE.) then
      self%energies = self%keVs
    else
      do c = 1, self%nchannels
        self%energies(c) = dble(c-1)
      end do
    end if

  case ('sum')
    do iy = -n, n
      do ix = -n, n
        self%nhsel(1,ix,iy) = sum(self%mLPNH(ix,iy,:))
        self%shsel(1,ix,iy) = sum(self%mLPSH(ix,iy,:))
      end do
    end do
    self%energies(1) = ieee_value(1.D0, ieee_quiet_nan)

  case ('weighted')
    if (size(self%eweights).ne.self%numEbins) then
      call Message%printError('selectEnergy', 'the number of energy weights does not match the number of energy bins')
    end if
    do iy = -n, n
      do ix = -n, n
        self%nhsel(1,ix,iy) = sum(self%eweights * self%mLPNH(ix,iy,:))
        self%shsel(1,ix,iy) = sum(self%eweights * self%mLPSH(ix,iy,:))
      end do
    end do
    if (self%haskeVs.eqv..TRUE.) then
      sw = sum(self%eweights)
      self%energies(1) = sum(self%eweights * self%keVs) / sw
    else
      self%energies(1) = ieee_value(1.D0, ieee_quiet_nan)
    end if

  case ('index')
    ib = self%energybin + 1
    if ((ib.lt.1).or.(ib.gt.self%numEbins)) then
      call Message%printError('selectEnergy', 'the requested energy bin lies outside the available range')
    end if
    do iy = -n, n
      do ix = -n, n
        self%nhsel(1,ix,iy) = self%mLPNH(ix,iy,ib)
        self%shsel(1,ix,iy) = self%mLPSH(ix,iy,ib)
      end do
    end do
    if (self%haskeVs.eqv..TRUE.) then
      self%energies(1) = self%keVs(ib)
    else
      self%energies(1) = ieee_value(1.D0, ieee_quiet_nan)
    end if
end select

end subroutine selectEnergy_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! sampling of the master pattern
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine chooseKernel_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: chooseKernel_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! bicubic only when the master pattern is barely finer than the vertex cells
!!
!! 63.435 degrees is the mean vertex spacing of the base icosahedron; it halves
!! with every subdivision level.

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self

real(kind=dbl)                  :: spacing

if (trim(self%interpmode).eq.'auto') then
  spacing = 63.435D0 / 2.D0**self%level
  if (spacing/self%pixelDegrees_() .lt. 4.D0) then
    self%kernel = 'bicubic'
  else
    self%kernel = 'bilinear'
  end if
else
  self%kernel = trim(self%interpmode)
end if

end subroutine chooseKernel_

!--------------------------------------------------------------------------
recursive subroutine cubicTaps_(t, w)
!DEC$ ATTRIBUTES DLLEXPORT :: cubicTaps_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! the four Catmull-Rom weights for a fractional offset t

IMPLICIT NONE

real(kind=dbl), INTENT(IN)    :: t
real(kind=dbl), INTENT(OUT)   :: w(4)

real(kind=dbl)                :: t2, t3

t2 = t*t
t3 = t2*t

w(1) = -0.5D0*t3 + t2 - 0.5D0*t
w(2) =  1.5D0*t3 - 2.5D0*t2 + 1.D0
w(3) = -1.5D0*t3 + 2.D0*t2 + 0.5D0*t
w(4) =  0.5D0*t3 - 0.5D0*t2

end subroutine cubicTaps_

!--------------------------------------------------------------------------
recursive subroutine sampleDirections_(self, nd, dirs, out)
!DEC$ ATTRIBUTES DLLEXPORT :: sampleDirections_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! interpolate the selected master pattern at a list of unit directions
!!
!! Directions with z >= 0 read the northern hemisphere, the others read the
!! southern one at |z|.  Both kernels are interpolating, which is what keeps the
!! equator seam exact.  selectEnergy_ must have been called first.
!!
!! This routine only reads from self and is safe to call from several threads.

use mod_Lambert

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
integer(kind=irg), INTENT(IN)   :: nd
real(kind=dbl), INTENT(IN)      :: dirs(nd,3)
real(kind=dbl), INTENT(OUT)     :: out(self%nchannels,nd)

type(Lambert_T)                 :: L
integer(kind=irg)               :: d, n, ix, iy, jx, jy, xx, yy, ntap, ierr, offs(4)
real(kind=dbl)                  :: dc(3), xy(2), px, py, fx, fy, wx(4), wy(4), wgt
real(kind=dbl)                  :: acc(self%nchannels)
logical                         :: north

n = self%npx
L = Lambert_T()

if (trim(self%kernel).eq.'bicubic') then
  ntap = 4
  offs = (/ -1, 0, 1, 2 /)
else
  ntap = 2
  offs = (/ 0, 1, 0, 0 /)
end if

do d = 1, nd
  dc = dirs(d,:)
  north = (dc(3).ge.0.D0)
  dc(3) = dabs(dc(3))

  call L%setxyzd(dc)
  ierr = L%LambertSphereToSquare(xy)

  px = (xy(1) + 1.D0) * dble(n)
  py = (xy(2) + 1.D0) * dble(n)
  ix = min(max(floor(px), 0), 2*n-1)
  iy = min(max(floor(py), 0), 2*n-1)
  fx = px - dble(ix)
  fy = py - dble(iy)

  if (ntap.eq.2) then
    wx(1) = 1.D0 - fx
    wx(2) = fx
    wy(1) = 1.D0 - fy
    wy(2) = fy
  else
    call cubicTaps_(fx, wx)
    call cubicTaps_(fy, wy)
  end if

  acc = 0.D0
  do jy = 1, ntap
    yy = min(max(iy+offs(jy), 0), 2*n) - n
    do jx = 1, ntap
      xx = min(max(ix+offs(jx), 0), 2*n) - n
      wgt = wx(jx) * wy(jy)
      if (north.eqv..TRUE.) then
        acc = acc + self%nhsel(:,xx,yy) * wgt
      else
        acc = acc + self%shsel(:,xx,yy) * wgt
      end if
    end do
  end do
  out(:,d) = acc
end do

end subroutine sampleDirections_

!--------------------------------------------------------------------------
recursive subroutine discQuadrature_(nrings, nazim, quad)
!DEC$ ATTRIBUTES DLLEXPORT :: discQuadrature_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! equal-area unit-disc quadrature points (rho, phi), all with the same weight
!!
!! The rings are staggered by the golden angle so that the points do not line up
!! radially.

IMPLICIT NONE

integer(kind=irg), INTENT(IN)   :: nrings
integer(kind=irg), INTENT(IN)   :: nazim
real(kind=dbl), INTENT(OUT)     :: quad(nrings*nazim,2)

integer(kind=irg)               :: r, a, k
real(kind=dbl)                  :: rho, golden

golden = cPi * (3.D0 - dsqrt(5.D0))

k = 0
do r = 0, nrings-1
  rho = dsqrt((dble(r) + 0.5D0)/dble(nrings))
  do a = 0, nazim-1
    k = k + 1
    quad(k,1) = rho
    quad(k,2) = 2.D0*cPi*dble(a)/dble(nazim) + golden*dble(r)
  end do
end do

end subroutine discQuadrature_

!--------------------------------------------------------------------------
recursive subroutine sampleCellAverage_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: sampleCellAverage_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! average the master pattern over each vertex cell instead of point sampling it
!!
!! The cell is approximated by a spherical cap of radius spacing/sqrt(3), the
!! circumradius of a hexagon.  Each vertex gets its own azimuthal offset, so that
!! whatever aliasing is left looks like noise rather than structure.

use mod_memory
use omp_lib

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self

type(memory_T)                  :: mem
real(kind=dbl), allocatable     :: quad(:,:), radius(:), stagger(:)
real(kind=dbl), allocatable     :: dirs(:,:), s(:,:)
integer(kind=irg)               :: nq, i, k, col, deg
real(kind=dbl)                  :: golden, twopi, dsum, dots
real(kind=dbl)                  :: v(3), seed(3), e1(3), e2(3), tang(3), th, ph

mem = memory_T()

nq = self%aarings * self%aaazim
golden = cPi * (3.D0 - dsqrt(5.D0))
twopi = 2.D0 * cPi

call mem%alloc(quad, (/ nq, 2 /), 'quad')
call discQuadrature_(self%aarings, self%aaazim, quad)

call mem%alloc(radius, (/ self%nverts-1 /), 'radius', startdims = (/ 0 /))
call mem%alloc(stagger, (/ self%nverts-1 /), 'stagger', startdims = (/ 0 /))

do i = 0, self%nverts-1
  dsum = 0.D0
  deg = 0
  do col = 1, 6
    if (self%neigh(i,col).ge.0) then
      dots = dot_product(self%verts(i,:), self%verts(self%neigh(i,col),:))
      dsum = dsum + dacos(min(max(dots, -1.D0), 1.D0))
      deg = deg + 1
    end if
  end do
  radius(i) = dsum / dble(deg) / dsqrt(3.D0)
  stagger(i) = modulo(dble(i) * golden, twopi)
end do

call mem%alloc(self%signal, (/ self%nchannels, self%nverts-1 /), 'signal', startdims = (/ 1, 0 /))

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, k, v, seed, e1, e2, tang, th, ph, dirs, s)
allocate(dirs(nq,3), s(self%nchannels,nq))

!$OMP DO SCHEDULE(DYNAMIC,256)
do i = 0, self%nverts-1
  v = self%verts(i,:)
! avoid a null cross product near the poles
  if (dabs(v(3)).gt.0.9D0) then
    seed = (/ 1.D0, 0.D0, 0.D0 /)
  else
    seed = (/ 0.D0, 0.D0, 1.D0 /)
  end if
  e1 = cross_(v, seed)
  e1 = e1 / dsqrt(sum(e1*e1))
  e2 = cross_(v, e1)

  do k = 1, nq
    th = radius(i) * quad(k,1)
    ph = quad(k,2) + stagger(i)
    tang = dcos(ph)*e1 + dsin(ph)*e2
    dirs(k,:) = dcos(th)*v + dsin(th)*tang
  end do

  call self%sampleDirections_(nq, dirs, s)
  self%signal(:,i) = sum(s, 2) / dble(nq)
end do
!$OMP END DO

deallocate(dirs, s)
!$OMP END PARALLEL

call mem%dealloc(quad, 'quad')
call mem%dealloc(radius, 'radius')
call mem%dealloc(stagger, 'stagger')

end subroutine sampleCellAverage_

!--------------------------------------------------------------------------
recursive subroutine normalizeSignal_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: normalizeSignal_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! per-channel normalization of the vertex signal

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self

integer(kind=irg)               :: c, nv
real(kind=dbl)                  :: mu, sd, lo, hi

nv = self%nverts

select case (trim(self%normmode))
  case ('none')
    return

  case ('zscore')
    do c = 1, self%nchannels
      mu = sum(self%signal(c,:)) / dble(nv)
      sd = dsqrt(sum((self%signal(c,:) - mu)**2) / dble(nv))
      if (sd.le.0.D0) sd = 1.D0
      self%signal(c,:) = (self%signal(c,:) - mu) / sd
    end do

  case ('minmax')
    do c = 1, self%nchannels
      lo = minval(self%signal(c,:))
      hi = maxval(self%signal(c,:))
      if (hi.gt.lo) then
        self%signal(c,:) = (self%signal(c,:) - lo) / (hi - lo)
      else
        self%signal(c,:) = self%signal(c,:) - lo
      end if
    end do

  case ('mean1')
    do c = 1, self%nchannels
      mu = sum(self%signal(c,:)) / dble(nv)
      if (dabs(mu).gt.0.D0) self%signal(c,:) = self%signal(c,:) / mu
    end do
end select

end subroutine normalizeSignal_

!--------------------------------------------------------------------------
recursive subroutine projectMaster_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: projectMaster_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! the full pipeline: icosphere, energy selection, sampling, normalization and
!! the derived chart representation

use mod_io
use mod_memory

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self

type(IO_T)                      :: Message
type(memory_T)                  :: mem
integer(kind=irg), allocatable  :: used(:)
integer(kind=irg)               :: i, r, c, hh, ww, np, k

mem = memory_T()

if (self%hasmaster.eqv..FALSE.) then
  call Message%printError('projectMaster', 'no master pattern has been read')
end if

if (self%nverts.eq.0) call self%buildIcosphere_()

call self%selectEnergy_()
call self%chooseKernel_()

if (self%antialias.eqv..TRUE.) then
  call self%sampleCellAverage_()
else
  call mem%alloc(self%signal, (/ self%nchannels, self%nverts-1 /), 'signal', startdims = (/ 1, 0 /))
  call self%sampleDirections_(self%nverts, self%verts, self%signal)
end if

call self%normalizeSignal_()

! the charts include the border shared with the next chart; dropping the first
! row and the last column leaves a disjoint tiling of everything but the poles
hh = self%hchart - 1
ww = self%wchart - 1

call mem%alloc(self%chartsignal, (/ self%nchannels, 5, hh, ww /), 'chartsignal')
call mem%alloc(used, (/ self%nverts-1 /), 'used', initval = 0, startdims = (/ 0 /))

do i = 0, 4
  do r = 1, hh
    do c = 1, ww
      k = self%chrt(i,r,c-1)
      self%chartsignal(:,i+1,r,c) = self%signal(:,k)
      used(k) = used(k) + 1
    end do
  end do
end do

np = count(used.eq.0)
call mem%alloc(self%poleindex, (/ np /), 'poleindex')
k = 0
do i = 0, self%nverts-1
  if (used(i).eq.0) then
    k = k + 1
    self%poleindex(k) = i
  end if
end do
call mem%dealloc(used, 'used')

if (self%dolaplacian.eqv..TRUE.) call self%cotangentLaplacian_()

end subroutine projectMaster_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! output
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive function i2s_(i) result(s)
!DEC$ ATTRIBUTES DLLEXPORT :: i2s_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! shortest decimal representation of an integer

IMPLICIT NONE

integer(kind=irg), INTENT(IN)   :: i
character(len=:), allocatable   :: s

character(25)                   :: line

write (line,"(I0)") i
s = trim(line)

end function i2s_

!--------------------------------------------------------------------------
recursive function d2s_(x, frm) result(s)
!DEC$ ATTRIBUTES DLLEXPORT :: d2s_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! decimal representation of a double, for the metadata string

IMPLICIT NONE

real(kind=dbl), INTENT(IN)            :: x
character(*), INTENT(IN), OPTIONAL    :: frm
character(len=:), allocatable         :: s

character(40)                         :: line

if (present(frm)) then
  write (line,frm) x
else
  write (line,"(ES23.15E3)") x
end if
s = trim(adjustl(line))

end function d2s_

!--------------------------------------------------------------------------
recursive function metaJSON_(self) result(js)
!DEC$ ATTRIBUTES DLLEXPORT :: metaJSON_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! assemble the metadata string that is stored as meta_json in the .npz file
!!
!! The keys and their order match the python original; the numeric formatting is
!! Fortran's, so the string is not byte-identical to what json.dumps produces.

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
character(len=:), allocatable   :: js

character(len=:), allocatable   :: nl, emode, wts
integer(kind=irg)               :: i

nl = achar(10)

if (trim(self%energymode).eq.'index') then
  emode = i2s_(self%energybin)
else
  emode = trim(self%energymode)
end if

if (allocated(self%eweights).and.(trim(self%energymode).eq.'weighted')) then
  wts = '['
  do i = 1, size(self%eweights)
    wts = wts//d2s_(self%eweights(i))
    if (i.lt.size(self%eweights)) wts = wts//', '
  end do
  wts = wts//']'
else
  wts = 'null'
end if

js = '{'//nl// &
     ' "generator": "mod_icoMP.f90 v'//icoMPversion//'",'//nl// &
     ' "source": "'//trim(self%masterfile)//'",'//nl// &
     ' "master_npx": '//i2s_(self%npx)//','//nl// &
     ' "master_n_energy_bins": '//i2s_(self%numEbins)//','//nl// &
     ' "equator_mismatch": '//d2s_(self%eqmismatch)//','//nl// &
     ' "level": '//i2s_(self%level)//','//nl// &
     ' "n_vertices": '//i2s_(self%nverts)//','//nl// &
     ' "n_faces": '//i2s_(self%nfaces)//','//nl// &
     ' "chart_shape": [5, '//i2s_(self%hchart-1)//', '//i2s_(self%wchart-1)//'],'//nl// &
     ' "energy_mode": "'//emode//'",'//nl// &
     ' "energy_weights": '//wts//','//nl// &
     ' "normalize": "'//trim(self%normmode)//'",'//nl

if (self%antialias.eqv..TRUE.) then
  js = js//' "sampling": "cell average, '//i2s_(self%aarings)//'x'//i2s_(self%aaazim)// &
       ' equal-area disc quadrature",'//nl
else
  js = js//' "sampling": "point sample",'//nl
end if

js = js// &
     ' "interpolation": "'//trim(self%kernel)//'",'//nl// &
     ' "master_pixels_per_cell": '// &
     d2s_(63.435D0/2.D0**self%level/self%pixelDegrees_(), "(F0.2)")//','//nl// &
     ' "convention": "square Lambert (EMsoft); z>=0 from mLPNH, z<0 from mLPSH at |z|",'//nl// &
     ' "chart_slice": "charts[:, 1:, :-1]; the 2 poles are excluded"'//nl// &
     '}'

end function metaJSON_

!--------------------------------------------------------------------------
recursive subroutine writeNPZ_(self, npzfile)
!DEC$ ATTRIBUTES DLLEXPORT :: writeNPZ_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! write all output arrays to a numpy .npz file
!!
!! Note that add_npz drives the external zip program and drops a temporary .npy
!! file in the current working directory for each array, so zip must be on the
!! path and two runs must not share a working directory.  The entries are stored
!! rather than deflated, so the file is larger than np.savez_compressed output
!! but reads back identically.

use mod_io
use mod_memory
use mod_npz

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
character(fnlen), INTENT(IN)    :: npzfile

type(IO_T)                      :: Message
type(memory_T)                  :: mem
real(kind=sgl), allocatable     :: s2(:,:), s4(:,:,:,:)
character(len=:), allocatable   :: js
character(fnlen)                :: fname
logical                         :: f_exists
integer(kind=irg)               :: iunit

mem = memory_T()
fname = trim(npzfile)

! zip appends to an existing archive, so an old file has to go first
inquire(file=trim(fname), exist=f_exists)
if (f_exists.eqv..TRUE.) then
  open(newunit=iunit, file=trim(fname), status='old')
  close(unit=iunit, status='delete')
end if

if (self%doubleprec.eqv..TRUE.) then
  call add_npz(trim(fname), 'signal', self%signal)
  call add_npz(trim(fname), 'chart_signal', self%chartsignal)
  call add_npz(trim(fname), 'vertices', self%verts)
else
  call mem%alloc(s2, (/ self%nchannels, self%nverts /), 's2')
  s2 = real(self%signal, sgl)
  call add_npz(trim(fname), 'signal', s2)
  call mem%dealloc(s2, 's2')

  call mem%alloc(s4, (/ self%nchannels, 5, self%hchart-1, self%wchart-1 /), 's4')
  s4 = real(self%chartsignal, sgl)
  call add_npz(trim(fname), 'chart_signal', s4)
  call mem%dealloc(s4, 's4')

  call mem%alloc(s2, (/ self%nverts, 3 /), 's2')
  s2 = real(self%verts, sgl)
  call add_npz(trim(fname), 'vertices', s2)
  call mem%dealloc(s2, 's2')
end if

call add_npz(trim(fname), 'faces', self%faces)
call add_npz(trim(fname), 'neighbours', self%neigh)
call add_npz(trim(fname), 'charts', self%chrt)
call add_npz(trim(fname), 'level_vertex_counts', self%counts)
call add_npz(trim(fname), 'pole_index', self%poleindex)
call add_npz(trim(fname), 'energies_keV', self%energies)
call add_npz(trim(fname), 'level', self%level)

js = self%metaJSON_()
call add_npz(trim(fname), 'meta_json', js)

if (self%dolaplacian.eqv..TRUE.) then
  call add_npz(trim(fname), 'laplacian_row', self%laprow)
  call add_npz(trim(fname), 'laplacian_col', self%lapcol)
  call add_npz(trim(fname), 'laplacian_val', self%lapval)
  call add_npz(trim(fname), 'vertex_area', self%vertexarea)
end if

call Message%printMessage(' -> wrote '//trim(fname))

end subroutine writeNPZ_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! verification
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive function edgeExists_(adj, adjcnt, i, j) result(yes)
!DEC$ ATTRIBUTES DLLEXPORT :: edgeExists_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! is (i,j) an edge of the adjacency structure?

IMPLICIT NONE

integer(kind=irg), INTENT(IN)   :: adj(0:,:)
integer(kind=irg), INTENT(IN)   :: adjcnt(0:)
integer(kind=irg), INTENT(IN)   :: i
integer(kind=irg), INTENT(IN)   :: j
logical                         :: yes

integer(kind=irg)               :: p

yes = .FALSE.
do p = 1, adjcnt(i)
  if (adj(i,p).eq.j) then
    yes = .TRUE.
    return
  end if
end do

end function edgeExists_

!--------------------------------------------------------------------------
recursive subroutine selfTest_(self, maxlevel)
!DEC$ ATTRIBUTES DLLEXPORT :: selfTest_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! geometry self tests; these need no master pattern
!!
!! This is the Fortran counterpart of the python --self-test option.

use mod_io
use mod_memory
use mod_Lambert

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)             :: self
integer(kind=irg), INTENT(IN), OPTIONAL   :: maxlevel

type(IO_T)                                :: Message
type(memory_T)                            :: mem
type(icoMP_T)                             :: t, tprev
type(Lambert_T)                           :: L
integer(kind=irg), allocatable            :: fadj(:,:), fcnt(:), used(:)
real(kind=dbl), allocatable               :: LM(:,:)
integer(kind=irg)                         :: mx, lvl, i, j, k, p, r, c, nv, nf, deg, degn, n5, n6
integer(kind=irg)                         :: ns(icoMPmaxdegree), ms(icoMPmaxdegree), nEf, nEn, hh, ww
integer(kind=irg)                         :: ierr, ntest
real(kind=dbl)                            :: dmax, angmin, angmax, ang, dots, golden, z, th, ph
real(kind=dbl)                            :: dc(3), dc2(3), xy(2), asum
character(fnlen)                          :: charline

mem = memory_T()
mx = 4
if (present(maxlevel)) mx = maxlevel

call Message%printMessage(' ')
call Message%printMessage(' mod_icoMP self test')
call Message%printMessage(' -------------------')

!--------------------------------------------------------------------------
! Lambert round trip
!--------------------------------------------------------------------------
ntest = 20000
golden = cPi * (3.D0 - dsqrt(5.D0))
L = Lambert_T()
dmax = 0.D0
do i = 1, ntest
  z = 2.D0*(dble(i)-0.5D0)/dble(ntest) - 1.D0
  th = dacos(min(max(z, -1.D0), 1.D0))
  ph = golden * dble(i)
  dc = (/ dsin(th)*dcos(ph), dsin(th)*dsin(ph), dabs(z) /)
  dc = dc / dsqrt(sum(dc*dc))
  call L%setxyzd(dc)
  ierr = L%LambertSphereToSquare(xy)
  call L%setxyd(xy)
  ierr = L%LambertSquareToSphere(dc2)
  dmax = max(dmax, maxval(dabs(dc2 - dc)))
end do
write (charline,"(' Lambert round trip: max err ',ES12.4)") dmax
call Message%printMessage(trim(charline))
if (dmax.ge.1.0D-12) call Message%printError('selfTest', 'Lambert round trip error is too large')

!--------------------------------------------------------------------------
! icosphere geometry, level by level
!--------------------------------------------------------------------------
do lvl = 0, mx
  t = icoMP_T( lvl )
  call t%buildIcosphere_()
  nv = t%nverts
  nf = t%nfaces
  hh = t%hchart
  ww = t%wchart

  if (nv.ne.10*4**lvl+2) call Message%printError('selfTest', 'wrong number of vertices')
  if (nf.ne.20*4**lvl) call Message%printError('selfTest', 'wrong number of faces')

! all vertices on the unit sphere
  dmax = 0.D0
  do i = 0, nv-1
    dmax = max(dmax, dabs(dsqrt(sum(t%verts(i,:)**2)) - 1.D0))
  end do
  if (dmax.ge.1.0D-12) call Message%printError('selfTest', 'vertices are not on the unit sphere')

! the adjacency implied by the face list
  call mem%alloc(fadj, (/ nv-1, icoMPmaxdegree /), 'fadj', startdims = (/ 0, 1 /))
  call mem%alloc(fcnt, (/ nv-1 /), 'fcnt', initval = 0, startdims = (/ 0 /))
  do i = 0, nf-1
    call addNeighbour_(fadj, fcnt, t%faces(i,1), t%faces(i,2))
    call addNeighbour_(fadj, fcnt, t%faces(i,1), t%faces(i,3))
    call addNeighbour_(fadj, fcnt, t%faces(i,2), t%faces(i,1))
    call addNeighbour_(fadj, fcnt, t%faces(i,2), t%faces(i,3))
    call addNeighbour_(fadj, fcnt, t%faces(i,3), t%faces(i,1))
    call addNeighbour_(fadj, fcnt, t%faces(i,3), t%faces(i,2))
  end do

! twelve pentagons and the rest hexagons, and the neighbour table must describe
! exactly the same edges as the face list
  n5 = 0
  n6 = 0
  nEf = 0
  nEn = 0
  do i = 0, nv-1
    deg = fcnt(i)
    degn = count(t%neigh(i,:).ge.0)
    if (deg.ne.degn) call Message%printError('selfTest', 'neighbour table and face list disagree on the degree')
    if (deg.eq.5) n5 = n5 + 1
    if (deg.eq.6) n6 = n6 + 1
! sort the CCW neighbours so that the two sets can be compared entry by entry
    ns(1:deg) = t%neigh(i,1:deg)
    do p = 2, deg
      k = ns(p)
      j = p - 1
      do while (j.ge.1)
        if (ns(j).le.k) exit
        ns(j+1) = ns(j)
        j = j - 1
      end do
      ns(j+1) = k
    end do
    ms(1:deg) = fadj(i,1:deg)
    if (any(ns(1:deg).ne.ms(1:deg))) then
      call Message%printError('selfTest', 'neighbour table and face list describe different edges')
    end if
    do p = 1, deg
      if (ns(p).gt.i) nEn = nEn + 1
      if (fadj(i,p).gt.i) nEf = nEf + 1
    end do
  end do
  if (n5.ne.12) call Message%printError('selfTest', 'there should be exactly twelve degree-five vertices')
  if (n6.ne.nv-12) call Message%printError('selfTest', 'all other vertices should have degree six')
  if (nEn.ne.nEf) call Message%printError('selfTest', 'the two edge sets have different sizes')

! every row, column and anti-diagonal step of a chart must be a real edge
  do i = 0, 4
    do r = 0, hh-1
      do c = 0, ww-2
        if (.not.edgeExists_(fadj, fcnt, t%chrt(i,r,c), t%chrt(i,r,c+1))) then
          call Message%printError('selfTest', 'a chart row step is not an edge')
        end if
      end do
    end do
    do r = 0, hh-2
      do c = 0, ww-1
        if (.not.edgeExists_(fadj, fcnt, t%chrt(i,r,c), t%chrt(i,r+1,c))) then
          call Message%printError('selfTest', 'a chart column step is not an edge')
        end if
      end do
    end do
    do r = 0, hh-2
      do c = 0, ww-2
        if (.not.edgeExists_(fadj, fcnt, t%chrt(i,r,c+1), t%chrt(i,r+1,c))) then
          call Message%printError('selfTest', 'a chart anti-diagonal step is not an edge')
        end if
      end do
    end do
  end do

! the charts tile everything except the two poles, exactly once
  call mem%alloc(used, (/ nv-1 /), 'used', initval = 0, startdims = (/ 0 /))
  do i = 0, 4
    do r = 1, hh-1
      do c = 0, ww-2
        used(t%chrt(i,r,c)) = used(t%chrt(i,r,c)) + 1
      end do
    end do
  end do
  if (maxval(used).ne.1) call Message%printError('selfTest', 'the charts do not tile the sphere exactly once')
  if (count(used.eq.0).ne.2) call Message%printError('selfTest', 'the charts should leave exactly two vertices out')
  if ((used(0).ne.0).or.(used(11).ne.0)) then
    call Message%printError('selfTest', 'the two vertices left out should be the poles 0 and 11')
  end if
  call mem%dealloc(used, 'used')

! coarse levels are prefixes of fine ones, so pooling is a slice
  if (lvl.gt.0) then
    tprev = icoMP_T( lvl-1 )
    call tprev%buildIcosphere_()
    if (t%counts(lvl-1).ne.tprev%nverts) then
      call Message%printError('selfTest', 'the level vertex counts are inconsistent')
    end if
    dmax = 0.D0
    do i = 0, tprev%nverts-1
      dmax = max(dmax, maxval(dabs(t%verts(i,:) - tprev%verts(i,:))))
    end do
    if (dmax.ge.1.0D-12) call Message%printError('selfTest', 'the coarse level is not a prefix of the fine one')
  end if

! edge lengths
  angmin = 1.0D30
  angmax = -1.0D30
  do i = 0, nf-1
    dots = dot_product(t%verts(t%faces(i,1),:), t%verts(t%faces(i,2),:))
    ang = dacos(min(max(dots, -1.D0), 1.D0)) * 180.D0/cPi
    angmin = min(angmin, ang)
    angmax = max(angmax, ang)
    if (lvl.eq.0) then
      if (dabs(dots - 1.D0/dsqrt(5.D0)).ge.1.0D-12) then
        call Message%printError('selfTest', 'the base icosahedron edges are not all equal')
      end if
    end if
  end do

  write (charline,"(' level ',I1,': ',I7,' verts, charts (',I4,',',I4,'), edge ',F7.3,'-',F7.3,' deg  OK')") &
        lvl, nv, hh-1, ww-1, angmin, angmax
  call Message%printMessage(trim(charline))

  call mem%dealloc(fadj, 'fadj')
  call mem%dealloc(fcnt, 'fcnt')
end do

!--------------------------------------------------------------------------
! the cotangent Laplacian
!--------------------------------------------------------------------------
t = icoMP_T( 3 )
call t%buildIcosphere_()
call t%cotangentLaplacian_()
nv = t%nverts

call mem%alloc(LM, (/ nv-1, nv-1 /), 'LM', initval = 0.D0, startdims = (/ 0, 0 /))
do i = 0, size(t%lapval)-1
  LM(t%laprow(i),t%lapcol(i)) = LM(t%laprow(i),t%lapcol(i)) + t%lapval(i)
end do

dmax = 0.D0
do i = 0, nv-1
  dmax = max(dmax, dabs(sum(LM(i,:))))
end do
if (dmax.ge.1.0D-9) call Message%printError('selfTest', 'the Laplacian does not annihilate the constant vector')

dmax = 0.D0
do i = 0, nv-1
  do j = 0, nv-1
    dmax = max(dmax, dabs(LM(i,j) - LM(j,i)))
  end do
end do
if (dmax.ge.1.0D-9) call Message%printError('selfTest', 'the Laplacian is not symmetric')

asum = sum(t%vertexarea)
if (dabs(asum - 4.D0*cPi).ge.0.05D0*4.D0*cPi) then
  call Message%printError('selfTest', 'the vertex areas do not add up to about 4 pi')
end if
call mem%dealloc(LM, 'LM')

call Message%printMessage(' Laplacian: symmetric, constant null vector, area about 4 pi  OK')
call Message%printMessage(' ')
call Message%printMessage(' all self tests passed')
call Message%printMessage(' ')

end subroutine selfTest_

!--------------------------------------------------------------------------
recursive function verifySampling_(self) result(ok)
!DEC$ ATTRIBUTES DLLEXPORT :: verifySampling_
!! author: MDG
!! version: 1.0
!! date: 08/07/26
!!
!! sampling at the master pattern's own grid nodes must return the stored values,
!! and the two sides of the equator must agree
!!
!! This is the Fortran counterpart of the python --verify option.  It always runs
!! on the energy sum, and leaves the energy mode and kernel settings as it found
!! them.

use mod_io
use mod_memory
use mod_Lambert

IMPLICIT NONE

class(icoMP_T), INTENT(INOUT)   :: self
logical                         :: ok

type(IO_T)                      :: Message
type(memory_T)                  :: mem
type(Lambert_T)                 :: L
real(kind=dbl), allocatable     :: nodes(:,:), got(:,:), ref(:), seam(:,:), sgot(:,:), sref(:)
integer(kind=irg)               :: n, ix, iy, k, ierr, kk, nseam, nnodes
real(kind=dbl)                  :: ab(2), dc(3), err, jump, rng, twopi
character(fnlen)                :: charline, savemode, savekernel

mem = memory_T()
ok = .TRUE.
n = self%npx
nnodes = (2*n+1)*(2*n+1)
nseam = 20000
twopi = 2.D0*cPi

savemode = trim(self%energymode)
savekernel = trim(self%kernel)
self%energymode = 'sum'
call self%selectEnergy_()

call mem%alloc(nodes, (/ nnodes, 3 /), 'nodes')
call mem%alloc(ref, (/ nnodes /), 'ref')
call mem%alloc(got, (/ 1, nnodes /), 'got')
call mem%alloc(seam, (/ nseam, 3 /), 'seam')
call mem%alloc(sref, (/ nseam /), 'sref')
call mem%alloc(sgot, (/ 1, nseam /), 'sgot')

! the master pattern's own grid nodes, in the order that numpy would ravel them
L = Lambert_T()
k = 0
do iy = 0, 2*n
  do ix = 0, 2*n
    k = k + 1
    ab = (/ dble(ix)/dble(n) - 1.D0, dble(iy)/dble(n) - 1.D0 /)
    call L%setxyd(ab)
    ierr = L%LambertSquareToSphere(dc)
    nodes(k,:) = dc
  end do
end do

do k = 1, nseam
  ab(1) = twopi*dble(k-1)/dble(nseam)
  seam(k,:) = (/ dcos(ab(1)), dsin(ab(1)), 1.0D-9 /)
end do

do kk = 1, 2
  if (kk.eq.1) then
    self%kernel = 'bilinear'
  else
    self%kernel = 'bicubic'
  end if

! northern hemisphere nodes
  k = 0
  do iy = 0, 2*n
    do ix = 0, 2*n
      k = k + 1
      ref(k) = self%nhsel(1,ix-n,iy-n)
    end do
  end do
  call self%sampleDirections_(nnodes, nodes, got)
  rng = max(maxval(ref) - minval(ref), 1.0D-30)
  err = maxval(dabs(got(1,:) - ref)) / rng
  if (err.ge.1.0D-9) ok = .FALSE.
  write (charline,"('  ',A8,' north nodes:  max rel err  ',ES12.4)") trim(self%kernel), err
  call Message%printMessage(trim(charline))

! southern hemisphere nodes
  k = 0
  do iy = 0, 2*n
    do ix = 0, 2*n
      k = k + 1
      ref(k) = self%shsel(1,ix-n,iy-n)
    end do
  end do
  nodes(:,3) = -nodes(:,3)
  call self%sampleDirections_(nnodes, nodes, got)
  nodes(:,3) = -nodes(:,3)
  rng = max(maxval(ref) - minval(ref), 1.0D-30)
  err = maxval(dabs(got(1,:) - ref)) / rng
  if (err.ge.1.0D-9) ok = .FALSE.
  write (charline,"('  ',A8,' south nodes:  max rel err  ',ES12.4)") trim(self%kernel), err
  call Message%printMessage(trim(charline))

! the equator seam
  call self%sampleDirections_(nseam, seam, sgot)
  sref = sgot(1,:)
  seam(:,3) = -seam(:,3)
  call self%sampleDirections_(nseam, seam, sgot)
  seam(:,3) = -seam(:,3)
  rng = max(maxval(self%nhsel(1,:,:)) - minval(self%nhsel(1,:,:)), 1.0D-30)
  jump = maxval(dabs(sgot(1,:) - sref)) / rng
  if (jump.ge.1.0D-9) ok = .FALSE.
  write (charline,"('  ',A8,' equator seam: max rel jump ',ES12.4)") trim(self%kernel), jump
  call Message%printMessage(trim(charline))
end do

call mem%dealloc(nodes, 'nodes')
call mem%dealloc(ref, 'ref')
call mem%dealloc(got, 'got')
call mem%dealloc(seam, 'seam')
call mem%dealloc(sref, 'sref')
call mem%dealloc(sgot, 'sgot')

self%energymode = trim(savemode)
self%kernel = trim(savekernel)

end function verifySampling_

end module mod_icoMP
