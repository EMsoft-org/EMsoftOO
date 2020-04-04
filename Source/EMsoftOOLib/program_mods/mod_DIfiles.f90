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

module mod_DIfiles
  !! author: MDG 
  !! version: 1.0 
  !! date: 04/03/20
  !!
  !! class definition for writing a Dictionary Indexing output file; for now
  !! we only write the TSL format (h5ebsd), but we apply it both to EBSD and TKD.
  !! ECP will be added at a later time, and when Oxford finally decides to release
  !! their HDF file we will likely add that format as well. 
  !! In the original EMsoft package, this was the EMh5ebsd module.
  !!
  !! Also contains everything dealing with dictionary indexing namelists.

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
  integer(kind=irg)  :: binning
  integer(kind=irg)  :: nthreads
  integer(kind=irg)  :: devid
  integer(kind=irg)  :: usenumd
  integer(kind=irg)  :: multidevid(8)
  integer(kind=irg)  :: platid
  integer(kind=irg)  :: nregions
  integer(kind=irg)  :: nlines
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
  character(1)       :: maskpattern
  character(3)       :: scalingmode
  character(3)       :: Notify
  character(1)       :: keeptmpfile
  character(fnlen)   :: anglefile
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
end type DictionaryIndexingNameListType

type, public, extends(DictionaryIndexingNameListType) :: EBSDDINameListType
end type EBSDDINameListType

type, public, extends(DictionaryIndexingNameListType) :: ECPDINameListType
end type ECPDINameListType

type, public, extends(DictionaryIndexingNameListType) :: TKDDINameListType
end type TKDDINameListType


type, public :: MPdataType

end type MPdataType

! class definition
type, public :: DIfile_T
private 
  character(fnlen)                              :: DIfile
  type(MPdatatype)                              :: DIDT
  character(fnlen)                              :: modality = 'unknown'
  type(DictionaryIndexingNameListType), public  :: nml
contains
private 

  procedure, pass(self) :: get_filename_
  procedure, pass(self) :: set_filename_
  procedure, pass(self) :: get_Modality_
  procedure, pass(self) :: set_Modality_
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: writeHDFNameList_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: get_filename => get_filename_
  generic, public :: set_filename => set_filename_
  generic, public :: get_Modality => get_Modality_
  generic, public :: set_Modality => set_Modality_
  generic, public :: writeHDFNameList => writeHDFNameList_

end type DIfile_T

!DEC$ ATTRIBUTES DLLEXPORT :: getNameList
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList
!DEC$ ATTRIBUTES DLLEXPORT :: get_filename
!DEC$ ATTRIBUTES DLLEXPORT :: set_filename
!DEC$ ATTRIBUTES DLLEXPORT :: get_Modality
!DEC$ ATTRIBUTES DLLEXPORT :: set_Modality
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList

! the constructor routine for this class 
interface DIfile_T
  module procedure DIfile_constructor
end interface DIfile_T

contains

!--------------------------------------------------------------------------
type(DIfile_T) function DIfile_constructor( nmlfile, fname ) result(DIfile)
!! author: MDG 
!! version: 1.0 
!! date: 04/03/20
!!
!! constructor for the DIfile_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 
character(fnlen), OPTIONAL   :: fname

if (present(nmlfile)) call DIfile%readNameList(nmlfile)
if (present(fname)) call DIfile%set_filename(fname)

end function DIfile_constructor

!--------------------------------------------------------------------------
subroutine DIfile_destructor(self) 
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
function get_filename_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 04/03/20
!!
!! get filename from the DIfile_T class

IMPLICIT NONE 

class(DIfile_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = self%DIfile

end function get_filename_

!--------------------------------------------------------------------------
subroutine set_filename_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 04/03/20
!!
!! set filename in the DIfile_T class

IMPLICIT NONE 

class(DIfile_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%DIfile = inp

end subroutine set_filename_

!--------------------------------------------------------------------------
function get_Modality_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/31/20
!!
!! get Modality from the DPfile_T class

IMPLICIT NONE 

class(DIfile_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = self%Modality

end function get_Modality_

!--------------------------------------------------------------------------
subroutine set_Modality_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/31/20
!!
!! set Modality in the DPfile_T class

IMPLICIT NONE 

class(DIfile_T), INTENT(INOUT)     :: self
character(*), INTENT(IN)           :: inp

self%Modality = inp

end subroutine set_Modality_


!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!! author: MDG 
!! version: 1.0 
!! date: 03/31/20
!!
!! read the namelist from an nml file for the DI_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(DIfile_T), INTENT(INOUT)      :: self
character(fnlen),INTENT(IN)         :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)         :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                      :: EMsoft 
type(IO_T)                          :: Message       
logical                             :: skipread = .FALSE.

integer(kind=irg)  :: numsx
integer(kind=irg)  :: numsy
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
character(1)       :: maskpattern
character(1)       :: keeptmpfile
character(3)       :: scalingmode
character(3)       :: Notify
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

! define the IO namelist to facilitate passing variables to the program.
namelist  / EBSDIndexingdata / thetac, delta, numsx, numsy, xpc, ypc, masterfile, devid, platid, inputtype, &
                               beamcurrent, dwelltime, binning, gammavalue, energymin, nregions, nlines, maskfile, &
                               scalingmode, maskpattern, L, omega, nthreads, energymax, datafile, angfile, ctffile, &
                               ncubochoric, numexptsingle, numdictsingle, ipf_ht, ipf_wd, nnk, nnav, exptfile, maskradius, &
                               dictfile, indexingmode, hipassw, stepX, stepY, tmpfile, avctffile, nosm, eulerfile, Notify, &
                               HDFstrings, ROI, keeptmpfile, multidevid, usenumd, nism, isangle, refinementNMLfile

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
ROI             = (/ 0, 0, 0, 0 /)  ! Region of interest (/ x0, y0, w, h /)
maskradius      = 240
binning         = 1             ! binning mode  (1, 2, 4, or 8)
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
keeptmpfile     = 'n'
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
inputtype       = 'Binary'    ! Binary, EMEBSD, TSLHDF, TSLup2, OxfordHDF, OxfordBinary, BrukerHDF 
HDFstrings      = ''

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=EBSDIndexingdata)
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

    if (numsx.eq.0) then 
        call Message%printError('readNameList:',' pattern size numsx is zero in '//nmlfile)
    end if

    if (numsy.eq.0) then 
        call Message%printError('readNameList:',' pattern size numsy is zero in '//nmlfile)
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
self%nml%exptfile      = exptfile
self%nml%nnk           = nnk
self%nml%nnav          = nnav
self%nml%nosm          = nosm
self%nml%nism          = nism
self%nml%isangle       = isangle
self%nml%ipf_ht        = ipf_ht
self%nml%ipf_wd        = ipf_wd
self%nml%nthreads      = nthreads
self%nml%datafile      = datafile
self%nml%tmpfile       = tmpfile
self%nml%ctffile       = ctffile
self%nml%avctffile     = avctffile
self%nml%angfile       = angfile
self%nml%eulerfile     = eulerfile
self%nml%maskradius    = maskradius
self%nml%numdictsingle = numdictsingle
self%nml%numexptsingle = numexptsingle
self%nml%hipassw       = hipassw
self%nml%masterfile    = masterfile
self%nml%energyfile    = masterfile
self%nml%maskfile      = maskfile
self%nml%StepX         = stepX
self%nml%StepY         = stepY
self%nml%indexingmode  = trim(indexingmode)
self%nml%Notify        = Notify
self%nml%inputtype     = inputtype
self%nml%HDFstrings    = HDFstrings
self%nml%L             = L
self%nml%numsx         = numsx
self%nml%numsy         = numsy
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
self%nml%dictfile      = dictfile 
self%nml%refinementNMLfile = refinementNMLfile

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
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
integer(kind=irg)                                   :: n_int, n_real
integer(kind=irg)                                   :: hdferr
integer(kind=irg),allocatable                       :: io_int(:)
real(kind=sgl),allocatable                          :: io_real(:)
character(20),allocatable                           :: intlist(:), reallist(:)
character(fnlen)                                    :: dataset, sval(1),groupname
character(fnlen,kind=c_char)                        :: line2(1), line10(10)
logical                                             :: g_exists, overwrite=.TRUE., isEBSD=.FALSE., &
                                                       isECP=.FALSE., isTKD=.FALSE., isEBSDSHT=.FALSE.

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

select type(emnl) 
  type is (EBSDDINameListType)
    isEBSD = .TRUE.
    n_int = 19
    n_real = 15
    allocate( io_int(n_int), intlist(n_int), io_real(n_real), reallist(n_real) )
  type is (ECPDINameListType)
    isECP = .TRUE.
    n_int = 19
    n_real = 15
    allocate( io_int(n_int), intlist(n_int), io_real(n_real), reallist(n_real) )
  class default 
    call Message%printError('writeHDFNameList', 'unknown name list type requested')
end select

! write all the single integers
io_int = (/ emnl%ncubochoric, emnl%numexptsingle, emnl%numdictsingle, emnl%ipf_ht, &
            emnl%ipf_wd, emnl%nnk, emnl%maskradius, emnl%numsx, emnl%numsy, emnl%binning, &
            emnl%nthreads, emnl%devid, emnl%platid, emnl%nregions, emnl%nnav, &
            emnl%nosm, emnl%nlines, emnl%usenumd, emnl%nism /)
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
call HDF%writeNMLintegers(io_int, intlist, n_int)

io_real = (/ emnl%L, emnl%thetac, emnl%delta, emnl%omega, emnl%xpc, &
             emnl%ypc, emnl%energymin, emnl%energymax, emnl%gammavalue, emnl%StepX, &
             emnl%stepY, emnl%isangle, emnl%beamcurrent, emnl%dwelltime, emnl%hipassw /)
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

dataset = SC_anglefile
line2(1) = emnl%anglefile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create anglefile dataset', hdferr)

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

! and pop this group off the stack
call HDF%pop()

end subroutine writeHDFNameList_


end module mod_DIfiles