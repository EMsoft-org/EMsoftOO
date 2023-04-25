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

module mod_Kosselmaster
  !! author: MDG
  !! version: 1.0
  !! date: 03/25/20
  !!
  !! class definition for the EMKosselmaster program

use mod_kinds
use mod_global

IMPLICIT NONE

private :: CalcKthick, CalcKint

! namelist for the EMKosselmaster program
type, public :: KosselmasterNameListType
  integer(kind=irg)       :: numthick
  integer(kind=irg)       :: npx
  integer(kind=irg)       :: nthreads
  real(kind=sgl)          :: voltage
  real(kind=sgl)          :: dmin
  real(kind=sgl)          :: startthick
  real(kind=sgl)          :: thickinc
  real(kind=sgl)          :: tfraction
  character(6)            :: Kosselmode
  character(fnlen)        :: xtalname
  character(fnlen)        :: outname
  character(fnlen)        :: BetheParametersFile
end type KosselmasterNameListType

! class definition
type, public :: Kosselmaster_T
private
  character(fnlen)                :: nmldeffile = 'EMKosselmaster.nml'
  type(KosselmasterNameListType)  :: nml

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: Kosselmaster_
  procedure, pass(self) :: get_numthick_
  procedure, pass(self) :: get_npx_
  procedure, pass(self) :: get_nthreads_
  procedure, pass(self) :: get_voltage_
  procedure, pass(self) :: get_dmin_
  procedure, pass(self) :: get_startthick_
  procedure, pass(self) :: get_thickinc_
  procedure, pass(self) :: get_tfraction_
  procedure, pass(self) :: get_Kosselmode_
  procedure, pass(self) :: get_xtalname_
  procedure, pass(self) :: get_outname_
  procedure, pass(self) :: get_BetheParametersFile_
  procedure, pass(self) :: set_numthick_
  procedure, pass(self) :: set_npx_
  procedure, pass(self) :: set_nthreads_
  procedure, pass(self) :: set_voltage_
  procedure, pass(self) :: set_dmin_
  procedure, pass(self) :: set_startthick_
  procedure, pass(self) :: set_thickinc_
  procedure, pass(self) :: set_tfraction_
  procedure, pass(self) :: set_Kosselmode_
  procedure, pass(self) :: set_xtalname_
  procedure, pass(self) :: set_outname_
  procedure, pass(self) :: set_BetheParametersFile_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: Kosselmaster => Kosselmaster_
  generic, public :: get_numthick => get_numthick_
  generic, public :: get_npx => get_npx_
  generic, public :: get_nthreads => get_nthreads_
  generic, public :: get_voltage => get_voltage_
  generic, public :: get_dmin => get_dmin_
  generic, public :: get_startthick => get_startthick_
  generic, public :: get_thickinc => get_thickinc_
  generic, public :: get_tfraction => get_tfraction_
  generic, public :: get_Kosselmode => get_Kosselmode_
  generic, public :: get_xtalname => get_xtalname_
  generic, public :: get_outname => get_outname_
  generic, public :: get_BetheParametersFile => get_BetheParametersFile_
  generic, public :: set_numthick => set_numthick_
  generic, public :: set_npx => set_npx_
  generic, public :: set_nthreads => set_nthreads_
  generic, public :: set_voltage => set_voltage_
  generic, public :: set_dmin => set_dmin_
  generic, public :: set_startthick => set_startthick_
  generic, public :: set_thickinc => set_thickinc_
  generic, public :: set_tfraction => set_tfraction_
  generic, public :: set_Kosselmode => set_Kosselmode_
  generic, public :: set_xtalname => set_xtalname_
  generic, public :: set_outname => set_outname_
  generic, public :: set_BetheParametersFile => set_BetheParametersFile_
end type Kosselmaster_T

! the constructor routine for this class
interface Kosselmaster_T
  module procedure Kosselmaster_constructor
end interface Kosselmaster_T

contains

!--------------------------------------------------------------------------
type(Kosselmaster_T) function Kosselmaster_constructor( nmlfile ) result(Kosselmaster)
!DEC$ ATTRIBUTES DLLEXPORT :: Kosselmaster_constructor
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! constructor for the Kosselmaster_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call Kosselmaster%readNameList(nmlfile)

end function Kosselmaster_constructor

!--------------------------------------------------------------------------
subroutine Kosselmaster_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT ::  Kosselmaster_destructor
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! destructor for the Kosselmaster_T Class

IMPLICIT NONE

type(Kosselmaster_T), INTENT(INOUT)  :: self

call reportDestructor('Kosselmaster_T')

end subroutine Kosselmaster_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! read the namelist from an nml file for the Kosselmaster_T Class

use mod_io
use mod_EMsoft

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT) :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft
type(IO_T)                           :: Message
logical                              :: skipread = .FALSE.

integer(kind=irg)       :: numthick
integer(kind=irg)       :: npx
integer(kind=irg)       :: nthreads
real(kind=sgl)          :: voltage
real(kind=sgl)          :: dmin
real(kind=sgl)          :: startthick
real(kind=sgl)          :: thickinc
real(kind=sgl)          :: tfraction
character(6)            :: Kosselmode
character(fnlen)        :: xtalname
character(fnlen)        :: outname
character(fnlen)        :: BetheParametersFile

namelist /Kosselmasterlist/ xtalname, voltage, dmin, nthreads, BetheParametersFile, &
                            startthick, thickinc, numthick, tfraction, outname, npx, Kosselmode

! set the input parameters to default values (except for xtalname, which must be present)
numthick = 10                   ! number of increments
npx = 256                       ! output arrays will have size npix x npix
nthreads = 4                    ! default number of threads for OpenMP
voltage = 20.0                  ! acceleration voltage [kV]
dmin = 0.05                     ! smallest d-spacing to include in dynamical matrix [nm]
startthick = 10.0               ! starting thickness [nm]
thickinc = 10.0                 ! thickness increment
xtalname = 'undefined'          ! initial value to check that the keyword is present in the nml file
outname = 'Kosselout.data'      ! output filename
BetheParametersFile = 'BetheParameters.nml' ! file name for Bethe parameters
Kosselmode = 'normal'           ! 'thicks' for thickness determination, 'normal' for normal plot
tfraction = 0.1                 ! thickness fraction for 'thicks' mode

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=Kosselmasterlist)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(xtalname).eq.'undefined') then
  call Message%printError('readNameList:',' structure file name is undefined in '//nmlfile)
 end if
end if

! if we get here, then all appears to be ok, and we need to fill in the knl fields
self%nml%numthick = numthick
self%nml%npx = npx
self%nml%nthreads = nthreads
self%nml%voltage = voltage
self%nml%dmin = dmin
self%nml%startthick = startthick
self%nml%thickinc = thickinc
self%nml%tfraction = tfraction
self%nml%Kosselmode = Kosselmode
self%nml%xtalname = xtalname
self%nml%outname = outname
self%nml%BetheParametersFile = BetheParametersFile

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! pass the namelist for the Kosselmaster_T Class to the calling program

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)          :: self
type(KosselmasterNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants

use ISO_C_BINDING

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)    :: self
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 3, n_real = 5
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( knl => self%nml )

! create the group for this namelist
groupname = trim(HDFnames%get_NMLlist())
hdferr = HDF%createGroup(groupname)

! write all the single integers
io_int = (/ knl%numthick, knl%npx, knl%nthreads /)
intlist(1) = 'numthick'
intlist(2) = 'npx'
intlist(3) = 'nthreads'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ knl%voltage, knl%dmin, knl%startthick, knl%thickinc, knl%tfraction /)
reallist(1) = 'voltage'
reallist(2) = 'dmin'
reallist(3) = 'startthick'
reallist(4) = 'thickinc'
reallist(5) = 'tfraction'
call HDF%writeNMLreals(io_real, reallist, n_real)

! write all the strings
dataset = SC_Kosselmode
line2(1) = knl%Kosselmode
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create Kosselmode dataset',hdferr)

dataset = SC_xtalname
line2(1) = knl%xtalname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create xtalname dataset',hdferr)

dataset = SC_outname
line2(1) = knl%outname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create outname dataset',hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_


!--------------------------------------------------------------------------
function get_numthick_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_numthick_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! get numthick from the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
integer(kind=irg)                        :: out

out = self%nml%numthick

end function get_numthick_

!--------------------------------------------------------------------------
subroutine set_numthick_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_numthick_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! set numthick in the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)            :: inp

self%nml%numthick = inp

end subroutine set_numthick_

!--------------------------------------------------------------------------
function get_npx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_npx_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! get npx from the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
integer(kind=irg)                        :: out

out = self%nml%npx

end function get_npx_

!--------------------------------------------------------------------------
subroutine set_npx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_npx_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! set npx in the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)            :: inp

self%nml%npx = inp

end subroutine set_npx_

!--------------------------------------------------------------------------
function get_nthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nthreads_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! get nthreads from the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
integer(kind=irg)                        :: out

out = self%nml%nthreads

end function get_nthreads_

!--------------------------------------------------------------------------
subroutine set_nthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nthreads_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! set nthreads in the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)            :: inp

self%nml%nthreads = inp

end subroutine set_nthreads_

!--------------------------------------------------------------------------
function get_voltage_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_voltage_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! get voltage from the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
real(kind=sgl)                           :: out

out = self%nml%voltage

end function get_voltage_

!--------------------------------------------------------------------------
subroutine set_voltage_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_voltage_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! set voltage in the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)               :: inp

self%nml%voltage = inp

end subroutine set_voltage_

!--------------------------------------------------------------------------
function get_dmin_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dmin_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! get dmin from the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
real(kind=sgl)                           :: out

out = self%nml%dmin

end function get_dmin_

!--------------------------------------------------------------------------
subroutine set_dmin_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_dmin_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! set dmin in the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)               :: inp

self%nml%dmin = inp

end subroutine set_dmin_

!--------------------------------------------------------------------------
function get_startthick_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_startthick_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! get startthick from the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
real(kind=sgl)                           :: out

out = self%nml%startthick

end function get_startthick_

!--------------------------------------------------------------------------
subroutine set_startthick_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_startthick_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! set startthick in the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)               :: inp

self%nml%startthick = inp

end subroutine set_startthick_

!--------------------------------------------------------------------------
function get_thickinc_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_thickinc_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! get thickinc from the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
real(kind=sgl)                           :: out

out = self%nml%thickinc

end function get_thickinc_

!--------------------------------------------------------------------------
subroutine set_thickinc_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_thickinc_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! set thickinc in the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)               :: inp

self%nml%thickinc = inp

end subroutine set_thickinc_

!--------------------------------------------------------------------------
function get_tfraction_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_tfraction_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! get tfraction from the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
real(kind=sgl)                           :: out

out = self%nml%tfraction

end function get_tfraction_

!--------------------------------------------------------------------------
subroutine set_tfraction_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_tfraction_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! set tfraction in the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)               :: inp

self%nml%tfraction = inp

end subroutine set_tfraction_

!--------------------------------------------------------------------------
function get_Kosselmode_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Kosselmode_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! get Kosselmode from the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
character(6)                             :: out

out = self%nml%Kosselmode

end function get_Kosselmode_

!--------------------------------------------------------------------------
subroutine set_Kosselmode_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Kosselmode_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! set Kosselmode in the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
character(6), INTENT(IN)                 :: inp

self%nml%Kosselmode = inp

end subroutine set_Kosselmode_

!--------------------------------------------------------------------------
function get_xtalname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_xtalname_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! get xtalname from the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
character(fnlen)                         :: out

out = self%nml%xtalname

end function get_xtalname_

!--------------------------------------------------------------------------
subroutine set_xtalname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_xtalname_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! set xtalname in the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)             :: inp

self%nml%xtalname = inp

end subroutine set_xtalname_

!--------------------------------------------------------------------------
function get_outname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_outname_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! get outname from the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
character(fnlen)                         :: out

out = self%nml%outname

end function get_outname_

!--------------------------------------------------------------------------
subroutine set_outname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_outname_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! set outname in the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)             :: inp

self%nml%outname = inp

end subroutine set_outname_

!--------------------------------------------------------------------------
function get_BetheParametersFile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_BetheParametersFile_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! get BetheParametersFile from the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
character(fnlen)                         :: out

out = self%nml%BetheParametersFile

end function get_BetheParametersFile_

!--------------------------------------------------------------------------
subroutine set_BetheParametersFile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_BetheParametersFile_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! set BetheParametersFile in the Kosselmaster_T class

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)             :: inp

self%nml%BetheParametersFile = inp

end subroutine set_BetheParametersFile_

!--------------------------------------------------------------------------
subroutine Kosselmaster_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: Kosselmaster_
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! perform the computations

use mod_EMsoft
use mod_initializers
use mod_symmetry
use mod_crystallography
use mod_gvectors
use mod_kvectors
use mod_io
use mod_diffraction
use mod_timing
use mod_Lambert
use mod_MPfiles
use mod_math
use HDF5
use mod_HDFsupport
use mod_HDFnames
use ISO_C_BINDING
use omp_lib
use mod_OMPsupport
use stringconstants

IMPLICIT NONE

class(Kosselmaster_T), INTENT(INOUT)    :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname

type(Cell_T)                    :: cell
type(DynType)                   :: Dyn
type(Timing_T)                  :: timer
type(IO_T)                      :: Message
type(Lambert_T)                 :: L
type(HDF_T)                     :: HDF
type(SpaceGroup_T)              :: SG
type(Diffraction_T)             :: Diff
type(MPfile_T)                  :: MPFT
type(kvectors_T)                :: kvec
type(gvectors_T)                :: reflist
type(HDFnames_T)                :: HDFnames

real(kind=dbl)                  :: ctmp(192,3), arg
integer(HSIZE_T)                :: dims3(3), cnt3(3), offset3(3)
integer(HSIZE_T)                :: dims2(2), cnt2(2), offset2(2)
integer(kind=irg)               :: isym,i,j,ik,npy,ipx,ipy,ipz,debug,izz, izzmax, iequiv(3,48), nequiv, num_el, MCnthreads, & ! counters
                                   SamplingType,numk,numthick,  nthreads, & ! number of independent incident beam directions
                                   ir,kk(3), npyhex, skip, ijmax, one, NUMTHREADS, TID, hdferr, tickstart, &
                                   n,ix,iy, io_int(6), nns, nnw, nref, nix, niy, nixp, niyp, ierr, &
                                   istat,gzero,ic,ip,ikk, totstrong, totweak     ! counters
real(kind=dbl)                  :: tpi,Znsq, kkl, DBWF, kin, xy(2), dc(3), edge, scl, tmp, dx, dxm, dy, dym, xyz(3), Radius !!
real(kind=sgl)                  :: io_real(5), selE, kn, FN(3), kkk(3), bp(4), tstop, xyzs(3), tav
complex(kind=dbl)               :: czero
real(kind=sgl),allocatable      :: mLPNH(:,:,:), mLPSH(:,:,:), Iz(:), thick(:), trange(:,:), masterSPNH(:,:,:), &
                                   masterSPSH(:,:,:), trangeSP(:,:)
real(kind=sgl),allocatable      :: auxNH(:,:,:), auxSH(:,:,:), auxtrange(:,:)
logical                         :: usehex, switchmirror, verbose, insert=.TRUE., overwrite=.TRUE., silent=.TRUE., g_exists
character(11)                   :: dstr
character(15)                   :: tstrb
character(15)                   :: tstre
character(fnlen)                :: oldprogname, groupname, energyfile, outname, attributename, datagroupname, HDF_FileVersion
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)
character(fnlen,kind=c_char)                     :: line2(1)

type(reflisttype), pointer      :: firstw
type(kvectorlist), pointer      :: ktmp
type(gnode),save                :: rlp
type(BetheParameterType)        :: BetheParameters
real(kind=sgl),allocatable      :: karray(:,:)
integer(kind=irg),allocatable   :: kij(:,:)
complex(kind=dbl),allocatable   :: DynMat(:,:)
character(fnlen)                :: dataset, instring

!$OMP THREADPRIVATE(rlp)

call openFortranHDFInterface()

associate( kmnl=>self%nml )

! set the HDF group names for this program
HDF = HDF_T()
HDFnames = HDFnames_T()

! initialize the timing routines
timer = Timing_T()
tstrb = timer%getTimeString()

tpi = 2.D0*cPi
czero = cmplx(0.D0,0.D0)
npy = kmnl%npx
gzero = 1  ! index of incident beam

!=============================================
!=============================================
! crystallography section
verbose = .TRUE.

call cell%setFileName(kmnl%xtalname)
call Diff%setrlpmethod('WK')

call Diff%setV(dble(kmnl%voltage))
call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, kmnl%dmin, verbose, useHDF=HDF)

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

! determine the point group number
 j=0
 do i=1,32
  if (SGPG(i).le.SG%getSpaceGroupNumber()) j=i
 end do
 isym = j

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
! this is where we determine the values for the thicknesses
if (kmnl%Kosselmode.eq.'thicks') then
  numthick = 1
else
  numthick = kmnl%numthick
end if

allocate(thick(numthick),Iz(numthick))

do i=1,numthick
  thick(i) = kmnl%startthick + float(i-1)*kmnl%thickinc
end do
!=============================================
!=============================================

!=============================================
!=============================================
! ---------- allocate memory for the master pattern
! we need to sample the stereographic projection Northern hemisphere or a portion
! thereoff, depending on the current symmetry.
if (kmnl%Kosselmode.eq.'normal') then
  allocate(mLPNH(-kmnl%npx:kmnl%npx,-npy:npy,1:kmnl%numthick),stat=istat)
  allocate(mLPSH(-kmnl%npx:kmnl%npx,-npy:npy,1:kmnl%numthick),stat=istat)
  allocate(masterSPNH(-kmnl%npx:kmnl%npx,-npy:npy,1:kmnl%numthick))
  allocate(masterSPSH(-kmnl%npx:kmnl%npx,-npy:npy,1:kmnl%numthick))

! set various arrays to zero
   mLPNH = 0.0
   mLPSH = 0.0
else
! Kosselmode must be 'thicks'
  allocate(trange(-kmnl%npx:kmnl%npx,-npy:npy),trangeSP(-kmnl%npx:kmnl%npx,-npy:npy),stat=istat)
  trange = 0.0
  trangeSP = 0.0
end if
! ---------- end allocate memory for the master pattern
!=============================================
!=============================================

! force dynamical matrix routine to read new Bethe parameters from file
  call Diff%SetBetheParameters(EMsoft, .FALSE., kmnl%BetheParametersFile)

!=============================================
! create the HDF5 output file
!=============================================
  HDF = HDF_T()
  call HDFnames%set_ProgramData(SC_Kosselmaster)
  call HDFnames%set_NMLlist(SC_KosselmasterNameList)
  call HDFnames%set_NMLfilename(SC_KosselmasterNML)

  outname = EMsoft%generateFilePath('EMdatapathname',trim(kmnl%outname))

! Create a new file using the default properties.
  hdferr =  HDF%createFile(outname)

! write the EMheader to the file
groupname = trim(HDFnames%get_ProgramData())
  call HDF%writeEMheader(EMsoft,dstr, tstrb, tstre, progname, groupname)

! add the CrystalData group at the top level of the file
call cell%addXtalDataGroup(SG, EMsoft, HDF)

! create a namelist group to write all the namelist files into
  hdferr = HDF%createGroup(HDFnames%get_NMLfiles())

! read the text file and write the array to the file
dataset = trim(HDFnames%get_NMLfilename())
  hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%nmldeffile)

! leave this group
  call HDF%pop()

! create a namelist group to write all the namelist files into
  hdferr = HDF%createGroup(HDFnames%get_NMLparameters())
  call self%writeHDFNameList(HDF, HDFnames)
  call Diff%writeBetheparameterNameList(HDF)

! leave this group
  call HDF%pop()

! then the remainder of the data in a EMData group
  hdferr = HDF%createGroup(HDFnames%get_EMData())
  hdferr = HDF%createGroup(HDFnames%get_ProgramData())

! add a HDF_FileVersion attribute
  HDF_FileVersion = '4.0'
  HDF_FileVersion = cstringify(HDF_FileVersion)
  attributename = SC_HDFFileVersion
  hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)

dataset = SC_xtalname
  allocate(stringarray(1))
  stringarray(1)= trim(kmnl%xtalname)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1, overwrite)
  else
    hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)
  end if

dataset = SC_npx
  hdferr = HDF%writeDatasetInteger(dataset, kmnl%npx)

dataset = SC_BetheParameters
  bp(1) = Diff%getBetheParameter('c1')
  bp(2) = Diff%getBetheParameter('c2')
  bp(3) = Diff%getBetheParameter('c3')
  bp(4) = Diff%getBetheParameter('sg')
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetFloatArray(dataset, bp, 4, overwrite)
  else
    hdferr = HDF%writeDatasetFloatArray(dataset, bp, 4)
  end if

dataset = SC_numthick
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetInteger(dataset, kmnl%numthick, overwrite)
  else
    hdferr = HDF%writeDatasetInteger(dataset, kmnl%numthick)
  end if

dataset = SC_startthick
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetFloat(dataset, kmnl%startthick, overwrite)
  else
    hdferr = HDF%writeDatasetFloat(dataset, kmnl%startthick)
  end if

dataset = SC_thickinc
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetFloat(dataset, kmnl%thickinc, overwrite)
  else
    hdferr = HDF%writeDatasetFloat(dataset, kmnl%thickinc)
  end if

dataset = SC_Kosselmode
  stringarray(1)= trim(kmnl%Kosselmode)
  hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

! create the hyperslabs and write zeroes to them for now
  if (kmnl%Kosselmode.eq.'normal') then
dataset = SC_mLPNH
    dims3 = (/  2*kmnl%npx+1, 2*kmnl%npx+1, numthick /)
    cnt3 = (/ 2*kmnl%npx+1, 2*kmnl%npx+1, numthick /)
    offset3 = (/ 0, 0, 0 /)
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims3, offset3, cnt3)

dataset = SC_mLPSH
    dims3 = (/  2*kmnl%npx+1, 2*kmnl%npx+1, numthick /)
    cnt3 = (/ 2*kmnl%npx+1, 2*kmnl%npx+1, numthick /)
    offset3 = (/ 0, 0, 0 /)
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims3, offset3, cnt3)
  else
dataset = SC_trange
    dims2 = (/  2*kmnl%npx+1, 2*kmnl%npx+1 /)
    cnt2 = (/ 2*kmnl%npx+1, 2*kmnl%npx+1 /)
    offset2 = (/ 0, 0 /)
    hdferr = HDF%writeHyperslabFloatArray(dataset, trange, dims2, offset2, cnt2)

dataset = SC_trangeSP
    dims2 = (/  2*kmnl%npx+1, 2*kmnl%npx+1 /)
    cnt2 = (/ 2*kmnl%npx+1, 2*kmnl%npx+1 /)
    offset2 = (/ 0, 0 /)
    hdferr = HDF%writeHyperslabFloatArray(dataset, trangeSP, dims2, offset2, cnt2)
  end if

  call HDF%popall()
!=============================================
! completes the HDF5 output file
!=============================================

!=============================================
!=============================================
! print a message to indicate we're starting the computation
call Message%printMessage('Starting computation', frm = "(/A)")

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
  call kvec%Calckvectors(cell, SG, Diff, (/ 0.D0, 0.D0, 0.D0 /),kmnl%npx,npy, ijmax,usehex)
else
  call kvec%Calckvectors(cell, SG, Diff, (/ 0.D0, 0.D0, 0.D0 /),kmnl%npx,npy, ijmax,usehex)
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
call timer%Time_tick(1)

! set the number of OpenMP threads
call OMP_setNThreads(kmnl%nthreads)

! use OpenMP to run on multiple cores ...
!$OMP PARALLEL COPYIN(rlp) &
!$OMP& PRIVATE(DynMat,ik,FN,TID,kn,ipx,ipy,ix,iequiv,nequiv,reflist,firstw) &
!$OMP& PRIVATE(kkk,nns,nnw,nref,io_int,Iz)

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

  call reflist%Initialize_ReflectionList(cell, SG, Diff, FN, kkk, self%nml%dmin, verbose)
  nref = reflist%get_nref()
! ---------- end of "create the master reflection list"
!=============================================

! determine strong and weak reflections
  nullify(firstw)
  nns = 0
  nnw = 0
  call reflist%Apply_BethePotentials(Diff, firstw, nns, nnw)

! generate the dynamical matrix
  if (allocated(DynMat)) deallocate(DynMat)
  allocate(DynMat(nns,nns))
  call reflist%GetDynMat(cell, Diff, firstw, DynMat, nns, nnw)
  totstrong = totstrong + nns
  totweak = totweak + nnw

! solve the dynamical eigenvalue equation for this beam direction
  kn = karray(4,ik)
  if (self%nml%Kosselmode.eq.'thicks') then
    call CalcKthick(DynMat,kn,nns,self%nml%tfraction,Iz)
  else
    call CalcKint(DynMat,kn,nns,numthick,thick,Iz)
  end if
  deallocate(DynMat)

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
  if (self%nml%Kosselmode.eq.'normal') then
    do ix=1,nequiv
      if (iequiv(3,ix).eq.-1) mLPSH(iequiv(1,ix),iequiv(2,ix),1:numthick) = Iz(1:numthick)
      if (iequiv(3,ix).eq.1) mLPNH(iequiv(1,ix),iequiv(2,ix),1:numthick) = Iz(1:numthick)
    end do
  else
    do ix=1,nequiv
      trange(iequiv(1,ix),iequiv(2,ix)) = Iz(1)
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

deallocate(karray, kij)

if (usehex) then
! and finally, we convert the hexagonally sampled array to a square Lambert projection which will be used
! for all EBSD pattern interpolations;  we need to do this for both the Northern and Southern hemispheres

! we begin by allocating auxiliary arrays to hold copies of the hexagonal data; the original arrays will
! then be overwritten with the newly interpolated data.
  if (self%nml%Kosselmode.eq.'normal') then
    allocate(auxNH(-kmnl%npx:kmnl%npx,-npy:npy,1:kmnl%numthick),stat=istat)
    allocate(auxSH(-kmnl%npx:kmnl%npx,-npy:npy,1:kmnl%numthick),stat=istat)
    auxNH = mLPNH
    auxSH = mLPSH
  else
    allocate(auxtrange(-kmnl%npx:kmnl%npx,-npy:npy),stat=istat)
    auxtrange = trange
  end if

  edge = 1.D0 / dble(kmnl%npx)
  scl = float(kmnl%npx)
  do i=-kmnl%npx,kmnl%npx
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
        if (nixp.gt.kmnl%npx) nixp = nix
        if (niyp.gt.kmnl%npx) niyp = niy
        dx = xy(1) - nix
        dy = xy(2) - niy
        dxm = 1.D0 - dx
        dym = 1.D0 - dy
        if (kmnl%Kosselmode.eq.'normal') then
          mLPNH(i,j,1:numthick) = auxNH(nix,niy,1:numthick)*dxm*dym + auxNH(nixp,niy,1:numthick)*dx*dym + &
                                  auxNH(nix,niyp,1:numthick)*dxm*dy + auxNH(nixp,niyp,1:numthick)*dx*dy
          mLPSH(i,j,1:numthick) = auxSH(nix,niy,1:numthick)*dxm*dym + auxSH(nixp,niy,1:numthick)*dx*dym + &
                                  auxSH(nix,niyp,1:numthick)*dxm*dy + auxSH(nixp,niyp,1:numthick)*dx*dy
        else
          trange(i,j) = auxtrange(nix,niy)*dxm*dym + auxtrange(nixp,niy)*dx*dym + &
                        auxtrange(nix,niyp)*dxm*dy + auxtrange(nixp,niyp)*dx*dy
        end if
      end if
    end do
  end do
  deallocate(auxNH, auxSH)
end if

! make sure that the outer pixel rim of the mLPSH patterns is identical to
! that of the mLPNH array.
if (kmnl%Kosselmode.eq.'normal') then
  mLPSH(-kmnl%npx,-kmnl%npx:kmnl%npx,1:numthick) = mLPNH(-kmnl%npx,-kmnl%npx:kmnl%npx,1:numthick)
  mLPSH( kmnl%npx,-kmnl%npx:kmnl%npx,1:numthick) = mLPNH( kmnl%npx,-kmnl%npx:kmnl%npx,1:numthick)
  mLPSH(-kmnl%npx:kmnl%npx,-kmnl%npx,1:numthick) = mLPNH(-kmnl%npx:kmnl%npx,-kmnl%npx,1:numthick)
  mLPSH(-kmnl%npx:kmnl%npx, kmnl%npx,1:numthick) = mLPNH(-kmnl%npx:kmnl%npx, kmnl%npx,1:numthick)

! get stereographic projections (summed over the atomic positions)
  Radius = 1.0
  do i=-kmnl%npx,kmnl%npx
    do j=-kmnl%npx,kmnl%npx
      L = Lambert_T( xyd = (/ dble(i), dble(j) /) / dble(kmnl%npx) )
      ierr = L%StereoGraphicInverse( xyz, Radius )
      xyz = xyz/vecnorm(xyz)
      if (ierr.ne.0) then
        masterSPNH(i,j,1) = 0.0
        masterSPSH(i,j,1) = 0.0
      else
        masterSPNH(i,j,1:numthick) = InterpolateLambert(xyz, mLPNH, kmnl%npx, numthick)
        masterSPSH(i,j,1:numthick) = InterpolateLambert(xyz, mLPSH, kmnl%npx, numthick)
      end if
    end do
  end do
else
  ! set the background to the average intensity level
  tav = sum(trange)/float(2*kmnl%npx+1)**2
  ! get stereographic projections (summed over the atomic positions)
  Radius = 1.0
  do i=-kmnl%npx,kmnl%npx
    do j=-kmnl%npx,kmnl%npx
      L = Lambert_T( xyd = (/ dble(i), dble(j) /) / dble(kmnl%npx) )
      ierr = L%StereoGraphicInverse( xyz, Radius )
      xyz = xyz/vecnorm(xyz)
      if (ierr.ne.0) then
        trangeSP(i,j) = tav
      else
        xyzs = sngl(xyz)
        trangeSP(i,j) = InterpolateLambert(xyzs, trange, kmnl%npx)
      end if
    end do
  end do
end if

io_int(1) = nint(float(totstrong)/float(numk))
call Message%WriteValue(' -> Average number of strong reflections = ',io_int, 1, "(I5)")
io_int(1) = nint(float(totweak)/float(numk))
call Message%WriteValue(' -> Average number of weak reflections   = ',io_int, 1, "(I5)")

call timer%makeTimeStamp()
dstr = timer%getDateString()
tstre = timer%getTimeString()

!===================================
!===================================
! and write the final output arrays to the HDF5 file
datagroupname = trim(HDFnames%get_ProgramData())

! open the existing file using the default properties.
hdferr =  HDF%openFile(outname)

! update the time string
hdferr = HDF%openGroup(HDFnames%get_EMheader())
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
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetFloat(dataset, tstop, overwrite)
  else
    hdferr = HDF%writeDatasetFloat(dataset, tstop)
  end if

  call HDF%pop()
  call HDF%pop()

  hdferr = HDF%openGroup(HDFnames%get_EMData())
  hdferr = HDF%openGroup(datagroupname)

! add data to the hyperslab
if (kmnl%Kosselmode.eq.'normal') then
  dataset = SC_mLPNH
    dims3 = (/  2*kmnl%npx+1, 2*kmnl%npx+1, numthick /)
    cnt3 = (/ 2*kmnl%npx+1, 2*kmnl%npx+1, numthick /)
    offset3 = (/ 0, 0, 0 /)
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims3, offset3, cnt3, insert)

  dataset = SC_mLPSH
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims3, offset3, cnt3, insert)

  dataset = SC_masterSPNH
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then
      hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPNH, dims3, offset3, cnt3, insert)
    else
      hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPNH, dims3, offset3, cnt3)
    end if

  dataset = SC_masterSPSH
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then
      hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPSH, dims3, offset3, cnt3, insert)
    else
      hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPSH, dims3, offset3, cnt3)
    end if

else
  dataset = SC_trange
    dims2 = (/  2*kmnl%npx+1, 2*kmnl%npx+1 /)
    cnt2 = (/ 2*kmnl%npx+1, 2*kmnl%npx+1 /)
    offset2 = (/ 0, 0 /)
    hdferr = HDF%writeHyperslabFloatArray(dataset, trange, dims2, offset2, cnt2, insert)

  dataset = SC_trangeSP
    dims2 = (/  2*kmnl%npx+1, 2*kmnl%npx+1 /)
    cnt2 = (/ 2*kmnl%npx+1, 2*kmnl%npx+1 /)
    offset2 = (/ 0, 0 /)
    hdferr = HDF%writeHyperslabFloatArray(dataset, trangeSP, dims2, offset2, cnt2, insert)
end if

call HDF%popall()

! and close the fortran hdf interface
call closeFortranHDFInterface()

call Message%printMessage('Final data stored in file '//trim(kmnl%outname), frm = "(A/)")

end associate

end subroutine Kosselmaster_

!--------------------------------------------------------------------------
recursive subroutine CalcKthick(DynMat,kn,nn,thresh,Iz)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcKthick
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! compute the thickness for which the Kossel intensity drops below a treshold

use mod_io
use mod_diffraction
use mod_kvectors
use mod_gvectors

IMPLICIT NONE

integer(kind=irg),INTENT(IN)    :: nn                   !< number of strong beams
complex(kind=dbl),INTENT(IN)    :: DynMat(nn,nn)
real(kind=sgl),INTENT(IN)       :: kn
real(kind=sgl),INTENT(IN)       :: thresh               !< thickness fraction parameter
real(kind=sgl),INTENT(INOUT)    :: Iz(1)                !< output (thickness)
!f2py intent(in,out) ::  Iz

type(diffraction_T)             :: Diff
integer(kind=irg)               :: j, IPIV(nn), k
complex(kind=dbl)               :: CGinv(nn,nn), Minp(nn,nn), Wloc(nn), lCG(nn,nn), lW(nn), lalpha(nn)
real(kind=dbl)                  :: s(nn), q(nn), t, ss

! compute the eigenvalues and eigenvectors
 Minp = DynMat
 IPIV = 0
 call Diff%BWsolve(Minp,Wloc,lCG,CGinv,nn,IPIV)

! the alpha coefficients are in the first column of the inverse matrix
 lW = cPi*Wloc/cmplx(kn,0.0)
 lalpha(1:nn) = CGinv(1:nn,1)

! make sure the alpha excitation coefficients are normalized
 ss = sum(abs(lalpha(1:nn))**2)
 if (ss.ne.1.D0) then
  ss = cmplx(1.D0/dsqrt(ss),0.D0)
  lalpha = lalpha*ss
 endif

! compute the thickness value in steps of 0.25 nm until less than thresh
 do j=1,nn
  q(j) = -4.D0*cPi*aimag(lW(j))
  s(j) = abs(lalpha(j))**2
 end do
 t = 0.D0
 do
   t = t+0.25D0
   ss = sum(s*dexp(t*q))
   if (ss.le.thresh) EXIT
 end do
 Iz(1) = t

end subroutine CalcKthick

!--------------------------------------------------------------------------
recursive subroutine CalcKint(DynMat,kn,nn,nt,thick,Iz)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcKint
!! author: MDG
!! version: 1.0
!! date: 03/25/20
!!
!! compute the Kossel intensities for a range of thicknesses

use mod_io
use mod_diffraction
use mod_kvectors
use mod_gvectors

IMPLICIT NONE

integer(kind=irg),INTENT(IN)    :: nn                   !< number of strong beams
complex(kind=dbl),INTENT(IN)    :: DynMat(nn,nn)
real(kind=sgl),INTENT(IN)       :: kn
integer(kind=irg),INTENT(IN)    :: nt                   !< number of thickness values
real(kind=sgl),INTENT(IN)       :: thick(nt)            !< thickness array
real(kind=sgl),INTENT(INOUT)    :: Iz(nt)               !< output intensities
!f2py intent(in,out) ::  Iz

type(diffraction_T)             :: Diff
integer(kind=irg)               :: j, IPIV(nn), k
complex(kind=dbl)               :: CGinv(nn,nn), Minp(nn,nn), Wloc(nn), lCG(nn,nn), lW(nn), lalpha(nn)
real(kind=dbl)                  :: s, q, t

! compute the eigenvalues and eigenvectors
 Minp = DynMat
 IPIV = 0
 call Diff%BWsolve(Minp,Wloc,lCG,CGinv,nn,IPIV)

! the alpha coefficients are in the first column of the inverse matrix
 lW = cPi*Wloc/cmplx(kn,0.0)
 lalpha(1:nn) = CGinv(1:nn,1)

! make sure the alpha excitation coefficients are normalized
 s = sum(abs(lalpha(1:nn))**2)
 if (s.ne.1.D0) then
  s = cmplx(1.D0/dsqrt(s),0.D0)
  lalpha = lalpha*s
 endif

! compute the thickness array
 Iz = 0.D0
 do j=1,nn
    q = -4.D0*cPi*aimag(lW(j))
    s = abs(lalpha(j))**2
    do k=1,nt
      t = q*thick(k)
      if (abs(t).lt.30.D0) Iz(k) = Iz(k) +  s * exp(t)
    end do
 end do

end subroutine CalcKint


end module mod_Kosselmaster
