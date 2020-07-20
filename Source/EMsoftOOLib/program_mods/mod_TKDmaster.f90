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

module mod_TKDmaster
  !! author: MDG
  !! version: 1.0
  !! date: 03/18/20
  !!
  !! class definition for the EMTKDmaster program

use mod_kinds
use mod_global
use mod_MPfiles

IMPLICIT NONE

! class definition
type, public :: TKDmaster_T
private
  character(fnlen)              :: nmldeffile = 'EMTKDmaster.nml'
  type(TKDmasterNameListType)   :: nml

contains
private
  procedure, pass(self) :: readNameList_
  ! procedure, pass(self) :: writeHDFNameList_    replaced by routine in mod_MPfiles
  procedure, pass(self) :: getNameList_
  ! procedure, pass(self) :: TKDmaster_
  procedure, pass(self) :: get_npx_
  procedure, pass(self) :: get_Esel_
  procedure, pass(self) :: get_nthreads_
  procedure, pass(self) :: get_dmin_
  procedure, pass(self) :: get_copyfromenergyfile_
  procedure, pass(self) :: get_h5copypath_
  procedure, pass(self) :: get_energyfile_
  procedure, pass(self) :: get_BetheParametersFile_
  procedure, pass(self) :: get_combinesites_
  procedure, pass(self) :: get_restart_
  procedure, pass(self) :: get_uniform_
  procedure, pass(self) :: get_Notify_
  procedure, pass(self) :: get_kinematical_
  procedure, pass(self) :: get_thickness_
  procedure, pass(self) :: set_npx_
  procedure, pass(self) :: set_Esel_
  procedure, pass(self) :: set_nthreads_
  procedure, pass(self) :: set_dmin_
  procedure, pass(self) :: set_copyfromenergyfile_
  procedure, pass(self) :: set_h5copypath_
  procedure, pass(self) :: set_energyfile_
  procedure, pass(self) :: set_BetheParametersFile_
  procedure, pass(self) :: set_combinesites_
  procedure, pass(self) :: set_restart_
  procedure, pass(self) :: set_uniform_
  procedure, pass(self) :: set_Notify_
  procedure, pass(self) :: set_kinematical_
  procedure, pass(self) :: set_thickness_

  generic, public :: getNameList => getNameList_
  ! generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  ! generic, public :: TKDmaster => TKDmaster_
  generic, public :: get_npx => get_npx_
  generic, public :: get_Esel => get_Esel_
  generic, public :: get_nthreads => get_nthreads_
  generic, public :: get_dmin => get_dmin_
  generic, public :: get_copyfromenergyfile => get_copyfromenergyfile_
  generic, public :: get_h5copypath => get_h5copypath_
  generic, public :: get_energyfile => get_energyfile_
  generic, public :: get_BetheParametersFile => get_BetheParametersFile_
  generic, public :: get_combinesites => get_combinesites_
  generic, public :: get_restart => get_restart_
  generic, public :: get_uniform => get_uniform_
  generic, public :: get_Notify => get_Notify_
  generic, public :: get_kinematical => get_kinematical_
  generic, public :: get_thickness => get_thickness_
  generic, public :: set_npx => set_npx_
  generic, public :: set_Esel => set_Esel_
  generic, public :: set_nthreads => set_nthreads_
  generic, public :: set_dmin => set_dmin_
  generic, public :: set_copyfromenergyfile => set_copyfromenergyfile_
  generic, public :: set_h5copypath => set_h5copypath_
  generic, public :: set_energyfile => set_energyfile_
  generic, public :: set_BetheParametersFile => set_BetheParametersFile_
  generic, public :: set_combinesites => set_combinesites_
  generic, public :: set_restart => set_restart_
  generic, public :: set_uniform => set_uniform_
  generic, public :: set_Notify => set_Notify_
  generic, public :: set_kinematical => set_kinematical_
  generic, public :: set_thickness => set_thickness_
end type TKDmaster_T

! the constructor routine for this class
interface TKDmaster_T
  module procedure TKDmaster_constructor
end interface TKDmaster_T

contains

!--------------------------------------------------------------------------
type(TKDmaster_T) function TKDmaster_constructor( nmlfile ) result(TKDmaster)
!DEC$ ATTRIBUTES DLLEXPORT :: TKDmaster_constructor
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! constructor for the TKDmaster_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call TKDmaster%readNameList(nmlfile)

end function TKDmaster_constructor

!--------------------------------------------------------------------------
subroutine TKDmaster_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: TKDmaster_destructor
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! destructor for the TKDmaster_T Class

IMPLICIT NONE

type(TKDmaster_T), INTENT(INOUT)  :: self

call reportDestructor('TKDmaster_T')

end subroutine TKDmaster_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! read the namelist from an nml file for the TKDmaster_T Class

use mod_io
use mod_EMsoft

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)    :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft
type(IO_T)                           :: Message
logical                              :: skipread = .FALSE.

integer(kind=irg)  :: npx
integer(kind=irg)  :: Esel
integer(kind=irg)  :: nthreads
real(kind=sgl)     :: dmin
character(fnlen)   :: copyfromenergyfile
character(fnlen)   :: h5copypath
character(fnlen)   :: energyfile
character(fnlen)   :: BetheParametersFile
logical            :: combinesites
logical            :: restart
logical            :: uniform
character(3)       :: Notify
logical            :: kinematical

! define the IO namelist to facilitate passing variables to the program.
namelist /TKDmastervars/ dmin,npx,nthreads,energyfile,Esel,restart,uniform,combinesites, &
                         copyfromenergyfile, h5copypath, BetheParametersFile, kinematical, &
                         Notify

! set the input parameters to default values (except for xtalname, which must be present)
npx = 500                       ! Nx pixels (total = 2Nx+1)
nthreads = 1
Esel = -1                       ! selected energy value for single energy run
dmin = 0.025                    ! smallest d-spacing to include in dynamical matrix [nm]
energyfile = 'undefined'        ! default filename for z_0(E_e) data from EMMC Monte Carlo simulations
combinesites = .FALSE.          ! keep asymmetric unit sites separate or not
restart = .FALSE.               ! when .TRUE. an existing file will be assumed
uniform = .FALSE.               ! when .TRUE., the output master patterns will contain 1.0 everywhere
copyfromenergyfile = 'undefined'
kinematical = .FALSE.
Notify = 'Off'
BetheParametersFile = 'BetheParameters.nml'
h5copypath = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=TKDmastervars)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(energyfile).eq.'undefined') then
  call Message%printError('EMTKDmaster:',' energy file name is undefined in '//nmlfile)
 end if
end if

! if we get here, then all appears to be ok, and we need to fill in the emnl fields
self%nml%npx = npx
self%nml%Esel = Esel
self%nml%nthreads = nthreads
self%nml%dmin = dmin
self%nml%energyfile = energyfile
self%nml%combinesites = combinesites
self%nml%restart = restart
self%nml%uniform = uniform
self%nml%BetheParametersFile = BetheParametersFile
self%nml%h5copypath = h5copypath
self%nml%Notify = Notify
self%nml%kinematical = kinematical
self%nml%copyfromenergyfile = copyfromenergyfile

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! pass the namelist for the TKDmaster_T Class to the calling program

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)          :: self
type(TKDmasterNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
function get_npx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_npx_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get npx from the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
integer(kind=irg)                     :: out

out = self%nml%npx

end function get_npx_

!--------------------------------------------------------------------------
subroutine set_npx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_npx_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set npx in the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%npx = inp

end subroutine set_npx_

!--------------------------------------------------------------------------
function get_Esel_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Esel_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get Esel from the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
integer(kind=irg)                     :: out

out = self%nml%Esel

end function get_Esel_

!--------------------------------------------------------------------------
subroutine set_Esel_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Esel_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set Esel in the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%Esel = inp

end subroutine set_Esel_

!--------------------------------------------------------------------------
function get_nthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nthreads_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get nthreads from the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
integer(kind=irg)                     :: out

out = self%nml%nthreads

end function get_nthreads_

!--------------------------------------------------------------------------
subroutine set_nthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nthreads_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set nthreads in the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%nthreads = inp

end subroutine set_nthreads_

!--------------------------------------------------------------------------
function get_dmin_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dmin_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get dmin from the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
real(kind=sgl)                        :: out

out = self%nml%dmin

end function get_dmin_

!--------------------------------------------------------------------------
subroutine set_dmin_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_dmin_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set dmin in the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)            :: inp

self%nml%dmin = inp

end subroutine set_dmin_

!--------------------------------------------------------------------------
function get_thickness_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_thickness_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get thickness from the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
real(kind=sgl)                        :: out

out = self%nml%thickness

end function get_thickness_

!--------------------------------------------------------------------------
subroutine set_thickness_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_thickness_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set thickness in the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)            :: inp

self%nml%thickness = inp

end subroutine set_thickness_

!--------------------------------------------------------------------------
function get_copyfromenergyfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_copyfromenergyfile_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get copyfromenergyfile from the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
character(fnlen)                      :: out

out = self%nml%copyfromenergyfile

end function get_copyfromenergyfile_

!--------------------------------------------------------------------------
subroutine set_copyfromenergyfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_copyfromenergyfile_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set copyfromenergyfile in the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%copyfromenergyfile = inp

end subroutine set_copyfromenergyfile_

!--------------------------------------------------------------------------
function get_h5copypath_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_h5copypath_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get h5copypath from the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
character(fnlen)                      :: out

out = self%nml%h5copypath

end function get_h5copypath_

!--------------------------------------------------------------------------
subroutine set_h5copypath_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_h5copypath_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set h5copypath in the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%h5copypath = inp

end subroutine set_h5copypath_

!--------------------------------------------------------------------------
function get_energyfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_energyfile_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get energyfile from the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
character(fnlen)                      :: out

out = self%nml%energyfile

end function get_energyfile_

!--------------------------------------------------------------------------
subroutine set_energyfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_energyfile_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set energyfile in the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%energyfile = inp

end subroutine set_energyfile_

!--------------------------------------------------------------------------
function get_BetheParametersFile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_BetheParametersFile_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get BetheParametersFile from the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
character(fnlen)                      :: out

out = self%nml%BetheParametersFile

end function get_BetheParametersFile_

!--------------------------------------------------------------------------
subroutine set_BetheParametersFile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_BetheParametersFile_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set BetheParametersFile in the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%BetheParametersFile = inp

end subroutine set_BetheParametersFile_

!--------------------------------------------------------------------------
function get_combinesites_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_combinesites_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get combinesites from the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
logical                               :: out

out = self%nml%combinesites

end function get_combinesites_

!--------------------------------------------------------------------------
subroutine set_combinesites_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_combinesites_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set combinesites in the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
logical, INTENT(IN)                   :: inp

self%nml%combinesites = inp

end subroutine set_combinesites_

!--------------------------------------------------------------------------
function get_restart_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_restart_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get restart from the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
logical                               :: out

out = self%nml%restart

end function get_restart_

!--------------------------------------------------------------------------
subroutine set_restart_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_restart_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set restart in the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
logical, INTENT(IN)                   :: inp

self%nml%restart = inp

end subroutine set_restart_

!--------------------------------------------------------------------------
function get_uniform_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_uniform_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get uniform from the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
logical                               :: out

out = self%nml%uniform

end function get_uniform_

!--------------------------------------------------------------------------
subroutine set_uniform_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_uniform_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set uniform in the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
logical, INTENT(IN)                   :: inp

self%nml%uniform = inp

end subroutine set_uniform_

!--------------------------------------------------------------------------
function get_Notify_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Notify_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get Notify from the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
character(3)                          :: out

out = self%nml%Notify

end function get_Notify_

!--------------------------------------------------------------------------
subroutine set_Notify_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Notify_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set Notify in the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
character(3), INTENT(IN)              :: inp

self%nml%Notify = inp

end subroutine set_Notify_

!--------------------------------------------------------------------------
function get_kinematical_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_kinematical_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get kinematical from the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
logical                               :: out

out = self%nml%kinematical

end function get_kinematical_

!--------------------------------------------------------------------------
subroutine set_kinematical_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_kinematical_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set kinematical in the TKDmaster_T class

IMPLICIT NONE

class(TKDmaster_T), INTENT(INOUT)     :: self
logical, INTENT(IN)                   :: inp

self%nml%kinematical = inp

end subroutine set_kinematical_

! !--------------------------------------------------------------------------
! subroutine TKDmaster_(self, EMsoft, progname)
! !! author: MDG
! !! version: 1.0
! !! date: 03/18/20
! !!
! !! perform the computations

! use mod_EMsoft
! use mod_initializers
! use mod_symmetry
! use mod_crystallography
! use mod_gvectors
! use mod_kvectors
! use mod_io
! use mod_math
! use mod_diffraction
! use mod_timing
! use mod_Lambert
! use HDF5
! use mod_HDFsupport
! use mod_HDFnames
! use ISO_C_BINDING
! use omp_lib
! use mod_notifications
! use stringconstants
! use mod_MCfiles

! IMPLICIT NONE

! class(TKDmaster_T), INTENT(INOUT)       :: self
! type(EMsoft_T), INTENT(INOUT)           :: EMsoft
! character(fnlen), INTENT(INOUT)         :: progname

! type(Cell_T)            :: cell
! type(DynType)           :: Dyn
! type(Timing_T)          :: timer
! type(IO_T)              :: Message
! type(Lambert_T)         :: L
! type(HDF_T)             :: HDF
! type(SpaceGroup_T)      :: SG
! type(Diffraction_T)     :: Diff
! type(MCfile_T)          :: MCFT
! type(MPfile_T)          :: MPFT
! type(kvectors_T)        :: kvec
! type(gvectors_T)        :: reflist
! type(HDFnames_T)        :: HDFnames

! real(kind=dbl)          :: ctmp(192,3), arg, Radius, xyz(3)
! integer(HSIZE_T)        :: dims4(4), cnt4(4), offset4(4)
! integer(HSIZE_T)        :: dims3(3), cnt3(3), offset3(3)
! integer(kind=irg)       :: isym,i,j,ik,npy,ipx,ipy,ipz,debug,iE,izz, izzmax, iequiv(3,48), nequiv, num_el, MCnthreads, & ! counters
!                            numk, timestart, timestop, numsites, nthreads, & ! number of independent incident beam directions
!                            ir,nat(100),kk(3), skip, ijmax, one, NUMTHREADS, TID, SamplingType, &
!                            numset,n,ix,iy,iz, io_int(6), nns, nnw, nref, Estart, &
!                            istat,gzero,ic,ip,ikk, totstrong, totweak, jh, ierr, nix, niy, nixp, niyp     ! counters
! real(kind=dbl)          :: tpi,Znsq, kkl, DBWF, kin, delta, h, lambda, omtl, srt, dc(3), xy(2), edge, scl, tmp, dx, dxm, dy, dym !
! real(kind=sgl)          :: io_real(5), selE, kn, FN(3), kkk(3), tstop, bp(4), nabsl, etotal, density, Ze, at_wt
! real(kind=sgl),allocatable      :: EkeVs(:), svals(:), auxNH(:,:,:,:), auxSH(:,:,:,:), Z2percent(:)  ! results
! real(kind=sgl),allocatable      :: mLPNH(:,:,:,:), mLPSH(:,:,:,:), masterSPNH(:,:,:), masterSPSH(:,:,:)
! real(kind=dbl),allocatable      :: LegendreArray(:), upd(:), diagonal(:)
! integer(kind=irg),allocatable   :: accum_z(:,:,:,:)
! complex(kind=dbl)               :: czero
! complex(kind=dbl),allocatable   :: Lgh(:,:), Sgh(:,:,:)
! logical                 :: usehex, switchmirror, verbose
! character(fnlen)        :: xtalname

! ! Monte Carlo derived quantities
! integer(kind=irg)       :: numEbins, nsx, nsy, hdferr, nlines, lastEnergy    ! variables used in MC energy file
! integer(kind=irg),allocatable :: thick(:)
! real(kind=sgl),allocatable :: lambdaE(:,:)
! character(fnlen)        :: oldprogname, groupname, energyfile, outname, datagroupname, attributename, HDF_FileVersion, fname
! character(8)            :: MCscversion
! character(11)           :: dstr
! character(15)           :: tstrb
! character(15)           :: tstre
! logical                 :: f_exists, readonly, overwrite=.TRUE., insert=.TRUE., stereog, g_exists, xtaldataread, FL, doLegendre
! character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)
! character(fnlen,kind=c_char)                     :: line2(1)

! type(gnode),save                :: rlp
! real(kind=sgl),allocatable      :: karray(:,:)
! integer(kind=irg),allocatable   :: kij(:,:)
! complex(kind=dbl),allocatable   :: DynMat(:,:)
! character(fnlen)                :: dataset, instring
! type(MCOpenCLNameListType)      :: mcnl
! type(kvectorlist), pointer      :: ktmp
! type(reflisttype), pointer      :: firstw

! character(fnlen),ALLOCATABLE    :: MessageLines(:)
! integer(kind=irg)               :: NumLines, info
! character(fnlen)                :: SlackUsername, exectime
! character(100)                  :: c


! !$OMP THREADPRIVATE(rlp)

! call openFortranHDFInterface()

! ! set the HDF group names for this program
! HDF = HDF_T()
! HDFnames = HDFnames_T()
! call MPFT%setModality('TKD')

! ! simplify the notation a little
! associate( emnl => self%nml )

! ! initialize the timing routines
! timer = Timing_T()
! tstrb = timer%getTimeString()

! ! if copyfromenergyfile is different from 'undefined', then we need to
! ! copy all the Monte Carlo data from that file into a new file, which
! ! will then be read from and written to by the ComputeMasterPattern routine.
! if (emnl%copyfromenergyfile.ne.'undefined') then
!   call MCFT%copyMCdata(EMsoft, HDF, emnl%copyfromenergyfile, emnl%energyfile, emnl%h5copypath)
! end if

! stereog = .TRUE.

! tpi = 2.D0*cPi
! czero = cmplx(0.D0,0.D0)

! ! is the master pattern used for spherical indexing only ?  If so, then we need to modifiy the k-vector sampling
! doLegendre = .FALSE.

! !=============================================
! !=============================================
! ! ---------- read Monte Carlo .h5 output file and extract necessary parameters
! ! set the HDF group names for reading the MC input file
! call HDFnames%set_ProgramData(SC_MCOpenCL)
! call HDFnames%set_NMLlist(SC_MCCLNameList)
! call HDFnames%set_NMLfilename(SC_MCOpenCLNML)
! fname = EMsoft%generateFilePath('EMdatapathname',trim(emnl%energyfile))
! call MCFT%setFileName(fname)
! call MCFT%readMCfile(HDF, HDFnames, getAccumz=.TRUE.)
! mcnl = MCFT%getnml()
! call MCFT%copyaccumz(accum_z)

! nsx = (mcnl%numsx - 1)/2
! nsy = nsx
! etotal = float(mcnl%totnum_el)

! io_int(1) = mcnl%totnum_el
! call Message%WriteValue(' --> total number of BSE electrons in MC data set ', io_int, 1)



! !=============================================
! ! should we create a new file or open an existing file?
! !=============================================
!   lastEnergy = -1
!   energyfile = trim(EMsoft%generateFilePath('EMdatapathname',emnl%energyfile))
!   outname = trim(energyfile)

! ! set the HDFnames to the correct strings for this program
! call HDFnames%set_ProgramData(SC_TKDmaster)
! call HDFnames%set_NMLlist(SC_TKDmasterNameList)
! call HDFnames%set_NMLfilename(SC_TKDmasterNML)
! call HDFnames%set_Variable(SC_MCOpenCL)

! if (emnl%restart.eqv..TRUE.) then
! ! in this case we need to check whether or not the file exists, then open
! ! it and read the value of the last energy level that was simulated and written
! ! to that file; if this level is different from the lowest energy level we
! ! know that there is at least one more level to be simulated.  If it is equal,
! ! then we can abort the program here.

!   inquire(file=trim(outname), exist=f_exists)
!   if (.not.f_exists) then
!     call Message%printError('TKDmaster','restart HDF5 file does not exist')
!   end if

! !=============================================
! ! open the existing HDF5 file
! !=============================================
!   datagroupname = trim(HDFnames%get_ProgramData())

! ! Create a new file using the default properties.
!   readonly = .TRUE.
!   hdferr =  HDF%openFile(outname, readonly)

! ! all we need to get from the file is the lastEnergy parameter
!   hdferr = HDF%openGroup(HDFnames%get_EMData())
!   hdferr = HDF%openGroup(datagroupname)

! dataset = SC_lastEnergy
!   call HDF%readDatasetInteger(dataset, hdferr, lastEnergy)

!   call HDF%pop(.TRUE.)
! end if
! !=============================================
! !=============================================
! ! crystallography section;
! verbose = .TRUE.

! call cell%setFileName(mcnl%xtalname)
! call Diff%setrlpmethod('WK')

! if (emnl%restart.eqv..TRUE.) then
!   call Diff%setV(dble(EkeVs(lastEnergy-1)))
!   call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, emnl%dmin, verbose, useHDF=HDF)
! else
!   call Diff%setV(dble(mcnl%EkeV))
!   call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, emnl%dmin, verbose, useHDF=HDF)
! end if

! ! check the crystal system and setting; abort the program for trigonal with rhombohedral setting with
! ! an explanation for the user

! if ((SG%getSpaceGroupXtalSystem().eq.5).and.(cell%getLatParm('b').eq.cell%getLatParm('c'))) then
!     call Message%printMessage( (/ &
!     '                                                                         ', &
!     ' ========Program Aborted========                                         ', &
!     ' The TKD master pattern simulation for rhombohedral/trigonal structures  ', &
!     ' requires that the structure be described using the hexagonal reference  ', &
!     ' frame.  Please re-enter the crystal structure in this setting and re-run', &
!     ' the Monte Carlo calculation and this master pattern program.            '/) )
!     stop
! end if

! ! then calculate density, average atomic number and average atomic weight
! call cell%calcDensity(Z2percent)
! io_real(1:3) = cell%getDensity()
! density = io_real(1)
! at_wt = io_real(2)
! Ze = io_real(3)
! call Message%WriteValue('Density, avA, avZ = ',io_real,3,"(2f10.5,',',f10.5)")

! ! allocate and compute the Sgh loop-up table
!  numset = cell%getNatomtype()
!  call Diff%Initialize_SghLUT(cell, SG, emnl%dmin, numset, nat, verbose)

! ! determine the point group number
!  j=0
!  do i=1,32
!   if (SGPG(i).le.SG%getSpaceGroupNumber()) j=i
!  end do
!  isym = j

! ! here is new code dealing with all the special cases (quite a few more compared to the
! ! Laue group case)...  isym is the point group number. Once the symmetry case has been
! ! fully determined (taking into account things like 31m and 3m1 an such), then the only places
! ! that symmetry is handled are the modified Calckvectors routine, and the filling of the modified
! ! Lambert projections after the dynamical simulation step.  We are also changing the name of the
! ! sr array (or srhex) to mLPNH and mLPSH (modified Lambert Projection Northern/Southern Hemisphere),
! ! and we change the output HDF5 file a little as well. We need to make sure that the EMEBSD program
! ! issues a warning when an old format HDF5 file is read.

! ! Here, we encode isym into a new number that describes the sampling scheme; the new schemes are
! ! described in detail in the EBSD manual pdf file.

! SamplingType = PGSamplingType(isym)

! ! next, intercept the special cases (hexagonal vs. rhombohedral cases that require special treatment)
! if ((SamplingType.eq.-1).or.(isym.eq.14).or.(isym.eq.26)) then
!   SamplingType = SG%getHexvsRho(isym)
! end if

! ! if the point group is trigonal or hexagonal, we need to switch usehex to .TRUE. so that
! ! the program will use the hexagonal sampling method
! usehex = .FALSE.
! if ((SG%getSpaceGroupXtalSystem().eq.4).or.(SG%getSpaceGroupXtalSystem().eq.5)) usehex = .TRUE.

! ! ---------- end of symmetry and crystallography section
! !=============================================
! !=============================================
! numEbins = MCFT%getnumEbins()

! allocate(EkeVs(numEbins),thick(numEbins))

! do i=1,numEbins
!   EkeVs(i) = mcnl%Ehistmin + float(i-1)*mcnl%Ebinsize
! end do

! ! then, for each energy determine the 95% histogram thickness
! izzmax = 0
! do iE = 1,numEbins
!  do ix=-nsx/10,nsx/10
!   do iy=-nsy/10,nsy/10
!    istat = sum(accum_z(iE,:,ix,iy))
!    izz = 1
!    do while (sum(accum_z(iE,1:izz,ix,iy)).lt.(0.99*istat))
!     izz = izz+1
!    end do
!    if (izz.gt.izzmax) izzmax = izz
!   end do
!  end do
!  thick(iE) = dble(izzmax) * mcnl%depthstep
! end do

! izz = nint(maxval(thick)/mcnl%depthstep)
! allocate(lambdaE(1:numEbins,1:izz),stat=istat)
! do iE=1,numEbins
!  call Diff%setV(dble(Ekevs(iE)))
!  call Diff%CalcUcg(cell,(/0,0,0/))
!  rlp = Diff%getrlp()
!  nabsl = rlp%xgp
!  do iz=1,izz
!   lambdaE(iE,iz) = float(sum(accum_z(iE,iz,-nsx/10:nsx/10,-nsy/10:nsy/10)))/etotal
!   lambdaE(iE,iz) = lambdaE(iE,iz) * exp(2.0*sngl(cPi)*(iz-1)*mcnl%depthstep/nabsl)
!  end do
! end do

! !=============================================
! !=============================================
! ! ---------- a couple of initializations
!    npy = emnl%npx
!    allocate(svals(numset),stat=istat)
!    gzero = 1  ! index of incident beam
! ! ----------
! !=============================================
! !=============================================





! end subroutine TKDmaster_

end module mod_TKDmaster
