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

module mod_EBSDmaster
  !! author: MDG
  !! version: 1.0
  !! date: 02/05/20
  !!
  !! class definition for the EMEBSDmaster program

use mod_kinds
use mod_global
use mod_MPfiles

IMPLICIT NONE
private


! class definition
type, public :: EBSDmaster_T
  private
    character(fnlen)              :: nmldeffile = 'EMEBSDmaster.nml'
    type(EBSDmasterNameListType)  :: nml

  contains
  private

    procedure, pass(self) :: readNameList_
    procedure, pass(self) :: getNameList_
    procedure, pass(self) :: EBSDmaster_
    procedure, pass(self) :: testEBSDmasterWrapper_
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
    ! procedure, pass(self) :: get_thickness_
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
    ! procedure, pass(self) :: set_thickness_


    generic, public :: getNameList => getNameList_
    generic, public :: readNameList => readNameList_
    generic, public :: EBSDmaster => EBSDmaster_
    generic, public :: testEBSDmasterWrapper => testEBSDmasterWrapper_
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
    ! generic, public :: get_thickness => get_thickness_
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
    ! generic, public :: set_thickness => set_thickness_
 end type EBSDmaster_T

! the constructor routine for this class
interface EBSDmaster_T
  module procedure EBSDmaster_constructor
end interface EBSDmaster_T

contains

!--------------------------------------------------------------------------
type(EBSDmaster_T) function EBSDmaster_constructor( nmlfile ) result(EBSDmaster)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDmaster_constructor
!! author: MDG
!! version: 1.0
!! date: 02/05/20
!!
!! constructor for the EBSDmaster_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call EBSDmaster%readNameList(nmlfile)

end function EBSDmaster_constructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 02/05/20
!!
!! read the namelist from an nml file for the EBSDmaster_T Class

use mod_io

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)                 :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)                 :: initonly
 !! fill in the default values only; do not read the file

type(IO_T)                                  :: Message
logical                                     :: skipread = .FALSE.

integer(kind=irg) :: npx
integer(kind=irg) :: Esel
integer(kind=irg) :: nthreads
real(kind=sgl)    :: dmin
character(3)      :: Notify
character(fnlen)  :: copyfromenergyfile
character(fnlen)  :: energyfile
character(fnlen)  :: BetheParametersFile
character(fnlen)  :: h5copypath
logical           :: combinesites
logical           :: restart
logical           :: uniform
logical           :: kinematical

! define the IO namelist to facilitate passing variables to the program.
namelist /EBSDmastervars/ dmin,npx,nthreads,copyfromenergyfile,energyfile,Esel,restart,uniform,Notify, &
                          combinesites,h5copypath,BetheParametersFile,kinematical

! set the input parameters to default values (except for xtalname, which must be present)
npx = 500                       ! Nx pixels (total = 2Nx+1)
nthreads = 1
Esel = -1                       ! selected energy value for single energy run
dmin = 0.025                    ! smallest d-spacing to include in dynamical matrix [nm]
Notify = 'Off'
copyfromenergyfile = 'undefined'! default filename for z_0(E_e) data from a different Monte Carlo simulation
h5copypath = 'undefined'
energyfile = 'undefined'        ! default filename for z_0(E_e) data from EMMC Monte Carlo simulations
BetheParametersFile='BetheParameters.nml'
combinesites = .FALSE.          ! combine all atom sites into one BSE yield or not
restart = .FALSE.               ! when .TRUE. an existing file will be assumed
uniform = .FALSE.               ! when .TRUE., the output master patterns will contain 1.0 everywhere
kinematical = .FALSE.           ! use the kinematical approximation if .TRUE.

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=EBSDmastervars)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(energyfile).eq.'undefined') then
  call Message%printError('readNameList:',' output (energy) file name is undefined in '//nmlfile)
 end if
end if

! if we get here, then all appears to be ok, and we need to fill in the nml fields
self%nml%npx = npx
self%nml%Esel = Esel
self%nml%nthreads = nthreads
self%nml%dmin = dmin
self%nml%copyfromenergyfile = copyfromenergyfile
self%nml%h5copypath = h5copypath
self%nml%energyfile = energyfile
self%nml%BetheParametersFile = BetheParametersFile
self%nml%Notify = Notify
self%nml%combinesites = combinesites
self%nml%restart = restart
self%nml%uniform = uniform
self%nml%kinematical = kinematical

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 02/05/20
!!
!! pass the namelist for the EBSDmaster_T Class to the calling program

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)          :: self
type(EBSDmasterNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
function get_npx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_npx_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get npx from the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
integer(kind=irg)                      :: out

out = self%nml%npx

end function get_npx_

!--------------------------------------------------------------------------
subroutine set_npx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_npx_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set npx in the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)          :: inp

self%nml%npx = inp

end subroutine set_npx_

!--------------------------------------------------------------------------
function get_Esel_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Esel_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get Esel from the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
integer(kind=irg)                      :: out

out = self%nml%Esel

end function get_Esel_

!--------------------------------------------------------------------------
subroutine set_Esel_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Esel_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set Esel in the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)          :: inp

self%nml%Esel = inp

end subroutine set_Esel_

!--------------------------------------------------------------------------
function get_nthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT ::get_nthreads_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get nthreads from the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
integer(kind=irg)                      :: out

out = self%nml%nthreads

end function get_nthreads_

!--------------------------------------------------------------------------
subroutine set_nthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT ::set_nthreads_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set nthreads in the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)          :: inp

self%nml%nthreads = inp

end subroutine set_nthreads_

!--------------------------------------------------------------------------
function get_dmin_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dmin_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get dmin from the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
real(kind=sgl)                         :: out

out = self%nml%dmin

end function get_dmin_

!--------------------------------------------------------------------------
subroutine set_dmin_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_dmin_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set dmin in the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)             :: inp

self%nml%dmin = inp

end subroutine set_dmin_

! !--------------------------------------------------------------------------
! function get_thickness_(self) result(out)
! !! author: MDG
! !! version: 1.0
! !! date: 03/18/20
! !!
! !! get thickness from the EBSDmaster_T class

! IMPLICIT NONE

! class(EBSDmaster_T), INTENT(INOUT)     :: self
! real(kind=sgl)                         :: out

! out = self%nml%thickness

! end function get_thickness_

! !--------------------------------------------------------------------------
! subroutine set_thickness_(self,inp)
! !! author: MDG
! !! version: 1.0
! !! date: 03/18/20
! !!
! !! set thickness in the EBSDmaster_T class

! IMPLICIT NONE

! class(EBSDmaster_T), INTENT(INOUT)     :: self
! real(kind=sgl), INTENT(IN)             :: inp

! self%nml%thickness = inp

! end subroutine set_thickness_

!--------------------------------------------------------------------------
function get_copyfromenergyfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_copyfromenergyfile_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get copyfromenergyfile from the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
character(fnlen)                       :: out

out = self%nml%copyfromenergyfile

end function get_copyfromenergyfile_

!--------------------------------------------------------------------------
subroutine set_copyfromenergyfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_copyfromenergyfile_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set copyfromenergyfile in the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)           :: inp

self%nml%copyfromenergyfile = inp

end subroutine set_copyfromenergyfile_

!--------------------------------------------------------------------------
function get_h5copypath_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_h5copypath_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get h5copypath from the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
character(fnlen)                       :: out

out = self%nml%h5copypath

end function get_h5copypath_

!--------------------------------------------------------------------------
subroutine set_h5copypath_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_h5copypath_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set h5copypath in the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)           :: inp

self%nml%h5copypath = inp

end subroutine set_h5copypath_

!--------------------------------------------------------------------------
function get_energyfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_energyfile_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get energyfile from the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
character(fnlen)                       :: out

out = self%nml%energyfile

end function get_energyfile_

!--------------------------------------------------------------------------
subroutine set_energyfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_energyfile_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set energyfile in the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)           :: inp

self%nml%energyfile = inp

end subroutine set_energyfile_

!--------------------------------------------------------------------------
function get_BetheParametersFile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_BetheParametersFile_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get BetheParametersFile from the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
character(fnlen)                       :: out

out = self%nml%BetheParametersFile

end function get_BetheParametersFile_

!--------------------------------------------------------------------------
subroutine set_BetheParametersFile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_BetheParametersFile_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set BetheParametersFile in the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)           :: inp

self%nml%BetheParametersFile = inp

end subroutine set_BetheParametersFile_

!--------------------------------------------------------------------------
function get_combinesites_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_combinesites_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get combinesites from the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
logical                                :: out

out = self%nml%combinesites

end function get_combinesites_

!--------------------------------------------------------------------------
subroutine set_combinesites_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_combinesites_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set combinesites in the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
logical, INTENT(IN)                    :: inp

self%nml%combinesites = inp

end subroutine set_combinesites_

!--------------------------------------------------------------------------
function get_restart_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_restart_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get restart from the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
logical                                :: out

out = self%nml%restart

end function get_restart_

!--------------------------------------------------------------------------
subroutine set_restart_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_restart_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set restart in the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
logical, INTENT(IN)                    :: inp

self%nml%restart = inp

end subroutine set_restart_

!--------------------------------------------------------------------------
function get_uniform_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_uniform_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get uniform from the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
logical                                :: out

out = self%nml%uniform

end function get_uniform_

!--------------------------------------------------------------------------
subroutine set_uniform_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_uniform_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set uniform in the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
logical, INTENT(IN)                    :: inp

self%nml%uniform = inp

end subroutine set_uniform_

!--------------------------------------------------------------------------
function get_Notify_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Notify_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get Notify from the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
character(3)                           :: out

out = self%nml%Notify

end function get_Notify_

!--------------------------------------------------------------------------
subroutine set_Notify_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Notify_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set Notify in the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
character(3), INTENT(IN)               :: inp

self%nml%Notify = inp

end subroutine set_Notify_

!--------------------------------------------------------------------------
function get_kinematical_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_kinematical_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! get kinematical from the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
logical                                :: out

out = self%nml%kinematical

end function get_kinematical_

!--------------------------------------------------------------------------
subroutine set_kinematical_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_kinematical_
!! author: MDG
!! version: 1.0
!! date: 03/18/20
!!
!! set kinematical in the EBSDmaster_T class

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)     :: self
logical, INTENT(IN)                    :: inp

self%nml%kinematical = inp

end subroutine set_kinematical_

!--------------------------------------------------------------------------
subroutine EBSDmaster_(self, EMsoft, progname, HDFnames, thickness)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDmaster_
!! author: MDG
!! version: 1.0
!! date: 02/05/20
!!
!! compute an EBSD or TKD master pattern as a function of energy
!!
!! The TKD and EBSD master pattern code is pretty much identical; it is just
!! the Monte Carlo input data that is different, and we use the HDFnames class
!! to differentiate the HDF5 output files (different group names)

use mod_EMsoft
use mod_initializers
use mod_symmetry
use mod_crystallography
use mod_gvectors
use mod_kvectors
use mod_io
use mod_math
use mod_diffraction
use mod_timing
use mod_Lambert
use HDF5
use mod_HDFsupport
use mod_HDFnames
use ISO_C_BINDING
use omp_lib
use mod_OMPsupport
use mod_memory
use mod_notifications
use stringconstants
use mod_MCfiles

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)  :: self
type(EMsoft_T), INTENT(INOUT)       :: EMsoft
character(fnlen),INTENT(IN)         :: progname
type(HDFnames_T),INTENT(INOUT)      :: HDFnames
real(kind=sgl),INTENT(IN),OPTIONAL  :: thickness  ! if present, then TKD mode

type(Cell_T)            :: cell
type(DynType)           :: Dyn
type(Timing_T)          :: timer
type(IO_T)              :: Message
type(Lambert_T)         :: L
type(HDF_T)             :: HDF
type(SpaceGroup_T)      :: SG
type(Diffraction_T)     :: Diff
type(MCfile_T)          :: MCFT
type(MPfile_T)          :: MPFT
type(kvectors_T)        :: kvec
type(gvectors_T)        :: reflist
type(HDFnames_T)        :: saveHDFnames
type(memory_T)          :: mem, memth

real(kind=dbl)          :: ctmp(192,3), arg, Radius, xyz(3)
integer(HSIZE_T)        :: dims4(4), cnt4(4), offset4(4)
integer(HSIZE_T)        :: dims3(3), cnt3(3), offset3(3)
integer(kind=irg)       :: isym,i,j,ik,npy,ipx,ipy,ipz,debug,iE,izz, izzmax, iequiv(3,48), nequiv, num_el, MCnthreads, & ! counters
                           numk, timestart, timestop, numsites, nthreads, & ! number of independent incident beam directions
                           ir,nat(maxpasym),kk(3), skip, ijmax, one, NUMTHREADS, TID, SamplingType, &
                           numset,n,ix,iy,iz, io_int(6), nns, nnw, nref, Estart, &
                           istat,gzero,ic,ip,ikk, totstrong, totweak, jh, ierr, nix, niy, nixp, niyp     ! counters
real(kind=dbl)          :: tpi,Znsq, kkl, DBWF, kin, delta, h, lambda, omtl, srt, dc(3), xy(2), edge, scl, tmp, dx, dxm, dy, dym !
real(kind=sgl)          :: io_real(5), selE, kn, FN(3), kkk(3), tstop, nabsl, etotal, density, Ze, at_wt, bp(4)
real(kind=sgl),allocatable      :: EkeVs(:), svals(:), auxNH(:,:,:,:), auxSH(:,:,:,:), Z2percent(:)  ! results
real(kind=sgl),allocatable      :: mLPNH(:,:,:,:), mLPSH(:,:,:,:), masterSPNH(:,:,:), masterSPSH(:,:,:)
real(kind=dbl),allocatable      :: LegendreArray(:), upd(:), diagonal(:)
integer(kind=irg),allocatable   :: accum_z(:,:,:,:)
complex(kind=dbl)               :: czero
complex(kind=dbl),allocatable   :: Lgh(:,:), Sgh(:,:,:)
logical                 :: usehex, switchmirror, verbose
character(fnlen)        :: xtalname

! Monte Carlo derived quantities
integer(kind=irg)       :: numEbins, nsx, nsy, hdferr, nlines, lastEnergy    ! variables used in MC energy file
integer(kind=irg),allocatable :: thick(:)
real(kind=sgl),allocatable :: lambdaE(:,:)
character(fnlen)        :: oldprogname, groupname, energyfile, outname, datagroupname, attributename, HDF_FileVersion, fname
character(8)            :: MCscversion
character(11)           :: dstr
character(15)           :: tstrb
character(15)           :: tstre
logical                 :: f_exists, readonly, overwrite=.TRUE., insert=.TRUE., stereog, g_exists, xtaldataread, FL, &
                           doLegendre, isTKD = .FALSE.
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)
character(fnlen,kind=c_char)                     :: line2(1)

type(gnode),save                :: rlp
real(kind=sgl),allocatable      :: karray(:,:)
integer(kind=irg),allocatable   :: kij(:,:)
complex(kind=dbl),allocatable   :: DynMat(:,:)
character(fnlen)                :: dataset, instring
type(MCOpenCLNameListType)      :: mcnl
type(kvectorlist), pointer      :: ktmp
type(reflisttype), pointer      :: firstw

character(fnlen),ALLOCATABLE    :: MessageLines(:)
integer(kind=irg)               :: NumLines, info
character(fnlen)                :: SlackUsername, exectime
character(100)                  :: c

!$OMP THREADPRIVATE(rlp)

if (present(thickness)) then
  isTKD = .TRUE.
  saveHDFnames = HDFnames
end if

call openFortranHDFInterface()

! set the HDF group names for this program
HDF = HDF_T()
HDFnames = HDFnames_T()
if (isTKD.eqv..TRUE.) then
  call MPFT%setModality('TKD')
else
  call MPFT%setModality('EBSD')
end if

! simplify the notation a little
associate( emnl => self%nml )

! initialize the timing routines
timer = Timing_T()
tstrb = timer%getTimeString()

! initialize the memory class 
mem = memory_T()

! if copyfromenergyfile is different from 'undefined', then we need to
! copy all the Monte Carlo data from that file into a new file, which
! will then be read from and written to by the ComputeMasterPattern routine.
if (emnl%copyfromenergyfile.ne.'undefined') then
  call MCFT%copyMCdata(EMsoft, HDF, emnl%copyfromenergyfile, emnl%energyfile, emnl%h5copypath)
end if

stereog = .TRUE.

tpi = 2.D0*cPi
czero = cmplx(0.D0,0.D0)

! is the master pattern used for spherical indexing only ?  If so, then we need to modifiy the k-vector sampling
doLegendre = .FALSE.

!=============================================
!=============================================
! ---------- read Monte Carlo .h5 output file and extract necessary parameters
! set the HDF group names for reading the MC input file
call HDFnames%set_ProgramData(SC_MCOpenCL)
call HDFnames%set_NMLlist(SC_MCCLNameList)
call HDFnames%set_NMLfilename(SC_MCOpenCLNML)
fname = EMsoft%generateFilePath('EMdatapathname',trim(emnl%energyfile))
call MCFT%setFileName(fname)
call MCFT%readMCfile(HDF, HDFnames, getAccumz=.TRUE.)
mcnl = MCFT%getnml()
call MCFT%copyaccumz(accum_z)

! set the HDFnames to the correct strings for this program
if (isTKD.eqv..TRUE.) then
  HDFnames = saveHDFnames
else
  call HDFnames%set_ProgramData(SC_EBSDmaster)
  call HDFnames%set_NMLlist(SC_EBSDmasterNameList)
  call HDFnames%set_NMLfilename(SC_EBSDmasterNML)
  call HDFnames%set_Variable(SC_MCOpenCL)
end if

nsx = (mcnl%numsx - 1)/2
nsy = nsx
etotal = float(mcnl%totnum_el)

io_int(1) = mcnl%totnum_el
call Message%WriteValue(' --> total number of BSE electrons in MC data set ', io_int, 1)
!=============================================
!=============================================
numEbins = MCFT%getnumEbins()

call mem%alloc( EkeVs, (/ numEbins /), 'EkeVs')
call mem%alloc( thick, (/ numEbins /), 'thick')

do i=1,numEbins
  EkeVs(i) = mcnl%Ehistmin + float(i-1)*mcnl%Ebinsize
end do

!=============================================
! should we create a new file or open an existing file?
!=============================================
  lastEnergy = -1
  energyfile = trim(EMsoft%generateFilePath('EMdatapathname',emnl%energyfile))
  outname = trim(energyfile)

if (emnl%restart.eqv..TRUE.) then
! in this case we need to check whether or not the file exists, then open
! it and read the value of the last energy level that was simulated and written
! to that file; if this level is different from the lowest energy level we
! know that there is at least one more level to be simulated.  If it is equal,
! then we can abort the program here.

  inquire(file=trim(outname), exist=f_exists)
  if (.not.f_exists) then
    call Message%printError('EBSDmaster','restart HDF5 file does not exist')
  end if

!=============================================
! open the existing HDF5 file
!=============================================
  datagroupname = 'EBSDmaster'
  HDF = HDF_T()

! Create a new file using the default properties.
  readonly = .TRUE.
  hdferr =  HDF%openFile(outname, readonly)

! all we need to get from the file is the lastEnergy parameter
groupname = SC_EMData
  hdferr = HDF%openGroup(groupname)
  hdferr = HDF%openGroup(datagroupname)

dataset = SC_lastEnergy
  call HDF%readDatasetInteger(dataset, hdferr, lastEnergy)

  call HDF%pop(.TRUE.)
end if
!=============================================
!=============================================
! crystallography section;
verbose = .TRUE.

call cell%setFileName(mcnl%xtalname)
call Diff%setrlpmethod('WK')

if (emnl%restart.eqv..TRUE.) then
  call Diff%setV(dble(EkeVs(lastEnergy-1)))
  call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, emnl%dmin, verbose, useHDF=HDF)
else
  call Diff%setV(dble(mcnl%EkeV))
  call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, emnl%dmin, verbose, useHDF=HDF)
end if

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

! then calculate density, average atomic number and average atomic weight
call cell%calcDensity(Z2percent)
io_real(1:3) = cell%getDensity()
density = io_real(1)
at_wt = io_real(2)
Ze = io_real(3)
call Message%WriteValue('Density, avA, avZ = ',io_real,3,"(2f10.5,',',f10.5)")

! allocate and compute the Sgh loop-up table
 numset = cell%getNatomtype()
 call Diff%Initialize_SghLUT(cell, SG, emnl%dmin, numset, nat, verbose)

! determine the point group number
 j=0
 do i=1,32
  if (SGPG(i).le.SG%getSpaceGroupNumber()) j=i
 end do
 isym = j

! here is new code dealing with all the special cases (quite a few more compared to the
! Laue group case)...  isym is the point group number. Once the symmetry case has been
! fully determined (taking into account things like 31m and 3m1 an such), then the only places
! that symmetry is handled are the modified Calckvectors routine, and the filling of the modified
! Lambert projections after the dynamical simulation step.  We are also changing the name of the
! sr array (or srhex) to mLPNH and mLPSH (modified Lambert Projection Northern/Southern Hemisphere),
! and we change the output HDF5 file a little as well. We need to make sure that the EMEBSD program
! issues a warning when an old format HDF5 file is read.

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

write (*,*) '========================'
write (*,*) 'isym = ',isym
write (*,*) 'SamplingType = ', SamplingType
write (*,*) '========================'

! ---------- end of symmetry and crystallography section
!=============================================
!=============================================

!=============================================
!=============================================
! this is where we determine the value for the thickness integration limit for the CalcLgh3 routine...


! then, for each energy determine the 95% histogram thickness
izzmax = 0
do iE = 1,numEbins
 do ix=-nsx/10,nsx/10
  do iy=-nsy/10,nsy/10
   istat = sum(accum_z(iE,:,ix,iy))
   izz = 1
   do while (sum(accum_z(iE,1:izz,ix,iy)).lt.(0.99*istat))
    izz = izz+1
   end do
   if (izz.gt.izzmax) izzmax = izz
  end do
 end do
 thick(iE) = dble(izzmax) * mcnl%depthstep
end do

izz = nint(maxval(thick)/mcnl%depthstep)
call mem%alloc(lambdaE, (/ numEbins, izz/), 'lambdaE')
do iE=1,numEbins
 call Diff%setV(dble(Ekevs(iE)))
 call Diff%CalcUcg(cell,(/0,0,0/))
 rlp = Diff%getrlp()
 nabsl = rlp%xgp
 do iz=1,izz
  lambdaE(iE,iz) = float(sum(accum_z(iE,iz,-nsx/10:nsx/10,-nsy/10:nsy/10)))/etotal
  lambdaE(iE,iz) = lambdaE(iE,iz) * exp(2.0*sngl(cPi)*(iz-1)*mcnl%depthstep/nabsl)
 end do
end do

! ---------- end of 'read Monte Carlo output file and extract necessary parameters' section
!=============================================
!=============================================

!=============================================
!=============================================
! ---------- a couple of initializations
   npy = emnl%npx
   
   gzero = 1  ! index of incident beam
   debug = 0  ! no longer used
! ----------
!=============================================
!=============================================

!=============================================
!=============================================
! if the combinesites parameter is .TRUE., then we only need to
! allocate a dimension of 1 in the master pattern array since we are adding
! together the master patterns for all sites in the asymmetric unit.
  if (emnl%combinesites.eqv..TRUE.) then
    numsites = 1
  else
    numsites = numset
  end if

! ---------- allocate memory for the master patterns
  call mem%alloc(mLPNH, (/ emnl%npx, npy, 1, numsites /), 'mLPNH', startdims=(/-emnl%npx,-npy,1,1/))
  call mem%alloc(mLPSH, (/ emnl%npx, npy, 1, numsites /), 'mLPNH', startdims=(/-emnl%npx,-npy,1,1/))
  call mem%alloc(masterSPNH, (/ emnl%npx, npy, 1 /), 'masterSPNH', startdims=(/-emnl%npx,-npy,1/))
  call mem%alloc(masterSPSH, (/ emnl%npx, npy, 1 /), 'masterSPSH', startdims=(/-emnl%npx,-npy,1/))

! set various arrays to zero
  if (emnl%uniform.eqv..TRUE.) then
   mLPNH = 1.0
   mLPSH = 1.0
   masterSPNH = 1.0
   masterSPSH = 1.0
  else
   mLPNH = 0.0
   mLPSH = 0.0
   masterSPNH = 0.0
   masterSPSH = 0.0
  end if
! ---------- end allocate memory for the master patterns
!=============================================
!=============================================

! force dynamical matrix routine to read new Bethe parameters from file
! this will all be changed with the new version of the Bethe potentials
  call Diff%SetBetheParameters(EMsoft, .FALSE., emnl%BetheParametersFile)

if (emnl%restart.eqv..FALSE.) then
!=============================================
! create or update the HDF5 output file
!=============================================
  HDF = HDF_T()

! Open an existing file or create a new file using the default properties.
  if (trim(energyfile).eq.trim(outname)) then
    hdferr =  HDF%openFile(outname)
  else
    hdferr =  HDF%createFile(outname)
  end if

! write the EMheader to the file
  datagroupname = trim(HDFnames%get_ProgramData())
  call HDF%writeEMheader(EMsoft,dstr, tstrb, tstre, progname, datagroupname)

! add the Duration field to the EMheader group
  hdferr = HDF%openGroup(HDFnames%get_EMheader())
  hdferr = HDF%openGroup(HDFnames%get_ProgramData())

dataset = SC_Duration
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  tstop = 0
  if (g_exists) then
    hdferr = HDF%writeDatasetFloat(dataset, tstop, overwrite)
  else
    hdferr = HDF%writeDatasetFloat(dataset, tstop)
  end if
  call HDF%pop()
  call HDF%pop()

! open or create a namelist group to write all the namelist files into
groupname = SC_NMLfiles
  hdferr = HDF%createGroup(HDFnames%get_NMLfiles())

! read the text file and write the array to the file
dataset = trim(HDFnames%get_NMLfilename())
  hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%nmldeffile)

! leave this group
  call HDF%pop()

! create a namelist group to write all the namelist files into
  hdferr = HDF%createGroup(HDFnames%get_NMLparameters())

  call MPFT%writeHDFNameList(HDF, HDFnames, emnl)
  call Diff%writeBetheparameterNameList(HDF)

! leave this group
  call HDF%pop()

! then the remainder of the data in a EMData group
  hdferr = HDF%openGroup(HDFnames%get_EMData())
  hdferr = HDF%createGroup(HDFnames%get_ProgramData())

! create the EBSDmaster group and add a HDF_FileVersion attribute to it
  HDF_FileVersion = '4.0'
  HDF_FileVersion = cstringify(HDF_FileVersion)
  attributename = SC_HDFFileVersion
  call H5Aexists_f(HDF%getobjectID(),trim(attributename),g_exists, hdferr)
  if (.not.g_exists) then
    hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)
  end if

! =====================================================
! The following write commands constitute HDF_FileVersion = 4.0
! =====================================================
dataset = SC_xtalname
  allocate(stringarray(1))
  stringarray(1)= trim(mcnl%xtalname)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1, overwrite)
  else
    hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)
  end if

dataset = SC_numset
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetInteger(dataset, numset, overwrite)
  else
    hdferr = HDF%writeDatasetInteger(dataset, numset)
  end if

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

dataset = SC_lastEnergy
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetInteger(dataset, lastEnergy, overwrite)
  else
    hdferr = HDF%writeDatasetInteger(dataset, lastEnergy)
  end if

dataset = 'Z2percent'  ! SC_Z2percent
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetFloatArray(dataset, Z2percent, cell%getNatomtype(), overwrite)
  else
    hdferr = HDF%writeDatasetFloatArray(dataset, Z2percent, cell%getNatomtype())
  end if

  if (emnl%Esel.eq.-1) then
dataset = SC_numEbins
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then
      hdferr = HDF%writeDatasetInteger(dataset, numEbins, overwrite)
    else
      hdferr = HDF%writeDatasetInteger(dataset, numEbins)
    end if

dataset = SC_EkeVs
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then
      hdferr = HDF%writeDatasetFloatArray(dataset, EkeVs, numEbins, overwrite)
    else
      hdferr = HDF%writeDatasetFloatArray(dataset, EkeVs, numEbins)
    end if
  else
dataset = SC_numEbins
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then
      hdferr = HDF%writeDatasetInteger(dataset, one, overwrite)
    else
      hdferr = HDF%writeDatasetInteger(dataset, one)
    end if

dataset = SC_selE
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then
      hdferr = HDF%writeDatasetFloat(dataset, selE, overwrite)
    else
      hdferr = HDF%writeDatasetFloat(dataset, selE)
    end if
  end if

! create the hyperslabs and write zeroes to them for now
dataset = SC_mLPNH
  dims4 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins, numsites /)
  cnt4 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1, numsites /)
  offset4 = (/ 0, 0, 0, 0 /)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims4, offset4, cnt4, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims4, offset4, cnt4)
  end if

dataset = SC_mLPSH
  dims4 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins, numsites /)
  cnt4 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1, numsites /)
  offset4 = (/ 0, 0, 0, 0 /)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims4, offset4, cnt4, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims4, offset4, cnt4)
  end if

dataset = SC_masterSPNH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1 /)
  offset3 = (/ 0, 0, 0 /)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPNH, dims3, offset3, cnt3, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPNH, dims3, offset3, cnt3)
  end if

dataset = SC_masterSPSH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1 /)
  offset3 = (/ 0, 0, 0 /)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPSH, dims3, offset3, cnt3, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPSH, dims3, offset3, cnt3)
  end if

! =====================================================
! end of HDF_FileVersion = 4.0 write statements
! =====================================================

  call HDF%pop(.TRUE.)
end if

!=============================================
!=============================================
! do we need to precompute the Legendre array for the new lattitudinal grid values?
if (doLegendre.eqv..TRUE.) then
  call Message%printMessage(' Computing Legendre lattitudinal grid values')
  call mem%alloc(diagonal, (/2*emnl%npx+1/), 'diagonal')
  call mem%alloc(upd, (/2*emnl%npx+1/), 'upd')
  diagonal = 0.D0
  upd = (/ (dble(i) / dsqrt(4.D0 * dble(i)**2 - 1.D0), i=1,2*emnl%npx+1) /)
  call dsterf(2*emnl%npx-1, diagonal, upd, info)
! the eigenvalues are stored from smallest to largest and we need them in the opposite direction
  call mem%alloc(LegendreArray, (/2*emnl%npx+1/), 'LegendreArray')
  LegendreArray(1:2*emnl%npx+1) = diagonal(2*emnl%npx+1:1:-1)
! set the center eigenvalue to 0
  LegendreArray(emnl%npx+1) = 0.D0
  call mem%dealloc(diagonal, 'diagonal')
  call mem%dealloc(upd, 'upd')
end if

!=============================================
!=============================================
! figure out what the start energy value is for the energyloop
if (lastEnergy.ne.-1) then
  Estart = lastEnergy-1
  if (Estart.eq.0) then
    call Message%printMessage('All energy levels are present in the HDF5 file')
    call Message%printMessage('No further computations needed.')
    stop 'Program halted.'
  end if
else
  Estart = numEbins
end if

!=============================================
!=============================================
! ---------- from here on, we need to repeat the entire computation for each energy value
! so this is where we could in principle implement an OpenMP approach; alternatively,
! we could do the inner loop over the incident beam directions in OpenMP (probably simpler)

! we use two times, one (1) for each individual energy level, the other (2) for the overall time
call timer%Time_tick(2)
reflist = gvectors_T()

! instantiate the memory class for the OpenMP section
memth = memory_T( nt = emnl%nthreads, silent=.TRUE. )

energyloop: do iE=Estart,1,-1
 if (emnl%uniform.eqv..FALSE.) then
! is this a single-energy run ?
   if (emnl%Esel.ne.-1) then
     if (emnl%Esel.ne.iE) CYCLE energyloop
   end if
   ! start the energy level timer
   call timer%Time_tick(1)

! is there any intensity in the Monte Carlo lambda function ?  If not, then skip this energylevel
   if (isTKD.eqv..TRUE.) then
     if (sum(lambdaE(iE,:)).lt.1.0e-7) then
       io_int(1) = iE
       call Message%WriteValue('There are very few electrons in energy bin ',io_int,1)
       call Message%printMessage('---> Skipping this energy level ')
       CYCLE energyloop
     end if
   end if

! print a message to indicate where we are in the computation
   io_int(1)=iE
   io_int(2)=Estart
   call Message%printMessage(' Starting computation for energy bin (in reverse order)', frm = "(/A)",advance="no")
   call Message%WriteValue(' ',io_int,2,"(I4,' of ',I4)",advance="no")
   io_real(1) = EkeVs(iE)
   call Message%WriteValue('; energy [keV] = ',io_real,1,"(F6.2/)")
   selE = EkeVs(iE)

! set the accelerating voltage
   call cell%setFileName(mcnl%xtalname)
   call Diff%setV(dble(EkeVs(iE)))
   call Diff%setrlpmethod('WK')
   if(iE .ne. Estart) then
    verbose = .TRUE.
!   call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, emnl%dmin, useHDF=HDF, noLUT=.TRUE., verbose=verbose)
    call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, emnl%dmin, useHDF=HDF, verbose=verbose)
   end if

!=============================================
! ---------- create the incident beam directions list
! determine all independent incident beam directions (use a linked list starting at khead)
! numk is the total number of k-vectors to be included in this computation;
! note that this needs to be redone for each energy, since the wave vector changes with energy
   kvec = kvectors_T()   ! initialize the wave vector list
   call kvec%set_kinp( (/ 0.D0, 0.D0, 1.D0 /) )
   call kvec%set_ktmax( 0.D0 )
   call kvec%set_SamplingType( SamplingType )

   if (doLegendre.eqv..FALSE.) then
    call kvec%set_mapmode('RoscaLambert')
    if (usehex) then
      call kvec%Calckvectors(cell, SG, Diff, (/ 0.D0, 0.D0, 0.D0 /),emnl%npx,npy, ijmax,usehex)
    else
      call kvec%Calckvectors(cell, SG, Diff, (/ 0.D0, 0.D0, 0.D0 /),emnl%npx,npy, ijmax,usehex)
    end if
   else
    call kvec%set_mapmode('RoscaLambertLegendre')
    if (usehex) then
      call kvec%Calckvectors(cell, SG, Diff, (/ 0.D0, 0.D0, 0.D0 /),emnl%npx,npy, ijmax,usehex, LegendreArray)
    else
      call kvec%Calckvectors(cell, SG, Diff, (/ 0.D0, 0.D0, 0.D0 /),emnl%npx,npy, ijmax,usehex, LegendreArray)
    end if
   end if
   numk = kvec%get_numk()
   io_int(1)=numk
   call Message%WriteValue('# independent beam directions to be considered = ', io_int, 1, "(I8)")

! convert part of the kvector linked list into arrays for OpenMP
  call mem%alloc(karray, (/4,numk/), 'karray')
  call mem%alloc(kij, (/3,numk/), 'kij')
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

! set the number of OpenMP threads
  call OMP_setNThreads(emnl%nthreads)

! use OpenMP to run on multiple cores ...
!$OMP PARALLEL COPYIN(rlp) &
!$OMP& PRIVATE(DynMat,Sgh,Lgh,ik,FN,TID,kn,ipx,ipy,ipz,ix,iequiv,nequiv,reflist,firstw) &
!$OMP& PRIVATE(kkk,nns,nnw,nref,svals,io_int)

  NUMTHREADS = OMP_GET_NUM_THREADS()
  TID = OMP_GET_THREAD_NUM()

  call memth%alloc(svals, (/ numset /), 'svals', TID=TID)

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

     if (self%nml%kinematical.eqv..FALSE.) then
! generate the dynamical matrix
       ! if (allocated(DynMat)) deallocate(DynMat)
       call memth%alloc(DynMat, (/nns,nns/), 'DynMat', TID=TID)
       call reflist%GetDynMat(cell, Diff, firstw, DynMat, nns, nnw)
       totstrong = totstrong + nns
       totweak = totweak + nnw
     else
! all reflections are strong, but they are not coupled to each other, only to the
! incident beam; all q_{g-g'} are zero except the ones with g'=0.  In addition, there
! is no anomalous absorption, only normal absorption.
       ! if (allocated(DynMat)) deallocate(DynMat)
       call memth%alloc(DynMat, (/nns,nns/), 'DynMat', TID=TID)
       call reflist%GetDynMatKin(cell, Diff, firstw, DynMat, nns)
       totstrong = totstrong + nns
       totweak = 0
     end if

! then we need to initialize the Sgh and Lgh arrays
     ! if (allocated(Sgh)) deallocate(Sgh)
     ! if (allocated(Lgh)) deallocate(Lgh)
     call memth%alloc(Sgh, (/ nns,nns,numset /), 'Sgh', TID=TID)
     call memth%alloc(Lgh, (/ nns,nns /), 'Lgh', TID=TID)
     Sgh = czero
     Lgh = czero
     call reflist%getSghfromLUT(Diff,nns,numset,Sgh)


! solve the dynamical eigenvalue equation for this beam direction
     kn = karray(4,ik)
     call reflist%CalcLgh(DynMat,Lgh,dble(thick(iE)),dble(kn),nns,gzero,mcnl%depthstep,lambdaE(iE,1:izzmax),izzmax)

! sum over the element-wise (Hadamard) product of the Lgh and Sgh arrays
     svals = 0.0
     do ix=1,numset
       svals(ix) = real(sum(Lgh(1:nns,1:nns)*Sgh(1:nns,1:nns,ix)))
     end do
     svals = svals/float(sum(nat(1:numset)))

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
  if (self%nml%combinesites.eqv..FALSE.) then
     do ix=1,nequiv
       if (iequiv(3,ix).eq.-1) mLPSH(iequiv(1,ix),iequiv(2,ix),1,1:numset) = svals(1:numset)
       if (iequiv(3,ix).eq.1) mLPNH(iequiv(1,ix),iequiv(2,ix),1,1:numset) = svals(1:numset)
     end do
  else
     do ix=1,nequiv
       if (iequiv(3,ix).eq.-1) mLPSH(iequiv(1,ix),iequiv(2,ix),1,1) = sum(svals)
       if (iequiv(3,ix).eq.1) mLPNH(iequiv(1,ix),iequiv(2,ix),1,1) = sum(svals)
     end do
  end if
!$OMP END CRITICAL
     if (mod(ik,5000).eq.0) then
       io_int(1) = ik
       io_int(2) = numk
       call Message%WriteValue('  completed beam direction ',io_int, 2, "(I8,' of ',I8)")
     end if

     call reflist%Delete_gvectorlist()
     call memth%dealloc(Sgh, 'Sgh', TID=TID)
     call memth%dealloc(Lgh, 'Sgh', TID=TID)
     call memth%dealloc(DynMat, 'DynMat', TID=TID)

    end do beamloop
    
  call memth%dealloc(svals, 'svals', TID=TID)

! end of OpenMP portion
!$OMP END PARALLEL

! deallocate arrays that will need to be re-allocated in the next cycle
  call mem%dealloc(karray, 'karray')
  call mem%dealloc(kij, 'kij')

  if (usehex) then
! and finally, we convert the hexagonally sampled array to a square Lambert projection which will be used
! for all EBSD pattern interpolations;  we need to do this for both the Northern and Southern hemispheres

! we begin by allocating auxiliary arrays to hold copies of the hexagonal data; the original arrays will
! then be overwritten with the newly interpolated data.
    call mem%alloc(auxNH, (/emnl%npx,npy,1,numsites/), 'auxNH', startdims=(/-emnl%npx,-npy,1,1/))
    call mem%alloc(auxSH, (/emnl%npx,npy,1,numsites/), 'auxSH', startdims=(/-emnl%npx,-npy,1,1/))
    ! allocate(auxSH(-emnl%npx:emnl%npx,-npy:npy,1,1:numsites),stat=istat)
    auxNH = mLPNH
    auxSH = mLPSH

    edge = 1.D0 / dble(emnl%npx)
    scl = float(emnl%npx)
    do i=-emnl%npx,emnl%npx
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
          if (nixp.gt.emnl%npx) nixp = nix
          if (niyp.gt.emnl%npx) niyp = niy
          dx = xy(1) - nix
          dy = xy(2) - niy
          dxm = 1.D0 - dx
          dym = 1.D0 - dy
          mLPNH(i,j,1,1:numsites) = auxNH(nix,niy,1,1:numsites)*dxm*dym + auxNH(nixp,niy,1,1:numsites)*dx*dym + &
                               auxNH(nix,niyp,1,1:numsites)*dxm*dy + auxNH(nixp,niyp,1,1:numsites)*dx*dy
          mLPSH(i,j,1,1:numsites) = auxSH(nix,niy,1,1:numsites)*dxm*dym + auxSH(nixp,niy,1,1:numsites)*dx*dym + &
                               auxSH(nix,niyp,1,1:numsites)*dxm*dy + auxSH(nixp,niyp,1,1:numsites)*dx*dy
        end if
      end do
    end do
    call mem%dealloc(auxNH, 'auxNH')
    call mem%dealloc(auxSH, 'auxSH')
  end if

! make sure that the outer pixel rim of the mLPSH patterns is identical to
! that of the mLPNH array.
  mLPSH(-emnl%npx,-emnl%npx:emnl%npx,1,1:numsites) = mLPNH(-emnl%npx,-emnl%npx:emnl%npx,1,1:numsites)
  mLPSH( emnl%npx,-emnl%npx:emnl%npx,1,1:numsites) = mLPNH( emnl%npx,-emnl%npx:emnl%npx,1,1:numsites)
  mLPSH(-emnl%npx:emnl%npx,-emnl%npx,1,1:numsites) = mLPNH(-emnl%npx:emnl%npx,-emnl%npx,1,1:numsites)
  mLPSH(-emnl%npx:emnl%npx, emnl%npx,1,1:numsites) = mLPNH(-emnl%npx:emnl%npx, emnl%npx,1,1:numsites)

! get stereographic projections (summed over the atomic positions)
  Radius = 1.0
  do i=-emnl%npx,emnl%npx
    do j=-emnl%npx,emnl%npx
      L = Lambert_T( xyd = (/ dble(i), dble(j) /) / dble(emnl%npx) )
      ierr = L%StereoGraphicInverse( xyz, Radius )
      xyz = xyz/vecnorm(xyz)
      if (ierr.ne.0) then
        masterSPNH(i,j,1) = 0.0
        masterSPSH(i,j,1) = 0.0
      else
        masterSPNH(i,j,1) = InterpolateLambert(xyz, mLPNH, emnl%npx, numsites)
        masterSPSH(i,j,1) = InterpolateLambert(xyz, mLPSH, emnl%npx, numsites)
      end if
    end do
  end do

! since these computations can take a long time, here we store
! all the output at the end of each pass through the energyloop.

  io_int(1) = nint(float(totstrong)/float(numk))
  call Message%WriteValue(' -> Average number of strong reflections = ',io_int, 1, "(I5)")
  io_int(1) = nint(float(totweak)/float(numk))
  call Message%WriteValue(' -> Average number of weak reflections   = ',io_int, 1, "(I5)")

  call timer%makeTimeStamp()
  dstr = timer%getDateString()
  tstre = timer%getTimeString()
 end if  ! (emnl%uniform.eqv..FALSE.)

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

! update the current energy level counter, so that the restart option will function
dataset = SC_lastEnergy
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDataSetInteger(dataset, iE, overwrite)
  else
    hdferr = HDF%writeDataSetInteger(dataset, iE)
  end if

! add data to the hyperslab
dataset = SC_mLPNH
  dims4 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins, numsites /)
  cnt4 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1, numsites /)
  offset4 = (/ 0, 0, iE-1, 0 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims4, offset4, cnt4, insert)

dataset = SC_mLPSH
  dims4 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins, numsites /)
  cnt4 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1, numsites /)
  offset4 = (/ 0, 0, iE-1, 0 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims4, offset4, cnt4, insert)

dataset = SC_masterSPNH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1 /)
  offset3 = (/ 0, 0, iE-1 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPNH, dims3, offset3, cnt3, insert)

dataset = SC_masterSPSH
  dims3 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins /)
  cnt3 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1 /)
  offset3 = (/ 0, 0, iE-1 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPSH, dims3, offset3, cnt3, insert)

  call HDF%pop(.TRUE.)

 if ((emnl%Esel.eq.-1).and.(iE.ne.1)) then
  call Message%printMessage(' Intermediate data stored in file '//trim(emnl%energyfile), frm = "(A/)")
 end if

 if ((emnl%Esel.eq.-1).and.(iE.eq.1)) then
  call Message%printMessage(' Final data stored in file '//trim(emnl%energyfile), frm = "(A/)")
 end if

end do energyloop

call timer%Time_tock(2)
io_int(1) = timer%getInterval(2)
call Message%WriteValue(' Total execution time [s] ',io_int,1)

! if requested, we notify the user that this program has completed its run
if (trim(EMsoft%getConfigParameter('EMNotify')).ne.'Off') then
  if (trim(emnl%Notify).eq.'On') then
    NumLines = 3
    allocate(MessageLines(NumLines))

    call hostnm(c)
    if (isTKD.eqv..TRUE.) then
      MessageLines(1) = 'EMTKDmaster program has ended successfully'
    else
      MessageLines(1) = 'EMEBSDmaster program has ended successfully'
    end if
    MessageLines(2) = 'Master pattern data stored in '//trim(outname)
    write (exectime,"(F10.4)") timer%getInterval(2)
    MessageLines(3) = 'Total execution time [s]: '//trim(exectime)
    SlackUsername = 'EMsoft on '//trim(c)
    i = PostMessage(EMsoft, MessageLines, NumLines, SlackUsername)
  end if
end if

end associate

call mem%dealloc(EkeVs, 'EkeVs')
call mem%dealloc(thick, 'thick')
call mem%dealloc(lambdaE, 'lambdaE')
call mem%dealloc(mLPNH, 'mLPNH')
call mem%dealloc(mLPSH, 'mLPSH')
call mem%dealloc(masterSPNH, 'masterSPNH')
call mem%dealloc(masterSPSH, 'masterSPSH')

! call mem%allocated_memory_use()
! call memth%thread_memory_use()

end subroutine EBSDmaster_


!--------------------------------------------------------------------------
subroutine testEBSDmasterWrapper_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: testEBSDmasterWrapper_
!! author: MDG
!! version: 1.0
!! date: 04/16/21
!!
!! bare bones testing routine for the C-callable EMsoftCgetEBSDmaster routine in mod_SEMwrappers
!!

use mod_EMsoft
use mod_initializers
use mod_symmetry
use mod_crystallography
use mod_gvectors
use mod_kvectors
use mod_io
use mod_math
use mod_diffraction
use mod_timing
use mod_Lambert
use HDF5
use mod_HDFsupport
use mod_HDFnames
use ISO_C_BINDING
use mod_memory
use mod_notifications
use stringconstants
use mod_MCfiles
use mod_SEMwrappers

IMPLICIT NONE

class(EBSDmaster_T), INTENT(INOUT)  :: self
type(EMsoft_T), INTENT(INOUT)       :: EMsoft
character(fnlen),INTENT(IN)         :: progname
type(HDFnames_T),INTENT(INOUT)      :: HDFnames

type(Cell_T)            :: cell
type(IO_T)              :: Message
type(DynType)           :: Dyn
type(HDF_T)             :: HDF
type(SpaceGroup_T)      :: SG
type(Diffraction_T)     :: Diff
type(MCfile_T)          :: MCFT
type(memory_T)          :: mem, memth

real(kind=dbl)          :: ctmp(192,3), arg, Radius, xyz(3)
integer(HSIZE_T)        :: dims4(4), cnt4(4), offset4(4)
integer(HSIZE_T)        :: dims3(3), cnt3(3), offset3(3)
integer(kind=irg)       :: isym,i,j,ik,npy,ipx,ipy,ipz,debug,iE,izz, izzmax, iequiv(3,48), nequiv, num_el, MCnthreads, & ! counters
                           numk, timestart, timestop, numsites, nthreads, & ! number of independent incident beam directions
                           ir,nat(maxpasym),kk(3), skip, ijmax, one, NUMTHREADS, TID, SamplingType, &
                           numset,n,ix,iy,iz, io_int(6), nns, nnw, nref, Estart, &
                           istat,gzero,ic,ip,ikk, totstrong, totweak, jh, ierr, nix, niy, nixp, niyp     ! counters
real(kind=dbl)          :: tpi,Znsq, kkl, DBWF, kin, delta, h, lambda, omtl, srt, dc(3), xy(2), edge, scl, tmp, dx, dxm, dy, dym !
real(kind=sgl)          :: io_real(5), selE, kn, FN(3), kkk(3), tstop, nabsl, etotal, density, Ze, at_wt, bp(4)
real(kind=sgl),allocatable      :: mLPNH(:,:,:,:), mLPSH(:,:,:,:)
integer(kind=irg),allocatable   :: accum_z(:,:,:,:)
complex(kind=dbl)               :: czero
logical                 :: usehex, switchmirror, verbose

! Monte Carlo derived quantities
integer(kind=irg)       :: numEbins, nsx, nsy, hdferr, nlines, lastEnergy    ! variables used in MC energy file
character(fnlen)        :: oldprogname, groupname, energyfile, outname, datagroupname, attributename, HDF_FileVersion, fname
character(8)            :: MCscversion
character(11)           :: dstr
character(15)           :: tstrb
character(15)           :: tstre
logical                 :: f_exists, readonly, overwrite=.TRUE., insert=.TRUE., stereog, g_exists, xtaldataread, FL, &
                           doLegendre, isTKD = .FALSE.

character(fnlen)                :: dataset, instring
type(MCOpenCLNameListType)      :: mcnl

real(kind=sgl), allocatable     :: atompos(:,:)
integer(kind=irg),allocatable   :: atomtypes(:)
real(kind=sgl),allocatable      :: atpos(:,:)
integer(kind=irg),allocatable   :: attp(:)
real(kind=sgl)                  :: latparm(6)
integer(c_int32_t)              :: ipar(wraparraysize)
real(kind=sgl)                  :: fpar(wraparraysize)
integer(c_size_t)               :: objAddress
character(len=1)                :: cancel

call openFortranHDFInterface()

! set the HDF group names for this program
HDF = HDF_T()
HDFnames = HDFnames_T()

! simplify the notation a little
associate( emnl => self%nml )

! initialize the memory class 
mem = memory_T()

! if copyfromenergyfile is different from 'undefined', then we need to
! copy all the Monte Carlo data from that file into a new file, which
! will then be read from and written to by the ComputeMasterPattern routine.
if (emnl%copyfromenergyfile.ne.'undefined') then
  call MCFT%copyMCdata(EMsoft, HDF, emnl%copyfromenergyfile, emnl%energyfile, emnl%h5copypath)
end if

energyfile = trim(EMsoft%generateFilePath('EMdatapathname',emnl%energyfile))
outname = trim(energyfile)

!=============================================
!=============================================
! ---------- read Monte Carlo .h5 output file and extract necessary parameters
! set the HDF group names for reading the MC input file
call HDFnames%set_ProgramData(SC_MCOpenCL)
call HDFnames%set_NMLlist(SC_MCCLNameList)
call HDFnames%set_NMLfilename(SC_MCOpenCLNML)
fname = EMsoft%generateFilePath('EMdatapathname',trim(emnl%energyfile))
call MCFT%setFileName(fname)
call MCFT%readMCfile(HDF, HDFnames, getAccumz=.TRUE.)
mcnl = MCFT%getnml()
call MCFT%copyaccumz(accum_z)

! set the HDFnames to the correct strings for this program
call HDFnames%set_ProgramData(SC_EBSDmaster)
call HDFnames%set_NMLlist(SC_EBSDmasterNameList)
call HDFnames%set_NMLfilename(SC_EBSDmasterNML)
call HDFnames%set_Variable(SC_MCOpenCL)

nsx = (mcnl%numsx - 1)/2
nsy = nsx
etotal = float(mcnl%totnum_el)
numEbins = MCFT%getnumEbins()

!=============================================
!=============================================
! crystallography section;
verbose = .TRUE.

call cell%setFileName(mcnl%xtalname)
call Diff%setrlpmethod('WK')

call Diff%setV(dble(mcnl%EkeV))
call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, emnl%dmin, verbose, useHDF=HDF)

! from here on we simply get all the parameters set properly to call the wrapper routine 
! ---------- a couple of initializations
npy = emnl%npx
numsites = cell%getNatomtype()

! ---------- allocate memory for the master patterns
call mem%alloc(mLPNH, (/ emnl%npx, npy, numEbins, numsites /), 'mLPNH', 0.0, startdims=(/-emnl%npx,-npy,1,1/))
call mem%alloc(mLPSH, (/ emnl%npx, npy, numEbins, numsites /), 'mLPNH', 0.0, startdims=(/-emnl%npx,-npy,1,1/))

! define the ipar and fpar parameter arrays 
ipar = 0
ipar(1) = nsx
ipar(2) = mcnl%globalworkgrpsz
ipar(3) = mcnl%num_el
ipar(4) = mcnl%totnum_el
ipar(5) = mcnl%multiplier
ipar(6) = mcnl%devid
ipar(7) = mcnl%platid
ipar(8) = SG%getSpaceGroupXtalSystem()
ipar(9) = numsites
ipar(10)= SG%getSpaceGroupNumber()
ipar(11)= SG%getSpaceGroupSetting()
ipar(12)= MCFT%getnumEbins()
ipar(13)= MCFT%getnumzbins()
ipar(14)= 1
ipar(15)= 1
ipar(16)= nsx/10
ipar(17)= emnl%npx
ipar(18)= emnl%nthreads
ipar(36)= 0

fpar = 0.0
fpar(1) = mcnl%sig
fpar(2) = mcnl%omega
fpar(3) = mcnl%EkeV
fpar(4) = mcnl%Ehistmin
fpar(5) = mcnl%Ebinsize
fpar(6) = mcnl%depthmax
fpar(7) = mcnl%depthstep
fpar(8) = mcnl%sigstart
fpar(9) = mcnl%sigend
fpar(10)= mcnl%sigstep
fpar(11)= emnl%dmin
fpar(12)= 4.0
fpar(13)= 8.0
fpar(14)= 50.0

! get the atom positions and lattice parameters into separate arrays 
atomtypes = cell%getAtomtype()
atompos = cell%getAsymPosData()
latparm = cell%getLatParm()
allocate(atpos(numsites,5), attp(numsites))
atpos = atompos(1:numsites,1:5)
attp = atomtypes(1:numsites) 

cancel = char(0)
objAddress = 0_c_size_t

! call the wrapper routine 
write (*,*) 'calling wrapper routine '
call EMsoftCgetEBSDmaster(ipar,fpar,atpos,attp,latparm,accum_z,mLPNH,mLPSH,C_NULL_FUNPTR,objAddress,cancel)
write (*,*) '   done.'

!=============================================
! create a simple HDF5 output file
!=============================================
HDF = HDF_T()

! Open an existing file or create a new file using the default properties.
hdferr =  HDF%openFile(outname)

! then the remainder of the data in a EMData group
hdferr = HDF%openGroup(HDFnames%get_EMData())
hdferr = HDF%createGroup(HDFnames%get_ProgramData())

! create the hyperslabs and write zeroes to them for now
dataset = SC_mLPNH
  dims4 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins, numsites /)
  cnt4 = (/ 2*emnl%npx+1, 2*emnl%npx+1, numEbins, numsites /)
  offset4 = (/ 0, 0, 0, 0 /)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims4, offset4, cnt4, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims4, offset4, cnt4)
  end if

dataset = SC_mLPSH
  dims4 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins, numsites /)
  cnt4 = (/ 2*emnl%npx+1, 2*emnl%npx+1, numEbins, numsites /)
  offset4 = (/ 0, 0, 0, 0 /)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims4, offset4, cnt4, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims4, offset4, cnt4)
  end if

  call HDF%pop(.TRUE.)

write (*,*) 'output written to '//trim(outname)

call mem%dealloc(mLPNH, 'mLPNH')
call mem%dealloc(mLPSH, 'mLPSH')

call closeFortranHDFInterface()

end associate

end subroutine testEBSDmasterWrapper_

end module mod_EBSDmaster
