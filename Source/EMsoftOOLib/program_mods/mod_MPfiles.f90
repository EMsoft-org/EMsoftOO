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

module mod_MPfiles
  !! author: MDG
  !! version: 1.0 
  !! date: 02/12/20
  !!
  !! generic class for Master Pattern file handling; programs must inherit this class
  !! and extend it, for instance for EBSD, ECP, or TKD master patterns; we pass the namelist-type 
  !! to the routines using the type(*) specification.

use mod_kinds
use mod_global
use stringconstants

IMPLICIT NONE 
private

type, public :: SEMmasterNameListType
  integer(kind=irg) :: npx
  integer(kind=irg) :: nthreads
  real(kind=sgl)    :: dmin
  character(3)      :: Notify
  character(fnlen)  :: copyfromenergyfile
  character(fnlen)  :: h5copypath
  character(fnlen)  :: energyfile
  character(fnlen)  :: BetheParametersFile
  logical           :: combinesites
  logical           :: uniform
  logical           :: kinematical 
end type SEMmasterNameListType

! inherit the SEMmasterNameListType and extend it with additional parameters
type, public, extends(SEMmasterNameListType) :: EBSDmasterNameListType
  integer(kind=irg) :: Esel
  logical           :: restart
  logical           :: useEnergyWeighting 
end type EBSDmasterNameListType

type, public, extends(SEMmasterNameListType) :: EBSDmasterSHTNameListType
  integer(kind=irg) :: Esel
  logical           :: restart
  logical           :: useEnergyWeighting 
end type EBSDmasterSHTNameListType

type, public, extends(SEMmasterNameListType) :: ECPmasterNameListType
end type ECPmasterNameListType

type, public, extends(SEMmasterNameListType) :: TKDmasterNameListType
  integer(kind=irg) :: Esel
  logical           :: restart
  real(kind=sgl)    :: thickness
end type TKDmasterNameListType


type, public :: MPdataType
  integer(kind=irg)               :: lastEnergy
  integer(kind=irg)               :: numEbins
  integer(kind=irg)               :: numset
  integer(kind=irg)               :: newPGnumber
  logical                         :: AveragedMP
  character(fnlen)                :: xtalname
  real(kind=sgl),allocatable      :: keVs(:)
  real(kind=sgl),allocatable      :: mLPNH4(:,:,:,:)
  real(kind=sgl),allocatable      :: mLPSH4(:,:,:,:)
  real(kind=sgl),allocatable      :: mLPNH(:,:,:)
  real(kind=sgl),allocatable      :: mLPSH(:,:,:)
  real(kind=sgl),allocatable      :: masterSPNH(:,:,:)
  real(kind=sgl),allocatable      :: masterSPSH(:,:,:)
end type MPdataType

type, public :: MPfile_T 
  private 
    type(MPdataType),public       :: MPDT
    character(fnlen)              :: MPfile
    character(fnlen)              :: modality = 'unknown'

  contains
  private 

    procedure, pass(self) :: readMPfile_
    procedure, pass(self) :: setFileName_
    ! procedure, pass(self) :: writeMPfile_
    procedure, pass(self) :: writeHDFNameList_
    ! procedure, pass(self) :: copynml_
    ! procedure, pass(self) :: getnml_
    procedure, pass(self) :: set_Modality_
    procedure, pass(self) :: get_Modality_
    procedure, pass(self) :: determine_Modality_
    procedure, pass(self) :: getlastEnergy_
    procedure, pass(self) :: getnumEbins_
    procedure, pass(self) :: getnumset_
    procedure, pass(self) :: getnewPGnumber_
    procedure, pass(self) :: getAveragedMP_
    procedure, pass(self) :: getxtalname_
    ! procedure, pass(self) :: copyBetheParameters_
    procedure, pass(self) :: copykeVs_
    procedure, pass(self) :: copymLPNH4_
    procedure, pass(self) :: copymLPSH4_
    procedure, pass(self) :: copymLPNH_
    procedure, pass(self) :: copymLPSH_
    procedure, pass(self) :: copysummLPNH_
    procedure, pass(self) :: copysummLPSH_
    procedure, pass(self) :: copymasterSPNH_
    procedure, pass(self) :: copymasterSPSH_
    procedure, pass(self) :: copyMPdata_
    procedure, pass(self) :: copyMPoverlapdata_ 

    generic, public :: readMPfile => readMPfile_
    generic, public :: setFileName => setFileName_
    generic, public :: setModality => set_Modality_
    generic, public :: getModality => get_Modality_
    generic, public :: determineModality => determine_Modality_
    ! generic, public :: writeMPfile => writeMPfile_
    generic, public :: writeHDFNameList => writeHDFNameList_
    ! generic, public :: copynml => copynml_
    ! generic, public :: getnml => getnml_
    generic, public :: getlastEnergy => getlastEnergy_
    generic, public :: getnumEbins => getnumEbins_
    generic, public :: getnumset => getnumset_
    generic, public :: getnewPGnumber => getnewPGnumber_
    generic, public :: getAveragedMP => getAveragedMP_
    generic, public :: getxtalname => getxtalname_
    ! generic, public :: copyBetheParameters => copyBetheParameters_
    generic, public :: copykeVs => copykeVs_
    generic, public :: copymLPNH4 => copymLPNH4_
    generic, public :: copymLPSH4 => copymLPSH4_
    generic, public :: copymLPNH => copymLPNH_
    generic, public :: copymLPSH => copymLPSH_
    generic, public :: copysummLPNH => copysummLPNH_
    generic, public :: copysummLPSH => copysummLPSH_
    generic, public :: copymasterSPNH => copymasterSPNH_
    generic, public :: copymasterSPSH => copymasterSPSH_
    generic, public :: copyMPdata => copyMPdata_
    generic, public :: copyMPoverlapdata => copyMPoverlapdata_

end type MPfile_T

! the constructor routine for this class 
interface MPfile_T
  module procedure MPfile_constructor
end interface MPfile_T

contains 

!--------------------------------------------------------------------------
type(MPfile_T) function MPfile_constructor( ) result(MPfile)
!! author: MDG 
!! version: 1.0 
!! date: 02/05/20
!!
!! constructor for the MPfile_T Class
 
IMPLICIT NONE

end function MPfile_constructor

! !--------------------------------------------------------------------------
! subroutine copynml_(self, nml)
! !! author: MDG 
! !! version: 1.0 
! !! date: 02/17/20
! !!
! !! copy the namelist into the MPfile_T class for writing to file

! IMPLICIT NONE 

! class(MPfile_T), INTENT(INOUT)            :: self
! type(EBSDmasterNameListType),INTENT(IN)   :: nml

! self%nml = nml

! end subroutine copynml_

!--------------------------------------------------------------------------
subroutine copykeVs_(self, acc, keep)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! copy the keVs array

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)                :: self
real(kind=sgl), allocatable, INTENT(OUT)      :: acc(:)
logical, INTENT(IN), OPTIONAL                 :: keep 

integer(kind=irg)                             :: s(1)

s = shape(self%MPDT%keVs)
allocate(acc(s(1)))

acc = self%MPDT%keVs 

if (.not.present(keep)) deallocate(self%MPDT%keVs)

end subroutine copykeVs_

!--------------------------------------------------------------------------
subroutine copymLPNH4_(self, acc, keep)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! copy the mLPNH4 array

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)                :: self
real(kind=sgl), allocatable, INTENT(OUT)      :: acc(:,:,:,:)
logical, INTENT(IN), OPTIONAL                 :: keep 

integer(kind=irg)                             :: s(4), nx

s = shape(self%MPDT%mLPNH4)
nx = (s(1)-1)/2
allocate(acc(-nx:nx,-nx:nx,s(3),s(4)))

acc = self%MPDT%mLPNH4 

if (.not.present(keep)) deallocate(self%MPDT%mLPNH4)

end subroutine copymLPNH4_

!--------------------------------------------------------------------------
subroutine copymLPSH4_(self, acc, keep)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! copy the mLPSH4 array

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)                :: self
real(kind=sgl), allocatable, INTENT(OUT)      :: acc(:,:,:,:)
logical, INTENT(IN), OPTIONAL                 :: keep 

integer(kind=irg)                             :: s(4), nx

s = shape(self%MPDT%mLPSH4)
nx = (s(1)-1)/2
allocate(acc(-nx:nx,-nx:nx,s(3),s(4)))

acc = self%MPDT%mLPSH4 

if (.not.present(keep)) deallocate(self%MPDT%mLPSH4)

end subroutine copymLPSH4_

!--------------------------------------------------------------------------
subroutine copymLPNH_(self, acc, keep)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! copy the mLPNH array

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)                :: self
real(kind=sgl), allocatable, INTENT(OUT)      :: acc(:,:,:)
logical, INTENT(IN), OPTIONAL                 :: keep 

integer(kind=irg)                             :: s(3), nx

s = shape(self%MPDT%mLPNH)
nx = (s(1)-1)/2
allocate(acc(-nx:nx,-nx:nx,s(3)))

acc = self%MPDT%mLPNH 

if (.not.present(keep)) deallocate(self%MPDT%mLPNH)

end subroutine copymLPNH_

!--------------------------------------------------------------------------
subroutine copymLPSH_(self, acc, keep)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! copy the mLPSH array

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)                :: self
real(kind=sgl), allocatable, INTENT(OUT)      :: acc(:,:,:)
logical, INTENT(IN), OPTIONAL                 :: keep 

integer(kind=irg)                             :: s(3), nx

s = shape(self%MPDT%mLPSH)
nx = (s(1)-1)/2
allocate(acc(-nx:nx,-nx:nx,s(3)))

acc = self%MPDT%mLPSH 

if (.not.present(keep)) deallocate(self%MPDT%mLPSH)

end subroutine copymLPSH_

!--------------------------------------------------------------------------
subroutine copysummLPNH_(self, acc, keep)
!! author: MDG 
!! version: 1.0 
!! date: 03/15/20
!!
!! copy the mLPNH array and return it summed over the last dimension

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)                :: self
real(kind=sgl), allocatable, INTENT(OUT)      :: acc(:,:)
logical, INTENT(IN), OPTIONAL                 :: keep 

integer(kind=irg)                             :: s(3), nx

s = shape(self%MPDT%mLPNH)
nx = (s(1)-1)/2
allocate( acc(-nx:nx, -nx:nx) )

acc = sum(self%MPDT%mLPNH, 3)

if (.not.present(keep)) deallocate(self%MPDT%mLPNH)

end subroutine copysummLPNH_

!--------------------------------------------------------------------------
subroutine copysummLPSH_(self, acc, keep)
!! author: MDG 
!! version: 1.0 
!! date: 03/15/20
!!
!! copy the mLPSH array and return it summed over the last dimension

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)                :: self
real(kind=sgl), allocatable, INTENT(OUT)      :: acc(:,:)
logical, INTENT(IN), OPTIONAL                 :: keep 

integer(kind=irg)                             :: s(3), nx

s = shape(self%MPDT%mLPSH)
nx = (s(1)-1)/2
allocate( acc(-nx:nx, -nx:nx) )

acc = sum(self%MPDT%mLPSH, 3)

if (.not.present(keep)) deallocate(self%MPDT%mLPSH)

end subroutine copysummLPSH_

!--------------------------------------------------------------------------
subroutine copymasterSPNH_(self, acc, keep)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! copy the masterSPNH array

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)                :: self
real(kind=sgl), allocatable, INTENT(OUT)      :: acc(:,:,:)
logical, INTENT(IN), OPTIONAL                 :: keep 

integer(kind=irg)                             :: s(3), nx

s = shape(self%MPDT%masterSPNH)
nx = (s(1)-1)/2
allocate(acc(-nx:nx,-nx:nx,s(3)))

acc = self%MPDT%masterSPNH 

if (.not.present(keep)) deallocate(self%MPDT%masterSPNH)

end subroutine copymasterSPNH_

!--------------------------------------------------------------------------
subroutine copymasterSPSH_(self, acc, keep)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! copy the masterSPSH array

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)                :: self
real(kind=sgl), allocatable, INTENT(OUT)      :: acc(:,:,:)
logical, INTENT(IN), OPTIONAL                 :: keep 

integer(kind=irg)                             :: s(3), nx

s = shape(self%MPDT%masterSPSH)
nx = (s(1)-1)/2
allocate(acc(-nx:nx,-nx:nx,s(3)))

acc = self%MPDT%masterSPSH 

if (.not.present(keep)) deallocate(self%MPDT%masterSPSH)

end subroutine copymasterSPSH_

!--------------------------------------------------------------------------
function get_Modality_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/03/20
!!
!! get Modality from the MPfile_T class

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = self%Modality

end function get_Modality_

!--------------------------------------------------------------------------
subroutine set_Modality_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/03/20
!!
!! set Modality in the MPfile_T class

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)     :: self
character(*), INTENT(IN)           :: inp

self%Modality = inp

end subroutine set_Modality_

!--------------------------------------------------------------------------
subroutine setFileName_(self, MPfile)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! set the Master Pattern file name 

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)          :: self
character(fnlen), INTENT(IN)            :: MPfile

self%MPfile = trim(MPfile)

end subroutine setFileName_

!--------------------------------------------------------------------------
subroutine determine_Modality_(self, HDF, MPfile)
!! author: MDG 
!! version: 1.0 
!! date: 03/24/20
!!
!! determine what type of Master Pattern file this is 

use HDF5 
use mod_HDFsupport 
use mod_io
use stringconstants 

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)    :: self
type(HDF_T), INTENT(INOUT)        :: HDF
character(fnlen), INTENT(IN)      :: MPfile 

type(IO_T)                        :: Message 
character(fnlen)                  :: groupname 
logical                           :: f_exists, g_exists, stat 
integer(kind=irg)                 :: hdferr

! we assume that MPfile contains the full path to the master pattern file 
inquire(file=trim(MPfile), exist=f_exists)

if (.not.f_exists) then
  call Message%printError('determineModality','Master Pattern file '//trim(MPfile)//' does not exist')
end if

! is this a proper HDF5 file ?
call h5fis_hdf5_f(trim(MPfile), stat, hdferr)

! open the file 
hdferr =  HDF%openFile(MPfile) 

! go to the NMLfiles group and see what's there ... 
groupname = SC_NMLfiles 
hdferr = HDF%opengroup(groupname)

groupname = SC_EBSDmasterNML
call H5Lexists_f(HDF%getobjectID(),trim(groupname),g_exists, hdferr)
if (g_exists) then
  call self%set_Modality_('EBSD')
else 
  groupname = SC_ECPmasterNML
  call H5Lexists_f(HDF%getobjectID(),trim(groupname),g_exists, hdferr)
  if (g_exists) then
    call self%set_Modality_('ECP')
  else 
    groupname = SC_TKDmasterNML
    call H5Lexists_f(HDF%getobjectID(),trim(groupname),g_exists, hdferr)
    if (g_exists) then
      call self%set_Modality_('TKD')
    else 
      groupname = SC_KosselmasterNML
      call H5Lexists_f(HDF%getobjectID(),trim(groupname),g_exists, hdferr)
      if (g_exists) then
        call self%set_Modality_('Kossel')
      else 
        call self%set_Modality_('unknown')
      end if 
    end if 
  end if 
end if 

! close the file 
call HDF%pop(.TRUE.)

end subroutine determine_Modality_

! !--------------------------------------------------------------------------
! function getnml_(self) result(nml)
! !! author: MDG 
! !! version: 1.0 
! !! date: 02/17/20
! !!
! !! get the namelist from the MPfile_T class

! IMPLICIT NONE 

! class(MPfile_T), INTENT(INOUT)            :: self
! class(SEMmasterNameListType)              :: nml

! select type(nml)
!   type is (EBSDmasterNameListType)
!     nml = self%EBSDnml
!   type is (EBSDmasterSHTNameListType)
!     nml = self%EBSDSHTnml
!   type is (ECPmasterNameListType)
!     nml = self%ECPnml
!   type is (TKDmasterNameListType)
!     nml = self%TKDnml
! end select

! end function getnml_

!--------------------------------------------------------------------------
function getnumEbins_(self) result(n)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! get numEbins from the MPfile_T class

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)    :: self
integer(kind=irg)                 :: n

n = self%MPDT%numEbins

end function getnumEbins_

!--------------------------------------------------------------------------
function getnumset_(self) result(n)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! get numset from the MPfile_T class

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)    :: self
integer(kind=irg)                 :: n

n = self%MPDT%numset

end function getnumset_

!--------------------------------------------------------------------------
function getlastEnergy_(self) result(n)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! get lastEnergy from the MPfile_T class

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)    :: self
integer(kind=irg)                 :: n

n = self%MPDT%lastEnergy

end function getlastEnergy_

!--------------------------------------------------------------------------
function getnewPGnumber_(self) result(n)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! get newPGnumber from the MPfile_T class

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)    :: self
integer(kind=irg)                 :: n

n = self%MPDT%newPGnumber

end function getnewPGnumber_

!--------------------------------------------------------------------------
function getAveragedMP_(self) result(n)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! get AveragedMP from the MPfile_T class

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)    :: self
logical                           :: n

n = self%MPDT%AveragedMP

end function getAveragedMP_

!--------------------------------------------------------------------------
function getxtalname_(self) result(n)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! get xtalname from the MPfile_T class

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)    :: self
character(fnlen)                  :: n

n = trim(self%MPDT%xtalname)

end function getxtalname_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames, emnl)
!! author: MDG 
!! version: 1.0 
!! date: 02/05/20
!!
!! write namelist to HDF file; this routine can handle EBSD, ECP, and TKD name lists

use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_IO
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(MPfile_T), INTENT(INOUT)          :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames
class(SEMmasterNameListType), INTENT(IN):: emnl 

type(IO_T)                              :: Message
integer(kind=irg)                       :: hdferr, restart, uniform, combinesites, Esel, &
                                           useEnergyWeighting, dokinematical, n_int
integer(kind=irg), allocatable          :: io_int(:)
character(20), allocatable              :: intlist(:)
character(fnlen)                        :: dataset, groupname
character(fnlen,kind=c_char)            :: line2(1)
logical                                 :: g_exists, overwrite=.TRUE., isEBSD=.FALSE., &
                                           isECP=.FALSE., isTKD=.FALSE., isEBSDSHT=.FALSE.

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! next determine what kind of namelist we are dealing with (EBSD, ECP, or TKD)
select type(emnl) 
  type is (EBSDmasterNameListType)
    isEBSD = .TRUE.
    if (emnl%useEnergyWeighting) then 
      useEnergyWeighting = 1
    else 
      useEnergyWeighting = 0
    end if
    if (emnl%restart) then 
      restart = 1
    else 
      restart = 0
    end if
    Esel = emnl%Esel
    n_int = 8
    allocate( io_int(n_int), intlist(n_int) )
  type is (EBSDmasterSHTNameListType)
    isEBSDSHT = .TRUE.
    if (emnl%useEnergyWeighting) then 
      useEnergyWeighting = 1
    else 
      useEnergyWeighting = 0
    end if
    if (emnl%restart) then 
      restart = 1
    else 
      restart = 0
    end if
    Esel = emnl%Esel
    n_int = 8
    allocate( io_int(n_int), intlist(n_int) )
  type is (ECPmasterNameListType)
    isECP = .TRUE.
    n_int = 5
    allocate( io_int(n_int), intlist(n_int) )
  type is (TKDmasterNameListType)
    isTKD = .TRUE.
    if (emnl%restart) then 
      restart = 1
    else 
      restart = 0
    end if
    Esel = emnl%Esel
    n_int = 7
    allocate( io_int(n_int), intlist(n_int) )
  class default 
    call Message%printError('writeHDFNameList', 'unknown name list type requested')
end select

! convert all logicals to integer 1 or 0
if (emnl%kinematical) then 
  dokinematical = 1
else 
  dokinematical = 0
end if 
if (emnl%combinesites) then 
  combinesites = 1
else 
  combinesites = 0
end if

if (emnl%uniform) then 
  uniform = 1
else 
  uniform = 0
end if

! these are the common integer parameters for EBSD, ECP, and TKD
io_int(1) = emnl%npx  
io_int(2) = emnl%nthreads
io_int(3) = combinesites
io_int(4) = dokinematical
io_int(5) = uniform

intlist(1) = 'npx'
intlist(2) = 'nthreads'
intlist(3) = 'combinesites'
intlist(4) = 'kinematical'
intlist(5) = 'uniform'

! and here are the specific integer parameters
if ((isEBSD.eqv..TRUE.) .or. (isEBSDSHT.eqv..TRUE.)) then
  io_int(6) = restart
  io_int(7) = useEnergyWeighting
  io_int(8) = Esel
  intlist(6) = 'restart'
  intlist(7) = 'useEnergyWeighting'
  intlist(8) = 'Esel'
end if 

if (isTKD.eqv..TRUE.) then
  io_int(6) = restart
  io_int(7) = Esel
  intlist(6) = 'restart'
  intlist(7) = 'Esel'
end if 

! and write them to the HDF file
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write a single real
dataset = SC_dmin
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then 
  hdferr = HDF%writeDatasetFloat(dataset, emnl%dmin, overwrite)
else
  hdferr = HDF%writeDatasetFloat(dataset, emnl%dmin)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create dmin dataset',hdferr)

if (isEBSDSHT.eqv..TRUE.) then 
  dataset = SC_latgridtype
  line2(1) = 'Lambert'
  line2(1) = cstringify(line2(1))
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
  if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create latgridtype dataset',hdferr)
end if

dataset = SC_copyfromenergyfile
line2(1) = emnl%copyfromenergyfile
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then 
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create copyfromenergyfile dataset',hdferr)

dataset = SC_energyfile
line2(1) = emnl%energyfile
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then 
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create energyfile dataset',hdferr)

dataset = SC_BetheParametersFile
line2(1) = emnl%BetheParametersFile
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then 
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create BetheParametersFile dataset',hdferr)

! and pop this group off the stack
call HDF%pop()

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
recursive subroutine readMPfile_(self, HDF, HDFnames, mpnl, getkeVs, getmLPNH, getmLPSH, &
                                 getmasterSPNH, getmasterSPSH, keep4, defectMP)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! read a Master Pattern file into the correct namelist and data structure

use HDF5 
use mod_HDFsupport
use mod_HDFnames
use mod_io
use stringconstants
use ISO_C_BINDING 

IMPLICIT NONE

class(MPfile_T), INTENT(INOUT)                   :: self 
type(HDF_T), INTENT(INOUT)                       :: HDF
type(HDFnames_T), INTENT(INOUT)                  :: HDFnames
class(SEMmasterNameListType), INTENT(INOUT)      :: mpnl
logical,INTENT(IN),OPTIONAL                      :: getkeVs
logical,INTENT(IN),OPTIONAL                      :: getmLPNH
logical,INTENT(IN),OPTIONAL                      :: getmLPSH
logical,INTENT(IN),OPTIONAL                      :: getmasterSPNH
logical,INTENT(IN),OPTIONAL                      :: getmasterSPSH
logical,INTENT(IN),OPTIONAL                      :: keep4
logical,INTENT(IN),OPTIONAL                      :: defectMP

type(IO_T)                                       :: Message
character(fnlen)                                 :: infile, groupname, datagroupname, dataset
logical                                          :: stat, readonly, g_exists, f_exists, FL, keepall, dfMP
integer(kind=irg)                                :: ii, nlines, restart, combinesites, uniform, istat, hdferr, &
                                                    dokinematical, useEnergyWeighting, Esel
integer(kind=irg),allocatable                    :: iarray(:)
real(kind=sgl),allocatable                       :: farray(:)
real(kind=sgl),allocatable                       :: mLPNH(:,:,:,:), mLPNH3(:,:,:) 
integer(HSIZE_T)                                 :: dims(1), dims2(2), dims3(3), offset3(3), dims4(4) 
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)
logical                                          :: isEBSD=.FALSE., isECP=.FALSE., isTKD=.FALSE.
character(7)                                     :: modality 

dfMP = .FALSE.
if (present(defectMP)) then
  if (defectMP.eqv..TRUE.) dfMP = .TRUE.
end if

keepall = .FALSE.
if (present(keep4)) then
  if (keep4.eqv..TRUE.) keepall = .TRUE.
end if

! next determine what kind of namelist we are dealing with (EBSD, ECP, or TKD)
select type (mpnl)
  class is (EBSDmasterNameListType)
    isEBSD = .TRUE. 
  class is (ECPmasterNameListType)  ! this class has a different MP array size !!!
    isECP = .TRUE. 
  class is (TKDmasterNameListType)
    isTKD = .TRUE. 
  class default 
    call Message%printError('readMPfile', 'unknown master pattern type requested')
end select

associate( MPDT => self%MPDT )

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

hdferr = HDF%openGroup(HDFnames%get_ProgramData())  ! SC_MCOpenCL
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
    call HDF%pop(.TRUE.)
    call Message%printError('readMPfile','this is not a valid Master Pattern file')
end if
call HDF%pop()

!====================================
! check if this is an overlap EBSD pattern or a regular one  [ added by MDG, 06/19/19 ]; revised [10/18/19]
!====================================
dataset = 'READMEFIRST'
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
! if the dataset exists, then this is an overlap EBSD pattern file 
MPDT%AveragedMP = .FALSE.
MPDT%newPGnumber = -1
if (g_exists.eqv..TRUE.) then
    MPDT%AveragedMP = .TRUE.
    call Message%printMessage(' This master pattern file contains an averaged master pattern (generated by EMEBSDoverlap)')
end if

if (MPDT%AveragedMP.eqv..TRUE.) then
  ! read the new point group number from the newpgnum data set in the EBSDoverlapNameList NMLparameters group
  hdferr = HDF%openGroup(HDFnames%get_NMLparameters())
  groupname = SC_EBSDoverlapNameList
      hdferr = HDF%openGroup(groupname)

  dataset = 'newpgnum'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if(g_exists) call HDF%readDatasetInteger(dataset, hdferr, MPDT%newPGnumber)

  call HDF%pop()
  call HDF%pop()
end if

!====================================
! read all NMLparameters group datasets
!====================================
hdferr = HDF%openGroup(HDFnames%get_NMLparameters())
hdferr = HDF%openGroup(HDFnames%get_NMLlist())

! we need to set the newPGnumber parameter to the correct value, to reflect the fact that 
! the symmetry of the overlap pattern will be different [ added by MDG, 06/20/19 ]
if (MPDT%AveragedMP.eqv..TRUE.) then 
  dataset = 'newpgnumber'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if(g_exists) call HDF%readDatasetInteger(dataset, hdferr, MPDT%newPGnumber)
end if

! we need to first read all the parameters that are common to EBSD, ECP, and TKD namelists
dataset = SC_combinesites
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
mpnl%combinesites = .FALSE.
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetInteger(dataset, hdferr, combinesites)
    if (combinesites.ne.0) mpnl%combinesites = .TRUE.
end if

dataset = SC_copyfromenergyfile
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetStringArray(dataset, nlines, hdferr, stringarray)
    mpnl%copyfromenergyfile = trim(stringarray(1))
    deallocate(stringarray)
else
    mpnl%copyfromenergyfile = 'n'
end if 

dataset = SC_dmin
    call HDF%readDatasetFloat(dataset, hdferr, mpnl%dmin)

dataset = SC_energyfile
    call HDF%readDatasetStringArray(dataset, nlines, hdferr, stringarray)
    mpnl%energyfile = trim(stringarray(1))
    deallocate(stringarray)

dataset = SC_npx
    call HDF%readDatasetInteger(dataset, hdferr, mpnl%npx)

dataset = SC_nthreads
    call HDF%readDatasetInteger(dataset, hdferr, mpnl%nthreads)


dataset = SC_kinematical
call H5Lexists_f(HDF%getobjectID(), trim(dataset), g_exists, hdferr)
if(g_exists) then
    call HDF%readDatasetInteger(dataset, hdferr, dokinematical)
    mpnl%kinematical= .FALSE.
    if (dokinematical.ne.0) mpnl%kinematical = .TRUE.
end if

dataset = SC_uniform
mpnl%uniform = .FALSE.
call H5Lexists_f(HDF%getobjectID(), trim(dataset), g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetInteger(dataset, hdferr, uniform)
    if (uniform.ne.0) mpnl%uniform = .TRUE.
end if 

! then we read parameters that are specific to the stated modality
! we need to do this inside a select type construct, since we need to access 
! members of inherited types...
select type (mpnl)
  class is (EBSDmasterNameListType)
    dataset = SC_Esel
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if(g_exists) call HDF%readDatasetInteger(dataset, hdferr, Esel)
    mpnl%Esel = Esel

    dataset = SC_restart
    call H5Lexists_f(HDF%getobjectID(), trim(dataset), g_exists, hdferr)
    if(g_exists) then
        call HDF%readDatasetInteger(dataset, hdferr, restart)
        mpnl%restart = .FALSE.
        if (restart.ne.0) mpnl%restart = .TRUE.
    end if

    dataset = SC_useEnergyWeighting
    call H5Lexists_f(HDF%getobjectID(), trim(dataset), g_exists, hdferr)
    if(g_exists) then
        call HDF%readDatasetInteger(dataset, hdferr, useEnergyWeighting)
        mpnl%useEnergyWeighting = .FALSE.
        if (useEnergyWeighting.ne.0) mpnl%useEnergyWeighting = .TRUE.
    end if
  class is (TKDmasterNameListType)
    dataset = SC_Esel
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if(g_exists) call HDF%readDatasetInteger(dataset, hdferr, mpnl%Esel)

    dataset = SC_restart
    call H5Lexists_f(HDF%getobjectID(), trim(dataset), g_exists, hdferr)
    if(g_exists) then
        call HDF%readDatasetInteger(dataset, hdferr, restart)
        mpnl%restart = .FALSE.
        if (restart.ne.0) mpnl%restart = .TRUE.
    end if
end select

! and close the NMLparameters group
call HDF%pop()
call HDF%pop()

!====================================
!====================================

! open the Master Pattern data group
hdferr = HDF%openGroup(HDFnames%get_EMData())
hdferr = HDF%openGroup(HDFnames%get_ProgramData())

! integers
dataset = SC_lastEnergy
call H5Lexists_f(HDF%getobjectID(), trim(dataset), g_exists, hdferr)
if(g_exists) then
    call HDF%readDatasetInteger(dataset, hdferr, MPDT%lastEnergy)
end if

MPDT%numEbins = 1
dataset = SC_numEbins
call H5Lexists_f(HDF%getobjectID(), trim(dataset), g_exists, hdferr)
if(g_exists) then
    call HDF%readDatasetInteger(dataset, hdferr, MPDT%numEbins)
end if 

dataset = SC_numset
    call HDF%readDatasetInteger(dataset, hdferr, MPDT%numset)

dataset = SC_xtalname
    call HDF%readDatasetStringArray(dataset, nlines, hdferr, stringarray)
    MPDT%xtalname = trim(stringarray(1))
    deallocate(stringarray)

! dataset = SC_BetheParameters
!     call HDF%readDatasetFloatArray(dataset, dims, hdferr, MPDT%BetheParameters)

! various optional arrays
if (present(getkeVs)) then 
  if (getkeVs.eqv..TRUE.) then
    dataset = SC_keVs
    call HDF%readDatasetFloatArray(dataset, dims, hdferr, MPDT%keVs)
  end if 
end if


if (present(getmLPNH)) then 
  if (getmLPNH.eqv..TRUE.) then
    dataset = SC_mLPNH
    if ((isEBSD.eqv..TRUE.).or.(isTKD.eqv..TRUE.)) then 
      if (dfMP.eqv..TRUE.) then 
        call HDF%readDatasetFloatArray(dataset, dims3, hdferr, mLPNH3)
        allocate(MPDT%mLPNH(-mpnl%npx:mpnl%npx,-mpnl%npx:mpnl%npx,dims3(3)),stat=istat)
        MPDT%mLPNH = mLPNH3
        deallocate(mLPNH3)
      else
        call HDF%readDatasetFloatArray(dataset, dims4, hdferr, mLPNH)
        if (keepall) then
          allocate(MPDT%mLPNH4(-mpnl%npx:mpnl%npx,-mpnl%npx:mpnl%npx,MPDT%numEbins, dims4(4)),stat=istat)
          MPDT%mLPNH4 = mLPNH
        else
          allocate(MPDT%mLPNH(-mpnl%npx:mpnl%npx,-mpnl%npx:mpnl%npx,MPDT%numEbins),stat=istat)
          MPDT%mLPNH(-mpnl%npx:mpnl%npx,-mpnl%npx:mpnl%npx,1:MPDT%numEbins) = sum(mLPNH,4)
        end if
        deallocate(mLPNH)
      end if
    end if 
    if (isECP.eqv..TRUE.) then 
        call HDF%readDatasetFloatArray(dataset, dims3, hdferr, mLPNH3)
        if (keepall) then 
          allocate(MPDT%mLPNH4(-mpnl%npx:mpnl%npx,-mpnl%npx:mpnl%npx,dims3(3),1),stat=istat)
          MPDT%mLPNH4(:,:,:,1) = mLPNH3(:,:,:)
        else
          allocate(MPDT%mLPNH(-mpnl%npx:mpnl%npx,-mpnl%npx:mpnl%npx,dims3(3)),stat=istat)
          MPDT%mLPNH = mLPNH3
        end if
        deallocate(mLPNH3)
    end if 
  end if 
end if

if (present(getmLPSH)) then 
  if (getmLPSH.eqv..TRUE.) then
    dataset = SC_mLPSH
    if ((isEBSD.eqv..TRUE.).or.(isTKD.eqv..TRUE.)) then 
      if (dfMP.eqv..TRUE.) then 
        call HDF%readDatasetFloatArray(dataset, dims3, hdferr, mLPNH3)
        allocate(MPDT%mLPSH(-mpnl%npx:mpnl%npx,-mpnl%npx:mpnl%npx,dims3(3)),stat=istat)
        MPDT%mLPSH = mLPNH3
        deallocate(mLPNH3)
      else
        call HDF%readDatasetFloatArray(dataset, dims4, hdferr, mLPNH)
        if (keepall) then
          allocate(MPDT%mLPSH4(-mpnl%npx:mpnl%npx,-mpnl%npx:mpnl%npx,MPDT%numEbins, dims4(4)),stat=istat)
          MPDT%mLPSH4 = mLPNH
        else
          allocate(MPDT%mLPSH(-mpnl%npx:mpnl%npx,-mpnl%npx:mpnl%npx,MPDT%numEbins),stat=istat)
          MPDT%mLPSH(-mpnl%npx:mpnl%npx,-mpnl%npx:mpnl%npx,1:MPDT%numEbins) = sum(mLPNH,4)
        end if
        deallocate(mLPNH)
      end if
    end if 
    if (isECP.eqv..TRUE.) then 
        call HDF%readDatasetFloatArray(dataset, dims3, hdferr, mLPNH3)
        if (keepall) then 
          allocate(MPDT%mLPSH4(-mpnl%npx:mpnl%npx,-mpnl%npx:mpnl%npx,dims3(3),1),stat=istat)
          MPDT%mLPSH4(:,:,:,1) = mLPNH3(:,:,:)
        else
          allocate(MPDT%mLPSH(-mpnl%npx:mpnl%npx,-mpnl%npx:mpnl%npx,dims3(3)),stat=istat)
          MPDT%mLPSH = mLPNH3
        end if
        deallocate(mLPNH3)
    end if 
  end if 
end if

if (present(getmasterSPNH)) then 
  if (getmasterSPNH.eqv..TRUE.) then
    dataset = SC_masterSPNH
    call HDF%readDatasetFloatArray(dataset, dims3, hdferr, MPDT%masterSPNH)
  end if 
end if

if (present(getmasterSPSH)) then 
  if (getmasterSPSH.eqv..TRUE.) then
    dataset = SC_masterSPSH
    call HDF%readDatasetFloatArray(dataset, dims3, hdferr, MPDT%masterSPSH)
  end if 
end if

! and close the HDF5 Master Pattern file
call HDF%pop(.TRUE.)

call Message%printMessage(' --> Completed reading master pattern data from '//trim(self%MPfile), frm = "(A/)")

end associate 

end subroutine readMPfile_


!--------------------------------------------------------------------------
recursive subroutine copyMPdata_(self, EMsoft, HDF, HDFnames, inputfile, outputfile, h5, skipCrystalData)
  !! author: MDG 
  !! version: 1.0 
  !! date: 03/20/20
  !!
  !! copy Master Pattern data from one file to a new file using h5copy

use mod_EMsoft
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_io
use stringconstants

IMPLICIT NONE

class(MPfile_T),INTENT(INOUT)     :: self
type(EMsoft_T),INTENT(INOUT)      :: EMsoft
type(HDF_T),INTENT(INOUT)         :: HDF 
type(HDFnames_T),INTENT(INOUT)    :: HDFnames
character(fnlen),INTENT(IN)       :: inputfile
character(fnlen),INTENT(IN)       :: outputfile
character(fnlen),INTENT(IN)       :: h5
logical,INTENT(IN),OPTIONAL       :: skipCrystalData

type(IO_T)                        :: Message
character(fnlen)                  :: infile, outfile, h5copypath, groupname
character(512)                    :: cmd, cmd2
logical                           :: f_exists, readonly, developer
integer(kind=irg)                 :: hdferr
character(fnlen)                  :: dev

! first we make sure that we actually have the h5copy program available
! check for EMDevelop parameter 
developer = .FALSE.
dev = EMsoft%getConfigParameter('Develop')
if (trim(dev).eq.'Yes') developer = .TRUE.

if (developer.eqv..TRUE.) then 
! if TRUE, use EMsoft_geth5copypath which is defined at configure time 
  h5copypath = trim(EMsoft%getConfigParameter('h5copypath'))//' -p -v '
  h5copypath = EMsoft%toNativePath(h5copypath)
else 
! if FALSE, check name list h5copypath parameter 
  if (trim(h5).ne.'undefined') then 
    h5copypath = trim(h5)//' -p -v '
    h5copypath = EMsoft%toNativePath(h5copypath)
  else 
! if undefined, then fail
    call Message%printError('copyMPdata','h5copypath must be set in the name list file ')
  end if
end if

call Message%printMessage(' Using '//trim(h5copypath)//' to copy Monte Carlo data to new file')

! first make sure that the input file exists and has MC data in it
infile = trim(EMsoft%generateFilePath('EMdatapathname',inputfile))
inquire(file=trim(infile), exist=f_exists)

outfile = trim(EMsoft%generateFilePath('EMdatapathname',outputfile))

! if the file does not exist, abort the program with an error message
if (f_exists.eqv..FALSE.) then 
  call Message%printError('copyMPdata','Master PAttern copyfromenergyfile does not exist: '//trim(infile))
end if

! make sure it has MCopenCL data in it; hdf open is done in the calling program
readonly = .TRUE.
hdferr =  HDF%openFile(infile, readonly)

hdferr = HDF%openGroup(HDFnames%get_EMData())
if (hdferr.eq.-1) then 
  call Message%printError('copyMPdata','EMData group does not exist in '//trim(infile))
end if

hdferr = HDF%openGroup(HDFnames%get_ProgramData())
if (hdferr.eq.-1) then 
  call Message%printError('copyMPdata','master group does not exist in '//trim(infile))
end if

call HDF%pop(.TRUE.)

! OK, if we get here, then the file does exist and it contains Master Pattern data, 
! so we let the user know
call Message%printMessage('--> Input file contains Master Pattern data')

! next, we copy the necessary groups into the new Monte Carlo file
cmd = trim(h5copypath)//' -i "'//trim(infile)
cmd = trim(cmd)//'" -o "'//trim(outfile)

if (.not.present(skipCrystalData)) then 
  cmd2 = trim(cmd)//'" -s "/CrystalData" -d "/CrystalData"'
  call system(trim(cmd2))
end if 

cmd2 = trim(cmd)//'" -s "/EMData" -d "/EMData"'
call system(trim(cmd2))

cmd2 = trim(cmd)//'" -s "/EMheader" -d "/EMheader"'
call system(trim(cmd2))

cmd2 = trim(cmd)//'" -s "/NMLfiles" -d "/NMLfiles"'
call system(trim(cmd2))

cmd2 = trim(cmd)//'" -s "/NMLparameters" -d "/NMLparameters"'
call system(trim(cmd2))

call Message%printMessage('--> Output file generated with Master Pattern data copied from '//trim(infile))

end subroutine copyMPdata_


!--------------------------------------------------------------------------
recursive subroutine copyMPoverlapdata_(self, EMsoft, HDF, HDFnames, inputfile, outputfile, h5, skipCrystalData)
  !! author: MDG 
  !! version: 1.0 
  !! date: 03/23/20
  !!
  !! copy Master Pattern data from one file to a new MP overlap file using h5copy

use mod_EMsoft
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_io
use stringconstants

IMPLICIT NONE

class(MPfile_T),INTENT(INOUT)     :: self
type(EMsoft_T),INTENT(INOUT)      :: EMsoft
type(HDF_T),INTENT(INOUT)         :: HDF 
type(HDFnames_T),INTENT(INOUT)    :: HDFnames
character(fnlen),INTENT(IN)       :: inputfile
character(fnlen),INTENT(IN)       :: outputfile
character(fnlen),INTENT(IN)       :: h5
logical,INTENT(IN),OPTIONAL       :: skipCrystalData

type(IO_T)                        :: Message
character(fnlen)                  :: infile, outfile, h5copypath, groupname
character(512)                    :: cmd, cmd2
logical                           :: f_exists, readonly, developer
integer(kind=irg)                 :: hdferr
character(fnlen)                  :: dev

! first we make sure that we actually have the h5copy program available
! check for EMDevelop parameter 
developer = .FALSE.
dev = EMsoft%getConfigParameter('Develop')
if (trim(dev).eq.'Yes') developer = .TRUE.

if (developer.eqv..TRUE.) then 
! if TRUE, use EMsoft_geth5copypath which is defined at configure time 
  h5copypath = trim(EMsoft%getConfigParameter('h5copypath'))//' -p -v '
  h5copypath = EMsoft%toNativePath(h5copypath)
else 
! if FALSE, check name list h5copypath parameter 
  if (trim(h5).ne.'undefined') then 
    h5copypath = trim(h5)//' -p -v '
    h5copypath = EMsoft%toNativePath(h5copypath)
  else 
! if undefined, then fail
    call Message%printError('copyMPdata','h5copypath must be set in the name list file ')
  end if
end if

call Message%printMessage(' Using '//trim(h5copypath)//' to copy Monte Carlo data to new file')

! first make sure that the input file exists and has MC data in it
infile = trim(EMsoft%generateFilePath('EMdatapathname',inputfile))
inquire(file=trim(infile), exist=f_exists)

outfile = trim(EMsoft%generateFilePath('EMdatapathname',outputfile))

! if the file does not exist, abort the program with an error message
if (f_exists.eqv..FALSE.) then 
  call Message%printError('copyMPdata','Master Pattern copyfromenergyfile does not exist: '//trim(infile))
end if

! make sure it has MCopenCL data in it; hdf open is done in the calling program
readonly = .TRUE.
hdferr =  HDF%openFile(infile, readonly)

hdferr = HDF%openGroup(HDFnames%get_EMData())
if (hdferr.eq.-1) then 
  call Message%printError('copyMPdata','EMData group does not exist in '//trim(infile))
end if

hdferr = HDF%openGroup(HDFnames%get_ProgramData())
if (hdferr.eq.-1) then 
  call Message%printError('copyMPdata','master group does not exist in '//trim(infile))
end if

call HDF%pop(.TRUE.)

! OK, if we get here, then the file does exist and it contains Master Pattern data, 
! so we let the user know
call Message%printMessage('--> Input file contains Master Pattern data')

! next, we copy the necessary groups into the new Monte Carlo file
cmd = trim(h5copypath)//' -i "'//trim(infile)
cmd = trim(cmd)//'" -o "'//trim(outfile)

if (.not.present(skipCrystalData)) then 
  cmd2 = trim(cmd)//'" -s "/CrystalData" -d "/CrystalData"'
  call system(trim(cmd2))
end if 

cmd2 = trim(cmd)//'" -s "/EMData/MCOpenCL" -d "/EMData/MCOpenCL"'
call system(trim(cmd2))

cmd2 = trim(cmd)//'" -s "/EMheader/MCOpenCL" -d "/EMheader/MCOpenCL"'
call system(trim(cmd2))

cmd2 = trim(cmd)//'" -s "/NMLfiles/MCOpenCLNML" -d "/NMLfiles/MCOpenCLNML"'
call system(trim(cmd2))

cmd2 = trim(cmd)//'" -s "/NMLparameters/MCCLNameList" -d "/NMLparameters/MCCLNameList"'
call system(trim(cmd2))

call Message%printMessage('--> Output file generated with Master Pattern data copied from '//trim(infile))

end subroutine copyMPoverlapdata_





end module mod_MPfiles