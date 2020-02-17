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
  !! class for Master Pattern file handling 

use mod_kinds
use mod_global
use stringconstants

IMPLICIT NONE 
private

! namelist for the EMEBSDmaster program
type, public :: EBSDmasterNameListType
  integer(kind=irg) :: stdout
  integer(kind=irg) :: npx
  integer(kind=irg) :: Esel
  integer(kind=irg) :: nthreads
  real(kind=sgl)    :: dmin
  character(3)      :: Notify
  character(fnlen)  :: copyfromenergyfile
  character(fnlen)  :: h5copypath
  character(fnlen)  :: energyfile
  character(fnlen)  :: outname
  character(fnlen)  :: BetheParametersFile
  logical           :: useEnergyWeighting
  logical           :: combinesites
  logical           :: restart
  logical           :: uniform
end type EBSDmasterNameListType


type, public :: MPdataType
        integer(kind=irg)               :: lastEnergy
        integer(kind=irg)               :: numEbins
        integer(kind=irg)               :: numset
        integer(kind=irg)               :: newPGnumber
        logical                         :: AveragedMP
        character(fnlen)                :: xtalname
        real(kind=sgl),allocatable      :: BetheParameters(:)
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
    type(MPdataType),public             :: MPDT
    type(EBSDmasterNameListType),public :: nml 
    character(fnlen)                    :: MPfile

  contains
  private 

    procedure, pass(self) :: readMPfile_
    procedure, pass(self) :: setFileName_
    ! procedure, pass(self) :: writeMPfile_
    procedure, pass(self) :: writeHDFNameList_
    procedure, pass(self) :: copynml_
    procedure, pass(self) :: getnml_
    procedure, pass(self) :: getlastEnergy_
    procedure, pass(self) :: getnumEbins_
    procedure, pass(self) :: getnumset_
    procedure, pass(self) :: getnewPGnumber_
    procedure, pass(self) :: getAveragedMP_
    procedure, pass(self) :: getxtalname_
    procedure, pass(self) :: copyBetheParameters_
    procedure, pass(self) :: copykeVs_
    procedure, pass(self) :: copymLPNH4_
    procedure, pass(self) :: copymLPSH4_
    procedure, pass(self) :: copymLPNH_
    procedure, pass(self) :: copymLPSH_
    procedure, pass(self) :: copymasterSPNH_
    procedure, pass(self) :: copymasterSPSH_
 

    generic, public :: readMPfile => readMPfile_
    generic, public :: setFileName => setFileName_
    ! generic, public :: writeMPfile => writeMPfile_
    generic, public :: writeHDFNameList => writeHDFNameList_
    generic, public :: copynml => copynml_
    generic, public :: getnml => getnml_
    generic, public :: getlastEnergy => getlastEnergy_
    generic, public :: getnumEbins => getnumEbins_
    generic, public :: getnumset => getnumset_
    generic, public :: getnewPGnumber => getnewPGnumber_
    generic, public :: getAveragedMP => getAveragedMP_
    generic, public :: getxtalname => getxtalname_
    generic, public :: copyBetheParameters => copyBetheParameters_
    generic, public :: copykeVs => copykeVs_
    generic, public :: copymLPNH4 => copymLPNH4_
    generic, public :: copymLPSH4 => copymLPSH4_
    generic, public :: copymLPNH => copymLPNH_
    generic, public :: copymLPSH => copymLPSH_
    generic, public :: copymasterSPNH => copymasterSPNH_
    generic, public :: copymasterSPSH => copymasterSPSH_

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

!--------------------------------------------------------------------------
subroutine copynml_(self, nml)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! copy the namelist into the MPfile_T class for writing to file

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)            :: self
type(EBSDmasterNameListType),INTENT(IN)   :: nml

self%nml = nml

end subroutine copynml_

!--------------------------------------------------------------------------
subroutine copyBetheParameters_(self, acc, keep)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! copy the BetheParameters array

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)                :: self
real(kind=sgl), allocatable, INTENT(OUT)      :: acc(:)
logical, INTENT(IN), OPTIONAL                 :: keep 

integer(kind=irg)                             :: s(1)

s = shape(self%MPDT%BetheParameters)
allocate(acc(s(1)))

acc = self%MPDT%BetheParameters 

if (.not.present(keep)) deallocate(self%MPDT%BetheParameters)

end subroutine copyBetheParameters_

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
function getnml_(self) result(nml)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! get the namelist from the MPfile_T class

IMPLICIT NONE 

class(MPfile_T), INTENT(INOUT)            :: self
type(EBSDmasterNameListType)              :: nml

nml = self%nml

end function getnml_

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
recursive subroutine writeHDFNameList_(self, HDF)
!! author: MDG 
!! version: 1.0 
!! date: 02/05/20
!!
!! write namelist to HDF file

use HDF5
use mod_HDFsupport
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(MPfile_T), INTENT(INOUT)          :: self 
type(HDF_T), INTENT(INOUT)              :: HDF

integer(kind=irg),parameter             :: n_int = 8, n_real = 1
integer(kind=irg)                       :: hdferr, io_int(n_int), restart, uniform, combinesites, &
                                           useEnergyWeighting
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, groupname
character(fnlen,kind=c_char)            :: line2(1)
logical                                 :: g_exists, overwrite=.TRUE.

associate( emnl => self%nml )

! create the group for this namelist
groupname = SC_EBSDMasterNameList
hdferr = HDF%createGroup(groupname)

! write all the single integers
if (emnl%combinesites) then 
  combinesites = 1
else 
  combinesites = 0
end if
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
if (emnl%uniform) then 
  uniform = 1
else 
  uniform = 0
end if
io_int = (/ emnl%stdout, emnl%npx, emnl%Esel, emnl%nthreads, combinesites, restart, uniform, useEnergyWeighting /)
intlist(1) = 'stdout'
intlist(2) = 'npx'
intlist(3) = 'Esel'
intlist(4) = 'nthreads'
intlist(5) = 'combinesites'
intlist(6) = 'restart'
intlist(7) = 'uniform'
intlist(8) = 'useEnergyWeighting'
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

dataset = SC_latgridtype
line2(1) = 'Lambert'
line2(1) = cstringify(line2(1))
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create latgridtype dataset',hdferr)

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

dataset = 'BetheParametersFile'
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

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
recursive subroutine readMPfile_(self, HDF, getkeVs, getmLPNH, getmLPSH, &
                                 getmasterSPNH, getmasterSPSH, keep4, defectMP)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! read an EBSD Master Pattern file into the correct namelist and data structure

use HDF5 
use mod_HDFsupport
use mod_io
use stringconstants
use ISO_C_BINDING 

IMPLICIT NONE

class(MPfile_T), INTENT(INOUT)                   :: self 
type(HDF_T), INTENT(INOUT)                       :: HDF
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
integer(kind=irg)                                :: ii, nlines, restart, combinesites, uniform, istat, hdferr
integer(kind=irg),allocatable                    :: iarray(:)
real(kind=sgl),allocatable                       :: farray(:)
real(kind=sgl),allocatable                       :: mLPNH(:,:,:,:), mLPNH3(:,:,:)
integer(HSIZE_T)                                 :: dims(1), dims2(2), dims3(3), offset3(3), dims4(4) 
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)

dfMP = .FALSE.
if (present(defectMP)) then
  if (defectMP.eqv..TRUE.) dfMP = .TRUE.
end if

keepall = .FALSE.
if (present(keep4)) then
  if (keep4.eqv..TRUE.) keepall = .TRUE.
end if

associate( mpnl => self%nml,  MPDT => self%MPDT )

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
if (dfMP.eqv..TRUE.) then 
  datagroupname = '/EMheader/EBSDdefectmaster'
else
  datagroupname = '/EMheader/EBSDmaster'
end if
call H5Lexists_f(HDF%getobjectID(),trim(datagroupname),g_exists, hdferr)
if (.not.g_exists) then
  call Message%printError('readMPfile','This HDF file does not contain Master Pattern header data')
end if

groupname = SC_EMheader
hdferr = HDF%openGroup(groupname)
groupname = SC_MCOpenCL
hdferr = HDF%openGroup(groupname)
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
groupname = SC_NMLfiles
    hdferr = HDF%openGroup(groupname)
dataset = 'EBSDmasterNML'
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists.eqv..FALSE.) then
    call HDF%pop(.TRUE.)
    call Message%printError('readMPfile','this is not an EBSD Master Pattern file')
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
  groupname = SC_NMLparameters
      hdferr = HDF%openGroup(groupname)
  groupname = 'EBSDoverlapNameList'
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
groupname = SC_NMLparameters
    hdferr = HDF%openGroup(groupname)
groupname = SC_EBSDMasterNameList
    hdferr = HDF%openGroup(groupname)

! we'll read these roughly in the order that the HDFView program displays them...
dataset = SC_Esel
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if(g_exists) call HDF%readDatasetInteger(dataset, hdferr, mpnl%Esel)

! we need to set the newPGnumber parameter to the correct value, to reflect the fact that 
! the symmetry of the overlap pattern will be different [ added by MDG, 06/20/19 ]
if (MPDT%AveragedMP.eqv..TRUE.) then 
  dataset = 'newpgnumber'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if(g_exists) call HDF%readDatasetInteger(dataset, hdferr, MPDT%newPGnumber)
end if

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

dataset = SC_restart
call H5Lexists_f(HDF%getobjectID(), trim(dataset), g_exists, hdferr)
if(g_exists) then
    call HDF%readDatasetInteger(dataset, hdferr, restart)
    mpnl%restart = .FALSE.
    if (restart.ne.0) mpnl%restart = .TRUE.
end if

dataset = SC_stdout
call H5Lexists_f(HDF%getobjectID(), trim(dataset), g_exists, hdferr)
if(g_exists) then
    call HDF%readDatasetInteger(dataset, hdferr, mpnl%stdout)
end if

dataset = SC_uniform
mpnl%uniform = .FALSE.
call H5Lexists_f(HDF%getobjectID(), trim(dataset), g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetInteger(dataset, hdferr, uniform)
    if (uniform.ne.0) mpnl%uniform = .TRUE.
end if 

! and close the NMLparameters group
    call HDF%pop()
    call HDF%pop()
!====================================
!====================================

! open the Monte Carlo data group
groupname = SC_EMData
    hdferr = HDF%openGroup(groupname)
if (dfMP.eqv..TRUE.) then 
  groupname = SC_EBSDdefectmaster
else
  groupname = SC_EBSDmaster
end if
    hdferr = HDF%openGroup(groupname)

! integers
dataset = SC_lastEnergy
call H5Lexists_f(HDF%getobjectID(), trim(dataset), g_exists, hdferr)
if(g_exists) then
    call HDF%readDatasetInteger(dataset, hdferr, MPDT%lastEnergy)
end if

dataset = SC_numEbins
    call HDF%readDatasetInteger(dataset, hdferr, MPDT%numEbins)

dataset = SC_numset
    call HDF%readDatasetInteger(dataset, hdferr, MPDT%numset)

dataset = SC_xtalname
    call HDF%readDatasetStringArray(dataset, nlines, hdferr, stringarray)
    MPDT%xtalname = trim(stringarray(1))
    deallocate(stringarray)

dataset = SC_BetheParameters
    call HDF%readDatasetFloatArray(dataset, dims, hdferr, MPDT%BetheParameters)

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
        MPDT%mLPNH = sum(mLPNH,4)
      end if
      deallocate(mLPNH)
    end if
  end if 
end if

if (present(getmLPSH)) then 
  if (getmLPSH.eqv..TRUE.) then
    dataset = SC_mLPSH
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
        MPDT%mLPSH = sum(mLPNH,4)
      end if
      deallocate(mLPNH)
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

call Message%printMessage(' -> completed reading master pattern data from '//trim(infile), frm = "(A/)")

end associate 

end subroutine readMPfile_




end module mod_MPfiles