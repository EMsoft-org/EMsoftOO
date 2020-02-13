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
    type(MPdataType)            :: MPDT
    type(EBSDmasterNameListType):: nml 
    character(fnlen)            :: MPfile

  contains
  private 

    ! procedure, pass(self) :: readMPfile_
    ! procedure, pass(self) :: writeMPfile_
    procedure, pass(self) :: writeHDFNameList_
    ! procedure, pass(self) :: copynml_
    ! procedure, pass(self) :: getnml_
    ! procedure, pass(self) :: setFileName_

    ! generic, public :: readMPfile => readMPfile_
    generic, public :: writeHDFNameList => writeHDFNameList_
    ! generic, public :: writeMPfile => writeMPfile_
    ! generic, public :: copynml => copynml_
    ! generic, public :: getnml => getnml_
    ! generic, public :: setFileName => setFileName_

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

class(MPfile_T), INTENT(INOUT)      :: self 
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




end module mod_MPfiles