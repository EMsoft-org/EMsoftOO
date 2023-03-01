! ###################################################################
! Copyright (c) 2013-2023, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_ISEDI
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/27/23
  !!
  !! class definition for the EMISEDI program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMISEDI program
type, public :: ISEDINameListType
  integer(kind=irg)         :: ipf_wd
  integer(kind=irg)         :: ipf_ht
  integer(kind=irg)         :: nnk
  integer(kind=irg)         :: nosm
  integer(kind=irg)         :: ncubochoric
  integer(kind=irg)         :: nsteps
  integer(kind=irg)         :: nbatch
  integer(kind=irg)         :: nthreads
  real(kind=sgl)            :: stepX
  real(kind=sgl)            :: stepY
  real(kind=sgl)            :: omega
  real(kind=sgl)            :: omega_step
  real(kind=sgl)            :: tiltaxis(3)
  character(fnlen)          :: HDFstrings(10)
  character(fnlen)          :: exptfile
  character(fnlen)          :: datafile
  character(fnlen)          :: ctffile
  character(fnlen)          :: angfile
  character(fnlen)          :: masterfile
end type ISEDINameListType

! class definition
type, public :: ISEDI_T
private 
  character(fnlen)       :: nmldeffile = 'EMISEDI.nml'
  type(ISEDINameListType)  :: nml 

contains
private 
  procedure, pass(self) :: setipf_wd_
  procedure, pass(self) :: getipf_wd_
  procedure, pass(self) :: setipf_ht_
  procedure, pass(self) :: getipf_ht_
  procedure, pass(self) :: setnnk_
  procedure, pass(self) :: getnnk_
  procedure, pass(self) :: setnosm_
  procedure, pass(self) :: getnosm_
  procedure, pass(self) :: setncubochoric_
  procedure, pass(self) :: getncubochoric_
  procedure, pass(self) :: setnsteps_
  procedure, pass(self) :: getnsteps_
  procedure, pass(self) :: setnbatch_
  procedure, pass(self) :: getnbatch_
  procedure, pass(self) :: setnthreads_
  procedure, pass(self) :: getnthreads_
  procedure, pass(self) :: setstepX_
  procedure, pass(self) :: getstepX_
  procedure, pass(self) :: setstepY_
  procedure, pass(self) :: getstepY_
  procedure, pass(self) :: settiltaxis_
  procedure, pass(self) :: gettiltaxis_
  procedure, pass(self) :: setomega_
  procedure, pass(self) :: getomega_
  procedure, pass(self) :: setomega_step_
  procedure, pass(self) :: getomega_step_
  procedure, pass(self) :: setHDFstrings_
  procedure, pass(self) :: getHDFstrings_
  procedure, pass(self) :: setexptfile_
  procedure, pass(self) :: getexptfile_
  procedure, pass(self) :: setdatafile_
  procedure, pass(self) :: getdatafile_
  procedure, pass(self) :: setctffile_
  procedure, pass(self) :: getctffile_
  procedure, pass(self) :: setangfile_
  procedure, pass(self) :: getangfile_
  procedure, pass(self) :: setmasterfile_
  procedure, pass(self) :: getmasterfile_
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: ISEDI_

  generic, public :: setipf_wd => setipf_wd_
  generic, public :: getipf_wd => getipf_wd_
  generic, public :: setipf_ht => setipf_ht_
  generic, public :: getipf_ht => getipf_ht_
  generic, public :: setnnk => setnnk_
  generic, public :: getnnk => getnnk_
  generic, public :: setnosm => setnosm_
  generic, public :: getnosm => getnosm_
  generic, public :: setncubochoric => setncubochoric_
  generic, public :: getncubochoric => getncubochoric_
  generic, public :: setnsteps => setnsteps_
  generic, public :: getnsteps => getnsteps_
  generic, public :: setnbatch => setnbatch_
  generic, public :: getnbatch => getnbatch_
  generic, public :: setnthreads => setnthreads_
  generic, public :: getnthreads => getnthreads_
  generic, public :: setstepX => setstepX_
  generic, public :: getstepX => getstepX_
  generic, public :: setstepY => setstepY_
  generic, public :: getstepY => getstepY_
  generic, public :: settiltaxis => settiltaxis_
  generic, public :: gettiltaxis => gettiltaxis_
  generic, public :: setomega => setomega_
  generic, public :: getomega => getomega_
  generic, public :: setomega_step => setomega_step_
  generic, public :: getomega_step => getomega_step_
  generic, public :: setHDFstrings => setHDFstrings_
  generic, public :: getHDFstrings => getHDFstrings_
  generic, public :: setexptfile => setexptfile_
  generic, public :: getexptfile => getexptfile_
  generic, public :: setdatafile => setdatafile_
  generic, public :: getdatafile => getdatafile_
  generic, public :: setctffile => setctffile_
  generic, public :: getctffile => getctffile_
  generic, public :: setangfile => setangfile_
  generic, public :: getangfile => getangfile_
  generic, public :: setmasterfile => setmasterfile_
  generic, public :: getmasterfile => getmasterfile_
  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: ISEDI => ISEDI_

end type ISEDI_T

! the constructor routine for this class 
interface ISEDI_T
  module procedure ISEDI_constructor
end interface ISEDI_T

contains

!--------------------------------------------------------------------------
type(ISEDI_T) function ISEDI_constructor( nmlfile ) result(ISEDI)
!! author: MDG 
!! version: 1.0 
!! date: 02/27/23
!!
!! constructor for the ISEDI_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call ISEDI%readNameList(nmlfile)

end function ISEDI_constructor

!--------------------------------------------------------------------------
subroutine ISEDI_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 02/27/23
!!
!! destructor for the ISEDI_T Class
 
IMPLICIT NONE

type(ISEDI_T), INTENT(INOUT)  :: self 

call reportDestructor('ISEDI_T')

end subroutine ISEDI_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/27/23
!!
!! read the namelist from an nml file for the ISEDI_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(ISEDI_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg)                    :: ipf_wd
integer(kind=irg)                    :: ipf_ht
integer(kind=irg)                    :: nnk
integer(kind=irg)                    :: nosm
integer(kind=irg)                    :: ncubochoric
integer(kind=irg)                    :: nsteps
integer(kind=irg)                    :: nbatch
integer(kind=irg)                    :: nthreads
real(kind=sgl)                       :: stepX
real(kind=sgl)                       :: stepY
real(kind=sgl)                       :: omega
real(kind=sgl)                       :: omega_step
real(kind=sgl)                       :: tiltaxis(3)
character(fnlen)                     :: HDFstrings(10)
character(fnlen)                     :: exptfile
character(fnlen)                     :: datafile
character(fnlen)                     :: ctffile
character(fnlen)                     :: angfile
character(fnlen)                     :: masterfile

namelist / ISEDIdata / ipf_wd, ipf_ht, nnk, nosm, ncubochoric, nsteps, nbatch, nthreads, &
                       stepX, stepY, omega, omega_step, HDFstrings, exptfile, datafile, &
                       ctffile, angfile, masterfile, tiltaxis

! set the input parameters to default values
ipf_wd = 100
ipf_ht = 100
stepX = 1.0
stepY = 1.0
nnk = 20
nosm = 20
ncubochoric = 100
omega = 0.0
omega_step = 1.0
tiltaxis = (/ 0.0, 1.0, 0.0 /)
nsteps = 4
HDFstrings = (/ '', '', '', '', '', '', '', '', '', '' /)
exptfile = 'undefined'
datafile = 'undefined'
ctffile = 'undefined'
angfile = 'undefined'
masterfile = 'undefined'
nbatch = 1024
nthreads = 1

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=ISEDIdata)
    close(UNIT=dataunit,STATUS='keep')

! check for required entries
    if (trim(exptfile).eq.'undefined') then
        call Message%printError('readNameList:',' experimental file name is undefined in '//nmlfile)
    end if

    if (trim(datafile).eq.'undefined') then
        call Message%printError('readNameList:',' datafile file name is undefined in '//nmlfile)
    end if

    if (trim(masterfile).eq.'undefined') then
        call Message%printError('readNameList:',' masterfile file name is undefined in '//nmlfile)
    end if
end if  

self%nml%ipf_wd = ipf_wd
self%nml%ipf_ht = ipf_ht
self%nml%nnk = nnk
self%nml%nosm = nosm
self%nml%ncubochoric = ncubochoric
self%nml%nsteps = nsteps
self%nml%nbatch = nbatch
self%nml%nthreads = nthreads
self%nml%tiltaxis = tiltaxis
self%nml%stepX = stepX
self%nml%stepY = stepY
self%nml%omega = omega
self%nml%omega_step = omega_step
self%nml%HDFstrings = HDFstrings
self%nml%exptfile = exptfile
self%nml%datafile = datafile
self%nml%ctffile = ctffile
self%nml%angfile = angfile
self%nml%masterfile = masterfile

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/27/23
!!
!! pass the namelist for the ISEDI_T Class to the calling program

IMPLICIT NONE 

class(ISEDI_T), INTENT(INOUT)          :: self
type(ISEDINameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/27/23
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)        :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 11, n_real = 9
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( mcnl => self%nml )

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine setipf_wd_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setipf_wd_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set ipf_wd in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%ipf_wd = inp

end subroutine setipf_wd_

!--------------------------------------------------------------------------
function getipf_wd_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getipf_wd_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get ipf_wd from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%ipf_wd

end function getipf_wd_

!--------------------------------------------------------------------------
subroutine setipf_ht_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setipf_ht_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set ipf_ht in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%ipf_ht = inp

end subroutine setipf_ht_

!--------------------------------------------------------------------------
function getipf_ht_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getipf_ht_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get ipf_ht from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%ipf_ht

end function getipf_ht_

!--------------------------------------------------------------------------
subroutine setnnk_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnnk_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set nnk in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%nnk = inp

end subroutine setnnk_

!--------------------------------------------------------------------------
function getnnk_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnnk_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get nnk from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%nnk

end function getnnk_

!--------------------------------------------------------------------------
subroutine setnosm_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnosm_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set nosm in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%nosm = inp

end subroutine setnosm_

!--------------------------------------------------------------------------
function getnosm_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnosm_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get nosm from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%nosm

end function getnosm_

!--------------------------------------------------------------------------
subroutine setncubochoric_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setncubochoric_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set ncubochoric in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%ncubochoric = inp

end subroutine setncubochoric_

!--------------------------------------------------------------------------
function getncubochoric_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getncubochoric_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get ncubochoric from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%ncubochoric

end function getncubochoric_

!--------------------------------------------------------------------------
subroutine setnsteps_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnsteps_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set nsteps in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%nsteps = inp

end subroutine setnsteps_

!--------------------------------------------------------------------------
function getnsteps_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnsteps_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get nsteps from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%nsteps

end function getnsteps_

!--------------------------------------------------------------------------
subroutine setnbatch_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnbatch_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set nbatch in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%nbatch = inp

end subroutine setnbatch_

!--------------------------------------------------------------------------
function getnbatch_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnbatch_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get nbatch from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%nbatch

end function getnbatch_

!--------------------------------------------------------------------------
subroutine setnthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnthreads_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set nthreads in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%nthreads = inp

end subroutine setnthreads_

!--------------------------------------------------------------------------
function getnthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnthreads_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get nthreads from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%nthreads

end function getnthreads_

!--------------------------------------------------------------------------
subroutine setstepX_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setstepX_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set stepX in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%stepX = inp

end subroutine setstepX_

!--------------------------------------------------------------------------
function getstepX_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getstepX_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get stepX from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%stepX

end function getstepX_

!--------------------------------------------------------------------------
subroutine setstepY_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setstepY_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set stepY in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%stepY = inp

end subroutine setstepY_

!--------------------------------------------------------------------------
function getstepY_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getstepY_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get stepY from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%stepY

end function getstepY_

!--------------------------------------------------------------------------
subroutine settiltaxis_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: settiltaxis_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set tiltaxis in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp(3)

self%nml%tiltaxis = inp

end subroutine settiltaxis_

!--------------------------------------------------------------------------
function gettiltaxis_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: gettiltaxis_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get tiltaxis from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out(3)

out = self%nml%tiltaxis

end function gettiltaxis_

!--------------------------------------------------------------------------
subroutine setomega_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setomega_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set omega in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%omega = inp

end subroutine setomega_

!--------------------------------------------------------------------------
function getomega_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getomega_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get omega from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%omega

end function getomega_

!--------------------------------------------------------------------------
subroutine setomega_step_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setomega_step_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set omega_step in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%omega_step = inp

end subroutine setomega_step_

!--------------------------------------------------------------------------
function getomega_step_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getomega_step_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get omega_step from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%omega_step

end function getomega_step_

!--------------------------------------------------------------------------
subroutine setHDFstrings_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setHDFstrings_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set HDFstrings in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)      :: inp(10)

integer(kind=irg)                 :: i 

do i=1,10
  self%nml%HDFstrings(i) = trim(inp(i))
end do 

end subroutine setHDFstrings_

!--------------------------------------------------------------------------
function getHDFstrings_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getHDFstrings_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get HDFstrings from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)   :: self
character(fnlen)                :: out(10)

integer(kind=irg)               :: i 

do i=1,10
  out(i) = trim(self%nml%HDFstrings(i))
end do

end function getHDFstrings_

!--------------------------------------------------------------------------
subroutine setexptfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setexptfile_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set exptfile in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%exptfile = trim(inp)

end subroutine setexptfile_

!--------------------------------------------------------------------------
function getexptfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getexptfile_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get exptfile from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%exptfile)

end function getexptfile_

!--------------------------------------------------------------------------
subroutine setdatafile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdatafile_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set datafile in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%datafile = trim(inp)

end subroutine setdatafile_

!--------------------------------------------------------------------------
function getdatafile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdatafile_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get datafile from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%datafile)

end function getdatafile_

!--------------------------------------------------------------------------
subroutine setctffile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setctffile_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set ctffile in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%ctffile = trim(inp)

end subroutine setctffile_

!--------------------------------------------------------------------------
function getctffile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getctffile_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get ctffile from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%ctffile)

end function getctffile_

!--------------------------------------------------------------------------
subroutine setangfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setangfile_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set angfile in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%angfile = trim(inp)

end subroutine setangfile_

!--------------------------------------------------------------------------
function getangfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getangfile_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get angfile from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%angfile)

end function getangfile_

!--------------------------------------------------------------------------
subroutine setmasterfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setmasterfile_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! set masterfile in the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%masterfile = trim(inp)

end subroutine setmasterfile_

!--------------------------------------------------------------------------
function getmasterfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getmasterfile_
!! author: MDG
!! version: 1.0
!! date: 02/27/23
!!
!! get masterfile from the ISEDI class

IMPLICIT NONE

class(ISEDI_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%masterfile)

end function getmasterfile_

!--------------------------------------------------------------------------
subroutine ISEDI_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: ISEDI_
!! author: MDG 
!! version: 1.0 
!! date: 02/27/23
!!
!! perform the computations

use mod_EMsoft
use mod_io
use mod_vendors
use mod_Lambert
use mod_others
use mod_others
use HDF5
use mod_memory
use mod_HDFsupport
use mod_HDFnames
use mod_quaternions
use mod_rotations
use mod_so3
use mod_ISEmaster
use mod_ISE
use stringconstants
use mod_crystallography
use mod_symmetry
use omp_lib
use ISO_C_BINDING
use mod_IPF
use mod_IPFsupport

IMPLICIT NONE 

class(ISEDI_T), INTENT(INOUT)       :: self
type(EMsoft_T), INTENT(INOUT)       :: EMsoft
character(fnlen), INTENT(INOUT)     :: progname 

type(IO_T)                          :: Message
type(HDF_T)                         :: HDF
type(Vendor_T)                      :: VT
type(HDFnames_T)                    :: HDFnames
type(so3_T)                         :: SO
type(ISE_T)                         :: ISE
type(Cell_T)                        :: cell
type(SpaceGroup_T)                  :: SG
type(memory_T)                      :: mem
type(e_T)                           :: eu 
type(a_T)                           :: ax
type(q_T)                           :: q
type(QuaternionArray_T)             :: Qartilt, qudictarray, qAR, sym, qAR2
type(Quaternion_T)                  :: quat, qu
type(IPF_T)                         :: IPF 
type(IPFmap_T)                      :: IPFmap 

type(ISEmasterNameListType)         :: mpnml

integer(kind=irg)                   :: io_int(2), pgnum, FZcnt, ii, jj, kk, ng, hdferr, npat, icnt, jjj, &
                                        ninbatch, nbatches, nremainder, pinbatch, pbatches, premainder, jstart, kstart
integer(kind=irg)                   :: nix, niy, nixp, niyp, TID, limit
integer(HSIZE_T)                    :: dims(3)
real(kind=sgl)                      :: dx, dy, dxm, dym, scl
real(kind=sgl)                      :: dc(3), nfactor
type(FZpointd),pointer              :: FZlist, FZtmp
real(kind=dbl), allocatable         :: ISEdict(:,:)
real(kind=sgl),allocatable          :: mLPNH(:,:), mLPSH(:,:), ISEimage(:,:,:), dp(:,:), minsortarr(:), resultarray(:), &
                                        dictblock(:,:), exptblock(:,:), dps(:,:), expt(:,:), maxsortarr(:), resultmain(:,:), &
                                        resulttmp(:,:)
integer(kind=irg),allocatable       :: dplabel(:,:), indexlist(:), indexarray(:), indexmain(:,:), indextmp(:,:)
real(kind=dbl)                      :: tiltaxis(3), angle, ddc(3), s, m, sd
character(fnlen)                    :: fname, groupname, dataset, IPFmapfile, IPFmode

call setRotationPrecision('d')

associate( nml => self%nml )

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()
HDFnames = HDFnames_T()

! 1. read ISE master pattern file (HDF format)
call HDFnames%set_ProgramData(SC_ISEmaster)
call HDFnames%set_NMLlist(SC_ISEmasterNameList)
call HDFnames%set_NMLfilename(SC_ISEmasterNML)
fname = EMsoft%generateFilePath('EMdatapathname',trim(nml%masterfile))
ISE = ISE_T()
call ISE%setISEMPfile(fname)
call ISE%readISEMPfile(HDF, HDFnames, mpnml, getmLPNH=.TRUE., getmLPSH=.TRUE.)
allocate(mLPNH(-mpnml%npx:mpnml%npx,-mpnml%npx:mpnml%npx))
allocate(mLPSH(-mpnml%npx:mpnml%npx,-mpnml%npx:mpnml%npx))
mLPNH = ISE%getmLPNH(mpnml%npx)
mLPSH = ISE%getmLPSH(mpnml%npx)

! 2. get the crystal structure information and point group number
call cell%getCrystalData(mpnml%xtalname, SG, EMsoft, verbose=.TRUE.)

! extract the point group number
pgnum = SG%getPGnumber()
io_int = pgnum
call Message%WriteValue(' Setting point group number to ',io_int,1)

! 3. generate an orientation sampling
SO = so3_T(pgnum, zerolist='FZ')
call SO%sampleRFZ(nml%ncubochoric)
FZcnt = SO%getListCount('FZ')
call SO%listtoQuaternionArray( qudictarray, 'FZ' )

! 4. pre-compute the entire dictionary
call Message%printMessage(' Computing dictionary of ISE intensity multiplets')

! for each orientation in the list, we generate the nsteps tilted orientations
! and sample the ISE intensities from the master pattern
mem = memory_T()
call mem%alloc(ISEdict, (/ nml%nsteps, FZcnt /), 'ISEdict')

! generate the quaternion list for the sample tilts 
Qartilt = QuaternionArray_T( n=nml%nsteps, s='d' )
tiltaxis = dble(nml%tiltaxis)
do ii = 1, nml%nsteps
  angle = nml%omega + dble(ii-1) * dble(nml%omega_step)
  ax = a_T( adinp = (/ tiltaxis(1), tiltaxis(2), tiltaxis(3), cvtoRadians(angle) /) )
  q = ax%aq()
  ! call q%q_print( 'quat : ' )
  quat = Quaternion_T( qd = q%q_copyd() )
  call Qartilt%insertQuatinArray( ii, quat )
end do

scl = float(mpnml%npx) 

! loop over all the dictionary orientations
call OMP_SET_NUM_THREADS(nml%nthreads)

!$OMP PARALLEL default(shared) private(ii,jj,s,quat,dc,nix,niy,nixp,niyp,dx,dy,dxm,dym)&
!$OMP& private(qu, ddc, io_int, m, sd ) 

TID = OMP_GET_THREAD_NUM()

!$OMP DO SCHEDULE(DYNAMIC)
do ii = 1, FZcnt
! get the orientation quaternion
  quat = qudictarray%getQuatfromArray( ii )

  do jj=1, nml%nsteps
    s = 0.D0
! get the pixel direction cosines 
    ddc = (/ 0.D0, 0.D0, 1.D0 /)
! apply the sample tilt to the detector direction cosines
    qu = Qartilt%getQuatfromArray( jj )
    ddc = qu%quat_Lp( ddc )
! apply the grain rotation to the detector direction cosines
    ddc = quat%quat_Lp( ddc )
    dc = sngl(ddc/sqrt(sum(ddc*ddc)))
! convert these direction cosines to interpolation coordinates in the Rosca-Lambert projection
    dc = dc / sqrt(sum(dc*dc))
    call LambertgetInterpolation(dc, scl, mpnml%npx, mpnml%npx, nix, niy, nixp, niyp, dx, dy, dxm, dym)

    if (dc(3) .ge. 0.0) then
          s = mLPNH(nix,niy) * dxm * dym + &
              mLPNH(nixp,niy) * dx * dym + &
              mLPNH(nix,niyp) * dxm * dy + &
              mLPNH(nixp,niyp) * dx * dy 
    else
          s = mLPSH(nix,niy) * dxm * dym + &
              mLPSH(nixp,niy) * dx * dym + &
              mLPSH(nix,niyp) * dxm * dy + &
              mLPSH(nixp,niyp) * dx * dy 
    end if
    ISEdict(jj,ii) = s 
  end do
  if (mod(ii,10000).eq.0) then 
    io_int(1) = ii
    call Message%WriteValue(' completed orientation ', io_int, 1)
  end if
! subtract mean and divide by standard deviation
  m = sum( ISEdict(:,ii) ) / dble(nml%nsteps)
  sd = sqrt( sum( (ISEdict(:,ii)-m)**2 ) / dble(nml%nsteps) )
  ISEdict(:,ii) = (ISEdict(:,ii)-m)/sd
end do 
!$OMP END DO
!$OMP END PARALLEL

call Message%printMessage('  --> Dictionary completed ')

! 5. load the experimental ISE images from their HDF5 file

! get the number of HDF groups to traverse 
ng = get_num_HDFgroups_(nml%HDFstrings)

fname = EMsoft%generateFilePath('EMdatapathname',trim(nml%exptfile))

hdferr =  HDF%openFile(fname, readonly=.TRUE.)
if (hdferr.ne.0) call HDF%error_check('openExpPatternFile:HDF%openFile', hdferr)
! open all the groups to the correct level of the data set
do ii=1,ng
    groupname = trim(nml%HDFstrings(ii))
    hdferr = HDF%openGroup(groupname)
    if (hdferr.ne.0) call HDF%error_check('openExpPatternFile:HDF%openGroup: groupname issue, check for typos.', hdferr)
end do

dataset = trim(nml%HDFstrings(ng+1))
 call HDF%readDatasetFloatArray(dataset, dims, hdferr, ISEimage)
 if (hdferr.ne.0) &
   call HDF%error_check('HDF%readDatasetIntegerArray: problem reading ISEimages array', hdferr)

call HDF%pop(.TRUE.)

! normalize the image intensities and reduce the array dimensions
npat =  nml%ipf_wd*nml%ipf_ht
call mem%alloc(expt, (/ nml%nsteps, npat /), 'expt')
icnt = 1

do ii=1,nml%ipf_ht
  do jj=1,nml%ipf_wd
! subtract mean and divide by standard deviation
    m = sum( ISEimage(jj, ii, :) ) / dble(nml%nsteps)
    sd = sqrt( sum( (ISEimage(jj,ii,:)-m)**2 ) / dble(nml%nsteps) )
    expt(:,icnt) = (ISEimage(jj,ii,:)-m)/sd
    icnt = icnt+1
  end do
end do

call Message%printMessage(' Read and normalized experimental ISE images ')

! 6. compute all the dot products and then rank them 
call mem%alloc(dp, (/ npat, nml%nnk /), 'dp')
call mem%alloc(dplabel, (/ npat, nml%nnk /), 'dplabel')
call mem%alloc(dictblock, (/ nml%nbatch, nml%nsteps /), 'dictblock')
call mem%alloc(exptblock, (/ nml%nsteps, nml%nbatch /), 'exptblock')
call mem%alloc(dps, (/ nml%nbatch, nml%nbatch /), 'exptblock')
call mem%alloc(indexlist, (/ FZcnt /), 'indexlist')
do ii=1,FZcnt 
  indexlist(ii) = ii 
end do
call mem%alloc(maxsortarr, (/ npat /), 'maxsortarr', initval = 0.0)
call mem%alloc(minsortarr, (/ npat /), 'minsortarr', initval =-2.0)
call mem%alloc(resultarray, (/ nml%nbatch /), 'resultarray')
call mem%alloc(indexarray, (/ nml%nbatch /), 'indexarray')
call mem%alloc(resultmain, (/ nml%nnk, npat /), 'resultmain', initval = -2.0)
call mem%alloc(indexmain, (/ nml%nnk,npat /), 'indexmain', initval = 0)
call mem%alloc(resulttmp, (/ 2*nml%nnk,npat /), 'resulttmp', initval = -2.0)
call mem%alloc(indextmp, (/ 2*nml%nnk, npat /), 'indextmp', initval = 0)

nfactor = 1.0/ float(nml%nsteps)

! we'll chunk the dictionary in batches of nbatch vectors with a remainder 
ninbatch = nml%nbatch 
nbatches = int(FZcnt/ninbatch)
nremainder = FZcnt - nbatches*ninbatch
pinbatch = nml%nbatch 
pbatches = int(npat/pinbatch)
premainder = npat - pbatches*pinbatch

io_int(1) = nbatches+1 
call Message%WriteValue(' Number of batches to index : ', io_int, 1)


outerloop: do ii=1,nbatches+1    ! loop over the dictionary
! fill the dictblock, taking into account that the last one is a different size
  dictblock = 0.0
  jstart = (ii-1)*ninbatch
  if (ii.ne.nbatches+1) then ! a regular dictblock
    do jj=1,ninbatch
      dictblock(jj,:) = ISEdict(:, jstart + jj)
    end do
  else   ! this is the remainder part 
    do jj=1,nremainder
      dictblock(jj,:) = ISEdict(:, jstart + jj)
    end do
  end if 
! next loop over all the exptblocks to compute the dps array
  innerloop: do kk=1,pbatches+1  ! loop over the experimental patterns
 ! fill the dictblock, taking into account that the last one is a different size
    exptblock = 0.0
    kstart = (kk-1)*pinbatch
    if (kk.ne.pbatches+1) then ! a regular exptblock
      do jj=1,pinbatch
        exptblock(:, jj) = expt(:, kstart + jj)
      end do
    else   ! this is the remainder part 
      if (premainder.ne.0) then
        do jj=1,premainder
          exptblock(:, jj) = expt(:, kstart + jj)
        end do
      end if 
    end if 
    ! write (*,*) 'min/max ',minval(dps), maxval(dps), maxval(indexlist)
! next, rank the dot products and keep the top nnk for each experimental pattern 
! the dps array has the dictionary patterns index as the first index, and the 
! experimental pattern index as the second index. 
    if (kk.eq.pbatches+1) then 
      limit = premainder
    else
      limit = pinbatch 
    end if 
    if (limit.ne.0) then 
! compute the dot product array
      dps = matmul(dictblock,exptblock) * nfactor 
      do jj=1,limit
        jjj = (kk-1)*pinbatch + jj
        maxsortarr(jjj) = maxval(dps(:,jj))
        if (maxsortarr(jjj).gt.minsortarr(jjj)) then ! only sort if the max falls in the range
          resultarray(1:ninbatch) = dps(1:ninbatch,jj)
          indexarray(1:ninbatch) = indexlist((ii-1)*ninbatch+1:ii*ninbatch)

          call SSORT(resultarray,indexarray,ninbatch,-2)

          resulttmp(nml%nnk+1:2*nml%nnk,jjj) = resultarray(1:nml%nnk)
          indextmp(nml%nnk+1:2*nml%nnk,jjj) = indexarray(1:nml%nnk)

          call SSORT(resulttmp(:,jjj),indextmp(:,jjj),2*nml%nnk,-2)

          resultmain(1:nml%nnk,jjj) = resulttmp(1:nml%nnk,jjj)
          indexmain(1:nml%nnk,jjj) = indextmp(1:nml%nnk,jjj)
          minsortarr(jjj) = resulttmp(nml%nnk,jjj)
        end if 
      end do
    end if 
  end do innerloop
  if (mod(ii,1).eq.0) then 
    io_int(1) = ii 
    io_int(2) = nbatches+1 
    call Message%WriteValue(' completed batch # ', io_int, 2, "(I4,' of ',I4)")
  end if
end do outerloop


qAR2 = QuaternionArray_T( n=1, s='d' )
call qAR2%QSym_Init(pgnum, sym)

qAR = QuaternionArray_T( n=npat, s='d' )
do icnt=1,npat
  quat = qudictarray%getQuatfromArray( indexmain(1,icnt) )
  call qAR%insertQuatinArray( icnt, quat )
end do 

IPF = IPF_T()
IPFmapfile = 'currentIPFZmap.tiff'
call IPF%set_IPFfilename(IPFmapfile)
call IPF%set_sampleDir( (/ 0, 0, 1 /) )
call IPF%set_nthreads(1)
IPFmode = 'TSL'
call IPF%set_IPFmode(IPFmode)
call IPF%updateIPFmap(EMsoft, progname, nml%ipf_wd, nml%ipf_ht, pgnum, IPFmapfile, qAR, sym) 


end associate

! and conclude with explicit memory deallocations 
call mem%dealloc(ISEdict, 'ISEdict')


end subroutine ISEDI_

end module mod_ISEDI