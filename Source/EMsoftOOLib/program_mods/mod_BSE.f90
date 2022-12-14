!! ###################################################################
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

module mod_BSE
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/04/22
  !!
  !! class definition for the EMBSE program

use mod_kinds
use mod_global

IMPLICIT NONE 

type, private :: BSEAccumType
  integer(kind=irg)               :: numdet
  integer(kind=irg)               :: ipf_wd
  integer(kind=irg)               :: ipf_ht
  integer(kind=irg),allocatable   :: accum_e(:,:,:),accum_z(:,:,:,:)
  real(kind=sgl),allocatable      :: accum_e_detector(:,:,:)
  real(kind=sgl),allocatable      :: mLPNH(:,:,:), mLPSH(:,:,:)
  real(kind=sgl),allocatable      :: EmLPNH(:,:), EmLPSH(:,:)
  real(kind=sgl),allocatable      :: rgx(:), rgy(:), rgz(:), corfactor(:)
  real(kind=sgl),allocatable      :: beamtiltq(:,:,:)
end type BSEAccumType

! namelist for the EMBSE program
type, public :: BSENameListType
  real(kind=sgl)    :: energymin
  real(kind=sgl)    :: energymax
  real(kind=sgl)    :: incidence
  real(kind=sgl)    :: beamcurrent
  real(kind=sgl)    :: dwelltime
  real(kind=sgl)    :: gammavalue
  real(kind=sgl)    :: workingdistance
  real(kind=sgl)    :: BSEdistance
  real(kind=sgl)    :: rin
  real(kind=sgl)    :: rout
  integer(kind=irg) :: NsqL
  integer(kind=irg) :: nthreads
  character(fnlen)  :: useangles
  character(3)      :: scalingmode
  character(6)      :: scanmode
  character(fnlen)  :: masterfile
  character(fnlen)  :: Kosselmasterfile
  character(fnlen)  :: datafile
  character(fnlen)  :: imagefile
end type BSENameListType



! class definition
type, public :: BSE_T
private 
  character(fnlen)        :: nmldeffile = 'EMBSE.nml'
  type(BSENameListType)   :: nml 
  type(BSEAccumType)      :: det

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: BSE_
  procedure, pass(self) :: setenergymin_
  procedure, pass(self) :: getenergymin_
  procedure, pass(self) :: setenergymax_
  procedure, pass(self) :: getenergymax_
  procedure, pass(self) :: setincidence_
  procedure, pass(self) :: getincidence_
  procedure, pass(self) :: setbeamcurrent_
  procedure, pass(self) :: getbeamcurrent_
  procedure, pass(self) :: setdwelltime_
  procedure, pass(self) :: getdwelltime_
  procedure, pass(self) :: setgammavalue_
  procedure, pass(self) :: getgammavalue_
  procedure, pass(self) :: setworkingdistance_
  procedure, pass(self) :: getworkingdistance_
  procedure, pass(self) :: setBSEdistance_
  procedure, pass(self) :: getBSEdistance_
  procedure, pass(self) :: setrin_
  procedure, pass(self) :: getrin_
  procedure, pass(self) :: setrout_
  procedure, pass(self) :: getrout_
  procedure, pass(self) :: setNsqL_
  procedure, pass(self) :: getNsqL_
  procedure, pass(self) :: setnthreads_
  procedure, pass(self) :: getnthreads_
  procedure, pass(self) :: setuseangles_
  procedure, pass(self) :: getuseangles_
  procedure, pass(self) :: setscalingmode_
  procedure, pass(self) :: getscalingmode_
  procedure, pass(self) :: setscanmode_
  procedure, pass(self) :: getscanmode_  
  procedure, pass(self) :: setmasterfile_
  procedure, pass(self) :: getmasterfile_
  procedure, pass(self) :: setKosselmasterfile_
  procedure, pass(self) :: getKosselmasterfile_
  procedure, pass(self) :: setdatafile_
  procedure, pass(self) :: getdatafile_
  procedure, pass(self) :: setimagefile_
  procedure, pass(self) :: getimagefile_
  procedure, pass(self) :: GenerateBSEDetector_
  procedure, pass(self) :: GenerateBSEbeamtiltquaternions_
  procedure, pass(self) :: ComputeBSEimage_
  procedure, pass(self) :: ComputeBSErings_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: BSE => BSE_
  generic, public :: setenergymin => setenergymin_
  generic, public :: getenergymin => getenergymin_
  generic, public :: setenergymax => setenergymax_
  generic, public :: getenergymax => getenergymax_
  generic, public :: setincidence => setincidence_
  generic, public :: getincidence => getincidence_
  generic, public :: setbeamcurrent => setbeamcurrent_
  generic, public :: getbeamcurrent => getbeamcurrent_
  generic, public :: setdwelltime => setdwelltime_
  generic, public :: getdwelltime => getdwelltime_
  generic, public :: setgammavalue => setgammavalue_
  generic, public :: getgammavalue => getgammavalue_
  generic, public :: setworkingdistance => setworkingdistance_
  generic, public :: getworkingdistance => getworkingdistance_
  generic, public :: setBSEdistance => setBSEdistance_
  generic, public :: getBSEdistance => getBSEdistance_
  generic, public :: setrin => setrin_
  generic, public :: getrin => getrin_
  generic, public :: setrout => setrout_
  generic, public :: getrout => getrout_
  generic, public :: setNsqL => setNsqL_
  generic, public :: getNsqL => getNsqL_
  generic, public :: setnthreads => setnthreads_
  generic, public :: getnthreads => getnthreads_
  generic, public :: setuseangles => setuseangles_
  generic, public :: getuseangles => getuseangles_
  generic, public :: setscalingmode => setscalingmode_
  generic, public :: getscalingmode => getscalingmode_
  generic, public :: setscanmode => setscanmode_
  generic, public :: getscanmode => getscanmode_  
  generic, public :: setmasterfile => setmasterfile_
  generic, public :: getmasterfile => getmasterfile_
  generic, public :: setKosselmasterfile => setKosselmasterfile_
  generic, public :: getKosselmasterfile => getKosselmasterfile_
  generic, public :: setdatafile => setdatafile_
  generic, public :: getdatafile => getdatafile_
  generic, public :: setimagefile => setimagefile_
  generic, public :: getimagefile => getimagefile_
  generic, public :: GenerateBSEDetector => GenerateBSEDetector_
  generic, public :: GenerateBSEbeamtiltquaternions => GenerateBSEbeamtiltquaternions_
  generic, public :: ComputeBSEimage => ComputeBSEimage_
  generic, public :: ComputeBSErings => ComputeBSErings_

end type BSE_T

! the constructor routine for this class 
interface BSE_T
  module procedure BSE_constructor
end interface BSE_T

contains

!--------------------------------------------------------------------------
type(BSE_T) function BSE_constructor( nmlfile ) result(BSE)
!! author: MDG 
!! version: 1.0 
!! date: 12/04/22
!!
!! constructor for the BSE_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call BSE%readNameList(nmlfile)

end function BSE_constructor

!--------------------------------------------------------------------------
subroutine BSE_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 12/04/22
!!
!! destructor for the BSE_T Class
 
IMPLICIT NONE

type(BSE_T), INTENT(INOUT)  :: self 

call reportDestructor('BSE_T')

end subroutine BSE_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 12/04/22
!!
!! read the namelist from an nml file for the BSE_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(BSE_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

real(kind=sgl)    :: energymin
real(kind=sgl)    :: energymax
real(kind=sgl)    :: incidence
real(kind=sgl)    :: beamcurrent
real(kind=sgl)    :: dwelltime
real(kind=sgl)    :: gammavalue
real(kind=sgl)    :: workingdistance
real(kind=sgl)    :: BSEdistance
real(kind=sgl)    :: rin
real(kind=sgl)    :: rout
integer(kind=irg) :: NsqL
integer(kind=irg) :: nthreads
character(fnlen)  :: useangles
character(3)      :: scalingmode
character(6)      :: scanmode
character(fnlen)  :: masterfile
character(fnlen)  :: Kosselmasterfile
character(fnlen)  :: datafile
character(fnlen)  :: imagefile

! define the IO namelist to facilitate passing variables to the program.
namelist  / BSEdata / energymin, energymax, incidence, beamcurrent, dwelltime, gammavalue, workingdistance, &
                      BSEdistance, rin, rout, NsqL, nthreads, useangles, scalingmode, masterfile, &
                      Kosselmasterfile, datafile, imagefile, scanmode

! set the input parameters to default values
energymin = 5.0
energymax = 20.0
masterfile = 'undefined'
Kosselmasterfile = 'undefined'
datafile = 'undefined'
useangles = 'original'
incidence = 0.0
beamcurrent = 150.0
dwelltime = 100.0
scalingmode = 'not'
scanmode = 'single'
gammavalue = 1.0
workingdistance = 10.0
BSEdistance = 9.5
rin = 5.0
rout = 12.0
NsqL = 40
imagefile = 'undefined'
nthreads = 1

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=BSEdata)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(masterfile).eq.'undefined') then
  call Message%printError('readNameList:',' master pattern file name is undefined in '//nmlfile)
 end if

 if (trim(datafile).eq.'undefined') then
  call Message%printError('readNameList:',' output file name is undefined in '//nmlfile)
 end if

 if (trim(imagefile).eq.'undefined') then
  call Message%printError('readNameList:',' image file name is undefined in '//nmlfile)
 end if
end if 

self%nml%energymin = energymin
self%nml%energymax = energymax
self%nml%incidence = incidence
self%nml%beamcurrent = beamcurrent
self%nml%dwelltime = dwelltime
self%nml%gammavalue = gammavalue
self%nml%workingdistance = workingdistance
self%nml%BSEdistance = BSEdistance
self%nml%rin = rin
self%nml%rout = rout
self%nml%NsqL = NsqL
self%nml%nthreads = nthreads
self%nml%useangles = useangles
self%nml%scalingmode = scalingmode
self%nml%scanmode = scanmode
self%nml%masterfile = masterfile
self%nml%Kosselmasterfile = Kosselmasterfile
self%nml%datafile = datafile
self%nml%imagefile = imagefile

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 12/04/22
!!
!! pass the namelist for the BSE_T Class to the calling program

IMPLICIT NONE 

class(BSE_T), INTENT(INOUT)          :: self
type(BSENameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 12/04/22
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)        :: self 
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
subroutine setenergymin_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setenergymin_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set energymin in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%energymin = inp

end subroutine setenergymin_

!--------------------------------------------------------------------------
function getenergymin_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getenergymin_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get energymin from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%energymin

end function getenergymin_

!--------------------------------------------------------------------------
subroutine setenergymax_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setenergymax_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set energymax in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%energymax = inp

end subroutine setenergymax_

!--------------------------------------------------------------------------
function getenergymax_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getenergymax_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get energymax from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%energymax

end function getenergymax_

!--------------------------------------------------------------------------
subroutine setincidence_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setincidence_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set incidence in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%incidence = inp

end subroutine setincidence_

!--------------------------------------------------------------------------
function getincidence_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getincidence_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get incidence from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%incidence

end function getincidence_

!--------------------------------------------------------------------------
subroutine setbeamcurrent_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setbeamcurrent_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set beamcurrent in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%beamcurrent = inp

end subroutine setbeamcurrent_

!--------------------------------------------------------------------------
function getbeamcurrent_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getbeamcurrent_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get beamcurrent from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%beamcurrent

end function getbeamcurrent_

!--------------------------------------------------------------------------
subroutine setdwelltime_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdwelltime_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set dwelltime in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%dwelltime = inp

end subroutine setdwelltime_

!--------------------------------------------------------------------------
function getdwelltime_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdwelltime_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get dwelltime from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%dwelltime

end function getdwelltime_

!--------------------------------------------------------------------------
subroutine setgammavalue_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setgammavalue_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set gammavalue in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%gammavalue = inp

end subroutine setgammavalue_

!--------------------------------------------------------------------------
function getgammavalue_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getgammavalue_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get gammavalue from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%gammavalue

end function getgammavalue_

!--------------------------------------------------------------------------
subroutine setworkingdistance_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setworkingdistance_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set workingdistance in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%workingdistance = inp

end subroutine setworkingdistance_

!--------------------------------------------------------------------------
function getworkingdistance_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getworkingdistance_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get workingdistance from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%workingdistance

end function getworkingdistance_

!--------------------------------------------------------------------------
subroutine setBSEdistance_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setBSEdistance_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set BSEdistance in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%BSEdistance = inp

end subroutine setBSEdistance_

!--------------------------------------------------------------------------
function getBSEdistance_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getBSEdistance_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get BSEdistance from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%BSEdistance

end function getBSEdistance_

!--------------------------------------------------------------------------
subroutine setrin_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setrin_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set rin in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%rin = inp

end subroutine setrin_

!--------------------------------------------------------------------------
function getrin_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getrin_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get rin from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%rin

end function getrin_

!--------------------------------------------------------------------------
subroutine setrout_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setrout_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set rout in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%rout = inp

end subroutine setrout_

!--------------------------------------------------------------------------
function getrout_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getrout_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get rout from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%rout

end function getrout_

!--------------------------------------------------------------------------
subroutine setNsqL_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setNsqL_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set NsqL in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%NsqL = inp

end subroutine setNsqL_

!--------------------------------------------------------------------------
function getNsqL_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getNsqL_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get NsqL from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%NsqL

end function getNsqL_

!--------------------------------------------------------------------------
subroutine setnthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnthreads_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set nthreads in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%nthreads = inp

end subroutine setnthreads_

!--------------------------------------------------------------------------
function getnthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnthreads_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get nthreads from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%nthreads

end function getnthreads_

!--------------------------------------------------------------------------
subroutine setuseangles_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setuseangles_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set useangles in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%useangles = trim(inp)

end subroutine setuseangles_

!--------------------------------------------------------------------------
function getuseangles_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getuseangles_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get useangles from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%useangles)

end function getuseangles_

!--------------------------------------------------------------------------
subroutine setscalingmode_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setscalingmode_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set scalingmode in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
character(3), INTENT(IN)       :: inp

self%nml%scalingmode = trim(inp)

end subroutine setscalingmode_

!--------------------------------------------------------------------------
function getscalingmode_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getscalingmode_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get scalingmode from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
character(3)                   :: out

out = trim(self%nml%scalingmode)

end function getscalingmode_

!--------------------------------------------------------------------------
subroutine setscanmode_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setscanmode_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set scanmode in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
character(3), INTENT(IN)       :: inp

self%nml%scanmode = trim(inp)

end subroutine setscanmode_

!--------------------------------------------------------------------------
function getscanmode_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getscanmode_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get scanmode from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
character(3)                   :: out

out = trim(self%nml%scanmode)

end function getscanmode_


!--------------------------------------------------------------------------
subroutine setmasterfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setmasterfile_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set masterfile in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%masterfile = trim(inp)

end subroutine setmasterfile_

!--------------------------------------------------------------------------
function getmasterfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getmasterfile_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get masterfile from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%masterfile)

end function getmasterfile_

!--------------------------------------------------------------------------
subroutine setKosselmasterfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setKosselmasterfile_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set Kosselmasterfile in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%Kosselmasterfile = trim(inp)

end subroutine setKosselmasterfile_

!--------------------------------------------------------------------------
function getKosselmasterfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getKosselmasterfile_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get Kosselmasterfile from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%Kosselmasterfile)

end function getKosselmasterfile_

!--------------------------------------------------------------------------
subroutine setdatafile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdatafile_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set datafile in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%datafile = trim(inp)

end subroutine setdatafile_

!--------------------------------------------------------------------------
function getdatafile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdatafile_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get datafile from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%datafile)

end function getdatafile_

!--------------------------------------------------------------------------
subroutine setimagefile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setimagefile_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! set imagefile in the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%imagefile = trim(inp)

end subroutine setimagefile_

!--------------------------------------------------------------------------
function getimagefile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getimagefile_
!! author: MDG
!! version: 1.0
!! date: 12/05/22
!!
!! get imagefile from the BSE_T class

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%imagefile)

end function getimagefile_

!--------------------------------------------------------------------------
subroutine GenerateBSEDetector_(self, sigma, verbose, SO2) 
!DEC$ ATTRIBUTES DLLEXPORT :: GenerateBSEDetector_
!! author: MDG 
!! version: 1.0 
!! date: 12/06/22
!!
!! generate BSE detector direction cosine arrays

use mod_so2
use mod_io
use mod_quaternions
use mod_rotations

class(BSE_T), INTENT(INOUT)             :: self
real(kind=dbl),INTENT(IN)               :: sigma
logical,INTENT(IN),OPTIONAL             :: verbose
type(SO2_T), OPTIONAL                   :: SO2

type(so2_T)                             :: SO
type(IO_T)                              :: Message
type(quaternion_T)                      :: qu

real(kind=dbl)                          :: BSEni(3), BSEno(3), theta, v(3), dtor, io_dble(3)
type(SO2pointd),pointer                 :: SO2list, SO2tmp
integer(kind=irg)                       :: i, io_int(1), SO2cnt

call setRotationPrecision('d')

associate( BSEdetector => self%det )

dtor = cPi/180.D0 

BSEni = (/ dble(self%nml%rin), 0.D0, dble(self%nml%BSEdistance) /)
BSEno = (/ dble(self%nml%rout), 0.D0, dble(self%nml%BSEdistance) /)
BSEni = BSEni/ sqrt(sum(BSEni*BSEni))
BSEno = BSEno/ sqrt(sum(BSEno*BSEno))

! create a linked list of unit vectors that point inside the annular BSE detector
if (self%nml%scanmode.eq.'single') then
  SO = SO2_T(self%nml%nsqL, BSEni(3), BSEno(3))
else ! scanmode='scan' leads to a square spiral ordering of the points in the linked list
  SO = SO2_T(self%nml%nsqL)
end if 

! put some stuff on the output
if (present(verbose)) then 
    if (verbose.eqv..TRUE.) then 
        io_int(1) = SO%getSO2cnt()
        call Message%WriteValue(' Number of BSE detector points : ',io_int,1)
        io_dble = BSEni
        call Message%WriteValue(' Inner detector vector ', io_dble, 3)
        io_dble = BSEno
        call Message%WriteValue(' Outer detector vector ', io_dble, 3)
    end if 
end if 

SO2cnt = SO%getSO2cnt()
BSEdetector%numdet = SO2cnt

if (allocated(BSEdetector%rgx)) deallocate(BSEdetector%rgx)
if (allocated(BSEdetector%rgy)) deallocate(BSEdetector%rgy)
if (allocated(BSEdetector%rgz)) deallocate(BSEdetector%rgz)
if (allocated(BSEdetector%corfactor)) deallocate(BSEdetector%corfactor)
allocate( BSEdetector%rgx(SO2cnt), BSEdetector%rgy(SO2cnt), BSEdetector%rgz(SO2cnt), BSEdetector%corfactor(SO2cnt) )

! loop over the linked list and rotate each vector onto a conical volume surrounding the incident beam
theta = - sigma * dtor * 0.5D0
qu = quaternion_T( qd = (/ cos(theta), sin(theta), 0.D0, 0.D0 /) )
SO2tmp => SO%getSO2listhead()
do i = 1, SO2cnt
    BSEdetector%corfactor(i) = 1.D0/SO2tmp%nvec(3)
! rotate the unit vector to the optical axis (opposite of sample tilt)
    v = qu%quat_Lp( SO2tmp%nvec ) 
    BSEdetector%rgx(i) = v(1)
    BSEdetector%rgy(i) = v(2)
    BSEdetector%rgz(i) = v(3)
! next point if it exists in the list     
    if (associated(SO2tmp)) then 
        SO2tmp => SO2tmp%next
    end if 
end do 

! normalize the area correction factor so that smallest value equals 1
BSEdetector%corfactor = BSEdetector%corfactor / minval(BSEdetector%corfactor)

! and clean up the linked list 
if (present(SO2)) then 
  SO2 = SO
else
  call SO%delete_SO2list()
end if

end associate

end subroutine GenerateBSEDetector_

!--------------------------------------------------------------------------
subroutine GenerateBSEbeamtiltquaternions_(self, dinl, verbose) 
!DEC$ ATTRIBUTES DLLEXPORT :: GenerateBSEbeamtiltquaternions_
!! author: MDG 
!! version: 1.0 
!! date: 12/06/22
!!
!! generate quaternions to describe the incident beam tilt in each point

use mod_io
use mod_DIfiles
use mod_quaternions
use mod_rotations

IMPLICIT NONE 

class(BSE_T), INTENT(INOUT)                 :: self
type(DictionaryIndexingNameListType),INTENT(IN)   :: dinl
logical,INTENT(IN),OPTIONAL                 :: verbose 

! type(Quaternion_T)                          :: qu  
type(IO_T)                                  :: Message

real(kind=sgl)                              :: ctrx, ctry, px, py, WD2, th, ct, st, p, maxth, io_real(4)
integer(kind=irg)                           :: ipf_wd, ipf_ht, pxstart, pystart, pxend, pyend, ix, iy 

call setRotationPrecision('d')

associate( enl => self%nml, BSEdetector => self%det )

WD2 = (enl%workingdistance*1000.0)**2

! get the size of the ROI as well as the central pixel of the complete IPF map
ctrx = BSEdetector%ipf_wd * 0.5    ! in pixels 
ctry = BSEdetector%ipf_ht * 0.5    ! in pixels

if (sum(dinl%ROI).eq.0) then 
    pxstart = 1
    pystart = 1
    pxend = BSEdetector%ipf_wd 
    pyend = BSEdetector%ipf_ht
    ipf_wd = BSEdetector%ipf_wd
    ipf_ht = BSEdetector%ipf_ht
else 
    pxstart = dinl%ROI(1)
    pystart = dinl%ROI(2)
    pxend = dinl%ROI(1) + dinl%ROI(3) - 1
    pyend = dinl%ROI(2) + dinl%ROI(4) - 1
    ipf_wd = dinl%ROI(3)
    ipf_ht = dinl%ROI(4)
end if

! allocate the beam tilt quaternion array 
if (allocated(BSEdetector%beamtiltq)) deallocate(BSEdetector%beamtiltq)
allocate(BSEdetector%beamtiltq(4, ipf_wd, ipf_ht))

! scan across the ROI and determine for each point the coordinates (px,py) in microns
! then convert that to a rotation quaternion for the beam tilt; we assume 
! that the beam pivots on the intersection point of the optical axis and the 
! objective lens exit plane (at distance equal to the working distance WD) 
! This is really a very small correction, in particular when the scan area is small
maxth = -1000.0
do ix=pxstart, pxend 
    px = ( ctrx - float(ix) ) * dinl%StepX
    do iy=pystart, pyend 
        py = -( ctry - float(iy) ) * dinl%StepY
        if ( (px.eq.0.0).and.(py.eq.0.0) ) then ! normal incidence in this point
            BSEdetector%beamtiltq(1:4, ix, iy) = (/ 1.0, 0.0 ,0.0 ,0.0 /)
        else  ! get the quaternion for this beam tilt 
            th = 0.5 * acos(enl%workingdistance*1000.0 / sqrt(px*px+py*py+WD2))  ! half the rotation angle 
            ! write(*,*) px, py, th
            if (th.gt.maxth) maxth = th
            ct = cos(th)
            st = sin(th)
            p = 1.0 / sqrt(px*px+py*py)
            BSEdetector%beamtiltq(1:4, ix, iy) = (/ ct, st * py * p, -st * px * p, 0.0 /)
        end if
    end do 
end do 

if (present(verbose)) then 
    if (verbose) then 
        io_real(1) = 2.0*maxth*180.0/sngl(cPi)
        call Message%WriteValue(' maximum beam tilt angle ', io_real, 1)
        io_real(1:4) = BSEdetector%beamtiltq( 1:4, pxstart, pystart)
        call Message%WriteValue(' (s,s) : ', io_real, 4) 
        io_real(1:4) = BSEdetector%beamtiltq( 1:4, pxend, pystart)
        call Message%WriteValue(' (e,s) : ', io_real, 4) 
        io_real(1:4) = BSEdetector%beamtiltq( 1:4, pxstart, pyend)
        call Message%WriteValue(' (s,e) : ', io_real, 4) 
        io_real(1:4) = BSEdetector%beamtiltq( 1:4, pxend, pyend)
        call Message%WriteValue(' (e,e) : ', io_real, 4) 
    end if 
end if 

end associate 

end subroutine GenerateBSEbeamtiltquaternions_

!--------------------------------------------------------------------------
subroutine ComputeBSEimage_(self, mcnl, mpnl, numang, Eangles, Emin, Emax, BSEimage)
!DEC$ ATTRIBUTES DLLEXPORT :: ComputeBSEimage_
!! author: MDG 
!! version: 1.0 
!! date: 12/06/22
!!
!! compute an energy-weighted BSE image for a given ROI

use mod_MCfiles
use mod_MPfiles
use mod_io
use mod_quaternions 
use mod_rotations
use mod_Lambert
use mod_image 
use omp_lib

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)                 :: self
type(MCOpenCLNameListType),INTENT(IN)       :: mcnl
type(EBSDMasterNameListType),INTENT(IN)     :: mpnl
integer(kind=irg),INTENT(IN)                :: numang 
real(kind=sgl),INTENT(IN)                   :: Eangles(3,numang)
integer(kind=irg),INTENT(IN)                :: Emin
integer(kind=irg),INTENT(IN)                :: Emax 
real(kind=sgl),INTENT(OUT),allocatable      :: BSEimage(:,:)

type(IO_T)                                  :: Message
type(Quaternion_T)                          :: quat, dquat 
type(e_T)                                   :: eu
type(q_T)                                   :: qu
type(Lambert_T)                             :: L

real(kind=sgl)                              :: s, mi, ma, ECPfactor, q(4)
integer(kind=irg)                           :: ix, iy, icnt, jd, sz(3), nxmc
real(kind=sgl)                              :: dc(3), avdc(3), newavdc(3), ixy(2), scl, sclmc
real(kind=dbl)                              :: ddc(3)
real(kind=sgl)                              :: dx, dy, dxm, dym, x, y, z
integer(kind=irg)                           :: ii, jj, kk, istat
integer(kind=irg)                           :: nix, niy, nixp, niyp, nixmc, niymc, TID

call setRotationPrecision('d')

associate( enl => self%nml, BSE => self%det )

allocate(BSEimage(BSE%ipf_wd,BSE%ipf_ht))
BSEimage = 0.0
scl = float(mpnl%npx) 
sz = shape(BSE%accum_e)
nxmc = (sz(2)-1)/2
sclmc = float(nxmc)

! THIS ALL NEEDS TO BE REVIEWED AND IMPROVED !!!

! loop over all the image pixels 
icnt = 0
call OMP_SET_NUM_THREADS(enl%nthreads)

!$OMP PARALLEL default(shared) private(ix,iy,s,icnt,qu,jd,kk,dc,nix,niy,nixp,niyp,dx,dy,dxm,dym,newavdc,ECPfactor)&
!$OMP& private(nixmc, niymc, eu, quat, dquat, ddc) 

TID = OMP_GET_THREAD_NUM()

!$OMP DO SCHEDULE(DYNAMIC)
do iy = 1, BSE%ipf_ht
    do ix = 1, BSE%ipf_wd
        s = 0.0
        icnt = BSE%ipf_wd*(iy-1) + ix
! get the orientation and determine the quaternion to be applied to all the 
! detector vectors
        eu = e_T( edinp = dble( Eangles(1:3, icnt) ) )
        qu = eu%eq()
        quat = Quaternion_T( qd = qu%q_copyd() )

! then perform the usual interpolation from the master pattern
        do jd = 1, BSE%numdet
! get the pixel direction cosines from the pre-computed array
            ddc = dble( (/ BSE%rgx(jd),BSE%rgy(jd),BSE%rgz(jd) /) )
! apply the beam tilt correction to this direction cosine
            dquat = Quaternion_T( qd = dble(BSE%beamtiltq(1:4, ix, iy)) )
            ddc = dquat%quat_Lp( ddc )
! apply the grain rotation to the detector direction cosines
            ddc = quat%quat_Lp( ddc )
            dc = sngl(ddc/sqrt(sum(ddc*ddc)))
! convert these direction cosines to interpolation coordinates in the Rosca-Lambert projection
            dc = dc / sqrt(sum(dc*dc))
            call LambertgetInterpolation(dc, scl, mpnl%npx, mpnl%npx, nix, niy, nixp, niyp, dx, dy, dxm, dym)

            if (dc(3) .ge. 0.0) then
                do kk = Emin, Emax
                    s = s +  BSE%corfactor(jd) * ( BSE%mLPNH(nix,niy,kk) * dxm * dym + &
                                                   BSE%mLPNH(nixp,niy,kk) * dx * dym + &
                                                   BSE%mLPNH(nix,niyp,kk) * dxm * dy + &
                                                   BSE%mLPNH(nixp,niyp,kk) * dx * dy )
                end do
              else
                do kk = Emin, Emax
                    s = s +  BSE%corfactor(jd) * ( BSE%mLPSH(nix,niy,kk) * dxm * dym + &
                                                   BSE%mLPSH(nixp,niy,kk) * dx * dym + &
                                                   BSE%mLPSH(nix,niyp,kk) * dxm * dy + &
                                                   BSE%mLPSH(nixp,niyp,kk) * dx * dy )
                end do
            end if
        end do
        BSEimage(ix,iy) = s 
    end do 
    if (mod(iy,10).eq.0) write (*,*) ' working on line ', iy
end do 
!$OMP END DO
!$OMP END PARALLEL

end associate 

end subroutine ComputeBSEimage_

!--------------------------------------------------------------------------
subroutine ComputeBSErings_(self, mcnl, mpnl, numang, Eangles, Emin, Emax, BSEimage, SO2)
!DEC$ ATTRIBUTES DLLEXPORT :: ComputeBSErings_
!! author: MDG 
!! version: 1.0 
!! date: 12/12/22
!!
!! compute an energy-weighted BSE image split into Lambert rings for a given ROI

use mod_MCfiles
use mod_MPfiles
use mod_io
use mod_so2
use mod_quaternions 
use mod_rotations
use mod_Lambert
use mod_image 
use omp_lib

IMPLICIT NONE

class(BSE_T), INTENT(INOUT)                 :: self
type(MCOpenCLNameListType),INTENT(IN)       :: mcnl
type(EBSDMasterNameListType),INTENT(IN)     :: mpnl
integer(kind=irg),INTENT(IN)                :: numang 
real(kind=sgl),INTENT(IN)                   :: Eangles(3,numang)
integer(kind=irg),INTENT(IN)                :: Emin
integer(kind=irg),INTENT(IN)                :: Emax 
real(kind=sgl),INTENT(OUT),allocatable      :: BSEimage(:,:,:)
type(SO2_T),INTENT(IN)                      :: SO2 

type(IO_T)                                  :: Message
type(Quaternion_T)                          :: quat, dquat 
type(e_T)                                   :: eu
type(q_T)                                   :: qu
type(Lambert_T)                             :: L

real(kind=sgl)                              :: s, mi, ma, ECPfactor, q(4)
integer(kind=irg)                           :: ix, iy, icnt, jd, sz(3), nxmc
real(kind=sgl)                              :: dc(3), avdc(3), newavdc(3), ixy(2), scl, sclmc
real(kind=dbl)                              :: ddc(3)
real(kind=sgl)                              :: dx, dy, dxm, dym, x, y, z
integer(kind=irg)                           :: ii, jj, kk, istat, kd
integer(kind=irg)                           :: nix, niy, nixp, niyp, nixmc, niymc, TID
integer(kind=irg),allocatable               :: SO2ringcount(:)
integer(kind=irg),allocatable               :: SO2ringstart(:)

call setRotationPrecision('d')

associate( enl => self%nml, BSE => self%det )

scl = float(mpnl%npx) 
sz = shape(BSE%accum_e)
nxmc = (sz(2)-1)/2
sclmc = float(nxmc)

! allocate all necessary arrays
allocate(BSEimage(enl%NsqL+1,BSE%ipf_wd,BSE%ipf_ht))
BSEimage = 0.0
SO2ringstart = SO2%getSO2ringstart()

! loop over all the image pixels 
icnt = 0
call OMP_SET_NUM_THREADS(enl%nthreads)

!$OMP PARALLEL default(shared) private(ix,iy,s,icnt,qu,jd,kk,dc,nix,niy,nixp,niyp,dx,dy,dxm,dym,newavdc,ECPfactor)&
!$OMP& private(kd, nixmc, niymc, eu, quat, dquat, ddc) 

TID = OMP_GET_THREAD_NUM()

!$OMP DO SCHEDULE(DYNAMIC)
do iy = 1, BSE%ipf_ht
    do ix = 1, BSE%ipf_wd
        icnt = BSE%ipf_wd*(iy-1) + ix
! get the orientation and determine the quaternion to be applied to all the 
! detector vectors
        eu = e_T( edinp = dble( Eangles(1:3, icnt) ) )
        qu = eu%eq()
        quat = Quaternion_T( qd = qu%q_copyd() )

! then perform the usual interpolation from the master pattern, but do this ring by ring
! in the list of Lambert points...
        do jd = 0, enl%NsqL-1
          do kd = SO2ringstart(jd+1),SO2ringstart(jd+2)-1 
            s = 0.0
! get the pixel direction cosines from the pre-computed array
            ddc = dble( (/ BSE%rgx(kd+1),BSE%rgy(kd+1),BSE%rgz(kd+1) /) )
! apply the beam tilt correction to this direction cosine
            dquat = Quaternion_T( qd = dble(BSE%beamtiltq(1:4, ix, iy)) )
            ddc = dquat%quat_Lp( ddc )
! apply the grain rotation to the detector direction cosines
            ddc = quat%quat_Lp( ddc )
            dc = sngl(ddc/sqrt(sum(ddc*ddc)))
! convert these direction cosines to interpolation coordinates in the Rosca-Lambert projection
            dc = dc / sqrt(sum(dc*dc))
            call LambertgetInterpolation(dc, scl, mpnl%npx, mpnl%npx, nix, niy, nixp, niyp, dx, dy, dxm, dym)

            if (dc(3) .ge. 0.0) then
                do kk = Emin, Emax
                    s = s +  BSE%corfactor(kd+1) * ( BSE%mLPNH(nix,niy,kk) * dxm * dym + &
                                                     BSE%mLPNH(nixp,niy,kk) * dx * dym + &
                                                     BSE%mLPNH(nix,niyp,kk) * dxm * dy + &
                                                     BSE%mLPNH(nixp,niyp,kk) * dx * dy )
                end do
              else
                do kk = Emin, Emax
                    s = s +  BSE%corfactor(kd+1) * ( BSE%mLPSH(nix,niy,kk) * dxm * dym + &
                                                     BSE%mLPSH(nixp,niy,kk) * dx * dym + &
                                                     BSE%mLPSH(nix,niyp,kk) * dxm * dy + &
                                                     BSE%mLPSH(nixp,niyp,kk) * dx * dy )
                end do
            end if
          end do
          BSEimage(jd+1,ix,iy) = s 
        end do
    end do 
    if (mod(iy,10).eq.0) write (*,*) ' working on line ', iy
end do 
!$OMP END DO
!$OMP END PARALLEL

end associate 

end subroutine ComputeBSErings_

!--------------------------------------------------------------------------
subroutine BSE_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: BSE_
!! author: MDG 
!! version: 1.0 
!! date: 12/05/22
!!
!! perform the computations

use mod_EMsoft
use mod_so2
use mod_quaternions
use mod_MCfiles
use mod_MPfiles
use mod_DIfiles
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_io
use mod_rotations
use stringconstants
use mod_vendors
use mod_memory
use mod_image 
use ISO_C_BINDING
use, intrinsic :: iso_fortran_env

IMPLICIT NONE 

class(BSE_T), INTENT(INOUT)         :: self
type(EMsoft_T), INTENT(INOUT)       :: EMsoft
character(fnlen), INTENT(INOUT)     :: progname 

type(MCfile_T)                      :: MCFT
type(MPfile_T)                      :: MPFT
type(HDF_T)                         :: HDF
type(HDFnames_T)                    :: HDFnames
type(so2_T)                         :: SO2
type(IO_T)                          :: Message
type(Quaternion_T)                  :: quat
type(QuaternionArray_T)             :: qAR
type(memory_T)                      :: mem
type(DIfile_T)                      :: DIFT
type(DictionaryIndexingNameListType):: dinl
type(Vendor_T)                      :: VT

type(EBSDmasterNameListType)        :: mpnl
type(MCOpenCLNameListType)          :: mcnl

integer(kind=irg)                   :: i, sz(3), nx, hdferr, resang, resctf
integer(kind=irg)                   :: Emin, Emax, dims3(3)      ! various parameters
character(fnlen)                    :: fname, DIfile
logical                             :: refined
real(kind=sgl)                      :: scl, mi, ma
real(kind=sgl),allocatable          :: Eangles(:,:), BSEimage(:,:), BSEimagescan(:,:,:)

! declare variables for use in object oriented image module
character(fnlen)                    :: TIFF_filename, dataset
integer                             :: iostat
character(len=128)                  :: iomsg
logical                             :: isInteger
type(image_t)                       :: im
integer(int8), allocatable          :: TIFF_image(:,:)
integer                             :: dim2(2)
integer(c_int32_t)                  :: result

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()
HDFnames = HDFnames_T()

associate( enl => self%nml, BSEdetector => self%det, EBSDMCdata => MCFT%MCDT, &
           EBSDMPdata => MPFT%MPDT, DIdata => DIFT%DIDT )

! 1. read the Monte Carlo data file (HDF format)
call HDFnames%set_ProgramData(SC_MCOpenCL)
call HDFnames%set_NMLlist(SC_MCCLNameList)
call HDFnames%set_NMLfilename(SC_MCOpenCLNML)
fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfile))
call MCFT%setFileName(fname)
write (*,*) 'File name = ', trim(fname)
call MCFT%readMCfile(HDF, HDFnames, getAccume=.TRUE.)
mcnl = MCFT%getnml()

! 2. read EBSD master pattern file (HDF format)
call MPFT%setModality('EBSD')
call HDFnames%set_ProgramData(SC_EBSDmaster)
call HDFnames%set_NMLlist(SC_EBSDmasterNameList)
call HDFnames%set_NMLfilename(SC_EBSDmasterNML)
call HDFnames%set_Variable(SC_MCOpenCL)
fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfile))
call MPFT%setFileName(fname)
call MPFT%readMPfile(HDF, HDFnames, mpnl, getmLPNH=.TRUE., getmLPSH=.TRUE.)

! 3. copy some of the arrays into the BSEdetector structure and delete the originals
sz = shape(EBSDMCdata%accum_e)
nx = (sz(2)-1)/2
allocate(BSEdetector%accum_e(sz(1),-nx:nx, -nx:nx))
BSEdetector%accum_e = EBSDMCdata%accum_e
deallocate(EBSDMCdata%accum_e)

allocate(BSEdetector%mLPNH(-mpnl%npx:mpnl%npx,-mpnl%npx:mpnl%npx,EBSDMPdata%numEbins))
BSEdetector%mLPNH = EBSDMPdata%mLPNH
deallocate(EBSDMPdata%mLPNH)

allocate(BSEdetector%mLPSH(-mpnl%npx:mpnl%npx,-mpnl%npx:mpnl%npx,EBSDMPdata%numEbins))
BSEdetector%mLPSH = EBSDMPdata%mLPSH
deallocate(EBSDMPdata%mLPSH)

Emin = nint((enl%energymin - mcnl%Ehistmin)/mcnl%Ebinsize) + 1
if (Emin.lt.1)  Emin=1
if (Emin.gt.EBSDMCdata%numEbins)  Emin=EBSDMCdata%numEbins

Emax = nint((enl%energymax - mcnl%Ehistmin)/mcnl%Ebinsize) + 1
if (Emax.lt.1)  Emax=1
if (Emax.gt.EBSDMCdata%numEbins)  Emax=EBSDMCdata%numEbins

! 4. generate BSE detector direction cosine arrays
if (enl%scanmode.eq.'single') then 
  call self%GenerateBSEDetector( dble(enl%incidence), verbose=.TRUE.)
else
  call self%GenerateBSEDetector( dble(enl%incidence), verbose=.TRUE., SO2=SO2)
end if 

! 5. read the angular arrays from the HDF5 file (DI only for now)
resang = index(enl%datafile, '.ang')
resctf = index(enl%datafile, '.ctf')
DIfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(enl%datafile)
if ( (resang.eq.0).and.(resctf.eq.0) ) then 
  refined = .FALSE.
  call HDFnames%set_NMLfiles(SC_NMLfiles)
  call HDFnames%set_NMLfilename(SC_DictionaryIndexingNML)
  call HDFnames%set_NMLparameters(SC_NMLparameters)
  call HDFnames%set_NMLlist(SC_DictionaryIndexingNameListType)
  if (trim(enl%useangles).eq.'refined') then
      call DIFT%readDotProductFile(EMsoft, HDF, HDFnames, DIfile, hdferr, &
                                   getRefinedEulerAngles=.TRUE.)
      refined = .TRUE.
  else
      call DIFT%readDotProductFile(EMsoft, HDF, HDFnames, DIfile, hdferr, &
                                   getPhi1=.TRUE., &
                                   getPhi=.TRUE., &
                                   getPhi2=.TRUE.)
  end if

  dinl = DIFT%getNameList()

  if (sum(dinl%ROI).eq.0) then 
      BSEdetector%ipf_wd = dinl%ipf_wd
      BSEdetector%ipf_ht = dinl%ipf_ht
  else 
      BSEdetector%ipf_wd = dinl%ROI(3)
      BSEdetector%ipf_ht = dinl%ROI(4)
  end if

  nx = BSEdetector%ipf_wd * BSEdetector%ipf_ht
  allocate(Eangles(3, nx))
  if (trim(enl%useangles).eq.'original') then
      do i=1,nx 
        Eangles(1:3,i) = (/ DIdata%Phi1(i), DIdata%Phi(i), DIdata%Phi2(i) /)
      end do
      deallocate(DIdata%Phi1, DIdata%Phi, DIdata%Phi2)
  else 
      Eangles = DIdata%RefinedEulerAngles
      deallocate(DIdata%RefinedEulerAngles)
  end if

  if (maxval(Eangles).gt.(2.D0*cPi)) Eangles = Eangles * dtor
else
  if (resang.ne.0) then ! we have an .ang file
    VT = Vendor_T()
    call VT%getAnglesfromANGfile(DIfile, dinl%ipf_wd, dinl%ipf_ht, dinl%StepX, dinl%StepY, Eangles)
  else  ! we must have a .ctf file
    VT = Vendor_T()
    call VT%getAnglesfromCTFfile(DIfile, dinl%ipf_wd, dinl%ipf_ht, dinl%StepX, dinl%StepY, Eangles)
  end if 
  BSEdetector%ipf_wd = dinl%ipf_wd
  BSEdetector%ipf_ht = dinl%ipf_ht
  dinl%ROI = (/ 0, 0, 0, 0 /)
end if 

! 6. generate the beam tilt quaternions for the ROI 
call self%GenerateBSEbeamtiltquaternions(dinl, verbose=.TRUE.)

! 7. save the detector pixel unit vectors for debugging purposes... this is in PoVray format
! to visualize the hemispherical ring of directions that fall onto the detector
! This part should be commented out once everything works properly
! open(dataunit,file='BSEpoints2.txt',status='unknown',form='formatted')
! write (dataunit,"(A)") '#declare gridpoints = '
! write (dataunit,"(A)") 'union{'
! do i=1,BSEdetector%numdet
!   write (dataunit,"('sphere{<',2(F10.6,','),F10.6,'>,0.01}')") BSEdetector%rgx(i), BSEdetector%rgy(i), BSEdetector%rgz(i) 
! end do 
! write (dataunit,"(A)") '};'

! scl = enl%workingdistance

! write (dataunit,"(A)") '#declare pixelrods= '
! write (dataunit,"(A)") 'union{'
! do i=1,BSEdetector%numdet
!   write (dataunit,"('cylinder{<0.0,0.0,0.0><',2(F10.6,','),F10.6,'>,0.005}')") &
!                      BSEdetector%rgx(i)/BSEdetector%rgz(i)*scl, BSEdetector%rgy(i)/BSEdetector%rgz(i)*scl, scl 
! end do 
! write (dataunit,"(A)") '};'

! close(dataunit,status='keep')

! 8. and finally perform the image computations
if (self%nml%scanmode.eq.'single') then 
  call self%ComputeBSEimage(mcnl, mpnl, nx, Eangles, Emin, Emax, BSEimage)

! and save the resulting BSE image to a tiff file
! output the ADP map as a tiff file 
  fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%imagefile))
  TIFF_filename = trim(fname)

! allocate memory for image
  allocate(TIFF_image(BSEdetector%ipf_wd,BSEdetector%ipf_ht))

! fill the image with whatever data you have (between 0 and 255)
  ma = maxval(BSEimage)
  mi = minval(BSEimage)

  TIFF_image = int(255 * (BSEimage-mi)/(ma-mi))

! set up the image_t structure
  im = image_t(TIFF_image)
  if(im%empty()) call Message%printMessage("EMBSE","failed to convert array to image")

! create the file
  call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
  if(0.ne.iostat) then
    call Message%printMessage("failed to write image to file : "//iomsg)
  else  
    call Message%printMessage('BSE image written to '//trim(TIFF_filename))
  end if 
  deallocate(TIFF_image)
else ! scanmode = 'scan' 
! this requires a different intensity computation, namely the intensity for each 
! square Lambert ring; since each ring corresponds to a particular theta value, storing
! the intensities by ring provides an easy way to study the angular dependence of 
! the BSE intensity for a given microstructure.  We'll store these in an HDF5 file
! so that we can do some post-processing in IDL or Matlab.
  call self%ComputeBSErings(mcnl, mpnl, nx, Eangles, Emin, Emax, BSEimagescan, SO2)

! and dump this into a simple hdf5 file called rings.h5
  fname = 'rings.h5'
  hdferr =  HDF%createFile(fname)
  dims3 = shape(BSEimagescan)
  write (*,*) 'shape(BSEimagescan) = ',dims3
  dataset = 'rings'
  hdferr = HDF%writeDatasetFloatArray(dataset, BSEimagescan, dims3(1), dims3(2), dims3(3) )
  call HDF%pop(.TRUE.)
end if 

end associate

! open the HDF interface
call closeFortranHDFInterface()

end subroutine BSE_



end module mod_BSE