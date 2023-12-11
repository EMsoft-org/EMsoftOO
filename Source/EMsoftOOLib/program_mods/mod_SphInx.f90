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

module mod_SphInx
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/08/23
  !!
  !! class definition for the EMSphInx program

use mod_kinds
use mod_global
use mod_SphInxSupport

IMPLICIT NONE 

! class definition
type, public :: SphInx_T
private 
  character(fnlen)          :: nmldeffile = 'EMSphInx.nml'
  type(SphInxNameListType)  :: nml 
  real(kind=sgl)            :: sig
  real(kind=sgl)            :: xpc
  real(kind=sgl)            :: ypc
  real(kind=sgl)            :: L

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: SphInx_
  procedure, pass(self) :: setbw_
  procedure, pass(self) :: getbw_
  procedure, pass(self) :: setnormed_
  procedure, pass(self) :: getnormed_
  procedure, pass(self) :: setrefine_
  procedure, pass(self) :: getrefine_
  procedure, pass(self) :: setflipy_
  procedure, pass(self) :: getflipy_
  procedure, pass(self) :: setROImask_
  procedure, pass(self) :: getROImask_
  procedure, pass(self) :: setROIfile_
  procedure, pass(self) :: getROIfile_
  procedure, pass(self) :: setnregions_
  procedure, pass(self) :: getnregions_
  procedure, pass(self) :: setnthread_
  procedure, pass(self) :: getnthread_
  procedure, pass(self) :: setbatchsize_
  procedure, pass(self) :: getbatchsize_
  procedure, pass(self) :: setscandims_
  procedure, pass(self) :: getscandims_
  procedure, pass(self) :: setpatdims_
  procedure, pass(self) :: getpatdims_
  procedure, pass(self) :: setdelta_
  procedure, pass(self) :: getdelta_
  procedure, pass(self) :: setpctr_
  procedure, pass(self) :: getpctr_
  procedure, pass(self) :: setvendor_
  procedure, pass(self) :: getvendor_
  procedure, pass(self) :: setthetac_
  procedure, pass(self) :: getthetac_
  procedure, pass(self) :: setbinning_
  procedure, pass(self) :: getbinning_
  procedure, pass(self) :: setcircmask_
  procedure, pass(self) :: getcircmask_
  procedure, pass(self) :: setmasterfile_
  procedure, pass(self) :: getmasterfile_
  procedure, pass(self) :: setpatfile_
  procedure, pass(self) :: getpatfile_
  procedure, pass(self) :: setHDFstrings_
  procedure, pass(self) :: getHDFstrings_
  procedure, pass(self) :: setinputtype_
  procedure, pass(self) :: getinputtype_
  procedure, pass(self) :: setdatafile_
  procedure, pass(self) :: getdatafile_
  procedure, pass(self) :: setctffile_
  procedure, pass(self) :: getctffile_
  procedure, pass(self) :: setangfile_
  procedure, pass(self) :: getangfile_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: SphInx => SphInx_
  generic, public :: setbw => setbw_
  generic, public :: getbw => getbw_
  generic, public :: setnormed => setnormed_
  generic, public :: getnormed => getnormed_
  generic, public :: setrefine => setrefine_
  generic, public :: getrefine => getrefine_
  generic, public :: setflipy => setflipy_
  generic, public :: getflipy => getflipy_
  generic, public :: setROImask => setROImask_
  generic, public :: getROImask => getROImask_
  generic, public :: setROIfile => setROIfile_
  generic, public :: getROIfile => getROIfile_
  generic, public :: setnregions => setnregions_
  generic, public :: getnregions => getnregions_
  generic, public :: setnthread => setnthread_
  generic, public :: getnthread => getnthread_
  generic, public :: setbatchsize => setbatchsize_
  generic, public :: getbatchsize => getbatchsize_
  generic, public :: setscandims => setscandims_
  generic, public :: getscandims => getscandims_
  generic, public :: setpatdims => setpatdims_
  generic, public :: getpatdims => getpatdims_
  generic, public :: setdelta => setdelta_
  generic, public :: getdelta => getdelta_
  generic, public :: setpctr => setpctr_
  generic, public :: getpctr => getpctr_
  generic, public :: setvendor => setvendor_
  generic, public :: getvendor => getvendor_
  generic, public :: setthetac => setthetac_
  generic, public :: getthetac => getthetac_
  generic, public :: setbinning => setbinning_
  generic, public :: getbinning => getbinning_
  generic, public :: setcircmask => setcircmask_
  generic, public :: getcircmask => getcircmask_
  generic, public :: setmasterfile => setmasterfile_
  generic, public :: getmasterfile => getmasterfile_
  generic, public :: setpatfile => setpatfile_
  generic, public :: getpatfile => getpatfile_
  generic, public :: setHDFstrings => setHDFstrings_
  generic, public :: getHDFstrings => getHDFstrings_
  generic, public :: setinputtype => setinputtype_
  generic, public :: getinputtype => getinputtype_
  generic, public :: setdatafile => setdatafile_
  generic, public :: getdatafile => getdatafile_
  generic, public :: setctffile => setctffile_
  generic, public :: getctffile => getctffile_
  generic, public :: setangfile => setangfile_
  generic, public :: getangfile => getangfile_

end type SphInx_T

! the constructor routine for this class 
interface SphInx_T
  module procedure SphInx_constructor
end interface SphInx_T

contains

!--------------------------------------------------------------------------
type(SphInx_T) function SphInx_constructor( nmlfile ) result(SphInx)
!! author: MDG 
!! version: 1.0 
!! date: 12/08/23
!!
!! constructor for the SphInx_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call SphInx%readNameList(nmlfile)

end function SphInx_constructor

!--------------------------------------------------------------------------
subroutine SphInx_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 12/08/23
!!
!! destructor for the SphInx_T Class
 
IMPLICIT NONE

type(SphInx_T), INTENT(INOUT)  :: self 

call reportDestructor('SphInx_T')

end subroutine SphInx_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 12/08/23
!!
!! read the namelist from an nml file for the SphInx_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(SphInx_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg)       :: bw
logical                 :: normed 
logical                 :: refine
logical                 :: flipy
integer(kind=irg)       :: ROImask(4)
character(fnlen)        :: ROIfile
integer(kind=irg)       :: nregions
integer(kind=irg)       :: nthread
integer(kind=irg)       :: batchsize
real(kind=sgl)          :: scandims(4)
integer(kind=sgl)       :: patdims(2)
real(kind=sgl)          :: delta
real(kind=sgl)          :: pctr(3)
character(fnlen)        :: vendor
real(kind=sgl)          :: thetac
integer(kind=irg)       :: binning
logical                 :: circmask
character(fnlen)        :: masterfile
character(fnlen)        :: patfile
character(fnlen)        :: HDFstrings(10)
character(fnlen)        :: inputtype
character(fnlen)        :: datafile
character(fnlen)        :: ctffile
character(fnlen)        :: angfile

namelist /SphInxNameList/ bw, normed, refine, ROImask, ROIfile, nregions, nthread, batchsize, thetac, delta, &
                          patdims, pctr, scandims, binning, circmask, masterfile, vendor, &
                          patfile, HDFstrings, inputtype, datafile, ctffile, angfile, flipy

! default values
bw = 88
normed = .TRUE.
refine = .TRUE.
flipy = .FALSE.
ROImask = (/ 0, 0, 0, 0 /)
ROIfile = 'undefined'
nregions = 10
nthread = 0
batchsize = 0
scandims = (/ 0.0, 0.0, 0.1, 0.1 /)
patdims = (/ 640, 480 /)
delta = 55.0
pctr = (/ 0.0, 0.0, 15000.0 /)
vendor = 'EMsoft' 
thetac = 10.0
binning = 1
circmask = .FALSE.
masterfile = 'undefined'
patfile = 'undefined'
inputtype = 'Binary'
HDFstrings = (/ '', '', '', '', '', '', '', '', '', '' /)
datafile = 'undefined'
ctffile = 'undefined' 
angfile = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=SphInxNameList)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(masterfile).eq.'undefined') then
  call Message%printError('readNameList:',' masterfile file name is undefined in '//nmlfile)
 end if

 if (trim(patfile).eq.'undefined') then
  call Message%printError('readNameList:',' patfile file name is undefined in '//nmlfile)
 end if

 if (trim(datafile).eq.'undefined') then
  call Message%printError('readNameList:',' datafile file name is undefined in '//nmlfile)
 end if

 if ((trim(ctffile).eq.'undefined').AND.(trim(angfile).eq.'undefined')) then 
  call Message%printError('readNameList:',' either ctffile or angfile must be defined in '//nmlfile)
 end if 
end if

! and set the actual values
self%nml%bw = bw
self%nml%normed = normed 
self%nml%refine = refine 
self%nml%flipy = flipy 
self%nml%ROImask = ROImask
self%nml%ROIfile = ROIfile
self%nml%nregions = nregions
self%nml%nthread = nthread
self%nml%batchsize = batchsize
self%nml%scandims = scandims
self%nml%patdims = patdims
self%nml%delta = delta
self%nml%pctr = pctr
self%nml%vendor = vendor
self%nml%thetac = thetac
self%nml%binning = binning
self%nml%circmask = circmask
self%nml%masterfile = masterfile
self%nml%patfile = patfile
self%nml%inputtype = inputtype
self%nml%HDFstrings = HDFstrings
self%nml%datafile = datafile
self%nml%ctffile = ctffile
self%nml%angfile = angfile

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 12/08/23
!!
!! pass the namelist for the SphInx_T Class to the calling program

IMPLICIT NONE 

class(SphInx_T), INTENT(INOUT)          :: self
type(SphInxNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine setbw_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setbw_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set bw in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%bw = inp

end subroutine setbw_

!--------------------------------------------------------------------------
function getbw_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getbw_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get bw from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%bw

end function getbw_

!--------------------------------------------------------------------------
subroutine setnormed_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnormed_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set normed in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
logical, INTENT(IN)       :: inp

self%nml%normed = inp

end subroutine setnormed_

!--------------------------------------------------------------------------
function getnormed_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnormed_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get normed from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
logical                   :: out

out = self%nml%normed

end function getnormed_

!--------------------------------------------------------------------------
subroutine setrefine_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setrefine_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set refine in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
logical, INTENT(IN)       :: inp

self%nml%refine = inp

end subroutine setrefine_

!--------------------------------------------------------------------------
function getrefine_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getrefine_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get refine from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
logical                   :: out

out = self%nml%refine

end function getrefine_

!--------------------------------------------------------------------------
subroutine setflipy_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setflipy_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set flipy in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
logical, INTENT(IN)       :: inp

self%nml%flipy = inp

end subroutine setflipy_

!--------------------------------------------------------------------------
function getflipy_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getflipy_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get flipy from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
logical                   :: out

out = self%nml%flipy

end function getflipy_

!--------------------------------------------------------------------------
subroutine setROImask_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setROImask_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set ROImask in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp(4)

self%nml%ROImask = inp

end subroutine setROImask_

!--------------------------------------------------------------------------
function getROImask_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getROImask_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get ROImask from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out(4)

out = self%nml%ROImask

end function getROImask_

!--------------------------------------------------------------------------
subroutine setROIfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setROIfile_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set ROIfile in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%ROIfile = trim(inp)

end subroutine setROIfile_

!--------------------------------------------------------------------------
function getROIfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getROIfile_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get ROIfile from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%ROIfile)

end function getROIfile_

!--------------------------------------------------------------------------
subroutine setnregions_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnregions_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set nregions in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%nregions = inp

end subroutine setnregions_

!--------------------------------------------------------------------------
function getnregions_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnregions_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get nregions from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%nregions

end function getnregions_

!--------------------------------------------------------------------------
subroutine setnthread_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnthread_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set nthread in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%nthread = inp

end subroutine setnthread_

!--------------------------------------------------------------------------
function getnthread_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnthread_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get nthread from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%nthread

end function getnthread_

!--------------------------------------------------------------------------
subroutine setbatchsize_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setbatchsize_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set batchsize in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%batchsize = inp

end subroutine setbatchsize_

!--------------------------------------------------------------------------
function getbatchsize_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getbatchsize_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get batchsize from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%batchsize

end function getbatchsize_

!--------------------------------------------------------------------------
subroutine setscandims_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setscandims_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set scandims in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp(4)

self%nml%scandims = inp

end subroutine setscandims_

!--------------------------------------------------------------------------
function getscandims_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getscandims_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get scandims from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out(4)

out = self%nml%scandims

end function getscandims_

!--------------------------------------------------------------------------
subroutine setpatdims_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setpatdims_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set patdims in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
integer(kind=sgl), INTENT(IN)      :: inp(2)

self%nml%patdims = inp

end subroutine setpatdims_

!--------------------------------------------------------------------------
function getpatdims_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getpatdims_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get patdims from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
integer(kind=sgl)                  :: out(2)

out = self%nml%patdims

end function getpatdims_

!--------------------------------------------------------------------------
subroutine setdelta_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdelta_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set delta in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%delta = inp

end subroutine setdelta_

!--------------------------------------------------------------------------
function getdelta_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdelta_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get delta from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%delta

end function getdelta_

!--------------------------------------------------------------------------
subroutine setpctr_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setpctr_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set pctr in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp(3)

self%nml%pctr = inp

end subroutine setpctr_

!--------------------------------------------------------------------------
function getpctr_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getpctr_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get pctr from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out(3)

out = self%nml%pctr

end function getpctr_

!--------------------------------------------------------------------------
subroutine setvendor_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setvendor_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set vendor in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%vendor = trim(inp)

end subroutine setvendor_

!--------------------------------------------------------------------------
function getvendor_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getvendor_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get vendor from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%vendor)

end function getvendor_

!--------------------------------------------------------------------------
subroutine setthetac_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setthetac_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set thetac in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%thetac = inp

end subroutine setthetac_

!--------------------------------------------------------------------------
function getthetac_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getthetac_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get thetac from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%thetac

end function getthetac_

!--------------------------------------------------------------------------
subroutine setbinning_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setbinning_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set binning in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%binning = inp

end subroutine setbinning_

!--------------------------------------------------------------------------
function getbinning_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getbinning_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get binning from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%binning

end function getbinning_

!--------------------------------------------------------------------------
subroutine setcircmask_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setcircmask_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set circmask in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
logical, INTENT(IN)       :: inp

self%nml%circmask = inp

end subroutine setcircmask_

!--------------------------------------------------------------------------
function getcircmask_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getcircmask_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get circmask from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
logical                   :: out

out = self%nml%circmask

end function getcircmask_

!--------------------------------------------------------------------------
subroutine setmasterfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setmasterfile_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set masterfile in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%masterfile = trim(inp)

end subroutine setmasterfile_

!--------------------------------------------------------------------------
function getmasterfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getmasterfile_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get masterfile from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%masterfile)

end function getmasterfile_

!--------------------------------------------------------------------------
subroutine setpatfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setpatfile_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set patfile in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%patfile = trim(inp)

end subroutine setpatfile_

!--------------------------------------------------------------------------
function getpatfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getpatfile_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get patfile from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%patfile)

end function getpatfile_

!--------------------------------------------------------------------------
subroutine setHDFstrings_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setHDFstrings_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set HDFstrings in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp(10)

integer(kind=irg)                  :: i 

do i=1,10
  self%nml%HDFstrings(i) = trim(inp(i))
end do

end subroutine setHDFstrings_

!--------------------------------------------------------------------------
function getHDFstrings_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getHDFstrings_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get HDFstrings from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out(10)

integer(kind=irg)                  :: i 

do i=1,10 
  out(i) = trim(self%nml%HDFstrings(i))
end do

end function getHDFstrings_

!--------------------------------------------------------------------------
subroutine setinputtype_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setinputtype_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set inputtype in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%inputtype = trim(inp)

end subroutine setinputtype_

!--------------------------------------------------------------------------
function getinputtype_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getinputtype_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get inputtype from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%inputtype)

end function getinputtype_

!--------------------------------------------------------------------------
subroutine setdatafile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdatafile_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set datafile in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%datafile = trim(inp)

end subroutine setdatafile_

!--------------------------------------------------------------------------
function getdatafile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdatafile_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get datafile from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%datafile)

end function getdatafile_

!--------------------------------------------------------------------------
subroutine setctffile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setctffile_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set ctffile in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%ctffile = trim(inp)

end subroutine setctffile_

!--------------------------------------------------------------------------
function getctffile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getctffile_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get ctffile from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%ctffile)

end function getctffile_

!--------------------------------------------------------------------------
subroutine setangfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setangfile_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! set angfile in the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%angfile = trim(inp)

end subroutine setangfile_

!--------------------------------------------------------------------------
function getangfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getangfile_
!! author: MDG
!! version: 1.0
!! date: 12/08/23
!!
!! get angfile from the SphInx_T class

IMPLICIT NONE

class(SphInx_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%angfile)

end function getangfile_

!--------------------------------------------------------------------------
subroutine SphInx_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: SphInx_
!! author: WL/MDG 
!! version: 1.0 
!! date: 12/08/23
!!
!! perform the computations

use mod_EMsoft
use mod_crystallography
use mod_symmetry
use mod_so3
use mod_quaternions
use mod_fft_wrap
use mod_fftw3
use mod_MCfiles
use mod_MPfiles
use mod_DIfiles
use mod_DIsupport
use mod_filters
use mod_patterns
use omp_lib
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_io
use mod_rotations
use stringconstants
use mod_memory
use mod_timing
use mod_SphereIndexer
use mod_vendors

IMPLICIT NONE 

class(SphInx_T), INTENT(INOUT)      :: self
type(EMsoft_T), INTENT(INOUT)       :: EMsoft
character(fnlen), INTENT(INOUT)     :: progname 
type(HDFnames_T), INTENT(INOUT)     :: HDFnames

type(MCfile_T)                      :: MCFT
type(MPfile_T)                      :: MPFT
type(DIfile_T)                      :: DIFT
type(HDF_T)                         :: HDF
type(HDFnames_T)                    :: saveHDFnames
type(EBSDmasterNameListType)        :: mpnl
type(MCOpenCLNameListType)          :: mcnl
type(SphereIndexer)                 :: threadindexer 
type(SpaceGroup_T)                  :: SG
type(Cell_T)                        :: cell
type(Timing_T)                      :: timer
type(IO_T)                          :: Message
type(Vendor_T)                      :: VT

real(kind=dbl),allocatable          :: mLPNH(:,:), mLPSH(:,:), weights(:)
integer(kind=irg)                   :: L, totnumexpt, hdferr, TID, NUMTHREADS, io_int(3), rcnt, ipos
integer(kind=irg)                   :: dims(2), dim2(2), dim3(3)
character(11)                       :: dstr 
character(15)                       :: tstrb, tstre
character(3)                        :: vendor
integer(kind=irg)                   :: i,j, std=6, kk, d, itype, n
integer(kind=irg)                   :: pgnum
! integer  (kind=irg                ),allocatable   :: indexmain(:,:)
real(kind=sgl)                      :: tstart, tstop, io_real(2), ma, mi, pcvals(3)
character(fnlen)                    :: xtalname, fname
integer(kind=irg)                   :: istart, iend, jstart, jend, ipf_wd, ipf_ht
logical                             :: ROIselected 
integer(kind=irg)                   :: ipar(4)
character(fnlen),allocatable        :: MessageLines(:)
integer(kind=irg)                   :: NumLines, numsx, numsy, recordsize, correctsize
character(fnlen)                    :: TitleMessage, exectime, wisdomFile
character(100)                      :: c
real(kind=dbl)                      :: w, Jres, fpar(5)

real(kind=dbl),allocatable          :: quats(:,:), xcorr(:), pat(:,:), ksqarray(:,:)
real(kind=sgl),allocatable          :: pat32(:), patrow(:), exptIQ(:), EBSDpat(:,:)
type(IdxRes)                        :: ires
integer(kind=irg)                   :: istat, iunitexpt = 41
integer(HSIZE_T)                    :: dims3(3), offset3(3)
type(C_PTR)                         :: planf

call setRotationPrecision('d')

call MPFT%setModality('EBSD')

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()

associate( sinl => self%nml, EBSDMCdata => MCFT%MCDT, EBSDMPdata => MPFT%MPDT )

! 1. read the Monte Carlo data file (HDF format)
call HDFnames%set_ProgramData(SC_MCOpenCL)
call HDFnames%set_NMLlist(SC_MCCLNameList)
call HDFnames%set_NMLfilename(SC_MCOpenCLNML)
fname = EMsoft%generateFilePath('EMdatapathname',trim(sinl%masterfile))
call MCFT%setFileName(fname)
call MCFT%readMCfile(HDF, HDFnames, getAccume=.TRUE.)
mcnl = MCFT%getnml()
xtalname = trim(mcnl%xtalname)
self%sig = mcnl%sig

! 2. read EBSD master pattern file (HDF format)
call HDFnames%set_ProgramData(SC_EBSDmaster)
call HDFnames%set_NMLlist(SC_EBSDmasterNameList)
call HDFnames%set_NMLfilename(SC_EBSDmasterNML)
call HDFnames%set_Variable(SC_MCOpenCL)

fname = EMsoft%generateFilePath('EMdatapathname',trim(sinl%masterfile))
call MPFT%setFileName(fname)
call MPFT%readMPfile(HDF, HDFnames, mpnl, getmLPNH=.TRUE., getmLPSH=.TRUE.)

! 3. build the energy averaged master pattern
! first make sure that we do some appropriate energy weighting to make the master patterns 2D instead of 3D
n = size(EBSDMCdata%accum_e, 1) ! get number of energy bins
allocate(weights(n)) ! allocate space for energy histogram
do i = 1, n
  weights(i) = sum(EBSDMCdata%accum_e(i,:,:)) ! this could be modified to sum over partial rectangle
enddo
weights = weights / sum(weights) ! this is currently weighted over the full square Lambert

! build energy weighted master pattern
d = mpnl%npx
allocate(mLPNH(-d:d,-d:d))
allocate(mLPSH(-d:d,-d:d))
mLPNH = 0.D0
mLPSH = 0.D0
do i = 1, n
  mLPNH = mLPNH + EBSDMPdata%mLPNH(:,:,i) * weights(i)
  mLPSH = mLPSH + EBSDMPdata%mLPSH(:,:,i) * weights(i)
enddo

! set the timer
timer = Timing_T()
dstr = timer%getDateString()
tstrb = timer%getTimeString()
tstre = ''

! convert some parameters from the name list to the regular EMsoft variables
numsx = sinl%patdims(1)
numsy = sinl%patdims(2)
ipf_wd = int(sinl%scandims(1)) 
ipf_ht = int(sinl%scandims(2))

pcvals = getEMsoftPCcoordinates(sinl%pctr, sinl%vendor, sinl%delta, numsx, numsy)
self%xpc = pcvals(1)
self%ypc = pcvals(2)
self%L   = pcvals(3)

! do we need a region of interest (ROI) ?
if (sum(sinl%ROImask).ne.0) then
  ROIselected = .TRUE.
  istart = sinl%ROImask(1)
  jstart = sinl%ROImask(2)
  iend   = istart+sinl%ROImask(3)-1
  jend   = jstart+sinl%ROImask(4)-1
else
  ROIselected = .FALSE.
  istart = 1
  jstart = 1
  iend   = ipf_wd
  jend   = ipf_ht
end if
if (ROIselected.eqv..TRUE.) then 
    totnumexpt = sinl%ROImask(3)*sinl%ROImask(4)
else
    totnumexpt = ipf_wd*ipf_ht
end if
L = numsx*numsy/sinl%binning**2

! make sure that correctsize is a multiple of 16; if not, make it so
if (mod(L,16) .ne. 0) then
    correctsize = 16*ceiling(float(L)/16.0)
else
    correctsize = L
end if

recordsize = correctsize*4

! we know that the master pattern file exists, and it also has all the
! crystallographic data in it, so we read that here instead of assuming
! that the actual .xtal file exists on this system ...
call cell%setFileName(xtalname)
call cell%readDataHDF(SG, EMsoft, useXtalName=fname)
! extract the point group number
pgnum = SG%getPGnumber()
io_int = pgnum
call Message%WriteValue(' Setting point group number to ',io_int,1)

!=====================================================
! open the pattern file for reading
!=====================================================
! set the vendor inputtype for the pattern file
VT = Vendor_T( sinl%inputtype )
itype = VT%get_itype()
call VT%set_filename(sinl%patfile)

!===================================================================================
! open the file with experimental patterns; depending on the inputtype parameter, this
! can be a regular binary file, as produced by a MatLab or IDL script (default); a
! pattern file produced by EMEBSD.f90 etc.; or a vendor binary or HDF5 file... in each case we need to
! open the file and leave it open, then use the getExpPatternRow() routine to read a row
! of patterns into the exppatarray variable ...  at the end, we use closeExpPatternFile() to
! properly close the experimental pattern file
if ( (itype.eq.4) .or. (itype.eq.6) .or. (itype.eq.7) .or. (itype.eq.8) ) then
  istat = VT%openExpPatternFile(EMsoft, ipf_wd, L, recordsize, sinl%HDFstrings, HDF)
else
  istat = VT%openExpPatternFile(EMsoft, ipf_wd, L, recordsize)
end if

if (sinl%flipy.eqv..TRUE.) then 
  call Message%printMessage(' EBSD patterns will be flipped vertically ')
else
  call Message%printMessage(' EBSD patterns will not be flipped vertically ')
end if

dims3 = (/ numsx, numsy, ipf_wd /)
offset3 = (/ 0, 0, 0 /)

!=====================================================
!=====================================================
!=====================================================
! indexing loop goes here.
allocate(quats (4,totnumexpt))
allocate(xcorr (  totnumexpt))

if (sinl%nthread.eq.0) then 
  call OMP_SET_NUM_THREADS(OMP_GET_MAX_THREADS())
else
  call OMP_SET_NUM_THREADS(sinl%nthread)
end if

!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(TID, threadindexer, i, j, offset3, patrow, pat32, mi, ma, pat, ires, rcnt, ipos)

NUMTHREADS = OMP_GET_NUM_THREADS()
TID = OMP_GET_THREAD_NUM()
!$OMP CRITICAL  
if (TID.eq.0) then 
  io_int(1) = NUMTHREADS
  call Message%WriteValue('  --> number of threads allocated : ',io_int,1)
  rcnt = 1
end if
!$OMP END CRITICAL  

allocate(pat   ( numsx ,numsy ))
allocate(patrow(    L * ipf_wd))
allocate(pat32 (             L))

!$OMP BARRIER
! this should be possible in a parallel fashion if we remove accessing the EMsoftConfig.json and fftw wisdom files 
! from the DiscreteSHTConstants_Init routine...
! 
if (TID.eq.0) then 
  io_int(1) = TID
  call Message%WriteValue(' loading wisdom for thread ',io_int,1)
  call FFTWisdom%load(EMsoft)
  call Message%printMessage(' initializing all other threads ... ',"(A/)")
end if
!$OMP BARRIER
!$OMP CRITICAL  
  call threadindexer%init(sinl%bw, dble(self%sig), dble(self%L), dble(sinl%thetac), dble(sinl%delta), numsx, numsy, &
                    mpnl%npx, mLPNH, mLPSH, sinl%circmask)
!$OMP END CRITICAL  

!$OMP BARRIER
!$OMP DO SCHEDULE(DYNAMIC)
do j = jstart, jend ! loop over rows
  if (TID.eq.0) then 
    io_int(1:3) = (/ rcnt, minval( (/ rcnt+NUMTHREADS-1, jend /) ), jend /)
    call Message%WriteValue(" reading/indexing rows ", io_int, 3,"(I4,' through ',I4,' of',I4)")
  end if

! get a pattern row (slightly faster than reading one pattern at a time); make sure threads wait in line to acces the file
  offset3 = (/ 0, 0, (j-1) * ipf_wd /) ! get hyperslab offset for this pattern row
!$OMP CRITICAL
  call VT%getExpPatternRow(j, ipf_wd, L, L, dims3, offset3, patrow, HDFstrings=sinl%HDFstrings, HDF=HDF, flipy=sinl%flipy) 
!$OMP END CRITICAL

  do i = istart, iend ! loop over columns
    ! get a pattern from the row
    pat32 = patrow(L*(i-1)+1:L*i)
    
    ! do adaptive histogram equalization if needed
    if(sinl%nregions.gt.1) then
      ma = maxval(pat32)
      mi = minval(pat32)
      pat32 = ( (pat32 - mi) / (ma - mi) ) * 255.0 ! rescale from 0-255
      pat = adhisteq(sinl%nregions, numsx, numsy, nint(pat32)) ! do AHE
    else
      pat = reshape(pat32, (/ numsx, numsy /) ) ! convert from vectorized to 2d (and float to double)
    endif

    ! index the pattern
    ires = threadindexer%index(pat, dble(self%xpc), dble(self%ypc), sinl%refine) ! for now just use a single pattern center
    if (ROIselected.eqv..TRUE.) then 
      ipos = i + (j-1) * sinl%ROImask(3)
    else
      ipos = i + (j-1) * ipf_wd
    end if 
    quats(:,ipos) = ires%qu
    xcorr(  ipos) = ires%xc
  enddo
  if (TID.eq.0) rcnt = rcnt + NUMTHREADS
enddo
!$OMP END DO
deallocate(pat32 )
deallocate(pat   )
deallocate(patrow)
!$OMP BARRIER
if (TID.eq.0) call FFTWisdom%save(EMsoft)
!$OMP END PARALLEL

!=====================================================
!=====================================================
!=====================================================

! perform some timing stuff
call CPU_TIME(tstop)
tstop = tstop - tstart
io_real(1) = float(sinl%nthread * totnumexpt) / tstop
call Message%WriteValue(' Number of experimental patterns indexed per second : ',io_real,1,"(/,F10.2,/)")

!=====================================================
!=====================================================
!=====================================================

! while the pattern file is still open, compute the pattern quality array exptIQ 
call Message%printMessage(' Computing pattern quality array ... ',"(A)",advance="no")
allocate(exptIQ(totnumexpt))

! prepare the fftw plan for this pattern size to compute pattern quality (pattern sharpness Q)
allocate(EBSDPat(numsx,numsy),stat=istat)
if (istat .ne. 0) stop '  Could not allocate arrays for EBSDPat filter'
EBSDPat = 0.0
allocate(ksqarray(numsx,numsy),stat=istat)
if (istat .ne. 0) stop '  Could not allocate ksqarray array'
Jres = 0.0
call init_getEBSDIQ(numsx, numsy, EBSDPat, ksqarray, Jres, planf)
deallocate(EBSDPat)

!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(TID, i, j, kk, offset3, EBSDpat, patrow, pat, ipos, rcnt)

NUMTHREADS = OMP_GET_NUM_THREADS()
TID = OMP_GET_THREAD_NUM()
allocate(EBSDPat(numsx, numsy))
allocate(patrow(   L * ipf_wd))

!$OMP DO SCHEDULE(DYNAMIC)
do j = jstart, jend ! loop over rows
! get a pattern row (slightly faster than reading one pattern at a time); make sure threads wait in line to acces the file
  offset3 = (/ 0, 0, (j-1) * ipf_wd /) ! get hyperslab offset for this pattern row
!$OMP CRITICAL
  call VT%getExpPatternRow(j, ipf_wd, L, L, dims3, offset3, patrow, HDFstrings=sinl%HDFstrings, HDF=HDF, flipy=sinl%flipy) 
!$OMP END CRITICAL

  do i = istart, iend ! loop over columns
    ! get a pattern from the row
    do kk=1,numsy
        EBSDPat(1:numsx,kk) = patrow((i-1)*L+(kk-1)*numsx+1:(i-1)*L+kk*numsx)
    end do
    ipos = i+(j-1)*ipf_wd
    exptIQ(ipos) = sngl(computeEBSDIQ(numsx, numsy, EBSDPat, ksqarray, Jres, planf))
  end do
enddo
!$OMP END DO
deallocate(EBSDpat, patrow)
!$OMP END PARALLEL
call Message%printMessage('   done.',"(A/)")

call VT%closeExpPatternFile()


!=====================================================
!=====================================================
!=====================================================

! ===================
! MAIN OUTPUT SECTION
! ===================

! fill the ipar array with integer parameters that are needed to write the h5ebsd file
! (anything other than what is already in the sinl structure)
ipar = 0
ipar(1) = totnumexpt
ipar(2) = pgnum
if (ROIselected.eqv..TRUE.) then
  ipar(3) = sinl%ROImask(3)
  ipar(4) = sinl%ROImask(4)
else
  ipar(3) = ipf_wd
  ipar(4) = ipf_ht
end if 

fpar = dble( (/ self%xpc, self%ypc, self%L, sinl%delta, self%sig /) )

if (sinl%datafile.ne.'undefined') then 
  vendor = 'TSL'
  fname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(sinl%datafile)
  call DIFT%setfilename(fname)
  call DIFT%sphh5ebsd_writeFile(EMsoft, HDF, HDFnames, vendor, sinl, mcnl%xtalname, dstr, tstrb, ipar, fpar, self%sig, &
                                xcorr, quats, exptIQ, progname, self%nmldeffile)
  call Message%printMessage(' Data stored in h5ebsd file : '//trim(sinl%datafile),"(A)")
end if

! we will need to also write .ang and .ctf files, if requested 
if (trim(sinl%ctffile).ne.'undefined') then 
  call VT%sphctfebsd_writeFile(EMsoft,cell,SG, sinl, mcnl%xtalname, ipar, sngl(mcnl%EkeV), sngl(mcnl%sig), xcorr, quats, exptIQ)
  call Message%printMessage(' Data stored in ctf file : '//trim(sinl%ctffile),"(/A)")
end if

if (trim(sinl%angfile).ne.'undefined') then 
    call VT%sphangebsd_writeFile(EMsoft,cell,SG, sinl,mcnl%xtalname,ipar,fpar, xcorr,quats, exptIQ)
    call Message%printMessage(' Data stored in ang file : '//trim(sinl%angfile),"(/A)")
end if


! close the fortran HDF5 interface
call closeFortranHDFInterface()


end associate

end subroutine SphInx_



end module mod_SphInx