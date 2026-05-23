! ###################################################################
! Copyright (c) 2013-2026, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_EBSDdefectIV
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/30/26
  !!
  !! class definition for the EMEBSDdefectIV program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMEBSDdefectIV program
type, public :: EBSDdefectIVNameListType
  integer(kind=irg)       :: numsx
  integer(kind=irg)       :: numsy
  integer(kind=irg)       :: binning
  integer(kind=irg)       :: nthreads
  real(kind=sgl)          :: thetac
  real(kind=sgl)          :: delta
  real(kind=sgl)          :: spotsize
  real(kind=sgl)          :: gammavalue
  real(kind=dbl)          :: beamcurrent
  real(kind=dbl)          :: dwelltime
  character(3)            :: scalingmode
  logical                 :: sampleInteractionVolume
  character(fnlen)        :: deformationfile
  character(fnlen)        :: ivolfile
  character(fnlen)        :: masterfile
  character(fnlen)        :: datafile
end type EBSDdefectIVNameListType

! class definition
type, public :: EBSDdefectIV_T
private 
  character(fnlen)                :: nmldeffile = 'EMEBSDdefectIV.nml'
  type(EBSDdefectIVNameListType)  :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: EBSDdefectIV_
  procedure, pass(self) :: setnumsx_
  procedure, pass(self) :: getnumsx_
  procedure, pass(self) :: setnumsy_
  procedure, pass(self) :: getnumsy_
  procedure, pass(self) :: setbinning_
  procedure, pass(self) :: getbinning_
  procedure, pass(self) :: setnthreads_
  procedure, pass(self) :: getnthreads_
  procedure, pass(self) :: setthetac_
  procedure, pass(self) :: getthetac_
  procedure, pass(self) :: setdelta_
  procedure, pass(self) :: getdelta_
  procedure, pass(self) :: setspotsize_
  procedure, pass(self) :: getspotsize_
  procedure, pass(self) :: setgammavalue_
  procedure, pass(self) :: getgammavalue_
  procedure, pass(self) :: setbeamcurrent_
  procedure, pass(self) :: getbeamcurrent_
  procedure, pass(self) :: setdwelltime_
  procedure, pass(self) :: getdwelltime_
  procedure, pass(self) :: setscalingmode_
  procedure, pass(self) :: getscalingmode_
  procedure, pass(self) :: setsampleInteractionVolume_
  procedure, pass(self) :: getsampleInteractionVolume_
  procedure, pass(self) :: setdeformationfile_
  procedure, pass(self) :: getdeformationfile_
  procedure, pass(self) :: setivolfile_
  procedure, pass(self) :: getivolfile_
  procedure, pass(self) :: setmasterfile_
  procedure, pass(self) :: getmasterfile_
  procedure, pass(self) :: setdatafile_
  procedure, pass(self) :: getdatafile_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: EBSDdefectIV => EBSDdefectIV_
  generic, public :: setnumsx => setnumsx_
  generic, public :: getnumsx => getnumsx_
  generic, public :: setnumsy => setnumsy_
  generic, public :: getnumsy => getnumsy_
  generic, public :: setbinning => setbinning_
  generic, public :: getbinning => getbinning_
  generic, public :: setnthreads => setnthreads_
  generic, public :: getnthreads => getnthreads_
  generic, public :: setthetac => setthetac_
  generic, public :: getthetac => getthetac_
  generic, public :: setdelta => setdelta_
  generic, public :: getdelta => getdelta_
  generic, public :: setspotsize => setspotsize_
  generic, public :: getspotsize => getspotsize_
  generic, public :: setgammavalue => setgammavalue_
  generic, public :: getgammavalue => getgammavalue_
  generic, public :: setbeamcurrent => setbeamcurrent_
  generic, public :: getbeamcurrent => getbeamcurrent_
  generic, public :: setdwelltime => setdwelltime_
  generic, public :: getdwelltime => getdwelltime_
  generic, public :: setscalingmode => setscalingmode_
  generic, public :: getscalingmode => getscalingmode_
  generic, public :: setsampleInteractionVolume => setsampleInteractionVolume_
  generic, public :: getsampleInteractionVolume => getsampleInteractionVolume_
  generic, public :: setdeformationfile => setdeformationfile_
  generic, public :: getdeformationfile => getdeformationfile_
  generic, public :: setivolfile => setivolfile_
  generic, public :: getivolfile => getivolfile_
  generic, public :: setmasterfile => setmasterfile_
  generic, public :: getmasterfile => getmasterfile_
  generic, public :: setdatafile => setdatafile_
  generic, public :: getdatafile => getdatafile_
end type EBSDdefectIV_T

! the constructor routine for this class 
interface EBSDdefectIV_T
  module procedure EBSDdefectIV_constructor
end interface EBSDdefectIV_T

contains

!--------------------------------------------------------------------------
type(EBSDdefectIV_T) function EBSDdefectIV_constructor( nmlfile ) result(EBSDdefectIV)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDdefectIV_constructor
!! author: MDG 
!! version: 1.0 
!! date: 01/30/26
!!
!! constructor for the EBSDdefectIV_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call EBSDdefectIV%readNameList(nmlfile)

end function EBSDdefectIV_constructor

!--------------------------------------------------------------------------
subroutine EBSDdefectIV_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 01/30/26
!!
!! destructor for the EBSDdefectIV_T Class
 
IMPLICIT NONE

type(EBSDdefectIV_T), INTENT(INOUT)  :: self 

call reportDestructor('EBSDdefectIV_T')

end subroutine EBSDdefectIV_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 01/30/26
!!
!! read the namelist from an nml file for the EBSDdefectIV_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(EBSDdefectIV_T), INTENT(INOUT) :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg)                    :: numsx
integer(kind=irg)                    :: numsy
integer(kind=irg)                    :: binning
integer(kind=irg)                    :: nthreads
real(kind=sgl)                       :: thetac
real(kind=sgl)                       :: delta
real(kind=sgl)                       :: spotsize
real(kind=sgl)                       :: gammavalue
real(kind=dbl)                       :: beamcurrent
real(kind=dbl)                       :: dwelltime
character(3)                         :: scalingmode
logical                              :: sampleInteractionVolume
character(fnlen)                     :: deformationfile
character(fnlen)                     :: ivolfile
character(fnlen)                     :: masterfile
character(fnlen)                     :: datafile

! define the IO namelist to facilitate passing variables to the program.
namelist /EBSDdefectIV/ numsx, numsy, binning, nthreads, thetac, delta, spotsize, gammavalue, beamcurrent, &
                        dwelltime, scalingmode, sampleInteractionVolume, deformationfile, ivolfile, masterfile, &
                        datafile

! note that the remaining pattern center parameters are read from the deformationfile...
numsx = 0
numsy = 0
binning = 1
nthreads = 1
thetac = 10.0
delta = 50.0
spotsize = 1.0
gammavalue = 0.3333
beamcurrent = 1.0
dwelltime = 1.0
scalingmode = 'gam'
sampleInteractionVolume = .TRUE.
deformationfile = 'undefined'
ivolfile = 'undefined'
masterfile = 'undefined'
datafile = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=EBSDdefectIV)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(datafile).eq.'undefined') then
  call Message%printError('readNameList:',' datafile file name is undefined in '//nmlfile)
 end if

  if (trim(masterfile).eq.'undefined') then
  call Message%printError('readNameList:',' masterfile file name is undefined in '//nmlfile)
 end if

 if (trim(ivolfile).eq.'undefined') then
  call Message%printError('readNameList:',' ivolfile file name is undefined in '//nmlfile)
 end if

 if (trim(deformationfile).eq.'undefined') then
  call Message%printError('readNameList:',' deformationfile file name is undefined in '//nmlfile)
 end if

end if

self%nml%numsx = numsx
self%nml%numsy = numsy
self%nml%binning = binning
self%nml%nthreads = nthreads
self%nml%thetac = thetac
self%nml%delta = delta
self%nml%spotsize = spotsize
self%nml%gammavalue = gammavalue
self%nml%beamcurrent = beamcurrent
self%nml%dwelltime = dwelltime
self%nml%scalingmode = scalingmode
self%nml%sampleInteractionVolume = sampleInteractionVolume
self%nml%deformationfile = deformationfile
self%nml%ivolfile = ivolfile
self%nml%masterfile = masterfile
self%nml%datafile = datafile

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 01/30/26
!!
!! pass the namelist for the EBSDdefectIV_T Class to the calling program

IMPLICIT NONE 

class(EBSDdefectIV_T), INTENT(INOUT)          :: self
type(EBSDdefectIVNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 01/30/26
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT)    :: self 
type(HDF%T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 5, n_real = 4, n_double = 2
integer(kind=irg)                       :: hdferr,  io_int(n_int), NB=0 
real(kind=sgl)                          :: io_real(n_real)
real(kind=dbl)                          :: io_real(n_double)
character(20)                           :: intlist(n_int), reallist(n_real), doublelist(n_double)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( enl => self%nml )

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

if (enl%sampleInteractionVolume.eqv..TRUE.) NB=1

io_int = (/ enl%numsx, enl%numsy, enl%binning, enl%nthreads, NB /)
intlist(1) = 'numsx'
intlist(2) = 'numsy'
intlist(3) = 'binning'
intlist(4) = 'nthreads'
intlist(5) = 'sampleInteractionVolume'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the double reals
io_real = (/ enl%beamcurrent, enl%dwelltime /)
doublelist(1) = 'beamcurrent'
doublelist(2) = 'dwelltime'
call HDF%writeNMLdbles(io_double, doublelist, n_double)

dataset = 'scalingmode'
line2(1) = trim(enl%scalingmode)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create scalingmode dataset', hdferr)

dataset ='deformationfile' 
line2(1) = trim(enl%deformationfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create deformationfile dataset', hdferr)

dataset = SC_masterfile
line2(1) = trim(enl%masterfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create masterfile dataset', hdferr)

dataset ='ivolfile' 
line2(1) = trim(enl%ivolfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create ivolfile dataset', hdferr)

dataset = 'datafile'
line2(1) = trim(enl%datafile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create datafile dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine setnumsx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnumsx_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! set numsx in the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)        :: inp

self%nml%numsx = inp

end subroutine setnumsx_

!--------------------------------------------------------------------------
function getnumsx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnumsx_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! get numsx from the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
integer(kind=irg)                    :: out

out = self%nml%numsx

end function getnumsx_

!--------------------------------------------------------------------------
subroutine setnumsy_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnumsy_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! set numsy in the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)        :: inp

self%nml%numsy = inp

end subroutine setnumsy_

!--------------------------------------------------------------------------
function getnumsy_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnumsy_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! get numsy from the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
integer(kind=irg)                    :: out

out = self%nml%numsy

end function getnumsy_

!--------------------------------------------------------------------------
subroutine setbinning_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setbinning_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! set binning in the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)        :: inp

self%nml%binning = inp

end subroutine setbinning_

!--------------------------------------------------------------------------
function getbinning_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getbinning_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! get binning from the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
integer(kind=irg)                    :: out

out = self%nml%binning

end function getbinning_

!--------------------------------------------------------------------------
subroutine setnthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnthreads_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! set nthreads in the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)        :: inp

self%nml%nthreads = inp

end subroutine setnthreads_

!--------------------------------------------------------------------------
function getnthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnthreads_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! get nthreads from the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
integer(kind=irg)                    :: out

out = self%nml%nthreads

end function getnthreads_

!--------------------------------------------------------------------------
subroutine setthetac_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setthetac_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! set thetac in the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
real(kind=sgl), INTENT(IN)           :: inp

self%nml%thetac = inp

end subroutine setthetac_

!--------------------------------------------------------------------------
function getthetac_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getthetac_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! get thetac from the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
real(kind=sgl)                       :: out

out = self%nml%thetac

end function getthetac_

!--------------------------------------------------------------------------
subroutine setdelta_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdelta_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! set delta in the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
real(kind=sgl), INTENT(IN)           :: inp

self%nml%delta = inp

end subroutine setdelta_

!--------------------------------------------------------------------------
function getdelta_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdelta_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! get delta from the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
real(kind=sgl)                       :: out

out = self%nml%delta

end function getdelta_

!--------------------------------------------------------------------------
subroutine setspotsize_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setspotsize_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! set spotsize in the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
real(kind=sgl), INTENT(IN)           :: inp

self%nml%spotsize = inp

end subroutine setspotsize_

!--------------------------------------------------------------------------
function getspotsize_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getspotsize_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! get spotsize from the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
real(kind=sgl)                       :: out

out = self%nml%spotsize

end function getspotsize_

!--------------------------------------------------------------------------
subroutine setgammavalue_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setgammavalue_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! set gammavalue in the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
real(kind=sgl), INTENT(IN)           :: inp

self%nml%gammavalue = inp

end subroutine setgammavalue_

!--------------------------------------------------------------------------
function getgammavalue_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getgammavalue_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! get gammavalue from the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
real(kind=sgl)                       :: out

out = self%nml%gammavalue

end function getgammavalue_

!--------------------------------------------------------------------------
subroutine setbeamcurrent_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setbeamcurrent_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! set beamcurrent in the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
real(kind=dbl), INTENT(IN)           :: inp

self%nml%beamcurrent = inp

end subroutine setbeamcurrent_

!--------------------------------------------------------------------------
function getbeamcurrent_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getbeamcurrent_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! get beamcurrent from the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
real(kind=dbl)                       :: out

out = self%nml%beamcurrent

end function getbeamcurrent_

!--------------------------------------------------------------------------
subroutine setdwelltime_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdwelltime_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! set dwelltime in the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
real(kind=dbl), INTENT(IN)           :: inp

self%nml%dwelltime = inp

end subroutine setdwelltime_

!--------------------------------------------------------------------------
function getdwelltime_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdwelltime_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! get dwelltime from the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
real(kind=dbl)                       :: out

out = self%nml%dwelltime

end function getdwelltime_

!--------------------------------------------------------------------------
subroutine setscalingmode_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setscalingmode_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! set scalingmode in the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
character(3), INTENT(IN)             :: inp

self%nml%scalingmode = trim(inp)

end subroutine setscalingmode_

!--------------------------------------------------------------------------
function getscalingmode_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getscalingmode_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! get scalingmode from the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
character(3)                         :: out

out = trim(self%nml%scalingmode)

end function getscalingmode_

!--------------------------------------------------------------------------
subroutine setsampleInteractionVolume_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setsampleInteractionVolume_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! set sampleInteractionVolume in the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
logical, INTENT(IN)                  :: inp

self%nml%sampleInteractionVolume = inp

end subroutine setsampleInteractionVolume_

!--------------------------------------------------------------------------
function getsampleInteractionVolume_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getsampleInteractionVolume_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! get sampleInteractionVolume from the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
logical                              :: out

out = self%nml%sampleInteractionVolume

end function getsampleInteractionVolume_

!--------------------------------------------------------------------------
subroutine setdeformationfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdeformationfile_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! set deformationfile in the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
character(fnlen), INTENT(IN)         :: inp

self%nml%deformationfile = trim(inp)

end subroutine setdeformationfile_

!--------------------------------------------------------------------------
function getdeformationfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdeformationfile_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! get deformationfile from the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
character(fnlen)                     :: out

out = trim(self%nml%deformationfile)

end function getdeformationfile_

!--------------------------------------------------------------------------
subroutine setivolfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setivolfile_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! set ivolfile in the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
character(fnlen), INTENT(IN)         :: inp

self%nml%ivolfile = trim(inp)

end subroutine setivolfile_

!--------------------------------------------------------------------------
function getivolfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getivolfile_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! get ivolfile from the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
character(fnlen)                     :: out

out = trim(self%nml%ivolfile)

end function getivolfile_

!--------------------------------------------------------------------------
subroutine setmasterfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setmasterfile_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! set masterfile in the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
character(fnlen), INTENT(IN)         :: inp

self%nml%masterfile = trim(inp)

end subroutine setmasterfile_

!--------------------------------------------------------------------------
function getmasterfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getmasterfile_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! get masterfile from the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
character(fnlen)                     :: out

out = trim(self%nml%masterfile)

end function getmasterfile_

!--------------------------------------------------------------------------
subroutine setdatafile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdatafile_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! set datafile in the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
character(fnlen), INTENT(IN)         :: inp

self%nml%datafile = trim(inp)

end subroutine setdatafile_

!--------------------------------------------------------------------------
function getdatafile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdatafile_
!! author: MDG
!! version: 1.0
!! date: 01/30/26
!!
!! get datafile from the EBSDdefectIV_T class

IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT) :: self
character(fnlen)                     :: out

out = trim(self%nml%datafile)

end function getdatafile_



!--------------------------------------------------------------------------
!
! SUBROUTINE:EBSDreadorpcdefHDF
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read angles, pattern centers, and deformation tensor field from an HDF5 file
!
!> @param enl EBSD name list structure
!> @param ipar integer parameter array 
!> @param fpar float parameter array 
!> @param orpcdef array of unit quaternions, pattern centers, and deformation tensors (output)
!
! file format description:
! HDF5 file with a single Group at the top level:  
!
! DeformationFieldInfo
!
! the following datasets should be inside this group
!
! integers:
! npix  : number of region-of-interest pixels along x
! npiy  : number of region-of-interest pixels along y
! npiz  : number of depth steps
!
! floats:
! stepx : stepsize in nm along x
! stepy : stepsize in nm along y
! stepz : stepsize in nm along z
! 
! eu(3) : (phi_1, Phi, phi_2) Euler angles for grain orientation (in degrees)
!
! ; pattern center coordinates in EMsoft convention
! pcx(npix,npiy) : pattern center x-coordinates for all points in ROI
! pcy(npix,npiy) : pattern center y-coordinates for all points in ROI
! L              : distance sample-scintillator in microns
!
! deftensor(9,npiz,npix,npiy) : deformation tensor components for each voxel in ROI volume
!    Note the order of the dimensions in this array!
!
!> @date 11/05/19 MDG 1.0 original
!--------------------------------------------------------------------------
recursive subroutine EBSDreadorpcdefHDF_(self, enl, ipar, fpar, orpcdef)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDreadorpcdefHDF_

use mod_EBSD
use mod_EBSDdefect
use mod_io
use mod_quaternions
use mod_rotations
use mod_HDFsupport


IMPLICIT NONE

class(EBSDdefectIV_T), INTENT(INOUT)        :: self
type(EBSDdefectNameListType),INTENT(INOUT)  :: enl
integer(kind=irg),INTENT(INOUT)             :: ipar(3)
real(kind=dbl),INTENT(INOUT)                :: fpar(4)
type(EBSDAnglePCDefType),INTENT(INOUT)      :: orpcdef

type(HDF%T)                                 :: HDF 
type(IO_T)                                  :: Message 

integer(kind=irg)                           :: io_int(1), i, hdferr, k, j, hdferr
integer(kind=irg)                           :: istat
character(fnlen)                            :: deformationfile, groupname, dataset
logical                                     :: g_exists, stat, readonly 
real(kind=dbl),allocatable                  :: pcxy(:,:), eu(:,:,:,:)
integer(HSIZE_T)                            :: dims1(1), dims2(2), dims3(3), dims4(4)

HDF= HDF%T()
deformationfile = EMsoft%generateFilePath('EMdatapathname', enl%deformationfile)

inquire(file=trim(deformationfile), exist=g_exists)

if (.not.g_exists) then
  call Message%printError('EBSDreadorpcdefHDF_','deformation field input file does not exist')
end if

! is this a proper HDF5 file ?
call h5fis_hdf5_f(trim(deformationfile), stat, hdferr)

if (stat.eqv..FALSE.) then ! the file exists, so let's open it an first make sure it is an EBSD dot product file
   call Message%printError('EBSDreadorpcdefHDF_','This is not a proper HDF5 file')
end if

! open the Monte Carlo file
readonly = .TRUE.
hdferr =  HDF%openFile(self%MCfile, readonly)

groupname = 'DeformationFieldInfo'
hdferr = HDF%openGroup(groupname)
call H5Lexists_f(HDF%getobjectID(),trim(groupname),g_exists, hdferr)
if (.not.g_exists) then
  call Message%printError('EBSDreadorpcdefHDF_','This HDF file does not contain deformation field data')
end if

! read the single integer parameters
ipar = 0
dataset = 'npix'  ! ipar(1)
call H5Lexists_f(HDF%head%next%objectID,trim(dataset),g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetInteger(dataset,  hdferr, ipar(1))
end if

dataset = 'npiy'  ! ipar(2)
call H5Lexists_f(HDF%head%next%objectID,trim(dataset),g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetInteger(dataset,  hdferr, ipar(2))
end if

dataset = 'npiz'  ! ipar(3)
call H5Lexists_f(HDF%head%next%objectID,trim(dataset),g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetInteger(dataset,  hdferr, ipar(3))
end if

! read single floats
fpar = 0.0
dataset = 'L'  ! fpar(1)
call H5Lexists_f(HDF%head%next%objectID,trim(dataset),g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetDouble(dataset,  hdferr, fpar(1))
end if

dataset = 'stepx'  ! fpar(2)
call H5Lexists_f(HDF%head%next%objectID,trim(dataset),g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetDouble(dataset,  hdferr, fpar(2))
end if

dataset = 'stepy'  ! fpar(3)
call H5Lexists_f(HDF%head%next%objectID,trim(dataset),g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetDouble(dataset,  hdferr, fpar(3))
end if

dataset = 'stepz'  ! fpar(4)
call H5Lexists_f(HDF%head%next%objectID,trim(dataset),g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetDouble(dataset,  hdferr, fpar(4))
end if

! read Euler angle triplet and convert to quaternion
! dataset = 'eu'
! call H5Lexists_f(HDF%head%next%objectID,trim(dataset),g_exists, hdferr)
! if (g_exists.eqv..TRUE.) then
!     call HDF%readDatasetDoubleArray2D(dataset, dims2,  hdferr, eu)
!     allocate(orpcdef%quatang(4,1),stat=istat)
!     eu(1,1) = eu(1,1) + 90.D0
!     orpcdef%quatang(1:4,1) = eu2qu(sngl(eu(1,1:3))*dtor)
! end if
dataset = 'eu'
call H5Lexists_f(HDF%head%next%objectID,trim(dataset),g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetDoubleArray(dataset, dims4, hdferr, eu)
    allocate(orpcdef%quatangfield(4,dims4(2),dims4(3),dims4(4)),stat=istat)
    do k=1,dims4(2)
        do i=1,dims4(3)
            do j=1,dims4(4)
                eu(1,k,i,j) = eu(1,k,i,j) + 90.D0
                orpcdef%quatangfield(1:4,k,i,j) = eu2qu(sngl(eu(1:3,k,i,j))*dtor)
             end do
        end do
    end do
end if


! read the pattern center data 
dataset = 'pcx'
call H5Lexists_f(HDF%head%next%objectID,trim(dataset),g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetDoubleArray(dataset, dims2,  hdferr, pcxy)
    allocate( orpcdef%pcfield(2,dims2(1),dims2(2)) )
    orpcdef%pcfield(1,:,:) = pcxy
    deallocate(pcxy)
end if

dataset = 'pcy'
call H5Lexists_f(HDF%head%next%objectID,trim(dataset),g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetDoubleArray(dataset, dims2,  hdferr, pcxy)
    orpcdef%pcfield(2,:,:) = pcxy
    deallocate(pcxy)
end if

! and finally, read the deformation field dataset
dataset = 'deftensor'
call H5Lexists_f(HDF%head%next%objectID,trim(dataset),g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetDoubleArray(dataset, dims4,  hdferr, orpcdef%deftensorfield)
end if

! close the group and file
call HDF%pop(.TRUE.)

call Message('')
call Message(' -> completed reading deformation field info from file '//trim(deformationfile))
call Message('')

end subroutine EBSDreadorpcdefHDF%

!--------------------------------------------------------------------------
subroutine EBSDdefectIV_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDdefectIV_
!! author: MDG 
!! version: 1.0 
!! date: 01/30/26
!!
!! perform the computations

use mod_EMsoft
use mod_HDFnames
use mod_io
use mod_HDFsupport
use mod_EBSD
use mod_MCfiles
use mod_MPfiles
use stringconstants

IMPLICIT NONE 

class(EBSDdefectIV_T), INTENT(INOUT)    :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

type(MCfile_T)                          :: MCFT
type(MPfile_T)                          :: MPFT
type(IO_T)                              :: Message

integer(kind=irg)                       :: res, error_cnt, hdferr, numangles, ipar(3), nx, ny, nz
integer(kind=irg)                       :: istat, sz(3), io_int(3)
real(kind=dbl)                          :: fpar(4)
real(kind=sgl)                          :: io_real(3)
logical                                 :: verbose
character(fnlen)                        :: writetofile

call openFortranHDFInterface()

associate( enl => self%nml, mcnl => MCFT%nml, &
           EBSDMCdata => MCFT%MCDT, EBSDMPdata => MPFT%MPDT, EBSDdetector => self%det )

! this program needs a lot of data
! 1. read the angle and deformation arrays from the deformation HDF file



EBSDreadorpcdefHDF













end subroutine EBSDdefectIV_



end module mod_EBSDdefectIV