! ###################################################################
! Copyright (c) 2013-2021, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_Laue
  !! author: MDG 
  !! version: 1.0 
  !! date: 05/11/21
  !!
  !! class definition for the EMLaue program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMLaue program
type, public :: LaueNameListType
        integer(kind=irg)       :: numpx
        integer(kind=irg)       :: numpy
        integer(kind=irg)       :: nthreads
        integer(kind=irg)       :: BPx
        real(kind=sgl)          :: spotw
        real(kind=sgl)          :: pixelsize
        real(kind=sgl)          :: maxVoltage
        real(kind=sgl)          :: minVoltage
        real(kind=sgl)          :: SDdistance
        real(kind=sgl)          :: gammavalue
        character(fnlen)        :: backprojection
        character(fnlen)        :: Lauemode
        character(fnlen)        :: orientationfile
        character(fnlen)        :: tiffprefix
        character(fnlen)        :: hdfname
        character(fnlen)        :: xtalname
end type LaueNameListType

! class definition
type, public :: Laue_T
private 
  character(fnlen)              :: nmldeffile = 'EMLaue.nml'
  type(LaueNameListType)        :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: ComputeLauePatterns_
  procedure, pass(self) :: get_numpx_
  procedure, pass(self) :: get_numpy_
  procedure, pass(self) :: get_nthreads_
  procedure, pass(self) :: get_BPx_
  procedure, pass(self) :: get_spotw_
  procedure, pass(self) :: get_pixelsize_
  procedure, pass(self) :: get_maxVoltage_
  procedure, pass(self) :: get_minVoltage_
  procedure, pass(self) :: get_SDdistance_
  procedure, pass(self) :: get_gammavalue_
  procedure, pass(self) :: get_backprojection_
  procedure, pass(self) :: get_Lauemode_
  procedure, pass(self) :: get_orientationfile_
  procedure, pass(self) :: get_tiffprefix_
  procedure, pass(self) :: get_hdfname_
  procedure, pass(self) :: get_xtalname_
  procedure, pass(self) :: set_numpx_
  procedure, pass(self) :: set_numpy_
  procedure, pass(self) :: set_nthreads_
  procedure, pass(self) :: set_BPx_
  procedure, pass(self) :: set_spotw_
  procedure, pass(self) :: set_pixelsize_
  procedure, pass(self) :: set_maxVoltage_
  procedure, pass(self) :: set_minVoltage_
  procedure, pass(self) :: set_SDdistance_
  procedure, pass(self) :: set_gammavalue_
  procedure, pass(self) :: set_backprojection_
  procedure, pass(self) :: set_Lauemode_
  procedure, pass(self) :: set_orientationfile_
  procedure, pass(self) :: set_tiffprefix_
  procedure, pass(self) :: set_hdfname_
  procedure, pass(self) :: set_xtalname_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: ComputeLauePatterns => ComputeLauePatterns_
  generic, public :: get_numpx => get_numpx_
  generic, public :: get_numpy => get_numpy_
  generic, public :: get_nthreads => get_nthreads_
  generic, public :: get_BPx => get_BPx_
  generic, public :: get_spotw => get_spotw_
  generic, public :: get_pixelsize => get_pixelsize_
  generic, public :: get_maxVoltage => get_maxVoltage_
  generic, public :: get_minVoltage => get_minVoltage_
  generic, public :: get_SDdistance => get_SDdistance_
  generic, public :: get_gammavalue => get_gammavalue_
  generic, public :: get_backprojection => get_backprojection_
  generic, public :: get_Lauemode => get_Lauemode_
  generic, public :: get_orientationfile => get_orientationfile_
  generic, public :: get_tiffprefix => get_tiffprefix_
  generic, public :: get_hdfname => get_hdfname_
  generic, public :: get_xtalname => get_xtalname_
  generic, public :: set_numpx => set_numpx_
  generic, public :: set_numpy => set_numpy_
  generic, public :: set_nthreads => set_nthreads_
  generic, public :: set_BPx => set_BPx_
  generic, public :: set_spotw => set_spotw_
  generic, public :: set_pixelsize => set_pixelsize_
  generic, public :: set_maxVoltage => set_maxVoltage_
  generic, public :: set_minVoltage => set_minVoltage_
  generic, public :: set_SDdistance => set_SDdistance_
  generic, public :: set_gammavalue => set_gammavalue_
  generic, public :: set_backprojection => set_backprojection_
  generic, public :: set_Lauemode => set_Lauemode_
  generic, public :: set_orientationfile => set_orientationfile_
  generic, public :: set_tiffprefix => set_tiffprefix_
  generic, public :: set_hdfname => set_hdfname_
  generic, public :: set_xtalname => set_xtalname_
end type Laue_T

! the constructor routine for this class 
interface Laue_T
  module procedure Laue_constructor
end interface Laue_T

contains

!--------------------------------------------------------------------------
type(Laue_T) function Laue_constructor( nmlfile ) result(Laue)
!DEC$ ATTRIBUTES DLLEXPORT :: Laue_constructor
!! author: MDG 
!! version: 1.0 
!! date: 05/11/21
!!
!! constructor for the Laue_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call Laue%readNameList(nmlfile)

end function Laue_constructor

!--------------------------------------------------------------------------
subroutine Laue_destructor(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: Laue_destructor
!! author: MDG 
!! version: 1.0 
!! date: 05/11/21
!!
!! destructor for the Laue_T Class
 
IMPLICIT NONE

type(Laue_T), INTENT(INOUT)  :: self 

call reportDestructor('Laue_T')

end subroutine Laue_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 05/11/21
!!
!! read the namelist from an nml file for the Laue_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)         :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg)       :: numpx
integer(kind=irg)       :: numpy
integer(kind=irg)       :: nthreads
integer(kind=irg)       :: BPx
real(kind=sgl)          :: spotw
real(kind=sgl)          :: pixelsize
real(kind=sgl)          :: maxVoltage
real(kind=sgl)          :: minVoltage
real(kind=sgl)          :: SDdistance
real(kind=sgl)          :: gammavalue
character(fnlen)        :: backprojection
character(fnlen)        :: Lauemode
character(fnlen)        :: orientationfile
character(fnlen)        :: tiffprefix
character(fnlen)        :: xtalname
character(fnlen)        :: hdfname

! define the IO namelist to facilitate passing variables to the program.
namelist  / LaueData / numpx, numpy, nthreads, spotw, pixelsize, maxVoltage, minVoltage, SDdistance, &
                       gammavalue, Lauemode, orientationfile, tiffprefix, xtalname, hdfname, BPx, &
                       backprojection

numpx = 1024                   ! detector x-size (pixels)
numpy = 768                    ! detector y-size (pixels)
nthreads = 1                   ! number of parallel threads for pattern computation
BPx = 300                      ! semi-edge length for back projection square Lambert maps
pixelsize = 50.0               ! micron
spotw = 0.1                    ! spot size weight factor (1/(2*sigma^2))
maxVoltage = 30.0              ! in kV
minVoltage = 15.0              ! in kV
SDdistance = 100.0             ! mm
gammavalue = 1.0               ! scaling factor for gamma intensity scaling
backprojection = 'No'          ! 'Yes' or 'No'; adds backprojections to output file
Lauemode = 'transmission'      ! 'transmission' or 'reflection'
orientationfile = 'undefined'  ! input file with orientation list 
tiffprefix = 'undefined'       ! prefix for tiff output files with individual patterns
xtalname = 'undefined'         ! structure file name
hdfname = 'undefined'          ! HDF output file name

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=LaueData)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(xtalname).eq.'undefined') then
  call Message%printError('GetLaueNameList:',' crystal structure file name is undefined in '//nmlfile)
 end if
 if (trim(hdfname).eq.'undefined') then
  call Message%printError('GetLaueNameList:',' output file name is undefined in '//nmlfile)
 end if
 if (trim(orientationfile).eq.'undefined') then
  call Message%printError('GetLaueNameList:',' orientation file name is undefined in '//nmlfile)
 end if
end if

self%nml%numpx = numpx
self%nml%numpy = numpy
self%nml%nthreads = nthreads
self%nml%pixelsize = pixelsize
self%nml%BPx = BPx
self%nml%spotw = spotw
self%nml%maxVoltage= maxVoltage
self%nml%minVoltage= minVoltage
self%nml%SDdistance = SDdistance  
self%nml%gammavalue = gammavalue
self%nml%backprojection = backprojection
self%nml%Lauemode = Lauemode
self%nml%orientationfile = orientationfile
self%nml%xtalname = xtalname
self%nml%hdfname = hdfname
self%nml%tiffprefix = tiffprefix

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 05/11/21
!!
!! pass the namelist for the Laue_T Class to the calling program

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)          :: self
type(LaueNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 05/11/21
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(Laue_T), INTENT(INOUT)            :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 4, n_real = 6
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( lnl => self%nml )

! create the group for this namelist
groupname = trim(HDFnames%get_NMLlist())
hdferr = HDF%createGroup(groupname)

! write all the single integers
io_int = (/ lnl%numpx, lnl%numpy, lnl%nthreads, lnl%BPx /)
intlist(1) = 'numpx'
intlist(2) = 'numpy'
intlist(3) = 'nthreads'
intlist(4) = 'BPx'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ lnl%spotw, lnl%pixelsize, lnl%maxVoltage, lnl%minVoltage, lnl%SDdistance, lnl%gammavalue /)
reallist(1) = 'spotw'
reallist(2) = 'pixelsize'
reallist(3) = 'maxVoltage'
reallist(4) = 'minVoltage'
reallist(5) = 'SDdistance'
reallist(6) = 'gammavalue'
call HDF%writeNMLreals(io_real, reallist, n_real)

! write all the strings
dataset = SC_xtalname
line2(1) = lnl%xtalname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create xtalname dataset',hdferr)

dataset = 'hdfname'
line2(1) = lnl%hdfname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create hdfname dataset',hdferr)

dataset = 'tiffprefix'
line2(1) = lnl%tiffprefix
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create tiffprefix dataset',hdferr)

dataset = 'Lauemode'
line2(1) = lnl%Lauemode
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create Lauemode dataset',hdferr)

dataset = 'backprojection'
line2(1) = lnl%backprojection
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create backprojection dataset',hdferr)

dataset = 'orientationfile'
line2(1) = lnl%orientationfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create orientationfile dataset',hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
function get_numpx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_numpx_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! get numpx from the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
integer(kind=irg)              :: out

out = self%nml%numpx

end function get_numpx_

!--------------------------------------------------------------------------
subroutine set_numpx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_numpx_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! set numpx in the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
integer(kind=irg), INTENT(IN)  :: inp

self%nml%numpx = inp

end subroutine set_numpx_

!--------------------------------------------------------------------------
function get_numpy_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_numpy_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! get numpy from the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
integer(kind=irg)              :: out

out = self%nml%numpy

end function get_numpy_

!--------------------------------------------------------------------------
subroutine set_numpy_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_numpy_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! set numpy in the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
integer(kind=irg), INTENT(IN)  :: inp

self%nml%numpy = inp

end subroutine set_numpy_

!--------------------------------------------------------------------------
function get_nthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nthreads_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! get nthreads from the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
integer(kind=irg)              :: out

out = self%nml%nthreads

end function get_nthreads_

!--------------------------------------------------------------------------
subroutine set_nthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nthreads_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! set nthreads in the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
integer(kind=irg), INTENT(IN)  :: inp

self%nml%nthreads = inp

end subroutine set_nthreads_

!--------------------------------------------------------------------------
function get_BPx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_BPx_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! get BPx from the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
integer(kind=irg)              :: out

out = self%nml%BPx

end function get_BPx_

!--------------------------------------------------------------------------
subroutine set_BPx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_BPx_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! set BPx in the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
integer(kind=irg), INTENT(IN)  :: inp

self%nml%BPx = inp

end subroutine set_BPx_

!--------------------------------------------------------------------------
function get_spotw_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_spotw_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! get spotw from the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
real(kind=sgl)                 :: out

out = self%nml%spotw

end function get_spotw_

!--------------------------------------------------------------------------
subroutine set_spotw_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_spotw_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! set spotw in the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
real(kind=sgl), INTENT(IN)     :: inp

self%nml%spotw = inp

end subroutine set_spotw_

!--------------------------------------------------------------------------
function get_pixelsize_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_pixelsize_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! get pixelsize from the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
real(kind=sgl)                 :: out

out = self%nml%pixelsize

end function get_pixelsize_

!--------------------------------------------------------------------------
subroutine set_pixelsize_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_pixelsize_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! set pixelsize in the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
real(kind=sgl), INTENT(IN)     :: inp

self%nml%pixelsize = inp

end subroutine set_pixelsize_

!--------------------------------------------------------------------------
function get_maxVoltage_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_maxVoltage_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! get maxVoltage from the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
real(kind=sgl)                 :: out

out = self%nml%maxVoltage

end function get_maxVoltage_

!--------------------------------------------------------------------------
subroutine set_maxVoltage_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_maxVoltage_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! set maxVoltage in the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
real(kind=sgl), INTENT(IN)     :: inp

self%nml%maxVoltage = inp

end subroutine set_maxVoltage_

!--------------------------------------------------------------------------
function get_minVoltage_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_minVoltage_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! get minVoltage from the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
real(kind=sgl)                 :: out

out = self%nml%minVoltage

end function get_minVoltage_

!--------------------------------------------------------------------------
subroutine set_minVoltage_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_minVoltage_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! set minVoltage in the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
real(kind=sgl), INTENT(IN)     :: inp

self%nml%minVoltage = inp

end subroutine set_minVoltage_

!--------------------------------------------------------------------------
function get_SDdistance_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_SDdistance_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! get SDdistance from the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
real(kind=sgl)                 :: out

out = self%nml%SDdistance

end function get_SDdistance_

!--------------------------------------------------------------------------
subroutine set_SDdistance_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_SDdistance_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! set SDdistance in the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
real(kind=sgl), INTENT(IN)     :: inp

self%nml%SDdistance = inp

end subroutine set_SDdistance_

!--------------------------------------------------------------------------
function get_gammavalue_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_gammavalue_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! get gammavalue from the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
real(kind=sgl)                 :: out

out = self%nml%gammavalue

end function get_gammavalue_

!--------------------------------------------------------------------------
subroutine set_gammavalue_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_gammavalue_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! set gammavalue in the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
real(kind=sgl), INTENT(IN)     :: inp

self%nml%gammavalue = inp

end subroutine set_gammavalue_

!--------------------------------------------------------------------------
function get_backprojection_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_backprojection_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! get backprojection from the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
character(fnlen)               :: out

out = self%nml%backprojection

end function get_backprojection_

!--------------------------------------------------------------------------
subroutine set_backprojection_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_backprojection_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! set backprojection in the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
character(fnlen), INTENT(IN)   :: inp

self%nml%backprojection = inp

end subroutine set_backprojection_

!--------------------------------------------------------------------------
function get_Lauemode_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Lauemode_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! get Lauemode from the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
character(fnlen)               :: out

out = self%nml%Lauemode

end function get_Lauemode_

!--------------------------------------------------------------------------
subroutine set_Lauemode_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Lauemode_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! set Lauemode in the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
character(fnlen), INTENT(IN)   :: inp

self%nml%Lauemode = inp

end subroutine set_Lauemode_

!--------------------------------------------------------------------------
function get_orientationfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_orientationfile_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! get orientationfile from the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
character(fnlen)               :: out

out = self%nml%orientationfile

end function get_orientationfile_

!--------------------------------------------------------------------------
subroutine set_orientationfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_orientationfile_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! set orientationfile in the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
character(fnlen), INTENT(IN)   :: inp

self%nml%orientationfile = inp

end subroutine set_orientationfile_

!--------------------------------------------------------------------------
function get_tiffprefix_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_tiffprefix_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! get tiffprefix from the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
character(fnlen)               :: out

out = self%nml%tiffprefix

end function get_tiffprefix_

!--------------------------------------------------------------------------
subroutine set_tiffprefix_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_tiffprefix_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! set tiffprefix in the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
character(fnlen), INTENT(IN)   :: inp

self%nml%tiffprefix = inp

end subroutine set_tiffprefix_

!--------------------------------------------------------------------------
function get_hdfname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_hdfname_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! get hdfname from the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
character(fnlen)               :: out

out = self%nml%hdfname

end function get_hdfname_

!--------------------------------------------------------------------------
subroutine set_hdfname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_hdfname_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! set hdfname in the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
character(fnlen), INTENT(IN)   :: inp

self%nml%hdfname = inp

end subroutine set_hdfname_

!--------------------------------------------------------------------------
function get_xtalname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_xtalname_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! get xtalname from the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
character(fnlen)               :: out

out = self%nml%xtalname

end function get_xtalname_

!--------------------------------------------------------------------------
subroutine set_xtalname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_xtalname_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! set xtalname in the Laue class

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)   :: self
character(fnlen), INTENT(IN)   :: inp

self%nml%xtalname = inp

end subroutine set_xtalname_

!--------------------------------------------------------------------------
subroutine ComputeLauePatterns_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: ComputeLauePatterns_
!! author: MDG 
!! version: 1.0 
!! date: 05/25/21
!!
!! perform the computations

use mod_EMsoft
! use initializersHDF
use mod_initializers
use mod_diffraction
use mod_memory
use mod_xrd
use mod_symmetry
use mod_crystallography
use mod_rotations
use mod_quaternions
use mod_LaueSupport
use mod_io
use mod_so3
use mod_timing
use mod_Lambert
use HDF5
use mod_HDFnames
use mod_HDFsupport
use ISO_C_BINDING
use omp_lib
use mod_notifications
use stringconstants
use mod_image
use, intrinsic :: iso_fortran_env

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)            :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

type(Cell_T)                            :: cell
type(Timing_T)                          :: timer
type(IO_T)                              :: Message
type(HDF_T)                             :: HDF
type(SpaceGroup_T)                      :: SG
type(Diffraction_T)                     :: Diff
type(DynType)                           :: Dyn
type(memory_T)                          :: mem, memth
type(so3_T)                             :: SO
type(QuaternionArray_T)                 :: qAR
type(Quaternion_T)                      :: quat
type(LaueReflist_T)                     :: LaueReflist
type(q_T)                               :: qu

type(LaueReflist_T),pointer             :: reflist 
integer(kind=irg)                       :: numangles, numbatches, remainder, ii, jj, pid, tickstart, gcnt
integer(kind=irg),allocatable           :: batchnumangles(:)
integer(kind=irg),parameter             :: batchsize = 100

integer(kind=irg)                       :: i, hdferr, npx, npy, refcnt, io_int(1), NUMTHREADS, TID, BPnpx, BPnpy, info
real(kind=sgl)                          :: kouter, kinner, tstart, tstop, mi, ma, dmin, lambdamin
real(kind=sgl),allocatable              :: pattern(:,:), patternbatch(:,:,:), bppatterns(:,:,:), bp(:,:)
real(kind=dbl)                          :: intfactor
real(kind=dbl),allocatable              :: LegendreArray(:), upd(:), diagonal(:)

logical                                 :: verbose, f_exists, g_exists, insert=.TRUE., overwrite=.TRUE.

character(fnlen)                        :: hdfname, groupname, datagroupname, attributename, dataset, fname, TIFF_filename
character(11)                           :: dstr
character(15)                           :: tstrb
character(15)                           :: tstre
character(4)                            :: pnum
character(fnlen)                        :: HDF_FileVersion
integer(HSIZE_T)                        :: dims3(3), cnt3(3), offset3(3)
character(fnlen,kind=c_char)            :: line2(1)

! declare variables for use in object oriented image module
integer                                 :: iostat, Lstart
character(len=128)                      :: iomsg
logical                                 :: isInteger
type(image_t)                           :: im
integer(int8)                           :: i8 (3,4)
integer(int8), allocatable              :: TIFF_image(:,:)


! initialize the HDF class 
call openFortranHDFInterface()
HDF = HDF_T()

! simplify the notation a little
associate( lnl => self%nml )

! initialize the timing routines
timer = Timing_T()
tstrb = timer%getTimeString()
call timer%Time_tick(1)

! initialize the memory class 
mem = memory_T()

! rotations in double precision
call setRotationPrecision('d')

! read the list of orientations and convert them all to quaternions if they are not already
fname = EMsoft%generateFilePath('EMdatapathname',trim(lnl%orientationfile))
call SO%nullifyList()
call SO%getOrientationsfromFile(fname)
numangles = SO%getListCount('FZ')
call SO%listtoQuaternionArray( qAR )
call SO%delete_FZlist()
io_int(1) = numangles
call Message%WriteValue(' Number of orientations read from file: ', io_int, 1)

! compute the limiting wave numbers for the outer and inner Ewald spheres
kouter = getXRDwavenumber(lnl%maxVoltage)
kinner = getXRDwavenumber(lnl%minVoltage)

write (*,*) 'wave numbers : ', kouter, kinner 

!=============================================
!=============================================
! crystallography section 
verbose = .TRUE.

call cell%setFileName(lnl%xtalname)
call Diff%setrlpmethod('XR')

dmin = 0.05
call Diff%setV(dble(1.0))
call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, dmin, verbose, useHDF=HDF, noLUT=.TRUE.)

if ((SG%getSpaceGroupXtalSystem().eq.5).and.(cell%getLatParm('b').eq.cell%getLatParm('c'))) then
    call Message%printMessage( (/ &
    '                                                                         ', &
    ' ========Program Aborted========                                         ', &
    ' The Laue pattern simulation for rhombohedral/trigonal structures        ', &
    ' requires that the structure be described using the hexagonal reference  ', &
    ' frame.  Please re-enter the crystal structure in this setting and re-run', &
    ' the program.                                                            '/) )
    stop
end if

!=============================================
!=============================================
! compute possible reflection list with kinematical intensities
lambdamin = 1.0/kouter
intfactor = 0.0001D0   ! default intensity cutoff factor (from EMLauemaster program)
LaueReflist = LaueReflist_T()
call LaueReflist%Init_Reflist(cell, SG, Diff, gcnt, lambdamin, intfactor, verbose)

!=============================================
!=============================================
! start creation of the output file, using a hyperslab approach for the Laue patterns 

! Open a new file
  hdfname = trim(EMsoft%generateFilePath('EMdatapathname',lnl%hdfname))

! but delete it first if it already exists
  inquire(file=trim(hdfname), exist=f_exists)
  if (f_exists) then
    open(unit=dataunit, file=trim(hdfname), status='old',form='unformatted')
    close(unit=dataunit, status='delete')
  end if

  hdferr =  HDF%createFile(hdfname)

! write the EMheader to the file
datagroupname = trim(HDFnames%get_ProgramData())
  call HDF%writeEMheader(EMsoft, dstr, tstrb, tstre, progname, datagroupname)

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

  call self%writeHDFNameList(HDF, HDFnames)

! leave this group
  call HDF%pop()

! then the remainder of the data in a EMData group
  hdferr = HDF%createGroup(HDFnames%get_EMData())
  hdferr = HDF%createGroup(HDFnames%get_ProgramData())

! create the Lauemaster group and add a HDF_FileVersion attribbute to it 
  HDF_FileVersion = '4.0'
  HDF_FileVersion = cstringify(HDF_FileVersion)
  attributename = SC_HDFFileVersion
  hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)

! finally, write all the necessary data:  orientations and simulated patterns along with geometrical parameters
dataset = 'kouter'
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then 
      hdferr = HDF%writeDatasetFloat(dataset, kouter, overwrite)
    else
      hdferr = HDF%writeDatasetFloat(dataset, kouter)
    end if

dataset = 'kinner'
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then 
      hdferr = HDF%writeDatasetFloat(dataset, kinner, overwrite)
    else
      hdferr = HDF%writeDatasetFloat(dataset, kinner)
    end if

! adjust array dimension to be multiple of 100
  npx = lnl%numpx
  npy = lnl%numpy
  if (numangles.lt.batchsize) then 
    numbatches = 1
    call mem%alloc(batchnumangles, (/ numbatches /), 'batchnumangles')
    batchnumangles(1) = numangles 
    remainder = 0
    call mem%alloc(patternbatch, (/ npx, npy, numangles /), 'patternbatch', initval = 0.0)
  else
    numbatches = numangles/batchsize+1
    call mem%alloc(batchnumangles, (/ numbatches /), 'batchnumangles')
    batchnumangles(1:numbatches-1) = batchsize
    remainder = mod(numangles,batchsize)
    batchnumangles(numbatches) = remainder
    allocate(patternbatch(npx, npy, batchnumangles(1)))
    call mem%alloc(patternbatch, (/ npx, npy, batchnumangles(1) /), 'patternbatch', initval = 0.0)
  end if

dataset = 'numangles'
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then 
      hdferr = HDF%writeDatasetInteger(dataset, numangles, overwrite)
    else
      hdferr = HDF%writeDatasetInteger(dataset, numangles)
    end if

! create the hyperslabs and write zeroes to them for now
dataset = 'LauePatterns'
  if (numangles.lt.batchsize) then 
    dims3 = (/ npx, npy, numangles /)
    cnt3 = (/ npx, npy, numangles /)
    offset3 = (/ 0, 0, 0 /)
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then 
      hdferr = HDF%writeHyperslabFloatArray(dataset, patternbatch, dims3, offset3, cnt3, insert)
    else
      hdferr = HDF%writeHyperslabFloatArray(dataset, patternbatch, dims3, offset3, cnt3)
    end if
  else 
    dims3 = (/ npx, npy, numangles /)
    cnt3 = (/ npx, npy, batchnumangles(1) /)
    offset3 = (/ 0, 0, 0 /)
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then 
      hdferr = HDF%writeHyperslabFloatArray(dataset, patternbatch, dims3, offset3, cnt3, insert)
    else
      hdferr = HDF%writeHyperslabFloatArray(dataset, patternbatch, dims3, offset3, cnt3)
    end if
  end if

! should we add the backprojections to the same file ?
if (trim(lnl%backprojection).eq.'Yes') then 
  BPnpx = 2*lnl%BPx+1
  BPnpy = 2*lnl%BPx+1
  if (numangles.lt.batchsize) then 
    call mem%alloc(bppatterns, (/ BPnpx, BPnpy, numangles /), 'bppatterns', initval=0.0)
  else
    call mem%alloc(bppatterns, (/ BPnpx, BPnpy, batchnumangles(1) /), 'bppatterns', initval=0.0)
  end if

dataset = 'backprojections'
  if (numangles.lt.batchsize) then 
    dims3 = (/ BPnpx, BPnpy, numangles /)
    cnt3 = (/ BPnpx, BPnpy, numangles /)
    offset3 = (/ 0, 0, 0 /)
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then 
      hdferr = HDF%writeHyperslabFloatArray(dataset, bppatterns, dims3, offset3, cnt3, insert)
    else
      hdferr = HDF%writeHyperslabFloatArray(dataset, bppatterns, dims3, offset3, cnt3)
    end if
  else 
    dims3 = (/ BPnpx, BPnpy, numangles /)
    cnt3 = (/ BPnpx, BPnpy, batchnumangles(1) /)
    offset3 = (/ 0, 0, 0 /)
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then 
      hdferr = HDF%writeHyperslabFloatArray(dataset, bppatterns, dims3, offset3, cnt3, insert)
    else
      hdferr = HDF%writeHyperslabFloatArray(dataset, bppatterns, dims3, offset3, cnt3)
    end if
  end if
end if

! leave the output HDF5 file open so that we can write the hyperslabs as they are completed 

!=============================================
!=============================================
! precompute the Legendre array for the new lattitudinal grid values
if (trim(lnl%backprojection).eq.'Yes') then 
  call Message%printMessage(' Computing Legendre lattitudinal grid values')
  call mem%alloc(diagonal, (/ BPnpx /), 'diagonal', initval=0.D0)
  call mem%alloc(upd, (/ BPnpx /), 'upd')
  upd = (/ (dble(i) / dsqrt(4.D0 * dble(i)**2 - 1.D0), i=1,BPnpx) /)
  call dsterf(BPnpx-2, diagonal, upd, info) 
  ! the eigenvalues are stored from smallest to largest and we need them in the opposite direction
  call mem%alloc(LegendreArray, (/ BPnpx-1 /), 'LegendreArray', startdims = (/ 0 /) )
  LegendreArray(0:BPnpx-1) = diagonal(BPnpx:1:-1)
  ! set the center eigenvalue to 0
  LegendreArray((BPnpx-1)/2) = 0.D0
  call mem%dealloc(diagonal, 'diagonal') 
  call mem%dealloc(upd, 'upd')
end if 

!=============================================
!=============================================
! here we perform the actual simulations; we compute batchsize patterns at a time 
! and write them to the HDF5 output file, along with (optionally) the pattern tiff files
! set the number of OpenMP threads 
  call OMP_SET_NUM_THREADS(lnl%nthreads)
  io_int(1) = lnl%nthreads
  call Message%WriteValue(' Attempting to set number of threads to ',io_int, 1, frm = "(I4)")

Lstart = 8
memth = Memory_T( nt = lnl%nthreads )

! outer loop ... 
  do ii = 1, numbatches

! use OpenMP to run on multiple cores ... 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID, pid, jj, pattern, bp, quat, qu)

  NUMTHREADS = OMP_GET_NUM_THREADS()
  TID = OMP_GET_THREAD_NUM()
  call memth%alloc(pattern, (/ npx,npy /), 'pattern', TID=TID)
  patternbatch = 0.0

!$OMP DO SCHEDULE(DYNAMIC)
    do jj = 1,batchnumangles(ii)
      pid = (ii-1) * batchnumangles(1) + jj  
      quat = conjg( qAR%getQuatfromArray( pid ) )
      ! qu = q_T( qdinp = quat%get_quatd() )
      pattern = LaueReflist%getLauePattern(lnl%Lauemode, lnl%SDdistance, lnl%Pixelsize, lnl%spotw, &
                                           quat, kouter, kinner, npx, npy)
      patternbatch(1:npx,1:npy,jj) = pattern**lnl%gammavalue
    end do 
!$OMP END DO
    if (TID.eq.0) then 
      io_int(1) = ii 
      call Message%WriteValue(' patterns completed for batch ', io_int, 1)
    end if 
    call memth%dealloc(pattern, 'pattern', TID=TID)

    if (trim(lnl%backprojection).eq.'Yes') then 
      call memth%alloc(bp, (/ BPnpx, BPnpy /), 'bp', initval=0.0, TID=TID)

!$OMP DO SCHEDULE(DYNAMIC)
      do jj = 1,batchnumangles(ii)
        bp = backprojectLauePattern( (/kouter, kinner/), lnl%pixelsize, lnl%SDdistance, Lstart, &
                                     (/npx, npy/), (/ lnl%BPx, lnl%BPx /), patternbatch(1:npx,1:npy,jj), &
                                     lnl%Lauemode, LegendreArray)
        bppatterns(1:BPnpx,1:BPnpy,jj) = bp
      end do 
!$OMP END DO
    if (TID.eq.0) then 
      io_int(1) = ii 
      call Message%WriteValue(' backprojections completed for batch ', io_int, 1)
    end if 
    call memth%dealloc(bp, 'bp', TID=TID)
    end if 

! end of OpenMP portion
!$OMP END PARALLEL

! write the hyperslab to the HDF5 file 
 if (numangles.lt.100) then 
  dims3 = (/ npx, npy, numangles /)
  cnt3 = (/ npx, npy, numangles /)
  offset3 = (/ 0, 0, 0 /)
    hdferr = HDF%writeHyperslabFloatArray(dataset, patternbatch, dims3, offset3, cnt3, insert)
  else 
  dims3 = (/ npx, npy, numangles /)
  cnt3 = (/ npx, npy, batchnumangles(ii) /)
  offset3 = (/ 0, 0, (ii-1) * batchnumangles(1) /)
    hdferr = HDF%writeHyperslabFloatArray(dataset, patternbatch(:,:,1:batchnumangles(ii)), dims3, offset3, &
                                          cnt3, insert)
  end if


  if (trim(lnl%backprojection).eq.'Yes') then 
    if (numangles.lt.100) then 
      dims3 = (/ BPnpx, BPnpy, numangles /)
      cnt3 = (/ BPnpx, BPnpy, numangles /)
      offset3 = (/ 0, 0, 0 /)
        hdferr = HDF%writeHyperslabFloatArray(dataset, bppatterns, dims3, offset3, cnt3, insert)
    else 
      dims3 = (/ BPnpx, BPnpy, numangles /)
      cnt3 = (/ BPnpx, BPnpy, batchnumangles(ii) /)
      offset3 = (/ 0, 0, (ii-1) * batchnumangles(1) /)
        hdferr = HDF%writeHyperslabFloatArray(dataset, bppatterns(:,:,1:batchnumangles(ii)), dims3, offset3, &
                                          cnt3, insert)
    end if
  end if

! optionally, write the individual tiff image files 

! output the ADP map as a tiff file 
  do jj=1,batchnumangles(ii)
    write (pnum,"(I4.4)") (ii-1) * batchnumangles(1) + jj
    fname = EMsoft%generateFilePath('EMdatapathname',trim(lnl%tiffprefix)//'_'//pnum//'.tiff')
    TIFF_filename = trim(fname)

! allocate memory for image
    allocate(TIFF_image(npx,npy))

! fill the image with whatever data you have (between 0 and 255)
    ma = maxval(patternbatch(:,:,jj))
    mi = minval(patternbatch(:,:,jj))

    TIFF_image = int(255 * (patternbatch(:,:,jj)-mi)/(ma-mi))

! set up the image_t structure
    im = image_t(TIFF_image)
    if(im%empty()) call Message%printMessage("ComputeLauePattern","failed to convert array to image")

! create the file
    call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
    if(0.ne.iostat) call Message%printMessage("failed to write image to file : "//iomsg)
    deallocate(TIFF_image)
  end do

 end do ! outer loop
! 

 call HDF%pop()
 call HDF%pop()
 call timer%makeTimeStamp()
 dstr = timer%getDateString()
 tstre = timer%getTimeString()

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

 call HDF%pop(.TRUE.)

! and close the fortran hdf interface
call closeFortranHDFInterface()

! clean up allocated memory 
call mem%dealloc(patternbatch, 'patternbatch')
call mem%dealloc(batchnumangles, 'batchnumangles')
if (trim(lnl%backprojection).eq.'Yes') then 
  call mem%dealloc(LegendreArray, 'LegendreArray')
end if 

! call mem%allocated_memory_use()
! call memth%thread_memory_use()

end associate


end subroutine ComputeLauePatterns_



end module mod_Laue
