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

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: ComputeLauePatterns => ComputeLauePatterns_

end type Laue_T

! the constructor routine for this class 
interface Laue_T
  module procedure Laue_constructor
end interface Laue_T

contains

!--------------------------------------------------------------------------
type(Laue_T) function Laue_constructor( nmlfile ) result(Laue)
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
call HDF%writeNMLintegers(HDF_head, io_int, intlist, n_int)

! write all the single reals
io_real = (/ lnl%spotw, lnl%pixelsize, lnl%maxVoltage, lnl%minVoltage, lnl%SDdistance, lnl%gammavalue /)
reallist(1) = 'spotw'
reallist(2) = 'pixelsize'
reallist(3) = 'maxVoltage'
reallist(4) = 'minVoltage'
reallist(5) = 'SDdistance'
reallist(6) = 'gammavalue'
call HDF%writeNMLreals(HDF_head, io_real, reallist, n_real)

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
subroutine ComputeLauePatterns_(self, EMsoft, progname)
!! author: MDG 
!! version: 1.0 
!! date: 05/11/21
!!
!! perform the computations

use mod_EMsoft

IMPLICIT NONE 

class(Laue_T), INTENT(INOUT)            :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 


end subroutine ComputeLauePatterns_



end module mod_Laue