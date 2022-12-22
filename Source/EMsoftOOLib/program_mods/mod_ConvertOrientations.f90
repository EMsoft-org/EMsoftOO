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

module mod_ConvertOrientations
  !! author: MDG
  !! version: 1.0
  !! date: 01/24/20
  !!
  !! class definition for the EMConvertOrientations program

use mod_kinds
use mod_global

IMPLICIT NONE

! namelist for the EMConvertOrientations program
type, public :: ConvertOrientationsNameListType
  integer(kind=irg)       :: reducetoRFZ
  character(fnlen)        :: xtalname
  character(fnlen)        :: inputfile
  character(fnlen)        :: outputfile
  character(2)            :: outputrepresentation
end type ConvertOrientationsNameListType

! class definition
type, public :: ConvertOrientations_T
private
  character(fnlen)                       :: nmldeffile = 'EMConvertOrientations.nml'
  type(ConvertOrientationsNameListType)  :: nml

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: ConvertOrientations_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: ConvertOrientations => ConvertOrientations_

end type ConvertOrientations_T

! the constructor routine for this class
interface ConvertOrientations_T
  module procedure ConvertOrientations_constructor
end interface ConvertOrientations_T

contains

!--------------------------------------------------------------------------
type(ConvertOrientations_T) function ConvertOrientations_constructor( nmlfile ) result(ConvertOrientations)
!DEC$ ATTRIBUTES DLLEXPORT :: ConvertOrientations_constructor
!! author: MDG
!! version: 1.0
!! date: 01/22/20
!!
!! constructor for the ConvertOrientations_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call ConvertOrientations%readNameList_(nmlfile)

end function ConvertOrientations_constructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 01/24/20
!!
!! read the namelist from an nml file for the ConvertOrientations_T Class

use mod_io

IMPLICIT NONE

class(ConvertOrientations_T), INTENT(INOUT) :: self
character(fnlen),INTENT(IN)                 :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)                 :: initonly
 !! fill in the default values only; do not read the file

type(IO_T)                                  :: Message
logical                                     :: skipread = .FALSE.
integer(kind=irg)                           :: cnt

integer(kind=irg)       :: reducetoRFZ
character(fnlen)        :: xtalname
character(fnlen)        :: inputfile
character(fnlen)        :: outputfile
character(2)            :: outputrepresentation

! define the IO namelist to facilitate passing variables to the program.
namelist  / EMConvertOrientations / inputfile, outputfile, outputrepresentation, &
                                    xtalname, reducetoRFZ

! initialize
reducetoRFZ = 1
inputfile = 'undefined'
outputfile = 'undefined'
outputrepresentation = 'xx'
xtalname = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=EMConvertOrientations)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(xtalname).eq.'undefined') then
  call Message%printError('readNameList:',' structure file name is undefined in '//nmlfile)
 end if
 if (trim(inputfile).eq.'undefined') then
  call Message%printError('readNameList:',' input file name is undefined in '//nmlfile)
 end if
 if (trim(outputfile).eq.'undefined') then
  call Message%printError('readNameList:',' output file name is undefined in '//nmlfile)
 end if
 if (outputrepresentation.eq.'xx') then
  call Message%printError('readNameList:',' output representation is undefined in '//nmlfile)
 end if
end if

self%nml%reducetoRFZ = reducetoRFZ
self%nml%inputfile = inputfile
self%nml%outputfile = outputfile
self%nml%outputrepresentation = outputrepresentation
self%nml%xtalname = xtalname

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 01/24/20
!!
!! pass the namelist for the ConvertOrientations_T Class to the calling program

IMPLICIT NONE

class(ConvertOrientations_T), INTENT(INOUT)          :: self
type(ConvertOrientationsNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine ConvertOrientations_(self, EMsoft)
!DEC$ ATTRIBUTES DLLEXPORT :: ConvertOrientations_
!! author: MDG
!! version: 1.0
!! date: 01/24/20
!!
!! perform the computations

use mod_EMsoft
use mod_so3
use mod_io
use mod_symmetry
use mod_crystallography
use mod_quaternions
use mod_rotations

IMPLICIT NONE

class(ConvertOrientations_T), INTENT(INOUT)       :: self
type(EMsoft_T), INTENT(INOUT)                     :: EMsoft

type(Cell_T)            :: cell
type(SpaceGroup_T)      :: SG
type(so3_T)             :: SO
type(QuaternionArray_T) :: QA, qsym, qAR
type(IO_T)              :: Message
type(q_T)               :: qu
type(Quaternion_T)      :: qq
type(r_T)               :: roFZ

integer(kind=irg)       :: i,j,k, ierr, io_int(3), pgnum
character(fnlen)        :: fname
integer(kind=irg)       :: FZtype, FZorder

! get the point group number from the crystal file
call cell%getCrystalData(self%nml%xtalname, SG, EMsoft)

! define the symmetry operators
pgnum = SG%getPGnumber()
call QA%QSym_init( pgnum, qsym )

! define the fundamental zone class
SO = so3_T( pgnum )

! print some output
call SO%getFZtypeandorder(FZtype, FZorder)
io_int = (/ pgnum, FZtype, FZorder /)
call Message%WriteValue('  Point group, type, and order: ',io_int,3)

! read the input orientations
fname = EMsoft%generateFilePath('EMdatapathname',self%nml%inputfile)
call SO%getOrientationsfromFile(fname)
call Message%printMessage(' - read orientation data from file '//trim(fname), frm="(/A)")

! apply the reduction to the RFZ ?
if (self%nml%reducetoRFZ.eq.1) then
  call SO%ReducelisttoRFZ( qsym )
end if

! and write the results to an output file
fname = EMsoft%generateFilePath('EMdatapathname',self%nml%outputfile)
call SO%writeOrientationstoFile( fname, self%nml%outputrepresentation, 'FZ' )

! print a final message
call Message%printMessage(' - wrote orientation data to file :'//trim(fname),"(/A/)")

end subroutine ConvertOrientations_

end module mod_ConvertOrientations
