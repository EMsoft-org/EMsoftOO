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

module mod_Disorientations
  !! author: MDG
  !! version: 1.0
  !! date: 01/22/20
  !!
  !! class definition for the EMDisorientations program

use mod_kinds
use mod_global

IMPLICIT NONE

! namelist for the EMDisorientations program
type, public :: DisorientationsNameListType
  integer(kind=irg)       :: pgnum
  integer(kind=irg)       :: pgnum2
  character(fnlen)        :: inputfile1
  character(fnlen)        :: inputfile2
  character(fnlen)        :: outputfile
end type DisorientationsNameListType

! class definition
type, public :: Disorientations_T
private
  character(fnlen)                    :: nmldeffile = 'EMDisorientations.nml'
  type(DisorientationsNameListType)   :: nml

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: Disorientations_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: Disorientations => Disorientations_

end type Disorientations_T

! the constructor routine for this class
interface Disorientations_T
  module procedure Disorientations_constructor
end interface Disorientations_T

contains

!--------------------------------------------------------------------------
type(Disorientations_T) function Disorientations_constructor( nmlfile ) result(Disorientations)
!DEC$ ATTRIBUTES DLLEXPORT :: Disorientations_constructor
!! author: MDG
!! version: 1.0
!! date: 01/22/20
!!
!! constructor for the Disorientations_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call Disorientations%readNameList(nmlfile)

end function Disorientations_constructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 01/22/20
!!
!! read the namelist from an nml file for the Disorientations_T Class

use mod_io

IMPLICIT NONE

class(Disorientations_T), INTENT(INOUT)     :: self
character(fnlen),INTENT(IN)                 :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)                 :: initonly
 !! fill in the default values only; do not read the file

type(IO_T)                                  :: Message
logical                                     :: skipread = .FALSE.

integer(kind=irg)                           :: pgnum
integer(kind=irg)                           :: pgnum2
character(fnlen)                            :: inputfile1
character(fnlen)                            :: inputfile2
character(fnlen)                            :: outputfile

! define the IO namelist to facilitate passing variables to the program.
namelist /Disorientations/ pgnum, pgnum2, inputfile1, inputfile2, outputfile

! set the input parameters to default values
pgnum = 32                  !
pgnum2 = 0                  !
inputfile1 = 'undefined'    ! default filename for input file
inputfile2 = 'undefined'    ! default filename for input file
outputfile = 'undefined'    ! default filename for input file

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=Disorientations)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(inputfile1).eq.'undefined') then
  call Message%printError(' EMDisorientations',' input file name 1 is undefined in '//nmlfile)
 end if
 if (trim(inputfile2).eq.'undefined') then
  call Message%printError(' EMDisorientations',' input file name 2 is undefined in '//nmlfile)
 end if
 if (trim(outputfile).eq.'undefined') then
  call Message%printError(' EMDisorientations',' output file name is undefined in '//nmlfile)
 end if
end if

! if we get here, then all appears to be ok, and we need to fill in the emnl fields
self%nml%pgnum = pgnum
self%nml%pgnum2 = pgnum2
self%nml%inputfile1 = inputfile1
self%nml%inputfile2 = inputfile2
self%nml%outputfile = outputfile

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 01/22/20
!!
!! pass the namelist for the Disorientations_T Class to the calling program

IMPLICIT NONE

class(Disorientations_T), INTENT(INOUT)          :: self
type(DisorientationsNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine Disorientations_(self, EMsoft)
!DEC$ ATTRIBUTES DLLEXPORT :: Disorientations_
!! author: MDG
!! version: 1.0
!! date: 01/24/20
!!
!! perform the computations

use mod_kinds
use mod_global
use mod_EMsoft
use mod_rotations
use mod_quaternions
use mod_io
use mod_so3

IMPLICIT NONE

class(Disorientations_T), INTENT(INOUT)       :: self
type(EMsoft_T), INTENT(INOUT)                 :: EMsoft

type(so3_T)                                   :: SO1, SO2
type(QuaternionArray_T)                       :: qdummy, qAR1, qAR2
type(a_T)                                     :: disax
type(IO_T)                                    :: Message
type(r_T)                                     :: ro1, ro2
logical                                       :: singlePhase
integer(kind=irg)                             :: i, Pmdims1, Pmdims2, FZt, FZo
real(kind=dbl)                                :: ax(4)
type(FZpointd),pointer                        :: FZ1, FZ2
character(fnlen)                              :: fname

! one or two phases ?
singlePhase = .FALSE.
if (self%nml%pgnum2.eq.0) singlePhase = .TRUE.

! initialize the so3_T classes and symmetry operators
if (singlePhase.eqv..TRUE.) then
  SO1 = so3_T( self%nml%pgnum )
  call qdummy%QSym_Init( self%nml%pgnum, qAR1 )
  Pmdims1 = qAR1%getQnumber()
  SO2 = SO1
  qAR2 = qAR1
  Pmdims1 = Pmdims2
else
  SO1 = so3_T( self%nml%pgnum )
  call qdummy%QSym_Init( self%nml%pgnum, qAR1 )
  Pmdims1 = qAR1%getQnumber()
  SO2 = so3_T( self%nml%pgnum2 )
  call qdummy%QSym_Init( self%nml%pgnum2, qAR2 )
  Pmdims2 = qAR2%getQnumber()
end if

write (*,*) ' point groups and numbers ', self%nml%pgnum, self%nml%pgnum2, Pmdims1, Pmdims2
call SO1%getFZtypeandOrder(FZt, FZo)
write(*,*) ' FZ1 : ', FZt, FZo
call SO2%getFZtypeandOrder(FZt, FZo)
write(*,*) ' FZ2 : ', FZt, FZo

! get the orientations from files
fname = EMsoft%generateFilePath('EMdatapathname', self%nml%inputfile1)
call SO1%getOrientationsfromFile( fname  )
call Message%printMessage(' - read orientation data from file '//trim(fname), frm="(/A)")
fname = EMsoft%generateFilePath('EMdatapathname', self%nml%inputfile2)
call SO2%getOrientationsfromFile( fname )
call Message%printMessage(' - read orientation data from file '//trim(fname), frm="(/A)")

! only proceed if the two lists have the same number of entries
if (SO1%getListCount('FZ').ne.SO2%getListCount('FZ')) then
  call Message%printError('Disorientations', 'input files have different number of orientations in them...')
end if

! for each point, make sure it lies in the fundamental zone
call Message%printMessage('Reducing orientations to Fundamental Zone')
call SO1%ReducelisttoRFZ(qAR1)
call SO2%ReducelisttoRFZ(qAR2)

! and here we compute the disorientation angle from a quaternion product...
call Message%printMessage('Computing disorientations')

! and we write the output axis-angle pair (angle in degrees) to the outputfile
open(unit=10,file=trim(EMsoft%generateFilePath('EMdatapathname', self%nml%outputfile)),status='unknown',form='formatted')
write (10, "(I7)") SO1%getListCount('FZ')

FZ1 => SO1%getListHead('FZ')
FZ2 => SO2%getListHead('FZ')
do i=1,SO1%getListCount('FZ')
  if (singlePhase.eqv..TRUE.) then
    call SO1%getDisorientation(qAR1, FZ1%rod, FZ2%rod, disax)
  else
    call SO1%getDisorientation(qAR1, qAR2, FZ1%rod, FZ2%rod, disax)
  endif
  ax = disax%a_copyd()
  write (10,"(4F12.6)") ax(1:3), ax(4)/dtor
! next list entry
  FZ1 => FZ1%next
  FZ2 => FZ2%next
end do
close(10,status='keep')

call Message%printMessage(' Done.')

end subroutine Disorientations_



end module mod_Disorientations
