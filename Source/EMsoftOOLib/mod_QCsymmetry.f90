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

module mod_QCsymmetry
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/30/22
  !!
  !! class definition for 2D (axial) and 3D (icosahedral) quasi-crystal symmetry

use mod_kinds
use mod_global

IMPLICIT NONE 

private 

!>  SGname all icosahedral space group names
! ICOSAHEDRAL SPACE GROUPS
  character(13), public, dimension(11) :: SGname= (/  &
          " P 5 3 2     ", " P 5_1 3 2   ", " I 5 3 2     ", " I 5_1 3 2   ", &
          " F 5 3 2     ", " F 5_1 3 2   ", " P -5 -3 2/m ", " P -5 -3 2/q ", &
          " I -5 -3 2/m ", " F -5 -3 2/m ", " F -5 -3 2/q " /)

!>  GL encoded generator strings
  character(57), public, dimension(11) :: GL= (/  &
        "2b000000c000000", "2bA00000c000000", &
        "3aBBBBBBb000000c000000", "3aBBBBBBbA00000c000000", &
        "7aBB0000a0BB000a00BB00a000BB0a0000BBb000000c000000", &
        "7aBB0000a0BB000a00BB00a000BB0a0000BBbA00000c000000", &
        "3b000000c000000d000000", "3b000000c000000d000000", &
        "4aBBBBBBb000000c000000d000000", &
        "8aBB0000a0BB000a00BB00a000BB0a0000BBb000000c000000d000000", &
        "3b000000c000000d000000" /)


! class definition
type, public :: QCspacegroup_T
private 
  integer(kind=irg)     :: nsg 
  integer(kind=irg)     :: SGnum

contains
private 
  procedure, pass(self) :: printSGtable_
  procedure, pass(self) :: GetQCSpaceGroup_
  procedure, pass(self) :: 

  generic, public :: printSGtable => printSGtable_
  generic, public :: GetQCSpaceGroup => GetQCSpaceGroup_
  generic, public :: 

end type QCspacegroup_T

! the constructor routine for this class 
interface QCspacegroup_T
  module procedure QCspacegroup_constructor
end interface QCspacegroup_T

contains

!--------------------------------------------------------------------------
type(QCspacegroup_T) function QCspacegroup_constructor( ) result(QCSG)
!DEC$ ATTRIBUTES DLLEXPORT :: QCspacegroup_constructor
!! author: MDG 
!! version: 1.0 
!! date: 01/30/22
!!
!! constructor for the QCspacegroup_T Class
 
IMPLICIT NONE

integer(kind=irg)         :: sz(1) 

! get the number of space groups in the SGname array
sz = shape(SGname)
QCSG%nsg = sz(1)

end function QCspacegroup_constructor

!--------------------------------------------------------------------------
subroutine QCspacegroup_destructor(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: QCspacegroup_destructor
!! author: MDG 
!! version: 1.0 
!! date: 01/30/22
!!
!! destructor for the QCcell_T Class
 
IMPLICIT NONE

type(QCspacegroup_T), INTENT(INOUT)  :: self 

call reportDestructor('QCspacegroup_T')

end subroutine QCspacegroup_destructor

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
recursive subroutine printSGtable_(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: printSGtable_
!! author: MDG 
!! version: 1.0 
!! date: 01/30/22
!!
!! list the space groups... 
!!
!! Daniel Rokhsar, David Wright, David Mermin, Scale equivalence of 
!! quasicrystallographic space groups, PHYS REV B, 37(14), 1987.
 
use mod_io 

IMPLICIT NONE

class(QCspacegroup_T), INTENT(INOUT)  :: self 

type(IO_T)                            :: Message 

character(25)                         :: sgname1, sgname2
character(fnlen)                      :: str, str1
integer(kind=irg)                     :: ii

if(mod(self%nsg,2) .ne. 0) then
  do ii = 1, floor(float(self%nsg)/2.0) + 1
    if(ii .lt. floor(float(self%nsg)/2.0) + 1) then
      write(str,'(I3,A)') 2*ii-1,':'
      write(str1,'(I3,A)') 2*ii,':'
      sgname1 = ''
      sgname2 = ''
      sgname1 = trim(SGname(2*ii-1))
      sgname2 = trim(SGname(2*ii))
      call Message(trim(str)//' '//sgname1//trim(str1)//' '//sgname2)
    else
      write(str,'(I3,A)')2*ii-1,':'
      sgname1 = ''
      sgname1 = trim(SGname(2*ii-1))
      call Message(trim(str)//sgname1)
    end if
  end do
else
  do ii = 1, self%nsg/2
    write(str,'(I3,A)')2*ii-1,':'
    write(str1,'(I3,A)')2*ii,': '
    sgname1 = ''
    sgname2 = ''
    sgname1 = trim(SGname(2*ii-1))
    sgname2 = trim(SGname(2*ii))
    call Message(trim(str)//' '//sgname1//trim(str1)//' '//sgname2)
  end do
end if

end subroutine printSGtable_

!--------------------------------------------------------------------------
recursive subroutine GetQCSpaceGroup_(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: GetQCSpaceGroup_
!! author: MDG 
!! version: 1.0 
!! date: 01/30/22

class(QCspacegroup_T), INTENT(INOUT)  :: self 

type(IO_T)                            :: Message 

integer(kind=irg)                     :: io_int(1) 

 ! get the space group number
 call Message('Enter space group number', frm = "(//A)")

 call ReadValue('    Space group number --> ', io_int, 1)
 self%SGnum = io_int(1)

end subroutine GetQCSpaceGroup_


end module mod_QCsymmetry