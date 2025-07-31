! ###################################################################
! Copyright (c) 2013-2025, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_CSL
  !! author: MDG 
  !! version: 1.0 
  !! date: 07/21/25
  !!
  !! class definition for CSL boundaries

use mod_kinds
use mod_global

IMPLICIT NONE 

! class definition
type, public :: CSL_T

integer(kind=irg)              :: CSLnum

character(3),allocatable       :: CSLlabels(:)

integer(kind=irg),allocatable  :: CSLintegers(:,:)
                                                 
contains
private 
  procedure, pass(self) :: getCSLrod_

  generic, public :: getCSLrod => getCSLrod_

end type CSL_T

! the constructor routine for this class 
interface CSL_T
  module procedure CSL_constructor
end interface CSL_T

contains

!--------------------------------------------------------------------------
type(CSL_T) function CSL_constructor( ) result(CSL)
!! author: MDG 
!! version: 1.0 
!! date: 07/21/25
!!
!! constructor for the CSL_T Class; reads the name list 
 
IMPLICIT NONE

! initialize the CSL data
CSL%CSLnum = 29

allocate(CSL%CSLlabels(CSL%CSLnum))
CSL%CSLlabels = (/ 'I  ', '3  ', '5  ', '7  ', '9  ', '11 ', '13a', '13b', '15 ', &
                   '17a', '17b', '19a', '19b', '21a', '21b', '23 ', '25a', '25b', &
                   '27a', '27b', '29a', '29b', '31a', '31b', '33a', '33b', '33c', &
                   '35a', '35b' /)

allocate(CSL%CSLintegers(6,CSL%CSLnum))
CSL%CSLintegers = reshape((/ 0,1,0,1,0,1, &
                             1,3,1,3,1,3, &
                             1,3,0,1,0,1, &
                             1,5,1,5,1,5, &
                             1,4,1,4,0,1, &
                             1,3,1,3,0,1, &
                             1,5,0,1,0,1, &
                             1,7,1,7,1,7, &
                             2,5,1,5,0,1, &
                             1,4,0,1,0,1, &
                             2,5,2,5,1,5, &
                             1,6,1,6,0,1, &
                             1,4,1,4,1,4, &
                             1,9,1,9,1,9, &
                             1,3,1,6,1,6, &
                             1,3,1,9,1,9, &
                             1,7,0,1,0,1, &
                             1,3,1,3,1,9, &
                             1,5,1,5,0,1, &
                             2,7,1,7,0,1, &
                             2,5,0,1,0,1, &
                             2,7,2,7,1,7, &
                             1,11,1,11,1,11, &
                             2,5,1,5,1,5, &
                             1,8,1,8,0,1, &
                             3,11,1,11,1,11, &
                             2,5,2,5,0,1, &
                             1,4,1,8,1,8, &
                             3,11,3,11,1,11 /), (/ 6, CSL%CSLnum /))

end function CSL_constructor

!--------------------------------------------------------------------------
subroutine CSL_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 07/21/25
!!
!! destructor for the CSL_T Class
 
IMPLICIT NONE

type(CSL_T), INTENT(INOUT)  :: self 

call reportDestructor('CSL_T')

end subroutine CSL_destructor

!--------------------------------------------------------------------------
recursive function getCSLrod_(self, CSLlabel, CSLnumber, quat )  result(rod)
!DEC$ ATTRIBUTES DLLEXPORT :: getCSLrod_

use mod_rotations 
use mod_io
use mod_math
use mod_quaternions

IMPLICIT NONE

class(CSL_T),INTENT(IN)                   :: self
character(3),INTENT(IN)                   :: CSLlabel
integer(kind=irg),INTENT(OUT)             :: CSLnumber
type(Quaternion_T),INTENT(OUT),OPTIONAL   :: quat
type(r_T)                                 :: rod

type(IO_T)                                :: Message 
type(q_T)                                 :: qu 

integer(kind=irg)                         :: i
real(kind=dbl)                            :: l, rd(4)

! first find the sequential number for this boundary 
CSLnumber = 0
do i=1, self%CSLnum
  if (trim(self%CSLlabels(i)).eq.trim(CSLlabel)) CSLnumber = i
end do

if (CSLnumber.eq.0) then
  call Message%printError('getCSLrod','requested CSL type not recognized')
end if

rd(1:3) = (/ dble(self%CSLintegers(1,CSLnumber))/dble(self%CSLintegers(2,CSLnumber)), &
             dble(self%CSLintegers(3,CSLnumber))/dble(self%CSLintegers(4,CSLnumber)), &
             dble(self%CSLintegers(5,CSLnumber))/dble(self%CSLintegers(6,CSLnumber)) /)

l = vecnorm(rd(1:3))
rd(1:3) = rd(1:3)/l
rd(4) = l

rod = r_T( rdinp = rd )

if (present(quat)) then 
  qu = rod%rq()
  quat = Quaternion_T( qd = qu%q_copyd() )
  call quat%quat_normalize()
end if 

end function getCSLrod_



end module mod_CSL