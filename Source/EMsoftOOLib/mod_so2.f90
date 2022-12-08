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

module mod_so2
  !! author: MDG
  !! version: 1.0
  !! date: 01/21/20
  !!
  !! everything that has to do with sampling on the 2-sphere
  !!
  !! mostly used to sample points for BSE intensity computations

use mod_kinds
use mod_global
use mod_rotations

IMPLICIT NONE
private


type, public :: SO2pointd
        real(kind=dbl)          :: sql(2)        ! point in square lambert grid 
        real(kind=dbl)          :: nvec(3)       ! corresponding point on unit 2-sphere
        type(SO2pointd),pointer :: next          ! link to next point
end type SO2pointd


type, public :: so2_T
  private
    type(SO2pointd),pointer     :: SO2list
    integer(kind=irg)           :: SO2cnt

  contains
  private
    procedure, pass(self) :: delete_SO2list_
    procedure, pass(self) :: getSO2cnt_
    procedure, pass(self) :: getSO2listhead_
    final :: so2_destructor

    generic, public :: delete_SO2list => delete_SO2list_
    generic, public :: getSO2cnt => getSO2cnt_
    generic, public :: getSO2listhead => getSO2listhead_
end type so2_T

! the constructor routine for this class
interface so2_T
  module procedure so2_constructor
end interface so2_T

contains

!--------------------------------------------------------------------------
type(so2_T) function so2_constructor( nsteps, niz, noz ) result(SO)
!DEC$ ATTRIBUTES DLLEXPORT :: so2_constructor
!! author: MDG
!! version: 1.0
!! date: 12/06/22
!!
!! constructor for the so2_T Class
!> nsteps number of points along semi-edge of 2D Lambert space
!> niz z-component of unit vector pointing to inner edge of BSE detector
!> noz z-component of unit vector pointing to outer edge of BSE detector

use mod_Lambert

IMPLICIT NONE

integer(kind=irg),INTENT(IN)     :: nsteps 
real(kind=dbl),INTENT(IN)        :: niz 
real(kind=dbl),INTENT(IN)        :: noz 

type(Lambert_T)                  :: L  

real(kind=dbl)                   :: ds, xy(2), xyz(3) 
integer(kind=irg)                :: ierr, ix, iy 
type(SO2pointd), pointer         :: SO2tmp, SO2tmp2 

! set the counter to zero
SO%SO2cnt = 0

! empty the list if it isn't
! call SO%delete_SO2list_()

ds = 1.D0/dble(nsteps)

allocate(SO%SO2list)
SO2tmp => SO%SO2list
nullify(SO2tmp%next)

! loop through the square Lambert grid and test each point to see whether or not 
! it falls inside the annular BSE detector; this is a simple test in which we 
! compute the angle between the unit vector on the 2-sphere and the optical axis.
! We are actually just comparing the third direction cosine.
do ix = -nsteps,nsteps
    do iy = -nsteps,nsteps
        xy = (/ dble(ix), dble(iy) /) * ds 
        L = Lambert_T( xyd = xy )
        ierr = L%LambertSquareToSphere( xyz )
! does the z-component fall in the correct range ?
        if ((xyz(3).le.niz).and.(xyz(3).ge.noz)) then 
          ! and set the correct values 
          SO2tmp%sql = xy 
          SO2tmp%nvec = xyz
          SO%SO2cnt = SO%SO2cnt + 1 
! add this point to the linked list 
          allocate(SO2tmp%next)
          SO2tmp => SO2tmp%next
          nullify(SO2tmp%next)
        end if
    end do 
end do 

end function so2_constructor

!--------------------------------------------------------------------------
subroutine so2_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: so2_destructor
!! author: MDG
!! version: 1.0
!! date: 12/06/22
!!
!! destructor for the so2_T Class

IMPLICIT NONE

type(so2_T), INTENT(INOUT)  :: self

call reportDestructor('so2_T')

! nothing to do for now...

end subroutine so2_destructor

!--------------------------------------------------------------------------
function getSO2cnt_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getSO2cnt_
!! author: MDG
!! version: 1.0
!! date: 12/06/22
!!
!! get SO2cnt from the SO2_T class

IMPLICIT NONE

class(so2_T), INTENT(IN)     :: self
integer(kind=irg)            :: out

out = self%SO2cnt

end function getSO2cnt_

!--------------------------------------------------------------------------
function getSO2listhead_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getSO2listhead_
!! author: MDG
!! version: 1.0
!! date: 12/06/22
!!
!! get SO2list head from the SO2_T class

IMPLICIT NONE

class(so2_T), INTENT(IN)     :: self
type(SO2pointd),pointer      :: out

out => self%SO2list

end function getSO2listhead_

!--------------------------------------------------------------------------
subroutine delete_SO2list_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: delete_SO2list_
!! author: MDG
!! version: 1.0
!! date: 12/06/22
!!
!! nullify the SO2 list

IMPLICIT NONE

class(so2_T), INTENT(INOUT)           :: self
 
type(SO2pointd),pointer               :: ltail, ltmp

if (associated(self%SO2list)) then
  ltail => self%SO2list
  if (associated(ltail%next)) ltmp => ltail % next
  do 
    deallocate(ltail)
    if (.not. associated(ltmp)) EXIT
    ltail => ltmp
    ltmp => ltail % next
  end do 

  nullify(self%SO2list)
  self%SO2cnt = 0
end if

end subroutine delete_SO2list_

end module mod_so2