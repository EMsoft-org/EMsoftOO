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

module mod_PGA3D
  !! author: MDG 
  !! version: 1.0 
  !! date: 07/20/21
  !!
  !! class definition for the PGA3D module (Projective Geometric Algebra in 3D) 
  !!
  !! This is based on the file pga3d.cpp from BiVector.net

use mod_kinds
use mod_global

IMPLICIT NONE 



private
character(5), parameter :: basis(16) = (/ "    1","   e0","   e1","   e2","   e3", &
                                          "  e01","  e02","  e03","  e12","  e31","  e23", &
                                          " e021"," e013"," e032"," e123","e0123" /)


intrinsic :: conjg 
public :: conjg

interface conjg
  procedure mvec_conjugate_
end interface conjg

! class definition
type, public :: PGA3D_T
private 
  real(kind=dbl)        :: mvec(0:15)

contains
private 
  procedure, pass(self) :: mvec_setcomp_
  procedure, pass(self) :: mvec_getcomp_
  procedure, pass(self) :: mvec_log_
  procedure, pass(self) :: mvec_conjugate_
  procedure, pass(self) :: mvec_involute_
  procedure, pass(self) :: mvec_mult_
  procedure, pass(self) :: mvec_wedge_
  procedure, pass(self) :: mvec_vee_
  procedure, pass(self) :: mvec_inner_
  procedure, pass(self) :: mvec_add_
  procedure, pass(self) :: mvec_subtract_
  procedure, pass(self) :: mvec_muls_
  procedure, pass(self) :: mvec_adds_
  ! procedure, pass(self) :: mvec_norm_
  ! procedure, pass(self) :: mvec_inorm_
  ! procedure, pass(self) :: mvec_normalized_
  procedure, pass(self) :: reverseOrder
  procedure, pass(self) :: PoincareDuality

  generic, public :: setcomp => mvec_setcomp_
  generic, public :: getcomp => mvec_getcomp_
  generic, public :: log => mvec_log_
  generic, public :: operator(.involute.) => mvec_involute_
  generic, public :: operator(*) => mvec_mult_
  generic, public :: operator(.wedge.) => mvec_wedge_
  generic, public :: operator(.vee.) => mvec_vee_
  generic, public :: operator(.inner.) => mvec_inner_
  generic, public :: operator(+) => mvec_add_, mvec_adds_
  generic, public :: operator(-) => mvec_subtract_
  generic, public :: operator(.muls.) => mvec_muls_ ! , mvec_smul_
  ! generic, public :: norm => mvec_norm_
  ! generic, public :: inorm => mvec_inorm_
  ! generic, public :: normalized => mvec_normalized_
  generic,public :: operator(.reverse.) => reverseOrder 
  generic,public :: operator(.dual.) => PoincareDuality

end type PGA3D_T

! interface operator(.smul.)
!   recursive function mvec_smul_(a, b) result(mvout)
!     use mod_global
!     real(kind=dbl), INTENT(IN)           :: a
!     type(PGA3D_T), INTENT(IN)            :: b 
!     type(PGA3D_T)                        :: mvout
!   end function mvec_smul_
! end interface 


! the constructor routine for this class 
interface PGA3D_T
  module procedure PGA3D_constructor
end interface PGA3D_T


contains

!--------------------------------------------------------------------------
type(PGA3D_T) function PGA3D_constructor( val, ind ) result(PGA3D)
!! author: MDG 
!! version: 1.0 
!! date: 07/20/21
!!
!! constructor for the PGA3D_T Class; initialize the multivector mvec
 
IMPLICIT NONE

real(kind=dbl),INTENT(IN),OPTIONAL      :: val
integer(kind=irg),INTENT(IN),OPTIONAL   :: ind 

PGA3D%mvec = 0.D0

if (present(val).and.present(ind)) then 
  PGA3D%mvec(ind) = val
end if 

end function PGA3D_constructor

!--------------------------------------------------------------------------
subroutine PGA3D_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 07/20/21
!!
!! destructor for the PGA3D_T Class
 
IMPLICIT NONE

type(PGA3D_T), INTENT(INOUT)  :: self 

call reportDestructor('PGA3D_T')

end subroutine PGA3D_destructor


!--------------------------------------------------------------------------
function reverseOrder( self ) result(mvout)
!DEC$ ATTRIBUTES DLLEXPORT :: reverseOrder
!! author: MDG 
!! version: 1.0 
!! date: 07/20/21
!!
!! reverse a multivector using the .reverse. operator

IMPLICIT NONE 

class(PGA3D_T),INTENT(IN)  :: self
type(PGA3D_T)              :: mvout 

mvout%mvec = self%mvec
mvout%mvec(5:14) = -mvout%mvec(5:14)

end function reverseOrder

!--------------------------------------------------------------------------
function PoincareDuality( self ) result(mvout)
!DEC$ ATTRIBUTES DLLEXPORT :: PoincareDuality
!! author: MDG 
!! version: 1.0 
!! date: 07/20/21
!!
!! Poincare duality operation 

IMPLICIT NONE 

class(PGA3D_T),INTENT(IN)   :: self
type(PGA3D_T)               :: mvout 

integer(kind=irg)           :: i

do i=0,15
  mvout%mvec(i) = self%mvec(16-i)
end do

end function PoincareDuality


!--------------------------------------------------------------------------
recursive function mvec_getcomp_(self, ind) result(c)
!DEC$ ATTRIBUTES DLLEXPORT :: mvec_getcomp_
!! author: MDG 
!! version: 1.0 
!! date: 07/20/21
!!
!! get a multivector component

IMPLICIT NONE 

class(PGA3D_T), INTENT(INOUT)         :: self
integer(kind=irg),INTENT(IN)          :: ind
real(kind=dbl)                        :: c 

c = self%mvec(ind)

end function mvec_getcomp_

!--------------------------------------------------------------------------
recursive subroutine mvec_setcomp_(self, val, ind) 
!DEC$ ATTRIBUTES DLLEXPORT :: mvec_setcomp_
!! author: MDG 
!! version: 1.0 
!! date: 07/20/21
!!
!! set a multivector component

IMPLICIT NONE 

class(PGA3D_T), INTENT(INOUT)         :: self
real(kind=dbl),INTENT(IN)             :: val 
integer(kind=irg),INTENT(IN)          :: ind

self%mvec(ind) = val

end subroutine mvec_setcomp_

!--------------------------------------------------------------------------
recursive function mvec_conjugate_(self) result(mvout)
!DEC$ ATTRIBUTES DLLEXPORT :: mvec_conjugate_
!! author: MDG 
!! version: 1.0 
!! date: 07/20/21
!!
!! conjugate a multivector

IMPLICIT NONE 

class(PGA3D_T), INTENT(INOUT)         :: self
type(PGA3D_T)                         :: mvout 

mvout%mvec = self%mvec
mvout%mvec(1:10) = -self%mvec(1:10)

end function mvec_conjugate_

!--------------------------------------------------------------------------
recursive function mvec_involute_(self) result(mvout)
!DEC$ ATTRIBUTES DLLEXPORT :: mvec_involute_
!! author: MDG 
!! version: 1.0 
!! date: 07/20/21
!!
!! involute a multivector

IMPLICIT NONE 

class(PGA3D_T), INTENT(IN)         :: self
type(PGA3D_T)                      :: mvout 

mvout%mvec = self%mvec
mvout%mvec(1:4) = -self%mvec(1:4)
mvout%mvec(11:14) = -self%mvec(11:14)

end function mvec_involute_

!--------------------------------------------------------------------------
recursive function mvec_mult_(self, mvb) result(mvout)
!DEC$ ATTRIBUTES DLLEXPORT :: mvec_mult_
!! author: MDG 
!! version: 1.0 
!! date: 07/20/21
!!
!! involute a multivector

IMPLICIT NONE 

class(PGA3D_T), INTENT(IN)          :: self
type(PGA3D_T), INTENT(IN)           :: mvb 
type(PGA3D_T)                       :: mvout

real(kind=dbl)                      :: a(0:15), b(0:15)

a = self%mvec
b = mvb%mvec

mvout%mvec = (/ b(0)*a(0)+b(2)*a(2)+b(3)*a(3)+b(4)*a(4)-b(8)*a(8)-b(9)*a(9)-b(10)*a(10)-b(14)*a(14), &
                b(1)*a(0)+b(0)*a(1)-b(5)*a(2)-b(6)*a(3)-b(7)*a(4)+b(2)*a(5)+b(3)*a(6)+b(4)*a(7) &
                +b(11)*a(8)+b(12)*a(9)+b(13)*a(10)+b(8)*a(11)+b(9)*a(12)+b(10)*a(13)+b(15)*a(14)-b(14)*a(15), &
                b(2)*a(0)+b(0)*a(2)-b(8)*a(3)+b(9)*a(4)+b(3)*a(8)-b(4)*a(9)-b(14)*a(10)-b(10)*a(14), &
                b(3)*a(0)+b(8)*a(2)+b(0)*a(3)-b(10)*a(4)-b(2)*a(8)-b(14)*a(9)+b(4)*a(10)-b(9)*a(14), &
                b(4)*a(0)-b(9)*a(2)+b(10)*a(3)+b(0)*a(4)-b(14)*a(8)+b(2)*a(9)-b(3)*a(10)-b(8)*a(14), &
                b(5)*a(0)+b(2)*a(1)-b(1)*a(2)-b(11)*a(3)+b(12)*a(4)+b(0)*a(5)-b(8)*a(6)+b(9)*a(7)+b(6)*a(8) &
                -b(7)*a(9)-b(15)*a(10)-b(3)*a(11)+b(4)*a(12)+b(14)*a(13)-b(13)*a(14)-b(10)*a(15), &
                b(6)*a(0)+b(3)*a(1)+b(11)*a(2)-b(1)*a(3)-b(13)*a(4)+b(8)*a(5)+b(0)*a(6)-b(10)*a(7)-b(5)*a(8) &
                -b(15)*a(9)+b(7)*a(10)+b(2)*a(11)+b(14)*a(12)-b(4)*a(13)-b(12)*a(14)-b(9)*a(15), &
                b(7)*a(0)+b(4)*a(1)-b(12)*a(2)+b(13)*a(3)-b(1)*a(4)-b(9)*a(5)+b(10)*a(6)+b(0)*a(7)-b(15)*a(8) &
                +b(5)*a(9)-b(6)*a(10)+b(14)*a(11)-b(2)*a(12)+b(3)*a(13)-b(11)*a(14)-b(8)*a(15), &
                b(8)*a(0)+b(3)*a(2)-b(2)*a(3)+b(14)*a(4)+b(0)*a(8)+b(10)*a(9)-b(9)*a(10)+b(4)*a(14), &
                b(9)*a(0)-b(4)*a(2)+b(14)*a(3)+b(2)*a(4)-b(10)*a(8)+b(0)*a(9)+b(8)*a(10)+b(3)*a(14), &
                b(10)*a(0)+b(14)*a(2)+b(4)*a(3)-b(3)*a(4)+b(9)*a(8)-b(8)*a(9)+b(0)*a(10)+b(2)*a(14), &
                b(11)*a(0)-b(8)*a(1)+b(6)*a(2)-b(5)*a(3)+b(15)*a(4)-b(3)*a(5)+b(2)*a(6)-b(14)*a(7)-b(1)*a(8) &
                +b(13)*a(9)-b(12)*a(10)+b(0)*a(11)+b(10)*a(12)-b(9)*a(13)+b(7)*a(14)-b(4)*a(15), &
                b(12)*a(0)-b(9)*a(1)-b(7)*a(2)+b(15)*a(3)+b(5)*a(4)+b(4)*a(5)-b(14)*a(6)-b(2)*a(7)-b(13)*a(8) &
                -b(1)*a(9)+b(11)*a(10)-b(10)*a(11)+b(0)*a(12)+b(8)*a(13)+b(6)*a(14)-b(3)*a(15), &
                b(13)*a(0)-b(10)*a(1)+b(15)*a(2)+b(7)*a(3)-b(6)*a(4)-b(14)*a(5)-b(4)*a(6)+b(3)*a(7)+b(12)*a(8) &
                -b(11)*a(9)-b(1)*a(10)+b(9)*a(11)-b(8)*a(12)+b(0)*a(13)+b(5)*a(14)-b(2)*a(15), &
                b(14)*a(0)+b(10)*a(2)+b(9)*a(3)+b(8)*a(4)+b(4)*a(8)+b(3)*a(9)+b(2)*a(10)+b(0)*a(14), &
                b(15)*a(0)+b(14)*a(1)+b(13)*a(2)+b(12)*a(3)+b(11)*a(4)+b(10)*a(5)+b(9)*a(6)+b(8)*a(7)+b(7)*a(8) &
                +b(6)*a(9)+b(5)*a(10)-b(4)*a(11)-b(3)*a(12)-b(2)*a(13)-b(1)*a(14)+b(0)*a(15) /)

end function mvec_mult_

!--------------------------------------------------------------------------
recursive function mvec_wedge_(self, mvb) result(mvout)
!DEC$ ATTRIBUTES DLLEXPORT :: mvec_wedge_
!! author: MDG 
!! version: 1.0 
!! date: 07/20/21
!!
!! wedge product of two multivectors

IMPLICIT NONE 

class(PGA3D_T), INTENT(IN)          :: self
type(PGA3D_T), INTENT(IN)           :: mvb 
type(PGA3D_T)                       :: mvout

real(kind=dbl)                      :: a(0:15), b(0:15)

a = self%mvec
b = mvb%mvec

mvout%mvec = (/ b(0)*a(0), b(1)*a(0)+b(0)*a(1), b(2)*a(0)+b(0)*a(2), b(3)*a(0)+b(0)*a(3), b(4)*a(0)+b(0)*a(4), &
                b(5)*a(0)+b(2)*a(1)-b(1)*a(2)+b(0)*a(5), b(6)*a(0)+b(3)*a(1)-b(1)*a(3)+b(0)*a(6), &
                b(7)*a(0)+b(4)*a(1)-b(1)*a(4)+b(0)*a(7), b(8)*a(0)+b(3)*a(2)-b(2)*a(3)+b(0)*a(8), &
                b(9)*a(0)-b(4)*a(2)+b(2)*a(4)+b(0)*a(9), b(10)*a(0)+b(4)*a(3)-b(3)*a(4)+b(0)*a(10), &
                b(11)*a(0)-b(8)*a(1)+b(6)*a(2)-b(5)*a(3)-b(3)*a(5)+b(2)*a(6)-b(1)*a(8)+b(0)*a(11), &
                b(12)*a(0)-b(9)*a(1)-b(7)*a(2)+b(5)*a(4)+b(4)*a(5)-b(2)*a(7)-b(1)*a(9)+b(0)*a(12), &
                b(13)*a(0)-b(10)*a(1)+b(7)*a(3)-b(6)*a(4)-b(4)*a(6)+b(3)*a(7)-b(1)*a(10)+b(0)*a(13), &
                b(14)*a(0)+b(10)*a(2)+b(9)*a(3)+b(8)*a(4)+b(4)*a(8)+b(3)*a(9)+b(2)*a(10)+b(0)*a(14), &
                b(15)*a(0)+b(14)*a(1)+b(13)*a(2)+b(12)*a(3)+b(11)*a(4)+b(10)*a(5)+b(9)*a(6)+b(8)*a(7) &
                +b(7)*a(8)+b(6)*a(9)+b(5)*a(10)-b(4)*a(11)-b(3)*a(12)-b(2)*a(13)-b(1)*a(14)+b(0)*a(15) /)

end function mvec_wedge_

!--------------------------------------------------------------------------
recursive function mvec_vee_(self, mvb) result(mvout)
!DEC$ ATTRIBUTES DLLEXPORT :: mvec_vee_
!! author: MDG 
!! version: 1.0 
!! date: 07/20/21
!!
!! Vee (or regressive) product of two multivectors

IMPLICIT NONE 

class(PGA3D_T), INTENT(IN)          :: self
type(PGA3D_T), INTENT(IN)           :: mvb 
type(PGA3D_T)                       :: mvout

real(kind=dbl)                      :: a(0:15), b(0:15), c(0:15)
integer(kind=irg)                   :: i

a = self%mvec
b = mvb%mvec

c(0:15) = (/ 1*(a(15)*b(15)), -1*(a(14)*-1*b(15)+a(15)*b(14)*-1), -1*(a(13)*-1*b(15)+a(15)*b(13)*-1), &
            -1*(a(12)*-1*b(15)+a(15)*b(12)*-1), &
            -1*(a(11)*-1*b(15)+a(15)*b(11)*-1), 1*(a(10)*b(15)+a(13)*-1*b(14)*-1-a(14)*-1*b(13)*-1+a(15)*b(10)), &
            1*(a(9)*b(15)+a(12)*-1*b(14)*-1-a(14)*-1*b(12)*-1+a(15)*b(9)), &
            1*(a(8)*b(15)+a(11)*-1*b(14)*-1-a(14)*-1*b(11)*-1+a(15)*b(8)), &
            1*(a(7)*b(15)+a(12)*-1*b(13)*-1-a(13)*-1*b(12)*-1+a(15)*b(7)), &
            1*(a(6)*b(15)-a(11)*-1*b(13)*-1+a(13)*-1*b(11)*-1+a(15)*b(6)), &
            1*(a(5)*b(15)+a(11)*-1*b(12)*-1-a(12)*-1*b(11)*-1+a(15)*b(5)), &
            1*(a(4)*b(15)-a(7)*b(14)*-1+a(9)*b(13)*-1-a(10)*b(12)*-1-a(12)*-1*b(10)+a(13)*-1*b(9)-a(14)*-1*b(7)+a(15)*b(4)), &
            1*(a(3)*b(15)-a(6)*b(14)*-1-a(8)*b(13)*-1+a(10)*b(11)*-1+a(11)*-1*b(10)-a(13)*-1*b(8)-a(14)*-1*b(6)+a(15)*b(3)), &
            1*(a(2)*b(15)-a(5)*b(14)*-1+a(8)*b(12)*-1-a(9)*b(11)*-1-a(11)*-1*b(9)+a(12)*-1*b(8)-a(14)*-1*b(5)+a(15)*b(2)), &
            1*(a(1)*b(15)+a(5)*b(13)*-1+a(6)*b(12)*-1+a(7)*b(11)*-1+a(11)*-1*b(7)+a(12)*-1*b(6)+a(13)*-1*b(5)+a(15)*b(1)), &
            1*(a(0)*b(15)+a(1)*b(14)*-1+a(2)*b(13)*-1+a(3)*b(12)*-1+a(4)*b(11)*-1+a(5)*b(10)+a(6)*b(9)+a(7)*b(8)+a(8)*b(7)&
            +a(9)*b(6)+a(10)*b(5)-a(11)*-1*b(4)-a(12)*-1*b(3)-a(13)*-1*b(2)-a(14)*(-1)*b(1)+a(15)*b(0)) /) 

do i=0,15 
  mvout%mvec(i) = c(15-i)
end do 

end function mvec_vee_

!--------------------------------------------------------------------------
recursive function mvec_inner_(self, mvb) result(mvout)
!DEC$ ATTRIBUTES DLLEXPORT :: mvec_inner_
!! author: MDG 
!! version: 1.0 
!! date: 07/20/21
!!
!! Vee (or regressive) product of two multivectors

IMPLICIT NONE 

class(PGA3D_T), INTENT(IN)          :: self
type(PGA3D_T), INTENT(IN)           :: mvb 
type(PGA3D_T)                       :: mvout

real(kind=dbl)                      :: a(0:15), b(0:15), c(0:15)
integer(kind=irg)                   :: i

a = self%mvec
b = mvb%mvec

mvout%mvec = (/ b(0)*a(0)+b(2)*a(2)+b(3)*a(3)+b(4)*a(4)-b(8)*a(8)-b(9)*a(9)-b(10)*a(10)-b(14)*a(14), &
                b(1)*a(0)+b(0)*a(1)-b(5)*a(2)-b(6)*a(3)-b(7)*a(4)+b(2)*a(5)+b(3)*a(6)+b(4)*a(7)+b(11)*a(8) &
                +b(12)*a(9)+b(13)*a(10)+b(8)*a(11)+b(9)*a(12)+b(10)*a(13)+b(15)*a(14)-b(14)*a(15), &
                b(2)*a(0)+b(0)*a(2)-b(8)*a(3)+b(9)*a(4)+b(3)*a(8)-b(4)*a(9)-b(14)*a(10)-b(10)*a(14), &
                b(3)*a(0)+b(8)*a(2)+b(0)*a(3)-b(10)*a(4)-b(2)*a(8)-b(14)*a(9)+b(4)*a(10)-b(9)*a(14), &
                b(4)*a(0)-b(9)*a(2)+b(10)*a(3)+b(0)*a(4)-b(14)*a(8)+b(2)*a(9)-b(3)*a(10)-b(8)*a(14), &
                b(5)*a(0)-b(11)*a(3)+b(12)*a(4)+b(0)*a(5)-b(15)*a(10)-b(3)*a(11)+b(4)*a(12)-b(10)*a(15), &
                b(6)*a(0)+b(11)*a(2)-b(13)*a(4)+b(0)*a(6)-b(15)*a(9)+b(2)*a(11)-b(4)*a(13)-b(9)*a(15), &
                b(7)*a(0)-b(12)*a(2)+b(13)*a(3)+b(0)*a(7)-b(15)*a(8)-b(2)*a(12)+b(3)*a(13)-b(8)*a(15), &
                b(8)*a(0)+b(14)*a(4)+b(0)*a(8)+b(4)*a(14), &
                b(9)*a(0)+b(14)*a(3)+b(0)*a(9)+b(3)*a(14), &
                b(10)*a(0)+b(14)*a(2)+b(0)*a(10)+b(2)*a(14), &
                b(11)*a(0)+b(15)*a(4)+b(0)*a(11)-b(4)*a(15), &
                b(12)*a(0)+b(15)*a(3)+b(0)*a(12)-b(3)*a(15), &
                b(13)*a(0)+b(15)*a(2)+b(0)*a(13)-b(2)*a(15), &
                b(14)*a(0)+b(0)*a(14), &
                b(15)*a(0)+b(0)*a(15) /)

mvout%mvec(0) = 0.D0

end function mvec_inner_

!--------------------------------------------------------------------------
recursive function mvec_add_(self, mvb) result(mvout)
!DEC$ ATTRIBUTES DLLEXPORT :: mvec_add_
!! author: MDG 
!! version: 1.0 
!! date: 07/20/21
!!
!! add two multivectors

IMPLICIT NONE 

class(PGA3D_T), INTENT(IN)          :: self
type(PGA3D_T), INTENT(IN)           :: mvb 
type(PGA3D_T)                       :: mvout


mvout%mvec = self%mvec + mvb%mvec

end function mvec_add_

!--------------------------------------------------------------------------
recursive function mvec_adds_(self, s) result(mvout)
!DEC$ ATTRIBUTES DLLEXPORT :: mvec_adds_
!! author: MDG 
!! version: 1.0 
!! date: 07/20/21
!!
!! add a scalar to a multivector

IMPLICIT NONE 

class(PGA3D_T), INTENT(IN)          :: self
real(kind=dbl), INTENT(IN)          :: s 
type(PGA3D_T)                       :: mvout

mvout%mvec = self%mvec 
mvout%mvec(0) = mvout%mvec(0) + s 

end function mvec_adds_

!--------------------------------------------------------------------------
recursive function mvec_subtract_(self, mvb) result(mvout)
!DEC$ ATTRIBUTES DLLEXPORT :: mvec_subtract_
!! author: MDG 
!! version: 1.0 
!! date: 07/20/21
!!
!! add two multivectors

IMPLICIT NONE 

class(PGA3D_T), INTENT(IN)          :: self
type(PGA3D_T), INTENT(IN)           :: mvb 
type(PGA3D_T)                       :: mvout


mvout%mvec = self%mvec - mvb%mvec

end function mvec_subtract_

!--------------------------------------------------------------------------
recursive function mvec_smul_(a, b) result(mvout)
!DEC$ ATTRIBUTES DLLEXPORT :: mvec_smul_
!! author: MDG 
!! version: 1.0 
!! date: 07/20/21
!!
!! scalar times multivector

IMPLICIT NONE 

real(kind=dbl), INTENT(IN)           :: a
type(PGA3D_T), INTENT(IN)            :: b 
type(PGA3D_T)                        :: mvout

mvout%mvec = a * b%mvec 

end function mvec_smul_

!--------------------------------------------------------------------------
recursive function mvec_muls_(self, b) result(mvout)
!DEC$ ATTRIBUTES DLLEXPORT :: mvec_muls_
!! author: MDG 
!! version: 1.0 
!! date: 07/20/21
!!
!! multivector times scalar

IMPLICIT NONE 

class(PGA3D_T), INTENT(IN)          :: self
real(kind=dbl), INTENT(IN)          :: b 
type(PGA3D_T)                       :: mvout

mvout%mvec = self%mvec * b

end function mvec_muls_

!--------------------------------------------------------------------------
recursive subroutine mvec_log_(self, pre) 
!DEC$ ATTRIBUTES DLLEXPORT :: mvec_log_
!! author: MDG 
!! version: 1.0 
!! date: 07/20/21
!!
!! output the components of a multivector if they are non-zero

use mod_IO 

IMPLICIT NONE 

class(PGA3D_T), INTENT(INOUT)         :: self
character(*), INTENT(IN), OPTIONAL    :: pre 

integer(kind=irg)                     :: n, i, j
real(kind=dbl)                        :: c
character(fnlen)                      :: str 
character(3),parameter                :: plus = ' + '
character(3),parameter                :: minus = ' - '
character(3)                          :: sgn 
character(30)                         :: entry
type(IO_T)                            :: Message 

str = ''
if (present(pre)) call Message%printMessage(pre, frm="(A,' : ',$)")
do i=0,15
  c = self%mvec(i)
  if (c.ne.0.0) then 
! get the sign
    sgn = plus
    if (c.lt.0.0) then 
      sgn = minus
      c = -c
    end if
! create the substring entry
    entry = ''
    write (entry,"(F14.6)") c 
    entry = sgn//trim(adjustl(entry))//' '//trim(adjustl(basis(i)))
! add to the output string
    call Message%printMessage(entry, frm="(A,$)")
  end if 
end do 

call Message%printMessage(str, frm="(A)")

end subroutine mvec_log_



end module mod_PGA3D