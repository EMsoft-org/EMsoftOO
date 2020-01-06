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

!--------------------------------------------------------------------------
! EMsoft:quaternions.f90
!--------------------------------------------------------------------------
!
! MODULE: quaternions
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief module with basic quaternion functions (some overloaded operators)
!                    
!> @details   [verified against Mathematica Quaternion package on 3/15/12]
!>
!> REMEMBER THAT QUATERNION MULTIPLICATION IS NON-COMMUTATIVE !!!!!
!>
!> quaternions are defined as arrays of 4 single or double precision reals;
!> the first entry is the scalar part, the remaining three form the vector part.
!>
!> If you want to try it out, here is an example test program\n
!>\n
!> program qtest\n
!>\n
!>use local\n
!>use quaternions\n
!>\n
!>IMPLICIT NONE\n
!>\n
!>! for single precision, use the following lines\n
!>!real(kind=sgl)  ::  u(4), v(4), w(4)\n
!>!real(kind=sgl) :: x, a=2\n
!>! define two quaternions (single)\n
!>!u = (/1.0,2.0,3.0,4.0/)\n
!>!v = (/5.0,6.0,7.0,8.0/)\n
!>\n
!>! for double precision, uncomment the next set and comment the previous set lines\n
!>real(kind=dbl)  ::  u(4), v(4), w(4)\n
!>real(kind=dbl) :: x, a=2.D0\n
!>! double\n
!>u = (/1.D0,2.D0,3.D0,4.D0/)\n
!>v = (/5.D0,6.D0,7.D0,8.D0/)\n
!>\n
!>\n
!>write (stdout,*) ' quaternion u '\n
!>call quaternion_print(u)\n
!>write (stdout,*) ' quaternion v '\n
!>call quaternion_print(v)\n
!>\n
!>! next, do all the operations to make sure that they are correct\n
!>\n
!>write (stdout,*) '   addition u+v '\n
!>w = u+v\n
!>call quaternion_print(w)\n
!>write (stdout,*) ' '\n
!>\n
!>write (stdout,*) '   subtraction u-v '\n
!>w = u-v\n
!>call quaternion_print(w)\n
!>write (stdout,*) ' '\n
!>\n
!>write (stdout,*) '   scalar multiplication (both orderings)  '\n
!>w = a*u\n
!>call quaternion_print(w)\n
!>w = u*a\n
!>call quaternion_print(w)\n
!>write (stdout,*) ' '\n
!>\n
!>write (stdout,*) '   multiplication uv '\n
!>w = quat_mult(u,v)\n
!>call quaternion_print(w)\n
!>write (stdout,*) ' '\n
!>\n
!write (stdout,*) '   conjugate u '\n
!>w = conjg(u)\n
!>call quaternion_print(w)\n
!>write (stdout,*) ' '\n
!>\n
!>write (stdout,*) '   norm(u) '\n
!>x = cabs(u)\n
!>write (stdout,*) x\n
!>write (stdout,*) ' '\n
!>\n
!>write (stdout,*) '   division u/v '\n
!>w = quat_div(u,v)\n
!>call quaternion_print(w)\n
!>write (stdout,*) ' '\n
!>\n
!>end program qtest\n
!>
!
!> @date 03/15/12   MDG 1.0 original
!> @date 08/04/13   MDG 1.1 moved rotation conversion functions to rotations.f90
!> @date 08/12/13   MDG 2.0 re-defined quaternions to be arrays of 4 reals rather than by letter ... 
!> @date 02/06/15   MDG 2.1 added quat_slerp interpolation routines
!> @date 03/11/15   MDG 2.2 renamed quaternion vector rotation routines and split into active and passive
!> @date 03/11/15   MDG 2.3 redefined quaternion product using epsijk constant (see rotations tutorial paper)
!> @date 03/14/15   MDG 2.4 sign change in quaternion product; removed quat_Lpstar routines
!--------------------------------------------------------------------------
module mod_quaternions
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/03/20
  !!
  !! Quaternion arithmetic class 
  !!
  !! Quaternions are defined with the scalar part in position 1, and the vector part in positions 2:4.
  !!


use mod_global
use mod_kinds
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit

IMPLICIT NONE

  interface conjg 
    module procedure :: quatconjg
  end interface 

  interface abs
    module procedure :: quatnorm
  end interface 

private

  type, public :: T_QuaternionClass 
    !! Quaternion Class definition
    private
      real(kind=sgl), dimension(4) :: q
      real(kind=dbl), dimension(4) :: qd
      character(1)                 :: s

    contains
    private 
      procedure, pass(self) :: quatprint
      procedure, pass(self) :: quatadd
      procedure, pass(self) :: quatsubtract
      procedure, pass(self) :: quatmult
      procedure, pass(self) :: quatsmult
      procedure, pass(self) :: quatdiv
      procedure, pass(self) :: quatsdiv
      procedure, pass(self) :: quatinnerproduct
      procedure, pass(self) :: quatangle
      procedure, pass(self) :: quatLp
      procedure, pass(self) :: quatslerp
      
      procedure, pass(self), public :: quatconjg
      procedure, pass(self), public :: quatnorm


      generic, public :: quat_print => quatprint
      generic, public :: operator(+) => quatadd
      generic, public :: operator(-) => quatsubtract
      generic, public :: operator(*) => quatmult, quatsmult
      generic, public :: operator(/) => quatdiv, quatsdiv
      generic, public :: quat_innerproduct => quatinnerproduct
      generic, public :: quat_angle => quatangle
      generic, public :: quat_Lp => quatLp
      generic, public :: quat_slerp => quatslerp

  end type T_QuaternionClass 

! the constructor routine for this class 
  interface T_QuaternionClass
    module procedure :: Quaternion_constructor
  end interface T_QuaternionClass

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! We begin with the functions/subroutines that are public in this class
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! FUNCTION: Quaternion_constructor
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief initialize the Quaternion class 
!
!> @date  01/03/20 MDG 1.0 new function
!--------------------------------------------------------------------------
type(T_QuaternionClass) function Quaternion_constructor( q, qd ) result(Quat)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/03/20
  !!
  !! constructor for the Quaternion Class 
  
IMPLICIT NONE

  real(kind=sgl), INTENT(IN), OPTIONAL      :: q(4)
  real(kind=dbl), INTENT(IN), OPTIONAL      :: qd(4)

! fill in one or the other quaternion 
  if (present(q)) then 
    Quat % q = q 
    Quat % qd = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
    Quat % s = 's'
  end if 

  if (present(qd)) then 
    Quat % q = (/ 0.0, 0.0, 0.0, 0.0 /)
    Quat % qd = qd 
    Quat % s = 'd'
  end if 

end function Quaternion_constructor

!--------------------------------------------------------------------------
!
! SUBROUTINE: quatprint
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief print a quaternion (for debugging purposes mostly)
!
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
recursive subroutine quatprint(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/03/20
  !!
  !! print a single precision quaternion (for debugging purposes mostly)

use mod_io

IMPLICIT NONE 

  class(T_QuaternionClass),intent(in)    :: self
   !! input quaternion 

  type(T_IOClass)                        :: Message 

  if (self%s.eq.'s') then 
    call Message % WriteValue('', self%q, 4, frm="('(',4f12.6,')')")
  else 
    call Message % WriteValue('', self%qd, 4, frm="('(',4f12.6,')')")
  end if 

end subroutine quatprint

!--------------------------------------------------------------------------
!
! FUNCTION: quatadd
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief quaternion addition (single precision)
!
!> @date 1/05/20   MDG 1.0 original
!--------------------------------------------------------------------------
pure recursive function quatadd(self, y) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/05/20
  !!
  !! quaternion addition (single precision)

IMPLICIT NONE

  class(T_QuaternionClass),intent(in) :: self, y
  type(T_QuaternionClass)             :: qres 

  if (self%s.eq.'s') then
    qres % q = self%q + y%q 
    qres%s = 's'
  else
    qres % qd = self%qd + y%qd 
    qres%s = 'd'
  end if 

end function quatadd

!--------------------------------------------------------------------------
!
! FUNCTION: quatsubtract
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief quaternion subtraction (single precision)
!
!> @date 1/05/20   MDG 1.0 original
!--------------------------------------------------------------------------
recursive function quatsubtract(self, y) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/05/20
  !!
  !! quaternion subtraction (single precision)

IMPLICIT NONE

  class(T_QuaternionClass),intent(in) :: self, y
  type(T_QuaternionClass)             :: qres 

  if (self%s.eq.'s') then 
    qres % q = self%q - y%q 
    qres%s = 's'
  else
    qres % qd = self%qd - y%qd 
    qres%s = 'd'
  end if

end function quatsubtract

!--------------------------------------------------------------------------
!
! FUNCTION: quatmult
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief quaternion multiplication   (single precision)
!
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!> @date 3/11/15   MDG 3.0 redefined quaternion product (see rotations tutorial paper)
!--------------------------------------------------------------------------
pure recursive function quatmult(self, y) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/05/20
  !!
  !! quaternion multiplication   (single precision)

IMPLICIT NONE

  class(T_QuaternionClass),intent(in) :: self, y
   !! input quaternions
  type(T_QuaternionClass)             :: qres 
   !! output quaternion 

! the following is a way to reduce the number of multiplications
! needs to be tested and merged with epsijk approach
!
! QuatMul(QUAT *q1, QUAT *q2, QUAT *res){
! float A, B, C, D, E, F, G, H;
! A = (q1->w + q1->x)*(q2->w + q2->x);
! B = (q1->z - q1->y)*(q2->y - q2->z);
! C = (q1->w - q1->x)*(q2->y + q2->z); 
! D = (q1->y + q1->z)*(q2->w - q2->x);
! E = (q1->x + q1->z)*(q2->x + q2->y);
! F = (q1->x - q1->z)*(q2->x - q2->y);
! G = (q1->w + q1->y)*(q2->w - q2->z);
! H = (q1->w - q1->y)*(q2->w + q2->z);
! res->w = B + (-E - F + G + H) /2;
! res->x = A - (E + F + G + H)/2; 
! res->y = C + (E - F + G - H)/2; 
! res->z = D + (E - F - G + H)/2;
! }

  if (self%s.eq.'s') then 
    qres%q = (/ self%q(1)*y%q(1) - self%q(2)*y%q(2) -          ( self%q(3)*y%q(3) + self%q(4)*y%q(4) ), &
                self%q(1)*y%q(2) + self%q(2)*y%q(1) + epsijk * ( self%q(3)*y%q(4) - self%q(4)*y%q(3) ), &
                self%q(1)*y%q(3) + self%q(3)*y%q(1) + epsijk * ( self%q(4)*y%q(2) - self%q(2)*y%q(4) ), &
                self%q(1)*y%q(4) + self%q(4)*y%q(1) + epsijk * ( self%q(2)*y%q(3) - self%q(3)*y%q(2) ) /)
    qres%s = 's'
  else 
    qres%qd = (/ self%qd(1)*y%qd(1) - self%qd(2)*y%qd(2) -           ( self%qd(3)*y%qd(3) + self%qd(4)*y%qd(4) ), &
                 self%qd(1)*y%qd(2) + self%qd(2)*y%qd(1) + epsijkd * ( self%qd(3)*y%qd(4) - self%qd(4)*y%qd(3) ), &
                 self%qd(1)*y%qd(3) + self%qd(3)*y%qd(1) + epsijkd * ( self%qd(4)*y%qd(2) - self%qd(2)*y%qd(4) ), &
                 self%qd(1)*y%qd(4) + self%qd(4)*y%qd(1) + epsijkd * ( self%qd(2)*y%qd(3) - self%qd(3)*y%qd(2) ) /) 
    qres%s = 'd'
  end if
      
end function quatmult

!--------------------------------------------------------------------------
!
! FUNCTION: quatsmult
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief quaternion multiplication   (single precision)
!
!> @date 01/03/20  MDG 1.0 original
!--------------------------------------------------------------------------
pure recursive function quatsmult(self, s) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/05/20
  !!
  !! scalar quaternion multiplication   (single precision)

IMPLICIT NONE

  class(T_QuaternionClass),intent(in)   :: self 
   !! input quaternion
  real(kind=sgl), INTENT(IN)            :: s
   !! scalar input 
  type(T_QuaternionClass)               :: qres 
   !! output quaternion

  qres%q = (/ s*self%q(1), s*self%q(2), s*self%q(3), s*self%q(4) /) 
  qres%s = 's'

end function quatsmult

!--------------------------------------------------------------------------
!
! FUNCTION: quatsmultd
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief quaternion multiplication   (single precision)
!
!> @date 01/03/20  MDG 1.0 original
!--------------------------------------------------------------------------
pure recursive function quatsmultd(self, s) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/05/20
  !!
  !! scalar quaternion multiplication   (double precision)

IMPLICIT NONE

  class(T_QuaternionClass),intent(in)   :: self 
   !! input quaternion
  real(kind=dbl), INTENT(IN)            :: s
   !! scalar input 
  type(T_QuaternionClass)               :: qres 
   !! output quaternion

  qres%qd = (/ s*self%qd(1), s*self%qd(2), s*self%qd(3), s*self%qd(4) /) 
  qres%s = 'd'

end function quatsmultd


!--------------------------------------------------------------------------
!
! FUNCTION: quatconjg
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief quaternion complex conjugation (extends intrinsic routine conjg)
!
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
pure recursive function quatconjg(self) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/05/20
  !!
  !! quaternion conjugation (extends intrinsic routine conjg)

IMPLICIT NONE 

  class(T_QuaternionClass),intent(in)    :: self
   !! input quaternion
  type(T_QuaternionClass)                :: qres
   !! output quaternion

  if (self%s.eq.'s') then 
    qres%q = (/ self%q(1), -self%q(2), -self%q(3), -self%q(4) /) 
    qres%s = 's'
  else
    qres%qd = (/ self%qd(1), -self%qd(2), -self%qd(3), -self%qd(4) /) 
    qres%s = 'd'
  end if 

end function quatconjg

!--------------------------------------------------------------------------
!
! FUNCTION: quatnorm
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief quaternion norm (extends intrinsic routine cabs)
!
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
pure recursive function quatnorm(self) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/05/20
  !!
  !! quaternion norm (extends intrinsic routine abs)

IMPLICIT NONE 

  class(T_QuaternionClass),intent(in) :: self
   !! input quaternion
  real(kind=sgl)                      :: res
   !! output norm
  
  if (self%s.eq.'s') then 
    res = sqrt( self%q(1)*self%q(1) + self%q(2)*self%q(2) + self%q(3)*self%q(3) + self%q(4)*self%q(4) )
  else
    res = sqrt( self%qd(1)*self%qd(1) + self%qd(2)*self%qd(2) + self%qd(3)*self%qd(3) + self%qd(4)*self%qd(4) )
  end if

end function quatnorm

!--------------------------------------------------------------------------
!
! FUNCTION: quatdiv
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief quaternion division (single precision)
!
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
recursive function quatdiv(self, y) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/05/20
  !!
  !! quaternion division (single precision)

IMPLICIT NONE 

  class(T_QuaternionClass),intent(in)    :: self, y
   !! input quaternions
  type(T_QuaternionClass)                :: qres
   !! output quaternion 

  type(T_QuaternionClass)                :: p, cy
  real(kind=sgl)                         :: q
  real(kind=dbl)                         :: qd

  if (self%s.eq.'s') then 
      q = quatnorm(y)
      cy = quatconjg(y)
      p = quatsdiv( cy, q*q )
      qres = quatmult(self,p)
      qres%s = self%s
  else
      qd = quatnorm(y)
      cy = quatconjg(y)
      p = quatsdivd( cy, qd*qd )
      qres = quatmult(self,p)
      qres%s = self%s
  end if 


end function quatdiv

!--------------------------------------------------------------------------
!
! FUNCTION: quatsdiv
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief quaternion division by scalar (single precision)
!
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
recursive function quatsdiv(self, s) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/05/20
  !!
  !! quaternion division (single precision)

use mod_io 

IMPLICIT NONE 

  class(T_QuaternionClass),intent(in)    :: self
   !! input quaternion (numerator)
  real(kind=sgl), INTENT(IN)             :: s            
   !! input quaternion (denominator)

  type(T_QuaternionClass)                :: qres
  type(T_IOClass)                        :: Message

  if (s.ne.0.0) then 
      qres%q = (/ self%q(1)/s, self%q(2)/s, self%q(3)/s, self%q(4)/s /) 
      qres%s = 's'
  else 
    call Message % printWarning('quatsdiv', (/ 'Attempting to divide quaternion by zero; skipping operation ...' /) )
    qres = self
  end if

end function quatsdiv

!--------------------------------------------------------------------------
!
! FUNCTION: quatsdivd
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief quaternion division (double precision)
!
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
recursive function quatsdivd(self, s) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/05/20
  !!
  !! quaternion division (doubgle precision)

use mod_io 

IMPLICIT NONE 

  class(T_QuaternionClass),intent(in)    :: self
   !! input quaternion (numerator)
  real(kind=dbl), INTENT(IN)             :: s            
   !! input quaternion (denominator)

  type(T_QuaternionClass)                :: qres
  type(T_IOClass)                        :: Message

  if (s.ne.0.0) then 
    qres%qd = (/ self%qd(1)/s, self%qd(2)/s, self%qd(3)/s, self%qd(4)/s /) 
    qres%s = 'd'
  else 
    call Message % printWarning('quatsdivd', (/ 'Attempting to divide quaternion by zero; skipping operation ...' /) )
    qres = self
  end if

end function quatsdivd

!--------------------------------------------------------------------------
!
! FUNCTION: quatinnerproduct
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  quaternion inner product  (single precision)
!
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
pure recursive function quatinnerproduct(self, y) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/05/20
  !!
  !! quaternion inner product (single precision)

IMPLICIT NONE 

  class(T_QuaternionClass),intent(in)    :: self, y
   !! input quaternions
  real(kind=sgl)                         :: res
   !! inner product 

  res = self%q(1) * y%q(1) + self%q(2) * y%q(2) + self%q(3) * y%q(3) + self%q(4) * y%q(4)

end function quatinnerproduct

!--------------------------------------------------------------------------
!
! FUNCTION: quatinnerproductd
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  quaternion inner product  (double precision)
!
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
pure recursive function quatinnerproductd(self, y) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/05/20
  !!
  !! quaternion inner product (double precision)

IMPLICIT NONE 

  class(T_QuaternionClass),intent(in)    :: self, y
   !! input quaternions
  real(kind=sgl)                         :: res
   !! inner product 

  res = self%qd(1) * y%qd(1) + self%qd(2) * y%qd(2) + self%qd(3) * y%qd(3) + self%qd(4) * y%qd(4)

end function quatinnerproductd

!--------------------------------------------------------------------------
!
! FUNCTION: quatangle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief   interquaternion angle   (single precision)
!
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------!
pure recursive function quatangle(self, y) result(res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/05/20
  !!
  !! interquaternion angle   (single precision)

IMPLICIT NONE 

  class(T_QuaternionClass),intent(in)    :: self, y
   !! input quaternions
  real(kind=sgl)                         :: res
   !! angle (radians)

  real(kind=sgl)                         :: q

  q = quatinnerproduct(self,y)
  res = acos( 2.0*q*q - 1.0 )

end function quatangle

!--------------------------------------------------------------------------
!
! FUNCTION: quatangled
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief   interquaternion angle   (double precision)
!
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------!
pure recursive function quatangled(self, y) result(res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/03/20
  !!
  !! interquaternion angle   (double precision)

IMPLICIT NONE 

   class(T_QuaternionClass),intent(in)   :: self, y
   !! input quaternions
  real(kind=sgl)                         :: res
   !! angle (radians)

  real(kind=dbl)                         :: q

  q = quatinnerproductd(self,y)
  res = dacos( 2.D0*q*q - 1.D0 )

end function quatangled

!--------------------------------------------------------------------------
!
! FUNCTION: quatLp
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  actively rotate a unit vector by a unit quaternion, L_p = p v p*
!
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!> @date 3/11/15   MDG 3.0 name change, to be compatible with rotations tutorial paper
!--------------------------------------------------------------------------!
pure recursive function quatLp(self, v) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/05/20
  !!
  !! actively rotate a unit vector by a unit quaternion, L_p = p v p* (single precision)

IMPLICIT NONE 

  class(T_QuaternionClass),intent(in)    :: self
   !! input quaternion
  real(kind=sgl),intent(in)              :: v(3)         
   !! input vector to be rotated 
  real(kind=sgl)                         :: res(3)
   !! output vector

  type(T_QuaternionClass)                :: qv, rqv, cq

  qv%q = (/ 0.0, v(1), v(2), v(3) /) 
  cq = quatconjg(self)
  rqv = quatmult(self, quatmult(qv, cq) )
  res(1:3) = rqv%q(2:4)

end function quatLp

!--------------------------------------------------------------------------
!
! FUNCTION: quatLpd
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief   actively rotate a unit vector by a unit quaternion, L_p = p v p*
!
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!> @date 3/11/15   MDG 3.0 name change, to be compatible with rotations tutorial paper
!--------------------------------------------------------------------------!
pure recursive function quatLpd(self, v) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/03/20
  !!
  !! actively rotate a unit vector by a unit quaternion, L_p = p v p* (double precision)

IMPLICIT NONE 

  class(T_QuaternionClass),intent(in)    :: self
   !! input quaternion
  real(kind=sgl),intent(in)              :: v(3)         
   !! input vector to be rotated 
  real(kind=dbl)                         :: res(3)
   !! output vector

  type(T_QuaternionClass)                :: qv, rqv, cq

  qv%q = (/ 0.0, v(1), v(2), v(3) /) 
  cq = quatconjg(self)
  rqv = quatmult(self, quatmult(qv, cq) )
  res(1:3) = rqv%q(2:4)

end function quatLpd

!--------------------------------------------------------------------------
!
! FUNCTION: quatslerp
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief return an array of interpolated quaternions
!
!> @date 02/06/15   MDG 1.0 original
!--------------------------------------------------------------------------!
pure recursive function quatslerp(self, qb, n) result(res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/05/20
  !!
  !! return an array of interpolated quaternions

IMPLICIT NONE 

  class(T_QuaternionClass),intent(in)         :: self   ! = qa 
   !! input quaternion (start)
  class(T_QuaternionClass),intent(in)         :: qb
   !! input quaternion (end)
  integer(kind=irg),intent(in)                :: n            
   !! number of steps in the interpolation
  type(T_QuaternionClass)                     :: res(n)
   !! output interpolated quaternion list

  type(T_QuaternionClass)                     :: cqa
  real(kind=sgl)                              :: theta, phi, dphi, s
  real(kind=dbl)                              :: thetad, phid, dphid, sd
  integer(kind=irg)                           :: i

  if (self%s.eq.'s') then 
      do i=1,n 
        res(i)%q = (/ 0.0, 0.0, 0.0, 0.0 /) 
        res(i)%s = 's'
      end do
      cqa = quatconjg(self)
      theta = acos( quatinnerproduct(qb, cqa ))
      if (theta.ne.0.0) then
        s = 1.0/sin(theta)
        dphi = theta/real(n-1)

        do i=1,n
          phi = real(i-1)*dphi
          res(i) = quatsmult(self, sin(theta-phi)*s ) + quatsmult(qb, sin(phi)*s )
        end do
      else
        do i=1,n
          res(i) = self
        end do
      end if
  else
      do i=1,n 
        res(i)%qd = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        res(i)%s = 'd'
      end do
      cqa = quatconjg(self)
      thetad = acos( quatinnerproduct(qb, cqa ))
      if (thetad.ne.0.D0) then
        sd = 1.D0/sin(theta)
        dphi = theta/dble(n-1)

        do i=1,n
          phi = dble(i-1)*dphi
          res(i) = quatsmultd(self, sin(theta-phi)*sd ) + quatsmultd(qb, sin(phi)*sd )
        end do
      else
        do i=1,n
          res(i) = self
        end do
      end if
  end if 

end function quatslerp

!--------------------------------------------------------------------------
!
! FUNCTION: quat_Marsaglia
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief return a single random quaternion using the Marsaglia approach
!
!> @param seed seed number 
! 
!> @date 04/23/18   MDG 1.0 original
!--------------------------------------------------------------------------!
recursive function quat_Marsaglia(seed) result(q)
!DEC$ ATTRIBUTES DLLEXPORT :: quatMarsagliad

use mod_rng 

IMPLICIT NONE

type(rng_t),INTENT(INOUT)           :: seed 
!f2py intent(in,out) ::  seed 
type(T_QuaternionClass)             :: q 

real(kind=dbl)                      :: x1,x2,y1,y2,s1,s2

x1 = 2.D0*rng_uniform(seed)-1.D0
y1 = 2.D0*rng_uniform(seed)-1.D0
s1 = x1*x1+y1*y1
if (s1.gt.1.D0) then 
  do while (s1.gt.1.D0) 
    x1 = 2.D0*rng_uniform(seed)-1.D0
    x2 = 2.D0*rng_uniform(seed)-1.D0
    s1 = x1*x1+y1*y1
  end do
end if 

x2 = 2.D0*rng_uniform(seed)-1.D0
y2 = 2.D0*rng_uniform(seed)-1.D0
s2 = x2*x2+y2*y2
if (s2.gt.1.D0) then 
  do while (s2.gt.1.D0) 
    x2 = 2.D0*rng_uniform(seed)-1.D0
    y2 = 2.D0*rng_uniform(seed)-1.D0
    s2 = x2*x2+y2*y2
  end do
end if 

s1 = sqrt( (1.D0-s2)/s2 )

q%qd = (/ x1, y1, x2*s1, y2*s1 /) 
q%s = 'd'

end function quat_Marsaglia


end module mod_quaternions
