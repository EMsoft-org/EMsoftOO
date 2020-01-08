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

module mod_quaternions
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/03/20
  !!
  !! Quaternion and Quaternion Array arithmetic class 
  !!
  !! Quaternions are defined with the scalar part in position 1, and the vector part in positions 2:4.
  !!
  !! There are two class definitions in this file, one for single quaternions, the other for 
  !! quaternion array operations (using OpenMP threads). The program MODQuaternionsTest.f90 
  !! can be used as part of ctest to run a test program on this module.


use mod_global
use mod_kinds
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit

IMPLICIT NONE
  private

! we overload the conjg, cabs, and .eq. intrinsics
    intrinsic :: conjg, cabs
    public :: conjg, cabs

    interface conjg
      procedure quatconjg
      procedure quatarrayconjg
    end interface conjg

    interface cabs
      procedure quatnorm
      procedure quatarraynorm
    end interface cabs

    interface operator(.eq.)
      procedure quatsequal
    end interface 

! definition of the quaternion class 
  type, public :: Quaternion_T
    !! Quaternion Class definition
    private
      real(kind=sgl), dimension(4) :: q
       !! single precision quaternion
      real(kind=dbl), dimension(4) :: qd
       !! double precision quaternion
      character(1)                 :: s
       !! precision indicator ('s' or 'd')

    contains
    private 
! quaternion IO routines
      procedure, pass(self) :: quatprint
! quaternion arithmetic routines 
      procedure, pass(self) :: quatadd
      procedure, pass(self) :: quatsubtract
      procedure, pass(self) :: quatmult
      procedure, pass(self) :: quatsmult
      procedure, pass(self) :: quatsmultd
      procedure, pass(self) :: quatdiv
      procedure, pass(self) :: quatsdiv
      procedure, pass(self) :: quatsdivd
      procedure, pass(self) :: quatconjg
      procedure, pass(self) :: quatnorm
      procedure, pass(self) :: quatnormalize
! quaternion-based transformations
      procedure, pass(self) :: quatLp
      procedure, pass(self) :: quatLpd
! routines with two or more input quaternions
      procedure, pass(self) :: quatinnerproduct
      procedure, pass(self) :: quatangle
      procedure, pass(self) :: quatslerp
! miscellaneous routines 
      procedure, pass(self), public :: quatsequal

      generic, public :: quat_print => quatprint
      generic, public :: operator(+) => quatadd
      generic, public :: operator(-) => quatsubtract
      generic, public :: operator(*) => quatmult 
      generic, public :: operator(*) => quatsmult, quatsmultd
      generic, public :: operator(/) => quatdiv 
      generic, public :: operator(/) => quatsdiv, quatsdivd
      generic, public :: quat_normalize => quatnormalize
      generic, public :: quat_Lp => quatLp, quatLpd
      generic, public :: quat_innerproduct => quatinnerproduct
      generic, public :: quat_angle => quatangle
      generic, public :: quat_slerp => quatslerp

  end type Quaternion_T 

! the constructor routine for this class 
  interface Quaternion_T
    module procedure Quaternion_constructor
  end interface Quaternion_T



! next we define the quaternion array class 
  type, public :: QuaternionArray_T
    !! Quaternion Class definition
    private
      integer(kind=irg)            :: n 
      integer(kind=irg)            :: nthreads
      real(kind=sgl), allocatable  :: q(:,:)
       !! single precision quaternion
      real(kind=dbl), allocatable  :: qd(:,:)
       !! double precision quaternion
      character(1)                 :: s
       !! precision indicator ('s' or 'd')

    contains
    private 
! quaternion arithmetic routines 
      procedure, pass(self) :: quatarrayadd
      procedure, pass(self) :: quatarraysubtract
      procedure, pass(self) :: quatarraymult
      procedure, pass(self) :: quatarraysmult
      procedure, pass(self) :: quatarraysmultd
      procedure, pass(self) :: quatarraydiv
      procedure, pass(self) :: quatarraysdiv
      procedure, pass(self) :: quatarrayconjg
      procedure, pass(self) :: quatarraynorm
      procedure, pass(self) :: quatarraynormalize
! quaternion-based transformations
      procedure, pass(self) :: quatarrayLp
      procedure, pass(self) :: quatarrayLpd
! routines with two or more input quaternion arrays
      procedure, pass(self) :: quatarrayinnerproduct
      procedure, pass(self) :: quatarrayangle
! miscellaneous routines 
      procedure, pass(self) :: extractfromQuaternionArray

! generics
      generic, public :: operator(+) => quatarrayadd
      generic, public :: operator(-) => quatarraysubtract
      generic, public :: operator(*) => quatarraymult 
      generic, public :: operator(*) => quatarraysmult, quatarraysmultd
      generic, public :: operator(/) => quatarraydiv 
      generic, public :: operator(/) => quatarraysdiv
      generic, public :: quat_normalize => quatarraynormalize
      generic, public :: quat_Lp => quatarrayLp, quatarrayLpd
      generic, public :: quat_innerproduct => quatarrayinnerproduct
      generic, public :: quat_angle => quatarrayangle
      generic, public :: getQuatfromArray => extractfromQuaternionArray

  end type QuaternionArray_T 

! the constructor routine for this class 
  interface QuaternionArray_T
    module procedure QuaternionArray_constructor
  end interface QuaternionArray_T

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! We begin with the functions/subroutines that are public in the 
! Quaternion_T class and pair up functions for individual and quaternion 
! arrays for easier module maintenance.
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
type(Quaternion_T) function Quaternion_constructor( q, qd ) result(Quat)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! constructor for the Quaternion Class 
  
IMPLICIT NONE

  real(kind=sgl), INTENT(IN), OPTIONAL      :: q(4)
  real(kind=dbl), INTENT(IN), OPTIONAL      :: qd(4)

! fill in one or the other quaternion 
  if ((.not.present(q)).and.(.not.present(qd))) then 
    Quat % q = (/ 0.0, 0.0, 0.0, 0.0 /)
    Quat % qd = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
    Quat % s = 's'
  else 
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
  end if

end function Quaternion_constructor

!--------------------------------------------------------------------------
type(QuaternionArray_T) function QuaternionArray_constructor( n, nthreads, q, qd, s ) result(QuatArray)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! constructor for the QuaternionArray Class 
  !!
  !! either call with parameters n and s 
  !! or with n and either one of q or qd 
  
use mod_io 

IMPLICIT NONE

  integer(kind=irg), INTENT(IN)             :: n
  integer(kind=irg), INTENT(IN), OPTIONAL   :: nthreads
  real(kind=sgl), INTENT(IN), OPTIONAL      :: q(4,n)
  real(kind=dbl), INTENT(IN), OPTIONAL      :: qd(4,n)
  character(1), INTENT(IN), OPTIONAL        :: s

  type(IO_T)                                :: Message

! OpenMP threads
  QuatArray % nthreads = 0
  if (present(nthreads)) QuatArray % nthreads = nthreads

! are we declaring just an empty variable with no entries, but with a given precision ?
  if ( present(s) .and. (.not.present(q)) .and. (.not.present(qd)) ) then 
    QuatArray % n = n 
    QuatArray % s = s 
    if (s.eq.'s') then 
      allocate(QuatArray % q(4,n))  
      QuatArray % q = 0.0
    else 
      allocate(QuatArray % qd(4,n))  
      QuatArray % qd = 0.D0
    end if 
    return 
  end if  
    
! number of quaternions in array
  QuatArray % n = n 

! single precision
  if (present(q)) then 
    allocate(QuatArray % q(4,n))  
    QuatArray % q = q
    QuatArray % s = 's'
  end if 

! double precision
  if (present(qd)) then 
    allocate(QuatArray % qd(4,n))  
    QuatArray % qd = qd
    QuatArray % s = 'd'
  end if 

end function QuaternionArray_constructor


!--------------------------------------------------------------------------
recursive subroutine quatprint(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! print a single precision quaternion (for debugging purposes mostly)

use mod_io

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion 

  type(IO_T)                        :: Message 

  if (self%s.eq.'s') then 
    call Message % WriteValue('', self%q, 4, frm="('(',4f12.6,'); precision: '$)")
    call Message % WriteValue('',self%s)
  else 
    call Message % WriteValue('', self%qd, 4, frm="('(',4f20.14,'); precision: '$)")
    call Message % WriteValue('',self%s)
  end if 

end subroutine quatprint

!--------------------------------------------------------------------------
pure recursive function quatadd(self, y) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion addition (single/double precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in) :: self, y
  type(Quaternion_T)             :: qres 

  if (self%s.eq.'s') then
    qres%q = self%q + y%q 
    qres%s = 's'
  else
    qres%qd = self%qd + y%qd 
    qres%s = 'd'
  end if 

end function quatadd

!--------------------------------------------------------------------------
recursive function quatarrayadd(self, y) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! quaternion array addition (single/double precision)

use mod_io

IMPLICIT NONE

  class(QuaternionArray_T),intent(in) :: self, y
  type(QuaternionArray_T)             :: qres 

  type(IO_T)                          :: Message 
  integer(kind=irg)                   :: sz(2)

! test to make sure that both arrays have the same number of quaternions 
  if (self%n.ne.y%n) then 
    call Message%printError('quatarrayadd','input arrays must have the same number of quaternions')
  end if 

  if (self%s.ne.y%s) then 
    call Message%printError('quatarrayadd','input arrays must have the same precision')
  end if 

  qres%n = self%n
  qres%s = self%s
  qres%nthreads = self%nthreads

  if (self%s.eq.'s') then
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%q)) then 
      sz = shape(qres%q)
      if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%q)
    end if 
    allocate(qres%q(4,self%n))
    qres%q = self%q + y%q 
  else
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%qd)) then 
      sz = shape(qres%qd)
      if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%qd)
    end if 
    allocate(qres%qd(4,self%n))
    qres%qd = self%qd + y%qd 
  end if 

end function quatarrayadd

!--------------------------------------------------------------------------
recursive function quatsubtract(self, y) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion subtraction (single/double precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in) :: self, y
  type(Quaternion_T)             :: qres 

  if (self%s.eq.'s') then 
    qres%q = self%q - y%q 
    qres%s = 's'
  else
    qres%qd = self%qd - y%qd 
    qres%s = 'd'
  end if

end function quatsubtract

!--------------------------------------------------------------------------
recursive function quatarraysubtract(self, y) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! quaternion array subtraction (single/double precision)

use mod_io 

IMPLICIT NONE

  class(QuaternionArray_T),intent(in) :: self, y
  type(QuaternionArray_T)             :: qres 

  type(IO_T)                          :: Message 
  integer(kind=irg)                   :: sz(2)

! test to make sure that both arrays have the same number of quaternions 
  if (self%n.ne.y%n) then 
    call Message%printError('quatarrayadd','input arrays must have the same number of quaternions')
  end if 

  if (self%s.ne.y%s) then 
    call Message%printError('quatarrayadd','input arrays must have the same precision')
  end if 

  qres%n = self%n
  qres%s = self%s
  qres%nthreads = self%nthreads

   if (self%s.eq.'s') then
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%q)) then 
      sz = shape(qres%q)
      if ( (sz(1).ne.4).or.(sz(2).ne.self%n) ) deallocate(qres%q)
    end if 
    allocate(qres%q(4,self%n))
    qres%q = self%q - y%q 
  else
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%qd)) then 
      sz = shape(qres%qd)
      if ( (sz(1).ne.4).or.(sz(2).ne.self%n) ) deallocate(qres%qd)
    end if 
    allocate(qres%qd(4,self%n))
    qres%qd = self%qd - y%qd 
  end if 

end function quatarraysubtract

!--------------------------------------------------------------------------
pure recursive function quatmult(self, y) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion multiplication   (single/double precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in) :: self, y
   !! input quaternions
  type(Quaternion_T)             :: qres 
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
recursive function quatarraymult(self, y) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! quaternion array multiplication   (single/double precision)

use mod_io
use omp_lib

IMPLICIT NONE

  class(QuaternionArray_T),intent(in) :: self, y
   !! input quaternion arrays
  type(QuaternionArray_T)             :: qres 
   !! output quaternion array

  type(IO_T)                          :: Message 
  integer(kind=irg)                   :: i, maxthreads, ompthreads, sz(2) 

! test to make sure that both arrays have the same number of quaternions 
  if (self%n.ne.y%n) then 
    call Message%printError('quatarrayadd','input arrays must have the same number of quaternions')
  end if 

  if (self%s.ne.y%s) then 
    call Message%printError('quatarrayadd','input arrays must have the same precision')
  end if 

  qres%n = self%n
  qres%s = self%s
  qres%nthreads = self%nthreads

! set the number of OpenMP threads
  maxthreads = OMP_GET_MAX_THREADS()
  if (self%nthreads.eq.0) then ! use the maximum number of threads unless self%n is less than that
    ompthreads = maxthreads
    if (self%n.lt.ompthreads) ompthreads = self%n 
  else
    ompthreads = self%nthreads
    if (self%nthreads.gt.maxthreads) ompthreads = maxthreads
  end if 
  call OMP_SET_NUM_THREADS(ompthreads)

  if (self%s.eq.'s') then 
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%q)) then 
      sz = shape(qres%q)
      if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%q)
    end if 
    allocate(qres%q(4,self%n))

!$OMP PARALLEL DEFAULT(shared)
!$OMP DO SCHEDULE(DYNAMIC,1) 
    do i=1,self%n
      qres%q(1:4,i) = (/ &
         self%q(1,i)*y%q(1,i) - self%q(2,i)*y%q(2,i) -          ( self%q(3,i)*y%q(3,i) + self%q(4,i)*y%q(4,i) ), &
         self%q(1,i)*y%q(2,i) + self%q(2,i)*y%q(1,i) + epsijk * ( self%q(3,i)*y%q(4,i) - self%q(4,i)*y%q(3,i) ), &
         self%q(1,i)*y%q(3,i) + self%q(3,i)*y%q(1,i) + epsijk * ( self%q(4,i)*y%q(2,i) - self%q(2,i)*y%q(4,i) ), &
         self%q(1,i)*y%q(4,i) + self%q(4,i)*y%q(1,i) + epsijk * ( self%q(2,i)*y%q(3,i) - self%q(3,i)*y%q(2,i) ) /)
    end do 
!$OMP END DO 
!$OMP END PARALLEL    
  else 
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%qd)) then 
      sz = shape(qres%qd)
      if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%qd)
    end if 
    allocate(qres%qd(4,self%n))

!$OMP PARALLEL DEFAULT(shared)
!$OMP DO SCHEDULE(DYNAMIC,1) 
    do i=1,self%n
    qres%qd(1:4,i) = (/ &
         self%qd(1,i)*y%qd(1,i) - self%qd(2,i)*y%qd(2,i) -           ( self%qd(3,i)*y%qd(3,i) + self%qd(4,i)*y%qd(4,i) ), &
         self%qd(1,i)*y%qd(2,i) + self%qd(2,i)*y%qd(1,i) + epsijkd * ( self%qd(3,i)*y%qd(4,i) - self%qd(4,i)*y%qd(3,i) ), &
         self%qd(1,i)*y%qd(3,i) + self%qd(3,i)*y%qd(1,i) + epsijkd * ( self%qd(4,i)*y%qd(2,i) - self%qd(2,i)*y%qd(4,i) ), &
         self%qd(1,i)*y%qd(4,i) + self%qd(4,i)*y%qd(1,i) + epsijkd * ( self%qd(2,i)*y%qd(3,i) - self%qd(3,i)*y%qd(2,i) ) /) 
    end do 
!$OMP END DO 
!$OMP END PARALLEL    
  end if
      
end function quatarraymult

!--------------------------------------------------------------------------
pure recursive function quatsmult(self, s) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! scalar quaternion multiplication   (single precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in)   :: self 
   !! input quaternion
  real(kind=sgl), INTENT(IN)       :: s
   !! scalar input 
  type(Quaternion_T)               :: qres 
   !! output quaternion

  qres%q = (/ s*self%q(1), s*self%q(2), s*self%q(3), s*self%q(4) /) 
  qres%s = 's'

end function quatsmult

!--------------------------------------------------------------------------
pure recursive function quatarraysmult(self, s) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! scalar quaternion array multiplication   (single precision)

IMPLICIT NONE

  class(QuaternionArray_T),intent(in)   :: self 
   !! input quaternion
  real(kind=sgl), INTENT(IN)            :: s
   !! scalar input 
  type(QuaternionArray_T)               :: qres 
   !! output quaternion

  integer(kind=irg)                     :: sz(2) 

! if the quaternion array is already allocated, check to make sure it has the right dimensions
  if (allocated(qres%q)) then 
    sz = shape(qres%q)
    if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%q)
  end if 
  allocate(qres%q(4,self%n))

  qres%q = s*self%q
  qres%s = self%s
  qres%n = self%n
  qres%nthreads = self%nthreads

end function quatarraysmult

!--------------------------------------------------------------------------
pure recursive function quatsmultd(self, s) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! scalar quaternion multiplication   (double precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in)   :: self 
   !! input quaternion
  real(kind=dbl), INTENT(IN)       :: s
   !! scalar input 
  type(Quaternion_T)               :: qres 
   !! output quaternion

  qres%qd = (/ s*self%qd(1), s*self%qd(2), s*self%qd(3), s*self%qd(4) /) 
  qres%s = 'd'

end function quatsmultd

!--------------------------------------------------------------------------
pure recursive function quatarraysmultd(self, s) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! scalar quaternion array multiplication   (double precision)

IMPLICIT NONE

  class(QuaternionArray_T),intent(in)   :: self 
   !! input quaternion
  real(kind=dbl), INTENT(IN)            :: s
   !! scalar input 
  type(QuaternionArray_T)               :: qres 
   !! output quaternion

  integer(kind=irg)                     :: sz(2) 

! if the quaternion array is already allocated, check to make sure it has the right dimensions
  if (allocated(qres%qd)) then 
    sz = shape(qres%qd)
    if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%qd)
  end if 
  allocate(qres%qd(4,self%n))

  qres%qd = s*self%qd
  qres%s = self%s
  qres%n = self%n
  qres%nthreads = self%nthreads

end function quatarraysmultd

!--------------------------------------------------------------------------
pure recursive function quatconjg(self) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion conjugation (extends intrinsic routine conjg)

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion
  type(Quaternion_T)                :: qres
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
pure recursive function quatarrayconjg(self) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! quaternion array conjugation (extends intrinsic routine conjg)

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)    :: self
   !! input quaternion
  type(QuaternionArray_T)                :: qres
   !! output quaternion

  integer(kind=irg)                      :: sz(2) 

  if (self%s.eq.'s') then 
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%q)) then 
      sz = shape(qres%q)
      if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%q)
    end if 
    allocate(qres%q(4,self%n))

    qres%q(1,:)   =  self%q(1,:)
    qres%q(2:4,:) = -self%q(2:4,:)
  else
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%qd)) then 
      sz = shape(qres%qd)
      if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%qd)
    end if 
    allocate(qres%qd(4,self%n))

    qres%qd(1,:)   =  self%qd(1,:)
    qres%qd(2:4,:) = -self%qd(2:4,:)
  end if 

  qres%s = self%s
  qres%n = self%n
  qres%nthreads = self%nthreads

end function quatarrayconjg

!--------------------------------------------------------------------------
pure recursive function quatnorm(self) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion norm (extends intrinsic routine abs)

IMPLICIT NONE 

  class(Quaternion_T),intent(in) :: self
   !! input quaternion
  real(kind=dbl)                 :: res
   !! output norm

  real(kind=sgl)                 :: n
  real(kind=dbl)                 :: nd, resd

  if (self%s.eq.'s') then
    n = self%q(1)**2 + self%q(2)**2 + self%q(3)**2 + self%q(4)**2
    resd = dsqrt( dble(n) )
    res = dble(sngl(resd))
  else
    nd = self%qd(1)**2 + self%qd(2)**2 + self%qd(3)**2 + self%qd(4)**2
    res = dsqrt( nd )
  end if 

end function quatnorm

!--------------------------------------------------------------------------
pure recursive function quatarraynorm(self) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! quaternion array norm (extends intrinsic routine abs)
  !!
  !! this routine requires the output array to be allocated in the calling program.

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in) :: self
   !! input quaternion
  real(kind=dbl)                      :: res(self%n)
   !! output norm

  real(kind=sgl),allocatable          :: n(:)
  real(kind=dbl),allocatable          :: nd(:), resd(:)

  if (self%s.eq.'s') then
    allocate(n(self%n), resd(self%n))
    n(:) = self%q(1,:)**2 + self%q(2,:)**2 + self%q(3,:)**2 + self%q(4,:)**2
    resd = dsqrt( dble(n) )
    res = dble(sngl(resd))
    deallocate(n, resd)
  else
    allocate(nd(self%n)) 
    nd(:) = self%qd(1,:)**2 + self%qd(2,:)**2 + self%qd(3,:)**2 + self%qd(4,:)**2
    res = dsqrt( nd )
    deallocate(nd)
  end if 

end function quatarraynorm

!--------------------------------------------------------------------------
recursive subroutine quatnormalize(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/07/20
  !!
  !! normalize the input quaternion 

IMPLICIT NONE 

  class(Quaternion_T),intent(inout) :: self
   !! input quaternion

  type(Quaternion_T)                :: q
  real(kind=sgl)                    :: n
  real(kind=dbl)                    :: nd

  if (self%s.eq.'s') then
    n = self%q(1)**2 + self%q(2)**2 + self%q(3)**2 + self%q(4)**2
    n = sqrt( n )
    q = self%quatsdiv(n)
    self%q = q%q
  else
    nd = self%qd(1)**2 + self%qd(2)**2 + self%qd(3)**2 + self%qd(4)**2
    nd = sqrt( nd )
    q = self%quatsdivd(nd)
    self%qd = q%qd
  end if 

end subroutine quatnormalize

!--------------------------------------------------------------------------
recursive subroutine quatarraynormalize(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/07/20
  !!
  !! normalize the input quaternion 

IMPLICIT NONE 

  class(QuaternionArray_T),intent(inout) :: self
   !! input quaternion

  real(kind=sgl),allocatable             :: n(:)
  real(kind=dbl),allocatable             :: nd(:)
  integer(kind=irg)                      :: i

  if (self%s.eq.'s') then
    allocate(n(self%n))
    n(:) = self%q(1,:)**2 + self%q(2,:)**2 + self%q(3,:)**2 + self%q(4,:)**2
    n = sqrt( n )
    do i=1,self%n
      self%q(:,i) = self%q(:,i)/n(i)
    end do
    deallocate(n)
  else
    allocate(nd(self%n))
    nd(:) = self%qd(1,:)**2 + self%qd(2,:)**2 + self%qd(3,:)**2 + self%qd(4,:)**2
    nd = sqrt( nd )
    do i=1,self%n
      self%qd(:,i) = self%qd(:,i)/nd(i)
    end do
    deallocate(nd)
  end if 

end subroutine quatarraynormalize

!--------------------------------------------------------------------------
recursive function quatdiv(self, y) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion division (single/double precision)

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self, y
   !! input quaternions
  type(Quaternion_T)                :: qres
   !! output quaternion 

  type(Quaternion_T)                :: p, cy
  real(kind=sgl)                    :: q
  real(kind=dbl)                    :: qd

  if (self%s.eq.'s') then 
      q = quatnorm(y)
      cy = quatconjg(y)
      p = quatsdiv( cy, q*q )
      qres = quatmult(self,p)
      qres%s = 's'
  else
      qd = quatnorm(y)
      cy = quatconjg(y)
      p = quatsdivd( cy, qd*qd )
      qres = quatmult(self,p)
      qres%s = 'd'
  end if 


end function quatdiv

!--------------------------------------------------------------------------
recursive function quatarraydiv(self, y) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion array division (single/double precision)

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)    :: self, y
   !! input quaternions
  type(QuaternionArray_T)                :: qres
   !! output quaternion 

  type(QuaternionArray_T)                :: p, cy
  real(kind=sgl),allocatable             :: q(:)
  real(kind=dbl),allocatable             :: qd(:)
  integer(kind=irg)                      :: i

  if (self%s.eq.'s') then 
      q = quatarraynorm(y)
      cy = quatarrayconjg(y)
      do i=1,self%n 
        cy%q(:,i) = cy%q(:,i) / q(i)**2
      end do
      qres = quatarraymult(self,cy)
  else
      qd = quatarraynorm(y)
      cy = quatarrayconjg(y)
      do i=1,self%n 
        cy%qd(:,i) = cy%qd(:,i) / qd(i)**2
      end do
      qres = quatarraymult(self,cy)
  end if 

end function quatarraydiv

!--------------------------------------------------------------------------
recursive function quatsdiv(self, s) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion division (single precision)

use mod_io 

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion (numerator)
  real(kind=sgl), INTENT(IN)        :: s            
   !! input quaternion (denominator)

  type(Quaternion_T)                :: qres
  type(IO_T)                        :: Message

  if (s.ne.0.0) then 
      qres%q = (/ self%q(1)/s, self%q(2)/s, self%q(3)/s, self%q(4)/s /) 
      qres%s = 's'
  else 
    call Message % printWarning('quatsdiv', (/ 'Attempting to divide quaternion by zero; skipping operation ...' /) )
    qres = self
  end if

end function quatsdiv

!--------------------------------------------------------------------------
recursive function quatarraysdiv(self, s) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! quaternion array scalar division (single precision)

use mod_io 

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)    :: self
   !! input quaternion (numerator)
  real(kind=sgl), INTENT(IN)             :: s            
   !! input quaternion (denominator)

  type(QuaternionArray_T)                :: qres
  type(IO_T)                             :: Message
  integer(kind=irg)                      :: sz(2)

  if (s.ne.0.0) then 
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%q)) then 
      sz = shape(qres%q)
      if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%q)
    end if 
    allocate(qres%q(4,self%n))

    qres%q = self%q/s
    qres%s = self%s
    qres%n = self%n
    qres%nthreads = self%nthreads
  else 
    call Message % printError('quatarraysdiv', 'Attempting to divide quaternion aray by zero' )
  end if

end function quatarraysdiv

!--------------------------------------------------------------------------
recursive function quatsdivd(self, s) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion division (doubgle precision)

use mod_io 

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion (numerator)
  real(kind=dbl), INTENT(IN)        :: s            
   !! input quaternion (denominator)

  type(Quaternion_T)                :: qres
  type(IO_T)                        :: Message

  if (s.ne.0.0) then 
    qres%qd = (/ self%qd(1)/s, self%qd(2)/s, self%qd(3)/s, self%qd(4)/s /) 
    qres%s = 'd'
  else 
    call Message % printWarning('quatsdivd', (/ 'Attempting to divide quaternion by zero; skipping operation ...' /) )
    qres = self
  end if

end function quatsdivd

!--------------------------------------------------------------------------
recursive function quatarraysdivd(self, s) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! quaternion array scalar division (double precision)

use mod_io 

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)    :: self
   !! input quaternion (numerator)
  real(kind=dbl), INTENT(IN)             :: s            
   !! input quaternion (denominator)

  type(QuaternionArray_T)                :: qres
  type(IO_T)                             :: Message
  integer(kind=irg)                      :: sz(2)

  if (s.ne.0.D0) then 
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%qd)) then 
      sz = shape(qres%qd)
      if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%qd)
    end if 
    allocate(qres%qd(4,self%n))

    qres%qd = self%qd/s
    qres%s = self%s
    qres%n = self%n
    qres%nthreads = self%nthreads
  else 
    call Message % printError('quatarraysdivd', 'Attempting to divide quaternion aray by zero' )
  end if

end function quatarraysdivd

!--------------------------------------------------------------------------
pure recursive function quatinnerproduct(self, y) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion inner product (single precision)

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self, y
   !! input quaternions
  real(kind=dbl)                    :: res
   !! inner product 

  if (self%s.eq.'s') then 
    res = dble(self%q(1) * y%q(1) + self%q(2) * y%q(2) + self%q(3) * y%q(3) + self%q(4) * y%q(4))
  else
    res = self%qd(1) * y%qd(1) + self%qd(2) * y%qd(2) + self%qd(3) * y%qd(3) + self%qd(4) * y%qd(4)
  end if

end function quatinnerproduct

!--------------------------------------------------------------------------
recursive function quatarrayinnerproduct(self, y) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! quaternion array inner product (single precision)
  !!
  !! calling program must allocate the output array

use mod_io 

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)    :: self, y
   !! input quaternions
  real(kind=dbl)                         :: res(self%n)
   !! inner product 

   type(IO_T)                            :: Message

! test to make sure that both arrays have the same number of quaternions 
  if (self%n.ne.y%n) then 
    call Message%printError('quatarrayinnerproduct','input arrays must have the same number of quaternions')
  end if 

  if (self%s.ne.y%s) then 
    call Message%printError('quatarrayinnerproduct','input arrays must have the same precision')
  end if 
  
  if (self%s.eq.'s') then 
    res(:) = dble(self%q(1,:) * y%q(1,:) + self%q(2,:) * y%q(2,:) + self%q(3,:) * y%q(3,:) + self%q(4,:) * y%q(4,:))
  else
    res(:) = self%qd(1,:) * y%qd(1,:) + self%qd(2,:) * y%qd(2,:) + self%qd(3,:) * y%qd(3,:) + self%qd(4,:) * y%qd(4,:)
  end if

end function quatarrayinnerproduct

!--------------------------------------------------------------------------!
pure recursive function quatangle(self, y) result(res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! interquaternion angle   (single/double precision)
  !! this only has meaning for a unit quaternion, so we test first and return -10000.0
  !! if either of the quaternions is not a unit quaternion.

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self, y
   !! input quaternions
  real(kind=dbl)                    :: res
   !! angle (radians)

  real(kind=sgl)                    :: q, nself, ny
  real(kind=dbl)                    :: qd, nselfd, nyd

  if (self%s.eq.'s') then 
      nself = self%quatnorm()
      ny = y%quatnorm()
      q = self%quat_innerproduct(y)
      res = dble(acos( q/(nself * ny) ))
  else 
      nselfd = self%quatnorm()
      nyd = y%quatnorm()
      qd = self%quat_innerproduct(y)
      res = dacos( qd/(nselfd * nyd) )
  end if

end function quatangle

!--------------------------------------------------------------------------!
recursive function quatarrayangle(self, y) result(res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! interquaternion angle for a pair of quaternion arrays (single/double precision)
  !! this only has meaning for a unit quaternion, so we test first and return -10000.0
  !! if either of the quaternions is not a unit quaternion.

use mod_io

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)     :: self, y
   !! input quaternions
  real(kind=dbl)                          :: res(self%n)
   !! angle (radians)

  type(IO_T)                              :: Message
  real(kind=sgl), allocatable             :: q(:), nself(:), ny(:)
  real(kind=dbl), allocatable             :: qd(:), nselfd(:), nyd(:)

! test to make sure that both arrays have the same number of quaternions 
  if (self%n.ne.y%n) then 
    call Message%printError('quatarrayangle','input arrays must have the same number of quaternions')
  end if 

  if (self%s.ne.y%s) then 
    call Message%printError('quatarrayangle','input arrays must have the same precision')
  end if 

  if (self%s.eq.'s') then 
      allocate(nself(self%n), ny(self%n), q(self%n)) 
      nself = self%quatarraynorm()
      ny = y%quatarraynorm()
      q = self%quatarrayinnerproduct(y)
      res = dble(acos( q/(nself * ny) ))
      deallocate(nself, ny, q)
  else 
      allocate(nselfd(self%n), nyd(self%n), qd(self%n)) 
      nselfd = self%quatarraynorm()
      nyd = y%quatarraynorm()
      qd = self%quatarrayinnerproduct(y)
      res = dacos( qd/(nselfd * nyd) )
      deallocate(nselfd, nyd, qd)
  end if

end function quatarrayangle


!--------------------------------------------------------------------------!
! pure recursive function quatLp(self, v) result (res)
recursive function quatLp(self, v) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! actively rotate a unit vector by a unit quaternion, L_p = p v p* (single precision)

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion
  real(kind=sgl),intent(in)         :: v(3)         
   !! input vector to be rotated 
  real(kind=sgl)                    :: res(3)
   !! output vector

  type(Quaternion_T)                :: qv, rqv, cq

  qv%q = (/ 0.0, v(1), v(2), v(3) /) 
  qv%s = 's'
  cq = quatconjg(self)
  rqv = quatmult(self, quatmult(qv, cq) )
  res(1:3) = rqv%q(2:4)

end function quatLp

!--------------------------------------------------------------------------!
recursive function quatarrayLp(self, v) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! actively rotate a single vector by an array of unit quaternions, L_p = p v p* (single precision)

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)    :: self
   !! input quaternion
  real(kind=sgl),intent(in)              :: v(3)         
   !! input vector to be rotated 
  real(kind=sgl)                         :: res(3,self%n)
   !! output vector

  type(Quaternion_T)                     :: qv, rqv, q, cq 
  integer(kind=irg)                      :: i

  qv%q = (/ 0.0, v(1), v(2), v(3) /) 
  qv%s = 's'
  do i=1,self%n
    q = extractfromQuaternionArray(self, i)
    cq = q%quatconjg()
    rqv = quatmult(q, quatmult(qv, cq) )
    res(1:3,i) = rqv%q(2:4)
  end do 

end function quatarrayLp

!--------------------------------------------------------------------------!
! pure recursive function quatLpd(self, v) result (res)
recursive function quatLpd(self, v) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/03/20
  !!
  !! actively rotate a unit vector by a unit quaternion, L_p = p v p* (double precision)

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion
  real(kind=dbl),intent(in)         :: v(3)         
   !! input vector to be rotated 
  real(kind=dbl)                    :: res(3)
   !! output vector

  type(Quaternion_T)                :: qv, rqv, cq

  qv%qd = (/ 0.D0, v(1), v(2), v(3) /) 
  qv%s = 'd'
  cq = quatconjg(self)
  rqv = quatmult(self, quatmult(qv, cq) )
  res(1:3) = rqv%qd(2:4)

end function quatLpd

!--------------------------------------------------------------------------!
recursive function quatArrayLpd(self, v) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/03/20
  !!
  !! actively rotate a unit vector by an array of unit quaternions, L_p = p v p* (double precision)

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)    :: self
   !! input quaternion
  real(kind=dbl),intent(in)              :: v(3)         
   !! input vector to be rotated 
  real(kind=dbl)                         :: res(3,self%n)
   !! output vector

  type(Quaternion_T)                     :: qv, rqv, q, cq 
  integer(kind=irg)                      :: i

  qv%qd = (/ 0.D0, v(1), v(2), v(3) /) 
  qv%s = 'd'
  do i=1,self%n
    q = extractfromQuaternionArray(self, i)
    cq = q%quatconjg()
    rqv = quatmult(q, quatmult(qv, cq) )
    res(1:3,i) = rqv%qd(2:4)
  end do 

end function quatArrayLpd

!--------------------------------------------------------------------------!
recursive function extractfromQuaternionArray(self, i) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! actively rotate a single vector by an array of unit quaternions, L_p = p v p* (single precision)

use mod_io 

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)   :: self
   !! input quaternion array
  integer(kind=irg), intent(in)         :: i 
   !! quaternion to be extracted 
  type(Quaternion_T)                    :: res
   !! extracted quaternion 

  type(IO_T)                            :: Message 

  if (i.le.self%n) then 
    res%s = self%s 
    if (self%s.eq.'s') then 
      res%q(:) = self%q(:,i) 
    else
      res%qd(:) = self%qd(:,i) 
    end if 
  else 
    call Message%printWarning('extractfromQuaternionArray: requested quaternion index larger than array size', &
                              (/'   ---> returning empty quaternion'/) )
    if (self%s.eq.'s') then 
      res = Quaternion_T( q = (/ 0.0, 0.0, 0.0, 0.0 /) )
    else
      res = Quaternion_T( qd = (/ 0.D0, 0.D0, 0.D0, 0.D0 /) )
    end if
  end if 

end function extractfromQuaternionArray

!--------------------------------------------------------------------------!
! pure recursive function quatslerp(self, qb, n) result(res)
recursive function quatslerp(self, qb, n) result(res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! return an array of interpolated quaternions

IMPLICIT NONE 

  class(Quaternion_T),intent(in)         :: self   ! = qa 
   !! input quaternion (start)
  class(Quaternion_T),intent(in)         :: qb
   !! input quaternion (end)
  integer(kind=irg),intent(in)           :: n            
   !! number of steps in the interpolation
  type(Quaternion_T)                     :: res(n)
   !! output interpolated quaternion list

  type(Quaternion_T)                     :: cqa
  real(kind=sgl)                         :: theta, phi, dphi, s
  real(kind=dbl)                         :: thetad, phid, dphid, sd
  integer(kind=irg)                      :: i

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
recursive function quatsequal(self, qb) result(res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion division (doubgle precision)

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self, qb
   !! input quaternions 
  logical                           :: res

  type(Quaternion_T)                :: diff

  real(kind=sgl)                    :: d, eps=1.0e-7
  real(kind=dbl)                    :: dd, epsd=1.0e-12

  res = .TRUE.
  diff = self - qb 

  if (self%s.eq.'s') then 
    d = maxval( abs( diff%q(:) ) )
    if (d.gt.eps) res = .FALSE.
  else 
    dd = maxval( abs( diff%qd(:) ) )
    if (dd.gt.epsd) res = .FALSE.
  end if

end function quatsequal













! this needs to be moved into another module (sampling module?)
!--------------------------------------------------------------------------!
recursive function quat_Marsaglia(seed) result(q)
!DEC$ ATTRIBUTES DLLEXPORT :: quatMarsaglia

use mod_rng 

IMPLICIT NONE

type(rng_t),INTENT(INOUT)           :: seed 
!f2py intent(in,out) ::  seed 
type(Quaternion_T)                  :: q 

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
