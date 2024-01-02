! ###################################################################
! Copyright (c) 2013-2024, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_dualquaternions
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! Dual Quaternion and Dual Quaternion Array arithmetic classes
  !!
  !! Dual Quaternions are defined with the scalar part in position 1, and the vector part in positions 2:4;
  !! the transaltion part is in positions 5-8.
  !!
  !! There are two class definitions in this file, one for single dual quaternions, the other for
  !! dual quaternion array operations (some using OpenMP threads). The program MODDualQuaternionsTest.f90
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
      procedure dualquatconjg
      procedure dualquatarrayconjg
    end interface conjg

    interface cabs
      procedure dualquatnorm
      procedure dualquatarraynorm
    end interface cabs

    ! interface operator(.eq.)
    !   procedure dualquatsequal
    ! end interface

! definition of the quaternion class
  type, public :: DualQuaternion_T
    !! DualQuaternion Class definition
    private
      real(kind=sgl), dimension(8) :: q
       !! single precision quaternion
      real(kind=dbl), dimension(8) :: qd
       !! double precision quaternion
      character(1)                 :: s
       !! precision indicator ('s' or 'd')
      real(kind=dbl)               :: sc

    contains
    private
! dual quaternion IO routines
      procedure, pass(self) :: dualquatprint
      procedure, pass(self) :: getdualquats
      procedure, pass(self) :: getdualquatd
      procedure, pass(self) :: setdualquats
      procedure, pass(self) :: setdualquatd
! dual quaternion arithmetic routines
      procedure, pass(self) :: dualquatflip
      procedure, pass(self) :: dualquatpos
      procedure, pass(self) :: dualquatadd
      procedure, pass(self) :: dualquatsubtract
      procedure, pass(self) :: dualquatmult
      procedure, pass(self) :: dualquatsmult
      procedure, pass(self) :: dualquatsmultd
      procedure, pass(self) :: dualquatdiv
      procedure, pass(self) :: dualquatsdiv
      procedure, pass(self) :: dualquatsdivd
      procedure, pass(self) :: dualquatconjg
      procedure, pass(self) :: dualquatnorm
      procedure, pass(self) :: dualquatnormalize
      procedure, pass(self) :: generatedualquat 
! quaternion-based transformations
      procedure, pass(self) :: dualquatLp
      procedure, pass(self) :: dualquatLpd
      ! procedure, pass(self) :: dualquatLp_vecarray
      ! procedure, pass(self) :: dualquatLpd_vecarray
! routines with two or more input quaternions
      procedure, pass(self) :: dualquatinnerproduct
      ! procedure, pass(self) :: quatslerp
! miscellaneous routines
      ! procedure, pass(self), public :: dualquatsequal

      generic, public :: dualquat_print => dualquatprint
      generic, public :: dualquat_flip => dualquatflip
      generic, public :: dualquat_pos => dualquatpos
      generic, public :: get_dualquats => getdualquats
      generic, public :: get_dualquatd => getdualquatd
      generic, public :: set_dualquats => setdualquats
      generic, public :: set_dualquatd => setdualquatd
      generic, public :: dualquat_norm => dualquatnorm
      generic, public :: operator(+) => dualquatadd
      generic, public :: operator(-) => dualquatsubtract
      generic, public :: operator(*) => dualquatmult
      generic, public :: operator(*) => dualquatsmult, dualquatsmultd
      generic, public :: operator(/) => dualquatdiv
      generic, public :: operator(/) => dualquatsdiv, dualquatsdivd
      generic, public :: dualquat_normalize => dualquatnormalize
      generic, public :: generate_dualquat => generatedualquat
      generic, public :: dualquat_Lp => dualquatLp, dualquatLpd
      ! generic, public :: dualquat_Lp_vecarray => dualquatLp_vecarray, dualquatLpd_vecarray
      generic, public :: dualquat_innerproduct => dualquatinnerproduct
      ! generic, public :: quat_slerp => quatslerp

  end type DualQuaternion_T

! next we define the dualquaternion array class
  type, public :: DualQuaternionArray_T
    !! DualQuaternion Class definition
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
! quaternion IO routines
      procedure, pass(self) :: dualquatarrayprint
! quaternion arithmetic routines
      procedure, pass(self) :: dualquatarrayadd
      procedure, pass(self) :: dualquatarraysubtract
      procedure, pass(self) :: dualquatarraymult
      procedure, pass(self) :: dualquatarraysmult
      procedure, pass(self) :: dualquatarraysmultd
      procedure, pass(self) :: dualquatarraydiv
      procedure, pass(self) :: dualquatarraysdiv
      procedure, pass(self) :: dualquatarrayconjg
      procedure, pass(self) :: dualquatarraynorm
      procedure, pass(self) :: dualquatarraynormalize
! quaternion-based transformations
      ! procedure, pass(self) :: dualquatarrayLp
      ! procedure, pass(self) :: dualquatarrayLpd
! routines with two or more input quaternion arrays
      procedure, pass(self) :: dualquatarrayinnerproduct
      ! procedure, pass(self) :: dualquatarrayangle
! miscellaneous routines
      ! procedure, pass(self) :: extractfromDualQuaternionArray
      ! procedure, pass(self) :: insertDualQuatintoArray
      ! procedure, pass(self) :: deleteArray_

! generics
      generic, public :: dualquat_print => dualquatarrayprint
      generic, public :: operator(+) => dualquatarrayadd
      generic, public :: operator(-) => dualquatarraysubtract
      generic, public :: operator(*) => dualquatarraymult
      generic, public :: operator(*) => dualquatarraysmult, dualquatarraysmultd
      generic, public :: operator(/) => dualquatarraydiv
      generic, public :: operator(/) => dualquatarraysdiv
      generic, public :: dualquat_normalize => dualquatarraynormalize
      ! generic, public :: dualquat_Lp => dualquatarrayLp, dualquatarrayLpd
      generic, public :: dualquat_innerproduct => dualquatarrayinnerproduct
      ! generic, public :: dualquat_angle => dualquatarrayangle
      ! generic, public :: getDualQuatfromArray => extractfromDualQuaternionArray
      ! generic, public :: insertDualQuatinArray => insertDualQuatintoArray
      ! generic, public :: deleteArray => deleteArray_

  end type DualQuaternionArray_T

! the constructor routines for these classes
  interface DualQuaternion_T
    module procedure DualQuaternion_constructor
  end interface DualQuaternion_T

  interface DualQuaternionArray_T
    module procedure DualQuaternionArray_constructor
  end interface DualQuaternionArray_T

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! We begin with the functions/subroutines that are public in the
! two classes and pair up functions for individual and quaternion
! arrays for easier module maintenance.
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
type(DualQuaternion_T) function DualQuaternion_constructor( q, qd ) result(DualQuat)
!DEC$ ATTRIBUTES DLLEXPORT :: DualQuaternion_constructor
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! constructor for the DualQuaternion Class

IMPLICIT NONE

  real(kind=sgl), INTENT(IN), OPTIONAL      :: q(8)
  real(kind=dbl), INTENT(IN), OPTIONAL      :: qd(8)

! fill in one or the other quaternion
  if ((.not.present(q)).and.(.not.present(qd))) then
    DualQuat % q = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
    DualQuat % qd = (/ 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /)
    DualQuat % s = 's'
  else
      if (present(q)) then
        DualQuat % q = q
        DualQuat % qd = (/ 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /)
        DualQuat % s = 's'
      end if

      if (present(qd)) then
        DualQuat % q = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
        DualQuat % qd = qd
        DualQuat % s = 'd'
      end if
  end if

end function DualQuaternion_constructor

!--------------------------------------------------------------------------
subroutine DualQuaternion_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: DualQuaternion_destructor
!! author: MDG
!! version: 1.0
!! date: 02/02/20
!!
!! destructor for the DualQuaternion_T Class

IMPLICIT NONE

type(DualQuaternion_T), INTENT(INOUT)     :: self

call reportDestructor('DualQuaternion_T')

end subroutine DualQuaternion_destructor

!--------------------------------------------------------------------------
type(DualQuaternionArray_T) function DualQuaternionArray_constructor( n, nthreads, q, qd, s ) result(DualQuatArray)
!DEC$ ATTRIBUTES DLLEXPORT :: DualQuaternionArray_constructor
  !! author: MDG
  !! version: 1.0
  !! date: 10/13/22
  !!
  !! constructor for the DualQuaternionArray Class
  !!
  !! either call with parameters n and s
  !! or with n and either one of q or qd

IMPLICIT NONE

  integer(kind=irg), INTENT(IN)             :: n
  integer(kind=irg), INTENT(IN), OPTIONAL   :: nthreads
  real(kind=sgl), INTENT(IN), OPTIONAL      :: q(8,n)
  real(kind=dbl), INTENT(IN), OPTIONAL      :: qd(8,n)
  character(1), INTENT(IN), OPTIONAL        :: s

! OpenMP threads
  DualQuatArray % nthreads = 0
  if (present(nthreads)) DualQuatArray % nthreads = nthreads

! are we declaring just an empty variable with no entries, but with a given precision ?
  if ( present(s) .and. (.not.present(q)) .and. (.not.present(qd)) ) then
    DualQuatArray % n = n
    DualQuatArray % s = s
    if (s.eq.'s') then
      allocate(DualQuatArray % q(8,n))
      DualQuatArray % q = 0.0
    else
      allocate(DualQuatArray % qd(8,n))
      DualQuatArray % qd = 0.D0
    end if
    return
  end if

! single precision
  if (present(q)) then
    allocate(DualQuatArray % q(8,n))
    DualQuatArray % n = n
    DualQuatArray % q = q
    DualQuatArray % s = 's'
  end if

! double precision
  if (present(qd)) then
    allocate(DualQuatArray % qd(8,n))
    DualQuatArray % n = n
    DualQuatArray % qd = qd
    DualQuatArray % s = 'd'
  end if

end function DualQuaternionArray_constructor

!--------------------------------------------------------------------------
subroutine DualQuaternionArray_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: DualQuaternionArray_destructor
!! author: MDG
!! version: 1.0
!! date: 02/02/20
!!
!! destructor for the DualQuaternionArray_T Class

IMPLICIT NONE

type(DualQuaternionArray_T), INTENT(INOUT)     :: self

call reportDestructor('DualQuaternionArray_T')

if (allocated(self%q)) deallocate(self%q)
if (allocated(self%qd)) deallocate(self%qd)

end subroutine DualQuaternionArray_destructor

!--------------------------------------------------------------------------
recursive subroutine dualquatprint(self)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatprint
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! print a dual quaternion

use mod_io

IMPLICIT NONE

  class(DualQuaternion_T),intent(in)    :: self
   !! input dual quaternion

  type(IO_T)                        :: Message

  if (self%s.eq.'s') then
    call Message % WriteValue('', self%q, 8, frm="('(',8f12.6,'); precision: '$)")
    call Message % WriteValue('',self%s)
  else
    call Message % WriteValue('', self%qd, 8, frm="('(',8f20.14,'); precision: '$)")
    call Message % WriteValue('',self%s)
  end if

end subroutine dualquatprint

!--------------------------------------------------------------------------
recursive function getdualquats(self) result(qs)
!DEC$ ATTRIBUTES DLLEXPORT :: getdualquats
  !! author: MDG
  !! version: 1.0
  !! date: 01/22/20
  !!
  !! return a dual quaternion

IMPLICIT NONE

class(DualQuaternion_T),intent(in)    :: self
 !! input quaternion
real(kind=sgl)                        :: qs(8)

qs = self%q

end function getdualquats

!--------------------------------------------------------------------------
recursive function getdualquatd(self) result(qd)
!DEC$ ATTRIBUTES DLLEXPORT :: getdualquatd
  !! author: MDG
  !! version: 1.0
  !! date: 01/22/20
  !!
  !! return a dual quaternion

IMPLICIT NONE

class(DualQuaternion_T),intent(in)    :: self
 !! input dual quaternion
real(kind=dbl)                        :: qd(8)

qd = self%qd

end function getdualquatd

!--------------------------------------------------------------------------
recursive subroutine setdualquats(self, qs)
!DEC$ ATTRIBUTES DLLEXPORT :: setdualquats
  !! author: MDG
  !! version: 1.0
  !! date: 02/18/20
  !!
  !! set a dual quaternion

IMPLICIT NONE

class(DualQuaternion_T),intent(inout)    :: self
real(kind=sgl),intent(in)                :: qs(8)
 !! input dual quaternion

self%q = qs
self%s = 's'

end subroutine setdualquats

!--------------------------------------------------------------------------
recursive subroutine setdualquatd(self, qd)
!DEC$ ATTRIBUTES DLLEXPORT :: setdualquatd
  !! author: MDG
  !! version: 1.0
  !! date: 01/22/20
  !!
  !! set a dual quaternion

IMPLICIT NONE

class(DualQuaternion_T),intent(inout)    :: self
real(kind=dbl),intent(in)                :: qd(8)
 !! input dual quaternion

self%qd = qd
self%s = 'd'

end subroutine setdualquatd

!--------------------------------------------------------------------------
recursive subroutine dualquatarrayprint(self, listN)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatarrayprint
  !! author: MDG 
  !! version: 1.0 
  !! date: 10/13/22
  !!
  !! print an array of dual quaternions

use mod_io

IMPLICIT NONE 

  class(DualQuaternionArray_T),intent(in)   :: self
   !! input dual quaternion 
  integer(kind=irg),INTENT(IN),OPTIONAL     :: listN

  type(IO_T)                                :: Message 
  integer(kind=irg)                         :: i, n

  if (present(listN)) then 
    n = listN
  else 
    n = self%n 
  end if
  if (self%s.eq.'s') then 
    do i=1,n
      call Message % WriteValue('', self%q(:,i), 8, frm="('(',8f12.6,')')")
    end do
  else 
    do i=1,n
      call Message % WriteValue('', self%qd(:,i), 8, frm="('(',8f20.14,')')")
    end do
  end if 

end subroutine dualquatarrayprint

!--------------------------------------------------------------------------
pure recursive subroutine dualquatflip(self)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatflip
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! change the sign of the complete dual quaternion

IMPLICIT NONE

class(DualQuaternion_T),intent(inout) :: self

if (self%s.eq.'s') then
  self%q = -self%q
else
  self%qd = -self%qd
end if

end subroutine dualquatflip

!--------------------------------------------------------------------------
pure recursive subroutine dualquatpos(self)
!DEC$ ATTRIBUTES DLLEXPORT :: quatpos
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! convert a quaternion to one with a positive scalar part, if it is negative

IMPLICIT NONE

class(DualQuaternion_T),intent(inout) :: self

if (self%s.eq.'s') then
  if (self%q(1).lt.0.0) self%q = -self%q
else
  if (self%qd(1).lt.0.D0) self%qd = -self%qd
end if

end subroutine dualquatpos

!--------------------------------------------------------------------------
pure recursive function dualquatadd(self, y) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatadd
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! dual quaternion addition (single/double precision)

IMPLICIT NONE

  class(DualQuaternion_T),intent(in) :: self, y
  type(DualQuaternion_T)             :: qres

  if (self%s.eq.'s') then
    qres%q = self%q + y%q
    qres%s = 's'
  else
    qres%qd = self%qd + y%qd
    qres%s = 'd'
  end if

end function dualquatadd

!--------------------------------------------------------------------------
recursive function dualquatarrayadd(self, y) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatarrayadd
  !! author: MDG
  !! version: 1.0
  !! date: 10/13/22
  !!
  !! dual quaternion array addition (single/double precision)

use mod_io

IMPLICIT NONE

  class(DualQuaternionArray_T),intent(in) :: self, y
  type(DualQuaternionArray_T)             :: qres

  type(IO_T)                              :: Message
  integer(kind=irg)                       :: sz(2)

! test to make sure that both arrays have the same number of quaternions
  if (self%n.ne.y%n) then
    call Message%printError('dualquatarrayadd','input arrays must have the same number of quaternions')
  end if

  if (self%s.ne.y%s) then
    call Message%printError('dualquatarrayadd','input arrays must have the same precision')
  end if

  qres%n = self%n
  qres%s = self%s
  qres%nthreads = self%nthreads

  if (self%s.eq.'s') then
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%q)) then
      sz = shape(qres%q)
      if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(qres%q)
    end if
    allocate(qres%q(8,self%n))
    qres%q = self%q + y%q
  else
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%qd)) then
      sz = shape(qres%qd)
      if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(qres%qd)
    end if
    allocate(qres%qd(8,self%n))
    qres%qd = self%qd + y%qd
  end if

end function dualquatarrayadd

!--------------------------------------------------------------------------
recursive function dualquatsubtract(self, y) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatsubtract
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! dual quaternion subtraction (single/double precision)

IMPLICIT NONE

  class(DualQuaternion_T),intent(in) :: self, y
  type(DualQuaternion_T)             :: qres

  if (self%s.eq.'s') then
    qres%q = self%q - y%q
    qres%s = 's'
  else
    qres%qd = self%qd - y%qd
    qres%s = 'd'
  end if

end function dualquatsubtract

!--------------------------------------------------------------------------
recursive function dualquatarraysubtract(self, y) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatarraysubtract
  !! author: MDG
  !! version: 1.0
  !! date: 10/13/22
  !!
  !! quaternion array subtraction (single/double precision)

use mod_io

IMPLICIT NONE

  class(DualQuaternionArray_T),intent(in) :: self, y
  type(DualQuaternionArray_T)             :: qres

  type(IO_T)                              :: Message
  integer(kind=irg)                       :: sz(2)

! test to make sure that both arrays have the same number of quaternions
  if (self%n.ne.y%n) then
    call Message%printError('dualquatarraysubtract','input arrays must have the same number of quaternions')
  end if

  if (self%s.ne.y%s) then
    call Message%printError('dualquatarraysubtract','input arrays must have the same precision')
  end if

  qres%n = self%n
  qres%s = self%s
  qres%nthreads = self%nthreads

   if (self%s.eq.'s') then
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%q)) then
      sz = shape(qres%q)
      if ( (sz(1).ne.8).or.(sz(2).ne.self%n) ) deallocate(qres%q)
    end if
    allocate(qres%q(8,self%n))
    qres%q = self%q - y%q
  else
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%qd)) then
      sz = shape(qres%qd)
      if ( (sz(1).ne.8).or.(sz(2).ne.self%n) ) deallocate(qres%qd)
    end if
    allocate(qres%qd(8,self%n))
    qres%qd = self%qd - y%qd
  end if

end function dualquatarraysubtract

!--------------------------------------------------------------------------
recursive function dualquatmult(self, y) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatmult
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! dual quaternion multiplication   (single/double precision)
  !! q1 q2 = qr1 qr2 +(qr1 qd2+qd1 qr2) epsilon

use mod_quaternions

IMPLICIT NONE

  class(DualQuaternion_T),intent(in) :: self, y
   !! input dual quaternions
  type(DualQuaternion_T)             :: qres
   !! output dual quaternion

!! intermediate quaternions
  type(Quaternion_T)                 :: qr1, qd1, qr2, qd2, qp

  if (self%s.eq.'s') then
    qr1 = Quaternion_T( q = self%q(1:4) )
    qd1 = Quaternion_T( q = self%q(5:8) )
    qr2 = Quaternion_T( q = y%q(1:4) )
    qd2 = Quaternion_T( q = y%q(5:8) )
!! the real part
    qp = qr1 * qr2
    qres%q(1:4) = qp%get_quats()
!! the dual part
    qp = qr1 * qd2 + qd1 * qr2
    qres%q(5:8) = qp%get_quats()
!! and set the type
    qres%s = 's'
  else
    qr1 = Quaternion_T( qd = self%qd(1:4) )
    qd1 = Quaternion_T( qd = self%qd(5:8) )
    qr2 = Quaternion_T( qd = y%qd(1:4) )
    qd2 = Quaternion_T( qd = y%qd(5:8) )
!! the real part
    qp = qr1 * qr2
    qres%qd(1:4) = qp%get_quatd()
!! the dual part
    qp = qr1 * qd2 + qd1 * qr2
    qres%qd(5:8) = qp%get_quatd()
!! and set the type
    qres%s = 'd'
  end if

end function dualquatmult

!--------------------------------------------------------------------------
recursive function dualquatarraymult(self, y) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatarraymult
  !! author: MDG
  !! version: 1.0
  !! date: 10/13/22
  !!
  !! dual quaternion array multiplication   (single/double precision)

use mod_io
use omp_lib
use mod_OMPsupport

IMPLICIT NONE

  class(DualQuaternionArray_T),intent(in) :: self, y
   !! input quaternion arrays
  type(DualQuaternionArray_T)             :: qres
   !! output quaternion array

  type(DualQuaternion_T)                  :: q1, q2, qp
  type(IO_T)                              :: Message
  integer(kind=irg)                       :: i, sz(2)

! test to make sure that both arrays have the same number of quaternions
  if (self%n.ne.y%n) then
    call Message%printError('dualquatarraymult','input arrays must have the same number of quaternions')
  end if

  if (self%s.ne.y%s) then
    call Message%printError('dualquatarraymult','input arrays must have the same precision')
  end if

  qres%n = self%n
  qres%s = self%s
  qres%nthreads = self%nthreads

! set the number of OpenMP threads
  call OMP_setNThreads(qres%nthreads)

  if (self%s.eq.'s') then
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%q)) then
      sz = shape(qres%q)
      if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(qres%q)
    end if
    allocate(qres%q(8,self%n))

!$OMP PARALLEL DEFAULT(shared)
!$OMP DO SCHEDULE(DYNAMIC,1)
    do i=1,self%n
      q1 = DualQuaternion_T( q = self%q(1:8,i) )
      q2 = DualQuaternion_T( q = y%q(1:8,i) )
      qp = q1 * q2
      qres%q(1:8,i) = qp%get_dualquats()
    end do
!$OMP END DO
!$OMP END PARALLEL
  else
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%qd)) then
      sz = shape(qres%qd)
      if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(qres%qd)
    end if
    allocate(qres%qd(8,self%n))

!$OMP PARALLEL DEFAULT(shared)
!$OMP DO SCHEDULE(DYNAMIC,1)
    do i=1,self%n
      q1 = DualQuaternion_T( qd = self%qd(1:8,i) )
      q2 = DualQuaternion_T( qd = y%qd(1:8,i) )
      qp = q1 * q2
      qres%qd(1:8,i) = qp%get_dualquatd()
    end do
!$OMP END DO
!$OMP END PARALLEL
  end if

end function dualquatarraymult

!--------------------------------------------------------------------------
pure recursive function dualquatsmult(self, s) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatsmult
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! scalar dual quaternion multiplication   (single precision)

IMPLICIT NONE

  class(DualQuaternion_T),intent(in)   :: self
   !! input dual quaternion
  real(kind=sgl), INTENT(IN)           :: s
   !! scalar input
  type(DualQuaternion_T)               :: qres
   !! output dual quaternion

  qres%q(1:8) = s*self%q(1:8)
  qres%s = 's'

end function dualquatsmult

!--------------------------------------------------------------------------
pure recursive function dualquatarraysmult(self, s) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatarraysmult
  !! author: MDG
  !! version: 1.0
  !! date: 10/13/22
  !!
  !! scalar dual quaternion array multiplication   (single precision)

IMPLICIT NONE

  class(DualQuaternionArray_T),intent(in)   :: self
   !! input quaternion
  real(kind=sgl), INTENT(IN)                :: s
   !! scalar input
  type(DualQuaternionArray_T)               :: qres
   !! output quaternion

  integer(kind=irg)                         :: sz(2)

! if the quaternion array is already allocated, check to make sure it has the right dimensions
  if (allocated(qres%q)) then
    sz = shape(qres%q)
    if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(qres%q)
  end if
  allocate(qres%q(8,self%n))

  qres%q = s*self%q
  qres%s = self%s
  qres%n = self%n
  qres%nthreads = self%nthreads

end function dualquatarraysmult

!--------------------------------------------------------------------------
pure recursive function dualquatsmultd(self, s) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatsmultd
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! scalar dual quaternion multiplication   (double precision)

IMPLICIT NONE

  class(DualQuaternion_T),intent(in)   :: self
   !! input quaternion
  real(kind=dbl), INTENT(IN)           :: s
   !! scalar input
  type(DualQuaternion_T)               :: qres
   !! output quaternion

  qres%qd(1:8) = s*self%qd(1:8)
  qres%s = 'd'

end function dualquatsmultd

!--------------------------------------------------------------------------
pure recursive function dualquatarraysmultd(self, s) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatarraysmultd
  !! author: MDG
  !! version: 1.0
  !! date: 10/13/22
  !!
  !! scalar dual quaternion array multiplication   (double precision)

IMPLICIT NONE

  class(DualQuaternionArray_T),intent(in)   :: self
   !! input dual quaternion array
  real(kind=dbl), INTENT(IN)                :: s
   !! scalar input
  type(DualQuaternionArray_T)               :: qres
   !! output dual quaternion array

  integer(kind=irg)                         :: sz(2)

! if the quaternion array is already allocated, check to make sure it has the right dimensions
  if (allocated(qres%qd)) then
    sz = shape(qres%qd)
    if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(qres%qd)
  end if
  allocate(qres%qd(8,self%n))

  qres%qd = s*self%qd
  qres%s = self%s
  qres%n = self%n
  qres%nthreads = self%nthreads

end function dualquatarraysmultd

!--------------------------------------------------------------------------
pure recursive function dualquatconjg(self) result (qres)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatconjg
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! dual quaternion conjugation (extends intrinsic routine conjg)
  !! There are three possible definitions for the conjugate of a dual quaternion:
  !!
  !! ++++----
  !! +---+---
  !! +----+++
  !! 
  !! this routine uses the last one, indicated with a diamond in section 5 of the following paper:
  !!
  !! <https://faculty.sites.iastate.edu/jia/files/inline-files/dual-quaternion.pdf>


IMPLICIT NONE

  class(DualQuaternion_T),intent(in)    :: self
   !! input dual quaternion
  type(DualQuaternion_T)                :: qres
   !! output dual quaternion

  if (self%s.eq.'s') then
    qres%q = (/ self%q(1), -self%q(2), -self%q(3), -self%q(4), -self%q(5), self%q(6), self%q(7), self%q(8) /)
    qres%s = 's'
  else
    qres%qd = (/ self%qd(1), -self%qd(2), -self%qd(3), -self%qd(4), -self%qd(5), self%qd(6), self%qd(7), self%qd(8) /)
    qres%s = 'd'
  end if

end function dualquatconjg

!--------------------------------------------------------------------------
pure recursive function dualquatarrayconjg(self) result (qres)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatarrayconjg
  !! author: MDG
  !! version: 1.0
  !! date: 10/13/22
  !!
  !! dual quaternion array conjugation (extends intrinsic routine conjg)

IMPLICIT NONE

  class(DualQuaternionArray_T),intent(in)    :: self
   !! input dual quaternion array
  type(DualQuaternionArray_T)                :: qres
   !! output dual quaternion array

  integer(kind=irg)                          :: sz(2)

  if (self%s.eq.'s') then
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%q)) then
      sz = shape(qres%q)
      if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(qres%q)
    end if
    allocate(qres%q(8,self%n))

    qres%q(1,:)   =  self%q(1,:)
    qres%q(2:4,:) = -self%q(2:4,:)
    qres%q(5,:)   =  self%q(5,:)
    qres%q(6:8,:) = -self%q(6:8,:)
  else
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%qd)) then
      sz = shape(qres%qd)
      if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(qres%qd)
    end if
    allocate(qres%qd(8,self%n))

    qres%qd(1,:)   =  self%qd(1,:)
    qres%qd(2:4,:) = -self%qd(2:4,:)
    qres%qd(5,:)   =  self%qd(5,:)
    qres%qd(6:8,:) = -self%qd(6:8,:)
  end if

  qres%s = self%s
  qres%n = self%n
  qres%nthreads = self%nthreads

end function dualquatarrayconjg

!--------------------------------------------------------------------------
pure recursive function dualquatnorm(self) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatnorm
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! dual quaternion norm (extends intrinsic routine abs)

IMPLICIT NONE

  class(DualQuaternion_T),intent(in) :: self
   !! input quaternion
  real(kind=dbl)                     :: res
   !! output norm

  real(kind=sgl)                     :: n
  real(kind=dbl)                     :: nd, resd

  if (self%s.eq.'s') then
    n = sum(self%q(1:8)**2)
    resd = dsqrt( dble(n) )
    res = dble(sngl(resd))
  else
    nd = sum(self%qd(1:8)**2)
    res = dsqrt( nd )
  end if

end function dualquatnorm

!--------------------------------------------------------------------------
pure recursive function dualquatarraynorm(self) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatarraynorm
  !! author: MDG
  !! version: 1.0
  !! date: 10/13/22
  !!
  !! dual quaternion array norm (extends intrinsic routine abs)
  !!
  !! this routine requires the output array to be allocated in the calling program.

IMPLICIT NONE

  class(DualQuaternionArray_T),intent(in) :: self
   !! input quaternion
  real(kind=dbl)                          :: res(self%n)
   !! output norm

  real(kind=sgl),allocatable              :: n(:)
  real(kind=dbl),allocatable              :: nd(:), resd(:)
  integer(kind=irg)                       :: i

  if (self%s.eq.'s') then
    allocate(n(self%n), resd(self%n))
    do i=1,self%n
      n(i) = sum(self%q(1:8,i)**2)
    end do
    resd = dsqrt( dble(n) )
    res = dble(sngl(resd))
    deallocate(n, resd)
  else
    allocate(nd(self%n))
    do i=1,self%n
      nd(i) = sum(self%qd(1:8,i)**2)
    end do
    res = dsqrt( nd )
    deallocate(nd)
  end if

end function dualquatarraynorm

!--------------------------------------------------------------------------
recursive subroutine dualquatnormalize(self)
!DEC$ ATTRIBUTES DLLEXPORT :: quatnormalize
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! normalize the input dual quaternion

IMPLICIT NONE

  class(DualQuaternion_T),intent(inout) :: self
   !! input quaternion

  type(DualQuaternion_T)                :: q
  real(kind=sgl)                        :: n
  real(kind=dbl)                        :: nd

  if (self%s.eq.'s') then
    n = sum(self%q(1:8)**2) 
    n = sqrt( n )
    q = self%dualquatsdiv(n)
    self%q = q%q
  else
    nd = sum(self%qd(1:8)**2)
    nd = sqrt( nd )
    q = self%dualquatsdivd(nd)
    self%qd = q%qd
  end if

end subroutine dualquatnormalize

!--------------------------------------------------------------------------
recursive subroutine dualquatarraynormalize(self)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatarraynormalize
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! normalize the input dual quaternion array

IMPLICIT NONE

  class(DualQuaternionArray_T),intent(inout) :: self
   !! input quaternion

  real(kind=sgl),allocatable                 :: n(:)
  real(kind=dbl),allocatable                 :: nd(:)
  integer(kind=irg)                          :: i

  if (self%s.eq.'s') then
    allocate(n(self%n))
    do i=1,self%n
      n(i) = sum(self%q(1:8,i)**2)
    end do 
    n = sqrt( n )
    do i=1,self%n
      self%q(:,i) = self%q(:,i)/n(i)
    end do
    deallocate(n)
  else
    allocate(nd(self%n))
    do i=1,self%n
      nd(i) = sum(self%qd(1:8,i)**2)
    end do 
    nd = sqrt( nd )
    do i=1,self%n
      self%qd(:,i) = self%qd(:,i)/nd(i)
    end do
    deallocate(nd)
  end if

end subroutine dualquatarraynormalize

!--------------------------------------------------------------------------
recursive function dualquatdiv(self, y) result (qres)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatdiv
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! dual quaternion division (single/double precision)

IMPLICIT NONE

  class(DualQuaternion_T),intent(in)    :: self, y
   !! input quaternions
  type(DualQuaternion_T)                :: qres
   !! output quaternion

  type(DualQuaternion_T)                :: p, cy
  real(kind=sgl)                        :: q
  real(kind=dbl)                        :: qd

  if (self%s.eq.'s') then
      q = dualquatnorm(y)
      cy = dualquatconjg(y)
      p = dualquatsdiv( cy, q*q )
      qres = dualquatmult(self,p)
      qres%s = 's'
  else
      qd = dualquatnorm(y)
      cy = dualquatconjg(y)
      p = dualquatsdivd( cy, qd*qd )
      qres = dualquatmult(self,p)
      qres%s = 'd'
  end if


end function dualquatdiv

!--------------------------------------------------------------------------
recursive function dualquatarraydiv(self, y) result (qres)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatarraydiv
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! dual quaternion array division (single/double precision)

IMPLICIT NONE

  class(DualQuaternionArray_T),intent(in)    :: self, y
   !! input quaternions
  type(DualQuaternionArray_T)                :: qres
   !! output quaternion

  type(DualQuaternionArray_T)                :: p, cy
  real(kind=sgl),allocatable                 :: q(:)
  real(kind=dbl),allocatable                 :: qd(:)
  integer(kind=irg)                          :: i

  if (self%s.eq.'s') then
      q = dualquatarraynorm(y)
      cy = dualquatarrayconjg(y)
      do i=1,self%n
        cy%q(:,i) = cy%q(:,i) / q(i)**2
      end do
      qres = dualquatarraymult(self,cy)
  else
      qd = dualquatarraynorm(y)
      cy = dualquatarrayconjg(y)
      do i=1,self%n
        cy%qd(:,i) = cy%qd(:,i) / qd(i)**2
      end do
      qres = dualquatarraymult(self,cy)
  end if

end function dualquatarraydiv

!--------------------------------------------------------------------------
recursive function dualquatsdiv(self, s) result (qres)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatsdiv
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! dual quaternion division (single precision)

use mod_io

IMPLICIT NONE

  class(DualQuaternion_T),intent(in)    :: self
   !! input quaternion (numerator)
  real(kind=sgl), INTENT(IN)            :: s
   !! input quaternion (denominator)

  type(DualQuaternion_T)                :: qres
  type(IO_T)                            :: Message

  if (s.ne.0.0) then
      qres%q = self%q/s
      qres%s = 's'
  else
    call Message % printWarning('dualquatsdiv', (/ 'Attempting to divide dual quaternion by zero; skipping operation ...' /) )
    qres = self
  end if

end function dualquatsdiv

!--------------------------------------------------------------------------
recursive function dualquatarraysdiv(self, s) result (qres)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatarraysdiv
  !! author: MDG
  !! version: 1.0
  !! date: 10/13/22
  !!
  !! dual quaternion array scalar division (single precision)

use mod_io

IMPLICIT NONE

  class(DualQuaternionArray_T),intent(in)    :: self
   !! input quaternion (numerator)
  real(kind=sgl), INTENT(IN)                 :: s
   !! input quaternion (denominator)

  type(DualQuaternionArray_T)                :: qres
  type(IO_T)                                 :: Message
  integer(kind=irg)                          :: sz(2)

  if (s.ne.0.0) then
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%q)) then
      sz = shape(qres%q)
      if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(qres%q)
    end if
    allocate(qres%q(8,self%n))

    qres%q = self%q/s
    qres%s = self%s
    qres%n = self%n
    qres%nthreads = self%nthreads
  else
    call Message % printError('dualquatarraysdiv', 'Attempting to divide dual quaternion aray by zero' )
  end if

end function dualquatarraysdiv

!--------------------------------------------------------------------------
recursive function dualquatsdivd(self, s) result (qres)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatsdivd
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! dual quaternion division (double precision)

use mod_io

IMPLICIT NONE

  class(DualQuaternion_T),intent(in)    :: self
   !! input quaternion (numerator)
  real(kind=dbl), INTENT(IN)            :: s
   !! input quaternion (denominator)

  type(DualQuaternion_T)                :: qres
  type(IO_T)                            :: Message

  if (s.ne.0.0) then
    qres%qd = self%qd/s
    qres%s = 'd'
  else
    call Message % printWarning('dualquatsdivd', (/ 'Attempting to divide dual quaternion by zero; skipping operation ...' /) )
    qres = self
  end if

end function dualquatsdivd

!--------------------------------------------------------------------------
recursive function dualquatarraysdivd(self, s) result (qres)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatarraysdivd
  !! author: MDG
  !! version: 1.0
  !! date: 10/13/22
  !!
  !! dual quaternion array scalar division (double precision)

use mod_io

IMPLICIT NONE

  class(DualQuaternionArray_T),intent(in)    :: self
   !! input quaternion (numerator)
  real(kind=dbl), INTENT(IN)                 :: s
   !! input quaternion (denominator)

  type(DualQuaternionArray_T)                :: qres
  type(IO_T)                                 :: Message
  integer(kind=irg)                          :: sz(2)

  if (s.ne.0.D0) then
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%qd)) then
      sz = shape(qres%qd)
      if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(qres%qd)
    end if
    allocate(qres%qd(8,self%n))

    qres%qd = self%qd/s
    qres%s = self%s
    qres%n = self%n
    qres%nthreads = self%nthreads
  else
    call Message % printError('dualquatarraysdivd', 'Attempting to divide dual quaternion aray by zero' )
  end if

end function dualquatarraysdivd

!--------------------------------------------------------------------------
pure recursive function dualquatinnerproduct(self, y) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatinnerproduct
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! dual quaternion inner product (single precision)

IMPLICIT NONE

  class(DualQuaternion_T),intent(in)    :: self, y
   !! input quaternions
  real(kind=dbl)                        :: res
   !! inner product

  if (self%s.eq.'s') then
    res = self%q(1) * y%q(1) + self%q(2) * y%q(2) + self%q(3) * y%q(3) + self%q(4) * y%q(4)
    res = res + self%q(5) * y%q(5) + self%q(6) * y%q(6) + self%q(7) * y%q(7) + self%q(8) * y%q(8)
  else
    res = self%qd(1) * y%qd(1) + self%qd(2) * y%qd(2) + self%qd(3) * y%qd(3) + self%qd(4) * y%qd(4)
    res = res + self%qd(5) * y%qd(5) + self%qd(6) * y%qd(6) + self%qd(7) * y%qd(7) + self%qd(8) * y%qd(8)
  end if

end function dualquatinnerproduct

!--------------------------------------------------------------------------
recursive function dualquatarrayinnerproduct(self, y) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatarrayinnerproduct
  !! author: MDG
  !! version: 1.0
  !! date: 10/13/22
  !!
  !! dual quaternion array inner product (single precision)
  !!
  !! calling program must allocate the output array

use mod_io

IMPLICIT NONE

  class(DualQuaternionArray_T),intent(in)    :: self, y
   !! input quaternions
  real(kind=dbl)                             :: res(self%n)
   !! inner product

   type(IO_T)                                :: Message
   integer(kind=irg)                         :: i

! test to make sure that both arrays have the same number of quaternions
  if (self%n.ne.y%n) then
    call Message%printError('dualquatarrayinnerproduct','input arrays must have the same number of quaternions')
  end if

  if (self%s.ne.y%s) then
    call Message%printError('dualquatarrayinnerproduct','input arrays must have the same precision')
  end if

  if (self%s.eq.'s') then
    do i=1,self%n 
      res(i) = sum(self%q(1:8,i) * y%q(1:8,i))
    end do
  else
    do i=1,self%n 
      res(i) = sum(self%qd(1:8,i) * y%qd(1:8,i))
    end do
  end if

end function dualquatarrayinnerproduct

!--------------------------------------------------------------------------!
recursive subroutine generatedualquat(self, ang, ax, tr, sc) 
!DEC$ ATTRIBUTES DLLEXPORT :: generatedualquat
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! create a dual quaternion that represents a rotation followed by a translation
  !! assumes angle in radians, unit axis vector, and arbitrary translation vector
  !! the isotropic scaling factor is given by sc 

use mod_quaternions

IMPLICIT NONE

  class(DualQuaternion_T),intent(inout) :: self
  real(kind=dbl),intent(in)             :: ang
  real(kind=dbl),intent(in)             :: ax(3)
  real(kind=dbl),intent(in)             :: tr(3)
  real(kind=dbl),intent(in)             :: sc

  type(DualQuaternion_T)                :: qr, qt, qp
  real(kind=sgl)                        :: ss 
  real(kind=dbl)                        :: sd, nax(3) 

! normalize the axis vector
  nax = ax/sqrt(sum(ax**2))
  self%sc = sc

  if (self%s.eq.'s') then 
    ss = sin(real(ang)*0.5)
    qr = DualQuaternion_T( q = (/ cos(real(ang)*0.5), ss*real(nax(1)), ss*real(nax(2)), ss*real(nax(3)), 0.0, 0.0, 0.0, 0.0 /) )
    qt = DualQuaternion_T( q = (/ 1.0, 0.0, 0.0, 0.0, 0.0, 0.5*real(tr(1)), 0.5*real(tr(2)), 0.5*real(tr(3)) /) )
    qp = qt * qr 
    self%q(1:8) = qp%get_dualquats()
  else
    sd = sin(ang*0.5D0)
    qr = DualQuaternion_T( qd = (/ cos(ang*0.5D0), sd*nax(1), sd*nax(2), sd*nax(3), 0.D0, 0.D0, 0.D0, 0.D0 /) )
    qt = DualQuaternion_T( qd = (/ 1.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.5D0*tr(1), 0.5D0*tr(2), 0.5D0*tr(3) /) )
    qp = qt * qr 
    self%qd(1:8) = qp%get_dualquatd()
  end if

end subroutine generatedualquat

!--------------------------------------------------------------------------!
recursive function dualquatLp(self, v) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatLp
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! actively rotate and translate a unit vector by a dual quaternion, L_p = p v p* (single precision)

IMPLICIT NONE

  class(DualQuaternion_T),intent(in)    :: self
   !! input quaternion
  real(kind=sgl),intent(in)             :: v(3)
   !! input vector to be rotated
  real(kind=sgl)                        :: res(3)
   !! output vector

  type(DualQuaternion_T)                :: qv, rqv, cq

  qv%q = (/ 1.0,  0.0, 0.0, 0.0, 0.0, real(self%sc)*v(1), real(self%sc)*v(2), real(self%sc)*v(3) /)
  qv%s = 's'
  cq = dualquatconjg(self)
  rqv = dualquatmult(self, dualquatmult(qv, cq) )
  res(1:3) = rqv%q(6:8)

end function dualquatLp

!--------------------------------------------------------------------------!
recursive function dualquatLpd(self, v) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: dualquatLpd
  !! author: MDG
  !! version: 1.0
  !! date: 10/11/22
  !!
  !! actively rotate and translate a unit vector by a dual quaternion, L_p = p v p* (double precision)

IMPLICIT NONE

  class(DualQuaternion_T),intent(in)    :: self
   !! input quaternion
  real(kind=dbl),intent(in)             :: v(3)
   !! input vector to be rotated
  real(kind=dbl)                        :: res(3)
   !! output vector

  type(DualQuaternion_T)                :: qv, rqv, cq

  qv%qd = (/ 1.D0, 0.D0, 0.D0, 0.D0, 0.D0, self%sc*v(1), self%sc*v(2), self%sc*v(3) /)
  qv%s = 'd'
  cq = dualquatconjg(self)
  rqv = dualquatmult(self, dualquatmult(qv, cq) )
  res(1:3) = rqv%qd(6:8)

end function dualquatLpd


! !--------------------------------------------------------------------------!
! recursive function quatLp_vecarray(self, N, v) result (res)
! !DEC$ ATTRIBUTES DLLEXPORT :: quatLp_vecarray
!   !! author: MDG
!   !! version: 1.0
!   !! date: 06/02/21
!   !!
!   !! actively rotate an array of unit vectors by a unit quaternion, L_p = p v p* (single precision)

! IMPLICIT NONE

!   class(Quaternion_T),intent(in)    :: self
!    !! input quaternion
!   integer(kind=irg),intent(in)      :: N
!   real(kind=sgl),intent(in)         :: v(3, N)
!    !! input vector to be rotated
!   real(kind=sgl)                    :: res(3, N)
!    !! output vector

!   type(Quaternion_T)                :: qv, rqv, cq
!   integer(kind=irg)                 :: i

!   do i=1,N 
!     qv%q = (/ 0.0, v(1, i), v(2, i), v(3, i) /)
!     qv%s = 's'
!     cq = quatconjg(self)
!     rqv = quatmult(self, quatmult(qv, cq) )
!     res(1:3, i) = rqv%q(2:4)
!   end do

! end function quatLp_vecarray

! !--------------------------------------------------------------------------!
! recursive function quatarrayLp(self, v) result (res)
! !DEC$ ATTRIBUTES DLLEXPORT :: quatarrayLp
!   !! author: MDG
!   !! version: 1.0
!   !! date: 10/11/22
!   !!
!   !! actively rotate a single vector by an array of unit quaternions, L_p = p v p* (single precision)

! IMPLICIT NONE

!   class(QuaternionArray_T),intent(in)    :: self
!    !! input quaternion
!   real(kind=sgl),intent(in)              :: v(3)
!    !! input vector to be rotated
!   real(kind=sgl)                         :: res(3,self%n)
!    !! output vector

!   type(Quaternion_T)                     :: qv, rqv, q, cq
!   integer(kind=irg)                      :: i

!   qv%q = (/ 0.0, v(1), v(2), v(3) /)
!   qv%s = 's'
!   do i=1,self%n
!     q = extractfromQuaternionArray(self, i)
!     cq = q%quatconjg()
!     rqv = quatmult(q, quatmult(qv, cq) )
!     res(1:3,i) = rqv%q(2:4)
!   end do

! end function quatarrayLp

! !--------------------------------------------------------------------------!
! ! pure recursive function quatLpd(self, v) result (res)
! recursive function quatLpd(self, v) result (res)
! !DEC$ ATTRIBUTES DLLEXPORT :: quatLpd
!   !! author: MDG
!   !! version: 1.0
!   !! date: 01/03/20
!   !!
!   !! actively rotate a unit vector by a unit quaternion, L_p = p v p* (double precision)

! IMPLICIT NONE

!   class(Quaternion_T),intent(in)    :: self
!    !! input quaternion
!   real(kind=dbl),intent(in)         :: v(3)
!    !! input vector to be rotated
!   real(kind=dbl)                    :: res(3)
!    !! output vector

!   type(Quaternion_T)                :: qv, rqv, cq

!   qv%qd = (/ 0.D0, v(1), v(2), v(3) /)
!   qv%s = 'd'
!   cq = quatconjg(self)
!   rqv = quatmult(self, quatmult(qv, cq) )
!   res(1:3) = rqv%qd(2:4)

! end function quatLpd

! !--------------------------------------------------------------------------!
! recursive function quatLpd_vecarray(self, N, v) result (res)
! !DEC$ ATTRIBUTES DLLEXPORT :: quatLpd_vecarray
!   !! author: MDG
!   !! version: 1.0
!   !! date: 06/02/21
!   !!
!   !! actively rotate an array of unit vectors by a unit quaternion, L_p = p v p* (double precision)

! IMPLICIT NONE

!   class(Quaternion_T),intent(in)    :: self
!    !! input quaternion
!   integer(kind=irg),intent(in)      :: N
!   real(kind=dbl),intent(in)         :: v(3, N)
!    !! input vector to be rotated
!   real(kind=dbl)                    :: res(3, N)
!    !! output vector

!   type(Quaternion_T)                :: qv, rqv, cq
!   integer(kind=irg)                 :: i

!   do i=1,N 
!     qv%qd = (/ 0.D0, v(1, i), v(2, i), v(3, i) /)
!     qv%s = 'd'
!     cq = quatconjg(self)
!     rqv = quatmult(self, quatmult(qv, cq) )
!     res(1:3, i) = rqv%qd(2:4)
!   end do

! end function quatLpd_vecarray



! !--------------------------------------------------------------------------!
! recursive function quatArrayLpd(self, v) result (res)
! !DEC$ ATTRIBUTES DLLEXPORT :: quatArrayLpd
!   !! author: MDG
!   !! version: 1.0
!   !! date: 01/03/20
!   !!
!   !! actively rotate a unit vector by an array of unit quaternions, L_p = p v p* (double precision)

! IMPLICIT NONE

!   class(QuaternionArray_T),intent(in)    :: self
!    !! input quaternion
!   real(kind=dbl),intent(in)              :: v(3)
!    !! input vector to be rotated
!   real(kind=dbl)                         :: res(3,self%n)
!    !! output vector

!   type(Quaternion_T)                     :: qv, rqv, q, cq
!   integer(kind=irg)                      :: i

!   qv%qd = (/ 0.D0, v(1), v(2), v(3) /)
!   qv%s = 'd'
!   do i=1,self%n
!     q = extractfromQuaternionArray(self, i)
!     cq = q%quatconjg()
!     rqv = quatmult(q, quatmult(qv, cq) )
!     res(1:3,i) = rqv%qd(2:4)
!   end do

! end function quatArrayLpd

! !--------------------------------------------------------------------------!
! recursive function extractfromQuaternionArray(self, i) result (res)
! !DEC$ ATTRIBUTES DLLEXPORT :: extractfromQuaternionArray
!   !! author: MDG
!   !! version: 1.0
!   !! date: 10/13/22
!   !!
!   !! extract a quaternion from an array of quaternions

! use mod_io

! IMPLICIT NONE

!   class(QuaternionArray_T),intent(in)   :: self
!    !! input quaternion array
!   integer(kind=irg), intent(in)         :: i
!    !! quaternion to be extracted
!   type(Quaternion_T)                    :: res
!    !! extracted quaternion

!   type(IO_T)                            :: Message

!   if (i.le.self%n) then
!     res%s = self%s
!     if (self%s.eq.'s') then
!       res%q(:) = self%q(:,i)
!     else
!       res%qd(:) = self%qd(:,i)
!     end if
!   else
!     call Message%printWarning('extractfromQuaternionArray: requested quaternion index larger than array size', &
!                               (/'   ---> returning empty quaternion'/) )
!     if (self%s.eq.'s') then
!       res = Quaternion_T( q = (/ 0.0, 0.0, 0.0, 0.0 /) )
!     else
!       res = Quaternion_T( qd = (/ 0.D0, 0.D0, 0.D0, 0.D0 /) )
!     end if
!   end if

! end function extractfromQuaternionArray


! !--------------------------------------------------------------------------!
! recursive subroutine insertQuatintoArray(self, i, q)
! !DEC$ ATTRIBUTES DLLEXPORT :: insertQuatintoArray
!   !! author: MDG
!   !! version: 1.0
!   !! date: 01/23/20
!   !!
!   !! insert a quaternion into an array of quaternions

! use mod_io

! IMPLICIT NONE

!   class(QuaternionArray_T),intent(inout):: self
!    !! input quaternion array
!   integer(kind=irg), intent(in)         :: i
!    !! quaternion to be extracted
!   type(Quaternion_T), intent(in)        :: q
!    !! extracted quaternion

!   type(IO_T)                            :: Message

!   if (i.le.self%n) then
!     if (self%s.eq.'s') then
!       self%q(:,i) = q%get_quats()
!     else
!       self%qd(:,i) = q%get_quatd()
!     end if
!   else
!     call Message%printWarning('extractfromQuaternionArray: requested quaternion index larger than array size', &
!                               (/'   ---> no quaternion inserted'/) )
!   end if

! end subroutine insertQuatintoArray

! !--------------------------------------------------------------------------!
! recursive function extractfromQuaternion3DArray(self, i) result (res)
! !DEC$ ATTRIBUTES DLLEXPORT :: extractfromQuaternion3DArray
!   !! author: MDG
!   !! version: 1.0
!   !! date: 06/03/21
!   !!
!   !! extract a quaternion from an array of quaternions

! use mod_io

! IMPLICIT NONE

!   class(Quaternion3DArray_T),intent(in) :: self
!    !! input quaternion array
!   integer(kind=irg), intent(in)         :: i(3)
!    !! quaternion to be extracted
!   type(Quaternion_T)                    :: res
!    !! extracted quaternion

!   type(IO_T)                            :: Message

!   if ( (i(1).le.self%n(1)) .and. (i(2).le.self%n(2)) .and. (i(3).le.self%n(3)) ) then
!     res%s = self%s
!     if (self%s.eq.'s') then
!       res%q(:) = self%q(:,i(1),i(2),i(3))
!     else
!       res%qd(:) = self%qd(:,i(1),i(2),i(3))
!     end if
!   else
!     call Message%printWarning('extractfromQuaternion3DArray: requested quaternion index larger than array size', &
!                               (/'   ---> returning empty quaternion'/) )
!     if (self%s.eq.'s') then
!       res = Quaternion_T( q = (/ 0.0, 0.0, 0.0, 0.0 /) )
!     else
!       res = Quaternion_T( qd = (/ 0.D0, 0.D0, 0.D0, 0.D0 /) )
!     end if
!   end if

! end function extractfromQuaternion3DArray


! !--------------------------------------------------------------------------!
! recursive subroutine insertQuatinto3DArray(self, i, q)
! !DEC$ ATTRIBUTES DLLEXPORT :: insertQuatinto3DArray
!   !! author: MDG
!   !! version: 1.0
!   !! date: 06/03/21
!   !!
!   !! insert a quaternion into an array of quaternions

! use mod_io

! IMPLICIT NONE

!   class(Quaternion3DArray_T),intent(inout):: self
!    !! input quaternion array
!   integer(kind=irg), intent(in)           :: i(3)
!    !! quaternion to be extracted
!   type(Quaternion_T), intent(in)          :: q
!    !! extracted quaternion

!   type(IO_T)                              :: Message

!   if ( (i(1).le.self%n(1)) .and. (i(2).le.self%n(2)) .and. (i(3).le.self%n(3)) ) then
!     if (self%s.eq.'s') then
!       self%q(:,i(1),i(2),i(3)) = q%get_quats()
!     else
!       self%qd(:,i(1),i(2),i(3)) = q%get_quatd()
!     end if
!   else
!     call Message%printWarning('extractfromQuaternion3DArray: requested quaternion index larger than array size', &
!                               (/'   ---> no quaternion inserted'/) )
!   end if

! end subroutine insertQuatinto3DArray

! !--------------------------------------------------------------------------!
! ! pure recursive function quatslerp(self, qb, n) result(res)
! recursive function quatslerp(self, qb, n) result(res)
! !DEC$ ATTRIBUTES DLLEXPORT :: quatslerp
!   !! author: MDG
!   !! version: 1.0
!   !! date: 10/11/22
!   !!
!   !! return an array of interpolated quaternions

! IMPLICIT NONE

!   class(Quaternion_T),intent(in)         :: self   ! = qa
!    !! input quaternion (start)
!   class(Quaternion_T),intent(in)         :: qb
!    !! input quaternion (end)
!   integer(kind=irg),intent(in)           :: n
!    !! number of steps in the interpolation
!   type(Quaternion_T)                     :: res(n)
!    !! output interpolated quaternion list

!   type(Quaternion_T)                     :: cqa
!   real(kind=sgl)                         :: theta, phi, dphi, s
!   real(kind=dbl)                         :: thetad, phid, dphid, sd
!   integer(kind=irg)                      :: i

!   if (self%s.eq.'s') then
!       do i=1,n
!         res(i)%q = (/ 0.0, 0.0, 0.0, 0.0 /)
!         res(i)%s = 's'
!       end do
!       cqa = quatconjg(self)
!       theta = acos( quatinnerproduct(qb, cqa ))
!       if (theta.ne.0.0) then
!         s = 1.0/sin(theta)
!         dphi = theta/real(n-1)

!         do i=1,n
!           phi = real(i-1)*dphi
!           res(i) = quatsmult(self, sin(theta-phi)*s ) + quatsmult(qb, sin(phi)*s )
!         end do
!       else
!         do i=1,n
!           res(i) = self
!         end do
!       end if
!   else
!       do i=1,n
!         res(i)%qd = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
!         res(i)%s = 'd'
!       end do
!       cqa = quatconjg(self)
!       thetad = acos( quatinnerproduct(qb, cqa ))
!       if (thetad.ne.0.D0) then
!         sd = 1.D0/sin(theta)
!         dphi = theta/dble(n-1)

!         do i=1,n
!           phi = dble(i-1)*dphi
!           res(i) = quatsmultd(self, sin(theta-phi)*sd ) + quatsmultd(qb, sin(phi)*sd )
!         end do
!       else
!         do i=1,n
!           res(i) = self
!         end do
!       end if
!   end if

! end function quatslerp

! !--------------------------------------------------------------------------
! recursive function quatsequal(self, qb) result(res)
! !DEC$ ATTRIBUTES DLLEXPORT :: quatsequal
!   !! author: MDG
!   !! version: 1.0
!   !! date: 10/11/22
!   !!
!   !! quaternion comparison (double precision)

! IMPLICIT NONE

!   class(Quaternion_T),intent(in)    :: self, qb
!    !! input quaternions
!   logical                           :: res

!   type(Quaternion_T)                :: diff

!   real(kind=sgl)                    :: d, eps=1.0e-7
!   real(kind=dbl)                    :: dd, epsd=1.0e-12

!   res = .TRUE.
!   diff = self - qb

!   if (self%s.eq.'s') then
!     d = maxval( abs( diff%q(:) ) )
!     if (d.gt.eps) res = .FALSE.
!   else
!     dd = maxval( abs( diff%qd(:) ) )
!     if (dd.gt.epsd) res = .FALSE.
!   end if

! end function quatsequal


! !--------------------------------------------------------------------------
! recursive function generateRandomArray(n, s, seed, northern) result(res)
! !DEC$ ATTRIBUTES DLLEXPORT :: generateRandomArray
!   !! author: MDG
!   !! version: 1.0
!   !! date: 10/13/22
!   !!
!   !! generate an array of random unit quaternions (single/double precision)

! use mod_rng

!   integer(kind=irg), intent(in)       :: n
!    !! number of unut quaternions to be generated
!   character(1), intent(in)            :: s
!    !! single 's' or double 'd' precision ?
!   type(rng_t), intent(inout)          :: seed
!    !! a seed number for the Marsaglia random number generator
!   logical, intent(in), OPTIONAL       :: northern

!   type(QuaternionArray_T)             :: res
!   type(Quaternion_T)                  :: q
!   integer(kind=irg)                   :: i
!   logical                             :: pos

!   pos = .FALSE.
!   if (present(northern)) then
!     if (northern.eqv..TRUE.) pos = .TRUE.
!   end if

!   res%s = s
!   if (s.eq.'s') then
!     allocate(res%q(4,n))
!   else
!     allocate(res%qd(4,n))
!   end if
!   res%n = n

!   do i=1,n
!     q = quat_Marsaglia(seed)
!     if (pos.eqv..TRUE.) then
!       if (q%qd(1).lt.0.0) q%qd = -q%qd
!     end if
!     if (s.eq.'s') then
!       res%q(:,i) = sngl(q%qd(:))
!     else
!       res%qd(:,i) = q%qd(:)
!     end if
!   end do

! end function generateRandomArray

! !--------------------------------------------------------------------------!
! recursive function quat_Marsaglia(seed) result(q)
! !DEC$ ATTRIBUTES DLLEXPORT :: quat_Marsaglia

! use mod_rng

! IMPLICIT NONE

! type(rng_t),INTENT(INOUT)           :: seed
! !f2py intent(in,out) ::  seed
! type(Quaternion_T)                  :: q

! real(kind=dbl)                      :: x1,x2,y1,y2,s1,s2

! x1 = 2.D0*rng_uniform(seed)-1.D0
! y1 = 2.D0*rng_uniform(seed)-1.D0
! s1 = x1*x1+y1*y1
! if (s1.gt.1.D0) then
!   do while (s1.gt.1.D0)
!     x1 = 2.D0*rng_uniform(seed)-1.D0
!     x2 = 2.D0*rng_uniform(seed)-1.D0
!     s1 = x1*x1+y1*y1
!   end do
! end if

! x2 = 2.D0*rng_uniform(seed)-1.D0
! y2 = 2.D0*rng_uniform(seed)-1.D0
! s2 = x2*x2+y2*y2
! if (s2.gt.1.D0) then
!   do while (s2.gt.1.D0)
!     x2 = 2.D0*rng_uniform(seed)-1.D0
!     y2 = 2.D0*rng_uniform(seed)-1.D0
!     s2 = x2*x2+y2*y2
!   end do
! end if

! s1 = sqrt( (1.D0-s2)/s2 )

! q%qd = (/ x1, y1, x2*s1, y2*s1 /)
! q%s = 'd'

! end function quat_Marsaglia

! !--------------------------------------------------------------------------
! recursive subroutine QSym_Init_(self, pgnum, qsym)
! !DEC$ ATTRIBUTES DLLEXPORT :: QSym_Init_
!   !! author: MDG
!   !! version: 1.0
!   !! date: 10/13/22
!   !!
!   !! initialize the quaternion symmetry operators for a given rotational point group

! use mod_symmetry
! use mod_io

! IMPLICIT NONE

! class(QuaternionArray_T), INTENT(INOUT)   :: self
! integer(kind=irg), INTENT(IN)             :: pgnum
! type(QuaternionArray_T), INTENT(OUT)      :: qsym

! type(IO_T)                                :: Message
! integer(kind=irg)                         :: i, Nqsym, prot
! real(kind=dbl), allocatable               :: Pm(:,:)

! ! here we need to analyze the rotational symmetry group, and copy the appropriate
! ! quaternion symmetry operators into a QuaternionArray_T class.

! ! first get the number of the rotational point group that corresponds to the crystal point group
! prot = PGrot(pgnum)
! ! possible values for prot are: (/1,3,6,9,12,16,18,21,24,28,30/)
! ! corresponding to the point groups 1, 2, 222, 4, 422, 3, 32, 6, 622, 23, 432 and 532 respectively

! !------------
! ! IMPORTANT NOTE: the original von Mises-Fischer (VMF) approach requires that q and -q are considered to
! ! be separate quaternions, so the original Matlab code included the negatives of all quaternion symmetry operators
! ! as well, leading to a cardinality of twice the rotational point group order.  It appears that we do not have to
! ! do so if we replace the exponential in the VMF by a hyperbolic cosine function, which would account directly
! ! for the q, -q duplicity... Alternatively, one can use the axial Watson distribution.
! !------------

! ! identity operator is part of all point groups

! ! select statement for each individual rotational point group (see typedefs.f90 for SYM_Qsymop definitions)
! select case (prot)
!         case(1)         ! 1 (no additional symmetry elements)
!                 allocate(Pm(4,1))
!                 Pm(1:4,1) = SYM_Qsymop(1:4,1)
!                 Nqsym = 1

!         case(3)         ! 2  (we'll assume that the two-fold axis lies along the e_y-axis)
!                 allocate(Pm(4,2))
!                 Pm(1:4,1) = SYM_Qsymop(1:4,1)
!                 Nqsym = 2
!                 Pm(1:4,2) = SYM_Qsymop(1:4,3)

!         case(6)         ! 222
!                 allocate(Pm(4,4))
!                 Pm(1:4,1) = SYM_Qsymop(1:4,1)
!                 Nqsym = 4
!                 do i=2,4
!                   Pm(1:4,i) = SYM_Qsymop(1:4,i)
!                 end do

!         case(9)         ! 4
!                 allocate(Pm(4,4))
!                 Pm(1:4,1) = SYM_Qsymop(1:4,1)
!                 Nqsym = 4
!                 Pm(1:4,2) = SYM_Qsymop(1:4,4)
!                 Pm(1:4,3) = SYM_Qsymop(1:4,7)
!                 Pm(1:4,4) = SYM_Qsymop(1:4,10)

!         case(12)        ! 422
!                 allocate(Pm(4,8))
!                 Pm(1:4,1) = SYM_Qsymop(1:4,1)
!                 Nqsym = 8
!                 Pm(1:4,2) = SYM_Qsymop(1:4,4)
!                 Pm(1:4,3) = SYM_Qsymop(1:4,7)
!                 Pm(1:4,4) = SYM_Qsymop(1:4,10)
!                 Pm(1:4,5) = SYM_Qsymop(1:4,2)
!                 Pm(1:4,6) = SYM_Qsymop(1:4,3)
!                 Pm(1:4,7) = SYM_Qsymop(1:4,11)
!                 Pm(1:4,8) = SYM_Qsymop(1:4,12)

!         case(16)        ! 3
!                 allocate(Pm(4,3))
!                 Pm(1:4,1) = SYM_Qsymop(1:4,1)
!                 Nqsym = 3
!                 Pm(1:4,2) = SYM_Qsymop(1:4,26)
!                 Pm(1:4,3) = SYM_Qsymop(1:4,28)

!         case(18)        ! 32 (needs special handling)
!                 allocate(Pm(4,6))
!                 Pm(1:4,1) = SYM_Qsymop(1:4,1)
!                 Nqsym = 6
!                 Pm(1:4,2) = SYM_Qsymop(1:4,26)
!                 Pm(1:4,3) = SYM_Qsymop(1:4,28)
!                 Pm(1:4,4) = SYM_Qsymop(1:4,30)
!                 Pm(1:4,5) = SYM_Qsymop(1:4,32)
!                 Pm(1:4,6) = SYM_Qsymop(1:4,34)

!         case(21)        ! 6
!                 allocate(Pm(4,6))
!                 Pm(1:4,1) = SYM_Qsymop(1:4,1)
!                 Nqsym = 6
!                 do i=25,29
!                   Pm(1:4,i-23) = SYM_Qsymop(1:4,i)
!                 end do

!         case(24)        ! 622
!                 allocate(Pm(4,6))
!                 Pm(1:4,1) = SYM_Qsymop(1:4,1)
!                 Nqsym = 12
!                 do i=25,35
!                   Pm(1:4,i-23) = SYM_Qsymop(1:4,i)
!                 end do

!         case(28)        ! 23
!                 allocate(Pm(4,12))
!                 Pm(1:4,1) = SYM_Qsymop(1:4,1)
!                 Nqsym = 12
!                 do i=2,4
!                   Pm(1:4,i) = SYM_Qsymop(1:4,i)
!                 end do
!                 do i=17,24
!                   Pm(1:4,4+(i-16)) = SYM_Qsymop(1:4,i)
!                 end do

!         case(30)        ! 432
!                 allocate(Pm(4,24))
!                 Pm(1:4,1) = SYM_Qsymop(1:4,1)
!                 Nqsym = 24
!                 do i=2,24
!                   Pm(1:4,i) = SYM_Qsymop(1:4,i)
!                 end do
!         case(33)        ! 532
!                 allocate(Pm(4,60))
!                 Pm(1:4,1) = SYM_Qsymop(1:4,1)
!                 Nqsym = 60
!                 do i=2,60
!                   Pm(1:4,i) = SYM_Qsymop(1:4,35+i)
!                 end do
!         case(34)  ! 822
!                 allocate(Pm(4,16))
!                 Pm(1:4,1) = SYM_Qsymop(1:4,1)
!                 Nqsym = 16
!                 do i = 1,15
!                   Pm(1:4,i+1) = SYM_Qsymop(1:4,95+i)
!                 end do
!         case(35)  ! 1022
!                 allocate(Pm(4,20))
!                 Pm(1:4,1) = SYM_Qsymop(1:4,1)
!                 Nqsym = 20
!                 do i = 1,19
!                   Pm(1:4,i+1) = SYM_Qsymop(1:4,110+i)
!                 end do
!         case(36)  ! 1222
!                 allocate(Pm(4,24))
!                 Pm(1:4,1) = SYM_Qsymop(1:4,1)
!                 Nqsym = 24
!                 do i = 1,23
!                   Pm(1:4,i+1) = SYM_Qsymop(1:4,129+i)
!                 end do

!         case default    ! this should never happen ...
!                 write (*,*) 'requested rotational point group ', prot
!                 call Message%printError('QSym_Init','unknown rotational point group number')
! end select

! ! and initialize the quatenrion array class
! qsym = QuaternionArray_T( n = Nqsym, nthreads = 1, qd = Pm )

! end subroutine QSym_Init_

! !--------------------------------------------------------------------------
! recursive function getQnumber_(self) result(num)
! !DEC$ ATTRIBUTES DLLEXPORT :: getQnumber_
!   !! author: MDG
!   !! version: 1.0
!   !! date: 10/13/22
!   !!
!   !! returns the number of quaternions in the QuaternionArray_T class

! IMPLICIT NONE

! class(QuaternionArray_T), INTENT(INOUT)   :: self
! integer(kind=irg)                         :: num

! num = self%n

! end function getQnumber_

! !--------------------------------------------------------------------------
! recursive function get3DQnumber_(self) result(num)
! !DEC$ ATTRIBUTES DLLEXPORT :: get3DQnumber_
!   !! author: MDG
!   !! version: 1.0
!   !! date: 07/18/21
!   !!
!   !! returns the number of quaternions in the Quaternion3DArray_T class

! IMPLICIT NONE

! class(Quaternion3DArray_T), INTENT(INOUT)    :: self
! integer(kind=irg)                            :: num(3)

! num = self%n

! end function get3DQnumber_

! !--------------------------------------------------------------------------
! recursive subroutine deleteArray_(self)
! !DEC$ ATTRIBUTES DLLEXPORT :: deleteArray_
!   !! author: MDG
!   !! version: 1.0
!   !! date: 07/16/21
!   !!
!   !! deletes the current array of quaternions in this class

! IMPLICIT NONE

! class(QuaternionArray_T), INTENT(INOUT)   :: self

! if (self%s.eq.'s') then 
!   if (allocated(self%q)) deallocate(self%q)
! else 
!   if (allocated(self%qd)) deallocate(self%qd)
! end if

! self%n = 0

! end subroutine deleteArray_

end module mod_dualquaternions
