! ###################################################################
! Copyright (c) 2015-2026, Marc De Graef Research Group/Carnegie Mellon University
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

module c_quaternions
  !! author: MDG
  !! version: 1.0
  !! date: 03/27/26
  !!
  !! C-interop wrappers for Quaternion_T and QuaternionArray_T.
  !! All functions use bind(c) and opaque handles for Fortran objects.

use iso_c_binding
use mod_kinds
use mod_quaternions

IMPLICIT NONE

private

contains

!--------------------------------------------------------------------------
! Quaternion_T constructors and destructors
!--------------------------------------------------------------------------

function c_quat_create(qd) result(handle) bind(c, name='emsoft_quat_create')
  !! Create a Quaternion_T from a double-precision 4-vector [w, x, y, z].
  real(c_double), INTENT(IN) :: qd(4)
  type(c_ptr)                :: handle
  type(Quaternion_T), pointer :: obj

  allocate(obj)
  obj = Quaternion_T( qd = qd )
  handle = c_loc(obj)

end function c_quat_create

!--------------------------------------------------------------------------
function c_quat_create_identity() result(handle) bind(c, name='emsoft_quat_create_identity')
  !! Create the identity quaternion [1, 0, 0, 0].
  type(c_ptr)                :: handle
  type(Quaternion_T), pointer :: obj

  allocate(obj)
  obj = Quaternion_T( qd = (/ 1.D0, 0.D0, 0.D0, 0.D0 /) )
  handle = c_loc(obj)

end function c_quat_create_identity

!--------------------------------------------------------------------------
subroutine c_quat_destroy(handle) bind(c, name='emsoft_quat_destroy')
  !! Destroy a Quaternion_T and free its memory.
  type(c_ptr), value, INTENT(IN) :: handle
  type(Quaternion_T), pointer    :: obj

  call c_f_pointer(handle, obj)
  deallocate(obj)

end subroutine c_quat_destroy

!--------------------------------------------------------------------------
! Getters and setters
!--------------------------------------------------------------------------

subroutine c_quat_get(handle, qd) bind(c, name='emsoft_quat_get')
  !! Get the quaternion components as double precision.
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double), INTENT(OUT)    :: qd(4)
  type(Quaternion_T), pointer    :: obj

  call c_f_pointer(handle, obj)
  qd = obj%get_quatd()

end subroutine c_quat_get

!--------------------------------------------------------------------------
subroutine c_quat_set(handle, qd) bind(c, name='emsoft_quat_set')
  !! Set the quaternion components from double precision.
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double), INTENT(IN)     :: qd(4)
  type(Quaternion_T), pointer    :: obj

  call c_f_pointer(handle, obj)
  call obj%set_quatd(qd)

end subroutine c_quat_set

!--------------------------------------------------------------------------
! Arithmetic operations
!--------------------------------------------------------------------------

function c_quat_add(h1, h2) result(handle) bind(c, name='emsoft_quat_add')
  !! Add two quaternions: result = h1 + h2.
  type(c_ptr), value, INTENT(IN) :: h1, h2
  type(c_ptr)                    :: handle
  type(Quaternion_T), pointer    :: q1, q2, qres

  call c_f_pointer(h1, q1)
  call c_f_pointer(h2, q2)
  allocate(qres)
  qres = q1 + q2
  handle = c_loc(qres)

end function c_quat_add

!--------------------------------------------------------------------------
function c_quat_subtract(h1, h2) result(handle) bind(c, name='emsoft_quat_subtract')
  !! Subtract two quaternions: result = h1 - h2.
  type(c_ptr), value, INTENT(IN) :: h1, h2
  type(c_ptr)                    :: handle
  type(Quaternion_T), pointer    :: q1, q2, qres

  call c_f_pointer(h1, q1)
  call c_f_pointer(h2, q2)
  allocate(qres)
  qres = q1 - q2
  handle = c_loc(qres)

end function c_quat_subtract

!--------------------------------------------------------------------------
function c_quat_multiply(h1, h2) result(handle) bind(c, name='emsoft_quat_multiply')
  !! Multiply two quaternions: result = h1 * h2.
  type(c_ptr), value, INTENT(IN) :: h1, h2
  type(c_ptr)                    :: handle
  type(Quaternion_T), pointer    :: q1, q2, qres

  call c_f_pointer(h1, q1)
  call c_f_pointer(h2, q2)
  allocate(qres)
  qres = q1 * q2
  handle = c_loc(qres)

end function c_quat_multiply

!--------------------------------------------------------------------------
function c_quat_divide(h1, h2) result(handle) bind(c, name='emsoft_quat_divide')
  !! Divide two quaternions: result = h1 / h2.
  type(c_ptr), value, INTENT(IN) :: h1, h2
  type(c_ptr)                    :: handle
  type(Quaternion_T), pointer    :: q1, q2, qres

  call c_f_pointer(h1, q1)
  call c_f_pointer(h2, q2)
  allocate(qres)
  qres = q1 / q2
  handle = c_loc(qres)

end function c_quat_divide

!--------------------------------------------------------------------------
function c_quat_scale(h, s) result(handle) bind(c, name='emsoft_quat_scale')
  !! Scalar multiplication: result = h * s.
  type(c_ptr), value, INTENT(IN) :: h
  real(c_double), value, INTENT(IN) :: s
  type(c_ptr)                    :: handle
  type(Quaternion_T), pointer    :: q, qres

  call c_f_pointer(h, q)
  allocate(qres)
  qres = q * s
  handle = c_loc(qres)

end function c_quat_scale

!--------------------------------------------------------------------------
! Conjugate, norm, and normalization
!--------------------------------------------------------------------------

function c_quat_conjugate(h) result(handle) bind(c, name='emsoft_quat_conjugate')
  !! Return the conjugate of a quaternion.
  type(c_ptr), value, INTENT(IN) :: h
  type(c_ptr)                    :: handle
  type(Quaternion_T), pointer    :: q, qres

  call c_f_pointer(h, q)
  allocate(qres)
  qres = conjg(q)
  handle = c_loc(qres)

end function c_quat_conjugate

!--------------------------------------------------------------------------
function c_quat_norm(h) result(res) bind(c, name='emsoft_quat_norm')
  !! Return the norm of a quaternion.
  type(c_ptr), value, INTENT(IN) :: h
  real(c_double)                 :: res
  type(Quaternion_T), pointer    :: q

  call c_f_pointer(h, q)
  res = cabs(q)

end function c_quat_norm

!--------------------------------------------------------------------------
subroutine c_quat_normalize(h) bind(c, name='emsoft_quat_normalize')
  !! Normalize a quaternion in-place.
  type(c_ptr), value, INTENT(IN) :: h
  type(Quaternion_T), pointer    :: q

  call c_f_pointer(h, q)
  call q%quat_normalize()

end subroutine c_quat_normalize

!--------------------------------------------------------------------------
subroutine c_quat_flip(h) bind(c, name='emsoft_quat_flip')
  !! Negate all components of a quaternion in-place.
  type(c_ptr), value, INTENT(IN) :: h
  type(Quaternion_T), pointer    :: q

  call c_f_pointer(h, q)
  call q%quat_flip()

end subroutine c_quat_flip

!--------------------------------------------------------------------------
subroutine c_quat_pos(h) bind(c, name='emsoft_quat_pos')
  !! Make scalar part positive in-place (q or -q convention).
  type(c_ptr), value, INTENT(IN) :: h
  type(Quaternion_T), pointer    :: q

  call c_f_pointer(h, q)
  call q%quat_pos()

end subroutine c_quat_pos

!--------------------------------------------------------------------------
! Geometric operations
!--------------------------------------------------------------------------

function c_quat_innerproduct(h1, h2) result(res) bind(c, name='emsoft_quat_innerproduct')
  !! Compute the inner product of two quaternions.
  type(c_ptr), value, INTENT(IN) :: h1, h2
  real(c_double)                 :: res
  type(Quaternion_T), pointer    :: q1, q2

  call c_f_pointer(h1, q1)
  call c_f_pointer(h2, q2)
  res = q1%quat_innerproduct(q2)

end function c_quat_innerproduct

!--------------------------------------------------------------------------
function c_quat_angle(h1, h2) result(res) bind(c, name='emsoft_quat_angle')
  !! Compute the angle (in radians) between two unit quaternions.
  type(c_ptr), value, INTENT(IN) :: h1, h2
  real(c_double)                 :: res
  type(Quaternion_T), pointer    :: q1, q2

  call c_f_pointer(h1, q1)
  call c_f_pointer(h2, q2)
  res = q1%quat_angle(q2)

end function c_quat_angle

!--------------------------------------------------------------------------
function c_quat_equal(h1, h2) result(res) bind(c, name='emsoft_quat_equal')
  !! Test equality of two quaternions.
  type(c_ptr), value, INTENT(IN) :: h1, h2
  logical(c_bool)                :: res
  type(Quaternion_T), pointer    :: q1, q2

  call c_f_pointer(h1, q1)
  call c_f_pointer(h2, q2)
  res = logical(q1%quatsequal(q2), c_bool)

end function c_quat_equal

!--------------------------------------------------------------------------
! Vector rotation
!--------------------------------------------------------------------------

subroutine c_quat_rotate_vector(h, v, vout) bind(c, name='emsoft_quat_rotate_vector')
  !! Rotate a 3-vector by a unit quaternion using the L_p operation: v' = q v q*.
  type(c_ptr), value, INTENT(IN) :: h
  real(c_double), INTENT(IN)     :: v(3)
  real(c_double), INTENT(OUT)    :: vout(3)
  type(Quaternion_T), pointer    :: q

  call c_f_pointer(h, q)
  vout = q%quat_Lp(v)

end subroutine c_quat_rotate_vector

!--------------------------------------------------------------------------
subroutine c_quat_rotate_vecarray(h, n, v, vout) bind(c, name='emsoft_quat_rotate_vecarray')
  !! Rotate an array of n 3-vectors by a unit quaternion.
  type(c_ptr), value, INTENT(IN)    :: h
  integer(c_int), value, INTENT(IN) :: n
  real(c_double), INTENT(IN)        :: v(3,n)
  real(c_double), INTENT(OUT)       :: vout(3,n)
  type(Quaternion_T), pointer       :: q

  call c_f_pointer(h, q)
  vout = q%quat_Lp_vecarray(n, v)

end subroutine c_quat_rotate_vecarray

!--------------------------------------------------------------------------
! QuaternionArray_T constructors and destructors
!--------------------------------------------------------------------------

function c_quatarray_create(n, qd) result(handle) bind(c, name='emsoft_quatarray_create')
  !! Create a QuaternionArray_T from n double-precision quaternions stored as (4,n).
  integer(c_int), value, INTENT(IN) :: n
  real(c_double), INTENT(IN)        :: qd(4,n)
  type(c_ptr)                        :: handle
  type(QuaternionArray_T), pointer   :: obj

  allocate(obj)
  obj = QuaternionArray_T( n = n, qd = qd, s = 'd' )
  handle = c_loc(obj)

end function c_quatarray_create

!--------------------------------------------------------------------------
function c_quatarray_create_empty(n) result(handle) bind(c, name='emsoft_quatarray_create_empty')
  !! Create an empty QuaternionArray_T with n slots.
  integer(c_int), value, INTENT(IN) :: n
  type(c_ptr)                        :: handle
  type(QuaternionArray_T), pointer   :: obj

  allocate(obj)
  obj = QuaternionArray_T( n = n, s = 'd' )
  handle = c_loc(obj)

end function c_quatarray_create_empty

!--------------------------------------------------------------------------
subroutine c_quatarray_destroy(handle) bind(c, name='emsoft_quatarray_destroy')
  !! Destroy a QuaternionArray_T and free its memory.
  type(c_ptr), value, INTENT(IN)   :: handle
  type(QuaternionArray_T), pointer :: obj

  call c_f_pointer(handle, obj)
  call obj%deleteArray()
  deallocate(obj)

end subroutine c_quatarray_destroy

!--------------------------------------------------------------------------
! QuaternionArray_T getters
!--------------------------------------------------------------------------

function c_quatarray_size(handle) result(n) bind(c, name='emsoft_quatarray_size')
  !! Return the number of quaternions in the array.
  type(c_ptr), value, INTENT(IN)   :: handle
  integer(c_int)                   :: n
  type(QuaternionArray_T), pointer :: obj

  call c_f_pointer(handle, obj)
  n = obj%getQnumber()

end function c_quatarray_size

!--------------------------------------------------------------------------
subroutine c_quatarray_get_element(handle, i, qd) bind(c, name='emsoft_quatarray_get_element')
  !! Get the i-th quaternion (1-based index) as a double-precision 4-vector.
  type(c_ptr), value, INTENT(IN)    :: handle
  integer(c_int), value, INTENT(IN) :: i
  real(c_double), INTENT(OUT)       :: qd(4)
  type(QuaternionArray_T), pointer  :: obj
  type(Quaternion_T)                :: qt

  call c_f_pointer(handle, obj)
  qt = obj%getQuatfromArray(i)
  qd = qt%get_quatd()

end subroutine c_quatarray_get_element

!--------------------------------------------------------------------------
subroutine c_quatarray_set_element(handle, i, qd) bind(c, name='emsoft_quatarray_set_element')
  !! Set the i-th quaternion (1-based index) from a double-precision 4-vector.
  type(c_ptr), value, INTENT(IN)    :: handle
  integer(c_int), value, INTENT(IN) :: i
  real(c_double), INTENT(IN)        :: qd(4)
  type(QuaternionArray_T), pointer  :: obj
  type(Quaternion_T)                :: qt

  call c_f_pointer(handle, obj)
  qt = Quaternion_T( qd = qd )
  call obj%insertQuatinArray(i, qt)

end subroutine c_quatarray_set_element

!--------------------------------------------------------------------------
! QuaternionArray_T arithmetic
!--------------------------------------------------------------------------

function c_quatarray_multiply(h1, h2) result(handle) bind(c, name='emsoft_quatarray_multiply')
  !! Element-wise multiplication of two quaternion arrays.
  type(c_ptr), value, INTENT(IN)   :: h1, h2
  type(c_ptr)                      :: handle
  type(QuaternionArray_T), pointer :: a1, a2, ares

  call c_f_pointer(h1, a1)
  call c_f_pointer(h2, a2)
  allocate(ares)
  ares = a1 * a2
  handle = c_loc(ares)

end function c_quatarray_multiply

!--------------------------------------------------------------------------
subroutine c_quatarray_normalize(handle) bind(c, name='emsoft_quatarray_normalize')
  !! Normalize all quaternions in the array in-place.
  type(c_ptr), value, INTENT(IN)   :: handle
  type(QuaternionArray_T), pointer :: obj

  call c_f_pointer(handle, obj)
  call obj%quat_normalize()

end subroutine c_quatarray_normalize

!--------------------------------------------------------------------------
subroutine c_quatarray_rotate_vector(handle, v, vout, n) bind(c, name='emsoft_quatarray_rotate_vector')
  !! Rotate a single vector by each quaternion in the array.
  !! vout is (3,n) where n = number of quaternions.
  type(c_ptr), value, INTENT(IN)    :: handle
  real(c_double), INTENT(IN)        :: v(3)
  integer(c_int), value, INTENT(IN) :: n
  real(c_double), INTENT(OUT)       :: vout(3,n)
  type(QuaternionArray_T), pointer  :: obj

  call c_f_pointer(handle, obj)
  vout = obj%quat_Lp(v)

end subroutine c_quatarray_rotate_vector

end module c_quaternions
