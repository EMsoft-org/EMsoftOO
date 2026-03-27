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

module c_rotations
  !! author: MDG
  !! version: 1.0
  !! date: 03/27/26
  !!
  !! C-interop wrappers for rotation representations using Orientation_T.
  !! Provides a unified interface: create from any representation,
  !! extract any representation.

use iso_c_binding
use mod_kinds
use mod_rotations

IMPLICIT NONE

private

contains

!--------------------------------------------------------------------------
! Orientation_T constructors — one per representation
!--------------------------------------------------------------------------

function c_rot_from_euler(eu) result(handle) bind(c, name='emsoft_rot_from_euler')
  !! Create from Euler angles [phi1, Phi, phi2] in radians.
  real(c_double), INTENT(IN) :: eu(3)
  type(c_ptr)                :: handle
  type(orientation_T), pointer :: obj
  type(e_T)                  :: e

  call setRotationPrecision('d')
  allocate(obj)
  e = e_T( edinp = eu )
  obj = orientation_T( e )
  handle = c_loc(obj)

end function c_rot_from_euler

!--------------------------------------------------------------------------
function c_rot_from_quaternion(qu) result(handle) bind(c, name='emsoft_rot_from_quaternion')
  !! Create from unit quaternion [w, x, y, z].
  real(c_double), INTENT(IN) :: qu(4)
  type(c_ptr)                :: handle
  type(orientation_T), pointer :: obj
  type(q_T)                  :: q

  call setRotationPrecision('d')
  allocate(obj)
  q = q_T( qdinp = qu )
  obj = orientation_T( q )
  handle = c_loc(obj)

end function c_rot_from_quaternion

!--------------------------------------------------------------------------
function c_rot_from_matrix(om) result(handle) bind(c, name='emsoft_rot_from_matrix')
  !! Create from 3x3 rotation matrix.
  real(c_double), INTENT(IN) :: om(3,3)
  type(c_ptr)                :: handle
  type(orientation_T), pointer :: obj
  type(o_T)                  :: o

  call setRotationPrecision('d')
  allocate(obj)
  o = o_T( odinp = om )
  obj = orientation_T( o )
  handle = c_loc(obj)

end function c_rot_from_matrix

!--------------------------------------------------------------------------
function c_rot_from_axisangle(ax) result(handle) bind(c, name='emsoft_rot_from_axisangle')
  !! Create from axis-angle pair [n1, n2, n3, angle] (angle in radians).
  real(c_double), INTENT(IN) :: ax(4)
  type(c_ptr)                :: handle
  type(orientation_T), pointer :: obj
  type(a_T)                  :: a

  call setRotationPrecision('d')
  allocate(obj)
  a = a_T( adinp = ax )
  obj = orientation_T( a )
  handle = c_loc(obj)

end function c_rot_from_axisangle

!--------------------------------------------------------------------------
function c_rot_from_rodrigues(ro) result(handle) bind(c, name='emsoft_rot_from_rodrigues')
  !! Create from Rodrigues vector [n1, n2, n3, tan(angle/2)].
  real(c_double), INTENT(IN) :: ro(4)
  type(c_ptr)                :: handle
  type(orientation_T), pointer :: obj
  type(r_T)                  :: r

  call setRotationPrecision('d')
  allocate(obj)
  r = r_T( rdinp = ro )
  obj = orientation_T( r )
  handle = c_loc(obj)

end function c_rot_from_rodrigues

!--------------------------------------------------------------------------
function c_rot_from_homochoric(ho) result(handle) bind(c, name='emsoft_rot_from_homochoric')
  !! Create from homochoric vector [h1, h2, h3].
  real(c_double), INTENT(IN) :: ho(3)
  type(c_ptr)                :: handle
  type(orientation_T), pointer :: obj
  type(h_T)                  :: h

  call setRotationPrecision('d')
  allocate(obj)
  h = h_T( hdinp = ho )
  obj = orientation_T( h )
  handle = c_loc(obj)

end function c_rot_from_homochoric

!--------------------------------------------------------------------------
function c_rot_from_cubochoric(cu) result(handle) bind(c, name='emsoft_rot_from_cubochoric')
  !! Create from cubochoric vector [c1, c2, c3].
  real(c_double), INTENT(IN) :: cu(3)
  type(c_ptr)                :: handle
  type(orientation_T), pointer :: obj
  type(c_T)                  :: cc

  call setRotationPrecision('d')
  allocate(obj)
  cc = c_T( cdinp = cu )
  obj = orientation_T( cc )
  handle = c_loc(obj)

end function c_rot_from_cubochoric

!--------------------------------------------------------------------------
function c_rot_from_stereographic(st) result(handle) bind(c, name='emsoft_rot_from_stereographic')
  !! Create from stereographic vector [s1, s2, s3].
  real(c_double), INTENT(IN) :: st(3)
  type(c_ptr)                :: handle
  type(orientation_T), pointer :: obj
  type(s_T)                  :: s

  call setRotationPrecision('d')
  allocate(obj)
  s = s_T( sdinp = st )
  obj = orientation_T( s )
  handle = c_loc(obj)

end function c_rot_from_stereographic

!--------------------------------------------------------------------------
function c_rot_from_rotvec(rv) result(handle) bind(c, name='emsoft_rot_from_rotvec')
  !! Create from rotation vector [v1, v2, v3].
  real(c_double), INTENT(IN) :: rv(3)
  type(c_ptr)                :: handle
  type(orientation_T), pointer :: obj
  type(v_T)                  :: v

  call setRotationPrecision('d')
  allocate(obj)
  v = v_T( vdinp = rv )
  obj = orientation_T( v )
  handle = c_loc(obj)

end function c_rot_from_rotvec

!--------------------------------------------------------------------------
! Destructor
!--------------------------------------------------------------------------

subroutine c_rot_destroy(handle) bind(c, name='emsoft_rot_destroy')
  type(c_ptr), value, INTENT(IN)   :: handle
  type(orientation_T), pointer     :: obj

  call c_f_pointer(handle, obj)
  deallocate(obj)

end subroutine c_rot_destroy

!--------------------------------------------------------------------------
! Representation extractors — all double precision
!--------------------------------------------------------------------------

subroutine c_rot_to_euler(handle, eu) bind(c, name='emsoft_rot_to_euler')
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double), INTENT(OUT)    :: eu(3)
  type(orientation_T), pointer   :: obj

  call c_f_pointer(handle, obj)
  eu = obj%get_ed()

end subroutine c_rot_to_euler

!--------------------------------------------------------------------------
subroutine c_rot_to_quaternion(handle, qu) bind(c, name='emsoft_rot_to_quaternion')
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double), INTENT(OUT)    :: qu(4)
  type(orientation_T), pointer   :: obj

  call c_f_pointer(handle, obj)
  qu = obj%get_qd()

end subroutine c_rot_to_quaternion

!--------------------------------------------------------------------------
subroutine c_rot_to_matrix(handle, om) bind(c, name='emsoft_rot_to_matrix')
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double), INTENT(OUT)    :: om(3,3)
  type(orientation_T), pointer   :: obj

  call c_f_pointer(handle, obj)
  om = obj%get_od()

end subroutine c_rot_to_matrix

!--------------------------------------------------------------------------
subroutine c_rot_to_axisangle(handle, ax) bind(c, name='emsoft_rot_to_axisangle')
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double), INTENT(OUT)    :: ax(4)
  type(orientation_T), pointer   :: obj

  call c_f_pointer(handle, obj)
  ax = obj%get_ad()

end subroutine c_rot_to_axisangle

!--------------------------------------------------------------------------
subroutine c_rot_to_rodrigues(handle, ro) bind(c, name='emsoft_rot_to_rodrigues')
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double), INTENT(OUT)    :: ro(4)
  type(orientation_T), pointer   :: obj

  call c_f_pointer(handle, obj)
  ro = obj%get_rd()

end subroutine c_rot_to_rodrigues

!--------------------------------------------------------------------------
subroutine c_rot_to_homochoric(handle, ho) bind(c, name='emsoft_rot_to_homochoric')
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double), INTENT(OUT)    :: ho(3)
  type(orientation_T), pointer   :: obj

  call c_f_pointer(handle, obj)
  ho = obj%get_hd()

end subroutine c_rot_to_homochoric

!--------------------------------------------------------------------------
subroutine c_rot_to_cubochoric(handle, cu) bind(c, name='emsoft_rot_to_cubochoric')
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double), INTENT(OUT)    :: cu(3)
  type(orientation_T), pointer   :: obj

  call c_f_pointer(handle, obj)
  cu = obj%get_cd()

end subroutine c_rot_to_cubochoric

!--------------------------------------------------------------------------
subroutine c_rot_to_stereographic(handle, st) bind(c, name='emsoft_rot_to_stereographic')
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double), INTENT(OUT)    :: st(3)
  type(orientation_T), pointer   :: obj

  call c_f_pointer(handle, obj)
  st = obj%get_sd()

end subroutine c_rot_to_stereographic

!--------------------------------------------------------------------------
subroutine c_rot_to_rotvec(handle, rv) bind(c, name='emsoft_rot_to_rotvec')
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double), INTENT(OUT)    :: rv(3)
  type(orientation_T), pointer   :: obj

  call c_f_pointer(handle, obj)
  rv = obj%get_vd()

end subroutine c_rot_to_rotvec

end module c_rotations
