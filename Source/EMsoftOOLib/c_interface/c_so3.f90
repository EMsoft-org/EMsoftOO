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

module c_so3
  !! author: MDG
  !! version: 1.0
  !! date: 03/27/26
  !!
  !! C-interop wrappers for so3_T (SO(3) sampling and fundamental zones).

use iso_c_binding
use mod_kinds
use mod_global
use mod_so3
use mod_rotations

IMPLICIT NONE

private

contains

!--------------------------------------------------------------------------
! Constructor / destructor
!--------------------------------------------------------------------------

function c_so3_create(pgnum) result(handle) bind(c, name='emsoft_so3_create')
  !! Create a so3_T for the given point group number (1-32).
  integer(c_int), value, INTENT(IN) :: pgnum
  type(c_ptr)                        :: handle
  type(so3_T), pointer               :: obj

  allocate(obj)
  obj = so3_T(pgnum)
  handle = c_loc(obj)

end function c_so3_create

!--------------------------------------------------------------------------
subroutine c_so3_destroy(handle) bind(c, name='emsoft_so3_destroy')
  type(c_ptr), value, INTENT(IN) :: handle
  type(so3_T), pointer           :: obj

  call c_f_pointer(handle, obj)
  deallocate(obj)

end subroutine c_so3_destroy

!--------------------------------------------------------------------------
! Fundamental zone properties
!--------------------------------------------------------------------------

subroutine c_so3_get_fz_type_order(handle, fztype, fzorder) &
    bind(c, name='emsoft_so3_get_fz_type_order')
  !! Get the fundamental zone type and order.
  !! Types: 0=none, 1=cyclic, 2=dihedral, 3=tetrahedral, 4=octahedral
  type(c_ptr), value, INTENT(IN) :: handle
  integer(c_int), INTENT(OUT)   :: fztype
  integer(c_int), INTENT(OUT)   :: fzorder
  type(so3_T), pointer           :: obj

  call c_f_pointer(handle, obj)
  call obj%getFZtypeandorder(fztype, fzorder)

end subroutine c_so3_get_fz_type_order

!--------------------------------------------------------------------------
! Fundamental zone membership test
!--------------------------------------------------------------------------

function c_so3_is_inside_fz(handle, rod) result(inside) &
    bind(c, name='emsoft_so3_is_inside_fz')
  !! Test if a Rodrigues vector [n1, n2, n3, tan(angle/2)] is inside the
  !! fundamental zone for this point group.
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double), INTENT(IN)    :: rod(4)
  logical(c_bool)               :: inside
  type(so3_T), pointer           :: obj
  type(r_T)                     :: r

  call c_f_pointer(handle, obj)
  call setRotationPrecision('d')
  r = r_T( rdinp = rod )
  inside = logical(obj%IsinsideFZ(r), c_bool)

end function c_so3_is_inside_fz

!--------------------------------------------------------------------------
! MacKenzie distribution
!--------------------------------------------------------------------------

subroutine c_so3_mackenzie(handle, nsteps, misor, mk) &
    bind(c, name='emsoft_so3_mackenzie')
  !! Compute the theoretical MacKenzie misorientation distribution.
  !! misor(0:nsteps) contains angle values in radians.
  !! mk(0:nsteps) is filled with the distribution values.
  type(c_ptr), value, INTENT(IN)    :: handle
  integer(c_int), value, INTENT(IN) :: nsteps
  real(c_double), INTENT(IN)        :: misor(0:nsteps)
  real(c_double), INTENT(INOUT)     :: mk(0:nsteps)
  type(so3_T), pointer              :: obj

  call c_f_pointer(handle, obj)
  call obj%getMacKenzieDistribution(nsteps, misor, mk)

end subroutine c_so3_mackenzie

end module c_so3
