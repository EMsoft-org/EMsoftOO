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

module c_symmetry
  !! author: MDG
  !! version: 1.0
  !! date: 03/27/26
  !!
  !! C-interop wrappers for SpaceGroup_T (symmetry operations).

use iso_c_binding
use mod_kinds
use mod_global
use mod_symmetry

IMPLICIT NONE

private

contains

!--------------------------------------------------------------------------
! Constructor / destructor
!--------------------------------------------------------------------------

function c_sg_create(sgnum) result(handle) bind(c, name='emsoft_sg_create')
  !! Create a SpaceGroup_T from a space group number (1-230).
  !! Generates the symmetry matrices including point group operators.
  integer(c_int), value, INTENT(IN) :: sgnum
  type(c_ptr)                        :: handle
  type(SpaceGroup_T), pointer        :: obj

  allocate(obj)
  obj = SpaceGroup_T( SGnumber = sgnum )
  handle = c_loc(obj)

end function c_sg_create

!--------------------------------------------------------------------------
subroutine c_sg_destroy(handle) bind(c, name='emsoft_sg_destroy')
  type(c_ptr), value, INTENT(IN) :: handle
  type(SpaceGroup_T), pointer    :: obj

  call c_f_pointer(handle, obj)
  deallocate(obj)

end subroutine c_sg_destroy

!--------------------------------------------------------------------------
! Property getters
!--------------------------------------------------------------------------

function c_sg_get_number(handle) result(n) bind(c, name='emsoft_sg_get_number')
  type(c_ptr), value, INTENT(IN) :: handle
  integer(c_int)                 :: n
  type(SpaceGroup_T), pointer    :: obj

  call c_f_pointer(handle, obj)
  n = obj%getSpaceGroupNumber()

end function c_sg_get_number

!--------------------------------------------------------------------------
function c_sg_get_order(handle) result(n) bind(c, name='emsoft_sg_get_order')
  !! Get the order of the space group (number of symmetry operations).
  type(c_ptr), value, INTENT(IN) :: handle
  integer(c_int)                 :: n
  type(SpaceGroup_T), pointer    :: obj

  call c_f_pointer(handle, obj)
  n = obj%getSpaceGroupOrder()

end function c_sg_get_order

!--------------------------------------------------------------------------
function c_sg_get_matnum(handle) result(n) bind(c, name='emsoft_sg_get_matnum')
  !! Get the number of symmetry matrices.
  type(c_ptr), value, INTENT(IN) :: handle
  integer(c_int)                 :: n
  type(SpaceGroup_T), pointer    :: obj

  call c_f_pointer(handle, obj)
  n = obj%getSpaceGroupMATnum()

end function c_sg_get_matnum

!--------------------------------------------------------------------------
function c_sg_get_numpt(handle) result(n) bind(c, name='emsoft_sg_get_numpt')
  !! Get the number of point group operators.
  type(c_ptr), value, INTENT(IN) :: handle
  integer(c_int)                 :: n
  type(SpaceGroup_T), pointer    :: obj

  call c_f_pointer(handle, obj)
  n = obj%getSpaceGroupNUMpt()

end function c_sg_get_numpt

!--------------------------------------------------------------------------
function c_sg_get_xtal_system(handle) result(n) bind(c, name='emsoft_sg_get_xtal_system')
  !! Get crystal system number (1-7).
  type(c_ptr), value, INTENT(IN) :: handle
  integer(c_int)                 :: n
  type(SpaceGroup_T), pointer    :: obj

  call c_f_pointer(handle, obj)
  n = obj%getSpaceGroupXtalSystem()

end function c_sg_get_xtal_system

!--------------------------------------------------------------------------
function c_sg_get_centro(handle) result(c) bind(c, name='emsoft_sg_get_centro')
  !! Check if the space group is centrosymmetric.
  type(c_ptr), value, INTENT(IN) :: handle
  logical(c_bool)                :: c
  type(SpaceGroup_T), pointer    :: obj

  call c_f_pointer(handle, obj)
  c = logical(obj%getSpaceGroupCentro(), c_bool)

end function c_sg_get_centro

!--------------------------------------------------------------------------
function c_sg_get_symmorphic(handle) result(s) bind(c, name='emsoft_sg_get_symmorphic')
  !! Check if the space group is symmorphic.
  type(c_ptr), value, INTENT(IN) :: handle
  logical(c_bool)                :: s
  type(SpaceGroup_T), pointer    :: obj

  call c_f_pointer(handle, obj)
  s = logical(obj%getSpaceGroupSymmorphic(), c_bool)

end function c_sg_get_symmorphic

!--------------------------------------------------------------------------
function c_sg_is_g_allowed(handle, g) result(allowed) bind(c, name='emsoft_sg_is_g_allowed')
  !! Check if a reflection g = [h, k, l] is allowed (not extinct).
  type(c_ptr), value, INTENT(IN) :: handle
  integer(c_int), INTENT(IN)    :: g(3)
  logical(c_bool)               :: allowed
  type(SpaceGroup_T), pointer   :: obj

  call c_f_pointer(handle, obj)
  allowed = logical(obj%IsGAllowed(g), c_bool)

end function c_sg_is_g_allowed

!--------------------------------------------------------------------------
! Symmetry computations
!--------------------------------------------------------------------------

subroutine c_sg_calc_orbit(handle, site, n, ctmp, maxn) &
    bind(c, name='emsoft_sg_calc_orbit')
  !! Compute the orbit of a position. Returns n equivalent positions in ctmp.
  !! ctmp must be pre-allocated to at least (maxn, 3).
  type(c_ptr), value, INTENT(IN)    :: handle
  real(c_double), INTENT(IN)        :: site(3)
  integer(c_int), INTENT(OUT)       :: n
  integer(c_int), value, INTENT(IN) :: maxn
  real(c_double), INTENT(OUT)       :: ctmp(maxn, 3)
  type(SpaceGroup_T), pointer       :: obj
  real(kind=dbl), allocatable       :: work(:,:)
  integer(kind=irg)                 :: nout, i

  call c_f_pointer(handle, obj)
  call obj%CalcOrbit(site, nout, work)
  n = nout
  ctmp = 0.D0
  do i = 1, min(nout, maxn)
    ctmp(i, 1:3) = work(i, 1:3)
  end do
  if (allocated(work)) deallocate(work)

end subroutine c_sg_calc_orbit

!--------------------------------------------------------------------------
subroutine c_sg_calc_star(handle, kk, n, stmp, space, maxn) &
    bind(c, name='emsoft_sg_calc_star')
  !! Compute the star of a reciprocal/direct vector.
  !! stmp must be pre-allocated to at least (maxn, 3).
  type(c_ptr), value, INTENT(IN)       :: handle
  real(c_double), INTENT(IN)           :: kk(3)
  integer(c_int), INTENT(OUT)          :: n
  integer(c_int), value, INTENT(IN)    :: maxn
  real(c_double), INTENT(OUT)          :: stmp(maxn, 3)
  character(c_char), value, INTENT(IN) :: space
  type(SpaceGroup_T), pointer          :: obj
  real(kind=dbl), allocatable          :: work(:,:)
  integer(kind=irg)                    :: nout, i

  call c_f_pointer(handle, obj)
  call obj%CalcStar(kk, nout, work, space)
  n = nout
  stmp = 0.D0
  do i = 1, min(nout, maxn)
    stmp(i, 1:3) = work(i, 1:3)
  end do
  if (allocated(work)) deallocate(work)

end subroutine c_sg_calc_star

!--------------------------------------------------------------------------
subroutine c_sg_calc_family(handle, ind, num, itmp, space, maxn) &
    bind(c, name='emsoft_sg_calc_family')
  !! Compute the family of symmetry-equivalent planes/directions.
  !! itmp must be pre-allocated to at least (maxn, 3).
  type(c_ptr), value, INTENT(IN)       :: handle
  integer(c_int), INTENT(IN)           :: ind(3)
  integer(c_int), INTENT(OUT)          :: num
  integer(c_int), value, INTENT(IN)    :: maxn
  integer(c_int), INTENT(OUT)          :: itmp(maxn, 3)
  character(c_char), value, INTENT(IN) :: space
  type(SpaceGroup_T), pointer          :: obj
  integer(kind=irg), allocatable       :: work(:,:)
  integer(kind=irg)                    :: nout, i

  call c_f_pointer(handle, obj)
  call obj%CalcFamily(ind, nout, space, work)
  num = nout
  itmp = 0
  do i = 1, min(nout, maxn)
    itmp(i, 1:3) = work(i, 1:3)
  end do
  if (allocated(work)) deallocate(work)

end subroutine c_sg_calc_family

end module c_symmetry
