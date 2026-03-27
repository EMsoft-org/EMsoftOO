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

module c_crystallography
  !! author: MDG
  !! version: 1.0
  !! date: 03/27/26
  !!
  !! C-interop wrappers for Cell_T (crystallography).

use iso_c_binding
use mod_kinds
use mod_global
use mod_crystallography
use mod_symmetry

IMPLICIT NONE

private

contains

!--------------------------------------------------------------------------
! Constructor / destructor
!--------------------------------------------------------------------------

function c_cell_create(latparm) result(handle) bind(c, name='emsoft_cell_create')
  !! Create a Cell_T from lattice parameters [a, b, c, alpha, beta, gamma].
  !! Lengths in nm, angles in degrees. Automatically computes metric tensors.
  real(c_double), INTENT(IN) :: latparm(6)
  type(c_ptr)                :: handle
  type(Cell_T), pointer      :: obj

  allocate(obj)
  obj = Cell_T( latparm = latparm )
  call obj%CalcMatrices()
  handle = c_loc(obj)

end function c_cell_create

!--------------------------------------------------------------------------
subroutine c_cell_destroy(handle) bind(c, name='emsoft_cell_destroy')
  type(c_ptr), value, INTENT(IN) :: handle
  type(Cell_T), pointer          :: obj

  call c_f_pointer(handle, obj)
  deallocate(obj)

end subroutine c_cell_destroy

!--------------------------------------------------------------------------
! Lattice parameter access
!--------------------------------------------------------------------------

subroutine c_cell_get_latparm(handle, latparm) bind(c, name='emsoft_cell_get_latparm')
  !! Get lattice parameters [a, b, c, alpha, beta, gamma].
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double), INTENT(OUT)    :: latparm(6)
  type(Cell_T), pointer          :: obj

  call c_f_pointer(handle, obj)
  latparm = obj%getLatParm()

end subroutine c_cell_get_latparm

!--------------------------------------------------------------------------
function c_cell_get_volume(handle) result(vol) bind(c, name='emsoft_cell_get_volume')
  !! Get unit cell volume in nm^3.
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double)                 :: vol
  type(Cell_T), pointer          :: obj

  call c_f_pointer(handle, obj)
  vol = obj%getVolume()

end function c_cell_get_volume

!--------------------------------------------------------------------------
! Metric tensors and structure matrices
!--------------------------------------------------------------------------

subroutine c_cell_get_dmt(handle, dmt) bind(c, name='emsoft_cell_get_dmt')
  !! Get direct metric tensor (3x3).
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double), INTENT(OUT)    :: dmt(3,3)
  type(Cell_T), pointer          :: obj

  call c_f_pointer(handle, obj)
  dmt = obj%getdmt()

end subroutine c_cell_get_dmt

!--------------------------------------------------------------------------
subroutine c_cell_get_rmt(handle, rmt) bind(c, name='emsoft_cell_get_rmt')
  !! Get reciprocal metric tensor (3x3).
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double), INTENT(OUT)    :: rmt(3,3)
  type(Cell_T), pointer          :: obj

  call c_f_pointer(handle, obj)
  rmt = obj%getrmt()

end subroutine c_cell_get_rmt

!--------------------------------------------------------------------------
subroutine c_cell_get_dsm(handle, dsm) bind(c, name='emsoft_cell_get_dsm')
  !! Get direct structure matrix (3x3).
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double), INTENT(OUT)    :: dsm(3,3)
  type(Cell_T), pointer          :: obj

  call c_f_pointer(handle, obj)
  dsm = obj%getdsm()

end subroutine c_cell_get_dsm

!--------------------------------------------------------------------------
subroutine c_cell_get_rsm(handle, rsm) bind(c, name='emsoft_cell_get_rsm')
  !! Get reciprocal structure matrix (3x3).
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double), INTENT(OUT)    :: rsm(3,3)
  type(Cell_T), pointer          :: obj

  call c_f_pointer(handle, obj)
  rsm = obj%getrsm()

end subroutine c_cell_get_rsm

!--------------------------------------------------------------------------
! Crystallographic computations
!--------------------------------------------------------------------------

subroutine c_cell_trans_space(handle, t, d, inspace, outspace) &
    bind(c, name='emsoft_cell_trans_space')
  !! Transform a 3-vector between coordinate systems.
  !! Spaces: 'd' (direct), 'r' (reciprocal), 'c' (Cartesian).
  type(c_ptr), value, INTENT(IN)  :: handle
  real(c_double), INTENT(IN)      :: t(3)
  real(c_double), INTENT(OUT)     :: d(3)
  character(c_char), value, INTENT(IN) :: inspace
  character(c_char), value, INTENT(IN) :: outspace
  type(Cell_T), pointer           :: obj

  call c_f_pointer(handle, obj)
  call obj%TransSpace(t, d, inspace, outspace)

end subroutine c_cell_trans_space

!--------------------------------------------------------------------------
function c_cell_calc_dot(handle, p, q, space) result(cdot) &
    bind(c, name='emsoft_cell_calc_dot')
  !! Compute dot product of two vectors in the given space.
  type(c_ptr), value, INTENT(IN)  :: handle
  real(c_double), INTENT(IN)      :: p(3), q(3)
  character(c_char), value, INTENT(IN) :: space
  real(c_double)                  :: cdot
  type(Cell_T), pointer           :: obj

  call c_f_pointer(handle, obj)
  cdot = obj%calcDot(p, q, space)

end function c_cell_calc_dot

!--------------------------------------------------------------------------
function c_cell_calc_length(handle, p, space) result(length) &
    bind(c, name='emsoft_cell_calc_length')
  !! Compute length of a vector in the given space.
  type(c_ptr), value, INTENT(IN)  :: handle
  real(c_double), INTENT(IN)      :: p(3)
  character(c_char), value, INTENT(IN) :: space
  real(c_double)                  :: length
  type(Cell_T), pointer           :: obj

  call c_f_pointer(handle, obj)
  length = obj%calcLength(p, space)

end function c_cell_calc_length

!--------------------------------------------------------------------------
function c_cell_calc_angle(handle, p, q, space) result(angle) &
    bind(c, name='emsoft_cell_calc_angle')
  !! Compute angle (radians) between two vectors in the given space.
  type(c_ptr), value, INTENT(IN)  :: handle
  real(c_double), INTENT(IN)      :: p(3), q(3)
  character(c_char), value, INTENT(IN) :: space
  real(c_double)                  :: angle
  type(Cell_T), pointer           :: obj

  call c_f_pointer(handle, obj)
  angle = obj%calcAngle(p, q, space)

end function c_cell_calc_angle

!--------------------------------------------------------------------------
subroutine c_cell_norm_vec(handle, p, space) bind(c, name='emsoft_cell_norm_vec')
  !! Normalize a vector in-place in the given space.
  type(c_ptr), value, INTENT(IN)  :: handle
  real(c_double), INTENT(INOUT)   :: p(3)
  character(c_char), value, INTENT(IN) :: space
  type(Cell_T), pointer           :: obj

  call c_f_pointer(handle, obj)
  call obj%NormVec(p, space)

end subroutine c_cell_norm_vec

!--------------------------------------------------------------------------
subroutine c_cell_calc_cross(handle, p, q, r, inspace, outspace) &
    bind(c, name='emsoft_cell_calc_cross')
  !! Compute cross product of two vectors.
  type(c_ptr), value, INTENT(IN)  :: handle
  real(c_double), INTENT(IN)      :: p(3), q(3)
  real(c_double), INTENT(OUT)     :: r(3)
  character(c_char), value, INTENT(IN) :: inspace, outspace
  type(Cell_T), pointer           :: obj

  call c_f_pointer(handle, obj)
  call obj%calcCross(p, q, r, inspace, outspace, 0)

end subroutine c_cell_calc_cross

!--------------------------------------------------------------------------
! Atom setup (programmatic crystal definition)
!--------------------------------------------------------------------------

subroutine c_cell_set_xtal_system(handle, xs) bind(c, name='emsoft_cell_set_xtal_system')
  !! Set the crystal system number (1=cubic..7=triclinic).
  type(c_ptr), value, INTENT(IN)    :: handle
  integer(c_int), value, INTENT(IN) :: xs
  type(Cell_T), pointer             :: obj

  call c_f_pointer(handle, obj)
  call obj%setXtalSystem(xs)

end subroutine c_cell_set_xtal_system

!--------------------------------------------------------------------------
subroutine c_cell_set_natomtype(handle, n) bind(c, name='emsoft_cell_set_natomtype')
  !! Set the number of atom types in the asymmetric unit.
  type(c_ptr), value, INTENT(IN)    :: handle
  integer(c_int), value, INTENT(IN) :: n
  type(Cell_T), pointer             :: obj

  call c_f_pointer(handle, obj)
  call obj%setNatomtype(n)

end subroutine c_cell_set_natomtype

!--------------------------------------------------------------------------
subroutine c_cell_setup_atoms(handle, sg_handle, natom, atomtypes, atomdata) &
    bind(c, name='emsoft_cell_setup_atoms')
  !! Set up all atoms in the asymmetric unit and compute equivalent positions.
  !! atomtypes(natom): atomic numbers (e.g. 28 for Ni).
  !! atomdata(natom, 5): each row is [x, y, z, occupancy, Debye-Waller].
  !! Calls calcPositions internally so the cell is ready for diffraction.
  type(c_ptr), value, INTENT(IN)    :: handle
  type(c_ptr), value, INTENT(IN)    :: sg_handle
  integer(c_int), value, INTENT(IN) :: natom
  integer(c_int), INTENT(IN)        :: atomtypes(natom)
  real(c_double), INTENT(IN)        :: atomdata(natom, 5)
  type(Cell_T), pointer             :: obj
  type(SpaceGroup_T), pointer       :: sg
  real(kind=dbl)                    :: pos(maxpasym, 5)
  integer(kind=irg)                 :: i

  call c_f_pointer(handle, obj)
  call c_f_pointer(sg_handle, sg)

  call obj%setNatomtype(natom)
  call obj%setXtalSystem(sg%getSpaceGroupXtalSystem())

  ! Set atom types
  call obj%setAtomtype(atomtypes(1:natom))

  ! Set atom positions via the full array interface
  pos = 0.D0
  do i = 1, natom
    pos(i, 1:5) = atomdata(i, 1:5)
  end do
  call obj%setAtomPos(pos)

  ! Generate all equivalent positions
  call obj%calcPositions(sg, 'v')

end subroutine c_cell_setup_atoms

end module c_crystallography
