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

module c_diffraction
  !! author: MDG
  !! version: 1.0
  !! date: 03/27/26
  !!
  !! C-interop wrappers for Diffraction_T.
  !! Provides voltage/wavelength setup and structure factor calculations.

use iso_c_binding
use mod_kinds
use mod_global
use mod_diffraction
use mod_crystallography

IMPLICIT NONE

private

contains

!--------------------------------------------------------------------------
! Constructor / destructor
!--------------------------------------------------------------------------

function c_diff_create(voltage_kv) result(handle) bind(c, name='emsoft_diff_create')
  !! Create a Diffraction_T with the given accelerating voltage (keV).
  !! Computes relativistic correction, wavelength, and interaction constant.
  real(c_double), value, INTENT(IN) :: voltage_kv
  type(c_ptr)                       :: handle
  type(Diffraction_T), pointer      :: obj

  allocate(obj)
  obj = Diffraction_T()
  call obj%setV(voltage_kv)
  handle = c_loc(obj)

end function c_diff_create

!--------------------------------------------------------------------------
subroutine c_diff_destroy(handle) bind(c, name='emsoft_diff_destroy')
  type(c_ptr), value, INTENT(IN) :: handle
  type(Diffraction_T), pointer   :: obj

  call c_f_pointer(handle, obj)
  deallocate(obj)

end subroutine c_diff_destroy

!--------------------------------------------------------------------------
! Compute wavelength (requires a Cell_T for V0mod correction)
!--------------------------------------------------------------------------

subroutine c_diff_calc_wavelength(handle, cell_handle) &
    bind(c, name='emsoft_diff_calc_wavelength')
  !! Compute the relativistic electron wavelength using the cell's mean
  !! inner potential for refraction correction.
  type(c_ptr), value, INTENT(IN) :: handle
  type(c_ptr), value, INTENT(IN) :: cell_handle
  type(Diffraction_T), pointer   :: obj
  type(Cell_T), pointer          :: cell

  call c_f_pointer(handle, obj)
  call c_f_pointer(cell_handle, cell)
  call obj%CalcWaveLength(cell)

end subroutine c_diff_calc_wavelength

!--------------------------------------------------------------------------
! Property getters
!--------------------------------------------------------------------------

function c_diff_get_voltage(handle) result(v) bind(c, name='emsoft_diff_get_voltage')
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double)                 :: v
  type(Diffraction_T), pointer   :: obj

  call c_f_pointer(handle, obj)
  v = obj%getV()

end function c_diff_get_voltage

!--------------------------------------------------------------------------
function c_diff_get_wavelength(handle) result(v) bind(c, name='emsoft_diff_get_wavelength')
  !! Get electron wavelength in nm.
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double)                 :: v
  type(Diffraction_T), pointer   :: obj

  call c_f_pointer(handle, obj)
  v = obj%getWaveLength()

end function c_diff_get_wavelength

!--------------------------------------------------------------------------
function c_diff_get_relcor(handle) result(v) bind(c, name='emsoft_diff_get_relcor')
  !! Get relativistic correction factor (gamma).
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double)                 :: v
  type(Diffraction_T), pointer   :: obj

  call c_f_pointer(handle, obj)
  v = obj%getRelcor()

end function c_diff_get_relcor

!--------------------------------------------------------------------------
function c_diff_get_sigma(handle) result(v) bind(c, name='emsoft_diff_get_sigma')
  !! Get interaction constant (V^-1 nm^-1).
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double)                 :: v
  type(Diffraction_T), pointer   :: obj

  call c_f_pointer(handle, obj)
  v = obj%getSigma()

end function c_diff_get_sigma

!--------------------------------------------------------------------------
function c_diff_get_psihat(handle) result(v) bind(c, name='emsoft_diff_get_psihat')
  !! Get relativistically corrected accelerating potential (V).
  type(c_ptr), value, INTENT(IN) :: handle
  real(c_double)                 :: v
  type(Diffraction_T), pointer   :: obj

  call c_f_pointer(handle, obj)
  v = obj%getPsihat()

end function c_diff_get_psihat

!--------------------------------------------------------------------------
! Scattering method
!--------------------------------------------------------------------------

subroutine c_diff_set_method(handle, m1, m2) bind(c, name='emsoft_diff_set_method')
  !! Set scattering factor method: 'WK', 'DT', or 'XR'.
  type(c_ptr), value, INTENT(IN) :: handle
  character(c_char), value, INTENT(IN) :: m1, m2
  type(Diffraction_T), pointer   :: obj
  character(2)                   :: method

  call c_f_pointer(handle, obj)
  method(1:1) = m1
  method(2:2) = m2
  call obj%setrlpmethod(method)

end subroutine c_diff_set_method

!--------------------------------------------------------------------------
! Structure factor calculation
!--------------------------------------------------------------------------

subroutine c_diff_calc_ucg(handle, cell_handle, hkl, xg, xgp, ucg_r, ucg_i, vphase) &
    bind(c, name='emsoft_diff_calc_ucg')
  !! Compute structure factor for reflection hkl.
  !! Cell must have atom positions calculated (via calcPositions).
  !! Returns extinction distance, absorption length, Ucg components, and phase.
  type(c_ptr), value, INTENT(IN)    :: handle
  type(c_ptr), value, INTENT(IN)    :: cell_handle
  integer(c_int), INTENT(IN)        :: hkl(3)
  real(c_double), INTENT(OUT)       :: xg
  real(c_double), INTENT(OUT)       :: xgp
  real(c_double), INTENT(OUT)       :: ucg_r, ucg_i
  real(c_double), INTENT(OUT)       :: vphase
  type(Diffraction_T), pointer      :: obj
  type(Cell_T), pointer             :: cell
  type(gnode)                       :: rlp

  call c_f_pointer(handle, obj)
  call c_f_pointer(cell_handle, cell)
  call obj%CalcUcg(cell, hkl)
  rlp = obj%getrlp()
  xg = dble(rlp%xg)
  xgp = dble(rlp%xgp)
  ucg_r = dble(real(rlp%Ucg))
  ucg_i = dble(aimag(rlp%Ucg))
  vphase = dble(rlp%Vphase)

end subroutine c_diff_calc_ucg

end module c_diffraction
