! ###################################################################
! Copyright (c) 2014-2026, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_KRtetrahedral
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! Homochoric to Tetrahedral FZ mapping module
  !! 
  !!
  !! This module provides vectorized mapping functions for batch processing.
  !! For OpenMP parallelization, use the batch functions with appropriate
  !! compiler flags: -fopenmp (GCC) or -openmp (Intel)

  use mod_kinds
  use mod_global
  use mod_KRsupport
  
  use, intrinsic :: iso_fortran_env, only: real64

  IMPLICIT NONE
  
  private

  public :: KRtetrahedral

  interface KRtetrahedral
    procedure KRtetrahedral_single
    procedure KRtetrahedral_batch
  end interface KRtetrahedral
  
  ! Padé numerator coefficients (n=12)
  real(real64), parameter :: A_COEFFS(13) = [ &
    0.36781788207117677_real64, 0.3947242675359163_real64, 0.025552232221354085_real64, &
    -0.0019925177097101683_real64, -0.0006826585228387496_real64, -3.461937624041875e-05_real64, &
    1.0711964436465566e-05_real64, 8.266471280731193e-07_real64, 4.637247599777716e-07_real64, &
    6.556999583536072e-07_real64, -1.4219967302931709e-08_real64, 4.576707508498134e-09_real64, &
    1.3484067046460952e-10_real64]
  
  ! Padé denominator coefficients (m=10, b[0]=1 fixed)
  real(real64), parameter :: B_COEFFS(11) = [ &
    1.0_real64, 2.9141418502709427e-06_real64, -1.5269094632305715e-06_real64, &
    6.416353453852009e-07_real64, 1.6615935995742397e-06_real64, -1.5174845503286455e-06_real64, &
    -1.8299029656639242e-06_real64, -6.052492204511186e-06_real64, 4.967030741589562e-06_real64, &
    -5.04089587435592e-07_real64, 6.408342665323226e-08_real64]
  
contains

!--------------------------------------------------------------------------
pure function c_T(theta, phi) result(c)
!DEC$ ATTRIBUTES DLLEXPORT :: c_T
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! T support function

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: theta, phi
  real(real64)              :: c
  
  c = cos(theta) + sin(theta) * (cos(phi) + sin(phi))

end function c_T

!--------------------------------------------------------------------------
function A_phi_T(phi) result(A)
!DEC$ ATTRIBUTES DLLEXPORT :: A_phi_T
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! Column area A(φ) = ∫_0^{π/2} ρ^3(c_T(θ, φ)) sinθ dθ

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: phi
  real(real64)              :: A
  integer                   :: i
  real(real64)              :: xa, wa, th_eval, f_val
  
  A = 0.0_real64
  do i = 1, 16
    xa = 0.5_real64 * (GL16_X(i) * HALF_PI + HALF_PI)
    wa = 0.5_real64 * HALF_PI * GL16_W(i)
    th_eval = xa
    f_val = rho3_from_c(c_T(th_eval, phi)) * sin(th_eval)
    A = A + wa * f_val
  end do
  A = max(A, EPS_SLOPE)
end function A_phi_T

!--------------------------------------------------------------------------
function G_theta_phi_T(theta, phi) result(G)
!DEC$ ATTRIBUTES DLLEXPORT :: G_theta_phi_T
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! G(θ; φ) = ∫_0^{θ} ρ^3(c_T(τ, φ)) sinτ dτ

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: theta, phi
  real(real64)              :: G
  integer                   :: i
  real(real64)              :: xa, wa, th_eval, f_val
  
  G = 0.0_real64
  do i = 1, 16
    xa = 0.5_real64 * (GL16_X(i) * theta + theta)
    wa = 0.5_real64 * theta * GL16_W(i)
    th_eval = xa
    f_val = rho3_from_c(c_T(th_eval, phi)) * sin(th_eval)
    G = G + wa * f_val
  end do

end function G_theta_phi_T

  ! Padé inverse CDF for T
!--------------------------------------------------------------------------
pure function phi_inv_T(u) result(phi)
!DEC$ ATTRIBUTES DLLEXPORT :: omega_max_from_cos
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! Padé inverse CDF for T

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: u
  real(real64)              :: phi
  real(real64)              :: t, P, Q
  real(real64)              :: d_coeffs(size(B_COEFFS))
  
  t = max(-1.0_real64, min(1.0_real64, 2.0_real64 * u - 1.0_real64))
  P = chebval_T(t, A_COEFFS)
  
  if (size(B_COEFFS) == 1) then
    Q = 1.0_real64
  else
    d_coeffs = B_COEFFS
    d_coeffs(1) = 0.0_real64
    Q = 1.0_real64 + chebval_T(t, d_coeffs)
  end if
  
  phi = max(0.0_real64, min(PHI_MAX_T, P / Q))

end function phi_inv_T

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Public API functions
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
subroutine KRtetrahedral_single(ho, mapped, newton_iters)
 !DEC$ ATTRIBUTES DLLEXPORT :: KRtetrahedral_single
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! Main mapping function

IMPLICIT NONE

real(real64), INTENT(IN)      :: ho(3)
real(real64), INTENT(OUT)     :: mapped(3)
integer, INTENT(IN), OPTIONAL :: newton_iters

real(real64)                  :: v(3), v_sw(3)
real(real64)                  :: sgn(3)
logical                       :: mask_xy
integer                       :: niter, i
real(real64)                  :: r, th, ph
real(real64)                  :: u, ph_p
real(real64)                  :: y, A, th_p
real(real64)                  :: G, Gn, f_theta, fn, step, th_new
real(real64)                  :: big_mask
real(real64)                  :: rho, r_p

niter = 8
if (present(newton_iters)) niter = newton_iters

v = ho

! Fold to first octant and enforce x >= y
sgn = sign(1.0_real64, v)
where (sgn == 0.0_real64) sgn = 1.0_real64
v = abs(v)

mask_xy = v(2) > v(1)
if (mask_xy) then
  v_sw = v
  v_sw(1) = v(2)
  v_sw(2) = v(1)
  v = v_sw
end if

! Spherical
call sph_from_cart_vec(v, r, th, ph)

! Azimuth reparameterization via Padé inverse CDF
u = max(0.0_real64, min(1.0_real64, ph / PHI_MAX_T))
ph_p = phi_inv_T(u)

! Polar via safeguarded Newton on y = (G/A)(θ'; φ')
y = 1.0_real64 - cos(th)
A = A_phi_T(ph_p)
th_p = max(0.0_real64, min(th, HALF_PI))

do i = 1, niter
  G = G_theta_phi_T(th_p, ph_p)
  Gn = G / A
  
  ! Slope at θ': f(θ') / A
  f_theta = max(rho3_from_c(c_T(th_p, ph_p)) * sin(th_p), EPS_SLOPE)
  fn = f_theta / A
  step = (Gn - y) / fn
  th_new = max(0.0_real64, min(HALF_PI, th_p - step))
  
  ! Damping for very large steps
  big_mask = merge(1.0_real64, 0.0_real64, abs(step) > 0.25_real64)
  if (big_mask > 0.5_real64) then
    th_new = 0.75_real64 * th_p + 0.25_real64 * th_new
  end if
  
  ! Early break if converged
  if (abs(th_new - th_p) < 1.0e-12_real64) then
    th_p = th_new
    exit
  end if
  th_p = th_new
end do

! Radial
rho = R_from_c(c_T(th_p, ph_p))
r_p = rho * (r / H_MAX)
mapped = cart_from_sph(r_p, th_p, ph_p)

! Undo y/x swap and signs
if (mask_xy) then
  v_sw = mapped
  v_sw(1) = mapped(2)
  v_sw(2) = mapped(1)
  mapped = v_sw
end if
mapped = mapped * sgn

end subroutine KRtetrahedral_single

!--------------------------------------------------------------------------
subroutine KRtetrahedral_batch(h_in, h_out, n, newton_iters)
!DEC$ ATTRIBUTES DLLEXPORT :: KRtetrahedral_batch
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! Batch processing function with OpenMP support

  IMPLICIT NONE

  real(real64), INTENT(IN)      :: h_in(3, n)
  real(real64), INTENT(OUT)     :: h_out(3, n)
  integer, INTENT(IN)           :: n
  integer, INTENT(IN), OPTIONAL :: newton_iters
  integer                       :: i
  
!$omp parallel do default(none) shared(h_in, h_out, n, newton_iters) private(i)
  do i = 1, n
    if (present(newton_iters)) then
      call KRtetrahedral_single(h_in(:, i), h_out(:, i), newton_iters)
    else
      call KRtetrahedral_single(h_in(:, i), h_out(:, i))
    end if
  end do
!$omp end parallel do

end subroutine KRtetrahedral_batch

end module mod_KRtetrahedral

