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

  module mod_KRicosahedral
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/09/26
  !!
  !! Homochoric to Octahedral FZ mapping module
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

  public :: KRicosahedral

  interface KRicosahedral
    procedure KRicosahedral_single
    procedure KRicosahedral_batch
  end interface KRicosahedral
  
  ! 12 five-fold axes (unit, pre-normalized)
  real(real64), parameter :: AXES_5F(3, 12) = reshape([ &
    -0.723606797749979_real64, -0.5257311121191336_real64,  0.4472135954999579_real64, &
    -0.723606797749979_real64,  0.5257311121191336_real64,  0.4472135954999579_real64, &
     0.0_real64,                0.0_real64,                 1.0_real64,              &
     0.276393202250021_real64, -0.8506508083520400_real64,  0.4472135954999579_real64, &
     0.276393202250021_real64,  0.8506508083520400_real64,  0.4472135954999579_real64, &
     0.894427190999916_real64,  0.0_real64,                 0.4472135954999579_real64, &
     0.723606797749979_real64, -0.5257311121191336_real64, -0.4472135954999579_real64, &
    -0.276393202250021_real64,  0.8506508083520400_real64, -0.4472135954999579_real64, &
     0.0_real64,                0.0_real64,                -1.0_real64,              &
    -0.276393202250021_real64, -0.8506508083520400_real64, -0.4472135954999579_real64, &
    -0.894427190999916_real64,  0.0_real64,                -0.4472135954999579_real64, &
     0.723606797749979_real64,  0.5257311121191336_real64, -0.4472135954999579_real64  &
  ], [3, 12])
  
  ! Hard-coded Chebyshev ψ'(u) fit coefficients
  real(real64), parameter :: CHEB_NUM(18) = [ &
     0.3304924856175057_real64, &
     0.31583416591403574_real64, &
    -0.016688974639899247_real64, &
    -0.0016721216544865944_real64, &
     0.0003619870280989233_real64, &
    -3.4011730639492566e-06_real64, &
    -6.302515376491188e-06_real64, &
     6.421149861038183e-07_real64, &
     6.949651306457666e-08_real64, &
    -2.0250514424156133e-08_real64, &
     4.204414792689178e-10_real64, &
     4.1271900706295387e-10_real64, &
    -4.990132935908946e-11_real64, &
    -4.6492593956442474e-12_real64, &
     1.6306230055684583e-12_real64, &
    -5.1885682094719584e-14_real64, &
    -3.3872343958599216e-14_real64, &
     4.631902188019415e-15_real64  &
  ]
  real(real64), parameter :: PSI_MAX = cPi / 5.0_real64  

contains

!--------------------------------------------------------------------------
pure subroutine project_to_tangent(v, a, vt)
 !DEC$ ATTRIBUTES DLLEXPORT :: project_to_tangent
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! Project vector to tangent plane

  IMPLICIT NONE 

  real(real64), INTENT(IN)  :: v(3), a(3)
  real(real64), INTENT(OUT) :: vt(3)
  real(real64)              :: dot, n
  real(real64)              :: vtmp(3)
  
  dot = dot_product(v, a)
  vtmp = v - dot * a
  n = norm2(vtmp)
  if (n > 1.0e-15_real64) then
    vt = vtmp / n
  else
    vt = vtmp
  end if

end subroutine project_to_tangent

!--------------------------------------------------------------------------
pure function theta_max(psi) result(th_max)
!DEC$ ATTRIBUTES DLLEXPORT :: theta_max
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! KR geometry: theta_max

  IMPLICIT NONE 

  real(real64), INTENT(IN)  :: psi
  real(real64)              :: th_max
  real(real64)              :: num, den
  
  num = 1.0_real64 - COS_BETA
  den = SIN_BETA * max(cos(psi), EPS)
  th_max = atan(num / den)

end function theta_max

!--------------------------------------------------------------------------
pure function R3_of_theta(theta) result(R3)
!DEC$ ATTRIBUTES DLLEXPORT :: R3_of_theta
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! KR geometry: R3_of_theta

  real(real64), INTENT(IN)  :: theta
  real(real64)              :: R3
  real(real64)              :: c, term1, term2
  
  c = max(-1.0_real64 + EPS_EPS, min(1.0_real64 - EPS_EPS, cos(theta)))
  term1 = atan(1.0_real64 / (P_CONST * c))
  term2 = (P_CONST * c) / (1.0_real64 + (P_CONST * c) ** 2)
  R3 = 1.5_real64 * (term1 - term2)

end function R3_of_theta

!--------------------------------------------------------------------------
pure function F_cap(theta) result(F)
!DEC$ ATTRIBUTES DLLEXPORT :: F_cap
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! KR geometry: F_cap

  IMPLICIT NONE 

  real(real64), INTENT(IN)  :: theta
  real(real64)              :: F
  real(real64)              :: c
  
  c = max(-1.0_real64 + EPS_EPS, min(1.0_real64 - EPS_EPS, cos(theta)))
  F = 1.5_real64 * (0.5_real64 * cPi * (1.0_real64 - c) - atan(P_CONST) + c * atan(P_CONST * c))

end function F_cap

!--------------------------------------------------------------------------
pure function psi_prime_from_u(u) result(psi)
!DEC$ ATTRIBUTES DLLEXPORT :: psi_prime_from_u
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! ψ'(u) from Chebyshev fit

  IMPLICIT NONE 

  real(real64), INTENT(IN)  :: u
  real(real64)              :: psi
  real(real64)              :: t
  
  t = 2.0_real64 * max(0.0_real64, min(1.0_real64, u)) - 1.0_real64
  psi = chebval_T(t, CHEB_NUM)
  psi = max(0.0_real64, min(PSI_MAX, psi))

end function psi_prime_from_u

!--------------------------------------------------------------------------
pure subroutine fold_azimuth(psi, psi_delta, mirror_val, sector_idx)
!DEC$ ATTRIBUTES DLLEXPORT :: fold_azimuth
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! Azimuth folding

  IMPLICIT NONE 

  real(real64), INTENT(IN)    :: psi
  real(real64), INTENT(OUT)   :: psi_delta, mirror_val
  integer, INTENT(OUT)        :: sector_idx
  real(real64)                :: psi_mod, psi72
  integer                     :: k
  
  psi_mod = modulo(psi, 2.0_real64 * cPi)
  k = int(psi_mod / WEDGE_72)
  psi72 = psi_mod - real(k, real64) * WEDGE_72
  
  if (psi72 > WEDGE_36) then
    psi_delta = WEDGE_72 - psi72
    mirror_val = 1.0_real64
  else
    psi_delta = psi72
    mirror_val = 0.0_real64
  end if
  sector_idx = k

end subroutine fold_azimuth

!--------------------------------------------------------------------------
pure function unfold_azimuth(psi_prime, mirror_val, sector_idx) result(psi)
!DEC$ ATTRIBUTES DLLEXPORT :: unfold_azimuth
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! Azimuth unfolding

  IMPLICIT NONE 

  real(real64), INTENT(IN)  :: psi_prime, mirror_val
  integer, INTENT(IN)       :: sector_idx
  real(real64)              :: psi
  real(real64)              :: psi72, psi_full
  
  if (mirror_val > 0.5_real64) then
    psi72 = WEDGE_72 - psi_prime
  else
    psi72 = psi_prime
  end if
  psi_full = psi72 + real(sector_idx, real64) * WEDGE_72
  psi = modulo(psi_full, 2.0_real64 * cPi)

end function unfold_azimuth

!--------------------------------------------------------------------------
subroutine polar_solve(y, th_max, th_out)
!DEC$ ATTRIBUTES DLLEXPORT :: polar_solve
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! Polar solve using bisection

  IMPLICIT NONE 

  real(real64), INTENT(IN)    :: y, th_max
  real(real64), INTENT(OUT)   :: th_out
  real(real64)                :: lo, hi, mid, val, denom
  integer                     :: i
  
  lo = 0.0_real64
  hi = th_max
  denom = max(F_cap(th_max), EPS)
  
  do i = 1, 100
    mid = 0.5_real64 * (lo + hi)
    val = F_cap(mid) / denom - y
    if (val > 0.0_real64) then
      hi = mid
    else
      lo = mid
    end if
    if ((hi - lo) < 1.0e-14_real64) exit
  end do
  
  th_out = 0.5_real64 * (lo + hi)

end subroutine polar_solve

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Public API functions
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
subroutine KRicosahedral_single(h, h_out)
!DEC$ ATTRIBUTES DLLEXPORT :: KRicosahedral_single
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! Main mapping function

  IMPLICIT NONE 

  real(real64), INTENT(IN)  :: h(3)
  real(real64)              :: h_out(3)
  
  real(real64)              :: rho, uhat(3)
  integer                   :: idx, i
  real(real64)              :: a(3), az(3), e0(3), e1(3)
  real(real64)              :: ua, theta, ue0, ue1, psi
  real(real64)              :: psi_delta, mirror_val
  integer                   :: sector_idx
  real(real64)              :: u, psi_prime
  real(real64)              :: th_ceiling, y_src, theta_prime
  real(real64)              :: R3, R_cap, rho_prime
  real(real64)              :: psi_full, sin_t, xloc, yloc, zloc
  real(real64)              :: dots(12), max_dot
  logical                   :: az_posZ, az_negZ
  real(real64)              :: sign_z
  
  ! Compute rho and unit vector
  rho = norm2(h)
  if (rho > EPS) then
    uhat = h / rho
  else
    uhat = 0.0_real64
    h_out = 0.0_real64
    return
  end if
  
  ! Find nearest face axis (max dot product)
  max_dot = -huge(1.0_real64)
  idx = 1
  do i = 1, 12
    dots(i) = dot_product(uhat, AXES_5F(:, i))
    if (dots(i) > max_dot) then
      max_dot = dots(i)
      idx = i
    end if
  end do
  a = AXES_5F(:, idx)
  
  ! Azimuth-0 direction
  az_posZ = (1.0_real64 - a(3)) < 1.0e-12_real64  ! a ≈ +Z
  az_negZ = (1.0_real64 + a(3)) < 1.0e-12_real64  ! a ≈ -Z
  
  if (az_posZ) then
    az = [1.0_real64, 0.0_real64, 0.0_real64]  ! +Z face → +X
  else if (az_negZ) then
    az = [-1.0_real64, 0.0_real64, 0.0_real64]  ! -Z face → -X
  else
    sign_z = sign(1.0_real64, a(3))
    az = sign_z * [0.0_real64, 0.0_real64, 1.0_real64]  ! z>0 → +Z, z<0 → -Z
  end if
  
  ! Tangent basis on face: e0 = proj(az, ⟂ a), e1 = a × e0
  call project_to_tangent(az, a, e0)
  e1(1) = a(2) * e0(3) - a(3) * e0(2)
  e1(2) = a(3) * e0(1) - a(1) * e0(3)
  e1(3) = a(1) * e0(2) - a(2) * e0(1)
  
  ! Normalize e1
  e1 = e1 / norm2(e1)
  
  ! Local spherical: theta = arccos(uhat·a), psi = atan2(uhat·e1, uhat·e0)
  ua = max(-1.0_real64, min(1.0_real64, dot_product(uhat, a)))
  theta = acos(ua)
  ue0 = dot_product(uhat, e0)
  ue1 = dot_product(uhat, e1)
  psi = atan2(ue1, ue0)
  psi = modulo(psi, 2.0_real64 * cPi)
  
  ! KR azimuth
  call fold_azimuth(psi, psi_delta, mirror_val, sector_idx)
  u = psi_delta / WEDGE_36
  psi_prime = psi_prime_from_u(u)
  
  ! KR polar
  th_ceiling = theta_max(psi_prime)
  y_src = (1.0_real64 - cos(theta)) / max(1.0_real64 - cos(th_ceiling), EPS)
  y_src = max(0.0_real64, min(1.0_real64, y_src))
  
  call polar_solve(y_src, th_ceiling, theta_prime)
  
  ! Radial KR scale; slight retract
  R3 = R3_of_theta(theta_prime)
  R_cap = max(R3, 0.0_real64) ** (1.0_real64 / 3.0_real64)
  rho_prime = (1.0_real64 - 1.0e-12_real64) * rho * (R_cap / H_MAX)
  
  ! Rebuild homochoric in global xyz using (e0,e1,a) local frame
  psi_full = unfold_azimuth(psi_prime, mirror_val, sector_idx)
  sin_t = sin(theta_prime)
  xloc = rho_prime * sin_t * cos(psi_full)
  yloc = rho_prime * sin_t * sin(psi_full)
  zloc = rho_prime * cos(theta_prime)
  h_out = xloc * e0 + yloc * e1 + zloc * a

end subroutine KRicosahedral_single

! Batch processing function with OpenMP support
!--------------------------------------------------------------------------
subroutine KRicosahedral_batch(h_in, h_out, n)
!DEC$ ATTRIBUTES DLLEXPORT :: KRicosahedral_batch
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! Batch processing function with OpenMP support

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: h_in(3, n)
  real(real64), INTENT(OUT) :: h_out(3, n)
  integer, INTENT(IN)       :: n
  integer                   :: i
  
  !$omp parallel do default(none) shared(h_in, h_out, n) private(i)
  do i = 1, n
    call KRicosahedral_single( h_in(:, i), h_out(:, i) )
  end do
  !$omp end parallel do

end subroutine KRicosahedral_batch

end module mod_KRicosahedral

