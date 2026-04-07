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

module mod_KRdihedral
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! Homochoric to Dihedral FZ mapping module
  !! Supports D2, D3, D4, D6 symmetries
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

  public :: KRdihedral

  interface KRdihedral
    procedure KRdihedral_single
    procedure KRdihedral_batch
  end interface KRdihedral
  
  ! Chebyshev coefficients for phi'(u) for different k
  real(real64), parameter :: CHEB_NUM_2(15) = [ &
    0.39627987780106344_real64, 0.39338452161224874_real64, -0.003544966435596199_real64, &
    -0.0006823304169040724_real64, -3.604149307300239e-05_real64, -3.1605464484494636e-06_real64, &
    2.0432777448919986e-07_real64, 5.053742895615855e-08_real64, 7.523441384495686e-09_real64, &
    5.248732909091277e-10_real64, -2.2865932040122213e-11_real64, -1.2314903279531874e-11_real64, &
    -2.0242943434689516e-12_real64, -1.6294199721876062e-13_real64, 3.866040574806242e-15_real64]
  
  real(real64), parameter :: CHEB_NUM_3(13) = [ &
    0.26372079493664446_real64, 0.2621214219744662_real64, -0.0019209446532809178_real64, &
    -0.00032215106358461057_real64, -5.042526775859698e-07_real64, 1.1411079542434066e-07_real64, &
    4.183849390993462e-08_real64, 2.789976730130904e-09_real64, -6.884027295928439e-11_real64, &
    -1.2483951830916921e-11_real64, -1.1964284939950267e-12_real64, -2.0890865865872803e-14_real64, &
    6.360346416223058e-15_real64]
  
  real(real64), parameter :: CHEB_NUM_4(11) = [ &
    0.19740759198647356_real64, 0.1965230199482321_real64, -0.0010592628708664318_real64, &
    -0.0001736206447833958_real64, 1.2067879204577828e-06_real64, 1.413140039530892e-07_real64, &
    4.969033409020158e-09_real64, 2.3319140405071373e-10_real64, -2.3209109765292826e-11_real64, &
    -1.2843902376725979e-12_real64, 1.0043923509978181e-14_real64]
  
  real(real64), parameter :: CHEB_NUM_6(10) = [ &
    0.13130528077331868_real64, 0.13096612621072035_real64, -0.00040606488851938056_real64, &
    -6.647953265805797e-05_real64, 4.781749066832917e-07_real64, 4.723681792250029e-08_real64, &
    -1.5943494316860916e-10_real64, -1.5278450598504546e-11_real64, -6.980104484333946e-13_real64, &
    -2.7112824751072206e-14_real64]
  
  real(real64), parameter :: PHI_MAX_2 = 0.7853981633974483_real64
  real(real64), parameter :: PHI_MAX_3 = 0.5235987755982988_real64
  real(real64), parameter :: PHI_MAX_4 = 0.39269908169872414_real64
  real(real64), parameter :: PHI_MAX_6 = 0.2617993877991494_real64

contains

!--------------------------------------------------------------------------
pure function a_polar(k) result(a)
  !DEC$ ATTRIBUTES DLLEXPORT :: a_polar
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! 

  IMPLICIT NONE

  integer, INTENT(IN) :: k
  real(real64)        :: a
  
  a = 1.0_real64 / tan(cPi / (2.0_real64 * real(k, real64)))

end function a_polar

!--------------------------------------------------------------------------
pure function rho_bound(thp, php, k) result(rho)
  !DEC$ ATTRIBUTES DLLEXPORT :: a_polar
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! 

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: thp, php
  integer, INTENT(IN)       :: k
  real(real64)              :: rho
  real(real64)              :: a_val, b, c, rho3
  
  a_val = a_polar(k)
  b = cos(php)
  c = max(a_val * cos(thp), b * sin(thp))
  c = max(c, EPS)
  rho3 = 1.5_real64 * (atan(1.0_real64 / c) - c / (1.0_real64 + c * c))
  rho = max(rho3, 0.0_real64) ** (1.0_real64 / 3.0_real64)

end function rho_bound

!--------------------------------------------------------------------------
pure function F_pc(theta, a) result(F)
  !DEC$ ATTRIBUTES DLLEXPORT :: a_polar
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! 

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: theta, a
  real(real64)              :: F
  real(real64)              :: cth
  
  cth = cos(theta)
  F = 1.5_real64 * (0.5_real64 * cPi * (1.0_real64 - cth) - atan(a) + cth * atan(a * cth))

end function F_pc

!--------------------------------------------------------------------------
pure function F_lat(theta, b) result(F)
  !DEC$ ATTRIBUTES DLLEXPORT :: F_lat
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! 

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: theta, b
  real(real64)              :: F
  real(real64)              :: bs, term1, sqrt1pb2, term2
  
  bs = b * sin(theta)
  term1 = 0.5_real64 * cPi - atan(bs)
  sqrt1pb2 = sqrt(1.0_real64 + b * b)
  term2 = (b / sqrt1pb2) * atan2(sqrt1pb2 * sin(theta), cos(theta))
  F = 1.5_real64 * (0.5_real64 * cPi - cos(theta) * term1 - term2)

end function F_lat

!--------------------------------------------------------------------------
pure function theta_switch(phi_p, a) result(ths)
  !DEC$ ATTRIBUTES DLLEXPORT :: theta_switch
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! 

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: phi_p, a
  real(real64)              :: ths
  real(real64)              :: b
  
  b = cos(phi_p)
  ths = atan2(a, b)

end function theta_switch

!--------------------------------------------------------------------------
pure function A_phi(phi_p, k) result(A)
  !DEC$ ATTRIBUTES DLLEXPORT :: A_phi
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! 

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: phi_p
  integer, INTENT(IN)       :: k
  real(real64)              :: A
  real(real64)              :: a_val, ths, b
  
  a_val = a_polar(k)
  ths = theta_switch(phi_p, a_val)
  b = cos(phi_p)
  A = F_pc(ths, a_val) + (F_lat(0.5_real64 * cPi, b) - F_lat(ths, b))

end function A_phi

!--------------------------------------------------------------------------
pure function G_theta(theta, phi_p, k) result(G)
  !DEC$ ATTRIBUTES DLLEXPORT :: G_theta
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! 

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: theta, phi_p
  integer, INTENT(IN)       :: k
  real(real64)              :: G
  real(real64)              :: a_val, ths, b
  
  a_val = a_polar(k)
  ths = theta_switch(phi_p, a_val)
  b = cos(phi_p)
  
  if (theta <= ths) then
    G = F_pc(theta, a_val)
  else
    G = F_pc(ths, a_val) + (F_lat(theta, b) - F_lat(ths, b))
  end if

end function G_theta

!--------------------------------------------------------------------------
pure function f_theta(theta, phi_p, k) result(f)
  !DEC$ ATTRIBUTES DLLEXPORT :: f_theta
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! 

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: theta, phi_p
  integer, INTENT(IN)       :: k
  real(real64)              :: f
  real(real64)              :: a_val, b, c, rho3
  
  a_val = a_polar(k)
  b = cos(phi_p)
  c = max(a_val * cos(theta), b * sin(theta))
  c = max(c, EPS)
  rho3 = 1.5_real64 * (atan(1.0_real64 / c) - c / (1.0_real64 + c * c))
  f = rho3 * sin(theta)

end function f_theta

!--------------------------------------------------------------------------
pure function phi_from_u_pade(u, k) result(phi)
  !DEC$ ATTRIBUTES DLLEXPORT :: phi_from_u_pade
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! 

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: u
  integer, INTENT(IN)       :: k
  real(real64)              :: phi
  real(real64)              :: t
  real(real64)              :: coeffs(size(transfer(k, [1.0_real64])))
  real(real64), parameter   :: PHI_MAX_VALUES(4) = [PHI_MAX_2, PHI_MAX_3, PHI_MAX_4, PHI_MAX_6]
  
  t = max(-1.0_real64, min(1.0_real64, 2.0_real64 * u - 1.0_real64))
  
  select case (k)
  case (2)
    phi = chebval_T(t, CHEB_NUM_2)
  case (3)
    phi = chebval_T(t, CHEB_NUM_3)
  case (4)
    phi = chebval_T(t, CHEB_NUM_4)
  case (6)
    phi = chebval_T(t, CHEB_NUM_6)
  case default
    phi = 0.0_real64
  end select
  
  select case (k)
  case (2)
    phi = max(0.0_real64, min(PHI_MAX_2, phi))
  case (3)
    phi = max(0.0_real64, min(PHI_MAX_3, phi))
  case (4)
    phi = max(0.0_real64, min(PHI_MAX_4, phi))
  case (6)
    phi = max(0.0_real64, min(PHI_MAX_6, phi))
  end select

end function phi_from_u_pade

!--------------------------------------------------------------------------
pure function invert_theta_newton(y, phi_p, k, max_iter) result(theta)
  !DEC$ ATTRIBUTES DLLEXPORT :: invert_theta_newton
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! 

  IMPLICIT NONE

  real(real64), INTENT(IN)      :: y, phi_p
  integer, INTENT(IN)           :: k
  integer, INTENT(IN), OPTIONAL :: max_iter
  real(real64)                  :: theta
  integer                       :: iter, niter
  real(real64)                  :: A, Gn, fn, step
  
  niter = 11
  if (present(max_iter)) niter = max_iter
  
  A = max(A_phi(phi_p, k), EPS)
  theta = acos(max(-1.0_real64, min(1.0_real64, 1.0_real64 - y)))
  
  do iter = 1, niter
    Gn = G_theta(theta, phi_p, k) / A
    fn = max(f_theta(theta, phi_p, k) / A, EPS_EPS)
    step = (Gn - y) / fn
    theta = max(0.0_real64, min(HALF_PI, theta - step))
  end do

end function invert_theta_newton

!--------------------------------------------------------------------------
pure function rotate_xy(v, ang) result(vr)
  !DEC$ ATTRIBUTES DLLEXPORT :: rotate_xy
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! 

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: v(3), ang
  real(real64)              :: vr(3)
  real(real64)              :: ca, sa
  
  ca = cos(ang)
  sa = sin(ang)
  vr(1) = v(1) * ca - v(2) * sa
  vr(2) = v(1) * sa + v(2) * ca
  vr(3) = v(3)

end function rotate_xy

!--------------------------------------------------------------------------
pure function reflect_about_axis(v, alpha, mask) result(vr)
  !DEC$ ATTRIBUTES DLLEXPORT :: reflect_about_axis
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! 

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: v(3), alpha
  logical, INTENT(IN)       :: mask
  real(real64)              :: vr(3)
  real(real64)              :: ca2, sa2
  
  vr = v
  if (mask) then
    ca2 = cos(2.0_real64 * alpha)
    sa2 = sin(2.0_real64 * alpha)
    vr(1) = ca2 * v(1) + sa2 * v(2)
    vr(2) = sa2 * v(1) - ca2 * v(2)
  end if

end function reflect_about_axis

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Public API functions
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
subroutine KRdihedral_single(ho, mapped, k, eps)
  !DEC$ ATTRIBUTES DLLEXPORT :: KRdihedral_single
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! Core mapper

  IMPLICIT NONE

  real(real64), INTENT(IN)            :: ho(3)
  integer, INTENT(IN)                 :: k
  real(real64), INTENT(OUT)           :: mapped(3)
  real(real64), INTENT(IN), OPTIONAL  :: eps
  
  real(real64)                        :: v(3), v_temp(3)
  real(real64)                        :: sgn_z, sector_half, wedge_max
  real(real64)                        :: sph(3), phi0, krot, phi1, phi2
  logical                             :: mask_neg, mask_over
  real(real64)                        :: r, th, ph, u, ph_p, y, th_p, rho, r_p
  real(real64)                        :: eps_val
  
  eps_val = 1.0e-8_real64
  if (present(eps)) eps_val = eps
  
  v = ho
  
  ! Fold to theta in [0, cPi/2]
  if (v(3) >= 0.0_real64) then
    sgn_z = 1.0_real64
  else
    sgn_z = -1.0_real64
  end if
  v(3) = abs(v(3))
  
  sector_half = cPi / real(k, real64)
  wedge_max = cPi / (2.0_real64 * real(k, real64))
  
  ! Rotate to sector
  sph = sph_from_cart(v)
  phi0 = sph(3)
  krot = floor((phi0 + sector_half) / (2.0_real64 * sector_half))
  v = rotate_xy(v, -krot * (2.0_real64 * sector_half))
  
  ! Reflect into [0, wedge_max]
  sph = sph_from_cart(v)
  phi1 = sph(3)
  mask_neg = phi1 < 0.0_real64
  if (mask_neg) then
    v = reflect_about_axis(v, 0.0_real64, .true.)
  end if
  
  sph = sph_from_cart(v)
  phi2 = sph(3)
  mask_over = phi2 > wedge_max
  if (mask_over) then
    v = reflect_about_axis(v, wedge_max, .true.)
  end if
  
  ! KR inside wedge
  sph = sph_from_cart(v)
  r = sph(1)
  th = sph(2)
  ph = sph(3)
  
  u = max(0.0_real64, min(1.0_real64, ph / wedge_max))
  ph_p = phi_from_u_pade(u, k)
  y = 1.0_real64 - cos(th)
  th_p = invert_theta_newton(y, ph_p, k)
  rho = rho_bound(th_p, ph_p, k)
  r_p = (1.0_real64 - eps_val) * rho * (r / H_MAX)
  mapped = cart_from_sph(r_p, th_p, ph_p)
  
  ! Undo reflections/rotations, restore sign
  if (mask_over) then
    mapped = reflect_about_axis(mapped, wedge_max, .true.)
  end if
  if (mask_neg) then
    mapped = reflect_about_axis(mapped, 0.0_real64, .true.)
  end if
  mapped = rotate_xy(mapped, krot * (2.0_real64 * sector_half))
  mapped(3) = sign(mapped(3), sgn_z)

end subroutine KRdihedral_single

!--------------------------------------------------------------------------
subroutine KRdihedral_batch(h_in, h_out, n, k)
  !DEC$ ATTRIBUTES DLLEXPORT :: KRdihedral_batch
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! batch conversion with OpenMP

IMPLICIT NONE

  real(real64), INTENT(IN)      :: h_in(3,n)
  real(real64), INTENT(OUT)     :: h_out(3,n)
  integer(kind=irg), INTENT(IN) :: n, k

  integer(kind=irg)             :: i
  
!$omp parallel do default(none) shared(h_in, h_out, n, k)
  do i = 1, n
    call KRdihedral_single(h_in(:,i), h_out(:,i), k)
  end do
!$omp end parallel do

end subroutine KRdihedral_batch

end module mod_KRdihedral

