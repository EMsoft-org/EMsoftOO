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

module mod_KRsupport
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! Support routines for KR sampling modules (to avoid code repetition)
  !! 
  use mod_kinds
  use mod_global

  use, intrinsic :: iso_fortran_env, only: real64

  IMPLICIT NONE
  private

  public :: chebval_T, sph_from_cart, cart_from_sph, rho3_from_c, R_from_c, sph_from_cart_vec
  
  ! Constants
  real(real64), parameter, public :: BALL_RADIUS = 1.3306700394914688_real64
  real(real64), parameter, public :: EPS = 1.0e-15_real64
  real(real64), parameter, public :: EPS_EPS = 1.0e-14_real64
  real(real64), parameter, public :: EPS_C = 1.0e-15_real64
  real(real64), parameter, public :: EPS_SLOPE = 1.0e-15_real64
  real(real64), parameter, public :: PHI_LO_O = 0.25_real64 * cPi
  real(real64), parameter, public :: PHI_HI_O = 0.50_real64 * cPi
  real(real64), parameter, public :: KAPPA = sqrt(2.0_real64) - 1.0_real64
  real(real64), parameter, public :: H_MAX = (3.0_real64 * cPi / 4.0_real64) ** (1.0_real64 / 3.0_real64)
  real(real64), parameter, public :: PHI_MAX_T = 0.25_real64 * cPi
  real(real64), parameter, public :: HALF_PI = 0.5_real64 * cPi

  ! Gauss-Legendre 16-point quadrature nodes
  real(real64), parameter, public :: GL16_X(16) = [ &
    -0.9894009349916499_real64, -0.9445750230732326_real64, -0.8656312023878317_real64, &
    -0.7554044083550031_real64, -0.6178762444026437_real64, -0.4580167776572274_real64, &
    -0.2816035507792589_real64, -0.09501250983763744_real64, 0.09501250983763744_real64, &
    0.2816035507792589_real64, 0.4580167776572274_real64, 0.6178762444026437_real64, &
    0.7554044083550031_real64, 0.8656312023878317_real64, 0.9445750230732326_real64, &
    0.9894009349916499_real64]
  
  ! Gauss-Legendre 16-point quadrature weights
  real(real64), parameter, public :: GL16_W(16) = [ &
    0.027152459411754095_real64, 0.06225352393864789_real64, 0.09515851168249278_real64, &
    0.12462897125553387_real64, 0.14959598881657673_real64, 0.16915651939500254_real64, &
    0.1826034150449236_real64, 0.1894506104550685_real64, 0.1894506104550685_real64, &
    0.1826034150449236_real64, 0.16915651939500254_real64, 0.14959598881657673_real64, &
    0.12462897125553387_real64, 0.09515851168249278_real64, 0.06225352393864789_real64, &
    0.027152459411754095_real64]

contains

!--------------------------------------------------------------------------
pure function chebval_T(t, coeffs) result(val)
!DEC$ ATTRIBUTES DLLEXPORT :: chebval_T
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !1 Evaluate Chebyshev polynomial sum_k coeffs(k) * T_k(t) using Clenshaw algorithm
  !! constructor for the SpaceGroup Class

IMPLICIT NONE

  real(real64), INTENT(IN)  :: t
  real(real64), INTENT(IN)  :: coeffs(:)
  real(real64)              :: val
  
  integer(kind=irg)         :: k, n
  real(real64)              :: b_k, b_kp1, b_kp2
  
  n = size(coeffs)
  if (n == 1) then
    val = coeffs(1)
    return
  end if
  
  b_kp2 = 0.0_real64
  b_kp1 = 0.0_real64
  
  do k = n, 2, -1
    b_k = 2.0_real64 * t * b_kp1 - b_kp2 + coeffs(k)
    b_kp2 = b_kp1
    b_kp1 = b_k
  end do
  
  val = coeffs(1) + t * b_kp1 - b_kp2
end function chebval_T

!--------------------------------------------------------------------------
pure function sph_from_cart(v) result(sph)
!DEC$ ATTRIBUTES DLLEXPORT :: sph_from_cart
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! Spherical to cartesian conversion

IMPLICIT NONE  

  real(real64), INTENT(IN)  :: v(3)
  real(real64)              :: sph(3)  ! [r, theta, phi]
  real(real64)              :: r, ct
  
  r = norm2(v)
  if (r > 0.0_real64) then
    ct = max(-1.0_real64, min(1.0_real64, v(3) / max(r, EPS)))
    sph(2) = acos(ct)
  else
    sph(2) = 0.0_real64
  end if
  sph(1) = r
  sph(3) = atan2(v(2), v(1))

end function sph_from_cart

!--------------------------------------------------------------------------
pure function cart_from_sph(r, th, ph) result(cart)
!DEC$ ATTRIBUTES DLLEXPORT :: cart_from_sph
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! Cartesian to spherical conversion

IMPLICIT NONE

  real(real64), INTENT(IN)  :: r, th, ph
  real(real64)              :: cart(3)
  real(real64)              :: st, ct, cp, sp
  
  st = sin(th)
  ct = cos(th)
  cp = cos(ph)
  sp = sin(ph)
  cart(1) = r * st * cp
  cart(2) = r * st * sp
  cart(3) = r * ct

end function cart_from_sph

!--------------------------------------------------------------------------
pure function rho3_from_c(c) result(rho3)
!DEC$ ATTRIBUTES DLLEXPORT :: rho3_from_c
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! Homochoric radial law

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: c
  real(real64)              :: rho3
  real(real64)              :: c_safe
  
  c_safe = max(c, EPS)
  rho3 = 1.5_real64 * (atan(1.0_real64 / c_safe) - c_safe / (1.0_real64 + c_safe * c_safe))

end function rho3_from_c

!--------------------------------------------------------------------------
pure function R_from_c(c) result(R)
!DEC$ ATTRIBUTES DLLEXPORT :: rho3_from_c
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    

  IMPLICIT NONE

  real(real64), intent(in)  :: c
  real(real64)              :: R
  
  R = max(rho3_from_c(c), 0.0_real64) ** (1.0_real64 / 3.0_real64)

end function R_from_c

!--------------------------------------------------------------------------
pure subroutine sph_from_cart_vec(v, r, theta, phi)
!DEC$ ATTRIBUTES DLLEXPORT :: sph_from_cart_vec
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    

  IMPLICIT NONE

  real(real64), intent(in)  :: v(3)
  real(real64), intent(out) :: r, theta, phi
  real(real64)              :: ct
  
  r = norm2(v)
  if (r > 0.0_real64) then
    ct = max(-1.0_real64, min(1.0_real64, v(3) / r))
    theta = acos(ct)
  else
    theta = 0.0_real64
  end if
  phi = atan2(v(2), v(1))

end subroutine sph_from_cart_vec

end module mod_KRsupport

