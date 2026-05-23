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

  module mod_KRoctahedral
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
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

  public :: KRoctahedral

  interface KRoctahedral
    procedure KRoctahedral_single
    procedure KRoctahedral_batch
  end interface KRoctahedral
  
  
  ! Padé numerator coefficients (n=15)
  real(real64), parameter :: A_COEFFS(16) = [ &
    1.1822546284090303e+00_real64, 6.2667883943387037e-01_real64, 6.0044992184703982e-01_real64, &
    3.8819631496309054e-01_real64, -4.4145326792933237e-01_real64, -3.0362113191584633e-01_real64, &
    1.8863124699476139e-01_real64, 7.6897095197285464e-02_real64, 9.8597202925534602e-02_real64, &
    2.7635698327477154e-01_real64, -2.0057747872002395e-02_real64, 6.1319952952553401e-02_real64, &
    -1.3776491479530030e-03_real64, 8.2840323424582154e-04_real64, 2.1685353523857393e-04_real64, &
    1.4705470527970237e-05_real64]
  
  ! Padé denominator coefficients (m=13, b[0]=1 fixed)
  real(real64), parameter :: B_COEFFS(14) = [ &
    1.0000000000000000e+00_real64, 1.2427500771143379e-01_real64, 4.2410008089771739e-01_real64, &
    3.3639100888253853e-01_real64, -4.0743628876474214e-01_real64, -2.3247960338991300e-01_real64, &
    2.0383650417779056e-01_real64, 2.5279219590521870e-02_real64, 3.7732606982845002e-02_real64, &
    2.4437987986578913e-01_real64, -7.0431650497779760e-02_real64, 6.4426522575716655e-02_real64, &
    -1.1610440319712067e-02_real64, 2.0125442894998288e-03_real64]
  
contains

!--------------------------------------------------------------------------
pure function theta_max_O(phi) result(th_max)
!DEC$ ATTRIBUTES DLLEXPORT :: theta_max_O
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!

  IMPLICIT NONE 

  real(real64), INTENT(IN)  :: phi
  real(real64)              :: th_max
  real(real64)              :: s
  
  s = max(sin(phi), EPS_C)
  th_max = atan(1.0_real64 / s)

end function theta_max_O

!--------------------------------------------------------------------------
pure function theta_switch_O(phi) result(ths)
!DEC$ ATTRIBUTES DLLEXPORT :: theta_switch_O
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!

IMPLICIT NONE 

  real(real64), INTENT(IN)   :: phi
  real(real64)              :: ths
  real(real64)              :: denom
  
  denom = max(cos(phi) + sin(phi), EPS_C)
  ths = atan(sqrt(2.0_real64) / denom)

end function theta_switch_O

!--------------------------------------------------------------------------
pure subroutine c_O_components(theta, phi, c_max, c_sum)
!DEC$ ATTRIBUTES DLLEXPORT :: c_O_components
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!

  IMPLICIT NONE 
  
  real(real64), INTENT(IN)  :: theta, phi
  real(real64), INTENT(OUT) :: c_max, c_sum
  
  c_max = (1.0_real64 / KAPPA) * cos(theta)
  c_sum = cos(theta) + sin(theta) * (cos(phi) + sin(phi))

end subroutine c_O_components

!--------------------------------------------------------------------------
 function A_phi_O(phi) result(A)
!DEC$ ATTRIBUTES DLLEXPORT :: A_phi_O
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! Column area A(φ) using GL16 integration

  IMPLICIT NONE 
  
  real(real64), INTENT(IN)  :: phi
  real(real64)              :: A
  real(real64)              :: ths, thh
  integer                   :: i
  real(real64)              :: xa, wa, sum_pc, sum_sum
  real(real64)              :: th_eval, f_val
  
  ths = theta_switch_O(phi)
  thh = theta_max_O(phi)
  
  if (ths >= thh) then
    ! Only polar cap branch
    sum_pc = 0.0_real64
    do i = 1, 16
      xa = 0.5_real64 * (GL16_X(i) * thh + thh)
      wa = 0.5_real64 * thh * GL16_W(i)
      th_eval = xa
      f_val = rho3_from_c((1.0_real64 / KAPPA) * cos(th_eval)) * sin(th_eval)
      sum_pc = sum_pc + wa * f_val
    end do
    A = sum_pc
  else
    ! Polar cap part
    sum_pc = 0.0_real64
    do i = 1, 16
      xa = 0.5_real64 * (GL16_X(i) * ths + ths)
      wa = 0.5_real64 * ths * GL16_W(i)
      th_eval = xa
      f_val = rho3_from_c((1.0_real64 / KAPPA) * cos(th_eval)) * sin(th_eval)
      sum_pc = sum_pc + wa * f_val
    end do
    
    ! Sum branch part
    sum_sum = 0.0_real64
    do i = 1, 16
      xa = 0.5_real64 * (GL16_X(i) * (thh - ths) + (ths + thh))
      wa = 0.5_real64 * (thh - ths) * GL16_W(i)
      th_eval = xa
      f_val = rho3_from_c(cos(th_eval) + sin(th_eval) * (cos(phi) + sin(phi))) * sin(th_eval)
      sum_sum = sum_sum + wa * f_val
    end do
    
    A = sum_pc + sum_sum
  end if
  
  A = max(A, EPS_SLOPE)

end function A_phi_O

!--------------------------------------------------------------------------
function C_src_phi_O(phi, phi_lo, phi_hi) result(C)
 !DEC$ ATTRIBUTES DLLEXPORT :: C_src_phi_O
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! Source azimuthal CDF

  IMPLICIT NONE 

  real(real64), INTENT(IN)  :: phi, phi_lo, phi_hi
  real(real64)              :: C
  integer                   :: i
  real(real64)              :: xa, wa, psi, thh, sum_val
  
  sum_val = 0.0_real64
  do i = 1, 16
    xa = 0.5_real64 * (GL16_X(i) * (phi - phi_lo) + (phi_lo + phi))
    wa = 0.5_real64 * (phi - phi_lo) * GL16_W(i)
    psi = xa
    thh = theta_max_O(psi)
    sum_val = sum_val + wa * (1.0_real64 - cos(thh))
  end do
  C = sum_val

end function C_src_phi_O

!--------------------------------------------------------------------------
pure function phi_inv_O(u) result(phi)
!DEC$ ATTRIBUTES DLLEXPORT :: phi_inv_O
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! Padé inverse CDF for O

  IMPLICIT NONE 
  
  real(real64), INTENT(IN)  :: u
  real(real64)              :: phi
  real(real64)              :: t, P, Q
  real(real64)              :: d_coeffs(size(B_COEFFS))
  
  t = max(-1.0_real64, min(1.0_real64, 2.0_real64 * u - 1.0_real64))
  P = chebval_T(t, A_COEFFS)
  d_coeffs = B_COEFFS
  d_coeffs(1) = 0.0_real64
  Q = 1.0_real64 + chebval_T(t, d_coeffs)
  phi = max(PHI_LO_O, min(PHI_HI_O, P / Q))

end function phi_inv_O

!--------------------------------------------------------------------------
 pure subroutine swap_int(a, b)
!DEC$ ATTRIBUTES DLLEXPORT :: swap_int
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!

  IMPLICIT NONE 
  
  integer, INTENT(INOUT)  :: a, b
  integer                 :: tmp

  tmp = a
  a = b
  b = tmp

end subroutine swap_int

!--------------------------------------------------------------------------
function gl16_integrate_0_to_b(b, phi, mode) result(val)
!DEC$ ATTRIBUTES DLLEXPORT :: gl16_integrate_0_to_b
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! GL16 integration helper for [0, b]

  IMPLICIT NONE 
  
  real(real64), INTENT(IN)  :: b, phi
  integer, INTENT(IN)       :: mode  ! 0: polar cap, 1: sum branch
  real(real64)              :: val
  integer                   :: i
  real(real64)              :: xa, wa, th_eval, f_val
  real(real64)              :: trig_phi
  
  val = 0.0_real64
  trig_phi = cos(phi) + sin(phi)
  
  do i = 1, 16
    xa = 0.5_real64 * (GL16_X(i) * b + b)
    wa = 0.5_real64 * b * GL16_W(i)
    th_eval = xa
    
    if (mode == 0) then
      f_val = rho3_from_c((1.0_real64 / KAPPA) * cos(th_eval)) * sin(th_eval)
    else
      f_val = rho3_from_c(cos(th_eval) + sin(th_eval) * trig_phi) * sin(th_eval)
    end if
    
    val = val + wa * f_val
  end do

end function gl16_integrate_0_to_b

!--------------------------------------------------------------------------
function gl16_integrate_a_to_b(a, b, phi, mode) result(val)
!DEC$ ATTRIBUTES DLLEXPORT :: gl16_integrate_a_to_b
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! GL16 integration helper for [a, b]

  IMPLICIT NONE 
  
  real(real64), INTENT(IN)  :: a, b, phi
  integer, INTENT(IN)       :: mode  ! 0: polar cap, 1: sum branch
  real(real64)              :: val
  integer                   :: i
  real(real64)              :: xa, wa, th_eval, f_val
  real(real64)              :: trig_phi, half, mid
  
  val = 0.0_real64
  trig_phi = cos(phi) + sin(phi)
  half = 0.5_real64 * (b - a)
  mid = 0.5_real64 * (b + a)
  
  do i = 1, 16
    xa = half * GL16_X(i) + mid
    wa = half * GL16_W(i)
    th_eval = xa
    
    if (mode == 0) then
      f_val = rho3_from_c((1.0_real64 / KAPPA) * cos(th_eval)) * sin(th_eval)
    else
      f_val = rho3_from_c(cos(th_eval) + sin(th_eval) * trig_phi) * sin(th_eval)
    end if
    
    val = val + wa * f_val
  end do

end function gl16_integrate_a_to_b

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Public API functions
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
subroutine KRoctahedral_single(ho, mapped, newton_iters)
!DEC$ ATTRIBUTES DLLEXPORT :: KRoctahedral_single
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! Main mapping function

  IMPLICIT NONE 
  
  real(real64), INTENT(IN)      :: ho(3)
  real(real64), INTENT(OUT)     :: mapped(3)
  integer, INTENT(IN), OPTIONAL :: newton_iters
  
  real(real64)                  :: v0(3), v(3), v_sorted(3)
  real(real64)                  :: sgn(3)
  integer                       :: idxs(3), inv_idx(3)
  integer                       :: i, niter
  real(real64)                  :: r, th, ph
  real(real64)                  :: C_hi, u_src, ph_p
  real(real64)                  :: th_hi_src, denom, y_src
  real(real64)                  :: A_tar, th_hi_tar, th_p
  real(real64)                  :: ths, only_pc_mask
  real(real64)                  :: G, Gn, fn, step, th_new
  real(real64)                  :: c_max, c_sum, c_act, rho, r_p
  real(real64)                  :: big_mask
  integer                       :: j
  
  ! Initialize
  niter = 8
  if (present(newton_iters)) niter = newton_iters
  
  v0 = ho
  
  ! Fold to first octant and enforce x ≤ y ≤ z (ordered simplex)
  sgn = sign(1.0_real64, v0)
  where (sgn == 0.0_real64) sgn = 1.0_real64
  v = abs(v0)
  
  ! Sort to ordered simplex
  idxs = [1, 2, 3]
  do i = 1, 2
    do j = i + 1, 3
      if (v(idxs(j)) < v(idxs(i))) then
        call swap_int(idxs(i), idxs(j))
      end if
    end do
  end do
  v_sorted = [v(idxs(1)), v(idxs(2)), v(idxs(3))]
  
  ! Build inverse permutation
  inv_idx(idxs(1)) = 1
  inv_idx(idxs(2)) = 2
  inv_idx(idxs(3)) = 3
  
  ! Spherical on source sector
  call sph_from_cart_vec(v_sorted, r, th, ph)
  
  ! Outer KR: source CDF then Padé inverse to φ'
  C_hi = C_src_phi_O(PHI_HI_O, PHI_LO_O, PHI_HI_O)
  u_src = max(0.0_real64, min(1.0_real64, C_src_phi_O(ph, PHI_LO_O, PHI_HI_O) / C_hi))
  ph_p = phi_inv_O(u_src)
  
  ! Inner KR: y_src normalized by source ceiling
  th_hi_src = theta_max_O(ph)
  denom = max(1.0_real64 - cos(th_hi_src), 1.0e-12_real64)
  y_src = max(0.0_real64, min(1.0_real64, (1.0_real64 - cos(th)) / denom))
  
  ! Target column area A_tar(φ')
  A_tar = A_phi_O(ph_p)
  
  ! Initial guess
  th_hi_tar = theta_max_O(ph_p)
  th_p = acos(max(-1.0_real64, min(1.0_real64, 1.0_real64 - y_src * (1.0_real64 - cos(th_hi_tar)))))
  
  ! Safeguarded Newton
  do i = 1, niter
    th_p = max(0.0_real64, min(th_p, th_hi_tar))
    
    ths = theta_switch_O(ph_p)
    only_pc_mask = merge(1.0_real64, 0.0_real64, ths >= th_hi_tar)
    
    if (only_pc_mask > 0.5_real64) then
      ! Entirely polar-cap
      G = gl16_integrate_0_to_b(th_p, ph_p, 0)
    else
      ! Column that switches
      if (th_p <= ths) then
        G = gl16_integrate_0_to_b(th_p, ph_p, 0)
      else
        G = gl16_integrate_0_to_b(ths, ph_p, 0) + &
            gl16_integrate_a_to_b(ths, th_p, ph_p, 1)
      end if
    end if
    
    Gn = max(0.0_real64, min(1.0_real64, G / A_tar))
    
    ! Slope f_tar/A_tar at θ'
    call c_O_components(th_p, ph_p, c_max, c_sum)
    c_act = max(c_max, c_sum)
    c_act = max(c_act, EPS_C)
    fn = max(rho3_from_c(c_act) * sin(th_p) / A_tar, EPS_SLOPE)
    
    step = (Gn - y_src) / fn
    th_new = th_p - step
    big_mask = merge(1.0_real64, 0.0_real64, abs(step) > 0.25_real64)
    if (big_mask > 0.5_real64) then
      th_new = 0.75_real64 * th_p + 0.25_real64 * th_new
    end if
    th_new = max(0.0_real64, min(th_new, th_hi_tar))
    
    if (abs(th_new - th_p) < 1.0e-12_real64) then
      th_p = th_new
      exit
    end if
    th_p = th_new
  end do
  
  ! Radial scaling
  call c_O_components(th_p, ph_p, c_max, c_sum)
  rho = R_from_c(max(c_max, c_sum))
  r_p = rho * (r / H_MAX)
  mapped = cart_from_sph(r_p, th_p, ph_p)
  
  ! Undo permutation and signs
  mapped = [mapped(inv_idx(1)), mapped(inv_idx(2)), mapped(inv_idx(3))]
  mapped = mapped * sgn

end subroutine KRoctahedral_single
 
!--------------------------------------------------------------------------
subroutine KRoctahedral_batch(h_in, h_out, n, newton_iters)
!DEC$ ATTRIBUTES DLLEXPORT :: KRoctahedral_batch
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
      call KRoctahedral_single(h_in(:, i), h_out(:, i), newton_iters)
    else
      call KRoctahedral_single(h_in(:, i), h_out(:, i))
    end if
  end do
!$omp end parallel do

end subroutine KRoctahedral_batch

end module mod_KRoctahedral

