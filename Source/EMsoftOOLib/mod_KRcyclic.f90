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

module mod_KRcyclic
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!
  !! Homochoric to cyclic FZ mapping module
  !! Supports C2, C3, C4, and C6 symmetries
  !! 
  !! This module provides vectorized mapping functions for batch processing.
  !! For OpenMP parallelization, use the batch functions with appropriate
  !! compiler flags: -fopenmp (GCC) or -openmp (Intel)
  
  use mod_kinds
  use mod_global
  use mod_KRsupport
  use omp_lib

  use, intrinsic :: iso_fortran_env, only: real64

  IMPLICIT NONE

  private

  public :: KRcyclic

  interface KRcyclic
    procedure KRcyclic_single
    procedure KRcyclic_batch
  end interface KRcyclic

  ! Coefficients for C2
  real(real64), parameter :: tau_C2 = 1.0_real64
  real(real64), parameter :: num_C2(11) = [ &
    1.7737848743975086e00_real64, 4.8055082084235462e-01_real64, -4.2901615601460463e-01_real64, &
    3.5067723530850670e-01_real64, -2.5078415551551952e-01_real64, 1.4245817838614297e-02_real64, &
    2.0364985682640563e-01_real64, 7.3494423671352616e-02_real64, 7.4055714754984337e-03_real64, &
    -2.6791200292731122e-06_real64, -1.3705227262528180e-05_real64]
  real(real64), parameter :: den_C2(11) = [ &
    1.0000000000000000e00_real64, 4.1494079856907912e-01_real64, -2.1735212942622864e-01_real64, &
    1.5858436443785079e-01_real64, -1.1731723315019479e-01_real64, 1.6329335772084991e-03_real64, &
    1.1851852745686715e-01_real64, 5.0425562752368529e-02_real64, 6.3676661436295858e-03_real64, &
    5.5532866871394670e-05_real64, -1.8803246198828434e-05_real64]
  
  ! Coefficients for C3
  real(real64), parameter :: tau_C3 = 1.0_real64 / sqrt(3.0_real64)
  real(real64), parameter :: num_C3(11) = [ &
    1.8258508913178009e00_real64, 1.0204210503954347e00_real64, -6.8990113854728552e-01_real64, &
    2.8365299613149925e-01_real64, 1.7567952865394379e-01_real64, -3.1873401624280262e-01_real64, &
    -2.2715481819120306e-01_real64, -5.2585665924886747e-02_real64, -4.3354111038882591e-03_real64, &
    -3.5816570541181354e-05_real64, 3.6250774202656139e-06_real64]
  real(real64), parameter :: den_C3(11) = [ &
    1.0000000000000000e00_real64, 7.2739743922667077e-01_real64, -2.8657823943183158e-01_real64, &
    1.0022667454700472e-01_real64, 9.6960768233829570e-02_real64, -1.6807702882407449e-01_real64, &
    -1.4360331365310397e-01_real64, -4.0424831299393671e-02_real64, -4.3922845051649244e-03_real64, &
    -8.8411935275716294e-05_real64, 6.4825070887252555e-06_real64]
  
  ! Coefficients for C4
  real(real64), parameter :: tau_C4 = sqrt(2.0_real64) - 1.0_real64
  real(real64), parameter :: num_C4(21) = [ &
    2.0638830467685660e00_real64, 2.6021401239522907e-01_real64, -1.2186320073474909e-01_real64, &
    5.5221717032306668e-02_real64, -4.2520303912717834e-02_real64, 3.0236893037655817e-02_real64, &
    -2.9012807237718201e-02_real64, 1.6415416344340801e-02_real64, -4.3515734269519143e-02_real64, &
    2.2695699364858511e-02_real64, -2.2550554009399114e-02_real64, 1.0886164310905866e-03_real64, &
    1.7816049007855050e-02_real64, -3.0923682555146758e-04_real64, -2.8394264832619254e-02_real64, &
    -5.3085979265811899e-02_real64, -3.4727668250554601e-02_real64, -9.1350981843724040e-03_real64, &
    -8.3214638244406738e-04_real64, 9.5344026331649555e-07_real64, 1.2203423886000405e-06_real64]
  real(real64), parameter :: den_C4(21) = [ &
    1.0000000000000000e00_real64, 5.0285292503462342e-01_real64, -9.7682933036777386e-02_real64, &
    2.6524209648882916e-02_real64, -1.5967794739660854e-02_real64, 8.2740600918967833e-03_real64, &
    -8.4219496026539543e-03_real64, -3.3279979553858394e-05_real64, -1.5859151539722664e-02_real64, &
    3.9761685229014994e-03_real64, -7.7305903613928651e-03_real64, -6.8929659704462581e-04_real64, &
    9.3824990360368652e-03_real64, -6.2969309144329875e-04_real64, -1.8185327299804658e-02_real64, &
    -3.0889735741043020e-02_real64, -2.1869520857486984e-02_real64, -6.9416967896511351e-03_real64, &
    -8.8759920534314980e-04_real64, -1.8630446847254365e-05_real64, 2.0084648627427544e-06_real64]
  
  ! Coefficients for C6
  real(real64), parameter :: tau_C6 = 2.0_real64 - sqrt(3.0_real64)
  real(real64), parameter :: num_C6(21) = [ &
    2.0791971633618513e00_real64, 5.4572868416551912e-01_real64, -3.3878248222116070e-01_real64, &
    2.5194939508555830e-01_real64, -2.1640338620750996e-01_real64, 1.9307970107561104e-01_real64, &
    -1.6750994994148374e-01_real64, 1.2930791464816124e-01_real64, -7.2433309310857660e-02_real64, &
    7.9094380324021239e-03_real64, 6.2273851101229023e-02_real64, -1.0405187160310844e-01_real64, &
    9.8843104227216927e-02_real64, -2.1019188535004646e-02_real64, -8.4291747319314506e-02_real64, &
    7.9501326731339361e-02_real64, 8.4591828839222991e-02_real64, 1.7150398056128016e-02_real64, &
    -1.3735157685798671e-03_real64, -5.2349932177729337e-04_real64, -2.0176275817922699e-05_real64]
  real(real64), parameter :: den_C6(21) = [ &
    1.0000000000000000e00_real64, 6.4478059179981684e-01_real64, -1.3942100975789928e-01_real64, &
    6.1704157630504450e-02_real64, -4.3641053539025823e-02_real64, 3.6876921485187543e-02_real64, &
    -3.2528514191822989e-02_real64, 2.6018436502876999e-02_real64, -1.3762001764068915e-02_real64, &
    2.8421427339296650e-04_real64, 1.8762998314057977e-02_real64, -3.0576085328707852e-02_real64, &
    3.1881795401240434e-02_real64, -6.6283972799896819e-03_real64, -3.5726723266283830e-02_real64, &
    3.6686209833398514e-02_real64, 4.9742917859944558e-02_real64, 1.4894263218343567e-02_real64, &
    1.2815175100273083e-04_real64, -4.3855328224429044e-04_real64, -3.5502169839994451e-05_real64]

contains

!--------------------------------------------------------------------------
pure function omega_max_from_cos(c, tau) result(omega)
!DEC$ ATTRIBUTES DLLEXPORT :: omega_max_from_cos
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: c, tau
  real(real64)              :: omega
  real(real64)              :: c_safe
  
  c_safe = max(c, EPS)
  omega = 2.0_real64 * atan(tau / c_safe)

end function omega_max_from_cos

!--------------------------------------------------------------------------
pure function sin_omega_from_cos_theta(c, tau) result(sin_omega)
!DEC$ ATTRIBUTES DLLEXPORT :: sin_omega_from_cos_theta
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: c, tau
  real(real64)              :: sin_omega
  
  sin_omega = (2.0_real64 * tau * c) / (c * c + tau * tau)

end function sin_omega_from_cos_theta

!--------------------------------------------------------------------------
pure function R_of_theta(theta, tau) result(R)
!DEC$ ATTRIBUTES DLLEXPORT :: R_of_theta
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    

  IMPLICIT NONE

  real(real64), intent(in)  :: theta, tau
  real(real64)              :: R
  real(real64)              :: c, omega, sin_omega, R3
  
  c = cos(theta)
  omega = omega_max_from_cos(c, tau)
  sin_omega = sin_omega_from_cos_theta(c, tau)
  R3 = 0.75_real64 * (omega - sin_omega)
  R = max(R3, 0.0_real64) ** (1.0_real64 / 3.0_real64)

end function R_of_theta

!--------------------------------------------------------------------------
pure function theta_fz_from_rational(theta, num, den) result(theta_fz)
 !DEC$ ATTRIBUTES DLLEXPORT :: theta_fz_from_rational
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: theta
  real(real64), INTENT(IN)  :: num(:), den(:)
  real(real64)              :: theta_fz
  real(real64)              :: y, t, num_val, den_val, g, y_clamped, sqrt_y
  real(real64)              :: den_coeffs(size(den))
  integer(kind=irg)         :: i
  
  y = 1.0_real64 - cos(theta)
  t = 2.0_real64 * y - 1.0_real64
  
  num_val = chebval_T(t, num)
  
  ! den encodes Q(t) = 1 + sum_{k>=1} b_k T_k(t)
  den_coeffs = 0.0_real64
  if (size(den) > 1) then
    do i = 2, size(den)
      den_coeffs(i) = den(i)
    end do
  end if
  den_val = 1.0_real64 + chebval_T(t, den_coeffs)
  
  g = num_val / den_val
  
  y_clamped = max(y, EPS)
  sqrt_y = sqrt(y_clamped)
  theta_fz = sqrt_y * g

end function theta_fz_from_rational

!--------------------------------------------------------------------------
subroutine map_ho_to_fz(h, tau, num, den, out)
 !DEC$ ATTRIBUTES DLLEXPORT :: map_ho_to_fz
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! Core mapper: homochoric h -> mapped h_FZ

  IMPLICIT NONE

  real(real64), INTENT(IN)  :: h(:), tau
  real(real64), INTENT(IN)  :: num(:), den(:)
  real(real64), INTENT(OUT) :: out(3)
  
  real(real64)              :: x, y, z, zsign, za, rho, xy, theta, theta_fz, R, rho_p
  real(real64)              :: az, s, c
  
  x = h(1)
  y = h(2)
  z = h(3)
  
  if (z >= 0.0_real64) then
    zsign = 1.0_real64
  else
    zsign = -1.0_real64
  end if
  za = abs(z)
  
  rho = norm2([x, y, z])
  xy = hypot(x, y)
  
  theta = atan2(xy, za)
  theta_fz = theta_fz_from_rational(theta, num, den)
  R = R_of_theta(theta_fz, tau)
  rho_p = rho * (R / BALL_RADIUS)
  
  az = atan2(y, x)
  s = sin(theta_fz)
  c = cos(theta_fz)
  
  out(1) = rho_p * s * cos(az)
  out(2) = rho_p * s * sin(az)
  out(3) = rho_p * c * zsign

end subroutine map_ho_to_fz

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Public API functions
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
subroutine KRcyclic_single(h, out, k)
  !DEC$ ATTRIBUTES DLLEXPORT :: KRcyclic_single
  !! author: Z. Varley, adapted for EMsoftOO by MDG
  !! version: 1.0
  !! date: 04/05/26
  !!    
  !! single conversion

IMPLICIT NONE

  real(real64), INTENT(IN)      :: h(3)
  real(real64), INTENT(OUT)     :: out(3)
  integer(kind=irg), INTENT(IN) :: k
  
  select case(k)
    case(2) 
      call map_ho_to_fz(h, tau_C2, num_C2, den_C2, out)

    case(3)
      call map_ho_to_fz(h, tau_C3, num_C3, den_C3, out)
    
    case(4)
      call map_ho_to_fz(h, tau_C4, num_C4, den_C4, out)

    case(6)
      call map_ho_to_fz(h, tau_C6, num_C6, den_C6, out)

    case default 
  end select  
end subroutine KRcyclic_single

!--------------------------------------------------------------------------
subroutine KRcyclic_batch(h_in, h_out, n, k)
  !DEC$ ATTRIBUTES DLLEXPORT :: KRcyclic_batch
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
    call KRcyclic_single(h_in(:, i), h_out(:, i), k)
  end do
!$omp end parallel do

end subroutine KRcyclic_batch

end module mod_KRcyclic

