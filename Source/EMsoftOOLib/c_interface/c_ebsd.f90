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

module c_ebsd
  !! author: MDG
  !! version: 1.0
  !! date: 03/27/26
  !!
  !! C-interop wrappers for EBSD pattern computation from master patterns.
  !! Computes simulated EBSD patterns by interpolating a master pattern
  !! on a modified Lambert projection for a given detector geometry
  !! and crystal orientation.

use iso_c_binding
use mod_kinds
use mod_global
use mod_quaternions
use mod_Lambert

IMPLICIT NONE

private

contains

!--------------------------------------------------------------------------
subroutine c_ebsd_compute_detector(numsx, numsy, xpc, ypc, delta, thetac, &
    omega, L, rgx, rgy, rgz) bind(c, name='emsoft_ebsd_compute_detector')
  !! Compute detector direction cosines for each pixel.
  !! Call once per detector geometry, then reuse for all orientations.
  !!
  !! Output arrays rgx, rgy, rgz must be pre-allocated to (numsx, numsy).
  !! Direction cosines are normalized unit vectors from sample to each pixel.
  integer(c_int), value, INTENT(IN)  :: numsx, numsy
  real(c_double), value, INTENT(IN)  :: xpc, ypc
  real(c_double), value, INTENT(IN)  :: delta
  real(c_double), value, INTENT(IN)  :: thetac
  real(c_double), value, INTENT(IN)  :: omega
  real(c_double), value, INTENT(IN)  :: L
  real(c_double), INTENT(OUT)        :: rgx(numsx, numsy)
  real(c_double), INTENT(OUT)        :: rgy(numsx, numsy)
  real(c_double), INTENT(OUT)        :: rgz(numsx, numsy)

  integer(kind=irg)                  :: i, j
  real(kind=dbl)                     :: dc(3), scin, alp, ca, sa, pcx, pcy
  real(kind=dbl)                     :: Ls, x, y, z, rr

  ! Detector tilt angle
  alp = 0.5D0 * cPi - (thetac - omega) * dtor
  ca = dcos(alp)
  sa = dsin(alp)

  ! Pixel coordinates relative to pattern center
  ! L is in microns (delta is in microns), convert consistently
  Ls = L * 1000.D0  ! convert mm to microns

  do j = 1, numsy
    do i = 1, numsx
      ! Pixel position relative to pattern center
      pcx = (dble(i) - 0.5D0 - dble(numsx) * 0.5D0 - xpc) * delta
      pcy = (dble(j) - 0.5D0 - dble(numsy) * 0.5D0 - ypc) * delta

      ! Direction cosine from sample to pixel, accounting for detector tilt
      x = pcx
      y = Ls * ca + pcy * sa
      z = -Ls * sa + pcy * ca

      ! Normalize
      rr = dsqrt(x*x + y*y + z*z)
      rgx(i, j) = x / rr
      rgy(i, j) = y / rr
      rgz(i, j) = z / rr
    end do
  end do

end subroutine c_ebsd_compute_detector

!--------------------------------------------------------------------------
subroutine c_ebsd_compute_pattern(numsx, numsy, npx, rgx, rgy, rgz, &
    mLPNH, mLPSH, quat, pattern) bind(c, name='emsoft_ebsd_compute_pattern')
  !! Compute a single EBSD pattern by interpolating the master pattern.
  !!
  !! For each detector pixel, the direction cosine is rotated by the
  !! orientation quaternion, projected onto the Lambert square, and
  !! the master pattern is bilinearly interpolated.
  !!
  !! Master patterns (mLPNH, mLPSH) should be 2D arrays of size
  !! (2*npx+1, 2*npx+1), already summed over energy bins if applicable.
  integer(c_int), value, INTENT(IN)  :: numsx, numsy
  integer(c_int), value, INTENT(IN)  :: npx
  real(c_double), INTENT(IN)         :: rgx(numsx, numsy)
  real(c_double), INTENT(IN)         :: rgy(numsx, numsy)
  real(c_double), INTENT(IN)         :: rgz(numsx, numsy)
  real(c_double), INTENT(IN)         :: mLPNH(2*npx+1, 2*npx+1)
  real(c_double), INTENT(IN)         :: mLPSH(2*npx+1, 2*npx+1)
  real(c_double), INTENT(IN)         :: quat(4)
  real(c_double), INTENT(OUT)        :: pattern(numsx, numsy)

  type(Quaternion_T)                 :: qu
  integer(kind=irg)                  :: i, j, nix, niy, nixp, niyp, npxi
  real(kind=dbl)                     :: dc(3), dcr(3), rr, scl
  real(kind=dbl)                     :: dx, dy, dxm, dym, xy(2)
  real(kind=dbl)                     :: sq2pi
  integer(kind=irg)                  :: ierr

  ! Set up quaternion for rotation
  qu = Quaternion_T( qd = quat )

  ! Scale factor for Lambert indexing
  scl = dble(npx)
  sq2pi = dsqrt(cPi * 0.5D0)
  npxi = npx + 1  ! 1-based index offset (array goes from 1 to 2*npx+1)

  do j = 1, numsy
    do i = 1, numsx
      ! Get detector direction cosine
      dc = (/ rgx(i,j), rgy(i,j), rgz(i,j) /)

      ! Rotate into crystal reference frame
      dcr = qu%quat_Lp(dc)

      ! Normalize
      rr = dsqrt(dcr(1)**2 + dcr(2)**2 + dcr(3)**2)
      dcr = dcr / rr

      ! Lambert projection: sphere to square
      if (dabs(dcr(3)) .ge. 1.D0) then
        ! At a pole
        nix = 0
        niy = 0
        nixp = 0
        niyp = 0
        dx = 0.D0
        dy = 0.D0
        dxm = 1.D0
        dym = 1.D0
      else
        ! Compute Lambert square coordinates
        if (dabs(dcr(2)) .le. dabs(dcr(1))) then
          xy(1) = dsign(1.D0, dcr(1)) * dsqrt(2.D0 * (1.D0 - dabs(dcr(3)))) / sq2pi
          xy(2) = xy(1) * datan2(dcr(2), dcr(1)) / (cPi * 0.25D0)
        else
          xy(2) = dsign(1.D0, dcr(2)) * dsqrt(2.D0 * (1.D0 - dabs(dcr(3)))) / sq2pi
          xy(1) = xy(2) * datan2(dcr(1), dcr(2)) / (cPi * 0.25D0)
        end if

        ! Scale to grid indices
        xy = xy * scl

        ! Bilinear interpolation indices
        nix = int(xy(1) + scl) - npx
        niy = int(xy(2) + scl) - npx
        nixp = nix + 1
        niyp = niy + 1
        dx = xy(1) - dble(nix)
        dy = xy(2) - dble(niy)
        dxm = 1.D0 - dx
        dym = 1.D0 - dy

        ! Clamp to array bounds
        nix = max(1, min(2*npx+1, nix + npxi))
        niy = max(1, min(2*npx+1, niy + npxi))
        nixp = max(1, min(2*npx+1, nixp + npxi))
        niyp = max(1, min(2*npx+1, niyp + npxi))
      end if

      ! Bilinear interpolation
      if (dcr(3) .ge. 0.D0) then
        ! Northern hemisphere
        pattern(i,j) = mLPNH(nix,niy) * dxm * dym &
                     + mLPNH(nixp,niy) * dx * dym &
                     + mLPNH(nix,niyp) * dxm * dy &
                     + mLPNH(nixp,niyp) * dx * dy
      else
        ! Southern hemisphere
        pattern(i,j) = mLPSH(nix,niy) * dxm * dym &
                     + mLPSH(nixp,niy) * dx * dym &
                     + mLPSH(nix,niyp) * dxm * dy &
                     + mLPSH(nixp,niyp) * dx * dy
      end if
    end do
  end do

end subroutine c_ebsd_compute_pattern

!--------------------------------------------------------------------------
subroutine c_ebsd_compute_patterns(numsx, numsy, npx, rgx, rgy, rgz, &
    mLPNH, mLPSH, quats, nquats, patterns) &
    bind(c, name='emsoft_ebsd_compute_patterns')
  !! Compute multiple EBSD patterns for an array of orientations.
  !! patterns must be pre-allocated to (numsx, numsy, nquats).
  integer(c_int), value, INTENT(IN)  :: numsx, numsy, npx, nquats
  real(c_double), INTENT(IN)         :: rgx(numsx, numsy)
  real(c_double), INTENT(IN)         :: rgy(numsx, numsy)
  real(c_double), INTENT(IN)         :: rgz(numsx, numsy)
  real(c_double), INTENT(IN)         :: mLPNH(2*npx+1, 2*npx+1)
  real(c_double), INTENT(IN)         :: mLPSH(2*npx+1, 2*npx+1)
  real(c_double), INTENT(IN)         :: quats(4, nquats)
  real(c_double), INTENT(OUT)        :: patterns(numsx, numsy, nquats)

  integer(kind=irg)                  :: iq

  do iq = 1, nquats
    call c_ebsd_compute_pattern(numsx, numsy, npx, rgx, rgy, rgz, &
        mLPNH, mLPSH, quats(:, iq), patterns(:, :, iq))
  end do

end subroutine c_ebsd_compute_patterns

end module c_ebsd
