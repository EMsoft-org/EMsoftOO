! ###################################################################
! Copyright (c) 2013-2024, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_Lambert
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! everything that has to do with the modified Lambert projections.
  !!
  !! This module contains a number of projection functions for the modified
  !! Lambert projection between square lattice and 2D hemisphere, hexagonal lattice
  !! and 2D hemisphere, as well as the more complex mapping between a 3D cubic grid
  !! and the unit quaternion hemisphere with positive scalar component.  In addition, there
  !! are some other projections, such as the stereographic one.  Each function is named
  !! by the projection, the dimensionality of the starting grid, and the forward or inverse
  !! character.  For each function, there is also a single precision and a double precision
  !! version, but we use the generic function formalism to have only a single call.  The Forward
  !! mapping is taken to be the one from the simple grid to the curved grid.  Since the module
  !! deals with various grids, we also add a few functions/subroutines that apply symmetry
  !! operations on those grids.  Finally, since Lambert grid interpolation is very common in
  !! EMsoft programs, we provide several versions that work with a variety of master pattern
  !! types.
  !!
  !! original modification history
  !! @date 07/10/13   MDG 1.0 original
  !! @date 07/12/13   MDG 1.1 added forward cube to ball to quaternion mappings
  !! @date 08/01/13   MDG 1.2 added standard Lambert projection
  !! @date 08/12/13   MDG 1.3 added inverse Lambert projections for Ball to Cube
  !! @date 09/20/13   MDG 1.4 added ApplyLaueSymmetry
  !! @date 08/29/15   MDG 1.5 small changes to hexagonal mapping routines; coordinate swap inside routines
  !! @date 01/20/18   MDG 1.6 added Lambert interpolation routines

use mod_kinds
use mod_global
use mod_quaternions

IMPLICIT NONE
  private

  public :: LambertgetInterpolation
  interface LambertgetInterpolation
          module procedure LambertgetInterpolationSingle
          module procedure LambertgetInterpolationDouble
  end interface

  public :: InterpolateLambert
  interface InterpolateLambert
          module procedure InterpolationLambert2DSingle
          module procedure InterpolationLambert2DDouble
          module procedure InterpolationLambert3DSingle
          module procedure InterpolationLambert3DInteger
          module procedure InterpolationLambert4DSingle
          module procedure InterpolationLambert4DDouble4b4
  end interface

! auxiliary private functions for the hexagonal mappings
  private :: GetSextantSingle, GetSextantDouble, GetPyramidSingle, GetPyramidDouble

  type, public :: Lambert_T
    !! Lambert class definition
    private
      real(kind=sgl), dimension(2)    :: xy
       !! single precision 2D coordinate pair
      real(kind=sgl), dimension(3)    :: xyz
       !! single precision 3D coordinate triplet
      real(kind=dbl), dimension(2)    :: xyd
       !! double precision 2D coordinate pair
      real(kind=dbl), dimension(3)    :: xyzd
       !! double precision 3D coordinate triplet
      character(1)                    :: s
       !! precision indicator ('s' or 'd')

    contains
      private
        procedure, pass(self) :: Lambert2DSquareForwardSingle
        procedure, pass(self) :: Lambert2DSquareForwardDouble
        procedure, pass(self) :: Lambert2DSquareInverseSingle
        procedure, pass(self) :: Lambert2DSquareInverseDouble
        procedure, pass(self) :: Lambert2DHexForwardSingle
        procedure, pass(self) :: Lambert2DHexForwardDouble
        procedure, pass(self) :: Lambert2DHexInverseSingle
        procedure, pass(self) :: Lambert2DHexInverseDouble
        procedure, pass(self) :: Lambert3DCubeForwardSingle
        procedure, pass(self) :: Lambert3DCubeForwardDouble
        procedure, pass(self) :: Lambert3DCubeInverseSingle
        procedure, pass(self) :: Lambert3DCubeInverseDouble
        procedure, pass(self) :: Lambert3DBallToQuaternion
!       procedure, pass(self) :: Lambert3DBallToQuaternionSingleInverse
!       procedure, pass(self) :: Lambert3DBallToQuaternionDoubleInverse
        procedure, pass(self) :: Lambert3DCubeToQuaternion
!       procedure, pass(self) :: Lambert3DCubeToQuaternionSingleInverse
!       procedure, pass(self) :: Lambert3DCubeToQuaternionDoubleInverse
        procedure, pass(self) :: StereoGraphicForwardSingle
        procedure, pass(self) :: StereoGraphicForwardDouble
        procedure, pass(self) :: StereoGraphicInverseSingle
        procedure, pass(self) :: StereoGraphicInverseDouble
        procedure, pass(self) :: LambertForwardSingle
        procedure, pass(self) :: LambertForwardDouble
        procedure, pass(self) :: LambertInverseSingle
        procedure, pass(self) :: LambertInverseDouble
        procedure, pass(self) :: Apply3DPGSymmetry_
        procedure, pass(self) :: sampleVMF_
        procedure, pass(self) :: HemiCheck_
        procedure, pass(self) :: setxy_
        procedure, pass(self) :: setxyd_
        procedure, pass(self) :: setxyz_
        procedure, pass(self) :: setxyzd_
        final :: Lambert_destructor

! mappings from 2D square grid to the Northern hemisphere of a 2D sphere
        generic, public :: LambertSquareToSphere =>  Lambert2DSquareForwardSingle, Lambert2DSquareForwardDouble
        generic, public :: LambertSphereToSquare => Lambert2DSquareInverseSingle, Lambert2DSquareInverseDouble
! mappings from 2D hexagonal grid to the Northern hemisphere of a 2D sphere
        generic, public :: LambertHexToSphere => Lambert2DHexForwardSingle, Lambert2DHexForwardDouble
        generic, public :: LambertSphereToHex => Lambert2DHexInverseSingle, Lambert2DHexInverseDouble
! mappings from the 3D cubic grid to the 3D spherical grid
        generic, public :: LambertCubeToBall => Lambert3DCubeForwardSingle, Lambert3DCubeForwardDouble
        generic, public :: LambertBallToCube => Lambert3DCubeInverseSingle, Lambert3DCubeInverseDouble
! mappings from the 3D spherical grid to the unit quaternion sphere
        generic, public :: LambertBallToQuaternion => Lambert3DBallToQuaternion
        ! generic, public :: LambertQuaternionToBall => Lambert3DBallToQuaternionSingleInverse, &
                                                        ! Lambert3DBallToQuaternionDoubleInverse
! !DEC$ ATTRIBUTES DLLEXPORT :: LambertQuaternionToBall
! mappings from the 3D cube grid to the unit quaternion sphere
        generic, public :: LambertCubeToQuaternion => Lambert3DCubeToQuaternion
        ! generic, public :: LambertQuaternionToCube => Lambert3DCubeToQuaternionSingleInverse, &
                                                        ! Lambert3DCubeToQuaternionDoubleInverse
! !DEC$ ATTRIBUTES DLLEXPORT :: LambertQuaternionToCube
! simple stereographic projection
        generic, public :: StereoGraphicForward => StereoGraphicForwardSingle, StereoGraphicForwardDouble
        generic, public :: StereoGraphicInverse => StereoGraphicInverseSingle, StereoGraphicInverseDouble
! simple Lambert projection
        generic, public :: LambertForward => LambertForwardSingle, LambertForwardDouble
        generic, public :: LambertInverse => LambertInverseSingle, LambertInverseDouble
        generic, public :: Apply3DPGSymmetry => Apply3DPGSymmetry_
        generic, public :: sampleVMF => sampleVMF_
        generic, public :: HemiCheck => HemiCheck_
        generic, public :: setxy => setxy_
        generic, public :: setxyz => setxyz_
        generic, public :: setxyd => setxyd_
        generic, public :: setxyzd => setxyzd_

  end type Lambert_T

  ! the constructor routine for this class
  interface Lambert_T
    module procedure Lambert_constructor
  end interface Lambert_T

contains

!--------------------------------------------------------------------------
type(Lambert_T) function Lambert_constructor( xy, xyz, xyd, xyzd ) result(Lambert)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert_constructor
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! constructor for the Lambert Class

IMPLICIT NONE

  real(kind=sgl), INTENT(IN), OPTIONAL      :: xy(2)
  real(kind=sgl), INTENT(IN), OPTIONAL      :: xyz(3)
  real(kind=dbl), INTENT(IN), OPTIONAL      :: xyd(2)
  real(kind=dbl), INTENT(IN), OPTIONAL      :: xyzd(3)

! fill in one or the other coordinate set
  if ( ((.not.present(xy)).and.(.not.present(xyd))) .and. ((.not.present(xyz)).and.(.not.present(xyzd)))) then
    Lambert % xy = (/  0.0, 0.0 /)
    Lambert % xyz = (/  0.0, 0.0, 0.0 /)
    Lambert % xyd = (/  0.D0, 0.D0 /)
    Lambert % xyzd = (/  0.D0, 0.D0, 0.D0 /)
    Lambert % s = 's'
  else
      if (present(xy)) then
        Lambert % xy = xy
        Lambert % xyd = (/ 0.D0, 0.D0 /)
        Lambert % xyz = (/ 0.0, 0.0, 0.0 /)
        Lambert % xyzd = (/ 0.D0, 0.D0, 0.D0 /)
        Lambert % s = 's'
      end if
      if (present(xyd)) then
        Lambert % xy =  (/ 0.0, 0.0 /)
        Lambert % xyd = xyd
        Lambert % xyz = (/ 0.0, 0.0, 0.0 /)
        Lambert % xyzd = (/ 0.D0, 0.D0, 0.D0 /)
        Lambert % s = 'd'
      end if
      if (present(xyz)) then
        Lambert % xy = (/ 0.0, 0.0 /)
        Lambert % xyd = (/ 0.D0, 0.D0 /)
        Lambert % xyz = xyz
        Lambert % xyzd = (/ 0.D0, 0.D0, 0.D0 /)
        Lambert % s = 's'
      end if
      if (present(xyzd)) then
        Lambert % xy = (/ 0.0, 0.0 /)
        Lambert % xyd = (/ 0.D0, 0.D0 /)
        Lambert % xyz = (/ 0.0, 0.0, 0.0 /)
        Lambert % xyzd = xyzd
        Lambert % s = 'd'
      end if
  end if

end function Lambert_constructor

!--------------------------------------------------------------------------
subroutine Lambert_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert_destructor
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! Lambert class destructor

IMPLICIT NONE

type(Lambert_T), INTENT(INOUT)  :: self

! this produces too many calls so we disable the output
! call reportDestructor('Lambert_T')

end subroutine Lambert_destructor

!--------------------------------------------------------------------------
subroutine setxy_(self, xy)
!DEC$ ATTRIBUTES DLLEXPORT :: setxy_
  !! author: MDG
  !! version: 1.0
  !! date: 02/17/20
  !!
  !! set a coordinate

IMPLICIT NONE

class(Lambert_T), INTENT(INOUT)  :: self
real(kind=sgl), INTENT(IN)       :: xy(2)

self%xy = xy
self%s = 's'

end subroutine setxy_

!--------------------------------------------------------------------------
subroutine setxyz_(self, xyz)
!DEC$ ATTRIBUTES DLLEXPORT :: setxyz_
  !! author: MDG
  !! version: 1.0
  !! date: 02/17/20
  !!
  !! set a coordinate

IMPLICIT NONE

class(Lambert_T), INTENT(INOUT)  :: self
real(kind=sgl), INTENT(IN)       :: xyz(3)

self%xyz = xyz
self%s = 's'

end subroutine setxyz_

!--------------------------------------------------------------------------
subroutine setxyd_(self, xyd)
!DEC$ ATTRIBUTES DLLEXPORT :: setxyd_
  !! author: MDG
  !! version: 1.0
  !! date: 02/17/20
  !!
  !! set a coordinate

IMPLICIT NONE

class(Lambert_T), INTENT(INOUT)  :: self
real(kind=dbl), INTENT(IN)       :: xyd(2)

self%xyd = xyd
self%s = 'd'

end subroutine setxyd_

!--------------------------------------------------------------------------
subroutine setxyzd_(self, xyzd)
!DEC$ ATTRIBUTES DLLEXPORT :: setxyzd_

  !! author: MDG
  !! version: 1.0
  !! date: 02/17/20
  !!
  !! set a coordinate

IMPLICIT NONE

class(Lambert_T), INTENT(INOUT)  :: self
real(kind=dbl), INTENT(IN)       :: xyzd(3)

self%xyzd = xyzd
self%s = 'd'

end subroutine setxyzd_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! We begin with the functions/subroutines that are public in this class
! Note that the result and error parameters have been swapped in all routines
! below compared to the original Lambert module !!!
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive function Lambert2DSquareForwardSingle(self, res) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert2DSquareForwardSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! forward projection from 2D square to 3D hemisphere, single precision

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=sgl),intent(out)              :: res(3)
   !! transformation result
  integer(kind=irg)                       :: ierr
   !f2py intent(in,out) ::  ierr

  real(kind=sgl)                          :: q, qq, xy2(2)

xy2 = self%xy * sngl(LPs%sPio2)

ierr = 0
! check to make sure that the input point lies inside the square of edge length 2 sqrt(pi/2)
if (maxval(abs(xy2)).gt.LPs%sPio2) then
  res = (/ 0.0, 0.0, 0.0 /)
  ierr = 1
  return
else
! Forward projection from square grid to Northern hemisphere.
! Equations (8) and (9) from D. Rosca, "New uniform grids on the sphere,"
! Astronomy & Astrophysics, 520, A63 (2010)

! deal with the origin first:
 if (maxval(abs(xy2)).eq.0.0) then
   res = (/ 0.0, 0.0, 1.0 /)
 else
  if (abs(xy2(1)).le.abs(xy2(2))) then
   q = 2.0*xy2(2)*LPs%iPi*sqrt(LPs%Pi-xy2(2)*xy2(2))
   qq = xy2(1)*LPs%Pi*0.25/xy2(2)
   res = (/ q*sin(qq), q*cos(qq), 1.0-2.0*xy2(2)*xy2(2)*sngl(LPs%iPi) /)
  else
   q = 2.0*xy2(1)*LPs%iPi*sqrt(LPs%Pi-xy2(1)*xy2(1))
   qq = xy2(2)*LPs%Pi*0.25/xy2(1)
   res = (/ q*cos(qq), q*sin(qq), 1.0-2.0*xy2(1)*xy2(1)*sngl(LPs%iPi) /)
  end if
  res = res/sqrt(sum(res*res))
 end if
end if

end function Lambert2DSquareForwardSingle

!--------------------------------------------------------------------------
recursive function Lambert2DSquareForwardDouble(self, res) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert2DSquareForwardDouble
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! forward projection from 2D square to 3D hemisphere, double precision

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=dbl),intent(out)              :: res(3)
   !! transformation result
  integer(kind=irg)                       :: ierr
   !f2py intent(in,out) ::  ierr

  real(kind=dbl)                          :: q, qq, xy2(2)

xy2 = self%xyd * LPs%sPio2

ierr = 0
! check to make sure that the input point lies inside the square of edge length 2 sqrt(pi)
if (maxval(dabs(xy2)).gt.LPs%sPio2) then
  res = (/ 0.D0, 0.D0, 0.D0 /)
  ierr = 1   ! input point does not lie inside square with edge length 2 sqrt(pi/2)
  return
else
! Forward projection from square grid to Northern hemisphere.
! Equations (8) and (9) from D. Rosca, "New uniform grids on the sphere,"
! Astronomy & Astrophysics, 520, A63 (2010)

! deal with the origin first:
 if (maxval(abs(xy2)).eq.0.0) then
   res = (/ 0.D0, 0.D0, 1.D0 /)
 else
  if (dabs(xy2(1)).le.dabs(xy2(2))) then
   q = 2.D0*xy2(2)*LPs%iPi*dsqrt(LPs%Pi-xy2(2)*xy2(2))
   qq = xy2(1)*LPs%Pi*0.25D0/xy2(2)
   res = (/ q*dsin(qq), q*dcos(qq), 1.D0-2.D0*xy2(2)*xy2(2)*LPs%iPi /)
  else
   q = 2.D0*xy2(1)*LPs%iPi*dsqrt(LPs%Pi-xy2(1)*xy2(1))
   qq = xy2(2)*LPs%Pi*0.25D0/xy2(1)
   res = (/ q*dcos(qq), q*dsin(qq), 1.D0-2.D0*xy2(1)*xy2(1)*LPs%iPi /)
  end if
  res = res/dsqrt(sum(res*res))
 end if
end if

end function Lambert2DSquareForwardDouble

!--------------------------------------------------------------------------
recursive function Lambert2DSquareInverseSingle(self, res) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert2DSquareInverseSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! inverse projection from 3D hemisphere to 2D square, single precision
  !!
  !! IMPORTANT NOTE: the calling routine must keep track of the sign of xyz(3);
  !! this routine only uses the absolute value.

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=sgl),intent(out)              :: res(2)
   !! output coordinate pair
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=sgl)                          :: q
  real(kind=sgl),parameter                :: eps = 1.0E-6

ierr = 0
! check to make sure that the input point lies on the unit sphere
if (abs(1.0-sum(self%xyz**2)).gt.eps) then
  res = (/ 0.0, 0.0 /)
  ierr = 1
else
! intercept the points (0,0,+-1)
  if (abs(self%xyz(3)).eq.1.0) then
    res = (/ 0.0, 0.0 /)
  else
    if (abs(self%xyz(2)).le.abs(self%xyz(1))) then
      q = abs(self%xyz(1))/self%xyz(1) * sqrt(2.0*(1.0-abs(self%xyz(3))))
      res = (/ q * sngl(LPs%sPi2), q * atan(self%xyz(2)/self%xyz(1))/sngl(LPs%sPi2) /)
    else
      q = abs(self%xyz(2))/self%xyz(2) * sqrt(2.0*(1.0-abs(self%xyz(3))))
      res = (/  q * atan(self%xyz(1)/self%xyz(2))/sngl(LPs%sPi2), q * sngl(LPs%sPi2) /)
    end if
  end if
end if

res = res / sngl(LPs%sPio2)

end function Lambert2DSquareInverseSingle

!--------------------------------------------------------------------------
recursive function Lambert2DSquareInverseDouble(self, res) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert2DSquareInverseDouble
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! inverse projection from 3D hemisphere to 2D square, double precision
  !!
  !! IMPORTANT NOTE: the calling routine must keep track of the sign of xyz(3);
  !! this routine only uses the absolute value.

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=dbl),intent(out)              :: res(2)
   !! output coordinate pair
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=dbl)                          :: q
  real(kind=dbl),parameter                :: eps = 1.0D-12

ierr = 0
! check to make sure that the input point lies on the unit sphere
if (dabs(1.D0-sum(self%xyzd**2)).gt.eps) then
  res = (/ 0.D0, 0.D0 /)
  ierr = 1
else
! intercept the points (0,0,+-1)
  if (dabs(self%xyzd(3)).eq.1.D0) then
    res = (/ 0.D0, 0.D0 /)
  else
    if (dabs(self%xyzd(2)).le.dabs(self%xyzd(1))) then
      q = dabs(self%xyzd(1))/self%xyzd(1) * dsqrt(2.D0*(1.D0-dabs(self%xyzd(3))))
      res = (/ q * LPs%sPi2, q * datan(self%xyzd(2)/self%xyzd(1))/LPs%sPi2 /)
    else
      q = dabs(self%xyzd(2))/self%xyzd(2) * dsqrt(2.D0*(1.D0-dabs(self%xyzd(3))))
      res = (/  q * datan(self%xyzd(1)/self%xyzd(2))/LPs%sPi2, q * LPs%sPi2 /)
    end if
  end if
end if

res = res / LPs%sPio2

end function Lambert2DSquareInverseDouble

!--------------------------------------------------------------------------
! the functions below deal with the hexagonal to 2D hemisphere projection
!
! all derivations and equations can be found in
!
! D. Rosca and M. De Graef, "Area-preserving projections from hexagonal and triangular
! domains to the sphere and applications to electron back-scatter diffraction pattern simulations,"
! Modelling Simul. Mater. Sci. Eng. 21 (2013) 055021.
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive function Lambert2DHexForwardSingle(self, res) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert2DHexForwardSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! forward projection from 2D hexagon to 3D hemisphere, single precision

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=sgl),intent(out)              :: res(3)
   !! output coordinate triplet
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=sgl)                          :: q, XX, YY, xp, yp, XY2(2), xyc(2)
  integer(kind=irg)                       :: ks

 ierr = 0

 xyc = (/ self%xy(1) - 0.5 * self%xy(2), self%xy(2) * sngl(LPS%srt) /) * sngl(LPs%preg)

 if (maxval(abs(xyc)).eq.0.0) then
  res = (/ 0.0, 0.0, 1.0 /)
 else
! flip coordinates
  XY2 = (/ xyc(2), xyc(1) /)
! determine in which sextant this point lies
  ks = GetSextantSingle(XY2)

  select case (ks)
    case (0,3)
        q = XY2(2)*LPs%prec/XY2(1)
        XX = LPs%preb*XY2(1)*cos(q)
        YY = LPs%preb*XY2(1)*sin(q)
    case (1,4)
        xp = XY2(1)+LPs%rtt*XY2(2)
        yp = XY2(1)*LPs%pred/xp
        XX = LPs%prea*xp*sin(yp)
        YY = LPs%prea*xp*cos(yp)
    case (2,5)
        xp = XY2(1)-LPs%rtt*XY2(2)
        yp = XY2(1)*LPs%pred/xp
        XX = LPs%prea*xp*sin(yp)
        YY = -LPs%prea*xp*cos(yp)
  end select
  q = XX**2+YY**2
! does the point lie outside the hexagon ?
  if (q.gt.4.0) then
    res = (/ 0.0, 0.0, 0.0 /)
    ierr = 1
  else
    res = (/ 0.50*XX*sqrt(4.0-q), 0.50*YY*sqrt(4.0-q), 1.0-0.5*q /)
  end if

! and flip the x and y coordinates
  xp = res(1)
  res(1) = res(2)
  res(2) = xp
 end if

end function Lambert2DHexForwardSingle

!--------------------------------------------------------------------------
recursive function Lambert2DHexForwardDouble(self, res) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert2DHexForwardDouble
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! forward projection from 2D hexagon to 3D hemisphere, double precision

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=dbl),intent(out)              :: res(3)
   !! output coordinate triplet
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=dbl)                          :: q, XX, YY, xp, yp, XY2(2), xyc(2)
  integer(kind=irg)                       :: ks

  ierr = 0

 xyc = (/ self%xyd(1) - 0.5D0 * self%xyd(2), self%xyd(2) * LPS%srt /) * LPs%preg

 if (maxval(dabs(xyc)).eq.0.D0) then
  res = (/ 0.D0, 0.D0, 1.D0 /)
 else
! flip coordinates
  XY2 = (/ xyc(2), xyc(1) /)
! determine in which sextant this point lies
  ks = GetSextantDouble(XY2)

  select case (ks)
    case (0,3)
        q = XY2(2)*LPs%prec/XY2(1)
        XX = LPs%preb*XY2(1)*dcos(q)
        YY = LPs%preb*XY2(1)*dsin(q)
    case (1,4)
        xp = XY2(1)+LPs%rtt*XY2(2)
        yp = XY2(1)*LPs%pred/xp
        XX = LPs%prea*xp*dsin(yp)
        YY = LPs%prea*xp*dcos(yp)
    case (2,5)
        xp = XY2(1)-LPs%rtt*XY2(2)
        yp = XY2(1)*LPs%pred/xp
        XX = LPs%prea*xp*dsin(yp)
        YY = -LPs%prea*xp*dcos(yp)
  end select
  q = XX**2+YY**2

! does the point lie outside the hexagon ?
  if (q.gt.4.D0) then
    res = (/ 0.D0, 0.D0, 0.D0 /)
    ierr = 1
  else
    res = (/ 0.5D0*XX*dsqrt(4.D0-q), 0.5D0*YY*dsqrt(4.D0-q), 1.D0-0.5D0*q /)
  end if

! and flip the x and y coordinates
  xp = res(1)
  res(1) = res(2)
  res(2) = xp
 end if

end function Lambert2DHexForwardDouble

!--------------------------------------------------------------------------
recursive function Lambert2DHexInverseSingle(self, res) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert2DHexInverseSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! forward projection from 2D hexagon to 3D hemisphere, single precision
  !!
  !! IMPORTANT NOTE: The calling program must keep track of the sign of xyz(3), since this
  !! routine will take the absolute value |xyz(3)| for the z-component of the input vector!

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=sgl),intent(out)              :: res(2)
   !! output coordinate triplet
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=sgl)                          :: q, qq, XX, YY, xxx, yyy, sgnX, XYZ2(3), xy(2)
  integer(kind=irg)                       :: ks
  real(kind=sgl),parameter                :: eps = 1.0E-7, eps2 = 1.0E-4

ierr = 0
! check to make sure that the input point lies on the unit sphere
if (abs(1.0-sum(self%xyz**2)).gt.eps) then
  res = (/ 0.0, 0.0 /)
  ierr = 1
else
 if (abs(self%xyz(3)).eq.1.0) then
  res = (/ 0.0, 0.0 /)
 else
! flip x and y components, and take the | | of the third component.
  XYZ2 = (/ self%xyz(2), self%xyz(1), abs(self%xyz(3)) /)

! first do the Lambert projection
  q = sqrt(2.0/(1.0+XYZ2(3)))
  XX = q * XYZ2(1)+eps2
  YY = q * XYZ2(2)+eps2

! determine in which sextant this point lies
  ks = GetSextantSingle( (/ XX, YY /) )
  sgnX = sqrt(XX**2+YY**2)
  if (XX.lt.0.0) sgnX=-sgnX

! then perform the inverse to the hexagonal grid
  select case (ks)
    case (0,3)
        q = LPs%pree * sgnX
        xxx = q * LPs%sPio2
        if (XX.eq.0.0) then
          yyy = LPs%pref * LPs%Pi * 0.5
        else
          yyy = q * LPs%pref * atan(YY/XX)
        end if
    case (1,4)
        q = LPs%prea * sgnX
        qq= atan((YY-LPs%rtt*XX)/(XX+LPs%rtt*YY))
        xxx = q * LPs%rtt *( LPs%Pi/6.0 - qq )
        yyy = q * ( 0.5*LPs%Pi + qq )
    case (2,5)
        q = LPs%prea * sgnX
        qq= atan((YY+LPs%rtt*XX)/(XX-LPs%rtt*YY))
        xxx = q * LPs%rtt *( LPs%Pi/6.0 + qq )
        yyy = q * ( -0.5*LPs%Pi + qq )
  end select
! and flip the coordinates
  res = (/ yyy, xxx /)
 end if
end if

! and finally, transform the coordinates back to the hexagonal grid
xy = res
res = (/ xy(1) + xy(2) * sngl(LPS%isrt), xy(2) * 2.0 * sngl(LPS%isrt) /) / LPs%preg

end function Lambert2DHexInverseSingle

!--------------------------------------------------------------------------
recursive function Lambert2DHexInverseDouble(self, res) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert2DHexInverseDouble
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! forward projection from 2D hexagon to 3D hemisphere, single precision
  !!
  !! IMPORTANT NOTE: The calling program must keep track of the sign of xyz(3), since this
  !! routine will take the absolute value |xyz(3)| for the z-component of the input vector!

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=dbl),intent(out)              :: res(2)
   !! output coordinate triplet
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=dbl)                          :: q, qq, XX, YY, xxx, yyy, sgnX, XYZ2(3), xy(2)
  integer(kind=irg)                       :: ks
  real(kind=dbl),parameter                :: eps = 1.0E-12, eps2 = 1.0E-4

ierr = 0
! check to make sure that the input point lies on the unit sphere
if (dabs(1.0-sum(self%xyzd**2)).gt.eps) then
  res = (/ 0.D0, 0.D0 /)
  ierr = 1
else
 if (dabs(self%xyzd(3)).eq.1.D0) then
  res = (/ 0.D0, 0.D0 /)
 else
! flip x and y components, and take the | | of the third component.
  XYZ2 = (/ self%xyzd(2), self%xyzd(1), dabs(self%xyzd(3)) /)

! first do the Lambert projection
  q = dsqrt(2.D0/(1.D0+XYZ2(3)))
  XX = q * XYZ2(1)+eps2
  YY = q * XYZ2(2)+eps2

! determine in which sextant this point lies
  ks = GetSextantDouble( (/ XX, YY /) )
  sgnX = dsqrt(XX**2+YY**2)
  if (XX.lt.0.D0) sgnX=-sgnX

! then perform the inverse to the hexagonal grid
  select case (ks)
    case (0,3)
        q = LPs%pree * sgnX
        xxx = q * LPs%sPio2
        if (XX.eq.0.0) then
          yyy = LPs%pref * LPs%Pi * 0.5D0
        else
          yyy = q * LPs%pref * datan(YY/XX)
        end if
    case (1,4)
        q = LPs%prea * sgnX
        qq= datan((YY-LPs%rtt*XX)/(XX+LPs%rtt*YY))
        xxx = q * LPs%rtt *( LPs%Pi/6.D0 - qq )
        yyy = q * ( 0.5D0*LPs%Pi + qq )
    case (2,5)
        q = LPs%prea * sgnX
        qq= datan((YY+LPs%rtt*XX)/(XX-LPs%rtt*YY))
        xxx = q * LPs%rtt *( LPs%Pi/6.D0 + qq )
        yyy = q * ( -0.5D0*LPs%Pi + qq )
  end select
! and flip the coordinates back
    res = (/ yyy, xxx /)
end if
end if

! and finally, transform the coordinates back to the hexagonal grid
xy = res
res = (/ xy(1) + xy(2) * LPS%isrt, xy(2) * 2.D0 * LPS%isrt /) / LPs%preg

end function Lambert2DHexInverseDouble

!--------------------------------------------------------------------------
recursive function GetSextantSingle(xy) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: GetSextantSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! determine to which sextant a point in hexagonal coordinates belongs, single precision
  !! private function

IMPLICIT NONE

real(kind=sgl),INTENT(IN)       :: xy(2)
integer(kind=irg)               :: res

real(kind=sgl)                  :: xx

xx = abs(xy(1)*LPs%isrt)        ! |x| / sqrt(3)

if (xy(1).ge.0.0) then
  if (abs(xy(2)).le.xx) then
    res = 0
  else
    if (xy(2).gt.xx) then
      res = 1
    else
      res = 5
    end if
  end if
else
  if (abs(xy(2)).le.xx) then
    res = 3
  else
    if (xy(2).gt.xx) then
      res = 2
    else
      res = 4
    end if
  end if
end if

end function GetSextantSingle

!--------------------------------------------------------------------------
recursive function GetSextantDouble(xy) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: GetSextantDouble
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! determine to which sextant a point in hexagonal coordinates belongs, double precision
  !! private function

IMPLICIT NONE

real(kind=dbl),INTENT(IN)       :: xy(2)
integer(kind=irg)               :: res

real(kind=dbl)                  :: xx

xx = dabs(xy(1)*LPs%isrt)       ! |x| / sqrt(3)

if (xy(1).ge.0.D0) then
  if (dabs(xy(2)).le.xx) then
    res = 0
  else
    if (xy(2).gt.xx) then
      res = 1
    else
      res = 5
    end if
  end if
else
  if (dabs(xy(2)).le.xx) then
    res = 3
  else
    if (xy(2).gt.xx) then
      res = 2
    else
      res = 4
    end if
  end if
end if

end function GetSextantDouble


!--------------------------------------------------------------------------
! the functions below deal with the cubic grid to the 3D ball, and then from
! the 3D ball to the unit quaternion hemisphere projection
!
! all derivations and equations can be found in
!
! D. Rosca, A. Morawiec, and M. De Graef. “A new method of constructing a grid in the space of
! 3D rotations and its applications to texture analysis”. Modeling and Simulations in Materials
! Science and Engineering 22, 075013 (2014)
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive function Lambert3DCubeForwardSingle(self, res) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert3DCubeForwardSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! map from 3D cubic grid to 3D ball, single precision

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=sgl),intent(out)              :: res(3)
   !! output coordinate triplet
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=sgl)                          :: XYZ(3), sXYZ(3), T1, T2, c, s, q, LamXYZ(3), eps
  integer(kind=irg)                       :: p

eps = 1.0D-7

ierr = 0
if (maxval(abs(self%xyz)).gt.LPs%ap/2) then
  res = (/ 0.0, 0.0, 0.0 /)
  ierr = 1
  return
end if
if (maxval(abs(self%xyz)).lt.eps) then
  res = (/ 0.0, 0.0, 0.0 /)
  return
end if

! determine which pyramid pair the point lies in and copy coordinates in correct order (see paper)
p = GetPyramidSingle(self%xyz)
select case (p)
 case (1,2)
  sXYZ = self%xyz
 case (3,4)
  sXYZ = (/ self%xyz(2), self%xyz(3), self%xyz(1) /)
 case (5,6)
  sXYZ = (/ self%xyz(3), self%xyz(1), self%xyz(2) /)
end select

! scale by grid parameter ratio sc
XYZ = LPs%sc * sXYZ

! transform to the sphere grid via the curved square, and intercept the zero point
if (maxval(abs(XYZ)).eq.0.0) then
  LamXYZ = (/ 0.0, 0.0, 0.0 /)
else
! intercept all the points along the z-axis
  if (maxval(abs(XYZ(1:2))).eq.0.0) then
    LamXYZ = (/ 0.0, 0.0, sngl(LPs%pref) * XYZ(3) /)
  else  ! this is a general grid point
    if (abs(XYZ(2)).le.abs(XYZ(1))) then
      q = LPs%pi12 * XYZ(2)/XYZ(1)
      c = cos(q)
      s = sin(q)
      q = LPs%prek * XYZ(1) / sqrt(LPs%r2-c)
      T1 = (LPs%r2*c - 1.0) * q
      T2 = LPs%r2 * s * q
    else
      q = LPs%pi12 * XYZ(1)/XYZ(2)
      c = cos(q)
      s = sin(q)
      q = LPs%prek * XYZ(2) / sqrt(LPs%r2-c)
      T1 = LPs%r2 * s * q
      T2 = (LPs%r2*c - 1.0) * q
    end if

! transform to sphere grid (inverse Lambert)
! [note that there is no need to worry about dividing by zero, since XYZ(3) can not become zero]
    c = T1**2+T2**2
    s = LPs%Pi * c/(24.0*XYZ(3)**2)
    c = LPs%sPi * c / LPs%r24 / XYZ(3)
    q = sqrt( 1.0 - s )
    LamXYZ = (/ T1 * q, T2 * q, sngl(LPs%pref) * XYZ(3) - c /)
  end if
end if

! reverse the coordinates back to the regular order according to the original pyramid number
select case (p)
 case (1,2)
  res = LamXYZ
 case (3,4)
  res = (/ LamXYZ(3), LamXYZ(1), LamXYZ(2) /)
 case (5,6)
  res = (/ LamXYZ(2), LamXYZ(3), LamXYZ(1) /)
end select

end function Lambert3DCubeForwardSingle

!--------------------------------------------------------------------------
recursive function Lambert3DCubeForwardDouble(self, res) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert3DCubeForwardDouble
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! map from 3D cubic grid to 3D ball, double precision

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=dbl),intent(out)              :: res(3)
   !! output coordinate triplet
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=dbl)                          :: XYZ(3), sXYZ(3), T1, T2, c, s, q, LamXYZ(3), eps
  integer(kind=irg)                       :: p

eps = 1.0D-12

ierr = 0
if (maxval(dabs(self%xyzd)).gt.(LPs%ap/2.D0+eps)) then
  res = (/ 0.D0, 0.D0, 0.D0 /)
  ierr = 1
  return
end if

if (maxval(dabs(self%xyzd)).lt.eps) then
  res = (/ 0.D0, 0.D0, 0.D0 /)
  return
end if

! determine which pyramid pair the point lies in and copy coordinates in correct order (see paper)
p = GetPyramidDouble(self%xyzd)
select case (p)
 case (1,2)
  sXYZ = self%xyzd
 case (3,4)
  sXYZ = (/ self%xyzd(2), self%xyzd(3), self%xyzd(1) /)
 case (5,6)
  sXYZ = (/ self%xyzd(3), self%xyzd(1), self%xyzd(2) /)
end select

! scale by grid parameter ratio sc
XYZ = LPs%sc * sXYZ

! transform to the sphere grid via the curved square, and intercept the zero point
if (maxval(dabs(XYZ)).eq.0.D0) then
  LamXYZ = (/ 0.D0, 0.D0, 0.D0 /)
else
! intercept all the points along the z-axis
  if (maxval(dabs(XYZ(1:2))).eq.0.D0) then
    LamXYZ = (/ 0.D0, 0.D0, LPs%pref * XYZ(3) /)
  else  ! this is a general grid point
    if (dabs(XYZ(2)).le.dabs(XYZ(1))) then
      q = LPs%pi12 * XYZ(2)/XYZ(1)
      c = dcos(q)
      s = dsin(q)
      q = LPs%prek * XYZ(1) / dsqrt(LPs%r2-c)
      T1 = (LPs%r2*c - 1.D0) * q
      T2 = LPs%r2 * s * q
    else
      q = LPs%pi12 * XYZ(1)/XYZ(2)
      c = dcos(q)
      s = dsin(q)
      q = LPs%prek * XYZ(2) / dsqrt(LPs%r2-c)
      T1 = LPs%r2 * s * q
      T2 = (LPs%r2*c - 1.D0) * q
    end if

! transform to sphere grid (inverse Lambert)
! [note that there is no need to worry about dividing by zero, since XYZ(3) can not become zero]
    c = T1**2+T2**2
    s = LPs%Pi * c/(24.D0*XYZ(3)**2)
    c = LPs%sPi * c / LPs%r24 / XYZ(3)
    q = dsqrt( 1.0 - s )
    LamXYZ = (/ T1 * q, T2 * q, LPs%pref * XYZ(3) - c /)
  end if
end if

! reverse the coordinates back to the regular order according to the original pyramid number
select case (p)
 case (1,2)
  res = LamXYZ
 case (3,4)
  res = (/ LamXYZ(3), LamXYZ(1), LamXYZ(2) /)
 case (5,6)
  res = (/ LamXYZ(2), LamXYZ(3), LamXYZ(1) /)
end select

end function Lambert3DCubeForwardDouble

!--------------------------------------------------------------------------
recursive function Lambert3DCubeInverseSingle(self, res) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert3DCubeInverseSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! map from 3D ball to 3D cubic grid, single/double precision

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=sgl),intent(out)              :: res(3)
   !! output coordinate triplet
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=sgl)                          :: rs, xyz3(3), xyz2(3), qxy, q2xy, sq2xy, q, ac, T1inv, T2inv, &
                                             xyz1(3), sx, sy, qx2y, sqx2y, tt
  integer(kind=irg)                       :: p

ierr = 0

rs = sqrt(sum(self%xyz*self%xyz))
if (rs.gt.LPs%R1) then
  res = (/ 0.0,0.0,0.0 /)
  ierr = 1
  return
end if

if (maxval(abs(self%xyz)).eq.0.0) then
  res = (/ 0.0,0.0,0.0 /)
  ierr = 2
  return
end if

! determine pyramid
p = GetPyramidSingle(self%xyz)
select case (p)
 case (1,2)
  xyz3 = self%xyz
 case (3,4)
  xyz3 = (/ self%xyz(2), self%xyz(3), self%xyz(1) /)
 case (5,6)
  xyz3 = (/ self%xyz(3), self%xyz(1), self%xyz(2) /)
end select

! inverse M_3
q = sqrt( 2.0*rs/(rs+abs(xyz3(3))) )
xyz2 = (/ xyz3(1) * q, xyz3(2) * q, (abs(xyz3(3))/xyz3(3)) * rs / sngl(LPs%pref) /)

! inverse M_2
qxy = xyz2(1)*xyz2(1)+xyz2(2)*xyz2(2)
sx = 1.0
if (xyz2(1).ne.0.0)  sx = abs(xyz2(1))/xyz2(1)
sy = 1.0
if (xyz2(2).ne.0.0)  sy = abs(xyz2(2))/xyz2(2)

if (qxy.ne.0.0) then
 if (abs(xyz2(2)).le.abs(xyz2(1))) then
  q2xy = qxy + xyz2(1)*xyz2(1)
  sq2xy = sqrt(q2xy)
  q = (LPs%beta/LPs%r2/LPs%R1) * sqrt(q2xy*qxy/(q2xy-abs(xyz2(1))*sq2xy))
  tt = (xyz2(2)*xyz2(2)+abs(xyz2(1))*sq2xy)/LPs%r2/qxy
  if (tt.gt.1.0) tt = 1.0
  if (tt.lt.-1.0) tt = -1.0
  ac = acos(tt)
  T1inv = q * sx
  T2inv = q * sy * ac/LPs%pi12
 else
  qx2y = qxy + xyz2(2)*xyz2(2)
  sqx2y = sqrt(qx2y)
  q = (LPs%beta/LPs%r2/LPs%R1) * sqrt(qx2y*qxy/(qx2y-abs(xyz2(2))*sqx2y))
  tt = (xyz2(1)*xyz2(1)+abs(xyz2(2))*sqx2y)/LPs%r2/qxy
  if (tt.gt.1.0) tt = 1.0
  if (tt.lt.-1.0) tt = -1.0
  ac = acos(tt)
  T1inv = q * sx * ac/LPs%pi12
  T2inv = q * sy
 end if
else
 T1inv = 0.0
 T2inv = 0.0
end if

xyz1 = (/ T1inv, T2inv, xyz2(3) /)

! inverse M_1
xyz1 = xyz1 / LPs%sc

! reverst the coordinates back to the regular order according to the original pyramid number
select case (p)
 case (1,2)
  res = xyz1
 case (3,4)
  res = (/ xyz1(3), xyz1(1), xyz1(2) /)
 case (5,6)
  res = (/ xyz1(2), xyz1(3), xyz1(1) /)
end select

end function Lambert3DCubeInverseSingle

!--------------------------------------------------------------------------
recursive function Lambert3DCubeInverseDouble(self, res) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert3DCubeInverseDouble
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! map from 3D ball to 3D cubic grid, double precision

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=dbl),intent(out)              :: res(3)
   !! output coordinate triplet
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=dbl)                          :: rs, xyz3(3), xyz2(3), qxy, q2xy, sq2xy, q, ac, T1inv, T2inv, &
                                             xyz1(3), sx, sy, qx2y, sqx2y, tt
  integer(kind=irg)                       :: p

ierr = 0

rs = dsqrt(sum(self%xyzd*self%xyzd))
if (rs.gt.LPs%R1) then
  res = (/ 0.D0,0.D0,0.D0 /)
  ierr = 1
  return
end if

if (maxval(dabs(self%xyzd)).eq.0.D0) then
  res = (/ 0.D0,0.D0,0.D0 /)
  return
end if

! determine pyramid
p = GetPyramidDouble(self%xyzd)
select case (p)
 case (1,2)
  xyz3 = self%xyzd
 case (3,4)
  xyz3 = (/ self%xyzd(2), self%xyzd(3), self%xyzd(1) /)
 case (5,6)
  xyz3 = (/ self%xyzd(3), self%xyzd(1), self%xyzd(2) /)
end select

! inverse M_3
q = dsqrt( 2.D0*rs/(rs+dabs(xyz3(3))) )
xyz2 = (/ xyz3(1) * q, xyz3(2) * q, (dabs(xyz3(3))/xyz3(3)) * rs / LPs%pref /)

! inverse M_2
qxy = xyz2(1)*xyz2(1)+xyz2(2)*xyz2(2)
sx = 1.D0
if (xyz2(1).ne.0.D0)  sx = dabs(xyz2(1))/xyz2(1)
sy = 1.D0
if (xyz2(2).ne.0.D0)  sy = dabs(xyz2(2))/xyz2(2)

if (qxy.ne.0.D0) then
 if (dabs(xyz2(2)).le.dabs(xyz2(1))) then
  q2xy = qxy + xyz2(1)*xyz2(1)
  sq2xy = dsqrt(q2xy)
  q = (LPs%beta/LPs%r2/LPs%R1) * dsqrt(q2xy*qxy/(q2xy-dabs(xyz2(1))*sq2xy))
  tt = (xyz2(2)*xyz2(2)+dabs(xyz2(1))*sq2xy)/LPs%r2/qxy
  if (tt.gt.1.D0) tt = 1.D0
  if (tt.lt.-1.D0) tt = -1.D0
  ac = dacos(tt)
  T1inv = q * sx
  T2inv = q * sy * ac/LPs%pi12
 else
  qx2y = qxy + xyz2(2)*xyz2(2)
  sqx2y = dsqrt(qx2y)
  q = (LPs%beta/LPs%r2/LPs%R1) * dsqrt(qx2y*qxy/(qx2y-dabs(xyz2(2))*sqx2y))
  tt = (xyz2(1)*xyz2(1)+dabs(xyz2(2))*sqx2y)/LPs%r2/qxy
  if (tt.gt.1.D0) tt = 1.D0
  if (tt.lt.-1.D0) tt = -1.D0
  ac = dacos(tt)
  T1inv = q * sx * ac/LPs%pi12
  T2inv = q * sy
 end if
else
  T1inv = 0.D0
  T2inv = 0.D0
end if
xyz1 = (/ T1inv, T2inv, xyz2(3) /)

! inverse M_1
xyz1 = xyz1 / LPs%sc

! reverse the coordinates back to the regular order according to the original pyramid number
select case (p)
 case (1,2)
  res = xyz1
 case (3,4)
  res = (/ xyz1(3), xyz1(1), xyz1(2) /)
 case (5,6)
  res = (/ xyz1(2), xyz1(3), xyz1(1) /)
end select

end function Lambert3DCubeInverseDouble

!--------------------------------------------------------------------------
recursive function GetPyramidSingle(xyz) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: GetPyramidSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! determine to which pyramid a point in a cubic grid belongs, single precision
  !! private function

IMPLICIT NONE

real(kind=sgl),INTENT(IN)       :: xyz(3)
integer(kind=irg)               :: res, p

logical                         :: next

next = .TRUE.
if ((abs(xyz(1)).le.xyz(3)).and.(abs(xyz(2)).le.xyz(3))) then
  p = 1                         ! pyramid 1
  next = .FALSE.
end if
if (next) then
 if ((abs(xyz(1)).le.-xyz(3)).and.(abs(xyz(2)).le.-xyz(3))) then
  p = 2                         ! pyramid 2
  next = .FALSE.
 end if
end if

if (next) then
 if ((abs(xyz(3)).le.xyz(1)).and.(abs(xyz(2)).le.xyz(1))) then
  p = 3                         ! pyramid 3
  next = .FALSE.
 end if
end if
if (next) then
 if ((abs(xyz(3)).le.-xyz(1)).and.(abs(xyz(2)).le.-xyz(1))) then
  p = 4                         ! pyramid 4
  next = .FALSE.
 end if
end if

if (next) then
 if ((abs(xyz(1)).le.xyz(2)).and.(abs(xyz(3)).le.xyz(2))) then
  p = 5                         ! pyramid 5
  next = .FALSE.
 end if
end if
if (next) then
 if ((abs(xyz(1)).le.-xyz(2)).and.(abs(xyz(3)).le.-xyz(2))) then
  p = 6                         ! pyramid 6
  next = .FALSE.
 end if
end if
res = p

end function GetPyramidSingle

!--------------------------------------------------------------------------
recursive function GetPyramidDouble(xyz) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: GetPyramidDouble
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! determine to which pyramid a point in a cubic grid belongs, double precision
  !! private function

IMPLICIT NONE

real(kind=dbl),INTENT(IN)       :: xyz(3)
integer(kind=irg)               :: res, p

logical                         :: next

next = .TRUE.
if ((dabs(xyz(1)).le.xyz(3)).and.(dabs(xyz(2)).le.xyz(3))) then
  p = 1                         ! pyramid 1
  next = .FALSE.
end if
if (next) then
 if ((dabs(xyz(1)).le.-xyz(3)).and.(dabs(xyz(2)).le.-xyz(3))) then
  p = 2                         ! pyramid 2
  next = .FALSE.
 end if
end if

if (next) then
 if ((dabs(xyz(3)).le.xyz(1)).and.(dabs(xyz(2)).le.xyz(1))) then
  p = 3                         ! pyramid 3
  next = .FALSE.
 end if
end if
if (next) then
 if ((dabs(xyz(3)).le.-xyz(1)).and.(dabs(xyz(2)).le.-xyz(1))) then
  p = 4                         ! pyramid 4
  next = .FALSE.
 end if
end if

if (next) then
 if ((dabs(xyz(1)).le.xyz(2)).and.(dabs(xyz(3)).le.xyz(2))) then
  p = 5                         ! pyramid 5
  next = .FALSE.
 end if
end if
if (next) then
 if ((dabs(xyz(1)).le.-xyz(2)).and.(dabs(xyz(3)).le.-xyz(2))) then
  p = 6                         ! pyramid 6
  next = .FALSE.
 end if
end if

res = p

end function GetPyramidDouble

!---------------
! and finally, here are the 3D cube and ball to 4D unit quaternion sphere mappings
!---------------

!--------------------------------------------------------------------------
recursive function Lambert3DBallToQuaternion(self, res) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert3DBallToQuaternion
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! map from 3D ball to unit quaternion sphere, single/double precision

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  type(Quaternion_T),intent(out)          :: res
   !! output coordinate triplet
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=sgl)                          :: q, x(21), ft, t
  real(kind=dbl)                          :: qd, xd(21), ftd, td
  integer(kind=irg)                       :: j

ierr = 0

if (self%s.eq.'s') then

  q = sqrt(sum(self%xyz**2))
  ! is the input point inside the ball that maps onto the unit quaternion sphere ?
  if (q.gt.LPs%R1) then
    res = Quaternion_T( q = (/ 0.0,0.0,0.0,0.0 /) )
    ierr = 1
    return
  end if

  if (maxval(abs(self%xyz)).eq.0.0) then
    res = Quaternion_T( q = (/ 1.0,0.0,0.0,0.0 /) )
  else
  ! get the value of t
    x(1) = 1.0
    x(2) = q**2
    do j=3,21
      x(j) = x(j-1) * x(1)
    end do
    t = sum( x * LPs%tfit )

  ! and get f(t)
    q = sqrt(1.0-t**2)
    ft = (1.5*(acos(t)-t*q))**(1.0/3.0) / q

  ! and finally create the quaternion
    res = Quaternion_T( q = (/ t, self%xyz(1)/ft, self%xyz(2)/ft, self%xyz(3)/ft /) )
  end if

else
  qd = dsqrt(sum(self%xyzd**2))
  ! is the input point inside the ball that maps onto the unit quaternion sphere ?
  if (qd.gt.LPs%R1) then
    res = Quaternion_T( qd = (/ 0.D0,0.D0,0.D0,0.D0 /) )
    ierr = 1
    return
  end if

  if (maxval(dabs(self%xyzd)).eq.0.D0) then
    res = Quaternion_T( qd = (/ 1.D0,0.D0,0.D0,0.D0 /) )
  else
  ! get the value of t
    xd(1) = 1.D0
    xd(2) = q**2
    do j=3,21
      xd(j) = xd(j-1) * xd(1)
    end do
    td = sum( xd * LPs%tfit )

  ! and get f(t)
    qd = dsqrt(1.D0-td**2)
    ftd = (1.5D0*(dacos(td)-td*qd))**(1.D0/3.D0) / qd

  ! and finally create the quaternion
    res = Quaternion_T( qd = (/ td, self%xyzd(1)/ftd, self%xyzd(2)/ftd, self%xyzd(3)/ftd /) )
  end if

end if

end function Lambert3DBallToQuaternion

!--------------------------------------------------------------------------
recursive function Lambert3DCubeToQuaternion(self, res) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert3DCubeToQuaternion
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! map from 3D cube to unit quaternion sphere, single/double precision

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  type(Quaternion_T),intent(out)          :: res
   !! output coordinate triplet
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=sgl)                          :: q(3)
  real(kind=dbl)                          :: qd(3)
  type(Lambert_T)                         :: s

  if (self%s.eq.'s') then
    ierr = Lambert3DCubeForwardSingle(self, q)
    s = Lambert_T( xyz = q )
    ierr = Lambert3DBallToQuaternion(s,res)
  else
    ierr = Lambert3DCubeForwardDouble(self, qd)
    s = Lambert_T( xyzd = qd )
    ierr = Lambert3DBallToQuaternion(s,res)
 end if

end function Lambert3DCubeToQuaternion

!---------
! the next couple of routines implement the forward and inverse stereographic projection
!---------

!--------------------------------------------------------------------------
recursive function StereoGraphicForwardSingle(self, res, Radius) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: StereoGraphicForwardSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! Forward stereographic projection (from unit sphere to plane), single precision

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=sgl),intent(out)              :: res(2)
   !! output coordinate triplet
  real(kind=sgl),INTENT(IN),OPTIONAL      :: Radius
   !! optional sphere radius
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=sgl),parameter                :: eps = 1.E-7

! input point must be on unit sphere
ierr = 0
if ( abs(1.0-sum(self%xyz**2)).gt.eps) then
  res = (/ 0.0, 0.0 /)
  ierr = 1
  return
end if

! projection
res = (/ self%xyz(1)/(1.0+self%xyz(3)), self%xyz(2)/(1.0+self%xyz(3)) /)

! scale if necessary
if (present(Radius)) res = Radius * res

end function StereoGraphicForwardSingle

!--------------------------------------------------------------------------
recursive function StereoGraphicForwardDouble(self, res, Radius) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: StereoGraphicForwardDouble

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=dbl),intent(out)              :: res(2)
   !! output coordinate triplet
  real(kind=dbl),INTENT(IN),OPTIONAL      :: Radius
   !! optional sphere radius
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=dbl),parameter                :: eps = 1.E-12

! input point must be on unit sphere
ierr = 0
if ( dabs(1.D0-sum(self%xyzd**2)).gt.eps) then
  res = (/ 0.D0, 0.D0 /)
  ierr = 1
  return
end if

! projection
res = (/ self%xyzd(1)/(1.D0+self%xyzd(3)), self%xyzd(2)/(1.D0+self%xyzd(3)) /)

! scale if necessary
if (present(Radius)) res = Radius * res

end function StereoGraphicForwardDouble

!--------------------------------------------------------------------------
recursive function StereoGraphicInverseSingle(self, res, Radius, quat) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: StereoGraphicInverseSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! Inverse stereographic projection (from plane to upper unit hemisphere), single precision

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=sgl),intent(out)              :: res(3)
   !! output coordinate triplet
  real(kind=sgl),INTENT(IN)               :: Radius
   !! projection sphere radius
  type(Quaternion_T),INTENT(IN),OPTIONAL  :: quat
   !! optional rotation quaternion
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=sgl)                          :: q, qq
  logical                                                 :: torot

torot = .FALSE.
if(present(quat)) torot = .TRUE.

ierr = 0
if (maxval(abs(self%xy)).eq.0.0) then
  res = (/ 0.0, 0.0, 1.0 /)
  if (torot.eqv..TRUE.) res = quat%quat_Lp(res)
else
  qq = sum(self%xy**2)
  if (qq.gt.Radius**2) then
    res = (/ 0.0, 0.0, 0.0 /)
    ierr = 1
  else
    q = 1.0/(Radius**2 + qq)
    res = (/ 2.0*self%xy(1), 2.0*self%xy(2), 1.0-qq /)
    res = res * q
    if (torot.eqv..TRUE.) res = quat%quat_Lp(res)
  end if
end if

end function StereoGraphicInverseSingle

!--------------------------------------------------------------------------
recursive function StereoGraphicInverseDouble(self, res, Radius, quat) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: StereoGraphicInverseDouble
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! Inverse stereographic projection (from plane to upper unit hemisphere), double precision

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=dbl),intent(out)              :: res(3)
   !! output coordinate triplet
  real(kind=dbl),INTENT(IN)               :: Radius
   !! projection sphere radius
  type(Quaternion_T),INTENT(IN),OPTIONAL  :: quat
   !! optional rotation quaternion
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=dbl)                          :: q, qq
  logical                                 :: torot

torot = .FALSE.
if(present(quat)) torot = .TRUE.

ierr = 0
if (maxval(dabs(self%xyd)).eq.0.D0) then
  res = (/ 0.D0, 0.D0, 1.D0 /)
  if (torot.eqv..TRUE.) res = quat%quat_Lp(res)
else
  qq = sum(self%xyd**2)
  if (qq.gt.Radius**2) then
    res = (/ 0.D0, 0.D0, 0.D0 /)
    ierr = 1
  else
    q = 1.D0/(Radius**2 + qq)
    res = (/ 2.D0*self%xyd(1), 2.D0*self%xyd(2), 1.D0-qq /)
    res = res * q
    if (torot.eqv..TRUE.) res = quat%quat_Lp(res)
  end if
end if

end function StereoGraphicInverseDouble

!---------
! the next couple of routines implement the standard forward and inverse Lambert projection
!---------

!--------------------------------------------------------------------------
recursive function LambertForwardSingle(self, res, Radius) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: LambertForwardSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! Forward Lambert projection (from sphere South pole to plane), single precision

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=sgl),intent(out)              :: res(2)
   !! output coordinate triplet
  real(kind=sgl),INTENT(IN),OPTIONAL      :: Radius
   !! projection sphere radius
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=sgl)                          :: q
  real(kind=sgl),parameter                :: eps = 1.E-7

! input point must be on sphere with radius Radius
ierr = 0
if ( abs(Radius**2-sum(self%xyz**2)).gt.eps) then
  res = (/ 0.0, 0.0 /)
  ierr = 1
  return
end if

! if the point is the North pole, then we have a degenerate projection
! onto the entire circle of radius 2*Radius; we'll signal that case
! with ierr=2
if (self%xyz(3).eq.Radius) then
  res = (/ 0.0, 0.0 /)
  ierr = 2
  return
end if

! otherwise we apply the forward projection
q = sqrt(2.0*Radius/(Radius - self%xyz(3)))
res = (/ q * self%xyz(1), q * self%xyz(2) /)

end function LambertForwardSingle

!--------------------------------------------------------------------------
recursive function LambertForwardDouble(self, res, Radius) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT ::  LambertForwardDouble
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! Forward Lambert projection (from sphere South pole to plane), double precision

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=dbl),intent(out)              :: res(2)
   !! output coordinate triplet
  real(kind=dbl),INTENT(IN),OPTIONAL      :: Radius
   !! projection sphere radius
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=dbl)                          :: q
  real(kind=dbl),parameter                :: eps = 1.E-12

! input point must be on sphere with radius Radius
ierr = 0
if ( dabs(Radius**2-sum(self%xyzd**2)).gt.eps) then
  res = (/ 0.D0, 0.D0 /)
  ierr = 1
  return
end if

! if the point is the North pole, then we have a degenerate projection
! onto the entire circle of radius 2*Radius; we'll signal that case
! with ierr=2
if (self%xyzd(3).eq.Radius) then
  res = (/ 0.D0, 0.D0 /)
  ierr = 2
  return
end if

! otherwise we'll apply the forward projection
q = dsqrt(2.D0*Radius/(Radius - self%xyzd(3)))
res = (/ q * self%xyzd(1), q * self%xyzd(2) /)

end function LambertForwardDouble

!--------------------------------------------------------------------------
recursive function LambertInverseSingle(self, res, Radius) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: LambertInverseSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! Inverse Lambert projection (from plane to sphere tangent at South pole), single precision

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=sgl),intent(out)              :: res(3)
   !! output coordinate triplet
  real(kind=sgl),INTENT(IN),OPTIONAL      :: Radius
   !! projection sphere radius
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=sgl)                          :: q, tr

! input point must be inside circle with radius 2*Radius
ierr = 0
q = sum(self%xy**2)
tr = 4.0*Radius**2
if (q.gt.tr) then
  res = (/ 0.0, 0.0, 0.0 /)
  ierr = 1
  return
end if

! projection (automatically takes care of projection from degenerate outer circle)
tr = sqrt(1.0 - q/tr)
res = (/ tr * self%xy(1), tr * self%xy(2), -Radius + q/(2.0*Radius) /)

end function LambertInverseSingle

!--------------------------------------------------------------------------
recursive function LambertInverseDouble(self, res, Radius) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: LambertInverseDouble
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! Inverse Lambert projection (from plane to sphere tangent at South pole), double precision

IMPLICIT NONE

  class(Lambert_T),intent(in)             :: self
   !! input Lambert class
  real(kind=dbl),intent(out)              :: res(3)
   !! output coordinate triplet
  real(kind=dbl),INTENT(IN),OPTIONAL      :: Radius
   !! projection sphere radius
  integer(kind=irg)                       :: ierr
  !f2py intent(in,out) ::  ierr

  real(kind=dbl)                          :: q, tr

! input point must be inside circle with radius 2*Radius
ierr = 0
q = sum(self%xyd**2)
tr = 4.D0*Radius**2
if (q.gt.tr) then
  res = (/ 0.D0, 0.D0, 0.D0 /)
  ierr = 1
  return
end if

! projection (automatically takes care of projection from degenerate outer circle)
tr = dsqrt(1.D0 - q/tr)
res = (/ tr * self%xyd(1), tr * self%xyd(2), -Radius + q/(2.D0*Radius) /)

end function LambertInverseDouble

!--------------------------------------------------------------------------
recursive subroutine LambertgetInterpolationSingle(dc, scl, npx, npy, nix, niy, nixp, niyp, dx, dy, dxm, dym, swap)
!DEC$ ATTRIBUTES DLLEXPORT :: LambertgetInterpolationSingle

  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! take direction cosines and return all parameters for square Lambert interpolation, single precision
  !! this piece of code replaces code that occurred many times in various programs

use mod_io

real(kind=sgl),INTENT(IN)           :: dc(3)
 !! direction cosines
real(kind=sgl),INTENT(IN)           :: scl
 !! scale parameter for square Lambert projection
integer(kind=irg),INTENT(IN)        :: npx
 !! number of pixels along square semi-edge
integer(kind=irg),INTENT(IN)        :: npy
 !! should be the same as npx
integer(kind=irg),INTENT(OUT)       :: nix
 !! x-coordinate of point
integer(kind=irg),INTENT(OUT)       :: niy
 !! y-coordinate of point
integer(kind=irg),INTENT(OUT)       :: nixp
 !! x of neighboring point
integer(kind=irg),INTENT(OUT)       :: niyp
 !! y of neighboring point
real(kind=sgl),INTENT(OUT)          :: dx
 !! x interpolation weight factors
real(kind=sgl),INTENT(OUT)          :: dy
 !! y interpolation weight factors
real(kind=sgl),INTENT(OUT)          :: dxm
 !! 1-x interpolation weight factors
real(kind=sgl),INTENT(OUT)          :: dym
 !! 1-y interpolation weight factors
logical,INTENT(IN),OPTIONAL         :: swap
 !! sometimes we need to swap the x and y coordinates (OPTIONAL)

type(Lambert_T)                     :: Lambert
type(IO_T)                          :: Message
real(kind=sgl)                      :: xy(2), x, io_real(3)
integer(kind=irg)                   :: istat

! Lambert sphere to square transformation
Lambert = Lambert_T( xyz=dc )
istat = Lambert2DSquareInverseSingle(Lambert, xy)
if (istat .ne. 0) then
  io_real(1:3) = dc
  call Message%WriteValue('input direction cosines : ', io_real, 3)
  io_real(1) = scl
  call Message%WriteValue('input scale factor      : ', io_real, 1)
  call Message%printWarning( 'LambertgetInterpolationSingle', (/'Something went wrong during interpolation...'/) )
end if
xy = scl * xy

if (present(swap)) then
  if (swap.eqv..TRUE.) then
    x = xy(1)
    xy(1) = xy(2)
    xy(2) = -x
  end if
end if

! four-point interpolation (bi-quadratic)
nix = int(npx+xy(1))-npx
niy = int(npy+xy(2))-npy
nixp = nix+1
niyp = niy+1
if (nixp.gt.npx) nixp = nix
if (niyp.gt.npy) niyp = niy
if (nix.lt.-npx) nix = nixp
if (niy.lt.-npy) niy = niyp
dx = xy(1)-nix
dy = xy(2)-niy
dxm = 1.0-dx
dym = 1.0-dy

end subroutine LambertgetInterpolationSingle

!--------------------------------------------------------------------------
recursive subroutine LambertgetInterpolationDouble(dc, scl, npx, npy, nix, niy, nixp, niyp, dx, dy, dxm, dym, swap)
!DEC$ ATTRIBUTES DLLEXPORT :: LambertgetInterpolationDouble
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! take direction cosines and return all parameters for square Lambert interpolation, double precision
  !! this piece of code replaces code that occurred many times in various programs

use mod_io

real(kind=dbl),INTENT(IN)           :: dc(3)
 !! direction cosines
real(kind=dbl),INTENT(IN)           :: scl
 !! scale parameter for square Lambert projection
integer(kind=irg),INTENT(IN)        :: npx
 !! number of pixels along square semi-edge
integer(kind=irg),INTENT(IN)        :: npy
 !! should be the same as npx
integer(kind=irg),INTENT(OUT)       :: nix
 !! x-coordinate of point
integer(kind=irg),INTENT(OUT)       :: niy
 !! y-coordinate of point
integer(kind=irg),INTENT(OUT)       :: nixp
 !! x of neighboring point
integer(kind=irg),INTENT(OUT)       :: niyp
 !! y of neighboring point
real(kind=dbl),INTENT(OUT)          :: dx
 !! x interpolation weight factors
real(kind=dbl),INTENT(OUT)          :: dy
 !! y interpolation weight factors
real(kind=dbl),INTENT(OUT)          :: dxm
 !! 1-x interpolation weight factors
real(kind=dbl),INTENT(OUT)          :: dym
 !! 1-y interpolation weight factors
logical,INTENT(IN),OPTIONAL         :: swap
 !! sometimes we need to swap the x and y coordinates (OPTIONAL)

type(Lambert_T)                     :: Lambert
type(IO_T)                          :: Message
real(kind=dbl)                      :: xy(2), x, io_real(3)
integer(kind=irg)                   :: istat

! Lambert sphere to square transformation
Lambert = Lambert_T( xyzd=dc )
istat = Lambert2DSquareInverseDouble(Lambert, xy)
if (istat .ne. 0) then
  io_real(1:3) = dc
  call Message%WriteValue('input direction cosines : ', io_real, 3)
  io_real(1) = scl
  call Message%WriteValue('input scale factor      : ', io_real, 1)
  call Message%printWarning( 'LambertgetInterpolationDouble', (/'Something went wrong during interpolation...'/) )
end if
xy = scl * xy

if (present(swap)) then
  if (swap.eqv..TRUE.) then
    x = xy(1)
    xy(1) = xy(2)
    xy(2) = -x
  end if
end if

! four-point interpolation (bi-quadratic)
nix = int(npx+xy(1))-npx
niy = int(npy+xy(2))-npy
nixp = nix+1
niyp = niy+1
if (nixp.gt.npx) nixp = nix
if (niyp.gt.npy) niyp = niy
if (nix.lt.-npx) nix = nixp
if (niy.lt.-npy) niy = niyp
dx = xy(1)-nix
dy = xy(2)-niy
dxm = 1.D0-dx
dym = 1.D0-dy

end subroutine LambertgetInterpolationDouble

!--------------------------------------------------------------------------
recursive function InterpolationLambert2DSingle(dc, m, npx) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: InterpolationLambert2DSingle

  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! perform a Lambert interpolation

IMPLICIT NONE

real(kind=sgl),INTENT(INOUT)            :: dc(3)
 !! direction cosines
!f2py intent(in,out) ::  dc
integer(kind=irg),INTENT(IN)            :: npx
 !! number of pixels along master patter nsemi-edge
real(kind=sgl),INTENT(IN)               :: m(-npx:npx,-npx:npx)
 !! master pattern
real(kind=sgl)                          :: res
 !! interpolated intensity

integer(kind=irg)                       :: nix, niy, nixp, niyp, istat
real(kind=sgl)                          :: xy(2), dx, dy, dxm, dym, scl

scl = float(npx)

if (dc(3).lt.0.0) dc = -dc

! convert direction cosines to lambert projections
call LambertgetInterpolation(dc, scl, npx, npx, nix, niy, nixp, niyp, dx, dy, dxm, dym)

res = m(nix,niy)*dxm*dym + m(nixp,niy)*dx*dym + m(nix,niyp)*dxm*dy + m(nixp,niyp)*dx*dy

end function InterpolationLambert2DSingle

!--------------------------------------------------------------------------
recursive function InterpolationLambert2DDouble(dc, m, npx) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: InterpolationLambert2DDouble

  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! perform a Lambert interpolation

IMPLICIT NONE

real(kind=dbl),INTENT(INOUT)            :: dc(3)
 !! direction cosines
!f2py intent(in,out) ::  dc
integer(kind=irg),INTENT(IN)            :: npx
 !! number of pixels along master pattern semi-edge
real(kind=dbl),INTENT(IN)               :: m(-npx:npx,-npx:npx)
 !! master pattern
real(kind=dbl)                          :: res
 !! interpolated intensity

integer(kind=irg)                       :: nix, niy, nixp, niyp, istat
real(kind=dbl)                          :: xy(2), dx, dy, dxm, dym, scl

scl = dble(npx)

if (dc(3).lt.0.0) dc = -dc

! convert direction cosines to lambert projections
call LambertgetInterpolation(dc, scl, npx, npx, nix, niy, nixp, niyp, dx, dy, dxm, dym)

res = m(nix,niy)*dxm*dym + m(nixp,niy)*dx*dym + m(nix,niyp)*dxm*dy + m(nixp,niyp)*dx*dy

end function InterpolationLambert2DDouble

!--------------------------------------------------------------------------
recursive function InterpolationLambert3DSingle(dc, m, npx, nn, s) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: InterpolationLambert3DSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! perform a Lambert interpolation

IMPLICIT NONE

real(kind=dbl),INTENT(INOUT)            :: dc(3)
 !! direction cosines
!f2py intent(in,out) ::  dc
integer(kind=irg),INTENT(IN)            :: npx
 !! number of pixels along master pattern semi-edge
integer(kind=irg),INTENT(IN)            :: nn
 !! number of entries in output array
real(kind=sgl),INTENT(IN)               :: m(-npx:npx,-npx:npx, nn)
 !! master pattern
real(kind=sgl),INTENT(INOUT),OPTIONAL   :: s(nn)
real(kind=sgl)                          :: res
 !! output intensity array

integer(kind=irg)                       :: nix, niy, nixp, niyp, istat
real(kind=sgl)                          :: xy(2), dx, dy, dxm, dym, scl, locals(nn)

scl = float(npx)

if (dc(3).lt.0.0) dc = -dc

! convert direction cosines to lambert projections
call LambertgetInterpolation(sngl(dc), scl, npx, npx, nix, niy, nixp, niyp, dx, dy, dxm, dym)

locals(1:nn) = m(nix,niy,1:nn)*dxm*dym + m(nixp,niy,1:nn)*dx*dym + &
               m(nix,niyp,1:nn)*dxm*dy + m(nixp,niyp,1:nn)*dx*dy

if (present(s)) s = locals 

res = sum(locals)

end function InterpolationLambert3DSingle

!--------------------------------------------------------------------------
recursive function InterpolationLambert3DInteger(dc, m, npx, nn) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: InterpolationLambert3DInteger
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! perform a Lambert interpolation

IMPLICIT NONE

real(kind=dbl),INTENT(INOUT)            :: dc(3)
 !! direction cosines
!f2py intent(in,out) ::  dc
integer(kind=irg),INTENT(IN)            :: npx
 !! number of pixels along master pattern semi-edge
integer(kind=irg),INTENT(IN)            :: nn
 !! number of entries in output array
integer(kind=irg),INTENT(IN)            :: m(nn,-npx:npx,-npx:npx)
 !! master pattern
real(kind=sgl)                          :: res(nn)
 !! output array

integer(kind=irg)                       :: nix, niy, nixp, niyp, istat
real(kind=sgl)                          :: xy(2), dx, dy, dxm, dym, scl

scl = float(npx)

if (dc(3).lt.0.0) dc = -dc

! convert direction cosines to lambert projections
call LambertgetInterpolation(sngl(dc), scl, npx, npx, nix, niy, nixp, niyp, dx, dy, dxm, dym)

res(1:nn) = m(1:nn,nix,niy)*dxm*dym + m(1:nn,nixp,niy)*dx*dym + &
            m(1:nn,nix,niyp)*dxm*dy + m(1:nn,nixp,niyp)*dx*dy

end function InterpolationLambert3DInteger

!--------------------------------------------------------------------------
recursive function InterpolationLambert4DSingle(dc, m, npx, nn) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: InterpolationLambert4DSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! perform a Lambert interpolation

IMPLICIT NONE

real(kind=dbl),INTENT(INOUT)            :: dc(3)
 !! direction cosines
!f2py intent(in,out) ::  dc
integer(kind=irg),INTENT(IN)            :: npx
 !! number of pixels along master pattern semi-edge
integer(kind=irg),INTENT(IN)            :: nn
 !! number oentries in output array
real(kind=sgl),INTENT(IN)               :: m(-npx:npx,-npx:npx, 1, nn)
 !! master pattern
real(kind=sgl)                          :: res
 !! output array

integer(kind=irg)                       :: nix, niy, nixp, niyp, istat
real(kind=sgl)                          :: xy(2), dx, dy, dxm, dym, scl, resarray(nn)

scl = float(npx)

if (dc(3).lt.0.0) dc = -dc

! convert direction cosines to lambert projections
call LambertgetInterpolation(sngl(dc), scl, npx, npx, nix, niy, nixp, niyp, dx, dy, dxm, dym)

resarray(1:nn) = m(nix,niy,1,1:nn)*dxm*dym + m(nixp,niy,1,1:nn)*dx*dym + &
                 m(nix,niyp,1,1:nn)*dxm*dy + m(nixp,niyp,1,1:nn)*dx*dy

res = sum(resarray)

end function InterpolationLambert4DSingle

!--------------------------------------------------------------------------
recursive function InterpolationLambert4DDouble4b4(dc, m, npx) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: InterpolationLambert4DDouble4b4
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! perform a Lambert interpolation

IMPLICIT NONE

real(kind=dbl),INTENT(INOUT)            :: dc(3)
 !! direction cosines
!f2py intent(in,out) ::  dc
integer(kind=irg),INTENT(IN)            :: npx
 !! number of pixels along master pattern semi edge
real(kind=dbl),INTENT(IN)               :: m(4,4,-npx:npx,-npx:npx)
 !! master pattern
real(kind=dbl)                          :: res(4,4)
 !! output array

integer(kind=irg)                       :: nix, niy, nixp, niyp, istat
real(kind=dbl)                          :: xy(2), dx, dy, dxm, dym, scl

scl = dble(npx)

if (dc(3).lt.0.0) dc = -dc

! convert direction cosines to lambert projections
call LambertgetInterpolation(dc, scl, npx, npx, nix, niy, nixp, niyp, dx, dy, dxm, dym)

res(1:4,1:4) = m(1:4,1:4,nix,niy)*dxm*dym + m(1:4,1:4,nixp,niy)*dx*dym + &
               m(1:4,1:4,nix,niyp)*dxm*dy + m(1:4,1:4,nixp,niyp)*dx*dy

end function InterpolationLambert4DDouble4b4

!--------------------------------------------------------------------------
recursive subroutine Apply3DPGSymmetry_(self,cell,SG,ipx,ipy,ipz,npx,iequiv,nequiv,usehex,stereographic,cubictype)
!DEC$ ATTRIBUTES DLLEXPORT :: Apply3DPGSymmetry_
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! Apply the 3D point group symmetry to a pair of coordinates on a Lambert grid

use mod_crystallography
use mod_symmetry
use mod_io

IMPLICIT NONE

class(Lambert_T), INTENT(INOUT) :: self
type(Cell_T),INTENT(INOUT)      :: cell
type(SpaceGroup_T),INTENT(INOUT):: SG
integer(kind=irg),INTENT(IN)    :: ipx
integer(kind=irg),INTENT(IN)    :: ipy
integer(kind=irg),INTENT(IN)    :: ipz
integer(kind=irg),INTENT(IN)    :: npx
integer(kind=irg),INTENT(OUT)   :: iequiv(3,48)
integer(kind=irg),INTENT(OUT)   :: nequiv
logical,INTENT(IN),OPTIONAL     :: usehex
logical,INTENT(IN),OPTIONAL     :: stereographic
integer(kind=irg),INTENT(IN),OPTIONAL :: cubictype

type(Lambert_T)                 :: L
type(IO_T)                      :: Message
real(kind=dbl)                  :: xy(2), xyz(3), kstar(3)
real(kind=dbl),parameter        :: neps = -0.0001D0
integer(kind=irg)               :: ierr, i, ix, iy
real(kind=dbl),allocatable      :: stmp(:,:)            !< output array with equivalent vectors
integer(kind=irg)               :: n                    !< number of entries in equivalent vector array
character(1)                    :: space                !< 'd' or 'r'


! for the cubic groups, we need to apply a lower symmetry group due to the fact that we
! do not use interpolations to apply the three-fold axes.  So, for all space groups with
! number .ge. 195, we must apply a lower symmetry group.  To avoid unnecessary repetitions,
! we pre-compute the relevant symmetry operations and store them in an extra SG%SYM_extra array
! when we apply the symmetry for the first time.  Then, we test the space group number in
! the CalcStar routine to see which set of arrays to use.  Since the cubic symmetry is used a lot,
! we hard-coded the symmetry operations to make things go slightly faster.


! we have points on the square/hexagonal Lambert projection and we need to determine the
! set of equivalent points; we can use the CalcStar routine to do this, but first we
! need to convert the 2D coordinates into a 3D vector in reciprocal space.
xy = (/ dble(ipx), dble(ipy) /) / dble(npx)
L = Lambert_T( xyd = xy )
if (present(usehex)) then
  ierr = L%LambertHexToSphere(xyz)
else
  ierr = L%LambertSquareToSphere(xyz)
end if
if (ipz.lt.0) xyz(3) = -xyz(3)
! convert to reciprocal space
call cell%NormVec(xyz, 'c')
call cell%TransSpace(xyz, kstar, 'c', 'r')

iequiv = 0
nequiv = 0

! apply the 3D point group to get the complete star of kstar
if (present(cubictype)) then
  select case (cubictype)
    case(3)
        n = 2
        allocate(stmp(n,3))
        stmp(1,1:3) = kstar(1:3)
        stmp(2,1:3) = (/ -kstar(1), kstar(2), -kstar(3) /)

    case(6)
        n = 8
        allocate(stmp(n,3))
        stmp(1,1:3) = kstar(1:3)
        stmp(2,1:3) = (/ -kstar(1), -kstar(2),  kstar(3) /)
        stmp(3,1:3) = (/ -kstar(1),  kstar(2), -kstar(3) /)
        stmp(4,1:3) = (/ -kstar(1), -kstar(2), -kstar(3) /)
        stmp(5,1:3) = (/  kstar(1), -kstar(2), -kstar(3) /)
        stmp(6,1:3) = (/  kstar(1),  kstar(2), -kstar(3) /)
        stmp(7,1:3) = (/  kstar(1), -kstar(2),  kstar(3) /)
        stmp(8,1:3) = (/ -kstar(1),  kstar(2),  kstar(3) /)

    case(8)
        n = 8
        allocate(stmp(n,3))
        stmp(1,1:3) = kstar(1:3)
        stmp(2,1:3) = (/ -kstar(1), -kstar(2),  kstar(3) /)
        stmp(3,1:3) = (/ -kstar(1),  kstar(2), -kstar(3) /)
        stmp(4,1:3) = (/  kstar(1), -kstar(2), -kstar(3) /)
        stmp(5,1:3) = (/  kstar(2), -kstar(1), -kstar(3) /)
        stmp(6,1:3) = (/ -kstar(2),  kstar(1), -kstar(3) /)
        stmp(7,1:3) = (/ -kstar(2), -kstar(1),  kstar(3) /)
        stmp(8,1:3) = (/  kstar(2),  kstar(1),  kstar(3) /)

    case(9)
        n = 16
        allocate(stmp(n,3))
        stmp(1,1:3) = kstar(1:3)
        stmp(2,1:3) =  (/ -kstar(1), -kstar(2),  kstar(3) /)
        stmp(3,1:3) =  (/ -kstar(1),  kstar(2), -kstar(3) /)
        stmp(4,1:3) =  (/ -kstar(1), -kstar(2), -kstar(3) /)
        stmp(5,1:3) =  (/  kstar(1), -kstar(2), -kstar(3) /)
        stmp(6,1:3) =  (/  kstar(1),  kstar(2), -kstar(3) /)
        stmp(7,1:3) =  (/  kstar(1), -kstar(2),  kstar(3) /)
        stmp(8,1:3) =  (/ -kstar(1),  kstar(2),  kstar(3) /)
        stmp(9,1:3) =  (/ -kstar(2),  kstar(1),  kstar(3) /)
        stmp(10,1:3) = (/  kstar(2), -kstar(1),  kstar(3) /)
        stmp(11,1:3) = (/  kstar(2),  kstar(1), -kstar(3) /)
        stmp(12,1:3) = (/  kstar(2), -kstar(1), -kstar(3) /)
        stmp(13,1:3) = (/ -kstar(2), -kstar(1), -kstar(3) /)
        stmp(14,1:3) = (/ -kstar(2),  kstar(1), -kstar(3) /)
        stmp(15,1:3) = (/  kstar(2),  kstar(1),  kstar(3) /)
        stmp(16,1:3) = (/ -kstar(2), -kstar(1),  kstar(3) /)

    case default
        call Message%printError('Apply3DPGSymmetry','unknown cubictype parameter [3, 6, 8, or 9]')
  end select

else
  space = 'r'
  call SG%CalcStar(kstar,n,stmp,space)
end if

! then convert the equivalent points back into 2D Lambert coordinates
do i=1,n
  call cell%TransSpace(stmp(i,1:3), xyz, 'r', 'c')
  call cell%NormVec(xyz, 'c')
  iequiv(3,i) = 1
  if (xyz(3).lt.neps) iequiv(3,i) = -1
  if (present(stereographic)) then
! export stereographic coordinates
   if (iequiv(3,i).eq.1) then
    ix = int(dble(npx)*xyz(1)/(1.D0+xyz(3)))
    iy = int(dble(npx)*xyz(2)/(1.D0+xyz(3)))
   else
    ix = int(dble(npx)*xyz(1)/(1.D0-xyz(3)))
    iy = int(dble(npx)*xyz(2)/(1.D0-xyz(3)))
   end if
   iequiv(1,i) = ix
   iequiv(2,i) = iy
  else
    L = Lambert_T( xyzd = xyz )
    if (present(usehex)) then
      ierr = L%LambertSphereToHex(xy)
    else
      ierr = L%LambertSphereToSquare(xy)
    end if
    xy = xy * dble(npx)
    iequiv(1,i) = nint(xy(1))
    iequiv(2,i) = nint(xy(2))
  end if
end do
nequiv = n

end subroutine Apply3DPGSymmetry_

recursive subroutine sampleVMF_(self, mu, kappa, VMFscale, inten, npx, nix, niy, w, mLPNH, mLPSH, LegendreArray)
!DEC$ ATTRIBUTES DLLEXPORT :: sampleVMF_
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! sample a p=3 von Mises-Fisher distribution for a Lambert grid patch around a given direction

IMPLICIT NONE 


class(Lambert_T), INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)     :: mu(3)
real(kind=dbl),INTENT(IN)     :: kappa
real(kind=dbl),INTENT(IN)     :: VMFscale
real(kind=dbl),INTENT(IN)     :: inten
integer(kind=irg),INTENT(IN)  :: npx
integer(kind=irg),INTENT(IN)  :: nix
integer(kind=irg),INTENT(IN)  :: niy
integer(kind=irg),INTENT(IN)  :: w
real(kind=sgl),INTENT(INOUT)  :: mLPNH(-npx:npx, -npx:npx)
!f2py intent(in,out) ::  mLPNH
real(kind=sgl),INTENT(INOUT)  :: mLPSH(-npx:npx, -npx:npx)
!f2py intent(in,out) ::  mLPSH
real(kind=dbl),INTENT(IN)     :: LegendreArray(0:2*npx)

real(kind=sgl)                :: xyz(3), vmf , LegendreLattitude, p
integer(kind=irg)             :: i, j, ix, iy
logical                       :: North, xN, yN  

North = .TRUE.
if (mu(3).lt.0.0) North = .FALSE.

do i=-w, w 
  ix = nix + i 
  do j=-w, w 
    iy = niy + j  
! check the hemisphere and properly wrap where needed
    xyz = self%HemiCheck_(ix, iy, npx, North)
! correct the angle to the Legendre lattitude 
    LegendreLattitude = sngl(LegendreArray( maxval( abs( (/ ix, iy /) ) )) )
! the factor p rescales the x and y components of kstar to maintain a unit vector
    p = sqrt((1.D0-LegendreLattitude**2)/(1.D0-xyz(3)**2))
    xyz = (/ p*xyz(1), p*xyz(2), LegendreLattitude /)
! compute the VMF value
    vmf = (-1.D0 + sum(mu*xyz)) * kappa + Log(inten) + VMFscale
! put this value in the correct array location
!    if (xyz(3).ge.0.0) then
      mLPNH(ix, iy) = mLPNH(ix, iy) + exp(vmf) 
      mLPSH(-ix, -iy) = mLPSH(-ix, -iy) + exp(vmf)
!    end if
  end do 
end do 

end subroutine sampleVMF_

!--------------------------------------------------------------------------
recursive function HemiCheck_(self, ix, iy, npx, North) result(xyz) 
!DEC$ ATTRIBUTES DLLEXPORT :: HemiCheck_
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! get the xyz belonging to an integer input point that may be Northern or Southern hemishpere...
IMPLICIT NONE 

class(Lambert_T), INTENT(INOUT) :: self
integer(kind=irg),INTENT(INOUT)     :: ix
!f2py intent(in,out) ::  ix
integer(kind=irg),INTENT(INOUT)     :: iy
!f2py intent(in,out) ::  iy
integer(kind=irg),INTENT(IN)        :: npx
logical,INTENT(IN)                  :: North
real(kind=sgl)                      :: xyz(3)

type(Lambert_T)                     :: L

integer(kind=irg)                   :: ierr 

if ((abs(ix).le.npx).and.(abs(iy).le.npx)) then
! regular case
  L = Lambert_T( xy = (float( (/ ix, iy/) ) / float(npx)))
  ierr = L%LambertSquareToSphere(xyz)
else 
    if ((abs(ix).gt.npx).and.(abs(iy).gt.npx)) then 
! corner case
        if (ix.gt.0) then 
          ix = 2*npx-ix
        else
          ix = -2*npx-ix
        end if 
        if (iy.gt.0) then 
          iy = 2*npx-iy
        else
          iy = -2*npx-iy
        end if 
        L = Lambert_T( xy = (float( (/ ix, iy/) ) / float(npx)))
        ierr = L%LambertSquareToSphere(xyz)
        if (North.eqv..TRUE.) xyz(3) = -xyz(3)
     else
! remaining cases
        if ((abs(ix).gt.npx).and.(abs(iy).le.npx)) then 
! edge case
          if (ix.gt.0) then 
            ix = 2*npx-ix
          else
            ix = -2*npx-ix
          end if 
          L = Lambert_T( xy = (float( (/ ix, iy/) ) / float(npx)))
          ierr = L%LambertSquareToSphere(xyz)
          if (North.eqv..TRUE.) xyz(3) = -xyz(3)
        end if
        if ((abs(ix).le.npx).and.(abs(iy).gt.npx)) then 
! edge case
          if (iy.gt.0) then 
            iy = 2*npx-iy
          else
            iy = -2*npx-iy
          end if 
          L = Lambert_T( xy = (float( (/ ix, iy/) ) / float(npx)))
          ierr = L%LambertSquareToSphere(xyz)
          if (North.eqv..TRUE.) xyz(3) = -xyz(3)
        end if
    end if
end if
end function HemiCheck_

end module mod_Lambert
