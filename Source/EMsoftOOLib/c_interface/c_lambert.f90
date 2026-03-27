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

module c_lambert
  !! author: MDG
  !! version: 1.0
  !! date: 03/27/26
  !!
  !! C-interop wrappers for Lambert projection utilities.
  !! Provides square<->sphere and cube<->ball mappings.

use iso_c_binding
use mod_kinds
use mod_global
use mod_Lambert

IMPLICIT NONE

private

contains

!--------------------------------------------------------------------------
! Square <-> Sphere (2D Lambert projection)
!--------------------------------------------------------------------------

subroutine c_lambert_square_to_sphere(xy, xyz, ierr) &
    bind(c, name='emsoft_lambert_square_to_sphere')
  !! Map a 2D point on the square [-1,1]^2 to a point on the unit hemisphere.
  real(c_double), INTENT(IN)     :: xy(2)
  real(c_double), INTENT(OUT)    :: xyz(3)
  integer(c_int), INTENT(OUT)    :: ierr
  type(Lambert_T)                :: L
  integer(kind=irg)              :: err

  L = Lambert_T( xyd = xy )
  xyz = L%LambertSquareToSphere(err)
  ierr = err

end subroutine c_lambert_square_to_sphere

!--------------------------------------------------------------------------
subroutine c_lambert_sphere_to_square(xyz, xy, ierr) &
    bind(c, name='emsoft_lambert_sphere_to_square')
  !! Map a point on the unit hemisphere to the square [-1,1]^2.
  real(c_double), INTENT(IN)     :: xyz(3)
  real(c_double), INTENT(OUT)    :: xy(2)
  integer(c_int), INTENT(OUT)    :: ierr
  type(Lambert_T)                :: L
  integer(kind=irg)              :: err

  L = Lambert_T( xyzd = xyz )
  xy = L%LambertSphereToSquare(err)
  ierr = err

end subroutine c_lambert_sphere_to_square

!--------------------------------------------------------------------------
! Cube <-> Ball (3D Lambert projection)
!--------------------------------------------------------------------------

subroutine c_lambert_cube_to_ball(cube, ball, ierr) &
    bind(c, name='emsoft_lambert_cube_to_ball')
  !! Map a 3D point in the cube to a point in the unit ball.
  real(c_double), INTENT(IN)     :: cube(3)
  real(c_double), INTENT(OUT)    :: ball(3)
  integer(c_int), INTENT(OUT)    :: ierr
  type(Lambert_T)                :: L
  integer(kind=irg)              :: err

  L = Lambert_T( xyzd = cube )
  ball = L%LambertCubeToBall(err)
  ierr = err

end subroutine c_lambert_cube_to_ball

!--------------------------------------------------------------------------
subroutine c_lambert_ball_to_cube(ball, cube, ierr) &
    bind(c, name='emsoft_lambert_ball_to_cube')
  !! Map a 3D point in the unit ball to a point in the cube.
  real(c_double), INTENT(IN)     :: ball(3)
  real(c_double), INTENT(OUT)    :: cube(3)
  integer(c_int), INTENT(OUT)    :: ierr
  type(Lambert_T)                :: L
  integer(kind=irg)              :: err

  L = Lambert_T( xyzd = ball )
  cube = L%LambertBallToCube(err)
  ierr = err

end subroutine c_lambert_ball_to_cube

!--------------------------------------------------------------------------
! Stereographic projection
!--------------------------------------------------------------------------

subroutine c_lambert_stereo_forward(xyz, xy, ierr) &
    bind(c, name='emsoft_lambert_stereo_forward')
  !! Forward stereographic projection: 3D unit sphere -> 2D plane.
  real(c_double), INTENT(IN)     :: xyz(3)
  real(c_double), INTENT(OUT)    :: xy(2)
  integer(c_int), INTENT(OUT)    :: ierr
  type(Lambert_T)                :: L
  integer(kind=irg)              :: err

  L = Lambert_T( xyzd = xyz )
  xy = L%StereoGraphicForward(err)
  ierr = err

end subroutine c_lambert_stereo_forward

!--------------------------------------------------------------------------
subroutine c_lambert_stereo_inverse(xy, xyz, ierr) &
    bind(c, name='emsoft_lambert_stereo_inverse')
  !! Inverse stereographic projection: 2D plane -> 3D unit sphere.
  real(c_double), INTENT(IN)     :: xy(2)
  real(c_double), INTENT(OUT)    :: xyz(3)
  integer(c_int), INTENT(OUT)    :: ierr
  type(Lambert_T)                :: L
  integer(kind=irg)              :: err

  L = Lambert_T( xyd = xy )
  xyz = L%StereoGraphicInverse(err)
  ierr = err

end subroutine c_lambert_stereo_inverse

end module c_lambert
