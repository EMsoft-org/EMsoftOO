! ###################################################################
! Copyright (c) 2013-2022, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_IPFsupport
  !! author: MDG 
  !! version: 1.0 
  !! date: 09/03/21
  !!
  !! module for inverse pole figure maps

use mod_kinds
use mod_global

IMPLICIT NONE 

type, private :: sPatch
  real(kind=dbl),allocatable    :: rx(:)
  real(kind=dbl),allocatable    :: ry(:)
  real(kind=dbl),allocatable    :: center(:)
  real(kind=dbl),allocatable    :: normals(:,:)
  real(kind=dbl),allocatable    :: cutoffs(:)
  real(kind=dbl),allocatable    :: coeffs(:,:)
  real(kind=dbl),allocatable    :: cumAngles(:)
  real(kind=dbl),allocatable    :: omega(:)
end type sPatch 


! class definition
type, public :: IPFmap_T
private 
  integer(kind=irg)   :: ipf_wd
  integer(kind=irg)   :: ipf_ht
  integer(kind=irg)   :: ipf_LaueClass
  integer(kind=irg)   :: ipf_nthreads
  character(fnlen)    :: ipf_mode
  character(fnlen)    :: ipf_filename
  type(sPatch)        :: sphericalPatch

contains
private 
  procedure, pass(self) :: get_ipf_wd_
  procedure, pass(self) :: get_ipf_ht_
  procedure, pass(self) :: get_ipf_mode_
  procedure, pass(self) :: get_ipf_filename_
  procedure, pass(self) :: get_ipf_LaueClass_
  procedure, pass(self) :: get_ipf_nthreads_
  procedure, pass(self) :: set_ipf_wd_
  procedure, pass(self) :: set_ipf_ht_
  procedure, pass(self) :: set_ipf_mode_
  procedure, pass(self) :: set_ipf_filename_
  procedure, pass(self) :: set_ipf_LaueClass_
  procedure, pass(self) :: set_ipf_nthreads_
  ! procedure, pass(self) :: SphericalTriangle_
  ! procedure, pass(self) :: SphericalWedge_
  ! procedure, pass(self) :: init_sPatch_
  ! procedure, pass(self) :: toHemi_
  ! procedure, pass(self) :: toColor_ 
  ! procedure, pass(self) :: build_
  procedure, pass(self) :: get_ipf_RGB_
  procedure, pass(self) :: get_IPFMap_

  generic, public :: get_ipf_wd => get_ipf_wd_
  generic, public :: get_ipf_ht => get_ipf_ht_
  generic, public :: get_ipf_mode => get_ipf_mode_
  generic, public :: get_ipf_filename => get_ipf_filename_
  generic, public :: get_ipf_LaueClass => get_ipf_LaueClass_
  generic, public :: get_ipf_nthreads => get_ipf_nthreads_
  generic, public :: set_ipf_wd => set_ipf_wd_
  generic, public :: set_ipf_ht => set_ipf_ht_
  generic, public :: set_ipf_mode => set_ipf_mode_
  generic, public :: set_ipf_filename => set_ipf_filename_
  generic, public :: set_ipf_LaueClass => set_ipf_LaueClass_
  generic, public :: set_ipf_nthreads => set_ipf_nthreads_
  generic, public :: get_ipf_RGB => get_ipf_RGB_
  generic, public :: get_IPFMap => get_IPFMap_

end type IPFmap_T

! the constructor routine for this class 
interface IPFmap_T
  module procedure IPF_constructor
end interface IPFmap_T

contains

!--------------------------------------------------------------------------
type(IPFmap_T) function IPF_constructor( ) result(IPF)
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! constructor for the IPFmap_T Class; reads the name list 
 
IMPLICIT NONE


end function IPF_constructor

!--------------------------------------------------------------------------
subroutine IPF_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! destructor for the IPFmap_T Class
 
IMPLICIT NONE

type(IPFmap_T), INTENT(INOUT)  :: self 

call reportDestructor('IPFmap_T')

end subroutine IPF_destructor


!--------------------------------------------------------------------------
function get_ipf_wd_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_wd_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get ipf_wd from the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out

out = self%ipf_wd

end function get_ipf_wd_

!--------------------------------------------------------------------------
subroutine set_ipf_wd_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_wd_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set ipf_wd in the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp

self%ipf_wd = inp

end subroutine set_ipf_wd_

!--------------------------------------------------------------------------
function get_ipf_ht_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_ht_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get ipf_ht from the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out

out = self%ipf_ht

end function get_ipf_ht_

!--------------------------------------------------------------------------
subroutine set_ipf_ht_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_ht_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set ipf_ht in the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp

self%ipf_ht = inp

end subroutine set_ipf_ht_

!--------------------------------------------------------------------------
function get_ipf_nthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_nthreads_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get ipf_nthreads from the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out

out = self%ipf_nthreads

end function get_ipf_nthreads_

!--------------------------------------------------------------------------
subroutine set_ipf_nthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_nthreads_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set ipf_nthreads in the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp

self%ipf_nthreads = inp

end subroutine set_ipf_nthreads_

!--------------------------------------------------------------------------
function get_ipf_filename_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_filename_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get ipf_filename from the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)  :: self
character(fnlen)                :: out

out = self%ipf_filename

end function get_ipf_filename_

!--------------------------------------------------------------------------
subroutine set_ipf_filename_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_filename_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set ipf_filename in the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)  :: self
character(fnlen), INTENT(IN)    :: inp

self%ipf_filename = inp

end subroutine set_ipf_filename_

!--------------------------------------------------------------------------
function get_ipf_mode_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_mode_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get ipf_mode from the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)  :: self
character(fnlen)                :: out

out = self%ipf_mode

end function get_ipf_mode_

!--------------------------------------------------------------------------
subroutine set_ipf_mode_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_mode_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set ipf_mode in the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)  :: self
character(fnlen), INTENT(IN)    :: inp

self%ipf_mode = inp

end subroutine set_ipf_mode_

!--------------------------------------------------------------------------
function get_ipf_LaueClass_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_LaueClass_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get ipf_LaueClass from the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)  :: self
integer(kind=irg)               :: out

out = self%ipf_LaueClass

end function get_ipf_LaueClass_

!--------------------------------------------------------------------------
subroutine set_ipf_LaueClass_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_LaueClass_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set ipf_LaueClass in the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)  :: self
integer(kind=irg), INTENT(IN)   :: inp

self%ipf_LaueClass = inp

end subroutine set_ipf_LaueClass_

! !--------------------------------------------------------------------------
! !---------based on Will Lenthe's sphere_sector.hpp module------------------
! !--------------------------------------------------------------------------
! recursive subroutine SphericalTriangle_(self, nRed, nGreen, nBlue)
! !DEC$ ATTRIBUTES DLLEXPORT :: SphericalTriangle_
! !! author: MDG (based on Will Lenthe's sphere_sector.hpp module)
! !! version: 1.0 
! !! date: 09/29/21
! !!
! !! construct a spherical triangle patch for IPF coloring

! use mod_io

! IMPLICIT NONE 

! class(IPFmap_T), INTENT(INOUT)  :: self
! real(kind=dbl), INTENT(IN)      :: nRed(0:2)
! real(kind=dbl), INTENT(IN)      :: nGreen(0:2)
! real(kind=dbl), INTENT(IN)      :: nBlue(0:2)

! type(IO_T)                      :: Message 
! real(kind=dbl)                  :: verts(0:2,0:2), det, ctr(0:2), eps=1.0D-12 
! integer(kind=irg)               :: i

! verts(0,0:2) = nRed 
! verts(1,0:2) = nGreen
! verts(2,0:2) = nBlue 

! det = verts(0,0) * verts(1,1) * verts(2,2) &
!     + verts(0,1) * verts(1,2) * verts(2,0) &
!     + verts(0,2) * verts(1,0) * verts(2,1) &
!     - verts(0,0) * verts(1,2) * verts(2,1) &
!     - verts(0,1) * verts(1,0) * verts(2,2) &
!     - verts(0,2) * verts(1,1) * verts(2,0) 

! if (det.lt.eps) call Message%printError('SphericalTriangle_','spherical triangle must be within single hemisphere')

! ! compute center of spherical triangle (assumes triangle covers less than hemisphere)
! do i=0,2
!   ctr(i) = verts(0,i) + verts(1,i) + verts(2,i)
! end do 
! ctr = ctr/sqrt(sum(ctr*ctr))

! ! now build a general spherical patch (this fills the omega lookup table)
! call self%build_(3, verts, ctr)

! end subroutine SphericalTriangle_

! !--------------------------------------------------------------------------
! recursive subroutine SphericalWedge_(self, nGreen, nBlue)
! !DEC$ ATTRIBUTES DLLEXPORT :: SphericalWedge_
! !! author: MDG (based on Will Lenthe's sphere_sector.hpp module)
! !! version: 1.0 
! !! date: 09/29/21
! !!
! !! construct a spherical triangle wedge for IPF coloring
! !! red is assumed to be at (0,0,1)

! use mod_io

! IMPLICIT NONE 

! class(IPFmap_T), INTENT(INOUT)  :: self
! real(kind=dbl), INTENT(IN)      :: nGreen(0:1)
! real(kind=dbl), INTENT(IN)      :: nBlue(0:1)

! type(IO_T)                      :: Message 
! real(kind=dbl)                  :: verts(0:3,0:2), dot, ctr(0:2), eps=1.0D-12 
! integer(kind=irg)               :: i

! verts = 0.D0
! verts(0,2) = 1.D0 
! verts(1,0:1) = nGreen
! verts(2,2) = -1.D0
! verts(3,0:1) = nBlue 

! ! compute center 
! ctr = (/ nGreen(0) + nBlue(0), nGreen(1) + nBlue(1), 0.D0 /)

! ! handle special case of antipodal blue/green
! dot = nGreen(0) * nBlue(0) + nGreen(1) * nBlue(1)
! if (abs(dot + 1.D0).lt.eps) then
!   ctr(0) = -nGreen(1)
!   ctr(1) =  nGreen(0)
! end if 

! ctr = ctr/sqrt(sum(ctr*ctr))

! ! now build a general spherical patch (this fills the omega lookup table)
! call self%build_(4, verts, ctr)

! end subroutine SphericalWedge_

! !--------------------------------------------------------------------------
! recursive subroutine init_sPatch_(self, N)
! !DEC$ ATTRIBUTES DLLEXPORT :: build_
! !! author: MDG (based on Will Lenthe's sphere_sector.hpp module)
! !! version: 1.0 
! !! date: 09/29/21
! !!

! IMPLICIT NONE 

! class(IPFmap_T), INTENT(INOUT)        :: self

! allocate(self%sphericalPatch%rx(0:N-1))
! allocate(self%sphericalPatch%ry(0:N-1))
! allocate(self%sphericalPatch%center(0:N-1))
! allocate(self%sphericalPatch%normals(0:N-1,0:2))
! allocate(self%sphericalPatch%cutoffs(0:3*N-1))
! allocate(self%sphericalPatch%coeffs(0:N-1,0:3))
! allocate(self%sphericalPatch%cumAngles(0:N))
! allocate(self%sphericalPatch%omega(0:256*N-1))

! end subroutine init_sPatch_

! !--------------------------------------------------------------------------
! recursive subroutine build_(self, N, verts, ctr, filFin)
! !DEC$ ATTRIBUTES DLLEXPORT :: build_
! !! author: MDG (based on Will Lenthe's sphere_sector.hpp module)
! !! version: 1.0 
! !! date: 09/29/21
! !!
! !! construct a spherical triangle patch to map to the unit hemisphere
! !! verts: vertices of spherical patch (maps to evenly spaced points to equator)
! !! ctr  : center of spherical patch (maps to north pole)
! !! filF : fillet fraction, should be 0 for original paper, ~0.05 for perceptually uniform, must be in (0,0.5)
! !!      : verts in CCW order

! use mod_io
! use mod_math

! IMPLICIT NONE 

! class(IPFmap_T), INTENT(INOUT)        :: self
! integer(kind=irg),INTENT(IN)          :: N
! real(kind=dbl),INTENT(IN)             :: verts(0:N-1,0:2)
! real(kind=dbl),INTENT(IN)             :: ctr(0:2)
! real(kind=dbl),INTENT(IN),OPTIONAL    :: filFin

! real(kind=dbl)                        :: center(0:2), vx(0:N-1,0:2), vy(0:N-1,0:2), angles(0:N-1), filF, &
!                                          deltas(0:N-1), delta, radii(0:2*N-1), dRadii(0:2*N-1), thetas(0:5), nn(0:2), &
!                                          r(0:5), c, s, x, y, z, m(0:2), nx, ny, nz, mx, my, mz, den, v1, v2, m1, &
!                                          m2, rhoN(0:N-1), v(0:2), irho(0:256*N-1), normxn(0:2), mag, sm
! integer(kind=irg)                     :: i, j, idx(0:N) 

! associate( sP => self%sphericalPatch )

! call self%init_sPatch_(N)

! center = ctr

! ! build orthogonal coordinate system for center -> each vertex
! do i=0,N-1
!   vy(i,0:2) = cross3( center(0:2), verts(i,0:2) )
!   vx(i,0:2) = cross3( vy(i,0:2), center(0:2) )
!   vx(i,0:2) = vx(i,0:2)/sqrt(sum(vx(i,0:2)**2))
!   vy(i,0:2) = vy(i,0:2)/sqrt(sum(vy(i,0:2)**2))
! end do 
! sP%rx(0:2) = vx(0,0:2)    ! red is the global x direction
! sP%ry(0:2) = vy(0,0:2)    ! global y is perpendicular to x and patch center (z)

! ! compute angles between successive verts
! do i=0, N-1
!   angles(i) = acos( DOT_PRODUCT( vx(i,0:2), vx(mod(i+1,N)) ) )
!   if (i.eq.0) then 
!     sP%cumAngles(0) = angles(0)
!   else
!     sP%cumAngles(i) = sP%cumAngles(i-1) + angles(i)  ! needs to be checked !!! std::partial_sum(angles, angles+N, cumAngles+1);
!   end if
! end do 

! ! compute normals of circles defining edges of domain
! do i=0, N-1
!   sP%normals(i,0:2) = cross3( verts(i,0:2), verts( mod(i,N), 0:2) )
!   sP%normals(i,0:2) = sP%normals(i,0:2) / sqrt( sum( sP%normals(i,0:2)**2 ) )
! end do

! ! compute cutoff angles for filleting
! filF = 0.D0 
! if (present(filFin)) filF = filFin
! deltas = 2.D0*cPi/dble(N)   ! for now just evenly space deltas
! do i=0, N-1                 ! loop over edges
!   do j=0, 2                 ! loop over verts
!     if (j.eq.2) then
!       delta = 0.D0 
!     else 
!       if (j.eq.0) then 
!         delta = filF * deltas(i)
!       else
!         delta = -filF * deltas(i)
!       end if
!     end if 
!     if (j.eq.0) then 
!       sP%cutoffs(i*3+j) = sP%cumAngles(i) + delta
!     else
!       sP%cutoffs(i*3+j) = sP%cumAngles(i+1) + delta
!     end if
!   end do 
! end do

! ! numerically compute r and dr/dtheta at transition points between linear and filleted regions
! do i= 0, N-1  ! loop over edges
! ! compute angles where radius needs to be calculated
!   dT = 0.10D0   ! angular offset from vertex -> transition points (for numerical derivative calculation), ~1/2 degree
!   hf = 0.010D0  ! fractional distance numerical derivative calculation (should probably be <= 1)
!   thetas = (/ sP%cutoffs(3*i+0) - hf * dT, & ! symmetric points for derivative calculation
!               sP%cutoffs(3*i+0)          , & ! first transition
!               sP%cutoffs(3*i+0) + hf * dT, & ! symmetric points for derivative calculation
!               sP%cutoffs(3*i+1) - hf * dT, & ! symmetric points for derivative calculation
!               sP%cutoffs(3*i+1)          , & ! second transition
!               sP%cutoffs(3*i+1) + hf * dT /) ! symmetric points for derivative calculation

! ! apply rotations and compute radius/angle at each point
!   do j= 0, 5
! ! compute normal of circle at desired angle (ry rotated about center)
!     c = cos(thetas(j) / 2.D0)
!     s = sin(thetas(j) / 2.D0)

! ! q * n (w == 0 since rotation axis is perpendicular to vector)
!     x = c * sP%ry(0) + s * (center(1) * sP%ry(2) - center(2) * sP%ry(1))
!     y = c * sP%ry(1) + s * (center(2) * sP%ry(0) - center(0) * sP%ry(2))
!     z = c * sP%ry(2) + s * (center(0) * sP%ry(1) - center(1) * sP%ry(0))

! ! q * n * q.conj() (normal of circle at desired angle)
!     m = (/  x * c + s * (z * center(1) - y * center(2)), &
!             y * c + s * (x * center(2) - z * center(0)), &
!             z * c + s * (y * center(0) - x * center(1)) /)

! ! now compute intersection of two unit circles at origin w/ normals v and normals(edge)
!     nx = sP%normals(i,0)
!     ny = sP%normals(i,1)
!     nz = sP%normals(i,2)
!     mx = m(0)
!     my = m(1)
!     mz = m(2)
!     den = sqrt( nx * nx * ( my * my + mz * mz ) + ny * ny * ( mz * mz + mx * mx ) + nz * nz * ( mx * mx + my * my ) &
!                - 2.D0 * ( nz * nx * mz * mx + nx * ny * mx * my + ny * nz * my * mz ) )

! ! intersection of two circles (point along edge i at angle thetas(j))
!     v = (/ (ny * mz - nz * my) / den, &
!            (nz * mx - nx * mz) / den, & 
!            (nx * my - ny * mx) / den /)

!     if( DOT_PRODUCT( v(0:2), center(0:2) ).lt.0.D0) v = -v  ! select intersection point closest to center
!     r(j) = acos( DOT_PRODUCT( v(0:2), center(0:2) ))        ! compute angle from center -> edge at this theta
!   end do 

! ! save radii and compute derivative
!   radii (i*2+0) = r(1)
!   radii (i*2+1) = r(4)
!   dRadii(i*2+0) = (r(2) - r(0)) / (hf * dT * 2.D0)
!   dRadii(i*2+1) = (r(5) - r(3)) / (hf * dT * 2.D0)
! end do

! ! compute polynomial coefficients to remove discontinuity in r
! do i= 0, N-1 ! loop over edge
!   j = mod(i+1, N)   ! index of next edge
!   v1 = radii (i*2+1)                    ! value of radius at transition point in edge i (near edge j)
!   v2 = radii (j*2+0)                    ! value of radius at transition point in edge j (near edge i)
!   m1 = dradii(i*2+1) * filf * angles(i) ! value of d(radius)/d(theta) at transition point (multiply by range of -1->0 to correct derivative scaling)
!   m2 = dradii(j*2+0) * filf * angles(j) ! value of d(radius)/d(theta) at transition point (multiply by range of  0->1 to correct derivative scaling)
!   sP%coeffs(i,0) = ( m1 + m2 + v1        - v2       ) / 4.D0
!   sP%coeffs(i,1) = (-m1 + m2                        ) / 4.D0
!   sP%coeffs(i,2) = (-m1 - m2 - v1 * 3.D0 + v2 * 3.D0) / 4.D0
!   sP%coeffs(i,3) = ( m1 - m2 + v1 * 2.D0 + v2 * 2.D0) / 4.D0
! end do

! ! build lookup table for nonlinear hue adjustment (to keep hue at vert i i/N)
! ! this is a numerical lookup table to solve A.6

! ! compute the fractional angle of each vertex
! do i=1,N-2
!   v(0:2) = verts(i,0:2) - center(0:2)
!   angle = atan2( DOT_PRODUCT( sP%ry, v ), DOT_PRODUCT( sP%rx, v ) )
!   if (angle.lt.0.D0) angle = angle + 2.D0*cPi  
!   rhoN(i-1) = angle / (2.D0*cPi)
! end do 
! rhoN(N-1) = 1.D0

! ! create evenly spaced list for angle from 0->1
! irho = (/ (dble(i)/dble(256.D0*N-1), i=0,256*N-1) /)

! ! compute the distance to the sector edge at each angle (in irho)
! sP%omega(0) = 0.D0
! do i=0, 256*N-1
! ! create vector normal to center at angle irho(i)
!   s = sin(2.D0*cPi * irho(i));
!   c = cos(2.D0*cPi * irho(i));
!   do j=0,2 
!     nn(j) = sP%rx(j) * s - sP%ry(j) * c   ! std::transform(rx, rx+3, ry, n, (s, c)(Real i, Real j){return i * s - j * c;});
!   end do

! ! determine which edge is closest and compute distance to edge
!   j = 0 
!   do while (rhoN(j).lt.irho(i)) 
!     j=j+1 
!   end do
!   normxn = cross3( SP%normals(j,0:2), nn(0:2) )
!   mag = sqrt(sum(normxn*normxn))
!   sP%omega(i+1) = acos( DOT_PRODUCT(normxn(0:2), center(0:2)) / mag)
! end do

! ! get the offset to the vertices
! idx(N) = 0
! do i=0,N-1
!   j = 0 
!   do while (irho(j).lt.rhoN(i)) 
!     j=j+1 
!   end do
!   idx(i+1) = j
! end do 

! ! normalize
! do i=0,N-1
!   sm = sum(sP%omega(idx(i):idx(i+1)))
!   sP%omega(idx(i):idx(i+1)) = sP%omega(idx(i):idx(i+1)) / sm
! end do

! ! integrate to obtain the final lookup table
! do i=1,N-1
!     sP%omega(i) = sP%omega(i) + sP%omega(i-1)
! end do

! end associate

! end subroutine build_

! !--------------------------------------------------------------------------
! recursive subroutine toHemi_(self, n, tht, phi)
! !DEC$ ATTRIBUTES DLLEXPORT :: toHemi_
! !! author: MDG (based on Will Lenthe's sphere_sector.hpp module)
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! n  : unit direction to color
! !! tht: fractional azimuthal angle (0,1) maps to (0,2*pi)
! !! phi: fractional polar angle (0,1) maps to (0,pi/2)

! IMPLICIT NONE 

! class(IPFmap_T), INTENT(INOUT)        :: self
! real(kind=dbl),INTENT(IN)             :: n(0:2)
! real(kind=dbl),INTENT(OUT)            :: tht
! real(kind=dbl),INTENT(OUT)            :: phi

! real(kind=dbl)                        :: v(0:2), n0(0:2), angle, eps=1.D-12, nxc(0:2), x, den
! integer(kind=irg)                     :: iOmg, idx, i, j 

! associate( sP => self%sphericalPatch )

! ! compute angle with red direction
! n0 = n      ! in case input and output overlap
! v = n0 - sP%center 
! angle = atan2( DOT_PRODUCT(sP%ry(0:2), v(0:2)), DOT_PRODUCT(sP%rx(0:2), v(0:2)) )
! if (angle.lt.0.D0) angle = angle + 2.D0*cPi 
! tht = angle / (2.D0*cPi)

! ! apply adaptive hue gradient from omega lookup table
! tht = tht * dble(self%N*256-1)      ! rescale from (0,1) to lookup table size
! iOmg = int(tht)                     ! get index of lower bound in lookup table
! if ((iOmg+1).lt.(self%N*256)) then  ! don't go out of bounds if tht == 1
! ! linearly interpolate from lookup table
!   tht = sP%omega(iOmg) + (sP%omega(iOmg+1) - sP%omega(iOmg)) * (dble(iOmg+1) - tht)
! else
!   tht = 1.D0
! end if 

! ! compute polar angle
! ! determine which region this angle falls in
! j=0
! do while(sP%cutoffs(j).lt.angle) j = j+1
! idx = j 
! i = idx/3  
! phi = acos( DOT_PRODUCT(n0(0:2), sP%center(0:2)) ) ! angle between center and point

! if (phi.gt.eps) then    ! avoid divide by zero issues
! ! normalize polar angle
!   select case(idx - i * 3) 
!         case(1) ! in linear region, normalize angle by max possible angle
!           nxc = cross3(n0(0:2), sP%center(0:2))     ! normal of arc through n/center
!           v = cross3(sP%normals(i), nxc)            ! intersection of two circles (edge of patch in direction tht)
!           v = v/sqrt(sum(v*v))
!           phi = phi/( 2.D0 * acos(DOT_PRODUCT(v(0:2), sP%center(0:2))) )  ! compute fractional progress ot edge

!         case(0) ! in first fillet
!           j = mod( (i+self%N-1), self%N)            ! get i-1 with periodic boundary conditions
!           x =  (angle - sP%cumAngles(i)) / (sP%cutoffs(idx) - sP%cumAngles(i))
!           den = sP%coeffs(j,0) * x * x * x + sP%coeffs(j,1) * x * x + sP%coeffs(j,2) * x + sP%coeffs(j,3)
!           phi = phi/ (2.D0 * maxval( (/phi, den/) ) ! normalize, clipping at 1/2

!         case(2) ! in second fillet
!           x = -(angle - sP%cumAngles(i+1)) / (sP%cutoffs(idx-1) - sP%cumAngles(i+1))
!           den = sP%coeffs(i,0) * x * x * x + sP%coeffs(i,1) * x * x + sP%coeffs(i,2) * x + sP%coeffs(i,3)
!           phi = phi/ (2.D0 * maxval( (/phi, den/) ) ! normalize, clipping at 1/2
        
!   end select 
! end if 

! end associate

! end subroutine toHemi_

! !--------------------------------------------------------------------------
! recursive function toColor_(self, n, h2r, wCn, nTh) result(rgb)
! !DEC$ ATTRIBUTES DLLEXPORT :: build_
! !! author: MDG (based on Will Lenthe's implementation in sphere_sector.hpp) 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! compute coloring for a unit direction in the fundamental sector
! !! n  : unit direction in fundamental sector (undefined behavior for directions outside of sector)
! !! rgb: location to write color (0,1)
! !! h2r: hsl2rgb like coloring function to use with h, s, and l in (0,1) and output as (0,1) rgb: void(Real const * const hsl, Real * const rgb)
! !! wCn: white/black center (north/south hemisphere)
! !! nTh: should theta be negated (to create enatiomorph / reverse color progression)

! IMPLICIT NONE 

! class(IPFmap_T), INTENT(INOUT)        :: self
! real(kind=dbl), INTENT(IN)            :: n(0:2)
! character(*),INTENT(IN)               :: h2r
! logical,INTENT(IN)                    :: wCn
! logical,INTENT(IN)                    :: nTh
! real(kind=dbl)                        :: rgb(0:2)

! call self%toHemi(n, rgb(0), rgb(2))
! rgb(1) = 1.D0   ! fully saturated
! if (nTh.eqv..TRUE.) rgb(0) = 1.D0 - rgb(0)
! if (wCn.eqv..TRUE.) rgb(2) = 1.D0 - rgb(2)

! ! transform to rgb using selected coloring function
! if (trim(h2r).eq.'hsl2rgb') then 
!   rgb = self%hsl2rgb_(rgb)
! end if 


! end function toColor_

! !--------------------------------------------------------------------------
! !----the following routines are converted from Will Lenthe's code base-----
! !--------------------------------------------------------------------------
! recursive function rot2d(xy, c, s) result(xyr)
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! rotate a 2d vector 
! !! xy: vector to rotate
! !! c : cosine of rotation angle
! !! s : sine of rotation angle

! IMPLICIT NONE 

! real(kind=dbl),INTENT(IN)       :: xy(0:1)
! real(kind=dbl),INTENT(IN)       :: c
! real(kind=dbl),INTENT(IN)       :: s
! real(kind=dbl)                  :: xyr(0:1)

! xyr(0) = c*xy(0) - s*xy(1) 
! xyr(1) = s*xy(0) + c*xy(1) 

! end function rot2d

! !-----------------------------------------------
! recursive function mir2d(xy, c, s) result(xyr)
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! mirror a 2d vector 
! !! xy: vector to mirror
! !! c : cosine of mirror plane angle
! !! s : sine of mirror plane angle

! IMPLICIT NONE 

! real(kind=dbl),INTENT(IN)       :: xy(0:1)
! real(kind=dbl),INTENT(IN)       :: c
! real(kind=dbl),INTENT(IN)       :: s
! real(kind=dbl)                  :: xyr(0:1)

! real(kind=dbl)                  :: d, pt(0:1)

! d = c*xy(0) + s*xy(1) 
! pt = (/ d*c, d*s /)

! xyr = pt - (xy - pt)

! end function mir2d

! !-----------------------------------------------
! recursive function r111(n) result(nr)
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! move a point in the first octant to the maximum z rotation possible by 120 @ 111
! !! n: vector to rotate
! !! return: true if z was already maximize, false otherwise

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: nr

! integer(kind=irg)               :: idx(1)

! nr = .FALSE.
! idx = maxloc(n)
! if (n(2).eq.n(idx(1))) nr = .TRUE.
! if (nr.eq..FALSE.) n = cshift(n, 2-idx(1))

! end function r111

! !-----------------------------------------------
! recursive function b1   (n)  result(res)   ! -1
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! reduce a direction to the -1 fundamental sector
! !! n : unit direction to reduce (in place)
! !! return  : true/false if n was/wasn't already in the fundamental sector

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 

! res = .FALSE.
! if (n(2).lt.0.D0) then 
!   res = .TRUE.
!   n = -n
! end if 
! res = .not.res

! end function b1

! !-----------------------------------------------
! recursive function _211 (n)  result(res)   ! 211
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! reduce a direction to the 211 fundamental sector
! !! n : unit direction to reduce (in place)
! !! return  : true/false if n was/wasn't already in the fundamental sector

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 

! res = .FALSE.
! if (n(2).lt.0.D0) then 
!   res = .TRUE.
!   n(1:2) = -n(1:2)
! end if 
! res = .not.res

! end function _211

! !-----------------------------------------------
! recursive function _211r(n)  result(res)   ! 211 rotated 45 @ 2 (2 fold @ xy)
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! reduce a direction to the 211 fundamental sector rotated 45 @ z
! !! n : unit direction to reduce (in place)
! !! return  : true/false if n was/wasn't already in the fundamental sector

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 

! real(kind=dbl)                  :: t

! res = .FALSE.
! if (n(2).lt.0.D0) then 
!   res = .TRUE.
!   t = n(0)
!   n(0) = n(1)
!   n(1) = t
!   n(2) = -n(2)
! end if 
! res = .not.res

! end function _211r

! !-----------------------------------------------
! recursive function _121 (n)  result(res)   ! 121
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! reduce a direction to the 121 fundamental sector
! !! n : unit direction to reduce (in place)
! !! return  : true/false if n was/wasn't already in the fundamental sector

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 

! res = .FALSE.
! if (n(2).lt.0.D0) then 
!   res = .TRUE.
!   n(0) = -n(0)
!   n(2) = -n(2)
! end if 
! res = .not.res

! end function _121

! !-----------------------------------------------
! recursive function _112 (n)  result(res)   ! 112
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! reduce a direction to the 112 fundamental sector
! !! n : unit direction to reduce (in place)
! !! return  : true/false if n was/wasn't already in the fundamental sector

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 

! res = .FALSE.
! if (n(1).lt.0.D0) then 
!   res = .TRUE.
!   n(0) = -n(0)
!   n(1) = -n(1)
! end if 
! res = .not.res

! end function _112

! !-----------------------------------------------
! recursive function _m11 (n)  result(res)   ! m11
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! reduce a direction to the m11 fundamental sector
! !! n : unit direction to reduce (in place)
! !! return  : true/false if n was/wasn't already in the fundamental sector

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 

! res = .FALSE.
! if (n(0).lt.0.D0) then 
!   res = .TRUE.
!   n(0) = -n(0)
! end if 
! res = .not.res

! end function _m11

! !-----------------------------------------------
! recursive function _1m1 (n)  result(res)   ! 1m1
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! reduce a direction to the 1m1 fundamental sector
! !! n : unit direction to reduce (in place)
! !! return  : true/false if n was/wasn't already in the fundamental sector

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 

! res = .FALSE.
! if (n(1).lt.0.D0) then 
!   res = .TRUE.
!   n(1) = -n(1)
! end if 
! res = .not.res

! end function _1m1

! !-----------------------------------------------
! recursive function _11m (n)  result(res)   ! 11m
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! reduce a direction to the 11m fundamental sector
! !! n : unit direction to reduce (in place)
! !! return  : true/false if n was/wasn't already in the fundamental sector

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 

! res = .FALSE.
! if (n(2).lt.0.D0) then 
!   res = .TRUE.
!   n(2) = -n(2)
! end if 
! res = .not.res

! end function _11m

! !-----------------------------------------------
! recursive function _222r(n)  result(res)   ! 222r
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! reduce a direction to the 222r fundamental sector
! !! n : unit direction to reduce (in place)
! !! return  : true/false if n was/wasn't already in the fundamental sector

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 

! logical                         :: ret 

! res = _211r(n)        ! first bring to +z hemisphere with rotation about xy
! if(n(1).lt.n(0)) then ! are we below the y==x line
!   ret = _112(n)
!   res = .FALSE.
! end if 

! end function _222r

! !-----------------------------------------------
! recursive function mm2r (n)  result(res)   ! mm2r
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! reduce a direction to the mm2r fundamental sector
! !! n : unit direction to reduce (in place)
! !! return  : true/false if n was/wasn't already in the fundamental sector

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 

! real(kind=dbl)                  :: ax 

! res = _112(n)     ! first bring to +y hemisphere with rotation about z
! ax = abs(n(0))
! if(ax.gt.n(1)) then
!   if (n(0).lt.0.D0) then 
!     n(0) = -n(1) 
!   else 
!     n(0) = n(1) 
!   end if 
!   n(1) = ax
!   res = .FALSE.
! end if 

! end function mm2r

! !-----------------------------------------------
! recursive function _4   (n)  result(res)   ! 4
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! reduce a direction to the 4 fundamental sector
! !! n : unit direction to reduce (in place)
! !! return  : true/false if n was/wasn't already in the fundamental sector

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 

! real(kind=dbl)                  :: t 

! res = _112(n)           ! first move to +y with 2 fold
! if(n(0).lt.0.D0) then   ! need to move to +x with 4 fold
!   n(0) = -n(0)
!   t = n(0) 
!   n(0) = n(1) 
!   n(1) = t
!   res = .FALSE.
! end if 

! end function _4

! !-----------------------------------------------
! recursive function b4   (n)  result(res)   ! -4
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! reduce a direction to the -4 fundamental sector
! !! n : unit direction to reduce (in place)
! !! return  : true/false if n was/wasn't already in the fundamental sector

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 

! real(kind=dbl)                  :: t 

! res = .TRUE.
! if(n(2).lt.0.D0) then   ! need to move to +z with -4
!   n(0) = -n(0)
!   t = n(0)
!   n(0) = n(1)
!   n(1) = t
!   n(2) = -n(2)
!   res = .FALSE.
! end if 
 
! if(n(1).lt.0.D0) then   ! need to move to +y with 2
!   n(0:1) = -n(0:1)
!   ret = false;
!   res = .FALSE.
! end if 

! end function b4

! !-----------------------------------------------
! recursive function _4mm (n)  result(res)   ! 4mm
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! reduce a direction to the 4mm fundamental sector
! !! n : unit direction to reduce (in place)
! !! return  : true/false if n was/wasn't already in the fundamental sector

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 

! real(kind=dbl)                  :: t 

! res = _4(n)               ! first move to first quadrant with 4 fold
! if(n(1).gt.n(0)) then     ! use xy mirror if needed
!   t = n(0)
!   n(0) = n(1)
!   n(1) = t 
!   res = .FALSE.
! end if 

! end function _4mm

! !-----------------------------------------------
! recursive function _3   (n)  result(res)   ! 3
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! reduce a direction to the 3 fundamental sector
! !! n : unit direction to reduce (in place)
! !! return  : true/false if n was/wasn't already in the fundamental sector

! use mod_IO 

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 

! type(IO_T)                      :: Message 

! real(kind=dbl)                  :: t 
! integer(kind=irg)               :: idx

! ! determine sector
! t = abs(n(1) / n(0))      ! |tan(theta)|
! ! 0, 1, or 2 for first, second, or third sector, initialize with easy cases
! idx = 0
! if (n(1).lt.0D0) idx = 2
! if( (n(0).lt.0.D0).and.(t.le.sqrt(3.D0)) ) idx = 1    ! check for second sector

! ! apply rotation
! select case(idx) 
!   case(0)     ! in first sector, we're done
!     res = .TRUE.
!   case(1)     ! in second sector, rotate -120 @ z
!     n = rot2d(n, -0.50, -sqrt(0.75D0)) 
!     res = .FALSE. 
!   case(2)     ! in third sector, rotate 120 @ z
!     n = rot2d(n, -0.50,  sqrt(0.75D0)) 
!     res = .FALSE. 
!   case default 
!     call Message%printError('_3','unhandled 3 fold case')
! end select

! end function _3

! !-----------------------------------------------
! recursive function _3r  (n)  result(res)   ! 3 rotated 30 degrees (use for other sectors)
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! reduce a direction to the 3 fundamental sector rotated 30 degrees about z
! !! n : unit direction to reduce [in place]
! !! return  : true/false if n was/wasn't already in the fundamental sector

! use mod_IO 

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 

! type(IO_T)                      :: Message 

! real(kind=dbl)                  :: t, eps = 1.D-12 
! integer(kind=irg)               :: idx

! ! determine sector
! t = abs(n(1) / maxval( (/ abs(n(0)), eps /) )      ! |tan(theta)|
! ! 0, 1, or 2 for first, second, or third sector, initialize with easy cases
! idx = 2
! if (n(0).lt.0D0) idx = 1
! if( (n(1).lt.0.D0).and.(t.ge.1.D0/sqrt(3.D0)) ) idx = 0    ! check for first sector

! ! apply rotation
! select case(idx) 
!   case(0)     ! in first sector, we're done
!     res = .TRUE.
!   case(1)     ! in second sector, rotate -120 @ z
!     n = rot2d(n, -0.50, -sqrt(0.75D0)) 
!     res = .FALSE. 
!   case(2)     ! in third sector, rotate 120 @ z
!     n = rot2d(n, -0.50,  sqrt(0.75D0)) 
!     res = .FALSE. 
!   case default 
!     call Message%printError('_3r','unhandled 3 fold case')
! end select

! end function _3r

! !-----------------------------------------------
! recursive function _31m (n)  result(res)   ! 31m
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! reduce a direction to the 31m fundamental sector
! !! n : unit direction to reduce [in place]
! !! return  : true/false if n was/wasn't already in the fundamental sector

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 

! real(kind=dbl)                  :: t, eps = 1.D-12 

! res = _3(n)       ! first reduce to 3 fold fz
! t = n(1) / n(0)
! ! check if we're above 60 degrees (< eps instead of signbit to handle divide by 0)
! if ( (n(1).lt.eps) .or. (.not( (0.D0.le.t).and.(t.le.sqrt(3.D0))) ) ) then 
!   n = mir2d(n, 0.5D0, sqrt(0.75D0))
!   res = .FALSE.
! end if 

! end function _31m

! !-----------------------------------------------
! recursive function _6   (n)  result(res)   ! 6
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!
! !! reduce a direction to the 6 fundamental sector
! !! n : unit direction to reduce [in place]
! !! return  : true/false if n was/wasn't already in the fundamental sector

! use mod_IO

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 

! type(IO_T)                      :: Message 

! real(kind=dbl)                  :: t, eps = 1.D-12
! integer(kind=irg)               :: idx 

! res = _112(n)       ! first bring to +y hemisphere

! ! determine sector
! t = abs(n(1)) / maxval( (/ abs(n(0)), eps /) ) ! |tan(theta)| if in +y hemisphere
! idx = 0
! if (n(0).lt.0.D0) idx = 2 
! if(t.gt.sqrt(3.D0)) idx = 1                    ! check for second sector

! ! apply rotation
! select case(idx) 
!   case(0)     ! in first sector, we're done
!     res = .TRUE.
!   case(1)     ! in second sector, rotate -120 @ z
!     n = rot2d(n,  0.50, -sqrt(0.75D0)) 
!     res = .FALSE. 
!   case(2)     ! in third sector, rotate 120 @ z
!     n = rot2d(n, -0.50, -sqrt(0.75D0)) 
!     res = .FALSE. 
!   case default 
!     call Message%printError('_6','unhandled 6 fold case')
! end select
  
! end function _6

! HERE


! !-----------------------------------------------
! recursive function _12m1(n)  result(res)   ! 12/m1
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function _12m1

! !-----------------------------------------------
! recursive function _112m(n)  result(res)   ! 112/m
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function _112m

! !-----------------------------------------------
! recursive function _222 (n)  result(res)   ! 222
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function _222



! !-----------------------------------------------
! recursive function mm2  (n)  result(res)   ! mm2
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function mm2



! !-----------------------------------------------
! recursive function mmm  (n)  result(res)   ! mmm
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function mmm

! !-----------------------------------------------
! recursive function mmmr (n)  result(res)   ! mmmr
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function mmmr



! !-----------------------------------------------
! recursive function _4m  (n)  result(res)   ! 4/m
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function _4m

! !-----------------------------------------------
! recursive function _422 (n)  result(res)   ! 422
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function _422



! !-----------------------------------------------
! recursive function b42m (n)  result(res)   ! -42m
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function b42m

! !-----------------------------------------------
! recursive function b4m2 (n)  result(res)   ! -4m2
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function b4m2

! !-----------------------------------------------
! recursive function _4mmm(n)  result(res)   ! 4/mmm
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function _4mm

! !-----------------------------------------------
! recursive function _3   (n)  result(res)   ! 3
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function _3

! !-----------------------------------------------
! recursive function _3r  (n)  result(res)   ! 3 rotated 30 degrees (use for other sectors)
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function _3r

! !-----------------------------------------------
! recursive function b3   (n)  result(res)   ! -3
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function b3

! !-----------------------------------------------
! recursive function _321 (n)  result(res)   ! 321
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function _321

! !-----------------------------------------------
! recursive function _312 (n)  result(res)   ! 312
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function _312

! !-----------------------------------------------
! recursive function _3m1 (n)  result(res)   ! 3m1
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function _3m1

! !-----------------------------------------------
! recursive function _31m (n)  result(res)   ! 31m
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function _31m

! !-----------------------------------------------
! recursive function b3m1 (n)  result(res)   ! -3m1
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function b3m1

! !-----------------------------------------------
! recursive function b31m (n)  result(res)   ! -31m
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function b31m



! !-----------------------------------------------
! recursive function b6   (n)  result(res)   ! -6
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function b6

! !-----------------------------------------------
! recursive function _6m  (n)  result(res)   ! 6/m
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function _6m

! !-----------------------------------------------
! recursive function _622 (n)  result(res)   ! 622
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function _622

! !-----------------------------------------------
! recursive function _6mm (n)  result(res)   ! 6mm
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function _6mm

! !-----------------------------------------------
! recursive function b6m2 (n)  result(res)   ! -6m2
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function b6m2

! !-----------------------------------------------
! recursive function b62m (n)  result(res)   ! -62m
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function b62m

! !-----------------------------------------------
! recursive function _6mmm(n)  result(res)   ! 6/mmm
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function _6mmm

! !-----------------------------------------------
! recursive function _23  (n)  result(res)   ! 23
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function _23

! !-----------------------------------------------
! recursive function mb3  (n)  result(res)   ! m3
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function mb3

! !-----------------------------------------------
! recursive function _432 (n)  result(res)   ! 432
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function _432

! !-----------------------------------------------
! recursive function b43m (n)  result(res)   ! -43m
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function b43m

! !-----------------------------------------------
! recursive function mb3m (n)  result(res)   ! m3m
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/01/21
! !!

! IMPLICIT NONE 

! real(kind=dbl),INTENT(INOUT)    :: n(0:2)
! logical                         :: res 


! end function mb3m






!--------------------------------------------------------------------------
!----------the following routines are based on the DREAM.3D code base------
!--------------------------------------------------------------------------
recursive function in_generic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMin, chiMax ) result(inside)
!! author: MDG 
!! version: 1.0 
!! date: 09/06/21
!!
!! generic unit triangle test 

IMPLICIT NONE 

real(kind=dbl),INTENT(IN)       :: etaDeg
real(kind=dbl),INTENT(IN)       :: chiDeg
real(kind=dbl),INTENT(IN)       :: etaMin 
real(kind=dbl),INTENT(IN)       :: etaMax 
real(kind=dbl),INTENT(IN)       :: chiMin 
real(kind=dbl),INTENT(IN)       :: chiMax 
logical                         :: inside 

real(kind=dbl),parameter        :: eps = 1.0D-12

inside = .FALSE.

if ( (etaDeg.ge.etaMin) .and. (etaDeg.le.(etaMax+eps)) ) then 
  if ( (chiDeg.ge.chiMin) .and. (chiDeg.le.(chiMax+eps)) )  then 
    inside = .TRUE.
  end if 
end if 

end function in_generic_unit_triangle

!--------------------------------------------------------------------------
recursive function in_cubic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMax ) result(inside)
!! author: MDG 
!! version: 1.0 
!! date: 09/06/21
!!
!! cubic unit triangle test 

IMPLICIT NONE 

real(kind=dbl),INTENT(IN)       :: etaDeg
real(kind=dbl),INTENT(IN)       :: chiDeg
real(kind=dbl),INTENT(IN)       :: etaMin
real(kind=dbl),INTENT(IN)       :: etaMax 
real(kind=dbl),INTENT(OUT)      :: chiMax 
logical                         :: inside 

real(kind=dbl),parameter        :: eps = 1.0D-12

inside = .FALSE.

! does this lie in the cubic high unit triangle ? 
if (etaDeg.ge.45D0) then 
  chiMax = dsqrt(1.D0 / (2.D0 + tan((90.D0 - etaDeg)*dtor)**2 ))
else 
  chiMax = dsqrt(1.D0 / (2.D0 + tan(etaDeg*dtor)**2 ))
end if 
if (chiMax.gt.1.D0) chiMax = 1.D0
chiMax = dacos(chiMax) * rtod
if ( (etaDeg.ge.0.D0) .and. (etaDeg.le.(etaMax + eps)) ) then 
  if ( (chiDeg.ge.0.D0) .and. (chiDeg.le.((chiMax + eps ) * rtod) ) )  then 
    inside = .TRUE.
  end if 
end if 

end function in_cubic_unit_triangle

!--------------------------------------------------------------------------
recursive function get_ipf_RGB_(self, sampleDir, qu, sym, Pm, clr ) result(RGB)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_RGB_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! returns an RGB triplet for the specific Laue Class, sample direction, and grain orientation
!! This code is based on the LaueOps routines from DREAM.3D, but combines all Laue classes into
!! a single routine.
!!
!! At present, only Euler coloring and standard Laue Class RGB coloring are available.
!! The more extensive color schemes proposed in (Nolze, G., & Hielscher, R. (2016). Orientations–perfectly colored. 
!! Journal of Applied Crystallography, 49(5), 1786-1802.), which have minimal color discontinuities
!! will be implemented in a future version of this module.

use mod_quaternions
use mod_rotations 
use mod_io
use mod_colorspace

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)        :: self
real(kind=dbl), INTENT(IN)            :: sampleDir(3)
type(Quaternion_T), INTENT(IN)        :: qu
type(QuaternionArray_T), INTENT(IN)   :: sym 
integer(kind=irg), INTENT(IN)         :: Pm
integer(kind=ish)                     :: RGB(3)
type(colorspace_T),INTENT(IN),OPTIONAL:: clr

type(e_T)                             :: eu
type(q_T)                             :: qq
type(IO_T)                            :: Message
type(Quaternion_T)                    :: qs, qc, qtest 

real(kind=dbl)                        :: e(3), tp, hsl(0:2) 
real(kind=dbl)                        :: refDir(3), chi, eta, etaDeg, chiDeg, RGBd(3), ma 
real(kind=dbl)                        :: chiMin, chiMax, etaMin, etaMax
integer(kind=irg)                     :: i 
logical                               :: inside

! do a simple test 
! qtest = Quaternion_T( qd = (/ 1.D0, 0.D0, 0.D0, 0.D0 /) )
! RGB = self%get_RGB_cubic_high_( sampleDir, qtest, sym, Pm )
! write (*,*) ' RGB (001) = ', RGB 
! qtest = Quaternion_T( qd = (/ cos(0.5D0*45.D0*dtor), 0.D0, sin(0.5D0*45.D0*dtor), 0.D0 /) )
! RGB = self%get_RGB_cubic_high_( sampleDir, qtest, sym, Pm )
! write (*,*) ' RGB (101) = ', RGB 
! stop

! use a simple Euler angle color scheme, mapping all of Euler space linearly onto RGB color space
if (trim(self%ipf_mode).eq.'Euler') then 
  tp = 2.D0 * cPi 
! scale the full Euler space and map onto the RGB space
  qq = q_T( qdinp = qu%get_quatd() )
  eu = qq%qe()
  e = eu%e_copyd()
  RGB = (/ int(e(1) * 255/ tp), int(e(2) * 255 / cPi), int(e(3) * 255/ tp) /) 
end if

! use the standard IPF colors based on the TSL Euler angles
if ( (trim(self%ipf_mode).ne.'Euler')  ) then 
  RGB = (/ 0, 0, 0 /)

  findloop: do i=1,Pm 
! apply the symmetry operation
    qs = sym%getQuatfromArray(i)
    qc = qs * qu
    call qc%quat_pos()
! convert the sampleDir to refDir
    refDir = qc%quat_Lp( sampleDir )
    if (refDir(3).lt.0.D0) refDir = -refDir
! convert to stereographic projection angles 
    chi = acos(refDir(3))
    eta = atan2(refDir(2),refDir(1))
    etaDeg = eta * rtod
    chiDeg = chi * rtod
! does this lie in the unit triangle ? 
    select case (self%ipf_LaueClass)
      case(1,2)     ! triclinic, monoclinic
          etaMin=0.D0
          etaMax=180.D0
          chiMin=0.D0
          chiMax=90.D0
          inside = in_generic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMin, chiMax )
          inside = .TRUE.
      case(3,4)     ! orthorhombic, tetragonal-low
          etaMin=0.D0
          etaMax=90.D0
          chiMin=0.D0
          chiMax=90.D0
          inside = in_generic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMin, chiMax )
      case(5)       ! tetragonal-high
          etaMin=0.D0
          etaMax=45.D0
          chiMin=0.D0
          chiMax=90.D0
          inside = in_generic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMin, chiMax )
      case(6)       ! trigonal-low
          etaMin=-120.D0
          etaMax=0.D0
          chiMin=0.D0
          chiMax=90.D0
          inside = in_generic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMin, chiMax )
      case(7)       ! trigonal-high
          etaMin=-90.D0
          etaMax=-30.D0
          chiMin=0.D0
          chiMax=90.D0
          inside = in_generic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMin, chiMax )
      case(8)       ! hexagonal-low
          etaMin=0.D0
          etaMax=60.D0
          chiMin=0.D0
          chiMax=90.D0
          inside = in_generic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMin, chiMax )
      case(9)       ! hexagonal-high
          etaMin=0.D0
          etaMax=30.D0
          chiMin=0.D0
          chiMax=90.D0
          inside = in_generic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMin, chiMax )
      case(10)      ! cubic-low
          etaMin=0.D0
          etaMax=90.D0
          inside = in_cubic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMax )
      case(11)      ! cubic-high
          etaMin=0.D0
          etaMax=45.D0
          inside = in_cubic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMax )
      case default 
        call Message%printError('get_ipf_RGB_',' Unknown Laue class')
    end select 

    if (inside.eqv..TRUE.) then ! if .TRUE., then get the correct color and exit the loop
      chiDeg = chiDeg / chiMax
      if (chiDeg.gt.1.D0) chiDeg = 1.D0
      RGBd(1) = 1.D0 - chiDeg 
      RGBd(3) = abs(etaDeg - etaMin) / (etaMax-etaMin)
      RGBd(2) = (1.D0 - RGBd(3)) * chiDeg
      RGBd(3) = RGBd(3) * chiDeg
      RGBd = sqrt(RGBd) 
      ma = maxval(RGBd)
      RGBd = RGBd / ma 
      if (trim(self%ipf_mode).eq.'OPC') then 
        hsl = clr%rgb2hsl(RGBd)
        RGBd = clr%sph2rgb(hsl)
      end if 
      if (trim(self%ipf_mode).eq.'PUC') then 
        hsl = clr%rgb2hsl(RGBd)
        if (hsl(2).gt.0.5D0) then 
          RGBd = clr%sphere2rgb(hsl(2), hsl(0), w0=.FALSE. )
        else
          RGBd = clr%sphere2rgb(hsl(2), hsl(0), w0=.TRUE.)
        end if 
      end if 
      RGB = int(RGBd * 255)
      exit findloop
    end if 
  end do findloop
end if 

end function get_ipf_RGB_

!--------------------------------------------------------------------------
subroutine get_IPFMap_( self, EMsoft, sampleDir, Orientations, sym, cDir ) 
!DEC$ ATTRIBUTES DLLEXPORT :: get_IPFMap_
!! author: MDG
!! version: 1.0
!! date: 09/03/21
!!
!! This takes an orientation list in quaternion representation and turns it into 
!! an inverse pole figure map (using one of a few color schemes) that is then stored 
!! in a .tiff file.
!!
!! This routine assumes that the IPF map size and filename as well as the Laue Class 
!! have been set by the calling program using the appropriate methods; the number of orientations
!! in the Orientations array must equal the product self%ipf_wd * self%ipf_ht.

use mod_EMsoft
use mod_quaternions
use mod_OMPsupport
use omp_lib
use mod_image
use mod_io
use mod_colorspace
use ISO_C_BINDING
use, intrinsic :: iso_fortran_env

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)              :: self
type(EMsoft_T), INTENT(INOUT)               :: EMsoft
integer(kind=irg), INTENT(IN)               :: sampleDir(3)
type(QuaternionArray_T), INTENT(IN)         :: Orientations
type(QuaternionArray_T), INTENT(INOUT)      :: sym
logical,INTENT(IN),OPTIONAL                 :: cDir

integer(kind=irg),allocatable               :: IPFmap(:,:,:)
integer(kind=irg)                           :: ix, iy, iq, TID, RGB(0:2) 
type(Quaternion_T)                          :: qu
type(IO_T)                                  :: Message
type(colorspace_T)                          :: clr

character(fnlen)                            :: fname, TIFF_filename
real(kind=dbl)                              :: sDir(3)

! declare variables for use in object oriented image module
integer                                     :: iostat, io_int(2)
character(len=128)                          :: iomsg
logical                                     :: isInteger, OPC, PUC
type(image_t)                               :: im
integer(int8), allocatable                  :: TIFF_image(:,:)
integer                                     :: dim2(2), Pm
integer(c_int32_t)                          :: result

allocate(IPFmap(3, self%ipf_wd, self%ipf_ht))

! make sure sampleDir is normalized 
sDir = dble(sampleDir)/NORM2(dble(sampleDir))

if (.not.present(cDir)) then
  io_int = (/ self%ipf_wd, self%ipf_ht /)
  call Message%printMessage(' ')
  call Message%WriteValue(' Generating IPF map of dimensions : ', io_int,2,frm="(I5,' x ',I5)")
end if 

Pm = sym%getQnumber()

OPC = .FALSE.
if (trim(self%ipf_mode).eq.'OPC') OPC = .TRUE.

PUC = .FALSE.
if (trim(self%ipf_mode).eq.'PUC') PUC = .TRUE.

! io_int(1) = Pm
! call Message%WriteValue(' # symops = ', io_int, 1) 

! io_int(1) = self%get_ipf_nthreads_()
! call Message%WriteValue(' # threads = ', io_int, 1)

if (present(cDir)) then  ! don't do a parallel run ...

  if (OPC.eqv..TRUE.) clr = colorspace_T()
  if (PUC.eqv..TRUE.) clr = colorspace_T( Nfold = 4 )

  do iy = 1, self%ipf_ht 
    do ix = 1, self%ipf_wd 
      iq = (iy-1) * self%ipf_wd + ix 
      qu = Orientations%getQuatfromArray(iq)
      call qu%quat_pos()
      if ((OPC.eqv..TRUE.).or.(PUC.eqv..TRUE.)) then 
        IPFmap(1:3, ix, iy) = self%get_ipf_RGB_( sDir, qu, sym, Pm, clr )
      else 
        IPFmap(1:3, ix, iy) = self%get_ipf_RGB_( sDir, qu, sym, Pm )
      end if 
    end do 
  end do
else
  call OMP_SET_NUM_THREADS( self%get_ipf_nthreads_() )

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID, iy, iq, qu, clr)

  TID = OMP_GET_THREAD_NUM()

  if (OPC.eqv..TRUE.) clr = colorspace_T()
  if (PUC.eqv..TRUE.) clr = colorspace_T( Nfold = 4 )

  !$OMP DO SCHEDULE(DYNAMIC)
    do iy = 1, self%ipf_ht 
      do ix = 1, self%ipf_wd 
        iq = (iy-1) * self%ipf_wd + ix 
        qu = Orientations%getQuatfromArray(iq)
        call qu%quat_pos()
        if ((OPC.eqv..TRUE.).or.(PUC.eqv..TRUE.)) then 
          IPFmap(1:3, ix, iy) = self%get_ipf_RGB_( sDir, qu, sym, Pm, clr )
        else 
          IPFmap(1:3, ix, iy) = self%get_ipf_RGB_( sDir, qu, sym, Pm )
        end if 
      end do 
    end do
  !$OMP END DO

  !$OMP END PARALLEL
end if

! finally, store the IPFmap in a color TIFF file.
if (present(cDir)) then 
  fname = trim(self%get_ipf_filename_())
else 
  fname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(self%get_ipf_filename_())
end if
TIFF_filename = trim(fname)

! allocate memory for a color image; each pixel has 3 bytes (RGB)
allocate(TIFF_image(3*self%ipf_wd,self%ipf_ht))
TIFF_image = reshape( IPFmap, (/ 3*self%ipf_wd, self%ipf_ht /) )

! set up the image_t structure
im = image_t(TIFF_image)
im%dims = (/ self%ipf_wd, self%ipf_ht /)
im%samplesPerPixel = 3
if(im%empty()) call Message%printMessage("EMdpmerge: failed to convert array to rgb image")

! create the file
call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage(" Failed to write image to file : "//iomsg)
else
  call Message%printMessage(' IPF map written to '//trim(TIFF_filename),"(A)")
end if

deallocate(TIFF_image, IPFmap)

end subroutine get_IPFMap_

end module mod_IPFsupport
