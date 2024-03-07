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

module mod_Ylm
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/07/23
  !!
  !! class definition for the Ylm Spherical Harmonics
  !!
  !! much of this is based on https://arxiv.org/pdf/1410.1748.pdf

use mod_kinds
use mod_global

IMPLICIT NONE 

type,public :: SH_Coefficients
    integer(kind=irg)                   :: maxL
    real(kind=dbl),allocatable          :: Alm(:)
    real(kind=dbl),allocatable          :: Blm(:)
    real(kind=dbl),allocatable          :: Plm(:)
    real(kind=dbl),allocatable          :: Ylm(:)
    real(kind=dbl)                      :: lastPhi
    real(kind=dbl)                      :: lastTheta
end type SH_Coefficients

! only used in SH_buildHarmonic routine, so commented out for now...
! type,public :: SH_Mode 
!   integer  :: l
!   integer  :: m 
!   complex  :: weight
! end type SH_Mode

! class definition
type, public :: Ylm_T
private 
  type(SH_Coefficients) :: SHcoeff
 
contains
private 
  procedure, pass(self) :: SH_PT_
  procedure, pass(self) :: SH_YR_
  procedure, pass(self) :: SH_initialize_
  procedure, pass(self) :: SH_ComputeP_
  procedure, pass(self) :: SH_ComputeY_
  procedure, pass(self) :: SH_ComputeYlm_
  ! procedure, pass(self) :: SH_buildHarmonic

  generic, public :: SH_PT => SH_PT_
  generic, public :: SH_YR => SH_YR_
  generic, public :: SH_initialize => SH_initialize_
  generic, public :: SH_ComputeP => SH_ComputeP_
  generic, public :: SH_ComputeY => SH_ComputeY_
  generic, public :: SH_ComputeYlm => SH_ComputeYlm_
  ! generic, public :: SH_buildHarmonic => SH_buildHarmonic_

end type Ylm_T

! the constructor routine for this class 
interface Ylm_T
  module procedure Ylm_constructor
end interface Ylm_T

contains

!--------------------------------------------------------------------------
type(Ylm_T) function Ylm_constructor( ) result(Ylm)
!! author: MDG 
!! version: 1.0 
!! date: 12/07/23
!!
!! constructor for the Ylm_T Class; 
 
IMPLICIT NONE

end function Ylm_constructor

!--------------------------------------------------------------------------
subroutine Ylm_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 12/07/23
!!
!! destructor for the Ylm_T Class
 
IMPLICIT NONE

type(Ylm_T), INTENT(INOUT)  :: self 

call reportDestructor('Ylm_T')

end subroutine Ylm_destructor

!--------------------------------------------------------------------------
recursive function SH_PT_(self, l, m) result(PT)
!DEC$ ATTRIBUTES DLLEXPORT :: SH_PT_
!! author: WL/MDG 
!! version: 1.0 
!! date: 12/07/23
!!
!! combine l and m into an array index for Spherical Harmonics

IMPLICIT NONE

class(Ylm_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: l
integer(kind=irg), INTENT(IN)   :: m
integer(kind=irg)               :: PT

PT = m + l * (l+1) / 2

end function SH_PT_

!--------------------------------------------------------------------------
recursive function SH_YR_(self, l, m) result(PT)
!DEC$ ATTRIBUTES DLLEXPORT :: SH_YR_
!! author: WL/MDG 
!! version: 1.0 
!! date: 12/07/23
!!
!! combine l and m into an array index for Spherical Harmonics

IMPLICIT NONE

class(Ylm_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: l
integer(kind=irg), INTENT(IN)   :: m
integer(kind=irg)               :: PT

PT = m + l * (l+1) 

end function SH_YR_

!--------------------------------------------------------------------------
recursive subroutine SH_initialize_(self, maxL)
!DEC$ ATTRIBUTES DLLEXPORT :: SH_initialize_
!! author: WL/MDG 
!! version: 1.0 
!! date: 12/07/23
!!
!! initialize Spherical Harmonic auxiliary arrays

IMPLICIT NONE

class(Ylm_T), INTENT(INOUT)           :: self
integer(kind=irg), INTENT(IN)         :: maxL

integer(kind=irg)                     :: l, m, ls, lm1s, ms

if (allocated(self%SHcoeff%Alm).eqv..TRUE.) deallocate(self%SHcoeff%Alm)
if (allocated(self%SHcoeff%Blm).eqv..TRUE.) deallocate(self%SHcoeff%Blm)
if (allocated(self%SHcoeff%Plm).eqv..TRUE.) deallocate(self%SHcoeff%Plm)
if (allocated(self%SHcoeff%Ylm).eqv..TRUE.) deallocate(self%SHcoeff%Ylm)
 
allocate(self%SHcoeff%Alm(0:self%SH_PT_(maxL,maxL)+1))
allocate(self%SHcoeff%Blm(0:self%SH_PT_(maxL,maxL)+1))
allocate(self%SHcoeff%Plm(0:self%SH_PT_(maxL,maxL)+1))
allocate(self%SHcoeff%Ylm(0:self%SH_YR_(maxL,maxL)+1))

self%SHcoeff%Alm = 0.D0
self%SHcoeff%Blm = 0.D0
self%SHcoeff%Plm = 0.D0
self%SHcoeff%Ylm = 0.D0

do l=2,maxL
  ls = l*l
  lm1s = (l-1)*(l-1)
  do m=0,l-2
    ms = m*m
    self%SHcoeff%Alm(self%SH_PT_(l,m)) = sqrt((4.D0*dble(ls)-1.D0)/dble(ls-ms))
    self%SHcoeff%Blm(self%SH_PT_(l,m)) = -sqrt((dble(lm1s-ms))/dble(4*lm1s-1))
  end do
end do

self%SHcoeff%maxL = maxL

! set these to random large numbers
self%SHcoeff%lastTheta = 213424324.D0
self%SHcoeff%lastPhi = 54325154325.D0

end subroutine SH_initialize_

!--------------------------------------------------------------------------
recursive subroutine SH_ComputeP_(self, x)
!DEC$ ATTRIBUTES DLLEXPORT :: SH_ComputeP_
!! author: WL/MDG 
!! version: 1.0 
!! date: 12/07/23
!!
!! compute a set of Plm for a given value of x

IMPLICIT NONE

class(Ylm_T), INTENT(INOUT)           :: self
real(kind=dbl), INTENT(IN)            :: x

real(kind=dbl), parameter             :: sqrt3 = 1.7320508075688772935D0, sqrt3div2 = -1.2247448713915890491D0
real(kind=dbl)                        :: temp, sintheta
integer(kind=irg)                     :: L, ll, m, PTlm

if (x.ne.self%SHcoeff%lastTheta) then 
    self%SHcoeff%lastTheta = x

    L = self%SHcoeff%maxL
    temp = 0.39894228040143267794D0
    sintheta = sqrt(1.D0-x*x)

    self%SHcoeff%Plm(self%SH_PT_(0,0)) = temp

    if (L.gt.0) then 
      self%SHcoeff%Plm(self%SH_PT_(1,0)) = x * sqrt3 * temp
      temp = sqrt3div2*sintheta*temp
      self%SHcoeff%Plm(self%SH_PT_(1,1)) = temp

      do ll=2,L
        do m=0,l-2
          PTlm = self%SH_PT_(ll,m)
          self%SHcoeff%Plm(PTlm) = self%SHcoeff%Alm(PTlm) * ( x * self%SHcoeff%Plm(self%SH_PT_(ll-1,m)) &
                                                    + self%SHcoeff%Blm(PTlm) * self%SHcoeff%Plm(self%SH_PT_(ll-2,m)) )
        end do
        self%SHcoeff%Plm(self%SH_PT_(ll,ll-1)) = x * sqrt(dble(2*(ll-1)+3)) * temp
        temp = -sqrt(1.D0+0.5D0/dble(ll)) * sintheta * temp
        self%SHcoeff%Plm(self%SH_PT_(ll,ll)) = temp
      end do
    end if 
end if

end subroutine SH_ComputeP_

!--------------------------------------------------------------------------
recursive subroutine SH_ComputeY_(self, theta, phi)
!DEC$ ATTRIBUTES DLLEXPORT :: SH_ComputeY_
!! author: WL/MDG 
!! version: 1.0 
!! date: 12/07/23
!!
!! precompute a set of Ylm for a given value of x

IMPLICIT NONE

class(Ylm_T), INTENT(INOUT)           :: self
real(kind=dbl), INTENT(IN)            :: theta
real(kind=dbl), INTENT(IN)            :: phi

real(kind=dbl)                        :: c1, c2, s1, s2, tc, s, c, tt, ctheta
integer(kind=irg)                     :: L, ll, m 

! this routine assumes that phi is the fastest changing variable !!!

if (phi.ne.self%SHcoeff%lastPhi) then  ! has phi changed ?  if so, then recompute the arrays
    ctheta = cos(theta)
    if (ctheta.ne.self%SHcoeff%lastTheta) then  ! has theta changed ? 
        call self%SH_ComputeP_(ctheta)
    end if

    L = self%SHcoeff%maxL

    m = 0
    do ll=0,L
      self%SHcoeff%Ylm(self%SH_YR_(ll,m)) = self%SHcoeff%Plm(self%SH_PT_(ll,m)) * sqrt(0.5D0)
    end do 

    c1 = 1.D0
    c2 = cos(phi)
    s1 = 0.D0
    s2 = -sin(phi)
    tc = 2.D0 * c2
    do m=1,L
      s = tc * s1 - s2
      c = tc * c1 - c2
      s2 = s1
      s1 = s
      c2 = c1
      c1 = c
      do ll=m,L 
        tt = self%SHcoeff%Plm(self%SH_PT_(ll,m))
        self%SHcoeff%Ylm(self%SH_YR_(ll,-m)) = tt * s
        self%SHcoeff%Ylm(self%SH_YR_(ll,m)) = tt * c
      end do
    end do
    self%SHcoeff%lastPhi = phi
end if

end subroutine SH_ComputeY_

!--------------------------------------------------------------------------
recursive function SH_ComputeYlm_(self, l, m, theta, phi) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: SH_ComputeYlm_
!! author: WL/MDG 
!! version: 1.0 
!! date: 12/07/23
!!
!! compute a set of Ylm for a given value of l, m, theta, phi

IMPLICIT NONE

class(Ylm_T), INTENT(INOUT)           :: self
integer(kind=irg), INTENT(IN)         :: l
integer(kind=irg), INTENT(IN)         :: m
real(kind=dbl), INTENT(IN)            :: theta
real(kind=dbl), INTENT(IN)            :: phi
complex(kind=dbl)                     :: res

real(kind=dbl),parameter              :: s = 0.707106781186547524D0   ! sqrt(1/2)

call self%SH_ComputeY( theta, phi)

if (m.lt.0) then
    res = cmplx(self%SHcoeff%Ylm(self%SH_YR_(l,-m))*s,-self%SHcoeff%Ylm(self%SH_YR_(l,m))*s)
    if (mod(m,2).eq.1) res = -res 
else if (m.eq.0) then
    res = cmplx(self%SHcoeff%Ylm(self%SH_YR_(l,0)),0.D0)
else
    res = cmplx(self%SHcoeff%Ylm(self%SH_YR_(l,m))*s,self%SHcoeff%Ylm(self%SH_YR_(l,-m))*s)
end if

end function SH_ComputeYlm_

!--------------------------------------------------------------------------
! this routine is apparently not used anywhere...so we leave it commented out
! recursive subroutine SH_buildHarmonic_(SHcoeff, dim, modes, data) 
! !DEC$ ATTRIBUTES DLLEXPORT :: SH_buildHarmonic_
! !! author: WL/MDG 
! !! version: 1.0 
! !! date: 12/07/23
! !!
! !! build a set of Spherical Harmonics

! use mod_Lambert

! IMPLICIT NONE

! type(SH_Coefficients),INTENT(INOUT)   :: SHcoeff
! integer(kind=irg),INTENT(IN)          :: dim
! type(SH_Mode),INTENT(IN)              :: modes(:)
! real(kind=dbl),INTENT(INOUT)          :: data(2*dim*dim)

! integer(kind=irg)                     :: i, j, l, m, sz, ierr
! real(kind=dbl)                        :: xyz(3), XY(2), Phi, thetaN, thetaS
! complex(kind=dbl)                     :: nh, sh 

! ! determine the maximum mode
! SHcoeff%maxL = 0
! sz = size(modes)
! do m=1,sz 
!   if (modes(m)%l.gt.SHcoeff%maxL) SHcoeff%maxL=modes(m)%l
! end do

! ! initialize the spherical harmonic coefficient arrays
! call SH_initialize(SHcoeff, SHcoeff%maxL)

! ! loop over the Lambert square
! do j=0,dim-1
!   XY(2) = dble(j)/dble(dim-1)
!   do i=0,dim-1
!     XY(1) = dble(i)/dble(dim-1)
!     nh = cmplx(0.D0,0.D0)
!     sh = cmplx(0.D0,0.D0)

! ! transform to a unit vector on the sphere
!     xyz = LambertSquareToSphere(XY, ierr)
! ! we may need to convert to the Legendre grid at this point ...

! ! get the angles
!     Phi = atan2(xyz(2),xyz(1))
!     thetaN = acos(xyz(3))
!     thetaS = acos(-xyz(3))

! ! compute signal in NH
!     call SH_ComputeY(SHcoeff, thetaN, Phi)
!     do m=1,sz
!       nh = nh + SH_ComputeYlm(SHcoeff, modes(m)%l, modes(m)%m, thetaN, Phi) * modes(m)%weight
!       if (modes(m)%m.ne.0) then
!         if (mod(modes(m)%m,2).eq.0) then 
!           nh = nh + SH_ComputeYlm(SHcoeff, modes(m)%l, -modes(m)%m, thetaN, Phi) * conjg(modes(m)%weight)
!         else
!           nh = nh - SH_ComputeYlm(SHcoeff, modes(m)%l, -modes(m)%m, thetaN, Phi) * conjg(modes(m)%weight)
!         end if
!       end if
!     end do

! ! compute signal in SH
!     call SH_ComputeY(SHcoeff, thetaS, Phi)
!     do m=1,sz
!       sh = sh + SH_ComputeYlm(SHcoeff, modes(m)%l, modes(m)%m, thetaN, Phi) * modes(m)%weight
!       if (modes(m)%m.ne.0) then
!         if (mod(modes(m)%m,2).eq.0) then 
!           sh = sh + SH_ComputeYlm(SHcoeff, modes(m)%l, -modes(m)%m, thetaN, Phi) * conjg(modes(m)%weight)
!         else
!           sh = sh - SH_ComputeYlm(SHcoeff, modes(m)%l, -modes(m)%m, thetaN, Phi) * conjg(modes(m)%weight)
!         end if
!       end if
!     end do

!     data(j*dim+i) = real(nh)
!     data(dim*dim + j*dim+i) = real(sh)

!   end do 
! end do

! end subroutine SH_buildHarmonic_



end module mod_Ylm