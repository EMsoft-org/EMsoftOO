! ###################################################################
! Copyright (c) 2016-2020, Marc De Graef Research Group/Carnegie Mellon University
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

module MODQuaternionsTest
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/05/20
  !!
  !! perform a series of unit tests on the quaternion module 
  !!
  !! TODO: add test for quaternion slerp

contains 

subroutine MODQuaternionsExecuteTest(res) &
           bind(c, name='MODQuaternionsExecuteTest')    ! this routine is callable from a C/C++ program
!DEC$ ATTRIBUTES DLLEXPORT :: MODQuaternionsExecuteTest

use,INTRINSIC :: ISO_C_BINDING
use mod_kinds
use mod_quaternions

IMPLICIT NONE

integer(C_INT32_T),INTENT(OUT)  :: res

type(Quaternion_T)      :: u, v, w
real(kind=dbl)          :: x, a=2.D0, diffd
real(kind=sgl)          :: y, b=2.0, diff

! threshold values 
real(kind=dbl),parameter:: epsd = 1.0D-12
real(kind=sgl),parameter:: eps  = 1.0E-7

!===================================================
! set the reference values (verified with Mathematica scripts)
type(Quaternion_T) :: resdzero 
type(Quaternion_T) :: resdupv 
type(Quaternion_T) :: resdumv 
type(Quaternion_T) :: resduta 
type(Quaternion_T) :: resdutv 
type(Quaternion_T) :: resduc 
type(Quaternion_T) :: resdudv 
real(kind=dbl), parameter  :: absvald = 5.477225575051661D0
real(kind=dbl),parameter   :: ipd = 70.D0
real(kind=dbl),parameter   :: qand = 1.570796326794897D0
real(kind=dbl)             :: vecd(3) = (/ 1.D0, 0.D0, 0.D0 /), rvecd(3)

type(Quaternion_T) :: resszero 
type(Quaternion_T) :: ressupv 
type(Quaternion_T) :: ressumv 
type(Quaternion_T) :: ressuta 
type(Quaternion_T) :: ressutv 
type(Quaternion_T) :: ressuc 
type(Quaternion_T) :: ressudv 
real(kind=sgl), parameter  :: absvals = 5.4772256
real(kind=sgl),parameter   :: ips = 70.0
real(kind=sgl),parameter   :: qans = 1.5707963
real(kind=sgl)             :: vecs(3) = (/ 1.0, 0.0, 0.0 /), rvecs(3)

resdzero = Quaternion_T( qd = (/ 0.D0, 0.D0, 0.D0, 0.D0/) )
resdupv  = Quaternion_T( qd = (/ 6.D0, 8.D0, 10.D0, 12.D0/) )
resdumv  = Quaternion_T( qd = (/-4.D0, -4.D0, -4.D0, -4.D0/) )
resduta  = Quaternion_T( qd = (/2.D0, 4.D0, 6.D0, 8.D0/) )
resdutv  = Quaternion_T( qd = (/-60.D0, 12.D0, 30.D0, 24.D0/) )
resduc   = Quaternion_T( qd = (/1.D0,-2.D0,-3.D0,-4.D0/) )
resdudv  = Quaternion_T( qd = (/0.40229885057471D0, 0.045977011494252D0, 0.D0, 0.09195402298850575D0/) )

resszero = Quaternion_T( q = (/ 0.0, 0.0, 0.0, 0.0/))
ressupv  = Quaternion_T( q = (/6.0, 8.0, 10.0, 12.0/) )
ressumv  = Quaternion_T( q = (/-4.0, -4.0, -4.0, -4.0/) )
ressuta  = Quaternion_T( q = (/2.0, 4.0, 6.0, 8.0/) )
ressutv  = Quaternion_T( q = (/-60.0, 12.0, 30.0, 24.0/) )
ressuc   = Quaternion_T( q = (/1.0,-2.0,-3.0,-4.0/) )
ressudv  = Quaternion_T( q = (/0.40229885, 0.045977011, 0.0, 0.091954023/) )

! initialize the error identifier to zero (should remain zero upon successful exit)
res = 0

!===================================================
!=============Double Precision Tests================
!===================================================

!===================================================
! initialize zero quaternion 
u = Quaternion_T( qd = (/0.D0, 0.D0, 0.D0, 0.D0/) )
diffd = cabs(u)
if (diffd.gt.epsd) then 
  res = 1
  write (*,"('double precision zero initialization test failed = ',D18.10)") diff
  return
end if

if (.not.(u%quatsequal(resdzero))) then 
  res = 2
  write (*,"('double precision zero comparison test failed = ',D18.10)") diff
  return
end if

! arithmetic tests 
! ! u = Quaternion_T( qd=(/ 1.D0,2.D0,3.D0,4.D0/) )
! v = Quaternion_T( qd=(/ 5.D0,6.D0,7.D0,8.D0/) )
u = Quaternion_T( qd=(/ 1.0_dbl,2.0_dbl,3.0_dbl,4.0_dbl/) )
v = Quaternion_T( qd=(/ 5.0_dbl,6.0_dbl,7.0_dbl,8.0_dbl/) )

w = u+v
if (.not.(w%quatsequal(resdupv))) then 
  res = 3
  write (*,"('double precision addition test failed = ')") 
  return
end if

w = u-v
if (.not.(w%quatsequal(resdumv))) then 
  res = 4
  write (*,"('double precision subtraction test failed = ')")
  return
end if

w = u*a
if (.not.(w%quatsequal(resduta))) then 
  res = 5
  write (*,"('double precision scalar multiplication test failed = ')")
  return
end if

w = u*v
if (.not.(w%quatsequal(resdutv))) then 
  res = 6
  write (*,"('double precision quaternion multiplication test failed = ')") 
  return
end if

w = u/v
if (.not.(w%quatsequal(resdudv))) then 
  res = 7
  write (*,"('double precision division test failed = ')")
  return
end if

w = conjg(u)
if (.not.(w%quatsequal(resduc))) then 
  res = 8
  write (*,"('double precision conjugation test failed = ')") 
  return
end if

diffd = cabs(u) - absvald
if (diffd.gt.epsd) then 
  res = 9
  write (*,"('double precision norm test failed = ',D18.10)") diffd
  return
end if

diffd = u%quat_innerproduct(v) - ipd 
if (diffd.gt.epsd) then 
  res = 10
  write (*,"('double precision innerproduct test failed = ',D18.10)") diffd
  return
end if

diffd = u%quat_angle(v) - qand
if (diffd.gt.epsd) then 
  res = 11
  write (*,"('double precision quaternion angle test failed = ',D18.10)") diffd
  return
end if

! from here on we work with unit quaternions
call u%quat_normalize()
call v%quat_normalize()

diffd = cabs(u) - 1.D0 
if (diffd.gt.epsd) then 
  res = 12
  write (*,"('double precision normalization test failed = ',D18.10)") diffd
  return
end if

! make u a rotation of 120° around the [111] axis 
u = Quaternion_T( qd=(/ 0.5_dbl, 0.5_dbl, 0.5_dbl, 0.5_dbl /) )
rvecd = u%quat_Lp(vecd)
diffd = sum( rvecd - (/ 0.D0, 1.D0, 0.D0/) )
if (diffd.gt.epsd) then 
  res = 13
  write (*,"('double precision vector rotation test failed = ',D18.10)") diffd
  return
end if



!===================================================
!=============Single Precision Tests================
!===================================================

!===================================================
! initialize zero quaternion 
u = Quaternion_T( q = (/0.0, 0.0, 0.0, 0.0/) )
diff = cabs( u )
if (diff.gt.eps) then 
  res = 20
  write (*,"('single precision zero initialization test failed = ',F12.8)") diff
  return
end if

if (.not.(u%quatsequal(resszero))) then 
  res = 21
  write (*,"('single precision zero comparison test failed = ')") 
  return
end if

! arithmetic tests 
u = Quaternion_T( q=(/1.0,2.0,3.0,4.0/) )
v = Quaternion_T( q=(/5.0,6.0,7.0,8.0/) )

w = u+v
if (.not.(w%quatsequal(ressupv))) then 
  res = 22
  write (*,"('single precision addition test failed = ')") 
  return
end if

w = u-v
if (.not.(w%quatsequal(ressumv))) then 
  res = 23
  write (*,"('single precision subtraction test failed = ')") 
  return
end if

w = u*b
if (.not.(w%quatsequal(ressuta))) then 
  res = 24
  write (*,"('single precision scalar multiplication test failed = ')")
  return
end if

w = u*v
if (.not.(w%quatsequal(ressutv))) then 
  res = 25 
  write (*,"('single precision quaternion multiplication test failed = ')")
  return
end if

w = u/v
if (.not.(w%quatsequal(ressudv))) then 
  res = 26
  write (*,"('single precision division test failed = ')") 
  return
end if

w = conjg(u)
if (.not.(w%quatsequal(ressuc))) then 
  res = 27
  write (*,"('single precision conjugation test failed = ')") 
  return
end if

diff = cabs(u) - absvals
if (diff.gt.eps) then 
  res = 28
  write (*,"('single precision norm test failed = ',F12.8)") diff
  return
end if

diff = u%quat_innerproduct(v) - ips 
if (diff.gt.eps) then 
  res = 29
  write (*,"('single precision innerproduct test failed = ',F12.8)") diff
  return
end if

diff = u%quat_angle(v) - qans
if (diff.gt.eps) then 
  res = 30
  write (*,"('single precision quaternion angle test failed = ',F12.8)") diff
  return
end if

! from here on we work with unit quaternions
call u%quat_normalize()
call v%quat_normalize()

diff = cabs(u) - 1.0 
if (diff.gt.eps) then 
  res = 31
  write (*,"('single precision normalization test failed = ',F12.8)") diff
  return
end if

! make u a rotation of 120° around the [111] axis 
u = Quaternion_T( q=(/ 0.5_sgl, 0.5_sgl, 0.5_sgl, 0.5_sgl /) )
rvecs = u%quat_Lp(vecs)
diff = sum( rvecs - (/ 0.0, 1.0, 0.0/) )
if (diff.gt.eps) then 
  res = 32
  write (*,"('single precision vector rotation test failed = ',D18.10)") diff
  return
end if



end subroutine MODQuaternionsExecuteTest


end module MODQuaternionsTest
