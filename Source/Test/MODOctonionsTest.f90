! ###################################################################
! Copyright (c) 2016-2022, Marc De Graef Research Group/Carnegie Mellon University
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

module MODOctonionsTest
  !! author: MDG 
  !! version: 1.0 
  !! date: 10/17/22
  !!
  !! perform a series of unit tests on the octonion module 
  !!

contains 

subroutine MODOctonionsExecuteTest(res) &
           bind(c, name='MODOctonionsExecuteTest')    ! this routine is callable from a C/C++ program
!DEC$ ATTRIBUTES DLLEXPORT :: MODOctonionsExecuteTest

use,INTRINSIC :: ISO_C_BINDING
use mod_kinds
use mod_octonions

IMPLICIT NONE

integer(C_INT32_T),INTENT(OUT)  :: res

type(Octonion_T)        :: a, b, c, d, u
real(kind=dbl)          :: diffd
real(kind=sgl)          :: diff

! threshold values 
real(kind=dbl),parameter:: epsd = 1.0D-12
real(kind=sgl),parameter:: eps  = 1.0E-7

! various parameters
integer(kind=irg)       :: i, errcnt 

type(Octonion_T) :: resdsum
type(Octonion_T) :: resdsub
type(Octonion_T) :: resdmult 
type(Octonion_T) :: resdsmult 
type(Octonion_T) :: resddiv 
type(Octonion_T) :: resdconjg 
type(Octonion_T) :: resdinv 
type(Octonion_T) :: resdzero
real(kind=dbl)   :: resdabs

type(Octonion_T) :: resssum
type(Octonion_T) :: resssub
type(Octonion_T) :: ressmult 
type(Octonion_T) :: resssmult 
type(Octonion_T) :: ressdiv 
type(Octonion_T) :: ressconjg 
type(Octonion_T) :: ressinv 
type(Octonion_T) :: resszero
real(kind=sgl)   :: ressabs
!===================================================
! set the reference values (verified against values on <https://pypi.org/project/pyoctonion/#description>)

! correct answers for the Octonion_T class
resdzero = Octonion_T()
resdsum = Octonion_T( od = (/  2.0D0,   5.0D0,   8.0D0,  11.0D0,  14.0D0,   8.0D0,  11.0D0,  14.0D0 /) )
resdsub = Octonion_T( od = (/  0.0D0,   -1.0D0,   -2.0D0,   -3.0D0,   -4.0D0,    4.0D0,    3.0D0,    2.0D0 /) )
resdsmult = Octonion_T( od = (/ 1.41421356237310D0, 2.82842712474619D0, 4.24264068711929D0, 5.65685424949238D0, &
                                7.07106781186548D0, 8.48528137423857D0, 9.89949493661167D0,11.31370849898476D0 /) )
resdmult = Octonion_T( od = (/-181.0D0,  -48.0D0,  -17.0D0,  -40.0D0,   83.0D0,    0.0D0,   35.0D0,    4.0D0 /) )
resddiv = Octonion_T( od = (/  0.82805429864253D0,    0.23529411764706D0,    0.10407239819005D0,    0.21719457013575D0,  &
                              -0.33031674208145D0,    0.05429864253394D0,   -0.09502262443439D0,    0.05429864253394D0 /) )
resdabs = 14.282856857085701D0
resdconjg = Octonion_T( od = (/  1.0D0,   -2.0D0,   -3.0D0,   -4.0D0,   -5.0D0,   -6.0D0,   -7.0D0,   -8.0D0 /) )
resdinv = Octonion_T( od = (/  0.00490196078431D0,  -0.00980392156863D0,  -0.01470588235294D0,  -0.01960784313725D0,  &
                              -0.02450980392157D0,  -0.02941176470588D0,  -0.03431372549020D0, -0.03921568627451D0 /) )

resszero = Octonion_T(smode = 's')
resssum = Octonion_T( o = (/  2.00,   5.00,   8.00,  11.00,  14.00,   8.00,  11.00,  14.00 /) )
resssub = Octonion_T( o = (/  0.00,   -1.00,   -2.00,   -3.00,   -4.00,    4.00,    3.00,    2.00 /) )
resssmult = Octonion_T( o = (/ 1.414214, 2.828427, 4.242640, 5.656854, 7.071068, 8.485281, 9.899495, 11.313708 /) )
ressmult = Octonion_T( o = (/-181.00,  -48.00,  -17.00,  -40.00,   83.00,    0.00,   35.00,    4.00 /) )
ressdiv = Octonion_T( o = (/  0.828054,   0.235294,   0.104072,   0.217195,  -0.330317,   0.054299,  -0.095023,   0.054299 /) )
ressabs  = 14.2828569
ressconjg = Octonion_T( o = (/  1.0,   -2.0,   -3.0,   -4.0,   -5.0,   -6.0,   -7.0,   -8.0 /) )
ressinv = Octonion_T( o = (/  0.004902,  -0.009804,  -0.014706,  -0.019608,  -0.024510,  -0.029412,  -0.034314,  -0.039216 /) )

! correct answers for the QuaternionArray_T class 
! type(QuaternionArray_T) :: resarrdadd

! resarrdadd = QuaternionArray_T( n=4, )

! initialize the error identifier to zero (should remain zero upon successful exit)
res = 0

!===================================================
!=============Double Precision Tests================
!===================================================

!===================================================
! initialize zero quaternion 
u = Octonion_T( od = (/ 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /) )
diffd = cabs(u)
if (diffd.gt.epsd) then 
  res = 1
  write (*,"('double precision zero initialization test failed = ',D18.10)") diff
  return
end if

if (.not.(u%octsequal(resdzero))) then 
  res = 2
  write (*,"('double precision zero comparison test failed = ',D18.10)") diff
  return
end if

! arithmetic tests 
a = Octonion_T( od = (/ 1.D0, 2.D0, 3.D0, 4.D0, 5.D0, 6.D0, 7.D0, 8.D0 /) )
b = Octonion_T( od = (/ 1.D0, 3.D0, 5.D0, 7.D0, 9.D0, 2.D0, 4.D0, 6.D0 /) )
c = Octonion_T( od = (/ 8.D0, 7.D0, 6.D0, 5.D0, 4.D0, 3.D0, 2.D0, 1.D0 /) )

d = a+b
if (.not.(d%octsequal(resdsum))) then 
  res = 3
  write (*,"('double precision addition test failed = ')") 
  return
end if

d = a-b
if (.not.(d%octsequal(resdsub))) then 
  res = 4
  write (*,"('double precision subtraction test failed = ')")
  return
end if

d = a*sqrt(2.D0)
if (.not.(d%octsequal(resdsmult))) then 
  res = 5
  write (*,"('double precision scalar multiplication test failed = ')")
  return
end if

d = a*b
if (.not.(d%octsequal(resdmult))) then 
  res = 6
  write (*,"('double precision quaternion multiplication test failed = ')") 
  return
end if

d = a/b
if (.not.(d%octsequal(resddiv))) then 
  res = 7
  write (*,"('double precision division test failed = ')")
  return
end if

d = conjg(a)
if (.not.(d%octsequal(resdconjg))) then 
  res = 8
  write (*,"('double precision conjugation test failed = ')") 
  return
end if

diffd = abs(cabs(a) - resdabs)
if (diffd.gt.epsd) then 
  res = 9
  write (*,"('double precision norm test failed = ',D18.10)") diffd
  return
end if

d = a%octinverse()
if (.not.(d%octsequal(resdinv))) then 
  res = 10
  write (*,"('double precision inverse test failed = ')") 
  return
end if

! from here on we work with unit quaternions
call a%octnormalize()
call a%oct_print(' regular normalization : ')
diffd = abs(cabs(a) - 1.D0)
if (diffd.gt.epsd) then 
  res = 11
  write (*,"('double precision normalization test failed = ',D18.10)") diffd
  return
end if

! make u a rotation of 120° around the [111] axis 
call a%set_GBmode(.TRUE.)
call a%octnormalize()
call a%oct_print(' GBOM normalization : ')
diffd = abs(cabs(a) - 1.D0)
if (diffd.gt.epsd) then 
  res = 12
  write (*,"('double precision normalization test failed = ',D18.10)") diffd
  return
end if
call a%set_GBmode(.FALSE.)

!===================================================
!=============Single Precision Tests================
!===================================================

!===================================================
! initialize zero quaternion 
u = Octonion_T( o = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /) )
diffd = cabs(u)
if (diffd.gt.eps) then 
  res = 13
  write (*,"('single precision zero initialization test failed = ',D18.10)") diff
  return
end if

if (.not.(u%octsequal(resszero))) then 
  res = 14
  write (*,"('single precision zero comparison test failed = ',D18.10)") diff
  return
end if

! arithmetic tests 
a = Octonion_T( o = (/ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 /) )
b = Octonion_T( o = (/ 1.0, 3.0, 5.0, 7.0, 9.0, 2.0, 4.0, 6.0 /) )
c = Octonion_T( o = (/ 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 /) )
d = Octonion_T( smode = 's' )

d = a+b
if (.not.(d%octsequal(resssum))) then 
  res = 15
  write (*,"('single precision addition test failed = ')") 
  return
end if

d = a-b
if (.not.(d%octsequal(resssub))) then 
  res = 16
  write (*,"('single precision subtraction test failed = ')")
  return
end if

d = a*sqrt(2.D0)
call d%oct_print('d = ')
call resssmult%oct_print('resssmult = ')
if (.not.(d%octsequal(resssmult))) then 
  res = 17
  write (*,"('single precision scalar multiplication test failed = ')")
  return
end if

d = a*b
if (.not.(d%octsequal(ressmult))) then 
  res = 18
  write (*,"('single precision quaternion multiplication test failed = ')") 
  return
end if

d = a/b
if (.not.(d%octsequal(ressdiv))) then 
  res = 19
  write (*,"('single precision division test failed = ')")
  return
end if

d = conjg(a)
if (.not.(d%octsequal(ressconjg))) then 
  res = 20
  write (*,"('single precision conjugation test failed = ')") 
  return
end if

diff = abs(cabs(a) - ressabs)
if (diff.gt.eps) then 
  res = 21
  write (*,"('single precision norm test failed = ',D18.10)") diffd
  return
end if

d = a%octinverse()
if (.not.(d%octsequal(ressinv))) then 
  res = 22
  write (*,"('single precision inverse test failed = ')") 
  return
end if

! from here on we work with unit quaternions
call a%octnormalize()
call a%oct_print(' regular normalization : ')
diff = abs(cabs(a) - 1.0)
if (diff.gt.eps) then 
  res = 23
  write (*,"('single precision normalization test failed = ',D18.10)") diffd
  return
end if

! make u a rotation of 120° around the [111] axis 
call a%set_GBmode(.TRUE.)
call a%octnormalize()
call a%oct_print(' GBOM normalization : ')
diff = abs(cabs(a) - 1.0)
if (diff.gt.eps) then 
  res = 24
  write (*,"('single precision normalization test failed = ',D18.10)") diffd
  return
end if
call a%set_GBmode(.FALSE.)




! !===================================================
! !==========Double Precision Array Tests=============
! !===================================================

! !===================================================
! ! declare two arrays with identical quaternions to simplify the testing ... 
! uA = QuaternionArray_T( n=4, qd = reshape( (/ 1.0_dbl,2.0_dbl,3.0_dbl,4.0_dbl, &
!                                               1.0_dbl,2.0_dbl,3.0_dbl,4.0_dbl, &
!                                               1.0_dbl,2.0_dbl,3.0_dbl,4.0_dbl, &
!                                               1.0_dbl,2.0_dbl,3.0_dbl,4.0_dbl /), (/ 4, 4 /) ) )

! vA = QuaternionArray_T( n=4, qd = reshape( (/ 5.0_dbl,6.0_dbl,7.0_dbl,8.0_dbl, &
!                                               5.0_dbl,6.0_dbl,7.0_dbl,8.0_dbl, &
!                                               5.0_dbl,6.0_dbl,7.0_dbl,8.0_dbl, &
!                                               5.0_dbl,6.0_dbl,7.0_dbl,8.0_dbl /), (/ 4, 4 /) ) )

! wA = uA + vA 
! errcnt = 0
! do i=1,4
!   u = wA%getQuatfromArray(i)
!   if (.not.(u%octsequal(resdupv))) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 40
!   write (*,"('double precision quaternion array addition test failed = ')")
!   return
! end if

! wA = uA - vA 
! errcnt = 0
! do i=1,4
!   u = wA%getQuatfromArray(i)
!   if (.not.(u%octsequal(resdumv))) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 41
!   write (*,"('double precision quaternion array subtraction test failed = ')")
!   return
! end if

! wA = uA * vA 
! errcnt = 0
! do i=1,4
!   u = wA%getQuatfromArray(i)
!   if (.not.(u%octsequal(resdutv))) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 42
!   write (*,"('double precision quaternion array multiplication test failed = ')")
!   return
! end if

! wA = uA / vA 
! errcnt = 0
! do i=1,4
!   u = wA%getQuatfromArray(i)
!   if (.not.(u%octsequal(resdudv))) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 43
!   write (*,"('double precision quaternion array division test failed = ')")
!   return
! end if

! wA = uA * a 
! errcnt = 0
! do i=1,4
!   u = wA%getQuatfromArray(i)
!   if (.not.(u%octsequal(resduta))) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 44
!   write (*,"('double precision quaternion array scalar multiplication test failed = ')")
!   return
! end if

! wA = conjg(uA) 
! errcnt = 0
! do i=1,4
!   u = wA%getQuatfromArray(i)
!   if (.not.(u%octsequal(resduc))) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 45
!   write (*,"('double precision quaternion array conjugation test failed = ')")
!   return
! end if

! dd = cabs(uA)
! errcnt = 0
! do i=1,4
!   if (abs(dd(i)-absvald).gt.epsd) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 46
!   write (*,"('double precision quaternion array abs test failed = ')")
!   return
! end if

! dd = uA%quat_innerproduct(vA)
! errcnt = 0
! do i=1,4
!   if (abs(dd(i)-ipd).gt.epsd) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 47
!   write (*,"('double precision quaternion array inner product test failed = ')")
!   return
! end if

! dd = uA%quat_angle(vA)
! errcnt = 0
! do i=1,4
!   if (abs(dd(i)-qand).gt.epsd) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 48
!   write (*,"('double precision quaternion array interquaternion angle test failed ')")
!   return
! end if

! call uA%quat_normalize()
! errcnt = 0
! dd = cabs(uA)
! do i=1,4
!   if (abs(dd(i)-1.D0).gt.epsd) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 49
!   write (*,"('double precision normalization test failed ')") 
!   return
! end if

! ! make u a rotation of 120° around the [111] axis 
! uA = QuaternionArray_T( n=4,  qd=reshape( (/ 0.5_dbl, 0.5_dbl, 0.5_dbl, 0.5_dbl, &
!                                              0.5_dbl, 0.5_dbl, 0.5_dbl, 0.5_dbl, &
!                                              0.5_dbl, 0.5_dbl, 0.5_dbl, 0.5_dbl, &
!                                              0.5_dbl, 0.5_dbl, 0.5_dbl, 0.5_dbl /), (/ 4, 4 /) ) )
! rvecarrd = uA%quat_Lp(vecd)
! errcnt = 0
! do i=1,4 
!   diffd = sum( abs(rvecarrd(:,i) - (/ 0.D0, 1.D0, 0.D0/)) )
!   if (diffd.gt.epsd) errcnt = errcnt+1
! end do
! if (errcnt.gt.0) then 
!   res = 50
!   write (*,"('double precision quaternion array vector rotation test failed ')")
!   return
! end if


! !===================================================
! !==========Single Precision Array Tests=============
! !===================================================

! !===================================================
! ! declare two arrays with identical quaternions to simplify the testing ... 
! uA = QuaternionArray_T( n=4, q = reshape( (/ 1.0_sgl,2.0_sgl,3.0_sgl,4.0_sgl, &
!                                              1.0_sgl,2.0_sgl,3.0_sgl,4.0_sgl, &
!                                              1.0_sgl,2.0_sgl,3.0_sgl,4.0_sgl, &
!                                              1.0_sgl,2.0_sgl,3.0_sgl,4.0_sgl /), (/ 4, 4 /) ) )

! vA = QuaternionArray_T( n=4, q = reshape( (/ 5.0_sgl,6.0_sgl,7.0_sgl,8.0_sgl, &
!                                              5.0_sgl,6.0_sgl,7.0_sgl,8.0_sgl, &
!                                              5.0_sgl,6.0_sgl,7.0_sgl,8.0_sgl, &
!                                              5.0_sgl,6.0_sgl,7.0_sgl,8.0_sgl /), (/ 4, 4 /) ) )

! wA = uA + vA 
! errcnt = 0
! do i=1,4
!   u = wA%getQuatfromArray(i)
!   if (.not.(u%octsequal(ressupv))) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 51
!   write (*,"('single precision quaternion array addition test failed = ')")
!   return
! end if

! wA = uA - vA 
! errcnt = 0
! do i=1,4
!   u = wA%getQuatfromArray(i)
!   if (.not.(u%octsequal(ressumv))) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 52
!   write (*,"('single precision quaternion array subtraction test failed = ')")
!   return
! end if

! wA = uA * vA 
! errcnt = 0
! do i=1,4
!   u = wA%getQuatfromArray(i)
!   if (.not.(u%octsequal(ressutv))) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 53
!   write (*,"('single precision quaternion array multiplication test failed = ')")
!   return
! end if

! wA = uA / vA 
! errcnt = 0
! do i=1,4
!   u = wA%getQuatfromArray(i)
!   if (.not.(u%octsequal(ressudv))) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 54
!   write (*,"('single precision quaternion array division test failed = ')")
!   return
! end if

! wA = uA * b 
! errcnt = 0
! do i=1,4
!   u = wA%getQuatfromArray(i)
!   if (.not.(u%octsequal(ressuta))) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 55
!   write (*,"('single precision quaternion array scalar multiplication test failed = ')")
!   return
! end if

! wA = conjg(uA) 
! errcnt = 0
! do i=1,4
!   u = wA%getQuatfromArray(i)
!   if (.not.(u%octsequal(ressuc))) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 56
!   write (*,"('single precision quaternion array conjugation test failed = ')")
!   return
! end if

! d = cabs(uA)
! errcnt = 0
! do i=1,4
!   if (abs(d(i)-absvals).gt.eps) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 57
!   write (*,"('single precision quaternion array abs test failed = ')")
!   return
! end if

! d = uA%quat_innerproduct(vA)
! errcnt = 0
! do i=1,4
!   if (abs(d(i)-ips).gt.eps) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 58
!   write (*,"('single precision quaternion array inner product test failed = ')")
!   return
! end if

! d = uA%quat_angle(vA)
! errcnt = 0
! do i=1,4
!   if (abs(d(i)-qans).gt.eps) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 59
!   write (*,"('single precision quaternion array interquaternion angle test failed ')")
!   return
! end if

! call uA%quat_normalize()
! errcnt = 0
! d = cabs(uA)
! do i=1,4
!   if (abs(d(i)-1.0).gt.eps) errcnt = errcnt + 1
! end do 
! if (errcnt.ne.0) then 
!   res = 60
!   write (*,"('single precision normalization test failed ')") 
!   return
! end if

! ! make u a rotation of 120° around the [111] axis 
! uA = QuaternionArray_T( n=4,  q=reshape( (/ 0.5_sgl, 0.5_sgl, 0.5_sgl, 0.5_sgl, &
!                                             0.5_sgl, 0.5_sgl, 0.5_sgl, 0.5_sgl, &
!                                             0.5_sgl, 0.5_sgl, 0.5_sgl, 0.5_sgl, &
!                                             0.5_sgl, 0.5_sgl, 0.5_sgl, 0.5_sgl /), (/ 4, 4 /) ) )
! rvecarrs = uA%quat_Lp(vecs)
! errcnt = 0
! do i=1,4 
!   diff = sum( abs(rvecarrs(:,i) - (/ 0.0, 1.0, 0.0/)) )
!   if (diff.gt.eps) errcnt = errcnt+1
! end do
! if (errcnt.gt.0) then 
!   res = 61
!   write (*,"('single precision quaternion array vector rotation test failed ')")
!   return
! end if


end subroutine MODOctonionsExecuteTest


end module MODOctonionsTest
