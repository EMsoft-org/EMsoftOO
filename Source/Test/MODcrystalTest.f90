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

module MODcrystalTest
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/13/20
  !!
  !! test of functions in the crystallography module 

use mod_kinds
use mod_global

contains 

subroutine MODcrystalExecuteTest(res) &
           bind(c, name='MODcrystalExecuteTest')    ! this routine is callable from a C/C++ program
!DEC$ ATTRIBUTES DLLEXPORT :: MODcrystalExecuteTest

use,INTRINSIC :: ISO_C_BINDING
use mod_symmetry
use mod_crystallography

IMPLICIT NONE

integer(C_INT32_T),INTENT(OUT)  :: res

type(Cell_T)        :: cell
type(SpaceGroup_T)  :: SG
integer(kind=irg)        :: i, j, k
real(kind=dbl)            :: diff
real(kind=sgl)            :: diffs
character(fnlen)    :: xtalname 

! reference values
real(kind=dbl),parameter:: eps = 1.0D-10
real(kind=sgl),parameter:: epss = 1.0E-7
real(kind=dbl)           :: dmt_ref(3,3), rmt_ref(3,3), dsm_ref(3,3), rsm_ref(3,3), vol_ref, alpha(3,3), &
                       d(3), d_ref(3), t(3), dp_ref, dp 
real(kind=sgl)            :: trsps_ref(27)
real(kind=dbl)            :: trspd_ref(27)
real(kind=sgl)            :: tins(3), touts(3)
real(kind=dbl)            :: tind(3), toutd(3)
character(1)        :: inspace(3), outspace(3)
real(kind=dbl)      :: ATOMpos(maxpasym,5)

!===================================================
! set the reference values (verified with Mathematica scripts)
dmt_ref = reshape( (/ 0.360000000000D+00,0.729322346201D-01,0.164169668796D+00, &
                      0.729322346201D-01,0.490000000000D+00,0.280000000000D+00, &
                      0.164169668796D+00,0.280000000000D+00,0.640000000000D+00/), (/ 3, 3 /))
rmt_ref = reshape( (/ 0.314579509017D+01,-.948450390289D-02,-.802794495252D+00, &
                      -.948450390289D-02,0.272111703095D+01,-.118805578250D+01, &
                      -.802794495252D+00,-.118805578250D+01,0.228820332109D+01/), (/ 3, 3 /))
dsm_ref = reshape( (/ 0.600000000000D+00,0.121553724367D+00,0.273616114661D+00, &
                      0.000000000000D+00,0.689365427109D+00,0.357924741383D+00, &
                      0.000000000000D+00,0.000000000000D+00,0.661077984283D+00/), (/ 3, 3 /))
rsm_ref = reshape( (/ 0.166666666667D+01,0.000000000000D+00,0.000000000000D+00, &
                      -.293878301181D+00,0.145060944555D+01,0.000000000000D+00, &
                      -.530709766715D+00,-.785397521912D+00,0.151268083914D+01/), (/ 3, 3 /))
vol_ref = 0.273434584192D+00

trsps_ref = (/  2.00000000,  3.00000000, -4.00000000, &
                0.28211802,  0.49586448, -1.39166069, &
                0.47019672,  0.63639730, -2.64431190, &
                9.47431469, 12.89660549,-14.32256985, &
                2.00000000,  3.00000000, -4.00000000, &
                3.33333325,  3.76407170, -9.46833515, &
                3.33333325,  3.76407170, -9.46833515, &
                1.20000005,  2.31120372, -1.02330554, &
                2.00000000,  3.00000000, -4.00000000 /)


trspd_ref =  (/ 2.000000000000D+00,3.000000000000D+00,-4.000000000000D+00, &
                0.282118028675D+00,0.495864469240D+00,-.139166066241D+01, &
                0.470196714458D+00,0.636397315794D+00,-.264431193713D+01, &
                0.947431464964D+01,0.128966052150D+02,-.143225696224D+02, &
                2.000000000000D+00,3.000000000000D+00,-4.000000000000D+00, &
                0.333333333333D+01,0.376407173429D+01,-.946833545571D+01, &
                0.333333333333D+01,0.376407173429D+01,-.946833545571D+01, &
                0.120000000000D+01,0.231120373006D+01,-.102330548366D+01, &
                2.000000000000D+00,3.000000000000D+00,-4.000000000000D+00 /)

d_ref = (/ 0.3660254037844386D0, 0.D0, 1.366025403784439D0 /)

dp_ref = 1.812474868352421D0

res = 0
!===================================================

!===================================================
! in this test, we generate a crystal structure and then we test all the distance and angle 
! routines in the crystal.f90 module.
! The test crystal structure is a hypothetical triclinic structure with 2 different atoms in the 
! unit cell and no symmetry:
!
! a = 0.6
! b = 0.7
! c = 0.8
! alpha = 60.0
! beta = 70.0
! gamma = 80.0
!
! Al @ (0.2, 0.3, 0.4)
! Ti @ (0.5, 0.6, 0.7)
!
xtalname = 'hypothetical.xtal'
call cell%setFileName(xtalname)
! incorrectly set the a parameter, then correct it 
call cell%setLatParm( (/ 0.3_dbl, 0.7_dbl, 0.8_dbl, 60.0_dbl, 70.0_dbl, 80.0_dbl /) )
call cell%setLatParm('a', 0.6_dbl)
call cell%setnumat(2_irg)
call cell%setAtomtype( (/ 13_irg, 22_irg /) )
call SG%setSpaceGroupNumber(1_irg)
call SG%setSpaceGroupXtalSystem(7_irg)
call cell%setXtalSystem(7_irg)
call SG%setSpaceGroupsetting(1_irg)
ATOMpos = 0.D0
ATOMpos(1,1:5) = (/ 0.2_dbl, 0.3_dbl, 0.4_dbl, 1.0_dbl, 0.004_dbl /)
ATOMpos(2,1:5) = (/ 0.5_dbl, 0.6_dbl, 0.7_dbl, 1.0_dbl, 0.006_dbl /)
call cell%setAtomPos(ATOMpos)

! generate the symmetry and all relevant matrices 
call cell%CalcMatrices()
SG = SpaceGroup_T( SGnumber = 1_irg )
! call SG%GenerateSymmetry(.TRUE.)
!===================================================

!===================================================
! compute the metric tensors, unit cell volume, and structure matrices

! direct metric tensor
diff = maxval(abs(cell%getdmt() - dmt_ref))
if (diff.gt.eps) then
  res = 1
  write (*,"('direct metric tensor difference = ',D18.10)") diff
  return
end if

! reciprocal metric tensor
diff = maxval(abs(cell%getrmt() - rmt_ref))
if (diff.gt.eps) then
  res = 2
  write (*,"('reciprocal metric tensor difference = ',D18.10)") diff
  return
end if

! direct structure matrix
diff = maxval(abs(cell%getdsm() - transpose(dsm_ref)))
if (diff.gt.eps) then
  res = 3
  write (*,"('direct structure matrix difference = ',D18.10)") diff
  return
end if

! reciprocal structure matrix
diff = maxval(abs(cell%getrsm() - transpose(rsm_ref)))
if (diff.gt.eps) then
  res = 4
  write (*,"('reciprocal structure matrix difference = ',D18.10)") diff
  return
end if

! cell volume
diff = abs(cell%getVolume() - vol_ref)
if (diff.gt.eps) then
  res = 5
  write (*,"('unit cell volume difference = ',D18.10)") diff
  return
end if
!===================================================

!===================================================
! TransSpace (s/d)
inspace = (/ 'd', 'r', 'c' /)
outspace = (/ 'd', 'r', 'c' /)
tins = (/ 2.0, 3.0, -4.0 /)
tind = (/ 2.D0, 3.D0, -4.D0 /)

! return code starts at 5
do i=1,3
  do j=1,3
    call cell%TransSpace(tins,touts,inspace(i),outspace(j))
    k = (i-1)*9 + (j-1)*3 + 1
    diffs = maxval(abs(touts(1:3) - trsps_ref(k:k+2)))
    if (diffs.gt.epss) then
       res = 5 + (i-1)*3+j
       write (*,*) 'unequal vectors:',touts,trsps_ref(k:k+2)
       return
    end if
  end do
end do

! return code starts at 15
do i=1,3
  do j=1,3
    call cell%TransSpace(tind,toutd,inspace(i),outspace(j))
    k = (i-1)*9 + (j-1)*3 + 1
    diff = maxval(dabs(toutd(1:3) - trspd_ref(k:k+2)))
    if (diff.gt.eps) then
       res = 15 + (i-1)*3+j
       write (*,*) 'unequal vectors:',toutd,trspd_ref(k:k+2)
       return
    end if
  end do
end do
!===================================================

! transformCoordinates, using a passive +30° rotation around the y-axis
alpha = 0.D0 
alpha(1,1) = cos( 30.D0 * dtor )
alpha(2,2) = 1.D0
alpha(3,3) = cos( 30.D0 * dtor )
alpha(3,1) = sin( 30.D0 * dtor )
alpha(1,3) =-sin( 30.D0 * dtor )
d = (/ 1.D0, 0.D0, 1.D0 /)

call cell%transCoor(d,t,alpha,'d','on')
diff = sum(abs(t - d_ref))
if (diff.gt.eps) then
  res = 25
  write (*,"('error in transCoor test = ',D18.10)") diff
  return
end if

call cell%transCoor(d,t,alpha,'r','no')
diff = sum(abs(t - d_ref(3:1:-1)))
if (diff.gt.eps) then
  res = 26
  write (*,"('error in transCoor test = ',D18.10)") diff
  return
end if


! CalcDot 
d = (/ 1.D0, 1.D0, 1.D0 /)
t = (/ -1.D0, -1.D0, 3.D0 /)

dp = cell%CalcDot(d,t,'d')
diff = abs( dp - dp_ref )
if (diff.gt.eps) then
  res = 27
  write (*,"('error in dot product test = ',D18.10)") diff
  return
end if

dp = cell%CalcDot(d,t,'c')
diff = abs( dp - 1.D0 )
if (diff.gt.eps) then
  res = 28
  write (*,"('error in dot product test = ',D18.10)") diff
  return
end if

! NormVec 
d = (/ 1.D0, 1.D0, 1.D0 /)
call cell%NormVec(d,'d')
diff = abs( d(1) - 0.6294160183482335D0 )
if (diff.gt.eps) then
  res = 29
  write (*,"('error in vector normalization test = ',D18.10)") diff
  return
end if

d = (/ 1.D0, 1.D0, 1.D0 /)
call cell%NormVec(d,'c')
diff = abs( d(1) - 1.D0/sqrt(3.D0) )
if (diff.gt.eps) then
  res = 30
  write (*,"('error in vector normalization test = ',D18.10)") diff
  return
end if

! CalcLength
d = (/ 1.D0, 1.D0, 1.D0 /)
dp = cell%CalcLength(d,'d')
diff = abs( dp - 1.588774309596194D0 )
if (diff.gt.eps) then
  res = 31
  write (*,"('error in vector length test = ',D18.10)") diff
  return
end if

! CalcAngle (s/d)
d = (/ 1.D0, 1.D0, 1.D0 /)
t = (/ -1.D0, -1.D0, 3.D0 /)
dp = cell%CalcAngle(d,t,'d')
diff = abs( dp - 0.9715366701557355D0 )
if (diff.gt.eps) then
  res = 32
  write (*,"('error in angle test = ',D18.10)") diff
  return
end if

dp = cell%CalcAngle(d,t,'c')
diff = abs( dp - 1.39582728112920776D0 )
if (diff.gt.eps) then
  res = 33
  write (*,"('error in angle test = ',D18.10)") diff
  return
end if

! more to come...



end subroutine MODcrystalExecuteTest


end module MODcrystalTest
