! ###################################################################
! Copyright (c) 2016-2024, Marc De Graef Research Group/Carnegie Mellon University
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

module MODsymmetryTest
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/11/20
  !!
  !! perform a series of unit tests on the symmetry module 
  !!
  !! - construct a space group 
  !! - test the direct space matrices for the point group and the space group 
  !! - test the group name, order, crystal system, number of generators, number of matrices,
  !!   centrosymmetry or not...
  !! - then test a family of planes and directions plus an orbit
  !! - multiplicity, IsGAllowed, BF symmetry, pattern symmetry
  !!
  !! TODO: 

contains 

subroutine MODsymmetryExecuteTest(res) &
           bind(c, name='MODsymmetryExecuteTest')    ! this routine is callable from a C/C++ program
!DEC$ ATTRIBUTES DLLEXPORT :: MODsymmetryExecuteTest

use,INTRINSIC :: ISO_C_BINDING
use mod_kinds
use mod_symmetry

IMPLICIT NONE

integer(C_INT32_T),INTENT(OUT)  :: res

type(SpaceGroup_T)              :: SG, SG2, SG3 

real(kind=dbl), allocatable     :: SGdata(:,:,:), SGdirec(:,:,:), SGrecip(:,:,:)
integer(kind=irg)               :: sz(3), i, g(3), nn, sz2(2)
logical                         :: centro, symmorphic, allowed
real(kind=dbl)                  :: m4b4(4,4), m3b3(3,3), o = 1.D0, z = 0.D0, diff, sd(3), diffd, &
                                   cref(3), cref2(3), cref3(3), kk(3)
real(kind=dbl), parameter       :: epsd = 1.0D-12 
integer(kind=sgl),allocatable   :: itmp(:,:)
real(kind=dbl),allocatable      :: ctmp(:,:)
real(kind=sgl)                  :: s(3)

!===================================================
! set the reference values (verified with Mathematica scripts)

! correct answers for the SpaceGroup_T class
m4b4 = reshape( (/-o, z, z, z,  z, -o, z, z,  z, z, o, z,  z, z, z, o /), (/4,4/) )
m3b3 = reshape( (/-o, z, z,  z, -o, z,  z, z, o /), (/3,3/) )
cref = (/ 0.3D0, -0.2D0, -0.4D0 /)
cref2= (/-0.1D0, -0.1D0,  0.8D0 /)
cref3= (/ 0.2D0,  0.3D0,  0.9D0 /)

! initialize the error identifier to zero (should remain zero upon successful exit)
res = 0

!===================================================
!===================Start Tests=====================
!===================================================

! do a test with space group # 32 (Pba2) of order 4 
SG = SpaceGroup_T( SGnumber = 32, setting = 1 )

! test for a couple of parameters 
if (SG%getSpaceGroupName().ne.' P b a 2   ') then 
  res = 1
  write (*,"('space group name test failed ')")
  return
end if

if (SG%getSpaceGroupOrder().ne.4) then 
  res = 2
  write (*,"('space group order test failed ')")
  return
end if

if (SG%getSpaceGroupNumber().ne.32) then 
  res = 3
  write (*,"('space group number test failed ')")
  return
end if

if (SG%getSpaceGroupSetting().ne.1) then 
  res = 4
  write (*,"('space group setting test failed ')")
  return
end if

! this space group is not centrosymmetric ... 
centro = SG%getSpaceGroupCentro()
if (centro.eqv..TRUE.) then 
  res = 5
  write (*,"('space group centrosymmetry test failed ')")
  return
end if

if (SG%getSpaceGroupXtalSystem().ne.3) then 
  res = 6
  write (*,"('space group crystal system test failed ')")
  return
end if

if (SG%getSpaceGroupSymmorphic().eqv..TRUE.) then 
  res = 7
  write (*,"('space group non-symmorphic test failed ')")
  return
end if

! next, check some of the symmetry matrices in the data, direc, and recip allocatable arrays
! as well as their dimensions 

SGdata = SG%getSpaceGroupDataMatrices()

! check the shape
sz = shape(SGdata)
if ((sz(1).ne.4).or.(sz(2).ne.4).or.(sz(3).ne.4)) then 
  res = 8
  write (*,"('space group data array size test failed ')")
  return
end if

! and the values of the second matrix in the array
diff = sum(abs( SGdata(2,:,:) - m4b4(:,:) ) )
if (diff.gt.epsd) then 
  res = 9
  write (*,"('space group data array values test failed ')")
  return
end if

SGdirec = SG%getSpaceGroupPGdirecMatrices()

! check the shape
sz = shape(SGdirec)
if ((sz(1).ne.4).or.(sz(2).ne.3).or.(sz(3).ne.3)) then 
  res = 10
  write (*,"('space group direc array size test failed ')")
  return
end if

! and the values of the second matrix in the array
diff = sum(abs( SGdirec(2,:,:) - m3b3(:,:) ) )
if (diff.gt.epsd) then 
  res = 11
  write (*,"('space group direc array values test failed ')")
  return
end if

! next, generate two other space groups and do some similar tests 

! do a test with space group # 114 (P-42_1c) of order 8 
SG2 = SpaceGroup_T( SGnumber = 114, setting = 1 )

if (SG2%getSpaceGroupSymmorphic().eqv..TRUE.) then 
  res = 12
  write (*,"('space group non-symmorphic test failed ')")
  return
end if

if (SG2%getSpaceGroupOrder().ne.8) then 
  res = 13
  write (*,"('space group order test failed ')")
  return
end if

deallocate(SGdata)
SGdata = SG2%getSpaceGroupDataMatrices()

! check the shape
sz = shape(SGdata)
if ((sz(1).ne.8).or.(sz(2).ne.4).or.(sz(3).ne.4)) then 
  res = 14 
  write (*,"('space group data array size test failed ')")
  return
end if

! do a test with space group # 206 (Ia-3) of order 48 
SG3 = SpaceGroup_T( SGnumber = 206, setting = 1 )

if (SG3%getSpaceGroupOrder().ne.48) then 
  res = 15
  write (*,"('space group order test failed ')")
  return
end if

! this space group is not centrosymmetric ... 
centro = SG3%getSpaceGroupCentro()
if (centro.eqv..FALSE.) then 
  res = 16
  write (*,"('space group centrosymmetry test failed ')")
  return
end if

deallocate(SGdata)
SGdata = SG3%getSpaceGroupDataMatrices()

! check the shape
sz = shape(SGdata)
if ((sz(1).ne.48).or.(sz(2).ne.4).or.(sz(3).ne.4)) then 
  res = 17 
  write (*,"('space group data array size test failed ')")
  return
end if

! that concludes the symmetry checks; next we check a few other routines ...
g = (/ 1,0,0 /)
allowed = SG%IsGAllowed( g )
if (allowed.eqv..FALSE.) then 
  res = 18
  write (*,"('space group IsGAllowed test failed ')")
  return
end if

allowed = SG2%IsGAllowed( g )
if (allowed.eqv..FALSE.) then 
  res = 19
  write (*,"('space group IsGAllowed test failed ')")
  return
end if

g = (/ 1,0,0 /)
allowed = SG3%IsGAllowed( g )
if (allowed.eqv..TRUE.) then 
  res = 20
  write (*,"('space group IsGAllowed test failed ')")
  return
end if

! reset space group 3 and test 
call SG3%resetSpaceGroup()

if (SG3%getSpaceGroupNumber().ne.0) then 
  res = 21
  write (*,"('space group number test failed (after space group reset) ')")
  return
end if

! let's check a family
g = (/ 1, 2, 3 /)
call SG%CalcFamily(g, nn, 'd', itmp)

sz2 = shape(itmp)
if ((sz2(1).ne.4).or.(sz2(2).ne.3)) then 
  res = 22 
  write (*,"('space group CalcFamily size test failed ')")
  return
end if

if ((itmp(2,1).ne.-1).or.(itmp(2,2).ne.-2)) then 
  res = 23 
  write (*,"('space group CalcFamily test failed ')")
  return
end if

! check an orbit for the second space group
sd = (/ 0.2D0, 0.3D0, 0.4D0 /)
call SG2%CalcOrbit(sd, nn, ctmp)
sz2 = shape(ctmp)

if ((sz2(1).ne.8).or.(sz2(2).ne.3)) then 
  res = 24 
  write (*,"('space group CalcOrbit size test failed ')")
  return
end if

diffd = sum( abs( ctmp(6,:) - cref3(:)))
if (diffd.gt.epsd) then 
  res = 25 
  write (*,"('space group CalcOrbit test failed ')")
  return
end if

! check a star in SG ... 
deallocate(ctmp)
kk = (/  0.1D0, 0.1D0, 0.8D0 /)
call SG%CalcStar(kk, nn, ctmp, 'd')
sz2 = shape(ctmp)

if ((sz2(1).ne.4).or.(sz2(2).ne.3)) then 
  res = 26 
  write (*,"('space group CalcOrbit size test failed ')")
  return
end if

diffd = sum( abs( ctmp(2,:) - cref2(:)))
if (diffd.gt.epsd) then 
  res = 27 
  write (*,"('space group CalcOrbit test failed ')")
  return
end if





end subroutine MODsymmetryExecuteTest

end module MODsymmetryTest
