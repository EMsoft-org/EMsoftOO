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

module mod_BFDF
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/21/24
  !!
  !! auxiliary routines for the EMTBBFDF program from the CTEM text book
  !! [this replaces the original BFDF.routines file]

use mod_kinds
use mod_global
use mod_io

IMPLICIT NONE 

contains

!--------------------------------------------------------------------------
recursive subroutine BFDF_straight_wedge(jdim, z, exe, xig)
!DEC$ ATTRIBUTES DLLEXPORT :: BFDF_straight_wedge
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! generate arrays for a straight wedge foil shape

integer(kind=irg),INTENT(IN)      :: jdim 
real(kind=sgl),INTENT(INOUT)      :: z(jdim,jdim)
real(kind=sgl),INTENT(INOUT)      :: exe(jdim,jdim)
real(kind=sgl),INTENT(IN)         :: xig

type(IO_T)                        :: Message 

real(kind=sgl)                    :: io_real(2), zmin, zmax, dz, ss 
integer(kind=irg)                 :: i, j

! create the thickness array 
call Message%ReadValue(' Enter minimum,maximum thickness [nm]  = ', io_real, 2)
zmin = io_real(1)
zmax = io_real(2)
dz = (zmax-zmin)/float(jdim)
do i=1,jdim
  do j=1,jdim
    z(i,j) = float(j)*dz
  end do
end do

! create the excitation error array 
call Message%ReadValue(' Enter w-parameter (s * xi)  = ', io_real, 1)
ss = io_real(1)
ss = ss/xig
do i=1,jdim
  do j=1,jdim
    exe(i,j) = ss
  end do
end do

end subroutine BFDF_straight_wedge

!--------------------------------------------------------------------------
recursive subroutine BFDF_bent_wedge(jdim, z, exe, xig)
!DEC$ ATTRIBUTES DLLEXPORT :: BFDF_bent_wedge
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! generate arrays for a bent wedge foil shape

integer(kind=irg),INTENT(IN)      :: jdim 
real(kind=sgl),INTENT(INOUT)      :: z(jdim,jdim)
real(kind=sgl),INTENT(INOUT)      :: exe(jdim,jdim)
real(kind=sgl),INTENT(IN)         :: xig

type(IO_T)                        :: Message 

real(kind=sgl)                    :: io_real(2), zmin, zmax, dz, ss, sm 
integer(kind=irg)                 :: i, j

! create the thickness array 
call Message%ReadValue(' Enter minimum,maximum thickness [nm]  = ', io_real, 2)
zmin = io_real(1)
zmax = io_real(2)
dz = (zmax-zmin)/float(jdim)
do i=1,jdim
  do j=1,jdim
    z(i,j) = float(j)*dz
  end do
end do

! create the excitation error array 
call Message%ReadValue(' Enter w-parameter (s * xi)  = ', io_real, 1)
sm = io_real(1)
do i=1,jdim
  ss = (-sm+2.0*sm*float(i-1)/float(jdim))/xig
  do j=1,jdim
    exe(i,j) = ss
  end do
end do

end subroutine BFDF_bent_wedge

!--------------------------------------------------------------------------
recursive subroutine BFDF_bent_hole(jdim, z, exe, xig)
!DEC$ ATTRIBUTES DLLEXPORT :: BFDF_bent_hole
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! generate arrays for a bent wedge foil shape with a hole

integer(kind=irg),INTENT(IN)      :: jdim 
real(kind=sgl),INTENT(INOUT)      :: z(jdim,jdim)
real(kind=sgl),INTENT(INOUT)      :: exe(jdim,jdim)
real(kind=sgl),INTENT(IN)         :: xig

type(IO_T)                        :: Message 

real(kind=sgl)                    :: io_real(2), zmin, zmax, dz, ss, sm, f1, dd 
integer(kind=irg)                 :: i, j, ihole, jhole, irad, jrad

! ------
call Message%ReadValue(' Enter minimum,maximum thickness [nm]  = ', io_real, 2)
zmin = io_real(1)
zmax = io_real(2)
dz = (zmax-zmin)/float(jdim)
ihole = 4*jdim/6
jhole = jdim/2
! outer radius
irad = 150
! inner radius
jrad = 40
f1 = 1.0/float(irad-jrad)
do i=1,jdim
  do j=1,jdim
    dd = sqrt(float(i-ihole)**2+float(j-jhole)**2)
! inside or outside the hole ?
    z(i,j) = float(j)*dz
    if (dd.lt.jrad) then
      z(i,j) = 0.D0
    else
      if (dd.lt.irad) then
        z(i,j)=z(i,j)*f1*(dd-float(jrad))
      end if
    end if
  end do
end do

! create the excitation error array 
call Message%ReadValue(' Enter w-parameter (s * xi)  = ', io_real, 1)
sm = io_real(1)
do i=1,jdim
  ss = (-sm+2.0*sm*float(i-1)/float(jdim))/xig
  do j=1,jdim
    exe(i,j) = ss
  end do
end do

end subroutine BFDF_bent_hole

!--------------------------------------------------------------------------
recursive subroutine BFDF_random_thickness(jdim, z, exe, xig)
!DEC$ ATTRIBUTES DLLEXPORT :: BFDF_random_thickness
!! author: MDG
!! version: 1.0
!! date: 02/21/24
!!
!! generate arrays for a bent random thickness foil shape with (an) optional hole(s)

integer(kind=irg),INTENT(IN)      :: jdim 
real(kind=sgl),INTENT(INOUT)      :: z(jdim,jdim)
real(kind=sgl),INTENT(INOUT)      :: exe(jdim,jdim)
real(kind=sgl),INTENT(IN)         :: xig

type(IO_T)                        :: Message 

real(kind=sgl)                    :: io_real(2), zav, zdev, ff, T0, T1, r(200), p(200), fi, fj, &
                                     zmin, zmax, ztot, q, ss 
integer(kind=irg)                 :: i, j, k, kk, nw, io_int(1)
integer(kind=irg)                 :: values(1:8)
integer, allocatable              :: seed(:)
real(kind=dbl)                    :: zz

T0 = SECNDS(0.0)
call Message%printMessage(' Average thickness, minimum thickness [nm]  = ')
call Message%ReadValue('   (minimum thickness may be <0, to create holes) ', io_real, 2)
zav = io_real(1)
zdev = io_real(2)
call Message%ReadValue(' How many random waves ? (<100) = ', io_int, 1)
nw = io_int(1)
nw = min(nw,100)
ff = cPi/float(jdim)

  ! superimpose a bunch of waves with random periodicities and phases
call date_and_time(values=values)
kk = 1
call random_seed(size=kk)
allocate(seed(1:kk))
seed(:) = values(8)
call random_seed(put=seed)
do k=1,2*nw
 call random_number(zz)
 r(k) = zz*10.0
 call random_number(zz)
 p(k) = zz*3.0
end do

do i=1,jdim
  fi = float(i)*ff
  do j=1,jdim
    fj=float(j)*ff
    z(i,j) = 0.0 
    do k=1,nw
      z(i,j) = z(i,j) + cos(r(2*k)*fi+p(2*k))*cos(r(2*k-1)*(fi+fj)+p(2*k-1)) + &
                        sin(r(2*k-1)*fj+p(2*k))*sin(r(2*k)*(fi-fj)+p(2*k-1))
    end do
  end do
end do
zmin = minval(z)
zmax = maxval(z)
ztot = sum(z)/float(jdim)**2

! create the excitation error array 
q = zdev/(ztot-zmin)
do i=1,jdim
  ss = (-0.5+1.0*float(i-1)/float(jdim))/xig
  do j=1,jdim
    exe(i,j) = ss
    z(i,j) = zav-zdev+q*(z(i,j)-zmin)
    if (z(i,j).lt.0.0) z(i,j)=0.0
  end do
end do
zmin = minval(z)
zmax = maxval(z)
ztot = sum(z)/float(jdim)**2
io_real(1) = ztot
call Message%WriteValue('Actual average thickness = ', io_real, 1)

end subroutine BFDF_random_thickness

end module mod_BFDF