! ###################################################################
! Copyright (c) 2016-2026, Marc De Graef Research Group/Carnegie Mellon University
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

module MODproductgridTest
  !! author: ZTV
  !! version: 1.0
  !! date: 11/02/26
  !!
  !! test of the product grids of mod_productgrid and so3_T%sample_PRD against the point counts and the
  !! closest pairs of the reference implementation, and against the symmetry operators of QSym_Init

use mod_kinds
use mod_global

contains

subroutine MODproductgridExecuteTest(res) &
           bind(c, name='MODproductgridExecuteTest')    ! this routine is callable from a C/C++ program
!DEC$ ATTRIBUTES DLLEXPORT :: MODproductgridExecuteTest

use,INTRINSIC :: ISO_C_BINDING
use mod_math
use mod_quaternions
use mod_rotations
use mod_so3
use mod_productgrid

IMPLICIT NONE

integer(C_INT32_T),INTENT(OUT)  :: res

integer(kind=irg), parameter    :: ntests = 10
! point group number, frequency, number of orientations
integer(kind=irg), parameter    :: tests(3,ntests) = reshape( (/ 32,6,73, 33,4,54, 27,6,364, 29,3,115, 15,6,219, &
                                                                 20,3,230, 39,3,230, 40,3,345, 8,3,345, 5,3,690 /), &
                                                              (/ 3,ntests /))
! largest |<p, g q>| over all distinct pairs of the extended set, from the reference implementation
real(kind=dbl), parameter       :: closest(ntests) = (/ 0.984427733721423D0, 0.987898528638932D0, 0.991444861373811D0, &
                                                        0.978147600733806D0, 0.984427733721423D0, 0.978147600733806D0, &
                                                        0.978147600733806D0, 0.978147600733806D0, 0.978147600733806D0, &
                                                        0.978147600733806D0 /)
! the operators in the table of SYM_Qsymop are good to about 8 digits for 532
real(kind=dbl), parameter       :: tol = 1.D-7
! every point group with a product grid
integer(kind=irg), parameter    :: pgs(38) = (/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28, &
                                                29,30,31,32,33,37,38,39,40,41 /)

type(so3_T)                     :: SO
type(QuaternionArray_T)         :: Pm, dummy
type(Quaternion_T)              :: mu, cq
type(FZpointd), pointer         :: p
type(q_T)                       :: qin
type(r_T)                       :: roFZ
real(kind=dbl), allocatable     :: q(:,:)
real(kind=dbl)                  :: x(4), y(4), best
integer(kind=irg)               :: t, i, j, g, n, nop, ft, fo

res = 0

! for every point group the grid is in the zone of IsinsideFZ, and ReduceOrientationtoRFZ leaves it unchanged
do t=1,size(pgs)
  SO = so3_T( pgs(t), zerolist='SF' )
  call SO%sample_PRD( 6 )
  dummy = QuaternionArray_T( n=1, s='d' )
  call dummy%QSym_Init( pgs(t), Pm )
  p => SO%getListHead('SF')
  do i=1,SO%getListCount('SF')
    if (.not.SO%IsinsideFZ(p%rod)) then
      res = 400 + pgs(t)
      write (*,*) 'orientation outside IsinsideFZ for point group ', pgs(t), ' : ', i
      return
    end if
    qin = p%qu
    call SO%ReduceOrientationtoRFZ(qin, Pm, roFZ)
    x = roFZ%r_copyd()
    y = p%rod%r_copyd()
    if (maxval(abs(x-y)).gt.tol) then
      if ((y(4).ne.inftyd()).or.(maxval(abs(x(1:3)-y(1:3))).gt.tol)) then
        res = 500 + pgs(t)
        write (*,*) 'orientation changed by ReduceOrientationtoRFZ for point group ', pgs(t), ' : ', i
        return
      end if
    end if
    p => p%next
  end do
end do

do t=1,ntests
  SO = so3_T( tests(1,t), zerolist='SF' )
  call SO%sample_PRD( tests(2,t) )
  n = SO%getListCount('SF')

! the number of orientations
  if (n.ne.tests(3,t)) then
    res = t
    write (*,*) 'wrong number of orientations for test ', t, ' : ', n
    return
  end if

  call SO%getFZtypeandorder(ft, fo)
  if (productGridCount(ft, fo, merge(-1, 1, fo.lt.0), tests(2,t)).ne.int(n,ill)) then
    res = 100 + t
    write (*,*) 'wrong productGridCount for test ', t
    return
  end if

  allocate(q(4,n))
  p => SO%getListHead('SF')
  do i=1,n
    q(1:4,i) = p%qu%q_copyd()
    p => p%next
  end do

  dummy = QuaternionArray_T( n=1, s='d' )
  call dummy%QSym_Init( tests(1,t), Pm )
  nop = Pm%getQnumber()

  best = 0.D0
  do i=1,n
    mu = Quaternion_T( qd = q(1:4,i) )
    do g=1,nop
      cq = Pm%getQuatfromArray(g) * mu
      x = cq%get_quatd()

! the orientation is the symmetry equivalent closest to the identity, so it is in the zone
      if (abs(x(1)).gt.abs(q(1,i))+tol) then
        res = 200 + t
        write (*,*) 'orientation outside the zone for test ', t, ' : ', i
        return
      end if

! the closest pair of the extended set, leaving out each orientation paired with itself
      do j=1,n
        if ((j.eq.i).and.(g.eq.1)) cycle
        best = max(best, abs(dot_product(q(1:4,j), x)))
      end do
    end do
  end do

  if (abs(best-closest(t)).gt.tol) then
    res = 300 + t
    write (*,*) 'wrong closest pair for test ', t, ' : ', best, closest(t)
    return
  end if
  deallocate(q)
end do

end subroutine MODproductgridExecuteTest

end module MODproductgridTest
