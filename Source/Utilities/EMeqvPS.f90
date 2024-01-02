! ###################################################################
! Copyright (c) 2015-2024, Marc De Graef Research Group/Carnegie Mellon University
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

program EMeqvPS
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/25/20
  !!
  !! For a given crystal symmetry and a pseudo-symmetry operation, list the independent equivalent operators

use mod_kinds
use mod_global
use mod_EMsoft
use mod_rotations 
use mod_quaternions
use mod_so3
use mod_io 
use mod_math

IMPLICIT NONE

character(fnlen)              :: progname = 'EMeqvPS.f90'
character(fnlen)              :: progdesc = 'List equivalent pseudo-symmetric rotations'
  
type(EMsoft_T)                :: EMsoft
type(IO_T)                    :: Message 
type(so3_T)                   :: SO
type(QuaternionArray_T)       :: Pm, dummy 
type(a_T)                     :: ax 
type(q_T)                     :: qu, qu1, qu2 
type(r_T)                     :: ro1, roFZ 
type(e_T)                     :: eu, eu2
type(Quaternion_T)            :: qm, qq, q2 

type(a_T),allocatable         :: axlist(:)
type(e_T),allocatable         :: eulist(:)
integer(kind=irg),allocatable :: unique(:)

integer(kind=irg)             :: pgnum, io_int(1), num, k, FZtype, FZorder, i, j
real(kind=dbl)                :: ro(4), rod(3), io_dbl(4), qus(4), io_dbl3(3), &
                                 eu1(3),  diff, a(4)


! program header and command line argument handling
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 910 /) )

! ask for point group number
call Message%ReadValue('Enter the point group number ',io_int)
pgnum = io_int(1)

! define the symmetry matrices 
SO = so3_T( pgnum )
call dummy%QSym_Init( pgnum, Pm )
num = Pm%getQnumber()

! read the pseudo-symmetric axis angle pair
call Message%ReadValue('Enter the pseudo-symmetry axis-angle pair (axis, angle°) ',io_dbl,4)
io_dbl(4) = io_dbl(4) * dtor
io_dbl(1:3) = io_dbl(1:3) / vecnorm(io_dbl(1:3))
ax = a_T( adinp = io_dbl )

! convert to quaternion
qu1 = ax%aq()

! for all symmetry quaternions Pm, compute the quaternion corresponding to the rotated rotation axis of ax;
! the rotation angle remains the same for all of them
allocate(axlist(num))
do j=1,num
  qm = Pm%getQuatfromArray(j)
  axlist(j) = a_T( adinp = (/ qm%quat_Lp( io_dbl(1:3) ) , io_dbl(4) /) )
end do

! take a random orientation, reduce it to the fundamental zone, and convert to quaternion q2
eu2 = e_T( edinp = (/ 94.43D0, 102.02D0, 228.77D0 /) * dtor )
call SO%ReduceOrientationtoRFZ(eu2, Pm, ro1)
qu2 = ro1%rq()
q2 = Quaternion_T( qd = qu2%q_copyd() )
call q2%quat_pos()

! apply this symmetrized pseudo-symmetry operators to the random orientation, and store them in an Euler list
allocate(eulist(num), unique(num))
unique = 1
do j=1,num
    qu = axlist(j)%aq() 
    qm = Quaternion_T( qd = qu%q_copyd() )
    call qm%quat_pos()
    qm = qm * q2
    qu = q_T( qdinp = qm%get_quatd() )
    eu = qu%qe()
    call SO%ReduceOrientationtoRFZ(eu, Pm, roFZ )
    eulist(j) = roFZ%re()
end do

! analyze this list for uniqueness, keeping only a single copy of each unique operator
do i=1,num-1 
  do j=i+1, num
    if (unique(i).eq.1) then 
      diff = sum(abs(eulist(i)%e_copyd()-eulist(j)%e_copyd()))
      if (diff.lt.0.01) unique(j) = 0
    end if 
  end do 
end do 

! and generate the list of pseudo-symmetric variants
call Message%printMessage(' ')
call Message%printMessage('Unique pseudo-symmetry operators (axis,angle°):')
do i=1,num
  if (unique(i).eq.1) then 
    a = axlist(i)%a_copyd()
    a(4) = a(4) / dtor
    io_dbl(1:4) = a(1:4)
    call Message%WriteValue('',io_dbl,4,"(4(F12.8,'  '))")
  end if 
end do 

end program EMeqvPS
