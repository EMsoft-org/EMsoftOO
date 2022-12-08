! ###################################################################
! Copyright (c) 2015-2023, Marc De Graef Research Group/Carnegie Mellon University
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

program EMeqvrot
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/25/20
  !!
  !! For a given crystal symmetry (point group number), compute equivalent orientations 

use mod_kinds
use mod_global
use mod_EMsoft
use mod_symmetry 
use mod_so3 
use mod_rotations 
use mod_quaternions 
use mod_io
use mod_math

IMPLICIT NONE

character(fnlen)            :: progname = 'EMeqvrot.f90'
character(fnlen)            :: progdesc = 'List equivalent rotations'

type(EMsoft_T)              :: EMsoft 
type(IO_T)                  :: Message 
type(so3_T)                 :: SO
type(SpaceGroup_T)          :: SG 
type(QuaternionArray_T)     :: qdummy, Pm
type(r_T)                   :: ro 
type(q_T)                   :: qu, qq 
type(a_T)                   :: ax
type(e_T)                   :: eu
type(Quaternion_T)          :: qm, qus 
type(Orientation_T)         :: ot, saveot

integer(kind=irg)           :: pgnum, Rtype, Otype, io_int(1), num, k
real(kind=dbl)              :: rod(3), io_dbl(4), io_dbl3(3), a(4), qqd(4)
logical                     :: next

! print some information
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 911 /) )

! ask for point group number
call SG%ListPointGroups()
call Message%ReadValue('Enter the desired point group number: ',io_int)
pgnum = io_int(1)

! what format are we using for the input?
call Message%ReadValue('Input representation: Rodrigues (1), axis-angle (2), quaternion (3), or Euler (4) ',io_int)
Rtype = io_int(1)
call Message%ReadValue('Output representation: Rodrigues (1), axis-angle (2), quaternion (3), or Euler (4) ',io_int)
Otype = io_int(1)

! get the fundamental zone parameters
SO = so3_T(pgnum)
call qdummy%QSym_Init( pgnum, Pm )
num = Pm%getQnumber()

! then in a loop ask for an orientation in Rodrigues form or axis angle
! and list all the equivalent orientations as axis angle pairs
next = .TRUE.
do while (next) 
  select case (Rtype)
    case (1) 
      call Message%ReadValue('Rodrigues vector (nx, ny, nz, |R|) :', io_dbl,4)
      ro = r_T( rdinp = io_dbl )
      qu = ro%rq()
    case (2) 
      call Message%ReadValue('Axis-angle pair (nx, ny, nz, omega°) :', io_dbl,4)
      io_dbl(4) = io_dbl(4) * dtor 
      io_dbl(1:3) = io_dbl(1:3) / vecnorm( io_dbl(1:3) )
      ax = a_T( adinp = io_dbl )
      qu = ax%aq()
    case (3) 
      call Message%ReadValue('Quaternion (q0, q1, q2, q3) :', io_dbl,4)
      qu = q_T( qdinp = io_dbl )
    case (4) 
      call Message%ReadValue('Euler (phi1°, Phi°, phi2°) :', io_dbl3,3)
      eu = e_T( edinp = io_dbl3 * dtor )
      qu = eu%eq()
    case default
  end select
  call Message%printMessage(' ')

  do k=1,num
    qm = Pm%getQuatfromArray(k)
    qus = qm * Quaternion_T( qd = qu%q_copyd() )

    qq = q_T( qdinp = qus%get_quatd() )
    
    qqd = qq%q_copyd()

    if (qqd(1).lt.0.0) then
      qqd = -qqd
      qus = Quaternion_T(qd = qqd)
      qq = q_T( qdinp = qus%get_quatd() )
    end if

    ot = Orientation_T( qq )
    ro = ot%getClass_r()
    ax = ot%getClass_a()
    eu = ot%getClass_e()

    if (SO%IsinsideFZ(ro).eqv..TRUE.) then
      saveot = ot
      select case(Otype)
        case(1)
          call Message%WriteValue('FZ  ', ro%r_copyd(), num=4)
        case(2)
          a = ax%a_copyd()
          a(4) = a(4) / dtor
          call Message%WriteValue('FZ  ', a, num=4)
        case(3)
          call Message%WriteValue('FZ  ', qq%q_copyd(), num=4)
        case(4)
          call Message%WriteValue('FZ  ', eu%e_copyd()/dtor, num=3)
        case default
      end select
    else
      select case(Otype)
        case(1)
          call Message%WriteValue('    ', ro%r_copyd(), num=4)
        case(2)
          a = ax%a_copyd()
          a(4) = a(4) / dtor
          call Message%WriteValue('    ', a, num=4)
        case(3)
          call Message%WriteValue('    ', qq%q_copyd(), num=4)
        case(4)
          call Message%WriteValue('    ', eu%e_copyd()/dtor, num=3)
        case default
      end select
    end if
  end do

! write all the equivalent representations for the one in the FZ
  call Message%printMessage(' ')
  call Message%printMessage('Equivalent representations for the rotation inside the FZ:')
  call saveot%print_orientation('d')
  call Message%printMessage(' ')
  call Message%ReadValue('---> Another one ? (1/0) ',io_int)
  if (io_int(1).eq.0) next=.FALSE.
end do

end program EMeqvrot
