! ###################################################################
! Copyright (c) 2013-2020, Marc De Graef Research Group/Carnegie Mellon University
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

program t 
use mod_kinds
use mod_global
use mod_math
use mod_rotations
use mod_Lambert
use,INTRINSIC :: ISO_C_BINDING

IMPLICIT NONE 

type(r_T)     :: r, rA, rB, rC  
type(s_T)     :: s, ss 
type(q_T)     :: q 
type(a_T)     :: a 
type(e_t)     :: e
type(orientation_T)             :: oo

integer(kind=irg) :: ierr 
real(kind=dbl)    :: x =  1.D0/dsqrt(3.D0)

type(e_T)                       :: ine, oute, int1e, int2e
type(o_T)                       :: ino, outo, int1o, int2o
type(a_T)                       :: ina, outa, int1a, int2a
type(r_T)                       :: inr, outr, int1r, int2r
type(q_T)                       :: inq, outq, int1q, int2q
type(h_T)                       :: inh, outh, int1h, int2h
type(c_T)                       :: inc, outc, int1c, int2c
type(s_T)                       :: ins, outs, int1s, int2s
type(v_T)                       :: inv, outv, int1v, int2v
type(orientation_T)             :: ot

integer(C_INT32_T)              :: res

real(kind=dbl)                  :: iom(3,3), oom(3,3), omm(3,3), qd(4), qd2(4), ad(4), ad2(4), rd(4), rd2(4), &
                                   hd(3), ccd(3), sd(3), vd(3), diff, diffmax, aux, ivec(3)
real(kind=dbl),parameter        :: maxerr = 1.0D-9
integer(kind=irg)               :: tcnt, i,  numarg, testcounter, testsfailed
integer(kind=irg),parameter     :: rcnt = 75
character(fnlen)                :: arg
logical                         :: verbose


call setRotationPrecision('Double')

ine = e_T( edinp = cvtoRadians( (/ 180.D0, 90.D0, 0.D0 /) ) )
ot = orientation_T( ine )
call ot%print_orientation('d')

inv = ot%getClass_v()
ina = inv%va()
call ina%a_print('ina ')
outv = ina%av()
call outv%v_print('outv ')

inv = v_T( vdinp = (/ 0.0000000000D0,   -(dsqrt(2.D0)*0.5D0)*cPi,   -(dsqrt(2.D0)*0.5D0)*cPi /) )
ot = orientation_T( inv )
call ot%print_orientation('d')

stop
  inc = ot%getClass_c()
  call inc%c_print('inc ')
  int1e = inc%ce()
  ! call int1e%e_print('int1e ')
  ! int1a = int1e%ea()
  ! call int1a%a_print('int1a ')
  ! int1h = int1a%ah()
  ! call int1h%h_print('int1h ')
  ! outc = int1h%hc()

  outc = int1e%ec()
  call outc%c_print('outc ')
  ccd = outc%c_copyd()
  write (*,*) 'ccd = ', ccd 
  aux = maxval(abs(ccd))
  if (abs(aux-LPs%ap*0.5).lt.1.D-08) then
    diff = maxval(abs(ccd)-abs(ot%get_cd()))
  else
    diff = maxval(abs(ccd-ot%get_cd()))
  end if 
  write (*,*) 'c%ce() max cu difference = ', diff
 
  ! outc = int1h%hc()
  ! call inc%c_print('inc ')
  ! call int1e%e_print('int1e ')
  ! call int1h%h_print('int1h ')
  ! call outc%c_print('outc ')
  ! outc = int1e%ec()
  ! call outc%c_print('outc ')







stop
a = a_T( adinp = (/ 0.D0, 0.D0, -1.D0, cPi/4.D0 /) )

oo = orientation_T( a )

call oo%print_orientation('d')

call oo%print_orientation('r')


! let's do a 45° rotation around [111], and start with an axis angle pair, also print it

a = a_T( adinp = (/ x, x, x, cvtoRadians(45.D0) /) )
call a%a_print('Input axis angle pair : ')

! convert this to a quaternion and print it
q = a%aq()
call q%q_print('Quaternion            : ')

! convert this to a Rodrigues vector 
r = q%qr() 
call r%r_print('Rodrigues vector      : ')

! and back to an axis angle pair 
a = r%ra() 
call a%a_print('Resulting axis-angle  : ')

write (*,*) ' ----- '
write (*,*) 'setting precision to s'
call setRotationPrecision('s')
write (*,*) '    precision = ', getRotationPrecision()

q = s%sq()
s = r%rs()
r = q%qr()

write (*,*) 'setting precision to d'
call setRotationPrecision('Double')
write (*,*) '    precision = ', getRotationPrecision()

q = s%sq()
s = r%rs()
r = q%qr()

s = s_T( sdinp = (/ .1D0, .2D0, .3D0 /) )
call s%s_print('stereographic vector')
call setRotationPrecision('s')

ss = s_T( sinp = (/ .10, .20, .30 /) )
call ss%s_print('stereographic vector')

call setRotationPrecision('d')
rA = r_T( rdinp = (/ 0.D0, 0.D0, 1.D0, dtan(0.5D0*cvtoRadians(45.D0)) /) )
rB = r_T( rdinp = (/ 1.D0, 0.D0, 0.D0, dtan(0.5D0*cvtoRadians(60.D0)) /) )

rC = RodriguesProduct(rA, rB)
call rA%r_print('rA               : ')
call rB%r_print('rB               : ')
call rC%r_print('resulting vector : ')


end program t








