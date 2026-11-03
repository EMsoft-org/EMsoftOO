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

module MODso3Test

use, intrinsic :: iso_c_binding, only: c_int
use mod_kinds, only: dbl, irg
use mod_so3, only: so3_T
use mod_rotations, only: q_T, r_T, setRotationPrecision
use mod_quaternions, only: Quaternion_T, QuaternionArray_T
use mod_symmetry, only: PGTHD

implicit none
private
public :: MODso3ExecuteTest

contains

subroutine MODso3ExecuteTest(res) bind(c, name='MODso3ExecuteTest')
!DEC$ ATTRIBUTES DLLEXPORT :: MODso3ExecuteTest
integer(c_int), intent(out) :: res
real(kind=dbl) :: q(4), q1(4), q2(4), reduced1(4), reduced2(4), theta, distance
integer(kind=irg) :: i
integer(kind=irg), parameter :: cyclicPG(3) = (/21, 16, 9/)
integer(kind=irg), parameter :: cyclicCode(3) = (/201, 221, 231/)
type(QuaternionArray_T) :: qdummy, Pm
! PGTHD in mod_symmetry: 2 -> 3, 3 -> 16, 4 -> 9, 6 -> 21, 422 -> 12, 432 -> 30.
! q_T and Quaternion_T are scalar first: (w,x,y,z). mod_global sets epsijkd=+1,
! giving the standard Hamilton product. EMsoft interprets rotations passively;
! symmetry acts on the left, Pm(j)*q, as in EMeqvrot and ReduceOrientationtoRFZ.
! Rodrigues storage is (unit axis, tan(angle/2)), not a Cartesian 3-vector.
! Run every assertion; res retains the first distinct failing assertion code.
res = 0
call setRotationPrecision('Double')
q = (/5.3028761936245346D-17, 0.86602540378443871D0, &
      0.50000000000000011D0, 3.0616169978683836D-17/)
call checkReduction(21, q, 101, 'basal half-turn with roundoff', &
  'Length > tol=1e5 is clamped to infinity; insideCyclicFZ then requires x(3) == 0 exactly. '// &
  'Without the clamp, abs(x(3)*x(4)) ~ 3e-17*1.9e16 ~ 0.5 > tan(pi/12). '// &
  'Here q_T%qr already returns infinity for w < 1e-10, before the reducer clamp.', res)

q = (/3.5355346428289641D-10, 0.70710678118654757D0, &
      0.70710678118654757D0, 3.5355346428289641D-10/)
q = q / sqrt(sum(q*q))
do i=1,size(cyclicPG)
  call checkReduction(cyclicPG(i), q, cyclicCode(i), 'near half-turn', &
    'A finite symmetry image lies in the cyclic slab; the 1e5 clamp makes its length infinite '// &
    'and the exact zero-axis test rejects it (6-fold: r_z ~ tan(15 deg), length ~ 2e9).', res)
end do

! On develop c978e1f, cyclic 2 passes the same input and is omitted.
! The supplied 422 boundary input and all eight symmetry images of
! 45 degrees minus 1e-12 radians about [001] also pass; no 422 defect is claimed.

call qdummy%QSym_Init(30, Pm)
theta = acos(-1.D0)/4.D0
q1 = (/cos(theta/2.D0), 0.D0, 0.D0, sin(theta/2.D0)/)
q2 = (/cos(theta/2.D0), 0.D0, 0.D0, -sin(theta/2.D0)/)
write(*,'(/,A)') '432 boundary uniqueness: inputs are +45 and -45 degrees about [001].'
write(*,'(A,4ES25.16)') '  q1 (w,x,y,z) = ', q1
write(*,'(A,4ES25.16)') '  q2 (w,x,y,z) = ', q2
distance = symmetryDistance(q1, q2, Pm)
write(*,'(A,ES25.16)') '  input symmetry chord distance = ', distance
call check(distance < 1.D-12, 401, '432 (PG 30): inputs must be symmetry-equivalent', res)
call checkReduction(30, q1, 402, '+45-degree boundary', 'Returned rotation must lie in FZ and preserve symmetry.', &
                    res, reduced1)
call checkReduction(30, q2, 404, '-45-degree boundary', 'Returned rotation must lie in FZ and preserve symmetry.', &
                    res, reduced2)
distance = chordDistance(reduced1, reduced2)
write(*,'(A,ES25.16)') '  reduced representative chord distance = ', distance
call check(distance < 1.D-12, 406, &
  '432 (PG 30): non-unique FZ representative on the boundary (first symmetry operator inside the FZ wins)', res)
write(*,'(/,A,I0)') 'MODso3 first failing assertion code (0 means all passed): ', res
end subroutine MODso3ExecuteTest

subroutine checkReduction(pgnum, input, code, label, cause, res, reduced)
integer(kind=irg), intent(in) :: pgnum, code
real(kind=dbl), intent(in) :: input(4)
character(*), intent(in) :: label, cause
integer(c_int), intent(inout) :: res
real(kind=dbl), optional, intent(out) :: reduced(4)
type(so3_T) :: SO
type(QuaternionArray_T) :: qdummy, Pm
type(q_T) :: qu, returned
type(r_T) :: roFZ
real(kind=dbl) :: output(4), distance
logical :: inside

SO = so3_T(pgnum)
call qdummy%QSym_Init(pgnum, Pm)
qu = q_T(qdinp=input)
call SO%ReduceOrientationtoRFZ(qu, Pm, roFZ)
inside = SO%IsinsideFZ(roFZ)
returned = roFZ%rq()
output = returned%q_copyd()
distance = symmetryDistance(input, output, Pm)
write(*,'(/,A,A,A,I0,2A)') 'Point group ', trim(adjustl(PGTHD(pgnum))), ' (PG ', pgnum, '): ', label
write(*,'(A,4ES25.16)') '  input (w,x,y,z) = ', input
write(*,'(A,4ES25.16)') '  returned (axis,length) = ', roFZ%r_copyd()
write(*,'(A,4ES25.16)') '  returned (w,x,y,z) = ', output
write(*,'(A,L1,A,ES25.16)') '  IsinsideFZ = ', inside, '; min symmetry chord distance = ', distance
call check(inside, code, 'returned rotation must be inside FZ', res)
call check(distance < 1.D-9, code+1, 'returned rotation must be symmetry-equivalent to input (<1e-9)', res)
if (.not.inside.or..not.(distance < 1.D-9)) write(*,'(A)') '  Cause: '//cause
if (present(reduced)) reduced = output
end subroutine checkReduction

! The oracle enumerates the supplied group, without calling either FZ method.
! Euclidean chord distance avoids acos(dot) loss of precision near coincidence.
function symmetryDistance(input, output, Pm) result(distance)
real(kind=dbl), intent(in) :: input(4), output(4)
type(QuaternionArray_T), intent(inout) :: Pm
real(kind=dbl) :: distance
integer(kind=irg) :: j
type(Quaternion_T) :: image

distance = huge(1.D0)
do j=1,Pm%getQnumber()
  image = Pm%getQuatfromArray(j)*Quaternion_T(qd=input)
  distance = min(distance, chordDistance(image%get_quatd(), output))
end do
end function symmetryDistance

function chordDistance(a, b) result(distance)
real(kind=dbl), intent(in) :: a(4), b(4)
real(kind=dbl) :: distance

distance = min(sqrt(sum((a-b)**2)), sqrt(sum((a+b)**2)))
end function chordDistance

subroutine check(passed, code, message, res)
logical, intent(in) :: passed
integer(kind=irg), intent(in) :: code
character(*), intent(in) :: message
integer(c_int), intent(inout) :: res

if (passed) then
  write(*,'(A,I0,2A)') '  PASS [', code, '] ', message
else
  write(*,'(A,I0,2A)') '  FAIL [', code, '] ', message
  if (res == 0) res = code
end if
end subroutine check

end module MODso3Test
