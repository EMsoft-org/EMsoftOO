 program qtest

use mod_kinds
use mod_quaternions
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit

IMPLICIT NONE

! for single precision, use the following lines
!real(kind=sgl)  ::  u(4), v(4), w(4)
!real(kind=sgl) :: x, a=2
! define two quaternions (single)
!u = (/1.0,2.0,3.0,4.0/)
!v = (/5.0,6.0,7.0,8.0/)

! for double precision, uncomment the next set and comment the previous set lines
type(Quaternion_T)  ::  u, v, w
real(kind=dbl)  :: x, a=2.D0

write (stdout,*) 'no initialization'
u = Quaternion_T()
call u%quat_print 


! double
u = Quaternion_T( qd=(/1.D0,2.D0,3.D0,4.D0/) )
v = Quaternion_T( qd=(/5.D0,6.D0,7.D0,8.D0/) )


write (stdout,*) ' quaternion u '
call u%quat_print 
write (stdout,*) ' quaternion v '
call v%quat_print

! next, do all the operations to make sure that they are correct

write (stdout,*) '   addition u+v '
w = u+v
call w%quat_print
write (stdout,*) ' '

write (stdout,*) '   subtraction u-v '
w = u-v
call w%quat_print
write (stdout,*) ' '

write (stdout,*) '   scalar multiplication (both orderings)  '
!w = a*u
!call quat_print(w)
w = u*a
call w%quat_print
write (stdout,*) ' '

write (stdout,*) '   multiplication uv '
w = u*v
call w%quat_print
write (stdout,*) ' '

write (stdout,*) '   conjugate u '
w = conjg(u)
call w%quat_print
write (stdout,*) ' '

write (stdout,*) '   norm(u) '
x = cabs(u)
write (stdout,*) x
write (stdout,*) ' '

write (stdout,*) '   division u/v '
w = u/v
call w%quat_print
write (stdout,*) ' '

end program qtest