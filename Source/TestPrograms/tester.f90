program tester 

use mod_global 
use mod_kinds
! use mod_symmetry
! use mod_crystallography
! use mod_QCsymmetry
! use mod_QCcrystallography
use mod_io
! use mod_dualquaternions
use mod_octonions

IMPLICIT NONE 

! type(SpaceGroup_T)          :: SG 
! type(QCSpaceGroup_T)        :: QCSG 
! type(IO_T)                  :: Message
! integer(kind=irg)           :: i, j, isg, SamplingType, isym

! type(DualQuaternion_T)      :: dq1, dq2, dqp, tmp
! real(kind=dbl)              :: v(3)

type(Octonion_T)                :: a, b, c, d 
real(kind=dbl)                  :: onorm
real(kind=sgl)                  :: onorms

a = Octonion_T( od = (/ 1.D0, 2.D0, 3.D0, 4.D0, 5.D0, 6.D0, 7.D0, 8.D0 /) )
b = Octonion_T( od = (/ 1.D0, 3.D0, 5.D0, 7.D0, 9.D0, 2.D0, 4.D0, 6.D0 /) )
c = Octonion_T( od = (/ 8.D0, 7.D0, 6.D0, 5.D0, 4.D0, 3.D0, 2.D0, 1.D0 /) )

d = a+b 
call d%oct_print('a+b = ')

d = a-b 
call d%oct_print('a-b = ')

d = a*sqrt(2.D0)
call d%oct_print('a*s = ')

d = a*b
call d%oct_print('a*b = ')

d = a/b 
call d%oct_print('a/b = ')

onorm = cabs(a)
write (*,*) 'a%norm = ', onorm

d = conjg(a)
call d%oct_print('conj(a) = ')

d = a%octinverse()
call d%oct_print('a%inv = ')

a = Octonion_T( o = (/ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 /) )
b = Octonion_T( o = (/ 1.0, 3.0, 5.0, 7.0, 9.0, 2.0, 4.0, 6.0 /) )
c = Octonion_T( o = (/ 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 /) )

d = a+b 
call d%oct_print('a+b = ')

d = a-b 
call d%oct_print('a-b = ')

d = a*sqrt(2.D0)
call d%oct_print('a*s = ')

d = a*b 
call d%oct_print('a*b = ')

d = a/b 
call d%oct_print('a/b = ')

onorms = real(cabs(a))
write (*,*) 'a%norm = ', onorm

d = conjg(a)
call d%oct_print('conj(a) = ')

d = a%octinverse()
call d%oct_print('a%inv = ')



! ! simple test of the dual quaternion package 
! dq1 = DualQuaternion_T( qd = (/ 1.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /) )
! dq2 = DualQuaternion_T( qd = (/ 1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0 /) )

! dqp = dq1 + dq2 
! call dqp%dualquat_print()

! dqp = dq1 * dq2 
! call dqp%dualquat_print()

! dqp = dq1 / dq2 
! call dqp%dualquat_print()

! dqp = DualQuaternion_T( qd = (/ 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /) )
! call dqp%generate_dualquat( cPi/2.D0, (/ 0.D0, 0.D0, 1.D0/), (/1.D0, 0.D0, 0.D0/), 5.D0 )
! ! dq2 = DualQuaternion_T( qd = (/ sqrt(3.D0), 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /) )
! ! dqp = dq2 * dqp
! call dqp%dualquat_print()
! v = (/ 10.D0, 0.D0, 0.D0 /)
! v = dqp%dualquat_Lp( v )
! write (*,*) v



! ! test for the 230 space groups to ensure that all symmetry codes are correct
! ! both for point groups and k-vector sampling 

! open(unit=dataunit, file='sym-new.txt', status='unknown',form='formatted')

! do isg=1,230 
!   call SG%setSpaceGroupNumber(isg)
!   call SG%setSpaceGrouptrigonal(.FALSE.)
!   if ((isg.ge.143).and.(isg.le.167)) call SG%setSpaceGrouptrigonal(.TRUE.)

!   j=0
!   do i=1,32
!     if (SGPG(i).le.isg) j=i
!   end do
!   isym = j
!   SamplingType = PGSamplingType(isym)

! ! next, intercept the special cases (hexagonal vs. rhombohedral cases that require special treatment)
!   if ((SamplingType.eq.-1).or.(isym.eq.14).or.(isym.eq.26)) then
!     SamplingType = SG%getHexvsRho(isym)
!   end if
!   write (dataunit,"(3I6)") isg, isym, SamplingType
! end do 

! close(unit=dataunit, status='keep')


! stop 

! do i=6,11 
!   if (i.ne.7) then 
!     QCSG = QCspacegroup_T( nD = 3, QCtype = 'Ico')
!     call QCSG%setSGnum(i)
! ! call QCcell_icosahedral%setMetricParametersQC()
!     call QCSG%GenerateQCSymmetry(dopg=.FALSE.) 
!     write (*,*) ' order of icosahedral space group ', i, ' equals ',QCSG%getnsym()
!   end if 
! end do

end program