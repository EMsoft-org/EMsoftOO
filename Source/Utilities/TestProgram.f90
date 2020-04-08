
module modA 

use mod_kinds
use mod_global 

IMPLICIT NONE 

type, public :: type1 
    integer(kind=irg)   :: a 

    contains
        private 
        procedure, pass(self) :: writeHDFnamelist_
        procedure, pass(self) :: p_

        generic, public :: writeHDFnamelist => writeHDFnamelist_
        generic, public :: p => p_

end type type1

contains 

subroutine writeHDFnamelist_(self, n)

class(type1), INTENT(INOUT) :: self 
integer(kind=irg),INTENT(IN):: n 

write (*,*) '================'
write (*,*) 'type1:writeHDFnamelist: n = ', n 
write (*,*) 'type1:writeHDFnamelist: a = ', self%a 
write (*,*) '================'

end subroutine writeHDFnamelist_

subroutine p_(self, n)

class(type1), INTENT(INOUT) :: self 
integer(kind=irg),INTENT(IN):: n 

write (*,"(' ====>', I1,'<==== ')") n 

end subroutine p_

end module modA

! ===================================================
module modB 

use mod_kinds
use mod_global 
use modA

IMPLICIT NONE 

type, public, extends(type1) :: type2 
    integer(kind=irg)   :: b 

    contains
        procedure, pass(self) :: writeHDFnamelist_ => writelist

end type type2

contains 

subroutine writelist(self, n)

class(type2), INTENT(INOUT) :: self 
integer(kind=irg),INTENT(IN):: n 

write (*,*) '================'
write (*,*) 'type2:writeHDFnamelist: n = ', n 
write (*,*) 'type2:writeHDFnamelist: a = ', self%a 
write (*,*) 'type2:writeHDFnamelist: b = ', self%b 
write (*,*) '================'

end subroutine writelist

end module modB



program TestProgram

use mod_kinds
use mod_global
use modA 
use modB 

! use mod_EMsoft
! use mod_io
! use mod_timing
! use mod_quaternions 
! use mod_rng
use mod_HDFsupport
use ISO_C_BINDING
IMPLICIT NONE

type(type1)  :: t1 
type(type2)  :: t2 

character(fnlen)        :: progname = 'this is the program name'
character(kind=c_char)  :: Cprogname(fnlen)
integer(kind=irg)       :: i, slen
! character(fnlen)        :: progdesc = 'and this is the descriptor'
! character(fnlen)        :: m

! type(EMsoft_T)     :: EMsoft 
! type(IO_T)         :: Message
! type(Timing_T)     :: Timing

! integer(kind=irg)       :: io_int(2), status
! type(rng_t)             :: seed 
! type(QuaternionArray_T) :: qra , qrb


write (*,*) '->'//trim(progname)//'<-', len_trim(progname)

Cprogname = carstringify(progname)

do i=1,len_trim(progname) 
    write (*,"(A1$)") Cprogname(i:i)
end do 
write (*,*) ''

progname = ''
i=1
do while(Cprogname(i).ne.C_NULL_CHAR) 
  progname(i:i) = Cprogname(i)
  i = i+1
end do

write (*,*) '->'//trim(progname)//'<-', len_trim(progname)





stop

t1%a = 10
call t1%writeHDFnamelist(40)
call t1%p(1)

t2%a = 20 
t2%b = 30
call t2%writeHDFnamelist(60)
call t2%p(2)




! Timing = Timing_T( showDateTime = .TRUE. )

! call Timing % Time_tick()

! EMsoft = EMsoft_T(progname, progdesc, showconfig=.TRUE., makeconfig=.FALSE.) 

! call EMsoft%setConfigParameter('WyckoffPositionsfilename','no idea where this file is located ')

! write (*,*) 'modified string = ',trim(EMsoft%getConfigParameter('WyckoffPositionsfilename'))




! Message = IO_T()

! m = 'Revision = '//trim( EMsoft % getConfigParameter('EMsoftRevision') )
! call Message % printMessage(m)

! m = 'Completed path = '//trim( EMsoft % generateFilePath('EMsofttestpath','test.h5'))
! call Message % printMessage(m)

! m = 'Completed path = '//trim( EMsoft % generateFilePath('Templatecodefilename','')) 
! call Message % printMessage(m)


! io_int = (/ 10, 20 /)
! call Message % WriteValue('test message ', io_int, 2)
! call Message % WriteValue('test message ', io_int, 2, "(I10,'f----f',I10)")

! call Message % printWarning('hmmm',(/'why not ... '/))


!     Message = IO_T()


!          call Message % printWarning('EMdatapathname was not defined in the json file', &
!                      (/ 'EMDATAPATHNAME environment variable was NOT defined as a backup.', &
!                         '----> using absolute path convention                            '/) )


!     status = 999001
!     call Message % printError('EMsoftpathname was not defined in the json file', &
!                    status, (/ 'EMSOFTPATHNAME environment variable was NOT defined as a backup.' /) )



! call Message % printError('blah','more blah')

! call Message % printMessage( 'Time : '//Timing % getTimeString() )
! call Message % printMessage( 'Date : '//Timing % getDateString() )

! ! call Timing % printTimeStamp(redirect=10)

! call sleep(2) 

! call Timing % Time_tick(2)

! call sleep(1) 
! call Timing % Time_tock(2)
! call Timing % Time_tock(1)

! write (*,*) 'Interval 1 : ', Timing % getInterval(1)
! write (*,*) 'Interval 2 : ', Timing % getInterval(2)

! call Timing % Time_reset(2)

! write (*,*) 'Interval 1 : ', Timing % getInterval(1)
! write (*,*) 'Interval 2 : ', Timing % getInterval(2)

! write (*,*) 'calling Marsaglia generator for 10 random quaternions '

! qra = quat_randomArray(10, 'd', seed, northern=.TRUE.)
! call qra%quat_print() 

! qrb = quat_randomArray(10, 's', seed, northern=.TRUE.)
! call qrb%quat_print() 



end program TestProgram
