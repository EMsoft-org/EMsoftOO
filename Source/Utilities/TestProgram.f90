program TestProgram

use mod_kinds
use mod_global
use mod_EMsoft
use mod_io
use mod_timing
use mod_quaternions 
use mod_rng

IMPLICIT NONE

character(fnlen)        :: progname = 'this is the program name'
character(fnlen)        :: progdesc = 'and this is the descriptor'
character(fnlen)        :: m

type(EMsoft_T)     :: EMsoft 
type(IO_T)         :: Message
type(Timing_T)     :: Timing

integer(kind=irg)       :: io_int(2), status
type(rng_t)             :: seed 
type(QuaternionArray_T) :: qra , qrb


Timing = Timing_T( showDateTime = .TRUE. )

call Timing % Time_tick()

EMsoft = EMsoft_T(progname, progdesc, showconfig=.TRUE., makeconfig=.FALSE.) 

call EMsoft%setConfigParameter('WyckoffPositionsfilename','no idea where this file is located ')

write (*,*) 'modified string = ',trim(EMsoft%getConfigParameter('WyckoffPositionsfilename'))




Message = IO_T()

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

call Message % printMessage( 'Time : '//Timing % getTimeString() )
call Message % printMessage( 'Date : '//Timing % getDateString() )

! call Timing % printTimeStamp(redirect=10)

call sleep(2) 

call Timing % Time_tick(2)

call sleep(1) 
call Timing % Time_tock(2)
call Timing % Time_tock(1)

write (*,*) 'Interval 1 : ', Timing % getInterval(1)
write (*,*) 'Interval 2 : ', Timing % getInterval(2)

call Timing % Time_reset(2)

write (*,*) 'Interval 1 : ', Timing % getInterval(1)
write (*,*) 'Interval 2 : ', Timing % getInterval(2)

write (*,*) 'calling Marsaglia generator for 10 random quaternions '

qra = quat_randomArray(10, 'd', seed, northern=.TRUE.)
call qra%quat_print() 

qrb = quat_randomArray(10, 's', seed, northern=.TRUE.)
call qrb%quat_print() 



end program TestProgram