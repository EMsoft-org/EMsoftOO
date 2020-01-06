program EMmkxtal

use mod_kinds
use mod_global
use mod_EMsoft
use mod_io
use mod_timing

IMPLICIT NONE

character(fnlen)        :: progname = 'this is the program name'
character(fnlen)        :: progdesc = 'and this is the descriptor'
character(fnlen)        :: m

type(EMsoft_T)     :: EMsoft 
type(T_IOClass)         :: Message
type(T_TimingClass)     :: Timing

integer(kind=irg)       :: io_int(2), status

Timing = T_TimingClass( showDateTime = .TRUE. )

call Timing % Time_tick()

EMsoft = EMsoft_T(progname, progdesc, showconfig=.TRUE., makeconfig=.FALSE.) 

Message = T_IOClass()

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


!     Message = T_IOClass()


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



end program EMmkxtal
