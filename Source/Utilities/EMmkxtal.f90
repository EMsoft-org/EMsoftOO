program EMmkxtal

use mod_kinds
use mod_global
use mod_EMsoft
use mod_io

IMPLICIT NONE

character(fnlen)        :: progname = 'this is the program name'
character(fnlen)        :: progdesc = 'and this is the descriptor'

type(T_EMsoftClass)     :: EMsoft 
type(T_IOClass)         :: Message
character(fnlen)        :: m

integer(kind=irg)       :: io_int(2)


EMsoft = T_EMsoftClass(progname, progdesc, showconfig=.TRUE.) 
Message = T_IOClass()

write (*,*) 'Revision = ', trim( EMsoft % getConfigParameter('EMsoftRevision') )

m = trim( EMsoft % generateFilePath('EMsofttestpath','test.h5'))
write (*,*) 'Completed path = ', trim(m)
! m = trim( EMsoft % generateFilePath('blahblah','test.h5')) 
! write (*,*) 'Completed path = ', trim(m)
m = trim( EMsoft % generateFilePath('Templatecodefilename','')) 
write (*,*) 'Completed path = ', trim(m)


io_int = (/ 10, 20 /)
call Message % WriteValue('test message ', io_int, 2)
call Message % WriteValue('test message ', io_int, 2, "(I10,'f----f',I10)")

call Message % printWarning('hmmm','why not ... ')
call Message % printError('blah','more blah')

end program EMmkxtal
