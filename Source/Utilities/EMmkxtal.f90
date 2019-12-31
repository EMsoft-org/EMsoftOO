program EMmkxtal

use mod_kinds
use mod_EMsoft

IMPLICIT NONE

character(fnlen)        :: progname, progdesc 

type(T_EMsoftClass)     :: EMsoft 

progname = 'this is the program name'
progdesc = 'and this is the descriptor'

EMsoft = T_EMsoftClass(progname, progdesc) 

end program EMmkxtal
