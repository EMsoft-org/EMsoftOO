program tester 

use mod_global 
use mod_kinds
use mod_QCsymmetry
use mod_QCcrystallography
use mod_io

IMPLICIT NONE 

type(QCSpaceGroup_T)        :: QCSG 
type(IO_T)                  :: Message
integer(kind=irg)           :: i



do i=6,11 
  if (i.ne.7) then 
    QCSG = QCspacegroup_T( nD = 3, QCtype = 'Ico')
    call QCSG%setSGnum(i)
! call QCcell_icosahedral%setMetricParametersQC()
    call QCSG%GenerateQCSymmetry(dopg=.FALSE.) 
    write (*,*) ' order of icosahedral space group ', i, ' equals ',QCSG%getnsym()
  end if 
end do

end program