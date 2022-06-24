program tester 

use mod_global 
use mod_kinds
use mod_symmetry
use mod_crystallography
use mod_QCsymmetry
use mod_QCcrystallography
use mod_io

IMPLICIT NONE 

type(SpaceGroup_T)          :: SG 
type(QCSpaceGroup_T)        :: QCSG 
type(IO_T)                  :: Message
integer(kind=irg)           :: i, j, isg, SamplingType, isym


! test for the 230 space groups to ensure that all symmetry codes are correct
! both for point groups and k-vector sampling 

open(unit=dataunit, file='sym-new.txt', status='unknown',form='formatted')

do isg=1,230 
  call SG%setSpaceGroupNumber(isg)
  call SG%setSpaceGrouptrigonal(.FALSE.)
  if ((isg.ge.143).and.(isg.le.167)) call SG%setSpaceGrouptrigonal(.TRUE.)

  j=0
  do i=1,32
    if (SGPG(i).le.isg) j=i
  end do
  isym = j
  SamplingType = PGSamplingType(isym)

! next, intercept the special cases (hexagonal vs. rhombohedral cases that require special treatment)
  if ((SamplingType.eq.-1).or.(isym.eq.14).or.(isym.eq.26)) then
    SamplingType = SG%getHexvsRho(isym)
  end if
  write (dataunit,"(3I6)") isg, isym, SamplingType
end do 

close(unit=dataunit, status='keep')


stop 

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