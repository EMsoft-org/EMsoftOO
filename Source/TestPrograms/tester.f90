program tester 

use mod_global 
use mod_kinds
use mod_EMsoft
use mod_symmetry
! use mod_crystallography
! use mod_QCsymmetry
! use mod_QCcrystallography
use mod_io
! use mod_dualquaternions
! use mod_octonions
use mod_HDFsupport
use HDF5
use mod_vendors
use mod_HallSG
use mod_math
use mod_platformsupport


IMPLICIT NONE 

type(EMsoft_T)          :: EMsoft
type(HDF_T)             :: HDF
type(Vendor_T)          :: VT
type(HallSG_T)          :: HSG 
type(SpaceGroup_T)      :: SSG

character(fnlen)        :: fname, groupname, inputtype, progname, progdesc, HDFstrings(10) 
integer(kind=irg)       :: hdferr, itype, istat, ipf_wd, ipf_ht, sz(3), L, recordsize, &
                           patsz, i, j, numsx, numsy, correctsize, s1, s2,HSGn, info, status
real(kind=sgl),allocatable   :: exppatarray(:), tot(:), totold(:)
real(kind=dbl),allocatable   :: SG(:,:,:)
integer(HSIZE_T)        :: dims3(3), offset3(3)
character(16)           :: HS
real(kind=dbl),allocatable          :: SGdirec(:,:,:)
real(kind=dbl)          :: z(11,11), fit(5), mp1, mp2, sig1, sig2 


status = system_hostnm(fname)
write (*,*) 'output of subroutine : ', trim(fname)

info = system_hostnm(fname)
write (*,*) 'output of function : ', trim(fname), info


stop 


mp1 = 5.D0 
mp2 = 5.D0 
sig1 = 0.5D0 
sig2 = 0.5D0 

z = reshape( (/  0.00001051,  0.00004493,  0.00012878,  0.00024743,  0.00031868,  0.00027512,  0.00015921,  0.00006176,  &
  0.00001606,  0.00000280,  0.00000033,&
  0.00010879,  0.00046518,  0.00133337,  0.00256190,  0.00329956,  0.00284861,  0.00164851,  0.00063949, &
  0.00016629,  0.00002898,  0.00000339,&
  0.00075503,  0.00322858,  0.00925424,  0.01778085,  0.02290057,  0.01977072,  0.01144144,  0.00443835, & 
  0.00115410,  0.00020116,  0.00002350,&
  0.00351266,  0.01502047,  0.04305395,  0.08272269,  0.10654143,  0.09198025,  0.05322956,  0.02064873, & 
  0.00536928,  0.00093588,  0.00010935,&
  0.01095445,  0.04684227,  0.13426640,  0.25797583,  0.33225605,  0.28684612,  0.16599969,  0.06439435, & 
  0.01674443,  0.00291861,  0.00034101,&
  0.02289956,  0.09792069,  0.28067513,  0.53928161,  0.69455954,  0.59963306,  0.34701149,  0.13461218, & 
  0.03500314,  0.00610115,  0.00071285,&
  0.03208825,  0.13721237,  0.39329890,  0.75567387,  0.97325867,  0.84024196,  0.48625340,  0.18862669, & 
  0.04904851,  0.00854930,  0.00099889,&
  0.03014026,  0.12888260,  0.36942284,  0.70979904,  0.91417489,  0.78923324,  0.45673433,  0.17717570, & 
  0.04607092,  0.00803030,  0.00093825,&
  0.01897712,  0.08114794,  0.23259852,  0.44690851,  0.57558902,  0.49692241,  0.28757218,  0.11155457, & 
  0.02900748,  0.00505609,  0.00059075,&
  0.00800932,  0.03424861,  0.09816855,  0.18861839,  0.24292819,  0.20972683,  0.12137026,  0.04708177, & 
  0.01224265,  0.00213393,  0.00024933,&
  0.00226591,  0.00968926,  0.02777282,  0.05336195,  0.06872671,  0.05933373,  0.03433681,  0.01331988, & 
   0.00346356,  0.00060371,  0.00007054 /), (/11,11/) )

z = transpose(z)

fit = (/ z(6,6), mp1, mp2, sig1, sig2 /)
call GaussianFit( 11, z, fit, info)

write (*,*) 'fit = ', fit(1), fit(2) - 5.D0, fit(3)-5.D0, fit(4), fit(5) 



stop

HS = List_Hall_Symbols(62, HSGn)

write (*,*) ' Hall Space Group number = ', HSGn, trim(HS)

SSG = SpaceGroup_T( SGnumber = 62, useHall=.TRUE., HallSGnumber=HSGn )

SGdirec = SSG%getSpaceGroupPGdirecMatrices()

sz = shape(SGdirec)

do i=1,sz(1) 
  write(*,*) 'pg matrix ',i
  do j=1,3
    write (*,*) SGdirec(i,j,1:3)
  end do 
end do

stop
HSG = HallSG_T( HS )

numsx = HSG%get_NHallgenerators()
write (*,*) 'number of generators in '//trim(HS)//' : ', numsx
allocate(SG(4,4,numsx))

SG = HSG%get_Hall_SeitzGenerators()

do i=1,numsx 
  do j=1,4
    write (*,*)  SG(j,1:4,i)
  end do 
  write (*,*) '-----'
end do 


stop
! HSG = HallSG_T( '-P 1', verbose=.TRUE. )

! HSG = HallSG_T( '-I 2xb', verbose=.TRUE. )
! HSG = HallSG_T( '-I 2xb (0 0 1)', verbose=.TRUE. )

HSG = HallSG_T( '-P 31 2c', verbose=.TRUE. )
HSG = HallSG_T( 'P 31 2c (0 0 1)', verbose=.TRUE. )

stop


HSG = HallSG_T( 'P 2 2ab -1ab', verbose=.TRUE. )
HSG = HallSG_T( 'P 4ab 2ab -1ab', verbose=.TRUE. )
HSG = HallSG_T( '-F 4 2 3', verbose=.TRUE. )
HSG = HallSG_T( 'F 4d 2 3 -1cd', verbose=.TRUE. )

stop
stop
progname = 'tester'
progdesc = 'test program to read problematic HDF5 file'
EMsoft = EMsoft_T( progname, progdesc)

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()
fname = 'playarea/Oxford/Al-large.h5'
fname = EMsoft%generateFilePath('EMdatapathname',trim(fname))

inputtype = 'TSLHDF'
VT = Vendor_T( inputtype )
itype = VT%get_itype()
call VT%set_filename(fname)

ipf_wd = 750
ipf_ht = 500
numsx = 156 
numsy = 128
L = numsx*numsy 
correctsize = 16*ceiling(float(L)/16.0)
recordsize = correctsize*4
patsz = correctsize

HDFstrings = ''
HDFstrings(1) = '1'
HDFstrings(2) = 'EBSD'
HDFstrings(3) = 'Data'
HDFstrings(4) = 'Processed Patterns'
! open the pattern file
istat = VT%openExpPatternFile(EMsoft, ipf_wd, L, recordsize, HDFstrings, HDF)

allocate(exppatarray(patsz * ipf_wd), tot(ipf_wd), totold(ipf_wd),stat=istat)

dims3 = (/ numsx, numsy, ipf_wd /)

write (*,*) L, correctsize, recordsize, patsz 

do i=250,275 ! ipf_ht
  exppatarray = 0.0
  offset3 = (/ 0, 0, (i-1)*ipf_wd /)
  call VT%getExpPatternRow(i, ipf_wd, patsz, L, dims3, offset3, exppatarray, &
                                     HDFstrings=HDFstrings, HDF=HDF)
  write (*,*) 'row number = ', i
  tot = 0.0
  do j=1,ipf_wd
    s1 = (j-1)*patsz+1
    s2 = j*patsz
    tot(j) = sum(exppatarray( s1:s2 ))
  end do
  if (i.eq.262) then 
    write (*,*) tot(19960:patsz) 
  else 
    write (*,*) tot(19960:patsz) - totold(19960:patsz)
  end if 
  write (*,*) '--------------' 
  totold = tot
end do 

call VT%closeExpPatternFile()





! type(SpaceGroup_T)          :: SG 
! type(QCSpaceGroup_T)        :: QCSG 
! type(IO_T)                  :: Message
! integer(kind=irg)           :: i, j, isg, SamplingType, isym

! type(DualQuaternion_T)      :: dq1, dq2, dqp, tmp
! real(kind=dbl)              :: v(3)

! type(Octonion_T)                :: a, b, c, d 
! real(kind=dbl)                  :: onorm
! real(kind=sgl)                  :: onorms

! a = Octonion_T( od = (/ 1.D0, 2.D0, 3.D0, 4.D0, 5.D0, 6.D0, 7.D0, 8.D0 /) )
! b = Octonion_T( od = (/ 1.D0, 3.D0, 5.D0, 7.D0, 9.D0, 2.D0, 4.D0, 6.D0 /) )
! c = Octonion_T( od = (/ 8.D0, 7.D0, 6.D0, 5.D0, 4.D0, 3.D0, 2.D0, 1.D0 /) )

! d = a+b 
! call d%oct_print('a+b = ')

! d = a-b 
! call d%oct_print('a-b = ')

! d = a*sqrt(2.D0)
! call d%oct_print('a*s = ')

! d = a*b
! call d%oct_print('a*b = ')

! d = a/b 
! call d%oct_print('a/b = ')

! onorm = cabs(a)
! write (*,*) 'a%norm = ', onorm

! d = conjg(a)
! call d%oct_print('conj(a) = ')

! d = a%octinverse()
! call d%oct_print('a%inv = ')

! a = Octonion_T( o = (/ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 /) )
! b = Octonion_T( o = (/ 1.0, 3.0, 5.0, 7.0, 9.0, 2.0, 4.0, 6.0 /) )
! c = Octonion_T( o = (/ 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 /) )

! d = a+b 
! call d%oct_print('a+b = ')

! d = a-b 
! call d%oct_print('a-b = ')

! d = a*sqrt(2.D0)
! call d%oct_print('a*s = ')

! d = a*b 
! call d%oct_print('a*b = ')

! d = a/b 
! call d%oct_print('a/b = ')

! onorms = real(cabs(a))
! write (*,*) 'a%norm = ', onorm

! d = conjg(a)
! call d%oct_print('conj(a) = ')

! d = a%octinverse()
! call d%oct_print('a%inv = ')



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