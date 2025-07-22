! ###################################################################
! Copyright (c) 2013-2025, Marc De Graef Research Group/Carnegie Mellon University
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
!     - Redistributions of source code must retain the above copyright notice, this list 
!        of conditions and the following disclaimer.
!     - Redistributions in binary form must reproduce the above copyright notice, this 
!        list of conditions and the following disclaimer in the documentation and/or 
!        other materials provided with the distribution.
!     - Neither the names of Marc De Graef, Carnegie Mellon University nor the names 
!        of its contributors may be used to endorse or promote products derived from 
!        this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
! USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! ###################################################################

module mod_GBO
  !! author: MDG 
  !! version: 1.0 
  !! date: 07/21/25
  !!
  !! class definition for the EMGBO program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMGBO program
type, public :: GBONameListType
  integer(kind=irg)       :: pgnum
  integer(kind=irg)       :: numsamples
  integer(kind=irg)       :: numbins
  integer(kind=irg)       :: nthreads
  integer(kind=irg)       :: seed
  character(3)            :: CSLtype
  logical                 :: fixedAB
  character(fnlen)        :: outname
  character(fnlen)        :: octonions
end type GBONameListType

! class definition
type, public :: GBO_T
private 
  character(fnlen)       :: nmldeffile = 'EMGBO.nml'
  type(GBONameListType)  :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: GBO_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: GBO => GBO_

end type GBO_T

! the constructor routine for this class 
interface GBO_T
  module procedure GBO_constructor
end interface GBO_T

contains

!--------------------------------------------------------------------------
type(GBO_T) function GBO_constructor( nmlfile ) result(GBO)
!! author: MDG 
!! version: 1.0 
!! date: 07/21/25
!!
!! constructor for the GBO_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call GBO%readNameList(nmlfile)

end function GBO_constructor

!--------------------------------------------------------------------------
subroutine GBO_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 07/21/25
!!
!! destructor for the GBO_T Class
 
IMPLICIT NONE

type(GBO_T), INTENT(INOUT)  :: self 

call reportDestructor('GBO_T')

end subroutine GBO_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 07/21/25
!!
!! read the namelist from an nml file for the GBO_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(GBO_T), INTENT(INOUT)           :: self
character(fnlen),INTENT(IN)           :: nmlfile
logical,OPTIONAL,INTENT(IN)           :: initonly

type(EMsoft_T)                        :: EMsoft 
type(IO_T)                            :: Message       
logical                               :: skipread = .FALSE.

integer(kind=irg)                     :: pgnum
integer(kind=irg)                     :: numsamples
integer(kind=irg)                     :: numbins
integer(kind=irg)                     :: nthreads
integer(kind=irg)                     :: seed
character(3)                          :: CSLtype
logical                               :: fixedAB
character(fnlen)                      :: outname
character(fnlen)                      :: octonions

namelist /GBOlist/ pgnum, numsamples, numbins, outname, nthreads, CSLtype, fixedAB, octonions, seed

nthreads = 1
outname = 'undefined' 
octonions = 'random'
pgnum = 32
numsamples = 100000
numbins = 180
CSLtype = ''
fixedAB = .FALSE.
seed = 543254

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=GBOlist)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(outname).eq.'undefined') then
  call Message%printError('EMGBO:',' output file name is undefined in '//nmlfile)
 end if
end if

self%nml%nthreads = nthreads
self%nml%pgnum = pgnum
self%nml%numsamples = numsamples
self%nml%numbins = numbins
self%nml%outname = trim(outname)
self%nml%octonions = trim(octonions)
self%nml%CSLtype = trim(CSLtype)
self%nml%fixedAB = fixedAB
self%nml%seed = seed

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 07/21/25
!!
!! pass the namelist for the GBO_T Class to the calling program

IMPLICIT NONE 

class(GBO_T), INTENT(INOUT)          :: self
type(GBONameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine GBO_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: GBO_
!! author: MDG 
!! version: 1.0 
!! date: 07/21/25
!!
!! perform the computations

use mod_EMsoft
use mod_io 
use mod_quaternions
use mod_GBoctonions
use mod_rotations
use mod_OMPsupport
use mod_math
use mod_dirstats
use mod_CSL

IMPLICIT NONE 

class(GBO_T), INTENT(INOUT)             :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 

type(IO_T)                              :: Message
type(Octonion_T)                        :: o
type(OctonionArray_T)                   :: octarray 
type(GBoctonion_T)                      :: oct1, oct2
type(Quaternion_T)                      :: qu, q1, q2, q3, q4
type(QuaternionArray_T)                 :: qsym
type(DirStat_T)                         :: DS
type(CSL_T)                             :: CSL
type(r_T)                               :: ro 
type(q_T)                               :: qa, qb, qc, qd
type(c_T)                               :: cu

integer(kind=irg)                       :: numb, ip, TID, i, j, io_int(2), CSLnumber, numoct, Nqsym
integer(kind=irg),allocatable           :: histogramSYM(:), histogramNBSYM(:)
integer(kind=irg)                       :: seed, myseed
real(kind=dbl)                          :: aa(3),bb(3),cc(3),dd(3),pp,tt, x(12), y(6), acube,&
                                           x1,x2,y1,y2,s1,s2, eu1(3), eu2(3), scale, Sqa(4), Sqc(4), Sqb(4), Sqd(4), &
                                           oac, obd, oo(8), qq(4), oct(8)
character(fnlen)                        :: fname 
character(4)                            :: mode
logical                                 :: f_exists  

associate( nml => self%nml )


! what is the octonion generator mode?  'random' or a user-provided text file
if (trim(nml%octonions).eq.'random') then 
  mode = 'rand'
  numoct = nml%numsamples
  call Message%printMessage(' Octonion mode:  random generation based on cubochoric sampling')
else
  mode = 'file'
  fname = EMsoft%generateFilePath('EMdatapathname',trim(nml%octonions))
! does this file exist ?
  inquire(file=trim(fname), exist=f_exists)
  if (.not.f_exists) then 
    call Message%printError('EMGBO','input octonion file does not exist')
  end if
  call Message%printMessage(' Octonion mode:  reading octonion pairs from '//trim(fname))
  open(unit=dataunit,file=trim(fname),status='old',form='formatted')
  read(dataunit,*) numoct
  octarray = OctonionArray_T( n = 2*numoct, s='d', nthreads=nml%nthreads )
  do i=1,2*numoct
    read(dataunit,*) oct
    o = Octonion_T( od = oct )
    call octarray%insertOctinArray(i, o)
  end do
  close(dataunit,status='keep')
  if (trim(nml%CSLtype).ne.'') numoct = 2*numoct
end if

! set the number of bins for the interval [0,180]
numb = nml%numbins
scale = (180.D0/cPi) * (float(numb)/180.D0)

! generate the histograms
allocate(histogramSYM(0:numb), histogramNBSYM(0:numb))
histogramSYM = 0
histogramNBSYM = 0

acube = 0.5D0 * cPi**0.666666666

! initialize the directional statistics class
DS = DirStat_T( pgnum = nml%pgnum )
qsym = DS%getQuatArray(slot='qsym')
Nqsym = qsym%getQnumber()

io_int(1) = Nqsym
call Message%WriteValue( ' Number of symmetry operators ', io_int, 1)

! initialize the CSL class 
CSL = CSL_T()

if (trim(nml%CSLtype).ne.'') then
    call Message%printMessage(' Available CSL Boundary quaternions : ')
    do i=1,CSL%CSLnum
      ro = CSL%getCSLrod( CSL%CSLlabels(i), CSLnumber, qu )
      write (*,"(I3,'  ',A3,'  ',4F10.6)") CSLnumber, CSL%CSLlabels(i), qu%get_quatd()
    end do
end if

! initialize the random number generator
seed = nml%seed

if (trim(nml%CSLtype).eq.'') then
!open(unit=20,file='testresults.txt',status='unknown',form='formatted')
  call OMP_SET_NUM_THREADS(nml%nthreads)
  io_int(1) = nml%nthreads
  call Message%WriteValue(' -> Number of threads set to ',io_int,1,"(I3)")

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID,aa,bb,cc,dd,pp,qa,qb,qc,qd,o,qq,tt,ip,i,j,myseed,x,x1,x2,y1,y2,s1,s2,oct1,oct2) 

  TID = OMP_GET_THREAD_NUM()

  ! make sure each thread has a different initial seed ... 
  myseed = abs( seed * (TID+1) ) 

!$OMP DO SCHEDULE(DYNAMIC)
  do i=1,numoct
    if (mode.eq.'rand') then 
! get four random cubochoric points
      call r8vec_uniform_01(12,myseed,x)
      aa = acube * (/ 2.D0*x(1)-1.D0, 2.D0*x(8)-1.D0, 2.D0*x(9)-1.D0 /)
      bb = acube * (/ 2.D0*x(2)-1.D0, 2.D0*x(7)-1.D0, 2.D0*x(10)-1.D0 /)
      cc = acube * (/ 2.D0*x(3)-1.D0, 2.D0*x(6)-1.D0, 2.D0*x(11)-1.D0 /)
      dd = acube * (/ 2.D0*x(4)-1.D0, 2.D0*x(5)-1.D0, 2.D0*x(12)-1.D0 /)
      cu = c_T( cdinp = aa )
      qa = cu%cq()
      cu = c_T( cdinp = bb )
      qb = cu%cq()
      cu = c_T( cdinp = cc )
      qc = cu%cq()
      cu = c_T( cdinp = dd )
      qd = cu%cq()
      oct1 = GBoctonion_T( Quaternion_T( qd=qa%q_copyd() ), Quaternion_T( qd=qb%q_copyd() ) )
      oct2 = GBoctonion_T( Quaternion_T( qd=qc%q_copyd() ), Quaternion_T( qd=qd%q_copyd() ) )
    else
! get 2 octonions from the octarray
      j = 2*(i-1)+1
      o = octarray%getOctfromArray(j) 
      oct1 = GBoctonion_T( oct = o)
      o = octarray%getOctfromArray(j+1) 
      oct2 = GBoctonion_T( oct = o)
    end if

! GBO geodesic distance with symmetry
    if (nml%fixedAB.eqv..TRUE.) then
      tt = oct1%GBO_Omega_symmetric(oct2, DS, single=.TRUE.)
    else
      tt = oct1%GBO_Omega_symmetric(oct2, DS)
    end if
    ip = nint(tt*scale)
    if (ip.lt.numb) histogramSYM(ip) = histogramSYM(ip) + 1

! No-Boundary GBO geodesic distance
    tt = oct1%GBO_Omega_symmetric_NB(oct2, Nqsym, qsym)
    ip = nint(tt*scale)
    if (ip.lt.numb) histogramNBSYM(ip) = histogramNBSYM(ip) + 1

    if (mod(i,5000).eq.0) then
      io_int(1) = i 
      io_int(2) = nml%numsamples
      call Message%WriteValue('completed ',io_int,2,"(I14,' out of ',I14)")
    end if
  end do
!$OMP END DO
!$OMP END PARALLEL
else
! we're doing a CSL boundary so let's get the correct quaternion for it ... 
  q1 = Quaternion_T( qd=(/ 1.D0, 0.D0, 0.D0, 0.D0 /) )
  ro = CSL%getCSLrod( nml%CSLtype, CSLnumber, qu )
  oct1 = GBoctonion_T( Quaternion_T( qd=q1%get_quatd() ), Quaternion_T( qd=qu%get_quatd() ) )
  
  call OMP_SET_NUM_THREADS(nml%nthreads)
  io_int(1) = nml%nthreads
  call Message%WriteValue(' -> Number of threads set to ',io_int,1,"(I3)")

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID, cc, dd, pp, qc, qd, qq, tt, ip, i,j, x1,x2,y1,y2,s1,s2,myseed,oct2) 

  TID = OMP_GET_THREAD_NUM()

! make sure each thread has a different initial seed ... 
  myseed = abs( seed * (TID+1) ) 

  ! generate random quartets of quaternions
!$OMP DO SCHEDULE(DYNAMIC)
  do i=1,numoct
    if (mode.eq.'rand') then 
      call r8vec_uniform_01(6,myseed,y)
! get four random cubochoric points
      cc = acube * (/ 2.D0*y(1)-1.D0, 2.D0*y(4)-1.D0, 2.D0*y(5)-1.D0 /)
      dd = acube * (/ 2.D0*y(2)-1.D0, 2.D0*y(3)-1.D0, 2.D0*y(6)-1.D0 /)
! generate the unit quaternions with positive scalar part
      qb = cu%cq()
      cu = c_T( cdinp = cc )
      qc = cu%cq()
      cu = c_T( cdinp = dd )
      oct2 = GBoctonion_T( Quaternion_T( qd=qc%q_copyd() ), Quaternion_T( qd=qd%q_copyd() ) )
    else 
      o = octarray%getOctfromArray(i) 
      oct2 = GBoctonion_T( oct = o)
    end if 

  ! GBO geodesic distance with symmetry
   if (nml%fixedAB.eqv..TRUE.) then
      tt = oct1%GBO_Omega_symmetric(oct2, DS, single=.TRUE.)
    else
      tt = oct1%GBO_Omega_symmetric(oct2, DS)
    end if
    ip = nint(tt*scale)
    if (ip.lt.numb) histogramSYM(ip) = histogramSYM(ip) + 1

  ! No-Boundary GBO geodesic distance
    tt = oct1%GBO_Omega_symmetric_NB(oct2, Nqsym, qsym)
    ip = nint(tt*scale)
    if (ip.lt.numb) histogramNBSYM(ip) = histogramNBSYM(ip) + 1

    if (mod(i,5000).eq.0) then
      io_int(1) = i 
      io_int(2) = nml%numsamples
      call Message%WriteValue('completed ',io_int,2,"(I14,' out of ',I14)")
    end if
  end do
!$OMP END DO
!$OMP END PARALLEL
end if

! output results
fname = EMsoft%generateFilePath('EMdatapathname',trim(nml%outname))

open(unit=dataunit,file=trim(fname),status='unknown',form='formatted')
write (dataunit,"(I6)") numb
do i=0,numb
  write (dataunit,"(I10,' ',F12.8,' ',F12.8)") i, real(histogramSYM(i))/real(numoct), &
                                               real(histogramNBSYM(i))/real(numoct)
end do 
close(unit=dataunit,status='keep')


end associate

end subroutine GBO_



end module mod_GBO