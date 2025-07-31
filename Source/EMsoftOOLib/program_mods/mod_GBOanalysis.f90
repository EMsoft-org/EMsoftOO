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

module mod_GBOanalysis
  !! author: MDG 
  !! version: 1.0 
  !! date: 07/22/25
  !!
  !! class definition for the EMGBOanalysis program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMGBOanalysis program
type, public :: GBOanalysisNameListType
  integer(kind=irg)     :: pgnum 
  character(3)          :: Component
  real(kind=dbl)        :: oct(8)
end type GBOanalysisNameListType

! class definition
type, public :: GBOanalysis_T
private 
  character(fnlen)       :: nmldeffile = 'EMGBOanalysis.nml'
  type(GBOanalysisNameListType)  :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: GBOanalysis_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: GBOanalysis => GBOanalysis_

end type GBOanalysis_T

! the constructor routine for this class 
interface GBOanalysis_T
  module procedure GBOanalysis_constructor
end interface GBOanalysis_T

contains

!--------------------------------------------------------------------------
type(GBOanalysis_T) function GBOanalysis_constructor( nmlfile ) result(GBOanalysis)
!! author: MDG 
!! version: 1.0 
!! date: 07/22/25
!!
!! constructor for the GBOanalysis_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call GBOanalysis%readNameList(nmlfile)

end function GBOanalysis_constructor

!--------------------------------------------------------------------------
subroutine GBOanalysis_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 07/22/25
!!
!! destructor for the GBOanalysis_T Class
 
IMPLICIT NONE

type(GBOanalysis_T), INTENT(INOUT)  :: self 

call reportDestructor('GBOanalysis_T')

end subroutine GBOanalysis_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 07/22/25
!!
!! read the namelist from an nml file for the GBOanalysis_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(GBOanalysis_T), INTENT(INOUT)   :: self
character(fnlen),INTENT(IN)           :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)           :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                        :: EMsoft 
type(IO_T)                            :: Message       
logical                               :: skipread = .FALSE.

integer(kind=irg)                     :: pgnum 
character(3)                          :: Component
real(kind=dbl)                        :: oct(8)

namelist /GBOanalysis/ pgnum, Component, oct

pgnum = 32 
Component = 'All'
oct = (/ 1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0 /)

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=GBOanalysis)
 close(UNIT=dataunit,STATUS='keep')
end if 

self%nml%pgnum = pgnum 
self%nml%Component = Component
self%nml%oct = oct 

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 07/22/25
!!
!! pass the namelist for the GBOanalysis_T Class to the calling program

IMPLICIT NONE 

class(GBOanalysis_T), INTENT(INOUT)          :: self
type(GBOanalysisNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine GBOanalysis_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: GBOanalysis_
!! author: MDG 
!! version: 1.0 
!! date: 07/22/25
!!
!! perform the computations

use mod_EMsoft
use mod_octonions
use mod_GBoctonions
use mod_dirstats
use mod_quaternions
use mod_io
use mod_symmetry 

IMPLICIT NONE 

class(GBOanalysis_T), INTENT(INOUT)     :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 

type(DirStat_T)                         :: DS
type(octonion_T)                        :: oct 
type(GBoctonion_T)                      :: GBoct 
type(GBoctonionArray_T)                 :: GBO_equiv 
type(QuaternionArray_T)                 :: qsym 
type(IO_T)                              :: Message
type(Quaternion_T)                      :: qa, qb

integer(kind=irg)                       :: Nqsym, io_int(1), NBflag
logical                                 :: enantiomorphic, centrosymmetric
real(kind=dbl)                          :: diff, io_real(1), epsd=1.0D-12  


associate( nml => self%nml )

! turn the input octonion into a properly normalized GBoctonion_T class 
oct = octonion_T( od = nml%oct )
GBoct = GBoctonion_T( oct = oct )

! first make sure that the two member quaternions are different; if 
! they are the same, then there is effectively no grain boundary and 
! we treat this case separately.
qa = GBoct%GBO_get_q(1)
qb = GBoct%GBO_get_q(2)
diff = sum( abs(qa%get_quatd() - qb%get_quatd()) )
if (diff.lt.epsd) then 
  io_real(1) = diff
  call Message%WriteValue(' difference between member quaternions : ', io_real, 1)
  call qa%quat_print(' qa: ')
  call qb%quat_print(' qb: ')
  call Message%printMessage(' so there is effectively no grain boundary here ... ')
  call Message%printMessage(' ')
  stop 'All is well that ends well ... [Shakespeare, 1623]'
end if 
 
! initialize the directional statistics class
DS = DirStat_T( pgnum = nml%pgnum )
qsym = DS%getQuatArray(slot='qsym')
Nqsym = qsym%getQnumber()

io_int(1) = Nqsym
call Message%WriteValue( ' Number of symmetry operators ', io_int, 1)

! we'll need to distinguish between the purely rotation (enantiomorphic) point
! groups and the centrosymmetric point groups 

! is this an enantiomorphic point group ?
enantiomorphic = .FALSE.
if (PGrot(nml%pgnum).eq.nml%pgnum) then 
  enantiomorphic = .TRUE.
  call Message%printMessage(' This point group ('//trim(adjustl(PGTHD(nml%pgnum)))//') is enantiomorphic. ')
else
  call Message%printMessage(' This point group ('//trim(adjustl(PGTHD(nml%pgnum)))//') is not enantiomorphic. ')
end if 

! is this a centrosymmetric point group ?
centrosymmetric = .FALSE.
if (PGLaue(nml%pgnum).eq.nml%pgnum) then 
  centrosymmetric = .TRUE.
  call Message%printMessage(' This point group ('//trim(adjustl(PGTHD(nml%pgnum)))//') is centrosymmetric. ')
else
  call Message%printMessage(' This point group ('//trim(adjustl(PGTHD(nml%pgnum)))//') is not centrosymmetric. ')
end if 

! check for the 10 point groups that have neither property and abort the program 
! for them since they are, as far as we know, incompatible with the octonion framework
if ( (.not.centrosymmetric).and.(.not.enantiomorphic)) then 
  call Message%printMessage( (/ ' This point group belongs to the set of 10 point groups that cannot be', &
                                ' handled by means of octonions (as far as we know).  This approach can', &
                                ' only be carried out using symmetry matrices, not quaternions nor     ', &
                                ' octonions. Hence, there is no point in continuing this program...    ' /)) 

  stop 'All is well that ends well ... [Shakespeare, 1623]'
end if 

! get the symmetrically equivalent GB octonions 
GBO_equiv = GBoct%GBO_get_equivalent(qsym, nthreads=1, NBflag=NBflag)
io_int(1) = GBO_equiv%getOnumber()
call Message%WriteValue(' Total number of unique equivalent GB octonions :', io_int, 1)
call GBO_equiv%oct_arrayprint()

if (NBflag.ne.0) then 
  io_int(1) = NBflag
  call Message%WriteValue(' The following octonion is a No Boundary octonion :', io_int, 1)
  call Message%printMessage(' so basically all of them are No Boundary octonions... ')
  call Message%printMessage(' Nothing further to do ... ')
  call Message%printMessage(' ')
  stop 'All is well that ends well ... [Shakespeare, 1623]'
end if 


end associate

end subroutine GBOanalysis_



end module mod_GBOanalysis