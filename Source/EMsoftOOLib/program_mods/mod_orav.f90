! ###################################################################
! Copyright (c) 2013-2024, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_orav
  !! author: MDG 
  !! version: 1.0 
  !! date: 08/19/22
  !!
  !! class definition for the EMorav program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMorav program
type, public :: oravNameListType
  character(fnlen)        :: orientationfilename
  integer(kind=irg)       :: pgnum
end type oravNameListType

! class definition
type, public :: orav_T
private 
  character(fnlen)        :: nmldeffile = 'EMorav.nml'
  type(oravNameListType)  :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: orav_
  procedure, pass(self) :: get_orientationfilename_
  procedure, pass(self) :: get_pgnum_
  procedure, pass(self) :: set_orientationfilename_
  procedure, pass(self) :: set_pgnum_
 
  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: orav => orav_
  generic, public :: get_orientationfilename => get_orientationfilename_
  generic, public :: get_pgnum => get_pgnum_
  generic, public :: set_orientationfilename => set_orientationfilename_
  generic, public :: set_pgnum => set_pgnum_


end type orav_T

! the constructor routine for this class 
interface orav_T
  module procedure orav_constructor
end interface orav_T

contains

!--------------------------------------------------------------------------
type(orav_T) function orav_constructor( nmlfile ) result(orav)
!DEC$ ATTRIBUTES DLLEXPORT :: orav_constructor
!! author: MDG 
!! version: 1.0 
!! date: 08/19/22
!!
!! constructor for the orav_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call orav%readNameList(nmlfile)

end function orav_constructor

!--------------------------------------------------------------------------
subroutine orav_destructor(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: orav_destructor
!! author: MDG 
!! version: 1.0 
!! date: 08/19/22
!!
!! destructor for the orav_T Class
 
IMPLICIT NONE

type(orav_T), INTENT(INOUT)  :: self 

call reportDestructor('orav_T')

end subroutine orav_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 08/19/22
!!
!! read the namelist from an nml file for the orav_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(orav_T), INTENT(INOUT)         :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

character(fnlen)        :: orientationfilename
integer(kind=irg)       :: pgnum

namelist / orav /  orientationfilename, pgnum 

orientationfilename = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=orav)
    close(UNIT=dataunit,STATUS='keep')

    if (trim(orientationfilename).eq.'undefined') then
        call Message%printError('readNameList:',' orientationfilename is undefined in '//nmlfile)
    end if
end if

self%nml%orientationfilename = trim(orientationfilename)
self%nml%pgnum = pgnum 

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 08/19/22
!!
!! pass the namelist for the orav_T Class to the calling program

IMPLICIT NONE 

class(orav_T), INTENT(INOUT)          :: self
type(oravNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
function get_orientationfilename_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_orientationfilename_
!! author: MDG 
!! version: 1.0 
!! date: 08/19/22
!!
!! get orientationfilename from the orav_T class

IMPLICIT NONE 

class(orav_T), INTENT(INOUT)     :: self
character(fnlen)                 :: out

out = self%nml%orientationfilename

end function get_orientationfilename_

!--------------------------------------------------------------------------
subroutine set_orientationfilename_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_orientationfilename_
!! author: MDG 
!! version: 1.0 
!! date: 08/19/22
!!
!! set orientationfilename in the orav_T class

IMPLICIT NONE 

class(orav_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)     :: inp

self%nml%orientationfilename = inp

end subroutine set_orientationfilename_

!--------------------------------------------------------------------------
function get_pgnum_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_pgnum_
!! author: MDG 
!! version: 1.0 
!! date: 08/19/22
!!
!! get pgnum from the orav_T class

IMPLICIT NONE 

class(orav_T), INTENT(INOUT)     :: self
integer(kind=irg)                :: out

out = self%nml%pgnum

end function get_pgnum_

!--------------------------------------------------------------------------
subroutine set_pgnum_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_pgnum_
!! author: MDG 
!! version: 1.0 
!! date: 08/19/22
!!
!! set pgnum in the orav_T class

IMPLICIT NONE 

class(orav_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)    :: inp

self%nml%pgnum = inp

end subroutine set_pgnum_

!--------------------------------------------------------------------------
subroutine orav_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: orav_
!! author: MDG 
!! version: 1.0 
!! date: 08/19/22
!!
!! perform the computations

use mod_EMsoft
use mod_HDFnames
use mod_dirstats 
use mod_quaternions
use mod_so3
use mod_rotations
use mod_memory
use mod_io

IMPLICIT NONE 

class(orav_T), INTENT(INOUT)            :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 

type(DirStat_T)                         :: dictVMF, dictWAT
type(so3_T)                             :: SO
type(memory_T)                          :: mem
type(IO_T)                              :: Message
type(QuaternionArray_T)                 :: qAR, QA, qsym
type(Quaternion_T)                      :: muhat
type(r_T)                               :: rod 
type(q_T)                               :: quat

integer(kind=irg)                       :: pgnum, numors, ii, io_int(1), seed
type(FZpointd),pointer                  :: FZlist, FZtmp
real(kind=dbl),allocatable              :: samples(:,:)
real(kind=dbl)                          :: qav(4), qstdev(4), io_dbl(4), kappahat
character(fnlen)                        :: fname

associate( nml => self%nml )

! rotation precision
call setRotationPrecision('Double')

! memory management class 
mem = memory_T()

! handle the point group number
pgnum = self%get_pgnum() 
if (nml%pgnum.eq.0) pgnum = 1 

! get the symmetry quaternions 
call QA%QSym_init( pgnum, qsym )

! get the input orientations
SO = so3_T( pgnum )
call SO%setFZtypeandorder( pgnum )
fname = EMsoft%generateFilePath('EMdatapathname',trim(nml%orientationfilename))
call SO%getOrientationsfromFile(fname)
if (pgnum.gt.1) call SO%ReducelisttoRFZ( qsym ) ! reduce to the correct fundamental zone
numors = SO%getListCount('FZ')
FZlist => SO%getListHead('FZ')
io_int(1) = numors
call Message%WriteValue(' Number of orientations read : ', io_int, 1, "(I8)")

! copy the orientations into a quaternion array qAR
call mem%alloc(samples, (/ 4, numors /), 'samples')
FZtmp => Fzlist
do ii = 1,numors
  rod = r_T( rdinp = FZtmp%rod%r_copyd() )
  quat = rod%rq()
  samples(1:4,ii) = quat%q_copyd()
  FZtmp => FZtmp%next
end do
qAR = QuaternionArray_T( n=numors, qd=samples )
! we can now delete the linked list since we have the quaternions in an array
call SO%delete_FZlist()

! initialize the directional statistics classes for von Mises-Fisher and Watson distributions
dictVMF = DirStat_T( DStype='VMF', pgnum = pgnum)
call dictVMF%setNumEM(15)
call dictVMF%setNumIter(40)

dictWAT = DirStat_T( DStype='WAT', pgnum = pgnum)
call dictWAT%setNumEM(15)
call dictWAT%setNumIter(40)

! next we need to store the qAR array into both dict classes and clean up qAR 
call dictVMF%setQuatArray( qAR )
call dictWAT%setQuatArray( qAR )
! if (allocated(qAR%qd)) deallocate(qAR%qd)

! now we are ready to compute the orientation average using four different approaches:
! - quaternion logarithm
! - quaternion eigenvalue (T-matrix)
! - von Mises-Fisher averaging
! - Watson averaging 

! quaternion logarithm
qav = quat_average( samples, numors, qstdev )
call Message%printMessage(' Quaternion logarithmic average [CAUTION]')
io_dbl = qav 
call Message%WriteValue(' <q>      ', io_dbl, 4)
io_dbl = qstdev
call Message%WriteValue(' st. dev. ', io_dbl, 4)

! quaternion T-matrix 
qav = quat_average( samples, numors, qstdev, Tmatrix=.TRUE. )
call Message%printMessage(' ----------- ')
call Message%printMessage(' Quaternion T-matrix average [CAUTION]')
io_dbl = qav 
call Message%WriteValue(' <q> ', io_dbl, 4)
io_dbl = qstdev
call Message%WriteValue(' st. dev. ', io_dbl, 4)

! VMF averaging
seed = 43514
muhat = Quaternion_T( qd=(/ 1.D0, 0.D0, 0.D0,0.D0 /) )
call dictVMF%EMforDS( seed, muhat, kappahat )
call Message%printMessage(' ----------- ')
call Message%printMessage(' Quaternion von Mises-Fisher average')
io_dbl = muhat%get_quatd()
call Message%WriteValue(' <q>      ', io_dbl, 4)
io_dbl(1) = kappahat
call Message%WriteValue(' kappa    ', io_dbl, 1)
io_dbl(1) = 180.D0*dacos(1.D0-1.D0/kappahat)/cPi 
call Message%WriteValue(' eq. deg. ', io_dbl, 1)

! WAT averaging
seed = 43514
muhat = Quaternion_T( qd=(/ 1.D0, 0.D0, 0.D0,0.D0 /) )
call dictWAT%EMforDS( seed, muhat, kappahat )
call Message%printMessage(' ----------- ')
call Message%printMessage(' Quaternion Watson average')
io_dbl = muhat%get_quatd()
call Message%WriteValue(' <q>      ', io_dbl, 4)
io_dbl(1) = kappahat
call Message%WriteValue(' kappa    ', io_dbl, 1)
io_dbl(1) = 180.D0*dacos(1.D0-1.D0/kappahat)/cPi 
call Message%WriteValue(' eq. deg. ', io_dbl, 1)

call Message%printMessage(' ')
call Message%printMessage(' ----------- ')
call Message%printMessage(' ')
call Message%printMessage((/ ' CAUTION: if the input orientation list straddles the fundamental zone ', & 
                             '          boundary, then the quaternion averaging procedure will not   ', &
                             '          generate the correct result.  In that case, only averaging on', &
                             '          the quaternion sphere will produce the correct average.      ', &
                             '                                                                       ' /))

! delete samples
call mem%dealloc(samples, 'samples')

end associate 

end subroutine orav_



end module mod_orav