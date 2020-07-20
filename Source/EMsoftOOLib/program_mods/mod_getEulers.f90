! ###################################################################
! Copyright (c) 2013-2020, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_getEulers
  !! author: MDG
  !! version: 1.0
  !! date: 04/08/20
  !!
  !! class definition for the EMgetEulers program

use mod_kinds
use mod_global

IMPLICIT NONE

! namelist for the EMgetEulers program
type, public :: getEulersNameListType
  character(8)      :: angledataset
  character(3)      :: raddeg
  character(fnlen)  :: txtfile
  character(fnlen)  :: datafile
  character(fnlen)  :: EMEBSDnmlfile
  character(fnlen)  :: dotproductfile
end type getEulersNameListType

! class definition
type, public :: getEulers_T
private
  character(fnlen)       :: nmldeffile = 'EMgetEulers.nml'
  type(getEulersNameListType)  :: nml

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: getEulers_
  procedure, pass(self) :: get_angledataset_
  procedure, pass(self) :: get_raddeg_
  procedure, pass(self) :: get_txtfile_
  procedure, pass(self) :: get_datafile_
  procedure, pass(self) :: get_EMEBSDnmlfile_
  procedure, pass(self) :: get_dotproductfile_
  procedure, pass(self) :: set_angledataset_
  procedure, pass(self) :: set_raddeg_
  procedure, pass(self) :: set_txtfile_
  procedure, pass(self) :: set_datafile_
  procedure, pass(self) :: set_EMEBSDnmlfile_
  procedure, pass(self) :: set_dotproductfile_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: getEulers => getEulers_
  generic, public :: get_angledataset => get_angledataset_
  generic, public :: get_raddeg => get_raddeg_
  generic, public :: get_txtfile => get_txtfile_
  generic, public :: get_datafile => get_datafile_
  generic, public :: get_EMEBSDnmlfile => get_EMEBSDnmlfile_
  generic, public :: get_dotproductfile => get_dotproductfile_
  generic, public :: set_angledataset => set_angledataset_
  generic, public :: set_raddeg => set_raddeg_
  generic, public :: set_txtfile => set_txtfile_
  generic, public :: set_datafile => set_datafile_
  generic, public :: set_EMEBSDnmlfile => set_EMEBSDnmlfile_
  generic, public :: set_dotproductfile => set_dotproductfile_

end type getEulers_T

! the constructor routine for this class
interface getEulers_T
  module procedure getEulers_constructor
end interface getEulers_T

contains

!--------------------------------------------------------------------------
type(getEulers_T) function getEulers_constructor( nmlfile ) result(getEulers)
!DEC$ ATTRIBUTES DLLEXPORT :: getEulers_constructor
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! constructor for the getEulers_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call getEulers%readNameList(nmlfile)

end function getEulers_constructor

!--------------------------------------------------------------------------
subroutine getEulers_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: getEulers_destructor
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! destructor for the getEulers_T Class

IMPLICIT NONE

type(getEulers_T), INTENT(INOUT)  :: self

call reportDestructor('getEulers_T')

end subroutine getEulers_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! read the namelist from an nml file for the getEulers_T Class

use mod_io
use mod_EMsoft

IMPLICIT NONE

class(getEulers_T), INTENT(INOUT)    :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft
type(IO_T)                           :: Message
logical                              :: skipread = .FALSE.

character(8)            :: angledataset   ! 'original' or 'refined'
character(3)            :: raddeg         ! 'rad' or 'deg'
character(fnlen)        :: txtfile
character(fnlen)        :: datafile
character(fnlen)        :: EMEBSDnmlfile
character(fnlen)        :: dotproductfile

namelist /Eulerslist/ datafile, txtfile, angledataset, dotproductfile, EMEBSDnmlfile, raddeg

dotproductfile = 'undefined'
txtfile = 'undefined'
datafile = 'undefined'
EMEBSDnmlfile = 'undefined'
angledataset = 'original'
raddeg = 'deg'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=Eulerslist)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(dotproductfile).eq.'undefined') then
  call Message%printError('readNameList:',' dotproductfile name is undefined in '//nmlfile)
 end if

 if (trim(txtfile).eq.'undefined') then
  call Message%printError('readNameList:',' txtfile name is undefined in '//nmlfile)
 end if

 if (trim(datafile).eq.'undefined') then
  call Message%printError('readNameList:',' datafile name is undefined in '//nmlfile)
 end if
end if

self%nml%dotproductfile = dotproductfile
self%nml%txtfile = txtfile
self%nml%datafile = datafile
self%nml%EMEBSDnmlfile = EMEBSDnmlfile
self%nml%angledataset = trim(angledataset)
self%nml%raddeg = raddeg

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! pass the namelist for the getEulers_T Class to the calling program

IMPLICIT NONE

class(getEulers_T), INTENT(INOUT)          :: self
type(getEulersNameListType)                :: nml

nml = self%nml

end function getNameList_


!--------------------------------------------------------------------------
function get_angledataset_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_angledataset_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! get angledataset from the getEulers_T class

IMPLICIT NONE

class(getEulers_T), INTENT(INOUT)     :: self
character(8)                          :: out

out = self%nml%angledataset

end function get_angledataset_

!--------------------------------------------------------------------------
subroutine set_angledataset_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_angledataset_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! set angledataset in the getEulers_T class

IMPLICIT NONE

class(getEulers_T), INTENT(INOUT)     :: self
character(8), INTENT(IN)              :: inp

self%nml%angledataset = inp

end subroutine set_angledataset_

!--------------------------------------------------------------------------
function get_raddeg_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_raddeg_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! get raddeg from the getEulers_T class

IMPLICIT NONE

class(getEulers_T), INTENT(INOUT)     :: self
character(3)                          :: out

out = self%nml%raddeg

end function get_raddeg_

!--------------------------------------------------------------------------
subroutine set_raddeg_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_raddeg_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! set raddeg in the getEulers_T class

IMPLICIT NONE

class(getEulers_T), INTENT(INOUT)     :: self
character(3), INTENT(IN)              :: inp

self%nml%raddeg = inp

end subroutine set_raddeg_

!--------------------------------------------------------------------------
function get_txtfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_txtfile_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! get txtfile from the getEulers_T class

IMPLICIT NONE

class(getEulers_T), INTENT(INOUT)     :: self
character(fnlen)                      :: out

out = self%nml%txtfile

end function get_txtfile_

!--------------------------------------------------------------------------
subroutine set_txtfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_txtfile_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! set txtfile in the getEulers_T class

IMPLICIT NONE

class(getEulers_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%txtfile = inp

end subroutine set_txtfile_

!--------------------------------------------------------------------------
function get_datafile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_datafile_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! get datafile from the getEulers_T class

IMPLICIT NONE

class(getEulers_T), INTENT(INOUT)     :: self
character(fnlen)                      :: out

out = self%nml%datafile

end function get_datafile_

!--------------------------------------------------------------------------
subroutine set_datafile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_datafile_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! set datafile in the getEulers_T class

IMPLICIT NONE

class(getEulers_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%datafile = inp

end subroutine set_datafile_

!--------------------------------------------------------------------------
function get_EMEBSDnmlfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_EMEBSDnmlfile_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! get EMEBSDnmlfile from the getEulers_T class

IMPLICIT NONE

class(getEulers_T), INTENT(INOUT)     :: self
character(fnlen)                      :: out

out = self%nml%EMEBSDnmlfile

end function get_EMEBSDnmlfile_

!--------------------------------------------------------------------------
subroutine set_EMEBSDnmlfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_EMEBSDnmlfile_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! set EMEBSDnmlfile in the getEulers_T class

IMPLICIT NONE

class(getEulers_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%EMEBSDnmlfile = inp

end subroutine set_EMEBSDnmlfile_

!--------------------------------------------------------------------------
function get_dotproductfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dotproductfile_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! get dotproductfile from the getEulers_T class

IMPLICIT NONE

class(getEulers_T), INTENT(INOUT)     :: self
character(fnlen)                      :: out

out = self%nml%dotproductfile

end function get_dotproductfile_

!--------------------------------------------------------------------------
subroutine set_dotproductfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_dotproductfile_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! set dotproductfile in the getEulers_T class

IMPLICIT NONE

class(getEulers_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%dotproductfile = inp

end subroutine set_dotproductfile_

!--------------------------------------------------------------------------
subroutine getEulers_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: getEulers_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! perform the computations

use mod_EMsoft
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_io
use stringconstants
use mod_DIfiles
use mod_crystallography
use mod_symmetry

IMPLICIT NONE

class(getEulers_T), INTENT(INOUT)       :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname

type(HDF_T)                             :: HDF
type(HDFnames_T)                        :: HDFnames
type(IO_T)                              :: Message
type(DIfile_T)                          :: DIFT
type(cell_T)                            :: cell
type(SpaceGroup_T)                      :: SG
type(DictionaryIndexingNameListType)    :: dinl

logical                                 :: refined
character(fnlen)                        :: ename, nmlname, DIfile
integer(kind=irg)                       :: hdferr, j, istat, Nexp

real(kind=sgl),allocatable              :: euler_best(:,:)

associate(eunl=>self%nml, DIDT=>DIFT%DIDT)

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()

HDFnames = HDFnames_T()

call HDFnames%set_NMLfiles(SC_NMLfiles)
call HDFnames%set_NMLfilename(SC_DictionaryIndexingNML)
call HDFnames%set_NMLparameters(SC_NMLparameters)
call HDFnames%set_NMLlist(SC_DictionaryIndexingNameListType)

refined = .FALSE.
DIfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(eunl%dotproductfile)
if (trim(eunl%angledataset).eq.'refined') then
    call DIFT%readDotProductFile(EMsoft, HDF, HDFnames, DIfile, hdferr, &
                                getRefinedEulerAngles=.TRUE.)
    refined = .TRUE.
else
    call DIFT%readDotProductFile(EMsoft, HDF, HDFnames, DIfile, hdferr, &
                                getEulerAngles=.TRUE.)
end if

Nexp = DIDT%Nexp
allocate(euler_best(3,Nexp),stat=istat)
if (istat .ne. 0) then
    call Message%printError('EMgetEulers','Failed to allocate euler angle array')
end if
euler_best = 0.0
if (refined.eqv..TRUE.) then
  if (eunl%raddeg.eq.'deg') then
    euler_best(1:3,1:Nexp) = DIDT%RefinedEulerAngles(1:3,1:Nexp)*180.0/cPi
  else
    euler_best(1:3,1:Nexp) = DIDT%RefinedEulerAngles(1:3,1:Nexp)
  end if
  deallocate(DIDT%RefinedEulerAngles)
  call Message%printMessage(' Extracting refined Euler angles from dot product file')
else
  if (eunl%raddeg.eq.'deg') then
    euler_best(1:3,1:Nexp) = DIDT%EulerAngles(1:3,1:Nexp)*180.0/cPi
  else
    euler_best(1:3,1:Nexp) = DIDT%EulerAngles(1:3,1:Nexp)
  end if
  deallocate(DIDT%EulerAngles)
  call Message%printMessage(' Extracting Euler angles from dot product file')
end if

call Message%printMessage('  --> dot product EBSD HDF5 file read')

!==============================
! and prepare the Euler angle text file
ename = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(eunl%txtfile)

open(unit=dataunit,file=trim(ename),status='unknown',action='write')
write (dataunit,*) 'eu'
write (dataunit,*) Nexp

do j=1,Nexp
  write(dataunit,*) euler_best(1:3,j)
end do

close(unit=dataunit,status='keep')
call Message%printMessage('  --> closed Euler angle text file '//trim(ename))

!==============================
! if requested, also prepare an EMEBSD.nml file.
if (trim(eunl%EMEBSDnmlfile).ne.'undefined') then
  nmlname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(eunl%EMEBSDnmlfile)

! write all the relevant parameters to the EMEBSD.nml file
  open(UNIT=dataunit,FILE=trim(nmlname),form='formatted',STATUS='unknown')
  write (dataunit,"(A11)") '&EBSDdata'
  write (dataunit,"(' L = ',F12.4,',')") dinl%L
  write (dataunit,"(' thetac = ',F12.4,',')") dinl%thetac
  write (dataunit,"(' delta = ',F12.4,',')") dinl%delta
  write (dataunit,"(' numsx = ',I4,',')") dinl%numsx
  write (dataunit,"(' numsy = ',I4,',')") dinl%numsy
  write (dataunit,"(' xpc = ',F12.4,',')") dinl%xpc
  write (dataunit,"(' ypc = ',F12.4,',')") dinl%ypc
  write (dataunit,"(' omega = ',F12.4,',')") dinl%omega
  write (dataunit,"(' alphaBD = ',F12.4,',')") 0.0
  write (dataunit,"(' energymin = ',F12.4,',')") dinl%energymin
  write (dataunit,"(' energymax = ',F12.4,',')") dinl%energymax
  write (dataunit,"(' includebackground = ''',A1,''',')") 'y'
  write (dataunit,"(' anglefile = ''',A,''',')") trim(eunl%txtfile)
  write (dataunit,"(' anglefiletype = ''',A12,''',')") 'orientations'
  write (dataunit,"(' eulerconvention = ''',A3,''',')") 'tsl'
  write (dataunit,"(' makedictionary = ''',A1,''',')") 'n'
  write (dataunit,"(' masterfile = ''',A,''',')") trim(dinl%masterfile)
  write (dataunit,"(' datafile = ''',A,''',')") trim(eunl%datafile)
  write (dataunit,"(' bitdepth = ''',A4,''',')") '8bit'
  write (dataunit,"(' beamcurrent = ',F12.4,',')") dinl%beamcurrent
  write (dataunit,"(' dwelltime = ',F12.4,',')") dinl%dwelltime
  write (dataunit,"(' poisson = ''',A1,''',')") 'n'
  write (dataunit,"(' binning = ',I2,',')") dinl%binning
  write (dataunit,"(' applyDeformation = ''',A1,''',')") 'n'
  write (dataunit,"(' scalingmode = ''',A3,''',')") dinl%scalingmode
  write (dataunit,"(' gammavalue = ',F12.4,',')") dinl%gammavalue
  write (dataunit,"(' nthreads = ',I2,',')") dinl%nthreads
  write (dataunit,"(A2)") ' /'
  close(UNIT=dataunit,STATUS='keep')
  call Message%printMessage('')
  call Message%printMessage('      You may need to edit the namelist file to make sure that the')
  call Message%printMessage('      file paths are correct, or to change any other parameters before')
  call Message%printMessage('      running the EMEBSD program with this input file.')
  call Message%printMessage('')
  call Message%printMessage('  --> closed EMEBSD name list file '//trim(nmlname))
end if

! close the fortran HDF interface
call closeFortranHDFInterface()

end associate

end subroutine getEulers_

end module mod_getEulers
