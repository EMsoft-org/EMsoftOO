! ###################################################################
! Copyright (c) 2013-2022, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_FitOrientation
  !! author: MDG
  !! version: 1.0
  !! date: 04/08/20
  !!
  !! class definition for the EMFitOrientation program

use mod_kinds
use mod_global

IMPLICIT NONE

! namelist for the EMFitOrientation program
type, public :: FitOrientationNameListType
  integer(kind=irg) :: nthreads
  integer(kind=irg) :: matchdepth
  character(fnlen)  :: dotproductfile
  character(fnlen)  :: ctffile
  character(fnlen)  :: angfile
  character(fnlen)  :: tmpfile
  character(fnlen)  :: PSvariantfile
  character(fnlen)  :: method
  character(4)      :: modality
  logical           :: inRAM
  real(kind=sgl)    :: step
  integer(kind=irg) :: nmis
  integer(kind=irg) :: niter
  integer(kind=irg) :: initialx
  integer(kind=irg) :: initialy
  character(fnlen)  :: PCcorrection
  real(kind=sgl)    :: truedelta
end type FitOrientationNameListType

! class definition
type, public :: FitOrientation_T
private
  character(fnlen)                  :: nmldeffile = 'EMFitOrientation.nml'
  type(FitOrientationNameListType)  :: nml

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeFitHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: FitOrientation_
  procedure, pass(self) :: get_nthreads_
  procedure, pass(self) :: get_matchdepth_
  procedure, pass(self) :: get_dotproductfile_
  procedure, pass(self) :: get_ctffile_
  procedure, pass(self) :: get_angfile_
  procedure, pass(self) :: get_tmpfile_
  procedure, pass(self) :: get_PSvariantfile_
  procedure, pass(self) :: get_method_
  procedure, pass(self) :: get_modality_
  procedure, pass(self) :: get_inRAM_
  procedure, pass(self) :: get_step_
  procedure, pass(self) :: get_nmis_
  procedure, pass(self) :: get_niter_
  procedure, pass(self) :: set_nthreads_
  procedure, pass(self) :: set_matchdepth_
  procedure, pass(self) :: set_dotproductfile_
  procedure, pass(self) :: set_ctffile_
  procedure, pass(self) :: set_angfile_
  procedure, pass(self) :: set_tmpfile_
  procedure, pass(self) :: set_PSvariantfile_
  procedure, pass(self) :: set_method_
  procedure, pass(self) :: set_modality_
  procedure, pass(self) :: set_inRAM_
  procedure, pass(self) :: set_step_
  procedure, pass(self) :: set_nmis_
  procedure, pass(self) :: set_niter_

  generic, public :: getNameList => getNameList_
  ! generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: FitOrientation => FitOrientation_
  generic, public :: get_nthreads => get_nthreads_
  generic, public :: get_matchdepth => get_matchdepth_
  generic, public :: get_dotproductfile => get_dotproductfile_
  generic, public :: get_ctffile => get_ctffile_
  generic, public :: get_angfile => get_angfile_
  generic, public :: get_tmpfile => get_tmpfile_
  generic, public :: get_PSvariantfile => get_PSvariantfile_
  generic, public :: get_method => get_method_
  generic, public :: get_modality => get_modality_
  generic, public :: get_inRAM => get_inRAM_
  generic, public :: get_step => get_step_
  generic, public :: get_nmis => get_nmis_
  generic, public :: get_niter => get_niter_
  generic, public :: set_nthreads => set_nthreads_
  generic, public :: set_matchdepth => set_matchdepth_
  generic, public :: set_dotproductfile => set_dotproductfile_
  generic, public :: set_ctffile => set_ctffile_
  generic, public :: set_angfile => set_angfile_
  generic, public :: set_tmpfile => set_tmpfile_
  generic, public :: set_PSvariantfile => set_PSvariantfile_
  generic, public :: set_method => set_method_
  generic, public :: set_modality => set_modality_
  generic, public :: set_inRAM => set_inRAM_
  generic, public :: set_step => set_step_
  generic, public :: set_nmis => set_nmis_
  generic, public :: set_niter => set_niter_
end type FitOrientation_T

! the constructor routine for this class
interface FitOrientation_T
  module procedure FitOrientation_constructor
end interface FitOrientation_T

contains

!--------------------------------------------------------------------------
type(FitOrientation_T) function FitOrientation_constructor( nmlfile ) result(FitOrientation)
!DEC$ ATTRIBUTES DLLEXPORT :: FitOrientation_constructor
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! constructor for the FitOrientation_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call FitOrientation%readNameList(nmlfile)

end function FitOrientation_constructor

!--------------------------------------------------------------------------
subroutine FitOrientation_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: FitOrientation_destructor
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! destructor for the FitOrientation_T Class

IMPLICIT NONE

type(FitOrientation_T), INTENT(INOUT)  :: self

call reportDestructor('FitOrientation_T')

end subroutine FitOrientation_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! read the namelist from an nml file for the FitOrientation_T Class

use mod_io
use mod_EMsoft

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft
type(IO_T)                           :: Message
logical                              :: skipread = .FALSE.


integer(kind=irg)   :: nthreads
integer(kind=irg)   :: matchdepth
character(fnlen)    :: dotproductfile
character(fnlen)    :: ctffile
character(fnlen)    :: angfile
character(fnlen)    :: tmpfile
character(fnlen)    :: PSvariantfile
character(fnlen)    :: method
character(4)        :: modality
logical             :: inRAM
integer(kind=irg)   :: nmis
integer(kind=irg)   :: niter
real(kind=sgl)      :: step
integer(kind=irg)   :: initialx
integer(kind=irg)   :: initialy
character(fnlen)    :: PCcorrection
real(kind=sgl)      :: truedelta

namelist / RefineOrientations / nthreads, dotproductfile, ctffile, modality, nmis, niter, step, inRAM, method, &
                                matchdepth, PSvariantfile, tmpfile, angfile, initialx, initialy, PCcorrection, truedelta

nthreads = 1
matchdepth = 1
dotproductfile = 'undefined'
ctffile = 'undefined'
angfile = 'undefined'
tmpfile = 'undefined'
PSvariantfile = 'undefined'
method = 'FIT'
inRAM = .FALSE.
nmis = 1
niter = 1
step = 1.0
modality = 'EBSD'
initialx = 0
initialy = 0
PCcorrection = 'off'
truedelta = 50.0

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=RefineOrientations)
    close(UNIT=dataunit,STATUS='keep')

! check for required entries
    if (trim(dotproductfile).eq.'undefined') then
        call Message%printError('readNameList:',' dotproduct file name is undefined in '//nmlfile)
    end if

    if (trim(ctffile).eq.'undefined') then
        call Message%printError('readNameList:',' ctf file name is undefined in '//nmlfile)
    end if

    if (trim(tmpfile).eq.'undefined') then
        call Message%printError('readNameList:',' tmp file name is undefined in '//nmlfile)
    end if
end if

self%nml%nthreads = nthreads
self%nml%matchdepth = matchdepth
self%nml%dotproductfile = dotproductfile
self%nml%ctffile = ctffile
self%nml%angfile = angfile
self%nml%tmpfile = tmpfile
self%nml%PSvariantfile = PSvariantfile
self%nml%method = method
self%nml%inRAM = inRAM
self%nml%nmis = nmis
self%nml%niter = niter
self%nml%step = step
self%nml%modality = modality
self%nml%initialx = initialx 
self%nml%initialy = initialy
self%nml%PCcorrection = PCcorrection
self%nml%truedelta = truedelta 

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! pass the namelist for the FitOrientation_T Class to the calling program

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)          :: self
type(FitOrientationNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeFitHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeFitHDFNameList_
!! author: MDG
!! version: 1.0
!! date: 04/22/20
!!
!! write namelist to HDF file

use HDF5
use mod_HDFsupport
use mod_HDFnames
use stringconstants

use ISO_C_BINDING

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)  :: self
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 7, n_real = 1
integer(kind=irg)                       :: hdferr,  io_int(n_int), inRAMi
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)
logical                                 :: g_exists, overwrite = .TRUE.

associate( ronl => self%nml )

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! integers 
io_int = (/ ronl%nthreads, ronl%matchdepth, ronl%nmis, ronl%niter, 0 , ronl%initialx, ronl%initialy/)
if (ronl%inRAM.eqv..TRUE.) then 
  io_int(5) = 1
end if
intlist(1) = 'nthreads'
intlist(2) = 'matchdepth'
intlist(3) = 'nmis'
intlist(4) = 'niter'
intlist(5) = 'inRAM'
intlist(6) = 'initialx'
intlist(7) = 'initialy'

! and write them to the HDF file
call HDF%writeNMLintegers(io_int, intlist, n_int)

! floats
reallist = (/ 'step' /)
io_real(1) = ronl%step
call HDF%writeNMLreals(io_real, reallist, n_real)

! strings

dataset = 'PCcorrection'
line2(1) = ronl%PCcorrection
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then 
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create PCcorrection dataset',hdferr)

dataset = 'modality'
line2(1) = ronl%modality
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create modality dataset',hdferr)

dataset = 'dotproductfile'
line2(1) = ronl%dotproductfile
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create dotproductfile dataset',hdferr)

dataset = 'ctffile'
line2(1) = ronl%ctffile
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create ctffile dataset',hdferr)

dataset = 'angfile'
line2(1) = ronl%angfile
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create angfile dataset',hdferr)

dataset = 'tmpfile'
line2(1) = ronl%tmpfile
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create tmpfile dataset',hdferr)

dataset = 'PSvariantfile'
line2(1) = ronl%PSvariantfile
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create PSvariantfile dataset',hdferr)

dataset = 'method'
line2(1) = ronl%method
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
end if
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create method dataset',hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeFitHDFNameList_


!--------------------------------------------------------------------------
function get_nthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nthreads_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! get nthreads from the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
integer(kind=irg)                          :: out

out = self%nml%nthreads

end function get_nthreads_

!--------------------------------------------------------------------------
subroutine set_nthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nthreads_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! set nthreads in the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)              :: inp

self%nml%nthreads = inp

end subroutine set_nthreads_

!--------------------------------------------------------------------------
function get_matchdepth_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_matchdepth_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! get matchdepth from the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
integer(kind=irg)                          :: out

out = self%nml%matchdepth

end function get_matchdepth_

!--------------------------------------------------------------------------
subroutine set_matchdepth_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_matchdepth_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! set matchdepth in the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)              :: inp

self%nml%matchdepth = inp

end subroutine set_matchdepth_

!--------------------------------------------------------------------------
function get_dotproductfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dotproductfile_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! get dotproductfile from the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
character(fnlen)                           :: out

out = self%nml%dotproductfile

end function get_dotproductfile_

!--------------------------------------------------------------------------
subroutine set_dotproductfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_dotproductfile_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! set dotproductfile in the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)               :: inp

self%nml%dotproductfile = inp

end subroutine set_dotproductfile_

!--------------------------------------------------------------------------
function get_ctffile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ctffile_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! get ctffile from the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
character(fnlen)                           :: out

out = self%nml%ctffile

end function get_ctffile_

!--------------------------------------------------------------------------
subroutine set_ctffile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ctffile_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! set ctffile in the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)               :: inp

self%nml%ctffile = inp

end subroutine set_ctffile_

!--------------------------------------------------------------------------
function get_angfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_angfile_
!! author: MDG
!! version: 1.0
!! date: 04/20/20
!!
!! get angfile from the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
character(fnlen)                           :: out

out = self%nml%angfile

end function get_angfile_

!--------------------------------------------------------------------------
subroutine set_angfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_angfile_
!! author: MDG
!! version: 1.0
!! date: 04/20/20
!!
!! set angfile in the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)               :: inp

self%nml%angfile = inp

end subroutine set_angfile_

!--------------------------------------------------------------------------
function get_tmpfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_tmpfile_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! get tmpfile from the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
character(fnlen)                           :: out

out = self%nml%tmpfile

end function get_tmpfile_

!--------------------------------------------------------------------------
subroutine set_tmpfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_tmpfile_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! set tmpfile in the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)               :: inp

self%nml%tmpfile = inp

end subroutine set_tmpfile_

!--------------------------------------------------------------------------
function get_PSvariantfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_PSvariantfile_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! get PSvariantfile from the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
character(fnlen)                           :: out

out = self%nml%PSvariantfile

end function get_PSvariantfile_

!--------------------------------------------------------------------------
subroutine set_PSvariantfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_PSvariantfile_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! set PSvariantfile in the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)               :: inp

self%nml%PSvariantfile = inp

end subroutine set_PSvariantfile_

!--------------------------------------------------------------------------
function get_method_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_method_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! get method from the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
character(fnlen)                           :: out

out = self%nml%method

end function get_method_

!--------------------------------------------------------------------------
subroutine set_method_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_method_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! set method in the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)               :: inp

self%nml%method = inp

end subroutine set_method_

!--------------------------------------------------------------------------
function get_modality_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_modality_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! get modality from the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
character(4)                               :: out

out = self%nml%modality

end function get_modality_

!--------------------------------------------------------------------------
subroutine set_modality_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_modality_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! set modality in the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
character(4), INTENT(IN)                   :: inp

self%nml%modality = inp

end subroutine set_modality_

!--------------------------------------------------------------------------
function get_inRAM_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_inRAM_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! get inRAM from the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
logical                                    :: out

out = self%nml%inRAM

end function get_inRAM_

!--------------------------------------------------------------------------
subroutine set_inRAM_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_inRAM_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! set inRAM in the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
logical, INTENT(IN)                        :: inp

self%nml%inRAM = inp

end subroutine set_inRAM_

!--------------------------------------------------------------------------
function get_step_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_step_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! get step from the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
real(kind=sgl)                             :: out

out = self%nml%step

end function get_step_

!--------------------------------------------------------------------------
subroutine set_step_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_step_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! set step in the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)                 :: inp

self%nml%step = inp

end subroutine set_step_

!--------------------------------------------------------------------------
function get_nmis_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nmis_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! get nmis from the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
integer(kind=irg)                          :: out

out = self%nml%nmis

end function get_nmis_

!--------------------------------------------------------------------------
subroutine set_nmis_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nmis_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! set nmis in the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)              :: inp

self%nml%nmis = inp

end subroutine set_nmis_

!--------------------------------------------------------------------------
function get_niter_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_niter_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! get niter from the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
integer(kind=irg)                          :: out

out = self%nml%niter

end function get_niter_

!--------------------------------------------------------------------------
subroutine set_niter_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_niter_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! set niter in the FitOrientation_T class

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)              :: inp

self%nml%niter = inp

end subroutine set_niter_

!--------------------------------------------------------------------------
subroutine FitOrientation_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: FitOrientation_
!! author: MDG
!! version: 1.0
!! date: 04/08/20
!!
!! perform the computations

use mod_EMsoft
use mod_crystallography
use mod_symmetry
use mod_patterns
use mod_vendors
use HDF5
use h5im
use h5lt
use mod_HDFsupport
use mod_HDFnames
use mod_math
use mod_MCfiles
use mod_MPfiles
use mod_DIfiles
use mod_EBSD
use mod_ECP
use mod_timing
use mod_rotations
use mod_quaternions
use mod_so3
use omp_lib
use mod_filters
use mod_timing
use mod_memory
use mod_io
use omp_lib
use mod_OMPsupport
use mod_bobyqa_refinement,only:bobyqa
use mod_FitOrientations,only:EMFitOrientationcalfunEBSD
use stringconstants

IMPLICIT NONE

class(FitOrientation_T), INTENT(INOUT)  :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname

type(IO_T)                              :: Message
type(HDF_T)                             :: HDF
type(HDFnames_T)                        :: HDFnames
type(Cell_T)                            :: cell
type(SpaceGroup_T)                      :: SG
type(q_T)                               :: qu, q, myqu
type(Quaternion_T)                      :: qq, quat, quat2, qq2, qquat2, qqq, myquat
type(QuaternionArray_T)                 :: quPS
type(e_T)                               :: eu, myeu, euinp2
type(a_T)                               :: ax
type(r_T)                               :: rfz
type(h_T)                               :: ho, myho
type(c_T)                               :: cu
type(MCfile_T)                          :: MCFT
type(MPfile_T)                          :: MPFT
type(DIfile_T)                          :: DIFT
type(so3_T)                             :: SO
type(QuaternionArray_T)                 :: qdummy, qAR
type(Timing_T)                          :: timer
type(EBSD_T)                            :: EBSD, myEBSD
type(ECP_T)                             :: ECP
type(Vendor_T)                          :: VT
type(memory_T)                          :: mem, memth

! type(EBSDIndexingNameListType)          :: dinl
type(MCOpenCLNameListType)              :: mcnl
type(EBSDMasterNameListType)            :: mpnl
! type(EBSDNameListType)                  :: ebsdnl


logical                                 :: stat, readonly, noindex, ROIselected
character(fnlen)                        :: dpfile, masterfile, energyfile
integer(kind=irg)                       :: hdferr, ii, jj, kk, iii, istat, npy, jjj, iparecp(4)

real(kind=dbl)                          :: misang       ! desired misorientation angle (degrees)
integer(kind=irg)                       :: Nmis         ! desired number of sampling points along cube edge
integer(kind=irg)                       :: CMcnt        ! number of entries in linked list
type(FZpointd),pointer                  :: CMlist, CMtmp       ! pointer to start of linked list and temporary one
real(kind=dbl)                          :: rhozero(4), hipassw

real(kind=sgl),allocatable              :: euPS(:,:), euler_bestmatch(:,:,:), CIlist(:), CMarray(:,:,:)
integer(kind=irg),allocatable           :: indexmain(:,:)
real(kind=sgl),allocatable              :: resultmain(:,:), DPCX(:), DPCY(:), DPCL(:)
integer(HSIZE_T)                        :: dims(1),dims2D(2),dims3(3),offset3(3)

character(fnlen, KIND=c_char),allocatable,TARGET    :: stringarray(:)
character(fnlen)                        :: dataset, groupname
character(fnlen)                        :: ename, fname
character(2)                            :: anglemode
real(kind=dbl),parameter                :: nAmpere = 6.241D+18   ! Coulomb per second

integer(c_size_t),allocatable           :: IPAR2(:)
real(kind=dbl),allocatable              :: X(:), XL(:), XU(:)
real(kind=sgl),allocatable              :: INITMEANVAL(:)
real(kind=dbl)                          :: RHOBEG, RHOEND
integer(kind=irg)                       :: NPT, N, IPRINT, NSTEP, NINIT
integer(kind=irg),parameter             :: MAXFUN = 10000
logical                                 :: verbose

logical                                 :: f_exists, init, g_exists, overwrite, isEBSD=.FALSE., isTKD=.FALSE., &
                                           isECP=.FALSE., switchwfoff
integer(kind=irg),parameter             :: iunitexpt = 41, itmpexpt = 42
integer(kind=irg)                       :: binx, biny, recordsize, pos(2), nsig, numk, FZt, FZo 
real(kind=sgl),allocatable              :: tmpimageexpt(:), EBSDPattern(:,:), mask(:,:), masklin(:), imageexpt(:)
real(kind=sgl),allocatable              :: imagedictflt(:), exppatarray(:)
real(kind=sgl),allocatable              :: EBSDpatternintd(:,:), binned(:,:), euler_best(:,:)
real(kind=sgl),allocatable              :: epatterns(:,:)
real(kind=sgl),allocatable              :: anglewf(:)
integer(kind=irg),allocatable           :: kij(:,:)
integer(kind=irg),allocatable           :: EBSDpatterninteger(:,:), EBSDpatternad(:,:), EBSDpint(:,:), Rvar(:)
type(C_PTR)                             :: planf, HPplanf, HPplanb
integer(kind=8)                         :: size_in_bytes_dict,size_in_bytes_expt
type(IncidentListECP),pointer           :: ktmp
real(kind=dbl),allocatable              :: klist(:,:)

real(kind=dbl)                          :: w, Jres, nel, emult, MCsig
real(kind=dbl)                          :: stpsz, cu0(3)
real(kind=dbl),allocatable              :: cubneighbor(:,:)
real(kind=sgl)                          :: ma, mi, dp, tstart, tstop, io_real(4), tmp, &
                                           vlen, avec(3), fpar1(1), fpar2(2), WD, alpha, ca, sa, c2a, s2a, nn(3), omega, dx, dy, rho
integer(kind=irg)                       :: ipar(10), Emin, Emax, nthreads, TID, io_int(2), tickstart, ierr, L, nvar, niter,i,j, &
                                           samplex, sampley, maxeindex, unchanged
integer(kind=irg)                         :: ll, mm, jpar(7), Nexp, pgnum, FZcnt, nlines, dims2(2), correctsize, totnumexpt, mystat

real(kind=dbl)                          :: prefactor, F, angleaxis(4)
real(kind=sgl),allocatable              :: axPS(:,:), dpPS(:,:), eulerPS(:,:,:)
real(kind=dbl)                          :: ratioE, euinp(3)
integer(kind=irg)                       :: cratioE, fratioE, eindex, jjend, iiistart, iiiend, Nd, Ne, patsz, pp
integer(kind=irg),allocatable           :: ppendE(:)
real(kind=sgl),allocatable              :: exptpatterns(:,:)

real(kind=sgl),allocatable              :: STEPSIZE(:), OSMmap(:,:), IQmap(:)
character(fnlen)                        :: modalityname, DIfile, xtalname


call setRotationPrecision('d')
call OMP_showAvailableThreads()

associate(ronl=>self%nml, dinl=>DIFT%nml, DIDT=>DIFT%DIDT, MCDT=>MCFT%MCDT, &
          MPDT=>MPFT%MPDT, det=>EBSD%det, mydet=> myEBSD%det, enl=>EBSD%nml, ecpnl=>ECP%nml)

init = .TRUE.
overwrite = .TRUE.
verbose = .FALSE.

!====================================
! read the relevant fields from the dot product HDF5 file

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()

! first we need to get the DIModality from the dot product file; this then
! determines a number of other parameters, including the HDFnames as well
! as the various arrays that we can read from the DI file
DIfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(ronl%dotproductfile)
DIFT = DIfile_T()
call DIFT%readDIModality(HDF, DIfile)
modalityname = DIFT%getModality()
! maybe this is an old dot product file (pre-6.0) so we use the modality switch in the 
! namelist for this program to set the modality 
if (trim(modalityname).eq.'unknown') then    
  modalityname = trim(ronl%modality)
end if

HDFnames = HDFnames_T()

call HDFnames%set_NMLfiles(SC_NMLfiles)
call HDFnames%set_NMLfilename(SC_DictionaryIndexingNML)
call HDFnames%set_NMLparameters(SC_NMLparameters)
call HDFnames%set_NMLlist(SC_DictionaryIndexingNameListType)

! allocate memory handling classes
mem = memory_T()
! call mem%toggle_verbose()
memth = memory_T( nt = ronl%nthreads )
! call memth%toggle_verbose()

!====================================
! read the relevant fields from the dot product HDF5 file
!====================================
if ( (trim(modalityname) .eq. 'EBSD').or.(trim(modalityname) .eq. 'TKD') )  then
  if ( (ronl%matchdepth.eq.1).or.(trim(ronl%method).eq.'SUB') ) then 
    call DIFT%readDotProductFile(EMsoft, HDF, HDFnames, DIfile, hdferr, &
                                 getCI=.TRUE., &
                                 getIQ=.TRUE., &
                                 getOSM=.TRUE., &
                                 getPhi1=.TRUE., &
                                 getPhi=.TRUE., &
                                 getPhi2=.TRUE.)

    w = dinl%hipassw
    Nexp = DIDT%Nexp
    call mem%alloc(euler_bestmatch, (/ 3,1,Nexp /),'euler_bestmatch',0.0)
    call mem%alloc(euler_best, (/ 3,Nexp /), 'euler_best', 0.0)
    call mem%alloc(CIlist, (/ Nexp /), 'CIlist', 0.0)
    euler_bestmatch(1,1,1:Nexp) = DIDT%Phi1(1:Nexp)
    euler_bestmatch(2,1,1:Nexp) = DIDT%Phi(1:Nexp)
    euler_bestmatch(3,1,1:Nexp) = DIDT%Phi2(1:Nexp)
    deallocate(DIDT%Phi1,DIDT%Phi,DIDT%Phi2)
  else
    call DIFT%readDotProductFile(EMsoft, HDF, HDFnames, DIfile, hdferr, &
                                 getCI=.TRUE., &
                                 getIQ=.TRUE., &
                                 getOSM=.TRUE., &
                                 getEulerAngles=.TRUE., &
                                 getDictionaryEulerAngles=.TRUE., &
                                 getTopMatchIndices=.TRUE.)

    w = dinl%hipassw
    Nexp = DIDT%Nexp
    call mem%alloc(euler_bestmatch, (/ 3,ronl%matchdepth,Nexp /),'euler_bestmatch',0.0)
    call mem%alloc(euler_best, (/ 3,Nexp /), 'euler_best', 0.0)
    call mem%alloc(CIlist, (/ Nexp /), 'CIlist', 0.0)
! read the appropriate set of Euler angles from the array of near matches
    do ii=1,ronl%matchdepth
      do jj=1,Nexp
        euler_bestmatch(1:3,ii,jj) = DIDT%DictionaryEulerAngles(1:3,DIDT%TopMatchIndices(ii,jj))
      end do
    end do
    euler_bestmatch = euler_bestmatch * dtor
    deallocate(DIDT%DictionaryEulerAngles, DIDT%TopMatchIndices)
  end if

  CIlist(1:Nexp) = DIDT%CI(1:Nexp)

  deallocate(DIDT%CI)

! the following arrays are kept in the EBSDDIdata structure
!   OSMmap = EBSDDIdata%OSM
!   IQmap = EBSDDIdata%IQ

else  ! ECP modality (to be written)
   call Message%printError('FitOrientation','This program only handles EBSPs; use EMECPDIrefine for ECPs') 
end if

call Message%printMessage(' -> completed reading of dot product file')

!=========================

Ne = dinl%numexptsingle
Nd = dinl%numdictsingle
nthreads = ronl%nthreads
dinl%nthreads = ronl%nthreads

if (sum(dinl%ROI).ne.0) then
  ROIselected = .TRUE.
  iiistart = dinl%ROI(2)
  iiiend = dinl%ROI(2)+dinl%ROI(4)-1
  jjend = dinl%ROI(3)
else
  ROIselected = .FALSE.
  iiistart = 1
  iiiend = dinl%ipf_ht
  jjend = dinl%ipf_wd
end if

if (ROIselected.eqv..TRUE.) then
    totnumexpt = dinl%ROI(3)*dinl%ROI(4)
else
    totnumexpt = dinl%ipf_wd*dinl%ipf_ht
end if

!===================================================================================
!===============READ MASTER AND MC FILE=============================================
!===================================================================================
!
! determine the modality from the master pattern file, and also set it in the dinl name list
fname = EMsoft%generateFilePath('EMdatapathname',trim(dinl%masterfile))
call MPFT%determineModality(HDF, fname)
call Message%printMessage(' Master Pattern modality : '//trim(MPFT%getModality()))
call DIFT%setModality(MPFT%getModality())

if (trim(MPFT%getModality()).eq.'EBSD') then
  isEBSD = .TRUE.
else if (trim(MPFT%getModality()).eq.'TKD') then
  isTKD = .TRUE.
else if (trim(MPFT%getModality()).eq.'ECP') then
  isECP = .TRUE.
  end if

! 1. read the Monte Carlo data file
call HDFnames%set_ProgramData(SC_MCOpenCL)
call HDFnames%set_NMLlist(SC_MCCLNameList)
call HDFnames%set_NMLfilename(SC_MCOpenCLNML)
fname = EMsoft%generateFilePath('EMdatapathname',trim(dinl%masterfile))
call MCFT%setFileName(fname)
call MCFT%readMCfile(HDF, HDFnames, getAccume=.TRUE.)
mcnl = MCFT%getnml()
xtalname = trim(mcnl%xtalname)

! 2. read the master pattern file
if (isTKD.eqv..TRUE.) then
  call HDFnames%set_ProgramData(SC_TKDmaster)
  call HDFnames%set_NMLlist(SC_TKDmasterNameList)
  call HDFnames%set_NMLfilename(SC_TKDmasterNML)
end if
if (isEBSD.eqv..TRUE.) then
  call HDFnames%set_ProgramData(SC_EBSDmaster)
  call HDFnames%set_NMLlist(SC_EBSDmasterNameList)
  call HDFnames%set_NMLfilename(SC_EBSDmasterNML)
end if
if (isECP.eqv..TRUE.) then
  call HDFnames%set_ProgramData(SC_ECPmaster)
  call HDFnames%set_NMLlist(SC_ECPmasterNameList)
  call HDFnames%set_NMLfilename(SC_ECPmasterNML)
end if
call HDFnames%set_Variable(SC_MCOpenCL)

fname = EMsoft%generateFilePath('EMdatapathname',trim(dinl%masterfile))
call MPFT%setFileName(fname)
call MPFT%readMPfile(HDF, HDFnames, mpnl, getmLPNH=.TRUE., getmLPSH=.TRUE.)

! set the HDFnames for the current program (same for all modalities)
call HDFnames%set_ProgramData(SC_EMDI)
call HDFnames%set_NMLlist(SC_EMDINameList)
call HDFnames%set_NMLfilename(SC_EMDI)

! we know that the master pattern file exists, and it also has all the
! crystallographic data in it, so we read that here instead of assuming
! that the actual .xtal file exists on this system ...
call cell%setFileName(xtalname)
call cell%readDataHDF(SG, EMsoft, useXtalName=fname)
! extract the point group number
pgnum = SG%getPGnumber()
io_int = pgnum
call Message%WriteValue(' Setting point group number to ',io_int,1)

! 3. for EBSD/TKD copy a few parameters from dinl to enl
! and then generate the detector arrays
if ( (isEBSD.eqv..TRUE.) .or. (isTKD.eqv..TRUE.)) then
  call mem%alloc(det%rgx, (/ dinl%numsx,dinl%numsy /), 'det%rgx', 0.0)
  call mem%alloc(det%rgy, (/ dinl%numsx,dinl%numsy /), 'det%rgy', 0.0)
  call mem%alloc(det%rgz, (/ dinl%numsx,dinl%numsy /), 'det%rgz', 0.0)
  call mem%alloc(det%accum_e_detector, (/ MCDT%numEbins,dinl%numsx,dinl%numsy /), 'det%accum_e_detector')
  enl%numsx = dinl%numsx
  enl%numsy = dinl%numsy
  enl%xpc = dinl%xpc
  enl%ypc = dinl%ypc
  enl%delta = dinl%delta
  enl%thetac = dinl%thetac
  enl%L = dinl%L
  enl%energymin = dinl%energymin
  enl%energymax = dinl%energymax

  if (isTKD.eqv..TRUE.) then
    call EBSD%GenerateDetector(MCFT, verbose, isTKD)
  end if
  if (isEBSD.eqv..TRUE.) then
    call EBSD%GenerateDetector(MCFT, verbose)
  end if
else  ! this must be an ECP indexing run so we initialize the appropriate detector arrays
  if (isECP.eqv..TRUE.) then
    ECP = ECP_T()
  ! copy a few parameters
    ecpnl%conesemiangle = dinl%conesemiangle
    ecpnl%sampletilt = dinl%sampletilt
    ecpnl%npix = dinl%npix
    ecpnl%workingdistance = dinl%workingdistance
    ecpnl%Rin = dinl%Rin
    ecpnl%Rout = dinl%Rout

    call ECP%ECPGenerateDetector(verbose=.TRUE.)
    nsig = nint((ecpnl%conesemiangle) + abs(ecpnl%sampletilt)) + 1
    call mem%alloc(anglewf, (/ nsig /), 'anglewf', 0.0)

    call Message%printMessage(' -> Calculating weight factors', frm = "(A)" )
    call ECP%ECPGetWeightFactors(mcnl, MCFT, anglewf, nsig, verbose=.TRUE.)

    !=================================================================
    ! check if there are enough angles in MC for detector geometry
    !=================================================================
    if (mcnl%sigend .lt. (abs(ecpnl%sampletilt) + ecpnl%conesemiangle)) then
      call Message%printMessage('Not enough angles in Monte carlo file...interpolation will be done without &
      appropriate weight factors',frm = "(A)")
      switchwfoff = .TRUE.
    end if

    if ((-mcnl%sigend .gt. (ecpnl%conesemiangle - abs(ecpnl%sampletilt))) .and. (switchwfoff .eqv. .FALSE.)) then
      call Message%printMessage('Not enough angles in Monte carlo file...interpolation will be done without &
      appropriate weight factors',frm = "(A)")
      switchwfoff = .TRUE.
    end if

    !=================================================================
    ! generate list of incident vectors
    !=================================================================
    numk = 0
    call ECP%GetVectorsCone()
    numk = ECP%get_numk()
    call mem%alloc(kij, (/ 2,numk /), 'kij', 0)
    call mem%alloc(klist, (/ 3,numk /), 'klist', 0.D0)

    io_int(1) = numk
    call Message%WriteValue('Number of beams for which interpolation will be done = ',io_int,1)

    ktmp => ECP%get_ListHead()
    ! converting to array for OpenMP parallelization
    do i = 1,numk
       klist(1:3,i) = ktmp%k(1:3)
       kij(1:2,i) = (/ktmp%i,ktmp%j/)
       ktmp => ktmp%next
    end do
    iparecp(1) = nsig
    iparecp(2) = numk
    iparecp(3) = ecpnl%npix
    iparecp(4) = mpnl%npx
  end if
end if

! also copy the sample tilt angle into the correct variable for writing to the dot product file
MCsig = mcnl%sig

call Message%printMessage(' Completed reading all MC/MP input data; generated detector ')

! reset the HDFnames fields to the correct values for the present program
call HDFnames%set_NMLfiles(SC_NMLfiles)
call HDFnames%set_NMLfilename(SC_FitOrientationNML)
call HDFnames%set_NMLparameters(SC_NMLparameters)
call HDFnames%set_NMLlist(SC_FitOrientationNameListType)

!=====================================================
! get the indices of the minimum and maximum energy
!=====================================================
Emin = nint((dinl%energymin - mcnl%Ehistmin)/mcnl%Ebinsize) +1
if (Emin.lt.1)  Emin=1
if (Emin.gt.MCDT%numEbins)  Emin=MCDT%numEbins

Emax = nint((dinl%energymax - mcnl%Ehistmin)/mcnl%Ebinsize) + 1
if (Emax .lt. 1) Emax = 1
if (Emax .gt. MCDT%numEbins) Emax = MCDT%numEbins

!=====================================================
!==========fill important parameters in namelist======
!=====================================================

binx = dinl%numsx/dinl%binning
biny = dinl%numsy/dinl%binning
recordsize = binx*biny*4
L = binx*biny
npy = mpnl%npx

! make sure that correctsize is a multiple of 16; if not, make it so
if (mod(L,16) .ne. 0) then
    correctsize = 16*ceiling(float(L)/16.0)
else
    correctsize = L
end if

! determine the experimental and dictionary sizes in bytes
size_in_bytes_dict = Nd*correctsize*sizeof(correctsize)
size_in_bytes_expt = Ne*correctsize*sizeof(correctsize)
patsz              = correctsize

allocate(IPAR2(10))

! define the jpar array
jpar(1) = dinl%binning
jpar(2) = dinl%numsx
jpar(3) = dinl%numsy
jpar(4) = mpnl%npx
jpar(5) = npy
jpar(6) = MCDT%numEbins
jpar(7) = MCDT%numEbins
!jpar(7) = dinl%nE

IPAR2(1:7) = jpar(1:7)
IPAR2(8) = Emin
IPAR2(9) = Emax
IPAR2(10)= dinl%nregions

dims2 = (/binx, biny/)

! intensity prefactor
! modified by MDG, 03/26/18
nel = float(mcnl%totnum_el) * float(mcnl%multiplier)
emult = nAmpere * 1e-9 / nel  ! multiplicative factor to convert MC data to an equivalent incident beam of 1 nanoCoulomb
! intensity prefactor  (redefined by MDG, 3/23/18)
prefactor = emult * dinl%beamcurrent * dinl%dwelltime * 1.0D-6

!=====================================================
!=====================================================
! set up the equivalent pseudo-symmetric orientation quaternions
!=====================================================
!=====================================================
! do we have a PS variant file ?
if (trim(ronl%PSvariantfile).ne.'undefined') then
    ! we need to get the direct structure matrix to convert the axis-angle pair
    ! to the correct reference frame, so we'll read the complete xtal file
    ! and initialize all matrices as usual
    ! allocate (cell)
    call cell%getCrystalData(xtalname, SG, EMsoft)

    call Message%printMessage('Reading pseudo-symmetry variant operators: ')
    dpfile = trim(EMsoft%getConfigParameter('EMdatapathname'))//trim(ronl%PSvariantfile)

    ! this is a simple text file, similar to an euler angle file; the input should
    ! be in quaternion format, so abort when the file does not have quaternions...
    open(unit=53,file=trim(dpfile),status='old',action='read')
    read (53,*) anglemode
    if ((anglemode.ne.'ax').and.(anglemode.ne.'eu')) call Message%printError('EMFitOrientationPS','angle type must be eu or ax')
    read (53,*) nvar
    nvar = nvar + 1     ! identity operation is first entry

    quPS = QuaternionArray_T(nvar, s='d')

    if (anglemode.eq.'ax') then
    ! allocate some arrays
        call mem%alloc(axPS, (/ 4,nvar /), 'axPS', 0.0)
        axPS(1:4,1) = (/ 0.0, 0.0, 1.0, 0.0 /)

        do ii = 2,nvar
            read(53,"(4F12.9)") axPS(1:4,ii)
        end do
    ! the axis should be given in crystal coordinates as a direction, so
    ! we need to transform the axis first to the crystal cartesian frame
    ! using the direct structure matrix, and then normalize it...
        call Message%printMessage(' -> Converting operators to cartesian reference frame...')
        do ii=1,nvar
            call cell%TransSpace(axPS(1:3,ii),avec,'d','c')
            avec = avec/vecnorm(avec)

            axPS(1:3,ii) = avec(1:3)
            axPS(4,ii) = axPS(4,ii) * dtor
        end do

        call Message%printMessage(' -> Final pseudo-symmetric quaternion operator(s): ')
        do ii = 1,nvar
          ax = a_T( adinp = dble(axPS(1:4,ii)) )
          qu = ax%aq()
          qq = Quaternion_T( qd = qu%q_copyd() )
          call qq%quat_pos()
          call quPS%insertQuatinArray(ii, qq)
          io_real(1:4) = qq%get_quatd()
          call Message%WriteValue('',io_real,4)
        end do
        call mem%dealloc(axPS, 'axPS')
    else
        ! allocate some arrays
        call mem%alloc(euPS, (/ 3,nvar /), 'euPS', 0.0)
        euPS(1:3,1) = (/ 0.0, 0.0, 0.0 /)

        do ii = 2,nvar
            read(53,"(3F12.9)") euPS(1:3,ii)
        end do

        call Message%printMessage(' -> Final pseudo-symmetric quaternion operator(s): ')
        do ii = 1,nvar
          eu = e_T( edinp = dble(euPS(1:3,ii)) )
          qu = eu%eq()
          qq = Quaternion_T( qd = qu%q_copyd() )
          call qq%quat_pos()
          call quPS%insertQuatinArray(ii, qq)
          io_real(1:4) = qq%get_quatd()
          call Message%WriteValue('',io_real,4)
        end do
        call mem%dealloc(euPS, 'euPS')
    end if
    close(52,status='keep')
    call Message%printMessage('--------')
else  ! there are no pseudo-symmetric variants in this run
  nvar = 1     ! identity operation is the only entry
  quPS = QuaternionArray_T(nvar, s='d')
  call quPS%insertQuatinArray(1, Quaternion_T( qd = (/ 1.D0, 0.D0, 0.D0, 0.D0 /) ))
end if

!=====================================================
!=====================================================
! set up the correct fundamental zone and symmetry operators
SO = so3_T(pgnum, zerolist='FZ')
call SO%getFZtypeandorder(FZt, FZo)
io_int(1:2) = (/ FZt, FZo /)
call Message%WriteValue('Fundamental Zone parameters FZt, FZo : ',io_int, 2)
call qdummy%QSym_Init( pgnum, qAR )
! print out the symmetry array qAR
call qAR%quat_print()

!=====================================================
!==========ALLOCATE ALL ARRAYS HERE===================
!=====================================================
call mem%alloc(mask, (/ binx,biny /), 'mask', 1.0)
call mem%alloc(masklin, (/ binx*biny /), 'masklin', 0.0)

!===============================================================
! define the circular mask if necessary and convert to 1D vector
!===============================================================
if (dinl%maskpattern.eq.'y') then
  do ii = 1,biny
      do jj = 1,binx
          if((ii-biny/2)**2 + (jj-binx/2)**2 .ge. dinl%maskradius**2) then
              mask(jj,ii) = 0.0
          end if
      end do
  end do
end if

do ii = 1,biny
    do jj = 1,binx
        masklin((ii-1)*binx+jj) = mask(jj,ii)
    end do
end do

!===============================================================
!======== pre-process the experimental patterns=================
!===============================================================
! is the output to a temporary file or will it be kept in memory?
if (ronl%inRAM.eqv..TRUE.) then
! allocate the array that will hold all the processed experimental patterns
  call mem%alloc(epatterns, (/ correctsize,totnumexpt /), 'epatterns', 0.0)
  call PreProcessPatterns(EMsoft, HDF, .TRUE., dinl, binx, biny, masklin, correctsize, totnumexpt, epatterns=epatterns)
  io_real(1) = minval(epatterns)
  io_real(2) = maxval(epatterns)
  call Message%WriteValue(' --> preprocessed patterns intensity range (kept in RAM) = ',io_real,2)
else
  ! get the tmp file name from the input name list instead of the dot product file
  ! to allow for multiple instantiations of this program to run simultaneously
  dinl%tmpfile = trim(ronl%tmpfile)
  call PreProcessPatterns(EMsoft, HDF, .FALSE., dinl, binx, biny, masklin, correctsize, totnumexpt)
end if

!===============================================================
!========Pattern center correction parameters===================
!===============================================================

if (trim(ronl%PCcorrection).eq.'on') then 
  alpha = 0.5 * sngl(cPi) - (mcnl%sig - dinl%thetac) * dtor  
  ca = cos(alpha)
  c2a = cos(2.0*alpha)
  sa = sin(alpha)
  s2a = sin(2.0*alpha)
! determine the shift vector for each sampling point (on the sample!) with respect to the 
! (initialx, initialy) position
  if (ROIselected.eqv..TRUE.) then 
    call mem%alloc(DPCX, (/ dinl%ROI(3) /), 'DPCX', 0.0)
    call mem%alloc(DPCY, (/ dinl%ROI(4) /), 'DPCY', 0.0)
    call mem%alloc(DPCL, (/ dinl%ROI(4) /), 'DPCL', 0.0)
    do i=1,dinl%ROI(3)
      DPCX(i) = - ( ronl%initialx - (dinl%ROI(1)+(i-1)) ) * dinl%StepX
    end do 
    do j=1,dinl%ROI(4)
      DPCY(j) = - ( ronl%initialy - (dinl%ROI(2)+(j-1)) ) * dinl%StepY
    end do 
  else
    call mem%alloc(DPCX, (/ dinl%ipf_wd /), 'DPCX', 0.0)
    call mem%alloc(DPCY, (/ dinl%ipf_ht /), 'DPCY', 0.0)
    call mem%alloc(DPCL, (/ dinl%ipf_ht /), 'DPCL', 0.0)
    do i=1,dinl%ipf_wd
      DPCX(i) = - ( ronl%initialx - i ) * dinl%StepX
    end do 
    do j=1,dinl%ipf_ht
      DPCY(j) = - ( ronl%initialy - j ) * dinl%StepY
    end do 
  end if
! convert these shifts to shifts in the detector reference frame 
! and put them in units of the detector pixel size 
  DPCX = - DPCX / dinl%delta
  DPCL = - DPCY * sa 
  DPCY = - DPCY * ca / dinl%delta
end if  

!=================================
!========LOOP VARIABLES===========
!=================================

ratioE = float(Nexp)/float(dinl%numexptsingle)
cratioE = ceiling(ratioE)
fratioE = floor(ratioE)

ppendE = (/ (dinl%numexptsingle, ii=1,cratioE) /)
if (fratioE.lt.cratioE) then
  ppendE(cratioE) = MODULO(Nexp,dinl%numexptsingle)
end if

!===================================================================================
! method = 'SUB' ... define necessary parameters
!===================================================================================
if (ronl%method.eq.'SUB') then
    Nmis = ronl%nmis
    niter = ronl%niter
    call mem%alloc(cubneighbor, (/ 3, (2*Nmis + 1)**3 /), 'cubneighbor', 0.D0)
end if

!===================================================================================
!===============BOBYQA VARIABLES====================================================
!===================================================================================
N = 3

RHOBEG = 0.1D0
RHOEND = 0.0001D0
IPRINT = 0
NPT = N + 6

verbose = .FALSE.

if (ronl%inRAM.eqv..FALSE.) then
   fname = trim(EMsoft%getConfigParameter('EMtmppathname'))//trim(dinl%tmpfile)
   open(unit=itmpexpt,file=trim(fname),&
   status='unknown',form='unformatted',access='direct',recl=correctsize*4,iostat=ierr)
end if

!===================================================================================
!===============MAIN COMPUTATION LOOP===============================================
!===================================================================================
call OMP_setNThreads(ronl%nthreads)

if (ROIselected.eqv..TRUE.) then 
  maxeindex = dinl%ROI(3) * dinl%ROI(4)
else 
  maxeindex = dinl%ipf_wd * dinl%ipf_ht
end if 

timer = Timing_T()
call timer%Time_tick()

call mem%alloc(exptpatterns, (/ binx*biny, dinl%numexptsingle /), 'exptpatterns', 0.0)

unchanged = 0 

! call mem%allocated_memory_use()

! parameters for orientation correction
! if (trim(ronl%PCcorrection).eq.'on') then 
!   alpha = cPi/2.D0 - (ebsdnl%MCsig - ebsdnl%thetac)*cPi/180.0
!   ca = cos(alpha)
!   sa = sin(alpha)
!   drd = ebsdnl%stepY*ca
!   dtd = ebsdnl%stepX
!   zs = ebsdnl%L/ebsdnl%numsx/ebsdnl%delta
!   xs = 0.D0
!   ys = 0.D0
!   rho = sqrt(xs**2 + ys**2 + zs**2)
!   !r = (/(ys*ca + zs*sa)/rho, -xs/rho, (-ys*sa + zs*ca)/rho /)
!   r = (/sa, 0.D0, ca/)
! end if 

! depending on the ronl%method, we perform the optimization with different routines...
if (ronl%method.eq.'FIT') then

    call Message%printMessage(' --> Starting regular refinement loop')

    do iii = 1,cratioE
        call mem%alloc(tmpimageexpt, (/ binx*biny /), 'tmpimageexpt', 0.0)
        if (ronl%inRAM.eqv..FALSE.) then
            do jj = 1,ppendE(iii)
                eindex = (iii - 1)*Ne + jj
                read(itmpexpt,rec=eindex) tmpimageexpt
                exptpatterns(1:binx*biny,jj) = tmpimageexpt(1:binx*biny)
            end do
        end if
        call mem%dealloc(tmpimageexpt, 'tmpimageexpt')

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID,ii,tmpimageexpt,jj,quat,quat2,binned,ma,mi,eindex) &
!$OMP& PRIVATE(EBSDpatternintd,EBSDpatterninteger,EBSDpatternad,imagedictflt,kk,ll,mm, myEBSD) &
!$OMP& PRIVATE(X,INITMEANVAL,XL,XU,STEPSIZE,dpPS,eulerPS,rfz,euinp,pos, q, qu, qq2, qq, eu, ho, mystat)

          TID = OMP_GET_THREAD_NUM()

! allocate all private arrays; one thread at a time to prevent conflicts accessing the memth class
!$OMP CRITICAL
          call memth%alloc(X, (/ N /), 'X', 0.D0, TID=TID)
          call memth%alloc(XL, (/ N /), 'XL', 0.D0, TID=TID)
          call memth%alloc(XU, (/ N /), 'XU', 1.D0, TID=TID)
          call memth%alloc(INITMEANVAL, (/ N /), 'INITMEANVAL', 0.0, TID=TID)
          call memth%alloc(STEPSIZE, (/ N /), 'STEPSIZE', ronl%step, TID=TID)

          call memth%alloc(dpPS, (/ ronl%matchdepth, nvar /), 'dpPS', 0.0, TID=TID)
          call memth%alloc(eulerPS, (/ 3, ronl%matchdepth, nvar /), 'eulerPS', 0.0, TID=TID)

          call memth%alloc(tmpimageexpt, (/ binx*biny /), 'tmpimageexpt', 0.0, TID=TID)
          call memth%alloc(binned, (/ binx,biny /), 'binned', 0.0, TID=TID)
          
          call memth%alloc(EBSDpatternintd, (/ binx,biny /), 'EBSDpatternintd', 0.0, TID=TID)
          call memth%alloc(EBSDpatterninteger, (/ binx,biny /), 'EBSDpatterninteger', 0, TID=TID)
          call memth%alloc(EBSDpatternad, (/ binx,biny /), 'EBSDpatternad', 0, TID=TID)
          call memth%alloc(imagedictflt, (/ binx*biny /), 'imagedictflt', 0.0, TID=TID)

          if (trim(ronl%PCcorrection).eq.'on') then 
! allocate the necessary arrays 
            call memth%alloc(myEBSD%det%rgx, (/ dinl%numsx, dinl%numsy /), 'myEBSD%det%rgx', 0.0, TID=TID)
            call memth%alloc(myEBSD%det%rgy, (/ dinl%numsx, dinl%numsy /), 'myEBSD%det%rgy', 0.0, TID=TID)
            call memth%alloc(myEBSD%det%rgz, (/ dinl%numsx, dinl%numsy /), 'myEBSD%det%rgz', 0.0, TID=TID)
            call memth%alloc(myEBSD%det%accum_e_detector, (/ MCDT%numEbins, dinl%numsx, dinl%numsy /), &
                             'mydet%accum_e_detector', 0.0, TID=TID)
            myEBSD%det%accum_e_detector = EBSD%det%accum_e_detector
          end if 

! call memth%thread_memory_use()
!$OMP END CRITICAL 
!$OMP BARRIER

!$OMP DO SCHEDULE(DYNAMIC)
          do jj = 1,ppendE(iii)

            eindex = (iii - 1)*Ne + jj
            if (eindex.gt.maxeindex) CYCLE
            if (self%nml%inRAM.eqv..TRUE.) then
                tmpimageexpt(1:correctsize) = epatterns(1:correctsize,eindex)
                tmpimageexpt = tmpimageexpt/vecnorm(tmpimageexpt)
            else
                tmpimageexpt = exptpatterns(:,jj)
            end if

    ! calculate the dot product for each of the orientations in the neighborhood of the best match
    ! including the pseudosymmetric variant; do this for all the selected top matches (matchdepth)
            do kk = 1, self%nml%matchdepth

                eu = e_T( edinp = dble(euler_bestmatch(1:3,kk,eindex)) )
                qu = eu%eq()
                quat = Quaternion_T( qd = qu%q_copyd() )
                call quat%quat_normalize()

                do ll = 1,nvar
                    qq2 = quPS%getQuatfromArray(ll)
                    quat2 = qq2 * quat
                    call quat2%quat_normalize()
                    q = q_T( qdinp = quat2%get_quatd())  ! RFZ reduction requires q_T class
                    call SO%ReduceOrientationtoRFZ(q, qAR, rfz)

                    ho = rfz%rh()

                    INITMEANVAL(1:3) = ho%h_copyd()
                    X(1:3) = 0.5D0
                    call bobyqa (IPAR2, INITMEANVAL, tmpimageexpt, N, NPT, X, XL, &
                                 XU, RHOBEG, RHOEND, IPRINT, MAXFUN, EMFitOrientationcalfunEBSD, EBSD%det%accum_e_detector,&
                                 MPDT%mLPNH, MPDT%mLPSH, mask, prefactor, EBSD%det%rgx, EBSD%det%rgy, &
                                 EBSD%det%rgz, STEPSIZE, DIFT%nml%gammavalue, verbose)
                    ho = h_T( hdinp = dble(X*2.0*STEPSIZE - STEPSIZE + INITMEANVAL) )
                    eu = ho%he()
                    eulerPS(1:3,kk,ll) = eu%e_copyd()

                    call EMFitOrientationcalfunEBSD(IPAR2, INITMEANVAL, tmpimageexpt, EBSD%det%accum_e_detector, &
                                                    MPDT%mLPNH, MPDT%mLPSH, N, X, F, mask, prefactor, EBSD%det%rgx, &
                                                    EBSD%det%rgy, EBSD%det%rgz, STEPSIZE, DIFT%nml%gammavalue, verbose)

                    dpPS(kk,ll) = 1.D0 - F

! do we need to perform a pattern center correction ?  This would be necessary for large area
! scans.  First, we apply the equivalent rotation to the refined orientation, then we create a 
! new set of detector arrays for this pattern center location, and we do another refinement step
! to get the final corrected orientation.  At the end, we make sure the new orientation falls in 
! the appropriate RFZ.
                    if ( (trim(ronl%PCcorrection).eq.'on') .and. (eindex.le.maxeindex) ) then  
                      ! get the corrected pattern center coordinates
                      ! first undo the pattern center shift by an equivalent rotation (see J. Appl. Cryst. (2017). 50, 1664–1676)
                      
                      ! generate the new detector arrays 

                      if (ROIselected.eqv..TRUE.) then 
                        samplex = mod(eindex-1, dinl%ROI(3))+1
                        sampley = (eindex-1)/dinl%ROI(3)+1
                      else 
                        samplex = mod(eindex-1, dinl%ipf_wd)+1
                        sampley = (eindex-1)/dinl%ipf_wd+1
                      end if 
                      myEBSD%nml = enl 
                      dx = DPCX(samplex)
                      dy = DPCY(sampley)
                      myEBSD%nml%xpc = enl%xpc - dx
                      myEBSD%nml%ypc = enl%ypc - dy
                      myEBSD%nml%L = enl%L - DPCL(sampley)
                      call EBSD%GeneratemyEBSDDetector(MCFT, dinl%numsx, dinl%numsy, MCDT%numEbins, myEBSD%det%rgx, &
                      myEBSD%det%rgy, myEBSD%det%rgz, myEBSD%det%accum_e_detector, &
                      (/ myEBSD%nml%xpc, myEBSD%nml%ypc, myEBSD%nml%L /))
                      
                      ! first undo the pattern center shift by an equivalent rotation (see J. Appl. Cryst. (2017). 50, 1664–1676, eq.15)
                      if ((dx.ne.0.0).or.(dy.ne.0.0)) then 
                        myeu = e_T(edinp = dble(eulerPS(1:3,kk,ll)))
                        myqu = myeu%eq()
                        rho = dx**2+dy**2
                        nn = -(/dx*ca,-dy,-dx*sa/)/sqrt(rho)
                        omega = acos(myEBSD%nml%L/sqrt(myEBSD%nml%L**2 + dinl%delta**2 * rho))
                        qqq = Quaternion_T( qd = dble((/ cos(omega*0.5), sin(omega*0.5) * nn, 0.0, 0.0 /)) ) 
                        myquat = Quaternion_T(qd = myqu%q_copyd())

                        qquat2 = myquat * qqq

                        myqu = q_T( qdinp = qquat2%get_quatd())
                        myho = myqu%qh()
                        INITMEANVAL(1:3) = sngl(myho%h_copyd())
                      else
                        myeu = e_T(edinp = dble(eulerPS(1:3,kk,ll)))
                        myho = myeu%eh()
                        INITMEANVAL(1:3) = sngl(myho%h_copyd())
                      end if 

                      ! refine the orientation using the new detector array and initial orientation 
                      ! INITMEANVAL(1:3) = ho%h_copyd()
                      X = 0.5D0
                      call bobyqa (IPAR2, INITMEANVAL, tmpimageexpt, N, NPT, X, XL, &
                                XU, RHOBEG, RHOEND, IPRINT, MAXFUN, EMFitOrientationcalfunEBSD, &
                                myEBSD%det%accum_e_detector,MPDT%mLPNH, MPDT%mLPSH, mask, prefactor, &
                                myEBSD%det%rgx, myEBSD%det%rgy, myEBSD%det%rgz, STEPSIZE, DIFT%nml%gammavalue, &
                                verbose)
                      
                      ho = h_T( hdinp = dble(X*2.0*STEPSIZE - STEPSIZE + INITMEANVAL) )
                      eu = ho%he()
                      eulerPS(1:3,kk,ll) = eu%e_copyd()

                      call EMFitOrientationcalfunEBSD(IPAR2, INITMEANVAL, tmpimageexpt, &
                           myEBSD%det%accum_e_detector, MPDT%mLPNH, MPDT%mLPSH, N, X, F, mask, &
                           prefactor, myEBSD%det%rgx, myEBSD%det%rgy, myEBSD%det%rgz, STEPSIZE, &
                           DIFT%nml%gammavalue, verbose)

                      dpPS(kk,ll) = 1.D0 - F
              
                      ! and return this orientation to the RFZ
                      euinp(1:3) = eulerPS(1:3,kk,ll)
                      euinp2 = e_T( edinp = euinp(1:3) )
                      q = euinp2%eq() ! RFZ reduction requires q_T class
                      call SO%ReduceOrientationtoRFZ(q, qAR, rfz)
                      eu = rfz%re()
                      eulerPS(1:3,kk,ll) = sngl(eu%e_copyd())
                    end if 
                      
                end do
            end do

! updating the confidence index only if a better match is found
            dp = maxval(dpPS)
            if (dp .gt. CIlist(eindex)) then
              CIlist(eindex) = dp
              pos = maxloc(dpPS)
              euler_best(1:3,eindex) = eulerPS(1:3,pos(1),pos(2))
            else
              euler_best(1:3,eindex) = euler_bestmatch(1:3,1,eindex)
              !$OMP CRITICAL
              unchanged = unchanged + 1 
              !$OMP END CRITICAL
            end if

            if (mod(eindex,250) .eq. 0) then
                io_int(1) = eindex
                io_int(2) = totnumexpt
                call Message%Writevalue('      completed refining pattern #',io_int,2,'(I8,'' of '',I8)')
            end if

        end do
    !$OMP END DO

!$OMP CRITICAL
        call memth%dealloc(X, 'X', TID=TID)
        call memth%dealloc(XL, 'XL', TID=TID)
        call memth%dealloc(XU, 'XU', TID=TID)
        call memth%dealloc(INITMEANVAL, 'INITMEANVAL', TID=TID)
        call memth%dealloc(STEPSIZE, 'STEPSIZE', TID=TID)
        call memth%dealloc(dpPS, 'dpPS', TID=TID)
        call memth%dealloc(eulerPS, 'eulerPS', TID=TID)
        call memth%dealloc(tmpimageexpt, 'tmpimageexpt', TID=TID)
        call memth%dealloc(binned, 'binned', TID=TID)
        call memth%dealloc(EBSDpatternintd, 'EBSDpatternintd', TID=TID)
        call memth%dealloc(EBSDpatterninteger, 'EBSDpatterninteger', TID=TID)
        call memth%dealloc(EBSDpatternad, 'EBSDpatternad', TID=TID)
        call memth%dealloc(imagedictflt, 'imagedictflt', TID=TID)
        if (trim(ronl%PCcorrection).eq.'on') then
          call memth%dealloc(myEBSD%det%rgx, 'myEBSD%det%rgx', TID=TID)
          call memth%dealloc(myEBSD%det%rgy, 'myEBSD%det%rgy', TID=TID)
          call memth%dealloc(myEBSD%det%rgz, 'myEBSD%det%rgz', TID=TID)
          call memth%dealloc(myEBSD%det%accum_e_detector, 'myEBSD%det%accum_e_detector', TID=TID)
        end if
!$OMP END CRITICAL

    !$OMP BARRIER    
    !$OMP END PARALLEL

    end do
else  ! sub-divide the cubochoric grid in half steps and determine for which gridpoint the dot product is largest
    do iii = 1,cratioE

        exptpatterns = 0.0
        stpsz = LPs%ap/2.D0/DIFT%nml%ncubochoric/2.D0

        if (self%nml%inRAM.eqv..FALSE.) then
            call mem%alloc(tmpimageexpt, (/ binx*biny /), 'tmpimageexpt', 0.0)
            do jj = 1,ppendE(iii)
                eindex = (iii - 1)*DIFT%nml%numexptsingle + jj
                read(itmpexpt,rec=eindex) tmpimageexpt
                exptpatterns(1:binx*biny,jj) = tmpimageexpt(1:binx*biny)
            end do
            call mem%dealloc(tmpimageexpt, 'tmpimageexpt')
        end if

        do kk = 1,niter

            io_int(1) = kk
            io_int(2) = niter
            call Message%Writevalue(' --> Starting cubochoric grid refinement ',io_int,2,'(I3,'' of '',I3)')

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,tmpimageexpt,jj,quat,binned,ma,mi,eindex) &
!$OMP& PRIVATE(EBSDpatternintd,EBSDpatterninteger,EBSDpatternad,imagedictflt,ll,mm,dp) &
!$OMP& PRIVATE(cubneighbor,cu0, cu, eu)

!$OMP CRITICAL
          call memth%alloc(tmpimageexpt, (/ binx*biny /), 'tmpimageexpt', 0.0, TID=TID)
          call memth%alloc(binned, (/ binx,biny /), 'binned', 0.0, TID=TID)
          call memth%alloc(EBSDpatternintd, (/ binx,biny /), 'EBSDpatternintd', 0.0, TID=TID)
          call memth%alloc(EBSDpatterninteger, (/ binx,biny /), 'EBSDpatterninteger', 0, TID=TID)
          call memth%alloc(EBSDpatternad, (/ binx,biny /), 'EBSDpatternad', 0, TID=TID)
          call memth%alloc(imagedictflt, (/ binx*biny /), 'imagedictflt', 0.0, TID=TID)
!$OMP END CRITICAL

!$OMP DO SCHEDULE(DYNAMIC)      
            do ii = 1,ppendE(iii)

               eindex = (iii - 1)*DIFT%nml%numexptsingle + ii
               if (self%nml%inRAM.eqv..TRUE.) then
                    tmpimageexpt(1:correctsize) = epatterns(1:correctsize,eindex)
                    tmpimageexpt = tmpimageexpt/vecnorm(tmpimageexpt)
                else
                    tmpimageexpt = exptpatterns(:, ii)
                end if
                eu = e_T( edinp = dble(euler_bestmatch(1:3,1,eindex)) )
                cu = eu%ec()
                cu0 = cu%c_copyd()

                call SO%CubochoricNeighbors(cubneighbor,Nmis,cu0,stpsz)

    ! calculate the dot product for each of the orientations in the neighborhood of the best match
                do jj = 1,(2*Nmis + 1)**3
                    cu = c_T( cdinp = dble(cubneighbor(1:3,jj)) )
                    qu = cu%cq()
                    quat = Quaternion_T( qd = qu%q_copyd() )

                    call EBSD%CalcEBSDPatternSingleFull(jpar,quat,det%accum_e_detector,MPDT%mLPNH,MPDT%mLPSH,&
                                                        det%rgx, det%rgy,det%rgz,binned,Emin,Emax,mask,&
                                                        prefactor)

                    ma = maxval(binned)
                    mi = minval(binned)

                    EBSDpatternintd = ((binned - mi)/ (ma-mi))
                    EBSDpatterninteger = nint(EBSDpatternintd*255.0)
                    EBSDpatternad =  adhisteq(DIFT%nml%nregions,binx,biny,EBSDpatterninteger)
                    binned = float(EBSDpatternad)

                    imagedictflt = 0.0

                    do ll = 1,biny
                        do mm = 1,binx
                            imagedictflt((ll-1)*binx+mm) = binned(mm,ll)
                        end do
                    end do

                    imagedictflt = imagedictflt/vecnorm(imagedictflt)

                    dp = DOT_PRODUCT(tmpimageexpt,imagedictflt)

    ! updating the confidence index if a better match is found
                    if(dp .gt. CIlist(eindex)) then
                        CIlist(eindex) = dp
                        qu = q_T( qdinp = quat%get_quatd() )
                        eu = qu%qe()
                        euler_best(1:3,eindex) = eu%e_copyd()
                    end if

                end do

                if (mod(eindex,250) .eq. 0) then
                    io_int(1) = eindex
                    io_int(2) = totnumexpt
                    call Message%Writevalue('      completed refining pattern #',io_int,2,'(I8,'' of '',I8)')
                end if
            end do

!$OMP END DO

!$OMP CRITICAL
        call memth%dealloc(tmpimageexpt, 'tmpimageexpt', TID=TID)
        call memth%dealloc(binned, 'binned', TID=TID)
        call memth%dealloc(EBSDpatternintd, 'EBSDpatternintd', TID=TID)
        call memth%dealloc(EBSDpatterninteger, 'EBSDpatterninteger', TID=TID)
        call memth%dealloc(EBSDpatternad, 'EBSDpatternad', TID=TID)
        call memth%dealloc(imagedictflt, 'imagedictflt', TID=TID)
!$OMP END CRITICAL
!$OMP BARRIER 
!$OMP END PARALLEL

        stpsz = stpsz/2.D0
        end do

    end do

end if

if (ronl%inRAM.eqv..FALSE.) then
   close(unit=itmpexpt, status='delete')
end if

!===========================================
write (*,*) 'total number of unchanged points : ', unchanged,' out of ', maxeindex
!===========================================

!===========================================
! output section
!===========================================
! add fitted dot product values to HDF5 file
! open the fortran HDF interface
dpfile = trim(EMsoft%getConfigParameter('EMdatapathname'))//trim(ronl%dotproductfile)

! open the fortran HDF interface
hdferr =  HDF%openFile(dpfile)

! add the name list file and the parsed name list for the refinement program
groupname = trim(HDFnames%get_NMLfiles())
hdferr = HDF%createGroup(groupname)

dataset = trim(HDFnames%get_NMLfilename())
hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%nmldeffile)

! leave this group
call HDF%pop()

! open the NMLparameters group to write all the namelist parameters into
hdferr = HDF%createGroup(HDFnames%get_NMLparameters())
call self%writeFitHDFNameList_(HDF, HDFnames)
call HDF%pop()

! open the Scan 1/EBSD/Data group
groupname = 'Scan 1'
hdferr = HDF%openGroup(groupname)
groupname = SC_EBSD
hdferr = HDF%openGroup(groupname)
groupname = SC_Data
hdferr = HDF%openGroup(groupname)

dataset = SC_RefinedDotProducts
call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetFloatArray(dataset, CIlist, Nexp, overwrite)
else
  hdferr = HDF%writeDatasetFloatArray(dataset, CIlist, Nexp)
end if

dataset = SC_RefinedEulerAngles
call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetFloatArray(dataset, sngl(euler_best), 3, Nexp, overwrite)
else
  hdferr = HDF%writeDatasetFloatArray(dataset, sngl(euler_best), 3, Nexp)
end if

call HDF%pop(.TRUE.)

!===========================================
! and generate the ctf/ang output file as well...
dinl%ctffile = ronl%ctffile
dinl%angfile = ronl%angfile

ipar = 0
ipar(1) = 1
ipar(2) = Nexp
ipar(3) = Nexp
ipar(4) = Nexp
ipar(5) = DIDT%FZcnt
ipar(6) = pgnum
if (sum(dinl%ROI).ne.0) then
    ipar(7) = dinl%ROI(3)
    ipar(8) = dinl%ROI(4)
else
    ipar(7) = dinl%ipf_wd
    ipar(8) = dinl%ipf_ht
end if

call mem%alloc(indexmain, (/ ipar(1), ipar(2) /), 'indexmain', 0)
call mem%alloc(resultmain, (/ ipar(1), ipar(2) /), 'resultmain', 0.0)
resultmain(1,1:ipar(2)) = CIlist(1:Nexp)

VT = Vendor_T()
call VT%set_Modality(MPFT%getModality())
if (ronl%ctffile.ne.'undefined') then
  fpar2(1) = mcnl%EkeV
  fpar2(2) = MCsig
  call VT%ctf_writeFile(EMsoft,cell,SG,dinl,ipar,fpar2,indexmain,euler_best*sngl(rtod),resultmain, &
                        DIDT%OSM, DIDT%IQ, noindex=.TRUE.)
  call Message%printMessage(' Data stored in ctf file : '//trim(ronl%ctffile))
end if

if (ronl%angfile.ne.'undefined') then
    fpar1(1) = WD
    call VT%ang_writeFile(EMsoft,cell,SG,dinl,ipar,fpar1,indexmain,euler_best,resultmain,DIDT%IQ)
    call Message%printMessage(' Data stored in ang file : '//trim(ronl%angfile))
end if

! close the fortran HDF5 interface
call closeFortranHDFInterface()

call timer%Time_tock()
tstop = timer%getInterval()

io_real(1) = tstop
call Message%WriteValue('Execution time [system_clock()] = ',io_real,1,"(F14.6,' [s]')")

call mem%allocated_memory_use()

end associate

end subroutine FitOrientation_

end module mod_FitOrientation
