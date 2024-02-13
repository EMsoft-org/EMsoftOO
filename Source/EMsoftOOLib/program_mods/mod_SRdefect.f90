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

module mod_SRdefect
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/13/24
  !!
  !! class definition for the EMSRdefect program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMSRdefect program
type, public :: SRdefectNameListType
  integer(kind=irg) :: DF_npix
  integer(kind=irg) :: DF_npiy
  integer(kind=irg) :: dinfo
  integer(kind=irg) :: t_interval
  integer(kind=irg) :: nthreads
  integer(kind=irg) :: SRG(3)
  integer(kind=irg) :: SRF(3)
  integer(kind=irg) :: Grange
  real(kind=sgl)    :: voltage
  real(kind=sgl)    :: GLaue
  real(kind=sgl)    :: DF_L
  real(kind=sgl)    :: DF_slice
  real(kind=sgl)    :: dmin
  character(4)      :: progmode
  character(3)      :: dispmode
  character(fnlen)  :: outname
  character(fnlen)  :: dispfile
  character(fnlen)  :: xtalname
  character(fnlen)  :: STEMnmlfile
  character(fnlen)  :: defectjsonfile
end type SRdefectNameListType

! class definition
type, public :: SRdefect_T
private 
  character(fnlen)              :: nmldeffile = 'EMSRdefect.nml'
  type(SRdefectNameListType)    :: nml 

contains
private 
  procedure, pass(self) :: setDF_npix_
  procedure, pass(self) :: getDF_npix_
  procedure, pass(self) :: setDF_npiy_
  procedure, pass(self) :: getDF_npiy_
  procedure, pass(self) :: setdinfo_
  procedure, pass(self) :: getdinfo_
  procedure, pass(self) :: sett_interval_
  procedure, pass(self) :: gett_interval_
  procedure, pass(self) :: setnthreads_
  procedure, pass(self) :: getnthreads_
  procedure, pass(self) :: setSRG_
  procedure, pass(self) :: getSRG_
  procedure, pass(self) :: setSRF_
  procedure, pass(self) :: getSRF_
  procedure, pass(self) :: setGrange_
  procedure, pass(self) :: getGrange_
  procedure, pass(self) :: setvoltage_
  procedure, pass(self) :: getvoltage_
  procedure, pass(self) :: setGLaue_
  procedure, pass(self) :: getGLaue_
  procedure, pass(self) :: setDF_L_
  procedure, pass(self) :: getDF_L_
  procedure, pass(self) :: setDF_slice_
  procedure, pass(self) :: getDF_slice_
  procedure, pass(self) :: setdmin_
  procedure, pass(self) :: getdmin_
  procedure, pass(self) :: setprogmode_
  procedure, pass(self) :: getprogmode_
  procedure, pass(self) :: setdispmode_
  procedure, pass(self) :: getdispmode_
  procedure, pass(self) :: setoutname_
  procedure, pass(self) :: getoutname_
  procedure, pass(self) :: setdispfile_
  procedure, pass(self) :: getdispfile_
  procedure, pass(self) :: setxtalname_
  procedure, pass(self) :: getxtalname_
  procedure, pass(self) :: setSTEMnmlfile_
  procedure, pass(self) :: getSTEMnmlfile_
  procedure, pass(self) :: setdefectjsonfile_
  procedure, pass(self) :: getdefectjsonfile_
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: SRdefect_

  generic, public :: setDF_npix => setDF_npix_
  generic, public :: getDF_npix => getDF_npix_
  generic, public :: setDF_npiy => setDF_npiy_
  generic, public :: getDF_npiy => getDF_npiy_
  generic, public :: setdinfo => setdinfo_
  generic, public :: getdinfo => getdinfo_
  generic, public :: sett_interval => sett_interval_
  generic, public :: gett_interval => gett_interval_
  generic, public :: setnthreads => setnthreads_
  generic, public :: getnthreads => getnthreads_
  generic, public :: setSRG => setSRG_
  generic, public :: getSRG => getSRG_
  generic, public :: setSRF => setSRF_
  generic, public :: getSRF => getSRF_
  generic, public :: setGrange => setGrange_
  generic, public :: getGrange => getGrange_
  generic, public :: setvoltage => setvoltage_
  generic, public :: getvoltage => getvoltage_
  generic, public :: setGLaue => setGLaue_
  generic, public :: getGLaue => getGLaue_
  generic, public :: setDF_L => setDF_L_
  generic, public :: getDF_L => getDF_L_
  generic, public :: setDF_slice => setDF_slice_
  generic, public :: getDF_slice => getDF_slice_
  generic, public :: setdmin => setdmin_
  generic, public :: getdmin => getdmin_
  generic, public :: setprogmode => setprogmode_
  generic, public :: getprogmode => getprogmode_
  generic, public :: setdispmode => setdispmode_
  generic, public :: getdispmode => getdispmode_
  generic, public :: setoutname => setoutname_
  generic, public :: getoutname => getoutname_
  generic, public :: setdispfile => setdispfile_
  generic, public :: getdispfile => getdispfile_
  generic, public :: setxtalname => setxtalname_
  generic, public :: getxtalname => getxtalname_
  generic, public :: setSTEMnmlfile => setSTEMnmlfile_
  generic, public :: getSTEMnmlfile => getSTEMnmlfile_
  generic, public :: setdefectjsonfile => setdefectjsonfile_
  generic, public :: getdefectjsonfile => getdefectjsonfile_
  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: SRdefect => SRdefect_

end type SRdefect_T

! the constructor routine for this class 
interface SRdefect_T
  module procedure SRdefect_constructor
end interface SRdefect_T

contains

!--------------------------------------------------------------------------
type(SRdefect_T) function SRdefect_constructor( nmlfile ) result(SRdefect)
!! author: MDG 
!! version: 1.0 
!! date: 02/13/24
!!
!! constructor for the SRdefect_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen),INTENT(IN)                     :: nmlfile 

call SRdefect%readNameList(nmlfile)


end function SRdefect_constructor

!--------------------------------------------------------------------------
subroutine SRdefect_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 02/13/24
!!
!! destructor for the SRdefect_T Class
 
IMPLICIT NONE

type(SRdefect_T), INTENT(INOUT)  :: self 

call reportDestructor('SRdefect_T')

end subroutine SRdefect_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/13/24
!!
!! read the namelist from an nml file for the SRdefect_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(SRdefect_T), INTENT(INOUT)     :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg) :: DF_npix
integer(kind=irg) :: DF_npiy
integer(kind=irg) :: dinfo
integer(kind=irg) :: t_interval
integer(kind=irg) :: nthreads
integer(kind=irg) :: SRG(3)
integer(kind=irg) :: SRF(3)
integer(kind=irg) :: Grange
real(kind=sgl)    :: voltage
real(kind=sgl)    :: GLaue
real(kind=sgl)    :: DF_L
real(kind=sgl)    :: DF_slice
real(kind=sgl)    :: dmin
character(4)      :: progmode
character(3)      :: dispmode
character(fnlen)  :: outname
character(fnlen)  :: dispfile
character(fnlen)  :: xtalname
character(fnlen)  :: STEMnmlfile
character(fnlen)  :: defectjsonfile

! define the IO namelist to facilitate passing variables to the program.
namelist /SRdeflist/ DF_L, DF_npix, DF_npiy, DF_slice, dmin, progmode,&
                    dinfo, outname, t_interval, dispfile, SRF, &
                    nthreads, dispmode, xtalname, voltage, SRG, Grange, GLaue, &
                    STEMnmlfile, defectjsonfile

nthreads = 1
voltage = 200
progmode = 'STEM'
SRG = (/ 1, 0, 0 /)
SRF = (/ 0, 0, 1 /)
Grange = 4
GLaue = 0.5
dispmode = 'not'
dinfo = 0 ! 1 is verbose
t_interval = 5 ! update every x steps
DF_L = 1.0
DF_npix = 256
DF_npiy = 256
DF_slice = 1.0 ! slice thickness for scattering matrix approach (nmu)
dmin = 0.03
xtalname = 'undefined'
defectjsonfile = 'undefined'
outname = 'undefined'
dispfile = 'undefined'
STEMnmlfile = 'STEM_rundata.nml'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=SRdeflist)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(xtalname).eq.'undefined') then
  call Message%printError('readNameList:',' xtalname file name is undefined in '//nmlfile)
 end if

  if (trim(defectjsonfile).eq.'undefined') then
  call Message%printError('readNameList:',' defectjsonfile file name is undefined in '//nmlfile)
 end if

 if (trim(outname).eq.'undefined') then
  call Message%printError('readNameList:',' outname file name is undefined in '//nmlfile)
 end if

 if (trim(dispfile).eq.'undefined') then
  call Message%printError('readNameList:',' dispfile file name is undefined in '//nmlfile)
 end if

 if (trim(STEMnmlfile).eq.'undefined') then
  call Message%printError('readNameList:',' STEMnmlfilefile name is undefined in '//nmlfile)
 end if
end if

self%nml%DF_L = DF_L
self%nml%DF_npix = DF_npix
self%nml%DF_npiy = DF_npiy
self%nml%DF_slice = DF_slice
self%nml%dmin = dmin
self%nml%progmode = progmode
self%nml%dinfo = dinfo
self%nml%outname = outname
self%nml%t_interval = t_interval
self%nml%dispfile = dispfile
self%nml%nthreads = nthreads
self%nml%dispmode = dispmode
self%nml%xtalname = xtalname
self%nml%voltage = voltage
self%nml%SRG = SRG
self%nml%SRF = SRF
self%nml%Grange = Grange
self%nml%GLaue = GLaue
self%nml%STEMnmlfile = STEMnmlfile
self%nml%defectjsonfile = defectjsonfile

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/13/24
!!
!! pass the namelist for the SRdefect_T Class to the calling program

IMPLICIT NONE 

class(SRdefect_T), INTENT(INOUT)          :: self
type(SRdefectNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/13/24
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)        :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 6, n_real = 5
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( enl => self%nml )


! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! write all the single integers
io_int = (/ enl%DF_npix, enl%DF_npiy, enl%dinfo, enl%t_interval, enl%nthreads, enl%Grange/)
intlist(1) = 'DF_npix'
intlist(2) = 'DF_npiy'
intlist(3) = 'dinfo'
intlist(4) = 't_interval'
intlist(5) = 'nthreads'
intlist(6) = 'Grange'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ enl%voltage, enl%GLaue, enl%DF_L, enl%DF_slice, enl%dmin /)
reallist(1) = 'voltage'
reallist(2) = 'GLaue'
reallist(3) = 'DF_L'
reallist(4) = 'DF_slice'
reallist(5) = 'dmin'
call HDF%writeNMLreals(io_real, reallist, n_real)

! a 3-vector
dataset = 'SRG'
hdferr = HDF%writeDatasetIntegerArray(dataset, enl%SRG, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create SRG dataset', hdferr)

dataset = 'SRF'
hdferr = HDF%writeDatasetIntegerArray(dataset, enl%SRF, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create SRF dataset', hdferr)

! write all the strings
dataset = 'progmode'
line2(1) = trim(enl%progmode)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create progmode dataset', hdferr)

dataset = 'dispmode'
line2(1) = trim(enl%dispmode)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create dispmode dataset', hdferr)

dataset = 'outname'
line2(1) = trim(enl%outname)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create outname dataset', hdferr)

dataset = 'dispfile'
line2(1) = trim(enl%dispfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create dispfile dataset', hdferr)

dataset = 'xtalname'
line2(1) = trim(enl%xtalname)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create xtalname dataset', hdferr)

dataset = 'STEMnmlfile'
line2(1) = trim(enl%STEMnmlfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create STEMnmlfile dataset', hdferr)

dataset = 'defectjsonfile'
line2(1) = trim(enl%defectjsonfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create defectjsonfile dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine setDF_npix_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setDF_npix_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set DF_npix in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%DF_npix = inp

end subroutine setDF_npix_

!--------------------------------------------------------------------------
function getDF_npix_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getDF_npix_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get DF_npix from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%DF_npix

end function getDF_npix_

!--------------------------------------------------------------------------
subroutine setDF_npiy_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setDF_npiy_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set DF_npiy in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%DF_npiy = inp

end subroutine setDF_npiy_

!--------------------------------------------------------------------------
function getDF_npiy_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getDF_npiy_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get DF_npiy from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%DF_npiy

end function getDF_npiy_

!--------------------------------------------------------------------------
subroutine setdinfo_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdinfo_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set dinfo in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%dinfo = inp

end subroutine setdinfo_

!--------------------------------------------------------------------------
function getdinfo_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdinfo_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get dinfo from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%dinfo

end function getdinfo_

!--------------------------------------------------------------------------
subroutine sett_interval_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: sett_interval_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set t_interval in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%t_interval = inp

end subroutine sett_interval_

!--------------------------------------------------------------------------
function gett_interval_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: gett_interval_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get t_interval from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%t_interval

end function gett_interval_

!--------------------------------------------------------------------------
subroutine setnthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnthreads_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set nthreads in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%nthreads = inp

end subroutine setnthreads_

!--------------------------------------------------------------------------
function getnthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnthreads_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get nthreads from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%nthreads

end function getnthreads_

!--------------------------------------------------------------------------
subroutine setSRG_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setSRG(3)_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set SRG in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp(3)

self%nml%SRG = inp

end subroutine setSRG_

!--------------------------------------------------------------------------
function getSRG_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getSRG_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get SRG from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out(3)

out = self%nml%SRG

end function getSRG_

!--------------------------------------------------------------------------
subroutine setSRF_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setSRF(3)_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set SRF in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp(3)

self%nml%SRF = inp

end subroutine setSRF_

!--------------------------------------------------------------------------
function getSRF_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getSRF_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get SRF from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out(3)

out = self%nml%SRF

end function getSRF_

!--------------------------------------------------------------------------
subroutine setGrange_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setGrange_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set Grange in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%Grange = inp

end subroutine setGrange_

!--------------------------------------------------------------------------
function getGrange_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getGrange_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get Grange from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%Grange

end function getGrange_

!--------------------------------------------------------------------------
subroutine setvoltage_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setvoltage_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set voltage in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%voltage = inp

end subroutine setvoltage_

!--------------------------------------------------------------------------
function getvoltage_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getvoltage_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get voltage from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%voltage

end function getvoltage_

!--------------------------------------------------------------------------
subroutine setGLaue_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setGLaue_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set GLaue in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%GLaue = inp

end subroutine setGLaue_

!--------------------------------------------------------------------------
function getGLaue_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getGLaue_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get GLaue from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%GLaue

end function getGLaue_

!--------------------------------------------------------------------------
subroutine setDF_L_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setDF_L_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set DF_L in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%DF_L = inp

end subroutine setDF_L_

!--------------------------------------------------------------------------
function getDF_L_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getDF_L_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get DF_L from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%DF_L

end function getDF_L_

!--------------------------------------------------------------------------
subroutine setDF_slice_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setDF_slice_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set DF_slice in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%DF_slice = inp

end subroutine setDF_slice_

!--------------------------------------------------------------------------
function getDF_slice_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getDF_slice_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get DF_slice from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%DF_slice

end function getDF_slice_

!--------------------------------------------------------------------------
subroutine setdmin_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdmin_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set dmin in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%dmin = inp

end subroutine setdmin_

!--------------------------------------------------------------------------
function getdmin_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdmin_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get dmin from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%dmin

end function getdmin_

!--------------------------------------------------------------------------
subroutine setprogmode_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setprogmode_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set progmode in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
character(4), INTENT(IN)       :: inp

self%nml%progmode = trim(inp)

end subroutine setprogmode_

!--------------------------------------------------------------------------
function getprogmode_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getprogmode_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get progmode from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
character(4)                   :: out

out = trim(self%nml%progmode)

end function getprogmode_

!--------------------------------------------------------------------------
subroutine setdispmode_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdispmode_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set dispmode in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
character(3), INTENT(IN)       :: inp

self%nml%dispmode = trim(inp)

end subroutine setdispmode_

!--------------------------------------------------------------------------
function getdispmode_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdispmode_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get dispmode from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
character(3)                   :: out

out = trim(self%nml%dispmode)

end function getdispmode_

!--------------------------------------------------------------------------
subroutine setoutname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setoutname_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set outname in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%outname = trim(inp)

end subroutine setoutname_

!--------------------------------------------------------------------------
function getoutname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getoutname_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get outname from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%outname)

end function getoutname_

!--------------------------------------------------------------------------
subroutine setdispfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdispfile_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set dispfile in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%dispfile = trim(inp)

end subroutine setdispfile_

!--------------------------------------------------------------------------
function getdispfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdispfile_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get dispfile from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%dispfile)

end function getdispfile_

!--------------------------------------------------------------------------
subroutine setxtalname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setxtalname_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set xtalname in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%xtalname = trim(inp)

end subroutine setxtalname_

!--------------------------------------------------------------------------
function getxtalname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getxtalname_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get xtalname from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%xtalname)

end function getxtalname_

!--------------------------------------------------------------------------
subroutine setSTEMnmlfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setSTEMnmlfile_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set STEMnmlfile in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%STEMnmlfile = trim(inp)

end subroutine setSTEMnmlfile_

!--------------------------------------------------------------------------
function getSTEMnmlfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getSTEMnmlfile_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get STEMnmlfile from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%STEMnmlfile)

end function getSTEMnmlfile_

!--------------------------------------------------------------------------
subroutine setdefectjsonfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdefectjsonfile_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set defectjsonfile in the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%defectjsonfile = trim(inp)

end subroutine setdefectjsonfile_

!--------------------------------------------------------------------------
function getdefectjsonfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdefectjsonfile_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get defectjsonfile from the SRdefect_T class

IMPLICIT NONE

class(SRdefect_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%defectjsonfile)

end function getdefectjsonfile_

!--------------------------------------------------------------------------
subroutine SRdefect_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: SRdefect_
!! author: MDG 
!! version: 1.0 
!! date: 02/13/24
!!
!! perform the computations

use mod_EMsoft
use mod_HDFnames
use mod_STEM 
use mod_initializers
use mod_crystallography
use mod_symmetry
use mod_io 
use mod_defect
use HDF5
use mod_HDFsupport
use ISO_C_BINDING
use omp_lib
use mod_math
use mod_rotations
use mod_gvectors
use mod_kvectors
use mod_diffraction
use mod_timing
use mod_memory
use stringconstants

IMPLICIT NONE 

class(SRdefect_T), INTENT(INOUT) :: self
type(EMsoft_T), INTENT(INOUT)    :: EMsoft
character(fnlen), INTENT(INOUT)  :: progname 
type(HDFnames_T), INTENT(INOUT)  :: HDFnames

type(STEM_T)                     :: STEM 
type(gvectors_T)                 :: gvec 
type(kvectors_T)                 :: kvec
type(SpaceGroup_T)               :: SG
type(Cell_T)                     :: cell
type(Defect_T)                   :: defects
type(Diffraction_T)              :: Diff
type(IO_T)                       :: Message
type(gnode)                      :: rlp
type(Timing_T)                   :: timer
type(DynType)                    :: Dyn
type(memory_T)                   :: mem, memth
type(HDF_T)                      :: HDF

integer(kind=irg)                :: nn,izero,i,j,k,npix,npiy,ii,jj,numvoids,numdisl,numYdisl,numsf, skip, &
                                    numinc,dinfo,t_interval,DF_nums_new, io_int(3), Grange, hdferr, &
                                    DF_npix_new,DF_npiy_new, numstart,numstop, isg, TID, NTHR, numCL, iCL, &
                                    SRG(3), SETNTHR, iSTEM, ic, gg(3), ir, DF_npix,DF_npiy, DF_nums, & 
                                    error_cnt, ier, isym, Nmat
real(kind=dbl)                   :: DynFN(3), dgr, arg, kv(3), qv(3), lambda
real(kind=sgl)                   :: ind(3),hkl(3),thick, qx, qy, c(3), att,xgp,DF_gf(3),io_real(1), voltage, &
                                    frac, GLaue, gx(3),gy(3),imat, DF_gc(3), DF_gstar(3), DF_slice, gdotR
character(fnlen)                 :: dataname,sgname, dispfile,xtalname, foilnmlfile, STEMnmlfile, dataset, datagroupname
character(4)                     :: outputformat, dispmode, progmode
complex(kind=dbl)                :: czero=cmplx(0.0,0.0,dbl),cone=cmplx(1.0,0.0,dbl)
logical                          :: verbose = .TRUE., f_exists
character(11)                    :: dstr
character(15)                    :: tstrb
character(15)                    :: tstre

complex(kind=dbl),allocatable    :: DHWM(:,:),Afirst(:,:),DHWMvoid(:,:),Azz(:,:)
complex(kind=dbl),allocatable    :: amp(:),amp2(:)
real(kind=sgl),allocatable       :: weights(:), inten(:), STEMimages(:,:,:,:), sgarray(:,:), qxy(:,:)
integer(kind=irg),allocatable    :: disparray(:,:,:), ggg(:,:)
real(kind=sgl),allocatable       :: BFweightsarray(:,:,:),ADFweightsarray(:,:,:)
real(kind=sgl),allocatable       :: DF_R(:,:)
complex(kind=dbl),allocatable    :: DF_Sarray(:,:,:), theta(:), DF_Svoid(:,:), DHWMz(:,:)
complex(kind=dbl),allocatable    :: Az(:,:)

! initialize the STEM geometry namelist as well 
STEM = STEM_T( self%getSTEMnmlfile() )

call openFortranHDFInterface()
HDF = HDF_T()

associate( nml => self%nml )

! initialize the timing routines
timer = Timing_T()
tstrb = timer%getTimeString()

call cell%setFileName(nml%xtalname)

call Diff%setrlpmethod('WK')
call Diff%setV(dble(nml%voltage))

call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, nml%dmin, verbose, useHDF=HDF, noLUT=.TRUE.)
call Diff%Printrlp()
lambda = Diff%getWaveLength()

! determine the point group number
j=0
do i=1,32
  if (SGPG(i).le.SG%getSpaceGroupNumber()) j=i
end do
call SG%BFsymmetry(self%getSRG(),j,isym,ir)

SRG = self%getSRG()
! use systematic row vector to compute G* as in eq. 8.28.
DF_gf = float(SRG)
DF_gstar = DF_gf/cell%CalcLength(DF_gf,'r')**2    ! define G* such that G.G* = 1
call cell%TransSpace(DF_gf,DF_gc,'r','c')         ! convert to Cartesian reference frame

! convert some parameters from namelist to local variables. 
DF_npix = self%getDF_npix()
DF_npiy = self%getDF_npiy()
npix = DF_npix
npiy = DF_npiy
Grange = self%getGrange()
progmode = self%getprogmode()
SETNTHR = self%getnthreads()
defects%DF_slice = self%getDF_slice()
defects%DF_npix = DF_npix
defects%DF_npiy = DF_npiy
defects%DF_L = self%getDF_L()
kv = dble(self%getSRF())
qv = dble(SRG)
DynFN = self%getSRF() 
DF_slice = self%getDF_slice()

! read foil and defect data 
call defects%InitializeDefects(EMsoft,cell,self%getdefectjsonfile(),npix,npiy, &
                               self%getDF_L(),DF_gf,kv,qv,error_cnt,verbose)

! compute total number of beams
nn = 2*Grange+1  ! total number of beams
izero = (nn+1)/2 ! id number of g=0 beam 

! initialize all STEM related arrays for the systematic row case
call STEM%init_STEM(cell,Diff,nn,SRG)  ! used to be a call to read_STEM_data

! allocate and initialize various arrays
mem = memory_T()
call mem%alloc(DF_Sarray, (/ Nmat-1,nn,nn /), 'DF_Sarray', initval=czero, startdims=(/0,1,1/) )
call mem%alloc(theta, (/nn /), 'theta',initval=czero, startdims=(/-nn/) )
call mem%alloc(DF_Svoid, (/ nn,nn /), 'DF_Svoid', initval=czero )
call mem%alloc(DHWMz, (/ nn,nn /), 'DF_SvDHWMzoid', initval=czero )
call mem%alloc(DHWM, (/ nn,nn /), 'DHWM', initval=czero )
call mem%alloc(DHWMvoid, (/ nn,nn /), 'DHWMvoid', initval=czero )

! Compute the off-diagonal part of the complex DHW matrix (factor i is included)
! We can precompute those because they will not change at all during the run
!       (these lines implement the equations on page 476 of the EM book)
ind = float(SRG)
do i=1,nn
 do j=1,nn
  hkl=(-Grange+i-1)*ind-(-Grange+j-1)*ind     ! difference vector
  if (i.ne.j) then
   call Diff%CalcUcg(cell,int(hkl), applyqgshift=.TRUE.)    ! compute the interaction parameters
   rlp = Diff%getrlp()
   DHWMz(i,j) = cPi*cmplx(-aimag(rlp%qg),real(rlp%qg),dbl)  ! and initalize the off-diagonal matrix element
  else
   DHWMz(i,j) = czero                         ! for now at least; this will be filled in later
  endif
 end do
end do
call Message%printMessage(' Reference Darwin-Howie-Whelan matrix initialized')

! display the diffraction information for the fundamental reflection of the systematic row
call Diff%CalcUcg(cell,SRG)
rlp = Diff%getrlp()
io_real(1) = rlp%xg
call Message%WriteValue('Extinction distance for g : ', io_real, 1, "(F10.5)")
io_real(1) = rlp%xgp
call Message%WriteValue('Anomalous absorption length for g : ', io_real, 1, "(F10.5)")
io_real(1) = rlp%xgp/rlp%xg
call Message%WriteValue('Absorption Ratio : ', io_real, 1, "(F10.5)")
! compute the normal absorption factor xgp
call Diff%CalcUcg(cell,(/0,0,0/))
rlp = Diff%getrlp()
xgp = aimag(rlp%qg)
io_real(1) = 1.0/xgp
call Message%WriteValue('Normal absorption length : ', io_real, 1, "(F10.5/)")

! define the foil thickness, attenuation, and number slices per column
thick = defects%foil%zb    ! this is the same everywhere for this version; needs to be updated in the next version
att = exp(-2.0*cPi*thick*xgp)  ! this is the global attenuation factor; remove the factor 2.0 if amplitudes are needed
DF_nums = nint(thick/defects%DF_slice)  ! this is the number of slices for this particular column
defects%DF_nums = DF_nums

! setup the defect displacement field parameters
if ((self%nml%dispmode.eq.'new').or.(self%nml%dispmode.eq.'not')) then
  call Message%printMessage(' Starting Displacement Field Computation (multi-threaded)')

! precompute ALL the defect columns and, if needed, store them in dispfile
! this portion should be carried out in multi-threaded mode as much as possible
  call mem%alloc(disparray, (/ DF_nums,DF_npix,DF_npiy /), 'disparray', initval=0  )

!!==========================================================================
!!==========================================================================

! initiate multi-threaded segment
!!$OMP     PARALLEL PRIVATE(TID,DF_R,imatvals,gdotR,i,j,k,imat) &
!!$OMP&   SHARED(NTHR,DF_npix,DF_npiy,defects%DF_nums,defects%numvoids,defects%numdisl,defects%numsf,defects%numinc,disparray,t_interval)
!  NTHR = OMP_GET_NUM_THREADS()
!  TID = OMP_GET_THREAD_NUM()
!  write (*,*) TID,': entering parallel region'
  TID = 0
  NTHR = 1
!  if (TID.eq.0) then
! do time reporting only in the master thread
!    call Time_report(TT,t_interval*float(NTHR)/float(DF_npix))
!    call Time_start(TT)
!  end if
  call mem%alloc(DF_R, (/ DF_nums,3 /), 'DF_R')     ! each thread has its own DF_R array
  call mem%alloc(defects%DF_R, (/ DF_nums,3 /), 'defects%DF_R')     ! each thread has its own DF_R array  

!!$OMP DO SCHEDULE (GUIDED)
  ! write(*,*) TID,': starting Do Schedule'
  do i=1,DF_npix  
    do j=1,DF_npiy
      DF_R = 0.0
! compute the displacement vectors DF_R for all points in the column
      call defects%CalcR(cell,i,j)
! loop over the fixed thickness slices
      do k=1,DF_nums
! then convert to the dot-product 
       if (defects%DF_R(k,1).eq.-10000.0) then  ! this is point inside a void
         disparray(k,i,j) = -10000
       else  ! it is not a void, so use the full dot product g.R (all vectors must be Cartesian !)
! use gac and gbc to get the two dot products and store both of them as integers mapped onto the 0..numd range
         gdotR = Dot_Product(DF_gc,defects%DF_R(k,1:3))
         imat = nint(float(Nmat)*amod(sngl(gdotR+10000.D0),1.0))
         if (imat.lt.0) imat = 0
         if (imat.gt.Nmat) imat = Nmat - 1
       end if
       disparray(k,i,j) = imat
     end do ! k loop
   end do
!  if ((mod(i,t_interval).eq.0).and.(TID.eq.0)) call Time_remaining(TT,i,defects%DF_npix)
end do
!!$OMP END DO 
! if (TID.eq.0) call Time_stop(TT,defects%DF_npix*defects%DF_npiy)
!!$OMP END PARALLEL
end if
call Message%printMessage('   -> done.')

!!==========================================================================
!!==========================================================================
! to be reviewed and modified if necessary

! and, if needed, store the defect displacement field for re-runs
! if (self%nml%dispmode.ne.'not') then
!   if ((self%nml%dispfile.ne.'none').and.(self%nml%dispmode.eq.'new')) then 
!     call Message%printMessage('Displacement field data stored in file '//self%nml%dispfile)
!     open(unit=dataunit,file=dispfile,status='new',action='write',form='unformatted')
!     write (dataunit) DF_nums,DF_npix,DF_npiy
!     write (dataunit) disparray
!     close(unit=dataunit,status='keep',iostat=ier)
!   endif
!   if ((self%nml%dispfile.ne.'none').and.(self%nml%dispmode.eq.'old')) then ! there is a pre-computed defect file, so let's load it
!    allocate(disparray(DF_nums,DF_npix,DF_npiy))
!    disparray = 0
!    write (*,*) 'Opening old dispfile '//trim(self%nml%dispfile)
!    write (*,*) 'shape(disparray) = ',shape(disparray)
 
!    open(unit=dataunit,file=trim(dispfile),status='old',action='read',form='unformatted')
!    read (dataunit) DF_nums_new,DF_npix_new,DF_npiy_new
! ! check to make sure that these dimensions are the same as the ones used in the current run of the program
!    if ((DF_nums_new.ne.DF_nums).or.(DF_npix_new.ne.DF_npix).or.(DF_npiy_new.ne.DF_npiy)) then
!     io_int(1) = DF_nums_new; io_int(2) = DF_npix_new; io_int(3) = DF_npiy_new
!     call WriteValue('The dimensions of the defect array in the file are : ', io_int, 3, "(3I5)")
!     io_int(1) = DF_nums; io_int(2) = DF_npix; io_int(3) = DF_npiy
!     call WriteValue('The dimensions in the SRdef_rundata file do not agree : ', io_int, 3, "(3I5)")
!     call Message('Terminating program run')
!     stop
!   end if
! ! ok, we're good, so read the actual data...  
!   read (dataunit) disparray
!   close(unit=dataunit,status='keep')
!  end if
! end if

! next the STEMimages array; for consistency with EMZAdefect, we've changed the index ordering)
if (self%nml%progmode.eq.'STEM') then
  call mem%alloc(STEMimages, (/ npix,npiy,nn,STEM%getnumberofsvalues() /), 'STEMimages', initval=0.0)
else
  call mem%alloc(STEMimages, (/ 2,npix,npiy,STEM%getnumCL() /), 'STEMimages', initval=0.0)
end if

! allocate and initialize auxiliary variables 
call mem%alloc(Afirst, (/ nn,nn /), 'Afirst', initval=czero)
  
! initialize the timer
numstart = 1
numstop = STEM%getnumberofsvalues()    

! before we start the main computational loop, we need to store the 
! necessary variables for the reconstruction of STEM BF/HAADF images
! using the post-processing IDL routine.  This file should contain all
! the information needed to recreate the full CBED patterns at each image
! pixel.  This is then followed by the actual data.

! the file format is identical to that of the EMZAdefect program in STEM mode,
! so that the STEMDisplay visualization program can be used for both ZA and SR files.

! In the present version we change that format to HDF5

call Message%printMessage(' Storing data for IDL visualization program in '//self%nml%outname)

! get the filename; if it already exists, then delete it and create a new one
dataname = EMsoft%generateFilePath('EMdatapathname', self%nml%outname)
inquire(file=trim(dataname), exist=f_exists)

if (f_exists) then
  open(unit=dataunit, file=trim(dataname), status='old',form='unformatted')
  close(unit=dataunit, status='delete')
end if

! Create a new file using the default properties.
hdferr = HDF%createFile(dataname)

! write the EMheader to the file
datagroupname = trim(HDFnames%get_ProgramData()) 
call HDF%writeEMheader(EMsoft, dstr, tstrb, tstre, progname, datagroupname)

! add the CrystalData group at the top level of the file
call cell%addXtalDataGroup(SG, EMsoft, HDF)

! create a namelist group to write all the namelist files into
hdferr = HDF%createGroup(HDFnames%get_NMLfiles())

! read the text file and write the array to the file
dataset = trim(HDFnames%get_NMLfilename())
hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%nmldeffile)

! we also need to include the defect JSON file here in this group 
dataset = 'defectJSONfile'
hdferr = HDF%writeDatasetTextFile(dataset, self%nml%defectjsonfile)

! leave this group
call HDF%pop()

! create a namelist group to write all the namelist files into
hdferr = HDF%createGroup(HDFnames%get_NMLparameters())
call self%writeHDFNameList(HDF, HDFnames)

! we also need to write the other namelist file for the STEM_T class ... 
call HDFnames%set_NMLlist(SC_STEMGeometryNameList)
call STEM%writeHDFNameList(HDF, HDFnames)
call HDFnames%set_NMLlist(SC_SRdefectNameList)

! leave this group
call HDF%pop()

! then the remainder of the data in the EMData group
hdferr = HDF%createGroup(HDFnames%get_EMData())
hdferr = HDF%createGroup(HDFnames%get_ProgramData())

! bragg angle for the SRG reflection
dataset = 'Bragg_angle_SRG'
  hdferr = HDF%writeDatasetFloat(dataset, Diff%CalcDiffAngle(cell, SRG)*0.5 )

! number of reflections, and associated information (hkl, ...)
  call cell%TransSpace(sngl(defects%foil%F),c,'d','c')
  call cell%NormVec(c,'c')
! then make ga the x-axis
  call cell%TransSpace(float(SRG),gx,'r','c')
  call cell%NormVec(gx,'c')
! compute the cross product between k and gx; this is the y-axis
  call cell%CalcCross(c,gx,gy,'c','c',0)
  
! number of reflections in systematic row
dataset = 'numref'
  hdferr = HDF%writeDatasetInteger(dataset, nn)

! create some temporary arrays
call mem%alloc(ggg, (/ 3, Grange /), 'ggg', initval=0, startdims=(/ 1, -Grange /) )
call mem%alloc(qxy, (/ 2, Grange /), 'qxy', initval=0.0, startdims=(/ 1, -Grange /) )
do ic = -Grange,Grange
  ggg(1:3,ic) = ic * SRG(1:3)
  call cell%TransSpace(float(ggg(1:3,ic)),c,'r','c')
  qxy(1:2,ic) = (/ cell%CalcDot(c, gx, 'c'), 0.0 /)
end do
dataset = 'gg'
  hdferr = HDF%writeDatasetIntegerArray(dataset, ggg, 3, 2*Grange+1 )

dataset = 'qxy'
  hdferr = HDF%writeDatasetFloatArray(dataset, qxy, 2, 2*Grange+1 )

call mem%dealloc(ggg, 'ggg')
call mem%dealloc(qxy, 'qxy')

! CLarray
dataset = 'CLarray'
  hdferr = HDF%writeDatasetFloatArray(dataset, STEM%getCLarray(), 20)

! ok, that's it for the output for now ... let's do the major loop
! first we copy the weights arrays since they are private to the STEM class
BFweightsarray = STEM%getBFweightsarray()
ADFweightsarray = STEM%getADFweightsarray()
sgarray = STEM%getsgarray()
numCL = STEM%getnumCL()

call mem%alloc(Az, (/ nn,nn /), 'Az', initval=czero)
iSTEM = 0

mainloop: do isg = numstart,numstop   ! this is the main computational loop
! get the correct excitation errors for this beam orientation
! fill the diagonal of the reference dynamical matrix and the void matrix
  do i=1,nn
   DHWMz(i,i)=2.0*cPi*cmplx(0.0,sgarray(i,isg))    ! initialize the diagonal element of the dynamical matrix
   DHWMvoid(i,i) = DHWMz(i,i)
  end do

! compute the first slice scattering matrix from the exponential Taylor expansion  (eq. 5.27)
! for the void (i.e., sort of a vacuum propagator)
  Afirst = czero
  do i=1,nn
    Afirst(i,i) = cone    ! initialize Afirst to be the identity matrix
  end do
  call MatrixExponential(Afirst, DF_Svoid, dble(DF_slice), 'Pade', nn)  

! and for the "non-void" 
! main loop for the array of Nmat scattering matrices
  do k=0,Nmat-1  
! initialize theta array
    do i=-nn,nn 
      arg = dble(i)*dble(k)*dgr
      theta(i) = cmplx(dcos(arg),-dsin(arg))
    end do
! then multiply DHWMz by appropriate theta values  
    do i=-Grange,Grange
      do j=-Grange,Grange
        DHWM(i+Grange+1,j+Grange+1) = DHWMz(i+Grange+1,j+Grange+1)*theta(i-j)
      end do
    end do

    call MatrixExponential(DHWM, Az, dble(DF_slice), 'Pade', nn)  

! and store in the main array
    DF_Sarray(k,1:nn,1:nn) = Az(1:nn,1:nn)
end do  ! main loop

call Message%printMessage(' Scattering matrices computed')

!------------------------------------------------!
! Finally, here it is: the actual (threaded!) image computation  !
!------------------------------------------------!
! NTHR = OMP_GET_NUM_THREADS()
call OMP_SET_NUM_THREADS(self%getnthreads())
memth = Memory_T( nt = self%getnthreads() )

!$OMP    PARALLEL PRIVATE(TID,i,j,k,ii,jj,iCL,amp,amp2,Azz,inten,weights,ic) &
!$OMP&   SHARED(NTHR,npix,npiy,DF_nums,disparray,DF_Sarray,DF_Svoid,progmode,Nmat, &
!$OMP&   BFweightsarray,ADFweightsarray,att,STEMimages,t_interval,nn,numCL,izero,isg)
!$OMP     

TID = OMP_GET_THREAD_NUM()   
call memth%alloc(Azz, (/ nn,nn /),'Azz',initval = czero,TID=TID)
call memth%alloc(amp, (/ nn /),'amp',initval = czero,TID=TID)
! call memth%alloc(amp2, (/ nn /),'amp2',initval = czero,TID=TID)
call memth%alloc(weights, (/ nn /),'weights',initval = 0.0,TID=TID)
call memth%alloc(inten, (/ nn /),'inten',initval=0.0,TID=TID)

!$OMP DO SCHEDULE (GUIDED)
donpix: do i=1,npix
  donpiy:   do j=1,npiy

! initialize the wave function for this pixel with (1.0,0.0) for the incident beam izero
    amp = czero
    amp(izero) = cmplx(1.D0,0.D0)

! loop over the fixed thickness slices
      doslices: do k=1,DF_nums
! select the appropriate scattering matrix to propagate with (see section 8.3.3 in the book)
       if (disparray(k,i,j).eq.-10000) then  ! this is a point inside a void
         Azz = DF_Svoid    ! so we use the void propagator matrix
       else  ! it is not a void
         Azz = DF_Sarray(disparray(k,i,j),1:nn,1:nn)
       end if
! and multiply with this matrix
       amp = matmul(Az,amp)
! alternatively we do the computation explicitly
       ! amp2 = czero
       ! do ii=1,nn
       !  do jj=1,nn
       !   amp2(ii) = amp2(ii) + Azz(ii,jj) * amp(jj)
       !  end do
       ! end do
       ! amp = amp2
      end do doslices ! loop over slices
   
! compute the (attenuated) intensities for the EM and STEM images and store
      inten(1:nn) = att*abs(amp(1:nn))**2
      if (progmode.eq.'STEM') then 
          STEMimages(i,j,1:nn,isg) = inten
      end if
      if (progmode.eq.'BFDF') then 
          do iCL=1,numCL
! BF detector
          weights(1:nn) = BFweightsarray(1:nn,isg,iCL)
          STEMimages(1,i,j,iCL) = STEMimages(1,i,j,iCL) + sum(inten*weights)
! HAADF detector
          weights(1:nn) = ADFweightsarray(1:nn,isg,iCL)
          STEMimages(2,i,j,iCL) = STEMimages(2,i,j,iCL) + sum(inten*weights)
        end do
      end if
    end do donpiy
  end do donpix
!$OMP END DO 
call memth%dealloc(Azz, 'Azz', TID=TID)
call memth%dealloc(amp, 'amp', TID=TID)
! call memth%dealloc(amp2, 'amp2', TID=TID)
call memth%dealloc(weights, 'weights', TID=TID)
call memth%dealloc(inten, 'inten', TID=TID)

!$OMP barrier
!$OMP END PARALLEL

  if ((float(isg)/float(numstop) .gt. frac).and.(TID.eq.0)) then
!     call Time_remaining(isg,numstop)
     frac = frac + 0.05
  end if  

end do mainloop

! write the output to the still open HDF5 file
dataset = 'STEMimages'
if (self%nml%progmode.eq.'STEM') then
  hdferr = HDF%writeDatasetFloatArray(dataset, STEMimages, npix,npiy,nn,STEM%getnumberofsvalues() )
else
  hdferr = HDF%writeDatasetFloatArray(dataset, STEMimages, 2,npix,npiy,STEM%getnumCL() )
end if

call HDF%popall()

call Message%printMessage(' All data stored in file'//trim(self%nml%outname) )

end associate 

call closeFortranHDFInterface()

end subroutine SRdefect_



end module mod_SRdefect