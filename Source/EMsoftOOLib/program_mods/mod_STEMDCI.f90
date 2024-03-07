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

module mod_STEMDCI
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/14/24
  !!
  !! class definition for the EMSTEMDCI program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMSTEMDCI program
type, public :: STEMDCINameListType
  integer(kind=irg) :: DF_npix
  integer(kind=irg) :: DF_npiy
  integer(kind=irg) :: dinfo
  integer(kind=irg) :: nthreads
  integer(kind=irg) :: t_interval
  integer(kind=irg) :: kk(3)
  real(kind=sgl)    :: voltage
  real(kind=sgl)    :: lauec(2)
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
end type STEMDCINameListType

! class definition
type, public :: STEMDCI_T
private 
  character(fnlen)            :: nmldeffile = 'EMSTEMDCI.nml'
  type(STEMDCINameListType)   :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: STEMDCI_
  procedure, pass(self) :: setDF_npix_
  procedure, pass(self) :: getDF_npix_
  procedure, pass(self) :: setDF_npiy_
  procedure, pass(self) :: getDF_npiy_
  procedure, pass(self) :: setdinfo_
  procedure, pass(self) :: getdinfo_
  procedure, pass(self) :: setnthreads_
  procedure, pass(self) :: getnthreads_
  procedure, pass(self) :: sett_interval_
  procedure, pass(self) :: gett_interval_
  procedure, pass(self) :: setkk_
  procedure, pass(self) :: getkk_
  procedure, pass(self) :: setvoltage_
  procedure, pass(self) :: getvoltage_
  procedure, pass(self) :: setlauec_
  procedure, pass(self) :: getlauec_
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

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: STEMDCI => STEMDCI_
  generic, public :: setDF_npix => setDF_npix_
  generic, public :: getDF_npix => getDF_npix_
  generic, public :: setDF_npiy => setDF_npiy_
  generic, public :: getDF_npiy => getDF_npiy_
  generic, public :: setdinfo => setdinfo_
  generic, public :: getdinfo => getdinfo_
  generic, public :: setnthreads => setnthreads_
  generic, public :: getnthreads => getnthreads_
  generic, public :: sett_interval => sett_interval_
  generic, public :: gett_interval => gett_interval_
  generic, public :: setkk => setkk_
  generic, public :: getkk => getkk_
  generic, public :: setvoltage => setvoltage_
  generic, public :: getvoltage => getvoltage_
  generic, public :: setlauec => setlauec_
  generic, public :: getlauec => getlauec_
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
end type STEMDCI_T

! the constructor routine for this class 
interface STEMDCI_T
  module procedure STEMDCI_constructor
end interface STEMDCI_T

contains

!--------------------------------------------------------------------------
type(STEMDCI_T) function STEMDCI_constructor( nmlfile ) result(STEMDCI)
!! author: MDG 
!! version: 1.0 
!! date: 02/14/24
!!
!! constructor for the STEMDCI_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call STEMDCI%readNameList(nmlfile)

end function STEMDCI_constructor

!--------------------------------------------------------------------------
subroutine STEMDCI_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 02/14/24
!!
!! destructor for the STEMDCI_T Class
 
IMPLICIT NONE

type(STEMDCI_T), INTENT(INOUT)  :: self 

call reportDestructor('STEMDCI_T')

end subroutine STEMDCI_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/14/24
!!
!! read the namelist from an nml file for the STEMDCI_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(STEMDCI_T), INTENT(INOUT)          :: self
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
integer(kind=irg) :: nthreads
integer(kind=irg) :: t_interval
integer(kind=irg) :: kk(3)
real(kind=sgl)    :: voltage
real(kind=sgl)    :: lauec(2)
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
namelist / STEMDCIdata / nthreads, voltage, progmode, xtalname, kk, lauec, STEMnmlfile, &
                         outname, defectjsonfile, dispmode, dispfile, dinfo, t_interval, DF_L, &
                         DF_npix, DF_npiy, DF_slice, dmin

! SET INPUT PARAMETERS TO DEFAULT VALUES (EXCEPT XTALNAME, WHICH MUST BE PRESENT)
nthreads = 6
voltage = 200
progmode = 'STEM'
xtalname = 'undefined'
kk = (/ 0, 0, 1 /)
lauec = (/ 0.0, 0.0 /)
STEMnmlfile = 'undefined'
defectjsonfile = 'undefined'
outname = 'undefined'
dispfile = 'undefined'
dispmode = 'not'
dinfo = 0 ! 1 is verbose
t_interval = 5 ! update every x steps
DF_L = 1.0
DF_npix = 256
DF_npiy = 256
DF_slice = 1.0 ! slice thickness for scattering matrix approach (nmu)
dmin = 0.03

if (.not.skipread) then
! read the namelist file
  open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
  read(UNIT=dataunit,NML=STEMDCIdata)
  close(UNIT=dataunit,STATUS='keep')

! check for required entries
  if (trim(xtalname).eq.'undefined') then
    call Message%printError('STEMDCI:',' crystal structure file name is undefined in '//nmlfile)
  end if

  if (trim(defectjsonfile).eq.'undefined') then
    call Message%printError('STEMDCI:',' defectjsonfile file name is undefined in '//nmlfile)
  end if

  if (trim(outname).eq.'undefined') then
    call Message%printError('STEMDCI:',' outname file name is undefined in '//nmlfile)
  end if

  if ( (trim(dispfile).eq.'undefined').and.(dispmode.ne.'not') ) then
    call Message%printError('STEMDCI:',' dispfile file name is undefined in '//nmlfile)
  end if

  if (trim(STEMnmlfile).eq.'undefined') then
    call Message%printError('STEMDCI:',' STEMnmlfilefile name is undefined in '//nmlfile)
  end if
end if

self%nml%nthreads = nthreads
self%nml%voltage = voltage
self%nml%progmode = progmode
self%nml%xtalname = xtalname
self%nml%kk = kk
self%nml%lauec = lauec
self%nml%STEMnmlfile = STEMnmlfile
self%nml%outname = outname
self%nml%defectjsonfile = defectjsonfile
self%nml%dispmode = dispmode
self%nml%dispfile = dispfile
self%nml%dinfo = dinfo
self%nml%t_interval = t_interval
self%nml%DF_L = DF_L
self%nml%DF_npix = DF_npix
self%nml%DF_npiy = DF_npiy
self%nml%DF_slice = DF_slice
self%nml%dmin = dmin

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/14/24
!!
!! pass the namelist for the STEMDCI_T Class to the calling program

IMPLICIT NONE 

class(STEMDCI_T), INTENT(INOUT)          :: self
type(STEMDCINameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/14/24
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)         :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 5, n_real = 4
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( enl => self%nml )


! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! write all the single integers
io_int = (/ enl%DF_npix, enl%DF_npiy, enl%dinfo, enl%t_interval, enl%nthreads/)
intlist(1) = 'DF_npix'
intlist(2) = 'DF_npiy'
intlist(3) = 'dinfo'
intlist(4) = 't_interval'
intlist(5) = 'nthreads'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ enl%voltage, enl%DF_L, enl%DF_slice, enl%dmin /)
reallist(1) = 'voltage'
reallist(2) = 'DF_L'
reallist(3) = 'DF_slice'
reallist(4) = 'dmin'
call HDF%writeNMLreals(io_real, reallist, n_real)

! a 3-vector
dataset = 'kk'
hdferr = HDF%writeDatasetIntegerArray(dataset, enl%kk, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create kk dataset', hdferr)

! a 2-vector
dataset = 'lauec'
hdferr = HDF%writeDatasetFloatArray(dataset, enl%lauec, 2)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create lauec dataset', hdferr)

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
!! date: 02/14/24
!!
!! set DF_npix in the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%DF_npix = inp

end subroutine setDF_npix_

!--------------------------------------------------------------------------
function getDF_npix_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getDF_npix_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! get DF_npix from the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%DF_npix

end function getDF_npix_

!--------------------------------------------------------------------------
subroutine setDF_npiy_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setDF_npiy_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! set DF_npiy in the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%DF_npiy = inp

end subroutine setDF_npiy_

!--------------------------------------------------------------------------
function getDF_npiy_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getDF_npiy_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! get DF_npiy from the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%DF_npiy

end function getDF_npiy_

!--------------------------------------------------------------------------
subroutine setdinfo_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdinfo_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! set dinfo in the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%dinfo = inp

end subroutine setdinfo_

!--------------------------------------------------------------------------
function getdinfo_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdinfo_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! get dinfo from the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%dinfo

end function getdinfo_

!--------------------------------------------------------------------------
subroutine setnthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnthreads_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! set nthreads in the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%nthreads = inp

end subroutine setnthreads_

!--------------------------------------------------------------------------
function getnthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnthreads_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! get nthreads from the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%nthreads

end function getnthreads_

!--------------------------------------------------------------------------
subroutine sett_interval_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: sett_interval_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! set t_interval in the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%t_interval = inp

end subroutine sett_interval_

!--------------------------------------------------------------------------
function gett_interval_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: gett_interval_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! get t_interval from the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%t_interval

end function gett_interval_

!--------------------------------------------------------------------------
subroutine setkk_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setkk_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! set kk in the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp(3)

self%nml%kk = inp

end subroutine setkk_

!--------------------------------------------------------------------------
function getkk_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getkk_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! get kk from the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out(3)

out = self%nml%kk

end function getkk_

!--------------------------------------------------------------------------
subroutine setvoltage_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setvoltage_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! set voltage in the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%voltage = inp

end subroutine setvoltage_

!--------------------------------------------------------------------------
function getvoltage_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getvoltage_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! get voltage from the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%voltage

end function getvoltage_

!--------------------------------------------------------------------------
subroutine setlauec_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setlauec_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! set lauec in the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp(2)

self%nml%lauec = inp

end subroutine setlauec_

!--------------------------------------------------------------------------
function getlauec_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getlauec_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! get lauec from the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out(2)

out = self%nml%lauec

end function getlauec_

!--------------------------------------------------------------------------
subroutine setDF_L_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setDF_L_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! set DF_L in the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%DF_L = inp

end subroutine setDF_L_

!--------------------------------------------------------------------------
function getDF_L_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getDF_L_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! get DF_L from the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%DF_L

end function getDF_L_

!--------------------------------------------------------------------------
subroutine setDF_slice_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setDF_slice_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! set DF_slice in the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%DF_slice = inp

end subroutine setDF_slice_

!--------------------------------------------------------------------------
function getDF_slice_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getDF_slice_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! get DF_slice from the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%DF_slice

end function getDF_slice_

!--------------------------------------------------------------------------
subroutine setdmin_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdmin_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! set dmin in the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%dmin = inp

end subroutine setdmin_

!--------------------------------------------------------------------------
function getdmin_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdmin_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! get dmin from the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%dmin

end function getdmin_

!--------------------------------------------------------------------------
subroutine setprogmode_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setprogmode_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! set progmode in the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
character(4), INTENT(IN)       :: inp

self%nml%progmode = trim(inp)

end subroutine setprogmode_

!--------------------------------------------------------------------------
function getprogmode_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getprogmode_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! get progmode from the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
character(4)                   :: out

out = trim(self%nml%progmode)

end function getprogmode_

!--------------------------------------------------------------------------
subroutine setdispmode_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdispmode_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! set dispmode in the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
character(3), INTENT(IN)       :: inp

self%nml%dispmode = trim(inp)

end subroutine setdispmode_

!--------------------------------------------------------------------------
function getdispmode_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdispmode_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! get dispmode from the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
character(3)                   :: out

out = trim(self%nml%dispmode)

end function getdispmode_

!--------------------------------------------------------------------------
subroutine setoutname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setoutname_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! set outname in the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%outname = trim(inp)

end subroutine setoutname_

!--------------------------------------------------------------------------
function getoutname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getoutname_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! get outname from the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%outname)

end function getoutname_

!--------------------------------------------------------------------------
subroutine setdispfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdispfile_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! set dispfile in the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%dispfile = trim(inp)

end subroutine setdispfile_

!--------------------------------------------------------------------------
function getdispfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdispfile_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! get dispfile from the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%dispfile)

end function getdispfile_

!--------------------------------------------------------------------------
subroutine setxtalname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setxtalname_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! set xtalname in the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%xtalname = trim(inp)

end subroutine setxtalname_

!--------------------------------------------------------------------------
function getxtalname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getxtalname_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! get xtalname from the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%xtalname)

end function getxtalname_

!--------------------------------------------------------------------------
subroutine setSTEMnmlfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setSTEMnmlfile_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! set STEMnmlfile in the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%STEMnmlfile = trim(inp)

end subroutine setSTEMnmlfile_

!--------------------------------------------------------------------------
function getSTEMnmlfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getSTEMnmlfile_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! get STEMnmlfile from the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%STEMnmlfile)

end function getSTEMnmlfile_

!--------------------------------------------------------------------------
subroutine setdefectjsonfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdefectjsonfile_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! set defectjsonfile in the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%defectjsonfile = trim(inp)

end subroutine setdefectjsonfile_

!--------------------------------------------------------------------------
function getdefectjsonfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdefectjsonfile_
!! author: MDG
!! version: 1.0
!! date: 02/14/24
!!
!! get defectjsonfile from the STEMDCI_T class

IMPLICIT NONE

class(STEMDCI_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%defectjsonfile)

end function getdefectjsonfile_

!--------------------------------------------------------------------------
subroutine STEMDCI_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: STEMDCI_
!! author: MDG 
!! version: 1.0 
!! date: 02/14/24
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

class(STEMDCI_T), INTENT(INOUT)  :: self
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
type(reflisttype),pointer        :: reflist,firstw, rltmpa, rltmpb
type(kvectorlist),pointer        :: khead, ktmp
type(STEMGeometryNameListType)   :: stemnl

integer(kind=irg)                :: nn,i,j,k,npix,npiy,ii,jj,numvoids,numdisl, numYdisl, &
                                    numsf,nCL,numinc,dinfo,t_interval, DF_nums, dgn, pgnum, &
                                    DF_nums_new,DF_npix_new,DF_npiy_new, numstart,numstop, isg, TID, &
                                    NTHR, isym, ir, ga(3), gb(3),ic,g,numd,ix,iy, &
                                    numk,ixp,iyp,SETNTHR, io_int(6), skip, gg(3), iSTEM, hdferr, tick, tock
!                                OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
integer(kind=irg),parameter      :: numdd=180
real(kind=sgl)                   :: glen,exer,arg,thick, X(2), dmin, kk(3), &
                                    lauec(2), g3(3), gdotR,att,xgp,DF_gf(3), &
                                    DM(2,2), DD, H,FNr(3),ll(3),lpg(3),gplen,LC3, c(3), gx(3), gy(3), &
                                    sgdenom, gac(3), gbc(3),zmax, beamdiv, ktmax, io_real(2), kt, qx, qy
character(fnlen)                 :: dataname,sgname,dispfile,xtalname,foilnmlfile, STEMnmlfile
character(4)                     :: dispmode, progmode
character(11)                    :: dstr
character(15)                    :: tstrb
character(15)                    :: tstre
complex(kind=dbl),allocatable    :: DHWM(:,:),DHWMvoid(:,:),DDD(:,:),Sarray(:,:,:,:)
complex(kind=dbl),allocatable    :: amp(:),amp2(:),Azz(:,:)
complex(kind=dbl)                :: czero,cone
complex(kind=dbl)                :: para(0:numdd),dx,dy,dxm,dym
real(kind=sgl),allocatable       :: inten(:), sgarray(:,:), qxy(:,:)
real(kind=sgl),allocatable       :: disparray(:,:,:,:),imatvals(:,:), ZAimages(:,:,:,:)
integer(kind=irg),allocatable    :: BFweightsarray(:,:,:),ADFweightsarray(:,:,:), hkl(:,:), ggg(:,:)
integer(kind=sgl),allocatable    :: expval(:,:,:)
real(kind=sgl),allocatable       :: karray(:,:)
integer(kind=irg),allocatable    :: kij(:,:)
complex(kind=dbl),allocatable    :: DynMat(:,:)
character(fnlen)                 :: dataset, instring, datafile, datagroupname, groupname
integer(HSIZE_T), dimension(1:4) :: hdims, offset 
integer(HSIZE_T)                 :: dims4(4)
logical                          :: usehex, verbose, f_exists, insert=.TRUE., overwrite=.TRUE.
integer(kind=irg)                :: error_cnt, ier,numsval, nbeams
complex(kind=dbl),allocatable    :: DF_Svoid(:,:), DHWMz(:,:)!, DF_Sarray(:,:,:,:)
integer(kind=irg)                :: DF_npix,DF_npiy
real(kind=dbl)                   :: DynFN(3), xx, DF_slice, lambda
real(kind=sgl)                   :: ijmax, tstart, tstop
integer(kind=irg)                :: intijmax
character(fnlen)                 :: outname
character(fnlen,kind=c_char)     :: line2(1)

! initialize the STEM geometry namelist as well 
datafile = EMsoft%generateFilePath('EMdatapathname',trim(self%getSTEMnmlfile()))
STEM = STEM_T( datafile )
stemnl = STEM%getNameList()

call openFortranHDFInterface()
HDF = HDF_T()

associate( enl => self%nml )

! initialize the timing routines
timer = Timing_T()
tstrb = timer%getTimeString()
call timer%Time_tick()

! first we define some default values
czero=cmplx(0.0,0.0,dbl)
cone=cmplx(1.0,0.0,dbl)
numd = numdd
usehex = .FALSE.

! convert some parameters from namelist to local variables. 
DF_npix = self%getDF_npix()
DF_npiy = self%getDF_npiy()
npix = DF_npix
npiy = DF_npiy
DF_slice = self%getDF_slice()
dmin = self%getdmin()
DynFN = float(enl%kk)

! set the defect parameters
defects%DF_slice = self%getDF_slice()
defects%DF_npix = DF_npix
defects%DF_npiy = DF_npiy
defects%DF_L = self%getDF_L()

!=============================================
! crystallography section
verbose = .TRUE.

call cell%setFileName(enl%xtalname)
call Diff%setrlpmethod('WK')
call Diff%setV(dble(enl%voltage))

call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, enl%dmin, verbose, useHDF=HDF)
lambda = Diff%getWaveLength()

! determine the point group number
j=0
do i=1,32
if (SGPG(i).le.SG%getSpaceGroupNumber()) j=i
end do

dgn = SG%GetPatternSymmetry(enl%kk,j,.TRUE.)
pgnum = j
isym = WPPG(dgn) ! WPPG lists the whole pattern point group numbers vs. diffraction group numbers

! determine the shortest reciprocal lattice points for this zone
call cell%ShortestG(SG, enl%kk, ga, gb, isym)
io_int(1:3) = enl%kk(1:3)
call Message%WriteValue('', io_int, 3,  "(//,' ','[',3I2,'] has Bright Field symmetry ')",advance="no")
call Message%printMessage(PGTWD(isym),"(' ',A,', ')",advance="no")
io_int(1) = ir
call Message%WriteValue(' order = ', io_int, 1, "(I4/)")
io_int(1:3)=ga(1:3)
io_int(4:6)=gb(1:3)
call Message%WriteValue(' Reciprocal lattice vectors : ', io_int, 6,"('(',3I3,') and (',3I3,')',/)")
call Message%printMessage('  (the first lattice vector is horizontal in the CBED pattern)')
call Message%printMessage(' ')

! determine the cartesian components of ga
defects%DF_gf = float(ga)
defects%DF_gstar = defects%DF_gf/cell%CalcLength(defects%DF_gf,'r')**2    ! define G* such that G.G* = 1
call cell%TransSpace(defects%DF_gf,defects%DF_gc,'r','c')                 ! convert to Cartesian reference frame

! read foil and defect data 
if (self%getdinfo_().eq.0) then 
  call defects%setdinfo(.FALSE.)
else
  call defects%setdinfo(.TRUE.)
end if
call defects%InitializeDefects(EMsoft,cell,self%getdefectjsonfile(),npix,npiy, &
                               self%getDF_L(),DF_gf,dble(enl%kk),dble(ga),error_cnt)

! set the Bethe parameters
call Diff%SetBetheParameters(EMsoft, .TRUE.)

kvec = kvectors_T()
gvec = gvectors_T()

call kvec%set_mapmode('Conical')
call kvec%set_kinp( dble(enl%kk) )

if (enl%progmode.ne.'EM') then
  call Message%printMessage('Progmode = '//enl%progmode)
  beamdiv = stemnl%beamconvergence
  kt = stemnl%kt
  numsval = stemnl%numberofsvalues
  ktmax = 2.0*sin(beamdiv/2000.)/lambda/cell%CalcLength(float(ga),'r')     ! ktmax in units of |ga|
  call kvec%set_ktmax( dble(ktmax) )
  ijmax = float(numsval)**2
  intijmax=int(ijmax)
  call kvec%Calckvectors(cell,SG,Diff,dble(ga),numsval,numsval,intijmax,usehex)
  ! here we figure out how many beams there are
  xx = 1.0/cell%CalcLength(float(enl%kk),'d')
  kk = xx*float(enl%kk)/lambda  
  call gvec%Initialize_ReflectionList(cell, SG, Diff, sngl(DynFN), kk, enl%dmin, verbose) 
else  ! progmode = EM
  ijmax = 0.0
  intijmax=int(ijmax)
  call kvec%Calckvectors(cell,SG,Diff,dble(ga),numsval,numsval,intijmax,usehex)
  ! here we figure out how many beams there are
  xx = 1.0/cell%CalcLength(float(enl%kk),'d')
  kk = xx*float(enl%kk)/lambda  
  call gvec%Initialize_ReflectionList(cell, SG, Diff, sngl(DynFN), kk, enl%dmin, verbose)
end if
call Message%printMessage(' --> done initializing reflectionlist')

call STEM%setnumk(kvec%get_numk())
numk = kvec%get_numk()
nn = gvec%get_nref()
khead => kvec%get_ListHead()
reflist => gvec%get_ListHead()

if (enl%progmode.ne.'EM') then
  call STEM%init_STEM_ZA(stemnl, cell, DynFN, Diff, khead, reflist, nn)
end if

! allocate and initialize various arrays
mem = memory_T()
call mem%alloc(DF_Svoid, (/ nn,nn /), 'DF_Svoid', initval=czero )
call mem%alloc(DHWMz, (/ nn,nn /), 'DF_SvDHWMzoid', initval=czero )
call mem%alloc(DHWM, (/ nn,nn /), 'DHWM', initval=czero )
call mem%alloc(DHWMvoid, (/ nn,nn /), 'DHWMvoid', initval=czero )

! Compute the off-diagonal part of the complex DHW matrix (factor i is included)
! We can precompute those because they will not change at all during the run
!       (these lines implement the equations on page 476 of the EM book)
! In this program, the reflections are stored using linked lists, which does not lend itself to
! OpenMP acceleration;  so, we'll have to re-write this at some point in the future...
!
! this is also where we compute the decomposition of the reflection indices w.r.t. ga and gb,
! and store them in the nab(1:2) field of the linked list; these are used to compute the 
! defect contributions to the dynamical matrix in the displace routine of the MEmath module
!
DM(1,1) = cell%CalcDot(float(gb),float(gb),'c')
DM(1,2) = -cell%CalcDot(float(ga),float(gb),'c')
DM(2,1) = DM(1,2)
DM(2,2) = cell%CalcDot(float(ga),float(ga),'c')
DD = DM(1,1)*DM(2,2) - DM(1,2)*DM(2,1)

rltmpa => reflist%next    ! point to the front of the list
! ir is the row index
do ir=1,nn
 rltmpb => reflist%next   ! point to the front of the list
! ic is the column index
 do ic=1,nn
  if (ic.ne.ir) then  ! exclude the diagonal
! compute Fourier coefficient of electrostatic lattice potential 
   call Diff%CalcUcg(cell,rltmpa%hkl - rltmpb%hkl)
   rlp = Diff%getrlp()
   DHWMz(ir,ic) = cPi*cmplx(-aimag(rlp%qg),real(rlp%qg),dbl)  ! and initialize the off-diagonal matrix element (including i)
  end if
  rltmpb => rltmpb%next  ! move to next column-entry
 end do
! decompose this point w.r.t ga and gb
 X(1) = cell%CalcDot(float(rltmpa%hkl),float(ga),'c')
 X(2) = cell%CalcDot(float(rltmpa%hkl),float(gb),'c')
 X = matmul(DM,X)/DD
 rltmpa%nab(1:2) = int(X(1:2))
 rltmpa => rltmpa%next   ! move to next row-entry
end do
call Message%printMessage(' Reference Darwin-Howie-Whelan matrix initialized')

! compute the normal absorption factor xgp
call Diff%CalcUcg(cell,(/0,0,0/))
rlp = Diff%getrlp()
xgp = aimag(rlp%qg)
io_real(1) = 1.0/xgp
call Message%WriteValue('Normal absorption length : ', io_real, 1, "(F10.5/)")

! Set up the excitation errors for EM illumination mode;
! distance between consecutive HOLZ layers in nm-1
if (progmode.eq.'EM') then  
  H = 1.0/cell%CalcLength(float(enl%kk),'d')
! g3 basis vector, properly scaled
  call cell%CalcCross(float(ga),float(gb),g3,'r','r',1)
  call cell%NormVec(g3,'r')
  g3 = H * g3
! unit foil normal in reciprocal space  
  call cell%TransSpace(sngl(defects%foil%F),FNr,'d','r')
  call cell%NormVec(FNr,'r')
 
! parallel illumination, so the excitation errors need to be computed only once
! fill the diagonal of the reference dynamical matrix and the void matrix
  rltmpa => reflist%next   ! point to the front of the list
! ir is the row index
  do ir=1,nn
   glen = cell%CalcLength(float(rltmpa%hkl),'r')
   if (glen.eq.0.0) then
    DHWMz(ir,ir) = czero
   else  ! compute the excitation error
    ll = lauec(1)*ga + lauec(2)*gb   ! Laue center vector
    lpg = ll + rltmpa%hkl                ! Laue + g
    gplen = cell%CalcLength(lpg,'r')
    LC3 = sqrt(1.0-lambda**2*(cell%CalcLength(ll,'r')**2))   ! to ensure proper normalization of wave vector
    if (gplen.eq.0.0) then
      exer=-lambda*cell%CalcDot(float(rltmpa%hkl),ll+lpg,'r')/2.0*LC3*cos(cell%CalcAngle(dble(enl%kk),defects%foil%F,'d'))        
    else
      sgdenom=2.0*LC3*cos(cell%CalcAngle(dble(enl%kk),defects%foil%F,'d'))-2.0*lambda*gplen*cos(cell%CalcAngle(lpg,FNr,'r'))
      exer=-(lambda*cell%CalcDot(float(rltmpa%hkl),ll+lpg,'r')-2.0*LC3*cell%CalcDot(g3,lpg,'r'))/sgdenom
    end if
    rltmpa%sg = exer
    DHWMz(ir,ir) = cmplx(0.0,2.D0*cPi*exer,dbl)
    ! call CalcUcg(cell,rlp,rltmpa%hkl)
   end if 
   DHWMvoid(ir,ir) = DHWMz(ir,ir)
   rltmpa => rltmpa%next   ! move to next row-entry
  end do
end if

! define the foil thickness, attenuation, and number slices per column
thick = defects%foil%zb    ! this is the same everywhere for this version; needs to be updated in the next version
att = exp(-2.0*cPi*thick*xgp)  ! this is the global intensity attenuation factor; remove the factor 2.0 if amplitudes are needed
defects%DF_nums = nint(thick/enl%DF_slice)  ! this is the number of slices for this particular column
DF_nums = defects%DF_nums


! setup the defect displacement field parameters
if ((self%nml%dispmode.eq.'new').or.(self%nml%dispmode.eq.'not')) then
  call Message%printMessage(' Starting Displacement Field Computation')

! precompute ALL the defect columns and, if needed, store them in dispfile
! this portion should be carried out in multi-threaded mode as much as possible
  call mem%alloc(disparray, (/ 2, DF_nums, DF_npix, DF_npiy /), 'disparray', initval=0.0 )
  call mem%alloc(imatvals, (/ 2, DF_nums /), 'imatvals', initval=0.0  )

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
  ! call mem%alloc(DF_R, (/ DF_nums,3 /), 'DF_R')     ! each thread has its own DF_R array
  call mem%alloc(defects%DF_R, (/ DF_nums,3 /), 'defects%DF_R')     ! each thread has its own DF_R array  

 call cell%TransSpace(float(ga),gac,'r','c')
 call cell%TransSpace(float(gb),gbc,'r','c')

!!$OMP DO SCHEDULE (GUIDED)
  ! write(*,*) TID,': starting Do Schedule'
  do i=1,DF_npix  
    do j=1,DF_npiy
      defects%DF_R = 0.0
! compute the displacement vectors DF_R for all points in the column
      call defects%CalcR(cell,i,j)
! loop over the fixed thickness slices
      do k=1,DF_nums
! then convert to the dot-product 
        if (defects%DF_R(k,1).eq.-10000.0) then  ! this is point inside a void
          imatvals(1:2,k) = -10000
        else  ! it is not a void, so use the full dot product g.R (all vectors must be Cartesian !)
! use gac and gbc to get the two dot products and store both of them as integers mapped onto the 0..numd range
          gdotR = Dot_Product(gac,defects%DF_R(k,1:3))
          imatvals(1,k) = numd*amod(gdotR+1000.0,1.0)
          gdotR = Dot_Product(gbc,defects%DF_R(k,1:3))
          imatvals(2,k) = numd*amod(gdotR+1000.0,1.0)
        end if
        disparray(1:2,1:defects%DF_nums,i,j) = imatvals(1:2,1:defects%DF_nums)
     end do ! k loop
   end do
!  if ((mod(i,t_interval).eq.0).and.(TID.eq.0)) call Time_remaining(TT,i,defects%DF_npix)
end do
!!$OMP END DO 
! if (TID.eq.0) call Time_stop(TT,defects%DF_npix*defects%DF_npiy)
!!$OMP END PARALLEL
end if
call Message%printMessage('   -> done.')

! and, if needed, store the defect displacement field for re-runs
if (self%nml%dispmode.ne.'not') then
  if ((self%nml%dispfile.ne.'none').and.(self%nml%dispmode.eq.'new')) then 
    dataname = EMsoft%generateFilePath('EMdatapathname',trim(self%nml%dispfile))
    hdferr = HDF%createFile(dataname)

    dataset = 'dispfield'
      hdferr = HDF%writeDatasetFloatArray(dataset, disparray, 2, DF_nums, DF_npix, DF_npiy )

    call HDF%popall()
    call Message%printMessage('Displacement field data stored in file '//self%nml%dispfile)
  endif
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
end if

! ok, all the set up is now complete;! next, we prepare for the actual image simulation
! there are three types of simulation: EM, BF/HAADF, or STEM with storage of CBED patterns
if (trim(enl%progmode).eq.'EM') call mem%alloc(ZAimages, (/ nn,npix,npiy,1 /), 'ZAimages', initval=0.0)
if (enl%progmode.eq.'BFDF') call mem%alloc(ZAimages, (/ 2,npix,npiy,stemnl%numCL /), 'ZAimages', initval=0.0)
! for STEM mode; we'll store data in blocks of 10 k-vectors
if (enl%progmode.eq.'STEM') call mem%alloc(ZAimages, (/ npix,npiy,nn,10 /), 'ZAimages', initval=0.0)

! copy the weights arrays since they are private to the STEM class
if (enl%progmode.ne.'EM') then
  BFweightsarray = STEM%getZABFweightsarray()
  ADFweightsarray = STEM%getZAADFweightsarray()
  sgarray = STEM%getsgarray()
end if

! loop over all reflections to get the appropriate powers
call mem%alloc(expval, (/ 2,nn,nn /), 'expval')
rltmpa => reflist%next    ! point to the front of the list
! ir is the row index
do ir=1,nn
   rltmpb => reflist%next   ! point to the front of the list
! ic is the column index
   do ic=1,nn
     if (ic.ne.ir) then  ! exclude the diagonal
       expval(1,ir,ic) = rltmpa%nab(1)-rltmpb%nab(1) 
       expval(2,ir,ic) = rltmpa%nab(2)-rltmpb%nab(2)
     end if
     rltmpb => rltmpb%next  ! move to next column-entry
  end do
  rltmpa => rltmpa%next   ! move to next row-entry
end do
call Message%printMessage(' initialized expval')

if (trim(enl%progmode).eq.'EM') then
  numstart = 1
  numstop = 1
else
  numstart = 1
  numstop = kvec%get_numk()
  io_int(1) = numstop
  call Message%WriteValue(' STEM number of beam directions =  ', io_int, 1, "(I5)")
end if

nCL = stemnl%numCL     ! set the number of camera length values

! define the numd complex defect parameters
do i=0,numd
  arg = 2.D0*cPi*float(i)/dble(numd)
  para(i) = cmplx(cos(arg),sin(arg),dbl)
end do

! before we start the main computational loop, we need to store the 
! necessary variables for the reconstruction of STEM BF/HAADF images
! using the post-processing IDL routine.  This file should contain all
! the information needed to recreate the full CBED patterns at each image
! pixel.  This is then followed by the actual data.

! the file format is identical to that of the EMSTEMDCI program in STEM mode,
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

dataset = 'STEMGeometry'
hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%generateFilePath('EMdatapathname', self%nml%STEMnmlfile) )

! we also need to include the defect JSON file here in this group 
dataset = 'defectJSONfile'
hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%generateFilePath('EMdatapathname', self%nml%defectjsonfile) )

! and the foil JSON file
dataset = 'foilJSONfile'
hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%generateFilePath('EMdatapathname', defects%foilname) )

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
dataset = SC_BraggAngle
  hdferr = HDF%writeDatasetFloat(dataset, Diff%CalcDiffAngle(cell, ga)*0.5)

dataset = SC_wavelength
  hdferr = HDF%writeDatasetFloat(dataset, sngl(lambda) )

dataset = SC_pixelsize
  hdferr = HDF%writeDatasetFloat(dataset, defects%DF_L)

dataset = SC_numreflections
  hdferr = HDF%writeDatasetInteger(dataset, nn)

! number of reflections, and associated information (hkl, ...)
  call cell%TransSpace(float(enl%kk),c,'d','c')
  call cell%NormVec(c,'c')
! then make ga the x-axis
  call cell%TransSpace(float(ga),gx,'r','c')
  call cell%NormVec(gx,'c')
! compute the cross product between k and gx; this is the y-axis
  call cell%CalcCross(c,gx,gy,'c','c',0)
  
! create some temporary arrays
call mem%alloc(ggg, (/ 3, nn /), 'ggg', initval=0 )
call mem%alloc(qxy, (/ 2, nn /), 'qxy', initval=0.0 )
rltmpa => reflist%next
do ic = 1, nn
  ggg(1:3,ic) = rltmpa%hkl
  call cell%TransSpace(float(rltmpa%hkl),c,'r','c')
  qxy(1:2,ic) = (/ cell%CalcDot(c, gx, 'c'), cell%CalcDot(c, gy, 'c') /)
  rltmpa => rltmpa%next
end do

dataset = SC_hkl
  hdferr = HDF%writeDatasetIntegerArray(dataset, ggg, 3, nn )

dataset = SC_qxy
  hdferr = HDF%writeDatasetFloatArray(dataset, qxy, 2, nn )

call mem%dealloc(ggg, 'ggg')
call mem%dealloc(qxy, 'qxy')

dataset = SC_numk
  hdferr = HDF%writeDatasetInteger(dataset, numk)

! number of wave vectors, tangential components, etc...
call mem%alloc(kij, (/ 2,numk /), 'kij')
ktmp => khead
do ic=1,numk
  kij(1:2,ic) = (/ ktmp%i, ktmp%j /)
  ktmp => ktmp%next
end do
dataset = SC_kij
  hdferr = HDF%writeDatasetIntegerArray(dataset, kij, 2, numk)
call mem%dealloc(kij, 'kij')

! finally, we need to create the ZAimages output array in such a way that 
! we can hyperslab it for the STEM mode, and just write the whole thing 
! directly for the other modes.  Either way, we create the complete array
! here and then write to it when we need to...

dataset = SC_ZAimages
 offset = (/ 0, 0, 0, 0/)

! EM mode
 if (trim(enl%progmode).eq.'EM') then
   hdims = (/ nn, npix, npiy, 1 /)
   dims4 = (/ nn, npix, npiy, 1 /)
 end if

! BFDF mode
 if (enl%progmode.eq.'BFDF') then
   hdims = (/ 2, npix, npiy, stemnl%numCL /)
   dims4 = (/ 2, npix, npiy, stemnl%numCL /)
 end if

! STEM mode
 if (enl%progmode.eq.'STEM') then
   hdims = (/ npix, npiy, nn, numk /)
   dims4 = (/ npix, npiy, nn, 10 /)
 end if

 hdferr = HDF%writeHyperslabFloatArray(dataset, ZAimages, hdims, offset, dims4)
! we leave the HDF5 file open for further writing of data sets
call Message%printMessage(' Initialized output HDF file; starting main computational loop')

call mem%alloc(Sarray, (/ nn,nn,numd,numd/), 'Sarray', startdims=(/1,1,0,0/) )
iSTEM = 0

memth = Memory_T( nt = self%getnthreads() )
!--------------------------------------------------------------
!--------------------------------------------------------------
!--------------------------------------------------------------
mainloop: do isg = numstart,numstop   ! this is the main computational loop
 iSTEM = iSTEM+1
!--------------------------------------------------------------
! here we precompute an array of scattering matrices that can 
! then be used, either directly, or via bi-linear interpolation,
! by the image computation portion of this program.
!
! For starters, we'll subdivide the range of possible alpha values
! in 180 segments (2 degrees each), with a copy of the last one
! (i.e., 181x181 = 32761 entries or 240 Mb for 31 beams)
!
! this part is essentially the same as the threaded section of the older
! STEMDCI.all.f90 program (which was a test program).
!

! get the correct excitation errors for this beam orientation (in STEM mode);
! no need to do anything in EM mode
if (trim(enl%progmode).ne.'EM') then  
! fill the diagonal of the reference dynamical matrix and the void matrix
  do i=1,nn
   DHWMz(i,i)=2.0*cPi*cmplx(0.0,sgarray(i,isg))    ! initialize the diagonal elements of the dynamical matrix
   DHWMvoid(i,i) = DHWMz(i,i)
  end do
end if

NTHR = self%getnthreads()

call OMP_SET_NUM_THREADS(NTHR)

!$OMP  PARALLEL PRIVATE(TID,i,j,k,ii,jj,ic,ir,g,Azz,DDD,zmax)
TID = OMP_GET_THREAD_NUM() 

! these are private variables, so each thread must allocate them !
call memth%alloc(Azz, (/ nn,nn /), 'Azz', TID=TID)
call memth%alloc(DDD, (/ nn,nn /), 'DDD', TID=TID)   

!$OMP DO SCHEDULE(STATIC)
do j=0,numd
 do i=0,numd 
! loop over all reflections in the array DD using the information in expval
! ir is the row index
  do ir=1,nn
! ic is the column index
   do ic=1,nn
    if (ic.ne.ir) then  ! exclude the diagonal
     DDD(ir,ic) = DHWMz(ir,ic) * para(i)**expval(1,ir,ic) * para(j)**expval(2,ir,ic)
    else
     DDD(ir,ic) = DHWMz(ir,ic) 
    end if
   end do
  end do
    
  call MatrixExponential(DDD, Azz, dble(enl%DF_slice), 'Pade', nn)  
    
  Sarray(1:nn,1:nn,i,j) = Azz(1:nn,1:nn)
 end do
 if ((TID.eq.0).and.(mod(j,10).eq.0)) then
    call Message%printMessage('.',"(A1)",advance="no")
 end if
end do
!$OMP END DO

call memth%dealloc(Azz, 'Azz', TID=TID)
call memth%dealloc(DDD, 'DDD', TID=TID)
!$OMP END PARALLEL

call Message%printMessage(' 181x181 scattering matrices precomputed ',"(A,' ')",advance="no")

!----------------------------------------------------!
! Finally, here it is: the actual image computation  !
!----------------------------------------------------!
call OMP_SET_NUM_THREADS(NTHR)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID,i,j,k,ii,Azz,amp,amp2,ix,iy,dx,dy,dxm,dym,inten,ixp,iyp)
TID = OMP_GET_THREAD_NUM() 
  
  ! these are private variables, so each thread must allocate them !
call memth%alloc(Azz, (/ nn,nn /), 'Azz', TID=TID)
call memth%alloc(amp, (/ nn /), 'amp', TID=TID)   
call memth%alloc(amp2, (/ nn /), 'amp2', TID=TID)   
call memth%alloc(inten, (/ nn /), 'inten', TID=TID)   

!$OMP DO SCHEDULE (STATIC)
 donpix: do i=1,npix
 if ((TID.eq.0).and.(mod(i,10).eq.0)) then
   call Message%printMessage('.',"(A1)",advance="no")
 end if
 donpiy:   do j=1,npiy
! initialize the wave function for this pixel with (1.0,0.0) for the incident beam
    amp = czero
    amp(1) = cone
    doslices: do k=1,defects%DF_nums    ! loop over the fixed thickness slices
! compute the appropriate scattering matrix to propagate with (see section 8.3.3 in the book)
       if (disparray(1,k,i,j).eq.-10000) then  ! this is point inside a void
         Azz = DF_Svoid    ! so we use the void propagator matrix
       else  ! it is not a void
! in this version, we use complex bi-variate interpolation to select the appropriate Azz 
! matrix from the pre-computed Sarray.
         ix = int(disparray(1,k,i,j))
         ixp = ix+1
         if (ix.eq.numd) ixp=1
         iy = int(disparray(2,k,i,j))
         iyp = iy+1
         if (iy.eq.numd) iyp=1
         dx = cmplx(amod(disparray(1,k,i,j),1.0),0.0)
         dy = cmplx(amod(disparray(2,k,i,j),1.0),0.0)
         dxm = cone-dx
         dym = cone-dy
         Azz = dxm*dym*Sarray(1:nn,1:nn,ix,iy)+dx*dym*Sarray(1:nn,1:nn,ixp,iy)+ &
                    dxm*dy*Sarray(1:nn,1:nn,ix,iyp)+dx*dy*Sarray(1:nn,1:nn,ixp,iyp)
       end if
! and multiply with this matrix
       amp2 = matmul(Azz, amp)
       amp = amp2
    end do doslices ! loop over slices 
      
! compute the intensities and (for STEM mode) use the proper weight-factor
    inten(1:nn) = att*abs(amp(1:nn))**2

! and store the intensities in the ZAimages array
    if (enl%progmode.eq.'STEM') then
        ZAimages(i,j,1:nn,iSTEM) = inten(1:nn)
    end if
    if (enl%progmode.eq.'BFDF') then
        do ii=1,nn
         do jj=1,nCL
          if (BFweightsarray(ii,isg,jj).eq.1) ZAimages(1,i,j,jj) = ZAimages(1,i,j,jj) + inten(ii)
          if (ADFweightsarray(ii,isg,jj).eq.1) ZAimages(2,i,j,jj) = ZAimages(2,i,j,jj) + inten(ii)
         end do
        end do
      end if
    if (trim(enl%progmode).eq.'EM') then
      ZAimages(1:nn,i,j,1) = inten(1:nn)
    end if

    end do donpiy
end do donpix
!$OMP END DO
call memth%dealloc(Azz, 'Azz', TID=TID)
call memth%dealloc(amp, 'amp', TID=TID)
call memth%dealloc(amp2, 'amp2', TID=TID)
call memth%dealloc(inten, 'inten', TID=TID)
!$OMP END PARALLEL

io_int(1) = isg 
io_int(2) = numstop
call Message%WriteValue('... completed ', io_int, 2, "(I4,'/',I4)")
  
  if (enl%progmode.eq.'STEM') then
    if (mod(isg,10).eq.0) then
       call Message%printMessage('Storing block of 10 k-vector results')
       offset = (/ 0, 0, 0, isg - 10 /)
       hdims = (/ npix, npiy, nn, numk /)
       dims4 = (/ npix, npiy, nn, 10 /)
       hdferr = HDF%writeHyperslabFloatArray(dataset, ZAimages, hdims, offset, dims4, insert)
       iSTEM = 0
       ZAimages = 0.0
    end if
  end if

200 end do mainloop

call mem%dealloc(Sarray, 'Sarray')

! and write the final data arrays to the HDF5 file in the appropriate format
if (trim(enl%progmode).eq.'EM') then
   offset = (/ 0, 0, 0, 0/)
   hdims = (/ nn, npix, npiy, 1 /)
   dims4 = (/ nn, npix, npiy, 1 /)
end if 
if (enl%progmode.eq.'BFDF') then
   offset = (/ 0, 0, 0, 0/)
   hdims = (/ 2, npix, npiy, stemnl%numCL /)
   dims4 = (/ 2, npix, npiy, stemnl%numCL /)
end if 
! flush the final ZAimages array if necessary (depends on value of iSTEM)
if ((enl%progmode.eq.'STEM').AND.(iSTEM.ne.0)) then
   offset = (/ 0, 0, 0, numk - iSTEM /)
   hdims = (/ npix, npiy, nn, numk /)
   dims4 = (/ npix, npiy, nn, iSTEM /)
end if 

dataset = SC_ZAimages
   hdferr = HDF%writeHyperslabFloatArray(dataset, ZAimages, hdims, offset, dims4, insert)

! close the EMData group
call HDF%pop()
call HDF%pop()

! and update the end time
tstre = timer%getTimeString()
call timer%Time_tock()

groupname = SC_EMheader
  hdferr = HDF%openGroup(groupname)

datagroupname = "STEMDCI"
  hdferr = HDF%openGroup(datagroupname)

! stop time /EMheader/StopTime 'character'
dataset = SC_StopTime
line2(1) = dstr//', '//tstre
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)

dataset = SC_Duration
  hdferr = HDF%writeDatasetFloat(dataset, timer%getInterval())

! close the datafile
call HDF%popall()

call closeFortranHDFInterface

call Message%printMessage('Data stored in file '//trim(enl%outname),"(/A/)")

end associate

end subroutine STEMDCI_



end module mod_STEMDCI