! ###################################################################
! Copyright (c) 2013-2023, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_HREBSD
  !! author: MDG (OO version)/ Chaoyi Zhu (original version)
  !! version: 1.0 
  !! date: 06/28/23
  !!
  !! class definition for the EMHREBSD program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMHREBSD program
type, public :: HREBSDNameListType
   integer(kind=irg)    :: numsx
   integer(kind=irg)    :: numsy
   integer(kind=irg)    :: ipf_wd
   integer(kind=irg)    :: ipf_ht
   integer(kind=irg)    :: patx
   integer(kind=irg)    :: paty
   integer(kind=irg)    :: N_ROI
   integer(kind=irg)    :: size_ROI
   integer(kind=irg)    :: roi_distance
   integer(kind=irg)    :: nthreads
   real(kind=sgl)       :: step_size
   real(kind=sgl)       :: delta
   real(kind=sgl)       :: thetac
   real(kind=sgl)       :: highpass
   real(kind=sgl)       :: lowpass
   real(kind=sgl)       :: C11
   real(kind=sgl)       :: C12
   real(kind=sgl)       :: C44
   real(kind=sgl)       :: C13
   real(kind=sgl)       :: C33
   character(fnlen)     :: datafile
   character(fnlen)     :: masterfile
   character(fnlen)     :: exptfile
   character(fnlen)     :: inputtype
   character(fnlen)     :: HDFstrings(10)
   character(1)         :: Remap
   character(1)         :: PCrefine
   character(3)         :: crystal
end type HREBSDNameListType

! other type declarations
type EBSDSEMArray
        real(kind=sgl)                  :: Step_X
        real(kind=sgl)                  :: Step_Y
        integer(kind=irg),allocatable   :: SEM_X(:)
        integer(kind=irg),allocatable   :: SEM_Y(:)
        integer(kind=irg),allocatable   :: GrainID(:)
        real(kind=sgl),allocatable      :: GrainID_X(:)
        real(kind=sgl),allocatable      :: GrainID_Y(:)
end type

! class definition
type, public :: HREBSD_T
private 
  character(fnlen)          :: nmldeffile = 'EMHREBSD.nml'
  type(HREBSDNameListType)  :: nml 
  type(EBSDSEMArray)        :: SEM
  real(kind=sgl)            :: PC(3)

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: HREBSD_
  procedure, pass(self) :: setnumsx_
  procedure, pass(self) :: getnumsx_
  procedure, pass(self) :: setnumsy_
  procedure, pass(self) :: getnumsy_
  procedure, pass(self) :: setipf_wd_
  procedure, pass(self) :: getipf_wd_
  procedure, pass(self) :: setipf_ht_
  procedure, pass(self) :: getipf_ht_
  procedure, pass(self) :: setpatx_
  procedure, pass(self) :: getpatx_
  procedure, pass(self) :: setpaty_
  procedure, pass(self) :: getpaty_
  procedure, pass(self) :: setN_ROI_
  procedure, pass(self) :: getN_ROI_
  procedure, pass(self) :: setsize_ROI_
  procedure, pass(self) :: getsize_ROI_
  procedure, pass(self) :: setroi_distance_
  procedure, pass(self) :: getroi_distance_
  procedure, pass(self) :: setnthreads_
  procedure, pass(self) :: getnthreads_
  procedure, pass(self) :: setstep_size_
  procedure, pass(self) :: getstep_size_
  procedure, pass(self) :: setdelta_
  procedure, pass(self) :: getdelta_
  procedure, pass(self) :: setthetac_
  procedure, pass(self) :: getthetac_
  procedure, pass(self) :: sethighpass_
  procedure, pass(self) :: gethighpass_
  procedure, pass(self) :: setlowpass_
  procedure, pass(self) :: getlowpass_
  procedure, pass(self) :: setC11_
  procedure, pass(self) :: getC11_
  procedure, pass(self) :: setC12_
  procedure, pass(self) :: getC12_
  procedure, pass(self) :: setC44_
  procedure, pass(self) :: getC44_
  procedure, pass(self) :: setC13_
  procedure, pass(self) :: getC13_
  procedure, pass(self) :: setC33_
  procedure, pass(self) :: getC33_
  procedure, pass(self) :: setdatafile_
  procedure, pass(self) :: getdatafile_
  procedure, pass(self) :: setmasterfile_
  procedure, pass(self) :: getmasterfile_
  procedure, pass(self) :: setexptfile_
  procedure, pass(self) :: getexptfile_
  procedure, pass(self) :: setinputtype_
  procedure, pass(self) :: getinputtype_
  procedure, pass(self) :: setHDFstrings_
  procedure, pass(self) :: getHDFstrings_
  procedure, pass(self) :: setRemap_
  procedure, pass(self) :: getRemap_
  procedure, pass(self) :: setPCrefine_
  procedure, pass(self) :: getPCrefine_
  procedure, pass(self) :: setcrystal_
  procedure, pass(self) :: getcrystal_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: HREBSD => HREBSD_
  generic, public :: setnumsx => setnumsx_
  generic, public :: getnumsx => getnumsx_
  generic, public :: setnumsy => setnumsy_
  generic, public :: getnumsy => getnumsy_
  generic, public :: setipf_wd => setipf_wd_
  generic, public :: getipf_wd => getipf_wd_
  generic, public :: setipf_ht => setipf_ht_
  generic, public :: getipf_ht => getipf_ht_
  generic, public :: setpatx => setpatx_
  generic, public :: getpatx => getpatx_
  generic, public :: setpaty => setpaty_
  generic, public :: getpaty => getpaty_
  generic, public :: setN_ROI => setN_ROI_
  generic, public :: getN_ROI => getN_ROI_
  generic, public :: setsize_ROI => setsize_ROI_
  generic, public :: getsize_ROI => getsize_ROI_
  generic, public :: setroi_distance => setroi_distance_
  generic, public :: getroi_distance => getroi_distance_
  generic, public :: setnthreads => setnthreads_
  generic, public :: getnthreads => getnthreads_
  generic, public :: setstep_size => setstep_size_
  generic, public :: getstep_size => getstep_size_
  generic, public :: setdelta => setdelta_
  generic, public :: getdelta => getdelta_
  generic, public :: setthetac => setthetac_
  generic, public :: getthetac => getthetac_
  generic, public :: sethighpass => sethighpass_
  generic, public :: gethighpass => gethighpass_
  generic, public :: setlowpass => setlowpass_
  generic, public :: getlowpass => getlowpass_
  generic, public :: setC11 => setC11_
  generic, public :: getC11 => getC11_
  generic, public :: setC12 => setC12_
  generic, public :: getC12 => getC12_
  generic, public :: setC44 => setC44_
  generic, public :: getC44 => getC44_
  generic, public :: setC13 => setC13_
  generic, public :: getC13 => getC13_
  generic, public :: setC33 => setC33_
  generic, public :: getC33 => getC33_
  generic, public :: setdatafile => setdatafile_
  generic, public :: getdatafile => getdatafile_
  generic, public :: setmasterfile => setmasterfile_
  generic, public :: getmasterfile => getmasterfile_
  generic, public :: setexptfile => setexptfile_
  generic, public :: getexptfile => getexptfile_
  generic, public :: setinputtype => setinputtype_
  generic, public :: getinputtype => getinputtype_
  generic, public :: setHDFstrings => setHDFstrings_
  generic, public :: getHDFstrings => getHDFstrings_
  generic, public :: setRemap => setRemap_
  generic, public :: getRemap => getRemap_
  generic, public :: setPCrefine => setPCrefine_
  generic, public :: getPCrefine => getPCrefine_
  generic, public :: setcrystal => setcrystal_
  generic, public :: getcrystal => getcrystal_
end type HREBSD_T

! the constructor routine for this class 
interface HREBSD_T
  module procedure HREBSD_constructor
end interface HREBSD_T

contains

!--------------------------------------------------------------------------
type(HREBSD_T) function HREBSD_constructor( nmlfile ) result(HREBSD)
!! author: MDG 
!! version: 1.0 
!! date: 06/28/23
!!
!! constructor for the HREBSD_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call HREBSD%readNameList(nmlfile)

end function HREBSD_constructor

!--------------------------------------------------------------------------
subroutine HREBSD_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 06/28/23
!!
!! destructor for the HREBSD_T Class
 
IMPLICIT NONE

type(HREBSD_T), INTENT(INOUT)  :: self 

call reportDestructor('HREBSD_T')

end subroutine HREBSD_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 06/28/23
!!
!! read the namelist from an nml file for the HREBSD_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(HREBSD_T), INTENT(INOUT)       :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg)                    :: numsx
integer(kind=irg)                    :: numsy
integer(kind=irg)                    :: ipf_wd
integer(kind=irg)                    :: ipf_ht
integer(kind=irg)                    :: patx
integer(kind=irg)                    :: paty
integer(kind=irg)                    :: N_ROI
integer(kind=irg)                    :: size_ROI
integer(kind=irg)                    :: roi_distance
integer(kind=irg)                    :: nthreads
real(kind=sgl)                       :: step_size
real(kind=sgl)                       :: delta
real(kind=sgl)                       :: thetac
real(kind=sgl)                       :: highpass
real(kind=sgl)                       :: lowpass
real(kind=sgl)                       :: C11
real(kind=sgl)                       :: C12
real(kind=sgl)                       :: C44
real(kind=sgl)                       :: C13
real(kind=sgl)                       :: C33
character(fnlen)                     :: datafile
character(fnlen)                     :: masterfile
character(fnlen)                     :: exptfile
character(fnlen)                     :: inputtype
character(fnlen)                     :: HDFstrings(10)
character(1)                         :: Remap
character(1)                         :: PCrefine
character(3)                         :: crystal

namelist / HREBSDdata / numsx, numsy, ipf_wd, ipf_ht, patx, paty, N_ROI, size_ROI, roi_distance, &
                        nthreads, step_size, delta, thetac, highpass, lowpass, C11, C12, C44, C13, &
                        C33, datafile, masterfile, exptfile, inputtype, HDFstrings, Remap, &
                        PCrefine, crystal

numsx = 1244
numsy = 1024
ipf_wd = 1
ipf_ht = 1
patx = 0
paty = 0
N_ROI = 21
size_ROI = 8 
roi_distance = 250
nthreads = 1
step_size = 1.0
delta = 50.0
thetac = 10.0
highpass = 0.05
lowpass = 0.3
C11 = 276.0
C12 = 159.0
C44 = 132.0
C13 = 0.0
C33 = 0.0
datafile = 'undefined'
masterfile = 'undefined'
exptfile = 'undefined'
inputtype = 'Binary'
HDFstrings = (/ '', '', '', '', '', '', '', '', '', '' /)
Remap = 'n'
PCrefine = 'y'
crystal = 'cub'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=HREBSDdata)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries

! we no longer require the energyfile parameter, but for backwards compatibility
! we still allow the user to include it (it doesn't do anything though)
! if (trim(energyfile).eq.'undefined') then
!  call FatalError('GetEBSDNameList:',' energy file name is undefined in '//nmlfile)
! end if

 if (trim(exptfile).eq.'undefined') then
  call Message%printError('readNameList:',' exptfile file name is undefined in '//nmlfile)
 end if

 if (trim(masterfile).eq.'undefined') then
  call Message%printError('readNameList:',' master pattern file name is undefined in '//nmlfile)
 end if

 if (trim(datafile).eq.'undefined') then
  call Message%printError('readNameList:',' output file name is undefined in '//nmlfile)
 end if
end if

self%nml%numsx = numsx
self%nml%numsy = numsy
self%nml%ipf_wd = ipf_wd
self%nml%ipf_ht = ipf_ht
self%nml%patx = patx
self%nml%paty = paty
self%nml%N_ROI = N_ROI
self%nml%size_ROI = size_ROI
self%nml%roi_distance = roi_distance
self%nml%nthreads = nthreads
self%nml%step_size = step_size
self%nml%delta = delta
self%nml%thetac = thetac
self%nml%highpass = highpass
self%nml%lowpass = lowpass
self%nml%C11 = C11
self%nml%C12 = C12
self%nml%C44 = C44
self%nml%C13 = C13
self%nml%C33 = C33
self%nml%datafile = datafile
self%nml%masterfile = masterfile
self%nml%exptfile = exptfile
self%nml%inputtype = inputtype
self%nml%HDFstrings = HDFstrings
self%nml%Remap = Remap
self%nml%PCrefine = PCrefine
self%nml%crystal = crystal

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 06/28/23
!!
!! pass the namelist for the HREBSD_T Class to the calling program

IMPLICIT NONE 

class(HREBSD_T), INTENT(INOUT)          :: self
type(HREBSDNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 06/28/23
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)          :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 10, n_real = 10
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1),line10(10)

associate( enl => self%nml )

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! write all the single integers
io_int = (/ enl%numsx, enl%numsy, enl%ipf_wd, enl%ipf_ht, enl%nthreads, enl%patx, enl%paty, &
            enl%N_ROI, enl%size_ROI, enl%roi_distance /)
intlist(1) = 'numsx'
intlist(2) = 'numsy'
intlist(3) = 'ipf_wd'
intlist(4) = 'ipf_ht'
intlist(5) = 'nthreads'
intlist(6) = 'patx'
intlist(7) = 'paty'
intlist(8) = 'N_ROI'
intlist(9) = 'size_ROI'
intlist(10) = 'roi_distance'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ enl%step_size, enl%delta, enl%thetac, enl%highpass, enl%lowpass, &
             enl%C11, enl%C12, enl%C44, enl%C13, enl%C33 /)
reallist(1) = 'step_size'
reallist(2) = 'delta'
reallist(3) = 'thetac'
reallist(4) = 'highpass'
reallist(5) = 'lowpass'
reallist(6) = 'C11'
reallist(7) = 'C12'
reallist(8) = 'C44'
reallist(9) = 'C13'
reallist(10)= 'C33'
call HDF%writeNMLreals(io_real, reallist, n_real)

! write all the strings
dataset = SC_datafile
line2(1) = trim(enl%datafile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create datafile dataset', hdferr)

dataset = SC_masterfile
line2(1) = trim(enl%masterfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create masterfile dataset', hdferr)

dataset = SC_exptfile
line2(1) = trim(enl%exptfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create exptfile dataset', hdferr)

dataset = SC_inputtype
line2(1) = trim(enl%inputtype)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create inputtype dataset', hdferr)

dataset = 'Remap'
line2(1) = trim(enl%Remap)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create Remap dataset', hdferr)

dataset = 'PCrefine'
line2(1) = trim(enl%PCrefine)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create PCrefine dataset', hdferr)

dataset = 'crystal'
line2(1) = trim(enl%crystal)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create PCrefine crystal', hdferr)

dataset = SC_HDFstrings
line10 = enl%HDFstrings
hdferr = HDF%writeDatasetStringArray(dataset, line10, 10)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create HDFstrings dataset', hdferr)

! pop this group off the stack
call HDF%pop()
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine setnumsx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnumsx_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set numsx in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%numsx = inp

end subroutine setnumsx_

!--------------------------------------------------------------------------
function getnumsx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnumsx_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get numsx from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%numsx

end function getnumsx_

!--------------------------------------------------------------------------
subroutine setnumsy_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnumsy_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set numsy in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%numsy = inp

end subroutine setnumsy_

!--------------------------------------------------------------------------
function getnumsy_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnumsy_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get numsy from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%numsy

end function getnumsy_

!--------------------------------------------------------------------------
subroutine setipf_wd_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setipf_wd_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set ipf_wd in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%ipf_wd = inp

end subroutine setipf_wd_

!--------------------------------------------------------------------------
function getipf_wd_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getipf_wd_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get ipf_wd from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%ipf_wd

end function getipf_wd_

!--------------------------------------------------------------------------
subroutine setipf_ht_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setipf_ht_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set ipf_ht in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%ipf_ht = inp

end subroutine setipf_ht_

!--------------------------------------------------------------------------
function getipf_ht_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getipf_ht_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get ipf_ht from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%ipf_ht

end function getipf_ht_

!--------------------------------------------------------------------------
subroutine setpatx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setpatx_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set patx in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%patx = inp

end subroutine setpatx_

!--------------------------------------------------------------------------
function getpatx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getpatx_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get patx from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%patx

end function getpatx_

!--------------------------------------------------------------------------
subroutine setpaty_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setpaty_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set paty in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%paty = inp

end subroutine setpaty_

!--------------------------------------------------------------------------
function getpaty_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getpaty_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get paty from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%paty

end function getpaty_

!--------------------------------------------------------------------------
subroutine setN_ROI_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setN_ROI_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set N_ROI in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%N_ROI = inp

end subroutine setN_ROI_

!--------------------------------------------------------------------------
function getN_ROI_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getN_ROI_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get N_ROI from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%N_ROI

end function getN_ROI_

!--------------------------------------------------------------------------
subroutine setsize_ROI_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setsize_ROI_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set size_ROI in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%size_ROI = inp

end subroutine setsize_ROI_

!--------------------------------------------------------------------------
function getsize_ROI_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getsize_ROI_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get size_ROI from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%size_ROI

end function getsize_ROI_

!--------------------------------------------------------------------------
subroutine setroi_distance_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setroi_distance_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set roi_distance in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%roi_distance = inp

end subroutine setroi_distance_

!--------------------------------------------------------------------------
function getroi_distance_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getroi_distance_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get roi_distance from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%roi_distance

end function getroi_distance_

!--------------------------------------------------------------------------
subroutine setnthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnthreads_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set nthreads in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%nthreads = inp

end subroutine setnthreads_

!--------------------------------------------------------------------------
function getnthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnthreads_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get nthreads from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%nthreads

end function getnthreads_

!--------------------------------------------------------------------------
subroutine setstep_size_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setstep_size_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set step_size in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%step_size = inp

end subroutine setstep_size_

!--------------------------------------------------------------------------
function getstep_size_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getstep_size_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get step_size from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%step_size

end function getstep_size_

!--------------------------------------------------------------------------
subroutine setdelta_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdelta_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set delta in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%delta = inp

end subroutine setdelta_

!--------------------------------------------------------------------------
function getdelta_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdelta_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get delta from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%delta

end function getdelta_

!--------------------------------------------------------------------------
subroutine setthetac_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setthetac_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set thetac in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%thetac = inp

end subroutine setthetac_

!--------------------------------------------------------------------------
function getthetac_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getthetac_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get thetac from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%thetac

end function getthetac_

!--------------------------------------------------------------------------
subroutine sethighpass_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: sethighpass_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set highpass in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%highpass = inp

end subroutine sethighpass_

!--------------------------------------------------------------------------
function gethighpass_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: gethighpass_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get highpass from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%highpass

end function gethighpass_

!--------------------------------------------------------------------------
subroutine setlowpass_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setlowpass_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set lowpass in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%lowpass = inp

end subroutine setlowpass_

!--------------------------------------------------------------------------
function getlowpass_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getlowpass_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get lowpass from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%lowpass

end function getlowpass_

!--------------------------------------------------------------------------
subroutine setC11_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setC11_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set C11 in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%C11 = inp

end subroutine setC11_

!--------------------------------------------------------------------------
function getC11_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getC11_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get C11 from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%C11

end function getC11_

!--------------------------------------------------------------------------
subroutine setC12_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setC12_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set C12 in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%C12 = inp

end subroutine setC12_

!--------------------------------------------------------------------------
function getC12_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getC12_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get C12 from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%C12

end function getC12_

!--------------------------------------------------------------------------
subroutine setC44_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setC44_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set C44 in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%C44 = inp

end subroutine setC44_

!--------------------------------------------------------------------------
function getC44_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getC44_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get C44 from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%C44

end function getC44_

!--------------------------------------------------------------------------
subroutine setC13_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setC13_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set C13 in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%C13 = inp

end subroutine setC13_

!--------------------------------------------------------------------------
function getC13_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getC13_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get C13 from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%C13

end function getC13_

!--------------------------------------------------------------------------
subroutine setC33_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setC33_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set C33 in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%C33 = inp

end subroutine setC33_

!--------------------------------------------------------------------------
function getC33_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getC33_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get C33 from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%C33

end function getC33_

!--------------------------------------------------------------------------
subroutine setdatafile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdatafile_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set datafile in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%datafile = trim(inp)

end subroutine setdatafile_

!--------------------------------------------------------------------------
function getdatafile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdatafile_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get datafile from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%datafile)

end function getdatafile_

!--------------------------------------------------------------------------
subroutine setmasterfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setmasterfile_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set masterfile in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%masterfile = trim(inp)

end subroutine setmasterfile_

!--------------------------------------------------------------------------
function getmasterfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getmasterfile_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get masterfile from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%masterfile)

end function getmasterfile_

!--------------------------------------------------------------------------
subroutine setexptfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setexptfile_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set exptfile in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%exptfile = trim(inp)

end subroutine setexptfile_

!--------------------------------------------------------------------------
function getexptfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getexptfile_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get exptfile from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%exptfile)

end function getexptfile_

!--------------------------------------------------------------------------
subroutine setinputtype_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setinputtype_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set inputtype in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%inputtype = trim(inp)

end subroutine setinputtype_

!--------------------------------------------------------------------------
function getinputtype_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getinputtype_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get inputtype from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%inputtype)

end function getinputtype_

!--------------------------------------------------------------------------
subroutine setHDFstrings_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setHDFstrings_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set HDFstrings in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp(10)

integer(kind=irg)                  :: i

do i=1,10
  self%nml%HDFstrings(i) = trim(inp(i))
end do 

end subroutine setHDFstrings_

!--------------------------------------------------------------------------
function getHDFstrings_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getHDFstrings_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get HDFstrings from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out(10)

integer(kind=irg)                  :: i

do i=1,10
  out(i) = trim(self%nml%HDFstrings(i))
end do 

end function getHDFstrings_

!--------------------------------------------------------------------------
subroutine setRemap_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setRemap_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set Remap in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
character(1), INTENT(IN)       :: inp

self%nml%Remap = trim(inp)

end subroutine setRemap_

!--------------------------------------------------------------------------
function getRemap_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getRemap_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get Remap from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
character(1)                   :: out

out = trim(self%nml%Remap)

end function getRemap_

!--------------------------------------------------------------------------
subroutine setPCrefine_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setPCrefine_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set PCrefine in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
character(1), INTENT(IN)       :: inp

self%nml%PCrefine = trim(inp)

end subroutine setPCrefine_

!--------------------------------------------------------------------------
function getPCrefine_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getPCrefine_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get PCrefine from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
character(1)                   :: out

out = trim(self%nml%PCrefine)

end function getPCrefine_

!--------------------------------------------------------------------------
subroutine setcrystal_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setcrystal_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! set crystal in the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
character(3), INTENT(IN)       :: inp

self%nml%crystal = trim(inp)

end subroutine setcrystal_

!--------------------------------------------------------------------------
function getcrystal_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getcrystal_
!! author: MDG
!! version: 1.0
!! date: 06/28/23
!!
!! get crystal from the HREBSD_T class

IMPLICIT NONE

class(HREBSD_T), INTENT(INOUT)     :: self
character(3)                   :: out

out = trim(self%nml%crystal)

end function getcrystal_

!--------------------------------------------------------------------------
subroutine HREBSD_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: HREBSD_
!! author: MDG (OO version) / Chaoyi Zhu (original version) 
!! version: 1.0 
!! date: 06/28/23
!!
!! perform the computations

use mod_EMsoft
use mod_HDFnames
use mod_quaternions
use mod_rotations 
use mod_io
use mod_vendors
use mod_HDFsupport
use HDF5
use ISO_C_BINDING

IMPLICIT NONE 

class(HREBSD_T), INTENT(INOUT)          :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

type(Vendor_T)                          :: VT
type(QuaternionArray_T)                 :: qAR
type(IO_T)                              :: Message

character(fnlen)                        :: inpfile, HDFstring
real(kind=sgl)                          :: stepsizes(3), fpar(3)
real(kind=sgl),allocatable              :: PC(:,:)
integer(kind=irg)                       :: hdferr


associate( enl => self%nml )

! this program reads a pattern file, including all the experimental parameters' so 
! this is essentially limited to TSLHDF, EDAXH5, OxfordHDF, and BrukerHDF at the moment.  We can also 
! read EMEBSD output, which is useful for debugging and testing purposes. 

call setRotationPrecision('d')

inpfile = EMsoft%generateFilePath('EMdatapathname',trim(enl%exptfile))
fpar = (/ real(enl%numsx), real(enl%numsy), real(enl%delta) /)
HDFstring = trim(enl%HDFstrings(1))

call openFortranHDFInterface()

VT = Vendor_T( enl%inputtype ) 

select case (trim(enl%inputtype)) 
  case ('TSLHDF')
    call VT%getTSLmetadata(inpfile, HDFstring, stepsizes, qAR, PC, fpar)

  case ('EDAXH5')
    call VT%getEDAXH5metadata(inpfile, enl%HDFstrings, stepsizes, qAR, PC, fpar)

  case ('OxfordHDF')
    call VT%getOxfordmetadata(inpfile, HDFstring, stepsizes, qAR, PC, fpar)

  case ('BrukerHDF')
    call VT%getBrukermetadata(inpfile, HDFstring, stepsizes, qAR, PC, fpar, self%SEM%SEM_X, self%SEM%SEM_Y)

  case ('EMEBSD')
    call VT%getEMsoftmetadata(inpfile, HDFstring, stepsizes, qAR, PC, fpar)

  case default 
    call Message%printError('HREBSD_','input format currently not supported in EMHREBSD')
end select

call Message%printMessage(' Completed reading orientation data')




call closeFortranHDFInterface()

end associate 

end subroutine HREBSD_



end module mod_HREBSD