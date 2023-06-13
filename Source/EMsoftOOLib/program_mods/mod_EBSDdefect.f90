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

module mod_EBSDdefect
  !! author: MDG 
  !! version: 1.0 
  !! date: 06/05/23
  !!
  !! class definition for the EMEBSDdefect program

use mod_kinds
use mod_global

IMPLICIT NONE 

type, public :: EBSDdefectNameListType
   real(kind=sgl)         :: L
   real(kind=sgl)         :: thetac
   real(kind=sgl)         :: delta
   real(kind=sgl)         :: xpc
   real(kind=sgl)         :: ypc
   real(kind=sgl)         :: gammavalue
   real(kind=sgl)         :: DF_L
   real(kind=dbl)         :: k(3)
   real(kind=dbl)         :: q(3)
   integer(kind=irg)      :: numsx
   integer(kind=irg)      :: numsy
   integer(kind=irg)      :: DF_npix
   integer(kind=irg)      :: DF_npiy
   integer(kind=irg)      :: nthreads
   character(3)           :: scalingmode
   character(fnlen)       :: masterfile
   character(fnlen)       :: datafile
   character(fnlen)       :: defectfilename
end type EBSDdefectNameListType

! class definition
type, public :: EBSDdefect_T
private 
  character(fnlen)       :: nmldeffile = 'EMEBSDdefect.nml'
  type(EBSDdefectNameListType)  :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: EBSDdefect_
  procedure, pass(self) :: setL_
  procedure, pass(self) :: getL_
  procedure, pass(self) :: setthetac_
  procedure, pass(self) :: getthetac_
  procedure, pass(self) :: setdelta_
  procedure, pass(self) :: getdelta_
  procedure, pass(self) :: setxpc_
  procedure, pass(self) :: getxpc_
  procedure, pass(self) :: setypc_
  procedure, pass(self) :: getypc_
  procedure, pass(self) :: setgammavalue_
  procedure, pass(self) :: getgammavalue_
  procedure, pass(self) :: setDF_L_
  procedure, pass(self) :: getDF_L_
  procedure, pass(self) :: setnumsx_
  procedure, pass(self) :: getnumsx_
  procedure, pass(self) :: setnumsy_
  procedure, pass(self) :: getnumsy_
  procedure, pass(self) :: setk_
  procedure, pass(self) :: getk_
  procedure, pass(self) :: setq_
  procedure, pass(self) :: getq_
  procedure, pass(self) :: setDF_npix_
  procedure, pass(self) :: getDF_npix_
  procedure, pass(self) :: setDF_npiy_
  procedure, pass(self) :: getDF_npiy_
  procedure, pass(self) :: setnthreads_
  procedure, pass(self) :: getnthreads_
  procedure, pass(self) :: setscalingmode_
  procedure, pass(self) :: getscalingmode_
  procedure, pass(self) :: setmasterfile_
  procedure, pass(self) :: getmasterfile_
  procedure, pass(self) :: setdatafile_
  procedure, pass(self) :: getdatafile_
  procedure, pass(self) :: setdefectfilename_
  procedure, pass(self) :: getdefectfilename_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: EBSDdefect => EBSDdefect_
  generic, public :: setL => setL_
  generic, public :: getL => getL_
  generic, public :: setthetac => setthetac_
  generic, public :: getthetac => getthetac_
  generic, public :: setdelta => setdelta_
  generic, public :: getdelta => getdelta_
  generic, public :: setxpc => setxpc_
  generic, public :: getxpc => getxpc_
  generic, public :: setypc => setypc_
  generic, public :: getypc => getypc_
  generic, public :: setgammavalue => setgammavalue_
  generic, public :: getgammavalue => getgammavalue_
  generic, public :: setDF_L => setDF_L_
  generic, public :: getDF_L => getDF_L_
  generic, public :: setnumsx => setnumsx_
  generic, public :: getnumsx => getnumsx_
  generic, public :: setnumsy => setnumsy_
  generic, public :: getnumsy => getnumsy_
  generic, public :: setk => setk_
  generic, public :: getk => getk_
  generic, public :: setq => setq_
  generic, public :: getq => getq_
  generic, public :: setDF_npix => setDF_npix_
  generic, public :: getDF_npix => getDF_npix_
  generic, public :: setDF_npiy => setDF_npiy_
  generic, public :: getDF_npiy => getDF_npiy_
  generic, public :: setnthreads => setnthreads_
  generic, public :: getnthreads => getnthreads_
  generic, public :: setscalingmode => setscalingmode_
  generic, public :: getscalingmode => getscalingmode_
  generic, public :: setmasterfile => setmasterfile_
  generic, public :: getmasterfile => getmasterfile_
  generic, public :: setdatafile => setdatafile_
  generic, public :: getdatafile => getdatafile_
  generic, public :: setdefectfilename => setdefectfilename_
  generic, public :: getdefectfilename => getdefectfilename_

end type EBSDdefect_T

! the constructor routine for this class 
interface EBSDdefect_T
  module procedure EBSDdefect_constructor
end interface EBSDdefect_T

contains

!--------------------------------------------------------------------------
type(EBSDdefect_T) function EBSDdefect_constructor( nmlfile ) result(EBSDdefect)
!! author: MDG 
!! version: 1.0 
!! date: 06/05/23
!!
!! constructor for the EBSDdefect_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call EBSDdefect%readNameList(nmlfile)

end function EBSDdefect_constructor

!--------------------------------------------------------------------------
subroutine EBSDdefect_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 06/05/23
!!
!! destructor for the EBSDdefect_T Class
 
IMPLICIT NONE

type(EBSDdefect_T), INTENT(INOUT)  :: self 

call reportDestructor('EBSDdefect_T')

end subroutine EBSDdefect_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 06/05/23
!!
!! read the namelist from an nml file for the EBSDdefect_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(EBSDdefect_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

real(kind=sgl)         :: L
real(kind=sgl)         :: thetac
real(kind=sgl)         :: delta
real(kind=sgl)         :: xpc
real(kind=sgl)         :: ypc
real(kind=sgl)         :: gammavalue
real(kind=sgl)         :: rotang
real(kind=sgl)         :: DF_L
real(kind=dbl)         :: k(3)
real(kind=dbl)         :: q(3)
integer(kind=irg)      :: numsx
integer(kind=irg)      :: numsy
integer(kind=irg)      :: DF_npix
integer(kind=irg)      :: DF_npiy
integer(kind=irg)      :: nthreads
character(3)           :: scalingmode
character(fnlen)       :: masterfile
character(fnlen)       :: datafile
character(fnlen)       :: defectfilename

namelist / EBSDdefectdata / L, thetac, delta, xpc, ypc, gammavalue, rotang, DF_L, &
                            numsx, numsy, k, q, DF_npix, DF_npiy, nthreads, scalingmode, &
                            masterfile, datafile, defectfilename

! set the input parameters to default values (except for xtalname, which must be present)
L = 15000.0
thetac = 10.0
delta = 50.0
numsx = 0
numsy = 0
xpc = 0.0
ypc = 0.0
scalingmode = 'not'
gammavalue = 1.0
k = (/ 0.D0,0.D0,1.D0 /)
q = (/ 1.D0,0.D0,0.D0 /)
rotang = 0.0
DF_L = 1.0
DF_npix = 256
DF_npiy = 256
masterfile = 'master.h5'
datafile = 'EBSDdefectout.h5'
defectfilename = 'EMdefect.json'
nthreads = 1

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=EBSDdefectdata)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries

 if (trim(masterfile).eq.'undefined') then
  call Message%printError('readNameList:',' masterfile file name is undefined in '//nmlfile)
 end if

 if (trim(datafile).eq.'undefined') then
  call Message%printError('readNameList:',' datafile file name is undefined in '//nmlfile)
 end if

 if (trim(defectfilename).eq.'undefined') then
  call Message%printError('readNameList:',' defectfilename file name is undefined in '//nmlfile)
 end if
end if 

self%nml%L = L
self%nml%thetac = thetac 
self%nml%delta = delta
self%nml%xpc = xpc
self%nml%ypc = ypc
self%nml%gammavalue = gammavalue
self%nml%DF_L = DF_L
self%nml%numsx = numsx
self%nml%numsy = numsy
self%nml%k = k
self%nml%q = q
self%nml%DF_npix = DF_npix
self%nml%DF_npiy = DF_npiy
self%nml%nthreads = nthreads
self%nml%scalingmode = scalingmode
self%nml%masterfile = masterfile
self%nml%datafile = datafile
self%nml%defectfilename = defectfilename

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 06/05/23
!!
!! pass the namelist for the EBSDdefect_T Class to the calling program

IMPLICIT NONE 

class(EBSDdefect_T), INTENT(INOUT)          :: self
type(EBSDdefectNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 06/05/23
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)        :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 5, n_real = 7
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( enl => self%nml )

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! write all the single integers
io_int = (/ enl%numsx, enl%numsy, enl%nthreads, enl%DF_npix, enl%DF_npix /)
intlist(1) = 'numsx'
intlist(2) = 'numsy'
intlist(3) = 'nthreads'
intlist(4) = 'DF_npix'
intlist(5) = 'DF_npiy'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ enl%L, enl%thetac, enl%delta, enl%xpc, enl%ypc, enl%gammavalue, &
             enl%DF_L /)
reallist(1) = 'L'
reallist(2) = 'thetac'
reallist(3) = 'delta'
reallist(4) = 'xpc'
reallist(5) = 'ypc'
reallist(6) = 'gammavalue'
reallist(7) = 'DF_L'
call HDF%writeNMLreals(io_real, reallist, n_real)

! a 3-vector
dataset = SC_k
hdferr = HDF%writeDatasetDoubleArray(dataset, enl%k, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create k dataset', hdferr)

dataset = 'q'
hdferr = HDF%writeDatasetDoubleArray(dataset, enl%q, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create q dataset', hdferr)

! write all the strings
dataset = SC_scalingmode
line2(1) = trim(enl%scalingmode)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create scalingmode dataset', hdferr)

dataset = SC_masterfile
line2(1) = trim(enl%masterfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create masterfile dataset', hdferr)

dataset = SC_datafile
line2(1) = trim(enl%datafile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create datafile dataset', hdferr)

dataset = SC_defectfilename
line2(1) = trim(enl%defectfilename)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create defectfilename dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine setL_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setL_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! set L in the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%L = inp

end subroutine setL_

!--------------------------------------------------------------------------
function getL_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getL_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! get L from the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%L

end function getL_

!--------------------------------------------------------------------------
subroutine setthetac_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setthetac_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! set thetac in the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%thetac = inp

end subroutine setthetac_

!--------------------------------------------------------------------------
function getthetac_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getthetac_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! get thetac from the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%thetac

end function getthetac_

!--------------------------------------------------------------------------
subroutine setdelta_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdelta_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! set delta in the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%delta = inp

end subroutine setdelta_

!--------------------------------------------------------------------------
function getdelta_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdelta_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! get delta from the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%delta

end function getdelta_

!--------------------------------------------------------------------------
subroutine setxpc_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setxpc_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! set xpc in the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%xpc = inp

end subroutine setxpc_

!--------------------------------------------------------------------------
function getxpc_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getxpc_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! get xpc from the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%xpc

end function getxpc_

!--------------------------------------------------------------------------
subroutine setypc_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setypc_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! set ypc in the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%ypc = inp

end subroutine setypc_

!--------------------------------------------------------------------------
function getypc_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getypc_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! get ypc from the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%ypc

end function getypc_

!--------------------------------------------------------------------------
subroutine setgammavalue_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setgammavalue_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! set gammavalue in the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%gammavalue = inp

end subroutine setgammavalue_

!--------------------------------------------------------------------------
function getgammavalue_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getgammavalue_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! get gammavalue from the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%gammavalue

end function getgammavalue_

!--------------------------------------------------------------------------
subroutine setDF_L_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setDF_L_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! set DF_L in the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%DF_L = inp

end subroutine setDF_L_

!--------------------------------------------------------------------------
function getDF_L_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getDF_L_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! get DF_L from the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%DF_L

end function getDF_L_

!--------------------------------------------------------------------------
subroutine setnumsx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnumsx_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! set numsx in the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%numsx = inp

end subroutine setnumsx_

!--------------------------------------------------------------------------
function getnumsx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnumsx_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! get numsx from the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%numsx

end function getnumsx_

!--------------------------------------------------------------------------
subroutine setnumsy_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnumsy_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! set numsy in the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%numsy = inp

end subroutine setnumsy_

!--------------------------------------------------------------------------
function getnumsy_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnumsy_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! get numsy from the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%numsy

end function getnumsy_

!--------------------------------------------------------------------------
subroutine setk_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setk_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! set k in the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)    :: self
integer(kind=irg), INTENT(IN)       :: inp(3)

self%nml%k = inp

end subroutine setk_

!--------------------------------------------------------------------------
function getk_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getk_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! get k from the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)    :: self
integer(kind=irg)                   :: out(3)

out = self%nml%k

end function getk_

!--------------------------------------------------------------------------
subroutine setq_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setq_
!! author: MDG
!! version: 1.0
!! date: 06/07/23
!!
!! set q in the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)    :: self
integer(kind=irg), INTENT(IN)       :: inp(3)

self%nml%q = inp

end subroutine setq_

!--------------------------------------------------------------------------
function getq_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getq_
!! author: MDG
!! version: 1.0
!! date: 06/07/23
!!
!! get q from the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)    :: self
integer(kind=irg)                   :: out(3)

out = self%nml%q

end function getq_

!--------------------------------------------------------------------------
subroutine setDF_npix_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setDF_npix_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! set DF_npix in the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%DF_npix = inp

end subroutine setDF_npix_

!--------------------------------------------------------------------------
function getDF_npix_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getDF_npix_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! get DF_npix from the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%DF_npix

end function getDF_npix_

!--------------------------------------------------------------------------
subroutine setDF_npiy_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setDF_npiy_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! set DF_npiy in the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%DF_npiy = inp

end subroutine setDF_npiy_

!--------------------------------------------------------------------------
function getDF_npiy_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getDF_npiy_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! get DF_npiy from the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%DF_npiy

end function getDF_npiy_

!--------------------------------------------------------------------------
subroutine setnthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnthreads_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! set nthreads in the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%nthreads = inp

end subroutine setnthreads_

!--------------------------------------------------------------------------
function getnthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnthreads_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! get nthreads from the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%nthreads

end function getnthreads_

!--------------------------------------------------------------------------
subroutine setscalingmode_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setscalingmode_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! set scalingmode in the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
character(3), INTENT(IN)       :: inp

self%nml%scalingmode = trim(inp)

end subroutine setscalingmode_

!--------------------------------------------------------------------------
function getscalingmode_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getscalingmode_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! get scalingmode from the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
character(3)                   :: out

out = trim(self%nml%scalingmode)

end function getscalingmode_

!--------------------------------------------------------------------------
subroutine setmasterfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setmasterfile_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! set masterfile in the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%masterfile = trim(inp)

end subroutine setmasterfile_

!--------------------------------------------------------------------------
function getmasterfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getmasterfile_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! get masterfile from the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%masterfile)

end function getmasterfile_

!--------------------------------------------------------------------------
subroutine setdatafile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdatafile_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! set datafile in the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%datafile = trim(inp)

end subroutine setdatafile_

!--------------------------------------------------------------------------
function getdatafile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdatafile_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! get datafile from the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%datafile)

end function getdatafile_

!--------------------------------------------------------------------------
subroutine setdefectfilename_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdefectfilename_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! set defectfilename in the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%defectfilename = trim(inp)

end subroutine setdefectfilename_

!--------------------------------------------------------------------------
function getdefectfilename_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdefectfilename_
!! author: MDG
!! version: 1.0
!! date: 06/05/23
!!
!! get defectfilename from the EBSDdefect class

IMPLICIT NONE

class(EBSDdefect_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%defectfilename)

end function getdefectfilename_

!--------------------------------------------------------------------------
subroutine EBSDdefect_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDdefect_
!! author: MDG 
!! version: 1.0 
!! date: 06/05/23
!!
!! perform the computations

use mod_EMsoft
use mod_so3
use mod_quaternions
use mod_MCfiles
use mod_MPfiles
use HDF5
use mod_HDFsupport
use mod_HDFnames
use omp_lib
use mod_io
use mod_rotations
use stringconstants
use mod_memory
use mod_symmetry
use mod_math
use mod_crystallography
use mod_timing
use mod_defect
use ISO_C_BINDING
use h5im
use h5lt
use mod_image

use, intrinsic :: iso_fortran_env

IMPLICIT NONE 

class(EBSDdefect_T), INTENT(INOUT)      :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

type(MCfile_T)                          :: MCFT
type(MPfile_T)                          :: MPFT
type(HDF_T)                             :: HDF
type(HDFnames_T)                        :: saveHDFnames
type(so3_T)                             :: SO
type(IO_T)                              :: Message
type(Quaternion_T)                      :: quat
type(QuaternionArray_T)                 :: qAR
type(memory_T)                          :: mem
type(SpaceGroup_T)                      :: SG
type(Cell_T)                            :: cell
type(q_T)                               :: q
type(o_T)                               :: o
type(e_T)                               :: eu
type(Timing_T)                          :: timer
type(Defect_T)                          :: Defects

type(SEMmasterNameListType)             :: mpnl
type(MCOpenCLNameListType)              :: mcnl

logical                                 :: verbose, insert = .TRUE., overwrite = .TRUE., g_exists
character(fnlen)                        :: fname, nmldeffile, datafile
integer(kind=irg)                       :: numangles, istat, i, j, k, error_cnt, npix, npiy, hdferr, ix, iy
type(FZpointd),pointer                  :: FZtmp
type(r_T)                               :: rr
integer(kind=irg)                       :: ga(3), gb(3), io_int(6), sh(3), ipar(8), iipar(8)
real(kind=dbl)                          :: kc(3), gac(3), gbc(3), FF(3,3), FF_inv(3,3), prefactor
real(kind=dbl)                          :: om(3,3), pctr(3)
real(kind=sgl)                          :: io_real(1), xpos, ypos, mi, ma
real(kind=sgl),allocatable              :: patarray(:,:,:,:), trial(:,:,:,:), quarray(:,:), binned(:,:)
real(kind=dbl),allocatable              :: tFij(:,:,:), Fmatrix(:,:,:)
real(kind=sgl),allocatable              :: tmLPNH(:,:,:) , tmLPSH(:,:,:)
real(kind=sgl),allocatable              :: trgx(:,:), trgy(:,:), trgz(:,:) ! auxiliary detector arrays needed for interpolation
integer(kind=irg)                       :: NUMTHREADS, TID   ! number of allocated threads, thread ID
integer(kind=irg)                       :: nthreads

character(fnlen)                        :: groupname, dataset, datagroupname, attributename, HDF_FileVersion, TIFF_filename
character(11)                           :: dstr
character(15)                           :: tstrb
character(15)                           :: tstre
integer(HSIZE_T)                        :: dims4(4), cnt4(4), offset4(4)
character(fnlen,kind=c_char)            :: line2(1)
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)

! declare variables for use in object oriented image module
integer                                 :: iostat
character(len=128)                      :: iomsg
logical                                 :: isInteger
type(image_t)                           :: im
integer(int8)                           :: i8 (3,4)
integer(int8), allocatable              :: TIFF_image(:,:)

associate( enl => self%nml, EBSDMCdata => MCFT%MCDT )

! make sure that k and q are orthogonal
if (abs(DOT_PRODUCT(enl%k,enl%q)).gt.1.0e-8) then 
  call Message%printError('EMEBSDdefect','k and q vectors must be orthogonal')
end if 

call MPFT%setModality('EBSD')
nmldeffile = trim(EMsoft%nmldeffile)

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()
mem = memory_T()

! set the HDF group names for this program
HDFnames = HDFnames_T()

verbose = .TRUE.

call setRotationPrecision('d')

! 1. read EBSD defect master pattern file (HDF format)
! first the Monte Carlo data (mainly we need the sample tilt angle)
call HDFnames%set_ProgramData(SC_MCOpenCL)
call HDFnames%set_NMLlist(SC_MCCLNameList)
call HDFnames%set_NMLfilename(SC_MCOpenCLNML)
fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfile))
call MCFT%setFileName(fname)
call MCFT%readMCfile(HDF, HDFnames)
mcnl = MCFT%getnml()

! then the master pattern arrays
call HDFnames%set_ProgramData(SC_EBSDdepthmaster)
call HDFnames%set_NMLlist(SC_EBSDdepthmasterNameList)
call HDFnames%set_NMLfilename(SC_EBSDdepthmasterNML)
call HDFnames%set_Variable(SC_MCOpenCL)

fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfile))
call MPFT%setFileName(fname)
call MPFT%readMPfile(HDF, HDFnames, mpnl, defectMP=.TRUE., getmLPNH=.TRUE., getmLPSH=.TRUE.)

sh = shape(MPFT%MPDT%mLPNH)

! 2. set HDFnames for output
call HDFnames%set_ProgramData(SC_EBSDdefect)
call HDFnames%set_NMLlist(SC_EBSDdefectNameList)
call HDFnames%set_NMLfilename(SC_EBSDdefectNML)
call HDFnames%set_Variable(SC_MCOpenCL)


! 3. create the output HDF file 
! Create a new file using the default properties.
timer = Timing_T()
tstrb = timer%getTimeString()
dstr = timer%getDateString()

call Message%printMessage(' Creating HDF5 output file ')
io_real(1) = real(enl%numsx) * real(enl%numsy) * real(enl%DF_npix) * real(enl%DF_npiy) * 4.0
io_real(1) = io_real(1) / (1024.0)**3
call Message%WriteValue('   size of pattern output array (Gb): ', io_real, 1,"(f10.2)")

datafile = EMsoft%generateFilePath('EMdatapathname', enl%datafile)

hdferr =  HDF%createFile(datafile)
if (hdferr.ne.0) call HDF%error_check('HDF_createFile ', hdferr)

!====================================
! new in Release 4.3: add a Manufacturer string (null terminated)
dataset = SC_Manufacturer
line2(1) = 'EMsoftOO'
line2(1) = cstringify(line2(1))
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
!====================================

! write the EMheader to the file
datagroupname = trim(HDFnames%get_ProgramData()) ! 'EBSDdefect'
call HDF%writeEMheader(EMsoft, dstr, tstrb, tstre, progname, datagroupname)

! add the CrystalData group at the top level of the file
call cell%addXtalDataGroup(SG, EMsoft, HDF)

! create a namelist group to write all the namelist files into
hdferr = HDF%createGroup(HDFnames%get_NMLfiles())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup NMLfiles', hdferr)

! read the text file and write the array to the file
dataset = SC_EMEBSDdefectNML
hdferr = HDF%writeDatasetTextFile(dataset, nmldeffile)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetTextFile EMEBSDdefectNML', hdferr)

call HDF%pop()

! create a NMLparameters group to write all the namelist entries into
hdferr = HDF%createGroup(HDFnames%get_NMLparameters())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup NMLparameters', hdferr)

call self%writeHDFNameList_(HDF, HDFnames)

! and leave this group
call HDF%pop()

! then the remainder of the data in a EMData group
hdferr = HDF%createGroup(HDFnames%get_EMData())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup EMData', hdferr)

! create the EBSD group and add a HDF_FileVersion attribute to it
hdferr = HDF%createGroup(datagroupname)
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup EBSDdefect', hdferr)

HDF_FileVersion = '4.1'
attributename = SC_HDFFileVersion
hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)

! =====================================================
dataset = SC_xtalname
call mem%alloc(stringarray, (/ 1 /), 'stringarray')
stringarray(1)= trim(MPFT%MPDT%xtalname)
hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetStringArray xtalname', hdferr)

dataset = 'TransformationQuaternion'
hdferr = HDF%writeDatasetDoubleArray(dataset, quat%get_quatd(), 4 )
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetStringArray quat', hdferr)

! generate a 4D hyperslab array to store the individual EBSD patterns
dataset = SC_EBSDpatterns
  call mem%alloc(patarray, (/ enl%numsx, enl%numsy, enl%DF_npix, enl%DF_npiy /), 'patarray', initval = 0.0)
  dims4 = (/  enl%numsx, enl%numsy, enl%DF_npix, enl%DF_npiy /)
  cnt4 = (/ enl%numsx, enl%numsy, enl%DF_npix, 1 /)
  offset4 = (/ 0, 0, 0, 0 /)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeHyperslabFloatArray(dataset, patarray, dims4, offset4, cnt4, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, patarray, dims4, offset4, cnt4)
  end if
  call mem%dealloc(patarray, 'patarray')

!=============================================
!=============================================
! 4. Defects initialization section 
! get the crystal structure information to perform crystallographic computations
call cell%getCrystalData(MPFT%MPDT%xtalname, SG, EMsoft, useHDF=HDF)

! copy some of the namelist parameters into the defects structure
Defects%DF_npix = enl%DF_npix
Defects%DF_npiy = enl%DF_npiy
Defects%DF_L = enl%DF_L
Defects%DF_nums = sh(3)
Defects%DF_slice = mcnl%depthstep  ! integration slice thickness
Defects%APD%stepsize = 0.01D0 

! If there is no displacement field file we compute displacement field
call Message%printMessage(' --> initializing defects ')
  
! determine the cartesian components of ga (which should really be enl%q)
Defects%DF_gf = sngl(enl%q)  ! float(ga)
Defects%DF_gstar = Defects%DF_gf/cell%CalcLength(Defects%DF_gf,'r')**2    ! define G* such that G.G* = 1
call cell%TransSpace(Defects%DF_gf,Defects%DF_gc,'r','c')         ! convert to Cartesian reference frame

! next, we read all the foil and defect data using the new InitializeDefects routine in defectmodule.f90
verbose = .FALSE.
call Defects%InitializeDefects(EMsoft,cell,enl%defectfilename,enl%DF_npix,enl%DF_npiy,enl%DF_L,Defects%DF_gf, &
                               enl%k,enl%q,error_cnt,verbose)

! ok, all the set up is now complete;
npix = Defects%DF_npix
npiy = Defects%DF_npiy

ipar(1) = 1
ipar(2) = enl%numsx
ipar(3) = enl%numsy
ipar(4) = mpnl%npx
ipar(5) = mpnl%npx
ipar(6) = sh(3)
ipar(7) = 0
ipar(8) = 0
prefactor = 1.D0

! this array contains a row of EBSD patterns
call mem%alloc(patarray, (/ enl%numsx, enl%numsy, enl%DF_npix, 1 /), 'patarray', initval = 0.0)

call Message%printMessage(' --> starting parallel computation ...')
! use OpenMP to run on multiple cores ... 
!$OMP PARALLEL default(shared)  PRIVATE(TID, NUMTHREADS, i, j, binned, tmLPNH, tmLPSH, ix, pctr, iipar)

  NUMTHREADS = OMP_GET_NUM_THREADS()
  TID = OMP_GET_THREAD_NUM()

! each thread has its own copy of the ipar array
  iipar = ipar

! each thread needs a private copy of the master arrays; not having those can produce poor scaling...
  allocate(tmLPNH(-mpnl%npx:mpnl%npx,-mpnl%npx:mpnl%npx,ipar(7)), tmLPSH(-mpnl%npx:mpnl%npx,-mpnl%npx:mpnl%npx,ipar(7)))
  allocate(binned(enl%numsx,enl%numsy))

! and copy the data in
  tmLPNH = MPFT%MPDT%mLPNH
  tmLPSH = MPFT%MPDT%mLPSH

! we do this one row at a time
  do iy=1, enl%DF_npiy
!$OMP DO SCHEDULE(DYNAMIC,1)  
    do ix=1, enl%DF_npix   ! loop over the entries along rows in the output ROI

      iipar(7) = ix
      iipar(8) = iy
      binned = 0.0
      ! for now, assume that the ROI is much smaller than a detector pixel
      pctr = (/ enl%xpc, enl%ypc, enl%L /)
      call CalcEBSDPatternDefect(iipar,cell,Defects,tmLPNH,tmLPSH,prefactor,enl,mcnl,pctr,binned)

! and put the pattern in the correct spot in the batch array   
      patarray(1:enl%numsx, 1:enl%numsy, ix, 1) = binned(1:enl%numsx, 1:enl%numsy)
   end do  ! (ix loop)
!$OMP END DO

! we wait here for all the threads to come together, and then we write the patterns
! for this line to the HDF5 file using thread 0
!$OMP BARRIER
   if (TID.eq.0) then 
     io_int(1) = iy 
     io_int(2) = enl%DF_npiy
     call Message%WriteValue(' -> completed row/total rows ',io_int,2,"(I5,'/',I5)")
     dataset = SC_EBSDpatterns
        dims4 = (/ enl%numsx, enl%numsy, enl%DF_npix, enl%DF_npiy /)
        cnt4 = (/ enl%numsx, enl%numsy, enl%DF_npix, 1 /)
        offset4 = (/ 0, 0, 0, enl%DF_npiy-iy /)
        call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
        if (g_exists) then 
          hdferr = HDF%writeHyperslabFloatArray(dataset, patarray, dims4, offset4, cnt4, overwrite)
        end if
! and store the first pattern in this row in a tiff file (mostly useful for debugging and
! checking the orientation)
      TIFF_filename = 'currentpattern.tiff'
      mi = minval(patarray(:,:,1,1))
      ma = maxval(patarray(:,:,1,1))
      allocate(TIFF_image(enl%numsx,enl%numsy))
      do i=1,enl%numsx 
        do j=1,enl%numsy
          TIFF_image(i,enl%numsy+1-j) = int(255 * (patarray(i,j,1,1)-mi)/(ma-mi))
        end do 
      end do 
      im = image_t(TIFF_image)
      call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
      deallocate(TIFF_image)
   end if

  end do ! end of iy loop
!$OMP END PARALLEL

call HDF%popall()


end associate 

end subroutine EBSDdefect_

!--------------------------------------------------------------------------
recursive subroutine CalcEBSDPatternDefect(ipar,cell,Defects,mLPNH,mLPSH,prefactor,enl,mcnl,patcntr,binned)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcEBSDPatternDefect
!! author: MDG 
!! version: 1.0 
!! date: 06/08/23
!!
!! This routine computes an EBSD pattern for a pixel in an ROI; the detector geometry 
!! is also created here, since we need a different detetcor for each image pixel; so,
!! this routine combines a few different routines into a single one.

use mod_Lambert
use mod_quaternions
use mod_rotations
use mod_defect
use mod_crystallography
use mod_MCfiles
use mod_math

IMPLICIT NONE

integer(kind=irg),INTENT(IN)                    :: ipar(8)
type(Cell_T),INTENT(IN)                         :: cell
type(Defect_T),INTENT(INOUT)                    :: Defects
real(kind=sgl),INTENT(IN)                       :: mLPNH(-ipar(4):ipar(4),-ipar(5):ipar(5),ipar(6))
real(kind=sgl),INTENT(IN)                       :: mLPSH(-ipar(4):ipar(4),-ipar(5):ipar(5),ipar(6))
real(kind=dbl),INTENT(IN)                       :: prefactor
type(EBSDdefectNameListType),INTENT(INOUT)      :: enl
type(MCOpenCLNameListType),INTENT(INOUT)        :: mcnl
real(kind=dbl),INTENT(IN)                       :: patcntr(3)
real(kind=sgl),INTENT(INOUT)                    :: binned(ipar(2),ipar(3))

type(Quaternion_T)                              :: quat, qq, dquat
type(q_T)                                       :: qu
type(o_T)                                       :: om

real(kind=dbl)                                  :: dc(3),dcnew(3), omfc(3,3), omfct(3,3)
real(kind=dbl)                                  :: tFij(3,3,ipar(6)), FF(3,3), FF_inv(3,3)
real(kind=sgl)                                  :: dx,dy,dxm,dym,ixy(2),scl
integer(kind=irg)                               :: ii,jj,kk
integer(kind=irg)                               :: nix,niy,nixp,niyp

real(kind=sgl)                                  :: scin_x, scin_y, rgx, rgy, rgz, z  ! scintillator coordinate ararays [microns]
real(kind=sgl)                                  :: alp, ca, sa, cw, sw, xpos, ypos
real(kind=sgl)                                  :: L2, Ls, Lc, calpha     ! distances
integer(kind=irg)                               :: binx, biny , i, j, istat, k, ipx, ipy, nx, ny, elp   
real(kind=sgl)                                  :: alpha, theta, g, pcvec(3), s, dp           ! direction cosine array
real(kind=sgl)                                  :: sx, rhos, x, bindx, xpc, ypc, L         ! various parameters


! get the quaternion that transforms the foil coordinates into the crystal reference frame
qq = Quaternion_T( qd = (/ 1.D0/sqrt(2.D0),0.D0,0.D0, 1.D0/sqrt(2.D0) /) )
dquat = conjg(Defects%foil%a_fc) 
quat = qq * dquat * conjg(qq)
! also get the transformation matrix which is need to put the deformation tensor into
! the crystal reference frame.
qu = q_T( qdinp = dquat%get_quatd())
om = qu%qo()
omfc = om%o_copyd()
omfct = transpose(omfc)

! ipar(1) = not used 
! ipar(2) = ebsdnl%numsx
! ipar(3) = ebsdnl%numsy
! ipar(4) = ebsdnl%npx
! ipar(5) = ebsdnl%npy
! ipar(6) = number of depth steps
! ipar(7) = i  [image coordinate]
! ipar(8) = j  [image coordinate]

binned = 0.0
scl = float(ipar(4))

!====================================
! ------ generate the detector parameters
!====================================
xpc = sngl(patcntr(1))
ypc = sngl(patcntr(2))
L = sngl(patcntr(3))

! auxiliary angle to rotate between reference frames
alp = 0.5 * cPi - (mcnl%sig - enl%thetac) * dtor
ca = cos(alp)
sa = sin(alp)

cw = cos(mcnl%omega * dtor)
sw = sin(mcnl%omega * dtor)
L2 = L * L

do ii = 1,ipar(2)
    do jj = 1,ipar(3)
! this portion comes from the detector generator routine; the following lines
! determine the direction cosines of the detector pixel with respect to the 
! image pixel (illumination point on the sample)
      scin_x = - ( -xpc - ( 1.0 - ipar(2)) * 0.5 - real(ii) ) * enl%delta
      scin_y =   ( ypc - ( 1.0 - ipar(3) ) * 0.5 - real(jj) ) * enl%delta
      ! sx = L2 + scin_x * scin_x
      Ls = -sw * scin_x + L*cw
      Lc = cw * scin_x + L*sw
      ! rhos = 1.0/sqrt(sx + scin_y**2)  ! the vector will be normalized anyway...
      rgx = (scin_y * ca + sa * Ls) ! * rhos
      rgy = Lc ! * rhos
      rgz = (-sa * scin_y + ca * Ls) ! * rhos
      z = 1.0/sqrt(rgx*rgx+rgy*rgy+rgz*rgz)
      rgx = rgx*z
      rgy = rgy*z
      rgz = rgz*z
! get the pixel direction cosines; these are also the direction cosines that 
! represent the orientation of the integration column in the sample.  We need to
! also transform the deformation tensor into the crystal reference frame before
! carrying out the integration. 
      dc = dble((/ rgx,rgy,rgz /))
      if (rgz.lt.0.0) dc = -dc

! get all the deformation tensors in the sample reference frame 
      call Defects%CalcFcolumn(cell, dc, tFij, ipar(7), ipar(8) )

! loop over the depth instead of the energy; get the transformed deformation tensor at each depth 
! to determine the direction cosines of the sampling unit vector.        
      do kk = 1, ipar(6)
        FF = tFij(:,:,kk)
        FF = matmul(omfc,matmul(FF,omfct))
        call mInvert(FF, FF_inv, .FALSE.)
! apply the grain rotation 
        dcnew = quat%quat_Lp(dc)
! apply the deformation
        dcnew = matmul(FF_inv, dcnew)
! and normalize the direction cosines (to remove any rounding errors)
        dcnew = dcnew/sqrt(sum(dcnew**2))

! convert these direction cosines to interpolation coordinates in the Rosca-Lambert projection
          call LambertgetInterpolation(sngl(dcnew), scl, ipar(4), ipar(5), nix, niy, nixp, niyp, dx, dy, dxm, dym)

! interpolate the intensity
          if (dcnew(3) .ge. 0.0) then
                binned(ii,jj) = binned(ii,jj) + ( mLPNH(nix,niy,kk) * dxm * dym + &
                                                 mLPNH(nixp,niy,kk) * dx * dym + mLPNH(nix,niyp,kk) * dxm * dy + &
                                                 mLPNH(nixp,niyp,kk) * dx * dy )
          else
                binned(ii,jj) = binned(ii,jj) + ( mLPSH(nix,niy,kk) * dxm * dym + &
                                                 mLPSH(nixp,niy,kk) * dx * dym + mLPSH(nix,niyp,kk) * dxm * dy + &
                                                 mLPSH(nixp,niyp,kk) * dx * dy )
          end if
        end do 
    end do
end do

binned = prefactor * binned

end subroutine CalcEBSDPatternDefect


end module mod_EBSDdefect
