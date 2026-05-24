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

module mod_EBSDGlobalOpt
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/06/23
  !!
  !! class definition for the EMEBSDGlobalOpt program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMEBSDGlobalOpt program
type, public :: EBSDGlobalOptNameListType
  integer(kind=irg)        :: NP
  integer(kind=irg)        :: itermax
  integer(kind=irg)        :: strategy 
  integer(kind=irg)        :: refresh
  integer(kind=irg)        :: iwrite
  integer(kind=irg)        :: method(3)
  integer(kind=irg)        :: GrainID
  real(kind=sgl)           :: VTR 
  real(kind=sgl)           :: CR_XC
  real(kind=sgl)           :: F_XC
  real(kind=sgl)           :: F_CR
  real(kind=sgl)           :: bound(3)
  real(kind=sgl)           :: w
  real(kind=sgl)           :: w_damp
  real(kind=sgl)           :: c1 
  real(kind=sgl)           :: c2 
  integer(kind=irg)        :: objective
  character(fnlen)         :: outputfile
  character(fnlen)         :: EBSDnmlfile
  character(1)             :: hybrid
  character(2)             :: globalopt
  character(1)             :: single_opt
  character(1)             :: single_grain
end type EBSDGlobalOptNameListType

! class definition
type, public :: EBSDGlobalOpt_T
private 
  character(fnlen)       :: nmldeffile = 'EMEBSDGlobalOpt.nml'
  type(EBSDGlobalOptNameListType)  :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: EBSDGlobalOpt_
  procedure, pass(self) :: setNP_
  procedure, pass(self) :: getNP_
  procedure, pass(self) :: setitermax_
  procedure, pass(self) :: getitermax_
  procedure, pass(self) :: setstrategy_
  procedure, pass(self) :: getstrategy_
  procedure, pass(self) :: setrefresh_
  procedure, pass(self) :: getrefresh_
  procedure, pass(self) :: setiwrite_
  procedure, pass(self) :: getiwrite_
  procedure, pass(self) :: setmethod_
  procedure, pass(self) :: getmethod_
  procedure, pass(self) :: setGrainID_
  procedure, pass(self) :: getGrainID_
  procedure, pass(self) :: setVTR_
  procedure, pass(self) :: getVTR_
  procedure, pass(self) :: setCR_XC_
  procedure, pass(self) :: getCR_XC_
  procedure, pass(self) :: setF_XC_
  procedure, pass(self) :: getF_XC_
  procedure, pass(self) :: setF_CR_
  procedure, pass(self) :: getF_CR_
  procedure, pass(self) :: setbound_
  procedure, pass(self) :: getbound_
  procedure, pass(self) :: setw_
  procedure, pass(self) :: getw_
  procedure, pass(self) :: setw_damp_
  procedure, pass(self) :: getw_damp_
  procedure, pass(self) :: setc1_
  procedure, pass(self) :: getc1_
  procedure, pass(self) :: setc2_
  procedure, pass(self) :: getc2_
  procedure, pass(self) :: setobjective_
  procedure, pass(self) :: getobjective_
  procedure, pass(self) :: setoutputfile_
  procedure, pass(self) :: getoutputfile_
  procedure, pass(self) :: setEBSDnmlfile_
  procedure, pass(self) :: getEBSDnmlfile_
  procedure, pass(self) :: sethybrid_
  procedure, pass(self) :: gethybrid_
  procedure, pass(self) :: setglobalopt_
  procedure, pass(self) :: getglobalopt_
  procedure, pass(self) :: setsingle_opt_
  procedure, pass(self) :: getsingle_opt_
  procedure, pass(self) :: setsingle_grain_
  procedure, pass(self) :: getsingle_grain_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: EBSDGlobalOpt => EBSDGlobalOpt_
  generic, public :: setNP => setNP_
  generic, public :: getNP => getNP_
  generic, public :: setitermax => setitermax_
  generic, public :: getitermax => getitermax_
  generic, public :: setstrategy => setstrategy_
  generic, public :: getstrategy => getstrategy_
  generic, public :: setrefresh => setrefresh_
  generic, public :: getrefresh => getrefresh_
  generic, public :: setiwrite => setiwrite_
  generic, public :: getiwrite => getiwrite_
  generic, public :: setmethod => setmethod_
  generic, public :: getmethod => getmethod_
  generic, public :: setGrainID => setGrainID_
  generic, public :: getGrainID => getGrainID_
  generic, public :: setVTR => setVTR_
  generic, public :: getVTR => getVTR_
  generic, public :: setCR_XC => setCR_XC_
  generic, public :: getCR_XC => getCR_XC_
  generic, public :: setF_XC => setF_XC_
  generic, public :: getF_XC => getF_XC_
  generic, public :: setF_CR => setF_CR_
  generic, public :: getF_CR => getF_CR_
  generic, public :: setbound => setbound_
  generic, public :: getbound => getbound_
  generic, public :: setw => setw_
  generic, public :: getw => getw_
  generic, public :: setw_damp => setw_damp_
  generic, public :: getw_damp => getw_damp_
  generic, public :: setc1 => setc1_
  generic, public :: getc1 => getc1_
  generic, public :: setc2 => setc2_
  generic, public :: getc2 => getc2_
  generic, public :: setobjective => setobjective_
  generic, public :: getobjective => getobjective_
  generic, public :: setoutputfile => setoutputfile_
  generic, public :: getoutputfile => getoutputfile_
  generic, public :: setEBSDnmlfile => setEBSDnmlfile_
  generic, public :: getEBSDnmlfile => getEBSDnmlfile_
  generic, public :: sethybrid => sethybrid_
  generic, public :: gethybrid => gethybrid_
  generic, public :: setglobalopt => setglobalopt_
  generic, public :: getglobalopt => getglobalopt_
  generic, public :: setsingle_opt => setsingle_opt_
  generic, public :: getsingle_opt => getsingle_opt_
  generic, public :: setsingle_grain => setsingle_grain_
  generic, public :: getsingle_grain => getsingle_grain_

end type EBSDGlobalOpt_T

! the constructor routine for this class 
interface EBSDGlobalOpt_T
  module procedure EBSDGlobalOpt_constructor
end interface EBSDGlobalOpt_T

contains

!--------------------------------------------------------------------------
type(EBSDGlobalOpt_T) function EBSDGlobalOpt_constructor( nmlfile ) result(EBSDGlobalOpt)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDGlobalOpt_constructor
!! author: MDG 
!! version: 1.0 
!! date: 12/06/23
!!
!! constructor for the EBSDGlobalOpt_T Class; optionally reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

if (present(nmlfile)) then 
  call EBSDGlobalOpt%readNameList(nmlfile)
end if

end function EBSDGlobalOpt_constructor

!--------------------------------------------------------------------------
subroutine EBSDGlobalOpt_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 12/06/23
!!
!! destructor for the EBSDGlobalOpt_T Class
 
IMPLICIT NONE

type(EBSDGlobalOpt_T), INTENT(INOUT)  :: self 

call reportDestructor('EBSDGlobalOpt_T')

end subroutine EBSDGlobalOpt_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 12/06/23
!!
!! read the namelist from an nml file for the EBSDGlobalOpt_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
character(fnlen),INTENT(IN)           :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)           :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                        :: EMsoft 
type(IO_T)                            :: Message       
logical                               :: skipread = .FALSE.

integer(kind=irg)                     :: NP
integer(kind=irg)                     :: itermax
integer(kind=irg)                     :: strategy 
integer(kind=irg)                     :: refresh
integer(kind=irg)                     :: iwrite
integer(kind=irg)                     :: method(3)
integer(kind=irg)                     :: GrainID
real(kind=sgl)                        :: VTR 
real(kind=sgl)                        :: CR_XC
real(kind=sgl)                        :: F_XC
real(kind=sgl)                        :: F_CR
real(kind=sgl)                        :: bound(3)
real(kind=sgl)                        :: w
real(kind=sgl)                        :: w_damp
real(kind=sgl)                        :: c1 
real(kind=sgl)                        :: c2 
integer(kind=irg)                     :: objective
character(fnlen)                      :: outputfile
character(fnlen)                      :: EBSDnmlfile
character(1)                          :: hybrid
character(2)                          :: globalopt
character(1)                          :: single_opt
character(1)                          :: single_grain

! define the IO namelist to facilitate passing variables to the program.
namelist  / EBSDDEdata / NP, itermax, strategy, refresh, iwrite, method, VTR, CR_XC, F_XC, F_CR, bound, hybrid, globalopt, &
                         w, w_damp, c1, c2, GrainID, objective, outputfile, EBSDnmlfile, single_opt, single_grain

! set the input parameters to default values 
NP            = 60
itermax       = 100 
strategy      = 2
refresh       = 10 
! iwrite        = 
! method(3)     =
GrainID       = 0
VTR           = -1
CR_XC         = 0.9
F_XC          = 0.5
F_CR          = 0.5
bound         = (/ 0.001, 2.0, 2.0 /)
w             = 1.0
w_damp        = 0.9
c1            = 2.0 
c2            = 2.0 
objective     = 1  
outputfile    = 'undefined'
EBSDnmlfile   = 'undefined'
hybrid        = 'n'
globalopt     = 'DE'
! single_opt    = 
! single_grain  =

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=EBSDDEdata)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(outputfile).eq.'undefined') then
  call Message%printError('readNameList:',' output file name is undefined in '//nmlfile)
 end if

 if (trim(EBSDnmlfile).eq.'undefined') then
  call Message%printError('readNameList:',' EBSDnmlfile file name is undefined in '//nmlfile)
 end if
end if 

self%nml%NP = NP
self%nml%itermax = itermax 
self%nml%strategy = strategy
self%nml%refresh = refresh
self%nml%iwrite = iwrite
self%nml%method = method
self%nml%GrainID = GrainID
self%nml%VTR = VTR
self%nml%CR_XC = CR_XC
self%nml%F_XC = F_XC
self%nml%F_CR = F_CR
self%nml%bound = bound
self%nml%w = w
self%nml%w_damp = w_damp
self%nml%c1 = c1
self%nml%c2 = c2
self%nml%objective = objective
self%nml%outputfile = outputfile
self%nml%EBSDnmlfile = EBSDnmlfile
self%nml%hybrid = hybrid
self%nml%globalopt = globalopt
self%nml%single_opt = single_opt
self%nml%single_grain = single_grain

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 12/06/23
!!
!! pass the namelist for the EBSDGlobalOpt_T Class to the calling program

IMPLICIT NONE 

class(EBSDGlobalOpt_T), INTENT(INOUT)          :: self
type(EBSDGlobalOptNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 12/06/23
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT)   :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 7, n_real = 8
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( enl => self%nml )

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! write all the single integers
io_int = (/ enl%NP, enl%itermax, enl%strategy, enl%refresh, enl%iwrite, enl%GrainID, enl%objective /)
intlist(1) = 'NP'
intlist(2) = 'itermax'
intlist(3) = 'strategy'
intlist(4) = 'refresh'
intlist(5) = 'iwrite'
intlist(6) = 'GrainID'
intlist(7) = 'objective'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ enl%VTR, enl%CR_XC, enl%F_XC, enl%F_CR, enl%w, enl%w_damp, enl%c1, enl%c2 /)
reallist(1) = 'VTR'
reallist(2) = 'CR'
reallist(3) = 'F_XC'
reallist(4) = 'F_CR'
reallist(5) = 'w'
reallist(6) = 'w_damp'
reallist(7) = 'c1'
reallist(8) = 'c2'
call HDF%writeNMLreals(io_real, reallist, n_real)

! a 3-vector
dataset = 'method'
hdferr = HDF%writeDatasetIntegerArray(dataset, enl%method, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create method dataset', hdferr)

! a 4-vector
dataset = 'bound'
hdferr = HDF%writeDatasetFloatArray(dataset, enl%bound, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create bound dataset', hdferr)

! write all the strings
dataset = SC_outputfile
line2(1) = trim(enl%outputfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create outputfile dataset', hdferr)

dataset = 'EBSDnmlfile'
line2(1) = trim(enl%EBSDnmlfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create maskpattEBSDnmlfileern dataset', hdferr)

dataset = 'hybrid'
line2(1) = trim(enl%hybrid)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create hybrid dataset', hdferr)

dataset = 'globalopt'
line2(1) = trim(enl%globalopt)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create globalopt dataset', hdferr)

dataset = 'single_opt'
line2(1) = trim(enl%single_opt)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create single_opt dataset', hdferr)

dataset = 'single_grain'
line2(1) = trim(enl%single_grain)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create single_grain dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine setNP_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setNP_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set NP in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%NP = inp

end subroutine setNP_

!--------------------------------------------------------------------------
function getNP_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getNP_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get NP from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
integer(kind=irg)                     :: out

out = self%nml%NP

end function getNP_

!--------------------------------------------------------------------------
subroutine setitermax_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setitermax_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set itermax in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%itermax = inp

end subroutine setitermax_

!--------------------------------------------------------------------------
function getitermax_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getitermax_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get itermax from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
integer(kind=irg)                     :: out

out = self%nml%itermax

end function getitermax_

!--------------------------------------------------------------------------
subroutine setstrategy_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setstrategy_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set strategy in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%strategy = inp

end subroutine setstrategy_

!--------------------------------------------------------------------------
function getstrategy_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getstrategy_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get strategy from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
integer(kind=irg)                     :: out

out = self%nml%strategy

end function getstrategy_

!--------------------------------------------------------------------------
subroutine setrefresh_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setrefresh_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set refresh in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%refresh = inp

end subroutine setrefresh_

!--------------------------------------------------------------------------
function getrefresh_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getrefresh_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get refresh from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
integer(kind=irg)                     :: out

out = self%nml%refresh

end function getrefresh_

!--------------------------------------------------------------------------
subroutine setiwrite_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setiwrite_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set iwrite in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%iwrite = inp

end subroutine setiwrite_

!--------------------------------------------------------------------------
function getiwrite_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getiwrite_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get iwrite from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
integer(kind=irg)                     :: out

out = self%nml%iwrite

end function getiwrite_

!--------------------------------------------------------------------------
subroutine setmethod_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setmethod_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set method in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)             :: inp(3)

self%nml%method = inp

end subroutine setmethod_

!--------------------------------------------------------------------------
function getmethod_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getmethod_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get method from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT)     :: self
integer(kind=irg)                         :: out(3)

out = self%nml%method

end function getmethod_

!--------------------------------------------------------------------------
subroutine setGrainID_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setGrainID_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set GrainID in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%GrainID = inp

end subroutine setGrainID_

!--------------------------------------------------------------------------
function getGrainID_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getGrainID_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get GrainID from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
integer(kind=irg)                     :: out

out = self%nml%GrainID

end function getGrainID_

!--------------------------------------------------------------------------
subroutine setVTR_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setVTR_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set VTR in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
real(kind=sgl), INTENT(IN)            :: inp

self%nml%VTR = inp

end subroutine setVTR_

!--------------------------------------------------------------------------
function getVTR_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getVTR_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get VTR from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
real(kind=sgl)                        :: out

out = self%nml%VTR

end function getVTR_

!--------------------------------------------------------------------------
subroutine setCR_XC_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setCR_XC_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set CR_XC in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
real(kind=sgl), INTENT(IN)            :: inp

self%nml%CR_XC = inp

end subroutine setCR_XC_

!--------------------------------------------------------------------------
function getCR_XC_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getCR_XC_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get CR_XC from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
real(kind=sgl)                        :: out

out = self%nml%CR_XC

end function getCR_XC_

!--------------------------------------------------------------------------
subroutine setF_XC_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setF_XC_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set F_XC in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
real(kind=sgl), INTENT(IN)            :: inp

self%nml%F_XC = inp

end subroutine setF_XC_

!--------------------------------------------------------------------------
function getF_XC_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getF_XC_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get F_XC from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
real(kind=sgl)                        :: out

out = self%nml%F_XC

end function getF_XC_

!--------------------------------------------------------------------------
subroutine setF_CR_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setF_CR_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set F_CR in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
real(kind=sgl), INTENT(IN)            :: inp

self%nml%F_CR = inp

end subroutine setF_CR_

!--------------------------------------------------------------------------
function getF_CR_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getF_CR_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get F_CR from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
real(kind=sgl)                        :: out

out = self%nml%F_CR

end function getF_CR_

!--------------------------------------------------------------------------
subroutine setbound_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setbound_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set bound in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)                :: inp(3)

self%nml%bound = inp

end subroutine setbound_

!--------------------------------------------------------------------------
function getbound_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getbound_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get bound from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT)     :: self
real(kind=sgl)                            :: out(3)

out = self%nml%bound

end function getbound_

!--------------------------------------------------------------------------
subroutine setw_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setw_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set w in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
real(kind=sgl), INTENT(IN)            :: inp

self%nml%w = inp

end subroutine setw_

!--------------------------------------------------------------------------
function getw_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getw_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get w from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
real(kind=sgl)                        :: out

out = self%nml%w

end function getw_

!--------------------------------------------------------------------------
subroutine setw_damp_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setw_damp_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set w_damp in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
real(kind=sgl), INTENT(IN)            :: inp

self%nml%w_damp = inp

end subroutine setw_damp_

!--------------------------------------------------------------------------
function getw_damp_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getw_damp_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get w_damp from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
real(kind=sgl)                        :: out

out = self%nml%w_damp

end function getw_damp_

!--------------------------------------------------------------------------
subroutine setc1_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setc1_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set c1 in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
real(kind=sgl), INTENT(IN)            :: inp

self%nml%c1 = inp

end subroutine setc1_

!--------------------------------------------------------------------------
function getc1_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getc1_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get c1 from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
real(kind=sgl)                        :: out

out = self%nml%c1

end function getc1_

!--------------------------------------------------------------------------
subroutine setc2_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setc2_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set c2 in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
real(kind=sgl), INTENT(IN)            :: inp

self%nml%c2 = inp

end subroutine setc2_

!--------------------------------------------------------------------------
function getc2_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getc2_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get c2 from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
real(kind=sgl)                        :: out

out = self%nml%c2

end function getc2_

!--------------------------------------------------------------------------
subroutine setobjective_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setobjective_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set objective in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%objective = inp

end subroutine setobjective_

!--------------------------------------------------------------------------
function getobjective_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getobjective_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get objective from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
integer(kind=irg)                     :: out

out = self%nml%objective

end function getobjective_

!--------------------------------------------------------------------------
subroutine setoutputfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setoutputfile_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set outputfile in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%outputfile = trim(inp)

end subroutine setoutputfile_

!--------------------------------------------------------------------------
function getoutputfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getoutputfile_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get outputfile from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
character(fnlen)                      :: out

out = trim(self%nml%outputfile)

end function getoutputfile_

!--------------------------------------------------------------------------
subroutine setEBSDnmlfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setEBSDnmlfile_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set EBSDnmlfile in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%EBSDnmlfile = trim(inp)

end subroutine setEBSDnmlfile_

!--------------------------------------------------------------------------
function getEBSDnmlfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getEBSDnmlfile_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get EBSDnmlfile from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
character(fnlen)                      :: out

out = trim(self%nml%EBSDnmlfile)

end function getEBSDnmlfile_

!--------------------------------------------------------------------------
subroutine sethybrid_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: sethybrid_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set hybrid in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
character(1), INTENT(IN)              :: inp

self%nml%hybrid = trim(inp)

end subroutine sethybrid_

!--------------------------------------------------------------------------
function gethybrid_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: gethybrid_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get hybrid from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
character(1)                          :: out

out = trim(self%nml%hybrid)

end function gethybrid_

!--------------------------------------------------------------------------
subroutine setglobalopt_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setglobalopt_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set globalopt in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
character(2), INTENT(IN)              :: inp

self%nml%globalopt = trim(inp)

end subroutine setglobalopt_

!--------------------------------------------------------------------------
function getglobalopt_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getglobalopt_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get globalopt from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
character(2)                          :: out

out = trim(self%nml%globalopt)

end function getglobalopt_

!--------------------------------------------------------------------------
subroutine setsingle_opt_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setsingle_opt_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set single_opt in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
character(1), INTENT(IN)              :: inp

self%nml%single_opt = trim(inp)

end subroutine setsingle_opt_

!--------------------------------------------------------------------------
function getsingle_opt_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getsingle_opt_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get single_opt from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
character(1)                          :: out

out = trim(self%nml%single_opt)

end function getsingle_opt_

!--------------------------------------------------------------------------
subroutine setsingle_grain_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setsingle_grain_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! set single_grain in the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
character(1), INTENT(IN)              :: inp

self%nml%single_grain = trim(inp)

end subroutine setsingle_grain_

!--------------------------------------------------------------------------
function getsingle_grain_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getsingle_grain_
!! author: MDG
!! version: 1.0
!! date: 12/06/23
!!
!! get single_grain from the EBSDGlobalOpt_T class

IMPLICIT NONE

class(EBSDGlobalOpt_T), INTENT(INOUT) :: self
character(1)                          :: out

out = trim(self%nml%single_grain)

end function getsingle_grain_

!--------------------------------------------------------------------------
subroutine EBSDGlobalOpt_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDGlobalOpt_
!! author: MDG 
!! version: 1.0 
!! date: 12/06/23
!!
!! perform the computations

use mod_EMsoft
use mod_HDFnames

IMPLICIT NONE 

class(EBSDGlobalOpt_T), INTENT(INOUT)   :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 
type(HDFnames_T), INTENT(INOUT)         :: HDFnames



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

mem = Memory_T()















end subroutine EBSDGlobalOpt_



end module mod_EBSDGlobalOpt