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
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 07/22/25
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(GBOanalysis_T), INTENT(INOUT)        :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 11, n_real = 9
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( mcnl => self%nml )

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine GBOanalysis_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: GBOanalysis_
!! author: MDG 
!! version: 1.0 
!! date: 07/22/25
!!
!! perform the computations

use mod_EMsoft
use mod_HDFnames

IMPLICIT NONE 

class(GBOanalysis_T), INTENT(INOUT)       :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

end subroutine GBOanalysis_



end module mod_GBOanalysis