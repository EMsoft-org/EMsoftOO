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

module mod_DIpreview
  !! author: MDG 
  !! version: 1.0 
  !! date: 04/08/20
  !!
  !! class definition for the EMDIpreview program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMDIpreview program
type, public :: DIpreviewNameListType
  integer(kind=irg) :: numsx
  integer(kind=irg) :: numsy
  integer(kind=irg) :: hipasswnsteps
  integer(kind=irg) :: nregionsmin
  integer(kind=irg) :: nregionsmax
  integer(kind=irg) :: nregionsstepsize
  integer(kind=irg) :: patx
  integer(kind=irg) :: paty
  integer(kind=irg) :: ipf_wd
  integer(kind=irg) :: ipf_ht
  integer(kind=irg) :: numav
  real(kind=sgl)    :: hipasswmax
  character(fnlen)  :: patternfile
  character(fnlen)  :: tifffile
  character(fnlen)  :: exptfile
  character(fnlen)  :: inputtype
  character(fnlen)  :: HDFstrings(10)
end type DIpreviewNameListType

! class definition
type, public :: DIpreview_T
private 
  character(fnlen)       :: nmldeffile = 'EMDIpreview.nml'
  type(DIpreviewNameListType)  :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: DIpreview_
  procedure, pass(self) :: get_numsx_
  procedure, pass(self) :: get_numsy_
  procedure, pass(self) :: get_hipasswnsteps_
  procedure, pass(self) :: get_nregionsmin_
  procedure, pass(self) :: get_nregionsmax_
  procedure, pass(self) :: get_nregionsstepsize_
  procedure, pass(self) :: get_patx_
  procedure, pass(self) :: get_paty_
  procedure, pass(self) :: get_ipf_wd_
  procedure, pass(self) :: get_ipf_ht_
  procedure, pass(self) :: get_numav_
  procedure, pass(self) :: get_hipasswmax_
  procedure, pass(self) :: get_patternfile_
  procedure, pass(self) :: get_tifffile_
  procedure, pass(self) :: get_exptfile_
  procedure, pass(self) :: get_inputtype_
  procedure, pass(self) :: get_HDFstrings_
  procedure, pass(self) :: set_numsx_
  procedure, pass(self) :: set_numsy_
  procedure, pass(self) :: set_hipasswnsteps_
  procedure, pass(self) :: set_nregionsmin_
  procedure, pass(self) :: set_nregionsmax_
  procedure, pass(self) :: set_nregionsstepsize_
  procedure, pass(self) :: set_patx_
  procedure, pass(self) :: set_paty_
  procedure, pass(self) :: set_ipf_wd_
  procedure, pass(self) :: set_ipf_ht_
  procedure, pass(self) :: set_numav_
  procedure, pass(self) :: set_hipasswmax_
  procedure, pass(self) :: set_patternfile_
  procedure, pass(self) :: set_tifffile_
  procedure, pass(self) :: set_exptfile_
  procedure, pass(self) :: set_inputtype_
  procedure, pass(self) :: set_HDFstrings_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: DIpreview => DIpreview_
  generic, public :: get_numsx => get_numsx_
  generic, public :: get_numsy => get_numsy_
  generic, public :: get_hipasswnsteps => get_hipasswnsteps_
  generic, public :: get_nregionsmin => get_nregionsmin_
  generic, public :: get_nregionsmax => get_nregionsmax_
  generic, public :: get_nregionsstepsize => get_nregionsstepsize_
  generic, public :: get_patx => get_patx_
  generic, public :: get_paty => get_paty_
  generic, public :: get_ipf_wd => get_ipf_wd_
  generic, public :: get_ipf_ht => get_ipf_ht_
  generic, public :: get_numav => get_numav_
  generic, public :: get_hipasswmax => get_hipasswmax_
  generic, public :: get_patternfile => get_patternfile_
  generic, public :: get_tifffile => get_tifffile_
  generic, public :: get_exptfile => get_exptfile_
  generic, public :: get_inputtype => get_inputtype_
  generic, public :: get_HDFstrings => get_HDFstrings_
  generic, public :: set_numsx => set_numsx_
  generic, public :: set_numsy => set_numsy_
  generic, public :: set_hipasswnsteps => set_hipasswnsteps_
  generic, public :: set_nregionsmin => set_nregionsmin_
  generic, public :: set_nregionsmax => set_nregionsmax_
  generic, public :: set_nregionsstepsize => set_nregionsstepsize_
  generic, public :: set_patx => set_patx_
  generic, public :: set_paty => set_paty_
  generic, public :: set_ipf_wd => set_ipf_wd_
  generic, public :: set_ipf_ht => set_ipf_ht_
  generic, public :: set_numav => set_numav_
  generic, public :: set_hipasswmax => set_hipasswmax_
  generic, public :: set_patternfile => set_patternfile_
  generic, public :: set_tifffile => set_tifffile_
  generic, public :: set_exptfile => set_exptfile_
  generic, public :: set_inputtype => set_inputtype_
  generic, public :: set_HDFstrings => set_HDFstrings_
 
end type DIpreview_T

!DEC$ ATTRIBUTES DLLEXPORT :: getNameList
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList
!DEC$ ATTRIBUTES DLLEXPORT :: DIpreview
!DEC$ ATTRIBUTES DLLEXPORT :: get_numsx
!DEC$ ATTRIBUTES DLLEXPORT :: set_numsx
!DEC$ ATTRIBUTES DLLEXPORT :: get_numsy
!DEC$ ATTRIBUTES DLLEXPORT :: set_numsy
!DEC$ ATTRIBUTES DLLEXPORT :: get_hipasswnsteps
!DEC$ ATTRIBUTES DLLEXPORT :: set_hipasswnsteps
!DEC$ ATTRIBUTES DLLEXPORT :: get_nregionsmin
!DEC$ ATTRIBUTES DLLEXPORT :: set_nregionsmin
!DEC$ ATTRIBUTES DLLEXPORT :: get_nregionsmax
!DEC$ ATTRIBUTES DLLEXPORT :: set_nregionsmax
!DEC$ ATTRIBUTES DLLEXPORT :: get_nregionsstepsize
!DEC$ ATTRIBUTES DLLEXPORT :: set_nregionsstepsize
!DEC$ ATTRIBUTES DLLEXPORT :: get_patx
!DEC$ ATTRIBUTES DLLEXPORT :: set_patx
!DEC$ ATTRIBUTES DLLEXPORT :: get_paty
!DEC$ ATTRIBUTES DLLEXPORT :: set_paty
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_wd
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_wd
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_ht
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_ht
!DEC$ ATTRIBUTES DLLEXPORT :: get_numav
!DEC$ ATTRIBUTES DLLEXPORT :: set_numav
!DEC$ ATTRIBUTES DLLEXPORT :: get_hipasswmax
!DEC$ ATTRIBUTES DLLEXPORT :: set_hipasswmax
!DEC$ ATTRIBUTES DLLEXPORT :: get_patternfile
!DEC$ ATTRIBUTES DLLEXPORT :: set_patternfile
!DEC$ ATTRIBUTES DLLEXPORT :: get_tifffile
!DEC$ ATTRIBUTES DLLEXPORT :: set_tifffile
!DEC$ ATTRIBUTES DLLEXPORT :: get_exptfile
!DEC$ ATTRIBUTES DLLEXPORT :: set_exptfile
!DEC$ ATTRIBUTES DLLEXPORT :: get_inputtype
!DEC$ ATTRIBUTES DLLEXPORT :: set_inputtype
!DEC$ ATTRIBUTES DLLEXPORT :: get_HDFstrings
!DEC$ ATTRIBUTES DLLEXPORT :: set_HDFstrings

! the constructor routine for this class 
interface DIpreview_T
  module procedure DIpreview_constructor
end interface DIpreview_T

contains

!--------------------------------------------------------------------------
type(DIpreview_T) function DIpreview_constructor( nmlfile ) result(DIpreview)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! constructor for the DIpreview_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call DIpreview%readNameList(nmlfile)

end function DIpreview_constructor

!--------------------------------------------------------------------------
subroutine DIpreview_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! destructor for the DIpreview_T Class
 
IMPLICIT NONE

type(DIpreview_T), INTENT(INOUT)  :: self 

call reportDestructor('DIpreview_T')

end subroutine DIpreview_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! read the namelist from an nml file for the DIpreview_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.


integer(kind=irg) :: numsx
integer(kind=irg) :: numsy
integer(kind=irg) :: hipasswnsteps
integer(kind=irg) :: nregionsmin
integer(kind=irg) :: nregionsmax
integer(kind=irg) :: nregionsstepsize
integer(kind=irg) :: patx
integer(kind=irg) :: paty
integer(kind=irg) :: ipf_wd
integer(kind=irg) :: ipf_ht
integer(kind=irg) :: numav
real(kind=sgl)    :: hipasswmax
character(fnlen)  :: patternfile
character(fnlen)  :: tifffile
character(fnlen)  :: exptfile
character(fnlen)  :: inputtype
character(fnlen)  :: HDFstrings(10)

namelist / EBSDDIpreviewdata / numsx, numsy, hipasswmax, hipasswnsteps, nregionsstepsize, &
          nregionsmax, nregionsmin, patx, paty, tifffile, exptfile, inputtype, HDFstrings, ipf_wd, &
          ipf_ht, patternfile, numav

! set the input parameters to default values
numsx = 0
numsy = 0
hipasswmax = 0.5
hipasswnsteps = 10
nregionsmin = 1
nregionsmax = 10
nregionsstepsize = 1
patx = 1
paty = 1
ipf_wd = 100
ipf_ht = 100
numav = 0
patternfile = 'undefined'
tifffile = 'undefined'
exptfile = 'undefined'
inputtype = 'Binary'
HDFstrings = ''

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=EBSDDIpreviewdata)
    close(UNIT=dataunit,STATUS='keep')

! check for required entries
        
    if (trim(exptfile).eq.'undefined') then
        call Message%printError('readNameList:',' experimental file name is undefined in '//nmlfile)
    end if

    if (trim(tifffile).eq.'undefined') then
        call Message%printError('readNameList:',' TIFF file name is undefined in '//nmlfile)
    end if

    if (numsx.eq.0) then 
        call Message%printError('readNameList:',' pattern size numsx is zero in '//nmlfile)
    end if

    if (numsy.eq.0) then 
        call Message%printError('readNameList:',' pattern size numsy is zero in '//nmlfile)
    end if
end if

self%nml%numsx = numsx
self%nml%numsy = numsy
self%nml%hipasswnsteps = hipasswnsteps
self%nml%nregionsmin = nregionsmin
self%nml%nregionsmax = nregionsmax
self%nml%nregionsstepsize = nregionsstepsize
self%nml%patx = patx
self%nml%paty = paty
self%nml%ipf_wd = ipf_wd
self%nml%ipf_ht = ipf_ht
self%nml%numav = numav
self%nml%patternfile = patternfile
self%nml%hipasswmax = hipasswmax
self%nml%tifffile = tifffile
self%nml%exptfile = exptfile
self%nml%inputtype = inputtype
self%nml%HDFstrings = HDFstrings

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! pass the namelist for the DIpreview_T Class to the calling program

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)          :: self
type(DIpreviewNameListType)                :: nml

nml = self%nml

end function getNameList_


!--------------------------------------------------------------------------
function get_numsx_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! get numsx from the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg)                     :: out

out = self%nml%numsx

end function get_numsx_

!--------------------------------------------------------------------------
subroutine set_numsx_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! set numsx in the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%numsx = inp

end subroutine set_numsx_

!--------------------------------------------------------------------------
function get_numsy_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! get numsy from the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg)                     :: out

out = self%nml%numsy

end function get_numsy_

!--------------------------------------------------------------------------
subroutine set_numsy_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! set numsy in the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%numsy = inp

end subroutine set_numsy_

!--------------------------------------------------------------------------
function get_hipasswnsteps_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! get hipasswnsteps from the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg)                     :: out

out = self%nml%hipasswnsteps

end function get_hipasswnsteps_

!--------------------------------------------------------------------------
subroutine set_hipasswnsteps_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! set hipasswnsteps in the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%hipasswnsteps = inp

end subroutine set_hipasswnsteps_

!--------------------------------------------------------------------------
function get_nregionsmin_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! get nregionsmin from the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg)                     :: out

out = self%nml%nregionsmin

end function get_nregionsmin_

!--------------------------------------------------------------------------
subroutine set_nregionsmin_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! set nregionsmin in the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%nregionsmin = inp

end subroutine set_nregionsmin_

!--------------------------------------------------------------------------
function get_nregionsmax_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! get nregionsmax from the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg)                     :: out

out = self%nml%nregionsmax

end function get_nregionsmax_

!--------------------------------------------------------------------------
subroutine set_nregionsmax_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! set nregionsmax in the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%nregionsmax = inp

end subroutine set_nregionsmax_

!--------------------------------------------------------------------------
function get_nregionsstepsize_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! get nregionsstepsize from the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg)                     :: out

out = self%nml%nregionsstepsize

end function get_nregionsstepsize_

!--------------------------------------------------------------------------
subroutine set_nregionsstepsize_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! set nregionsstepsize in the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%nregionsstepsize = inp

end subroutine set_nregionsstepsize_

!--------------------------------------------------------------------------
function get_patx_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! get patx from the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg)                     :: out

out = self%nml%patx

end function get_patx_

!--------------------------------------------------------------------------
subroutine set_patx_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! set patx in the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%patx = inp

end subroutine set_patx_

!--------------------------------------------------------------------------
function get_paty_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! get paty from the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg)                     :: out

out = self%nml%paty

end function get_paty_

!--------------------------------------------------------------------------
subroutine set_paty_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! set paty in the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%paty = inp

end subroutine set_paty_

!--------------------------------------------------------------------------
function get_ipf_wd_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! get ipf_wd from the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg)                     :: out

out = self%nml%ipf_wd

end function get_ipf_wd_

!--------------------------------------------------------------------------
subroutine set_ipf_wd_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! set ipf_wd in the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%ipf_wd = inp

end subroutine set_ipf_wd_

!--------------------------------------------------------------------------
function get_ipf_ht_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! get ipf_ht from the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg)                     :: out

out = self%nml%ipf_ht

end function get_ipf_ht_

!--------------------------------------------------------------------------
subroutine set_ipf_ht_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! set ipf_ht in the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%ipf_ht = inp

end subroutine set_ipf_ht_

!--------------------------------------------------------------------------
function get_numav_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! get numav from the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg)                     :: out

out = self%nml%numav

end function get_numav_

!--------------------------------------------------------------------------
subroutine set_numav_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! set numav in the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%numav = inp

end subroutine set_numav_

!--------------------------------------------------------------------------
function get_hipasswmax_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! get hipasswmax from the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
real(kind=sgl)                        :: out

out = self%nml%hipasswmax

end function get_hipasswmax_

!--------------------------------------------------------------------------
subroutine set_hipasswmax_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! set hipasswmax in the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)            :: inp

self%nml%hipasswmax = inp

end subroutine set_hipasswmax_

!--------------------------------------------------------------------------
function get_patternfile_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! get patternfile from the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
character(fnlen)                      :: out

out = self%nml%patternfile

end function get_patternfile_

!--------------------------------------------------------------------------
subroutine set_patternfile_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! set patternfile in the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%patternfile = inp

end subroutine set_patternfile_

!--------------------------------------------------------------------------
function get_tifffile_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! get tifffile from the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
character(fnlen)                      :: out

out = self%nml%tifffile

end function get_tifffile_

!--------------------------------------------------------------------------
subroutine set_tifffile_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! set tifffile in the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%tifffile = inp

end subroutine set_tifffile_

!--------------------------------------------------------------------------
function get_exptfile_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! get exptfile from the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
character(fnlen)                      :: out

out = self%nml%exptfile

end function get_exptfile_

!--------------------------------------------------------------------------
subroutine set_exptfile_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! set exptfile in the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%exptfile = inp

end subroutine set_exptfile_

!--------------------------------------------------------------------------
function get_inputtype_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! get inputtype from the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
character(fnlen)                      :: out

out = self%nml%inputtype

end function get_inputtype_

!--------------------------------------------------------------------------
subroutine set_inputtype_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! set inputtype in the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%inputtype = inp

end subroutine set_inputtype_

!--------------------------------------------------------------------------
function get_HDFstrings_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! get HDFstrings from the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
character(fnlen)                      :: out(10)

out = self%nml%HDFstrings

end function get_HDFstrings_

!--------------------------------------------------------------------------
subroutine set_HDFstrings_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! set HDFstrings in the DIpreview_T class

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)          :: inp(10)

self%nml%HDFstrings = inp

end subroutine set_HDFstrings_

!--------------------------------------------------------------------------
subroutine DIpreview_(self, EMsoft, progname)
!! author: MDG 
!! version: 1.0 
!! date: 04/08/20
!!
!! perform the computations

use mod_EMsoft
use mod_image
use, intrinsic :: iso_fortran_env
use mod_io
use mod_filters
use mod_vendors
use HDF5
use mod_HDFsupport

IMPLICIT NONE 

class(DIpreview_T), INTENT(INOUT)       :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 

type(Vendor_T)                          :: VT 
type(IO_T)                              :: Message
type(HDF_T)                             :: HDF 

character(fnlen)                        :: ename, image_filename, fname
integer(kind=irg)                       :: iunitexpt, recordsize, ierr, kk, ii, jj, i, j, numr, numw, binx, biny, &
                                           xoffset, yoffset, io_int(2), istat, L, patsz , hdferr, nx, ny
integer(HSIZE_T)                        :: dims3(3), offset3(3)
logical                                 :: f_exists
real(kind=sgl)                          :: mi, ma, io_real(1)
real(kind=dbl)                          :: x, y, val, v2
real(kind=sgl),allocatable              :: expt(:), pattern(:,:), pcopy(:,:), hpvals(:), sumexpt(:)
integer(kind=irg),allocatable           :: nrvals(:), pint(:,:), ppp(:,:)
type(C_PTR)                             :: HPplanf, HPplanb
complex(kind=dbl),allocatable           :: hpmask(:,:)
complex(C_DOUBLE_COMPLEX),allocatable   :: inp(:,:), outp(:,:)
real(kind=dbl),allocatable              :: rrdata(:,:), ffdata(:,:), ksqarray(:,:)

! declare variables for use in object oriented image module
integer                                 :: iostat
character(len=128)                      :: iomsg
logical                                 :: isInteger
type(image_t)                           :: im, im2
integer(int8)                           :: i8 (3,4), int8val
integer(int8), allocatable              :: output_image(:,:)


! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()

associate(enl=>self%nml)

binx = enl%numsx
biny = enl%numsy
L = binx * biny
recordsize = 4 * L
patsz = L

! open the file with experimental patterns; depending on the inputtype parameter, this
! can be a regular binary file, as produced by a MatLab or IDL script (default); a 
! pattern file produced by EMEBSD.f90; or a vendor binary or HDF5 file... in each case we need to 
! open the file and leave it open, then use the getSingleExpPattern() routine to read a 
! pattern into the expt variable ...  at the end, we use closeExpPatternFile() to
! properly close the experimental pattern file
VT = Vendor_T(enl%inputtype)
fname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(enl%exptfile)
call VT%set_filename(fname)
istat = VT%openExpPatternFile(EMsoft, enl%ipf_wd, L, recordsize, enl%HDFstrings, HDF)
if (istat.ne.0) then
    call Message%printError("DIpreview:", "Fatal error handling experimental pattern file")
end if

! should we average patterns locally ?
allocate(expt(patsz))
dims3 = (/ binx, biny, 1 /)
if (enl%numav.ge.0) then
  io_int(1) = 2*enl%numav+1
  io_int(2) = 2*enl%numav+1
  call Message%WriteValue(' Averaging patterns over ', io_int, 2, "(I3,' by ',I3,' area')")
  allocate(sumexpt(patsz))
  sumexpt = 0.0
  jj = 0 
  do i=-enl%numav,enl%numav
    if ((enl%patx+i.gt.0).and.(enl%patx+i.lt.enl%ipf_wd)) then
      do j=-enl%numav,enl%numav
        if ((enl%paty+j.gt.0).and.(enl%paty+j.lt.enl%ipf_ht)) then
          offset3 = (/ 0, 0, (enl%paty+j) * enl%ipf_wd + (enl%patx+i) /)
          call VT%getSingleExpPattern(enl%paty, enl%ipf_wd, patsz, L, dims3, offset3, expt, enl%HDFstrings, HDF)
          sumexpt = sumexpt + expt
          jj = jj+1
        end if 
      end do 
    end if 
  end do
  sumexpt = sumexpt / float(jj)
end if

! and read the center pattern (again)
offset3 = (/ 0, 0, enl%paty * enl%ipf_wd + enl%patx /)
call VT%getSingleExpPattern(enl%paty, enl%ipf_wd, patsz, L, dims3, offset3, expt, enl%HDFstrings, HDF)

! and close the pattern file
call VT%closeExpPatternFile(HDF)

io_real(1) = maxval(expt) 
call Message%WriteValue('maximum intensity in pattern ',io_real,1)

! turn it into a 2D pattern
allocate(pattern(binx, biny), pcopy(binx, biny), pint(binx,biny), ppp(binx,biny), stat=ierr)
if (enl%numav.gt.0) then 
  do kk=1,biny
    pcopy(1:binx,kk) = sumexpt((kk-1)*binx+1:kk*binx)
  end do
else
  do kk=1,biny
    pcopy(1:binx,kk) = expt((kk-1)*binx+1:kk*binx)
  end do
end if

! do we need to extract this pattern from the file and store it as an image file ?
if (trim(enl%patternfile).ne.'undefined') then
! allocate a byte array for the final output TIFF image that will contain all individual images
  allocate(output_image(binx,biny))
  image_filename = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(enl%patternfile)

  ma = maxval(pcopy)
  mi = minval(pcopy)

  do i=1,binx
    do j=1,biny
     int8val = int(255.0*(pcopy(i,biny-j+1)-mi)/(ma-mi))
     output_image(i,j) = int8val
    end do
  end do

 ! set up the image_t structure
  im = image_t(output_image)
  if(im%empty()) call Message%printMessage("EMEBSDDIpreview","failed to convert array to image")

 ! create the file
  call im%write(trim(image_filename), iostat, iomsg) ! format automatically detected from extension
  if(0.ne.iostat) then
    call Message%printMessage("failed to write image to file : "//iomsg)
  else  
    call Message%printMessage('  Selected pattern written to '//trim(image_filename))
  end if 
  deallocate(output_image)
end if

! define the high-pass filter width array and the nregions array
numr = (enl%nregionsmax - enl%nregionsmin) / enl%nregionsstepsize + 1
allocate(nrvals(numr))
nrvals = enl%nregionsmin + (/ ((i-1)*enl%nregionsstepsize, i=1,numr) /)

! the array for the hi pass filter parameter is a non-linear progression
numw = enl%hipasswnsteps
allocate(hpvals(numw))
do ii=1,numw
    hpvals(ii) = 2.0**(float(ii-1-numw)) * 2.0 * enl%hipasswmax
end do

! allocate a byte array for the final output TIFF image that will contain all individual images
nx = numw * binx
ny = numr * biny
allocate(output_image(nx,ny))
image_filename = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(enl%tifffile)

! next we need to set up the high-pass filter fftw plans
allocate(hpmask(binx,biny),inp(binx,biny),outp(binx,biny),stat=istat)
if (istat .ne. 0) stop 'could not allocate hpmask, inp, outp arrays'
allocate(rrdata(binx,biny),ffdata(binx,biny),stat=istat)
if (istat .ne. 0) stop 'could not allocate rrdata, ffdata arrays'
call init_HiPassFilter(dble(hpvals(1)), (/enl%numsx, enl%numsy /), hpmask, inp, outp, HPplanf, HPplanb) 

! the outer loop goes over the hipass filter width and is displayed horizontally in the final image
do ii=1,numw
! Hi-Pass filter
    pattern = pcopy
    rrdata = dble(pattern)
    ffdata = applyHiPassFilter(rrdata, (/ binx, biny /), dble(hpvals(ii)), hpmask, inp, outp, HPplanf, HPplanb)
    pattern = sngl(ffdata)

    ma = maxval(pattern)
    mi = minval(pattern)

    pint = nint(((pattern - mi) / (ma-mi))*255.0)
    xoffset = (ii-1) * binx + 1
    do jj=1,numr
! adaptive histogram equalization
        if (nrvals(jj).eq.0) then
            ppp = pint
        else
            ppp = adhisteq(nrvals(jj),binx,biny,pint)
        end if 

! and store the pattern in the correct spot in the output_image array (flipped upside down !!!)
        yoffset =  (numr-jj) * biny + 1
        do i=1,binx
          do j=1,biny
           output_image(xoffset+i-1, yoffset+j-1) = ppp(i,biny-j+1)
          end do
        end do
    end do

! regenerate the complex inverted Gaussian mask with the next value of the mask width
    hpmask = cmplx(1.D0,0.D0)
    do i=1,binx/2 
      x = dble(i)**2
      do j=1,biny/2
        y = dble(j)**2
        v2 = hpvals(ii) * ( x+y )
        if (v2.lt.30.D0) then
          val = 1.D0-dexp(-v2)
          hpmask(i,j) = cmplx(val, 0.D0)
          hpmask(binx+1-i,j) = cmplx(val, 0.D0)
          hpmask(i,biny+1-j) = cmplx(val, 0.D0)
          hpmask(binx+1-i,biny+1-j) = cmplx(val, 0.D0)
        end if
      end do
    end do
end do

! set up the image_t structure
im2 = image_t(output_image)
if(im2%empty()) call Message%printMessage("DIpreview","failed to convert array to image")

! create the file
call im2%write(trim(image_filename), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage(" Failed to write image to file : "//iomsg)
else  
  call Message%printMessage('  Preprocessed pattern array written to '//trim(image_filename))
end if 
deallocate(output_image)

call Message%printMessage('')
call Message%printMessage(' High-pass filter parameter values along horizontal axis (L to R) :')
do ii=1,numw
    io_real(1) = hpvals(ii)
    call Message%WriteValue('',io_real,1,"(F10.6)")
end do

call Message%printMessage('')
call Message%printMessage(' nregions values along vertical axis (B to T):')
do ii=1,numr
    io_int(1) = nrvals(ii)
    call Message%WriteValue('',io_int,1)
end do

! close the fortran HDF interface
call closeFortranHDFInterface()

end associate

end subroutine DIpreview_



end module mod_DIpreview