! ###################################################################
! Copyright (c) 2013-2026, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_ppEBSD
  !! author: MDG
  !! version: 1.0
  !! date: 12/06/24
  !!
  !! class definition for the EMppEBSD program

use mod_kinds
use mod_global

IMPLICIT NONE

! namelist for the EMgetppEBSD program
type, public :: ppEBSDNameListType
 integer(kind=irg)  :: ipf_ht
 integer(kind=irg)  :: ipf_wd
 integer(kind=irg)  :: maskradius
 integer(kind=irg)  :: numsx
 integer(kind=irg)  :: numsy
 integer(kind=irg)  :: nthreads
 integer(kind=irg)  :: nregions
 integer(kind=irg)  :: logparam
 integer(kind=irg)  :: ROI(4)
 real(kind=dbl)     :: hipassw
 character(3)       :: filtertype
 character(1)       :: maskpattern
 character(1)       :: filterpattern
 character(fnlen)   :: exptfile
 character(fnlen)   :: tmpfile
 character(fnlen)   :: tiffname
 character(fnlen)   :: maskfile
 character(fnlen)   :: inputtype
 character(fnlen)   :: HDFstrings(10)
end type ppEBSDNameListType

! class definition
type, public :: ppEBSD_T
private
  character(fnlen)          :: nmldeffile = 'EMppEBSD.nml'
  type(ppEBSDNameListType)  :: nml

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: ppEBSD_
  procedure, pass(self) :: get_nml_
  procedure, pass(self) :: get_ipf_ht_
  procedure, pass(self) :: get_ipf_wd_
  procedure, pass(self) :: get_maskradius_
  procedure, pass(self) :: get_numsx_
  procedure, pass(self) :: get_numsy_
  procedure, pass(self) :: get_nthreads_
  procedure, pass(self) :: get_nregions_
  procedure, pass(self) :: get_logparam_
  procedure, pass(self) :: get_ROI_
  procedure, pass(self) :: get_hipassw_
  procedure, pass(self) :: get_filtertype_
  procedure, pass(self) :: get_maskpattern_
  procedure, pass(self) :: get_exptfile_
  procedure, pass(self) :: get_tmpfile_
  procedure, pass(self) :: get_maskfile_
  procedure, pass(self) :: get_inputtype_
  procedure, pass(self) :: get_HDFstrings_
  procedure, pass(self) :: set_ipf_ht_
  procedure, pass(self) :: set_ipf_wd_
  procedure, pass(self) :: set_maskradius_
  procedure, pass(self) :: set_numsx_
  procedure, pass(self) :: set_numsy_
  procedure, pass(self) :: set_nthreads_
  procedure, pass(self) :: set_nregions_
  procedure, pass(self) :: set_ROI_
  procedure, pass(self) :: set_hipassw_
  procedure, pass(self) :: set_filtertype_
  procedure, pass(self) :: set_maskpattern_
  procedure, pass(self) :: set_exptfile_
  procedure, pass(self) :: set_tmpfile_
  procedure, pass(self) :: set_maskfile_
  procedure, pass(self) :: set_inputtype_
  procedure, pass(self) :: set_HDFstrings_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: ppEBSD => ppEBSD_
  generic, public :: get_nml => get_nml_
  generic, public :: get_ipf_ht => get_ipf_ht_
  generic, public :: get_ipf_wd => get_ipf_wd_
  generic, public :: get_maskradius => get_maskradius_
  generic, public :: get_numsx => get_numsx_
  generic, public :: get_numsy => get_numsy_
  generic, public :: get_nthreads => get_nthreads_
  generic, public :: get_nregions => get_nregions_
  generic, public :: get_logparam => get_logparam_
  generic, public :: get_ROI => get_ROI_
  generic, public :: get_hipassw => get_hipassw_
  generic, public :: get_filtertype => get_filtertype_
  generic, public :: get_maskpattern => get_maskpattern_
  generic, public :: get_exptfile => get_exptfile_
  generic, public :: get_tmpfile => get_tmpfile_
  generic, public :: get_maskfile => get_maskfile_
  generic, public :: get_inputtype => get_inputtype_
  generic, public :: get_HDFstrings => get_HDFstrings_
  generic, public :: set_ipf_ht => set_ipf_ht_
  generic, public :: set_ipf_wd => set_ipf_wd_
  generic, public :: set_maskradius => set_maskradius_
  generic, public :: set_numsx => set_numsx_
  generic, public :: set_numsy => set_numsy_
  generic, public :: set_nthreads => set_nthreads_
  generic, public :: set_nregions => set_nregions_
  generic, public :: set_ROI => set_ROI_
  generic, public :: set_hipassw => set_hipassw_
  generic, public :: set_filtertype => set_filtertype_
  generic, public :: set_maskpattern => set_maskpattern_
  generic, public :: set_exptfile => set_exptfile_
  generic, public :: set_tmpfile => set_tmpfile_
  generic, public :: set_maskfile => set_maskfile_
  generic, public :: set_inputtype => set_inputtype_
  generic, public :: set_HDFstrings => set_HDFstrings_
end type ppEBSD_T

! the constructor routine for this class
interface ppEBSD_T
  module procedure ppEBSD_constructor
end interface ppEBSD_T

contains

!--------------------------------------------------------------------------
type(ppEBSD_T) function ppEBSD_constructor( nmlfile ) result(ppEBSD)
!DEC$ ATTRIBUTES DLLEXPORT :: ppEBSD_constructor
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! constructor for the ppEBSD_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call ppEBSD%readNameList(nmlfile)

end function ppEBSD_constructor

!--------------------------------------------------------------------------
subroutine ppEBSD_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: ppEBSD_destructor
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! destructor for the ppEBSD_T Class

IMPLICIT NONE

type(ppEBSD_T), INTENT(INOUT)  :: self

call reportDestructor('ppEBSD_T')

end subroutine ppEBSD_destructor

!--------------------------------------------------------------------------
function get_nml_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nml_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! get nml from the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT)    :: self
type(ppEBSDNameListType)          :: out

out = self%nml

end function get_nml_

!--------------------------------------------------------------------------
function get_ipf_ht_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_ht_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! get ipf_ht from the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
integer(kind=irg)              :: out

out = self%nml%ipf_ht

end function get_ipf_ht_

!--------------------------------------------------------------------------
subroutine set_ipf_ht_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_ht_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! set ipf_ht in the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)  :: inp

self%nml%ipf_ht = inp

end subroutine set_ipf_ht_

!--------------------------------------------------------------------------
function get_ipf_wd_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_wd_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! get ipf_wd from the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
integer(kind=irg)              :: out

out = self%nml%ipf_wd

end function get_ipf_wd_

!--------------------------------------------------------------------------
subroutine set_ipf_wd_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_wd_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! set ipf_wd in the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)  :: inp

self%nml%ipf_wd = inp

end subroutine set_ipf_wd_

!--------------------------------------------------------------------------
function get_maskradius_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_maskradius_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! get maskradius from the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
integer(kind=irg)              :: out

out = self%nml%maskradius

end function get_maskradius_

!--------------------------------------------------------------------------
subroutine set_maskradius_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_maskradius_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! set maskradius in the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)  :: inp

self%nml%maskradius = inp

end subroutine set_maskradius_

!--------------------------------------------------------------------------
function get_numsx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_numsx_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! get numsx from the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
integer(kind=irg)              :: out

out = self%nml%numsx

end function get_numsx_

!--------------------------------------------------------------------------
subroutine set_numsx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_numsx_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! set numsx in the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)  :: inp

self%nml%numsx = inp

end subroutine set_numsx_

!--------------------------------------------------------------------------
function get_numsy_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_numsy_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! get numsy from the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
integer(kind=irg)              :: out

out = self%nml%numsy

end function get_numsy_

!--------------------------------------------------------------------------
subroutine set_numsy_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_numsy_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! set numsy in the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)  :: inp

self%nml%numsy = inp

end subroutine set_numsy_

!--------------------------------------------------------------------------
function get_nthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nthreads_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! get nthreads from the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
integer(kind=irg)              :: out

out = self%nml%nthreads

end function get_nthreads_

!--------------------------------------------------------------------------
subroutine set_nthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nthreads_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! set nthreads in the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)  :: inp

self%nml%nthreads = inp

end subroutine set_nthreads_

!--------------------------------------------------------------------------
function get_nregions_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nregions_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! get nregions from the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
integer(kind=irg)              :: out

out = self%nml%nregions

end function get_nregions_

!--------------------------------------------------------------------------
subroutine set_nregions_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nregions_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! set nregions in the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)  :: inp

self%nml%nregions = inp

end subroutine set_nregions_

!--------------------------------------------------------------------------
function get_logparam_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_logparam_
!! author: MDG
!! version: 1.0
!! date: 01/23/25
!!
!! get logparam from the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
integer(kind=irg)              :: out

out = self%nml%logparam

end function get_logparam_

!--------------------------------------------------------------------------
subroutine set_logparam_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_logparam_
!! author: MDG
!! version: 1.0
!! date: 01/23/25
!!
!! set logparam in the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)  :: inp

self%nml%logparam = inp

end subroutine set_logparam_

!--------------------------------------------------------------------------
function get_ROI_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ROI_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! get ROI from the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
integer(kind=irg)              :: out(4)

out = self%nml%ROI

end function get_ROI_

!--------------------------------------------------------------------------
subroutine set_ROI_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ROI_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! set ROI in the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)  :: inp(4)

self%nml%ROI = inp

end subroutine set_ROI_

!--------------------------------------------------------------------------
function get_hipassw_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_hipassw_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! get hipassw from the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
real(kind=dbl)                 :: out

out = self%nml%hipassw

end function get_hipassw_

!--------------------------------------------------------------------------
subroutine set_hipassw_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_hipassw_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! set hipassw in the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
real(kind=dbl), INTENT(IN)     :: inp

self%nml%hipassw = inp

end subroutine set_hipassw_

!--------------------------------------------------------------------------
function get_maskpattern_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_maskpattern_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! get maskpattern from the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
character(1)                   :: out

out = self%nml%maskpattern

end function get_maskpattern_

!--------------------------------------------------------------------------
subroutine set_maskpattern_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_maskpattern_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! set maskpattern in the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
character(1), INTENT(IN)       :: inp

self%nml%maskpattern = inp

end subroutine set_maskpattern_

!--------------------------------------------------------------------------
function get_exptfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_exptfile_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! get exptfile from the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
character(fnlen)               :: out

out = self%nml%exptfile

end function get_exptfile_

!--------------------------------------------------------------------------
subroutine set_exptfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_exptfile_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! set exptfile in the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
character(fnlen), INTENT(IN)   :: inp

self%nml%exptfile = inp

end subroutine set_exptfile_

!--------------------------------------------------------------------------
function get_filtertype_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_filtertype_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! get filtertype from the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
character(fnlen)               :: out

out = self%nml%filtertype

end function get_filtertype_

!--------------------------------------------------------------------------
subroutine set_filtertype_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_filtertype_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! set filtertype in the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
character(fnlen), INTENT(IN)   :: inp

self%nml%filtertype = inp

end subroutine set_filtertype_

!--------------------------------------------------------------------------
function get_tmpfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_tmpfile_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! get tmpfile from the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
character(fnlen)               :: out

out = self%nml%tmpfile

end function get_tmpfile_

!--------------------------------------------------------------------------
subroutine set_tmpfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_tmpfile_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! set tmpfile in the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
character(fnlen), INTENT(IN)   :: inp

self%nml%tmpfile = inp

end subroutine set_tmpfile_

!--------------------------------------------------------------------------
function get_maskfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_maskfile_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! get maskfile from the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
character(fnlen)               :: out

out = self%nml%maskfile

end function get_maskfile_

!--------------------------------------------------------------------------
subroutine set_maskfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_maskfile_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! set maskfile in the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
character(fnlen), INTENT(IN)   :: inp

self%nml%maskfile = inp

end subroutine set_maskfile_

!--------------------------------------------------------------------------
function get_inputtype_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_inputtype_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! get inputtype from the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
character(fnlen)               :: out

out = self%nml%inputtype

end function get_inputtype_

!--------------------------------------------------------------------------
subroutine set_inputtype_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_inputtype_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! set inputtype in the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
character(fnlen), INTENT(IN)   :: inp

self%nml%inputtype = inp

end subroutine set_inputtype_

!--------------------------------------------------------------------------
function get_HDFstrings_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_HDFstrings_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! get HDFstrings from the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
character(fnlen)               :: out(10)

out = self%nml%HDFstrings

end function get_HDFstrings_

!--------------------------------------------------------------------------
subroutine set_HDFstrings_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_HDFstrings_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! set HDFstrings in the ppEBSD_T class

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
character(fnlen), INTENT(IN)   :: inp(10)

self%nml%HDFstrings = inp

end subroutine set_HDFstrings_

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! read the namelist from an nml file for the ppEBSD_T Class

use mod_io

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT) :: self
character(fnlen),INTENT(IN)    :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)    :: initonly
 !! fill in the default values only; do not read the file

type(IO_T)                     :: Message
logical                        :: skipread = .FALSE.

integer(kind=irg)              :: ipf_ht
integer(kind=irg)              :: ipf_wd
integer(kind=irg)              :: maskradius
integer(kind=irg)              :: numsx
integer(kind=irg)              :: numsy
integer(kind=irg)              :: nthreads
integer(kind=irg)              :: nregions
integer(kind=irg)              :: ROI(4)
integer(kind=irg)              :: logparam
real(kind=dbl)                 :: hipassw
character(1)                   :: maskpattern
character(3)                   :: filtertype
character(fnlen)               :: exptfile
character(fnlen)               :: tmpfile
character(fnlen)               :: maskfile
character(fnlen)               :: inputtype
character(fnlen)               :: HDFstrings(10)

! define the IO namelist to facilitate passing variables to the program.
namelist  / ppEBSDdata / numsx, numsy, nregions, maskpattern, nthreads, ipf_ht, ipf_wd, exptfile, maskradius, inputtype, &
                         filtertype, tmpfile, maskfile, HDFstrings, hipassw, ROI, logparam 

! set the input parameters to default values
 ipf_ht = 100
 ipf_wd = 100
 maskfile = 'undefined'
 maskpattern = 'n'
 maskradius = 240
 hipassw = 0.05
 nregions = 10
 logparam = 10
 numsx = 0
 numsy = 0
 ROI = (/ 0, 0, 0, 0 /)
 filtertype = 'fft'      ! or 'log'
 exptfile = 'undefined'
 inputtype = 'Binary'
 HDFstrings = ''
 tmpfile = 'EMEBSDDict_tmp.data'
 nthreads = 1

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=ppEBSDdata)
    close(UNIT=dataunit,STATUS='keep')

! check for required entries
    if (trim(exptfile).eq.'undefined') then
        call Message%printError('readNameList:',' experimental file name is undefined in '//nmlfile)
    end if

    if (numsx.eq.0) then
        call Message%printError('readNameList:',' patterns size numsx is zero in '//nmlfile)
    end if

    if (numsy.eq.0) then
        call Message%printError('readNameList:',' patterns size numsy is zero in '//nmlfile)
    end if
 end if

! if we get here, then all appears to be ok, and we need to fill in the enl fields
self%nml%ipf_ht = ipf_ht
self%nml%ipf_wd = ipf_wd
self%nml%maskradius = maskradius
self%nml%numsx = numsx
self%nml%numsy = numsy
self%nml%nthreads = nthreads
self%nml%nregions = nregions
self%nml%logparam = logparam
self%nml%ROI = ROI
self%nml%hipassw = hipassw
self%nml%filtertype = filtertype
self%nml%maskpattern = maskpattern
self%nml%exptfile = exptfile
self%nml%tmpfile = tmpfile
self%nml%maskfile = maskfile
self%nml%inputtype = inputtype
self%nml%HDFstrings = HDFstrings

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! pass the namelist for the ppEBSD_T Class to the calling program

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT)          :: self
type(ppEBSDNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine ppEBSD_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: ppEBSD_
!! author: MDG
!! version: 1.0
!! date: 12/06/24
!!
!! perform the computations

use mod_EMsoft
use mod_patterns
use mod_filters
use mod_io
use HDF5
use mod_HDFsupport
use mod_DIfiles
use mod_memory

use, intrinsic                  :: iso_fortran_env

IMPLICIT NONE

class(ppEBSD_T), INTENT(INOUT)  :: self
type(EMsoft_T), INTENT(INOUT)   :: EMsoft
character(fnlen), INTENT(INOUT) :: progname

type(IO_T)                      :: Message
type(HDF_T)                     :: HDF
type(memory_T)                  :: mem 

integer(kind=irg)               :: num,ierr,irec,istat
integer(kind=irg)               :: L,totnumexpt,imght,imgwd,nnk, recordsize, iii, hdferr,&
                                         recordsize_correct, patsz
real(kind=sgl),allocatable            :: mask(:,:),masklin(:)

integer(kind=irg)                     :: i,j,ii,jj,kk,ll,mm,pp,qq, io_int(3), iiistart, iiiend, jjend
real(kind=sgl)                        :: mi, ma, io_real(2), tstart, tmp, vlen, tstop
integer(kind=irg)                     :: binx,biny,TID,nthreads
integer(kind=irg)                     :: correctsize
logical                               :: f_exists, ROIselected
character(1000)                       :: charline
type(EBSDDINameListType)              :: dinl
character(fnlen)                      :: fname

call openFortranHDFInterface()

associate(ppnl=>self%nml)

! memory class 
mem = memory_T()

! copy various constants from the namelist
L = ppnl%numsx*ppnl%numsy
totnumexpt = ppnl%ipf_wd*ppnl%ipf_ht
imght = ppnl%numsx
imgwd = ppnl%numsy
recordsize = L*4
binx = ppnl%numsx
biny = ppnl%numsy

! make sure that correctsize is a multiple of 16; if not, make it so
! (probably not necessary in this program but kep in place for consistency)
if (mod(L,16) .ne. 0) then
    correctsize = 16*ceiling(float(L)/16.0)
else
    correctsize = L
end if

! determine the experimental and dictionary sizes in bytes
recordsize_correct = correctsize*4
patsz              = correctsize

ROIselected = .FALSE.
if (sum(ppnl%ROI).ne.0) ROIselected = .TRUE.

!=========================================
! ALLOCATION AND INITIALIZATION OF ARRAYS
!=========================================
call mem%alloc(mask, (/ binx,biny /), 'mask', 1.0)
call mem%alloc(masklin, (/ L /), 'masklin', 0.0)

!=====================================================
! define the circular mask if necessary and convert to 1D vector
!=====================================================

if (trim(ppnl%maskfile).ne.'undefined') then
! read the mask from file; the mask can be defined by a 2D array of 0 and 1 values
! that is stored in row form as strings, e.g.
!    0000001110000000
!    0000011111000000
! ... etc
!
    f_exists = .FALSE.
    fname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(ppnl%maskfile)
    inquire(file=trim(fname), exist=f_exists)
    if (f_exists.eqv..TRUE.) then
      mask = 0.0
      open(unit=dataunit,file=trim(fname),status='old',form='formatted')
      do jj=biny,1,-1
        read(dataunit,"(A)") charline
        do ii=1,binx
          if (charline(ii:ii).eq.'1') mask(ii,jj) = 1.0
        end do
      end do
      close(unit=dataunit,status='keep')
    else
      call Message%printError('ppEBSD','maskfile '//trim(fname)//' does not exist')
    end if
else
    if (ppnl%maskpattern.eq.'y') then
      do ii = 1,biny
          do jj = 1,binx
              if((ii-biny/2)**2 + (jj-binx/2)**2 .ge. ppnl%maskradius**2) then
                  mask(jj,ii) = 0.0
              end if
          end do
      end do
    end if
end if

! convert the mask to a linear (1D) array
do ii = 1,biny
    do jj = 1,binx
        masklin((ii-1)*binx+jj) = mask(jj,ii)
    end do
end do


!=====================================================
! Preprocess all the experimental patterns and store
! them in a temporary file as vectors
!=====================================================

! first, make sure that this file does not already exist
f_exists = .FALSE.
fname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(ppnl%tmpfile)
inquire(file=trim(fname), exist=f_exists)

! if the file exists, delete and recreate it
if (f_exists) then
    open(unit=dataunit,file=trim(fname),&
      status='unknown',form='unformatted',access='direct',recl=recordsize_correct,iostat=ierr)
    close(unit=dataunit,status='delete')
end if

! copy the relevant ppnl parameters into the dinl structure
dinl%nthreads = ppnl%nthreads
dinl%hipassw = ppnl%hipassw
dinl%ROI = ppnl%ROI
dinl%ipf_wd = ppnl%ipf_wd
dinl%ipf_ht = ppnl%ipf_ht
dinl%tmpfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(ppnl%tmpfile)
dinl%exptfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(ppnl%exptfile)
dinl%inputtype = ppnl%inputtype
dinl%HDFstrings = ppnl%HDFstrings
dinl%nregions = ppnl%nregions
dinl%DIModality = 'EBSD'

HDF = HDF_T() 

if (ppnl%filtertype.eq.'fft') then 
! standard hipass and adaptive histogram equalization pre-processing step
  call PreProcessPatterns(EMsoft, HDF, .FALSE., dinl, binx, biny, masklin, correctsize, totnumexpt)
else
  ! logarithmic high-pass filter only
  call PreProcessPatterns(EMsoft, HDF, .FALSE., dinl, binx, biny, masklin, correctsize, totnumexpt, &
                          log=.TRUE., logparam=ppnl%logparam)
end if 

!=====================================================
call Message%printMessage(' -> created pre-processed EBSD pattern file '//trim(fname))

! close the fortran HDF5 interface
call closeFortranHDFInterface()

end associate

call mem%dealloc(mask, 'mask')
call mem%dealloc(masklin, 'masklin')

! call mem%allocated_memory_use(' from end of ppEBSD_ subroutine ')

end subroutine ppEBSD_



end module mod_ppEBSD
