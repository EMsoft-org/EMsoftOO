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

module mod_ADP
  !! author: MDG
  !! version: 1.0
  !! date: 04/05/20
  !!
  !! class definition for the EMgetADP program

use mod_kinds
use mod_global

IMPLICIT NONE

! namelist for the EMgetADP program
type, public :: ADPNameListType
 integer(kind=irg)  :: ipf_ht
 integer(kind=irg)  :: ipf_wd
 integer(kind=irg)  :: maskradius
 integer(kind=irg)  :: numsx
 integer(kind=irg)  :: numsy
 integer(kind=irg)  :: nthreads
 integer(kind=irg)  :: nregions
 integer(kind=irg)  :: ROI(4)
 real(kind=dbl)     :: hipassw
 character(1)       :: maskpattern
 character(1)       :: filterpattern
 character(1)       :: keeptmpfile
 character(1)       :: usetmpfile
 character(fnlen)   :: exptfile
 character(fnlen)   :: tmpfile
 character(fnlen)   :: tiffname
 character(fnlen)   :: maskfile
 character(fnlen)   :: inputtype
 character(fnlen)   :: HDFstrings(10)
end type ADPNameListType

! class definition
type, public :: ADP_T
private
  character(fnlen)       :: nmldeffile = 'EMgetADP.nml'
  type(ADPNameListType)  :: nml

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: ADP_
  procedure, pass(self) :: get_ipf_ht_
  procedure, pass(self) :: get_ipf_wd_
  procedure, pass(self) :: get_maskradius_
  procedure, pass(self) :: get_numsx_
  procedure, pass(self) :: get_numsy_
  procedure, pass(self) :: get_nthreads_
  procedure, pass(self) :: get_nregions_
  procedure, pass(self) :: get_ROI_
  procedure, pass(self) :: get_hipassw_
  procedure, pass(self) :: get_maskpattern_
  procedure, pass(self) :: get_filterpattern_
  procedure, pass(self) :: get_keeptmpfile_
  procedure, pass(self) :: get_usetmpfile_
  procedure, pass(self) :: get_exptfile_
  procedure, pass(self) :: get_tmpfile_
  procedure, pass(self) :: get_tiffname_
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
  procedure, pass(self) :: set_maskpattern_
  procedure, pass(self) :: set_filterpattern_
  procedure, pass(self) :: set_keeptmpfile_
  procedure, pass(self) :: set_usetmpfile_
  procedure, pass(self) :: set_exptfile_
  procedure, pass(self) :: set_tmpfile_
  procedure, pass(self) :: set_tiffname_
  procedure, pass(self) :: set_maskfile_
  procedure, pass(self) :: set_inputtype_
  procedure, pass(self) :: set_HDFstrings_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: ADP => ADP_
  generic, public :: get_ipf_ht => get_ipf_ht_
  generic, public :: get_ipf_wd => get_ipf_wd_
  generic, public :: get_maskradius => get_maskradius_
  generic, public :: get_numsx => get_numsx_
  generic, public :: get_numsy => get_numsy_
  generic, public :: get_nthreads => get_nthreads_
  generic, public :: get_nregions => get_nregions_
  generic, public :: get_ROI => get_ROI_
  generic, public :: get_hipassw => get_hipassw_
  generic, public :: get_maskpattern => get_maskpattern_
  generic, public :: get_filterpattern => get_filterpattern_
  generic, public :: get_keeptmpfile => get_keeptmpfile_
  generic, public :: get_usetmpfile => get_usetmpfile_
  generic, public :: get_exptfile => get_exptfile_
  generic, public :: get_tmpfile => get_tmpfile_
  generic, public :: get_tiffname => get_tiffname_
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
  generic, public :: set_maskpattern => set_maskpattern_
  generic, public :: set_filterpattern => set_filterpattern_
  generic, public :: set_keeptmpfile => set_keeptmpfile_
  generic, public :: set_usetmpfile => set_usetmpfile_
  generic, public :: set_exptfile => set_exptfile_
  generic, public :: set_tmpfile => set_tmpfile_
  generic, public :: set_tiffname => set_tiffname_
  generic, public :: set_maskfile => set_maskfile_
  generic, public :: set_inputtype => set_inputtype_
  generic, public :: set_HDFstrings => set_HDFstrings_
end type ADP_T

! the constructor routine for this class
interface ADP_T
  module procedure ADP_constructor
end interface ADP_T

contains

!--------------------------------------------------------------------------
type(ADP_T) function ADP_constructor( nmlfile ) result(ADP)
!DEC$ ATTRIBUTES DLLEXPORT :: ADP_constructor
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! constructor for the ADP_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call ADP%readNameList(nmlfile)

end function ADP_constructor

!--------------------------------------------------------------------------
subroutine ADP_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: ADP_destructor
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! destructor for the ADP_T Class

IMPLICIT NONE

type(ADP_T), INTENT(INOUT)  :: self

call reportDestructor('ADP_T')

end subroutine ADP_destructor


!--------------------------------------------------------------------------
function get_ipf_ht_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_ht_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get ipf_ht from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out

out = self%nml%ipf_ht

end function get_ipf_ht_

!--------------------------------------------------------------------------
subroutine set_ipf_ht_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_ht_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set ipf_ht in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp

self%nml%ipf_ht = inp

end subroutine set_ipf_ht_

!--------------------------------------------------------------------------
function get_ipf_wd_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_wd_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get ipf_wd from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out

out = self%nml%ipf_wd

end function get_ipf_wd_

!--------------------------------------------------------------------------
subroutine set_ipf_wd_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_wd_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set ipf_wd in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp

self%nml%ipf_wd = inp

end subroutine set_ipf_wd_

!--------------------------------------------------------------------------
function get_maskradius_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_maskradius_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get maskradius from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out

out = self%nml%maskradius

end function get_maskradius_

!--------------------------------------------------------------------------
subroutine set_maskradius_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_maskradius_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set maskradius in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp

self%nml%maskradius = inp

end subroutine set_maskradius_

!--------------------------------------------------------------------------
function get_numsx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_numsx_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get numsx from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out

out = self%nml%numsx

end function get_numsx_

!--------------------------------------------------------------------------
subroutine set_numsx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_numsx_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set numsx in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp

self%nml%numsx = inp

end subroutine set_numsx_

!--------------------------------------------------------------------------
function get_numsy_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_numsy_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get numsy from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out

out = self%nml%numsy

end function get_numsy_

!--------------------------------------------------------------------------
subroutine set_numsy_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_numsy_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set numsy in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp

self%nml%numsy = inp

end subroutine set_numsy_

!--------------------------------------------------------------------------
function get_nthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nthreads_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get nthreads from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out

out = self%nml%nthreads

end function get_nthreads_

!--------------------------------------------------------------------------
subroutine set_nthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nthreads_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set nthreads in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp

self%nml%nthreads = inp

end subroutine set_nthreads_

!--------------------------------------------------------------------------
function get_nregions_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nregions_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get nregions from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out

out = self%nml%nregions

end function get_nregions_

!--------------------------------------------------------------------------
subroutine set_nregions_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nregions_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set nregions in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp

self%nml%nregions = inp

end subroutine set_nregions_

!--------------------------------------------------------------------------
function get_ROI_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ROI_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get ROI from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out(4)

out = self%nml%ROI

end function get_ROI_

!--------------------------------------------------------------------------
subroutine set_ROI_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ROI_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set ROI in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp(4)

self%nml%ROI = inp

end subroutine set_ROI_

!--------------------------------------------------------------------------
function get_hipassw_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_hipassw_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get hipassw from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
real(kind=dbl)                  :: out

out = self%nml%hipassw

end function get_hipassw_

!--------------------------------------------------------------------------
subroutine set_hipassw_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_hipassw_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set hipassw in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)      :: inp

self%nml%hipassw = inp

end subroutine set_hipassw_

!--------------------------------------------------------------------------
function get_maskpattern_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_maskpattern_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get maskpattern from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(1)                    :: out

out = self%nml%maskpattern

end function get_maskpattern_

!--------------------------------------------------------------------------
subroutine set_maskpattern_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_maskpattern_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set maskpattern in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(1), INTENT(IN)        :: inp

self%nml%maskpattern = inp

end subroutine set_maskpattern_

!--------------------------------------------------------------------------
function get_filterpattern_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_filterpattern_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get filterpattern from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(1)                    :: out

out = self%nml%filterpattern

end function get_filterpattern_

!--------------------------------------------------------------------------
subroutine set_filterpattern_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_filterpattern_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set filterpattern in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(1), INTENT(IN)        :: inp

self%nml%filterpattern = inp

end subroutine set_filterpattern_

!--------------------------------------------------------------------------
function get_keeptmpfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_keeptmpfile_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get keeptmpfile from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(1)                    :: out

out = self%nml%keeptmpfile

end function get_keeptmpfile_

!--------------------------------------------------------------------------
subroutine set_keeptmpfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_keeptmpfile_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set keeptmpfile in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(1), INTENT(IN)        :: inp

self%nml%keeptmpfile = inp

end subroutine set_keeptmpfile_

!--------------------------------------------------------------------------
function get_usetmpfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_usetmpfile_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get usetmpfile from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(1)                    :: out

out = self%nml%usetmpfile

end function get_usetmpfile_

!--------------------------------------------------------------------------
subroutine set_usetmpfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_usetmpfile_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set usetmpfile in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(1), INTENT(IN)        :: inp

self%nml%usetmpfile = inp

end subroutine set_usetmpfile_

!--------------------------------------------------------------------------
function get_exptfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_exptfile_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get exptfile from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(fnlen)                :: out

out = self%nml%exptfile

end function get_exptfile_

!--------------------------------------------------------------------------
subroutine set_exptfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_exptfile_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set exptfile in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)    :: inp

self%nml%exptfile = inp

end subroutine set_exptfile_

!--------------------------------------------------------------------------
function get_tmpfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_tmpfile_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get tmpfile from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(fnlen)                :: out

out = self%nml%tmpfile

end function get_tmpfile_

!--------------------------------------------------------------------------
subroutine set_tmpfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_tmpfile_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set tmpfile in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)    :: inp

self%nml%tmpfile = inp

end subroutine set_tmpfile_

!--------------------------------------------------------------------------
function get_tiffname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_tiffname_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get tiffname from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(fnlen)                :: out

out = self%nml%tiffname

end function get_tiffname_

!--------------------------------------------------------------------------
subroutine set_tiffname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_tiffname_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set tiffname in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)    :: inp

self%nml%tiffname = inp

end subroutine set_tiffname_

!--------------------------------------------------------------------------
function get_maskfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_maskfile_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get maskfile from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(fnlen)                :: out

out = self%nml%maskfile

end function get_maskfile_

!--------------------------------------------------------------------------
subroutine set_maskfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_maskfile_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set maskfile in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)    :: inp

self%nml%maskfile = inp

end subroutine set_maskfile_

!--------------------------------------------------------------------------
function get_inputtype_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_inputtype_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get inputtype from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(fnlen)                :: out

out = self%nml%inputtype

end function get_inputtype_

!--------------------------------------------------------------------------
subroutine set_inputtype_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_inputtype_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set inputtype in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)    :: inp

self%nml%inputtype = inp

end subroutine set_inputtype_

!--------------------------------------------------------------------------
function get_HDFstrings_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_HDFstrings_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! get HDFstrings from the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(fnlen)                :: out(10)

out = self%nml%HDFstrings

end function get_HDFstrings_

!--------------------------------------------------------------------------
subroutine set_HDFstrings_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_HDFstrings_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! set HDFstrings in the ADP_T class

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)    :: inp(10)

self%nml%HDFstrings = inp

end subroutine set_HDFstrings_

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! read the namelist from an nml file for the ADP_T Class

use mod_io

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(IO_T)                           :: Message
logical                              :: skipread = .FALSE.

integer(kind=irg)       :: ipf_ht
integer(kind=irg)       :: ipf_wd
integer(kind=irg)       :: maskradius
integer(kind=irg)       :: numsx
integer(kind=irg)       :: numsy
integer(kind=irg)       :: nthreads
integer(kind=irg)       :: nregions
integer(kind=irg)       :: ROI(4)
real(kind=dbl)          :: hipassw
character(1)            :: maskpattern
character(1)            :: filterpattern
character(1)            :: keeptmpfile
character(1)            :: usetmpfile
character(fnlen)        :: exptfile
character(fnlen)        :: tmpfile
character(fnlen)        :: tiffname
character(fnlen)        :: maskfile
character(fnlen)        :: inputtype
character(fnlen)        :: HDFstrings(10)

! define the IO namelist to facilitate passing variables to the program.
namelist  / getADP / numsx, numsy, nregions, maskpattern, nthreads, ipf_ht, ipf_wd, exptfile, maskradius, inputtype, &
                     tmpfile, maskfile, HDFstrings, hipassw, tiffname, filterpattern, keeptmpfile, usetmpfile, ROI

! set the input parameters to default values
 ipf_ht = 100
 ipf_wd = 100
 maskfile = 'undefined'
 filterpattern = 'y'
 maskpattern = 'n'
 keeptmpfile = 'n'
 usetmpfile = 'n'
 maskradius = 240
 hipassw = 0.05
 nregions = 10
 numsx = 0
 numsy = 0
 ROI = (/ 0, 0, 0, 0 /)
 exptfile = 'undefined'
 inputtype = 'Binary'
 HDFstrings = ''
 tmpfile = 'EMEBSDDict_tmp.data'
 tiffname = 'undefined'
 nthreads = 1

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=getADP)
    close(UNIT=dataunit,STATUS='keep')

! check for required entries
    if (trim(exptfile).eq.'undefined') then
        call Message%printError('readNameList:',' experimental file name is undefined in '//nmlfile)
    end if

    if (trim(tiffname).eq.'undefined') then
        call Message%printError('readNameList:',' output tiff file name is undefined in '//nmlfile)
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
self%nml%ROI = ROI
self%nml%hipassw = hipassw
self%nml%maskpattern = maskpattern
self%nml%filterpattern = filterpattern
self%nml%keeptmpfile = keeptmpfile
self%nml%usetmpfile = usetmpfile
self%nml%exptfile = exptfile
self%nml%tmpfile = tmpfile
self%nml%tiffname = tiffname
self%nml%maskfile = maskfile
self%nml%inputtype = inputtype
self%nml%HDFstrings = HDFstrings

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! pass the namelist for the ADP_T Class to the calling program

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)          :: self
type(ADPNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine ADP_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: ADP_
!! author: MDG
!! version: 1.0
!! date: 04/05/20
!!
!! perform the computations

use mod_EMsoft
use mod_patterns
use mod_filters
use mod_io
use mod_math
use mod_timing
use omp_lib
use HDF5
use h5im
use h5lt
use mod_HDFsupport
use ISO_C_BINDING
use mod_DIfiles
use mod_image

use, intrinsic :: iso_fortran_env

IMPLICIT NONE

class(ADP_T), INTENT(INOUT)                         :: self
type(EMsoft_T), INTENT(INOUT)                       :: EMsoft
character(fnlen), INTENT(INOUT)                     :: progname

type(IO_T)                                          :: Message
type(HDF_T)                                         :: HDF
type(timing_T)                                      :: timer

integer(kind=irg)                                   :: num,ierr,irec,istat
integer(kind=irg),parameter                         :: iunit = 40
integer(kind=irg),parameter                         :: iunitexpt = 41
integer(kind=irg),parameter                         :: iunitdict = 42
integer(kind=irg),parameter                         :: itmpexpt = 43

integer(kind=irg)                                   :: L,totnumexpt,imght,imgwd,nnk, recordsize, iii, hdferr,&
                                                       recordsize_correct, patsz, TIFF_nx, TIFF_ny
real(kind=sgl),allocatable                          :: imageexpt(:), mask(:,:),masklin(:), exppatarray(:), tmpexppatarray(:)
real(kind=sgl),allocatable                          :: imageexptflt(:),binned(:,:),imagedictflt(:),imagedictfltflip(:), &
                                                       tmpimageexpt(:)
real(kind=sgl),allocatable                          :: EBSDpattern(:,:), dpmap(:)
real(kind=sgl),allocatable                          :: EBSDpatternintd(:,:), EBSDpat(:,:)
integer(kind=irg),allocatable                       :: EBSDpatterninteger(:,:), EBSDpatternad(:,:), EBSDpint(:,:)
real(kind=dbl),allocatable                          :: rdata(:,:), fdata(:,:), rrdata(:,:), ffdata(:,:)
complex(kind=dbl),allocatable                       :: hpmask(:,:)
complex(C_DOUBLE_COMPLEX),allocatable               :: inp(:,:), outp(:,:)
character(11)                                       :: dstr
character(15)                                       :: tstrb
character(15)                                       :: tstre
real(kind=dbl)                                      :: w, Jres
integer(kind=irg)                                   :: dims(2)
character(fnlen, KIND=c_char),allocatable,TARGET    :: stringarray(:)
character(fnlen)                                    :: groupname, dataset, fname, ename, sourcefile, TIFF_filename
integer(hsize_t)                                    :: expwidth, expheight
integer(c_size_t),target                            :: slength
integer(c_int)                                      :: numd, nump
type(C_PTR)                                         :: planf, HPplanf, HPplanb
integer(HSIZE_T)                                    :: dims3(3), offset3(3)

integer(kind=irg)                                   :: i,j,ii,jj,kk,ll,mm,pp,qq, io_int(3), iiistart, iiiend, jjend
real(kind=sgl)                                      :: mi, ma, io_real(2), tstart, tmp, vlen, tstop
integer(kind=irg)                                   :: binx,biny,TID,nthreads
integer(kind=irg)                                   :: correctsize
logical                                             :: f_exists, ROIselected
character(1000)                                     :: charline
! type(DictionaryIndexingNameListType)                :: dinl
type(EBSDDINameListType)                            :: dinl

! declare variables for use in object oriented image module
integer                                             :: iostat
character(len=128)                                  :: iomsg
logical                                             :: isInteger
type(image_t)                                       :: im
integer(int8)                                       :: i8 (3,4)
integer(int8), allocatable                          :: TIFF_image(:,:)

call openFortranHDFInterface()

associate(adpnl=>self%nml)

! initialize the timing routines
timer = Timing_T()
tstrb = timer%getTimeString()

! copy various constants from the namelist
L = adpnl%numsx*adpnl%numsy
totnumexpt = adpnl%ipf_wd*adpnl%ipf_ht
imght = adpnl%numsx
imgwd = adpnl%numsy
recordsize = L*4
dims = (/imght, imgwd/)
w = adpnl%hipassw
binx = adpnl%numsx
biny = adpnl%numsy

! make sure that correctsize is a multiple of 16; if not, make it so
if (mod(L,16) .ne. 0) then
    correctsize = 16*ceiling(float(L)/16.0)
else
    correctsize = L
end if

! determine the experimental and dictionary sizes in bytes
recordsize_correct = correctsize*4
patsz              = correctsize

ROIselected = .FALSE.
if (sum(adpnl%ROI).ne.0) ROIselected = .TRUE.

!=========================================
! ALLOCATION AND INITIALIZATION OF ARRAYS
!=========================================
allocate(mask(binx,biny),masklin(L),stat=istat)
if (istat .ne. 0) stop 'Could not allocate arrays for masks'
mask = 1.0
masklin = 0.0

allocate(imageexpt(L),imageexptflt(correctsize),imagedictflt(correctsize),imagedictfltflip(correctsize),stat=istat)
allocate(tmpimageexpt(correctsize),stat=istat)
if (istat .ne. 0) stop 'Could not allocate array for reading experimental image patterns'
imageexpt = 0.0
imageexptflt = 0.0

allocate(EBSDpattern(binx,biny),binned(binx,biny),stat=istat)
if (istat .ne. 0) stop 'Could not allocate array for EBSD pattern'
EBSDpattern = 0.0
binned = 0.0

allocate(rdata(binx,biny),fdata(binx,biny),stat=istat)
if (istat .ne. 0) stop 'could not allocate arrays for Hi-Pass filter'
rdata = 0.D0
fdata = 0.D0

!=====================================================
! define the circular mask if necessary and convert to 1D vector
!=====================================================

if (trim(adpnl%maskfile).ne.'undefined') then
! read the mask from file; the mask can be defined by a 2D array of 0 and 1 values
! that is stored in row form as strings, e.g.
!    0000001110000000
!    0000011111000000
! ... etc
!
    f_exists = .FALSE.
    fname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(adpnl%maskfile)
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
      call Message%printError('ADP','maskfile '//trim(fname)//' does not exist')
    end if
else
    if (adpnl%maskpattern.eq.'y') then
      do ii = 1,biny
          do jj = 1,binx
              if((ii-biny/2)**2 + (jj-binx/2)**2 .ge. adpnl%maskradius**2) then
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
! them in a temporary file as vectors; also, create
! an average dot product map to be stored in the h5ebsd output file
!
! this could become a separate routine in the EMEBSDmod module ...
!=====================================================

! first, make sure that this file does not already exist
f_exists = .FALSE.
fname = trim(EMsoft%generateFilePath('EMtmppathname'))//trim(adpnl%tmpfile)
inquire(file=trim(fname), exist=f_exists)

call Message%WriteValue('Checking for temporary file ',trim(fname))

if ((adpnl%usetmpfile.eq.'y').and.(.not.f_exists)) then
  call Message%printError('ADP','tmp file does not exist ...')
end if

! if the file exists, and we do not want to keep it, then delete and recreate it
if (f_exists) then
  if (adpnl%usetmpfile.eq.'n') then   ! delete the file and open a new one
    open(unit=itmpexpt,file=trim(fname),&
      status='unknown',form='unformatted',access='direct',recl=recordsize_correct,iostat=ierr)
    close(unit=itmpexpt,status='delete')
  end if
end if

if (adpnl%usetmpfile.eq.'n') then

! copy the relevant adpnl parameters into the dinl structure
  dinl%nthreads = adpnl%nthreads
  dinl%hipassw = adpnl%hipassw
  dinl%ROI = adpnl%ROI
  dinl%ipf_wd = adpnl%ipf_wd
  dinl%ipf_ht = adpnl%ipf_ht
  dinl%tmpfile = adpnl%tmpfile
  dinl%exptfile = adpnl%exptfile
  dinl%inputtype = adpnl%inputtype
  dinl%HDFstrings = adpnl%HDFstrings
  dinl%nregions = adpnl%nregions

  call PreProcessPatterns(EMsoft, HDF, .FALSE., dinl, binx, biny, masklin, correctsize, totnumexpt)
end if

!=====================================================
call Message%printMessage(' -> computing Average Dot Product map (ADP)')
call Message%printMessage(' ')

! re-open the temporary file
open(unit=itmpexpt,file=trim(fname),&
     status='old',form='unformatted',access='direct',recl=recordsize_correct,iostat=ierr)

! use the getADPmap routine in the filters module
if (ROIselected.eqv..TRUE.) then
  allocate(dpmap(adpnl%ROI(3)*adpnl%ROI(4)))
  call getADPmap(itmpexpt, adpnl%ROI(3)*adpnl%ROI(4), L, adpnl%ROI(3), adpnl%ROI(4), dpmap)
  TIFF_nx = adpnl%ROI(3)
  TIFF_ny = adpnl%ROI(4)
else
  allocate(dpmap(totnumexpt))
  call getADPmap(itmpexpt, totnumexpt, L, adpnl%ipf_wd, adpnl%ipf_ht, dpmap)
  TIFF_nx = adpnl%ipf_wd
  TIFF_ny = adpnl%ipf_ht
end if

! output the ADP map as a tiff file
fname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(adpnl%tiffname)//'_ADP.tiff'
TIFF_filename = trim(fname)

! allocate memory for image
allocate(TIFF_image(TIFF_nx,TIFF_ny))

! fill the image with whatever data you have (between 0 and 255)
ma = maxval(dpmap)
mi = minval(dpmap)

do j=1,TIFF_ny
 do i=1,TIFF_nx
  ii = (j-1) * TIFF_nx + i
  TIFF_image(i,j) = int(255 * (dpmap(ii)-mi)/(ma-mi))
 end do
end do

! set up the image_t structure
im = image_t(TIFF_image)
if(im%empty()) call Message%printMessage("EMgetADP","failed to convert array to image")

! create the file
call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage("failed to write image to file : "//iomsg)
else
  call Message%printMessage('ADP map written to '//trim(TIFF_filename))
end if
deallocate(TIFF_image)

if (adpnl%keeptmpfile.eq.'n') then
  close(unit=itmpexpt, status = 'delete')
  call Message%printMessage(' -> tmp file deleted')
else
  call Message%printMessage(' -> keeping tmp file')
end if

! close the fortran HDF5 interface
call closeFortranHDFInterface()

end associate

end subroutine ADP_



end module mod_ADP
