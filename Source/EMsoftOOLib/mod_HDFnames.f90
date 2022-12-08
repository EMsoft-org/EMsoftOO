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

module mod_HDFnames
 !! author: MDG
 !! version: 1.0
 !! date: 03/03/20
 !!
 !! this module contains a class that defines HDF group names that are frequently used; all
 !! HDF files generated by EMsoft have the same internal layout to make it easy to generate
 !! readers in other packages, e.g., Matlab or python

use mod_kinds
use mod_global

IMPLICIT NONE

private

type, public :: HDFnames_T
  private
    character(fnlen) :: EMheader
    character(fnlen) :: EMData
    character(fnlen) :: NMLfiles
    character(fnlen) :: NMLfilename
    character(fnlen) :: NMLparameters
    character(fnlen) :: NMLlist
    character(fnlen) :: ProgramData
    character(fnlen) :: CrystalData
    character(fnlen) :: Variable

  contains
  private

    procedure, pass(self) :: set_EMheader_
    procedure, pass(self) :: set_EMData_
    procedure, pass(self) :: set_NMLfiles_
    procedure, pass(self) :: set_NMLfilename_
    procedure, pass(self) :: set_NMLparameters_
    procedure, pass(self) :: set_NMLlist_
    procedure, pass(self) :: set_ProgramData_
    procedure, pass(self) :: set_CrystalData_
    procedure, pass(self) :: set_Variable_
    procedure, pass(self) :: get_EMheader_
    procedure, pass(self) :: get_EMData_
    procedure, pass(self) :: get_NMLfiles_
    procedure, pass(self) :: get_NMLfilename_
    procedure, pass(self) :: get_NMLparameters_
    procedure, pass(self) :: get_NMLlist_
    procedure, pass(self) :: get_ProgramData_
    procedure, pass(self) :: get_CrystalData_
    procedure, pass(self) :: get_Variable_
    procedure, pass(self) :: get_AllNames_
    final :: HDFnames_destructor

    generic, public :: set_EMheader => set_EMheader_
    generic, public :: set_EMData => set_EMData_
    generic, public :: set_NMLfiles => set_NMLfiles_
    generic, public :: set_NMLfilename => set_NMLfilename_
    generic, public :: set_NMLparameters => set_NMLparameters_
    generic, public :: set_NMLlist => set_NMLlist_
    generic, public :: set_ProgramData => set_ProgramData_
    generic, public :: set_CrystalData => set_CrystalData_
    generic, public :: set_Variable => set_Variable_
    generic, public :: get_EMheader => get_EMheader_
    generic, public :: get_EMData => get_EMData_
    generic, public :: get_NMLfiles => get_NMLfiles_
    generic, public :: get_NMLfilename => get_NMLfilename_
    generic, public :: get_NMLparameters => get_NMLparameters_
    generic, public :: get_NMLlist => get_NMLlist_
    generic, public :: get_ProgramData => get_ProgramData_
    generic, public :: get_CrystalData => get_CrystalData_
    generic, public :: get_Variable => get_Variable_
    generic, public :: get_AllNames => get_AllNames_

end type HDFnames_T

! the constructor routine for this class
interface HDFnames_T
  module procedure HDFnames_constructor
end interface HDFnames_T

contains

!--------------------------------------------------------------------------
type(HDFnames_T) function HDFnames_constructor( ) result(HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: HDFnames_constructor
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! constructor for the HDFnames_T Class; sets a couple default group names

use stringconstants

IMPLICIT NONE

HDFnames%EMheader = SC_EMheader
HDFnames%CrystalData = SC_CrystalData
HDFnames%EMData = SC_EMData
HDFnames%NMLfiles = SC_NMLfiles
HDFnames%NMLparameters = SC_NMLparameters

HDFnames%NMLfilename = ''
HDFnames%NMLlist = ''
HDFnames%ProgramData = ''
HDFnames%Variable = ''

end function HDFnames_constructor

!--------------------------------------------------------------------------
subroutine HDFnames_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: HDFnames_destructor
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! destructor for the HDFnames_T Class

IMPLICIT NONE

type(HDFnames_T), INTENT(INOUT)  :: self

call reportDestructor('HDFnames_T')

end subroutine HDFnames_destructor

!--------------------------------------------------------------------------
function get_EMheader_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_EMheader_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! get EMheader from the HDFnames_T class

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)     :: self
character(fnlen)                     :: out

out = self%EMheader

end function get_EMheader_

!--------------------------------------------------------------------------
subroutine set_EMheader_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_EMheader_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! set EMheader in the HDFnames_T class

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)     :: self
character(*), INTENT(IN)             :: inp

self%EMheader = inp

end subroutine set_EMheader_

!--------------------------------------------------------------------------
function get_NMLparameters_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_NMLparameters_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! get NMLparameters from the HDFnames_T class

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)     :: self
character(fnlen)                     :: out

out = self%NMLparameters

end function get_NMLparameters_

!--------------------------------------------------------------------------
subroutine set_NMLparameters_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_NMLparameters_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! set NMLparameters in the HDFnames_T class

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)     :: self
character(*), INTENT(IN)             :: inp

self%NMLparameters = inp

end subroutine set_NMLparameters_

!--------------------------------------------------------------------------
function get_EMData_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_EMData_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! get EMData from the HDFnames_T class

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)     :: self
character(fnlen)                     :: out

out = self%EMData

end function get_EMData_

!--------------------------------------------------------------------------
subroutine set_EMData_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_EMData_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! set EMData in the HDFnames_T class

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)     :: self
character(*), INTENT(IN)             :: inp

self%EMData = inp

end subroutine set_EMData_

!--------------------------------------------------------------------------
function get_NMLfiles_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_NMLfiles_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! get NMLfiles from the HDFnames_T class

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)     :: self
character(fnlen)                     :: out

out = self%NMLfiles

end function get_NMLfiles_

!--------------------------------------------------------------------------
subroutine set_NMLfiles_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_NMLfiles_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! set NMLfiles in the HDFnames_T class

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)     :: self
character(*), INTENT(IN)             :: inp

self%NMLfiles = inp

end subroutine set_NMLfiles_

!--------------------------------------------------------------------------
function get_NMLlist_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_NMLlist_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! get NMLlist from the HDFnames_T class

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)     :: self
character(fnlen)                     :: out

out = self%NMLlist

end function get_NMLlist_

!--------------------------------------------------------------------------
subroutine set_NMLlist_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_NMLlist_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! set NMLlist in the HDFnames_T class

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)     :: self
character(*), INTENT(IN)             :: inp

self%NMLlist = inp

end subroutine set_NMLlist_

!--------------------------------------------------------------------------
function get_ProgramData_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ProgramData_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! get ProgramData from the HDFnames_T class

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)     :: self
character(fnlen)                     :: out

out = self%ProgramData

end function get_ProgramData_

!--------------------------------------------------------------------------
subroutine set_ProgramData_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ProgramData_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! set ProgramData in the HDFnames_T class

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)     :: self
character(*), INTENT(IN)             :: inp

self%ProgramData = inp

end subroutine set_ProgramData_

!--------------------------------------------------------------------------
function get_CrystalData_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_CrystalData_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! get CrystalData from the HDFnames_T class

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)     :: self
character(fnlen)                     :: out

out = self%CrystalData

end function get_CrystalData_

!--------------------------------------------------------------------------
subroutine set_CrystalData_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_CrystalData_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! set CrystalData in the HDFnames_T class

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)     :: self
character(*), INTENT(IN)             :: inp

self%CrystalData = inp

end subroutine set_CrystalData_

!--------------------------------------------------------------------------
function get_NMLfilename_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_NMLfilename_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! get NMLfilename from the HDFnames_T class

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)     :: self
character(fnlen)                     :: out

out = self%NMLfilename

end function get_NMLfilename_

!--------------------------------------------------------------------------
subroutine set_NMLfilename_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_NMLfilename_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! set NMLfilename in the HDFnames_T class

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)     :: self
character(*), INTENT(IN)             :: inp

self%NMLfilename = inp

end subroutine set_NMLfilename_

!--------------------------------------------------------------------------
function get_Variable_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Variable_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! get Variable from the HDFnames_T class

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)     :: self
character(fnlen)                     :: out

out = self%Variable

end function get_Variable_

!--------------------------------------------------------------------------
subroutine set_Variable_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Variable_
!! author: MDG
!! version: 1.0
!! date: 03/03/20
!!
!! set Variable in the HDFnames_T class

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)     :: self
character(*), INTENT(IN)             :: inp

self%Variable = inp

end subroutine set_Variable_

!--------------------------------------------------------------------------
subroutine get_AllNames_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: get_AllNames_
!! author: MDG
!! version: 1.0
!! date: 03/17/20
!!
!! print all variables in the HDFnames_T class

use mod_io

IMPLICIT NONE

class(HDFnames_T), INTENT(INOUT)  :: self

type(IO_T)                        :: Message

call Message%printMessage(' Current components of the HDFnames class')
call Message%printMessage(' EMheader      : '//trim(self%get_EMheader()))
call Message%printMessage(' EMData        : '//trim(self%get_EMData()))
call Message%printMessage(' NMLfiles      : '//trim(self%get_NMLfiles()))
call Message%printMessage(' NMLfilename   : '//trim(self%get_NMLfilename()))
call Message%printMessage(' NMLparameters : '//trim(self%get_NMLparameters()))
call Message%printMessage(' NMLlist       : '//trim(self%get_NMLlist()))
call Message%printMessage(' ProgramData   : '//trim(self%get_ProgramData()))
call Message%printMessage(' CrystalData   : '//trim(self%get_CrystalData()))
call Message%printMessage(' Variable      : '//trim(self%get_Variable()))

end subroutine get_AllNames_




end module mod_HDFnames
