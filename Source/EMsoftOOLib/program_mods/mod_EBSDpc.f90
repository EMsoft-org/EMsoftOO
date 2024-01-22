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

module mod_EBSDpc
  !! author: MDG 
  !! version: 1.0 
  !! date: 09/09/22
  !!
  !! class definition for the EMEBSDpc program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMEBSDpc program
type, public :: EBSDpcNameListType
  real(kind=sgl)    :: xstar 
  real(kind=sgl)    :: ystar 
  real(kind=sgl)    :: zstar 
  real(kind=sgl)    :: delta
  integer(kind=irg) :: Nx
  integer(kind=irg) :: Ny
  integer(kind=irg) :: binning
  character(fnlen)  :: convention
end type EBSDpcNameListType

! class definition
type, public :: EBSDpc_T
private 
  character(fnlen)          :: nmldeffile = 'EMEBSDpc.nml'
  type(EBSDpcNameListType)  :: nml 
  real(kind=sgl)            :: EMsoft_pc(3)
  real(kind=sgl)            :: TSL_pc(3)
  real(kind=sgl)            :: Oxford_pc(3)
  real(kind=sgl)            :: Bruker_pc(3)

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: EBSDpc_
  procedure, pass(self) :: get_xstar_
  procedure, pass(self) :: get_ystar_
  procedure, pass(self) :: get_zstar_
  procedure, pass(self) :: get_delta_
  procedure, pass(self) :: get_Nx_
  procedure, pass(self) :: get_Ny_
  procedure, pass(self) :: get_binning_
  procedure, pass(self) :: get_convention_
  procedure, pass(self) :: set_xstar_
  procedure, pass(self) :: set_ystar_
  procedure, pass(self) :: set_zstar_
  procedure, pass(self) :: set_delta_
  procedure, pass(self) :: set_Nx_
  procedure, pass(self) :: set_Ny_
  procedure, pass(self) :: set_binning_
  procedure, pass(self) :: set_convention_
  procedure, pass(self) :: EMsoft2TSL_
  procedure, pass(self) :: EMsoft2Oxford_
  procedure, pass(self) :: EMsoft2Bruker_
  procedure, pass(self) :: TSL2EMsoft_
  procedure, pass(self) :: Oxford2EMsoft_
  procedure, pass(self) :: Bruker2EMsoft_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: EBSDpc => EBSDpc_
  generic, public :: get_xstar => get_xstar_
  generic, public :: get_ystar => get_ystar_
  generic, public :: get_zstar => get_zstar_
  generic, public :: get_delta => get_delta_
  generic, public :: get_Nx => get_Nx_
  generic, public :: get_Ny => get_Ny_
  generic, public :: get_binning => get_binning_
  generic, public :: get_convention => get_convention_
  generic, public :: set_xstar => set_xstar_
  generic, public :: set_ystar => set_ystar_
  generic, public :: set_zstar => set_zstar_
  generic, public :: set_delta => set_delta_
  generic, public :: set_Nx => set_Nx_
  generic, public :: set_Ny => set_Ny_
  generic, public :: set_binning => set_binning_
  generic, public :: set_convention => set_convention_
end type EBSDpc_T

! the constructor routine for this class 
interface EBSDpc_T
  module procedure EBSDpc_constructor
end interface EBSDpc_T

contains

!--------------------------------------------------------------------------
type(EBSDpc_T) function EBSDpc_constructor( nmlfile ) result(EBSDpc)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDpc_constructor
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! constructor for the EBSDpc_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call EBSDpc%readNameList(nmlfile)

end function EBSDpc_constructor

!--------------------------------------------------------------------------
subroutine EBSDpc_destructor(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDpc_destructor
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! destructor for the EBSDpc_T Class
 
IMPLICIT NONE

type(EBSDpc_T), INTENT(INOUT)  :: self 

call reportDestructor('EBSDpc_T')

end subroutine EBSDpc_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! read the namelist from an nml file for the EBSDpc_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)       :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

real(kind=sgl)    :: xstar 
real(kind=sgl)    :: ystar 
real(kind=sgl)    :: zstar 
real(kind=sgl)    :: delta
integer(kind=irg) :: Nx
integer(kind=irg) :: Ny
integer(kind=irg) :: binning
character(fnlen)  :: convention

namelist /EBSDpc/ xstar, ystar, zstar, delta, Nx, Ny, binning, convention

xstar  = 0.0
ystar  = 0.0
zstar  = 0.0
delta = 50.0
Nx = 640
Ny = 480
binning = 1
convention = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=EBSDpc)
    close(UNIT=dataunit,STATUS='keep')

    if (trim(convention).eq.'undefined') then
      call Message%printError('readNameList:',' convention is undefined in '//nmlfile)
    end if
end if

self%nml%xstar = xstar
self%nml%ystar = ystar
self%nml%zstar = zstar
self%nml%delta = delta
self%nml%Nx = Nx
self%nml%Ny = Ny
self%nml%binning = binning
self%nml%convention = convention

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! pass the namelist for the EBSDpc_T Class to the calling program

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)          :: self
type(EBSDpcNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
function get_xstar_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_xstar_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! get xstar from the EBSDpc_T class

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self
real(kind=sgl)                     :: out

out = self%nml%xstar

end function get_xstar_

!--------------------------------------------------------------------------
subroutine set_xstar_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_xstar_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! set xstar in the EBSDpc_T class

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)         :: inp

self%nml%xstar = inp

end subroutine set_xstar_

!--------------------------------------------------------------------------
function get_ystar_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ystar_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! get ystar from the EBSDpc_T class

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self
real(kind=sgl)                     :: out

out = self%nml%ystar

end function get_ystar_

!--------------------------------------------------------------------------
subroutine set_ystar_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ystar_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! set ystar in the EBSDpc_T class

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)         :: inp

self%nml%ystar = inp

end subroutine set_ystar_

!--------------------------------------------------------------------------
function get_zstar_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_zstar_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! get zstar from the EBSDpc_T class

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self
real(kind=sgl)                     :: out

out = self%nml%zstar

end function get_zstar_

!--------------------------------------------------------------------------
subroutine set_zstar_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_zstar_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! set zstar in the EBSDpc_T class

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)         :: inp

self%nml%zstar = inp

end subroutine set_zstar_

!--------------------------------------------------------------------------
function get_delta_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_delta_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! get delta from the EBSDpc_T class

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self
real(kind=sgl)                     :: out

out = self%nml%delta

end function get_delta_

!--------------------------------------------------------------------------
subroutine set_delta_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_delta_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! set delta in the EBSDpc_T class

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)         :: inp

self%nml%delta = inp

end subroutine set_delta_

!--------------------------------------------------------------------------
function get_Nx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Nx_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! get Nx from the EBSDpc_T class

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self
integer(kind=irg)                  :: out

out = self%nml%Nx

end function get_Nx_

!--------------------------------------------------------------------------
subroutine set_Nx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Nx_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! set Nx in the EBSDpc_T class

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)      :: inp

self%nml%Nx = inp

end subroutine set_Nx_

!--------------------------------------------------------------------------
function get_Ny_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Ny_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! get Ny from the EBSDpc_T class

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self
integer(kind=irg)                  :: out

out = self%nml%Ny

end function get_Ny_

!--------------------------------------------------------------------------
subroutine set_Ny_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Ny_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! set Ny in the EBSDpc_T class

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)      :: inp

self%nml%Ny = inp

end subroutine set_Ny_

!--------------------------------------------------------------------------
function get_binning_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_binning_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! get binning from the EBSDpc_T class

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self
integer(kind=irg)                  :: out

out = self%nml%binning

end function get_binning_

!--------------------------------------------------------------------------
subroutine set_binning_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_binning_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! set binning in the EBSDpc_T class

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)      :: inp

self%nml%binning = inp

end subroutine set_binning_

!--------------------------------------------------------------------------
function get_convention_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_convention_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! get convention from the EBSDpc_T class

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = self%nml%convention

end function get_convention_

!--------------------------------------------------------------------------
subroutine set_convention_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_convention_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! set convention in the EBSDpc_T class

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%convention = inp

end subroutine set_convention_

!--------------------------------------------------------------------------
!-routines to convert to and from EMsoft convention and the others...     -
!--------------------------------------------------------------------------
subroutine EMsoft2TSL_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: EMsoft2TSL_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! convert pattern center from EMsoft to TSL convention 

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self

self%TSL_pc(1) = self%EMsoft_pc(1)/real(self%nml%Nx) + 0.5 
self%TSL_pc(2) = (self%EMsoft_pc(2)+real(self%nml%Ny)/2.0)/real(self%nml%Nx) 
self%TSL_pc(3) = self%EMsoft_pc(3)/real(self%nml%Nx)/self%nml%delta 

end subroutine EMsoft2TSL_

!--------------------------------------------------------------------------
subroutine TSL2EMsoft_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: TSL2EMsoft_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! convert pattern center from TSL to EMsoft convention 

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self

self%EMsoft_pc(1) = (self%TSL_pc(1)-0.5)*real(self%nml%Nx)
self%EMsoft_pc(2) = real(self%nml%Nx)*self%TSL_pc(2)-real(self%nml%Ny)/2.0
self%EMsoft_pc(3) = self%TSL_pc(3)*real(self%nml%Nx)*self%nml%delta 

end subroutine TSL2EMsoft_

!--------------------------------------------------------------------------
subroutine EMsoft2Oxford_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: EMsoft2Oxford_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! convert pattern center from EMsoft to Oxford convention 

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self

self%Oxford_pc(1) = self%EMsoft_pc(1)/real(self%nml%Nx) + 0.5 
self%Oxford_pc(2) = self%EMsoft_pc(2)/real(self%nml%Ny) + 0.5
self%Oxford_pc(3) = self%EMsoft_pc(3)/real(self%nml%Nx)/self%nml%delta 

end subroutine EMsoft2Oxford_

!--------------------------------------------------------------------------
subroutine Oxford2EMsoft_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: Oxford2EMsoft_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! convert pattern center from Oxford to EMsoft convention 

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self

self%EMsoft_pc(1) = (self%Oxford_pc(1)-0.5)*real(self%nml%Nx)
self%EMsoft_pc(2) = (self%Oxford_pc(2)-0.5)*real(self%nml%Ny)
self%EMsoft_pc(3) = self%Oxford_pc(3)*real(self%nml%Nx)*self%nml%delta 

end subroutine Oxford2EMsoft_

!--------------------------------------------------------------------------
subroutine EMsoft2Bruker_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: EMsoft2Bruker_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! convert pattern center from EMsoft to Bruker convention 

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self

self%Bruker_pc(1) = self%EMsoft_pc(1)/real(self%nml%Nx) + 0.5 
self%Bruker_pc(2) = 0.5-self%EMsoft_pc(2)/real(self%nml%Ny)
self%Bruker_pc(3) = self%EMsoft_pc(3)/real(self%nml%Nx)/self%nml%delta 

end subroutine EMsoft2Bruker_

!--------------------------------------------------------------------------
subroutine Bruker2EMsoft_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: Bruker2EMsoft_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! convert pattern center from Bruker to EMsoft convention 

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)     :: self

self%EMsoft_pc(1) = (self%Bruker_pc(1)-0.5)*real(self%nml%Nx)
self%EMsoft_pc(2) = (0.5-self%Bruker_pc(2))*real(self%nml%Ny)
self%EMsoft_pc(3) = self%Bruker_pc(3)*real(self%nml%Nx)*self%nml%delta 

end subroutine Bruker2EMsoft_

!--------------------------------------------------------------------------
subroutine EBSDpc_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDpc_
!! author: MDG 
!! version: 1.0 
!! date: 09/09/22
!!
!! perform the computations

use mod_EMsoft
use mod_io

IMPLICIT NONE 

class(EBSDpc_T), INTENT(INOUT)          :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 

type(IO_T)                              :: Message 

associate( nml => self%nml )

! do the pattern center conversions
select case(trim(nml%convention)) 
  case('EMsoft')
    self%EMsoft_pc = (/ nml%xstar, nml%ystar, nml%zstar /)
    call self%EMsoft2TSL_()
    call self%EMsoft2Oxford_()
    call self%EMsoft2Bruker_()
  case('TSL')
    self%TSL_pc = (/ nml%xstar, nml%ystar, nml%zstar /)
    call self%TSL2EMsoft_()
    call self%EMsoft2Oxford_()
    call self%EMsoft2Bruker_()
  case('Oxford')
    self%Oxford_pc = (/ nml%xstar, nml%ystar, nml%zstar /)
    call self%Oxford2EMsoft_()
    call self%EMsoft2TSL_()
    call self%EMsoft2Bruker_()
  case('Bruker')
    self%Bruker_pc = (/ nml%xstar, nml%ystar, nml%zstar /)
    call self%Bruker2EMsoft_()
    call self%EMsoft2TSL_()
    call self%EMsoft2Oxford_()
end select

! print the result to the terminal
call Message%printMessage(' EMsoft pattern center coordinates ')
call Message%WriteValue(' [xpc, ypc, L] = ', self%EMsoft_pc, 3)
call Message%printMessage(' --- ')

call Message%printMessage(' TSL pattern center coordinates ')
call Message%WriteValue(' [x*, y*, z*] = ', self%TSL_pc, 3)
call Message%printMessage(' --- ')

call Message%printMessage(' Oxford pattern center coordinates ')
call Message%WriteValue(' [x*, y*, z*] = ', self%Oxford_pc, 3)
call Message%printMessage(' --- ')

call Message%printMessage(' Bruker pattern center coordinates ')
call Message%WriteValue(' [x*, y*, z*] = ', self%Bruker_pc, 3)
call Message%printMessage(' --- ')

end associate

end subroutine EBSDpc_

end module mod_EBSDpc