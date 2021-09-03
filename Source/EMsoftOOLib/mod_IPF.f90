! ###################################################################
! Copyright (c) 2013-2021, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_IPF
  !! author: MDG 
  !! version: 1.0 
  !! date: 09/03/21
  !!
  !! module for inverse pole figure maps

use mod_kinds
use mod_global

IMPLICIT NONE 

! class definition
type, public :: IPF_T
private 
  integer(kind=irg)   :: ipf_wd
  integer(kind=irg)   :: ipf_ht
  character(fnlen)    :: ipf_filename

contains
private 
  procedure, pass(self) :: get_ipf_wd_
  procedure, pass(self) :: get_ipf_ht_
  procedure, pass(self) :: get_ipf_filename_
  procedure, pass(self) :: set_ipf_wd_
  procedure, pass(self) :: set_ipf_ht_
  procedure, pass(self) :: set_ipf_filename_
  procedure, pass(self) :: make_IPFMap_
 
  generic, public :: get_ipf_wd => get_ipf_wd_
  generic, public :: get_ipf_ht => get_ipf_ht_
  generic, public :: get_ipf_filename => get_ipf_filename_
  generic, public :: set_ipf_wd => set_ipf_wd_
  generic, public :: set_ipf_ht => set_ipf_ht_
  generic, public :: set_ipf_filename => set_ipf_filename_
  generic, public :: make_IPFMap => make_IPFMap_

end type IPF_T

! the constructor routine for this class 
interface IPF_T
  module procedure IPF_constructor
end interface IPF_T


contains

!--------------------------------------------------------------------------
function get_ipf_wd_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_wd_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get ipf_wd from the IPF_T class

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out

out = self%nml%ipf_wd

end function get_ipf_wd_

!--------------------------------------------------------------------------
subroutine set_ipf_wd_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_wd_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set ipf_wd in the IPF_T class

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp

self%nml%ipf_wd = inp

end subroutine set_ipf_wd_

!--------------------------------------------------------------------------
function get_ipf_ht_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_ht_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get ipf_ht from the IPF_T class

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out

out = self%nml%ipf_ht

end function get_ipf_ht_

!--------------------------------------------------------------------------
subroutine set_ipf_ht_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_ht_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set ipf_ht in the IPF_T class

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp

self%nml%ipf_ht = inp

end subroutine set_ipf_ht_

!--------------------------------------------------------------------------
function get_ipf_filename_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_filename_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get ipf_filename from the IPF_T class

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)     :: self
character(fnlen)                :: out

out = self%nml%ipf_filename

end function get_ipf_filename_

!--------------------------------------------------------------------------
subroutine set_ipf_filename_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_filename_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set ipf_filename in the IPF_T class

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)    :: inp

self%nml%ipf_filename = inp

end subroutine set_ipf_filename_


!--------------------------------------------------------------------------
subroutine make_IPFMap_() 
!DEC$ ATTRIBUTES DLLEXPORT :: make_IPFMap_
!! author: MDG
!! version: 1.0
!! date: 09/03/21
!!
!! This takes an orientation list in quaternion representation and turns it into 
!! an inverse pole figure map (using one of a few color schemes) that is then stored 
!! in a .tiff or .bmp file.
!!

IMPLICIT NONE 


end subroutine make_IPFMap_



end module mod_IPF