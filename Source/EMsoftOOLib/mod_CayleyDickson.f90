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

module mod_CayleyDickson
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/29/24
  !!
  !! class definition for the Cayley-Dickson construction

use mod_kinds
use mod_global

IMPLICIT NONE 

! class definition
type, public :: CayleyDickson_T
private 

contains
private 
  procedure, pass(self) :: CDconjgs_
  procedure, pass(self) :: CDconjgd_
  procedure, pass(self) :: CDmults_
  procedure, pass(self) :: CDmultd_

  generic, public :: CDconjgs => CDconjgs_, CDconjgd_
  generic, public :: CDmults => CDmults_, CDmultd_

end type CayleyDickson_T

! the constructor routine for this class 
interface CayleyDickson_T
  module procedure CayleyDickson_constructor
end interface CayleyDickson_T

contains

!--------------------------------------------------------------------------
type(CayleyDickson_T) function CayleyDickson_constructor( ) result(CayleyDickson)
!! author: MDG 
!! version: 1.0 
!! date: 02/29/24
!!
!! constructor for the CayleyDickson_T Class; reads the name list 
 
IMPLICIT NONE

end function CayleyDickson_constructor

!--------------------------------------------------------------------------
subroutine CayleyDickson_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 02/29/24
!!
!! destructor for the CayleyDickson_T Class
 
IMPLICIT NONE

type(CayleyDickson_T), INTENT(INOUT)  :: self 

call reportDestructor('CayleyDickson_T')

end subroutine CayleyDickson_destructor

!--------------------------------------------------------------------------
recursive function CDconjgs_(self, n1, ns) result(nout)
!DEC$ ATTRIBUTES DLLEXPORT :: CDconjgs_
!! author: MDG
!! version: 1.0
!! date: 02/29/24
!!
!! conjugate of a higher dimensional number (single precision)

IMPLICIT NONE

class(CayleyDickson_T),INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)          :: ns
real(kind=sgl),INTENT(IN)             :: n1(ns)
real(kind=sgl)                        :: nout(ns)

real(kind=sgl)                        :: a(ns/2), b(ns/2)
integer(kind=irg)                     :: nh

if (ns.eq.1) then
  nout = n1
else
  nh = ns/2
  a = self%CDconjgs_(n1(1:nh),nh)
  b = -n1(nh+1:ns)
  nout = (/ a, b /)
end if

end function CDconjgs_

!--------------------------------------------------------------------------
recursive function CDconjgd_(self, n1, ns) result(nout)
!DEC$ ATTRIBUTES DLLEXPORT :: CDconjgd_
!! author: MDG
!! version: 1.0
!! date: 02/29/24
!!
!! conjugate of a higher dimensional number (double precision)

IMPLICIT NONE

class(CayleyDickson_T),INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)          :: ns
real(kind=dbl),INTENT(IN)             :: n1(ns)
real(kind=dbl)                        :: nout(ns)

real(kind=dbl)                        :: a(ns/2), b(ns/2)
integer(kind=irg)                     :: nh

if (ns.eq.1) then
  nout = n1
else
  nh = ns/2
  a = self%CDconjgd_(n1(1:nh),nh)
  b = -n1(nh+1:ns)
  nout = (/ a, b /)
end if

end function CDconjgd_

!--------------------------------------------------------------------------
recursive function CDmults_(self, n1, n2, ns, split) result(nout)
!DEC$ ATTRIBUTES DLLEXPORT :: CDmults_
!! author: MDG
!! version: 1.0
!! date: 02/29/24
!!
!! multiplication of two higher dimensional numbers using the Cayley-Dickson construction
!! routine can handle split numbers (single precision)

IMPLICIT NONE

class(CayleyDickson_T),INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)          :: ns
real(kind=sgl),INTENT(IN)             :: n1(ns)
real(kind=sgl),INTENT(IN)             :: n2(ns)
logical,INTENT(IN),OPTIONAL           :: split
real(kind=sgl)                        :: nout(ns)

real(kind=sgl)                        :: a(ns/2), b(ns/2), c(ns/2), d(ns/2), dc(ns/2), cc(ns/2)
integer(kind=irg)                     :: nh

nout = 0.0
if (ns.eq.1) then
! here we just multiply two scalars
  nout(1) = n1(1) * n2(1)
else
! we still have a higher dimensional number so we use the Cayley-Dickson formula
  nh = ns/2
  a = n1(1:nh)
  b = n1(nh+1:ns)
  c = n2(1:nh)
  d = n2(nh+1:ns)
  cc = self%CDconjgs_(c,nh)
  dc = self%CDconjgs_(d,nh)
  if (.not.present(split)) then
      nout = (/ self%CDmults_(a,c,nh) - self%CDmults_(dc,b,nh), self%CDmults_(d,a,nh) + self%CDmults_(b,cc,nh) /)
  else if (split.eqv..TRUE.) then
      nout = (/ self%CDmults_(a,c,nh, .TRUE.) + self%CDmults_(dc,b,nh, .TRUE.), & 
                self%CDmults_(d,a,nh, .TRUE.) + self%CDmults_(b,cc,nh, .TRUE.) /)
  end if
end if 

end function CDmults_

!--------------------------------------------------------------------------
! Function: CDmultd
!
!> @author Marc De Graef, Carnegie Mellon University 
!
!> @brief multiplication of two higher dimensional numbers using the Caley-Dickson construction
!
!> @param n1 first input number
!> @param n2 second input number
!> @param n1 dimensionality of the number
!> @param split if present and .TRUE., then use the split multiplication rule
!> @param nout product of the two input numbers 
!
!> @date 03/14/18 MDG 1.0 original
!--------------------------------------------------------------------------
recursive function CDmultd_(self, n1, n2, ns, split) result(nout)
!DEC$ ATTRIBUTES DLLEXPORT :: CDmultd_
!! author: MDG
!! version: 1.0
!! date: 02/29/24
!!
!! multiplication of two higher dimensional numbers using the Cayley-Dickson construction
!! routine can handle split numbers (double precision)

IMPLICIT NONE

class(CayleyDickson_T),INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)          :: ns
real(kind=dbl),INTENT(IN)             :: n1(ns)
real(kind=dbl),INTENT(IN)             :: n2(ns)
logical,INTENT(IN),OPTIONAL           :: split
real(kind=dbl)                        :: nout(ns)

real(kind=dbl)                        :: a(ns/2), b(ns/2), c(ns/2), d(ns/2), dc(ns/2), cc(ns/2)
integer(kind=irg)                     :: nh

nout = 0.D0
if (ns.eq.1) then
! here we just multiply two scalars
  nout(1) = n1(1) * n2(1)
else
! we still have a higher dimensional number so we use the Caley-Dickson formula
  nh = ns/2
  a = n1(1:nh)
  b = n1(nh+1:ns)
  c = n2(1:nh)
  d = n2(nh+1:ns)
  cc = self%CDconjgd_(c,nh)
  dc = self%CDconjgd_(d,nh)
  if (.not.present(split)) then
      nout = (/ self%CDmultd_(a,c,nh) - self%CDmultd_(dc,b,nh), self%CDmultd_(d,a,nh) + self%CDmultd_(b,cc,nh) /)
  else if (split.eqv..TRUE.) then
      nout = (/ self%CDmultd_(a,c,nh, .TRUE.) + self%CDmultd_(dc,b,nh, .TRUE.), & 
                self%CDmultd_(d,a,nh, .TRUE.) + self%CDmultd_(b,cc,nh, .TRUE.) /)
  end if
end if 

end function CDmultd_

end module mod_CayleyDickson