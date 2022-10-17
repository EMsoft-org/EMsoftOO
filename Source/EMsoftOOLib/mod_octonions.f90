! ###################################################################
! Copyright (c) 2013-2022, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_octonions
  !! author: MDG 
  !! version: 1.0 
  !! date: 10/16/22
  !!
  !! class definition for the octonions module
  !!
  !! octonion multiplication is carried out in terms of the member quaternions
  !! both member quaternions are either single or double precision

use mod_kinds
use mod_global
use mod_quaternions, only:Quaternion_T

IMPLICIT NONE 

intrinsic :: conjg, cabs
public :: conjg, cabs

interface conjg
  procedure octconjg_
end interface conjg

interface cabs
  procedure octnorm_
end interface cabs

! class definition
type, public :: octonion_T
private 
  type(Quaternion_T)          :: o(2)
  logical                     :: GBmode = .FALSE.  ! normalizations are different when in grain boundary mode 

contains
private 
  procedure, pass(self) :: octprint_
  procedure, pass(self) :: octadd_
  procedure, pass(self) :: octsubtract_
  procedure, pass(self) :: octmult_
  procedure, pass(self) :: octsmult_
  procedure, pass(self) :: octdivide_
  procedure, pass(self) :: octsdivide_
  procedure, pass(self) :: octinverse_
  procedure, pass(self) :: octconjg_
  procedure, pass(self) :: octnorm_
  procedure, pass(self) :: octnormalize_
  procedure, pass(self) :: octsequal_
  procedure, pass(self) :: getGBmode_
  procedure, pass(self) :: getocts_
  procedure, pass(self) :: getoctd_
  procedure, pass(self) :: setGBmode_
  procedure, pass(self) :: setocts_
  procedure, pass(self) :: setoctd_

  generic, public :: oct_print => octprint_
  generic, public :: operator(+) => octadd_
  generic, public :: operator(-) => octsubtract_
  generic, public :: operator(*) => octmult_
  generic, public :: operator(*) => octsmult_
  generic, public :: operator(/) => octdivide_
  generic, public :: operator(/) => octsdivide_
  generic, public :: octinverse => octinverse_
  generic, public :: octnormalize => octnormalize_
  generic, public :: octsequal => octsequal_
  generic, public :: get_GBmode => getGBmode_
  generic, public :: get_octs => getocts_
  generic, public :: get_octd => getoctd_
  generic, public :: set_GBmode => setGBmode_
  generic, public :: set_octs => setocts_
  generic, public :: set_octd => setoctd_

end type octonion_T

! the constructor routine for this class 
interface octonion_T
  module procedure Octonion_constructor
end interface octonion_T

contains

!--------------------------------------------------------------------------
type(octonion_T) function Octonion_constructor( o, od, GBmode, smode ) result(octonion)
!DEC$ ATTRIBUTES DLLEXPORT :: Octonion_constructor
!! author: MDG 
!! version: 1.0 
!! date: 10/16/22
!!
!! constructor for the octonions_T Class
 
use mod_quaternions

IMPLICIT NONE

real(kind=sgl),OPTIONAL   :: o(8)  
real(kind=dbl),OPTIONAL   :: od(8)  
logical, OPTIONAL         :: GBmode
character(1), OPTIONAL    :: smode

if ((.not.present(o)).and.(.not.present(od))) then
  ! by default we construct a double precision zero octonion
  if (present(smode)) then 
    if (smode.eq.'s') then 
      octonion%o(1) = Quaternion_T( q = (/ 0.0, 0.0, 0.0, 0.0 /) )
      octonion%o(2) = Quaternion_T( q = (/ 0.0, 0.0, 0.0, 0.0 /) )
    else
      octonion%o(1) = Quaternion_T( qd = (/ 0.D0, 0.D0, 0.D0, 0.D0 /) )
      octonion%o(2) = Quaternion_T( qd = (/ 0.D0, 0.D0, 0.D0, 0.D0 /) )
    end if 
  else
    octonion%o(1) = Quaternion_T( qd = (/ 0.D0, 0.D0, 0.D0, 0.D0 /) )
    octonion%o(2) = Quaternion_T( qd = (/ 0.D0, 0.D0, 0.D0, 0.D0 /) )
  end if
else
  if (present(o)) then 
    octonion%o(1) = Quaternion_T( q = o(1:4) )
    octonion%o(2) = Quaternion_T( q = o(5:8) )
  end if 

  if (present(od)) then 
    octonion%o(1) = Quaternion_T( qd = od(1:4) )
    octonion%o(2) = Quaternion_T( qd = od(5:8) )
  end if 
end if

if (present(GBmode)) octonion%GBmode = GBmode

end function Octonion_constructor

!--------------------------------------------------------------------------
subroutine Octonion_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 10/16/22
!!
!! destructor for the octonions_T Class
 
IMPLICIT NONE

type(octonion_T), INTENT(INOUT)  :: self 

call reportDestructor('octonion_T')

end subroutine octonion_destructor

!--------------------------------------------------------------------------
function getGBmode_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getGBmode_
!! author: MDG 
!! version: 1.0 
!! date: 10/16/22
!!
!! get GBmode from the Octonion_T class

IMPLICIT NONE 

class(Octonion_T), INTENT(INOUT)     :: self
logical                              :: out

out = self%GBmode

end function getGBmode_

!--------------------------------------------------------------------------
subroutine setGBmode_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setGBmode_
!! author: MDG 
!! version: 1.0 
!! date: 10/16/22
!!
!! set GBmode in the Octonion_T class

IMPLICIT NONE 

class(Octonion_T), INTENT(INOUT)     :: self
logical, INTENT(IN)                  :: inp

self%GBmode = inp

end subroutine setGBmode_

!--------------------------------------------------------------------------
recursive subroutine octprint_(self, m)
!DEC$ ATTRIBUTES DLLEXPORT :: octprint_
  !! author: MDG
  !! version: 1.0
  !! date: 10/16/22
  !!
  !! print an octonion
 
use mod_quaternions
use mod_io

IMPLICIT NONE

class(octonion_T),intent(in)    :: self
character(*),intent(in)         :: m

type(IO_T)                      :: Message

if (self%o(1)%quat_getprecision().eq.'s') then
  call Message % WriteValue(trim(m), (/ self%o(1)%get_quats(), self%o(2)%get_quats() /), 8, frm="('(',8f16.6,'); precision: '$)")
  call Message % WriteValue('',self%o(1)%quat_getprecision())
else
  call Message % WriteValue(trim(m), (/ self%o(1)%get_quatd(), self%o(2)%get_quatd() /), 8, frm="('(',8f24.14,'); precision: '$)")
  call Message % WriteValue('',self%o(1)%quat_getprecision())
end if

end subroutine octprint_

!--------------------------------------------------------------------------
recursive function getocts_(self) result(qs)
!DEC$ ATTRIBUTES DLLEXPORT :: getocts_
  !! author: MDG
  !! version: 1.0
  !! date: 10/16/22
  !!
  !! return an octonion

use mod_quaternions

IMPLICIT NONE

class(Octonion_T),intent(in)    :: self

real(kind=sgl)                  :: qs(8)

qs = (/ self%o(1)%get_quats(), self%o(2)%get_quats() /)

end function getocts_

!--------------------------------------------------------------------------
recursive function getoctd_(self) result(qd)
!DEC$ ATTRIBUTES DLLEXPORT :: getoctd_
  !! author: MDG
  !! version: 1.0
  !! date: 10/16/22
  !!
  !! return an octonion

use mod_quaternions

IMPLICIT NONE

class(Octonion_T),intent(in)    :: self

real(kind=dbl)                  :: qd(8)

qd = (/ self%o(1)%get_quatd(), self%o(2)%get_quatd() /)

end function getoctd_

!--------------------------------------------------------------------------
recursive subroutine setocts_(self, qs)
!DEC$ ATTRIBUTES DLLEXPORT :: setocts_
  !! author: MDG
  !! version: 1.0
  !! date: 10/16/22
  !!
  !! set an octonion

use mod_quaternions

IMPLICIT NONE

class(Octonion_T),intent(inout)    :: self
real(kind=sgl),intent(in)          :: qs(8)

self%o(1) = Quaternion_T( q = qs(1:4) )
self%o(2) = Quaternion_T( q = qs(5:8) )

end subroutine setocts_

!--------------------------------------------------------------------------
recursive subroutine setoctd_(self, qd)
!DEC$ ATTRIBUTES DLLEXPORT :: setoctd_
  !! author: MDG
  !! version: 1.0
  !! date: 10/16/22
  !!
  !! set an octonion

use mod_quaternions

IMPLICIT NONE

class(Octonion_T),intent(inout)    :: self
real(kind=dbl),intent(in)          :: qd(8)

self%o(1) = Quaternion_T( qd = qd(1:4) )
self%o(2) = Quaternion_T( qd = qd(5:8) )

end subroutine setoctd_

!--------------------------------------------------------------------------
recursive function octsequal_(self, qb) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: octsequal_
  !! author: MDG
  !! version: 1.0
  !! date: 10/17/22
  !!
  !! octonion comparison

IMPLICIT NONE

  class(Octonion_T),intent(in)    :: self, qb
  logical                         :: res

  type(Octonion_T)                :: diff

  real(kind=sgl)                  :: d, eps=1.0e-6
  real(kind=dbl)                  :: dd, epsd=1.0e-12

  res = .TRUE.
  diff = self - qb

  if (self%o(1)%quat_getprecision().eq.'s') then
    d = maxval( abs( diff%getocts_() ) )
    if (d.gt.eps) res = .FALSE.
  else
    dd = maxval( abs( diff%getoctd_() ) )
    if (dd.gt.epsd) res = .FALSE.
  end if

end function octsequal_

!--------------------------------------------------------------------------
recursive function octconjg_(self) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: octconjg_
  !! author: MDG
  !! version: 1.0
  !! date: 10/16/22
  !!
  !! octonion conjugation (single/double precision)

IMPLICIT NONE

class(Octonion_T),intent(in) :: self
type(Octonion_T)             :: qres

real(kind=sgl)               :: a(8)
real(kind=dbl)               :: b(8)

if (self%o(1)%quat_getprecision().eq.'s') then
  a = self%get_octs()
  a(2:8) = -a(2:8)
  call qres%set_octs( a )
else
  b = self%get_octd()
  b(2:8) = -b(2:8)
  call qres%set_octd( b )
end if

end function octconjg_

!--------------------------------------------------------------------------
recursive function octadd_(self, y) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: octadd_
  !! author: MDG
  !! version: 1.0
  !! date: 10/16/22
  !!
  !! octonion addition (single/double precision)

IMPLICIT NONE

class(Octonion_T),intent(in) :: self, y
type(Octonion_T)             :: qres

if (self%o(1)%quat_getprecision().eq.'s') then
  call qres%set_octs( self%get_octs() + y%get_octs() )
else
  call qres%set_octd( self%get_octd() + y%get_octd() )
end if

end function octadd_

!--------------------------------------------------------------------------
recursive function octsubtract_(self, y) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: octsubtract_
  !! author: MDG
  !! version: 1.0
  !! date: 10/16/22
  !!
  !! octonion subtraction (single/double precision)

IMPLICIT NONE

class(Octonion_T),intent(in) :: self, y
type(Octonion_T)             :: qres

if (self%o(1)%quat_getprecision().eq.'s') then
  call qres%set_octs( self%get_octs() - y%get_octs() )
else
  call qres%set_octd( self%get_octd() - y%get_octd() )
end if

end function octsubtract_

!--------------------------------------------------------------------------
recursive function octmult_(self, y) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: octmult_
  !! author: MDG
  !! version: 1.0
  !! date: 10/16/22
  !!
  !! octonion multiplication   (single/double precision)

use mod_quaternions

IMPLICIT NONE

class(Octonion_T),intent(in) :: self, y
type(Octonion_T)             :: qres

type(Quaternion_T)           :: qa, qb, qc, qd, a, b  

! extract the quaternions
qa = self%o(1)
qb = self%o(2)
qc = y%o(1)
qd = y%o(2)

! carry out the quaternion multiplications according to the Cayley-Dickson scheme 
a = qa * qc - qd%qconjg() * qb
b = qd * qa + qb * qc%qconjg()

! and collect the results in the output octonion
if (self%o(1)%quat_getprecision().eq.'s') then
  call qres%set_octs( (/ a%get_quats(), b%get_quats() /) )
else
  call qres%set_octd( (/ a%get_quatd(), b%get_quatd() /) )
end if

end function octmult_

!--------------------------------------------------------------------------
recursive function octsmult_(self, y) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: octsmult_
  !! author: MDG
  !! version: 1.0
  !! date: 10/17/22
  !!
  !! scalar-octonion multiplication   (single/double precision)

use mod_quaternions

IMPLICIT NONE

class(Octonion_T),intent(in) :: self
real(kind=dbl),intent(in)    :: y
type(Octonion_T)             :: qres

type(Quaternion_T)           :: qa, qb

! extract the quaternions and multiply by the scalar
! and collect the results in the output octonion
if (self%o(1)%quat_getprecision().eq.'s') then
  qa = self%o(1) * real(y)
  qb = self%o(2) * real(y)
  call qres%set_octs( (/ qa%get_quats(), qb%get_quats() /) )
else
  qa = self%o(1) * y
  qb = self%o(2) * y
  call qres%set_octd( (/ qa%get_quatd(), qb%get_quatd() /) )
end if

end function octsmult_

!--------------------------------------------------------------------------
recursive function octsdivide_(self, y) result(qres)
!DEC$ ATTRIBUTES DLLEXPORT :: octsdivide_
  !! author: MDG
  !! version: 1.0
  !! date: 10/17/22
  !!
  !! scalar-octonion multiplication   (single/double precision)

use mod_quaternions
use mod_io 

IMPLICIT NONE

class(Octonion_T),intent(in) :: self
real(kind=dbl),intent(in)    :: y
type(Octonion_T)             :: qres

type(Quaternion_T)           :: qa, qb
type(IO_T)                   :: Message

if (y.eq.0.D0) then 
  call Message%printError('octsdivide','attmpting to divide by zero')
end if 

! extract the quaternions and divide by the scalar
! and collect the results in the output octonion
if (self%o(1)%quat_getprecision().eq.'s') then
  qa = self%o(1) / real(y)
  qb = self%o(2) / real(y)
  call qres%set_octs( (/ qa%get_quats(), qb%get_quats() /) )
else
  qa = self%o(1) / y
  qb = self%o(2) / y
  call qres%set_octd( (/ qa%get_quatd(), qb%get_quatd() /) )
end if

end function octsdivide_

!--------------------------------------------------------------------------
recursive function octnorm_(self) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: octnorm_
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! octonion norm (extends intrinsic routine abs)

IMPLICIT NONE

  class(Octonion_T),intent(in)   :: self
  real(kind=dbl)                 :: res

  real(kind=sgl)                 :: n, a(8)
  real(kind=dbl)                 :: nd, resd, b(8)

if (self%o(1)%quat_getprecision().eq.'s') then
  a = self%getocts_()
  n = sum(a*a)
  resd = dsqrt( dble(n) )
  res = dble(sngl(resd))
else
  b = self%getoctd_()
  nd = sum(b*b)
  res = dsqrt( nd )
end if

end function octnorm_

!--------------------------------------------------------------------------
recursive subroutine octnormalize_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: octnormalize_
  !! author: MDG
  !! version: 1.0
  !! date: 10/17/22
  !!
  !! octonion normalization (depends on GBmode)

IMPLICIT NONE

  class(Octonion_T),intent(inout):: self

  real(kind=dbl)                 :: a(8), b
  real(kind=dbl),parameter       :: s2 = 1.4142135623730951D0  ! sqrt(2)

b = self%octnorm_()

if (self%GBmode.eqv..TRUE.) then  ! we have a Grain Boundary octionion with two unit quaternions
  if (self%o(1)%quat_getprecision().eq.'s') then
    a = self%getocts_() 
    b = sqrt(sum(a(1:4)**2))
    a(1:4) = a(1:4)/b
    b = sqrt(sum(a(5:8)**2))
    a(5:8) = a(5:8)/b
    a = a / s2
    call self%setocts_( real(a) )
  else 
    a = self%getoctd_()
    b = sqrt(sum(a(1:4)**2))
    a(1:4) = a(1:4)/b
    b = sqrt(sum(a(5:8)**2))
    a(5:8) = a(5:8)/b
    a = a / s2
    call self%setoctd_( a )
  end if 
else ! a regular octonion
  if (self%o(1)%quat_getprecision().eq.'s') then
    a = self%getocts_() / b
    call self%setocts_( real(a) )
  else 
    a = self%getoctd_() / b
    call self%setoctd_( a )
  end if 
end if

end subroutine octnormalize_

!--------------------------------------------------------------------------
recursive function octinverse_(self) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: octinverse_ 
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! octonion norm (extends intrinsic routine abs)

IMPLICIT NONE

  class(Octonion_T),intent(in)   :: self
  type(Octonion_T)               :: res

  real(kind=dbl)                 :: b
  type(Octonion_T)               :: o

b = self%octnorm_()
b = b*b
o = self%octconjg_()
res = o / b 

end function octinverse_

!--------------------------------------------------------------------------
recursive function octdivide_(self,y) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: octdivide_ 
  !! author: MDG
  !! version: 1.0
  !! date: 01/06/20
  !!
  !! octonion division

IMPLICIT NONE

  class(Octonion_T),intent(in)   :: self, y
  type(Octonion_T)               :: res

res = self * y%octinverse_()

end function octdivide_





end module mod_octonions