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

module mod_octonions
  !! author: MDG 
  !! version: 1.0 
  !! date: 07/16/25
  !!
  !! basic octonion arithmetic (very loosely based on a version created by ChatGPT4, 
  !! but object oriented this time); operations are defined for single octonions and 
  !! for arrays of octonions. This class largely parallels the Quaternion_T class.

use mod_kinds
use mod_global

IMPLICIT NONE 

! we overload the conjg, cabs, and .eq. intrinsics
intrinsic :: conjg, cabs
public :: conjg, cabs

interface conjg
  procedure o_conjugate_
  procedure o_arrayconjugate_
end interface conjg

interface cabs
  procedure o_norm_
  procedure o_arraynorm_
end interface cabs

interface operator(.eq.)
  procedure octsequal
end interface

type, public :: Octonion_T
  real(kind=sgl)  :: o(8)
  real(kind=dbl)  :: od(8)
  character(1)    :: s

contains
private

! octonion arithmetic routines
  procedure, pass(self) :: o_add_
  procedure, pass(self) :: o_subtract_
  procedure, pass(self) :: o_multiply_
  procedure, pass(self) :: o_smultiply_
  procedure, pass(self) :: o_divide_
  procedure, pass(self) :: o_sdivide_
  procedure, pass(self) :: o_conjugate_
  procedure, pass(self) :: o_norm_
  procedure, pass(self),public :: o_normalize_
  procedure, pass(self) :: o_exp_
  procedure, pass(self) :: o_log_
  procedure, pass(self) :: o_innerproduct_
  procedure, pass(self) :: getocts_
  procedure, pass(self) :: getoctd_
  procedure, pass(self),public :: setocts_
  procedure, pass(self),public :: setoctd_
  procedure, pass(self),public :: octsequal
  procedure, pass(self) :: octprint_
 

! overload the basic operators
  generic, public :: operator(+) => o_add_
  generic, public :: operator(-) => o_subtract_
  generic, public :: operator(*) => o_multiply_
  generic, public :: operator(*) => o_smultiply_
  generic, public :: operator(/) => o_divide_
  generic, public :: o_sdivide => o_sdivide_
  generic, public :: o_normalize => o_normalize_
  generic, public :: o_exp => o_exp_
  generic, public :: o_log => o_log_
  generic, public :: o_innerproduct => o_innerproduct_
  generic, public :: get_octs => getocts_
  generic, public :: get_octd => getoctd_
  generic, public :: set_octs => setocts_
  generic, public :: set_octd => setoctd_
  generic, public :: oct_print => octprint_
 
end type Octonion_T

! the constructor routine for this class 
interface Octonion_T
  module procedure Octonion_constructor
end interface Octonion_T


! next we define the Octonion Array class
type, public :: OctonionArray_T
    integer(kind=irg)            :: n
    integer(kind=irg)            :: nthreads
    real(kind=sgl), allocatable  :: o(:,:)
    real(kind=dbl), allocatable  :: od(:,:)
    character(1)                 :: s

  contains
  private
! quaternion IO routines
    procedure, pass(self) :: o_arrayprint_
! quaternion arithmetic routines
    procedure, pass(self) :: o_arrayadd_
    procedure, pass(self) :: o_arraysubtract_
    procedure, pass(self) :: o_arraymult_
    procedure, pass(self) :: o_arraysmult_
    procedure, pass(self) :: o_arraydiv_
    procedure, pass(self) :: o_arrayconjugate_
    procedure, pass(self) :: o_arraynorm_
    procedure, pass(self) :: o_arraynormalize_
! routines with two or more input quaternion arrays
    procedure, pass(self) :: o_arrayinnerproduct_
! miscellaneous routines
    procedure, pass(self) :: extractfromOctonionArray_
    procedure, pass(self) :: insertOctintoArray_
    procedure, pass(self) :: getOnumber_
    procedure, pass(self) :: deleteArray_
    procedure, pass(self) :: writeArraytoFile_

! generics
    generic, public :: oct_arrayprint => o_arrayprint_
    generic, public :: operator(+) => o_arrayadd_
    generic, public :: operator(-) => o_arraysubtract_
    generic, public :: operator(*) => o_arraymult_
    generic, public :: operator(*) => o_arraysmult_
    generic, public :: operator(/) => o_arraydiv_
    generic, public :: oct_normalize => o_arraynormalize_
    generic, public :: oct_innerproduct => o_arrayinnerproduct_
    generic, public :: getOctfromArray => extractfromOctonionArray_
    generic, public :: insertOctinArray => insertOctintoArray_
    generic, public :: getOnumber => getOnumber_
    generic, public :: deleteArray => deleteArray_
    generic, public :: writeArraytoFile => writeArraytoFile_

end type OctonionArray_T

PRIVATE:: insertOctintoArray_

! the constructor routine for this class 
interface OctonionArray_T
  module procedure OctonionArray_constructor
end interface OctonionArray_T

contains

!--------------------------------------------------------------------------
type(Octonion_T) function Octonion_constructor( o, od ) result(Oct)
!DEC$ ATTRIBUTES DLLEXPORT :: Octonion_constructor
  !! author: MDG
  !! version: 1.0
  !! date: 07/16/25
  !!
  !! constructor for the Octonion_T Class

IMPLICIT NONE

real(kind=sgl), INTENT(IN), OPTIONAL      :: o(8)
real(kind=dbl), INTENT(IN), OPTIONAL      :: od(8)

! fill in one or the other octonion
if ((.not.present(o)).and.(.not.present(od))) then
  Oct % o(1:8) = 0.0
  Oct % od(1:8) = 0.D0
  Oct % s = 's'
else
    if (present(o)) then
      Oct % o = o
      Oct % od(1:8) = 0.D0
      Oct % s = 's'
    end if

    if (present(od)) then
      Oct % o(1:8) = 0.0
      Oct % od = od
      Oct % s = 'd'
    end if
end if

end function Octonion_constructor

!--------------------------------------------------------------------------
subroutine Octonion_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: Octonion_destructor
!! author: MDG
!! version: 1.0
!! date: 07/17/25
!!
!! destructor for the Octonion_T Class

IMPLICIT NONE

type(Octonion_T), INTENT(INOUT)     :: self

call reportDestructor('Octonion_T')

end subroutine Octonion_destructor

!--------------------------------------------------------------------------
type(OctonionArray_T) function OctonionArray_constructor( n, nthreads, o, od, s ) result(OctArray)
!DEC$ ATTRIBUTES DLLEXPORT :: OctonionArray_constructor
  !! author: MDG
  !! version: 1.0
  !! date: 07/17/25
  !!
  !! constructor for the OctonionArray Class
  !!
  !! either call with parameters n and s
  !! or with n and either one of q or qd

IMPLICIT NONE

  integer(kind=irg), INTENT(IN)             :: n
  integer(kind=irg), INTENT(IN), OPTIONAL   :: nthreads
  real(kind=sgl), INTENT(IN), OPTIONAL      :: o(4,n)
  real(kind=dbl), INTENT(IN), OPTIONAL      :: od(4,n)
  character(1), INTENT(IN), OPTIONAL        :: s

! OpenMP threads
  OctArray % nthreads = 0
  if (present(nthreads)) OctArray % nthreads = nthreads

! are we declaring just an empty variable with no entries, but with a given precision ?
  if ( present(s) .and. (.not.present(o)) .and. (.not.present(od)) ) then
    OctArray % n = n
    OctArray % s = s
    if (s.eq.'s') then
      allocate(OctArray % o(8,n))
      OctArray % o = 0.0
    else
      allocate(OctArray % od(8,n))
      OctArray % od = 0.D0
    end if
    return
  end if

! single precision
  if (present(o)) then
    allocate(OctArray % o(8,n))
    OctArray % n = n
    OctArray % o = o
    OctArray % s = 's'
  end if

! double precision
  if (present(od)) then
    allocate(OctArray % od(8,n))
    OctArray % n = n
    OctArray % od = od
    OctArray % s = 'd'
  end if

end function OctonionArray_constructor

!--------------------------------------------------------------------------
subroutine OctonionArray_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: OctonionArray_destructor
!! author: MDG
!! version: 1.0
!! date: 07/17/25
!!
!! destructor for the OctonionArray_T Class

IMPLICIT NONE

type(OctonionArray_T), INTENT(INOUT)     :: self

call reportDestructor('OctonionArray_T')

if (allocated(self%o)) deallocate(self%o)
if (allocated(self%od)) deallocate(self%od)

end subroutine OctonionArray_destructor

!--------------------------------------------------------------------------
recursive function o_add_(self, b) result(c)
!DEC$ ATTRIBUTES DLLEXPORT :: o_add_
!! author: MDG 
!! version: 1.0 
!! date: 07/16/25
!!
!! basic octonion addition

class(Octonion_T), INTENT(IN)       :: self
type(Octonion_T), INTENT(IN)        :: b
type(Octonion_T)                    :: c

if (self%s.eq.'s') then 
  c%o = self%o + b%o
  c%s = 's'
else
  c%od = self%od + b%od
  c%s = 'd'
end if

end function o_add_

!--------------------------------------------------------------------------
recursive function o_arrayadd_(self, y) result(oct)
!DEC$ ATTRIBUTES DLLEXPORT :: o_arrayadd_
  !! author: MDG
  !! version: 1.0
  !! date: 07/17/25
  !!
  !! octonion array addition (single/double precision)

use mod_io

IMPLICIT NONE

class(OctonionArray_T),intent(in)   :: self, y
type(OctonionArray_T)               :: oct

type(IO_T)                          :: Message
integer(kind=irg)                   :: sz(2)

! test to make sure that both arrays have the same number of quaternions
if (self%n.ne.y%n) then
  call Message%printError('o_arrayadd_','input arrays must have the same number of quaternions')
end if

if (self%s.ne.y%s) then
  call Message%printError('o_arrayadd_','input arrays must have the same precision')
end if

oct%n = self%n
oct%s = self%s
oct%nthreads = self%nthreads

if (self%s.eq.'s') then
! if the quaternion array is already allocated, check to make sure it has the right dimensions
  if (allocated(oct%o)) then
    sz = shape(oct%o)
    if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(oct%o)
  end if
  allocate(oct%o(8,self%n))
  oct%o = self%o + y%o
else
! if the quaternion array is already allocated, check to make sure it has the right dimensions
  if (allocated(oct%od)) then
    sz = shape(oct%od)
    if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(oct%od)
  end if
  allocate(oct%od(8,self%n))
  oct%od = self%od + y%od
end if

end function o_arrayadd_

!--------------------------------------------------------------------------
recursive function o_subtract_(self, b) result(c)
!DEC$ ATTRIBUTES DLLEXPORT :: o_subtract_
!! author: MDG 
!! version: 1.0 
!! date: 07/16/25
!!
!!  subtract two octonions 

class(Octonion_T), INTENT(IN)     :: self
type(Octonion_T), INTENT(IN)      :: b
type(Octonion_T)                  :: c

if (self%s.eq.'s') then 
  c%o = self%o - b%o
  c%s = 's'
else
  c%od = self%od - b%od
  c%s = 'd'
end if

end function o_subtract_

!--------------------------------------------------------------------------
recursive function o_arraysubtract_(self, y) result(oct)
!DEC$ ATTRIBUTES DLLEXPORT :: o_arraysubtract_
  !! author: MDG
  !! version: 1.0
  !! date: 07/17/25
  !!
  !! octonion array subtraction (single/double precision)

use mod_io

IMPLICIT NONE

class(OctonionArray_T),intent(in)   :: self, y
type(OctonionArray_T)               :: oct

type(IO_T)                          :: Message
integer(kind=irg)                   :: sz(2)

! test to make sure that both arrays have the same number of quaternions
if (self%n.ne.y%n) then
  call Message%printError('o_arraysubtract_','input arrays must have the same number of quaternions')
end if

if (self%s.ne.y%s) then
  call Message%printError('o_arraysubtract_','input arrays must have the same precision')
end if

oct%n = self%n
oct%s = self%s
oct%nthreads = self%nthreads

if (self%s.eq.'s') then
! if the quaternion array is already allocated, check to make sure it has the right dimensions
  if (allocated(oct%o)) then
    sz = shape(oct%o)
    if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(oct%o)
  end if
  allocate(oct%o(8,self%n))
  oct%o = self%o - y%o
else
! if the quaternion array is already allocated, check to make sure it has the right dimensions
  if (allocated(oct%od)) then
    sz = shape(oct%od)
    if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(oct%od)
  end if
  allocate(oct%od(8,self%n))
  oct%od = self%od - y%od
end if

end function o_arraysubtract_

!--------------------------------------------------------------------------
recursive function o_multiply_(self, z) result(c)
!DEC$ ATTRIBUTES DLLEXPORT :: o_multiply_
!! author: MDG 
!! version: 1.0 
!! date: 07/16/25
!!
!!  octonion multiplication (computation performed in double precision)

class(Octonion_T), INTENT(IN)   :: self
type(Octonion_T), INTENT(IN)    :: z
type(Octonion_T)                :: c

real(kind=dbl)                  :: A(8), B(8), D(8)

if (self%s.eq.'s') then 
  A = dble(self%o)
  B = dble(z%o)
else
  A = self%od
  B = z%od
end if 

D(1) = A(1)*B(1) - sum(A(2:8)*B(2:8))
D(2) = A(1)*B(2) + A(2)*B(1) + A(3)*B(4) - A(4)*B(3) + A(5)*B(6) - A(6)*B(5) - A(7)*B(8) + A(8)*B(7)
D(3) = A(1)*B(3) - A(2)*B(4) + A(3)*B(1) + A(4)*B(2) + A(5)*B(7) + A(6)*B(8) - A(7)*B(5) - A(8)*B(6)
D(4) = A(1)*B(4) + A(2)*B(3) - A(3)*B(2) + A(4)*B(1) + A(5)*B(8) - A(6)*B(7) + A(7)*B(6) - A(8)*B(5)
D(5) = A(1)*B(5) - A(2)*B(6) - A(3)*B(7) - A(4)*B(8) + A(5)*B(1) + A(6)*B(2) + A(7)*B(3) + A(8)*B(4)
D(6) = A(1)*B(6) + A(2)*B(5) - A(3)*B(8) + A(4)*B(7) - A(5)*B(2) + A(6)*B(1) - A(7)*B(4) + A(8)*B(3)
D(7) = A(1)*B(7) + A(2)*B(8) + A(3)*B(5) - A(4)*B(6) - A(5)*B(3) + A(6)*B(4) + A(7)*B(1) - A(8)*B(2)
D(8) = A(1)*B(8) - A(2)*B(7) + A(3)*B(6) + A(4)*B(5) - A(5)*B(4) - A(6)*B(3) + A(7)*B(2) + A(8)*B(1)

if (self%s.eq.'s') then 
  c%o = real(D)
  c%s = 's'
else
  c%od = D 
  c%s = 'd'
end if 

end function o_multiply_

!--------------------------------------------------------------------------
recursive function o_arraymult_(self, y) result(oct)
!DEC$ ATTRIBUTES DLLEXPORT :: o_arraymult_
  !! author: MDG
  !! version: 1.0
  !! date: 07/17/25
  !!
  !! octonion array multiplication (single/double precision)

use mod_io

IMPLICIT NONE

class(OctonionArray_T),intent(in)   :: self, y
type(OctonionArray_T)               :: oct

type(IO_T)                          :: Message

type(Octonion_T)                    :: o1, o2 
integer(kind=irg)                   :: sz(2), i

! test to make sure that both arrays have the same number of quaternions
if (self%n.ne.y%n) then
  call Message%printError('o_arraymult_','input arrays must have the same number of quaternions')
end if

if (self%s.ne.y%s) then
  call Message%printError('o_arraymult_','input arrays must have the same precision')
end if

oct%n = self%n
oct%s = self%s
oct%nthreads = self%nthreads

if (self%s.eq.'s') then
! if the quaternion array is already allocated, check to make sure it has the right dimensions
  if (allocated(oct%o)) then
    sz = shape(oct%o)
    if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(oct%o)
  end if
  allocate(oct%o(8,self%n))
else
! if the quaternion array is already allocated, check to make sure it has the right dimensions
  if (allocated(oct%od)) then
    sz = shape(oct%od)
    if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(oct%od)
  end if
  allocate(oct%od(8,self%n))
end if

do i=1,self%n 
  o1 = self%extractfromOctonionArray_(i)
  o2 = y%extractfromOctonionArray_(i)
  o1 = o1 * o2
  call oct%insertOctintoArray_(i, o1 )
end do

end function o_arraymult_

!--------------------------------------------------------------------------
recursive function o_smultiply_(self, s) result(c)
!DEC$ ATTRIBUTES DLLEXPORT :: o_smultiply_
!! author: MDG 
!! version: 1.0 
!! date: 07/16/25
!!
!!  octonion scalar multiplication (computation performed in double precision)

class(Octonion_T), INTENT(IN)   :: self
real(kind=dbl), INTENT(IN)      :: s  
type(Octonion_T)                :: c

if (self%s.eq.'s') then 
  c%o = real(s * self%o)
  c%s = 's'
else
  c%od = s * self%od 
  c%s = 'd'
end if 

end function o_smultiply_

!--------------------------------------------------------------------------
recursive function o_arraysmult_(self, s) result(oct)
!DEC$ ATTRIBUTES DLLEXPORT :: o_arraysmult_
  !! author: MDG
  !! version: 1.0
  !! date: 07/17/25
  !!
  !! octonion array scalar multiplication (single/double precision)

use mod_io

IMPLICIT NONE

class(OctonionArray_T),INTENT(IN)   :: self
real(kind=dbl),INTENT(IN)           :: s
type(OctonionArray_T)               :: oct

type(IO_T)                          :: Message

type(Octonion_T)                    :: o1, o2
integer(kind=irg)                   :: sz(2), i

oct%n = self%n
oct%s = self%s
oct%nthreads = self%nthreads

if (self%s.eq.'s') then
! if the quaternion array is already allocated, check to make sure it has the right dimensions
  if (allocated(oct%o)) then
    sz = shape(oct%o)
    if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(oct%o)
  end if
  allocate(oct%o(8,self%n))
else
! if the quaternion array is already allocated, check to make sure it has the right dimensions
  if (allocated(oct%od)) then
    sz = shape(oct%od)
    if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(oct%od)
  end if
  allocate(oct%od(8,self%n))
end if

do i=1,self%n 
  o1 = self%extractfromOctonionArray_(i)
  o2 = o1 * s
  call oct%insertOctintoArray_(i, o2)
end do

end function o_arraysmult_

!--------------------------------------------------------------------------
function o_divide_(self, b) result(c)
!DEC$ ATTRIBUTES DLLEXPORT :: o_divide_
!! author: MDG 
!! version: 1.0 
!! date: 07/16/25
!!
!!  octonion divison

use mod_io 

class(Octonion_T), INTENT(IN)     :: self 
type(Octonion_T), INTENT(IN)      :: b
type(Octonion_T)                  :: c

type(IO_T)                        :: Message 

type(Octonion_T)                  :: b_conj
real(kind=dbl)                    :: b_norm2 = 1.D0
real(kind=sgl)                    :: bs_norm2 = 1.0

b_conj = o_conjugate_(b)
if (b%s.eq.'s') then
  bs_norm2 = sum(b%o**2)
else
  b_norm2 = sum(b%od**2)
end if 

if (b_norm2 == 0.0D0) call Message%printError('o_divide_',' Division by zero octonion')
if (bs_norm2 == 0.0) call Message%printError('o_divide_',' Division by zero octonion')

c = self%o_multiply_( b_conj )

if (b%s.eq.'s') then
  c%o = c%o / bs_norm2
  c%s = 's'
else
  c%od = c%od / b_norm2
  c%s = 'd'
end if

end function o_divide_

!--------------------------------------------------------------------------
recursive function o_arraydiv_(self, y) result(oct)
!DEC$ ATTRIBUTES DLLEXPORT :: o_arraydiv_
  !! author: MDG
  !! version: 1.0
  !! date: 07/17/25
  !!
  !! octonion array division (single/double precision)

use mod_io

IMPLICIT NONE

class(OctonionArray_T),intent(in)   :: self, y
type(OctonionArray_T)               :: oct

type(IO_T)                          :: Message

type(Octonion_T)                    :: o1, o2 
integer(kind=irg)                   :: sz(2), i

! test to make sure that both arrays have the same number of quaternions
if (self%n.ne.y%n) then
  call Message%printError('o_arraydiv_','input arrays must have the same number of quaternions')
end if

if (self%s.ne.y%s) then
  call Message%printError('o_arraydiv_','input arrays must have the same precision')
end if

oct%n = self%n
oct%s = self%s
oct%nthreads = self%nthreads

if (self%s.eq.'s') then
! if the quaternion array is already allocated, check to make sure it has the right dimensions
  if (allocated(oct%o)) then
    sz = shape(oct%o)
    if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(oct%o)
  end if
  allocate(oct%o(8,self%n))
else
! if the quaternion array is already allocated, check to make sure it has the right dimensions
  if (allocated(oct%od)) then
    sz = shape(oct%od)
    if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(oct%od)
  end if
  allocate(oct%od(8,self%n))
end if

do i=1,self%n 
  o1 = self%extractfromOctonionArray_(i)
  o2 = y%extractfromOctonionArray_(i)
  o1 = o1 / o2
  call oct%insertOctintoArray_(i, o1)
end do

end function o_arraydiv_

!--------------------------------------------------------------------------
subroutine o_sdivide_(self, s)
!DEC$ ATTRIBUTES DLLEXPORT :: o_sdivide_
!! author: MDG 
!! version: 1.0 
!! date: 07/16/25
!!
!!  octonion scalar divison (expects a double precision scalar)

use mod_io 

class(Octonion_T), INTENT(INOUT)    :: self
real(kind=dbl), INTENT(IN)          :: s

type(IO_T)                          :: Message

if (s.eq.0.D0) call Message%printError('o_sdivide_',' can not divide octonion by zero ')

if (self%s.eq.'s') then 
  self%o = self%o / real(s)
else
  self%od = self%od / s
end if

end subroutine o_sdivide_

!--------------------------------------------------------------------------
recursive function octsequal(self, b) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: octsequal
  !! author: MDG
  !! version: 1.0
  !! date: 07/16/25
  !!
  !! octonion comparison (double precision)

IMPLICIT NONE

class(Octonion_T),INTENT(IN)      :: self
type(Octonion_T), INTENT(IN)      :: b
logical                           :: res

type(Octonion_T)                  :: diff

real(kind=sgl)                    :: d, eps=1.0e-6
real(kind=dbl)                    :: dd, epsd=1.0e-12

res = .TRUE.
diff = self%o_subtract_(b)

if (self%s.eq.'s') then
  d = maxval( abs( diff%o(:) ) )
  if (d.gt.eps) res = .FALSE.
else
  dd = maxval( abs( diff%od(:) ) )
  if (dd.gt.epsd) res = .FALSE.
end if

end function octsequal

!--------------------------------------------------------------------------
function o_conjugate_(self) result(c)
!DEC$ ATTRIBUTES DLLEXPORT :: o_conjugate_
!! author: MDG 
!! version: 1.0 
!! date: 07/16/25
!!
!! conjugate of an octonion 

class(Octonion_T), INTENT(IN)     :: self
type(Octonion_T)                  :: c

if (self%s.eq.'s') then 
  c%o(1) = self%o(1)
  c%o(2:8) = -self%o(2:8)
  c%s = 's'
else
  c%od(1) = self%od(1)
  c%od(2:8) = -self%od(2:8)
  c%s = 'd'
end if

end function o_conjugate_

!--------------------------------------------------------------------------
function o_arrayconjugate_(self) result(c)
!DEC$ ATTRIBUTES DLLEXPORT :: o_arrayconjugate_
!! author: MDG 
!! version: 1.0 
!! date: 07/16/25
!!
!! conjugate an octonion array

class(OctonionArray_T), INTENT(IN)     :: self
type(OctonionArray_T)                  :: c

type(Octonion_T)                       :: o1, o2
integer(kind=irg)                      :: sz(2), i 

c%n = self%n
c%s = self%s
c%nthreads = self%nthreads

if (self%s.eq.'s') then
! if the quaternion array is already allocated, check to make sure it has the right dimensions
  if (allocated(c%o)) then
    sz = shape(c%o)
    if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(c%o)
  end if
  allocate(c%o(8,self%n))
else
! if the quaternion array is already allocated, check to make sure it has the right dimensions
  if (allocated(c%od)) then
    sz = shape(c%od)
    if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(c%od)
  end if
  allocate(c%od(8,self%n))
end if

do i=1,self%n 
  o1 = self%extractfromOctonionArray_(i)
  o2 = conjg(o1)
  call c%insertOctintoArray_(i, o2 )
end do

end function o_arrayconjugate_

!--------------------------------------------------------------------------
function o_norm_(self) result(n)
!DEC$ ATTRIBUTES DLLEXPORT :: o_norm_
!! author: MDG 
!! version: 1.0 
!! date: 07/16/25
!!
!! norm of an octonion always returned as a double precision parameter but 
!! adjusted to single precision when necessary

class(Octonion_T), INTENT(IN)   :: self
real(kind=dbl)                  :: n

if (self%s.eq.'s') then 
  n = dble(real(sqrt(sum(self%o**2))))
else
  n = sqrt(sum(self%od**2))
end if

end function o_norm_

!--------------------------------------------------------------------------
function o_arraynorm_(self) result(c)
!DEC$ ATTRIBUTES DLLEXPORT :: o_arraynorm_
!! author: MDG 
!! version: 1.0 
!! date: 07/16/25
!!
!! conjugate an octonion array

class(OctonionArray_T), INTENT(IN)     :: self
real(kind=dbl),allocatable             :: c(:)

type(Octonion_T)                       :: o1
integer(kind=irg)                      :: i 

if (allocated(c)) then
  deallocate(c)
  allocate(c(self%n))
end if

do i=1,self%n 
  o1 = self%extractfromOctonionArray_(i)
  c(i) = cabs(o1)
end do

end function o_arraynorm_

!--------------------------------------------------------------------------
subroutine o_normalize_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: o_normalize_
!! author: MDG 
!! version: 1.0 
!! date: 07/16/25
!!
!! normalize an octonion 

use mod_io 

class(Octonion_T), INTENT(INOUT)    :: self

type(IO_T)                          :: Message 

real(kind=dbl)                      :: n

n = self%o_norm_()

if (n == 0.0D0) then 
  call self%octprint_()
  call Message%printError('o_normalize',' Cannot normalize zero octonion') 
end if 

call self%o_sdivide_(n)

end subroutine o_normalize_

!--------------------------------------------------------------------------
subroutine o_arraynormalize_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: o_arraynormalize_
!! author: MDG 
!! version: 1.0 
!! date: 07/16/25
!!
!! normalize an octonion array

use mod_io 

class(OctonionArray_T), INTENT(INOUT)     :: self

type(IO_T)                                :: Message 

type(Octonion_T)                          :: o1
real(kind=dbl)                            :: c
integer(kind=irg)                         :: i, io_int(1)

do i=1,self%n 
  o1 = self%extractfromOctonionArray_(i)
  c = cabs(o1)
  if (c == 0.0D0) then 
    io_int(1) = i
    call Message%writeValue('o_arraynormalize Warning: Cannot normalize octonion ', io_int, 1) 
  else
    call o1%o_sdivide_(c)
    call self%insertOctintoArray_(i, o1)
  end if 
end do

end subroutine o_arraynormalize_

!--------------------------------------------------------------------------
function o_exp_(self) result(c)
!DEC$ ATTRIBUTES DLLEXPORT :: o_exp
!! author: MDG 
!! version: 1.0 
!! date: 07/16/25
!!
!!  exponential of an octonion (carried out in double precision)

class(Octonion_T), INTENT(IN)   :: self
type(Octonion_T)                :: c

real(kind=dbl)                  :: scalar, vnorm, exp_scalar, cos_vnorm, sin_vnorm
real(kind=dbl)                  :: v(7)

if (self%s.eq.'s') then 
  scalar = dble(self%o(1))
  v = dble(self%o(2:8))
else
  scalar = self%od(1)
  v = self%od(2:8)
end if 

vnorm = sqrt(sum(v**2))

exp_scalar = exp(scalar)

if (vnorm == 0.0D0) then
  if (self%s.eq.'s') then 
    c%o(1) = real(exp_scalar)
    c%o(2:8) = 0.0
    c%s = 's'
  else
    c%od(1) = exp_scalar
    c%od(2:8) = 0.0D0
    c%s = 'd'
  end if 
else
  cos_vnorm = cos(vnorm)
  sin_vnorm = sin(vnorm) / vnorm
  if (self%s.eq.'s') then 
    c%o(1) = real(exp_scalar * cos_vnorm)
    c%o(2:8) = real(exp_scalar * sin_vnorm * v)
    c%s = 's'
  else
    c%od(1) = exp_scalar * cos_vnorm
    c%od(2:8) = exp_scalar * sin_vnorm * v
    c%s = 'd'
  end if 
end if

end function o_exp_

!--------------------------------------------------------------------------
function o_log_(self) result(c)
!DEC$ ATTRIBUTES DLLEXPORT :: o_log_
!! author: MDG 
!! version: 1.0 
!! date: 07/16/25
!!
!!  logarithm of an octonion (carried out in double precision)

use mod_io

class(Octonion_T), INTENT(IN)   :: self
type(Octonion_T)                :: c

type(IO_T)                      :: Message 

real(kind=dbl)                  :: norm_a, vnorm, theta, scalar
real(kind=dbl)                  :: v(7)

if (self%s.eq.'s') then 
  scalar = dble(self%o(1))
  v = dble(self%o(2:8))
else
  scalar = self%od(1)
  v = self%od(2:8)
end if 

norm_a = self%o_norm_()
vnorm = sqrt(sum(v**2))

if (norm_a == 0.0D0) call Message%printError('o_log','Logarithm of zero octonion is undefined')

if (vnorm == 0.0D0) then
  if (self%s.eq.'s') then 
    c%o(1) = real(log(norm_a))
    c%o(2:8) = 0.0
    c%s = 's'
  else 
    c%od(1) = log(norm_a)
    c%od(2:8) = 0.D0
    c%s = 'd'
  end if 
else
    theta = acos(scalar / norm_a)
    if (self%s.eq.'s') then
      c%o(1) = real(log(norm_a))
      c%o(2:8) = real(theta * (v / vnorm))
      c%s = 's'
    else
      c%od(1) = log(norm_a)
      c%od(2:8) = theta * (v / vnorm)
      c%s = 'd'
    end if 
end if

end function o_log_

!--------------------------------------------------------------------------
function o_innerproduct_(self, b) result(p)
!DEC$ ATTRIBUTES DLLEXPORT :: o_innerproduct_
!! author: MDG 
!! version: 1.0 
!! date: 07/17/25
!!
!! inner product of two octonions

class(Octonion_T),INTENT(IN)    :: self
type(Octonion_T),INTENT(IN)     :: b 
real(kind=dbl)                  :: p 

real(kind=sgl)                  :: ps 

if (self%s.eq.'s') then 
  ps = sum( self%o(1:8) * b%o(1:8) )
  p = dble(ps)
else
  p = sum( self%od(1:8) * b%od(1:8) )
end if 

end function o_innerproduct_

!--------------------------------------------------------------------------
function o_arrayinnerproduct_(self, b) result(p)
!DEC$ ATTRIBUTES DLLEXPORT :: o_arrayinnerproduct_
!! author: MDG 
!! version: 1.0 
!! date: 07/17/25
!!
!! inner product of two octonions

use mod_io 

class(OctonionArray_T),INTENT(IN)    :: self
type(OctonionArray_T),INTENT(IN)     :: b 
real(kind=dbl),allocatable           :: p(:) 

type(IO_T)                           :: Message

type(Octonion_T)                     :: o1, o2
real(kind=sgl)                       :: ip 
integer(kind=irg)                    :: i

! test to make sure that both arrays have the same number of quaternions
if (self%n.ne.b%n) then
  call Message%printError('o_arrayinnerproduct_','input arrays must have the same number of quaternions')
end if

if (self%s.ne.b%s) then
  call Message%printError('o_arrayinnerproduct_','input arrays must have the same precision')
end if

if (allocated(p)) deallocate(p)
allocate(p(self%n))

do i=1,self%n 
  o1 = self%extractfromOctonionArray_(i)
  o2 = b%extractfromOctonionArray_(i)
  p(i) = o1%o_innerproduct_(o2)
end do

end function o_arrayinnerproduct_

!--------------------------------------------------------------------------
recursive function getocts_(self) result(os)
!DEC$ ATTRIBUTES DLLEXPORT :: getocts_
  !! author: MDG
  !! version: 1.0
  !! date: 07/16/25
  !!
  !! return an octonion 

IMPLICIT NONE

class(Octonion_T),INTENT(IN)    :: self

real(kind=sgl)                  :: os(8)

os = self%o

end function getocts_

!--------------------------------------------------------------------------
recursive function getoctd_(self) result(od)
!DEC$ ATTRIBUTES DLLEXPORT :: getoctd_
  !! author: MDG
  !! version: 1.0
  !! date: 07/16/25
  !!
  !! return an octonion 

IMPLICIT NONE

class(Octonion_T),INTENT(IN)    :: self

real(kind=dbl)                  :: od(8)

od = self%od

end function getoctd_

!--------------------------------------------------------------------------
recursive subroutine setocts_(self, os)
!DEC$ ATTRIBUTES DLLEXPORT :: setocts_
  !! author: MDG
  !! version: 1.0
  !! date: 07/16/25
  !!
  !! set an octonion

IMPLICIT NONE

class(Octonion_T),INTENT(INOUT)    :: self
real(kind=sgl),INTENT(IN)          :: os(8)

self%o = os
self%s = 's'

end subroutine setocts_

!--------------------------------------------------------------------------
recursive subroutine setoctd_(self, od)
!DEC$ ATTRIBUTES DLLEXPORT :: setoctd_
  !! author: MDG
  !! version: 1.0
  !! date: 07/16/25
  !!
  !! set an octonion

IMPLICIT NONE

class(Octonion_T),INTENT(INOUT)    :: self
real(kind=dbl),INTENT(IN)          :: od(8)

self%od= od
self%s = 'd'

end subroutine setoctd_

!--------------------------------------------------------------------------
recursive function extractfromOctonionArray_(self, i) result(o)
!DEC$ ATTRIBUTES DLLEXPORT :: extractfromOctonionArray_
  !! author: MDG
  !! version: 1.0
  !! date: 07/16/25
  !!
  !! extract an octonion from an array 

IMPLICIT NONE

class(OctonionArray_T),INTENT(IN)     :: self
integer(kind=irg),INTENT(IN)          :: i
type(Octonion_T)                      :: o

if (self%s.eq.'s') then 
  o = Octonion_T( o = self%o(1:8,i) )
else
  o = Octonion_T( od = self%od(1:8,i) )
end if

end function extractfromOctonionArray_

!--------------------------------------------------------------------------
recursive subroutine insertOctintoArray_(self, i, o)
!DEC$ ATTRIBUTES DLLEXPORT :: insertOctintoArray_
  !! author: MDG
  !! version: 1.0
  !! date: 07/16/25
  !!
  !! insert an octonion in an existing array 

use mod_io 

IMPLICIT NONE

class(OctonionArray_T),INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)          :: i
type(Octonion_T),INTENT(INOUT)        :: o

type(IO_T)                            :: Message

! make sure that the index i is within the appropriate range 
if (i.gt.self%n) call Message%printError('insertOctintoArray_',' index too large for octonion array')

if (self%s.eq.'s') then 
  self%o(1:8,i) = o%getocts_()
else
  self%od(1:8,i) = o%getoctd_()
end if

end subroutine insertOctintoArray_

!--------------------------------------------------------------------------
recursive subroutine octprint_(self, st, onum)
!DEC$ ATTRIBUTES DLLEXPORT :: octprint_
  !! author: MDG
  !! version: 1.0
  !! date: 07/16/25
  !!
  !! print an octonion 

use mod_io

IMPLICIT NONE

class(Octonion_T),INTENT(IN)            :: self
character(*),INTENT(IN),OPTIONAL        :: st
integer(kind=irg),INTENT(IN),OPTIONAL   :: onum

type(IO_T)                              :: Message
integer(kind=irg)                       :: io_int(1)

if (present(st)) call Message%printMessage( trim(st), frm='(A,$)' )

if (present(onum)) then 
  io_int(1) = onum 
  call Message%WriteValue('',io_int,1,frm="(I4,':',$)")
  if (self%s.eq.'s') then
    call Message % WriteValue('', self%o, 8, frm="('(',8f12.6,')')")
  else
    call Message % WriteValue('', self%od, 8, frm="('(',8f20.14,')')")
  end if
else
  if (self%s.eq.'s') then
    call Message % WriteValue('', self%o, 8, frm="('(',8f12.6,'); precision: '$)")
    call Message % WriteValue('',self%s)
  else
    call Message % WriteValue('', self%od, 8, frm="('(',8f20.14,'); precision: '$)")
    call Message % WriteValue('',self%s)
  end if
end if 

end subroutine octprint_

!--------------------------------------------------------------------------
recursive subroutine o_arrayprint_(self, listN)
!DEC$ ATTRIBUTES DLLEXPORT :: o_arrayprint_
  !! author: MDG 
  !! version: 1.0 
  !! date: 07/17/25
  !!
  !! print an array of listN octonions

IMPLICIT NONE 

class(OctonionArray_T),intent(in)     :: self
integer(kind=irg),INTENT(IN),OPTIONAL :: listN

type(Octonion_T)                      :: oct
integer(kind=irg)                     :: i, n

if (present(listN)) then 
  n = listN
else 
  n = self%n 
end if

do i=1,n
  oct = self%getOctfromArray(i)
  call oct%octprint_(onum=i)
end do

end subroutine o_arrayprint_

!--------------------------------------------------------------------------
recursive function getOnumber_(self) result(num)
!DEC$ ATTRIBUTES DLLEXPORT :: getOnumber_
  !! author: MDG
  !! version: 1.0
  !! date: 07/17/25
  !!
  !! returns the number of octonions in the QuaternionArray_T class

IMPLICIT NONE

class(OctonionArray_T), INTENT(INOUT)     :: self
integer(kind=irg)                         :: num

num = self%n

end function getOnumber_

!--------------------------------------------------------------------------
recursive subroutine deleteArray_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: deleteArray_
  !! author: MDG
  !! version: 1.0
  !! date: 07/17/25
  !!
  !! deletes the current array of quaternions in this class

IMPLICIT NONE

class(OctonionArray_T), INTENT(INOUT)   :: self

if (self%s.eq.'s') then 
  if (allocated(self%o)) deallocate(self%o)
else 
  if (allocated(self%od)) deallocate(self%od)
end if

self%n = 0

end subroutine deleteArray_

!--------------------------------------------------------------------------
recursive subroutine writeArraytoFile_(self, filename)
!DEC$ ATTRIBUTES DLLEXPORT :: writeArraytoFile_
  !! author: MDG
  !! version: 1.0
  !! date: 07/17/25
  !!
  !! write the current array to a text file (mostly used for debugging purposes)

IMPLICIT NONE

class(OctonionArray_T), INTENT(INOUT)     :: self
character(fnlen),INTENT(IN)               :: filename

integer(kind=irg)                         :: i 

open(dataunit2,file=trim(filename),status='unknown',form='formatted')
write (dataunit2,"(A)") 'oc'
write (dataunit2,"(I6)") self%n
do i=1,self%n
  if (self%s.eq.'s') then 
    write (dataunit2,"(8(F10.8,' '))") self%o(1:8,i)
  else 
    write (dataunit2,"(8(F10.8,' '))") real(self%od(1:8,i))
  end if
end do

close(dataunit2,status='keep')

end subroutine writeArraytoFile_


end module mod_octonions
