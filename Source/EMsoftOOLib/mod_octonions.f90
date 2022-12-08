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

module mod_octonions
  !! author: MDG 
  !! version: 1.0 
  !! date: 10/16/22
  !!
  !! class definition for the octonions module
  !!
  !! octonion multiplication is carried out in terms of the member quaternions
  !! both member quaternions are either single or double precision
  !!
  !! there is a stand alone octonion class and an octonion array class
  !!
  !! there are two ways of octonion normalization, the normal way, and the 
  !! grain boundary way (combination of two unit quaternions and sqrt(2))
  !!
  !! see the following publication for details:
  !!
  !! T. Francis, I. Chesser, S. Singh, E.A. Holm and M. De Graef. 
  !! "A Geodesic Octonion Metric for Grain Boundaries". 
  !! Acta Materialia, 166:135-147 (2019)
  !! DOI: https://doi.org/10.1016/j.actamat.2018.12.034

use mod_kinds
use mod_global
use mod_quaternions, only:Quaternion_T
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit

IMPLICIT NONE 

! This parameter determines which normalization procedure is followed for the octonions
! either regular normalization (divide by square root of sum of squares) or
! the Grain Boundary Octonion normalization in which each quaternion is normalized 
! and then the whole is divided by sqrt(2).  Since most octonion computations in EMsoft
! are likely about grain boundaries, we set this parameter to .TRUE. as the default.
logical, private :: octonionGBmode = .TRUE.

! set the precision for all computations to 'd' for double (default) or 's' for single precision
character(1), private :: octonionprecision= 'd'

! overload the conjg and cabs routines
intrinsic :: conjg, cabs
public :: conjg, cabs, get_octonionGBmode, set_octonionGBmode, set_octonionprecision, get_octonionprecision

interface conjg
  procedure octconjg_
  procedure octarrayconjg_
end interface conjg

interface cabs
  procedure octnorm_
  procedure octarraynorm_
end interface cabs

! class definition
type, public :: octonion_T
private 
  type(Quaternion_T)          :: o(2)

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
  procedure, pass(self) :: getocts_
  procedure, pass(self) :: getoctd_
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
  generic, public :: get_octs => getocts_
  generic, public :: get_octd => getoctd_
  generic, public :: set_octs => setocts_
  generic, public :: set_octd => setoctd_

end type octonion_T

! next we define the quaternion array class
type, public :: OctonionArray_T
  !! OctonionArray_T Class definition
  private
    integer(kind=irg)            :: n
    integer(kind=irg)            :: nthreads
    real(kind=sgl), allocatable  :: o(:,:)
    real(kind=dbl), allocatable  :: od(:,:)

  contains
  private
! quaternion IO routines
    procedure, pass(self) :: octarrayprint_
! quaternion arithmetic routines
    procedure, pass(self) :: octarrayadd_
    procedure, pass(self) :: octarraysubtract_
    procedure, pass(self) :: octarraymult_
    procedure, pass(self) :: octarraysmult_
    procedure, pass(self) :: octarrayinverse_
    procedure, pass(self) :: octarraydivide_
    procedure, pass(self) :: octarraysdiv_
    procedure, pass(self) :: octarrayconjg_
    procedure, pass(self) :: octarraynorm_
    procedure, pass(self) :: octarraynormalize_
! miscellaneous routines
    procedure, pass(self) :: extractfromOctArray_
    procedure, pass(self) :: insertOctintoArray_
    procedure, pass(self) :: getOnumber_
    procedure, pass(self) :: deleteArray_

! generics
    generic, public :: octarray_print => octarrayprint_
    generic, public :: operator(+) => octarrayadd_
    generic, public :: operator(-) => octarraysubtract_
    generic, public :: operator(*) => octarraymult_
    generic, public :: operator(*) => octarraysmult_
    generic, public :: operator(/) => octarraydivide_
    generic, public :: operator(/) => octarraysdiv_
    generic, public :: octarray_normalize => octarraynormalize_
    generic, public :: octarray_inverse => octarrayinverse_
    generic, public :: getOctfromArray => extractfromOctArray_
    generic, public :: insertOctinArray => insertOctintoArray_
    generic, public :: getOnumber => getOnumber_
    generic, public :: deleteArray => deleteArray_

  end type OctonionArray_T

! the constructor routines for these classes 
interface octonion_T
  module procedure Octonion_constructor
end interface octonion_T

interface octonionArray_T
  module procedure OctonionArray_constructor
end interface octonionArray_T

contains

!--------------------------------------------------------------------------
type(octonion_T) function Octonion_constructor( o, od, smode ) result(octonion)
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
type(OctonionArray_T) function OctonionArray_constructor( n, nthreads, o, od, s ) result(OctArray)
!DEC$ ATTRIBUTES DLLEXPORT :: OctonionArray_constructor
  !! author: MDG
  !! version: 1.0
  !! date: 10/18/22
  !!
  !! constructor for the OctonionArray Class
  !!
  !! either call with parameters n and s
  !! or with n and either one of o or od

IMPLICIT NONE

  integer(kind=irg), INTENT(IN)             :: n
  integer(kind=irg), INTENT(IN), OPTIONAL   :: nthreads
  real(kind=sgl), INTENT(IN), OPTIONAL      :: o(8,n)
  real(kind=dbl), INTENT(IN), OPTIONAL      :: od(8,n)
  character(1), INTENT(IN), OPTIONAL        :: s

! OpenMP threads
  OctArray % nthreads = 0
  if (present(nthreads)) OctArray % nthreads = nthreads

! are we declaring just an empty variable with no entries, but with a given precision ?
  if ( present(s) .and. (.not.present(o)) .and. (.not.present(od)) ) then
    OctArray % n = n
    if (octonionprecision.eq.'s') then
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
  end if

! double precision
  if (present(od)) then
    allocate(OctArray % od(8,n))
    OctArray % n = n
    OctArray % od = od
  end if

end function OctonionArray_constructor

!--------------------------------------------------------------------------
subroutine OctonionArray_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: OctonionArray_destructor
!! author: MDG
!! version: 1.0
!! date: 10/18/22
!!
!! destructor for the OctonionArray_T Class

IMPLICIT NONE

type(OctonionArray_T), INTENT(INOUT)     :: self

call reportDestructor('OctonionArray_T')

if (allocated(self%o)) deallocate(self%o)
if (allocated(self%od)) deallocate(self%od)

end subroutine OctonionArray_destructor

!--------------------------------------------------------------------------
function get_octonionGBmode() result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_octonionGBmode
!! author: MDG 
!! version: 1.0 
!! date: 10/19/22
!!
!! get the global octonionGBmode parameter

IMPLICIT NONE 

logical                              :: out

out = octonionGBmode

end function get_octonionGBmode

!--------------------------------------------------------------------------
subroutine set_octonionGBmode(inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_octonionGBmode
!! author: MDG 
!! version: 1.0 
!! date: 10/19/22
!!
!! set the octonionGBmode parameter

IMPLICIT NONE 

logical, INTENT(IN)                  :: inp

octonionGBmode = inp

end subroutine set_octonionGBmode

!--------------------------------------------------------------------------
function get_octonionprecision() result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_octonionprecision
!! author: MDG 
!! version: 1.0 
!! date: 10/19/22
!!
!! get the global octonionGBmode parameter

IMPLICIT NONE 

character(1)                    :: out

out = octonionprecision

end function get_octonionprecision

!--------------------------------------------------------------------------
subroutine set_octonionprecision(inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_octonionprecision
!! author: MDG 
!! version: 1.0 
!! date: 10/19/22
!!
!! set the octonionGBmode parameter

IMPLICIT NONE 

character(1), INTENT(IN)        :: inp

octonionprecision = inp

end subroutine set_octonionprecision

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

if (octonionprecision.eq.'s') then
  call Message % WriteValue(trim(m), (/ self%o(1)%get_quats(), self%o(2)%get_quats() /), 8, frm="('(',8f16.6,'); precision: '$)")
  call Message % WriteValue('',octonionprecision)
else
  call Message % WriteValue(trim(m), (/ self%o(1)%get_quatd(), self%o(2)%get_quatd() /), 8, frm="('(',8f24.14,'); precision: '$)")
  call Message % WriteValue('',octonionprecision)
end if

end subroutine octprint_

!--------------------------------------------------------------------------
recursive subroutine octarrayprint_(self, listN)
!DEC$ ATTRIBUTES DLLEXPORT :: octarrayprint_
  !! author: MDG 
  !! version: 1.0 
  !! date: 10/18/22
  !!
  !! print an array of octonions

use mod_io

IMPLICIT NONE 

  class(OctonionArray_T),intent(in)     :: self
  integer(kind=irg),INTENT(IN),OPTIONAL :: listN

  type(IO_T)                            :: Message 
  integer(kind=irg)                     :: i, n

  if (present(listN)) then 
    n = listN
  else 
    n = self%n 
  end if
  if (octonionprecision.eq.'s') then 
    do i=1,n
      call Message % WriteValue('', self%o(:,i), 8, frm="('(',8f12.6,')')")
    end do
  else 
    do i=1,n
      call Message % WriteValue('', self%od(:,i), 8, frm="('(',8f20.14,')')")
    end do
  end if 

end subroutine octarrayprint_

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

  if (octonionprecision.eq.'s') then
    d = maxval( abs( diff%getocts_() ) )
    if (d.gt.eps) res = .FALSE.
  else
    dd = maxval( abs( diff%getoctd_() ) )
    if (dd.gt.epsd) res = .FALSE.
  end if

end function octsequal_

!--------------------------------------------------------------------------
recursive function octconjg_(self) result(ores)
!DEC$ ATTRIBUTES DLLEXPORT :: octconjg_
  !! author: MDG
  !! version: 1.0
  !! date: 10/16/22
  !!
  !! octonion conjugation (single/double precision)

IMPLICIT NONE

class(Octonion_T),intent(in) :: self
type(Octonion_T)             :: ores

real(kind=sgl)               :: a(8)
real(kind=dbl)               :: b(8)

if (octonionprecision.eq.'s') then
  a = self%get_octs()
  a(2:8) = -a(2:8)
  call ores%set_octs( a )
else
  b = self%get_octd()
  b(2:8) = -b(2:8)
  call ores%set_octd( b )
end if

end function octconjg_

!--------------------------------------------------------------------------
recursive function octarrayconjg_(self) result(ores)
!DEC$ ATTRIBUTES DLLEXPORT :: octarrayconjg_
  !! author: MDG
  !! version: 1.0
  !! date: 10/16/22
  !!
  !! octonion conjugation (single/double precision)

IMPLICIT NONE

class(OctonionArray_T),intent(in) :: self
type(OctonionArray_T)             :: ores

integer(kind=irg)                 :: i 
type(Octonion_T)                  :: o 

ores = self 
do i=1,self%n 
  o = self%extractfromOctArray_(i)
  call ores%insertOctintoArray_(i, conjg(o))
end do

end function octarrayconjg_

!--------------------------------------------------------------------------
recursive function octadd_(self, y) result(ores)
!DEC$ ATTRIBUTES DLLEXPORT :: octadd_
  !! author: MDG
  !! version: 1.0
  !! date: 10/16/22
  !!
  !! octonion addition (single/double precision)

IMPLICIT NONE

class(Octonion_T),intent(in) :: self, y
type(Octonion_T)             :: ores

if (octonionprecision.eq.'s') then
  call ores%set_octs( self%get_octs() + y%get_octs() )
else
  call ores%set_octd( self%get_octd() + y%get_octd() )
end if

end function octadd_

!--------------------------------------------------------------------------
recursive function octarrayadd_(self, y) result(ores)
!DEC$ ATTRIBUTES DLLEXPORT :: octarrayadd_
  !! author: MDG
  !! version: 1.0
  !! date: 10/18/22
  !!
  !! octonion array addition (single/double precision)

use mod_io

IMPLICIT NONE

  class(OctonionArray_T),intent(in) :: self, y
  type(OctonionArray_T)             :: ores

  type(IO_T)                        :: Message
  integer(kind=irg)                 :: sz(2)

! test to make sure that both arrays have the same number of quaternions
  if (self%n.ne.y%n) then
    call Message%printError('octarrayadd','input arrays must have the same number of octonions')
  end if

  ores%n = self%n
  ores%nthreads = self%nthreads

  if (octonionprecision.eq.'s') then
! if the octonion array is already allocated, check to make sure it has the right dimensions
    if (allocated(ores%o)) then
      sz = shape(ores%o)
      if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(ores%o)
    end if
    allocate(ores%o(8,self%n))
    ores%o = self%o + y%o
  else
! if the octonion array is already allocated, check to make sure it has the right dimensions
    if (allocated(ores%od)) then
      sz = shape(ores%od)
      if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(ores%od)
    end if
    allocate(ores%od(4,self%n))
    ores%od = self%od + y%od
  end if

end function octarrayadd_

!--------------------------------------------------------------------------
recursive function octsubtract_(self, y) result(ores)
!DEC$ ATTRIBUTES DLLEXPORT :: octsubtract_
  !! author: MDG
  !! version: 1.0
  !! date: 10/16/22
  !!
  !! octonion subtraction (single/double precision)

IMPLICIT NONE

class(Octonion_T),intent(in) :: self, y
type(Octonion_T)             :: ores

if (octonionprecision.eq.'s') then
  call ores%set_octs( self%get_octs() - y%get_octs() )
else
  call ores%set_octd( self%get_octd() - y%get_octd() )
end if

end function octsubtract_

!--------------------------------------------------------------------------
recursive function octarraysubtract_(self, y) result(ores)
!DEC$ ATTRIBUTES DLLEXPORT :: octarraysubtract_
  !! author: MDG
  !! version: 1.0
  !! date: 10/18/22
  !!
  !! octonion array subtraction (single/double precision)

use mod_io

IMPLICIT NONE

  class(OctonionArray_T),intent(in) :: self, y
  type(OctonionArray_T)             :: ores

  type(IO_T)                        :: Message
  integer(kind=irg)                 :: sz(2)

! test to make sure that both arrays have the same number of octonions
  if (self%n.ne.y%n) then
    call Message%printError('octarraysubtract','input arrays must have the same number of octonions')
  end if

  ores%n = self%n
  ores%nthreads = self%nthreads

   if (octonionprecision.eq.'s') then
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(ores%o)) then
      sz = shape(ores%o)
      if ( (sz(1).ne.8).or.(sz(2).ne.self%n) ) deallocate(ores%o)
    end if
    allocate(ores%o(8,self%n))
    ores%o = self%o - y%o
  else
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(ores%od)) then
      sz = shape(ores%od)
      if ( (sz(1).ne.8).or.(sz(2).ne.self%n) ) deallocate(ores%od)
    end if
    allocate(ores%od(8,self%n))
    ores%od = self%od - y%od
  end if

end function octarraysubtract_

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
if (octonionprecision.eq.'s') then
  call qres%set_octs( (/ a%get_quats(), b%get_quats() /) )
else
  call qres%set_octd( (/ a%get_quatd(), b%get_quatd() /) )
end if

end function octmult_

!--------------------------------------------------------------------------
recursive function octarraymult_(self, y) result(ores)
!DEC$ ATTRIBUTES DLLEXPORT :: octarraymult_
  !! author: MDG
  !! version: 1.0
  !! date: 10/18/22
  !!
  !! octonion array multiplication   (single/double precision)

use mod_io
use omp_lib
use mod_OMPsupport

IMPLICIT NONE

  class(OctonionArray_T),intent(in) :: self, y
  type(OctonionArray_T)             :: ores

  type(Octonion_T)                  :: oct
  type(IO_T)                        :: Message
  integer(kind=irg)                 :: i, sz(2)

! test to make sure that both arrays have the same number of octonions
  if (self%n.ne.y%n) then
    call Message%printError('octarraymult','input arrays must have the same number of octonions')
  end if

  ores%n = self%n
  ores%nthreads = self%nthreads

! set the number of OpenMP threads
  call OMP_setNThreads(ores%nthreads)

  if (octonionprecision.eq.'s') then
! if the octonion array is already allocated, check to make sure it has the right dimensions
    if (allocated(ores%o)) then
      sz = shape(ores%o)
      if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(ores%o)
    end if
    allocate(ores%o(8,self%n))

!$OMP PARALLEL DEFAULT(shared) PRIVATE(oct)
!$OMP DO SCHEDULE(DYNAMIC,1)
    do i=1,self%n
      oct = self%extractfromOctArray_(i) * y%extractfromOctArray_(i)
      call ores%insertOctintoArray_(i, oct)
    end do
!$OMP END DO
!$OMP END PARALLEL
  else
! if the octonion array is already allocated, check to make sure it has the right dimensions
    if (allocated(ores%od)) then
      sz = shape(ores%od)
      if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(ores%od)
    end if
    allocate(ores%od(8,self%n))

!$OMP PARALLEL DEFAULT(shared) PRIVATE(oct)
!$OMP DO SCHEDULE(DYNAMIC,1)
    do i=1,self%n
      oct = self%extractfromOctArray_(i) * y%extractfromOctArray_(i)
      call ores%insertOctintoArray_(i, oct)
    end do
!$OMP END DO
!$OMP END PARALLEL
  end if

end function octarraymult_

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
if (octonionprecision.eq.'s') then
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
recursive function octarraysmult_(self, s) result(ores)
!DEC$ ATTRIBUTES DLLEXPORT :: octarraysmult_
  !! author: MDG
  !! version: 1.0
  !! date: 01/08/20
  !!
  !! scalar octonion array multiplication   (single/double precision)

IMPLICIT NONE

  class(OctonionArray_T),intent(in)   :: self
  real(kind=dbl), INTENT(IN)          :: s
  type(OctonionArray_T)               :: ores

  integer(kind=irg)                   :: sz(2)

! if the octonion array is already allocated, check to make sure it has the right dimensions
  if (octonionprecision.eq.'s') then
    if (allocated(ores%o)) then
      sz = shape(ores%o)
      if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(ores%o)
    end if
    allocate(ores%o(8,self%n))

    ores%o = self%o * s
  else
    if (allocated(ores%od)) then
      sz = shape(ores%od)
      if ((sz(1).ne.8).or.(sz(2).ne.self%n)) deallocate(ores%od)
    end if
    allocate(ores%od(8,self%n))

    ores%od = self%od * s
  end if 

  ores%n = self%n
  ores%nthreads = self%nthreads

end function octarraysmult_

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
if (octonionprecision.eq.'s') then
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

if (octonionprecision.eq.'s') then
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
recursive function octarraynorm_(self) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: octarraynorm_
  !! author: MDG
  !! version: 1.0
  !! date: 10/18/22
  !!
  !! octonion array norm (extends intrinsic routine abs)
  !!
  !! this routine requires the output array to be allocated in the calling program.

IMPLICIT NONE

  class(OctonionArray_T),intent(in)   :: self
  real(kind=dbl)                      :: res(self%n)

  real(kind=dbl),allocatable          :: nd(:)
  type(Octonion_T)                    :: a
  integer(kind=irg)                   :: i

  allocate(nd(self%n))
  do i=1,self%n 
    a = self%extractfromOctArray_(i)
    nd(i) = a%octnorm_()
  end do 

  if (octonionprecision.eq.'s') then
    res = dble(sngl(nd))
  else 
    res = nd
  end if
  deallocate(nd)

end function octarraynorm_

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

if (octonionGBmode.eqv..TRUE.) then  ! we have a Grain Boundary octionion with two unit quaternions
  if (octonionprecision.eq.'s') then
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
  if (octonionprecision.eq.'s') then
    a = self%getocts_() / b
    call self%setocts_( real(a) )
  else 
    a = self%getoctd_() / b
    call self%setoctd_( a )
  end if 
end if

end subroutine octnormalize_

!--------------------------------------------------------------------------
recursive subroutine octarraynormalize_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: octarraynormalize_
  !! author: MDG
  !! version: 1.0
  !! date: 10/19/22
  !!
  !! normalize the input octonions

IMPLICIT NONE

  class(OctonionArray_T),intent(inout)   :: self

  integer(kind=irg)                      :: i
  type(Octonion_T)                       :: o 

do i=1,self%n 
  o = self%extractfromOctArray_(i)
  call o%octnormalize_()
  call self%insertOctintoArray_(i, o)
end do   

end subroutine octarraynormalize_

!--------------------------------------------------------------------------
recursive function octinverse_(self) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: octinverse_ 
  !! author: MDG
  !! version: 1.0
  !! date: 10/18/22
  !!
  !! octonion inverse

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
recursive function octarrayinverse_(self) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: octarrayinverse_ 
  !! author: MDG
  !! version: 1.0
  !! date: 10/19/22
  !!
  !! octonion array inverse

IMPLICIT NONE

  class(OctonionArray_T),intent(in)   :: self
  type(OctonionArray_T)               :: res

  type(Octonion_T)                    :: o, oinv
  integer(kind=irg)                   :: i 

res = self 
do i=1,self%n 
  o = self%extractfromOctArray_(i)
  oinv = o%octinverse_()
  call res%insertOctintoArray_(i, oinv)
end do

end function octarrayinverse_

!--------------------------------------------------------------------------
recursive function octdivide_(self,y) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: octdivide_ 
  !! author: MDG
  !! version: 1.0
  !! date: 10/18/22
  !!
  !! octonion division

IMPLICIT NONE

  class(Octonion_T),intent(in)   :: self, y
  type(Octonion_T)               :: res

res = self * y%octinverse_()

end function octdivide_

!--------------------------------------------------------------------------
recursive function octarraydivide_(self, y) result (ores)
!DEC$ ATTRIBUTES DLLEXPORT :: octarraydivide_
  !! author: MDG
  !! version: 1.0
  !! date: 10/19/22
  !!
  !! octonion array division (single/double precision)

use mod_io 

IMPLICIT NONE

  class(OctonionArray_T),intent(in)    :: self, y
  type(OctonionArray_T)                :: ores

  type(OctonionArray_T)                :: yinv
  type(IO_T)                           :: Message 

! test to make sure that both arrays have the same number of octonions
if (self%n.ne.y%n) then
  call Message%printError('octarraydivide_','input arrays must have the same number of octonions')
end if

yinv = y%octarrayinverse_() 
ores = self * yinv 

end function octarraydivide_

!--------------------------------------------------------------------------
recursive function octarraysdiv_(self, s) result(ores)
!DEC$ ATTRIBUTES DLLEXPORT :: octarraysdiv_
  !! author: MDG
  !! version: 1.0
  !! date: 10/19/22
  !!
  !! octonion array division by a scalar (single/double precision)

IMPLICIT NONE

  class(OctonionArray_T),intent(in)    :: self
  real(kind=dbl),intent(in)            :: s
  type(OctonionArray_T)                :: ores

  integer(kind=irg)                    :: i
  type(Octonion_T)                     :: o

ores = self
do i=1,self%n
  o = self%extractfromOctArray_(i)
  o = o / s 
  call ores%insertOctintoArray_(i,o)
end do 

end function octarraysdiv_

!--------------------------------------------------------------------------!
recursive function extractfromOctArray_(self, i) result (res)
!DEC$ ATTRIBUTES DLLEXPORT :: extractfromOctArray_
  !! author: MDG
  !! version: 1.0
  !! date: 10/18/22
  !!
  !! extract an octonion from an array of octonions

use mod_io

IMPLICIT NONE

  class(OctonionArray_T),intent(in)   :: self
  integer(kind=irg), intent(in)       :: i
  type(Octonion_T)                    :: res

  type(IO_T)                          :: Message

  if (i.le.self%n) then
    if (octonionprecision.eq.'s') then 
      res = Octonion_T( o = self%o(1:8,i) )
    else
      res = Octonion_T( od = self%od(1:8,i) )
    end if 
  else
    call Message%printWarning('extractfromOctonionArray_: requested octonion index larger than array size', &
                              (/'   ---> returning empty octonion'/) )
    if (octonionprecision.eq.'s') then
      res = Octonion_T( smode='s' )
    else
      res = Octonion_T( )
    end if
  end if

end function extractfromOctArray_

!--------------------------------------------------------------------------!
recursive subroutine insertOctintoArray_(self, i, o)
!DEC$ ATTRIBUTES DLLEXPORT :: insertOctintoArray_
  !! author: MDG
  !! version: 1.0
  !! date: 01/23/20
  !!
  !! insert an octonion into an array of octonions

use mod_io

IMPLICIT NONE

  class(OctonionArray_T),intent(inout):: self
  integer(kind=irg), intent(in)       :: i
  type(Octonion_T), intent(in)        :: o

  type(IO_T)                            :: Message

  if (i.le.self%n) then
    if (octonionprecision.eq.'s') then 
      self%o(1:8,i) = o%get_octs()
    else
      self%od(1:8,i) = o%get_octd()
    end if
  else
    call Message%printWarning('insertOctintoArray: requested octonion index larger than array size', &
                              (/'   ---> no octonion inserted'/) )
  end if

end subroutine insertOctintoArray_

!--------------------------------------------------------------------------
recursive subroutine deleteArray_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: deleteArray_
  !! author: MDG
  !! version: 1.0
  !! date: 10/18/22
  !!
  !! deletes the current array of octonions in this class

IMPLICIT NONE

class(OctonionArray_T), INTENT(INOUT)   :: self

if (octonionprecision.eq.'s') then 
  if (allocated(self%o)) deallocate(self%o)
else 
  if (allocated(self%od)) deallocate(self%od)
end if

self%n = 0

end subroutine deleteArray_

!--------------------------------------------------------------------------
recursive function getOnumber_(self) result(num)
!DEC$ ATTRIBUTES DLLEXPORT :: getOnumber_
  !! author: MDG
  !! version: 1.0
  !! date: 10/18/22
  !!
  !! returns the number of octonions in the OctonionArray_T class

IMPLICIT NONE

class(OctonionArray_T), INTENT(INOUT)   :: self
integer(kind=irg)                       :: num

num = self%n

end function getOnumber_



end module mod_octonions