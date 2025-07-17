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

module mod_GBoctonions
  !! author: MDG 
  !! version: 1.0 
  !! date: 07/16/25
  !!
  !! class definition for the Grain Boundary Octonions module
  !!
  !! there is a stand alone Octonion_T class and an OctonionArray_T class
  !! and we inherit those from the mod_octonions module
  !!
  !! What turns an octonion into a GBoctonion is the fact that it is made up
  !! of two unit quaternions, so it requires a normalization factor of sqrt(2) 
  !!
  !! see the following publication for details:
  !!
  !! T. Francis, I. Chesser, S. Singh, E.A. Holm and M. De Graef. 
  !! "A Geodesic Octonion Metric for Grain Boundaries". 
  !! Acta Materialia, 166:135-147 (2019)
  !! DOI: https://doi.org/10.1016/j.actamat.2018.12.034

use mod_kinds
use mod_global
use mod_octonions
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit

IMPLICIT NONE 

! If we use this module, which extends the Octonion_T class, then by definition
! we are using grain boundary octonions, which are constructed from two unit quaternions.
! Therefore, they use an octonion normalization by a factor of sqrt(2), but this is 
! handled transparently by the parent Octonion_T class.
! Other than that, there really are not many differences between the two
! classes.  There are of course grain boundary specific operations which are defined in this module.

! class definition
type, public, extends(Octonion_T) :: GBoctonion_T
private 

contains
private 


end type GBoctonion_T

! ! next we define the quaternion array class
! type, public, extends(OctonionArray_T) :: GBOctonionArray_T
! private
!     integer(kind=irg)            :: n
!     integer(kind=irg)            :: nthreads
!     real(kind=sgl), allocatable  :: o(:,:)
!     real(kind=dbl), allocatable  :: od(:,:)

!   contains
!   private
! ! quaternion IO routines
!     procedure, pass(self) :: octarrayprint_
! ! quaternion arithmetic routines
!     procedure, pass(self) :: octarrayadd_
!     procedure, pass(self) :: octarraysubtract_
!     procedure, pass(self) :: octarraymult_
!     procedure, pass(self) :: octarraysmult_
!     procedure, pass(self) :: octarrayinverse_
!     procedure, pass(self) :: octarraydivide_
!     procedure, pass(self) :: octarraysdiv_
!     procedure, pass(self) :: octarrayconjg_
!     procedure, pass(self) :: octarraynorm_
!     procedure, pass(self) :: octarraynormalize_
! ! miscellaneous routines
!     procedure, pass(self) :: extractfromOctArray_
!     procedure, pass(self) :: insertOctintoArray_
!     procedure, pass(self) :: getOnumber_
!     procedure, pass(self) :: deleteArray_

! ! generics
!     generic, public :: octarray_print => octarrayprint_
!     generic, public :: operator(+) => octarrayadd_
!     generic, public :: operator(-) => octarraysubtract_
!     generic, public :: operator(*) => octarraymult_
!     generic, public :: operator(*) => octarraysmult_
!     generic, public :: operator(/) => octarraydivide_
!     generic, public :: operator(/) => octarraysdiv_
!     generic, public :: octarray_normalize => octarraynormalize_
!     generic, public :: octarray_inverse => octarrayinverse_
!     generic, public :: getOctfromArray => extractfromOctArray_
!     generic, public :: insertOctinArray => insertOctintoArray_
!     generic, public :: getOnumber => getOnumber_
!     generic, public :: deleteArray => deleteArray_

!   end type GBOctonionArray_T

! the constructor routines for these classes 
interface GBoctonion_T
  module procedure GBoctonion_constructor
end interface GBoctonion_T

! interface octonionArray_T
!   module procedure OctonionArray_constructor
! end interface octonionArray_T

contains

!--------------------------------------------------------------------------
type(GBoctonion_T) function GBoctonion_constructor( qu1, qu2 ) result(GBoctonion)
!DEC$ ATTRIBUTES DLLEXPORT :: GBoctonion_constructor
!! author: MDG 
!! version: 1.0 
!! date: 10/16/22
!!
!! constructor for the GBoctonions_T Class
 
use mod_quaternions

IMPLICIT NONE

type(Quaternion_T), INTENT(IN)    :: qu1
type(Quaternion_T), INTENT(IN)    :: qu2

if (qu1%quat_getprecision().eq.'s') then
  GBoctonion%o = (/ qu1%get_quats(), qu2%get_quats() /)
  GBoctonion%s = 's'
else 
  GBoctonion%od = (/ qu1%get_quatd(), qu2%get_quatd() /)
  GBoctonion%s = 'd'
end if 

! this normalization involves sqrt(2) due to the two unit quaternions, but this 
! is correctly handled by the parent class Octonion_T
call GBoctonion%o_normalize()

end function GBoctonion_constructor

! !--------------------------------------------------------------------------
! subroutine Octonion_destructor(self) 
! !! author: MDG 
! !! version: 1.0 
! !! date: 10/16/22
! !!
! !! destructor for the octonions_T Class
 
! IMPLICIT NONE

! type(GBoctonion_T), INTENT(INOUT)  :: self 

! call reportDestructor('GBoctonion_T')

! end subroutine octonion_destructor

! !--------------------------------------------------------------------------
! type(GBOctonionArray_T) function OctonionArray_constructor( n, nthreads, o, od, s ) result(OctArray)
! !DEC$ ATTRIBUTES DLLEXPORT :: OctonionArray_constructor
!   !! author: MDG
!   !! version: 1.0
!   !! date: 10/18/22
!   !!
!   !! constructor for the OctonionArray Class
!   !!
!   !! either call with parameters n and s
!   !! or with n and either one of o or od

! IMPLICIT NONE

!   integer(kind=irg), INTENT(IN)             :: n
!   integer(kind=irg), INTENT(IN), OPTIONAL   :: nthreads
!   real(kind=sgl), INTENT(IN), OPTIONAL      :: o(8,n)
!   real(kind=dbl), INTENT(IN), OPTIONAL      :: od(8,n)
!   character(1), INTENT(IN), OPTIONAL        :: s

! ! OpenMP threads
!   OctArray % nthreads = 0
!   if (present(nthreads)) OctArray % nthreads = nthreads

! ! are we declaring just an empty variable with no entries, but with a given precision ?
!   if ( present(s) .and. (.not.present(o)) .and. (.not.present(od)) ) then
!     OctArray % n = n
!     if (octonionprecision.eq.'s') then
!       allocate(OctArray % o(8,n))
!       OctArray % o = 0.0
!     else
!       allocate(OctArray % od(8,n))
!       OctArray % od = 0.D0
!     end if
!     return
!   end if

! ! single precision
!   if (present(o)) then
!     allocate(OctArray % o(8,n))
!     OctArray % n = n
!     OctArray % o = o
!   end if

! ! double precision
!   if (present(od)) then
!     allocate(OctArray % od(8,n))
!     OctArray % n = n
!     OctArray % od = od
!   end if

! end function OctonionArray_constructor

! !--------------------------------------------------------------------------
! subroutine OctonionArray_destructor(self)
! !DEC$ ATTRIBUTES DLLEXPORT :: OctonionArray_destructor
! !! author: MDG
! !! version: 1.0
! !! date: 10/18/22
! !!
! !! destructor for the GBOctonionArray_T Class

! IMPLICIT NONE

! type(GBOctonionArray_T), INTENT(INOUT)     :: self

! call reportDestructor('GBOctonionArray_T')

! if (allocated(self%o)) deallocate(self%o)
! if (allocated(self%od)) deallocate(self%od)

! end subroutine OctonionArray_destructor

! !--------------------------------------------------------------------------
! recursive subroutine octarraynormalize_(self)
! !DEC$ ATTRIBUTES DLLEXPORT :: octarraynormalize_
!   !! author: MDG
!   !! version: 1.0
!   !! date: 10/19/22
!   !!
!   !! normalize the input octonions

! IMPLICIT NONE

!   class(GBOctonionArray_T),intent(inout)   :: self

!   integer(kind=irg)                      :: i
!   type(Octonion_T)                       :: o 

! do i=1,self%n 
!   o = self%extractfromOctArray_(i)
!   call o%octnormalize_()
!   call self%insertOctintoArray_(i, o)
! end do   

! end subroutine octarraynormalize_


! !--------------------------------------------------------------------------!
! recursive function extractfromOctArray_(self, i) result (res)
! !DEC$ ATTRIBUTES DLLEXPORT :: extractfromOctArray_
!   !! author: MDG
!   !! version: 1.0
!   !! date: 10/18/22
!   !!
!   !! extract an octonion from an array of octonions

! use mod_io

! IMPLICIT NONE

!   class(GBOctonionArray_T),intent(in)   :: self
!   integer(kind=irg), intent(in)       :: i
!   type(Octonion_T)                    :: res

!   type(IO_T)                          :: Message

!   if (i.le.self%n) then
!     if (octonionprecision.eq.'s') then 
!       res = Octonion_T( o = self%o(1:8,i) )
!     else
!       res = Octonion_T( od = self%od(1:8,i) )
!     end if 
!   else
!     call Message%printWarning('extractfromOctonionArray_: requested octonion index larger than array size', &
!                               (/'   ---> returning empty octonion'/) )
!     if (octonionprecision.eq.'s') then
!       res = Octonion_T( smode='s' )
!     else
!       res = Octonion_T( )
!     end if
!   end if

! end function extractfromOctArray_

! !--------------------------------------------------------------------------!
! recursive subroutine insertOctintoArray_(self, i, o)
! !DEC$ ATTRIBUTES DLLEXPORT :: insertOctintoArray_
!   !! author: MDG
!   !! version: 1.0
!   !! date: 01/23/20
!   !!
!   !! insert an octonion into an array of octonions

! use mod_io

! IMPLICIT NONE

!   class(GBOctonionArray_T),intent(inout):: self
!   integer(kind=irg), intent(in)       :: i
!   type(Octonion_T), intent(in)        :: o

!   type(IO_T)                            :: Message

!   if (i.le.self%n) then
!     if (octonionprecision.eq.'s') then 
!       self%o(1:8,i) = o%get_octs()
!     else
!       self%od(1:8,i) = o%get_octd()
!     end if
!   else
!     call Message%printWarning('insertOctintoArray: requested octonion index larger than array size', &
!                               (/'   ---> no octonion inserted'/) )
!   end if

! end subroutine insertOctintoArray_

! !--------------------------------------------------------------------------
! recursive subroutine deleteArray_(self)
! !DEC$ ATTRIBUTES DLLEXPORT :: deleteArray_
!   !! author: MDG
!   !! version: 1.0
!   !! date: 10/18/22
!   !!
!   !! deletes the current array of octonions in this class

! IMPLICIT NONE

! class(GBOctonionArray_T), INTENT(INOUT)   :: self

! if (octonionprecision.eq.'s') then 
!   if (allocated(self%o)) deallocate(self%o)
! else 
!   if (allocated(self%od)) deallocate(self%od)
! end if

! self%n = 0

! end subroutine deleteArray_

! !--------------------------------------------------------------------------
! recursive function getOnumber_(self) result(num)
! !DEC$ ATTRIBUTES DLLEXPORT :: getOnumber_
!   !! author: MDG
!   !! version: 1.0
!   !! date: 10/18/22
!   !!
!   !! returns the number of octonions in the GBOctonionArray_T class

! IMPLICIT NONE

! class(GBOctonionArray_T), INTENT(INOUT)   :: self
! integer(kind=irg)                       :: num

! num = self%n

! end function getOnumber_



end module mod_GBoctonions