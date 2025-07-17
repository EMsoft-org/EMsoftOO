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
use mod_quaternions

IMPLICIT NONE 

! If we use this module, which extends the Octonion_T class, then by definition
! we are using grain boundary octonions, which are constructed from two unit quaternions.
! Therefore, they use an octonion normalization by a factor of sqrt(2), but this is 
! handled transparently by the parent Octonion_T class.
! Other than that, there really are not many differences between the two
! classes.  There are of course grain boundary specific operations which are defined in this module.

! class definition for the Grain Boundary Octonion
type, public, extends(Octonion_T) :: GBoctonion_T
private 

contains
private 

end type GBoctonion_T

! class definition for the Grain Boundary Octonion Array
type, public, extends(OctonionArray_T) :: GBOctonionArray_T
private

contains
private

   procedure, pass(self) :: insertGBOctintoArray_
   generic, public :: insertGBOctinArray => insertGBOctintoArray_

end type GBOctonionArray_T

private:: insertGBOctintoArray_

! the constructor routines for these classes 
interface GBoctonion_T
  module procedure GBoctonion_constructor
end interface GBoctonion_T

interface GBoctonionArray_T
  module procedure GBOctonionArray_constructor
end interface GBoctonionArray_T

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

!--------------------------------------------------------------------------
subroutine GBOctonion_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 07/17/25
!!
!! destructor for the GBoctonion_T Class
 
IMPLICIT NONE

type(GBoctonion_T), INTENT(INOUT)  :: self 

call reportDestructor('GBoctonion_T')

end subroutine GBOctonion_destructor

!--------------------------------------------------------------------------
type(GBOctonionArray_T) function GBOctonionArray_constructor( qAr1, qAr2 ) result(OctArray)
!DEC$ ATTRIBUTES DLLEXPORT :: GBOctonionArray_constructor
  !! author: MDG
  !! version: 1.0
  !! date: 07/17/25
  !!
  !! constructor for the GBOctonionArray Class
  !!
  !! this constructor takes two QuaternionArrays and merges them into a GBOctonionArray
  !! 

use mod_io 

IMPLICIT NONE

  type(QuaternionArray_T),INTENT(INOUT)     :: qAr1 
  type(QuaternionArray_T),INTENT(INOUT)     :: qAr2 

  type(IO_T)                                :: Message 

  integer(kind=irg)                         :: i 
  type(Quaternion_T)                        :: q1, q2
  type(GBoctonion_T)                        :: gboct

! make sure the arrays have the same size
  if (qAr1%getQnumber().ne.qAr2%getQnumber()) then 
    call Message%printError('GBOctonionArray_constructor',' input quaternion arrays have different size')
  end if 

! inherit quaternion array parameters
  OctArray%nthreads = qAr1%getnthreads()
  OctArray%n = qAr1%getQnumber()
  OctArray%s = qAr1%getprecision()
  
! allocate the GBO array
  if (OctArray%s.eq.'s') then
    if (allocated(OctArray%o)) deallocate(OctArray%o)
    allocate( OctArray%o(8,OctArray%n) ) 
  else
    if (allocated(OctArray%od)) deallocate(OctArray%od)
    allocate( OctArray%od(8,OctArray%n) ) 
  end if 

  do i=1,OctArray%n
    q1 = qAr1%getQuatfromArray(i)
    q2 = qAr2%getQuatfromArray(i)
    gboct = GBoctonion_T( q1, q2 )
    call OctArray%insertGBOctinArray(i, gboct)
  end do

end function GBOctonionArray_constructor

!--------------------------------------------------------------------------
subroutine GBOctonionArray_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: GBOctonionArray_destructor
!! author: MDG
!! version: 1.0
!! date: 07/17/25
!!
!! destructor for the GBOctonionArray_T Class

IMPLICIT NONE

type(GBOctonionArray_T), INTENT(INOUT)     :: self

call reportDestructor('GBOctonionArray_T')

if (allocated(self%o)) deallocate(self%o)
if (allocated(self%od)) deallocate(self%od)

end subroutine GBOctonionArray_destructor

!--------------------------------------------------------------------------
recursive subroutine insertGBOctintoArray_(self, i, o)
!DEC$ ATTRIBUTES DLLEXPORT :: insertGBOctintoArray_
  !! author: MDG
  !! version: 1.0
  !! date: 07/16/25
  !!
  !! insert a GBoctonion in an existing array (overrides mod_octonions)

use mod_io 

IMPLICIT NONE

class(GBOctonionArray_T),INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)            :: i
type(GBOctonion_T),INTENT(INOUT)        :: o

type(IO_T)                              :: Message

! make sure that the index i is within the appropriate range 
if (i.gt.self%n) call Message%printError('insertGBOctintoArray_',' index too large for octonion array')

if (self%s.eq.'s') then 
  self%o(1:8,i) = o%get_octs()
else
  self%od(1:8,i) = o%get_octd()
end if

end subroutine insertGBOctintoArray_


end module mod_GBoctonions