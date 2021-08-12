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

module mod_STL
  !! author: MDG 
  !! version: 1.0 
  !! date: 08/12/21
  !!
  !! class definition for the STL file format; uses the Marching Cubes module mod_MCA for the input format

use mod_kinds
use mod_global
use,intrinsic :: iso_c_binding
use mod_MCA 

IMPLICIT NONE 

private 

type STLtriangle
  real(c_float)           :: nv(3)
  real(c_float)           :: v1(3)
  real(c_float)           :: v2(3)
  real(c_float)           :: v3(3)
  integer(c_int16_t)      :: attbytecnt = 0
end type STLtriangle

! class definition
type, public :: STL_T
private 
  type(STLtriangle),allocatable   :: STLtriangles 
  character(fnlen)                :: STLfilename
  integer(kind=irg)               :: Ntriangles
  integer(kind=irg)               :: STLunit = 128
  integer(kind=irg)               :: RecLength = 50  ! bytes

contains
private 
  procedure, pass(self) :: STLtest_         ! creates a very simple test file
  procedure, pass(self) :: writeSTLfile_    ! write the actual file 

  generic, public :: STLtest => STLtest_
  generic, public :: writeSTLfile => writeSTLfile_

end type STL_T

! the constructor routine for this class 
interface STL_T
  module procedure STL_constructor
end interface STL_T

contains

!--------------------------------------------------------------------------
type(STL_T) function STL_constructor( STLfilename, STLheader, Ntriangles, doSTLtest, MCAlist ) result(STL)
!DEC$ ATTRIBUTES DLLEXPORT :: STL_constructor
!! author: MDG 
!! version: 1.0 
!! date: 08/12/21
!!
!! constructor for the STL_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen),INTENT(IN)                   :: STLfilename 
character(kind=c_char,len=80),INTENT(IN)      :: STLheader
integer(kind=irg),INTENT(IN)                  :: Ntriangles
logical,INTENT(IN),OPTIONAL                   :: doSTLtest
type(MCAtriangle),INTENT(IN),OPTIONAL,pointer :: MCAlist

integer(c_int32_t)                            :: nt 

if (present(doSTLtest)) then 
  if (doSTLtest.eqv..TRUE.) then 
    STL%STLfilename = 'octahedron.stl'

  ! open the file; this must be a little-endian file !!!
    open(unit = STL%STLunit, file=trim(STL%STLfilename), access = "STREAM", action="WRITE", &
         status = 'REPLACE', form="UNFORMATTED", convert="LITTLE_ENDIAN")

    write (STL%STLunit) 'STL test file using a simple octahedron as object                               '
    STL%Ntriangles = 8
    write (STL%STLunit) STL%Ntriangles
    call STL%STLtest_()
    close(STL%STLunit, status = 'keep' )
  end if 
else 
  STL%STLfilename = trim(STLfilename)
  STL%Ntriangles = Ntriangles
  ! open the file; this must be a little-endian file !!!
  open(unit = STL%STLunit, file=trim(STL%STLfilename), access = "STREAM", action="WRITE", &
       status = 'REPLACE', form="UNFORMATTED", convert="LITTLE_ENDIAN")

  write (STL%STLunit) STLheader
  write (STL%STLunit) STL%Ntriangles
  call STL%writeSTLfile_(MCAlist)
  close(STL%STLunit, status = 'keep' )
end if 

end function STL_constructor

!--------------------------------------------------------------------------
subroutine STL_destructor(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: STL_constructor
!! author: MDG 
!! version: 1.0 
!! date: 08/12/21
!!
!! destructor for the STL_T Class
 
IMPLICIT NONE

type(STL_T), INTENT(INOUT)  :: self 

call reportDestructor('STL_T')

end subroutine STL_destructor


!--------------------------------------------------------------------------
subroutine writeSTLfile_(self, MCAlist)
!DEC$ ATTRIBUTES DLLEXPORT :: writeSTLfile_
!! author: MDG 
!! version: 1.0 
!! date: 08/12/21
!!
!! write the triangles from the linked list to the STL file 

IMPLICIT NONE 

class(STL_T), INTENT(INOUT)       :: self
type(MCAtriangle), INTENT(IN), pointer     :: MCAlist

type(STLtriangle)                 :: tr 
integer(kind=irg)                 :: i 
real(c_float)                     :: nv(3)
type(MCAtriangle),pointer         :: tmp

tmp => MCAlist

do i=1,self%Ntriangles
  nv(:) = tmp%v1(:) + tmp%v2(:) + tmp%v3(:)
  nv = nv/sqrt(sum(nv*nv))
  tr%nv = real(nv, c_float)
  tr%v1(:) = real(tmp%v1(:), c_float)
  tr%v2(:) = real(tmp%v2(:), c_float)
  tr%v3(:) = real(tmp%v3(:), c_float)
  write (*,*) i, tmp%v1, tmp%v2, tmp%v3 
  write(self%STLunit) tr
  if (associated(tmp%next)) tmp => tmp%next
end do

end subroutine writeSTLfile_

!--------------------------------------------------------------------------
subroutine STLtest_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: STLtest_
!! author: MDG 
!! version: 1.0 
!! date: 08/12/21
!!
!! write the triangles of a simple octahedron to a file 'octahedron.stl' for testing

use mod_EMsoft

IMPLICIT NONE 

class(STL_T), INTENT(INOUT)       :: self

type(STLtriangle)                 :: tr 
real(kind=sgl),parameter          :: v(3,6) = reshape( (/ 1.0,0.0,0.0, 0.0,1.0,0.0, -1.0,0.0,0.0, &
                                                          0.0,-1.0,0.0, 0.0,0.0,1.0, 0.0,0.0,-1.0 /), (/3,6/) )
integer(kind=irg),parameter       :: nf = 8 
integer(kind=irg),parameter       :: f(3, nf) = reshape( (/ 1,2,5, 2,3,5, 3,4,5, 4,1,5, 2,1,6, 3,2,6, 4,3,6, 1,4,6 /), &
                                                         (/ 3,nf /) )
integer(kind=irg)                 :: i 
real(c_float)                     :: nv(3)

do i=1,nf 
  nv(:) = v(:,f(1,i)) + v(:,f(2,i)) + v(:,f(3,i))
  nv = nv/sqrt(sum(nv*nv))
  tr%nv = real(nv, c_float)
  tr%v1(:) = real(v(:,f(1,i)), c_float)
  tr%v2(:) = real(v(:,f(2,i)), c_float)
  tr%v3(:) = real(v(:,f(3,i)), c_float)
  write(self%STLunit) tr
end do

end subroutine STLtest_



end module mod_STL