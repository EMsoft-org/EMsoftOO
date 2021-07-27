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

module mod_shapes
  !! author: MDG 
  !! version: 1.0 
  !! date: 07/27/21
  !!
  !! class definition for a polyhedron shape 
  !!
  !! This class can be used to initialize a 3D polyhedron shape, compute surface area 
  !! and volume using the 3D Projective Geometric Algebra module.  This module also 
  !! provides a method to compute the shape amplitude using Komrska's formula or a 
  !! simple 3D FFT. 

use mod_kinds
use mod_global
use mod_PGA3D
use mod_PGA3Dsupport

IMPLICIT NONE 

type, public :: face_T 
  integer(kind=irg)               :: nv 
  integer(kind=irg),allocatable   :: verts(:)
end type face_T 

! class definition
type, public :: shape_T
private 
  integer(kind=irg)               :: nvertices 
  integer(kind=irg)               :: nfaces
  integer(kind=irg)               :: nedges
  type(face_T),allocatable        :: faces(:)
  type(PGA3D_T), allocatable      :: facecenter(:)
  type(PGA3D_T), allocatable      :: vertex(:)
  type(PGA3D_T), allocatable      :: facenormal(:)
  type(PGA3D_T), allocatable      :: edge(:)
  real(kind=dbl)                  :: area 
  real(kind=dbl)                  :: volume
  character(fnlen)                :: shapename 
  character(fnlen)                :: shapefile 

contains
private 
  procedure, pass(self) :: shape_init_
  procedure, pass(self) :: shape_info_
  procedure, pass(self) :: build_shape_

  generic, public :: initialize_shape => shape_init_
  generic, public :: shape_info => shape_info_
  generic, public :: build_shape => build_shape_

end type shape_T

! the constructor routine for this class 
interface shape_T
  module procedure shape_constructor
end interface shape_T

contains

!--------------------------------------------------------------------------
type(shape_T) function shape_constructor( shapename, shapefile ) result(shape)
!! author: MDG 
!! version: 1.0 
!! date: 07/27/21
!!
!! constructor for the shape_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), INTENT(IN)            :: shapename 
character(fnlen), INTENT(IN),OPTIONAL   :: shapefile

if (present(shapefile)) then 
  call shape%shape_init_( shapename, shapefile )
else
  write (*,*) 'calling shape_init_()'
  call shape%shape_init_( shapename )
end if 

end function shape_constructor

!--------------------------------------------------------------------------
subroutine shape_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 07/27/21
!!
!! destructor for the shape_T Class
 
IMPLICIT NONE

type(shape_T), INTENT(INOUT)  :: self 

call reportDestructor('shape_T')

end subroutine shape_destructor

!--------------------------------------------------------------------------
subroutine shape_info_( self )
!DEC$ ATTRIBUTES DLLEXPORT :: shape_info_
!! author: MDG 
!! version: 1.0 
!! date: 07/27/21
!!
!! print information about the shape 
!!

use mod_IO 

IMPLICIT NONE

class(shape_T), INTENT(INOUT)  :: self 

type(IO_T)                    :: Message 
integer(kind=irg)             :: io_int(20)
real(kind=dbl)                :: io_real(10)

call Message%printMessage(' Shape Name : '//trim(self%shapename))
call Message%printMessage(' Shape File : '//trim(self%shapefile),frm="(A/)")

io_int(1) = self%nvertices 
call Message%WriteValue(' # vertices  : ', io_int, 1)
io_int(1) = self%nfaces 
call Message%WriteValue(' # faces     : ', io_int, 1)
io_int(1) = self%nedges
call Message%WriteValue(' # edges     : ', io_int, 1)
io_int(1) = self%nvertices + self%nfaces - self%nedges
call Message%WriteValue(' V+F-E (should be 2) ', io_int, 1, frm="(I3/)")

io_real(1:2) = (/ self%area, self%volume /)
call Message%WriteValue(' Area/Volume : ', io_real, 2)

end subroutine shape_info_


!--------------------------------------------------------------------------
subroutine build_shape_( self )
!DEC$ ATTRIBUTES DLLEXPORT :: build_shape_
!! author: MDG 
!! version: 1.0 
!! date: 07/27/21
!!
!! build the shape starting from the vertices; to do this we use PGA3D commands 

use mod_PGA3D

IMPLICIT NONE 

class(shape_T), INTENT(INOUT)   :: self 

integer(kind=irg)               :: i, j, icnt, i1, i2 
real(kind=dbl)                  :: s
type(PGA3D_T)                   :: mv, mv2, vsum 

! first we normalize all the vertices
do i=1,self%nvertices 
  self%vertex(i) = self%vertex(i)%normalized()
end do 

! find the face centers 
do j=1,self%nfaces
  mv = self%vertex(self%faces(j)%verts(1)) 
  do i=2,self%faces(j)%nv
    mv = mv + self%vertex(self%faces(j)%verts(i))
  end do 
  s = 1.D0/dble(self%faces(j)%nv)
  mv = mv.muls.s
  self%facecenter(j) = mv%normalized()
  ! call self%facecenter(j)%log()
end do

! then we build all the edges
icnt = 1
do i=1,self%nfaces
  do j=1,self%faces(i)%nv
    self%edge(icnt) = self%vertex(self%faces(i)%verts(j)).vee.self%vertex(self%faces(i)%verts(mod(j+1,self%faces(i)%nv)+1))
    ! call self%edge(icnt)%log()
    icnt = icnt + 1
  end do 
end do 

! and we get the area and volume
self%area = 0.D0 
self%volume = 0.D0 
do i=1,self%nfaces 
  do j=1,self%faces(i)%nv
    i1 = self%faces(i)%verts(j)
    if (j.eq.self%faces(i)%nv) then 
      i2 = self%faces(i)%verts(1)
    else
      i2 = self%faces(i)%verts(j+1)
    end if
    mv = self%vertex(i1).vee.self%vertex(i2).vee.self%facecenter(i)
    if ((i.eq.1).and.(j.eq.1)) then 
      vsum = mv
    else 
      vsum = vsum + mv
    end if  
    call vsum%log()
    self%area = self%area + mv%norm()
  end do 
end do 
self%area = self%area /  2.D0
self%volume = vsum%getcomp(1) / 6.D0 

end subroutine build_shape_

!--------------------------------------------------------------------------
subroutine shape_init_( self, shapename, shapefile )
!DEC$ ATTRIBUTES DLLEXPORT :: shape_init_
!! author: MDG 
!! version: 1.0 
!! date: 07/27/21
!!
!! populate the class with shape information
!!
!! We use 3D Projective Geometric Algebra to compute the face normals, the total area 
!! and the volume.  The routine then stores the parameters needed to compute the shape
!! amplitude using the Komrska formula.  A shape file must have the following format:
!!
! cube
! vertex coordinates
! 8
!  0.5D0,-0.5D0,-0.5D0
!  0.5D0, 0.5D0,-0.5D0
!  0.5D0, 0.5D0, 0.5D0
!  0.5D0,-0.5D0, 0.5D0
! -0.5D0,-0.5D0,-0.5D0
! -0.5D0, 0.5D0,-0.5D0
! -0.5D0, 0.5D0, 0.5D0
! -0.5D0,-0.5D0, 0.5D0
! vertices per face
! 6
! 4, 1,2,3,4
! 4, 1,5,6,2
! 4, 2,6,7,3
! 4, 3,7,8,4
! 4, 1,4,8,5
! 4, 5,8,7,6


use mod_EMsoft
use mod_IO 

IMPLICIT NONE 

class(shape_T), INTENT(INOUT)           :: self
character(fnlen), INTENT(IN)            :: shapename 
character(fnlen), INTENT(IN),OPTIONAL   :: shapefile

type(EMsoft_T)                          :: EMsoft
type(IO_T)                              :: Message 
type(PGA3D_T)                           :: pt 

character(fnlen)                        :: fname, line, pn
logical                                 :: fexists 
real(kind=dbl)                          :: x, y, z 
integer(kind=irg)                       :: i 

pn = ''
EMsoft = EMsoft_T(pn, pn, silent=.TRUE.)

if (.not.present(shapefile)) then 
  fname = EMsoft%generateFilePath('Resourcepathname','ShapeFiles')
  fname = trim(fname)//'/'//trim(shapename)
  fname = EMsoft%toNativePath(fname)
else
  fname = trim(shapefile)
end if 

write(*,*) 'looking for file '//trim(fname)

self%shapename = trim(shapename)
self%shapefile = trim(fname) 

inquire(file=trim(fname),exist=fexists)

if (fexists.eqv..TRUE.) then 
  open(unit=dataunit, file=trim(fname), status='old', form='formatted')
  read (dataunit,"(A)") line   ! shape name 
  read (dataunit,"(A)") line   ! vertex coordinates 
  read (dataunit,*) self%nvertices
  allocate(self%vertex(self%nvertices))
  do i=1,self%nvertices 
    read (dataunit,*) x, y, z 
    pt = point(x, y, z)
    self%vertex(i) = pt
    ! call pt%log()
  end do 
  read (dataunit,"(A)") line   ! vertex coordinates 
  read (dataunit,*) self%nfaces
  allocate(self%faces(self%nfaces))
  allocate(self%facecenter(self%nfaces))
  do i=1,self%nfaces
    line = ''
    read (dataunit,"(A)") line 
    read (line,*) self%faces(i)%nv 
    allocate(self%faces(i)%verts(self%faces(i)%nv))
    read (line,*) self%faces(i)%nv, self%faces(i)%verts(1:self%faces(i)%nv)
  end do 
  close(unit=dataunit, status='keep')
  self%nedges = sum(self%faces(1:self%nfaces)%nv) / 2
  allocate(self%edge(2*self%nedges))  ! each edge will need to be counted twice !!!
  call self%build_shape_()
else 
  call Message%printError('shape_init_','Shape file '//trim(fname)//' not found')
end if 

end subroutine shape_init_



end module mod_shapes