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

module mod_polyhedra
  !! author: MDG 
  !! version: 1.0 
  !! date: 07/27/21
  !!
  !! class definition for a 3D polyhedral shape (general) 
  !!
  !! This class can be used to initialize a 3D polyhedron shape, compute surface area 
  !! and volume using the 3D Projective Geometric Algebra module.  This module also 
  !! provides a method to compute the shape amplitude using Komrska's formula or provides 
  !! the shape function. The 5 Platonic and 13 Archimedian solids are available in the resources/ShapeFiles 
  !! folder and can be initialized directly; other shapes can be defined in the same 
  !! text file format:
  !!
  !! cube                           <- a generic shape name 
  !! vertex coordinates             <- starts the block of vertex coordinates 
  !! 8                              <- number of vertices 
  !!  0.5D0,-0.5D0,-0.5D0           <- x, y, z vertex coordinates in double precision
  !!  0.5D0, 0.5D0,-0.5D0
  !!  0.5D0, 0.5D0, 0.5D0
  !!  0.5D0,-0.5D0, 0.5D0
  !! -0.5D0,-0.5D0,-0.5D0
  !! -0.5D0, 0.5D0,-0.5D0
  !! -0.5D0, 0.5D0, 0.5D0
  !! -0.5D0,-0.5D0, 0.5D0
  !! vertices per face              <- starts the block of face descriptions
  !! 6                              <- number of faces 
  !! 4, 1,2,3,4                     <- # vertices, followed by the vertex labels
  !! 4, 1,5,6,2                        using a counterclockwise winding
  !! 4, 2,6,7,3
  !! 4, 3,7,8,4
  !! 4, 1,4,8,5
  !! 4, 5,8,7,6
  !!
  !! The faces MUST have the correct winding orientation in order for the 
  !! PGA3D routines to correctly generate the shape (and distinguish inside 
  !! from outside).
  !!

use mod_kinds
use mod_global
use mod_PGA3D
use mod_PGA3Dsupport

IMPLICIT NONE 

type, public :: face_T 
  integer(kind=irg)               :: nv             ! number of vertices
  integer(kind=irg),allocatable   :: verts(:)       ! vertex indices
  real(kind=dbl),allocatable      :: edgeL(:)       ! edge lengths
  real(kind=dbl)                  :: area           ! face area 
  real(kind=dbl)                  :: dorigin        ! shortest distance to origin
  integer(kind=irg),allocatable   :: otherface(:)   ! adjoining face for each edge
  type(PGA3D_T)                   :: fmv            ! face multivector
  type(PGA3D_T), allocatable      :: edgecenter(:)  ! edge center multivector
  type(PGA3D_T), allocatable      :: edge(:)        ! edge multivector
  type(PGA3D_T), allocatable      :: edgenormal(:)  ! unit normal to edge lying in face
end type face_T 

! class definition
type, public :: polyhedron_T
private 
  integer(kind=irg)               :: nvertices 
  integer(kind=irg)               :: nfaces
  integer(kind=irg)               :: nedges
  type(face_T),allocatable        :: faces(:)
  type(PGA3D_T), allocatable      :: face(:)
  type(PGA3D_T), allocatable      :: facecenter(:)
  type(PGA3D_T), allocatable      :: vertex(:)
  type(PGA3D_T), allocatable      :: facenormal(:)
  real(kind=dbl)                  :: area 
  real(kind=dbl)                  :: volume
  character(fnlen)                :: shapename 
  character(fnlen)                :: shapefile 

contains
private 
  procedure, pass(self) :: get_nvertices_
  procedure, pass(self) :: get_nfaces_
  procedure, pass(self) :: get_nedges_
  procedure, pass(self) :: get_area_
  procedure, pass(self) :: get_volume_
  procedure, pass(self) :: get_shapename_
  procedure, pass(self) :: get_shapefile_
  procedure, pass(self) :: set_nvertices_
  procedure, pass(self) :: set_nfaces_
  procedure, pass(self) :: set_nedges_
  procedure, pass(self) :: set_area_
  procedure, pass(self) :: set_volume_
  procedure, pass(self) :: set_shapename_
  procedure, pass(self) :: set_shapefile_
  procedure, pass(self) :: polyhedron_init_
  procedure, pass(self) :: polyhedron_info_
  procedure, pass(self) :: build_polyhedron_
  procedure, pass(self) :: polyhedron_shapefunction_
  procedure, pass(self) :: polyhedron_shapeamplitude_

  generic, public :: get_nvertices => get_nvertices_
  generic, public :: get_nfaces => get_nfaces_
  generic, public :: get_nedges => get_nedges_
  generic, public :: get_area => get_area_
  generic, public :: get_volume => get_volume_
  generic, public :: get_shapename => get_shapename_
  generic, public :: get_shapefile => get_shapefile_
  generic, public :: set_nvertices => set_nvertices_
  generic, public :: set_nfaces => set_nfaces_
  generic, public :: set_nedges => set_nedges_
  generic, public :: set_area => set_area_
  generic, public :: set_volume => set_volume_
  generic, public :: set_shapename => set_shapename_
  generic, public :: set_shapefile => set_shapefile_
  generic, public :: initialize_polyhedron => polyhedron_init_
  generic, public :: polyhedron_info => polyhedron_info_
  generic, public :: build_polyhedron => build_polyhedron_
  generic, public :: polyhedron_shapefunction => polyhedron_shapefunction_
  generic, public :: polyhedron_shapeamplitude => polyhedron_shapeamplitude_

end type polyhedron_T

! the constructor routine for this class 
interface polyhedron_T
  module procedure polyhedron_constructor
end interface polyhedron_T

contains

!--------------------------------------------------------------------------
type(polyhedron_T) function polyhedron_constructor( shapename, Ledge, shapefile ) result(shape)
!! author: MDG 
!! version: 1.0 
!! date: 07/27/21
!!
!! constructor for the polyhedron_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), INTENT(IN)            :: shapename 
real(kind=dbl),INTENT(IN)               :: Ledge 
character(fnlen), INTENT(IN),OPTIONAL   :: shapefile

if (present(shapefile)) then 
  call shape%polyhedron_init_( shapename, Ledge, shapefile )
else
  call shape%polyhedron_init_( shapename, Ledge )
end if 

end function polyhedron_constructor

!--------------------------------------------------------------------------
subroutine polyhedron_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 07/27/21
!!
!! destructor for the polyhedron_T Class
 
IMPLICIT NONE

type(polyhedron_T), INTENT(INOUT)  :: self 

call reportDestructor('polyhedron_T')

end subroutine polyhedron_destructor

!--------------------------------------------------------------------------
function get_nvertices_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nvertices_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get nvertices from the polyhedron_T class

IMPLICIT NONE 

class(polyhedron_T), INTENT(INOUT)     :: self
integer(kind=irg)                      :: out

out = self%nvertices

end function get_nvertices_

!--------------------------------------------------------------------------
subroutine set_nvertices_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nvertices_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set nvertices in the polyhedron_T class

IMPLICIT NONE 

class(polyhedron_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)          :: inp

self%nvertices = inp

end subroutine set_nvertices_

!--------------------------------------------------------------------------
function get_nfaces_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nfaces_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get nfaces from the polyhedron_T class

IMPLICIT NONE 

class(polyhedron_T), INTENT(INOUT)     :: self
integer(kind=irg)                      :: out

out = self%nfaces

end function get_nfaces_

!--------------------------------------------------------------------------
subroutine set_nfaces_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nfaces_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set nfaces in the polyhedron_T class

IMPLICIT NONE 

class(polyhedron_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)          :: inp

self%nfaces = inp

end subroutine set_nfaces_

!--------------------------------------------------------------------------
function get_nedges_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nedges_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get nedges from the polyhedron_T class

IMPLICIT NONE 

class(polyhedron_T), INTENT(INOUT)     :: self
integer(kind=irg)                      :: out

out = self%nedges

end function get_nedges_

!--------------------------------------------------------------------------
subroutine set_nedges_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nedges_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set nedges in the polyhedron_T class

IMPLICIT NONE 

class(polyhedron_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)          :: inp

self%nedges = inp

end subroutine set_nedges_

!--------------------------------------------------------------------------
function get_area_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_area_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get area from the polyhedron_T class

IMPLICIT NONE 

class(polyhedron_T), INTENT(INOUT)     :: self
real(kind=dbl)                         :: out

out = self%area

end function get_area_

!--------------------------------------------------------------------------
subroutine set_area_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_area_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set area in the polyhedron_T class

IMPLICIT NONE 

class(polyhedron_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)             :: inp

self%area = inp

end subroutine set_area_

!--------------------------------------------------------------------------
function get_volume_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_volume_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get volume from the polyhedron_T class

IMPLICIT NONE 

class(polyhedron_T), INTENT(INOUT)     :: self
real(kind=dbl)                         :: out

out = self%volume

end function get_volume_

!--------------------------------------------------------------------------
subroutine set_volume_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_volume_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set volume in the polyhedron_T class

IMPLICIT NONE 

class(polyhedron_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)             :: inp

self%volume = inp

end subroutine set_volume_

!--------------------------------------------------------------------------
function get_shapename_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_shapename_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get shapename from the polyhedron_T class

IMPLICIT NONE 

class(polyhedron_T), INTENT(INOUT)     :: self
character(fnlen)                       :: out

out = self%shapename

end function get_shapename_

!--------------------------------------------------------------------------
subroutine set_shapename_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_shapename_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set shapename in the polyhedron_T class

IMPLICIT NONE 

class(polyhedron_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)           :: inp

self%shapename = inp

end subroutine set_shapename_

!--------------------------------------------------------------------------
function get_shapefile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_shapefile_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! get shapefile from the polyhedron_T class

IMPLICIT NONE 

class(polyhedron_T), INTENT(INOUT)     :: self
character(fnlen)                       :: out

out = self%shapefile

end function get_shapefile_

!--------------------------------------------------------------------------
subroutine set_shapefile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_shapefile_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/21
!!
!! set shapefile in the polyhedron_T class

IMPLICIT NONE 

class(polyhedron_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)           :: inp

self%shapefile = inp

end subroutine set_shapefile_

!--------------------------------------------------------------------------
subroutine polyhedron_info_( self )
!DEC$ ATTRIBUTES DLLEXPORT :: polyhedron_info_
!! author: MDG 
!! version: 1.0 
!! date: 07/27/21
!!
!! print information about the shape 
!!

use mod_IO 

IMPLICIT NONE

class(polyhedron_T), INTENT(INOUT)  :: self 

type(IO_T)                          :: Message 
integer(kind=irg)                   :: io_int(20), i, j
real(kind=dbl)                      :: io_real(10), L
type(PGA3D_T)                       :: mv

mv = self%vertex(self%faces(1)%verts(1)).vee.self%vertex(self%faces(1)%verts(2))
L = mv%norm()

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
io_real(1) = L
call Message%WriteValue(' Edge length : ', io_real, 1)
io_real(1:2) = (/ self%area, self%volume /)
call Message%WriteValue(' Area/Volume : ', io_real, 2)
io_real(1:2) = (/ self%area/L**2, self%volume/L**3 /)
call Message%WriteValue('   for L=1   : ', io_real, 2)

call Message%printMessage('')

call Message%printMessage('Vertex tri-vectors:',frm="(/A)")
do i=1,self%nvertices 
  call self%vertex(i)%log()
end do

call Message%printMessage('Edge bi-vectors:',frm="(/A)")
do i=1,self%nfaces
  do j=1,self%faces(i)%nv
    call self%faces(i)%edge(j)%log()
  end do
end do

call Message%printMessage('Edge center tri-vectors:',frm="(/A)")
do i=1,self%nfaces
  do j=1,self%faces(i)%nv
    call self%faces(i)%edgecenter(j)%log()
  end do
end do

call Message%printMessage('Face vectors:',frm="(/A)")
do i=1,self%nfaces
  call self%faces(i)%fmv%log()
end do

call Message%printMessage('Face normal vectors:',frm="(/A)")
do i=1,self%nfaces
  mv = self%face(i) * E0123
  mv = mv.muls.(1.D0/mv%inorm())
  call mv%log()
end do

call Message%printMessage('Face edge normal bi-vectors:',frm="(/A)")
do i=1,self%nfaces
  do j=1,self%faces(i)%nv
    call self%faces(i)%edgenormal(j)%log()
  end do
end do

call Message%printMessage('Face areas & distance to origin :',frm="(/A)")
do i=1,self%nfaces
  io_int(1) = i 
  io_real(1:2) = (/ self%faces(i)%area*0.5D0, self%faces(i)%dorigin /)
  call Message%WriteValue(' ', io_int, 1, frm="(I3,$)")
  call Message%WriteValue(' ', io_real, 2)
end do

end subroutine polyhedron_info_

!--------------------------------------------------------------------------
subroutine build_polyhedron_( self, Ledge )
!DEC$ ATTRIBUTES DLLEXPORT :: build_polyhedron_
!! author: MDG 
!! version: 1.0 
!! date: 07/27/21
!!
!! build the shape parameters; to do this we use PGA3D commands 

use mod_PGA3D
use mod_IO

IMPLICIT NONE 

class(polyhedron_T), INTENT(INOUT)    :: self 
real(kind=dbl),INTENT(IN)             :: Ledge 

integer(kind=irg)                     :: i, j, i1, i2, j1, j2, k1, k2, io_int(13) 
real(kind=dbl)                        :: s, pp(3), x, y, z
character(5)                          :: str
type(PGA3D_T)                         :: mv, mv2, vsum, pt, lf 
type(IO_T)                            :: Message 

! first scale the vertices so that the first edge length of the first face equal Ledge 
mv = self%vertex(self%faces(1)%verts(1)).vee.self%vertex(self%faces(1)%verts(2))
s = Ledge / mv%norm()   ! scale factor to apply to all vertices 
do i=1,self%nvertices 
  call getpoint(self%vertex(i), x, y, z)
  self%vertex(i) = point(x*s, y*s, z*s)
end do

! build all the edges and edge centers 
do i=1,self%nfaces
  do j=1,self%faces(i)%nv
    i1 = self%faces(i)%verts(j)
    if (j.eq.self%faces(i)%nv) then 
      i2 = self%faces(i)%verts(1)
    else
      i2 = self%faces(i)%verts(j+1)
    end if
    self%faces(i)%edge(j) = self%vertex(i1).vee.self%vertex(i2)
  end do 
  do j=1,self%faces(i)%nv
    i1 = self%faces(i)%verts(j)
    if (j.eq.self%faces(i)%nv) then 
      i2 = self%faces(i)%verts(1)
    else
      i2 = self%faces(i)%verts(j+1)
    end if
    self%faces(i)%edgecenter(j) = self%vertex(i1) + self%vertex(i2)
    self%faces(i)%edgecenter(j) = self%faces(i)%edgecenter(j).muls.0.5D0
  end do 
end do 

! build the face multivectors 
do i=1,self%nfaces
  mv = self%vertex(self%faces(i)%verts(1)).vee.self%vertex(self%faces(i)%verts(2))
  self%faces(i)%fmv = mv.vee.self%vertex(self%faces(i)%verts(3))
  self%faces(i)%fmv = self%faces(i)%fmv.muls.(-1.D0)
end do 

pt = point(0.D0,0.D0,0.D0)
! find the face centers 
do j=1,self%nfaces
  mv = self%vertex(self%faces(j)%verts(1)) 
  do i=2,self%faces(j)%nv
    mv = mv + self%vertex(self%faces(j)%verts(i))
  end do 
  s = 1.D0/dble(self%faces(j)%nv)
  mv = mv.muls.s
  self%facecenter(j) = mv%normalized()
  mv2 = pt .vee. self%facecenter(j)
  self%facenormal(j) = mv2%normalized()
  self%face(j) = mv2.inner.self%facecenter(j)
end do

! and we get the edge lengths, areas, and volume
self%area = 0.D0 
self%volume = 0.D0 
do i=1,self%nfaces 
  self%faces(i)%area = 0.D0
  do j=1,self%faces(i)%nv
    i1 = self%faces(i)%verts(j)
    if (j.eq.self%faces(i)%nv) then 
      i2 = self%faces(i)%verts(1)
    else
      i2 = self%faces(i)%verts(j+1)
    end if
    ! edge length 
    mv = self%vertex(i1).vee.self%vertex(i2)
    self%faces(i)%edgeL(j) = mv%norm()

    ! unit normal to edge lying in the face 
    ! lf is the plane normal to the face p and containg the edge line ell (p . ell)
    lf = self%faces(i)%fmv.inner.self%faces(i)%edge(j)
    ! the line normal to the plane lf and containing the vertex P (lf . P)
    self%faces(i)%edgenormal(j) = lf.inner.self%vertex(self%faces(i)%verts(j))
    self%faces(i)%edgenormal(j) = self%faces(i)%edgenormal(j)%normalized()

    ! volume contribution
    mv = self%vertex(i1).vee.self%vertex(i2).vee.self%facecenter(i)
    if ((i.eq.1).and.(j.eq.1)) then 
      vsum = mv
    else 
      vsum = vsum + mv
    end if  

    ! area contribution
    self%faces(i)%area = self%faces(i)%area + mv%norm()
  end do 
  ! get the closest distance from the face to the origin 
  mv = self%face(i)%normalized()
  self%faces(i)%dorigin = mv%inorm()
end do 
self%area = sum(self%faces(:)%area) /  2.D0
self%volume = vsum%getcomp(1) / 6.D0 

! finally, for each edge in the self%faces structure, determine which other face it shares;
! this is a little tricky, since this edge will have the opposite orientation for the other face,
! assuming that the winding is consistent.  If any otherface number remains zero at the end, then there 
! must be an error in the winding of at least one of the faces...
do i=1,self%nfaces 
  self%faces(i)%otherface = 0
  do j=1,self%faces(i)%nv
    i1 = self%faces(i)%verts(j)
    if (j.eq.self%faces(i)%nv) then 
      i2 = self%faces(i)%verts(1)
    else
      i2 = self%faces(i)%verts(j+1)
    end if
! (i1, i2) is an oriented edge; next we look for the edge (i2,i1) in all other faces
    do j1=1,self%nfaces 
      do j2=1,self%faces(j1)%nv
        if ( (j1.ne.i).or.(j2.ne.j) ) then ! don't consider the current edge 
          k1 = self%faces(j1)%verts(j2)
          if (j2.eq.self%faces(j1)%nv) then 
            k2 = self%faces(j1)%verts(1)
          else
            k2 = self%faces(j1)%verts(j2+1)
          end if
          if ( (k1.eq.i2).and.(k2.eq.i1) ) then 
            self%faces(i)%otherface(j) = j1
          end if 
        end if 
      end do 
    end do
  end do 
  ! io_int(1) = i 
  ! io_int(2:self%faces(i)%nv+1) = self%faces(i)%otherface
  ! call Message%writeValue('',io_int,self%faces(i)%nv+1,frm="('face ',I3,' -> other faces ',12I3)")
end do

! test that the winding is consistent for all the faces 
do i=1,self%nfaces 
  if (minval(self%faces(i)%otherface).eq.0) then 
    write (str,"(I5)") i 
    call Message%printMessage('WARNING: face '//trim(str)//' has inconsistent winding !!! ')
  end if 
end do 

end subroutine build_polyhedron_

!--------------------------------------------------------------------------
subroutine polyhedron_init_( self, shapename, Ledge, shapefile )
!DEC$ ATTRIBUTES DLLEXPORT :: polyhedron_init_
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

use mod_EMsoft
use mod_IO 

IMPLICIT NONE 

class(polyhedron_T), INTENT(INOUT)      :: self
character(fnlen), INTENT(IN)            :: shapename 
real(kind=dbl),INTENT(IN)               :: Ledge        ! target edge length for first edge (allows for scaling)
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
  allocate(self%face(self%nfaces))
  allocate(self%facenormal(self%nfaces))
  allocate(self%facecenter(self%nfaces))
  do i=1,self%nfaces
    line = ''
    read (dataunit,"(A)") line 
    read (line,*) self%faces(i)%nv 
    allocate(self%faces(i)%verts(self%faces(i)%nv))
    read (line,*) self%faces(i)%nv, self%faces(i)%verts(1:self%faces(i)%nv)
    allocate(self%faces(i)%otherface(self%faces(i)%nv))
    allocate(self%faces(i)%edgenormal(self%faces(i)%nv))
    allocate(self%faces(i)%edgeL(self%faces(i)%nv))
    allocate(self%faces(i)%edge(self%faces(i)%nv))
    allocate(self%faces(i)%edgecenter(self%faces(i)%nv))
  end do 
  close(unit=dataunit, status='keep')
  self%nedges = sum(self%faces(1:self%nfaces)%nv) / 2
  call self%build_polyhedron_( Ledge )
else 
  call Message%printError('polyhedron_init_','Shape file '//trim(fname)//' not found')
end if 

end subroutine polyhedron_init_

!--------------------------------------------------------------------------
subroutine polyhedron_shapefunction_( self, shapearray, dims, dxyz)
!DEC$ ATTRIBUTES DLLEXPORT :: polyhedron_shapefunction_
!! author: MDG 
!! version: 1.0 
!! date: 07/28/21
!!
!! generate a discrete shape function array using pixelwise comparisons in PGA3D
!! (this is also known as the characteristic function, 1 inside and 0 outside)

IMPLICIT NONE 

class(polyhedron_T), INTENT(INOUT)      :: self
integer(kind=irg),INTENT(IN)            :: dims(3)
real(kind=sgl), INTENT(INOUT)           :: shapearray(-dims(1):dims(1),-dims(2):dims(2),-dims(3):dims(3))
real(kind=dbl),INTENT(IN)               :: dxyz 

type(PGA3D_T)                           :: pt, mv 
integer(kind=irg)                       :: i, j, k, f 
real(kind=dbl)                          :: sgn, sgnsum 

pt = point(0.D0,0.D0,0.D0)
mv = self%face(1) .wedge. pt 
sgn = dble(self%nfaces)
if (mv%getcomp(15).lt.0.D0) sgn = -sgn

shapearray = 0.D0
do i=-dims(1),dims(1)
  do j=-dims(2),dims(2)
    do k=-dims(3),dims(3)
      pt = point(dble(i)*dxyz,dble(j)*dxyz,dble(k)*dxyz)
      sgnsum = 0.D0
      do f=1,self%nfaces 
        mv = self%face(f) .wedge. pt 
        if (mv%getcomp(15).gt.0.D0) then 
          sgnsum = sgnsum + 1.D0 
        else
          sgnsum = sgnsum - 1.D0 
        end if
      end do 
      if (sgnsum.eq.sgn) shapearray(i,j,k) = 1.0
    end do 
  end do 
end do 

end subroutine polyhedron_shapefunction_

!--------------------------------------------------------------------------
subroutine polyhedron_shapeamplitude_( self, shamp, dims, dk, nthr )
!DEC$ ATTRIBUTES DLLEXPORT :: polyhedron_shapeamplitude_
!! author: MDG 
!! version: 1.0 
!! date: 08/06/21
!!
!! generate a shape amplitude using the Komrska algorithm
!! 
! @article{komrska1987a,
!   author = {Komrska, J.},
!   journal = {Optik},
!   pages = {171-183},
!   title = {Algebraic {E}xpressions of {S}hape {A}mplitudes of {P}olygons and {P}olyhedra},
!   volume = 80,
!   year = 1987}
!! 
!! This is based on a 1999 f77 version of the algorithm, written for the CTEM book illustrations.
!! The algorithm has been completely rewritten using 3D projective geometric algebra which
!! provides more flexibility and rigour to the computations (in particular in terms of 
!! directional quantities).

use omp_lib
use mod_io

IMPLICIT NONE 

class(polyhedron_T), INTENT(INOUT)      :: self
integer(kind=irg),INTENT(IN)            :: dims(3)
complex(kind=dbl),INTENT(INOUT)         :: shamp(-dims(1):dims(1)-1,-dims(2):dims(2)-1,-dims(3):dims(3)-1)
real(kind=dbl),INTENT(IN)               :: dk  ! step size in shape amplitude array
integer(kind=irg),INTENT(IN)            :: nthr

type(IO_T)                              :: Message
integer(kind=irg)                       :: i, j, k, l, f, e, vn, TID, io_int(2), zcnt
complex(kind=dbl)                       :: p, esum, ff 
real(kind=dbl)                          :: scl, knf, qq, pp, r, arg, d, ratio, dd, x, y, z
type(PGA3D_T)                           :: kvec, qn, mv, pt

zcnt = 0 

!$OMP PARALLEL NUM_THREADS(nthr) DEFAULT(SHARED) PRIVATE(TID,j,k,kvec,p,qq,f,e,qn,knf,pp,r,arg,ff,esum,mv,d,ratio,x,y,z)

TID = OMP_GET_THREAD_NUM()

!$OMP DO SCHEDULE (STATIC)
do i = -dims(1),dims(1)-1
  do j = -dims(2),dims(2)-1
    do k = -dims(3),dims(3)-1
      kvec = line(dk*dble(i), dk*dble(j), dk*dble(k))  ! this is the Fourier space frequency vector k
      if (kvec%norm().eq.0.D0) then 
        p = cmplx(self%volume,0.D0)
      else
        p = cmplx(0.D0, 0.D0)
  ! length squared of q
        qq = kvec%norm()**2
  ! loop over all faces
        do f = 1,self%nfaces
  ! first make sure that the second denominator in the Komrska equation is not zero
          qn = kvec.inner.self%facenormal(f)
          knf = qn%getcomp(0)
          pp = qq-knf**2

  ! it is zero, so use the special equation for this face 
          if (dabs(pp).lt.1.0D-8) then 
            r = 0.5D0*knf*self%faces(f)%area/qq 
            arg = knf*self%faces(f)%dorigin 
            ff = cmplx(r*dsin(arg), r*dcos(arg))
          else 
  ! its not zero, so do the general equation for this face
            esum = cmplx(0.D0, 0.D0)
  ! loop over all edges for this face
            do e = 1,self%faces(f)%nv
              mv = kvec.inner.self%faces(f)%edge(e)%normalized()
              d = mv%getcomp(0)
              arg = d*self%faces(f)%edgeL(e)*0.5D0 ! *cPi
              if (dabs(arg).lt.1.0D-8) then 
               ratio = 1.D0
              else 
               ratio = dsin(arg)/arg
              endif
              mv = kvec.inner.self%faces(f)%edgenormal(e)%normalized()
              dd = mv%getcomp(0)
              ratio = ratio*self%faces(f)%edgeL(e)*dd
              mv = point(0.D0,0.D0,0.D0).vee.self%faces(f)%edgecenter(e) 
              mv = kvec.inner.mv
              arg = mv%getcomp(0)
              esum = esum + cmplx(ratio*cos(arg),-ratio*sin(arg)) 
            end do 
            ff = -knf*esum/(pp*qq)
          endif
          p = p+ff
        end do
      end if 
      shamp(i, j, k) = p
    end do 
  end do 
  if (TID.eq.0) then 
    zcnt = zcnt+1 
    io_int = (/ nthr*zcnt, 2*dims(1)+1 /) 
    call Message%writeValue('   completed plane ', io_int, 2, frm="(I5,' out of ',I5)")
  end if 
end do 
!$OMP END DO
!$OMP END PARALLEL

call Message%printMessage(' --> done.')

end subroutine polyhedron_shapeamplitude_

end module mod_polyhedra