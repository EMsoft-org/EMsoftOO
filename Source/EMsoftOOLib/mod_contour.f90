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
!--------------------------------------------------------------------------
! EMsoft00:mod_contour.f90
!--------------------------------------------------------------------------
!
! MODULE: mod_contour
!
!> @author Marc De Graef, Carnegie Mellon University (but see below)
!
!> @brief module to draw contour plots
!
!> @date   ?/?/7? PVH 1.0 original f77 version by P. Van Houtte (KULeuven)
!> @date   ?/?/85 KOM 2.0 pascal version using pointers by Koen Mols
!> @date  10/10/01 MDG 3.0 f90 version (after a lot of cursing at pointers...)
!> @date  11/27/01 MDG 3.1 added kind support
!> @date   2/22/24 MDG 4.0 conversion to EMsoftOO 
!--------------------------------------------------------------------------

module mod_contour
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/22/24
  !!
  !! class definition for contour plots

use mod_kinds
use mod_global
use mod_postscript

IMPLICIT NONE 

! all variable names start with a capital C

! user defined pointer types
type Cnode 
 integer(kind=irg)      :: ix,iy
 real(kind=sgl)         :: f
 logical                :: gt
end type
       
type Cvector 
 real(kind=sgl)          :: sx,sy
 type(Ctriangle),pointer :: t1,t2
 integer(kind=irg)       :: iv
 logical                 :: border
 type(Cvector),pointer   :: next
end type
                   
type Ctriangle
 type(Cvector),pointer   :: v1,v2 
 logical                 :: to_plot
 integer(kind=irg)       :: it
 type(Ctriangle),pointer :: next
end type

type Crow_type
 type(Ctriangle),pointer :: rw
end type

! parameters specified in calling program
type Cparameters
 integer(kind=irg)       :: nx,ny 
 real(kind=sgl)          :: level,x0,y0,dx,dy
 real(kind=sgl)          :: distort(2,2)
end type


! class definition
type, public :: contour_T
! basic pointers and "pointer array"
  type(Cvector),pointer           :: Cv_root
  type(Ctriangle),pointer         :: Ct_root
  type(Crow_type),dimension(1024) :: Crow_t
  type(Cparameters)               :: Cparam
  integer(kind=irg)               :: Civ,Cit

! array (1D) of values to be contoured
  real(kind=sgl),allocatable      :: Cdata(:)

! Postscript class 
  type(Postscript_T)              :: PS 

contains
private 
  procedure, pass(self) :: contour_
  procedure, pass(self) :: make_node_
  procedure, pass(self) :: make_central_node_
  procedure, pass(self) :: make_triangle_
  procedure, pass(self) :: make_vector_
  procedure, pass(self) :: plot_contour_
  procedure, pass(self) :: plot_line_

  generic, public :: contour => contour_

end type contour_T

! the constructor routine for this class 
interface contour_T
  module procedure contour_constructor
end interface contour_T

contains

!--------------------------------------------------------------------------
type(contour_T) function contour_constructor( PS ) result(contour)
!! author: MDG 
!! version: 1.0 
!! date: 02/22/24
!!
!! constructor for the contour_T Class
 
IMPLICIT NONE

type(PostScript_T),INTENT(INOUT)    :: PS

contour%PS = PS 

end function contour_constructor

!--------------------------------------------------------------------------
subroutine contour_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 02/22/24
!!
!! destructor for the contour_T Class
 
IMPLICIT NONE

type(contour_T), INTENT(INOUT)  :: self 

call reportDestructor('contour_T')

end subroutine contour_destructor


!---------------------------------------------
recursive subroutine make_node_(self,i,j,nd)

IMPLICIT NONE

class(contour_T), INTENT(INOUT)   :: self
integer(kind=irg), INTENT(IN)     :: i,j
type(Cnode), INTENT(INOUT)        :: nd

! fill the node vector with coordinates, function value, and 
! a logical indicating whether or not the function value is
! greater than the level requested for the contour
 nd%ix=2*i
 nd%iy=2*j
 nd%f =self%Cdata(i+j*self%Cparam%nx+1)
 nd%gt=(nd%f.gt.self%Cparam%level)
end subroutine make_node_

!---------------------------------------------
recursive subroutine make_central_node_(self,i,j,nd)

IMPLICIT NONE

class(contour_T), INTENT(INOUT)   :: self
integer(kind=irg), INTENT(IN)     :: i,j
type(Cnode), INTENT(INOUT)        :: nd

! same a make_node but for the center of the grid square;
! uses an average value for the function
 nd%ix=2*i-1
 nd%iy=2*j-1
 nd%f =(self%Cdata(i+j*self%Cparam%nx+1)+self%Cdata((i-1)+j*self%Cparam%nx+1)+ &
        self%Cdata(i+(j-1)*self%Cparam%nx+1)+self%Cdata((i-1)+(j-1)*self%Cparam%nx+1))/4.0
 nd%gt=(nd%f.gt.self%Cparam%level)

end subroutine make_central_node_

!---------------------------------------------
recursive subroutine make_triangle_(self,n1,n2,n3,t)

IMPLICIT NONE

class(contour_T), INTENT(INOUT)       :: self
type(Cnode),INTENT(IN)                :: n1,n2,n3
type(Ctriangle),INTENT(INOUT),pointer :: t

! if the three corner points of the triangle are not all
! greater than or all smaller than the contour level, then
! create a triangle and insert it in the linked list, otherwise
! return a nil value.
 if ((n1%gt.eqv.n2%gt).and.(n1%gt.eqv.n3%gt)) then 
  nullify(t)
 else 
  self%Cit = self%Cit+1
  allocate(t)
  t%to_plot=.TRUE.
  nullify(t%v1)
  nullify(t%v2)
  t%it = self%Cit
  t%next => self%Ct_root
  self%Ct_root => t
 end if

end subroutine make_triangle_

!---------------------------------------------
recursive subroutine make_vector_(self,n1,n2,t1,t2)

IMPLICIT NONE

class(contour_T), INTENT(INOUT)         :: self
type(Cnode), INTENT(IN)                 :: n1,n2
type(Ctriangle), INTENT(INOUT),pointer  :: t1,t2
    
type(Cvector),pointer                   :: v 
real(kind=sgl)                          :: sx,sy

! Find the intersection point of the contour with the line
! common to the two triangles t1 and t2 (linear interpolation);
! Connect this vector to the two triangles
 if (.not.(n1%gt.eqv.n2%gt)) then 
  allocate(v)
  self%Civ=self%Civ+1
  sx=0.50*self%Cparam%dx*((self%Cparam%level-n1%f)*float(n2%ix-n1%ix)/(n2%f-n1%f)+float(n1%ix))
  sy=0.50*self%Cparam%dy*((self%Cparam%level-n1%f)*float(n2%iy-n1%iy)/(n2%f-n1%f)+float(n1%iy))
  v%sx = self%Cparam%x0+self%Cparam%distort(1,1)*sx+self%Cparam%distort(1,2)*sy
  v%sy = self%Cparam%y0+self%Cparam%distort(2,1)*sx+self%Cparam%distort(2,2)*sy
  v%border=.FALSE.
  v%t1=>t1
  v%t2=>t2
  v%iv=self%Civ
! connect t1 with vector if it is not a border point
  if (.not.associated(t1)) then 
   v%border=.TRUE.
  else 
   if (.not.associated(t1%v1)) then 
    t1%v1=>v
   else 
    t1%v2=>v
   end if
  end if
! connect t2 with vector if it is not a border point
  if (.not.associated(t2)) then 
   v%border=.TRUE.
  else 
   if (.not.associated(t2%v1)) then 
    t2%v1=>v
   else 
    t2%v2=>v
   end if
  end if
  v%next=>self%Cv_root
  self%Cv_root=>v
 end if

end subroutine make_vector_

!---------------------------------------------

recursive subroutine plot_line_(self,t,v)

use mod_postscript 

IMPLICIT NONE

class(contour_T), INTENT(INOUT)         :: self
type(Ctriangle), INTENT(INOUT),pointer  :: t
type(Cvector), INTENT(INOUT),pointer    :: v

! follow the linked list and draw all points along a single line
! Since this line may close on itself, we will have to do this
! many times, until there are no more triangles left.
 do while (associated(t)) 
   t%to_plot=.FALSE.
   if (t%v1%iv.eq.v%iv) then 
    v=>t%v2 
   else 
    v=>t%v1
   end if
   call self%PS%draw(v%sx,v%sy)
   self%Cit=self%Cit-1
   v%border=.FALSE.
   if (v%t1%it.eq.t%it) then 
    t=>v%t2 
   else 
    t=>v%t1
   end if
 end do 

end subroutine plot_line_

!---------------------------------------------
recursive subroutine plot_contour_(self)

use mod_postscript 

IMPLICIT NONE

class(contour_T), INTENT(INOUT)         :: self

type(Cvector),pointer                   :: v
type(Ctriangle),pointer                 :: t
integer(kind=irg)                       :: i

! find and plot the open contour lines  (2 passes should suffice)
 do i=1,2
  v=>self%Cv_root
  do while (associated(v)) 
   if (v%border) then 
    v%border=.FALSE.
    call self%PS%stroke
    call self%PS%move(v%sx,v%sy)
    if (associated(v%t1)) then 
      call self%plot_line_(v%t1,v)
    else
      call self%plot_line_(v%t2,v)
    end if
   end if
   v=>v%next
  end do
 end do
! 
! find and plot the closed contour lines (several passes needed)
!
 do while (self%Cit.gt.0) 
  t=>self%Ct_root
  do while (associated(t))
   if (t%to_plot) then
     v=>t%v1
     call self%PS%stroke
     call self%PS%move(v%sx,v%sy)
     if (v%t1%it.eq.t%it) then 
      nullify(v%t1)
     else 
      nullify(v%t2)
     end if
     call self%plot_line_(t,v)
   end if
   t=>t%next
  end do
 end do
!
! dispose of both linked lists
!
 do while (associated(self%Ct_root)) 
  t=>self%Ct_root
  self%Ct_root=>t%next
  deallocate(t)
 end do
 do while (associated(self%Cv_root))
  v=>self%Cv_root
  self%Cv_root=>v%next
  deallocate(v)
 end do

end subroutine plot_contour_


!---------------------------------------------
recursive subroutine contour_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: contour_                      
!                                              
! definition of vectors, nodes and triangles
!                                              
!            ul       ur                       
!            |\       /                         
!            | \  tt /                          
!            |  4   5                           
!            |   \ /                            
!            1 tl c  tr                         
!            |   / \                            
!            |  3   6                           
!            | /  tb \ 
!             /___2___\
!            ll       lr                       
!                                              
! Each grid square is triangulated by four triangles
! which are labelled tl, tt, tb, and tr.  Each grid
! square has four corners (ll, ul, lr, ur) and a center (c).
! The algorithm computes for each of the six connecting
! vectors whether or not there is an intersection of the 
! requested contour.  If there is one, then the corresponding
! coordinates are added to a linked list which keeps track
! of the two neighbouring triangles that make up the vector.
!
! At the end the entire linked list is traversed to draw all
! contours, both the ones that intersect the border and the
! ones that are fully inside the border.
!
! There may be better ways to do this, but it definitely works.
!

use mod_io

IMPLICIT NONE 

class(contour_T), INTENT(INOUT)         :: self

type(IO_T)                              :: Message 

type(Cnode)                             :: nd_ll,nd_ul,nd_lr,nd_ur,nd_c
type(Ctriangle),pointer                 :: t_l,t_r,t_b,t_t,dummy
integer(kind=irg)                       :: Cit_save, i, j, oi_int(2)
real(kind=sgl)                          :: oi_real(1)

! initialize vector and triangle counters
 self%Civ = 0
 self%Cit = 0

! initialize vector and triangle linked lists 
 if (.not.(associated(self%Cv_root))) then
  allocate(self%Cv_root)
  nullify(self%Cv_root)
 end if
 if (.not.(associated(self%Ct_root))) then
  allocate(self%Ct_root)
  nullify(self%Ct_root)
 end if 
! allocate a dummy triangle
 allocate(dummy)
 nullify(dummy)     
! loop over all rows in the grid
 do j=1,self%Cparam%ny-1 
  call self%make_node_(0,j-1,nd_ll)                     ! lower left node
  call self%make_node_(0,j  ,nd_ul)                     ! upper left node
! there is no neighbouring triangle at the left 
  nullify(t_r)                  
! loop over all colums in the row
  do i=1,self%Cparam%nx-1
   call self%make_node_(i,j-1,nd_lr)                    ! lower right node
   call self%make_node_(i,j  ,nd_ur)                    ! upper right node
   call self%make_central_node_(i,j,nd_c)               ! central node
             
   call self%make_triangle_(nd_ll,nd_c,nd_ul,t_l)       ! left triangle
   call self%make_vector_(nd_ll,nd_ul,t_r,t_l)          ! vector 1 
             
   call self%make_triangle_(nd_ll,nd_c,nd_lr,t_b)       ! bottom triangle
   call self%make_triangle_(nd_ul,nd_c,nd_ur,t_t)       ! top triangle
   call self%make_triangle_(nd_lr,nd_c,nd_ur,t_r)       ! right triangle
             
   call self%make_vector_(nd_ll,nd_lr,self%Crow_t(i)%rw,t_b) ! vector 2 
   self%Crow_t(i)%rw=>t_t                              ! keep for the next row

   call self%make_vector_(nd_ll,nd_c,t_l,t_b)           ! vector 3
   call self%make_vector_(nd_ul,nd_c,t_l,t_t)           ! vector 4
   call self%make_vector_(nd_ur,nd_c,t_r,t_t)           ! vector 5
   call self%make_vector_(nd_lr,nd_c,t_r,t_b)           ! vector 6
             
   nd_ll=nd_lr                                    ! copy lower right to lower left
   nd_ul=nd_ur                                    ! copy upper right to upper left
  end do  ! end loop over i
  call self%make_vector_(nd_ll,nd_ul,t_r,dummy)         ! last vector of a row
 end do   ! end loop over j 
     
! close the upperside 
 call self%make_node_(0,self%Cparam%ny-1,nd_ll)
 do i=1,self%Cparam%nx-1
  call self%make_node_(i,self%Cparam%ny-1,nd_lr)
  call self%make_vector_(nd_ll,nd_lr,self%Crow_t(i)%rw,dummy)  
  nd_ll=nd_lr
 end do 

! and plot it 
 Cit_save = self%Cit 
 call self%plot_contour_()
 call self%PS%stroke  ! to make sure the last line is drawn

oi_real(1)= self%Cparam%level
call Message%WriteValue(' Level : ',oi_real, 1, "(F10.5,';')",advance="no")
oi_int(1) = self%Civ 
oi_int(2) = Cit_save
call Message%WriteValue('   vectors/triangles : ', oi_int, 2, "(I6,'/',I6)")

end subroutine contour_



end module mod_contour