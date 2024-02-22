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
!
! SUBROUTINE: axonometry
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief draw a wireframe or phong-shaded surface using an axonometric projection
!
!> @details This is based on a Pascal routine of the same name 
!> that was written in 1987 by MDG and K. Mols, while graduate
!> students at the Catholic University of Leuven, Belgium.
!>
!> Given that this is somewhat older code, for now we will not extensively 
!> document each subroutine.
!>
!> Private subroutines:  setmenu, drawing, initframe, axonometry
! 
!> @todo comment the axonometry subroutines
!
!> @date   10/20/87 MDG/KM 1.0 original
!> @date    5/21/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/25/13 MDG 3.0 modified IO 
!> @date   06/09/14 MDG 4.0 made AXO and PS as arguments instead of globals
!> @date   02/22/24 MDG 5.0 update to EMsoftOO
!--------------------------------------------------------------------------

module mod_axonometry
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/22/24
  !!
  !! class definition for the axonometric projections

use mod_kinds
use mod_global

IMPLICIT NONE 

! class definition
type, public :: axonometry_T
private 
 integer(kind=irg)      :: xi,yi,beta,xmod,ymod,countx,county
 real(kind=sgl)         :: grid,scle,vscle,xstart,ystart
 logical                :: visibility
 real(kind=sgl)         :: axw,xll,yll
 character(fnlen)       :: progdesc

contains
private 
  procedure, pass(self) :: axonometry_
  procedure, pass(self) :: setmenu_ 
  procedure, pass(self) :: drawing_ 
  procedure, pass(self) :: initframe_

  generic, public :: axonometry => axonometry_

end type axonometry_T

! the constructor routine for this class 
interface axonometry_T
  module procedure axonometry_constructor
end interface axonometry_T

contains

!--------------------------------------------------------------------------
type(axonometry_T) function axonometry_constructor( progdesc, axw, xll, yll ) result(axonometry)
!! author: MDG 
!! version: 1.0 
!! date: 02/22/24
!!
!! constructor for the axonometry_T Class
 
IMPLICIT NONE

character(fnlen),INTENT(IN)   :: progdesc 
real(kind=sgl),INTENT(IN)     :: axw
real(kind=sgl),INTENT(IN)     :: xll
real(kind=sgl),INTENT(IN)     :: yll

! set parameters for what used to bet the axistype parameters
axonometry%axw = axw
axonometry%xll = xll 
axonometry%yll = yll

! initialize a bunch of internal parameters
axonometry%xi=1
axonometry%yi=2
axonometry%beta=30
axonometry%visibility=.FALSE.
axonometry%scle=1.0
axonometry%xmod=1
axonometry%ymod=1
axonometry%vscle=1.0
axonometry%progdesc=trim(progdesc)

end function axonometry_constructor

!--------------------------------------------------------------------------
subroutine axonometry_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 02/22/24
!!
!! destructor for the axonometry_T Class
 
IMPLICIT NONE

type(axonometry_T), INTENT(INOUT)  :: self 

call reportDestructor('axonometry_T')

end subroutine axonometry_destructor



!###################################################################
!###################################################################
!###################################################################
! ###################################################################
recursive subroutine setmenu_(self,what)
!! author: MDG 
!! version: 1.0 
!! date: 02/22/24
!!
!! "interactive" menu routine in a terminal window

IMPLICIT NONE

class(axonometry_T),INTENT(INOUT)   :: self
character(*),INTENT(IN)             :: what

! clear the screen
! call system('clear')

 write (*,'(3(/))')
 if (what.eq.'change_xi') then 
  write (*,"('1- xi        : ')",advance="no") 
  read (*,"(I8)") self%xi
 else 
  write (*,"('1- xi        : ',i4)") self%xi
 end if

 if (what.eq.'change_yi') then 
  write (*,"('2- yi        : ')",advance="no") 
  read (*,"(I8)") self%yi
 else 
  write (*,"('2- yi        : ',i4)") self%yi
 end if

 if (what.eq.'change_beta') then 
  write (*,"('3- beta      : ')",advance="no") 
  read (*,"(I8)") self%beta
 else 
  write (*,"('3- beta      : ',i4)") self%beta
 end if

 if (what.eq.'change_vis') then 
  self%visibility = .not. self%visibility 
 end if
 write (*,"('4- draw mode :')",advance="no") 
 if (self%visibility) then 
  write (*,"(' Phong shading')")
 else 
  write (*,"(' Wireframe')")
 end if

 if (what.eq.'change_scale') then 
  write (*,"('5- scale     : ')",advance="no") 
  read (*,*) self%scle
 else 
  write (*,"('5- scale     : ',f10.4)") self%scle
 end if

 if (what.eq.'change_xmod') then 
  write (*,"('6- xmod      : ')",advance="no") 
  read (*,"(I8)") self%xmod
 else 
  write (*,"('6- xmod      : ',i4)") self%xmod
 end if

 if (what.eq.'change_ymod') then 
  write (*,"('7- ymod      : ')",advance="no") 
  read (*,"(I8)") self%ymod
 else 
  write (*,"('7- ymod      : ',i4)") self%ymod
 end if

 if (what.eq.'change_vscale') then 
  write (*,"('8- vscale    : ')",advance="no") 
  read (*,*) self%vscle
 else 
  write (*,"('8- vscale    : ',f10.4)") self%vscle
 end if

 write (*,"('D- make drawing')") 
 write (*,"('E- Export PostScript')") 
 write (*,"('Q- quit routine')") 
 write (*,'(3(/))')

end subroutine setmenu_

! ###################################################################

recursive subroutine drawing_(self,PS,EMsoft,zz,inten,nx,ny,dmode,axname)
!! author: MDG 
!! version: 1.0 
!! date: 02/22/24
!!
!! main drawing routine that does all the work

use mod_postscript
use mod_EMsoft
use mod_io

IMPLICIT NONE

class(axonometry_T),INTENT(INOUT)         :: self
type(PostScript_T),INTENT(INOUT)          :: PS
type(EMsoft_T),INTENT(INOUT)              :: EMsoft
integer(kind=irg),INTENT(IN)              :: nx
integer(kind=irg),INTENT(IN)              :: ny
real(kind=sgl),INTENT(IN)                 :: zz(nx,ny)
real(kind=sgl),INTENT(INOUT)              :: inten(nx,ny)
character(6),INTENT(IN)                   :: dmode
character(fnlen),INTENT(IN)               :: axname

type(IO_T)                                :: Message 

integer(kind=irg)   :: i,j,kk,ip,jp
logical             :: pensw,ipp
real(kind=sgl)      :: alfa,v1,v2,w1,w2,w3,n(3),e(3),l(3),h(3),pointx,pointy,xr,yr,xp,yp,xi,yi,u1,u2, &
                       tv(4,2),pp(4,2),k_a,k_d,k_s,In_a,In_p,inmax,pn,zero,nl,hn,sx,sy, &
                       xmin,xmax,ymin,ymax,s,t,u,bb
    
! define the geometrical transformation to the image reference frame
 xi = float(self%xi)
 yi = float(self%yi)
 if (self%xi.eq.0) then 
  alfa=sngl(cPi)/2.0
 else 
  alfa=atan(yi/xi)
 end if
 v1=-sin(alfa)
 v2=cos(alfa)
 u1=sin(sngl(cPi)*self%beta/180.0)
 u2=cos(sngl(cPi)*self%beta/180.0)
 w1=-u1*v2
 w2=u1*v1
 w3=u2
! initalize the Postscript file
! open a complete file for online mode, else
! open a temporary file for export mode.
 if (dmode.eq.'export') then 
  open(UNIT=PS%psunit,FILE=trim(EMsoft%generateFilePath('EMdatapathname',axname)),STATUS='UNKNOWN',FORM='FORMATTED')
 else
  PS%psname = 'tmp.ps'
  call PS%openfile(self%progdesc, EMsoft, .TRUE.)  ! PS, progdesc, imanum, dontask
  PS%pspage = 0
 end if
 write (PS%psunit,*) 'gsave'
! set the origin
 write (PS%psunit,"(f8.3,f8.3,' T')") self%xll-1.0, self%yll-1.0
! determine the viewing window boundaries
 sx = self%scle/20.0
 sy = self%scle/20.0
 write (PS%psunit,"(E14.6,' ',E14.6,' scale')") sx,sy
! clip the drawing, so that none of it appears outside of the square
 xmin = -self%axw/2.0/sx
 xmax =  self%axw/2.0/sx
 ymin = -self%axw/2.0/sy
 ymax =  self%axw/2.0/sy
 write (PS%psunit,*) '1.0 setgray'
 call PS%closepath
 call PS%move(xmin,ymin)
 call PS%draw(xmax,ymin)
 call PS%draw(xmax,ymax)
 call PS%draw(xmin,ymax)
 call PS%clippath
 write (PS%psunit,*) '0.0 setgray'
! the next line MUST be present, otherwise the clippath will be
! visible on the drawing !
 call PS%newpath
 if (self%visibility) then
! Phong shading
! first compute the bisector unit vector for Phong shading
! (see Computer Graphics: Systems and Concepts, by R. Salmon and
! M. Slater, Addison-Wesley, 1987, p. 418)
  call Message%printMessage('Computing bisector vector', frm = "(A)")
  zero = 0.0
! eye direction
  e(1) = float(self%xi)
  e(2) = float(self%yi)
  e(3) = tan(sngl(cPi)*self%beta/180.0)*sqrt(e(1)**2+e(2)**2)
! light source direction
  l(1) = 5.0
  l(2) = -10.0
  l(3) = 50.0
! sum and normalize
  do i=1,3
   h(i) = e(i)+l(i)
  end do
  s = 0.0
  t = 0.0
  u = 0.0
  do i=1,3
   s = s + h(i)**2
   t = t + l(i)**2
   u = u + e(i)**2
  end do
  s = 1.0/sqrt(s)
  t = 1.0/sqrt(t)
  u = 1.0/sqrt(u)
  h = h*s
  l = l*t
  e = e*u
! compute all the polygon normals (triangles)
  call Message%printMessage('Computing polygon normals and intensities', frm = "(A)")
  bb = self%vscle
  inmax = -10000.0
  In_a = 0.2
  In_p = 1.0
  k_a  = 0.4
  k_d  = 0.3
  k_s  = 1.0-k_a-k_d
  pn   = 2.0
  do i=1,self%countx-1
   ip = i+1
   do j=1,self%county-1
    jp = j+1
! average normal 
    n(1)= bb*(zz(i,j)-zz(ip,j)-zz(ip,jp)+zz(i,jp))
    n(2)= bb*(zz(i,j)+zz(ip,j)-zz(ip,jp)-zz(i,jp))
    n(3) = 2.0
! normalize normals
    s = n(1)**2+n(2)**2+n(3)**2
    s = 1.0/sqrt(s)
    n = n*s
! compute Phong shading 
!   I = I_ak_a + I_p[k_d N.L + k_s (H.N)^n]
    nl=0.0
    hn=0.0
    do kk=1,3
     nl = nl + n(kk)*l(kk)
     hn = hn + n(kk)*h(kk)
    end do
! take only positive angles
    nl = max(nl,zero)
    hn = max(hn,zero)
! and compute the intensity (keep track of the maximum)
    inten(i,j)=In_a*k_a +  In_p*(k_d*nl + k_s * hn**pn)
    inmax = max(inmax,inten(i,j))
   end do
  end do
! normalize scattered intensities to maximum
  inmax = 1.0/inmax
  do i=1,self%countx
   do j=1,self%county
    inten(i,j) = max(zero,inten(i,j)*inmax)
   end do
  end do
  call Message%printMessage('Producing Postscript output', frm = "(A)")
  call PS%newpath
  do i=1,self%countx-1
   ip = i+1
   do j=1,self%county-1
    jp = j+1
! square 
    tv(1,1)=self%grid*float(i-self%countx/2)
    tv(1,2)=self%grid*(j-self%county/2)
    tv(2,1)=self%grid*float(ip-self%countx/2)
    tv(2,2)=self%grid*(j-self%county/2)
    tv(3,1)=self%grid*float(ip-self%countx/2)
    tv(3,2)=self%grid*(jp-self%county/2)
    tv(4,1)=self%grid*float(i-self%countx/2)
    tv(4,2)=self%grid*(jp-self%county/2)
    pp(1,1)=v1*tv(1,1)+v2*tv(1,2)
    pp(1,2)=w1*tv(1,1)+w2*tv(1,2)+w3*self%vscle*zz(i,j)
    pp(2,1)=v1*tv(2,1)+v2*tv(2,2)
    pp(2,2)=w1*tv(2,1)+w2*tv(2,2)+w3*self%vscle*zz(ip,j)
    pp(3,1)=v1*tv(3,1)+v2*tv(3,2)
    pp(3,2)=w1*tv(3,1)+w2*tv(3,2)+w3*self%vscle*zz(ip,jp)
    pp(4,1)=v1*tv(4,1)+v2*tv(4,2)
    pp(4,2)=w1*tv(4,1)+w2*tv(4,2)+w3*self%vscle*zz(i,jp)
    write (PS%psunit,*) inten(i,j),' setgray'
    call PS%move(pp(1,1),pp(1,2))
    call PS%draw(pp(2,1),pp(2,2))
    call PS%draw(pp(3,1),pp(3,2))
    call PS%draw(pp(4,1),pp(4,2))
    write (PS%psunit,*) 'F'
   end do
  end do
 else
! wireframe model 
! plot the x-lines
  do i=1,self%countx 
   if (mod(i,self%xmod).eq.0) then 
    call PS%newpath
    xr=self%grid*float(i-self%countx/2)
    do j=1,self%county
     yr=self%grid*(j-self%county/2)
     xp=v1*xr+v2*yr
     yp=w1*xr+w2*yr+w3*self%vscle*zz(i,j)
     if (j.eq.1) then 
      pensw=.FALSE.
     else 
      pensw=.TRUE.
     end if
     if (.not.pensw) then 
      pointx=xp
      pointy=yp
      ipp=.TRUE.
     else 
      if (ipp) then 
       call PS%move(pointx,pointy)
       ipp=.FALSE.
      end if
      call PS%draw(xp,yp)
     end if
    end do
    write (PS%psunit,*) 'S'
   end if
  end do
  call Message%printMessage('x-grid completed', frm = "(A)")
! plot y lines 
  do j=1,self%county
   if (mod(j,self%ymod).eq.0) then 
    call PS%newpath
    yr=self%grid*(j-self%county/2)
    do i=1,self%countx
     xr=self%grid*(i-self%countx/2)
     xp=v1*xr+v2*yr
     yp=w1*xr+w2*yr+w3*self%vscle*zz(i,j)
     if (i.eq.1) then 
      pensw=.FALSE.
     else 
      pensw=.TRUE.
     end if
     if (.not.pensw) then 
      pointx=xp
      pointy=yp
      ipp=.TRUE.
     else
      if (ipp) then 
       call PS%move(pointx,pointy)
       ipp=.FALSE.
      end if
      call PS%draw(xp,yp)
     end if
    end do
    write (PS%psunit,*) 'S'
   end if
  end do
  call Message%printMessage('y-grid completed', frm = "(A)")
 end if
 if (dmode.eq.'export') then
  close (unit=PS%psunit,status='keep')
 else
  call self%initframe_(PS,'stop ',.FALSE.)
  call PS%closefile()
  call Message%printMessage('Use a postscript viewing program to display the file '//PS%psname, frm = ("A"))
 end if
end subroutine drawing_

! ###################################################################
recursive subroutine initframe_(self,PS,mode,db)
!! author: MDG 
!! version: 1.0 
!! date: 02/22/24
!!
!! viewing window initialization stuff

use mod_postscript

IMPLICIT NONE

class(axonometry_T),INTENT(INOUT) :: self 
type(PostScript_T),INTENT(INOUT)  :: PS
character(5),INTENT(IN)           :: mode
logical,INTENT(IN)                :: db
 
! initialize the viewing window (from (-20,-20) to (120,120))
! (these are user coordinates, chosen so that the actual drawing
! will go from 0 to 100 along both x and y)
! On the output this means that the entire drawing, including
! axis labels and such, will always fit inside a region that
! is axw*axw inches squared. 
!
! For a single drawing on an 8.5*11 inch page with 1 inch margins 
! left and right and top and bottom 2.25 inch margins requires
! axw = 6.5 and (xll,yll) = (1.0,2.25) 
!
 if (mode.eq.'start') then 
  write (PS%psunit,*) 'gsave'
  write (PS%psunit,"(F8.3,f8.3,' T')") self%xll-1.0,self%yll-1.0
  write (PS%psunit,"(E14.6,' dup scale')") self%axw/140.0
! set the origin
  write (PS%psunit,*) '20 20 T'
! define the font (default size is 2.0)
  call PS%setfont(PSfonts(1),2.0)
! draw the main rectangle
  if (db) call PS%drawrect(0.0,0.0,100.0,100.0)
 else
! reset the origin
  write (PS%psunit,*) 'grestore'
 end if
end subroutine initframe_

! ###################################################################
recursive subroutine axonometry_(self,PS,EMsoft,zz,nx,ny,g,axname)
!DEC$ ATTRIBUTES DLLEXPORT :: axonometry_
!! author: MDG 
!! version: 1.0 
!! date: 02/22/24
!!
!! driver routine for axonometric projections; see top of file for more info 

use mod_postscript
use mod_EMsoft
use mod_io

IMPLICIT NONE

class(axonometry_T),INTENT(INOUT)    :: self
type(PostScript_T),INTENT(INOUT)     :: PS
type(EMsoft_T),INTENT(INOUT)         :: EMsoft
integer(kind=irg),INTENT(INOUT)      :: nx
integer(kind=irg),INTENT(INOUT)      :: ny
real(kind=sgl),INTENT(INOUT)         :: zz(nx,ny)
real(kind=sgl),INTENT(IN)            :: g
character(fnlen),INTENT(IN)          :: axname

type(IO_T)                           :: Message 

real(kind=sgl)                       :: inten(1)
logical                              :: more
character(1)                         :: selection 
character(fnlen)                     :: dummyname
    
! initalize some parameters; the constructor call took care of the rest
 self%countx=nx
 self%county=ny
 self%grid=g

 more = .TRUE.
 do while (more) 
  call self%setmenu_('all')
  call Message%ReadValue('Enter selection : ', selection, frm = "(A1)")
  select case (selection)
   case('1')
    call self%setmenu_('change_xi')
   case('2')
    call self%setmenu_('change_yi')
   case('3')
    call self%setmenu_('change_beta')
   case('4')
    call self%setmenu_('change_vis')
   case('5')
    call self%setmenu_('change_scale')
   case('6')
    call self%setmenu_('change_xmod')
   case('7')
    call self%setmenu_('change_ymod')
   case('8')
    call self%setmenu_('change_vscale')
   case('d','D')
    dummyname = ''
    call self%drawing_(PS,EMsoft,zz,inten,nx,ny,'online',dummyname)
   case('e','E')
    call self%drawing_(PS,EMsoft,zz,inten,nx,ny,'export',axname)
   case('q','Q')
    more = .FALSE.
   case default
  end select
 end do

end subroutine axonometry_

end module mod_axonometry