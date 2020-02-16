! ###################################################################
! Copyright (c) 2013-2020, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_Kikuchi
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/30/20
  !!
  !! all kikuchi map related routines 

use mod_kinds
use mod_global

IMPLICIT NONE 

private

type, public :: kikuchireflection
  integer(kind=irg)                     :: hkl(3),rnum
  logical                               :: drawh,drawk,drawc
  real(kind=sgl)                        :: hlx(2),hly(2),klx(2),kly(2),clx(2),cly(2),beta,theta,Ig
  type(kikuchireflection),pointer       :: next
end type kikuchireflection


type, public :: Kikuchi_T
private 
  type(kikuchireflection), pointer     :: Kiklist

contains
private 
    procedure, pass(self) :: getListHead_ 
    procedure, pass(self) :: DeleteList_
    procedure, pass(self) :: Calcmap_
    procedure, pass(self) :: PlotBands_
    procedure, pass(self) :: KikmapPage_
    final :: Kikuchi_destructor 

    generic, public :: getListHead => getListHead_ 
    generic, public :: DeleteList => DeleteList_ 
    generic, public :: KikmapPage => KikmapPage_

!DEC$ ATTRIBUTES DLLEXPORT :: getListHead    
!DEC$ ATTRIBUTES DLLEXPORT :: DeleteList
!DEC$ ATTRIBUTES DLLEXPORT :: KikmapPage

end type Kikuchi_T

! the constructor routine for this class 
interface Kikuchi_T
  module procedure :: Kikuchi_constructor
end interface Kikuchi_T

contains

!--------------------------------------------------------------------------
type(Kikuchi_T) function Kikuchi_constructor( progdesc, EMsoft, PS ) result(Kik)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/30/20
  !!
  !! constructor for the Kikuchi_T Class 

use mod_postscript 
use mod_EMsoft

IMPLICIT NONE

character(fnlen), INTENT(IN)                    :: progdesc 
type(EMsoft_T), INTENT(INOUT)                   :: EMsoft
type(PostScript_T), INTENT(INOUT), OPTIONAL     :: PS 

integer(kind=irg)                               :: imanum
type(kikuchireflection), pointer                :: temp

if (present(PS)) then 
  imanum = 1
  PS = PostScript_T(progdesc, EMsoft, imanum)
  call PS%setpspage(0)
end if 

call Kik%DeleteList() 

allocate(Kik%Kiklist)
nullify(Kik%Kiklist%next)

end function Kikuchi_constructor

!--------------------------------------------------------------------------
subroutine Kikuchi_destructor(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/30/20
  !!
  !! clean up

IMPLICIT NONE 

type(Kikuchi_T),INTENT(INOUT)          :: self

call reportDestructor('Kikuchi_T')

call self%DeleteList()

end subroutine Kikuchi_destructor

!--------------------------------------------------------------------------
subroutine getListHead_(self, top)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/30/20
  !!
  !! return the pointer to the top of the reflection list 

IMPLICIT NONE

class(Kikuchi_T),INTENT(INOUT)                  :: self
type(kikuchireflection), pointer,INTENT(OUT)    :: top 

top => self%Kiklist

end subroutine getListHead_

!--------------------------------------------------------------------------
subroutine DeleteList_(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/30/20
  !!
  !! delete the linked list 

IMPLICIT NONE

class(Kikuchi_T),INTENT(INOUT)       :: self

type(kikuchireflection), pointer     :: temp

if (associated(self%Kiklist)) then 
    temp => self%Kiklist%next
    do while (associated(temp%next))
      deallocate(self%Kiklist)
      self%Kiklist => temp
      temp => self%Kiklist%next
    end do
    deallocate(self%Kiklist)
end if
nullify(self%Kiklist)

end subroutine DeleteList_

!--------------------------------------------------------------------------
subroutine KikmapPage_(self, cell, SG, PS, Diff, camlen)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/30/20
  !!
  !! draw Kikuchi map 

use mod_crystallography
use mod_symmetry
use mod_io
use mod_postscript
use mod_diffraction
use mod_misc

IMPLICIT NONE

class(Kikuchi_T),INTENT(INOUT)      :: self
type(Cell_T),INTENT(INOUT)          :: cell
type(SpaceGroup_T),INTENT(INOUT)    :: SG
type(PostScript_T),INTENT(INOUT)    :: PS
class(Diffraction_T),INTENT(INOUT)  :: Diff 
real(kind=sgl),INTENT(IN)           :: camlen

type(IO_T)                          :: Message 
logical                             :: again, first, newzone, nexttop
real(kind=sgl)                      :: negative(2),twopi,ggl,gg(3),igl,RR,thr, &
                                       RHOLZmax,RHOLZ(20),xo,yo,sc,pos(2),dy, oi_real(6), io_real(1)
real(kind=sgl),parameter            :: xoff(0:5)=(/0.0,3.3125,0.0,3.3125,0.0,3.3125/), &
                                       yoff(0:5)=(/6.0,6.0,3.0,3.0,0.0,0.0/), &
                                       eps = 1.0E-3
integer(kind=irg)                   :: g1i(3),g2i(3),i,numHOLZ,nref
real(kind=sgl),parameter            :: le=3.25,he=2.9375
character(1)                        :: z
    
real(kind=sgl)                      :: g1(3),g2(3),PX,PY,laL,Gmax,radius, Imax, psfh
integer(kind=irg)                   :: uvw(3),contrast

! set contrast to 1 for white lines on gray background,
! 0 for black lines on white background
 contrast = 0
 radius = 80.0
 thr = 1.E-4 
 twopi = 2.0*cPi
 Imax = 0.0

! camera length
 laL = sngl(Diff%getWaveLength()) * camlen
 oi_real(1) = sngl(Diff%getWaveLength())
 call Message%WriteValue('wavelength [nm] = ', oi_real, 1, "(F10.6)")

 oi_real(1) = camlen
 call Message%WriteValue(' L         [mm] = ', oi_real, 1, "(F10.2)")

 oi_real(1) = laL
 call Message%WriteValue('camera length lambda*L [mm nm] = ', oi_real, 1, "(F10.5)")

! what portion of reciprocal space is to be covered 
 call Message%ReadValue('Enter maximum spatial frequency to be considered [nm^-1] ', io_real, 1)
 Gmax = io_real(1)

! get the zone axis
 call GetIndex(SG%getSpaceGrouphexset(), uvw, 'd')

! get the basis vectors g1 and g2
 call cell%ShortestG(SG,uvw,g1i,g2i,i)
 g1 = float(g1i); g2 = float(g2i)
    
! get Laue center in terms of g1 and g2
  oi_real(1:3) = g1(1:3)
  oi_real(4:6) = g2(1:3)
  call Message%WriteValue('The basis vectors for this zone axis are ',oi_real, 6, "(/'g1 = ',3f10.5,/'g2 = ',3f10.5,/)") 

  psfh = PS%getpsfigheight()

! do page preamble stuff if this is a new page
! [This assumes that the PostScript file has already been opened]
  call PS%newpage(.FALSE.,'Kinematical Kikuchi Map')
  call PS%text(5.25,-0.05,'scale bar in reciprocal nm')
  call PS%textvar(5.25,psfh+0.02,'Camera Constant [nm mm]',laL)
  call PS%setfont(PSfonts(2),0.15)
  call PS%text(-0.25,psfh+0.02,'Structure File : '//trim(cell%getFileName()))
! draw frame and related stuff
  PX = 3.75
  PY = 4.00
! add other data lines to the upper left
  call PS%setfont(PSfonts(2),0.15)
  call PS%textvar(-0.25,psfh-0.18,'Acc. Voltage [kV] ',sngl(Diff%getV()))
  call PS%text(-0.25,psfh-0.38,'Zone axis ')
  call PS%PrintIndices('d',SG%getSpaceGrouphexset(),uvw(1),uvw(2),uvw(3),-0.25+1.5,psfh-0.38)
! scale bar (sc is the conversion factor from nm-1 to inches)
  sc = laL/25.4
  call PS%setlinewidth(0.020)
  call PS%line(xo+0.05,yo+0.06,xo+0.05+5.0*sc,yo+0.06)
  call PS%setfont(PSfonts(2),0.15)
  call PS%text(xo+0.05+2.5*sc,yo+0.10,'5 ')
! draw main circle graylevel 50%
  if (contrast.eq.1) then
    call PS%filledcircle(PX,PY,radius/25.4,0.5)
  else
    call PS%filledcircle(PX,PY,radius/25.4,1.0)
  end if
! plot origin of reciprocal space 
  call PS%filledcircle(PX,PY,0.03,0.0)

! compute the Kikuchimap
  call Calcmap_(self, cell, SG, PS, Diff, camlen, g1, g2, Gmax, radius, uvw, Imax)

! once that is done, we can determine the intensities to be drawn and 
! draw all reflections and lines.
    call PlotBands_(self, PS, radius, contrast, Imax)

end subroutine KikmapPage_

!--------------------------------------------------------------------------
subroutine Calcmap_(self, cell, SG, PS, Diff, camlen, g1, g2, Gmax, radius, uvw, Imax)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/30/20
  !!
  !! compute the actual map

use mod_crystallography
use mod_symmetry
use mod_io
use mod_postscript
use mod_diffraction
use mod_misc

IMPLICIT NONE

class(Kikuchi_T),INTENT(INOUT)      :: self
type(Cell_T),INTENT(INOUT)          :: cell
type(SpaceGroup_T),INTENT(INOUT)    :: SG
type(PostScript_T),INTENT(INOUT)    :: PS
class(Diffraction_T),INTENT(INOUT)  :: Diff 
real(kind=sgl),INTENT(IN)           :: camlen
real(kind=sgl),INTENT(IN)           :: g1(3), g2(3), Gmax, radius
integer(kind=irg),INTENT(IN)        :: uvw(3)
real(kind=sgl),INTENT(INOUT)        :: Imax

type(IO_T)                          :: Message 
integer(kind=irg)                   :: inmhkl(3),hc,i,j,nref,istat,inm,hh,kk,ll,ind(3),ih,ik,il,imin,hhh, oi_int(3)
real(kind=sgl)                      :: gg(3),Ig,x,alp,beta,theta,D,gpx,gpy,hlx(2),hly(2),phi,glen,gtoc(2,2),pxy(2), &
                                       da,db,dd,dec(2,2)
real(kind=dbl)                      :: rmt(3,3)
logical                             :: drawit
logical,allocatable                 :: z(:,:,:)
character(1)                        :: q
type(kikuchireflection),pointer     :: top, temp, bot
type(gnode)                         :: rlp 

 call self%getListHead(top)
 bot => top 
 nullify(bot%next)

 call Message%WriteValue('','Computing map',"(/A)")

! set the index boundaries
rmt = cell%getrmt()
 do i=1,3
   inmhkl(i) = int(1.2*Gmax/sqrt(rmt(i,i)))
 end do
 oi_int(1:3) = inmhkl(1:3)
 call Message%WriteValue('Index range ', oi_int, 3, "(3I4)")

! allocate logical array to keep track of reflections
 allocate(z(-inmhkl(1):inmhkl(1),-inmhkl(2):inmhkl(2),-inmhkl(3):inmhkl(3)))
 z = .FALSE.
 inm = maxval(inmhkl)

! initialize the geometrical parameters
 alp = atan(radius/camlen)
 Imax = 0.0
 phi = cell%Calcangle(g1,g2,'r')
 glen = cell%CalcLength(g2,'r')
 gtoc(1,1) = cell%CalcLength(g1,'r')
 gtoc(1,2) = glen*cos(phi)
 gtoc(2,1) = 0.0
 gtoc(2,2) = glen*sin(phi)
! matrix to decompose g w.r.t. g1 and g2
 da = cell%CalcDot(g1,g1,'r')
 db = cell%CalcDot(g1,g2,'r')
 dd = cell%CalcDot(g2,g2,'r')
 dec = reshape( (/dd, -db, -db, da/), (/2,2/)) / (da*dd-db**2)
 nref=0
! eliminate central spot
 z(0,0,0) = .TRUE.

! loop over all reflections
 do hh=-inmhkl(1),inmhkl(1)
  do kk=-inmhkl(2),inmhkl(2)
   do ll=-inmhkl(3),inmhkl(3)
! have we done this one already ? if so, skip
   if (z(hh,kk,ll).eqv..TRUE.) cycle
! reduce to lowest common denominator
   ind= (/ hh, kk, ll /)
   call IndexReduce(ind)
   if (z(ind(1),ind(2),ind(3)).eqv..TRUE.) cycle
! next, make sure that this reflection is allowed; if not, then
! take twice that g
   imin = 1
   if (SG%IsGAllowed(ind).eqv..FALSE.) then
! check if a multiple of g is allowed
     i=1
     imin=-1
     do while ((i.le.inm).and.(imin.eq.-1))
      if (SG%IsGAllowed(i*ind).eqv..TRUE.) then 
       imin = i
      end if
      i=i+1     
     end do
   end if
   if (imin.eq.-1) then  ! no multiple is allowed, so take next reflection
    cycle
   else
    ind = imin*ind
   end if
! make sure no higher indices get through
   if (maxval(abs(ind)).gt.inm) cycle
   if (z(ind(1),ind(2),ind(3)).eqv..TRUE.) cycle
! compute the angle between this vector and the incident beam direction
   beta = -cell%CalcAngle(float(ind),float(uvw),'c')+0.5*cPi
! if negative, take the opposite reflection (-g)
   if (beta.lt.0.0) then
     ind = -ind
     beta = -beta
   end if
! diffraction angle
   theta = Diff%CalcDiffAngle( cell, ind )
! if the entire band falls outside of the region of interest, then return
   if ((beta-theta).gt.alp) cycle
! store data in linked list, along with calculated coordinates of lines
   allocate(bot%next,stat=istat)
   nref = nref+1
   bot => bot%next
   nullify(bot%next)
   bot%theta = theta
   bot%beta  = beta
   bot%hkl = ind
! get intensity (kinematical)
   call Diff%CalcUcg(cell,ind)
   rlp = Diff%getrlp()
   bot%Ig = rlp%Vmod**2
   if (bot%Ig.gt.Imax) Imax = bot%Ig
! find the normalized Cartesian projection of g
   gpx = cell%CalcDot(float(ind),g1,'r')
   gpy = cell%CalcDot(float(ind),g2,'r')
   pxy = matmul(dec,(/gpx,gpy/))
! convert it to a Cartesian reference frame instead of g1,g2
   pxy = matmul(gtoc,pxy)
   D = sqrt(pxy(1)**2+pxy(2)**2)
   gpx = pxy(1)/D
   gpy = pxy(2)/D
! get the intercept coordinate of both lines with the projection of g
! first line
    x = radius*tan(beta+theta)/tan(alp)/25.4
    call CalcLine_(x,gpx,gpy,radius/25.4,hlx,hly,drawit)
    if (drawit.eqv..TRUE.) then
     bot%drawh=.TRUE.
     bot%hlx=hlx
     bot%hly=hly
    else
     bot%drawh = .FALSE.
    end if
! second line
    x = radius*tan(beta-theta)/tan(alp)/25.4
    call CalcLine_(x,gpx,gpy,radius/25.4,hlx,hly,drawit)
    if (drawit.eqv..TRUE.) then
     bot%drawk=.TRUE.
     bot%klx=hlx
     bot%kly=hly
    else
     bot%drawk = .FALSE.
    end if
! center line, to be drawn as a thin line
    x = radius*tan(beta)/tan(alp)/25.4
    call CalcLine_(x,gpx,gpy,radius/25.4,hlx,hly,drawit)
    if (drawit.eqv..TRUE.) then
     bot%drawc=.TRUE.
     bot%clx=hlx
     bot%cly=hly
    else
     bot%drawc = .FALSE.
    end if
! remove the multiples of those Miller indices from the list 
! so that we only keep the shortest vectors in any direction
   do hhh=-inm,inm
     ih=ind(1)*hhh
     ik=ind(2)*hhh
     il=ind(3)*hhh
     if (((abs(ih).le.inmhkl(1)).and.(abs(ik).le.inmhkl(2))).and.(abs(il).le.inmhkl(3))) then 
         z(ih,ik,il)=.TRUE.
     end if
   end do
! 
  end do
 end do
end do

oi_int(1) = nref
call Message%WriteValue('number of bands drawn : ',oi_int, 1, "(I6)")

end subroutine Calcmap_

!--------------------------------------------------------------------------
subroutine PlotBands_(self, PS, radius, contrast, Imax)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/30/20
  !!
  !! draw a Kikuchi band

use mod_io
use mod_postscript
use mod_diffraction
use mod_misc

IMPLICIT NONE

class(Kikuchi_T),INTENT(INOUT)      :: self
type(PostScript_T),INTENT(INOUT)    :: PS
real(kind=sgl),INTENT(IN)           :: radius
integer(kind=irg),INTENT(IN)        :: contrast
real(kind=sgl),INTENT(IN)           :: Imax 

type(IO_T)                          :: Message 
real(kind=sgl)                      :: V,qx,qy,CB,limit,PX,PY
character(12)                       :: txt
integer(kind=irg)                   :: i,nref, psunit
type(kikuchireflection),pointer     :: top, temp, bot

 nullify(top, temp, bot)
 call self%getListHead(top)
 psunit = PS%getpsunit() 

 call Message%printMessage('Plotting Kikuchi bands and labels',frm="(/A/)")
 PX = 3.75
 PY = 4.00

! point to the top of the linked list
 temp => top%next
 nref = 0
 limit = (1.001*radius)**2
 open(unit=30,file='temp.txt',status='unknown',form='formatted')
 if (contrast.eq.1) then 
   write (psunit,"(F12.7,' setgray')") 1.0
 else
   write (psunit,"(F12.7,' setgray')") 0.0
 end if
 call PS%setfont(PSfonts(4),0.10)
! move through the entire list 
 do while (associated(temp))
! first line
  ! V=0.015*(temp%Ig/Imax)**0.1
  V=0.01*(temp%Ig/Imax)
  write (*,*) V, temp%Ig, Imax, temp%hkl
  if (temp%drawh.eqv..TRUE.) then  
   call PS%setlinewidth(V)
   !call PS%line(PX+temp%hlx(1),PY+temp%hly(1),PX+temp%hlx(2),PY+temp%hly(2))
   call PS%line_gray(PX+temp%hlx(1),PY+temp%hly(1),PX+temp%hlx(2),PY+temp%hly(2),0.01-V)
  end if
! second line
  if (temp%drawk.eqv..TRUE.) then  
   call PS%setlinewidth(V)
   !call PS%line(PX+temp%klx(1),PY+temp%kly(1),PX+temp%klx(2),PY+temp%kly(2))
   call PS%line_gray(PX+temp%klx(1),PY+temp%kly(1),PX+temp%klx(2),PY+temp%kly(2),0.01-V)
  end if
! central line
  if (temp%drawc.eqv..TRUE.) then  
!  call PS%setlinewidth(0.004)
!  call PS%line(PX+temp%clx(1),PY+temp%cly(1),PX+temp%clx(2),PY+temp%cly(2))
! add indices along continuation of lines
   qx= PX+1.01*temp%clx(2)
   qy= PY+1.01*temp%cly(2)
   V = 180.0*atan2(temp%cly(2)-temp%cly(1),temp%clx(2)-temp%clx(1))/cPi
   write (30,"(1x,I3,1x,I3,1x,I3,1x,3(f10.5,1x))") (temp%hkl(i),i=1,3),qx,qy,V
   nref = nref+1
  end if
! move to the next reflection
  temp=>temp%next
 end do
! add indices along continuation of lines
 close(unit=30,status='keep')
 CB = (radius/25.4)**2
 write (psunit,"(F12.7,' setgray')") 0.0
 open(unit=30,file='temp.txt',status='old',form='formatted')
 do i=1,nref
   read (30,"(A12,1x,3(f10.5,1x))") txt,qx,qy,V
   if (((qx-PX)**2+(qy-PY)**2).gt.CB) then
    call PS%move(qx,qy)  ! just outside clipping ring
    write (psunit,"(1x,F10.5,' rotate')") V
    write (psunit,"(1x,'( ',A12,' ) show')") txt
    write (psunit,"(1x,F10.5,' rotate')") -V
   end if
 end do
 close(unit=30,status='delete')

end subroutine PlotBands_

!--------------------------------------------------------------------------
subroutine CalcLine_(x,px,py,rad,hlx,hly,drawit)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/30/20
  !!
  !! compute the parameters to draw a single line 

IMPLICIT NONE

real(kind=sgl)  :: x,px,py,rad,hlx(2),hly(2),det,tgm,y,qx,qy
logical         :: drawit

drawit = .FALSE.
if (abs(x).le.rad) then
    if (abs(px*py).gt.1.0e-6) then
        tgm = py/px
        y = atan2(py,px)
        qx = x*cos(y)
        qy = x*sin(y)
        det = 1.0-(1.0+tgm**2)*(1.0-(tgm*rad/(qx+tgm*qy))**2)
        if (det.gt.0.0) then  ! there is an intersection for this line so it should be drawn
             drawit = .TRUE.
             hlx(1) = (qx+tgm*qy)*(1.0-sqrt(det))/(1.0+tgm**2)
             hly(1) = qy-(hlx(1)-qx)/tgm
             hlx(2) = (qx+tgm*qy)*(1.0+sqrt(det))/(1.0+tgm**2)
             hly(2) = qy-(hlx(2)-qx)/tgm
        end if
    else  
        if (abs(px).lt.1.0e-6) then
! parallel to the x-axis 
             drawit = .TRUE.
             hlx(1) = sqrt(rad**2-x**2)
             hly(1) = x*py/abs(py)
             hlx(2) = -hlx(1)
             hly(2) = hly(1)       
        else
! parallel to the x-axis 
             drawit = .TRUE.
             hly(1) = sqrt(rad**2-x**2)
             hlx(1) = x*px/abs(px)
             hly(2) = -hly(1)
             hlx(2) = hlx(1)       
        end if
    end if
end if

end subroutine CalcLine_



end module mod_Kikuchi