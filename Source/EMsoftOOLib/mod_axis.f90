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
! SUBROUTINE: axis
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief draw an x-y plot
!
!> @details This is based on a Pascal routine of the same name 
!> that was written in 1987 by MDG and K. Mols, while graduate
!> students at the Catholic University of Leuven, Belgium.
!>
!> Given that this is somewhat older code, for now we will not extensively 
!> document each subroutine.
!
!> Original written in Pascal by:   
!>                   Koen Mols
!>                   Dept Metaalkunde en Toegepaste Materiaalkunde
!>                   KU-Leuven
!>                   Belgium
!                   
!>                   version unix 1.1      oct 1987
!    
!>  Converted to Fortran by Marc De Graef, May 2001    
!     
!>   AXIS generates linear and logarithmic x-y plots in PostScript.
!>   The figure can be scaled according to user-boundaries or, more 
!>   conveniently, in an automatic way.
!>   
!>   
!>   Syntax :  [several options of the Pascal version are incompatible
!>              with Fortran syntax, and have been removed;  axis
!>              now generates PostScript output instead of screen-directed
!>              graphics commands]
!   
!>  subroutine axis(points,xvec,yvec,xmin,xmax,ymin,ymax,
!>                  xautorange,yautorange,xmode,ymode,pmode,mark,
!>                  scalex,scaley,overplot,drawborder);
!>      
!>      
!>  item             type           description
!>
!>  points           INTEGER        number of points to draw
!>  xvec             ARRAY of REAL  x-component of the array of drawing postions
!>  yvec             ARRAY of REAL  y-component of the array of drawing postions
!>  xmin,xmax,
!>  ymin,ymax        REAL           minimum en maximum value of x- and y-component 
!>                                  upon return from axis these vars will contain 
!>                                  the actual drawing boundaries.
!>                                  e.g.   if xmin=0.3 and xmax=73 then a drawing
!>                                         will be generated from 0 to 80, thus
!>                                         xmin--> 0  and xmax--> 80
!>  xautorange,
!>  yautorange       LOGICAL        .TRUE. minimum en maximum values of each component 
!>                                         are internally defined
!>                                  .FALSE. the values passed are being used
!>                                         ( the figure will not be clipped )
!>  xmode,ymode      character*3    kind of representation (LIN,LOG)
!>  pmode            character*3    CON          : continuous line
!>                                  DOT          : markers put at each (x,y) position
!>                                  BAR          : bar-graph (vertical lines)
!>  mark             INTEGER        marker    ( for DOT )  see set_marker
!>                                  line_type ( for CON )  see set_line_style
!>  scalex           character*3    BOT,TOP,NON
!>  scaley           character*3    LEF,RIG,NON
!>  overplot         logical        .FALSE. -> draw frame
!>                                  .TRUE.  -> do not draw frame, but use same scaling
!>  drawborder       logical        .FALSE. -> do not draw the axes
!>                                  .TRUE.  -> draw the axes
!>
!>
!>  NOTES
!>      
!>      When the min and max values of both the coordinates are known (whether
!>      by passing them or by autoranging, further scaling depends on the 
!>      axis_mode :
!>      
!>      a)  linear scaling :
!>      
!>      These values are then used to trim to the actual limits,
!>      member of following set :
!>      
!>                       (0,1,1.5,2,2.5,3,4,6,8,10)
!>      
!>      b)  logarithmic scaling :
!>      
!>      The graph is always trimmed to decade boundaries.
!>  
!> Functions: stringl; power; omag; limit; border; getshift; determinestep;
!> Subroutines:  setticks; setexponent; drawborder; drawfigure; axis
!>
!> only the axis subroutine is public, everything else is private
! 
!> @todo comment the axis subroutines
!
!> @date   10/20/87 MDG/KM 1.0 original
!> @date    5/21/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/25/13 MDG 3.0 modified IO 
!> @date   02/22/24 MDG 4.0 rewrite for EMsoftOO
!-------------------------------------------------------------------------

module mod_axis
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/22/24
  !!
  !! class definition for 2D axis plots

use mod_kinds
use mod_global
use mod_postscript

IMPLICIT NONE 

! namelist for the axis type
type, public :: axistype
  real(kind=sgl)  :: axw
  real(kind=sgl)  :: xll 
  real(kind=sgl)  :: yll
end type axistype

! class definition
type, public :: axis_T
  type(Postscript_T)  :: PS
  real(kind=sgl)      :: axw,xll,yll

contains
private 
  procedure, pass(self) :: axis_
  procedure, pass(self) :: initframe_
  procedure, pass(self) :: drawborder_
  procedure, pass(self) :: setticks_
  procedure, pass(self) :: drawfigure_
  procedure, pass(self) :: setexponent_

  generic, public :: axis => axis_
end type axis_T

! the constructor routine for this class 
interface axis_T
  module procedure axis_constructor
end interface axis_T

contains

!--------------------------------------------------------------------------
! type(axis_T) function axis_constructor( PS, axw, xll, yll ) result(axis)
type(axis_T) function axis_constructor( axw, xll, yll ) result(axis)
!! author: MDG 
!! version: 1.0 
!! date: 02/22/24
!!
!! constructor for the axis_T Class
 
IMPLICIT NONE

! type(PostScript_T),INTENT(INOUT)  :: PS 
real(kind=sgl),INTENT(IN)         :: axw
real(kind=sgl),INTENT(IN)         :: xll
real(kind=sgl),INTENT(IN)         :: yll

! pass the PostScript class into the local PS
! axis%PS = PS

! set parameters for what used to bet the axistype parameters
axis%axw = axw
axis%xll = xll 
axis%yll = yll

end function axis_constructor

!--------------------------------------------------------------------------
subroutine axis_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 02/22/24
!!
!! destructor for the axis_T Class
 
IMPLICIT NONE

type(axis_T), INTENT(INOUT)  :: self 

call reportDestructor('axis_T')

end subroutine axis_destructor

!###################################################################
!############### Auxiliary Routines ################################
!###################################################################

!******************************************************************************
integer recursive function stringl(t)
  !! author: MDG
  !! version: 1.0
  !! date: 02/22/24
  !!
  !! 

IMPLICIT NONE

character(*)      :: t
integer(kind=irg) :: i

 i=len(t) 
 if (i.ne.1) then
  do while (t(i:i).eq.' ')
    i=i-1
  end do 
 end if
 stringl=i
end function stringl

!  *****************************************************************************
real recursive function power(n)
  !! author: MDG
  !! version: 1.0
  !! date: 02/22/24
  !!
  !!  raises 10 to the power n
  !!  this seems silly, but for negative powers this is a simple way to do it

IMPLICIT NONE

real(kind=dbl)       :: q
integer(kind=irg)    :: n
        
 if (n.lt.0) then 
  q = 1.D0/10**(-n)
 else
  q = 10**n
 endif
 power = sngl(q)
end function power

!  ******************************************************************************
integer recursive function omag(x)
!DEC$ ATTRIBUTES DLLEXPORT :: omag
  !! author: MDG
  !! version: 1.0
  !! date: 02/22/24
  !!
  !!  determines the order of magnitude of a real value x
  !!   
  !!                   n-1                n
  !!   O(r) = n :    10    < abs(r)  <= 10 
  !!   
  !! [this could also be done by taking the base 10 logarithm, but this routine is faster]

IMPLICIT NONE

integer(kind=irg)   :: n
real(kind=sgl)      :: x,magn,y 
 
 n=0
 magn=1.0
 y=abs(x)
 if (x.eq.0.0) then 
  omag=-1000
 else 
  if (y.lt.1) then 
   do while (y.le.magn*0.1)
    n=n-1
    magn=magn*0.1
   end do
   omag=n
  else 
   do while (magn.lt.y) 
    n=n+1
    magn=magn*10.0
   end do
   omag=n
  end if
 end if
end function omag

!  ******************************************************************************
real recursive function limit(l,n,hilo)
  !! author: MDG
  !! version: 1.0
  !! date: 02/22/24
  !!
  !! 

IMPLICIT NONE

integer(kind=irg) :: i,n
real(kind=sgl)    :: l,r,h
character(2)      :: hilo
real(kind=sgl),parameter    :: hl(10) = (/0.0,1.0,1.50,2.0,2.50,3.0,4.0,5.0,6.0,8.0/)

intent(IN)        :: n,hilo
intent(INOUT)     :: l

 l=abs(l)
 r=power(1-n)
 l=l*r
 if (hilo.eq.'hi') then 
  h=10.0 
  do i=10,2,-1
   if (l.le.hl(i)) h=hl(i)
  end do
 else
  h=8.0 
  do i=10,2,-1
   if (l.lt.hl(i)) h=hl(i-1)
  end do
 end if
 limit = h
end function limit

!  ******************************************************************************
recursive subroutine border(mode,xmin,xmax,n)
!DEC$ ATTRIBUTES DLLEXPORT :: border
!    
!     border determines the world coordinates of the upper and
!                                 n
!     lower bound in the form a*10
! 
!    LOG
!    
!    inputs are the min and max values of a row of reals
!    outputs are the standard borders when a log graph is used
!    
!    caution :  
!    suppose min=10
!    since O(min) = 1  normally the lower border would then be one less = 0
!                      (resulting in an empty decade)
!    therefore O(min*delta) is used, in which delta is a factor that will
!    change min to the first representable real value greater then min
!    so O(min*delta) will then be 2 and lowborder will be 1 
!    
!    LIN
!    
!    inputs are the min and max values of a row of reals
!    outputs are the standard borders that will be used in the
!    drawing
!    n is max(O(min),O(max))
! 

IMPLICIT NONE

character(3),INTENT(IN)         :: mode
real(kind=sgl),INTENT(INOUT)    :: xmin,xmax
integer(kind=irg),INTENT(INOUT) :: n

real(kind=sgl),parameter        :: delta=1.0000001
integer(kind=irg)               :: omin,omax
real(kind=sgl)                  :: lowborder,highborder

 if (mode.eq.'log') then
  lowborder=omag(xmin*delta)-1
  highborder=omag(xmax)
  n=0
 else 
  lowborder =0
  highborder=0
  if (abs(xmax).gt.abs(xmin)) then 
   omax=omag(xmax)
   omin=omag(xmin*delta)
  else 
   omax=omag(xmax*delta)
   omin=omag(xmin)
  end if
!  this multiplication with delta will force the limit (the closest to zero)
!  to jump to a higher order of magnitude if it is a border value
!  e.g.
!  min = 1000    max = 1200
!  normally O(min)=3 and O(max)=4)  resulting in a plot starting from zero
!  (different orders of magn)
!  but  O(min*delta) = 4 so same order of magnitude and other scale then
!  zero used for lower bound
!
  if (omax.ge.omin) then 
   n=omax 
  else 
   n=omin
  end if
! orders of magnitude are different 
  if (omax.ne.omin) then  
   if (xmin.ge.0) highborder=limit(xmax,n,'hi')
   if (xmax.le.0) lowborder=-limit(xmin,n,'hi')
   if (xmax*xmin.lt.0) then 
    if (omax.gt.omin) then 
     highborder=limit(xmax,n,'hi')
     lowborder=-1.0
    else 
     highborder=1.0
     lowborder=-limit(xmin,n,'hi')
    end if
   end if 
! orders of magnitude are the same 
  else 
   if (xmin.ge.0) then 
    lowborder=limit(xmin,n,'lo')
   else 
    lowborder=-limit(xmin,n,'hi')
   end if
   if (xmax.lt.0) then 
    highborder=-limit(xmax,n,'lo')
   else 
    highborder=limit(xmax,n,'hi')
   end if
  end if
! lowborder and highborder are between 0..100 instead of  0..1 
  n=n-1  
 end if
 xmin=lowborder
 xmax=highborder
end subroutine border

!  ******************************************************************************
recursive subroutine setticks_(self,cp,d,low,high,ts,cs,cw,ch,sh,ich,typ,sx,sy)
!DEC$ ATTRIBUTES DLLEXPORT :: setticks_

use mod_postscript

IMPLICIT NONE

class(axis_T), INTENT(INOUT)  :: self
real(kind=sgl)                :: tick2,cp,d,low,high,ts,cs,cw,ch,sh,sx,sy
integer(kind=irg)             :: i,ich
character(4)                  :: typ
real(kind=sgl),parameter      :: logtick(9) = (/0.000000000,0.301029996,0.477121255,0.602059991, &
                                                0.698970004,0.778151250,0.845098040,0.903089987,0.954242509/)

 call self%PS%setfont(PSfonts(2),2.2)
 call self%PS%setlinewidth(0.001)
! draw the tickmarks for linear scaling on the x-axis
 if (typ.eq.'xlin') then 
  call self%PS%move(cp,ts)
  call self%PS%draw(cp,ts+2)
  call self%PS%stroke
  write (self%PS%psunit,"(E14.6,E14.6,' scale')") 1.0/sx,1.0/sy
  write (self%PS%psunit,"(E14.6,' ',E14.6,' M ( ',F5.1,' ) show')") sx*(cp-1.5*sh),sy*cs,cp
  write (self%PS%psunit,"(E14.6,E14.6,' scale')") sx,sy
  tick2=cp+d/2
  if ((tick2.ge.low).and.(tick2.le.high)) then
   call self%PS%move(tick2,ts) 
   call self%PS%draw(tick2,ts+2) 
   call self%PS%stroke
   write (self%PS%psunit,"(E14.6,E14.6,' scale')") 1.0/sx,1.0/sy
   write (self%PS%psunit,"(E14.6,' ',E14.6,' M ( ',F5.1,' ) show')") sx*(tick2-1.5*sh),sy*cs,tick2
   write (self%PS%psunit,"(E14.6,E14.6,' scale')") sx,sy
  end if
 end if
! draw the tickmarks for logarithmic scaling on the x-axis
 if (typ.eq.'xlog') then 
  if (ts.eq.100) then
   call self%PS%move(cp,ts)
   call self%PS%draw(cp,ts+2)
  else
   call self%PS%move(cp,ts-1)
   call self%PS%draw(cp,ts+2)
  end if
  call self%PS%stroke
  write (self%PS%psunit,"(E14.6,E14.6,' scale')") 1.0/sx,1.0/sy
  call self%PS%text(sx*(cp-0.8*cw),sy*(cs-0.1*ch),'10')
  call self%PS%textint(sx*(cp-0.8*cw),sy*(cs+0.1*ch),'   ',int(cp))
  write (self%PS%psunit,"(E14.6,E14.6,' scale')") sx,sy
  do i=1,9
   tick2=cp+logtick(i)
   if ((tick2.ge.low).and.(tick2.le.high)) then
    call self%PS%move(tick2,ts) 
    call self%PS%draw(tick2,ts+2) 
    call self%PS%stroke
   end if
  end do
 end if
! draw the tickmarks for linear scaling on the y-axis
 if (typ.eq.'ylin') then 
  call self%PS%move(ts,cp)
  call self%PS%draw(ts+2,cp)
  call self%PS%stroke
  write (self%PS%psunit,"(E14.6,E14.6,' scale')") 1.0/sx,1.0/sy
  write (self%PS%psunit,"(E14.6,' ',E14.6,' M ( ',F5.1,' ) show')") sx*(cs-0.2*sh),sy*(cp-0.1*ch),cp
  write (self%PS%psunit,"(E14.6,E14.6,' scale')") sx,sy
  tick2=cp+d/2
  if ((tick2.ge.low).and.(tick2.le.high)) then
   call self%PS%move(ts,tick2) 
   call self%PS%draw(ts+2,tick2) 
   call self%PS%stroke
   write (self%PS%psunit,"(E14.6,E14.6,' scale')") 1.0/sx,1.0/sy
   write (self%PS%psunit,"(E14.6,' ',E14.6,' M ( ',F5.1,' ) show')") sx*(cs-0.2*sh),sy*(tick2-0.1*ch),tick2
   write (self%PS%psunit,"(E14.6,E14.6,' scale')") sx,sy
  end if
 end if
! draw the tickmarks for logarithmic scaling on the y-axis
 if (typ.eq.'ylog') then 
  if (ts.eq.100) then
   call self%PS%move(ts,cp)
   call self%PS%draw(ts+3,cp)
  else
   call self%PS%move(ts-1,cp)
   call self%PS%draw(ts+2,cp)
  end if
  call self%PS%stroke
  write (self%PS%psunit,"(E14.6,E14.6,' scale')") 1.0/sx,1.0/sy
  call self%PS%text(sx*(cs-ch),sy*(cp-0.1*ch),'10')
  call self%PS%textint(sx*(cs-ch),sy*(cp+0.1*ch),'  ',int(cp))
  write (self%PS%psunit,"(E14.6,E14.6,' scale')") sx,sy
  do i=1,9
   tick2=cp+logtick(i)
   if ((tick2.ge.low).and.(tick2.le.high)) then
    call self%PS%move(ts,tick2) 
    call self%PS%draw(ts+2,tick2) 
    call self%PS%stroke
   end if
  end do
 end if
end subroutine setticks_

!******************************************************************************)
recursive subroutine setexponent_(self, n,s)

use mod_postscript

IMPLICIT NONE

! draws the order of magnitude of the axis 
!
!              n
!            10        ( if n<>0 )
!
class(axis_T), INTENT(INOUT)        :: self
real(kind=sgl)                      :: xs,ys
integer(kind=irg)                   :: n
character(3)                        :: s

intent(IN)                          :: n,s

 if (s.eq.'BOT') then 
  xs=102
  ys=2
 end if
 if (s.eq.'TOP') then 
  xs=95
  ys=104
 end if
 if (s.eq.'LEF') then 
  xs=-5
  ys=102
 end if
 if (s.eq.'RIG') then 
  xs=108
  ys=98
 end if
 call self%PS%setfont(PSfonts(2),2.5)
 call self%PS%text(xs,ys,'*10')
 call self%PS%setfont(PSfonts(2),2.0)
 call self%PS%textint(xs+2.3,ys+1.2,' ',n)

end subroutine setexponent_

!******************************************************************************
real recursive function determinestep(m,range)

IMPLICIT NONE

character(3)              :: m
integer(kind=irg)         :: i
real(kind=sgl)            :: range
real(kind=sgl),parameter  :: x(4) = (/2.0,3.0,6.0,15.0/), y(4) = (/0.4,0.5,1.0,2.0/)
!    
!    from the border values (low and high) a step used for making
!    the numbered ticks is derived
!    
!    caution : additional ticks are made at half way 
! 
 determinestep=4.0
 if (m.eq.'log') then 
  determinestep=1
 else
  do i=1,4
   if (range.le.x(i)) determinestep=y(i)
  end do
 end if

end function determinestep

!  ******************************************************************************
integer recursive function getshift(low,high,m)
! 
!     determines the maximum length of the numbering at an axis
!     

IMPLICIT NONE

real(kind=sgl)    :: low,high
character(3)      :: m

 if (m.eq.'lin') then
  if (low.eq.-100) then 
   if (high.ne.100) then 
    getshift=1
   else 
    getshift=1
   end if
  else 
   if (low.lt.0) then 
    getshift=1
   else 
    if (high.ne.100) then 
     getshift=2
    else 
     getshift=1
    end if
   end if
  end if
! mode = 'log'
 else 
  if (omag(low).gt.omag(high)) then 
   getshift=omag(low)+2 
  else 
   getshift=omag(high)+2
  end if
 end if

end function getshift

!  ******************************************************************************)
recursive subroutine drawborder_(self,low,high,n,s,m)

use mod_postscript

IMPLICIT NONE

class(axis_T), INTENT(INOUT)        :: self
integer(kind=irg)                   :: ich,mlog,n
real(kind=sgl)                      :: sh,xl,xh,yl,yh,ts,cs,cp,ch,cw,d,low,high,q,sx,sy
character(3)                        :: s,m 
character(4)                        :: settick

 if (m.eq.'log') then
  mlog = 1
 else
  mlog = 0
 end if
! bottom or top
 if ((s.eq.'BOT').or.(s.eq.'TOP')) then
  if (m.eq.'lin') then
   settick='xlin'
  else 
   settick='xlog'
  end if
  xl=low-0.2*(high-low)
  xh=high+0.2*(high-low) 
  yl=-20
  yh=120
  if (s.eq.'BOT') then
   ts=0 
   cs=-4-4*mlog 
  else 
   ts=98 
   cs=101+mlog
  end if
 end if
! left or right
 if ((s.eq.'LEF').or.(s.eq.'RIG')) then
  if (m.eq.'lin') then
   settick='ylin'
  else 
   settick='ylog'
  end if
  yl=low-0.2*(high-low)
  yh=high+0.2*(high-low) 
  xl=-20
  xh=120
  if (s.eq.'LEF') then
   ts=0 
   cs=-5
  else 
   ts=98 
   cs=102+mlog
  end if
 end if
 if ((n.ne.0).and.(m.ne.'log')) call self%setexponent_(n,s)
! go to the new origin
 if ((s.eq.'BOT').or.(s.eq.'TOP')) then  
  q = -xl*140.0/(xh-xl) - 20.0
  write (self%PS%psunit,"(f12.6,f12.6,' T')") q,0.0
 end if
 if ((s.eq.'LEF').or.(s.eq.'RIG')) then  
  q = -yl*140.0/(yh-yl) - 20.0
  write (self%PS%psunit,"(f12.6,f12.6,' T')") 0.0,q
 end if
! and switch the scale to draw the tickmarks 
 sx = 140.0/(xh-xl)
 sy = 140.0/(yh-yl)
 write (self%PS%psunit,"(E14.6,E14.6,' scale')") sx,sy
 cw=0.02*(xh-xl)
 ch=0.06*(yh-yl)

 d=determinestep(m,abs(high-low))
 ich=getshift(low,high,m)
 sh=cw*(ich-4.0/9.0)
 if (s.eq.'RIG') sh=0.0
 if (low*high.ge.0.0) then 
  cp=low
 else 
  d=-d
  cp=0
  do while (nint(cp).ge.low) 
   call self%setticks_(cp,d,low,high,ts,cs,cw,ch,sh,ich,settick,sx,sy)
   cp=cp+d
  end do
  cp=0
  d=-d
 end if
 do while (nint(cp).le.high) 
  call self%setticks_(cp,d,low,high,ts,cs,cw,ch,sh,ich,settick,sx,sy)
  cp=cp+d
 end do
! set the scale back to what is was !!!
 sx = 1.0/sx
 sy = 1.0/sy
 write (self%PS%psunit,"(E14.6,E14.6,' scale')") sx,sy
! and return to the origin
 if ((s.eq.'BOT').or.(s.eq.'TOP')) then  
  write (self%PS%psunit,"(f12.6,f12.6,' T')") -q,0.0
 end if
 if ((s.eq.'LEF').or.(s.eq.'RIG')) then  
  write (self%PS%psunit,"(f12.6,f12.6,' T')") 0.0,-q
 end if
end subroutine drawborder_

!  ******************************************************************************
recursive subroutine drawfigure_(self,xmin,xmax,ymin,ymax,pmode,mark,points,xmode,ymode,xvec,yvec)

use mod_postscript
 
IMPLICIT NONE
 
class(axis_T), INTENT(INOUT)        :: self
integer(kind=irg)                   :: i,mark,points
real(kind=sgl)                      :: xdraw,ydraw,xmin,xmax,ymin,ymax,qx,qy,sx,sy,q
character(3)                        :: pmode,xmode,ymode
real(kind=sgl)                      :: xvec(points),yvec(points) 
  
! go to the new origin
 qx = -xmin*100.0/(xmax-xmin)
 qy = -ymin*100.0/(ymax-ymin)
 call self%PS%translate(qx,qy)
! switch the scale to draw figure
 sx = 140.0/(1.4*(xmax-xmin))
 sy = 140.0/(1.4*(ymax-ymin))
 q = 0.001/sy
 call self%PS%setlinewidth(q)
 write(self%PS%psunit,"(E14.6,' ',E14.6,' scale')") sx,sy
! clip the drawing, so that none of it appears outside of the square
! write (psunit,*) '1.0 setgray'
! call self%PS%closepath
! call self%PS%move(xmin,ymin)
! call self%PS%draw(xmax,ymin)
! call self%PS%draw(xmax,ymax)
! call self%PS%draw(xmin,ymax)
! call self%PS%clippath
! write (psunit,*) '0.0 setgray'
! the next line MUST be present, otherwise the clippath will be
! visible on the drawing !
! call self%PS%newpath
! make the drawing 
 do i=1,points 
  if (xmode.eq.'log') then 
   xdraw=log10(xvec(i))
  else 
   xdraw=xvec(i)
  end if
  if (ymode.eq.'log') then 
   ydraw=log10(yvec(i))
  else
   ydraw=yvec(i)
  end if
  if (pmode.eq.'CON') then
   if (i.eq.1) then
    call self%PS%move(xdraw,ydraw) 
   else 
    call self%PS%draw(xdraw,ydraw)
   end if
  end if
! the folowing lines are to be implemented
!         if (pmode.eq.'DOT') then
!    point(1)=xdraw
!           point(2)=ydraw
!           polymarker2d(point,1,0)
!         end if
  if (pmode.eq.'BAR') then
   call self%PS%move(xdraw,ymin)
   call self%PS%draw(xdraw,ydraw)
  end if
 end do 
 call self%PS%stroke
! set the scale back to what is was !!!
 write (self%PS%psunit,"(E14.6,' ',E14.6,' scale')") 1.0/sx,1.0/sy
 call self%PS%translate(-qx,-qy)

end subroutine drawfigure_

! ###################################################################
recursive subroutine initframe_(self,mode,db)
!! author: MDG 
!! version: 1.0 
!! date: 02/22/24
!!
!! viewing window initialization stuff

use mod_postscript

IMPLICIT NONE

class(axis_T),INTENT(INOUT)       :: self 
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
  write (self%PS%psunit,*) 'gsave'
  write (self%PS%psunit,"(F8.3,f8.3,' T')") self%xll-1.0,self%yll-1.0
  write (self%PS%psunit,"(E14.6,' dup scale')") self%axw/140.0
! set the origin
  write (self%PS%psunit,*) '20 20 T'
! define the font (default size is 2.0)
  call self%PS%setfont(PSfonts(1),2.0)
! draw the main rectangle
  if (db) call self%PS%drawrect(0.0,0.0,100.0,100.0)
 else
! reset the origin
  write (self%PS%psunit,*) 'grestore'
 end if
end subroutine initframe_

!  ******************************************************************************
recursive subroutine axis_(self,points,xvec,yvec,xmin,xmax,ymin,ymax,xautorange,yautorange, &
                          xmode,ymode,pmode,mark,scalex,scaley,overplot,db,title,xtitle,ytitle)
!DEC$ ATTRIBUTES DLLEXPORT :: axis_

use mod_postscript

IMPLICIT NONE

class(axis_T), INTENT(INOUT)   :: self
integer(kind=irg)              :: points,mark
integer(kind=irg)              :: nx,ny
real(kind=sgl)                 :: xvec(points), yvec(points), xmin, xmax, ymin, ymax,q,r
real(kind=sgl)                 :: sxmin,sxmax,symin,symax 
logical                        :: xautorange, yautorange,overplot,db
character(3)                   :: xmode, ymode, pmode, scalex, scaley
character(*)                   :: title,xtitle,ytitle

! determine the scale of the figure
 if (xautorange) then
  xmin = minval(xvec)
  xmax = maxval(xvec)
 end if
 sxmin=xmin
 sxmax=xmax
 if (yautorange) then
  ymin = minval(yvec)
  ymax = maxval(yvec)
 end if
 symin=ymin
 symax=ymax
! determine the border parameters 
 call border(xmode,xmin,xmax,nx)
 call border(ymode,ymin,ymax,ny)
! initialize graphics stuff  
 call self%initframe_('start',db)
 if (db) then 
! draw the borders if needed
  if (scalex.ne.'NON') then
   call self%drawborder_(xmin,xmax,nx,scalex,xmode)
  end if
  if (scaley.ne.'NON') then
   call self%drawborder_(ymin,ymax,ny,scaley,ymode)
  end if
 end if
 xmin=xmin*power(nx)
 xmax=xmax*power(nx)
 ymin=ymin*power(ny)
 ymax=ymax*power(ny)
! draw the curve 
 call self%drawfigure_(xmin,xmax,ymin,ymax,pmode,mark,points,xmode,ymode,xvec,yvec)
! draw the titles (this could also be done from the 
! calling program, immediately after the axis call)
! write (psunit,*) 'initclip'
 call self%PS%setfont(PSfonts(2),6.0)
 q=20.0-stringl(title)/2.0
 r=110.0
 call self%PS%text(q,r,title)
 call self%PS%setfont(PSfonts(2),4.0)
 q=50.0-stringl(xtitle)/2.0
 r=-10.0
 call self%PS%text(q,r,xtitle)
! rotate the text by 90 degrees
 write (self%PS%psunit,"(f12.6, f12.6,' T 90.0 rotate')") -10.0,50.0-stringl(ytitle)/2.0
 r=0
 call self%PS%text(r,r,ytitle)
 write (self%PS%psunit,"(' -90.0 rotate ',f12.6, f12.6,' T')") 10.0,-50.0+stringl(ytitle)/2.0
! return the proper min and max values
 if (xmode.eq.'log') then 
  xmin=power(nint(xmin))
  xmax=power(nint(xmax))
 end if
 if (ymode.eq.'log') then 
  ymin=power(nint(ymin))
  ymax=power(nint(ymax))
 end if
 call self%initframe_('stop ',.FALSE.)
end subroutine axis_

end module mod_axis