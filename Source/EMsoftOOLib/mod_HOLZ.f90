! ###################################################################
! Copyright (c) 2013-2022, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_HOLZ
  !! author: MDG
  !! version: 1.0
  !! date: 01/28/20
  !!
  !! all HOLZ-related routines

use mod_kinds
use mod_global

IMPLICIT NONE

private


! a structure that contains all the relevant HOLZ geometry information;
! it can be filled by calling the GetHOLZGeometry subroutine
type, public :: HOLZentries
  real(kind=sgl)        :: g1(3),g2(3),g3(3),gx(3),gy(3),LC1,LC2,H,FNr(3),FNg(3),gp(2),gtoc(2,2),gshort(3)
  integer(kind=irg)     :: uvw(3),FN(3)
end type HOLZentries

! type definitions for kinematical HOLZ patterns
type, public :: HOLZreflection
  integer(kind=irg)             :: hkl(3),n1,n2,N
  logical                       :: draw,dbdiff
  real(kind=sgl)                :: hlphi,hlx(2),hly(2),sg,Ig,pxy(2)
  type(HOLZreflection),pointer  :: next
end type HOLZreflection

type, public :: HOLZvartype
  real(kind=sgl)                :: g1(3),g2(3),g3(3),H,FNg(3),FNr(3),gshort(3),gp(3),LC1,LC2,thickness,rectangle(2), &
                                   PX,PY,thetac,laL,Gmax,Imax,gtoc(2,2),glen,phi,CBEDrad,CBEDsc
  integer(kind=irg)             :: uvw(3),FN(3)
end type HOLZvartype


! the main HOLZ class definition
type, public :: HOLZ_T
private
  type(HOLZentries)                 :: HOLZentry
  type(HOLZreflection), pointer     :: HOLZlist
  type(HOLZvartype)                 :: HOLZvar

contains
private
    procedure, pass(self) :: getListHead_
    procedure, pass(self) :: DeleteList_
    procedure, pass(self) :: setrectangle_
    procedure, pass(self) :: setlaL_
    procedure, pass(self) :: setImax_
    procedure, pass(self) :: setGmax_
    procedure, pass(self) :: getrectangle_
    procedure, pass(self) :: getlaL_
    procedure, pass(self) :: getImax_
    procedure, pass(self) :: getGmax_
    procedure, pass(self) :: HOLZPage_
    procedure, pass(self) :: ShortestGFOLZ_
    procedure, pass(self) :: CalcHOLZ_
    procedure, pass(self) :: ReCalcHOLZ_
    procedure, pass(self) :: PlotHOLZ_
    procedure, pass(self) :: PlotHOLZlines_
    procedure, pass(self) :: CalcsgHOLZ_
    procedure, pass(self) :: GetHOLZcoordinates_
    procedure, pass(self) :: GetHOLZGeometry_

    final :: HOLZ_destructor

    generic, public :: getListHead => getListHead_
    generic, public :: DeleteList => DeleteList_
    generic, public :: setrectangle => setrectangle_
    generic, public :: setlaL => setlaL_
    generic, public :: setImax => setImax_
    generic, public :: setGmax => setGmax_
    generic, public :: getrectangle => getrectangle_
    generic, public :: getlaL => getlaL_
    generic, public :: getImax => getImax_
    generic, public :: getGmax => getGmax_
    generic, public :: HOLZPage => HOLZPage_
    generic, public :: ShortestGFOLZ => ShortestGFOLZ_
    generic, public :: CalcHOLZ => CalcHOLZ_
    generic, public :: ReCalcHOLZ => ReCalcHOLZ_
    generic, public :: PlotHOLZ => PlotHOLZ_
    generic, public :: PlotHOLZlines => PlotHOLZlines_
    generic, public :: CalsgHOLZ => CalcHOLZ_
    generic, public :: GetHOLZcoordinates => GetHOLZcoordinates_
    generic, public :: GetHOLZGeometry => GetHOLZGeometry_
end type HOLZ_T

! the constructor routine for this class
interface HOLZ_T
  module procedure :: HOLZ_constructor
end interface HOLZ_T

contains

!--------------------------------------------------------------------------
type(HOLZ_T) function HOLZ_constructor( progdesc, EMsoft, PS ) result(HOLZ)
!DEC$ ATTRIBUTES DLLEXPORT :: HOLZ_constructor
  !! author: MDG
  !! version: 1.0
  !! date: 01/28/20
  !!
  !! constructor for the HOLZ_T Class

use mod_postscript
use mod_EMsoft

IMPLICIT NONE

character(fnlen), INTENT(IN), OPTIONAL          :: progdesc
type(EMsoft_T), INTENT(INOUT), OPTIONAL         :: EMsoft
type(PostScript_T), INTENT(INOUT), OPTIONAL     :: PS

integer(kind=irg)                               :: imanum
type(HOLZreflection), pointer                   :: temp

if (present(PS)) then
  imanum = 1
  PS = PostScript_T(progdesc, EMsoft, imanum)
  call PS%setpspage(0)
end if

! call HOLZ%DeleteList()

allocate(HOLZ%HOLZlist)
nullify(HOLZ%HOLZlist%next)

end function HOLZ_constructor

!--------------------------------------------------------------------------
subroutine HOLZ_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: HOLZ_destructor
!! author: MDG
!! version: 1.0
!! date: 02/02/20
!!
!! destructor for the HOLZ_T Class

IMPLICIT NONE

type(HOLZ_T), INTENT(INOUT)  :: self

call reportDestructor('HOLZ_T')

call self%DeleteList()

end subroutine HOLZ_destructor

!--------------------------------------------------------------------------
subroutine getListHead_(self, top)
!DEC$ ATTRIBUTES DLLEXPORT :: getListHead_
  !! author: MDG
  !! version: 1.0
  !! date: 01/28/20
  !!
  !! return the pointer to the top of the reflection list

IMPLICIT NONE

class(HOLZ_T),INTENT(INOUT)                 :: self
type(HOLZreflection), pointer,INTENT(OUT)   :: top

top => self%HOLZlist

end subroutine getListHead_

!--------------------------------------------------------------------------
subroutine DeleteList_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: DeleteList_
  !! author: MDG
  !! version: 1.0
  !! date: 01/28/20
  !!
  !! delete the linked list

IMPLICIT NONE

class(HOLZ_T),INTENT(INOUT)       :: self
type(HOLZreflection), pointer     :: temp

if (associated(self%HOLZlist)) then
    temp => self%HOLZlist ! %next
    do while (associated(temp%next))
      deallocate(self%HOLZlist)
      self%HOLZlist => temp
      temp => self%HOLZlist%next
    end do
    deallocate(self%HOLZlist)
end if
nullify(self%HOLZlist)

end subroutine DeleteList_

!--------------------------------------------------------------------------
function getrectangle_(self) result(t)
!DEC$ ATTRIBUTES DLLEXPORT :: getrectangle_
  !! author: MDG
  !! version: 1.0
  !! date: 01/28/20
  !!
  !! return rectangle parameters

IMPLICIT NONE

class(HOLZ_T),INTENT(INOUT)    :: self
real(kind=sgl)                 :: t(2)

t = self%HOLZvar%rectangle

end function getrectangle_

!--------------------------------------------------------------------------
function getlaL_(self) result(t)
!DEC$ ATTRIBUTES DLLEXPORT :: getlaL_
  !! author: MDG
  !! version: 1.0
  !! date: 01/28/20
  !!
  !! return  lambda * L

IMPLICIT NONE

class(HOLZ_T),INTENT(INOUT)    :: self
real(kind=sgl)                 :: t

t = self%HOLZvar%laL

end function getlaL_

!--------------------------------------------------------------------------
function getImax_(self) result(t)
!DEC$ ATTRIBUTES DLLEXPORT :: getImax_
  !! author: MDG
  !! version: 1.0
  !! date: 01/28/20
  !!
  !! return Imax value

IMPLICIT NONE

class(HOLZ_T),INTENT(INOUT)    :: self
real(kind=sgl)                 :: t

t = self%HOLZvar%Imax

end function getImax_

!--------------------------------------------------------------------------
function getGmax_(self) result(t)
!DEC$ ATTRIBUTES DLLEXPORT :: getGmax_
  !! author: MDG
  !! version: 1.0
  !! date: 01/28/20
  !!
  !! return Gmax value

IMPLICIT NONE

class(HOLZ_T),INTENT(INOUT)    :: self
real(kind=sgl)                 :: t

t = self%HOLZvar%Gmax

end function getGmax_

! !--------------------------------------------------------------------------
! function get(self) result(t)
!   !! author: MDG
!   !! version: 1.0
!   !! date: 01/28/20
!   !!
!   !! return

! IMPLICIT NONE

! class(HOLZ_T),INTENT(INOUT)    :: self
! real(kind=sgl)                 :: t

! t = self%HOLZvar%

! end function get

!--------------------------------------------------------------------------
subroutine setrectangle_(self, t)
!DEC$ ATTRIBUTES DLLEXPORT :: setrectangle_
  !! author: MDG
  !! version: 1.0
  !! date: 01/28/20
  !!
  !! set rectangle parameter

IMPLICIT NONE

class(HOLZ_T),INTENT(INOUT)    :: self
real(kind=sgl),INTENT(IN)      :: t(2)

self%HOLZvar%rectangle = t

end subroutine setrectangle_

!--------------------------------------------------------------------------
subroutine setlaL_(self, t)
!DEC$ ATTRIBUTES DLLEXPORT :: setlaL_
  !! author: MDG
  !! version: 1.0
  !! date: 01/28/20
  !!
  !! set laL

IMPLICIT NONE

class(HOLZ_T),INTENT(INOUT)    :: self
real(kind=sgl),INTENT(IN)      :: t

self%HOLZvar%laL = t

end subroutine setlaL_

!--------------------------------------------------------------------------
subroutine setImax_(self, t)
!DEC$ ATTRIBUTES DLLEXPORT :: setImax_
  !! author: MDG
  !! version: 1.0
  !! date: 01/28/20
  !!
  !! set Imax

IMPLICIT NONE

class(HOLZ_T),INTENT(INOUT)    :: self
real(kind=sgl),INTENT(IN)      :: t

self%HOLZvar%Imax = t

end subroutine setImax_

!--------------------------------------------------------------------------
subroutine setGmax_(self, t)
!DEC$ ATTRIBUTES DLLEXPORT :: setGmax_
  !! author: MDG
  !! version: 1.0
  !! date: 01/28/20
  !!
  !! set Gmax

IMPLICIT NONE

class(HOLZ_T),INTENT(INOUT)    :: self
real(kind=sgl),INTENT(IN)      :: t

self%HOLZvar%Gmax = t

end subroutine setGmax_

! !--------------------------------------------------------------------------
! subroutine set(self, t)
!   !! author: MDG
!   !! version: 1.0
!   !! date: 01/28/20
!   !!
!   !! return

! IMPLICIT NONE

! class(HOLZ_T),INTENT(INOUT)    :: self
! real(kind=sgl),INTENT(IN)      :: t

! self%HOLZvar% = t

! end subroutine set

!--------------------------------------------------------------------------
subroutine HOLZPage_(self, cell, SG, PS, Diff, camlen)
!DEC$ ATTRIBUTES DLLEXPORT :: HOLZPage_
  !! author: MDG
  !! version: 1.0
  !! date: 01/28/20
  !!
  !! draw zone axis HOLZ diffraction pattern and line pattern

use mod_crystallography
use mod_symmetry
use mod_io
use mod_diffraction
use mod_misc
use mod_postscript

IMPLICIT NONE

class(HOLZ_T),INTENT(INOUT)             :: self
type(Cell_T),INTENT(INOUT)              :: cell
type(SpaceGroup_T),INTENT(INOUT)        :: SG
type(PostScript_T),INTENT(INOUT)        :: PS
type(Diffraction_T),INTENT(INOUT)       :: Diff
real(kind=sgl),INTENT(IN)               :: camlen

type(Cell_T)                            :: savecell
type(IO_T)                              :: Message

type(HOLZreflection),pointer            :: temp, bot
type(HOLZvartype)                       :: HOLZv
type(HOLZentries)                       :: HOLZe

logical                                 :: again, first, newzone, nexttop
real(kind=sgl)                          :: negative(2),twopi,thr, oi_real(9), rect(2), &
                                           RHOLZmax,RHOLZ(20),xo,yo,sc,pos(2),dy
real(kind=sgl),parameter                :: xoff(0:5)=(/0.0,3.3125,0.0,3.3125,0.0,3.3125/), &
                                           yoff(0:5)=(/6.0,6.0,3.0,3.0,0.0,0.0/), &
                                           eps = 1.0E-3
integer(kind=irg)                       :: g1i(3),g2i(3),i,numHOLZ, io_int(1), psunit
real(kind=sgl),parameter                :: le=3.25,he=2.9375

! Historical note:
! this program was originally written for comparison with a TEM negative, hence
! the  strange pattern dimensions ... This may need to be redone in some other way at some point...

 associate(HOLZvar => self%HOLZvar)

 savecell = cell
 psunit = PS%getpsunit()

! dimensions of a standard TEM negative in inches:  3.9375 x 3.1875
 negative = (/ 3.9375, 3.1875 /)
 call self%setrectangle(negative*0.5)
 thr = 1.E-4
 twopi = 2.0*cPi
 call self%setImax(0.0)

! camera length
 call self%setlaL(sngl(Diff%getWaveLength()) * camlen)

! and print some information
 oi_real(1) = sngl(Diff%getWaveLength())
 call Message%WriteValue(' wavelength [nm] = ', oi_real, 1, "(F10.6)")
 oi_real(1) = camlen
 call Message%WriteValue(' L         [mm] = ', oi_real, 1, "(f10.2)")
 oi_real(1) = self%getlaL()
 call Message%WriteValue(' camera length lambda*L [mm nm] = ', oi_real, 1, "(f10.5)")

! what portion of reciprocal space is covered by the negative along horizontal direction ?
 rect = self%getrectangle()
 call self%setGmax(sqrt(rect(1)**2+rect(2)**2)*25.4/self%getlaL())

! this is also the maximum allowable HOLZ radius
 RHOLZmax = self%getGmax()

! next loop over multiple zone axis orientations or different lattice parameters
 again = .TRUE.
 first = .TRUE.
 do while (again)

! new zone axis or modify lattice parameters ?
  if (first.eqv..TRUE.) then
   newzone = .TRUE.

! get the zone axis
   call GetIndex(SG%getSpaceGrouphexset(),HOLZvar%uvw,'d')
   call Message%ReadValue(' Enter foil thickness [nm, R] (to get relrods of proper length) : ', oi_real, 1)
   HOLZvar%thickness = oi_real(1)
   first = .FALSE.

! get the linked list to store all reflections
   call self%getListHead( bot )
   if (.not.associated(bot)) then
    nullify(bot%next)
   end if
 else
! either get a new zone or change the lattice parameters and keep the zone
   call Message%ReadValue(' New zone (1) or change lattice parameters for present zone (2) ', io_int, 1)
   if (io_int(1).eq.1) then
     newzone = .TRUE.
     call GetIndex(SG%getSpaceGrouphexset(),HOLZvar%uvw,'d')
     cell = savecell
! deallocate the previous linked list and allocate a new one
     call self%DeleteList()
     call self%getListHead( bot )
     nullify(bot%next)
   else
! show the current lattice parameters and save the parameters
     newzone = .FALSE.
   end if
  end if

! it is a new zone, so draw the diffraction pattern on the top half of the page
  if (newzone.eqv..TRUE.) then
! get the basis vectors g1 and g2
    call cell%ShortestG(SG, HOLZvar%uvw,g1i,g2i,i)
    HOLZvar%g1 = float(g1i); HOLZvar%g2 = float(g2i)

! get the beam divergence angle to determine the diameter of the central disk
    oi_real(1) = 500.0*minval( (/ Diff%CalcDiffAngle(cell,(/ g1i(1),g1i(2),g1i(3) /)), &
                                  Diff%CalcDiffAngle(cell,(/ g2i(1),g2i(2),g2i(3) /)) /) )
    call Message%WriteValue(' Maximum disk diameter without overlap [mrad]= ', oi_real, 1, "(f10.4)")
    call Message%ReadValue(' Enter the beam divergence angle [mrad, R] ', oi_real, 1);
    HOLZvar%thetac = oi_real(1)*0.001

! distance between consecutive HOLZ layers in nm-1
    HOLZvar%H = 1.0/cell%CalcLength(float(HOLZvar%uvw),'d')

! determine g3 basis vector
    call cell%CalcCross(HOLZvar%g1,HOLZvar%g2,HOLZvar%g3,'r','r',1)
    call cell%NormVec(HOLZvar%g3,'r')
    HOLZvar%g3 = HOLZvar%H * HOLZvar%g3

! get foil normal
    call Message%printMessage('Enter Foil Normal F [real space indices]',"(A)")
    call GetIndex(SG%getSpaceGrouphexset(),HOLZvar%FN,'d')
! compute components of FN with respect to g1, g2, g3
    call cell%TransSpace(float(HOLZvar%FN),HOLZvar%FNr,'d','r')
    call cell%NormVec(HOLZvar%FNr,'r')
    HOLZvar%FNg = (/ cell%CalcDot(HOLZvar%FNr,HOLZvar%g1,'r'), &
                          cell%CalcDot(HOLZvar%FNr,HOLZvar%g2,'r'), &
                          cell%CalcDot(HOLZvar%FNr,HOLZvar%g3,'r') /)

! determine shortest vector of FOLZ layer
    call self%ShortestGFOLZ(cell)
    oi_real(1:3) = HOLZvar%gp(1)*HOLZvar%g1(1:3)+HOLZvar%gp(2)*HOLZvar%g2(1:3)
    call Message%WriteValue(' HOLZ shift vector = ', oi_real, 3, "(3f9.4)")

! get Laue center in terms of g1 and g2
    oi_real(1:3) = HOLZvar%g1(1:3)
    oi_real(4:6) = HOLZvar%g2(1:3)
    oi_real(7:9) = HOLZvar%g3(1:3)
    call Message%WriteValue( 'The new basis vectors for this zone axis are ', &
          oi_real, 9, "(/'g1 = ',3f10.5,/'g2 = ',3f10.5,/'g3 = ',3f10.5,/)")
    oi_real(1) = HOLZvar%H
    call Message%WriteValue(' reciprocal interplanar spacing H = ', oi_real, 1, "(F10.4,' nm^-1'/)")
    call Message%ReadValue(' Enter the coordinates of the Laue center with respect to g1 and g2', oi_real, 2)
    HOLZvar%LC1 = oi_real(1)
    HOLZvar%LC2 = oi_real(2)

! compute how many HOLZ layers need to be drawn.
! this follows from the camera length and the size of the micrograph
    i=1
    numHOLZ = 0
    do while(i.lt.20)
     RHOLZ(i) = sqrt(2.0*HOLZvar%H*float(i)/Diff%getWaveLength() - (float(i)*HOLZvar%H)**2)
     if (RHOLZ(i).lt.RHOLZmax) numHOLZ = numHOLZ+1
     i=i+1
    end do

! print the  number and radii of the HOLZ rings in the field of view
    call Message%printMessage('Number and radii of possible HOLZ rings inside field of view',"(A)")
    oi_real(1) = RHOLZmax
    call Message%WriteValue(' RHOLZmax = ', oi_real, 1, "(F10.5)")
    do i=1,numHOLZ
      oi_real(1)=float(i); oi_real(2)=RHOLZ(i)
      call Message%WriteValue(' Ring ', oi_real, 2, "(F5.0,3x,F10.5)")
    end do

! do page preamble stuff if this is a new page
! [This assumes that the PostScript file has already been opened]
    if (newzone) then
     call PS%newpage(.FALSE.,'Kinematical HOLZ Diffraction Patterns')
     call PS%text(5.25,-0.05,'scale bar in reciprocal nm')
     call PS%textvar(5.25,PS%getpsfigheight()+0.02,'Camera Constant [nm mm]',HOLZvar%laL)
     call PS%setfont(PSfonts(2),0.15)
     call PS%text(-0.25,PS%getpsfigheight()+0.02,'Structure File : '//trim(cell%getFileName()))
    end if

! draw frame and related stuff
    xo = 2.25
    yo = 5.00
    HOLZvar%PX = xo + HOLZvar%rectangle(1)
    HOLZvar%PY = yo + HOLZvar%rectangle(2)
    HOLZvar%CBEDrad = 1.5
    HOLZvar%CBEDsc = 1.3
    call PS%setlinewidth(0.012)
    call PS%balloon(xo,yo,negative(1),negative(2),0.0312)

! zone axis
    call PS%setfont(PSfonts(2),0.12)
    call PS%text(xo+0.05,yo+negative(2)+0.12,'Zone axis ')
    call PS%PrintIndices('d',SG%getSpaceGrouphexset(),HOLZvar%uvw(1),HOLZvar%uvw(2),HOLZvar%uvw(3), &
                      xo+0.6,yo+negative(2)+0.12)

! add other data lines to the upper left
    call PS%setfont(PSfonts(2),0.15)
    call PS%textvar(-0.25,PS%getpsfigheight()-0.18,'Acc. Voltage [kV] ',sngl(Diff%getV()))
    call PS%text(-0.25,PS%getpsfigheight()-0.38,'Foil normal ')
    call PS%PrintIndices('d',SG%getSpaceGrouphexset(),HOLZvar%FN(1),HOLZvar%FN(2),HOLZvar%FN(3),&
                      -0.25+1.5,PS%getpsfigheight()-0.38)
    call PS%textvar(-0.25,PS%getpsfigheight()-0.58,'Foil thickness [nm] ',HOLZvar%thickness)
    call PS%text(-0.25,PS%getpsfigheight()-0.78,'Laue center ')
    call PS%textvar(-0.25+1.1,PS%getpsfigheight()-0.78,'',HOLZvar%LC1)
    call PS%textvar(-0.25+1.6,PS%getpsfigheight()-0.78,'',HOLZvar%LC2)

! HOLZ ring radii
    call PS%text(xo-1.5,PS%getpsfigheight()-1.45,'HOLZ radii [nm-1] ')
    do i=1,numHOLZ
        call PS%textint(xo-1.5,PS%getpsfigheight()-1.5-float(i)*0.14,'',i)
        call PS%textvar(xo-1.3,PS%getpsfigheight()-1.5-float(i)*0.14,'',RHOLZ(i))
    end do

! CBED 000 disk text
    call PS%setfont(PSfonts(2),0.12)
    call PS%textvar(-0.25,0.5,'Convergence angle [mrad] ',HOLZvar%thetac*1000.0)

! lattice parameters
    call PS%setfont(PSfonts(4),0.14)
    call PS%text(-0.25,2.00,'a :')
    call PS%text(-0.25,1.84,'b :')
    call PS%text(-0.25,1.68,'c :')
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,2.00,cell%getLatParm('a')
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.84,cell%getLatParm('b')
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.68,cell%getLatParm('c')
    call PS%setfont(PSfonts(1),0.14)
    call PS%text(-0.25,1.52,'a :')
    call PS%text(-0.25,1.36,'b :')
    call PS%text(-0.25,1.20,'g :')
    call PS%setfont(PSfonts(4),0.14)
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.52,cell%getLatParm('alpha')
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.36,cell%getLatParm('beta ')
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.20,cell%getLatParm('gamma')

! scale bar (sc is the conversion factor from nm-1 to inches)
    sc = HOLZvar%laL/25.4
    call PS%setlinewidth(0.020)
    call PS%line(xo+0.05,yo+0.06,xo+0.05+5.0*sc,yo+0.06)
    call PS%setfont(PSfonts(2),0.15)
    call PS%text(xo+0.05+2.5*sc,yo+0.10,'5 ')

! plot origin of reciprocal space
    call PS%filledcircle(HOLZvar%PX,HOLZvar%PY,0.03,0.0)

! set clip path
    call PS%closepath
    call PS%gsave
    call PS%move(xo,yo)
    call PS%draw(xo,yo+negative(2))
    call PS%draw(xo+negative(1),yo+negative(2))
    call PS%draw(xo+negative(1),yo)
    call PS%clippath
    call PS%newpath

! now we are ready to compute the HOLZ layers and store them in the linked list
! first the ZOLZ layer
    call self%CalcHOLZ(cell,SG,Diff,PS,0)

! then the HOLZ layers
    do i=1,numHOLZ
      call PS%setlinewidth(0.005)
      call PS%circle(HOLZvar%PX,HOLZvar%PY,RHOLZ(i)*HOLZvar%laL/25.4)
      call self%CalcHOLZ(cell,SG,Diff,PS,i)
    end do

! once that is done, we can determine the intensities to be drawn and
! draw all reflections.
    call self%PlotHOLZ(PS)

! and eliminate the current clippath
    call PS%closepath
    call PS%grestore

! draw the vectors g1 and g2, and the projection gp of G
    call PS%setlinewidth(0.02)
    call PS%filledcircle(0.5,HOLZvar%PY-2.5,0.015,0.0)

! g1
    pos = matmul(HOLZvar%gtoc,(/1.0,0.0/) )
    HOLZvar%glen = 2.0*sqrt(pos(1)**2+pos(2)**2)
    pos = pos/HOLZvar%glen
    call PS%line(0.5,HOLZvar%PY-2.5,0.5+pos(1),HOLZvar%PY-2.5+pos(2))
    call PS%filledcircle(0.5+pos(1),HOLZvar%PY-2.5+pos(2),0.04,0.0)
    call PS%PrintIndices('r',SG%getSpaceGrouphexset(),int(HOLZvar%g1(1)),int(HOLZvar%g1(2)),int(HOLZvar%g1(3)), &
                      0.5+pos(1)+0.1, HOLZvar%PY-2.5+pos(2))

! g2
    pos = matmul(HOLZvar%gtoc,(/0.0,1.0/) )
    pos = pos/HOLZvar%glen
    call PS%line(0.5,HOLZvar%PY-2.5,0.5+pos(1),HOLZvar%PY-2.5+pos(2))
    call PS%filledcircle(0.5+pos(1),HOLZvar%PY-2.5+pos(2),0.04,0.0)
    call PS%PrintIndices('r',SG%getSpaceGrouphexset(),int(HOLZvar%g2(1)),int(HOLZvar%g2(2)),int(HOLZvar%g2(3)), &
                           0.5+pos(1)+0.1, HOLZvar%PY-2.5+pos(2))

! and then the projection of G onto g1,g2
    pos = matmul(HOLZvar%gtoc,(/HOLZvar%gp(1),HOLZvar%gp(2)/) )
    pos = pos/HOLZvar%glen
    call PS%setlinewidth(0.02)
    call PS%line(0.45+pos(1),HOLZvar%PY-2.5+pos(2),0.55+pos(1),HOLZvar%PY-2.5+pos(2))
    call PS%line(0.5+pos(1),HOLZvar%PY-2.45+pos(2),0.5+pos(1),HOLZvar%PY-2.55+pos(2))
    call PS%text(-0.5,3.2,'Basis vectors g1, g2,')
    call PS%text(-0.5,3.0,'and projection of G (cross)')
! draw fixed radius circle for bright field CBED disk
    call PS%setlinewidth(0.025)
    call PS%circle(HOLZvar%PX,yo-2.5,HOLZvar%CBEDrad)
! indicate center of pattern
    call PS%setlinewidth(0.01)
    call PS%line(HOLZvar%PX-0.05,yo-2.5,HOLZvar%PX+0.05,yo-2.5)
    call PS%line(HOLZvar%PX,yo-2.45,HOLZvar%PX,yo-2.55)

    call self%PlotHOLZlines(PS, 0.0)
    nexttop = .TRUE.
  else  ! this is not a new zone axis
! let the user define new lattice parameters
    savecell = cell
    oi_real(1) = cell%getLatParm('a'); oi_real(2) = cell%getLatParm('b'); oi_real(3) = cell%getLatParm('c')
    call Message%WriteValue(' Current lattice parameters [nm] ', oi_real, 3, "(/'a = ',f7.5,', b = ',f7.5,', c = ',f7.5)")
    oi_real(1) = cell%getLatParm('alpha'); oi_real(2) = cell%getLatParm('beta'); oi_real(3) = cell%getLatParm('gamma')
    call Message%WriteValue('', oi_real, 3, "(/'alpha = ',f7.2,', beta = ',f7.2,', gamma = ',f7.2)")

! ask for the new parameters (all must be entered) and recompute metric information
    call Message%ReadValue(' Enter new lattice parameters a, b, and c [nm] ', oi_real, 3)
    call cell%setLatParm('a', dble(oi_real(1)))
    call cell%setLatParm('b', dble(oi_real(2)))
    call cell%setLatParm('c', dble(oi_real(3)))
    call Message%ReadValue(' Enter new angles alpha, beta, and gamma [degrees] ', oi_real, 3)
    call cell%setLatParm('alpha', dble(oi_real(1)))
    call cell%setLatParm('beta', dble(oi_real(2)))
    call cell%setLatParm('gamma', dble(oi_real(3)))
    call cell%calcMatrices()

! redo the geometry
    call self%ReCalcHOLZ(cell, SG, Diff)

! move to top or bottom for next drawing ?
    if (nexttop.eqv..TRUE.) then
     call PS%newpage(.FALSE.,'Kinematical HOLZ Diffraction Patterns')
     dy = 4.25
     nexttop=.FALSE.
    else
     dy = -0.25
     nexttop=.TRUE.
    end if

! CBED 000 disk text
    call PS%setfont(PSfonts(2),0.12)
    call PS%textvar(-0.25,0.5+dy,'Convergence angle [mrad] ',HOLZvar%thetac*1000.0)

! lattice parameters
    call PS%setfont(PSfonts(4),0.14)
    call PS%text(-0.25,2.00+dy,'a :')
    call PS%text(-0.25,1.84+dy,'b :')
    call PS%text(-0.25,1.68+dy,'c :')
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,2.00+dy,cell%getLatParm('a')
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.84+dy,cell%getLatParm('b')
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.68+dy,cell%getLatParm('c')
    call PS%setfont(PSfonts(1),0.14)
    call PS%text(-0.25,1.52+dy,'a :')
    call PS%text(-0.25,1.36+dy,'b :')
    call PS%text(-0.25,1.20+dy,'g :')
    call PS%setfont(PSfonts(4),0.14)
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.52+dy,cell%getLatParm('alpha')
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.36+dy,cell%getLatParm('beta')
    write (psunit,"(1x,F12.7,' ',F12.7,' M (',F12.4,') show')") 0.0,1.20+dy,cell%getLatParm('gamma')

! draw fixed radius circle for bright field CBED disk
    call PS%setlinewidth(0.025)
    call PS%circle(HOLZvar%PX,yo-2.5+dy,HOLZvar%CBEDrad)

! indicate center of pattern
    call PS%setlinewidth(0.01)
    call PS%line(HOLZvar%PX-0.05,yo-2.5+dy,HOLZvar%PX+0.05,yo-2.5+dy)
    call PS%line(HOLZvar%PX,yo-2.45+dy,HOLZvar%PX,yo-2.55+dy)

    call self%PlotHOLZlines(PS, dy)
  end if ! newzone .eq. .TRUE.

  call Message%ReadValue(' Another pattern ? [1/0] ', io_int, 1)
  if (io_int(1).ne.1) again=.FALSE.
 end do  ! end of main loop

end associate

end subroutine HOLZPage_

!--------------------------------------------------------------------------
subroutine ShortestGFOLZ_(self, cell)
!DEC$ ATTRIBUTES DLLEXPORT :: ShortestGFOLZ_
 !! author: MDG
 !! version: 1.0
 !! date: 01/28/20
 !!
 !! find the shortest G vector

use mod_io
use mod_postscript
use mod_crystallography

IMPLICIT NONE

class(HOLZ_T), INTENT(INOUT)            :: self
type(cell_T),INTENT(INOUT)              :: cell

type(IO_T)                              :: Message
real(kind=sgl)                          :: gmin,gam11,gam12,gam22
integer(kind=irg),parameter             :: inm = 8
integer(kind=irg)                       :: ih,ik,il,NN, oi_int(1)

associate(HOLZvar => self%HOLZvar)

! look for the shortest reflection satisfying hu+kv+lw = 1
! This could be replaced by code from Jackson's paper (1987),
! but it does essentially the same thing.
 gmin = 100.0
 NN=1
 do while((gmin.eq.100.0).and.(NN.lt.4))
  do ih=-inm,inm
   do ik=-inm,inm
    do il=-inm,inm
! does this reflection lie in the plane NN ?
     if ((ih*HOLZvar%uvw(1)+ik*HOLZvar%uvw(2)+il*HOLZvar%uvw(3)).eq.NN) then
      HOLZvar%glen = cell%CalcLength(float((/ih,ik,il/)),'r')
      if (HOLZvar%glen.lt.gmin) then
       gmin = HOLZvar%glen
       HOLZvar%gshort = float( (/ ih,ik,il /) )
      end if
     end if
    end do
   end do
  end do
  oi_int(1) = NN
  call Message%WriteValue(' Could not find any reflections with hu+kv+lw = ', oi_int, 1, "(I2)")
  NN = NN+1
 end do
 if (gmin.eq.100.0) then ! for some reason there is no reflection with N<=3 ...
  call Message%printError('ShortestGFOLZ: ',' could not find any reflections with hu+kv+lw<=3 ...')
 end if

! projected components of G
 gam11 = cell%CalcDot(HOLZvar%g1,HOLZvar%g1,'r')
 gam12 = cell%CalcDot(HOLZvar%g1,HOLZvar%g2,'r')
 gam22 = cell%CalcDot(HOLZvar%g2,HOLZvar%g2,'r')
 gmin = 1.0/(gam11*gam22-gam12**2)
 HOLZvar%gp(1) = (cell%CalcDot(HOLZvar%gshort,HOLZvar%g1,'r')*gam22-&
                  cell%CalcDot(HOLZvar%gshort,HOLZvar%g2,'r')*gam12)*gmin
 HOLZvar%gp(2) = (cell%CalcDot(HOLZvar%gshort,HOLZvar%g2,'r')*gam11-&
                  cell%CalcDot(HOLZvar%gshort,HOLZvar%g1,'r')*gam12)*gmin

end associate

end subroutine ShortestGFOLZ_

!--------------------------------------------------------------------------
subroutine CalcHOLZ_(self, cell, SG, Diff, PS, N)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcHOLZ_
 !! author: MDG
 !! version: 1.0
 !! date: 01/28/20
 !!
 !! compute the HOLZ reflections for zone N

use mod_io
use mod_postscript
use mod_crystallography
use mod_symmetry
use mod_diffraction

IMPLICIT NONE

class(HOLZ_T), INTENT(INOUT)               :: self
type(Cell_T),INTENT(INOUT)                 :: cell
type(SpaceGroup_T),INTENT(INOUT)           :: SG
type(Diffraction_T),INTENT(INOUT)          :: Diff
type(PostScript_T),INTENT(INOUT)           :: PS
integer(kind=irg),INTENT(IN)               :: N

type(HOLZreflection),pointer               :: top, bot
type(IO_T)                                 :: Message
type(gnode)                                :: rlp
integer(kind=irg)                          :: inmhkl(2),i,j,nref,istat, oi_int(1)
real(kind=sgl)                             :: correction,gg(3),Ig,smax,gxy(2),pxy(2),exer,sgdenom,x,tgm,qx,qy,y,det,LC3, &
                                              ll(3),lpg(3),gplen
logical                                    :: a,dbdiff

associate(HOLZvar => self%HOLZvar)

 call Message%printMessage('Computing HOLZ reflection data',"(/A)")

 call self%getListHead( top )
 call self%getListHead( bot )

! set the index boundaries
 inmhkl(1) = int(1.1*HOLZvar%Gmax/cell%CalcLength(HOLZvar%g1,'r'))
 inmhkl(2) = int(1.1*HOLZvar%Gmax/cell%CalcLength(HOLZvar%g2,'r'))

! we will take g1 to be the x-axis of the pattern
! so the transformation matrix from g1,g2 to Cartesian is given by
 if (N.eq.0) then
  HOLZvar%phi = cell%CalcAngle(HOLZvar%g1,HOLZvar%g2,'r')
  HOLZvar%glen = cell%CalcLength(HOLZvar%g2,'r')
  HOLZvar%gtoc(1,1) = cell%CalcLength(HOLZvar%g1,'r')
  HOLZvar%gtoc(1,2) = HOLZvar%glen*cos(HOLZvar%phi)
  HOLZvar%gtoc(2,1) = 0.0
  HOLZvar%gtoc(2,2) = HOLZvar%glen*sin(HOLZvar%phi)
  HOLZvar%gtoc = HOLZvar%gtoc*HOLZvar%laL/25.4  ! nm^-1 to inches via the camera length
end if

! loop over all possible reflections in this layer (2D loop !!!)
 smax = 20.0/HOLZvar%thickness  ! the number 20 is arbitrary and could be changed at will
 nref=0
 do i=-inmhkl(1),inmhkl(1)
  do j=-inmhkl(2),inmhkl(2)
! make sure this reflection is close enough to the Ewald sphere to give some intensity;
! use crystal thickness to determine this according to the argument of the sinc function.
      gg = i*HOLZvar%g1 + j*HOLZvar%g2 + N*HOLZvar%gshort

! Do this only for those reflections that are allowed by
! the lattice centering !
      a = SG%IsGAllowed(int(gg))
      if (a) then

! compute excitation error, including Laue center, foil normal, and HOLZ reflection.
       HOLZvar%glen = cell%CalcLength(gg,'r')
       if (HOLZvar%glen.ne.0.0) then
         ll = HOLZvar%LC1*HOLZvar%g1 + HOLZvar%LC2*HOLZvar%g2
         lpg = ll + gg
         gplen = cell%CalcLength(lpg,'r')
         LC3 = sqrt(1.0-Diff%getWaveLength()**2*cell%CalcLength(ll,'r')**2)
         if (gplen.eq.0.0) then
           exer = -Diff%getWaveLength()*cell%CalcDot(gg,2.0*ll+gg,'r')/ &
                   (2.0*LC3*cos(cell%CalcAngle(float(HOLZvar%uvw),float(HOLZvar%FN),'d')) )
         else
           sgdenom = 2.0*LC3*cos(cell%CalcAngle(float(HOLZvar%uvw),float(HOLZvar%FN),'d'))- &
                   2.0*Diff%getWaveLength()*gplen*cos(cell%CalcAngle(lpg,HOLZvar%FNr,'r'))
           exer = -(Diff%getWaveLength()*cell%CalcDot(gg,2.0*ll+gg,'r')- &
                    2.0*LC3*gplen*cos(cell%CalcAngle(HOLZvar%g3,lpg,'r')))/sgdenom
         end if
       else
         exer = 10000.0
       end if

! exclude the 000 reflection
       if (abs(exer).le.smax) then

! OK, it is close enough.  Does it have any intensity ?
        call Diff%CalcUcg(cell,int(gg))
        rlp = Diff%getrlp()

! yes, it does.  get the scaled intensity using the sinc function
        Ig = rlp%Vmod**2
        if (Ig.lt.1.0e-16) then
          dbdiff = .TRUE.
        else
          dbdiff = .FALSE.
        end if
        if (abs(exer).ge.1.0e-5) then
         Ig = Ig * (sin(cPi*exer*HOLZvar%thickness)/(cPi*exer*HOLZvar%thickness))**2
        end if

! store maximum intensity
        if (Ig.gt.HOLZvar%Imax) HOLZvar%Imax = Ig

! next, determine the drawing coordinates, first in terms of g1 and g2
        correction = 1.0/(1.0-Diff%getWaveLength()*HOLZvar%H*(float(N)+exer*HOLZvar%FNg(3)))
        gxy = (/ (i+N*HOLZvar%gp(1)+exer*HOLZvar%FNg(1)), (j+N*HOLZvar%gp(2)+exer*HOLZvar%FNg(2))  /) * correction

! convert to Cartesian drawing coordinates
        pxy = matmul(HOLZvar%gtoc,gxy)

! and add the point to the linked list
        allocate(bot%next,stat=istat)
        if (istat.ne.0) then
          self%HOLZvar = HOLZvar
          call self%PlotHOLZ(PS)
          call PS%closefile()
          call Message%printError('CalcHOLZ: ',' unable to allocate memory for linked list')
        end if
        bot => bot%next
        nullify(bot%next)
        bot%hkl = gg
        bot%n1  = i
        bot%n2  = j
        bot%N   = N
        bot%sg  = exer
        bot%Ig  = Ig
        bot%pxy = pxy
        bot%dbdiff = dbdiff
        nref = nref+1

! would this point contribute to the HOLZ line drawing in the central disk ?
        HOLZvar%phi = asin(Diff%getWaveLength()*HOLZvar%glen*0.5) - asin(N*HOLZvar%H/HOLZvar%glen)
        bot%draw = .FALSE.
        if (abs(HOLZvar%phi).le.HOLZvar%thetac) then
         x = HOLZvar%phi/HOLZvar%thetac  *  HOLZvar%CBEDrad
         if (pxy(1).ne.0.0) then
           tgm = pxy(2)/pxy(1)
           y = atan2(pxy(2),pxy(1))
           qx = x*cos(y)
           qy = x*sin(y)
           det = 1.0-(1.0+tgm**2)*(1.0-(tgm*HOLZvar%CBEDrad*HOLZvar%CBEDsc/(qx+tgm*qy))**2)
           if (det.gt.0.0) then  ! there is an intersection for this line so it should be drawn
             bot%draw = .TRUE.
             bot%hlx(1) = (qx+tgm*qy)*(1.0-sqrt(det))/(1.0+tgm**2)
             bot%hly(1) = qy-(bot%hlx(1)-qx)/tgm
             bot%hlx(2) = (qx+tgm*qy)*(1.0+sqrt(det))/(1.0+tgm**2)
             bot%hly(2) = qy-(bot%hlx(2)-qx)/tgm
           end if
         else  ! parallel to the y-axis (easy to deal with)
             bot%draw = .TRUE.
             bot%hlx(1) = qx
             bot%hly(1) = sqrt((HOLZvar%CBEDrad*HOLZvar%CBEDsc)**2-qx**2)
             bot%hlx(2) = qx
             bot%hly(2) = -bot%hly(1)
         end if
        end if
       end if
    end if
   end do
  end do

  oi_int(1) = nref
  call Message%WriteValue(' number of reflections to be drawn : ', oi_int, 1, "(I6)")

end associate

end subroutine CalcHOLZ_

!--------------------------------------------------------------------------
subroutine ReCalcHOLZ_(self, cell, SG, Diff)
!DEC$ ATTRIBUTES DLLEXPORT :: ReCalcHOLZ_
 !! author: MDG
 !! version: 1.0
 !! date: 01/28/20
 !!
 !! compute the HOLZ reflections for zone N for new lattice parameters

use mod_crystallography
use mod_symmetry
use mod_diffraction
use mod_io

IMPLICIT NONE

class(HOLZ_T), INTENT(INOUT)            :: self
type(Cell_T),INTENT(INOUT)              :: cell
type(SpaceGroup_T),INTENT(INOUT)        :: SG
type(Diffraction_T),INTENT(INOUT)       :: Diff

type(IO_T)                              :: Message
type(HOLZreflection),pointer            :: temp

real(kind=sgl)                          :: correction,gg(3),gxy(2),pxy(2),exer,sgdenom,x,tgm,qx,qy, &
                                           y,det,LC3,ll(3),lpg(3),gplen

associate(HOLZvar => self%HOLZvar)

call Message%printMessage('Computing HOLZ reflection data',"(/A/)")

    call self%getListHead( temp )

    do while (associated(temp))
       gg = temp%hkl

! compute excitation error
       HOLZvar%glen = cell%CalcLength(gg,'r')
       ll = HOLZvar%LC1*HOLZvar%g1 + HOLZvar%LC2*HOLZvar%g2
       lpg = ll + gg
       gplen = cell%CalcLength(lpg,'r')
       LC3 = sqrt(1.0-Diff%getWaveLength()**2*cell%CalcLength(ll,'r')**2)
       if (gplen.eq.0.0) then
         exer = -Diff%getWaveLength()*cell%CalcDot(gg,2.0*ll+gg,'r')/2.0*LC3* &
                 cos(cell%CalcAngle(float(HOLZvar%uvw),float(HOLZvar%FN),'d'))
       else
         sgdenom = 2.0*LC3*cos(cell%CalcAngle(float(HOLZvar%uvw),float(HOLZvar%FN),'d'))- &
                 2.0*Diff%getWaveLength()*gplen*cos(cell%CalcAngle(lpg,HOLZvar%FNr,'r'))
         exer = -(Diff%getWaveLength()*cell%CalcDot(gg,2.0*ll+gg,'r')- &
                  2.0*LC3*gplen*cos(cell%CalcAngle(HOLZvar%g3,lpg,'r')))/sgdenom
       end if

! next, determine the drawing coordinates, first in terms of g1 and g2
       correction = 1.0/(1.0-Diff%getWaveLength()*HOLZvar%H*(float(temp%N)+exer*HOLZvar%FNg(3)))
       gxy = (/ (temp%n1+temp%N*HOLZvar%gp(1)+exer*HOLZvar%FNg(1)), &
                (temp%n2+temp%N*HOLZvar%gp(2)+exer*HOLZvar%FNg(2))  /) * correction

! convert to Cartesian drawing coordinates
       pxy = matmul(HOLZvar%gtoc,gxy)
       HOLZvar%phi = asin(Diff%getWaveLength()*HOLZvar%glen*0.5) - asin(temp%N*HOLZvar%H/HOLZvar%glen)
       temp % draw = .FALSE.
       if (abs(HOLZvar%phi).le.HOLZvar%thetac) then
        x = HOLZvar%phi/HOLZvar%thetac  *  HOLZvar%CBEDrad
        if (pxy(1).ne.0.0) then
         tgm = pxy(2)/pxy(1)
         y = atan2(pxy(2),pxy(1))
         qx = x*cos(y)
         qy = x*sin(y)
         det = 1.0-(1.0+tgm**2)*(1.0-(tgm*HOLZvar%CBEDrad*HOLZvar%CBEDsc/(qx+tgm*qy))**2)
         if (det.gt.0.0) then  ! there is an intersection for this line so it should be drawn
          temp%draw = .TRUE.
          temp%hlx(1) = (qx+tgm*qy)*(1.0-sqrt(det))/(1.0+tgm**2)
          temp%hly(1) = qy-(temp%hlx(1)-qx)/tgm
          temp%hlx(2) = (qx+tgm*qy)*(1.0+sqrt(det))/(1.0+tgm**2)
          temp%hly(2) = qy-(temp%hlx(2)-qx)/tgm
         end if
        else  ! parallel to the y-axis (easy to deal with)
          temp%draw = .TRUE.
          temp%hlx(1) = qx
          temp%hly(1) = sqrt((HOLZvar%CBEDrad*HOLZvar%CBEDsc)**2-qx**2)
          temp%hlx(2) = qx
          temp%hly(2) = -temp%hly(1)
        end if
       end if

! move to the next reflection
       temp=>temp%next
     end do

end associate

end subroutine ReCalcHOLZ_

!--------------------------------------------------------------------------
subroutine PlotHOLZ_(self, PS)
!DEC$ ATTRIBUTES DLLEXPORT :: PlotHOLZ_
 !! author: MDG
 !! version: 1.0
 !! date: 01/28/20
 !!
 !! draw a single HOLZ zone axis diffraction pattern

use mod_io
use mod_postscript

IMPLICIT NONE

class(HOLZ_T),INTENT(INOUT)              :: self
type(PostScript_T),INTENT(INOUT)         :: PS

type(HOLZreflection),pointer             :: temp
type(IO_T)                               :: Message

real(kind=sgl)                           :: V,qx,qy

associate(HOLZvar => self%HOLZvar)

 call Message%printMessage('Plotting HOLZ reflections',"(/A/)")

! point to the top of the linked list
 call self%getListHead( temp )

! move through the entire list and draw all reflections with exponentially scaled intensity;
 do while (associated(temp))
  V=0.05*(temp%Ig/HOLZvar%Imax)**0.1
  qx = temp%pxy(1)
  qy = temp%pxy(2)

! make sure the point is inside the rectangle
  if ((abs(qx).lt.HOLZvar%rectangle(1)).and.(abs(qy).lt.HOLZvar%rectangle(2))) then
   if (temp%dbdiff) then  ! potential double diffraction spot
    call PS%setlinewidth(0.005)
    call PS%square(HOLZvar%PX+qx,HOLZvar%PY+qy,0.04)
   else ! regular reflection
    call PS%filledcircle(HOLZvar%PX+qx,HOLZvar%PY+qy,V,0.0)
   end if
  end if

! move to the next reflection
  temp=>temp%next
 end do

end associate

end subroutine PlotHOLZ_

!--------------------------------------------------------------------------
subroutine PlotHOLZlines_(self, PS, dy)
!DEC$ ATTRIBUTES DLLEXPORT :: PlotHOLZlines_
 !! author: MDG
 !! version: 1.0
 !! date: 01/28/20
 !!
 !! draw a single HOLZ line convergent beam disk

use mod_io
use mod_postscript

IMPLICIT NONE

class(HOLZ_T),INTENT(INOUT)             :: self
type(PostScript_T),INTENT(INOUT)        :: PS
real(kind=sgl),INTENT(IN)               :: dy

type(HOLZreflection),pointer            :: temp

type(IO_T)                              :: Message
real(kind=sgl)                          :: V,V2,qx,qy,CB
character(12)                           :: txt
integer(kind=irg)                       :: i,nref, psunit


associate(HOLZvar => self%HOLZvar)

 psunit = PS%getpsunit()

 call Message%printMessage('Plotting HOLZ lines and labels',"(/A/)")

! point to the top of the linked list
 call self%getListHead(temp)
 nref = 0
 open(unit=30,file='temp.txt',status='unknown',form='formatted')
 call PS%setfont(PSfonts(4),0.08)

! move through the entire list and draw the corresponding HOLZ line across the 000 diffraction disk
 do while (associated(temp))

  if (temp%draw.eqv..TRUE.) then ! draw the HOLZ line on the second drawing
   V=0.015*(temp%Ig/HOLZvar%Imax)**0.1
   V2=1.0-1.0*(temp%Ig/HOLZvar%Imax)**0.1
   if (V2.lt.0.0) V2 = 0.0
   call PS%setlinewidth(V)

! make sure that all the lines will actually fit inside the region of interest
   if ((maxval(abs(temp%hlx)).le.HOLZvar%CBEDrad*HOLZvar%CBEDsc).and. &
       (maxval(abs(temp%hly)).le.HOLZvar%CBEDrad*HOLZvar%CBEDsc)) then
    call PS%line_gray(HOLZvar%PX+temp%hlx(1),2.5+temp%hly(1)+dy,HOLZvar%PX+temp%hlx(2),2.5+temp%hly(2)+dy,V2)

! add indices along continuation of lines
    qx= HOLZvar%PX+1.02*temp%hlx(2)
    qy= 2.5+1.02*temp%hly(2)+dy
    V = 180.0*atan2(temp%hly(2)-temp%hly(1),temp%hlx(2)-temp%hlx(1))/cPi
    write (30,"(1x,I3,1x,I3,1x,I3,1x,3(f10.5,1x))") (temp%hkl(i),i=1,3),qx,qy,V
    nref = nref+1
   end if
  end if

! move to the next reflection
  temp=>temp%next
 end do

! add indices along continuation of lines
 close(unit=30,status='keep')

 CB = HOLZvar%CBEDrad**2
 open(unit=30,file='temp.txt',status='old',form='formatted')
 do i=1,nref
   read (30,"(A12,1x,3(f10.5,1x))") txt,qx,qy,V
   if (((qx-HOLZvar%PX)**2+(qy-(2.5+dy))**2).gt.CB) then
    call PS%move(qx,qy)  ! just outside clipping ring
    write (psunit,*) V,' rotate'
    write (psunit,"(1x,'( ',A12,' ) show')") txt
    write (psunit,*) -V,' rotate'
   end if
 end do
 close(unit=30,status='delete')

end associate

end subroutine PlotHOLZlines_

!--------------------------------------------------------------------------
recursive function CalcsgHOLZ_(self,cell,HOLZdata,gg,kt,lambda) result(exer)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcsgHOLZ_
  !! author: MDG
  !! version: 1.0
  !! date: 01/26/22
  !!
  !! compute the excitation error including HOLZ and Laue Center information

use mod_crystallography

IMPLICIT NONE

class(HOLZ_T),INTENT(INOUT)             :: self
type(Cell_T),INTENT(INOUT)              :: cell
type(HOLZentries),INTENT(INOUT)         :: HOLZdata
!f2py intent(in,out) ::  HOLZdata
real(kind=sgl),INTENT(IN)               :: gg(3), kt(3), lambda

real(kind=sgl)                          :: exer, g1len, g2len
real(kind=sgl)                          :: ll(3), lpg(3), glen, gplen, LC1, LC2, LC3, sgdenom


glen = cell%CalcLength(gg,'r')
g1len = cell%CalcLength(HOLZdata%g1,'r')
g2len = cell%CalcLength(HOLZdata%g2,'r')
if (glen.ne.0.0) then
  LC1 = cell%CalcDot(kt,HOLZdata%g1,'r')/g1len
  LC2 = cell%CalcDot(kt,HOLZdata%g2,'r')/g2len
  ll = LC1*HOLZdata%g1 + LC2*HOLZdata%g2
  lpg = ll + gg
  gplen = cell%CalcLength(lpg,'r')
  LC3 = sqrt(1.0-lambda**2*cell%CalcLength(ll,'r')**2)
  if (gplen.eq.0.0) then
    exer = -lambda*cell%CalcDot(gg,2.0*ll+gg,'r')/(2.0*LC3*cell%CalcDot(HOLZdata%g3,HOLZdata%FNr,'r'))
  else
    sgdenom = 2.0*cell%CalcDot(LC3*HOLZdata%g3-lambda*lpg,HOLZdata%FNr,'r')
    exer = (cell%CalcDot(lpg,2.0*LC3*HOLZdata%g3-lambda*gg,'r')-lambda*cell%CalcDot(gg,ll,'r'))/sgdenom
  end if
else
  exer = 10000.0
end if

end function CalcsgHOLZ_

!--------------------------------------------------------------------------
recursive function GetHOLZcoordinates_(self,cell,HOLZdata,gg,kt,lambda) result(pxy)
!DEC$ ATTRIBUTES DLLEXPORT :: GetHOLZcoordinates_
  !! author: MDG
  !! version: 1.0
  !! date: 01/26/22
  !!
  !! find the projected coordinates of an arbitrary HOLZ g-vector

use mod_crystallography

IMPLICIT NONE

class(HOLZ_T)                           :: self
type(Cell_T),INTENT(INOUT)              :: cell          
type(HOLZentries),INTENT(INOUT)         :: HOLZdata
!f2py intent(in,out) ::  HOLZdata
real(kind=sgl),INTENT(IN)               :: gg(3), kt(3), lambda

real(kind=sgl)                          :: pxy(2), h1, h2, g11, g12, g22, z
real(kind=sgl)                          :: exer, correction, gxy(2), nx, ny, hh(3)
integer(kind=irg)                       :: N

! get the Laue zone number
  N = abs( HOLZdata%uvw(1)*gg(1) + HOLZdata%uvw(2)*gg(2) + HOLZdata%uvw(3)*gg(3) )

! get components of gg w.r.t. g1 and g2
  hh = gg - N * HOLZdata%gshort
  h1 = cell%CalcDot(hh,HOLZdata%g1,'c')
  h2 = cell%CalcDot(hh,HOLZdata%g2,'c')
  g11 = cell%CalcDot(HOLZdata%g1,HOLZdata%g1,'c')
  g12 = cell%CalcDot(HOLZdata%g1,HOLZdata%g2,'c')
  g22 = cell%CalcDot(HOLZdata%g2,HOLZdata%g2,'c')
  z = 1.0/(g12**2-g11*g22)
  nx = (g12*h2-g22*h1)*z
  ny = (g12*h1-g11*h2)*z

! compute excitation error, including Laue center, foil normal, and HOLZ reflection.
  exer = self%CalcsgHOLZ_(cell,HOLZdata,gg,kt,lambda)

! next, determine the drawing coordinates, first in terms of g1 and g2
  correction = 1.0/(1.0-lambda*HOLZdata%H*(float(N)+exer*HOLZdata%FNg(3)))
  gxy = (/ (nx+N*HOLZdata%gp(1)+exer*HOLZdata%FNg(1)), (ny+N*HOLZdata%gp(2)+exer*HOLZdata%FNg(2))  /) * correction

! convert to Cartesian drawing coordinates
  pxy = matmul(HOLZdata%gtoc,gxy)

end function GetHOLZcoordinates_

!--------------------------------------------------------------------------
recursive subroutine GetHOLZGeometry_(self,cell,HOLZdata,g1,g2,uvw,fn)
!DEC$ ATTRIBUTES DLLEXPORT :: GetHOLZGeometry_
  !! author: MDG
  !! version: 1.0
  !! date: 01/26/22
  !!
  !! initialize HOLZ geometrical data for a given zone axis

use mod_crystallography
use mod_io

IMPLICIT NONE

class(HOLZ_T),INTENT(INOUT)             :: self
type(Cell_T),INTENT(INOUT)              :: cell
type(HOLZentries),INTENT(INOUT)         :: HOLZdata
!f2py intent(in,out) ::  HOLZdata
integer(kind=irg),INTENT(IN)            :: uvw(3), fn(3)
real(kind=sgl),INTENT(IN)               :: g1(3), g2(3)

type(IO_T)                              :: Message 

real(kind=sgl)                          :: gmin,gam11,gam12,gam22, phi, glen, g3(3), c(3), gx(3), gy(3), gshort(3)
integer(kind=irg),parameter             :: inm = 8
integer(kind=irg)                       :: ih,ik,il,NN, oi_int(1)

! set some basic values
    HOLZdata%g1 = g1
    HOLZdata%g2 = g2
    HOLZdata%uvw = uvw
    HOLZdata%FN = fn

! distance between consecutive HOLZ layers in nm-1
    HOLZdata%H = 1.0/cell%CalcLength(float(uvw),'d')

! determine g3 basis vector
    call cell%CalcCross(HOLZdata%g1,HOLZdata%g2,g3,'r','r',1)
    call cell%NormVec(g3,'r')
    HOLZdata%g3 = HOLZdata%H * g3

! compute components of FN with respect to ga, gb, g3
    call cell%TransSpace(float(HOLZdata%FN),HOLZdata%FNr,'d','r')
    call cell%NormVec(HOLZdata%FNr,'r')
    HOLZdata%FNg = (/ cell%CalcDot(HOLZdata%FNr,HOLZdata%g1,'r'), cell%CalcDot(HOLZdata%FNr,HOLZdata%g2,'r'), &
                        cell%CalcDot(HOLZdata%FNr,g3,'r') /)

! look for the shortest reflection satisfying hu+kv+lw = 1
! This could be replaced by code from Jackson's paper (1987),
! but it does essentially the same thing.
 gmin = 100.0
 NN=1
 do while((gmin.eq.100.0).and.(NN.lt.4))
  do ih=-inm,inm
   do ik=-inm,inm
    do il=-inm,inm
! does this reflection lie in the plane NN ?
     if ((ih*uvw(1)+ik*uvw(2)+il*uvw(3)).eq.NN) then
      glen = cell%CalcLength(float((/ih,ik,il/)),'r')
      if (glen.lt.gmin) then
       gmin = glen
       gshort = float( (/ ih,ik,il /) )
      end if
     end if
    end do
   end do
  end do
  oi_int(1) = NN
  call Message%WriteValue(' Could not find any reflections with hu+kv+lw = ', oi_int, 1, frm = "(I2)")
  NN = NN+1
 end do
 if (gmin.eq.100.0) then ! for some reason there is no reflection with N<=3 ...
  call Message%printError('ShortestGFOLZ: ',' could not find any reflections with hu+kv+lw<=3 ...')
 end if
 HOLZdata%gshort = gshort

! projected components of G
 gam11 = cell%CalcDot(g1,g1,'r')
 gam12 = cell%CalcDot(g1,g2,'r')
 gam22 = cell%CalcDot(g2,g2,'r')
 gmin = 1.0/(gam11*gam22-gam12**2)
 HOLZdata%gp(1) = (cell%CalcDot(gshort,g1,'r')*gam22-cell%CalcDot(gshort,g2,'r')*gam12)*gmin
 HOLZdata%gp(2) = (cell%CalcDot(gshort,g2,'r')*gam11-cell%CalcDot(gshort,g1,'r')*gam12)*gmin

! coordinate transformation matrix for g1 along x (our standard orientation for all programs)
 phi = cell%CalcAngle(g1,g2,'r')
 glen = cell%CalcLength(g2,'r')
 HOLZdata%gtoc(1,1) = cell%CalcLength(g1,'r')
 HOLZdata%gtoc(1,2) = glen * cos(phi)
 HOLZdata%gtoc(2,1) = 0.0
 HOLZdata%gtoc(2,2) = glen * sin(phi)

! first normalize the zone axis in cartesian components; this is the z-axis
  call cell%TransSpace(float(uvw),c,'d','c')
  call cell%NormVec(c,'c')

! then make ga the x-axis
  call cell%TransSpace(g1,gx,'r','c')
  call cell%NormVec(gx,'c')
  HOLZdata%gx = gx

! compute the cross product between k and gx; this is the y-axis
  call cell%CalcCross(c,gx,gy,'c','c',0)
  HOLZdata%gy = gy


end subroutine GetHOLZGeometry_


end module mod_HOLZ
