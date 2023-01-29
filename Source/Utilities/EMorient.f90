! ###################################################################
! Copyright (c) 2013-2023, Marc De Graef Research Group/Carnegie Mellon University
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

program EMorient
  !! author: MDG
  !! version: 1.0 
  !! date: 01/26/20
  !!
  !! stereographic projection of a crystal orientation relation

use mod_kinds
use mod_global
use mod_EMsoft
use mod_symmetry
use mod_io
use mod_crystallography
use mod_rotations
use mod_quaternions
use mod_so3
use mod_misc 
use mod_postscript
use mod_HDFsupport

IMPLICIT NONE

character(fnlen)      :: progname = 'EMorient.f90'
character(fnlen)      :: progdesc = 'Stereographic projection of orientation relation'

type(EMsoft_T)        :: EMsoft
type(IO_T)            :: Message
type(Cell_T)          :: cellA, cellB 
type(SpaceGroup_T)    :: SGA, SGB 
type(PostScript_T)    :: PS 

character(1)              :: sp
logical                   :: nn,topbot
type(OrientationRelation) :: orel
real(kind=sgl)            :: rr(3),gg(3),g(3),r(3),M(3,3),negthresh,p(3),Ep(3,3),E(3,3),TT(3,3), io_real(3), &
                                 CX, CY, CRad, xst, yst
real(kind=dbl)            :: dE(3,3),dgg(3)
integer(kind=irg)         :: h,k,l,cr,hkl(3),iview(3),inm, i, ih, ik, il, imanum
character(fnlen)          :: xtalnameA, xtalnameB
character(17)         :: str
character(12)         :: instr


! print header information and handle command line arguments 
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 919 /) )
 
call openFortranHDFInterface()

! ask for the crystal structure file
call Message%ReadValue(' Enter xtal file name (phase A) : ', xtalnameA,"(A)")
call cellA%getCrystalData(xtalnameA, SGA, EMsoft)
call Message%ReadValue(' Enter xtal file name (phase B) : ', xtalnameB,"(A)")
call cellB%getCrystalData(xtalnameB, SGB, EMsoft)

inm=2
topbot=.FALSE.

! 20cm radius projection circle, centered on page [inches]
CRad = 3.937
CX = 3.25
CY = 3.5
negthresh=-0.0001

! get orientation relation
call GetOR(orel)
TT = ComputeOR(orel, cellA, cellB, 'AB')
 
! real space or reciprocal space?
call Message%ReadValue('Real Space (d) or reciprocal space (r) : ', sp,'(A1)')

! viewing direction (watch for hexagonal indices !)
call GetViewingDirection(SGA%getSpaceGrouphexset(), iview)

! create transformation matrix
call ProjectionMatrix(cellA,iview,M)

! open PostScript file
PS = PostScript_T(progdesc, EMsoft, imanum)

! write text and draw projection circle
call PS%newpage(.FALSE.,'Stereographic Projection')
call PS%setlinewidth(0.016)
call PS%circle(CX,CY,CRad)
call PS%setlinewidth(0.004)
call PS%line(CX-CRad,CY,CX+CRad,CY)
call PS%line(CX,CY-CRad,CX,CY+CRad)
call PS%setfont(PSfonts(2),0.08)
call PS%text(CX-CRad-0.07,CY-0.025,'A')
call PS%text(CX+CRad+0.03,CY-0.025,'B')
call PS%text(CX-0.03,CY-CRad-0.09,'M''')
call PS%text(CX-0.03,CY+CRad+0.07,'M"')
call PS%setfont(PSfonts(2),0.12/PS%getpsscale())
call PS%text(0.35,8.30,'Crystal A : '//trim(cellA%getFileName()))
call PS%filledcircle(0.0,8.30,0.015,0.0)
call PS%DumpIndices(SGA%getSpaceGrouphexset(),sp,0,0,0,1,0.0,8.30,.TRUE.)
call PS%setfont(PSfonts(2),0.12)
call PS%text(0.35,8.10,'Crystal B : '//trim(cellB%getFileName()))
call PS%filledsquare(0.0,8.10,0.035,0.0)
call PS%DumpIndices(SGB%getSpaceGrouphexset(),sp,0,0,0,2,0.0,8.10,.TRUE.)
call PS%setfont(PSfonts(2),0.12)
call IndexString(SGA%getSpaceGrouphexset(),instr,iview,'d')
call PS%text(0.0,7.90,'Viewing Direction '//instr//' [A]')

if (sp.eq.'d') then 
str='direct space'
else
str='reciprocal space'
endif
call PS%text(0.0,7.70,'Projection of '//str)

call PS%text(CX,8.20,'Orientation Relation ')
hkl(1:3)=int(orel % gA(1:3))

call IndexString(SGA%getSpaceGrouphexset(),instr,hkl,'r')
call PS%text(CX,8.00,'\(hkl\) : ')
call PS%text(CX+0.4,8.00,'A-'//instr)
hkl(1:3)=int(orel % gB(1:3))

call IndexString(SGB%getSpaceGrouphexset(),instr,hkl,'r')
call PS%text(CX+0.9,8.00,'|| B-'//instr)

! Space=.True.
hkl(1:3)=int(orel % tA(1:3))
call IndexString(SGA%getSpaceGrouphexset(),instr,hkl,'d')
call PS%text(CX,7.80,'[uvw] : ')
call PS%text(CX+0.4,7.80,'A-'//instr)
hkl(1:3)=int(orel % tB(1:3))
call IndexString(SGB%getSpaceGrouphexset(),instr,hkl,'d')
call PS%text(CX+0.9,7.80,'|| B-'//instr)

! loop over all planes or directions
 do h=-inm,inm
  do k=-inm,inm
   do l=-inm,inm
    ih=h
    ik=k
    il=l

! skip the origin
    if ((ih**2+ik**2+il**2).ne.0) then

! reduce to smallest integers to avoid overlap
! of indices, such as (111) and (222)
     hkl(1)=ih
     hkl(2)=ik
     hkl(3)=il
     call IndexReduce(hkl)

! transform to cartesian coordinates
     g(1)=float(hkl(1))
     g(2)=float(hkl(2))
     g(3)=float(hkl(3))
     ih = hkl(1)
     ik = hkl(2)
     il = hkl(3)

! crystal A
     call cellA%TransSpace(g,r,sp,'c')
     call cellA%NormVec(r,'c')

! apply viewing tansformation
     rr = matmul(M,r)

! compute stereographic projection coordinates
     xst=CX+CRad*rr(1)/(1.0+abs(rr(3)))
     yst=CY+CRad*rr(2)/(1.0+abs(rr(3)))
     cr=1
     if (rr(3).gt.negthresh) then
      call PS%filledcircle(xst,yst,0.015/PS%getpsscale(),0.0)
      nn = .TRUE.
      call PS%DumpIndices(SGA%getSpaceGrouphexset(),sp,ih,ik,il,cr,xst,yst,nn)
     else if (topbot) then
      call PS%circle(xst,yst,0.035/PS%getpsscale())
      nn = .FALSE.
      call PS%DumpIndices(SGA%getSpaceGrouphexset(),sp,ih,ik,il,cr,xst,yst,nn)
     end if

! crystal B
!    call TransCoor(cellB,dble(matmul(cellB % rsm,g)),dgg,dble(TT),sp,'on')
     call cellB%TransSpace(g,r,sp,'c')
     call cellB%TransCoor(dble(r),dgg,dble(TT),sp,'on')
     gg = sngl(dgg)
     call cellB%NormVec(gg,'c')

! apply viewing tansformation
     rr = matmul(M,gg)

! compute stereographic projection coordinates
     xst=CX+CRad*rr(1)/(1.0+abs(rr(3)))
     yst=CY+CRad*rr(2)/(1.0+abs(rr(3)))
     cr=2
     if (rr(3).gt.negthresh) then
      call PS%filledsquare(xst,yst,0.035/PS%getpsscale(),0.0)
      nn = .TRUE.
      call PS%DumpIndices(SGB%getSpaceGrouphexset(),sp,ih,ik,il,cr,xst,yst,nn)
     else if (topbot) then
      call PS%square(xst,yst,0.050/PS%getpsscale())
      nn = .FALSE.
      call PS%DumpIndices(SGB%getSpaceGrouphexset(),sp,ih,ik,il,cr,xst,yst,nn)
     end if
    end if
   end do 
  end do 
 end do 

! close Postscript file
 call PS%closefile()

end program EMorient
