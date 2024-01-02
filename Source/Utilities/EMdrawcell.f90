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

program EMdrawcell 
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/30/19
  !!
  !! Draw a unit cell (or a few ...)
  !!
  !! This is not meant to be a very good drawing program; it is just a simple illustration
  !! of how one can use crystallographic information to make structure drawings.

use mod_kinds 
use mod_global
use mod_EMsoft
use mod_io
use mod_crystallography
use mod_symmetry 
use mod_postscript
use mod_math 
use mod_misc
use mod_HDFsupport

IMPLICIT NONE 

character(fnlen)               :: progname = 'EMdrawcell.f90'
character(fnlen)               :: progdesc='Draw one or more unit cells in perspective mode'

type(EMsoft_T)                 :: EMsoft
type(IO_T)                     :: Message 
type(Cell_T)                   :: cell 
type(SpaceGroup_T)             :: SG 
type(PostScript_T)             :: PS

real(kind=dbl),allocatable     :: apos(:,:,:)
integer(kind=irg),allocatable  :: atype(:), numat(:)

integer(kind=irg),parameter    :: n=1000
character(1)                   :: sp
integer(kind=irg)              :: acol(n)
        
real(kind=sgl)                 :: p(4),q(4),x(n),y(n),z(n),x1,y1,z1,asize(n),M(4,4),VD,diam,io_real(3), CX, CY, &
                                  shx, shy, xmid, ymid, xmin, xmax, ymin, ymax
integer(kind=irg)              :: i, j, icnt, ier, idx(n), iview(3), iform, io_int(3), imanum
character(fnlen)               :: gname, xtalname 

! program header and command line argument handling
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 909 /) )

call openFortranHDFInterface()

! read crystal information
call Message%ReadValue(' Enter xtal file name : ', xtalname,"(A)")
call cell%getCrystalData(xtalname, SG, EMsoft)

! reduce the coordinates 
call SG%setSpaceGroupreduce(.TRUE.)

! create all atoms
call cell%calcPositions(SG, 'm')

! dimension parameters
CX = 7.0
CY = 7.0

! real space drawing (reciprocal not implemented at this time)
sp = 'd'

! Viewing Distance = 3 times CX = 6 times xmax
call Message%ReadValue(' Viewing distance [nm] : ',io_real, 1)
VD = io_real(1)

! viewing direction
call GetViewingDirection( SG%getSpaceGrouphexset(), iview)

! create transformation matrix
M = ComputeViewTrans(cell, iview, VD)

! open PostScript file
imanum = 1
PS = PostScript_T( progdesc, EMsoft, imanum )

! write text and draw box
call PS%DrawcellFrame(cell, iview, sp, CX, CY, SG%getSpaceGrouphexset() )

! draw unit cell outline first
! which radii should be used for the drawing ?
call Message%ReadValue(' Use ionic radii (1) or metallic (2) :', io_int, 1)
iform = io_int(1)

! then get all atom coordinates
apos = cell%getapos()
numat = cell%getnumat()
atype = cell%getatomtype()

icnt=0
xmin=100.0
xmax=-100.0
ymin=100.0
ymax=-100.0
do i=1,cell%getNatomtype()
 do j=1,numat(i)
   p=(/sngl(apos(i,j,1)),sngl(apos(i,j,2)),sngl(apos(i,j,3)),1.0/)
   q = matmul(p,M)
   x1 = VD*q(1)/q(4) 
   y1 = VD*q(2)/q(4) 
   z1 = VD*q(3)/q(4) 
   x1 = VD*x1/z1
   y1 = VD*y1/z1
   icnt = icnt+1
   x(icnt) = 0.5*CX+2.5*x1
   y(icnt) = 0.5*CY+2.5*y1
   if (x(icnt).lt.xmin) xmin = x(icnt)
   if (x(icnt).gt.xmax) xmax = x(icnt)
   if (y(icnt).lt.ymin) ymin = y(icnt)
   if (y(icnt).gt.ymax) ymax = y(icnt)
   z(icnt) = z1
   if (iform.eq.1) then 
    asize(icnt) = ATOM_SPradii( atype(i) )
   else
    asize(icnt) = ATOM_MTradii( atype(i) )
   endif
   acol(icnt) = atype(i)
  end do
 end do

! shift the drawing back to the center of the page
 xmid = (xmax+xmin)*0.5
 ymid = (ymax+ymin)*0.5
 shx = 0.5*CX - xmid
 shy = 0.5*CY - ymid
 do i=1,icnt
  x(i) = x(i) + shx
  y(i) = y(i) + shy
 end do

! rank according to distance from observer
! draw atoms/relpoints in reverse order (farthest first)
 call SPSORT(z,icnt,idx,-1,ier)

 do i=1,icnt
  j = idx(i) 
  diam = asize(j)
  call PS%sphere(x(j), y(j), diam, acol(j))
 end do 

! close PostScript file
 call PS%closefile()

end program EMdrawcell

