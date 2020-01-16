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

module mod_misc
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! a collection of miscellaneous useful routines (no class definitions, just routines)

use mod_kinds
use mod_global

contains

!--------------------------------------------------------------------------
recursive subroutine IndexReduce(hkl)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! Reduce an index triplet to smallest integers 

IMPLICIT NONE

integer(kind=irg),INTENT(INOUT)      :: hkl(3)      
 !! indices

integer(kind=irg)                        :: mi,i,j
real(kind=sgl)                           :: rhkl(3),ir

 mi=100
 do i=1,3
  if ((abs(hkl(i)).lt.mi).and.(hkl(i).ne.0)) mi=abs(hkl(i))
 end do
 
! then check if this index is a common divider of the others
 j = 0
 do i=1,3
  rhkl(i) = float(hkl(i))/float(mi)
  ir = abs(rhkl(i))-float(int(abs(rhkl(i))))
  if (ir.eq.0.0) j=j+1
 end do
 
 if (j.eq.3) mi=1
 hkl = int(rhkl*mi)

end subroutine IndexReduce

!--------------------------------------------------------------------------
recursive subroutine IndexReduceMB(hkl)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! Reduce a Miller-Bravais index quartet to smallest integers

IMPLICIT NONE


integer(kind=irg),INTENT(INOUT)      :: hkl(4)      
 !! indices

integer(kind=irg)                        :: mi,i,j
real(kind=sgl)                           :: rhkl(4),ir

 mi=100
 do i=1,4
  if ((abs(hkl(i)).lt.mi).and.(hkl(i).ne.0)) mi=abs(hkl(i))
 end do

! then check if this index is a common divider of the others
 j = 0
 do i=1,4
  rhkl(i) = float(hkl(i))/float(mi)
  ir = abs(rhkl(i))-float(int(abs(rhkl(i))))
  if (ir.eq.0.0) j=j+1
 end do

 if (j.eq.4) mi=1
 hkl = int(rhkl*mi)

end subroutine

!--------------------------------------------------------------------------
recursive subroutine IndexString(hexset,st,hkl,sp)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! Return a string of indices for printing (only deals with indices up to 9)

IMPLICIT NONE

logical,INTENT(IN)                      :: hexset
character(12),INTENT(OUT)                   :: st         
 !! output string
integer(kind=irg),INTENT(INOUT)         :: hkl(3)   
 !! index triplet
character(1),INTENT(IN)                     :: sp         
 !! space character 'd' or 'r'

integer(kind=irg)                             :: l,hkil(4),i
character(1),parameter                  :: numbers(0:9) = (/'0','1','2','3','4','5','6','7','8','9'/)

 do l=1,12
  st(l:l) = ' '
 end do
 l=1
 if (sp.eq.'d') then
  st(l:l)='['
  l=l+1
 else
  st(l:l)='\'  ! '
  l=l+1
  st(l:l)='('
  l=l+1
 end if

 if (hexset.eqv..FALSE.) then 
  do i=1,3
   if (hkl(i).lt.0) then
    st(l:l)='-'
   else
    st(l:l)=' '
   end if
   l=l+1
   st(l:l)=numbers(abs(hkl(i)))
   l=l+1
  end do
 else
  if (sp.eq.'d') then 
    call MilBrav(hkl,hkil,'34')
    call IndexReduceMB(hkil)
  else
    hkil(1:2) = hkl(1:2)
    hkil(3) = -(hkl(1)+hkl(2))
    hkil(4) = hkl(3)
  end if
  do i=1,4
   if ((hkl(i).lt.0).and.(i.ne.3)) then
    st(l:l)='-'
   else
    st(l:l)=' '
   end if
   l=l+1
   if (i.eq.3) then 
    st(l:l)='.'
   else
    st(l:l)=numbers(abs(hkil(i)))
   end if
   l=l+1
  end do
 end if
 if (sp.eq.'d') then
  st(l:l)=']'
  l=l+1
  st(l:l)=' '
  l=l+1
  st(l:l)=' '
 else
  st(l:l)='\'  ! '
  l=l+1
  st(l:l)=')'
 end if

end subroutine IndexString

!--------------------------------------------------------------------------
recursive subroutine MilBrav(p,q,d)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/07/20
  !!
  !! conversion from Miller to Miller-Bravais indices for
  !! directions.  The switch d is either '34' or '43'.
  !! implements equations 1.31 and 1.32, pages 24-25. 

IMPLICIT NONE

integer(kind=irg),INTENT(INOUT)         :: p(3)         
 !! input/output vector
!f2py intent(in,out) ::  p
integer(kind=irg),INTENT(INOUT)         :: q(4)         
 !! input/output vector
!f2py intent(in,out) ::  q
character(2),INTENT(IN)                 :: d            
 !! direction string ('34' or '43')

integer(kind=irg)                       :: i, j         
real(kind=sgl)                          :: r(4), rm, tmp(4)     

 if (d.eq.'43') then 
! equation 1.31
! these will always be integers, so no reduction is required
  p(1) = q(1)-q(3)
  p(2) = q(2)-q(3)
  p(3) = q(4)
 else
! equation 1.32
! there is no need to divide by 3, since that would be taken out 
! by the reduction to integers in the next step
  r(1) = float(2*p(1)-p(2))
  r(2) = float(2*p(2)-p(1))
  r(3) = -float(p(1)+p(2))
  r(4) = float(3*p(3))

! next reduce to common integers
! first, find the non-zero minimum index
  rm = 100.0
  do i=1,4 
   if ((abs(r(i)).lt.rm).and.(r(i).ne.0.0)) rm = abs(r(i))
  end do

! then check if this index is a common divider of the others
  j = 0
  do i=1,4
   tmp(i) = r(i)/rm
   if ( ( abs(tmp(i))-int(abs(tmp(i))) ).eq.0.0) j=j+1
  end do
  if (j.eq.4) then
    q = tmp
  else  
    q = r
  end if
 end if

end subroutine MilBrav

!--------------------------------------------------------------------------
recursive subroutine ProjectionMatrix(cell,iview,M)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! construct a projection matrix for a given direct space viewing direction

use mod_crystallography

IMPLICIT NONE

type(Cell_T),INTENT(INOUT)      :: cell
real(kind=sgl),INTENT(OUT)      :: M(3,3)               
 !! output transformation matrix
integer(kind=irg),INTENT(IN)    :: iview(3)             
 !! input viewing direction indices

real(kind=sgl)                  :: g(3),r(3),q(3),qmin  ! auxiliary variables
integer(kind=irg)               :: i,imin               ! auxiliary variables

 g=float(iview)
 q=(/ 0.0,0.0,0.0 /)
 call cell%TransSpace(g,r,'d','c')
 call cell%NormVec(r,'c')

! the direction with the smallest direction cosine
! will be put parallel to the horizontal axis of the projection
 qmin=1.0
 do i=1,3
  if (abs(r(i)).lt.qmin) then
   qmin=abs(r(i))
   imin=i
  endif
 end do
 q(imin)=1.0

! cross rxq to get the y-direction.
 call cell%CalcCross(r,q,g,'c','c',0)
 call cell%NormVec(g,'c')

! cross gxr to get the x-direction.
 call cell%CalcCross(g,r,q,'c','c',0)
 call cell%NormVec(q,'c')

! fill the projection matrix
 M(1,1:3)=q(1:3)
 M(2,1:3)=g(1:3)
 M(3,1:3)=r(1:3)

end subroutine ProjectionMatrix

!--------------------------------------------------------------------------
recursive subroutine GetViewingDirection(hexset,iview)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! get the direct space indices of the viewing direction 

use mod_io

IMPLICIT NONE

logical,INTENT(IN)              :: hexset
integer(kind=irg),INTENT(OUT)   :: iview(3) 
 !! direction indices in three index notation

type(IO_T)                      :: Message
integer(kind=irg)               :: vview(4)     ! to transform between 3 and 4 index notation
integer(kind=irg)               :: io_int(4)    ! used for IO

! viewing direction (watch for hexagonal indices !)
 if (hexset.eqv..FALSE.) then
  call Message%ReadValue('Enter viewing direction indices [uvw] : ', io_int,3)
  iview(1:3) = io_int(1:3)
 else
  call Message%ReadValue('Enter viewing direction indices [uvtw] : ', io_int, 4)
  vview(1:4) = io_int(1:4)
  call MilBrav(iview,vview,'43')
 endif
 call IndexReduce(iview)

end subroutine GetViewingDirection

!--------------------------------------------------------------------------
recursive subroutine GetDrawingSpace(sp)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! real space or reciprocal space drawing

use mod_io

IMPLICIT NONE

character(1),INTENT(OUT)    :: sp        
 !! character 'd' or 'r'

type(IO_T)                  :: Message 

call Message%ReadValue('Real Space (d) or reciprocal space (r) : ', sp,'(A1)')

end subroutine GetDrawingSpace

!--------------------------------------------------------------------------
recursive subroutine GetIndex(hexset,ind,sp)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! get the u,v,w or h,k,l indices

use mod_io

IMPLICIT NONE

logical,INTENT(IN)            :: hexset
integer(kind=irg),INTENT(OUT) :: ind(3) 
 !! indices
character(1),INTENT(IN)       :: sp         
 !! space 'd' or 'r'

type(IO_T)                    :: Message
integer(kind=irg)               :: jnd(4)

 if (sp.eq.'d') then 
  if (hexset.eqv..FALSE.) then
   call Message%ReadValue('Enter u,v,w :', ind,3)
  else
   call Message%ReadValue('Enter u,v,t,w :',jnd,4)
   call MilBrav(ind,jnd,'43')
   call IndexReduce(ind)
  end if
 else
  call Message%ReadValue('Enter h,k,l :', ind,3)
 endif

end subroutine GetIndex





end module mod_misc