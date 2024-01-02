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

program EMlistSG
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/14/20
  !!
  !! list general equivalent positions for a 
  !! given space group.  This is a simple program
  !! to illustrate how one can use the space group matrices.
  !!
  !! updated on 11/21/22 to include Hall space group symbols

use mod_EMsoft
use mod_global
use mod_io
use mod_symmetry
use mod_crystallography
use mod_HallSG

IMPLICIT NONE

type(EMsoft_T)              :: EMsoft 
type(IO_T)                  :: Message
type(Cell_T)                :: cell
type(SpaceGroup_T)          :: SG
      
character(fnlen)            :: progname = 'EMlistSG.f90'
character(fnlen)            :: progdesc = 'List equivalent positions for arbitrary space group'
      
character(3)                :: pos
integer(kind=irg)           :: p(4),ii,jj,i, io_int(1), HSGnum
real(kind=sgl)              :: ppp, io_real(1)
real(kind=dbl),allocatable  :: sgdata(:,:,:)
character(fnlen)            :: mess
character(16)               :: HS

 EMsoft = EMsoft_T(progname, progdesc, tpl = (/ 915 /) )

 call Message%printMessage('This program can use the standard space groups as well as') 
 call Message%printMessage('the Hall space group symbols; enter a negative space group') 
 call Message%printMessage('number to display the Hall symbols for the desired space group')
 call Message%printMessage('')
 call Message%ReadValue(' Enter Space Group number : ', io_int, 1) 
 if (io_int(1).lt.0) then 
   HS = List_Hall_Symbols(abs(io_int(1)),HSGnum) 
   SG = SpaceGroup_T( SGnumber = -io_int(1), useHall=.TRUE., HallSGnumber=HSGnum )
   call Message%printMessage('  Hall Space Group Symbol  : '//trim(HS), frm = "(A)")
 else
   SG = SpaceGroup_T( SGnumber = io_int(1) )
   call Message%printMessage('  Space Group Symbol       : '//trim(SG%getSpaceGroupName()), frm = "(A)")
 end if 

 pos = 'xyz'

 io_int(1) = SG%getSpaceGroupMATnum()
 call Message%WriteValue('  number of operators      : ', io_int, 1, "(I3)")

! get the symmetry matrices 
 sgdata = SG%getSpaceGroupDataMatrices()

! loop over all symmetry matrices
 do i=1, SG%getSpaceGroupMATnum()
  io_int(1) = i; 
  call Message%WriteValue(' ', io_int, 1,"(1x,i3,2x,'-> (')",advance="no")

! loop over all rows
  do ii=1,3

! loop over colums (get the numbers)
     p(1:3)=int(sgdata(i,ii,1:3)) 
     ppp = sngl(sgdata(i,ii,4))

! print each entry 
! first the part containing x, y, and/or z
     do jj=1,3
      if (p(jj).ne.0) then
       mess(1:1)='+'
       if (p(jj).eq.-1) mess(1:1)='-'
       mess(2:2) = pos(jj:jj)
       call Message%printMessage(mess, frm = "(A2)",advance="no")
      end if
     end do 

! if there is a translation component, print it
   if (ppp.ne.0.0) then
    mess(1:1)='+'
    if (ppp.lt.0.0) mess(1:1)='-'
    call Message%printMessage(mess, frm = "(A1)",advance="no");
    io_real(1) = abs(ppp); 
    call Message%WriteValue('', io_real, 1, "(f5.3)",advance="no")
   end if

! print a comma, or close the brackets and do a newline
   if (ii.ne.3) then 
     mess(1:1) = ','; call Message%printMessage(mess, frm = "(A1)",advance="no")
    else
     mess(1:1) = ')'; call Message%printMessage(mess, frm = "(A1)")
   end if

  end do
 end do 

end program EMlistSG
