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

program EMorbit
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/14/20
  !!
  !! print the orbit of a point (inside unit cell)

use mod_EMsoft
use mod_global
use mod_io
use mod_symmetry
use mod_crystallography

IMPLICIT NONE

type(EMsoft_T)              :: EMsoft 
type(IO_T)                  :: Message
type(Cell_T)                :: cell
type(SpaceGroup_T)          :: SG

character(fnlen)            :: progname = 'EMorbit.f90'
character(fnlen)            :: progdesc = 'List the orbit of a given position'

real(kind=dbl),allocatable  :: ctmp(:,:)
integer(kind=irg)           :: i,m,n,ans, io_int(1) 
real(kind=sgl)              :: io_real(3)
character(fnlen)            :: xtalname

 EMsoft = EMsoft_T(progname, progdesc, tpl = (/ 918 /) )
 
 call Message%ReadValue(' Enter xtal file name : ', xtalname,"(A)")
 call cell%getCrystalData(xtalname, SG, EMsoft)

 ans = 1
 do while (ans.eq.1)
  call Message%ReadValue(' Enter fractional coordinates of atom : ', io_real, 3)
  call SG%CalcOrbit(dble(io_real),n,ctmp)
  io_int(1) = n
  call Message%WriteValue('   # equivalent atoms in orbit        : ', io_int, 1, "(I3)")

! spit them out 
  do i=1,n
   io_int(1) = i
   io_real(1:3) = ctmp(i,1:3)
   call Message%WriteValue('', io_int, 1, "(1x,i3,' -> (')",advance="no")
   call Message%WriteValue('', io_real, 3, "(f7.4,',',f7.4,',',f7.4,');   ')",advance="no") 
   if (mod(i+1,2).eq.1) then
     call Message%printMessage(' ', frm = "(A)")
   end if
  end do
  call Message%ReadValue(' Another orbit ? (1/0) ', io_int, 1)
  ans = io_int(1)
  deallocate(ctmp)
 end do 

end program EMorbit
