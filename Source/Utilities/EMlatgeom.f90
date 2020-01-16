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

!--------------------------------------------------------------------------
! EMsoft:EMlatgeom.f90
!--------------------------------------------------------------------------
!
! PROGRAM: EMlatgeom 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief general lattice geometry computations
!
!> @details  not meant to be an all encompassing program; simply an illustration of how
!> to do crystallographic computations.
! 
!> @date 01/05/99 MDG 1.0 original
!> @date 05/21/01 MDG 2.0 f90
!> @date 04/16/13 MDG 3.0 rewrite
!> @date 06/13/14 MDG 4.0 rewrite without global variables
!--------------------------------------------------------------------------
program EMlatgeom
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! general lattice geometry computations
  !!
  !! @details  not meant to be an all encompassing program; simply an illustration of how
  !! to do crystallographic computations.

use mod_kinds
use mod_global
use mod_crystallography
use mod_io
use mod_symmetry
use mod_EMsoft

type(EMsoft_T)            :: EMsoft 
type(Cell_T)              :: cell 
type(SpaceGroup_T)        :: SG
type(IO_T)                :: Message 

character(fnlen)          :: progname = 'EMlatgeom.f90'
character(fnlen)          :: progdesc = 'Simple lattice geometry program'

integer(kind=irg)        	:: isel,another, oi_int(3)
real(kind=sgl)           	:: v1(3),v2(3),vc(3),p,q,r, oi_real(1)
character(1)             	:: sp,sp2
character(fnlen)		      :: xtalname

 EMsoft = EMsoft_T(progname, progdesc, tpl = (/ 914 /) )

! read crystal information
 call Message%ReadValue(' Enter xtal file name : ', xtalname,"(A)")
 call cell%getCrystalData(xtalname, SG, EMsoft)

! set up loop
 another=1
 do while (another.eq.1) 
  call PrintMenu(isel)
  sp='d'
  if (mod(isel,2).eq.0) sp='r'
  if (isel.lt.3) then
   call Message%ReadValue('    Enter vector components : ', v1, 3)
   oi_real(1) = cell%CalcLength(v1,sp)
   if (sp.eq.'d') then 
    call Message%WriteValue(' -> Length [nm] = ', oi_real, 1, "(2x,F10.6)")
   else
    call Message%WriteValue(' -> Length [nm-1] = ', oi_real, 1, "(2x,F10.6)")
   end if
  else
   call Message%ReadValue('    Enter first vector components : ', v1, 3)
   call Message%ReadValue('    Enter second vector components : ', v2, 3)
   if (isel.lt.5) then 
    p = cell%CalcLength(v1,sp)
    q = cell%CalcLength(v2,sp)
    r = cell%CalcDot(v1,v2,sp)
    oi_real(1) = cell%CalcAngle(v1,v2,sp)*180.0/cPi
    call Message%WriteValue(' -> Angle [deg] = ', oi_real, 1, "(2x,F8.4)")
   else
    if (sp.eq.'d') sp2='r'
    if (sp.eq.'r') sp2='d'
    call cell%CalcCross(v1,v2,vc,sp,sp2,0)
    oi_int(1:3)=int(vc(1:3))
    if (sp.eq.'d') then
     call Message%WriteValue(' ', oi_int, 3, "(' -> p x q = (',2(i3,1x),i3,')')")
    else
     call Message%WriteValue(' ', oi_int, 3, "(' -> p x q = [',2(i3,1x),i3,']')")
    end if
   end if
  end if
  call Message%ReadValue(' Another computation? (1/0) : ', oi_int, 1)
  another = oi_int(1)
 end do 
 
end program EMlatgeom

!--------------------------------------------------------------------------
subroutine PrintMenu(isel)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! Print a menu and return a selection

use mod_kinds
use mod_io
        
type(IO_T)                      :: Message
integer(kind=irg),INTENT(OUT)  	:: isel	! selection 
integer(kind=irg)		            :: io_int(1)

 call Message%printMessage(' Select from the following options: ', frm = "(A/)")
 call Message%printMessage(' [1] length of direct space vector', frm = "(A)")
 call Message%printMessage(' [2] length of reciprocal space vector', frm = "(A)")
 call Message%printMessage(' [3] angle between direct space vectors', frm = "(A)")
 call Message%printMessage(' [4] angle between reciprocal space vectors', frm = "(A)")
 call Message%printMessage(' [5] cross product, direct space vectors', frm = "(A)")
 call Message%printMessage(' [6] cross product, reciprocal space vectors', frm = "(A/)")
 call Message%ReadValue('    Enter selection: ', io_int, 1)
 isel = io_int(1)

end subroutine PrintMenu
