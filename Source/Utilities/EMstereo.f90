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

program EMstereo
  !! author: MDG
  !! version: 1.0 
  !! date: 01/26/20
  !!
  !! Standard stereographic projections

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

IMPLICIT NONE

character(fnlen)      :: progname = 'EMstereo.f90'
character(fnlen)      :: progdesc = 'Stereographic projections (direct/ reciprocal space)'

type(EMsoft_T)        :: EMsoft
type(IO_T)            :: Message
type(Cell_T)          :: cell
type(SpaceGroup_T)    :: SG
type(PostScript_T)    :: PS

character(1)          :: sp
logical               :: topbot
integer(kind=irg)     :: iview(3), io_int(3), imanum
character(fnlen)      :: xtalname

! print header information and handle command line arguments 
 EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 925 /) )
 
! ask for the crystal structure file
call Message%ReadValue(' Enter xtal file name : ', xtalname,"(A)")
call cell%getCrystalData(xtalname, SG, EMsoft)

topbot=.FALSE.

! real space or reciprocal space?
call Message%ReadValue('Real Space (d) or reciprocal space (r) : ', sp,'(A1)')

! viewing direction (watch for hexagonal indices !)
call GetViewingDirection(SG%getSpaceGrouphexset(), iview)

! open PostScript file
 imanum = 1
 PS = PostScript_T(progdesc, EMsoft, imanum)

! get index ranges
 call Message%printMessage( (/ &
     '  Enter the maximum index for h,k and l, or for ', &
     '  u,v, and w. For a hexagonal system, please use', &
     '  4-index notation [uv.w] or (hk.l) to determine', &
     '  the largest index.                            ' /), frm = "(A/)")
 call Message%ReadValue(' Enter maximum indices (h,k,l) : ', io_int, 3)

! call the drawing routine
 call PS%StereoProj(cell, SG, sp, iview, io_int(1), io_int(2), io_int(3), topbot)

! close Postscript file
 call PS%closefile()

end program EMstereo

