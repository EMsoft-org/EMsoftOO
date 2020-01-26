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

program EMZAgeom 
  !! author: MDG
  !! version: 1.0 
  !! date: 01/26/20
  !!
  !! EMZAgeom computes zone axis geometry & symmetry

use mod_kinds
use mod_global
use mod_EMsoft
use mod_crystallography
use mod_symmetry
use mod_io

IMPLICIT NONE

character(fnlen)        :: progname = 'EMZAgeom.f90'
character(fnlen)        :: progdesc = 'Zone axis geometry and symmetry'

type(EMsoft_T)          :: EMsoft 
type(IO_T)              :: Message 
type(Cell_T)            :: cell 
type(SpaceGroup_T)      :: SG

integer(kind=irg)       :: kk(3), ga(3), gb(3), io_int(6), j, i, dgn
character(fnlen)        :: xtalname

! display the standard program info
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 905 /) ) 

! ask for the crystal structure file
call Message%ReadValue(' Enter xtal file name : ', xtalname,"(A)")
call cell%getCrystalData(xtalname, SG, EMsoft)

! get the zone axis
call Message%ReadValue('Enter the zone axis indices : ',io_int,3)
kk(1:3) = io_int(1:3)  

! determine the point group number and get the ZAP 2-D symmetry
 j = SG%getPGnumber()

! use the new routine to get the whole pattern 2D symmetry group, since that
! is the one that determines the independent beam directions.
 dgn = SG%GetPatternSymmetry(kk, j, .TRUE.)
 
! determine and display the shortest reciprocal lattice vectors for this zone
 call cell%ShortestG(SG, kk, ga, gb, j)
 io_int(1:3) = ga(1:3)
 io_int(4:6) = gb(1:3)
 call Message%WriteValue(' Reciprocal lattice vectors : ', io_int, 6, "(' (',3I3,') and (',3I3,')',/)")

 call Message%printMessage( (/ & 
   ' -> Note that the first of these vectors is by default the horizontal direction in     ', &
   '    any diffraction pattern or image simulation. All reference frames are right-handed.' /), frm = "(A)")

end program EMZAgeom

