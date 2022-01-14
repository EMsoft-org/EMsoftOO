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

!--------------------------------------------------------------------------
! EMsoft:EMshowxtal.f90
!--------------------------------------------------------------------------
!
! PROGRAM: EMshowxtal 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Display crystal structure information
!
!> @date   07/31/18 MDG 1.0 original
!> @date   12/23/18 MDG 1.1 minor reorganization
!--------------------------------------------------------------------------
program EMshowxtal
  !! author: MDG
  !! version: 1.0 
  !! date: 01/23/20
  !!
  !! Display crystal structure information

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

IMPLICIT NONE

character(fnlen)                :: progname = 'EMshowxtal.f90'
character(fnlen)                :: progdesc = 'Display crystal structure information'

type(EMsoft_T)                  :: EMsoft 
type(IO_T)                      :: Message 
type(Cell_T)                    :: cell 
type(SpaceGroup_T)              :: SG

character(fnlen)                :: xtalname
logical                         :: verbose=.TRUE.
integer(kind=irg)                      :: i, j
character(1)                    :: yesno
real(kind=dbl),allocatable      :: data(:,:,:), direc(:,:,:)

! header and command line arguments, if any
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 921 /) ) 
 
! ask for the crystal structure file
call Message%ReadValue(' Enter xtal file name : ', xtalname,"(A)")
call cell%getCrystalData(xtalname, SG, EMsoft, verbose=.TRUE.)
 

call Message%ReadValue('Do you want to print the symmetry matrices as well ? (y/n) ', yesno)

if (yesno.eq.'y') then
  call Message%printMessage('Space group operators (last column = translation)')
  data = SG%getSpaceGroupDataMatrices()
  do i=1,SG%getSpaceGroupMATnum() 
     write (*,*) i,':'
     write (*,*) (data(i,1,j),j=1,4)
     write (*,*) (data(i,2,j),j=1,4)
     write (*,*) (data(i,3,j),j=1,4)
     write (*,*) ' '
  end do

  call Message%printMessage('Point group operators')
  direc = SG%getSpaceGroupPGdirecMatrices()
  do i=1,SG%getSpaceGroupNUMpt() 
     write (*,*) i,':'
     write (*,*) (direc(i,1,j),j=1,3)
     write (*,*) (direc(i,2,j),j=1,3)
     write (*,*) (direc(i,3,j),j=1,3)
     write (*,*) ' '
  end do
endif    

end program EMshowxtal
