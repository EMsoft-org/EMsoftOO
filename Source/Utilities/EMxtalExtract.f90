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

!--------------------------------------------------------------------------
! EMsoft:EMxtalExtract.f90
!--------------------------------------------------------------------------
!
! PROGRAM: EMmkxtalExtract 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief extract information from a crystal structure file and put in txt file for input redirect
!
!> @date   07/19/18 MDG 1.0 original
!--------------------------------------------------------------------------
program EMxtalExtract
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! extract information from a crystal structure file and put in txt file for input redirect

use mod_EMsoft
use mod_kinds
use mod_global 
use mod_symmetry 
use mod_io 
use mod_crystallography

IMPLICIT NONE

type(EMsoft_T)                  :: EMsoft 
type(Cell_T)                    :: cell 
type(SpaceGroup_T)              :: SG
type(IO_T)                      :: Message 

character(fnlen)                :: progname = 'EMxtalExtract.f90'
character(fnlen)                :: progdesc = 'Extract all information from an HDF crystal structure file &
                                               and dump it in regular or Wyckoff form to a text file'
 

character(fnlen)                :: xtalname, txtname
integer(kind=irg)               :: i, io_int(3)
real(kind=dbl)                  :: io_real(5)
real(kind=dbl), allocatable     :: asp(:,:)
integer(kind=irg), allocatable  :: atp(:)

 EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 926 /) )

! read crystal information
 call Message%ReadValue(' Enter xtal file name : ', xtalname,"(A)")
 call cell%getCrystalData(xtalname, SG, EMsoft)

! get a filename for the output 
 call Message%ReadValue('Enter text file name for output (*.txt) ', txtname)
 open (unit=dataunit, file=trim(txtname), status='unknown',form='formatted')

! 1. write the crystal system
 if (SG%getSpaceGrouptrigonal().eqv..TRUE.) then
    write(dataunit,"(I1)") 5
    if (SG%getSpaceGroupXtalSystem().eq.4) then
      write(dataunit,"(I1)") 1
    else
      write(dataunit,"(I1)") 0
    end if
 else
    write(dataunit,"(I1)") SG%getSpaceGroupXtalSystem()
 end if

! 2. lattice parameters (only the ones that are needed)
! a is always needed
 write(dataunit,"(F10.5)") cell%getLatParm('a')
 select case (SG%getSpaceGroupXtalSystem())
    case (1)
  ! tetragonal
    case (2)
     write(dataunit,"(F10.5)") cell%getLatParm('c')
  ! orthorhombic
    case (3)
     write(dataunit,"(F10.5)") cell%getLatParm('b')
     write(dataunit,"(F10.5)") cell%getLatParm('c')
  ! hexagonal
    case (4)
     write(dataunit,"(F10.5)") cell%getLatParm('c')
  ! rhombohedral 
    case (5)
     write(dataunit,"(F10.5)") cell%getLatParm('alpha')
  ! monoclinic   
    case (6)
     write(dataunit,"(F10.5)") cell%getLatParm('b')
     write(dataunit,"(F10.5)") cell%getLatParm('c')
     write(dataunit,"(F10.5)") cell%getLatParm('beta')
  ! triclinic    
    case (7) 
     write(dataunit,"(F10.5)") cell%getLatParm('b')
     write(dataunit,"(F10.5)") cell%getLatParm('c')
     write(dataunit,"(F10.5)") cell%getLatParm('alpha')
     write(dataunit,"(F10.5)") cell%getLatParm('beta')
     write(dataunit,"(F10.5)") cell%getLatParm('gamma')
 end select

! 3. space group number
 write(dataunit,"(I3)") SG%getSpaceGroupNumber()

! 4. some space groups have a second setting
 if (minval(tworig-SG%getSpaceGroupNumber()).eq.0) then 
   write(dataunit,"(I3)") SG%getSpaceGroupSetting()
 end if

! 5. the atom coordinates, site occupations, and Debye-Waller factors are output here...
 asp = cell%getAsymPosData()
 atp = cell%getAtomtype()
 do i=1,cell%getNatomtype()
   if (i.gt.1) write(dataunit,"('y')")
   write(dataunit,"(I2)") atp(i)
   write(dataunit,"(4(F10.5,','),F10.5)") asp(i,1:5)
 end do
 write(dataunit,"('n')")

! 6. xtal file name
 write(dataunit,"(A)") trim(cell%getFileName())

! 7. source, if defined
 if (trim(cell%getSource()).eq.'undefined') then
   write(dataunit,"(A)") ''''''
   call Message%printMessage('')
   call Message%printMessage('=========================')
   call Message%printMessage('The current version of the structure file has an empty Source field.')
   call Message%printMessage('Please edit the text file and add a citation describing where the structure data &
                              comes from (max 512 characters).')
   call Message%printMessage('=========================')
   call Message%printMessage('')
 else
   write(dataunit,"('''',A,'''')") trim(cell%getSource())
 end if 

 close(unit=dataunit,status='keep')

 call Message%printMessage('')
 call Message%printMessage('The text output file '//trim(txtname)//' has been created.')
 call Message%printMessage('You can edit this file with a text editor to make corrections,')
 call Message%printMessage('or to add a Source string on the last line (if not already present) between single quotes') 
 call Message%printMessage('Then, use the following command to update the actual structure file:')
 call Message%printMessage('')
 call Message%printMessage('EMmkxtal < '//trim(txtname))
 call Message%printMessage('')

end program EMxtalExtract
