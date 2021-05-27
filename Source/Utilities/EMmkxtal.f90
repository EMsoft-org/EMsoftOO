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
program EMmkxtal
  !! author: MDG
  !! version: 1.0
  !! date: 01/14/20
  !!
  !! generate a .xtal file; tested on 01/14/20.

use mod_kinds
use mod_global
use mod_EMsoft
use mod_io
use mod_symmetry
use mod_crystallography 
use mod_HDFsupport

IMPLICIT NONE

character(fnlen)        :: progname = 'EMmkxtal.f90'
character(fnlen)        :: progdesc = 'Create an HDF crystal structure file and place it in the XtalFolder'

type(EMsoft_T)          :: EMsoft
type(IO_T)              :: Message
type(SpaceGroup_T)      :: SG
type(Cell_T)            :: cell

character(fnlen)        :: flag   ! we need to test for the -w Wyckoff positions command line argument
character(fnlen)        :: fname, source
logical                 :: useWyckoff  = .FALSE.

! initialize the cell and IO classes
cell = Cell_T()
Message = IO_T()

! display the splash screen and analyze any non-standard command line arguments
EMsoft%flagset = '-w'
EMsoft = EMsoft_T(progname, progdesc, tpl = (/ 917 /))
if (trim(EMsoft%flagset).eq.'yes') then
    useWyckoff = .TRUE.
    call cell%setWyckoff(useWyckoff)
end if

! generate the space group class; this generates all symmetry arrays
SG = SpaceGroup_T()

! to avoid circular references, the cell class also needs to know the crystal system
call cell%setXtalSystem(SG%getSpaceGroupXtalSystem())

! get the lattice parameters
call cell%setLatParm(SG)

! get the atom positions, either using Wyckoff positions or the regular way
if (useWyckoff) then
  call cell%GetAsymPosWyckoff(SG, EMsoft)
else
  call cell%GetAsymPos()
end if

! ask for the .xtal file name
call Message%ReadValue('Enter output file name (*.xtal) ', fname)
call cell%setFileName(fname)

! and for the source of the crystal data
call Message%ReadValue('Enter the source for this data [max. 512 characters, quoted] ', source)
call cell%setSource(source)

! finally, save the file
call openFortranHDFInterface()
call cell%SaveDataHDF(SG, EMsoft)
call closeFortranHDFInterface()

end program EMmkxtal
