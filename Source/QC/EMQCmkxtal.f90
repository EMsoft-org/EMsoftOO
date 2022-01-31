! ###################################################################
! Copyright (c) 2016-2022, Marc De Graef Research Group/Carnegie Mellon University
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

program EMQCmkxtal
  !! author: MDG/SS
  !! version: 1.0 
  !! date: 01/30/22
  !!
  !! create a quasicrystal structure file

use mod_kinds
use mod_global
use mod_EMsoft
use mod_QCcrystallography
use mod_QCsymmetry
use mod_io

IMPLICIT NONE

character(fnlen)     :: progname = 'EMQCmkxtal.f90'
character(fnlen)     :: progdesc = 'Create an axial or icosahedral quasi-crystal structure file'

type(EMsoft_T)       :: EMsoft
type(QCcell_T)       :: QCcell
! type(TDQCcell_T)     :: TDQCcell
type(QCspacegroup_T) :: QCSG
type(IO_T)           :: Message

integer(kind=irg)    :: qcdim, ii, jj, io_int(1)
character(fnlen)     :: fname

! print the EMsoft header and handle any command line arguments  
EMsoft = EMsoft_T( progname, progdesc )

! Is this a 2D or 3D quasicrystal ?  [was GetQCType in EMsoft 5.0]
call Message%printMessage(' Select the quasicrystal dimensionality : ', frm = "(A)")
call Message%printMessage('  2-dimensional (axial) quasicrystal', frm = "(A/)")
call Message%printMessage('  3-dimensional (icosahedral) quasicrystal', frm = "(A//)")
call Message%ReadValue(' quasi-crystal dimensionality ---> ', io_int, 1)
qcdim = io_int(1) 

! initialize the appropriate QC cell class and construct the HDF5 .qxtal file
 select case (qcdim)
  case(2)
    ! 2D QC
    ! TDQCcell = TDQCcell_T()
    ! call GetQCLatParm(TDQCcell)
    ! call PrintSGTable(TDQCcell)
    ! call GetQCSpaceGroup(TDQCcell)
    ! call GetQCAsymPos(TDQCcell)
    ! call ReadValue('Enter output file name (*.qxtal) ', fname)
    ! TDQCcell%fname = fname
    ! call SaveQCDataHDF(TDQCcell)

  case(3)
    ! 3D QC
    QCcell = QCcell_T()
    QCSG = QCspacegroup_T()
    call QCcell%GetQCLatParm()
    call QCSG%PrintSGTable()
    call QCSG%GetQCSpaceGroup()
    call QCcell%GetQCAsymPos()
    call Message%ReadValue('Enter output file name (*.qxtal) ', fname)
    call QCcell%setfname(fname)
    call QCcell%SaveQCDataHDF(QCSG, EMsoft)

  case DEFAULT
    ! anything else
    call Message%printError('EMQCmkxtal:','unknown quasicrystal dimensionality (only 2 and 3 are implemented)')

 end select

! perform the computations
! call YYY%QCmkxtal(EMsoft, progname)

end program EMQCmkxtal
