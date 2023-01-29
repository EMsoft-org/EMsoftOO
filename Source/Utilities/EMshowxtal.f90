! ###################################################################
! Copyright (c) 2013-2023, Marc De Graef Research Group/Carnegie Mellon University
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
  !! Display (quasi-)crystal structure information

use mod_kinds
use mod_global
use mod_EMsoft
use mod_symmetry
use mod_io
use mod_HDFsupport
use HDF5
use mod_crystallography
use mod_rotations
use mod_quaternions
use mod_so3
use mod_misc 
use mod_QCsymmetry
use mod_QCcrystallography
use stringconstants

IMPLICIT NONE

character(fnlen)                :: progname = 'EMshowxtal.f90'
character(fnlen)                :: progdesc = 'Display crystal structure information'

type(EMsoft_T)                  :: EMsoft 
type(IO_T)                      :: Message 
type(Cell_T)                    :: cell 
type(SpaceGroup_T)              :: SG
type(QCcell_axial_T)            :: QCcell_axial 
type(QCcell_icosahedral_T)      :: QCcell_icosahedral
type(QCSpaceGroup_T)            :: QCSG
type(HDF_T)                     :: HDF

character(fnlen)                :: xtalname, fname, dataset, groupname
logical                         :: verbose=.TRUE., g_exists
integer(kind=irg)               :: i, j, hdferr, N_Axial
character(1)                    :: yesno
real(kind=dbl),allocatable      :: data(:,:,:), direc(:,:,:), recip(:,:,:)

! header and command line arguments, if any
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 921 /) ) 
 
! ask for the crystal structure file
call Message%ReadValue(' Enter xtal file name [*.xtal, *.qxtal]: ', xtalname,"(A)")

! is this a regular .xtal file or a quasi-crystal file .qxtal? 
i = INDEX(trim(xtalname), 'qxtal')

call openFortranHDFInterface()

if (i.eq.0) then ! regular crystal structure file
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
    call Message%printMessage(' Direct space           Reciprocal space')
    direc = SG%getSpaceGroupPGdirecMatrices()
    recip = SG%getSpaceGroupPGrecipMatrices()
    do i=1,SG%getSpaceGroupNUMpt() 
       write (*,*) i,':'
       write (*,*) (direc(i,1,j),j=1,3),'       ',(recip(i,1,j),j=1,3)
       write (*,*) (direc(i,2,j),j=1,3),'       ',(recip(i,2,j),j=1,3)
       write (*,*) (direc(i,3,j),j=1,3),'       ',(recip(i,3,j),j=1,3)
       write (*,*) ' '
    end do
  endif    
else  ! quasi-crystal structure file
  fname = trim(EMsoft%generateFilePath('EMXtalFolderpathname',trim(xtalname)))

! first we need to check if the N_Axial field exists in this file 
  HDF = HDF_T()

  hdferr = HDF%openFile(fname, readonly = .TRUE.)
  if (hdferr.ne.0) call HDF%error_check('ReadQCDataHDF:HDF%openFile:'//trim(fname), hdferr)

  groupname = SC_CrystalData
  hdferr = HDF%openGroup(groupname)
  if (hdferr.ne.0) call HDF%error_check('ReadQCDataHDF:HDF%openGroup:'//trim(groupname), hdferr)

  ! is this an axial or icosahedral structure?
  dataset = SC_AxialSymmetry
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)

  if (g_exists) then
    QCcell_axial = QCcell_axial_T()
    call HDF%readDatasetInteger(dataset, hdferr, N_Axial) 
    if (hdferr.ne.0) call HDF%error_check('ReadQCDataHDF:HDF%readDatasetInteger:'//trim(dataset), hdferr)
    call HDF%pop(.TRUE.)
    call closeFortranHDFInterface()
   
    if (N_Axial.eq.8)  QCSG = QCspacegroup_T( nD = 2, QCtype = 'Oct')
    if (N_Axial.eq.10) QCSG = QCspacegroup_T( nD = 2, QCtype = 'Dec')
    if (N_Axial.eq.12) QCSG = QCspacegroup_T( nD = 2, QCtype = 'DoD')

    call QCcell_axial%setfname(xtalname)
    call QCcell_axial%ReadQCDataHDF(QCSG, EMsoft)
    ! call QCSG%printSGtable()
    call QCcell_axial%setMetricParametersQC()
    call QCSG%GenerateQCSymmetry(dopg=.FALSE.)
    call QCcell_axial%DumpQXtalInfo(QCSG)
  else
    call HDF%pop(.TRUE.)
    call closeFortranHDFInterface()
    QCcell_icosahedral = QCcell_icosahedral_T()
    call QCcell_icosahedral%setfname(xtalname)
    QCSG = QCspacegroup_T( nD = 3, QCtype = 'Ico')
    call QCcell_icosahedral%ReadQCDataHDF(QCSG, EMsoft)
    ! call QCSG%printSGtable()
    call QCcell_icosahedral%setMetricParametersQC()
    call QCSG%GenerateQCSymmetry(dopg=.FALSE.)
    call QCcell_icosahedral%DumpQXtalInfo(QCSG)
  end if
end if

end program EMshowxtal
