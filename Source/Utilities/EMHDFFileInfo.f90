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

program EMHDFFileInfo
  !! author: MDG
  !! version: 1.0 
  !! date: 04/17/20
  !!
  !! program to print out information about an arbitrary EMsoft HDF file, including 
  !! .xtal files

use mod_kinds
use mod_global
use mod_EMsoft
use mod_crystallography
use mod_symmetry
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_io
use mod_HDFFileInfo
use ISO_C_BINDING 

IMPLICIT NONE

character(fnlen)                :: progname = 'EMHDFFileInfo.f90'
character(fnlen)                :: progdesc = 'Program to print out information about an arbitrary EMsoft HDF file'

type(EMsoft_T)                  :: EMsoft
type(IO_T)                      :: Message
type(HDF_T)                     :: HDF 
type(HDFnames_T)                :: HDFnames 
type(cell_T)                    :: cell 
type(SpaceGroup_T)              :: SG

integer(kind=irg)               :: numarg, hdferr, nmembers, i, j, otype, io_int(1), PNindex, nlines
character(fnlen)                :: name_buffer 
character(fnlen)                :: HDFname, groupname, dataset 
character(fnlen),allocatable    :: membernames(:), prognames(:)    
integer(kind=irg),allocatable   :: membertypes(:)
logical                         :: f_exists, g_exists, stat, singleprogram, crystalfile
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)

! print the EMsoft header and handle any command line arguments  
EMsoft = EMsoft_T( progname, progdesc )

! there should be a command line argument (an HDF5 filename)
numarg = command_argument_count()
if (numarg.eq.0) then 
    call Message%ReadValue(' Please enter an HDF5 filename : ', HDFname)
else 
    call get_command_argument(1,HDFname)
end if 

! make sure the file exists
call Message%printMessage(' Checking to see if this file exists ... ',"(A$)")
inquire(file=trim(HDFname), exist=f_exists)
if (f_exists.eqv..TRUE.) then 
    call Message%printMessage(' it does ! ')
else 
    call Message%printMessage(' hmmmm... ')
    call Message%printError('EMHDFFileInfo','file not found ...')
end if 

! is this a proper HDF5 file ?
call h5fis_hdf5_f(trim(HDFname), stat, hdferr)
if (stat.eqv..FALSE.) then 
    call Message%printError('EMHDFFileInfo','this is not an HDF5 file')
end if 

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T() 
HDFnames = HDFnames_T()

! open the file 
hdferr =  HDF%openFile(HDFname, readonly=.TRUE.)
if (hdferr.ne.0) then 
    call Message%printError('EMHDFFileInfo','error opening file in read-only mode')
end if 

! Which EMsoft program created this file ? 
! First look for the EMheader group
groupname = trim(HDFnames%get_EMheader())
crystalfile = .FALSE.
call Message%printMessage(' Looking for group '//trim(groupname)//' ... ',"(/A$)")
    call H5Lexists_f(HDF%getobjectID(), trim(groupname), g_exists, hdferr)
if (g_exists.eqv..FALSE.) then 
! this could be an .xtal file which does not have an EMheader group 
    call Message%printMessage(' not found ')
    groupname = 'CrystalData'
    call Message%printMessage(' Looking for group '//trim(groupname)//' ... ',"(/A$)")
    call H5Lexists_f(HDF%getobjectID(), trim(groupname), g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then 
        call Message%printMessage(' found it!',"(A/)")
        crystalfile = .TRUE.
    else 
        call Message%printError('EMHDFFileInfo','EMheader and CrystalData groups not found; is this really an EMsoft file?')
    end if 
else
    call Message%printMessage(' found it!',"(A/)")
end if 

! the HDFanalyze routine does all the work of determining which program generated 
! the data in the HDFname file; it also then calls the appropriate routine to 
! extract basic information from the HDF5 file and show it to the user.
! first open the EMheader/CrystalData group
if (crystalfile.eqv..TRUE.) then 
! close the file and call the regular crystal data display routine ...
    call HDF%pop(.TRUE.)
    call cell%readDataHDF(SG, EMsoft, useHDF=HDF, useXtalName=HDFname)
    call cell%dumpXtalInfo(SG)
else 
! ok, if we get here we have a genuine EMsoft file, so let's get a list of any and 
! all the groups in EMheader / CrystalData
    call h5gn_members_f(HDF%getobjectID(), trim(groupname), nmembers, hdferr)
    allocate( membernames(0:nmembers-1), membertypes(0:nmembers-1) )

! Get each group member's name and type
    do i=0,nmembers-1
        call h5gget_obj_info_idx_f(HDF%getobjectID(), trim(groupname), i, name_buffer, otype, &
                                hdferr)
        membernames(i) = trim(name_buffer)
        membertypes(i) = otype
    end do

! there are two possibilities: 
!   - EMheader has a ProgramName parameter in it 
!   - EMheader has one or more subgroups in it 

! look for ProgramName in the membernames array 
    PNindex = -1
    singleprogram = .FALSE.
    do i=0,nmembers-1 
        if (trim(membernames(i)).eq.'ProgramName') then 
            singleprogram = .TRUE.
            PNindex = i 
        end if  
    end do 

    groupname = trim(HDFnames%get_EMheader())
    hdferr = HDF%openGroup(groupname)

    if (singleprogram.eqv..TRUE.) then 
! get the content of the ProgramName data set 
        dataset = trim(membernames(PNindex))
        call HDF%readDatasetStringArray(dataset, nlines, hdferr, stringarray)
        progname = trim(stringarray(1))
        deallocate(stringarray)

        call HDFanalyze(HDFname, progname)
    else ! loop over all the members that are groups, 
         ! extract their ProgramName parameter, then close the file, then get all the information
        allocate(prognames(0:nmembers-1))
        j=0
        do i=0,nmembers-1
            if (membertypes(i).eq.H5G_GROUP_F) then  ! this is a group, so we need to analyze it
                groupname = trim(membernames(i))
                hdferr = HDF%openGroup(groupname)
    ! get the content of the ProgramName data set 
                dataset = 'ProgramName'
                call HDF%readDatasetStringArray(dataset, nlines, hdferr, stringarray)
                prognames(j) = trim(stringarray(1))
                deallocate(stringarray)
                call HDF%pop()
                j = j+1
            end if 
        end do
! close the file 
        call HDF%pop(.TRUE.)
        call Message%printMessage(' The following programs have contributed to this file:',"(/A)")
        do i=0,j-1
            call Message%printMessage(' - '//trim(prognames(i)))
        end do

        do i=0,j-1
            call HDFanalyze(HDFname, prognames(i))
        end do
    end if 
end if 

! close the HDF interface
call closeFortranHDFInterface()

end program EMHDFFileInfo
