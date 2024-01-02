! ###################################################################
! Copyright (c) 2016-2024, Marc De Graef Research Group/Carnegie Mellon University
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

program EMextractnml
  !! author: MDG
  !! version: 1.0 
  !! date: 08/18/22
  !!
  !! extract nml files from an EMsoft HDF5 file

use mod_kinds
use mod_global
use mod_EMsoft
use mod_extractnml
use mod_io

IMPLICIT NONE

character(fnlen)                :: progname = 'EMextractnml.f90'
character(fnlen)                :: progdesc = 'Extract nml files from an EMsoft HDF5 file'

type(EMsoft_T)                  :: EMsoft
type(extractnml_T)              :: exnml
type(IO_T)                      :: Message 

character(fnlen)                :: hdffile

! print the EMsoft header 
EMsoft = EMsoft_T( progname, progdesc )

! ask the user for the hdf filename 
call Message%ReadValue('Enter the full path to the HDF5 file from which to extract namelist files:', hdffile)

! instantiate the class
exnml = extractnml_T( hdffile )

! if we get here then we can do the extraction
call exnml%extractnml( EMsoft, progname)

end program EMextractnml
