! ###################################################################
! Copyright (c) 2013-2024, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_extractnml
  !! author: MDG 
  !! version: 1.0 
  !! date: 08/18/22
  !!
  !! class definition for the EMextractnml program

use mod_kinds
use mod_global

IMPLICIT NONE 

! class definition
type, public :: extractnml_T
private 
  character(fnlen)  :: hdfname 

contains
private 
  procedure, pass(self) :: extractnml_

  generic, public :: extractnml => extractnml_

end type extractnml_T

! the constructor routine for this class 
interface extractnml_T
  module procedure extractnml_constructor
end interface extractnml_T

contains

!--------------------------------------------------------------------------
type(extractnml_T) function extractnml_constructor( hdffile ) result(extractnml)
!DEC$ ATTRIBUTES DLLEXPORT :: extractnml_constructor
!! author: MDG 
!! version: 1.0 
!! date: 08/18/22
!!
!! constructor for the extractnml_T Class; checks for the existence of the hdffile
 
 use mod_io 

IMPLICIT NONE

character(fnlen)  :: hdffile 

type(io_T)        :: Message  
logical           :: f_exists

extractnml%hdfname = trim(hdffile)

! make sure the file exists, if not we exit
inquire(file=trim(hdffile), exist=f_exists)

if (f_exists.eqv..FALSE.) call Message%printError('extractnml_constructor','HDF5 file not found')

end function extractnml_constructor

!--------------------------------------------------------------------------
subroutine extractnml_destructor(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: extractnml_destructor
!! author: MDG 
!! version: 1.0 
!! date: 08/18/22
!!
!! destructor for the extractnml_T Class
 
IMPLICIT NONE

type(extractnml_T), INTENT(INOUT)  :: self 

call reportDestructor('extractnml_T')

end subroutine extractnml_destructor

!--------------------------------------------------------------------------
subroutine extractnml_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: extractnml_
!! author: MDG 
!! version: 1.0 
!! date: 08/18/22
!!
!! perform the computations

use mod_EMsoft
use HDF5
use mod_HDFsupport
use stringconstants
use mod_io

IMPLICIT NONE 

class(extractnml_T), INTENT(INOUT)      :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 

type(HDF_T)                             :: HDF 
type(IO_T)                              :: Message 

character(fnlen)                        :: groupname, name, dataset, nmlname
logical                                 :: g_exists
integer(kind=irg)                       :: hdferr, storage_type, nlinks, max_corder, io_int(1)
integer(HSIZE_T)                        :: n

! open the HDF5 interface
call openFortranHDFInterface()

! instantiate the class
HDF = HDF_T()

! open the file in read_only mode
hdferr = HDF%openFile( self%hdfname, readonly=.TRUE. )

! check for the existence of the NMLfiles group 
groupname = SC_NMLfiles
  call H5Lexists_f(HDF%getobjectID(), trim(groupname), g_exists, hdferr)

! if the NMLfiles group exists, then check how many datasets it has
  if (g_exists) then
    hdferr = HDF%openGroup(groupname)
    call h5gget_info_f(HDF%getobjectID(), storage_type, nlinks, max_corder, hdferr)
    io_int(1) = nlinks
    call Message%WriteValue(' Number of namelist files found : ', io_int, 1)
    do n=1, nlinks 
! get the data size (string length)
      call h5lget_name_by_idx_f(HDF%getobjectID(), '.', H5_INDEX_NAME_F, H5_ITER_NATIVE_F, &
                              n-1, name, hdferr )
      dataset = trim(name)
      nmlname = trim(name)//'_extracted.nml'
      call Message%printMessage( ' extracting '//trim(name)//' into file '//trim(nmlname ) )
      hdferr = HDF%extractDatasetTextfile( dataset, nmlname )
    end do 
  else 
    call Message%printError('extractnml_','This file does not have an NMLfiles group')
  end if 

  call HDF%popall()

call closeFortranHDFInterface()

end subroutine extractnml_



end module mod_extractnml