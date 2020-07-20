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

module mod_HDFFileInfo
  !! author: MDG
  !! version: 1.0
  !! date: 04/17/20
  !!
  !! helper routines for the EMHDFFileInfo program

use mod_kinds
use mod_global
use HDF5
use mod_HDFsupport

IMPLICIT NONE


contains


!--------------------------------------------------------------------------
subroutine HDFanalyze(HDFname, ProgramName)
!DEC$ ATTRIBUTES DLLEXPORT :: HDFanalyze
!! author: MDG
!! version: 1.0
!! date: 04/17/20
!!
!! perform the data set analysis; every EMsoft file class must have a getFileInfo method
!! that can be called by the code below to display information regarding the file content.

use mod_io
use mod_MPfiles
use mod_MCfiles
use mod_DIfiles

IMPLICIT NONE

character(fnlen), INTENT(INOUT)       :: HDFName
character(fnlen), INTENT(INOUT)       :: ProgramName

type(IO_T)                            :: Message
type(MPfile_T)                        :: MPFT
type(MCfile_T)                        :: MCFT
type(DIfile_T)                        :: DIFT


call Message%printMessage(' Starting analysis of '//trim(ProgramName)//' data','(/A)')
call Message%printMessage(' ============================================= ')

! we'll select the relevant routine through a simple select statement
select case(ProgramName)
  case('EMEBSDmaster.f90')
    MPFT = MPfile_T()
    call MPFT%setModality('EBSD')
    call MPFT%setFileName(HDFname)
    call MPFT%getFileInfo('EBSD')

  case('EMTKDmaster.f90')
    MPFT = MPfile_T()
    call MPFT%setModality('TKD')
    call MPFT%setFileName(HDFname)
    call MPFT%getFileInfo('TKD')

  case('EMECPmaster.f90')
    MPFT = MPfile_T()
    call MPFT%setModality('ECP')
    call MPFT%setFileName(HDFname)
    call MPFT%getFileInfo('ECP')

  case('EMMCOpenCL.f90')
    MCFT = MCfile_T()
    call MCFT%setFileName(HDFname)
    call MCFT%getFileInfo()

  case('EMDI.f90')
    ! DIFT = DIfile_T()
    ! call MCFT%setFileName(HDFname)
    ! call MCFT%getFileInfo()



  case default
    call Message%printMessage(' --> not yet implemented ')
end select

end subroutine HDFanalyze

end module mod_HDFFileInfo
