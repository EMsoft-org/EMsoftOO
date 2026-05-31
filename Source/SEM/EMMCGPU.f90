! ###################################################################
! Copyright (c) 2016-2026, Marc De Graef Research Group/Carnegie Mellon University
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

program EMMCGPU
  !! author: MDG
  !! version: 1.0
  !! date: 02/04/20
  !!
  !! Monte Carlo backscattered electron simulation (GPU-accelerated).
  !! Formerly EMMCOpenCL; renamed because the GPU compute now runs on either
  !! OpenCL or Apple Metal.  The old EMMCOpenCL command is still built as an
  !! alias, and the program identity embedded in the output HDF5 files is kept
  !! as 'EMMCOpenCL.f90' (see MCfileProgName below) so that existing Monte Carlo
  !! files remain compatible with the master-pattern programs that read them.

use mod_kinds
use mod_global
use mod_EMsoft
use mod_MCOpenCL

IMPLICIT NONE

character(fnlen)     :: progname = 'EMMCGPU.f90'
character(fnlen)     :: progdesc = 'Monte Carlo backscattered electron simulation'

! Program name written into the output HDF5 file header.  This MUST stay
! 'EMMCOpenCL.f90' for backward/forward compatibility: downstream readers
! (mod_HDFFileInfo, the EBSD/ECP/TKD master programs) match Monte Carlo files on
! this exact string, and existing files in the wild carry it.
character(fnlen)     :: MCfileProgName = 'EMMCOpenCL.f90'

type(EMsoft_T)       :: EMsoft
type(MCOpenCL_T)     :: MCCL

! print the EMsoft header and handle any command line arguments
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 42 /) )

! deal with the namelist stuff
MCCL = MCOpenCL_T(EMsoft%nmldeffile)

! perform the sampling algorithm
call MCCL%MCOpenCL(EMsoft, MCfileProgName)

end program EMMCGPU
