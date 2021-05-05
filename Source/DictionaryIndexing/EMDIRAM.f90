! ###################################################################
! Copyright (c) 2016-2020, Marc De Graef Research Group/Carnegie Mellon University
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

program EMDIRAM
  !! author: MDG
  !! version: 1.0 
  !! date: 04/30/21
  !!
  !! In-RAM Indexing of EBSD/ECP/TKD patterns using a dynamically generated pattern dictionary

use mod_kinds
use mod_global
use mod_EMsoft
use mod_DI
use mod_HDFsupport
use ISO_C_BINDING

IMPLICIT NONE

character(fnlen)        :: progname
character(fnlen)        :: progdesc = 'In-RAM Indexing of EBSD/ECP/TKD patterns using a dynamically calculated pattern dictionary'

type(EMsoft_T)          :: EMsoft

progname = 'EMDIRAM.f90'

! print the EMsoft header and handle any command line arguments  
! this program uses the same name list file as the EMDI program
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 105 /) )

! call the DIRAMdriver routine to take care of the entire indexing process 
call DIRAMdriver(EMsoft%nmldeffile, progname, progdesc)

end program EMDIRAM
