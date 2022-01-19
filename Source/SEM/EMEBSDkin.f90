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

program EMEBSDkin
  !! author: MDG
  !! version: 1.0 
  !! date: 11/18/21
  !!
  !! computes kinematical EBSD patterns (with or without deformation)

use mod_kinds
use mod_global
use mod_EMsoft
use mod_HDFnames
use stringconstants
use mod_EBSDkin

IMPLICIT NONE

character(fnlen) :: progname = 'EMEBSDkin.f90'
character(fnlen) :: progdesc = 'kinematical EBSD patterns, using precomputed master Lambert projections'

type(EMsoft_T)   :: EMsoft
type(EBSDkin_T)  :: EBSDkin 
type(HDFnames_T) :: HDFnames

! print the EMsoft header and handle any command line arguments  
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 310 /) )

! deal with the namelist stuff
EBSDkin = EBSDkin_T(EMsoft%nmldeffile)

HDFnames = HDFnames_T()

call HDFnames%set_ProgramData(SC_EBSDkin) 
call HDFnames%set_NMLlist(SC_EBSDkinNameList) 
call HDFnames%set_NMLfilename(SC_EBSDkinNML) 

! perform the pattern computations
call EBSDkin%EBSDkin(EMsoft, progname, HDFnames)

end program EMEBSDkin
