! ###################################################################
! Copyright (c) 2016-2023, Marc De Graef Research Group/Carnegie Mellon University
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

program EMDemag
  !! author: MDG
  !! version: 1.0 
  !! date: 08/17/21
  !!
  !! Demagnetization factor computations

use mod_kinds
use mod_global
use mod_EMsoft
use mod_Demag
use mod_HDFnames 
use stringconstants

IMPLICIT NONE

character(fnlen)  :: progname = 'EMDemag.f90'
character(fnlen)  :: progdesc = 'Computation of demagnetization tensor field and factors'

type(EMsoft_T)    :: EMsoft
type(Demag_T)     :: Demag 
type(HDFnames_T)  :: HDFnames 

! print the EMsoft header and handle any command line arguments  
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 301 /) )

! deal with the namelist stuff
Demag = Demag_T(EMsoft%nmldeffile)

HDFnames = HDFnames_T()

call HDFnames%set_ProgramData(SC_Demag)
call HDFnames%set_NMLlist(SC_DemagNML) 
call HDFnames%set_NMLfilename(SC_DemagNML) 
call HDFnames%set_Variable(SC_Demag)

! perform the computations
call Demag%Demag(EMsoft, progname, HDFnames)

end program EMDemag
