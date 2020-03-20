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

program EMTKD
  !! author: MDG
  !! version: 1.0 
  !! date: 03/20/20
  !!
  !! Dynamical TKD patterns, using precomputed MC and master Lambert projections
  !!
  !! we use the regular EBSD code to do this, since the only real difference is the 
  !! detector geometry, and the Monte Carlo input file 

use mod_kinds
use mod_global
use mod_EMsoft
use mod_HDFnames
use stringconstants
use mod_EBSD

IMPLICIT NONE

character(fnlen)    :: progname = 'EMTKD.f90'
character(fnlen)    :: progdesc = 'Dynamical TKD patterns, using precomputed MC and master Lambert projections'

type(HDFnames_T)    :: HDFnames
type(EMsoft_T)      :: EMsoft
type(EBSD_T)        :: TKD 
logical             :: isTKD = .TRUE.

! print the EMsoft header and handle any command line arguments  
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 75 /) )

! deal with the namelist stuff; we'll use the EBSD name list for this case
TKD = EBSD_T(EMsoft%nmldeffile, isTKD)

! set the HDFnames class to TKD mode  
HDFnames = HDFnames_T() 
call HDFnames%set_ProgramData(SC_TKD) 
call HDFnames%set_NMLlist(SC_TKDNameList) 
call HDFnames%set_NMLfilename(SC_TKDNML) 
call HDFnames%set_Variable(SC_MCOpenCL) 

! perform the computations
call TKD%EBSD(EMsoft, progname, HDFnames, TKD=.TRUE.)

end program EMTKD
