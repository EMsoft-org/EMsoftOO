! ###################################################################
! Copyright (c) 2014-2020, Marc De Graef Research Group/Carnegie Mellon University
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

program EMsampleRFZ
  !! author: MDG
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! Basic program to generate a uniform sampling of Rodrigues Fundamental Zone
  !!
  !! This program calls the SampleRFZ routine of the so3 module to generate
  !! an angle file of euler angles for points that uniformly sample an RFZ for a given
  !! crystal symmetry. 

use mod_kinds
use mod_global
use mod_EMsoft
use mod_sampleRFZ

IMPLICIT NONE

character(fnlen)    :: progname = 'EMsampleRFZ.f90'
character(fnlen)    :: progdesc = 'Create a uniform sampling of Rodrigues space and output angle list'
character(fnlen)    :: nmlpath

type(EMsoft_T)      :: EMsoft
type(sampleRFZ_T)   :: RFZ 

! print some information
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 60 /) )

! deal with the namelist stuff
RFZ = sampleRFZ_T(EMsoft%nmldeffile)

! perform the sampling algorithm
call RFZ%CreateSampling(EMsoft)

end program EMsampleRFZ
