! ###################################################################
! Copyright (c) 2013-2024, Marc De Graef/Carnegie Mellon University
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

!--------------------------------------------------------------------------
! EMsoft: EMSRCBED.f90
!--------------------------------------------------------------------------
!
! PROGRAM: EMSRCBED 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Systematic Row Convergent Beam Electron Diffraction 
!>        using scattering matrix approach
! 
!> @date 06/07/01 MDG 1.0 original (for the CTEMbook)
!> @date 04/18/13 MDG 2.0 rewrite
!> @date 02/16/13 MDG 3.0 complete rewrite for EMsoftOO
!--------------------------------------------------------------------------
program EMSRCBED

use mod_kinds
use mod_global
use mod_EMsoft
use mod_SRCBED

IMPLICIT NONE


character(fnlen)    :: progname = 'EMSRCBED.f90'
character(fnlen)    :: progdesc = 'CTEMbook: Systematic row convergent beam pattern using scattering matrix'

type(EMsoft_T)      :: EMsoft
type(SRCBED_T)      :: SRCBED


EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 103 /) )

! generate a set of systematic row CBED patterns
SRCBED = SRCBED_T( EMsoft%nmldeffile )

! perform the computations
call SRCBED%SRCBED( EMsoft, progname )
 
end program EMSRCBED
