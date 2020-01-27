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

module mod_TB
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/27/20
  !!
  !! Everything that has to do with analytical two-beam solutions 

use mod_kinds 
use mod_global
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit

IMPLICIT NONE

contains 

!--------------------------------------------------------------------------
recursive subroutine TBCalcSM(Ar,Ai,sg,z,xig,xigp,xizero,betag)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/27/20
  !!
  !! 2-beam scattering matrix implementation
  !!
  !! compute crystal scattering matrix (real and imaginary
  !! part) for a given set of parameters. Optimized to minimize the 
  !! total number of function evaluations and multiplications/divisions
  !!
!DEC$ ATTRIBUTES DLLEXPORT :: TBCalcSM


IMPLICIT NONE

real(kind=sgl),INTENT(IN)       :: sg                   
 !! excitation error
real(kind=sgl),INTENT(IN)       :: z                    
 !! thickness
real(kind=sgl),INTENT(IN)       :: xig                  
 !! extinction distance
real(kind=sgl),INTENT(IN)       :: xigp                 
 !! anomalous absorption length
real(kind=sgl),INTENT(IN)       :: xizero               
 !! normal absorption length
real(kind=sgl),INTENT(IN)       :: betag                
 !! phase parameter
real(kind=sgl),INTENT(OUT)      :: Ar(2,2)              
 !! real part of result 
real(kind=sgl),INTENT(OUT)      :: Ai(2,2)              
 !! imaginary part of result 

real(kind=sgl)  :: pr, pi, cs, ss, ch, sh, q, q1, q2, sgs, sr, si, o , p, sb, cb, e, r, sq, xigi, xigpi

! setup auxiliary variables 
 xigi = 1.00/xig
 xigpi = 1.00/xigp

! sigma squared
 cb = cos(betag)
 sb = sin(betag)
 q = sg**2+xigi**2-xigpi**2
 r = cb*xigi*xigpi
 sq = sqrt(q**2+4.0*r**2)

! real part of sigma
 sr = sqrt(0.5*(q+sq))

! imaginary part of sigma
 si = r/sr
 sq = 1.0/sq

! s_g divided by sigma squared
 sgs = sg*sq

! arguments of trigonometric and hyperbolic functions
 e = cPi*z
 pr = e*sr
 pi = e*si
 e = exp(-e/xizero)

! trigonometric and hyperbolic functions
 cs = cos(pr)
 ss = sin(pr)
 ch = cosh(pi)
 sh = sinh(pi)
 o = ss*ch
 p = cs*sh

! transmitted amplitude T including normal absorption
 q = e*sgs
 q1 = q*(si*o-sr*p)
 q2 = q*(sr*o+si*p)
 Ar(1,1) = e*cs*ch
 Ai(1,1) = -e*ss*sh
 Ar(2,2) = Ar(1,1)+q1
 Ai(2,2) = Ai(1,1)+q2
 Ar(1,1) = Ar(1,1)-q1
 Ai(1,1) = Ai(1,1)-q2

! scattered amplitude S including normal absorption
 q1 = e*sq*(si*xigi-(sr*cb+si*sb)*xigpi)
 q2 = e*sq*(sr*xigi-(sr*sb-si*cb)*xigpi)
 Ar(1,2) = q1*o-q2*p
 Ai(1,2) = q2*o+q1*p
 Ar(2,1) = Ar(1,2)
 Ai(2,1) = Ai(1,2)

end subroutine TBCalcSM

!--------------------------------------------------------------------------
recursive subroutine TBCalcInten(It,Is,sg,z,xig,xigp,xizero,betag)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/27/20
  !!
  !! 2-beam transmitted and scattered intensities
  !!
  !! compute transmitted and scattered intensities for
  !! the perfect crystal case.   This routine does not make use of 
  !! complex number arithmetic, but instead uses the analytical 
  !! expressions for the two-beam intensities derived in section 
  !! 6.3.3.4 on page 356-357.
  !!
!DEC$ ATTRIBUTES DLLEXPORT :: TBCalcInten

IMPLICIT NONE

real(kind=sgl),INTENT(IN)       :: sg                   
 !! excitation error
real(kind=sgl),INTENT(IN)       :: z                    
 !! thickness
real(kind=sgl),INTENT(IN)       :: xig                  
 !! extinction distance
real(kind=sgl),INTENT(IN)       :: xigp                 
 !! anomalous absorption length
real(kind=sgl),INTENT(IN)       :: xizero               
 !! normal absorption length
real(kind=sgl),INTENT(IN)       :: betag                
 !! phase parameter
real(kind=sgl),INTENT(OUT)      :: It                   
 !! real part of result 
real(kind=sgl),INTENT(OUT)      :: Is                   
 !! imaginary part of result 

real(kind=sgl) :: q, r, sq, qgsi, e, sr, si, cp, ch, pr, pi, xigi, xigpi, sgs
     
! setup auxiliary variables 
 xigi = 1.0/xig
 xigpi = 1.0/xigp

! sigma squared
 q = sg**2+xigi**2-xigpi**2
 r = cos(betag)*xigi*xigpi
 sq = sqrt(q**2+4.0*r**2)

! real part of sigma
 sr = sqrt(0.5*(q+sq))

! imaginary part of sigma
 si = r/sr
 sq = 1.0/sq

! reciprocal of q_g squared
 qgsi = xigi**2+xigpi**2-2.0*sin(betag)*xigi*xigpi

! s_g squared divided by sigma squared
 sgs = sg**2*sq

! arguments of trigonometric and hyperbolic functions
 e = 2.0*cPi*z
 pr = e*sr
 pi = e*si
 e = exp(-e/xizero)

! trigonometric functions 
 cp = cos(pr)
 ch = cosh(pi)

! transmitted intensity It
 It = 0.5*ch*(1.0+sgs)+sg*sq*(sr*sinh(pi)-si*sin(pr))+0.5*cp*(1.0-sgs)
 It = It*e

! scattered intensity Is
 Is = 0.5*qgsi*sq*e*(ch-cp)

end subroutine TBCalcInten

end module mod_TB