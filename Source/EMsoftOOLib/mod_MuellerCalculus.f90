! ###################################################################
! Copyright (c) 2013-2022, Marc De Graef selfearch Group/Carnegie Mellon University
! All rights selferved.
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
! AND ANY EXPselfS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
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
! EMsoft:MuellerCalculus.f90
!--------------------------------------------------------------------------
!
! PROGRAM: MuellerCalculus
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief routines to generate and handle Mueller matrices and Stokes vectors for polarized light microscopy
!
!> @details Most of the routines in this module are based on the book by Collett:
!> Polarized Light: Fundamentals and Applications, E. Collett, 1993 (M. Decker, Inc)
!
!> @date 02/12/17 MDG 1.0 initial version
!> @date 09/05/19 MDG 1.1 various corrections
!--------------------------------------------------------------------------

module mod_MuellerCalculus

use mod_kinds
use mod_global

IMPLICIT NONE

type, public :: MuellerCalculus_T
  !! MuellerMatrix class definition
  private
    real(kind=dbl)                        :: M(4,4)
    character(fnlen)                      :: descriptor

  contains
    private
      procedure, pass(self) :: get_basicMuellerMatrix_
      procedure, pass(self) :: get_UniaxialReflectivities_
      procedure, pass(self) :: get_SampleMuellerMatrix_
      final :: MuellerCalculus_destructor

      generic, public :: get_basicMuellerMatrix => get_basicMuellerMatrix_
      generic, public :: get_UniaxialReflectivities => get_UniaxialReflectivities_
      generic, public :: get_SampleMuellerMatrix => get_SampleMuellerMatrix_
end type MuellerCalculus_T

! the constructor routine for this class
interface MuellerCalculus_T
    module procedure MuellerCalculus_constructor
end interface MuellerCalculus_T

contains

!--------------------------------------------------------------------------
type(MuellerCalculus_T) function MuellerCalculus_constructor( ) result(MuellerCalculus)
!DEC$ ATTRIBUTES DLLEXPORT :: MuellerCalculus_constructor
!! author: MDG
!! version: 1.0
!! date: 02/02/20
!!
!! constructor for the MuellerCalculus_T Class

IMPLICIT NONE

end function MuellerCalculus_constructor

!--------------------------------------------------------------------------
subroutine MuellerCalculus_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: MuellerCalculus_destructor
!! author: MDG
!! version: 1.0
!! date: 02/02/20
!!
!! destructor for the MuellerCalculus_T Class
    
IMPLICIT NONE
    
type(MuellerCalculus_T), INTENT(INOUT)  :: self
    
call reportDestructor('MuellerCalculus_T')
    
end subroutine MuellerCalculus_destructor

  !--------------------------------------------------------------------------
recursive subroutine get_basicMuellerMatrix_(self, MMtype)
!DEC$ ATTRIBUTES DLLEXPORT :: get_basicMuellerMatrix_
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! returns a basic 4x4 Mueller matrix by type

use mod_io

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)      :: MMtype
type(IO_T)                        :: Message

select case (MMtype)
    case (0)
        call Message%printMessage('The following basic Mueller matrix types are available:')
        call Message%printMessage('1: linear horizontal polarizer (along x)')
        call Message%printMessage('2: linear vertical polarizer (along y)')
        call Message%printMessage('3: linear polarizer at +45°')
        call Message%printMessage('4: linear polarizer at -45°')
        call Message%printMessage('5: quarter-wave plate, fast axis vertical')
        call Message%printMessage('6: quarter-wave plate, fast axis horizontal')
        call Message%printMessage('7: circular polarizer, right-handed')
        call Message%printMessage('8: circular polarizer, left-handed')
        call Message%printMessage('9: ideal mirror')
    case (1)
        self%descriptor = 'linear horizontal polarizer'
        self%M(1,1:4) = (/ 1.D0, 1.D0, 0.D0, 0.D0 /)
        self%M(2,1:4) = (/ 1.D0, 1.D0, 0.D0, 0.D0 /)
        self%M(3,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        self%M(4,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        self%M = 0.5D0 * self%M
    case (2)
        self%descriptor = 'linear vertical polarizer'
        self%M(1,1:4) = (/ 1.D0,-1.D0, 0.D0, 0.D0 /)
        self%M(2,1:4) = (/-1.D0, 1.D0, 0.D0, 0.D0 /)
        self%M(3,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        self%M(4,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        self%M = 0.5D0 * self%M
    case (3)
        self%descriptor = 'linear polarizer at +45°'
        self%M(1,1:4) = (/ 1.D0, 0.D0, 1.D0, 0.D0 /)
        self%M(2,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        self%M(3,1:4) = (/ 1.D0, 0.D0, 1.D0, 0.D0 /)
        self%M(4,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        self%M = 0.5D0 * self%M
    case (4)
        self%descriptor = 'linear polarizer at -45°'
        self%M(1,1:4) = (/ 1.D0, 0.D0,-1.D0, 0.D0 /)
        self%M(2,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        self%M(3,1:4) = (/-1.D0, 0.D0, 1.D0, 0.D0 /)
        self%M(4,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        self%M = 0.5D0 * self%M
    case (5)
        self%descriptor = 'quarter-wave plate, fast axis vertical'
        self%M(1,1:4) = (/ 1.D0, 0.D0, 0.D0, 0.D0 /)
        self%M(2,1:4) = (/ 0.D0, 1.D0, 0.D0, 0.D0 /)
        self%M(3,1:4) = (/ 0.D0, 0.D0, 0.D0,-1.D0 /)
        self%M(4,1:4) = (/ 0.D0, 0.D0, 1.D0, 0.D0 /)
    case (6)
        self%descriptor = 'quarter-wave plate, fast axis horizontal'
        self%M(1,1:4) = (/ 1.D0, 0.D0, 0.D0, 0.D0 /)
        self%M(2,1:4) = (/ 0.D0, 1.D0, 0.D0, 0.D0 /)
        self%M(3,1:4) = (/ 0.D0, 0.D0, 0.D0, 1.D0 /)
        self%M(4,1:4) = (/ 0.D0, 0.D0,-1.D0, 0.D0 /)
    case (7)
        self%descriptor = 'circular polarizer, right-handed'
        self%M(1,1:4) = (/ 1.D0, 0.D0, 0.D0, 1.D0 /)
        self%M(2,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        self%M(3,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        self%M(4,1:4) = (/ 1.D0, 0.D0, 0.D0, 1.D0 /)
        self%M = 0.5D0 * self%M
    case (8)
        self%descriptor = 'circular polarizer, left-handed'
        self%M(1,1:4) = (/ 1.D0, 0.D0, 0.D0,-1.D0 /)
        self%M(2,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        self%M(3,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        self%M(4,1:4) = (/-1.D0, 0.D0, 0.D0, 1.D0 /)
        self%M = 0.5D0 * self%M
    case (9)
        self%descriptor = 'ideal mirror'
        self%M(1,1:4) = (/ 1.D0, 0.D0, 0.D0, 0.D0 /)
        self%M(2,1:4) = (/ 0.D0, 1.D0, 0.D0, 0.D0 /)
        self%M(3,1:4) = (/ 0.D0, 0.D0,-1.D0, 0.D0 /)
        self%M(4,1:4) = (/ 0.D0, 0.D0, 0.D0,-1.D0 /)
        !self%M = 0.5D0 * self%M
    case default
end select

end subroutine get_basicMuellerMatrix_

  !--------------------------------------------------------------------------
recursive function get_UniaxialReflectivities_(self, wl, epsac, nincident, dc, beamtilt) result(rvals)
!DEC$ ATTRIBUTES DLLEXPORT :: get_UniaxialReflectivities_
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! compute the reflectivities [rss,rsp,rps,rpp] for uniaxial symmetry

use mod_io

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)  :: self
type(IO_T)                          :: Message

real(kind=dbl),INTENT(IN)           :: wl
complex(kind=dbl),INTENT(IN)        :: epsac(2)
real(kind=dbl),INTENT(IN)           :: nincident
real(kind=dbl),INTENT(IN)           :: dc(3)
real(kind=dbl),INTENT(IN)           :: beamtilt
complex(kind=dbl)                   :: rvals(4)

real(kind=dbl)                      :: k, theta1, ct, st, tt, gamma, theta
complex(kind=dbl)                   :: eps0, Deps, eps1, epsgam, no, n1, ko, k1, KK, q1, qt, qo, qroot, qe, q
complex(kind=dbl)                   :: A, B, Ap, Bp, factor, cdc(3), ke, ro, re, nne, ngam, etaO, etaE, cone

rvals = cmplx(0.D0,0.D0)
cone = cmplx(1.D0,0.D0)

! get the incident light wave number [m^-1]
! k = 2.D0 * cPi / wl 

! turn direction cosines into complex numbers
! cdc(1) = cmplx(dc(1),0.D0)
! cdc(2) = cmplx(dc(2),0.D0)
! cdc(3) = cmplx(dc(3),0.D0)

gamma = atan2(dc(2),dc(1))
theta = acos(dc(3))

etaO = sqrt(epsac(1))
etaE = sqrt(epsac(1)*epsac(2) / (epsac(1)+(epsac(2)-epsac(1))*cos(theta)**2) )

! dielectric parameters and refractive indices (complex valued !)
! Deps = epsac(2)-epsac(1)
! epsgam = epsac(1) + cdc(3)**2 * Deps
! n1 = cmplx(nincident,0.D0)
! nne = sqrt(epsac(2))
! no = sqrt(epsac(1))
! ngam = sqrt(epsgam)

if (beamtilt.eq.0.D0) then   ! we'll use the simplified expressions for the reflection coefficients
    ! source:  J. Lekner, "Normal-incidence reflection and tramsission by uniaxial crystals and crystal plates"
    ! J. Phys.: Condens. Matter 4 (1992) 1387-1398
        ro = (etaO-cone)/(etaO+cone)  ! (n1 - no)/(n1 + no)
        re = (etaE-cone)/(etaE+cone)  ! (n1*ngam - nne*no)/(n1*ngam + nne*no)
     
        rvals(1) = - (re * sin(gamma)**2 + ro * cos(gamma)**2)
        rvals(2) = - (re - ro) * sin(gamma) * cos(gamma)
        rvals(3) = - rvals(2)
        rvals(4) = ro * sin(gamma)**2 + re * cos(gamma)**2
    ! else   ! if there is a beam tilt, then we need to employ the full expressions
    ! ! source: J. Lekner, "Reflection and refraction by uniaxial crystals"
    ! ! J. Phys.: Condens. Matter 3 (1991) 6121-6133
    
    ! ! beam tilt angle
    !     theta1 = beamtilt * cPi / 180.D0
    !     ct = cos(theta1)
    !     st = sin(theta1)
    !     tt = tan(theta1)
    
    ! ! various wave numbers and wave vector components 
    !     ko = cmplx(k,0.D0) * no
    !     k1 = cmplx(k,0.D0) * nincident
    !     KK = cmplx(st,0.D0) * k
    !     q1 = cmplx(ct,0.D0) * k1
    !     qt = q1 + KK * cmplx(tt,0.D0)
    !     qo = sqrt(-KK**2+ko**2)
    !     qroot = sqrt( (epsac(1)/epsgam**2) * (k**2 * epsac(2) * epsgam - KK**2*(epsac(2) - cdc(2)**2 * Deps)) )
    !     qe = - KK * cdc(1) * cdc(3) * Deps / epsgam + qroot
    
    ! ! reflection constants
    !     A  = (qo*cdc(1)-KK*cdc(3))*(cdc(1)*(qe*ko**2+qo**2*qt)-KK*cdc(3)*(ko**2+qe*qt))
    !     Ap = (qo*cdc(1)-KK*cdc(3))*(cdc(1)*(qe*ko**2-qo**2*qt)-KK*cdc(3)*(ko**2-qe*qt))
    !     B  = (ko*cdc(2))**2 * (ko**2+qo*qt)
    !     Bp = (ko*cdc(2))**2 * (ko**2-qo*qt)
    
    ! ! and finally the four reflection parameters rss, rsp, rps, rpp
    !     factor = A*(q1+qo)+B*(q1+qe)
    
    !     rvals(1) = (A*(q1-qo)+B*(q1-qe))
    !     rvals(2) = 2.0*cdc(2)*(qo-qe)*(qo*cdc(1)+KK*cdc(3))*k1*ko**2
    !     rvals(3) = 2.0*cdc(2)*(qo-qe)*(qo*cdc(1)-KK*cdc(3))*k1*ko**2
    !     rvals(4) = - (Ap*(q1+qo)+Bp*(q1+qe))
    !     rvals = rvals/factor
end if 

end function get_UniaxialReflectivities_

  !--------------------------------------------------------------------------
recursive function get_SampleMuellerMatrix_(self, rvals) result(MM)
!DEC$ ATTRIBUTES DLLEXPORT :: get_SampleMuellerMatrix_
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! compute the reflectivities [rss,rsp,rps,rpp] for uniaxial symmetry

use mod_io

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)  :: self

real(kind=dbl)                         :: MM(4,4)
complex(kind=dbl),INTENT(IN)           :: rvals(4)
complex(kind=dbl)                      :: rpp2, rsp2, rps2, rss2, &
                                          rssrsp, rssrps, rssrpp, rsprss, rsprps, rsprpp, &
                                          rpsrss, rpsrsp, rpsrpp, rpprss, rpprsp, rpprps

type(IO_T)                          :: Message

! this is a straighforward application of the definitions in the 2012 paper
! by Letnes et al.  
!
! Calculation of the Mueller matrix for scattering of light from two-dimensional rough surfaces
! PA Letnes, AA Maradudin, T Nordam, I Simonsen
! Physical Review A 86 (3), 031803, 2012


! first get all the constants
rss2 = rvals(1) * conjg(rvals(1))
rsp2 = rvals(2) * conjg(rvals(2))
rps2 = rvals(3) * conjg(rvals(3))
rpp2 = rvals(4) * conjg(rvals(4))

rssrsp = rvals(1) * conjg(rvals(2))
rssrps = rvals(1) * conjg(rvals(3))
rssrpp = rvals(1) * conjg(rvals(4))

rsprss = rvals(2) * conjg(rvals(1))
rsprps = rvals(2) * conjg(rvals(3))
rsprpp = rvals(2) * conjg(rvals(4))

rpsrss = rvals(3) * conjg(rvals(1))
rpsrsp = rvals(3) * conjg(rvals(2))
rpsrpp = rvals(3) * conjg(rvals(4))

rpprss = rvals(4) * conjg(rvals(1))
rpprsp = rvals(4) * conjg(rvals(2))
rpprps = rvals(4) * conjg(rvals(3))

MM(1,1) = real(rpp2 + rsp2 + rps2 + rss2)
MM(1,2) = real(rpp2 - rsp2 + rps2 - rss2)
MM(1,3) = real(rsprpp + rssrps + rpprsp + rpsrss)
MM(1,4) = real(cmplx(0.D0,1.D0) * (rsprpp + rssrps - rpprsp - rpsrss))

MM(2,1) = real(rpp2 + rsp2 - rps2 - rss2)
MM(2,2) = real(rpp2 - rsp2 - rps2 + rss2)
MM(2,3) = real(rsprpp - rssrps + rpprsp - rpsrss)
MM(2,4) = real(cmplx(0.D0,1.D0) * (rsprpp - rssrps - rpprsp + rpsrss))

MM(3,1) = real(rpsrpp + rpprps + rssrsp + rsprss)
MM(3,2) = real(rpsrpp + rpprps - rssrsp - rsprss)
MM(3,3) = real(rssrpp + rsprps + rpsrsp + rpprss)
MM(3,4) = real(cmplx(0.D0,1.D0) * (rssrpp + rsprps - rpsrsp - rpprss))

MM(4,1) = real(cmplx(0.D0,-1.D0) * (rpsrpp - rpprps + rssrsp - rsprss))
MM(4,2) = real(cmplx(0.D0,-1.D0) * (rpsrpp - rpprps - rssrsp + rsprss))
MM(4,3) = real(cmplx(0.D0,-1.D0) * (rssrpp - rsprps + rpsrsp - rpprss))
MM(4,4) = real(rssrpp - rsprps - rpsrsp + rpprss)

MM = MM*0.5D0

end function get_SampleMuellerMatrix_


end module mod_MuellerCalculus
