! ###################################################################
! Copyright (c) 2013-2024, Marc De Graef selfearch Group/Carnegie Mellon University
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
!> @date 08/10/23 MDG 2.0 complet update
!--------------------------------------------------------------------------

module mod_MuellerCalculus

use mod_kinds
use mod_global
use mod_io

IMPLICIT NONE

type :: MuellerMatrixType
  real(kind=dbl)                        :: M(4,4)
  character(fnlen)                      :: descriptor
end type 

type :: StokesVectorType 
  real(kind=dbl)                        :: S(0:3)
  character(fnlen)                      :: descriptor
end type 

type, public :: MuellerCalculus_T
  !! MuellerMatrix class definition
  private
  type(MuellerMatrixType)   :: MM

  contains
    private
      procedure, pass(self) :: get_basicMuellerMatrix_
      procedure, pass(self) :: get_diattenuator_
      procedure, pass(self) :: get_rotator_
      procedure, pass(self) :: get_retarder_
      procedure, pass(self) :: rotate_MuellerMatrix_
      procedure, pass(self) :: print_MuellerMatrix_
      procedure, pass(self) :: print_StokesVector_
      procedure, pass(self) :: propagateStokesVector_
      procedure, pass(self) :: concatenateMuellerMatrices_
      procedure, pass(self) :: get_EllipticityAngle_
      procedure, pass(self) :: get_AuxiliaryAngle_
      procedure, pass(self) :: get_OrientationAngle_
      procedure, pass(self) :: get_PhaseShiftAngle_
      procedure, pass(self) :: get_Polarization_
      procedure, pass(self) :: get_Stokes_EO_
      procedure, pass(self) :: get_Stokes_AD_
      procedure, pass(self) :: get_AD_from_EO_
      procedure, pass(self) :: get_EO_from_AD_
      procedure, pass(self) :: get_UniaxialReflectivities_
      procedure, pass(self) :: get_SampleMuellerMatrix_
      final :: MuellerCalculus_destructor

      generic, public :: get_basicMuellerMatrix => get_basicMuellerMatrix_
      generic, public :: get_diattenuator => get_diattenuator_
      generic, public :: get_rotator => get_rotator_
      generic, public :: get_retarder => get_retarder_
      generic, public :: rotate_MuellerMatrix => rotate_MuellerMatrix_
      generic, public :: print_MuellerMatrix => print_MuellerMatrix_
      generic, public :: print_StokesVector => print_StokesVector_
      generic, public :: propagateStokesVector => propagateStokesVector_
      generic, public :: concatenateMuellerMatrices => concatenateMuellerMatrices_
      generic, public :: get_EllipticityAngle => get_EllipticityAngle_
      generic, public :: get_AuxiliaryAngle => get_AuxiliaryAngle_
      generic, public :: get_OrientationAngle => get_OrientationAngle_
      generic, public :: get_PhaseShiftAngle => get_PhaseShiftAngle_
      generic, public :: get_Polarization => get_Polarization_
      generic, public :: get_Stokes_EO => get_Stokes_EO_
      generic, public :: get_Stokes_AD => get_Stokes_AD_
      generic, public :: get_AD_from_EO => get_AD_from_EO_
      generic, public :: get_EO_from_AD => get_EO_from_AD_
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
recursive function get_basicMuellerMatrix_(self, MMtype) result(MM)
!DEC$ ATTRIBUTES DLLEXPORT :: get_basicMuellerMatrix_
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! returns a basic 4x4 Mueller matrix by type

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)     :: self
integer(kind=irg),INTENT(IN)                :: MMtype
type(MuellerMatrixType)                     :: MM

type(IO_T)                                  :: Message

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
        MM%descriptor = 'linear horizontal polarizer'
        MM%M(1,1:4) = (/ 1.D0, 1.D0, 0.D0, 0.D0 /)
        MM%M(2,1:4) = (/ 1.D0, 1.D0, 0.D0, 0.D0 /)
        MM%M(3,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        MM%M(4,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        MM%M = 0.5D0 * MM%M
    case (2)
        MM%descriptor = 'linear vertical polarizer'
        MM%M(1,1:4) = (/ 1.D0,-1.D0, 0.D0, 0.D0 /)
        MM%M(2,1:4) = (/-1.D0, 1.D0, 0.D0, 0.D0 /)
        MM%M(3,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        MM%M(4,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        MM%M = 0.5D0 * MM%M
    case (3)
        MM%descriptor = 'linear polarizer at +45°'
        MM%M(1,1:4) = (/ 1.D0, 0.D0, 1.D0, 0.D0 /)
        MM%M(2,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        MM%M(3,1:4) = (/ 1.D0, 0.D0, 1.D0, 0.D0 /)
        MM%M(4,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        MM%M = 0.5D0 * MM%M
    case (4)
        MM%descriptor = 'linear polarizer at -45°'
        MM%M(1,1:4) = (/ 1.D0, 0.D0,-1.D0, 0.D0 /)
        MM%M(2,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        MM%M(3,1:4) = (/-1.D0, 0.D0, 1.D0, 0.D0 /)
        MM%M(4,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        MM%M = 0.5D0 * MM%M
    case (5)
        MM%descriptor = 'quarter-wave plate, fast axis vertical'
        MM%M(1,1:4) = (/ 1.D0, 0.D0, 0.D0, 0.D0 /)
        MM%M(2,1:4) = (/ 0.D0, 1.D0, 0.D0, 0.D0 /)
        MM%M(3,1:4) = (/ 0.D0, 0.D0, 0.D0,-1.D0 /)
        MM%M(4,1:4) = (/ 0.D0, 0.D0, 1.D0, 0.D0 /)
    case (6)
        MM%descriptor = 'quarter-wave plate, fast axis horizontal'
        MM%M(1,1:4) = (/ 1.D0, 0.D0, 0.D0, 0.D0 /)
        MM%M(2,1:4) = (/ 0.D0, 1.D0, 0.D0, 0.D0 /)
        MM%M(3,1:4) = (/ 0.D0, 0.D0, 0.D0, 1.D0 /)
        MM%M(4,1:4) = (/ 0.D0, 0.D0,-1.D0, 0.D0 /)
    case (7)
        MM%descriptor = 'circular polarizer, right-handed'
        MM%M(1,1:4) = (/ 1.D0, 0.D0, 0.D0, 1.D0 /)
        MM%M(2,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        MM%M(3,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        MM%M(4,1:4) = (/ 1.D0, 0.D0, 0.D0, 1.D0 /)
        MM%M = 0.5D0 * MM%M
    case (8)
        MM%descriptor = 'circular polarizer, left-handed'
        MM%M(1,1:4) = (/ 1.D0, 0.D0, 0.D0,-1.D0 /)
        MM%M(2,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        MM%M(3,1:4) = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        MM%M(4,1:4) = (/-1.D0, 0.D0, 0.D0, 1.D0 /)
        MM%M = 0.5D0 * MM%M
    case (9)
        MM%descriptor = 'ideal mirror'
        MM%M(1,1:4) = (/ 1.D0, 0.D0, 0.D0, 0.D0 /)
        MM%M(2,1:4) = (/ 0.D0, 1.D0, 0.D0, 0.D0 /)
        MM%M(3,1:4) = (/ 0.D0, 0.D0,-1.D0, 0.D0 /)
        MM%M(4,1:4) = (/ 0.D0, 0.D0, 0.D0,-1.D0 /)
    case default
end select

end function get_basicMuellerMatrix_

 !--------------------------------------------------------------------------
recursive function get_diattenuator_(self, px, py, polar) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: get_diattenuator_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! returns a di-attenuator 4x4 Mueller matrix

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)     :: self
real(kind=dbl),INTENT(IN)                   :: px
real(kind=dbl),INTENT(IN)                   :: py
logical,OPTIONAL,INTENT(IN)                 :: polar

type(IO_T)                                  :: Message

type(MuellerMatrixType)                     :: res
logical                                     :: usepolar

! initialize a Mueller matrix for a diattenuator
res%descriptor = 'diattenuator'

usepolar = .FALSE.
if (present(polar)) then
    if (polar.eqv..TRUE.) usepolar = .TRUE.
end if 

if (usepolar) then
    if ((px.lt.0.D0).or.(px.gt.1.D0)) then
        call Message%printError('MC_get_diattenuator','attenuation magnitude must lie in range [0,1]')
    end if 
    res%M(1,1:4) = (/ 1.D0, cos(2.D0*py), 0.D0, 0.D0 /)
    res%M(2,1:4) = (/ cos(2.D0*py), 1.D0, 0.D0, 0.D0 /)
    res%M(3,1:4) = (/ 0.D0, 0.D0, sin(2.D0*py), 0.D0 /)
    res%M(4,1:4) = (/ 0.D0, 0.D0, 0.D0, sin(2.D0*py) /)
    res%M = 0.5D0*px*px*res%M
else
    if ((minval((/ px, py /)).lt.0.D0).or.(maxval((/px, py/)).gt.1.D0)) then
        call Message%printError('MC_get_diattenuator','attenuation factors must lie in range [0,1]')
    end if 
    res%M(1,1:4) = (/ px*px+py*py, px*px-py*py, 0.D0, 0.D0 /)
    res%M(2,1:4) = (/ px*px-py*py, px*px+py*py, 0.D0, 0.D0 /)
    res%M(3,1:4) = (/ 0.D0, 0.D0, 2.D0*px*py, 0.D0 /)
    res%M(4,1:4) = (/ 0.D0, 0.D0, 0.D0, 2.D0*px*py /)
    res%M = 0.5D0*res%M
end if 
    
end function get_diattenuator_

!--------------------------------------------------------------------------
recursive function get_rotator_(self, theta) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: get_rotator_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! returns a rotator 4x4 Mueller matrix

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)     :: self
real(kind=dbl),INTENT(IN)                   :: theta
type(MuellerMatrixType)                     :: res

real(kind=dbl)                              :: ct, st

ct = cos(2.D0*theta)
st = sin(2.D0*theta)

! initialize a Mueller matrix for a rotator
res%descriptor = 'rotator'
res%M(1,1:4) = (/ 1.D0, 0.D0, 0.D0, 0.D0 /)
res%M(2,1:4) = (/ 0.D0, ct, st, 0.D0 /)
res%M(3,1:4) = (/ 0.D0, -st, ct, 0.D0 /)
res%M(4,1:4) = (/ 0.D0, 0.D0, 0.D0, 1.D0 /)

end function get_rotator_

!--------------------------------------------------------------------------
recursive function get_retarder_(self, phi) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: get_retarder_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! returns a retarder 4x4 Mueller matrix

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)     :: self
real(kind=dbl),INTENT(IN)                   :: phi
type(MuellerMatrixType)                     :: res

real(kind=dbl)                              :: cp, sp

cp = cos(phi)
sp = sin(phi)

! initialize a Mueller matrix for a retarder
res%descriptor = 'retarder'
res%M(1,1:4) = (/ 1.D0, 0.D0, 0.D0, 0.D0 /)
res%M(2,1:4) = (/ 0.D0, 1.D0, 0.D0, 0.D0 /)
res%M(3,1:4) = (/ 0.D0, 0.D0, cp, -sp /)
res%M(4,1:4) = (/ 0.D0, 0.D0, sp, cp /)

end function get_retarder_

!--------------------------------------------------------------------------
recursive function rotate_MuellerMatrix_(self, MM, theta, normalincidence) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: rotate_MuellerMatrix_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! rotate a 4x4 Mueller matrix

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)     :: self
type(MuellerMatrixType),INTENT(IN)          :: MM 
real(kind=dbl),INTENT(IN)                   :: theta
logical,INTENT(IN),OPTIONAL                 :: normalincidence
type(MuellerMatrixType)                     :: res

type(MuellerMatrixType)                     :: Mrot 
logical                                     :: normal 

normal = .FALSE.
if (present(normalincidence)) normal = .TRUE.

! initialize the output Mueller matrix descriptor
res%descriptor = trim(MM%descriptor)//'-rotated'

Mrot = self%get_rotator(theta)
if (normal.eqv..FALSE.) then 
    res%M = matmul(transpose(Mrot%M), matmul(MM%M, Mrot%M))
else 
    res%M = matmul(Mrot%M, matmul(MM%M, Mrot%M))
end if

end function rotate_MuellerMatrix_

!--------------------------------------------------------------------------
recursive subroutine print_MuellerMatrix_(self, MM)
!DEC$ ATTRIBUTES DLLEXPORT :: print_MuellerMatrix_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! print a 4x4 Mueller matrix

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)     :: self
type(MuellerMatrixType),INTENT(IN)          :: MM

type(IO_T)                                  :: Message

real(kind=dbl)                              :: io_double(4)
integer(kind=irg)                           :: i

call Message%printMessage('Mueller Matrix Type : '//trim(MM%descriptor))

do i=1,4 
    io_double(1:4) = MM%M(i,1:4)
    call Message%WriteValue('  --> ',io_double,4)
end do

end subroutine print_MuellerMatrix_

!--------------------------------------------------------------------------
recursive subroutine print_StokesVector_(self, S)
!DEC$ ATTRIBUTES DLLEXPORT :: print_StokesVector_
!! author: MDG 
!! version: 1.0 
!! date: 08/13/23
!!
!! print a Stokes Vector

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)     :: self
type(StokesVectorType),INTENT(IN)           :: S

type(IO_T)                                  :: Message

real(kind=dbl)                              :: io_double(4)
integer(kind=irg)                           :: i

call Message%printMessage('Stokes Vector descriptor : '//trim(S%descriptor))

io_double(1:4) = S%S(0:3)
call Message%WriteValue('  --> ',io_double,4)

end subroutine print_StokesVector_

!--------------------------------------------------------------------------
recursive function propagateStokesVector_(self, MM, SV, descriptor) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: propagateStokesVector_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! multiplies a Stokes vector by a Mueller matrix

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)     :: self
type(MuellerMatrixType),INTENT(IN)          :: MM
type(StokesVectorType),INTENT(IN)           :: SV
character(fnlen),INTENT(IN)                 :: descriptor
type(StokesVectorType)                      :: res

res%S = matmul(MM%M, SV%S)
res%descriptor = trim(descriptor)

end function propagateStokesVector_

!--------------------------------------------------------------------------
recursive function concatenateMuellerMatrices_(self, MM1, MM2) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: concatenateMuellerMatrices_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! multiplies a Mueller matrix M1 by M2, in the order  M2 x M1 
!!
!! note that MM1 is earlier in the optical path than MM2

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)     :: self
type(MuellerMatrixType),INTENT(IN)          :: MM1
type(MuellerMatrixType),INTENT(IN)          :: MM2
type(MuellerMatrixType)                     :: res

res%M = matmul(MM2%M, MM1%M)
res%descriptor = trim(MM1%descriptor)//'->'//trim(MM2%descriptor)

end function concatenateMuellerMatrices_

!--------------------------------------------------------------------------
recursive function get_EllipticityAngle_(self, SV) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: get_EllipticityAngle_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! extracts the ellipticity angle from a Stokes vector 

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)     :: self
type(StokesVectorType),INTENT(IN)           :: SV
real(kind=dbl)                              :: res

type(IO_T)                                  :: Message 

real(kind=dbl)                              :: p4, io_double(1)

p4 = cPi * 0.25D0

res = 0.5D0 * asin(SV%S(3)/SV%S(0))

if (abs(res).gt.p4) then
    io_double(1) = res
    call Message%WriteValue('Ellipticity angle = ',io_double,1)
    call Message%printError('MC_get_EllipticityAngle','Ellipticity angle does not lie in range [-pi/4,pi/4]')
end if

end function get_EllipticityAngle_

!--------------------------------------------------------------------------
recursive function get_OrientationAngle_(self, SV) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: get_OrientationAngle_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! extracts the polarization ellipse orientation angle from a Stokes vector 

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)     :: self
type(StokesVectorType),INTENT(IN)           :: SV
real(kind=dbl)                              :: res

res = 0.5D0 * atan2(SV%S(2),SV%S(1))

res = mod(res+2.D0*cPi,cPi)

end function get_OrientationAngle_

!--------------------------------------------------------------------------
recursive function get_AuxiliaryAngle_(self, SV) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: get_AuxiliaryAngle_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! extracts the auxiliary angle from a Stokes vector

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)     :: self
type(StokesVectorType),INTENT(IN)           :: SV
real(kind=dbl)                              :: res

real(kind=dbl)                              :: psi, chi, alpha, delta

chi = self%get_EllipticityAngle(SV)
psi = self%get_OrientationAngle(SV)
call self%get_AD_from_EO(chi, psi, alpha, delta)

res = alpha

end function get_AuxiliaryAngle_

!--------------------------------------------------------------------------
recursive function get_PhaseShiftAngle_(self, SV) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: get_PhaseShiftAngle_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! extracts the phase shift angle from a Stokes vector

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)     :: self
type(StokesVectorType),INTENT(IN)           :: SV
real(kind=dbl)                              :: res

real(kind=dbl)                              :: psi, chi, alpha, delta

chi = self%get_EllipticityAngle(SV)
psi = self%get_OrientationAngle(SV)
call self%get_AD_from_EO(chi, psi, alpha, delta)

res = delta

end function get_PhaseShiftAngle_

!--------------------------------------------------------------------------
recursive function get_Polarization_(self, SV) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Polarization_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! extracts the polarization from a Stokes vector

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)     :: self
type(StokesVectorType),INTENT(IN)           :: SV
real(kind=dbl)                              :: res

type(IO_T)                                  :: Message 

if (SV%S(0).eq.0.D0) then
    call Message%printError('MC_get_Polarization','Total intensity in Stokes Vector is zero')
end if

res = dsqrt(sum(SV%S(1:3)**2)) / SV%S(0)

end function get_Polarization_

!--------------------------------------------------------------------------
recursive function get_Stokes_EO_(self, chi, psi, descriptor) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Stokes_EO_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! generate a Stokes vector for a given Ellipticity and Orientation angle

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)     :: self
real(kind=dbl),INTENT(IN)                   :: chi
real(kind=dbl),INTENT(IN)                   :: psi
character(fnlen),INTENT(IN)                 :: descriptor
type(StokesVectorType)                      :: res

real(kind=dbl)                              :: cp, sp, cc, sc

cp = cos(2.D0*psi)
sp = sin(2.D0*psi)
cc = cos(2.D0*chi)
sc = sin(2.D0*chi)

res%descriptor = trim(descriptor)
res%S = (/ 1.D0, cc*cp, cc*sp, sc /)

end function get_Stokes_EO_

!--------------------------------------------------------------------------
recursive function get_Stokes_AD_(self, alpha, delta, descriptor) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Stokes_AD_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! generate a Stokes vector for a given auxiliary and phase shift angle

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)     :: self
real(kind=dbl),INTENT(IN)                   :: alpha
real(kind=dbl),INTENT(IN)                   :: delta 
character(fnlen),INTENT(IN)                 :: descriptor
type(StokesVectorType)                      :: res

real(kind=dbl)                              :: ca, sa, cd, sd

ca = cos(2.D0*alpha)
sa = sin(2.D0*alpha)
cd = cos(delta)
sd = sin(delta)

res%descriptor = trim(descriptor)
res%S = (/ 1.D0, ca, sa*cd, sa*sd /)

end function get_Stokes_AD_

!--------------------------------------------------------------------------
recursive subroutine get_AD_from_EO_(self, chi, psi, alpha, delta)
!DEC$ ATTRIBUTES DLLEXPORT :: get_AD_from_EO_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! convert auxiliary and phase shift angle to ellipticity and orientation angles

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)     :: self
real(kind=dbl),INTENT(IN)                   :: chi
real(kind=dbl),INTENT(IN)                   :: psi
real(kind=dbl),INTENT(OUT)                  :: alpha
real(kind=dbl),INTENT(OUT)                  :: delta 

real(kind=dbl)                              :: sc, tp, cp, cc, p2, p4, st, tt, ss, ct, sa

p2 = cPi*0.5D0
p4 = cPi*0.25D0

cc = cos(2.D0*chi)
sc = sin(2.D0*chi)
tp = tan(2.D0*psi)
cp = cos(2.D0*psi)

st = dsqrt(sc*sc+tp*tp)
tt = dsqrt(1.D0+tp*tp)
ss = dsqrt(1.D0-sc*sc)

ct = cos(2.D0*chi) * tan(2.D0*psi)
sa = sc/abs(cp)

! get alpha
if (abs(psi-p2).ge.p4) then
  alpha = 0.5D0 * atan2(st/tt,ss/tt)
else
  alpha = 0.5D0 * (cPi - atan2(st/tt,ss/tt))
end if

! get delta, such that there is only one cut in the delta surface for chi=0, psi<pi/2
if (abs(psi-p2).lt.p4) then
    delta = atan2(-sa/st,ct/st)-cPi
else
    delta = atan2(sa/st,ct/st)
    if (chi.gt.0.D0) delta = delta - 2.0D0*cPi
end if

end subroutine get_AD_from_EO_

!--------------------------------------------------------------------------
recursive subroutine get_EO_from_AD_(self, alpha, delta, chi, psi)
!DEC$ ATTRIBUTES DLLEXPORT :: get_EO_from_AD_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! convert ellipticity and orientation angles to auxiliary and phase shift angles

IMPLICIT NONE

class(MuellerCalculus_T), INTENT(INOUT)     :: self
real(kind=dbl),INTENT(IN)                   :: alpha
real(kind=dbl),INTENT(IN)                   :: delta 
real(kind=dbl),INTENT(OUT)                  :: chi
real(kind=dbl),INTENT(OUT)                  :: psi

type(IO_T)                                  :: Message 

real(kind=dbl)                              :: p2, p4

p2 = cPi * 0.5D0
p4 = cPi * 0.25D0

chi = 0.5D0 * asin ( sin(2.D0 * alpha) * sin(delta))

if (delta.le.p2) then
    psi = 0.5D0 * atan(cos(delta) * tan(2.D0 * alpha))
else
    psi = cPi - 0.5D0 * atan(cos(delta) * tan(2.D0 * alpha))
end if

! make sure chi falls in the range [-pi/4,pi/4]
if (abs(chi).gt.p4) then
    call Message%printError('MC_get_EO_from_AD','ellipticity angle must be in interval [-pi/4,pi/4]')
end if

! make sure psi falls in the range [0,pi]
if (psi.lt.0.D0) psi = psi + cPi

end subroutine get_EO_from_AD_

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
