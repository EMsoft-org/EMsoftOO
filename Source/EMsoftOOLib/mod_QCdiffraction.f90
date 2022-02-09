! ###################################################################
! Copyright (c) 2014-2022, Marc De Graef Research Group/Carnegie Mellon University
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
module mod_QCdiffraction

use mod_kinds
use mod_global

IMPLICIT NONE

private

! The parameters in gnode are computed by CalcQCUcg
type, public :: QCgnode
  character(2)        :: method   ! computation method (WK = Weickenmeier-Kohl, DT = Doyle-Turner/Smith-Burge, XR for XRD)
  logical             :: absorption ! is absorption included or not ?
  integer(kind=irg)   :: hkl(3)   ! Miller indices
  real(kind=sgl)      :: xg, &    ! extinction distance [nm]
                         xgp, &   ! absorption length [nm]
                         ar, &    ! aborption ratio
                         g, &     ! length of reciprocal lattice vectors [nm^-1]
                         Vmod,Vpmod, & ! modulus of Vg and Vgprime [V]
                         Umod,Upmod, & ! modulus of Ug and Ugprime [nm^-2]
                         Vphase,Vpphase ! phase factors of Vg and Vgprime [rad]
  complex(kind=sgl)   :: Ucg, &   ! scaled potential Fourier coefficient [nm^-2]
                         Vg, &    ! potential Fourier coefficient [V]
                         qg       ! interaction parameter for Darwin-Howie-Whelan equations [nm^-1]
end type QCgnode


type, public :: DynType
  real(kind=sgl)                   :: WV(3)           ! wave vector expressed in reciprocal frame
  real(kind=sgl)                   :: FN(3)           ! Foil normal in reciprocal frame
  real(kind=sgl)                   :: Upz             ! U'_0 normal absorption parameter
! complex(kind=dbl),allocatable    :: W(:), &         ! eigenvalue vector for Bloch wave method
!                                   CG(:,:), &        ! eigenvector matrix
!                                   alpha(:), &       ! excitation amplitude vector
!                                   DHWMz(:,:),&      ! Darwin-Howie-Whelan matrix
  complex(kind=dbl),allocatable    :: DynMat(:,:)     ! dynamical matrix
!                                   DynMat0(:,:), &   ! dynamical matrix (for programs that need two or more of them)
!                                   DynMat1(:,:), &   ! dynamical matrix (for programs that need two or more of them)
!                                   DynMat2(:,:), &   ! dynamical matrix (for programs that need two or more of them)
!                                   DynMat3(:,:), &   ! dynamical matrix (for programs that need two or more of them)
!                                   phiz(:),Az(:,:)   ! used for Taylor expansion of scattering matrix
end type DynType

type, public :: BetheParameterType
        real(kind=sgl)                 :: c1 = 8.0_sgl         ! changed from 8 and 12 for a test on 8/14/15
        real(kind=sgl)                 :: c2 = 50.0_sgl
        real(kind=sgl)                 :: c3 = 100.0_sgl
        real(kind=sgl)                 :: sgdbdiff = 1.00_sgl    ! changed from 0.05 on 08/14/15 by MDG
        real(kind=sgl)                 :: weakcutoff = 0.0_sgl
        real(kind=sgl)                 :: cutoff = 0.0_sgl
        real(kind=sgl)                 :: sgcutoff = 0.0_sgl
        integer(kind=irg)              :: nns
        integer(kind=irg)              :: nnw
        integer(kind=irg)              :: minweak
        integer(kind=irg)              :: minstrong
        integer(kind=irg)              :: maxweak
        integer(kind=irg)              :: maxstrong
        integer(kind=irg)              :: totweak
        integer(kind=irg)              :: totstrong
        integer(kind=irg),allocatable  :: weaklist(:)
        integer(kind=irg),allocatable  :: stronglist(:)
        integer(kind=irg),allocatable  :: weakhkl(:,:)
        integer(kind=irg),allocatable  :: stronghkl(:,:)
        real(kind=sgl),allocatable     :: weaksg(:)
        real(kind=sgl),allocatable     :: strongsg(:)
        integer(kind=irg),allocatable  :: strongID(:)
        integer(kind=sgl),allocatable  :: reflistindex(:)              ! used to map strong reflections onto the original reflist
        integer(kind=sgl),allocatable  :: weakreflistindex(:)          ! used to map weak reflections onto the original reflist
end type BetheParameterType

! the diffraction class that deals with the above parameter types
type, public :: QCDiffraction_T
private

    real(kind=dbl)                  :: voltage          ! voltage always in keV !
    real(kind=dbl)                  :: Lambda           ! wave length in nm
    real(kind=dbl)                  :: Relcor
    real(kind=dbl)                  :: Sigma
    real(kind=dbl)                  :: Psihat
    real(kind=dbl)                  :: V0mod = 0.D0     ! used to compute refraction correction
    type(QCgnode)                   :: QCrlp            ! default variable for reciprocal lattice point
    type(DynType)                   :: Dyn              ! dynamical scattering parameters
    type(BetheParameterType)        :: BetheParameters
    real(kind=sgl),allocatable      :: scatfacg(:)
    complex(kind=sgl),allocatable   :: scatfac(:,:)
    integer(kind=irg)               :: numscatfac
    complex(kind=sgl),allocatable   :: LUT(:), SghLUT(:,:,:,:)
    complex(kind=sgl),allocatable   :: LUTqg(:)
    logical,allocatable             :: dbdiff(:)

contains
private

    procedure, pass(self) :: CalcQCUcg_
    procedure, pass(self) :: QCCalcWaveLength_
    procedure, pass(self) :: CalcsgQC_
    procedure, pass(self) :: getLUT_
    procedure, pass(self) :: getLUTqg_
    procedure, pass(self) :: getBetheParameter_
    procedure, pass(self) :: getWaveLength_
    procedure, pass(self) :: getshapeLUT_
    procedure, pass(self) :: setV_
    procedure, pass(self) :: getV_
    procedure, pass(self) :: allocateLUT_
    procedure, pass(self) :: setLUT_
    procedure, pass(self) :: setLUTqg_
    procedure, pass(self) :: setdbdiff_

    ! procedure, pass(self) :: CalcDiffAngle_
    ! procedure, pass(self) :: PreCalcFSCATT_
    ! procedure, pass(self) :: CalcsgSingle_
    ! procedure, pass(self) :: Set_Bethe_Parameters_
    ! procedure, pass(self) :: writeBetheparameterNameList_
    ! procedure, pass(self) :: BWsolve_
    ! procedure, pass(self) :: getScatfac_
    ! procedure, pass(self) :: getV0mod_
    ! procedure, pass(self) :: getRelcor_
    ! procedure, pass(self) :: getSigma_
    ! procedure, pass(self) :: getPsihat_
    ! procedure, pass(self) :: getrlp_
    ! procedure, pass(self) :: getSghLUT_
    ! procedure, pass(self) :: getdbdiff_
    ! procedure, pass(self) :: setrlpmethod_
    ! procedure, pass(self) :: Initialize_SghLUT_
    ! procedure, pass(self) :: preCalcSgh_
    ! procedure, pass(self) :: preCalcSghECCI_
    ! procedure, pass(self) :: Printrlp_
    final :: QCDiffraction_destructor


    generic, public :: CalcQCUcg => CalcQCUcg_
    generic, public :: QCCalcWaveLength => QCCalcWaveLength_
    generic, public :: CalcsgQC => CalcsgQC_
    generic, public :: getLUT => getLUT_
    generic, public :: getLUTqg => getLUTqg_
    generic, public :: getBetheParameter => getBetheParameter_
    generic, public :: getWaveLength => getWaveLength_
    generic, public :: getshapeLUT => getshapeLUT_
    generic, public :: setV => setV_
    generic, public :: getV => getV_
    generic, public :: allocateLUT => allocateLUT_
    generic, public :: setLUT => setLUT_
    generic, public :: setLUTqg => setLUTqg_
    generic, public :: setdbdiff => setdbdiff_

    ! generic, public :: CalcDiffAngle => CalcDiffAngle_
    ! generic, public :: PreCalcFSCATT => PreCalcFSCATT_
    ! generic, public :: getScatfac => getScatfac_
    ! generic, public :: SetBetheParameters => Set_Bethe_Parameters_
    ! generic, public :: writeBetheparameterNameList => writeBetheparameterNameList_
    ! generic, public :: BWsolve => BWsolve_
    ! generic, public :: setrlpmethod => setrlpmethod_
    ! generic, public :: getV0mod => getV0mod_
    ! generic, public :: getRelcor => getRelcor_
    ! generic, public :: getSigma => getSigma_
    ! generic, public :: getPsihat => getPsihat_
    ! generic, public :: getrlp => getrlp_
    ! generic, public :: getSghLUT => getSghLUT_
    ! generic, public :: getdbdiff => getdbdiff_
    ! generic, public :: Initialize_SghLUT => Initialize_SghLUT_
    ! generic, public :: preCalcSgh => preCalcSgh_
    ! generic, public :: preCalcSghECCI => preCalcSghECCI_
    ! generic, public :: Printrlp => Printrlp_

end type QCDiffraction_T

! the constructor routine for this class
  interface QCDiffraction_T
    module procedure QCDiffraction_constructor
  end interface QCDiffraction_T

contains

!--------------------------------------------------------------------------
type(QCDiffraction_T) function QCDiffraction_constructor( voltage, QCcell, QCSG, verbose ) result(QCDiff)
!DEC$ ATTRIBUTES DLLEXPORT :: QCDiffraction_constructor
  !! author: MDG
  !! version: 1.0
  !! date: 02/06/22
  !!
  !! constructor for the QCDiffraction Class

use mod_QCcrystallography
use mod_QCsymmetry

IMPLICIT NONE

real(kind=dbl), INTENT(IN)          :: voltage    ! in kV
type(QCCell_T), INTENT(INOUT)       :: QCcell
type(QCspacegroup_T),INTENT(INOUT)  :: QCSG
logical, INTENT(IN), OPTIONAL       :: verbose

QCDiff%voltage = voltage
QCDiff%QCrlp%method = 'WK'

if (present(verbose)) then
  if (verbose.eqv..TRUE.) call QCDiff%QCCalcWaveLength_( QCcell, QCSG, verbose )
else
  call QCDiff%QCCalcWaveLength_( QCcell, QCSG )
end if

end function QCDiffraction_constructor

!--------------------------------------------------------------------------
subroutine QCDiffraction_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: QCDiffraction_destructor
!! author: MDG
!! version: 1.0
!! date: 02/06/22
!!
!! destructor for the QCDiffraction_T Class

IMPLICIT NONE

type(QCDiffraction_T), INTENT(INOUT)  :: self

! call reportDestructor('Diffraction_T')

if (allocated(self%scatfacg)) deallocate(self%scatfacg)
if (allocated(self%scatfac)) deallocate(self%scatfac)
if (allocated(self%LUT)) deallocate(self%LUT)
if (allocated(self%SghLUT)) deallocate(self%SghLUT)
if (allocated(self%LUTqg)) deallocate(self%LUTqg)
if (allocated(self%dbdiff)) deallocate(self%dbdiff)

end subroutine QCDiffraction_destructor

!--------------------------------------------------------------------------
recursive subroutine QCCalcWaveLength_(self, QCcell, QCSG, verbose)
!DEC$ ATTRIBUTES DLLEXPORT :: QCCalcWaveLength_
  !! author: MDG
  !! version: 1.0
  !! date: 02/06/22
  !!
  !! compute the relativistic electron wavelength
  !!
  !! These quantities are computed in double precision because of the
  !! wide range of magnitudes.  If a crystal structure has been defined
  !! then the gamma*V_0 term is added to correct for refraction.

use mod_QCsymmetry
use mod_QCcrystallography
use mod_io

IMPLICIT NONE

class(QCDiffraction_T),INTENT(INOUT)    :: self
class(QCCell_T),INTENT(INOUT)           :: QCcell
type(QCspacegroup_T),INTENT(INOUT)      :: QCSG
logical,INTENT(IN),OPTIONAL             :: verbose

type(IO_T)                              :: Message

real(kind=dbl)                          :: temp1,temp2, oi_real(1), Vmod, Vpmod, xig, xgp, io_real(1)
integer(kind=irg)                       :: io_int(1)
complex(kind=dbl)                       :: Ucg, qg
integer(kind=irg),allocatable           :: hkl(:)

temp1 = 1.0D+9*cPlanck/dsqrt(2.D0*cRestmass*cCharge)
temp2 = cCharge*0.5D0*self%voltage*1000.D0/cRestmass/(cLight**2)

! relativistic correction factor (known as gamma)
self%Relcor = 1.0D0+2.0D0*temp2

! relativistic acceleration voltage
self%Psihat = self%voltage*(1.D0+temp2)*1000.D0

! compute V_0 and add it to mPsihat (corrected by mRelcor)
! this should have been done in the calling program
!     call cell%CalcPositions('v')
select type (QCcell)
  class is (QCcell_icosahedral_T)
    allocate( hkl(6) )
    hkl=(/ 0, 0, 0, 0, 0, 0/)
    Ucg = self%calcQCUcg_(QCcell,QCSG,hkl,qg,Vmod,Vpmod,xig,xgp)
  class is (QCcell_axial_T)
    allocate( hkl(5) )
    hkl=(/ 0, 0, 0, 0, 0/)
    Ucg = self%calcQCUcg_(QCcell,QCSG,hkl,qg,Vmod,Vpmod,xig,xgp) 
end select
deallocate( hkl )

! wave length, corrected for refraction
self%Psihat = self%Psihat + dble(self%V0mod)
self%Lambda = temp1/dsqrt(self%Psihat)

! interaction constant sigma
self%Sigma = 2.D0*cPi*cRestmass*self%Relcor*cCharge*self%Lambda
self%Sigma = 1.0D-18*self%Sigma/cPlanck**2

if (present(verbose)) then
  if (verbose) then
    oi_real(1) = self%V0mod
    call Message%WriteValue('', oi_real, 1,"(/' Mean inner potential [V]                ',E10.4)")
    oi_real(1) = self%Relcor
    call Message%WriteValue('', oi_real, 1,"(' Relativistic correction factor [gamma]  ',E10.4)")
    oi_real(1) = self%Psihat
    call Message%WriteValue('', oi_real, 1,"(' Relativistic Accelerating Potential [V] ',E10.4)")
    oi_real(1) = self%Lambda
    call Message%WriteValue('', oi_real, 1,"(' Electron Wavelength [nm]                ',E10.4)")
    oi_real(1) = self%Sigma
    call Message%WriteValue('', oi_real, 1,"(' Interaction constant [V nm]^(-1)        ',E10.4)")
  end if
end if

end subroutine QCCalcWaveLength_
  
!--------------------------------------------------------------------------
recursive function CalcQCUcg_(self, QCcell, QCSG, hkl, qg, Vmod, Vpmod, xig, xgp) result(Ucg)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcQCUcg_
  !! author: MDG
  !! version: 1.0
  !! date: 01/26/20
  !!
  !! compute the interaction coupling constants Ucg and qg for this reflection
  !!
  !! For now, in this first version of the QC codes, we use the simplistic
  !! primitive hypercubic structure model from Veit Elser's paper; we'll need to expand 
  !! on this routine to be able to simulate more realistic QC structures.

use mod_QCcrystallography
use mod_QCsymmetry
use mod_others

IMPLICIT NONE

class(QCDiffraction_T),INTENT(INOUT)  :: self
class(QCCell_T),INTENT(INOUT)         :: QCcell
type(QCspacegroup_T),INTENT(INOUT)    :: QCSG
integer(kind=irg),INTENT(IN)          :: hkl(:)
complex(kind=dbl),INTENT(OUT)         :: qg
real(kind=dbl),INTENT(OUT)            :: Vmod
real(kind=dbl),INTENT(OUT)            :: Vpmod
real(kind=dbl),INTENT(OUT)            :: xig
real(kind=dbl),INTENT(OUT)            :: xgp
complex(kind=dbl)                     :: Ucg 

integer(kind=irg)                     :: i,j,absflg,m,ii
real(kind=sgl)                        :: s,twopi,arg,swk,dwwk,pref,ul,pre,preg,sct,fs,fsp, g, go, p, q, &
                                         Umod, Upmod, Vphase, Vpphase, multiplicity
complex(kind=sgl)                     :: ff,gg,sf,p1,pp
complex(kind=sgl)                     :: czero
logical                               :: accflg, dwflg
character(2)                          :: smb

twopi = sngl(2.D0*cPi)
czero = cmplx(0.0,0.0)

! we'll only use the Weickenmeier-Kohl scattering parameters here since this
! allows us to also get the absorption part of the Fourier coefficients
! some of this code is directly adapted from the CalcUcg routine in diffraction.f90

! first we need to get the scattering parameter
! compute the scattering parameter s^2=(g/2)^2
if (sum(hkl**2).eq.0) then 
  s = 0.0
  g = 0.0
  go = 0.0
else
  select type (QCcell)
    class is (QCcell_icosahedral_T)
      g = QCcell%getvectorLength( hkl, 'P', 'r' )
    class is (QCcell_axial_T)
      g = QCcell%getvectorLength( hkl, 'P', 'r' )
  end select
  ! s = (0.50*g)**2 ! done later with proper scaling for WK
end if

! The Weickenmeier-Kohl (WK) subroutine works in Angstrom, and also
! scales reciprocal space by a factor of 2*pi;  this scaling
! is accomplished by changing the g-value in nm^{-1} by a
! scaling factor swk = 2*pi/10, to go from book units to WK units.
!
! A similar scaling must be performed on the Debye Waller factor;
! the book defines it as exp(-Bs^2), with s in [nm^2]; WK define
! it in A^2, and with a scaled reciprocal space.  The conversion
! factor dwwk = 100.0/8*pi^2
!
! To go from the standard B factor in [nm^2] to ul^2 in A^2,
! ul = sqrt(B*dwwk)
swk = 0.1*twopi
dwwk = 100.0/(8.0*cPi**2)

! properly scale the scattering parameter
s = (0.5*g*swk)**2

! let fscatt perform the relativistic corrections for f_g and fprime_g
accflg = .TRUE.

! include absorption ?
absflg = 3

! always include Debye-Waller factor
dwflg  = .TRUE.

! compute the scaling prefactors
! pref contains A to nm conversion, and divides by 4pi
! WE NEED TO FIGURE OUT WHAT TO DO FOR THE UNIT CELL VOLUME !!!  (this may not be correct...)
pref = 0.04787801/(4.0*cPi)/QCcell%get_vol()

! preg is used to go from V to U, remembering that gamma is already
! included in the output from fscatt
preg = 2.0 * sngl(cRestmass*cCharge/cPlanck**2)*1.0E-18
pre = pref * preg

! initialize the real and imaginary parts of the structure factor
ff=czero
gg=czero

! loop over all atoms 
select type (QCcell)
  class is (QCcell_icosahedral_T)
    do m=1,QCcell%get_ATOM_ntype()
      ul = sqrt(QCcell%ATOM_pos(m,8)*dwwk + QCcell%ATOM_pos(m,9)*dwwk)
      j  = QCcell%ATOM_type(m)
      sf = FSCATT(s,ul,j,smb,sngl(self%voltage),absflg,accflg,dwflg) * QCcell%ATOM_pos(m,7)
      pp = QCcell%ShapeTransformTriacontahedron(QCSG, hkl, m)

    ! loop over all atoms in the orbit
      p1 = czero
      do j = 1,QCcell%numat(m)
        arg = twopi*sum(float(hkl(1:6))*QCcell%apos(m,j,1:6))
        p1  = p1 + exp(cmplx(0.0,-arg))
      end do

      ff = ff + p1*real(sf)  * abs(real(pp))
      gg = gg + p1*aimag(sf) * abs(real(pp))
    end do
  class is (QCcell_axial_T)
    do m = 1,QCcell%get_ATOM_ntype()
      ul = sqrt(QCcell%ATOM_pos(m,7)*dwwk + QCcell%ATOM_pos(m,8)*dwwk)
      j  = QCcell%ATOM_type(m)
      sf = FSCATT(s,ul,j,smb,sngl(self%voltage),absflg,accflg,dwflg) * QCcell%ATOM_pos(m,6)
      pp = QCcell%ShapeTransformPolygonCa(QCSG, hkl, m)
      sf = sf * cabs(pp)

      ! loop over all atoms in the orbit
      p1 = czero
      do j = 1,QCcell%numat(m)
       arg = twopi*sum(float(hkl(1:5))*QCcell%apos(m,j,1:5))
       p1  = p1 + exp(cmplx(0.0,-arg))
      end do

      ff = ff + p1*real(sf) 
      gg = gg + p1*aimag(sf)  
    end do
end select
!
! these are the modulus and phase of the real part of Vg
Vmod = pref * cabs(ff)
Vphase = atan2(aimag(ff),real(ff))

! modulus of U_g
Umod = preg*Vmod

! if absorption is included, also compute the imaginary part of Vg, i.e., Vprime_g
Vpmod = pref * cabs(gg)
Vpphase = atan2(aimag(gg),real(gg))

! modulus of Uprime_g
Upmod = preg*Vpmod

! complex Ucg = U_g + i Uprime_g = U_g,r-Uprime_g,i + i(U_g,i+Uprime_g,r)
Ucg = pre * cmplx(real(ff)-aimag(gg),aimag(ff)+real(gg))

! complex Vg 
if (abs(Umod).gt.0.0) then 
  xig = 1.0/abs(Umod)/self%Lambda
else
  xig = 1.0E+8
end if 

if (abs(Upmod).gt.0.0) then 
  xgp = 1.0/abs(Upmod)/self%Lambda
else
  xgp = 1.0E+8
end if 

arg = Vpphase - Vphase
qg  = cmplx(cos(Vphase)/xig-sin(Vpphase)/xgp,cos(Vpphase)/xgp+sin(Vphase)/xig)

end function CalcQCUcg_

!--------------------------------------------------------------------------
recursive function CalcsgQC_(self, QCcell, gg, kk, FN) result(sg)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcsgQC_
  !! author: MDG
  !! version: 1.0
  !! date: 02/06/22
  !!
  !! compute the excitation error for a given QC reflection 

use mod_QCcrystallography

IMPLICIT NONE

class(QCDiffraction_T),INTENT(INOUT)  :: self
class(QCCell_T),INTENT(INOUT)         :: QCcell
real(kind=dbl),INTENT(IN)             :: gg(:)
 !! reciprocal lattice point
real(kind=dbl),INTENT(IN)             :: kk(3)
 !! wave vector
real(kind=dbl),INTENT(IN)             :: FN(3)
 !! foil normal

real(kind=dbl)                        :: kpg(3),tkpg(3),xnom,xden,q1,q2,sg
real(kind=dbl)                        :: gvec(3)

! get the g-vector
select type (QCcell)
  class is (QCcell_icosahedral_T)
    gvec = QCcell%getGvector(dble(gg), 'P')
  class is (QCcell_axial_T)
    gvec = QCcell%getGvector(dble(gg), 'P')
end select 

! auxiliary vectors
kpg=kk+gvec
tkpg=2.0*kk+gvec

! use equation of Ewald sphere
xnom = -DOT_PRODUCT(gvec,tkpg)

! 2|k0+g|cos(alpha) = 2(k0+g).Foilnormal
q2 = DOT_PRODUCT(kpg,FN)
xden = 2.D0*q2
sg = xnom/xden

end function CalcsgQC_

!--------------------------------------------------------------------------
recursive function getLUT_(self, QCindex) result(c)
!DEC$ ATTRIBUTES DLLEXPORT :: getLUT_
  !! author: MDG
  !! version: 1.0
  !! date: 02/08/22
  !!
  !! return an entry from the LUT array

IMPLICIT NONE

class(QCDiffraction_T),INTENT(INOUT)  :: self
integer(kind=irg)                     :: QCindex
complex(kind=dbl)                     :: c

 c = self%LUT( QCindex )

end function getLUT_

!--------------------------------------------------------------------------
recursive function getLUTqg_(self, QCindex) result(c)
!DEC$ ATTRIBUTES DLLEXPORT :: getLUTqg_
  !! author: MDG
  !! version: 1.0
  !! date: 02/08/22
  !!
  !! return an entry from the LUTqg array

IMPLICIT NONE

class(QCDiffraction_T),INTENT(INOUT)  :: self
integer(kind=irg)                     :: QCindex
complex(kind=dbl)                     :: c

c = self%LUTqg( QCindex )

end function getLUTqg_

!--------------------------------------------------------------------------
recursive function getBetheParameter_(self, c) result(bp)
!DEC$ ATTRIBUTES DLLEXPORT :: getBetheParameter_
  !! author: MDG
  !! version: 1.0
  !! date: 02/02/20
  !!
  !! returns one of the Bethe parameters

use mod_io

IMPLICIT NONE

class(QCDiffraction_T),INTENT(INOUT)  :: self
character(2),INTENT(IN)               :: c
real(kind=dbl)                        :: bp

type(IO_T)                            :: Message

select case(c)
  case('c1')
    bp = self%BetheParameters%c1
  case('c2')
    bp = self%BetheParameters%c2
  case('c3')
    bp = self%BetheParameters%c3
  case('sg')
    bp = self%BetheParameters%sgdbdiff
  case default
    call Message%printWarning('getBetheParameter: unknown Bethe parameter '//c//'; returning 0.0')
    bp = 0.D0
end select

end function getBetheParameter_

!--------------------------------------------------------------------------
recursive function getWaveLength_(self) result(Lambda)
!DEC$ ATTRIBUTES DLLEXPORT :: getWaveLength_
  !! author: MDG
  !! version: 1.0
  !! date: 02/08/22
  !!
  !! returns the electron wavelength

IMPLICIT NONE

class(QCDiffraction_T),INTENT(INOUT)  :: self
real(kind=dbl)                        :: Lambda

Lambda = self%Lambda

end function getWaveLength_

!--------------------------------------------------------------------------
recursive function getshapeLUT_(self) result(s)
!DEC$ ATTRIBUTES DLLEXPORT :: getshapeLUT_
  !! author: MDG
  !! version: 1.0
  !! date: 02/08/22
  !!
  !! returns the shape of the LUT array

IMPLICIT NONE

class(QCDiffraction_T),INTENT(INOUT)  :: self
integer(kind=irg)                     :: s(1)

s = shape(self%LUT)

end function getshapeLUT_

!--------------------------------------------------------------------------
recursive function getV_(self) result(V)
!DEC$ ATTRIBUTES DLLEXPORT :: getV_
  !! author: MDG
  !! version: 1.0
  !! date: 01/28/20
  !!
  !! return the accelerating voltage

IMPLICIT NONE

class(QCDiffraction_T),INTENT(INOUT)  :: self
real(kind=dbl)                        :: V

V = self%voltage

end function getV_

!--------------------------------------------------------------------------
recursive subroutine setV_(self, V)
!DEC$ ATTRIBUTES DLLEXPORT :: setV_
  !! author: MDG
  !! version: 1.0
  !! date: 02/05/20
  !!
  !! set the accelerating voltage

IMPLICIT NONE

class(QCDiffraction_T),INTENT(INOUT)  :: self
real(kind=dbl),INTENT(IN)             :: V

self%voltage = V

end subroutine setV_

!--------------------------------------------------------------------------
recursive subroutine allocateLUT_(self, nLUT)
!DEC$ ATTRIBUTES DLLEXPORT :: allocateLUT_
  !! author: MDG
  !! version: 1.0
  !! date: 02/09/22 
  !!
  !! allocates all the LUT arrays

use mod_memory

IMPLICIT NONE

class(QCDiffraction_T),INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)          :: nLUT

type(memory_T)                        :: mem

mem = memory_T() 

! the LUT arrays store all the Fourier coefficients etc...
call mem%alloc(self%LUT, (/ nLUT /), 'self%LUT', initval = cmplx(0.D0,0.D0))
call mem%alloc(self%LUTqg, (/ nLUT /), 'self%LUTqg', initval = cmplx(0.D0,0.D0))
call mem%alloc(self%dbdiff, (/ nLUT /), 'self%dbdiff', initval = .FALSE.)
! call mem%alloc(self%inverseIndex, (/ nLUT, m /), 'self%inverseIndex', initval = 0)

end subroutine allocateLUT_

!--------------------------------------------------------------------------
recursive subroutine setLUT_(self, id, c)
!DEC$ ATTRIBUTES DLLEXPORT :: setLUT_
  !! author: MDG
  !! version: 1.0
  !! date: 02/04/20
  !!
  !! set an entry in the LUT array

IMPLICIT NONE

class(QCDiffraction_T),INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)          :: id
complex(kind=dbl),INTENT(IN)          :: c

self%LUT( id ) = c

end subroutine setLUT_

!--------------------------------------------------------------------------
recursive subroutine setLUTqg_(self, id, c)
!DEC$ ATTRIBUTES DLLEXPORT :: setLUTqg_
  !! author: MDG
  !! version: 1.0
  !! date: 02/04/20
  !!
  !! set an entry in the LUTqg array

IMPLICIT NONE

class(QCDiffraction_T),INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)          :: id
complex(kind=dbl),INTENT(IN)          :: c

self%LUTqg( id ) = c

end subroutine setLUTqg_

!--------------------------------------------------------------------------
recursive subroutine setdbdiff_(self, id, c)
!DEC$ ATTRIBUTES DLLEXPORT :: setdbdiff_
  !! author: MDG
  !! version: 1.0
  !! date: 02/04/20
  !!
  !! set an entry in the dbdiff array

IMPLICIT NONE

class(QCDiffraction_T),INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)          :: id
logical,INTENT(IN)                    :: c

self%dbdiff( id ) = c

end subroutine setdbdiff_

  ! !--------------------------------------------------------------------------
  ! recursive function getrlp_(self) result(rlp)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: getrlp_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 01/26/20
  !   !!
  !   !! ask for accelerating voltage, then call CalcWaveLength
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T),INTENT(INOUT)  :: self
  ! type(gnode)                         :: rlp
  
  ! rlp = self%rlp
  
  ! end function getrlp_
  
  ! !--------------------------------------------------------------------------
  ! recursive subroutine allocateLUT_(self, dims)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: allocateLUT_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 02/04/20, incorporated memory_T class on 04/15/21
  !   !!
  !   !! allocates the LUT arrays
  
  ! use mod_memory
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T),INTENT(INOUT)  :: self
  ! integer(kind=irg)                   :: dims(3)
  
  ! integer(kind=irg)                   :: imh, imk, iml
  ! type(memory_T)                      :: mem
  
  ! imh = 2*dims(1)
  ! imk = 2*dims(2)
  ! iml = 2*dims(3)
  
  ! mem = memory_T() 

  ! ! the LUT array stores all the Fourier coefficients, so that we only need to compute them once... i.e., here and now
  ! call mem%alloc(self%LUT, (/ imh, imk, iml /), 'self%LUT', cmplx(0.D0,0.D0), startdims=(/ -imh, -imk, -iml /))
  ! call mem%alloc(self%LUTqg, (/ imh, imk, iml /), 'self%LUTqg', cmplx(0.D0,0.D0), startdims=(/ -imh, -imk, -iml /))
  ! call mem%alloc(self%dbdiff, (/ imh, imk, iml /), 'self%dbdiff', .FALSE., startdims=(/ -imh, -imk, -iml /))
  
  ! end subroutine allocateLUT_
  

  
 
  
  ! !--------------------------------------------------------------------------
  ! recursive function getSghLUT_(self, n, hkl) result(c)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: getSghLUT_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 02/12/20
  !   !!
  !   !! return an entry from the LUT array
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T),INTENT(INOUT)  :: self
  ! integer(kind=irg), INTENT(IN)       :: n
  ! integer(kind=irg), INTENT(IN)       :: hkl(3)
  ! complex(kind=dbl)                   :: c(n)
  
  !  c(1:n) = self%SghLUT( 1:n, hkl(1), hkl(2), hkl(3) )
  
  ! end function getSghLUT_
  

  
  ! !--------------------------------------------------------------------------
  ! recursive function getdbdiff_(self, hkl) result(c)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: getdbdiff_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 02/04/20
  !   !!
  !   !! return an entry from the dbdiff array
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T),INTENT(INOUT)  :: self
  ! integer(kind=irg)                   :: hkl(3)
  ! logical                             :: c
  
  !  c = self%dbdiff( hkl(1), hkl(2), hkl(3) )
  
  ! end function getdbdiff_
  

  
  ! !--------------------------------------------------------------------------
  ! recursive function getrlpmethod_(self) result(rlpmethod)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: getrlpmethod_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 01/26/20
  !   !!
  !   !! retrieve the rlp method
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T),INTENT(INOUT)  :: self
  ! character(2)                        :: rlpmethod
  
  ! rlpmethod = self%rlp%method
  
  ! end function getrlpmethod_
  

  
  ! !--------------------------------------------------------------------------
  ! recursive function getV0mod_(self) result(V)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: getV0mod_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 02/14/20
  !   !!
  !   !! return V0mod
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T),INTENT(INOUT)  :: self
  ! real(kind=dbl)                      :: V
  
  ! V = self%V0mod
  
  ! end function getV0mod_
  
  ! !--------------------------------------------------------------------------
  ! recursive function getRelcor_(self) result(V)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: getRelcor_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 01/14/20
  !   !!
  !   !! return the relativistic correction factor
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T),INTENT(INOUT)  :: self
  ! real(kind=dbl)                      :: V
  
  ! V = self%Relcor
  
  ! end function getRelcor_
  
  ! !--------------------------------------------------------------------------
  ! recursive function getSigma_(self) result(V)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: getSigma_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 01/14/20
  !   !!
  !   !! return the interaction constant
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T),INTENT(INOUT)  :: self
  ! real(kind=dbl)                      :: V
  
  ! V = self%Sigma
  
  ! end function getSigma_
  
  ! !--------------------------------------------------------------------------
  ! recursive function getPsihat_(self) result(V)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: getPsihat_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 01/14/20
  !   !!
  !   !! return the relativistically corrected accelerating voltage
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T),INTENT(INOUT)  :: self
  ! real(kind=dbl)                      :: V
  
  ! V = self%Psihat
  
  ! end function getPsihat_
  
  ! !--------------------------------------------------------------------------
  ! recursive subroutine Printrlp_(self,first)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: Printrlp_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 01/09/20
  !   !!
  !   !! output the contents of the rlp structure
  
  ! use mod_io
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T),INTENT(INOUT)      :: self
  ! logical,optional,intent(INOUT)          :: first
  !  !! switch for long/short output
  
  ! type(IO_T)                              :: Message
  ! integer(kind=irg)                       :: oi_int(3)
  ! real(kind=sgl)                          :: oi_real(7)
  ! complex(kind=sgl)                       :: oi_cmplx(1)
  
  
  ! if (present(first)) then
  !  if (first) then
  !   call Message%printMessage('     Scattering factors : ', frm = "(/A)",advance="no")
  !   if (self%rlp%method.eq.'WK') then
  !    if (self%rlp%absorption.eqv..TRUE.) then
  !     call Message%printMessage(' Weickenmeier-Kohl (with absorption)', frm = "(A/)")
  !    else
  !     call Message%printMessage(' Weickenmeier-Kohl', frm = "(A/)")
  !    end if
  !   else
  !     call Message%printMessage(' Doyle-Turner/Smith-Burge', frm = "(A/)")
  !   end if
  
  !   if (self%rlp%absorption.eqv..TRUE.) then
  !     call Message%printMessage( &
  !         '   h  k  l    |g|    Ucg_r  Ucg_i   |Ug|    phase   |Ugp|   phase   xi_g   xi_gp    ratio  Re-1/q_g-Im', &
  !         frm = "(A)")
  !   else
  !     call Message%printMessage('   h  k  l    |g|    Ucg_r  |Ug|    phase    xi_g   1/q_g', frm = "(A)")
  !   end if
  !   first = .FALSE.
  !  end if
  ! end if
  
  ! if (self%rlp%absorption.eqv..TRUE.) then
  !  oi_int(1:3) = self%rlp%hkl(1:3)
  !  call Message%WriteValue('',oi_int, 3, "(1x,3I3,1x)",advance="no")
  !  oi_real(1) = self%rlp%g
  !  call Message%WriteValue('',oi_real, 1, "(F9.4)",advance="no")
  !  oi_cmplx(1) = self%rlp%Ucg
  !  call Message%WriteValue('',oi_cmplx, 1, "(2F7.3,1x)",advance="no")
  !  oi_real(1:7)  = (/ self%rlp%Umod,self%rlp%Vphase*180.0/sngl(cPi),self%rlp%Upmod,self%rlp%Vpphase*180.0/sngl(cPi),&
  !                     self%rlp%xg,self%rlp%xgp,self%rlp%ar /)
  !  call Message%WriteValue('',oi_real, 7, "(4F8.3,3F8.1)",advance="no")
  !  oi_cmplx(1) = self%rlp%qg
  !  call Message%WriteValue('',oi_cmplx, 1, "(2F8.3)")
  ! else
  !  oi_int(1:3) = self%rlp%hkl(1:3)
  !  call Message%WriteValue('',oi_int, 3, "(1x,3I3,1x)",advance="no")
  !  oi_real(1) = self%rlp%g
  !  call Message%WriteValue('',oi_real, 1, "(F9.4)",advance="no")
  !  oi_real(1) = real(self%rlp%Ucg)
  !  call Message%WriteValue('',oi_real, 1, "(F7.3,1x)",advance="no")
  !  oi_real(1:3)  = (/ self%rlp%Umod,self%rlp%Vphase*180.0/sngl(cPi),self%rlp%xg /)
  !  call Message%WriteValue('',oi_real, 3, "(2F8.3,F8.1)",advance="no")
  !  oi_cmplx(1) = self%rlp%qg
  !  call Message%WriteValue('',oi_cmplx, 1, "(2F8.3)")
  ! end if
  
  ! end subroutine Printrlp_
  
  ! !--------------------------------------------------------------------------
  ! recursive subroutine GetVoltage_(self, cell, verbose)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: GetVoltage_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 01/26/20
  !   !!
  !   !! ask for accelerating voltage, then call CalcWaveLength
  
  ! use mod_io
  ! use mod_crystallography
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T),INTENT(INOUT)  :: self
  ! type(Cell_T),INTENT(INOUT)          :: cell
  ! logical,INTENT(IN),OPTIONAL         :: verbose
  
  
  ! real(kind=sgl)                      :: io_real(1)
  ! type(IO_T)                          :: Message
  
  ! call Message%ReadValue(' Enter the microscope accelerating voltage [kV, R] : ', io_real, 1)
  ! self%voltage = dble(io_real(1))
  
  ! if (present(verbose)) then
  !     if (verbose.eqv..TRUE.) call self%CalcWaveLength(cell, verbose = .TRUE.)
  ! else
  !     call self%CalcWaveLength(cell)
  ! end if
  
  ! end subroutine GetVoltage_
  
  ! !--------------------------------------------------------------------------
  ! recursive function CalcDiffAngle_(self, cell, hkl) result(tt)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: CalcDiffAngle_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 01/26/20
  !   !!
  !   !! compute the diffraction angle 2theta in radians
  
  ! use mod_crystallography
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T),INTENT(INOUT)  :: self
  ! type(Cell_T),INTENT(IN)             :: cell
  ! integer(kind=irg),INTENT(IN)        :: hkl(3)
  !  !! Miller indices
  
  ! real(kind=sgl)                      :: tt
  
  ! tt = 2.0*asin(0.50*sngl(self%Lambda)*cell%CalcLength(float(hkl), 'r') )
  
  ! end function CalcDiffAngle_
  

  
  ! !--------------------------------------------------------------------------
  ! recursive subroutine PreCalcFSCATT_(self, cell, dmin, gstep)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: PreCalcFSCATT_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 01/27/20
  !   !!
  !   !! precompute the FSCATT values for interpolation purposes, to speed up STEM-DCI and other codes
  
  ! use mod_crystallography
  ! use mod_others
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T),INTENT(INOUT)  :: self
  ! type(Cell_T),INTENT(INOUT)          :: cell
  ! real(kind=sgl),INTENT(IN)           :: dmin
  !  !! smallest d-spacing to consider
  ! real(kind=sgl),INTENT(IN)           :: gstep
  !  !! step size in the cell%scatfacg array
  
  ! integer(kind=irg)                   :: j,m,ii,i
  ! real(kind=sgl)                      :: s,ul
  ! real(kind=sgl),parameter            :: swk = 0.628318530717959
  ! real(kind=sgl),parameter            :: dwwk = 1.26651479552922
  ! integer(kind=irg),parameter         :: absflg = 3
  ! logical                             :: accflg=.TRUE., dwflg=.TRUE.
  ! character(2)                        :: smb
  ! real(kind=dbl), allocatable         :: apos(:,:)
  ! integer(kind=irg), allocatable      :: atp(:)
  
  ! ! first generate the array of s-values for which the scattering factors need to be computed
  ! s = 2.0/dmin   ! maximum range in reciprocal space
  ! self%numscatfac = nint(s/gstep) + 2
  ! allocate(self%scatfacg(self%numscatfac))
  ! self%scatfacg = (/ (gstep * float(i-1),i=1,self%numscatfac) /)
  ! self%scatfacg = self%scatfacg * swk
  
  ! ! allocate the scattering factor interpolation array
  ! allocate( self%scatfac(self%numscatfac, cell%getNatomtype()) )
  
  ! ! The Weickenmeier-Kohl (WK) subroutine works in Angstrom, and also
  ! ! scales reciprocal space by a factor of 2*pi;  this scaling
  ! ! is accomplished by changing the g-value in nm^{-1} by a
  ! ! scaling factor swk = 2*pi/10, to go from book units to WK units.
  ! !
  ! ! A similar scaling must be performed on the Debye Waller factor;
  ! ! the book defines it as exp(-Bs^2), with s in [nm^2]; WK define
  ! ! it in A^2, and with a scaled reciprocal space.  The conversion
  ! ! factor dwwk = 100.0/8*pi^2
  ! !
  ! ! To go from the standard B factor in [nm^2] to ul^2 in A^2,
  ! ! ul = sqrt(B*dwwk)
  
  ! apos = cell%getAsymPosData()
  ! atp = cell%getatomtype()
  
  ! do i=1,self%numscatfac
  ! ! properly scale the scattering parameter
  !  s = self%scatfacg(i)
  
  ! ! loop over all atoms in the asymmetric unit
  !  do m=1,cell%getNatomtype()
  ! ! get the atomic scattering factor for this atom
  ! ! scale and include Debye-Waller factor and site occupation parameter
  !   ul = sqrt(apos(m,5)*dwwk)
  !   j = atp(m)
  !   self%scatfac(i,m) = FSCATT(s,ul,j,smb,sngl(self%voltage),absflg,accflg,dwflg)*apos(m,4)
  !  end do
  ! end do
  
  ! end subroutine PreCalcFSCATT_
  
  ! !--------------------------------------------------------------------------
  ! recursive subroutine getScatfac_(self, cell, s, sfarray, ntypes)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: getScatfac_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 01/27/20
  !   !!
  !   !! interpolate the precomputed FSCATT values
  
  ! use mod_crystallography
  ! use mod_others
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T),INTENT(INOUT)  :: self
  ! type(Cell_T),INTENT(INOUT)          :: cell
  ! real(kind=sgl),INTENT(IN)           :: s
  !  !! reciprocal distance value
  ! integer(kind=irg),INTENT(IN)        :: ntypes
  !  !! number of types
  ! complex(kind=sgl),INTENT(OUT)       :: sfarray(ntypes)
  !  !! returned scattering factor values
  
  ! integer(kind=irg)                   :: jj
  ! real(kind=sgl)                      :: dx
  
  ! if (s.eq.0.0) then
  !     sfarray(1:ntypes) = self%scatfac(1,1:ntypes)
  ! else
  !     jj = ifix(s/self%scatfacg(2))
  !     if (jj.ge.self%numscatfac) then
  !         sfarray(1:ntypes) = self%scatfac(self%numscatfac,1:ntypes)
  !     else
  !         dx = s/self%scatfacg(2) - float(jj)
  !         sfarray(1:ntypes) = self%scatfac(jj,1:ntypes)*(1.0-dx) + &
  !                                      self%scatfac(jj+1,1:ntypes)*dx
  !     end if
  ! end if
  
  ! end subroutine getScatfac_
  
 
  

  
  ! !--------------------------------------------------------------------------
  ! recursive subroutine Set_Bethe_Parameters_(self, EMsoft, silent, filename, usevalues)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: Set_Bethe_Parameters_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 02/02/20; added usevalues on 04/15/21
  !   !!
  !   !! Read the Bethe potential parameters from a file, if it exists; otherwise take defaults
  !   !!
  !   !! The parameters set in this routine determine the difference between strong and
  !   !! weak reflections.  The routine checks for the presence of the BetheParameters.nml file
  !   !! in the current folder.  If present, it will read the parameters, otherwise it will use
  !   !! defaults which have been determined to be reasonable based on dynamical EBSD runs.
  
  ! use mod_EMsoft
  ! use mod_io
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T), INTENT(INOUT)         :: self
  ! type(EMsoft_T), INTENT(INOUT),OPTIONAL      :: EMsoft
  ! logical,INTENT(IN),OPTIONAL                 :: silent
  ! character(fnlen),INTENT(IN),OPTIONAL        :: filename
  ! real(kind=sgl),INTENT(IN),OPTIONAL          :: usevalues(4)
  
  ! type(IO_T)                                  :: Message
  ! character(fnlen)                            :: Bethefilename, fname
  ! logical                                     :: fexist
  ! real(kind=sgl)                              :: c1, c2, c3, sgdbdiff
  
  ! namelist /BetheList/ c1, c2, c3, sgdbdiff
  
  ! if (present(usevalues)) then 
  !   self%BetheParameters%c1 = usevalues(1)
  !   self%BetheParameters%c2 = usevalues(2)
  !   self%BetheParameters%c3 = usevalues(3)
  !   self%BetheParameters%sgdbdiff = usevalues(4)
  ! else

  !   if (present(filename)) then
  !     Bethefilename = trim(filename)
  !   else
  !     Bethefilename = 'BetheParameters.nml'
  !   end if
    
  !   ! check for the presence of the namelist file in the current folder
  !   inquire(file=trim(Bethefilename),exist=fexist)
    
  !   ! set all default values (must be done here, since nml file may not contain all of them)
  !   c1 = 20.0_sgl           ! changed from 8 and 12 for a test on 8/14/15
  !   c2 = 40.0_sgl           !
  !   c3 = 200.0_sgl          !
  !   sgdbdiff = 1.00_sgl    !
    
  !   if (fexist) then ! check for the file in the local folder
  !   ! read the parameters from the file
  !    fname = EMsoft%toNativePath(Bethefilename)
  !    open(UNIT=dataunit,FILE=trim(fname),DELIM='APOSTROPHE')
  !    READ(UNIT=dataunit,NML=BetheList)
  !    close(UNIT=dataunit)
  !    if (.not.present(silent)) then
  !      call Message%printMessage('Read Bethe parameters from BetheParameters.nml', frm = "(A)")
  !      write (6,nml=BetheList)
  !    end if
  !   end if
    
  !   self%BetheParameters%c1 = c1
  !   self%BetheParameters%c2 = c2
  !   self%BetheParameters%c3 = c3
  !   self%BetheParameters%sgdbdiff = sgdbdiff
  ! end if 

  ! end subroutine Set_Bethe_Parameters_
  
  ! !--------------------------------------------------------------------------
  ! recursive subroutine writeBetheparameterNameList_(self, HDF)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: writeBetheparameterNameList_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 02/12/20
  !   !!
  !   !! write namelist to already opened HDF file
  
  ! use HDF5
  ! use mod_HDFsupport
  ! use ISO_C_BINDING
  ! use stringconstants
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T), INTENT(INOUT)                   :: self
  ! type(HDF_T), INTENT(INOUT)                            :: HDF
  
  ! integer(kind=irg),parameter                           :: n_int = 1, n_real = 4
  ! integer(kind=irg)                                     :: hdferr,  io_int(n_int), restart, uniform
  ! real(kind=sgl)                                        :: io_real(n_real)
  ! character(20)                                         :: intlist(n_int), reallist(n_real)
  ! character(fnlen)                                      :: dataset, groupname
  ! character(fnlen,kind=c_char)                          :: line2(1)
  ! logical                                               :: g_exists, overwrite=.TRUE.
  
  ! ! We are not writing the complete BetheParameters structure to the name list, only
  ! ! the parameters of importance for the dynamical simulations ...
  
  ! ! create the group for this namelist
  ! groupname = SC_BetheList
  ! hdferr = HDF%createGroup(groupname)
  
  ! ! write all the single doubles
  ! io_real = (/ self%BetheParameters%c1, self%BetheParameters%c2, self%Betheparameters%c3, self%BetheParameters%sgdbdiff /)
  
  ! reallist(1) = 'c1'
  ! reallist(2) = 'c2'
  ! reallist(3) = 'c3'
  ! reallist(4) = 'sgdbdiff'
  ! call HDF%writeNMLreals(io_real, reallist, n_real)
  
  ! ! and pop this group off the stack
  ! call HDF%pop()
  
  ! end subroutine writeBetheparameterNameList_
  
  
  ! !--------------------------------------------------------------------------
  ! recursive function CalcsgSingle_(self, cell, gg, kk, FN) result(sg)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: CalcsgSingle_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 01/27/20
  !   !!
  !   !! compute the excitation error for a given reflection (single precision)
  
  ! use mod_crystallography
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T),INTENT(INOUT)  :: self
  ! type(Cell_T),INTENT(INOUT)          :: cell
  ! real(kind=sgl),INTENT(IN)           :: gg(3)
  !  !! reciprocal lattice point
  ! real(kind=sgl),INTENT(IN)           :: kk(3)
  !  !! wave vector
  ! real(kind=sgl),INTENT(IN)           :: FN(3)
  !  !! foil normal
  
  ! real(kind=sgl)                      :: kpg(3),tkpg(3),xnom,xden,q1,q2,sg
  
  !  kpg=kk+gg
  !  tkpg=2.0*kk+gg
  
  ! ! use equation of Ewald sphere
  !  xnom = -cell%CalcDot(gg,tkpg,'r')
  
  ! ! 2|k0+g|cos(alpha) = 2(k0+g).Foilnormal
  !  q1 = cell%CalcLength(kpg,'r')
  !  q2 = cell%CalcAngle(kpg,FN,'r')
  !  xden = 2.0*q1*cos(q2)
  !  sg = xnom/xden
  
  ! end function CalcsgSingle_
  
 
  ! !--------------------------------------------------------------------------
  ! recursive subroutine BWsolve_(self, M, W, CGG, CGinv, nn, IPIV)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: BWsolve_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 01/27/20
  !   !!
  !   !! Solve the Bloch wave eigenvalue equation for the
  !   !! N-beam case, using the ZGEEV LAPACK 3.0 routine
  
  ! use mod_io
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T),INTENT(INOUT)  :: self
  ! integer(kind=irg),INTENT(IN)        :: nn
  !  ! number of beams
  ! complex(kind=dbl),INTENT(IN)        :: M(nn,nn)
  !  ! input dynamical matrix
  ! complex(kind=dbl),INTENT(OUT)       :: W(nn)
  !  ! Bloch eigenvalues
  ! complex(kind=dbl),INTENT(OUT)       :: CGG(nn,nn)
  !  ! Bloch eigenvectors
  ! complex(kind=dbl),INTENT(OUT)       :: CGinv(nn,nn)
  !  ! inverse of eigenvector array
  ! integer(kind=irg),INTENT(IN)        :: IPIV(nn)
  !  ! pivot array, currently unused
  
  ! type(IO_T)                          :: Message
  ! integer(kind=irg)                   :: INFO, LDA, LDVR, LDVL, LWORK, JPIV(nn),MILWORK, i, io_int(1)
  ! integer(kind=irg),parameter         :: LWMAX = 5000
  ! complex(kind=dbl)                   :: VL(nn,nn),  WORK(LWMAX), normsum
  ! real(kind=dbl)                      :: RWORK(2*nn), io_real(1)
  ! character                           :: JOBVL, JOBVR
  ! complex(kind=dbl),allocatable       :: MIWORK(:)
  
  ! ! set some initial LAPACK variables
  !  LDA = nn
  !  LDVL = nn
  !  LDVR = nn
  !  INFO = 0
  
  ! ! first initialize the parameters for the LAPACK ZGEEV, CGETRF, and CGETRI routines
  !  JOBVL = 'N'   ! do not compute the left eigenvectors
  !  JOBVR = 'V'   ! do compute the right eigenvectors
  !  LWORK = -1    ! so that we can ask the routine for the actually needed value
  
  ! ! call the routine to determine the optimal workspace size
  !   call zgeev(JOBVL,JOBVR,nn,M,LDA,W,VL,LDVL,CGG,LDVR,WORK,LWORK,RWORK,INFO)
  !   LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
  
  ! ! then call the eigenvalue solver
  !   call zgeev(JOBVL,JOBVR,nn,M,LDA,W,VL,LDVL,CGG,LDVR,WORK,LWORK,RWORK,INFO)
  !   if (INFO.ne.0) call Message%printError('Error in BWsolve: ','ZGEEV return not zero')
  
  ! ! it appears that the eigenvectors may not always be normalized ...
  ! ! so we renormalize them here...
  ! ! do i=1,nn
  ! !   normsum = sum(abs(CGG(1:nn,i))**2)
  ! !   normsum = cmplx(1.0,0.0,dbl)/sqrt(normsum)
  ! !   CGG(1:nn,i) = CGG(1:nn,i)*normsum
  ! ! end do
  
  ! ! make a copy of CG for the matrix inversion routines
  !  CGinv = CGG
  
  ! ! invert CGinv to get the Bloch wave excitation amplitudes
  !  LDA = nn
  !  call zgetrf(nn,nn,CGinv,LDA,JPIV,INFO)
  !  if (INFO.ne.0) then
  !   io_int(1) = INFO
  !   call Message%WriteValue('zgetrf error code: ',io_int,1,frm="(I5)")
  !   call Message%printError('Error in BWsolve: ','ZGETRF return not zero')
  !  end if
  
  !  MILWORK = 64*nn
  !  allocate(MIWORK(MILWORK))
  
  !  MIWORK = cmplx(0.0_dbl,0.0_dbl)
  !  call zgetri(nn,CGinv,LDA,JPIV,MIWORK,MILWORK,INFO)
  !  if (INFO.ne.0) then
  !   io_int(1) = INFO
  !   call Message%WriteValue('zgetri error code: ',io_int,1,frm="(I5)")
  !   call Message%printError('Error in BWsolve: ','ZGETRI return not zero')
  !  end if
  
  ! ! if ((abs(sum(matmul(CGG,CGinv)))-dble(nn)).gt.1.E-8) then
  ! !  call Message('Error in matrix inversion; continuing', frm = "(A)")
  ! !  io_real(1) = abs(sum(matmul(CGG,CGinv)))-dble(nn)
  ! !  call WriteValue('   Matrix inversion error; this number should be zero: ',io_real,1,"(F)")
  ! ! endif
  
  !  deallocate(MIWORK)
  
  ! end subroutine BWsolve_
  
  ! !--------------------------------------------------------------------------
  ! recursive subroutine Initialize_SghLUT_(self, cell, SG, dmin, numset, nat, verbose)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: Initialize_SghLUT_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 02/11/20
  !   !!
  !   !! initialize a look-up table with Sgh structure-factor-like entries
  !   !!
  
  ! use mod_crystallography
  ! use mod_symmetry
  ! use mod_io
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T), INTENT(INOUT)        :: self
  ! type(Cell_T), INTENT(INOUT)                :: cell
  ! type(SpaceGroup_T), INTENT(INOUT)          :: SG
  ! real(kind=sgl),INTENT(IN)                  :: dmin
  ! integer(kind=sgl),INTENT(IN)               :: numset
  ! integer(kind=sgl),INTENT(INOUT)            :: nat(maxpasym)
  ! logical,INTENT(IN),optional                :: verbose
  
  ! type(IO_T)                                 :: Message
  ! integer(kind=irg)                          :: istat, io_int(3), skip
  ! integer(kind=irg)                          :: imh, imk, iml, gg(3), ix, iy, iz
  ! real(kind=sgl)                             :: dhkl, io_real(3), ddt
  ! complex(kind=dbl)                          :: Sghvec(numset)
  
  ! ! compute the range of reflections for the lookup table and allocate the table
  ! ! The master list is easily created by brute force
  !  imh = 1
  !  do
  !    dhkl = 1.0/cell%CalcLength( (/float(imh) ,0.0_sgl,0.0_sgl/), 'r')
  !    if (dhkl.lt.dmin) EXIT
  !    imh = imh + 1
  !  end do
  !  imk = 1
  !  do
  !    dhkl = 1.0/cell%CalcLength( (/0.0_sgl,float(imk),0.0_sgl/), 'r')
  !    if (dhkl.lt.dmin) EXIT
  !    imk = imk + 1
  !  end do
  !  iml = 1
  !  do
  !    dhkl = 1.0/cell%CalcLength( (/0.0_sgl,0.0_sgl,float(iml)/), 'r')
  !    if (dhkl.lt.dmin) EXIT
  !    iml = iml + 1
  !  end do
  
  !  if (present(verbose)) then
  !   if (verbose) then
  !     io_int = (/ imh, imk, iml /)
  !     call Message%WriteValue(' Range of reflections along a*, b* and c* = ',io_int,3)
  !   end if
  !  end if
  
  ! ! the SghLUT array stores all the structure-factor-like quantities of the Sgh matrix needed for EBSD, ECP, etc simulations
  !  allocate(self%SghLUT(1:numset,-2*imh:2*imh,-2*imk:2*imk,-2*iml:2*iml),stat=istat)
  !  if (istat.ne.0) call Message%printError('Initialize_SghLUT:',' unable to allocate SghLUT array')
  !  self%SghLUT = cmplx(0.D0,0.D0)
  
  !  if (present(verbose)) then
  !   if (verbose) then
  !    call Message%printMessage(' Generating Sgh coefficient lookup table ... ', frm = "(/A)",advance="no")
  !   end if
  !  end if
  
  ! ! note that the lookup table must be twice as large as the list of participating reflections,
  ! ! since the Sgh matrix uses g-h as its index !!!
  ! izl: do iz=-2*iml,2*iml
  ! iyl:  do iy=-2*imk,2*imk
  ! ixl:   do ix=-2*imh,2*imh
  !         gg = (/ ix, iy, iz /)
  !         if (SG%IsGAllowed(gg)) then  ! is this reflection allowed by lattice centering ?
  ! ! add the reflection to the look up table
  !            call self%preCalcSgh(cell,SG,gg,numset,nat,Sghvec)
  !            self%SghLUT(1:numset, ix, iy, iz) = Sghvec(1:numset)
  !         end if ! IsGAllowed
  !        end do ixl
  !       end do iyl
  !     end do izl
  
  !   if (present(verbose)) then
  !    if (verbose) then
  !     call Message%printMessage('Done', frm = "(A/)")
  !    end if
  !   end if
  
  ! ! that's it
  ! end subroutine Initialize_SghLUT_
  
  ! !--------------------------------------------------------------------------
  ! recursive subroutine preCalcSgh_(self, cell, SG, kkk, numset, nat, Sghvec)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: preCalcSgh_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 02/11/20
  !   !!
  !   !! compute structure factor-like Sgh array entry for EBSD, ECCI and ECP simulations
  
  ! use mod_crystallography
  ! use mod_symmetry
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T), INTENT(INOUT)     :: self
  ! type(Cell_T), INTENT(INOUT)             :: cell
  ! type(SpaceGroup_T),INTENT(INOUT)        :: SG
  ! integer(kind=irg),INTENT(IN)            :: kkk(3)
  ! integer(kind=irg),INTENT(IN)            :: numset
  ! integer(kind=irg),INTENT(INOUT)         :: nat(maxpasym)
  ! complex(kind=dbl),INTENT(INOUT)         :: Sghvec(numset)
  
  ! integer(kind=irg)                       :: ip, ir, ic, ikk, n
  ! real(kind=sgl)                          :: Znsq, DBWF, kkl
  ! complex(kind=dbl)                       :: carg
  ! real(kind=dbl),allocatable              :: ctmp(:,:), apdata(:,:), apos(:,:,:)
  ! real(kind=dbl)                          :: arg, tpi
  ! integer(kind=irg), allocatable          :: numat(:), atp(:)
  
  !   tpi = 2.D0 * cPi
  !   Sghvec = cmplx(0.D0,0.D0)
  
  !   apdata = cell%getAsymPosData()
  !   numat = cell%getnumat()
  !   atp = cell%getatomtype()
  !   apos = cell%getapos()
  
  ! ! for each special position we need to compute its contribution to the Sgh array
  !   do ip=1,cell%getNatomtype()
  !     call SG%CalcOrbit(apdata(ip,1:3),n,ctmp)
  !     nat(ip) = numat(ip)
  ! ! get Zn-squared for this special position, and include the site occupation parameter as well
  !     Znsq = float(atp(ip))**2 * apdata(ip,4)
  ! ! We'll assume isotropic Debye-Waller factors for now ...
  ! ! That means we need the square of the length of s=kk^2/4
  !     kkl = 0.25 * cell%CalcLength(float(kkk),'r')**2
  ! ! Debye-Waller exponential times Z^2
  !     DBWF = Znsq * exp(-apdata(ip,5)*kkl)
  ! ! here is where we insert the proper weight factor, Z^2 exp[-M_{h-g}]
  ! ! and also the detector geometry...   For now, we do nothing with the detector
  ! ! geometry; the Rossouw et al 1994 paper lists a factor A that does not depend
  ! ! on anything in particular, so we assume it is 1.
  !     do ikk=1,n ! sum over all the atoms in this orbit
  ! ! get the argument of the complex exponential
  !       arg = tpi*sum(dble(kkk(1:3))*apos(ip,ikk,1:3))
  !       carg = cmplx(dcos(arg),dsin(arg))
  ! ! multiply with the prefactor and add
  !       Sghvec(ip) = Sghvec(ip) + carg * cmplx(DBWF,0.D0)
  !     end do
  !   end do
  
  ! end subroutine preCalcSgh_
  
  ! !--------------------------------------------------------------------------
  ! recursive subroutine preCalcSghECCI_(self, cell, SG, numset, nat, Sghvec)
  ! !DEC$ ATTRIBUTES DLLEXPORT :: preCalcSghECCI_
  !   !! author: MDG
  !   !! version: 1.0
  !   !! date: 02/11/20
  !   !!
  !   !! compute structure factor-like Sgh array entry for EBSD, ECCI and ECP simulations
  
  ! use mod_crystallography
  ! use mod_symmetry
  
  ! IMPLICIT NONE
  
  ! class(Diffraction_T), INTENT(INOUT)     :: self
  ! type(Cell_T), INTENT(INOUT)             :: cell
  ! type(SpaceGroup_T),INTENT(INOUT)        :: SG
  ! integer(kind=irg),INTENT(IN)            :: numset
  ! integer(kind=irg),INTENT(INOUT)         :: nat(maxpasym)
  ! complex(kind=dbl),INTENT(INOUT)         :: Sghvec(numset)
  
  ! integer(kind=irg)                       :: ip, ir, ic, ikk, n
  ! real(kind=sgl)                          :: Znsq, DBWF, kkl
  ! complex(kind=dbl)                       :: carg
  ! real(kind=dbl),allocatable              :: ctmp(:,:), apdata(:,:), apos(:,:,:)
  ! real(kind=dbl)                          :: arg, tpi
  ! integer(kind=irg), allocatable          :: numat(:), atp(:)
  
  !   tpi = 2.D0 * cPi
  !   Sghvec = cmplx(0.D0,0.D0)
  
  !   apdata = cell%getAsymPosData()
  !   numat = cell%getnumat()
  !   atp = cell%getatomtype()
  !   apos = cell%getapos()
  
  ! ! for each special position we need to compute its contribution to the Sgh array
  !   do ip=1,cell%getNatomtype()
  !     call SG%CalcOrbit(apdata(ip,1:3),n,ctmp)
  !     nat(ip) = numat(ip)
  ! ! get Zn-squared for this special position, and include the site occupation parameter as well
  !     Znsq = float(atp(ip))**2 * apdata(ip,4)
  
  ! ! multiply with the prefactor and add
  !   Sghvec = Sghvec + cmplx(n * Znsq, 0.D0)
  !   end do
  
  ! end subroutine preCalcSghECCI_
  
  
  
  
  ! ! !
  ! ! ! ###################################################################
  ! ! !
  ! ! !  subroutine CalcFresnelPropagator
  ! ! !
  ! ! !                                    created: 4/16/97
  ! ! !  Author: Marc De Graef
  ! ! !
  ! ! !  Description: compute the Fresnel propagator (for possibly inclined
  ! ! !               illumination) and store it in a file
  ! ! !
  ! ! !  History
  ! ! !
  ! ! !  modified by  rev reason
  ! ! !  -------- --- --- -----------
  ! ! !   4/16/97 MDG 1.0 original
  ! ! !   9/29/01 MDG 2.0 f90
  ! ! !  11/27/01 MDG 2.1 added kind support
  ! ! ! ###################################################################
  ! ! recursive subroutine CalcFresnelPropagator(beam,dimi,dimj,dz,scl,propname,lambda)
  ! ! !DEC$ ATTRIBUTES DLLEXPORT :: CalcFresnelPropagator
  
  ! ! use constants
  ! ! use io
  ! ! use files
  
  ! ! IMPLICIT NONE
  
  ! ! real(kind=sgl)                  :: beam(3),b,bm(2),dz,fidim,fjdim,prefac,scl, oi_real(2), lambda
  ! ! real(kind=sgl),allocatable      :: idimi(:),jdimj(:)
  ! ! complex(kind=sgl),allocatable   :: fr(:,:)
  ! ! integer(kind=irg)               :: dimi,dimj,i,ix,iy
  ! ! character(fnlen)                   :: propname
  
  ! ! INTENT(IN) :: beam,dimi,dimj,dz
  
  ! !   fidim = 1.0/float(dimi)
  ! !   fjdim = 1.0/float(dimj)
  ! !   prefac = scl*cPi*lambda*dz
  ! !   call Message('Computing Fresnel propagator', frm = "(A)")
  ! ! ! normalize the incident beam direction and rescale to the wavevector
  ! !   b = sqrt(sum(beam**2))
  ! !   bm= beam(1:2)/b/lambda
  ! !   oi_real(1:2) = bm(1:2)
  ! !   call WriteValue(' Laue center at ', oi_real, 2, "(F8.4,',',F8.4)")
  ! ! ! allocate variables
  ! !   allocate(fr(dimi,dimj))
  ! !   allocate(idimi(dimi),jdimj(dimj))
  ! !   idimi = float((/ (i, i=0,dimi-1) /))
  ! !   jdimj = float((/ (i, i=0,dimj-1) /))
  ! ! !
  ! !   where(idimi.ge.dimi/2) idimi = idimi-float(dimi)
  ! !   where(jdimj.ge.dimj/2) jdimj = jdimj-float(dimj)
  ! ! !
  ! !   idimi = prefac*idimi*fidim
  ! !   jdimj = prefac*jdimj*fjdim
  ! ! !
  ! !   idimi = idimi*(idimi + 2.0*bm(1))
  ! !   jdimj = jdimj*(jdimj + 2.0*bm(2))
  ! ! ! loop over y axis
  ! !   do iy=1,dimj
  ! ! ! loop over x-axis
  ! !     do ix=1,dimi
  ! !       fr(ix,iy)=cmplx(cos(idimi(ix)+jdimj(iy)),sin(idimi(ix)+jdimj(iy)))
  ! !     end do
  ! !   end do
  ! !   deallocate(idimi,jdimj)
  ! ! ! and store it in a file
  ! !   open(unit=dataunit,file=trim(EMsoft_toNativePath(propname)),form='unformatted')
  ! !   write (dataunit) dimi,dimj
  ! !   write (dataunit) fr
  ! !   close(unit=dataunit,status='keep')
  ! !   deallocate(fr)
  ! ! end subroutine
  
  
  ! ! !--------------------------------------------------------------------------
  ! ! !
  ! ! ! FUNCTION: LorentzPF
  ! ! !
  ! ! !> @author Marc De Graef, Carnegie Mellon University
  ! ! !
  ! ! !> @brief compute the Lorentz Polarization Factor Lp
  ! ! !
  ! ! !> @param theta scattering angle
  ! ! !> @param HEDM optional string to indicate HEDM mode
  ! ! !
  ! ! !> @date  03/26/13  MDG  1.0 added for HEDM project
  ! ! !--------------------------------------------------------------------------
  ! ! recursive function LorentzPF(theta,HEDM) result(tt)
  ! ! !DEC$ ATTRIBUTES DLLEXPORT :: LorentzPF
  
  ! ! use crystal
  
  ! ! IMPLICIT NONE
  
  ! ! real(kind=sgl),INTENT(IN)                       :: theta                !< scattering angle
  ! ! character(*),INTENT(IN),OPTIONAL                :: HEDM         !< for HEDM we have a different polarization factor
  ! ! real(kind=sgl)                                  :: tt
  
  ! ! if (present(HEDM)) then
  ! !   tt = (1.0+cos(2.0*theta)**2) / sin(theta)**2 / cos(theta)
  ! ! else
  ! !   tt = (1.0+cos(2.0*theta)**2) / sin(theta)**2 / cos(theta)
  ! ! end if
  
  ! ! end function
  
  
  
  end module mod_QCdiffraction
  
