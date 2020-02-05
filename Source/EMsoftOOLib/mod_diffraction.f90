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

!--------------------------------------------------------------------------
! EMsoft:diffraction.f90
!--------------------------------------------------------------------------
!
! MODULE: diffraction
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Anything related to dynamical diffraction
!
!> @todo general cleanup; merge with dynamical module?; add Private/Public 
! 
!> @date   10/13/98 MDG 1.0 original
!> @date    5/22/01 MDG 2.0 f90
!> @date   11/27/01 MDG 2.1 added kind support
!> @date    3/14/02 MDG 2.2 added CalcDynMat routine
!> @date   01/10/14 MDG 3.0 update with new cell type etc...
!> @date   12/01/14 MDG 3.1 removal of all global variables, including mLambda etc...
!--------------------------------------------------------------------------
module mod_diffraction

use mod_kinds
use mod_global 

IMPLICIT NONE

private 

! atomic scattering factor parametrization (Doyle-Turner, Smith-Burge)
! used only if absorption is not taken into account;  otherwise
! the Weickenmeier-Kohl routine is used.
real(kind=sgl),parameter :: scatfac(8,98) = reshape( (/ &
        0.202,0.244,0.082,0.000,0.30868,0.08544,0.01273,0.00000, &
        0.091,0.181,0.110,0.036,0.18183,0.06212,0.01803,0.00284, &
        1.611,1.246,0.326,0.099,1.07638,0.30480,0.04533,0.00495, &
        1.250,1.334,0.360,0.106,0.60804,0.18591,0.03653,0.00416, &
        0.945,1.312,0.419,0.116,0.46444,0.14178,0.03223,0.00377, &
        0.731,1.195,0.456,0.125,0.36995,0.11297,0.02814,0.00346, &
        0.572,1.043,0.465,0.131,0.28847,0.09054,0.02421,0.00317, &
        0.455,0.917,0.472,0.138,0.23780,0.07622,0.02144,0.00296, &
        0.387,0.811,0.475,0.146,0.20239,0.06609,0.01931,0.00279, &
        0.303,0.720,0.475,0.153,0.17640,0.05860,0.01762,0.00266, &
        2.241,1.333,0.907,0.286,1.08004,0.24505,0.03391,0.00435, &
        2.268,1.803,0.839,0.289,0.73670,0.20175,0.03013,0.00405, &
        2.276,2.428,0.858,0.317,0.72322,0.19773,0.03080,0.00408, &
        2.129,2.533,0.835,0.322,0.57775,0.16476,0.02880,0.00386, &
        1.888,2.469,0.805,0.320,0.44876,0.13538,0.02642,0.00361, &
        1.659,2.386,0.790,0.321,0.36650,0.11488,0.02469,0.00340, &
        1.452,2.292,0.787,0.322,0.30935,0.09980,0.02234,0.00323, &
        1.274,2.190,0.793,0.326,0.26682,0.08813,0.02219,0.00307, &
        3.951,2.545,1.980,0.482,1.37075,0.22402,0.04532,0.00434, &
        4.470,2.971,1.970,0.482,0.99523,0.22696,0.04195,0.00417, &
        3.966,2.917,1.925,0.480,0.88960,0.20606,0.03856,0.00399, &
        3.565,2.818,1.893,0.483,0.81982,0.19049,0.03590,0.00386, &
        3.245,2.698,1.860,0.486,0.76379,0.17726,0.03363,0.00374, &
        2.307,2.334,1.823,0.490,0.78405,0.15785,0.03157,0.00364, &
        2.747,2.456,1.792,0.498,0.67786,0.15674,0.03000,0.00357, &
        2.544,2.343,1.759,0.506,0.64424,0.14880,0.02854,0.00350, &
        2.367,2.236,1.724,0.515,0.61431,0.14180,0.02725,0.00344, &
        2.210,2.134,1.689,0.524,0.58727,0.13553,0.02609,0.00339, &
        1.579,1.820,1.658,0.532,0.62940,0.12453,0.02504,0.00333, &
        1.942,1.950,1.619,0.543,0.54162,0.12518,0.02416,0.00330, &
        2.321,2.486,1.688,0.599,0.65602,0.15458,0.02581,0.00351, &
        2.447,2.702,1.616,0.601,0.55893,0.14393,0.02446,0.00342, &
        2.399,2.790,1.529,0.594,0.45718,0.12817,0.02280,0.00328, &
        2.298,2.854,1.456,0.590,0.38830,0.11536,0.02146,0.00316, &
        2.166,2.904,1.395,0.589,0.33899,0.10497,0.02041,0.00307, &
        2.034,2.927,1.342,0.589,0.29999,0.09598,0.01952,0.00299, &
        4.776,3.859,2.234,0.868,1.40782,0.18991,0.03701,0.00419, &
        5.848,4.003,2.342,0.880,1.04972,0.19367,0.03737,0.00414, &
        4.129,3.012,1.179,0.000,0.27548,0.05088,0.00591,0.000,   &
        4.105,3.144,1.229,0.000,0.28492,0.05277,0.00601,0.000,   &
        4.237,3.105,1.234,0.000,0.27415,0.05074,0.00593,0.000,   &
        3.120,3.906,2.361,0.850,0.72464,0.14642,0.03237,0.00366, &
        4.318,3.270,1.287,0.000,0.28246,0.05148,0.00590,0.000,   &
        4.358,3.298,1.323,0.000,0.27881,0.05179,0.00594,0.000,   &
        4.431,3.343,1.345,0.000,0.27911,0.05153,0.00592,0.000,   &
        4.436,3.454,1.383,0.000,0.28670,0.05269,0.00595,0.000,   &
        2.036,3.272,2.511,0.837,0.61497,0.11824,0.02846,0.00327, &
        2.574,3.259,2.547,0.838,0.55675,0.11838,0.02784,0.00322, &
        3.153,3.557,2.818,0.884,0.66649,0.14449,0.02976,0.00335, &
        3.450,3.735,2.118,0.877,0.59104,0.14179,0.02855,0.00327, &
        3.564,3.844,2.687,0.864,0.50487,0.13316,0.02691,0.00316, &
        4.785,3.688,1.500,0.000,0.27999,0.05083,0.00581,0.000,   &
        3.473,4.060,2.522,0.840,0.39441,0.11816,0.02415,0.00298, &
        3.366,4.147,2.443,0.829,0.35509,0.11117,0.02294,0.00289, &
        6.062,5.986,3.303,1.096,1.55837,0.19695,0.03335,0.00379, &
        7.821,6.004,3.280,1.103,1.17657,0.18778,0.03263,0.00376, &
        4.940,3.968,1.663,0.000,0.28716,0.05245,0.00594,0.000,   &
        5.007,3.980,1.678,0.000,0.28283,0.05183,0.00589,0.000,   &
        5.085,4.043,1.684,0.000,0.28588,0.05143,0.00581,0.000,   &
        5.151,4.075,1.683,0.000,0.28304,0.05073,0.00571,0.000,   &
        5.201,4.094,1.719,0.000,0.28079,0.05081,0.00576,0.000,   &
        5.255,4.113,1.743,0.000,0.28016,0.05037,0.00577,0.000,   &
        6.267,4.844,3.202,1.200,1.00298,0.16066,0.02980,0.00367, &
        5.225,4.314,1.827,0.000,0.29158,0.05259,0.00586,0.000,   &
        5.272,4.347,1.844,0.000,0.29046,0.05226,0.00585,0.000,   &
        5.332,4.370,1.863,0.000,0.28888,0.05198,0.00581,0.000,   &
        5.376,4.403,1.884,0.000,0.28773,0.05174,0.00582,0.000,   &
        5.436,4.437,1.891,0.000,0.28655,0.05117,0.00577,0.000,   &
        5.441,4.510,1.956,0.000,0.29149,0.05264,0.00590,0.000,   &
        5.529,4.533,1.945,0.000,0.28927,0.05144,0.00578,0.000,   &
        5.553,4.580,1.969,0.000,0.28907,0.05160,0.00577,0.000,   &
        5.588,4.619,1.997,0.000,0.29001,0.05164,0.00579,0.000,   &
        5.659,4.630,2.014,0.000,0.28807,0.05114,0.00578,0.000,   &
        5.709,4.677,2.019,0.000,0.28782,0.05084,0.00572,0.000,   &
        5.695,4.740,2.064,0.000,0.28968,0.05156,0.00575,0.000,   &
        5.750,4.773,2.079,0.000,0.28933,0.05139,0.00573,0.000,   &
        5.754,4.851,2.096,0.000,0.29159,0.05152,0.00570,0.000,   &
        5.803,4.870,2.127,0.000,0.29016,0.05150,0.00572,0.000,   &
        2.388,4.226,2.689,1.255,0.42866,0.09743,0.02264,0.00307, &
        2.682,4.241,2.755,1.270,0.42822,0.09856,0.02295,0.00307, &
        5.932,4.972,2.195,0.000,0.29086,0.05126,0.00572,0.000,   &
        3.510,4.552,3.154,1.359,0.52914,0.11884,0.02571,0.00321, &
        3.841,4.679,3.192,1.363,0.50261,0.11999,0.02560,0.00318, &
        6.070,4.997,2.232,0.000,0.28075,0.04999,0.00563,0.000,   &
        6.133,5.031,2.239,0.000,0.28047,0.04957,0.00558,0.000,   &
        4.078,4.978,3.096,1.326,0.38406,0.11020,0.02355,0.00299, &
        6.201,5.121,2.275,0.000,0.28200,0.04954,0.00556,0.000,   &
        6.215,5.170,2.316,0.000,0.28382,0.05002,0.00562,0.000,   &
        6.278,5.195,2.321,0.000,0.28323,0.04949,0.00557,0.000,   &
        6.264,5.263,2.367,0.000,0.28651,0.05030,0.00563,0.000,   &
        6.306,5.303,2.386,0.000,0.28688,0.05026,0.00561,0.000,   &
        6.767,6.729,4.014,1.561,0.85951,0.15642,0.02936,0.00335, &
        6.323,5.414,2.453,0.000,0.29142,0.05096,0.00568,0.000,   &
        6.415,5.419,2.449,0.000,0.28836,0.05022,0.00561,0.000,   &
        6.378,5.495,2.495,0.000,0.29156,0.05102,0.00565,0.000,   &
        6.460,5.469,2.471,0.000,0.28396,0.04970,0.00554,0.000,   &
        6.502,5.478,2.510,0.000,0.28375,0.04975,0.00561,0.000,   &
        6.548,5.526,2.520,0.000,0.28461,0.04965,0.00557,0.000/), (/8,98/))


! The parameters in gnode are computed by CalcUcg 
type, public :: gnode
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
end type gnode


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
type, public :: Diffraction_T
private 

    real(kind=dbl)                  :: voltage          ! voltage always in keV !
    real(kind=dbl)                  :: Lambda           ! wave length in nm 
    real(kind=dbl)                  :: Relcor
    real(kind=dbl)                  :: Sigma
    real(kind=dbl)                  :: Psihat   
    real(kind=dbl)                  :: V0mod = 0.D0     ! used to compute refraction correction 
    type(gnode)                     :: rlp              ! default variable for reciprocal lattice point 
    type(DynType)                   :: Dyn              ! dynamical scattering parameters 
    type(BetheParameterType)        :: BetheParameters
    real(kind=sgl),allocatable      :: scatfacg(:)
    complex(kind=sgl),allocatable   :: scatfac(:,:) 
    integer(kind=irg)               :: numscatfac
    complex(kind=sgl),allocatable   :: LUT(:,:,:), SghLUT(:,:,:,:)
    complex(kind=sgl),allocatable   :: LUTqg(:,:,:)
    logical,allocatable             :: dbdiff(:,:,:)

contains 
private 

    procedure, pass(self) :: GetVoltage_
    procedure, pass(self) :: CalcWaveLength_
    procedure, pass(self) :: CalcDiffAngle_
    procedure, pass(self) :: CalcUcg_
    procedure, pass(self) :: PreCalcFSCATT_
    procedure, pass(self) :: CalcsgSingle_
    procedure, pass(self) :: CalcsgDouble_
    procedure, pass(self) :: getV_
    procedure, pass(self) :: getScatfac_
    procedure, pass(self) :: getWaveLength_
    procedure, pass(self) :: getBetheParameter_
    procedure, pass(self) :: Set_Bethe_Parameters_
    procedure, pass(self) :: BWsolve_
    procedure, pass(self) :: setrlpmethod_
    procedure, pass(self) :: getrlp_
    procedure, pass(self) :: allocateLUT_
    procedure, pass(self) :: getshapeLUT_
    procedure, pass(self) :: getLUT_
    procedure, pass(self) :: getLUTqg_
    procedure, pass(self) :: getdbdiff_
    procedure, pass(self) :: setLUT_
    procedure, pass(self) :: setLUTqg_
    procedure, pass(self) :: setdbdiff_
    procedure, pass(self) :: Printrlp_

    generic, public :: GetVoltage => GetVoltage_
    generic, public :: CalcWaveLength => CalcWaveLength_
    generic, public :: CalcDiffAngle => CalcDiffAngle_
    generic, public :: CalcUcg => CalcUcg_
    generic, public :: PreCalcFSCATT => PreCalcFSCATT_
    generic, public :: Calcsg => CalcsgSingle_, CalcsgDouble_
    generic, public :: getScatfac => getScatfac_
    generic, public :: getWaveLength => getWaveLength_
    generic, public :: getBetheParameter => getBetheParameter_
    generic, public :: SetBetheParameters => Set_Bethe_Parameters_
    generic, public :: BWsolve => BWsolve_
    generic, public :: setrlpmethod => setrlpmethod_
    generic, public :: getV => getV_
    generic, public :: getrlp => getrlp_
    generic, public :: allocateLUT => allocateLUT_
    generic, public :: getshapeLUT => getshapeLUT_
    generic, public :: getLUT => getLUT_
    generic, public :: getLUTqg => getLUTqg_
    generic, public :: getdbdiff => getdbdiff_
    generic, public :: setLUT => setLUT_
    generic, public :: setLUTqg => setLUTqg_
    generic, public :: setdbdiff => setdbdiff_
    generic, public :: Printrlp => Printrlp_

end type Diffraction_T

! the constructor routine for this class 
  interface Diffraction_T
    module procedure Diffraction_constructor
  end interface Diffraction_T

!DEC$ ATTRIBUTES DLLEXPORT :: GetVoltage
!DEC$ ATTRIBUTES DLLEXPORT :: CalcWaveLength
!DEC$ ATTRIBUTES DLLEXPORT :: CalcDiffAngle
!DEC$ ATTRIBUTES DLLEXPORT :: CalcUcg
!DEC$ ATTRIBUTES DLLEXPORT :: PreCalcFSCATT
!DEC$ ATTRIBUTES DLLEXPORT :: CalcsgSingle
!DEC$ ATTRIBUTES DLLEXPORT :: CalcsgDouble
!DEC$ ATTRIBUTES DLLEXPORT :: getV
!DEC$ ATTRIBUTES DLLEXPORT :: getScatfac
!DEC$ ATTRIBUTES DLLEXPORT :: getWaveLength
!DEC$ ATTRIBUTES DLLEXPORT :: getBetheParameter
!DEC$ ATTRIBUTES DLLEXPORT :: Set_Bethe_Parameters
!DEC$ ATTRIBUTES DLLEXPORT :: BWsolve
!DEC$ ATTRIBUTES DLLEXPORT :: setrlpmethod
!DEC$ ATTRIBUTES DLLEXPORT :: getrlp
!DEC$ ATTRIBUTES DLLEXPORT :: getshapeLUT
!DEC$ ATTRIBUTES DLLEXPORT :: allocateLUT
!DEC$ ATTRIBUTES DLLEXPORT :: getLUT
!DEC$ ATTRIBUTES DLLEXPORT :: getLUTqg
!DEC$ ATTRIBUTES DLLEXPORT :: getdbdiff
!DEC$ ATTRIBUTES DLLEXPORT :: setLUT
!DEC$ ATTRIBUTES DLLEXPORT :: setLUTqg
!DEC$ ATTRIBUTES DLLEXPORT :: setdbdiff
!DEC$ ATTRIBUTES DLLEXPORT :: Printrlp


! ! some of these go in the PostScript module 

!     procedure, pass(self) :: DiffPage
!     procedure, pass(self) :: DumpZAP
!     procedure, pass(self) :: DumpPP
!     procedure, pass(self) :: studylist

! ! Lorentz module ?
!     procedure, pass(self) :: LorentzPF



contains

!--------------------------------------------------------------------------
type(Diffraction_T) function Diffraction_constructor( voltage, cell, verbose ) result(Diff)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/26/20
  !!
  !! constructor for the Diffraction Class 
  
use mod_crystallography 

IMPLICIT NONE

real(kind=dbl), INTENT(IN)          :: voltage    ! in kV 
type(Cell_T), INTENT(INOUT)         :: cell 
logical, INTENT(IN), OPTIONAL       :: verbose

integer(kind=irg)                   :: hkl(3)

Diff%voltage = voltage 

hkl=(/0,0,0/)
call Diff%CalcUcg(cell, hkl) 
Diff%V0mod = Diff%rlp%Vmod

if (present(verbose)) then 
    if (verbose.eqv..TRUE.) call Diff%CalcWaveLength( cell, verbose )
else
    call Diff%CalcWaveLength( cell )
end if 

end function Diffraction_constructor

!--------------------------------------------------------------------------
recursive subroutine CalcWaveLength_(self, cell, verbose)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/26/20
  !!
  !! compute the relativistic electron wavelength
  !!
  !! These quantities are computed in double precision because of the 
  !! wide range of magnitudes.  If a crystal structure has been defined
  !! then the gamma*V_0 term is added to correct for refraction.

use mod_symmetry
use mod_crystallography
use mod_io

IMPLICIT NONE

class(Diffraction_T),INTENT(INOUT)      :: self
type(Cell_T),INTENT(INOUT)              :: cell
logical,INTENT(IN),OPTIONAL             :: verbose

type(IO_T)                              :: Message
real(kind=dbl)                          :: temp1, temp2, oi_real(1)


  temp1 = 1.0D+9*cPlanck/dsqrt(2.D0*cRestmass*cCharge)
  temp2 = cCharge*0.5D0*self%voltage*1000.D0/cRestmass/(cLight**2)

! relativistic correction factor (known as gamma)      
  self%Relcor = 1.0D0+2.0D0*temp2

! relativistic acceleration voltage
  self%Psihat = self%voltage*(1.D0+temp2)*1000.D0

! compute the electron wavelength in nm
! compute V_0 and add it to mPsihat (corrected by mRelcor)
! this should have been done in the calling program 
!     call cell%CalcPositions('v')
 if (self%V0mod.eq.0.D0) then 
    call self%CalcUcg(cell, (/0,0,0/) ) 
    self%V0mod = self%rlp%Vmod
 end if 

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

end subroutine CalcWaveLength_

!--------------------------------------------------------------------------
recursive subroutine setrlpmethod_(self, m, ask )
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/26/20
  !!
  !! which scattering factors should be used ?

use mod_io 

IMPLICIT NONE 

class(Diffraction_T),INTENT(INOUT)      :: self
character(2), INTENT(IN), OPTIONAL      :: m
logical, INTENT(IN), OPTIONAL           :: ask 

type(IO_T)                              :: Message
integer(kind=irg)                       :: io_int(1) 

if (present(m)) then    ! we ignore the 'ask' parameter 
   self%rlp%absorption = .FALSE.
   self%rlp%method=m 
   if (m.eq.'WK') self%rlp%absorption = .TRUE.
else 
  if (present(ask)) then
   call Message%printMessage(  (/ & 
            ' The following scattering factor sets are available :', &
            '  [1] Doyle-Turner/Smith-Burge (no absorption)       ', & 
            '  [2] Weickenmeier-Kohl (no absorption)              ', &
            '  [3] Weickenmeier-Kohl (with absorption)            ' /) )
   call Message%ReadValue(' Which set do you want to use [1/2/3] ? ', io_int,1)
   self%rlp%absorption = .FALSE.
   select case (io_int(1)) 
    case(1); self%rlp%method='DT'; 
    case(2); self%rlp%method='WK'; 
    case(3); self%rlp%method='WK'; self%rlp%absorption=.TRUE.
   end select
  end if
end if 

end subroutine setrlpmethod_

!--------------------------------------------------------------------------
recursive function getrlp_(self) result(rlp)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/26/20
  !!
  !! ask for accelerating voltage, then call CalcWaveLength

IMPLICIT NONE 

class(Diffraction_T),INTENT(INOUT)  :: self
type(gnode)                         :: rlp 

rlp = self%rlp 

end function getrlp_

!--------------------------------------------------------------------------
recursive subroutine allocateLUT_(self, dims)
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/04/20
  !!
  !! allocates the LUT arrays

use mod_io 

IMPLICIT NONE 

class(Diffraction_T),INTENT(INOUT)  :: self
integer(kind=irg)                   :: dims(3) 

integer(kind=irg)                   :: imh, imk, iml, istat
type(IO_T)                          :: Message

imh = dims(1)
imk = dims(2)
iml = dims(3)

! the LUT array stores all the Fourier coefficients, so that we only need to compute them once... i.e., here and now
allocate(self%LUT(-2*imh:2*imh,-2*imk:2*imk,-2*iml:2*iml),stat=istat)
if (istat.ne.0) call Message%printError('allocateLUT:',' unable to allocate LUT array')

allocate(self%LUTqg(-2*imh:2*imh,-2*imk:2*imk,-2*iml:2*iml),stat=istat)
if (istat.ne.0) call Message%printError('allocateLUT:',' unable to allocate LUTqg array')
 
! allocate an array that keeps track of potential double diffraction reflections
allocate(self%dbdiff(-2*imh:2*imh,-2*imk:2*imk,-2*iml:2*iml),stat=istat)
if (istat.ne.0) call Message%printError('allocateLUT:',' unable to allocate dbdiff array')

self%LUT = cmplx(0.D0,0.D0)
self%LUTqg = cmplx(0.D0,0.D0)
self%dbdiff = .FALSE.

end subroutine allocateLUT_

!--------------------------------------------------------------------------
recursive function getshapeLUT_(self) result(s)
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/04/20
  !!
  !! returns the shape of the LUT array

IMPLICIT NONE 

class(Diffraction_T),INTENT(INOUT)  :: self
integer(kind=irg)                   :: s(3) 

 s = shape(self%LUT)

end function getshapeLUT_

!--------------------------------------------------------------------------
recursive function getLUT_(self, hkl) result(c)
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/04/20
  !!
  !! return an entry from the LUT array

IMPLICIT NONE 

class(Diffraction_T),INTENT(INOUT)  :: self
integer(kind=irg)                   :: hkl(3)
complex(kind=dbl)                   :: c 

 c = self%LUT( hkl(1), hkl(2), hkl(3) )

end function getLUT_

!--------------------------------------------------------------------------
recursive function getLUTqg_(self, hkl) result(c)
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/04/20
  !!
  !! return an entry from the LUTqg array

IMPLICIT NONE 

class(Diffraction_T),INTENT(INOUT)  :: self
integer(kind=irg)                   :: hkl(3)
complex(kind=dbl)                   :: c 

 c = self%LUTqg( hkl(1), hkl(2), hkl(3) )

end function getLUTqg_

!--------------------------------------------------------------------------
recursive function getdbdiff_(self, hkl) result(c)
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/04/20
  !!
  !! return an entry from the dbdiff array

IMPLICIT NONE 

class(Diffraction_T),INTENT(INOUT)  :: self
integer(kind=irg)                   :: hkl(3)
logical                             :: c

 c = self%dbdiff( hkl(1), hkl(2), hkl(3) )

end function getdbdiff_

!--------------------------------------------------------------------------
recursive subroutine setLUT_(self, hkl, c)
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/04/20
  !!
  !! set an entry in the LUT array

IMPLICIT NONE 

class(Diffraction_T),INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)        :: hkl(3)
complex(kind=sgl),INTENT(IN)        :: c 

self%LUT( hkl(1), hkl(2), hkl(3) ) = c

end subroutine setLUT_

!--------------------------------------------------------------------------
recursive subroutine setLUTqg_(self, hkl, c)
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/04/20
  !!
  !! set an entry in the LUTqg array

IMPLICIT NONE 

class(Diffraction_T),INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)        :: hkl(3)
complex(kind=sgl),INTENT(IN)        :: c 

self%LUTqg( hkl(1), hkl(2), hkl(3) ) = c

end subroutine setLUTqg_

!--------------------------------------------------------------------------
recursive subroutine setdbdiff_(self, hkl, c)
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/04/20
  !!
  !! set an entry in the dbdiff array

IMPLICIT NONE 

class(Diffraction_T),INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)        :: hkl(3)
logical,INTENT(IN)                  :: c

self%dbdiff( hkl(1), hkl(2), hkl(3) ) = c

end subroutine setdbdiff_

!--------------------------------------------------------------------------
recursive function getrlpmethod_(self) result(rlpmethod)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/26/20
  !!
  !! retrieve the rlp method

IMPLICIT NONE 

class(Diffraction_T),INTENT(INOUT)  :: self
character(2)                        :: rlpmethod 

rlpmethod = self%rlp%method

end function getrlpmethod_

!--------------------------------------------------------------------------
recursive function getV_(self) result(V)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/28/20
  !!
  !! return the accelerating voltage

IMPLICIT NONE 

class(Diffraction_T),INTENT(INOUT)  :: self
real(kind=dbl)                      :: V

V = self%voltage

end function getV_

!--------------------------------------------------------------------------
recursive subroutine Printrlp_(self,first)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! output the contents of the rlp structure

use mod_io

IMPLICIT NONE

class(Diffraction_T),INTENT(INOUT)      :: self 
logical,optional,intent(INOUT)          :: first                
 !! switch for long/short output

type(IO_T)                              :: Message 
integer(kind=irg)                       :: oi_int(3)
real(kind=sgl)                          :: oi_real(7)
complex(kind=sgl)                       :: oi_cmplx(1)


if (present(first)) then
 if (first) then
  call Message%printMessage('     Scattering factors : ', frm = "(/A)",advance="no")
  if (self%rlp%method.eq.'WK') then 
   if (self%rlp%absorption.eqv..TRUE.) then 
    call Message%printMessage(' Weickenmeier-Kohl (with absorption)', frm = "(A/)")
   else
    call Message%printMessage(' Weickenmeier-Kohl', frm = "(A/)")
   end if
  else
    call Message%printMessage(' Doyle-Turner/Smith-Burge', frm = "(A/)")
  end if

  if (self%rlp%absorption.eqv..TRUE.) then
    call Message%printMessage( &
        '   h  k  l    |g|    Ucg_r  Ucg_i   |Ug|    phase   |Ugp|   phase   xi_g   xi_gp    ratio  Re-1/q_g-Im', &
        frm = "(A)")
  else
    call Message%printMessage('   h  k  l    |g|    Ucg_r  |Ug|    phase    xi_g   1/q_g', frm = "(A)")
  end if
  first = .FALSE.
 end if
end if

if (self%rlp%absorption.eqv..TRUE.) then
 oi_int(1:3) = self%rlp%hkl(1:3)
 call Message%WriteValue('',oi_int, 3, "(1x,3I3,1x)",advance="no")
 oi_real(1) = self%rlp%g
 call Message%WriteValue('',oi_real, 1, "(F9.4)",advance="no")
 oi_cmplx(1) = self%rlp%Ucg
 call Message%WriteValue('',oi_cmplx, 1, "(2F7.3,1x)",advance="no")
 oi_real(1:7)  = (/ self%rlp%Umod,self%rlp%Vphase*180.0/sngl(cPi),self%rlp%Upmod,self%rlp%Vpphase*180.0/sngl(cPi),&
                    self%rlp%xg,self%rlp%xgp,self%rlp%ar /)
 call Message%WriteValue('',oi_real, 7, "(4F8.3,3F8.1)",advance="no")
 oi_cmplx(1) = self%rlp%qg
 call Message%WriteValue('',oi_cmplx, 1, "(2F8.3)")
else
 oi_int(1:3) = self%rlp%hkl(1:3)
 call Message%WriteValue('',oi_int, 3, "(1x,3I3,1x)",advance="no")
 oi_real(1) = self%rlp%g
 call Message%WriteValue('',oi_real, 1, "(F9.4)",advance="no")
 oi_real(1) = real(self%rlp%Ucg)
 call Message%WriteValue('',oi_real, 1, "(F7.3,1x)",advance="no")
 oi_real(1:3)  = (/ self%rlp%Umod,self%rlp%Vphase*180.0/sngl(cPi),self%rlp%xg /)
 call Message%WriteValue('',oi_real, 3, "(2F8.3,F8.1)",advance="no")
 oi_cmplx(1) = self%rlp%qg
 call Message%WriteValue('',oi_cmplx, 1, "(2F8.3)")
end if

end subroutine Printrlp_

!--------------------------------------------------------------------------
recursive subroutine GetVoltage_(self, cell, verbose)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/26/20
  !!
  !! ask for accelerating voltage, then call CalcWaveLength

use mod_io
use mod_crystallography

IMPLICIT NONE

class(Diffraction_T),INTENT(INOUT)  :: self
type(Cell_T),INTENT(INOUT)          :: cell 
logical,INTENT(IN),OPTIONAL         :: verbose


real(kind=sgl)                      :: io_real(1)
type(IO_T)                          :: Message 

call Message%ReadValue(' Enter the microscope accelerating voltage [kV, R] : ', io_real, 1)
self%voltage = dble(io_real(1))

if (present(verbose)) then
    if (verbose.eqv..TRUE.) call self%CalcWaveLength(cell, verbose = .TRUE.)
else
    call self%CalcWaveLength(cell)
end if

end subroutine GetVoltage_

!--------------------------------------------------------------------------
recursive function CalcDiffAngle_(self, cell, hkl) result(tt)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/26/20
  !!
  !! compute the diffraction angle 2theta in radians

use mod_crystallography

IMPLICIT NONE

class(Diffraction_T),INTENT(INOUT)  :: self
type(Cell_T),INTENT(IN)             :: cell
integer(kind=irg),INTENT(IN)        :: hkl(3)                
 !! Miller indices

real(kind=sgl)                      :: tt

tt = 2.0*asin(0.50*sngl(self%Lambda)*cell%CalcLength(float(hkl), 'r') )

end function CalcDiffAngle_

!--------------------------------------------------------------------------
recursive subroutine CalcUcg_(self, cell, hkl, applyqgshift, interpolate)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/26/20
  !!
  !! compute the complex Structure Factor for a given g
  !!
  !! includes computation of extinction distance, absorption length, etc...
  !! This routine is probably the most important one in all the dynamical routines,
  !! because it computes all possible relevant parameters and stores them in the rlp variable.
  !! We've added the XR rlp%method parameter so that the same routine can be used for the
  !! computation of kinematical x-ray scattering; this was needed for the HEDM package.

use mod_crystallography
use mod_symmetry
use mod_others

IMPLICIT NONE

class(Diffraction_T),INTENT(INOUT)  :: self
type(Cell_T),INTENT(INOUT)          :: cell
integer(kind=irg),INTENT(IN)        :: hkl(3)               
 !! Miller indices
logical,OPTIONAL,INTENT(IN)         :: applyqgshift
 !! [optional] multiply qg by exp[i theta_g] if present and .TRUE.
logical,OPTIONAL,INTENT(IN)         :: interpolate          ! requires rlp%mode = 'IP'

integer(kind=irg)                   :: j,absflg,m,ii
real(kind=sgl)                      :: s,twopi,arg,swk,dwwk,pref,ul,pre,sct,fs,fsp
real(kind=sgl),parameter            :: preg = 0.664840340614319 ! = 2.0 * sngl(cRestmass*cCharge/cPlanck**2)*1.0E-18
complex(kind=sgl)                   :: ff,gg,sf,p1
complex(kind=sgl)                   :: czero
complex(kind=sgl),allocatable       :: sfarray(:)
integer(kind=irg),allocatable       :: atp(:)
real(kind=dbl), allocatable         :: apos(:,:), atpos(:,:,:)
integer(kind=irg), allocatable      :: numat(:)
logical                             :: accflg, dwflg, interp
character(2)                        :: smb

! interpolation is only used for the WK mode to pre-compute the scattering factor array
interp = .FALSE.
if (present(interpolate)) then
  if (interpolate.eqv..TRUE.) interp=.TRUE.
end if

twopi = sngl(2.D0*cPi)
czero = cmplx(0.0,0.0)
self%rlp%hkl = hkl

! compute the scattering parameter s^2=(g/2)^2
if (sum(hkl**2).eq.0) then 
  s = 0.0
  self%rlp%g = 0.0
else
  self%rlp%g = sngl(cell%CalcLength(dble(hkl),'r'))
  s = (0.50*self%rlp%g)**2
end if

atp = cell%getAtomtype()
apos = cell%getAsymPosData()
atpos = cell%getapos()
numat = cell%getnumat() 

!----------------------------------
! first the simplest case: kinematical X-ray scattering factors
! this option was added to accomodate the HEDM forward projector needs
if (self%rlp%method.eq.'XR') then 

! set the prefactor for the Doyle-Turner summation
  pref = 0.4178214

! initialize the real and imaginary parts of the structure factor
  sf = czero
  
! loop over all atoms in the asymmetric unit
 do m=1,cell%getNatomtype()

! get the atomic scattering factor for this atom
  sct = 0.0
  j = atp(m)
  do ii=1,4
   sct=sct+scatfac(ii,j)*exp(-scatfac(ii+4,j)*s)
  end do

! scale and include Debye-Waller factor and site occupation parameter
  fsp = pref * s * sct 
  fs = (float(atp(m)) - fsp ) * apos(m,4) * exp(-apos(m,5)*s)

! loop over all atoms in the orbit
  do j=1,numat(m)
   arg=twopi * sum(hkl(1:3)*atpos(m,j,1:3))
   sf = sf + fs * cmplx(cos(arg),-sin(arg))
  end do

 end do ! m loop

! and fill in just two entries of the rlp variable
 self%rlp%hkl = hkl
 self%rlp%Ucg = sf
 
end if

!----------------------------------
! for Doyle-Turner scattering factors
if (self%rlp%method.eq.'DT') then 
! initialize the real and imaginary parts of the structure factor
 sf = czero

! compute the prefactor (this also scales from Angstrom to nm)
 pref = 0.04787801/cell%getVolume()

! loop over all atoms in the asymmetric unit
 do m=1,cell%getNatomtype()

! get the atomic scattering factor for this atom
  sct=0.0
  j=atp(m)
  do ii=1,4
   sct=sct+scatfac(ii,j)*exp(-scatfac(ii+4,j)*s)
  end do

! scale and include Debye-Waller factor
! and site occupation parameter
  fs=pref*sct*exp(-apos(m,5)*s)*apos(m,4)

! loop over all atoms in the orbit
  do j=1,numat(m)
   arg=twopi*sum(hkl(1:3)*atpos(m,j,1:3))
   sf = sf + fs*exp(cmplx(0.0,-arg))
  end do
 end do ! m loop

! and fill in the entries of the rlp variable
 pre = 2.0*sngl(cRestmass*cCharge/cPlanck**2)*1.0E-18
 self%rlp%hkl = hkl
 self%rlp%Vmod = cabs(sf)*self%Relcor
 self%rlp%Vphase = atan2(aimag(sf),real(sf))
 self%rlp%Vpmod = 0.0
 self%rlp%Vpphase = 0.0
 if (self%rlp%Vmod.gt.0.0) then
  self%rlp%xg = 1.0/(pre*self%rlp%Vmod*self%Lambda)
 else
  self%rlp%xg = 1.0E+8
 end if
 self%rlp%xgp = 0.0
 self%rlp%ar = 0.0
 self%rlp%Vg = self%rlp%Vmod * exp(cmplx(0.0,self%rlp%Vphase))
 self%rlp%Ucg = pre*self%rlp%Vg
 self%rlp%qg = cmplx(1.0/self%rlp%xg,0.0)
end if


!----------------------------------
if (self%rlp%method.eq.'WK') then 
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
 s = self%rlp%g*swk

! let fscatt perform the relativistic corrections for f_g and fprime_g
 accflg = .TRUE.

! include absorption ?
 absflg = 0
 if (self%rlp%absorption.eqv..TRUE.) absflg = 3  ! include phonon and core contributions

! always include Debye-Waller factor
 dwflg  = .TRUE.

! compute the scaling prefactors
! pref contains A to nm conversion, and divides by 4pi
 pref = 0.04787801/cell%getVolume()/(4.0*cPi) 

! preg is used to go from V to U, remembering that gamma is already
! included in the output from fscatt
 pre = pref * preg

! initialize the real and imaginary parts of the structure factor
 ff=czero
 gg=czero

! loop over all atoms in the asymmetric unit
 do m=1,cell%getNatomtype()
! get the atomic scattering factor for this atom
! scale and include Debye-Waller factor and site occupation parameter
  ul = sqrt(apos(m,5)*dwwk)
  j = atp(m)
  sf = FSCATT(s,ul,j,smb,sngl(self%voltage),absflg,accflg,dwflg)*cmplx(apos(m,4),0.0)

! loop over all atoms in the orbit
  p1 = czero
  do j=1,numat(m)
   arg=twopi*sum(float(hkl(1:3))*atpos(m,j,1:3))
   p1 = p1 + exp(cmplx(0.0,-arg))
  end do

  ff = ff + p1*real(sf)
  gg = gg + p1*aimag(sf)

 end do
!
! fill in the entries of the rlp variable
 self%rlp%hkl = hkl

! these are the modulus and phase of the real part of Vg
 self%rlp%Vmod = pref * cabs(ff)
 self%rlp%Vphase = atan2(aimag(ff),real(ff))

! modulus of U_g
 self%rlp%Umod = preg*self%rlp%Vmod

! if absorption is included, also compute the imaginary part of Vg, i.e., Vprime_g
 if (self%rlp%absorption.eqv..TRUE.) then 
  self%rlp%Vpmod = pref * cabs(gg)
  self%rlp%Vpphase = atan2(aimag(gg),real(gg))

! modulus of Uprime_g
  self%rlp%Upmod = preg*self%rlp%Vpmod

! complex Ucg = U_g + i Uprime_g = U_g,r-Uprime_g,i + i(U_g,i+Uprime_g,r)
  self%rlp%Ucg = pre * cmplx(real(ff)-aimag(gg),aimag(ff)+real(gg))
 else ! set absorption parameters to zero
  self%rlp%Vpmod = 0.0
  self%rlp%Vpphase = 0.0
! Ucg = U_g (complex number)
  self%rlp%Ucg = pre * ff
 end if

! complex Vg 
 self%rlp%Vg = self%rlp%Ucg/preg
 if (abs(self%rlp%Umod).gt.0.0) then 
  self%rlp%xg = 1.0/abs(self%rlp%Umod)/self%Lambda
 else
  self%rlp%xg = 1.0E+8
 end if 

 if (abs(self%rlp%Upmod).gt.0.0) then 
  self%rlp%xgp = 1.0/abs(self%rlp%Upmod)/self%Lambda
 else
  self%rlp%xgp = 1.0E+8
 end if 

 if (self%rlp%absorption.eqv..TRUE.) then 
  self%rlp%ar = self%rlp%xgp/self%rlp%xg
  if (present(applyqgshift)) then
    if (applyqgshift.eqv..TRUE.) then
      self%rlp%qg = cmplx(cos(self%rlp%Vphase)/self%rlp%xg-sin(self%rlp%Vpphase)/self%rlp%xgp, &
                          cos(self%rlp%Vpphase)/self%rlp%xgp+sin(self%rlp%Vphase)/self%rlp%xg)
    end if
  else
    arg = self%rlp%Vpphase-self%rlp%Vphase
    self%rlp%qg = cmplx(1.0/self%rlp%xg-sin(arg)/self%rlp%xgp,cos(arg)/self%rlp%xgp)
  end if
 else
  self%rlp%ar = 0.0
  self%rlp%qg = cmplx(1.0/self%rlp%xg,0.0)
 end if

end if

!----------------------------------
if (self%rlp%method.eq.'IP') then 
 allocate(sfarray(cell%getNatomtype()))

! The Weickenmeier-Kohl (WK) scattering parameters have been pre-calculated 
! so all we need to do is linear interpolation to get the correct value
 swk = 0.1*twopi
  
! properly scale the scattering parameter
 s = self%rlp%g*swk

! get the atomic scattering factors for all atom types by linear interpolation
 call self%getScatfac(cell, s, sfarray, cell%getNatomtype())

! compute the scaling prefactors
! pref contains A to nm conversion, and divides by 4pi
 pref = 0.04787801/cell%getVolume()/(4.0*cPi) 

! preg is used to go from V to U, remembering that gamma is already
! included in the output from fscatt
 pre = pref * preg

! initialize the real and imaginary parts of the structure factor
 ff=czero
 gg=czero

! loop over all atoms in the asymmetric unit
 do m=1,cell%getNatomtype()
  sf = sfarray(m)

! loop over all atoms in the orbit
  p1 = czero
  do j=1,numat(m)
   arg=twopi*sum(float(hkl(1:3))*atpos(m,j,1:3))
   p1 = p1 + exp(cmplx(0.0,-arg))
  end do

  ff = ff + p1*real(sf)
  gg = gg + p1*aimag(sf)

 end do
!
! fill in the entries of the rlp variable
 self%rlp%hkl = hkl

! these are the modulus and phase of the real part of Vg
 self%rlp%Vmod = pref * cabs(ff)
 self%rlp%Vphase = atan2(aimag(ff),real(ff))

! modulus of U_g
 self%rlp%Umod = preg*self%rlp%Vmod

! if absorption is included, also compute the imaginary part of Vg, i.e., Vprime_g
 if (self%rlp%absorption.eqv..TRUE.) then 
  self%rlp%Vpmod = pref * cabs(gg)
  self%rlp%Vpphase = atan2(aimag(gg),real(gg))

! modulus of Uprime_g
  self%rlp%Upmod = preg*self%rlp%Vpmod

! complex Ucg = U_g + i Uprime_g = U_g,r-Uprime_g,i + i(U_g,i+Uprime_g,r)
  self%rlp%Ucg = pre * cmplx(real(ff)-aimag(gg),aimag(ff)+real(gg))
 else ! set absorption parameters to zero
  self%rlp%Vpmod = 0.0
  self%rlp%Vpphase = 0.0
! Ucg = U_g (complex number)
  self%rlp%Ucg = pre * ff
 end if

! complex Vg 
 self%rlp%Vg = self%rlp%Ucg/preg
 if (abs(self%rlp%Umod).gt.0.0) then 
  self%rlp%xg = 1.0/abs(self%rlp%Umod)/self%Lambda
 else
  self%rlp%xg = 1.0E+8
 end if 

 if (abs(self%rlp%Upmod).gt.0.0) then 
  self%rlp%xgp = 1.0/abs(self%rlp%Upmod)/self%Lambda
 else
  self%rlp%xgp = 1.0E+8
 end if 

 if (self%rlp%absorption.eqv..TRUE.) then 
  self%rlp%ar = self%rlp%xgp/self%rlp%xg
  if (present(applyqgshift)) then
    if (applyqgshift.eqv..TRUE.) then
      self%rlp%qg = cmplx(cos(self%rlp%Vphase)/self%rlp%xg-sin(self%rlp%Vpphase)/self%rlp%xgp, &
                          cos(self%rlp%Vpphase)/self%rlp%xgp+sin(self%rlp%Vphase)/self%rlp%xg)
    end if
  else
    arg = self%rlp%Vpphase-self%rlp%Vphase
    self%rlp%qg = cmplx(1.0/self%rlp%xg-sin(arg)/self%rlp%xgp,cos(arg)/self%rlp%xgp)
  end if
 else
  self%rlp%ar = 0.0
  self%rlp%qg = cmplx(1.0/self%rlp%xg,0.0)
 end if

 deallocate(sfarray)
end if

end subroutine CalcUcg_

!--------------------------------------------------------------------------
recursive subroutine PreCalcFSCATT_(self, cell, dmin, gstep)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/27/20
  !!
  !! precompute the FSCATT values for interpolation purposes, to speed up STEM-DCI and other codes

use mod_crystallography
use mod_others

IMPLICIT NONE

class(Diffraction_T),INTENT(INOUT)  :: self 
type(Cell_T),INTENT(INOUT)          :: cell
real(kind=sgl),INTENT(IN)           :: dmin  
 !! smallest d-spacing to consider
real(kind=sgl),INTENT(IN)           :: gstep 
 !! step size in the cell%scatfacg array
    
integer(kind=irg)                   :: j,m,ii,i
real(kind=sgl)                      :: s,ul
real(kind=sgl),parameter            :: swk = 0.628318530717959
real(kind=sgl),parameter            :: dwwk = 1.26651479552922
integer(kind=irg),parameter         :: absflg = 3
logical                             :: accflg=.TRUE., dwflg=.TRUE.
character(2)                        :: smb
real(kind=dbl), allocatable         :: apos(:,:)
integer(kind=irg), allocatable      :: atp(:)

! first generate the array of s-values for which the scattering factors need to be computed
s = 2.0/dmin   ! maximum range in reciprocal space
self%numscatfac = nint(s/gstep) + 2
allocate(self%scatfacg(self%numscatfac))
self%scatfacg = (/ (gstep * float(i-1),i=1,self%numscatfac) /)
self%scatfacg = self%scatfacg * swk

! allocate the scattering factor interpolation array
allocate( self%scatfac(self%numscatfac, cell%getNatomtype()) )

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
  
apos = cell%getAsymPosData() 
atp = cell%getatomtype() 

do i=1,self%numscatfac
! properly scale the scattering parameter
 s = self%scatfacg(i)

! loop over all atoms in the asymmetric unit
 do m=1,cell%getNatomtype()
! get the atomic scattering factor for this atom
! scale and include Debye-Waller factor and site occupation parameter
  ul = sqrt(apos(m,5)*dwwk)
  j = atp(m)
  self%scatfac(i,m) = FSCATT(s,ul,j,smb,sngl(self%voltage),absflg,accflg,dwflg)*apos(m,4)
 end do 
end do 

end subroutine PreCalcFSCATT_

!--------------------------------------------------------------------------
recursive subroutine getScatfac_(self, cell, s, sfarray, ntypes)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/27/20
  !!
  !! interpolate the precomputed FSCATT values

use mod_crystallography
use mod_others

IMPLICIT NONE

class(Diffraction_T),INTENT(INOUT)  :: self 
type(Cell_T),INTENT(INOUT)          :: cell
real(kind=sgl),INTENT(IN)           :: s       
 !! reciprocal distance value
integer(kind=irg),INTENT(IN)        :: ntypes  
 !! number of types 
complex(kind=sgl),INTENT(OUT)       :: sfarray(ntypes) 
 !! returned scattering factor values
    
integer(kind=irg)                   :: jj
real(kind=sgl)                      :: dx

if (s.eq.0.0) then 
    sfarray(1:ntypes) = self%scatfac(1,1:ntypes)
else
    jj = ifix(s/self%scatfacg(2))
    if (jj.ge.self%numscatfac) then
        sfarray(1:ntypes) = self%scatfac(self%numscatfac,1:ntypes)
    else
        dx = s/self%scatfacg(2) - float(jj)
        sfarray(1:ntypes) = self%scatfac(jj,1:ntypes)*(1.0-dx) + &
                                     self%scatfac(jj+1,1:ntypes)*dx
    end if
end if

end subroutine getScatfac_

!--------------------------------------------------------------------------
recursive function getWaveLength_(self) result(Lambda)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/28/20
  !!
  !! returns the electron wavelength 

IMPLICIT NONE 

class(Diffraction_T),INTENT(INOUT)  :: self 
real(kind=dbl)                      :: Lambda 

Lambda = self%Lambda 

end function getWaveLength_

!--------------------------------------------------------------------------
recursive function getBetheParameter_(self, c) result(bp)
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/02/20
  !!
  !! returns one of the Bethe parameters

use mod_io

IMPLICIT NONE 

class(Diffraction_T),INTENT(INOUT)  :: self 
character(2),INTENT(IN)             :: c
real(kind=dbl)                      :: bp

type(IO_T)                          :: Message

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
recursive subroutine Set_Bethe_Parameters_(self, EMsoft, silent, filename)
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/02/20
  !!
  !! Read the Bethe potential parameters from a file, if it exists; otherwise take defaults
  !!
  !! The parameters set in this routine determine the difference between strong and
  !! weak reflections.  The routine checks for the presence of the BetheParameters.nml file
  !! in the current folder.  If present, it will read the parameters, otherwise it will use 
  !! defaults which have been determined to be reasonable based on dynamical EBSD runs. 

use mod_EMsoft
use mod_io

IMPLICIT NONE

class(Diffraction_T), INTENT(INOUT)         :: self
type(EMsoft_T), INTENT(INOUT)               :: EMsoft
logical,INTENT(IN),OPTIONAL                 :: silent
character(fnlen),INTENT(IN),OPTIONAL        :: filename

type(IO_T)                                  :: Message
character(fnlen)                            :: Bethefilename, fname
logical                                     :: fexist
real(kind=sgl)                              :: c1, c2, c3, sgdbdiff

namelist /BetheList/ c1, c2, c3, sgdbdiff

if (present(filename)) then
  Bethefilename = trim(filename)
else 
  Bethefilename = 'BetheParameters.nml'
end if 

! check for the presence of the namelist file in the current folder
inquire(file=trim(Bethefilename),exist=fexist)

! set all default values (must be done here, since nml file may not contain all of them)
c1 = 4.0_sgl           ! changed from 8 and 12 for a test on 8/14/15
c2 = 8.0_sgl           !
c3 = 50.0_sgl          !
sgdbdiff = 1.00_sgl    !

if (fexist) then ! check for the file in the local folder
! read the parameters from the file
 fname = EMsoft%toNativePath(Bethefilename)
 open(UNIT=dataunit,FILE=trim(fname),DELIM='APOSTROPHE')
 READ(UNIT=dataunit,NML=BetheList)
 close(UNIT=dataunit)
 if (.not.present(silent)) then
   call Message%printMessage('Read Bethe parameters from BetheParameters.nml', frm = "(A)")
   write (6,nml=BetheList)
 end if
end if

self%BetheParameters%c1 = c1
self%BetheParameters%c2 = c2
self%BetheParameters%c3 = c3
self%BetheParameters%sgdbdiff = sgdbdiff

end subroutine Set_Bethe_Parameters_

!--------------------------------------------------------------------------
recursive function CalcsgSingle_(self, cell, gg, kk, FN) result(sg)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/27/20
  !!
  !! compute the excitation error for a given reflection (single precision)

use mod_crystallography

IMPLICIT NONE

class(Diffraction_T),INTENT(INOUT)  :: self 
type(Cell_T),INTENT(INOUT)          :: cell
real(kind=sgl),INTENT(IN)           :: gg(3)                
 !! reciprocal lattice point
real(kind=sgl),INTENT(IN)           :: kk(3)                
 !! wave vector
real(kind=sgl),INTENT(IN)           :: FN(3)                
 !! foil normal
    
real(kind=sgl)                      :: kpg(3),tkpg(3),xnom,xden,q1,q2,sg

 kpg=kk+gg
 tkpg=2.0*kk+gg

! use equation of Ewald sphere
 xnom = -cell%CalcDot(gg,tkpg,'r')

! 2|k0+g|cos(alpha) = 2(k0+g).Foilnormal
 q1 = cell%CalcLength(kpg,'r')
 q2 = cell%CalcAngle(kpg,FN,'r')
 xden = 2.0*q1*cos(q2)
 sg = xnom/xden

end function CalcsgSingle_

!--------------------------------------------------------------------------
recursive function CalcsgDouble_(self, cell, gg, kk, FN) result(sg)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/27/20
  !!
  !! compute the excitation error for a given reflection (double precision)

use mod_crystallography

IMPLICIT NONE

class(Diffraction_T),INTENT(INOUT)  :: self 
type(Cell_T),INTENT(INOUT)          :: cell
real(kind=dbl),INTENT(IN)           :: gg(3)                
 !! reciprocal lattice point
real(kind=dbl),INTENT(IN)           :: kk(3)                
 !! wave vector
real(kind=dbl),INTENT(IN)           :: FN(3)                
 !! foil normal
    
real(kind=dbl)                      :: kpg(3),tkpg(3),xnom,xden,q1,q2,sg

 kpg=kk+gg
 tkpg=2.D0*kk+gg

! use equation of Ewald sphere
 xnom = -cell%CalcDot(gg,tkpg,'r')

! 2|k0+g|cos(alpha) = 2(k0+g).Foilnormal
 q1 = cell%CalcLength(kpg,'r')
 q2 = cell%CalcAngle(kpg,FN,'r')
 xden = 2.D0*q1*dcos(q2)
 sg = xnom/xden

end function CalcsgDouble_

!--------------------------------------------------------------------------
recursive subroutine BWsolve_(self, M, W, CGG, CGinv, nn, IPIV)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/27/20
  !!
  !! Solve the Bloch wave eigenvalue equation for the 
  !! N-beam case, using the ZGEEV LAPACK 3.0 routine

use mod_io

IMPLICIT NONE

class(Diffraction_T),INTENT(INOUT)  :: self 
integer(kind=irg),INTENT(IN)        :: nn           
 ! number of beams
complex(kind=dbl),INTENT(IN)        :: M(nn,nn)     
 ! input dynamical matrix
complex(kind=dbl),INTENT(OUT)       :: W(nn)        
 ! Bloch eigenvalues
complex(kind=dbl),INTENT(OUT)       :: CGG(nn,nn)   
 ! Bloch eigenvectors
complex(kind=dbl),INTENT(OUT)       :: CGinv(nn,nn) 
 ! inverse of eigenvector array
integer(kind=irg),INTENT(IN)        :: IPIV(nn)     
 ! pivot array, currently unused
    
type(IO_T)                          :: Message
integer(kind=irg)                   :: INFO, LDA, LDVR, LDVL, LWORK, JPIV(nn),MILWORK, i, io_int(1)
integer(kind=irg),parameter         :: LWMAX = 5000 
complex(kind=dbl)                   :: VL(nn,nn),  WORK(LWMAX), normsum
real(kind=dbl)                      :: RWORK(2*nn), io_real(1)
character                           :: JOBVL, JOBVR
complex(kind=dbl),allocatable       :: MIWORK(:)

! set some initial LAPACK variables 
 LDA = nn
 LDVL = nn
 LDVR = nn
 INFO = 0
 
! first initialize the parameters for the LAPACK ZGEEV, CGETRF, and CGETRI routines
 JOBVL = 'N'   ! do not compute the left eigenvectors
 JOBVR = 'V'   ! do compute the right eigenvectors
 LWORK = -1    ! so that we can ask the routine for the actually needed value

! call the routine to determine the optimal workspace size
  call zgeev(JOBVL,JOBVR,nn,M,LDA,W,VL,LDVL,CGG,LDVR,WORK,LWORK,RWORK,INFO)
  LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

! then call the eigenvalue solver
  call zgeev(JOBVL,JOBVR,nn,M,LDA,W,VL,LDVL,CGG,LDVR,WORK,LWORK,RWORK,INFO)
  if (INFO.ne.0) call Message%printError('Error in BWsolve: ','ZGEEV return not zero')

! it appears that the eigenvectors may not always be normalized ...
! so we renormalize them here...
! do i=1,nn
!   normsum = sum(abs(CGG(1:nn,i))**2)
!   normsum = cmplx(1.0,0.0,dbl)/sqrt(normsum)
!   CGG(1:nn,i) = CGG(1:nn,i)*normsum
! end do

! make a copy of CG for the matrix inversion routines
 CGinv = CGG

! invert CGinv to get the Bloch wave excitation amplitudes 
 LDA = nn
 call zgetrf(nn,nn,CGinv,LDA,JPIV,INFO)
 if (INFO.ne.0) then
  io_int(1) = INFO
  call Message%WriteValue('zgetrf error code: ',io_int,1,frm="(I5)")
  call Message%printError('Error in BWsolve: ','ZGETRF return not zero')
 end if

 MILWORK = 64*nn 
 allocate(MIWORK(MILWORK))

 MIWORK = cmplx(0.0_dbl,0.0_dbl)
 call zgetri(nn,CGinv,LDA,JPIV,MIWORK,MILWORK,INFO)
 if (INFO.ne.0) then
  io_int(1) = INFO
  call Message%WriteValue('zgetri error code: ',io_int,1,frm="(I5)")
  call Message%printError('Error in BWsolve: ','ZGETRI return not zero')
 end if

! if ((abs(sum(matmul(CGG,CGinv)))-dble(nn)).gt.1.E-8) then
!  call Message('Error in matrix inversion; continuing', frm = "(A)")
!  io_real(1) = abs(sum(matmul(CGG,CGinv)))-dble(nn)
!  call WriteValue('   Matrix inversion error; this number should be zero: ',io_real,1,"(F)")
! endif
  
 deallocate(MIWORK)
 
end subroutine BWsolve_

! !
! ! ###################################################################
! ! 
! !  subroutine CalcFresnelPropagator
! !
! !                                    created: 4/16/97
! !  Author: Marc De Graef
! !  
! !  Description: compute the Fresnel propagator (for possibly inclined
! !               illumination) and store it in a file
! ! 
! !  History
! ! 
! !  modified by  rev reason
! !  -------- --- --- -----------
! !   4/16/97 MDG 1.0 original
! !   9/29/01 MDG 2.0 f90
! !  11/27/01 MDG 2.1 added kind support
! ! ###################################################################
! recursive subroutine CalcFresnelPropagator(beam,dimi,dimj,dz,scl,propname,lambda)
! !DEC$ ATTRIBUTES DLLEXPORT :: CalcFresnelPropagator

! use constants
! use io
! use files

! IMPLICIT NONE

! real(kind=sgl)                  :: beam(3),b,bm(2),dz,fidim,fjdim,prefac,scl, oi_real(2), lambda
! real(kind=sgl),allocatable      :: idimi(:),jdimj(:)
! complex(kind=sgl),allocatable   :: fr(:,:)
! integer(kind=irg)               :: dimi,dimj,i,ix,iy
! character(fnlen)                   :: propname

! INTENT(IN) :: beam,dimi,dimj,dz

!   fidim = 1.0/float(dimi)
!   fjdim = 1.0/float(dimj)
!   prefac = scl*cPi*lambda*dz
!   call Message('Computing Fresnel propagator', frm = "(A)")
! ! normalize the incident beam direction and rescale to the wavevector
!   b = sqrt(sum(beam**2))
!   bm= beam(1:2)/b/lambda
!   oi_real(1:2) = bm(1:2) 
!   call WriteValue(' Laue center at ', oi_real, 2, "(F8.4,',',F8.4)")
! ! allocate variables
!   allocate(fr(dimi,dimj))
!   allocate(idimi(dimi),jdimj(dimj))
!   idimi = float((/ (i, i=0,dimi-1) /))
!   jdimj = float((/ (i, i=0,dimj-1) /))
! !
!   where(idimi.ge.dimi/2) idimi = idimi-float(dimi)
!   where(jdimj.ge.dimj/2) jdimj = jdimj-float(dimj)
! !
!   idimi = prefac*idimi*fidim
!   jdimj = prefac*jdimj*fjdim
! !
!   idimi = idimi*(idimi + 2.0*bm(1))
!   jdimj = jdimj*(jdimj + 2.0*bm(2))
! ! loop over y axis  
!   do iy=1,dimj
! ! loop over x-axis
!     do ix=1,dimi
!       fr(ix,iy)=cmplx(cos(idimi(ix)+jdimj(iy)),sin(idimi(ix)+jdimj(iy)))
!     end do
!   end do
!   deallocate(idimi,jdimj)
! ! and store it in a file
!   open(unit=dataunit,file=trim(EMsoft_toNativePath(propname)),form='unformatted')
!   write (dataunit) dimi,dimj
!   write (dataunit) fr
!   close(unit=dataunit,status='keep')
!   deallocate(fr)
! end subroutine


! !--------------------------------------------------------------------------
! !
! ! FUNCTION: LorentzPF
! !
! !> @author Marc De Graef, Carnegie Mellon University
! !
! !> @brief compute the Lorentz Polarization Factor Lp
! !
! !> @param theta scattering angle  
! !> @param HEDM optional string to indicate HEDM mode
! !
! !> @date  03/26/13  MDG  1.0 added for HEDM project
! !--------------------------------------------------------------------------
! recursive function LorentzPF(theta,HEDM) result(tt)
! !DEC$ ATTRIBUTES DLLEXPORT :: LorentzPF

! use crystal

! IMPLICIT NONE

! real(kind=sgl),INTENT(IN)                       :: theta                !< scattering angle
! character(*),INTENT(IN),OPTIONAL                :: HEDM         !< for HEDM we have a different polarization factor
! real(kind=sgl)                                  :: tt

! if (present(HEDM)) then
!   tt = (1.0+cos(2.0*theta)**2) / sin(theta)**2 / cos(theta)
! else
!   tt = (1.0+cos(2.0*theta)**2) / sin(theta)**2 / cos(theta)
! end if

! end function



end module mod_diffraction
