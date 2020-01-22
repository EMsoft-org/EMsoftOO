! ###################################################################
! Copyright (c) 2014-2020 Marc De Graef/Carnegie Mellon University
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

module mod_global
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/30/19
  !!
  !! global constant definitions (except for strings which are in stringconstants.in.f90)

use mod_kinds
use,intrinsic :: ISO_C_BINDING

IMPLICIT NONE
public

!> @note This module must be "use"d by every program, subroutine, and function.
!> These are the only global variables used by the EMsoft package.

!> character type used for json routines
  integer,parameter                     :: jsonCK = selected_char_kind('DEFAULT')
!DEC$ ATTRIBUTES DLLEXPORT :: jsonCK

!> standard string length for filenames
  integer(kind=irg),parameter           :: fnlen = 512
!DEC$ ATTRIBUTES DLLEXPORT :: fnlen

!> counter for non-fatal error messages in handling of environment variables;
!> if this counter equals 0, then error warnings will be printed, otherwise not.
  integer(kind=irg)                     :: displayEMsoftWarningMessages = 0
  integer(kind=irg)                     :: displayConfigFileMissingMessage = 0
!DEC$ ATTRIBUTES DLLEXPORT :: displayEMsoftWarningMessages
!DEC$ ATTRIBUTES DLLEXPORT :: displayConfigFileMissingMessage

!> standard array size for all wrapper routine calls; applies to ipar, fpar, and spar arrays
  integer(c_int32_t),parameter          :: wraparraysize = 80
!DEC$ ATTRIBUTES DLLEXPORT :: wraparraysize

!> reserved IO unit identifiers for data (21-23)
  integer(kind=irg), parameter          :: dataunit = 21, dataunit2 = 22, dataunit3 = 23
!DEC$ ATTRIBUTES DLLEXPORT :: dataunit
!DEC$ ATTRIBUTES DLLEXPORT :: dataunit2
!DEC$ ATTRIBUTES DLLEXPORT :: dataunit3

!> Maximum number of positions in asymmetric unit
  integer(kind=irg), parameter          :: maxpasym = 250   
!DEC$ ATTRIBUTES DLLEXPORT :: maxpasym

! ****************************************************
! ****************************************************
! ****************************************************
! these two parameters are used to change the sign of the permutation symbol from Adam 
! Morawiec's book to the convention used for the EMsoft package.  If you want to use Adam's 
! convention, both of these parameters should be set to +1; -1 will change the sign everywhere
! for all representations that involve the unit vector.  The quaternion product is 
! also redefined to include the epsijk parameter.  Doing so guarantees that the 
! quat_Lp operator ALWAYS returns an active result, regardless of the choice of epsijk;
! quat_LPstar ALWAYS returns a passive result.

! Uncomment these for an alternative way of doing things
!real(kind=sgl), parameter :: epsijk = -1.0
!real(kind=dbl), parameter :: epsijkd = -1.D0

! uncomment these for the Morawiec version.
  real(kind=sgl), parameter :: epsijk = 1.0
  real(kind=dbl), parameter :: epsijkd = 1.D0
!DEC$ ATTRIBUTES DLLEXPORT :: epsijk
!DEC$ ATTRIBUTES DLLEXPORT :: epsijkd

! In the first case, epsijk=-1, the rotation 120@[111] will result in 
! an axis angle pair of [111], 2pi/3.  In the second case, the axis-angle 
! pair will be -[111], 2pi/3.  In all cases, the rotations are interpreted
! in the passive sense.  The case epsijk=+1 corresponds to the mathematically 
! consistent case, using the standard definition for the quaternion product; in
! the other case, epsijk=-1, one must redefine the quaternion product in order
! to produce consistent results.  This takes a lengthy explanation ... see the
! rotations tutorial paper for an in-depth explanation.  These changes propagate
! to a number of files, notably quaternions.f90, and everywhere else that quaternions
! and rotations in general are used (in particular the mod_rotations module).
!
! Reference:  D.J. Rowenhorst, A.D. Rollett, G.S. Roher, M.A. Groeber, M.A. Jackson, 
!  P.J. Konijnenberg, and M. De Graef. "Tutorial: consistent representations of and 
!  conversions between 3D rotations". Modeling and Simulations in Materials Science 
!  and Engineering, 23, 083501 (2015).
!
! ****************************************************
! ****************************************************
! ****************************************************


! various physical constants
!> cPi          = pi [dimensionless]
!> cLight       = velocity of light [m/s]
!> cPlanck      = Planck''s constant [Js]
!> cBoltzmann   = Boltmann constant [J/K]
!> cPermea      = permeability of vacuum [4pi 10^7 H/m]
!> cPermit      = permittivity of vacuum [F/m]
!> cCharge      = electron charge [C]
!> cRestmass    = electron rest mass [kg]
!> cMoment      = electron magnetic moment [J/T]
!> cJ2eV        = Joules per eV
!> cAvogadro    = Avogadro's constant [mol^-1]
!
! The values of several of these constants have been updated to the new SI 2019 exact values [MDG, 01/22/19]
! The exact values below are the ones for cLight, cPlanck, cBoltzmann, cCharge; the others are derived using 
! the standard relations in the 2019 SI units document.  In the derivation, we used 0.0072973525664D0 as the 
! value for the hyperfine structure constant alpha. 
!
  real(kind=dbl), parameter :: cPi=3.141592653589793238D0, cLight = 299792458.D0, &
                               cPlanck = 6.62607015D-34, cBoltzmann = 1.380649D-23,  &
                               cPermea = 1.2566370616D-6, cPermit = 8.8541878163D-12, &
                               cCharge = 1.602176634D-19, cRestmass = 9.1093837090D-31, &
                               cMoment = 9.2740100707D-24, cJ2eV = 1.602176565D-19, &
                               cAvogadro = 6.02214076D23
!DEC$ ATTRIBUTES DLLEXPORT :: cPi
!DEC$ ATTRIBUTES DLLEXPORT :: cPlanck
!DEC$ ATTRIBUTES DLLEXPORT :: cPermea
!DEC$ ATTRIBUTES DLLEXPORT :: cCharge
!DEC$ ATTRIBUTES DLLEXPORT :: cMoment
!DEC$ ATTRIBUTES DLLEXPORT :: cAvogadro
!DEC$ ATTRIBUTES DLLEXPORT :: cLight
!DEC$ ATTRIBUTES DLLEXPORT :: cBoltzmann
!DEC$ ATTRIBUTES DLLEXPORT :: cPermit
!DEC$ ATTRIBUTES DLLEXPORT :: cRestmass
!DEC$ ATTRIBUTES DLLEXPORT :: cJ2eV

  real(kind=dbl), parameter :: dtor = cPi/180.D0 
  real(kind=dbl), parameter :: rtod = 180.D0/cPi
!DEC$ ATTRIBUTES DLLEXPORT :: dtor
!DEC$ ATTRIBUTES DLLEXPORT :: rtod

!> element symbols (we'll do 1-98 for all parameter lists)
  character(2), dimension(98) :: ATOM_sym=(/' H','He','Li','Be',' B',' C',' N',' O',' F','Ne', &
                                            'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca', &
                                            'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
                                            'Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr', &
                                            'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
                                            'Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd', &
                                            'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
                                            'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg', &
                                            'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
                                            'Pa',' U','Np','Pu','Am','Cm','Bk','Cf'/)
!DEC$ ATTRIBUTES DLLEXPORT :: ATOM_sym

!> atomic weights for things like density computations (from NIST elemental data base)
  real(kind=sgl),dimension(98)    :: ATOM_weights(98) = (/1.00794, 4.002602, 6.941, 9.012182, 10.811, &
                                                          12.0107, 14.0067, 15.9994, 18.9984032, 20.1797, &
                                                          22.98976928, 24.3050, 26.9815386, 28.0855, 30.973762, &
                                                          32.065, 35.453, 39.948, 39.0983, 40.078, &
                                                          44.955912, 47.867, 50.9415, 51.9961, 54.938045, &
                                                          55.845, 58.933195, 58.6934, 63.546, 65.38, &
                                                          69.723, 72.64, 74.92160, 78.96, 79.904, &
                                                          83.798, 85.4678, 87.62, 88.90585, 91.224, &
                                                          92.90638, 95.96, 98.9062, 101.07, 102.90550, &
                                                          106.42, 107.8682, 112.411, 114.818, 118.710, &
                                                          121.760, 127.60, 126.90447, 131.293, 132.9054519, &
                                                          137.327, 138.90547, 140.116, 140.90765, 144.242, &
                                                          145.0, 150.36, 151.964, 157.25, 158.92535, &
                                                          162.500, 164.93032, 167.259, 168.93421, 173.054, &
                                                          174.9668, 178.49, 180.94788, 183.84, 186.207, &
                                                          190.23, 192.217, 195.084, 196.966569, 200.59, &
                                                          204.3833, 207.2, 208.98040, 209.0, 210.0, &
                                                          222.0, 223.0, 226.0, 227.0, 232.03806, &
                                                          231.03588, 238.02891, 237.0, 244.0, 243.0, &
                                                          247.0, 251.0, 252.0 /)
!DEC$ ATTRIBUTES DLLEXPORT :: ATOM_weights


! these are a bunch of constants used for Lambert and related projections; they are all in double precision
  type LambertParametersType
          real(kind=dbl)          :: Pi=3.141592653589793D0       !  pi
          real(kind=dbl)          :: iPi=0.318309886183791D0      !  1/pi
          real(kind=dbl)          :: sPi=1.772453850905516D0      !  sqrt(pi)
          real(kind=dbl)          :: sPio2=1.253314137315500D0    !  sqrt(pi/2)
          real(kind=dbl)          :: sPi2=0.886226925452758D0     !  sqrt(pi)/2
          real(kind=dbl)          :: srt=0.86602540378D0      !  sqrt(3)/2
          real(kind=dbl)          :: isrt=0.57735026919D0    !  1/sqrt(3)
          real(kind=dbl)          :: alpha=1.346773687088598D0   !  sqrt(pi)/3^(1/4)
          real(kind=dbl)          :: rtt=1.732050807568877D0      !  sqrt(3)
          real(kind=dbl)          :: prea=0.525037567904332D0    !  3^(1/4)/sqrt(2pi)
          real(kind=dbl)          :: preb=1.050075135808664D0     !  3^(1/4)sqrt(2/pi)
          real(kind=dbl)          :: prec=0.906899682117109D0    !  pi/2sqrt(3)
          real(kind=dbl)          :: pred=2.094395102393195D0     !  2pi/3
          real(kind=dbl)          :: pree=0.759835685651593D0     !  3^(-1/4)
          real(kind=dbl)          :: pref=1.381976597885342D0     !  sqrt(6/pi)
          real(kind=dbl)          :: preg=1.5551203015562141D0    ! 2sqrt(pi)/3^(3/4)
! the following constants are used for the cube to quaternion hemisphere mapping
          real(kind=dbl)          :: a=1.925749019958253D0        ! pi^(5/6)/6^(1/6)
          real(kind=dbl)          :: ap=2.145029397111025D0       ! pi^(2/3)
          real(kind=dbl)          :: sc=0.897772786961286D0       ! a/ap
          real(kind=dbl)          :: beta=0.962874509979126D0     ! pi^(5/6)/6^(1/6)/2
          real(kind=dbl)          :: R1=1.330670039491469D0       ! (3pi/4)^(1/3)
          real(kind=dbl)          :: r2=1.414213562373095D0       ! sqrt(2)
          real(kind=dbl)          :: r22=0.707106781186547D0      ! 1/sqrt(2)
          real(kind=dbl)          :: pi12=0.261799387799149D0     ! pi/12
          real(kind=dbl)          :: pi8=0.392699081698724D0      ! pi/8
          real(kind=dbl)          :: prek=1.643456402972504D0     ! R1 2^(1/4)/beta
          real(kind=dbl)          :: r24=4.898979485566356D0      ! sqrt(24)

          real(kind=dbl)          :: tfit(21) = (/ 0.9999999999999968D0, -0.49999999999986866D0,  &
                                                -0.025000000000632055D0, - 0.003928571496460683D0, &
                                                -0.0008164666077062752D0, - 0.00019411896443261646D0, &
                                                -0.00004985822229871769D0, - 0.000014164962366386031D0, &
                                                -1.9000248160936107D-6, - 5.72184549898506D-6, &
                                                7.772149920658778D-6, - 0.00001053483452909705D0, &
                                                9.528014229335313D-6, - 5.660288876265125D-6, &
                                                1.2844901692764126D-6, 1.1255185726258763D-6, &
                                                -1.3834391419956455D-6, 7.513691751164847D-7, &
                                                -2.401996891720091D-7, 4.386887017466388D-8, &
                                                -3.5917775353564864D-9 /)

          real(kind=dbl)          :: BP(12)= (/ 0.D0, 1.D0, 0.577350269189626D0, 0.414213562373095D0, 0.D0,  &
                                             0.267949192431123D0, 0.D0, 0.198912367379658D0, 0.D0, &
                                             0.158384440324536D0, 0.D0, 0.131652497587396D0/)       ! used for Fundamental Zone determination in so3 module
end type LambertParametersType

type(LambertParametersType)        :: LPs
!DEC$ ATTRIBUTES DLLEXPORT :: LPs


  real(kind=sgl),dimension(81)        :: Butterfly9x9 = (/-10.0, -15.0, -22.0, -22.0, -22.0, -22.0, -22.0, -15.0, -10.0, &
                                                         -1.0, -6.0, -13.0, -22.0, -22.0, -22.0, -13.0, -6.0, -1.0, &
                                                          3.0, 6.0, 4.0, -3.0, -22.0, -3.0, 4.0, 6.0, 3.0, & 
                                                          3.0, 11.0, 19.0, 28.0, 42.0, 28.0, 19.0, 11.0, 3.0, &
                                                          3.0, 11.0, 27.0, 42.0, 42.0, 42.0, 27.0, 11.0, 3.0, &
                                                          3.0, 11.0, 19.0, 28.0, 42.0, 28.0, 19.0, 11.0, 3.0, &
                                                          3.0, 6.0, 4.0, -3.0, -22.0, -3.0, 4.0, 6.0, 3.0, & 
                                                         -1.0, -6.0, -13.0, -22.0, -22.0, -22.0, -13.0, -6.0, -1.0, &
                                                         -10.0, -15.0, -22.0, -22.0, -22.0, -22.0, -22.0, -15.0, -10.0/)
!DEC$ ATTRIBUTES DLLEXPORT :: Butterfly9x9

! vertex coordinates of the icosahedron (normalized)
  real(kind=dbl),dimension(3,12)      :: IcoVertices = reshape( (/ 0D0,0.D0,1.D0, &
                                       0.89442719099991587856D0,0.D0,0.44721359549995793928D0, &
                                       0.27639320225002103036D0,0.85065080835203993218D0,0.44721359549995793928D0, &
                                      -0.72360679774997896964D0,0.52573111211913360603D0,0.44721359549995793928D0, &
                                      -0.72360679774997896964D0,-0.52573111211913360603D0,0.44721359549995793928D0, &
                                       0.27639320225002103036D0,-0.85065080835203993218D0,0.44721359549995793928D0, &
                                      -0.89442719099991587856D0,0.D0,-0.44721359549995793928D0, &
                                      -0.27639320225002103036D0,-0.85065080835203993218D0,-0.44721359549995793928D0, &
                                       0.72360679774997896964D0,-0.52573111211913360603D0,-0.44721359549995793928D0, &
                                       0.72360679774997896964D0,0.52573111211913360603D0,-0.44721359549995793928D0, &
                                      -0.27639320225002103036D0,0.85065080835203993218D0,-0.44721359549995793928D0, &
                                       0.D0,0.D0,-1.D0 /), (/3,12/))
!DEC$ ATTRIBUTES DLLEXPORT :: IcoVertices

  character(30),dimension(27) :: ConfigStructureNames = (/ "EMsoftpathname                ", &
                                                           "EMXtalFolderpathname          ", &
                                                           "EMdatapathname                ", &
                                                           "EMtmppathname                 ", &
                                                           "EMsoftLibraryLocation         ", &
                                                           "EMSlackWebHookURL             ", &
                                                           "EMSlackChannel                ", &
                                                           "UserName                      ", &
                                                           "UserLocation                  ", &
                                                           "UserEmail                     ", &
                                                           "EMNotify                      ", &
                                                           "Develop                       ", &
                                                           "Release                       ", &
                                                           "h5copypath                    ", &
                                                           "EMsoftplatform                ", &
                                                           "EMsofttestpath                ", &
                                                           "EMsoftTestingPath             ", &
                                                           "EMsoftversion                 ", &
                                                           "Configpath                    ", &
                                                           "Templatepathname              ", &
                                                           "Resourcepathname              ", &
                                                           "Homepathname                  ", &
                                                           "OpenCLpathname                ", &
                                                           "Templatecodefilename          ", &
                                                           "WyckoffPositionsfilename      ", &
                                                           "Randomseedfilename            ", &
                                                           "EMsoftnativedelimiter         " /)
!DEC$ ATTRIBUTES DLLEXPORT :: ConfigStructureNames



end module mod_global
