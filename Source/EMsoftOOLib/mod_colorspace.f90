! ###################################################################
! Copyright (c) 2013-2022, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_colorspace
  !! author: MDG 
  !! version: 1.0 
  !! date: 09/07/21
  !!
  !! class definition for colorspace conversions
  !!
  !! This module is a direct fortran translation of Will Lenthe's colorspace.hpp
  !! <https:! github.com/wlenthe/UniformBicone>
  !!
  !! rgb, hsv, and hsl are restricted to the range (0,1) with values outside representing imaginary colors
  !! the range (0,1) is used for hue in hsv/hsl (not (0,2*pi) or (0,360))
  !! all conversions are available using the shortest path in network below
  !!        .->hsv <---> rgb <---> xyz <---> luv
  !!        `->hsl <-'                   `-> lab
  !!        therefore the following direct transformations are implemented
  !!         -rgb2xyz, rgb2hsv, rgb2hsl
  !!         -xyz2rgb, xyz2luv, xyz2lab
  !!         -luv2xyz, lab2xyz
  !!         -hsv2rgb, hsv2hsl
  !!         -hsl2rgb, hsl2hsv
  !!        with indirect transformations requiring a combination of the above
  !! abc: values to convert from abc color space
  !! ijk: location to write values converted to ijk color space (can be the same as parameter abc)
  !! ill: standard illuminant as xyz (only required for conversions involving xyz<->luv or xyz<->lab, defaults to CIE illuminant D65 for 2 degree observer)
  !! : true/false if the values fall outside the ijk gamut for conversions to that pass through xyz2rgb (void for others)
  !!

use mod_kinds
use mod_global

IMPLICIT NONE 

! six color representations
type colorspacetype 
  real(kind=dbl) :: rgb(3) ! {r , g , b }: standard red, green, blue (sRGB)
  real(kind=dbl) :: xyz(3) ! {X , Y , Z }: CIE 1931 XYZ color space ('master' (original) perceptually uniform color space)
  real(kind=dbl) :: luv(3) ! {L*, u*, v*}: 1976 CIELUV color space (perceptually uniform space for computer displays)
  real(kind=dbl) :: lab(3) ! {L*, a*, b*}: 1976 CIELab color space (perceptually uniform space for print)
  real(kind=dbl) :: hsv(3) ! {h , s , v }: hue, saturation, value (cylindrical)
  real(kind=dbl) :: hsl(3) ! {h , s , l }: hue, saturation, lightness (cylindrical)
end type colorspacetype

! spline control points for 4- and 6-fold perceptually uniform color spaces
real(kind=dbl)   :: P4(3,15) = reshape( (/ 50.000000D0,  20.260871D0, -61.121717D0, &
                                           40.000000D0,  55.743174D0, -72.075429D0, &
                                           50.000000D0,  66.074726D0, -25.689258D0, &
                                           60.000000D0,  76.406278D0,  20.696912D0, &
                                           70.000000D0,  86.737831D0,  67.083083D0, &
                                           60.000000D0,  45.201449D0,  61.044479D0, &
                                           50.000000D0,  3.6650678D0,  55.005876D0, &
                                           40.000000D0, -37.871314D0,  48.967272D0, &
                                           50.000000D0, -42.148787D0,  19.573417D0, &
                                           60.000000D0, -46.426261D0, -9.8204386D0, &
                                           70.000000D0, -50.703734D0, -39.214294D0, &
                                           60.000000D0, -15.221432D0, -50.168006D0, &
                                           50.000000D0,  20.260871D0, -61.121717D0, &
                                           40.000000D0,  55.743174D0, -72.075429D0, &
                                           50.000000D0,  66.074726D0, -25.689258D0/), (/ 3, 15 /) )

real(kind=dbl)   :: P6(3,21) = reshape( (/ 50.000000D0,  110.55193D0,  -1.8333167D0, &
                                           40.000000d0,  131.49157d0,   28.367970d0, &
                                           50.000000d0,  96.706947d0,   43.686247d0, &
                                           60.000000d0,  61.922323d0,   59.004523d0, &
                                           70.000000d0,  27.137700d0,   74.322800d0, &
                                           60.000000d0,  5.4672700d0,   65.871937d0, &
                                           50.000000d0, -16.203160d0,   57.421073d0, &
                                           40.000000d0, -37.873590d0,   48.970210d0, &
                                           50.000000d0, -43.154747d0,   27.643087d0, &
                                           60.000000d0, -48.435903d0,   6.3159633d0, &
                                           70.000000d0, -53.717060d0,  -15.011160d0, &
                                           60.000000d0, -40.878273d0,  -54.308313d0, &
                                           50.000000d0, -28.039487d0,  -93.605467d0, &
                                           40.000000d0, -15.200700d0,  -132.90262d0, &
                                           50.000000d0,  12.757087d0,  -109.34704d0, &
                                           60.000000d0,  40.714873d0,  -85.791467d0, &
                                           70.000000d0,  68.672660d0,  -62.235890d0, &
                                           60.000000d0,  89.612297d0,  -32.034603d0, &
                                           50.000000d0,  110.55193d0,  -1.8333167d0, &
                                           40.000000d0,  131.49157d0,   28.367970d0, &
                                           50.000000d0,  96.706947d0,   43.686247d0/), (/ 3, 21 /) )

! constants for standard illuminants and common rgb gamuts
type colorstandardstype 
! standard illuminants as xyz (normalized XYZ)
    real(kind=dbl) :: A_2(0:2)    = (/ 0.44757D0, 0.40745D0, 0.14498D0 /)
    real(kind=dbl) :: A_10(0:2)   = (/ 0.45117D0, 0.40594D0, 0.14289D0 /)
    real(kind=dbl) :: B_2(0:2)    = (/ 0.34842D0, 0.35161D0, 0.29997D0 /)
    real(kind=dbl) :: B_10(0:2)   = (/ 0.34980D0, 0.35270D0, 0.29750D0 /)
    real(kind=dbl) :: C_2(0:2)    = (/ 0.31006D0, 0.31616D0, 0.37378D0 /)
    real(kind=dbl) :: C_10(0:2)   = (/ 0.31039D0, 0.31905D0, 0.37056D0 /)
    real(kind=dbl) :: D50_2(0:2)  = (/ 0.34567D0, 0.35850D0, 0.29583D0 /)
    real(kind=dbl) :: D50_10(0:2) = (/ 0.34773D0, 0.35952D0, 0.29275D0 /)
    real(kind=dbl) :: D55_2(0:2)  = (/ 0.33242D0, 0.34743D0, 0.32015D0 /)
    real(kind=dbl) :: D55_10(0:2) = (/ 0.33411D0, 0.34877D0, 0.31712D0 /)
    real(kind=dbl) :: D65_2(0:2)  = (/ 0.31271D0, 0.32902D0, 0.35827D0 /)
    real(kind=dbl) :: D65_10(0:2) = (/ 0.31382D0, 0.33100D0, 0.35518D0 /)
    real(kind=dbl) :: D75_2(0:2)  = (/ 0.29902D0, 0.31485D0, 0.38613D0 /)
    real(kind=dbl) :: D75_10(0:2) = (/ 0.29968D0, 0.31740D0, 0.38292D0 /)
    real(kind=dbl) :: E(0:2)      = (/ 1.D0/3.D0, 1.D0/3.D0, 1.D0/3.D0 /)

! RGB chromaticities as xyz (normalized XYZ)
    real(kind=dbl) :: sRGB(0:2,0:2)     = transpose(reshape( (/ 0.6400D0, 0.3300D0, 0.0300D0, &
                                                                0.3000D0, 0.6000D0, 0.1000D0, &
                                                                0.1500D0, 0.0600D0, 0.7900D0 /), (/ 3,3 /) ) )
    real(kind=dbl) :: cieRGB(0:2,0:2)   = transpose(reshape( (/ 0.7347D0, 0.2653D0, 0.0000D0, &
                                                                0.2738D0, 0.7174D0, 0.0088D0, &
                                                                0.1666D0, 0.0089D0, 0.8245D0 /), (/ 3,3 /) ) )
    real(kind=dbl) :: appleRGB(0:2,0:2) = transpose(reshape( (/ 0.6250D0, 0.3400D0, 0.0350D0, &
                                                                0.2800D0, 0.5950D0, 0.1250D0, &
                                                                0.1550D0, 0.0700D0, 0.7750D0 /), (/ 3,3 /) ) )
    real(kind=dbl) :: adobeRGB(0:2,0:2) = transpose(reshape( (/ 0.6400D0, 0.3300D0, 0.0300D0, &
                                                                0.2100D0, 0.7100D0, 0.0800D0, &
                                                                0.1500D0, 0.0600D0, 0.7900D0 /), (/ 3,3 /) ) )
    real(kind=dbl) :: palRGB(0:2,0:2)   = transpose(reshape( (/ 0.6400D0, 0.3300D0, 0.0300D0, &
                                                                0.2900D0, 0.6000D0, 0.1100D0, &
                                                                0.1500D0, 0.0600D0, 0.7900D0 /), (/ 3,3 /) ) )
    real(kind=dbl) :: ntscRGB(0:2,0:2)  = transpose(reshape( (/ 0.6300D0, 0.3400D0, 0.0300D0, &
                                                                0.3100D0, 0.5950D0, 0.0950D0, &
                                                                0.1550D0, 0.0700D0, 0.7750D0 /), (/ 3,3 /) ) )

! sRGB gamma correction constants
    real(kind=dbl) :: sA     = 0.055D0  
    real(kind=dbl) :: sGamma = 2.4D0    
    real(kind=dbl) :: sPhi   = 12.92D0  
    real(kind=dbl) :: sK0    = 0.04045D0

! sRGB conversion matrices
    real(kind=dbl)           :: sRGBmat(0:2,0:2)  !  = rgbMat(Standards%sRGB, Standards<T>::D65_2)
    real(kind=dbl)           :: sRGBmatInv(0:2,0:2)  !  = inv3x3(rgbMat(Standards%sRGB, Standards<T>::D65_2))

end type colorstandardstype

! class definition
type, public :: colorspace_T
private 
  type(colorstandardstype)    :: Standards
  real(kind=dbl),allocatable  :: P(:,:)
  integer(kind=irg)           :: Nfold 
  integer(kind=irg)           :: N

contains
private 
  procedure, pass(self) :: rgb2xyz_
  procedure, pass(self) :: rgb2luv_
  procedure, pass(self) :: rgb2lab_
  procedure, pass(self) :: rgb2hsv_
  procedure, pass(self) :: rgb2hsl_

  procedure, pass(self) :: xyz2rgb_
  procedure, pass(self) :: xyz2luv_
  procedure, pass(self) :: xyz2lab_
  procedure, pass(self) :: xyz2hsv_
  procedure, pass(self) :: xyz2hsl_
! 
  procedure, pass(self) :: luv2rgb_
  procedure, pass(self) :: luv2xyz_
  procedure, pass(self) :: luv2lab_
  procedure, pass(self) :: luv2hsv_
  procedure, pass(self) :: luv2hsl_
! 
  procedure, pass(self) :: lab2rgb_
  procedure, pass(self) :: lab2xyz_
  procedure, pass(self) :: lab2luv_
  procedure, pass(self) :: lab2hsv_
  procedure, pass(self) :: lab2hsl_

  procedure, pass(self) :: hsv2rgb_
  procedure, pass(self) :: hsv2xyz_
  procedure, pass(self) :: hsv2luv_
  procedure, pass(self) :: hsv2lab_
  procedure, pass(self) :: hsv2hsl_

  procedure, pass(self) :: hsl2rgb_
  procedure, pass(self) :: hsl2xyz_
  procedure, pass(self) :: hsl2luv_
  procedure, pass(self) :: hsl2lab_
  procedure, pass(self) :: hsl2hsv_

  procedure, pass(self) :: sph2rgb_
  procedure, pass(self) :: sphere2rgb_

  generic, public :: rgb2xyz => rgb2xyz_
  generic, public :: rgb2luv => rgb2luv_
  generic, public :: rgb2lab => rgb2lab_
  generic, public :: rgb2hsv => rgb2hsv_
  generic, public :: rgb2hsl => rgb2hsl_

  generic, public :: xyz2rgb => xyz2rgb_
  generic, public :: xyz2luv => xyz2luv_
  generic, public :: xyz2lab => xyz2lab_
  generic, public :: xyz2hsv => xyz2hsv_
  generic, public :: xyz2hsl => xyz2hsl_
! 
  generic, public :: luv2rgb => luv2rgb_
  generic, public :: luv2xyz => luv2xyz_
  generic, public :: luv2lab => luv2lab_
  generic, public :: luv2hsv => luv2hsv_
  generic, public :: luv2hsl => luv2hsl_
! 
  generic, public :: lab2rgb => lab2rgb_
  generic, public :: lab2xyz => lab2xyz_
  generic, public :: lab2luv => lab2luv_
  generic, public :: lab2hsv => lab2hsv_
  generic, public :: lab2hsl => lab2hsl_

  generic, public :: hsv2rgb => hsv2rgb_
  generic, public :: hsv2xyz => hsv2xyz_
  generic, public :: hsv2luv => hsv2luv_
  generic, public :: hsv2lab => hsv2lab_
  generic, public :: hsv2hsl => hsv2hsl_

  generic, public :: hsl2rgb => hsl2rgb_
  generic, public :: hsl2xyz => hsl2xyz_
  generic, public :: hsl2luv => hsl2luv_
  generic, public :: hsl2lab => hsl2lab_
  generic, public :: hsl2hsv => hsl2hsv_

  generic, public :: sph2rgb => sph2rgb_
  generic, public :: sphere2rgb  => sphere2rgb_
end type colorspace_T

! the constructor routine for this class 
interface colorspace_T
  module procedure colorspace_constructor
end interface colorspace_T

contains

!--------------------------------------------------------------------------
type(colorspace_T) function colorspace_constructor( Nfold ) result(colorspace)
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!
!! constructor for the colorspace_T Class

use mod_IO 

IMPLICIT NONE

integer(kind=irg),INTENT(IN),OPTIONAL     :: Nfold 

type(IO_T)                                :: Message 


colorspace%Standards%sRGBmat = rgbMat_(colorspace%Standards%sRGB, colorspace%Standards%D65_2)
colorspace%Standards%sRGBmatInv = inv3x3(colorspace%Standards%sRGBmat)

if (present(Nfold)) then 
  colorspace%Nfold = Nfold 
  select case(Nfold) 
    case(4) 
      allocate(colorspace%P(0:14,0:2))
      colorspace%N = 15
      colorspace%P = transpose(P4)
    case(6) 
      allocate(colorspace%P(0:20,0:2))
      colorspace%N = 21
      colorspace%P = transpose(P6)
    case default
      call Message%printError('colorspace_constructor','unknown Nfold value')
  end select
end if 

end function colorspace_constructor

!--------------------------------------------------------------------------
subroutine colorspace_destructor(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: colorspace_destructor
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!
!! destructor for the colorspace_T Class
 
IMPLICIT NONE

type(colorspace_T), INTENT(INOUT)  :: self 

call reportDestructor('colorspace_T')

end subroutine colorspace_destructor

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
recursive function inv3x3(mat) result(invmat)
!DEC$ ATTRIBUTES DLLEXPORT :: inv3x3
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!
!! invert a 3x3 matrix (uses the mod_math routine mInvert)

use mod_math

IMPLICIT NONE 

real(kind=dbl),INTENT(IN)     :: mat(3,3)
real(kind=dbl)                :: invmat(3,3)

call mInvert(mat, invmat, .FALSE.)

end function inv3x3

!--------------------------------------------------------------------------
recursive function rgbMat_(rgb, w) result(rgbmat)
!DEC$ ATTRIBUTES DLLEXPORT :: rgbMat_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!
!! compute a matrix to convert from XYZ --> rgb
!! rgb: chromaticity of {red point, green point, blue point} as XYZ
!! w: chromaticity of white point as XYZ

IMPLICIT NONE 

real(kind=dbl),INTENT(IN)     :: rgb(0:2,0:2)
real(kind=dbl),INTENT(IN)     :: w(0:2)
real(kind=dbl)                :: rgbmat(0:2,0:2)

real(kind=dbl)                :: WW(0:2)
real(kind=dbl)                :: invR(0:2,0:2), invG(0:2,0:2), invB(0:2,0:2)

WW = w/w(1)

! build and invert 3x3 matrices to solve for rows of conversion matrix
! matrix * (r, g, b) = {x, 0, 0}, {0, x, 0}, {0, x, 0} and matrix^-1 * {1,1,1} = w
invR = inv3x3( transpose(reshape( (/  WW(0) ,   WW(1) ,   WW(2) , rgb(1,0), rgb(1,1), rgb(1,2), &
                                    rgb(2,0), rgb(2,1), rgb(2,2)/), (/3,3/))))
invG = inv3x3( transpose(reshape( (/rgb(0,0), rgb(0,1), rgb(0,2),   WW(0) ,   WW(1) ,   WW(2) , &
                                    rgb(2,0), rgb(2,1), rgb(2,2)/), (/3,3/))))
invB = inv3x3( transpose(reshape( (/rgb(0,0), rgb(0,1), rgb(0,2), rgb(1,0), rgb(1,1), rgb(1,2), &
                                      WW(0) ,   WW(1) ,   WW(2) /), (/3,3/))))

! assemble matrix
rgbmat = transpose(reshape( (/invR(0,0),invR(1,0),invR(2,0),invG(0,1),invG(1,1),invG(2,1),invB(0,2),invB(1,2),invB(2,2)/), &
                            (/ 3,3 /) ) )

end function rgbMat_

!--------------------------------------------------------------------------
recursive function hcm2rgb_(h, c, m) result(rgb)
!DEC$ ATTRIBUTES DLLEXPORT :: hcm2rgb_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!
!! convert from hue, chroma, and minimum to rgb

IMPLICIT NONE 

real(kind=dbl),INTENT(IN)     :: h
real(kind=dbl),INTENT(IN)     :: c
real(kind=dbl),INTENT(IN)     :: m
real(kind=dbl)                :: rgb(0:2)

real(kind=dbl)                :: x
integer(kind=irg)             :: h6

h6 = int(h * 6)   ! remember that Hue is in the range (0,1) 
x = c * (1.D0 - abs(mod(h*6.0D0, 2.D0) - 1.D0))

select case (h6) 
  case (0) 
    rgb = (/ c+m, x+m,   m /)
  case (1) 
    rgb = (/ x+m, c+m,   m /)
  case (2) 
    rgb = (/   m, c+m, x+m /)
  case (3) 
    rgb = (/   m, x+m, c+m /)
  case (4) 
    rgb = (/ x+m,   m, c+m /)
  case (5) 
    rgb = (/ c+m,   m, x+m /)
  case default 
    rgb = (/ 0.D0, 0.D0, 0.D0 /)
end select 

end function hcm2rgb_

!--------------------------------------------------------------------------
!------------------Implementation of Direct Conversions--------------------
!--------------------------------------------------------------------------
recursive function xyz2rgb_(self, xyz) result(rgb)
!DEC$ ATTRIBUTES DLLEXPORT :: xyz2rgb_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!
!! convert from XYZ to sRGB
!! xyz      : XYZ (X, Y, Z) values to convert
!! rgb      : location to write sRGB (red, green, blue) values
!! return   : true/false if xyz falls outside/inside the sRGB color gamut (removed on 9/8/21)
!!
!! verified on 9/7/21

IMPLICIT NONE 

class(colorspace_T)             :: self
real(kind=dbl),INTENT(IN)       :: xyz(0:2)
real(kind=dbl)                  :: rgb(0:2)

real(kind=dbl)                  :: gammaInv
real(kind=dbl)                  :: k0Inv
real(kind=dbl)                  :: a1
real(kind=dbl)                  :: work(0:2)
integer(kind=irg)               :: i 
! logical                         :: clamped

gammaInv = 1.D0 / self%Standards%sGamma
k0Inv    = self%Standards%sK0 / self%Standards%sPhi
a1       = 1.D0 + self%Standards%sA

work = matmul(self%Standards%sRGBmat, xyz)

! gamma correction
do i=0,2
  if (work(i).le.k0Inv) then 
    rgb(i) = work(i) * self%Standards%sPhi
  else
    rgb(i) = a1 * work(i)**gammaInv - self%Standards%sA
  end if
end do 

! check if this value is outside the sRGB color gamut
! clamped = .FALSE.
! do i=0,2
!   if(rgb(i).lt.0.0) then
!     rgb(i) = 0.0
!     clamped = .TRUE.
!   else if(rgb(i).gt.1.0) then
!     rgb(i) = 1.0
!     clamped = .TRUE.
!   end if 
! end do

end function xyz2rgb_

!--------------------------------------------------------------------------
recursive function rgb2xyz_(self, rgb) result(xyz) 
!DEC$ ATTRIBUTES DLLEXPORT :: rgb2xyz_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!
!! convert from sRGB to XYZ
!! rgb      : sRGB (red, green, blue) values to convert
!! xyz      : XYZ (X, Y, Z) values 
!!
!! verified on 9/7/21

IMPLICIT NONE 

class(colorspace_T)             :: self
real(kind=dbl),INTENT(IN)       :: rgb(0:2)
real(kind=dbl)                  :: xyz(0:2)

real(kind=dbl)                  :: a1
integer(kind=irg)               :: i 

a1 = 1.D0 + self%Standards%sA

! undo the gamma correction
do i=0,2
  if (rgb(i).le.self%Standards%sK0) then 
    xyz(i) = rgb(i) / self%Standards%sPhi
  else
    xyz(i) = ( (rgb(i) + self%Standards%sA)/a1)**self%Standards%sGamma
  end if
end do 

xyz = matmul(self%Standards%sRGBmatInv, xyz)

end function rgb2xyz_

!--------------------------------------------------------------------------
recursive function xyz2lab_(self, xyz, ill) result(lab)
!DEC$ ATTRIBUTES DLLEXPORT :: xyz2lab_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!
!! convert from XYZ to Lab
!! xyz: XYZ (X, Y, Z) values to convert
!! lab: location to write Lab (L*, a*, b*) values
!! ill: Lab illuminant as XYZ (or absent to use illuminant D65 for a 2 degree observer)
!!
!! verified on 9/7/21

IMPLICIT NONE 

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: xyz(0:2)
real(kind=dbl),INTENT(IN),OPTIONAL  :: ill(0:2)
real(kind=dbl)                      :: lab(0:2)

real(kind=dbl)                      :: delta, d2, d3, k, fXYZ(0:2), illum(0:2), t, inv3
integer(kind=irg)                   :: i 

delta = 6.D0 / 29.D0
d2 = delta * delta * 3.D0
d3 = delta * delta * delta
k = 4.D0 / 29.D0
inv3 = 1.D0/3.D0

if (present(ill)) then 
  illum = ill
else
  illum = self%Standards%D65_2
end if 
illum = illum/illum(1)

! compute f(i/i_N)
do i=0,2
  t = xyz(i) / illum(i)  ! normalize with respect to white point
  if (t.gt.d3) then 
    fXYZ(i) = t**inv3 
  else  
    fXYZ(i) = t / d2 + k ! apply nonlinear scaling
  end if 
end do

! change basis and scale
lab = (/ fXYZ(1) * 116.D0 - 16.D0, (fXYZ(0) - fXYZ(1)) * 500.D0, (fXYZ(1) - fXYZ(2)) * 200.D0 /)

end function xyz2lab_

!--------------------------------------------------------------------------
recursive function lab2xyz_(self, lab, ill) result(xyz)
!DEC$ ATTRIBUTES DLLEXPORT :: lab2xyz_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!
!! convert from Lab to xyz
!! lab: Lab (L*, a*, b*) values
!! xyz: XYZ (X, Y, Z) values 
!! ill: Lab illuminant as XYZ (or absent to use illuminant D65 for a 2 degree observer)
!!
!! verified on 9/7/21

IMPLICIT NONE 

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: lab(0:2)
real(kind=dbl),INTENT(IN),OPTIONAL  :: ill(0:2)
real(kind=dbl)                      :: xyz(0:2)

real(kind=dbl)                      :: delta, d2, d3, k, fXYZ(0:2), illum(0:2), t
integer(kind=irg)                   :: i 

delta = 6.D0 / 29.D0
d2 = delta * delta * 3.D0
d3 = delta * delta * delta
k = 4.D0 / 29.D0

if (present(ill)) then 
  illum = ill
else
  illum = self%Standards%D65_2
end if 
illum = illum/illum(1)

! change basis and scale
t = (lab(0)+16.D0)/116.D0 
xyz = (/ t + lab(1)/500.D0, t, t-lab(2)/200.D0 /)

! compute f(i/i_N)
do i=0,2
  if (xyz(i).gt.delta) then 
    xyz(i) = xyz(i)**3 
  else  
    xyz(i) = d2 * ( xyz(i) - k ) ! remove nonlinear scaling
  end if 
end do
xyz = xyz * illum  ! remove normalization with respect to white point

end function lab2xyz_

!--------------------------------------------------------------------------
recursive function xyz2luv_(self, xyz, ill) result(luv)
!DEC$ ATTRIBUTES DLLEXPORT :: xyz2luv_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!
!! convert from XYZ to Luv
!! xyz: XYZ (X, Y, Z) values 
!! luv: luv (L*, u*, v*) values
!! ill: Lab illuminant as XYZ (or absent to use illuminant D65 for a 2 degree observer)
!!
!! verified on 9/7/21

IMPLICIT NONE 

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: xyz(0:2)
real(kind=dbl),INTENT(IN),OPTIONAL  :: ill(0:2)
real(kind=dbl)                      :: luv(0:2)

real(kind=dbl)                      :: d, d_8, denn, den, u, v, inv3, illum(0:2) 
logical                             :: zero  

d   = 216.D0 / 24389.D0 ! (6/29)^3
d_8 = 8.D0 / d         ! (29/3)^3
inv3 = 1.D0/3.D0 

if (present(ill)) then 
  illum = ill
else
  illum = self%Standards%D65_2
end if 

! compute normalized X, Y, and Z and denomentators and u/v
denn = (illum(0) + illum(1) * 15.D0 + illum(2) * 3.D0) / illum(1)
den  =  xyz(0) + xyz(1) * 15.D0 + xyz(2) * 3.D0
if (den.eq.0.D0) then 
  u =  - (illum(0) / illum(1)) / denn
  v =  - 1.D0 / denn
else 
  u = xyz(0) / den - (illum(0) / illum(1)) / denn
  v = xyz(1) / den - 1.D0 / denn
end if 

! compute Luv
if (xyz(1).le.d) then 
  luv(0) = xyz(1) * d_8
else 
  luv(0) = xyz(1)**inv3 * 116.D0 - 16.D0
end if 
luv(1) = luv(0) *  52.D0 * u
luv(2) = luv(0) * 117.D0 * v

end function xyz2luv_

!--------------------------------------------------------------------------
recursive function luv2xyz_(self, luv, ill) result(xyz)
!DEC$ ATTRIBUTES DLLEXPORT :: luv2xyz_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!
!! convert from Luv to xyz
!! luv: luv (L*, u*, v*) values
!! xyz: XYZ (X, Y, Z) values 
!! ill: Lab illuminant as XYZ (or absent to use illuminant D65 for a 2 degree observer)
!!
!! verified on 9/7/21

IMPLICIT NONE 

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: luv(0:2)
real(kind=dbl),INTENT(IN),OPTIONAL  :: ill(0:2)
real(kind=dbl)                      :: xyz(0:2)

real(kind=dbl)                      :: d, denn, up, vp, Lp, illum(0:2) 
logical                             :: zero  

if (luv(0).eq.0.D0) then 
  xyz = (/ 0.D0, 0.D0, 0.D0 /)
else
  d = 27.D0 / 24389.D0 
  if (present(ill)) then 
    illum = ill
  else
    illum = self%Standards%D65_2
  end if 

! compute u' and v'
  denn = (illum(0) + illum(1) * 15.D0 + illum(2) * 3.D0) / illum(1)
  up = (luv(1) / 13.D0 + luv(0) * (illum(0) / illum(1)) * 4.D0 / denn) * 3.D0
  vp = (luv(2) / 13.D0 + luv(0) * 9.D0 / denn) * 4.D0
  Lp = (luv(0) + 16.D0) / 116.D0

! compute X, Y, and Z
  if (luv(0).le.8.D0) then 
    xyz(1) = luv(0) * d 
  else 
    xyz(1) = Lp**3
  end if 
  xyz(2) = xyz(1) * (12.D0 * luv(0) - up - vp * 5.D0) / vp
  xyz(0) = xyz(1) * (up * 3.D0) / vp

end if

end function luv2xyz_

!--------------------------------------------------------------------------
recursive function hsv2rgb_(self, hsv) result(rgb)
!DEC$ ATTRIBUTES DLLEXPORT :: hsv2rgb_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!
!! convert from hsv to rgb
!! hsv: hsv (hue, saturation, value) values to convert
!! rgb: location to write rgb (red, green, blue) values
!!
!! verified on 9/7/21

IMPLICIT NONE 

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: hsv(0:2)
real(kind=dbl)                      :: rgb(0:2)

real(kind=dbl)                      :: c

c = hsv(1) * hsv(2)  ! chroma
rgb = hcm2rgb_( hsv(0), c, hsv(2) - c )

end function hsv2rgb_

!--------------------------------------------------------------------------
recursive function hsl2rgb_(self, hsl) result(rgb)
!DEC$ ATTRIBUTES DLLEXPORT :: hsl2rgb_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!
!! convert from hsl to rgb
!! hsv: hsl: hsl (hue, saturation, lightness) values to convert
!! rgb: location to write rgb (red, green, blue) values
!!
!! verified on 9/7/21

IMPLICIT NONE 

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: hsl(0:2)
real(kind=dbl)                      :: rgb(0:2)

real(kind=dbl)                      :: c

c = (1.D0 - abs(hsl(2) * 2.D0 - 1.D0)) * hsl(1)  ! chroma
rgb = hcm2rgb_( hsl(0), c, hsl(2) - c*0.5D0 )

end function hsl2rgb_

!--------------------------------------------------------------------------
recursive function hsl2hsv_(self, hsl) result(hsv)
!DEC$ ATTRIBUTES DLLEXPORT :: hsl2hsv_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!
!! convert from hsl to hsv
!! hsl: hsl (hue, saturation, lightness) values to convert
!! hsv: location to write hsv (hue, saturation, value) values
!!
!! verified on 9/7/21

IMPLICIT NONE 

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: hsl(0:2)
real(kind=dbl)                      :: hsv(0:2)

real(kind=dbl)                      :: s

if (hsl(2).lt.0.5D0) then 
  s = hsl(1) * hsl(2)
else 
  s = hsl(1) * ( 1.D0 - hsl(2)) 
end if 

hsv(0) = hsl(0)
hsv(2) = s + hsl(2)
if (s.eq.0.D0) then 
  hsv(1) = 0.D0
else 
  hsv(1) = 2.D0 * s / hsv(2)
end if 

end function hsl2hsv_

!--------------------------------------------------------------------------
recursive function hsv2hsl_(self, hsv) result(hsl)
!DEC$ ATTRIBUTES DLLEXPORT :: hsv2hsl_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!
!! convert from hsv to hsl
!! hsv: hsv (hue, saturation, lightness) values to convert
!! hsl: location to write hsl (hue, saturation, value) values
!!
!! verified on 9/7/21

IMPLICIT NONE 

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: hsv(0:2)
real(kind=dbl)                      :: hsl(0:2)

real(kind=dbl)                      :: sv, x

sv = hsv(1) * hsv(2)
x = 2.D0 * hsv(2) - sv

hsl(0) = hsv(0)
hsl(2) = hsv(2) - sv / 2.D0
if (sv.eq.0.D0) then 
  hsl(1) = 0.D0
else 
  if (x.lt.1.D0) then 
    hsl(1) = sv / x 
  else
    hsl(1) = sv / ( 2.D0 - x )
  end if 
end if 

end function hsv2hsl_

!--------------------------------------------------------------------------
recursive function rgb2hsv_(self, rgb) result(hsv)
!DEC$ ATTRIBUTES DLLEXPORT :: rgb2hsv_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!
!! convert rgb to hsv
!! rgb: sRGB (red, green, blue) values to convert
!! hsv: location to write hsv (hue, saturation, value) values
!!
!! verified on 9/7/21

IMPLICIT NONE 

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: rgb(0:2)
real(kind=dbl)                      :: hsv(0:2)

real(kind=dbl)                      :: vMin, vMax, delta
integer(kind=irg)                   :: iMax(1)

vMax = maxval(rgb)

if(vMax.eq.0.D0) then ! black
  hsv = 0.D0
else 
  vMin = minval(rgb)
  delta = vMax - vMin
  if (delta.eq.0.D0) then ! gray
    hsv(0:1) = 0.D0
  else 
    iMax = maxloc(rgb)
    select case (iMax(1))
      case(1)
        hsv(0) = mod( (rgb(1)-rgb(2))/delta + 12.D0, 6.D0 )
      case(2)
        hsv(0) = (rgb(2)-rgb(0))/delta + 2.D0
      case(3)
        hsv(0) = (rgb(0)-rgb(1))/delta + 4.D0
    end select
    if (hsv(0).lt.0.0) hsv(0) = 1.D0 + hsv(0)
    hsv(1) = delta / vMax
  end if 
  hsv(2) = vMax
end if 
hsv(0) = hsv(0) / 6.D0

end function rgb2hsv_

!--------------------------------------------------------------------------
recursive function rgb2hsl_(self, rgb) result(hsl)
!DEC$ ATTRIBUTES DLLEXPORT :: rgb2hsl_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!
!! convert rgb to hsv
!! rgb: sRGB (red, green, blue) values to convert
!! hsv: location to write hsv (hue, saturation, value) values
!!
!! verified on 9/7/21

IMPLICIT NONE 

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: rgb(0:2)
real(kind=dbl)                      :: hsl(0:2)

real(kind=dbl)                      :: vMin, vMax, delta, sigma
integer(kind=irg)                   :: iMax(1) 

vMax = maxval(rgb)

if(vMax.eq.0.D0) then ! black
  hsl = 0.D0
else 
  vMin = minval(rgb)
  delta = vMax - vMin
  sigma = vMax + vMin
  if (delta.eq.0.D0) then ! gray
    hsl(0:1) = 0.D0
  else 
    iMax = maxloc(rgb)
    select case (iMax(1))
      case(1)
        hsl(0) = mod( (rgb(1)-rgb(2))/delta + 12.D0, 6.D0 )
      case(2)
        hsl(0) = (rgb(2)-rgb(0))/delta + 2.D0
      case(3)
        hsl(0) = (rgb(0)-rgb(1))/delta + 4.D0
    end select
    if (hsl(0).lt.0.0) hsl(0) = 1.D0 + hsl(0)
    if (sigma.lt.1.D0) then 
      hsl(1) = delta / sigma
    else 
      hsl(1) = delta / ( 2.D0 - sigma )
    end if 
  end if 
  hsl(2) = 0.5D0 * sigma
end if 
hsl(0) = hsl(0) / 6.D0

end function rgb2hsl_

!--------------------------------------------------------------------------
!------------------Implementation of Indirect Conversions------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! illuminant free cie spaces -> hsl/hsv (using cie spaces -> rgb)
!--------------------------------------------------------------------------
recursive function xyz2hsv_( self, xyz) result(hsv) 
!DEC$ ATTRIBUTES DLLEXPORT :: xyz2hsv_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!

IMPLICIT NONE 

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: xyz(0:2)
real(kind=dbl)                      :: hsv(0:2)

hsv = self%xyz2rgb_(xyz) 
hsv = self%rgb2hsv_(hsv) ! xyz->rgb->hsv

end function xyz2hsv_

!--------------------------------------------------------------------------
recursive function xyz2hsl_( self, xyz) result(hsl) 
!DEC$ ATTRIBUTES DLLEXPORT :: xyz2hsl_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!

IMPLICIT NONE 

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: xyz(0:2)
real(kind=dbl)                      :: hsl(0:2)

hsl = self%xyz2rgb_(xyz)
hsl = self%rgb2hsl_(hsl) ! xyz->rgb->hsl

end function xyz2hsl_

!--------------------------------------------------------------------------
! illuminated cie spaces -> rgb/hsv/hsl (using cie spaces -> xyz -> rgb)
!--------------------------------------------------------------------------
recursive function luv2rgb_( self, luv, ill) result(rgb) 
!DEC$ ATTRIBUTES DLLEXPORT :: luv2rgb_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!

IMPLICIT NONE

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: luv(0:2)
real(kind=dbl),INTENT(IN),OPTIONAL  :: ill(0:2)
real(kind=dbl)                      :: rgb(0:2)

real(kind=dbl)                      :: illum(0:2) 

if (present(ill)) then 
  illum = ill
else
  illum = self%Standards%D65_2
end if 

rgb = self%luv2xyz_(luv, illum)
rgb = self%xyz2rgb_(rgb)  ! luv->xyz->rgb

end function luv2rgb_

!--------------------------------------------------------------------------
recursive function luv2hsv_( self, luv, ill) result(hsv) 
!DEC$ ATTRIBUTES DLLEXPORT :: luv2hsv_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!

IMPLICIT NONE

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: luv(0:2)
real(kind=dbl),INTENT(IN),OPTIONAL  :: ill(0:2)
real(kind=dbl)                      :: hsv(0:2)

real(kind=dbl)                      :: illum(0:2) 

if (present(ill)) then 
  illum = ill
else
  illum = self%Standards%D65_2
end if 

hsv = self%luv2xyz_(luv, illum)
hsv = self%xyz2hsv_(hsv)    ! luv->xyz->rgb->hsv

end function luv2hsv_

!--------------------------------------------------------------------------
recursive function luv2hsl_( self, luv, ill) result(hsl) 
!DEC$ ATTRIBUTES DLLEXPORT :: luv2hsl_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!

IMPLICIT NONE

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: luv(0:2)
real(kind=dbl),INTENT(IN),OPTIONAL  :: ill(0:2)
real(kind=dbl)                      :: hsl(0:2)

real(kind=dbl)                      :: illum(0:2) 

if (present(ill)) then 
  illum = ill
else
  illum = self%Standards%D65_2
end if 

hsl = self%luv2xyz_(luv, illum)
hsl = self%xyz2hsl_(hsl)    ! luv->xyz->rgb->hsl

end function luv2hsl_

!--------------------------------------------------------------------------
recursive function lab2rgb_( self, lab, ill) result(rgb) 
!DEC$ ATTRIBUTES DLLEXPORT :: lab2rgb_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!

IMPLICIT NONE

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: lab(0:2)
real(kind=dbl),INTENT(IN),OPTIONAL  :: ill(0:2)
real(kind=dbl)                      :: rgb(0:2)

real(kind=dbl)                      :: illum(0:2) 

if (present(ill)) then 
  illum = ill
else
  illum = self%Standards%D65_2
end if 

rgb = self%lab2xyz_(lab, illum)
rgb = self%xyz2rgb_(rgb)    !  lab->xyz->rgb

end function lab2rgb_

!--------------------------------------------------------------------------
recursive function lab2hsv_( self, lab, ill) result(hsv) 
!DEC$ ATTRIBUTES DLLEXPORT :: lab2hsv_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!

IMPLICIT NONE

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: lab(0:2)
real(kind=dbl),INTENT(IN),OPTIONAL  :: ill(0:2)
real(kind=dbl)                      :: hsv(0:2)

real(kind=dbl)                      :: illum(0:2) 

if (present(ill)) then 
  illum = ill
else
  illum = self%Standards%D65_2
end if 

hsv = self%lab2xyz_(lab, illum)
hsv = self%xyz2hsv_(hsv)     ! lab->xyz->rgb->hsv

end function lab2hsv_

!--------------------------------------------------------------------------
recursive function lab2hsl_( self, lab, ill) result(hsl) 
!DEC$ ATTRIBUTES DLLEXPORT :: lab2hsl_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!

IMPLICIT NONE

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: lab(0:2)
real(kind=dbl),INTENT(IN),OPTIONAL  :: ill(0:2)
real(kind=dbl)                      :: hsl(0:2)

real(kind=dbl)                      :: illum(0:2) 

if (present(ill)) then 
  illum = ill
else
  illum = self%Standards%D65_2
end if 

hsl = self%lab2xyz_(lab, illum)
hsl = self%xyz2hsl_(hsl)     !  lab->xyz->rgb->hsl

end function lab2hsl_

!--------------------------------------------------------------------------
! within cie spaces (through xyz)
!--------------------------------------------------------------------------
recursive function luv2lab_( self, luv, ill) result(lab) 
!DEC$ ATTRIBUTES DLLEXPORT :: luv2lab_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!

IMPLICIT NONE

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: luv(0:2)
real(kind=dbl),INTENT(IN),OPTIONAL  :: ill(0:2)
real(kind=dbl)                      :: lab(0:2)

real(kind=dbl)                      :: illum(0:2) 

if (present(ill)) then 
  illum = ill
else
  illum = self%Standards%D65_2
end if 

lab = self%luv2xyz_(luv, illum)
lab = self%xyz2lab_(lab, illum)  ! luv->xyz->lab

end function luv2lab_

!--------------------------------------------------------------------------
recursive function lab2luv_( self, lab, ill) result(luv) 
!DEC$ ATTRIBUTES DLLEXPORT :: lab2luv_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!

IMPLICIT NONE

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: lab(0:2)
real(kind=dbl),INTENT(IN),OPTIONAL  :: ill(0:2)
real(kind=dbl)                      :: luv(0:2)

real(kind=dbl)                      :: illum(0:2) 

if (present(ill)) then 
  illum = ill
else
  illum = self%Standards%D65_2
end if 

luv = self%lab2xyz_(lab, illum) 
luv = self%xyz2luv_(luv, illum)   ! lab->xyz->luv

end function lab2luv_

!--------------------------------------------------------------------------
! rgb/hsv/hsl to cie spaces (all go through rgb -> xyz)
!--------------------------------------------------------------------------
recursive function rgb2luv_( self, rgb, ill) result(luv) 
!DEC$ ATTRIBUTES DLLEXPORT :: rgb2luv_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!

IMPLICIT NONE

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: rgb(0:2)
real(kind=dbl),INTENT(IN),OPTIONAL  :: ill(0:2)
real(kind=dbl)                      :: luv(0:2)

real(kind=dbl)                      :: illum(0:2) 

if (present(ill)) then 
  illum = ill
else
  illum = self%Standards%D65_2
end if 

luv = self%rgb2xyz_(rgb) 
luv = self%xyz2luv_(luv, illum)  ! rgb->xyz->luv

end function rgb2luv_

!--------------------------------------------------------------------------
recursive function rgb2lab_( self, rgb, ill) result(lab) 
!DEC$ ATTRIBUTES DLLEXPORT :: rgb2lab_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!
 
IMPLICIT NONE

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: rgb(0:2)
real(kind=dbl),INTENT(IN),OPTIONAL  :: ill(0:2)
real(kind=dbl)                      :: lab(0:2)

real(kind=dbl)                      :: illum(0:2) 

if (present(ill)) then 
  illum = ill
else
  illum = self%Standards%D65_2
end if 

lab = self%rgb2xyz_(rgb) 
lab = self%xyz2lab_(lab, illum)   ! rgb->xyz->lab

end function rgb2lab_

!--------------------------------------------------------------------------
recursive function hsv2xyz_( self, hsv) result(xyz) 
!DEC$ ATTRIBUTES DLLEXPORT :: hsv2xyz_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!

IMPLICIT NONE
class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: hsv(0:2)
real(kind=dbl)                      :: xyz(0:2)

xyz = self%hsv2rgb_(hsv) 
xyz = self%rgb2xyz_(xyz)    ! hsv->rgb->xyz

end function hsv2xyz_

!--------------------------------------------------------------------------
recursive function hsv2luv_( self, hsv, ill) result(luv) 
!DEC$ ATTRIBUTES DLLEXPORT :: hsv2luv_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!

IMPLICIT NONE

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: hsv(0:2)
real(kind=dbl),INTENT(IN),OPTIONAL  :: ill(0:2)
real(kind=dbl)                      :: luv(0:2)

real(kind=dbl)                      :: illum(0:2) 

if (present(ill)) then 
  illum = ill
else
  illum = self%Standards%D65_2
end if 

luv = self%hsv2rgb_(hsv) 
luv = self%rgb2xyz_(luv) 
luv = self%xyz2luv_(luv, illum)    ! hsv->rgb->xyz->luv

end function hsv2luv_

!--------------------------------------------------------------------------
recursive function hsv2lab_( self, hsv, ill) result(lab) 
!DEC$ ATTRIBUTES DLLEXPORT :: hsv2lab_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!

IMPLICIT NONE

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: hsv(0:2)
real(kind=dbl),INTENT(IN),OPTIONAL  :: ill(0:2)
real(kind=dbl)                      :: lab(0:2)

real(kind=dbl)                      :: illum(0:2) 

if (present(ill)) then 
  illum = ill
else
  illum = self%Standards%D65_2
end if 

lab = self%hsv2rgb_(hsv) 
lab = self%rgb2xyz_(lab) 
lab = self%xyz2lab_(lab, illum)    ! hsv->rgb->xyz->lab

end function hsv2lab_

!--------------------------------------------------------------------------
recursive function hsl2xyz_( self, hsl) result(xyz) 
!DEC$ ATTRIBUTES DLLEXPORT :: hsl2xyz_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!

IMPLICIT NONE

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: hsl(0:2)
real(kind=dbl)                      :: xyz(0:2)

xyz = self%hsl2rgb_(hsl) 
xyz = self%rgb2xyz_(xyz)    ! hsl->rgb->xyz

end function hsl2xyz_

!--------------------------------------------------------------------------
recursive function hsl2luv_( self, hsl, ill) result(luv) 
!DEC$ ATTRIBUTES DLLEXPORT :: hsl2luv_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!

IMPLICIT NONE

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: hsl(0:2)
real(kind=dbl),INTENT(IN),OPTIONAL  :: ill(0:2)
real(kind=dbl)                      :: luv(0:2)

real(kind=dbl)                      :: illum(0:2) 

if (present(ill)) then 
  illum = ill
else
  illum = self%Standards%D65_2
end if 

luv = self%hsl2rgb_(hsl)
luv = self%rgb2xyz_(luv) 
luv = self%xyz2luv_(luv, illum)    ! hsl->rgb->xyz->luv

end function hsl2luv_

!--------------------------------------------------------------------------
recursive function hsl2lab_( self, hsl, ill) result(lab) 
!DEC$ ATTRIBUTES DLLEXPORT :: hsl2lab_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!

IMPLICIT NONE

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: hsl(0:2)
real(kind=dbl),INTENT(IN),OPTIONAL  :: ill(0:2)
real(kind=dbl)                      :: lab(0:2)

real(kind=dbl)                      :: illum(0:2) 

if (present(ill)) then 
  illum = ill
else
  illum = self%Standards%D65_2
end if 

lab = self%hsl2rgb_(hsl) 
lab = self%rgb2xyz_(lab) 
lab = self%xyz2lab_(lab, illum)    ! hsl->rgb->xyz->lab

end function hsl2lab_

!--------------------------------------------------------------------------
recursive function sph2rgb_( self, hsl) result(rgb) 
!DEC$ ATTRIBUTES DLLEXPORT :: sph2rgb_
!! author: MDG 
!! version: 1.0 
!! date: 09/07/21
!!

!! phenomenologically adjusted hsl2rgb (based on Will Lenthe's original code)

IMPLICIT NONE

class(colorspace_T)                 :: self
real(kind=dbl),INTENT(IN)           :: hsl(0:2)
real(kind=dbl)                      :: rgb(0:2)

logical                             :: whiteCenter, half 
real(kind=dbl)                      :: yL, yS, h3, h6, hNew, sP, th, gray, q, c, m, h, x 

!  1 / ( 1 + sqrt(2*pi) * std::erf( 5 * sqrt(2) / 3 ) * 0.3 );
real(kind=dbl), parameter           :: iDen = 0.570990316610288181236261564684297686279447800757942106831501845990856895D0
! sqrt(pi/2) / 10
real(kind=dbl), parameter           :: k1   = 0.125331413731550025120788264240552262650349337030496915831496178817114683D0
!10 * sqrt(2)
real(kind=dbl), parameter           :: k2   = 14.1421356237309504880168872420969807856967187537694807317667973799073248D0
! 1/3
real(kind=dbl), parameter           :: k1_3 = 0.333333333333333333333333333333333333333333333333333333333333333333333333D0
! 1/6
real(kind=dbl), parameter           :: k1_6 = 0.166666666666666666666666666666666666666666666666666666666666666666666667D0
! pi/2
real(kind=dbl), parameter           :: pi_2 = 1.57079632679489661923132169163975144209858469968755291048747229615390820D0

! get lightness and saturation rescaling parameters
whiteCenter = .FALSE.
if (hsl(2).ge.0.5D0) whiteCenter = .TRUE.
if (whiteCenter.eqv..TRUE.) then 
  yL = 0.25D0 
  yS = 0.20D0 
else 
  yL = 0.5D0 
  yS = 0.5D0 
end if 

! adjust hue gradient (A.5)
h3   = mod(hsl(0), k1_3)
half = .FALSE. 
if (h3.gt.k1_6) half = .TRUE. 
if (half.eqv..TRUE.) then 
  h6 = k1_3 - h3 
else 
  h6 = h3
end if 
hNew = (h6 + k1 * erf(k2 * h6)) * iDen
if (half.eqv..TRUE.) then  ! save adjusted hue
  rgb(0) = hsl(0) - h3 + k1_3 - hNew  
else 
  rgb(0) = hsl(0) - h3 + hNew
end if 

! adjust lightness gradient (A.9)
sP   = sin(hsl(2) * pi_2)
th   = yL * hsl(2) + (1.0D0 - yL) * sP * sP
gray = 1.0D0 - 2.0D0 * yS * abs(th - 0.5D0)
rgb(2) = (th - 0.5D0) * gray + 0.5D0  ! save adjusted lightness

! adjust saturation gradient (A.10)
q = 1.0D0 - abs( 2.0D0 * hsl(2) - 1.0D0 ) 
if (q.eq.0.D0) then 
  rgb(1) = 0.D0
else 
  rgb(1) = gray * ( 1.0D0 - abs( 2.0D0 * th - 1.0D0 ) ) / q  ! save adjusted saturation
end if 

! convert adjusted hsl to rgb
c = (1.D0 - abs(rgb(2) * 2.D0 - 1.D0)) * rgb(1)   ! compute chroma
m = rgb(2) - c/2.D0   ! m
h = rgb(0) * 6.D0     ! hue (0,1) -> (0,6)
x = c * (1.D0 - abs(mod(h, 2.D0) - 1.D0))

select case (int(h)) 
  case (0) 
    rgb = (/ c+m, x+m,   m /)
  case (1) 
    rgb = (/ x+m, c+m,   m /)
  case (2) 
    rgb = (/   m, c+m, x+m /)
  case (3) 
    rgb = (/   m, x+m, c+m /)
  case (4) 
    rgb = (/ x+m,   m, c+m /)
  case (5) 
    rgb = (/ c+m,   m, x+m /)
  case default 
    rgb = (/ 0.D0, 0.D0, 0.D0 /)
end select 

end function sph2rgb_

!--------------------------------------------------------------------------
recursive function sphere2rgb_( self, a, p, w0 ) result(rgb) 
!DEC$ ATTRIBUTES DLLEXPORT :: sphere2rgb_
!! author: MDG (converted to f90 from Will Lenthe's original C++ version) 
!! version: 1.0 
!! date: 09/28/21
!!
!! a  : fractional theta angle [0,1] (0.5 for the equator, 0 N, 1 S)
!! p  : fractional phi angle [0,1]
!! rgb: location to write rgb color
!! w0 : true/false for white/black @ center 

IMPLICIT NONE

class(colorspace_T)         :: self
real(kind=dbl),INTENT(IN)   :: a
real(kind=dbl),INTENT(IN)   :: p
logical,INTENT(IN)          :: w0
real(kind=dbl)              :: rgb(0:2)

logical                     :: sh, swap
real(kind=dbl)              :: h, l, t, tt, delta, work(0:2,0:2), w 
integer(kind=irg)           :: s, j 
real(kind=dbl)              :: tl, deltaS, deltaN, l0, l1, deltaL, x, hmt, c1, d1, fcfc, minL, maxL, &
                               a2, b2, c2, d2, ll, hpt, c4, d4, a3, b3, c3, d3, fc, ac, bc, cc, dc, mc, tc, luv(0:2)

h = p 
l = a 
if (w0.eqv..TRUE.) l = 1.D0-a 

! bring t to [0,1)
t = mod(h, 1.D0)
if (t.lt.0.D0) t = t+1.D0 

! remap t to knot domain and find segment t falls in
tt = t * (self%N-3) + 3
s = int(tt) 
delta = tt - s 

! unrolled de Boor's spline interpolation
! k = 0
w = delta / 3.D0
x = 1.D0 - w
do j = 0, 2
  work(2,j) = self%P(s, j) * w + self%P(s-1, j) * x
end do 

w = (delta + 1.D0) / 3.D0
x = 1.D0 - w
do j = 0, 2
  work(1,j) = self%P(s-1, j) * w + self%P(s-2, j) * x
end do 

w = (delta + 2.D0) / 3.D0
x = 1.D0 - w
do j = 0, 2
  work(0,j) = self%P(s-2, j) * w + self%P(s-3, j) * x
end do 

! k = 1
w = delta * 0.5D0
x = 1.D0 - w
do j = 0, 2
  work(2,j) = work(2, j) * w + work(1,j) * x
end do 

w = (delta + 1.D0) * 0.5D0
x = 1.D0 - w
do j = 0, 2
  work(1,j) = work(1, j) * w + work(0,j) * x
end do 

! k = 2
w = delta
x = 1.D0 - w
do j = 0, 2
  work(2,j) = work(2, j) * w + work(1,j) * x
end do 
 
! copy last row into luv
luv = (/ work(2,0), work(2,1), work(2,2) /) 

minL = 12.D0
maxL = 98.D0

! apply nonlinear rescaling to lightness interpolation to impose C2 continuity at equator
! We'll use a piecewise polynomial with linear portions near the poles to maximize truly perceptually uniform region
!         / f1 =                       c1 * x + d1 : 0       <= x <= 1/2 - t (linear from pole to first transition)
!  f(x) = | f2 = a2 * x^3 + b2 * x^2 + c2 * x + d2 : 1/2 - t <  x <= 1/2     (cubic from first transition to equator)
!         | f3 = a3 * x^3 + b3 * x^2 + c3 * x + d3 : 1/2     <  x <  1/2 + t (cubic from equator to second transition)
!         \ f4 =                       c4 * x + d4 : 1/2 + t <= x <= 1       (linear from second transition to pole)
! 12 unknowns solved for with the following 12 constraints
! -C2 continuity between f1/f2, f2/f3, and f3, f4 (9 constraints)
! -f(0) = minL, f(1/2) = luv(0), f(1) = maxL (3 constraints)
tl = 0.1D0 ! offset from equator to end of transition for C2 continuity of L* (0,0.5) 
! 0 -> true perceptually uniformity with visual discontinuity, 
! 0.5 -> largest deviation from perceptually uniformity but spreads discontinuity over largest area
sh = .FALSE. 
if (l.le.0.5D0) sh = .TRUE.           ! does this point fall below/above the equator (true/false)
deltaS = minL - luv(0)                ! delta luminance from equator to south pole
deltaN = maxL - luv(0)                ! delta luminance from equator to north pole

x = ( deltaS + deltaN ) / ( tl * 2.D0 - 3.D0 )    ! this value will be needed in all for segments
if (sh.eqv..TRUE.) then                   ! in f1 or f2, don't bother calculating coefficients for f3 or f4
  hmt = 0.5D0 - tl                        ! transition from linear -> cubic in southern cone
  c1 = (x * tl - l0) * 2.D0               ! compute slope of linear region in southern cone
  d1 = deltaS                             ! compute intercept of linear region in southern cone
  if (l.le.hmt) then                      ! in first linear region (f1)
    deltaL = c1 * l + d1                  ! compute f1(l)
  else                                    ! in first cubic region (f2)
    a2 = -x / (tl * tl)
    b2 = -a2 * hmt * 3.D0
    c2 =  c1 - b2 * hmt
    d2 =  d1 - a2 * hmt * hmt * hmt
    ll = l * l                            ! compute l^2 once
    deltaL = a2 * ll * l + b2 * ll + c2 * l + d2      ! compute f2(l)
  end if 
else                                      ! in f3 or f4, don't bother calculating coefficients for f1 or f2
  hpt = 0.5D0 + tl                        ! transition from cubic -> linear in northern cone
  c4 = (deltaN - x * tl) * 2.D0           ! compute slope of linear region in northern cone
  d4 = deltaN - c4                        ! compute intercept of linear region in northern cone
  if (l.ge.hpt) then                      ! in second linear region (f4)
    deltaL = c4 * l + d4                  ! compute f4(l)
  else                                    ! in second cubic region (f3)
    a3 =  x / (tl * tl)
    b3 = -a3 * hpt * 3.D0
    c3 =  c4 - b3 * hpt
    d3 =  d4 - a3 * hpt * hpt * hpt
    ll = l * l                            ! compute l^2 once
    deltaL = a3 * ll * l + b3 * ll + c3 * l + d3      ! compute f3(l)
  end if
end if

! interpolate luminance and compute scaling factor for chromaticity
luv(0) = luv(0) + deltaL
if (sh.eqv..TRUE.) then 
  fc = 1.D0 - deltaL / deltaS             ! 0->1 bicone scaling
else
  fc = 1.D0 - deltaL / deltaN             ! 0->1 bicone scaling
end if 

! adjust chromaticity scaling factor similarly to lightness to make chromaticity C1 continous at equator
!  f(x) = / f1 =                         x     : 0  < x <= tc (linear from pole to tc)
!         \ f2 = a * x^3 + b * x^2 + c * x + d : tc < x <= 1  (cubic from tc to equator)
! 4 unknowns solved for with C2 continuity between f1/f2 (3 constraints) and f'(1) = 0 (C1 continuity at equator)
tc = 0.8D0        ! this parameter is less sensitive than tl since local uniformity is more strongly L* dependant
ac = -1.D0 / ( (tc - 1.D0) * (tc - 1.D0) * 3.D0 )
bc = -tc * 3.D0 * ac
cc = (tc * 6.D0 - 3.D0) * ac
dc = -tc * tc * tc * ac
mc = 3.D0 / (tc + 2.D0)
if (fc.gt.tc) then 
  fcfc = fc * fc                                  ! compute fc^2 once
  fc = ac * fcfc * fc + bc * fcfc + cc * fc + dc  ! compute f2(fc)
end if 
luv(1) = luv(1) * fc 
luv(2) = luv(2) * fc

rgb = self%luv2rgb_(luv)   ! luv -> rgb

end function sphere2rgb_

end module mod_colorspace
