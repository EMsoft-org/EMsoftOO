! ###################################################################
! Copyright (c) 2016-2022, Marc De Graef Research Group/Carnegie Mellon University
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

module MODcolorspace
  !! author: MDG 
  !! version: 1.0 
  !! date: 09/08/21
  !!
  !! perform a series of unit tests on the colorspace module 

use stringconstants
use mod_global
use mod_kinds

contains 

subroutine MODcolorspaceExecuteTest(res) &
           bind(c, name='MODcolorspaceExecuteTest')    ! this routine is callable from a C/C++ program
!DEC$ ATTRIBUTES DLLEXPORT :: MODcolorspaceExecuteTest

use,INTRINSIC :: ISO_C_BINDING
use mod_kinds
use mod_colorspace
use mod_EMsoft

IMPLICIT NONE

integer(C_INT32_T),INTENT(OUT)  :: res

type(colorspace_T)              :: clr
type(EMsoft_T)                  :: EMsoft

character(fnlen)                :: progname = 'MODcolorspaceTest.f90'
character(fnlen)                :: progdesc = 'test program for the mod_colorspace module'

real(kind=dbl)                  :: xyzref(0:2), rgbref(0:2), labref(0:2), luvref(0:2), hsvref(0:2), hslref(0:2)
real(kind=dbl)                  :: d, dxyz, drgb, dlab, dluv, dhsv, dhsl
real(kind=dbl)                  :: xyz(0:2), rgb(0:2), lab(0:2), luv(0:2), hsv(0:2), hsl(0:2)
real(kind=dbl)                  :: xyz2(0:2), rgb2(0:2), lab2(0:2), luv2(0:2), hsv2(0:2), hsl2(0:2)
real(kind=dbl),parameter        :: eps = 1.D-12 

!===================================================
! set the reference values

rgbref = (/  0.49500000000000000D0,  0.49299999999999999D0,  0.72099999999999997D0 /)     
xyzref = (/  0.24695570762809679D0,  0.22755750948829487D0,  0.48352890683431898D0 /)     
labref = (/   54.819857172503987D0, 13.800164884434807D0,  -30.480481381104131D0 /)     
luvref = (/   54.819857172503987D0, -3.2440012157009277D0, -48.187749865274554D0 /)     
hslref = (/  0.66812865497076024D0,  0.29007633587786258D0,  0.60699999999999998D0 /)     
hsvref = (/  0.66812865497076024D0,  0.31622746185852979D0,  0.72099999999999997D0 /)     

! initialize the error identifier to zero (should remain zero upon successful exit)
res = 0

!===================================================
!===================Start Tests=====================
!===================================================

! execute the constructor
clr = colorspace_T() 

! initial tests 
rgb = (/ 0.495D0, 0.493D0, 0.721D0 /)
xyz = clr%rgb2xyz(rgb)
d = sqrt( sum( (xyzref-xyz)**2 ) )
if (d.ge.eps) then 
  res = 1
  write (*,"(' rgb2xyz  produces incorrect result ')")
end if 
lab = clr%rgb2lab(rgb)
d = sqrt( sum( (labref-lab)**2 ) )
if (d.ge.eps) then 
  res = 2
  write (*,"(' rgb2lab  produces incorrect result ')")
end if 
luv = clr%rgb2luv(rgb)
d = sqrt( sum( (luvref-luv)**2 ) )
if (d.ge.eps) then 
  res = 3
  write (*,"(' rgb2luv  produces incorrect result ')")
end if 
hsl = clr%rgb2hsl(rgb)
d = sqrt( sum( (hslref-hsl)**2 ) )
if (d.ge.eps) then 
  res = 4
  write (*,"(' rgb2hsl  produces incorrect result ')")
end if 
hsv = clr%rgb2hsv(rgb)
d = sqrt( sum( (hsvref-hsv)**2 ) )
if (d.ge.eps) then 
  res = 5
  write (*,"(' rgb2hsv  produces incorrect result ')")
end if 

! two-way tests 
xyz2 = clr%xyz2rgb(clr%rgb2xyz(rgb))
lab2 = clr%lab2rgb(clr%rgb2lab(rgb))
luv2 = clr%luv2rgb(clr%rgb2luv(rgb))
hsl2 = clr%hsl2rgb(clr%rgb2hsl(rgb))
hsv2 = clr%hsv2rgb(clr%rgb2hsv(rgb))
dxyz = sqrt( sum( (rgb-xyz2)**2 ) )
dlab = sqrt( sum( (rgb-lab2)**2 ) )
dluv = sqrt( sum( (rgb-luv2)**2 ) )
dhsl = sqrt( sum( (rgb-hsl2)**2 ) )
dhsv = sqrt( sum( (rgb-hsv2)**2 ) )

if (dxyz.ge.eps) then
  res = 6 
  write(*,"('clr%xyz2rgb(clr%rgb2xyz(rgb)) produces incorrect result')")
end if 
if (dlab.ge.eps) then
  res = 7 
  write(*,"('clr%lab2rgb(clr%rgb2lab(rgb)) produces incorrect result')")
end if 
if (dluv.ge.eps) then
  res = 8 
  write(*,"('clr%luv2rgb(clr%rgb2luv(rgb)) produces incorrect result')")
end if 
if (dhsl.ge.eps) then
  res = 9 
  write(*,"('clr%hsl2rgb(clr%rgb2hsl(rgb)) produces incorrect result')")
end if 
if (dhsv.ge.eps) then
  res = 10 
  write(*,"('clr%hsv2rgb(clr%rgb2hsv(rgb)) produces incorrect result')")
end if 


rgb2 = clr%rgb2xyz(clr%xyz2rgb(xyz))
lab2 = clr%lab2xyz(clr%xyz2lab(xyz))
luv2 = clr%luv2xyz(clr%xyz2luv(xyz))
hsl2 = clr%hsl2xyz(clr%xyz2hsl(xyz))
hsv2 = clr%hsv2xyz(clr%xyz2hsv(xyz))

drgb = sqrt( sum( (xyz-rgb2)**2 ) )
dlab = sqrt( sum( (xyz-lab2)**2 ) )
dluv = sqrt( sum( (xyz-luv2)**2 ) )
dhsl = sqrt( sum( (xyz-hsl2)**2 ) )
dhsv = sqrt( sum( (xyz-hsv2)**2 ) )
if (drgb.ge.eps) then 
  res = 11 
  write(*,"('clr%rgb2xyz(clr%xyz2rgb(xyz)) produces incorrect results')")
end if 
if (dlab.ge.eps) then 
  res = 12 
  write(*,"('clr%lab2xyz(clr%xyz2lab(xyz)) produces incorrect results')")
end if 
if (dluv.ge.eps) then 
  res = 13 
  write(*,"('clr%luv2xyz(clr%xyz2luv(xyz)) produces incorrect results')")
end if 
if (dhsl.ge.eps) then 
  res = 14 
  write(*,"('clr%hsl2xyz(clr%xyz2hsl(xyz)) produces incorrect results')")
end if 
if (dhsv.ge.eps) then 
  res = 15 
  write(*,"('clr%hsv2xyz(clr%xyz2hsv(xyz)) produces incorrect results')")
end if 

xyz2 = clr%xyz2lab(clr%lab2xyz(lab))
rgb2 = clr%rgb2lab(clr%lab2rgb(lab))
luv2 = clr%luv2lab(clr%lab2luv(lab))
hsl2 = clr%hsl2lab(clr%lab2hsl(lab))
hsv2 = clr%hsv2lab(clr%lab2hsv(lab))

dxyz = sqrt( sum( (lab-xyz2)**2 ) )
drgb = sqrt( sum( (lab-rgb2)**2 ) )
dluv = sqrt( sum( (lab-luv2)**2 ) )
dhsl = sqrt( sum( (lab-hsl2)**2 ) )
dhsv = sqrt( sum( (lab-hsv2)**2 ) )
if (dxyz.ge.eps) then 
  res = 16
  write(*,"('clr%xyz2lab(clr%lab2xyz(lab)) produces incorrect results')")
end if 
if (drgb.ge.eps) then 
  res = 17
  write(*,"('clr%rgb2lab(clr%lab2rgb(lab)) produces incorrect results')")
end if 
if (dluv.ge.eps) then 
  res = 18
  write(*,"('clr%luv2lab(clr%lab2luv(lab)) produces incorrect results')")
end if 
if (dhsl.ge.eps) then 
  res = 19
  write(*,"('clr%hsl2lab(clr%lab2hsl(lab)) produces incorrect results')")
end if 
if (dhsv.ge.eps) then 
  res = 20
  write(*,"('clr%hsv2lab(clr%lab2hsv(lab)) produces incorrect results')")
end if 

xyz2 = clr%xyz2luv(clr%luv2xyz(luv))
rgb2 = clr%rgb2luv(clr%luv2rgb(luv))
lab2 = clr%lab2luv(clr%luv2lab(luv))
hsl2 = clr%hsl2luv(clr%luv2hsl(luv))
hsv2 = clr%hsv2luv(clr%luv2hsv(luv))

dxyz = sqrt( sum( (luv-xyz2)**2 ) )
drgb = sqrt( sum( (luv-rgb2)**2 ) )
dlab = sqrt( sum( (luv-lab2)**2 ) )
dhsl = sqrt( sum( (luv-hsl2)**2 ) )
dhsv = sqrt( sum( (luv-hsv2)**2 ) )
if (dxyz.ge.eps) then 
  res = 21
  write(*,"('clr%xyz2luv(clr%luv2xyz(luv))')")
end if 
if (drgb.ge.eps) then 
  res = 22
  write(*,"('clr%rgb2luv(clr%luv2rgb(luv))')")
end if 
if (dlab.ge.eps) then 
  res = 23
  write(*,"('clr%lab2luv(clr%luv2lab(luv))')")
end if 
if (dhsl.ge.eps) then 
  res = 24
  write(*,"('clr%hsl2luv(clr%luv2hsl(luv))')")
end if 
if (dhsv.ge.eps) then 
  res = 25
  write(*,"('clr%hsv2luv(clr%luv2hsv(luv))')")
end if 

xyz2 = clr%xyz2hsl(clr%hsl2xyz(hsl))
rgb2 = clr%rgb2hsl(clr%hsl2rgb(hsl))
lab2 = clr%lab2hsl(clr%hsl2lab(hsl))
luv2 = clr%luv2hsl(clr%hsl2luv(hsl))
hsv2 = clr%hsv2hsl(clr%hsl2hsv(hsl))

dxyz = sqrt( sum( (hsl-xyz2)**2 ) )
drgb = sqrt( sum( (hsl-rgb2)**2 ) )
dlab = sqrt( sum( (hsl-lab2)**2 ) )
dluv = sqrt( sum( (hsl-luv2)**2 ) )
dhsv = sqrt( sum( (hsl-hsv2)**2 ) )
if (dxyz.ge.eps) then 
  res = 26
  write(*,"('clr%xyz2hsl(clr%hsl2xyz(hsl))')")
end if 
if (drgb.ge.eps) then 
  res = 27
  write(*,"('clr%rgb2hsl(clr%hsl2rgb(hsl))')")
end if 
if (dlab.ge.eps) then 
  res = 28
  write(*,"('clr%lab2hsl(clr%hsl2lab(hsl))')")
end if 
if (dluv.ge.eps) then 
  res = 29
  write(*,"('clr%luv2hsl(clr%hsl2luv(hsl))')")
end if 
if (dhsv.ge.eps) then 
  res = 30
  write(*,"('clr%hsv2hsl(clr%hsl2hsv(hsl))')")
end if 


xyz2 = clr%xyz2hsv(clr%hsv2xyz(hsv))
rgb2 = clr%rgb2hsv(clr%hsv2rgb(hsv))
lab2 = clr%lab2hsv(clr%hsv2lab(hsv))
luv2 = clr%luv2hsv(clr%hsv2luv(hsv))
hsl2 = clr%hsl2hsv(clr%hsv2hsl(hsv))

dxyz = sqrt( sum( (hsv-xyz2)**2 ) )
drgb = sqrt( sum( (hsv-rgb2)**2 ) )
dlab = sqrt( sum( (hsv-lab2)**2 ) )
dluv = sqrt( sum( (hsv-luv2)**2 ) )
dhsl = sqrt( sum( (hsv-hsl2)**2 ) )
if (dxyz.ge.eps) then
  res = 31
  write(*,"('clr%xyz2hsv(clr%hsv2xyz(hsv))')")
end if 
if (drgb.ge.eps) then
  res = 32
  write(*,"('clr%rgb2hsv(clr%hsv2rgb(hsv))')")
end if 
if (dlab.ge.eps) then
  res = 33
  write(*,"('clr%lab2hsv(clr%hsv2lab(hsv))')")
end if 
if (dluv.ge.eps) then
  res = 34
  write(*,"('clr%luv2hsv(clr%hsv2luv(hsv))')")
end if 
if (dhsl.ge.eps) then
  res = 35
  write(*,"('clr%hsl2hsv(clr%hsv2hsl(hsv))')")
end if 


end subroutine MODcolorspaceExecuteTest

end module MODcolorspace
