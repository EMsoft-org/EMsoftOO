program colortest

use mod_kinds
use mod_global
use mod_colorspace 

IMPLICIT NONE 

type(colorspace_T)      :: clr 
real(kind=sgl)          :: xyz(0:2), rgb(0:2), xyz2(0:2), lab(0:2), luv(0:2), hsv(0:2), hsl(0:2)
logical                 :: clamped 

! initialize the color space class 
clr = colorspace_T()

! rgb = (/ 1.0, 1.0, 1.0 /)
! write (*,*) 'rgb = ', rgb, ' ===> ', clr%rgb2hsl(rgb), clr%rgb2hsv(rgb)
! rgb = (/ 1.0, 0.0, 0.0 /)
! write (*,*) 'rgb = ', rgb, ' ===> ', clr%rgb2hsl(rgb), clr%rgb2hsv(rgb)
! rgb = (/ 0.75, 0.75, 0.0 /)
! write (*,*) 'rgb = ', rgb, ' ===> ', clr%rgb2hsl(rgb), clr%rgb2hsv(rgb)
! rgb = (/ 0.0, 0.50, 0.0 /)
! write (*,*) 'rgb = ', rgb, ' ===> ', clr%rgb2hsl(rgb), clr%rgb2hsv(rgb)
! rgb = (/ 0.50, 1.0, 1.0 /)
! write (*,*) 'rgb = ', rgb, ' ===> ', clr%rgb2hsl(rgb), clr%rgb2hsv(rgb)
! rgb = (/ 0.50, 0.50, 1.0 /)
! write (*,*) 'rgb = ', rgb, ' ===> ', clr%rgb2hsl(rgb), clr%rgb2hsv(rgb)
! rgb = (/ 0.75, 0.25, 0.75 /)
! write (*,*) 'rgb = ', rgb, ' ===> ', clr%rgb2hsl(rgb), clr%rgb2hsv(rgb)

rgb = (/ 0.495, 0.493, 0.721 /)
write (*,*) 'target  : ', rgb
write (*,*) 'rgb2xyz : ', clr%xyz2rgb(clr%rgb2xyz(rgb), clamped)
write (*,*) 'rgb2lab : ', clr%lab2rgb(clr%rgb2lab(rgb))
write (*,*) 'rgb2luv : ', clr%luv2rgb(clr%rgb2luv(rgb))
write (*,*) 'rgb2hsl : ', clr%hsl2rgb(clr%rgb2hsl(rgb))
write (*,*) 'rgb2hsv : ', clr%hsv2rgb(clr%rgb2hsv(rgb))

xyz = clr%rgb2xyz(rgb)
write (*,*) '---------'
write (*,*) 'target  : ', xyz 
write (*,*) 'xyz2rgb : ', clr%rgb2xyz(clr%xyz2rgb(xyz, clamped))
write (*,*) 'xyz2lab : ', clr%lab2xyz(clr%xyz2lab(xyz))
write (*,*) 'xyz2luv : ', clr%luv2xyz(clr%xyz2luv(xyz))
write (*,*) 'xyz2hsl : ', clr%hsl2xyz(clr%xyz2hsl(xyz, clamped))
write (*,*) 'xyz2hsv : ', clr%hsv2xyz(clr%xyz2hsv(xyz, clamped))

lab = clr%rgb2lab(rgb)
write (*,*) '---------'
write (*,*) 'target  : ', lab
write (*,*) 'lab2xyz : ', clr%xyz2lab(clr%lab2xyz(lab))
write (*,*) 'lab2rgb : ', clr%rgb2lab(clr%lab2rgb(lab))
write (*,*) 'lab2luv : ', clr%luv2lab(clr%lab2luv(lab))
write (*,*) 'lab2hsl : ', clr%hsl2lab(clr%lab2hsl(lab))
write (*,*) 'lab2hsv : ', clr%hsv2lab(clr%lab2hsv(lab))

luv = clr%rgb2luv(rgb)
write (*,*) '---------'
write (*,*) 'target  : ', luv
write (*,*) 'luv2xyz : ', clr%xyz2luv(clr%luv2xyz(luv))
write (*,*) 'luv2rgb : ', clr%rgb2luv(clr%luv2rgb(luv))
write (*,*) 'luv2lab : ', clr%lab2luv(clr%luv2lab(luv))
write (*,*) 'luv2hsl : ', clr%hsl2luv(clr%luv2hsl(luv))
write (*,*) 'luv2hsv : ', clr%hsv2luv(clr%luv2hsv(luv))

hsl = clr%rgb2hsl(rgb)
write (*,*) '---------'
write (*,*) 'target  : ', hsl
write (*,*) 'hsl2xyz : ', clr%xyz2hsl(clr%hsl2xyz(hsl), clamped)
write (*,*) 'hsl2rgb : ', clr%rgb2hsl(clr%hsl2rgb(hsl))
write (*,*) 'hsl2lab : ', clr%lab2hsl(clr%hsl2lab(hsl))
write (*,*) 'hsl2luv : ', clr%luv2hsl(clr%hsl2luv(hsl))
write (*,*) 'hsl2hsv : ', clr%hsv2hsl(clr%hsl2hsv(hsl))

hsv = clr%rgb2hsv(rgb)
write (*,*) '---------'
write (*,*) 'target  : ', hsv
write (*,*) 'hsv2xyz : ', clr%xyz2hsv(clr%hsv2xyz(hsv), clamped)
write (*,*) 'hsv2rgb : ', clr%rgb2hsv(clr%hsv2rgb(hsv))
write (*,*) 'hsv2lab : ', clr%lab2hsv(clr%hsv2lab(hsv))
write (*,*) 'hsv2luv : ', clr%luv2hsv(clr%hsv2luv(hsv))
write (*,*) 'hsv2hsl : ', clr%hsl2hsv(clr%hsv2hsl(hsv))


end program colortest