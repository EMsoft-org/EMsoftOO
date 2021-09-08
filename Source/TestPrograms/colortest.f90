program colortest

use mod_kinds
use mod_global
use mod_colorspace 

IMPLICIT NONE 

type(colorspace_T)      :: clr 
real(kind=dbl)          :: xyz(0:2), rgb(0:2), lab(0:2), luv(0:2), hsv(0:2), hsl(0:2)
real(kind=dbl)          :: xyz2(0:2), rgb2(0:2), lab2(0:2), luv2(0:2), hsv2(0:2), hsl2(0:2)

! initialize the color space class 
clr = colorspace_T()

rgb = (/ 0.495D0, 0.493D0, 0.721D0 /)
write (*,*) 'target  : ', rgb
xyz2 = clr%xyz2rgb(clr%rgb2xyz(rgb))
lab2 = clr%lab2rgb(clr%rgb2lab(rgb))
luv2 = clr%luv2rgb(clr%rgb2luv(rgb))
hsl2 = clr%hsl2rgb(clr%rgb2hsl(rgb))
hsv2 = clr%hsv2rgb(clr%rgb2hsv(rgb))
write (*,*) 'rgb2xyz : ', xyz2, sqrt( sum( (rgb-xyz2)**2 ) )
write (*,*) 'rgb2lab : ', lab2, sqrt( sum( (rgb-lab2)**2 ) )
write (*,*) 'rgb2luv : ', luv2, sqrt( sum( (rgb-luv2)**2 ) )
write (*,*) 'rgb2hsl : ', hsl2, sqrt( sum( (rgb-hsl2)**2 ) )
write (*,*) 'rgb2hsv : ', hsv2, sqrt( sum( (rgb-hsv2)**2 ) )

xyz = clr%rgb2xyz(rgb)
write (*,*) '---------'
write (*,*) 'target  : ', xyz 
rgb2 =clr%rgb2xyz(clr%xyz2rgb(xyz))
lab2 =clr%lab2xyz(clr%xyz2lab(xyz))
luv2 =clr%luv2xyz(clr%xyz2luv(xyz))
hsl2 =clr%hsl2xyz(clr%xyz2hsl(xyz))
hsv2 =clr%hsv2xyz(clr%xyz2hsv(xyz))

write (*,*) 'xyz2rgb : ', rgb2, sqrt( sum( (xyz-rgb2)**2 ) )
write (*,*) 'xyz2lab : ', lab2, sqrt( sum( (xyz-lab2)**2 ) )
write (*,*) 'xyz2luv : ', luv2, sqrt( sum( (xyz-luv2)**2 ) )
write (*,*) 'xyz2hsl : ', hsl2, sqrt( sum( (xyz-hsl2)**2 ) )
write (*,*) 'xyz2hsv : ', hsv2, sqrt( sum( (xyz-hsv2)**2 ) )

lab = clr%rgb2lab(rgb)
write (*,*) '---------'
write (*,*) 'target  : ', lab
xyz2 = clr%xyz2lab(clr%lab2xyz(lab))
rgb2 = clr%rgb2lab(clr%lab2rgb(lab))
luv2 = clr%luv2lab(clr%lab2luv(lab))
hsl2 = clr%hsl2lab(clr%lab2hsl(lab))
hsv2 = clr%hsv2lab(clr%lab2hsv(lab))

write (*,*) 'lab2xyz : ', xyz2, sqrt( sum( (lab-xyz2)**2 ) )
write (*,*) 'lab2rgb : ', rgb2, sqrt( sum( (lab-rgb2)**2 ) )
write (*,*) 'lab2luv : ', luv2, sqrt( sum( (lab-luv2)**2 ) )
write (*,*) 'lab2hsl : ', hsl2, sqrt( sum( (lab-hsl2)**2 ) )
write (*,*) 'lab2hsv : ', hsv2, sqrt( sum( (lab-hsv2)**2 ) )

luv = clr%rgb2luv(rgb)
write (*,*) '---------'
write (*,*) 'target  : ', luv
xyz2 = clr%xyz2luv(clr%luv2xyz(luv))
rgb2 = clr%rgb2luv(clr%luv2rgb(luv))
lab2 = clr%lab2luv(clr%luv2lab(luv))
hsl2 = clr%hsl2luv(clr%luv2hsl(luv))
hsv2 = clr%hsv2luv(clr%luv2hsv(luv))

write (*,*) 'luv2xyz : ', xyz2, sqrt( sum( (luv-xyz2)**2 ) )
write (*,*) 'luv2rgb : ', rgb2, sqrt( sum( (luv-rgb2)**2 ) )
write (*,*) 'luv2lab : ', lab2, sqrt( sum( (luv-lab2)**2 ) )
write (*,*) 'luv2hsl : ', hsl2, sqrt( sum( (luv-hsl2)**2 ) )
write (*,*) 'luv2hsv : ', hsv2, sqrt( sum( (luv-hsv2)**2 ) )

hsl = clr%rgb2hsl(rgb)
write (*,*) '---------'
write (*,*) 'target  : ', hsl
xyz2 = clr%xyz2hsl(clr%hsl2xyz(hsl))
rgb2 = clr%rgb2hsl(clr%hsl2rgb(hsl))
lab2 = clr%lab2hsl(clr%hsl2lab(hsl))
luv2 = clr%luv2hsl(clr%hsl2luv(hsl))
hsv2 = clr%hsv2hsl(clr%hsl2hsv(hsl))

write (*,*) 'hsl2xyz : ', xyz2, sqrt( sum( (hsl-xyz2)**2 ) )
write (*,*) 'hsl2rgb : ', rgb2, sqrt( sum( (hsl-rgb2)**2 ) )
write (*,*) 'hsl2lab : ', lab2, sqrt( sum( (hsl-lab2)**2 ) )
write (*,*) 'hsl2luv : ', luv2, sqrt( sum( (hsl-luv2)**2 ) )
write (*,*) 'hsl2hsv : ', hsv2, sqrt( sum( (hsl-hsv2)**2 ) )

hsv = clr%rgb2hsv(rgb)
write (*,*) '---------'
write (*,*) 'target  : ', hsv
xyz2 = clr%xyz2hsv(clr%hsv2xyz(hsv))
rgb2 = clr%rgb2hsv(clr%hsv2rgb(hsv))
lab2 = clr%lab2hsv(clr%hsv2lab(hsv))
luv2 = clr%luv2hsv(clr%hsv2luv(hsv))
hsl2 = clr%hsl2hsv(clr%hsv2hsl(hsv))

write (*,*) 'hsv2xyz : ', xyz2, sqrt( sum( (hsv-xyz2)**2 ) )
write (*,*) 'hsv2rgb : ', rgb2, sqrt( sum( (hsv-rgb2)**2 ) )
write (*,*) 'hsv2lab : ', lab2, sqrt( sum( (hsv-lab2)**2 ) )
write (*,*) 'hsv2luv : ', luv2, sqrt( sum( (hsv-luv2)**2 ) )
write (*,*) 'hsv2hsl : ', hsl2, sqrt( sum( (hsv-hsl2)**2 ) )

end program colortest