program colortest

use mod_EMsoft
use mod_kinds
use mod_global
use mod_colorspace 
use mod_image
use mod_io
use ISO_C_BINDING
use mod_io
use, intrinsic :: iso_fortran_env

IMPLICIT NONE 

type(colorspace_T)      :: clr 
real(kind=dbl)          :: xyz(0:2), rgb(0:2), lab(0:2), luv(0:2), hsv(0:2), hsl(0:2), a, p
real(kind=dbl)          :: xyz2(0:2), rgb2(0:2), lab2(0:2), luv2(0:2), hsv2(0:2), hsl2(0:2)

character(fnlen)        :: progname = 'colortest.f90'
character(fnlen)        :: progdesc = 'Standalone test program to generate a stereographic color projection '

type(IO_T)                                  :: Message
type(EMsoft_T)                              :: EMsoft

character(fnlen)                            :: fname, TIFF_filename 
real(kind=dbl),allocatable                  :: IPFmap(:,:,:)
integer(kind=irg)                           :: i, j, N, NN
real(kind=dbl)                              :: ix, iy, R, RR, X, Y, Z, phi, theta, mina, maxa, minp, maxp 

! declare variables for use in object oriented image module
integer                                     :: iostat, io_int(2)
character(len=128)                          :: iomsg
logical                                     :: isInteger, OPC
type(image_t)                               :: im
integer(int8), allocatable                  :: TIFF_image(:,:)
integer                                     :: dim2(2), Pm
integer(c_int32_t)                          :: result


EMsoft = EMsoft_T( progname, progdesc )


! initialize the color space class 
clr = colorspace_T( Nfold = 4 )


! write (*,*) '================'
write (*,*) 'Test of perceptually uniform color scheme code ' 

p = 0.05D0
do i = 0,11 
    a = ( dble(i) * cPi / 6.D0 ) / (2.D0 * cPi)
    luv = clr%sphere2rgb(a, p, w0=.TRUE.)
    ! write (*,*)  a, p, luv
end do 

N = 501
NN = (N-1)/2
R = dble(NN)
allocate( IPFmap(3, -NN:NN, -NN:NN) )

mina = 100.D0 
maxa = -100.D0 
minp = 100.D0 
maxp = -100.D0 

do i=-NN,NN
    ix = dble(i)/R
    do j=-NN,NN
      iy = dble(j)/R
      RR = ix*ix+iy*iy
      if (RR.le.1.D0) then ! convert to azimuth and polar angles on [0,1]
! coordinates on the sphere
        X = 2.D0*ix/(1.D0+RR)
        Y = -2.D0*iy/(1.D0+RR)
        Z = -(-1.D0+RR)/(1.D0+RR)
        a = acos(Z)/cPi 
        p = atan2(Y,X)
        if (p.lt.0.D0) p = p + 2.D0*cPi
        p = p/(2.D0*cPi)
        if (a.lt.mina) mina = a
        if (a.gt.maxa) maxa = a
        if (p.lt.minp) minp = p
        if (p.gt.maxp) maxp = p
        rgb = clr%sphere2rgb(a, p, w0=.TRUE.)
        IPFmap(1:3,i,j) = rgb(0:2)*255.D0
        ! IPFmap(1:3,i,j) = (/  a, p, 0.D0 /)*255.D0
      end if 
    end do 
end do 

write (*,*) ' range of azimuth/polar : ', mina, maxa, minp, maxp 

! finally, store the IPFmap in a color TIFF file.
fname = trim(EMsoft%generateFilePath('EMdatapathname'))//'../EMsoftOObuild/PUC_NH.tiff'
TIFF_filename = trim(fname)

! allocate memory for a color image; each pixel has 3 bytes [RGB]
allocate(TIFF_image(3*N,N))
TIFF_image = reshape( IPFmap, (/ 3*N, N /) )

! set up the image_t structure
im = image_t(TIFF_image)
im%dims = (/ N, N /)
im%samplesPerPixel = 3
if(im%empty()) call Message%printMessage("EMdpmerge: failed to convert array to rgb image")

! create the file
call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage(" Failed to write image to file : "//iomsg)
else
  call Message%printMessage(' IPF map written to '//trim(TIFF_filename),"(/A)")
end if

end program colortest