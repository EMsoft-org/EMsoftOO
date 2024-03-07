!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!*                                                                     *
!* Copyright (c) 2019-2024, De Graef Group, Carnegie Mellon University *
!* All rights reserved.                                                *
!*                                                                     *
!* Author: William C. Lenthe                                           *
!*                                                                     *
!* EMSphInx is available for academic or non-profit non-commercial     *
!* research use. Please, see the license.txt file in this distribution *
!* for further details.                                                *
!*                                                                     *
!* Interested in a commercial license? Contact:                        *
!*                                                                     *
!* Center for Technology Transfer and Enterprise Creation              *
!* 4615 Forbes Avenue, Suite 302                                       *
!* Pittsburgh, PA 15213                                                *
!*                                                                     *
!* phone. : 412.268.7393                                               *
!* email  : innovation@cmu.edu                                         *
!* website: https://www.cmu.edu/cttec/                                 *
!*                                                                     *
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
program testSHT

use mod_kinds
use mod_DSHT 
use mod_io 
use mod_fft_wrap 
use mod_EMsoft

! use Ylm
! use rng

IMPLICIT NONE 

integer(kind=irg)                 :: dim, maxGeneric, maxLegendre
integer(kind=irg)                 :: maxL 
character(fnlen)                  :: layout = 'legendre'
logical                           :: time = .TRUE.

type(EMsoft_T)                    :: EMsoft

complex(kind=dbl),allocatable     :: refSpectra(:,:), alm(:,:)
real(kind=dbl),allocatable        :: mLPNH(:,:), mLPSH(:,:)
type(DiscreteSHT)                 :: SHTC 
integer(kind=irg)                 :: coeffCount, i, m, l, lmax, mmax, d, limL
real(kind=dbl)                    :: maxErr, rmsErr, relErr, delta, rDelta
character(fnlen)                  :: fname, progname, progdesc

progname = 'testSHT'
progdesc = 'Test program for spherical harmonic transform'
EMsoft = EMsoft_T(progname, progdesc, silent=.TRUE.)

d = 13
call FFTWisdom%load(EMsoft)
dim = 2*d+1
maxGeneric = d/2
maxLegendre = d-2
limL = 0

! set maxL depending on layout
if (trim(layout).eq.'legendre') then
  maxL = maxLegendre
else
  maxL = maxGeneric
end if

! allocate all memory
allocate(refSpectra(0:maxL-1,0:maxL-1))
refSpectra = cmplx(0.D0, 0.D0)

allocate(mLPNH(-d:d,-d:d),mLPSH(-d:d,-d:d))
mLPNH = 0.D0
mLPSH = 0.D0

allocate(alm(0:maxL-1,0:maxL-1))
alm = cmplx(0.D0, 0.D0)

! build SHT calculator and select linear combination of spherical harmonics to construct signal from
call SHTC%init(d, maxL, layout)

refSpectra = cmplx(1.D0,2.D0)
refSpectra(0:maxL-1,0) = cmplx(3.D0,0.D0)

! round trip calculation (spectra -> real -> spectra)
! compute the signal
call SHTC%synthesize(refSpectra, mLPNH, mLPSH, limL)

fname = 'testmLPNH.txt'
call SH_printRow(d, 10, mLPNH, fname) 
fname = 'testmLPSH.txt'
call SH_printRow(d, 10, mLPSH, fname) 


! compute the coefficients of the signal
call SHTC%analyze(mLPNH, mLPSH, alm)

! compute max absolute, relative, and rms errors in coefficients
coeffCount = 0
lmax = 0
mmax = 0
maxErr = 0.D0
rmsErr = 0.D0
relErr = 0.D0
do m=0,maxL-1
  do l=m,maxL-1
     write (*,*) m, l, alm(l,m), refSpectra(l,m)
     delta = abs( alm(l,m) - refSpectra(l,m))
     rDelta = delta / abs(refSpectra(l,m))
     if (delta.gt.maxErr) then
      maxErr = delta
      lmax = l
      mmax = m
    end if
    if (rDelta.gt.relErr) relErr = rDelta
    rmsErr = rmsErr + delta*delta
    coeffCount = coeffCount + 1
  end do 
end do
rmsErr = sqrt(rmsErr / dble(coeffCount))
write (*,"(' MaxErr = ',F12.8,'; relErr (%) = ',F12.8,'; rmsErr = ',F12.8)") maxErr, relErr*100.D0, rmsErr
write (*,"(' mode with the highest absolute error (l,m) : ',I4,I4)") lmax, mmax 
write (*,*) alm(lmax,mmax)

call FFTWisdom%save(EMsoft)

end program testSHT
