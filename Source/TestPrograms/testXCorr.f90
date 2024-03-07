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
program testXCorr

use mod_kinds
use mod_EMsoft
use mod_DSHT 
use mod_Wigner 
use mod_SHcorrelator 
use mod_io 
use mod_rotations 
use mod_global
use mod_fft_wrap 

use mod_SphereIndexer 

IMPLICIT NONE 

integer(kind=irg)                 :: dim, maxGeneric, maxLegendre
integer(kind=irg)                 :: maxL 
logical                           :: time = .TRUE.
character(fnlen)                  :: layout = 'legendre'
integer                           :: u

type(EMsoft_T)                    :: EMsoft
type(q_T)                         :: quat
type(e_T)                         :: euler

complex(kind=dbl),allocatable     :: alm(:,:), almRot(:,:)
real(kind=dbl),allocatable        :: mLPNH   (:,:), mLPSH   (:,:)
real(kind=dbl),allocatable        :: mLPNHRot(:,:), mLPSHRot(:,:)
real(kind=dbl)                    :: qu(4), qr(4), eu(0:2)
type(DiscreteSHT)                 :: SHTC 
type(SphereCorrelator)            :: s2corr 
integer(kind=irg)                 :: coeffCount, i, m, l, lmax, mmax, d, limL
real(kind=dbl)                    :: maxErr, rmsErr, relErr, delta, rDelta, xc
character(fnlen)                  :: fname, progname, progdesc 
logical                           :: fMr = .FALSE.
logical                           :: ref = .TRUE.
integer(kind=irg)                 :: fNf = 1
real(kind=dbl)                    :: eps = 1.0D-2, tmp

progname = 'testXCorr'
progdesc = 'Test program for spherical cross correlation'
EMsoft = EMsoft_T(progname, progdesc, silent=.TRUE.)

call FFTWisdom%load(EMsoft)

call setRotationPrecision('d')

! define size etc
d    = 32
dim  = 2*d+1
maxL = dim-2
limL = 0

! generate a random quaternion
call random_number(qu)
qu(:) = qu(:) / sqrt(dot_product(qu, qu))

! create a random signal on the sphere
write (*,*) ' ---> generating random data'
allocate(mLPNH(-d:d,-d:d),mLPSH(-d:d,-d:d))
call random_number(mLPNH)
call random_number(mLPSH)

! make the equator the same for both hemispheres
mLPSH(-d, :) = mLPNH(-d, :)
mLPSH( d, :) = mLPNH( d, :)
mLPSH( :,-d) = mLPNH( :,-d)
mLPSH( :, d) = mLPNH( :, d)

! save applied rotation + spherical function
open(newunit = u, file = 'input.raw', access = 'stream')
write(u) qu
write(u) mLPNH
write(u) mLPSH
close(u)

! build SHT calculator
write (*,*) ' ---> building dicreteSHTConstructor'
call SHTC%init(d, maxL, layout)

! compute the SHT of the random signal
allocate(alm(0:maxL-1,0:maxL-1))
call SHTC%analyze(mLPNH, mLPSH, alm)

! rotate spectra
write (*,*) ' ---> rotating spectra'
allocate(almRot(0:maxL-1,0:maxL-1))
call Wigner_rotateHarmonics(maxL, maxL, alm, almRot, qu)

! compute maximum cross correlation rotation
write (*,*) ' ---> constructing spherical correlator'
call s2corr%init(maxL)
write (*,*) ' ---> computing correlation and finding maximum'
xc = s2corr%correlate(alm, almRot, fMr, fNf, eu, ref, eps)

! convert peak rotation from zyz euler angle to quaternion
eu(0) = eu(0) + cPi / 2
eu(2) = eu(2) - cPi / 2
euler = e_T( edinp = eu )
quat = euler%eq()
qr = quat%q_copyd()

! print out applied + registered orientations
write(*,*) qu
write(*,*) qr, xc ! + cross correlation at peak

! print out error in registration
write(*,*) acos(min(sum(qu*qr), 1.D0)) * 180.0 / cPi, "degree error"

call FFTWisdom%save(EMsoft)

end program testXCorr
