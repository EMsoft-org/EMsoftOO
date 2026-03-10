! ###################################################################
! Copyright (c) 2016-2024, Marc De Graef Research Group/Carnegie Mellon University
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

program EMplay
  !! author: MDG
  !! version: 1.0 
  !! date: 02/26/24
  !!
  !! test program 

use mod_kinds
use mod_global
use mod_EMsoft
use mod_rotations
use mod_quaternions
use mod_math
use mod_filters
use mod_povray
use mod_so3
use mod_io
use mod_symmetry
use mod_crystallography
use mod_DIC
use HDF5
use mod_HDFsupport
use mod_platformsupport
use bspline_module
use bspline_kinds_module, only: wp, ip

IMPLICIT NONE

character(fnlen)             :: progname = 'EMplay.f90'
character(fnlen)             :: progdesc = 'test program for HR-EBSD-DIC'

type(EMsoft_T)               :: EMsoft
type(DIC_T)                  :: DIC
type(q_T)                    :: qu 
type(o_T)                    :: om 
type(e_T)                    :: eu
type(h_T)                    :: ho
type(c_T)                    :: cu
type(QuaternionArray_T)      :: qAR, qsym
type(Quaternion_T)           :: q
type(Orientation_T)          :: ori
type(IO_T)                   :: Message

real(kind=wp), allocatable   :: refpat(:,:), defpat(:,:)
real(kind=sgl),allocatable   :: tmppat(:,:)
real(kind=wp)                :: DD, PCx, PCy, val, err, errmax, rnxi, rnyi, hg(8), W(3,3), gradCIC(8), Hessian(8,8), &
                                minx, miny, xi1max, xi2max, normdp, oldnorm, oldW(3,3), horiginal(8), CIC, sol(8), &
                                homographies(8,1000), hpartial(8), scalingfactor
real(kind=dbl)               :: Wnew(3,3), Winv(3,3), dx, dy, p2(3), Woriginal(3,3), alp, srt(3,3), srtrot(3,3), ac
real(kind=dbl)               :: qqq(4,-10:10), hhh(3,-10:10), sss(3), ccc(3), aaa(4), omval, dp, somval
integer(kind=irg)            :: nx, ny, nxy, nbx, nby, i, ii, j, NSR, cnt, nxSR, nySR, jj, recordsize, ierr, maxnumit  
real(wp)                     :: tol
integer(kind=4)              :: hnstat
character(fnlen)             :: fname, gname, hostname 

call setRotationPrecision('d')

alp = (3.D0*cPi/4.D0)**(1.D0/3.D0)
ac = cPi**(2.D0/3.D0)*0.5D0
nx = 10

do i=-nx,nx
  cu = c_T( cdinp = (/ 0.D0, 0.D0, dble(i)*ac/dble(nx) /) )
  ho = cu%ch()
  hhh(1:3,i) = ho%h_copyd()
  qu = cu%cq()
  qqq(1:4,i) = qu%q_copyd()
end do 

somval = 0.D0
do i=-nx,nx-1
  dp = sum(qqq(1:4,i)*qqq(1:4,i+1))
  omval = 2.D0*acos(dp)/dtor 
  somval = somval + omval
  write (*,"(I3,' & ', F7.4,' & ', F7.4,' & ', F7.4,' & ', F7.4,' & ', F7.4,' \\')") i, dble(i)*ac/dble(nx), &
        hhh(3,i), qqq(1,i), dp, omval
end do 

write (*,*) ' sum of angles : ', somval


! hnstat = system_hostnm(hostname)

! EMsoft = EMsoft_T( progname, progdesc )

! ! this program is a first conversion of the IDL HR-EBSD-DIC code into f90; we'll test
! ! everything here before turning it into a fullblown EMsoftOO program. We will use the IDL 
! ! convention of starting arrays at coordinate 0 for ease of comparison.  Also, keep in mind
! ! that fortran arrays are column-based whereas IDL is row-based...

! ! read the reference pattern from a binary file; the patterns have been 
! ! pre-processed and scaled between 0.D0 and 1.D0
! nx = 640
! ny = 480

! ! instantiate the DIC class
! ! this also initializes the x and y coordinate arrays
! ! DIC = DIC_T( nx, ny, normalize = .FALSE.) !, cross = (/ 310, 330, 230, 250 /) )
! DIC = DIC_T( nx, ny, normalize = .TRUE., cross = (/ 0, 0, 0, 0 /) )
! call DIC%setverbose(.FALSE.)
! ! call DIC%setverbose(.TRUE.)

! nxy = nx * ny
! recordsize = nxy * 4

! allocate(tmppat(nx,ny), refpat(0:nx-1,0:ny-1), defpat(0:nx-1,0:ny-1))
! ! open(20,file='/Users/mdg/playarea/DIC/verification-sq.data',status='old',form='unformatted')
! if (trim(hostname).eq.'Mac-Studio.fios-router.home') then 
!     fname = EMsoft%generateFilePath('EMdatapathname','DIC/refpat.data')
! else
!     fname = EMsoft%generateFilePath('EMdatapathname','playarea/DIC/refpat.data')
! end if
! open(20,file=trim(fname),status='old',form='unformatted')
! read(20) tmppat
! close(20,status='keep')
! ! refpat = dble( loghipass(tmppat, 5, nx, ny) ) 
! refpat = dble( tmppat ) 
! refpat = refpat - minval(refpat)
! refpat = refpat/maxval(refpat)
! call DIC%setpattern('r', refpat)

! write (*,*) ' range refpat ', minval(refpat), maxval(refpat)
! write (*,*) refpat(1:3,1:3)

! ! define some parameters for now...
! tol = 100 * epsilon(1.0_wp)
! PCx = 0.5_wp ! real(nx/2-1,wp)
! PCy = 0.5_wp ! real(ny/2-1,wp)
! ! PCx = real(nx/2-1,wp)
! ! PCy = real(ny/2-1,wp)
! DD = 15000.0_wp/50.0_wp 

! ! define the border widths nbx and nby for the subregion
! nbx = 30
! nby = 30
! call DIC%defineSR( nbx, nby, PCx, PCy)

! call DIC%getbsplines(refp=.TRUE., verify=.TRUE., grads=.TRUE.)

! !zero-mean and normalize the referenceSR and targetSR arrays
! call DIC%applyZMN(doreference=.TRUE.)

! ! compute the Hessian matrix
! call DIC%getHessian( Hessian )

! write (*,*) 'Hessian'
! do i=1,8
!     write (*,*) Hessian(i,1:8)
! end do

! if (trim(hostname).eq.'Mac-Studio.fios-router.home') then 
!     gname = EMsoft%generateFilePath('EMdatapathname','DIC/test/homographypatterns2.data')
! else
!     gname = EMsoft%generateFilePath('EMdatapathname','playarea/DIC/test/homographypatterns2.data')
! end if
! open(unit=dataunit,file=trim(gname),&
!      status='unknown',form='unformatted',access='direct',recl=recordsize,iostat=ierr)

! open(unit=28,file='output2-verify.txt',status='unknown',form='formatted')

! horiginal = (/ (0.0_wp, i=1,8) /)
! call DIC%applyHomography(horiginal, PCx, PCy)

! maxnumit = 50 

! do jj=1, 1000
!     ! call Message%printMessage(' ---------------------- ')
!     if (mod(jj,100).eq.0) write (*,*) 'starting pattern ', jj

! read(dataunit,rec=jj, iostat=ierr) tmppat
! ! defpat = dble( loghipass(tmppat, 5, nx, ny) )
! defpat = dble( tmppat)
! defpat = defpat - minval(defpat)
! defpat = defpat/maxval(defpat)
! call DIC%setpattern('d', defpat)

! Woriginal = DIC%getShapeFunction(horiginal)

! ! get the derivatives of the reference pattern to determine the Hessian
! call DIC%getbsplines(defp=.TRUE.)

! ! ! zero-mean and normalize the referenceSR and targetSR arrays
! call DIC%applyZMN(dotarget=.TRUE.)

! ! ! initial CIC function (to be minimized ...)
! call DIC%getresiduals( CIC )
! ! stop

! ! initialize the homography coefficients to zero
! ! in a more general implementation we might initialize them to the 
! ! values for one of the nearest neighbor patterns

! ! simple multiplicative term for the partial solutions
! ! maybe this will speed up convergence a little...
! scalingfactor = 1.45D0  

! ! and here we start the loop 
! do ii=1, maxnumit
!     ! write (*,*) ' iteration # ',ii
!     if (ii.eq.1) then  ! initialize to identity homography in first cycle
!         hpartial = (/ (0.0_wp, i=1,8) /)
!         W = DIC%getShapeFunction(hpartial)
!     end if

!     call DIC%applyHomography(hpartial, PCx, PCy, dotarget=.TRUE.)

!     ! zero-mean and normalize the deformed sub-region array wtarget
!     call DIC%applyZMN( dotarget = .TRUE. )

!     ! residuals
!     call DIC%getresiduals( CIC )

!     ! gradCIC and solve equation
!     call DIC%solveHessian(SOL, normdp)

!     ! convert to a shape function and take the inverse
!     Wnew = DIC%getShapeFunction(reshape(SOL, (/8/))*scalingfactor)
!     Winv = matrixInvert_wp( Wnew )
!     W = matmul( W, Winv )
!     W = W / W(3,3)
!     hg = DIC%getHomography(W)
!     hpartial = reshape(SOL, (/8/)) * scalingfactor

!     ! if (normdp.lt.0.0005D0) then
!     if (normdp.lt.0.000025D0) then
!     ! if (normdp.lt.0.001D0) then
!         EXIT
!     end if
!     ! call Message%printMessage(' ---------------------- ')
! end do 

! ! W = matrixInvert_wp( W )
! ! W = W / W(3,3)
! hg = DIC%getHomography(W)
! ! write (*,*) '' 
! ! write (*,*) ' Final homography : '
! ! write (*,*) DIC%getHomography(W)
! ! write (*,*) '' 
! ! write (*,*) ' Target homography : '
! ! write (*,*) horiginal
! ! write (*,*) '' 
! ! write (*,*) ' Differences : '
! ! write (*,*) horiginal-hg

! ! write results to data file (single precision because IDL has a bug for double precision)
! if (ii.eq.maxnumit+1) then 
!     write (unit=28,FMT='(9(F12.8,","),I4)') (/ (0.0, i=1,8) /), real(normdp), ii
! else
!     write (unit=28,FMT='(9(F12.8,","),I4)') real(hg), real(normdp), ii
! end if

! ! if (jj.eq.3) stop
!   call DIC%cleanup()

! ! end if 

! end do ! loop over jj 

! close(28,status='keep')
! ! close(20,status='keep')
! close(dataunit,status='keep')

end program EMplay
