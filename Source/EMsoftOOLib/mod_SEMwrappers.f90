! ###################################################################
! Copyright (c) 2013-2020, Marc De Graef Research Group/Carnegie Mellon University
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
! EMsoft:mod_SEMwrappers.f90
!--------------------------------------------------------------------------
!
! MODULE: mod_SEMwrappers
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief routines that can be called by external code; all routines requiring HDF are in EMdymodHDF.f90
!
!> @date  01/21/18 MDG 1.0 separated C/C++ callable SEM routines out from original EMdymod module
!--------------------------------------------------------------------------
!
! general information: the ipar and fpar arrays for all the routines that are C-callable
! are identical, so we document here their component definitions; to allow for future expansion, each
! array has 40 entries, of which about half are currently (April 2016) used.
!
! integer(kind=irg) :: ipar(40)  components 
! ipar(1) : nx  = (numsx-1)/2
! ipar(2) : globalworkgrpsz
! ipar(3) : num_el
! ipar(4) : totnum_el
! ipar(5) : multiplier
! ipar(6) : devid
! ipar(7) : platid
! ipar(8) : CrystalSystem
! ipar(9) : Natomtypes
! ipar(10): SpaceGroupNumber
! ipar(11): SpaceGroupSetting
! ipar(12): numEbins
! ipar(13): numzbins
! ipar(14): mcmode  ( 1 = 'full', 2 = 'bse1' )
! ipar(15): numangle
! ipar(16): nxten = nx/10
! the following are only used in the master routine
! ipar(17): npx
! ipar(18): nthreads
! the following are only used in the EBSD pattern routine
! ipar(19): numx of detector pixels
! ipar(20): numy of detector pixels
! ipar(21): number of orientation in quaternion set
! ipar(22): binning factor (0-3)
! ipar(23): binned x-dimension
! ipar(24): binned y-dimension
! ipar(25): anglemode  (0 for quaternions, 1 for Euler angles)
! ipar(26): already initialized 
! ipar(27:40) : 0 (unused for now)

! real(kind=dbl) :: fpar(40)  components
! fpar(1) : sig
! fpar(2) : omega
! fpar(3) : EkeV
! fpar(4) : Ehistmin
! fpar(5) : Ebinsize
! fpar(6) : depthmax
! fpar(7) : depthstep
! fpar(8) : sigstart
! fpar(9) : sigend
! fpar(10): sigstep
! parameters only used in the master pattern routine
! fpar(11) : dmin
! fpar(12) : Bethe  c1
! fpar(13) : Bethe  c2
! fpar(14) : Bethe  c3
! parameters only used in the EBSD pattern routine
! fpar(15): pattern center x
! fpar(16): pattern center y
! fpar(17): scintillator pixel size
! fpar(18): detector tilt angle
! fpar(19): sample-scintillator distance
! fpar(20): beam current
! fpar(21): dwelltime
! fpar(22): gamma value
! fpar(23:40): 0 (unused for now)

! newly added in version 3.2, to facilitate passing EMsoft configuration
! strings back and forth to C/C++ programs that call EMdymod routines...
! character(fnlen)  :: spar(40)   configuration string components
! spar(1): EMsoftpathname
! spar(2): EMXtalFolderpathname
! spar(3): EMdatapathname
! spar(4): EMtmppathname
! spar(5): EMsoftLibraryLocation
! spar(6): EMSlackWebHookURL
! spar(7): EMSlackChannel
! spar(8): UserName
! spar(9): UserLocation
! spar(10): UserEmail
! spar(11): EMNotify
! spar(12): Develop
! spar(13): Release
! spar(14): h5copypath
! spar(15): EMsoftplatform
! spar(16): EMsofttestpath
! spar(17): EMsoftTestingPath
! spar(18): EMsoftversion
! spar(19): Configpath
! spar(20): Templatepathname
! spar(21): Resourcepathname
! spar(22): Homepathname
! spar(23): OpenCLpathname
! spar(24): Templatecodefilename
! spar(25): WyckoffPositionsfilename
! spar(26): Randomseedfilename
! spar(27): EMsoftnativedelimiter
! spar(28:40): '' (unused for now)


!
module mod_SEMwrappers

    !--------------------------------------------------------------------------
    ! Callback routine(s) to communicate progress with DREAM.3D package
    
    ! Define interface of call-back routine
    ! arguments are:
    !  objAddress: unique 8-byte integer to identify the calling class in DREAM.3D
    !  patternCompleted: integer indicating the current pattern ID number
    !
    ABSTRACT INTERFACE
       SUBROUTINE ProgressCallBack(objAddress, patternCompleted) bind(C)
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(c_size_t),INTENT(IN), VALUE          :: objAddress
        INTEGER(KIND=4), INTENT(IN), VALUE           :: patternCompleted
       END SUBROUTINE ProgressCallBack
    END INTERFACE
    
    
    ! similar callback routine, with two integer arguments
    ABSTRACT INTERFACE
       SUBROUTINE ProgressCallBack2(objAddress, loopCompleted, totalLoops, bseYield) bind(C)
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(c_size_t),INTENT(IN), VALUE          :: objAddress
        INTEGER(KIND=4), INTENT(IN), VALUE           :: loopCompleted
        INTEGER(KIND=4), INTENT(IN), VALUE           :: totalLoops
        REAL(KIND=4),INTENT(IN), VALUE              :: bseYield
       END SUBROUTINE ProgressCallBack2
    END INTERFACE
    
    ! similar callback routine, with two integer arguments
    ABSTRACT INTERFACE
       SUBROUTINE ProgressCallBack3(objAddress, loopCompleted, totalLoops, EloopCompleted, totalEloops) bind(C)
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(c_size_t),INTENT(IN), VALUE          :: objAddress
        INTEGER(KIND=4), INTENT(IN), VALUE           :: loopCompleted
        INTEGER(KIND=4), INTENT(IN), VALUE           :: totalLoops
        INTEGER(KIND=4), INTENT(IN), VALUE           :: EloopCompleted
        INTEGER(KIND=4), INTENT(IN), VALUE           :: totalELoops
       END SUBROUTINE ProgressCallBack3
    END INTERFACE
    
    !--------------------------------------------------------------------------

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! the routines starting with EMsoftC are callable from C/C++
! programs and can handle progress callback and a cancel request.
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! SUBROUTINE:EMsoftCgetEBSDPatterns
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief This subroutine can be called by a C/C++ program as a standalone function to compute EBSD patterns
!
!> @details This subroutine provides a method to compute a series of EBSD patterns and
!> can be called from an external C/C++ program; the routine provides a callback mechanism to
!> update the calling program about computational progress, as well as a cancel option.
!> The routine is intended to be called form a C/C++ program, e.g., DREAM.3D.  This routine is a simplified version
!> of the core of the EMEBSD program. 
!>
!> @param ipar array with integer input parameters
!> @param fpar array with float input parameters
!> @param EBSDpattern output array
!> @param quats quaternion input array
!> @param accum_e array with Monte Carlo histogram
!> @param mLPNH Northern hemisphere master pattern
!> @param mLPSH Southern hemisphere master pattern
!> @param cproc pointer to a C-function for the callback process
!> @param objAddress unique integer identifying the calling class in DREAM.3D
!> @param cancel character defined by DREAM.3D; when not equal to NULL (i.e., char(0)), the computation should be halted
!
!> @date 10/16/15 MDG 1.0 original in EMsoft
!> @date 04/14/21 MDG 4.0 rewrite for EMsoftOO
!--------------------------------------------------------------------------
recursive subroutine EMsoftCgetEBSDPatterns(ipar, fpar, EBSDpattern, quats, accum_e, mLPNH, mLPSH, cproc, objAddress, cancel) &
           bind(c, name='EMsoftCgetEBSDPatterns')    ! this routine is callable from a C/C++ program
!DEC$ ATTRIBUTES DLLEXPORT :: EMsoftCgetEBSDPatterns

! the input parameters are all part of a ipar and fpar input arrays instead of the usual namelist structure to
! make this routine callable by external programs, such as DREAM.3D

! The following is the mapping for the ipar and fpar array components used in this routine:
!
! ipar(1)  = mcnsx
! ipar(9)  = numset
! ipar(12) = detnumEbins
! ipar(17) = mpnpx
! ipar(19) = detnumsx
! ipar(20) = detnumsy
! ipar(21) = numquats
! ipar(22) = binning
! ipar(23) = binned x-dimension
! ipar(24) = binned y-dimension
! ipar(25) = anglemode
! ipar(26) = already initialized

! fpar(1)  = enl%MCsig
! fpar(2)  = enl%omega
! fpar(15) = enl%xpc
! fpar(16) = enl%ypc
! fpar(17) = enl%delta
! fpar(18) = enl%thetac
! fpar(19) = enl%L
! fpar(20) = enl%beamcurrent
! fpar(21) = enl%dwelltime
! fpar(22) = gammavalue

use mod_kinds
use mod_EMsoft
use mod_global
use mod_Lambert
use mod_quaternions
use mod_rotations
use mod_memory
use,INTRINSIC :: ISO_C_BINDING

IMPLICIT NONE

integer(c_int32_t),INTENT(IN)           :: ipar(wraparraysize)
real(kind=sgl),INTENT(IN)               :: fpar(wraparraysize)
integer(c_int32_t),PARAMETER            :: nq=4
real(kind=sgl),INTENT(IN)               :: quats(nq,ipar(21))
integer(c_int32_t),INTENT(IN)           :: accum_e(ipar(12),-ipar(1):ipar(1),-ipar(1):ipar(1))
real(kind=sgl),INTENT(IN)               :: mLPNH(-ipar(17):ipar(17), -ipar(17):ipar(17), ipar(12), ipar(9))
real(kind=sgl),INTENT(IN)               :: mLPSH(-ipar(17):ipar(17), -ipar(17):ipar(17), ipar(12), ipar(9))
real(kind=sgl),INTENT(OUT)              :: EBSDpattern(ipar(23),ipar(24),ipar(21))
TYPE(C_FUNPTR), INTENT(IN), VALUE       :: cproc
integer(c_size_t),INTENT(IN), VALUE     :: objAddress
character(len=1),INTENT(IN)             :: cancel

type(memory_T)                          :: mem 
type(EMsoft_T)                          :: EMsoft
type(Quaternion_T)                      :: quat
type(e_T)                               :: eu
type(q_T)                               :: qu

! various variables and arrays
character(fnlen)                        :: pname, pdesc
real(kind=sgl)                          :: fullsizepattern(ipar(19),ipar(20)), binned(ipar(23),ipar(24))
real(kind=irg),allocatable,save         :: accum_e_detector(:,:,:)
real(kind=sgl),allocatable,save         :: rgx(:,:), rgy(:,:), rgz(:,:)
real(kind=sgl),allocatable,save         :: mLPNHsum(:,:,:), mLPSHsum(:,:,:)
real(kind=sgl),save                     :: prefactor
real(kind=sgl),allocatable              :: scin_x(:), scin_y(:)                 ! scintillator coordinate arrays [microns]
real(kind=sgl)                          :: alp, ca, sa, cw, sw
real(kind=sgl)                          :: L2, Ls, Lc     ! distances
integer(kind=irg)                       :: nix, niy, binx, biny,  nixp, niyp, i, j, Emin, Emax, istat, k, ip, dn, cn, & 
                                           ii, jj, binfac, ipx, ipy, epl      ! various parameters
real(kind=sgl)                          :: dc(3), scl, alpha, theta, gam, pcvec(3), dp, calpha           ! direction cosine array
real(kind=sgl)                          :: sx, dx, dxm, dy, dym, rhos, x, bindx         ! various parameters
real(kind=sgl)                          :: ixy(2)
real(kind=dbl),parameter                :: nAmpere = 6.241D+18 
PROCEDURE(ProgressCallBack), POINTER    :: proc

! link the proc procedure to the cproc argument
CALL C_F_PROCPOINTER (cproc, proc)

call setRotationPrecision('Single')

! binned pattern dimensions
binx = ipar(23)
biny = ipar(24)
binfac = 2**ipar(22)
bindx = 1.0/float(binfac)**2

if (ipar(26).eq.0) then 
!====================================
! ------ generate the detector rgx, rgy, rgz arrays (and a few others)
!====================================
  mem = memory_T()

  call mem%alloc(mLPNHsum, (/ ipar(17), ipar(17), ipar(12) /), 'mLPNHsum', 0.0, startdims = (/ -ipar(17), -ipar(17), 1/) )
  call mem%alloc(mLPSHsum, (/ ipar(17), ipar(17), ipar(12) /), 'mLPSHsum', 0.0, startdims = (/ -ipar(17), -ipar(17), 1/) )

  ! Stuart Wright: for some reason the following calls do not work on my Windows 10 computer, VS 2015
  ! so I unwrapped the code to perform this function explicitly [modified with platform check, MDG]
  pname = ''
  pdesc = ''
  EMsoft = EMsoft_T( pname, pdesc, silent=.TRUE.)
  if (trim(EMsoft%getConfigParameter('EMsoftplatform')).ne.'Windows') then 
     mLPNHsum = sum(mLPNH,4)
     mLPSHsum = sum(mLPSH,4)
  else
    do i=-ipar(17),ipar(17)
        do j=-ipar(17),ipar(17)
            do k=1,ipar(12)
                do ii=1,ipar(9)
                    mLPNHsum(i,j,k) = mLPNHsum(i,j,k) + mLPNH(i,j,k,ii)
                    mLPSHsum(i,j,k) = mLPSHsum(i,j,k) + mLPSH(i,j,k,ii)
                end do
            end do
        end do
    end do
  end if

  call mem%alloc(scin_x, (/ ipar(19) /), 'scin_x')
  call mem%alloc(scin_y, (/ ipar(20) /), 'scin_y') 

  scin_x = - ( -fpar(15) - ( 1.0 - float(ipar(19)) ) * 0.5 - (/ (i-1, i=1,ipar(19)) /) ) * fpar(17)
  scin_y = ( fpar(16) - ( 1.0 - float(ipar(20)) ) * 0.5 - (/ (i-1, i=1,ipar(20)) /) ) * fpar(17)

! auxiliary angle to rotate between reference frames
  alp = 0.5 * cPi - (fpar(1) - fpar(18)) * dtor
  ca = cos(alp)
  sa = sin(alp)

  cw = cos(fpar(2) * dtor)
  sw = sin(fpar(2) * dtor)

! compute auxilliary interpolation arrays
  call mem%alloc(rgx, (/ ipar(19),ipar(20) /), 'rgx')
  call mem%alloc(rgy, (/ ipar(19),ipar(20) /), 'rgy')
  call mem%alloc(rgz, (/ ipar(19),ipar(20) /), 'rgz')

  epl = ipar(20)+1
  L2 = fpar(19) * fpar(19)
  do j=1,ipar(19)
    sx = L2 + scin_x(j) * scin_x(j)
    Ls = -sw * scin_x(j) + fpar(19) * cw
    Lc = cw * scin_x(j) + fpar(19) * sw
    do i=1,ipar(20)
!   rhos = 1.0/sqrt(sx + scin_y(i)**2)
     rgx(j,epl-i) = (scin_y(i) * ca + sa * Ls) ! * rhos
     rgy(j,epl-i) = Lc ! * rhos
     rgz(j,epl-i) = (-sa * scin_y(i) + ca * Ls) ! * rhos
! make sure that these vectors are normalized !
     x = sqrt(rgx(j,epl-i)**2+rgy(j,epl-i)**2+rgz(j,epl-i)**2)
     rgx(j,epl-i) = rgx(j,epl-i) / x
     rgy(j,epl-i) = rgy(j,epl-i) / x
     rgz(j,epl-i) = rgz(j,epl-i) / x
    end do
  end do

! remove the auxiliary arrays scin_x and scin_y
  call mem%dealloc(scin_x, 'scin_x')
  call mem%dealloc(scin_y, 'scin_y')

!====================================
! ------ create the equivalent detector energy array
!====================================
! from the Monte Carlo energy data, we need to extract the relevant
! entries for the detector geometry defined above.  

! determine the scale factor for the Lambert interpolation; the square has
! an edge length of 2 x sqrt(pi/2)
  scl = float(ipar(1)) 

! energy summation will go over all energy bins
  Emin = 1
  Emax = ipar(12)

  call mem%alloc(accum_e_detector, (/ ipar(12),ipar(19),ipar(20) /), 'accum_e_detector')

! correction of change in effective pixel area compared to equal-area Lambert projection
  alpha = atan(fpar(17)/fpar(19)/sqrt(sngl(cPi)))
  ipx = ipar(19)/2 + nint(fpar(15))
  ipy = ipar(20)/2 + nint(fpar(16))
  if (ipx .gt. ipar(19)) ipx = ipar(19)
  if (ipx .lt. 1) ipx = 1
  if (ipy .gt. ipar(20)) ipy = ipar(20)
  if (ipy .lt. 1) ipy = 1
  pcvec = (/ rgx(ipx,ipy), rgy(ipx,ipy), rgz(ipx,ipy) /)
  calpha = cos(alpha)
! initialize the Lambert class 
  dc = (/ 0.0, 0.0, 0.0 /)
  do i=1,ipar(19)
    do j=1,ipar(20)
! do the coordinate transformation for this detector pixel
       dc = (/ rgx(i,j), rgy(i,j), rgz(i,j) /)
! make sure the third one is positive; if not, switch all 
       if (dc(3).lt.0.0) dc = -dc

! convert these direction cosines to coordinates in the Rosca-Lambert projection
       call LambertgetInterpolation(dc, scl, int(ipar(1)), int(ipar(1)), nix, niy, nixp, niyp, &
                                    dx, dy, dxm, dym, swap=.TRUE.)

! do the area correction for this detector pixel
        dp = dot_product(pcvec,dc)
        if ((i.eq.ipx).and.(j.eq.ipy)) then
          gam = 0.25 
        else
          theta = calpha*calpha + dp*dp - 1.0
          gam = theta**1.5/(calpha**3) * 0.25
          
        end if
! interpolate the intensity 
        do k= Emin, Emax
          accum_e_detector(k,i,j) = gam * (accum_e(k,nix,niy) * dxm * dym + &
                                           accum_e(k,nix+1,niy) * dx * dym + &
                                           accum_e(k,nix,niy+1) * dxm * dy + &
                                           accum_e(k,nix+1,niy+1) * dx * dy)
        end do
    end do
  end do 
  prefactor = 0.25D0 * nAmpere * fpar(20) * fpar(21)  * 1.0D-15 / sum(accum_e_detector)
  accum_e_detector = accum_e_detector * prefactor
end if  ! initialize detector arrays 

! from here on, we simply compute the EBSD patterns by interpolation, using the above arrays
! no intensity scaling or anything else...other than multiplication by pre-factor
! intensity scaling is left to the user of the calling program.

! define some parameters and initialize EBSDpattern
scl = float(ipar(17)) 
EBSDpattern = 0.0
fullsizepattern = 0.0
dn = nint(float(ipar(21))*0.01)
cn = dn

! here is the main loop over all quaternions
quat = Quaternion_T( q = (/ 0.0, 0.0, 0.0, 0.0 /) )
eu = e_T( )
qu = q_T( qinp = (/ 0.0, 0.0, 0.0, 0.0 /) )
quatloop: do ip=1,ipar(21)
  binned = 0.0
  fullsizepattern = 0.0
  if (ipar(25).eq.0) then 
    call quat%set_quats( quats(1:4,ip) )
  else
    call eu%e_set( quats(1:3,ip) )  ! this assumes that the input Euler angles are in radians
    qu = eu%eq()
    call quat%set_quats( qu%q_copy() ) 
  end if
  do i=1,ipar(19)
    do j=1,ipar(20)
! do the active coordinate transformation for this euler angle
      dc = quat%quat_Lp( (/ rgx(i,j), rgy(i,j), rgz(i,j) /) )
! normalize dc
      dc = dc/sqrt(sum(dc*dc))
! convert these direction cosines to coordinates in the Rosca-Lambert projection (always square projection !!!)
      call LambertgetInterpolation(dc, scl, int(ipar(17)), int(ipar(17)), nix, niy, nixp, niyp, dx, dy, dxm, dym)

      if (dc(3).gt.0.0) then ! we're in the Northern hemisphere
        do k=1,ipar(12) 
          fullsizepattern(i,j) = fullsizepattern(i,j) + accum_e_detector(k,i,j) * ( mLPNHsum(nix,niy,k) * dxm * dym +&
                                      mLPNHsum(nixp,niy,k) * dx * dym + mLPNHsum(nix,niyp,k) * dxm * dy + &
                                      mLPNHsum(nixp,niyp,k) * dx * dy )
        end do
      else                   ! we're in the Southern hemisphere
        do k=1,ipar(12) 
          fullsizepattern(i,j) = fullsizepattern(i,j) + accum_e_detector(k,i,j) * ( mLPSHsum(nix,niy,k) * dxm * dym +&
                                      mLPSHsum(nixp,niy,k) * dx * dym + mLPSHsum(nix,niyp,k) * dxm * dy + &
                                      mLPSHsum(nixp,niyp,k) * dx * dy )
        end do
      end if
    end do
  end do

! bin the pattern if necessary and apply the gamma scaling factor
  if (binx.ne.ipar(19)) then 
    do ii=1,ipar(19),binfac
        do jj=1,ipar(20),binfac
            binned(ii/binfac+1,jj/binfac+1) = &
            sum(fullsizepattern(ii:ii+binfac-1,jj:jj+binfac-1))
        end do
    end do
    EBSDpattern(1:binx,1:biny,ip) = (binned(1:binx,1:biny)* bindx)**fpar(22)
  else
    EBSDpattern(1:binx,1:biny,ip) = (fullsizepattern(1:binx,1:biny))**fpar(22)
  end if

! has the cancel flag been set by the calling program ?
  if(cancel.ne.char(0)) EXIT quatloop

! update the progress counter and report it to the calling program via the proc callback routine
  if(objAddress.ne.0) then
    if (ip.ge.cn) then
      cn = cn+dn
      call proc(objAddress, ip)
    end if
  end if

end do quatloop


end subroutine EMsoftCgetEBSDPatterns





    
end module mod_SEMwrappers
    