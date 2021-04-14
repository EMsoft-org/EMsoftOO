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
quat = Quaternion_T( q = (/ 1.0, 0.0, 0.0, 0.0 /) )
eu = e_T( )
qu = q_T( qinp = (/ 1.0, 0.0, 0.0, 0.0 /) )
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

! bin the pattern if necessary and apply the gamma scaling factor (except when calling from IDL)
  if (binx.ne.ipar(19)) then 
    do ii=1,ipar(19),binfac
        do jj=1,ipar(20),binfac
            binned(ii/binfac+1,jj/binfac+1) = &
            sum(fullsizepattern(ii:ii+binfac-1,jj:jj+binfac-1))
        end do
    end do
    if (fpar(22).gt.0.0) then 
      EBSDpattern(1:binx,1:biny,ip) = (binned(1:binx,1:biny) * bindx)**fpar(22)
    else
      EBSDpattern(1:binx,1:biny,ip) = binned(1:binx,1:biny) * bindx
    end if
  else
    if (fpar(22).gt.0.0) then 
      EBSDpattern(1:binx,1:biny,ip) = (fullsizepattern(1:binx,1:biny))**fpar(22)
    else
      EBSDpattern(1:binx,1:biny,ip) = fullsizepattern(1:binx,1:biny)
    end if 
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

!--------------------------------------------------------------------------
!
! SUBROUTINE:EMsoftCgetECPatterns
!
!> @author Saransh Singh/Marc De Graef, Carnegie Mellon University
!
!> @brief This subroutine can be called by a C/C++ program as a standalone function to compute ECPs
!
!> @details This subroutine provides a method to compute a series of ECPs and
!> can be called from an external C/C++ program; the routine provides a callback mechanism to
!> update the calling program about computational progress, as well as a cancel option.
!> The routine is intended to be called form a C/C++ program, e.g., DREAM.3D.  This routine is a simplified version
!> of the core of the EMECP program. 
!>
!> This routine will first compute the incident cone vectors etc. if necessary, and then perform
!> the usual interpolation from the square Lambert projection. The pattern will be a basic pattern,
!> without any intensity scaling or binning etc; the calling program should take care of those 
!> operations.
!
!> @param ipar array with integer input parameters
!> @param fpar array with float input parameters
!> @param ECPattern output array
!> @param quats array of quaternions
!> @param accum_e array with Monte Carlo histogram
!> @param mLPNH Northern hemisphere master pattern
!> @param mLPSH Southern hemisphere master pattern
!
!> @date 10/16/15  SS 1.0 original
!> @date 11/02/14 MDG 1.1 put all integer parameters inside ipar and fixed size of ipar/fpar
!> @date 11/04/15 MDG 1.2 added array of quaternions as input parameter
!> @date 01/14/16 MDG 2.0 forked from original SingleECPattern routine; SAVE atrribute removed; ipar redefined (ipar(1) removed)
!> @date 04/15/21 MDG 3.0 rewrite for EMsoftOO
!--------------------------------------------------------------------------
recursive subroutine EMsoftCgetECPatterns(ipar, fpar, ECpattern, quats, accum_e, mLPNH, mLPSH, cproc, objAddress, cancel) &
           bind(c, name='EMsoftCgetECPatterns')    ! this routine is callable from a C/C++ program
!DEC$ ATTRIBUTES DLLEXPORT :: EMsoftCgetECPatterns

! the input parameters are all part of a ipar and fpar input arrays instead of the usual namelist structures.
! The following is the mapping:
!
! ipar(1) = detnumpix
! ipar(2) = numangle
! ipar(3) = mcnsx
! ipar(4) = numset
! ipar(5) = mpnpx
! ipar(6) = numquats

! fpar(1) = ecpnl%thetac
! fpar(2) = ecpnl%sampletilt
! fpar(3) = ecpnl%workingdistance
! fpar(4) = ecpnl%Rin
! fpar(5) = ecpnl%Rout
! fpar(6) = ecpnl%sigstart
! fpar(7) = ecpnl%sigend
! fpar(8) = ecpnl%sigstep

use mod_kinds
use mod_EMsoft
use mod_global
use mod_Lambert
use mod_quaternions
use mod_rotations
use mod_memory
use,INTRINSIC :: ISO_C_BINDING

IMPLICIT NONE

integer(c_int32_t),PARAMETER            :: nipar=6
integer(c_int32_t),PARAMETER            :: nfpar=8
integer(c_int32_t),PARAMETER            :: nq=4
integer(c_int32_t),INTENT(IN)           :: ipar(nipar)
real(kind=sgl),INTENT(IN)               :: fpar(nfpar)
real(kind=sgl),INTENT(OUT)              :: ECpattern(ipar(1),ipar(1),ipar(6))
real(kind=sgl),INTENT(IN)               :: quats(nq,ipar(6))
real(kind=sgl),INTENT(IN)               :: accum_e(ipar(2),-ipar(3):ipar(3),-ipar(3):ipar(3))
real(kind=sgl),INTENT(IN)               :: mLPNH(-ipar(5):ipar(5), -ipar(5):ipar(5), ipar(4))
real(kind=sgl),INTENT(IN)               :: mLPSH(-ipar(5):ipar(5), -ipar(5):ipar(5), ipar(4))
TYPE(C_FUNPTR), INTENT(IN), VALUE       :: cproc
integer(c_size_t),INTENT(IN), VALUE     :: objAddress
character(len=1),INTENT(IN)             :: cancel

type(memory_T)                          :: mem 
type(EMsoft_T)                          :: EMsoft
type(Quaternion_T)                      :: quat
type(q_T)                               :: qu

real(kind=sgl),allocatable              :: klist(:,:,:), rgx(:,:), rgy(:,:), rgz(:,:), weightfact(:)
real(kind=sgl),allocatable              :: mLPNHsum(:,:), mLPSHsum(:,:)

real(kind=sgl)                          :: kk(3), thetacr, ktmax, delta, wf
integer(kind=irg)                       :: istat, imin, imax, jmin, jmax, ii ,jj, nazimuth, npolar, nsig, ip, dn, cn
integer(kind=irg)                       :: ipolar, iazimuth, isig, isampletilt, nix, niy, nixp, niyp, isigp
real(kind=sgl)                          :: thetain, thetaout, polar, azimuthal, delpolar, delazimuth, om(3,3)
real(kind=sgl)                          :: dc(3), scl, deltheta, acc_sum, MCangle, ixy(2), dx, dy, dxm, dym, dp
PROCEDURE(ProgressCallBack), POINTER    :: proc

! link the proc procedure to the cproc argument
CALL C_F_PROCPOINTER (cproc, proc)

call setRotationPrecision('Single')

!==================================================================================
! ------ generate the detector klist, rgx, rgy, rgz, weightfactors arrays 
!==================================================================================

imin = 1
imax = ipar(1)
jmin = 1
jmax = ipar(1)

mem = memory_T() 
call mem%alloc(mLPNHsum, (/ ipar(5), ipar(5) /), 'mLPNHsum', startdims = (/ -ipar(5), -ipar(5) /) )
call mem%alloc(mLPSHsum, (/ ipar(5), ipar(5) /), 'mLPSHsum', startdims = (/ -ipar(5), -ipar(5) /) )
mLPNHsum = sum(mLPNH,3)
mLPSHsum = sum(mLPSH,3)

call mem%alloc(klist, (/ 3, ipar(1), ipar(1)/), 'klist')
kk = (/0.0,0.0,1.0/)
thetacr = DtoR*fpar(1)
ktmax = tan(thetacr)
delta = 2.0*ktmax/dble(ipar(1)-1)
 
do ii = imin, imax
    do jj = jmin, jmax
        klist(1:3,ii,jj) = (/-ktmax+delta*(ii-1),-ktmax+delta*(jj-1),0.0/) + kk(1:3)
        klist(1:3,ii,jj) =  klist(1:3,ii,jj)/sqrt(sum( klist(1:3,ii,jj)**2))
    end do
end do

thetain = atan2(fpar(4),fpar(3))
thetaout = atan2(fpar(5),fpar(3))

om(1,:) = (/cos(fpar(2)*sngl(dtor)),0.0,sin(fpar(2)*sngl(dtor))/)
om(2,:) = (/0.0,1.0,0.0/)
om(3,:) = (/-sin(fpar(2)*sngl(dtor)),0.0,cos(fpar(2)*sngl(dtor))/)

npolar = nint((thetaout - thetain)*180.0/cPi) + 1
delpolar = (thetaout - thetain)/float(npolar-1)

nazimuth = 361
delazimuth = 2.0*cPi/float(nazimuth-1)

call mem%alloc(rgx, (/ npolar, nazimuth /), 'rgx')
call mem%alloc(rgy, (/ npolar, nazimuth /), 'rgy')
call mem%alloc(rgz, (/ npolar, nazimuth /), 'rgz')

do ipolar = 1,npolar
     polar = thetain + float(ipolar-1)*delpolar

     do iazimuth = 1,nazimuth
         azimuthal = float(iazimuth-1)*delazimuth

         dc(1) = cos(azimuthal)*sin(polar)
         dc(2) = sin(azimuthal)*sin(polar)
         dc(3) = cos(polar)

         dc = matmul(om,dc)

         rgx(ipolar,iazimuth) = dc(1)
         rgy(ipolar,iazimuth) = dc(2)
         rgz(ipolar,iazimuth) = dc(3)
    end do
end do

!===================================================================
! ------ generate the weight factors from the monte carlo histogram
!===================================================================

scl = float(ipar(3))
nsig = nint(fpar(1) + abs(fpar(2))) + 1
deltheta = (fpar(1)+abs(fpar(2)))/float(nsig-1)

call mem%alloc(weightfact, (/ nsig /), 'weightfac', 0.0) 

do isig = 1,nsig
    acc_sum = 0.0
    MCangle = (isig - 1)*deltheta
    isampletilt = nint((MCangle - fpar(6))/fpar(8))

    if (isampletilt .lt. 1) then
        isampletilt = abs(isampletilt) + 1
    else
        isampletilt = isampletilt + 1
    end if

    do ipolar = 1,npolar
        do iazimuth = 1,nazimuth
            dc(1:3) = (/rgx(ipolar,iazimuth), rgy(ipolar,iazimuth), rgz(ipolar,iazimuth)/)
! convert to Rosca-lambert projection
            call LambertgetInterpolation(dc, scl, int(ipar(3)), int(ipar(3)), nix, niy, nixp, niyp, dx, dy, dxm, dym)
        
            acc_sum = 0.25*(accum_e(isampletilt,nix,niy) * dxm * dym + &
                            accum_e(isampletilt,nixp,niy) * dx * dym + &
                            accum_e(isampletilt,nix,niyp) * dxm * dy + &
                            accum_e(isampletilt,nixp,niyp) * dx * dy)
         
            weightfact(isig) = weightfact(isig) + acc_sum

        end do
    end do
end do

weightfact(1:nsig) = weightfact(1:nsig)/weightfact(1)

!===================================================================
! ------ perform interpolation from square lambert map
!===================================================================
scl = float(ipar(5))
ECPattern = 0.0
dn = nint(float(ipar(6))*0.01)
cn = dn
quat = Quaternion_T( q = (/ 1.0, 0.0, 0.0, 0.0 /) )
quatloop: do ip=1,ipar(6)
  do ii = imin, imax
    do jj = jmin, jmax

        dc(1:3) = klist(1:3,ii,jj)
        dc = dc/sqrt(sum(dc*dc))
        
        dp = DOT_PRODUCT(dc(1:3),(/sin(fpar(2)*dtoR),0.D0,cos(fpar(2)*dtoR)/))      
        if (dp .gt. 1.D0) dp = 1.0
        MCangle = acos(dp)*Rtod
        isig = int(MCangle) + 1
        if (isig .gt. nsig) isig = nsig

        isigp = isig + 1
        if (isigp .gt. nsig) isigp = nsig
        dx = MCangle - int(MCangle)
        dxm =  1.0 - dx
        
        wf = weightfact(isig) * dxm + weightfact(isigp) * dx
        wf = 1.0
        call quat%set_quats( quats(1:4,ip) )
        dc = quat%quat_Lp(dc)
        dc = dc/sqrt(sum(dc*dc))

        call LambertgetInterpolation(dc, scl, int(ipar(5)), int(ipar(5)), nix, niy, nixp, niyp, dx, dy, dxm, dym)

        if (dc(3).ge.0.D0) then 
            ECpattern(ii,jj,ip) = wf * ( mLPNHsum(nix,niy) * dxm * dym + &
                         mLPNHsum(nixp,niy) * dx * dym + &
                         mLPNHsum(nix,niyp) * dxm * dy + &
                         mLPNHsum(nixp,niyp) * dx * dy )

        else
            ECpattern(ii,jj,ip) =  wf * ( mLPSHsum(nix,niy) * dxm * dym + &
                         mLPSHsum(nixp,niy) * dx * dym + &
                         mLPSHsum(nix,niyp) * dxm * dy + &
                         mLPSHsum(nixp,niyp) * dx * dy )
        end if

    end do
  end do

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

end subroutine EMsoftCgetECPatterns



! !--------------------------------------------------------------------------
! !
! ! SUBROUTINE:EMsoftCgetMCOpenCL
! !
! !> @author Marc De Graef, Carnegie Mellon University
! !
! !> @brief This subroutine can be called by a C/C++ program as a standalone routine to compute Monte Carlo data
! !
! !> @details This subroutine provides a method to compute a Monte Carlo data set, normally computed
! !> with the EMMCOpenCL.f90 program.  The routine can be called from an external C/C++ program; 
! !> the routine provides a callback mechanism to update the calling program about computational 
! !> progress, as well as a cancel option.
! !>
! !> The routine is intended to be called from a C/C++ program, e.g., DREAM.3D.  This routine is a 
! !> simplified version of the core of the EMMCOpenCL program. 
! !>
! !> Since the HDF5 library with fortran90 support can only be a static library on Mac OS X, we must
! !> have the calling program read the .xtal HDF5 file and pass the necessary information on to
! !> this routine.  This is a workaround until the HDF group fixes the static library issue; DREAM.3D
! !> requires a dynamical HDF5 library, so for DREAM.3D and EMsoft to properly work together, the 
! !> callable routines in this file may not depend on any HDF code at all, either directly or indirectly.
! !>
! !> @param ipar array with integer input parameters
! !> @param fpar array with float input parameters
! !> @param atdata atom coordinate array
! !> @param attypes atom type array
! !> @param latparm lattice parameter array
! !> @param accum_e output array with Monte Carlo energy histogram
! !> @param accum_z output array with Monte Carlo depth histogram
! !
! !> @date 03/08/16 MDG 1.0 original
! !> @date 03/19/16 MDG 1.1 corrections to a few variable types
! !> @date 04/13/16 MDG 1.2 correction to accum_z array size due to changes in calling DREAM.3D filter
! !> @date 04/18/16 MDG 1.3 increased number of entries in ipar, fpar for compatibility with EMsoftCgetEBSDmaster routine
! !> @date 04/28/16 MDG 1.4 corrected error in indexing of init_seeds array; caused DREAM.3D to crash randomly
! !> @date 11/09/17 MDG 2.0 added spar string array to pass EMsoft configuration strings into the routine
! !--------------------------------------------------------------------------
! recursive subroutine EMsoftCgetMCOpenCL(ipar, fpar, spar, atompos, atomtypes, latparm, accum_e, accum_z, cproc, &
! objAddress, cancel) bind(c, name='EMsoftCgetMCOpenCL')    ! this routine is callable from a C/C++ program
! !DEC$ ATTRIBUTES DLLEXPORT :: EMsoftCgetMCOpenCL

! ! ipar components
! ! ipar(1) : integer(kind=irg)       :: nx  = (numsx-1)/2
! ! ipar(2) : integer(kind=irg)       :: globalworkgrpsz
! ! ipar(3) : integer(kind=irg)       :: num_el
! ! ipar(4) : integer(kind=irg)       :: totnum_el
! ! ipar(5) : integer(kind=irg)       :: multiplier
! ! ipar(6) : integer(kind=irg)       :: devid
! ! ipar(7) : integer(kind=irg)       :: platid
! ! ipar(8) : integer(kind=irg)       :: CrystalSystem
! ! ipar(9) : integer(kind=irg)       :: Natomtypes
! ! ipar(10): integer(kind=irg)       :: SpaceGroupNumber
! ! ipar(11): integer(kind=irg)       :: SpaceGroupSetting
! ! ipar(12): integer(kind=irg)       :: numEbins
! ! ipar(13): integer(kind=irg)       :: numzbins
! ! ipar(14): integer(kind=irg)       :: mcmode  ( 1 = 'full', 2 = 'bse1' )
! ! ipar(15): integer(kind=irg)       :: numangle
! ! ipar(16): integer(kind=irg)       :: nxten = nx/10
! ! other entries are not used

! ! fpar components
! ! fpar(1) : real(kind=dbl)          :: sig
! ! fpar(2) : real(kind=dbl)          :: omega
! ! fpar(3) : real(kind=dbl)          :: EkeV
! ! fpar(4) : real(kind=dbl)          :: Ehistmin
! ! fpar(5) : real(kind=dbl)          :: Ebinsize
! ! fpar(6) : real(kind=dbl)          :: depthmax
! ! fpar(7) : real(kind=dbl)          :: depthstep
! ! fpar(8) : real(kind=dbl)          :: sigstart
! ! fpar(9) : real(kind=dbl)          :: sigend
! ! fpar(10): real(kind=dbl)          :: sigstep
! ! other entries are not used

! ! spar components
! ! this routine needs the following parameters to be set:
! ! spar(23): OpenCLpathname
! ! spar(26): Randomseedfilename
! ! 

! use local
! use error 
! use configmod
! use constants
! use crystal
! use constants
! use symmetry
! use io
! use typedefs
! use clfortran
! use CLsupport
! use timing
! use,INTRINSIC :: ISO_C_BINDING


! IMPLICIT NONE

! integer(c_int32_t),INTENT(IN)           :: ipar(wraparraysize)
! real(kind=sgl),INTENT(IN)               :: fpar(wraparraysize)
! character(kind=c_char, len=1), target, INTENT(IN) :: spar(wraparraysize*fnlen)
! real(kind=sgl),INTENT(IN)               :: atompos(ipar(9),5)
! integer(kind=irg),INTENT(IN)            :: atomtypes(ipar(9))
! real(kind=sgl),INTENT(IN)               :: latparm(6)
! integer(kind=irg),INTENT(OUT)           :: accum_e(ipar(12),-ipar(1):ipar(1),-ipar(1):ipar(1))
! integer(kind=irg),INTENT(OUT)           :: accum_z(ipar(12),ipar(13),-ipar(16):ipar(16),-ipar(16):ipar(16))
! TYPE(C_FUNPTR), INTENT(IN), VALUE       :: cproc
! integer(c_size_t),INTENT(IN), VALUE     :: objAddress
! character(fnlen)                        :: clPath=''
! character(len=1),INTENT(IN)             :: cancel

! ! local variables and parameters
! type(unitcell)                          :: cell
! character(4)                            :: mode
! integer(kind=ill)                       :: i=0, j=0, k=0, io_int(1)=0, num_max=0, totnum_el=0, ipg=0, isave=0, istat=0
! integer(kind=irg)                       :: nx=0, numEbins=0, numzbins=0, numangle=0, iang=0, cn=0, dn=0, totn=0
! integer(kind=irg),target                :: globalworkgrpsz=0, num_el=0, steps=0
! integer(kind=8),target                  :: globalsize(2)=0, localsize(2)=0
! integer(kind=8)                         :: size_in_bytes=0,size_in_bytes_seeds=0

! real(kind=sgl),target                   :: dens=0, avA=0, avZ=0, omega=0, EkeV=0, sig=0, bseyield=0, io_real(3)
! real(kind=4),target                     :: density=0, Ze=0, at_wt=0, delta=0
! real(kind=8),parameter                  :: dtoR = 0.01745329251D0  ! pi/180
! real(kind=4),allocatable, target        :: Lamresx(:), Lamresy(:), depthres(:), energyres(:)

! integer(kind=4),allocatable             :: rnseeds(:)
! integer(kind=4),allocatable,target      :: init_seeds(:)
! integer(kind=4)                         :: idxy(2), iE=0, px=0, py=0, iz=0, nseeds=0, hdferr=0, tstart=0, trimSpace=0 ! auxiliary variables
! real(kind=4)                            :: cxyz(3), edis=0, bse=0, xy(2), xs=0, ys=0, zs=0, sclf=0 ! auxiliary variables
! real(kind=8)                            :: rand=0
! logical                                 :: f_exists=.FALSE.


! ! OpenCL variables
! integer(c_intptr_t),allocatable, target :: platform(:)
! integer(c_intptr_t),allocatable, target :: device(:)
! integer(c_intptr_t),target              :: context=0
! integer(c_intptr_t),target              :: command_queue=0
! integer(c_intptr_t),target              :: prog=0
! integer(c_intptr_t),target              :: kernel=0
! integer(c_intptr_t),target              :: LamX=0, LamY=0, LamZ=0, depth=0, energy=0, seeds=0
! type(c_ptr)                             :: event
! integer(c_int32_t)                      :: ierr=0, pcnt=0
! integer(c_size_t),target                :: slength=0
! integer(c_intptr_t),target              :: ctx_props(3)
! character(2),target                     :: kernelname=''
! character(19),target                    :: progoptions=''
! character(fnlen),target                 :: info='' ! info about the GPU
! integer(c_int64_t)                      :: cmd_queue_props=0

! integer, parameter                      :: iunit = 10
! integer, parameter                      :: source_length = 50000
! character(len=source_length),target     :: source=''
! character(len=source_length, KIND=c_char),TARGET :: csource=''
! type(c_ptr), target                     :: psource
! integer(c_int)                          :: nump=0, numd=0, irec=0, val=0,val1=0 ! auxiliary variables
! integer(c_size_t)                       :: cnum=0, cnuminfo=0
! character(fnlen)                        :: instring='', dataname='', fname='', sourcefile=''
! PROCEDURE(ProgressCallBack2), POINTER   :: proc
! character(250),target                   :: currentDir=''
! character(fnlen)                        :: emmcPath='', outname=''
! character(fnlen)                        :: randomSeedPath=''
! integer(c_int32_t),target               :: filestat=0
! INTEGER(kind=irg)                       :: getcwd, status
! CHARACTER(LEN=30)                       :: Format=''

! ! parameters to deal with the input string array spar
! type(ConfigStructureType)               :: CS

! nullify(proc)

! ! link the proc procedure to the cproc argument
! CALL C_F_PROCPOINTER (cproc, proc)

! ! the calling program passes a c-string array spar that we need to convert to the 
! ! standard EMsoft config structure for use inside this routine
! call C2F_configuration_strings(C_LOC(spar), CS)

! ! the following is necessitated by the fact that none of this code may 
! ! depend on HDF5 routines, so we need to cut-and-paste from various 
! ! other library routines to set things up so that we can compute the 
! ! density, and the average atomic number and atomic mass...

! ! copy all the unit cell parameters into the proper fields and compute the 
! ! density parameters needed by the Monte Carlo routine; then discard the cell structure
! ! lattice parameters
! cell%a = dble(latparm(1))
! cell%b = dble(latparm(2))
! cell%c = dble(latparm(3))
! cell%alpha = dble(latparm(4))
! cell%beta = dble(latparm(5))
! cell%gamma = dble(latparm(6))
! ! symmetry parameters
! cell%xtal_system = ipar(8)
! cell%SYM_SGset = ipar(11)
! cell%SYM_SGnum = ipar(10)
! if ((cell%SYM_SGnum.ge.143).and.(cell%SYM_SGnum.le.167)) then
!   cell%SG%SYM_trigonal = .TRUE.
! else
!   cell%SG%SYM_trigonal = .FALSE.
! end if 
! ! atom type and coordinate parameters
! cell%ATOM_ntype = ipar(9)
! cell%ATOM_type(1:cell%ATOM_ntype) = atomtypes(1:cell%ATOM_ntype) 
! cell%ATOM_pos(1:cell%ATOM_ntype,1:5) = atompos(1:cell%ATOM_ntype,1:5) 
! ! generate the symmetry operations
! cell%hexset = .FALSE.
! if (cell%xtal_system.eq.4) cell%hexset = .TRUE.
! if ((cell%xtal_system.eq.5).AND.(cell%SYM_SGset.ne.2)) cell%hexset = .TRUE.
! ! compute the metric matrices
!  call CalcMatrices(cell)
! ! First generate the point symmetry matrices, then the actual space group.
! ! Get the symmorphic space group corresponding to the point group
! ! of the actual space group
!  ipg=0
!  do i=1,32
!   if (SGPG(i).le.cell%SYM_SGnum) ipg=i
!  end do
! ! if the actual group is also the symmorphic group, then both 
! ! steps can be done simultaneously, otherwise two calls to 
! ! GenerateSymmetry are needed.
!  if (SGPG(ipg).eq.cell%SYM_SGnum) then
!   call GenerateSymmetry(cell,.TRUE.)
!  else
!   isave = cell%SYM_SGnum
!   cell%SYM_SGnum = SGPG(ipg)
!   call GenerateSymmetry(cell,.TRUE.)
!   cell%SYM_SGnum = isave
!   call GenerateSymmetry(cell,.FALSE.)
!  end if
! ! next we get all the atom positions
! call CalcPositions(cell,'v')
! ! and now we have all we need to compute the density, average A and average Z
! call CalcDensity(cell, dens, avZ, avA)

! ! and copy these values into the desired variables
! density = dble(dens)
! Ze = dble(avZ)
! at_wt = dble(avA)

! ! define a number of parameters
! steps = 300
! mode = 'full'
! if (ipar(14).ne.1) mode = 'bse1'   
! EkeV = sngl(fpar(3))
! !sig = mcnl%sig*dtoR    ! this is defined later on and depends on the mode
! omega = sngl(fpar(2))*dtoR
! globalworkgrpsz = ipar(2)
! num_el = int(ipar(3))  ! no. of electron simulation by one work item
! num_max = globalworkgrpsz*globalworkgrpsz*num_el ! total simulation in one loop
! totnum_el = ipar(4) * ipar(5) ! total number of electrons to simulate
! globalsize = (/ globalworkgrpsz, globalworkgrpsz /)
! numEbins =  int(ipar(12))
! numzbins =  int(ipar(13))
! nx = int(ipar(1))
! delta = dble(nx)
! size_in_bytes = num_max*sizeof(EkeV)
! size_in_bytes_seeds = 4*globalworkgrpsz*globalworkgrpsz*sizeof(EkeV)
! numangle = int(ipar(15))

! ! next allocate and initialize a couple of arrays
! allocate(Lamresx(num_max), Lamresy(num_max), depthres(num_max), energyres(num_max), stat=istat)
! depthres = 0.0
! energyres = 0.0
! Lamresx = 0.0
! Lamresy = 0.0
! accum_e = 0
! accum_z = 0

! !======================
! ! OpenCL INITIALIZATION
! !======================
! call CLinit_PDCCQ(platform, nump, int(ipar(7)), device, numd, int(ipar(6)), info, context, command_queue)

! !=====================
! ! BUILD THE KERNEL
! !=====================
! ! read the source file
! sourcefile='/EMMC.cl'
! emmcPath=trim(CS%OpenCLpathname)//trim(sourcefile)
! emmcPath=EMsoft_toNativePath(emmcPath)

! ! sourcefile = 'EMMC.cl'
! call CLread_source_file_wrapper(emmcPath, csource, slength)
! ! io_int(1) = slength
! ! call WriteValue('Kernel source length (characters) : ',io_int,1)

! ! we disable all screen output; perhaps we can feed error messages back to the calling program...

! ! create the program
! pcnt = 1
! psource = C_LOC(csource)
! prog = clCreateProgramWithSource(context, pcnt, C_LOC(psource), C_LOC(slength), ierr)
!  ! if(ierr /= CL_SUCCESS) call FatalError("clCreateProgramWithSource: ",'Error: cannot create program from source.')

! ! build the program
! ierr = clBuildProgram(prog, numd, C_LOC(device), C_NULL_PTR, C_NULL_FUNPTR, C_NULL_PTR)
! if (ierr.le.0) then
!   ierr = clGetProgramBuildInfo(prog, device(ipar(6)), CL_PROGRAM_BUILD_LOG, sizeof(source), C_LOC(source), cnum)
!   ! if(len(trim(source)) > 0) call Message(trim(source(1:cnum)),frm='(A)')
! endif
!  ! if(ierr /= CL_SUCCESS) call FatalError("clBuildProgram: ",'Error: cannot build program.')

! ! get the compilation log
! ierr = clGetProgramBuildInfo(prog, device(ipar(6)), CL_PROGRAM_BUILD_LOG, sizeof(source), C_LOC(source), cnum)
!  ! if(len(trim(source)) > 0) call Message(trim(source(1:cnum)),frm='(A)')
!  ! if(ierr /= CL_SUCCESS) call FatalError("clGetProgramBuildInfo: ",'Error building program.')

! ! if we get here, then the program build was successful and we can proceed with the creation of the kernel
!  ! call Message('Program Build Successful... Creating kernel')

! ! finally get the kernel and release the program
! kernelname = 'MC'//CHAR(0)
! kernel = clCreateKernel(prog, C_LOC(kernelname), ierr)
!  ! if(ierr /= CL_SUCCESS) call FatalError("clCreateKernel: ",'Error creating kernel MC.')

! ierr = clReleaseProgram(prog)
!  ! if(ierr /= CL_SUCCESS) call FatalError("clReleaseProgram: ",'Error releasing program.')

! open(unit = iunit, file = trim(EMsoft_toNativePath(CS%Randomseedfilename)), form='unformatted', status='old')
! read(iunit) nseeds
! allocate(rnseeds(nseeds))
! read(iunit) rnseeds
! close(unit=iunit,status='keep')

! ! the next error needs to be checked in the calling program
!  if (globalworkgrpsz**2 .gt. nseeds) call FatalError('EMMCOpenCL:','insufficient prime numbers')

! allocate(init_seeds(4*globalworkgrpsz*globalworkgrpsz),stat=istat)
! init_seeds = 0
! do i = 1,globalworkgrpsz
!     do j = 1,globalworkgrpsz
!         do k = 1,4
!             init_seeds(4*((i-1)*globalworkgrpsz+(j-1))+k) = rnseeds(4*((i-1)*globalworkgrpsz+j)+k)
!         end do
!     end do
! end do

! ! create device memory buffers
! LamX = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, C_NULL_PTR, ierr)
!  ! if(ierr /= CL_SUCCESS) call FatalError('clCreateBuffer: ','cannot allocate device memory for LamX.')

! LamY = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, C_NULL_PTR, ierr)
!  ! if(ierr /= CL_SUCCESS) call FatalError('clCreateBuffer: ','cannot allocate device memory for LamY.')

! depth = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, C_NULL_PTR, ierr)
!    ! if(ierr /= CL_SUCCESS) call FatalError('clCreateBuffer: ','cannot allocate device memory for depth.')

! energy = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, C_NULL_PTR, ierr)
!    ! if(ierr /= CL_SUCCESS) call FatalError('clCreateBuffer: ','cannot allocate device memory for energy.')

! seeds = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, C_NULL_PTR, ierr)
!  ! if(ierr /= CL_SUCCESS) call FatalError('clCreateBuffer: ','cannot allocate device memory for seeds.')

! ierr = clEnqueueWriteBuffer(command_queue, seeds, CL_TRUE, 0_8, size_in_bytes_seeds, C_LOC(init_seeds(1)), &
!                             0, C_NULL_PTR, C_NULL_PTR)
!  ! if(ierr /= CL_SUCCESS) call FatalError('clEnqueueWriteBuffer: ','cannot Enqueue write buffer.')

! ! set the callback parameters
! dn = 1
! cn = dn
! totn = numangle * (totnum_el/num_max+1)

! call Time_tick(tstart)

! ! loop over angles (used for BSE1, single run for full)
! angleloop: do iang = 1,numangle

!   if (mode .eq. 'bse1') then
!     sig = (fpar(8) + (iang-1)*fpar(10))*dtoR
!   else 
!     sig = fpar(1)*dtoR
!   end if

!   mainloop: do i = 1,(totnum_el/num_max+1)

! ! set the kernel arguments
!     ierr = clSetKernelArg(kernel, 0, sizeof(LamX), C_LOC(LamX))
!        ! if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument.'

!     ierr = clSetKernelArg(kernel, 1, sizeof(LamY), C_LOC(LamY))
!        ! if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument.'

!     ierr = clSetKernelArg(kernel, 2, sizeof(EkeV), C_LOC(EkeV))
!        ! if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument.'

!     ierr = clSetKernelArg(kernel, 3, sizeof(globalworkgrpsz), C_LOC(globalworkgrpsz))
!        ! if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument.'

!     ierr = clSetKernelArg(kernel, 4, sizeof(Ze), C_LOC(Ze))
!        ! if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument.'

!     ierr = clSetKernelArg(kernel, 5, sizeof(density), C_LOC(density))
!        ! if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument.'

!     ierr = clSetKernelArg(kernel, 6, sizeof(at_wt), C_LOC(at_wt))
!        ! if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument.'

!     ierr = clSetKernelArg(kernel, 7, sizeof(num_el), C_LOC(num_el))
!        ! if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument.'

!     ierr = clSetKernelArg(kernel, 8, sizeof(seeds), C_LOC(seeds))
!        ! if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument.'

!     ierr = clSetKernelArg(kernel, 9, sizeof(sig), C_LOC(sig))
!        ! if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument.'

!     ierr = clSetKernelArg(kernel, 10, sizeof(omega), C_LOC(omega))
!        ! if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument.'

!     ierr = clSetKernelArg(kernel, 11, sizeof(depth), C_LOC(depth))
!        ! if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument.'

!     ierr = clSetKernelArg(kernel, 12, sizeof(energy), C_LOC(energy))
!        ! if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument.'

!     ierr = clSetKernelArg(kernel, 13, sizeof(steps), C_LOC(steps))
!        ! if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument.'

! ! execute the kernel
! !   ierr = clEnqueueNDRangeKernel(command_queue, kernel, 2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), &
! !                                 0, C_NULL_PTR, C_NULL_PTR)
!     ierr = clEnqueueNDRangeKernel(command_queue, kernel, 2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
!                                   0, C_NULL_PTR, C_NULL_PTR)
!     ! if(ierr /= CL_SUCCESS) stop 'Error: clEnqueueNDRangeKernel'
! ! wait for the commands to finish
!     ierr = clFinish(command_queue)
!     ! if(ierr /= CL_SUCCESS) stop 'Error: clFinish'

! ! read the resulting vector from device memory
!     ierr = clEnqueueReadBuffer(command_queue,LamX,CL_TRUE,0_8,size_in_bytes,C_LOC(Lamresx(1)),0,C_NULL_PTR,C_NULL_PTR)
!     ! if(ierr /= CL_SUCCESS) stop 'Error: clEnqueueReadBuffer LamX '
!     ierr = clEnqueueReadBuffer(command_queue,LamY,CL_TRUE,0_8,size_in_bytes,C_LOC(Lamresy(1)),0,C_NULL_PTR,C_NULL_PTR)
!     ! if(ierr /= CL_SUCCESS) stop 'Error: clEnqueueReadBuffer LamY '
!     ierr = clEnqueueReadBuffer(command_queue,depth,CL_TRUE,0_8,size_in_bytes,C_LOC(depthres(1)),0,C_NULL_PTR,C_NULL_PTR)
!     ! if(ierr /= CL_SUCCESS) stop 'Error: clEnqueueReadBuffer depth '
!     ierr = clEnqueueReadBuffer(command_queue,energy,CL_TRUE,0_8,size_in_bytes,C_LOC(energyres(1)),0,C_NULL_PTR,C_NULL_PTR)
!     ! if(ierr /= CL_SUCCESS) stop 'Error: clEnqueueReadBuffer energy'

!     if (mode .eq. 'full') then
!       val = 0
!       subloopfull: do j = 1, num_max
!         if ((Lamresx(j) .ne. -10.0) .and. (Lamresy(j) .ne. -10.0) &
!           .and. (depthres(j) .ne. 10.0) .and. (energyres(j) .ne. 0.0) &
!           .and. .not.isnan(Lamresx(j)) .and. .not.isnan(Lamresy(j))) then
! ! and get the nearest pixel [ take into account reversal of coordinate frame (x,y) -> (y,-x) ]
!              ! if ((nint(delta*Lamresy(j)) .eq. 0.0) .and. (nint(-delta*Lamresx(j)) .eq. 0.0)) then
!              !   val1 = val1 + 1
!              ! end if

!              idxy = (/ nint(delta*Lamresy(j)), nint(-delta*Lamresx(j)) /)

!              if (maxval(abs(idxy)).le.nx) then
! ! If Ec larger than Emin, then we should count this electron
!                if (energyres(j).gt.fpar(4)) then

!                  val = val + 1
!                  iE = nint((energyres(j)-fpar(4))/fpar(5))+1
! ! first add this electron to the correct exit distance vs. energy bin (coarser than the angular plot)
!                  edis = abs(depthres(j))  ! distance from last scattering point to surface along trajectory
!                  iz = nint(edis/fpar(7)) +1
!                  if ( (iz.gt.0).and.(iz.le.ipar(13)) ) then

!                    px = nint(idxy(1)/10.0)
!                    py = nint(idxy(2)/10.0)
!                    accum_z(iE,iz,px,py) = accum_z(iE,iz,px,py) + 1

!                  end if
! ! then add it to the modified Lambert accumulator array.
!                  accum_e(iE,idxy(1),idxy(2)) = accum_e(iE,idxy(1),idxy(2)) + 1
!                end if
!              end if
!         end if
!       end do subloopfull
!     end if

!     if (mode .eq. 'bse1') then
!       subloopbse1: do j = 1, num_max

!         if ((Lamresx(j) .ne. -10.0) .and. (Lamresy(j) .ne. -10.0) &
!           .and. (depthres(j) .ne. 10.0) .and. (energyres(j) .ne. 0.0) &
!           .and. .not.isnan(Lamresx(j)) .and. .not.isnan(Lamresy(j))) then
! ! and get the nearest pixel [ take into account reversal of coordinate frame (x,y) -> (y,-x) ]
!           if ((nint(delta*Lamresy(j)) .eq. 0.0) .and. (nint(-delta*Lamresx(j)) .eq. 0.0)) then
!             val1 = val1 + 1
!           end if

!           val = val + 1
!           idxy = (/ nint(delta*Lamresy(j)), nint(-delta*Lamresx(j)) /)

!           if (maxval(abs(idxy)).le.nx) then
! ! first add this electron to the correct exit distance vs. sigma (coarser than the angular plot)
!             edis = abs(depthres(j))  ! distance from last scattering point to surface along trajectory
!             iz = nint(edis/fpar(7)) +1
!             if ( (iz.gt.0).and.(iz.le.ipar(13)) ) then
!               px = nint(idxy(1)/10.0)
!               py = nint(idxy(2)/10.0)
!               accum_z(iang,iz,px,py) = accum_z(iang,iz,px,py) + 1

!             end if
! ! then add it to the modified Lambert accumulator array.
!             accum_e(iang,idxy(1),idxy(2)) = accum_e(iang,idxy(1),idxy(2)) + 1
!           end if
!         end if
!       end do subloopbse1
!     end if

! ! has the cancel flag been set by the calling program ?
!     if(cancel.ne.char(0)) then 
!       EXIT angleloop
!     end if

! ! update the progress counter and report it to the calling program via the proc callback routine
!     if(objAddress.ne.0) then
!       bseyield = 100.0*float(sum(accum_e))/float(i*num_max)
!       call proc(objAddress, cn, totn, bseyield)
!       cn = cn+dn
!     end if

!   end do mainloop
! end do angleloop 

! !=====================
! ! RELEASE EVERYTHING
! !=====================

! ierr = clReleaseKernel(kernel)
! ierr = clReleaseCommandQueue(command_queue)
! ierr = clReleaseContext(context)
! ierr = clReleaseMemObject(LamX)
! ierr = clReleaseMemObject(LamY)
! ierr = clReleaseMemObject(depth)
! ierr = clReleaseMemObject(energy)
! ierr = clReleaseMemObject(seeds)


! end subroutine EMsoftCgetMCOpenCL



    
end module mod_SEMwrappers
    