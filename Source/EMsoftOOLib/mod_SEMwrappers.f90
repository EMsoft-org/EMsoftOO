! ###################################################################
! Copyright (c) 2013-2024, Marc De Graef Research Group/Carnegie Mellon University
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
! The following is the mapping:  [NOTE THAT THESE NEED TO BE REVISED TO THE INDICES AT THE TOP OF THIS FILE !!!]
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


!--------------------------------------------------------------------------
!
! SUBROUTINE:EMsoftCgetEBSDmaster
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief This subroutine can be called by a C/C++ program as a standalone routine to compute EBSD master patterns
!
!> @details This subroutine provides a method to compute an EBSD master pattern for the northern and southern
!> hemispheres, i.e., it implements the EMEBSDmaster.f90 program.  The routine can be called from an external C/C++ program; 
!> the routine provides a callback mechanism to update the calling program about computational 
!> progress, as well as a cancel option.
!>
!> The routine is intended to be called from a C/C++ program, e.g., DREAM.3D.  This routine is a 
!> simplified version of the core of the EMEBSDmaster program. 
!>
!> Since the HDF5 library with fortran90 support can only be a static library on Mac OS X, we must
!> have the calling program read the .xtal HDF5 file and pass the necessary information on to
!> this routine.  This is a workaround until the HDF group fixes the static library issue; DREAM.3D
!> requires a dynamical HDF5 library, so for DREAM.3D and EMsoft to properly work together, the 
!> callable routines in this file may not depend on any HDF code at all, either directly or indirectly.
!>
!> @param ipar array with integer input parameters
!> @param fpar array with float input parameters
!> @param atdata atom coordinate array
!> @param attypes atom type array
!> @param latparm lattice parameter array
!> @param accum_z output array with Monte Carlo depth histogram
!> @param mLPNH modified Lambert projection northern hemisphere (output)
!> @param mLPSH modified Lambert projection southern hemisphere (output)
!
!> @date 04/17/16 MDG 1.0 original
!> @date 08/16/18 MDG 1.1 added 'uniform' parameter [ipar(36)] to generate only the background
!> @date 04/15/21 MDG 2.0 conversion to Object Oriented version
!--------------------------------------------------------------------------
recursive subroutine EMsoftCgetEBSDmaster(ipar,fpar,atompos,atomtypes,latparm,accum_z,mLPNH,mLPSH,cproc,objAddress,cancel) &
           bind(c, name='EMsoftCgetEBSDmaster')    ! this routine is callable from a C/C++ program
!DEC$ ATTRIBUTES DLLEXPORT :: EMsoftCgetEBSDmaster

! these are the same as in the EMsoftCgetMCOpenCL routine, with a few extras at the end.
! ipar components
! ipar(1) : integer(kind=irg)       :: nx  = (numsx-1)/2
! ipar(2) : integer(kind=irg)       :: globalworkgrpsz
! ipar(3) : integer(kind=irg)       :: num_el
! ipar(4) : integer(kind=irg)       :: totnum_el
! ipar(5) : integer(kind=irg)       :: multiplier
! ipar(6) : integer(kind=irg)       :: devid
! ipar(7) : integer(kind=irg)       :: platid
! ipar(8) : integer(kind=irg)       :: CrystalSystem
! ipar(9) : integer(kind=irg)       :: Natomtypes
! ipar(10): integer(kind=irg)       :: SpaceGroupNumber
! ipar(11): integer(kind=irg)       :: SpaceGroupSetting
! ipar(12): integer(kind=irg)       :: numEbins
! ipar(13): integer(kind=irg)       :: numzbins
! ipar(14): integer(kind=irg)       :: mcmode  ( 1 = 'full', 2 = 'bse1' )
! ipar(15): integer(kind=irg)       :: numangle
! ipar(16): integer(kind=irg)       :: nxten = nx/10
! the following are only used in this routine, not in the Monte Carlo routine
! ipar(17): integer(kind=irg)       :: npx
! ipar(18): integer(kind=irg)       :: nthreads
! ipar(36): integer(kind=irg)       :: uniform  ['1' = yes (background only), '0' = no ]

! fpar components
! fpar(1) : real(kind=dbl)          :: sig
! fpar(2) : real(kind=dbl)          :: omega
! fpar(3) : real(kind=dbl)          :: EkeV
! fpar(4) : real(kind=dbl)          :: Ehistmin
! fpar(5) : real(kind=dbl)          :: Ebinsize
! fpar(6) : real(kind=dbl)          :: depthmax
! fpar(7) : real(kind=dbl)          :: depthstep
! fpar(8) : real(kind=dbl)          :: sigstart
! fpar(9) : real(kind=dbl)          :: sigend
! fpar(10): real(kind=dbl)          :: sigstep
! parameters only used in this routine, this includes the Bethe Parameters !!!!
! fpar(11) : real(kind=dbl)         :: dmin
! fpar(12) : real(kind=dbl)         :: Bethe  c1
! fpar(13) : real(kind=dbl)         :: Bethe  c2
! fpar(14) : real(kind=dbl)         :: Bethe  c3

use mod_EMsoft
use mod_initializers
use mod_symmetry
use mod_crystallography
use mod_gvectors
use mod_kvectors
use mod_io
use mod_math
use mod_diffraction
use mod_timing
use mod_Lambert
use ISO_C_BINDING
use omp_lib
use mod_OMPsupport
use mod_memory
use mod_notifications
use stringconstants
use mod_MCfiles

IMPLICIT NONE

integer(c_int32_t),INTENT(IN)           :: ipar(wraparraysize)
real(kind=sgl),INTENT(IN)               :: fpar(wraparraysize)
real(kind=sgl),INTENT(IN)               :: atompos(ipar(9),5)
integer(kind=irg),INTENT(IN)            :: atomtypes(ipar(9))
real(kind=sgl),INTENT(IN)               :: latparm(6)
integer(kind=irg),INTENT(IN)            :: accum_z(ipar(12),ipar(13),-ipar(16):ipar(16),-ipar(16):ipar(16))
real(kind=sgl),INTENT(INOUT)            :: mLPNH(-ipar(17):ipar(17),-ipar(17):ipar(17),1:ipar(12),1:ipar(9))
real(kind=sgl),INTENT(INOUT)            :: mLPSH(-ipar(17):ipar(17),-ipar(17):ipar(17),1:ipar(12),1:ipar(9))
TYPE(C_FUNPTR), INTENT(IN), VALUE       :: cproc
integer(c_size_t),INTENT(IN), VALUE     :: objAddress
character(len=1),INTENT(IN)             :: cancel

type(Cell_T)            :: cell
type(DynType)           :: Dyn
type(Timing_T)          :: timer
type(IO_T)              :: Message
type(Lambert_T)         :: L
type(SpaceGroup_T)      :: SG
type(Diffraction_T)     :: Diff
type(kvectors_T)        :: kvec
type(gvectors_T)        :: reflist
type(memory_T)          :: mem, memth

real(kind=dbl)          :: ctmp(192,3), arg

integer(kind=irg)       :: isym,i,j,ik,npy,ipx,ipy,ipz,debug,iE,izz, izzmax, iequiv(3,48), nequiv, num_el, MCnthreads, & ! counters
                           numk, numsites, & 
                           ir,nat(maxpasym),kk(3), skip, ijmax, one, NUMTHREADS, TID, SamplingType, cancelerr, &
                           numset,n,ix,iy,iz, nns, nnw, nref, Estart, ipg, isave, npx, nthreads,  &
                           istat,gzero,ic,ip,ikk, totstrong, totweak, jh, ierr, nix, niy, nixp, niyp, nxten     ! counters
real(kind=dbl)          :: tpi,Znsq, kkl, DBWF, kin, delta, h, lambda, omtl, srt, dc(3), xy(2), edge, scl, tmp, &
                           dx, dxm, dy, dym, dmin !
real(kind=sgl)          :: io_real(5), selE, kn, FN(3), kkk(3), tstart, tstop, bp(4), nabsl
real(kind=sgl),allocatable      :: EkeVs(:), svals(:), auxNH(:,:,:,:), auxSH(:,:,:,:)  ! results
complex(kind=dbl)               :: czero
complex(kind=dbl),allocatable   :: Lgh(:,:), Sgh(:,:,:)
logical                 :: usehex, switchmirror, verbose
real(kind=dbl)          :: apositions(maxpasym,5)
! Monte Carlo derived quantities
integer(kind=irg)       :: numEbins, numzbins, nsx, nsy, hdferr, nlines, lastEnergy, cn, dn, totn, cn2, totn2 ! variables used in MC energy file
real(kind=dbl)          :: EkeV, Ehistmin, Ebinsize, depthmax, depthstep, etotal ! enery variables from MC program
integer(kind=irg),allocatable :: thick(:), acc_z(:,:,:,:)
real(kind=sgl),allocatable :: lambdaE(:,:)
logical                 :: f_exists, readonly, overwrite=.TRUE., insert=.TRUE., stereog, uniform
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)
character(fnlen,kind=c_char)                     :: line2(1)
character(fnlen)                                 :: pname, pdesc
integer(kind=irg)       :: imh, imk, iml, gg(3)
real(kind=sgl)          :: dhkl, ddt

type(gnode),save                :: rlp
type(kvectorlist), pointer      :: ktmp
type(reflisttype), pointer      :: firstw
real(kind=sgl),allocatable      :: karray(:,:)
integer(kind=irg),allocatable   :: kij(:,:)
complex(kind=dbl),allocatable   :: DynMat(:,:)
character(fnlen)                :: dataset, instring
PROCEDURE(ProgressCallBack3), POINTER   :: proc

!$OMP THREADPRIVATE(rlp) 

nullify(proc)

! link the proc procedure to the cproc argument
CALL C_F_PROCPOINTER (cproc, proc)

! initalize a few variables
tpi = 2.D0*cPi
czero = dcmplx(0.D0,0.D0)

! parameters that would normally be read from the MC HDF5 file
npx = ipar(17)
nxten = ipar(16)
EkeV = fpar(3)
Ehistmin = fpar(4)
Ebinsize = fpar(5)
depthmax = fpar(6)
depthstep = fpar(7)
numEbins = ipar(12)
Estart = numEbins
numzbins = ipar(13)
num_el = ipar(3)
dmin = fpar(11)
nthreads = ipar(18)
etotal = dble(ipar(4))*dble(ipar(5))
numsites = ipar(9)
numset = numsites
uniform = .FALSE.
if (ipar(36).eq.1) uniform = .TRUE.

nsx = ipar(1) ! (mcnl%numsx - 1)/2
nsy = nsx

!=============================================
!=============================================
! crystallography section

! the following is necessitated by the fact that none of this code may 
! depend on HDF5 routines, so we need to cut-and-paste from various 
! other library routines to set things up so that we can compute the 
! density, and the average atomic number and atomic mass...

! copy all the unit cell parameters into the proper fields and compute the 
! density parameters needed by the Monte Carlo routine; 
! initialize cell and SG classes 
cell = Cell_T( dble(latparm) )
SG = SpaceGroup_T( SGnumber = ipar(10), xtalSystem = ipar(8), setting = ipar(11), &
                   dmt = cell%getdmt(), rmt = cell%getrmt() )
! fill in additional symmetry and cell parameters
if ((ipar(10).ge.143).and.(ipar(10).le.167)) then
  call SG%setSpaceGrouptrigonal(.TRUE.)
else
  call SG%setSpaceGrouptrigonal(.FALSE.)
end if 
! atom type and coordinate parameters
call cell%setNatomtype( ipar(9) )
call cell%setAtomtype( atomtypes(1:ipar(9)) )
apositions(1:ipar(9), 1:5) = atompos
call cell%setAtomPos( apositions )
! generate the symmetry operations
call SG%setSpaceGrouphexset( .FALSE. )
if (ipar(8).eq.4) call SG%setSpaceGrouphexset( .TRUE. )
if ((ipar(8).eq.5).AND.(ipar(11).ne.2)) call SG%setSpaceGrouphexset( .TRUE. )
! Get the symmorphic space group corresponding to the point group
! of the actual space group
ipg = SG%getPGnumber()
! if the actual group is also the symmorphic group, then both 
! steps can be done simultaneously, otherwise two calls to 
! GenerateSymmetry are needed.
if (SGPG(ipg).eq.ipar(10)) then
  call SG%GenerateSymmetry( .TRUE., cell%getdmt(), cell%getrmt() )
else
  isave = SG%getSpaceGroupNumber()
  call SG%setSpaceGroupNumber( SGPG(ipg) )
  call SG%GenerateSymmetry( .TRUE., cell%getdmt(), cell%getrmt() )
  call SG%setSpaceGroupNumber( int(isave) )
  call SG%GenerateSymmetry( .FALSE., cell%getdmt(), cell%getrmt() )
end if
! next we get all the atom positions
call cell%CalcPositions( SG, 'v' )

! initialize the memory class 
mem = memory_T()

! set the method for Fourier coefficient computation to Weickenmeier-Kohl
call Diff%setrlpmethod( 'WK' )
call Diff%setV( EkeV )
! next we need to initialize the cell but without reading any HDF files ... 
call Initialize_Cell_NoHDF(cell, Diff, SG, Dyn, fpar(11))

! set the BetheParameters ... 
call Diff%SetBetheParameters( usevalues = (/ sngl(fpar(12)), sngl(fpar(13)), sngl(fpar(14)), 1.0_sgl /) )

! allocate and compute the Sgh loop-up table
numset = cell%getNatomtype()
nat = 0
call Diff%Initialize_SghLUT(cell, SG, fpar(11), numset, nat, verbose)

! determine the point group number
j=0
do i=1,32
if (SGPG(i).le.SG%getSpaceGroupNumber()) j=i
end do
isym = j

! here is new code dealing with all the special cases (quite a few more compared to the 
! Laue group case)...  isym is the point group number. Once the symmetry case has been
! fully determined (taking into account things like 31m and 3m1 an such), then the only places
! that symmetry is handled are the modified Calckvectors routine, and the filling of the modified
! Lambert projections after the dynamical simulation step.  We are also changing the name of the 
! sr array (or srhex) to mLPNH and mLPSH (modified Lambert Projection Northern/Southern Hemisphere).

! Here, we encode isym into a new number that describes the sampling scheme; the new schemes are 
! described in detail in the EBSD manual pdf file.
SamplingType = PGSamplingType(isym)

! next, intercept the special cases (hexagonal vs. rhombohedral cases that require special treatment)
if ((SamplingType.eq.-1).or.(isym.eq.14).or.(isym.eq.26)) then
  SamplingType = SG%getHexvsRho(isym)
end if

! if the point group is trigonal or hexagonal, we need to switch usehex to .TRUE. so that
! the program will use the hexagonal sampling method
usehex = .FALSE.
if ((SG%getSpaceGroupXtalSystem().eq.4).or.(SG%getSpaceGroupXtalSystem().eq.5)) usehex = .TRUE.

! ---------- end of symmetry and crystallography section
!=============================================
!=============================================

!=============================================
!=============================================
! this is where we determine the value for the thickness integration limit for the CalcLgh3 routine...
mem = memory_T()
call mem%alloc(EkeVs, (/ numEbins /), 'EkeVs')
call mem%alloc(thick, (/ numEbins /), 'thick')
do i=1,numEbins
  EkeVs(i) = Ehistmin + float(i-1)*Ebinsize
end do

! then, for each energy determine the 95% histogram thickness
izzmax = 0
do iE = 1,numEbins
 do ix=-nsx/10,nsx/10
  do iy=-nsy/10,nsy/10
   istat = sum(accum_z(iE,:,ix,iy))
   izz = 1
   do while (sum(accum_z(iE,1:izz,ix,iy)).lt.(0.99*istat))
    izz = izz+1
   end do
   if (izz.gt.izzmax) izzmax = izz
  end do
 end do
 thick(iE) = dble(izzmax) * fpar(7)
end do

izz = nint(maxval(thick)/fpar(7))
call mem%alloc(lambdaE, (/ numEbins, izz/), 'lambdaE')
do iE=1,numEbins
 call Diff%setV(dble(Ekevs(iE)))
 call Diff%CalcWaveLength( cell )
 call Diff%CalcUcg(cell,(/0,0,0/))
 rlp = Diff%getrlp()
 nabsl = rlp%xgp
 do iz=1,izz
  lambdaE(iE,iz) = float(sum(accum_z(iE,iz,-nsx/10:nsx/10,-nsy/10:nsy/10)))/etotal
  lambdaE(iE,iz) = lambdaE(iE,iz) * exp(2.0*sngl(cPi)*(iz-1)*fpar(7)/nabsl)
 end do
end do

! ---------- end of 'read Monte Carlo output file and extract necessary parameters' section
!=============================================
!=============================================

!=============================================
!=============================================
! ---------- a couple of initializations
   numset = cell%getNatomtype()
   npy = npx
   gzero = 1  ! index of incident beam
   debug = 0  ! no longer used
! ----------
!=============================================
!=============================================

!=============================================
!=============================================

! set various arrays to zero or 1, depending on uniform parameter
if (uniform.eqv..TRUE.) then
   mLPNH = 1.0
   mLPSH = 1.0
else
   mLPNH = 0.0
   mLPSH = 0.0
end if

! set the callback parameters
dn = 1
cn = dn
cn2 = 0
totn2 = Estart

!=============================================
!=============================================
! ---------- from here on, we need to repeat the entire computation for each energy value, assuming that uniform = .FALSE.
if (uniform.eqv..FALSE.) then
  reflist = gvectors_T()

! instantiate the memory class for the OpenMP section
  memth = memory_T( nt = ipar(18), silent=.TRUE. )

  cancelerr = 0
  energyloop: do iE=Estart,1,-1
  cn2 = cn2+dn 
! if this is not the first run, then set the accelerating voltage and recompute the Fourier coefficients
  if (iE.ne.Estart) then 
   call Diff%setV(dble(EkeVs(iE)))
   call Diff%setrlpmethod('WK')
   call Initialize_Cell_NoHDF(cell, Diff, SG, Dyn, fpar(11))
  end if 

!=============================================
! ---------- create the incident beam directions list
! determine all independent incident beam directions (use a linked list starting at khead)
! numk is the total number of k-vectors to be included in this computation;
! note that this needs to be redone for each energy, since the wave vector changes with energy
   kvec = kvectors_T()   ! initialize the wave vector list
   call kvec%set_kinp( (/ 0.D0, 0.D0, 1.D0 /) )
   call kvec%set_ktmax( 0.D0 )
   call kvec%set_SamplingType( SamplingType )
   call kvec%set_mapmode('RoscaLambert')
   if (usehex) then
     call kvec%Calckvectors(cell, SG, Diff, (/ 0.D0, 0.D0, 0.D0 /),npx,npy, ijmax,usehex)
   else
     call kvec%Calckvectors(cell, SG, Diff, (/ 0.D0, 0.D0, 0.D0 /),npx,npy, ijmax,usehex)
   end if
   numk = kvec%get_numk()
   totn = numk
   cn = dn

! convert part of the kvector linked list into arrays for OpenMP
   call mem%alloc(karray, (/ 4,numk /), 'karray')
   call mem%alloc(kij, (/ 3,numk /), 'kij')
! point to the first beam direction
   ktmp => kvec%get_ListHead()
! and loop through the list, keeping k, kn, and i,j
   karray(1:3,1) = sngl(ktmp%k(1:3))
   karray(4,1) = sngl(ktmp%kn)
   kij(1:3,1) = (/ ktmp%i, ktmp%j, ktmp%hs /)
   do ik=2,numk
     ktmp => ktmp%next
     karray(1:3,ik) = sngl(ktmp%k(1:3))
     karray(4,ik) = sngl(ktmp%kn)
     kij(1:3,ik) = (/ ktmp%i, ktmp%j, ktmp%hs /)
   end do
! and remove the linked list
   call kvec%Delete_kvectorlist()

   verbose = .FALSE.
   totstrong = 0
   totweak = 0

! ---------- end of "create the incident beam directions list"
!=============================================

! here's where we introduce the OpenMP calls, to spead up the overall calculations...

! set the number of OpenMP threads 
    call OMP_SET_NUM_THREADS( ipar(18) )

! use OpenMP to run on multiple cores ... 
!$OMP PARALLEL COPYIN(rlp) &
!$OMP& PRIVATE(DynMat,Sgh,Lgh,ik,FN,TID,kn,ipx,ipy,ix,iequiv,nequiv,reflist,firstw) &
!$OMP& PRIVATE(kkk,nns,nnw,nref,svals) SHARED(cancelerr)

    NUMTHREADS = OMP_GET_NUM_THREADS()
    TID = OMP_GET_THREAD_NUM()

  call memth%alloc(svals, (/ numset /), 'svals', TID=TID)

!$OMP DO SCHEDULE(DYNAMIC,100)    
! ---------- and here we start the beam direction loop
     beamloop:do ik = 1,numk

!=============================================
! ---------- create the master reflection list for this beam direction
! Then we must determine the masterlist of reflections (also a linked list);
! This list basically samples a large reciprocal space volume; it does not 
! distinguish between zero and higher order Laue zones, since that 
! distinction becomes meaningless when we consider the complete 
! reciprocal lattice.  
      reflist = gvectors_T()
      kkk = karray(1:3,ik)
      FN = kkk

      call reflist%Initialize_ReflectionList(cell, SG, Diff, FN, kkk, fpar(11))
      nref = reflist%get_nref()
! ---------- end of "create the master reflection list"
!=============================================

! determine strong and weak reflections
       nullify(firstw)
       nns = 0
       nnw = 0
       call reflist%Apply_BethePotentials(Diff, firstw, nns, nnw)

! generate the dynamical matrix
       call memth%alloc(DynMat, (/nns,nns/), 'DynMat', TID=TID)
       call reflist%GetDynMat(cell, Diff, firstw, DynMat, nns, nnw)
       totstrong = totstrong + nns
       totweak = totweak + nnw

! then we need to initialize the Sgh and Lgh arrays
       call memth%alloc(Sgh, (/ nns,nns,numset /), 'Sgh', TID=TID)
       call memth%alloc(Lgh, (/ nns,nns /), 'Lgh', TID=TID)
       Sgh = czero
       Lgh = czero
       call reflist%getSghfromLUT(Diff,nns,numset,Sgh)

! solve the dynamical eigenvalue equation for this beam direction  
       kn = karray(4,ik)
       call reflist%CalcLgh(DynMat,Lgh,dble(thick(iE)),dble(kn),nns,gzero,depthstep,lambdaE(iE,1:izzmax),izzmax)

! sum over the element-wise (Hadamard) product of the Lgh and Sgh arrays 
       svals = 0.0
       do ix=1,numset
         svals(ix) = real(sum(Lgh(1:nns,1:nns)*Sgh(1:nns,1:nns,ix)))
       end do
       svals = svals/float(sum(nat(1:numset)))

! and store the resulting svals values, applying point group symmetry where needed.
       ipx = kij(1,ik)
       ipy = kij(2,ik)
       ipz = kij(3,ik)
!
     if (usehex) then
       call L%Apply3DPGSymmetry(cell,SG,ipx,ipy,ipz,npx,iequiv,nequiv,usehex)
     else
       if ((SG%getSpaceGroupNumber().ge.195).and.(SG%getSpaceGroupNumber().le.230)) then
         call L%Apply3DPGSymmetry(cell,SG,ipx,ipy,ipz,npx,iequiv,nequiv,cubictype=SamplingType)
       else
         call L%Apply3DPGSymmetry(cell,SG,ipx,ipy,ipz,npx,iequiv,nequiv)
       end if
     end if
!$OMP CRITICAL
       do ix=1,nequiv
         if (iequiv(3,ix).eq.-1) mLPSH(iequiv(1,ix),iequiv(2,ix),iE,1:numset) = svals(1:numset)
         if (iequiv(3,ix).eq.1) mLPNH(iequiv(1,ix),iequiv(2,ix),iE,1:numset) = svals(1:numset)
       end do
!$OMP END CRITICAL
  
     call reflist%Delete_gvectorlist()
     call memth%dealloc(Sgh, 'Sgh', TID=TID)
     call memth%dealloc(Lgh, 'Sgh', TID=TID)
     call memth%dealloc(DynMat, 'DynMat', TID=TID)

! has the cancel flag been set by the calling program ?
!!!!$OMP CANCELLATION POINT
     if(cancel.ne.char(0)) then
!$OMP ATOMIC WRITE
         cancelerr = 1
!$OMP CANCEL DO
      end if 

! update the progress counter and report it to the calling program via the proc callback routine
!$OMP CRITICAL
     if(objAddress.ne.0) then
       cn = cn+dn
       if (mod(cn,1000).eq.0) then 
         call proc(objAddress, cn, totn, cn2, totn2)
       end if
     end if
!$OMP END CRITICAL

! for testing purposes ... the wrapper should not produce any command line output
    ! if (mod(ik,5000).eq.0) write (*,*) ik,' completed out of ',numk 

    end do beamloop

    call memth%dealloc(svals, 'svals', TID=TID)

! end of OpenMP portion
!$OMP END PARALLEL
  
! was the Cancel button pressed in the calling program?
    if(cancelerr.ne.0) EXIT energyloop

    call mem%dealloc(karray, 'karray')
    call mem%dealloc(kij, 'kij')

   if (usehex) then
! and finally, we convert the hexagonally sampled array to a square Lambert projection which will be used 
! for all EBSD pattern interpolations;  we need to do this for both the Northern and Southern hemispheres

! we begin by allocating auxiliary arrays to hold copies of the hexagonal data; the original arrays will
! then be overwritten with the newly interpolated data.
    call mem%alloc(auxNH, (/npx,npy,1,numsites/), 'auxNH', startdims=(/-npx,-npy,1,1/))
    call mem%alloc(auxSH, (/npx,npy,1,numsites/), 'auxSH', startdims=(/-npx,-npy,1,1/))
    auxNH = mLPNH
    auxSH = mLPSH

    edge = 1.D0 / dble(npx)
    scl = float(npx)
    do i=-npx,npx
      do j=-npy,npy
! determine the spherical direction for this point
        L = Lambert_T( xyd = (/ dble(i), dble(j) /) * edge )
        ierr = L%LambertSquareToSphere(dc)
! convert direction cosines to hexagonal Lambert projections
        L = Lambert_T( xyzd = dc )
        ierr = L%LambertSphereToHex(xy)
        xy = xy * scl
! interpolate intensity from the neighboring points
        if (ierr.eq.0) then
          nix = floor(xy(1))
          niy = floor(xy(2))
          nixp = nix+1
          niyp = niy+1
          if (nixp.gt.npx) nixp = nix
          if (niyp.gt.npx) niyp = niy
          dx = xy(1) - nix
          dy = xy(2) - niy
          dxm = 1.D0 - dx
          dym = 1.D0 - dy
          mLPNH(i,j,1,1:numsites) = auxNH(nix,niy,1,1:numsites)*dxm*dym + auxNH(nixp,niy,1,1:numsites)*dx*dym + &
                               auxNH(nix,niyp,1,1:numsites)*dxm*dy + auxNH(nixp,niyp,1,1:numsites)*dx*dy
          mLPSH(i,j,1,1:numsites) = auxSH(nix,niy,1,1:numsites)*dxm*dym + auxSH(nixp,niy,1,1:numsites)*dx*dym + &
                               auxSH(nix,niyp,1,1:numsites)*dxm*dy + auxSH(nixp,niyp,1,1:numsites)*dx*dy
        end if
      end do
    end do
    call mem%dealloc(auxNH, 'auxNH')
    call mem%dealloc(auxSH, 'auxSH')
   end if

! make sure that the outer pixel rim of the mLPSH patterns is identical to
! that of the mLPNH array.
   mLPSH(-npx,-npx:npx,iE,1:numset) = mLPNH(-npx,-npx:npx,iE,1:numset)
   mLPSH( npx,-npx:npx,iE,1:numset) = mLPNH( npx,-npx:npx,iE,1:numset)
   mLPSH(-npx:npx,-npx,iE,1:numset) = mLPNH(-npx:npx,-npx,iE,1:numset)
   mLPSH(-npx:npx, npx,iE,1:numset) = mLPNH(-npx:npx, npx,iE,1:numset)

  end do energyloop

end if ! (uniform.eqv..FALSE.)

! that's the end of it...

end subroutine EMsoftCgetEBSDmaster

!--------------------------------------------------------------------------
recursive subroutine getKosselPatterns(ipar, fpar, Kosselpattern, quats, mLPNH, mLPSH)
!DEC$ ATTRIBUTES DLLEXPORT :: getKosselPatterns

! the input parameters are all part of a ipar and fpar input arrays instead of the usual namelist structures.
! The following is the mapping:
!
! ipar(1) = 1 if GetVectorsCone detector arrays need to be computed, 0 if not (arrays will have save status)
! ipar(2) = detnumsx
! ipar(3) = mpnpx
! ipar(4) = numquats
! ipar(5) = numdepths
! ipar(6) = depthsel

! fpar(1) = ecpnl%thetac

use mod_EMsoft
use mod_kinds
use mod_Lambert
use mod_quaternions
use,INTRINSIC :: ISO_C_BINDING

IMPLICIT NONE

integer(c_size_t),PARAMETER             :: nipar=6
integer(c_size_t),PARAMETER             :: nfpar=1
integer(c_size_t),PARAMETER             :: nq=4
integer(c_size_t),INTENT(IN)            :: ipar(nipar)
real(kind=sgl),INTENT(IN)               :: fpar(nfpar)
real(kind=sgl),INTENT(OUT)              :: Kosselpattern(ipar(2),ipar(2),ipar(4))
real(kind=sgl),INTENT(IN)               :: quats(nq,ipar(4))
real(kind=sgl),INTENT(IN)               :: mLPNH(-ipar(3):ipar(3), -ipar(3):ipar(3),ipar(5))
real(kind=sgl),INTENT(IN)               :: mLPSH(-ipar(3):ipar(3), -ipar(3):ipar(3),ipar(5))

type(Quaternion_T)                      :: quat

real(kind=sgl),allocatable,save         :: klist(:,:,:)

real(kind=dbl),parameter                :: Rtod = 57.2957795131D0
real(kind=dbl),parameter                :: dtoR = 0.01745329251D0

real(kind=sgl)                          :: kk(3), thetacr, ktmax, delta
integer(kind=irg)                       :: istat, imin, imax, jmin, jmax, ii ,jj, nsig, ip
integer(kind=irg)                       :: isig, nix, niy, nixp, niyp, isigp
real(kind=sgl)                          :: dc(3), scl, ixy(2), dx, dy, dxm, dym, dp


!==================================================================================
! ------ generate the detector klist array if needed 
!------- (calling program must decide this via ipar(1))
!==================================================================================

if (ipar(1).ge.1) then
    if (allocated(klist)) deallocate(klist)
    allocate(klist(1:3,-ipar(2):ipar(2),-ipar(2):ipar(2)), stat=istat)
    kk = (/0.0,0.0,1.0/)
    thetacr = DtoR*fpar(1)
    ktmax = tan(thetacr)
    delta = 2.0*ktmax/dble(ipar(2)-1)

    imin = 1
    imax = ipar(2)
    jmin = 1
    jmax = ipar(2)
     
    do ii = imin, imax
        do jj = jmin, jmax
            klist(1:3,ii,jj) = (/-ktmax+delta*(ii-1),-ktmax+delta*(jj-1),0.0/) + kk(1:3)
            klist(1:3,ii,jj) =  klist(1:3,ii,jj)/sqrt(sum( klist(1:3,ii,jj)**2))
        end do
    end do
end if

!===================================================================
! ------ perform interpolation from square lambert map
!===================================================================

scl = float(ipar(3))
quat = Quaternion_T( q = (/ 1.0, 0.0, 0.0, 0.0 /) )

do ip=1,ipar(4)
  do ii = imin, imax
    do jj = jmin, jmax

        dc(1:3) = klist(1:3,ii,jj)
        call quat%set_quats( quats(1:4,ip) )
        dc = quat%quat_Lp(dc)
        dc = dc/sqrt(sum(dc*dc))

        call LambertgetInterpolation(dc, scl, int(ipar(3)), int(ipar(3)), nix, niy, nixp, niyp, dx, dy, dxm, dym)

        if (dc(3).gt.0.0) then 
            Kosselpattern(ii,jj,ip) =  mLPNH(nix,niy,ipar(6)) * dxm * dym + &
                         mLPNH(nixp,niy,ipar(6)) * dx * dym + &
                         mLPNH(nix,niyp,ipar(6)) * dxm * dy + &
                         mLPNH(nixp,niyp,ipar(6)) * dx * dy 

        else
            Kosselpattern(ii,jj,ip) =  mLPSH(nix,niy,ipar(6)) * dxm * dym + &
                         mLPSH(nixp,niy,ipar(6)) * dx * dym + &
                         mLPSH(nix,niyp,ipar(6)) * dxm * dy + &
                         mLPSH(nixp,niyp,ipar(6)) * dx * dy 
        end if
    end do
  end do
end do

end subroutine getKosselPatterns




    
end module mod_SEMwrappers
    
