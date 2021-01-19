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
! EMsoft:EMSEMwrappermod.f90
!--------------------------------------------------------------------------
!
! MODULE: EMSEMwrappermod
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
module EMSEMwrappermod

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
    
end module EMSEMwrappermod
    