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

module mod_wrappers
  !! author: MDG 
  !! version: 1.0 
  !! date: 04/14/21
  !!
  !! wrappers module; contains wrappers to other EMsoftOO routines so that 
  !! IDL programs (among others) can call them 

use mod_kinds
use mod_global

IMPLICIT NONE 

contains

!--------------------------------------------------------------------------
recursive function getEBSDPatternsWrapper(argc, argv) bind(c, name='getEBSDPatternsWrapper') 
!DEC$ ATTRIBUTES DLLEXPORT :: getEBSDPatternsWrapper
  !! author: MDG
  !! version: 1.0
  !! date: 04/14/21
  !!
  !! wrapper routine for EMsoftCgetEBSDPatterns 
  !!
  !! see example at https://groups.google.com/forum/#!topic/comp.lang.idl-pvwave/Gk0xxVFbW8E

use,INTRINSIC :: ISO_C_BINDING
use mod_SEMwrappers

IMPLICIT NONE

INTEGER(c_size_t), VALUE, INTENT(IN)            :: argc 
type(c_ptr), dimension(argc), INTENT(INOUT)     :: argv
!f2py intent(in,out) ::  argv
REAL(c_float)                                   :: getEBSDPatternsWrapper

! wrapper function dependent declarations; they are all pointers 
! since we pass everything by reference from IDL 
integer(c_int32_t)                              :: nipar, nfpar, nq
integer(c_int32_t),dimension(:), pointer        :: ipar
real(c_float), dimension(:), pointer            :: fpar
real(c_float), dimension(:,:), pointer          :: quats
real(c_float), dimension(:,:,:), pointer        :: EBSDpattern 
integer(c_int32_t),dimension(:,:,:), pointer    :: accum_e 
real(c_float), dimension(:,:,:,:),pointer       :: mLPNH, mLPSH

TYPE(C_FUNPTR)                                  :: cproc
integer(c_size_t)                               :: objAddress
character(len=1)                                :: cancel

! the following line just helps in identifying the correct order of the subroutine arguments...
!                             1      2      3           4         5       6     7
!subroutine getEBSDPatterns(ipar, fpar, EBSDpattern, quats, accum_e, mLPNH, mLPSH)
!

! the ipar and fpar arrays are defined in the calling routine and are expected
! to have the following entries (in this case from the IDL routine EBSDExecute.pro)
! ipar = lonarr(wraparraysize)   wraparraysize is set to 80 in EMsoftOO 6.0
! fpar = fltarr(wraparraysize)   
! 
! the following parameters are used in this routine 
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

! transform the C pointers above to fortran pointers, and use them in the regular function call
nipar = 80
nfpar = 80
nq = 4

call c_f_pointer(argv(1),ipar,(/nipar/)) 
call c_f_pointer(argv(2),fpar,(/nfpar/)) 
call c_f_pointer(argv(3),EBSDpattern,(/ipar(2),ipar(3),ipar(8)/))
call c_f_pointer(argv(4),quats,(/nq,ipar(8)/))
call c_f_pointer(argv(5),accum_e,(/ipar(4),2*ipar(5)+1,2*ipar(5)+1/))
call c_f_pointer(argv(6),mLPNH,(/2*ipar(6)+1, 2*ipar(6)+1, ipar(4), ipar(7)/))
call c_f_pointer(argv(7),mLPSH,(/2*ipar(6)+1, 2*ipar(6)+1, ipar(4), ipar(7)/))

objAddress = 0
cancel = char(0)
call EMsoftCgetEBSDPatterns(ipar, fpar, EBSDpattern, quats, accum_e, mLPNH, mLPSH, C_NULL_FUNPTR, objAddress, cancel)

getEBSDPatternsWrapper = 1._c_float
end function getEBSDPatternsWrapper


!--------------------------------------------------------------------------
!
! SUBROUTINE:getECPatternsWrapper
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief wrapper routine for SingleECPPattern; based on Marc's routine above
!>
!> see https://groups.google.com/forum/#!topic/comp.lang.idl-pvwave/Gk0xxVFbW8E
!>
!
!> @param argc number of argument
!> @param argv pointers to subroutine parameters
!
!> @date 10/28/15  SS 1.0 original
!> @date 11/02/15 MDG 1.1 simplified parameters
!--------------------------------------------------------------------------
recursive function getECPatternsWrapper(argc, argv) bind(c, name='getECPatternsWrapper') 
!DEC$ ATTRIBUTES DLLEXPORT :: getECPatternsWrapper
  !! author: MDG
  !! version: 1.0
  !! date: 04/14/21
  !!
  !! wrapper routine for EMsoftCgetECPatterns 
  !!
  !! see example at https://groups.google.com/forum/#!topic/comp.lang.idl-pvwave/Gk0xxVFbW8E

use,INTRINSIC :: ISO_C_BINDING
use mod_SEMwrappers

IMPLICIT NONE

INTEGER(c_size_t), VALUE, INTENT(IN)            :: argc 
type(c_ptr), dimension(argc), INTENT(INOUT)     :: argv
!f2py intent(in,out) ::  argv
REAL(c_float)                                   :: getECPatternsWrapper

! wrapper function dependent declarations; they are all pointers 
! since we pass everything by reference from IDL 
integer(c_int32_t)                              :: nipar, nfpar, nq
integer(c_int32_t),dimension(:), pointer        :: ipar
real(c_float), dimension(:), pointer            :: fpar
real(c_float), dimension(:,:,:), pointer        :: accum_e 
real(c_float), dimension(:,:,:), pointer        :: mLPNH, mLPSH, ECPattern
real(c_float), dimension(:,:), pointer          :: quats
TYPE(C_FUNPTR)                                  :: cproc
integer(c_size_t)                               :: objAddress
character(len=1)                                :: cancel

! the following line just helps in identifying the correct order of the subroutine arguments...
!                             1      2     3       4       5       6       7
!subroutine getECPatterns(ipar, fpar, ECPattern, quats, accum_e, mLPNH, mLPSH)
!
! transform the C pointers above to fortran pointers, and use them in the regular function call
nipar = 8
nfpar = 8
nq = 4
call c_f_pointer(argv(1),ipar,(/nipar/)) 
call c_f_pointer(argv(2),fpar,(/nfpar/)) 
call c_f_pointer(argv(3),ECpattern,(/ipar(2),ipar(3),ipar(8)/))
call c_f_pointer(argv(4),quats,(/nq,ipar(8)/))
call c_f_pointer(argv(5),accum_e,(/ipar(4),2*ipar(5)+1,2*ipar(5)+1/))
call c_f_pointer(argv(6),mLPNH,(/2*ipar(7)+1, 2*ipar(7)+1, ipar(6)/))
call c_f_pointer(argv(7),mLPSH,(/2*ipar(7)+1, 2*ipar(7)+1, ipar(6)/))

objAddress = 0
cancel = char(0)
call EMsoftCgetECPatterns(ipar, fpar, ECpattern, quats, accum_e, mLPNH, mLPSH, C_NULL_FUNPTR, objAddress, cancel)

getECPatternsWrapper = 1._c_float
end function getECPatternsWrapper

end module mod_wrappers