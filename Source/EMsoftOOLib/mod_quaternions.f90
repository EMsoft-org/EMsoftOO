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

module mod_quaternions
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/03/20
  !!
  !! Quaternion and Quaternion Array arithmetic class 
  !!
  !! Quaternions are defined with the scalar part in position 1, and the vector part in positions 2:4.
  !!
  !! There are two class definitions in this file, one for single quaternions, the other for 
  !! quaternion array operations (some using OpenMP threads). The program MODQuaternionsTest.f90 
  !! can be used as part of ctest to run a test program on this module.


use mod_global
use mod_kinds
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit

IMPLICIT NONE
  private

! following are used to define the quaternion symmetry operators
    real(kind=dbl), public, parameter  :: sq22=0.7071067811865475244D0 ! sqrt(2)/2
    real(kind=dbl), public, parameter  :: sq32=0.8660254037844386467D0 ! sqrt(3)/2
    real(kind=dbl), public, parameter  :: half=0.5D0                   ! 1/2

! We define the rotational crystal symmetry operators in terms of quaternions (q0, q1,q2,q3) with q0 the scalar part;
! these are used in the dictmod EBSD dictionary indexing module, and are defined with respect to the standard cartesian 
! reference frame.  Note that these are not defined as Quaternion_T type, just as arrays of 4-component doubles.
    real(kind=dbl), public, dimension(4,152) :: SYM_Qsymop = reshape( (/ &
                      1.D0, 0.D0, 0.D0, 0.D0, &       ! 1: identity operator
                      0.D0, 1.D0, 0.D0, 0.D0, &       ! 2: 180@[100]
                      0.D0, 0.D0, 1.D0, 0.D0, &       ! 3: 180@[010]
                      0.D0, 0.D0, 0.D0, 1.D0, &       ! 4: 180@[001]
                      sq22, sq22, 0.D0, 0.D0, &       ! 5: 90@[100]
                      sq22, 0.D0, sq22, 0.D0, &       ! 6: 90@[010]
                      sq22, 0.D0, 0.D0, sq22, &       ! 7: 90@[001]
                      sq22,-sq22, 0.D0, 0.D0, &       ! 8: 270@[100]
                      sq22, 0.D0,-sq22, 0.D0, &       ! 9: 270@[010]
                      sq22, 0.D0, 0.D0,-sq22, &       !10: 270@[001]
                      0.D0, sq22, sq22, 0.D0, &       !11: 180@[110]
                      0.D0,-sq22, sq22, 0.D0, &       !12: 180@[-110]
                      0.D0, 0.D0, sq22, sq22, &       !13: 180@[011]
                      0.D0, 0.D0,-sq22, sq22, &       !14: 180@[0-11]
                      0.D0, sq22, 0.D0, sq22, &       !15: 180@[101]
                      0.D0,-sq22, 0.D0, sq22, &       !16: 180@[-101]
                      half, half, half, half, &       !17: 120@[111]
                      half,-half,-half,-half, &       !18: 120@[-1-1-1]
                      half, half,-half, half, &       !19: 120@[1-11]
                      half,-half, half,-half, &       !20: 120@[-11-1]
                      half,-half, half, half, &       !21: 120@[-111]
                      half, half,-half,-half, &       !22: 120@[1-1-1]
                      half,-half,-half, half, &       !23: 120@[-1-11]
                      half, half, half,-half, &       !24: 120@[11-1]
                      sq32, 0.D0, 0.D0, half, &       !25:  60@[001] (hexagonal/trigonal operators start here)
                      half, 0.D0, 0.D0, sq32, &       !26: 120@[001]
                      0.D0, 0.D0, 0.D0, 1.D0, &       !27: 180@[001] (duplicate from above, but useful to keep it here)
                     -half, 0.D0, 0.D0, sq32, &       !28: 240@[001]
                     -sq32, 0.D0, 0.D0, half, &       !29: 300@[001]
                      0.D0, 1.D0, 0.D0, 0.D0, &       !30: 180@[100]
                      0.D0, sq32, half, 0.D0, &       !31: 180@[xxx]
                      0.D0, half, sq32, 0.D0, &       !32: 180@[xxx]
                      0.D0, 0.D0, 1.D0, 0.D0, &       !33: 180@[010]
                      0.D0,-half, sq32, 0.D0, &       !34: 180@[xxx]
                      0.D0,-sq32, half, 0.D0, &       !35: 180@[xxx]
                      1.0000000000000000D0, 0.00000000000000000D0, 0.00000000000000000D0, 0.00000000000000000D0, & ! icosahedral operators
                      0.0000000000000000D0, 0.68819093704223555D0, 0.50000000000000000D0, 0.52573108673095648D0, &
                      0.0000000000000000D0,-0.26286554336547824D0, 0.80901700258254916D0, 0.52573108673095648D0, &
                      0.0000000000000000D0,-0.85065078735351463D0, 0.00000000000000000D0, 0.52573108673095648D0, &
                      0.0000000000000000D0,-0.26286554336547824D0,-0.80901700258254916D0, 0.52573108673095648D0, &
                      0.0000000000000000D0, 0.68819093704223555D0,-0.50000000000000000D0, 0.52573108673095648D0, &
                      0.0000000000000000D0, 0.52573108673095648D0, 0.00000000000000000D0, 0.85065078735351463D0, &
                      0.0000000000000000D0, 0.16245985031127913D0, 0.50000000000000000D0, 0.85065078735351463D0, &
                      0.0000000000000000D0,-0.42532539367675731D0, 0.30901700258254972D0, 0.85065078735351463D0, &
                      0.0000000000000000D0,-0.42532539367675731D0,-0.30901700258254972D0, 0.85065078735351463D0, &
                      0.0000000000000000D0, 0.16245985031127913D0,-0.50000000000000000D0, 0.85065078735351463D0, &
                      0.0000000000000000D0, 0.95105654001235851D0,-0.30901700258254972D0, 0.00000000000000000D0, &
                      0.0000000000000000D0, 0.95105654001235851D0, 0.30901700258254972D0, 0.00000000000000000D0, &
                      0.0000000000000000D0, 0.58778524398803644D0, 0.80901700258254916D0, 0.00000000000000000D0, &
                      0.0000000000000000D0, 0.00000000000000000D0, 1.00000000000000000D0, 0.00000000000000000D0, &
                      0.0000000000000000D0,-0.58778524398803644D0, 0.80901700258254916D0, 0.00000000000000000D0, &
                     0.50000000000000000D0, 0.42532540417601997D0, 0.30901699437494742D0, 0.68819096193209561D0, &
                     0.50000000000000000D0,-0.42532540417601997D0,-0.30901699437494742D0,-0.68819096193209561D0, &
                     0.50000000000000000D0,-0.16245984737382999D0, 0.50000000000000000D0, 0.68819096193209561D0, &
                     0.50000000000000000D0, 0.16245984737382999D0,-0.50000000000000000D0,-0.68819096193209561D0, &
                     0.50000000000000000D0,-0.52573108874869778D0, 0.00000000000000000D0, 0.68819096193209561D0, &
                     0.50000000000000000D0, 0.52573108874869778D0, 0.00000000000000000D0,-0.68819096193209561D0, &
                     0.50000000000000000D0,-0.16245984737382999D0,-0.50000000000000000D0, 0.68819096193209561D0, &
                     0.50000000000000000D0, 0.16245984737382999D0, 0.50000000000000000D0,-0.68819096193209561D0, &
                     0.50000000000000000D0, 0.42532539174817890D0,-0.30901700049082287D0, 0.68819096193209561D0, &
                     0.50000000000000000D0,-0.42532539174817890D0, 0.30901700049082287D0,-0.68819096193209561D0, &
                     0.50000000000000000D0,-0.85065078349635781D0, 0.00000000000000000D0, 0.16245984737382999D0, &
                     0.50000000000000000D0, 0.85065078349635781D0, 0.00000000000000000D0,-0.16245984737382999D0, &
                     0.50000000000000000D0,-0.68819096193209561D0,-0.50000000000000000D0,-0.16245984737382999D0, &
                     0.50000000000000000D0, 0.68819096193209561D0, 0.50000000000000000D0, 0.16245984737382999D0, &
                     0.50000000000000000D0,-0.26286554437434889D0,-0.80901695670145712D0, 0.16245984737382999D0, &
                     0.50000000000000000D0, 0.26286554437434889D0, 0.80901695670145712D0,-0.16245984737382999D0, &
                     0.50000000000000000D0, 0.26286554437434889D0,-0.80901695670145712D0,-0.16245984737382999D0, &
                     0.50000000000000000D0,-0.26286554437434889D0, 0.80901695670145712D0, 0.16245984737382999D0, &
                     0.50000000000000000D0, 0.68819096193209561D0,-0.50000000000000000D0, 0.16245984737382999D0, &
                     0.50000000000000000D0,-0.68819096193209561D0, 0.50000000000000000D0,-0.16245984737382999D0, &
                     0.80901700537708732D0, 0.00000000000000000D0, 0.00000000000000000D0, 0.58778523714932640D0, &
                     0.30901702997862029D0, 0.00000000000000000D0, 0.00000000000000000D0, 0.95105650472681824D0, &
                     0.80901700537708732D0, 0.00000000000000000D0, 0.00000000000000000D0,-0.58778523714932640D0, &
                     0.30901702997862029D0, 0.00000000000000000D0, 0.00000000000000000D0,-0.95105650472681824D0, &
                     0.80901700537708732D0, 0.52573109227969150D0, 0.00000000000000000D0, 0.26286554613984575D0, &
                     0.30901702997862029D0, 0.85065078781948245D0, 0.00000000000000000D0, 0.42532539390974122D0, &
                     0.80901700537708732D0,-0.52573109227969150D0, 0.00000000000000000D0,-0.26286554613984575D0, &
                     0.30901702997862029D0,-0.85065078781948245D0, 0.00000000000000000D0,-0.42532539390974122D0, &
                     0.80901700537708732D0, 0.16245984550474032D0, 0.50000000000000000D0, 0.26286554613984575D0, &
                     0.30901702997862029D0, 0.26286555540853851D0, 0.80901696456355054D0, 0.42532539390974122D0, &
                     0.80901700537708732D0,-0.16245984550474032D0,-0.50000000000000000D0,-0.26286554613984575D0, &
                     0.30901702997862029D0,-0.26286555540853851D0,-0.80901696456355054D0,-0.42532539390974122D0, &
                     0.80901700537708732D0,-0.42532540916195122D0, 0.30901697149092866D0, 0.26286554613984575D0, &
                     0.30901702997862029D0,-0.68819097766197224D0, 0.50000000000000000D0, 0.42532539390974122D0, &
                     0.80901700537708732D0, 0.42532540916195122D0,-0.30901697149092866D0,-0.26286554613984575D0, &
                     0.30901702997862029D0, 0.68819097766197224D0,-0.50000000000000000D0,-0.42532539390974122D0, &
                     0.80901700537708732D0,-0.42532540916195122D0,-0.30901697149092866D0, 0.26286554613984575D0, &
                     0.30901702997862029D0,-0.68819097766197224D0,-0.50000000000000000D0, 0.42532539390974122D0, &
                     0.80901700537708732D0, 0.42532540916195122D0, 0.30901697149092866D0,-0.26286554613984575D0, &
                     0.30901702997862029D0, 0.68819097766197224D0, 0.50000000000000000D0,-0.42532539390974122D0, &
                     0.80901700537708732D0, 0.16245984550474032D0,-0.50000000000000000D0, 0.26286554613984575D0, &
                     0.30901702997862029D0, 0.26286555540853851D0,-0.80901696456355054D0, 0.42532539390974122D0, &
                     0.80901700537708732D0,-0.16245984550474032D0, 0.50000000000000000D0,-0.26286554613984575D0, &
                     0.30901702997862029D0,-0.26286555540853851D0, 0.80901696456355054D0,-0.42532539390974122D0, & 

                     ! octagonal QC group 822                              
                     0.923879532511287D0, 0.D0, 0.D0, 0.38268343236509D0, &       ! 45@[001]
                     0.707106781186547D0, 0.D0, 0.D0, 0.707106781186547D0, &      ! 90@[001]
                     0.38268343236509D0, 0.D0, 0.D0, 0.923879532511287D0, &       ! 135@[001]
                     0.D0, 0.D0, 0.D0, 1.D0, &                                    ! 180@[001]
                    -0.382683432365090D0, 0.D0, 0.D0, 0.923879532511287D0, &     ! 225@[001]
                    -0.707106781186547D0, 0.0D0, 0.0D0, 0.707106781186548D0, &   ! 270@[001]
                    -0.923879532511287D0, 0.0D0, 0.0D0, 0.382683432365090D0, &   ! 315@[001]
                     ! all 2 fold rotation axes
                     0.D0, 1.D0, 0.D0, 0.D0, &                                    ! 180@[100]
                     0.D0, 0.923879532511287D0, 0.38268343236509D0, 0.D0, &       ! 180@[cos(pi/8) sin(pi/8) 0]
                     0.0D0, 0.923879532511287D0, -0.382683432365090D0, 0.0D0, &
                     0.0D0, 0.707106781186548D0, -0.707106781186547D0, 0.0D0, &
                     0.0D0, 0.382683432365090D0, -0.923879532511287D0, 0.0D0, &
                     0.0D0, 0.0D0, -1.0D0, 0.0D0, &
                     0.0D0, -0.382683432365090D0, -0.923879532511287D0, 0.0D0, &
                     0.0D0, -0.707106781186547D0, -0.707106781186548D0, 0.0D0, &
                     
                     ! decagonal QC group 1022
                     0.951056516295154D0, 0.0D0, 0.0D0, 0.309016994374947D0, &    ! 36@[001]
                     0.809016994374947D0, 0.0D0, 0.0D0, 0.587785252292473D0, &    ! 72@[001]
                     0.587785252292473D0, 0.0D0, 0.0D0, 0.809016994374947D0, &    ! 108@[001]
                     0.309016994374947D0, 0.0D0, 0.0D0, 0.951056516295154D0, &    ! 144@[001]
                     0.0D0, 0.0D0, 0.0D0, 1.0D0, &                                ! 180@[001]
                    -0.309016994374947D0, 0.0D0, 0.0D0, 0.951056516295154D0, &    ! 216@[001]
                    -0.587785252292473D0, 0.0D0, 0.0D0, 0.809016994374947D0, &    ! 252@[001]
                    -0.809016994374947D0, 0.0D0, 0.0D0, 0.587785252292473D0, &    ! 288@[001]
                    -0.951056516295154D0, 0.0D0, 0.0D0, 0.309016994374948D0, &    ! 324@[001]  
                    ! all 2-fold rotation axis
                     0.D0, 1.D0, 0.D0, 0.D0, &                                    ! 180@[100]
                     0.D0, 0.951056516295154D0, 0.309016994374947D0, 0.D0,   &    ! 180@[cos(pi/10) sin(pi/10) 0]
                     0.0D0, 0.951056516295154D0, -0.309016994374947D0, 0.0D0, &
                     0.0D0, 0.809016994374947D0, -0.587785252292473D0, 0.0D0, &
                     0.0D0, 0.587785252292473D0, -0.809016994374947D0, 0.0D0, &
                     0.0D0, 0.309016994374947D0, -0.951056516295154D0, 0.0D0, &
                     0.0D0, 0.0D0, -1.0D0, 0.0D0, &
                     0.0D0, -0.309016994374947D0, -0.951056516295154D0, 0.0D0, &
                     0.0D0, -0.587785252292473D0, -0.809016994374947D0, 0.0D0, &
                     0.0D0, -0.809016994374947D0, -0.587785252292473D0, 0.0D0, &

                     ! dodecagonal QC group 1222
                     0.965925826289068D0, 0.0D0, 0.0D0, 0.258819045102521D0, &    ! 30@[001] 
                     0.866025403784439D0, 0.0D0, 0.0D0, 0.5D0, &                  ! 60@[001] 
                     0.707106781186548D0, 0.0D0, 0.0D0, 0.707106781186547D0, &    ! 90@[001] 
                     0.5D0, 0.0D0, 0.0D0, 0.866025403784439D0, &                  ! 120@[001] 
                     0.258819045102521D0, 0.0D0, 0.0D0, 0.965925826289068D0, &    ! 150@[001] 
                     0.0D0, 0.0D0, 0.0D0, 1.0D0, &                                    ! 180@[001] 
                    -0.258819045102521D0, 0.0D0, 0.0D0, 0.965925826289068D0, &      ! 210@[001] 
                    -0.5D0, 0.0D0, 0.0D0, 0.866025403784439D0, &                  ! 240@[001] 
                    -0.707106781186547D0, 0.0D0, 0.0D0, 0.707106781186548D0, &    ! 270@[001] 
                    -0.866025403784439D0, 0.0D0, 0.0D0, 0.5D0, &                  ! 300@[001] 
                    -0.965925826289068D0, 0.0D0, 0.0D0, 0.258819045102521D0, &    ! 330@[001]
                    ! all 2-fold rotation axes
                    0.0D0, 1.0D0, 0.0D0, 0.0D0, &                                    ! 180@[100]
                    0.0D0, 0.965925826289068D0, 0.258819045102521D0, 0.0D0, &      ! 180@[cos(pi/12) sin(pi/12) 0]
                    0.0D0, 0.965925826289068D0, -0.258819045102521D0, 0.0D0, &
                    0.0D0, 0.866025403784439D0, -0.5D0, 0.0D0, &
                    0.0D0, 0.707106781186548D0, -0.707106781186547D0, 0.0D0, &
                    0.0D0, 0.5D0, -0.866025403784439D0, 0.0D0, &
                    0.0D0, 0.258819045102521D0, -0.965925826289068D0, 0.0D0, &
                    0.0D0, 0.0D0, -1.0D0, 0.0D0, &
                    0.0D0, -0.258819045102521D0, -0.965925826289068D0, 0.0D0, &
                    0.0D0, -0.5D0, -0.866025403784439D0, 0.0D0, &
                    0.0D0, -0.707106781186547D0, -0.707106781186548D0, 0.0D0, &
                    0.0D0, -0.866025403784439D0, -0.5D0, 0.0D0 &
                    /), (/4,152/) )
!DEC$ ATTRIBUTES DLLEXPORT :: SYM_Qsymop
 
! we overload the conjg, cabs, and .eq. intrinsics
    intrinsic :: conjg, cabs
    public :: conjg, cabs

    interface conjg
      procedure quatconjg
      procedure quatarrayconjg
    end interface conjg

    interface cabs
      procedure quatnorm
      procedure quatarraynorm
    end interface cabs

    interface operator(.eq.)
      procedure quatsequal
    end interface 

! definition of the quaternion class 
  type, public :: Quaternion_T
    !! Quaternion Class definition
    private
      real(kind=sgl), dimension(4) :: q
       !! single precision quaternion
      real(kind=dbl), dimension(4) :: qd
       !! double precision quaternion
      character(1)                 :: s
       !! precision indicator ('s' or 'd')

    contains
    private 
! quaternion IO routines
      procedure, pass(self) :: quatprint
      procedure, pass(self) :: getquats
      procedure, pass(self) :: getquatd
      procedure, pass(self) :: setquats
      procedure, pass(self) :: setquatd
! quaternion arithmetic routines 
      procedure, pass(self) :: quatflip
      procedure, pass(self) :: quatpos
      procedure, pass(self) :: quatadd
      procedure, pass(self) :: quatsubtract
      procedure, pass(self) :: quatmult
      procedure, pass(self) :: quatsmult
      procedure, pass(self) :: quatsmultd
      procedure, pass(self) :: quatdiv
      procedure, pass(self) :: quatsdiv
      procedure, pass(self) :: quatsdivd
      procedure, pass(self) :: quatconjg
      procedure, pass(self) :: quatnorm
      procedure, pass(self) :: quatnormalize
! quaternion-based transformations
      procedure, pass(self) :: quatLp
      procedure, pass(self) :: quatLpd
! routines with two or more input quaternions
      procedure, pass(self) :: quatinnerproduct
      procedure, pass(self) :: quatangle
      procedure, pass(self) :: quatslerp
! miscellaneous routines 
      procedure, pass(self), public :: quatsequal

      generic, public :: quat_print => quatprint
      generic, public :: quat_flip => quatflip
      generic, public :: quat_pos => quatpos
      generic, public :: get_quats => getquats
      generic, public :: get_quatd => getquatd
      generic, public :: set_quats => setquats
      generic, public :: set_quatd => setquatd
      generic, public :: quat_norm => quatnorm
      generic, public :: operator(+) => quatadd
      generic, public :: operator(-) => quatsubtract
      generic, public :: operator(*) => quatmult 
      generic, public :: operator(*) => quatsmult, quatsmultd
      generic, public :: operator(/) => quatdiv 
      generic, public :: operator(/) => quatsdiv, quatsdivd
      generic, public :: quat_normalize => quatnormalize
      generic, public :: quat_Lp => quatLp, quatLpd
      generic, public :: quat_innerproduct => quatinnerproduct
      generic, public :: quat_angle => quatangle
      generic, public :: quat_slerp => quatslerp

  end type Quaternion_T 

! the constructor routine for this class 
  interface Quaternion_T
    module procedure Quaternion_constructor
  end interface Quaternion_T



! next we define the quaternion array class 
  type, public :: QuaternionArray_T
    !! Quaternion Class definition
    private
      integer(kind=irg)            :: n 
      integer(kind=irg)            :: nthreads
      real(kind=sgl), allocatable  :: q(:,:)
       !! single precision quaternion
      real(kind=dbl), allocatable  :: qd(:,:)
       !! double precision quaternion
      character(1)                 :: s
       !! precision indicator ('s' or 'd')

    contains
    private 
    ! quaternion IO routines
      procedure, pass(self) :: quatarrayprint
! quaternion arithmetic routines 
      procedure, pass(self) :: quatarrayadd
      procedure, pass(self) :: quatarraysubtract
      procedure, pass(self) :: quatarraymult
      procedure, pass(self) :: quatarraysmult
      procedure, pass(self) :: quatarraysmultd
      procedure, pass(self) :: quatarraydiv
      procedure, pass(self) :: quatarraysdiv
      procedure, pass(self) :: quatarrayconjg
      procedure, pass(self) :: quatarraynorm
      procedure, pass(self) :: quatarraynormalize
! quaternion-based transformations
      procedure, pass(self) :: quatarrayLp
      procedure, pass(self) :: quatarrayLpd
! routines with two or more input quaternion arrays
      procedure, pass(self) :: quatarrayinnerproduct
      procedure, pass(self) :: quatarrayangle
! miscellaneous routines 
      procedure, pass(self) :: extractfromQuaternionArray
      procedure, pass(self) :: insertQuatintoArray
      procedure, pass(self) :: QSym_Init_
      procedure, pass(self) :: getQnumber_

! generics
      generic, public :: quat_print => quatarrayprint
      generic, public :: operator(+) => quatarrayadd
      generic, public :: operator(-) => quatarraysubtract
      generic, public :: operator(*) => quatarraymult 
      generic, public :: operator(*) => quatarraysmult, quatarraysmultd
      generic, public :: operator(/) => quatarraydiv 
      generic, public :: operator(/) => quatarraysdiv
      generic, public :: quat_normalize => quatarraynormalize
      generic, public :: quat_Lp => quatarrayLp, quatarrayLpd
      generic, public :: quat_innerproduct => quatarrayinnerproduct
      generic, public :: quat_angle => quatarrayangle
      generic, public :: getQuatfromArray => extractfromQuaternionArray
      generic, public :: insertQuatinArray => insertQuatintoArray
      generic, public :: QSym_Init => QSym_Init_
      generic, public :: getQnumber => getQnumber_

  end type QuaternionArray_T 

! the constructor routine for this class 
  interface QuaternionArray_T
    module procedure QuaternionArray_constructor
  end interface QuaternionArray_T

  public :: quat_randomArray
  interface quat_randomArray 
    module procedure generateRandomArray
  end interface 

  private :: quat_Marsaglia 

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! We begin with the functions/subroutines that are public in the 
! two classes and pair up functions for individual and quaternion 
! arrays for easier module maintenance.
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
type(Quaternion_T) function Quaternion_constructor( q, qd ) result(Quat)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! constructor for the Quaternion Class 
  
IMPLICIT NONE

  real(kind=sgl), INTENT(IN), OPTIONAL      :: q(4)
  real(kind=dbl), INTENT(IN), OPTIONAL      :: qd(4)

! fill in one or the other quaternion 
  if ((.not.present(q)).and.(.not.present(qd))) then 
    Quat % q = (/ 0.0, 0.0, 0.0, 0.0 /)
    Quat % qd = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
    Quat % s = 's'
  else 
      if (present(q)) then 
        Quat % q = q 
        Quat % qd = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        Quat % s = 's'
      end if 

      if (present(qd)) then 
        Quat % q = (/ 0.0, 0.0, 0.0, 0.0 /)
        Quat % qd = qd 
        Quat % s = 'd'
      end if 
  end if

end function Quaternion_constructor

!--------------------------------------------------------------------------
subroutine Quaternion_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 02/02/20
!!
!! destructor for the Quaternion_T Class

IMPLICIT NONE

type(Quaternion_T), INTENT(INOUT)     :: self 

call reportDestructor('Quaternion_T')

end subroutine Quaternion_destructor

!--------------------------------------------------------------------------
type(QuaternionArray_T) function QuaternionArray_constructor( n, nthreads, q, qd, s ) result(QuatArray)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! constructor for the QuaternionArray Class 
  !!
  !! either call with parameters n and s 
  !! or with n and either one of q or qd 
  
use mod_io 

IMPLICIT NONE

  integer(kind=irg), INTENT(IN)             :: n
  integer(kind=irg), INTENT(IN), OPTIONAL   :: nthreads
  real(kind=sgl), INTENT(IN), OPTIONAL      :: q(4,n)
  real(kind=dbl), INTENT(IN), OPTIONAL      :: qd(4,n)
  character(1), INTENT(IN), OPTIONAL        :: s

  type(IO_T)                                :: Message

! OpenMP threads
  QuatArray % nthreads = 0
  if (present(nthreads)) QuatArray % nthreads = nthreads

! are we declaring just an empty variable with no entries, but with a given precision ?
  if ( present(s) .and. (.not.present(q)) .and. (.not.present(qd)) ) then 
    QuatArray % n = n 
    QuatArray % s = s 
    if (s.eq.'s') then 
      allocate(QuatArray % q(4,n))  
      QuatArray % q = 0.0
    else 
      allocate(QuatArray % qd(4,n))  
      QuatArray % qd = 0.D0
    end if 
    return 
  end if  
    
! number of quaternions in array
  QuatArray % n = n 

! single precision
  if (present(q)) then 
    allocate(QuatArray % q(4,n))  
    QuatArray % q = q
    QuatArray % s = 's'
  end if 

! double precision
  if (present(qd)) then 
    allocate(QuatArray % qd(4,n))  
    QuatArray % qd = qd
    QuatArray % s = 'd'
  end if 

end function QuaternionArray_constructor

!--------------------------------------------------------------------------
subroutine QuaternionArray_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 02/02/20
!!
!! destructor for the QuaternionArray_T Class

IMPLICIT NONE

type(QuaternionArray_T), INTENT(INOUT)     :: self 

call reportDestructor('QuaternionArray_T')

if (allocated(self%q)) deallocate(self%q)
if (allocated(self%qd)) deallocate(self%qd)

end subroutine QuaternionArray_destructor

!--------------------------------------------------------------------------
recursive subroutine quatprint(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! print a quaternion

use mod_io

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion 

  type(IO_T)                        :: Message 

  if (self%s.eq.'s') then 
    call Message % WriteValue('', self%q, 4, frm="('(',4f12.6,'); precision: '$)")
    call Message % WriteValue('',self%s)
  else 
    call Message % WriteValue('', self%qd, 4, frm="('(',4f20.14,'); precision: '$)")
    call Message % WriteValue('',self%s)
  end if 

end subroutine quatprint

!--------------------------------------------------------------------------
recursive function getquats(self) result(qs)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! return a quaternion

IMPLICIT NONE 

class(Quaternion_T),intent(in)    :: self
 !! input quaternion 
real(kind=sgl)                    :: qs(4) 

qs = self%q 

end function getquats

!--------------------------------------------------------------------------
recursive function getquatd(self) result(qd)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! return a quaternion

IMPLICIT NONE 

class(Quaternion_T),intent(in)    :: self
 !! input quaternion 
real(kind=dbl)                    :: qd(4) 

qd = self%qd 

end function getquatd

!--------------------------------------------------------------------------
recursive subroutine setquats(self, qs)
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/18/20
  !!
  !! set a quaternion

IMPLICIT NONE 

class(Quaternion_T),intent(inout)    :: self
real(kind=sgl),intent(in)            :: qs(4) 
 !! input quaternion 

self%q = qs
self%s = 's'

end subroutine setquats

!--------------------------------------------------------------------------
recursive subroutine setquatd(self, qd)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! set a quaternion

IMPLICIT NONE 

class(Quaternion_T),intent(inout)    :: self
real(kind=dbl),intent(in)            :: qd(4) 
 !! input quaternion 

self%qd = qd
self%s = 'd'

end subroutine setquatd

!--------------------------------------------------------------------------
recursive subroutine quatarrayprint(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! print an array of quaternions

use mod_io

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)   :: self
   !! input quaternion 

  type(IO_T)                            :: Message 
  integer(kind=irg)                     :: i

  if (self%s.eq.'s') then 
    do i=1,self%n
      call Message % WriteValue('', self%q(:,i), 4, frm="('(',4f12.6,')')")
    end do
  else 
    do i=1,self%n
      call Message % WriteValue('', self%qd(:,i), 4, frm="('(',4f20.14,')')")
    end do
  end if 

end subroutine quatarrayprint

!--------------------------------------------------------------------------
pure recursive subroutine quatflip(self) 
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! change the sign of the complete quaternion

IMPLICIT NONE

class(Quaternion_T),intent(inout) :: self

if (self%s.eq.'s') then 
  self%q = -self%q 
else 
  self%qd = -self%qd 
end if 

end subroutine quatflip 

!--------------------------------------------------------------------------
pure recursive subroutine quatpos(self) 
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! convert a quaternion to one with a positive scalar part, if it is negative

IMPLICIT NONE

class(Quaternion_T),intent(inout) :: self

if (self%s.eq.'s') then 
  if (self%q(1).lt.0.0) self%q = -self%q 
else 
  if (self%qd(1).lt.0.D0) self%qd = -self%qd 
end if 

end subroutine quatpos 

!--------------------------------------------------------------------------
pure recursive function quatadd(self, y) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion addition (single/double precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in) :: self, y
  type(Quaternion_T)             :: qres 

  if (self%s.eq.'s') then
    qres%q = self%q + y%q 
    qres%s = 's'
  else
    qres%qd = self%qd + y%qd 
    qres%s = 'd'
  end if 

end function quatadd

!--------------------------------------------------------------------------
recursive function quatarrayadd(self, y) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! quaternion array addition (single/double precision)

use mod_io

IMPLICIT NONE

  class(QuaternionArray_T),intent(in) :: self, y
  type(QuaternionArray_T)             :: qres 

  type(IO_T)                          :: Message 
  integer(kind=irg)                   :: sz(2)

! test to make sure that both arrays have the same number of quaternions 
  if (self%n.ne.y%n) then 
    call Message%printError('quatarrayadd','input arrays must have the same number of quaternions')
  end if 

  if (self%s.ne.y%s) then 
    call Message%printError('quatarrayadd','input arrays must have the same precision')
  end if 

  qres%n = self%n
  qres%s = self%s
  qres%nthreads = self%nthreads

  if (self%s.eq.'s') then
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%q)) then 
      sz = shape(qres%q)
      if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%q)
    end if 
    allocate(qres%q(4,self%n))
    qres%q = self%q + y%q 
  else
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%qd)) then 
      sz = shape(qres%qd)
      if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%qd)
    end if 
    allocate(qres%qd(4,self%n))
    qres%qd = self%qd + y%qd 
  end if 

end function quatarrayadd

!--------------------------------------------------------------------------
recursive function quatsubtract(self, y) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion subtraction (single/double precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in) :: self, y
  type(Quaternion_T)             :: qres 

  if (self%s.eq.'s') then 
    qres%q = self%q - y%q 
    qres%s = 's'
  else
    qres%qd = self%qd - y%qd 
    qres%s = 'd'
  end if

end function quatsubtract

!--------------------------------------------------------------------------
recursive function quatarraysubtract(self, y) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! quaternion array subtraction (single/double precision)

use mod_io 

IMPLICIT NONE

  class(QuaternionArray_T),intent(in) :: self, y
  type(QuaternionArray_T)             :: qres 

  type(IO_T)                          :: Message 
  integer(kind=irg)                   :: sz(2)

! test to make sure that both arrays have the same number of quaternions 
  if (self%n.ne.y%n) then 
    call Message%printError('quatarrayadd','input arrays must have the same number of quaternions')
  end if 

  if (self%s.ne.y%s) then 
    call Message%printError('quatarrayadd','input arrays must have the same precision')
  end if 

  qres%n = self%n
  qres%s = self%s
  qres%nthreads = self%nthreads

   if (self%s.eq.'s') then
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%q)) then 
      sz = shape(qres%q)
      if ( (sz(1).ne.4).or.(sz(2).ne.self%n) ) deallocate(qres%q)
    end if 
    allocate(qres%q(4,self%n))
    qres%q = self%q - y%q 
  else
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%qd)) then 
      sz = shape(qres%qd)
      if ( (sz(1).ne.4).or.(sz(2).ne.self%n) ) deallocate(qres%qd)
    end if 
    allocate(qres%qd(4,self%n))
    qres%qd = self%qd - y%qd 
  end if 

end function quatarraysubtract

!--------------------------------------------------------------------------
pure recursive function quatmult(self, y) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion multiplication   (single/double precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in) :: self, y
   !! input quaternions
  type(Quaternion_T)             :: qres 
   !! output quaternion 

! the following is a way to reduce the number of multiplications
! needs to be tested and merged with epsijk approach
!
! QuatMul(QUAT *q1, QUAT *q2, QUAT *res){
! float A, B, C, D, E, F, G, H;
! A = (q1->w + q1->x)*(q2->w + q2->x);
! B = (q1->z - q1->y)*(q2->y - q2->z);
! C = (q1->w - q1->x)*(q2->y + q2->z); 
! D = (q1->y + q1->z)*(q2->w - q2->x);
! E = (q1->x + q1->z)*(q2->x + q2->y);
! F = (q1->x - q1->z)*(q2->x - q2->y);
! G = (q1->w + q1->y)*(q2->w - q2->z);
! H = (q1->w - q1->y)*(q2->w + q2->z);
! res->w = B + (-E - F + G + H) /2;
! res->x = A - (E + F + G + H)/2; 
! res->y = C + (E - F + G - H)/2; 
! res->z = D + (E - F - G + H)/2;
! }

  if (self%s.eq.'s') then 
    qres%q = (/ self%q(1)*y%q(1) - self%q(2)*y%q(2) -          ( self%q(3)*y%q(3) + self%q(4)*y%q(4) ), &
                self%q(1)*y%q(2) + self%q(2)*y%q(1) + epsijk * ( self%q(3)*y%q(4) - self%q(4)*y%q(3) ), &
                self%q(1)*y%q(3) + self%q(3)*y%q(1) + epsijk * ( self%q(4)*y%q(2) - self%q(2)*y%q(4) ), &
                self%q(1)*y%q(4) + self%q(4)*y%q(1) + epsijk * ( self%q(2)*y%q(3) - self%q(3)*y%q(2) ) /)
    qres%s = 's'
  else 
    qres%qd = (/ self%qd(1)*y%qd(1) - self%qd(2)*y%qd(2) -           ( self%qd(3)*y%qd(3) + self%qd(4)*y%qd(4) ), &
                 self%qd(1)*y%qd(2) + self%qd(2)*y%qd(1) + epsijkd * ( self%qd(3)*y%qd(4) - self%qd(4)*y%qd(3) ), &
                 self%qd(1)*y%qd(3) + self%qd(3)*y%qd(1) + epsijkd * ( self%qd(4)*y%qd(2) - self%qd(2)*y%qd(4) ), &
                 self%qd(1)*y%qd(4) + self%qd(4)*y%qd(1) + epsijkd * ( self%qd(2)*y%qd(3) - self%qd(3)*y%qd(2) ) /) 
    qres%s = 'd'
  end if
      
end function quatmult

!--------------------------------------------------------------------------
recursive function quatarraymult(self, y) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! quaternion array multiplication   (single/double precision)

use mod_io
use omp_lib

IMPLICIT NONE

  class(QuaternionArray_T),intent(in) :: self, y
   !! input quaternion arrays
  type(QuaternionArray_T)             :: qres 
   !! output quaternion array

  type(IO_T)                          :: Message 
  integer(kind=irg)                   :: i, maxthreads, ompthreads, sz(2) 

! test to make sure that both arrays have the same number of quaternions 
  if (self%n.ne.y%n) then 
    call Message%printError('quatarrayadd','input arrays must have the same number of quaternions')
  end if 

  if (self%s.ne.y%s) then 
    call Message%printError('quatarrayadd','input arrays must have the same precision')
  end if 

  qres%n = self%n
  qres%s = self%s
  qres%nthreads = self%nthreads

! set the number of OpenMP threads
  maxthreads = OMP_GET_MAX_THREADS()
  if (self%nthreads.eq.0) then ! use the maximum number of threads unless self%n is less than that
    ompthreads = maxthreads
    if (self%n.lt.ompthreads) ompthreads = self%n 
  else
    ompthreads = self%nthreads
    if (self%nthreads.gt.maxthreads) ompthreads = maxthreads
  end if 
  call OMP_SET_NUM_THREADS(ompthreads)

  if (self%s.eq.'s') then 
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%q)) then 
      sz = shape(qres%q)
      if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%q)
    end if 
    allocate(qres%q(4,self%n))

!$OMP PARALLEL DEFAULT(shared)
!$OMP DO SCHEDULE(DYNAMIC,1) 
    do i=1,self%n
      qres%q(1:4,i) = (/ &
         self%q(1,i)*y%q(1,i) - self%q(2,i)*y%q(2,i) -          ( self%q(3,i)*y%q(3,i) + self%q(4,i)*y%q(4,i) ), &
         self%q(1,i)*y%q(2,i) + self%q(2,i)*y%q(1,i) + epsijk * ( self%q(3,i)*y%q(4,i) - self%q(4,i)*y%q(3,i) ), &
         self%q(1,i)*y%q(3,i) + self%q(3,i)*y%q(1,i) + epsijk * ( self%q(4,i)*y%q(2,i) - self%q(2,i)*y%q(4,i) ), &
         self%q(1,i)*y%q(4,i) + self%q(4,i)*y%q(1,i) + epsijk * ( self%q(2,i)*y%q(3,i) - self%q(3,i)*y%q(2,i) ) /)
    end do 
!$OMP END DO 
!$OMP END PARALLEL    
  else 
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%qd)) then 
      sz = shape(qres%qd)
      if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%qd)
    end if 
    allocate(qres%qd(4,self%n))

!$OMP PARALLEL DEFAULT(shared)
!$OMP DO SCHEDULE(DYNAMIC,1) 
    do i=1,self%n
    qres%qd(1:4,i) = (/ &
         self%qd(1,i)*y%qd(1,i) - self%qd(2,i)*y%qd(2,i) -           ( self%qd(3,i)*y%qd(3,i) + self%qd(4,i)*y%qd(4,i) ), &
         self%qd(1,i)*y%qd(2,i) + self%qd(2,i)*y%qd(1,i) + epsijkd * ( self%qd(3,i)*y%qd(4,i) - self%qd(4,i)*y%qd(3,i) ), &
         self%qd(1,i)*y%qd(3,i) + self%qd(3,i)*y%qd(1,i) + epsijkd * ( self%qd(4,i)*y%qd(2,i) - self%qd(2,i)*y%qd(4,i) ), &
         self%qd(1,i)*y%qd(4,i) + self%qd(4,i)*y%qd(1,i) + epsijkd * ( self%qd(2,i)*y%qd(3,i) - self%qd(3,i)*y%qd(2,i) ) /) 
    end do 
!$OMP END DO 
!$OMP END PARALLEL    
  end if
      
end function quatarraymult

!--------------------------------------------------------------------------
pure recursive function quatsmult(self, s) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! scalar quaternion multiplication   (single precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in)   :: self 
   !! input quaternion
  real(kind=sgl), INTENT(IN)       :: s
   !! scalar input 
  type(Quaternion_T)               :: qres 
   !! output quaternion

  qres%q = (/ s*self%q(1), s*self%q(2), s*self%q(3), s*self%q(4) /) 
  qres%s = 's'

end function quatsmult

!--------------------------------------------------------------------------
pure recursive function quatarraysmult(self, s) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! scalar quaternion array multiplication   (single precision)

IMPLICIT NONE

  class(QuaternionArray_T),intent(in)   :: self 
   !! input quaternion
  real(kind=sgl), INTENT(IN)            :: s
   !! scalar input 
  type(QuaternionArray_T)               :: qres 
   !! output quaternion

  integer(kind=irg)                     :: sz(2) 

! if the quaternion array is already allocated, check to make sure it has the right dimensions
  if (allocated(qres%q)) then 
    sz = shape(qres%q)
    if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%q)
  end if 
  allocate(qres%q(4,self%n))

  qres%q = s*self%q
  qres%s = self%s
  qres%n = self%n
  qres%nthreads = self%nthreads

end function quatarraysmult

!--------------------------------------------------------------------------
pure recursive function quatsmultd(self, s) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! scalar quaternion multiplication   (double precision)

IMPLICIT NONE

  class(Quaternion_T),intent(in)   :: self 
   !! input quaternion
  real(kind=dbl), INTENT(IN)       :: s
   !! scalar input 
  type(Quaternion_T)               :: qres 
   !! output quaternion

  qres%qd = (/ s*self%qd(1), s*self%qd(2), s*self%qd(3), s*self%qd(4) /) 
  qres%s = 'd'

end function quatsmultd

!--------------------------------------------------------------------------
pure recursive function quatarraysmultd(self, s) result(qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! scalar quaternion array multiplication   (double precision)

IMPLICIT NONE

  class(QuaternionArray_T),intent(in)   :: self 
   !! input quaternion
  real(kind=dbl), INTENT(IN)            :: s
   !! scalar input 
  type(QuaternionArray_T)               :: qres 
   !! output quaternion

  integer(kind=irg)                     :: sz(2) 

! if the quaternion array is already allocated, check to make sure it has the right dimensions
  if (allocated(qres%qd)) then 
    sz = shape(qres%qd)
    if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%qd)
  end if 
  allocate(qres%qd(4,self%n))

  qres%qd = s*self%qd
  qres%s = self%s
  qres%n = self%n
  qres%nthreads = self%nthreads

end function quatarraysmultd

!--------------------------------------------------------------------------
pure recursive function quatconjg(self) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion conjugation (extends intrinsic routine conjg)

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion
  type(Quaternion_T)                :: qres
   !! output quaternion

  if (self%s.eq.'s') then 
    qres%q = (/ self%q(1), -self%q(2), -self%q(3), -self%q(4) /) 
    qres%s = 's'
  else
    qres%qd = (/ self%qd(1), -self%qd(2), -self%qd(3), -self%qd(4) /) 
    qres%s = 'd'
  end if 

end function quatconjg

!--------------------------------------------------------------------------
pure recursive function quatarrayconjg(self) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! quaternion array conjugation (extends intrinsic routine conjg)

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)    :: self
   !! input quaternion
  type(QuaternionArray_T)                :: qres
   !! output quaternion

  integer(kind=irg)                      :: sz(2) 

  if (self%s.eq.'s') then 
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%q)) then 
      sz = shape(qres%q)
      if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%q)
    end if 
    allocate(qres%q(4,self%n))

    qres%q(1,:)   =  self%q(1,:)
    qres%q(2:4,:) = -self%q(2:4,:)
  else
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%qd)) then 
      sz = shape(qres%qd)
      if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%qd)
    end if 
    allocate(qres%qd(4,self%n))

    qres%qd(1,:)   =  self%qd(1,:)
    qres%qd(2:4,:) = -self%qd(2:4,:)
  end if 

  qres%s = self%s
  qres%n = self%n
  qres%nthreads = self%nthreads

end function quatarrayconjg

!--------------------------------------------------------------------------
pure recursive function quatnorm(self) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion norm (extends intrinsic routine abs)

IMPLICIT NONE 

  class(Quaternion_T),intent(in) :: self
   !! input quaternion
  real(kind=dbl)                 :: res
   !! output norm

  real(kind=sgl)                 :: n
  real(kind=dbl)                 :: nd, resd

  if (self%s.eq.'s') then
    n = self%q(1)**2 + self%q(2)**2 + self%q(3)**2 + self%q(4)**2
    resd = dsqrt( dble(n) )
    res = dble(sngl(resd))
  else
    nd = self%qd(1)**2 + self%qd(2)**2 + self%qd(3)**2 + self%qd(4)**2
    res = dsqrt( nd )
  end if 

end function quatnorm

!--------------------------------------------------------------------------
pure recursive function quatarraynorm(self) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! quaternion array norm (extends intrinsic routine abs)
  !!
  !! this routine requires the output array to be allocated in the calling program.

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in) :: self
   !! input quaternion
  real(kind=dbl)                      :: res(self%n)
   !! output norm

  real(kind=sgl),allocatable          :: n(:)
  real(kind=dbl),allocatable          :: nd(:), resd(:)

  if (self%s.eq.'s') then
    allocate(n(self%n), resd(self%n))
    n(:) = self%q(1,:)**2 + self%q(2,:)**2 + self%q(3,:)**2 + self%q(4,:)**2
    resd = dsqrt( dble(n) )
    res = dble(sngl(resd))
    deallocate(n, resd)
  else
    allocate(nd(self%n)) 
    nd(:) = self%qd(1,:)**2 + self%qd(2,:)**2 + self%qd(3,:)**2 + self%qd(4,:)**2
    res = dsqrt( nd )
    deallocate(nd)
  end if 

end function quatarraynorm

!--------------------------------------------------------------------------
recursive subroutine quatnormalize(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/07/20
  !!
  !! normalize the input quaternion 

IMPLICIT NONE 

  class(Quaternion_T),intent(inout) :: self
   !! input quaternion

  type(Quaternion_T)                :: q
  real(kind=sgl)                    :: n
  real(kind=dbl)                    :: nd

  if (self%s.eq.'s') then
    n = self%q(1)**2 + self%q(2)**2 + self%q(3)**2 + self%q(4)**2
    n = sqrt( n )
    q = self%quatsdiv(n)
    self%q = q%q
  else
    nd = self%qd(1)**2 + self%qd(2)**2 + self%qd(3)**2 + self%qd(4)**2
    nd = sqrt( nd )
    q = self%quatsdivd(nd)
    self%qd = q%qd
  end if 

end subroutine quatnormalize

!--------------------------------------------------------------------------
recursive subroutine quatarraynormalize(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/07/20
  !!
  !! normalize the input quaternion 

IMPLICIT NONE 

  class(QuaternionArray_T),intent(inout) :: self
   !! input quaternion

  real(kind=sgl),allocatable             :: n(:)
  real(kind=dbl),allocatable             :: nd(:)
  integer(kind=irg)                      :: i

  if (self%s.eq.'s') then
    allocate(n(self%n))
    n(:) = self%q(1,:)**2 + self%q(2,:)**2 + self%q(3,:)**2 + self%q(4,:)**2
    n = sqrt( n )
    do i=1,self%n
      self%q(:,i) = self%q(:,i)/n(i)
    end do
    deallocate(n)
  else
    allocate(nd(self%n))
    nd(:) = self%qd(1,:)**2 + self%qd(2,:)**2 + self%qd(3,:)**2 + self%qd(4,:)**2
    nd = sqrt( nd )
    do i=1,self%n
      self%qd(:,i) = self%qd(:,i)/nd(i)
    end do
    deallocate(nd)
  end if 

end subroutine quatarraynormalize

!--------------------------------------------------------------------------
recursive function quatdiv(self, y) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion division (single/double precision)

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self, y
   !! input quaternions
  type(Quaternion_T)                :: qres
   !! output quaternion 

  type(Quaternion_T)                :: p, cy
  real(kind=sgl)                    :: q
  real(kind=dbl)                    :: qd

  if (self%s.eq.'s') then 
      q = quatnorm(y)
      cy = quatconjg(y)
      p = quatsdiv( cy, q*q )
      qres = quatmult(self,p)
      qres%s = 's'
  else
      qd = quatnorm(y)
      cy = quatconjg(y)
      p = quatsdivd( cy, qd*qd )
      qres = quatmult(self,p)
      qres%s = 'd'
  end if 


end function quatdiv

!--------------------------------------------------------------------------
recursive function quatarraydiv(self, y) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion array division (single/double precision)

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)    :: self, y
   !! input quaternions
  type(QuaternionArray_T)                :: qres
   !! output quaternion 

  type(QuaternionArray_T)                :: p, cy
  real(kind=sgl),allocatable             :: q(:)
  real(kind=dbl),allocatable             :: qd(:)
  integer(kind=irg)                      :: i

  if (self%s.eq.'s') then 
      q = quatarraynorm(y)
      cy = quatarrayconjg(y)
      do i=1,self%n 
        cy%q(:,i) = cy%q(:,i) / q(i)**2
      end do
      qres = quatarraymult(self,cy)
  else
      qd = quatarraynorm(y)
      cy = quatarrayconjg(y)
      do i=1,self%n 
        cy%qd(:,i) = cy%qd(:,i) / qd(i)**2
      end do
      qres = quatarraymult(self,cy)
  end if 

end function quatarraydiv

!--------------------------------------------------------------------------
recursive function quatsdiv(self, s) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion division (single precision)

use mod_io 

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion (numerator)
  real(kind=sgl), INTENT(IN)        :: s            
   !! input quaternion (denominator)

  type(Quaternion_T)                :: qres
  type(IO_T)                        :: Message

  if (s.ne.0.0) then 
      qres%q = (/ self%q(1)/s, self%q(2)/s, self%q(3)/s, self%q(4)/s /) 
      qres%s = 's'
  else 
    call Message % printWarning('quatsdiv', (/ 'Attempting to divide quaternion by zero; skipping operation ...' /) )
    qres = self
  end if

end function quatsdiv

!--------------------------------------------------------------------------
recursive function quatarraysdiv(self, s) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! quaternion array scalar division (single precision)

use mod_io 

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)    :: self
   !! input quaternion (numerator)
  real(kind=sgl), INTENT(IN)             :: s            
   !! input quaternion (denominator)

  type(QuaternionArray_T)                :: qres
  type(IO_T)                             :: Message
  integer(kind=irg)                      :: sz(2)

  if (s.ne.0.0) then 
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%q)) then 
      sz = shape(qres%q)
      if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%q)
    end if 
    allocate(qres%q(4,self%n))

    qres%q = self%q/s
    qres%s = self%s
    qres%n = self%n
    qres%nthreads = self%nthreads
  else 
    call Message % printError('quatarraysdiv', 'Attempting to divide quaternion aray by zero' )
  end if

end function quatarraysdiv

!--------------------------------------------------------------------------
recursive function quatsdivd(self, s) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion division (doubgle precision)

use mod_io 

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion (numerator)
  real(kind=dbl), INTENT(IN)        :: s            
   !! input quaternion (denominator)

  type(Quaternion_T)                :: qres
  type(IO_T)                        :: Message

  if (s.ne.0.0) then 
    qres%qd = (/ self%qd(1)/s, self%qd(2)/s, self%qd(3)/s, self%qd(4)/s /) 
    qres%s = 'd'
  else 
    call Message % printWarning('quatsdivd', (/ 'Attempting to divide quaternion by zero; skipping operation ...' /) )
    qres = self
  end if

end function quatsdivd

!--------------------------------------------------------------------------
recursive function quatarraysdivd(self, s) result (qres)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! quaternion array scalar division (double precision)

use mod_io 

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)    :: self
   !! input quaternion (numerator)
  real(kind=dbl), INTENT(IN)             :: s            
   !! input quaternion (denominator)

  type(QuaternionArray_T)                :: qres
  type(IO_T)                             :: Message
  integer(kind=irg)                      :: sz(2)

  if (s.ne.0.D0) then 
! if the quaternion array is already allocated, check to make sure it has the right dimensions
    if (allocated(qres%qd)) then 
      sz = shape(qres%qd)
      if ((sz(1).ne.4).or.(sz(2).ne.self%n)) deallocate(qres%qd)
    end if 
    allocate(qres%qd(4,self%n))

    qres%qd = self%qd/s
    qres%s = self%s
    qres%n = self%n
    qres%nthreads = self%nthreads
  else 
    call Message % printError('quatarraysdivd', 'Attempting to divide quaternion aray by zero' )
  end if

end function quatarraysdivd

!--------------------------------------------------------------------------
pure recursive function quatinnerproduct(self, y) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion inner product (single precision)

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self, y
   !! input quaternions
  real(kind=dbl)                    :: res
   !! inner product 

  if (self%s.eq.'s') then 
    res = dble(self%q(1) * y%q(1) + self%q(2) * y%q(2) + self%q(3) * y%q(3) + self%q(4) * y%q(4))
  else
    res = self%qd(1) * y%qd(1) + self%qd(2) * y%qd(2) + self%qd(3) * y%qd(3) + self%qd(4) * y%qd(4)
  end if

end function quatinnerproduct

!--------------------------------------------------------------------------
recursive function quatarrayinnerproduct(self, y) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! quaternion array inner product (single precision)
  !!
  !! calling program must allocate the output array

use mod_io 

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)    :: self, y
   !! input quaternions
  real(kind=dbl)                         :: res(self%n)
   !! inner product 

   type(IO_T)                            :: Message

! test to make sure that both arrays have the same number of quaternions 
  if (self%n.ne.y%n) then 
    call Message%printError('quatarrayinnerproduct','input arrays must have the same number of quaternions')
  end if 

  if (self%s.ne.y%s) then 
    call Message%printError('quatarrayinnerproduct','input arrays must have the same precision')
  end if 
  
  if (self%s.eq.'s') then 
    res(:) = dble(self%q(1,:) * y%q(1,:) + self%q(2,:) * y%q(2,:) + self%q(3,:) * y%q(3,:) + self%q(4,:) * y%q(4,:))
  else
    res(:) = self%qd(1,:) * y%qd(1,:) + self%qd(2,:) * y%qd(2,:) + self%qd(3,:) * y%qd(3,:) + self%qd(4,:) * y%qd(4,:)
  end if

end function quatarrayinnerproduct

!--------------------------------------------------------------------------!
pure recursive function quatangle(self, y) result(res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! interquaternion angle   (single/double precision)
  !! this only has meaning for a unit quaternion, so we test first and return -10000.0
  !! if either of the quaternions is not a unit quaternion.

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self, y
   !! input quaternions
  real(kind=dbl)                    :: res
   !! angle (radians)

  real(kind=sgl)                    :: q, nself, ny
  real(kind=dbl)                    :: qd, nselfd, nyd

  if (self%s.eq.'s') then 
      nself = self%quatnorm()
      ny = y%quatnorm()
      q = self%quat_innerproduct(y)
      res = dble(acos( q/(nself * ny) ))
  else 
      nselfd = self%quatnorm()
      nyd = y%quatnorm()
      qd = self%quat_innerproduct(y)
      res = dacos( qd/(nselfd * nyd) )
  end if

end function quatangle

!--------------------------------------------------------------------------!
recursive function quatarrayangle(self, y) result(res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! interquaternion angle for a pair of quaternion arrays (single/double precision)
  !! this only has meaning for a unit quaternion, so we test first and return -10000.0
  !! if either of the quaternions is not a unit quaternion.

use mod_io

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)     :: self, y
   !! input quaternions
  real(kind=dbl)                          :: res(self%n)
   !! angle (radians)

  type(IO_T)                              :: Message
  real(kind=sgl), allocatable             :: q(:), nself(:), ny(:)
  real(kind=dbl), allocatable             :: qd(:), nselfd(:), nyd(:)

! test to make sure that both arrays have the same number of quaternions 
  if (self%n.ne.y%n) then 
    call Message%printError('quatarrayangle','input arrays must have the same number of quaternions')
  end if 

  if (self%s.ne.y%s) then 
    call Message%printError('quatarrayangle','input arrays must have the same precision')
  end if 

  if (self%s.eq.'s') then 
      allocate(nself(self%n), ny(self%n), q(self%n)) 
      nself = self%quatarraynorm()
      ny = y%quatarraynorm()
      q = self%quatarrayinnerproduct(y)
      res = dble(acos( q/(nself * ny) ))
      deallocate(nself, ny, q)
  else 
      allocate(nselfd(self%n), nyd(self%n), qd(self%n)) 
      nselfd = self%quatarraynorm()
      nyd = y%quatarraynorm()
      qd = self%quatarrayinnerproduct(y)
      res = dacos( qd/(nselfd * nyd) )
      deallocate(nselfd, nyd, qd)
  end if

end function quatarrayangle


!--------------------------------------------------------------------------!
! pure recursive function quatLp(self, v) result (res)
recursive function quatLp(self, v) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! actively rotate a unit vector by a unit quaternion, L_p = p v p* (single precision)

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion
  real(kind=sgl),intent(in)         :: v(3)         
   !! input vector to be rotated 
  real(kind=sgl)                    :: res(3)
   !! output vector

  type(Quaternion_T)                :: qv, rqv, cq

  qv%q = (/ 0.0, v(1), v(2), v(3) /) 
  qv%s = 's'
  cq = quatconjg(self)
  rqv = quatmult(self, quatmult(qv, cq) )
  res(1:3) = rqv%q(2:4)

end function quatLp

!--------------------------------------------------------------------------!
recursive function quatarrayLp(self, v) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! actively rotate a single vector by an array of unit quaternions, L_p = p v p* (single precision)

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)    :: self
   !! input quaternion
  real(kind=sgl),intent(in)              :: v(3)         
   !! input vector to be rotated 
  real(kind=sgl)                         :: res(3,self%n)
   !! output vector

  type(Quaternion_T)                     :: qv, rqv, q, cq 
  integer(kind=irg)                      :: i

  qv%q = (/ 0.0, v(1), v(2), v(3) /) 
  qv%s = 's'
  do i=1,self%n
    q = extractfromQuaternionArray(self, i)
    cq = q%quatconjg()
    rqv = quatmult(q, quatmult(qv, cq) )
    res(1:3,i) = rqv%q(2:4)
  end do 

end function quatarrayLp

!--------------------------------------------------------------------------!
! pure recursive function quatLpd(self, v) result (res)
recursive function quatLpd(self, v) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/03/20
  !!
  !! actively rotate a unit vector by a unit quaternion, L_p = p v p* (double precision)

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self
   !! input quaternion
  real(kind=dbl),intent(in)         :: v(3)         
   !! input vector to be rotated 
  real(kind=dbl)                    :: res(3)
   !! output vector

  type(Quaternion_T)                :: qv, rqv, cq

  qv%qd = (/ 0.D0, v(1), v(2), v(3) /) 
  qv%s = 'd'
  cq = quatconjg(self)
  rqv = quatmult(self, quatmult(qv, cq) )
  res(1:3) = rqv%qd(2:4)

end function quatLpd

!--------------------------------------------------------------------------!
recursive function quatArrayLpd(self, v) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/03/20
  !!
  !! actively rotate a unit vector by an array of unit quaternions, L_p = p v p* (double precision)

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)    :: self
   !! input quaternion
  real(kind=dbl),intent(in)              :: v(3)         
   !! input vector to be rotated 
  real(kind=dbl)                         :: res(3,self%n)
   !! output vector

  type(Quaternion_T)                     :: qv, rqv, q, cq 
  integer(kind=irg)                      :: i

  qv%qd = (/ 0.D0, v(1), v(2), v(3) /) 
  qv%s = 'd'
  do i=1,self%n
    q = extractfromQuaternionArray(self, i)
    cq = q%quatconjg()
    rqv = quatmult(q, quatmult(qv, cq) )
    res(1:3,i) = rqv%qd(2:4)
  end do 

end function quatArrayLpd

!--------------------------------------------------------------------------!
recursive function extractfromQuaternionArray(self, i) result (res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! extract a quaternion from an array of quaternions

use mod_io 

IMPLICIT NONE 

  class(QuaternionArray_T),intent(in)   :: self
   !! input quaternion array
  integer(kind=irg), intent(in)         :: i 
   !! quaternion to be extracted 
  type(Quaternion_T)                    :: res
   !! extracted quaternion 

  type(IO_T)                            :: Message 

  if (i.le.self%n) then 
    res%s = self%s 
    if (self%s.eq.'s') then 
      res%q(:) = self%q(:,i) 
    else
      res%qd(:) = self%qd(:,i) 
    end if 
  else 
    call Message%printWarning('extractfromQuaternionArray: requested quaternion index larger than array size', &
                              (/'   ---> returning empty quaternion'/) )
    if (self%s.eq.'s') then 
      res = Quaternion_T( q = (/ 0.0, 0.0, 0.0, 0.0 /) )
    else
      res = Quaternion_T( qd = (/ 0.D0, 0.D0, 0.D0, 0.D0 /) )
    end if
  end if 

end function extractfromQuaternionArray


!--------------------------------------------------------------------------!
recursive subroutine insertQuatintoArray(self, i, q)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/23/20
  !!
  !! insert a quaternion into an array of quaternions 

use mod_io 

IMPLICIT NONE 

  class(QuaternionArray_T),intent(inout):: self
   !! input quaternion array
  integer(kind=irg), intent(in)         :: i 
   !! quaternion to be extracted 
  type(Quaternion_T), intent(in)        :: q
   !! extracted quaternion 

  type(IO_T)                            :: Message 

  if (i.le.self%n) then 
    if (self%s.eq.'s') then 
      self%q(:,i) = q%get_quats()
    else
      self%qd(:,i) = q%get_quatd()
    end if 
  else 
    call Message%printWarning('extractfromQuaternionArray: requested quaternion index larger than array size', &
                              (/'   ---> no quaternion inserted'/) )
  end if 

end subroutine insertQuatintoArray

!--------------------------------------------------------------------------!
! pure recursive function quatslerp(self, qb, n) result(res)
recursive function quatslerp(self, qb, n) result(res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! return an array of interpolated quaternions

IMPLICIT NONE 

  class(Quaternion_T),intent(in)         :: self   ! = qa 
   !! input quaternion (start)
  class(Quaternion_T),intent(in)         :: qb
   !! input quaternion (end)
  integer(kind=irg),intent(in)           :: n            
   !! number of steps in the interpolation
  type(Quaternion_T)                     :: res(n)
   !! output interpolated quaternion list

  type(Quaternion_T)                     :: cqa
  real(kind=sgl)                         :: theta, phi, dphi, s
  real(kind=dbl)                         :: thetad, phid, dphid, sd
  integer(kind=irg)                      :: i

  if (self%s.eq.'s') then 
      do i=1,n 
        res(i)%q = (/ 0.0, 0.0, 0.0, 0.0 /) 
        res(i)%s = 's'
      end do
      cqa = quatconjg(self)
      theta = acos( quatinnerproduct(qb, cqa ))
      if (theta.ne.0.0) then
        s = 1.0/sin(theta)
        dphi = theta/real(n-1)

        do i=1,n
          phi = real(i-1)*dphi
          res(i) = quatsmult(self, sin(theta-phi)*s ) + quatsmult(qb, sin(phi)*s )
        end do
      else
        do i=1,n
          res(i) = self
        end do
      end if
  else
      do i=1,n 
        res(i)%qd = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
        res(i)%s = 'd'
      end do
      cqa = quatconjg(self)
      thetad = acos( quatinnerproduct(qb, cqa ))
      if (thetad.ne.0.D0) then
        sd = 1.D0/sin(theta)
        dphi = theta/dble(n-1)

        do i=1,n
          phi = dble(i-1)*dphi
          res(i) = quatsmultd(self, sin(theta-phi)*sd ) + quatsmultd(qb, sin(phi)*sd )
        end do
      else
        do i=1,n
          res(i) = self
        end do
      end if
  end if 

end function quatslerp

!--------------------------------------------------------------------------
recursive function quatsequal(self, qb) result(res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/06/20
  !!
  !! quaternion comparison (double precision)

IMPLICIT NONE 

  class(Quaternion_T),intent(in)    :: self, qb
   !! input quaternions 
  logical                           :: res

  type(Quaternion_T)                :: diff

  real(kind=sgl)                    :: d, eps=1.0e-7
  real(kind=dbl)                    :: dd, epsd=1.0e-12

  res = .TRUE.
  diff = self - qb 

  if (self%s.eq.'s') then 
    d = maxval( abs( diff%q(:) ) )
    if (d.gt.eps) res = .FALSE.
  else 
    dd = maxval( abs( diff%qd(:) ) )
    if (dd.gt.epsd) res = .FALSE.
  end if

end function quatsequal


!--------------------------------------------------------------------------
recursive function generateRandomArray(n, s, seed, northern) result(res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! generate an array of random unit quaternions (single/double precision)

use mod_rng

  integer(kind=irg), intent(in)       :: n
   !! number of unut quaternions to be generated 
  character(1), intent(in)            :: s 
   !! single 's' or double 'd' precision ?
  type(rng_t), intent(inout)          :: seed 
   !! a seed number for the Marsaglia random number generator 
  logical, intent(in), OPTIONAL       :: northern 

  type(QuaternionArray_T)             :: res 
  type(Quaternion_T)                  :: q 
  integer(kind=irg)                   :: i
  logical                             :: pos 

  pos = .FALSE.
  if (present(northern)) then 
    if (northern.eqv..TRUE.) pos = .TRUE.
  end if

  res%s = s 
  if (s.eq.'s') then 
    allocate(res%q(4,n))
  else
    allocate(res%qd(4,n))
  end if
  res%n = n

  do i=1,n 
    q = quat_Marsaglia(seed)
    if (pos.eqv..TRUE.) then
      if (q%qd(1).lt.0.0) q%qd = -q%qd 
    end if
    if (s.eq.'s') then 
      res%q(:,i) = sngl(q%qd(:))
    else 
      res%qd(:,i) = q%qd(:)
    end if
  end do 

end function generateRandomArray

!--------------------------------------------------------------------------!
recursive function quat_Marsaglia(seed) result(q)
!DEC$ ATTRIBUTES DLLEXPORT :: quatMarsaglia

use mod_rng 

IMPLICIT NONE

type(rng_t),INTENT(INOUT)           :: seed 
!f2py intent(in,out) ::  seed 
type(Quaternion_T)                  :: q 

real(kind=dbl)                      :: x1,x2,y1,y2,s1,s2

x1 = 2.D0*rng_uniform(seed)-1.D0
y1 = 2.D0*rng_uniform(seed)-1.D0
s1 = x1*x1+y1*y1
if (s1.gt.1.D0) then 
  do while (s1.gt.1.D0) 
    x1 = 2.D0*rng_uniform(seed)-1.D0
    x2 = 2.D0*rng_uniform(seed)-1.D0
    s1 = x1*x1+y1*y1
  end do
end if 

x2 = 2.D0*rng_uniform(seed)-1.D0
y2 = 2.D0*rng_uniform(seed)-1.D0
s2 = x2*x2+y2*y2
if (s2.gt.1.D0) then 
  do while (s2.gt.1.D0) 
    x2 = 2.D0*rng_uniform(seed)-1.D0
    y2 = 2.D0*rng_uniform(seed)-1.D0
    s2 = x2*x2+y2*y2
  end do
end if 

s1 = sqrt( (1.D0-s2)/s2 )

q%qd = (/ x1, y1, x2*s1, y2*s1 /) 
q%s = 'd'

end function quat_Marsaglia

!--------------------------------------------------------------------------
recursive subroutine QSym_Init_(self, pgnum, qsym) 
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! initialize the quaternion symmetry operators for a given rotational point group

use mod_symmetry
use mod_io

IMPLICIT NONE

class(QuaternionArray_T), INTENT(INOUT)   :: self
integer(kind=irg), INTENT(IN)             :: pgnum
type(QuaternionArray_T), INTENT(OUT)      :: qsym

type(IO_T)                                :: Message
integer(kind=irg)                         :: i, Nqsym, prot
real(kind=dbl), allocatable               :: Pm(:,:)

! here we need to analyze the rotational symmetry group, and copy the appropriate 
! quaternion symmetry operators into a QuaternionArray_T class.

! first get the number of the rotational point group that corresponds to the crystal point group
prot = PGrot(pgnum)
! possible values for prot are: (/1,3,6,9,12,16,18,21,24,28,30/)
! corresponding to the point groups 1, 2, 222, 4, 422, 3, 32, 6, 622, 23, 432 and 532 respectively

!------------
! IMPORTANT NOTE: the original von Mises-Fischer (VMF) approach requires that q and -q are considered to 
! be separate quaternions, so the original Matlab code included the negatives of all quaternion symmetry operators 
! as well, leading to a cardinality of twice the rotational point group order.  It appears that we do not have to
! do so if we replace the exponential in the VMF by a hyperbolic cosine function, which would account directly
! for the q, -q duplicity... Alternatively, one can use the axial Watson distribution.
!------------

! identity operator is part of all point groups

! select statement for each individual rotational point group (see typedefs.f90 for SYM_Qsymop definitions)
select case (prot) 
        case(1)         ! 1 (no additional symmetry elements)
                allocate(Pm(4,1))
                Pm(1:4,1) = SYM_Qsymop(1:4,1)
                Nqsym = 1

        case(3)         ! 2  (we'll assume that the two-fold axis lies along the e_y-axis)
                allocate(Pm(4,2))
                Pm(1:4,1) = SYM_Qsymop(1:4,1)
                Nqsym = 2
                Pm(1:4,2) = SYM_Qsymop(1:4,3)

        case(6)         ! 222
                allocate(Pm(4,4))
                Pm(1:4,1) = SYM_Qsymop(1:4,1)
                Nqsym = 4
                do i=2,4
                  Pm(1:4,i) = SYM_Qsymop(1:4,i)
                end do

        case(9)         ! 4
                allocate(Pm(4,4))
                Pm(1:4,1) = SYM_Qsymop(1:4,1)
                Nqsym = 4
                Pm(1:4,2) = SYM_Qsymop(1:4,4)
                Pm(1:4,3) = SYM_Qsymop(1:4,7)
                Pm(1:4,4) = SYM_Qsymop(1:4,10)

        case(12)        ! 422
                allocate(Pm(4,8))
                Pm(1:4,1) = SYM_Qsymop(1:4,1)
                Nqsym = 8
                Pm(1:4,2) = SYM_Qsymop(1:4,4)
                Pm(1:4,3) = SYM_Qsymop(1:4,7)
                Pm(1:4,4) = SYM_Qsymop(1:4,10)
                Pm(1:4,5) = SYM_Qsymop(1:4,2)
                Pm(1:4,6) = SYM_Qsymop(1:4,3)
                Pm(1:4,7) = SYM_Qsymop(1:4,11)
                Pm(1:4,8) = SYM_Qsymop(1:4,12)

        case(16)        ! 3
                allocate(Pm(4,3))
                Pm(1:4,1) = SYM_Qsymop(1:4,1)
                Nqsym = 3
                Pm(1:4,2) = SYM_Qsymop(1:4,26)
                Pm(1:4,3) = SYM_Qsymop(1:4,28)

        case(18)        ! 32 (needs special handling)
                allocate(Pm(4,6))
                Pm(1:4,1) = SYM_Qsymop(1:4,1)
                Nqsym = 6
                Pm(1:4,2) = SYM_Qsymop(1:4,26)
                Pm(1:4,3) = SYM_Qsymop(1:4,28)
                Pm(1:4,4) = SYM_Qsymop(1:4,30)
                Pm(1:4,5) = SYM_Qsymop(1:4,32)
                Pm(1:4,6) = SYM_Qsymop(1:4,34)

        case(21)        ! 6
                allocate(Pm(4,6))
                Pm(1:4,1) = SYM_Qsymop(1:4,1)
                Nqsym = 6
                do i=25,29
                  Pm(1:4,i-23) = SYM_Qsymop(1:4,i)
                end do

        case(24)        ! 622
                allocate(Pm(4,6))
                Pm(1:4,1) = SYM_Qsymop(1:4,1)
                Nqsym = 12
                do i=25,35
                  Pm(1:4,i-23) = SYM_Qsymop(1:4,i)
                end do

        case(28)        ! 23
                allocate(Pm(4,12))
                Pm(1:4,1) = SYM_Qsymop(1:4,1)
                Nqsym = 12
                do i=2,4
                  Pm(1:4,i) = SYM_Qsymop(1:4,i)
                end do
                do i=17,24
                  Pm(1:4,4+(i-16)) = SYM_Qsymop(1:4,i)
                end do

        case(30)        ! 432
                allocate(Pm(4,24))
                Pm(1:4,1) = SYM_Qsymop(1:4,1)
                Nqsym = 24
                do i=2,24
                  Pm(1:4,i) = SYM_Qsymop(1:4,i)
                end do
        case(33)        ! 532
                allocate(Pm(4,60))
                Pm(1:4,1) = SYM_Qsymop(1:4,1)
                Nqsym = 60
                do i=2,60
                  Pm(1:4,i) = SYM_Qsymop(1:4,35+i)
                end do
        case(34)  ! 822
                allocate(Pm(4,16))
                Pm(1:4,1) = SYM_Qsymop(1:4,1)
                Nqsym = 16
                do i = 1,15
                  Pm(1:4,i+1) = SYM_Qsymop(1:4,95+i)
                end do
        case(35)  ! 1022
                allocate(Pm(4,20))
                Pm(1:4,1) = SYM_Qsymop(1:4,1)
                Nqsym = 20
                do i = 1,19
                  Pm(1:4,i+1) = SYM_Qsymop(1:4,110+i)
                end do
        case(36)  ! 1222
                allocate(Pm(4,24))
                Pm(1:4,1) = SYM_Qsymop(1:4,1)
                Nqsym = 24
                do i = 1,23
                  Pm(1:4,i+1) = SYM_Qsymop(1:4,129+i)
                end do

        case default    ! this should never happen ...
                write (*,*) 'requested rotational point group ', prot
                call Message%printError('QSym_Init','unknown rotational point group number')
end select

! and initialize the quatenrion array class 
qsym = QuaternionArray_T( n = Nqsym, nthreads = 1, qd = Pm )

end subroutine QSym_Init_

!--------------------------------------------------------------------------
recursive function getQnumber_(self) result(num) 
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! returns the number of quaternions in the QuaternionArray_T class

IMPLICIT NONE

class(QuaternionArray_T), INTENT(INOUT)   :: self
integer(kind=irg)                         :: num

num = self%n

end function getQnumber_

end module mod_quaternions
