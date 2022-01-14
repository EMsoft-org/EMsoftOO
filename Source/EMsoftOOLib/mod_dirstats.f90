! ###################################################################
! Copyright (c) 2014-2022, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_dirstats
 !! author: MDG
 !! version: 1.0
 !! date: 01/23/20
 !!
 !! Directional Statistics routines
 !!
 !! This module contains all routines that deal with directional statistics
 !! using a symmetrized quaternion unit sphere with the modified von Mises-Fisher
 !! distribution, or the modified mixture of axial Watson distributions.  All
 !! the details for this approach can be found in two papers:
 !!
 !! Y.H. Chen, Park S.U., D. Wei, G. Newstadt, M. Jackson, J.P. Simmons, M. De Graef,
 !! and A.O. Hero. “A selfionary approach to the EBSD indexing problem”.
 !! Microsc. MicroAnal. 21, 739-752 (2015).
 !!
 !! Y.H. Chen, D. Wei, G. Newstadt, M. De Graef, J.P. Simmons, and A.O. Hero.
 !! “Parameter estimation in spherical symmetry groups”.
 !! IEEE. Sign. Proc. Lett. 22, 1152-1155 (2015).
 !!
 !! Here is an example program showing how the routines can be called:
 !!
 !! program t
 !!
 !! use mod_kinds
 !! use mod_global
 !! use mod_quaternions
 !! use mod_dirstats
 !!
 !! integer(kind=irg)               :: nums, seed
 !! real(kind=dbl),allocatable      :: samples(:,:)
 !! type(selftype)                  :: self
 !! real(kind=dbl)                  :: muhat(4), kappahat
 !!
 !! ! this is a test of the selfionary indexing portion that deals with the
 !! ! modified von Mises-Fisher distribution; the results must be the same
 !! ! as those produced by the original Matlab code...
 !!
 !! seed = 432514
 !! nums = 1000
 !! allocate(samples(4,nums))
 !!
 !! allocate(self)
 !! self%NumEM = 3
 !! self%NumIter = 30
 !! self%pgnum = 32
 !!
 !! ! read a bunch of quaternions from a file... and store them in samples
 !!
 !! call DI_Init(self,'VMF') ! replace 'VMF' by 'WAT' to use the axial Watson distribution
 !!
 !! do i=1,10
 !!   call DI_EMforDS(samples, self, nums, seed, muhat, kappahat,'VMF')  ! replace 'VMF' by 'WAT' to use the Watson distribution
 !!
 !!   write (*,*) '  '
 !!   write (*,*) 'mu    = ',muhat
 !!   write (*,*) 'kappa = ',kappahat
 !!   write (*,*) 'equivalent angular precision : ',180.D0*dacos(1.D0-1.D0/kappahat)/cPi
 !! end do
 !!
 !! end program
 !!
 !! original history of the dictmod.f90 module (EMsoft v. 5.X and earlier)
 !! @date 12/31/14 MDG 1.0 original (based on UMich Matlab code and IDL intermediate version)
 !! @date 01/02/15 MDG 1.1 debug of code; produces same result as Matlab code
 !! @date 01/04/15 MDG 1.2 trial implementation of model using hyperbolic functions instead of exponential
 !! @date 01/06/15 MDG 1.3 changed public routine names with DI_ in front
 !! @date 01/07/15 MDG 1.4 added VMF sampling routines
 !! @date 01/09/15 MDG 1.5 replaced several computations by numerically more stable versions
 !! @date 02/05/15 MDG 1.6 added sampling for axial Watson distribution
 !! @date 02/06/15 MDG 1.7 streamlined sampling code by removing duplications; VMF vs. Watson is now an argument
 !! @date 02/06/15 MDG 1.8 general rewrite and removal of duplications in indexing routines; some name changes
 !!

use mod_kinds
use mod_global
use mod_quaternions

IMPLICIT NONE

! these need to be moved elsewhere, for now they are commented out ...
 ! Similarity_Classifier
 ! CardIntersection

!--------------------------------------------------------------------------
! this is a portion of the old dicttype structure, but now turned into a class
! DirStat = Directional Statistics
!--------------------------------------------------------------------------
type, public :: DirStat_T
private
  real(kind=dbl),allocatable          :: xAp(:)    ! kappa array
  real(kind=dbl),allocatable          :: yAp(:)    ! A_4(u) lookup table
  integer(kind=irg)                   :: pgnum     ! point group number; needed to define symmetry operators
  integer(kind=irg)                   :: Apnum     ! number of entries in lookup table
  integer(kind=irg)                   :: NumEM     ! number of times that the EM algorithm needs to be carried out (set by user)
  integer(kind=irg)                   :: NumIter   ! number of iterations inside each EM call (set by user)
  character(3)                        :: DStype    ! 'VMF' or 'WAT', sets the distribution type
  type(QuaternionArray_T)             :: Xquats    ! quaternion array (needs routines to set and get)
  type(QuaternionArray_T)             :: Xaux      ! auxiliary quaternion array
  type(QuaternionArray_T)             :: qsym      ! quaternion symmetry operators
  integer(kind=irg)                   :: N         ! number of samples in Xquats (to keep things simple)
  integer(kind=irg)                   :: Pmdims    ! number of symmetry operators
  type(Quaternion_T)                  :: Mumean    ! mean quaternion direction
  real(kind=dbl)                      :: kappa     ! concentration parameter

contains
private

  procedure, pass(self) :: RotateToMu_
  procedure, pass(self) :: randUniformSphere_
  procedure, pass(self) :: SampleDS_
  procedure, pass(self) :: randDSMarginal_
  procedure, pass(self) :: getDSDensityLBM_
  procedure, pass(self) :: VMFMeanDirDensity_
  procedure, pass(self) :: WatsonMeanDirDensity_
  procedure, pass(self) :: Estep_
  procedure, pass(self) :: Mstep_
  procedure, pass(self) :: Density_
  procedure, pass(self) :: logCp_
  procedure, pass(self) :: EMforDS_
  procedure, pass(self) :: getQandL_
  procedure, pass(self) :: setPGnum_
  procedure, pass(self) :: setNumEM_
  procedure, pass(self) :: setNumIter_
  procedure, pass(self) :: getN_
  procedure, pass(self) :: getNumEM_
  procedure, pass(self) :: getNumIter_
  procedure, pass(self) :: setQuatArray_
  procedure, pass(self) :: getQuatArray_
  final :: DirStats_destructor

  generic, public :: RotateToMu => RotateToMu_
  generic, public :: randUniformSphere => randUniformSphere_
  generic, public :: SampleDS => SampleDS_
  generic, public :: getDSDensityLBM => getDSDensityLBM_
  generic, public :: VMFMeanDirDensity => VMFMeanDirDensity_
  generic, public :: WatsonMeanDirDensity => WatsonMeanDirDensity_
  generic, public :: EMforDS => EMforDS_
  generic, public :: setNumEM => setNumEM_
  generic, public :: getNumEM => getNumEM_
  generic, public :: getN => getN_
  generic, public :: setNumIter => setNumIter_
  generic, public :: getNumIter => getNumIter_
  generic, public :: setPGnum => setPGnum_
  generic, public :: setQuatArray => setQuatArray_
  generic, public :: getQuatArray => getQuatArray_

end type DirStat_T

! the constructor routine for this class
interface DirStat_T
  module procedure DirStat_constructor
end interface DirStat_T

contains

!--------------------------------------------------------------------------
type(DirStat_T) function DirStat_constructor( DStype, PGnum ) result(DS)
!DEC$ ATTRIBUTES DLLEXPORT :: DirStat_constructor
 !! author: MDG
 !! version: 1.0
 !! date: 01/23/20
 !!
 !! constructor for the DirStat_T Class;

use mod_math
use mod_so3
use mod_quaternions

IMPLICIT NONE

character(3),INTENT(IN),OPTIONAL        :: DStype
integer(kind=irg),INTENT(IN),OPTIONAL   :: PGnum

type(QuaternionArray_T)                 :: qsym
type(QuaternionArray_T)                 :: qq
integer(kind=irg)                       :: i
real(kind=dbl)                          :: y1, y2

if (present(DStype)) then 
  DS%DStype = DStype

! von Mises-Fisher mode: (DStype='VMF')
! the next part of the initial Matlab code computes a lookup table for the parameter Ap(u) (Appendix in paper)
! this lookup table is only used when the ratio of the BesselI functions is between 0 and 0.95; for the
! region between 0.95 and 1, we use an analytical approximation (see VMF_Mstep routine).
!
! Watson mode: (DStype='WAT')
! we've used a similar approach to create a lookup table for values of kappa that are smaller than 35, in
! which case we use the standard ratio of Kummer functions:  Kummer[3/2,3,k]/Kummer[1/2,2,k]/k.  For
! larger kappa values, we have an expansion using the large argument behavior of the modified Bessel functions.
!

! allocate the parameter arrays
  DS%Apnum = 35000
  allocate(DS%xAp(DS%Apnum), DS%yAp(DS%Apnum))

! define the xAp array
  DS%xAp = (/ (0.001D0+dble(i-1)*0.001D0,i=1,DS%Apnum)  /)

  if (DS%DStype.eq.'VMF') then ! von Mises-Fisher distribution
    do i=1,DS%Apnum
      DS%yAp(i) = BesselIn(DS%xAp(i), 2) / BesselI1(DS%xAp(i))
    end do
  end if

  if (DS%DStype.eq.'WAT') then ! Watson distribution
    do i=1,DS%Apnum
      y1 = BesselI1(DS%xAp(i)*0.5D0)
      y2 = BesselI0(DS%xAp(i)*0.5D0)
      DS%yAp(i) = y1 / (y2-y1) / DS%xAp(i)
    end do
  end if
end if 

! do we need to initialize the list of quaternion symmetry operators ?
if (present(PGnum)) then
  DS%pgnum = PGnum
  call qq%QSym_Init(PGnum, qsym)
  DS%qsym = qsym
  DS%Pmdims = DS%qsym%getQnumber()
end if

end function DirStat_constructor

!--------------------------------------------------------------------------
subroutine DirStats_destructor( self )
!DEC$ ATTRIBUTES DLLEXPORT :: DirStats_destructor
 !! author: MDG
 !! version: 1.0
 !! date: 01/24/20
 !!
 !! deallocate all arrays

IMPLICIT NONE

type(DirStat_T), INTENT(INOUT)      :: self

call reportDestructor('DirStat_T')

if (allocated(self%xAp)) deallocate(self%xAp)
if (allocated(self%yAp)) deallocate(self%yAp)

end subroutine DirStats_destructor

!--------------------------------------------------------------------------
function getNumEM_( self ) result(NumEM)
!DEC$ ATTRIBUTES DLLEXPORT :: getNumEM_
 !! author: MDG
 !! version: 1.0
 !! date: 01/24/20
 !!
 !! get NumEM parameter

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)       :: self
integer(kind=irg)                     :: NumEM

NumEM = self%NumEM

end function getNumEM_

!--------------------------------------------------------------------------
function getNumIter_( self ) result(NumIter)
!DEC$ ATTRIBUTES DLLEXPORT :: getNumIter_
 !! author: MDG
 !! version: 1.0
 !! date: 01/24/20
 !!
 !! get NumIter parameter

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)       :: self
integer(kind=irg)                     :: NumIter

NumIter = self%NumIter

end function getNumIter_

!--------------------------------------------------------------------------
function getN_( self ) result(N)
!DEC$ ATTRIBUTES DLLEXPORT :: getN_
 !! author: MDG
 !! version: 1.0
 !! date: 01/24/20
 !!
 !! get N parameter

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)       :: self
integer(kind=irg)                     :: N

N = self%N

end function getN_

!--------------------------------------------------------------------------
subroutine setNumEM_( self, NumEM )
!DEC$ ATTRIBUTES DLLEXPORT :: setNumEM_
 !! author: MDG
 !! version: 1.0
 !! date: 01/24/20
 !!
 !! set NumEM parameter

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)       :: self
integer(kind=irg),INTENT(IN)          :: NumEM

self%NumEM = NumEM

end subroutine setNumEM_

!--------------------------------------------------------------------------
subroutine setNumIter_( self, NumIter )
!DEC$ ATTRIBUTES DLLEXPORT :: setNumIter_
 !! author: MDG
 !! version: 1.0
 !! date: 01/24/20
 !!
 !! set NumIter parameter

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)       :: self
integer(kind=irg),INTENT(IN)          :: NumIter

self%NumIter = NumIter

end subroutine setNumIter_

!--------------------------------------------------------------------------
subroutine setPGnum_( self, PGnum )
!DEC$ ATTRIBUTES DLLEXPORT :: setPGnum_
 !! author: MDG
 !! version: 1.0
 !! date: 01/24/20
 !!
 !! set PGnum parameter

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)       :: self
integer(kind=irg),INTENT(IN)          :: PGnum

self%PGnum = PGnum

end subroutine setPGnum_

!--------------------------------------------------------------------------
subroutine setQuatArray_( self, qAR, slot )
!DEC$ ATTRIBUTES DLLEXPORT :: setQuatArray_
 !! author: MDG
 !! version: 1.0
 !! date: 01/24/20
 !!
 !! put a quaternion array in Xquats or Xaux (slot='aux')

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)       :: self
type(QuaternionArray_T), INTENT(INOUT):: qAR
character(*), INTENT(IN), OPTIONAL    :: slot

if (present(slot)) then
  if (trim(slot).eq.'aux') self%Xaux = qAR
else
  self%Xquats = qAR
end if

self%N = qAR%getQnumber()

end subroutine setQuatArray_

!--------------------------------------------------------------------------
function getQuatArray_( self, slot ) result(qAR)
!DEC$ ATTRIBUTES DLLEXPORT :: getQuatArray_
 !! author: MDG
 !! version: 1.0
 !! date: 01/24/20
 !!
 !! get a quaternion array from Xquats or Xaux (slot='aux')

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)       :: self
character(*), INTENT(IN), OPTIONAL    :: slot
type(QuaternionArray_T)               :: qAR

if (present(slot)) then
  if (trim(slot).eq.'aux') qAR = self%Xaux
else
  qAR = self%Xquats
end if

end function getQuatArray_

!--------------------------------------------------------------------------
recursive function RotateToMu_(self, lmu, y) result(ymu)
!DEC$ ATTRIBUTES DLLEXPORT :: RotateToMu_
 !! author: MDG
 !! version: 1.0
 !! date: 01/23/20
 !!
 !! Rotate an array of quaternions to an average direction lmu using the null space approach
 !!
 !! In Matlab, one uses the null() operator which returns the null space of the argument
 !! This is then inserted into a 4x4 rotation matrix and multiplied with the quaternions
 !! from the random sample.  The null space of the input quaternion can be computed with
 !! singular value decomposition, which is done with the dgesvd Lapack routine. The matrix
 !! returned as u is the desired rotation matrix, except that the numbers in the first
 !! column must have their signs reversed.

use mod_quaternions
use mod_rotations

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)         :: self
type(Quaternion_T), INTENT(IN)          :: lmu
type(QuaternionArray_T), INTENT(INOUT)  :: y
type(QuaternionArray_T)                 :: ymu

integer(kind=irg)                       :: i

! parameters for the singular value decomposition
integer(kind=irg)                       :: nr, LDA, LDU, LDVT, lwork, info, N
real(kind=dbl)                          :: mA(4,4), ss(4), u(4,4), vt, work(20)
type(Quaternion_T)                      :: qq

! Rotate the distribution along the desired mean direction mu

mA = 0.D0
mA(1:4,1) = lmu%get_quatd()
nr = 4
LDA = 4
LDVT = 1
LDU = 4
lwork = 20
call DGESVD('A','N',nr,nr,mA,LDA,ss, u, LDU, vt, LDVT, work, lwork, info)
u(1:4,1) = -u(1:4,1)

! next, apply this 4x4 rotation matrix to all of the generated quaternions to
! rotate them along the mean direction mu
N = y%getQnumber()
ymu = QuaternionArray_T( n = N, s = 'd' )
do i=1,N
  qq = y%getQuatfromArray(i)
  call ymu%insertQuatinArray(i, quaternion_T( qd = matmul(u, qq%get_quatd() ) ) )
end do

end function RotateToMu_

!--------------------------------------------------------------------------
recursive function randUniformSphere_(self, N, seed) result(ranSphere)
!DEC$ ATTRIBUTES DLLEXPORT :: randUniformSphere_
 !! author: MDG, based on 2015 Yu-Hui's Matlab code
 !! version: 1.0
 !! date: 01/23/20
 !!
 !! Return a set of random vectors on the sphere S^2 using normal random sampling

use mod_math

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)      :: self
integer(kind=irg),INTENT(IN)            :: N
integer(kind=irg),INTENT(INOUT)         :: seed
real(kind=dbl)                          :: ranSphere(3,N)

real(kind=dbl)                          :: nq, NR(N*3), randNorm(3,N)
integer(kind=irg)                       :: i

ranSphere = 0.D0
call R8VEC_normal_01(N*3, seed, NR)
randNorm = reshape( NR, (/ 3, N /) )

! and normalize the three-vectors
do i=1,N
  nq = vecnorm(randNorm(1:3,i))
  RanSphere(1:3,i) = randNorm(1:3,i)/nq
end do

end function randUniformSphere_

!--------------------------------------------------------------------------
recursive function SampleDS_(self, N, seed, mu, kappa) result(sDS)
!DEC$ ATTRIBUTES DLLEXPORT :: SampleDS_
 !! author: MDG, based on 2015 Yu-Hui's Matlab code
 !! version: 1.0
 !! date: 01/23/20
 !!
 !! Sample a directional distribution (VMF or WAT) on the quaternion unit sphere

use mod_io
use mod_quaternions

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)      :: self
integer(kind=irg),INTENT(IN)      :: N
 !! number of samples to return
integer(kind=irg),INTENT(INOUT)   :: seed
 !! random number generator seed
type(Quaternion_T),INTENT(INOUT)  :: mu
 !! mean direction (unit quaternion)
real(kind=dbl),INTENT(IN)         :: kappa
 !! concentration parameter
type(QuaternionArray_T)           :: sDS
 !! return array

type(IO_T)                        :: Message
real(kind=dbl)                    :: nq, tmpmu(4), RandSphere(3,N), t(N), RS(4,N), y(4,N)
type(Quaternion_T)                :: lmu
integer(kind=irg)                 :: i
type(QuaternionArray_T)           :: yAR

! make sure the input quaternion is normalized
nq = mu%quat_norm()
if (nq.eq.0.D0) call Message%printError('SampleDS','Input quaternion has zero length')
lmu = mu/nq

! initialize a bunch of parameters
tmpmu = (/ 1.D0, 0.D0, 0.D0, 0.D0 /)
RS = 0.D0

! get the t-parameter
t = randDSMarginal_(self, N, kappa, seed)

! and the distribution of random directions on the 2-sphere
RandSphere = self%randUniformSphere(N,seed)
RS(2:4,1:N) = RandSphere(1:3,1:N)

! merge these two parameters into the desired random variables
y = transpose( spread(t,DIM=2,NCOPIES=4) * spread(tmpmu,DIM=1,NCOPIES=N) + &
               spread(dsqrt(1.D0-t*t),DIM=2,NCOPIES=4) * transpose(RS) )

! copy the y array into a QuaternionArray
yAR = QuaternionArray_T( n=N, qd=y )

! Rotate the distribution along the desired mean direction mu
sDS = self%RotateToMu(lmu, yAR)

end function SampleDS_

!--------------------------------------------------------------------------
recursive function randDSMarginal_(self, N, k, seed) result(t)
!DEC$ ATTRIBUTES DLLEXPORT :: randDSMarginal_
 !! author: MDG, based on 2015 Yu-Hui's Matlab code
 !! version: 1.0
 !! date: 01/23/20
 !!
 !! rejection sampling for the t parameter for the marginal directional distribution
 !!
 !! This algorithm samples the parameter t from a marginal distribution f(t), which
 !! by itself is the VMF or Watson distribution integrated around the mean direction.  These
 !! expressions were derived by Yu-Hui and verified by MDG [02/05/15] and then implemented;
 !! there are some typographical errors in the literature, and the versions documented here
 !! are correct [numerical verification].

use mod_math

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)      :: self
integer(kind=irg),INTENT(IN)            :: N
real(kind=dbl),INTENT(IN)               :: k
integer(kind=irg),INTENT(INOUT)         :: seed
real(kind=dbl)                          :: t(N)

real(kind=dbl)                          :: LBM(2), h, x, C
integer(kind=irg)                       :: i

! Find the left bound and maximum needed for rejection sampling
C = 0.D0
LBM = self%getDSDensityLBM(k, C)

! apply the rejection sampling algorithm to either VMF or Watson marginal distributions
t = 0.D0
do i=1,N
  do
    x = r8_uniform_01(seed)*(1.D0-LBM(1))+LBM(1)
    if (self%DStype.eq.'VMF') h = self%VMFMeanDirDensity(x, k, C)
    if (self%DStype.eq.'WAT') h = self%WatsonMeanDirDensity(x, k, C)
    if (r8_uniform_01(seed)*LBM(2).le.h) EXIT
  end do
  t(i) = x
end do

end function randDSMarginal_

!--------------------------------------------------------------------------
recursive function getDSDensityLBM_(self, k, C) result(LBM)
!DEC$ ATTRIBUTES DLLEXPORT :: getDSDensityLBM_
 !! author: MDG, based on 2015 Yu-Hui's Matlab code
 !! version: 1.0
 !! date: 01/23/20
 !!
 !! determines the left bound and maximum for rejection sampling

use mod_io

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)      :: self
real(kind=dbl),INTENT(IN)               :: k
real(kind=dbl),INTENT(INOUT)            :: C
real(kind=dbl)                          :: LBM(2)

type(IO_T)                              :: Message
real(kind=dbl),parameter                :: min_thresh=0.00001D0
real(kind=dbl)                          :: s, x, start
integer(kind=irg)                       :: f

start = -1.D0
if (self%DStype.eq.'WAT') start = 0.D0

! first we look for the left bound
f = 1
do
 x = start+dble(f)*0.00001D0
 if (x.eq.1.D0) call Message%printError('getDSDensityLBM','reached +1 in leftbound determination')
 if (self%DStype.eq.'VMF') s = self%VMFMeanDirDensity(x,k,C)
 if (self%DStype.eq.'WAT') s = self%WatsonMeanDirDensity(x,k,C)
 if (s.ge.min_thresh) EXIT
 f = f+1
end do
!
LBM(1) =  start+dble(f)*0.00001D0

if (self%DStype.eq.'VMF') then
! for the simplified version of the density function, we have an analytical
! expression for where the maximum of the function occurs [convert the BesselI(3/2,x)
! to hyperbolic functions, then to exponential, and ignore the negative exponential
! which will be very small for reasonably sized k and t...]
  x = (-1.D0+dsqrt(1.D0+4.D0*k*k))/(2.D0*k)
  LBM(2) = self%VMFMeanDirDensity(x,k,C)
end if

if (self%DStype.eq.'WAT') then
! for the simplified version of the density function, we have an analytical
! expression for where the maximum of the function occurs
  x = dsqrt((2.D0*k-1.D0)/(2.D0*k))
  LBM(2) = self%WatsonMeanDirDensity(x,k,C)
end if

end function getDSDensityLBM_

!--------------------------------------------------------------------------
recursive function VMFMeanDirDensity_(self, x, k, C) result(y)
!DEC$ ATTRIBUTES DLLEXPORT :: VMFMeanDirDensity_
 !! author: MDG, based on 2015 Yu-Hui's Matlab code
 !! version: 1.0
 !! date: 01/23/20
 !!
 !! function to be sampled for VMF random sampling; we're using a close approximation

use mod_math
use mod_io

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)      :: self
real(kind=dbl),INTENT(IN)       :: x
real(kind=dbl),INTENT(IN)       :: k
real(kind=dbl),INTENT(INOUT)    :: C
real(kind=dbl)                  :: y

type(IO_T)                      :: Message

if (dabs(x).gt.1.D0) call Message%printError('VMFMeanDirDensity','argument must be in [-1,1]')

! explicit expression for p=4 (Gamma[3/2]Gamma[1/2] = pi/2)
! and the BesselI(3/2) function reduces to hyperbolic functions
! diverges for k->0, and becomes really small for large k
if (C.eq.0.D0) then
  C = 2.D0*k**(2.5D0)/dsqrt(2.D0*cPi)/(k-1.D0)
end if

! this is a close approximation, really good for larger values of k
! and numerically more stable than the original, which has problems for k>600 or so
y = C * dexp(k*(x-1.D0))*dsqrt(1.D0-x*x)

end function VMFMeanDirDensity_

!--------------------------------------------------------------------------
recursive function WatsonMeanDirDensity_(self, x, k, C) result(y)
!DEC$ ATTRIBUTES DLLEXPORT :: WatsonMeanDirDensity_
 !! author: MDG, based on 2015 Yu-Hui's Matlab code
 !! version: 1.0
 !! date: 01/23/20
 !!
 !! function to be sampled for Watson random sampling; we're using a close approximation

use mod_math
use mod_io

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)      :: self
real(kind=dbl),INTENT(IN)       :: x
real(kind=dbl),INTENT(IN)       :: k
real(kind=dbl),INTENT(INOUT)    :: C
real(kind=dbl)                  :: y

real(kind=dbl),parameter        :: CC = 144.43253338822560946D0  ! 256/sqrt(pi)
type(IO_T)                      :: Message

if (dabs(x).gt.1.D0) call Message%printError('WatsonMeanDirDensity','argument must be in [-1,1]')

! approximate expression for p=4
if (C.eq.0.D0) then
  C = CC*k**4.5D0/(525.D0+4.D0*k*(45.D0+8.D0*k*(3.D0+4.D0*k)))
end if

! this is a close approximation, really good for larger values of k
! and numerically more stable than the original, which has problems for k>600 or so
y = C * dexp(k*(x*x-1.D0))*dsqrt(1.D0-x*x)

end function WatsonMeanDirDensity_

!--------------------------------------------------------------------------
recursive subroutine EMforDS_(self, seed, muhat, kappahat)
!DEC$ ATTRIBUTES DLLEXPORT :: EMforDS_
 !! author: MDG, based on 2015 Chen's Matlab code, with simplifications
 !! version: 1.0
 !! date: 01/23/20
 !!
 !! Expectation maximization approach to maximum likelihood problem for mu and kappa
 !!
 !! this routine expects the input quaternion array to be stored in Xquats using the setQuatArray method

use mod_math, only:r8vec_normal_01, r4_uniform_01  ! array of normal random numbers
use mod_quaternions
use mod_rotations
use mod_so3

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)      :: self
integer(kind=irg),INTENT(INOUT)      :: seed
type(Quaternion_T),INTENT(OUT)       :: muhat
real(kind=dbl),INTENT(OUT)           :: kappahat

type(so3_T)                          :: SO
integer(kind=irg)                    :: i, j, N, Pmdims, init, dd, NumEM, NumIter
real(kind=dbl),allocatable           :: Mu_All(:,:), Kappa_All(:), R_All(:,:,:), L_All(:), &
                                        R(:,:), Q(:), L(:)
real(kind=dbl)                       :: MuKa(5), Qi, Li, rod(4), Kappa, x(4)
type(Quaternion_T)                   :: Mu, PmMu, qu
type(q_T)                            :: MuMu
type(r_T)                            :: roFZ

! In this routine, we perform the EM algorithm to obtain an estimate for the
! mean direction and concentration parameter of the modified von Mises-Fisher (mVMF)
! distribution that models the statistics of the orientation point cloud.

! array sizes
N = self%getN()
Pmdims = self%qsym%getQnumber()
NumEM = self%NumEM
NumIter = self%NumIter

! initialize some auxiliary arrays
allocate(Mu_All(NumEM,4), Kappa_All(NumEM), &
         R_All(N,Pmdims,NumEM),L_All(NumIter))
Mu_All = 0.D0
Kappa_All = 0.D0
R_All = 0.D0
L_All = 0.D0


! main loop (EM typically uses a few starting parameter sets to make sure we don't get stuck in a local maximum)
do init=1,NumEM
! create a vector to hold the results
  allocate(R(N,Pmdims))
  R = 0.D0

! generate a normal random vector and normalize it as a starting guess for Mu (i.e., a unit quaternion)
  call R8VEC_normal_01(4,seed,x)
  Mu = Quaternion_T( qd = x )
  call Mu%quat_normalize()
  call Mu%quat_pos()

! starting value for Kappa
  Kappa = 30.D0

! define the number of iterations and the Q and L function arrays
  allocate (Q(NumIter), L(NumIter))
  Q = 0.D0
  L = 0.D0

! and here we go with the EM iteration...
! we use quaternion multiplication throughout instead of the matrix version in the Matlab version
! quaternion multiplication has been verified against the 4x4 matrix multiplication of the Matlab code on 01/02/15
  iloop: do i=1,NumIter
! E-step
    R = Estep_(self, Mu,Kappa)
! M-step
    MuKa = Mstep_(self, R)
! calculate the Q and Likelihood function values
    call getQandL_(self, MuKa,R,Qi,Li)
    L(i) = Li
    Q(i) = Qi

! update the containers
    Mu_All(init,1:4) = MuKa(1:4)
    Kappa_All(init) = MuKa(5)
    R_All(1:N,1:Pmdims,init) = R(1:N,1:Pmdims)
    L_All(init) = L(i)
    Mu = Quaternion_T( qd = MuKa(1:4) )
    Kappa = MuKa(5)

! and terminate if necessary
    if (i.ge.2) then
      if (abs(Q(i)-Q(i-1)).lt.0.01) then
        EXIT iloop
      end if
    end if
  end do iloop
  deallocate(R,Q,L)
end do

dd = maxloc(L_All,1)
Mu = Quaternion_T( qd = Mu_all(dd,1:4) )
call Mu%quat_pos()
kappahat = Kappa_All(dd)

! the final step is to make sure that the resulting Mu lies in the fundamental zone.
! since we start the EM iterations from a random quaternion, there is no guarantee that the
! result lies in the same fundamental zone. Therefore, we cycle through all the
! equivalent quaternions, and stop as soon as we find one in the Rodrigues
! fundamental zone, which requires routines from the rotations and so3 modules.
SO = so3_T( self%pgnum )
MuMu = q_T( qdinp = Mu%get_quatd() )
call SO%ReduceOrientationtoRFZ( MuMu, self%qsym, roFZ )
MuMu = roFZ%rq()
muhat = Quaternion_T( qd = MuMu%q_copyd() )

deallocate(Mu_All, Kappa_All, R_All, L_All)

end subroutine EMforDS_

!--------------------------------------------------------------------------
recursive function Estep_(self, Mu, Kappa) result(R)
!DEC$ ATTRIBUTES DLLEXPORT :: Estep_
 !! author: MDG, based on 2015 Chen's Matlab code, with simplifications
 !! version: 1.0
 !! date: 01/24/20
 !!
 !! computes the E step of the EM process, verified against Matlab code on 01/02/15

use mod_quaternions

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)         :: self
type(Quaternion_T),INTENT(IN)           :: Mu
 !! current guess for mean quaternion
real(kind=dbl),INTENT(IN)               :: Kappa
 !! current guess for concentration parameter
real(kind=dbl)                          :: R(self%N,self%Pmdims)

integer(kind=irg)                       :: j
real(kind=dbl)                          :: Rdenom(self%N), C
type(Quaternion_T)                      :: PmMu

C = logCp_(self, kappa)

do j=1,self%Pmdims
  PmMu = Mu * self%qsym%getQuatfromArray(j)
  R(1:self%N,j) = Density_(self, PmMu%get_quatd(), Kappa, C)
end do
! and determine the normalization factors
Rdenom = 1.D0/sum(R,2)

do j=1,self%Pmdims
  R(1:self%N,j) = R(1:self%N,j)*Rdenom(1:self%N)
end do

end function Estep_

!--------------------------------------------------------------------------
recursive function Mstep_(self, R) result(MuKa)
!DEC$ ATTRIBUTES DLLEXPORT :: Mstep_
 !! author: MDG, based on 2015 Chen's Matlab code, with simplifications
 !! version: 1.0
 !! date: 01/24/20
 !!
 !! computes the M step of the EM process

use mod_quaternions
use mod_math

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)         :: self
real(kind=dbl),INTENT(IN)               :: R(self%N,self%Pmdims)
 !! weight factors from the E step
real(kind=dbl)                          :: MuKa(5)

type(Quaternion_T)                      :: qu
real(kind=dbl)                          :: tmpGamma(4), nGamma, diff(self%Apnum), y, Tscatt(4,4), tmp(4,4), x(4)
integer(kind=irg)                       :: minp, i, j

! variables needed for the dsyev Lapack eigenvalue routine
CHARACTER                               :: JOBZ, UPLO
INTEGER                                 :: INFO, LDA, LWORK, NN
DOUBLE PRECISION                        :: A( 4 , 4 ), W( 4 ), WORK( 20 )

if (self%DStype.eq.'VMF') then
! this is simplified from the Matlab routine and uses straight summations and
! quaternion multiplication instead of arrays
  tmpGamma = 0.D0
  do j=1,self%Pmdims
    do i=1,self%N
      qu = self%Xquats%getQuatfromArray(i) * conjg(self%qsym%getQuatfromArray(j) )
      tmpGamma = tmpGamma +  R(i,j) * qu%get_quatd()
    end do
  end do
  nGamma = vecnorm(tmpGamma)
  MuKa(1:4) = tmpGamma/nGamma
  y = nGamma/dble(self%N)
end if

if (self%DStype.eq.'WAT') then
! here, we compute the modified scattering matrix Tscatt and compute its largest eigenvalue
! under the assumption that kappa will always be positive for the types of problems that we
! need to consider; we need to use the outer product, implemented using spread calls
  Tscatt = 0.D0
  do j=1,self%Pmdims
    do i=1,self%N
      qu = self%Xquats%getQuatfromArray(i) * conjg(self%qsym%getQuatfromArray(j) )
      x = qu%get_quatd()
      tmp = spread(x,dim=2,ncopies=4)*spread(x,dim=1,ncopies=4)
      Tscatt = Tscatt + R(i,j) * tmp
    end do
  end do
  Tscatt = Tscatt/dble(self%N)

  JOBZ = 'V'
  UPLO = 'U'
  NN = 4
  LDA = 4
  LWORK = 20
  A = Tscatt
! DSYEV computes all eigenvalues and, optionally, eigenvectors of a real symmetric matrix A.
  call DSYEV( JOBZ, UPLO, NN, A, LDA, W, WORK, LWORK, INFO )
  x(1:4) = A(1:4,4)
  MuKa(1:4) = x(1:4)
  y = dot_product(x,matmul(Tscatt,x))
end if

! find kappa corresponding to this value of gamma (equation 17 in appendix of paper)
! we split this into two regionds: 0<=y<0.94, for which we use the look-up table
! approach, and 0.94<=y<=1, for which we have derived an analytical approximation
! that is pretty accurate in the relevant region of kappa>30.
if (y.ge.0.94D0) then
  if (self%DStype.eq.'VMF') MuKa(5) = (15.D0-3.D0*y+dsqrt(15.D0+90.D0*y+39.D0*y*y))/(16.D0*(1.0D0-y))
  if (self%DStype.eq.'WAT') MuKa(5) = (5.D0*y-11.D0-dsqrt(39.D0-12.D0*y+9.D0*y**2))/(8.D0*(y-1.D0))
else
  diff = dabs( y - self%yAp )
  minp = minloc( diff, 1 )
  if (minp.eq.1) minp = 2
  MuKa(5) = self%xAp(minp)
end if

end function Mstep_

!--------------------------------------------------------------------------
recursive subroutine getQandL_(self, MuKa, R, Q, L)
!DEC$ ATTRIBUTES DLLEXPORT :: getQandL_
 !! author: MDG, based on 2015 Chen's Matlab code, with simplifications
 !! version: 1.0
 !! date: 01/24/20
 !!
 !! Computes the Q array and the log-likelihood array

use mod_quaternions

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)         :: self
real(kind=dbl),INTENT(IN)               :: MuKa(5)
 !! vector with mu and kappa values
real(kind=dbl),INTENT(IN)               :: R(self%N,self%Pmdims)
 !! output from the E step
real(kind=dbl),INTENT(INOUT)            :: Q
 !! output Q
real(kind=dbl),INTENT(INOUT)            :: L
 !! output L

real(kind=dbl)                          :: Phi(self%N,self%Pmdims), C, oldQ, oldL
type(Quaternion_T)                      :: qu, PmMu
integer(kind=irg)                       :: j
real(kind=dbl),parameter                :: eps = 0.00001D0

  oldQ = Q
  oldL = L

  C = logCp_(self, MuKa(5))
  if (self%DStype.eq.'VMF') C = dexp(C)

! compute the auxiliary Phi array
  Phi = 0.D0
  qu = Quaternion_T( qd = MuKa(1:4) )
  do j=1,self%Pmdims
    PmMu = self%qsym%getQuatfromArray(j) * qu
    Phi(1:self%N,j) = Density_(self, PmMu%get_quatd(), MuKa(5), C)
  end do
  Phi = Phi/dble(self%Pmdims)
  if (minval(Phi).gt.0.D0) then
! and convert the array into the Q and L parameters.
   L = sum(dlog(sum(Phi,2)))
   Q = sum(R*dlog(Phi))
  else
   L = oldL
   Q = oldQ
  end if ! else, we reuse the old values
end subroutine getQandL_

!--------------------------------------------------------------------------
recursive function Density_(self,  mu, kappa, C) result(y)
!DEC$ ATTRIBUTES DLLEXPORT :: Density_
 !! author: MDG, based on 2015 Chen's Matlab code, with simplifications
 !! version: 1.0
 !! date: 01/24/20
 !!
 !! computes the VMF or Watson density function
 !!
 !! original in Matlab by Yu-Hui Chen, U. Michigan
 !! converted to IDL by MDG, 12/18/14, simplified arguments
 !! converted to f90 by MDG, 12/31/14, further simplifications
 !! output validated against Matlab output on 12/31/14

use mod_quaternions

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT)         :: self
real(kind=dbl),INTENT(IN)               :: mu(4)
 !! mean direction
real(kind=dbl),INTENT(IN)               :: kappa
 !! concentration parameter
real(kind=dbl),INTENT(IN)               :: C
 !! logCp(kappa) or exp(logCp(kappa) (precomputed in calling routine)
real(kind=dbl)                          :: y(self%N)

integer(kind=irg)                       :: j
type(Quaternion_T)                      :: q

if (self%DStype.eq.'VMF') then
 do j=1,self%N
  q = self%Xquats%getQuatfromArray(j)
  y(j) = dexp(C+kappa*dot_product(mu, q%get_quatd() ) )
 end do
end if

if (self%DStype.eq.'WAT') then
 do j=1,self%N
  q = self%Xquats%getQuatfromArray(j)
  y(j) = dexp(C+kappa*dot_product(mu, q%get_quatd() )**2)
 end do
end if

end function Density_

!--------------------------------------------------------------------------
!
! FUNCTION: logCp
!
!> @author Marc De Graef, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
!
!> @brief computes the logarithm of Cp for VMF and Watson distributions
!
!> @details
!> original in Matlab by Yu-Hui Chen, U. Michigan
!> converted to IDL by MDG, 12/18/14, simplified arguments
!> converted to f90 by MDG, 12/31/14, further simplifications
!> output validated against Matlab output on 12/31/14
!
!> @param kappa input parameter
!> @param Dtype 'VMF' or 'WAT'
!
!> @date 12/31/14 MDG 1.0 original
!> @date 01/09/14 MDG 1.1 introduced more accurate numerical approximation for Cp
!> @date 01/09/14 MDG 1.2 moved some of the constants in front of the logarithm
!> @date 02/06/15 MDG 1.3 merged VMF and WAT routines
!--------------------------------------------------------------------------
recursive function logCp_(self, kappa) result(lCp)
!DEC$ ATTRIBUTES DLLEXPORT :: logCp_
 !! author: MDG, based on 2015 Chen's Matlab code, with simplifications
 !! version: 1.0
 !! date: 01/24/20
 !!
 !! computes the logarithm of Cp for VMF and Watson distributions

use mod_math

IMPLICIT NONE

class(DirStat_T), INTENT(INOUT) :: self
real(kind=dbl),INTENT(IN)       :: kappa
real(kind=dbl)                  :: lCp

! pre-computed constants
real(kind=dbl),parameter        :: C=-3.675754132818690967D0    ! C = ln(1.D0/(2.D0*cPi)**2)
real(kind=dbl),parameter        :: C2=4.1746562059854348688D0   ! C2 = ln(512/sqrt(2)/pi^(3/2))
real(kind=dbl),parameter        :: C2W=5.4243952068443172530D0  ! C2 = ln(128*sqrt(pi))

if (self%DStype.eq.'VMF') then
! for arguments larger than kappa=30, we use a simple numerical approximation
  if (kappa.gt.30.D0) then
    lCp = kappa**4.5D0/(-105D0+8.D0*kappa*(-15.D0+16.D0*kappa*(-3.D0+8.D0*kappa)))
    lCp = C2 - kappa + dlog(lCp)
  else
    lCp = C + dlog( kappa / BesselI1(kappa) )
  end if
end if

if (self%DStype.eq.'WAT') then
! for arguments larger than kappa=20, we use a simple numerical approximation
  if (kappa.gt.20.D0) then
    lCp = kappa**4.5D0/(525.D0 + 4.D0*kappa*(45.D0 + 8.D0*kappa*(3.D0 + 4.D0*kappa)))
    lCp = C2W - kappa + dlog(lCp)
  else
    lCp = -kappa*0.5D0 - dlog( BesselI0(kappa*0.5D0) - BesselI1(kappa*0.5D0) )
  end if
end if

end function logCp_

! !--------------------------------------------------------------------------
! !--------------------------------------------------------------------------
! !--------------------------------------------------------------------------
! ! the following routines have to do with the neighborhood similarity analysis;
! ! however, they are currently not used anywhere in EMsoft, so we comment them out.
! !--------------------------------------------------------------------------
! !--------------------------------------------------------------------------
! !--------------------------------------------------------------------------
! !
! ! SUBROUTINE: DI_similarity_classifier
! !
! !> @author Saransh Singh, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
! !
! !> @brief classify the point as grain interior or anomalous point
! !
! !> @details takes the kNN neighbor information as input and returns the
! !  whether the point lies in the interior of the grain or lies on the
! !  grain boundary. Details in pg 11 of the selfionary indexing paper
! !
! !> @param array input array
! !> @param k number of top matches for each pixel
! !> @param npx number of pixels in the x direction
! !> @param npy number of pixels in the y direction
! !
! !> @date 01/05/15 SS 1.0 original
! !> @date 01/06/15 MDG 1.1 simplified summation loop and renamed routine
! !--------------------------------------------------------------------------
! recursive subroutine DI_Similarity_Classifier(array,k,npx,npy,returnarr)
! !DEC$ ATTRIBUTES DLLEXPORT :: DI_Similarity_Classifier

! use local

! IMPLICIT NONE

! integer(kind=irg),INTENT(IN)            :: k
! integer(kind=irg),INTENT(IN)            :: npx
! integer(kind=irg),INTENT(IN)            :: npy
! integer(kind=sgl),INTENT(IN)            :: array(npx,npy,k)
! real(kind=sgl),INTENT(OUT)              :: returnarr(npx,npy)

! integer(kind=irg)                       :: ii,jj,ki,kj,similarity_measure_sum,res
! real(kind=sgl)                          :: similarity_measure


! similarity_measure_sum = 0
! similarity_measure = 0.0
! returnarr = 0.0

! do ii = 2,npx-1
!     do jj = 2,npy-1
!       do ki = -1, 1
!         do kj = -1, 1
!           if ((abs(ki)+abs(kj)).ne.0) then
!             call CardIntersection(array(ii+ki,jj+kj,1:k),array(ii,jj,1:k),k,res)
!             similarity_measure_sum = similarity_measure_sum + res

!             similarity_measure = float(similarity_measure_sum)/float(8*k)
!             returnarr(ii,jj) = similarity_measure

!             similarity_measure = 0.0
!             similarity_measure_sum = 0
!           end if
!         end do
!       end do
!     end do
! end do

! end subroutine DI_Similarity_Classifier

! !--------------------------------------------------------------------------
! !
! ! SUBROUTINE: CardIntersection
! !
! !> @author Saransh Singh, Carnegie Mellon University / Yu-Hui Chen, U. Michigan
! !
! !> @brief calculate the cardinality of the intersection of two sets
! !
! !> @param set1
! !> @param set2
! !> @param k number of elements in each set
! !
! !> @date 01/05/15 SS 1.0 original
! !> @date MDG 1.1 changed types to integer
! !--------------------------------------------------------------------------
! recursive subroutine CardIntersection(set1,set2,k,res)
! !DEC$ ATTRIBUTES DLLEXPORT :: CardIntersection

! use local

! IMPLICIT NONE

! integer(kind=irg),INTENT(IN)            :: k
! integer(kind=irg),INTENT(IN)            :: set1(k)
! integer(kind=irg),INTENT(IN)            :: set2(k)
! integer(kind=irg),INTENT(OUT)           :: res

! integer(kind=irg)                       :: ii,jj
! jj = 1
! res = 0

! do ii = 1,k
!     do jj = 1,k
!         if (set1(ii) .eq. set2(jj)) then
!             res = res + 1
!             EXIT
!         end if
!     end do
! end do

! end subroutine CardIntersection


end module mod_dirstats
