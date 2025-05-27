! ###################################################################
! Copyright (c) 2013-2025, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_cluster
  !! author: MDG 
  !! version: 1.0 
  !! date: 05/22/25
  !!
  !! class definition for the cluster module

use mod_kinds
use mod_global
use mod_DIfiles

IMPLICIT NONE 

! class definition
type, public :: Cluster_T
  integer(kind=irg), allocatable  :: grainID(:,:) ! grain identifier
  integer(kind=irg), allocatable  :: npixels(:)   ! number of pixels in each grain
  integer(kind=irg), allocatable  :: ROI(:,:)     ! bounding box information for each grain
  integer(kind=irg)               :: nGrains      ! number of grains found
  integer(kind=irg)               :: GrainSize    ! used to control recursion
  integer(kind=irg)               :: TotSize      ! used to control recursion
  integer(kind=irg)               :: MaxSize      ! used to control recursion
  integer(kind=irg)               :: ipf_wd       ! used to control recursion
  integer(kind=irg)               :: ipf_ht       ! used to control recursion
  real(kind=sgl)                  :: gangle       ! threshold angle for clustering
  real(kind=dbl), allocatable     :: avor(:,:)    ! average orientation per grain
  real(kind=dbl), allocatable     :: kappa(:)     ! concentration parameters
  real(kind=sgl), allocatable     :: kam(:,:)     ! kernel average orientation map

contains
private 
  procedure, pass(self) :: ScanGrain_
  procedure, pass(self) :: getROI_
end type Cluster_T

! the constructor routine for this class 
interface Cluster_T
  module procedure cluster_constructor
end interface Cluster_T

contains

!--------------------------------------------------------------------------
type(Cluster_T) function cluster_constructor( DIFT, gangle ) result(cluster)
!! author: MDG 
!! version: 1.0 
!! date: 05/22/25
!!
!! constructor for the Cluster_T Class; this routine performs the entire computation
 
use mod_DIsupport
use mod_io
use mod_dirstats 
use mod_quaternions
use mod_so3
use mod_rotations

IMPLICIT NONE

type(DIfile_T), INTENT(INOUT)   :: DIFT
real(kind=sgl), INTENT(IN)      :: gangle

type(IO_T)                      :: Message 
type(QuaternionArray_T)         :: qAR, QA, qsym
type(Quaternion_T)              :: muhat, quat
type(r_T)                       :: rod 
type(q_T)                       :: qu
type(e_T)                       :: eu
type(DirStat_T)                 :: dictVMF

integer(kind=irg)               :: nt, ix, iy, io_int(2), seed, icnt
integer(kind=irg),allocatable   :: grainIDs(:)
real(kind=dbl)                  :: kappahat

associate(nml=>DIFT%nml)

! rotation precision
call setRotationPrecision('Double')

cluster%gangle = gangle 

! dimensions of region of interest
if (sum(nml%ROI).eq.0) then 
  cluster%ipf_wd = nml%ipf_wd
  cluster%ipf_ht = nml%ipf_ht
else
  cluster%ipf_wd = nml%ROI(3)
  cluster%ipf_ht = nml%ROI(4)
end if
io_int = (/ cluster%ipf_wd, cluster%ipf_ht /)
call Message%WriteValue(' ROI dimensions : ', io_int, 2)
nt = cluster%ipf_wd * cluster%ipf_ht
allocate( cluster%grainID(cluster%ipf_wd, cluster%ipf_ht), &
          cluster%kam(cluster%ipf_wd, cluster%ipf_ht), grainIDs(nt) )

! next we need to compute the misorientations w.r.t. neighbors which
! is essentially the KAM map without thresholding... 
call Message%printMessage(' Computing misorientation map (KAM)')
call getKAMMap(nt, DIFT%DIDT%RefinedEulerAngles, cluster%ipf_wd, cluster%ipf_ht, DIFT%DIDT%pgnum, cluster%kam)
cluster%kam = cluster%kam*rtod

! to prevent weird edge cases, put the kam edges to a large value
cluster%kam(1,1:cluster%ipf_ht) = 360.0
cluster%kam(cluster%ipf_wd,1:cluster%ipf_ht) = 360.0
cluster%kam(1:cluster%ipf_wd,1) = 360.0
cluster%kam(1:cluster%ipf_wd,cluster%ipf_ht) = 360.0

! initialize number of grains and grainID array
cluster%nGrains = 0 
cluster%grainID = 0

! next, we recursively scan for each grain ... 
cluster%MaxSize = nt       ! used to make sure that the recursion will actually end
cluster%TotSize = 0

do ix=1,cluster%ipf_wd
  do iy=1,cluster%ipf_ht
    if (cluster%grainID(ix,iy).eq.0) then 
      if (cluster%kam(ix,iy).le.cluster%gangle) then
        cluster%nGrains = cluster%nGrains+1
        cluster%GrainSize = 0
        call cluster%ScanGrain_(ix, iy)
        if (cluster%GrainSize.eq.1) then 
          cluster%grainID(ix,iy) = -1
          cluster%nGrains = cluster%nGrains-1
        end if
      else
        cluster%grainID(ix,iy) = -1
      end if 
    end if 
  end do 
end do 

! for all grains, find the 2D bounding box needed for the modified DI algorithm
allocate( cluster%ROI(4,cluster%nGrains) )
call cluster%getROI_()

! from here on, it will be more useful to have the grain IDs in a 1D array
grainIDs = reshape(cluster%grainID, (/ nt /) )

! determine the number of pixels in each grain
allocate( cluster%npixels(cluster%nGrains) )
do ix = 1, cluster%nGrains
  cluster%npixels(ix) = count(cluster%grainID==ix)
end do

! next we need to compute the average orientation for each grain; this will 
! become the center of the misorientation ball used for the modified DI approach.
! avor will contain those orientations in quaternion form.
allocate( cluster%avor(4,cluster%nGrains), cluster%kappa(cluster%nGrains) )

io_int(1) = DIFT%DIDT%pgnum
call Message%WriteValue(' Initializing von Mises-Fisher distribution for point group # ',io_int,1)

! get the symmetry quaternions 
call QA%QSym_init( DIFT%DIDT%pgnum, qsym )

! initialize the von Mises-Fisher distribution code
dictVMF = DirStat_T( DStype='VMF', pgnum = DIFT%DIDT%pgnum)
call dictVMF%setNumEM(15)
call dictVMF%setNumIter(40)

! for each grain, set up the orientation array
seed = 32890

do ix = 1, cluster%nGrains
  icnt = 1
  qAR = QuaternionArray_T( n=cluster%npixels(ix), s='d' )
  do iy = 1, nt 
    if (grainIDs(iy).eq.ix) then 
      eu = e_T( edinp = dble(DIFT%DIDT%RefinedEulerAngles(1:3,iy)) )  
      qu = eu%eq()
      quat = Quaternion_T( qd = qu%q_copyd() )
      call qAR%insertQuatinArray(icnt,quat)
      icnt = icnt+1
    end if 
  end do

! pass these orientations to the dictVMF class  
  call dictVMF%setQuatArray( qAR )
! and perform the averaging step
  muhat = Quaternion_T( qd=(/ 1.D0, 0.D0, 0.D0,0.D0 /) )
  call dictVMF%EMforDS( seed, muhat, kappahat, verbose=.FALSE. )
  if (kappahat.gt.1000.D0) then   ! we only keep the orientations if Watson converged
! store muhat in the cluster%avor array 
    cluster%avor(1:4, ix) = muhat%get_quatd()
    cluster%kappa(ix) = kappahat
  else
    cluster%avor(1:4, ix) = (/ 1.D0, 0.D0, 0.D0, 0.D0 /)
    cluster%kappa(ix) = -1.D0
  end if

! and get rid of the orientation array
  call qAR%deleteArray()
end do

end associate

end function cluster_constructor

! !--------------------------------------------------------------------------
recursive subroutine ScanGrain_(self, is, js)
!DEC$ ATTRIBUTES DLLEXPORT :: ScanGrain_
!! author: MDG 
!! version: 1.0 
!! date: 05/22/25
!!

! this routine starts at the point (is, js), and determines
! whether or not it lies inside the grain with label nGrains.
! if it does, it will make sure that this particle is followed in 2D
!
! This is a recursive routine, so we must make sure that it uses 
! a minimal amount of memory...  To avoid a run-away recursive event,
! we'll count the number of on voxels and make sure it remains 
! smaller than MaxSize (which is the maximum possible number of on voxels).
! This should never happen, of course, but it was useful while 
! debugging the routine.
!
IMPLICIT NONE

class(Cluster_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN)   :: is, js   ! IPF coordinates and size

! have we reached the end of the recursion?  if yes, then force a stop
if (self%TotSize.lt.self%MaxSize) then 
! is this a point that we have not yet considered?
  if ((self%grainID(is,js).eq.0)) then 
! check the kam value; if it is smaller then the threshold, include this point
    if (self%kam(is,js).le.self%gangle) then 
      self%grainID(is,js) = self%nGrains          ! set this pixel to the current grain
      self%TotSize = self%TotSize + 1             ! increment the pixel-on counter
      self%GrainSize = self%GrainSize + 1
! then recurse to the neighbors
      if (is+1.le.self%ipf_wd) call self%ScanGrain_(is+1, js)   ! move to (x+1,y)
      if (is.gt.0) call self%ScanGrain_(is-1, js)               ! move to (x-1,y)
      if (js+1.le.self%ipf_ht) call self%ScanGrain_(is, js+1)   ! move to (x,y+1)
      if (js.gt.0) call self%ScanGrain_(is, js-1)               ! move to (x,y-1)
    end if
  end if
end if

! and return to the main program or the previous recursive routine call
end subroutine ScanGrain_

!--------------------------------------------------------------------------
recursive subroutine getROI_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: getROI_
!! author: MDG 
!! version: 1.0 
!! date: 05/21/25
!!
!! extract Region-of-Interest coordinates for all grain IDs

IMPLICIT NONE 

class(Cluster_T),INTENT(INOUT)    :: self

integer(kind=irg)                 :: i, ir, ic, xmin, xmax, ymin, ymax 

! simple brute force approach to finding each grain bounding box
do i=1,self%nGrains
  xmin = self%ipf_wd
  xmax = 1
  ymin = self%ipf_ht
  ymax = 1
  do ic=1,self%ipf_wd 
    do ir=1,self%ipf_ht
      if (self%grainID(ic,ir).eq.i) then 
        if (ic.lt.xmin) xmin = ic
        if (ic.gt.xmax) xmax = ic
        if (ir.lt.ymin) ymin = ir
        if (ir.gt.ymax) ymax = ir
      end if 
    end do 
  end do 
  self%ROI(1:4, i) = (/ xmin, ymin, xmax-xmin+1, ymax-ymin+1 /)
end do 

end subroutine getROI_

end module mod_cluster