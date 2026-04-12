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

integer(kind=irg), parameter :: max_stack_size = 1000000

! class definition
type, public :: Cluster_T
  integer(kind=irg), allocatable  :: grainID(:,:) ! grain identifier
  integer(kind=irg), allocatable  :: npixels(:)   ! number of pixels in each grain
  integer(kind=irg), allocatable  :: grainROI(:,:)! bounding box information for each grain
  integer(kind=irg)               :: ROI(4)       ! bounding box information for complete data set
  integer(kind=irg)               :: nGrains      ! number of grains found
  integer(kind=irg)               :: GrainSize    ! used to control recursion
  integer(kind=irg)               :: TotSize      ! used to control recursion
  integer(kind=irg)               :: MaxSize      ! used to control recursion
  integer(kind=irg)               :: ipf_wd       ! ROI width
  integer(kind=irg)               :: ipf_ht       ! ROI height
  real(kind=sgl)                  :: gangle       ! threshold angle for clustering
  real(kind=dbl), allocatable     :: avor(:,:)    ! average orientation per grain
  real(kind=dbl), allocatable     :: kappa(:)     ! concentration parameters
  real(kind=sgl), allocatable     :: kam(:,:)     ! kernel average orientation map
  integer(kind=irg), allocatable  :: x_stack(:) 
  integer(kind=irg), allocatable  :: y_stack(:)
contains
private 
  procedure, pass(self) :: ScanGrain_
  procedure, pass(self) :: getROI_
  procedure, pass(self) :: grain_dilate_
  procedure, pass(self) :: grow_region_
  procedure, pass(self) :: grow_region_driver_
end type Cluster_T

! the constructor routine for this class 
interface Cluster_T
  module procedure cluster_constructor
end interface Cluster_T

contains

!--------------------------------------------------------------------------
type(Cluster_T) function cluster_constructor( DIFT, gangle, dilate, orav, numEM, numIter, debug ) result(cluster)
!DEC$ ATTRIBUTES DLLEXPORT :: cluster_constructor
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
logical,INTENT(IN)              :: dilate
character(fnlen), INTENT(IN)    :: orav
integer(kind=irg), INTENT(IN)   :: numEM 
integer(kind=irg), INTENT(IN)   :: numIter
logical,INTENT(IN),OPTIONAL     :: debug

type(IO_T)                      :: Message 
type(QuaternionArray_T)         :: qAR, QA
type(Quaternion_T)              :: muhat, quat
type(r_T)                       :: rod 
type(q_T)                       :: qu
type(e_T)                       :: eu
type(DirStat_T)                 :: dictVMF

integer(kind=irg)               :: nt, ix, iy, i, j, k, io_int(2), seed, icnt, values(8)
integer(kind=irg),allocatable   :: grainIDs(:)
real(kind=dbl)                  :: kappahat
real(kind=sgl)                  :: ma

character(fnlen)                :: outname
character(3)                    :: filenum

associate(nml=>DIFT%nml)

! rotation precision
call setRotationPrecision('Double')

cluster%gangle = gangle 

! Note: there are two different ROI variables...
! - the original ROI (cluster%ROI) defined in the DI namelist file and read from the dp file
! - a per grain ROI (cluster%grainROI) that defines the bounding box for each individual grain
!
! the first one is used to determine the location of the patterns in the pattern file ... 

! dimensions of region of interest
cluster%ROI = nml%ROI
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
write (*,*) ' shape = ', shape(DIFT%DIDT%RefinedEulerAngles), maxval(DIFT%DIDT%RefinedEulerAngles)
call getKAMMap(nt, DIFT%DIDT%RefinedEulerAngles, cluster%ipf_wd, cluster%ipf_ht, DIFT%DIDT%pgnum, cluster%kam)
cluster%kam = cluster%kam*rtod

! ma = 1.2 * maxval( cluster%kam(2:cluster%ipf_wd-1,2:cluster%ipf_ht-1) )

! to prevent weird edge cases, put the kam edges to a large value
! cluster%kam(1,1:cluster%ipf_ht) = ma
! cluster%kam(cluster%ipf_wd,1:cluster%ipf_ht) = ma
! cluster%kam(1:cluster%ipf_wd,1) = ma
! cluster%kam(1:cluster%ipf_wd,cluster%ipf_ht) = ma

! find the grains
call cluster%grow_region_driver_()

! reset all the -1 values to 0
where(cluster%grainID==-1)
  cluster%grainID = 0
end where

! this is old code that performed a recursive grain search; for very large datasets 
! this caused potential memory problems, so we replaced it by the grow_region_driver
! routine ... [MDG, 06/25/25]

! initialize number of grains and grainID array
! cluster%nGrains = 0 
! cluster%grainID = 0


! ! next, we recursively scan for each grain ... 
! cluster%MaxSize = nt       ! used to make sure that the recursion will actually end
! cluster%TotSize = 0

! do ix=1,cluster%ipf_wd
!   do iy=1,cluster%ipf_ht
!     if (cluster%grainID(ix,iy).eq.0) then 
!       if (cluster%kam(ix,iy).le.cluster%gangle) then
!         cluster%nGrains = cluster%nGrains+1
!         cluster%GrainSize = 0
!         call cluster%ScanGrain_(ix, iy)
!         if (cluster%GrainSize.eq.1) then 
!           cluster%grainID(ix,iy) = -1
!           cluster%nGrains = cluster%nGrains-1
!         end if
!       else
!         cluster%grainID(ix,iy) = -1
!       end if 
!     end if 
!   end do 
! end do 

! for all grains, find the 2D bounding box needed for the modified DI algorithm
allocate( cluster%grainROI(4,cluster%nGrains) )
call cluster%getROI_()

! this next step is experimental at the moment...
! next, we optionally dilate all grains to remove most of the empty space at the grain boundaries
! and then we need to re-run the ROI finding routine to update the box sizes
if (dilate.eqv..TRUE.) then
  call cluster%grain_dilate_()
  call cluster%getROI_()
end if 

! determine the number of pixels in each grain
allocate( cluster%npixels(cluster%nGrains) )
do ix = 1, cluster%nGrains
  cluster%npixels(ix) = count(cluster%grainID==ix)
end do

! next we need to compute the average orientation for each grain; this will 
! become the center of the misorientation ball used for the modified DI approach.
! avor will contain those orientations in quaternion form.
! There are three options to get the orientation:
! - orientation of the center pixel of the bounding box
! - von Mises-Fisher average
! - Watson average
allocate( cluster%avor(4,cluster%nGrains), cluster%kappa(cluster%nGrains) )
if (trim(orav).eq.'center') then 
  open(dataunit,file='center.txt',status='unknown',form='formatted')
  write (dataunit,"(I5)") cluster%nGrains
  do i=1,cluster%nGrains
    ix = cluster%grainROI(1,i) + cluster%grainROI(3,i)/2
    iy = cluster%grainROI(2,i) + cluster%grainROI(4,i)/2
    j = (iy-1)*cluster%ipf_wd + ix 
    eu = e_T( edinp = dble(DIFT%DIDT%RefinedEulerAngles(1:3,j)))
    qu = eu%eq()
    cluster%avor(1:4,i) = qu%q_copyd()
    cluster%kappa(i) = 1.D0
    write (dataunit,"(4(F10.6,','),F14.6,',',I5)") qu%q_copyd(), cluster%kappa(i), cluster%npixels(i)
  end do 
  close(dataunit,status='keep')
else
  io_int(1) = DIFT%DIDT%pgnum

  ! initialize the von Mises-Fisher or Watson distribution code
  if (trim(orav).eq.'averageVMF') then 
    call Message%WriteValue(' Initializing von Mises-Fisher distribution for point group # ',io_int,1)
    dictVMF = DirStat_T( DStype='VMF', PGnum = DIFT%DIDT%pgnum)
    open(dataunit,file='VMF.txt',status='unknown',form='formatted')
    write (dataunit,"(I5)") cluster%nGrains
  else
    if (trim(orav).eq.'averageWAT') then 
      call Message%WriteValue(' Initializing Watson distribution for point group # ',io_int,1)
      dictVMF = DirStat_T( DStype='WAT', PGnum = DIFT%DIDT%pgnum)
      open(dataunit,file='WAT.txt',status='unknown',form='formatted')
      write (dataunit,"(I5)") cluster%nGrains
   else
      call Message%printError('cluster_constructor: ', 'unknown orientation averaging procedure') 
    end if
  end if 

  call dictVMF%setNumEM(numEM)
  call dictVMF%setNumIter(numIter)

  ! for each grain, set up the orientation array
  call date_and_time(values=values)
  seed = values(8)

  do j = 1, cluster%nGrains
    icnt = 1
    qAR = QuaternionArray_T( n=cluster%npixels(j), s='d' )
    do iy = 1, cluster%ipf_ht 
      do ix = 1, cluster%ipf_wd
        if (cluster%grainID(ix,iy).eq.j) then 
          k = (iy-1)*cluster%ipf_wd+ix
          eu = e_T( edinp = dble(DIFT%DIDT%RefinedEulerAngles(1:3,k)) )  
          qu = eu%eq()
          quat = Quaternion_T( qd = qu%q_copyd() )
          call qAR%insertQuatinArray(icnt,quat)
          icnt = icnt+1
        end if 
      end do
    end do

    if (present(debug)) then 
      if (debug.eqv..TRUE.) then 
        write (filenum,"(I3.3)") j
        outname = 'qu_grain_'//filenum//'.txt'
        call Message%printMessage(' writing orientations to '//trim(outname))
        call qAR%writeArraytoFile(outname)
      end if 
    end if 

  ! pass these orientations to the dictVMF class  
    call dictVMF%setQuatArray( qAR )
  ! and perform the averaging step
    muhat = Quaternion_T( qd=(/ 1.D0, 0.D0, 0.D0,0.D0 /) )
    call dictVMF%EMforDS( seed, muhat, kappahat, verbose=.TRUE. )
    if (kappahat.gt.5.D0) then   ! we only keep the orientations if Watson converged
  ! store muhat in the cluster%avor array 
      cluster%avor(1:4, j) = muhat%get_quatd()
      cluster%kappa(j) = kappahat
    else
      cluster%avor(1:4, j) = (/ 1.D0, 0.D0, 0.D0, 0.D0 /)
      cluster%kappa(j) = -1.D0
    end if

    qu = q_T( qdinp = cluster%avor(1:4,j) )
    write (dataunit,"(4(F10.6,','),F14.6,',',I5)") qu%q_copyd(), cluster%kappa(j), cluster%npixels(j)

  ! and get rid of the orientation array
    call qAR%deleteArray()
  end do
  close(dataunit,status='keep')
end if

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
  self%grainROI(1:4, i) = (/ xmin, ymin, xmax-xmin+1, ymax-ymin+1 /)
end do 

end subroutine getROI_

!--------------------------------------------------------------------------
recursive subroutine grain_dilate_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: grain_dilate_
!! author: MDG 
!! version: 1.0 
!! date: 06/12/25
!!
!! dilate all grains by one pixel in all directions to avoid empty spaces along 
!! the grain boundaries; this is a very naive implementation, probably not very 
!! efficient... ideally, we would use a convolution with a kernel to do this

IMPLICIT NONE 

class(Cluster_T),INTENT(INOUT)    :: self

integer                           :: i, j, m, sub(3,3)
integer(kind=irg),allocatable     :: im_in(:,:), im_out(:,:)


! allocate local arrays
allocate(im_in(0:self%ipf_wd+1, 0:self%ipf_ht+1), &
         im_out(0:self%ipf_wd+1, 0:self%ipf_ht+1))

!copy the current grainID into this array with a one-pixel zero border
im_in = 0
im_in(1:self%ipf_wd,1:self%ipf_ht) = self%grainID
im_out = im_in

! Scan each pixel (excluding borders)
do i = 1, self%ipf_wd-1
  do j = 1, self%ipf_ht-1
    sub = im_in(i:i+2, j:j+2)
    m = maxval(sub)
    if (m.ne.0) im_out(i+1,j+1) = m
  end do
end do

! copy the dilated array back into the grainID array
self%grainID = im_out(1:self%ipf_wd,1:self%ipf_ht)

end subroutine grain_dilate_

!--------------------------------------------------------------------------
subroutine grow_region_driver_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: grow_region_driver_
!! author: MDG 
!! version: 1.0 
!! date: 06/25/25
!!
!! based on Chat-GPT suggested algorithm

IMPLICIT NONE

class(Cluster_T), INTENT(INOUT)   :: self

integer(kind=irg)                 :: i, j

allocate(self%x_stack(max_stack_size), self%y_stack(max_stack_size) )

! initialize the grainID map and set the first label
self%grainID= 0
self%nGrains = 1

! Find and label all connected regions
do j = 1, self%ipf_ht
  do i = 1, self%ipf_wd
    if (self%grainID(i,j) == 0) then
      call self%grow_region_(i, j)
      ! don't accept grains that are only 1 pixel large
      if (count(self%grainID==self%nGrains).eq.1) then 
        self%grainID(i,j) = -1
      else
        self%nGrains = self%nGrains + 1
      end if 
    end if
  end do
end do

! this routine apparently overcounts the number of grains by 1, so we subtract 1
self%nGrains = self%nGrains - 1

deallocate(self%x_stack, self%y_stack)

end subroutine grow_region_driver_

!--------------------------------------------------------------------------
subroutine grow_region_(self, x_seed, y_seed)
!DEC$ ATTRIBUTES DLLEXPORT :: grow_region_
!! author: MDG 
!! version: 1.0 
!! date: 06/25/25
!!
!! based on Chat-GPT suggested algorithm

use mod_IO 

IMPLICIT NONE 

! Input
class(Cluster_T),INTENT(INOUT)    :: self
integer(kind=irg), INTENT(IN)     :: x_seed, y_seed

type(IO_T)                        :: Message 

! Local variables
integer(kind=irg)                 :: top, x, y, xn, yn, i

! 4-connected neighbor directions; can be extended to 8 neighbors if needed
! integer(kind=irg), parameter      :: dx(4) = [0, 0, -1, 1]
! integer(kind=irg), parameter      :: dy(4) = [-1, 1, 0, 0]
integer(kind=irg), parameter      :: dx(8) = [0, 0, -1, 1, 1, 1,-1,-1]
integer(kind=irg), parameter      :: dy(8) = [-1, 1, 0, 0, 1,-1, 1,-1]

! Initialize stack
self%x_stack = 0
self%y_stack = 0

! is this grain already labeled ?
if (self%grainID(x_seed, y_seed) /= 0) then 
  self%nGrains = self%nGrains-1 
  return  ! already labeled
end if

! is the kam value higher than the threshold? If so, then this is not a new grain
if (self%kam(x_seed, y_seed).gt.self%gangle) then 
  self%nGrains = self%nGrains-1 
  return  ! likely a grain boundary pixel
end if

top = 1
self%x_stack(top) = x_seed
self%y_stack(top) = y_seed
self%grainID(x_seed, y_seed) = self%nGrains

do while (top > 0)
  x = self%x_stack(top)
  y = self%y_stack(top)
  top = top - 1

  do i = 1, 8
    xn = x + dx(i)
    yn = y + dy(i)

    if (xn >= 1 .and. xn <= self%ipf_wd .and. yn >= 1 .and. yn <= self%ipf_ht) then
      if (self%grainID(xn, yn) == 0) then
        if (abs(self%kam(xn,yn) - self%kam(x,y)) <= self%gangle) then
          if (top < max_stack_size) then
            top = top + 1
            self%x_stack(top) = xn
            self%y_stack(top) = yn
            self%grainID(xn, yn) = self%nGrains
          else
            call Message%printMessage(' grow_region Error: Stack overflow! ... returning ')
            return
          end if
        end if
      end if
    end if
  end do
end do
end subroutine grow_region_




end module mod_cluster
