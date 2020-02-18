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

module mod_so3
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/21/20
  !!
  !! everything that has to do with sampling of rotation space SO(3)

use mod_kinds
use mod_global
use mod_rotations

IMPLICIT NONE
private 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! the following table is used for two-phase disorientation fundamental zones
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! this table encodes Figure 1 of the paper  "Representation of Orientation and
! Disorientation data for Cubic, Hexagonal, Tetragonal, and Orthorhombic Crystals", A. Heinz
! and P. Neumann, Acta Cryst. A47, 780-789 (1991)
! The following conversions are used
! 0 -> x  (no symmetry)
! 1 -> a  mixed cubic-hexagonal FZ
! 2 -> b  mixed FZ
! 3 -> c  octahedral FZ
! 4 -> d  tetrahedral FZ
! 5 -> e  24-sided prismatic FZ
! 6 -> f  622 hexagonal dihedral FZ
! 7 -> g  422 octagonal dihedral FZ
! 8 -> h  32 trigonal dihedral FZ
! 9 -> i  222 dihedral FZ
! This table is used in the so3.f90 module to figure out which FZ should be used for a single phase
! or two phase FZ computation; all FZs are also available in the povray.f90 module for 3D visualization.
! The new routine getFZtypeandorder in so3.f90 will take two point group numbers, possibly identical,
! and return the class FZtype and FZorder parameters that are currently used already in other routines.  
integer(kind=irg), dimension(32,32) :: FZtypeTable = reshape( (/ &
 0, 0, 0, 0, 0, 9, 0, 9, 0, 0, 0, 7, 0, 9, 7, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 8, 6, 4, 4, 3, 4, 3, &
 0, 0, 0, 0, 0, 9, 0, 9, 0, 0, 0, 7, 0, 9, 7, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 8, 6, 4, 4, 3, 4, 3, &
 0, 0, 0, 0, 0, 9, 9, 9, 7, 9, 7, 7, 7, 7, 7, 8, 8, 6, 8, 6, 6, 8, 6, 6, 6, 6, 6, 4, 4, 3, 4, 3, &
 0, 0, 0, 0, 0, 9, 0, 9, 0, 0, 0, 7, 0, 9, 7, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 8, 6, 4, 4, 3, 4, 3, &
 0, 0, 0, 0, 0, 9, 9, 9, 7, 9, 7, 7, 7, 7, 7, 8, 8, 6, 8, 6, 6, 8, 6, 6, 6, 6, 6, 4, 4, 3, 4, 3, &
 9, 9, 9, 9, 9, 9, 9, 9, 7, 9, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 4, 4, 3, 4, 3, &
 0, 0, 9, 0, 9, 9, 0, 9, 0, 0, 0, 7, 0, 9, 7, 0, 0, 6, 0, 6, 0, 0, 0, 6, 0, 6, 6, 4, 4, 3, 4, 3, &
 9, 9, 9, 9, 9, 9, 9, 9, 7, 9, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 4, 4, 3, 4, 3, &
 0, 0, 7, 0, 7, 7, 0, 7, 0, 0, 0, 7, 0, 7, 7, 0, 0, 5, 0, 5, 0, 0, 0, 5, 0, 5, 5, 3, 3, 3, 3, 3, &
 0, 0, 9, 0, 9, 9, 0, 9, 0, 0, 0, 7, 0, 9, 7, 0, 0, 6, 0, 6, 0, 0, 0, 6, 0, 6, 6, 4, 4, 3, 4, 3, &
 0, 0, 7, 0, 7, 7, 0, 7, 0, 0, 0, 7, 0, 7, 7, 0, 0, 5, 0, 5, 0, 0, 0, 5, 0, 5, 5, 3, 3, 3, 3, 3, &
 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 3, 3, 3, 3, &
 0, 0, 7, 0, 7, 7, 0, 7, 0, 0, 0, 7, 0, 7, 7, 0, 0, 5, 0, 5, 0, 0, 0, 5, 0, 5, 5, 3, 3, 3, 3, 3, &
 9, 9, 7, 9, 7, 7, 9, 7, 7, 9, 7, 7, 7, 9, 7, 6, 6, 5, 6, 5, 6, 6, 6, 5, 6, 5, 5, 3, 3, 3, 3, 3, &
 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 3, 3, 3, 3, &
 0, 0, 8, 0, 8, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 8, 6, 2, 2, 1, 2, 1, &
 0, 0, 8, 0, 8, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 8, 6, 2, 2, 1, 2, 1, &
 8, 8, 6, 8, 6, 6, 6, 6, 5, 6, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 6, 8, 6, 6, 6, 8, 6, 2, 2, 1, 2, 1, &
 0, 0, 8, 0, 8, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 8, 6, 2, 2, 1, 2, 1, &
 8, 8, 6, 8, 6, 6, 6, 6, 5, 6, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 6, 8, 6, 6, 6, 8, 6, 2, 2, 1, 2, 1, &
 0, 0, 6, 0, 6, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 6, 0, 6, 0, 0, 0, 6, 0, 6, 6, 2, 2, 1, 2, 1, &
 0, 0, 8, 0, 8, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 6, 6, 2, 2, 1, 2, 1, &
 0, 0, 6, 0, 6, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 6, 0, 6, 0, 0, 0, 6, 0, 6, 6, 2, 2, 1, 2, 1, &
 6, 6, 6, 6, 6, 6, 6, 6, 5, 6, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 1, 2, 1, &
 0, 0, 6, 0, 6, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 6, 0, 6, 0, 0, 0, 6, 0, 6, 6, 2, 2, 1, 2, 1, &
 8, 8, 6, 8, 6, 6, 6, 6, 5, 6, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 6, 6, 6, 6, 6, 8, 6, 2, 2, 1, 2, 1, &
 6, 6, 6, 6, 6, 6, 6, 6, 5, 6, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 1, 2, 1, &
 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 3, 4, 3, &
 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 3, 4, 3, &
 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, &
 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 3, 4, 3, &
 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3 &
 /), (/ 32, 32/) )

! The following two arrays are used to determine the FZtype (FZtarray) and primary rotation axis order (FZoarray)
! for each of the 32 crystallographic point group symmetries (in the order of the International Tables)
!
!                                       '    1','   -1','    2','    m','  2/m','  222', &
!                                       '  mm2','  mmm','    4','   -4','  4/m','  422', &
!                                       '  4mm',' -42m','4/mmm','    3','   -3','   32', &
!                                       '   3m','  -3m','    6','   -6','  6/m','  622', &
!                                       '  6mm',' -6m2','6/mmm','   23','   m3','  432', &
!                                       ' -43m',' m-3m'/
!
! 1 (C1), -1 (Ci), [triclinic]
! 2 (C2), m (Cs), 2/m (C2h), [monoclinic]
! 222 (D2), mm2 (C2v), mmm (D2h), [orthorhombic]
! 4 (C4), -4 (S4), 4/m (C4h), 422 (D4), 4mm (C4v), -42m (D2d), 4/mmm (D4h), [tetragonal]
! 3 (C3), -3 (C3i), 32 (D3), 3m (C3v), -3m (D3d), [trigonal]
! 6 (C6), -6 (C3h), 6/m (C6h), 622 (D6), 6mm (C6v), -6m2 (D3h), 6/mmm (D6h), [hexagonal]
! 23 (T), m3 (Th), 432 (O), -43m (Td), m-3m (Oh) [cubic]
!
! FZtype
! 0        no symmetry at all
! 1        cyclic symmetry
! 2        dihedral symmetry
! 3        tetrahedral symmetry
! 4        octahedral symmetry
!
integer(kind=irg),dimension(36)     :: FZtarray = (/ 0,0,1,1,1,2,2,2,1,1,1,2,2,2,2,1,1,2, &
                                                     2,2,1,1,1,2,2,2,2,3,3,4,3,4,5,2,2,2 /)

integer(kind=irg),dimension(36)     :: FZoarray = (/ 0,0,2,2,2,2,2,2,4,4,4,4,4,4,4,3,3,3, &
                                                     3,3,6,6,6,6,6,6,6,0,0,0,0,0,0,8,10,12 /)



! public :: SampleRFZ, IsinsideFZ, CubochoricNeighbors

! logical functions to determine if point is inside specific FZ
!private :: insideCyclicFZ, insideDihedralFZ, insideCubicFZ
! public:: insideCyclicFZ, insideDihedralFZ, insideCubicFZ

type, public :: FZpointd
  type(r_T)               :: rod       ! Rodrigues-Frank vector [nx, ny, nz, tan(omega/2) ]
  type(r_T)               :: trod      ! second Rodrigues-Frank vector; can be used for coordinate transformations
  integer(kind=irg)       :: gridpt(3) ! coordinates of grid point ! added on 06/19/18 by SS
  type(FZpointd),pointer  :: next      ! link to next point
end type FZpointd


type, public :: so3_T
  private
    integer(kind=irg)       :: FZtype 
    integer(kind=irg)       :: FZ2type 
    integer(kind=irg)       :: FZorder
    integer(kind=irg)       :: MFZtype 
    integer(kind=irg)       :: MFZorder
    integer(kind=irg)       :: pgnum
    integer(kind=irg)       :: pgnum2
    integer(kind=irg)       :: gridtype
    integer(kind=irg)       :: FZcnt
    integer(kind=irg)       :: CMcnt
    integer(kind=irg)       :: COcnt
    integer(kind=irg)       :: FBcnt
    type(FZpointd),pointer  :: FZlist 
    type(FZpointd),pointer  :: CMlist  ! CM = Constant Misorientation
    type(FZpointd),pointer  :: COlist  ! CO = Cone sampling
    type(FZpointd),pointer  :: FBlist  ! FB = Fiber sampling
  contains
  private 

    procedure, pass(self) :: getFZtypeandorder_ 
    procedure, pass(self) :: setFZtypeandorder_ 
    procedure, pass(self) :: getMFZtypeandorder_ 
    procedure, pass(self) :: setMFZtypeandorder_ 
    procedure, pass(self) :: IsinsideFZ_
    procedure, pass(self) :: IsinsideMFZ_
    procedure, pass(self) :: insideCubicMFZ_
    procedure, pass(self) :: insideDihedralMFZ_
    procedure, pass(self) :: insideIcosahedralFZ_
    procedure, pass(self) :: insideCyclicFZ_
    procedure, pass(self) :: insideDihedralFZ_
    procedure, pass(self) :: insideCubicFZ_
    procedure, pass(self) :: insideCubeHexFZ_
    procedure, pass(self) :: listtoArray_
    procedure, pass(self) :: listtoQuaternionArray_
    procedure, pass(self) :: QuaternionArraytolist_
    procedure, pass(self) :: getListHead_
    procedure, pass(self) :: getListCount_
    procedure, pass(self) :: setGridType_

    procedure, pass(self) :: delete_FZlist_
    procedure, pass(self) :: nullifyList_
    procedure, pass(self) :: SampleRFZ_
    procedure, pass(self) :: CubochoricNeighbors_
    procedure, pass(self) :: sample_isoCube_
    procedure, pass(self) :: sample_isoCubeFilled_
    procedure, pass(self) :: sample_Cone_
    procedure, pass(self) :: sample_Fiber_
    procedure, pass(self) :: SampleIsoMisorientation_
    procedure, pass(self) :: getOrientationsfromFile_
    procedure, pass(self) :: writeOrientationstoFile_
    procedure, pass(self) :: getVertex_
    procedure, pass(self) :: getMacKenzieDistribution_
! some other related routines
    procedure, pass(self) :: ReducelisttoRFZ_
    procedure, pass(self) :: ReduceDisorientationtoMFZ_
    procedure, pass(self) :: ReduceOrientationtoCubicEFZ_
    procedure, pass(self) :: ReduceOrientationtoRFZ_
    procedure, pass(self) :: getDisorientation_
    procedure, pass(self) :: getDisorientationTwoPhases_
    procedure, pass(self) :: getAverageDisorientationMap_
    final :: so3_destructor

    generic, public :: getFZtypeandorder => getFZtypeandorder_
    generic, public :: setFZtypeandorder => setFZtypeandorder_
    generic, public :: getMFZtypeandorder => getMFZtypeandorder_
    generic, public :: setMFZtypeandorder => setMFZtypeandorder_
    generic, public :: IsinsideFZ => IsinsideFZ_
    generic, public :: IsinsideMFZ => IsinsideMFZ_
    generic, public :: insideCubicMFZ => insideCubicMFZ_
    generic, public :: insideDihedralMFZ => insideDihedralMFZ_
    generic, public :: insideIcosahedralFZ => insideIcosahedralFZ_
    generic, public :: insideCyclicFZ => insideCyclicFZ_
    generic, public :: insideDihedralFZ => insideDihedralFZ_
    generic, public :: insideCubicFZ => insideCubicFZ_
    generic, public :: insideCubeHexFZ => insideCubeHexFZ_
    generic, public :: listtoArray => listtoArray_
    generic, public :: listtoQuaternionArray => listtoQuaternionArray_
    generic, public :: QuaternionArraytolist => QuaternionArraytolist_
    generic, public :: getListHead => getListHead_
    generic, public :: getListCount => getListCount_
    generic, public :: setGridType => setGridType_

    generic, public :: delete_FZlist => delete_FZlist_
    generic, public :: nullifyList => nullifyList_
    generic, public :: SampleRFZ => SampleRFZ_
    generic, public :: CubochoricNeighbors => CubochoricNeighbors_
    generic, public :: sample_isoCube => sample_isoCube_
    generic, public :: sample_isoCubeFilled => sample_isoCubeFilled_
    generic, public :: sample_Cone => sample_Cone_
    generic, public :: sample_Fiber => sample_Fiber_
    generic, public :: SampleIsoMisorientation => SampleIsoMisorientation_
    generic, public :: getOrientationsfromFile => getOrientationsfromFile_
    generic, public :: writeOrientationstoFile => writeOrientationstoFile_
    generic, public :: getVertex => getVertex_
    generic, public :: getMacKenzieDistribution => getMacKenzieDistribution_
    generic, public :: ReducelisttoRFZ => ReducelisttoRFZ_

    generic, public :: ReduceDisorientationtoMFZ => ReduceDisorientationtoMFZ_
    generic, public :: ReduceOrientationtoCubicEFZ => ReduceOrientationtoCubicEFZ_
    generic, public :: ReduceOrientationtoRFZ => ReduceOrientationtoRFZ_
    generic, public :: getDisorientation => getDisorientation_, getDisorientationTwoPhases_
    generic, public :: getAverageDisorientationMap => getAverageDisorientationMap_

end type so3_T

! the constructor routine for this class 
interface so3_T
  module procedure so3_constructor
end interface so3_T

contains

!--------------------------------------------------------------------------
type(so3_T) function so3_constructor( pgnum, pgnum2, zerolist ) result(SO)
!! author: MDG 
!! version: 1.0 
!! date: 01/21/20
!!
!! constructor for the so3_T Class 
 
IMPLICIT NONE

integer(kind=irg), INTENT(IN)             :: pgnum 
 !! primary point group
integer(kind=irg), INTENT(IN), OPTIONAL   :: pgnum2 
 !! optional secondary point group
character(2), INTENT(IN), OPTIONAL        :: zerolist
 !! optional selector for linked list to be reset

if (present(pgnum2)) then 
  call SO%setFZtypeandorder(pgnum, pgnum2)
  SO%pgnum2 = pgnum2
else 
  call SO%setFZtypeandorder(pgnum)
end if 
SO%pgnum = pgnum

if (present(zerolist)) then 
  call SO%nullifyList(zerolist)
else 
  call SO%nullifyList()
end if 

end function so3_constructor

!--------------------------------------------------------------------------
subroutine so3_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 02/02/20
!!
!! destructor for the so3_T Class
 
IMPLICIT NONE

type(so3_T), INTENT(INOUT)  :: self 

call reportDestructor('so3_T')

! nothing to do for now...

end subroutine so3_destructor

!--------------------------------------------------------------------------
subroutine nullifyList_(self, zerolist) 
!! author: MDG 
!! version: 1.0 
!! date: 02/18/20
!!
!! nullify the selected list 
 
IMPLICIT NONE

class(so3_T), INTENT(INOUT)           :: self 
character(2), INTENT(IN), OPTIONAL    :: zerolist
 !! optional selector for linked list to be reset

if (present(zerolist)) then 
  select case(zerolist)
  case('FZ')
    nullify(self%FZlist)
    self%FZcnt = 0
  case('CM')
    nullify(self%CMlist)
    self%CMcnt = 0
  case('CO')
    nullify(self%COlist)
    self%COcnt = 0
  case('FB')
    nullify(self%FBlist)
    self%FBcnt = 0
  case default 
    nullify(self%FZlist)
    self%FZcnt = 0
  end select
else
  nullify(self%FZlist)
  nullify(self%CMlist)
  nullify(self%COlist)
  nullify(self%FBlist)
  self%FZcnt = 0
  self%CMcnt = 0
  self%COcnt = 0
  self%FBcnt = 0
end if 

end subroutine nullifyList_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Routine to return the FZtype and FZorder parameters for single or two-phase
! fundamental zone (FZ) computations; this includes all the FZ types from the 
! following paper:
!
! "Representation of Orientation and Disorientation data for Cubic, Hexagonal, 
! Tetragonal, and Orthorhombic Crystals", A. Heinz and P. Neumann, Acta Cryst. A47, 
! 780-789 (1991)
!
! this routine also allows for icosahedral symmetry, although this is not part 
! of the paper above.
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine setFZtypeandorder_(self, pgnum1, pgnum2) 
!! author: MDG 
!! version: 1.0 
!! date: 01/21/20
!!
!! set the point group number(s) and the fundamental zone type and order  

IMPLICIT NONE

class(so3_T),INTENT(INOUT)                :: self 
integer(kind=irg),INTENT(IN)              :: pgnum1
 !! point group number for first point group
integer(kind=irg),INTENT(IN),OPTIONAL     :: pgnum2
 !! point group number for optional second point group

integer(kind=irg)                         :: thisFZType
logical                                   :: twophase

twophase = .FALSE.
if (present(pgnum2)) twophase = .TRUE.

if (twophase.eqv..TRUE.) then
  thisFZtype = FZtypeTable(pgnum1,pgnum2)

  select case(thisFZtype)
    case (0)
      ! this needs some more work since we need to properly handle the cyclic groups of FZtype = 1 ...
      self%FZtype = 0
      self%FZorder = 0
    case (1)
      self%FZtype = 6
      self%FZorder = 0
    case (2)
      self%FZtype = 7
      self%FZorder = 0
    case (3)
      self%FZtype = 4
      self%FZorder = 0
    case (4)
      self%FZtype = 3
      self%FZorder = 0
    case (5)
      self%FZtype = 8
      self%FZorder = 0
    case (6)
      self%FZtype = 2
      self%FZorder = 6
    case (7)
      self%FZtype = 2
      self%FZorder = 4
    case (8)
      self%FZtype = 2
      self%FZorder = 3
    case (9)
      self%FZtype = 2
      self%FZorder = 2
  end select
else  ! single phase so use the old way of doing things...
  self%FZtype = FZtarray(pgnum1)
  self%FZorder = FZoarray(pgnum1)
end if

end subroutine setFZtypeandorder_

!--------------------------------------------------------------------------
recursive subroutine getFZtypeandorder_(self, FZtype, FZorder) 
!! author: MDG 
!! version: 1.0 
!! date: 01/21/20
!!
!! set the point group number(s) and the fundamental zone type and order  

IMPLICIT NONE

class(so3_T),INTENT(INOUT)                :: self 
integer(kind=irg), INTENT(OUT)            :: FZtype 
integer(kind=irg), INTENT(OUT)            :: FZorder 

FZtype = self%FZtype 
FZorder = self%FZorder 

end subroutine getFZtypeandorder_

!--------------------------------------------------------------------------
recursive subroutine setMFZtypeandorder_(self, pgnum1, pgnum2) 
!! author: MDG 
!! version: 1.0 
!! date: 01/21/20
!!
!! set the point group number(s) and the Mackenzie fundamental zone type and order  

IMPLICIT NONE

class(so3_T),INTENT(INOUT)                :: self 
integer(kind=irg),INTENT(IN)              :: pgnum1
 !! point group number for first point group
integer(kind=irg),INTENT(IN),OPTIONAL     :: pgnum2
 !! point group number for optional second point group

integer(kind=irg)                         :: thisFZType
logical                                   :: twophase

twophase = .FALSE.
if (present(pgnum2)) twophase = .TRUE.

if (twophase.eqv..TRUE.) then
  thisFZtype = FZtypeTable(pgnum1,pgnum2)

  select case(thisFZtype)
    case (0)
      ! this needs some more work since we need to properly handle the cyclic groups of FZtype = 1 ...
      self%MFZtype = 0
      self%MFZorder = 0
    case (1)
      self%MFZtype = 6
      self%MFZorder = 0
    case (2)
      self%MFZtype = 7
      self%MFZorder = 0
    case (3)
      self%MFZtype = 4
      self%MFZorder = 0
    case (4)
      self%MFZtype = 3
      self%MFZorder = 0
    case (5)
      self%MFZtype = 8
      self%MFZorder = 0
    case (6)
      self%MFZtype = 2
      self%MFZorder = 6
    case (7)
      self%MFZtype = 2
      self%MFZorder = 4
    case (8)
      self%MFZtype = 2
      self%MFZorder = 3
    case (9)
      self%MFZtype = 2
      self%MFZorder = 2
  end select
else  ! single phase so use the old way of doing things...
  self%MFZtype = FZtarray(pgnum1)
  self%MFZorder = FZoarray(pgnum1)
end if

end subroutine setMFZtypeandorder_

!--------------------------------------------------------------------------
recursive subroutine getMFZtypeandorder_(self, MFZtype, MFZorder) 
!! author: MDG 
!! version: 1.0 
!! date: 01/21/20
!!
!! set the point group number(s) and the Mackenzie fundamental zone type and order  

IMPLICIT NONE

class(so3_T),INTENT(INOUT)        :: self 
integer(kind=irg), INTENT(OUT)    :: MFZtype 
integer(kind=irg), INTENT(OUT)    :: MFZorder 

MFZtype = self%MFZtype 
MFZorder = self%MFZorder 

end subroutine getMFZtypeandorder_


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! We define a number of logical routines, that decide whether or not 
! a point in Rodrigues representation lies inside the fundamental zone (FZ)
! for a given crystal symmetry. This follows the Morawiec@Field paper:
!
! A. Morawiec & D. P. Field (1996) Rodrigues parameterization for orientation 
! and misorientation distributions, Philosophical Magazine A, 73:4, 1113-1130, 
! DOI: 10.1080/01418619608243708
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive function IsinsideFZ_(self, rod, qFZ) result(insideFZ)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/21/20
  !!
  !! does Rodrigues point lie inside the relevant FZ?

use mod_math
use mod_quaternions

IMPLICIT NONE

class(so3_T),INTENT(INOUT)              :: self 
type(r_T), INTENT(INOUT)                :: rod
 !! input Rodrigues vector 
type(q_T), INTENT(INOUT), OPTIONAL      :: qFZ
 !! quaternion that rotates the FZ into a new orientation (optional)

logical                                 :: insideFZ

type(r_T)                               :: newrod
type(q_T)                               :: qu
type(quaternion_T)                      :: qu1, qu2, qq
real(kind=dbl)                          :: x(4)

! do we need to rotate the FZ ? (we do this by rotating rod in the opposite way)
if (present(qFZ)) then 
  qu = rod%rq()
  qu1 = quaternion_T( qd = qu%q_copyd() )
  qu2 = quaternion_T( qd = qFZ%q_copyd() )
  qq = qu2 * ( qu1 * conjg(qu2) )
  qu = q_T( qdinp = qq%get_quatd() )
  newrod = qu%qr()
else
  newrod = rod
end if

insideFZ = .FALSE.

! dealing with 180 rotations is needed only for 
! FZtypes 0 and 1; the other FZs are always finite.
x = newrod%r_copyd()
  select case (self%FZtype)
    case (0)
      insideFZ = .TRUE.   ! all points are inside the FZ
    case (1)
      insideFZ = self%insideCyclicFZ(newrod)        ! infinity is checked inside this function
    case (2)
      if (x(4).ne.inftyd()) insideFZ = self%insideDihedralFZ(newrod, self%FZorder)
    case (3)
      if (x(4).ne.inftyd()) insideFZ = self%insideCubicFZ(newrod,'tet')
    case (4)
      if (x(4).ne.inftyd()) insideFZ = self%insideCubicFZ(newrod,'oct')
    case (5) ! icosahedral symmetry
      if (x(4).ne.inftyd()) insideFZ = self%insideIcosahedralFZ(newrod)
    case (6) ! cubic-hexagonal misorientation FZ
      if (x(4).ne.inftyd()) insideFZ = self%insideCubeHexFZ(newrod)
    case (7)
!     if (x(4).ne.inftyd) insideFZ = self%insideCubicFZ(rod,'oct')
    case (8)
!     if (x(4).ne.inftyd) insideFZ = self%insideCubicFZ(rod,'oct')
  end select

end function IsinsideFZ_

!--------------------------------------------------------------------------
recursive function insideIcosahedralFZ_(self, rod) result(res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/21/20
  !!
  !! does Rodrigues point lie inside icosahedral FZ?

IMPLICIT NONE 

class(so3_T),INTENT(INOUT)        :: self 
type(r_T), INTENT(INOUT)          :: rod

logical                           :: res

real(kind=dbl)                    :: dval, rv(3), r(4)
integer(kind=irg)                 :: i, j

res = .FALSE.

dval=0.32491969623290632616D0  ! sqrt(1-2/sqrt(5)))
r = rod%r_copyd()
rv(1:3) = r(1:3)*r(4)
j = 0
do i=1,12
  if (DOT_PRODUCT(IcoVertices(1:3,i),rv)+dval.ge.0.D0) j = j+1
end do
if (j.eq.12) res = .TRUE.

end function insideIcosahedralFZ_

!--------------------------------------------------------------------------
recursive function insideCyclicFZ_(self, rod, M) result(res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/21/20
  !!
  !! does Rodrigues point lie inside cyclic FZ (for 2, 3, 4, and 6-fold)?

use mod_math

IMPLICIT NONE 

class(so3_T),INTENT(INOUT)        :: self 
type(r_T), INTENT(INOUT)          :: rod
logical,INTENT(IN),OPTIONAL       :: M

logical                           :: res, doM
real(kind=dbl)                    :: x(4)

res = .FALSE.
doM = .FALSE.
if (present(M)) then 
  if (M.eqv..TRUE.) doM = .TRUE.
end if 

x = rod%r_copyd()
if (x(4).ne.inftyd()) then
  if (doM.eqv..TRUE.) then
    if ((self%MFZtype.eq.1.).and.(self%MFZorder.eq.2)) then
! check the y-component vs. tan(pi/2n)
      res = dabs(x(2)*x(4)).le.LPs%BP(self%MFZorder)
    else
! check the z-component vs. tan(pi/2n)
      res = dabs(x(3)*x(4)).le.LPs%BP(self%MFZorder)
    end if
  else 
    if ((self%FZtype.eq.1.).and.(self%FZorder.eq.2)) then
! check the y-component vs. tan(pi/2n)
      res = dabs(x(2)*x(4)).le.LPs%BP(self%FZorder)
    else
! check the z-component vs. tan(pi/2n)
      res = dabs(x(3)*x(4)).le.LPs%BP(self%FZorder)
    end if
  end if 
else
  if (doM.eqv..TRUE.) then
    if ((self%MFZtype.eq.1.).and.(self%MFZorder.eq.2)) then
      if(x(2) .eq. 0.D0) res = .TRUE.
    else
      if (x(3).eq.0.D0) res = .TRUE.
    end if
  else
    if ((self%FZtype.eq.1.).and.(self%FZorder.eq.2)) then
      if(x(2) .eq. 0.D0) res = .TRUE.
    else
      if (x(3).eq.0.D0) res = .TRUE.
    end if
  end if 
endif

end function insideCyclicFZ_

!--------------------------------------------------------------------------
recursive function insideDihedralFZ_(self, rod, order) result(res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/21/20
  !!
  !! does Rodrigues point lie inside cyclic FZ (for 2, 3, 4, and 6-fold)?

IMPLICIT NONE 

class(so3_T),INTENT(INOUT)        :: self 
type(r_T), INTENT(INOUT)          :: rod
integer(kind=irg), INTENT(IN)     :: order

logical                           :: res, c1, c2
real(kind=dbl)                    :: r(3), x(4)
real(kind=dbl),parameter          :: r1 = 1.00D0
real(kind=dbl),allocatable        :: polygonvertex(:,:)
integer(kind=irg)                 :: inout

x = rod%r_copyd()
if (x(4).gt.sqrt(3.D0)) then 
  res = .FALSE.
else
  r(1:3) = x(1:3) * x(4)

  ! first, check the z-component vs. tan(pi/2n)  (same as insideCyclicFZ)
  c1 = dabs(r(3)).le.LPs%BP(order)
  res = .FALSE.

  ! check the square boundary planes if c1=.TRUE.
  if (c1) then
    select case (order)
      case (2)
        c2 = (dabs(r(1)).le.r1).and.(dabs(r(2)).le.r1)
      case (3)
        c2 =          dabs( LPs%srt*r(1)+0.5D0*r(2)).le.r1
        c2 = c2.and.( dabs( LPs%srt*r(1)-0.5D0*r(2)).le.r1 )
        c2 = c2.and.( dabs(r(2)).le.r1 )
      case (4)
        c2 = (dabs(r(1)).le.r1).and.(dabs(r(2)).le.r1)
        c2 = c2.and.((LPs%r22*dabs(r(1)+r(2)).le.r1).and.(LPs%r22*dabs(r(1)-r(2)).le.r1))
      case (6)
        c2 =          dabs( 0.5D0*r(1)+LPs%srt*r(2)).le.r1
        c2 = c2.and.( dabs( LPs%srt*r(1)+0.5D0*r(2)).le.r1 )
        c2 = c2.and.( dabs( LPs%srt*r(1)-0.5D0*r(2)).le.r1 )
        c2 = c2.and.( dabs( 0.5D0*r(1)-LPs%srt*r(2)).le.r1 )
        c2 = c2.and.( dabs(r(2)).le.r1 )
        c2 = c2.and.( dabs(r(1)).le.r1 )

      ! add the 2-D quasi crystal type for 822, 1022, and 1222 rotational groups
      case (8)
        c2 = .FALSE.
        allocate(polygonvertex(order*2, 2))
        polygonvertex = 0.D0
        call self%getVertex(order, polygonvertex)
        inout = PNPOLY(r(1),r(2),polygonvertex(1:2*order,1),polygonvertex(1:2*order,2),order*2)
        if(inout .ge. 0) c2 = .TRUE.

      case (10)
        c2 = .FALSE.
        allocate(polygonvertex(order*2, 2))
        polygonvertex = 0.D0
        call self%getVertex(order, polygonvertex)
        inout = PNPOLY(r(1),r(2),polygonvertex(1:2*order,1),polygonvertex(1:2*order,2),order*2)
        if(inout .ge. 0) c2 = .TRUE.

      case(12)
        c2 = .FALSE.
        allocate(polygonvertex(order*2, 2))
        polygonvertex = 0.D0
        call self%getVertex(order, polygonvertex)
        inout = PNPOLY(r(1),r(2),polygonvertex(1:2*order,1),polygonvertex(1:2*order,2),order*2)
        if(inout .ge. 0) c2 = .TRUE.

    end select
    res = c2
  end if
end if

end function insideDihedralFZ_

!--------------------------------------------------------------------------
recursive function insideCubicFZ_(self, rod, ot) result(res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/21/20
  !!
  !! does Rodrigues point lie inside cubic FZ (octahedral or tetrahedral)?

IMPLICIT NONE 

class(so3_T),INTENT(INOUT)        :: self 
type(r_T), INTENT(INOUT)          :: rod
character(3), INTENT(IN)          :: ot

logical                           :: res, c1, c2
real(kind=dbl)                    :: r(3), x(4)
real(kind=dbl),parameter          :: r1  = 1.0D0
real(kind=dbl),parameter          :: eps = 1.0D-8

x = rod%r_copyd()
r(1:3) = x(1:3) * x(4)

res = .FALSE.

! primary cube planes (only needed for octahedral case)
if (ot.eq.'oct') then
  c1 = (maxval(dabs(r)) - LPS%BP(4) .le. eps) 
else 
  c1 = .TRUE.
end if

! octahedral truncation planes, both for tetrahedral and octahedral point groups
c2 = ((dabs(r(1))+dabs(r(2))+dabs(r(3))) - r1 .le. eps)

! if both c1 and c2, then the point is inside
if (c1.and.c2) res = .TRUE.

end function insideCubicFZ_

!--------------------------------------------------------------------------
recursive function insideCubeHexFZ_(self, rod) result(res)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/21/20
  !!
  !! does Rodrigues point lie inside combined cubic-hexagonal FZ?
  !!
  !! For details on this test, see section 8 in "Representation of Orientation and
  !! Disorientation data for Cubic, Hexagonal, Tetragonal, and Orthorhombic Crystals", A. Heinz
  !! and P. Neumann, Acta Cryst. A47, 780-789 (1991)

IMPLICIT NONE 

class(so3_T),INTENT(INOUT)        :: self 
type(r_T), INTENT(INOUT)          :: rod

logical                           :: res
real(kind=dbl)                    :: r(3), x(4)
real(kind=dbl),parameter          :: r1 = 0.414213562373095D0, r2 = 0.131652497587396D0, &
                                     alpha = 0.267949192431123D0, beta = 0.464101615137755D0
real(kind=dbl),parameter          :: eps = 1.0D-6

x = rod%r_copyd()
r(1:3) = x(1:3) * x(4)

res = .FALSE.

if ( (r(2).ge.0.D0).and.(r(3).ge.0.D0) ) then 
  if ( ((alpha * (r(1)+r(3)) + r(2)) - beta .le. eps).and.( (alpha * (r(2)-r(3)) + r(1)) - beta .le. eps) ) then 
    if ( (r(1) - r1 .le. eps) .and. (r(2) - r1 .le. eps) .and. (r(3) - r2 .le. eps) ) res = .TRUE.
  end if
end if 

end function insideCubeHexFZ_

!--------------------------------------------------------------------------
recursive subroutine delete_FZlist_(self, l)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/21/20
  !!
  !! delete a linked list of rodrigues vectors

class(so3_T),INTENT(INOUT)        :: self 
character(2), INTENT(IN),OPTIONAL :: l    

type(FZpointd),pointer            :: ltail, ltmp

if (present(l)) then 
  select case(l)
    case('CM')
      ltail => self%CMlist
      self%CMcnt = 0
    case('FZ')
      ltail => self%FZlist
      self%CMcnt = 0
    case('CO')
      ltail => self%COlist
      self%COcnt = 0
    case('FB')
      ltail => self%FBlist
      self%FBcnt = 0
    case default 
      ltail => self%FZlist
      self%FZcnt = 0
    end select
else
  ltail => self%FZlist
  self%FZcnt = 0
end if

! deallocate the entire linked list before returning, to prevent memory leaks
ltmp => ltail % next
do 
  if (associated(ltail)) deallocate(ltail)
  if (.not. associated(ltmp)) EXIT
  ltail => ltmp
  ltmp => ltail % next
end do

end subroutine delete_FZlist_

!--------------------------------------------------------------------------
recursive subroutine SampleRFZ_(self, nsteps, qFZ)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/21/20
  !!
  !! Generate a uniform sampling of a Rodriguess FZ
  !!
  !! This routine fills in a linked list FZlist of Rodrigues points that 
  !! are inside a specific fundamental zone determined by the sample point group;
  !! this list can then be further dealt with in the calling program.  
  !!
  !! Here's how you would use this routine in a main program:
  !!
  !!    use mod_so3
  !!    use mod_rotations
  !!   
  !!    IMPLICIT NONE 
  !!   
  !!    type(so3_T)             :: SO
  !!    integer(kind=irg)       :: i, pgnum, nsteps, FZcnt
  !!    type(FZpointd), pointer :: FZtmp
  !!    type(e_T)               :: eu
  !!    
  !!    pgnum = 32
  !!    SO = so3_T( pgnum )
  !!    
  !!    nsteps = 10
  !!    call SO%sampleRFZ(nsteps)
  !! 
  !! Then you can access all the entries in the list and, for instance, convert them to Euler angles...
  !!
  !!    FZtmp => SO%getListHead('FZ')          ! point to the top of the list
  !!    FZcnt = SO%getListCount('FZ')          ! get the number of entries in the list
  !!    do i = 1, FZcnt                        ! loop over all entries
  !!      eu = FZtmp%rod%re()                  ! convert to Euler angles (in radians by default)
  !!    !  do something with eu                ! for instance, write eu to a file
  !!      FZtmp => FZtmp%next                  ! point to the next entry
  !!    end do
  !!
  !! If you just want to look at the first 10 entries on the list and show all other orientation representations:
  !!
  !!    type(orientation_T) :: ot
  !!    
  !!    FZtmp => SO%getListHead('FZ')
  !!    do i = 1,10
  !!      ot = orientation_T( FZtmp%rod )
  !!      call ot%print_orientation('d')    ! the argument 'd' means angles will be in degrees ('r' for radians)
  !!      FZtmp => FZtmp%next
  !!    end do

IMPLICIT NONE

class(so3_T),INTENT(INOUT)           :: self 

integer(kind=irg), INTENT(IN)        :: nsteps
type(q_T),INTENT(INOUT),OPTIONAL     :: qFZ

type(r_T)                            :: rod
type(c_T)                            :: cu
real(kind=dbl)                       :: x, y, z, delta, shift, sedge, ztmp
type(FZpointd), pointer              :: FZtmp, FZtmp2
integer(kind=irg)                    :: i, j, k
logical                              :: b, rotateFZ = .FALSE.

if (present(qFZ)) rotateFZ = .TRUE.

! cube semi-edge length s = 0.5D0 * LPs%ap
! step size for sampling of grid; total number of samples = (2*nsteps+1)**3
sedge = 0.5D0 * LPs%ap
delta = sedge / dble(nsteps)
if (self%gridtype.eq.0) then
  shift = 0.0D0
else
  shift = 0.5D0
end if

! set the counter to zero
self%FZcnt = 0

! note that when FZtype is cyclic (1) and FZorder is 2, then we must rotate the 
! rotation axis to lie along the b (y) direction, not z !!!!

! loop over the cube of volume pi^2; note that we do not want to include
! the opposite edges/facets of the cube, to avoid double counting rotations
! with a rotation angle of 180 degrees.  This only affects the cyclic groups.

 do i=-nsteps+1,nsteps
  x = (dble(i)+shift)*delta
  do j=-nsteps+1,nsteps
   y = (dble(j)+shift)*delta
   do k=-nsteps+1,nsteps
    z = (dble(k)+shift)*delta
! make sure that this point lies inside the cubochoric cell
    if (maxval( (/ abs(x), abs(y), abs(z) /) ).le.sedge) then

! convert to Rodrigues representation
      cu = c_T( cdinp = (/ x, y, z /) )
      rod = cu%cr()

! If insideFZ=.TRUE., then add this point to the linked list FZlist and keep
! track of how many points there are on this list
       if (rotateFZ.eqv..TRUE.) then 
         b = self%IsinsideFZ(rod, qFZ)
       else
         b = self%IsinsideFZ(rod)
       end if
       if (b) then 
        if (.not.associated(self%FZlist)) then
          allocate(self%FZlist)
          FZtmp => self%FZlist
        else
          allocate(FZtmp%next)
          FZtmp => FZtmp%next
        end if
        nullify(FZtmp%next)
! if monoclinic, then reorder the components !!!
!        if ((FZtype.eq.1).and.(FZorder.eq.2)) then
!          ztmp = rod(3)
!          rod(3) = rod(1)
!          rod(1) = rod(2)
!          rod(2) = ztmp
!        end if
        FZtmp%rod = rod
        FZtmp%gridpt(1:3) = (/i, j, k/)
        self%FZcnt = self%FZcnt + 1
       end if
    end if
  end do
 end do
end do

end subroutine SampleRFZ_

!--------------------------------------------------------------------------
recursive subroutine CubochoricNeighbors_(self, cubneighbor, nn, cub, stepsize)
  !! author: Saransh Singh
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! find the nearest neighbors of a point in s03 space, given the point
  !! and the step size used in the previous meshing. to be used in multi resolution
  !! indexing programs, specifically the PED, ECP and EBSD indexing. we're not worrying
  !! about keeping the neighbors in the FZ. that can just be done later.

use mod_io 

IMPLICIT NONE

class(so3_T),INTENT(INOUT)           :: self 

integer(kind=irg),INTENT(IN)         :: nn ! number of nearest neighbor in each direction (should be an odd number for symmetric meshing)
real(kind=dbl),INTENT(OUT)           :: cubneighbor(3,(2*nn+1)**3)
real(kind=dbl),INTENT(IN)            :: cub(3)
real(kind=dbl),INTENT(IN)            :: stepsize

type(IO_T)                           :: Message
integer(kind=irg)                    :: ii,jj,kk,ll,idx

if (dabs(stepsize) .gt. LPs%ap) then
    call Message%printError('CubochoricNeighbors', 'Step size is larger than edge length of the cube')
end if

do ii = -nn,nn
    do jj = -nn,nn
        do kk = -nn,nn
            idx  = (ii + nn)*(2*nn + 1)**2 + (jj + nn)*(2*nn + 1) + (kk + nn + 1)
            cubneighbor(1:3,idx) = cub + stepsize/2.D0*(/ii,jj,kk/)
            do ll = 1,3
                if (cubneighbor(ll,idx) .lt.  -0.5D0 * LPs%ap) then
                    cubneighbor(ll,idx) = cubneighbor(ll,idx) + LPs%ap
                else if (cubneighbor(ll,idx) .gt.  0.5D0 * LPs%ap) then
                    cubneighbor(ll,idx) = cubneighbor(ll,idx) - LPs%ap
                end if
            end do
        end do
    end do
end do

end subroutine CubochoricNeighbors_

!--------------------------------------------------------------------------
recursive subroutine sample_isoCube_(self, misang, N) ! CM = Constant Misorientation
  !! author: MDG
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! sample a centered cube surface inside the Cubochoric cube for a given misorientation angle
  !!
  !! linked list (self%CMlist) will have a length of 6(N-1)^2+2 entries

IMPLICIT NONE

class(so3_T),INTENT(INOUT)              :: self 

real(kind=dbl),INTENT(IN)               :: misang       
 !! desired misorientation angle (degrees)
integer(kind=irg),INTENT(IN)            :: N            
 !! desired number of sampling points along cube edge

type(FZpointd),pointer                  :: CMtmp, CMtmp2
type(c_T)                               :: cu
real(kind=dbl)                          :: edge, misangr, dx, x, y, z
integer(kind=irg)                       :: i, j, k

! initialize parameters
self%CMcnt = 0

! convert the misorientation angle to radians
misangr = misang * dtor 

! make sure the linked list is empty
if (associated(self%CMlist)) call self%delete_FZlist('CM')

! allocate the linked list
allocate(self%CMlist)
CMtmp => self%CMlist

! set the cube edge length based on the misorientation angle
edge = (cPi * (misangr-sin(misangr)))**(1.D0/3.D0)  * 0.5D0
dx = edge / dble(N)

! and generate the linked list of surface points

! do the x-y bottom and top planes first (each have N^2 points)
do i=-N,N
  x = dble(i)*dx
  do j=-N,N
    y = dble(j)*dx
! add the point to the list 
    cu = c_T( cdinp = (/ x, y, edge /) )
    CMtmp%rod = cu%cr()
    self%CMcnt = self%CMcnt + 1
    allocate(CMtmp%next)
    CMtmp => CMtmp%next
    nullify(CMtmp%next)
! and its mirror image in the top plane
    cu = c_T( cdinp = (/ x, y, -edge /) )
    CMtmp%rod = cu%cr()
    self%CMcnt = self%CMcnt + 1
    allocate(CMtmp%next)
    CMtmp => CMtmp%next
    nullify(CMtmp%next)
  end do
end do

! then we do the y-z planes; each have N*(N-2) points
do j=-N,N
  y =  dble(j)*dx
  do k=-N+1,N-1
    z = dble(k)*dx
! add the point to the list 
    cu = c_T( cdinp = (/ edge, y, z /) )
    CMtmp%rod = cu%cr()
    self%CMcnt = self%CMcnt + 1
    allocate(CMtmp%next)
    CMtmp => CMtmp%next
    nullify(CMtmp%next)
! and its mirror image in the top plane
    cu = c_T( cdinp = (/ -edge, y, z /) )
    CMtmp%rod = cu%cr()
    self%CMcnt = self%CMcnt + 1
    allocate(CMtmp%next)
    CMtmp => CMtmp%next
    nullify(CMtmp%next)
  end do
end do

! and finally the x-z planes, with (N-2)^2 points each
do i=-N+1,N-1
  x = dble(i)*dx
  do k=-N+1,N-1
    z = dble(k)*dx
! add the point to the list 
    cu = c_T( cdinp = (/ x, edge, z /) )
    CMtmp%rod = cu%cr()
    self%CMcnt = self%CMcnt + 1
    allocate(CMtmp%next)
    CMtmp => CMtmp%next
    nullify(CMtmp%next)
! and its mirror image in the top plane
    cu = c_T( cdinp = (/ x, -edge, z /) )
    CMtmp%rod = cu%cr()
    self%CMcnt = self%CMcnt + 1
    allocate(CMtmp%next)
    CMtmp => CMtmp%next
    nullify(CMtmp%next)
  end do
end do

end subroutine sample_isoCube_

!--------------------------------------------------------------------------
recursive subroutine sample_isoCubeFilled_(self, misang, N) ! CM = Constant Misorientation
  !! author: MDG
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! sample a centered cube inside the Cubochoric cube for a given misorientation angle
  !!  
  !! This routine is different from the sample_isoCube routine in that it 
  !! generates ALL the points inside the centered cube instead of just the points on
  !! the outer surface.  This can be useful to uniformly sample a small volume of orientation
  !! space around some point out to a given misorientation angle.  Since the sampling has concentric
  !! cubes, all the samples can be subdivided into discrete misorientation classes.
  !! The linked list wil have a length of N^3

IMPLICIT NONE

class(so3_T),INTENT(INOUT)              :: self 

real(kind=dbl),INTENT(IN)               :: misang       
 !! desired misorientation angle (degrees)
integer(kind=irg),INTENT(IN)            :: N            
 !! desired number of sampling points along cube edge

type(FZpointd),pointer                  :: CMtmp, CMtmp2
type(c_T)                               :: cu
real(kind=dbl)                          :: edge, misangr, dx, x, y, z, xc, yc, zc
integer(kind=irg)                       :: i, j, k

! initialize parameters
self%CMcnt = 0

! convert the misorientation angle to radians
misangr = misang * dtor

! make sure the linked list is empty
if (associated(self%CMlist)) call self%delete_FZlist('CM') 

! allocate the linked list
allocate(self%CMlist)
CMtmp => self%CMlist

! set the cube edge length based on the misorientation angle
edge = (cPi * (misangr-sin(misangr)))**(1.D0/3.D0) * 0.5D0
dx = edge / dble(N)

! and generate the linked list of surface points
! loop over the (2N+1)^3 points
do i=-N,N
  x = dble(i)*dx
  do j=-N,N
    y = dble(j)*dx
    do k=-N,N
      z = dble(k)*dx
! add the point to the list 
      cu = c_T( cdinp = (/ x, y, z /) )
      CMtmp%rod = cu%cr()
      self%CMcnt = self%CMcnt + 1
      allocate(CMtmp%next)
      CMtmp => CMtmp%next
      nullify(CMtmp%next)
    end do
  end do
end do

end subroutine sample_isoCubeFilled_

!--------------------------------------------------------------------------
recursive subroutine sample_Cone_(self, unitvector, dpmin, N) 
  !! author: MDG
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! sample a cone centered on a unit vector with apex in the origin and given opening angle

IMPLICIT NONE

class(so3_T),INTENT(INOUT)              :: self 

real(kind=dbl),INTENT(IN)               :: unitvector(3)
 !! axis of cone
real(kind=dbl),INTENT(IN)               :: dpmin        
 !! maximum dot product
integer(kind=irg),INTENT(IN)            :: N            
 !! number of sampling points along cube semi edge

type(FZpointd),pointer                  :: tmp, tmp2
type(r_T)                               :: rod 
type(c_T)                               :: cu
real(kind=dbl)                          :: dx, x, y, z, s, delta, dp, r(3), xx(4)

! initialize parameters
self%COcnt = 0
! cube semi-edge length
s = 0.5D0 * LPs%ap
! step size for sampling of grid; total number of samples = (2*nsteps+1)**3
delta = s/dble(N)

! make sure the linked list is empty
if (associated(self%COlist)) call self%delete_FZlist('CO')

! allocate the linked list and insert the origin
allocate(self%COlist)
tmp => self%COlist
nullify(tmp%next)
tmp%rod = r_T( rdinp = (/ 0.D0, 0.D0, 0.D0, 0.D0 /) )
self%COcnt = self%COcnt + 1

! and generate the linked list of points inside the cone
x = -s
do while (x.lt.s)
  y = -s
  do while (y.lt.s)
    z = -s
    do while (z.lt.s)

     if ((x.ne.0.D0).and.(y.ne.0.D0).and.(z.ne.0.D0)) then
! convert to Rodrigues representation
      cu = c_T( cdinp = (/ x, y, z /) )
      rod = cu%cr()
      xx = rod%r_copyd()
      r = xx(1:3)/sqrt(sum(xx(1:3)**2))
! compute the dot product of this vector and the unitvector
      dp = unitvector(1)*r(1)+unitvector(2)*r(2)+unitvector(3)*r(3)
! conditionally add the point to the list if it lies inside the cone (dpmax <= dp)
      if ((dp.ge.dpmin).and.(self%IsinsideFZ(rod))) then
        allocate(tmp%next)
        tmp => tmp%next
        nullify(tmp%next)
        tmp%trod = rod
        self%COcnt = self%COcnt + 1
      end if
     end if

    z = z + delta
  end do
  y = y + delta
 end do
 x = x + delta
end do

end subroutine sample_Cone_

!--------------------------------------------------------------------------
recursive subroutine sample_Fiber_(self, itmp, num, dpmin, N) 
  !! author: MDG
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! sample a fiber texture with a given fiber axis and angular spread

use mod_quaternions 

IMPLICIT NONE

class(so3_T),INTENT(INOUT)              :: self 

integer(kind=irg),INTENT(IN)            :: num
real(kind=dbl),INTENT(IN)               :: itmp(num,3)   
 !! equivalent fiber axis unit vectors
real(kind=dbl),INTENT(IN)               :: dpmin        
 !! maximum dot product
integer(kind=irg),INTENT(IN)            :: N            
 !! number of sampling points along cube semi edge

type(FZpointd),pointer                  :: tmp, tmp2
type(quaternion_T)                      :: qu 
type(r_T)                               :: rod
type(c_T)                               :: cu
type(q_T)                               :: q
real(kind=dbl)                          :: dx, x, y, z, s, delta, dp, Fr(3), Fn(3)
integer(kind=irg)                       :: j

! initialize parameters
self%FBcnt = 0
! cube semi-edge length
s = 0.5D0 * LPs%ap
! step size for sampling of grid; total number of samples = (2*nsteps+1)**3
delta = s/dble(N)

! make sure the linked list is empty
if (associated(self%FBlist)) call self%delete_FZlist('FB')

! allocate the linked list and insert the origin
allocate(self%FBlist)
tmp => self%FBlist
nullify(tmp%next)

! and generate the linked list of points inside the cone
x = -s
do while (x.lt.s)
  y = -s
  do while (y.lt.s)
    z = -s
    do while (z.lt.s)

! convert to Rodrigues representation
      cu = c_T( cdinp = (/ x, y, z /) )
      q = cu%cq()
      rod = cu%cr()
      qu = conjg( quaternion_T( qd = rod%r_copyd() ) )

! loop over the equivalent fiber axis indices
      do j=1,num
        Fr = qu%quat_Lp( itmp(j,1:3) )

! conditionally add the point to the list if it lies inside the cone (dpmax <= dp)
        if ((Fr(3).ge.dpmin).and.(self%IsinsideFZ(rod))) then
          tmp%trod = rod
          allocate(tmp%next)
          tmp => tmp%next
          nullify(tmp%next)
          self%FBcnt = self%FBcnt + 1
        end if
      end do

    z = z + delta
  end do
  y = y + delta
 end do
 x = x + delta
end do

end subroutine sample_Fiber_

!--------------------------------------------------------------------------
recursive subroutine SampleIsoMisorientation_(self, rhozero, misang) 
  !! author: MDG
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! Constant Misorientation sampling routine; input list must be generated by sampleCubeSurface

use mod_math

IMPLICIT NONE

class(so3_T),INTENT(INOUT)              :: self 

type(r_T),INTENT(INOUT)                 :: rhozero  
 !! center Rodrigues vector
real(kind=dbl),INTENT(IN)               :: misang       
 !! desired misorientation angle (degrees)

type(FZpointd),pointer                  :: CMtmp
real(kind=dbl)                          :: rhovec(3), s, vv(3), x(4)
integer(kind=irg)                       :: i

! go through the list and transform all points to the spheroid misorientation surface
! the resulting Rodrigues vectors are stored in the trod(4) entry.

x = rhozero%r_copyd()
rhovec(1:3) = x(1:3) * x(4)

CMtmp => self%CMlist
do i=1,self%CMcnt
! get the actual Rodrigues vector
  x = CMtmp%rod%r_copyd()
  vv(1:3) = x(1:3) * x(4)
! apply the Rodrigues transformation formula
  vv = (-vv + rhovec + cross3(rhovec,vv))/(1.D0 + DOT_PRODUCT(vv,rhovec))
! and convert back to the 4-component format
  s = dsqrt(sum(vv*vv))
  if (s.gt.0.D0) then 
    CMtmp%trod = r_T( rdinp = (/ vv(1)/s, vv(2)/s, vv(3)/s, s /) )
  else
    CMtmp%trod = r_T( rdinp = (/ 0.D0, 0.D0, 1.D0, 0.D0 /) )
  end if
  CMtmp=>CMtmp%next
end do

end subroutine SampleIsoMisorientation_

!--------------------------------------------------------------------------
recursive subroutine getOrientationsfromFile_(self, filename)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! read a list of orientations from a text file and insert them in a linked list

IMPLICIT NONE

class(so3_T),INTENT(INOUT)              :: self 

character(fnlen),INTENT(IN)             :: filename
 !! complete path to input file name 

type(e_T)                               :: e 
type(o_T)                               :: o 
type(a_T)                               :: a 
type(r_T)                               :: r 
type(q_T)                               :: q 
type(h_T)                               :: h 
type(c_T)                               :: c 
type(s_T)                               :: s 
type(v_T)                               :: v 

character(2)                            :: anglemode
integer(kind=irg)                       :: numang, i
real(kind=dbl)                          :: x3(3), x4(4), x9(9)
type(FZpointd),pointer                  :: FZtmp

open(unit=53,file=trim(filename),status='old',action='read')
read (53,*) anglemode
read (53,*) numang

! make sure the linked list is empty
if (associated(self%FZlist)) call self%delete_FZlist('FZ')

! allocate the linked list
allocate(self%FZlist)
FZtmp => self%FZlist

select case(anglemode)
  case('eu') ! angles must be in degrees
    do i=1,numang
      read (53,*) x3(1:3)
      x3 = x3 * dtor
      e = e_T( edinp = x3 )
      FZtmp%rod = e%er()
      self%FZcnt = self%FZcnt + 1
      allocate(FZtmp%next)
      FZtmp => FZtmp%next
      nullify(FZtmp%next)
    end do
  case('ro')
    do i=1,numang
      read (53,*) x4(1:4)
      FZtmp%rod = r_T( rdinp = x4 )
      self%FZcnt = self%FZcnt + 1
      allocate(FZtmp%next)
      FZtmp => FZtmp%next
      nullify(FZtmp%next)
    end do
  case('qu') 
    do i=1,numang
      read (53,*) x4(1:4)
      q = q_T( qdinp = x4 )
      FZtmp%rod = q%qr()
      self%FZcnt = self%FZcnt + 1
      allocate(FZtmp%next)
      FZtmp => FZtmp%next
      nullify(FZtmp%next)
    end do
  case('ax') ! angle must be in degrees
    do i=1,numang
      read (53,*) x4(1:4)
      x4(4) = x4(4) * dtor
      a = a_T( adinp = x4 )
      FZtmp%rod = a%ar()
      self%FZcnt = self%FZcnt + 1
      allocate(FZtmp%next)
      FZtmp => FZtmp%next
      nullify(FZtmp%next)
    end do
  case('ho') 
    do i=1,numang
      read (53,*) x3(1:3)
      h = h_T( hdinp = x3 )
      FZtmp%rod = h%hr()
      self%FZcnt = self%FZcnt + 1
      allocate(FZtmp%next)
      FZtmp => FZtmp%next
      nullify(FZtmp%next)
    end do
  case('cu')
    do i=1,numang
      read (53,*) x3(1:3)
      c = c_T( cdinp = x3 )
      FZtmp%rod = c%cr()
      self%FZcnt = self%FZcnt + 1
      allocate(FZtmp%next)
      FZtmp => FZtmp%next
      nullify(FZtmp%next)
    end do
  case('st')
    do i=1,numang
      read (53,*) x3(1:3)
      s = s_T( sdinp = x3 )
      FZtmp%rod = s%sr()
      self%FZcnt = self%FZcnt + 1
      allocate(FZtmp%next)
      FZtmp => FZtmp%next
      nullify(FZtmp%next)
    end do
  case('om')
    do i=1,numang
      read (53,*) x9(1:9)
      o = o_T( odinp = reshape( x9, (/3,3/) ) )
      FZtmp%rod = o%or()
      self%FZcnt = self%FZcnt + 1
      allocate(FZtmp%next)
      FZtmp => FZtmp%next
      nullify(FZtmp%next)
    end do
  case('rv')
    do i=1,numang
      read (53,*) x3(1:3)
      v = v_T( vdinp = x3 )
      FZtmp%rod = v%vr()
      self%FZcnt = self%FZcnt + 1
      allocate(FZtmp%next)
      FZtmp => FZtmp%next
      nullify(FZtmp%next)
    end do
end select

close(unit=53,status='keep')

end subroutine getOrientationsfromFile_

!--------------------------------------------------------------------------
recursive subroutine writeOrientationstoFile_(self, filename, mode, list, trod)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! write a list of orientations from a linked list to a text file 

use mod_io 
use mod_math

IMPLICIT NONE

class(so3_T),INTENT(INOUT)              :: self 

character(fnlen),INTENT(IN)             :: filename
 !! complete path to output file name 
character(2), INTENT(IN)                :: mode
 !! output orientation representation  (eu, ro, ho, ...)
character(2), INTENT(IN), OPTIONAL      :: list 
 !! list from which to write 
logical, INTENT(IN), OPTIONAL           :: trod
 !! list from which to write 

type(e_T)                               :: e
type(o_T)                               :: o
type(q_T)                               :: q
type(s_T)                               :: s
type(v_T)                               :: v
type(h_T)                               :: h
type(c_T)                               :: c
type(r_T)                               :: r
type(a_T)                               :: a

type(IO_T)                              :: Message
type(FZpointd), pointer                 :: FZtmp 
integer(kind=irg)                       :: cnt, i
real(kind=dbl)                          :: io_real(9)
logical                                 :: dotrod 

dotrod = .FALSE.
if (present(trod)) then 
  if (trod.eqv..TRUE.) dotrod = .TRUE. 
end if

if (present(list)) then
  select case(list)
    case('FZ')
      FZtmp => self%FZlist
      cnt = self%FZcnt 
    case('CM')
      FZtmp => self%CMlist
      cnt = self%CMcnt 
    case('CO')
      FZtmp => self%COlist
      cnt = self%COcnt 
    case('FB')
      FZtmp => self%FBlist
      cnt = self%FBcnt 
    case default 
      FZtmp => self%FZlist
      cnt = self%FZcnt 
  end select 
else
  FZtmp => self%FZlist
  cnt = self%FZcnt
end if 

open(unit=53, file=trim(filename), status='unknown', form='formatted')
write (53,"(A2)") mode 
write (53,"(I6)") cnt 

if (dotrod.eqv..TRUE.) then 
  do i=1, cnt
    select case(mode)
      case('eu')
        e = FZtmp%trod%re()
        io_real(1:3) = e%e_copyd() / dtor 
        call Message%WriteValue('', io_real, 3, frm="(2(F14.8,' '),F14.8)",redirect=53)
      case('ro')
        io_real(1:4) = FZtmp%trod%r_copyd() 
        if (io_real(4).eq.inftyd()) then 
          call Message%WriteValue('', io_real, 3, frm="(3(F14.8,' '),'infinity')",redirect=53)
        else
          call Message%WriteValue('', io_real, 4, frm="(3(F14.8,' '),F14.8)",redirect=53)
        end if 
      case('om')
        o = FZtmp%trod%ro()
        io_real(1:9) = reshape(o%o_copyd(), (/ 9 /) ) 
        call Message%WriteValue('', io_real, 9, frm="(8(F14.8,' '),F14.8)",redirect=53)
      case('ho')
        h = FZtmp%trod%rh()
        io_real(1:3) = h%h_copyd()
        call Message%WriteValue('', io_real, 3, frm="(2(F14.8,' '),F14.8)",redirect=53)
      case('cu')
        c = FZtmp%trod%rc()
        io_real(1:3) = c%c_copyd()
        call Message%WriteValue('', io_real, 3, frm="(2(F14.8,' '),F14.8)",redirect=53)
      case('rv')
        v = FZtmp%trod%rv()
        io_real(1:3) = v%v_copyd() 
        call Message%WriteValue('', io_real, 3, frm="(2(F14.8,' '),F14.8)",redirect=53)
      case('st')
        s = FZtmp%trod%rs()
        io_real(1:3) = s%s_copyd()
        call Message%WriteValue('', io_real, 3, frm="(2(F14.8,' '),F14.8)",redirect=53)
      case('ax')
        a = FZtmp%trod%ra()
        io_real(1:4) = a%a_copyd()
        io_real(4) = io_real(4) / dtor 
        call Message%WriteValue('', io_real, 4, frm="(3(F14.8,' '),F14.8)",redirect=53)
      case('qu')
        q = FZtmp%trod%rq()
        io_real(1:4) = q%q_copyd() 
        call Message%WriteValue('', io_real, 4, frm="(3(F14.8,' '),F14.8)",redirect=53)
      case default
    end select
    FZtmp => FZtmp%next 
  end do 
else
  do i=1, cnt
    select case(mode)
      case('eu')
        e = FZtmp%rod%re()
        io_real(1:3) = e%e_copyd() / dtor 
        call Message%WriteValue('', io_real, 3, frm="(2(F14.8,' '),F14.8)",redirect=53)
      case('ro')
        io_real(1:4) = FZtmp%rod%r_copyd() 
        if (io_real(4).eq.inftyd()) then 
          call Message%WriteValue('', io_real, 3, frm="(3(F14.8,' '),'infinity')",redirect=53)
        else
          call Message%WriteValue('', io_real, 4, frm="(3(F14.8,' '),F14.8)",redirect=53)
        end if 
      case('om')
        o = FZtmp%rod%ro()
        io_real(1:9) = reshape(o%o_copyd(), (/ 9 /) ) 
        call Message%WriteValue('', io_real, 9, frm="(8(F14.8,' '),F14.8)",redirect=53)
      case('ho')
        h = FZtmp%rod%rh()
        io_real(1:3) = h%h_copyd()
        call Message%WriteValue('', io_real, 3, frm="(2(F14.8,' '),F14.8)",redirect=53)
      case('cu')
        c = FZtmp%rod%rc()
        io_real(1:3) = c%c_copyd()
        call Message%WriteValue('', io_real, 3, frm="(2(F14.8,' '),F14.8)",redirect=53)
      case('rv')
        v = FZtmp%rod%rv()
        io_real(1:3) = v%v_copyd() 
        call Message%WriteValue('', io_real, 3, frm="(2(F14.8,' '),F14.8)",redirect=53)
      case('st')
        s = FZtmp%rod%rs()
        io_real(1:3) = s%s_copyd()
        call Message%WriteValue('', io_real, 3, frm="(2(F14.8,' '),F14.8)",redirect=53)
      case('ax')
        a = FZtmp%rod%ra()
        io_real(1:4) = a%a_copyd()
        io_real(4) = io_real(4) / dtor 
        call Message%WriteValue('', io_real, 4, frm="(3(F14.8,' '),F14.8)",redirect=53)
      case('qu')
        q = FZtmp%rod%rq()
        io_real(1:4) = q%q_copyd() 
        call Message%WriteValue('', io_real, 4, frm="(3(F14.8,' '),F14.8)",redirect=53)
      case default
    end select
    FZtmp => FZtmp%next 
  end do 
end if  

close(unit=53, status = 'keep')

end subroutine writeOrientationstoFile_

!--------------------------------------------------------------------------
recursive subroutine listtoArray_(self, l, eAR, oAR, qAR, sAR, vAR, hAR, cAR, rAR, aAR) 
  !! author: MDG
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! convert a linked list into an array of orientations 

IMPLICIT NONE

class(so3_T),INTENT(INOUT)              :: self 
character(2), INTENT(IN), OPTIONAL      :: l 
type(e_T),allocatable,OPTIONAL          :: eAR(:)
type(o_T),allocatable,OPTIONAL          :: oAR(:)
type(q_T),allocatable,OPTIONAL          :: qAR(:)
type(s_T),allocatable,OPTIONAL          :: sAR(:)
type(v_T),allocatable,OPTIONAL          :: vAR(:)
type(h_T),allocatable,OPTIONAL          :: hAR(:)
type(c_T),allocatable,OPTIONAL          :: cAR(:)
type(r_T),allocatable,OPTIONAL          :: rAR(:)
type(a_T),allocatable,OPTIONAL          :: aAR(:)

integer(kind=irg)                       :: cnt, i
type(FZpointd), pointer                 :: FZtmp 

if (present(l)) then 
  select case(l)
    case('FZ')
      FZtmp => self%FZlist
      cnt = self%FZcnt 
    case('CM')
      FZtmp => self%CMlist
      cnt = self%CMcnt 
    case('CO')
      FZtmp => self%COlist
      cnt = self%COcnt 
    case('FB')
      FZtmp => self%FBlist
      cnt = self%FBcnt 
    case default 
      FZtmp => self%FZlist
      cnt = self%FZcnt 
  end select 
else 
  FZtmp => self%FZlist
  cnt = self%FZcnt 
end if 

if (present(eAR)) then 
  if (allocated(eAR)) deallocate(eAR)
  allocate(eAR(cnt))
end if
if (present(oAR)) then 
  if (allocated(oAR)) deallocate(oAR)
  allocate(oAR(cnt))
end if
if (present(qAR)) then 
  if (allocated(qAR)) deallocate(qAR)
  allocate(qAR(cnt))
end if
if (present(sAR)) then 
  if (allocated(sAR)) deallocate(sAR)
  allocate(sAR(cnt))
end if
if (present(vAR)) then 
  if (allocated(vAR)) deallocate(vAR)
  allocate(vAR(cnt))
end if
if (present(hAR)) then 
  if (allocated(hAR)) deallocate(hAR)
  allocate(hAR(cnt))
end if
if (present(cAR)) then 
  if (allocated(cAR)) deallocate(cAR)
  allocate(cAR(cnt))
end if
if (present(rAR)) then 
  if (allocated(rAR)) deallocate(rAR)
  allocate(rAR(cnt))
end if
if (present(aAR)) then 
  if (allocated(aAR)) deallocate(aAR)
  allocate(aAR(cnt))
end if

do i=1,cnt 
  if (present(eAR)) eAR(i) = FZtmp%rod%re() 
  if (present(oAR)) oAR(i) = FZtmp%rod%ro() 
  if (present(qAR)) qAR(i) = FZtmp%rod%rq() 
  if (present(sAR)) sAR(i) = FZtmp%rod%rs() 
  if (present(vAR)) vAR(i) = FZtmp%rod%rv() 
  if (present(hAR)) hAR(i) = FZtmp%rod%rh() 
  if (present(cAR)) cAR(i) = FZtmp%rod%rc() 
  if (present(rAR)) rAR(i) = FZtmp%rod 
  if (present(aAR)) aAR(i) = FZtmp%rod%ra() 
  FZtmp => FZtmp%next 
end do

end subroutine listtoArray_

!--------------------------------------------------------------------------
recursive subroutine listtoQuaternionArray_(self, qAR, l)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/24/20
  !!
  !! convert a linked list into a QuaternionArray_T object 

use mod_quaternions 
use mod_rotations

IMPLICIT NONE

class(so3_T),INTENT(INOUT)              :: self 
type(QuaternionArray_T), INTENT(OUT)    :: qAR
character(2), INTENT(IN), OPTIONAL      :: l 

type(FZpointd), pointer                 :: FZtmp
integer(kind=irg)                       :: i, cnt
type(q_T)                               :: qu
type(r_T)                               :: rod
type(Quaternion_T)                      :: qq

if (present(l)) then 
  select case(l)
    case('FZ')
      FZtmp => self%FZlist
      cnt = self%FZcnt 
    case('CM')
      FZtmp => self%CMlist
      cnt = self%CMcnt 
    case('CO')
      FZtmp => self%COlist
      cnt = self%COcnt 
    case('FB')
      FZtmp => self%FBlist
      cnt = self%FBcnt 
    case default 
      FZtmp => self%FZlist
      cnt = self%FZcnt 
  end select 
else 
  FZtmp => self%FZlist
  cnt = self%FZcnt 
end if 

qAR = QuaternionArray_T( n = cnt, s='d' )
qq = Quaternion_T()

do i=1,cnt 
  rod = FZtmp%rod
  qu = rod%rq()
  call qq%set_quatd(qu%q_copyd())
  call qAR%insertQuatinArray( i, qq )
  FZtmp => FZtmp%next 
end do

end subroutine listtoQuaternionArray_

!--------------------------------------------------------------------------
recursive subroutine QuaternionArraytolist_(self, qAR, l)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/24/20
  !!
  !! convert a QuaternionArray_T object into a linked list

use mod_quaternions 
use mod_rotations
use mod_io

IMPLICIT NONE

class(so3_T),INTENT(INOUT)              :: self 
type(QuaternionArray_T), INTENT(INOUT)  :: qAR
character(2), INTENT(IN), OPTIONAL      :: l 

type(IO_T)                              :: Message 
type(FZpointd), pointer                 :: FZtmp
integer(kind=irg)                       :: i, cnt, N
type(Quaternion_T)                      :: qu
type(q_T)                               :: qq

if (present(l)) then 
  select case(l)
    case('FZ')
      FZtmp => self%FZlist
      cnt = self%FZcnt 
    case('CM')
      FZtmp => self%CMlist
      cnt = self%CMcnt 
    case('CO')
      FZtmp => self%COlist
      cnt = self%COcnt 
    case('FB')
      FZtmp => self%FBlist
      cnt = self%FBcnt 
    case default 
      FZtmp => self%FZlist
      cnt = self%FZcnt 
  end select 
else 
  FZtmp => self%FZlist
  cnt = self%FZcnt 
end if 

N = qAR%getQnumber()

! if there are fewer than cnt items in the Array, then we only copy those into the list 
if (N.lt.cnt) then 
  cnt = N 
  call Message%printWarning('QuaternionArraytolist: there are fewer items in array than in list; only copying those...')
end if 

do i=1,cnt 
  qu = qAR%getQuatfromArray(i)
  qq = q_T( qdinp = qu%get_quatd() )
  FZtmp%rod = qq%qr()
  FZtmp => FZtmp%next 
end do

end subroutine QuaternionArraytolist_

!--------------------------------------------------------------------------
recursive function getListHead_(self, l) result(FZptr)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! return the pointer to the selected list  

IMPLICIT NONE

class(so3_T),INTENT(INOUT)    :: self
character(2), INTENT(IN)      :: l  

type(FZpointd),pointer        :: FZptr 

select case(l)
  case('FZ')
    FZptr => self%FZlist 
  case('CM')
    FZptr => self%CMlist 
  case('CO')
    FZptr => self%COlist 
  case('FB')
    FZptr => self%FBlist 
  case default
    FZptr => self%FZlist 
end select 

end function getListHead_

!--------------------------------------------------------------------------
recursive function getListCount_(self, l) result(cnt)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! return the counter for the selected list  

IMPLICIT NONE

class(so3_T),INTENT(INOUT)    :: self
character(2), INTENT(IN)      :: l  

integer(kind=irg)             :: cnt 

select case(l)
  case('FZ')
    cnt = self%FZcnt 
  case('CM')
    cnt = self%CMcnt 
  case('CO')
    cnt = self%COcnt 
  case('FB')
    cnt = self%FBcnt 
  case default
    cnt = self%FZcnt 
end select 

end function getListCount_

!--------------------------------------------------------------------------
recursive subroutine setGridType_(self, g)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! set the gridtype parameter

IMPLICIT NONE

class(so3_T),INTENT(INOUT)    :: self
integer(kind=irg)             :: g 

self%gridtype = g 

end subroutine setGridType_

!--------------------------------------------------------------------------
recursive function IsinsideMFZ_(self, rod) result(insideMFZ)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! does Rodrigues point lie inside the relevant Mackenzie (disorientation) FZ

use mod_math

IMPLICIT NONE 

class(so3_T),INTENT(INOUT)    :: self

type(r_T), INTENT(INOUT)      :: rod
logical                       :: insideMFZ

real(kind=dbl)                :: r(4)

r = rod%r_copyd()

select case (self%MFZtype)
  case (0)
    insideMFZ = .TRUE.   ! all points are inside the FZ
  case (1)
    insideMFZ = self%insideCyclicFZ(rod, M=.TRUE.) ! infinity is checked inside this function
  case (2)
    if (r(4).ne.inftyd()) insideMFZ = self%insideDihedralMFZ(rod)
  case (3)
    if (r(4).ne.inftyd()) insideMFZ = self%insideCubicMFZ(rod,'tet')
  case (4)
    if (r(4).ne.inftyd()) insideMFZ = self%insideCubicMFZ(rod,'oct')
end select

end function IsinsideMFZ_

!--------------------------------------------------------------------------
recursive function insideCubicMFZ_(self, rod, ot) result(res)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! does Rodrigues point lie inside cubic MacKenzie FZ (octahedral or tetrahedral)?

IMPLICIT NONE

class(so3_T),INTENT(INOUT)    :: self

type(r_T), INTENT(INOUT)      :: rod
character(3), INTENT(IN)      :: ot

logical                       :: res, c0, c1, c2, c3
real(kind=dbl)                :: r(3), x(4)

res = .FALSE.

! first of all, we need to be inside the regular FZ
c0 = self%insideCubicFz(rod,ot)

x = rod%r_copyd()
r(1:3) = x(1:3) * x(4)

if (ot.eq.'oct') then
  c1 = (c0.and.(r(3).ge.0.D0))
  c2 = (c1.and.(r(2).ge.r(3)))
  c3 = (c2.and.(r(1).ge.r(2))) 
else 
  c1 = (c0.and.(minval(r).ge.0.D0))  ! in the first octant
  if (r(1).le.r(2)) then
    c3 = (c1.and.(r(3).le.r(1)))
  else
    c3 = (c1.and.(r(3).le.r(2)))
  end if
end if

res = c3

end function insideCubicMFZ_

!--------------------------------------------------------------------------
recursive function insideDihedralMFZ_(self, rod) result(res)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! does Rodrigues point lie inside dihedral MacKenzie FZ

IMPLICIT NONE

class(so3_T),INTENT(INOUT)    :: self

type(r_T), INTENT(INOUT)         :: rod

logical                       :: res, c0, c1, c2, c3
real(kind=dbl)                :: r(3), x(4)
real(kind=dbl),parameter      :: v = 0.57735026918962584D0

res = .FALSE.

! first of all, we need to be inside the regular FZ
c0 = self%insideDihedralFZ(rod, self%MFZorder)

x = rod%r_copyd()
r(1:3) = x(1:3) * x(4)

if (c0) then
select case (self%MFZorder)
    case (2)
      c2 = (minval(r).ge.0.D0)
    case (3)
      c1 = (minval( (/ r(1), r(3) /) ).ge.0.D0)
      c2 = (c1.and.(r(1).ge.dabs(r(2))*v)) 
    case (4)
      c1 = (minval(r).ge.0.D0)
      c2 = (c1.and.(r(1).ge.r(2))) 
    case (6)
      c1 = (minval(r).ge.0.D0)
      c2 = (c1.and.(r(1).ge.r(2)*v)) 
  end select
end if

res = c2

end function insideDihedralMFZ_

!--------------------------------------------------------------------------
!
! SUBROUTINE: getVertex
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief get vertices of RFZ for quasicrystals (dihedral symmetries)
!
!> @param order name of the Euler angle file (with usual path handling)
!> @param vertex the number of components in the returned linked list
!
!> @date 06/18/18 SS 1.0 original
!--------------------------------------------------------------------------
recursive subroutine getVertex_(self, order, vertex)
  !! author: Saransh Singh
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! get vertices of RFZ for quasicrystals (dihedral symmetries)

IMPLICIT NONE

class(so3_T),INTENT(INOUT)    :: self

integer(kind=irg),INTENT(IN)  :: order
real(kind=dbl),INTENT(OUT)    :: vertex(2*order,2)

integer(kind=irg)             :: ii
real(kind=dbl)                :: th

do ii = 1,2*order
  th  = (dble(ii - 1)/dble(order) + 1.D0/2.D0/dble(order)) * cPi
  vertex(ii,1:2) = (/dcos(th), dsin(th)/)
end do

end subroutine getVertex_
                                                                
!                                                                      
!     ..................................................................
!                                                                       
!        SUBROUTINE PNPOLY                                              
!                                                                       
!        PURPOSE                                                        
!           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON            
!                                                                       
!        USAGE                                                          
!           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )                     
!                                                                       
!        DESCRIPTION OF THE PARAMETERS                                  
!           PX      - X-COORDINATE OF POINT IN QUESTION.                
!           PY      - Y-COORDINATE OF POINT IN QUESTION.                
!           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF         
!                     VERTICES OF POLYGON.                              
!           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF           
!                     VERTICES OF POLYGON.                              
!           N       - NUMBER OF VERTICES IN THE POLYGON.                
!           INOUT   - THE SIGNAL RETURNED:                              
!                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,        
!                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,     
!                      1 IF THE POINT IS INSIDE OF THE POLYGON.         
!                                                                       
!        REMARKS                                                        
!           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.      
!           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY           
!           OPTIONALLY BE INCREASED BY 1.                               
!           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING      
!           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX    
!           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING   
!           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.              
!           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.         
!           THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM      
!           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.   
!                                                                       
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  
!           NONE                                                        
!                                                                       
!        METHOD                                                         
!           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT  
!           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE        
!           POINT IS INSIDE OF THE POLYGON.                             
!                                                                       
!     ..................................................................
!                                                                       
      RECURSIVE FUNCTION PNPOLY(PX,PY,XX,YY,N) RESULT(INOUT)

      IMPLICIT NONE

      REAL(KIND=DBL) PX, PY
      INTEGER(KIND=IRG) N

      REAL(KIND=DBL) X(200),Y(200),XX(N),YY(N)                                    
      LOGICAL MX,MY,NX,NY                                         
      INTEGER O,INOUT,I,J,MAXDIM                                                         
!      OUTPUT UNIT FOR PRINTED MESSAGES                                 
      DATA O/6/                                                         
      MAXDIM=200                                                        
      IF(N.LE.MAXDIM)GO TO 6                                            
      WRITE(O,7)                                                        
7     FORMAT('0WARNING:',I5,' TOO GREAT FOR THIS VERSION OF PNPOLY. RESULTS INVALID')                                                 
      RETURN                                                            
6     DO 1 I=1,N                                                        
      X(I)=XX(I)-PX                                                     
1     Y(I)=YY(I)-PY                                                     
      INOUT=-1                                                          
      DO 2 I=1,N                                                        
      J=1+MOD(I,N)                                                      
      MX=X(I).GE.0.0                                                    
      NX=X(J).GE.0.0                                                    
      MY=Y(I).GE.0.0                                                    
      NY=Y(J).GE.0.0                                                    
      IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2       
      IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3  
      INOUT=-INOUT                                                      
      GO TO 2                                                           
3     IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5                       
4     INOUT=0                                                           
      RETURN                                                            
5     INOUT=-INOUT                                                      
2     CONTINUE                                                          
      RETURN                                                            
      END FUNCTION PNPOLY


!--------------------------------------------------------------------------
recursive function MKCC(a, b, c) result(CC)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! auxiliary function for MacKenzie distributions
  !!
  !! Uses the expressions from A. Morawiec, J.Appl.Cryst. (1995) 28:289-293
 
IMPLICIT NONE

real(kind=dbl),INTENT(IN)       :: a 
real(kind=dbl),INTENT(IN)       :: b 
real(kind=dbl),INTENT(IN)       :: c 
real(kind=dbl)                  :: CC 

CC =  acos( (cos(c)-cos(a)*cos(b)) / (sin(a)*sin(b)) )

end function MKCC

!--------------------------------------------------------------------------
recursive function MKS2(a, b, c) result(S2)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! auxiliary function for MacKenzie distributions
  !!
  !! Uses the expressions from A. Morawiec, J.Appl.Cryst. (1995) 28:289-293

IMPLICIT NONE

real(kind=dbl),INTENT(IN)       :: a 
real(kind=dbl),INTENT(IN)       :: b 
real(kind=dbl),INTENT(IN)       :: c 
real(kind=dbl)                  :: S2 

S2 = 2.D0 * ( cPi - MKCC(a,b,c) - cos(a) * MKCC(c,a,b) - cos(b) * MKCC(b,c,a) )

end function MKS2

!--------------------------------------------------------------------------
recursive subroutine getMacKenzieDistribution_(self, Nmisor, misor, MK)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! computes the theoretical MacKenzie distribution for a given rotational point group
  !!
  !! Uses the expressions from A. Morawiec, J.Appl.Cryst. (1995) 28:289-293

use mod_io
use mod_symmetry

IMPLICIT NONE

class(so3_T),INTENT(INOUT)    :: self

integer(kind=irg),INTENT(IN)  :: Nmisor
real(kind=dbl),INTENT(IN)     :: misor(0:Nmisor)
real(kind=dbl),INTENT(INOUT)  :: MK(0:Nmisor)

type(IO_T)                    :: Message 
real(kind=dbl),allocatable    :: pomega(:), chi(:), rr(:)
real(kind=dbl)                :: h, h2, h3
integer(kind=irg)             :: prot, pgrotOrder, i
real(kind=dbl),parameter      :: nn(4) = dble((/ 2, 4, 3, 6 /)), nnn(4) = dble((/ 4, 8, 6, 12 /))

! allocate array
allocate(pomega(Nmisor))

! first get the distribution pomega for the no-symmetry case 
prot = PGrot(self%pgnum)
pgrotOrder = PGTHDorder(prot)
pomega = (1.D0 / (2.D0 * cPi * cPi)) * sin(misor*0.5D0)**2

! if there is rotational symmetry, then compute the solid angle function chi
if (prot.gt.1) then 
  allocate(chi(Nmisor), rr(Nmisor))
  chi = 4.D0 * cPi
  rr = tan(misor*0.5D0)
  if ((prot.eq.3).or.(prot.eq.6)) h = nn(1)
  if ((prot.eq.9).or.(prot.eq.12)) h = nn(2)
  if ((prot.eq.16).or.(prot.eq.18)) h = nn(3)
  if ((prot.eq.21).or.(prot.eq.24)) h = nn(4)

  select case(prot) 
    case (3, 9, 16, 21)   ! cyclic point groups 
      chi = 4.D0 * cPi
      h2 = tan(cPi * 0.5D0 / h)
      do i = 0, Nmisor
        if (rr(i).ge.h2) then 
          chi(i) = chi(i) - 4.D0 * cPi * (1.D0 - cos( acos(h2/rr(i)) ) )
        end if
      end do
      MK = (cPi/180.D0) * h * pomega * chi

    case (6, 12, 18, 24)   ! dihedral groups 
      h2 = tan(cPi * 0.5D0 / h)
      h3 = 1.D0
      do i = 0, Nmisor
        if (rr(i).gt.h2) then 
          chi(i) = chi(i) - 4.D0 * cPi * (1.D0 - cos( acos(h2/rr(i)) ) )
        end if
        if (rr(i).gt.h3) then 
          chi(i) = chi(i) - 4.D0 * h * cPi * (1.D0 - cos( acos(h3/rr(i)) ) )
        end if
        if (rr(i).gt.sqrt(h3+h2*h2)) then 
          chi(i) = chi(i) + 4.D0 * h * MKS2( acos(h2/rr(i)), acos(h3/rr(i)), cPi/2.D0) &
                   + 2.D0 * h * MKS2( acos(h3/rr(i)), acos(h3/rr(i)), cPi/h)           
        end if
        if (rr(i).gt.sqrt(h3+2.D0*h2*h2)) then 
          chi(i) = 0.D0
        end if
      end do
      MK = (cPi/180.D0) * (2.D0 * h) * pomega * chi
 
     case (28)   ! tetrahedral group 
      h2 = 1.D0/sqrt(2.D0)
      h3 = 1.D0/sqrt(3.D0)
      do i = 0, Nmisor
        if (rr(i).gt.h3) then 
          chi(i) = chi(i) - 16.D0 * cPi * (1.D0 - cos( acos(h3/rr(i)) ) )
        end if
        if (rr(i).gt.h2) then 
          chi(i) = chi(i) + 12.D0 * MKS2( acos(h3/rr(i)),  acos(h3/rr(i)), acos(h3*h3) ) 
        end if
        if (rr(i).gt.1.D0) then 
          chi(i) = 0.D0
        end if
      end do
      MK = (cPi/180.D0) * 12.D0 * pomega * chi

     case (30)   ! octahedral group 
      h2 = sqrt(2.D0) - 1.D0
      h3 = 1.D0/sqrt(3.D0)
      do i = 0, Nmisor
        if (rr(i).gt.h2) then 
          chi(i) = chi(i) - 12.D0 * cPi * (1.D0 - cos( acos(h2/rr(i)) ) )
        end if
        if (rr(i).gt.h3) then 
          chi(i) = chi(i) - 16.D0 * cPi * (1.D0 - cos( acos(h3/rr(i)) ) )
        end if
        if (rr(i).gt.(2.D0-sqrt(2.D0))) then 
          chi(i) = chi(i) + 12.D0 * MKS2( acos(h2/rr(i)),  acos(h2/rr(i)), cPi*0.5D0 ) + &
                   24.D0 * MKS2( acos(h2/rr(i)),  acos(h3/rr(i)), acos(h3) ) 
        end if
        if (rr(i).gt.sqrt(23.0-16.D0*sqrt(2.D0))) then 
          chi(i) = 0.D0
        end if
      end do
      MK = (cPi/180.D0) * 24.D0 * pomega * chi

    case default
      call Message%printError('getMacKenzieDistribution',' non-existent rotational point group number')
  end select
  deallocate(chi, rr)
! force the last point to zero  
  MK(Nmisor) = 0.D0
else 
  MK = (cPi/180.D0) * (4.D0 * cPi) * pomega
  deallocate(pomega)
end if

end subroutine getMacKenzieDistribution_

!--------------------------------------------------------------------------
recursive subroutine ReduceDisorientationtoMFZ_(self, ro, SG, roMFZ)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/23/20
  !!
  !! Reduce a disorientation (Rodrigues) to the Mackenzie Fundamental Zone
  !!
  !! This requires the symmetry operations of the Rodrigues FZ, which 
  !! includes mirrors and inversion symmetry, i.e., the regular (full) point group
  !! symmetry of the shape of the RFZ.  We already have those implemented in the 
  !! regular symmetry module, and we assume that those symmetry matrices have been
  !! initialized and are present in the cell class.  Then we just call the regular
  !! CalcStar routine to generate the equivalents and pick the one that is inside
  !! the MFZ.

use mod_quaternions
use mod_symmetry

IMPLICIT NONE

class(so3_T),INTENT(INOUT)        :: self
type(r_T), INTENT(INOUT)          :: ro
type(SpaceGroup_T), INTENT(INOUT) :: SG
type(r_T), INTENT(OUT)            :: roMFZ

type(r_T)                         :: rod
real(kind=dbl)                    :: r(3), x(4), mag
integer(kind=irg)                 :: i, j, n
real(kind=dbl),allocatable        :: stmp(:,:)
logical                           :: inMFZ

roMFZ = r_T( rdinp = (/ 0.D0, 0.D0, 1.D0, 0.D0 /) )
x = ro%r_copyd()
r(1:3) = x(1:3)*x(4)
call SG%CalcStar(r,n,stmp,'d')

MFZloop: do j=1,n
  mag = dsqrt(sum(stmp(j,1:3)**2))
  rod = r_T( rdinp = (/ stmp(j,1:3)/mag, mag /) )
  inMFZ = self%IsinsideMFZ(rod)
  if (inMFZ) then
   roMFZ = rod
   EXIT MFZloop
  end if
! we really should never get to the following line ...
  if (j.eq.n) then
    roMFZ = r_T( rdinp = (/ 0.D0, 0.D0, 1.D0, 0.D0 /) )
  end if
end do MFZloop

end subroutine ReduceDisorientationtoMFZ_

!--------------------------------------------------------------------------
recursive subroutine ReduceOrientationtoCubicEFZ_(self, eu, euFZ)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/23/20
  !!
  !! Reduce an orientation (Euler angles) to the cubic FZ in Euler space
  !!
  !! [no longer uses dict module]

use mod_quaternions
use mod_symmetry
use mod_io 

IMPLICIT NONE

class(so3_T),INTENT(INOUT)    :: self
type(e_T), INTENT(INOUT)      :: eu 
type(e_T), INTENT(OUT)        :: euFZ

type(IO_T)                    :: Message
type(quaternion_T)            :: Mu, qu, qS 
type(q_T)                     :: qq
real(kind=dbl)                :: c, s, z, x(4), y(3)
integer(kind=irg)             :: i, j, Pmdims

euFZ = e_T( edinp = (/ 0.D0, 0.D0, 0.D0 /) )
Pmdims = 24  ! cubic symmetry has rotational group 432

qq = eu%eq()
x = qq%q_copyd()
Mu = quaternion_T( qd = x )
call Mu%quat_pos()

FZloop: do j=1,Pmdims
  qS = quaternion_T( qd = SYM_Qsymop(1:4,j) )
  qu = qS*Mu
  call qu%quat_pos()
  qq = q_T( qdinp = qu%get_quatd() )
  euFZ = qq%qe()
  y = euFZ%e_copyd()
  ! apply the cubic Euler FZ boundary conditions
  c = cos(y(3))
  s = sin(y(3))
  z = acos(minval( (/ c/sqrt(1+c*c), s/sqrt(1+s*s) /) ))
  if ((y(2).gt.z).and.(y(2).lt.cPi/2D0).and.(y(3).lt.cPi/2.D0)) EXIT FZloop
! we really should never get to the following line ...
  if (j.eq.Pmdims) call Message%printWarning('ReduceOrientationtoCubicEFZ: no orientation found in cubic EFZ')
end do FZloop

end subroutine ReduceOrientationtoCubicEFZ_

!--------------------------------------------------------------------------
recursive subroutine ReducelisttoRFZ_(self, Pm)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/25/20
  !!
  !! takes a linked list l and reduces the entire list to the selected Rodrigues FZ

use mod_quaternions 
use mod_rotations

IMPLICIT NONE 

class(so3_T), INTENT(INOUT)             :: self 
type(QuaternionArray_T),INTENT(INOUT)   :: Pm

integer(kind=irg)                       :: i 
type(QuaternionArray_T)                 :: qAR
type(Quaternion_T)                      :: qq
type(q_T)                               :: qu
type(r_T)                               :: roFZ

! first, convert the linked list to a QuaternionArray_T object 
call self%listtoQuaternionArray( qAR, 'FZ' )
! then reduce the orientations to the RFZ 
do i = 1, qAR%getQnumber()
  qq = qAR%getQuatfromArray(i)
  qu = q_T( qdinp = qq%get_quatd() )
  call self%ReduceOrientationtoRFZ_( qu, Pm, roFZ )
  qu = roFZ%rq()
  call qAR%insertQuatinArray(i, Quaternion_T( qd = qu%q_copyd() ) )
end do
! and convert the array back to the linked list 
call self%QuaternionArraytolist( qAR, 'FZ')

end subroutine ReducelisttoRFZ_

!--------------------------------------------------------------------------
recursive subroutine ReduceOrientationtoRFZ_(self, rot, Pm, roFZ, MFZ)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/23/20
  !!
  !! Reduce an orientation (any of nine ... ) to the Rodrigues Fundamental Zone
  !!
  !! [no longer uses dict module but requires that the calling program 
  !! initializes the Pm QuaternionArray_T class for the correct point group
  !! using the QSym_Init method in mod_quaternions]

use mod_quaternions
use mod_symmetry
use mod_io 
use mod_math

IMPLICIT NONE

class(so3_T),INTENT(INOUT)             :: self
class(*), INTENT(INOUT)                :: rot 
type(QuaternionArray_T),INTENT(INOUT)  :: Pm 
type(r_T), INTENT(OUT)                 :: roFZ
logical,OPTIONAL,INTENT(IN)            :: MFZ

type(IO_T)                             :: Message
type(Quaternion_T)                     :: Mu, qu, qS, pp 
type(q_T)                              :: qq
type(r_T)                              :: rod      
real(kind=dbl)                         :: x(4), y(3)
integer(kind=irg)                      :: i, j, Pmdims
logical                                :: useMFZ
real(kind=dbl)                         :: tol

tol = 1.0D+5

useMFZ = .FALSE.
if (present(MFZ)) then 
  if (MFZ.eqv..TRUE.) useMFZ = .TRUE.
endif

roFZ = r_T( rdinp = (/ 0.D0, 0.D0, 1.D0, 0.D0 /) )
Pmdims = Pm%getQnumber()

! convert the input rot to quaternion form 
select type (rot)
  class is (e_T)
    qq = rot%eq()
  class is (o_T)
    qq = rot%oq()
  class is (a_T)
    qq = rot%aq()
  class is (h_T)
    qq = rot%hq()
  class is (s_T)
    qq = rot%sq()
  class is (c_T)
    qq = rot%cq()
  class is (q_T)
    qq = rot
  class is (v_T)
    qq = rot%vq()
  class is (r_T)
    qq = rot%rq()
  class default
end select

Mu = Quaternion_T( qd = qq%q_copyd() )
call Mu%quat_pos()

FZloop: do j=1,Pmdims
  qu = Pm%getQuatfromArray(j)*Mu
  x = qu%get_quatd()
  if (x(1).lt.0D0) x = -x 
  qq = q_T( qdinp = x )
  rod = qq%qr()
  x = rod%r_copyd()
  if(abs(x(4)) .gt. tol) rod = r_T( rdinp = (/ x(1:3), inftyd() /) )
  
  if (useMFZ.eqv..TRUE.) then
    if (self%IsinsideMFZ(rod)) EXIT FZloop
  else
    if (self%IsinsideFZ(rod))  EXIT FZloop
  end if
  ! we really should never get to the following line ...
  if (j.eq.Pmdims) call Message%printWarning( 'ReduceOrientationtoRFZ: no solution found')
end do FZloop

roFZ = rod

end subroutine ReduceOrientationtoRFZ_

!--------------------------------------------------------------------------
recursive subroutine getDisorientation_(self, Pm, or1, or2, disax)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/23/20
  !!
  !! Determine the disorientation angle between two orientations (in radians)
  !!
  !! This routine is now generalized to take arbitrary orientation classes as input
  !! (they need not have the same class!), and generates the complete disorientation as an
  !! axis-angle pair class. The calling program can then chose which part to use, the
  !! disorientation angle or the axis (or both).

use mod_quaternions
use mod_symmetry
use mod_rotations
use mod_math

IMPLICIT NONE

class(so3_T),INTENT(INOUT)             :: self
type(QuaternionArray_T),INTENT(INOUT)  :: Pm 
class(*), INTENT(INOUT)                :: or1 
class(*), INTENT(INOUT)                :: or2 
type(a_T), INTENT(OUT)                 :: disax

type(Quaternion_T)                     :: qu1, qu2, Mu, qu, Mus, qus, p 
real(kind=dbl)                         :: a, ac, ax(3), x(4)
integer(kind=irg)                      :: j, k, Pmdims

Pmdims = Pm%getQnumber() 

! first we need to convert the input orientations into Quaternion_T objects;
! since we do not know ahead of time which representations will be provided
! as input, we use a select type block in the getQfromClass method in 
! the rotations module.
qu1 = getQfromClass(or1)
qu2 = getQfromClass(or2)

Mu = Quaternion_T( qd = qu1%get_quatd() )
call Mu%quat_pos()

qu = Quaternion_T( qd = qu2%get_quatd() )
call qu%quat_pos()

ac = 1000.D0
do j=1,Pmdims ! loop over the symmetric equivalents of Mu
  Mus = Pm%getQuatfromArray(j)*Mu
  call Mus%quat_pos()
  do k=1,Pmdims
    qus = Pm%getQuatfromArray(k)*qu
    call qus%quat_pos()
    p = Mus*conjg(qus)
    call p%quat_pos()
    x = p%get_quatd()
    a = 2.0*acos(x(1))
    if (a.lt.ac) then
      ac = a
      ax(1:3) = x(2:4)/vecnorm(x(2:4))
    end if
    p = qus*conjg(Mus)
    call p%quat_pos()
    x = p%get_quatd()
    a = 2.0*acos(x(1))
    if (a.lt.ac) then
      ac = a
      ax(1:3) = x(2:4)/vecnorm(x(2:4))
    end if
  end do
end do 

disax = a_T( adinp = (/ ax(1:3), ac /) )

end subroutine getDisorientation_

!--------------------------------------------------------------------------
recursive subroutine getDisorientationTwoPhases_(self, Pm1, Pm2, or1, or2, disax)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/23/20
  !!
  !! Determine the disorientation angle between two orientations (in radians) for 
  !! two different phases.
  !!
  !! This routine is now generalized to take arbitrary orientation classes as input
  !! (they need not have the same class!), and generates the complete disorientation as an
  !! axis-angle pair class. The calling program can then chose which part to use, the
  !! disorientation angle or the axis (or both).

use mod_quaternions
use mod_symmetry
use mod_rotations
use mod_math

IMPLICIT NONE

class(so3_T),INTENT(INOUT)              :: self
type(QuaternionArray_T),INTENT(INOUT)   :: Pm1
type(QuaternionArray_T),INTENT(INOUT)   :: Pm2
class(*), INTENT(INOUT)                 :: or1 
class(*), INTENT(INOUT)                 :: or2 
type(a_T), INTENT(OUT)                  :: disax

type(Quaternion_T)                      :: qu1, qu2, Mu, qu, Mus, qus, p 
real(kind=dbl)                          :: a, ac, ax(3), x(4)
integer(kind=irg)                       :: j, k, Pm1dims, Pm2dims

Pm1dims = Pm1%getQnumber() 
Pm2dims = Pm2%getQnumber() 

! first we need to convert the input orientations into Quaternion_T objects;
! since we do not know ahead of time which representations will be provided
! as input, we use a select type block in the getQfromClass method in 
! the rotations module.
qu1 = getQfromClass(or1)
qu2 = getQfromClass(or2)

Mu = Quaternion_T( qd = qu1%get_quatd() )
call Mu%quat_pos()

qu = Quaternion_T( qd = qu2%get_quatd() )
call qu%quat_pos()

ac = 1000.D0
do j=1,Pm1dims ! loop over the symmetric equivalents of Mu
  Mus = Pm1%getQuatfromArray(j)*Mu
  call Mus%quat_pos()
  do k=1,Pm2dims
    qus = Pm2%getQuatfromArray(k)*qu
    call qus%quat_pos()
    p = Mus*conjg(qus)
    call p%quat_pos()
    x = p%get_quatd()
    a = 2.0*acos(x(1))
    if (a.lt.ac) then
      ac = a
      ax(1:3) = x(2:4)/vecnorm(x(2:4))
    end if
  end do
end do 

disax = a_T( adinp = (/ ax(1:3), ac /) )

end subroutine getDisorientationTwoPhases_

!--------------------------------------------------------------------------
recursive subroutine getAverageDisorientationMap_(self, quats, Pm, wd, ht, ADMap) 
  !! author: MDG
  !! version: 1.0 
  !! date: 01/23/20
  !!
  !! Determine the average disorientation map (in degrees)
  !!
  !! the calling program must initialize the symmetry elements in Pm, as well 
  !! as generate the input quaternion array in a QuaternionArray_T class.

use mod_quaternions
use mod_rotations

IMPLICIT NONE

class(so3_T),INTENT(INOUT)              :: self
type(QuaternionArray_T),INTENT(IN)      :: quats
 !! orientation array
type(QuaternionArray_T),INTENT(INOUT)   :: Pm
 !! symmetry operators in quaternion form
integer(kind=irg),INTENT(IN)            :: wd
 !! map width (in pixels)
integer(kind=irg),INTENT(IN)            :: ht
 !! map height
real(kind=dbl),INTENT(OUT)              :: ADMap(wd,ht)
 !! output average disorientation angle map

integer(kind=irg)                       :: i, j, ic, icr, ict
real(kind=dbl)                          :: misor(4,wd,ht), denom(wd, ht), x(4) 
type(a_T)                               :: disax
type(q_T)                               :: q1, q2 
type(Quaternion_T)                      :: qu1, qu2 

ADMap = 0.D0
denom = 4.D0
misor = 0.D0   ! contains four misorientation angles in the order r, t, l, b 

! correct the multiplicities along the edges
denom(2:wd-1,1) = 3.D0
denom(2:wd-1,ht) = 3.D0
denom(1,2:ht-1) = 3.D0
denom(wd,2:ht-1) = 3.D0
! correct the multiplicities in the corners
denom(1,1) = 2.D0
denom(1,ht) = 2.D0
denom(wd,1) = 2.D0
denom(wd,ht) = 2.D0

! we'll do this line by line (horizontally)
do j=1,ht
  do i=1,wd
    ic = wd*(j-1)+i
    icr = wd*(j-1)+i+1
    ict = wd*j + i

! right neighbor (also includes left one)
    if (i.lt.wd) then 
      qu1 = quats%getQuatfromArray(ic)
      q1 = q_T( qdinp = qu1%get_quatd() )
      qu2 = quats%getQuatfromArray(icr)
      q2 = q_T( qdinp = qu2%get_quatd() )
      call self%getDisorientation(Pm, q1, q2, disax)
      x = disax%a_copyd()
      misor(1,i,j) = x(4)
      misor(3,i+1,j) = x(4)
    end if

! top neighbor
    if (j.lt.ht) then
      qu1 = quats%getQuatfromArray(ic)
      q1 = q_T( qdinp = qu1%get_quatd() )
      qu2 = quats%getQuatfromArray(ict)
      q2 = q_T( qdinp = qu2%get_quatd() )
      call self%getDisorientation(Pm, q1, q2, disax)
      x = disax%a_copyd()
      misor(2,i,j) = x(4)
      misor(4,i,j+1) = x(4)
    end if
  end do 
end do 

! then take the average
ADMap = sum(misor,1)/denom

! and convert to degrees
ADMap = ADMap / dtor

end subroutine getAverageDisorientationMap_

end module mod_so3
