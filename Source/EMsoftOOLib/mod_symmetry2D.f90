! ###################################################################
! Copyright (c) 2014-2024, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_symmetry2D
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/10/20
  !!
  !! all 2D symmetry-related constants, derived types and methods

use mod_kinds
use mod_global

IMPLICIT NONE


! for many diffraction calculations we need the 2D planar point groups; 
! the maximum order of such a group is 12, and there are only 10 of them, with
! two settings for one of them (3m1 and 31m).
type symdata2D
  integer(kind=irg)     :: SYM_pgnum                    !< 2D point group number
  integer(kind=irg)     :: SYM_MATnum                   !< number of non-zero symmetry matrices (order)
  integer(kind=irg)     :: SYM_direc(12,2,2)            !< point group matrices (filled in by Generate2DSymmetry)
end type symdata2D

contains 

!--------------------------------------------------------------------------
recursive subroutine Generate2DSymmetry(TDPG,pgn)
!DEC$ ATTRIBUTES DLLEXPORT :: Generate2DSymmetry
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/10/20
  !!
  !! generate the symmetry matrices for one of the 2D planar point groups.
  !!
  !! The output of this routine is stored in the TDPG structure that 
  !! is defined at the end of the symmetryvars module.  Routine verified on 10/2/13.

use mod_io

IMPLICIT NONE

type (symdata2D),INTENT(INOUT)  :: TDPG
!f2py intent(in,out) ::  TDPG
integer(kind=irg),INTENT(IN)    :: pgn          !< point group number
type(IO_T)                      :: Message 

! here we define all 8 possible 2x2 matrices (remember column major order !)
integer(kind=irg),parameter     :: mI(2,2) = reshape( (/ 1, 0, 0, 1 /), (/2,2/))        ! identity
integer(kind=irg),parameter     :: mA(2,2) = reshape( (/-1, 0, 0, 1 /), (/2,2/))        ! mirror
integer(kind=irg),parameter     :: mB(2,2) = reshape( (/ 0, 1,-1, 0 /), (/2,2/))        ! 4-fold
integer(kind=irg),parameter     :: mC(2,2) = reshape( (/ 0, 1, 1, 0 /), (/2,2/))        ! diagonal mirror
integer(kind=irg),parameter     :: mD(2,2) = reshape( (/-1, 1,-1, 0 /), (/2,2/))        ! 3-fold
integer(kind=irg),parameter     :: mE(2,2) = reshape( (/ 0,-1, 1,-1 /), (/2,2/))        ! 3-fold twice
integer(kind=irg),parameter     :: mF(2,2) = reshape( (/-1, 0,-1, 1 /), (/2,2/))        ! 6-fold
integer(kind=irg),parameter     :: mG(2,2) = reshape( (/ 1,-1, 0,-1 /), (/2,2/))        ! 6-fold twice


write (*,*) 'Generate2DSymmetry : ', pgn 

TDPG%SYM_pgnum = pgn

select case (pgn)
  case (1) ! 1
        TDPG%SYM_MATnum = 1
        TDPG%SYM_direc(1,1:2,1:2) = mI
!------------
  case (2) ! 2
        TDPG%SYM_MATnum = 2
        TDPG%SYM_direc(1,1:2,1:2) = mI
        TDPG%SYM_direc(2,1:2,1:2) = -mI
!------------
  case (3) ! m
        TDPG%SYM_MATnum = 2
        TDPG%SYM_direc(1,1:2,1:2) = mI
        TDPG%SYM_direc(2,1:2,1:2) = mA
!------------
  case (4) ! 2mm
        TDPG%SYM_MATnum = 4
        TDPG%SYM_direc(1,1:2,1:2) = mI
        TDPG%SYM_direc(2,1:2,1:2) = -mI
        TDPG%SYM_direc(3,1:2,1:2) = mA
        TDPG%SYM_direc(4,1:2,1:2) = -mA
!------------
  case (5) ! 4
        TDPG%SYM_MATnum = 4
        TDPG%SYM_direc(1,1:2,1:2) = mI
        TDPG%SYM_direc(2,1:2,1:2) = -mI
        TDPG%SYM_direc(3,1:2,1:2) = mB
        TDPG%SYM_direc(4,1:2,1:2) = -mB
!------------
  case (6) ! 4mm
        TDPG%SYM_MATnum = 8
        TDPG%SYM_direc(1,1:2,1:2) = mI
        TDPG%SYM_direc(2,1:2,1:2) = -mI
        TDPG%SYM_direc(3,1:2,1:2) = mA
        TDPG%SYM_direc(4,1:2,1:2) = -mA
        TDPG%SYM_direc(5,1:2,1:2) = mB
        TDPG%SYM_direc(6,1:2,1:2) = -mB
        TDPG%SYM_direc(7,1:2,1:2) = mC
        TDPG%SYM_direc(8,1:2,1:2) = -mC
!------------
  case (7) ! 3
        TDPG%SYM_MATnum = 3
        TDPG%SYM_direc(1,1:2,1:2) = mI
        TDPG%SYM_direc(2,1:2,1:2) = mD
        TDPG%SYM_direc(3,1:2,1:2) = mE
!------------
  case (8) ! 3m1
        TDPG%SYM_MATnum = 6
        TDPG%SYM_direc(1,1:2,1:2) = mI
        TDPG%SYM_direc(2,1:2,1:2) = mD
        TDPG%SYM_direc(3,1:2,1:2) = mE
        TDPG%SYM_direc(4,1:2,1:2) = -mC
        TDPG%SYM_direc(5,1:2,1:2) = -mF
        TDPG%SYM_direc(6,1:2,1:2) = -mG
!------------
  case (9) ! 6
        TDPG%SYM_MATnum = 6
        TDPG%SYM_direc(1,1:2,1:2) = mI
        TDPG%SYM_direc(2,1:2,1:2) = -mI
        TDPG%SYM_direc(3,1:2,1:2) = mD
        TDPG%SYM_direc(4,1:2,1:2) = -mD
        TDPG%SYM_direc(5,1:2,1:2) = mE
        TDPG%SYM_direc(6,1:2,1:2) = -mE
!------------
  case (10) ! 6mm
        TDPG%SYM_MATnum = 12
        TDPG%SYM_direc(1,1:2,1:2) = mI
        TDPG%SYM_direc(2,1:2,1:2) = -mI
        TDPG%SYM_direc(3,1:2,1:2) = mD
        TDPG%SYM_direc(4,1:2,1:2) = -mD
        TDPG%SYM_direc(5,1:2,1:2) = mE
        TDPG%SYM_direc(6,1:2,1:2) = -mE
        TDPG%SYM_direc(7,1:2,1:2) = mF
        TDPG%SYM_direc(8,1:2,1:2) = -mF
        TDPG%SYM_direc(9,1:2,1:2) = mG
        TDPG%SYM_direc(10,1:2,1:2) = -mG
        TDPG%SYM_direc(11,1:2,1:2) = mC
        TDPG%SYM_direc(12,1:2,1:2) = -mC        
!------------
  case (11) ! 31m
        TDPG%SYM_MATnum = 6
        TDPG%SYM_direc(1,1:2,1:2) = mI
        TDPG%SYM_direc(2,1:2,1:2) = mC
        TDPG%SYM_direc(3,1:2,1:2) = mD
        TDPG%SYM_direc(4,1:2,1:2) = mE
        TDPG%SYM_direc(5,1:2,1:2) = mF
        TDPG%SYM_direc(6,1:2,1:2) = mG
!------------
  case default
        call Message%printError('Generate2DSymmetry',' unknown 2D point group number')
        stop
end select

end subroutine Generate2DSymmetry

!--------------------------------------------------------------------------
recursive subroutine Apply2DPGSymmetry(TDPG,ipx,ipy,isym,iequiv,nequiv)
!DEC$ ATTRIBUTES DLLEXPORT :: Apply2DPGSymmetry
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/07/20
  !!
  !! Apply the 2D point group symmetry to a pair of coordinates 
  !!
  !! This routine returns a set of unique equivalent coordinates for the
  !! input point.  Note that for special points, the number nequiv must be reduced
  !! from its point group order value.

use mod_io

IMPLICIT NONE

type(symdata2D),INTENT(INOUT)      :: TDPG
!f2py intent(in,out) ::  TDPG   
integer(kind=irg),INTENT(IN)       :: ipx
integer(kind=irg),INTENT(IN)       :: ipy
integer(kind=irg),INTENT(IN)       :: isym
integer(kind=irg),INTENT(OUT)      :: iequiv(2,12)
integer(kind=irg),INTENT(OUT)      :: nequiv
    
integer(kind=irg)                  :: i, j, pequiv(2,12), mequiv
real(kind=sgl),parameter           :: eps = 1.0E-6 
real(kind=sgl)                     :: diff
logical                            :: newp

! make sure that the symmetry matrices have been predefined; if not, then
! compute them first
if (TDPG%SYM_pgnum.ne.isym) call Generate2DSymmetry(TDPG,isym)

! set the order;  note that this may need to reduced for special points
mequiv = TDPG%SYM_MATnum
iequiv = 0
pequiv = 0

! compute the transformed coordinates
do i=1,mequiv
  pequiv(1:2,i) = matmul( (/ ipx, ipy /), TDPG%SYM_direc(i,1:2,1:2) )
end do

! the first point is always unique, so simply copy it
nequiv = 1
iequiv(1:2,1) = pequiv(1:2,1)

! next, identify double entries and remove them from the list
do i=2,mequiv
  newp = .TRUE.
  do j=1,nequiv
   diff= sum(abs(pequiv(1:2,i)-iequiv(1:2,j)))
   if (diff.le.eps) then
     newp  = .FALSE.
   endif
  end do

! yes, it is a new point
  if (newp) then
   nequiv=nequiv+1
   iequiv(1:2,nequiv) = pequiv(1:2,i)
  end if
end do

! that's it.
end subroutine Apply2DPGSymmetry



end module mod_symmetry2D
