! ###################################################################
! Copyright (c) 2013-2022, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_QCcrystallography
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/30/22
  !!
  !! class definition for the 6D quasi-crystal cell

use mod_kinds
use mod_global

IMPLICIT NONE 

! class definition
type, public :: QCcell_T
private 

  integer(kind=irg)                     :: atno
  integer(kind=irg)                     :: imax
  integer(kind=irg)                     :: numindices
  integer(kind=irg),allocatable         :: facts(:,:)
  integer(kind=irg),allocatable         :: Ucgindex(:)
  logical,allocatable                   :: Ucgcalc(:)
  integer(kind=irg),allocatable         :: inverseIndex(:,:)
  real(kind=dbl)                        :: epvec(3,6), epar(6,3)
  real(kind=dbl)                        :: eovec(3,6), eperp(6,3)
  real(kind=dbl)                        :: Mp(6,6), Picos(6,6)
  real(kind=dbl)                        :: Mo(6,6), Qicos(6,6)
  real(kind=dbl)                        :: dsm(6,6), rsm(6,6)
  real(kind=dbl)                        :: dmt(6,6), rmt(6,6)
  real(kind=dbl)                        :: scaling(6,6)
  real(kind=dbl)                        :: SYM_icos(6,6,120)      ! 532 rotational group in matrix representation
  real(kind=dbl)                        :: QClatparm, alphaij, alphastarij
  real(kind=dbl)                        :: dmin
  real(kind=dbl)                        :: vol
  real(kind=dbl)                        :: gmax_orth
  real(kind=dbl)                        :: DWF
  real(kind=dbl)                        :: voltage
  real(kind=dbl)                        :: mRelCor
  real(kind=dbl)                        :: mSigma
  real(kind=dbl)                        :: mPsihat
  real(kind=dbl)                        :: mLambda
  real(kind=dbl)                        :: Upzero
  real(kind=dbl)                        :: xizerop
  real(kind=dbl)                        :: multiplicity
  character(1)                          :: centering   ! 'P','I','F'
  complex(kind=dbl),allocatable         :: LUT(:)
  complex(kind=dbl),allocatable         :: LUTqg(:)
  logical, allocatable                  :: dbdiff(:)
  character(fnlen)                      :: SGname(11), QCtype, fname
  integer(kind=irg)                     :: ATOM_ntype, ATOM_type(maxpasym), numat(maxpasym)
  real(kind=sgl),allocatable            :: apos(:,:,:)
  real(kind=sgl)                        :: ATOM_pos(maxpasym,10)

contains
private 
  procedure, pass(self) :: get6Dindex_
  procedure, pass(self) :: invert6Dindex_
  procedure, pass(self) :: GetQCLatParm_
  procedure, pass(self) :: 
  procedure, pass(self) :: 

  generic, public :: get6Dindex => get6Dindex_
  generic, public :: invert6Dindex => invert6Dindex_
  generic, public :: GetQCLatParm => GetQCLatParm_
  generic, public :: 
  generic, public :: 

end type QCcell_T

! the constructor routine for this class 
interface QCcell_T
  module procedure QCcell_constructor
end interface QCcell_T

contains

!--------------------------------------------------------------------------
type(QCcell_T) function QCcell_constructor( ) result(QCcell)
!! author: MDG 
!! version: 1.0 
!! date: 01/30/22
!!
!! constructor for the QCcell_T Class
 
IMPLICIT NONE

allocate( self%inversindex(nLUT, 6) )




end function QCcell_constructor

!--------------------------------------------------------------------------
subroutine QCcell_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 01/30/22
!!
!! destructor for the QCcell_T Class
 
IMPLICIT NONE

type(QCcrystallography_T), INTENT(INOUT)  :: self 

call reportDestructor('QCcell_T')

end subroutine QCcell_destructor

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
recursive function get6Dindex_(self, QCindex ) result(gindex)
!DEC$ ATTRIBUTES DLLEXPORT :: get6Dindex_
!! author: MDG/SS 
!! version: 1.0 
!! date: 01/30/22
!!
!! Convert a 6-component Miller index to a single lookup index
 
IMPLICIT NONE

class(QCcell_T),INTENT(INOUT)     :: self
integer(kind=irg),INTENT(IN)      :: QCindex(6)
integer(kind=irg)                 :: gindex, imax, isize, g(6)

imax   = self%imax
isize  = 2 * imax + 1
g      = QCindex + imax 
gindex = g(1)*isize**5 + g(2)*isize**4 + g(3)*isize**3 + g(4)*isize**2 + g(5)*isize + g(6) + 1

end function get6Dindex_

!--------------------------------------------------------------------------
recursive function invert6Dindex_(self, gindex ) result(QCindex)
!DEC$ ATTRIBUTES DLLEXPORT :: invert6Dindex_
!! author: MDG/SS 
!! version: 1.0 
!! date: 01/30/22
!!
!! Convert a single lookup index into the corresponding 6-component Miller index 
 
IMPLICIT NONE

class(QCcell_T),INTENT(INOUT)     :: self
integer(kind=irg),INTENT(IN)      :: gindex
integer(kind=irg)                 :: QCindex(6)

QCindex(1:6) = self%inverseIndex(gindex,1:6)

end function invert6Dindex_

!--------------------------------------------------------------------------
recursive subroutine Get3DQCLatParm_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: Get3DQCLatParm_
!! author: MDG/SS 
!! version: 1.0 
!! date: 01/30/22
!!
!! input of lattice parameters for 3D (icosahedral) quasi-crystal 

use mod_io

IMPLICIT NONE

class(QCcell_T)          :: self

type(IO_T)               :: Message
real(kind=dbl)           :: io_real(1)   !< double precision real input array
integer(kind=irg)        :: std

 self%QCtype         = 'Ico'
 self%alphaij        =  90.D0
 self%alphastarij    =  90.D0

 call Message(' Using the higher dimensional cut-and-project approach:', frm = "(A)")
 call Message(' -------------------------------', frm = "(A)")
 
  ! get the lattice parameters
 call Message('Enter lattice parameters', frm = "(//A)")

  ! in-plane lattice parameters
 call ReadValue('    a_i | i = {1,2,3,4,5,6} (hyper-cube) [nm] = ', io_real, 1)
 self%QClatparm = io_real(1)
  
end subroutine Get3DQCLatParm_



end module mod_QCcrystallography