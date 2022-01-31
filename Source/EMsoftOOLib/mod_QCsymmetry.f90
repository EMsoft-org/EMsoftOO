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

module mod_QCsymmetry
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/30/22
  !!
  !! class definition for 2D (axial) and 3D (icosahedral) quasi-crystal symmetry
  !!
  !! For the 3D QCs, we use an approach similar to that of the regular space groups
  !! with generator strings; for the 2D axial QCs, that has not been implemented yet 
  !! and the approach is more direct in terms of the definition of the generator matrices.

use mod_kinds
use mod_global

IMPLICIT NONE 

private 

!>  SGname all icosahedral space group names
! ICOSAHEDRAL SPACE GROUPS
  character(13), public, dimension(11) :: SGname_ico = (/  &
          " P 5 3 2     ", " P 5_1 3 2   ", " I 5 3 2     ", " I 5_1 3 2   ", &
          " F 5 3 2     ", " F 5_1 3 2   ", " P -5 -3 2/m ", " P -5 -3 2/q ", &
          " I -5 -3 2/m ", " F -5 -3 2/m ", " F -5 -3 2/q " /)

!>  GL encoded generator strings
  character(57), public, dimension(11) :: GL_ico = (/  &
        "2b000000c000000                                          ", &
        "2bA00000c000000                                          ", &
        "3aBBBBBBb000000c000000                                   ", &
        "3aBBBBBBbA00000c000000                                   ", &
        "7aBB0000a0BB000a00BB00a000BB0a0000BBb000000c000000       ", &
        "7aBB0000a0BB000a00BB00a000BB0a0000BBbA00000c000000       ", &
        "3b000000c000000d000000                                   ", &
        "3b000000c000000d000000                                   ", &
        "4aBBBBBBb000000c000000d000000                            ", &
        "8aBB0000a0BB000a00BB00a000BB0a0000BBb000000c000000d000000", &
        "3b000000c000000d000000                                   " /)


! class definition
type, public :: QCspacegroup_T
private 
  integer(kind=irg)           :: nsg 
  integer(kind=irg)           :: SGnum
  integer(kind=irg)           :: GENnum             !< number of generator matrices
  integer(kind=irg)           :: MATnum             !< number of non-zero symmetry matrices
  integer(kind=irg)           :: NUMpt              !< number of point group operators
  integer(kind=irg)           :: N_Axial
  logical                     :: reduce             !< switch to enable/disable reduction to fundamental cell
  real(kind=dbl),allocatable  :: data(:,:,:)        !< all symmetry matrices for a given spacegroup
  real(kind=dbl),allocatable  :: direc(:,:,:)       !< direct space point group matrices
  real(kind=dbl),allocatable  :: recip(:,:,:)       !< reciprocal space point group matrices
  real(kind=dbl),allocatable  :: c(:,:)             !< dummy 6x6 matrix used for various computations
  character(11)               :: name
  character(fnlen),allocatable:: SGname(:)
  character(3)                :: QCtype
  character(57),allocatable   :: GL(:)

contains
private 
  procedure, pass(self) :: printSGtable_
  procedure, pass(self) :: GetQCSpaceGroup_
  procedure, pass(self) :: setnsg_
  procedure, pass(self) :: setSGnum_
  procedure, pass(self) :: setQCtype_
  procedure, pass(self) :: getnsg_
  procedure, pass(self) :: getSGnum_
  procedure, pass(self) :: getQCtype_
  procedure, pass(self) :: getAxialGroupNames_
  ! procedure, pass(self) :: 

  generic, public :: printSGtable => printSGtable_
  generic, public :: GetQCSpaceGroup => GetQCSpaceGroup_
  generic, public :: setnsg => setnsg_
  generic, public :: setSGnum => setSGnum_
  generic, public :: setQCtype => setQCtype_
  generic, public :: getnsg => getnsg_
  generic, public :: getSGnum => getSGnum_
  generic, public :: getQCtype => getQCtype_
  ! generic, public :: 

end type QCspacegroup_T

! the constructor routine for this class 
interface QCspacegroup_T
  module procedure QCspacegroup_constructor
end interface QCspacegroup_T

contains

!--------------------------------------------------------------------------
type(QCspacegroup_T) function QCspacegroup_constructor( nD, QCtype ) result(QCSG)
!DEC$ ATTRIBUTES DLLEXPORT :: QCspacegroup_constructor
!! author: MDG 
!! version: 1.0 
!! date: 01/30/22
!!
!! constructor for the QCspacegroup_T Class
 
IMPLICIT NONE

integer(kind=irg), INTENT(IN)   :: nD
character(3),INTENT(IN),OPTIONAL:: QCtype 

if (nD.eq.3) then 
! get the number of space groups in the SGname array
  QCSG%nsg = 11
  allocate( QCSG%data(5000,7,7)  )         !< all symmetry matrices for a given spacegroup
  allocate( QCSG%direc(5000,6,6) )         !< direct space point group matrices
  allocate( QCSG%recip(5000,6,6) )         !< reciprocal space point group matrices
  allocate( QCSG%c(7,7) )                  !< dummy matrix used for various computations
  allocate( QCSG%GL(QCSG%nsg), QCSG%SGname(QCSG%nsg) )
  QCSG%GL = GL_ico
  QCSG%SGname = SGname_ico
end if 

if (nD.eq.2) then ! axial space groups
  call QCSG%setQCtype(QCtype)
  allocate( QCSG%data(200,6,6)  )         !< all symmetry matrices for a given spacegroup
  allocate( QCSG%direc(100,5,5) )         !< direct space point group matrices
  allocate( QCSG%recip(100,5,5) )         !< reciprocal space point group matrices
  allocate( QCSG%c(6,6) )                  !< dummy matrix used for various computations
! get the number of space groups in the SGname array; those names are constructed on the fly...
  call QCSG%getAxialGroupNames_()
end if 

end function QCspacegroup_constructor

!--------------------------------------------------------------------------
subroutine QCspacegroup_destructor(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: QCspacegroup_destructor
!! author: MDG 
!! version: 1.0 
!! date: 01/30/22
!!
!! destructor for the QCcell_T Class
 
IMPLICIT NONE

type(QCspacegroup_T), INTENT(INOUT)  :: self 

call reportDestructor('QCspacegroup_T')

end subroutine QCspacegroup_destructor

!--------------------------------------------------------------------------
recursive subroutine setnsg_(self, nsg)
!DEC$ ATTRIBUTES DLLEXPORT :: setnsg_
  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20
  !!
  !! set nsg

IMPLICIT NONE

class(QCspacegroup_T), INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)          :: nsg

self%nsg= nsg

end subroutine setnsg_

!--------------------------------------------------------------------------
recursive function getnsg_(self) result(nsg)
!DEC$ ATTRIBUTES DLLEXPORT :: getnsg_
  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20
  !!
  !! get nsg

IMPLICIT NONE

class(QCspacegroup_T), INTENT(INOUT)  :: self
integer(kind=irg)                     :: nsg

nsg = self%nsg

end function getnsg_

!--------------------------------------------------------------------------
recursive subroutine setSGnum_(self, SGnum)
!DEC$ ATTRIBUTES DLLEXPORT :: setSGnum_
  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20
  !!
  !! set SGnum

IMPLICIT NONE

class(QCspacegroup_T), INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)          :: SGnum

self%SGnum= SGnum

end subroutine setSGnum_

!--------------------------------------------------------------------------
recursive function getSGnum_(self) result(SGnum)
!DEC$ ATTRIBUTES DLLEXPORT :: getSGnum_
  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20
  !!
  !! get SGnum

IMPLICIT NONE

class(QCspacegroup_T), INTENT(INOUT)  :: self
integer(kind=irg)                     :: SGnum

SGnum = self%SGnum

end function getSGnum_

!--------------------------------------------------------------------------
recursive subroutine setQCtype_(self, QCtype)
!DEC$ ATTRIBUTES DLLEXPORT :: setQCtype_
  !! author: MDG
  !! version: 1.0
  !! date: 01/31/22
  !!
  !! set QCtype

IMPLICIT NONE

class(QCspacegroup_T), INTENT(INOUT)    :: self
character(3),INTENT(IN)                 :: QCtype

self%QCtype = QCtype

end subroutine setQCtype_

!--------------------------------------------------------------------------
recursive function getQCtype_(self) result(QCtype)
!DEC$ ATTRIBUTES DLLEXPORT :: getQCtype_
  !! author: MDG
  !! version: 1.0
  !! date: 01/31/22
  !!
  !! get QCtype

IMPLICIT NONE

class(QCspacegroup_T), INTENT(INOUT)    :: self
character(3)                            :: QCtype

QCtype = self%QCtype

end function getQCtype_


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
recursive subroutine printSGtable_(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: printSGtable_
!! author: MDG 
!! version: 1.0 
!! date: 01/30/22
!!
!! list the space groups... 
!!
!! Daniel Rokhsar, David Wright, David Mermin, Scale equivalence of 
!! quasicrystallographic space groups, PHYS REV B, 37(14), 1987.
 
use mod_io 

IMPLICIT NONE

class(QCspacegroup_T), INTENT(INOUT)  :: self 

type(IO_T)                            :: Message 

character(25)                         :: sgname1, sgname2
character(fnlen)                      :: str, str1
integer(kind=irg)                     :: ii, nsg

nsg = self%getnsg_()

if(mod(nsg,2) .ne. 0) then
  do ii = 1, floor(float(nsg)/2.0) + 1
    if(ii .lt. floor(float(nsg)/2.0) + 1) then
      write(str,'(I3,A)') 2*ii-1,':'
      write(str1,'(I3,A)') 2*ii,':'
      sgname1 = ''
      sgname2 = ''
      sgname1 = trim(self%SGname(2*ii-1))
      sgname2 = trim(self%SGname(2*ii))
      call Message%printMessage(trim(str)//' '//sgname1//trim(str1)//' '//sgname2)
    else
      write(str,'(I3,A)')2*ii-1,':'
      sgname1 = ''
      sgname1 = trim(self%SGname(2*ii-1))
      call Message%printMessage(trim(str)//sgname1)
    end if
  end do
else
  do ii = 1, nsg/2
    write(str,'(I3,A)')2*ii-1,':'
    write(str1,'(I3,A)')2*ii,': '
    sgname1 = ''
    sgname2 = ''
    sgname1 = trim(self%SGname(2*ii-1))
    sgname2 = trim(self%SGname(2*ii))
    call Message%printMessage(trim(str)//' '//sgname1//trim(str1)//' '//sgname2)
  end do
end if

end subroutine printSGtable_

!--------------------------------------------------------------------------
recursive subroutine getAxialGroupNames_(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: getAxialGroupNames_
!! author: MDG 
!! version: 1.0 
!! date: 01/31/22
!!
!! generate the space group symbols for the axial groups of order 8, 10, or 12 

use mod_io 

class(QCspacegroup_T), INTENT(INOUT)  :: self 

type(IO_T)                            :: Message 

integer(kind=irg)                     :: n      ! axial rotational symmetry
character(1)                          :: lat    ! primitive or centered
character(fnlen)                      :: str, str1, str2    ! concatanate various symbols here
character(25)                         :: bspace, sgname1, sgname2
integer(kind=irg)                     :: nsg, styp, p, io_int(1), ii
logical                               :: pflag

! n is the rotational symmetry of the highest symmetry axis
! styp is the categorization of n based on the following conditions
! 1. styp == 1 for n power of odd prime e.g. 3-fold, 7-fold, 9-fold etc.
! 2. styp == 2 for n power of two e.g. 2-fold, 4-fold, 8-fold etc.
! 3. styp == 3 for n twice the power of odd prime e.g. 6-fold, 10-fold etc.
! 4. styp == 4 for n even but not twice a prime power e.g. 12-fold, 20-fold etc.
! 5. styp == 5 for n odd but not a prime power e.g. 15-fold, 21-fold etc.

if(trim(self%getQCtype_()) .eq. 'Oct') then
  n     = 8
  styp  = 2
else if(trim(self%getQCtype_()) .eq. 'Dec') then
  n     = 10
  styp  = 3
else if(trim(self%getQCtype_()) .eq. 'DoD') then
  n     = 12
  styp  = 4 
else
  call Message%printError('getAxialGroupNames_:','only 8, 10 and 12 fold symmetries have been implemented thus far.')
end if

select case(styp)
  case(2)
    nsg = 57 + 3*(n - 1) + 2*(n/2 - 1)
    call self%setnsg_(nsg)
    io_int(1) = nsg
    write(str,'(I2)') n
    call Message%WriteValue('Number of space groups for axial symmetry of '//trim(str)//' = ',io_int,1,'(I3,//)')
    if(.not.(allocated(self%SGname))) allocate(self%SGname(nsg))
    self%SGname = ''

  case(3)
    nsg = 17 + 2*(n - 1)
    call self%setnsg_(nsg)
    io_int(1) = nsg
    write(str,'(I2)') n
    call Message%WriteValue('Number of space groups for axial symmetry of '//trim(str)//' = ',io_int,1,'(I3,//)')
    if(.not.(allocated(self%SGname))) allocate(self%SGname(nsg))

    ! fill out the names and print on screen
    write(str1,'(I6)') n
    str1 = adjustl(str1)
    write(str,'(A,A)')'P',trim(str1)
    self%SGname(1) = trim(str)

    write(str1,'(I6)') -n
    str1 = adjustl(str1)
    write(str,'(A,A)')'P',trim(str1)
    self%SGname(2) = trim(str)

    do ii = 3, n+1
      write(str1,'(I6)') ii-2
      str1 = adjustl(str1)
      str1 = '_'//trim(str1)
      write(str2,'(I6,A)') n,trim(str1)
      str2 = adjustl(str2)
      write(str,'(A,A)')'P',trim(str2)
      self%SGname(ii) = trim(str)
    end do

    write(str1,'(I6)') n
    str1 = adjustl(str1)
    str1 = trim(str1)//' 2 2'
    write(str,'(A,A)')'P',trim(str1)
    self%SGname(n+2) = trim(str)

    do ii = n+3, 2*n+1
      write(str1,'(I6)') ii-n-2
      str1 = adjustl(str1)
      str1 = '_'//trim(str1)//' 2 2'
      write(str2,'(I6,A)') n,trim(str1)
      str2 = adjustl(str2)
      write(str,'(A,A)')'P',trim(str2)
      self%SGname(ii) = trim(str)
    end do

    write(str1,'(I6)') -n
    str1 = adjustl(str1)
    self%SGname(2*n+2) = 'P'//trim(str1)//' m 2'
    self%SGname(2*n+3) = 'P'//trim(str1)//' c 2'
    self%SGname(2*n+4) = 'P'//trim(str1)//' 2 m'
    self%SGname(2*n+5) = 'P'//trim(str1)//' 2 c'

    write(str1,'(I6)') n
    str1 = adjustl(str1)
    write(str2,'(I6)') n/2
    str2 = adjustl(str2)

    self%SGname(2*n+6) = 'P'//trim(str1)//' m m'
    self%SGname(2*n+7) = 'P'//trim(str1)//' c c'
    self%SGname(2*n+8) = 'P'//trim(str1)//'_'//trim(str2)//' m c'
    self%SGname(2*n+9) = 'P'//trim(str1)//'_'//trim(str2)//' c m'

    self%SGname(2*n+10) = 'P'//trim(str1)//' / m'
    self%SGname(2*n+11) = 'P'//trim(str1)//'_'//trim(str2)//' / m'
    self%SGname(2*n+12) = 'P'//trim(str1)//' / m 2 / m 2 / m'
    self%SGname(2*n+13) = 'P'//trim(str1)//'/ m 2 / c 2 / c'
    self%SGname(2*n+14) = 'P'//trim(str1)//'_'//trim(str2)//' / m 2 / m 2 / c'
    self%SGname(2*n+15) = 'P'//trim(str1)//'_'//trim(str2)//' / m 2 / c 2 / m'

  case(4)
    nsg = 13 + 2*(n - 1)
    call self%setnsg_(nsg)
    io_int(1) = nsg
    write(str,'(I2)') n
    call Message%WriteValue('Number of space groups for axial symmetry of '//trim(str)//' = ',io_int,1,'(I3,//)')
    if(.not.(allocated(self%SGname))) allocate(self%SGname(nsg))

    write(str1,'(I6)') n
    str1 = adjustl(str1)
    write(str,'(A,A)')'P',trim(str1)
    self%SGname(1) = trim(str)

    write(str1,'(I6)') -n
    str1 = adjustl(str1)
    write(str,'(A,A)')'P',trim(str1)
    self%SGname(2) = trim(str)

    do ii = 3, n+1
      write(str1,'(I6)') ii-2
      str1 = adjustl(str1)
      str1 = '_'//trim(str1)
      write(str2,'(I6,A)') n,trim(str1)
      str2 = adjustl(str2)
      write(str,'(A,A)')'P',trim(str2)
      self%SGname(ii) = trim(str)
    end do

    write(str1,'(I6)') n
    str1 = adjustl(str1)
    str1 = trim(str1)//' 2 2'
    write(str,'(A,A)')'P',trim(str1)
    self%SGname(n+2) = trim(str)

    do ii = n+3, 2*n+1
      write(str1,'(I6)') ii-n-2
      str1 = adjustl(str1)
      str1 = '_'//trim(str1)//' 2 2'
      write(str2,'(I6,A)') n,trim(str1)
      str2 = adjustl(str2)
      write(str,'(A,A)')'P',trim(str2)
      self%SGname(ii) = trim(str)
    end do

    write(str1,'(I6)') -n
    str1 = adjustl(str1)
    self%SGname(2*n+2) = 'P'//trim(str1)//' 2 m'
    self%SGname(2*n+3) = 'P'//trim(str1)//' 2 c'

    write(str1,'(I6)') n
    str1 = adjustl(str1)
    write(str2,'(I6)') n/2
    str2 = adjustl(str2)

    self%SGname(2*n+4) = 'P'//trim(str1)//' m m'
    self%SGname(2*n+5) = 'P'//trim(str1)//' c c'
    self%SGname(2*n+6) = 'P'//trim(str1)//'_'//trim(str2)//' c m'

    self%SGname(2*n+7) = 'P'//trim(str1)//' / m'
    self%SGname(2*n+8) = 'P'//trim(str1)//'_'//trim(str2)//' / m'
    self%SGname(2*n+9) = 'P'//trim(str1)//' / m 2 / m 2 / m'
    self%SGname(2*n+10) = 'P'//trim(str1)//'/ m 2 / c 2 / c'
    self%SGname(2*n+11) = 'P'//trim(str1)//'_'//trim(str2)//' / m 2 / c 2 / m'

  case DEFAULT
    call Message%printError('getAxialGroupNames_:','n-fold axis couldn''t be categorized into one of the 5 categories. &
      see: The space groups of axial crystals and quasicrystals, REVIEW OF MODERN PHYSICS, 63(3), 1991');
end select

end subroutine getAxialGroupNames_

!--------------------------------------------------------------------------
recursive subroutine GetQCSpaceGroup_(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: GetQCSpaceGroup_
!! author: MDG 
!! version: 1.0 
!! date: 01/30/22

use mod_io 

class(QCspacegroup_T), INTENT(INOUT)  :: self 

type(IO_T)                            :: Message 

integer(kind=irg)                     :: io_int(1) 

 ! get the space group number
 call Message%printMessage('Enter space group number', frm = "(//A)")

 call Message%ReadValue('    Space group number --> ', io_int, 1)
 self%SGnum = io_int(1)

end subroutine GetQCSpaceGroup_




end module mod_QCsymmetry