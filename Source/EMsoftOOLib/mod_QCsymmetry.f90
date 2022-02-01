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
  integer(kind=irg)           :: nsym 
  integer(kind=irg)           :: SGnum
  integer(kind=irg)           :: GENnum             !< number of generator matrices
  integer(kind=irg)           :: MATnum             !< number of non-zero symmetry matrices
  integer(kind=irg)           :: NUMpt              !< number of point group operators
  integer(kind=irg)           :: N_Axial
  logical                     :: reduce             !< switch to enable/disable reduction to fundamental cell
  real(kind=dbl),allocatable, public  :: data(:,:,:)        !< all symmetry matrices for a given spacegroup
  real(kind=dbl),allocatable, public  :: direc(:,:,:)       !< direct space point group matrices
  real(kind=dbl),allocatable, public  :: recip(:,:,:)       !< reciprocal space point group matrices
  real(kind=dbl),allocatable, public  :: icos(:,:,:)        !< dummy 6x6 matrix used for various computations
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
  procedure, pass(self) :: setNUMpt_
  procedure, pass(self) :: setMATnum_
  procedure, pass(self) :: setQCtype_
  procedure, pass(self) :: getnsg_
  procedure, pass(self) :: getSGnum_
  procedure, pass(self) :: getNUMpt_
  procedure, pass(self) :: getMATnum_
  procedure, pass(self) :: getQCtype_
  procedure, pass(self) :: getAxialGroupNames_
  procedure, pass(self) :: fillgen_QC_
  procedure, pass(self) :: MakeQCGenerators_
  procedure, pass(self) :: GenerateQCSymmetry_
  procedure, pass(self) :: matrixmult_
  procedure, pass(self) :: isnew_  
  procedure, pass(self) :: isitnew_  
  procedure, pass(self) :: GetSymmetryOperators_
  ! procedure, pass(self) :: 
  ! procedure, pass(self) :: 
  ! procedure, pass(self) :: 

  generic, public :: printSGtable => printSGtable_
  generic, public :: GetQCSpaceGroup => GetQCSpaceGroup_
  generic, public :: setnsg => setnsg_
  generic, public :: setSGnum => setSGnum_
  generic, public :: setNUMpt => setNUMpt_
  generic, public :: setMATnum => setMATnum_
  generic, public :: setQCtype => setQCtype_
  generic, public :: getnsg => getnsg_
  generic, public :: getSGnum => getSGnum_
  generic, public :: getNUMpt => getNUMpt_
  generic, public :: getMATnum => getMATnum_
  generic, public :: getQCtype => getQCtype_
  generic, public :: fillgen_QC => fillgen_QC_
  generic, public :: MakeQCGenerators => MakeQCGenerators_
  generic, public :: GenerateQCSymmetry => GenerateQCSymmetry_
  generic, public :: matrixmult => matrixmult_
  generic, public :: isnew => isnew_
  generic, public :: isitnew => isitnew_
  generic, public :: GetSymmetryOperators => GetSymmetryOperators_
  ! generic, public :: 
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
  allocate( QCSG%icos(6,6,120))
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
  allocate( QCSG%icos(5,5,40))
  allocate( QCSG%c(6,6) )                 !< dummy matrix used for various computations
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
recursive subroutine setNUMpt_(self, NUMpt)
!DEC$ ATTRIBUTES DLLEXPORT :: setNUMpt_
  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20
  !!
  !! set NUMpt

IMPLICIT NONE

class(QCspacegroup_T), INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)          :: NUMpt

self%NUMpt= NUMpt

end subroutine setNUMpt_

!--------------------------------------------------------------------------
recursive function getNUMpt_(self) result(NUMpt)
!DEC$ ATTRIBUTES DLLEXPORT :: getNUMpt_
  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20
  !!
  !! get NUMpt

IMPLICIT NONE

class(QCspacegroup_T), INTENT(INOUT)  :: self
integer(kind=irg)                     :: NUMpt

NUMpt = self%NUMpt

end function getNUMpt_

!--------------------------------------------------------------------------
recursive subroutine setMATnum_(self, MATnum)
!DEC$ ATTRIBUTES DLLEXPORT :: setMATnum_
  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20
  !!
  !! set MATnum

IMPLICIT NONE

class(QCspacegroup_T), INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)          :: MATnum

self%MATnum= MATnum

end subroutine setMATnum_

!--------------------------------------------------------------------------
recursive function getMATnum_(self) result(MATnum)
!DEC$ ATTRIBUTES DLLEXPORT :: getMATnum_
  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20
  !!
  !! get MATnum

IMPLICIT NONE

class(QCspacegroup_T), INTENT(INOUT)  :: self
integer(kind=irg)                     :: MATnum

MATnum = self%MATnum

end function getMATnum_

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
!! author: MDG/SS 
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
!! author: MDG/SS 
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
!! author: MDG/SS 
!! version: 1.0 
!! date: 01/30/22
!!
!! ask the user for space group number

use mod_io 

class(QCspacegroup_T), INTENT(INOUT)  :: self 

type(IO_T)                            :: Message 

integer(kind=irg)                     :: io_int(1) 

 ! get the space group number
 call Message%printMessage('Enter space group number', frm = "(//A)")

 call Message%ReadValue('    Space group number --> ', io_int, 1)
 self%SGnum = io_int(1)

end subroutine GetQCSpaceGroup_

!--------------------------------------------------------------------------
recursive subroutine fillgen_QC_(self, t)
!DEC$ ATTRIBUTES DLLEXPORT :: fillgen_QC_
!! author: MDG/SS
!! version: 1.0 
!! date: 02/01/22
!!
!! fills in a generator matrix based on an input string

IMPLICIT NONE

class(QCspacegroup_T),INTENT(INOUT)         :: self
character(1),INTENT(IN)                     :: t(7) !< 7-character input string
integer(kind=irg)                           :: j    !< auxiliary variable

! first fill the array with zeroes and a 1 at 4,4
 self%c(1:7,1:7) = 0.0D0
 self%c(7,7)     = 1.0D0

! then check for the particular matrix type
 select case (t(1))
 case('a')
  self%c(1,1) = 1.0D0; self%c(2,2) = 1.0D0; self%c(3,3) = 1.0D0
  self%c(4,4) = 1.0D0; self%c(5,5) = 1.0D0; self%c(6,6) = 1.0D0
 case('b')
  self%c(1,1) = 1.0D0; self%c(6,2) = 1.0D0; self%c(2,3) = 1.0D0
  self%c(3,4) = 1.0D0; self%c(4,5) = 1.0D0; self%c(5,6) = 1.0D0
 case('c')
  self%c(2,1) =  1.0D0; self%c(6,2) = 1.0D0; self%c(4,3) = -1.0D0
  self%c(5,4) = -1.0D0; self%c(3,5) = 1.0D0; self%c(1,6) =  1.0D0
 case('d')
  self%c(1,1) = -1.0D0; self%c(2,2) = -1.0D0; self%c(3,3) = -1.0D0
  self%c(4,4) = -1.0D0; self%c(5,5) = -1.0D0; self%c(6,6) = -1.0D0
 end select

! then fill in the translational component
 do j=1,6 
  select case (t(j+1))
   case('A'); self%c(j,7) = 1.0D0/5.0D0
   case('B'); self%c(j,7) = 1.0D0/2.0D0
  end select
 end do

end subroutine fillgen_QC_

!--------------------------------------------------------------------------
recursive subroutine GenerateQCSymmetry_(self, dopg)
!DEC$ ATTRIBUTES DLLEXPORT :: GenerateQCSymmetry_
!! author: MDG/SS
!! version: 1.0 
!! date: 02/01/22
!!
!! compute all relevant symmetry operators

use mod_io

IMPLICIT NONE

class(QCspacegroup_T),INTENT(INOUT)       :: self
logical,INTENT(IN)                        :: dopg

type(IO_T)                                :: Message

integer(kind=irg)                         :: i,j,k,nsym,k1,k2,l1,l2,sg(4) !< loop counters (mostly)
real(kind=dbl)                            :: q,sm, eps                    !< auxiliary variables.

eps = 0.1D0

! create the space group generator matrices
call self%MakeQCGenerators()

sg = (/5,6,10,11/)
do k = 1,4
 if(self%SGnum .eq. sg(k) ) then
  call Message%printMessage('--> Face-centered Icosahedral quasicrystal detected. Generating space group symmetry. ')
 end if
end do

nsym = self%GENnum

! generate new elements from the squares of the generators 
 do k=1,self%GENnum 
  call self%matrixmult(k,k)
  if (self%isitnew_(nsym).eqv..TRUE.) then 
   nsym=nsym+1
   self%data(nsym,:,:) = self%c(:,:)
  end if
 end do

! generate the remainder of the factorgroup
 k1=1
 do while (k1.le.nsym) 
  k2=k1+1
  do while (k2.le.nsym)
   call self%matrixmult(k2,k1)
   if (self%isitnew_(nsym).eqv..TRUE.) then 
    nsym=nsym+1
    self%data(nsym,:,:) = self%c(:,:)
   end if
   k2=k2+1
  end do
  k1=k1+1
 end do
 self%MATnum = nsym

! reduce the translation operators to the fundamental unit cell
 do i=1,self%MATnum
   self%data(i,6,7)=mod( self%data(i,6,7),1.0_dbl)
 end do 

! generate point group symmetry if flag is passed
if(dopg.eqv..TRUE.) then
  self%NUMpt = 0
  do i = 1,self%MATnum
    sm = 0.D0
    do j = 1,6
      sm = sm + self%data(i,j,7)**2
    end do
    if(sm .lt. eps) then
      self%NUMpt = self%NUMpt + 1
      self%direc(self%NUMpt,1:6,1:6) = self%data(i,1:6,1:6)
    end if
  end do
end if

end subroutine GenerateQCSymmetry_

!--------------------------------------------------------------------------
recursive subroutine MakeQCGenerators_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: MakeQCGenerators_
!! author: MDG/SS
!! version: 1.0 
!! date: 02/01/22
!!
!! get all the QC generators

IMPLICIT NONE

class(QCspacegroup_T),INTENT(INOUT)       :: self 

real(kind=dbl)                            :: g1(6,6), g2(6,6), g3(6,6)
character(100)                            :: genst                        !< full generator string
character(1)                              :: t(7), ngen
integer(kind=irg)                         :: i, k, l

genst = trim(self%GL(self%SGnum))

ngen = genst(1:1)

read(ngen,*) self%GENnum

! create the generator matrices 
do i = 1,self%GENnum
  do  k = 1,7
      l = 1 + 7*(i-1) + k
      t(k) = genst(l:l)
  end do
  call self%fillgen_QC(t)
  self%data(i,:,:) = self%c(:,:)
end do

end subroutine MakeQCGenerators_

!--------------------------------------------------------------------------
recursive subroutine matrixmult_(self, k1, k2)
!DEC$ ATTRIBUTES DLLEXPORT :: matrixmult_
!! author: MDG/SS
!! version: 1.0 
!! date: 02/01/22
!!
!! multiplies symmetry matrices together

IMPLICIT NONE

class(QCspacegroup_T),INTENT(INOUT)   :: self
integer(kind=irg),INTENT(IN)          :: k1              !< index of first 6x6 input matrix
integer(kind=irg),INTENT(IN)          :: k2              !< index of second 6x6 input matrix

integer(kind=irg)                     :: i,j,k,nn        !< loop counters
real(kind=dbl),parameter              :: eps=0.0005_dbl  !< truncation constant

if (self%getQCtype_().eq.'Ico') then 
  nn = 7 
else
  nn = 6
end if 

do i=1,nn
  do j=1,nn
    self%c(i,j) = 0.0_dbl
    do k=1,nn
      self%c(i,j)=self%c(i,j)+self%data(k1,i,k)*self%data(k2,k,j)
    end do
  end do
end do

! bring the translational part of the matrix back to
! the first unit cell and correct possible rounding errors
if (self%getQCtype_().eq.'Ico') then 
  do k = 1,6
    if (abs(self%c(k,7)).lt.eps) self%c(k,7)=0.0_dbl
    if (self%c(k,7).lt.0.0_dbl) then 
     do while (self%c(k,7).lt.0.0_dbl) 
      self%c(k,7)=self%c(k,7)+1.0_dbl
     end do
    end if
    if (self%c(k,7).gt.1.0_dbl) then 
     do while (self%c(k,7).gt.1.0_dbl) 
      self%c(k,7)=self%c(k,7)-1.0_dbl
     end do
    end if
    if (abs(self%c(k,7)-1.0_dbl).lt.eps) self%c(k,7)=0.0_dbl
  end do
else
! only a_5 is checked since the other are always 0
  if (abs(self%c(5,6)).lt.eps) self%c(5,6)=0.0_dbl
  if (self%c(5,6).lt.0.0_dbl) then 
   do while (self%c(5,6).lt.0.0_dbl) 
    self%c(5,6)=self%c(5,6)+1.0_dbl
   end do
  end if
  if (self%c(5,6).gt.1.0_dbl) then 
   do while (self%c(5,6).gt.1.0_dbl) 
    self%c(5,6)=self%c(5,6)-1.0_dbl
   end do
  end if
  if (abs(self%c(5,6)-1.0_dbl).lt.eps) self%c(5,6)=0.0_dbl
end if 

end subroutine matrixmult_

!--------------------------------------------------------------------------
recursive function MatrixPower_(s, A, n) result(B)
!DEC$ ATTRIBUTES DLLEXPORT :: MatrixPower_
!! author: MDG/SS
!! version: 1.0 
!! date: 02/01/22
!!
!! multiplies symmetry matrices together

IMPLICIT NONE

integer(kind=irg),INTENT(IN)          :: s
real(kind=dbl),INTENT(IN)             :: A(s,s)
integer(kind=irg),INTENT(IN)          :: n

real(kind=dbl)                        :: B(s,s), tmp(s,s)
integer(kind=irg)                     :: i, m

if (n.eq.0) then
  B = 0.D0
  do i = 1,m
    B(i,i) = 1.D0
  end do
else if (n.eq.1) then
  B = A
else if (n.gt.1) then
  tmp = A
  do i = 1,n-1
    tmp = matmul(A, tmp)
  end do
  B = tmp
end if

end function MatrixPower_

!--------------------------------------------------------------------------
recursive subroutine GetSymmetryOperators_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: GetSymmetryOperators_
!! author: MDG/SS
!! version: 1.0 
!! date: 02/01/22
!!
!! generate all symmetry operators using the generators in 6 index notation. 
!! The generators are borrowed from:
!! F. Gahler, CRYSTALLOGRAPHY OF DODECAGONAL QUASICRYSTALS  (Dodecahedral, Decahedral, Octahedral)
!! T. Janssen, Crystallography of Quasi-Crystals, Acta Cryst. (1986). A42, 261-271 (Icosahedral)

IMPLICIT NONE

class(QCspacegroup_T),INTENT(INOUT)   :: self

real(kind=dbl),allocatable            :: g1(:,:), g2(:,:), g3(:,:), O1(:,:,:), O2(:,:,:), O3(:,:,:), O(:,:)
integer(kind=irg)                     :: ii, jj, kk, nsym, nsym2



if (self%getQCtype_().eq.'DoD') then 
  allocate(g1(5,5), g2(5,5), g3(5,5),O1(5,5,12),O2(5,5,2), O3(5,5,2), O(5,5))

  g1(1,1:5) = (/0.D0, 1.D0, 0.D0, -1.D0, 0.D0/)
  g1(2,1:5) = (/1.D0, 0.D0, 0.D0, 0.D0, 0.D0/)
  g1(3,1:5) = (/0.D0, 1.D0, 0.D0, 0.D0, 0.D0/)
  g1(4,1:5) = (/0.D0, 0.D0, 1.D0, 0.D0, 0.D0/)
  g1(5,1:5) = (/0.D0, 0.D0, 0.D0, 0.D0, 1.D0/)

  g2(1,1:5) = (/0.D0, 0.D0, 0.D0, 1.D0, 0.D0/)
  g2(2,1:5) = (/0.D0, 0.D0, 1.D0, 0.D0, 0.D0/)
  g2(3,1:5) = (/0.D0, 1.D0, 0.D0, 0.D0, 0.D0/)
  g2(4,1:5) = (/1.D0, 0.D0, 0.D0, 0.D0, 0.D0/)
  g2(5,1:5) = (/0.D0, 0.D0, 0.D0, 0.D0, 1.D0/)

! generate all powers of the matrix
  do ii = 1,12
    O1(:,:,ii) = MatrixPower_(5,g1,ii)
    if(ii .le. 2) then
      O2(:,:,ii) = MatrixPower_(5,g2,ii)
    end if
  end do

! combine the rotations and mirrors to generate all 24 operators
! START WITH IDENTITY OPERATOR
  nsym = 1
  self%icos = 0.D0
  do ii = 1,5
    self%icos(ii,ii,1) = 1.D0
  end do

  do ii = 1,12
    O = O1(:,:,ii)
    if(self%isnew_(5, O, nsym)) then
      nsym = nsym + 1
      self%icos(:,:,nsym) = O(1:5,1:5)
    end if
  end do

  do ii = 1,2
    O = O2(:,:,ii)
    if(self%isnew_(5, O, nsym)) then
      nsym = nsym + 1
      self%icos(:,:,nsym) = O(1:5,1:5)
    end if
  end do

  ! go through all combinations two at a time (2 possible configurations)
  do ii = 1,12
    do jj = 1,2
      O = matmul(O1(:,:,ii), O2(:,:,jj))
      if(self%isnew_(5, O, nsym)) then
        nsym = nsym + 1
        self%icos(:,:,nsym) = O(1:5,1:5)
      end if

      O = matmul(O2(:,:,jj), O1(:,:,ii))
      if(self%isnew_(5, O, nsym)) then
        nsym = nsym + 1
        self%icos(:,:,nsym) = O(1:5,1:5)
      end if

    end do
  end do

  self%nsym = nsym
end if 

if (self%getQCtype_().eq.'Dec') then 
  allocate(g1(5,5), g2(5,5), g3(5,5),O1(5,5,10),O2(5,5,2), O3(5,5,2), O(5,5))

  g1(1,:) = (/0.D0,-1.D0,1.D0,0.D0,0.D0/)
  g1(2,:) = (/0.D0,-1.D0,0.D0,1.D0,0.D0/)
  g1(3,:) = (/0.D0,-1.D0,0.D0,0.D0,0.D0/)
  g1(4,:) = (/1.D0,-1.D0,0.D0,0.D0,0.D0/)
  g1(5,:) = (/0.D0,0.D0,0.D0,0.D0,1.D0/)

  ! generate all powers of the matrix
  do ii = 1,4
    O1(:,:,ii) = MatrixPower_(5,g1,ii)
  end do

! combine the rotations and mirrors to generate all 24 operators
! START WITH IDENTITY OPERATOR
  nsym = 1
  self%icos = 0.D0
  do ii = 1,5
    self%icos(ii,ii,1) = 1.D0
  end do

  do ii = 1,4
    O = O1(:,:,ii)
    if(self%isnew_(5, O, nsym)) then
      nsym = nsym + 1
      self%icos(:,:,nsym) = O(1:5,1:5)
    end if
  end do
  self%nsym = nsym
end if

if (self%getQCtype_().eq.'Oct') then 
  allocate(g1(5,5), g2(5,5), g3(5,5),O1(5,5,8),O2(5,5,2), O3(5,5,2), O(5,5))

  g1(1,1:5) = (/0.D0, 0.D0, 0.D0, -1.D0, 0.D0/)
  g1(2,1:5) = (/1.D0, 0.D0, 0.D0, 0.D0, 0.D0/)
  g1(3,1:5) = (/0.D0, 1.D0, 0.D0, 0.D0, 0.D0/)
  g1(4,1:5) = (/0.D0, 0.D0, 1.D0, 0.D0, 0.D0/)
  g1(5,1:5) = (/0.D0, 0.D0, 0.D0, 0.D0, 1.D0/)

  g2(1,1:5) = (/0.D0, 0.D0, 0.D0, 1.D0, 0.D0/)
  g2(2,1:5) = (/0.D0, 0.D0, 1.D0, 0.D0, 0.D0/)
  g2(3,1:5) = (/0.D0, 1.D0, 0.D0, 0.D0, 0.D0/)
  g2(4,1:5) = (/1.D0, 0.D0, 0.D0, 0.D0, 0.D0/)
  g2(5,1:5) = (/0.D0, 0.D0, 0.D0, 0.D0, 1.D0/)

  ! generate all powers of the matrix
  do ii = 1,8
    O1(:,:,ii) = MatrixPower_(5,g1,ii)
    if(ii .le. 2) then
      O2(:,:,ii) = MatrixPower_(5,g2,ii)
    end if
  end do

  ! combine the rotations and mirrors to generate all 24 operators
  ! START WITH IDENTITY OPERATOR
  nsym = 1
  self%icos = 0.D0
  do ii = 1,5
    self%icos(ii,ii,1) = 1.D0
  end do

  do ii = 1,8
    O = O1(:,:,ii)
    if(self%isnew_(5, O, nsym)) then
      nsym = nsym + 1
      self%icos(:,:,nsym) = O(1:5,1:5)
    end if
  end do

  do ii = 1,2
    O = O2(:,:,ii)
    if(self%isnew_(5, O, nsym)) then
      nsym = nsym + 1
      self%icos(:,:,nsym) = O(1:5,1:5)
    end if
  end do

  ! go through all combinations two at a time (2 possible configurations)
  do ii = 1,8
    do jj = 1,2
      O = matmul(O1(:,:,ii), O2(:,:,jj))
      if(self%isnew_(5, O, nsym)) then
        nsym = nsym + 1
        self%icos(:,:,nsym) = O(1:5,1:5)
      end if

      O = matmul(O2(:,:,jj), O1(:,:,ii))
      if(self%isnew_(5, O, nsym)) then
        nsym = nsym + 1
        self%icos(:,:,nsym) = O(1:5,1:5)
      end if

    end do
  end do

  self%nsym = nsym
end if

if (self%getQCtype_().eq.'Ico') then 
  allocate(g1(6,6), g2(6,6), g3(6,6),O1(6,6,5),O2(6,6,6), O3(6,6,2), O(6,6))

  g1(1,:) = (/1.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0/)
  g1(2,:) = (/0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0/)
  g1(3,:) = (/0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0/)
  g1(4,:) = (/0.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0/)
  g1(5,:) = (/0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 1.D0/)
  g1(6,:) = (/0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 0.D0/)

  g2(1,:) = (/0.D0, 0.D0, 0.D0,  0.D0, 0.D0, 1.D0/)
  g2(2,:) = (/1.D0, 0.D0, 0.D0,  0.D0, 0.D0, 0.D0/)
  g2(3,:) = (/0.D0, 0.D0, 0.D0,  0.D0, 1.D0, 0.D0/)
  g2(4,:) = (/0.D0, 0.D0, -1.D0, 0.D0, 0.D0, 0.D0/)
  g2(5,:) = (/0.D0, 0.D0, 0.D0, -1.D0, 0.D0, 0.D0/)
  g2(6,:) = (/0.D0, 1.D0, 0.D0,  0.D0, 0.D0, 0.D0/)

  ! Inversion operator
  g3 = 0.D0 
  do ii = 1,6
    g3(ii,ii) = -1.D0
  end do

  do ii = 1,5
    O1(:,:,ii) = MatrixPower_(6, g1, ii)
  end do

  do ii = 1,3
    O2(:,:,ii)  = MatrixPower_(6, g2, ii)
  end do

  do ii = 1,2
    O3(:,:,ii) = MatrixPower_(6, g3, ii)
  end do

  ! START WITH IDENTITY OPERATOR
  nsym = 1
  self%icos = 0.D0
  do ii = 1,6
    self%icos(ii,ii,1) = 1.D0
  end do

  ! go through all combinations one at a time (2 possible configurations)

  do ii = 1,5
    O = O1(:,:,ii)
    if(self%isnew_(6, O, nsym)) then
      nsym = nsym + 1
      self%icos(:,:,nsym) = O(1:6,1:6)
    end if
  end do

  do ii = 1,3
    O = O2(:,:,ii)
    if(self%isnew_(6, O, nsym)) then
      nsym = nsym + 1
      self%icos(:,:,nsym) = O(1:6,1:6)
    end if
  end do

  ! go through all combinations two at a time (2 possible configurations)
  do ii = 1,5
    do jj = 1,3
      O = matmul(O1(:,:,ii), O2(:,:,jj))
      if(self%isnew_(6, O, nsym)) then
        nsym = nsym + 1
        self%icos(:,:,nsym) = O(1:6,1:6)
      end if

      O = matmul(O2(:,:,jj), O1(:,:,ii))
      if(self%isnew_(6, O, nsym)) then
        nsym = nsym + 1
        self%icos(:,:,nsym) = O(1:6,1:6)
      end if
    end do
  end do


  ! now multiply all the genrators with the already existing symmetry operators
  ! in the SYM_icos matrix
  nsym2 = nsym
  do ii = 1,5
    do jj = 1,nsym2
      O = matmul(O1(:,:,ii),self%icos(:,:,jj))
      if(self%isnew_(6, O, nsym)) then
          nsym = nsym + 1
          self%icos(:,:,nsym) = O(1:6,1:6)
        end if

        O = matmul(self%icos(:,:,jj), O1(:,:,ii))
        if(self%isnew_(6, O, nsym)) then
          nsym = nsym + 1
          self%icos(:,:,nsym) = O(1:6,1:6)
        end if

    end do
  end do

  nsym2 = nsym
  do ii = 1,3
    do jj = 1,nsym2
      O = matmul(O2(:,:,ii),self%icos(:,:,jj))
      if(self%isnew_(6, O, nsym)) then
          nsym = nsym + 1
          self%icos(:,:,nsym) = O(1:6,1:6)
        end if

        O = matmul(self%icos(:,:,jj), O2(:,:,ii))
        if(self%isnew_(6, O, nsym)) then
          nsym = nsym + 1
          self%icos(:,:,nsym) = O(1:6,1:6)
        end if

    end do
  end do

  nsym2 = nsym
  do ii = 1,5
    do jj = 1,nsym2
      O = matmul(O1(:,:,ii),self%icos(:,:,jj))
      if(self%isnew_(6, O, nsym)) then
          nsym = nsym + 1
          self%icos(:,:,nsym) = O(1:6,1:6)
        end if

        O = matmul(self%icos(:,:,jj), O1(:,:,ii))
        if(self%isnew_(6, O, nsym)) then
          nsym = nsym + 1
          self%icos(:,:,nsym) = O(1:6,1:6)
        end if

    end do
  end do

  nsym2 = nsym
  do ii = 1,3
    do jj = 1,nsym2
      O = matmul(O2(:,:,ii),self%icos(:,:,jj))
      if(self%isnew_(6, O, nsym)) then
          nsym = nsym + 1
          self%icos(:,:,nsym) = O(1:6,1:6)
        end if

        O = matmul(self%icos(:,:,jj), O2(:,:,ii))
        if(self%isnew_(6, O, nsym)) then
          nsym = nsym + 1
          self%icos(:,:,nsym) = O(1:6,1:6)
        end if

    end do
  end do

  nsym2 = nsym
  do ii = 1,nsym2
    O = matmul(g3,self%icos(:,:,ii))
    if(self%isnew_(6, O, nsym)) then
        nsym = nsym + 1
        self%icos(:,:,nsym) = O(1:6,1:6)
    end if
  end do

  self%nsym = nsym
end if

deallocate(g1, g2, g3, O1, O2, O3, O)

end subroutine GetSymmetryOperators_

!--------------------------------------------------------------------------
recursive function isnew_(self, s, sym, nsym) result(new_sym)
!DEC$ ATTRIBUTES DLLEXPORT :: isnew_
!! author: MDG/SS
!! version: 1.0 
!! date: 02/01/22
!!
!! is this a new matrix ?

IMPLICIT NONE

class(QCspacegroup_T),INTENT(INOUT)   :: self
integer(kind=irg),INTENT(IN)          :: s
real(kind=dbl),INTENT(IN)             :: sym(s,s)
integer(kind=irg),INTENT(IN)          :: nsym

integer(kind=irg)                     :: ii, jj
real(kind=dbl)                        :: eps = 1.0D-4
logical                               :: new_sym

new_sym = .TRUE.

do ii = 1,nsym
  if(sum(abs(self%icos(:,:,ii) - sym)) .lt. eps) then
    new_sym = .FALSE.
    EXIT
  end if
end do

end function isnew_

!--------------------------------------------------------------------------
logical recursive function isitnew_(self, nsym)
!DEC$ ATTRIBUTES DLLEXPORT :: isitnew_
!! author: MDG/SS
!! version: 1.0 
!! date: 02/01/22
!!
!! is this a new matrix ?

IMPLICIT NONE

class(QCspacegroup_T),INTENT(INOUT)     :: self
integer(kind=irg),INTENT(IN)            :: nsym                 !< index of matrix to be compared

integer(kind=irg)                       :: i,j,k,n,im,jm,nm     !< loop counters
real(kind=dbl),parameter                :: eps=0.0005_dbl       !< comparison threshold

if (self%getQCtype_().eq.'Ico') then 
  im = 6 
  jm = 7 
  nm = 42
else
  im = 5 
  jm = 6 
  nm = 42
end if

k=0
n=0
do while ((k.le.nsym).and.(n.ne.nm))
  n=0
  k=k+1
  do i=1,im
    do j=1,jm
      if (abs(self%c(i,j)- self%data(k,i,j)).lt.eps) n=n+1
    end do
  end do
end do
 
 if (n.ne.nm) then 
  isitnew_ = .TRUE.
 else
  isitnew_ = .FALSE.
 end if

end function isitnew_



end module mod_QCsymmetry