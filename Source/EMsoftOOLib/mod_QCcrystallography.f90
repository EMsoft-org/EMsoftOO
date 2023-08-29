! ###################################################################
! Copyright (c) 2013-2023, Marc De Graef Research Group/Carnegie Mellon University
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
  !! class definition for the 5D and 6D quasi-crystal cell
  !!
  !! we use polymorphism to implement both the axial and icosahedral QCs, so that 
  !! we don't have to use two different routines for every single computation.

use mod_kinds
use mod_global

IMPLICIT NONE 

private :: QCextractposition 

! class definition
type, public :: QCcell_T
private 
  integer(kind=irg)                     :: atno
  integer(kind=irg)                     :: numindices
  integer(kind=irg)                     :: N_Axial 
  integer(kind=irg),allocatable         :: facts(:,:)
  integer(kind=irg),allocatable         :: Ucgindex(:)
  integer(kind=irg)                     :: imax_qc, imax_p
  integer(kind=irg)                     :: imaxz
  integer(kind=irg)                     :: imax
  logical,allocatable                   :: Ucgcalc(:)
  integer(kind=irg),allocatable         :: inverseIndex(:,:)
  real(kind=dbl)                        :: vol
  real(kind=dbl)                        :: dmin
  real(kind=dbl)                        :: dmin_qc, dmin_p  
  real(kind=dbl)                        :: gmax_orth
  real(kind=dbl)                        :: DWF
  real(kind=dbl)                        :: multiplicity
  character(1)                          :: centering   ! 'P','I','F'
  logical                               :: reduce
  character(3)                          :: QCtype
  character(fnlen)                      :: fname
  integer(kind=irg)                     :: ATOM_ntype 
  integer(kind=irg),public              :: ATOM_type(maxpasym)
  integer(kind=irg),public              :: numat(maxpasym)  
  real(kind=sgl),public                 :: ATOM_pos(maxpasym,10)
  real(kind=sgl),allocatable,public     :: apos(:,:,:)

contains
private 
  procedure, pass(self) :: getnDindex_
  procedure, pass(self) :: invertnDindex_
  procedure, pass(self) :: GetQCLatParm_
  procedure, pass(self) :: GetQCAsymPos_
  procedure, pass(self) :: SaveQCDataHDF_
  procedure, pass(self) :: ReadQCDataHDF_
  procedure, pass(self) :: DumpQXtalInfo_
  procedure, pass(self) :: setdmin_
  procedure, pass(self) :: setfname_
  procedure, pass(self) :: getfname_
  procedure, pass(self) :: setQCtype_
  procedure, pass(self) :: getQCtype_
  procedure, pass(self) :: displayPeriodicTable
  procedure, pass(self) :: CalcQCPositions_
  procedure, pass(self) :: GetQCOrbit_
  procedure, pass(self) :: isnewvector_
  procedure, pass(self) :: setMetricParametersQC_
  procedure, pass(self) :: get_atno_
  procedure, pass(self) :: get_numindices_
  procedure, pass(self) :: get_N_Axial_
  procedure, pass(self) :: get_vol_
  procedure, pass(self) :: get_gmax_orth_
  procedure, pass(self) :: get_multiplicity_
  procedure, pass(self) :: get_centering_
  procedure, pass(self) :: get_ATOM_ntype_
  procedure, pass(self) :: set_atno_
  procedure, pass(self) :: set_numindices_
  procedure, pass(self) :: set_N_Axial_
  procedure, pass(self) :: set_vol_
  procedure, pass(self) :: set_gmax_orth_
  procedure, pass(self) :: set_multiplicity_
  procedure, pass(self) :: set_centering_
  procedure, pass(self) :: set_ATOM_ntype_
  procedure, pass(self) :: setinverseIndex_
  procedure, pass(self) :: set_imaxQC5_
  procedure, pass(self) :: set_imaxQC6_

  generic, public :: getnDindex => getnDindex_
  generic, public :: invertnDindex => invertnDindex_
  generic, public :: GetQCLatParm => GetQCLatParm_
  generic, public :: GetQCAsymPos => GetQCAsymPos_
  generic, public :: SaveQCDataHDF => SaveQCDataHDF_
  generic, public :: ReadQCDataHDF => ReadQCDataHDF_
  generic, public :: DumpQXtalInfo => DumpQXtalInfo_
  generic, public :: setdmin => setdmin_
  generic, public :: setfname => setfname_
  generic, public :: getfname => getfname_
  generic, public :: setQCtype => setQCtype_
  generic, public :: getQCtype => getQCtype_
  generic, public :: CalcQCPositions => CalcQCPositions_
  generic, public :: GetQCOrbit => GetQCOrbit_
  generic, public :: isnewvector => isnewvector_
  generic, public :: setMetricParametersQC => setMetricParametersQC_
  generic, public :: get_atno => get_atno_
  generic, public :: get_numindices => get_numindices_
  generic, public :: get_N_Axial => get_N_Axial_
  generic, public :: get_vol => get_vol_
  generic, public :: get_gmax_orth => get_gmax_orth_
  generic, public :: get_multiplicity => get_multiplicity_
  generic, public :: get_centering => get_centering_
  generic, public :: get_ATOM_ntype => get_ATOM_ntype_
  generic, public :: set_atno => set_atno_
  generic, public :: set_numindices => set_numindices_
  generic, public :: set_N_Axial => set_N_Axial_
  generic, public :: set_vol => set_vol_
  generic, public :: set_gmax_orth => set_gmax_orth_
  generic, public :: set_multiplicity => set_multiplicity_
  generic, public :: set_centering => set_centering_
  generic, public :: set_ATOM_ntype => set_ATOM_ntype_
  generic, public :: setinverseIndex => setinverseIndex_
  generic, public :: set_imax => set_imaxQC5_, set_imaxQC6_
end type QCcell_T

type, public, extends(QCcell_T) :: QCcell_axial_T
private 
  real(kind=dbl)                        :: epvec(3,5), epar(5,3), scaling(5,5), scalingfact
  real(kind=dbl)                        :: dsm(5,5), rsm(5,5)
  real(kind=dbl)                        :: rmt(5,5), dmt(5,5)
  real(kind=dbl)                        :: QClatparm_a
  real(kind=dbl)                        :: QClatparm_c
  real(kind=dbl)                        :: alphaij, alphai5, alphastarij

contains
  private 
    procedure, pass(self) :: get_imaxQC5_
    procedure, pass(self) :: allocateinverseIndexQC5_
    procedure, pass(self) :: TransSpaceQC5_
    procedure, pass(self) :: getGvectorQC5_
    procedure, pass(self) :: getvectorLengthQC5_
    procedure, pass(self) :: ShapeTransformTriangle_
    procedure, pass(self) :: ShapeTransformPolygonCa_
    procedure, pass(self) :: getinverseIndexQC5_

    generic, public :: get_imax => get_imaxQC5_
    generic, public :: allocateinverseIndex => allocateinverseIndexQC5_
    generic, public :: getinverseIndex => getinverseIndexQC5_
    generic, public :: TransSpace => TransSpaceQC5_
    generic, public :: getGvector => getGvectorQC5_
    generic, public :: getvectorLength => getvectorLengthQC5_
    generic, public :: ShapeTransformTriangle => ShapeTransformTriangle_
    generic, public :: ShapeTransformPolygonCa => ShapeTransformPolygonCa_
end type QCcell_axial_T 

type, public, extends(QCcell_T) :: QCcell_icosahedral_T
private 
  real(kind=dbl)                        :: epvec(3,6), epar(6,3)
  real(kind=dbl)                        :: eovec(3,6), eperp(6,3)
  real(kind=dbl)                        :: Mp(6,6), Picos(6,6)
  real(kind=dbl)                        :: Mo(6,6), Qicos(6,6)
  real(kind=dbl)                        :: dsm(6,6), rsm(6,6)
  real(kind=dbl)                        :: dmt(6,6), rmt(6,6)
  real(kind=dbl)                        :: scaling(6,6)
  real(kind=dbl)                        :: QClatparm
  real(kind=dbl)                        :: alphaij, alphastarij

contains
  private 
    procedure, pass(self) :: get_imaxQC6_
    procedure, pass(self) :: allocateinverseIndexQC6_
    procedure, pass(self) :: TransSpaceQC6_
    procedure, pass(self) :: getGvectorQC6_
    procedure, pass(self) :: getvectorLengthQC6_
    procedure, pass(self) :: ShapeTransformPyramid_
    procedure, pass(self) :: ShapeTransformTriacontahedron_
    procedure, pass(self) :: getinverseIndexQC6_

    generic, public :: get_imax => get_imaxQC6_
    generic, public :: allocateinverseIndex => allocateinverseIndexQC6_
    generic, public :: getinverseIndex => getinverseIndexQC6_
    generic, public :: TransSpace => TransSpaceQC6_
    generic, public :: getGvector => getGvectorQC6_
    generic, public :: getvectorLength => getvectorLengthQC6_
    generic, public :: ShapeTransformPyramid => ShapeTransformPyramid_
    generic, public :: ShapeTransformTriacontahedron => ShapeTransformTriacontahedron_
end type QCcell_icosahedral_T 

! the constructor routine for this class 
interface QCcell_T
  module procedure QCcell_T_constructor
end interface QCcell_T

interface QCcell_axial_T
  module procedure QCcell_T_axial_constructor
end interface QCcell_axial_T

interface QCcell_icosahedral_T
  module procedure QCcell_T_icosahedral_constructor
end interface QCcell_icosahedral_T


contains

!--------------------------------------------------------------------------
type(QCcell_T) function QCcell_T_constructor( ) result(self)
!DEC$ ATTRIBUTES DLLEXPORT :: QCcell_T_constructor
!! author: MDG 
!! version: 1.0 
!! date: 01/30/22
!!
!! constructor for the QCcell_T Class
 
IMPLICIT NONE

end function QCcell_T_constructor

!--------------------------------------------------------------------------
type(QCcell_axial_T) function QCcell_T_axial_constructor( ) result(self)
!DEC$ ATTRIBUTES DLLEXPORT :: QCcell_T_axial_constructor
!! author: MDG 
!! version: 1.0 
!! date: 01/30/22
!!
!! constructor for the QCcell_axial_T Class
 
IMPLICIT NONE

end function QCcell_T_axial_constructor

!--------------------------------------------------------------------------
type(QCcell_icosahedral_T) function QCcell_T_icosahedral_constructor( ) result(self)
!DEC$ ATTRIBUTES DLLEXPORT :: QCcell_T_icosahedral_constructor
!! author: MDG 
!! version: 1.0 
!! date: 01/30/22
!!
!! constructor for the QCcell_icosahedral_T Class
 
IMPLICIT NONE

end function QCcell_T_icosahedral_constructor

!--------------------------------------------------------------------------
subroutine QCcell_T_destructor(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: QCcell_T_destructor
!! author: MDG 
!! version: 1.0 
!! date: 01/30/22
!!
!! destructor for the QCcell_T Class
 
IMPLICIT NONE

type(QCcell_T), INTENT(INOUT)  :: self 

call reportDestructor('QCcell_T')

end subroutine QCcell_T_destructor

!--------------------------------------------------------------------------
recursive subroutine setdmin_(self, dmin1, dmin2)
!DEC$ ATTRIBUTES DLLEXPORT :: setdmin_
  !! author: MDG
  !! version: 1.0
  !! date: 01/31/22
  !!
  !! set fname

IMPLICIT NONE

class(QCcell_T), INTENT(INOUT)       :: self
real(kind=sgl),INTENT(IN)            :: dmin1 
real(kind=sgl),INTENT(IN),OPTIONAL   :: dmin2

if (present(dmin2)) then 
  self%dmin_qc = dmin1
  self%dmin_p = dmin2 
else
  self%dmin = dmin1
end if

end subroutine setdmin_

!--------------------------------------------------------------------------
recursive subroutine setfname_(self, fname)
!DEC$ ATTRIBUTES DLLEXPORT :: setfname_
  !! author: MDG
  !! version: 1.0
  !! date: 01/31/22
  !!
  !! set fname

IMPLICIT NONE

class(QCcell_T), INTENT(INOUT)  :: self
character(fnlen),INTENT(IN)     :: fname

self%fname = trim(fname)

end subroutine setfname_

!--------------------------------------------------------------------------
recursive function getfname_(self) result(fname)
!DEC$ ATTRIBUTES DLLEXPORT :: getfname_
  !! author: MDG
  !! version: 1.0
  !! date: 01/31/22
  !!
  !! get fname

IMPLICIT NONE

class(QCcell_T), INTENT(INOUT)  :: self
character(fnlen)                :: fname

fname = trim(self%fname)

end function getfname_

!--------------------------------------------------------------------------
recursive subroutine setQCtype_(self, QCtype)
!DEC$ ATTRIBUTES DLLEXPORT :: setQCtype_
  !! author: MDG
  !! version: 1.0
  !! date: 01/31/22
  !!
  !! set QCtype

IMPLICIT NONE

class(QCcell_T), INTENT(INOUT)  :: self
character(3),INTENT(IN)         :: QCtype

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

class(QCcell_T), INTENT(INOUT)  :: self
character(3)                    :: QCtype

QCtype = self%QCtype

end function getQCtype_

!--------------------------------------------------------------------------
function get_atno_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_atno_
!! author: MDG 
!! version: 1.0 
!! date: 02/06/22
!!
!! get atno from the QCcell_T class

IMPLICIT NONE 

class(QCcell_T), INTENT(INOUT)     :: self
integer(kind=irg)                  :: out

out = self%atno

end function get_atno_

!--------------------------------------------------------------------------
subroutine set_atno_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_atno_
!! author: MDG 
!! version: 1.0 
!! date: 02/06/22
!!
!! set atno in the QCcell_T class

IMPLICIT NONE 

class(QCcell_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)      :: inp

self%atno = inp

end subroutine set_atno_

!--------------------------------------------------------------------------
function get_numindices_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_numindices_
!! author: MDG 
!! version: 1.0 
!! date: 02/06/22
!!
!! get numindices from the QCcell_T class

IMPLICIT NONE 

class(QCcell_T), INTENT(INOUT)     :: self
integer(kind=irg)                  :: out

out = self%numindices

end function get_numindices_

!--------------------------------------------------------------------------
subroutine set_numindices_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_numindices_
!! author: MDG 
!! version: 1.0 
!! date: 02/06/22
!!
!! set numindices in the QCcell_T class

IMPLICIT NONE 

class(QCcell_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)      :: inp

self%numindices = inp

end subroutine set_numindices_

!--------------------------------------------------------------------------
function get_N_Axial_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_N_Axial_
!! author: MDG 
!! version: 1.0 
!! date: 02/06/22
!!
!! get N_Axial from the QCcell_T class

IMPLICIT NONE 

class(QCcell_T), INTENT(INOUT)     :: self
integer(kind=irg)                  :: out

out = self%N_Axial

end function get_N_Axial_

!--------------------------------------------------------------------------
subroutine set_N_Axial_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_N_Axial_
!! author: MDG 
!! version: 1.0 
!! date: 02/06/22
!!
!! set N_Axial in the QCcell_T class

IMPLICIT NONE 

class(QCcell_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)      :: inp

self%N_Axial = inp

end subroutine set_N_Axial_

!--------------------------------------------------------------------------
function get_vol_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_vol_
!! author: MDG 
!! version: 1.0 
!! date: 02/06/22
!!
!! get vol from the QCcell_T class

IMPLICIT NONE 

class(QCcell_T), INTENT(INOUT)     :: self
real(kind=dbl)                     :: out

out = self%vol

end function get_vol_

!--------------------------------------------------------------------------
subroutine set_vol_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_vol_
!! author: MDG 
!! version: 1.0 
!! date: 02/06/22
!!
!! set vol in the QCcell_T class

IMPLICIT NONE 

class(QCcell_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)         :: inp

self%vol = inp

end subroutine set_vol_

!--------------------------------------------------------------------------
function get_gmax_orth_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_gmax_orth_
!! author: MDG 
!! version: 1.0 
!! date: 02/06/22
!!
!! get gmax_orth from the QCcell_T class

IMPLICIT NONE 

class(QCcell_T), INTENT(INOUT)     :: self
real(kind=dbl)                     :: out

out = self%gmax_orth

end function get_gmax_orth_

!--------------------------------------------------------------------------
subroutine set_gmax_orth_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_gmax_orth_
!! author: MDG 
!! version: 1.0 
!! date: 02/06/22
!!
!! set gmax_orth in the QCcell_T class

IMPLICIT NONE 

class(QCcell_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)         :: inp

self%gmax_orth = inp

end subroutine set_gmax_orth_

!--------------------------------------------------------------------------
function get_multiplicity_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_multiplicity_
!! author: MDG 
!! version: 1.0 
!! date: 02/06/22
!!
!! get multiplicity from the QCcell_T class

IMPLICIT NONE 

class(QCcell_T), INTENT(INOUT)     :: self
real(kind=dbl)                     :: out

out = self%multiplicity

end function get_multiplicity_

!--------------------------------------------------------------------------
subroutine set_multiplicity_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_multiplicity_
!! author: MDG 
!! version: 1.0 
!! date: 02/06/22
!!
!! set multiplicity in the QCcell_T class

IMPLICIT NONE 

class(QCcell_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)         :: inp

self%multiplicity = inp

end subroutine set_multiplicity_

!--------------------------------------------------------------------------
function get_centering_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_centering_
!! author: MDG 
!! version: 1.0 
!! date: 02/06/22
!!
!! get centering from the QCcell_T class

IMPLICIT NONE 

class(QCcell_T), INTENT(INOUT)     :: self
character(1)                       :: out

out = self%centering

end function get_centering_

!--------------------------------------------------------------------------
subroutine set_centering_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_centering_
!! author: MDG 
!! version: 1.0 
!! date: 02/06/22
!!
!! set centering in the QCcell_T class

IMPLICIT NONE 

class(QCcell_T), INTENT(INOUT)     :: self
character(1), INTENT(IN)           :: inp

self%centering = inp

end subroutine set_centering_

!--------------------------------------------------------------------------
function get_ATOM_ntype_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ATOM_ntype_
!! author: MDG 
!! version: 1.0 
!! date: 02/06/22
!!
!! get ATOM_ntype from the QCcell_T class

IMPLICIT NONE 

class(QCcell_T), INTENT(INOUT)     :: self
integer(kind=irg)                  :: out

out = self%ATOM_ntype

end function get_ATOM_ntype_

!--------------------------------------------------------------------------
subroutine set_ATOM_ntype_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ATOM_ntype_
!! author: MDG 
!! version: 1.0 
!! date: 02/06/22
!!
!! set ATOM_ntype in the QCcell_T class

IMPLICIT NONE 

class(QCcell_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)      :: inp

self%ATOM_ntype = inp

end subroutine set_ATOM_ntype_

!--------------------------------------------------------------------------
recursive subroutine set_imaxQC5_(self, imax_qc, imax_p)
!DEC$ ATTRIBUTES DLLEXPORT :: set_imaxQC5_
!! author: MDG 
!! version: 1.0 
!! date: 02/07/22
!!
!! set imax parameters in the QCcell_axial_T class

IMPLICIT NONE 

class(QCcell_T), INTENT(INOUT)        :: self
integer(kind=irg), INTENT(IN)         :: imax_qc
integer(kind=irg), INTENT(IN)         :: imax_p

self%imax_qc = imax_qc
self%imax_p  = imax_p

end subroutine set_imaxQC5_

!--------------------------------------------------------------------------
recursive subroutine set_imaxQC6_(self, imax)
!DEC$ ATTRIBUTES DLLEXPORT :: set_imaxQC6_
!! author: MDG 
!! version: 1.0 
!! date: 02/07/22
!!
!! set imax parameters in the QCcell_icosahedral_T class

IMPLICIT NONE 

class(QCcell_T), INTENT(INOUT)           :: self
integer(kind=irg), INTENT(IN)            :: imax

self%imax = imax

end subroutine set_imaxQC6_

!--------------------------------------------------------------------------
recursive subroutine get_imaxQC5_(self, imax_qc, imax_p)
!DEC$ ATTRIBUTES DLLEXPORT :: get_imaxQC5_
!! author: MDG 
!! version: 1.0 
!! date: 02/07/22
!!
!! get imax parameters in the QCcell_axial_T class

IMPLICIT NONE 

class(QCcell_axial_T), INTENT(INOUT)  :: self
integer(kind=irg), INTENT(INOUT)      :: imax_qc
integer(kind=irg), INTENT(INOUT)      :: imax_p

imax_qc = self%imax_qc
imax_p  = self%imax_p

end subroutine get_imaxQC5_

!--------------------------------------------------------------------------
recursive function get_imaxQC6_(self) result(imax)
!DEC$ ATTRIBUTES DLLEXPORT :: get_imaxQC6_
!! author: MDG 
!! version: 1.0 
!! date: 02/07/22
!!
!! get imax parameters in the QCcell_icosahedral_T class

IMPLICIT NONE 

class(QCcell_icosahedral_T), INTENT(INOUT)  :: self
integer(kind=irg)                           :: imax

imax = self%imax 

end function get_imaxQC6_

!--------------------------------------------------------------------------
recursive subroutine allocateinverseIndexQC5_(self, nLUT)
!DEC$ ATTRIBUTES DLLEXPORT :: allocateinverseIndexQC5_
  !! author: MDG
  !! version: 1.0
  !! date: 02/09/22 
  !!
  !! allocates the inverseIndex array

use mod_memory

IMPLICIT NONE

class(QCcell_axial_T),INTENT(INOUT)         :: self
integer(kind=irg),INTENT(IN)                :: nLUT

type(memory_T)                              :: mem

mem = memory_T() 

! the LUT arrays store all the Fourier coefficients etc...
call mem%alloc(self%inverseIndex, (/ nLUT, 5 /), 'self%inverseIndex', initval = 0)

end subroutine allocateinverseIndexQC5_

!--------------------------------------------------------------------------
recursive subroutine allocateinverseIndexQC6_(self, nLUT)
!DEC$ ATTRIBUTES DLLEXPORT :: allocateinverseIndexQC6_
  !! author: MDG
  !! version: 1.0
  !! date: 02/09/22 
  !!
  !! allocates the inverseIndex array

use mod_memory

IMPLICIT NONE

class(QCcell_icosahedral_T),INTENT(INOUT)   :: self
integer(kind=irg),INTENT(IN)                :: nLUT

type(memory_T)                              :: mem

mem = memory_T() 

! the LUT arrays store all the Fourier coefficients etc...
call mem%alloc(self%inverseIndex, (/ nLUT, 6 /), 'self%inverseIndex', initval = 0)

end subroutine allocateinverseIndexQC6_

!--------------------------------------------------------------------------
recursive subroutine setinverseIndex_(self, id, g)
!DEC$ ATTRIBUTES DLLEXPORT :: setinverseIndex_
  !! author: MDG
  !! version: 1.0
  !! date: 02/09/22
  !!
  !! set an entry in the inverseIndex array

IMPLICIT NONE

class(QCcell_T),INTENT(INOUT)         :: self
integer(kind=irg),INTENT(IN)          :: id
integer(kind=irg),INTENT(IN)          :: g(:)

self%inverseIndex( id, : ) = g(:)

end subroutine setinverseIndex_

!--------------------------------------------------------------------------
recursive function getinverseIndexQC5_(self, id) result(g)
!DEC$ ATTRIBUTES DLLEXPORT :: getinverseIndexQC5_
  !! author: MDG
  !! version: 1.0
  !! date: 02/09/22
  !!
  !! get an entry in the inverseIndex array

IMPLICIT NONE

class(QCcell_axial_T),INTENT(INOUT)   :: self
integer(kind=irg),INTENT(IN)          :: id
integer(kind=irg)                     :: g(5)

g(1:5) = self%inverseIndex( id, 1:5 )

end function getinverseIndexQC5_

!--------------------------------------------------------------------------
recursive function getinverseIndexQC6_(self, id) result(g)
!DEC$ ATTRIBUTES DLLEXPORT :: getinverseIndexQC6_
  !! author: MDG
  !! version: 1.0
  !! date: 02/09/22
  !!
  !! get an entry in the inverseIndex array

IMPLICIT NONE

class(QCcell_icosahedral_T),INTENT(INOUT)   :: self
integer(kind=irg),INTENT(IN)                :: id
integer(kind=irg)                           :: g(6)

g(1:6) = self%inverseIndex( id, 1:6 )

end function getinverseIndexQC6_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
recursive function getnDindex_(self, QCindex ) result(gindex)
!DEC$ ATTRIBUTES DLLEXPORT :: getnDindex_
!! author: MDG/SS 
!! version: 1.0 
!! date: 01/31/22
!!
!! Convert a 5 or 6-component Miller index to a single lookup index
 
IMPLICIT NONE

class(QCcell_T),INTENT(INOUT)     :: self
integer(kind=irg),INTENT(IN)      :: QCindex(*)
integer(kind=irg)                 :: gindex

integer(kind=irg)                 :: isize, isize_qc, isize_p 
integer(kind=irg), allocatable    :: g(:)

select type (self)
  class is (QCcell_icosahedral_T)
    allocate(g(6))
    isize  = 2 * self%imax + 1
    g(1:6) = QCindex(1:6) + self%imax 
    gindex = g(1)*isize**5 + g(2)*isize**4 + g(3)*isize**3 + g(4)*isize**2 + g(5)*isize + g(6) + 1
  class is (QCcell_axial_T)
    allocate(g(5))
    isize_qc = 2 * self%imax_qc + 1
    isize_p  = 2 * self%imax_p  + 1
    g(1:4)   = QCindex(1:4) + self%imax_qc
    g(5)    = QCindex(5)   + self%imax_p 
    gindex    = g(1) * isize_qc**3 * isize_p + g(2)*isize_qc**2 * isize_p + g(3)*isize_qc * isize_p + &
              g(4)*isize_p + g(5) + 1
  class default
end select 

end function getnDindex_

!--------------------------------------------------------------------------
recursive function invertnDindex_(self, gindex ) result(QCindex)
!DEC$ ATTRIBUTES DLLEXPORT :: invertnDindex_
!! author: MDG/SS 
!! version: 1.0 
!! date: 01/31/22
!!
!! Convert a single lookup index into the corresponding a 5 or 6-component Miller index 
 
IMPLICIT NONE

class(QCcell_T),INTENT(INOUT)     :: self
integer(kind=irg),INTENT(IN)      :: gindex
integer(kind=irg),allocatable     :: QCindex(:)

select type (self)
  class is (QCcell_icosahedral_T)
    allocate(QCindex(6))
    QCindex(1:6) = self%inverseIndex(gindex,1:6)
  class is (QCcell_axial_T)
    allocate(QCindex(5))
    QCindex(1:5) = self%inverseIndex(gindex,1:5)
  class default
end select 

end function invertnDindex_

!--------------------------------------------------------------------------
recursive subroutine GetQCLatParm_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: GetQCLatParm_
!! author: MDG/SS 
!! version: 1.0 
!! date: 01/31/22
!!
!! input of lattice parameters for 3D (icosahedral) and 2D (axial) quasi-crystal 

use mod_io

IMPLICIT NONE

class(QCcell_T),INTENT(INOUT)     :: self

type(IO_T)                        :: Message
real(kind=dbl)                    :: io_real(1)   !< double precision real input array
integer(kind=irg)                 :: std, io_int(1)

select type (self)
  class is (QCcell_icosahedral_T)
    self%QCtype         = 'Ico'
    self%alphaij        =  90.D0
    self%alphastarij    =  90.D0

    call Message%printMessage(' Using the higher dimensional cut-and-project approach', frm = "(A/)")

    ! get the lattice parameters
    call Message%printMessage('Enter lattice parameters', frm = "(//A)")

    ! in-plane lattice parameters
    call Message%ReadValue('    a_i | i = {1,2,3,4,5,6} (hyper-cube) [nm] = ', io_real, 1)
    self%QClatparm = io_real(1)
  class is (QCcell_axial_T)
    call Message%printMessage(' Select the 2D quasicrystal type (axial symmetry) : ', frm = "(//A)")
    call Message%printMessage('  1. octagonal (8-fold) quasicrystal', frm = "(A)")
    call Message%printMessage('  2. decagonal (10-fold) quasicrystal', frm = "(A)")
    call Message%printMessage('  3. dodecagonal (12-fold) quasicrystal', frm = "(A/)")
    call Message%ReadValue(' quasi-crystal type ---> ', io_int, 1)

    select case (io_int(1))
      case(1)
    ! 8-fold
        self%QCtype         = 'Oct'
        self%N_Axial        =  8
        self%alphaij        =  90.D0
        self%alphastarij    =  90.D0
        self%alphai5        =  90.D0

      case(2)
    ! 10-fold
        self%QCtype         = 'Dec'
        self%N_Axial        =  10
        self%alphaij        =  60.D0
        self%alphastarij    =  104.5D0
        self%alphai5        =  90.D0

      case(3)
    ! 12-fold
        self%QCtype         = 'DoD'
        self%N_Axial        =  12
        self%alphaij        =  120.D0
        self%alphastarij    =  60.D0
        self%alphai5        =  90.D0
    end select

    call Message%printMessage(' Using the higher dimensional cut and project approach', frm = "(A/)")


    ! get the lattice parameters
    call Message%printMessage('Enter lattice parameters', frm = "(/A)")

    ! in-plane lattice parameters
    call Message%ReadValue('    a_i | i = {1,2,3,4} (quasicrystal plane) [nm] = ', io_real, 1)
    self%QClatparm_a = io_real(1)

    ! out of plane lattice parameters
    call Message%ReadValue('    a_5 (axial direction) [nm] = ', io_real, 1)
    self%QClatparm_c = io_real(1)

    call Message%printMessage('', frm = "(/A)")

  class default
end select

end subroutine GetQCLatParm_

!--------------------------------------------------------------------------
recursive subroutine GetQCAsymPos_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: GetQCAsymPos_
!! author: MDG/SS 
!! version: 1.0 
!! date: 01/31/22
!!
!! read the atom coordinates from standard input

use mod_io

class(QCcell_T)                 :: self

type(IO_T)                      :: Message

logical                         :: more                 !< logical to determine if more atoms need to be entered
character(1)                    :: ans, list(256)       !< used for IO
real(kind=sgl)                  :: pt(10), out_real(10)   !< used to read and write asymmetric position data
integer(kind=irg)               :: j, io_int(1) , std   !< auxiliary variables

more=.TRUE.
self%ATOM_ntype = 0
call Message%printMessage(' Enter atoms in asymmetric unit ', frm = "(/A)")
call self%displayPeriodicTable()

do while (more)
  self%ATOM_ntype = self%ATOM_ntype + 1

  ! atomic number
  call Message%ReadValue(' ->  Atomic number : ', io_int, 1)
  self%ATOM_type(self%ATOM_ntype) = io_int(1)

  ! general atom coordinate
  list = (/ (' ',j=1,256) /)
  select type(self)
    class is (QCcell_icosahedral_T)
      call Message%printMessage(' ->  Fractional coordinates, site occupation, Bpar [nm^2],'&
      'Bperp [nm^2], radial atomic size (fraction of a_i) : ',&
      frm = "(A,' ')",advance="no")
    class is (QCcell_axial_T)
      call Message%printMessage(' ->  Fractional coordinates, site occupation, Bpar_11 [nm^2],'&
      'Bperp_33 [nm^2], Bperp [nm^2], radial atomic size (fraction of a_i) : ',&
      frm = "(A,' ')",advance="no")
  end select
  read (*,"(256A)") list

  ! interpret this string and extract coordinates and such ...
  call QCextractposition(list,pt,iQC=.TRUE.) 

  ! store in the appropriate component of the cell variable  
  self%ATOM_pos(self%ATOM_ntype,1:10) = pt(1:10)

  ! and write the coordinate back to the terminal  
  out_real = (/ (self%ATOM_pos(self%ATOM_ntype,j),j=1,10) /)
  call Message%WriteValue('    -> ', out_real, 10, frm = "(1x,6(F10.7,2x),3(F10.7,2x),F10.7)") 

  call Message%ReadValue(' ->  Another atom ? (y/n) ', ans, frm = "(A1)")
  if ((ans.eq.'y').or.(ans.eq.'Y')) then 
   more=.TRUE.
  else
   more=.FALSE.
  end if 
end do

end subroutine GetQCAsymPos_

!--------------------------------------------------------------------------
recursive subroutine CalcQCPositions_(self, QCSG)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcQCPositions_
!! author: MDG/SS
!! version: 1.0 
!! date: 02/01/22
!!
!! compute atom positions, always reduced to unit cell

use mod_QCsymmetry
use mod_io

IMPLICIT NONE

class(QCcell_T),INTENT(INOUT)         :: self
type(QCspacegroup_T),INTENT(INOUT)    :: QCSG

type(IO_T)                            :: Message

logical                               :: inside         !< auxiliary logical
integer(kind=irg)                     :: i,j,k,l,m,mm,icnt,ncells,n,kk,ier, io_int(3)  !< various auxiliary variables
real(kind=dbl),allocatable            :: ctmp(:,:)      !< auxiliary variables  

! make sure all coordinates are reduced to the fundamental unit cell
 self%reduce=.TRUE.

if (QCSG%getQCtype().eq.'Ico') then 
  m = 6
else
  m = 5
end if 
allocate(ctmp(QCSG%getMATnum(),m))

! main loop
! first allocate the apos variable (contains CARTESIAN coordinates
 if (allocated(self%apos)) deallocate(self%apos)
 allocate (self%apos(self%ATOM_ntype, QCSG%getMATnum(), m),stat=ier)
 if (ier.ne.0) call Message%printError('CalcQCPositions',' unable to allocate memory for array apos')
 ctmp = 0.D0
 do i=1,self%ATOM_ntype
! for each atom in the asymmetric unit
  call self%GetQCOrbit_(QCSG, ctmp, i, n)
  self%apos(i,1:n,1:m) = ctmp(1:n,1:m)
  self%numat(i)        = n
 end do  

end subroutine CalcQCPositions_

!--------------------------------------------------------------------------
recursive subroutine GetQCOrbit_(self, QCSG, orbit, mm, nn)
!DEC$ ATTRIBUTES DLLEXPORT ::GetQCOrbit_
!! author: MDG/SS
!! version: 1.0 
!! date: 02/01/22
!!
!! generate all symmetrically equivalent direct space vectors

use mod_QCsymmetry

IMPLICIT NONE

class(QCcell_T),INTENT(INOUT)         :: self
type(QCspacegroup_T),INTENT(INOUT)    :: QCSG
integer(kind=irg),INTENT(IN)          :: mm
real(kind=dbl),allocatable,INTENT(OUT):: orbit(:,:)
integer(kind=irg),INTENT(OUT)         :: nn


integer(kind=irg)                     :: ii, jj, kk, Pmdims, m
real(kind=dbl),allocatable            :: r(:), s(:)

Pmdims = QCSG%getMATnum()
if (QCSG%getQCtype().eq.'Ico') then 
  m = 6
else
  m = 5
end if 
allocate(orbit(Pmdims,m), r(m), s(m))

nn = 1
orbit(1,1:m) = self%ATOM_pos(mm,1:m)
r(1:m)       = orbit(1,1:m)

! get all the equivalent atom positions
do ii = 1,Pmdims
  do jj = 1,m
    s(jj) = QCSG%data(ii,jj,m+1)
    do kk = 1,m
      s(jj) = s(jj) + QCSG%data(ii,jj,kk) * r(kk)
    end do
  end do

  do jj = 1,m
    if(s(jj) .lt. 0.D0) then
      do while(s(jj) .lt. 0.D0) 
        s(jj) = s(jj) + 1.D0
      end do
    else if(s(jj) .gt. 0.D0) then
      s(jj) = mod(s(jj),1.D0)
    end if
  end do

  if(self%isnewvector(Pmdims, m, s, orbit, nn)) then
    nn = nn + 1
    orbit(nn,1:m) = s
  end if
end do

end subroutine GetQCOrbit_

!--------------------------------------------------------------------------
recursive function isnewvector_(self, Pmdims, m, hkl, orbit, nn) result(isnew)
!DEC$ ATTRIBUTES DLLEXPORT :: isnewvector_

use mod_QCsymmetry

IMPLICIT NONE

class(QCcell_T),INTENT(INOUT)         :: self
integer(kind=irg),INTENT(IN)          :: Pmdims
integer(kind=irg),INTENT(IN)          :: m
real(kind=dbl),INTENT(IN)             :: hkl(m)
real(kind=dbl),INTENT(IN)             :: orbit(Pmdims,m)
integer(kind=irg),INTENT(IN)          :: nn
logical                               :: isnew

integer(kind=irg)                     :: ii
real(kind=dbl),parameter              :: eps = 1.0D-6

isnew = .TRUE.

do ii = 1,nn
  if(sum(abs(hkl - orbit(ii,1:m))) .lt. eps) then
    isnew = .FALSE.
    EXIT
  end if
end do

end function isnewvector_

!--------------------------------------------------------------------------
recursive subroutine SaveQCDataHDF_(self, QCSG, EMsoft)
!DEC$ ATTRIBUTES DLLEXPORT :: SaveQCDataHDF_
!! author: MDG/SS 
!! version: 1.0 
!! date: 01/31/22
!!
!! save 2D or 3D QCcrystal structure data to an HDF file

use mod_EMsoft
use mod_QCsymmetry
use mod_io
use HDF5
use mod_HDFsupport
use mod_timing
use stringconstants
use ISO_C_BINDING

IMPLICIT NONE

class(QCcell_T),INTENT(INOUT)           :: self
type(QCspacegroup_T),INTENT(INOUT)      :: QCSG
type(EMsoft_T),INTENT(INOUT)            :: EMsoft

type(HDF_T)                             :: HDF
type(IO_T)                              :: Message
type(Timing_T)                          :: timer

character(11)                           :: dstr
character(15)                           :: tstr
character(fnlen)                        :: progname = 'EMQCmkxtal.f90', groupname, dataset, fname
integer(kind=irg)                       :: hdferr
real(kind=dbl)                          :: cellparams(3), cellparams5(5)
integer(kind=irg),allocatable           :: atomtypes(:)
real(kind=sgl),allocatable              :: atompos(:,:)
integer(kind=irg)                       :: naxial
character(fnlen,kind=c_char)            :: line2(1)


timer = Timing_T()
tstr = timer%getTimeString()
dstr = timer%getDateString()

! Initialize FORTRAN interface if needed.
call openFortranHDFInterface()

! generate the filepath and open the file for writing
fname = EMsoft%generateFilePath('EMXtalFolderpathname',self%fname)
hdferr = HDF%createFile(fname)

groupname = SC_CrystalData
hdferr = HDF%createGroup(groupname)
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%createGroup:'//trim(groupname), hdferr)

dataset = SC_ProgramName
line2(1) = trim(progname)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetStringArray:'//trim(dataset), hdferr)

dataset = SC_CreationDate
line2(1) = dstr
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetStringArray:'//trim(dataset), hdferr)

dataset = SC_CreationTime
line2(1) = tstr
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetStringArray:'//trim(dataset), hdferr)

dataset = SC_Creator
line2(1) = trim(EMsoft%getConfigParameter('UserName'))
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetStringArray:'//trim(dataset), hdferr)

select type (self)
  class is (QCcell_icosahedral_T)
    dataset = SC_LatticeParameters
    cellparams = (/ self%QClatparm, self%alphaij, self%alphastarij /)
    hdferr = HDF%writeDatasetDoubleArray(dataset, cellparams, 3)
    if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetDoubleArray:'//trim(dataset), hdferr)

    dataset = SC_AxialSymmetry 
    naxial = 0
    hdferr = HDF%writeDatasetInteger(dataset, naxial)

  class is (QCcell_axial_T)
    dataset = SC_LatticeParameters
    cellparams5 = (/ self%QClatparm_a, self%QClatparm_c, self%alphaij, self%alphai5, self%alphastarij /)
    hdferr = HDF%writeDatasetDoubleArray(dataset, cellparams5, 5)
    if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetDoubleArray:'//trim(dataset), hdferr)

    dataset = SC_AxialSymmetry !'Axial Symmetry'
    hdferr = HDF%writeDatasetInteger(dataset, self%N_Axial)

  class default
end select

dataset = SC_SpaceGroupNumber
hdferr = HDF%writeDatasetInteger(dataset, QCSG%getSGnum())
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetInteger:'//trim(dataset), hdferr)

dataset = SC_Natomtypes
hdferr = HDF%writeDatasetInteger(dataset, self%ATOM_ntype)
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetInteger:'//trim(dataset), hdferr)

allocate(atomtypes(self%ATOM_ntype))
atomtypes(1:self%ATOM_ntype) = self%ATOM_type(1:self%ATOM_ntype)
dataset = SC_Atomtypes
hdferr = HDF%writeDatasetIntegerArray(dataset, atomtypes, self%ATOM_ntype)
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetIntegerArray:'//trim(dataset), hdferr)
deallocate(atomtypes)

allocate(atompos(self%ATOM_ntype,10))
atompos(1:self%ATOM_ntype,1:10) = self%ATOM_pos(1:self%ATOM_ntype,1:10)
dataset = SC_AtomData
hdferr = HDF%writeDatasetFloatArray(dataset, atompos, self%ATOM_ntype, 10)
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetFloatArray:'//trim(dataset), hdferr)
deallocate(atompos)

call HDF%popall()
call closeFortranHDFInterface()

end subroutine SaveQCDataHDF_


!--------------------------------------------------------------------------
recursive subroutine ReadQCDataHDF_(self, QCSG, EMsoft)
!DEC$ ATTRIBUTES DLLEXPORT :: ReadQCDataHDF_
!! author: MDG/SS 
!! version: 1.0 
!! date: 01/31/22
!!
!! read 2D or 3D QCcrystal structure data from an HDF file

use mod_EMsoft
use mod_io
use HDF5
use mod_HDFsupport
use mod_QCsymmetry
use stringconstants
 
IMPLICIT NONE

class(QCcell_T),INTENT(INOUT)           :: self
type(QCspacegroup_T),INTENT(INOUT)      :: QCSG
type(EMsoft_T),INTENT(INOUT)            :: EMsoft

type(HDF_T)                             :: HDF
type(IO_T)                              :: Message

character(fnlen)                        :: dataset, groupname, fname
integer(HSIZE_T)                        :: dims(1), dims2(2)
integer(kind=irg)                       :: hdferr, SGnum
real(kind=dbl),allocatable              :: cellparams(:)
integer(kind=irg),allocatable           :: atomtypes(:)
real(kind=sgl),allocatable              :: atompos(:,:)
character(fnlen)                        :: pp

! Initialize FORTRAN interface if needed.
call openFortranHDFInterface()

HDF = HDF_T()

fname = EMsoft%generateFilePath('EMXtalFolderpathname',self%fname)
hdferr = HDF%openFile(fname, readonly = .TRUE.)
if (hdferr.ne.0) call HDF%error_check('ReadQCDataHDF:HDF%openFile:'//trim(fname), hdferr)

groupname = SC_CrystalData
hdferr = HDF%openGroup(groupname)
if (hdferr.ne.0) call HDF%error_check('ReadQCDataHDF:HDF%openGroup:'//trim(groupname), hdferr)

dataset = SC_LatticeParameters
call HDF%readDatasetDoubleArray(dataset, dims, hdferr, cellparams)
if (hdferr.ne.0) call HDF%error_check('ReadQCDataHDF:HDF%readDatasetDoubleArray1D:'//trim(dataset), hdferr)

select type (self)
  class is (QCcell_axial_T)
    self%QClatparm_a = cellparams(1)
    self%QClatparm_c = cellparams(2)
    self%alphaij     = cellparams(3)
    self%alphai5     = cellparams(4)
    self%alphastarij = cellparams(5)

    dataset = SC_AxialSymmetry
    call HDF%readDatasetInteger(dataset, hdferr, self%N_Axial) 
    if (hdferr.ne.0) call HDF%error_check('ReadQCDataHDF:HDF%readDatasetInteger:'//trim(dataset), hdferr)

    if(self%N_Axial .eq. 8) then
      call self%setQCType('Oct')
      QCSG = QCspacegroup_T( nD = 2, QCtype = 'Oct' )
    else if(self%N_Axial .eq. 10) then
      call self%setQCType('Dec')
      QCSG = QCspacegroup_T( nD = 2, QCtype = 'Dec' )
    else if(self%N_Axial .eq. 12) then
      call self%setQCType('DoD')
      QCSG = QCspacegroup_T( nD = 2, QCtype = 'DoD' )
    else
      call Message%printError('ReadQCDataHDF',&
      'The axial symmetry is not one of the implemented ones (only 8, 10 and 12 fold implemented.)')
    end if


  class is (QCcell_icosahedral_T)
    self%QClatparm   = cellparams(1)
    self%alphaij     = cellparams(2)
    self%alphastarij = cellparams(3)
    call self%setQCType('Ico')
    QCSG = QCspacegroup_T( nD = 3 )
end select

dataset = SC_SpaceGroupNumber
call HDF%readDatasetInteger(dataset, hdferr, SGnum) 
call QCSG%setSGnum(SGnum)
if (hdferr.ne.0) call HDF%error_check('ReadQCDataHDF:HDF%readDatasetInteger:'//trim(dataset), hdferr)

dataset = SC_Natomtypes
call HDF%readDatasetInteger(dataset, hdferr, self%ATOM_ntype)
if (hdferr.ne.0) call HDF%error_check('ReadDataHDF:HDF%readDatasetInteger:'//trim(dataset), hdferr)

dataset = SC_Atomtypes
call HDF%readDatasetIntegerArray(dataset, dims, hdferr, atomtypes)
if (hdferr.ne.0) call HDF%error_check('ReadDataHDF:HDF%readDatasetIntegerArray1D:'//trim(dataset), hdferr)

self%ATOM_type(1:self%ATOM_ntype) = atomtypes(1:self%ATOM_ntype) 
deallocate(atomtypes)

dataset = SC_AtomData
call HDF%readDatasetFloatArray(dataset, dims2, hdferr, atompos)
if (hdferr.ne.0) call HDF%error_check('ReadDataHDF:HDF%readDatasetFloatArray2D:'//trim(dataset), hdferr)

self%ATOM_pos(1:self%ATOM_ntype,1:10) = atompos(1:self%ATOM_ntype,1:10) 
deallocate(atompos)

call HDF%popall()

end subroutine ReadQCDataHDF_

!--------------------------------------------------------------------------
recursive subroutine DumpQXtalInfo_(self, QCSG)    
!DEC$ ATTRIBUTES DLLEXPORT :: DumpQXtalInfo_

use mod_io
use mod_QCsymmetry

IMPLICIT NONE

class(QCcell_T),INTENT(INOUT)           :: self
type(QCspacegroup_T),INTENT(INOUT)      :: QCSG

type(IO_T)                              :: Message

integer(kind=irg)                       :: i, j, oi_int(3)
real(kind=dbl)                          :: oi_real(10)

select type (self)
  class is (QCcell_axial_T)
     call Message%printMessage('', frm = "(A/)")
     call Message%printMessage('Quasicrystal Structure Information', frm = "('-->',A,'<--')")
     oi_real(1) = self%QClatparm_a
     call Message%WriteValue('  a_i | i = {1,2,3,4} [nm]             : ', oi_real, 1, "(F9.5)")
     oi_real(1) = self%QClatparm_c
     call Message%WriteValue('  a_5 [nm]                             : ', oi_real, 1, "(F9.5)")
     oi_real(1) = self%alphaij
     call Message%WriteValue('  alpha_ij | (i,j) = {1,2,3,4} [deg]   : ', oi_real, 1, "(F10.5)")
      oi_real(1) = self%alphai5
     call Message%WriteValue('  alpha_i5 | (i)   = {1,2,3,4} [deg]   : ', oi_real, 1, "(F10.5)")
     oi_int(1)  = self%N_Axial
     call Message%WriteValue('  Highest axial rotational symmetry    : ', oi_int, 1, "(I4)")
     oi_real(1) = self%vol
     call Message%WriteValue('  Volume [nm^5]                        : ', oi_real, 1, "(F12.8)")
     oi_int(1) = QCSG%getSGnum()
     call Message%WriteValue('  Space group #                        : ', oi_int, 1, "(1x,I3)")
     call Message%WriteValue('  Space group symbol                   : ', '  '//trim(QCSG%SGname(oi_int(1))) )
     oi_int(1) = QCSG%getnsym()
     call Message%WriteValue('  Number of symmetry operators         : ', oi_int, 1, "(1x,I4)")
     
    ! generate atom positions and dump output  
     call Message%printMessage('', frm = "(A/)")
     call self%CalcQCPositions(QCSG)
     oi_int(1) = self%ATOM_ntype
     call Message%WriteValue('  Number of asymmetric atom positions ', oi_int, 1)
     do i=1,self%ATOM_ntype
      oi_int(1:3) = (/i, self%ATOM_type(i), self%numat(i)/)
      call Message%WriteValue('  General position / atomic number / multiplicity :', oi_int, 3,"(1x,I3,'/',I2,'/',I3)",advance="no")
      call Message%printMessage(' ('//ATOM_sym(self%ATOM_type(i))//')', frm = "(A)")
      call Message%printMessage('   Equivalent positions  (a_1 a_2 a_3 a_4 a_5  occ  Bpar_11 Bpar_33 Bperp lambda_k) ', frm = "(A)")
      do j=1,self%numat(i)
        oi_real(1:10) = (/dble(self%apos(i, j,1:5)),dble(self%ATOM_pos(i,6:10))/)
        call Message%WriteValue('         > ', oi_real, 10,"(2x,6(F9.5,','),4F9.5)")
      end do
    end do
    call Message%printMessage('', frm = "(A/)")
  class is (QCcell_icosahedral_T)
     call Message%printMessage('', frm = "(A/)")
     call Message%printMessage('Quasicrystal Structure Information', frm = "('-->',A,'<--')")
     oi_real(1) = self%QClatparm
     call Message%WriteValue('  a_i | i = {1,2,3,4,5} [nm]           : ', oi_real, 1, "(F9.5)")
     oi_real(1) = self%alphaij
     call Message%WriteValue('  alpha_ij  [deg]                      : ', oi_real, 1, "(F10.5)")
     oi_real(1) = self%vol
     call Message%WriteValue('  Volume [nm^5]                        : ', oi_real, 1, "(F12.8)")
     oi_int(1) = QCSG%getSGnum()
     call Message%WriteValue('  Space group #                        : ', oi_int, 1, "(1x,I3)")
     call Message%WriteValue('  Space group symbol                   : ', '  '//trim(QCSG%SGname(oi_int(1))) )
     call Message%WriteValue('  Space group generator string         : ', '  '//trim(QCSG%GL(oi_int(1))) )
     
    ! generate atom positions and dump output  
     call Message%printMessage('', frm = "(A/)")
     call self%CalcQCPositions(QCSG)
     oi_int(1) = self%ATOM_ntype
     call Message%WriteValue('  Number of asymmetric atom positions ', oi_int, 1)
     do i=1,self%ATOM_ntype
      oi_int(1:3) = (/i, self%ATOM_type(i), self%numat(i)/)
      call Message%WriteValue('  General position / atomic number / multiplicity :', oi_int, 3,"(1x,I3,'/',I2,'/',I3)",advance="no")
      call Message%printMessage(' ('//ATOM_sym(self%ATOM_type(i))//')', frm = "(A)")
      call Message%printMessage('   Equivalent positions  (a_1 a_2 a_3 a_4 a_5 a_6  occ  Bpar Bperp lambda_k) ', frm = "(A)")
      do j=1,self%numat(i)
        oi_real(1:10) = (/dble(self%apos(i, j,1:6)),dble(self%ATOM_pos(i,7:10))/)
        call Message%WriteValue('         > ', oi_real, 10,"(2x,6(F9.5,','),4F9.5)")
      end do
    end do
    call Message%printMessage('', frm = "(A/)")
end select

end subroutine DumpQXtalInfo_


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
recursive subroutine setMetricParametersQC_(self)
!DEC$ ATTRIBUTES DLLEXPORT ::setMetricParametersQC_ 

use mod_math

IMPLICIT NONE

class(QCcell_T),INTENT(INOUT)       :: self

real(kind=dbl)                    :: QClatparm_a, QClatparm_c, QClatparm
real(kind=dbl),allocatable        :: s(:), c(:), s2(:), c2(:)
real(kind=dbl)                    :: tmp(5,3), c3(4), s3(4), c4(4), s4(4), A, B, CC
integer(kind=irg)                 :: i, j
real(kind=dbl),parameter          :: ct = 1.D0/dsqrt(5.D0),&
                                     st = dsqrt(1.D0 - ct * ct)

select type (self) 
  class is (QCcell_icosahedral_T)
    allocate(s(5), c(5), c2(5), s2(5))
    QClatparm = self%QClatparm
  class is (QCcell_axial_T)
    allocate(s(4), c(4), c2(4), s2(4))
    QClatparm_a = self%QClatparm_a
    QClatparm_c = self%QClatparm_c
end select

if(self%getQCtype() .eq. 'DoD') then
  do i = 1,4
    c3(i) = dcos(2.D0*cPi*dble(i-1)/12.D0)
    s3(i) = dsin(2.D0*cPi*dble(i-1)/12.D0)
    c4(i) = -dcos(10.D0*cPi*dble(i-1)/12.D0)
    s4(i) = -dsin(10.D0*cPi*dble(i-1)/12.D0)
  end do

  do i = 0,1
    c(i+1)  = dcos(2.D0*cPi*dble(i-1)/12.D0)
    s(i+1)  = dsin(2.D0*cPi*dble(i-1)/12.D0)
    c2(i+1) = dcos(10.D0*cPi*dble(i-1)/12.D0)
    s2(i+1) = dsin(10.D0*cPi*dble(i-1)/12.D0)    
  end do

  do i = 2,3
    c(i+1)  = dcos(2.D0*cPi*dble(i+1)/12.D0)
    s(i+1)  = dsin(2.D0*cPi*dble(i+1)/12.D0)
    c2(i+1) = dcos(10.D0*cPi*dble(i+1)/12.D0)
    s2(i+1) = dsin(10.D0*cPi*dble(i+1)/12.D0)
  end do

  select type (self) 
    class is (QCcell_axial_T)
      self%scaling(1,1:5) = (/1.D0, 1.D0, 0.D0, -1.D0, 0.D0/)
      self%scaling(2,1:5) = (/1.D0, 1.D0, 1.D0, 0.D0, 0.D0/)
      self%scaling(3,1:5) = (/0.D0, 1.D0, 1.D0, 1.D0, 0.D0/)
      self%scaling(4,1:5) = (/-1.D0, 0.D0, 1.D0, 1.D0, 0.D0/)
      self%scaling(5,1:5) = (/0.D0, 0.D0, 0.D0, 0.D0, 1.D0/)

      self%dsm(1,1:5) = (/c(1), s(1), c2(1), s2(1), 0.D0/)
      self%dsm(2,1:5) = (/c(2), s(2), c2(2), s2(2), 0.D0/)
      self%dsm(3,1:5) = (/c(3), s(3), c2(3), s2(3), 0.D0/)
      self%dsm(4,1:5) = (/c(4), s(4), c2(4), s2(4), 0.D0/)
      self%dsm(5,1:5) = (/0.D0, 0.D0, 0.D0, 0.D0, sqrt(3.D0)*QClatparm_c/(sqrt(2.D0)*QClatparm_a)/) ! 

      self%dsm(1:5,1:5) = (sqrt(2.D0)*self%QClatparm_a/sqrt(3.D0)) * self%dsm(1:5,1:5)

      self%rsm(1,1:5) = (/c3(1), s3(1), c4(1), s4(1), 0.D0/)
      self%rsm(2,1:5) = (/c3(2), s3(2), c4(2), s4(2), 0.D0/)
      self%rsm(3,1:5) = (/c3(3), s3(3), c4(3), s4(3), 0.D0/)
      self%rsm(4,1:5) = (/c3(4), s3(4), c4(4), s4(4), 0.D0/)
      self%rsm(5,1:5) = (/0.D0, 0.D0, 0.D0, 0.D0, sqrt(2.D0)*QClatparm_a/QClatparm_c/) 

      self%rsm(1:5,1:5) = (1.D0/QClatparm_a/dsqrt(2.D0)) * self%rsm(1:5,1:5) 

      self%dsm = transpose(self%dsm)
      self%rsm = transpose(self%rsm)

      A   = (sqrt(2.D0) * QClatparm_a / dsqrt(3.D0)) ** 2
      B   = QClatparm_c ** 2
      CC  = A * dcos(self%alphaij * cPi / 180.D0)

      self%dmt(1,1:5) = (/A, 0.D0, CC, 0.D0, 0.D0/)
      self%dmt(2,1:5) = (/0.D0, A, 0.D0, CC, 0.D0/)
      self%dmt(3,1:5) = (/CC, 0.D0, A, 0.D0, 0.D0/)
      self%dmt(4,1:5) = (/0.D0, CC, 0.D0, A, 0.D0/)
      self%dmt(5,1:5) = (/0.D0, 0.D0, 0.D0, 0.D0, B/)

      self%vol  = sqrt(CalcDeterminant(self%dmt, 5, 5))

      A   = 2.D0 * (1.D0 / QClatparm_a) ** 2
      B   = (1.D0 / QClatparm_c) ** 2
      CC  = A * dcos(self%alphastarij * cPi / 180.D0)

      self%rmt(1,1:5) = (/A, 0.D0, CC, 0.D0, 0.D0/)
      self%rmt(2,1:5) = (/0.D0, A, 0.D0, CC, 0.D0/)
      self%rmt(3,1:5) = (/CC, 0.D0, A, 0.D0, 0.D0/)
      self%rmt(4,1:5) = (/0.D0, CC, 0.D0, A, 0.D0/)
      self%rmt(5,1:5) = (/0.D0, 0.D0, 0.D0, 0.D0, B/)
    end select
end if

if(self%getQCtype() .eq. 'Dec') then
  do i = 1,4
    c(i)  = dcos(2.D0*cPi*dble(i)/5.D0) 
    s(i)  = dsin(2.D0*cPi*dble(i)/5.D0) 
    c2(i) = dcos(6.D0*cPi*dble(i)/5.D0)
    s2(i) = dsin(6.D0*cPi*dble(i)/5.D0)
  end do

  select type (self)
    class is (QCcell_axial_T)
      self%scaling(1,1:5)    = (/0.D0, 1.D0, 0.D0, -1.D0, 0.D0/)
      self%scaling(2,1:5)    = (/0.D0, 1.D0, 1.D0, -1.D0, 0.D0/)
      self%scaling(3,1:5)    = (/-1.D0, 1.D0, 1.D0, 0.D0, 0.D0/)
      self%scaling(4,1:5)    = (/-1.D0, 0.D0, 1.D0, -1.D0, 0.D0/)
      self%scaling(5,1:5)    = (/0.D0, 0.D0, 0.D0, 0.D0, 1.D0/)

      self%dsm(1,1:5) = (/c(1)-1.D0, s(1), c2(1)-1.D0, s2(1), 0.D0/)
      self%dsm(2,1:5) = (/c(2)-1.D0, s(2), c2(2)-1.D0, s2(2), 0.D0/)
      self%dsm(3,1:5) = (/c(3)-1.D0, s(3), c2(3)-1.D0, s2(3), 0.D0/)
      self%dsm(4,1:5) = (/c(4)-1.D0, s(4), c2(4)-1.D0, s2(4), 0.D0/)
      self%dsm(5,1:5) = (/0.D0, 0.D0, 0.D0, 0.D0, 5.D0 * QClatparm_c/(2.D0*QClatparm_a)/) ! 

      self%dsm(1:5,1:5) = (2.D0*QClatparm_a/5.D0) * self%dsm(1:5,1:5) 

      self%rsm(1,1:5) = (/c(1), s(1), c2(1), s2(1), 0.D0/)
      self%rsm(2,1:5) = (/c(2), s(2), c2(2), s2(2), 0.D0/)
      self%rsm(3,1:5) = (/c(3), s(3), c2(3), s2(3), 0.D0/)
      self%rsm(4,1:5) = (/c(4), s(4), c2(4), s2(4), 0.D0/)
      self%rsm(5,1:5) = (/0.D0, 0.D0, 0.D0, 0.D0, QClatparm_a/QClatparm_c/) 

      self%rsm(1:5,1:5) = (1.D0/QClatparm_a) * self%rsm(1:5,1:5) 

      self%dsm = transpose(self%dsm)
      self%rsm = transpose(self%rsm)

      A   = (2.D0 * QClatparm_a / dsqrt(5.D0)) ** 2
      B   = QClatparm_c ** 2
      CC  = A * dcos(self%alphaij * cPi / 180.D0)

      self%dmt(1,1:5) = (/A, CC, CC, CC, 0.D0/)
      self%dmt(2,1:5) = (/CC, A, CC, CC, 0.D0/)
      self%dmt(3,1:5) = (/CC, CC, A, CC, 0.D0/)
      self%dmt(4,1:5) = (/CC, CC, CC, A, 0.D0/)
      self%dmt(5,1:5) = (/0.D0, 0.D0, 0.D0, 0.D0, B/)

      self%vol  = sqrt(CalcDeterminant(self%dmt, 5, 5))

      A   = 2.D0 * (1.D0 / QClatparm_a) ** 2
      B   = (1.D0 / QClatparm_c) ** 2
      CC  = A * dcos(self%alphastarij * cPi / 180.D0)

      self%rmt(1,1:5) = (/A, CC, CC, CC, 0.D0/)
      self%rmt(2,1:5) = (/CC, A, CC, CC, 0.D0/)
      self%rmt(3,1:5) = (/CC, CC, A, CC, 0.D0/)
      self%rmt(4,1:5) = (/CC, CC, CC, A, 0.D0/)
      self%rmt(5,1:5) = (/0.D0, 0.D0, 0.D0, 0.D0, B/)
    end select
end if 

if(self%getQCtype() .eq. 'Oct') then
  do i = 1,4
    c(i)  = dcos(2.D0*cPi*dble(i)/8.D0)
    s(i)  = dsin(2.D0*cPi*dble(i)/8.D0)
    c2(i) = dcos(6.D0*cPi*dble(i)/8.D0)
    s2(i) = dsin(6.D0*cPi*dble(i)/8.D0)
    c3(i) = dcos(2.D0*cPi*dble(i-1)/8.D0)
    s3(i) = dsin(2.D0*cPi*dble(i-1)/8.D0)
    c4(i) = dcos(6.D0*cPi*dble(i-1)/8.D0)
    s4(i) = dsin(6.D0*cPi*dble(i-1)/8.D0)
  end do

  select type (self)
    class is (QCcell_axial_T)
      self%scaling(1,1:5)    = (/1.D0, 1.D0, 0.D0, -1.D0, 0.D0/)
      self%scaling(2,1:5)    = (/1.D0, 1.D0, 1.D0, 0.D0, 0.D0/)
      self%scaling(3,1:5)    = (/0.D0, 1.D0, 1.D0, 1.D0, 0.D0/)
      self%scaling(4,1:5)    = (/-1.D0, 0.D0, 1.D0, 1.D0, 0.D0/)
      self%scaling(5,1:5)    = (/0.D0, 0.D0, 0.D0, 0.D0, 1.D0/)

      self%dsm(1,1:5) = (/c(1), s(1), c2(1), s2(1), 0.D0/)
      self%dsm(2,1:5) = (/c(2), s(2), c2(2), s2(2), 0.D0/)
      self%dsm(3,1:5) = (/c(3), s(3), c2(3), s2(3), 0.D0/)
      self%dsm(4,1:5) = (/c(4), s(4), c2(4), s2(4), 0.D0/)
      self%dsm(5,1:5) = (/0.D0, 0.D0, 0.D0, 0.D0, dsqrt(2.D0)*QClatparm_c/QClatparm_a/) ! 

      self%dsm(1:5,1:5) = (QClatparm_a / dsqrt(2.D0))*self%dsm(1:5,1:5) 

      self%rsm(1,1:5) = (/c(1), s(1), c2(1), s2(1), 0.D0/)
      self%rsm(2,1:5) = (/c(2), s(2), c2(2), s2(2), 0.D0/)
      self%rsm(3,1:5) = (/c(3), s(3), c2(3), s2(3), 0.D0/)
      self%rsm(4,1:5) = (/c(4), s(4), c2(4), s2(4), 0.D0/)
      self%rsm(5,1:5) = (/0.D0, 0.D0, 0.D0, 0.D0, dsqrt(2.D0)*QClatparm_a/QClatparm_c/) 

      self%rsm(1:5,1:5) = (1.D0 / QClatparm_a / dsqrt(2.D0)) * self%rsm(1:5,1:5) 

      self%dsm = transpose(self%dsm)
      self%rsm = transpose(self%rsm)

      A   = (QClatparm_a / dsqrt(2.D0)) ** 2
      B   = QClatparm_c ** 2

      self%dmt(1,1:5) = (/A, 0.D0, 0.D0, 0.D0, 0.D0/)
      self%dmt(2,1:5) = (/0.D0, A, 0.D0, 0.D0, 0.D0/)
      self%dmt(3,1:5) = (/0.D0, 0.D0, A, 0.D0, 0.D0/)
      self%dmt(4,1:5) = (/0.D0, 0.D0, 0.D0, A, 0.D0/)
      self%dmt(5,1:5) = (/0.D0, 0.D0, 0.D0, 0.D0, B/)

      self%vol  = sqrt(CalcDeterminant(self%dmt, 5, 5))

      A   = 2.D0 * (1.D0 / QClatparm_a) ** 2
      B   = (1.D0 / QClatparm_c) ** 2

      self%rmt(1,1:5) = (/A, 0.D0, 0.D0, 0.D0, 0.D0/)
      self%rmt(2,1:5) = (/0.D0, A, 0.D0, 0.D0, 0.D0/)
      self%rmt(3,1:5) = (/0.D0, 0.D0, A, 0.D0, 0.D0/)
      self%rmt(4,1:5) = (/0.D0, 0.D0, 0.D0, A, 0.D0/)
      self%rmt(5,1:5) = (/0.D0, 0.D0, 0.D0, 0.D0, B/)
    end select
end if 

if(self%getQCtype() .eq. 'Ico') then
  do i = 2,6
    c(i-1)  = dcos(2.D0*cPi*dble(i)/5.D0)
    s(i-1)  = dsin(2.D0*cPi*dble(i)/5.D0)
    c2(i-1) = dcos(4.D0*cPi*dble(i)/5.D0)
    s2(i-1) = dsin(4.D0*cPi*dble(i)/5.D0)
  end do  

  select type (self)
    class is (QCcell_icosahedral_T)
      self%scaling(1,1:6) = (/1.D0, 1.D0, 1.D0, 1.D0, 1.D0, 1.D0/)
      self%scaling(2,1:6) = (/1.D0, 1.D0, 1.D0, -1.D0, -1.D0, 1.D0/)
      self%scaling(3,1:6) = (/1.D0, 1.D0, 1.D0, 1.D0, -1.D0, -1.D0/)
      self%scaling(4,1:6) = (/1.D0, -1.D0, 1.D0, 1.D0, 1.D0, -1.D0/)
      self%scaling(5,1:6) = (/1.D0, -1.D0, -1.D0, 1.D0, 1.D0, 1.D0/)
      self%scaling(6,1:6) = (/1.D0, 1.D0, -1.D0, -1.D0, 1.D0, 1.D0/)

      self%scaling      = 0.5D0 * self%scaling

      self%dsm(1,1:6) = (/0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 1.D0/)
      self%dsm(2,1:6) = (/c(1)*st, s(1)*st, ct, c2(1)*st, s2(1)*st, -ct/)
      self%dsm(3,1:6) = (/c(2)*st, s(2)*st, ct, c2(2)*st, s2(2)*st, -ct/)
      self%dsm(4,1:6) = (/c(3)*st, s(3)*st, ct, c2(3)*st, s2(3)*st, -ct/)
      self%dsm(5,1:6) = (/c(4)*st, s(4)*st, ct, c2(4)*st, s2(4)*st, -ct/)
      self%dsm(6,1:6) = (/c(5)*st, s(5)*st, ct, c2(5)*st, s2(5)*st, -ct/)

      self%dsm = self%dsm * QClatparm

      self%rsm(1,1:6) = (/0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 1.D0/)
      self%rsm(2,1:6) = (/c(1)*st, s(1)*st, ct, c2(1)*st, s2(1)*st, -ct/)
      self%rsm(3,1:6) = (/c(2)*st, s(2)*st, ct, c2(2)*st, s2(2)*st, -ct/)
      self%rsm(4,1:6) = (/c(3)*st, s(3)*st, ct, c2(3)*st, s2(3)*st, -ct/)
      self%rsm(5,1:6) = (/c(4)*st, s(4)*st, ct, c2(4)*st, s2(4)*st, -ct/)
      self%rsm(6,1:6) = (/c(5)*st, s(5)*st, ct, c2(5)*st, s2(5)*st, -ct/)

      self%rsm = ( self%rsm / (QClatparm * 2.D0) ) 

      self%dsm = transpose(self%dsm)
      self%rsm = transpose(self%rsm)

      A   = QClatparm ** 2

      self%dmt(1,1:6) = (/A, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0/)
      self%dmt(2,1:6) = (/0.D0, A, 0.D0, 0.D0, 0.D0, 0.D0/)
      self%dmt(3,1:6) = (/0.D0, 0.D0, A, 0.D0, 0.D0, 0.D0/)
      self%dmt(4,1:6) = (/0.D0, 0.D0, 0.D0, A, 0.D0, 0.D0/)
      self%dmt(5,1:6) = (/0.D0, 0.D0, 0.D0, 0.D0, A, 0.D0/)
      self%dmt(6,1:6) = (/0.D0, 0.D0, 0.D0, 0.D0, 0.D0, A/)

      self%vol  = sqrt(CalcDeterminant(self%dmt, 6, 6))

      A   = (1.D0 / QClatparm) ** 2

      self%rmt(1,1:6) = (/A, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0/)
      self%rmt(2,1:6) = (/0.D0, A, 0.D0, 0.D0, 0.D0, 0.D0/)
      self%rmt(3,1:6) = (/0.D0, 0.D0, A, 0.D0, 0.D0, 0.D0/)
      self%rmt(4,1:6) = (/0.D0, 0.D0, 0.D0, A, 0.D0, 0.D0/)
      self%rmt(5,1:6) = (/0.D0, 0.D0, 0.D0, 0.D0, A, 0.D0/)
      self%rmt(6,1:6) = (/0.D0, 0.D0, 0.D0, 0.D0, 0.D0, A/)
    end select
end if 

end subroutine setMetricParametersQC_

!--------------------------------------------------------------------------
recursive subroutine TransSpaceQC5_(self, t, d, inspace, outspace)
!DEC$ ATTRIBUTES DLLEXPORT ::TransSpaceQC5_
  !! author: MDG/SS
  !! version: 1.0
  !! date: 02/01/22
  !!
  !! convert between different vector spaces 'd', 'r', 'c'

IMPLICIT NONE

class(QCcell_axial_T),intent(inout)   :: self
real(kind=dbl),INTENT(IN)             :: t(5)
real(kind=dbl),INTENT(OUT)            :: d(5)
character(1),INTENT(IN)               :: inspace
character(1),INTENT(IN)               :: outspace

! intercept the case where inspace and outspace are the same 
if (inspace.eq.outspace) then
  d = t
  return
end if

if (inspace.eq.'d') then
! direct to Cartesian (pre-multiplication)
  if (outspace.eq.'c') then
   d = matmul(self%dsm,t)
   return
  end if
! direct to reciprocal (post-multiplication)
  if (outspace.eq.'r') then
   d = matmul(t,self%dmt)
   return
  end if
 end if

 if (inspace.eq.'r') then
! reciprocal to Cartesian (pre-multiplication)
  if (outspace.eq.'c') then
   d = matmul(self%rsm,t)
   return
  end if
! reciprocal to direct (post-multiplication)
  if (outspace.eq.'d') then
   d = matmul(t,self%rmt)
   return
  end if
 end if

 if (inspace.eq.'c') then
! Cartesian to direct (post-multiplication)
  if (outspace.eq.'d') then
   d = matmul(t,self%rsm)
   return
  end if
! Cartesian to reciprocal (post-multiplication)
  if (outspace.eq.'r') then
   d = matmul(t,self%dsm)
   return
  end if
 end if

end subroutine TransSpaceQC5_

!--------------------------------------------------------------------------
recursive subroutine TransSpaceQC6_(self, t, d, inspace, outspace)
!DEC$ ATTRIBUTES DLLEXPORT ::TransSpaceQC6_
  !! author: MDG/SS
  !! version: 1.0
  !! date: 02/01/22
  !!
  !! convert between different vector spaces 'd', 'r', 'c'

IMPLICIT NONE

class(QCcell_icosahedral_T),intent(inout)   :: self
real(kind=dbl),INTENT(IN)                   :: t(6)
real(kind=dbl),INTENT(OUT)                  :: d(6)
character(1),INTENT(IN)                     :: inspace
character(1),INTENT(IN)                     :: outspace

! intercept the case where inspace and outspace are the same 
if (inspace.eq.outspace) then
  d = t
  return
end if

if (inspace.eq.'d') then
! direct to Cartesian (pre-multiplication)
  if (outspace.eq.'c') then
   d = matmul(self%dsm,t)
   return
  end if
! direct to reciprocal (post-multiplication)
  if (outspace.eq.'r') then
   d = matmul(t,self%dmt)
   return
  end if
 end if

 if (inspace.eq.'r') then
! reciprocal to Cartesian (pre-multiplication)
  if (outspace.eq.'c') then
   d = matmul(self%rsm,t)
   return
  end if
! reciprocal to direct (post-multiplication)
  if (outspace.eq.'d') then
   d = matmul(t,self%rmt)
   return
  end if
 end if

 if (inspace.eq.'c') then
! Cartesian to direct (post-multiplication)
  if (outspace.eq.'d') then
   d = matmul(t,self%rsm)
   return
  end if
! Cartesian to reciprocal (post-multiplication)
  if (outspace.eq.'r') then
   d = matmul(t,self%dsm)
   return
  end if
 end if

end subroutine TransSpaceQC6_

!--------------------------------------------------------------------------
recursive function getGvectorQC5_(self, QCindex, OP) result(gvector)
!DEC$ ATTRIBUTES DLLEXPORT ::getGvectorQC5_ 
  !! author: MDG/SS
  !! version: 1.0
  !! date: 02/01/22
  !!
  !! Generate the parallel/orthogonal reciprocal lattice vector for a set of indices

use mod_io

IMPLICIT NONE

class(QCcell_axial_T),INTENT(INOUT)    :: self
real(kind=dbl),INTENT(IN)              :: QCindex(5)
character(1),INTENT(IN)                :: OP

real(kind=dbl)                         :: gvector(3), dvec(5)

type(IO_T)                             :: Message

call self%TransSpace(QCindex,dvec,'r','c')
if (OP .eq. 'O') then
  gvector = (/dvec(3), dvec(4), 0.D0/)
else if (OP .eq. 'P') then
  gvector  = (/dvec(1), dvec(2), dvec(5)/)
else
  call Message%printError('getGvectorQC5','Unknown OP parameter passed to function')
end if

end function getGvectorQC5_

!--------------------------------------------------------------------------
recursive function getGvectorQC6_(self, QCindex, OP) result(gvector)
!DEC$ ATTRIBUTES DLLEXPORT ::getGvectorQC6_ 
  !! author: MDG/SS
  !! version: 1.0
  !! date: 02/01/22
  !!
  !! Generate the parallel/orthogonal reciprocal lattice vector for a set of indices

use mod_io

IMPLICIT NONE

class(QCcell_icosahedral_T),INTENT(INOUT) :: self
real(kind=dbl),INTENT(IN)                 :: QCindex(6)
character(1),INTENT(IN)                   :: OP

real(kind=dbl)                            :: gvector(3), dvec(6)

type(IO_T)                                :: Message

call self%TransSpace(QCindex,dvec,'r','c')
if (OP .eq. 'O') then
  gvector = (/dvec(3), dvec(4), dvec(6)/)
else if (OP .eq. 'P') then
  gvector  = (/dvec(1), dvec(2), dvec(3)/)
else
  call Message%printError('getGvectorQC6','Unknown OP parameter passed to function')
end if

end function getGvectorQC6_

!--------------------------------------------------------------------------
recursive function getvectorLengthQC5_(self, QCindex, OP, space) result(gl)
!DEC$ ATTRIBUTES DLLEXPORT ::getvectorLengthQC5_
 !! author: MDG/SS
 !! version: 1.0
 !! date: 02/01/22
 !!
 !! Get the length of a parallel/orthogonal reciprocal lattice vector for a set of indices

use mod_io

IMPLICIT NONE

class(QCcell_axial_T),INTENT(INOUT)   :: self
integer(kind=irg),INTENT(IN)          :: QCindex(5)
character(1),INTENT(IN)               :: OP
character(1),INTENT(IN)               :: space  ! 'r' or 'd'

type(IO_T)                            :: Message

real(kind=dbl)                        :: gl
real(kind=dbl)                        :: dvec(5), dvec2(5)

dvec = dble(QCindex)

if(space .eq. 'r') then
  call self%TransSpace(dvec,dvec2,'r','c')
  if (OP.eq.'O') then
    !gl = dsqrt(DOT_PRODUCT(dvec, matmul(QCcell%rmto, dvec)))
    gl = sqrt(dvec2(3)**2 + dvec2(4)**2)
  else if (OP .eq. 'P') then
    !gl = dsqrt(DOT_PRODUCT(dvec, matmul(QCcell%rmtp, dvec)))
    gl = sqrt(dvec2(1)**2 + dvec2(2)**2 + dvec2(5)**2)
  else
    call Message%printError('QC_getvectorLength : ','Unknown OP parameter passed to function')
  end if
else if(space .eq. 'd') then
  !dvec2 = matmul(dvec,QCcell%Mdirect)
  call self%TransSpace(dvec,dvec2,'d','c')
  if (OP.eq.'O') then
    !gl = dsqrt(DOT_PRODUCT(dvec, matmul(QCcell%dmto, dvec)))
    gl = sqrt(dvec2(3)**2 + dvec2(4)**2)
  else if (OP .eq. 'P') then
    !gl = dsqrt(DOT_PRODUCT(dvec, matmul(QCcell%dmtp, dvec)))
    gl = sqrt(dvec2(1)**2 + dvec2(2)**2 + dvec2(5)**2)
  else
    call Message%printError('QC_getvectorLength : ','Unknown OP parameter passed to function')
  end if
else
  call Message%printError('QC_getvectorLength: ', 'Unkown space input.')
end if

end function getvectorLengthQC5_

!--------------------------------------------------------------------------
recursive function getvectorLengthQC6_(self, QCindex, OP, space) result(gl)
!DEC$ ATTRIBUTES DLLEXPORT ::getvectorLengthQC6_
 !! author: MDG/SS
 !! version: 1.0
 !! date: 02/01/22
 !!
 !! Get the length of a parallel/orthogonal reciprocal lattice vector for a set of indices

use mod_io

IMPLICIT NONE

class(QCcell_icosahedral_T),INTENT(INOUT)   :: self
integer(kind=irg),INTENT(IN)                :: QCindex(6)
character(1),INTENT(IN)                     :: OP
character(1),INTENT(IN)                     :: space  ! 'r' or 'd'

type(IO_T)                                  :: Message

real(kind=dbl)                              :: gl
real(kind=dbl)                              :: dvec(6), dvec2(6)

dvec = dble(QCindex)

if(space .eq. 'r') then
  call self%TransSpace(dvec,dvec2,'r','c')
  if (OP.eq.'O') then
    gl = sqrt(dvec2(4)**2 + dvec2(5)**2 + dvec2(6)**2)
  else if (OP .eq. 'P') then
    gl = sqrt(dvec2(1)**2 + dvec2(2)**2 + dvec2(3)**2)
  else
    call Message%printError('getvectorLengthQC6 : ','Unknown OP parameter passed to function')
  end if
else if(space .eq. 'd') then
  call self%TransSpace(dvec,dvec2,'d','c')
  if (OP.eq.'O') then
    gl = sqrt(dvec2(4)**2 + dvec2(5)**2 + dvec2(6)**2)
  else if (OP .eq. 'P') then
    gl = sqrt(dvec2(1)**2 + dvec2(2)**2 + dvec2(3)**2)
  else
    call Message%printError('getvectorLengthQC6 : ','Unknown OP parameter passed to function')
  end if
else
  call Message%printError('getvectorLengthQC6 : ', 'Unkown space input.')
end if

end function getvectorLengthQC6_

!--------------------------------------------------------------------------
recursive function ShapeTransformTriangle_(self, ap, g) result(stt)
!DEC$ ATTRIBUTES DLLEXPORT :: ShapeTransformTriangle_
 !! author: MDG/SS
 !! version: 1.0
 !! date: 02/01/22
 !!
 !! calculate shape transform of the a triangle, which
 !! is the building block for any polygon occupation domain. this result is known
 !! analytically. we will simply sum it up for the 24 triangles. the
 !! limits of the expression was calculated using mathematica

use mod_io

IMPLICIT NONE

class(QCcell_axial_T),INTENT(INOUT)   :: self
real(kind=dbl),INTENT(IN)             :: g(5)
real(kind=dbl),INTENT(IN)             :: ap

type(IO_T)                            :: Message 

complex(kind=dbl)                     :: stt, imagj
real(kind=dbl)                        :: e1(5), e2(5), etmp(5), vlen, dvec1(2), dvec2(2), dvec3(2)
real(kind=dbl)                        :: a1, a2, a3, Ar, prod, t
complex(kind=dbl)                     :: t1, t2
integer(kind=irg)                     :: icase
real(kind=dbl),parameter              :: eps = 1.0D-12
real(kind=dbl),parameter              :: tau = (1.D0 + dsqrt(5.D0))/2.D0

imagj = cmplx(0.D0, 1.D0)

if(self%getQCtype() .eq. 'DoD') then

  !e1    = (/0.D0, -1.D0, 0.D0, 1.D0, 0.D0/)
  e1 = (/1.D0, -1.D0, 0.D0, 1.D0, 0.D0/)
  call self%TransSpace(e1, etmp, 'd', 'c')
  dvec1 = (/etmp(3), etmp(4)/)

  !e2    = (/0.D0, 0.D0, 0.D0, 1.D0, 0.D0/)!/dsqrt(3.D0)
  e2    = (/0.D0, -1.D0, -1.D0, 1.D0, 0.D0/)
  call self%TransSpace(e2, etmp, 'd', 'c')
  dvec2 = (/etmp(3), etmp(4)/)

  Ar    = ap**2 * dsin(cPi/6.D0)

else if(self%getQCtype() .eq. 'Dec') then

  !e1    = (/1.D0, 0.D0, 0.D0, 1.D0, 0.D0/)
  e1    = (/1.D0, 0.D0, 0.D0, 0.D0, 0.D0/) - (/1.D0, 1.D0, 1.D0, 1.D0, 0.D0/)/5.D0
  call self%TransSpace(e1, etmp, 'd', 'c')
  dvec1 = (/etmp(3), etmp(4)/)

  e2    = (/0.D0, 0.D0, 1.D0, 0.D0, 0.D0/) - (/1.D0, 1.D0, 1.D0, 1.D0, 0.D0/)/5.D0
  !e2    = (/0.D0, 1.D0, 0.D0, 1.D0, 0.D0/) - (/1.D0, 1.D0, 1.D0, 1.D0, 0.D0/)/5.D0
  call self%TransSpace(e2, etmp, 'd', 'c')
  dvec2 = (/etmp(3), etmp(4)/)
  
  Ar    = ap**2 * dsin(cPi/5.D0)

else if(self%getQCtype() .eq. 'Oct') then

  e1    = (/1.D0, -1.D0, 0.D0, 1.D0, 0.D0/)/2.D0
  call self%TransSpace(e1, etmp, 'd', 'c')
  dvec1 = (/etmp(3), etmp(4)/)

  e2    = (/1.D0, -1.D0, -1.D0, 1.D0, 0.D0/)/2.D0
  call self%TransSpace(e2, etmp, 'd', 'c')
  dvec2 = (/etmp(3), etmp(4)/)

  Ar    = ap**2 * dsin(cPi/4.D0)

else
  call Message%printError('ShapeTransformPolygon','unknown symmetry type.')
end if

call self%TransSpace(g, etmp, 'r', 'c')
dvec3 = (/etmp(3), etmp(4)/)

a1  = 2.D0 * cPi * DOT_PRODUCT(dvec1, dvec3)
a2  = 2.D0 * cPi * DOT_PRODUCT(dvec2, dvec3)

a3 = a1 - a2

icase = 0

! (a1 -> 0) or (a2 -> 0) or (a1 -> a2 ~-> 0) or (a1 -> a2 -> 0) 
if((abs(a1) .lt. eps) .and. (abs(a2) .gt. eps) .and. (abs(a3) .gt. eps)) then
  icase = 1
else if((abs(a1) .gt. eps) .and. (abs(a2) .lt. eps) .and. (abs(a3) .gt. eps)) then
  icase = 2
else if((abs(a1) .gt. eps) .and. (abs(a2) .gt. eps) .and. (abs(a3) .lt. eps)) then
  icase = 3
else if((abs(a1) .lt. eps) .and. (abs(a2) .lt. eps)) then
  icase = 4
end if

select case (icase)
  case(1)
    t1    = 1 + imagj * a2 - exp(imagj * a2) 
    prod  = a2 * a2
    stt   = Ar * t1 / prod 

  case(2)
    t1    = 1 + imagj * a1 - exp(imagj * a1) 
    prod  = a1 * a1
    stt   = Ar * t1 / prod

  case(3)
    t1    = imagj * a2 * exp(imagj * a2) - exp(imagj * a2) + cmplx(1.D0, 0.D0)
    prod  = a2 * a2
    stt   = -Ar * t1 / prod

  case(4)
    stt = cmplx(Ar/2, 0.D0)

  case DEFAULT
    t1    = a2 * exp(imagj * a1)
    t2    = a1 * exp(imagj * a2)
    prod  = a1 * a2 * a3
    stt   = -Ar * (t1 - t2 + a3)/prod
end select

end function ShapeTransformTriangle_

!--------------------------------------------------------------------------
recursive function ShapeTransformPolygonCa_(self, QCSG, hkl, asite) result(stp)
!DEC$ ATTRIBUTES DLLEXPORT :: ShapeTransformPolygonCa_
 !! author: MDG/SS
 !! version: 1.0
 !! date: 02/01/22
 !!
 !! calculate shape transform of the polygon Ca (see Franz gahler, 
 !! Crystallography of Dodecagonal Quasicrystals)

use mod_io
use mod_QCsymmetry

IMPLICIT NONE

class(QCcell_axial_T),INTENT(INOUT)   :: self
type(QCspacegroup_T),INTENT(IN)       :: QCSG
integer(kind=irg),INTENT(IN)          :: hkl(5)
integer(kind=irg),INTENT(IN)          :: asite

type(IO_T)                            :: Message 

real(kind=dbl)                        :: ar
real(kind=dbl),parameter              :: tau = (1.D0 + dsqrt(5.D0))/2.D0
complex(kind=dbl)                     :: stp
integer(kind=irg)                     :: Pmdims, ii
real(kind=dbl)                        :: A_Polygon, hkl2(5), mat(5,5), mat2(5,5)

stp = cmplx(0.D0,0.D0)

if(self%getQCtype() .eq. 'DoD') then
  ar        = self%ATOM_pos(asite, 10) * dsqrt(2.D0/3.D0) * (self%QClatparm_a)
  Pmdims    = 24

else if(self%getQCtype() .eq. 'Dec') then
  ar        = self%ATOM_pos(asite, 10)  * self%QClatparm_a !* 2.D0/ sqrt(5.D0)
  Pmdims    = 5

else if(self%getQCtype() .eq. 'Oct') then
  ar        = self%ATOM_pos(asite, 10) * self%QClatparm_a/dsqrt(2.D0)
  Pmdims    = 16 

else
  call Message%printError('ShapeTransformPolygonCa_','unknown symmetry type.')

end if

do ii = 1,Pmdims
  hkl2 = matmul(QCSG%icos(:,:,ii),dble(hkl))
  stp  = stp + self%ShapeTransformTriangle(ar, hkl2)
end do

!stp = stp / A_Polygon

end function ShapeTransformPolygonCa_

!--------------------------------------------------------------------------
recursive function ShapeTransformPyramid_(self, ap, g) result(stp)
!DEC$ ATTRIBUTES DLLEXPORT :: ShapeTransformPyramid_
 !! author: MDG/SS
 !! version: 1.0
 !! date: 02/01/22
 !!
 !! calculate shape transform of the a triangular pyramid, which
 !! is the building block for the Triacontahedron. this result is known
 !! analytically. we will simply sum it up for the 120 pyramids. the
 !! limits of the expression was calculated using mathematica

use mod_io

IMPLICIT NONE

class(QCcell_icosahedral_T),INTENT(INOUT)   :: self
real(kind=dbl),INTENT(IN)                   :: g(6)
real(kind=dbl),INTENT(IN)                   :: ap
complex(kind=dbl)                           :: stp

type(IO_T)                                  :: Message

real(kind=dbl)                              :: a1, a2, a3, a4, a5, a6, prod, Vr
complex(kind=dbl)                           :: t1, t2, t3, t4, imagj
real(kind=dbl)                              :: e1(6), e2(6), e3(6), etmp(6)              ! basis vectors
real(kind=dbl)                              :: e1_c(3), e2_c(3), e3_c(3), gperp(3), vlen
real(kind=dbl),parameter                    :: eps = 1.0D-12
integer(kind=irg)                           :: icase

imagj = cmplx(0.D0, 1.D0)

! edges of tetrahedron in perp space
e1    = (/1.D0, -1.D0, -1.D0, -1.D0, -1.D0, -1.D0/) * 0.5D0 
call self%TransSpace(e1, etmp, 'd', 'c')
e1_c  = etmp(4:6) 

e2    =  (/1.D0, 1.D0, -1.D0, -1.D0, -1.D0, -1.D0/) * 0.5D0 
call self%TransSpace(e2, etmp, 'd', 'c')
e2_c  =  etmp(4:6) 

e3    =  (/1.D0, 0.D0, -1.D0, -1.D0, 0.D0, -1.D0/) * 0.5D0 
call self%TransSpace(e3, etmp, 'd', 'c')
e3_c  =  etmp(4:6) 

gperp =  self%getGvector(g, 'O')

a1    =  2.D0 * cPi * DOT_PRODUCT(e1_c, gperp)
a2    =  2.D0 * cPi * DOT_PRODUCT(e2_c, gperp)
a3    =  2.D0 * cPi * DOT_PRODUCT(e3_c, gperp)

a4    =  a2 - a3
a5    =  a3 - a1
a6    =  a1 - a2

! volume of tetrahedron
Vr    =  ap**3 * DOT_PRODUCT(e1_c,(/e2_c(2)*e3_c(3) - e2_c(3)*e3_c(2), e2_c(3)*e3_c(1) - e2_c(1)*e3_c(3), &
                          e2_c(1)*e3_c(2) - e2_c(2)*e3_c(1)/))

! define different cases with limits (icase variable)
! the limits were calculated in Mathematica

icase = 0

! case 1: a1 -> 0, a2 -> 0, a3 -> 0
if((abs(a1) .lt. eps) .and. (abs(a2) .lt. eps) .and. (abs(a3) .lt. eps)) then
  icase = 1
! case 2: a1 -> 0, a2 -> 0
else if((abs(a1) .lt. eps) .and. (abs(a2) .lt. eps) .and. (abs(a3) .gt. eps)) then
  icase = 2
! case3: a1 -> 0, a3 -> 0
else if((abs(a1) .lt. eps) .and. (abs(a3) .lt. eps) .and. (abs(a2) .gt. eps)) then
  icase = 3
! case 4: a2 -> 0, a3 -> 0
else if((abs(a2) .lt. eps) .and. (abs(a3) .lt. eps) .and. (abs(a1) .gt. eps)) then
  icase = 4
! case 5/6: a1 -> 0
else if((abs(a1) .lt. eps) .and. (abs(a2) .gt. eps) .and. (abs(a3) .gt. eps)) then

  if(abs(a4) .lt. eps) then
    icase = 5
  else
    icase = 6
  end if

! case 7/8: a2 -> 0
else if((abs(a2) .lt. eps) .and. (abs(a1) .gt. eps) .and. (abs(a3) .gt. eps)) then

  if(abs(a5) .lt. eps) then
    icase = 7
  else
    icase = 8
  end if

! case 9/10: a3 -> 0
else if((abs(a3) .lt. eps) .and. (abs(a2) .gt. eps) .and. (abs(a1) .gt. eps)) then

  if(abs(a6) .lt. eps) then
    icase = 9
  else
    icase = 10
  end if

! case 11/14: a3 -> 0
else if((abs(a1) .gt. eps) .and. (abs(a2) .gt. eps) .and. (abs(a3) .gt. eps)) then

  if((abs(a4) .lt. eps) .and. (abs(a5) .gt. eps) .and. (abs(a6) .gt. eps)) then
    icase = 11
  else if((abs(a5) .lt. eps) .and. (abs(a4) .gt. eps) .and. (abs(a6) .gt. eps)) then
    icase = 12
  else if((abs(a6) .lt. eps) .and. (abs(a4) .gt. eps) .and. (abs(a5) .gt. eps)) then
    icase = 13
  else if((abs(a4) .lt. eps) .and. (abs(a5) .lt. eps)) then
    icase = 14
  end if
end if

select case (icase)
  case(1)
    stp = cmplx(Vr/6.D0, 0.D0)
    
  case(2)
    t1    = 2.D0 * exp(imagj * a3)
    t2    = a3 * (2.D0 * imagj - a3)
    t3    = 2.D0 * a3**3.D0
    stp   = -imagj * Vr * (2.D0 - t1 + t2)/t3
    
  case(3)
    t1    = 2.D0 * exp(imagj * a2)
    t2    = a2 * (2.D0 * imagj - a2)
    t3    = 2.D0 * a2**3.D0
    stp   = -imagj * Vr * (2.D0 - t1 + t2)/t3
    
  case(4)
    t1    = 2.D0 * exp(imagj * a1)
    t2    = a1 * (2.D0 * imagj - a1)
    t3    = 2.D0 * a1**3.D0
    stp   = -imagj * Vr * (2.D0 - t1 + t2)/t3
    
  case(5)
    t1    = (2.D0 - imagj * a3) * exp(imagj * a3)
    t2    = -2.D0 - imagj * a3
    t3    = a3**3
    stp   = -imagj * Vr * (t1 + t2)/t3
    
  case(6)
    t1    = (-1.D0 + exp(imagj * a3) - imagj*a3) * a2**2
    t2    = (1.D0 - exp(imagj * a2)) * a3**2
    t3    = imagj * a2 * a3**2
    t4    = a4 * (a2 * a3)**2
    stp   = -imagj * Vr * (t1 + t2 + t3)/t4
    
  case(7)
    t1    = (2.D0 - imagj * a1) * exp(imagj * a1)
    t2    = -2.D0 - imagj * a1
    t3    = a1**3
    stp   = -imagj * Vr * (t1 + t2)/t3
    
  case(8)
    t1    = (-1.D0 + exp(imagj * a3) - imagj*a3) * a1**2
    t2    = (1.D0 - exp(imagj * a1)) * a3**2
    t3    = imagj * a1 * a3**2
    t4    = -a5 * (a1 * a3)**2
    stp   = -imagj * Vr * (t1 + t2 + t3)/t4
    
  case(9)
    t1    = (2.D0 - imagj * a2) * exp(imagj * a2)
    t2    = -2.D0 - imagj * a2
    t3    = a2**3
    stp   = -imagj * Vr * (t1 + t2)/t3
    
  case(10)
    t1    = (-1.D0 + exp(imagj * a2) - imagj*a2) * a1**2
    t2    = (1.D0 - exp(imagj * a1)) * a2**2
    t3    = imagj * a1 * a2**2
    t4    = a6 * (a1 * a2)**2
    stp   = -imagj * Vr * (t1 + t2 + t3)/t4
    
  case(11)
    t1    = a3 * a1 * (-2.D0 + exp(imagj * a3) * (2.D0 - imagj * a3))
    t2    = (1.D0 - exp(imagj * a1)) * a3**2
    t3    = (1.D0 + imagj * exp(imagj * a3) * (imagj + a3)) * a1**2
    t4    = a1 * (a3 * a5)**2
    stp   = -imagj * Vr * (t1 + t2 + t3)/t4
    
  case(12)
    t1    = a3 * a2 * (-2.D0 + exp(imagj * a3) * (2.D0 - imagj * a3))
    t2    = (1.D0 - exp(imagj * a2)) * a3**2
    t3    = (1.D0 + imagj * exp(imagj * a3) * (imagj + a3)) * a2**2
    t4    = a2 * (a3 * a4)**2
    stp   = -imagj * Vr * (t1 + t2 + t3)/t4

  case(13)
    t1    = a3 * a2 * (-2.D0 + exp(imagj * a2) * (2.D0 + imagj * a3))
    t2    = (1.D0 - exp(imagj * a2)) * a3**2
    t3    = (1.D0 - imagj * exp(imagj * a2) * a3 - exp(imagj * a3) ) * a2**2
    t4    = a1 * (a3 * a5)**2
    stp   = -imagj * Vr * (t1 + t2 + t3)/t4

  case(14) 
    t1    =  (-2.D0 + a2 * (2.D0 * imagj + a2)) * exp(imagj * a2)
    t2    =  2.D0 * a2**3
    stp   =  -imagj * Vr * (2.D0 + t1)/t2

  case DEFAULT
    t1    =  a2 * a3 * a4 * exp(imagj * a1)
    t2    =  a1 * a3 * a5 * exp(imagj * a2)
    t3    =  a1 * a2 * a6 * exp(imagj * a3)
    t4    =  a4 * a5 * a6
    prod  =  a1 * a2 * a3 * a4 * a5 * a6
    stp   =  -imagj * Vr * (t1 + t2 + t3 + t4) / prod
end select 

end function ShapeTransformPyramid_

!--------------------------------------------------------------------------
recursive function ShapeTransformTriacontahedron_(self, QCSG, hkl, asite) result(stt)
!DEC$ ATTRIBUTES DLLEXPORT :: ShapeTransformTriacontahedron_

use mod_QCsymmetry

IMPLICIT NONE

class(QCcell_icosahedral_T),INTENT(INOUT)   :: self
type(QCspacegroup_T),INTENT(INOUT)          :: QCSG
integer(kind=irg),INTENT(IN)                :: hkl(6)
integer(kind=irg),INTENT(IN)                :: asite
complex(kind=dbl)                           :: stt

integer(kind=irg)                           :: Pmdims, ii, jj, nn
real(kind=dbl)                              :: V_Triacontahedron, hkl2(6), ar


ar = self%ATOM_pos(asite,10) * self%QClatparm

stt = cmplx(0.D0,0.D0)

Pmdims = QCSG%getNUMpt()

do ii = 1,Pmdims
  hkl2 = matmul(QCSG%direc(ii,:,:),dble(hkl))
  stt  = stt + self%ShapeTransformPyramid(ar, hkl2)
end do

end function ShapeTransformTriacontahedron_

!--------------------------------------------------------------------------
recursive subroutine displayPeriodicTable(self)
!DEC$ ATTRIBUTES DLLEXPORT :: displayPeriodicTable
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! display the periodic table so that the user can look up the atomic number

use mod_io

IMPLICIT NONE

class(QCcell_T),intent(inout)           :: self
type (IO_T)                             :: Message

call Message%printMessage( (/ &
   ' ------------------------------------ Periodic Table of the Elements --------------------------------------', &
   '1:H                                                                                                    2:He', &
   '3:Li  4:Be                                                               5:B   6:C   7:N   8:O   9:F  10:Ne', &
   '11:Na 12:Mg                                                             13:Al 14:Si 15:P  16:S  17:Cl 18:Ar', &
   '19:K  20:Ca 21:Sc 22:Ti 23:V  24:Cr 25:Mn 26:Fe 27:Co 28:Ni 29:Cu 30:Zn 31:Ga 32:Ge 33:As 34:Se 35:Br 36:Kr', &
   '37:Rb 38:Sr 39:Y  40:Zr 41:Nb 42:Mo 43:Tc 44:Ru 45:Rh 46:Pd 47:Ag 48:Cd 49:In 50:Sn 51:Sb 52:Te 53: I 54:Xe', &
   '55:Cs 56:Ba ----- 72:Hf 73:Ta 74:W  75:Re 76:Os 77:Ir 78:Pt 79:Au 80:Hg 81:Tl 82:Pb 83:Bi 84:Po 85:At 86:Rn', &
   '87:Fr 88:Ra 89:Ac 90:Th 91:Pa 92:U                                                                         ', &
   '57:La 58:Ce 59:Pr 60:Nd 61:Pm 62:Sm 63:Eu 64:Gd 65:Tb 66:Dy 67:Ho 68:Er 69:Tm 70:Yb 71:Lu                  ', &
   ' ----------------------------------------------------------------------------------------------------------' /) )

end subroutine displayPeriodicTable

!--------------------------------------------------------------------------
recursive subroutine QCextractposition(list,pt,iQC)
!DEC$ ATTRIBUTES DLLEXPORT :: QCextractposition
!! author: MDG/SS 
!! version: 1.0 
!! date: 01/30/22
!!
!! extract quasicrystal atom position data from a string

IMPLICIT NONE

character(1),INTENT(IN)                 :: list(256)                              !< input string
real(kind=sgl),INTENT(OUT)              :: pt(10)                                !< output real array
logical,INTENT(IN),OPTIONAL             :: iQC

integer(kind=irg)                       :: comma(12),slash(12),period(12), &
                                           ccnt,scnt,pcnt,pp,i,j,hcnt, &
                                           ip,ipt,icnt,nd,n,k,ns,nocc                !< auxiliary variables
integer(kind=irg),parameter             :: nmb(48:57)=(/0,1,2,3,4,5,6,7,8,9/)   !< list of numbers
real(kind=dbl)                          :: nominator,denominator,x              !< used for fraction interpretation
logical                                 :: hasperiod                            !< used for decimal interpretation

nocc = 6
if(present(iQC)) then
  if(iQC) then
    nocc = 7
  end if
end if

! initalize a few variables
 comma(1:6) = 0
 slash(1:5) = 0
 period(1:5) = 0
 ccnt = 0
 scnt = 0
 pcnt = 0
 j = 0
 hcnt = 0
 
! count characters and search for , . and /
 ccnt = ccnt+1
 comma(ccnt) = 0
 do i=1,256
  if (list(i)(1:1).ne.' ') j=j+1
  if (list(i)(1:1).eq.',') then 
   ccnt = ccnt+1
   comma(ccnt)=i
  end if
  if (list(i)(1:1).eq.'/') then 
   scnt = scnt+1
   slash(scnt)=i
  end if
  if (list(i)(1:1).eq.'.') then 
   pcnt = pcnt+1
   period(pcnt)=i
  end if
 end do 
 ccnt = ccnt+1
 comma(ccnt) = j+1
 do while (ccnt.lt.8) 
  ccnt = ccnt+1
  comma(ccnt) = comma(ccnt-1)+1
 end do

! interpret the string
 j = 1
 ip = 1
 icnt = 0
 ipt = 1
 pp = 1
 do i=1,ccnt-1
! is it a real number or a fraction ?
  if (((slash(j).lt.comma(i+1)).and.(scnt.gt.0)).and.(j.le.scnt)) then
! it is a fraction;  get the nominator
   nd = slash(j)-ip
   n = 0
   do k=0,nd-1
    n = 10*n+nmb(ichar(list(ip+k)(1:1)))
   end do
   nominator = dble(n)
   ip = slash(j)+1
! and then the denominator
   nd = comma(i+1)-ip
   n = 0
   do k=0,nd-1
    n = 10*n+nmb(ichar(list(ip+k)(1:1)))
   end do
   denominator = dble(n)
! and fill in the entire range
   pt(ipt) = sngl(nominator/denominator)
   ipt = ipt+1
   ip = comma(i+1)+1
   j=j+1
  else
! no, it is a real number, possibly without a period
! is there a period in this number ?
   if ((period(pp).gt.comma(i)).and.(period(pp).lt.comma(i+1))) then
     hasperiod = .TRUE.
   else
     hasperiod = .FALSE.
   endif
   nd = comma(i+1)-ip
   if (hasperiod) then 
    if (period(pp).eq.comma(i)+1) then
     x = 0.D0
     ns = 2
    else
      x = dble(nmb(ichar(list(ip)(1:1))))
      ns = 3
    end if 
    do k=ns,nd
     x = x + 10.D0**(ns-k-1)*dble(nmb(ichar(list(ip+k-1)(1:1))))
    end do
    pt(ipt)= sngl(x)
    ipt=ipt+1
    ip = comma(i+1)+1
    pp = pp+1
   else
    nd = comma(i+1)-ip
    n = 0
    do k=0,nd-1
     n = 10*n+nmb(ichar(list(ip+k)(1:1)))
    end do
    pt(ipt) = float(n)
    ipt=ipt+1
    ip = comma(i+1)+1
   end if
  end if
 end do 

! set default values
 if (pt(nocc).eq.0.0) pt(nocc) = 1.0

end subroutine QCextractposition

end module mod_QCcrystallography