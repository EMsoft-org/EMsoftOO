! ###################################################################
! Copyright (c) 2014-2023, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_crystallography
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! Everything that has to do with crystallographic computations (used to be crystal.f90 module)

use mod_kinds
use mod_global
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit

IMPLICIT NONE
  private

  type, public :: Cell_T
    !! declaration of the unitcell Cell_T class   (formerly in typedefs.f90 module)
      private

! user provided data (this all goes into the .xtal file)
      real(kind=dbl)                       :: a, b, c, alpha, beta, gamma
       !! lattice parameters [nanometers and degrees]
      integer(kind=irg)                    :: xtal_system
       !! crystal system number; duplicated from the entry in the SpaceGroup class (mod_symmetry.f90)
      integer(kind=irg)                    :: ATOM_type(maxpasym)
       !! atom types in asymetric unit
      integer(kind=irg)                    :: ATOM_ntype
       !! number of atom types in asymmetric unit
      real(kind=sgl)                       :: ATOM_pos(maxpasym,5)
       !! atom coordinates, site occupations, and Debye-Waller factors for atoms in asymmetric unit
      character(fnlen)                     :: xtalname
       !! crystal structure file name
      character(fnlen)                     :: source
       !! string describing the citation for the crystallographic data
! derived metric quantities
      real(kind=dbl)                       :: dmt(3,3)
       !! dirct metric tensor
      real(kind=dbl)                       :: rmt(3,3)
       !! reciprocal metric tensor
      real(kind=dbl)                       :: dsm(3,3)
       !! direct structure matrix
      real(kind=dbl)                       :: rsm(3,3)
       !! reciprocal structure matrix
      real(kind=dbl)                       :: trigmat(3,3)
       !! direct structure matrix for the trigonal/rhombohedral case
      real(kind=dbl)                       :: vol
       !! unit cell volume [nm^3]
! other useful derived quantities
      integer(kind=irg)                    :: numat(maxpasym)
       !! number of atoms of each type in asymmetric unit
      real(kind=dbl),allocatable, public   :: apos(:,:,:)
       !! array with atom coordinates
      real(kind=dbl)                       :: density
       !! theoretical density
      real(kind=dbl)                       :: avZ
       !! average atomic number
      real(kind=dbl)                       :: avA
       !! average atomic number
      logical                              :: Wyckoff = .FALSE.
       !! should we use Wyckoff positions when creating a crystal data file ?

      contains
        ! public
        private
! routines for crystallographic computations
          procedure, pass(self) :: computeMatrices
          procedure, pass(self) :: TransSpaceSingle
          procedure, pass(self) :: TransSpaceDouble
          procedure, pass(self) :: transformCoordinates
          procedure, pass(self) :: calcDotSingle
          procedure, pass(self) :: calcDotDouble
          procedure, pass(self) :: normVecSingle
          procedure, pass(self) :: normVecDouble
          procedure, pass(self) :: calcLengthSingle
          procedure, pass(self) :: calcLengthDouble
          procedure, pass(self) :: calcAngleSingle
          procedure, pass(self) :: calcAngleDouble
          procedure, pass(self) :: calcCrossSingle
          procedure, pass(self) :: calcCrossDouble
          procedure, pass(self) :: calcPositions_
          procedure, pass(self) :: convertfromRtoH_
! routines to set class parameters
          procedure, pass(self) :: getFileName_
          procedure, pass(self) :: getSource_
          procedure, pass(self) :: getVolume_
          procedure, pass(self) :: getNatomtype_
          procedure, pass(self) :: getnumat_
          procedure, pass(self) :: getapos_
          procedure, pass(self) :: getatomtype_
          procedure, pass(self) :: getDensity_
          procedure, pass(self) :: setWyckoff_
          procedure, pass(self) :: setFileName_
          procedure, pass(self) :: setSource_
          procedure, pass(self) :: setnumat_
          procedure, pass(self) :: setXtalSystem_
          procedure, pass(self) :: setAtomPos_
          procedure, pass(self) :: setatomtype_
          procedure, pass(self) :: setAtomPosAll_
          procedure, pass(self) :: setNatomtype_
          procedure, pass(self) :: requestLatticeParameters
          procedure, pass(self) :: setLatticeParameterSingle
          procedure, pass(self) :: getLatticeParameterSingle
          procedure, pass(self) :: getLatticeParametersAll
          procedure, pass(self) :: setLatticeParametersAll
          procedure, pass(self) :: getAsymmetricPosition
          procedure, pass(self) :: getDirectStructureMatrix
          procedure, pass(self) :: getReciprocalStructureMatrix
          procedure, pass(self) :: getDirectMetricTensor
          procedure, pass(self) :: getReciprocalMetricTensor
          procedure, pass(self) :: getAsymPosArray
          ! procedure, pass(self) :: extractAtomPositionData
          procedure, pass(self) :: calcTheoreticalDensity
          procedure, pass(self) :: displayPeriodicTable
! routines to read/write .xtal files
          procedure, pass(self) :: readDataHDF_
          procedure, pass(self) :: saveDataHDF_
          procedure, pass(self) :: addXtalDataGroup_
          procedure, pass(self) :: getCrystalData_
          procedure, pass(self) :: dumpXtalInfo_
! miscellaneous routines
          procedure, pass(self) :: resetUnitCell
          procedure, pass(self) :: ShortestG_
          procedure, pass(self) :: GetAsymPosWyckoff_
          final :: Cell_destructor

          ! procedure, pass(self) :: Convert_kgs_to_Substrate
          ! procedure, pass(self) :: setEMsoftXtalSystem
          ! procedure, pass(self) :: getEMsoftXtalSystem

          generic, public :: resetCell => resetUnitCell

          generic, public :: calcMatrices => computeMatrices
          generic, public :: transSpace => TransSpaceDouble, TransSpaceSingle
          generic, public :: transCoor => transformCoordinates
          generic, public :: calcDot => calcDotSingle, calcDotDouble
          generic, public :: normVec => normVecSingle, normVecDouble
          generic, public :: calcLength => calcLengthSingle, calcLengthDouble
          generic, public :: calcAngle => calcAngleSingle, calcAngleDouble
          generic, public :: calcCross => calcCrossSingle, calcCrossDouble
          generic, public :: calcPositions => calcPositions_
          generic, public :: convertfromRtoH => convertfromRtoH_
          generic, public :: getLatParm => getLatticeParameterSingle, getLatticeParametersAll
          generic, public :: setLatParm => setLatticeParameterSingle, setLatticeParametersAll
          generic, public :: setLatParm => requestLatticeParameters
          generic, public :: setXtalSystem => setXtalSystem_
          generic, public :: getAsymPos => getAsymmetricPosition
          generic, public :: getAsymPosData => getAsymPosArray
          generic, public :: getdsm => getDirectStructureMatrix
          generic, public :: getrsm => getReciprocalStructureMatrix
          generic, public :: getdmt => getDirectMetricTensor
          generic, public :: getrmt => getReciprocalMetricTensor
          generic, public :: calcDensity => calcTheoreticalDensity
          generic, public :: getDensity => getDensity_
          generic, public :: GetAsymPosWyckoff => GetAsymPosWyckoff_
          generic, public :: getnumat => getnumat_
          generic, public :: setnumat => setnumat_
          generic, public :: getapos => getapos_
          generic, public :: setWyckoff => setWyckoff_
          generic, public :: setFileName => setFileName_
          generic, public :: getFileName => getFileName_
          generic, public :: setSource => setSource_
          generic, public :: getSource => getSource_
          generic, public :: getVolume => getVolume_
          generic, public :: setNatomtype => setNatomtype_
          generic, public :: getNatomtype => getNatomtype_
          generic, public :: getAtomtype => getatomtype_

          generic, public :: setAtomtype => setatomtype_
          generic, public :: readDataHDF => readDataHDF_
          generic, public :: saveDataHDF => saveDataHDF_
          generic, public :: addXtalDataGroup => addXtalDataGroup_
          generic, public :: setAtomPos => setAtomPos_, setAtomPosAll_
          generic, public :: getCrystalData => getCrystalData_
          generic, public :: dumpXtalInfo => dumpXtalInfo_
          generic, public :: ShortestG => ShortestG_

  end type Cell_T

! the following parameters and arrays need to be moved to somewhere else ...
      ! integer(kind=irg)                    :: DynNbeams, DynNbeamsLinked, nns, nnw, numscatfac

  ! the constructor routine for this class
  interface Cell_T
    module procedure Cell_constructor
  end interface Cell_T

contains

!--------------------------------------------------------------------------
type(Cell_T) function Cell_constructor( latparm ) result(cell)
!DEC$ ATTRIBUTES DLLEXPORT :: Cell_constructor
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! constructor for the Cell Class

IMPLICIT NONE

real(kind=dbl), INTENT(IN), OPTIONAL :: latparm(6)
 !! optional lattice parameter array; cell will be reset to zero if not present

if (present(latparm)) then
! set the lattice parameters
  cell%a = latparm(1)
  cell%b = latparm(2)
  cell%c = latparm(3)
  cell%alpha = latparm(4)
  cell%beta  = latparm(5)
  cell%gamma = latparm(6)

! initialize all the relevant geometry matrices
  call cell%calcMatrices()
else 
  call cell%resetUnitCell()
end if

end function Cell_constructor

!--------------------------------------------------------------------------
subroutine Cell_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: Cell_destructor
!! author: MDG
!! version: 1.0
!! date: 02/02/20
!!
!! destructor for the Cell_T Class

IMPLICIT NONE

type(Cell_T), INTENT(INOUT)  :: self

call reportDestructor('Cell_T')

! and deallocate any arrays
if (allocated(self%apos)) deallocate(self%apos)

end subroutine Cell_destructor

!--------------------------------------------------------------------------
recursive subroutine resetUnitCell(self)
!DEC$ ATTRIBUTES DLLEXPORT :: resetUnitCell
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! reset all unit cell variables to zero

IMPLICIT NONE

class(Cell_T), intent(inout)    :: self

! initialize cell to zero
 self%a = 0.0_dbl
 self%b = 0.0_dbl
 self%c = 0.0_dbl
 self%alpha = 0.0_dbl
 self%beta  = 0.0_dbl
 self%gamma = 0.0_dbl
 self%vol   = 0.0_dbl
 self%dmt = 0.0_dbl
 self%rmt = 0.0_dbl
 self%dsm = 0.0_dbl
 self%rsm = 0.0_dbl
 self%trigmat = 0.0_dbl
 self%xtal_system = 0
 self%ATOM_type = 0_irg
 self%ATOM_ntype = 0_irg
 self%ATOM_pos = 0.0_dbl
 self%xtalname = ''
 self%source = ''

! and deallocate any arrays
 if (allocated(self%apos)) deallocate(self%apos)

end subroutine resetUnitCell

!--------------------------------------------------------------------------
recursive function getLatticeParameterSingle( self, p ) result(lp)
!DEC$ ATTRIBUTES DLLEXPORT :: getLatticeParameterSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! return a single lattice parameter

use mod_global

IMPLICIT NONE

class(Cell_T),intent(inout)    :: self
character(*), intent(in)       :: p

real(kind=dbl)                 :: lp

select case(p)
  case('a')
    lp = self%a
  case('b')
    lp = self%b
  case('c')
    lp = self%c
  case('alpha')
    lp = self%alpha
  case('beta')
    lp = self%beta
  case('gamma')
    lp = self%gamma
  end select

end function getLatticeParameterSingle

!--------------------------------------------------------------------------
recursive subroutine setLatticeParameterSingle( self, p, lp )
!DEC$ ATTRIBUTES DLLEXPORT :: setLatticeParameterSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/28/20
  !!
  !! set a single lattice parameter

use mod_global

IMPLICIT NONE

class(Cell_T),intent(inout)    :: self
character(*), intent(in)       :: p
real(kind=dbl), intent(in)     :: lp

select case(p)
  case('a')
    self%a = lp
  case('b')
    self%b = lp
  case('c')
    self%c = lp
  case('alpha')
    self%alpha = lp
  case('beta')
    self%beta = lp
  case('gamma')
    self%gamma = lp
  end select

end subroutine setLatticeParameterSingle

!--------------------------------------------------------------------------
recursive function getLatticeParametersAll( self ) result(lp)
!DEC$ ATTRIBUTES DLLEXPORT :: getLatticeParametersAll
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! return all lattice parameters

use mod_global

IMPLICIT NONE

class(Cell_T),intent(inout)    :: self

real(kind=dbl)                 :: lp(6)

lp = (/ self%a, self%b, self%c, self%alpha, self%beta, self%gamma /)

end function getLatticeParametersAll

!--------------------------------------------------------------------------
recursive subroutine setLatticeParametersAll( self, lp)
!DEC$ ATTRIBUTES DLLEXPORT :: setLatticeParametersAll
  !! author: MDG
  !! version: 1.0
  !! date: 01/28/20
  !!
  !! set all lattice parameters

use mod_global

IMPLICIT NONE

class(Cell_T),intent(inout)    :: self
real(kind=dbl),intent(in)      :: lp(6)

self%a = lp(1)
self%b = lp(2)
self%c = lp(3)
self%alpha = lp(4)
self%beta = lp(5)
self%gamma = lp(6)

end subroutine setLatticeParametersAll

!--------------------------------------------------------------------------
recursive subroutine computeMatrices(self)
!DEC$ ATTRIBUTES DLLEXPORT :: computeMatrices
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! computes important crystallographic matrices

use mod_io
use mod_global

IMPLICIT NONE

class(Cell_T),intent(inout)   :: self

type(IO_T)                    :: Message

! auxiliary variables for geometric computation
real(kind=dbl)                :: det,ca,cb,cg,sa,sb,sg,tg,pirad
real(kind=dbl)                :: Mx(3,3), My(3,3), x, y

! auxiliary variables for the various tensors
 pirad = cPi/180.0_dbl
 ca = dcos(pirad*self%alpha)
 cb = dcos(pirad*self%beta)
 cg = dcos(pirad*self%gamma)
 sa = dsin(pirad*self%alpha)
 sb = dsin(pirad*self%beta)
 sg = dsin(pirad*self%gamma)
 tg = dtan(pirad*self%gamma)

 if (sg.eq.0.0_dbl) call Message%printError('mod_crystallography:computeMatrices',' Invalid gamma angle')

! compute the direct metric tensor [equation 1.5, page 6]
 self%dmt(1,1) = self%a**2
 self%dmt(2,2) = self%b**2
 self%dmt(3,3) = self%c**2
 self%dmt(1,2) = self%a*self%b*cg
 self%dmt(2,1) = self%dmt(1,2)
 self%dmt(1,3) = self%a*self%c*cb
 self%dmt(3,1) = self%dmt(1,3)
 self%dmt(2,3) = self%b*self%c*ca
 self%dmt(3,2) = self%dmt(2,3)
! cell volume via the determinant of dmt
 det = (self%a*self%b*self%c)**2*(1.D0-ca**2-cb**2-cg**2+2.D0*ca*cb*cg)
 self%vol = dsqrt(det)

 if (self%vol.lt.1D-6) call Message%printError('mod_crystallography:computeMatrices',' Unit cell volume zero or suspiciously small')

! compute the reciprocal metric tensor as the inverse of the direct
! metric tensor
 self%rmt(1,1) = (self%b*self%c*sa)**2
 self%rmt(2,2) = (self%a*self%c*sb)**2
 self%rmt(3,3) = (self%a*self%b*sg)**2
 self%rmt(1,2) = self%a*self%b*self%c**2*(ca*cb-cg)
 self%rmt(2,1) = self%rmt(1,2)
 self%rmt(1,3) = self%a*self%b**2*self%c*(cg*ca-cb)
 self%rmt(3,1) = self%rmt(1,3)
 self%rmt(2,3) = self%a**2*self%b*self%c*(cb*cg-ca)
 self%rmt(3,2) = self%rmt(2,3)
 self%rmt = self%rmt/det

! compute the direct structure matrix [equation 1.64, page 57]
 self%dsm(1,1) = self%a
 self%dsm(1,2) = self%b*cg
 self%dsm(1,3) = self%c*cb
 self%dsm(2,1) = 0.0_dbl
 self%dsm(2,2) = self%b*sg
 self%dsm(2,3) = -self%c*(cb*cg-ca)/sg
 self%dsm(3,1) = 0.0_dbl
 self%dsm(3,2) = 0.0_dbl
 self%dsm(3,3) = self%vol/(self%a*self%b*sg)

! compute the reciprocal structure matrix [equation 1.65, page 58]
 self%rsm(1,1) = 1.0_dbl/self%a
 self%rsm(1,2) = 0.0_dbl
 self%rsm(1,3) = 0.0_dbl
 self%rsm(2,1) = -1.0_dbl/(self%a*tg)
 self%rsm(2,2) = 1.0_dbl/(self%b*sg)
 self%rsm(2,3) = 0.0_dbl
 self%rsm(3,1) = self%b*self%c*(cg*ca-cb)/(self%vol*sg)
 self%rsm(3,2) = self%a*self%c*(cb*cg-ca)/(self%vol*sg)
 self%rsm(3,3) = (self%a*self%b*sg)/self%vol

! finally, if we have the trigonal/rhombohedral case, we need a second direct structure matrix
! that can transform the three-fold [111] axis to the z-axis of a cartesian reference frame,
! while keeping the rhombohedral [100] direction in the x-z plane.  This is used for the k-vector
! sampling routines for master pattern calculations.
!
! added by MDG on 08/30/15; computations in Mathematica notebook rhombohedral.nb in manuals folder, validated.
if (self%xtal_system.eq.5) then
! Mx matrix
  x = 0.5D0/dcos(pirad*0.5D0*self%alpha)
  y = dsqrt(1.D0-x*x)
  Mx = transpose(reshape( (/ 1.D0, 0.D0, 0.D0,  0.D0, x, -y,  0.D0, y, x /), (/ 3,3 /) ))
! My matrix
  x = 2.0D0*dsin(pirad*0.5D0*self%alpha)/dsqrt(3.D0)
  y = dsqrt(1.D0-x*x)
  My = transpose(reshape( (/ x, 0.D0, -y,  0.D0, 1.0D0, 0.D0,  y, 0.D0, x /), (/ 3,3 /) ))
  self%trigmat = matmul(matmul(My,Mx),self%dsm)
end if

end subroutine computeMatrices


!--------------------------------------------------------------------------
recursive subroutine dumpXtalInfo_(self, SG)
!DEC$ ATTRIBUTES DLLEXPORT :: dumpXtalInfo_
  !! author: MDG
  !! version: 1.0
  !! date: 01/12/20
  !!
  !! a brief summary of the crystal structure

use mod_io
use mod_symmetry
use mod_HallSG

IMPLICIT NONE

class(Cell_T), intent(inout)            :: self
type(SpaceGroup_T), intent(inout)       :: SG

type(IO_T)                              :: Message

integer(kind=irg)                       :: i, j, oi_int(3)
real(kind=dbl)                          :: oi_real(5)


 call Message%printMessage('', frm = "(A/)")
 call Message%printMessage('Crystal Structure Information', frm = "('-->',A,'<--')")
 oi_real(1) = self%getLatParm('a')
 call Message%WriteValue('  a [nm]             : ', oi_real, 1, "(F9.5)")
 oi_real(1) = self%getLatParm('b')
 call Message%WriteValue('  b [nm]             : ', oi_real, 1, "(F9.5)")
 oi_real(1) = self%getLatParm('c')
 call Message%WriteValue('  c [nm]             : ', oi_real, 1, "(F9.5)")
 oi_real(1) = self%getLatParm('alpha')
 call Message%WriteValue('  alpha [deg]        : ', oi_real, 1, "(F9.5)")
 oi_real(1) = self%getLatParm('beta')
 call Message%WriteValue('  beta  [deg]        : ', oi_real, 1, "(F9.5)")
 oi_real(1) = self%getLatParm('gamma')
 call Message%WriteValue('  gamma [deg]        : ', oi_real, 1, "(F9.5)")
 oi_real(1) = self%vol
 call Message%WriteValue('  Volume [nm^3]      : ', oi_real, 1, "(F12.8)")
 oi_int(1) = SG%getSpaceGroupNumber()
 call Message%WriteValue('  Space group #      : ', oi_int, 1, "(1x,I3)")
 if (SG%getuseHallSG().eqv..TRUE.) then 
   oi_int(1) = SG%getHallSpaceGroupNumber()
   call Message%WriteValue('  Hall Space group # : ', oi_int, 1, "(1x,I3)")
   call Message%WriteValue('  Hall space group   : ', trim(get_HallString(SG%getHallSpaceGroupNumber())) )
 else
   call Message%WriteValue('  Space group symbol : ', trim(SYM_SGname(SG%getSpaceGroupNumber())) )
   call Message%WriteValue('  Generator String   : ',  trim(SYM_GL(SG%getSpaceGroupNumber())) )
   if ((SG%getSpaceGroupSetting().eq.2).AND.(SG%getSpaceGroupXtalSystem().ne.5)) then
    call Message%printMessage('   Using second origin setting', frm = "(A)")
   endif
   if ((SG%getSpaceGroupSetting().eq.2).AND.(SG%getSpaceGroupXtalSystem().eq.5)) then
    call Message%printMessage('   Using rhombohedral parameters', frm = "(A)")
   endif
 end if 
 if (SG%getSpaceGroupCentro()) then
    call Message%printMessage('   Structure is centrosymmetric', frm = "(A)")
 else
   call Message%printMessage('   Structure is non-centrosymmetric', frm = "(A)")
 end if

! space group and point group information
 if (SG%getuseHallSG().eqv..TRUE.) then 
   oi_int(1) = SG%HallSG%get_NHallgenerators()
 else
   oi_int(1) = SG%getSpaceGroupGENnum()
 end if 
 call Message%WriteValue('  # generators       : ', oi_int, 1, "(1x,I3)")
 oi_int(1) = SG%getSpaceGroupMATnum()
 call Message%WriteValue('  # symmetry matrices: ', oi_int, 1, "(1x,I3)")
 oi_int(1) = SG%getSpaceGroupNUMpt()
 call Message%WriteValue('  # point sym. matr. : ', oi_int, 1, "(1x,I3)")

! generate atom positions and dump output
 call Message%printMessage('', frm = "(A/)")
 call self%calcPositions(SG, 'v')
 oi_int(1) = self%ATOM_ntype
 call Message%WriteValue('  Number of asymmetric atom positions ', oi_int, 1)
 do i=1,self%ATOM_ntype
  oi_int(1:3) = (/i, self%ATOM_type(i), self%numat(i)/)
  call Message%WriteValue('  General position / atomic number / multiplicity :', oi_int, 3,"(1x,I3,'/',I2,'/',I3)",advance="no")
  call Message%printMessage(' ('//ATOM_sym(self%ATOM_type(i))//')', frm = "(A)")
  call Message%printMessage('   Equivalent positions  (x y z  occ  DWF) ', frm = "(A)")
  do j=1,self%numat(i)
    oi_real(1:5) = (/self%apos(i, j,1:3),dble(self%ATOM_pos(i,4:5))/)
    call Message%WriteValue('         > ', oi_real, 5,"(2x,4(F9.5,','),F9.5)")
  end do
end do
call Message%printMessage('', frm = "(A/)")

call self%calcDensity()
oi_real(1:3) = self%getDensity()
call Message%WriteValue(' Density, avA, avZ = ',oi_real,3,"(2f10.5,',',f10.5)")

end subroutine dumpXtalInfo_

!--------------------------------------------------------------------------
recursive subroutine setWyckoff_(self, useWyckoff)
!DEC$ ATTRIBUTES DLLEXPORT :: setWyckoff_
  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20
  !!
  !! set the Wyckoff parameter

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
logical, INTENT(IN)             :: useWyckoff

self%Wyckoff = useWyckoff

end subroutine setWyckoff_

!--------------------------------------------------------------------------
recursive subroutine setFileName_(self, xtalname)
!DEC$ ATTRIBUTES DLLEXPORT :: setFileName_

  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20
  !!
  !! set the output filename

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
character(fnlen),INTENT(IN)     :: xtalname

self%xtalname = trim(xtalname)

end subroutine setFileName_

!--------------------------------------------------------------------------
recursive function getFileName_(self) result(xtalname)
!DEC$ ATTRIBUTES DLLEXPORT :: getFileName_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! get the .xtal filename

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
character(fnlen)                :: xtalname

xtalname = trim(self%xtalname)

end function getFileName_

!--------------------------------------------------------------------------
recursive function getDirectStructureMatrix(self) result(dsm)
!DEC$ ATTRIBUTES DLLEXPORT :: getDirectStructureMatrix
  !! author: MDG
  !! version: 1.0
  !! date: 01/26/20
  !!
  !! get the direct structure matrix

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
real(kind=dbl)                  :: dsm(3,3)

dsm = self%dsm

end function getDirectStructureMatrix

!--------------------------------------------------------------------------
recursive function getReciprocalStructureMatrix(self) result(rsm)
!DEC$ ATTRIBUTES DLLEXPORT :: getReciprocalStructureMatrix

  !! author: MDG
  !! version: 1.0
  !! date: 01/26/20
  !!
  !! get the reciprocal structure matrix

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
real(kind=dbl)                  :: rsm(3,3)

rsm = self%rsm

end function getReciprocalStructureMatrix

!--------------------------------------------------------------------------
recursive function getDirectMetricTensor(self) result(dmt)
!DEC$ ATTRIBUTES DLLEXPORT :: getDirectMetricTensor
  !! author: MDG
  !! version: 1.0
  !! date: 01/27/20
  !!
  !! get the direct metric tensor

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
real(kind=dbl)                  :: dmt(3,3)

dmt = self%dmt

end function getDirectMetricTensor

!--------------------------------------------------------------------------
recursive function getReciprocalMetricTensor(self) result(rmt)
!DEC$ ATTRIBUTES DLLEXPORT :: getReciprocalMetricTensor
  !! author: MDG
  !! version: 1.0
  !! date: 01/27/20
  !!
  !! get the reciprocal metric tensor

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
real(kind=dbl)                  :: rmt(3,3)

rmt = self%rmt

end function getReciprocalMetricTensor

!--------------------------------------------------------------------------
recursive subroutine setSource_(self, source)
!DEC$ ATTRIBUTES DLLEXPORT :: setSource_
  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20
  !!
  !! set the Source

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
character(fnlen),INTENT(IN)     :: source

self%source = trim(source)

end subroutine setSource_

!--------------------------------------------------------------------------
recursive function getSource_(self) result(source)
!DEC$ ATTRIBUTES DLLEXPORT :: getSource_
  !! author: MDG
  !! version: 1.0
  !! date: 01/16/20
  !!
  !! get the Source

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)   :: self
character(fnlen)               :: source

source = trim(self%source)

end function getSource_

!--------------------------------------------------------------------------
recursive function getVolume_(self) result(vol)
!DEC$ ATTRIBUTES DLLEXPORT :: getVolume_
  !! author: MDG
  !! version: 1.0
  !! date: 01/27/20
  !!
  !! get the cell volume

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)   :: self
real(kind=dbl)                 :: vol

vol = self%vol

end function getVolume_

!--------------------------------------------------------------------------
recursive subroutine setNatomtype_(self, n)
!DEC$ ATTRIBUTES DLLEXPORT :: setNatomtype_
  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20
  !!
  !! set the Source

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
integer(kind=irg),INTENT(IN)    :: n

self%ATOM_ntype = n

end subroutine setNatomtype_

!--------------------------------------------------------------------------
recursive function getNatomtype_(self) result(n)
!DEC$ ATTRIBUTES DLLEXPORT :: getNatomtype_
  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20
  !!
  !! set the Source

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
integer(kind=irg)               :: n

n = self%ATOM_ntype

end function getNatomtype_

!--------------------------------------------------------------------------
recursive subroutine setXtalSystem_(self, xs)
!DEC$ ATTRIBUTES DLLEXPORT :: setXtalSystem_
  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20
  !!
  !! set the crystal system (necessary to avoid circular dependency with space group class)

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
integer(kind=irg),INTENT(IN)    :: xs

self%xtal_system = xs

end subroutine setXtalSystem_

!--------------------------------------------------------------------------
recursive subroutine setAtomPos_(self, pos)
!DEC$ ATTRIBUTES DLLEXPORT :: setAtomPos_
  !! author: MDG
  !! version: 1.0
  !! date: 01/14/20
  !!
  !! set a single atom position

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
real(kind=sgl),INTENT(IN)       :: pos(3)

self%ATOM_pos = 0.0
self%ATOM_pos(1,1:3) = pos(1:3)

end subroutine setAtomPos_

!--------------------------------------------------------------------------
recursive subroutine setAtomPosAll_(self, pos)
!DEC$ ATTRIBUTES DLLEXPORT :: setAtomPosAll_
  !! author: MDG
  !! version: 1.0
  !! date: 01/14/20
  !!
  !! set a single atom position

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
real(kind=dbl),INTENT(IN)       :: pos(maxpasym,5)

self%ATOM_pos = pos

end subroutine setAtomPosAll_

!--------------------------------------------------------------------------
recursive function getatomtype_(self) result(atp)
!DEC$ ATTRIBUTES DLLEXPORT :: getatomtype_
  !! author: MDG
  !! version: 1.0
  !! date: 01/16/20
  !!
  !! returns the ATOM_type array to the calling program

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
integer(kind=irg), allocatable  :: atp(:)

integer(kind=irg)               :: sz(1)

sz = shape(self%ATOM_type)
allocate(atp(sz(1)))

atp = self%ATOM_type

end function getatomtype_

!--------------------------------------------------------------------------
recursive subroutine setatomtype_(self, atp)
!DEC$ ATTRIBUTES DLLEXPORT :: setatomtype_
  !! author: MDG
  !! version: 1.0
  !! date: 02/14/20
  !!
  !! sets the ATOM_type array

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
integer(kind=irg), INTENT(IN)   :: atp(1:self%ATOM_ntype)

self%ATOM_type(1:self%ATOM_ntype) = atp(1:self%ATOM_ntype)

end subroutine setatomtype_

!--------------------------------------------------------------------------
recursive function getnumat_(self) result(numat)
!DEC$ ATTRIBUTES DLLEXPORT :: getnumat_
  !! author: MDG
  !! version: 1.0
  !! date: 01/25/20
  !!
  !! returns the apos array to the calling program

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
integer(kind=irg), allocatable  :: numat(:)

integer(kind=irg)               :: sz(3)

allocate( numat( self%ATOM_ntype ) )

numat = self%numat(1:self%ATOM_ntype)

end function getnumat_

!--------------------------------------------------------------------------
recursive subroutine setnumat_(self, numat)
!DEC$ ATTRIBUTES DLLEXPORT :: setnumat_
  !! author: MDG
  !! version: 1.0
  !! date: 02/14/20
  !!
  !! set the number of atoms in the asymmetric unit

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
integer(kind=irg), INTENT(IN)   :: numat

self%ATOM_ntype = numat

end subroutine setnumat_

!--------------------------------------------------------------------------
recursive function getapos_(self) result(apos)
!DEC$ ATTRIBUTES DLLEXPORT :: getapos_
  !! author: MDG
  !! version: 1.0
  !! date: 01/25/20
  !!
  !! returns the apos array to the calling program

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
real(kind=dbl), allocatable     :: apos(:,:,:)

integer(kind=irg)               :: sz(3)

sz = shape(self%apos)
allocate(apos(sz(1),sz(2),sz(3)))

apos = self%apos

end function getapos_

!--------------------------------------------------------------------------
recursive function getAsymPosArray(self) result(asp)
!DEC$ ATTRIBUTES DLLEXPORT :: getAsymPosArray
  !! author: MDG
  !! version: 1.0
  !! date: 01/16/20
  !!
  !! returns the AsymPos array to the calling program

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)    :: self
real(kind=dbl), allocatable     :: asp(:,:)

integer(kind=irg)               :: sz(2)

sz = shape(self%ATOM_pos)
allocate(asp(sz(1),sz(2)))

asp = self%ATOM_pos

end function getAsymPosArray

!--------------------------------------------------------------------------
recursive subroutine getAsymmetricPosition(self)
!DEC$ ATTRIBUTES DLLEXPORT :: getAsymmetricPosition
  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20, updated 4/1/21 to new string parsing code in extractAtomPositionData
  !!
  !! ask the user for the atom type, coordinates, site occupation parameter
  !! and Debye-Waller parameter for each atom type.

use mod_io

IMPLICIT NONE

class(Cell_T),INTENT(INOUT)             :: self

type(IO_T)                              :: Message
logical                                 :: more                 ! logical to determine if more atoms need to be entered
character(1)                            :: ans, list(256)       ! used for IO
real(kind=sgl)                          :: pt(5), out_real(5)   ! used to read and write asymmetric position data
integer(kind=irg)                       :: i, j, io_int(1), std, sl   ! auxiliary variables
character(fnlen)                        :: instring

 more=.TRUE.
 self%ATOM_ntype = 0
 call Message%printMessage(' Enter atoms in the asymmetric unit ', frm = "(/A)")
 call self%displayPeriodicTable()

 do while (more)
  self%ATOM_ntype = self%ATOM_ntype + 1

! atomic number
  call Message%ReadValue(' ->  Atomic number : ', io_int, 1)
  self%ATOM_type(self%ATOM_ntype) = io_int(1)

! general atom coordinate
  list = (/ (' ',j=1,256) /)
  call Message%printMessage(' ->  Fractional coordinates, site occupation, and Debye-Waller Factor [nm^2] : ', &
                            frm = "(A,' ')",advance="no")
  read (stdin,"(A)") instring

! we can simply read the numbers from the string unless there is a / character which requires
! more detailed parsing using extractposition() ... so we look for potential / characters first 
  if (index(instring,'/').ne.0) then 
    call extractAtomPositionData(instring,pt)   
  else  
    read (instring,*) pt
  end if

! check for negative values and convert them to positive by adding 10.0 and using mod; 
! i.e., reduction to interval [0..1[
  do i=1,3
    if (pt(i).lt.0.0) pt(i) = mod(dble(pt(i))+10.D0,1.D0)
  end do

! default values 
  if (pt(4).eq.0.0) pt(4) = 1.0   ! this may not be a good assumption...but why define an atom with zero occupation factor ???
  if (pt(5).eq.0.0) pt(5) = 0.004 ! we always need a non-zero Debye-Waller factor ! 

! store in the appropriate component of the cell variable
  self%ATOM_pos(self%ATOM_ntype,1:5) = pt(1:5)

! and write the coordinate back to the terminal
  out_real = (/ (self%ATOM_pos(self%ATOM_ntype,j),j=1,5) /)
  call Message%WriteValue('    -> ', out_real, 5, frm = "(1x,4(F10.7,2x),F10.7)")

  call Message%ReadValue(' ->  Another atom ? (y/n) ', ans, frm = "(A1)")
  if ((ans.eq.'y').or.(ans.eq.'Y')) then
   more=.TRUE.
  else
   more=.FALSE.
  end if

 end do

end subroutine getAsymmetricPosition

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! routines for crystallographic computations
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine TransSpaceSingle(self, t, d, inspace, outspace)
!DEC$ ATTRIBUTES DLLEXPORT :: TransSpaceSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! convert vector components from inspace to outspace ('d', 'r', 'c') (single precision)

IMPLICIT NONE

class(Cell_T), intent(in)       :: self
real(kind=sgl),INTENT(IN)       :: t(3)
 !! input vector in inspace reference frame
real(kind=sgl),INTENT(OUT)      :: d(3)
 !! output vector in outspace reference frame
character(1),INTENT(IN)         :: inspace
 !! characters to label input space (d, r, or c)
character(1),INTENT(IN)         :: outspace
 !! characters to label output space (d, r, or c)

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

end subroutine TransSpaceSingle

!--------------------------------------------------------------------------
recursive subroutine TransSpaceDouble(self,t,d,inspace,outspace)
!DEC$ ATTRIBUTES DLLEXPORT :: TransSpaceDouble
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! convert vector components from inspace to outspace ('d', 'r', 'c') (double precision)

IMPLICIT NONE

class(Cell_T), intent(in)       :: self
real(kind=dbl),INTENT(IN)       :: t(3)
 !! input vector in inspace reference frame
real(kind=dbl),INTENT(OUT)      :: d(3)
 !! output vector in outspace reference frame
character(1),INTENT(IN)         :: inspace
 !! characters to label input space (d, r, or c)
character(1),INTENT(IN)         :: outspace
 !! characters to label output space (d, r, or c)

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

end subroutine TransSpaceDouble

!--------------------------------------------------------------------------
recursive subroutine transformCoordinates(self, t, d, talpha, space, direction)
!DEC$ ATTRIBUTES DLLEXPORT :: transformCoordinates
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! convert vector components from one reference frame
  !! to another; this is a general coordinate transformation using the
  !! old-to-new matrix alpha.  The details of this routine are summarized in
  !! Table 1.6, page 51, of the textbook. The direction of the
  !! transformation is 'on' (old-to-new) or 'no' (new-to-old).

IMPLICIT NONE

class(Cell_T), intent(in)       :: self
real(kind=dbl),INTENT(IN)       :: t(3)
 !! input vector w.r.t. input space reference frame
real(kind=dbl),INTENT(OUT)      :: d(3)
 !! transformed vector components
real(kind=dbl),INTENT(IN)       :: talpha(3,3)
 !! transformation matrix
real(kind=dbl)                  :: alinv(3,3)
 !! inverse of transformation matrix
character(1),INTENT(IN)         :: space
 !! space in which to perform transformation ('d', 'r', 'c')
character(2),INTENT(IN)         :: direction
 !! transformation direction (no=new-to-old, on=old-to-new)

! these matrices are typically unitary, so inverse is simply the transpose
 if (space.eq.'d') then
  if (direction.eq.'on') then
   d = matmul(t,transpose(talpha))
  else
   d = matmul(t,talpha)
  end if
 else
  if (direction.eq.'on') then
   d = matmul(talpha,t)
  else
   d = matmul(transpose(talpha),t)
  end if
 end if

end subroutine transformCoordinates

!--------------------------------------------------------------------------
recursive function CalcDotSingle(self, p, q, space) result(cdot)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcDotSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! computes the dot product between two vectors in
  !! real, reciprocal, or Cartesian space; implements
  !! equations 1.6 (page 7), and 1.16 (page 15).  (single precision)

IMPLICIT NONE

class(Cell_T),intent(in)        :: self
real(kind=sgl),INTENT(IN)       :: p(3)
 !! first input vector in space reference frame
real(kind=sgl),INTENT(IN)       :: q(3)
 !! second input vector
character(1),INTENT(IN)         :: space
 !! space in which to compute product ('d', 'r', or 'c')
real(kind=sgl)                  :: cdot
 !! dot product p.q

 cdot = 0.0_sgl
 if (space.eq.'d') cdot = dot_product(p,matmul(self%dmt,q))
 if (space.eq.'r') cdot = dot_product(p,matmul(self%rmt,q))
 if (space.eq.'c') cdot = dot_product(p,q)

end function CalcDotSingle

!--------------------------------------------------------------------------
recursive function CalcDotDouble(self, p, q, space) result(cdot)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcDotDouble
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! computes the dot product between two vectors in
  !! real, reciprocal, or Cartesian space; implements
  !! equations 1.6 (page 7), and 1.16 (page 15).  (double precision)

IMPLICIT NONE

class(Cell_T),intent(in)        :: self
real(kind=dbl),INTENT(IN)       :: p(3)
 !! first input vector in space reference frame
real(kind=dbl),INTENT(IN)       :: q(3)
 !! second input vector
character(1),INTENT(IN)         :: space
 !! space in which to compute product ('d', 'r', or 'c')
real(kind=dbl)                  :: cdot
 !! dot product p.q

 cdot = 0.0_dbl
 if (space.eq.'d') cdot = dot_product(p,matmul(self%dmt,q))
 if (space.eq.'r') cdot = dot_product(p,matmul(self%rmt,q))
 if (space.eq.'c') cdot = dot_product(p,q)

end function CalcDotDouble

!--------------------------------------------------------------------------
recursive subroutine NormVecSingle(self, p, space)
!DEC$ ATTRIBUTES DLLEXPORT :: NormVecSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! single precision vector normalization
  !!
  !! vector normalization in arbitrary space.  Note that it is not
  !! necessarily so that the Cartesian length of a vector [sqrt(x^2+y^2+z^2)]
  !! becomes unity after this normalization.  That is only the case for
  !! the cartesian metric; for all other metrics, the length of a normalized
  !! vector is in general not equal to 1.0 !

IMPLICIT NONE

class(Cell_T),intent(in)                :: self
real(kind=sgl),INTENT(INOUT)            :: p(3)
 !! input/output vector components
!f2py intent(in,out) ::  p
character(1),INTENT(IN)                 :: space
 !! space character ('d', 'r', or 'c')
real(kind=sgl)                          :: x
 !! auxiliary variable

 x = self%CalcLength(p, space)
 if (x.ne.0.0) then
   p = p/x
 else
   p = (/0.0,0.0,0.0/)
 end if

end subroutine NormVecSingle

!--------------------------------------------------------------------------
recursive subroutine NormVecDouble(self, p, space)
!DEC$ ATTRIBUTES DLLEXPORT :: NormVecDouble
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! double precision vector normalization
  !!
  !! vector normalization in arbitrary space.  Note that it is not
  !! necessarily so that the Cartesian length of a vector [sqrt(x^2+y^2+z^2)]
  !! becomes unity after this normalization.  That is only the case for
  !! the cartesian metric; for all other metrics, the length of a normalized
  !! vector is in general not equal to 1.0 !

IMPLICIT NONE

class(Cell_T),intent(in)                :: self
real(kind=dbl),INTENT(INOUT)            :: p(3)
 !! input/output vector components
!f2py intent(in,out) ::  p
character(1),INTENT(IN)                 :: space
 !! space character ('d', 'r', or 'c')
real(kind=dbl)                          :: x
 !! auxiliary variable

 x = self%CalcLength(p,space)
 if (x.ne.0.D0) then
   p = p/x
 else
   p = (/0.D0,0.D0,0.D0/)
 end if

end subroutine NormVecDouble

!--------------------------------------------------------------------------
recursive function CalcLengthSingle(self, p, space) result(x)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcLengthSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! single precision vector length

IMPLICIT NONE

class(Cell_T),intent(in)                :: self
real(kind=sgl),INTENT(IN)               :: p(3)
 !! input/output vector components
character(1),INTENT(IN)                 :: space
 !! space character ('d', 'r', or 'c')
real(kind=sgl)                          :: x
 !! auxiliary variable

 x = sqrt(self%CalcDot(p, p, space))

end function CalcLengthSingle

!--------------------------------------------------------------------------
recursive function CalcLengthDouble(self, p, space) result(x)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcLengthDouble
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! double precision vector length

IMPLICIT NONE

class(Cell_T),intent(in)                :: self
real(kind=dbl),INTENT(IN)               :: p(3)
 !! input/output vector components
character(1),INTENT(IN)                 :: space
 !! space character ('d', 'r', or 'c')
real(kind=dbl)                          :: x
 !! auxiliary variable

 x = dsqrt(self%CalcDot(p, p, space))

end function CalcLengthDouble

!--------------------------------------------------------------------------
recursive function CalcAngleSingle(self, p, q, space) result(a)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcAngleSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! single precision angle in arbitrary space

use mod_io
use mod_global

IMPLICIT NONE

class(Cell_T),intent(in)                :: self
real(kind=sgl),INTENT(IN)               :: p(3)
  !! first vector components
real(kind=sgl),INTENT(IN)               :: q(3)
  !! second vector components
character(1),INTENT(IN)                 :: space
  !! space of the computation ('d', 'r', 'c')
real(kind=sgl)                          :: a
  !! angle (radians)

real(kind=sgl)                          :: x, y, z, t
type(IO_T)                              :: Message

 x = self%CalcDot(p,q,space)
 y = self%CalcLength(p,space)
 z = self%CalcLength(q,space)

 if ((y.eq.0.0_sgl).or.(z.eq.0.0_sgl)) then
  call Message%printError('CalcAngleSingle',' vector of zero length specified')
 end if

 t = x/(y*z)
 if (t.ge.1.0_sgl) then
  a = 0.0_sgl
 else
  if (t.le.-1.0_sgl) then
   a = sngl(cPi)
  else
   a = acos(t)
  end if
 end if

end function CalcAngleSingle

!--------------------------------------------------------------------------
recursive function CalcAngleDouble(self,p,q,space) result(a)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcAngleDouble
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! double precision angle in arbitrary space

use mod_io
use mod_global

IMPLICIT NONE

class(Cell_T),intent(in)                :: self
real(kind=dbl),INTENT(IN)               :: p(3)
 !! first vector components
real(kind=dbl),INTENT(IN)               :: q(3)
 !! second vector components
character(1),INTENT(IN)                 :: space
 !! space of the computation ('d', 'r', 'c')
real(kind=dbl)                          :: a
 !! angle in radians

type(IO_T)                              :: Message
real(kind=dbl)                          :: x, y, z, t

 x = self%CalcDot(p,q,space)
 y = self%CalcLength(p,space)
 z = self%CalcLength(q,space)

 if ((y.eq.0.0_dbl).or.(z.eq.0.0_dbl)) then
  call Message%printError('CalcAngleDouble',' vector of zero length specified')
 end if

 t = x/(y*z)
 if (t.ge.1.0_dbl) then
  a = 0.0_dbl
 else
  if (t.le.-1.0_dbl) then
   a = cPi
  else
   a = dacos(t)
  end if
 end if

end function CalcAngleDouble

!--------------------------------------------------------------------------
recursive subroutine CalcCrossSingle(self,p,q,r,inspace,outspace,iv)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcCrossSingle
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! vector cross product in arbitrary space (single precision)
  !!
  !! computes the cross product between two vectors and
  !! expresses it in either real space or reciprocal space.
  !! The output can also be expressed in the standard
  !! Cartesian reference frame.  The switch iv indicates
  !! whether the result should be scaled by the unit cell
  !! volume. More information in section 1.3.5, page 18.


IMPLICIT NONE

class(Cell_T),intent(in)                :: self
real(kind=sgl),INTENT(IN)               :: p(3)
 !! first input vector (order is important here !)
real(kind=sgl),INTENT(IN)               :: q(3)
 !! second input vector
real(kind=sgl),INTENT(OUT)              :: r(3)
 !! output vector
character(1),INTENT(IN)                 :: inspace
 !! inspace character ('d','r','c')
character(1),INTENT(IN)                 :: outspace
 !! outspace character
integer(kind=irg),INTENT(IN)            :: iv
 !! volume division switch

real(kind=sgl)                          :: x(3), vl

! divide by volume?
 if (iv.eq.1) then
  vl = sngl(self%vol)
 else
  vl = 1.0_sgl
 endif

! in direct space
 if (inspace.eq.'d') then               ! so the output is in reciprocal space !
  r(1) = vl*(p(2)*q(3)-p(3)*q(2))
  r(2) = vl*(p(3)*q(1)-p(1)*q(3))
  r(3) = vl*(p(1)*q(2)-p(2)*q(1))
  if (outspace.eq.'d') then             ! output in direct space
   x = matmul(r,self%rmt)
   r = x
  end if
  if (outspace.eq.'c') then             ! output in cartesian frame
   x = matmul(self%rsm,r)
   r = x
  end if
 end if

! in reciprocal space
 if (inspace.eq.'r') then               ! so the output is in direct space !
  r(1) = (p(2)*q(3)-p(3)*q(2))/vl
  r(2) = (p(3)*q(1)-p(1)*q(3))/vl
  r(3) = (p(1)*q(2)-p(2)*q(1))/vl
  if (outspace.eq.'r') then             ! output in reciprocal space
   x = matmul(r,self%dmt)
   r = x
  end if
  if (outspace.eq.'c') then             ! output in cartesian frame
   x = matmul(self%dsm,r)
   r = x
  end if
 end if

! in  cartesian
 if (inspace.eq.'c') then               ! so no conversion needed.
  r(1) = p(2)*q(3)-p(3)*q(2)
  r(2) = p(3)*q(1)-p(1)*q(3)
  r(3) = p(1)*q(2)-p(2)*q(1)
 end if

end subroutine CalcCrossSingle

!--------------------------------------------------------------------------
recursive subroutine CalcCrossDouble(self,p,q,r,inspace,outspace,iv)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcCrossDouble
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! vector cross product in arbitrary space (single precision)
  !!
  !! computes the cross product between two vectors and
  !! expresses it in either real space or reciprocal space.
  !! The output can also be expressed in the standard
  !! Cartesian reference frame.  The switch iv indicates
  !! whether the result should be scaled by the unit cell
  !! volume. More information in section 1.3.5, page 18.

IMPLICIT NONE

class(Cell_T),intent(in)                :: self
real(kind=dbl),INTENT(IN)               :: p(3)
 !! first input vector (order is important here !)
real(kind=dbl),INTENT(IN)               :: q(3)
 !! second input vector
real(kind=dbl),INTENT(OUT)              :: r(3)
 !! output vector
character(1),INTENT(IN)                 :: inspace
 !! inspace character ('d','r','c')
character(1),INTENT(IN)                 :: outspace
 !! outspace character
integer(kind=irg),INTENT(IN)            :: iv
 !! volume division switch

real(kind=dbl)                          :: x(3), vl

 if (iv.eq.1) then
  vl = self%vol
 else
  vl = 1.0_dbl
 endif

! in direct space
 if (inspace.eq.'d') then               ! so the output is in reciprocal space !
  r(1) = vl*(p(2)*q(3)-p(3)*q(2))
  r(2) = vl*(p(3)*q(1)-p(1)*q(3))
  r(3) = vl*(p(1)*q(2)-p(2)*q(1))
  if (outspace.eq.'d') then             ! output in direct space
   x = matmul(r,self%rmt)
   r = x
  end if
  if (outspace.eq.'c') then             ! output in cartesian frame
   x = matmul(self%rsm,r)
   r = x
  end if
 end if

! in reciprocal space
 if (inspace.eq.'r') then               ! so the output is in direct space !
  r(1) = (p(2)*q(3)-p(3)*q(2))/vl
  r(2) = (p(3)*q(1)-p(1)*q(3))/vl
  r(3) = (p(1)*q(2)-p(2)*q(1))/vl
  if (outspace.eq.'r') then             ! output in reciprocal space
   x = matmul(r,self%dmt)
   r = x
  end if
  if (outspace.eq.'c') then             ! output in cartesian frame
   x = matmul(self%dsm,r)
   r = x
  end if
 end if

! in  cartesian
 if (inspace.eq.'c') then               ! so no conversion needed.
  r(1) = p(2)*q(3)-p(3)*q(2)
  r(2) = p(3)*q(1)-p(1)*q(3)
  r(3) = p(1)*q(2)-p(2)*q(1)
 end if

end subroutine CalcCrossDouble

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! routines for crystallographic IO
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine requestLatticeParameters(self, SG)
!DEC$ ATTRIBUTES DLLEXPORT :: requestLatticeParameters
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! Input of crystal system followed by the appropriate set of lattice
  !! parameters; all are stored in the cell class.
  !!
  !! For the trigonal and hexagonal crystal systems, the following switch settings are to be used:
  !!              xtal_system   hexset    SYM_trigonal   SYM_second
  !! hexagonal          4         T            F              F
  !! trig/hex           5         T            T              F
  !! trig/rhomb         5         F            T              T

use mod_io
use mod_symmetry
use mod_HallSG

IMPLICIT NONE

class(Cell_T),intent(inout)             :: self
type(SpaceGroup_T), intent(inout)       :: SG

type(IO_T)                              :: Message
integer(kind=irg)                       :: io_int(1), TRIG(7)
real(kind=dbl)                          :: io_real(1)
character(8)                            :: SGshort 
character(16)                           :: HS

  ! make sure the symmetry operations will be reduced to the
  ! fundamental unit cell
   call SG%setSpaceGroupreduce(.TRUE.)
   call SG%setSpaceGrouphexset(.FALSE.)

  ! deal with the rhombohedral vs. hexagonal setting in the trigonal crystal system
  ! (the rhombohedral axes are considered as the second setting)
   call SG%setSpaceGrouptrigonal(.FALSE.)
   call SG%setSpaceGroupsecond(.FALSE.)

   if (self%xtal_system.eq.5) then
    call SG%setSpaceGrouptrigonal(.TRUE.)
 
  ! this is a bit different if we have the Hall symbols vs. the regular space groups
    if (SG%getuseHallSG().eqv..TRUE.) then 
      TRIG = (/ 146,148,155,160,161,166,167 /)
      if (minval(abs(TRIG-SG%getSpaceGroupNumber())).eq.0) then ! we have one of these cases
       HS = get_HallString( SG%getHallSpaceGroupNumber(), SGshort )
       if (scan(SGshort,'H').ne.0) then  ! this is the hexagonal setting
          self%xtal_system=4   ! this is set to 4 so that we ask for the correct lattice parameters below
       else ! and this is the rhombohedral setting
          call SG%setSpaceGroupsecond(.TRUE.)
       end if         
     end if         
    else 
      call Message%printMessage('Enter 1 for rhombohedral setting ,', frm = "(/A)")
      call Message%ReadValue('0 for hexagonal setting : ', io_int, 1)
      if (io_int(1).eq.0) then
       self%xtal_system=4   ! this is set to 4 so that we ask for the correct lattice parameters below
      else
       call SG%setSpaceGroupsecond(.TRUE.)
      end if
    end if

   end if 

  ! get the lattice parameters
   call Message%printMessage('Enter lattice parameters', frm = "(//A)")

  ! put default values based on cubic symmetry, then change them later
   call Message%ReadValue('    a [nm] = ', io_real, 1)
   self%a = io_real(1)
   self%b = self%a
   self%c = self%a
   self%alpha = 90.0_dbl
   self%beta = 90.0_dbl
   self%gamma = 90.0_dbl

  ! now get the proper lattice parameters
   select case (self%xtal_system)
    case (1)
  ! tetragonal
    case (2)
     call Message%ReadValue('    c [nm] = ', io_real, 1)
     self%c = io_real(1)
  ! orthorhombic
    case (3)
     call Message%ReadValue('    b [nm] = ', io_real, 1)
     self%b = io_real(1)
     call Message%ReadValue('    c [nm] = ', io_real, 1)
     self%c = io_real(1)
  ! hexagonal
    case (4)
     call Message%ReadValue('    c [nm] = ', io_real, 1)
     self%c = io_real(1)
     self%gamma=120.0_dbl
  ! rhombohedral
    case (5)
     call Message%ReadValue('    alpha [deg] = ', io_real, 1)
     self%alpha = io_real(1)
     self%beta = self%alpha
     self%gamma = self%alpha
  ! monoclinic
    case (6)
     call Message%ReadValue('    b [nm] = ', io_real, 1)
     self%b = io_real(1)
     call Message%ReadValue('    c [nm] = ', io_real, 1)
     self%c = io_real(1)
     if (SG%getuseHallSG().eqv..TRUE.) then
! if we are using the Hall space group symbols, then for the monoclinic case we need 
! to distinguish between the three possible choices for the unique axis
       SGshort = SG%HallSG%get_HallSGlabel( SG%getHallSpaceGroupNumber() )
       if (scan(SGshort,'b').ne.0) then 
         call Message%ReadValue('    beta  [deg] = ', io_real, 1)
         self%beta = io_real(1)
       else if (scan(SGshort,'c').ne.0) then 
         call Message%ReadValue('    gamma  [deg] = ', io_real, 1)
         self%gamma = io_real(1)
       else ! this must be the "a" unique axis
         call Message%ReadValue('    alpha  [deg] = ', io_real, 1)
         self%alpha = io_real(1)
       end if
     else  ! standard space group setting so we ask for the beta angle
       call Message%ReadValue('    beta  [deg] = ', io_real, 1)
       self%beta = io_real(1)
     end if
  ! triclinic
    case (7)
     call Message%ReadValue('    b [nm] = ', io_real, 1)
     self%b = io_real(1)
     call Message%ReadValue('    c [nm] = ', io_real, 1)
     self%c = io_real(1)
     call Message%ReadValue('    alpha [deg] = ', io_real, 1)
     self%alpha = io_real(1)
     call Message%ReadValue('    beta  [deg] = ', io_real, 1)
     self%beta = io_real(1)
     call Message%ReadValue('    gamma [deg] = ', io_real, 1)
     self%gamma = io_real(1)
   end select

  ! if trigonal symmetry was selected in the first setting,
  ! then the xtal_system must be reset to 5
   if (SG%getSpaceGrouptrigonal()) then
     self%xtal_system=5
   end if

  ! if hexagonal setting is used, then Miller-Bravais indices must be enabled
   if ((self%xtal_system.eq.4).OR.((self%xtal_system.eq.5).AND.(.not.SG%getSpaceGroupsecond() ) )) then
     call SG%setSpaceGrouphexset(.TRUE.)
   else
     call SG%setSpaceGrouphexset(.FALSE.)
   end if
!end if

end subroutine requestLatticeParameters

! !--------------------------------------------------------------------------
! recursive subroutine GetAsymmetricPosition(self)
!   !! author: MDG
!   !! version: 1.0
!   !! date: 01/07/20
!   !!
!   !! ask the user for the atom type, coordinates, site occupation parameter
!   !! and Debye-Waller parameter for each atom type.

! use mod_io

! IMPLICIT NONE

! class(Cell_T),intent(inout)             :: self
! !f2py intent(in,out) ::  cell

! type(IO_T)                              :: Message
! logical                                 :: more
! character(1)                            :: ans, list(256)
! real(kind=sgl)                          :: pt(5), out_real(5)
! integer(kind=irg)                       :: i, j, io_int(1), std, sl
! character(fnlen)                        :: instring

!  more=.TRUE.
!  self%ATOM_ntype = 0
!  call Message%printMessage(' Enter atoms in asymmetric unit ', frm = "(/A)")
!  call self%DisplayElements()

!  do while (more)
!   self%ATOM_ntype = cell%ATOM_ntype + 1

! ! atomic number
!   call Message%ReadValue(' ->  Atomic number : ', io_int, 1)
!   self%ATOM_type(cell%ATOM_ntype) = io_int(1)

! ! general atom coordinate
!   list = (/ (' ',j=1,256) /)
!   call Message%printMessage(' ->  Fractional coordinates, site occupation, and Debye-Waller Factor [nm^2] : ', &
!                             frm = "(A,' ')",advance="no")
!   call Message%ReadValue('', instring)
!   sl = len(trim(instring))
!   j = 0
!   do i=1,sl
!     if (instring(i:i).ne.' ') then
!       j = j+1
!       list(j) = instring(i:i)
!     end if
!   end do

! ! interpret this string and extract coordinates and such ...
!   call extractAtomPositionData(list,pt)

! ! store in the appropriate component of the cell variable
!   self%ATOM_pos(self%ATOM_ntype,1:5) = pt(1:5)

! ! and write the coordinate back to the terminal
!   out_real = (/ (self%ATOM_pos(self%ATOM_ntype,j),j=1,5) /)
!   call Message%WriteValue('    -> ', out_real, 5, frm = "(1x,4(F10.7,2x),F10.7)")

!   call Message%ReadValue(' ->  Another atom ? (y/n) ', ans, frm = "(A1)")
!   if ((ans.eq.'y').or.(ans.eq.'Y')) then
!    more=.TRUE.
!   else
!    more=.FALSE.
!   end if

!  end do

! end subroutine GetAsymmetricPosition

!--------------------------------------------------------------------------
recursive subroutine GetAsymPosWyckoff_(self, SG, EMsoft)
!DEC$ ATTRIBUTES DLLEXPORT :: GetAsymPosWyckoff_
  !! author: MDG
  !! version: 1.0
  !! date: 01/10/20
  !!
  !! read the atom coordinates from standard input in Wyckoff format (uses mod_symmetry routines)

use mod_io
use mod_symmetry
use mod_EMsoft

IMPLICIT NONE

class(cell_T),INTENT(INOUT)             :: self
type(SpaceGroup_T),INTENT(INOUT)        :: SG
type(EMsoft_T),INTENT(INOUT)            :: EMsoft

type(IO_T)                              :: Message
logical                                 :: more, found              ! logical to determine if more atoms need to be entered
character(1)                            :: ans                      ! used for IO
real(kind=sgl)                          :: pt(3), out_real(5)       ! used to read and write asymmetric position data
integer(kind=irg)                       :: i, j, ipos, std, numsp   ! auxiliary variables
real(kind=sgl)                          :: io_real(3)
character(fnlen)                        :: wpstring, WPfile
character(3)                            :: Wyckoffpos
character(6)                            :: Wyckoffstring
character(6)                            :: WyckoffList(27)
character(6)                            :: list

 more=.TRUE.
 self%ATOM_ntype = 0
 call Message%printMessage(' Enter atoms in asymmetric unit using Wyckoff positions', frm = "(/A)")
 call self%displayPeriodicTable()
 WPfile = trim(EMsoft%getConfigParameter('WyckoffPositionsfilename'))
 call SG%printWyckoffPositions(wpstring, WPfile, WyckoffList)

! number of special positions encoded in string
 numsp = (len(trim(wpstring))-1)/4+1

 call Message%printMessage('WPstring : '//trim(wpstring))

 do while (more)
  self%ATOM_ntype = self%ATOM_ntype + 1

! atomic number
  call Message%ReadValue(' ->  Atomic number, site occupation, Debye-Waller factor : ', io_real, 3)
  self%ATOM_type(self%ATOM_ntype) = int(io_real(1))
  self%ATOM_pos(self%ATOM_ntype,4:5) = io_real(2:3)

! ask for the Wyckoff position and make sure it actually exists in the list
  found = .FALSE.
  do while (found.eqv..FALSE.)
! ask for the Wyckoff position
    call Message%printMessage(' ->  Wyckoff position : ', frm = "(A,' ')",advance="no")
    list = '      '
    call Message%ReadValue(' ', list, frm="(A)")

! find the corresponding encoded triplet
    do i=1,6
      Wyckoffstring(i:i) = list(i:i)
    end do
    do i=1,numsp
      if (trim(Wyckoffstring).eq.trim(WyckoffList(i))) then
        found = .TRUE.
        if (i.eq.numsp) then
          Wyckoffpos = 'xyz'
        else
          do j=1,3
            ipos = (i-1)*4+1+j
            Wyckoffpos(j:j) = wpstring(ipos:ipos)
          end do
        end if
      end if
    end do
    if (found.eqv..FALSE.) then
      call Message%printMessage(' incorrect Wyckoff position; please try again ', frm = "(A,' ')",advance="no")
!   else
!     write (*,*) 'Found Wyckoff position '//Wyckoffpos
    end if
  end do

! interpret this encoded string and extract coordinates and such ...
  call SG%extractWyckoffposition(Wyckoffpos, pt)

! store in the appropriate component of the cell variable
  self%ATOM_pos(self%ATOM_ntype,1:3) = pt(1:3)

! and write the coordinate back to the terminal
  out_real = (/ (self%ATOM_pos(self%ATOM_ntype,j),j=1,5) /)
  call Message%WriteValue('    -> ', out_real, 5, frm = "(1x,4(F10.7,2x),F10.7)")

  call Message%ReadValue(' ->  Another atom ? (y/n) ', ans, frm = "(A1)")
  if ((ans.eq.'y').or.(ans.eq.'Y')) then
   more=.TRUE.
  else
   more=.FALSE.
  end if

 end do

end subroutine GetAsymPosWyckoff_

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

class(Cell_T),intent(inout)             :: self
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
recursive subroutine findinline(list,sl,ch,pos,j,rep)
!DEC$ ATTRIBUTES DLLEXPORT :: findinline
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20

IMPLICIT NONE

integer(kind=irg), INTENT(IN)           :: sl
character(1),INTENT(INOUT)              :: list(sl)  
character(1), INTENT(IN)                :: ch 
integer(kind=irg), INTENT(INOUT)        :: pos(sl)
integer(kind=irg),INTENT(INOUT)         :: j
character(1), INTENT(IN), OPTIONAL      :: rep

integer(kind=irg)                       :: i

pos = -1
j=0
do i=1,sl
  if (list(i)(1:1).eq.ch) then 
    j = j+1 
    pos(j) = i 
    if (present(rep)) list(i)(1:1) = rep
  end if 
end do 

end subroutine findinline

!--------------------------------------------------------------------------
recursive subroutine extractAtomPositionData(instring,pt)
!DEC$ ATTRIBUTES DLLEXPORT :: extractAtomPositionData
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20; rewritten 04/01/21 based on EMsoft 5.0 version changes 
  !!
  !! Extract the coordinates, site occupation, and DW factor from the
  !! input string; note that coordinates can be entered in decimal or in fractional
  !! notation, hence the somewhat convoluted way of interpreting this string...

IMPLICIT NONE

character(fnlen),INTENT(IN)             :: instring
 !! input string
real(kind=sgl),INTENT(OUT)              :: pt(5)
 !! output real array

integer(kind=irg),allocatable           :: slash(:)
integer(kind=irg)                       :: scnt, ccnt, spcnt, i, j, sl, ps, vs             !< auxiliary variables
character(1),allocatable                :: list(:)       
integer(kind=irg),parameter             :: nmb(48:57)=(/0,1,2,3,4,5,6,7,8,9/)   !< list of numbers
real(kind=dbl)                          :: nominator,denominator,x              !< used for fraction interpretation
character(fnlen)                        :: tmpstr, newstr 
real(kind=dbl),allocatable              :: vals(:)
character(1)                            :: cors   ! comma separated or space separated ? 

! make a copy of the input string and convert it to an array of characters
 tmpstr = trim(instring)
 sl = len(tmpstr)
 j = 0
 allocate(list(sl))
 do i=1,sl
   list(i) = tmpstr(i:i)
 end do
 allocate(slash(sl))

! search for all / locations; findinline changes those to spaces
 call findinline(list, sl, '/', slash, scnt, rep=' ')
! 4 commas/spaces plus 1 plus number of slashes is the number of reals that need to extracted
 vs = 5 + scnt
 allocate(vals(vs))

! next we generate a new string of more than 5 numbers 
 do i=1,sl 
   newstr(i:i) = list(i)
 end do

! read all numbers from this newstr into an array of floats
 read(newstr,*) vals

! for those numbers that are part of fractions, they must be larger than or equal to 1.0, so we look for them 
 j=1  ! index into vals array
 i=1  ! index into pt array
 do while (i.lt.6) 
   if ((abs(vals(j)).ge.1.0).and.(j.lt.vs-2)) then 
     pt(i) = vals(j) / vals(j+1)
     j = j+2
   else 
     pt(i) = vals(j)
     j = j+1
   end if 
   i = i+1 
 end do

! and clean up
 deallocate(slash, vals)

end subroutine extractAtomPositionData

!--------------------------------------------------------------------------
recursive function getDensity_(self) result(dza)
!DEC$ ATTRIBUTES DLLEXPORT :: getDensity_

IMPLICIT NONE

class(Cell_T), INTENT(INOUT)  :: self
real(kind=dbl)                :: dza(3)

dza = (/ self%density, self%avA, self%avZ /)

end function getDensity_

!--------------------------------------------------------------------------
recursive subroutine calcTheoreticalDensity(self, Z2percent)
!DEC$ ATTRIBUTES DLLEXPORT :: calcTheoreticalDensity

  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! compute the theoretical density as well as average Z and A

use mod_global

IMPLICIT NONE

class(Cell_T),intent(inout)             :: self

real(kind=sgl),OPTIONAL,INTENT(OUT),allocatable  :: Z2percent(:)
 !! optional array with percentages of Z^2 values for each general atom type

real(kind=sgl),allocatable              :: Z2list(:)
real(kind=dbl)                          :: AW, Z
integer(kind=irg)                       :: i

! compute the total atomic weight for the unit cell (g/mol)
! also compute the total atomic number
AW = 0.D0
Z = 0.D0
do i = 1, self % ATOM_ntype
  AW = AW + self%numat(i) * ATOM_weights(self % ATOM_type(i)) * self % ATOM_pos(i,4)
  Z = Z + self%numat(i) * float(self % ATOM_type(i))
end do
self%avA = AW/sum( self%numat(1:self % ATOM_ntype) * self%ATOM_pos(1:self%ATOM_ntype,4) )
self%avZ = Z/dble(sum( self%numat(1:self % ATOM_ntype) ))

! and compute the density in gram/centimeter^3
self%density = AW / (self % vol * 1.D-21 * cAvogadro)

! do we need to fill the Z2percent array?
! this estimates the percentage contribution of each atom type to the
! Rutherford scattering process by simply taking Z^2 for each atom in the unit cell
if (present(Z2percent)) then
  allocate(Z2percent(self%ATOM_ntype),Z2list(self%ATOM_ntype))
  do i=1, self%ATOM_ntype
    Z2list(i) = self%numat(i) * self%ATOM_pos(i,4) * self%ATOM_type(i)**2
  end do
  Z = sum(Z2list)
  Z2percent = 100.0 * Z2list / Z
end if

end subroutine calcTheoreticalDensity


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! HDF support routines for writing and reading .xtal files
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine saveDataHDF_(self, SG, EMsoft)
!DEC$ ATTRIBUTES DLLEXPORT :: saveDataHDF_
  !! author: MDG
  !! version: 1.0
  !! date: 01/14/20
  !!
  !! save crystal structure data to an HDF file

use mod_io
use mod_global
use mod_EMsoft
use HDF5
use mod_HDFsupport
use mod_timing
use mod_symmetry
use stringconstants
use mod_HallSG

IMPLICIT NONE

class(Cell_T),INTENT(INOUT)             :: self
type(SpaceGroup_T),INTENT(INOUT)        :: SG
type(EMsoft_T),INTENT(INOUT)            :: EMsoft

type(IO_T)                              :: Message
type(Timing_T)                          :: Timer

type(HDF_T)                             :: me  ! local HDF class

character(11)                           :: dstr
character(15)                           :: tstr
character(fnlen)                        :: progname = 'EMmkxtal.f90', groupname, dataset, xtalname, strings(1)
integer(kind=irg)                       :: hdferr, setting
real(kind=dbl)                          :: cellparams(6)
integer(kind=irg),allocatable           :: atomtypes(:)
real(kind=sgl),allocatable              :: atompos(:,:)
logical                                 :: openHDFfile

Message = IO_T()
Timer = Timing_T()

dstr = Timer%getDateString()
tstr = Timer%getTimeString()

me = HDF_T()

  xtalname = trim(EMsoft%generateFilePath('EMXtalFolderpathname',self%xtalname))
  hdferr =  me%createFile(xtalname)
  call me%error_check('SaveDataHDF:HDF_createFile:'//trim(xtalname), hdferr)

  groupname = SC_CrystalData
  hdferr = me%createGroup(groupname)
  call me%error_check( 'SaveDataHDF:HDF_createGroup:'//trim(groupname), hdferr)

  dataset = SC_ProgramName
  strings(1) = trim(progname)
  hdferr = me%writeDatasetStringArray(dataset, strings, 1)
  call me%error_check( 'SaveDataHDF:writeDatasetStringArray:'//trim(dataset), hdferr)

  dataset = SC_CreationDate
  strings(1) = dstr
  hdferr = me%writeDatasetStringArray(dataset, strings, 1)
  call me%error_check( 'SaveDataHDF:writeDatasetStringArray:'//trim(dataset), hdferr)

  dataset = SC_CreationTime
  strings(1) = tstr
  hdferr = me%writeDatasetStringArray(dataset, strings, 1)
  call me%error_check( 'SaveDataHDF:writeDatasetStringArray:'//trim(dataset), hdferr)

  dataset = SC_Creator
  strings(1) = trim(EMsoft%getConfigParameter('UserName'))
  hdferr = me%writeDatasetStringArray(dataset, strings, 1)
  call me%error_check( 'SaveDataHDF:writeDatasetStringArray:'//trim(dataset), hdferr)

  dataset = SC_CrystalSystem
  hdferr = me%writeDatasetInteger(dataset, self%xtal_system)
  call me%error_check( 'SaveDataHDF:writeDatasetStringArray:'//trim(dataset), hdferr)

  dataset = SC_LatticeParameters
  cellparams = (/ self%a, self%b, self%c, self%alpha, self%beta, self%gamma /)
  hdferr = me%writeDatasetDoubleArray(dataset, cellparams, 6)
  call me%error_check( 'SaveDataHDF:writeDatasetDoubleArray:'//trim(dataset), hdferr)

  dataset = SC_SpaceGroupNumber
  hdferr = me%writeDatasetInteger(dataset, SG%getSpaceGroupNumber())
  call me%error_check( 'SaveDataHDF:writeDatasetInteger:'//trim(dataset), hdferr)

  if (SG%getuseHallSG().eqv..TRUE.) then 
    dataset = SC_HallSGnumber
    hdferr = me%writeDatasetInteger(dataset, SG%getHallSpaceGroupNumber())
    call me%error_check( 'SaveDataHDF:writeDatasetInteger:'//trim(dataset), hdferr)

    dataset = SC_HallSG
    strings(1) = trim(get_HallString( SG%getHallSpaceGroupNumber() ) )
    hdferr = me%writeDatasetStringArray(dataset, strings, 1)
    call me%error_check( 'SaveDataHDF:writeDatasetStringArray:'//trim(dataset), hdferr)
  end if 

  ! make sure we do not write a '0' for the SGset variable; it must be either 1 or 2
  if (SG%getSpaceGroupSetting().eq.0) then
    setting = 1
  else
    setting = SG%getSpaceGroupSetting()
  end if
  dataset = SC_SpaceGroupSetting
  hdferr = me%writeDatasetInteger(dataset, setting)
  call me%error_check( 'SaveDataHDF:writeDatasetInteger:'//trim(dataset), hdferr)

  dataset = SC_Natomtypes
  hdferr = me%writeDatasetInteger(dataset, self%ATOM_ntype)
  call me%error_check( 'SaveDataHDF:writeDatasetInteger:'//trim(dataset), hdferr)

  allocate(atomtypes(self%ATOM_ntype))
  atomtypes(1:self%ATOM_ntype) = self%ATOM_type(1:self%ATOM_ntype)
  dataset = SC_Atomtypes
  hdferr = me%writeDatasetIntegerArray(dataset, atomtypes, self%ATOM_ntype)
  call me%error_check( 'SaveDataHDF:writeDatasetIntegerArray:'//trim(dataset), hdferr)
  deallocate(atomtypes)

  allocate(atompos(self%ATOM_ntype,5))
  atompos(1:self%ATOM_ntype,1:5) = self%ATOM_pos(1:self%ATOM_ntype,1:5)
  dataset = SC_AtomData
  hdferr = me%writeDatasetFloatArray(dataset, atompos, self%ATOM_ntype, 5)
  call me%error_check( 'SaveDataHDF:writeDatasetFloatArray:'//trim(dataset), hdferr)
  deallocate(atompos)

  dataset = SC_Source
  strings(1) = trim(self%source)
  hdferr = me%writeDatasetStringArray(dataset, strings, 1)
  call me%error_check( 'SaveDataHDF:writeDatasetStringArray:'//trim(dataset), hdferr)

  call me%popall()

end subroutine saveDataHDF_


!--------------------------------------------------------------------------
recursive subroutine addXtalDataGroup_(self, SG, EMsoft, useHDF)
!DEC$ ATTRIBUTES DLLEXPORT :: addXtalDataGroup_
  !! author: MDG
  !! version: 1.0
  !! date: 01/14/20
  !!
  !! save crystal structure data to an HDF file

use mod_io
use mod_global
use mod_EMsoft
use HDF5
use mod_HDFsupport
use mod_timing
use mod_symmetry
use stringconstants
use mod_HallSG

IMPLICIT NONE

class(Cell_T),INTENT(INOUT)             :: self
type(SpaceGroup_T),INTENT(INOUT)        :: SG
type(EMsoft_T),INTENT(INOUT)            :: EMsoft
type(HDF_T),INTENT(INOUT)               :: useHDF

type(IO_T)                              :: Message
type(Timing_T)                          :: Timer

character(11)                           :: dstr
character(15)                           :: tstr
character(fnlen)                        :: progname = 'EMmkxtal.f90', groupname, dataset, xtalname, strings(1)
integer(kind=irg)                       :: hdferr, setting, HSGn
real(kind=dbl)                          :: cellparams(6)
integer(kind=irg),allocatable           :: atomtypes(:)
real(kind=sgl),allocatable              :: atompos(:,:)
logical                                 :: openHDFfile

Message = IO_T()
Timer = Timing_T()

dstr = Timer%getDateString()
tstr = Timer%getTimeString()

  groupname = SC_CrystalData
  hdferr = useHDF%createGroup(groupname)
  call useHDF%error_check( 'SaveDataHDF:HDF_createGroup:'//trim(groupname), hdferr)

  dataset = SC_ProgramName
  strings(1) = trim(progname)
  hdferr = useHDF%writeDatasetStringArray(dataset, strings, 1)
  call useHDF%error_check( 'SaveDataHDF:writeDatasetStringArray:'//trim(dataset), hdferr)

  dataset = SC_CreationDate
  strings(1) = dstr
  hdferr = useHDF%writeDatasetStringArray(dataset, strings, 1)
  call useHDF%error_check( 'SaveDataHDF:writeDatasetStringArray:'//trim(dataset), hdferr)

  dataset = SC_CreationTime
  strings(1) = tstr
  hdferr = useHDF%writeDatasetStringArray(dataset, strings, 1)
  call useHDF%error_check( 'SaveDataHDF:writeDatasetStringArray:'//trim(dataset), hdferr)

  dataset = SC_Creator
  strings(1) = trim(EMsoft%getConfigParameter('UserName'))
  hdferr = useHDF%writeDatasetStringArray(dataset, strings, 1)
  call useHDF%error_check( 'SaveDataHDF:writeDatasetStringArray:'//trim(dataset), hdferr)

  dataset = SC_CrystalSystem
  hdferr = useHDF%writeDatasetInteger(dataset, self%xtal_system)
  call useHDF%error_check( 'SaveDataHDF:writeDatasetStringArray:'//trim(dataset), hdferr)

  dataset = SC_LatticeParameters
  cellparams = (/ self%a, self%b, self%c, self%alpha, self%beta, self%gamma /)
  hdferr = useHDF%writeDatasetDoubleArray(dataset, cellparams, 6)
  call useHDF%error_check( 'SaveDataHDF:writeDatasetDoubleArray:'//trim(dataset), hdferr)

  dataset = SC_SpaceGroupNumber
  hdferr = useHDF%writeDatasetInteger(dataset, SG%getSpaceGroupNumber())
  call useHDF%error_check( 'SaveDataHDF:writeDatasetInteger:'//trim(dataset), hdferr)

  ! make sure we do not write a '0' for the SGset variable; it must be either 1 or 2
  if (SG%getSpaceGroupSetting().eq.0) then
    setting = 1
  else
    setting = SG%getSpaceGroupSetting()
  end if
  dataset = SC_SpaceGroupSetting
  hdferr = useHDF%writeDatasetInteger(dataset, setting)
  call useHDF%error_check( 'SaveDataHDF:writeDatasetInteger:'//trim(dataset), hdferr)

  if (SG%getuseHallSG().eqv..TRUE.) then 
    dataset = SC_HallSGnumber
    HSGn = SG%getHallSpaceGroupNumber()
    hdferr = useHDF%writeDatasetInteger(dataset, HSGn )
    call useHDF%error_check( 'SaveDataHDF:writeDatasetInteger:'//trim(dataset), hdferr)

    dataset = SC_HallSG
    strings(1) = trim(get_HallString( HSGn ) )
    hdferr = useHDF%writeDatasetStringArray(dataset, strings, 1)
    call useHDF%error_check( 'SaveDataHDF:writeDatasetStringArray:'//trim(dataset), hdferr)
  end if 

  dataset = SC_Natomtypes
  hdferr = useHDF%writeDatasetInteger(dataset, self%ATOM_ntype)
  call useHDF%error_check( 'SaveDataHDF:writeDatasetInteger:'//trim(dataset), hdferr)

  allocate(atomtypes(self%ATOM_ntype))
  atomtypes(1:self%ATOM_ntype) = self%ATOM_type(1:self%ATOM_ntype)
  dataset = SC_Atomtypes
  hdferr = useHDF%writeDatasetIntegerArray(dataset, atomtypes, self%ATOM_ntype)
  call useHDF%error_check( 'SaveDataHDF:writeDatasetIntegerArray:'//trim(dataset), hdferr)
  deallocate(atomtypes)

  allocate(atompos(self%ATOM_ntype,5))
  atompos(1:self%ATOM_ntype,1:5) = self%ATOM_pos(1:self%ATOM_ntype,1:5)
  dataset = SC_AtomData
  hdferr = useHDF%writeDatasetFloatArray(dataset, atompos, self%ATOM_ntype, 5)
  call useHDF%error_check( 'SaveDataHDF:writeDatasetFloatArray:'//trim(dataset), hdferr)
  deallocate(atompos)

  dataset = SC_Source
  strings(1) = trim(self%source)
  hdferr = useHDF%writeDatasetStringArray(dataset, strings, 1)
  call useHDF%error_check( 'SaveDataHDF:writeDatasetStringArray:'//trim(dataset), hdferr)

  call useHDF%pop()

end subroutine addXtalDataGroup_

!--------------------------------------------------------------------------
recursive subroutine getCrystalData_(self, xtalname, SG, EMsoft, verbose, useHDF)
!DEC$ ATTRIBUTES DLLEXPORT :: getCrystalData_
  !! author: MDG
  !! version: 1.0
  !! date: 01/09/20
  !!
  !! load or generate crystal data

use mod_EMsoft
use mod_io
use mod_symmetry
use mod_HDFsupport
use mod_HallSG

IMPLICIT NONE

class(Cell_T),INTENT(INOUT)             :: self
character(fnlen),INTENT(IN)             :: xtalname
type(SpaceGroup_T),INTENT(INOUT)        :: SG
type(EMsoft_T),INTENT(INOUT)            :: EMsoft
logical,INTENT(IN),OPTIONAL             :: verbose
type(HDF_T),OPTIONAL,INTENT(INOUT)      :: useHDF

type(HDF_T)                             :: me
type(IO_T)                              :: Message
integer(kind=irg)                       :: i, ipg, SGnum

self%xtalname = trim(xtalname)
if (present(useHDF)) then
  call self%readDataHDF(SG, EMsoft, useHDF)
else
  call self%readDataHDF(SG, EMsoft)
end if

call SG%setSpaceGrouphexset(.FALSE.)
if (SG%getSpaceGroupXtalSystem().eq.4)  call SG%setSpaceGrouphexset(.TRUE.)
if ((SG%getSpaceGroupXtalSystem().eq.5).AND.(SG%getSpaceGroupSetting().ne.2)) call SG%setSpaceGrouphexset(.TRUE.)

call self%calcPositions(SG,'s')

! and print the information on the screen
if (present(verbose)) then
 if (verbose) then
   call self%dumpXtalInfo(SG)
 end if
end if

end subroutine getCrystalData_

!--------------------------------------------------------------------------
recursive subroutine readDataHDF_(self, SG, EMsoft, useHDF, useXtalName)
!DEC$ ATTRIBUTES DLLEXPORT :: readDataHDF_
  !! author: MDG
  !! version: 1.0
  !! date: 01/09/20
  !!
  !! Read crystal structure data from an HDF file

use mod_EMsoft
use mod_io
use mod_symmetry
use HDF5
use mod_HDFsupport
use ISO_C_BINDING
use stringconstants
use mod_HallSG

IMPLICIT NONE

class(Cell_T),INTENT(INOUT)             :: self
type(SpaceGroup_T),INTENT(INOUT)        :: SG
type(EMsoft_T),INTENT(INOUT)            :: EMsoft
type(HDF_T),OPTIONAL,INTENT(INOUT)      :: useHDF
character(fnlen),OPTIONAL,INTENT(IN)    :: useXtalName

type(HDF_T)                             :: me

type(IO_T)                              :: Message
character(fnlen)                        :: dataset, groupname, xtalname
integer(HSIZE_T)                        :: dims(1), dims2(2)
integer(kind=irg)                       :: hdferr, nlines, xtal_system, SGnum, setting, &
                                           ATOM_ntype, i, HallSGnum
real(kind=dbl),allocatable              :: cellparams(:)
integer(kind=irg),allocatable           :: atomtypes(:)
real(kind=sgl),allocatable              :: atompos(:,:)
character(fnlen)                        :: pp
logical                                 :: openHDFfile, d_exists
character(fnlen, KIND=c_char),allocatable,TARGET    :: stringarray(:)
character(16)                           :: Hallstring 

openHDFfile = .TRUE.

! call openFortranHDFInterface()

if (present(useHDF)) then
  openHDFfile = .FALSE.
  me = useHDF
else
  me = HDF_T()
end if

if (present(useXtalName)) then
  xtalname = trim(useXtalName)
else
  xtalname = trim(EMsoft%generateFilePath('EMXtalFolderpathname',self%xtalname))
end if
hdferr =  me%openFile(xtalname)

groupname = SC_CrystalData
hdferr = me%openGroup(groupname)
call me%error_check('readDataHDF:HDF_openGroup:'//trim(groupname), hdferr)

dataset = SC_CrystalSystem
call me%readDatasetInteger(dataset, hdferr, xtal_system)
call SG%setSpaceGroupXtalSystem(xtal_system)
self%xtal_system = xtal_system
call me%error_check('readDataHDF:readDatasetInteger:'//trim(dataset), hdferr)

dataset = SC_LatticeParameters
call me%readDatasetDoubleArray(dataset, dims, hdferr, cellparams)
call me%error_check('readDataHDF:readDatasetDoubleArray:'//trim(dataset), hdferr)

self%a = cellparams(1)
self%b = cellparams(2)
self%c = cellparams(3)
self%alpha = cellparams(4)
self%beta = cellparams(5)
self%gamma = cellparams(6)

! compute the metric matrices
call self%CalcMatrices()

! check for the Hall Space Group parameter; if it is present in the .xtal file,
! then we must use the Hall generators instead of the ones created from the 
! regular space group number.  
dataset = SC_SpaceGroupNumber
call me%readDatasetInteger(dataset, hdferr, SGnum)
call me%error_check('readDataHDF:readDatasetInteger:'//trim(dataset), hdferr)

dataset = SC_HallSGnumber
call H5Lexists_f(me%getobjectID(),trim(dataset),d_exists, hdferr)
if (d_exists) then
  call me%readDatasetInteger(dataset, hdferr, HallSGnum)
  call me%error_check('readDataHDF:readDatasetInteger:'//trim(dataset), hdferr)

  SG = SpaceGroup_T( SGnumber = SGnum, xtalSystem = xtal_system, setting = setting, &
                     dmt=self%dmt, rmt=self%rmt, useHall=.TRUE., HallSGnumber=HallSGnum )
! read the Hall Space Group symbol and figure out whether or not there is inversion symmetry
  dataset = SC_HallSG
  call me%readDatasetStringArray(dataset, nlines, hdferr, stringarray)
  Hallstring = trim(stringarray(1))
  deallocate(stringarray)
  if (Hallstring(1:1).eq.'-') call SG%setSpaceGroupCentro(.TRUE.)
else 
  dataset = SC_SpaceGroupSetting
  call me%readDatasetInteger(dataset, hdferr, setting)
  call me%error_check('readDataHDF:readDatasetInteger:'//trim(dataset), hdferr)

  if (setting.eq.0) setting = 1
  SG = SpaceGroup_T( SGnumber = SGnum, xtalSystem = xtal_system, setting = setting, &
                     dmt=self%dmt, rmt=self%rmt )
end if 

! here we also set the point group number, which is needed by many routines
if (SGnum.ge.221) then
  i = 32
else
  i=0
  do while (SGPG(i+1).le.SGnum)
    i = i+1
  end do
end if
call SG%setPGnumber(i)

dataset = SC_Natomtypes
call me%readDatasetInteger(dataset, hdferr, ATOM_ntype)
self%ATOM_ntype = ATOM_ntype
call me%error_check('readDataHDF:readDatasetInteger:'//trim(dataset), hdferr)

dataset = SC_Atomtypes
call me%readDatasetIntegerArray(dataset, dims, hdferr, atomtypes)
call me%error_check('readDataHDF:readDatasetIntegerArray:'//trim(dataset), hdferr)
self%ATOM_type(1:self%ATOM_ntype) = atomtypes(1:self%ATOM_ntype)
deallocate(atomtypes)

dataset = SC_AtomData
call me%readDatasetFloatArray(dataset, dims2, hdferr, atompos)
call me%error_check('readDataHDF:readDatasetFloatArray:'//trim(dataset), hdferr)
self%ATOM_pos(1:self%ATOM_ntype,1:5) = atompos(1:self%ATOM_ntype,1:5)
deallocate(atompos)

! the following data set does not exist in older .xtal files so we test for its presence
dataset = SC_Source
call H5Lexists_f(me%getobjectID(),trim(dataset),d_exists, hdferr)
if (d_exists) then
  call me%readDatasetStringArray(dataset, nlines, hdferr, stringarray)
  self%source = trim(stringarray(1))
  deallocate(stringarray)
else
  self%source = 'undefined'
  call Message%printMessage('readDataHDF Warning: there is no Source dataset in '//trim(xtalname))
end if

if (openHDFfile) then
  call me%popall()
else ! just close this group, but not the file
  call me%pop()
end if

! for trigonal space groups we need to set SYM_trigonal to .TRUE.
if ((SGnum.ge.143).and.(SGnum.le.167)) then
  call SG%setSpaceGrouptrigonal(.TRUE.)
else
  call SG%setSpaceGrouptrigonal(.FALSE.)
end if

! we have not yet implemented the rhombohedral setting of the trigonal
! space groups, so this needs to remain at .FALSE. always.
call SG%setSpaceGroupsecond(.FALSE.)

! if (openHDFfile) call me%closer()

end subroutine readDataHDF_

!--------------------------------------------------------------------------
recursive subroutine calcPositions_(self, SG, switch, numcells)
!DEC$ ATTRIBUTES DLLEXPORT :: calcPositions_
  !! author: MDG
  !! version: 1.0
  !! date: 01/09/20
  !!
  !! Compute all atom positions in the fundamental unit
  !! cell and translate to neighbouring cells if needed
  !! (used for structure drawings and structure factor computations)

use mod_io
use mod_symmetry

IMPLICIT NONE

class(Cell_T),INTENT(INOUT)     :: self
type(SpaceGroup_T),INTENT(INOUT):: SG
character(1),INTENT(IN)         :: switch ! if switch='m', then multiple unit cells, otherwise single cell
integer(kind=irg),INTENT(IN),OPTIONAL :: numcells(3)

type(IO_T)                      :: Message
logical                         :: inside                       ! auxiliary logical
integer(kind=irg)               :: i,j,k,l,mm,icnt,celln(3),ncells,n,kk,ier, io_int(3), &
                                   jstart, kstart, lstart       ! various auxiliary variables
real(kind=dbl)                  :: ff(3),sh(3)                  ! auxiliary variables
real(kind=sgl)                  :: r(3),g(3)                    ! auxiliary variables
real(kind=dbl),allocatable      :: ctmp(:,:)


! multiple cells ?
 if (present(numcells)) then 
   call SG%setSpaceGroupreduce(.FALSE.)
   celln = numcells
   sh = 0.5_dbl*celln+1.0_dbl
   ncells = (2*celln(1)+1)*(2*celln(2)+1)*(2*celln(3)+1)
   jstart = -celln(1)
   kstart = -celln(2)
   lstart = -celln(3)
 else
! make sure all coordinates are reduced to the fundamental unit cell
   call SG%setSpaceGroupreduce(.TRUE.)
   jstart = 1
   kstart = 1
   lstart = 1
   if (switch.eq.'m') then
    call Message%ReadValue('Number of unit cells in a, b and c direction ?: ', io_int,3)
    celln = io_int
    sh = 0.5_dbl*celln+1.0_dbl
    ncells = celln(1)*celln(2)*celln(3)
   else
  ! no, just one cell
    celln=0
    sh=1.0_dbl
    ncells = 1
   end if
 end if

! main loop
! first allocate the apos variable (contains CARTESIAN coordinates
! if switch is 'm', crystal coordinates otherwise)
 if (allocated(self%apos)) deallocate(self%apos)
 allocate (self%apos(self%ATOM_ntype, ncells * SG%getSpaceGroupMATnum(), 3),stat=ier)
 if (ier.ne.0) call Message%printError('calcPositions',' unable to allocate memory for array cell%apos')

 do i=1,self%ATOM_ntype

! for each atom in the asymmetric unit
  call SG%CalcOrbit(dble(self%ATOM_pos(i,1:3)),n,ctmp)
  self%numat(i)=n
  icnt=1

! replicate in all cells
  do j=jstart,celln(1)+1
   ff(1)=dble(j)
   do k=kstart,celln(2)+1
    ff(2)=dble(k)
    do l=lstart,celln(3)+1
     ff(3)=dble(l)
     do kk=1,self%numat(i)
      do mm=1,3
       r(mm)=ctmp(kk,mm)+ff(mm)-sh(mm)
      end do
      if (switch.eq.'m') then
! make sure the atom is actually inside the block of unit
! cells, or on one of the edges/faces/corners
       inside=.TRUE.
       do mm=1,3
        if (abs(r(mm)+sh(mm)).gt.(celln(mm)+1.0)) inside=.FALSE.
       end do
       if (inside) then
        call self%TransSpace(r,g,'d','c')
        do mm=1,3
         self%apos(i,icnt,mm)=g(mm)
        end do
        icnt=icnt+1
       end if
      else ! switch
! prepare for structure factor computation
       do mm=1,3
        self%apos(i,icnt,mm)=r(mm)
       end do
       icnt=icnt+1
      end if  ! switch
     end do ! kk
    end do ! l
   end do ! k
  end do ! j
  self%numat(i)=icnt-1
  deallocate(ctmp)

 end do  ! self%ATOM_type

end subroutine calcPositions_

!--------------------------------------------------------------------------
recursive subroutine convertfromRtoH_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: convertfromRtoH_
  !! author: MDG
  !! version: 1.0
  !! date: 12/14/22
  !!
  !! Convert lattice parameters and atom coordinates to the hexagonal unit cell.
  !! We use the obverse setting in all cases.

use mod_io

IMPLICIT NONE

class(Cell_T),INTENT(INOUT)  :: self

integer(kind=irg)            :: i
real(kind=dbl)               :: ar, calpha, pos(3)
real(kind=dbl),parameter     :: M(3,3) = reshape( (/ 2.D0/3.D0,  1.D0/3.D0,  1.D0/3.D0, &
                                                    -1.D0/3.D0,  1.D0/3.D0,  1.D0/3.D0, &
                                                    -1.D0/3.D0, -2.D0/3.D0,  1.D0/3.D0 /), (/3,3/) ) 

! change the lattice parameters to hexagonal
calpha = cos( self%alpha * dtor )
ar = self%a 

self%alpha = 90.D0 
self%beta = 90.D0
self%gamma = 120.D0
self%a = ar * sqrt( 2.D0 - 2.D0 * calpha )
self%b = self%a 
self%c = ar * sqrt( 3.D0 + 6.D0 * calpha )

! and then the atom coordinates 
do i=1,self%ATOM_ntype
  pos = dble(self%ATOM_pos(i,1:3))
  pos = matmul( M, pos )
  self%ATOM_pos(i,1:3) = sngl(pos)
end do 

end subroutine convertfromRtoH_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! miscellaneous crystallographic computations
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine ShortestG_(self, SG, k, gone, gtwo, isym)
!DEC$ ATTRIBUTES DLLEXPORT :: ShortestG_
  !! author: MDG
  !! version: 1.0
  !! date: 01/26/20
  !!
  !! determine the pair of shortest reciprocal lattice
  !! vectors for a given zone axis; used to draw diffraction
  !! patterns and in a bunch of other routines
  !!
  !! we look for 4 integer numbers which transform ga and gb
  !! simultaneously into two shortest possible vectors of the same zone;
  !! a range of -10:10 should be sufficient as a search space
  !!
  !! If g1 and g2 are those shortest vectors, then we have
  !!
  !!    ga  =  na*g1 + ma*g2\n
  !!    gb  =  nb*g1 + mb*g2\n
  !!
  !! Inversion of this relation gives
  !!
  !!    g1  =  (mb*ga - ma*gb)/D\n
  !!    g2  =  (-nb*ga + na*gb)/D\n
  !!
  !! with D = na*mb - nb*ma.
  !!
  !! The procedure below searches for the combination of
  !! (ma,mb,na,nb) which simultaneously minimizes the
  !! length of g1 and g2, and makes sure that both g1 and g2
  !! are integer linear combinations of the reciprocal basis
  !! vectors.

use mod_io
use mod_symmetry

IMPLICIT NONE

class(Cell_T),INTENT(INOUT)             :: self
type(SpaceGroup_T),INTENT(INOUT)        :: SG
integer(kind=irg),INTENT(INOUT)         :: isym
 !! used to resolve some potential ambiguities for 3-fold and 6-fold symmetries
integer(kind=irg),INTENT(IN)            :: k(3)
 !! input zone axis indices
integer(kind=irg),INTENT(OUT)           :: gone(3)
 !! output first vector
integer(kind=irg),INTENT(OUT)           :: gtwo(3)
 !! output second vector

type(IO_T)                              :: Message
integer(kind=irg)                       :: ga(3),gb(3),nzero(3),u,v,w,snz,ml(4),igsave(3)
integer(kind=irg)                       :: ima,imb,ina,inb,el(6),denom,minsum,inm,il(48),jcnt,num
real(kind=sgl)                          :: fel(6),fit,gsave(3)
integer(kind=irg),allocatable           :: ifit(:,:,:,:)
integer(kind=irg),allocatable           :: itmp(:,:)

 u = k(1)
 v = k(2)
 w = k(3)

! determine two arbitrary vectors normal to k
! first count the zeroes in k
 nzero = (/0,0,0/)
 where (k.eq.0) nzero = 1
 snz = sum(nzero)

 if (snz.eq.0) then  ! make sure ga x gb is parallel to k
   ga = (/v,-u,0/)
   gb = (/w,0,-u/)
 else
  select case (snz)
  case(1);  ga = nzero; gb = 1 - nzero
            if (nzero(1).eq.1) gb = gb * (/0,w,-v/)
            if (nzero(2).eq.1) gb = gb * (/w,0,-u/)
            if (nzero(3).eq.1) gb = gb * (/v,-u,0/)
  case(2);  if ((nzero(1).eq.1).and.(nzero(2).eq.1)) then
             ga = (/1,0,0/); gb = (/0,1,0/)
            endif
            if ((nzero(1).eq.1).and.(nzero(3).eq.1)) then
             ga = (/0,0,1/); gb = (/1,0,0/)
            endif
            if ((nzero(2).eq.1).and.(nzero(3).eq.1)) then
             ga = (/0,1,0/); gb = (/0,0,1/)
            endif
  case(3); call Message%printError('ShortestG',' beam direction cannot be [0,0,0]')
  end select
 end if

! check linear combinations to see if there are any shorter ones
 inm = 10
 allocate(ifit(-inm:inm,-inm:inm,-inm:inm,-inm:inm))
 do ima=-inm,inm
  do imb=-inm,inm
   do ina=-inm,inm
    do inb=-inm,inm
     el(1) = imb*ga(1)-ima*gb(1)
     el(2) = imb*ga(2)-ima*gb(2)
     el(3) = imb*ga(3)-ima*gb(3)
     el(4) = ina*gb(1)-inb*ga(1)
     el(5) = ina*gb(2)-inb*ga(2)
     el(6) = ina*gb(3)-inb*ga(3)
     denom = ina*imb-inb*ima
     ifit(ima,imb,ina,inb)=100
     if (denom.ne.0) then
      fel = float(el)/float(denom)
      fit = sum(abs(float(int(fel))-fel))
! here is where we only keep the integer combinations
      if (fit.eq.0.0) then
        gone(1:3) = int(fel(1:3))
        gtwo(1:3) = int(fel(4:6))
! keep the sum of the squares of the lengths
       ifit(ima,imb,ina,inb)=sum(gone**2)+sum(gtwo**2)
      end if
     end if
    end do
   end do
  end do
 end do
 minsum = 50

! look for the minimum of ifit with the smallest and most
! positive coefficients; store them in ml
! [minloc does not work here because there may be multiple minima]
 do ima=-inm,inm
  do imb=-inm,inm
   do ina=-inm,inm
    do inb=-inm,inm
     if (ifit(ima,imb,ina,inb).le.minsum) then
      minsum = ifit(ima,imb,ina,inb)
      ml(1) = ima
      ml(2) = imb
      ml(3) = ina
      ml(4) = inb
     end if
    end do
   end do
  end do
 end do
 deallocate(ifit)

! transform ga and gb into g1 and g2
 gone = (ml(2)*ga-ml(1)*gb)/(ml(3)*ml(2)-ml(4)*ml(1))
 gtwo = (ml(3)*gb-ml(4)*ga)/(ml(3)*ml(2)-ml(4)*ml(1))

! next rank these two vectors so that their cross product is along +k
 call self%CalcCross(float(gone),float(gtwo),gsave,'r','r',0)
 fit = self%CalcDot(gsave,float(k),'r')
 if (fit.lt.0.0) then
  igsave = gone
  gone = gtwo
  gtwo = igsave
 end if

! finally, if isym.ne.0 make sure that the selection of the
! basis vectors for the 3-fold and 6-fold 2D point groups is
! correctly done.
!
! For isym < 7:   90 degrees between vectors (should not be a problem)
! For isym = 7:  120 degrees between vectors (should be ok too)
! distinguish between 3 coming from a cubic group
! vs. the same symmetries originating from a hexagonal setting
 if ((isym.eq.7).and.(self%gamma.eq.120.0)) then
   isym=11
   fit = self%CalcAngle(float(gone),float(gtwo),'r')/dtor
   if (abs(fit-120.0).lt.1.0) then
     gtwo=gone+gtwo
   end if
 end if

! For isym = 8:  here we should distinguish between the settings 3m1 and 31m !!!
!                The angle should always be 120 degrees, so we must check that
!                this is the case for the selected gone and gtwo.
 if ((isym.eq.8).and.(self%gamma.eq.120.0)) then
   isym=12
   fit = self%CalcAngle(float(gone),float(gtwo),'r')/dtor
   if (abs(fit-120.0).lt.1.0) then
     gtwo=gone+gtwo
   end if
 end if
!
 if (isym.eq.8) then
   fit = self%CalcAngle(float(gone),float(gtwo),'r')/dtor
   if (abs(fit-120.0).gt.1.0) then
     gtwo=gtwo-gone
   end if

! we assume it is the 31m setting;  if the order of gone is 6, then that is not true
   call SG%CalcFamily(gone,num,'r',itmp)
   call SG%GetOrderinZone(float(k),il,num,jcnt,itmp)
   if (jcnt.eq.6) then  ! it is the 3m1 setting
     isym = 13
   end if
 end if

! it could the 3m1 setting for the 3m hexagonal case
 if (isym.eq.12) then
! we assume it is the 31m setting;  if the order of gone is 6, then that is not true
   call SG%CalcFamily(gone,num,'r',itmp)
   call SG%GetOrderinZone(float(k),il,num,jcnt,itmp)
   if (jcnt.eq.6) then  ! it is the 3m1 setting
     isym = 14
   end if
 end if

! For isym = 9 or 10:   60 degrees between vectors (may not be the case)
 if ((isym.eq.9).or.(isym.eq.10)) then
   fit = self%CalcAngle(float(gone),float(gtwo),'r')/dtor
   if (abs(fit-120.0).lt.1.0) then
     gtwo=gone+gtwo
   end if
 end if

end subroutine ShortestG_

! !--------------------------------------------------------------------------
! !
! ! FUNCTION:Convert_kgs_to_Substrate
! !
! !> @author Marc De Graef, Carnegie Mellon University
! !
! !> @brief convert a film wave vector into substrate coordinates
! !
! !> @param cell film unit cell pointer
! !> @param cellS substrate unit cell pointer
! !> @param kg wave vector components in film
! !> @param TTinv coordinate transformation matrix
! !> @param lambdaS electron wavelength in substrate, corrected for refraction
! !> @param FN foil normal in substrate frame
! !
! !> @date 11/25/14 MDG 1.0 original
! !> @date 08/19/15 SS  1.1 replaced FN by r in evaluation of tangential
! !> component and corrected calculation of kgS from p1
! !> @date 09/23/15 commented out a few lines
! !--------------------------------------------------------------------------
! recursive function Convert_kgs_to_Substrate(cell, cellS, kg, TTinv, FN) result(kgS)
! !DEC$ ATTRIBUTES DLLEXPORT :: Convert_kgs_to_Substrate

! use local
! use typedefs

! IMPLICIT NONE

! type(unitcell)                          :: cell, cellS
! real(kind=sgl),INTENT(IN)               :: kg(3)
! real(kind=sgl),INTENT(IN)               :: TTinv(3,3)
! real(kind=sgl),INTENT(IN)               :: FN(3)
! real(kind=sgl)                          :: kgS(3)
! real(kind=sgl)                          :: dp,normal
! real(kind=sgl)                          :: tangential(3)

! real(kind=sgl)                          :: p(3), r(3), p1(3)

! ! convert to direct space for transforming to substrate frame
! call TransSpace(cell,kg,p,'r','d')

! ! convert to substrate frame
! p = matmul(TTinv,p)

! ! convert to cartesian frame and get correct length of wavevector
! call TransSpace(cellS,p,p1,'d','c')

! !call NormVec(cellS,p1,'c')
! !p1 = p1/cell%mLambda

! ! convert foil normal to cartesian frame and get normal component of wavevector
! call TransSpace(cellS,FN,r,'r','c')
! call NormVec(cellS,r,'c')
! dp = CalcDot(cellS,p1,r,'c')

! ! subtract out the normal component to get tangential component; this is conserved across the interface
! tangential = p1 - dp*r

! ! get magnitude of normal component
! normal = sqrt((1.0/cellS%mLambda)**2 - sum(tangential**2))

! p1 = normal*r + tangential

! call NormVec(cellS,p1,'c')

! p1 = p1/cellS%mLambda

! call TransSpace(cellS,p1,kgS,'c','r')

! end function Convert_kgs_to_Substrate

end module mod_crystallography
