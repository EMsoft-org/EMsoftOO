! ###################################################################
! Copyright (c) 2014-2020, Marc De Graef Research Group/Carnegie Mellon University
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
      character(fnlen)                     :: fname
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
      real(kind=dbl),allocatable           :: apos(:,:,:)
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
        private
! routines for crystallographic computations
          procedure, pass(self) :: computeMatrices
          procedure, pass(self) :: transSpaceSingle
          procedure, pass(self) :: transSpaceDouble
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
! routines to set class parameters 
          procedure, pass(self) :: setWyckoff_
          procedure, pass(self) :: setFileName_
          procedure, pass(self) :: setSource_
          procedure, pass(self) :: setXtalSystem_
          procedure, pass(self) :: requestLatticeParameters
          procedure, pass(self) :: getLatticeParameterSingle
          procedure, pass(self) :: getLatticeParametersAll
          procedure, pass(self) :: getAsymmetricPosition
          procedure, pass(self) :: displayPeriodicTable
          ! procedure, pass(self) :: extractAtomPositionData
          procedure, pass(self) :: calcTheoreticalDensity
! routines to get class parameters 

! routines to read/write .xtal files 

! miscellaneous routines 
          procedure, pass(self) :: resetUnitCell
          procedure, pass(self) :: GetAsymPosWyckoff_


          ! procedure, pass(self) :: GetOR
          ! procedure, pass(self) :: ComputeOR
          ! procedure, pass(self) :: CalcsgHOLZ
          ! procedure, pass(self) :: GetHOLZGeometry
          ! procedure, pass(self) :: GetHOLZcoordinates
          ! procedure, pass(self) :: Convert_kgs_to_Substrate
          ! procedure, pass(self) :: setEMsoftXtalSystem
          ! procedure, pass(self) :: getEMsoftXtalSystem

          generic, public :: resetCell => resetUnitCell
          !DEC$ ATTRIBUTES DLLEXPORT :: resetCell
          generic, public :: calcMatrices => computeMatrices
          !DEC$ ATTRIBUTES DLLEXPORT :: calcMatrices 
          generic, public :: transSpace => transSpaceDouble, transSpaceSingle
          !DEC$ ATTRIBUTES DLLEXPORT :: transSpace 
          generic, public :: transCoor => transformCoordinates
          !DEC$ ATTRIBUTES DLLEXPORT :: transCoor 
          generic, public :: calcDot => calcDotSingle, calcDotDouble
          !DEC$ ATTRIBUTES DLLEXPORT :: calcDot 
          generic, public :: normVec => normVecSingle, normVecDouble
          !DEC$ ATTRIBUTES DLLEXPORT :: normVec 
          generic, public :: calcLength => calcLengthSingle, calcLengthDouble
          !DEC$ ATTRIBUTES DLLEXPORT :: calcLength 
          generic, public :: calcAngle => calcAngleSingle, calcAngleDouble
          !DEC$ ATTRIBUTES DLLEXPORT :: calcAngle 
          generic, public :: calcCross => calcCrossSingle, calcCrossDouble
          !DEC$ ATTRIBUTES DLLEXPORT :: calcCross 
          ! generic, public :: MilBrav => convertMilBrav
          !DEC$ ATTRIBUTES DLLEXPORT :: MilBrav 
          generic, public :: getLatParm => getLatticeParameterSingle, getLatticeParametersAll
          !DEC$ ATTRIBUTES DLLEXPORT :: getLatParm 
          generic, public :: setLatParm => requestLatticeParameters
          !DEC$ ATTRIBUTES DLLEXPORT :: setLatParm 
          generic, public :: setXtalSystem => setXtalSystem_
          !DEC$ ATTRIBUTES DLLEXPORT :: setXtalSystem 
          generic, public :: getAsymPos => getAsymmetricPosition
          !DEC$ ATTRIBUTES DLLEXPORT :: getAsymPos 
          ! generic, public :: displayElements => displayPeriodicTable 
          !DEC$ ATTRIBUTES DLLEXPORT :: displayElements 
          ! generic, public :: extractposition => extractAtomPositionData 
          !DEC$ ATTRIBUTES DLLEXPORT :: extractposition 
          generic, public :: calcDensity => calcTheoreticalDensity 
          !DEC$ ATTRIBUTES DLLEXPORT :: calcDensity 
          generic, public :: GetAsymPosWyckoff => GetAsymPosWyckoff_
          !DEC$ ATTRIBUTES DLLEXPORT :: GetAsymPosWyckoff 
          generic, public :: setWyckoff => setWyckoff_
          !DEC$ ATTRIBUTES DLLEXPORT :: setWyckoff 
          generic, public :: setFileName => setFileName_
          !DEC$ ATTRIBUTES DLLEXPORT :: setFileName 
          generic, public :: setSource => setSource_
          !DEC$ ATTRIBUTES DLLEXPORT :: setSource 

  end type Cell_T

! the following parameters and arrays need to be moved to somewhere else ... 
      ! real(kind=sgl),allocatable           :: scatfacg(:)
      ! complex(kind=sgl),allocatable        :: scatfac(:,:) 
      ! complex(kind=dbl),allocatable        :: LUT(:,:,:), SghLUT(:,:,:,:)
      ! complex(kind=dbl),allocatable        :: LUTqg(:,:,:)
      ! logical,allocatable                  :: dbdiff(:,:,:)
      ! type(symdata)                        :: SG
      ! type(reflisttype),pointer            :: reflist
      ! type(reflisttype),pointer            :: firstw                ! connection to first weak entry in linked list
      ! integer(kind=irg)                    :: DynNbeams, DynNbeamsLinked, nns, nnw, numscatfac
      ! real(kind=dbl)                       :: voltage, mLambda, mRelcor, mSigma, mPsihat   ! voltage always in keV !

  ! the constructor routine for this class 
  interface Cell_T
    module procedure Cell_constructor
  end interface Cell_T

contains

!--------------------------------------------------------------------------
type(Cell_T) function Cell_constructor( latparm ) result(cell)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/07/20
  !!
  !! constructor for the Cell Class 
  
IMPLICIT NONE

real(kind=dbl), INTENT(IN), OPTIONAL :: latparm(6)
 !! optional lattice parameter array; cell will be reset to zero if not present 

call cell%resetUnitCell()

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
end if

end function Cell_constructor

!--------------------------------------------------------------------------
recursive subroutine resetUnitCell(self)
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
 self%fname = ''
 self%source = ''

! and deallocate any arrays
 if (allocated(self%apos)) deallocate(self%apos)
 
end subroutine resetUnitCell

!--------------------------------------------------------------------------
recursive function getLatticeParameterSingle( self, p ) result(lp)
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
recursive function getLatticeParametersAll( self ) result(lp)
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
recursive subroutine computeMatrices(self)
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
recursive subroutine DumpXtalInfo(self, SG)    
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/12/20
  !!
  !! a brief summary of the crystal structure

use mod_io
use mod_symmetry

IMPLICIT NONE

type(Cell_T), intent(inout)             :: self
type(SpaceGroup_T), intent(inout)       :: SG 

type(IO_T)                              :: Message 

integer(kind=irg)                       :: i, j, oi_int(3)
real(kind=dbl)                          :: oi_real(5)


 call Message%printMessage('', frm = "(A/)")
 call Message%printMessage('Crystal Structure Information', frm = "('-->',A,'<--')")
 oi_real(1) = self%getLatParm('a')
 call WriteValue('  a [nm]             : ', oi_real, 1, "(F9.5)")
 oi_real(1) = self%getLatParm('b')
 call WriteValue('  b [nm]             : ', oi_real, 1, "(F9.5)")
 oi_real(1) = self%getLatParm('c')
 call WriteValue('  c [nm]             : ', oi_real, 1, "(F9.5)")
 oi_real(1) = self%getLatParm('alpha')
 call WriteValue('  alpha [deg]        : ', oi_real, 1, "(F9.5)")
 oi_real(1) = self%getLatParm('beta')
 call WriteValue('  beta  [deg]        : ', oi_real, 1, "(F9.5)")
 oi_real(1) = self%getLatParm('gamma')
 call WriteValue('  gamma [deg]        : ', oi_real, 1, "(F9.5)")
 oi_real(1) = self%vol
 call WriteValue('  Volume [nm^3]      : ', oi_real, 1, "(F12.8)")
 oi_int(1) = SG%getSpaceGroupNumber()
 call WriteValue('  Space group #      : ', oi_int, 1, "(1x,I3)")
 call WriteValue('  Space group symbol : ', trim(SYM_SGname(SG%getSpaceGroupNumber())) )
 call WriteValue('  Generator String   : ',  trim(SYM_GL(SG%getSpaceGroupNumber())) )
 if ((SG%getSpaceGroupSetting().eq.2).AND.(SG%getSpaceGroupXtalSystem().ne.5)) then 
  call Message%printMessage('   Using second origin setting', frm = "(A)")
 endif
 if ((SG%getSpaceGroupSetting().eq.2).AND.(SG%getSpaceGroupXtalSystem().eq.5)) then 
  call Message%printMessage('   Using rhombohedral parameters', frm = "(A)")
 endif
  if (SG%getSpaceGroupCentro()) then 
    call Message%printMessage('   Structure is centrosymmetric', frm = "(A)")
 else 
   call Message%printMessage('   Structure is non-centrosymmetric', frm = "(A)")
 end if

! space group and point group information
 oi_int(1) = SG%getSpaceGroupGENnum()
 call WriteValue('  # generators       : ', oi_int, 1, "(1x,I3)")
 oi_int(1) = SG%getSpaceGroupMATnum()
 call WriteValue('  # symmetry matrices: ', oi_int, 1, "(1x,I3)")
 oi_int(1) = SG%getSpaceGroupNUMpt()
 call WriteValue('  # point sym. matr. : ', oi_int, 1, "(1x,I3)")

! generate atom positions and dump output  
 call Message%printMessage('', frm = "(A/)")
 call CalcPositions(self,'v')
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

end subroutine DumpXtalInfo

!--------------------------------------------------------------------------
recursive subroutine setWyckoff_(self, useWyckoff)
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
recursive subroutine setFileName_(self, fname)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/13/20
  !!
  !! set the output filename 

IMPLICIT NONE 

class(Cell_T), INTENT(INOUT)    :: self 
character(fnlen),INTENT(IN)     :: fname 

self%fname = trim(fname)

end subroutine setFileName_

!--------------------------------------------------------------------------
recursive subroutine setSource_(self, source)
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
recursive subroutine setXtalSystem_(self, xs)
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
recursive subroutine getAsymmetricPosition(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/13/20
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
  call Message%ReadValue('', instring)
  sl = len(trim(instring))
  j = 0
  do i=1,sl
    if (instring(i:i).ne.' ') then
      j = j+1
      list(j) = instring(i:i)
    end if
  end do

! interpret this string and extract coordinates and such ...
  call extractAtomPositionData(list,pt) 
  
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
   d = matmul(self%rsm,t)
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
   d = matmul(self%rsm,t)
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
recursive subroutine MilBrav(p,q,d)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/07/20
  !!
  !! conversion from Miller to Miller-Bravais indices for
  !! directions.  The switch d is either '34' or '43'.
  !! implements equations 1.31 and 1.32, pages 24-25. 

IMPLICIT NONE

integer(kind=irg),INTENT(INOUT)         :: p(3)         
 !! input/output vector
!f2py intent(in,out) ::  p
integer(kind=irg),INTENT(INOUT)         :: q(4)         
 !! input/output vector
!f2py intent(in,out) ::  q
character(2),INTENT(IN)                 :: d            
 !! direction string ('34' or '43')

integer(kind=irg)                       :: i, j         
real(kind=sgl)                          :: r(4), rm, tmp(4)     

 if (d.eq.'43') then 
! equation 1.31
! these will always be integers, so no reduction is required
  p(1) = q(1)-q(3)
  p(2) = q(2)-q(3)
  p(3) = q(4)
 else
! equation 1.32
! there is no need to divide by 3, since that would be taken out 
! by the reduction to integers in the next step
  r(1) = float(2*p(1)-p(2))
  r(2) = float(2*p(2)-p(1))
  r(3) = -float(p(1)+p(2))
  r(4) = float(3*p(3))

! next reduce to common integers
! first, find the non-zero minimum index
  rm = 100.0
  do i=1,4 
   if ((abs(r(i)).lt.rm).and.(r(i).ne.0.0)) rm = abs(r(i))
  end do

! then check if this index is a common divider of the others
  j = 0
  do i=1,4
   tmp(i) = r(i)/rm
   if ( ( abs(tmp(i))-int(abs(tmp(i))) ).eq.0.0) j=j+1
  end do
  if (j.eq.4) then
    q = tmp
  else  
    q = r
  end if
 end if

end subroutine MilBrav

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! routines for crystallographic IO 
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine requestLatticeParameters(self, SG)
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

IMPLICIT NONE

class(Cell_T),intent(inout)             :: self
type(SpaceGroup_T), intent(inout)       :: SG 

type(IO_T)                              :: Message
integer(kind=irg)                       :: io_int(1)
real(kind=dbl)                          :: io_real(1)

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
    call Message%printMessage('Enter 1 for rhombohedral setting ,', frm = "(/A)")
    call Message%ReadValue('0 for hexagonal setting : ', io_int, 1)
    if (io_int(1).eq.0) then
     self%xtal_system=4   ! this is set to 4 so that we ask for the correct lattice parameters below
    else
     call SG%setSpaceGroupsecond(.TRUE.)
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
     call Message%ReadValue('    beta  [deg] = ', io_real, 1)
     self%beta = io_real(1)
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
recursive subroutine extractAtomPositionData(list,pt)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/07/20
  !!
  !! Extract the coordinates, site occupation, and DW factor from the 
  !! input string; note that coordinates can be entered in decimal or in fractional
  !! notation, hence the somewhat convoluted way of interpreting this string...

IMPLICIT NONE

character(1),INTENT(IN)                 :: list(256)                            
 !! input string
real(kind=sgl),INTENT(OUT)              :: pt(5)                                
 !! output real array

integer(kind=irg)                       :: comma(6),slash(5),period(5), &
                                           ccnt,scnt,pcnt,pp,i,j,hcnt, &
                                           ip,ipt,icnt,nd,n,k,ns                
integer(kind=irg),parameter             :: nmb(48:57)=(/0,1,2,3,4,5,6,7,8,9/)   
real(kind=dbl)                          :: nominator,denominator,x              
logical                                 :: hasperiod                            

! first, make sure all the spaces are removed from the list array

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
 do while (ccnt.lt.6) 
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
 if (pt(4).eq.0.0) pt(4) = 1.0

end subroutine extractAtomPositionData

!--------------------------------------------------------------------------
recursive subroutine calcTheoreticalDensity(self, Z2percent)
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


! !--------------------------------------------------------------------------
! !
! ! SUBROUTINE: CalcPositions
! !
! !> @author Marc De Graef, Carnegie Mellon University
! !
! !> @brief  compute atom positions in one or more unit cells
! !
! !> @details Compute all atom positions in the fundamental unit
! !> cell and translate to neighbouring cells if needed
! !> (used for structure drawings and structure factor computations)
! !!
! !> @param cell unit cell pointer
! !> @param switch input vector 
! !
! !> @date  10/13/98 MDG 1.0 original
! !> @date   5/19/01 MDG 2.0 f90
! !> @date  11/27/01 MDG 2.1 added kind support
! !> @date  03/21/13 MDG 3.0 clean up and updated IO
! !> @date  01/10/14 MDG 4.0 SG is now part of the unitcell type
! !> @date  06/05/14 MDG 4.1 made cell an argument instead of global variable; replaced itmp by argument
! !--------------------------------------------------------------------------
! recursive subroutine CalcPositions(cell,switch)
! !DEC$ ATTRIBUTES DLLEXPORT :: CalcPositions

! use mod_io

! IMPLICIT NONE

! type(unitcell)                  :: cell
! character(1),INTENT(IN)         :: switch                       !< if switch='m', then multiple unit cells, otherwise single cell

! logical                         :: inside                       !< auxiliary logical
! integer(kind=irg)               :: i,j,k,l,mm,icnt,celln(3),ncells,n,kk,ier, io_int(3)  !< various auxiliary variables
! real(kind=dbl)                  :: ctmp(192,3),ff(3),sh(3)      !< auxiliary variables  
! real(kind=sgl)                  :: r(3),g(3)                    !< auxiliary variables  

! ! make sure all coordinates are reduced to the fundamental unit cell
!  cell%SG%SYM_reduce=.TRUE.

! ! multiple cells ?
!  if (switch.eq.'m') then 
!   call ReadValue('Number of unit cells in a, b and c direction ?: ', io_int,3)
!   do j=1,3
!    celln(j) = io_int(j)
!    sh(j) = 0.5_dbl*celln(j)+1.0_dbl
!   end do
!   ncells = celln(1)*celln(2)*celln(3)
!  else
! ! no, just one cell
!   do j=1,3
!    celln(j)=0
!    sh(j)=1.0_dbl
!   end do
!   ncells = 1
!  end if

! ! main loop
! ! first allocate the apos variable (contains CARTESIAN coordinates
! ! if switch is 'm', crystal coordinates otherwise)
!  if (allocated(cell%apos)) deallocate(cell%apos)
!  allocate (cell%apos(cell%ATOM_ntype, ncells * cell%SG%SYM_MATnum, 3),stat=ier)
!  if (ier.ne.0) call FatalError('CalcPositions',' unable to allocate memory for array cell%apos')

!  do i=1,cell%ATOM_ntype

! ! for each atom in the asymmetric unit
!   call CalcOrbit(cell,i,n,ctmp)
!   cell%numat(i)=n
!   icnt=1

! ! replicate in all cells
!   do j=1,celln(1)+1
!    ff(1)=dble(j)
!    do k=1,celln(2)+1
!     ff(2)=dble(k)
!     do l=1,celln(3)+1
!      ff(3)=dble(l)
!      do kk=1,cell%numat(i)
!       do mm=1,3
!        r(mm)=ctmp(kk,mm)+ff(mm)-sh(mm)
!       end do 
!       if (switch.eq.'m') then
! ! make sure the atom is actually inside the block of unit
! ! cells, or on one of the edges/faces/corners
!        inside=.TRUE.
!        do mm=1,3
!         if ((r(mm)+sh(mm)).gt.(celln(mm)+1.0)) inside=.FALSE.
!        end do
!        if (inside) then
!         call TransSpace(cell,r,g,'d','c')
!         do mm=1,3
!          cell%apos(i,icnt,mm)=g(mm)
!         end do
!         icnt=icnt+1
!        end if
!       else ! switch

! ! prepare for structure factor computation
!        do mm=1,3
!         cell%apos(i,icnt,mm)=r(mm)
!        end do
!        icnt=icnt+1
!       end if  ! switch
!      end do ! kk
!     end do ! l 
!    end do ! k
!   end do ! j
!   cell%numat(i)=icnt-1
!  end do  ! cell%ATOM_type

! end subroutine CalcPositions

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! miscellaneous crystallographic computations
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

! !--------------------------------------------------------------------------
! recursive subroutine ShortestG_(self, k, gone, gtwo, isym)
!   !! author: MDG 
!   !! version: 1.0 
!   !! date: 01/10/20
!   !!
!   !! determine the pair of shortest reciprocal lattice 
!   !! vectors for a given zone axis; used to draw diffraction 
!   !! patterns and in a bunch of other routines
!   !!
!   !! we look for 4 integer numbers which transform ga and gb
!   !! simultaneously into two shortest possible vectors of the same zone;
!   !! a range of -10:10 should be sufficient as a search space
!   !! 
!   !! If g1 and g2 are those shortest vectors, then we have
!   !! 
!   !!    ga  =  na*g1 + ma*g2\n
!   !!    gb  =  nb*g1 + mb*g2\n
!   !! 
!   !! Inversion of this relation gives
!   !! 
!   !!    g1  =  (mb*ga - ma*gb)/D\n
!   !!    g2  =  (-nb*ga + na*gb)/D\n
!   !! 
!   !! with D = na*mb - nb*ma.
!   !! 
!   !! The procedure below searches for the combination of 
!   !! (ma,mb,na,nb) which simultaneously minimizes the 
!   !! length of g1 and g2, and makes sure that both g1 and g2
!   !! are integer linear combinations of the reciprocal basis
!   !! vectors.


! use error
! use crystal
! use constants

! IMPLICIT NONE

! class(SpaceGroup_T),INTENT(INOUT)       :: self
! integer(kind=irg),INTENT(INOUT)         :: isym                 
!  !! used to resolve some potential ambiguities for 3-fold and 6-fold symmetries
! integer(kind=irg),INTENT(IN)            :: k(3)                 
!  !! input zone axis indices
! integer(kind=irg),INTENT(OUT)           :: gone(3)              
!  !! output first vector
! integer(kind=irg),INTENT(OUT)           :: gtwo(3)              
!  !! output second vector

! integer(kind=irg)                       :: ga(3),gb(3),nzero(3),u,v,w,snz,ml(4),igsave(3)       ! auxiliary variables
! integer(kind=irg)                       :: ima,imb,ina,inb,el(6),denom,minsum,inm,il(48),jcnt,num ! auxiliary variables
! real(kind=sgl)                          :: fel(6),fit,gsave(3)  ! auxiliary variables
! integer(kind=irg),allocatable           :: ifit(:,:,:,:)        ! array used to search
! integer(kind=irg)                       :: itmp(48,3)           ! array used for family computations etc

!  u = k(1)
!  v = k(2)
!  w = k(3)

! ! determine two arbitrary vectors normal to k 
! ! first count the zeroes in k
!  nzero = (/0,0,0/)
!  where (k.eq.0) nzero = 1
!  snz = sum(nzero)
 
!  if (snz.eq.0) then  ! make sure ga x gb is parallel to k
!    ga = (/v,-u,0/)
!    gb = (/w,0,-u/)
!  else
!   select case (snz)
!   case(1);  ga = nzero; gb = 1 - nzero
!             if (nzero(1).eq.1) gb = gb * (/0,w,-v/)
!             if (nzero(2).eq.1) gb = gb * (/w,0,-u/)
!             if (nzero(3).eq.1) gb = gb * (/v,-u,0/)
!   case(2);  if ((nzero(1).eq.1).and.(nzero(2).eq.1)) then
!              ga = (/1,0,0/); gb = (/0,1,0/)
!             endif
!             if ((nzero(1).eq.1).and.(nzero(3).eq.1)) then
!              ga = (/0,0,1/); gb = (/1,0,0/)
!             endif
!             if ((nzero(2).eq.1).and.(nzero(3).eq.1)) then
!              ga = (/0,1,0/); gb = (/0,0,1/)
!             endif
!   case(3); call FatalError('ShortestG',' beam direction cannot be [0,0,0]')
!   end select
!  end if 

! ! check linear combinations to see if there are any shorter ones
!  inm = 10
!  allocate(ifit(-inm:inm,-inm:inm,-inm:inm,-inm:inm))
!  do ima=-inm,inm
!   do imb=-inm,inm
!    do ina=-inm,inm
!     do inb=-inm,inm
!      el(1) = imb*ga(1)-ima*gb(1) 
!      el(2) = imb*ga(2)-ima*gb(2) 
!      el(3) = imb*ga(3)-ima*gb(3) 
!      el(4) = ina*gb(1)-inb*ga(1) 
!      el(5) = ina*gb(2)-inb*ga(2) 
!      el(6) = ina*gb(3)-inb*ga(3) 
!      denom = ina*imb-inb*ima
!      ifit(ima,imb,ina,inb)=100
!      if (denom.ne.0) then
!       fel = float(el)/float(denom)
!       fit = sum(abs(float(int(fel))-fel))
! ! here is where we only keep the integer combinations
!       if (fit.eq.0.0) then
!         gone(1:3) = int(fel(1:3))
!         gtwo(1:3) = int(fel(4:6))
! ! keep the sum of the squares of the lengths 
!        ifit(ima,imb,ina,inb)=sum(gone**2)+sum(gtwo**2) 
!       end if
!      end if
!     end do
!    end do
!   end do
!  end do
!  minsum = 50

! ! look for the minimum of ifit with the smallest and most
! ! positive coefficients; store them in ml
! ! [minloc does not work here because there may be multiple minima]
!  do ima=-inm,inm
!   do imb=-inm,inm
!    do ina=-inm,inm
!     do inb=-inm,inm
!      if (ifit(ima,imb,ina,inb).le.minsum) then
!       minsum = ifit(ima,imb,ina,inb)
!       ml(1) = ima
!       ml(2) = imb
!       ml(3) = ina
!       ml(4) = inb
!      end if
!     end do
!    end do
!   end do
!  end do
!  deallocate(ifit)

! ! transform ga and gb into g1 and g2 
!  gone = (ml(2)*ga-ml(1)*gb)/(ml(3)*ml(2)-ml(4)*ml(1))
!  gtwo = (ml(3)*gb-ml(4)*ga)/(ml(3)*ml(2)-ml(4)*ml(1))

! ! next rank these two vectors so that their cross product is along +k
!  call CalcCross(cell,float(gone),float(gtwo),gsave,'r','r',0)
!  fit = CalcDot(cell,gsave,float(k),'r')
!  if (fit.lt.0.0) then
!   igsave = gone
!   gone = gtwo
!   gtwo = igsave
!  end if

! ! finally, if isym.ne.0 make sure that the selection of the 
! ! basis vectors for the 3-fold and 6-fold 2D point groups is
! ! correctly done.
! !
! ! For isym < 7:   90 degrees between vectors (should not be a problem)
! ! For isym = 7:  120 degrees between vectors (should be ok too)
! ! distinguish between 3 coming from a cubic group
! ! vs. the same symmetries originating from a hexagonal setting
!  if ((isym.eq.7).and.(cell%gamma.eq.120.0)) then
!    isym=11
!    fit = CalcAngle(cell,float(gone),float(gtwo),'r')*180.0/cPi 
!    if (abs(fit-120.0).lt.1.0) then
!      gtwo=gone+gtwo
!    end if
!  end if

! ! For isym = 8:  here we should distinguish between the settings 3m1 and 31m !!!
! !                The angle should always be 120 degrees, so we must check that
! !                this is the case for the selected gone and gtwo.
!  if ((isym.eq.8).and.(cell%gamma.eq.120.0)) then
!    isym=12
!    fit = CalcAngle(cell,float(gone),float(gtwo),'r')*180.0/cPi 
!    if (abs(fit-120.0).lt.1.0) then
!      gtwo=gone+gtwo
!    end if
!  end if
! !
!  if (isym.eq.8) then
!    fit = CalcAngle(cell,float(gone),float(gtwo),'r')*180.0/cPi 
!    if (abs(fit-120.0).gt.1.0) then
!      gtwo=gtwo-gone
!    end if

! ! we assume it is the 31m setting;  if the order of gone is 6, then that is not true
!    call CalcFamily(cell,gone,num,'r',itmp)
!    call GetOrder(float(k),il,num,jcnt,itmp)
!    if (jcnt.eq.6) then  ! it is the 3m1 setting
!      isym = 13
!    end if
!  end if

! ! it could the 3m1 setting for the 3m hexagonal case
!  if (isym.eq.12) then
! ! we assume it is the 31m setting;  if the order of gone is 6, then that is not true
!    call CalcFamily(cell,gone,num,'r',itmp)
!    call GetOrder(float(k),il,num,jcnt,itmp)
!    if (jcnt.eq.6) then  ! it is the 3m1 setting
!      isym = 14
!    end if
!  end if

! ! For isym = 9 or 10:   60 degrees between vectors (may not be the case)
!  if ((isym.eq.9).or.(isym.eq.10)) then
!    fit = CalcAngle(cell,float(gone),float(gtwo),'r')*180.0/cPi 
!    if (abs(fit-120.0).lt.1.0) then
!      gtwo=gone+gtwo
!    end if
!  end if

! end subroutine ShortestG_



! !--------------------------------------------------------------------------
! !
! ! SUBROUTINE: GetOR
! !
! !> @author Marc De Graef, Carnegie Mellon University
! !
! !> @brief ask for orientation relation between two crystals
! !
! !> @details ask for orientation relation between two crystals in terms of parallel planes
! !> and parallel directions; 
! !
! !> @param orel output variable of type orientation
! !
! !> @todo Is this routine really necessary ? It is not called very often.
! !
! !> @date 10/13/98 MDG 1.0 original
! !> @date    5/19/01 MDG 2.0 f90 version
! !> @date   11/27/01 MDG 2.1 added kind support
! !> @date   03/19/13 MDG 3.0 interface support
! !> @date   01/10/14 MDG 4.0 checked for changes to unitcell type
! !> @date   06/05/14 MDG 4.1 added stdout 
! !> @date   03/29/18 MDG 4.2 removed stdout 
! !--------------------------------------------------------------------------
! recursive subroutine GetOR(orel)
! !DEC$ ATTRIBUTES DLLEXPORT :: GetOR

! use io

! IMPLICIT NONE

! type(orientation),INTENT(OUT)                   :: orel                 !< orientation relation type

! real(kind=sgl)                                  :: c1,c2                !< auxiliary variables
! integer(kind=irg)                               :: io_int(6)            !< used for IO

 
!  c1 = 1.0_sgl
!  c2 = 1.0_sgl
!  do while ((c1.ne.0.0_sgl).or.(c2.ne.0.0_sgl))
!   call Message('Enter orientation relation in following form:', frm = "(A)")
!   call Message('planes:     h_A,k_A,l_A,h_B,k_B,l_B ', frm = "(A)")
!   call Message('directions: u_A,v_A,w_A,u_B,v_B,w_B ', frm = "(A)")
!   call ReadValue('Plane normals :', io_int, 6) 
!   orel%gA(1:3) = float(io_int(1:3))
!   orel%gB(1:3) = float(io_int(4:6))
!   call ReadValue('Directions    :', io_int, 6) 
!   orel%tA(1:3) = float(io_int(1:3))
!   orel%tB(1:3) = float(io_int(4:6))

! ! check for orthonormality using zone equation
!   c1=sum(orel%tA*orel%gA)
!   if (c1.ne.0.0_sgl) then
!    call Message('Plane does not contain direction (crystal A)', frm ="(A)")
!   end if
!   c2=sum(orel%tB*orel%gB)
!   if (c2.ne.0.0_sgl) then
!    call Message('Plane does not contain direction (crystal B)', frm = "(A)")
!   end if
!  end do

! end subroutine GetOR


! !--------------------------------------------------------------------------
! !
! ! SUBROUTINE: ComputeOR
! !
! !> @author Marc De Graef, Carnegie Mellon University
! !
! !> @brief compute the orientation relation transformation matrix
! !
! !> @param orel output variable of type orientation
! !> @param cellA unit cell A pointer
! !> @param cellB unit cell B pointer
! !> @param direction 'AB' for A new, B old; 'BA' for B new, A old
! !
! !> @date   12/20/13 MDG 1.0 first version, used for EMoverlap and EMorient
! !> @date   01/10/14 MDG 4.0 checked for changes to unitcell type
! !> @date   06/05/14 MDG 4.1 modification for cell pointers
! !> @date   09/21/15 SS  4.2 correction in final step (buggy lines commented out and new lines marked)
! !--------------------------------------------------------------------------
! recursive function ComputeOR(orel, cellA, cellB, direction) result(TT)
! !DEC$ ATTRIBUTES DLLEXPORT :: ComputeOR

! use math
! use io
! use error

! IMPLICIT NONE

! type(orientation),INTENT(INOUT) :: orel         !< orientation relation type
! !f2py intent(in,out) ::  orel         !< orientation relation type
! type(unitcell)                  :: cellA, cellB 
! character(2),INTENT(IN)         :: direction  !< direction of transformation (AB or BA)
! real(kind=sgl)                  :: TT(3,3)

! real(kind=sgl)                  :: r(3), p(3), Ep(3,3), E(3,3)
! real(kind=dbl)                  :: dE(3,3)
! integer(kind=irg)               :: i


! ! compute E matrix  [page 74]
!  call TransSpace(cellA, orel % gA,r,'r','d')
!  call NormVec(cellA, r,'d')
!  call NormVec(cellA, orel % tA,'d')
!  call CalcCross(cellA, orel % tA,r,p,'d','d',0)
!  call NormVec(cellA, p,'d')
!  E(1,1:3)=r(1:3)
!  E(2,1:3)=p(1:3)
!  E(3,1:3)=orel % tA(1:3)
!  !if (direction.eq.'AB') then
!  !  call mInvert(dble(E),dE,.FALSE.)
!  !  E = sngl(dE)
!  !end if

! ! compute E-prime matrix 
!  call TransSpace(cellB, orel % gB,r,'r','d')
!  call NormVec(cellB, r,'d')
!  call NormVec(cellB, orel % tB,'d')
!  call CalcCross(cellB, orel % tB,r,p,'d','d',0)
!  call NormVec(cellB, p,'d')
!  Ep(1,1:3)=r(1:3)
!  Ep(2,1:3)=p(1:3)
!  Ep(3,1:3)=orel % tB(1:3)
!  !if (direction.eq.'BA') then
!  !  call mInvert(dble(Ep),dE,.FALSE.)
!  !  Ep = sngl(dE)
!  !end if

! ! and multiply E^(-1)Ep or Ep^(-1)E matrices to get transformation matrix M
!  if (direction.eq.'BA') then
! !=============NEW==================
!    call mInvert(dble(Ep),dE,.FALSE.)
!    Ep = sngl(dE) 
! !==================================
!    TT = matmul(Ep,E)
!  else if (direction .eq. 'AB') then
! !=============NEW==================
!    call mInvert(dble(E),dE,.FALSE.)
!    E = sngl(dE)
! !==================================   
!    TT = matmul(E,Ep)
!  else
!    call FatalError('ComputeOR','Unknown direction specified')
!  end if

! end function ComputeOR

! !--------------------------------------------------------------------------
! ! 
! ! FUNCTION:CalcsgHOLZ
! !
! !> @author Marc De Graef, Carnegie Mellon University
! !
! !> @brief compute the excitation error including HOLZ and Laue Center information
! !
! !> @details  see chapter 3
! !
! !> @param cell unit cell pointer
! !> @param HOLZdata HOLZ data structure
! !> @param gg input g vector
! !> @param kt tangential components of wave vector
! !> @param lambda electron wavelength
! ! 
! !> @date   10/16/13 MDG 1.0 new version, includes HOLZ stuff
! !> @date   01/10/14 MDG 4.0 checked for changes to unitcell type
! !> @date   06/05/14 MDG 4.1 added unit cell pointer argument and HOLZdata argument
! !--------------------------------------------------------------------------
! recursive function CalcsgHOLZ(cell,HOLZdata,gg,kt,lambda) result(exer)
! !DEC$ ATTRIBUTES DLLEXPORT :: CalcsgHOLZ

! IMPLICIT NONE

! type(unitcell)                          :: cell
! type(HOLZentries),INTENT(INOUT)         :: HOLZdata
! !f2py intent(in,out) ::  HOLZdata
! real(kind=sgl),INTENT(IN)               :: gg(3), kt(3), lambda

! real(kind=sgl)                          :: exer, g1len, g2len
! real(kind=sgl)                          :: ll(3), lpg(3), glen, gplen, LC1, LC2, LC3, sgdenom


! glen = CalcLength(cell,gg,'r')
! g1len = CalcLength(cell,HOLZdata%g1,'r')
! g2len = CalcLength(cell,HOLZdata%g2,'r')
! if (glen.ne.0.0) then
!   LC1 = CalcDot(cell,kt,HOLZdata%g1,'r')/g1len
!   LC2 = CalcDot(cell,kt,HOLZdata%g2,'r')/g2len
!   ll = LC1*HOLZdata%g1 + LC2*HOLZdata%g2
!   lpg = ll + gg
!   gplen = CalcLength(cell,lpg,'r')
!   LC3 = sqrt(1.0-lambda**2*CalcLength(cell,ll,'r')**2)
!   if (gplen.eq.0.0) then
!     exer = -lambda*CalcDot(cell,gg,2.0*ll+gg,'r')/(2.0*LC3*CalcDot(cell,HOLZdata%g3,HOLZdata%FNr,'r'))
!   else
!     sgdenom = 2.0*CalcDot(cell,LC3*HOLZdata%g3-lambda*lpg,HOLZdata%FNr,'r')
!     exer = (CalcDot(cell,lpg,2.0*LC3*HOLZdata%g3-lambda*gg,'r')-lambda*CalcDot(cell,gg,ll,'r'))/sgdenom
!   end if
! else
!   exer = 10000.0
! end if

! end function CalcsgHOLZ




! !--------------------------------------------------------------------------
! ! 
! ! SUBROUTINE:GetHOLZGeometry
! !
! !> @author Marc De Graef, Carnegie Mellon University
! !
! !> @brief initialize HOLZ geometrical data for a given zone axis
! !
! !> @details  see chapter 3
! !
! !> @param cell unit cell pointer
! !> @param HOLZdata HOLZ data structure
! !> @param g1 first ZOLZ vector
! !> @param g2 second ZOLZ vector
! !> @param uvw zone axis
! !> @param FN foil normal
! ! 
! !> @date 10/17/13 MDG 1.0 original
! !> @date   01/10/14 MDG 4.0 checked for changes to unitcell type
! !> @date   06/05/14 MDG 4.1 added unit cell pointer argument and HOLZdata argument
! !--------------------------------------------------------------------------
! recursive subroutine GetHOLZGeometry(cell,HOLZdata,g1,g2,uvw,fn)
! !DEC$ ATTRIBUTES DLLEXPORT :: GetHOLZGeometry

! use io
! use error

! IMPLICIT NONE

! type(unitcell)                          :: cell
! type(HOLZentries),INTENT(INOUT)         :: HOLZdata
! !f2py intent(in,out) ::  HOLZdata
! integer(kind=irg),INTENT(IN)            :: uvw(3), fn(3)
! real(kind=sgl),INTENT(IN)               :: g1(3), g2(3)

! real(kind=sgl)                          :: gmin,gam11,gam12,gam22, phi, glen, g3(3), c(3), gx(3), gy(3), gshort(3)
! integer(kind=irg),parameter             :: inm = 8
! integer(kind=irg)                       :: ih,ik,il,NN, oi_int(1)

! ! set some basic values
!     HOLZdata%g1 = g1
!     HOLZdata%g2 = g2
!     HOLZdata%uvw = uvw
!     HOLZdata%FN = fn
    
! ! distance between consecutive HOLZ layers in nm-1
!     HOLZdata%H = 1.0/CalcLength(cell,float(uvw),'d')

! ! determine g3 basis vector
!     call CalcCross(cell,HOLZdata%g1,HOLZdata%g2,g3,'r','r',1)
!     call NormVec(cell,g3,'r')
!     HOLZdata%g3 = HOLZdata%H * g3

! ! compute components of FN with respect to ga, gb, g3
!     call TransSpace(cell,float(HOLZdata%FN),HOLZdata%FNr,'d','r')
!     call NormVec(cell,HOLZdata%FNr,'r')
!     HOLZdata%FNg = (/ CalcDot(cell,HOLZdata%FNr,HOLZdata%g1,'r'), CalcDot(cell,HOLZdata%FNr,HOLZdata%g2,'r'), &
!                         CalcDot(cell,HOLZdata%FNr,g3,'r') /)

! ! look for the shortest reflection satisfying hu+kv+lw = 1
! ! This could be replaced by code from Jackson's paper (1987),
! ! but it does essentially the same thing.
!  gmin = 100.0
!  NN=1
!  do while((gmin.eq.100.0).and.(NN.lt.4))
!   do ih=-inm,inm
!    do ik=-inm,inm
!     do il=-inm,inm
! ! does this reflection lie in the plane NN ?
!      if ((ih*uvw(1)+ik*uvw(2)+il*uvw(3)).eq.NN) then
!       glen = CalcLength(cell,float((/ih,ik,il/)),'r')
!       if (glen.lt.gmin) then
!        gmin = glen
!        gshort = float( (/ ih,ik,il /) )
!       end if
!      end if
!     end do
!    end do
!   end do
!   oi_int(1) = NN
!   call WriteValue(' Could not find any reflections with hu+kv+lw = ', oi_int, 1, frm = "(I2)")
!   NN = NN+1
!  end do
!  if (gmin.eq.100.0) then ! for some reason there is no reflection with N<=3 ...
!   call FatalError('ShortestGFOLZ: ',' could not find any reflections with hu+kv+lw<=3 ...')
!  end if
!  HOLZdata%gshort = gshort

! ! projected components of G
!  gam11 = CalcDot(cell,g1,g1,'r')
!  gam12 = CalcDot(cell,g1,g2,'r')
!  gam22 = CalcDot(cell,g2,g2,'r')
!  gmin = 1.0/(gam11*gam22-gam12**2)
!  HOLZdata%gp(1) = (CalcDot(cell,gshort,g1,'r')*gam22-CalcDot(cell,gshort,g2,'r')*gam12)*gmin
!  HOLZdata%gp(2) = (CalcDot(cell,gshort,g2,'r')*gam11-CalcDot(cell,gshort,g1,'r')*gam12)*gmin

! ! coordinate transformation matrix for g1 along x (our standard orientation for all programs)
!  phi = CalcAngle(cell,g1,g2,'r')
!  glen = CalcLength(cell,g2,'r')
!  HOLZdata%gtoc(1,1) = CalcLength(cell,g1,'r')
!  HOLZdata%gtoc(1,2) = glen * cos(phi)
!  HOLZdata%gtoc(2,1) = 0.0
!  HOLZdata%gtoc(2,2) = glen * sin(phi)

! ! first normalize the zone axis in cartesian components; this is the z-axis
!   call TransSpace(cell,float(uvw),c,'d','c')
!   call NormVec(cell,c,'c')

! ! then make ga the x-axis
!   call TransSpace(cell,g1,gx,'r','c')
!   call NormVec(cell,gx,'c')
!   HOLZdata%gx = gx

! ! compute the cross product between k and gx; this is the y-axis
!   call CalcCross(cell,c,gx,gy,'c','c',0)
!   HOLZdata%gy = gy


! end subroutine GetHOLZGeometry

! !--------------------------------------------------------------------------
! ! 
! ! FUNCTION:GetHOLZcoordinates
! !
! !> @author Marc De Graef, Carnegie Mellon University
! !
! !> @brief find the projected coordinates of an arbitrary HOLZ g-vector
! !
! !> @details  see chapter 3
! !
! !> @param gg input g vector
! !> @param kt tangential wave vector component
! !> @param lambda electron wavelength
! ! 
! !> @date 1/29/02  MDG 1.0 original
! !> @date 04/08/13 MDG 2.0 rewrite
! !> @date 10/16/13 MDG 3.0 incorporation into LACBED code
! !> @date 01/10/14 MDG 4.0 checked for changes to unitcell type
! !> @date   06/05/14 MDG 4.1 added unit cell pointer argument and HOLZdata argument
! !--------------------------------------------------------------------------
! recursive function GetHOLZcoordinates(cell,HOLZdata,gg,kt,lambda) result(pxy)
! !DEC$ ATTRIBUTES DLLEXPORT :: GetHOLZcoordinates

! IMPLICIT NONE

! type(unitcell)                          :: cell
! type(HOLZentries),INTENT(INOUT)         :: HOLZdata
! !f2py intent(in,out) ::  HOLZdata
! real(kind=sgl),INTENT(IN)               :: gg(3), kt(3), lambda

! real(kind=sgl)                          :: pxy(2), h1, h2, g11, g12, g22, z
! real(kind=sgl)                          :: exer, correction, gxy(2), nx, ny, hh(3)
! integer(kind=irg)                       :: N

! ! get the Laue zone number
!         N = abs( HOLZdata%uvw(1)*gg(1) + HOLZdata%uvw(2)*gg(2) + HOLZdata%uvw(3)*gg(3) )

! ! get components of gg w.r.t. g1 and g2
!         hh = gg - N * HOLZdata%gshort 
!         h1 = CalcDot(cell,hh,HOLZdata%g1,'c')
!         h2 = CalcDot(cell,hh,HOLZdata%g2,'c')
!         g11 = CalcDot(cell,HOLZdata%g1,HOLZdata%g1,'c')
!         g12 = CalcDot(cell,HOLZdata%g1,HOLZdata%g2,'c')
!         g22 = CalcDot(cell,HOLZdata%g2,HOLZdata%g2,'c')
!         z = 1.0/(g12**2-g11*g22)
!         nx = (g12*h2-g22*h1)*z
!         ny = (g12*h1-g11*h2)*z

! ! compute excitation error, including Laue center, foil normal, and HOLZ reflection.
!         exer = CalcsgHOLZ(cell,HOLZdata,gg,kt,lambda)

! ! next, determine the drawing coordinates, first in terms of g1 and g2
!         correction = 1.0/(1.0-lambda*HOLZdata%H*(float(N)+exer*HOLZdata%FNg(3)))
!         gxy = (/ (nx+N*HOLZdata%gp(1)+exer*HOLZdata%FNg(1)), (ny+N*HOLZdata%gp(2)+exer*HOLZdata%FNg(2))  /) * correction

! ! convert to Cartesian drawing coordinates
!         pxy = matmul(HOLZdata%gtoc,gxy)

! end function GetHOLZcoordinates
 
 
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
