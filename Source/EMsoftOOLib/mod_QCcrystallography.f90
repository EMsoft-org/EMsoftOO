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
  procedure, pass(self) :: GetQCAsymPos_
  procedure, pass(self) :: SaveQCDataHDF_
  procedure, pass(self) :: setfname_
  procedure, pass(self) :: getfname_
  procedure, pass(self) :: displayPeriodicTable

  generic, public :: get6Dindex => get6Dindex_
  generic, public :: invert6Dindex => invert6Dindex_
  generic, public :: GetQCLatParm => GetQCLatParm_
  generic, public :: GetQCAsymPos => GetQCAsymPos_
  generic, public :: SaveQCDataHDF => SaveQCDataHDF_
  generic, public :: setfname => setfname_
  generic, public :: getfname => getfname_

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

! allocate( QCcell%inverseIndex(nLUT, 6) )

end function QCcell_constructor

!--------------------------------------------------------------------------
subroutine QCcell_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 01/30/22
!!
!! destructor for the QCcell_T Class
 
IMPLICIT NONE

type(QCcell_T), INTENT(INOUT)  :: self 

call reportDestructor('QCcell_T')

end subroutine QCcell_destructor

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
recursive subroutine GetQCLatParm_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: GetQCLatParm_
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

 call Message%printMessage(' Using the higher dimensional cut-and-project approach:', frm = "(A)")
 call Message%printMessage(' -------------------------------', frm = "(A)")
 
  ! get the lattice parameters
 call Message%printMessage('Enter lattice parameters', frm = "(//A)")

  ! in-plane lattice parameters
 call Message%ReadValue('    a_i | i = {1,2,3,4,5,6} (hyper-cube) [nm] = ', io_real, 1)
 self%QClatparm = io_real(1)
  
end subroutine GetQCLatParm_

!--------------------------------------------------------------------------
recursive subroutine GetQCAsymPos_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: GetQCAsymPos_
!! author: MDG/SS 
!! version: 1.0 
!! date: 01/30/22
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
  call Message%printMessage(' ->  Fractional coordinates, site occupation, Bpar [nm^2],'&
  'Bperp [nm^2], radial atomic size (fraction of a_i) : ',&
  frm = "(A,' ')",advance="no")
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
recursive subroutine SaveQCDataHDF_(self, QCSG, EMsoft)
!DEC$ ATTRIBUTES DLLEXPORT :: SaveQCDataHDF_
!! author: MDG/SS 
!! version: 1.0 
!! date: 01/31/22
!!
!! save 3D QCcrystal structure data to an HDF file

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
real(kind=dbl)                          :: cellparams(3)
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

write (*,*) '-----SC_CrystalData----'
call HDF%stackdump()


groupname = SC_CrystalData
hdferr = HDF%createGroup(groupname)
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%createGroup:'//trim(groupname), hdferr)

write (*,*) '-----SC_ProgramName----'
call HDF%stackdump()

dataset = SC_ProgramName
line2(1) = trim(progname)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetStringArray:'//trim(dataset), hdferr)

write (*,*) '-----SC_CreationDate----'
call HDF%stackdump()

dataset = SC_CreationDate
line2(1) = dstr
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetStringArray:'//trim(dataset), hdferr)

write (*,*) '-----SC_CreationTime----'
call HDF%stackdump()

dataset = SC_CreationTime
line2(1) = tstr
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetStringArray:'//trim(dataset), hdferr)

write (*,*) '-----SC_Creator----'
call HDF%stackdump()

dataset = SC_Creator
line2(1) = trim(EMsoft%getConfigParameter('UserName'))
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetStringArray:'//trim(dataset), hdferr)

write (*,*) '-----SC_LatticeParameters----'
call HDF%stackdump()

dataset = SC_LatticeParameters
cellparams = (/ self%QClatparm, self%alphaij, self%alphastarij /)
hdferr = HDF%writeDatasetDoubleArray(dataset, cellparams, 3)
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetDoubleArray:'//trim(dataset), hdferr)

write (*,*) '-----SC_SpaceGroupNumber----'
call HDF%stackdump()

dataset = SC_SpaceGroupNumber
hdferr = HDF%writeDatasetInteger(dataset, QCSG%getSGnum())
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetInteger:'//trim(dataset), hdferr)

! the axial symmetry will be used to differentiate between QCs and regular crystals
! for the 2D QCs, the axial symmetry will be equal to the axial symmetry along the 
! periodic direction. for icosahedral QCs, this will be set to 0, which will be 
! reserved for this case. this is done for minimum code change for dictionary indexing
! of iQCs.
write (*,*) '-----naxial----'
call HDF%stackdump()

dataset = SC_AxialSymmetry 
naxial = 0
hdferr = HDF%writeDatasetInteger(dataset, naxial)

write (*,*) '-----ATOM_ntype----'
call HDF%stackdump()

dataset = SC_Natomtypes
hdferr = HDF%writeDatasetInteger(dataset, self%ATOM_ntype)
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetInteger:'//trim(dataset), hdferr)

write (*,*) '-----atomtypes----'
call HDF%stackdump()

allocate(atomtypes(self%ATOM_ntype))
atomtypes(1:self%ATOM_ntype) = self%ATOM_type(1:self%ATOM_ntype)
dataset = SC_Atomtypes
hdferr = HDF%writeDatasetIntegerArray(dataset, atomtypes, self%ATOM_ntype)
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetIntegerArray:'//trim(dataset), hdferr)
deallocate(atomtypes)

write (*,*) '-----atompos----'
call HDF%stackdump()

allocate(atompos(self%ATOM_ntype,10))
atompos(1:self%ATOM_ntype,1:10) = self%ATOM_pos(1:self%ATOM_ntype,1:10)
dataset = SC_AtomData
hdferr = HDF%writeDatasetFloatArray(dataset, atompos, self%ATOM_ntype, 10)
if (hdferr.ne.0) call HDF%error_check('SaveDataHDF:HDF%writeDatasetFloatArray:'//trim(dataset), hdferr)
deallocate(atompos)

call HDF%stackdump()



call HDF%pop(.TRUE.)
call closeFortranHDFInterface()

end subroutine SaveQCDataHDF_


















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

class(QCCell_T),intent(inout)           :: self
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