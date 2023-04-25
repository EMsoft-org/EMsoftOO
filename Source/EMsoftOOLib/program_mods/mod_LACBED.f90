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

module mod_LACBED
    !! author: MDG 
    !! version: 1.0 
    !! date: 01/22/20
    !!
    !! class definition for the EMLACBED program
  
  use mod_kinds
  use mod_global
  
  IMPLICIT NONE 
  
  ! namelist for the EMLACBED program
  type, public :: LACBEDNameListType
    integer(kind=irg)       :: k(3)
    integer(kind=irg)       :: fn(3)
    integer(kind=irg)       :: maxHOLZ
    integer(kind=irg)       :: numthick
    integer(kind=irg)       :: npix
    integer(kind=irg)       :: nthreads
    real(kind=sgl)          :: voltage
    real(kind=sgl)          :: dmin
    real(kind=sgl)          :: convergence
    real(kind=sgl)          :: startthick
    real(kind=sgl)          :: thickinc
    real(kind=sgl)          :: minten
    character(fnlen)        :: xtalname
    character(fnlen)        :: outname
    character(fnlen)        :: BetheParametersFile
  end type LACBEDNameListType
  
  ! class definition
  type, public :: LACBED_T
  private 
    character(fnlen)       :: nmldeffile = 'EMLACBED.nml'
    type(LACBEDNameListType)  :: nml 
  
  contains
  private 
    procedure, pass(self) :: readNameList_
    procedure, pass(self) :: writeHDFNameList_
    procedure, pass(self) :: getNameList_
    procedure, pass(self) :: Calc2DFamily_
    procedure, pass(self) :: CheckPatternSymmetry_
    procedure, pass(self) :: LACBED_
  
    generic, public :: getNameList => getNameList_
    generic, public :: writeHDFNameList => writeHDFNameList_
    generic, public :: readNameList => readNameList_
    generic, public :: Calc2DFamily => Calc2DFamily_
    generic, public :: CheckPatternSymmetry => CheckPatternSymmetry_
    generic, public :: LACBED => LACBED_
  
  end type LACBED_T
  
  ! the constructor routine for this class 
  interface LACBED_T
    module procedure LACBED_constructor
  end interface LACBED_T
  
  contains
  
  !--------------------------------------------------------------------------
  type(LACBED_T) function LACBED_constructor( nmlfile ) result(LACBED)
  !DEC$ ATTRIBUTES DLLEXPORT :: LACBED_constructor
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! constructor for the LACBED_T Class; reads the name list 
   
  IMPLICIT NONE
  
  character(fnlen), OPTIONAL   :: nmlfile 
  
  call LACBED%readNameList(nmlfile)
  
  end function LACBED_constructor
  
  !--------------------------------------------------------------------------
  subroutine LACBED_destructor(self)
  !DEC$ ATTRIBUTES DLLEXPORT :: LACBED_destructor
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! destructor for the LACBED_T Class
   
  IMPLICIT NONE
  
  type(LACBED_T), INTENT(INOUT)  :: self 
  
  call reportDestructor('LACBED_T')
  
  end subroutine LACBED_destructor

  !--------------------------------------------------------------------------

subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! read the namelist from an nml file for the LACBED_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(LACBED_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg)       :: k(3)
integer(kind=irg)       :: fn(3)
integer(kind=irg)       :: maxHOLZ
integer(kind=irg)       :: numthick
integer(kind=irg)       :: npix
integer(kind=irg)       :: nthreads
real(kind=sgl)          :: voltage
real(kind=sgl)          :: dmin
real(kind=sgl)          :: convergence
real(kind=sgl)          :: startthick
real(kind=sgl)          :: thickinc
real(kind=sgl)          :: minten
character(fnlen)        :: xtalname
character(fnlen)        :: outname

namelist /LACBEDlist/ xtalname, voltage, k, fn, dmin, convergence, minten, &
                              nthreads, startthick, thickinc, numthick, outname, npix, maxHOLZ

k = (/ 0, 0, 1 /)               ! beam direction [direction indices]
fn = (/ 0, 0, 1 /)              ! foil normal [direction indices]
maxHOLZ = 2                     ! maximum HOLZ layer index to be used for the output file; note that his number
                                ! does not affect the actual computations; it only determines which reflection 
                                ! families will end up in the output file
numthick = 10                   ! number of increments
npix = 256                      ! output arrays will have size npix x npix
nthreads = 1                    ! number of computational threads
voltage = 200000.0              ! acceleration voltage [V]
dmin = 0.025                    ! smallest d-spacing to include in dynamical matrix [nm]
convergence = 25.0              ! beam convergence angle [mrad]
startthick = 10.0               ! starting thickness [nm]
thickinc = 10.0                 ! thickness increment
minten = 1.0E-6                 ! minimum intensity in diffraction disk to make it into the output file
xtalname = 'undefined'          ! initial value to check that the keyword is present in the nml file
outname = 'lacbedout.data'      ! output filename

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=LACBEDlist)
close(UNIT=dataunit,STATUS='keep')

! check for required entries
if (trim(xtalname).eq.'undefined') then
  call Message%printError('EMLACBED:',' structure file name is undefined in '//nmlfile)
end if
end if

self%nml%k = k
self%nml%fn = fn
self%nml%maxHOLZ = maxHOLZ
self%nml%numthick = numthick
self%nml%npix = npix
self%nml%nthreads = nthreads
self%nml%voltage = voltage
self%nml%dmin = dmin
self%nml%convergence = convergence
self%nml%startthick = startthick
self%nml%thickinc = thickinc
self%nml%minten = minten
self%nml%xtalname = xtalname
self%nml%outname = outname


end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! pass the namelist for the LACBED_T Class to the calling program

IMPLICIT NONE 

class(LACBED_T), INTENT(INOUT)          :: self
type(LACBEDNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(LACBED_T), INTENT(INOUT)          :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 4, n_real = 6
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( lacbednl => self%nml )

! create the group for this namelist
groupname = trim(HDFnames%get_NMLlist())
hdferr = HDF%createGroup(groupname)

! write all the single integers
io_int = (/ lacbednl%maxHOLZ, lacbednl%numthick, lacbednl%npix, lacbednl%nthreads /)
intlist(1) = 'maxHOLZ'
intlist(2) = 'numthick'
intlist(3) = 'npix'
intlist(4) = 'nthreads'

call HDF%writeNMLintegers(io_int, intlist, n_int)

! vectors
dataset = SC_k
hdferr = HDF%writeDatasetIntegerArray(dataset, lacbednl%k, 3)
if (hdferr.ne.0) call HDF%error_check('HDFwriteLACBEDNameList: unable to create k dataset', hdferr)

dataset = SC_fn
hdferr = HDF%writeDatasetIntegerArray(dataset, lacbednl%fn, 3)
if (hdferr.ne.0) call HDF%error_check('HDFwriteLACBEDNameList: unable to create fn dataset', hdferr)

! write all the single reals
io_real = (/ lacbednl%voltage, lacbednl%dmin, lacbednl%convergence, lacbednl%startthick, lacbednl%thickinc, lacbednl%minten/)
reallist(1) = 'voltage'
reallist(2) = 'dmin'
reallist(3) = 'convergence'
reallist(4) = 'startthick'
reallist(5) = 'thickinc'
reallist(6) = 'minten'
call HDF%writeNMLreals(io_real, reallist, n_real)

! write all the strings
dataset = SC_outname
line2(1) = lacbednl%outname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteLACBEDNameList: unable to create outname dataset', hdferr)

dataset = SC_xtalname
line2(1) = lacbednl%xtalname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteLACBEDNameList: unable to create xtalname dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

recursive subroutine Calc2DFamily_(self,cell,SG,ind,ksame,numksame,nunique,itmp)
!DEC$ ATTRIBUTES DLLEXPORT :: Calc2DFamily_
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! compute the indices of equivalent planes w.r.t. 2D symmetry and store them in the itmp array

use mod_global
use mod_io
use mod_crystallography
use mod_symmetry

IMPLICIT NONE


class(LACBED_T), INTENT(INOUT)          :: self
type(Cell_T)                            :: cell
type(SpaceGroup_T)                      :: SG
integer(kind=irg),INTENT(IN)            :: ind(3)                       !< input triplet
logical,INTENT(IN)                      :: ksame(*)                     !< list of symmetry operators
integer(kind=irg),INTENT(IN)            :: numksame                     !< number on the input list
integer(kind=irg),INTENT(OUT)           :: nunique                      !< number of equivalent entries generated
integer(kind=irg),INTENT(OUT)           :: itmp(48,3)                   !< array used for family computations etc

integer(kind=irg)                       :: m,i,j                        !< loop counters and such
real(kind=sgl)                          :: h,k,l,ih,ik,il,idiff !< auxiliary variables
logical                                 :: newpoint                     !< is this a new point ?
real,parameter                          :: eps=0.0001_sgl               !< comparison threshold

real(kind=dbl),allocatable              :: SGdirec(:,:,:)
integer(kind=irg)                       :: NUMpt

! first take the identity
 j=1
 itmp(j,1:3)=ind(1:3)
 h=float(ind(1))
 k=float(ind(2))
 l=float(ind(3))

 SGdirec = SG%getSpaceGroupPGdirecMatrices()
 NUMpt = SG%getSpaceGroupNUMpt()

! multiply with all point group elements that have the value .TRUE. in ksame
 do i=2,NUMpt 
  if (ksame(i)) then 
    ih=SGdirec(i,1,1)*h+SGdirec(i,1,2)*k+SGdirec(i,1,3)*l
    ik=SGdirec(i,2,1)*h+SGdirec(i,2,2)*k+SGdirec(i,2,3)*l
    il=SGdirec(i,3,1)*h+SGdirec(i,3,2)*k+SGdirec(i,3,3)*l

! is this a new point ?
   newpoint=.TRUE.
   do m=1,j+1
    idiff=(itmp(m,1)-ih)**2+(itmp(m,2)-ik)**2+(itmp(m,3)-il)**2
    if (idiff.lt.eps) newpoint=.FALSE.
   end do

   if (newpoint) then 
    j=j+1
    itmp(j,1)=nint(ih)
    itmp(j,2)=nint(ik)
    itmp(j,3)=nint(il)
   end if
  end if
 end do 
 nunique=j

end subroutine Calc2DFamily_

recursive subroutine CheckPatternSymmetry_(self,cell,SG,k,ga,isym,thetam)
!DEC$ ATTRIBUTES DLLEXPORT :: CheckPatternSymmetry_

use mod_crystallography
use mod_symmetry
use mod_io

IMPLICIT NONE
class(LACBED_T),INTENT(INOUT)           :: self
type(Cell_T)                            :: cell
type(SpaceGroup_T)                      :: SG
type(IO_T)                              :: Message  
integer(kind=irg),INTENT(IN)            :: k(3)         !< zone axis
integer(kind=irg),INTENT(IN)            :: ga(3)        !< g-vector
integer(kind=irg),INTENT(INOUT) :: isym         !< 2D point group number

real(kind=sgl),INTENT(OUT)              :: thetam       !< rotation angle (degrees, CCW)

integer(kind=irg)                       :: num
real(kind=sgl)                          :: io_real(1)
integer(kind=irg), allocatable          :: itmp(:,:)   !< array used for family computations etc

! no action is needed for the following 2D point groups: 1, 2, 2mm, 3, 4, 4mm, 6, 6mm
thetam = 0.0

! for the group m (isym=3), we need to determine the cardinality of the 
! family of ga; if equal to 1, then ga lies in the mirror plane and all
! is correct.  If the cardinality is 2, then we compute the angle between
! ga and ga' and set thetam to half of this angle.

if (isym.eq.3) then 
  call SG%CalcFamily(ga, num, 'r',itmp)
  if (num.ne.1) then
    thetam = 0.5 * cell%CalcAngle(float(ga),float(itmp(2,1:3)),'r') *180.0/cPi
  end if  
  io_real(1) = thetam
  call Message%WriteValue('  --> Pattern symmetry m correction; point group rotation angle [deg]',io_real, 1, "(F6.3/)")
end if

! for the groups 3m1 and 31m, we need to check which one we have

if ((isym.eq.8).or.(isym.eq.11)) then 
  call SG%CalcFamily(ga, num, 'r',itmp)
  if (num.eq.3) then
    isym = 11
  else
    isym = 8
  end if  
  call Message%printMessage('  --> Pattern symmetry verified to be '//PGTWD(isym), frm = "(A/)")
end if

end subroutine CheckPatternSymmetry_

!--------------------------------------------------------------------------
subroutine LACBED_(self, EMsoft, progname, nmldeffile)
!DEC$ ATTRIBUTES DLLEXPORT :: LACBED_
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! perform the computations

use mod_EMsoft
use mod_io
use mod_crystallography
use mod_symmetry
use mod_diffraction
use mod_math
use mod_HDFsupport
use mod_HDFnames
use HDF5
use mod_quaternions
use mod_rotations
use mod_initializers
use mod_gvectors
use mod_kvectors
use mod_OMPsupport
use mod_timing
use mod_HOLZ
use stringconstants

use ISO_C_BINDING

use, intrinsic :: iso_fortran_env

IMPLICIT NONE 

class(LACBED_T), INTENT(INOUT)          :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname
character(fnlen), INTENT(INOUT)         :: nmldeffile

type(LACBEDNameListType)                  :: nml
type(Cell_T)                              :: cell
type(reflisttype),pointer                 :: firstw, nexts, rl, rltmpa, rltmpb
type(gvectors_T)                          :: reflist
type(BetheParameterType)                  :: BetheParameters
type(gnode),save                          :: rlp
type(DynType),save                        :: Dyn
type(Diffraction_T)                       :: Diff, lDiff
type(HDF_T)                               :: HDF
type(HDFnames_T)                          :: HDFnames
type(kvectorlist),pointer                 :: klist, ktmp
type(kvectors_T)                          :: kvec
type(Timing_T)                            :: timer
type(SpaceGroup_T)                        :: SG
type(IO_T)                                :: Message
type(HOLZ_T)                              :: HOLZ
logical                                   :: verbose

integer(kind=irg)                         :: sLUT, i, ii, jj, ik, ithick, parity, hkl(6,23), Pmdims, h(6), gindex, io_int(6)
complex(kind=dbl)                         :: Ucg, qg, Ucg2, qg0
real(kind=dbl)                            :: Vmod, Vpmod, xig, xgp
real(kind=dbl)                            :: lgpar, lgperp, st, nfact
integer(HSIZE_T)                          :: dims4(4), cnt4(4), offset4(4)

real(kind=sgl)                            :: FN(3), kk(3), dmin, kp(3), ku(3), io_real(3), qu(4), gx(3), cen, delta, gac(3), sc,&
                                             ktmax, rad, theta, thetacr, voltage, ma, mi, kkk(3), bragg, galen, minten, thetam
real(kind=dbl)                            :: epar(6,3), mLambda
integer(kind=irg)                         :: nref, nns, nnw, counter, numk, ierr, ga(3), gb(3), imax, imin, ir, istat, isym
real(kind=sgl),allocatable                :: klistarray(:,:), image(:,:), familytwotheta(:), diskoffset(:,:), &
                                             rfamilytwotheta(:), rdiskoffset(:,:)
integer(kind=irg),allocatable             :: kpix(:,:), hklarray(:,:), ranking(:), rranking(:), refctr(:)
integer(kind=irg),allocatable             :: familymult(:), familyhkl(:,:), whichHOLZ(:), gequiv(:,:), rfamilymult(:), &
                                             rfamilyhkl(:,:), rwhichHOLZ(:)
integer(kind=irg)                         :: itmp(48,3), famhkl(3), famnum, ghkl(3), numksame, nunique
logical,allocatable                       :: ksame(:)
real(kind=dbl)                            :: s(3)

real(kind=dbl),allocatable                :: SGdirec(:,:,:)
integer(kind=irg)                         :: NUMpt

complex(kind=dbl),allocatable             :: DynMat(:,:)
type(HOLZentries)                         :: HOLZdata

logical                                   :: f_exists, firstloop, insert=.TRUE.
real(kind=sgl),allocatable                :: intensity(:,:,:), thick(:), inten(:,:), slice(:,:,:,:)
integer(kind=irg)                         :: nthreads, TID, j, jmax, jmin, NUMTHREADS, it, TIFF_nx, TIFF_ny, dgn, icnt, ih, numir

character(fnlen)                          :: dataset, instring, outname, groupname, fname, TIFF_filename, tmpfile
character(fnlen)                          :: datagroupname
character(fnlen)                          :: mode, xtalname
integer(HSIZE_T)                          :: dims3(3), cnt3(3) 
integer(kind=irg)                         :: hdferr, d(3), refcnt
integer(kind=irg)                         :: tstart, tstop, clock_rate, nn, nt, numt, pgnum, ijmax, Tnref, ifamily
real(kind=sgl)                            :: exec_time, kc(3), pxy(2)
character(11)                             :: dstr
character(15)                             :: tstrb
character(15)                             :: tstre
character(2)                              :: str

! simplify the notation a little
associate( cbednl => self%nml )

nml = self%getNameList_()

call openFortranHDFInterface()
! set the HDF group names for this program
HDF = HDF_T()
HDFnames = HDFnames_T()

! initialize the timing routines
timer = Timing_T()
tstrb = timer%getTimeString()

call timer%Time_tick(1)

voltage         = cbednl%voltage
minten          = cbednl%minten
dmin            = cbednl%dmin
xtalname        = trim(cbednl%xtalname)
nthreads        = cbednl%nthreads
!===========================================================================================
! CRYSTALLOGRAPHY
!===========================================================================================
verbose = .TRUE.

call cell%setFileName(xtalname)

call Diff%setrlpmethod('WK')
call Diff%setV(dble(voltage))

call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, dmin, verbose, useHDF=HDF)

! call Diff%Printrlp()

! determine the point group number
j=0
do i=1,32
    if (SGPG(i).le.SG%getSpaceGroupNumber()) j=i
end do
isym = j

! use the new routine to get the whole pattern 2D symmetry group, since that
! is the one that determines the independent beam directions.
dgn = SG%GetPatternSymmetry(cbednl%k,j,.TRUE.)
pgnum = j
isym = WPPG(dgn) ! WPPG lists the whole pattern point group numbers vs. diffraction group numbers


! to do all this requires knowledge of the subset of 3D symmetry operators that keeps 
! the incident beam direction invariant, so we determine this subset first.  This 
! code comes from the CalcStar routine, but we only need a portion of it.
SGdirec = SG%getSpaceGroupPGdirecMatrices()
NUMpt = SG%getSpaceGroupNUMpt()

allocate(ksame(NUMpt))
ksame = .FALSE.
numksame = 1
ksame(1) = .TRUE.

! get all the symmetry operation IDs that leave the zone axis invariant (skip the identity operation)
do i=2,NUMpt
  s = matmul(SGdirec(i,1:3,1:3),dble(cbednl%k)) 
  if (sum(abs(s-dble(cbednl%k))).lt.1.0D-10) then 
    ksame(i) = .TRUE.
    numksame = numksame+1
  end if
end do

! and output the number of 3D symmetry operators that leave k invariant; they form
! the subset that we are after... This should really coincide with the 2D whole pattern symmetry group,
! but with 3D operators instead of 2D operators; so now we have the actual symmetry matrices.
io_int(1) = numksame
call Message%WriteValue(' Number of 3D symmetry operators that leave k invariant : ',io_int, 1, "(I3)")
io_int(1) = PGTWDorder(WPPG(dgn))
call Message%WriteValue(' Order of Whole Pattern point group : ',io_int, 1, "(I3,/)")
allocate(gequiv(numksame,3))

!determine the shortest reciprocal lattice points for this zone
call cell%ShortestG(SG,cbednl%k,ga,gb,isym)
io_int(1:3)=ga(1:3)
io_int(4:6)=gb(1:3)
call Message%WriteValue(' Reciprocal lattice vectors : ', io_int, 6,"('(',3I3,') and (',3I3,')')")
call Message%printMessage('  (the first lattice vector is horizontal in the CBED pattern)')
call Message%printMessage(' ')
galen = cell%CalcLength(float(ga), 'r')

! get number of thicknesses for which to compute the CBED disk images
numt = cbednl%numthick
allocate(thick(numt),stat=istat)
thick = cbednl%startthick + cbednl%thickinc* (/ (float(i),i=0,numt-1) /)

! set up the incident wave vector list
nullify(klist)
thetacr   = cbednl%convergence / 1000.0
bragg     = Diff%CalcDiffAngle(cell,(/ga(1),ga(2),ga(3)/)) * 0.5
ktmax     = thetacr / bragg
ijmax     = cbednl%npix**2

kvec = kvectors_T()   ! initialize the wave vector list
call kvec%set_kinp(dble(cbednl%k))
call kvec%set_ktmax(dble(ktmax))
call kvec%set_mapmode('Conical')
call kvec%Calckvectors(cell,SG,Diff,dble(ga),cbednl%npix, cbednl%npix,ijmax,usehex=.FALSE.)

numk = kvec%get_numk() 
io_int(1) = numk
call Message%WriteValue(' # independent beam directions to be considered = ', io_int, 1, "(I8)")

! copy some of these values into regular arrays for easier access in parallel section
allocate(klistarray(4,numk), kpix(2,numk))
klistarray = 0.0
kpix       = 0

!point to the first beam direction
ktmp => kvec%get_ListHead()
do ii = 1,numk
  klistarray(1:3,ii)  =   ktmp%k(1:3)
  klistarray(4,ii)    =   ktmp%kn
  kpix(1:2,ii)        =   (/ktmp%i, ktmp%j/)
  ktmp                =>  ktmp%next
end do

! get rid of the linked list of k-vectors
call kvec%Delete_kvectorlist()

!===========================================================================================
!force dynamical matrix routine to read new Bethe parameters from file
! call Diff%SetBetheParameters(EMsoft, .FALSE., cbednl%BetheParametersFile)
call Diff%SetBetheParameters(EMsoft, .TRUE., cbednl%BetheParametersFile)
! nullify(mainreflist)
kkk = klistarray(1:3,1)
FN = kkk

! some parameters required for simulation
qg0 = cmplx(cPi,0.0) * Diff%getLUTqg((/0,0,0/))
qg0 = cmplx(0.D0,aimag(qg0))

!====================================
! in this computation, we are going to stick to the same set of reflections, and for each 
! incident beam direction we will update the excitation errors, apply the Bethe potentials,
! and then compute the intensities; the sequence of nref reflections is thus fixed for the 
! entire computation, and weak reflections are handled with the Bethe potentials when 
! appropriate.

reflist = gvectors_T()
call reflist%Initialize_ReflectionList(cell, SG, Diff, FN, kkk, cbednl%dmin, verbose)
nref = reflist%get_nref()

allocate(hklarray(3,nref))
hklarray = 0

nexts => reflist%Get_ListHead()
nexts => nexts%next
do ii = 1,nref
  hklarray(1:3,ii)  =   nexts%hkl
  nexts             =>  nexts%next 
end do

! get rid of the linked list of k-vectors
!call reflist%Delete_gvectorlist()

!======================================
! allocate intensity array
allocate(intensity(nref,numk,numt))
intensity  = 0.0

! !======================================
! ! for some unclear reason (Heisenbug...) this program produces incorrect results unless there is 
! ! a write statement immediately following the GetDynMat routine.  We'll write some stuff to a temp file
! ! and then delete the file afterwards until we can figure out what the underlying issue is [MDG, 11/28/18] 
! tmpfile= trim(EMsoft%getEMtmppathname())//'tmpoutput.txt'
! tmpfile= EMsoft%toNativePath(tmpfile)
! open(unit=600,file=trim(tmpfile),status='unknown',form='formatted')

!======================================
! initialize the OpenMP parallel processing environment
call OMP_setNThreads(cbednl%nthreads)

io_int(1) = cbednl%nthreads
call Message%WriteValue(' Attempting to set number of threads to ',io_int, 1, frm = "(I4)")

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(kk, firstloop, NUMTHREADS, TID, reflist, rl, nns, nnw) &
!$OMP& PRIVATE(DynMat, Tnref, inten, jj, ii, io_real, firstw, lDiff)

NUMTHREADS = OMP_GET_NUM_THREADS()
TID = OMP_GET_THREAD_NUM()
firstloop = .TRUE.
lDiff = Diff

! perform the main loop over the incident beam directions in parallel
!$OMP DO SCHEDULE(DYNAMIC)
do ik = 1,numk
! local wave vector
  kk  = klistarray(1:3,ik)

! we will have each thread do its own reflist initialization to allow for the lists to be used
! in combination with the Bethe potential approach; this needs to be done only once for each thread
  if (firstloop.eqv..TRUE.) then 
    reflist = gvectors_T()
    call reflist%Initialize_ReflectionList(cell, SG, lDiff, FN, kkk, nml%dmin, verbose=.False.)
    firstloop = .FALSE.
  end if

! go through the reflist and reset all strong and weak links to null pointers and .FALSE. identifiers;
! at the same time, compute the excitation errors for this incident beam direction kk
  rl => reflist%Get_ListHead()
  rl => rl%next
  do
    if (.not.associated(rl)) EXIT
    nullify(rl%nexts, rl%nextw)
    rl%strong = .FALSE.
    rl%weak = .FALSE.
    rl%sg = Diff%Calcsg(cell,float(rl%hkl(1:3)),kk,FN)
    rl => rl%next
  end do
  
! determine strong and weak reflections
  nullify(firstw)
  nns = 0
  nnw = 0
  call reflist%Apply_BethePotentials(lDiff, firstw, nns, nnw)

! generate the dynamical matrix
  allocate(DynMat(nns,nns))
  call reflist%GetDynMat(cell, lDiff, firstw, DynMat, nns, nnw)
  Tnref = reflist%get_nref()
!========================
! this line is here as a workaround for some unidentified bug in the GetDynMat routine ...
! without this write, the remainder of this loop occasionally produces NaN results
  ! if (mod(ik,1000) .eq. 0) write (600,"(4I10,F12.4)") ik, TID, nns, nnw, maxval(abs(DynMat))
!========================

! ! solve the dynamical eigenvalue equation for this beam direction  
  allocate(inten(numt,nns))
  call reflist%CalcCBEDint(cell,DynMat,klistarray(4,ik),nns,numt,thick,inten)

! and copy the intensities into the correct slot for the strong reflections only ...
  rl => reflist%Get_ListHead()
  jj = 1
  intensity(1,ik,1:numt) = inten(1:numt,jj)    ! first reflection is always BF 
  rl => rl%next%next
  do ii = 2, Tnref
    if (rl%strong.eqv..TRUE.) then
      jj = jj+1
      intensity(ii,ik,1:numt) = inten(1:numt,jj)
    end if
    rl => rl%next
  end do

! every now and then, print a progress update
  if(mod(ik,5000) .eq. 0) then
    io_real(1) = 100.D0 * dble(ik)/dble(numk)
    call Message%WriteValue('  completed ',io_real,1,'(F10.2,"% of beam directions")')
  end if
   deallocate(DynMat, inten)
end do
!$OMP END DO
!$OMP END PARALLEL

write (*,*) 'maximum intensity value = ', maxval(intensity)


! call timestamp(datestring=dstr, timestring=tstre)

! !===============================
! ! get rid of the bug workaround file
 close(unit=600,status='delete')
! !===============================



!===============================
!===============================
! next we need to compute a number of arrays needed by the visualization routine as well as the 
! EMCBEDpattern program.  The arrays have to do with reflector identifications...
!
! We need to go through the reflection list and identify for each 
! reflection the multiplicity with respect to the family generated by the
! 2D whole pattern point group [see manual for IDL visualization program
! for an explicit example].  In other words, we need to tag each reflection
! that will need to end up in the final output file for the LACBED program.
! For instance, for the [112] Copper zone axis orientation, (1 1 -1) and
! (-1 -1 1) do not belong to the same family (since they lie on the mirror
! plane that makes up the 2D Whole Pattern group m.  Therefore, these must 
! be tagged as separate families with appropriate multiplicities, otherwise
! the visualization program will only know of one reflection (1 1 -1) and
! will not be able to generate the other one (-1 -1 1) {which need not be 
! identical to begin with, due to the m symmetry].
! At the same time we need to determine how many independent families there 
! are, as well as their multiplicities.
 
! set the scale parameter for a default camera length of 1000 mm.
mLambda = Diff%getWaveLength()
sc = mLambda * 1000.0 * 300.0 / 25.4  ! the absolute value does not matter and is derived from legacy Postscript code
! The original code used 300 dpi (hence 300/25.4) which was convenient for Postscript output; in the current case, we
! do not actually use the true value, but in the IDL visualization program, we scale the user defined camera length by
! 1000.0, and use this ratio to scale the diskoffset coordinates.  So, there's no absolute length scale, only a relative scale.



! to be safe, we need to reset the family number for each of the reflections in the current list
! to zero, to reflect the fact that we haven't considered that reflection yet in terms of 2D symmetry.
  rltmpa => reflist%Get_ListHead()
  rltmpa => rltmpa%next
  do while (associated(rltmpa))
    rltmpa%famnum = 0
    if (.not.associated(rltmpa%next)) EXIT
    rltmpa => rltmpa%next
  end do

! next we go through the entire list again, this time with a double loop, and we
! determine for each hkl all the 2D equivalent ones and tag them as belonging to
! the same family in terms of 2D symmetry.  
  rltmpa => reflist%Get_ListHead()
  rltmpa%famnum = 1             ! reflection at origin is always its own family
  ifamily = 1
  rltmpa => rltmpa%next%next
    whileloop1: do while (associated(rltmpa))
    ! only look at points that haven't been tagged yet
    if (rltmpa%famnum.eq.0) then
        ghkl = rltmpa%hkl           ! get the reflection
        call self%Calc2DFamily(cell,SG,ghkl,ksame,numksame,nunique,itmp)
    ! we add this point to a new famnum value
        ifamily = ifamily + 1
        rltmpa%famnum = ifamily
        rltmpa%famhkl = ghkl
        famhkl = ghkl
        if (nunique.gt.1) then
    ! the order is larger than 1, so we need to go through the list and tag all the equivalent ones;
    ! we need only consider those that have the same original famhkl.
        do i=2,nunique
            rltmpb => rltmpa%next
    whileloop2: do while (associated(rltmpb))
    ! look for this reflection on the list
            if (sum(abs(itmp(i,1:3) - rltmpb%hkl(1:3))).eq.0) then  ! found it!
            rltmpb%famnum = ifamily
            rltmpb%famhkl = rltmpa%famhkl
            EXIT whileloop2
            end if
            if (.not.associated(rltmpb%next)) EXIT whileloop2
            rltmpb => rltmpb%next
            end do whileloop2
        end do
        end if
    end if
    if (.not.associated(rltmpa%next)) EXIT whileloop1
    rltmpa => rltmpa%next
    end do whileloop1

    ! ok, so there are ifamily families; next we need to store the corresponding
    ! hkl, and multiplicity, as well as the diffraction angle and the position of 
    ! the diffraction disk center for a standard camera length.
    allocate(familyhkl(3,ifamily), familymult(ifamily), familytwotheta(ifamily), diskoffset(2,ifamily), ranking(ifamily))

    ! redo the above loop, sort of, but now fill in the actual data
    ! we no longer need to keep the famnum entries in the linked list, so we
    ! can reset those to zero to keep track of points already visited.
    familyhkl(1:3,1) = (/ 0, 0, 0 /)
    familymult(1) = 1
    diskoffset(1:2,1) = (/ 0.0, 0.0 /)
    familytwotheta(1) = 0.0
    ranking(1) = 1
    rltmpa => reflist%Get_ListHead()
    rltmpa => rltmpa%next%next%next
    ifamily = 1
  
    ! for some of the 2D point groups, the standard orientation of the group according to ITC vol A
    ! may not be the orientation that we have here, so we need to determine by how much the 2D point
    ! group is rotated (CCW) with respect to the standard setting...
    call self%CheckPatternSymmetry(cell,SG,cbednl%k,ga,isym,thetam)

    ! initialize the HOLZ geometry type
    HOLZ = HOLZ_T()
    call HOLZ%GetHOLZGeometry(cell,HOLZdata,float(ga),float(gb),cbednl%k,cbednl%fn) 

! keep track of the reflection identifier in ranking
refcnt = 1
outerloop2: do while (associated(rltmpa))
  refcnt = refcnt+1
  if (rltmpa%famnum.ne.0) then 
    ifamily = ifamily+1
    famnum = rltmpa%famnum
    rltmpa%famnum = 0
    famhkl = rltmpa%famhkl
    familyhkl(1:3,ifamily) = famhkl(1:3)
    familytwotheta(ifamily) = Diff%CalcDiffAngle(cell,(/famhkl(1),famhkl(2),famhkl(3)/))*1000.0
    familymult(ifamily) = 1
    ranking(ifamily) = refcnt
! get the disk offset parameters
    pxy = sc * HOLZ%GetHOLZcoordinates(cell,HOLZdata,float(famhkl), (/ 0.0, 0.0, 0.0 /), sngl(mLambda))
    diskoffset(1:2,ifamily) = pxy
  
! and remove the equivalent reflections from the list
    rltmpb => rltmpa%next
    whileloop3: do while (associated(rltmpb))
      if (rltmpb%famnum.eq.famnum) then
         familymult(ifamily) = familymult(ifamily) + 1
         rltmpb%famnum = 0
      end if
      if (.not.associated(rltmpb%next)) EXIT whileloop3
      rltmpb => rltmpb%next
    end do whileloop3
   end if
   if (.not.associated(rltmpa%next)) EXIT outerloop2
   rltmpa => rltmpa%next
  end do outerloop2

  io_int(1) = ifamily
  call Message%WriteValue(' Maximum number of unique families in output = ', io_int, 1, "(I5)")

!===============================
!===============================
! before we write the intensities, we need to determine which reflection families 
! need to be written to the file; this is determined by the maxHOLZ parameter.  So,
! first we determine to which HOLZ layer each family belongs by using the zone 
! equation.   [check special case of hexagonal indices !!!!]
! Along the way, we count the ones up to order maxHOLZ.
! Also, to reduce the size of the output file a bit, we ignore those reflections that
! have a maximum intensity less than minten for the initial thickness value. On the
! other hand, if a reflection is due to double diffraction, then we include it, always.
  allocate(whichHOLZ(ifamily))
  icnt = 0
  do ir=1,ifamily
    whichHOLZ(ir) = iabs(cbednl%k(1)*familyhkl(1,ir)+cbednl%k(2)*familyhkl(2,ir)+cbednl%k(3)*familyhkl(3,ir))
    if (whichHOLZ(ir).le.cbednl%maxHOLZ) then 
      if ((maxval(intensity(ranking(ir),:,:)).ge.minten).or. &
          (Diff%getdbdiff((/familyhkl(1,ir),familyhkl(2,ir),familyhkl(3,ir)/)))) then
        icnt = icnt+1
      else  ! just change the HOLZ value to some large value to make sure it does not get written to the file
        whichHOLZ(ir) = 100
      end if
    end if
  end do    
  io_int(1) = icnt
  call Message%WriteValue(' Actual number of unique families in output = ', io_int, 1, "(I5)")

! ok, so the whichHOLZ array tells us which reflections need to be included, but they are likely not in the 
! same order as in the familyhkl array, so we have some reordering to do... this can be done by defining new 
! familyhkl etc arrays that will then be written to the file with the correct reflection ordering.
numir = icnt
allocate(rfamilyhkl(3,numir), rfamilymult(numir), rfamilytwotheta(numir), rdiskoffset(2,numir), rwhichHOLZ(numir))
allocate(rranking(numir))

! fill the new arrays with the correctly ranked data
! first the central reflection
rfamilyhkl(1:3,1) = familyhkl(1:3,1)
rfamilymult(1) = familymult(1)
rfamilytwotheta(1) = familytwotheta(1)
rdiskoffset(1:2,1) = diskoffset(1:2,1)
rwhichHOLZ(1) = whichHOLZ(1)
rranking(1) = ranking(1)
icnt = 1
! then the others
do ih = 0,cbednl%maxHOLZ
 do ir = ifamily,2,-1
  if (whichHOLZ(ir).eq.ih) then
    icnt = icnt +1 
    rfamilyhkl(1:3,icnt) = familyhkl(1:3,ir)
    rfamilymult(icnt) = familymult(ir)
    rfamilytwotheta(icnt) = familytwotheta(ir)
    rdiskoffset(1:2,icnt) = diskoffset(1:2,ir)
    rwhichHOLZ(icnt) = whichHOLZ(ir)
    rranking(icnt) = ranking(ir)
  end if
 end do 
end do 

! then report a final count of the number of reflections in the output for each HOLZ layer
allocate(refctr(0:cbednl%maxHOLZ))
refctr = 0
do ir = 1,numir
  refctr(rwhichHOLZ(ir)) =  refctr(rwhichHOLZ(ir)) + 1
end do

call Message%PrintMessage(' ')
do ir=0,cbednl%maxHOLZ
  io_int(1) = ir
  io_int(2) = refctr(ir)
  call Message%WriteValue('  Number of unique reflections in Laue Zone ',io_int, 2,"(I3,' : ', I4)") 
end do
call Message%PrintMessage(' ')


! !===============================
! ! HDF5 I/O
! ! write out the data to the file

! ! Create a new file using the default properties.
 outname = EMsoft%generateFilePath('EMdatapathname',trim(cbednl%outname))

 !=============================================
! create or update the HDF5 output file
!=============================================
 call HDFnames%set_ProgramData(SC_LACBED)
 call HDFnames%set_NMLlist(SC_LACBEDNameList)
 call HDFnames%set_NMLfilename(SC_LACBEDNML)

! Open an existing file or create a new file using the default properties.
 hdferr =  HDF%createFile(outname)

 ! write the EMheader to the file
 datagroupname = trim(HDFnames%get_ProgramData())
 call HDF%writeEMheader(EMsoft,dstr, tstrb, tstre, progname, datagroupname)

   ! create a namelist group to write all the namelist files into
 groupname = SC_NMLfiles
 hdferr = HDF%createGroup(groupname)

! read the text file and write the array to the file
 dataset = SC_LACBEDNameList
 hdferr = HDF%writeDatasetTextFile(dataset, nmldeffile)

 ! leave this group
 call HDF%pop()
  
 ! create a namelist group to write all the namelist files into
 hdferr = HDF%createGroup(HDFnames%get_NMLparameters())
 call self%writeHDFNameList(HDF, HDFnames)

 ! leave this group
 call HDF%pop()

! then the remainder of the data in a EMData group
! [this is the list of variables needed for the IDL visualization program CBEDDisplay.pro and the EMCBEDpattern.f90 program]
 groupname = SC_EMData
 hdferr = HDF%createGroup(groupname)

!========================================
! the following parameters are already part of either the header group or the namelist group ... 
! progname, EMsoftversion, npix, numt, xtalname, voltage, convergence, k, fn, dmin, maxHOLZ,
! startthick, thickinc

! to be written: icnt, numk, ga, galen, minten, pgnum, PGLaue(pgnum), dgn, PDG(dgn), BFPG(dgn), 
! WPPG(dgn), DFGN(dgn), DFSP(dgn), thetam, familyhkl, familymult, familytwotheta, 
! diskoffset, whichHOLZ, and all the diffraction disks, one per family... [using hyperslabs]
!========================================

! write integers 
dataset = 'numk'
hdferr = HDF%writeDatasetInteger(dataset, numk)

dataset = 'ifamily'
hdferr = HDF%writeDatasetInteger(dataset, ifamily)

dataset = 'ga'
hdferr = HDF%writeDatasetIntegerArray(dataset, ga, 3)

! integer arrays
dataset = 'diffgroup'
hdferr = HDF%writeDatasetIntegerArray(dataset, (/ pgnum, PGLaue(pgnum), dgn, PDG(dgn), BFPG(dgn), &
                                        WPPG(dgn), DFGN(dgn), DFSP(dgn) /), 8)

dataset = SC_hkl
hdferr = HDF%writeDatasetIntegerArray(dataset, hklarray, 3, nref)

dataset = SC_PixelLocation
hdferr = HDF%writeDatasetIntegerArray(dataset, kpix, 2, numk)

dataset = 'familymult'
hdferr = HDF%writeDatasetIntegerArray(dataset, rfamilymult, numir)

dataset = 'familyhkl'
hdferr = HDF%writeDatasetIntegerArray(dataset, rfamilyhkl, 3, numir)

! floats
dataset = 'galen'
hdferr = HDF%writeDatasetFloat(dataset, galen)

dataset = 'minten'
hdferr = HDF%writeDatasetFloat(dataset, minten)

dataset = 'thetam'
hdferr = HDF%writeDatasetFloat(dataset, thetam)

! float arrays
dataset = 'familytwotheta'
hdferr = HDF%writeDatasetFloatArray(dataset, rfamilytwotheta, numir)

dataset = 'diskoffset'
hdferr = HDF%writeDatasetFloatArray(dataset, rdiskoffset, 2, numir)

dataset = SC_klist  
hdferr = HDF%writeDatasetFloatArray(dataset, klistarray, 4, numk)

! before we write the intensities, we need to reorganize the array into a 4D array of diffraction disks 
! dataset = SC_Intensities
! hdferr = HDF%writeDatasetFloatArray(dataset, intensity, nref, numk, numt)


dataset = 'whichHOLZ'
hdferr = HDF%writeDatasetIntegerArray(dataset, rwhichHOLZ, numir)


! The last data to be written is the set of diffraction disks, one for each family... 
! We also need to define the hyperslab parameters so that we can write each of the 
! diffraction disks for the complete thickness range; the disks array is a 4D array, but 
! we write 3D slabs into it
! create the hyperslabs and write zeroes to them for now

dataset = 'disks'
  allocate(slice(2*cbednl%npix+1, 2*cbednl%npix+1, numt, 1))
  slice = 0.0
  dims4 = (/  2*cbednl%npix+1, 2*cbednl%npix+1, numt, numir /)
  cnt4 = (/ 2*cbednl%npix+1, 2*cbednl%npix+1, numt, 1 /)
  offset4 = (/ 0, 0, 0, 0 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, slice, dims4, offset4, cnt4)

! Finally, write the correct diffraction disks to the hyperslab array
! Note that in the original EMLACBED program, the first subscript into the diffraction disks was reversed,
! so we will do this here as well... (before we do the offset to positive subscripts!)
  kpix(1,:) = -kpix(1,:)
  kpix = kpix + cbednl%npix + 1
  do ir = 1,numir
    do ik=1,numk
      slice(kpix(1,ik),kpix(2,ik),1:numt,1) = intensity(rranking(ir),ik,1:numt)
    end do 
    offset4 = (/ 0, 0, 0, ir-1 /)
    hdferr = HDF%writeHyperslabFloatArray(dataset, slice, dims4, offset4, cnt4, insert)
  end do

! write the icnt number to the file
dataset = 'icnt'
  hdferr = HDF%writeDatasetInteger(dataset, numir)

! leave this group and close the file
call HDF%popall()

call Message%printMessage(' Output data stored in '//trim(outname))
call Message%printMessage(' ')

! timing data output
call timer%Time_tock(1)
io_real(1) = timer%getInterval(1)
call Message%WriteValue('Total run time [s] : ', io_real, 1)

end associate

end subroutine LACBED_


end module mod_LACBED
  
