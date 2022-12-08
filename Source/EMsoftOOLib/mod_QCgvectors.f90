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
module mod_QCgvectors
  !! author: MDG
  !! version: 1.0
  !! date: 02/07/22
  !!
  !! variables and types needed to determine lists of reciprocal lattice vectors
  !! for icosahedral and axial quasi-crystals.
  !!
  !! Due to some complicated module interdependencies the CalcDynMatQC routine is in
  !! this module rather than in mod_QCdiffraction, where it would logically belong.  We may need
  !! to figure out how to change that.

use mod_kinds
use mod_global

IMPLICIT NONE
private

! linked list of reflections
type, public :: QCreflisttype
  integer(kind=irg)             :: num, &        ! sequential number
                                   strongnum,&   ! sequential number for strong beams
                                   weaknum,&     ! sequential number for weak beams
                                   hkl5(5),&     ! 5 component Miller index
                                   hkl6(6)       ! 6 component Miller index
  real(kind=dbl)                :: sg, &         ! excitation error
                                   xg, &         ! extinction distance
                                   glen          ! length of reciprocal lattice vector
  logical                       :: strong, weak  ! is this a strong beam or not; both .FALSE. means 'do not consider'
  complex(kind=dbl)             :: Ucg           ! Fourier coefficient, copied from cell%LUT
  complex(kind=dbl)             :: qg            ! scaled Fourier coefficient, copied from cell%LUTqg
  type(QCreflisttype),pointer   :: next          ! connection to next entry in master linked list
  type(QCreflisttype),pointer   :: nexts         ! connection to next strong entry in linked list
  type(QCreflisttype),pointer   :: nextw         ! connection to next weak entry in linked list
end type QCreflisttype

! define the actual gvectors class; closely follows mod_gvectors.f90
type, public :: QCgvectors_T
private
  type(QCreflisttype), pointer  :: reflist
  type(QCreflisttype), pointer  :: rltail
  integer(kind=irg)             :: nref

contains
private
  procedure, pass(self) :: MakeQCRefList_
  procedure, pass(self) :: AddQCReflection_
  procedure, pass(self) :: Apply_BethePotentialsQC_
  procedure, pass(self) :: Delete_QCgvectorlist_
  procedure, pass(self) :: Initialize_QCReflectionList_
  ! procedure, pass(self) :: GetDynMat_
  procedure, pass(self) :: get_QCnref_
  procedure, pass(self) :: Get_QCListHead_

  final :: QCgvectors_destructor

  generic, public :: MakeQCRefList => MakeQCRefList_
  generic, public :: AddQCReflection => AddQCReflection_
  generic, public :: Apply_BethePotentialsQC => Apply_BethePotentialsQC_
  generic, public :: Delete_QCgvectorlist => Delete_QCgvectorlist_
  generic, public :: Initialize_QCReflectionList => Initialize_QCReflectionList_
  ! generic, public :: GetDynMat => GetDynMat_
  generic, public :: get_QCnref => get_QCnref_
  generic, public :: Get_QCListHead => Get_QCListHead_
end type QCgvectors_T

! the constructor routine for this class
interface QCgvectors_T
  module procedure QCgvectors_constructor
end interface QCgvectors_T

contains

!--------------------------------------------------------------------------
type(QCgvectors_T) function QCgvectors_constructor( ) result(GVec)
!DEC$ ATTRIBUTES DLLEXPORT :: QCgvectors_constructor
!! author: MDG
!! version: 1.0
!! date: 02/05/22
!!
!! constructor for the QCgvectors_constructor Class

IMPLICIT NONE

! initialize the reflist
nullify(GVec%reflist)
GVec%nref = 0
call GVec%MakeQCRefList()

end function QCgvectors_constructor

!--------------------------------------------------------------------------
subroutine QCgvectors_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: QCgvectors_destructor
!! author: MDG
!! version: 1.0
!! date: 02/05/22
!!
!! destructor for the gvectors_T Class

IMPLICIT NONE

type(QCgvectors_T), INTENT(INOUT)  :: self

call reportDestructor('gvectors_T')

end subroutine QCgvectors_destructor

!--------------------------------------------------------------------------
recursive subroutine MakeQCRefList_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: MakeQCRefList_
  !! author: MDG
  !! version: 1.0
  !! date: 02/05/22
  !!
  !! allocate and initialize the linked reflection list

use mod_io

IMPLICIT NONE

class(QCgvectors_T), INTENT(INOUT)  :: self

type(IO_T)                          :: Message
integer(kind=irg)                   :: istat

! create it if it does not already exist
if (.not.associated(self%reflist)) then
  allocate(self%reflist,stat=istat)
  if (istat.ne.0) call Message%printError('MakeQCRefList:',' unable to allocate pointer')
  self%rltail => self%reflist           ! tail points to new value
  nullify(self%rltail%next)             ! nullify next in new value
end if

end subroutine MakeQCRefList_

!--------------------------------------------------------------------------
recursive subroutine Delete_QCgvectorlist_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: Delete_QCgvectorlist_
  !! author: MDG
  !! version: 1.0
  !! date: 02/05/22
  !!
  !! delete the entire linked list

IMPLICIT NONE

class(QCgvectors_T), INTENT(INOUT)  :: self

type(QCreflisttype),pointer         :: rltmpa

! deallocate the entire linked list before returning, to prevent memory leaks
if (associated(self%reflist)) then
  self%rltail => self%reflist
  if (associated(self%rltail%next)) then
    rltmpa => self%rltail%next
    do
      if (associated(self%rltail)) deallocate(self%rltail)
      if (.not. associated(rltmpa)) EXIT
      self%rltail => rltmpa
      rltmpa => self%rltail%next
    end do
  end if
end if

nullify(self%reflist)
nullify(self%rltail)

self%nref = 0

end subroutine Delete_QCgvectorlist_

!--------------------------------------------------------------------------
recursive function get_QCnref_(self) result(nref)
!DEC$ ATTRIBUTES DLLEXPORT :: get_QCnref_
!! author: MDG
!! version: 1.0
!! date: 02/05/22
!!
!! return the number of QC g-vectors

IMPLICIT NONE

class(QCgvectors_T), INTENT(INOUT)  :: self
integer(kind=irg)                   :: nref

nref = self%nref

end function get_QCnref_

!--------------------------------------------------------------------------
recursive subroutine AddQCReflection_(self, QCcell, QCDiff, hkl)
!DEC$ ATTRIBUTES DLLEXPORT :: AddQCReflection_
  !! author: MDG
  !! version: 1.0
  !! date: 02/05/22
  !!
  !! add a reflection to the linked reflection list

use mod_io
use mod_QCcrystallography
use mod_QCdiffraction

IMPLICIT NONE

class(QCgvectors_T),INTENT(INOUT)                 :: self
class(QCcell_T),INTENT(INOUT)                     :: QCcell
type(QCDiffraction_T),INTENT(INOUT)               :: QCDiff
integer(kind=irg),INTENT(IN)                      :: hkl(:)

type(IO_T)                                        :: Message
integer(kind=irg)                                 :: istat, QCindex

! create linked list if it does not already exist
 if (.not.associated(self%reflist)) then
   call self%MakeQCRefList()
 end if 

! create a new entry
 if (.not.associated(self%rltail%next)) allocate(self%rltail%next,stat=istat)               ! allocate new value
 if (istat.ne.0) call Message%printError('AddQCReflection',' unable to add new reflection')

 self%rltail => self%rltail%next                     ! tail points to new value
 nullify(self%rltail%next)                           ! nullify next in new value

 self%nref = self%nref + 1                           ! update reflection counter
 self%rltail%num = self%nref                         ! store reflection number

 select type (QCcell)
   class is (QCcell_Icosahedral_T)
     self%rltail%hkl6 = hkl                           ! store 6-component Miller indices
     QCindex = QCcell%getnDIndex(hkl)
     self%rltail%glen = QCcell%getvectorLength(hkl, 'P', 'r')
   class is (QCcell_Axial_T)
     self%rltail%hkl5 = hkl                           ! store 5-component Miller indices
     QCindex = QCcell%getnDIndex(hkl)
     self%rltail%glen = QCcell%getvectorLength(hkl, 'P', 'r')
 end select

 self%rltail%Ucg = QCDiff%getLUT( QCindex )          ! store Ucg  in the list
 self%rltail%qg = QCDiff%getLUTqg( QCindex )         ! store pi/qg  in the list

 nullify(self%rltail%nextw)
 nullify(self%rltail%nexts)

end subroutine AddQCReflection_

!--------------------------------------------------------------------------
recursive function Get_QCListHead_(self) result(listrootw)
!DEC$ ATTRIBUTES DLLEXPORT :: Get_QCListHead_
  !! author: MDG
  !! version: 1.0
  !! date: 02/05/22
  !!
  !! return the head of the linked list

IMPLICIT NONE

class(QCgvectors_T),INTENT(INOUT)     :: self
type(QCreflisttype),pointer           :: listrootw

listrootw => self%reflist

end function Get_QCListHead_

!--------------------------------------------------------------------------
recursive subroutine Initialize_QCReflectionList_(self, QCcell, QCSG, QCDiff, FN, k, nref, verbose)
!DEC$ ATTRIBUTES DLLEXPORT :: Initialize_QCReflectionList_
  !! author: MDG
  !! version: 1.0
  !! date: 02/04/20
  !!
  !! initialize the potential reflection list for a given wave vector

use mod_io
use mod_QCcrystallography
use mod_QCdiffraction
use mod_QCsymmetry

IMPLICIT NONE

class(QCgvectors_T),INTENT(INOUT)     :: self
class(QCCell_T),INTENT(INOUT)         :: QCcell
type(QCSpaceGroup_T),INTENT(INOUT)    :: QCSG
type(QCDiffraction_T),INTENT(INOUT)   :: QCDiff
real(kind=dbl),INTENT(IN)             :: FN(3)
real(kind=dbl),INTENT(IN)             :: k(3)
integer(kind=irg),INTENT(IN)          :: nref
logical,INTENT(IN),OPTIONAL           :: verbose

type(QCgvectors_T)                    :: rltail
type(IO_T)                            :: Message

integer(kind=irg)                     :: imh, imhz, m, i, minholz, RHOLZ, im, istat, N, &
                                        ig, numr, ir, irsel, i1, i2, i3, i4, i5, i6, QCindex
integer(kind=irg),allocatable         :: orbit(:,:), gg(:)
real(kind=sgl)                        :: dhkl, io_real(9), H, g3(3), g3n(3), FNg(3), ddt, s, kr(3), exer, &
                                        rBethe_i, rBethe_d, sgp, r_g, la, dval
integer(kind=irg)                     :: io_int(3), gshort(3), isym, nn
complex(kind=dbl)                     :: Ucg, qg
real(kind=dbl)                        :: Vmod, Vpmod, xig, xgp
logical                               :: isnew

! set the truncation parameters
  rBethe_i = QCDiff%getBetheParameter('c3')   ! if larger than this value, we ignore the reflection completely
  rBethe_d = QCDiff%getBetheParameter('sg')   ! excitation error cutoff for double diffraction reflections
  la = 1.0/QCDiff%getWaveLength()

select type (QCcell)
  class is (QCcell_axial_T)
    m = 5
    allocate(orbit(m,QCSG%getnsym()), gg(m))
    orbit = 0
    gg = (/ 0,0,0,0,0 /)
    call self%AddQCReflection(QCcell, QCDiff, gg)   ! this guarantees that 000 is always the first reflection
    self%rltail%sg = 0.0

    call QCcell%get_imax(imh, imhz)
    imh = imh/2
    imhz = imhz/2

! now compute |sg|/|U_g|/lambda for the other allowed reflections; if this parameter is less than
! the threshhold, rBethe_i, then add the reflection to the list of potential reflections
    i1l: do i1=-imh,imh
     i2l: do i2=-imh,imh
      i3l: do i3=-imh,imh
       i4l: do i4=-imh,imh
        i5l: do i5=-imhz,imhz

          if ((abs(i1)+abs(i2)+abs(i3)+abs(i4)+abs(i5)).ne.0) then  ! avoid double counting the origin
            gg      = (/ i1, i2, i3, i4, i5/)
            QCindex = QCcell%getnDindex(gg)
            sgp     = QCDiff%CalcsgQC(QCcell,dble(gg),k,FN)
            Ucg     = QCDiff%getLUT(QCindex)   !QCcell%LUT(gg(1),gg(2),gg(3),gg(4),gg(5))
            r_g     = la * abs(sgp)/abs(Ucg)
            
            if (r_g.le.rBethe_i) then 
              call self%AddQCReflection( QCcell, QCDiff, gg )
              self%rltail%sg   = sgp
              self%rltail%glen = QCcell%getvectorLength(gg, 'P', 'r')
            end if
          end if

        end do i5l
       end do i4l
      end do i3l
     end do i2l
    end do i1l

  class is (QCcell_icosahedral_T)
    m = 6
    allocate(orbit(m,QCSG%getnsym()), gg(m))
    orbit = 0
    gg = (/ 0,0,0,0,0 /)
    call self%AddQCReflection(QCcell, QCDiff, gg)   ! this guarantees that 000 is always the first reflection
    self%rltail%sg = 0.0

    imh = QCcell%get_imax()
    imh = imh/2

! now compute |sg|/|U_g|/lambda for the other allowed reflections; if this parameter is less than
! the threshhold, rBethe_i, then add the reflection to the list of potential reflections
    i1r: do i1=-imh,imh
     i2r: do i2=-imh,imh
      i3r: do i3=-imh,imh
       i4r: do i4=-imh,imh
        i5r: do i5=-imh,imh
         i6r: do i6=-imh,imh

          if ((abs(i1)+abs(i2)+abs(i3)+abs(i4)+abs(i5)+abs(i6)).ne.0) then  ! avoid double counting the origin
            gg      = (/ i1, i2, i3, i4, i5, i6/)
            QCindex = QCcell%getnDindex(gg)
            sgp     = QCDiff%CalcsgQC(QCcell,dble(gg),k,FN)
            Ucg     = QCDiff%getLUT(QCindex)   !QCcell%LUT(gg(1),gg(2),gg(3),gg(4),gg(5))
            r_g     = la * abs(sgp)/abs(Ucg)
            
            if (r_g.le.rBethe_i) then 
              call self%AddQCReflection( QCcell, QCDiff, gg )
              self%rltail%sg   = sgp
              self%rltail%glen = QCcell%getvectorLength(gg, 'P', 'r')
            end if
          end if

         end do i6r
        end do i5r
       end do i4r
      end do i3r
     end do i2r
    end do i1r

  end select

  if (present(verbose)) then
    if (verbose) then
      io_int(1) = self%nref
      call Message%WriteValue(' Length of the master list of reflections : ', io_int, 1, "(I8)")
    end if
  end if

end subroutine Initialize_QCReflectionList_


!--------------------------------------------------------------------------
recursive subroutine Apply_BethePotentialsQC_(self, QCcell, QCDiff, listrootw, nns, nnw)
!DEC$ ATTRIBUTES DLLEXPORT :: Apply_BethePotentialsQC_
  !! author: MDG
  !! version: 1.0
  !! date: 02/08/22
  !!
  !! tag weak and strong reflections in self%reflist
  !!
  !! This routine steps through the reflist linked list and
  !! determines for each reflection whether it is strong or weak or should be
  !! ignored.  Strong and weak reflections are then linked in a new list via
  !! the nexts and nextw pointers, along with the nns and nnw counters.
  !! This routine makes use of the BetheParameter variable in the QCDiff_T class.

use mod_io
use mod_QCcrystallography
use mod_QCdiffraction

IMPLICIT NONE

class(QCgvectors_T), INTENT(INOUT)             :: self
class(QCcell_T),INTENT(INOUT)                  :: QCcell
type(QCDiffraction_T), INTENT(INOUT)           :: QCDiff
type(QCreflisttype),pointer                    :: listrootw
integer(kind=irg),INTENT(OUT)                  :: nns
integer(kind=irg),INTENT(OUT)                  :: nnw

type(IO_T)                                     :: Message
type(QCreflisttype),pointer                    :: rl, lastw, lasts

integer(kind=irg),allocatable                  :: glist(:,:), gmh(:)
real(kind=dbl),allocatable                     :: rh(:)
integer(kind=irg)                              :: icnt, istat, ir, ih, m, QCindex
real(kind=dbl)                                 :: sgp, la, mm

complex(kind=dbl)                              :: Ugh, qg
real(kind=dbl)                                 :: Vmod, Vpmod, xig, xgp

nullify(lasts)
nullify(lastw)
nullify(rl)

! first we extract the list of g-vectors from reflist, so that we can compute
! all the g-h difference vectors
select type (QCcell)
  class is (QCcell_axial_T)
    m = 5
    allocate(glist(m,self%nref),rh(self%nref),stat=istat)
  class is (QCcell_icosahedral_T)
    m = 6
    allocate(glist(m,self%nref),rh(self%nref),stat=istat)
end select

rl => self%reflist%next
icnt = 0
do
  if (.not.associated(rl)) EXIT
  icnt = icnt+1
  if (m.eq.5) then 
    glist(1:m,icnt) = rl%hkl5(1:m)
  else
    glist(1:m,icnt) = rl%hkl6(1:m)
  end if
  rl => rl%next
end do

! initialize the strong and weak reflection counters
nns = 1
nnw = 0

! the first reflection is always strong
rl => self%reflist%next
rl%strong = .TRUE.
rl%weak = .FALSE.
lasts => rl
nullify(lasts%nextw)

la = 1.D0/QCDiff%getWaveLength()

! next we need to iterate through all reflections in glist and
! determine which category the reflection belongs to: strong, weak, ignore
irloop: do ir = 2,icnt
  rl => rl%next
  rh = 0.D0
  sgp = la * abs(rl%sg)
  do ih = 1,icnt
    gmh(1:m) = glist(1:m,ir) - glist(1:m,ih)
    QCindex = QCcell%getnDIndex( gmh )
    Ugh = QCDiff%getLUT(QCindex)
    if (abs(Ugh).eq.0.D0) then 
      rh(ih) = 10000.D0
    else
      rh(ih) = sgp/abs(Ugh)
    end if
  end do

! which category does reflection ir belong to ?
  mm = minval(rh)

! m > c2 => ignore this reflection
  if (mm.gt.QCDiff%getBetheParameter('c2')) then
    rl%weak = .FALSE.
    rl%strong = .FALSE.
    CYCLE irloop
  end if

! c1 < m < c2 => weak reflection
  if ((QCDiff%getBetheParameter('c1').lt.mm).and.(mm.le.QCDiff%getBetheParameter('c2'))) then
    if (nnw.eq.0) then
      listrootw => rl
      lastw => rl
    else
      lastw%nextw => rl
      lastw => rl
      nullify(lastw%nexts)
    end if
    rl%weak = .TRUE.
    rl%strong = .FALSE.
    nnw = nnw + 1
    CYCLE irloop
  end if

! m < c1 => strong
  if (mm.le.QCDiff%getBetheParameter('c1')) then
    lasts%nexts => rl
    nullify(lasts%nextw)
    lasts => rl
    rl%weak = .FALSE.
    rl%strong = .TRUE.
    nns = nns + 1
  end if
end do irloop

deallocate(glist, rh)

end subroutine Apply_BethePotentialsQC_



! !--------------------------------------------------------------------------
! recursive subroutine GetDynMat_(self, cell, Diff, listrootw, DynMat, nns, nnw, BlochMode, noNormAbs)
! !DEC$ ATTRIBUTES DLLEXPORT :: GetDynMat_
!   !! author: MDG
!   !! version: 1.0
!   !! date: 02/12/20
!   !!
!   !! compute the dynamical matrix, including Bethe potentials
!   !!
!   !! We compute the dynamical matrix as the structure matrix A, with
!   !! the q_g elements along the off-diagonal; the reason for this is the fact
!   !! that this approach leads to a dynamical matrix that is shift invariant.
!   !! A conversion to the Bloch wave dynamical matrix can be obtained by setting
!   !! the optional keyword BlochMode

! use mod_crystallography
! use mod_diffraction
! use mod_kvectors

! IMPLICIT NONE

! class(gvectors_T), INTENT(INOUT) :: self
! type(Cell_T),INTENT(INOUT)       :: cell
! type(Diffraction_T),INTENT(INOUT):: Diff
! type(reflisttype),pointer        :: listrootw
! integer(kind=irg),INTENT(IN)     :: nns
! complex(kind=dbl),INTENT(INOUT)  :: DynMat(nns,nns)
! integer(kind=irg),INTENT(IN)     :: nnw
! character(5),INTENT(IN),OPTIONAL :: BlochMode   ! 'Bloch' or 'Struc'
! logical,INTENT(IN),OPTIONAL      :: noNormAbs

! type(gnode)                      :: rlp
! complex(kind=dbl)                :: czero, ughp, uhph, weaksum, cv, Agh, Ahgp, Ahmgp, Ahg, weakdiagsum, pq0, Ahh, Agpgp, ccpi
! real(kind=dbl)                   :: weaksgsum, tpi, Pioxgp, mlambda
! real(kind=sgl)                   :: Upz
! integer(kind=sgl)                :: ir, ic, ll(3), istat, wc
! type(reflisttype),pointer        :: listroot, rlr, rlc, rlw
! character(1)                     :: AorD

! czero = cmplx(0.0,0.0,dbl)      ! complex zero
! tpi = 2.D0 * cPi
! ccpi = cmplx(cPi,0.0D0,dbl)
! mLambda = Diff%getWaveLength()

! nullify(rlr)
! nullify(rlc)
! nullify(rlw)

! ! if BlochMode is absent, then we compute the Bloch dynamical matrix D directly
! ! if Blochmode = Struc, we compute the structure matrix A directly
! ! if Blochmode = Bloch, we do compute the structure matrix A and convert it to D
! ! [in the absence of the BlochMode keyword, the dynamical matrix D
! ! will not be invariant to origin shifts; A, on the other hand, is always shift
! ! invariant, so that D derived from A will also be shift invariant]

! call Diff%setrlpmethod('WK')
! listroot => self%reflist

! AorD = 'D'
! if (present(Blochmode)) AorD = 'A'

! ! Standard Bloch wave mode
! if (AorD.eq.'D') then

!         DynMat = czero
!         call Diff%CalcUcg(cell, (/0,0,0/) )
!         rlp = Diff%getrlp()
!         Upz = rlp%Upmod
!         if (present(noNormAbs)) then
!           if (noNormAbs.eqv..TRUE.) then
!             Upz = 0.0
!           end if
!         end if
!         !Pioxgp = cPi/rlp%xgp

!         rlr => listroot%next
!         ir = 1
!         do
!           if (.not.associated(rlr)) EXIT
!           rlc => listroot%next
!           ic = 1
!           do
!           if (.not.associated(rlc)) EXIT
!           if (ic.ne.ir) then  ! not a diagonal entry
! ! here we need to do the Bethe corrections if necessary
!             if (nnw.ne.0) then
!               weaksum = czero
!               rlw => listrootw
!               do
!                if (.not.associated(rlw)) EXIT
!                ll = rlr%hkl - rlw%hkl
!                ughp = Diff%getLUT( ll )
!                ll = rlw%hkl - rlc%hkl
!                uhph = Diff%getLUT( ll )
!                weaksum = weaksum +  ughp * uhph *cmplx(1.D0/rlw%sg,0.0,dbl)
!                rlw => rlw%nextw
!               end do
! !        ! and correct the dynamical matrix element to become a Bethe potential coefficient
!               ll = rlr%hkl - rlc%hkl
!               DynMat(ir,ic) = Diff%getLUT( ll ) - cmplx(0.5D0*mlambda,0.0D0,dbl)*weaksum
!              else
!               ll = rlr%hkl - rlc%hkl
!               DynMat(ir,ic) = Diff%getLUT( ll )
!             end if
!           else  ! it is a diagonal entry, so we need the excitation error and the absorption length
! ! determine the total contribution of the weak beams
!             if (nnw.ne.0) then
!               weaksgsum = 0.D0
!               rlw => listrootw
!               do
!                if (.not.associated(rlw)) EXIT
!                 ll = rlr%hkl - rlw%hkl
!                 ughp = Diff%getLUT( ll )
!                 weaksgsum = weaksgsum +  abs(ughp)**2/rlw%sg
!                 rlw => rlw%nextw
!               end do
!               weaksgsum = weaksgsum * mLambda/2.D0
!               DynMat(ir,ir) = cmplx(2.D0*rlr%sg/mLambda-weaksgsum,Upz,dbl)
!             else
!               DynMat(ir,ir) = cmplx(2.D0*rlr%sg/mLambda,Upz,dbl)
!             end if

!            end if
!            rlc => rlc%nexts
!            ic = ic + 1
!           end do
!           rlr => rlr%nexts
!           ir = ir+1
!         end do

! else ! AorD = 'A' so we need to compute the structure matrix using LUTqg ...

! ! note that the factor of i pi is added in at the end...
!         DynMat = czero
!         call Diff%CalcUcg(cell, (/0,0,0/), applyqgshift=.FALSE. )
!         rlp = Diff%getrlp()
!         pq0 = cmplx(0.D0,1.D0/rlp%xgp,dbl)

!         rlr => listroot%next
!         ir = 1
!         do
!           if (.not.associated(rlr)) EXIT
!           rlc => listroot%next
!           ic = 1
!           do
!           if (.not.associated(rlc)) EXIT
!           if (ic.ne.ir) then  ! not a diagonal entry
! ! here we need to do the Bethe corrections if necessary
!             if (nnw.ne.0) then
!               weaksum = czero
!               rlw => listrootw
!               do
!                if (.not.associated(rlw)) EXIT
!                ll = rlr%hkl - rlw%hkl
!                Agh = Diff%getLUTqg( ll )
!                ll = rlw%hkl - rlc%hkl
!                Ahgp = Diff%getLUTqg( ll )
! ! denominator Ahh - Ag'g'
!                Ahh = cmplx(2.D0 * rlw%sg,0.D0,dbl) + pq0
!                Agpgp = cmplx(2.D0 * rlc%sg,0.D0,dbl) + pq0
!                weaksum = weaksum +  Agh * Ahgp / (Ahh - Agpgp)
!                rlw => rlw%nextw
!               end do
! ! and correct the dynamical matrix element to become a Bethe potential coefficient
!               ll = rlr%hkl - rlc%hkl
!               DynMat(ir,ic) = Diff%getLUTqg( ll )  -  weaksum
!              else
!               ll = rlr%hkl - rlc%hkl
!               DynMat(ir,ic) = Diff%getLUTqg( ll )
!             end if
!           else  ! it is a diagonal entry, so we need the excitation error and the absorption length
! ! determine the total contribution of the weak beams
!             if (nnw.ne.0) then
!               weakdiagsum = 0.D0
!               rlw => listrootw
!               do
!                if (.not.associated(rlw)) EXIT
!                 ll = rlr%hkl - rlw%hkl
!                 Agh = Diff%getLUTqg( ll )
!                 Ahg = Diff%getLUTqg(-ll )
! ! denominator Ahh - Agg
!                 Ahh = cmplx(2.D0 * rlw%sg,0.D0,dbl) + pq0
!                 Agpgp = cmplx(2.D0 * rlr%sg,0.D0,dbl) + pq0
!                 weakdiagsum = weakdiagsum +  Agh * Ahg  / (Ahh - Agpgp)
!                 rlw => rlw%nextw
!               end do
!               DynMat(ir,ir) = cmplx( 2.D0 * rlr%sg, 0.D0, dbl) + pq0 - weakdiagsum
!             else
!               DynMat(ir,ir) = cmplx( 2.D0 * rlr%sg, 0.D0,dbl) + pq0
!             end if

!            end if
!            rlc => rlc%nexts
!            ic = ic + 1
!           end do
!           rlr => rlr%nexts
!           ir = ir+1
!         end do
!         DynMat = DynMat * ccpi ! cmplx(cPi, 0.D0)
!         !DynMat = DynMat * cmplx(0.D0,cPi)
! end if
! !cv = cmplx(1.D0/cPi/mLambda,0.D0)
! !DynMat = DynMat * cv

! if (present(BlochMode)) then
!   if (BlochMode.eq.'Bloch') then
!     cv = cmplx(1.D0/cPi/mLambda,0.D0)
!     DynMat = DynMat * cv
!   end if
! end if

! end subroutine GetDynMat_

! !--------------------------------------------------------------------------
! recursive subroutine GetDynMatFull_(self, cell, Diff, listrootw, DynMat, nref)
! !DEC$ ATTRIBUTES DLLEXPORT :: GetDynMatFull_
!   !! author: MDG
!   !! version: 1.0
!   !! date: 03/01/20
!   !!
!   !! compute the Bloch wave dynamical matrix for the full dynamical case
!   !! [no Bethe perturbations]

! use mod_crystallography
! use mod_diffraction
! use mod_kvectors

! IMPLICIT NONE

! class(gvectors_T), INTENT(INOUT) :: self
! type(Cell_T),INTENT(INOUT)       :: cell
! type(Diffraction_T),INTENT(INOUT):: Diff
! type(reflisttype),pointer        :: listrootw
! integer(kind=irg),INTENT(IN)     :: nref
! complex(kind=dbl),INTENT(INOUT)  :: DynMat(nref,nref)

! type(gnode)                      :: rlp
! complex(kind=dbl)                :: czero, weaksum, qg0, pq0
! real(kind=dbl)                   :: mlambda, weaksgsum
! real(kind=sgl)                   :: Upz
! integer(kind=sgl)                :: ir, ic, ll(3), istat, wc
! type(reflisttype),pointer        :: listroot, rlr, rlc

! czero = cmplx(0.0,0.0,dbl)      ! complex zero

! mLambda = Diff%getWaveLength()

! nullify(rlr)
! nullify(rlc)

! call Diff%setrlpmethod('WK')
! listroot => self%reflist

! DynMat = czero
! call Diff%CalcUcg(cell, (/0,0,0/) )
! rlp = Diff%getrlp()
! Upz = rlp%Upmod

! rlr => listroot%next
! ir = 1
!     do
!         if (.not.associated(rlr)) EXIT
!         rlc => listroot%next
!         ic = 1
!         do
!             if (.not.associated(rlc)) EXIT
!             if (ic.ne.ir) then  ! not a diagonal entry
!                 ll = rlr%hkl - rlc%hkl
!                 DynMat(ir,ic) = Diff%getLUT(ll)
!             else
!                 DynMat(ir,ir) = cmplx(2.D0*rlr%sg/mLambda,Upz,dbl)
!             end if
!             rlc => rlc%nexts
!             ic = ic + 1
!         end do
!         rlr => rlr%nexts
!         ir = ir+1
!     end do

! !DynMat = DynMat * cmplx(cPi,0.D0)

! end subroutine GetDynMatFull_

! !--------------------------------------------------------------------------
! recursive subroutine GetDynMatKin_(self, cell, Diff, listrootw, DynMat, nns)
! !DEC$ ATTRIBUTES DLLEXPORT :: GetDynMatKin_
!   !! author: MDG
!   !! version: 1.0
!   !! date: 03/01/20
!   !!
!   !! compute the Bloch wave dynamical matrix for the kinematical case
!   !! [no Bethe perturbations]

! use mod_crystallography
! use mod_diffraction
! use mod_kvectors

! IMPLICIT NONE

! class(gvectors_T), INTENT(INOUT) :: self
! type(Cell_T),INTENT(INOUT)       :: cell
! type(Diffraction_T),INTENT(INOUT):: Diff
! type(reflisttype),pointer        :: listrootw
! integer(kind=irg),INTENT(IN)     :: nns
! complex(kind=dbl),INTENT(INOUT)  :: DynMat(nns,nns)

! type(gnode)                      :: rlp
! complex(kind=dbl)                :: czero
! real(kind=dbl)                   :: mlambda
! real(kind=sgl)                   :: Upz
! integer(kind=sgl)                :: ir, ll(3)
! type(reflisttype),pointer        :: listroot, rlr

! czero = cmplx(0.0,0.0,dbl)      ! complex zero
! mLambda = Diff%getWaveLength()

! nullify(rlr)

! call Diff%setrlpmethod('WK')
! listroot => self%reflist

! DynMat = czero
! call Diff%CalcUcg(cell, (/0,0,0/) )
! rlp = Diff%getrlp()
! Upz = rlp%Upmod

! rlr => listroot%next
! ir = 1
! do
!   if (.not.associated(rlr)) EXIT
!   if (ir.ne.1) then  ! not a diagonal entry
! ! we only need to fill the first column in the kinematical approximation ...
!     ll = rlr%hkl
!     DynMat(ir,1) = Diff%getLUT( ll )
!   end if
! ! add the diagonal entry
!   DynMat(ir,ir) = cmplx(2.D0*rlr%sg/mLambda,Upz,dbl)
!   rlr => rlr%nexts
!   ir = ir+1
! end do

! end subroutine GetDynMatKin_


! !--------------------------------------------------------------------------
! recursive subroutine GetDynMatDHW_(self, cell, Diff, listroot, rltmpa, rltmpb, DynMat, nn, DM, ga, gb)
! !DEC$ ATTRIBUTES DLLEXPORT :: GetDynMatDHW_
!   !! author: MDG
!   !! version: 1.0
!   !! date: 03/01/20
!   !!
!   !! compute the Darwin-Howie-Whelan dynamical matrix
!   !! [no Bethe perturbations]

! use mod_crystallography
! use mod_diffraction
! use mod_kvectors
! use mod_io

! IMPLICIT NONE

! class(gvectors_T), INTENT(INOUT) :: self
! type(Cell_T),INTENT(INOUT)       :: cell
! type(Diffraction_T),INTENT(INOUT):: Diff
! type(reflisttype),pointer        :: listroot, rltmpa, rltmpb
! integer(kind=irg),INTENT(IN)     :: nn, ga(3), gb(3)
! real(kind=sgl),INTENT(IN)        :: DM(2,2)
! complex(kind=dbl),INTENT(INOUT)  :: DynMat(nn,nn)
! type(IO_T)                       :: Message

! type(gnode)                      :: rlp
! complex(kind=dbl)                :: cone
! real(kind=dbl)                   :: mlambda
! real(kind=sgl)                   :: Upz, DD
! integer(kind=sgl)                :: ir, ll(3), ic
! type(reflisttype),pointer        :: rlr
! real(kind=sgl)                   :: X(2)

!  cone = cmplx(0.D0,cPi)

!  mLambda = Diff%getWaveLength()

!  DD = DM(1,1)*DM(2,2) - DM(1,2)*DM(2,1)

!  nullify(rlr)

!  !call Diff%setrlpmethod('WK')
!  rlp = Diff%getrlp()
!  listroot => self%reflist
!  rltmpa => listroot%next

!  ! ir is the row index
!    do ir=1,nn
!     rltmpb => listroot%next   ! point to the front of the list
!  ! ic is the column index
!     do ic=1,nn
!      if (ic.ne.ir) then  ! exclude the diagonal
!  ! compute Fourier coefficient of electrostatic lattice potential
!       call Diff%CalcUcg(cell,rltmpa%hkl - rltmpb%hkl)
!       rlp = Diff%getrlp()
!       DynMat(ir,ic) = rlp%qg * cone
!  !cmplx(- cPi * aimag(rlp%qg), cPi * real(rlp%qg),dbl)  ! and initialize the off-diagonal matrix element (including i Pi)
!      end if
!      rltmpb => rltmpb%next  ! move to next column-entry
!     end do
!  ! decompose this point w.r.t ga and gb
!     X(1) = cell%CalcDot(float(rltmpa%hkl),float(ga),'c')
!     X(2) = cell%CalcDot(float(rltmpa%hkl),float(gb),'c')
!     X = matmul(DM,X)/DD

!     rltmpa%nab(1:2) = int(X(1:2))
!     rltmpa => rltmpa%next   ! move to next row-entry
!    end do

!    call Message%printMessage(' --> Reference Darwin-Howie-Whelan matrix initialized',"(A/)")

! end subroutine GetDynMatDHW_




end module mod_QCgvectors
