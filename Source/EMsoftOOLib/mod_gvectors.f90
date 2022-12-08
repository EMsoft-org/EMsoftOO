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
module mod_gvectors
  !! author: MDG
  !! version: 1.0
  !! date: 02/02/20
  !!
  !! variables and types needed to determine lists of reciprocal lattice vectors.
  !!
  !! This was part of the dynamical module, but was moved into a separate
  !! module.  The new module includes a routine to delete the linked list, and also routines
  !! to allocate linked lists, which are used by almost all dynamical scattering codes.
  !!
  !! Due to some complicated module interdependencies the CalcDynMat routine is in
  !! this module rather than in diffraction, where it would logically belong.  We may need
  !! to figure out how to change that.

use mod_kinds
use mod_global

IMPLICIT NONE
private

! [note added on 1/10/14]
! first we define the reflisttype (all the information needed for a given reciprocal lattice point or rlp).
! In this linked list, we want to keep everything that might be needed to perform rlp-related
! simulations, except for the Fourier coefficient of the lattice potential, wich is kept in
! a lookup table.  Anything that can easily be derived from the LUT does not need to be stored.
! [end note]
!
! linked list of reflections
type, public :: reflisttype
  integer(kind=irg)             :: num, &        ! sequential number
                                   hkl(3),&      ! Miller indices
                                   famhkl(3),&   ! family representative Miller indices
                                   HOLZN,&       ! belongs to this HOLZ layer
                                   strongnum,&   ! sequential number for strong beams
                                   weaknum,&     ! sequential number for weak beams
                                   famnum,&      ! family number
                                   nab(2)        ! decomposition with respect to ga and gb
  logical                       :: dbdiff        ! double diffraction reflection ?
  real(kind=dbl)                :: sg, &         ! excitation error
                                   xg, &         ! extinction distance
! removed 1/10/14                 Ucgmod, &      ! modulus of Fourier coefficient
                                   sangle, &     ! scattering angle (mrad)
                                   thetag        ! phase angle, needed for ECCI simulations
  logical                       :: strong, weak  ! is this a strong beam or not; both .FALSE. means 'do not consider'
  integer(kind=irg)             :: variant       ! one of the four variants of a superlattice reflection in fcc
  complex(kind=dbl)             :: Ucg           ! Fourier coefficient, copied from cell%LUT
  complex(kind=dbl)             :: qg            ! scaled Fourier coefficient, copied from cell%LUTqg
  type(reflisttype),pointer     :: next          ! connection to next entry in master linked list
  type(reflisttype),pointer     :: nexts         ! connection to next strong entry in linked list
  type(reflisttype),pointer     :: nextw         ! connection to next weak entry in linked list
end type reflisttype

type, public :: gvectors_T
private
  type(reflisttype), pointer :: reflist
  type(reflisttype), pointer :: rltail
  integer(kind=irg)          :: nref

contains
private
  procedure, pass(self) :: MakeRefList_
  procedure, pass(self) :: AddReflection_
  procedure, pass(self) :: Apply_BethePotentials_
  procedure, pass(self) :: GetSubRefList_
  procedure, pass(self) :: Delete_gvectorlist_
  procedure, pass(self) :: Compute_ReflectionListZoneAxis_
  procedure, pass(self) :: Initialize_ReflectionList_
  procedure, pass(self) :: Initialize_ReflectionList_PED_
  procedure, pass(self) :: Initialize_ReflectionList_EwaldSweep_
  procedure, pass(self) :: GetDynMat_
  procedure, pass(self) :: GetDynMatKin_
  procedure, pass(self) :: GetDynMatDHW_
  procedure, pass(self) :: GetDynMatFull_
  procedure, pass(self) :: GetExpthetagh_
  procedure, pass(self) :: get_nref_
  procedure, pass(self) :: CalcLgh_
  procedure, pass(self) :: CalcLghdepth_
  procedure, pass(self) :: CalcCBEDint_
  procedure, pass(self) :: CalcPEDint_
  procedure, pass(self) :: getSghfromLUT_
  procedure, pass(self) :: getSghfromLUTsum_
  procedure, pass(self) :: GetSgArray_ECCI_
  procedure, pass(self) :: GetExpval_ECCI_
  procedure, pass(self) :: Get_ListHead_
  procedure, pass(self) :: PrintRefList_

  final :: gvectors_destructor

  generic, public :: MakeRefList => MakeRefList_
  generic, public :: AddReflection => AddReflection_
  generic, public :: Apply_BethePotentials => Apply_BethePotentials_
  generic, public :: GetSubRefList => GetSubRefList_
  generic, public :: Delete_gvectorlist => Delete_gvectorlist_
  generic, public :: Compute_ReflectionListZoneAxis => Compute_ReflectionListZoneAxis_
  generic, public :: Initialize_ReflectionList => Initialize_ReflectionList_, &
                                                  Initialize_ReflectionList_PED_,&
                                                  Initialize_ReflectionList_EwaldSweep_
  generic, public :: GetDynMat => GetDynMat_
  generic, public :: GetDynMatKin => GetDynMatKin_
  generic, public :: GetDynMatDHW => GetDynMatDHW_
  generic, public :: GetDynMatFull => GetDynMatFull_
  generic, public :: GetExpthetagh => GetExpthetagh_
  generic, public :: get_nref => get_nref_
  generic, public :: CalcLgh => CalcLgh_
  generic, public :: CalcLghdepth => CalcLghdepth_
  generic, public :: CalcCBEDint => CalcCBEDint_
  generic, public :: CalcPEDint => CalcPEDint_
  generic, public :: getSghfromLUT => getSghfromLUT_
  generic, public :: getSghfromLUTsum => getSghfromLUTsum_
  generic, public :: GetSgArray_ECCI => GetSgArray_ECCI_
  generic, public :: GetExpval_ECCI => GetExpval_ECCI_
  generic, public :: Get_ListHead => Get_ListHead_
  generic, public :: PrintRefList => PrintRefList_

end type gvectors_T


! the constructor routine for this class
interface gvectors_T
  module procedure gvectors_constructor
end interface gvectors_T

contains

!--------------------------------------------------------------------------
type(gvectors_T) function gvectors_constructor( ) result(GVec)
!DEC$ ATTRIBUTES DLLEXPORT :: gvectors_constructor
!! author: MDG
!! version: 1.0
!! date: 02/02/20
!!
!! constructor for the gvectors_T Class

IMPLICIT NONE

! the calling program must make sure that the reflist is empty ...
! initialize the reflist
nullify(GVec%reflist)
GVec%nref = 0
call GVec%MakeRefList()

end function gvectors_constructor

!--------------------------------------------------------------------------
subroutine gvectors_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: gvectors_destructor
!! author: MDG
!! version: 1.0
!! date: 02/02/20
!!
!! destructor for the gvectors_T Class

IMPLICIT NONE

type(gvectors_T), INTENT(INOUT)  :: self

type(reflisttype),pointer        :: rltmpa

call reportDestructor('gvectors_T')

end subroutine gvectors_destructor

!--------------------------------------------------------------------------
recursive subroutine MakeRefList_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: MakeRefList_
  !! author: MDG
  !! version: 1.0
  !! date: 02/02/20
  !!
  !! allocate and initialize the linked reflection list

use mod_io

IMPLICIT NONE

class(gvectors_T), INTENT(INOUT)  :: self

type(IO_T)                        :: Message
integer(kind=irg)                 :: istat

! create it if it does not already exist
if (.not.associated(self%reflist)) then
  allocate(self%reflist,stat=istat)
  if (istat.ne.0) call Message%printError('MakeRefList:',' unable to allocate pointer')
  self%rltail => self%reflist           ! tail points to new value
  nullify(self%rltail%next)             ! nullify next in new value
end if

end subroutine MakeRefList_

!--------------------------------------------------------------------------
recursive subroutine Delete_gvectorlist_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: Delete_gvectorlist_
  !! author: MDG
  !! version: 1.0
  !! date: 02/02/20
  !!
  !! delete the entire linked list

IMPLICIT NONE

class(gvectors_T), INTENT(INOUT)  :: self

type(reflisttype),pointer         :: rltmpa

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

end subroutine Delete_gvectorlist_

!--------------------------------------------------------------------------
recursive function get_nref_(self) result(nref)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nref_
!! author: MDG
!! version: 1.0
!! date: 02/16/20
!!
!! return the number of g-vectors

IMPLICIT NONE

class(gvectors_T), INTENT(INOUT)  :: self
integer(kind=irg)                 :: nref

nref = self%nref

end function get_nref_

!--------------------------------------------------------------------------
recursive subroutine AddReflection_(self, Diff, hkl)
!DEC$ ATTRIBUTES DLLEXPORT :: AddReflection_
  !! author: MDG
  !! version: 1.0
  !! date: 02/02/20
  !!
  !! add a reflection to the linked reflection list

use mod_io
use mod_crystallography
use mod_diffraction

IMPLICIT NONE

class(gvectors_T),INTENT(INOUT)    :: self
type(Diffraction_T),INTENT(INOUT)  :: Diff
integer(kind=irg),INTENT(IN)       :: hkl(3)
 !! Miller indices of reflection to be added to list

type(IO_T)                         :: Message
integer(kind=irg)                  :: istat

! create linked list if it does not already exist
 if (.not.associated(self%reflist)) then
   call self%MakeRefList()
   self%rltail => self%reflist
   nullify(self%rltail%next)
   write (*,*) '.not.associated(self%reflist)'
 end if

! create a new entry
 if (.not.associated(self%rltail%next)) allocate(self%rltail%next,stat=istat)               ! allocate new value
 if (istat.ne.0) call Message%printError('AddReflection',' unable to add new reflection')

 self%rltail => self%rltail%next                          ! tail points to new value
 nullify(self%rltail%next)                           ! nullify next in new value

 self%nref = self%nref + 1                      ! update reflection counter
 self%rltail%num = self%nref                              ! store reflection number
 self%rltail%hkl = hkl                               ! store Miller indices
 self%rltail%Ucg = Diff%getLUT( hkl )                ! store Ucg  in the list
 self%rltail%qg = Diff%getLUTqg( hkl )               ! store pi/qg  in the list
 self%rltail%famnum = 0                              ! init this value for Prune_ReflectionList
! self%rltail%Ucgmod = cabs(rlp%Ucg)                 ! added on 2/29/2012 for Bethe potential computations
! self%rltail%sangle = 1000.0*dble(CalcDiffAngle(hkl(1),hkl(2),hkl(3)))    ! added 4/18/2012 for EIC project HAADF/BF tomography simulations
 self%rltail%thetag = Diff%getVphase( hkl )          ! added 12/14/2013 for EMECCI program
 nullify(self%rltail%nextw)
 nullify(self%rltail%nexts)

end subroutine AddReflection_

!--------------------------------------------------------------------------
recursive subroutine GetSubReflist_(self, cell, Diff, FN, kk, listrootw, nns, nnw, first)
!DEC$ ATTRIBUTES DLLEXPORT :: GetSubReflist_
  !! author: MDG
  !! version: 1.0
  !! date: 02/02/20
  !!
  !! used to subselect a series of reflections from an existing reflist

use mod_crystallography
use mod_diffraction

IMPLICIT NONE

class(gvectors_T),INTENT(INOUT)                :: self
type(Cell_T),INTENT(INOUT)                     :: cell
type(Diffraction_T),INTENT(INOUT)              :: Diff
real(kind=sgl),INTENT(IN)                      :: FN(3)
real(kind=sgl),INTENT(IN)                      :: kk(3)
type(reflisttype),pointer                      :: listrootw
integer(kind=irg),INTENT(OUT)                  :: nns
integer(kind=irg),INTENT(OUT)                  :: nnw
logical,INTENT(IN),OPTIONAL                    :: first

integer(kind=irg),allocatable,SAVE             :: glist(:,:)
real(kind=dbl),allocatable,SAVE                :: rh(:)
integer(kind=irg),SAVE                         :: icnt
real(kind=dbl),SAVE                            :: la

type(reflisttype),pointer                      :: rl, lastw, lasts
integer(kind=irg)                              :: istat, gmh(3), ir, ih
real(kind=dbl)                                 :: sgp, m, sg, gg(3)


! this routine must go through the entire reflist and apply the Bethe criteria for a
! given beam and foil normal pair. Hence, it is similar to the Apply_BethePotentials
! routine, except that the excitation errors are recomputed each time this routine
! is called, to account for the different incident beam directions.

nullify(lasts)
nullify(lastw)
nullify(rl)

! first we extract the list of g-vectors from reflist, so that we can compute
! all the g-h difference vectors; we only need to do this the first time we call
! this routine, since this list never changes (hence the SAVE qualifier for glist)
if (PRESENT(first)) then
  if (first) then
    la = 1.D0/Diff%getWaveLength()
    allocate(glist(3,self%nref),rh(self%nref),stat=istat)
    rl => self%reflist%next
    icnt = 0
    do
      if (.not.associated(rl)) EXIT
      icnt = icnt+1
      glist(1:3,icnt) = rl%hkl(1:3)
      rl => rl%next
    end do
  end if
else
end if

! initialize the strong and weak reflection counters
nns = 1
nnw = 0

! the first reflection is always strong
rl => self%reflist%next
rl%strong = .TRUE.
rl%weak = .FALSE.
lasts => rl
nullify(lasts%nextw)


! next we need to iterate through all reflections in glist and
! determine which category the reflection belongs to: strong, weak, ignore
irloop: do ir = 2,icnt
  rl => rl%next
  rh = 0.D0
! here we need to actually compute the excitation error sg, since kk and FN are different each time the routine is called
  gg = dble(rl%hkl)
  sg = Diff%Calcsg(cell,gg,dble(kk),dble(FN))
  sgp = la * abs(sg)

  do ih = 1,icnt
   gmh(1:3) = glist(1:3,ir) - glist(1:3,ih)
   if (Diff%getdbdiff( gmh )) then  ! it is a double diffraction reflection with |U|=0
! to be written
     rh(ih) = 10000.D0
   else
     rh(ih) = sgp/abs( Diff%getLUT( gmh ) )
   end if
  end do

! which category does reflection ir belong to ?
  m = minval(rh)

! m > c2 => ignore this reflection
  if (m.gt.Diff%getBetheParameter('c2')) then
    rl%weak = .FALSE.
    rl%strong = .FALSE.
    CYCLE irloop
  end if

! c1 < m < c2 => weak reflection
  if ((Diff%getBetheParameter('c1').lt.m).and.(m.le.Diff%getBetheParameter('c2'))) then
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
  if (m.le.Diff%getBetheParameter('c1')) then
    lasts%nexts => rl
    nullify(lasts%nextw)
    lasts => rl
    rl%weak = .FALSE.
    rl%strong = .TRUE.
    nns = nns + 1
  end if
end do irloop

end subroutine GetSubRefList_

!--------------------------------------------------------------------------
recursive function Get_ListHead_(self) result(listroot)
!DEC$ ATTRIBUTES DLLEXPORT :: Get_ListHead_
  !! author: MDG
  !! version: 1.0
  !! date: 02/02/20
  !!
  !! return the head of the linked list


IMPLICIT NONE

class(gvectors_T),INTENT(INOUT)                :: self
type(reflisttype),pointer                      :: listroot

listroot => self%reflist

end function Get_ListHead_

!--------------------------------------------------------------------------
recursive subroutine PrintRefList_(self, swall)
!DEC$ ATTRIBUTES DLLEXPORT :: PrintRefList_
  !! author: MDG
  !! version: 1.0
  !! date: 07/12/22
  !!
  !! print the reflection list 

use mod_io

IMPLICIT NONE

class(gvectors_T),INTENT(INOUT)               :: self
character(1),INTENT(IN)                       :: swall

type(IO_T)                                    :: Message 

type(reflisttype),pointer                     :: rltmp
integer(kind=irg)                             :: io_int(3)

rltmp => self%reflist%next

do 
  if (.not.associated(rltmp)) EXIT 
  io_int(1:3) = rltmp%hkl(1:3)
  call Message%WriteValue('',io_int,3)
  select case (swall) 
    case('a')  
      rltmp => rltmp%next 
    case('s')  
      rltmp => rltmp%nexts
    case('w')  
      rltmp => rltmp%nextw 
  end select
end do 

end subroutine PrintRefList_

!--------------------------------------------------------------------------
recursive subroutine Apply_BethePotentials_(self, Diff, listrootw, nns, nnw, verbose)
!DEC$ ATTRIBUTES DLLEXPORT :: Apply_BethePotentials_
  !! author: MDG
  !! version: 1.0
  !! date: 02/02/20
  !!
  !! tag weak and strong reflections in self%reflist
  !!
  !! This routine steps through the reflist linked list and
  !! determines for each reflection whether it is strong or weak or should be
  !! ignored.  Strong and weak reflections are then linked in a new list via
  !! the nexts and nextw pointers, along with the nns and nnw counters.
  !! This routine makes use of the BetheParameter variable in the Diff_T class.

use mod_io
use mod_diffraction

IMPLICIT NONE

class(gvectors_T), INTENT(INOUT)               :: self
type(Diffraction_T), INTENT(INOUT)             :: Diff
type(reflisttype),pointer                      :: listrootw
integer(kind=irg),INTENT(OUT)                  :: nns
integer(kind=irg),INTENT(OUT)                  :: nnw
logical,INTENT(IN),OPTIONAL                    :: verbose

type(IO_T)                                     :: Message
integer(kind=irg),allocatable                  :: glist(:,:)
real(kind=dbl),allocatable                     :: rh(:)
type(reflisttype),pointer                      :: rl, lastw, lasts
integer(kind=irg)                              :: icnt, istat, gmh(3), ir, ih, io_int(1)
real(kind=dbl)                                 :: sgp, la, m

nullify(lasts)
nullify(lastw)
nullify(rl)

! first we extract the list of g-vectors from reflist, so that we can compute
! all the g-h difference vectors
allocate(glist(3,self%nref),rh(self%nref),stat=istat)
rl => self%reflist%next
icnt = 0
do
  if (.not.associated(rl)) EXIT
  icnt = icnt+1
  glist(1:3,icnt) = rl%hkl(1:3)
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

la = 1.D0/Diff%getWaveLength()

! next we need to iterate through all reflections in glist and
! determine which category the reflection belongs to: strong, weak, ignore
irloop: do ir = 2,icnt
  rl => rl%next
  rh = 0.D0
  sgp = la * abs(rl%sg)
  do ih = 1,icnt
   gmh(1:3) = glist(1:3,ir) - glist(1:3,ih)
   if (Diff%getdbdiff( gmh )) then  ! it is a double diffraction reflection with |U|=0
! to be written
     rh(ih) = 10000.D0
   else
     rh(ih) = sgp/abs( Diff%getLUT( gmh ) )
   end if
  end do

! which category does reflection ir belong to ?
  m = minval(rh)

! m > c2 => ignore this reflection
  if (m.gt.Diff%getBetheParameter('c2')) then
    rl%weak = .FALSE.
    rl%strong = .FALSE.
    CYCLE irloop
  end if

! c1 < m < c2 => weak reflection
  if ((Diff%getBetheParameter('c1').lt.m).and.(m.le.Diff%getBetheParameter('c2'))) then
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
  if (m.le.Diff%getBetheParameter('c1')) then
    lasts%nexts => rl
    nullify(lasts%nextw)
    lasts => rl
    rl%weak = .FALSE.
    rl%strong = .TRUE.
    nns = nns + 1
  end if
end do irloop

deallocate(glist, rh)

if (present(verbose)) then 
  if (verbose.eqv..TRUE.) then 
    io_int(1) = nns 
    call Message%WriteValue(' Number of strong beams : ', io_int, 1)
  end if 
end if 

end subroutine Apply_BethePotentials_

!--------------------------------------------------------------------------
recursive subroutine Compute_ReflectionListZoneAxis_(self,cell,SG,Diff,FN,dmin,k,ga,gb)
!DEC$ ATTRIBUTES DLLEXPORT :: Compute_ReflectionListZoneAxis_
  !! author: MDG
  !! version: 1.0
  !! date: 02/02/20
  !!
  !! compute the entire reflection list for general conditions (including HOLZ)

use mod_io
use mod_crystallography
use mod_diffraction
use mod_symmetry

IMPLICIT NONE

class(gvectors_T),INTENT(INOUT)   :: self
type(Cell_T),INTENT(INOUT)        :: cell
type(SpaceGroup_T),INTENT(INOUT)  :: SG
type(Diffraction_T),INTENT(INOUT) :: Diff
real(kind=sgl),INTENT(IN)         :: FN(3)
real(kind=sgl),INTENT(IN)         :: dmin
real(kind=sgl),INTENT(IN)         :: k(3)
integer(kind=irg),INTENT(IN)      :: ga(3)
integer(kind=irg),INTENT(IN)      :: gb(3)

type(IO_T)                        :: Message
integer(kind=irg)                 :: imh, imk, iml, gg(3), ix, iy, iz, i, minholz, RHOLZ, im, istat, N, &
                                     ig, numr, ir, irsel
real(kind=sgl)                    :: dhkl, io_real(6), H, g3(3), g3n(3), FNg(3), ddt, s, kr(3), exer, &
                                     rBethe_i, rBethe_d, sgp, r_g, la, dval
integer(kind=irg)                 :: io_int(3), gshort(3), gp(3)

! set the truncation parameters
rBethe_i = Diff%getBetheParameter('c3')  ! if larger than this value, we ignore the reflection completely
rBethe_d = Diff%getBetheParameter('sg')  ! excitation error cutoff for double diffraction reflections
la = 1.0/sngl(Diff%getWaveLength())

! get the size of the lookup table
gp = Diff%getShapeLUT()
imh = (gp(1)-1)/4
imk = (gp(2)-1)/4
iml = (gp(3)-1)/4

! nullify(self%reflist)
! nullify(self%rltail)

gg = (/0,0,0/)
call self%AddReflection(Diff, gg)  ! this guarantees that 000 is always the first reflection

self%rltail%sg = 0.0
! now compute |sg|/|U_g|/lambda for the other allowed reflections; if this parameter is less than
! the threshhold, rBethe_i, then add the reflection to the list of potential reflections
! note that this uses the older form of the Bethe Potential truncation parameters for now

do ix=-imh,imh
  do iy=-imk,imk
    gg = ix*ga + iy*gb
    if ((abs(gg(1))+abs(gg(2))+abs(gg(3))).ne.0) then  ! avoid double counting the origin
      dval = 1.0/cell%CalcLength(float(gg), 'r' )
      if ((SG%IsGAllowed(gg)).AND.(dval .gt. dmin)) then
        sgp = Diff%Calcsg(cell, float(gg),k,FN)
        if  ((abs(gg(1)).le.imh).and.(abs(gg(2)).le.imk).and.(abs(gg(3)).le.iml) ) then
          if (Diff%getdbdiff( gg )) then ! potential double diffraction reflection
            if (abs(sgp).le.rBethe_d) then
              call self%AddReflection(Diff, gg)
              self%rltail%sg = sgp
              self%rltail%dbdiff = .TRUE.
            end if
          else
            r_g = la * abs(sgp)/abs(Diff%getLUT( gg ))
            if (r_g.le.rBethe_i) then
              call self%AddReflection(Diff, gg )
              self%rltail%sg = sgp
              self%rltail%dbdiff = .FALSE.
            end if
          end if
        end if
      end if
    end if
  end do
end do
io_int(1) = self%nref
call Message%WriteValue(' Length of the master list of reflections : ', io_int, 1, "(I5,/)")

end subroutine Compute_ReflectionListZoneAxis_

!--------------------------------------------------------------------------
recursive subroutine Initialize_ReflectionList_(self, cell, SG, Diff, FN, k, dmin, verbose)
!DEC$ ATTRIBUTES DLLEXPORT :: Initialize_ReflectionList_
  !! author: MDG
  !! version: 1.0
  !! date: 02/04/20
  !!
  !! initialize the potential reflection list for a given wave vector

use mod_io
use mod_crystallography
use mod_diffraction
use mod_symmetry

IMPLICIT NONE

class(gvectors_T),INTENT(INOUT)     :: self
type(Cell_T),INTENT(INOUT)          :: cell
type(SpaceGroup_T),INTENT(INOUT)    :: SG
type(Diffraction_T),INTENT(INOUT)   :: Diff
real(kind=sgl),INTENT(IN)           :: FN(3)
real(kind=sgl),INTENT(IN)           :: k(3)
real(kind=sgl),INTENT(IN)           :: dmin
logical,INTENT(IN),OPTIONAL         :: verbose

type(IO_T)                          :: Message
integer(kind=irg)                   :: imh, imk, iml, gg(3), ix, iy, iz, i, minholz, RHOLZ, im, istat, N, &
                                       ig, numr, ir, irsel
real(kind=sgl)                      :: dhkl, io_real(9), H, g3(3), g3n(3), FNg(3), ddt, s, kr(3), exer, &
                                       rBethe_i, rBethe_d, sgp, r_g, la, dval
integer(kind=irg)                   :: io_int(3), gshort(3), gp(3)

! set the truncation parameters
  rBethe_i = Diff%getBetheParameter('c3')   ! if larger than this value, we ignore the reflection completely
  rBethe_d = Diff%getBetheParameter('sg')   ! excitation error cutoff for double diffraction reflections
  la = 1.0/Diff%getWaveLength()

! get the size of the lookup table
  gp = Diff%getshapeLUT()
  imh = (gp(1)-1)/4
  imk = (gp(2)-1)/4
  iml = (gp(3)-1)/4

! transmitted beam has excitation error zero
  gg = (/ 0,0,0 /)
  call self%AddReflection(Diff, gg )   ! this guarantees that 000 is always the first reflection
  self%rltail%sg = 0.0

! now compute |sg|/|U_g|/lambda for the other allowed reflections; if this parameter is less than
! the threshhold, rBethe_i, then add the reflection to the list of potential reflections
ixl: do ix=-imh,imh
iyl:  do iy=-imk,imk
izl:   do iz=-iml,iml
         if ((abs(ix)+abs(iy)+abs(iz)).ne.0) then  ! avoid double counting the origin
          gg = (/ ix, iy, iz /)
          dval = 1.0/cell%CalcLength(float(gg), 'r' )
          if ((SG%IsGAllowed(gg)).and.(dval.gt.dmin)) then ! allowed by the lattice centering, if any
            sgp = Diff%Calcsg(cell,float(gg),k,FN)
            if (Diff%getdbdiff( gg )) then ! potential double diffraction reflection
              if (abs(sgp).le.rBethe_d) then
                call self%AddReflection(Diff, gg )
                self%rltail%sg = sgp
                self%rltail%dbdiff = .TRUE.
              end if
            else
              r_g = la * abs(sgp)/abs(Diff%getLUT(gg))
              if (r_g.le.rBethe_i) then
                call self%AddReflection(Diff, gg )
                self%rltail%sg = sgp
                self%rltail%dbdiff = .FALSE.
              end if
            end if
          end if ! IsGAllowed
         end if
       end do izl
      end do iyl
    end do ixl

  if (present(verbose)) then
    if (verbose) then
      io_int(1) = self%nref
      call Message%WriteValue(' Length of the master list of reflections : ', io_int, 1, "(I8)")
    end if
  end if

end subroutine Initialize_ReflectionList_

!--------------------------------------------------------------------------
recursive subroutine Initialize_ReflectionList_PED_(self, cell, SG, Diff, gmax, Igmax, verbose)
!DEC$ ATTRIBUTES DLLEXPORT :: Initialize_ReflectionList_PED_
  !! author: MDG
  !! version: 1.0
  !! date: 02/04/20
  !!
  !! initialize the potential reflection list for a given wave vector

use mod_io
use mod_crystallography
use mod_diffraction
use mod_symmetry

IMPLICIT NONE

class(gvectors_T),INTENT(INOUT)     :: self
type(Cell_T),INTENT(INOUT)          :: cell
type(SpaceGroup_T),INTENT(INOUT)    :: SG
type(Diffraction_T),INTENT(INOUT)   :: Diff
real(kind=sgl),INTENT(IN)           :: gmax
real(kind=sgl),INTENT(OUT)          :: Igmax
logical,INTENT(IN),OPTIONAL         :: verbose

type(IO_T)                          :: Message
integer(kind=irg)                   :: imh, imk, iml, gg(3), ix, iy, iz, i, minholz, RHOLZ, im, istat, N, &
                                       ig, numr, ir, irsel
real(kind=sgl)                      :: dhkl, io_real(9), H, g3(3), g3n(3), FNg(3), ddt, s, kr(3), exer, &
                                       rBethe_i, rBethe_d, sgp, r_g, la, dval, xgmin, glen
integer(kind=irg)                   :: io_int(3), gshort(3), gp(3)


! get the size of the lookup table
  gp = Diff%getshapeLUT()
  imh = (gp(1)-1)/4
  imk = (gp(2)-1)/4
  iml = (gp(3)-1)/4

! transmitted beam has excitation error zero
  gg = (/ 0,0,0 /)
  call self%AddReflection(Diff, gg )   ! this guarantees that 000 is always the first reflection
  self%rltail%xg = 0.0
  xgmin = 100000.0
  Igmax = 0.0

! now compute |sg|/|U_g|/lambda for the other allowed reflections; if this parameter is less than
! the threshhold, rBethe_i, then add the reflection to the list of potential reflections
ixl: do ix=-imh,imh
iyl:  do iy=-imk,imk
izl:   do iz=-iml,iml
         if ((abs(ix)+abs(iy)+abs(iz)).ne.0) then  ! avoid double counting the origin
           gg = (/ ix, iy, iz /)
           glen = cell%CalcLength(float(gg), 'r' )
  
  ! find all reflections, ignoring double diffraction spots
           if ((SG%IsGAllowed(gg)).and.(glen.le.gmax)) then ! allowed by the lattice centering, if any
              call self%AddReflection(Diff, gg)
  ! we'll use the sangle field of the rltail structure to store |Ug|^2; we will also need the extinction distance
              self%rltail%sangle = abs(Diff%getLUT((/ix, iy, iz/)))**2
              ! print*,Diff%getLUT((/ix, iy, iz/))
              if (self%rltail%sangle.gt.Igmax) Igmax = self%rltail%sangle
              self%rltail%xg = 1.0/(abs(Diff%getLUT((/ix, iy, iz/)))*Diff%getWaveLength())
              if ( self%rltail%xg.lt.xgmin) xgmin =  self%rltail%xg
           end if ! IsGAllowed
          end if
       end do izl
      end do iyl
    end do ixl

  if (present(verbose)) then
    if (verbose) then
      io_int(1) = self%nref
      call Message%WriteValue(' Length of the master list of reflections : ', io_int, 1, "(I8)")
    end if
  end if

end subroutine Initialize_ReflectionList_PED_

!--------------------------------------------------------------------------
recursive subroutine Initialize_ReflectionList_EwaldSweep_(self,cell,SG,Diff,FN,k,pedangle, goffset, verbose)
!DEC$ ATTRIBUTES DLLEXPORT :: Initialize_ReflectionList_EwaldSweep_
  !! author: MDG
  !! version: 1.0
  !! date: 02/04/20
  !!
  !! initialize the potential reflection list for a given precession electron diffraction geometry

use mod_io
use mod_crystallography
use mod_diffraction
use mod_symmetry

IMPLICIT NONE

class(gvectors_T),INTENT(INOUT)     :: self
type(Cell_T),INTENT(INOUT)          :: cell
type(SpaceGroup_T),INTENT(INOUT)    :: SG
type(Diffraction_T),INTENT(INOUT)   :: Diff
real(kind=sgl),INTENT(IN)           :: FN(3)
real(kind=sgl),INTENT(IN)           :: k(3)
real(kind=sgl),INTENT(IN)           :: pedangle
real(kind=sgl),INTENT(IN)           :: goffset
logical,INTENT(IN),OPTIONAL         :: verbose

type(IO_T)                          :: Message
integer(kind=irg)                   :: imh, imk, iml, gg(3), ix, iy, iz, io_int(3), gp(3)
real(kind=sgl)                      :: FNg(3), c, s, kr(3), sgp, la, kstar(3), gperp(3), gpara(3), bup, blo, y, z, &
                                       gdk, glen, gplen

! init a couple of parameters
  la = 1.0/Diff%getWaveLength()
  c = la * cos(pedangle/1000.0)
  s = 2.0 * la * sin(pedangle/1000.0)

! reciprocal space wave vector
  call cell%TransSpace(k, kstar, 'd', 'r')
  call cell%NormVec(kstar, 'r')
  kstar = la * kstar

! get the size of the lookup table
  gp = Diff%getshapeLUT()
  imh = (gp(1)-1)/4
  imk = (gp(2)-1)/4
  iml = (gp(3)-1)/4

  if (associated(self%reflist)) then
    call self%Delete_gvectorlist()
  end if

! transmitted beam has excitation error zero, and set xg to zero; xg will store the accumulated intensity for each reflection
  gg = (/ 0,0,0 /)
  call self%AddReflection(Diff, gg )   ! this guarantees that 000 is always the first reflection
  self%rltail%sg = 0.0
  self%rltail%xg = 0.0

! scan through the reciprocal lattice volume corresponding to the dmin value
ixl: do ix=-imh,imh
iyl:  do iy=-imk,imk
izl:   do iz=-iml,iml
        if ((abs(ix)+abs(iy)+abs(iz)).ne.0) then  ! avoid double counting the origin
         gg = (/ ix, iy, iz /)
         if (SG%IsGAllowed(gg)) then ! allowed by the lattice centering, if any
! first we need to determine the parallel and perpendicular components of this g vector with respect to the beam direction in reciprocal space
          gdk = cell%CalcDot(float(gg),kstar,'r')       ! projection of gg onto k*
          gpara = gdk * kstar
          gperp = float(gg) - gpara
! then get the length of the perpendicular and parallel components, including sign of parallel component
          glen = cell%CalcLength(gperp, 'r')
          gplen = cell%CalcLength(gpara, 'r')
        ! sign of length depends on dot product sign of gg onto k*
          if (gdk.le.0.0) gplen = -gplen
! compute the upper and lower bounds for this value of glen
          y = glen*s
          z = c*c-glen*glen
          bup = goffset+c-sqrt(z-y)
          blo = -goffset+c-sqrt(z+y)
! and check whether or not this point should be taken into account
          if ((blo.le.gplen).and.(gplen.le.bup)) then
            sgp = Diff%Calcsg(cell,float(gg),k,FN)
! note that we are not applying any Bethe parameter conditions here since those will be applied for each beam orientation separately
            if (Diff%getdbdiff( gg )) then ! potential double diffraction reflection
                call self%AddReflection(Diff, gg )
                self%rltail%sg = sgp
                self%rltail%xg = 0.0
                self%rltail%dbdiff = .TRUE.
            else
                call self%AddReflection(Diff, gg )
                self%rltail%sg = sgp
                self%rltail%xg = 0.0
                self%rltail%dbdiff = .FALSE.
            end if
          end if ! reflection inside precession-swept Ewald sphere volume
         end if ! IsGAllowed
        end if ! not the origin
       end do izl
      end do iyl
    end do ixl

  if (present(verbose)) then
    if (verbose) then
      io_int(1) = self%nref
      call Message%WriteValue(' Length of the master list of reflections : ', io_int, 1, "(I8)")
    end if
  end if

end subroutine Initialize_Reflectionlist_EwaldSweep_

!--------------------------------------------------------------------------
recursive subroutine GetDynMat_(self, cell, Diff, listrootw, DynMat, nns, nnw, MatrixType, noNormAbs, noDiagonal)
!DEC$ ATTRIBUTES DLLEXPORT :: GetDynMat_
  !! author: MDG
  !! version: 1.0
  !! date: 02/12/20
  !!
  !! compute the dynamical matrix, including Bethe potentials
  !!
  !! We compute the dynamical matrix as the structure matrix A, with
  !! the q_g elements along the off-diagonal; the reason for this is the fact
  !! that this approach leads to a dynamical matrix that is shift invariant.
  !! A conversion to the Bloch wave dynamical matrix can be obtained by setting
  !! the optional keyword BlochMode

use mod_crystallography
use mod_diffraction
use mod_kvectors

IMPLICIT NONE

class(gvectors_T), INTENT(INOUT) :: self
type(Cell_T),INTENT(INOUT)       :: cell
type(Diffraction_T),INTENT(INOUT):: Diff
type(reflisttype),pointer        :: listrootw
integer(kind=irg),INTENT(IN)     :: nns
complex(kind=dbl),INTENT(INOUT)  :: DynMat(nns,nns)
integer(kind=irg),INTENT(IN)     :: nnw
character(5),INTENT(IN),OPTIONAL :: MatrixType   ! 'Bloch' or 'Struc'
logical,INTENT(IN),OPTIONAL      :: noNormAbs
logical,INTENT(IN),OPTIONAL      :: noDiagonal

type(gnode)                      :: rlp
complex(kind=dbl)                :: czero, ughp, uhph, weaksum, cv, Agh, Ahgp, Ahmgp, Ahg, weakdiagsum, pq0, Ahh, Agpgp, ccpi
real(kind=dbl)                   :: weaksgsum, tpi, Pioxgp, mlambda
real(kind=dbl)                   :: Upz
integer(kind=sgl)                :: ir, ic, ll(3), istat, wc
type(reflisttype),pointer        :: listroot, rlr, rlc, rlw
character(1)                     :: AorD

czero = cmplx(0.0,0.0,dbl)      ! complex zero
tpi = 2.D0 * cPi
ccpi = cmplx(cPi,0.0D0,dbl)
mLambda = Diff%getWaveLength()

nullify(rlr)
nullify(rlc)
nullify(rlw)

! if MatrixType is absent, then we compute the Bloch dynamical matrix D directly
! if MatrixType = Struc, we compute the structure matrix A directly

call Diff%setrlpmethod('WK')
listroot => self%reflist

AorD = 'D'
if (present(MatrixType)) then
  if (MatrixType.eq.'Struc') AorD = 'A'
end if 

! Standard Bloch wave mode
if (AorD.eq.'D') then
        DynMat = czero
        call Diff%CalcUcg(cell, (/0,0,0/) )
        rlp = Diff%getrlp()
        Upz = rlp%Upmod
        if (present(noNormAbs)) then
          if (noNormAbs.eqv..TRUE.) then
            Upz = 0.0
          end if
        end if
        !Pioxgp = cPi/rlp%xgp

        rlr => listroot%next
        ir = 1
        do
          if (.not.associated(rlr)) EXIT
          rlc => listroot%next
          ic = 1
          do
          if (.not.associated(rlc)) EXIT
          if (ic.ne.ir) then  ! not a diagonal entry
! here we need to do the Bethe corrections if necessary
            if (nnw.ne.0) then
              weaksum = czero
              rlw => listrootw
              do
               if (.not.associated(rlw)) EXIT
               ll = rlr%hkl - rlw%hkl
               ughp = Diff%getLUT( ll )
               ll = rlw%hkl - rlc%hkl
               uhph = Diff%getLUT( ll )
               weaksum = weaksum +  ughp * uhph *cmplx(1.D0/rlw%sg,0.0,dbl)
               rlw => rlw%nextw
              end do
!        ! and correct the dynamical matrix element to become a Bethe potential coefficient
              ll = rlr%hkl - rlc%hkl
              DynMat(ir,ic) = Diff%getLUT( ll ) - cmplx(0.5D0*mlambda,0.0D0,dbl)*weaksum
             else
              ll = rlr%hkl - rlc%hkl
              DynMat(ir,ic) = Diff%getLUT( ll )
            end if
          else  ! it is a diagonal entry, so we need the excitation error and the absorption length
! determine the total contribution of the weak beams
            if (nnw.ne.0) then
              weaksgsum = 0.D0
              rlw => listrootw
              do
               if (.not.associated(rlw)) EXIT
                ll = rlr%hkl - rlw%hkl
                ughp = Diff%getLUT( ll )
                weaksgsum = weaksgsum +  abs(ughp)**2/rlw%sg
                rlw => rlw%nextw
              end do
              weaksgsum = weaksgsum * mLambda/2.D0
              DynMat(ir,ir) = cmplx(2.D0*rlr%sg/mLambda-weaksgsum,Upz,dbl)
            else
              DynMat(ir,ir) = cmplx(2.D0*rlr%sg/mLambda,Upz,dbl)
            end if

           end if
           rlc => rlc%nexts
           ic = ic + 1
          end do
          rlr => rlr%nexts
          ir = ir+1
        end do

else ! AorD = 'A' so we need to compute the structure matrix using LUTqg ...

! note that the factor of i pi is added in at the end...
        DynMat = czero
        call Diff%CalcUcg(cell, (/0,0,0/), applyqgshift=.FALSE. )
        rlp = Diff%getrlp()
        pq0 = cmplx(0.D0,1.D0/rlp%xgp,dbl)

        rlr => listroot%next
        ir = 1
        do
          if (.not.associated(rlr)) EXIT
          rlc => listroot%next
          ic = 1
          do
          if (.not.associated(rlc)) EXIT
          if (ic.ne.ir) then  ! not a diagonal entry
! here we need to do the Bethe corrections if necessary
            if (nnw.ne.0) then
              weaksum = czero
              rlw => listrootw
              do
               if (.not.associated(rlw)) EXIT
               ll = rlr%hkl - rlw%hkl
               Agh = Diff%getLUTqg( ll )
               ll = rlw%hkl - rlc%hkl
               Ahgp = Diff%getLUTqg( ll )
! denominator Ahh - Ag'g'
               Ahh = cmplx(2.D0 * rlw%sg,0.D0,dbl) + pq0
               Agpgp = cmplx(2.D0 * rlc%sg,0.D0,dbl) + pq0
               weaksum = weaksum +  Agh * Ahgp / (Ahh - Agpgp)
               rlw => rlw%nextw
              end do
! and correct the dynamical matrix element to become a Bethe potential coefficient
              ll = rlr%hkl - rlc%hkl
              DynMat(ir,ic) = Diff%getLUTqg( ll )  -  weaksum
             else
              ll = rlr%hkl - rlc%hkl
              DynMat(ir,ic) = Diff%getLUTqg( ll )
            end if
          else  ! it is a diagonal entry, so we need the excitation error and the absorption length
! determine the total contribution of the weak beams
            if (nnw.ne.0) then
              weakdiagsum = 0.D0
              rlw => listrootw
              do
               if (.not.associated(rlw)) EXIT
                ll = rlr%hkl - rlw%hkl
                Agh = Diff%getLUTqg( ll )
                Ahg = Diff%getLUTqg(-ll )
! denominator Ahh - Agg
                Ahh = cmplx(2.D0 * rlw%sg,0.D0,dbl) + pq0
                Agpgp = cmplx(2.D0 * rlr%sg,0.D0,dbl) + pq0
                weakdiagsum = weakdiagsum +  Agh * Ahg  / (Ahh - Agpgp)
                rlw => rlw%nextw
              end do
              DynMat(ir,ir) = cmplx( 2.D0 * rlr%sg, 0.D0, dbl) + pq0 - weakdiagsum
            else
              DynMat(ir,ir) = cmplx( 2.D0 * rlr%sg, 0.D0,dbl) + pq0
            end if

           end if
           rlc => rlc%nexts
           ic = ic + 1
          end do
          rlr => rlr%nexts
          ir = ir+1
        end do
        DynMat = DynMat * ccpi ! cmplx(cPi, 0.D0)
end if

if (present(noDiagonal)) then 
  if (noDiagonal.eqv..TRUE.) then 
    do ic=1,nns 
      DynMat(ic,ic) = czero 
    end do 
  end if 
end if 

end subroutine GetDynMat_

!--------------------------------------------------------------------------
recursive subroutine GetDynMatFull_(self, cell, Diff, listrootw, DynMat, nref)
!DEC$ ATTRIBUTES DLLEXPORT :: GetDynMatFull_
  !! author: MDG
  !! version: 1.0
  !! date: 03/01/20
  !!
  !! compute the Bloch wave dynamical matrix for the full dynamical case
  !! [no Bethe perturbations]

use mod_crystallography
use mod_diffraction
use mod_kvectors

IMPLICIT NONE

class(gvectors_T), INTENT(INOUT) :: self
type(Cell_T),INTENT(INOUT)       :: cell
type(Diffraction_T),INTENT(INOUT):: Diff
type(reflisttype),pointer        :: listrootw
integer(kind=irg),INTENT(IN)     :: nref
complex(kind=dbl),INTENT(INOUT)  :: DynMat(nref,nref)

type(gnode)                      :: rlp
complex(kind=dbl)                :: czero, weaksum, qg0, pq0
real(kind=dbl)                   :: mlambda, weaksgsum
real(kind=sgl)                   :: Upz
integer(kind=sgl)                :: ir, ic, ll(3), istat, wc
type(reflisttype),pointer        :: listroot, rlr, rlc

czero = cmplx(0.0,0.0,dbl)      ! complex zero

mLambda = Diff%getWaveLength()

nullify(rlr)
nullify(rlc)

call Diff%setrlpmethod('WK')
listroot => self%reflist

DynMat = czero
call Diff%CalcUcg(cell, (/0,0,0/) )
rlp = Diff%getrlp()
Upz = rlp%Upmod

rlr => listroot%next
ir = 1
    do
        if (.not.associated(rlr)) EXIT
        rlc => listroot%next
        ic = 1
        do
            if (.not.associated(rlc)) EXIT
            if (ic.ne.ir) then  ! not a diagonal entry
                ll = rlr%hkl - rlc%hkl
                DynMat(ir,ic) = Diff%getLUT(ll)
            else
                DynMat(ir,ir) = cmplx(2.D0*rlr%sg/mLambda,Upz,dbl)
            end if
            rlc => rlc%nexts
            ic = ic + 1
        end do
        rlr => rlr%nexts
        ir = ir+1
    end do

!DynMat = DynMat * cmplx(cPi,0.D0)

end subroutine GetDynMatFull_

!--------------------------------------------------------------------------
recursive subroutine GetDynMatKin_(self, cell, Diff, listrootw, DynMat, nns)
!DEC$ ATTRIBUTES DLLEXPORT :: GetDynMatKin_
  !! author: MDG
  !! version: 1.0
  !! date: 03/01/20
  !!
  !! compute the Bloch wave dynamical matrix for the kinematical case
  !! [no Bethe perturbations]

use mod_crystallography
use mod_diffraction
use mod_kvectors

IMPLICIT NONE

class(gvectors_T), INTENT(INOUT) :: self
type(Cell_T),INTENT(INOUT)       :: cell
type(Diffraction_T),INTENT(INOUT):: Diff
type(reflisttype),pointer        :: listrootw
integer(kind=irg),INTENT(IN)     :: nns
complex(kind=dbl),INTENT(INOUT)  :: DynMat(nns,nns)

type(gnode)                      :: rlp
complex(kind=dbl)                :: czero
real(kind=dbl)                   :: mlambda
real(kind=sgl)                   :: Upz
integer(kind=sgl)                :: ir, ll(3)
type(reflisttype),pointer        :: listroot, rlr

czero = cmplx(0.0,0.0,dbl)      ! complex zero
mLambda = Diff%getWaveLength()

nullify(rlr)

call Diff%setrlpmethod('WK')
listroot => self%reflist

DynMat = czero
call Diff%CalcUcg(cell, (/0,0,0/) )
rlp = Diff%getrlp()
Upz = rlp%Upmod

rlr => listroot%next
ir = 1
do
  if (.not.associated(rlr)) EXIT
  if (ir.ne.1) then  ! not a diagonal entry
! we only need to fill the first column in the kinematical approximation ...
    ll = rlr%hkl
    DynMat(ir,1) = Diff%getLUT( ll )
  end if
! add the diagonal entry
  DynMat(ir,ir) = cmplx(2.D0*rlr%sg/mLambda,Upz,dbl)
  rlr => rlr%nexts
  ir = ir+1
end do

end subroutine GetDynMatKin_


!--------------------------------------------------------------------------
recursive subroutine GetDynMatDHW_(self, cell, Diff, DynMat, nn, DM, ga, gb)
!DEC$ ATTRIBUTES DLLEXPORT :: GetDynMatDHW_
  !! author: MDG
  !! version: 1.0
  !! date: 03/01/20
  !!
  !! compute the Darwin-Howie-Whelan dynamical matrix
  !! [no Bethe perturbations]

use mod_crystallography
use mod_diffraction
use mod_kvectors
use mod_io

IMPLICIT NONE

class(gvectors_T), INTENT(INOUT) :: self
type(Cell_T),INTENT(INOUT)       :: cell
type(Diffraction_T),INTENT(INOUT):: Diff
integer(kind=irg),INTENT(IN)     :: nn, ga(3), gb(3)
real(kind=sgl),INTENT(IN)        :: DM(2,2)
complex(kind=dbl),INTENT(INOUT)  :: DynMat(nn,nn)
type(IO_T)                       :: Message

type(reflisttype),pointer        :: listroot, rltmpa, rltmpb
type(gnode)                      :: rlp
complex(kind=dbl)                :: cone
real(kind=dbl)                   :: mlambda
real(kind=sgl)                   :: Upz, DD
integer(kind=sgl)                :: ir, ll(3), ic
type(reflisttype),pointer        :: rlr
real(kind=sgl)                   :: X(2)

 cone = cmplx(0.D0,cPi)

 mLambda = Diff%getWaveLength()

 DD = DM(1,1)*DM(2,2) - DM(1,2)*DM(2,1)

 nullify(rlr)

 !call Diff%setrlpmethod('WK')
 rlp = Diff%getrlp()
 listroot => self%reflist
 rltmpa => listroot%next

 ! ir is the row index
   do ir=1,nn
    rltmpb => listroot%next   ! point to the front of the list
 ! ic is the column index
    do ic=1,nn
     if (ic.ne.ir) then  ! exclude the diagonal
 ! compute Fourier coefficient of electrostatic lattice potential
      call Diff%CalcUcg(cell,rltmpa%hkl - rltmpb%hkl)
      rlp = Diff%getrlp()
      DynMat(ir,ic) = rlp%qg * cone
 !cmplx(- cPi * aimag(rlp%qg), cPi * real(rlp%qg),dbl)  ! and initialize the off-diagonal matrix element (including i Pi)
     end if
     rltmpb => rltmpb%next  ! move to next column-entry
    end do
 ! decompose this point w.r.t ga and gb
    X(1) = cell%CalcDot(float(rltmpa%hkl),float(ga),'c')
    X(2) = cell%CalcDot(float(rltmpa%hkl),float(gb),'c')
    X = matmul(DM,X)/DD

    rltmpa%nab(1:2) = int(X(1:2))
    rltmpa => rltmpa%next   ! move to next row-entry
   end do

   call Message%printMessage(' --> Reference Darwin-Howie-Whelan matrix initialized',"(A/)")

end subroutine GetDynMatDHW_


!--------------------------------------------------------------------------
recursive subroutine GetSgArray_ECCI_(self, cell, sg, kl, FN, Diff, nns)
!DEC$ ATTRIBUTES DLLEXPORT :: GetSgArray_ECCI_
  !! author: MDG
  !! version: 1.0
  !! date: 03/01/20
  !!
  !! compute the Sg array for ECCI

use mod_crystallography
use mod_diffraction
use mod_kvectors
use mod_io

IMPLICIT NONE

class(gvectors_T), INTENT(INOUT) :: self
type(Cell_T),INTENT(INOUT)       :: cell
type(Diffraction_T),INTENT(INOUT):: Diff
integer(kind=irg),INTENT(IN)     :: nns
real(kind=sgl),INTENT(INOUT)     :: sg(nns)
real(kind=sgl),INTENT(IN)        :: kl(3)
real(kind=sgl),INTENT(IN)        :: FN(3)

type(reflisttype),pointer        :: rltmpa
type(IO_T)                       :: Message
integer(kind=irg)                :: ig, gg(3)

  nullify(rltmpa)
  rltmpa => self%reflist
  rltmpa => rltmpa%next
  reflectionloopCL: do ig=1,nns
    gg = float(rltmpa%hkl)
    sg(ig) = Diff%Calcsg(cell,float(gg),kl,FN)
! ! and we move to the next strong reflection in the list
    rltmpa => rltmpa%nexts
  end do reflectionloopCL

end subroutine GetSgArray_ECCI_

!--------------------------------------------------------------------------
recursive subroutine GetExpval_ECCI_(self, cell, expval, Diff, nn, DM, ga, gb)
!DEC$ ATTRIBUTES DLLEXPORT :: GetExpval_ECCI_
  !! author: MDG
  !! version: 1.0
  !! date: 03/01/20
  !!
  !! compute the exponential values array used in ECCI computation

use mod_crystallography
use mod_diffraction
use mod_kvectors
use mod_io

IMPLICIT NONE

class(gvectors_T), INTENT(INOUT) :: self
type(Cell_T),INTENT(INOUT)       :: cell
type(Diffraction_T),INTENT(INOUT):: Diff
integer(kind=irg),INTENT(IN)     :: nn, ga(3), gb(3)
integer(kind=sgl),INTENT(INOUT)  :: expval(2, nn, nn)
real(kind=sgl),INTENT(IN)        :: DM(2,2)

type(IO_T)                       :: Message

type(reflisttype),pointer        :: listroot, rltmpa, rltmpb
integer(kind=irg)                :: ir, ic
real(kind=sgl)                   :: X(2)
real(kind=sgl)                   :: DD

  DD = DM(1,1)*DM(2,2) - DM(1,2)*DM(2,1)

  expval = 0
  listroot => self%reflist
  rltmpa => listroot%next

! this part has already been done in the Apply_BethePotentials routine ...
  ! ir is the row index
    do ir=1,nn
     ! rltmpb => listroot%next   ! point to the front of the list
  ! ic is the column index
     ! do ic=1,nn
     !  rltmpb => rltmpb%next  ! move to next column-entry
     ! end do
  ! decompose this point w.r.t ga and gb
     X(1) = cell%CalcDot(float(rltmpa%hkl),float(ga),'c')
     X(2) = cell%CalcDot(float(rltmpa%hkl),float(gb),'c')
     X = matmul(DM,X)/DD

     rltmpa%nab(1:2) = int(X(1:2))
     rltmpa => rltmpa%nexts   ! move to next row-entry
    end do
  nullify(rltmpa)
  rltmpa => listroot%next

  ! ir is the row index
    do ir=1,nn
      rltmpb => listroot%next   ! point to the front of the list
  ! ic is the column index
       do ic=1,nn
         if (ic.ne.ir) then  ! exclude the diagonal
           expval(1,ir,ic) = rltmpb%nab(1)-rltmpa%nab(1)
           expval(2,ir,ic) = rltmpb%nab(2)-rltmpa%nab(2)
         end if
         rltmpb => rltmpb%nexts ! move to next column-entry
      end do
      rltmpa => rltmpa%nexts  ! move to next row-entry
    end do
end subroutine GetExpval_ECCI_


!--------------------------------------------------------------------------
recursive subroutine getSghfromLUT_(self,Diff,nns,numset,Sgh)
!DEC$ ATTRIBUTES DLLEXPORT :: getSghfromLUT_
  !! author: MDG
  !! version: 1.0
  !! date: 02/12/20
  !!
  !! compute structure factor-like Sgh array entry for EBSD, ECCI and ECP simulations

use mod_crystallography
use mod_diffraction
! use symmetry

IMPLICIT NONE

class(gvectors_T),INTENT(INOUT)         :: self
type(Diffraction_T),INTENT(INOUT)       :: Diff
integer(kind=irg),INTENT(IN)            :: nns
integer(kind=irg),INTENT(IN)            :: numset
complex(kind=dbl),INTENT(INOUT)         :: Sgh(nns,nns,numset)

type(reflisttype),pointer               :: rltmpa, rltmpb
integer(kind=irg)                       :: ir, ic, kkk(3)

associate( reflist => self%reflist )
! loop over all contributing reflections
! ir is the row index
    rltmpa => reflist%next    ! point to the front of the list
    do ir=1,nns
! ic is the column index
      rltmpb => reflist%next    ! point to the front of the list
      do ic=1,nns
        kkk = rltmpb%hkl - rltmpa%hkl
        Sgh(ir,ic,1:numset) = Diff%getSghLUT( numset, kkk )
        ! Sgh(ic,ir,1:numset) = Diff%getSghLUT( numset, kkk )
        rltmpb => rltmpb%nexts  ! move to next column-entry
      end do
     rltmpa => rltmpa%nexts  ! move to next row-entry
   end do
 end associate

end subroutine getSghfromLUT_

!--------------------------------------------------------------------------
recursive subroutine getSghfromLUTsum_(self,Diff,nns,numset,Sgh)
!DEC$ ATTRIBUTES DLLEXPORT :: getSghfromLUTsum_
  !! author: MDG
  !! version: 1.0
  !! date: 02/12/20
  !!
  !! compute structure factor-like Sgh array entry for EBSD, ECCI and ECP simulations

use mod_crystallography
use mod_diffraction
! use symmetry

IMPLICIT NONE

class(gvectors_T),INTENT(INOUT)         :: self
type(Diffraction_T),INTENT(INOUT)       :: Diff
integer(kind=irg),INTENT(IN)            :: nns
integer(kind=irg),INTENT(IN)            :: numset
complex(kind=dbl),INTENT(INOUT)         :: Sgh(nns,nns)

type(reflisttype),pointer               :: rltmpa, rltmpb
integer(kind=irg)                       :: ir, ic, kkk(3)
real(kind=dbl)                          :: thdiff

associate( reflist => self%reflist )

! loop over all contributing reflections
! ir is the row index
    rltmpa => reflist%next    ! point to the front of the list
    do ir=1,nns
! ic is the column index
      rltmpb => reflist%next    ! point to the front of the list
      do ic=1,nns
        kkk = rltmpb%hkl - rltmpa%hkl
        thdiff = rltmpb%thetag - rltmpa%thetag
        Sgh(ir,ic) = cmplx( cos(thdiff), sin(thdiff) ) * sum(Diff%getSghLUT( numset, kkk ))
        rltmpb => rltmpb%nexts  ! move to next column-entry
      end do
     rltmpa => rltmpa%nexts  ! move to next row-entry
   end do
 end associate

end subroutine getSghfromLUTsum_

!--------------------------------------------------------------------------
recursive subroutine GetExpthetagh_(self,Diff,nns,thetagh)
!DEC$ ATTRIBUTES DLLEXPORT :: GetExpthetagh_
  !! author: MDG
  !! version: 1.0
  !! date: 04/20/22
  !!
  !! compute phase factor thetagh array for scattering matrix ECCI simulations

use mod_crystallography
use mod_diffraction
! use symmetry

IMPLICIT NONE

class(gvectors_T),INTENT(INOUT)         :: self
type(Diffraction_T),INTENT(INOUT)       :: Diff
integer(kind=irg),INTENT(IN)            :: nns
complex(kind=dbl),INTENT(INOUT)         :: thetagh(nns,nns)

type(reflisttype),pointer               :: rltmpa, rltmpb
integer(kind=irg)                       :: ir, ic 
real(kind=dbl)                          :: thdiff

associate( reflist => self%reflist )

! loop over all contributing reflections
! ir is the row index
    rltmpa => reflist%next    ! point to the front of the list
    do ir=1,nns
! ic is the column index
      rltmpb => reflist%next    ! point to the front of the list
      do ic=1,nns
        thdiff = rltmpb%thetag - rltmpa%thetag
        thetagh(ir,ic) = cmplx( cos(thdiff), sin(thdiff) )
        rltmpb => rltmpb%nexts  ! move to next strong beam column-entry
      end do
     rltmpa => rltmpa%nexts  ! move to next strong beam row-entry
   end do
 end associate

end subroutine GetExpthetagh_

!--------------------------------------------------------------------------
recursive subroutine CalcCBEDint_(self,cell,DynMat,kn,nn,nt,thick,inten)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcCBEDint_
  !! author: MDG
  !! version: 1.0
  !! date: 02/12/20
  !!
  !! compute the scattered intensities for a range of thicknesses (specifically for CBED mode)

use mod_io
use mod_crystallography
use mod_diffraction
use mod_kvectors
use mod_global

IMPLICIT NONE

class(gvectors_T),INTENT(INOUT) :: self
integer(kind=irg),INTENT(IN)    :: nn                   !< number of strong beams
complex(kind=dbl),INTENT(IN)    :: DynMat(nn,nn)
type(Cell_T)                    :: cell
real(kind=sgl),INTENT(IN)       :: kn
integer(kind=irg),INTENT(IN)    :: nt                   !< number of thickness values
real(kind=sgl),INTENT(IN)       :: thick(nt)            !< thickness array
real(kind=sgl),INTENT(INOUT)    :: inten(nt,nn)         !< output intensities (both strong and weak)

type(Diffraction_T)             :: Diff
integer(kind=irg)               :: i,j,IPIV(nn), ll(3), jp
complex(kind=dbl)               :: CGinv(nn,nn), Minp(nn,nn),diag(nn),Wloc(nn), lCG(nn,nn), lW(nn), &
                                   lalpha(nn), delta(nn,nn) 
real(kind=sgl)                  :: th

! compute the eigenvalues and eigenvectors
 Minp = DynMat
 IPIV = 0
 call Diff%BWsolve(Minp,Wloc,lCG,CGinv,nn,IPIV)

! the alpha coefficients are in the first column of the inverse matrix
! the minus sign in W(i) stems from the fact that k_n is in the direction
! opposite to the foil normal
 lW = cPi*Wloc/cmplx(kn,0.0)
 do i=1,nn
  lalpha(i) = CGinv(i,1)
 end do

! we are going to ignore weak beam intensities for now...

! compute the strong beam intensities, stored in the first nn slots of inten 
! we could also compute the weak beams, since they make use of the same diag(1:nn) expression
! as the strong beams, plus a few other factors (excitation error, wave length, Fourier coefficients)
! that part would need to rewritten entirely
 do i=1,nt
  th = thick(i)
  diag(1:nn)=exp(-th*aimag(lW(1:nn)))*cmplx(cos(th*real(lW(1:nn))),sin(th*real(lW(1:nn))))*lalpha(1:nn)
! the delta array is common to the strong and weak beam intensity computation, so we compute it first
  do j=1,nn
   delta(j,1:nn) = lCG(j,1:nn)*diag(1:nn)
  end do
! strong beams
  do j=1,nn
   inten(i,j) = abs(sum(delta(j,1:nn)))**2
  end do 
 end do

end subroutine CalcCBEDint_

recursive subroutine CalcPEDint_(self,cell,DynMat,kn,nn,nt,thick,inten)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcPEDint_
  !! author: MDG
  !! version: 1.0
  !! date: 02/12/20
  !!
  !! compute the scattered intensities for a range of thicknesses (specifically for CBED mode)

use mod_io
use mod_crystallography
use mod_diffraction
use mod_kvectors
use mod_global

IMPLICIT NONE

class(gvectors_T),INTENT(INOUT) :: self
integer(kind=irg),INTENT(IN)    :: nn                   !< number of strong beams
complex(kind=dbl),INTENT(IN)    :: DynMat(nn,nn)
type(Diffraction_T) :: Diff
type(Cell_T)                    :: cell
real(kind=sgl),INTENT(IN)       :: kn

integer(kind=irg),INTENT(IN)    :: nt                   !< number of thickness values
real(kind=sgl),INTENT(IN)       :: thick(nt)            !< thickness array
real(kind=sgl),INTENT(INOUT)    :: inten(nt,nn)         !< output intensities (both strong and weak)

integer(kind=irg)               :: i,j,IPIV(nn), ll(3), jp
complex(kind=dbl)               :: CGinv(nn,nn), Minp(nn,nn),diag(nn),Wloc(nn), lCG(nn,nn), lW(nn), &
                                   lalpha(nn), delta(nn,nn) 
real(kind=sgl)                  :: th

! compute the eigenvalues and eigenvectors
Minp = DynMat
IPIV = 0
call Diff%BWsolve(Minp,Wloc,lCG,CGinv,nn,IPIV)

! the alpha coefficients are in the first column of the inverse matrix
lW = cPi*Wloc/cmplx(kn,0.0)
lalpha(1:nn) = CGinv(1:nn,1)

! make sure the alpha excitation coefficients are normalized 
! s = sum(abs(lalpha(1:nn))**2)
! if (s.ne.1.D0) then
!  s = cmplx(1.D0/dsqrt(s),0.D0)
!  lalpha = lalpha*s
! endif 

! compute the strong beam intensities, stored in the first nn slots of inten 
do i=1,nt
 th = thick(i)
 diag(1:nn)=exp(-th*aimag(lW(1:nn)))*cmplx(cos(th*real(lW(1:nn))),sin(th*real(lW(1:nn))))*lalpha(1:nn)
! the delta array is common to the strong and weak beam intensity computation, so we compute it first
 do j=1,nn
  delta(j,1:nn) = lCG(j,1:nn)*diag(1:nn)
 end do
! strong beams
 do j=1,nn
  inten(i,j) = abs(sum(delta(j,1:nn)))**2
 end do 
end do
end subroutine CalcPEDint_

!--------------------------------------------------------------------------
recursive subroutine CalcLgh_(self,DMat,Lgh,thick,kn,nn,gzero,depthstep,lambdaE,izz)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcLgh_
  !! author: MDG
  !! version: 1.0
  !! date: 02/12/20
  !!
  !! compute the Lgh matrix for EBSD, ECCI, ECP, etc simulations

use mod_io

IMPLICIT NONE

class(gvectors_T),INTENT(INOUT)     :: self
integer(kind=irg),INTENT(IN)        :: nn
complex(kind=dbl),INTENT(IN)        :: DMat(nn,nn)
complex(kind=dbl),INTENT(OUT)       :: Lgh(nn,nn)
real(kind=dbl),INTENT(IN)           :: thick
real(kind=dbl),INTENT(IN)           :: kn
integer(kind=irg),INTENT(IN)        :: gzero
real(kind=dbl),INTENT(IN)           :: depthstep
integer(kind=irg),INTENT(IN)        :: izz
real(kind=sgl),INTENT(IN)           :: lambdaE(izz)

type(IO_T)                          :: Message
integer                             :: i,j,k, iz
complex(kind=dbl)                   :: CGinv(nn,nn), Minp(nn,nn), tmp3(nn,nn)

real(kind=dbl)                      :: tpi, dzt
complex(kind=dbl)                   :: Ijk(nn,nn), q, getMIWORK(1), qold

integer(kind=irg)                   :: INFO, LDA, LDVR, LDVL,  JPIV(nn), MILWORK
complex(kind=dbl)                   :: CGG(nn,nn), W(nn)
complex(kind=dbl),allocatable       :: MIWORK(:)

integer(kind=irg),parameter         :: LWMAX = 5000
complex(kind=dbl)                   :: VL(nn,nn),  WORK(LWMAX)
real(kind=dbl)                      :: RWORK(2*nn)
character                           :: JOBVL, JOBVR
integer(kind=sgl)                   :: LWORK
integer(kind=sgl)                   :: ilaenv

! compute the eigenvalues and eigenvectors
! using the LAPACK CGEEV, CGETRF, and CGETRI routines
!
 Minp = DMat

! set some initial LAPACK variables
 LDA = nn
 LDVL = nn
 LDVR = nn
 INFO = 0

 ! first initialize the parameters for the LAPACK ZGEEV, CGETRF, and CGETRI routines
 JOBVL = 'N'   ! do not compute the left eigenvectors
 JOBVR = 'V'   ! do compute the right eigenvectors
 LWORK = -1 ! so that we can ask the routine for the actually needed value

! call the routine to determine the optimal workspace size
 LDA = nn
  call zgeev(JOBVL,JOBVR,nn,Minp,LDA,W,VL,LDVL,CGG,LDVR,WORK,LWORK,RWORK,INFO)
  LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

! then call the eigenvalue solver
  LDA = nn
  call zgeev(JOBVL,JOBVR,nn,Minp,LDA,W,VL,LDVL,CGG,LDVR,WORK,LWORK,RWORK,INFO)
  if (INFO.ne.0) call Message%printError('Error in CalcLgh: ','ZGEEV return not zero')

 CGinv = CGG

 LDA=nn
 call zgetrf(nn,nn,CGinv,LDA,JPIV,INFO)
 MILWORK = -1
 LDA=nn
! 
! MILWORK=ilaenv( 1, 'ZGETRI', ' ', nn, -1, -1, -1 )
! if (.not.allocated(MIWORK)) allocate(MIWORK(nn*MILWORK))
 
 call zgetri(nn,CGinv,LDA,JPIV,getMIWORK,MILWORK,INFO)
 MILWORK =  INT(real(getMIWORK(1)))
 if (.not.allocated(MIWORK)) allocate(MIWORK(MILWORK))
 MIWORK = cmplx(0.D0,0.D0)
 LDA=nn
 call zgetri(nn,CGinv,LDA,JPIV,MIWORK,MILWORK,INFO)
 deallocate(MIWORK)

! in all the time that we've used these routines, we haven't
! had a single problem with the matrix inversion, so we don't
! really need to do this test:
!
! if ((cabs(sum(matmul(CGG,CGinv)))-dble(nn)).gt.1.E-8) write (*,*) 'Error in matrix inversion; continuing'


! then compute the integrated intensity matrix
 W = W/cmplx(2.0*kn,0.0)

! recall that alpha(1:nn) = CGinv(1:nn,gzero)

! first the Ijk matrix (this is Winkelmann's B^{ij}(t) matrix)
! combined with numerical integration over [0, z0] interval,
! taking into account depth profiles from Monte Carlo simulations ...
! the depth profile lambdaE must be added to the absorption
! components of the Bloch wave eigenvalues.

tpi = 2.D0*cPi*depthstep
dzt = depthstep/thick
 do j=1,nn
  do k=1,nn
     q =  cmplx(0.D0,0.D0)
     qold = cmplx(tpi*(aimag(W(j))+aimag(W(k))),tpi*(real(W(j))-real(W(k))))
     if(real(qold) .lt. 0.0) qold = -qold
     do iz = 1,izz
       q = q + dble(lambdaE(iz)) * exp( - qold * dble(iz) )
     end do
     Ijk(j,k) = conjg(CGinv(j,gzero)) * q * CGinv(k,gzero)
  end do
 end do
Ijk = Ijk * dzt

! then the matrix multiplications to obtain Lgh
tmp3 = matmul(conjg(CGG),Ijk)
Lgh = matmul(tmp3,transpose(CGG))

end subroutine CalcLgh_

!--------------------------------------------------------------------------
recursive subroutine CalcLghdepth_(self,DMat,Lgh,thick,kn,nn,gzero,depthstep,lambdaE,izz)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcLghdepth_
  !! author: MDG
  !! version: 1.0
  !! date: 02/12/20
  !!
  !! compute the Lgh matrix for EBSD, ECCI, ECP, etc simulations

use mod_io

IMPLICIT NONE

class(gvectors_T),INTENT(INOUT)     :: self
integer(kind=irg),INTENT(IN)        :: nn
complex(kind=dbl),INTENT(IN)        :: DMat(nn,nn)
complex(kind=dbl),INTENT(OUT)       :: Lgh(nn,nn,izz)
real(kind=dbl),INTENT(IN)           :: thick
real(kind=dbl),INTENT(IN)           :: kn
integer(kind=irg),INTENT(IN)        :: gzero
real(kind=dbl),INTENT(IN)           :: depthstep
integer(kind=irg),INTENT(IN)        :: izz
real(kind=sgl),INTENT(IN)           :: lambdaE(izz)

type(IO_T)                          :: Message

integer                             :: i,j,k, iz
complex(kind=dbl)                   :: CGinv(nn,nn), Minp(nn,nn), tmp3(nn,nn,izz)

real(kind=dbl)                      :: tpi, dzt
complex(kind=dbl)                   :: Ijk(nn,nn,izz), q, getMIWORK, qold

integer(kind=irg)                   :: INFO, LDA, LDVR, LDVL, JPIV(nn), MILWORK
complex(kind=dbl)                   :: CGG(nn,nn), W(nn)
complex(kind=dbl),allocatable       :: MIWORK(:)

integer(kind=irg),parameter         :: LWMAX = 5000 
complex(kind=dbl)                   :: VL(nn,nn),  WORK(LWMAX)
real(kind=dbl)                      :: RWORK(2*nn)
character                           :: JOBVL, JOBVR
integer(kind=sgl)                   :: LWORK
integer(kind=sgl)                   :: ilaenv

! compute the eigenvalues and eigenvectors
! using the LAPACK CGEEV, CGETRF, and CGETRI routines
!
 Minp = DMat

! set some initial LAPACK variables
 LDA = nn
 LDVL = nn
 LDVR = nn
 INFO = 0

 ! first initialize the parameters for the LAPACK ZGEEV, CGETRF, and CGETRI routines
 JOBVL = 'N'   ! do not compute the left eigenvectors
 JOBVR = 'V'   ! do compute the right eigenvectors
 LWORK = -1 ! so that we can ask the routine for the actually needed value

! call the routine to determine the optimal workspace size
 LDA = nn
  call zgeev(JOBVL,JOBVR,nn,Minp,LDA,W,VL,LDVL,CGG,LDVR,WORK,LWORK,RWORK,INFO)
  LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

! then call the eigenvalue solver
  LDA = nn
  call zgeev(JOBVL,JOBVR,nn,Minp,LDA,W,VL,LDVL,CGG,LDVR,WORK,LWORK,RWORK,INFO)
  if (INFO.ne.0) call Message%printError('Error in CalcLgh: ','ZGEEV return not zero')

 CGinv = CGG

 LDA=nn
 call zgetrf(nn,nn,CGinv,LDA,JPIV,INFO)
 MILWORK = -1
 LDA=nn

!  call zgetri(nn,CGinv,LDA,JPIV,getMIWORK,MILWORK,INFO)
!  MILWORK =  INT(real(getMIWORK))
 MILWORK=ilaenv( 1, 'ZGETRI', ' ', nn, -1, -1, -1 )
 if (.not.allocated(MIWORK)) allocate(MIWORK(nn*MILWORK))
 MIWORK = cmplx(0.D0,0.D0)
 LDA=nn
 call zgetri(nn,CGinv,LDA,JPIV,MIWORK,MILWORK,INFO)
 deallocate(MIWORK)

! in all the time that we've used these routines, we haven't
! had a single problem with the matrix inversion, so we don't
! really need to do this test:
!
! if ((cabs(sum(matmul(CGG,CGinv)))-dble(nn)).gt.1.E-8) write (*,*) 'Error in matrix inversion; continuing'


! then compute the integrated intensity matrix
 W = W/cmplx(2.0*kn,0.0)

! recall that alpha(1:nn) = CGinv(1:nn,gzero)

! first the Ijk matrix (this is Winkelmann's B^{ij}(t) matrix)
! combined with numerical integration over [0, z0] interval,
! taking into account depth profiles from Monte Carlo simulations ...
! the depth profile lambdaE must be added to the absorption
! components of the Bloch wave eigenvalues.

 tpi = 2.D0*cPi*depthstep
 dzt = depthstep/thick
  do j=1,nn
   do k=1,nn
      q =  cmplx(0.D0,0.D0)
      qold = cmplx(tpi*(aimag(W(j))+aimag(W(k))),tpi*(real(W(j))-real(W(k))))
      if(real(qold) .lt. 0.0) qold = -qold
      do iz = 1,izz
        q = dble(lambdaE(iz)) * exp( - qold * dble(iz) )
        Ijk(j,k,iz) = conjg(CGinv(j,gzero)) * q * CGinv(k,gzero)
      end do
   end do
  end do
 ! finally, multiply by dz/z_0 for the discrete integration (so we 
 ! don't forget to do that in the calling program...) 
 Ijk = Ijk * dzt
 
 ! then the matrix multiplications to obtain Lgh(z) 
 do iz=1,izz
   tmp3(1:nn,1:nn,iz) = matmul(conjg(CGG),Ijk(1:nn,1:nn,iz)) 
 end do
 do iz=1,izz
   Lgh(1:nn,1:nn,iz) = matmul(tmp3(1:nn,1:nn,iz),transpose(CGG))
 end do

end subroutine CalcLghdepth_

! ! this may need to be moved to the kvectors_T class ...
! ! this routine is currently not used anywhere !

! !> @param cell unit cell pointer
! !> @param khead start of reflisttype linked list
! !> @param reflist reflection linked list
! !> @param Dyn dynamical scattering structure
! !> @param BetheParameter Bethe parameter structure
! !> @param numk number of wave vectors to consider
! !> @param nbeams total number of unique beams
! !--------------------------------------------------------------------------
! recursive subroutine Prune_ReflectionList_(self, cell, Diff, khead, Dyn, numk, nbeams)
!   !! author: MDG
!   !! version: 1.0
!   !! date: 02/02/20
!   !!
!   !! select from the reflection list those g-vectors that will be counted in an LACBED computation
!   !!
!   !! This routine basicaly repeats a section from the Compute_DynMat routine
!   !! without actually computing the dynamical matrix; it simply keeps track of all the
!   !! beams that are at one point or another regarded as strong or weak beams.
!   !! We'll use the famnum field in the rlp linked list to flag the strong reflections.
!   !! Linked list entries that are not used are instantly removed from the linked list.

! use mod_io
! use mod_crystallography
! use mod_kvectors
! use mod_diffraction

! IMPLICIT NONE

! class(gvectors_T), INTENT(INOUT)        :: self
! type(Cell_T),INTENT(INOUT)              :: cell
! type(kvectorlist),pointer               :: khead
! type(DynType),INTENT(INOUT)             :: Dyn
! integer(kind=irg),INTENT(IN)            :: numk
! integer(kind=irg),INTENT(OUT)           :: nbeams

! type(IO_T)                              :: Message
! integer(kind=irg)                       :: ik, ig, istrong, curfam(3), newfam(3)
! real(kind=sgl)                          :: sgp, lUg, cut1, cut2
! !integer(kind=irg),allocatable          :: strongreflections(:,:)
! type(kvectorlist),pointer               :: ktmp
! type(reflisttype),pointer               :: rltmpa, rltmpb

! ! reset the value of DynNbeams in case it was modified in a previous call
! cell%DynNbeams = cell%DynNbeamsLinked

! nbeams = 0

! ! reset the reflection linked list
!   rltmpa => cell%reflist%next

! ! pick the first reflection since that is the transmitted beam (only on the first time)
!   rltmpa%famnum = 1
!   nbeams = nbeams + 1

! ! loop over all reflections in the linked list
! !!!! this will all need to be changed with the new Bethe potential criteria ...
!   rltmpa => rltmpa%next
!   reflectionloop: do ig=2,cell%DynNbeamsLinked
!     lUg = abs(rltmpa%Ucg) * cell%mLambda
!     cut1 = BetheParameter%cutoff * lUg
!     cut2 = BetheParameter%weakcutoff * lUg

! ! loop over all the incident beam directions
!     ktmp => khead
! ! loop over all beam orientations, selecting them from the linked list
!     kvectorloop: do ik = 1,numk
! ! We compare |sg| with two multiples of lambda |Ug|
! !
! !  |sg|>cutoff lambda |Ug|   ->  don't count reflection
! !  cutoff lambda |Ug| > |sg| > weakcutoff lambda |Ug|  -> weak reflection
! !  weakcutoff lambda |Ug| > |sg|  -> strong reflection
! !
! !       sgp = abs(CalcsgHOLZ(float(rltmpa%hkl),sngl(ktmp%kt),sngl(cell%mLambda)))
! !       write (*,*) rltmpa%hkl,CalcsgHOLZ(float(rltmpa%hkl),sngl(ktmp%kt), &
! !                       sngl(cell%mLambda)),Calcsg(float(rltmpa%hkl),sngl(ktmp%k),DynFN)
!         sgp = abs(Calcsg(cell,float(rltmpa%hkl),sngl(ktmp%k),Dyn%FN))
! ! we have to deal separately with double diffraction reflections, since
! ! they have a zero potential coefficient !
!         if ( cell%dbdiff(rltmpa%hkl(1),rltmpa%hkl(2),rltmpa%hkl(3)) ) then  ! it is a double diffraction reflection
!           if (sgp.le.BetheParameter%sgcutoff) then
!             nbeams = nbeams + 1
!             rltmpa%famnum = 1
!             EXIT kvectorloop    ! this beam did contribute, so we no longer need to consider it
!           end if
!         else   ! it is not a double diffraction reflection
!           if (sgp.le.cut1) then  ! count this beam, whether it is weak or strong
!             nbeams = nbeams + 1
!             rltmpa%famnum = 1
!             EXIT kvectorloop    ! this beam did contribute, so we no longer need to consider it
!           end if
!         end if

! ! go to the next incident beam direction
!        if (ik.ne.numk) ktmp => ktmp%next
!      end do kvectorloop  ! ik loop

! ! go to the next beam in the list
!    rltmpa => rltmpa%next
!   end do reflectionloop

!   call Message(' Renumbering reflections', frm = "(A)")

! ! change the following with the new next2 pointer in the reflist type !!!

! ! ok, now that we have the list, we'll go through it again to set sequential numbers instead of 1's
! ! at the same time, we'll deallocate those entries that are no longer needed.
!   rltmpa => reflist%next
!   rltmpb => rltmpa
!   rltmpa => rltmpa%next ! we keep the first entry, always.
!   istrong = 1
!   reflectionloop2: do ig=2,cell%DynNbeamsLinked
!     if (rltmpa%famnum.eq.1) then
!         istrong = istrong + 1
!         rltmpa%famnum = istrong
!         rltmpa => rltmpa%next
!         rltmpb => rltmpb%next
!     else   ! remove this entry from the linked list
!         rltmpb%next => rltmpa%next
!         deallocate(rltmpa)
!         rltmpa => rltmpb%next
!     endif
! ! go to the next beam in the list
!   end do reflectionloop2

! ! reset the number of beams to the newly obtained number
!   cell%DynNbeamsLinked = nbeams
!   cell%DynNbeams = nbeams

! ! go through the entire list once again to correct the famhkl
! ! entries, which may be incorrect now; famhkl is supposed to be one of the
! ! reflections on the current list, but that might not be the case since
! ! famhkl was first initialized when there were additional reflections on
! ! the list... so we set famhkl to be the same as the first hkl in each family.
!   rltmpa => reflist%next%next  ! no need to check the first one
! reflectionloop3:  do while (associated(rltmpa))
!     curfam = rltmpa%famhkl
!     if (sum(abs(curfam-rltmpa%hkl)).ne.0) then
!       newfam = rltmpa%hkl
!       do while (sum(abs(rltmpa%famhkl-curfam)).eq.0)
!         rltmpa%famhkl = newfam
!         rltmpa => rltmpa%next
!         if ( .not.associated(rltmpa) ) EXIT reflectionloop3
!       end do
!     else
!       do while (sum(abs(rltmpa%famhkl-curfam)).eq.0)
!         rltmpa => rltmpa%next
!         if ( .not.associated(rltmpa) ) EXIT reflectionloop3
!       end do
!     end if
! ! go to the next beam in the list
!   end do reflectionloop3

! end subroutine Prune_ReflectionList

! ! this is an old and currently unused routine ...
! !--------------------------------------------------------------------------
! !
! ! SUBROUTINE: Compute_DynMat
! !
! !> @author Marc De Graef, Carnegie Mellon University
! !
! !> @brief compute the entire dynamical matrix, including HOLZ and Bethe potentials
! !
! !> @details This is a complicated routine because it forms the core of the dynamical
! !> scattering simulations.  The routine must be capable of setting up the dynamical
! !> matrix for systematic row and zone axis, with or without HOLZ reflections, and
! !> with Bethe potentials for the Bloch wave case.  The routine must also be able to
! !> decide which reflections are weak, and which are strong (again for the Bloch wave
! !> case, but potentially also for other cases? Further research needed...).
! !
! !> @param cell unit cell pointer
! !> @param reflist reflection list pointer
! !> @param Dyn dynamical scattering structure
! !> @param BetheParameter Bethe parameter structure
! !> @param calcmode string that describes the particular matrix mode
! !> @param kk incident wave vector
! !> @param kt tangential component of incident wave vector (encodes the Laue Center)
! !> @param IgnoreFoilNormal switch for foil normal inclusion in sg computation
! !> @param IncludeSecondOrder (optional) switch to include second order correction to Bethe potentials
! !
! !> @date   05/06/13 MDG 1.0 original
! !> @date   08/30/13 MDG 1.1 correction of effective excitation error
! !> @date   09/20/13 MDG 1.2 added second order Bethe potential correction switch
! !> @date   06/09/14 MDG 2.0 added cell, reflist and BetheParameter as arguments
! !--------------------------------------------------------------------------
! recursive subroutine Compute_DynMat(cell,reflist,Dyn,BetheParameter,calcmode,kk,kt,IgnoreFoilNormal,IncludeSecondOrder)
! !DEC$ ATTRIBUTES DLLEXPORT :: Compute_DynMat
!   !! author: MDG
!   !! version: 1.0
!   !! date: 02/02/20
!   !!
!   !!

! use error
! use constants
! use crystal
! use diffraction
! use io

! IMPLICIT NONE

! type(unitcell)                          :: cell
! type(reflisttype),pointer               :: reflist
! type(DynType),INTENT(INOUT)            :: Dyn
! !f2py intent(in,out) ::  Dyn
! type(BetheParameterType),INTENT(INOUT) :: BetheParameter
! !f2py intent(in,out) ::  BetheParameter
! character(*),INTENT(IN)         :: calcmode             !< computation mode
! real(kind=dbl),INTENT(IN)               :: kk(3),kt(3)          !< incident wave vector and tangential component
! logical,INTENT(IN)                      :: IgnoreFoilNormal     !< how to deal with the foil normal
! logical,INTENT(IN),OPTIONAL             :: IncludeSecondOrder   !< second order Bethe potential correction switch

! complex(kind=dbl)                       :: czero,pre, weaksum, ughp, uhph
! integer(kind=irg)                       :: istat,ir,ic,nn, iweak, istrong, iw, ig, ll(3), gh(3), nnn, nweak, io_int(1)
! real(kind=sgl)                          :: glen,exer,gg(3), kpg(3), gplen, sgp, lUg, cut1, cut2, io_real(3)
! real(kind=dbl)                          :: lsfour, weaksgsum
! logical                                 :: AddSecondOrder
! type(gnode)                           :: rlp
! type(reflisttype),pointer               :: rltmpa, rltmpb

! AddSecondOrder = .FALSE.
! if (present(IncludeSecondOrder)) AddSecondOrder = .TRUE.

! ! has the list of reflections been allocated ?
! if (.not.associated(reflist)) call FatalError('Compute_DynMat',' reflection list has not been allocated')

! ! if the dynamical matrix has already been allocated, deallocate it first
! ! this is partially so that no program will allocate DynMat itself; it must be done
! ! via this routine only.
! if (allocated(Dyn%DynMat)) deallocate(Dyn%DynMat)

! ! initialize some parameters
! czero = cmplx(0.0,0.0,dbl)      ! complex zero
! pre = cmplx(0.0,cPi,dbl)                ! i times pi

! if (calcmode.ne.'BLOCHBETHE') then

! ! allocate DynMat
!           allocate(Dyn%DynMat(cell%DynNbeams,cell%DynNbeams),stat=istat)
!           Dyn%DynMat = czero
! ! get the absorption coefficient
!           call CalcUcg(cell, rlp, (/0,0,0/),applyqgshift=.TRUE. )
!           Dyn%Upz = rlp%Vpmod

! ! are we supposed to fill the off-diagonal part ?
!          if ((calcmode.eq.'D-H-W').or.(calcmode.eq.'BLOCH')) then
!           rltmpa => cell%reflist%next    ! point to the front of the list
! ! ir is the row index
!           do ir=1,cell%DynNbeams
!            rltmpb => cell%reflist%next   ! point to the front of the list
! ! ic is the column index
!            do ic=1,cell%DynNbeams
!             if (ic.ne.ir) then  ! exclude the diagonal
! ! compute Fourier coefficient of electrostatic lattice potential
!              gh = rltmpa%hkl - rltmpb%hkl
!              if (calcmode.eq.'D-H-W') then
!               call CalcUcg(cell, rlp,gh,applyqgshift=.TRUE.)
!               Dyn%DynMat(ir,ic) = pre*rlp%qg
!              else
!               Dyn%DynMat(ir,ic) = cell%LUT(gh(1),gh(2),gh(3))
!              end if
!             end if
!             rltmpb => rltmpb%next  ! move to next column-entry
!            end do
!            rltmpa => rltmpa%next   ! move to next row-entry
!           end do
!          end if

! ! or the diagonal part ?
!          if ((calcmode.eq.'DIAGH').or.(calcmode.eq.'DIAGB')) then
!           rltmpa => reflist%next   ! point to the front of the list
! ! ir is the row index
!           do ir=1,cell%DynNbeams
!            glen = CalcLength(cell,float(rltmpa%hkl),'r')
!            if (glen.eq.0.0) then
!             Dyn%DynMat(ir,ir) = cmplx(0.0,Dyn%Upz,dbl)
!            else  ! compute the excitation error
!             exer = Calcsg(cell,float(rltmpa%hkl),sngl(kk),Dyn%FN)

!             rltmpa%sg = exer
!             if (calcmode.eq.'DIAGH') then  !
!              Dyn%DynMat(ir,ir) = cmplx(0.0,2.D0*cPi*exer,dbl)
!             else
!              Dyn%DynMat(ir,ir) = cmplx(2.D0*exer/cell%mLambda,Dyn%Upz,dbl)
!             end if
!            endif
!            rltmpa => rltmpa%next   ! move to next row-entry
!           end do
!          end if

! else  ! this is the Bloch wave + Bethe potentials initialization (originally implemented in the EBSD programs)

! ! we don't know yet how many strong reflections there are so we'll need to determine this first
! ! this number depends on some externally supplied parameters, which we will get from a namelist
! ! file (which should be read only once by the Set_Bethe_Parameters routine), or from default values
! ! if there is no namelist file in the folder.
!         if (BetheParameter%cutoff.eq.0.0) call Set_Bethe_Parameters(BetheParameter)


! ! reset the value of DynNbeams in case it was modified in a previous call
!         cell%DynNbeams = cell%DynNbeamsLinked

! ! precompute lambda^2/4
!         lsfour = cell%mLambda**2*0.25D0

! ! first, for the input beam direction, determine the excitation errors of
! ! all the reflections in the master list, and count the ones that are
! ! needed for the dynamical matrix (weak as well as strong)
!         if (.not.allocated(BetheParameter%weaklist)) allocate(BetheParameter%weaklist(cell%DynNbeams))
!         if (.not.allocated(BetheParameter%stronglist)) allocate(BetheParameter%stronglist(cell%DynNbeams))
!         if (.not.allocated(BetheParameter%reflistindex)) allocate(BetheParameter%reflistindex(cell%DynNbeams))
!         if (.not.allocated(BetheParameter%weakreflistindex)) allocate(BetheParameter%weakreflistindex(cell%DynNbeams))

!         BetheParameter%weaklist = 0
!         BetheParameter%stronglist = 0
!         BetheParameter%reflistindex = 0
!         BetheParameter%weakreflistindex = 0

!         rltmpa => cell%reflist%next

! ! deal with the transmitted beam first
!     nn = 1              ! nn counts all the scattered beams that satisfy the cutoff condition
!     nnn = 1             ! nnn counts only the strong beams
!     nweak = 0           ! counts only the weak beams
!     BetheParameter%stronglist(nn) = 1   ! make sure that the transmitted beam is always a strong beam ...
!     BetheParameter%weaklist(nn) = 0
!     BetheParameter%reflistindex(nn) = 1

!     rltmpa%sg = 0.D0
! ! write (*,*) 'DynNbeamsLinked = ',DynNbeamsLinked

! ! loop over all reflections in the linked list
!     rltmpa => rltmpa%next
!     reflectionloop: do ig=2,cell%DynNbeamsLinked
!       gg = float(rltmpa%hkl)                    ! this is the reciprocal lattice vector

! ! deal with the foil normal; if IgnoreFoilNormal is .TRUE., then assume it is parallel to the beam direction
!      if (IgnoreFoilNormal) then
! ! we're taking the foil normal to be parallel to the incident beam direction at each point of
! ! the standard stereographic triangle, so cos(alpha) = 1 always in eqn. 5.11 of EM
!         kpg = kk+gg                             ! k0 + g (vectors)
!         gplen = CalcLength(cell,kpg,'r')        ! |k0+g|
!         rltmpa%sg = (1.0/cell%mLambda**2 - gplen**2)*0.5/gplen
!      else
!         rltmpa%sg = Calcsg(cell,gg,sngl(kk),Dyn%FN)
! ! here we need to determine the components of the Laue Center w.r.t. the g1 and g2 vectors
! ! and then pass those on to the routine;
! !       rltmpa%sg = CalcsgHOLZ(gg,sngl(kt),sngl(cell%mLambda))
! ! write (*,*) gg, Calcsg(gg,sngl(kk),DynFN), CalcsgHOLZ(gg,sngl(kt),sngl(cell%mLambda))
!      end if

! ! use the reflection num entry to indicate whether or not this
! ! reflection should be used for the dynamical matrix
! ! We compare |sg| with two multiples of lambda |Ug|
! !
! !  |sg|>cutoff lambda |Ug|   ->  don't count reflection
! !  cutoff lambda |Ug| > |sg| > weakcutoff lambda |Ug|  -> weak reflection
! !  weakcutoff lambda |Ug| > |sg|  -> strong reflection
! !
!         sgp = abs(rltmpa%sg)
!         lUg = abs(rltmpa%Ucg) * cell%mLambda
!         cut1 = BetheParameter%cutoff * lUg
!         cut2 = BetheParameter%weakcutoff * lUg

! ! we have to deal separately with double diffraction reflections, since
! ! they have a zero potential coefficient !
!         if ( cell%dbdiff(rltmpa%hkl(1),rltmpa%hkl(2),rltmpa%hkl(3)) ) then  ! it is a double diffraction reflection
!           if (sgp.le.BetheParameter%sgcutoff) then
!                 nn = nn+1
!                 nnn = nnn+1
!                 BetheParameter%stronglist(ig) = 1
!                 BetheParameter%reflistindex(ig) = nnn
!           end if
!         else   ! it is not a double diffraction reflection
!           if (sgp.le.cut1) then  ! count this beam
!                 nn = nn+1
! ! is this a weak or a strong reflection (in terms of Bethe potentials)?
!                 if (sgp.le.cut2) then ! it's a strong beam
!                         nnn = nnn+1
!                         BetheParameter%stronglist(ig) = 1
!                         BetheParameter%reflistindex(ig) = nnn
!                 else ! it's a weak beam
!                         nweak = nweak+1
!                         BetheParameter%weaklist(ig) = 1
!                         BetheParameter%weakreflistindex(ig) = nweak
!                 end if
!           end if
!         end if
! ! go to the next beam in the list
!       rltmpa => rltmpa%next
!     end do reflectionloop

! ! if we don't have any beams in this list (unlikely, but possible if the cutoff and
! ! weakcutoff parameters have unreasonable values) then we abort the run
! ! and we report some numbers to the user
!          if (nn.eq.0) then
!            call Message(' no beams found for the following parameters:', frm = "(A)")
!            io_real(1:3) = kk(1:3)
!            call WriteValue(' wave vector = ', io_real,3)
!            io_int(1) = nn
!            call WriteValue('  -> number of beams = ', io_int, 1)
!            call Message( '   -> check cutoff and weakcutoff parameters for reasonableness', frm = "(A)")
!            call FatalError('Compute_DynMat','No beams in list')
!         end if

! ! next, we define nns to be the number of strong beams, and nnw the number of weak beams.
!          BetheParameter%nns = sum(BetheParameter%stronglist)
!          BetheParameter%nnw = sum(BetheParameter%weaklist)

! ! add nns to the weakreflistindex to offset it; this is used for plotting reflections on CBED patterns
!         do ig=2,cell%DynNbeamsLinked
!           if (BetheParameter%weakreflistindex(ig).ne.0) then
!             BetheParameter%weakreflistindex(ig) = BetheParameter%weakreflistindex(ig) + BetheParameter%nns
!           end if
!         end do

! ! We may want to keep track of the total and average numbers of strong and weak beams
!          BetheParameter%totweak = BetheParameter%totweak + BetheParameter%nnw
!          BetheParameter%totstrong = BetheParameter%totstrong + BetheParameter%nns
!          if (BetheParameter%nnw.lt.BetheParameter%minweak) BetheParameter%minweak=BetheParameter%nnw
!          if (BetheParameter%nnw.gt.BetheParameter%maxweak) BetheParameter%maxweak=BetheParameter%nnw
!          if (BetheParameter%nns.lt.BetheParameter%minstrong) BetheParameter%minstrong=BetheParameter%nns
!          if (BetheParameter%nns.gt.BetheParameter%maxstrong) BetheParameter%maxstrong=BetheParameter%nns

! ! allocate arrays for weak and strong beam information
!         if (allocated(BetheParameter%weakhkl)) deallocate(BetheParameter%weakhkl)
!         if (allocated(BetheParameter%weaksg)) deallocate(BetheParameter%weaksg)
!         if (allocated(BetheParameter%stronghkl)) deallocate(BetheParameter%stronghkl)
!         if (allocated(BetheParameter%strongsg)) deallocate(BetheParameter%strongsg)
!         if (allocated(BetheParameter%strongID)) deallocate(BetheParameter%strongID)
!         allocate(BetheParameter%weakhkl(3,BetheParameter%nnw),BetheParameter%weaksg(BetheParameter%nnw))
!         allocate(BetheParameter%stronghkl(3,BetheParameter%nns),BetheParameter%strongsg(BetheParameter%nns))
!         allocate(BetheParameter%strongID(BetheParameter%nns))

! ! here's where we extract the relevant information from the linked list (much faster
! ! than traversing the list each time...)
!         rltmpa => cell%reflist%next    ! reset the a list
!         iweak = 0
!         istrong = 0
!         do ir=1,cell%DynNbeamsLinked
!              if (BetheParameter%weaklist(ir).eq.1) then
!                 iweak = iweak+1
!                 BetheParameter%weakhkl(1:3,iweak) = rltmpa%hkl(1:3)
!                 BetheParameter%weaksg(iweak) = rltmpa%sg
!              end if
!              if (BetheParameter%stronglist(ir).eq.1) then
!                 istrong = istrong+1
!                 BetheParameter%stronghkl(1:3,istrong) = rltmpa%hkl(1:3)
!                 BetheParameter%strongsg(istrong) = rltmpa%sg
! ! make an inverse index list
!                 BetheParameter%strongID(istrong) = ir
!              end if
!            rltmpa => rltmpa%next
!         end do

! ! now we are ready to create the dynamical matrix
!         cell%DynNbeams = BetheParameter%nns

! ! allocate DynMat if it hasn't already been allocated and set to complex zero
!           if (allocated(Dyn%DynMat)) deallocate(Dyn%DynMat)
!           allocate(Dyn%DynMat(cell%DynNbeams,cell%DynNbeams),stat=istat)
!           Dyn%DynMat = czero

! ! get the absorption coefficient
!           call CalcUcg(cell, rlp, (/0,0,0/) )
!           Dyn%Upz = rlp%Vpmod

! ! ir is the row index
!        do ir=1,BetheParameter%nns
! ! ic is the column index
!           do ic=1,BetheParameter%nns
! ! compute the Bethe Fourier coefficient of the electrostatic lattice potential
!               if (ic.ne.ir) then  ! not a diagonal entry
!                  ll = BetheParameter%stronghkl(1:3,ir) - BetheParameter%stronghkl(1:3,ic)
!                  Dyn%DynMat(ir,ic) = cell%LUT(ll(1),ll(2),ll(3))
!         ! and subtract from this the total contribution of the weak beams
!          weaksum = czero
!          do iw=1,BetheParameter%nnw
!               ll = BetheParameter%stronghkl(1:3,ir) - BetheParameter%weakhkl(1:3,iw)
!               ughp = cell%LUT(ll(1),ll(2),ll(3))
!               ll = BetheParameter%weakhkl(1:3,iw) - BetheParameter%stronghkl(1:3,ic)
!               uhph = cell%LUT(ll(1),ll(2),ll(3))
!               weaksum = weaksum +  ughp * uhph *cmplx(1.D0/BetheParameter%weaksg(iw),0.0,dbl)
!          end do
!         ! and correct the dynamical matrix element to become a Bethe potential coefficient
!          Dyn%DynMat(ir,ic) = Dyn%DynMat(ir,ic) - cmplx(0.5D0*cell%mLambda,0.0D0,dbl)*weaksum
! ! do we need to add the second order corrections ?
!                   if (AddSecondOrder) then
!                     weaksum = czero
!                   end if
!               else  ! it is a diagonal entry, so we need the excitation error and the absorption length
! ! determine the total contribution of the weak beams
!                  weaksgsum = 0.D0
!                   do iw=1,BetheParameter%nnw
!                       ll = BetheParameter%stronghkl(1:3,ir) - BetheParameter%weakhkl(1:3,iw)
!                       ughp = cell%LUT(ll(1),ll(2),ll(3))
!                       weaksgsum = weaksgsum +  abs(ughp)**2/BetheParameter%weaksg(iw)
!                  end do
!                  weaksgsum = weaksgsum * cell%mLambda/2.D0
!                  Dyn%DynMat(ir,ir) = cmplx(2.D0*BetheParameter%strongsg(ir)/cell%mLambda-weaksgsum,Dyn%Upz,dbl)
! ! do we need to add the second order corrections ?
!                   if (AddSecondOrder) then
!                     weaksum = czero
!                   end if
!                end if
!           end do
!         end do
! ! that should do it for the initialization of the dynamical matrix

! end if   ! Bethe potential initialization

! end subroutine Compute_DynMat



end module mod_gvectors
