! ###################################################################
! Copyright (c) 2013-2020, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_LaueSupport
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! class definition for the Laue modality
  
  use mod_kinds
  use mod_global
  
  IMPLICIT NONE 
  private

  ! linked list for Laue XRD computations.
  type, public :: Laue_g_list  
    integer(kind=irg)   :: hkl(3)       ! Miller indices
    real(kind=dbl)      :: xyz(3)       ! Cartesian components of the plane normal
    real(kind=dbl)      :: tt           ! 2theta value
    real(kind=dbl)      :: polar        ! polarization factor
    real(kind=dbl)      :: sfs          ! |structure factor|^2
    type(Laue_g_list),pointer :: next   ! connection to next reflector 
  end type Laue_g_list

  ! class definition
  type, public :: Laue_T
  private 
    type(Laue_g_list), pointer :: reflist
    type(Laue_g_list), pointer :: rltail
    integer(kind=irg)          :: nref
  
  contains
  private 
    procedure, pass(self) :: MakeRefList_
    procedure, pass(self) :: get_ListHead_
    procedure, pass(self) :: Laue_Init_Reflist_
  
    generic, public :: MakeRefList => MakeRefList_
    generic, public :: get_ListHead => get_ListHead_
    generic, public :: Init_Reflist => Laue_Init_Reflist_
  
  end type Laue_T
  
  ! the constructor routine for this class 
  interface Laue_T
    module procedure Laue_constructor
  end interface Laue_T
  
  contains
  
  !--------------------------------------------------------------------------
  type(Laue_T) function Laue_constructor( ) result(GVec)
  !DEC$ ATTRIBUTES DLLEXPORT :: Laue_constructor
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! constructor for the Laue_T Class; reads the name list 
   
  IMPLICIT NONE
  
  ! the calling program must make sure that the reflist is empty ...
  ! initialize the reflist
  nullify(GVec%reflist)
  GVec%nref = 0
  call GVec%MakeRefList()

  end function Laue_constructor
  
  !--------------------------------------------------------------------------
  subroutine Laue_destructor(self) 
  !DEC$ ATTRIBUTES DLLEXPORT :: Laue_destructor
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! destructor for the Laue_T Class
   
  IMPLICIT NONE
  
  type(Laue_T), INTENT(INOUT)  :: self 
  
  call reportDestructor('Laue_T')
  
  end subroutine Laue_destructor

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

  class(Laue_T), INTENT(INOUT)  :: self

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
  recursive function get_ListHead_(self) result(glist)
  !DEC$ ATTRIBUTES DLLEXPORT :: get_ListHead_
  !! author: MDG
  !! version: 1.0
  !! date: 02/12/20
  !!
  !! return the number of k-vectors

  IMPLICIT NONE

  class(Laue_T), INTENT(INOUT)  :: self
  type(Laue_g_list), pointer    :: glist

  glist => self%reflist

  end function get_ListHead_
  
  !--------------------------------------------------------------------------
  subroutine Laue_Init_Reflist_(self, cell, SG, Diff, gcnt, lambdamin, intfactor, verbose)
  !DEC$ ATTRIBUTES DLLEXPORT :: Laue_Init_Reflist_
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

  IMPLICIT NONE 
  
  class(Laue_T), INTENT(INOUT)       :: self
  type(IO_T)                         :: Message
  type(SpaceGroup_T), INTENT(INOUT)  :: SG
  type(Diffraction_T), INTENT(INOUT) :: Diff
  type(Cell_T), INTENT(INOUT)        :: cell
  integer(kind=irg),INTENT(OUT)      :: gcnt
  real(kind=sgl), INTENT(IN)         :: lambdamin
  real(kind=dbl) , INTENT(IN)        :: intfactor
  logical,INTENT(IN),OPTIONAL        :: verbose

  type(Laue_g_list),pointer         :: gtmp, gtail                ! linked list for allowed g-vector search 
  type(gnode)                       :: rlp
  real(kind=sgl)                    :: gmax                       !< diameter of limiting sphere
  real(kind=sgl)                    :: ghkl                       !< length of a reciprocal lattice vector
  integer(kind=irg)                 :: imh, imk, iml              !< maximum index along a*, b*, and c*
  real(kind=sgl)                    :: g(3), tt                   !< g-vector and 2theta
  integer(kind=irg)                 :: io_int(3)                  !< io variable
  real(kind=sgl)                    :: io_real(1)                 !< io variable
  integer(kind=irg)                 :: istat, h, k, l, icnt       !< status variables and such
  !real(kind=sgl),parameter          :: tdtr = 114.5915590262      !< 2 * 180.0 / pi
  real(kind=sgl)                    :: threshold, th, sfs         !< threshold for discarding allowed reflections, |F|^2


  ! first get the range of Miller indices based on the lattice parameters and the xray wave length
  gmax = 2.0 / lambdamin      ! radius of the limiting sphere for smallest wave length  [nm^-1]
  imh = 1
  do   ! a* direction
    imh = imh + 1
    ghkl = cell%CalcLength((/float(imh) ,0.0_sgl,0.0_sgl/), 'r')
    if (ghkl.gt.gmax) EXIT
  end do
  imk = 1
  do  ! b* direction
    imk = imk + 1
    ghkl = cell%CalcLength((/0.0_sgl,float(imk),0.0_sgl/), 'r')
    if (ghkl.gt.gmax) EXIT
  end do
  iml = 1
  do  ! c* direction
    iml = iml + 1
    ghkl = cell%CalcLength((/0.0_sgl,0.0_sgl,float(iml)/), 'r')
    if (ghkl.gt.gmax) EXIT
  end do

  ! output range
  if (present(verbose)) then 
    if (verbose) then
      io_int = (/ imh, imk ,iml /)
      call Message%WriteValue('Range of reflections along a*, b* and c* = ', io_int, 3)
    end if
   end if 

  ! next we make a list of all rlp's that satisfy the following conditions
  !  - rlp must be inside the limiting sphere;
  !  - rlp must not be a systematic extinction;
  !  - rlp must not be a symmetry-induced extinction (since everything is kinematical)
  ! Since we don't know a-priori how many rlps will satisfy all
  ! three conditions, we'll first create a linked list and then copy
  ! that list into an allocatable array reflist
  if (.not.associated(self%reflist)) then ! allocate the head and tail of the linked list
   allocate(self%reflist,stat=istat)         ! allocate new value
   if (istat.ne.0) call Message%printError('Laue_Init_Reflist', 'unable to allocate reflist pointer')
   gtail => self%reflist                     ! tail points to new value
   nullify(gtail%next)                  ! nullify next in new value
  end if
  gtail => self%reflist
  call Diff%setrlpmethod('XR')

  ! compute the intensity threshold parameter as a fraction of |F(000)|^2 
  call Diff%CalcUcg(cell,(/0, 0, 0/))
  rlp = Diff%getrlp()
  threshold = cabs(rlp%Ucg)**2
  if (present(verbose)) then
    if (verbose) then
      io_real(1) = threshold 
      call Message%WriteValue(' Intensity Threshold value : ', io_real, 1)
    end if
  end if
  
  ! now loop over all g-vectors inside the imh, imk, iml range
  gcnt = 0
  icnt = 0
  th = sngl(intfactor) * threshold
  do h=-imh,imh
    do k=-imk,imk
      do l=-iml,iml
        icnt = icnt + 1
        g =float( (/ h, k, l /) )
! first of all, is this reflection inside the limiting sphere? (CalcLength)
        ghkl = cell%CalcLength(g,'r')
        if ((ghkl.le.gmax).and.(ghkl.gt.0.0)) then 
! then see if the reflection is allowed by systematic extinction (IsGAllowed)
          if ( SG%IsGAllowed( (/ h,k,l /) ) ) then    ! this is not a systematic extinction
! does this reflection have a non-zero structure factor larger than the threshold?
            call Diff%CalcUcg(cell, (/ h, k, l /) )
            sfs = cabs(rlp%Ucg)**2 
            if (sfs.ge.th) then   ! count this reflection 
              gcnt = gcnt + 1
! fill in the values
              gtail%hkl = (/ h, k, l /)
              write(*,*) gtail%hkl
              call cell%TransSpace(dble(gtail%hkl), gtail%xyz, 'r', 'c')
!            call NormVec(cell, gtail%xyz, 'c')    ! removed by MDG, 07/30/19 for EMLaue program
              gtail % tt = Diff%CalcDiffAngle(cell, (/h, k, l/))
              gtail % polar = (1.D0+ cos(2.D0*gtail%tt)**2)*0.5D0
              gtail % sfs = sfs / threshold
! and add it to the linked list
              allocate(gtail%next,stat=istat)    ! allocate new value
              if (istat.ne.0) call Message%printError('Laue_Init_Reflist', 'unable to allocate new entry in ghead linked list')
              gtail => gtail%next              ! gtail points to new value
              nullify(gtail%next)              ! nullify next in new value
            end if
          end if
        end if
      end do
    end do
  end do

  if (present(verbose)) then
    if (verbose) then
      io_int(1) = gcnt
      call Message%WriteValue(' Total number of reflections = ', io_int, 1)
    end if
  end if 


  end subroutine Laue_Init_Reflist_
  

  end module mod_LaueSupport

! these routines need to go in the Laue_T Class
! !--------------------------------------------------------------------------
! !
! ! SUBROUTINE: Laue_Init_Reflist
! !
! !> @author Marc De Graef, Carnegie Melon University
! !
! !> @brief compute the list of all possible reciprocal lattice points for Laue XRD
! !
! !> @param verbose print output when .TRUE.
! !
! !> @date 03/14/19 MDG 1.0 original
! !--------------------------------------------------------------------------
! subroutine Laue_Init_Reflist(cell, lmnl, reflist, gcnt, verbose)
!   !! author: MDG 
!   !! version: 1.0 
!   !! date: 02/02/20
!   !!
!   !!

! use local
! use io
! use crystal
! use error
! use symmetry
! use diffraction 
! use NameListTypedefs

! IMPLICIT NONE

! type(unitcell)                    :: cell
! type(LaueMasterNameListType),INTENT(INOUT) :: lmnl
! !f2py intent(in,out) ::  lmnl
! type(Laue_g_list),pointer         :: reflist                    ! linked list for allowed g-vector search 
! integer(kind=irg),INTENT(OUT)     :: gcnt
! logical,OPTIONAL,INTENT(IN)       :: verbose                    ! print output or not ?

! type(Laue_g_list),pointer         :: gtmp, gtail                ! linked list for allowed g-vector search 
! type(gnode)                       :: rlp

! real(kind=sgl)                    :: gmax                       !< diameter of limiting sphere
! real(kind=sgl)                    :: ghkl                       !< length of a reciprocal lattice vector
! integer(kind=irg)                 :: imh, imk, iml              !< maximum index along a*, b*, and c*
! real(kind=sgl)                    :: g(3), tt                   !< g-vector and 2theta
! integer(kind=irg)                 :: io_int(3)                  !< io variable
! real(kind=sgl)                    :: io_real(1)                 !< io variable
! integer(kind=irg)                 :: istat, h, k, l, icnt       !< status variables and such
! !real(kind=sgl),parameter          :: tdtr = 114.5915590262      !< 2 * 180.0 / pi
! real(kind=sgl)                    :: threshold, th, sfs         !< threshold for discarding allowed reflections, |F|^2

! ! first get the range of Miller indices based on the lattice parameters and the xray wave length
!  gmax = 2.0 / lmnl%lambdamin      ! radius of the limiting sphere for smallest wave length  [nm^-1]
!  imh = 1
!  do   ! a* direction
!    imh = imh + 1
!    ghkl = CalcLength(cell,  (/float(imh) ,0.0_sgl,0.0_sgl/), 'r')
!    if (ghkl.gt.gmax) EXIT
!  end do
!  imk = 1
!  do  ! b* direction
!    imk = imk + 1
!    ghkl = CalcLength(cell, (/0.0_sgl,float(imk),0.0_sgl/), 'r')
!    if (ghkl.gt.gmax) EXIT
!  end do
!  iml = 1
!  do  ! c* direction
!    iml = iml + 1
!    ghkl = CalcLength(cell, (/0.0_sgl,0.0_sgl,float(iml)/), 'r')
!    if (ghkl.gt.gmax) EXIT
!  end do

! ! output range
! if (present(verbose)) then 
!   if (verbose) then
!     io_int = (/ imh, imk ,iml /)
!     call WriteValue('Range of reflections along a*, b* and c* = ', io_int, 3)
!   end if
! end if 
 
! ! next we make a list of all rlp's that satisfy the following conditions
! !  - rlp must be inside the limiting sphere;
! !  - rlp must not be a systematic extinction;
! !  - rlp must not be a symmetry-induced extinction (since everything is kinematical)
! ! Since we don't know a-priori how many rlps will satisfy all
! ! three conditions, we'll first create a linked list and then copy
! ! that list into an allocatable array reflist
!  if (.not.associated(reflist)) then ! allocate the head and tail of the linked list
!    allocate(reflist,stat=istat)         ! allocate new value
!    if (istat.ne.0) call FatalError('Laue_Init_Reflist', 'unable to allocate reflist pointer')
!    gtail => reflist                     ! tail points to new value
!    nullify(gtail%next)                  ! nullify next in new value
!  end if

! ! initialize the computation mode X-Ray
! rlp%method = 'XR'

! ! compute the intensity threshold parameter as a fraction of |F(000)|^2 
!  call CalcUcg(cell, rlp, (/0, 0, 0/) )
!  threshold = cabs(rlp%Ucg)**2
! if (present(verbose)) then
!   if (verbose) then
!     io_real(1) = threshold 
!     call WriteValue(' Intensity Threshold value : ', io_real, 1)
!   end if
! end if

! ! now loop over all g-vectors inside the imh, imk, iml range
! gcnt = 0
! icnt = 0
! th = lmnl%intfactor * threshold
! do h=-imh,imh
!  do k=-imk,imk
!   do l=-iml,iml
!    icnt = icnt + 1
!     g =float( (/ h, k, l /) )
! ! first of all, is this reflection inside the limiting sphere? (CalcLength)
!     ghkl = CalcLength(cell,g,'r')
!     if ((ghkl.le.gmax).and.(ghkl.gt.0.0)) then 
! ! then see if the reflection is allowed by systematic extinction (IsGAllowed)
!        if ( IsGAllowed(cell, (/ h,k,l /) ) ) then    ! this is not a systematic extinction
! ! does this reflection have a non-zero structure factor larger than the threshold?
!            call CalcUcg(cell, rlp, (/ h, k, l /) )
!            sfs = cabs(rlp % Ucg)**2 
!            if (sfs.ge.th) then   ! count this reflection 
!              gcnt = gcnt + 1
! ! fill in the values
!              gtail % hkl = (/ h, k, l /)
!              call TransSpace(cell, dble(gtail % hkl), gtail % xyz, 'r', 'c')
! !             call NormVec(cell, gtail%xyz, 'c')    ! removed by MDG, 07/30/19 for EMLaue program
!              gtail % tt = CalcDiffAngle(cell,h,k,l)
!              gtail % polar = (1.D0+ cos(2.D0*gtail%tt)**2)*0.5D0
!              gtail % sfs = sfs / threshold
! ! and add it to the linked list
!              allocate(gtail%next,stat=istat)    ! allocate new value
!              if (istat.ne.0) call FatalError('Laue_Init_Reflist', 'unable to allocate new entry in ghead linked list')
!              gtail => gtail%next              ! gtail points to new value
!              nullify(gtail%next)              ! nullify next in new value
!            end if
!        end if
!     end if
!   end do
!  end do
! end do

! if (present(verbose)) then
!   if (verbose) then
!     io_int(1) = gcnt
!     call WriteValue(' Total number of reflections = ', io_int, 1)
!   end if
! end if 

! end subroutine Laue_Init_Reflist



! !--------------------------------------------------------------------------
! !
! ! SUBROUTINE: Laue_Init_Unit_Reflist
! !
! !> @author Marc De Graef, Carnegie Melon University
! !
! !> @brief compute the list of all possible reciprocal lattice points for Laue XRD;
! !> in this particular version, we return only unit length g-vectors, and then only
! !> one Friedel pair for each lattice row.  We also keep the structure factors for
! !> entire row, since we will need to evaluate which multiple of g actually causes the
! !> reflection; we will also need the d-spacings for each of them.  The linked list 
! !> generated by this routine has a type that is different from the regular Laue linked list
! !
! !> @param verbose print output when .TRUE.
! !
! !> @date 01/29/20 MDG 1.0 original
! !--------------------------------------------------------------------------
! subroutine Laue_Init_Unit_Reflist(cell, lmnl, reflist, gcnt, verbose)
! !DEC$ ATTRIBUTES DLLEXPORT :: Laue_Init_Unit_Reflist

! use local
! use io
! use crystal
! use error
! use symmetry
! use diffraction 
! use postscript
! use NameListTypedefs

! IMPLICIT NONE

! type(unitcell)                    :: cell
! type(LaueMasterNameListType),INTENT(INOUT) :: lmnl
! !f2py intent(in,out) ::  lmnl
! type(Laue_grow_list),pointer      :: reflist                    ! linked list for allowed g-vector search 
! integer(kind=irg),INTENT(OUT)     :: gcnt
! logical,OPTIONAL,INTENT(IN)       :: verbose                    ! print output or not ?

! type(Laue_grow_list),pointer      :: gtmp, gtail                ! linked list for allowed g-vector search 
! type(gnode)                       :: rlp

! logical,allocatable               :: z(:,:,:)

! real(kind=sgl)                    :: gmax                       !< diameter of limiting sphere
! real(kind=sgl)                    :: ghkl                       !< length of a reciprocal lattice vector
! integer(kind=irg)                 :: imh, imk, iml              !< maximum index along a*, b*, and c*
! real(kind=sgl)                    :: tt                         !< 2theta
! integer(kind=irg)                 :: io_int(3)                  !< io variable
! real(kind=sgl)                    :: io_real(1)                 !< io variable
! integer(kind=irg)                 :: i, istat, h, k, l, icnt, g(3), gr(3), rf, lcnt       !< status variables and such
! !real(kind=sgl),parameter          :: tdtr = 114.5915590262      !< 2 * 180.0 / pi
! real(kind=sgl)                    :: threshold, th, sfs         !< threshold for discarding allowed reflections, |F|^2

! ! first get the range of Miller indices based on the lattice parameters and the xray wave length
!  gmax = 2.0 / lmnl%lambdamin      ! radius of the limiting sphere for smallest wave length  [nm^-1]
!  imh = 1
!  do   ! a* direction
!    imh = imh + 1
!    ghkl = CalcLength(cell,  (/float(imh) ,0.0_sgl,0.0_sgl/), 'r')
!    if (ghkl.gt.gmax) EXIT
!  end do
!  imk = 1
!  do  ! b* direction
!    imk = imk + 1
!    ghkl = CalcLength(cell, (/0.0_sgl,float(imk),0.0_sgl/), 'r')
!    if (ghkl.gt.gmax) EXIT
!  end do
!  iml = 1
!  do  ! c* direction
!    iml = iml + 1
!    ghkl = CalcLength(cell, (/0.0_sgl,0.0_sgl,float(iml)/), 'r')
!    if (ghkl.gt.gmax) EXIT
!  end do

! ! output range
! if (present(verbose)) then 
!   if (verbose) then
!     io_int = (/ imh, imk ,iml /)
!     call WriteValue('Range of reflections along a*, b* and c* = ', io_int, 3)
!   end if
! end if 
 
! ! logical array to keep track of reflections that we have already dealt with (.TRUE.)
! allocate(z(-imh:imh, -imk:imk, -iml:iml))
! z = .FALSE.   ! all are set to .FALSE. to start the triple loop through the reciprocal lattice 

! ! next we make a list of all rlp's that satisfy the following conditions
! !  - rlp must be inside the limiting sphere;
! !  - rlp must not be a systematic extinction;
! !  - rlp must not be a symmetry-induced extinction (since everything is kinematical)
! ! Since we don't know a-priori how many rlps will satisfy all
! ! three conditions, we'll first create a linked list and then copy
! ! that list into an allocatable array reflist
!  if (.not.associated(reflist)) then ! allocate the head and tail of the linked list
!    allocate(reflist,stat=istat)         ! allocate new value
!    if (istat.ne.0) call FatalError('Laue_Init_Reflist', 'unable to allocate reflist pointer')
!    gtail => reflist                     ! tail points to new value
!    nullify(gtail%next)                  ! nullify next in new value
!  end if

! ! initialize the computation mode X-Ray
! rlp%method = 'XR'

! ! compute the intensity threshold parameter as a fraction of |F(000)|^2 
!  call CalcUcg(cell, rlp, (/0, 0, 0/) )
!  threshold = cabs(rlp%Ucg)**2
! if (present(verbose)) then
!   if (verbose) then
!     io_real(1) = threshold 
!     call WriteValue(' Intensity Threshold value : ', io_real, 1)
!   end if
! end if

! ! now loop over all g-vectors inside the imh, imk, iml range
! ! we keep only the unit normal, and we also keep the structure factors 
! ! and d-spacings for a series of positive multiples of g (i.e., we keep Friedel
! ! pairs apart).  This will then allow for an efficient scan through the Ewald
! ! volume in the pattern generation module. 

! ! to ensure that we get the correct ranges for the sfs and dspacing arrays, we 
! ! first reduce each hkl to the smallest common denominator, which sets the Nentries
! ! parameter; then we compute all the sfs and dspacing values, as well as the unit 
! ! g-vector (in the Cartesian crystal reference frame) and we set all the points along
! ! the g row to .TRUE. in the z logical array
! gcnt = 0
! icnt = 0
! lcnt = 0
! th = lmnl%intfactor * threshold

! ! loop over the entire reciprocal space sublattice
! do h=-imh,imh
!  do k=-imk,imk
!   do l=-iml,iml
! ! skip the origin !  (transmitted beam will be handled separately)
!     if (maxval(abs( (/ h, k, l /) ) ).eq.0) CYCLE
! ! have we already dealt with this reflection?
!     if (z(h,k,l).eqv..TRUE.) CYCLE 
! ! is this reflection forbidden by lattice centering ?
!     if ( .not.IsGAllowed(cell, (/ h, k, l /) ) ) then 
!       z(h,k,l) = .TRUE.
!       CYCLE
!     end if 
! ! we haven't covered this one yet so let's reduce the Miller indices to the lowest common denominator
!     g = (/ h, k, l /)
!     gr = g
!     call IndexReduce( gr )
! ! the reduction factor is ...
!     do i=1,3
!       if (gr(i).ne.0) then 
!         rf = g(i)/gr(i)
!         EXIT 
!       end if 
!     end do 
! ! have we done this one yet ?  If not, then we fill in the linked list entry, and 
! ! allocate the next one 
!     if (z(gr(1),gr(2),gr(3)).eqv..FALSE.) then
! ! set the entire systematic row to .TRUE. to avoid visiting them again 
!       do rf=1,100
!         ghkl = CalcLength(cell,float(rf*gr),'r')
!         if (ghkl.gt.gmax) EXIT 
!       end do
!       rf = rf-1
!       do i=1,rf 
!         z(i*gr(1), i*gr(2), i*gr(3)) = .TRUE.
!       end do
! ! the reduction factor is also the Nentries parameter for the linked list, so we create a 
! ! new entry in the list and generate the proper arrays sfs and dspacing
!       gtail % hkl = gr 
! ! convert the shortest g-vector to a unit cartesian vector
!       call TransSpace(cell, dble(gr), gtail % xyz, 'r', 'c')
!       call NormVec(cell, gtail%xyz, 'c')
! ! then deal with the intensities and d-spacings
!       gtail % Nentries = rf 
!       allocate( gtail%sfs(rf), gtail%dspacing(rf) )
!       gtail % sfs = 0.0
!       gtail % dspacing = 0.0
!       do i=1,rf
! ! is this reflection inside the limiting sphere? (CalcLength)
!         ghkl = CalcLength(cell,float(i*gr),'r')
!         if ((ghkl.le.gmax).and.(ghkl.gt.0.0)) then 
!           if ( IsGAllowed(cell, i*gr ) ) then ! allowed reflection, so compute the entries
!             call CalcUcg(cell, rlp, i*gr )
!             gtail % sfs(i) = cabs(rlp % Ucg)**2 / threshold 
!             if (gtail%sfs(i).gt.lmnl%intfactor) then 
!               gtail % dspacing(i) = 1.0/CalcLength(cell, float(i*gr), 'r')
!             else 
!               gtail % dspacing(i) = 0.0
!               gtail % sfs(i) = 0.0
!             end if 
!           end if
!         end if
!       end do
!       if (sum(gtail%sfs).eq.0.0) then 
!         deallocate(gtail%sfs, gtail%dspacing)
!       else
! ! extend the linked list
!         allocate(gtail%next,stat=istat)    ! allocate new value
!         if (istat.ne.0) call FatalError('Laue_Init_Unit_Reflist', 'unable to allocate new entry in linked list')
!         gtail => gtail%next              ! gtail points to new value
!         nullify(gtail%next)              ! nullify next in new value
!         gcnt = gcnt + 1
!       end if
!     end if
!     icnt = icnt + 1
!   end do
!  end do
! end do

! if (present(verbose)) then
!   if (verbose) then
!     io_int(1:2) = (/ gcnt, icnt /)
!     call WriteValue(' Total number of reflections accepted/tested = ', io_int, 2, frm="(I8,'/',I8)")
!   end if
! end if 

! end subroutine Laue_Init_Unit_Reflist
