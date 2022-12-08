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

module mod_kvectors
  !! author: MDG
  !! version: 1.0
  !! date: 02/02/20
  !!
  !! variables and types needed to determine lists of wave vectors
 
 use mod_kinds
 use mod_global
 
 IMPLICIT NONE
 
 private
 
 ! linked list of wave vectors (used by all diffraction programs)
 type, public :: kvectorlist
   integer(kind=irg)             :: i,j,hs       ! image coordinates
   real(kind=dbl)                :: kt(3)        ! tangential component of wavevector
   real(kind=dbl)                :: kn           ! normal component
   real(kind=dbl)                :: k(3)         ! full wave vector
   type(kvectorlist),pointer     :: next         ! connection to next wave vector
 end type kvectorlist
 
 type, public :: kvectors_T
   private
     type(kvectorlist), pointer  :: klist
     real(kind=dbl)              :: kinp(3)
     real(kind=dbl)              :: ktmax
     integer(kind=irg)           :: numk
     integer(kind=irg)           :: isym
     character(fnlen)            :: mapmode
     real(kind=dbl)              :: delta
     real(kind=dbl)              :: gan(3)
     real(kind=dbl)              :: gperp(3)
     real(kind=dbl)              :: kstar(3)
 
   contains
   private
 
    procedure, pass(self) :: MakeRefList_
    procedure, pass(self) :: Calckvectors_
    procedure, pass(self) :: CalckvectorsSymmetry_
    procedure, pass(self) :: get_ListHead_
    procedure, pass(self) :: check_mapmode_
    procedure, pass(self) :: get_numk_
    procedure, pass(self) :: get_mapmode_
    procedure, pass(self) :: set_kinp_
    procedure, pass(self) :: set_ktmax_
    procedure, pass(self) :: set_SamplingType_
    procedure, pass(self) :: set_mapmode_
    procedure, pass(self) :: Add_knode_
    procedure, pass(self) :: AddkVector_
    procedure, pass(self) :: Delete_kvectorlist_
    procedure, pass(self) :: CalckvectorsECP_
    procedure, pass(self) :: CalckvectorsPrecession_
    procedure, pass(self) :: CalckvectorsGPU_
    procedure, pass(self) :: Calckvectorcircle_
    procedure, pass(self) :: Calckvectorcone_
    procedure, pass(self) :: Calckvectortrace_

    final :: kvectors_destructor
 
    generic, public :: MakeRefList => MakeRefList_
    generic, public :: Calckvectors => Calckvectors_
    generic, public :: CalckvectorsSymmetry => CalckvectorsSymmetry_
    generic, public :: get_ListHead => get_ListHead_
    generic, public :: get_numk => get_numk_
    generic, public :: set_kinp => set_kinp_
    generic, public :: set_ktmax => set_ktmax_
    generic, public :: set_SamplingType => set_SamplingType_
    generic, public :: check_mapmode => check_mapmode_
    generic, public :: set_mapmode => set_mapmode_
    generic, public :: get_mapmode => get_mapmode_
    generic, public :: Add_knode => Add_knode_
    generic, public :: AddkVector => AddkVector_
    generic, public :: Delete_kvectorlist => Delete_kvectorlist_
    generic, public :: CalckvectorsECP => CalckvectorsECP_
    generic, public :: CalckvectorsPrecession => CalckvectorsPrecession_
    generic, public :: CalckvectorsGPU => CalckvectorsGPU_
    generic, public :: Calckvectorcircle => Calckvectorcircle_
    generic, public :: Calckvectorcone => Calckvectorcone_
    generic, public :: Calckvectortrace => Calckvectortrace_
 
 end type kvectors_T
 
 ! the constructor routine for this class
 interface kvectors_T
   module procedure kvectors_constructor
 end interface kvectors_T
 
 contains
 
 !--------------------------------------------------------------------------
 type(kvectors_T) function kvectors_constructor( ) result(KVec)
 !DEC$ ATTRIBUTES DLLEXPORT :: kvectors_constructor
 !! author: MDG
 !! version: 1.0
 !! date: 02/02/20
 !!
 !! constructor for the kvectors_T Class
 
 IMPLICIT NONE
 
 integer(kind=irg)     :: nref
 
 ! simply initialize the reflist; nref will be 0 but is not needed in calling program
 
 nullify(KVec%klist)
 KVec%numk = 0
 
 call KVec%MakeRefList(nref)
 
 end function kvectors_constructor
 
 !--------------------------------------------------------------------------
 subroutine kvectors_destructor(self)
 !DEC$ ATTRIBUTES DLLEXPORT :: kvectors_destructor
 !! author: MDG
 !! version: 1.0
 !! date: 02/02/20
 !!
 !! destructor for the kvectors_T Class
 
 IMPLICIT NONE
 
 type(kvectors_T), INTENT(INOUT)  :: self
 
 call reportDestructor('kvectors_T')
 ! 02-06-2020 CLément Lafond : cause fortran error 157 on Windows, commented for the moment has it not modify program behaviour
 !call self%Delete_kvectorlist()
 
 end subroutine kvectors_destructor
 
 !--------------------------------------------------------------------------
 recursive subroutine Delete_kvectorlist_(self)
 !DEC$ ATTRIBUTES DLLEXPORT :: Delete_kvectorlist_
 !! author: MDG
 !! version: 1.0
 !! date: 02/02/20
 !!
 !! delete the entire linked list
 
 IMPLICIT NONE
 
 class(kvectors_T), INTENT(INOUT) :: self
 
 type(kvectorlist),pointer        :: ktmp, ktail
 
 ! deallocate the entire linked list before returning, to prevent memory leaks
 if (associated(self%klist)) then
   ktail => self%klist
   if (associated(ktail%next)) then
     ktmp => ktail % next
     do
       if (associated(ktail)) deallocate(ktail)
       if (.not. associated(ktmp)) EXIT
       ktail => ktmp
       ktmp => ktail % next
     end do
   end if
 end if
 
 nullify(self%klist)
 self%numk = 0
 end subroutine Delete_kvectorlist_
 
 !--------------------------------------------------------------------------
 recursive subroutine MakeRefList_(self, nref)
 !DEC$ ATTRIBUTES DLLEXPORT :: MakeRefList_
   !! author: MDG
   !! version: 1.0
   !! date: 02/04/20
   !!
   !! allocate and initialize the linked reflection list
 
 use mod_io
 
 IMPLICIT NONE
 
 class(kvectors_T), INTENT(INOUT)  :: self
 integer(kind=irg),INTENT(INOUT)   :: nref
 
 type(IO_T)                        :: Message
 type(kvectorlist),pointer         :: rltail
 integer(kind=irg)                 :: istat
 
 if (associated(self%klist)) then
   call self%Delete_kvectorlist()
 end if
 
 ! create it if it does not already exist
 if (.not.associated(self%klist)) then
   nref = 0
   allocate(self%klist,stat=istat)
   if (istat.ne.0) call Message%printError('MakeRefList:',' unable to allocate pointer')
   rltail => self%klist           ! tail points to new value
   nullify(rltail%next)           ! nullify next in new value
 end if
 
 end subroutine MakeRefList_
 
 !--------------------------------------------------------------------------
 recursive function Kdelta(i,j) result(res)
 !DEC$ ATTRIBUTES DLLEXPORT :: Kdelta
 !! author: MDG
 !! version: 1.0
 !! date: 02/02/20
 !!
 !! Kronecker delta function, returns 1 or 0
 
 IMPLICIT NONE
 
 integer(kind=irg),INTENT(IN)    :: i,j
 integer(kind=irg)               :: res
 
  if (i.eq.j) then
    res = 1
  else
    res = 0
  end if
 
 end function Kdelta
 
 !--------------------------------------------------------------------------
 recursive function check_mapmode_(self, mp) result(ok)
 !DEC$ ATTRIBUTES DLLEXPORT :: check_mapmode_
 !! author: MDG
 !! version: 1.0
 !! date: 02/02/20
 !!
 !! check whether or not the requested mapmode exists
 
 IMPLICIT NONE
 
 class(kvectors_T),INTENT(INOUT)       :: self
 character(fnlen),INTENT(IN),OPTIONAL  :: mp
 logical                               :: ok
 
 integer(kind=irg) :: i
 character(20)     :: modes(6) = (/ 'Conical             ', &
                                    'ECCI                ', &
                                    'Standard            ', &
                                    'StandardConical     ', &
                                    'RoscaLambert        ', &
                                    'RoscaLambertLegendre' /)
 
 ok = .FALSE.
 if (present(mp)) then
   do i = 1, 5
     if (trim(modes(i)).eq.trim(mp)) ok = .TRUE.
   end do
 else
   do i = 1, 5
     if (trim(modes(i)).eq.trim(self%mapmode)) ok = .TRUE.
   end do
 end if
 
 end function check_mapmode_
 
 !--------------------------------------------------------------------------
 recursive subroutine set_mapmode_(self, mp)
 !DEC$ ATTRIBUTES DLLEXPORT :: set_mapmode_
 !! author: MDG
 !! version: 1.0
 !! date: 02/12/20
 !!
 !! set (and check) the map mode
 
 use mod_io
 
 IMPLICIT NONE
 
 class(kvectors_T), INTENT(INOUT)  :: self
 character(*), INTENT(IN)          :: mp
 
 type(IO_T)                        :: Message
 
 self%mapmode = trim(mp)
 
 if (.not.self%check_mapmode()) then
   call Message%printError('set_mapmode','kvector mapping mode '//trim(mp)//' not known')
 end if
 
 end subroutine set_mapmode_
 
 !--------------------------------------------------------------------------
 recursive function get_mapmode_(self) result(mp)
 !DEC$ ATTRIBUTES DLLEXPORT :: get_mapmode_
 !! author: MDG
 !! version: 1.0
 !! date: 02/13/20
 !!
 !! get the map mode
 
 IMPLICIT NONE
 
 class(kvectors_T), INTENT(INOUT)  :: self
 character(fnlen)                  :: mp
 
 mp = trim(self%mapmode)
 
 end function get_mapmode_
 
 !--------------------------------------------------------------------------
 recursive subroutine set_kinp_(self, k)
 !DEC$ ATTRIBUTES DLLEXPORT :: set_kinp_
 !! author: MDG
 !! version: 1.0
 !! date: 02/12/20
 !!
 !! set the input wave vector
 
 IMPLICIT NONE
 
 class(kvectors_T), INTENT(INOUT)  :: self
 real(kind=dbl), INTENT(IN)        :: k(3)
 
 self%kinp = k
 
 end subroutine set_kinp_
 
 !--------------------------------------------------------------------------
 recursive subroutine set_ktmax_(self, k)
 !DEC$ ATTRIBUTES DLLEXPORT :: set_ktmax_
 !! author: MDG
 !! version: 1.0
 !! date: 02/12/20
 !!
 !! set the max tangential component
 
 IMPLICIT NONE
 
 class(kvectors_T), INTENT(INOUT)  :: self
 real(kind=dbl), INTENT(IN)        :: k
 
 self%ktmax = k
 
 end subroutine set_ktmax_
 
 !--------------------------------------------------------------------------
 recursive subroutine set_SamplingType_(self, i)
 !DEC$ ATTRIBUTES DLLEXPORT :: set_SamplingType_
 !! author: MDG
 !! version: 1.0
 !! date: 02/12/20
 !!
 !! set the the sampling type parameter isym
 
 IMPLICIT NONE
 
 class(kvectors_T), INTENT(INOUT)  :: self
 integer(kind=irg), INTENT(IN)     :: i
 
 self%isym = i
 
 end subroutine set_SamplingType_
 
 !--------------------------------------------------------------------------
 recursive function get_numk_(self) result(numk)
 !DEC$ ATTRIBUTES DLLEXPORT :: get_numk_
 !! author: MDG
 !! version: 1.0
 !! date: 02/12/20
 !!
 !! return the number of k-vectors
 
 IMPLICIT NONE
 
 class(kvectors_T), INTENT(INOUT)  :: self
 integer(kind=irg)                 :: numk
 
 numk = self%numk
 
 end function get_numk_
 
 !--------------------------------------------------------------------------
 recursive function get_ListHead_(self) result(klist)
 !DEC$ ATTRIBUTES DLLEXPORT :: get_ListHead_
 !! author: MDG
 !! version: 1.0
 !! date: 02/12/20
 !!
 !! return the number of k-vectors
 
 IMPLICIT NONE
 
 class(kvectors_T), INTENT(INOUT)  :: self
 type(kvectorlist), pointer        :: klist
 
 klist => self%klist
 
 end function get_ListHead_
 
 !--------------------------------------------------------------------------
 recursive subroutine Calckvectors_(self,cell,SG,Diff,ga,npx,npy,ijmax,usehex,LegendreArray)
 !DEC$ ATTRIBUTES DLLEXPORT :: Calckvectors_
  !! author: MDG
  !! version: 1.0
  !! date: 02/02/20
  !!
  !! create a linked list of wave vectors
  !!
  !! This is a new version that combines several older routines.  The most important
  !! aspects of this routine are a) linked list can use regular mapping or modified Lambert mapping;
  !! b) list makes use of crystal symmetry (although that feature can be turned off); c) routine
  !! has been cleaned up, and there is now a Delete_kvectorlist function as well.  This is a very
  !! complex routine so make sure you fully understand it before you attempt to modify or add anything!
  !!
  !! todo: The Standard and RoscaLambert mapmodes have different considerations of the
  !! Laue groups; this needs to be verified and, if necessary, simplified to a single set of
  !! conditions.  This might also allow Addkvector and Add_knode to become a single routine.
 
 use mod_io
 use mod_diffraction
 use mod_crystallography
 use mod_symmetry
 
 IMPLICIT NONE
 
 class(kvectors_T),INTENT(INOUT)         :: self
 type(Cell_T),INTENT(INOUT)              :: cell
 type(SpaceGroup_T),INTENT(INOUT)        :: SG
 type(Diffraction_T),INTENT(INOUT)       :: Diff
 real(kind=dbl),INTENT(IN)               :: ga(3)
  !! "horizontal" reciprocal lattice vector
 integer(kind=irg),INTENT(IN)            :: npx
  !! number of kvectors along x
 integer(kind=irg),INTENT(IN)            :: npy
  !! number of kvectors along y
 integer(kind=irg),INTENT(INOUT)         :: ijmax
  !! max parameter used for Conical and StandardConical modes
 real(kind=dbl),INTENT(IN),OPTIONAL      :: LegendreArray(0:2*npx)
  !! Legendre lattitude grid points for spherical indexing
 logical,INTENT(IN),OPTIONAL             :: usehex
  !! hexagonal mode for RoscaLambert mapmode
 
 type(IO_T)                              :: Message
 integer(kind=irg)                       :: istat,i,j,istart,iend,jstart,jend, imin, imax, jmin, jmax, ii, jj, sqring
 real(kind=dbl)                          :: glen, xy(2), xx, yy, eps
 logical                                 :: hexgrid = .FALSE., yes = .TRUE., flip = .TRUE., check
 character(3)                            :: grid
 type(kvectorlist),pointer               :: ktail, ktmp
 real(kind=sgl)                          :: xytest(2), xxtest, yytest
 
 ! first, if self%klist already exists, delete it
  if (associated(self%klist)) then  ! deallocate the entire linked list
    call self%Delete_kvectorlist()
  end if
 
 ! do we know this mapmode ?
 if ( self%check_mapmode().eqv..FALSE.) then
   call Message%printError('Calckvectors','mapmode unknown')
 end if
 
 if (trim(self%mapmode).eq.'ECCI') then ! used for ECCI without symmetry application, including EMZAdefect
   ! compute geometrical factors
    glen = cell%CalcLength(ga,'r')                         ! length of ga
    self%gan = ga/glen                                     ! normalized ga
    !self%delta = 2.0*self%ktmax*glen/(2.0*float(npx)+1.0)
    self%delta = self%ktmax*glen/dble(npx)
    ! print*,self%delta, self%gan   ! grid step size in nm-1
    call cell%TransSpace(self%kinp,self%kstar,'d','r')     ! transform incident direction to reciprocal space
    call cell%CalcCross(ga,self%kstar,self%gperp,'r','r',0)! compute g_perp = ga x k
    call cell%NormVec(self%gperp,'r')                      ! normalize g_perp
    call cell%NormVec(self%kstar,'r')                      ! normalize reciprocal beam vector
   
   ! allocate the head and tail of the linked list
    allocate(self%klist,stat=istat)                        ! allocate new value
    if (istat.ne.0) call Message%printError('Calckvectors','unable to allocate self%klist pointer')
    ktail => self%klist                                    ! tail points to new value
    nullify(ktail%next)                                    ! nullify next in new value
    self%numk = 1                                          ! keep track of number of k-vectors so far
    ktail%i = 0                                            ! i-index of beam
    ktail%j = 0                                            ! j-index of beam
    ktail%hs = 0                                           ! spare index
    ktail%kt = (/0.0,0.0,0.0/)                             ! no tangential component for central beam direction
    ktail%k = self%kstar/Diff%getWaveLength()              ! divide by wavelength
    ktail%kn = cell%CalcDot(ktail%k,self%kstar,'r')        ! normal component
   
   ! set the loop limits
    imin = -npx; imax = npx; jmin = -npy; jmax = npy;
   
   ! and loop over the entire range (without symmetry considerations
    do i=imin,imax
     do j=jmin,jmax
      if (.not.((i.eq.0).and.(j.eq.0))) then               ! the point (0,0) has already been taken care of
       if ((i**2+j**2).le.ijmax) then                      ! only directions inside the incident cone
        call self%Add_knode(cell,Diff,ktail,i,j,(/ 0.0,0.0/))   ! add k-vector to linked list
       end if
      end if
     end do
    end do
 end if  ! mapmode = ECCI
   
 
 if (trim(self%mapmode).eq.'Conical') then ! used for CBED without symmetry application, including EMZAdefect
 ! compute geometrical factors
  glen = cell%CalcLength(ga,'r')                         ! length of ga
  self%gan = ga/glen                                     ! normalized ga
  self%delta = 2.0*self%ktmax*glen/(2.0*float(npx)+1.0)  ! grid step size in nm-1
  call cell%TransSpace(self%kinp,self%kstar,'d','r')     ! transform incident direction to reciprocal space
  call cell%CalcCross(ga,self%kstar,self%gperp,'r','r',0)! compute g_perp = ga x k
  call cell%NormVec(self%gperp,'r')                      ! normalize g_perp
  call cell%NormVec(self%kstar,'r')                      ! normalize reciprocal beam vector
 
 ! allocate the head and tail of the linked list
  allocate(self%klist,stat=istat)                        ! allocate new value
  if (istat.ne.0) call Message%printError('Calckvectors','unable to allocate self%klist pointer')
  ktail => self%klist                                    ! tail points to new value
  nullify(ktail%next)                                    ! nullify next in new value
  self%numk = 1                                          ! keep track of number of k-vectors so far
  ktail%i = 0                                            ! i-index of beam
  ktail%j = 0                                            ! j-index of beam
  ktail%kt = (/0.0,0.0,0.0/)                             ! no tangential component for central beam direction
  ktail%k = self%kstar/Diff%getWaveLength()              ! divide by wavelength
  ktail%kn = cell%CalcDot(ktail%k,self%kstar,'r')        ! normal component
 
 ! set the loop limits
  imin = -npx; imax = npx; jmin = -npy; jmax = npy;
 
 ! and loop over the entire range (without symmetry considerations
  do i=imin,imax
   do j=jmin,jmax
    if (.not.((i.eq.0).and.(j.eq.0))) then               ! the point (0,0) has already been taken care of
     if ((i**2+j**2).le.ijmax) then                      ! only directions inside the incident cone
      call self%Add_knode(cell,Diff,ktail,i,j,(/ 0.0,0.0/))   ! add k-vector to linked list
     end if
    end if
   end do
  end do
 end if  ! mapmode = Conical
 
 
 ! standard or standard-conical kvector list, as used by CBED and other programs
 if ( (self%mapmode.eq.'Standard').or.(self%mapmode.eq.'StandardConical') ) then
 
 ! for standard mode, we want to make sure that ijmax, which need not be defined by
 ! the calling program for this mode, is set to a large value
   if (self%mapmode.eq.'Standard') then
     ijmax = (5*npx)**2
   end if
 
 ! compute geometrical factors
  glen = cell%CalcLength(ga,'r')                         ! length of ga
  self%gan = ga/glen                                     ! normalized ga
  self%delta = 2.0*self%ktmax*glen/(2.0*float(npx)+1.0)  ! grid step size in nm-1
  call cell%TransSpace(self%kinp,self%kstar,'d','r')     ! transform incident direction to reciprocal space
  call cell%CalcCross(ga,self%kstar,self%gperp,'r','r',0)! compute g_perp = ga x k
  call cell%NormVec(self%gperp,'r')                      ! normalize g_perp
  call cell%NormVec(self%kstar,'r')                      ! normalize reciprocal beam vector
 
 ! allocate the head and tail of the linked list
  allocate(self%klist,stat=istat)                        ! allocate new value
  if (istat.ne.0) call Message%printError('Calckvectors','unable to allocate self%klist pointer')
  ktail => self%klist                                    ! tail points to new value
  nullify(ktail%next)                                    ! nullify next in new value
  self%numk = 1                                               ! keep track of number of k-vectors so far
  ktail%i = 0                                            ! i-index of beam
  ktail%j = 0                                            ! j-index of beam
  ktail%kt = (/0.0,0.0,0.0/)                             ! no tangential component for central beam direction
  ktail%k = self%kstar/Diff%getWaveLength()              ! divide by wavelength
  ktail%kn = cell%CalcDot(ktail%k,self%kstar,'r')        ! normal component
 
 ! implement symmetry Table 7.3 from EM book
   select case(self%isym)  ! negative values -> systematic row; positive -> zone axis
    case(-1)  ! centrosymmetric systematic row
      imin = 0; imax = npx; grid = 'srw'
    case(-2)  ! non-centrosymmetric systematic row
      imin = -npx; imax = npx; grid = 'srw'
    case(1)  ! 2D Group 1
      imin = -npx; imax = npx; jmin = -npy; jmax = npy; grid = 'sqa'
    case(2)  ! 2D Group 2
      imin = -npx; imax = npx; jmin = 0; jmax = npy; grid = 'sqb'
    case(3)  ! 2D Group m
      imin = -npx; imax = npx; jmin = 0; jmax = npy; grid = 'sqa'
    case(4)  ! 2D Group 2mm
      imin = 0; imax = npx; jmin = 0; jmax = npy; grid = 'sqa'
    case(5)  ! 2D Group 4
      imin = 1; imax = npx; jmin = 0; jmax = npy; grid = 'sqa'
    case(6)  ! 2D Group 4mm
      imin = 0; imax = npx; jmin = 0; jmax = npy; grid = 'sqc'
    case(7)  ! 2D Group 3   (cubic version)
      grid = 'hxa'; hexgrid=.TRUE.
    case(8)  ! 2D Group 31m  (cubic version)
      grid = 'hxb'; hexgrid=.TRUE.
    case(9)  ! 2D Group 6
      grid = 'hxe'; hexgrid=.TRUE.
    case(10)  ! 2D Group 6mm
      grid = 'hxf'; hexgrid=.TRUE.
    case(11)  ! 2D Group 3  (hexagonal setting)
      grid = 'hxc'; hexgrid=.TRUE.
    case(12)  ! 2D Group 31m (hexagonal setting)
      grid = 'hxd'; hexgrid=.TRUE.
    case(13)  ! 2D Group 3m1 (cubic setting)
      grid = 'hxg'; hexgrid=.TRUE.
    case(14)  ! 2D Group 3m1 (hexagonal setting)
      grid = 'hxh'; hexgrid=.TRUE.
    case default   ! we should never get here
      call Message%printError('Calckvectors','unknown isym value')
   end select
 
 ! now do the real work for standard sets of wave vectors
   select case(grid)
 
    case('srw')          ! systematic row incident beam orientations
      do i=imin,imax
       if (i.ne.0) then                                  ! the point (0,0) has already been taken care of
        call self%Add_knode(cell,Diff,ktail,i,0,(/ 0.0,0.0/))
       end if
      end do
 
    case('sqa')          ! from here on, all orientations are zone axis cases for all Laue groups
      do i=imin,imax
       do j=jmin,jmax
        if (.not.((i.eq.0).and.(j.eq.0))) then   ! the point (0,0) has already been taken care of
         if (i**2+j**2.le.ijmax) then
          call self%Add_knode(cell,Diff,ktail,i,j,(/ 0.0,0.0/))
         end if
        end if
       end do
      end do
 
    case('sqb')
      do i=imin,imax
       jloop_sqb:  do j=jmin,jmax
        if ((j.eq.0).and.(i.lt.0)) cycle jloop_sqb       ! skip the points  (i<0,0)
        if (.not.((i.eq.0).and.(j.eq.0))) then   ! the point (0,0) has already been taken care of
         if (i**2+j**2.le.ijmax) then
          call self%Add_knode(cell,Diff,ktail,i,j,(/ 0.0,0.0/))
         end if
        end if
       end do jloop_sqb
      end do
 
    case('sqc')
      do j=0,jmax
       do i=j,imax
        if (.not.((i.eq.0).and.(j.eq.0))) then   ! the point (0,0) has already been taken care of
         if (i**2+j**2 .le. ijmax) then
          call self%Add_knode(cell,Diff,ktail,i,j,(/ 0.0,0.0/))
         end if
        end if
       end do
      end do
 
    case('hxa')
      do j=0,npy
       do i=1-Kdelta(j,0),npx
        if (.not.((i.eq.0).and.(j.eq.0))) then   ! the point (0,0) has already been taken care of
         if (i**2+j**2.le.ijmax) then
          call self%Add_knode(cell,Diff,ktail,i,j,(/ 0.0,0.0/),hexgrid)
         end if
        end if
       end do
     end do
 
    case('hxb')
      do j=0,npy
       do i=j,npx
        if (.not.((i.eq.0).and.(j.eq.0))) then   ! the point (0,0) has already been taken care of
         if (i**2+j**2.le.ijmax) then
          call self%Add_knode(cell,Diff,ktail,i,j,(/ 0.0,0.0/),hexgrid)
         end if
        end if
       end do
     end do
 
    case('hxc')
      do j=0,npy
       do i=1-Kdelta(j,0)-j,npx-j
        if (.not.((i.eq.0).and.(j.eq.0))) then   ! the point (0,0) has already been taken care of
         if (i**2+j**2.le.ijmax) then
          call self%Add_knode(cell,Diff,ktail,i,j,(/ 0.0,0.0/),hexgrid)
         end if
        end if
       end do
      end do
 
    case('hxd')
      do j=0,npy
       do i=0,npx-j
        if (.not.((i.eq.0).and.(j.eq.0))) then   ! the point (0,0) has already been taken care of
         if (i**2+j**2.le.ijmax) then
          call self%Add_knode(cell,Diff,ktail,i,j,(/ 0.0,0.0/),hexgrid)
         end if
        end if
       end do
      end do
 
    case('hxe')
      do j=0,npy-1
       do i=1-Kdelta(j,0),npx-j
        if (.not.((i.eq.0).and.(j.eq.0))) then   ! the point (0,0) has already been taken care of
         if (i**2+j**2.le.ijmax) then
          call self%Add_knode(cell,Diff,ktail,i,j,(/ 0.0,0.0/),hexgrid)
        end if
       end if
       end do
      end do
 
    case('hxf')
      do j=0,npy/2
       do i=j,npx-j
        if (.not.((i.eq.0).and.(j.eq.0))) then   ! the point (0,0) has already been taken care of
         if (i**2+j**2.le.ijmax) then
          call self%Add_knode(cell,Diff,ktail,i,j,(/ 0.0,0.0/),hexgrid)
         end if
        end if
       end do
      end do
 
    case('hxg')
      do j=0,npy
       do i=j/2,min(2*j,npy)
        if (.not.((i.eq.0).and.(j.eq.0))) then   ! the point (0,0) has already been taken care of
         if (i**2+j**2.le.ijmax) then
          call self%Add_knode(cell,Diff,ktail,i,j,(/ 0.0,0.0/),hexgrid)
         end if
        end if
       end do
      end do
 
    case('hxh')
      do j=0,npy
       do i=-j/2,min(j,npy-1)
        if (.not.((i.eq.0).and.(j.eq.0))) then   ! the point (0,0) has already been taken care of
         if (i**2+j**2.le.ijmax) then
          call self%Add_knode(cell,Diff,ktail,i,j,(/ 0.0,0.0/),hexgrid)
         end if
        end if
       end do
      end do
 
   case default  ! we should never get here
     call Message%printError('Calckvectors:','unknown grid type value')
 
   end select  ! grid value
 
 end if ! mapmode.eq.'Standard' or 'StandardConical'
 
 
 ! the next type of grid is the one used for the modified Lambert maps in the dynamical EBSD
 ! programs; this requires some special care, since these mappings are a little trickier than
 ! those of the standard mapmode.  While it is possible to use a plain Lambert projection as
 ! well, here we only allow for the RoscaLambert mode.
 
 if (self%mapmode.eq.'RoscaLambert') then
    self%delta =  1.D0 / dble(npx)
    hexgrid = .FALSE.
    if (present(usehex)) then             ! hexagonal grid if .TRUE.
      if (usehex.eqv..TRUE.) hexgrid = .TRUE.
    end if
 
 ! allocate the head of the linked list
    allocate(self%klist,stat=istat)              ! allocate new value
    if (istat.ne.0) call Message%printError('Calckvectors',' unable to allocate self%klist pointer')
    ktail => self%klist                          ! tail points to new value
    nullify(ktail%next)                          ! nullify next in new value
    self%numk = 1                                     ! keep track of number of k-vectors so far
    ktail%hs = 1                                 ! this lies in the Northern Hemisphere
    ktail%i = 0                                  ! i-index of beam
    ktail%j = 0                                  ! j-index of beam
    self%kstar = (/ 0.D0, 0.D0, 1.D0 /)          ! we always use c* as the center of the RoscaLambert projection
    call cell%NormVec(self%kstar,'c')            ! normalize incident direction
    self%kstar = self%kstar/Diff%getWaveLength() ! divide by wavelength
 ! and transform to reciprocal crystal space using the direct structure matrix
    ktail%k = matmul(transpose(cell%getdsm()),self%kstar)
    ktail%kn = 1.0/Diff%getWaveLength()
 
 ! MDG: as of 8/25/15, we no longer use the Laue groups to determine the set of independent wave vectors,
 ! but instead we use the complete point group symmetry, as it should be.  Upon reflection, using
 ! the Laue groups was equivalent to implicitly using Friedel's law, which makes all diffraction patterns
 ! centrosymmetric, and that is not correct for EBSD. So, the symbol isym now encodes the full point group,
 ! not the Laue group.  This will require a modification in each calling program as well.
 
 ! in addition, the modified Lambert projection will now require two hemispheres (NH and SH). We can handle this
 ! by means of an optional argument to the AddkVector routine; when the argument is present, the -k_z version
 ! of the direction is also added to the list.
 
 ! deal with each point group symmetry separately or in sets, depending on the value of isym
  select case (self%isym)
 
    case (1)  ! triclinic 1
         istart = -npx
         iend = npx
         jstart = -npy
         jend = npy
           do j=jstart,jend
             do i=istart,iend   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 call self%AddkVector(cell,Diff,ktail,xy,i,j,addSH = yes)
             end do
           end do
 
    case (2)  ! triclinic -1
         istart = -npx
         iend = npx
         jstart = -npy
         jend = npy
           do j=jstart,jend
             do i=istart,iend   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 call self%AddkVector(cell,Diff,ktail,xy,i,j)
             end do
           end do
 
   case (3)   !  monoclinic 2
         istart = 0
         iend = npx
         jstart = -npy
         jend = npy
           do j=jstart,jend
            do i=istart,iend   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 call self%AddkVector(cell,Diff,ktail,xy,i,j, addSH = yes)
            end do
           end do
 
   case (4)   !  monoclinic m
         istart = -npx
         iend = npx
         jstart = 0
         jend = npy
           do j=jstart,jend
            do i=istart,iend   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 call self%AddkVector(cell,Diff,ktail,xy,i,j, addSH = yes)
            end do
           end do
 
   case (5)  ! monoclinic 2/m, orthorhombic 222, mm2, tetragonal 4, -4
         istart = 0
         iend = npx
         jstart = 0
         jend = npy
           do j=jstart,jend
            do i=istart,iend   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 call self%AddkVector(cell,Diff,ktail,xy,i,j, addSH = yes)
            end do
           end do
 
   case (6)  ! orthorhombic mmm, tetragonal 4/m, 422, -4m2, cubic m-3, 432 (for now)
         istart = 0
         iend = npx
         jstart = 0
         jend = npy
           do j=jstart,jend
            do i=istart,iend   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 call self%AddkVector(cell,Diff,ktail,xy,i,j)
            end do
           end do
 
   case (7)  ! tetragonal 4mm
         istart = 0
         iend = npx
         jstart = 0
         jend = npx
           do i=istart,iend
            do j=jstart,i   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 call self%AddkVector(cell,Diff,ktail,xy,i,j, addSH = yes)
            end do
           end do
 
   case (8)  ! tetragonal -42m, cubic -43m (for now)
         istart = 0
         iend = npx
         jstart = -npx
         jend = npx
           do i=istart,iend
            do j=-i, i   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 call self%AddkVector(cell,Diff,ktail,xy,i,j)
            end do
           end do
 
   case (9)  ! tetragonal 4/mmm, cubic m-3m (for now)
         istart = 0
         iend = npx
         jstart = 0
         jend = npx
           do i=istart,iend
            do j=jstart,i   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 call self%AddkVector(cell,Diff,ktail,xy,i,j)
            end do
           end do
 
 ! cases 10 through 19 are all on a hexagonal grid...
 ! for now (08/31/15), we have not yet implemented the rhombohedral setting of the trigonal space groups;
 ! this case is truly a pain in the neck to implement...
 
   case (10)   ! hexagonal 3
         istart = 0
         iend = npx
         jstart = 0
         jend = npx
           do j=jstart,jend
             do i=istart,iend   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 if (InsideHexGrid(xy)) call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid, addSH = yes)
             end do
           end do
 
   case (11)   ! rhombohedral 3
         call Message%printError('Calckvectors: ','rhombohedral setting not yet implemented, use hexagonal setting instead')
 !       istart = 0
 !       iend = npx
 !       jstart = 0
 !       jend = npx
 !         do j=jstart,jend
 !           do i=istart,iend   !
 !               ii = 2*j-i
 !               jj = j-2*i
 !               xy = (/ dble(ii), dble(jj) /) * delta * LPs%isrt
 !               if (InsideHexGrid(xy)) call self%AddkVector(cell,Diff,ktail,xy,ii,jj,hexgrid, addSH = yes)
 !           end do
 !         end do
 
   case (12)   ! hexagonal -3, 321, -6; [not implemented: rhombohedral 32]
         if ((SG%getSpaceGrouptrigonal()).and.(SG%getSpaceGroupsecond())) then
           call Message%printError('Calckvectors: ','rhombohedral setting not yet implemented, use hexagonal setting instead')
         else
           istart = 0
           iend = npx
           jstart = 0
           jend = npx
             do j=jstart,jend
               do i=istart,iend   !
                   xy = (/ dble(i), dble(j) /) * self%delta
                   if (InsideHexGrid(xy)) call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid)
               end do
             end do
         end if
 
   case (13)   ! [not implemented: rhombohedral -3], hexagonal 312  [ modified 7/31/18, MDG ]
         if ((SG%getSpaceGrouptrigonal()).and.(SG%getSpaceGroupsecond())) then
           call Message%printError('Calckvectors: ','rhombohedral setting not yet implemented, use hexagonal setting instead')
         else
           istart = 0
           iend = npx
           jstart = 0
           jend = -npx
             do j=jstart,jend,-1
               do i=istart+j/2,iend   !
                 xy = (/ dble(i),  dble(j) /) * self%delta
                 if (InsideHexGrid(xy)) then
                   ! call self%AddkVector(cell,Diff,ktail,xy,-i,-j,hexgrid)
                   ! ktail%k(2) = -ktail%k(2)
                  call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid)
                 end if
               end do
             end do
           istart = 0
           iend = npx
           jstart = 0
           jend = npx
             do i=istart,iend   !
               do j=jstart,i/2
                 xy = (/ dble(i),  dble(j) /) * self%delta
                 if (InsideHexGrid(xy)) then
                   call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid)
                  ! call self%AddkVector(cell,Diff,ktail,xy,-i,-j,hexgrid)
                  !  ktail%k(2) = -ktail%k(2)
                 end if
               end do
             end do
         end if
 
   case (14)   ! hexagonal 3m1, [not implemented: rhombohedral 3m]
         if ((SG%getSpaceGrouptrigonal()).and.(SG%getSpaceGroupsecond())) then
           call Message%printError('Calckvectors: ','rhombohedral setting not yet implemented, use hexagonal setting instead')
         else
           istart = 0
           iend = npx
           jstart = 0
           jend = npx
             do j=jstart,jend
               do i=istart+(j-1)/2,2*j
                 xy = (/ dble(i),  dble(j) /) * self%delta
                 if (InsideHexGrid(xy)) then
                   call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid, addSH = yes)
                 end if
               end do
             end do
         end if
 
   case (15)   ! hexagonal 31m, 6
           istart = 0
           iend = npx
           jstart = 0
           jend = npx
             do j=jstart,jend
               do i=istart+j,jend
                 xy = (/ dble(i),  dble(j) /) * self%delta
                 if (InsideHexGrid(xy)) then
                   call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid, addSH = yes)
                 end if
               end do
             end do
 
   case (16)   ! hexagonal -3m1, 622, -6m2 [not implemented: rhombohedral -3m]
         if ((SG%getSpaceGrouptrigonal()).and.(SG%getSpaceGroupsecond())) then
           call Message%printError('Calckvectors: ','rhombohedral setting not yet implemented, use hexagonal setting instead')
         else
           istart = 0
           iend = npx
           jstart = 0
           jend = npx
           eps = 1.0D-4
             do j=jstart,jend
               do i=istart,iend
                   xy = (/ dble(i), dble(j) /) * self%delta
                   xx = dble(i)-dble(j)/2.D0
                   yy = dble(j)*LPs%srt
                   check = .TRUE.
                   if (xx.lt.0.D0) then
                     check = .FALSE.
                   else
                     if (xx.ge.0.D0) then
                       yy = datan2(yy,xx)
                       if (yy .lt. (LPs%Pi/6.D0-eps)) check = .FALSE.
                     end if
                   end if
                   if (InsideHexGrid(xy).and.(check)) call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid)
               end do
             end do
         end if
 
   case (17)   ! hexagonal -31m, 6/m, -62m
         istart = 0
         iend = npx
         jstart = 0
         jend = npx
         eps = 1.0D-4
           do j=jstart,jend
             do i=istart,iend
                 xy = (/ dble(i), dble(j) /) * self%delta
                 xx = dble(i)-dble(j)/2.D0
                 yy = dble(j)*LPs%srt
                 check = .TRUE.
                 if (xx.lt.0.D0) then
                    check = .FALSE.
                 else
                    if (xx.ge.0.D0) then
                      yy = datan2(yy,xx)
                      if (yy.gt.(cPi/3.D0+eps)) check = .FALSE.
                    end if
                 end if
                 if (InsideHexGrid(xy).and.(check)) call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid)
             end do
           end do
 
   case (18)   ! hexagonal 6mm
         istart = 0
         iend = npx
         jstart = 0
         jend = npx
         eps = 1.0D-4
           do j=jstart,jend
             do i=istart,iend
                 xy = (/ dble(i), dble(j) /) * self%delta
                 xx = dble(i)-dble(j)/2.D0
                 yy = dble(j)*LPs%srt
                 check = .TRUE.
                 if (xx.lt.0.D0) then
                    check = .FALSE.
                 else
                    if (xx.ge.0.D0) then
                      yy = datan2(yy,xx)
                      if (yy.gt.(cPi/6.D0+eps)) check = .FALSE.
                    end if
                 end if
                 if (InsideHexGrid(xy).and.(check)) call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid, addSH = yes)
             end do
           end do
 
   case (19)   ! hexagonal 6/mmm
         istart = 0
         iend = npx
         jstart = 0
         jend = npx
         eps = 1.0D-4
           do j=jstart,jend
             do i=istart,iend
                 xy = (/ dble(i), dble(j) /) * self%delta
                 xx = dble(i)-dble(j)/2.D0
                 yy = dble(j)*LPs%srt
 
                 check = .TRUE.
                 if (xx.lt.0.D0) then
                    check = .FALSE.
                 else
                    if (xx.ge.0.D0) then
                      yy = datan2(yy, xx)
                      if (yy.gt.(LPs%Pi/6.D0+eps)) check = .FALSE.
                    end if
 
                 end if
                 if (InsideHexGrid(xy).and.(check)) call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid)
 
             end do
           end do
 
  end select
 
 end if
 
 
 ! the next type of grid is the one used for the modified Lambert maps with Legendre lattitude grid values in the
 ! dynamical EBSD programs; this requires some special care, since these mappings are a little
 ! trickier than those of the standard mapmode.  While it is possible to use a plain Lambert
 ! projection as well, here we only allow for the RoscaLambert mode with modified lattitudinal angles.
 
 if (self%mapmode.eq.'RoscaLambertLegendre') then
    self%delta =  1.D0 / dble(npx)
    if (usehex) then             ! hexagonal grid
       hexgrid = .TRUE.
    else                         ! square grid
       hexgrid = .FALSE.
    end if
 
 ! allocate the head of the linked list
    allocate(self%klist,stat=istat)                   ! allocate new value
    if (istat.ne.0) call Message%printError('Calckvectors',' unable to allocate self%klist pointer')
    ktail => self%klist                               ! tail points to new value
    nullify(ktail%next)                          ! nullify next in new value
    self%numk = 1                                     ! keep track of number of k-vectors so far
    ktail%hs = 1                                 ! this lies in the Northern Hemisphere
    ktail%i = 0                                  ! i-index of beam
    ktail%j = 0                                  ! j-index of beam
    self%kstar = (/ 0.0, 0.0, 1.0 /)             ! we always use c* as the center of the RoscaLambert projection
    call cell%NormVec(self%kstar,'c')            ! normalize incident direction
    self%kstar = self%kstar/Diff%getWaveLength() ! divide by wavelength
 ! and transform to reciprocal crystal space using the structure matrix
    ktail%k = matmul(transpose(cell%getdsm()),self%kstar)
    ktail%kn = 1.0/Diff%getWaveLength()
 
 ! MDG: as of 8/25/15, we no longer use the Laue groups to determine the set of independent wave vectors,
 ! but instead we use the complete point group symmetry, as it should be.  Upon reflection, using
 ! the Laue groups was equivalent to implicitly using Friedel's law, which makes all diffraction patterns
 ! centrosymmetric, and that is not correct for EBSD. So, the symbol isym now encodes the full point group,
 ! not the Laue group.  This will require a modification in each calling program as well.
 
 ! in addition, the modified Lambert projection will now require two hemispheres (NH and SH). We can handle this
 ! by means of an optional argument to the AddkVector routine; when the argument is present, the -k_z version
 ! of the direction is also added to the list.
 
 ! The main difference with the regular RoscaLambert case is the fact that the lattitudinal value of the direction
 ! cosines needs to be replaced by the one from the Legendre array, and then the in-plane dorection cosines need to
 ! be properly scaled; this type of master pattern is then used for Spherical Indexing in the EMSphInx program.
 
 ! deal with each point group symmetry separately or in sets, depending on the value of isym
  select case (self%isym)
 
    case (1)  ! triclinic 1
         istart = -npx
         iend = npx
         jstart = -npy
         jend = npy
           do j=jstart,jend
             do i=istart,iend   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 sqring = maxval( (/ abs(i), abs(j) /) )
                 call self%AddkVector(cell,Diff,ktail,xy,i,j,addSH = yes,LegendreLattitude=LegendreArray(sqring))
             end do
           end do
 
    case (2)  ! triclinic -1
         istart = -npx
         iend = npx
         jstart = -npy
         jend = npy
           do j=jstart,jend
             do i=istart,iend   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 sqring = maxval( (/ abs(i), abs(j) /) )
                 call self%AddkVector(cell,Diff,ktail,xy,i,j,LegendreLattitude=LegendreArray(sqring))
             end do
           end do
 
   case (3)   !  monoclinic 2
         istart = 0
         iend = npx
         jstart = -npy
         jend = npy
           do j=jstart,jend
            do i=istart,iend   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 sqring = maxval( (/ abs(i), abs(j) /) )
                 call self%AddkVector(cell,Diff,ktail,xy,i,j, addSH = yes,LegendreLattitude=LegendreArray(sqring))
            end do
           end do
 
   case (4)   !  monoclinic m
         istart = -npx
         iend = npx
         jstart = 0
         jend = npy
           do j=jstart,jend
            do i=istart,iend   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 sqring = maxval( (/ abs(i), abs(j) /) )
                 call self%AddkVector(cell,Diff,ktail,xy,i,j, addSH = yes,LegendreLattitude=LegendreArray(sqring))
            end do
           end do
 
   case (5)  ! monoclinic 2/m, orthorhombic 222, mm2, tetragonal 4, -4
         istart = 0
         iend = npx
         jstart = 0
         jend = npy
           do j=jstart,jend
            do i=istart,iend   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 sqring = maxval( (/ abs(i), abs(j) /) )
                 call self%AddkVector(cell,Diff,ktail,xy,i,j, addSH = yes,LegendreLattitude=LegendreArray(sqring))
            end do
           end do
 
   case (6)  ! orthorhombic mmm, tetragonal 4/m, 422, -4m2, cubic m-3, 432 (for now)
         istart = 0
         iend = npx
         jstart = 0
         jend = npy
           do j=jstart,jend
            do i=istart,iend   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 sqring = maxval( (/ abs(i), abs(j) /) )
                 call self%AddkVector(cell,Diff,ktail,xy,i,j,LegendreLattitude=LegendreArray(sqring))
            end do
           end do
 
   case (7)  ! tetragonal 4mm
         istart = 0
         iend = npx
         jstart = 0
         jend = npx
           do i=istart,iend
            do j=jstart,i   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 sqring = maxval( (/ abs(i), abs(j) /) )
                 call self%AddkVector(cell,Diff,ktail,xy,i,j, addSH = yes,LegendreLattitude=LegendreArray(sqring))
            end do
           end do
 
   case (8)  ! tetragonal -42m, cubic -43m (for now)
         istart = 0
         iend = npx
         jstart = -npx
         jend = npx
           do i=istart,iend
            do j=-i, i   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 sqring = maxval( (/ abs(i), abs(j) /) )
                 call self%AddkVector(cell,Diff,ktail,xy,i,j,LegendreLattitude=LegendreArray(sqring))
            end do
           end do
 
   case (9)  ! tetragonal 4/mmm, cubic m-3m (for now)
         istart = 0
         iend = npx
         jstart = 0
         jend = npx
           do i=istart,iend
            do j=jstart,i   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                 sqring = maxval( (/ abs(i), abs(j) /) )
                 call self%AddkVector(cell,Diff,ktail,xy,i,j,LegendreLattitude=LegendreArray(sqring))
            end do
           end do
 
 ! cases 10 through 19 are all on a hexagonal grid...
 ! for now (08/31/15), we have not yet implemented the rhombohedral setting of the trigonal space groups;
 ! this case is truly a pain in the neck to implement...
 
   case (10)   ! hexagonal 3
         istart = 0
         iend = npx
         jstart = 0
         jend = npx
           do j=jstart,jend
             do i=istart,iend   !
                 xy = (/ dble(i), dble(j) /) * self%delta
                  sqring = maxval( (/ abs(i), abs(j) /) )
                if (InsideHexGrid(xy)) call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid, addSH = yes, &
                                                        LegendreLattitude=LegendreArray(sqring))
             end do
           end do
 
   case (11)   ! rhombohedral 3
         call Message%printError('Calckvectors: ','rhombohedral setting not yet implemented, use hexagonal setting instead')
 !       istart = 0
 !       iend = npx
 !       jstart = 0
 !       jend = npx
 !         do j=jstart,jend
 !           do i=istart,iend   !
 !               ii = 2*j-i
 !               jj = j-2*i
 !               xy = (/ dble(ii), dble(jj) /) * delta * LPs%isrt
 !               if (InsideHexGrid(xy)) call self%AddkVector(cell,Diff,ktail,xy,ii,jj,hexgrid, addSH = yes)
 !           end do
 !         end do
 
   case (12)   ! hexagonal -3, 321, -6; [not implemented: rhombohedral 32]
         if ((SG%getSpaceGrouptrigonal()).and.(SG%getSpaceGroupsecond())) then
           call Message%printError('Calckvectors: ','rhombohedral setting not yet implemented, use hexagonal setting instead')
         else
           istart = 0
           iend = npx
           jstart = 0
           jend = npx
             do j=jstart,jend
               do i=istart,iend   !
                   xy = (/ dble(i), dble(j) /) * self%delta
                 sqring = maxval( (/ abs(i), abs(j) /) )
                   if (InsideHexGrid(xy)) call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid, &
                                                               LegendreLattitude=LegendreArray(sqring))
               end do
             end do
         end if
 
   case (13)   ! [not implemented: rhombohedral -3], hexagonal 312  [ modified 7/31/18, MDG ]
         if ((SG%getSpaceGrouptrigonal()).and.(SG%getSpaceGroupsecond())) then
           call Message%printError('Calckvectors: ','rhombohedral setting not yet implemented, use hexagonal setting instead')
         else
           istart = 0
           iend = npx
           jstart = 0
           jend = -npx
             do j=jstart,jend,-1
               do i=istart+j/2,iend   !
                 xy = (/ dble(i),  dble(j) /) * self%delta
                 if (InsideHexGrid(xy)) then
                   ! call self%AddkVector(cell,Diff,ktail,xy,-i,-j,hexgrid)
                   ! ktail%k(2) = -ktail%k(2)
                 sqring = maxval( (/ abs(i), abs(j) /) )
                  call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid,LegendreLattitude=LegendreArray(sqring))
                 end if
               end do
             end do
           istart = 0
           iend = npx
           jstart = 0
           jend = npx
             do i=istart,iend   !
               do j=jstart,i/2
                 xy = (/ dble(i),  dble(j) /) * self%delta
                 if (InsideHexGrid(xy)) then
                 sqring = maxval( (/ abs(i), abs(j) /) )
                   call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid,LegendreLattitude=LegendreArray(sqring))
                  ! call self%AddkVector(cell,Diff,ktail,xy,-i,-j,hexgrid)
                  !  ktail%k(2) = -ktail%k(2)
                 end if
               end do
             end do
         end if
 
   case (14)   ! hexagonal 3m1, [not implemented: rhombohedral 3m]
         if ((SG%getSpaceGrouptrigonal()).and.(SG%getSpaceGroupsecond())) then
           call Message%printError('Calckvectors: ','rhombohedral setting not yet implemented, use hexagonal setting instead')
         else
           istart = 1
           iend = npx
           jstart = 1
           jend = npx
             do j=jstart,jend
               do i=istart+(j-1)/2,2*j
                 xy = (/ dble(i),  dble(j) /) * self%delta
                 sqring = maxval( (/ abs(i), abs(j) /) )
                 if (InsideHexGrid(xy)) then
                   call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid, addSH = yes,LegendreLattitude=LegendreArray(sqring))
                 end if
               end do
             end do
         end if
 
   case (15)   ! hexagonal 31m, 6
           istart = 0
           iend = npx
           jstart = 1
           jend = npx
             do j=jstart,jend
               do i=istart+j,jend
                 xy = (/ dble(i),  dble(j) /) * self%delta
                 sqring = maxval( (/ abs(i), abs(j) /) )
                 if (InsideHexGrid(xy)) then
                   call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid, addSH = yes,LegendreLattitude=LegendreArray(sqring))
                 end if
               end do
             end do
 
   case (16)   ! hexagonal -3m1, 622, -6m2 [not implemented: rhombohedral -3m]
         if ((SG%getSpaceGrouptrigonal()).and.(SG%getSpaceGroupsecond())) then
           call Message%printError('Calckvectors: ','rhombohedral setting not yet implemented, use hexagonal setting instead')
         else
           istart = 0
           iend = npx
           jstart = 0
           jend = npx
           eps = 1.0D-4
             do j=jstart,jend
               do i=istart,iend
                   xy = (/ dble(i), dble(j) /) * self%delta
                   xx = dble(i)-dble(j)/2.D0
                   yy = dble(j)*LPs%srt
                   check = .TRUE.
                   if (xx.lt.0.D0) then
                     check = .FALSE.
                   else
                     if (xx.ge.0.D0) then
                       yy = datan2(yy,xx)
                       if (yy .lt. (LPs%Pi/6.D0-eps)) check = .FALSE.
                     end if
                   end if
                 sqring = maxval( (/ abs(i), abs(j) /) )
                   if (InsideHexGrid(xy).and.(check)) call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid,&
                                                                      LegendreLattitude=LegendreArray(sqring))
               end do
             end do
         end if
 
   case (17)   ! hexagonal -31m, 6/m, -62m
         istart = 0
         iend = npx
         jstart = 0
         jend = npx
         eps = 1.0D-4
           do j=jstart,jend
             do i=istart,iend
                 xy = (/ dble(i), dble(j) /) * self%delta
                 xx = dble(i)-dble(j)/2.D0
                 yy = dble(j)*LPs%srt
                 check = .TRUE.
                 if (xx.lt.0.D0) then
                    check = .FALSE.
                 else
                    if (xx.ge.0.D0) then
                      yy = datan2(yy,xx)
                      if (yy.gt.(cPi/3.D0+eps)) check = .FALSE.
                    end if
                 end if
                 sqring = maxval( (/ abs(i), abs(j) /) )
                 if (InsideHexGrid(xy).and.(check)) call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid, &
                                                                    LegendreLattitude=LegendreArray(sqring))
             end do
           end do
 
   case (18)   ! hexagonal 6mm
         istart = 0
         iend = npx
         jstart = 0
         jend = npx
         eps = 1.0D-4
           do j=jstart,jend
             do i=istart,iend
                 xy = (/ dble(i), dble(j) /) * self%delta
                 xx = dble(i)-dble(j)/2.D0
                 yy = dble(j)*LPs%srt
                 check = .TRUE.
                 if (xx.lt.0.D0) then
                    check = .FALSE.
                 else
                    if (xx.ge.0.D0) then
                      yy = datan2(yy,xx)
                      if (yy.gt.(cPi/6.D0+eps)) check = .FALSE.
                    end if
                 end if
                 sqring = maxval( (/ abs(i), abs(j) /) )
                 if (InsideHexGrid(xy).and.(check)) call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid, &
                                                                    addSH = yes,LegendreLattitude=LegendreArray(sqring))
             end do
           end do
 
   case (19)   ! hexagonal 6/mmm
         istart = 0
         iend = npx
         jstart = 0
         jend = npx
         eps = 1.0D-4
           do j=jstart,jend
             do i=istart,iend
                 xy = (/ dble(i), dble(j) /) * self%delta
                 xx = dble(i)-dble(j)/2.D0
                 yy = dble(j)*LPs%srt
 
                 check = .TRUE.
                 if (xx.lt.0.D0) then
                    check = .FALSE.
                 else
                    if (xx.ge.0.D0) then
                      yy = datan2(yy, xx)
                      if (yy.gt.(LPs%Pi/6.D0+eps)) check = .FALSE.
                    end if
 
                 end if
                 sqring = maxval( (/ abs(i), abs(j) /) )
                 if (InsideHexGrid(xy).and.(check)) call self%AddkVector(cell,Diff,ktail,xy,i,j,hexgrid,&
                                                                    LegendreLattitude=LegendreArray(sqring))
 
             end do
           end do
 
  end select
 
 end if
 
 end subroutine Calckvectors_
 
 recursive subroutine CalckvectorsPrecession_(self,cell,Diff,k,ga,precangle,prechalfwidth,precsample,precazimuthal,numk)
 !DEC$ ATTRIBUTES DLLEXPORT :: CalckvectorsPrecession_
 !> @author Marc De Graef, Carnegie Mellon University
 !
 !> @brief create a linked list of wave vectors for precession electron diffraction
 !
 !> @details This is a new version to test whether or not we can use the whole pattern
 !> symmetry to determine the relevant list of incident wave vectors; this should be a 
 !> general routine, so that we do not need to consider each symmetry case separately.
 !> This will require a floating point version of the Apply2DPGSymmetry routine in symmetry.f90.
 !
 !> @param khead wave vector list pointer
 !> @param cell unit cell pointer
 !> @param k central wave vector
 !> @param ga reciprocal lattice vector normal to k
 !> @param precangle main precession angle [mrad]
 !> @param prechalfwidth precession beam half width [mrad]
 !> @param precsample number of samples in half width 
 !> @param precazimuthal number of samples around each precession circle
 !> @param numk total number of wave vectors in list
 !
 !> @date  03/05/14 MDG 1.0 original
 !> @date  11/28/14 MDG 2.0 rewrite without global variables
 !--------------------------------------------------------------------------
 use mod_global
 use mod_io
 use mod_diffraction
 use mod_crystallography
 use mod_Lambert
 
 IMPLICIT NONE
 
 class(kvectors_T),INTENT(INOUT)         :: self
 type(Cell_T)        ,INTENT(IN)         :: cell
 type(Diffraction_T),INTENT(INOUT)       :: Diff
 type(io_T)                              :: Message
 
 real(kind=dbl),INTENT(IN)               :: k(3)         !< initial wave vector
 real(kind=dbl),INTENT(IN)               :: ga(3)        !< "horizontal" reciprocal lattice vector
 real(kind=sgl),INTENT(IN)               :: precangle    !< precession angle in [mrad]
 real(kind=sgl),INTENT(IN)               :: prechalfwidth        !< halfwidth of tilted beam [mrad]
 integer(kind=irg),INTENT(IN)            :: precsample           !< number of kvectors along beam tilt
 integer(kind=irg),INTENT(IN)            :: precazimuthal        !< number of kvectors along circumference
 integer(kind=irg),INTENT(OUT)           :: numk         !< total number of kvectors in linked list
 
 type(kvectorlist),pointer               :: ktail, khead
 integer(kind=irg)                       :: istat,i,j, iequiv(2,12), nequiv, jj, nx, ny, il, ith
 real(kind=dbl)                          :: gp, dgp, glen, gan(3), gperp(3), kstar(3), dth
 real(kind=dbl),allocatable              :: gw(:), ct(:), st(:), th(:), mLambda
 logical                                 :: hexgrid = .FALSE.
 real(kind=sgl)                          :: kt(3),kr(3)
 real(kind=sgl)                          :: ktlen
 
 
 ! compute geometrical factors 
  glen = cell%CalcLength(ga,'r')                 ! length of ga
  gan = ga/glen
  mLambda = Diff%getWaveLength()                                ! normalized ga
  gp = 2.0*sin(precangle/1000.0)*mLambda         ! precession angle converted to reciprocal length gp in units of glen
  dgp = 0.0
  if (precsample.gt.0) then
    dgp = 2.0*sin(0.001*(precangle-prechalfwidth))/mLambda/glen/float(precsample)        ! half width step size converted to reciprocal length dgp in units of glen
  end if
  allocate(gw(2*precsample+1))                           ! sampling radii
  gw = gp + dgp * (/ (i,i=-precsample,precsample) /)     ! sampling radii
 
 ! pre-compute cosines and sines
  allocate(ct(precazimuthal),st(precazimuthal), th(precazimuthal))
  dth = 2.D0*cPi / dble(precazimuthal)
  th = (/ (i-1,i=1,precazimuthal) /) * dth
  ct = cos(th)
  st = sin(th)
  
  call cell%TransSpace(k,kstar,'d','r')               ! transform incident direction to reciprocal space
  call cell%CalcCross(ga,kstar,gperp,'r','r',0)       ! compute g_perp = ga x k
  call cell%NormVec(gperp,'r')                        ! normalize g_perp
  call cell%NormVec(kstar,'r')                        ! normalize reciprocal beam vector
 
  ! first, if self%klist already exists, delete it
  if (associated(self%klist)) then  ! deallocate the entire linked list
   call self%Delete_kvectorlist()
  end if
 
  ! allocate the head of the linked list
  allocate(self%klist,stat=istat)                   ! allocate new value
  if (istat.ne.0) call Message%printError('Calckvectors',' unable to allocate self%klist pointer')
  ktail => self%klist                               ! tail points to new value
  nullify(ktail%next)                          ! nullify next in new value
  self%numk = 0                                     ! keep track of number of k-vectors so far
  
 ! next loop around each of the precession circles
  do il = 1,2*precsample+1                               ! number of concentric circles
   do ith = 1,precazimuthal                              ! number of points along each circle
 ! make a new one in the list, except for the first one
    if (self%numk.ne.0) then
      allocate(ktail%next,stat=istat)                    ! allocate new value
      if (istat.ne.0) call Message%printError('Add_knode:',' unable to allocate pointer')
      ktail => ktail%next                        ! tail points to new value
      nullify(ktail%next)                                ! nullify next in new value
    end if
 ! and populate the fields   
    kt = - gw(il)*ct(ith)*gan - gw(il)*st(ith)*gperp     ! tangential component of k
    ktail%kt = kt                                        ! store tangential component of k
    ktlen = cell%CalcLength(kt,'r')**2                   ! squared length of tangential component
 
    kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar         ! complete wave vector
    ktail%k = kr                                         ! store in pointer list
    ktail%kn = cell%CalcDot(ktail%k,kstar,'r')           ! normal component of k
    self%numk = self%numk + 1
   end do
  end do
 
 end subroutine CalckvectorsPrecession_
 
 !--------------------------------------------------------------------------
 recursive function InsideHexGrid(xy) result(res)
 !DEC$ ATTRIBUTES DLLEXPORT :: InsideHexGrid
 !! author: MDG
 !! version: 1.0
 !! date: 02/02/20
 !!
 !! determines whether or not a point lies inside the standard hexagonal Lambert grid
 !!
 !! The hexagon has unit edge length, with one vertex at (1,0). Take the
 !! absolute values of x and y and check if this point lies inside a box with edge
 !! lengths (1, sqrt(3)/2).  If not, exit; else, check on which side of the line
 !! |x|+|y|/sqrt(3)=1 the point lies.
 
 IMPLICIT NONE
 
 real(kind=dbl),INTENT(IN)       :: xy(2)
 logical                         :: res
 real(kind=dbl)                  :: ax, ay
 
 ! assume it is a good point
 res = .TRUE.
 
 ! first of all, take the absolute values and see if the transformed point lies inside the
 ! rectangular box with edge lengths (1,sqrt(3)/2)
 ax = abs(xy(1)-0.5D0*xy(2))
 ay = abs(xy(2)*LPs%srt)
 
 if ((ax.gt.1.D0).or.(ay.gt.LPs%srt)) res = .FALSE.
 ! then check for the inclined edge
 if (res) then
   if (ax+ay*LPs%isrt .gt. 1.D0) res = .FALSE.
 end if
 
 end function InsideHexGrid
 
 !--------------------------------------------------------------------------
 recursive subroutine CalckvectorsSymmetry_(self,cell,Diff,TDPG,ga,npx,npy,ijmax,klaue,debug)
 !DEC$ ATTRIBUTES DLLEXPORT :: CalckvectorsSymmetry_
 
 use mod_io
 use mod_diffraction
 use mod_crystallography
 use mod_Lambert
 use mod_symmetry2D
 
 IMPLICIT NONE
 
 class(kvectors_T),INTENT(INOUT)         :: self
 type(Cell_T),INTENT(INOUT)              :: cell
 type(Diffraction_T),INTENT(INOUT)       :: Diff
 type(symdata2D),INTENT(INOUT)           :: TDPG
 !f2py intent(in,out) ::  TDPG
 real(kind=dbl),INTENT(IN)               :: ga(3)
  !! "horizontal" reciprocal lattice vector
 integer(kind=irg),INTENT(IN)            :: npx
  !! number of kvectors along x
 integer(kind=irg),INTENT(IN)            :: npy
  !! number of kvectors along y
 integer(kind=irg),INTENT(INOUT)         :: ijmax
  !! max parameter used for Conical and StandardConical modes
 real(kind=sgl),INTENT(IN)               :: klaue(2)
  !! fractional Laue center coordinates
 logical,INTENT(IN),OPTIONAL             :: debug
 
 type(IO_T)                              :: Message
 integer(kind=irg),allocatable           :: kselected(:,:)       ! keeps track of which k-vectors have already been considered
 integer(kind=irg)                       :: istat,i,j, iequiv(2,12), nequiv, jj, nx, ny
 real(kind=dbl)                          :: glen, Lauexy(2)
 logical                                 :: hexgrid = .FALSE.
 real(kind=sgl)                          :: kt(3),kr(3)
 real(kind=sgl)                          :: ktlen
 type(kvectorlist),pointer               :: ktail
 
 nx = 2*npx
 ny = 2*npy
 allocate(kselected(-nx:nx,-ny:ny))
 
 ! initialize the kselected array to 0
 kselected = 0
 
 ! compute geometrical factors
  glen = cell%CalcLength(ga,'r')                         ! length of ga
  Lauexy = glen * klaue                                  ! scaled Laue center coordinates
  self%gan = ga/glen                                     ! normalized ga
  self%delta = 2.0*self%ktmax*glen/(2.0*float(npx)+1.0)  ! grid step size in nm-1
  call cell%TransSpace(self%kinp,self%kstar,'d','r')     ! transform incident direction to reciprocal space
  call cell%CalcCross(ga,self%kstar,self%gperp,'r','r',0)! compute g_perp = ga x k
  call cell%NormVec(self%gperp,'r')                      ! normalize g_perp
  call cell%NormVec(self%kstar,'r')                      ! normalize reciprocal beam vector
 
 ! allocate the head and tail of the linked list
  allocate(self%klist,stat=istat)                        ! allocate new value
  if (istat.ne.0) call Message%printError('Calckvectors','unable to allocate self%klist pointer')
  ktail => self%klist                                    ! tail points to new value
  nullify(ktail%next)                                    ! nullify next in new value
  self%numk = 1                                          ! keep track of number of k-vectors so far
  ktail%i = 0                                            ! i-index of beam
  ktail%j = 0                                            ! j-index of beam
 
 ! use the Laue center coordinates to define the tangential component of the incident wave vector
  kt = - Lauexy(1)*self%gan - Lauexy(2)*self%gperp       ! tangential component of k
  ktail%kt = kt                                          ! store tangential component of k
  ktlen = cell%CalcLength(kt,'r')**2                     ! squared length of tangential component
 
  kr = kt + sqrt(1.0/Diff%getWaveLength()**2 - ktlen)*self%kstar  ! complete wave vector
  ktail%k = kr                                           ! store in pointer list
  ktail%kn = cell%CalcDot(ktail%k,self%kstar,'r')        ! normal component of k
 
  kselected(0,0) = 2
 
  if (maxval(abs(klaue)).eq.0.0) then                    ! zone axis orientation, so we should use symmetry
 ! we scan over the entire range of potential beam directions, defined by npx and npy along with
 ! the conical truncation parameter ijmax; for each point we check whether or not it has been considered
 ! before; it it has, we move on, if it hasn't, then we add this point to the linked list in the usual way.
 ! we do this by computing the equivalent (i,j) using the Whole Pattern symmetry.
      do i=-nx,nx
       do j=-ny,ny
         if (kselected(i,j).eq.0) then
           if ((i*i+j*j).le.ijmax) then
 ! first of all, add the present point to the linked list
            call self%Add_knode(cell,Diff,ktail,i,j,(/ 0.0,0.0/))
 ! then compute the equivalent points and flag all of them in kselected
            call Apply2DPGSymmetry(TDPG,i,j,self%isym,iequiv,nequiv)
            kselected(iequiv(1,1),iequiv(2,1)) = 2
            if (nequiv.gt.1) then
             do jj=2,nequiv
              kselected(iequiv(1,jj),iequiv(2,jj)) = 1
             end do
            end if
         end if
        end if
       end do
      end do
   else                                                  ! not a zone axis, so no symmmetry
      do i=-nx,nx
       do j=-ny,ny
         if (kselected(i,j).eq.0) then
           if ((i*i+j*j).le.ijmax) then
 ! first of all, add the present point to the linked list
             call self%Add_knode(cell,Diff,ktail,i,j,sngl(Lauexy))
             kselected(i,j) = 2
         end if
        end if
       end do
      end do
   end if
 
 ! and clean up the kselected array
 deallocate(kselected)
 
 end subroutine CalckvectorsSymmetry_
 
 !--------------------------------------------------------------------------
 recursive subroutine Add_knode_(self,cell,Diff,ktail,i,j,klaue,hexgrid)
 !DEC$ ATTRIBUTES DLLEXPORT :: Add_knode_
 !! author: MDG
 !! version: 1.0
 !! date: 02/02/20
 !!
 !! add one entry to the linked wave vector list (standard mode)
 
 use mod_io
 use mod_crystallography
 use mod_diffraction
 
 IMPLICIT NONE
 
 class(kvectors_T),INTENT(INOUT)         :: self
 type(Cell_T),INTENT(INOUT)              :: cell
 type(Diffraction_T),INTENT(INOUT)       :: Diff
 type(kvectorlist),pointer               :: ktail
 integer(kind=irg),INTENT(IN)            :: i
 integer(kind=irg),INTENT(IN)            :: j
 real(kind=sgl),INTENT(IN)               :: klaue(2)
 logical,INTENT(IN),OPTIONAL             :: hexgrid
 
 type(IO_T)                              :: Message
 real(kind=sgl)                          :: kt(3),kr(3)
 real(kind=sgl)                          :: ktlen
 integer(kind=irg)                       :: istat
 
 
 allocate(ktail%next,stat=istat)                 ! allocate new value
 if (istat.ne.0) call Message%printError('Add_knode:',' unable to allocate pointer')
 ktail => ktail%next                             ! tail points to new value
 nullify(ktail%next)                             ! nullify next in new value
 self%numk = self%numk + 1                                 ! keep track of number of k-vectors so far
 ktail%i = i                                     ! i-index of beam
 ktail%j = j                                     ! j-index of beam
 
 ! is it a square or hexagonal grid ?
 if (present(hexgrid)) then
   kt = -(float(i)-float(j)*0.5)*self%delta*self%gan - float(j)*self%delta*self%gperp*0.5*sqrt(3.0)  ! tangential component of k
 else
   kt = -(klaue(1)+float(i)*self%delta)*self%gan - (klaue(2)+float(j)*self%delta)*self%gperp  ! tangential component of k
 end if
 
 ktail%kt = kt                                   ! store tangential component of k
 ktlen = cell%CalcLength(kt,'r')**2              ! squared length of tangential component
 
 kr = kt + sqrt(1.0/Diff%getWaveLength()**2 - ktlen)*self%kstar  ! complete wave vector
 ktail%k = kr                                    ! store in pointer list
 ktail%kn = cell%CalcDot(ktail%k,self%kstar,'r') ! normal component of k
 
 end subroutine Add_knode_
 
 !--------------------------------------------------------------------------
 recursive function GetSextant(x,y) result(res)
 !DEC$ ATTRIBUTES DLLEXPORT :: GetSextant
 !! author: MDG
 !! version: 1.0
 !! date: 02/02/20
 !!
 !! determines in which sextant a point (x,y) is located (used for RoscaLambert mapmode)
 
 IMPLICIT NONE
 
 real(kind=dbl),INTENT(IN):: x, y
 
 real(kind=dbl),parameter        :: srt = 1.732050808   ! sqrt(3.D0)
 integer(kind=irg)               :: res
 real(kind=dbl)                  :: xx
 
 xx = dabs(x*srt)        ! |x| sqrt(3)
 
 if (y.ge.0) then
   if (y.ge.xx) then
         res = 0
   else
         if (x.gt.0.D0) then
           res = 1
         else
           res = 5
         end if
   end if
 else
   if (dabs(y).ge.xx) then
         res = 3
   else
         if (x.gt.0.D0) then
           res = 2
         else
           res = 4
         end if
   end if
 end if
 
 end function GetSextant
 
 !--------------------------------------------------------------------------
 recursive subroutine AddkVector_(self,cell,Diff,ktail,xyval,i,j,usehex,addSH,LegendreLattitude)
 !DEC$ ATTRIBUTES DLLEXPORT :: AddkVector_
 !! author: MDG
 !! version: 1.0
 !! date: 02/02/20
 !!
 !! add a k-vector for square or hexagonal grid sampling mode (used for RoscaLambert mapmode)
 
 use mod_io
 use mod_diffraction
 use mod_crystallography
 use mod_Lambert
 
 IMPLICIT NONE
 
 class(kvectors_T),INTENT(INOUT)         :: self
 type(Cell_T),INTENT(INOUT)              :: cell
 type(Diffraction_T),INTENT(INOUT)       :: Diff
 type(kvectorlist),pointer               :: ktail
 real(kind=dbl),INTENT(IN)               :: xyval(2)
 integer(kind=irg),INTENT(IN)            :: i
 integer(kind=irg),INTENT(IN)            :: j
 logical,INTENT(IN),OPTIONAL             :: usehex
 logical,INTENT(IN),OPTIONAL             :: addSH
 real(kind=dbl),INTENT(IN),OPTIONAL      :: LegendreLattitude
 
 type(Lambert_T)                         :: Lambert
 type(IO_T)                              :: Message
 integer(kind=irg)                       :: istat, ks, ii, ierr
 real(kind=dbl)                          :: kstar(3), p
 logical                                 :: hex
 
 ! project the coordinate up to the sphere, to get a unit 3D vector kstar in cartesian space
 if (present(usehex)) then
   Lambert = Lambert_T( xyd = xyval )
   ierr = Lambert%LambertHexToSphere(kstar)
   hex = .TRUE.
 else
   Lambert = Lambert_T( xyd = xyval )
   ierr = Lambert%LambertSquareToSphere(kstar)
   hex = .FALSE.
 end if
 
 ! do we need to modify the direction cosines to coincide with the Legendre lattitudinal grid values?
 if (present(LegendreLattitude)) then
   if (kstar(3).ne.1.D0) then
 ! the factor p rescales the x and y components of kstar to maintain a unit vector
     p = sqrt((1.D0-LegendreLattitude**2)/(1.D0-kstar(3)**2))
     kstar = (/ p*kstar(1), p*kstar(2), LegendreLattitude /)
   end if
 end if
 
 ! add this vector to the linked list
      allocate(ktail%next,stat=istat)                    ! allocate new value
      if (istat.ne.0) call Message%printError('Addkvector:',' unable to allocate ktail pointer')
      ktail => ktail%next                                ! tail points to new value
      nullify(ktail%next)                                ! nullify next in new value
      self%numk = self%numk + 1                          ! keep track of number of k-vectors so far
      ktail%hs = 1                                       ! which hemisphere (Northern = 1, Southern = -1)
      ! if (iv) then                                      ! transform the hex coordinates to square-array coordinates
      !   ktail%i = i - j/2+mod(j,2)/2                     ! i-index of beam
      !   ktail%j = j                                      ! j-index of beam
      ! else                                               ! leave the square coordinates unchanged
        ktail%i = i                                      ! i-index of beam
        ktail%j = j                                      ! j-index of beam
      ! end if
      call cell%NormVec(kstar,'c')                       ! normalize incident direction in cartesian space
      kstar = kstar/Diff%getWaveLength()                 ! divide by wavelength
 ! and transform to reciprocal crystal space using the direct structure matrix
      call cell%TransSpace(kstar, ktail%k, 'c', 'r')
      ktail%kn = 1.0/Diff%getWaveLength()
 
 ! do we also need to add a Southern Hemisphere vector ?
      if (present(addSH)) then
        if (addSH.eqv..TRUE.) then
          allocate(ktail%next,stat=istat)                    ! allocate new value
          if (istat.ne.0) call Message%printError('Addkvector:',' unable to allocate ktail pointer')
          ktail => ktail%next                                ! tail points to new value
          nullify(ktail%next)                                ! nullify next in new value
          self%numk = self%numk + 1                          ! keep track of number of k-vectors so far
          ktail%hs = -1                                      ! which hemisphere (Northern = 1, Southern = -1)
          ! if (iv) then                                      ! transform the hex coordinates to square-array coordinates
          !   ktail%i = i - j/2+mod(j,2)/2                     ! i-index of beam
          !   ktail%j = j                                      ! j-index of beam
          ! else                                               ! leave the square coordinates unchanged
            ktail%i = i                                      ! i-index of beam
            ktail%j = j                                      ! j-index of beam
          ! end if
 ! get the Southern hemisphere version of kstar
          kstar(3) = -kstar(3)
 ! and transform to reciprocal crystal space using the direct structure matrix
          call cell%TransSpace(kstar, ktail%k, 'c', 'r')
          ktail%kn = 1.0/Diff%getWaveLength()
        end if
      end if
 
 end subroutine AddkVector_
 
 !--------------------------------------------------------------------------
 recursive subroutine CalckvectorsECP_(self,cell,Diff,rotmat,thetac,npx,npy,FN)
 !DEC$ ATTRIBUTES DLLEXPORT :: CalckvectorsECP_
 !! author: Saransh Singh
 !! version: 1.0
 !! date: 02/03/20
 !!
 !! create a linked list of wave vectors for conical incidence in ECP
 
 use mod_io
 use mod_diffraction
 use mod_crystallography
 use mod_Lambert
 
 IMPLICIT NONE
 
 class(kvectors_T),INTENT(INOUT)         :: self
 type(Cell_T),INTENT(INOUT)              :: cell
 type(Diffraction_T),INTENT(INOUT)       :: Diff
 real(kind=sgl),INTENT(IN)               :: rotmat(3,3)  !< initial wave vector
 real(kind=sgl),INTENT(IN)               :: thetac       !< half angle of cone of incident beam directions in degrees
 integer(kind=irg),INTENT(IN)            :: npx          !< number of kvectors along x
 integer(kind=irg),INTENT(IN)            :: npy          !< number of kvectors along y
 real(kind=sgl),INTENT(IN)               :: FN(3)        ! foil normal in reciprocal frame
 
 type(IO_T)                              :: Message
 type(kvectorlist),pointer               :: ktail
 real(kind=sgl)                          :: thetacr
 integer(kind=irg)                       :: i,j,imin,imax,jmin,jmax,ijmax,istat
 real(kind=sgl)                          :: kk(3),krec(3)
 real(kind=sgl)                          :: k(3),kcart(3),kkk
 
 if (associated(self%klist)) then   ! deallocate the entire linked list
     call self%Delete_kvectorlist()
 end if
 
 k = (/0.0,0.0,1.0/)
 k = k/Diff%getWaveLength()
 thetacr = dtor*thetac
 kk = k
 
 self%ktmax = tan(thetacr)*sqrt(sum(kk**2))
 self%delta = 2.0*self%ktmax/(2.0*float(npx)+1.0)
 imin = -npx
 imax = npx
 jmin = -npy
 jmax = npy
 
 ijmax = npx**2
 
 ! allocate the head and tail of the linked list
 allocate(self%klist,stat=istat)                         ! allocate new value
 if (istat.ne.0) call Message%printError('Calckvectors','unable to allocate self%klist pointer')
 
 ktail => self%klist                                     ! tail points to new value
 nullify(ktail%next)                                     ! nullify next in new value
 self%numk = 0
 ktail%i = 0                                             ! i-index of beam
 ktail%j = 0                                             ! j-index of beam
 ktail%kt = (/0.0,0.0,0.0/)                              ! no tangential component for central beam direction
 ktail%k = matmul(rotmat,kk)
 k = ktail%k
 call cell%TransSpace(k,krec,'c','r')
 ktail%k = krec
 
 kkk = cell%CalcDot(sngl(ktail%k),FN,'r')                ! normal component
 ktail%kn = kkk
 
 do i = imin,imax
     do j = jmin,jmax
             allocate(ktail%next)
             ktail => ktail%next
             nullify(ktail%next)
             self%numk = self%numk + 1
             ktail%i = i
             ktail%j = j
             ktail%kt = (/self%delta*i,self%delta*j,0.D0/)
             ktail%k = ktail%kt + kk
             call cell%NormVec(ktail%k,'c')
             ktail%k = ktail%k/Diff%getWaveLength()
             ktail%kt = matmul(rotmat,ktail%kt)
             ktail%k = matmul(rotmat,ktail%k)
             k = ktail%k
             call cell%TransSpace(k,krec,'c','r')
             ktail%k = krec
             k = ktail%kt
             call cell%TransSpace(k,krec,'c','r')
             ktail%kt = krec
             kkk = cell%CalcDot(sngl(ktail%k),FN,'r')
             ktail%kn = kkk
     end do
 end do
 
 end subroutine CalckvectorsECP_
 
 !--------------------------------------------------------------------------
 recursive subroutine CalckvectorsGPU_(self,cell,Diff,npx,npix,centralpix,usehex)
 !DEC$ ATTRIBUTES DLLEXPORT :: CalckvectorsGPU_
  !! author: Saransh Singh
  !! version: 1.0
  !! date: 02/03/20
  !!
  !! create a linked list of wave vectors for approximate master pattern calculation
  !!
  !! this subroutine calculates the list of k vectors in a small patch of
  !! size npix*npix in lambert space. This central beam is used to calculate the list
  !! of g vectors for the whole patch, and the k vectors are used to calculate the
  !! diagonal components i.e. the sg values.
 
 use mod_io
 use mod_diffraction
 use mod_crystallography
 use mod_Lambert
 
 IMPLICIT NONE
 
 class(kvectors_T),INTENT(INOUT)         :: self
 type(cell_T),INTENT(INOUT)              :: cell
 type(Diffraction_T),INTENT(INOUT)       :: Diff
 integer(kind=irg),INTENT(IN)            :: npx ! 2*npx+1 is size of master pattern
 integer(kind=irg),INTENT(IN)            :: npix ! the small patch will be 4*npix*npix
 integer(kind=irg),INTENT(IN)            :: centralpix(2)
 logical,INTENT(IN),OPTIONAL             :: usehex
 
 type(IO_T)                              :: Message
 type(kvectorlist),pointer               :: ktail
 logical                                 :: verbose,switchmirror,hex
 integer(kind=irg)                       :: i,j,isym,pgnum,ks,istat
 integer(kind=irg)                       :: istart,iend,jstart,jend
 integer(kind=irg),parameter             :: LaueTest(11) = (/ 149, 151, 153, 156, 158, 160, 161, 164, 165, 166, 167 /)  ! space groups with 2 or mirror at 30 degrees
 real(kind=dbl)                          :: x,y,q,XX,YY,xp,yp,rr, xy(2)
 
 
 allocate(self%klist,stat=istat)                   ! allocate new value
 if (istat.ne.0) call Message%printError('Calckvectors',' unable to allocate self%klist pointer')
 ktail => self%klist                       ! tail points to new value
 nullify(ktail%next)
 self%numk = 1
 ! we deal with the symmetry of the master pattern in the main subroutine
 istart = centralpix(1)-npix
 iend = centralpix(1)+npix-1
 jstart = centralpix(2)-npix
 jend = centralpix(2)+npix-1
 
 if (present(usehex)) then
     if (usehex) then
 ! hexagonal grid step size
         self%delta =  1.D0/dble(npx)
         ktail%i = centralpix(1)
         ktail%j = centralpix(2)
         x = (i - j*0.5)*self%delta
         y = j*self%delta*LPs%srt
         rr = x*x+y*y
         hex = .TRUE.
         ks = GetSextant(x,y)
         select case (ks)
         case (0,3)
             XX = LPs%preb*y*dcos(x*LPs%prec/y)
             YY = LPs%preb*y*dsin(x*LPs%prec/y)
         case (1,4)
             xp = y+LPs%rtt*x
             yp = y*LPs%pred/xp
             XX = LPs%prea*xp*dsin(yp)
             YY = LPs%prea*xp*dcos(yp)
         case (2,5)
             xp = y-LPs%rtt*x
             yp = y*LPs%pred/xp
             XX = LPs%prea*xp*dsin(yp)
             YY = -LPs%prea*xp*dcos(yp)
         end select
 
         q = XX**2+YY**2
         self%kstar = (/ 0.5D0*XX*dsqrt(4.D0-q), 0.5D0*YY*dsqrt(4.D0-q),1.D0-0.5D0*q /)
         ktail%i = centralpix(1) - centralpix(2)/2+mod(centralpix(2),2)/2                     ! i-index of beam
         ktail%j = centralpix(2)   ! j-index of beam
         call cell%NormVec(self%kstar,'c')
         self%kstar = self%kstar/Diff%getWaveLength()
 ! and transform to reciprocal crystal space using the direct structure matrix
         ktail%k = matmul(transpose(cell%getdsm()),self%kstar)
         ktail%kn = 1.0/Diff%getWaveLength()
     else
         self%delta = 1.D0/dble(npx)
         ktail%i = centralpix(1)
         ktail%j = centralpix(2)
         x = centralpix(1)*self%delta
         y = centralpix(2)*self%delta
         rr = x*x+y*y
         hex = .FALSE.
         if (maxval(abs((/x,y/))).eq.0.0) then
             self%kstar = (/0.D0,0.D0,1.D0/)
         else
             if (dabs(x).le.dabs(y)) then
                 q = 2.D0*y*LPs%iPi*dsqrt(cPi-y*y)
                 self%kstar = (/ q*dsin(x*LPs%Pi*0.25D0/y), q*dcos(x*LPs%Pi*0.25D0/y), 1.D0-2.D0*y*y*LPs%iPi /)
             else
                 q = 2.D0*x*LPs%iPi*dsqrt(cPi-x*x)
                 self%kstar = (/ q*dcos(y*LPs%Pi*0.25D0/x), q*dsin(y*LPs%Pi*0.25D0/x), 1.D0-2.D0*x*x*LPs%iPi /)
             end if
         end if
 
         ktail%i = centralpix(1)
         ktail%j = centralpix(2)
         call cell%NormVec(self%kstar,'c')
         self%kstar = self%kstar/Diff%getWaveLength()
 ! and transform to reciprocal crystal space using the direct structure matrix
         ktail%k = matmul(transpose(cell%getdsm()),self%kstar)
         ktail%kn = 1.0/Diff%getWaveLength()
     end if
 else
     self%delta = LPs%ap/dble(npx)
     ktail%i = centralpix(1)
     ktail%j = centralpix(2)
     x = centralpix(1)*self%delta
     y = centralpix(2)*self%delta
     rr = x*x+y*y
     hex = .FALSE.
     if (maxval(abs((/x,y/))).eq.0.0) then
         self%kstar = (/0.D0,0.D0,1.D0/)
     else
 
         if (dabs(x).le.dabs(y)) then
             q = 2.D0*y*LPs%iPi*dsqrt(cPi-y*y)
             self%kstar = (/ q*dsin(x*LPs%Pi*0.25D0/y), q*dcos(x*LPs%Pi*0.25D0/y), 1.D0-2.D0*y*y*LPs%iPi /)
         else
             q = 2.D0*x*LPs%iPi*dsqrt(cPi-x*x)
             self%kstar = (/ q*dcos(y*LPs%Pi*0.25D0/x), q*dsin(y*LPs%Pi*0.25D0/x), 1.D0-2.D0*x*x*LPs%iPi /)
         end if
     end if
     ktail%i = centralpix(1)
     ktail%j = centralpix(2)
     call cell%NormVec(self%kstar,'c')
     self%kstar = self%kstar/Diff%getWaveLength()
 ! and transform to reciprocal crystal space using the direct structure matrix
     ktail%k = matmul(transpose(cell%getdsm()),self%kstar)
     ktail%kn = 1.0/Diff%getWaveLength()
 end if
 
 do j=jstart,jend
     do i=istart,iend
         if (.not.((i .eq. centralpix(1)) .and. (j .eq. centralpix(2)))) then ! central pixel already taken care of
             if (present(usehex)) then
                 xy = (/ dble(i), dble(j) /) * self%delta
                 call self%AddkVector(cell,Diff,ktail,xy,i,j,usehex)
             else
                 xy = (/ dble(i), dble(j) /) * self%delta
                 call self%AddkVector(cell,Diff,ktail,xy,i,j)
             end if
         end if
     end do
 end do
 
 end subroutine CalckvectorsGPU_
 
 !--------------------------------------------------------------------------
!
! SUBROUTINE: Calckvectorcone
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute a set of incident beam directions for single image ECCI mode
!
!> @param cell unit cell pointer
!> @param k incident wave vector (zone axis)
!> @param ga principal g vector
!> @param ktxy tangential components
!> @param ktrad cone opening angle
!> @param ktstep number of steps along cone radius
!
!> @date 11/29/01 MDG 1.0 original
!> @date 12/05/13 MDG 2.0 adapted for ECCI simulations
!> @date 12/01/15 MDG 2.1 simplification of input parameters
!--------------------------------------------------------------------------
recursive subroutine Calckvectorcone_(self, cell, mLambda, k, ga, ktxy, ktrad, ktstep)
!DEC$ ATTRIBUTES DLLEXPORT :: Calckvectorcone_

use mod_io
use mod_diffraction
use mod_crystallography

IMPLICIT NONE

class(kvectors_T), INTENT(INOUT)    :: self
type(Cell_T),INTENT(IN)             :: cell
type(IO_T)                          :: Message
real(kind=dbl), INTENT(IN)          :: mlambda
integer(kind=irg),INTENT(IN)        :: k(3)
integer(kind=irg),INTENT(IN)        :: ga(3)
real(kind=sgl),INTENT(IN)           :: ktxy(2)
real(kind=sgl),INTENT(IN)           :: ktrad
integer(kind=irg),INTENT(IN)        :: ktstep

type(kvectorlist),pointer           :: ktmp,ktail
integer                             :: istat,imin,imax,jmin,jmax,ijmax,i,j,ic,jc,ki
real                                :: kr(3),glen,delta,kstar(3),kt(3),gan(3),gperp(3),ktlen, dkt

! compute geometrical factors
 glen = cell%CalcLength(float(ga),'r')              ! length of ga
 self%gan = ga/glen                                 ! normalized ga
 self%delta = ktrad*glen/float(ktstep)              ! grid step size in nm-1
 dkt = ktrad/float(ktstep)
 call cell%TransSpace(dble(k),self%kstar,'d','r')       ! transform incident direction to reciprocal space
 call cell%CalcCross(dble(ga),self%kstar,self%gperp,'r','r',0)! compute g_perp = ga x k
 call cell%NormVec(self%gperp,'r')                       ! normalize g_perp
 call cell%NormVec(self%kstar,'r')                       ! normalize reciprocal beam vector

! deal only with the incident beam (parallel illumination)
if (ktstep.eq.0) then
 ! if (.not.associated(self%klist)) then     ! allocate the head and ktail of the linked list
   ! allocate(self%klist,stat=istat)         ! allocate new value
   ! if (istat.ne.0) call Message%printError('Calckvectorcone: unable to allocate head pointer',' ')
   ktail => self%klist                     ! ktail points to new value
   nullify(ktail%next)                     ! nullify next in new value
   self%numk = 1                           ! keep track of number of k-vectors so far
 ! this should be the center vector of the illumination cone !!!
   kt = - glen * (ktxy(1)*self%gan + ktxy(2) * self%gperp)
   ktail%kt = kt                           ! store tangential component of k
   ktlen = glen**2*(ktxy(1)**2+ktxy(2)**2) ! squared length of tangential component

   kr = kt + sqrt(1.0/mLambda**2 - ktlen)*self%kstar ! complete wave vector
   ktail%k = kr                            ! store in pointer list
   ktail%kn = cell%CalcDot(ktail%k,dble(self%kstar),'r')    ! normal component of k
 ! end if
else
! next, put the center of the cone in units of (i,j) (original ECP "screen" coordinates)
  ic = int(ktxy(1)*glen/self%delta)
  jc = int(ktxy(2)*glen/self%delta)
  ki = ktstep

 ! if (.not.associated(self%klist)) then     ! allocate the head and ktail of the linked list
   ! allocate(self%klist,stat=istat)         ! allocate new value
   ! if (istat.ne.0) call Message%printError('Calckvectorcone: unable to allocate head pointer',' ')
   ktail => self%klist                      ! ktail points to new value
   nullify(ktail%next)                ! nullify next in new value
   self%numk = 1                          ! keep track of number of k-vectors so far
 ! this should be the center vector of the illumination cone !!!
   ktail%i = ic                            ! i-index of beam
   ktail%j = jc                            ! j-index of beam
   kt = -float(ktail%i)*self%delta*self%gan - float(ktail%j)*self%delta*self%gperp  ! tangential component of k
   ktail%kt = kt                           ! store tangential component of k
   ktlen = self%delta**2*(ktail%i**2+ktail%j**2)         ! squared length of tangential component

   kr = kt + sqrt(1.0/mLambda**2 - ktlen)*self%kstar ! complete wave vector
   ktail%k = kr                            ! store in pointer list
   ktail%kn = cell%CalcDot(ktail%k,dble(self%kstar),'r')    ! normal component of k
 ! else
   ! call Message%printError('Calckvectorcone: pointer head already allocated',' ')
 ! end if

! the following lines are quite different if symmetry is taken into account;
! check the MBsym.f90 program to determine how that can be done.
  imin =  -ki; imax = ki; jmin = -ki; jmax = ki;
  ijmax = ki**2
! now do the real work
  do i=imin,imax
   do j=jmin,jmax
    if (.not.((i.eq.0).and.(j.eq.0))) then  ! the point (0,0) has already been taken care of
     if ((i**2+j**2).le.ijmax) then   ! is point inside the incident cone ?
      allocate(ktail%next,stat=istat)  ! allocate new value
      if (istat.ne.0) call Message%printError('Calckvectorcone: unable to allocate pointer',' ')
      ktail => ktail%next               ! ktail points to new value
      nullify(ktail%next)              ! nullify next in new value
      self%numk = self%numk + 1                 ! keep track of number of k-vectors so far
      ktail%i = ic+i                   ! i-index of beam
      ktail%j = jc+j                   ! j-index of beam
      kt = - float(ktail%i)*self%delta*self%gan - float(ktail%j)*self%delta*self%gperp  ! tangential component of k
      ktail%kt = kt                    ! store tangential component of k
      ! ktlen = self%delta**2*(ktail%i**2+ktail%j**2)         ! squared length of tangential component
      ktlen = cell%CalcLength(kt,'r')**2

      kr = kt + sqrt(1.0/mLambda**2 - ktlen)*self%kstar ! complete wave vector
      ktail%k = kr                     ! store in pointer list
      ktail%kn = cell%CalcDot(ktail%k,dble(self%kstar),'r')    ! normal component of k
     end if
    end if
   end do
  end do
end if

end subroutine Calckvectorcone_



!--------------------------------------------------------------------------
!
! SUBROUTINE: Calckvectortrace
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute a set of incident beam directions for line scan ECCI mode
!
!> @param cell unit cell pointer
!> @param k incident wave vector (zone axis)
!> @param ga principal g vector
!> @param ktxy tangential components of trace start point
!> @param ktxy2 tangential components of trace end point
!> @param ktrad cone opening angle
!> @param ktstep number of steps along cone radius
!
!> @date 12/08/13 MDG 1.0 original
!> @date 12/01/15 MDG 1.1 simplifcation of input variables
!--------------------------------------------------------------------------
recursive subroutine Calckvectortrace_(self,cell,mlambda,k,ga,ktxy,ktxy2,ktrad,ktstep)
!DEC$ ATTRIBUTES DLLEXPORT :: Calckvectortrace_

use mod_io
use mod_diffraction
use mod_crystallography

IMPLICIT NONE

class(kvectors_T), INTENT(INOUT)    :: self
type(Cell_T),INTENT(IN)             :: cell
type(IO_T)                          :: Message
real(kind=dbl), INTENT(IN)          :: mlambda
integer(kind=irg),INTENT(IN)        :: k(3)
integer(kind=irg),INTENT(IN)        :: ga(3)
real(kind=sgl),INTENT(IN)           :: ktxy(2)
real(kind=sgl),INTENT(IN)           :: ktxy2(2)
real(kind=sgl),INTENT(IN)           :: ktrad
integer(kind=irg),INTENT(IN)        :: ktstep

type(kvectorlist),pointer           :: ktail
integer                             :: istat,j
real                                :: kr(3),glen,delta,kstar(3),kt(3),gan(3),gperp(3),ktlen, dktx, dkty

! compute geometrical factors
 glen = cell%CalcLength(float(ga),'r')              ! length of ga
 self%gan = ga/glen                                 ! normalized ga
 self%delta = 2.0*ktrad*glen/float(2*ktstep+1)      ! grid step size in nm-1
 call cell%TransSpace(dble(k),self%kstar,'d','r')       ! transform incident direction to reciprocal space
 call cell%CalcCross(dble(ga),self%kstar,self%gperp,'r','r',0)! compute g_perp = ga x k
 call cell%NormVec(self%gperp,'r')                       ! normalize g_perp
 call cell%NormVec(self%kstar,'r')                       ! normalize reciprocal beam vector

 j = 0
 if (.not.associated(self%klist)) then     ! allocate the head and ktail of the linked list
   allocate(self%klist,stat=istat)         ! allocate new value
   if (istat.ne.0) call Message%printError('Calckvectortrace: unable to allocate head pointer',' ')
   ktail => self%klist                     ! ktail points to new value
   nullify(ktail%next)                ! nullify next in new value
   self%numk = 1                           ! keep track of number of k-vectors so far
! this should be the starting point of the line trace
!   kt = - glen * ( ktxy(1)*gan + ktxy(2) * gperp)
   kt = - glen * ( ktxy(1)*self%gan - ktxy(2) * self%gperp)
   ktail%kt = kt                           ! store tangential component of k
   ktlen = glen**2*(ktxy(1)**2+ktxy(2)**2)         ! squared length of tangential component
   kr = kt + sqrt(1.0/mLambda**2 - ktlen)*self%kstar ! complete wave vector
   ktail%k = kr                            ! store in pointer list
   ktail%kn = cell%CalcDot(ktail%k,dble(self%kstar),'r')    ! normal component of k
 end if

 dktx = (ktxy2(1) - ktxy(1))/float(ktstep-1)
 dkty = (ktxy2(2) - ktxy(2))/float(ktstep-1)

 do j=1,ktstep-1
      allocate(ktail%next,stat=istat)  ! allocate new value
      if (istat.ne.0) call Message%printError('Calckvectortrace: unable to allocate pointer',' ')
      ktail => ktail%next              ! ktail points to new value
      nullify(ktail%next)              ! nullify next in new value
      self%numk = self%numk + 1                  ! keep track of number of k-vectors so far
!     kt = - glen * ( (ktxy(1)+float(j)*dktx)*gan + (ktxy(2)+float(j)*dkty) * gperp) ! tangential component of k
      kt = - glen * ( (ktxy(1)+float(j)*dktx)*self%gan - (ktxy(2)+float(j)*dkty) * self%gperp) ! tangential component of k
      ktail%kt = kt                    ! store tangential component of k
      ktlen = self%delta**2*(ktail%i**2+ktail%j**2)         ! squared length of tangential component
      kr = kt + sqrt(1.0/mLambda**2 - ktlen)*self%kstar ! complete wave vector
      ktail%k = kr                     ! store in pointer list
      ktail%kn = cell%CalcDot(ktail%k,dble(self%kstar),'r')    ! normal component of k
 end do

end subroutine Calckvectortrace_


!--------------------------------------------------------------------------
!
! SUBROUTINE: Calckvectorcone
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute a set of incident beam directions for single image ECCI mode
!
!> @param cell unit cell pointer
!> @param khead head of kvector linked list
!> @param k incident wave vector (zone axis)
!> @param ga principal g vector
!> @param ktxy tangential components
!> @param ktrad cone opening angle
!> @param ktstep number of steps along cone radius
!> @param numk resulting number of incident beam directions
!
!> @date 11/29/01 MDG 1.0 original
!> @date 12/05/13 MDG 2.0 adapted for ECCI simulations
!> @date 12/01/15 MDG 2.1 simplification of input parameters
!--------------------------------------------------------------------------
recursive subroutine Calckvectorcircle_(self, cell, mLambda, k, ga, ktxy, ktrad, ktstep)
!DEC$ ATTRIBUTES DLLEXPORT :: Calckvectorcircle_

use mod_io
use mod_diffraction
use mod_crystallography

IMPLICIT NONE

class(kvectors_T), INTENT(INOUT)    :: self
type(Cell_T),INTENT(IN)             :: cell
type(IO_T)                          :: Message
real(kind=dbl), INTENT(IN)          :: mlambda
integer(kind=irg),INTENT(IN)        :: k(3)
integer(kind=irg),INTENT(IN)        :: ga(3)
real(kind=sgl),INTENT(IN)           :: ktxy(2)
real(kind=sgl),INTENT(IN)           :: ktrad
integer(kind=irg),INTENT(IN)        :: ktstep

type(kvectorlist),pointer           :: ktmp,ktail
integer                             :: istat,imin,imax,jmin,jmax,ijmax,i,j,ic,jc, iang
real                                :: kr(3),glen,delta, dtang, kstar(3),kt(3),gan(3),gperp(3),ktlen, dkt, ii, &
                                       jj, ii1, ii2, jj1, jj2
real(kind=dbl)                      :: ki, iangrad
! compute geometrical factors
 glen = cell%CalcLength(float(ga),'r')         ! length of ga
 self%gan = ga/glen                                 ! normalized ga
 self%delta = ktrad*glen/float(ktstep)              ! grid step size in nm-1
 dkt = ktrad/float(ktstep)
 call cell%TransSpace(dble(k),self%kstar,'d','r')       ! transform incident direction to reciprocal space
 call cell%CalcCross(dble(ga),self%kstar,self%gperp,'r','r',0)! compute g_perp = ga x k
 call cell%NormVec(self%gperp,'r')                       ! normalize g_perp
 call cell%NormVec(self%kstar,'r')                       ! normalize reciprocal beam vector

! deal only with the incident beam (parallel illumination)
if (ktstep.eq.0) then
 if (.not.associated(self%klist)) then     ! allocate the head and ktail of the linked list
   allocate(self%klist,stat=istat)         ! allocate new value
   if (istat.ne.0) call Message%printError('Calckvectorcone: unable to allocate head pointer',' ')
   ktail => self%klist                      ! ktail points to new value
   nullify(ktail%next)                ! nullify next in new value
   self%numk = 1                          ! keep track of number of k-vectors so far
 ! this should be the center vector of the illumination cone !!!
   kt = - glen * (ktxy(1)*self%gan + ktxy(2) * self%gperp)
   ktail%kt = kt                           ! store tangential component of k
   ktlen = glen**2*(ktxy(1)**2+ktxy(2)**2)         ! squared length of tangential component

   kr = kt + sqrt(1.0/mLambda**2 - ktlen)*self%kstar ! complete wave vector
   ktail%k = kr                            ! store in pointer list
   ktail%kn = cell%CalcDot(ktail%k,dble(self%kstar),'r')    ! normal component of k
 end if
else
! next, put the center of the cone in units of (i,j) (original ECP "screen" coordinates)
  ic = int(ktxy(1)*glen/self%delta)
  jc = int(ktxy(2)*glen/self%delta)

  imin =  -ki; imax = ki; jmin = -ki; jmax = ki;
  ijmax = (ki)**2

  if (.not.associated(self%klist)) then 
    allocate(self%klist,stat=istat)
    ktail => self%klist
  end if

  imax = 360
  ki = dble(ktrad)*dtor
  self%numk = 0

  iangrad = 1.0*dtor
  ii1 = ki*(dcos(iangrad))
  jj1 = ki*(dsin(iangrad))
  ii1 = dble(ii1)*rtod
  jj1 = dble(jj1)*rtod
  iangrad = 3.0*dtor
  ii2 = ki*(dcos(iangrad))
  jj2 = ki*(dsin(iangrad))
  ii2 = dble(ii2)*rtod
  jj2 = dble(jj2)*rtod
  dtang = sqrt(((ii2-ii1)**2+(jj2-jj1)**2))
  
  do iang = 1, 360, 1
    iangrad = iang*dtor
    ii = ki*(dcos(iangrad))
    jj = ki*(dsin(iangrad))

    ii = dble(ii)*rtod
    jj = dble(jj)*rtod

    if (.not.((ii.eq.0).and.(jj.eq.0))) then  ! the point (0,0) has already been taken care of
      allocate(ktail%next,stat=istat)  ! allocate new value
      if (istat.ne.0) call Message%printError('Calckvectorcone: unable to allocate pointer',' ')
      
      self%numk = self%numk + 1                 ! keep track of number of k-vectors so far
      kt = - (ii)*self%gan*self%delta- (jj)*self%gperp*self%delta  ! tangential component of k
      ktail%kt = kt                    ! store tangential component of k
      ktlen = self%delta**2*(ii**2+jj**2)         ! squared length of tangential component

      kr = kt + sqrt(1.0/mLambda**2 - ktlen)*self%kstar ! complete wave vector
      ktail%k = kr                     ! store in pointer list
      ktail%kn = cell%CalcDot(ktail%k,dble(self%kstar),'r')
      ktail => ktail%next               ! ktail points to new value
      nullify(ktail%next)              ! nullify next in new value
    end if
  end do
 end if

end subroutine Calckvectorcircle_



 end module mod_kvectors
 
