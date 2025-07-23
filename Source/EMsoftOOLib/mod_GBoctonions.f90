! ###################################################################
! Copyright (c) 2013-2025, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_GBoctonions
  !! author: MDG 
  !! version: 1.0 
  !! date: 07/16/25
  !!
  !! class definition for the Grain Boundary Octonions module
  !!
  !! there is a stand alone Octonion_T class and an OctonionArray_T class
  !! and we inherit those from the mod_octonions module
  !!
  !! What turns an octonion into a GBoctonion is the fact that it is made up
  !! of two unit quaternions, so it requires a normalization factor of sqrt(2) 
  !!
  !! see the following publication for details:
  !!
  !! T. Francis, I. Chesser, S. Singh, E.A. Holm and M. De Graef. 
  !! "A Geodesic Octonion Metric for Grain Boundaries". 
  !! Acta Materialia, 166:135-147 (2019)
  !! DOI: https://doi.org/10.1016/j.actamat.2018.12.034

use mod_kinds
use mod_global
use mod_octonions
use mod_quaternions

IMPLICIT NONE 

! If we use this module, which extends the Octonion_T class, then by definition
! we are using grain boundary octonions, which are constructed from two unit quaternions.
! Therefore, they use an octonion normalization by a factor of sqrt(2), but this is 
! handled transparently by the parent Octonion_T class.
! Other than that, there really are not many differences between the two
! classes.  There are of course grain boundary specific operations which are 
! defined in this module.

! This module also implements the findings of Prof. Oliver Johnson from BYU.
! In particular, his model introduces the concept of "simulated inversion" in 
! addition to a correct grain exchange symmetry operation (the one in the octonion paper
! is actually incorrect).  


! class definition for the Grain Boundary Octonion
type, public, extends(Octonion_T) :: GBoctonion_T
private 

contains
private 

   procedure, pass(self) :: GBO_get_q_
   procedure, pass(self) :: GBO_get_equivalent_
   procedure, pass(self) :: GBO_Omega_
   procedure, pass(self) :: GBO_Omega_NB_
   procedure, pass(self) :: GBO_SLERP_
   procedure, pass(self) :: GBO_Omega_Refine_
   procedure, pass(self) :: GBO_Omega_Symmetric_
   procedure, pass(self) :: GBO_Omega_symmetric_NB_

   generic, public :: GBO_get_q => GBO_get_q_
   generic, public :: GBO_get_equivalent => GBO_get_equivalent_
   generic, public :: GBO_Omega => GBO_Omega_
   generic, public :: GBO_Omega_NB => GBO_Omega_NB_
   generic, public :: GBO_SLERP => GBO_SLERP_
   generic, public :: GBO_Omega_Refine => GBO_Omega_Refine_
   generic, public :: GBO_Omega_Symmetric => GBO_Omega_Symmetric_
   generic, public :: GBO_Omega_Symmetric_NB => GBO_Omega_Symmetric_NB_

end type GBoctonion_T

! class definition for the Grain Boundary Octonion Array
type, public, extends(OctonionArray_T) :: GBOctonionArray_T
private

contains
private

   procedure, pass(self) :: insertGBOctintoArray_

   generic, public :: insertGBOctinArray => insertGBOctintoArray_

end type GBOctonionArray_T

private:: insertGBOctintoArray_, GBO_Omega_, GBO_Omega_NB_, GBO_SLERP_
public :: GBO_minimal_U1_angle, GBO_minimize_U1_angle

! the constructor routines for these classes 
interface GBoctonion_T
  module procedure GBoctonion_constructor
end interface GBoctonion_T

interface GBoctonionArray_T
  module procedure GBOctonionArray_constructor
end interface GBoctonionArray_T

contains

!--------------------------------------------------------------------------
type(GBoctonion_T) function GBoctonion_constructor( qu1, qu2, oct ) result(GBoctonion)
!DEC$ ATTRIBUTES DLLEXPORT :: GBoctonion_constructor
!! author: MDG 
!! version: 1.0 
!! date: 10/16/22
!!
!! constructor for the GBoctonions_T Class
 
use mod_quaternions

IMPLICIT NONE

type(Quaternion_T), INTENT(IN),OPTIONAL    :: qu1
type(Quaternion_T), INTENT(IN),OPTIONAL    :: qu2
type(Octonion_T),INTENT(IN),OPTIONAL       :: oct

if (present(qu1)) then 
  if (qu1%quat_getprecision().eq.'s') then
    GBoctonion%o = (/ qu1%get_quats(), qu2%get_quats() /)
    GBoctonion%s = 's'
  else 
    GBoctonion%od = (/ qu1%get_quatd(), qu2%get_quatd() /)
    GBoctonion%s = 'd'
  end if 
end if 

if (present(oct)) then 
  if (oct%s.eq.'s') then 
    GBoctonion%o = oct%get_octs()
    GBoctonion%s = 's'
  else
    GBoctonion%od = oct%get_octd()
    GBoctonion%s = 'd'
  end if 
end if

! this normalization involves sqrt(2) due to the two unit quaternions, but this 
! is correctly handled by the parent class Octonion_T
call GBoctonion%o_normalize()

end function GBoctonion_constructor

!--------------------------------------------------------------------------
subroutine GBOctonion_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 07/17/25
!!
!! destructor for the GBoctonion_T Class
 
IMPLICIT NONE

type(GBoctonion_T), INTENT(INOUT)  :: self 

call reportDestructor('GBoctonion_T')

end subroutine GBOctonion_destructor

!--------------------------------------------------------------------------
type(GBOctonionArray_T) function GBOctonionArray_constructor( qAr1, qAr2, n, nthreads, s ) result(OctArray)
!DEC$ ATTRIBUTES DLLEXPORT :: GBOctonionArray_constructor
  !! author: MDG
  !! version: 1.0
  !! date: 07/17/25
  !!
  !! constructor for the GBOctonionArray Class
  !!
  !! this constructor takes two QuaternionArrays and merges them into a GBOctonionArray
  !! Alternatively, it initializes the arrays to a requested size

use mod_io 

IMPLICIT NONE

  type(QuaternionArray_T),INTENT(INOUT),OPTIONAL      :: qAr1 
  type(QuaternionArray_T),INTENT(INOUT),OPTIONAL      :: qAr2 
  integer(kind=irg),INTENT(IN),OPTIONAL               :: n 
  integer(kind=irg),INTENT(IN),OPTIONAL               :: nthreads
  character(1),INTENT(IN),OPTIONAL                    :: s 

  type(IO_T)                                          :: Message 

  integer(kind=irg)                                   :: i 
  type(Quaternion_T)                                  :: q1, q2
  type(GBoctonion_T)                                  :: gboct

if (present(qAr1)) then
! make sure the arrays have the same size
  if (qAr1%getQnumber().ne.qAr2%getQnumber()) then 
    call Message%printError('GBOctonionArray_constructor',' input quaternion arrays have different size')
  end if 

! inherit quaternion array parameters
  OctArray%nthreads = qAr1%getnthreads()
  OctArray%n = qAr1%getQnumber()
  OctArray%s = qAr1%getprecision()
  
! allocate the GBO array
  if (OctArray%s.eq.'s') then
    if (allocated(OctArray%o)) deallocate(OctArray%o)
    allocate( OctArray%o(8,OctArray%n) ) 
  else
    if (allocated(OctArray%od)) deallocate(OctArray%od)
    allocate( OctArray%od(8,OctArray%n) ) 
  end if 

  do i=1,OctArray%n
    q1 = qAr1%getQuatfromArray(i)
    q2 = qAr2%getQuatfromArray(i)
    gboct = GBoctonion_T( q1, q2 )
    call OctArray%insertGBOctinArray(i, gboct)
  end do
else
! set array parameters
  OctArray%nthreads = nthreads
  OctArray%n = n
  OctArray%s = s

  ! allocate the GBO array
  if (s.eq.'s') then
    if (allocated(OctArray%o)) deallocate(OctArray%o)
    allocate( OctArray%o(8,OctArray%n) ) 
    Octarray%o = 0.0
    Octarray%o(1,:) = 1.0
  else
    if (allocated(OctArray%od)) deallocate(OctArray%od)
    allocate( OctArray%od(8,OctArray%n) ) 
    Octarray%od = 0.D0
    Octarray%od(1,:) = 1.D0
  end if 
end if 

end function GBOctonionArray_constructor

!--------------------------------------------------------------------------
subroutine GBOctonionArray_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: GBOctonionArray_destructor
!! author: MDG
!! version: 1.0
!! date: 07/17/25
!!
!! destructor for the GBOctonionArray_T Class

IMPLICIT NONE

type(GBOctonionArray_T), INTENT(INOUT)     :: self

call reportDestructor('GBOctonionArray_T')

if (allocated(self%o)) deallocate(self%o)
if (allocated(self%od)) deallocate(self%od)

end subroutine GBOctonionArray_destructor

!--------------------------------------------------------------------------
recursive subroutine insertGBOctintoArray_(self, i, o)
!DEC$ ATTRIBUTES DLLEXPORT :: insertGBOctintoArray_
  !! author: MDG
  !! version: 1.0
  !! date: 07/16/25
  !!
  !! insert a GBoctonion in an existing array (overrides mod_octonions)

use mod_io 

IMPLICIT NONE

class(GBOctonionArray_T),INTENT(INOUT)  :: self
integer(kind=irg),INTENT(IN)            :: i
type(GBOctonion_T),INTENT(INOUT)        :: o

type(IO_T)                              :: Message

! make sure that the index i is within the appropriate range 
if (i.gt.self%n) call Message%printError('insertGBOctintoArray_',' index too large for octonion array')

if (self%s.eq.'s') then 
  self%o(1:8,i) = o%get_octs()
else
  self%od(1:8,i) = o%get_octd()
end if

end subroutine insertGBOctintoArray_

!--------------------------------------------------------------------------
recursive function GBO_get_q_(self, n)  result(qu)
!DEC$ ATTRIBUTES DLLEXPORT :: GBO_get_q_
!! author: MDG
!! version: 1.0
!! date: 07/19/25
!!
!! extracts quaternion n from a grain boundary octonion 

class(GBoctonion_T),INTENT(INOUT)                 :: self   ! this is oct1 
integer(kind=irg),INTENT(IN)                      :: n
type(Quaternion_T)                                :: qu

if (self%s.eq.'s') then 
  if (n.eq.1) then 
    qu = Quaternion_T( q = self%o(1:4) )
  else
    qu = Quaternion_T( q = self%o(5:8) )
  end if 
else
  if (n.eq.1) then 
    qu = Quaternion_T( qd = self%od(1:4) )
   else
    qu = Quaternion_T( qd = self%od(5:8) )
  end if 
end if 

! make sure that this is a properly normalized quaternion 
call qu%quat_normalize()

end function GBO_get_q_

!--------------------------------------------------------------------------
recursive function GBO_Omega_symmetric_(self, oct2, DS, solution, arclengths, single, noU1, &
                                       metric, refine)  result(Omega)
!DEC$ ATTRIBUTES DLLEXPORT :: GBO_Omega_symmetric_
!! author: MDG
!! version: 1.0
!! date: 07/16/25
!!
!! Compute the S^7 geodesic arc length for a GBO pair (U(1) symmetry, grain exchange and crystal symmetry)

use mod_dirstats
use mod_IO

IMPLICIT NONE

class(GBoctonion_T),INTENT(INOUT)                 :: self   ! this is oct1 
type(GBoctonion_T),INTENT(INOUT)                  :: oct2
type(DirStat_T),INTENT(INOUT)                     :: DS
type(GBOctonionArray_T),INTENT(OUT),OPTIONAL      :: solution
real(kind=dbl),allocatable,INTENT(OUT),OPTIONAL   :: arclengths(:,:)  ! (/ Nqsym**2, Nqsym**2 /)
logical,INTENT(IN),OPTIONAL                       :: single
logical,INTENT(IN),OPTIONAL                       :: noU1
character(fnlen),INTENT(IN),OPTIONAL              :: metric
logical,INTENT(IN),OPTIONAL                       :: refine
real(kind=dbl)                                    :: Omega

type(GBoctonion_T)                                :: GBab, GBcd 
type(QuaternionArray_T)                           :: qsym, qAr
type(Quaternion_T)                                :: Sqa, Sqb, Sqc, Sqd, qa, qb, qc, qd, qu
type(IO_T)                                        :: Message

integer(kind=irg)                                 :: i, j, k, l, Nqsym, io_int(1), ss
logical                                           :: keep, arcs, skipU1, dorefine
real(kind=dbl)                                    :: smallest, x
character(fnlen)                                  :: usemetric

! handle the optional input parameters
usemetric = 'octonion'
if (present(metric)) usemetric = trim(metric)
! call Message%printMessage(' Metric to be used for computation : '//trim(usemetric))

dorefine = .FALSE. 
if (present(refine)) then 
  if (refine.eqv..TRUE.) dorefine=.TRUE.
end if 

skipU1 = .FALSE.
if (present(noU1)) then
  if (noU1.eqv..TRUE.) skipU1 = .TRUE.
end if

keep = .FALSE.
if (present(solution)) then
  if (keep.eqv..TRUE.) then
    keep = .TRUE.
    qAr = QuaternionArray_T( n = 2, s = self%s )
    solution = GBOctonionArray_T( qAr, qAr )
  end if 
end if 

qsym = DS%getQuatArray(slot='qsym')
Nqsym = qsym%getQnumber()
io_int(1) = Nqsym
! call Message%WriteValue(' Number of symmetry operators generated : ', io_int,1)
arcs = .FALSE.
if (present(arclengths)) then
  if (arcs.eqv..TRUE.) then
    arcs = .TRUE.
    allocate( arclengths(Nqsym**2, Nqsym**2) )
  end if 
end if 
smallest = 100.D0

if (Nqsym.eq.1) then
  if (skipU1.eqv..TRUE.) then
    smallest = self%GBO_Omega_(oct2,noU1=.TRUE.,metric=usemetric)
  else
    if (dorefine.eqv..TRUE.) then
      smallest = self%GBO_Omega_Refine_(oct2,metric=usemetric)
    else
      smallest = self%GBO_Omega_(oct2,metric=usemetric)
    end if
  end if
  if (keep.eqv..TRUE.) then 
    call solution%insertGBOctinArray(1, self)
    call solution%insertGBOctinArray(2, oct2)
  end if 
  if (arcs.eqv..TRUE.) then
    arclengths(1,1) = smallest
  end if
else
  if (present(single)) then 
    if (single.eqv..TRUE.) then
      do k=1,Nqsym
        qc = oct2%GBO_get_q(1)
        Sqc = qsym%getQuatfromArray(k) * qc
        call Sqc%quat_pos()
        do l=1,Nqsym
          qd = oct2%GBO_get_q(2)
          Sqd = qsym%getQuatfromArray(l) * qd
          call Sqc%quat_pos()
          GBcd = GBoctonion_T( Sqc, Sqd )
          if (skipU1.eqv..TRUE.) then  
            x = self%GBO_Omega_( GBcd, noU1=.TRUE., metric=usemetric)
          else
            if (dorefine.eqv..TRUE.) then 
              if ((k+l).eq.2) then    ! first time we need to initialize some arrays
                x = GBab%GBO_Omega_Refine(GBcd,metric=usemetric,init=.TRUE.)
              end if 
              x = self%GBO_Omega_Refine_( GBcd, metric=usemetric)
            else
              x = self%GBO_Omega_( GBcd, metric=usemetric)
            end if
          end if
          if (arcs.eqv..TRUE.) then 
            arclengths((i-1)*Nqsym+j,(k-1)*Nqsym+l) = x
          end if
          if (x.lt.smallest) then 
            smallest = x
            if (keep) then 
              call solution%insertGBOctinArray(1, self)
              call solution%insertGBOctinArray(2, GBcd)
            end if 
          end if
        end do
      end do
    end if 
  else
    qa = self%GBO_get_q(1)
    qb = self%GBO_get_q(2)
    qc = oct2%GBO_get_q(1)
    qd = oct2%GBO_get_q(2)
    ! write (*,*) 'working with the following quaternions:'
    ! call qa%quat_print(' qa : ')
    ! call qb%quat_print(' qb : ')
    ! call qc%quat_print(' qc : ')
    ! call qd%quat_print(' qd : ')
    ss = 0
    do i=1,Nqsym
      Sqa = qsym%getQuatfromArray(i) * qa
      do j=1,Nqsym
        Sqb = qsym%getQuatfromArray(j) * qb
        GBab = GBoctonion_T( Sqa, Sqb )
        do k=1,Nqsym
          Sqc = qsym%getQuatfromArray(k) * qc
          do l=1,Nqsym
            ss = ss+1
            Sqd = qsym%getQuatfromArray(l) * qd
            GBcd = GBoctonion_T( Sqc, Sqd )
            if (skipU1.eqv..TRUE.) then  
              x = GBab%GBO_Omega(GBcd,noU1=.TRUE.,metric=usemetric)
            else
              if (dorefine.eqv..TRUE.) then 
                if ((i+j+k+l).eq.4) then    ! first time we need to initialize some arrays
                  x = GBab%GBO_Omega_Refine(GBcd,metric=usemetric,init=.TRUE.)
                end if 
                x = GBab%GBO_Omega_Refine(GBcd,metric=usemetric)
              else
                x = GBab%GBO_Omega(GBcd,metric=usemetric)
              end if
            end if
            if (arcs) then 
              arclengths((i-1)*Nqsym+j,(k-1)*Nqsym+l) = x
            end if
            if (x.lt.smallest) then 
              smallest = x
              if (keep) then 
                call solution%insertGBOctinArray(1, GBab)
                call solution%insertGBOctinArray(2, GBcd)
              end if 
            end if
          end do
        end do
      end do
    end do
  end if 
end if 

Omega = smallest

end function GBO_Omega_symmetric_

!--------------------------------------------------------------------------
recursive function GBO_Omega_symmetric_NB_(self, oct2, Nqsym, qsym)  result(Omega)
!DEC$ ATTRIBUTES DLLEXPORT :: GBO_Omega_symmetric_NB_

use mod_quaternions

IMPLICIT NONE

class(GBoctonion_T),INTENT(INOUT) :: self
type(GBoctonion_T),INTENT(INOUT)  :: oct2
integer(kind=irg),INTENT(IN)      :: Nqsym 
type(QuaternionArray_T),INTENT(IN):: qsym
real(kind=dbl)                    :: Omega

integer(kind=irg)                 :: i, k
real(kind=dbl)                    :: smallest, x
type(Quaternion_T)                :: qa, qc, Sqa, Sqc

qa = self%GBO_get_q_(1)
qc = oct2%GBO_get_q_(1)

if (Nqsym.eq.1) then
  smallest = self%GBO_Omega_NB_(qa,qc)
else
  smallest = 1000.D0
  do i=1,Nqsym
    Sqa = qsym%getQuatfromArray(i) * qa
    call Sqa%quat_pos()
    do k=1,Nqsym
      Sqc = qsym%getQuatfromArray(k) * qc
      call Sqc%quat_pos()
      x = self%GBO_Omega_NB_(Sqa,Sqc)
      if (x.lt.smallest) smallest = x
    end do
  end do
end if 

Omega = smallest

end function GBO_Omega_symmetric_NB_

!--------------------------------------------------------------------------
recursive function GBO_Omega_(self,oct2,metric,noU1)  result(Omega)
!DEC$ ATTRIBUTES DLLEXPORT :: GBO_Omega_
!! author: MDG
!! version: 1.0
!! date: 07/16/25
!!
!! Compute the angle that will minimize a geodesic quaternion arc length with respect to U(1) symmetry

IMPLICIT NONE

class(GBoctonion_T),INTENT(INOUT) :: self 
type(GBoctonion_T),INTENT(INOUT)  :: oct2
character(fnlen),INTENT(IN)       :: metric
logical,INTENT(IN),OPTIONAL       :: noU1
real(kind=dbl)                    :: Omega

type(Quaternion_T)                :: qu

real(kind=dbl)                    :: qq1(4), qq2(4), zeta, sigma, cac, cbd, cbc, cad, cz, sz, cs, ss, &
                                     sum1, sum2, sum3, sum4, sums(4), smax, qa(4), qb(4), qc(4), qd(4)
integer(kind=irg)                 :: isum(1), m
real(kind=dbl),parameter          :: srt = 1.D0/sqrt(2.D0)     

! extract the quaternions as regular 4-component arrays
qu = self%GBO_get_q(1)
qa = qu%get_quatd()
qu = self%GBO_get_q(2)
qb = qu%get_quatd()
qu = oct2%GBO_get_q(1)
qc = qu%get_quatd()
qu = oct2%GBO_get_q(2)
qd = qu%get_quatd()

m = 1
if (trim(metric).eq.'Olmsted') then 
  m=2
else if (trim(metric).eq.'Riemannian') then 
       m=3
     end if

if (present(noU1)) then
  if (noU1.eqv..TRUE.) then
    cac = sum(qa*qc)
    cbd = sum(qb*qd)
    cbc = sum(qb*qc)
    cad = sum(qa*qd)
   
    select case(m)
    case(1) 
      sum1 = 0.5D0 * maxval( abs( (/ cac+cbd, cac-cbd /) ) )
      sum2 = 0.5D0 * maxval( abs( (/ cbc+cad, cbc-cad /) ) )
      Omega = 2.0 * minval( (/ acos(sum1), acos(sum2) /) )
    case(2)
      sum1 = sqrt(4.0D0 * ( 2.D0 - cac*cac - cbd*cbd ))
      sum2 = sqrt(4.0D0 * ( 2.D0 - cbc*cbc - cad*cad ))
      Omega = minval( (/ sum1, sum2 /) )
    case(3)
      cac = 2.D0 * acos(cac)
      cbd = 2.D0 * acos(cbd)
      cbc = 2.D0 * acos(cbc)
      cad = 2.D0 * acos(cad)
      sum1 = sqrt(cac*cac + cbd*cbd )
      sum2 = sqrt(cbc*cbc + cad*cad )
      Omega = minval( (/ sum1, sum2 /) )
    end select

  end if
else
! determine the minimal U(1) angle for the (a,b) - (c,d) boundary pair
  zeta = GBO_minimal_U1_angle(qa,qb,qc,qd)
  ! write (*,*) ' (a,b) - (c,d) ',zeta/dtor
  cz = cos(zeta*0.5D0)
  sz = sin(zeta*0.5D0)
  qq1 = (/ qc(1)*cz-qc(4)*sz, cz*qc(2)+sz*qc(3), cz*qc(3)-sz*qc(2), cz*qc(4)+sz*qc(1) /)
  qq2 = (/ qd(1)*cz-qd(4)*sz, cz*qd(2)+sz*qd(3), cz*qd(3)-sz*qd(2), cz*qd(4)+sz*qd(1) /)
  cac = sum(qa*qq1)
  cbd = sum(qb*qq2)
  
  select case(m)
    case(1) 
      sum1 = 0.5D0 * maxval( abs( (/ cac+cbd, cac-cbd /) ) )
    case(2)
      sum1 = sqrt(4.0D0 * ( 2.D0 - cac*cac - cbd*cbd ))
    case(3)
      cac = 2.D0 * acos(cac)
      cbd = 2.D0 * acos(cbd)
      sum1 = sqrt(cac*cac + cbd*cbd)
  end select

! determine the minimal U(1) angle for the (a,-b) - (c,d) boundary pair
  zeta = GBO_minimal_U1_angle(qa,-qb,qc,qd)
  ! write (*,*) ' (a,-b) - (c,d) ',zeta/dtor
  cz = cos(zeta*0.5D0)
  sz = sin(zeta*0.5D0)
  qq1 = (/ qc(1)*cz-qc(4)*sz, cz*qc(2)+sz*qc(3), cz*qc(3)-sz*qc(2), cz*qc(4)+sz*qc(1) /)
  qq2 = (/ qd(1)*cz-qd(4)*sz, cz*qd(2)+sz*qd(3), cz*qd(3)-sz*qd(2), cz*qd(4)+sz*qd(1) /)
  cac = sum(qa*qq1)
  cbd = sum(-qb*qq2)
  
  select case(m)
    case(1) 
      sum3 = 0.5D0 * maxval( abs( (/ cac+cbd, cac-cbd /) ) )
    case(2)
      sum3 = sqrt(4.0D0 * ( 2.D0 - cac*cac - cbd*cbd ))
    case(3)
      cac = 2.D0 * acos(cac)
      cbd = 2.D0 * acos(cbd)
      sum3 = sqrt(cac*cac + cbd*cbd)
  end select

! determine the minimal U(1) angle for the (b,a) - (c,d) boundary pair
  sigma = GBO_minimal_U1_angle(qa,qb,qc,qd,exchange=.TRUE.)
  ! write (*,*) ' (b, a) - (c,d) ',sigma/dtor
  cs = cos(sigma*0.5D0)
  ss = sin(sigma*0.5D0)
  qq1 = (/ qc(1)*cs-qc(4)*ss, cs*qc(2)+ss*qc(3), cs*qc(3)-ss*qc(2), cs*qc(4)+ss*qc(1) /)
  qq2 = (/ qd(1)*cs-qd(4)*ss, cs*qd(2)+ss*qd(3), cs*qd(3)-ss*qd(2), cs*qd(4)+ss*qd(1) /)
  cbc = sum(qb*qq1)
  cad = sum(qa*qq2)

  select case(m)
    case(1) 
      sum2 = 0.5D0 * maxval( abs( (/ cbc+cad, cbc-cad /) ) )
    case(2)
      sum2 = sqrt(4.0D0 * ( 2.D0 - cbc*cbc - cad*cad ))
    case(3)
      cbc = 2.D0 * acos(cbc)
      cad = 2.D0 * acos(cad)
      sum2 = sqrt(cbc*cbc + cad*cad)
  end select

! determine the minimal U(1) angle for the (b,-a) - (c,d) boundary pair
  sigma = GBO_minimal_U1_angle(-qa,qb,qc,qd,exchange=.TRUE.)
  ! write (*,*) ' (b,-a) - (c,d) ',sigma/dtor
  cs = cos(sigma*0.5D0)
  ss = sin(sigma*0.5D0)
  qq1 = (/ qc(1)*cs-qc(4)*ss, cs*qc(2)+ss*qc(3), cs*qc(3)-ss*qc(2), cs*qc(4)+ss*qc(1) /)
  qq2 = (/ qd(1)*cs-qd(4)*ss, cs*qd(2)+ss*qd(3), cs*qd(3)-ss*qd(2), cs*qd(4)+ss*qd(1) /)
  cbc = sum(qb*qq1)
  cad = sum(-qa*qq2)

  select case(m)
    case(1) 
      sum4 = 0.5D0 * maxval( abs( (/ cbc+cad, cbc-cad /) ) )
    case(2)
      sum4 = sqrt(4.0D0 * ( 2.D0 - cbc*cbc - cad*cad ))
    case(3)
      cbc = 2.D0 * acos(cbc)
      cad = 2.D0 * acos(cad)
      sum4 = sqrt(cbc*cbc + cad*cad)
  end select

  sums = (/ sum1, sum2, sum3, sum4 /)
  smax = maxval(sums)
  isum = maxloc(sums)

  ! write (*,*) ' sums : ', sums 

! and determine the smallest geodesic distance on S^7
  select case(m)
    case(1)
      Omega = 2.0 * acos(smax)
    case(2,3)
      Omega = minval(sums)
      Omega = Omega * srt
  end select

end if

end function GBO_Omega_

!--------------------------------------------------------------------------
recursive function GBO_Omega_NB_(self, qa, qc)  result(Omega)
!DEC$ ATTRIBUTES DLLEXPORT :: GBO_Omega_NB_

use mod_quaternions

IMPLICIT NONE

class(GBoctonion_T),INTENT(IN)    :: self
type(Quaternion_T),INTENT(INOUT)  :: qa
type(Quaternion_T),INTENT(INOUT)  :: qc
real(kind=dbl)                    :: Omega

type(Quaternion_T)                :: qq, pp 
real(kind=dbl)                    :: q(4), p(4)

qq = qa * conjg(qc)
pp = qc * conjg(qa)
q = qq%get_quatd()
p = pp%get_quatd()
Omega = minval( (/ 2.D0 * acos(abs(q)), 2.D0*acos(abs(p)) /) )

end function GBO_Omega_NB_

!--------------------------------------------------------------------------
recursive function GBO_minimal_U1_angle(qa,qb,qc,qd,exchange)  result(zeta)
!DEC$ ATTRIBUTES DLLEXPORT :: GBO_minimal_U1_angle
!! author: MDG
!! version: 1.0
!! date: 07/16/25
!!
!! Compute the angle that will minimize a geodesic quaternion arc length with respect to U(1) symmetry
!!
!! this was rewritten with Oliver Johnson's new solution; in particular the grain exchange
!! expression.

use mod_rotations
use mod_quaternions

IMPLICIT NONE

real(kind=dbl),INTENT(IN)         :: qa(4)
real(kind=dbl),INTENT(IN)         :: qb(4)
real(kind=dbl),INTENT(IN)         :: qc(4)
real(kind=dbl),INTENT(IN)         :: qd(4)
logical,INTENT(IN),OPTIONAL       :: exchange
real(kind=dbl)                    :: zeta

type(Quaternion_T)                :: qpi_x 
type(Quaternion_T)                :: qua_x, qub_x

real(kind=dbl)                    :: nom, denom, mu, v1, v4, qax(4), qbx(4)

zeta = 0.D0

if (present(exchange)) then 
  if (exchange.eqv..TRUE.) then 
    qpi_x = Quaternion_T( qd = (/ 0.D0, 1.D0, 0.D0, 0.D0 /) ) 
    qua_x = Quaternion_T( qd = qb ) * qpi_x
    qub_x = Quaternion_T( qd = qa ) * qpi_x
    qax = qua_x%get_quatd()
    qbx = qub_x%get_quatd()
    v1 = sum(qax*qc) + sum(qbx*qd)
    v4 = (qax(4)*qc(1)-qax(1)*qc(4)) - (qax(3)*qc(2)-qax(2)*qc(3)) + &
         (qbx(4)*qd(1)-qbx(1)*qd(4)) - (qbx(3)*qd(2)-qbx(2)*qd(3))
    zeta = 2.D0 * atan2(v4, v1)
  end if
else
! using Oliver Johnson's new results
  v1 = sum(qa*qc) + sum(qb*qd)
  v4 = (qa(4)*qc(1)-qa(1)*qc(4)) - (qa(3)*qc(2)-qa(2)*qc(3)) + &
       (qb(4)*qd(1)-qb(1)*qd(4)) - (qb(3)*qd(2)-qb(2)*qd(3))
  zeta = 2.D0 * atan2(v4, v1)
end if 

if (zeta.lt.0.D0) zeta = 4.D0*cPi + zeta

end function GBO_minimal_U1_angle

!--------------------------------------------------------------------------
recursive function GBO_minimize_U1_angle(qa,qb,qc,qd,numz,z,czs,szs,m,exchange)  result(zval)
!DEC$ ATTRIBUTES DLLEXPORT :: GBO_minimize_U1_angle
!! author: MDG
!! version: 1.0
!! date: 07/16/25
!!

IMPLICIT NONE

real(kind=dbl),INTENT(IN)         :: qa(4)
real(kind=dbl),INTENT(IN)         :: qb(4)
real(kind=dbl),INTENT(IN)         :: qc(4)
real(kind=dbl),INTENT(IN)         :: qd(4)
integer(kind=irg),INTENT(IN)      :: numz 
real(kind=dbl),INTENT(IN)         :: z(numz)
real(kind=dbl),INTENT(IN)         :: czs(numz)
real(kind=dbl),INTENT(IN)         :: szs(numz)
integer(kind=irg),INTENT(IN)      :: m
logical,INTENT(IN),OPTIONAL       :: exchange
real(kind=dbl)                    :: zval

real(kind=dbl)                    :: nom, denom, mu, qq1(4,numz), qq2(4,numz), cac(numz), cbd(numz), &
                                     sm(numz), mval, x1, x2, x3, y1, y2, y3, A, B, C
integer(kind=irg)                 :: i, mpos(1)

zval = 0.D0

qq1(1,:) = czs(:)*qc(1)-szs(:)*qc(4)
qq1(2,:) = czs(:)*qc(2)+szs(:)*qc(3)
qq1(3,:) = czs(:)*qc(3)-szs(:)*qc(2)
qq1(4,:) = czs(:)*qc(4)+szs(:)*qc(1)

qq2(1,:) = czs(:)*qd(1)-szs(:)*qd(4)
qq2(2,:) = czs(:)*qd(2)+szs(:)*qd(3)
qq2(3,:) = czs(:)*qd(3)-szs(:)*qd(2)
qq2(4,:) = czs(:)*qd(4)+szs(:)*qd(1)

if (present(exchange)) then 
  if (exchange.eqv..TRUE.) then 
    do i=1,numz
      cac(i) = sum(qa(:)*qq1(:,i))
      cbd(i) = sum(qb(:)*qq2(:,i))
    end do
  else
    do i=1,numz
      cac(i) = sum(qb(:)*qq1(:,i))
      cbd(i) = sum(qa(:)*qq2(:,i))
    end do
  end if 
end if 

if (m.eq.2) then 
  sm = sqrt(4.0D0 * ( 2.D0 - cac*cac - cbd*cbd ))
else if (m.eq.3) then
        cac = 2.D0 * acos(cac)
        cbd = 2.D0 * acos(cbd)
        sm = sqrt(cac*cac + cbd*cbd)
     end if

mval = minval(sm)
mpos = minloc(sm)

if ((mpos(1).ne.1).and.(mpos(1).ne.numz)) then 
  x1 = z(mpos(1)-1)
  x2 = z(mpos(1))
  x3 = z(mpos(1)+1)

  y1 = sm(mpos(1)-1)
  y2 = sm(mpos(1))
  y3 = sm(mpos(1)+1)
else if (mpos(1).eq.1) then 
        x1 = z(numz)
        x2 = z(1)
        x3 = z(2)

        y1 = sm(numz)
        y2 = sm(1)
        y3 = sm(2)
     else 
        x1 = z(numz-1)
        x2 = z(numz)
        x3 = z(1)

        y1 = sm(numz-1)
        y2 = sm(numz)
        y3 = sm(1)
     end if

! simply fit a parabola through three points and determine the location of the minimum.
! denom = (x1 - x2) * (x1 - x3) * (x2 - x3)
A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) ! / denom;
B = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) !  / denom;
! we don't need the value at the minimum, so no need to compute C
! C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

zval = -B / (2.D0*A)

end function GBO_minimize_U1_angle

!--------------------------------------------------------------------------
recursive function GBO_Omega_Refine_(self ,oct2,metric,init)  result(Omega)
!DEC$ ATTRIBUTES DLLEXPORT :: GBO_Omega_Refine_
!! author: MDG
!! version: 1.0
!! date: 07/16/25
!!


IMPLICIT NONE

class(GBoctonion_T),INTENT(INOUT) :: self 
type(GBoctonion_T),INTENT(INOUT)  :: oct2
character(fnlen),INTENT(IN)       :: metric
logical,INTENT(IN),OPTIONAL       :: init 
real(kind=dbl)                    :: Omega

type(Quaternion_T)                :: qu 

real(kind=dbl)                    :: qq1(4), qq2(4), zeta, sigma, cac, cbd, cbc, cad, cz, sz, cs, ss, &
                                     sum1, sum2, sum3, sum4, sums(4), smax, qa(4), qb(4), qc(4), qd(4)
integer(kind=irg)                 :: isum(1), m, i
real(kind=dbl),parameter          :: srt = 1.D0/sqrt(2.D0)     
integer(kind=irg),parameter       :: numz = 180
real(kind=dbl),save               :: czs(numz), szs(numz), z(numz)


! extract the quaternions as regular 4-component arrays
qu = self%GBO_get_q(1)
qa = qu%get_quatd()
qu = self%GBO_get_q(2)
qb = qu%get_quatd()
qu = oct2%GBO_get_q(1)
qc = qu%get_quatd()
qu = oct2%GBO_get_q(2)
qd = qu%get_quatd()

if (present(init)) then 
  if (init.eqv..TRUE.) then 
    z = (/ (i-1, i=1,numz) /) * 4.D0 * cPi / dble(numz) 
    czs = cos( z * 0.5D0 ) 
    szs = sin( z * 0.5D0 ) 
    Omega = 0.D0
    return
  end if 
end if 

m = 1
if (trim(metric).eq.'Olmsted') then 
  m=2
else if (trim(metric).eq.'Riemannian') then 
       m=3
     end if

! determine the minimal U(1) angle for the (a,b) - (c,d) boundary pair
  if (m.eq.1) then 
    zeta = GBO_minimal_U1_angle(qa,qb,qc,qd)
  else 
    zeta = GBO_minimize_U1_angle(qa,qb,qc,qd,numz,z,czs,szs,m)
  end if
  cz = cos(zeta*0.5D0)
  sz = sin(zeta*0.5D0)
  qq1 = (/ qc(1)*cz-qc(4)*sz, cz*qc(2)+sz*qc(3), cz*qc(3)-sz*qc(2), cz*qc(4)+sz*qc(1) /)
  qq2 = (/ qd(1)*cz-qd(4)*sz, cz*qd(2)+sz*qd(3), cz*qd(3)-sz*qd(2), cz*qd(4)+sz*qd(1) /)
  cac = sum(qa*qq1)
  cbd = sum(qb*qq2)

  select case(m)
    case(1) 
      sum1 = 0.5D0 * maxval( abs( (/ cac+cbd, cac-cbd /) ) )
    case(2)
      sum1 = sqrt(4.0D0 * ( 2.D0 - cac*cac - cbd*cbd ))
    case(3)
      cac = 2.D0 * acos(cac)
      cbd = 2.D0 * acos(cbd)
      sum1 = sqrt(cac*cac + cbd*cbd)
  end select

! determine the minimal U(1) angle for the (a,-b) - (c,d) boundary pair
  if (m.eq.1) then 
    zeta = GBO_minimal_U1_angle(qa,-qb,qc,qd)
  else 
    zeta = GBO_minimize_U1_angle(qa,-qb,qc,qd,numz,z,czs,szs,m)
  end if
  cz = cos(zeta*0.5D0)
  sz = sin(zeta*0.5D0)
  qq1 = (/ qc(1)*cz-qc(4)*sz, cz*qc(2)+sz*qc(3), cz*qc(3)-sz*qc(2), cz*qc(4)+sz*qc(1) /)
  qq2 = (/ qd(1)*cz-qd(4)*sz, cz*qd(2)+sz*qd(3), cz*qd(3)-sz*qd(2), cz*qd(4)+sz*qd(1) /)
  cac = sum(qa*qq1)
  cbd = sum(-qb*qq2)
  
  select case(m)
    case(1) 
      sum3 = 0.5D0 * maxval( abs( (/ cac+cbd, cac-cbd /) ) )
    case(2)
      sum3 = sqrt(4.0D0 * ( 2.D0 - cac*cac - cbd*cbd ))
    case(3)
      cac = 2.D0 * acos(cac)
      cbd = 2.D0 * acos(cbd)
      sum3 = sqrt(cac*cac + cbd*cbd)
  end select

! determine the minimal U(1) angle for the (b,a) - (c,d) boundary pair
 if (m.eq.1) then 
    sigma = GBO_minimal_U1_angle(qa,qb,qc,qd,exchange=.TRUE.)
  else 
    sigma = GBO_minimize_U1_angle(qa,qb,qc,qd,numz,z,czs,szs,m,exchange=.TRUE.)
  end if
  cs = cos(sigma*0.5D0)
  ss = sin(sigma*0.5D0)
  qq1 = (/ qc(1)*cs-qc(4)*ss, cs*qc(2)+ss*qc(3), cs*qc(3)-ss*qc(2), cs*qc(4)+ss*qc(1) /)
  qq2 = (/ qd(1)*cs-qd(4)*ss, cs*qd(2)+ss*qd(3), cs*qd(3)-ss*qd(2), cs*qd(4)+ss*qd(1) /)
  cbc = sum(qb*qq1)
  cad = sum(qa*qq2)

  select case(m)
    case(1) 
      sum2 = 0.5D0 * maxval( abs( (/ cbc+cad, cbc-cad /) ) )
    case(2)
      sum2 = sqrt(4.0D0 * ( 2.D0 - cbc*cbc - cad*cad ))
    case(3)
      cbc = 2.D0 * acos(cbc)
      cad = 2.D0 * acos(cad)
      sum2 = sqrt(cbc*cbc + cad*cad)
  end select

! determine the minimal U(1) angle for the (b,-a) - (c,d) boundary pair
 if (m.eq.1) then 
    sigma = GBO_minimal_U1_angle(-qa,qb,qc,qd,exchange=.TRUE.)
  else 
    sigma = GBO_minimize_U1_angle(-qa,qb,qc,qd,numz,z,czs,szs,m,exchange=.TRUE.)
  end if
  cs = cos(sigma*0.5D0)
  ss = sin(sigma*0.5D0)
  qq1 = (/ qc(1)*cs-qc(4)*ss, cs*qc(2)+ss*qc(3), cs*qc(3)-ss*qc(2), cs*qc(4)+ss*qc(1) /)
  qq2 = (/ qd(1)*cs-qd(4)*ss, cs*qd(2)+ss*qd(3), cs*qd(3)-ss*qd(2), cs*qd(4)+ss*qd(1) /)
  cbc = sum(qb*qq1)
  cad = sum(-qa*qq2)
 
  select case(m)
    case(1) 
      sum4 = 0.5D0 * maxval( abs( (/ cbc+cad, cbc-cad /) ) )
    case(2)
      sum4 = sqrt(4.0D0 * ( 2.D0 - cbc*cbc - cad*cad ))
    case(3)
      cbc = 2.D0 * acos(cbc)
      cad = 2.D0 * acos(cad)
      sum4 = sqrt(cbc*cbc + cad*cad)
  end select

  sums = (/ sum1, sum2, sum3, sum4 /)
  smax = maxval(sums)
  isum = maxloc(sums)

! and determine the smallest geodesic distance on S^7
  select case(m)
    case(1)
      Omega = 2.0 * acos(smax)
    case(2,3)
      Omega = minval(sums)
      Omega = Omega * srt
  end select

end function GBO_Omega_Refine_

!--------------------------------------------------------------------------
recursive function GBO_SLERP_(self, hcn1, hcn2, Omega, t, n)  result(hcnt)
!DEC$ ATTRIBUTES DLLEXPORT :: GBO_SLERP_
!! author: MDG
!! version: 1.0
!! date: 07/16/25
!!


IMPLICIT NONE

class(GBoctonion_T),INTENT(INOUT)     :: self
integer(kind=irg),INTENT(IN)          :: n
real(kind=dbl),INTENT(IN)             :: hcn1(n)
real(kind=dbl),INTENT(IN)             :: hcn2(n)
real(kind=dbl),INTENT(IN)             :: Omega
real(kind=dbl),INTENT(IN)             :: t
real(kind=dbl)                        :: hcnt(n)

real(kind=dbl)                        :: st, sp, sm, theta
 
theta = Omega*0.5D0

st = sin(theta)
sp = sin(t*theta)
sm = sin((1.D0-t)*theta)

hcnt = hcn1 * sm/st + hcn2 * sp/st

end function GBO_SLERP_

!--------------------------------------------------------------------------
recursive function GBO_get_equivalent_(self, qsym, nthreads, NBflag)  result(GBO_equiv)
!DEC$ ATTRIBUTES DLLEXPORT :: GBO_get_equivalent_
!! author: MDG
!! version: 1.0
!! date: 07/23/25
!!


IMPLICIT NONE

class(GBoctonion_T),INTENT(INOUT)     :: self
type(QuaternionArray_T),INTENT(INOUT) :: qsym
integer(kind=irg),INTENT(IN)          :: nthreads
integer(kind=irg),INTENT(INOUT)       :: NBflag     ! no boundary flag 
type(GBOctonionArray_T)               :: GBO_equiv

type(GBOctonionArray_T)               :: GBO_temp
type(Quaternion_T)                    :: qa, qb, Sqa, Sqb, qpi_x, qua_x, qub_x 
type(GBoctonion_T)                    :: GBab
type(Octonion_T)                      :: o

integer(kind=irg)                     :: icnt, j, k, l, Nqsym, Nlist 
real(kind=dbl)                        :: diff, epsd=1.0D-12, octo(8) 
real(kind=dbl),allocatable            :: olist(:,:)
logical                               :: newoct

qa = self%GBO_get_q(1)
qb = self%GBO_get_q(2)

Nqsym = qsym%getQnumber()
GBO_temp = GBOctonionArray_T( n = 2*Nqsym**2, s='d', nthreads = nthreads )
qpi_x = Quaternion_T( qd = (/ 0.D0, 1.D0, 0.D0, 0.D0 /) ) 

icnt = 0
NBflag = 0
do k=1,Nqsym
  Sqa = qsym%getQuatfromArray(k) * qa
  call Sqa%quat_pos()
  do l=1,Nqsym
    Sqb = qsym%getQuatfromArray(l) * qb
    call Sqb%quat_pos()
! this is regular crystallographic symmetry
    GBab = GBoctonion_T( Sqa, Sqb )
    diff = sum( abs(Sqa%get_quatd() - Sqb%get_quatd()) )
    icnt = icnt + 1
    if (diff.lt.epsd) NBflag = icnt
    call GBO_temp%insertGBOctintoArray_(icnt, GBab)
! next we do grain exchange symmetry
    qua_x = Sqb * qpi_x
    qub_x = Sqa * qpi_x
    GBab = GBoctonion_T( qua_x, qub_x )
    diff = sum( abs(qua_x%get_quatd() - qub_x%get_quatd()) )
    icnt = icnt + 1
    if (diff.lt.epsd) NBflag = icnt
    call GBO_temp%insertGBOctintoArray_(icnt, GBab)
  end do 
end do

if (NBflag.eq.0) then 
! next we need to determine how many unique octonions there are and return
! only those... that requires a two-step process...  for now we just return 
! the temp class
  Nlist = 1
  allocate( olist(8, icnt) )
  o = GBO_temp%getOctfromArray(1)
  olist(1:8,1) = o%get_octd()
  do k=2,icnt 
    o = GBO_temp%getOctfromArray(k)
    octo = o%get_octd()
    newoct = .TRUE.
    do j=1,k-1
      diff = sum( abs( olist(1:8,j) - octo(1:8) ) )
      if (diff.lt.epsd) then 
        newoct = .FALSE.
        EXIT 
      end if 
    end do 
    if (newoct.eqv..TRUE.) then 
      Nlist = Nlist + 1
      olist(1:8,Nlist) = octo(1:8)
    end if 
  end do

  if (Nlist.eq.icnt) then ! simply return the original list since they are all unique
    GBO_equiv = GBO_temp
  else  ! make a new (shorter) list
    GBO_equiv = GBOctonionArray_T( n = Nlist, s='d', nthreads = nthreads )
    do k = 1, Nlist
      o = Octonion_T( od = olist(1:8,k) )
      GBab = GBoctonion_T( oct = o  )
      call GBO_equiv%insertGBOctintoArray_(k, GBab)
    end do 
  end if 
else
    GBO_equiv = GBO_temp
end if 

end function GBO_get_equivalent_


end module mod_GBoctonions