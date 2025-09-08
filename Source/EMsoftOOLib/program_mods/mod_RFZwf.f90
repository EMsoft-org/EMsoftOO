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

module mod_RFZwf
  !! author: MDG 
  !! version: 1.0 
  !! date: 09/08/25
  !!
  !! class definition for the EMRFZwf program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMRFZwf program
type, public :: RFZwfNameListType
  character(fnlen)  :: hdfname
end type RFZwfNameListType

! class definition
type, public :: RFZwf_T
private 
  character(fnlen)                :: nmldeffile = 'EMRFZwf.nml'
  type(RFZwfNameListType)         :: nml 
  integer(kind=irg),dimension(10) :: nLaue = (/ 3, 16, 9, 21, 6, 18, 12, 24, 28, 30 /)
  character(3),dimension(10)      :: Laue = (/ '2  ', '3  ', '4  ', '6  ', '222', &
                                               '32 ', '422', '622', '23 ', '432' /)
  ! group names in the HDF output file
  character(13)                   :: ortype(5) = (/ 'Cubochoric   ','Homochoric   ','Stereographic', &
                                                    'Rodrigues    ','Euler        '/)

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: RFZwf_
  ! procedure, pass(self) :: doCubochoric_
  ! procedure, pass(self) :: doHomochoric_
  ! procedure, pass(self) :: doStereographic_
  ! procedure, pass(self) :: doRodrigues_
  ! procedure, pass(self) :: doEuler_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: RFZwf => RFZwf_
  ! generic, public :: doCubochoric => doCubochoric_
  ! generic, public :: doHomochoric => doHomochoric_
  ! generic, public :: doStereographic => doStereographic_
  ! generic, public :: doRodrigues => doRodrigues_
  ! generic, public :: doEuler => doEuler_

end type RFZwf_T

! the constructor routine for this class 
interface RFZwf_T
  module procedure RFZwf_constructor
end interface RFZwf_T

contains

!--------------------------------------------------------------------------
type(RFZwf_T) function RFZwf_constructor( nmlfile ) result(RFZwf)
!! author: MDG 
!! version: 1.0 
!! date: 09/08/25
!!
!! constructor for the RFZwf_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call RFZwf%readNameList(nmlfile)

end function RFZwf_constructor

!--------------------------------------------------------------------------
subroutine RFZwf_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 09/08/25
!!
!! destructor for the RFZwf_T Class
 
IMPLICIT NONE

type(RFZwf_T), INTENT(INOUT)  :: self 

call reportDestructor('RFZwf_T')

end subroutine RFZwf_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 09/08/25
!!
!! read the namelist from an nml file for the RFZwf_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(RFZwf_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

character(fnlen)  :: hdfname

namelist / EMRFZwf / hdfname 

hdfname = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
  open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
  read(UNIT=dataunit,NML=EMRFZwf)
  close(UNIT=dataunit,STATUS='keep')

  ! if (trim(hdfname).eq.'undefined') then
  !     call Message%printError('readNameList:',' hdfname is undefined in '//nmlfile)
  ! end if
end if 

self%nml%hdfname = hdfname 

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 09/08/25
!!
!! pass the namelist for the RFZwf_T Class to the calling program

IMPLICIT NONE 

class(RFZwf_T), INTENT(INOUT)          :: self
type(RFZwfNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine do_(self, HDF)
!DEC$ ATTRIBUTES DLLEXPORT :: do_
!! author: MDG 
!! version: 1.0 
!! date: 09/08/25
!!
!! handle the   representation

use HDF5
use mod_HDFsupport

class(RFZwf_T), INTENT(INOUT)   :: self
type(HDF_T),INTENT(INOUT)       :: HDF



end subroutine do_

!--------------------------------------------------------------------------
subroutine RFZwf_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: RFZwf_
!! author: MDG 
!! version: 1.0 
!! date: 09/08/25
!!
!! make lists of all the vertices needed to drawn an RFZ in any of 5 orientation representations.
!!
!! this routine calls code that is similar to that in the mod_povray module, but without
!! actually calling any of the povray commands...

use mod_EMsoft
use mod_rotations
use mod_quaternions
use mod_io
use HDF5
use mod_HDFsupport

IMPLICIT NONE 

class(RFZwf_T), INTENT(INOUT)       :: self
type(EMsoft_T), INTENT(INOUT)       :: EMsoft
character(fnlen), INTENT(INOUT)     :: progname 

type(HDF_T)                         :: HDF
type(IO_T)                          :: Message

real(kind=dbl),allocatable          :: ropos(:,:), sppos(:,:), cupos(:,:), hopos(:,:), eupos(:,:)

integer(kind=irg)                   :: FZorder, sz(2), i, hdferr, iL, FZcyclic(4)
character(fnlen)                    :: datafile, dataset, groupname


FZcyclic = (/ 2, 3, 4, 6 /)

call openFortranHDFInterface()

! open the output HDF5 file
HDF = HDF_T()

datafile = EMsoft%generateFilePath('EMdatapathname', self%nml%hdfname)

hdferr =  HDF%createFile(datafile)
if (hdferr.ne.0) call HDF%error_check('HDF_createFile ', hdferr)

! create all the orientation representation groups
do i=1,5
  groupname = trim(self%ortype(i))
  hdferr = HDF%createGroup(groupname)
  call HDF%pop()
end do 

! each one of these groups will contain 10 datasets that represent the 
! entire wireframe for each of the Laue groups in each of the representations.

! loop over the cyclic Laue groups
do i=1,4
  FZorder = FZcyclic(i)
  if (allocated(ropos)) deallocate(ropos)
! get the Rodrigues wireframe and store it in the Rodrigues group
  call initFZCyclic_(FZorder, ropos)
  sz = shape(ropos)
  groupname = 'Rodrigues'
  hdferr = HDF%openGroup(groupname)
  dataset = 'Laue_'//trim(self%Laue(i))
  hdferr = HDF%writeDatasetDoubleArray(dataset, ropos, sz(1), sz(2))
  if (hdferr.ne.0) call HDF%error_check('writeDatasetDoubleArray ropos', hdferr)
  call HDF%pop()

! then do the other representations
  call convert_orep_(ropos, sz, sppos, cupos, hopos, eupos)

! Stereographic
  groupname = 'Stereographic'
  hdferr = HDF%openGroup(groupname)
  dataset = 'Laue_'//trim(self%Laue(i))
  hdferr = HDF%writeDatasetDoubleArray(dataset, sppos, sz(1), sz(2))
  if (hdferr.ne.0) call HDF%error_check('writeDatasetDoubleArray sppos', hdferr)
  deallocate(sppos)
  call HDF%pop()

! Cubochoric
  groupname = 'Cubochoric'
  hdferr = HDF%openGroup(groupname)
  dataset = 'Laue_'//trim(self%Laue(i))
  hdferr = HDF%writeDatasetDoubleArray(dataset, cupos, sz(1), sz(2))
  if (hdferr.ne.0) call HDF%error_check('writeDatasetDoubleArray cupos', hdferr)
  deallocate(cupos)
  call HDF%pop()

! Homochoric
  groupname = 'Homochoric'
  hdferr = HDF%openGroup(groupname)
  dataset = 'Laue_'//trim(self%Laue(i))
  hdferr = HDF%writeDatasetDoubleArray(dataset, hopos, sz(1), sz(2))
  if (hdferr.ne.0) call HDF%error_check('writeDatasetDoubleArray hopos', hdferr)
  deallocate(hopos)
  call HDF%pop()

! Euler
  groupname = 'Euler'
  hdferr = HDF%openGroup(groupname)
  dataset = 'Laue_'//trim(self%Laue(i))
  hdferr = HDF%writeDatasetDoubleArray(dataset, eupos, sz(1), sz(2))
  if (hdferr.ne.0) call HDF%error_check('writeDatasetDoubleArray eupos', hdferr)
  deallocate(eupos)
  call HDF%pop()
end do 



call HDF%popall()


call closeFortranHDFInterface()

end subroutine RFZwf_



!--------------------------------------------------------------------------
recursive subroutine convert_orep_(ropos, sz, sppos, cupos, hopos, eupos)
!DEC$ ATTRIBUTES DLLEXPORT :: convert_orep_
!! author: MDG
!! version: 1.0
!! date: 09/08/25
!!
!! convert the input Rodrigues vector arrays to the other representations

use mod_rotations
use mod_io

IMPLICIT NONE

integer(kind=irg),INTENT(IN)          :: sz(2)
real(kind=dbl),INTENT(IN)             :: ropos(sz(1),sz(2))
real(kind=dbl),INTENT(OUT),allocatable:: sppos(:,:)
real(kind=dbl),INTENT(OUT),allocatable:: cupos(:,:)
real(kind=dbl),INTENT(OUT),allocatable:: hopos(:,:)
real(kind=dbl),INTENT(OUT),allocatable:: eupos(:,:)

type(r_T)                             :: r
type(s_T)                             :: s
type(c_T)                             :: c
type(h_T)                             :: h
type(e_T)                             :: e

integer(kind=irg)                     :: i 
real(kind=dbl)                        :: x

call setRotationPrecision('d')

allocate( sppos(sz(1),sz(2)), cupos(sz(1),sz(2)), hopos(sz(1),sz(2)), eupos(sz(1),sz(2)))

sppos = 0.D0 
cupos = 0.D0 
hopos = 0.D0 
eupos = 0.D0 

do i=1,sz(2)
  if (sum(abs(ropos(1:3,i))).ne.0.D0) then 
    x = sqrt(sum(ropos(1:3,i)**2))
    r = r_T( rdinp = (/ ropos(1:3,i)/x, x /) )
    s = r%rs()
    sppos(1:3,i) = s%s_copyd() 
    c = r%rc()
    cupos(1:3,i) = c%c_copyd() 
    h = r%rh()
    hopos(1:3,i) = h%h_copyd() 
    e = r%re()
    eupos(1:3,i) = e%e_copyd() 
  end if 
end do 

end subroutine convert_orep_

!--------------------------------------------------------------------------
recursive subroutine initFZCyclic_(FZorder, ropos)
!DEC$ ATTRIBUTES DLLEXPORT :: initFZCyclic_
!! author: MDG
!! version: 1.0
!! date: 09/08/25
!!
!! generate the coordinates of the wireframe for the cyclic rotational groups
!! these are stored as 3-component Rodrigues vectors, with (0,0,0) entries
!! separating the major line segments.

use mod_rotations
use mod_io

IMPLICIT NONE

integer(kind=irg),INTENT(IN)          :: FZorder
real(kind=dbl),INTENT(OUT),allocatable:: ropos(:,:)

type(e_T)                             :: eul, eu, euld, eulast
type(r_T)                             :: ro1, ro2, ro, rolast, ron
type(q_T)                             :: qu
type(s_T)                             :: sp, splast
type(h_T)                             :: h, ho, holast
type(o_T)                             :: om
type(c_T)                             :: cu, culast
type(a_T)                             :: axang
type(orientation_T)                   :: ot
type(IO_T)                            :: Message

real(kind=dbl)                        :: rmax, dx, r, xmax, x, y, z, zsmall, ac, sh(3), xx, &
                                         tpi, hpi, aux(4), aux4a(4), aux4b(4)

integer(kind=irg)                     :: i,j,k, icnt, imax, nt, ns, idpos, icpos, ihedge
integer(kind=irg),allocatable         :: h_edge(:,:)
real(kind=dbl),allocatable            :: cpos(:,:), dpos(:)
! parameters that depend on the cyclic group
real(kind=dbl)                        :: a, b, c, dt, ds, d, dd, zz, oo, c2, tmp

select case(FZorder)
  case(2) ! define the coordinates of the monoclinic C2 (2) FZ in Rodrigues Space
    a = 570.289922125538D0
    b = 1.0D0
    c = 1.0D0
    dt = 114.57984425107713D0
    ds = 2.0D0
    d = 1.7320508075688772D0
    zz = 0.D0
    oo = 1.D0
    idpos = 200
    icpos = 200
    ihedge = 100
    allocate(cpos(3,icpos), h_edge(2,ihedge), dpos(idpos) )

    do i=-12,12
      if (abs(i).ne.12) then
        cpos(1:3, 13+i) = (/ -a, dtan(dble(i)*15.D0*dtor*0.5D0) ,  c /)
      else
        cpos(1:3, 13+i) = (/ -a, a ,  c /)
        if (i.lt.0) cpos(2,13+i) = -cpos(2,13+i)
      end if
    end do

    do i=1,25
      cpos(1:3,25+i) = cpos(1:3,i)
      cpos(1,25+i) = -cpos(1,25+i)
    end do

    do i=1,50
      cpos(1:3,50+i) = cpos(1:3,i)
      tmp = cpos(1,50+i)
      cpos(1,50+i) = cpos(2,50+i)
      cpos(2,50+i) = tmp
    end do

    do i=1,100
      cpos(1:3,100+i) = cpos(1:3,i)
      cpos(3,100+i) = -cpos(3,100+i)
    end do

    ! and normalize
    do i=1,200
      dpos(i) = dsqrt(sum(cpos(1:3,i)*cpos(1:3,i)))
    end do

! this FZ must be rotated so that the two-fold axis falls along the monoclinic b-axis.
    do i=1,200
       tmp = cpos(2,i)
       cpos(2,i) = cpos(3,i)
       cpos(3,i) = tmp
    end do

    ns = 2000
    dx = dt/float(ns-1)

    ! define the connectivity of all the edges
    do i=1,25
      h_edge(1:2,i)    = (/    i, 25+i /)
      h_edge(1:2,25+i) = (/ 50+i, 75+i /)
      h_edge(1:2,50+i) = (/100+i,125+i /)
      h_edge(1:2,75+i) = (/150+i,175+i /)
    end do
  case(3) ! define the coordinates of the trigonal C3 (3) FZ in Rodrigues Space
    a = 57.289922125538D0
    b = 1.0D0
    c = 0.577350269120D0
    dt = 114.57984425107713D0
    ds = 2.0D0
    d = 1.7320508075688772D0
    zz = 0.D0
    oo = 1.D0
    c2 = c ! 1.7320508075688767D0
    icpos = 208
    idpos = 208
    ihedge = 104
    allocate(cpos(3,icpos), h_edge(2,ihedge), dpos(idpos) )

    do i=-6,6
      if (abs(i).ne.6) then
        cpos(1:3, 7+i) = (/ -a, dtan(dble(i)*30.D0*dtor*0.5D0) ,  c /)
      else
        cpos(1:3, 7+i) = (/ -a, a ,  c /)
        if (i.lt.0) cpos(2,7+i) = -cpos(2,7+i)
      end if
    end do
    do i=1,13
      cpos(1:3,13+i) = cpos(1:3,i)
      cpos(1,13+i) = -cpos(1,13+i)
    end do

    do i=1,26
      cpos(1:3,26+i) = cpos(1:3,i)
      tmp = cpos(1,26+i)
      cpos(1,26+i) = cpos(2,26+i)
      cpos(2,26+i) = tmp
    end do

    do i=1,52
      cpos(1:3,52+i) = cpos(1:3,i)
      cpos(3,52+i) = -cpos(3,52+i)
    end do

    do i=1,104
      cpos(1:3,104+i) = cpos(1:3,i)
      if (cpos(3,104+i).lt.0.D0) then
        cpos(3,104+i) = -c2
      else
        cpos(3,104+i) = c2
      end if
    end do

    ! and normalize
    do i=1,104
      dpos(i) = dsqrt(sum(cpos(1:3,i)*cpos(1:3,i)))
    end do

    ns = 200
    dx = dt/float(ns-1)

    ! define the connectivity of all the edges
    do i=1,13
      h_edge(1:2,i)    = (/    i, 13+i /)
      h_edge(1:2,13+i) = (/ 26+i, 39+i /)
      h_edge(1:2,26+i) = (/ 52+i, 65+i /)
      h_edge(1:2,39+i) = (/ 78+i, 91+i /)

      h_edge(1:2,52+i) = (/104+i,117+i /)
      h_edge(1:2,65+i) = (/130+i,143+i /)
      h_edge(1:2,78+i) = (/156+i,169+i /)
      h_edge(1:2,91+i) = (/182+i,195+i /)
    end do
  case(4) ! define the coordinates of the tetragonal C4 (4) FZ in Rodrigues Space

    ! NEEDS TO BE FIXED !!!
    a = 57.289922125538D0
    b = 1.0D0
    c = 0.414213562D0
    dt = 114.57984425107713D0
    ds = 2.0D0
    d = 1.7320508075688772D0
    zz = 0.D0
    oo = 1.D0
    icpos = 104
    idpos = 104
    ihedge = 52
    allocate(cpos(3,icpos), h_edge(2,ihedge), dpos(idpos) )

    do i=-6,6
      if (abs(i).ne.6) then
        cpos(1:3, 7+i) = (/ -a, dtan(dble(i)*30.D0*dtor*0.5D0) ,  c /)
      else
        cpos(1:3, 7+i) = (/ -a, a ,  c /)
        if (i.lt.0) cpos(2,7+i) = -cpos(2,7+i)
      end if
    end do
    do i=1,13
      cpos(1:3,13+i) = cpos(1:3,i)
      cpos(1,13+i) = -cpos(1,13+i)
    end do

    do i=1,26
      cpos(1:3,26+i) = cpos(1:3,i)
      tmp = cpos(1,26+i)
      cpos(1,26+i) = cpos(2,26+i)
      cpos(2,26+i) = tmp
    end do

    do i=1,52
      cpos(1:3,52+i) = cpos(1:3,i)
      cpos(3,52+i) = -cpos(3,52+i)
    end do

    ! and normalize
    do i=1,104
      dpos(i) = dsqrt(sum(cpos(1:3,i)*cpos(1:3,i)))
    end do

    ns = 200
    dx = dt/float(ns-1)

    ! define the connectivity of all the edges
    do i=1,13
      h_edge(1:2,i)    = (/    i, 13+i /)
      h_edge(1:2,13+i) = (/ 26+i, 39+i /)
      h_edge(1:2,26+i) = (/ 52+i, 65+i /)
      h_edge(1:2,39+i) = (/ 78+i, 91+i /)
    end do
  case(6) ! define the coordinates of the hexagonal C6 (6) FZ in Rodrigues Space
    a = 57.289922125538D0
    b = 1.0D0
    c = 0.2679491924D0
    dt = 114.57984425107713D0
    ds = 2.0D0
    d = 1.7320508075688772D0
    zz = 0.D0
    oo = 1.D0
    icpos = 104
    idpos = 104
    ihedge = 52
    allocate(cpos(3,icpos), h_edge(2,ihedge), dpos(idpos) )

    do i=-6,6
      if (abs(i).ne.6) then
        cpos(1:3, 7+i) = (/ -a, dtan(dble(i)*30.D0*dtor*0.5D0) ,  c /)
      else
        cpos(1:3, 7+i) = (/ -a, a ,  c /)
        if (i.lt.0) cpos(2,7+i) = -cpos(2,7+i)
      end if
    end do
    do i=1,13
      cpos(1:3,13+i) = cpos(1:3,i)
      cpos(1,13+i) = -cpos(1,13+i)
    end do

    do i=1,26
      cpos(1:3,26+i) = cpos(1:3,i)
      tmp = cpos(1,26+i)
      cpos(1,26+i) = cpos(2,26+i)
      cpos(2,26+i) = tmp
    end do

    do i=1,52
      cpos(1:3,52+i) = cpos(1:3,i)
      cpos(3,52+i) = -cpos(3,52+i)
    end do

    ! and normalize
    do i=1,104
      dpos(i) = dsqrt(sum(cpos(1:3,i)*cpos(1:3,i)))
    end do

    ns = 200
    dx = dt/float(ns-1)

    ! define the connectivity of all the edges
    do i=1,13
      h_edge(1:2,i)    = (/    i, 13+i /)
      h_edge(1:2,13+i) = (/ 26+i, 39+i /)
      h_edge(1:2,26+i) = (/ 52+i, 65+i /)
      h_edge(1:2,39+i) = (/ 78+i, 91+i /)
    end do
  case default
end select

allocate( ropos( 3, (ns+2)*ihedge) )
ropos = 0.D0

icnt = 1
dx = 1.D0/dble(ns)
do i=1, ihedge
  ro1 = r_T( rdinp = (/ cpos(1:3,h_edge(1,i))/dpos(i), dpos(i) /) )
  ro2 = r_T( rdinp = (/ cpos(1:3,h_edge(2,i))/dpos(i), dpos(i) /) )
  rolast = ro1
  ropos(1:3,icnt) = (/ 0.D0, 0.D0, 0.D0 /) 
  icnt = icnt+1
  do j=1,ns+1
    ! if (j.eq.1) then ! make sure to get the starting point
    !   aux4a = rolast%r_copyd()
    !   ropos(1:3,icnt) = aux4a(1:3)*aux4a(4)
    !   icnt = icnt+1
    ! end if 
    aux = dpos(i)*ro1%r_copyd() + dpos(i)*(ro2%r_copyd() - ro1%r_copyd()) * j * dx
    xx = dsqrt( sum (aux(1:3)**2) )
    ro = r_T( rdinp = (/ aux(1:3)/xx, xx /) )
    aux4b = ro%r_copyd()
    ropos(1:3,icnt) = aux4b(1:3)*aux4b(4)
    rolast = ro
    icnt = icnt+1
  end do
end do

end subroutine initFZCyclic_



end module mod_RFZwf