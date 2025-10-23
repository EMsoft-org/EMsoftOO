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
  !!
  !! updated on 10/21/25 to cover the case of rotated RFZs (e.g. 321 vs. 312)

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMRFZwf program
type, public :: RFZwfNameListType
  character(fnlen)  :: hdfname
  character(fnlen)  :: prefix
  character(fnlen)  :: PVexec             ! path to PoVray executable
  character(fnlen)  :: PVincludepath      ! path to PoVray include files
end type RFZwfNameListType

! class definition
type, public :: RFZwf_T
private 
  character(fnlen)                :: nmldeffile = 'EMRFZwf.nml'
  type(RFZwfNameListType)         :: nml 
  integer(kind=irg),dimension(12) :: nLaue = (/ 3, 16, 9, 21, 6, 18, 12, 24, 28, 30, 37, 40 /)
  character(4),dimension(12)      :: Laue = (/ '2   ', '3   ', '4   ', '6   ', '222 ', &
                                               '32  ', '422 ', '622 ', '23  ', '432 ', &
                                               '32R ', '222R' /)
  ! group names in the HDF output file
  character(13)                   :: ortype(5) = (/ 'Cubochoric   ','Homochoric   ','Stereographic', &
                                                    'Rodrigues    ','Euler        '/)

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: RFZwf_
  procedure, pass(self) :: verify_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: RFZwf => RFZwf_
  generic, public :: verify => verify_

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
character(fnlen)  :: prefix
character(fnlen)  :: PVexec             ! path to PoVray executable
character(fnlen)  :: PVincludepath      ! path to PoVray include files

namelist / EMRFZwf / hdfname, prefix, PVexec, PVincludepath

hdfname = 'undefined'
prefix  = 'undefined'
PVexec  = 'undefined'
PVincludepath = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
  open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
  read(UNIT=dataunit,NML=EMRFZwf)
  close(UNIT=dataunit,STATUS='keep')

  if (trim(hdfname).eq.'undefined') then
      call Message%printError('readNameList:',' hdfname is undefined in '//nmlfile)
  end if
  if (PVexec.ne.'undefined') then
    if (trim(prefix).eq.'undefined') then
        call Message%printError('readNameList:',' prefix is undefined in '//nmlfile)
    end if
    if (trim(PVincludepath).eq.'undefined') then
        call Message%printError('readNameList:',' PVincludepath is undefined in '//nmlfile)
    end if
    self%nml%PVexec = PVexec
    self%nml%PVincludepath = PVincludepath
  end if 
end if 

self%nml%hdfname = hdfname 
self%nml%prefix = prefix 
self%nml%PVexec = PVexec
self%nml%PVincludepath = PVincludepath

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
use mod_so3
use HDF5
use mod_HDFsupport

IMPLICIT NONE 

class(RFZwf_T), INTENT(INOUT)       :: self
type(EMsoft_T), INTENT(INOUT)       :: EMsoft
character(fnlen), INTENT(INOUT)     :: progname 

type(HDF_T)                         :: HDF
type(so3_T)                         :: SO
type(IO_T)                          :: Message

real(kind=dbl),allocatable          :: ropos(:,:), sppos(:,:), cupos(:,:), hopos(:,:), eupos(:,:)

integer(kind=irg)                   :: FZorder, FZtype, sz(2), i, hdferr, iL, FZcyclic(4), FZdihedral(4)
character(fnlen)                    :: datafile, dataset, groupname


FZcyclic = (/ 2, 3, 4, 6 /)
FZdihedral = (/ 2, 3, 4, 6 /)

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

! loop over the cyclic Laue groups and do all but the Euler representations
do i=1,12
  if (allocated(ropos)) deallocate(ropos)
  call Message%printMessage(' --> starting on point group '//trim(self%Laue(i)))

! get the Rodrigues wireframe and store it in the Rodrigues group
! pay attention to the rotated RFZs !!! (i>10)
  if (i.le.4) then
    call initFZCyclic_(FZcyclic(i), ropos)
  else
    SO = so3_T( self%nLaue(i), zerolist='FZ')
    call SO%getFZtypeandorder(FZtype, FZorder) 
    write (*,*) self%nLaue(i), FZtype, FZorder 
    if (FZorder.lt.0) then 
      call initFZother_(abs(FZorder), FZtype, ropos, rotorder=abs(FZorder) ) 
    else
      call initFZother_(FZorder, FZtype, ropos)
    end if
  end if
  sz = shape(ropos)
  groupname = 'Rodrigues'
  hdferr = HDF%openGroup(groupname)
  dataset = 'Laue_'//trim(self%Laue(i))
  hdferr = HDF%writeDatasetDoubleArray(dataset, ropos, sz(1), sz(2))
  if (hdferr.ne.0) call HDF%error_check('writeDatasetDoubleArray ropos', hdferr)
  call HDF%pop()

! then do the other representations (except for Euler which is handled separately)
  call convert_orep_(ropos, sz, sppos, cupos, hopos)

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
end do 


! Euler: this is a bit different from the others since there are extra bits
!        to be drawn...in addition, the cyclic RFZs are basically just prisms
! we have to be careful here because the extra lines that were drawn for the 
! dihedral groups, to show the curvature of top and bottom faces more clearly,
! causes issues with the Euler plots; so for the Euler plots, we recalculate
! the ropos array
groupname = 'Euler'
hdferr = HDF%openGroup(groupname)
call Message%printMessage(' Starting on Euler RFZs')
do i=1,12
  if (allocated(ropos)) deallocate(ropos)
  if (allocated(eupos)) deallocate(eupos)
  call Message%printMessage(' --> starting on point group '//trim(self%Laue(i)))

! get the Euler wireframe 
  SO = so3_T( self%nLaue(i), zerolist='FZ')
  call SO%getFZtypeandorder(FZtype, FZorder) 
! for the non-cyclic groups, we need to first get the ropos array
  if (i.gt.4) then 
    if (FZorder.lt.0) then 
      call initFZother_(abs(FZorder), FZtype, ropos, rotorder=abs(FZorder), euler=.TRUE.)
    else
      call initFZother_(FZorder, FZtype, ropos, euler=.TRUE.)
    end if
    sz = shape(ropos)
    call EulerinitFZ_(FZorder, FZtype, eupos, ropos)
  else
    call EulerinitFZ_(FZorder, FZtype, eupos)
  end if
  sz = shape(eupos)

  dataset = 'Laue_'//trim(self%Laue(i))
  hdferr = HDF%writeDatasetDoubleArray(dataset, eupos, sz(1), sz(2))
  if (hdferr.ne.0) call HDF%error_check('writeDatasetDoubleArray eupos', hdferr)
end do 

! and close the HDF5 file
call HDF%popall()

! next we generate all the RFZdrawings for all the point groups to make sure they
! are correct.  This requires a routine that employs the PoVray module to generate
! the correct drawing volumes, followed by reading the wire frame from the hdf file
! that we just generated.
if (self%nml%prefix.ne.'undefined') call self%verify_(EMsoft, HDF)

call closeFortranHDFInterface()

end subroutine RFZwf_

!--------------------------------------------------------------------------
recursive subroutine convert_orep_(ropos, sz, sppos, cupos, hopos)
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

type(r_T)                             :: r
type(s_T)                             :: s
type(c_T)                             :: c
type(h_T)                             :: h

integer(kind=irg)                     :: i 
real(kind=dbl)                        :: x

call setRotationPrecision('d')

allocate( sppos(sz(1),sz(2)), cupos(sz(1),sz(2)), hopos(sz(1),sz(2)) )

sppos = 0.D0 
cupos = 0.D0 
hopos = 0.D0 

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
                                         tPi, hpi, aux(4), aux4a(4), aux4b(4)

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

!--------------------------------------------------------------------------
recursive subroutine EulerinitFZ_(FZorder, FZtype, eupos, ropos)
!DEC$ ATTRIBUTES DLLEXPORT :: EulerinitFZ_
!! author: MDG
!! version: 1.0
!! date: 09/09/25
!!
!! generate the coordinates of the wireframe for the Euler representation.

use mod_rotations
use mod_io

IMPLICIT NONE

integer(kind=irg),INTENT(IN)          :: FZorder
integer(kind=irg),INTENT(IN)          :: FZtype
real(kind=dbl),INTENT(OUT),allocatable:: eupos(:,:)
real(kind=dbl),INTENT(IN),OPTIONAL    :: ropos(:,:)

type(e_T)                             :: eu, eulast
type(r_T)                             :: ro
type(IO_T)                            :: Message

real(kind=dbl)                        :: rmax, dx, r, xmax, x, y, z, zsmall, ac, sh(3), xx, &
                                         tPi, aux(4), aux4a(4), aux4b(4), hPi

integer(kind=irg)                     :: i,j,k, icnt, imax, nt, ns, idpos, icpos, ihedge, rosz(2)
integer(kind=irg),allocatable         :: h_edge(:,:)
real(kind=dbl),allocatable            :: cpos(:,:), dpos(:)
! parameters that depend on the cyclic group
real(kind=dbl)                        :: a, b, c, dt, ds, d, dd, zz, oo, c2, tmp

tPi = 2.D0 * cPi
hPi = 0.5D0 * cPi
sh = (/ cPi, cPi/2.D0, cPi /)

if (present(ropos)) then 
  rosz = shape(ropos)
end if

if (FZtype.eq.1) then   ! these are the cyclic groups
  allocate( eupos(3,48) )
  eupos = 0.D0
  xx = cPi/dble(FZorder)
  ! draw four diagonal lines
  icnt = 1
  eupos(1:3, icnt+1) = (/ 0.D0, 0.D0, xx /) - sh
  eupos(1:3, icnt+2) = (/ xx, 0.D0, 0.D0 /) - sh
  icnt = icnt + 3 
  eupos(1:3, icnt+1) = (/ 0.D0, 0.D0, tPi-xx /) - sh
  eupos(1:3, icnt+2) = (/ tPi-xx, 0.D0, 0.D0 /) - sh
  icnt = icnt + 3 
  eupos(1:3, icnt+1) = (/ xx, 0.D0, tPi /) - sh
  eupos(1:3, icnt+2) = (/ tPi, 0.D0, xx /) - sh
  icnt = icnt + 3 
  eupos(1:3, icnt+1) = (/ tPi-xx, 0.D0, tPi /) - sh
  eupos(1:3, icnt+2) = (/ tPi, 0.D0, tPi-xx /) - sh
  icnt = icnt + 3 
! for cyclic groups we also need to draw the diagonals in the top surface
! and the vertical lines connecting bottom and top planes
! top plane
  eupos(1:3, icnt+1) = (/ 0.D0, cPi, xx /) - sh
  eupos(1:3, icnt+2) = (/ xx, cPi, 0.D0 /) - sh
  icnt = icnt+3
  eupos(1:3, icnt+1) = (/ 0.D0, cPi, tPi-xx /) - sh
  eupos(1:3, icnt+2) = (/ tPi-xx, cPi, 0.D0 /) - sh
  icnt = icnt+3
  eupos(1:3, icnt+1) = (/ xx, cPi, tPi /) - sh
  eupos(1:3, icnt+2) = (/ tPi, cPi, xx /) - sh
  icnt = icnt+3
  eupos(1:3, icnt+1) = (/ tPi-xx, cPi, tPi /) - sh
  eupos(1:3, icnt+2) = (/ tPi, cPi, tPi-xx /) - sh
  icnt = icnt+3
  ! verticals
  eupos(1:3, icnt+1) = (/ xx, cPi, 0.D0 /) - sh
  eupos(1:3, icnt+2) = (/ xx, 0.D0, 0.D0 /) - sh
  icnt = icnt+3
  eupos(1:3, icnt+1) = (/ 0.D0, cPi, xx /) - sh 
  eupos(1:3, icnt+2) = (/ 0.D0, 0.D0, xx /) - sh
  icnt = icnt+3
  eupos(1:3, icnt+1) = (/ tPi-xx, cPi, 0.D0 /) - sh
  eupos(1:3, icnt+2) = (/ tPi-xx, 0.D0, 0.D0 /) - sh
  icnt = icnt+3
  eupos(1:3, icnt+1) = (/ 0.D0, cPi, tPi-xx /) - sh
  eupos(1:3, icnt+2) = (/ 0.D0, 0.D0, tPi-xx /) - sh
  icnt = icnt+3
  eupos(1:3, icnt+1) = (/ tPi, cPi, xx /) - sh
  eupos(1:3, icnt+2) = (/ tPi, 0.D0, xx /) - sh
  icnt = icnt+3
  eupos(1:3, icnt+1) = (/ xx, cPi, tPi /) - sh
  eupos(1:3, icnt+2) = (/ xx, 0.D0, tPi /) - sh
  icnt = icnt+3
  eupos(1:3, icnt+1) = (/ tPi, cPi, tPi-xx /) - sh
  eupos(1:3, icnt+2) = (/ tPi, 0.D0, tPi-xx /) - sh
  icnt = icnt+3
  eupos(1:3, icnt+1) = (/ tPi-xx, cPi, tPi /) - sh
  eupos(1:3, icnt+2) = (/ tPi-xx, 0.D0, tPi /) - sh

! and return to the calling routine
  RETURN
end if 

if (FZtype.eq.2) then   ! these are the dihedral groups
! first we need to convert the existing ropos entries to Euler space
  allocate( eupos(3,rosz(2) + 36) )
  eupos = 0.D0
  do i=1,rosz(2)
    xx = sum(abs(ropos(1:3,i)))
    if (xx.ne.0.D0) then
      xx = sqrt(sum(ropos(1:3,i)**2))
      ro = r_T( rdinp = (/ ropos(1:3,i)/xx, xx /) )
      eu = ro%re()
      eupos(1:3,i) = eu%e_copyd() 
      eupos(1,i) = mod(eupos(1,i)+10.D0*cPi,2.D0*cPi)
      eupos(2,i) = mod(eupos(2,i)+10.D0*cPi,cPi)
      eupos(3,i) = mod(eupos(3,i)+10.D0*cPi,2.D0*cPi)
      eupos(1:3,i) = eupos(1:3,i) - sh
    end if
  end do   

  ! draw four diagonal lines
  xx = cPi/dble(FZorder)
  icnt = rosz(2) + 1
  eupos(1:3, icnt+1) = (/ 0.D0, 0.D0, xx /) - sh
  eupos(1:3, icnt+2) = (/ xx, 0.D0, 0.D0 /) - sh
  icnt = icnt + 3 
  eupos(1:3, icnt+1) = (/ 0.D0, 0.D0, tPi-xx /) - sh
  eupos(1:3, icnt+2) = (/ tPi-xx, 0.D0, 0.D0 /) - sh
  icnt = icnt + 3 
  eupos(1:3, icnt+1) = (/ xx, 0.D0, tPi /) - sh
  eupos(1:3, icnt+2) = (/ tPi, 0.D0, xx /) - sh
  icnt = icnt + 3 
  eupos(1:3, icnt+1) = (/ tPi-xx, 0.D0, tPi /) - sh
  eupos(1:3, icnt+2) = (/ tPi, 0.D0, tPi-xx /) - sh
  icnt = icnt + 3 

! the verticals need to be drawn but only up to the level of the FZ surface
  eupos(1:3, icnt+1) = (/ xx, hPi, 0.D0 /) - sh
  eupos(1:3, icnt+2) = (/ xx, 0.D0, 0.D0 /) - sh
  icnt = icnt+3
  eupos(1:3, icnt+1) = (/ 0.D0, hPi, xx /) - sh
  eupos(1:3, icnt+2) = (/ 0.D0, 0.D0, xx /) - sh
  icnt = icnt+3
  eupos(1:3, icnt+1) = (/ tPi-xx, hPi, 0.D0 /) - sh
  eupos(1:3, icnt+2) = (/ tPi-xx, 0.D0, 0.D0 /) - sh
  icnt = icnt+3
  eupos(1:3, icnt+1) = (/ 0.D0, hPi, tPi-xx /) - sh
  eupos(1:3, icnt+2) = (/ 0.D0, 0.D0, tPi-xx /) - sh
  icnt = icnt+3
  eupos(1:3, icnt+1) = (/ tPi, hPi, xx /) - sh
  eupos(1:3, icnt+2) = (/ tPi, 0.D0, xx /) - sh
  icnt = icnt+3
  eupos(1:3, icnt+1) = (/ xx, hPi, tPi /) - sh
  eupos(1:3, icnt+2) = (/ xx, 0.D0, tPi /) - sh
  icnt = icnt+3
  eupos(1:3, icnt+1) = (/ tPi, hPi, tPi-xx /) - sh
  eupos(1:3, icnt+2) = (/ tPi, 0.D0, tPi-xx /) - sh
  icnt = icnt+3
  eupos(1:3, icnt+1) = (/ tPi-xx, hPi, tPi /) - sh
  eupos(1:3, icnt+2) = (/ tPi-xx, 0.D0, tPi /) - sh

! and return to the calling routine
  RETURN
end if

if (FZtype.eq.3) then   ! this is the tetrahedral group
! first we need to convert the existing ropos entries to Euler space
  allocate( eupos(3,rosz(2) + 48) )
  eupos = 0.D0
  do i=1,rosz(2)
    xx = sum(abs(ropos(1:3,i)))
    if (xx.ne.0.D0) then
    xx = sqrt(sum(ropos(1:3,i)**2))
      ro = r_T( rdinp = (/ ropos(1:3,i)/xx, xx /) )
      eu = ro%re()
      eupos(1:3,i) = eu%e_copyd() 
      eupos(1,i) = mod(eupos(1,i)+10.D0*cPi,2.D0*cPi)
      eupos(2,i) = mod(eupos(2,i)+10.D0*cPi,cPi)
      eupos(3,i) = mod(eupos(3,i)+10.D0*cPi,2.D0*cPi)
      eupos(1:3,i) = eupos(1:3,i) - sh
    end if 
  end do   

  ! draw four diagonal lines
  icnt = rosz(2) + 1
  xx = cPi/dble(2)
  ! draw four diagonal lines
  eupos(1:3, icnt+1) = (/ 0.D0, 0.D0, xx /) - sh
  eupos(1:3, icnt+2) = (/ xx, 0.D0, 0.D0 /) - sh
  icnt = icnt + 3
  eupos(1:3, icnt+1) = (/ 0.D0, 0.D0, tpi-xx /) - sh
  eupos(1:3, icnt+2) = (/ tpi-xx, 0.D0, 0.D0 /) - sh
  icnt = icnt + 3
  eupos(1:3, icnt+1) = (/ xx, 0.D0, tpi /) - sh
  eupos(1:3, icnt+2) = (/ tpi, 0.D0, xx /) - sh
  icnt = icnt + 3
  eupos(1:3, icnt+1) = (/ tpi-xx, 0.D0, tpi /) - sh
  eupos(1:3, icnt+2) = (/ tpi, 0.D0, tpi-xx /) - sh
  icnt = icnt + 3

! and finally the corner posts
  eupos(1:3, icnt+1) = (/ hPi, 0.D0, 0.D0 /) - sh 
  eupos(1:3, icnt+2) = (/ 0.D0, 0.D0, 0.D0 /) - sh
  icnt = icnt + 3
  eupos(1:3, icnt+1) = (/ 0.D0, hPi, 0.D0 /) - sh
  eupos(1:3, icnt+2) = (/ 0.D0, 0.D0, 0.D0 /) - sh
  icnt = icnt + 3
  eupos(1:3, icnt+1) = (/ 0.D0, 0.D0, hPi /) - sh
  eupos(1:3, icnt+2) = (/ 0.D0, 0.D0, 0.D0 /) - sh
  icnt = icnt + 3

  eupos(1:3, icnt+1) = (/ tpi - hPi, 0.D0, 0.D0 /) - sh
  eupos(1:3, icnt+2) = (/ tpi, 0.D0, 0.D0 /) - sh
  icnt = icnt + 3
  eupos(1:3, icnt+1) = (/ tpi, hPi, 0.D0 /) - sh
  eupos(1:3, icnt+2) = (/ tpi, 0.D0, 0.D0 /) - sh
  icnt = icnt + 3
  eupos(1:3, icnt+1) = (/ tpi, 0.D0, hPi /) - sh 
  eupos(1:3, icnt+2) = (/ tpi, 0.D0, 0.D0 /) - sh
  icnt = icnt + 3

  eupos(1:3, icnt+1) = (/ tpi - hPi, 0.D0, tpi /) - sh
  eupos(1:3, icnt+2) = (/ tpi, 0.D0, tpi /) - sh
  icnt = icnt + 3
  eupos(1:3, icnt+1) = (/ tpi, hPi, tpi /) - sh
  eupos(1:3, icnt+2) = (/ tpi, 0.D0, tpi /) - sh
  icnt = icnt + 3
  eupos(1:3, icnt+1) = (/ tpi, 0.D0, tpi-hPi /) - sh
  eupos(1:3, icnt+2) = (/ tpi, 0.D0, tpi /) - sh
  icnt = icnt + 3

  eupos(1:3, icnt+1) = (/ 0.D0 + hPi, 0.D0, tpi /) - sh
  eupos(1:3, icnt+2) = (/ 0.D0, 0.D0, tpi /) - sh
  icnt = icnt + 3
  eupos(1:3, icnt+1) = (/ 0.D0, hPi, tpi /) - sh 
  eupos(1:3, icnt+2) = (/ 0.D0, 0.D0, tpi /) - sh 
  icnt = icnt + 3
  eupos(1:3, icnt+1) = (/ 0.D0, 0.D0, tpi-hPi /) - sh
  eupos(1:3, icnt+2) = (/ 0.D0, 0.D0, tpi /) - sh 

! and return to the calling routine
  RETURN
end if

if (FZtype.eq.4) then   ! this is the octahedral group
! first we need to convert the existing ropos entries to Euler space
  allocate( eupos(3,rosz(2) + 96) )
  eupos = 0.D0
  do i=1,rosz(2)
    xx = sqrt(sum(ropos(1:3,i)**2))
    if (xx.ne.0.D0) then
      ro = r_T( rdinp = (/ ropos(1:3,i)/xx, xx /) )
      eu = ro%re()
      eupos(1:3,i) = eu%e_copyd() 
      eupos(1,i) = mod(eupos(1,i)+10.D0*cPi,2.D0*cPi)
      eupos(2,i) = mod(eupos(2,i)+10.D0*cPi,cPi)
      eupos(3,i) = mod(eupos(3,i)+10.D0*cPi,2.D0*cPi)
      eupos(1:3,i) = eupos(1:3,i) - sh
    end if 
  end do   

  ! draw four diagonal lines
  icnt = rosz(2) + 1

  xx = cPi/dble(4)
! draw four diagonal lines
  eu = e_T( edinp = (/ xx, 0.D0, 0.D0 /) - sh )
  eulast = e_T( edinp = (/ 0.D0, 0.D0, xx /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ tpi-xx, 0.D0, 0.D0 /) - sh )
  eulast = e_T( edinp = (/ 0.D0, 0.D0, tpi-xx /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ tpi, 0.D0, xx /) - sh )
  eulast = e_T( edinp = (/ xx, 0.D0, tpi /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ tpi, 0.D0, tpi-xx /) - sh )
  eulast = e_T( edinp = (/ tpi-xx, 0.D0, tpi /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
! and verticals
  hPi = hPi * 0.5D0
  eu = e_T( edinp = (/ xx, 0.D0, 0.D0 /) - sh )
  eulast = e_T( edinp = (/ xx, hPi, 0.D0 /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ 0.D0, 0.D0, xx /) - sh )
  eulast = e_T( edinp = (/ 0.D0, hPi, xx /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ tpi-xx, 0.D0, 0.D0 /) - sh )
  eulast = e_T( edinp = (/ tpi-xx, hPi, 0.D0 /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ 0.D0, 0.D0, tpi-xx /) - sh )
  eulast = e_T( edinp = (/ 0.D0, hPi, tpi-xx /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ tpi, 0.D0, xx /) - sh )
  eulast = e_T( edinp = (/ tpi, hPi, xx /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ xx, 0.D0, tpi /) - sh )
  eulast = e_T( edinp = (/ xx, hPi, tpi /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ tpi, 0.D0, tpi-xx /) - sh )
  eulast = e_T( edinp = (/ tpi, hPi, tpi-xx /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ tpi-xx, 0.D0, tpi /) - sh )
  eulast = e_T( edinp = (/ tpi-xx, hPi, tpi /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
! and the closing segments
  eu = e_T( edinp = (/ xx, hPi, 0.D0 /) - sh )
  eulast = e_T( edinp = (/ 0.D0, hPi, 0.D0 /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ 0.D0, hPi, xx /) - sh )
  eulast = e_T( edinp = (/ 0.D0, hPi, 0.D0 /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ tpi-xx, hPi, 0.D0 /) - sh )
  eulast = e_T( edinp = (/ tpi, hPi, 0.D0 /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ 0.D0, hPi, tpi-xx /) - sh )
  eulast = e_T( edinp = (/ 0.D0, hPi, tpi /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ tpi, hPi, xx /) - sh )
  eulast = e_T( edinp = (/ tpi, hPi, 0.D0 /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ xx, hPi, tpi /) - sh )
  eulast = e_T( edinp = (/ 0.D0, hPi, tpi /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ tpi, hPi, tpi-xx /) - sh )
  eulast = e_T( edinp = (/ tpi, hPi, tpi /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ tpi-xx, hPi, tpi /) - sh )
  eulast = e_T( edinp = (/ tpi, hPi, tpi /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
! and finally the corner posts
  eu = e_T( edinp = (/ 0.D0, 0.D0, 0.D0 /) - sh )
  eulast = e_T( edinp = (/ hPi, 0.D0, 0.D0 /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ 0.D0, 0.D0, 0.D0 /) - sh )
  eulast = e_T( edinp = (/ 0.D0, hPi, 0.D0 /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ 0.D0, 0.D0, 0.D0 /) - sh )
  eulast = e_T( edinp = (/ 0.D0, 0.D0, hPi /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3

  eu = e_T( edinp = (/ tpi, 0.D0, 0.D0 /) - sh )
  eulast = e_T( edinp = (/ tpi - hPi, 0.D0, 0.D0 /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ tpi, 0.D0, 0.D0 /) - sh )
  eulast = e_T( edinp = (/ tpi, hPi, 0.D0 /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ tpi, 0.D0, 0.D0 /) - sh )
  eulast = e_T( edinp = (/ tpi, 0.D0, hPi /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3

  eu = e_T( edinp = (/ tpi, 0.D0, tpi /) - sh )
  eulast = e_T( edinp = (/ tpi - hPi, 0.D0, tpi /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ tpi, 0.D0, tpi /) - sh )
  eulast = e_T( edinp = (/ tpi, hPi, tpi /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ tpi, 0.D0, tpi /) - sh )
  eulast = e_T( edinp = (/ tpi, 0.D0, tpi-hPi /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3

  eu = e_T( edinp = (/ 0.D0, 0.D0, tpi /) - sh )
  eulast = e_T( edinp = (/ 0.D0 + hPi, 0.D0, tpi /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ 0.D0, 0.D0, tpi /) - sh )
  eulast = e_T( edinp = (/ 0.D0, hPi, tpi /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()
  icnt = icnt+3
  eu = e_T( edinp = (/ 0.D0, 0.D0, tpi /) - sh )
  eulast = e_T( edinp = (/ 0.D0, 0.D0, tpi-hPi /) - sh )
  eupos(1:3, icnt+1) = eulast%e_copyd()
  eupos(1:3, icnt+2) = eu%e_copyd()

! and return to the calling routine
  RETURN
end if

end subroutine EulerinitFZ_

!--------------------------------------------------------------------------
recursive subroutine initFZother_(FZorder, FZtype, ropos, euler, rotorder)
!DEC$ ATTRIBUTES DLLEXPORT :: initFZother_
!! author: MDG
!! version: 1.0
!! date: 09/08/25
!!
!! generate the coordinates of the wireframe for the cyclic rotational groups
!! these are stored as 3-component Rodrigues vectors, with (0,0,0) entries
!! separating the major line segments.
!!
!! 10/21/25: added option to rotate the coordinates for second settings of 
!! some of the point groups (e.g., 32 vs. 312)

use mod_rotations
use mod_povray

IMPLICIT NONE

integer(kind=irg),INTENT(IN)          :: FZorder
integer(kind=irg),INTENT(IN)          :: FZtype
real(kind=dbl),INTENT(OUT),allocatable:: ropos(:,:)
logical,INTENT(IN),OPTIONAL           :: euler
integer(kind=irg),INTENT(IN),OPTIONAL :: rotorder

type(r_T)                             :: ro1, ro2, rolast, ro
type(PoVRay_T)                        :: PoV

real(kind=dbl)                        :: aux4b(4), d, aux(3), xx, dx

integer(kind=irg)                     :: i,j,k, icnt, dims(3), nt, ns, rotate
integer(kind=irg),allocatable         :: s_edge(:,:), t_edge(:,:)
real(kind=dbl),allocatable            :: cpos(:,:)
logical                               :: twostep

rotate = 0
if (present(rotorder)) rotate = rotorder

! use the PoVray routines to get all the coordinates and connectivities
if (FZtype.eq.2) then
    if (FZorder.eq.6) then
        twostep = .TRUE.
        if (present(euler)) then
          dims = (/ 24, 24, 12 /)
          allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
          call PoV%getpos_FZ622(dims, cpos, s_edge, t_edge, ns, d, nt, rotate, euler=.TRUE.)
        else
          dims = (/ 24, 24, 24 /)
          allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
          call PoV%getpos_FZ622(dims, cpos, s_edge, t_edge, ns, d, nt, rotate)
        end if
    end if
    if (FZorder.eq.4) then
        twostep = .TRUE.
        if (present(euler)) then
          dims = (/ 16, 16, 8 /)
          allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
          call PoV%getpos_FZ422(dims, cpos, s_edge, t_edge, ns, d, nt, rotate, euler=.TRUE.)
        else
          dims = (/ 16, 16, 16 /)
          allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
          call PoV%getpos_FZ422(dims, cpos, s_edge, t_edge, ns, d, nt, rotate)
        end if
    end if
    if (FZorder.eq.3) then
        twostep = .TRUE.
        if (present(euler)) then
          dims = (/ 12, 12, 6 /)
          allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
          call PoV%getpos_FZ32(dims, cpos, s_edge, t_edge, ns, d, nt, rotate, euler=.TRUE.)
        else
          dims = (/ 12, 12, 12 /)
          allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
          call PoV%getpos_FZ32(dims, cpos, s_edge, t_edge, ns, d, nt, rotate)
        end if
    end if
    if (FZorder.eq.2) then
        twostep = .TRUE.
        if (present(euler)) then
          dims = (/ 8, 8, 4 /)
          allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
          call PoV%getpos_FZ222(dims, cpos, s_edge, t_edge, ns, d, nt, rotate, euler=.TRUE.)
        else
          dims = (/ 8, 8, 16 /)
          allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
          call PoV%getpos_FZ222(dims, cpos, s_edge, t_edge, ns, d, nt, rotate)
        end if
    end if
end if

if (FZtype.eq.3) then
! rotational group 23
      twostep = .FALSE.
      dims = (/ 6, 12, 1 /)
      allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
      call PoV%getpos_FZ23(dims, cpos, s_edge, t_edge, ns, d, nt)
      nt = 0  ! to avoid extra point appearing in drawings...
end if

if (FZtype.eq.4) then
! rotational group 432
      twostep = .TRUE.
      dims = (/ 24, 12, 24 /)
      allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
      call PoV%getpos_FZ432(dims, cpos, s_edge, t_edge, ns, d, nt)
end if

! allocate( ropos(3, (ns+2)*dims(2) + (nt+2)*dims(3) ) )
! allocate( ropos(3, (ns+1)*dims(2) + (nt+1)*dims(3) ) )
allocate( ropos(3, ns*dims(2) + nt*dims(3) ) )
ropos = 0.D0

! next we determine the actual Rodrigues coordinates that will go into the ropos array 
icnt = 1
 dx = 1.D0/dble(ns)
 do i=1,dims(2)
  ro1 = r_T( rdinp = (/ cpos(1:3,s_edge(1,i)), d /) )
  ro2 = r_T( rdinp = (/ cpos(1:3,s_edge(2,i)), d /) )
  rolast = ro1
  ropos(1:3,icnt) = (/ 0.D0, 0.D0, 0.D0 /) 
  icnt = icnt+1
  do j=1,ns-1!+1
    aux = d*ro1%r_copyd() + d*(ro2%r_copyd() - ro1%r_copyd()) * j * dx
    xx = dsqrt( sum (aux(1:3)**2) )
    ro = r_T( rdinp = (/ aux(1:3)/xx, xx /) )
    aux4b = ro%r_copyd()
    ropos(1:3,icnt) = aux4b(1:3)*aux4b(4)
    rolast = ro
    icnt = icnt+1
  end do
 end do

 if (twostep) then
   dx = 1.D0/dble(nt)
   do i=1,dims(3)
    ro1 = r_T( rdinp = (/ cpos(1:3,t_edge(1,i)), d /) )
    ro2 = r_T( rdinp = (/ cpos(1:3,t_edge(2,i)), d /) )
    rolast = ro1
    ropos(1:3,icnt) = (/ 0.D0, 0.D0, 0.D0 /) 
    icnt = icnt+1
    do j=1,nt-1!+1
      aux = d*ro1%r_copyd() + d*(ro2%r_copyd() - ro1%r_copyd()) * j * dx
      xx = dsqrt( sum (aux(1:3)**2) )
      ro = r_T( rdinp = (/ aux(1:3)/xx, xx /) )
      aux4b = ro%r_copyd()
      ropos(1:3,icnt) = aux4b(1:3)*aux4b(4)
      rolast = ro
      icnt = icnt+1
    end do
   end do
 end if

end subroutine initFZother_

!--------------------------------------------------------------------------
recursive subroutine verify_(self, EMsoft, HDF)
!DEC$ ATTRIBUTES DLLEXPORT :: verify_
!! author: MDG
!! version: 1.0
!! date: 09/09/25

use mod_EMsoft
use HDF5
use mod_HDFsupport
use mod_povray
use ISO_C_BINDING

class(RFZwf_T), INTENT(INOUT)       :: self
type(EMsoft_T), INTENT(INOUT)       :: EMsoft
type(HDF_T),INTENT(INOUT)           :: HDF

type(PoVRay_T)                      :: PoV

real(real_kind_15),allocatable      :: wireframe(:,:)
integer(HSIZE_T)                    :: dims(2)

integer(kind=irg)                   :: sz(2), i, j, iG, iL, hdferr, icnt
character(fnlen)                    :: povname, groupname, dataset
character(fnlen)                    :: locationline, locline2, datafile, pvcmd, str, skyline
real(kind=sgl)                      :: eyepos(3), dd, dis(5)
real(kind=dbl)                      :: cylr, ac, va(5)
character(9)                        :: px, py, pz, pd
character(21)                       :: p1, p2
logical                             :: readonly = .TRUE., fexists

! initialize PoVray parameters
cylr = 0.005D0
dis = (/ 4.0, 4.0, 2.5, 3.5, 11.0 /)
va = (/ 15.D0, 15.D0, 15.D0, 15.D0, 155.D0 /)
locline2 = "location < "
eyepos = (/ 0.911259, 0.0, 0.112 /)
eyepos = eyepos/sqrt( sum( eyepos*eyepos))
write (px,"(F9.3)") eyepos(1)
write (py,"(F9.3)") eyepos(2)
write (pz,"(F9.3)") eyepos(3)

p1 = "*cos(clck*0.0174533)"
p2 = "*sin(clck*0.0174533)"

locline2 = trim(locline2)//px//p1//"-"//py//p2//","//px//p2//"+"//py//p1//","//pz//">*"

datafile = EMsoft%generateFilePath('EMdatapathname', self%nml%hdfname)

hdferr = HDF%openFile(datafile, readonly)

! loop over the orientation representations
! 1(cubochoric)|2(homochoric)|3(stereographic)|4(Rodrigues)|5(Euler)
 do iG=1,5 
  groupname = trim(self%ortype(iG))
  hdferr = HDF%openGroup(groupname)
  do iL=1,12
! read the dataset
    dataset = 'Laue_'//trim(self%Laue(iL))
    call HDF%readDatasetDoubleArray(dataset, dims, hdferr, wireframe)
    sz = dims
! set up the PoVray output file
    datafile = EMsoft%generateFilePath('EMdatapathname', self%nml%prefix)
    datafile = trim(datafile)//'_'//trim(groupname)//'_'//trim(self%Laue(iL))
    povname = trim(datafile)//'.pov'
    dd = dis(iG)
    write (pd,"(F9.3)") dd 
    locationline = trim(locline2)//pd
    if (iG.lt.5) then 
      skyline = 'sky <0.0, 0.0, 1.0>'
    else
      skyline = 'sky <0.0, 1.0, 0.0>'
    end if
    PoV = PoVRay_T( EMsoft, povname, locationline=locationline, skyline=skyline, viewangle = va(iG) )
    write (*,*) ' Creating '//trim(povname), dims, hdferr
! add the reference frame and any necessary wireframes
    if (iG.eq.1) then
      ac = 0.5D0 * LPs%ap
      call PoV%addReferenceFrame(ac, cylr)
      call PoV%addCubochoricCube()
      call PoV%addOrigin( (/ 0.0, 0.0, 0.0 /) )
    end if
    if (iG.eq.2) then
      ac = 1.33067D0
      call PoV%addReferenceFrame(ac, cylr)
      call PoV%addWireFrameSphere(ac)
      call PoV%addOrigin( (/ 0.0, 0.0, 0.0 /) )
    end if
    if (iG.eq.3) then
      ac = 1.0D0
      call PoV%addReferenceFrame(ac, cylr)
      call PoV%addWireFrameSphere(ac)
      call PoV%addOrigin( (/ 0.0, 0.0, 0.0 /) )
    end if
    if (iG.eq.4) then
      ac = 1.0D0
      call PoV%addReferenceFrame(ac, cylr)
      call PoV%addOrigin( (/ 0.0, 0.0, 0.0 /) )
    end if
    if (iG.eq.5) then
      call PoV%addEulerBox()
      call PoV%addOrigin( (/ -3.141593,-1.570796,-3.141593 /) )
    end if
! add the current wireframe array as rendered cylinders
    icnt = 1
    do while (icnt.lt.sz(2)) 
      if (sum(abs(wireframe(1:3,icnt))).eq.0.D0) then 
        icnt = icnt+1
      end if
      if (sum(abs(wireframe(1:3,icnt+1))).ne.0.D0) then 
        call PoV%addCylinder(wireframe(1:3,icnt),wireframe(1:3,icnt+1),cylr,(/ 1.0, 0.0, 0.0 /))
      end if 
      icnt = icnt+1
    end do
    call PoV%closeFile()

! and run the rendering using a povray.ini file
!---------------------------------------------------------------------
! next we generate the PoVRay.ini file with the rendering instructions
    open(unit=dataunit,file='povray.ini',status='unknown',form='formatted')
    write(dataunit,"('Input_File_Name=',A)") trim(povname)
    write(dataunit,"('Output_File_Name=',A)") trim(datafile)//'.png'

    str = '+L'//trim(self%nml%PVincludepath)
    write(dataunit,"(A)") trim(str)
    write(dataunit,"('+W',I4,' +H',I4)") 1024, 1024
    write(dataunit,"('Initial_Clock=1')")
    write(dataunit,"('Initial_Frame=1')")
    write(dataunit,"('Final_Clock=1')")
    write(dataunit,"('Final_Frame=1')")
    write(dataunit,"('Work_Threads=6')")
    close(unit=dataunit,status='keep')

    pvcmd = trim(self%nml%PVexec)
    inquire(file=trim(pvcmd),exist=fexists)
    if (fexists.eqv..TRUE.) then
      call system(trim(pvcmd)//' povray.ini >/dev/null 2>/dev/null')
    end if
  end do 
  call HDF%pop()
end do 

call HDF%popall()

end subroutine verify_

end module mod_RFZwf