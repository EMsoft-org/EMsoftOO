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

module mod_kinematical
  !! author: MDG 
  !! version: 1.0 
  !! date: 03/26/20
  !!
  !! class definition for the EMkinematical program

use mod_kinds
use mod_global

IMPLICIT NONE 

private :: AntiAlias, KinematicalLines, KinematicalBands 

! namelist for the EMkinematical program
type, public :: kinematicalNameListType
  real(kind=sgl)   :: dmin
  real(kind=sgl)   :: thr
  real(kind=sgl)   :: voltage
  integer(kind=irg):: nx
  character(fnlen) :: xtalname
  character(fnlen) :: datafile
  character(5)     :: mode
end type kinematicalNameListType

! class definition
type, public :: kinematical_T
private 
  character(fnlen)       :: nmldeffile = 'EMkinematical.nml'
  type(kinematicalNameListType)  :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: kinematical_
  procedure, pass(self) :: get_dmin_
  procedure, pass(self) :: get_thr_
  procedure, pass(self) :: get_voltage_
  procedure, pass(self) :: get_nx_
  procedure, pass(self) :: get_xtalname_
  procedure, pass(self) :: get_datafile_
  procedure, pass(self) :: get_mode_
  procedure, pass(self) :: set_dmin_
  procedure, pass(self) :: set_thr_
  procedure, pass(self) :: set_voltage_
  procedure, pass(self) :: set_nx_
  procedure, pass(self) :: set_xtalname_
  procedure, pass(self) :: set_datafile_
  procedure, pass(self) :: set_mode_
  procedure, pass(self) :: KinematicalLines
  procedure, pass(self) :: KinematicalBands

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: kinematical => kinematical_
  generic, public :: get_dmin => get_dmin_
  generic, public :: get_thr => get_thr_
  generic, public :: get_voltage => get_voltage_
  generic, public :: get_nx => get_nx_
  generic, public :: get_xtalname => get_xtalname_
  generic, public :: get_datafile => get_datafile_
  generic, public :: get_mode => get_mode_
  generic, public :: set_dmin => set_dmin_
  generic, public :: set_thr => set_thr_
  generic, public :: set_voltage => set_voltage_
  generic, public :: set_nx => set_nx_
  generic, public :: set_xtalname => set_xtalname_
  generic, public :: set_datafile => set_datafile_
  generic, public :: set_mode => set_mode_
end type kinematical_T

!DEC$ ATTRIBUTES DLLEXPORT :: getNameList
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList
!DEC$ ATTRIBUTES DLLEXPORT :: kinematical
!DEC$ ATTRIBUTES DLLEXPORT :: get_dmin
!DEC$ ATTRIBUTES DLLEXPORT :: set_dmin
!DEC$ ATTRIBUTES DLLEXPORT :: get_thr
!DEC$ ATTRIBUTES DLLEXPORT :: set_thr
!DEC$ ATTRIBUTES DLLEXPORT :: get_voltage
!DEC$ ATTRIBUTES DLLEXPORT :: set_voltage
!DEC$ ATTRIBUTES DLLEXPORT :: get_nx
!DEC$ ATTRIBUTES DLLEXPORT :: set_nx
!DEC$ ATTRIBUTES DLLEXPORT :: get_xtalname
!DEC$ ATTRIBUTES DLLEXPORT :: set_xtalname
!DEC$ ATTRIBUTES DLLEXPORT :: get_datafile
!DEC$ ATTRIBUTES DLLEXPORT :: set_datafile
!DEC$ ATTRIBUTES DLLEXPORT :: get_mode
!DEC$ ATTRIBUTES DLLEXPORT :: set_mode

! the constructor routine for this class 
interface kinematical_T
  module procedure kinematical_constructor
end interface kinematical_T

contains

!--------------------------------------------------------------------------
type(kinematical_T) function kinematical_constructor( nmlfile ) result(kinematical)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! constructor for the kinematical_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call kinematical%readNameList(nmlfile)

end function kinematical_constructor

!--------------------------------------------------------------------------
subroutine kinematical_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! destructor for the kinematical_T Class
 
IMPLICIT NONE

type(kinematical_T), INTENT(INOUT)  :: self 

call reportDestructor('kinematical_T')

end subroutine kinematical_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! read the namelist from an nml file for the kinematical_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(kinematical_T), INTENT(INOUT)  :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

real(kind=sgl)   :: dmin
real(kind=sgl)   :: thr
real(kind=sgl)   :: voltage
integer(kind=irg):: nx
character(fnlen) :: xtalname
character(fnlen) :: datafile
character(5)     :: mode

! define the IO namelist to facilitate passing variables to the program.
namelist /EMkinematical/ dmin, voltage, thr, xtalname, datafile, mode, nx

! set the input parameters to default values (except for xtalname, which must be present)
dmin = 0.05                    ! smallest d-spacing to include in dynamical matrix [nm]
thr = 1.0                      ! smallest |structurefactor|^2 to include
voltage = 30000.0              ! microscope voltage [V]
nx = 500                       ! semi edge length of square Lambert projection
datafile = 'undefined'         ! output file name
xtalname = 'undefined'         ! structure file name
mode = 'lines'                 ! default plot mode

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
  open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
  read(UNIT=dataunit,NML=EMkinematical)
  close(UNIT=dataunit,STATUS='keep')

! check for required entries
  if (trim(xtalname).eq.'undefined') then
    call Message%printError('readNameList:','  crystal structure file name is undefined in '//nmlfile)
  end if
  if (trim(datafile).eq.'undefined') then
    call Message%printError('readNameList:',' output file name is undefined in '//nmlfile)
  end if
end if

! if we get here, then all appears to be ok, and we need to fill in the emnl fields
self%nml%dmin = dmin
self%nml%thr = thr
self%nml%voltage = voltage
self%nml%nx = nx
self%nml%xtalname = xtalname
self%nml%datafile = datafile
self%nml%mode = mode 

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! pass the namelist for the kinematical_T Class to the calling program

IMPLICIT NONE 

class(kinematical_T), INTENT(INOUT)          :: self
type(kinematicalNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(kinematical_T), INTENT(INOUT)     :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T),INTENT(INOUT)          :: HDFnames

integer(kind=irg),parameter             :: n_int = 1, n_real = 3
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( knl => self%nml )

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! integers 
io_int = (/ knl%nx /)
intlist(1) = 'nx'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ knl%voltage, knl%dmin, knl%thr /)
reallist(1) = 'voltage'
reallist(2) = 'dmin'
reallist(3) = 'thr'
call HDF%writeNMLreals(io_real, reallist, n_real)

! write all the strings
dataset = SC_xtalname
line2(1) = knl%xtalname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create xtalname dataset',hdferr)

dataset = SC_datafile
line2(1) = knl%datafile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create datafile dataset',hdferr)

dataset = 'mode'
line2(1) = knl%mode
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create mode dataset',hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_


!--------------------------------------------------------------------------
function get_dmin_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! get dmin from the kinematical_T class

IMPLICIT NONE 

class(kinematical_T), INTENT(INOUT)     :: self
real(kind=sgl)                          :: out

out = self%nml%dmin

end function get_dmin_

!--------------------------------------------------------------------------
subroutine set_dmin_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! set dmin in the kinematical_T class

IMPLICIT NONE 

class(kinematical_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)              :: inp

self%nml%dmin = inp

end subroutine set_dmin_

!--------------------------------------------------------------------------
function get_thr_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! get thr from the kinematical_T class

IMPLICIT NONE 

class(kinematical_T), INTENT(INOUT)     :: self
real(kind=sgl)                          :: out

out = self%nml%thr

end function get_thr_

!--------------------------------------------------------------------------
subroutine set_thr_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! set thr in the kinematical_T class

IMPLICIT NONE 

class(kinematical_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)              :: inp

self%nml%thr = inp

end subroutine set_thr_

!--------------------------------------------------------------------------
function get_voltage_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! get voltage from the kinematical_T class

IMPLICIT NONE 

class(kinematical_T), INTENT(INOUT)     :: self
real(kind=sgl)                          :: out

out = self%nml%voltage

end function get_voltage_

!--------------------------------------------------------------------------
subroutine set_voltage_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! set voltage in the kinematical_T class

IMPLICIT NONE 

class(kinematical_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)              :: inp

self%nml%voltage = inp

end subroutine set_voltage_

!--------------------------------------------------------------------------
function get_nx_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! get nx from the kinematical_T class

IMPLICIT NONE 

class(kinematical_T), INTENT(INOUT)     :: self
integer(kind=irg)                       :: out

out = self%nml%nx

end function get_nx_

!--------------------------------------------------------------------------
subroutine set_nx_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! set nx in the kinematical_T class

IMPLICIT NONE 

class(kinematical_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)           :: inp

self%nml%nx = inp

end subroutine set_nx_

!--------------------------------------------------------------------------
function get_xtalname_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! get xtalname from the kinematical_T class

IMPLICIT NONE 

class(kinematical_T), INTENT(INOUT)     :: self
character(fnlen)                        :: out

out = self%nml%xtalname

end function get_xtalname_

!--------------------------------------------------------------------------
subroutine set_xtalname_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! set xtalname in the kinematical_T class

IMPLICIT NONE 

class(kinematical_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)            :: inp

self%nml%xtalname = inp

end subroutine set_xtalname_

!--------------------------------------------------------------------------
function get_datafile_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! get datafile from the kinematical_T class

IMPLICIT NONE 

class(kinematical_T), INTENT(INOUT)     :: self
character(fnlen)                        :: out

out = self%nml%datafile

end function get_datafile_

!--------------------------------------------------------------------------
subroutine set_datafile_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! set datafile in the kinematical_T class

IMPLICIT NONE 

class(kinematical_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)            :: inp

self%nml%datafile = inp

end subroutine set_datafile_

!--------------------------------------------------------------------------
function get_mode_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! get mode from the kinematical_T class

IMPLICIT NONE 

class(kinematical_T), INTENT(INOUT)     :: self
character(5)                            :: out

out = self%nml%mode

end function get_mode_

!--------------------------------------------------------------------------
subroutine set_mode_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! set mode in the kinematical_T class

IMPLICIT NONE 

class(kinematical_T), INTENT(INOUT)     :: self
character(5), INTENT(IN)                :: inp

self%nml%mode = inp

end subroutine set_mode_

!--------------------------------------------------------------------------
subroutine kinematical_(self, EMsoft, progname)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! perform the computations

use mod_EMsoft

IMPLICIT NONE 

class(kinematical_T), INTENT(INOUT)     :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 

associate( knl=>self%nml )

! call the main routine 
if (knl%mode.eq.'lines') call self%KinematicalLines(EMsoft, progname)

if (knl%mode.eq.'bands') call self%KinematicalBands(EMsoft, progname)

end associate

end subroutine kinematical_

!--------------------------------------------------------------------------
subroutine KinematicalLines(self, EMsoft, progname)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! generate a master pattern file with Kossel lines

use mod_EMsoft
use mod_initializers
use mod_crystallography
use mod_diffraction
use mod_symmetry
use mod_io
use mod_rotations
use mod_quaternions
use mod_timing
use mod_Lambert
use mod_math
use HDF5
use mod_HDFsupport
use mod_HDFnames
use ISO_C_BINDING
use stringconstants

IMPLICIT NONE 

class(kinematical_T),INTENT(INOUT)  :: self
type(EMsoft_T),INTENT(INOUT)        :: EMsoft
character(fnlen),INTENT(IN)         :: progname

type(HDF_T)                         :: HDF 
type(HDFnames_T)                    :: HDFnames 
type(IO_T)                          :: Message 
type(cell_T)                        :: cell 
type(SpaceGroup_T)                  :: SG 
type(diffraction_T)                 :: Diff 
type(DynType)                       :: Dyn
type(Timing_T)                      :: timer
type(Quaternion_T)                  :: qu
type(q_T)                           :: qq 
type(a_T)                           :: ax
type(Lambert_T)                     :: Lambert

character(fnlen)                    :: datafile, groupname, dataset, xtalname
logical                             :: f_exists, readonly, verbose
integer(kind=irg)                   :: hdferr, nlines, i, istat, ix, iy, nx
integer(HSIZE_T)                    :: dims3(3), dims4(4)
real(kind=dbl)                      :: EkeV, mLambda
real(kind=sgl)                      :: m, xyzs(3)

integer(kind=irg)                   :: imh, imk, iml, ii, j, num, nums, mhkl, numphi, ierr, io_int(3)
integer(kind=irg),allocatable       :: gvec(:,:)
integer(kind=irg)                   :: h,k,l,totfam,ind(3),icnt, oi_int(1), itmp(48,3), g1(3), g2(3)
logical                             :: first
real(kind=sgl)                      :: g(3), thr, dphi, gc(3), gax(3), gz(3), v(3), x, gg, sgn, &
                                       xy(2), xyz(3)
real(kind=sgl),allocatable          :: Vgg(:),th(:), unitvec(:,:), scaledVgg(:)
character(1)                        :: space
real(kind=sgl)                      :: dhkl, Radius, tav
real(kind=sgl)                      :: ixy(2),scl
real(kind=sgl)                      :: dx,dy,dxm,dym, sa, ca
integer(kind=irg)                   :: jj,kk
integer(kind=irg)                   :: nix,niy,nixp,niyp
logical,allocatable                 :: keep(:)
character(11)                       :: dstr
character(15)                       :: tstrb
character(15)                       :: tstre

integer(kind=irg),allocatable       :: acc_e(:,:,:)
real(kind=sgl),allocatable          :: Eweights(:)
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)

! arrays for master pattern and stereographic projections
real(kind=sgl),allocatable          :: masterNH(:,:), masterSH(:,:), stereoNH(:,:), stereoSH(:,:), &
                                       phi(:), cp(:), sp(:), dc(:,:), mLPNH(:,:,:), mLPSH(:,:,:)
type(gnode),save                    :: rlp

!$OMP THREADPRIVATE(rlp) 

verbose = .TRUE.
space = 'r'

call openFortranHDFInterface()

call setRotationPrecision('single')

associate( knl=>self%nml )

! set the HDF group names for this program
HDF = HDF_T() 
HDFnames = HDFnames_T() 

! initialize the timing routines
timer = Timing_T()
tstrb = timer%getTimeString()
dstr  = timer%getDateString()

! initialize the crystal structure and compute a list of potential reflectors 
! get the crystal structure from the *.xtal file
call cell%setFileName(knl%xtalname)
call Diff%setrlpmethod('WK')

call Diff%setV(dble(knl%voltage))
call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, knl%dmin, verbose, useHDF=HDF)
mLambda = Diff%getWaveLength()

! generate a list of hkl indices 

! compute the range of reflections for the lookup table and allocate the table
! The master list is easily created by brute force
 imh = 1
 do 
   dhkl = 1.0/cell%CalcLength((/float(imh) ,0.0_sgl,0.0_sgl/), 'r')
   if (dhkl.lt.knl%dmin) EXIT
   imh = imh + 1
 end do
 imk = 1
 do 
   dhkl = 1.0/cell%CalcLength((/0.0_sgl,float(imk),0.0_sgl/), 'r')
   if (dhkl.lt.knl%dmin) EXIT
   imk = imk + 1
 end do
 iml = 1
 do 
   dhkl = 1.0/cell%CalcLength((/0.0_sgl,0.0_sgl,float(iml)/), 'r')
   if (dhkl.lt.knl%dmin) EXIT
   iml = iml + 1
 end do

io_int(1:3) = (/ imh, imk, iml /)
call Message%WriteValue(' reflection range : ',io_int,3)

! allocate all arrays
 ii = (2*imh+1)*(2*imk+1)*(2*iml+1)
 allocate(gvec(3,ii))
 allocate(Vgg(ii))
 allocate(unitvec(3,ii))
 allocate(th(ii))

! loop through all (hkl) values
 first = .TRUE.
 icnt = 1
 do h=-imh,imh
  ind(1)=-h
  do k=-imk,imk
   ind(2)=-k
   do l=-iml,iml
     ind(3)=-l

! compute the Fourier coefficient
     call Diff%CalcUcg(cell, ind)
     rlp = Diff%getrlp() 

! ignore the reciprocal lattice point if Vgg is small
     if (rlp%Vmod.ge.knl%thr) then 

! store the Fourier coefficient of the lattice potential
      if (abs(h)+abs(k)+abs(l).eq.0) then
        Vgg(icnt) = 0.0
      else
        Vgg(icnt)=rlp%Vmod**2
      end if
! g-vector (Miller indices)
      gvec(1:3,icnt) = (/ h, k, l /)
! unit vector in cartesian basis
      g(1:3)=float( (/ h, k, l /) )
      gg=cell%CalcLength(g,'r')
      call cell%TransSpace(g,gc,'r','c')
      gc = gc/sqrt(sum(gc*gc))
      unitvec(1:3,icnt) = gc(1:3) 
! and Bragg angle
      th(icnt)=asin(0.5*mLambda*gg)

! increment counter
      icnt=icnt+1
    end if
   end do
  end do
 end do

 icnt=icnt-1
 oi_int(1)=icnt
 call Message%WriteValue(' Total number of entries found = ', oi_int, 1, "(I6)")

 allocate(scaledVgg(icnt))
do i=1,icnt
  if (Vgg(i).ne.0.0) scaledVgg(i) = (Vgg(i))**0.25
end do

! sort these points by increasing Bragg angle


! next, we need to generate a master pattern that has a line with Gaussian inverted
! profile for each of the Kossel lines.  We'll use the standard square Lambert projection
! to create this pattern, and then convert it also to stereographic projections.
nx = self%get_nx()
allocate(masterNH(-nx:nx,-nx:nx),stat=istat)
allocate(masterSH(-nx:nx,-nx:nx),stat=istat)
masterNH = 1.0
masterSH = 1.0

numphi = 360*16
allocate(phi(numphi), cp(numphi), sp(numphi))
allocate( dc(3,numphi) )
dphi = 2.0*cPi/dble(numphi)

phi = (/ (float(i-1)*dphi,i=1,numphi) /)
cp = cos(phi)
sp = sin(phi)
gz = (/ 0.0, 0.0, 1.0 /)
scl = float(nx)

do kk=1,icnt   ! ignore the last point
  k = kk ! ksort(kk)   ! pick them in the correct order

! get the unrotated direction cosines of the sampling points on the sphere; this generate a circle on the sphere (Kossel cone trace)
  ca = cos( th(k) )
  sa = sin( th(k) )
  v = unitvec(1:3,k) 
  sgn = 1.0
  x = cell%CalcDot(gz,v,'c')
  if (x.lt.0.0) sgn = -1.0
  do i=1,numphi
      dc(1,i) = ca * cp(i)
      dc(2,i) = ca * sp(i)
      dc(3,i) = sa 
  end do

! then determine the rotation quaternion to bring the z axis onto the g direction (cartesian)
  if (x.ne.1.0) then   ! the cross product exists
    call cell%CalcCross(v,gz,gax,'c','c',0) ! gax is the rotation axis
    gax = gax/sqrt(sum(gax*gax))
    x = acos(x)
    if (x.lt.0.0) then
      ax = a_T( ainp = (/-gax(1),-gax(2),-gax(3), -x /) )
    else
      ax = a_T( ainp = (/ gax(1), gax(2), gax(3), x /) )
    end if
    qq = ax%aq()
    call qu%set_quats( qq%q_copy() )
    qu = conjg(qu)
    do i=1,numphi
      v = qu%quat_Lp( dc(:,i) )
      v = v / sqrt(sum(v*v))
      dc(:,i) = v(:)
    end do
  end if

! ok, so we have rotated this set of directions; next we need to add some intensity to each 
! corresponding point in the master pattern
  do i=1,numphi
! convert these direction cosines to coordinates in the Rosca-Lambert projection
    call Lambert%setxyz( dc(:,i) )
    istat = Lambert%LambertSphereToSquare( ixy )
    if (istat .ne. 0) stop 'Something went wrong during interpolation...'
    ixy = ixy * scl
    nix = int(nx+ixy(1))-nx
    niy = int(nx+ixy(2))-nx
    if (dc(3,i) .ge. 0.0) then
      call AntiAlias(masterNH,ixy,nix,niy,nx,scaledVgg(k))
    else
      call AntiAlias(masterSH,ixy,nix,niy,nx,scaledVgg(k))
    end if
  end do
end do

! make sure the master pattern intensities are positive
x = minval(masterNH)
ca = minval(masterSH)
masterNH = masterNH - minval( (/x, ca/) )
masterNH = masterNH/maxval(masterNH)
masterSH = masterSH - minval( (/x, ca/) )
masterSH = masterSH/maxval(masterSH)

! get stereographic projections
allocate(stereoNH(-nx:nx,-nx:nx),stat=istat)
allocate(stereoSH(-nx:nx,-nx:nx),stat=istat)
tav = sum(masterNH)/float(2*nx+1)**2
Radius = 1.0
do i=-nx,nx 
  do j=-nx,nx 
    Lambert = Lambert_T( xy = (/ float(i), float(j) /) / float(knl%nx) )
    ierr = Lambert%StereoGraphicInverse( xyz, Radius )
    xyz = xyz/vecnorm(xyz)
    if (ierr.ne.0) then 
      stereoNH(i,j) = tav
      stereoSH(i,j) = tav
    else
      stereoNH(i,j) = InterpolateLambert(xyz, masterNH, nx)
      stereoSH(i,j) = InterpolateLambert(xyz, masterNH, nx)
    end if
  end do
end do

! prepare the output
call timer%makeTimeStamp()
dstr = timer%getDateString() 
tstre = timer%getTimeString()

!------------------------------
! write the output to an HDF5 file
!------------------------------
call Message%printMessage('opening '//trim(knl%datafile), frm = "(A)" )
call HDFnames%set_ProgramData(SC_EMkinematical) 
call HDFnames%set_NMLlist(SC_EMkinematicalNameList) 
call HDFnames%set_NMLfilename(SC_EMkinematicalNML) 

! open the output file
datafile = EMsoft%generateFilePath('EMdatapathname',trim(knl%datafile))

! open the file using the default properties.
hdferr =  HDF%createFile(datafile)

! write the EMheader to the file
groupname = trim(HDFnames%get_ProgramData())
call HDF%writeEMheader(dstr, tstrb, tstre, progname, groupname)

! add the CrystalData group at the top level of the file
call cell%addXtalDataGroup(SG, EMsoft, HDF)

! open or create a namelist group to write all the namelist files into
hdferr = HDF%createGroup(HDFnames%get_NMLfiles())

! read the text file and write the array to the file
dataset = trim(HDFnames%get_NMLfilename())
hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%nmldeffile)

! leave this group
call HDF%pop()

! create a namelist group to write all the namelist files into
hdferr = HDF%createGroup(HDFnames%get_NMLparameters())

call self%writeHDFNameList(HDF, HDFnames)

! leave this group
call HDF%pop()

! then the remainder of the data in a EMData group
hdferr = HDF%createGroup(HDFnames%get_EMData())
hdferr = HDF%createGroup(HDFnames%get_ProgramData())

dataset = SC_numg
hdferr = HDF%writeDatasetInteger(dataset, icnt)

dataset = SC_hklmax
hdferr = HDF%writeDatasetIntegerArray(dataset, (/ imh, imk, iml /), 3)

dataset = SC_hklarray
hdferr = HDF%writeDatasetIntegerArray(dataset, gvec(1:icnt,1:3), icnt, 3)

dataset = SC_modFsquared
hdferr = HDF%writeDatasetFloatArray(dataset, Vgg(1:icnt), icnt)

dataset = SC_unitGvectors
hdferr = HDF%writeDatasetFloatArray(dataset, unitvec(1:icnt, 1:3), icnt, 3)

dataset = SC_BraggAngle
hdferr = HDF%writeDatasetFloatArray(dataset, th(1:icnt), icnt)

dataset = SC_masterNH
hdferr = HDF%writeDatasetFloatArray(dataset, masterNH, 2*nx+1, 2*nx+1)

dataset = SC_masterSH
hdferr = HDF%writeDatasetFloatArray(dataset, masterSH, 2*nx+1, 2*nx+1)

dataset = SC_stereoNH
hdferr = HDF%writeDatasetFloatArray(dataset, stereoNH, 2*nx+1, 2*nx+1)

dataset = SC_stereoSH
hdferr = HDF%writeDatasetFloatArray(dataset, stereoSH, 2*nx+1, 2*nx+1)

! and close everything
call HDF%pop(.TRUE.)

call Message%printMessage(' Output data stored in file '//trim(datafile))

end associate

end subroutine KinematicalLines

!--------------------------------------------------------------------------
subroutine KinematicalBands(self, EMsoft, progname)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! generate a master pattern file with Kossel bands

use mod_EMsoft
use mod_initializers
use mod_crystallography
use mod_diffraction
use mod_symmetry
use mod_io
use mod_rotations
use mod_quaternions
use mod_timing
use mod_Lambert
use mod_math
use HDF5
use mod_HDFsupport
use mod_HDFnames
use ISO_C_BINDING
use stringconstants

IMPLICIT NONE 

class(kinematical_T),INTENT(INOUT)  :: self
type(EMsoft_T),INTENT(INOUT)        :: EMsoft
character(fnlen),INTENT(IN)         :: progname

type(HDF_T)                         :: HDF 
type(HDFnames_T)                    :: HDFnames 
type(IO_T)                          :: Message 
type(cell_T)                        :: cell 
type(SpaceGroup_T)                  :: SG 
type(diffraction_T)                 :: Diff 
type(DynType)                       :: Dyn
type(Timing_T)                      :: timer
type(Quaternion_T)                  :: qu
type(q_T)                           :: qq 
type(a_T)                           :: ax
type(Lambert_T)                     :: Lambert  

character(fnlen)                    :: listfile, masterfile, groupname, dataset, xtalname, outputfile, infile, datafile
logical                             :: f_exists, readonly, verbose
integer(kind=irg)                   :: hdferr, nlines, i, istat, ix, iy, nx, io_int(3), nkeep, ierr
integer(HSIZE_T)                    :: dims3(3), dims4(3)
real(kind=dbl)                      :: EkeV, mLambda
real(kind=sgl)                      :: m

integer(kind=irg)                   :: imh, imk, iml, ii, j, num, nums, mhkl, valpos, numphi,numtheta,iequiv(3,48),nequiv
integer(kind=irg),allocatable       :: family(:,:,:),numfam(:),idx(:), idx2(:), sfi(:), itmp(:,:)
integer(kind=irg)                   :: h,k,l,totfam,ind(3),icnt, oi_int(1), g1(3), g2(3), NUMTHREADS, TID
logical                             :: first
logical,allocatable                 :: z(:,:,:)
real(kind=sgl)                      :: g(3), thr, dphi, gc(3), gax(3), gz(3), v(3), x, val1, val2, &
                                       ang, valmax, edge, dc(3), Radius, xy(2), xyz(3)
real(kind=sgl),allocatable          :: Vgg(:),ddg(:),gg(:),th(:), gcart(:,:,:), gcrys(:,:), cp(:), sp(:), ca(:), sa(:), &
                                       phi(:), theta(:), Vg(:), VgX(:), VggX(:), cosnorm(:)
character(1)                        :: space
real(kind=sgl)                      :: dhkl, incrad, glen, sd, corr
real(kind=sgl)                      :: ixy(2),scl, cp2, dp
real(kind=sgl)                      :: dx,dy,dxm,dym, tav
integer(kind=irg)                   :: jj,kk
integer(kind=irg)                   :: nix,niy,nixp,niyp
logical,allocatable                 :: keep(:)
character(11)                       :: dstr
character(15)                       :: tstrb
character(15)                       :: tstre
integer(kind=irg),allocatable       :: acc_e(:,:,:)
real(kind=sgl),allocatable          :: Eweights(:)
real(kind=sgl),allocatable          :: srtmp(:,:,:,:), mLPNH(:,:,:), mLPSH(:,:,:), masterNH(:,:), masterSH(:,:), &
                                       KBI(:), kinmasterNH(:,:), kinmasterSH(:,:), kinNH(:,:), kinSH(:,:), &
                                       stereoNH(:,:), stereoSH(:,:)
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)
type(gnode),save                    :: rlp

!$OMP THREADPRIVATE(rlp) 

verbose = .TRUE.
space = 'r'

call openFortranHDFInterface()

call setRotationPrecision('single')

associate( knl=>self%nml )

! set the HDF group names for this program
HDF = HDF_T() 
HDFnames = HDFnames_T() 

! initialize the timing routines
timer = Timing_T()
tstrb = timer%getTimeString()
dstr  = timer%getDateString()

! initialize the crystal structure and compute a list of potential reflectors 
! get the crystal structure from the *.xtal file
call cell%setFileName(knl%xtalname)
call Diff%setrlpmethod('WK')

call Diff%setV(dble(knl%voltage))
call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, knl%dmin, verbose, useHDF=HDF)
mLambda = Diff%getWaveLength()

! compute the range of reflections for the lookup table and allocate the table
! The master list is easily created by brute force
 imh = 1
 do 
   dhkl = 1.0/cell%CalcLength((/float(imh) ,0.0_sgl,0.0_sgl/), 'r')
   if (dhkl.lt.knl%dmin) EXIT
   imh = imh + 1
 end do
 imk = 1
 do 
   dhkl = 1.0/cell%CalcLength((/0.0_sgl,float(imk),0.0_sgl/), 'r')
   if (dhkl.lt.knl%dmin) EXIT
   imk = imk + 1
 end do
 iml = 1
 do 
   dhkl = 1.0/cell%CalcLength((/0.0_sgl,0.0_sgl,float(iml)/), 'r')
   if (dhkl.lt.knl%dmin) EXIT
   iml = iml + 1
 end do

io_int(1:3) = (/ imh, imk, iml /)
call Message%WriteValue(' reflection range : ',io_int,3)

! allocate all arrays
 allocate(z(-2*imh:2*imh,-2*imk:2*imk,-2*iml:2*iml))
 ii = (2*imh+1)*(2*imk+1)*(2*iml+1)
 allocate(family(ii,48,3))
 allocate(numfam(ii))
 allocate(Vgg(ii))
 allocate(VggX(ii))
 allocate(ddg(ii))
 allocate(gg(ii))
 allocate(th(ii))
! determine the families of reflections with (hkl)<=(imh,imk,iml)
! first initialize the boolean array z
 z = .FALSE.
! then loop through all (hkl) values
 first = .TRUE.
 icnt = 1
 totfam=0
 rlp%method= 'DT'   ! we're computing simple Doyle-Turner or Smith-Burge scattering factors to get the list of reflectors
 do h=-imh,imh
  ind(1)=-h
  do k=-imk,imk
   ind(2)=-k
   do l=-iml,iml
    ind(3)=-l

! make sure we have not already done this one in another family
    if (.not.z(-h,-k,-l)) then

! if it is a new one, then determine the entire family
      call Diff%CalcUcg(cell, ind)
     rlp = Diff%getrlp() 

! but ignore the reciprocal lattice point if Vgg is small
     if (abs(rlp%Ucg).ge.thr) then 
      Vgg(icnt) = (abs(rlp%Ucg))**2

! copy family in array and label all its members and multiples in z-array
      call SG%CalcFamily(ind,num,space,itmp)
      do i=1,num
        family(icnt,i,1:3)=itmp(i,1:3)
        z(itmp(i,1),itmp(i,2),itmp(i,3))=.TRUE.
      end do

! increment family counter
      numfam(icnt)=num
      totfam=totfam+num-1
      icnt=icnt+1
     end if
    end if
   end do
  end do
 end do

 icnt=icnt-1
 oi_int(1)=icnt
 call Message%WriteValue(' Total number of families        = ', oi_int, 1, "(I6)")
 oi_int(1)=totfam
 call Message%WriteValue(' Total number of family members  = ', oi_int, 1, "(I6)")

 icnt = icnt-1   ! to eliminate the origin ... 
! compute d-spacings, g-spacings, theta
 allocate(gcart(3,icnt,48))
 do k=1,icnt
  g(1:3)=float(family(k,1,1:3))
  gg(k)=cell%CalcLength(g,'r')
  th(k)=asin(0.5*mLambda*gg(k))
  do j=1,numfam(k)
    g(1:3)=float(family(k,j,1:3))
    gg(k)=cell%CalcLength(g,'r')
    call cell%TransSpace(g,gc,'r','c')
    gc = gc/sqrt(sum(gc*gc))
    gcart(1:3,k,j) = gc(1:3) 
  end do
 end do

! we're taking the dot product between the plane normals and the direction cosine
! vectors, so we will need the complement of the theta angles ... 
cp2 = sngl(cPi/2.D0)
! th = cp2 - th 

! next, we scan over the entire Lambert square and, for each point, 
! determine which Kikuchi bands it belongs to; then we add all those 
! intensities together and place that value in the array. This should
! automatically take care of overlapping bands etc...
nx = self%get_nx()
allocate(kinmasterNH(-nx:nx,-nx:nx), kinmasterSH(-nx:nx,-nx:nx))
kinmasterNH = 0.0 
kinmasterSH = 0.0 

edge = 1.0 / float(nx)

do ix=-nx,nx 
  do iy=-nx,nx 
! get the direction cosines for this point 
    call Lambert%setxy( xy = (/ float(ix), float(iy) /) * edge )
    ierr = Lambert%LambertSquareToSphere( dc )

! loop over all family members 
    do k=1,icnt
      do j=1,numfam(k) 
        dp = DOT_PRODUCT( gcart(1:3,k,j), dc(1:3) )
          ang = acos( dp ) - cp2
          if  (abs(ang).le.th(k)) then 
            kinmasterNH(ix,iy) = kinmasterNH(ix,iy) + Vgg(k)
            kinmasterSH(-ix,-iy) = kinmasterSH(-ix,-iy) + Vgg(k)
          end if 
      end do 
    end do 
  end do 
end do 

! we need to correct the edge intensities because of some strange double counting 
! that occurs for the perimeter of the square Lambert projection ... 
kinmasterNH(-nx,:) = kinmasterNH(-nx+1,:) 
kinmasterNH( nx,:) = kinmasterNH( nx-1,:) 
kinmasterNH( -nx+1:nx-1,-nx) = kinmasterNH(-nx+1:nx-1,-nx+1) 
kinmasterNH( -nx+1:nx-1, nx) = kinmasterNH(-nx+1:nx-1, nx-1)
kinmasterNH(-nx,-nx) = kinmasterNH(-nx+1,-nx+1)
kinmasterNH( nx,-nx) = kinmasterNH( nx-1,-nx+1)
kinmasterNH(-nx, nx) = kinmasterNH(-nx+1, nx-1)
kinmasterNH( nx, nx) = kinmasterNH( nx-1, nx-1)

kinmasterSH(-nx,:) = kinmasterSH(-nx+1,:)
kinmasterSH( nx,:) = kinmasterSH( nx-1,:)
kinmasterSH( -nx+1:nx-1,-nx) = kinmasterSH(-nx+1:nx-1,-nx+1) 
kinmasterSH( -nx+1:nx-1, nx) = kinmasterSH(-nx+1:nx-1, nx-1)
kinmasterSH(-nx,-nx) = kinmasterSH(-nx+1,-nx+1)
kinmasterSH( nx,-nx) = kinmasterSH( nx-1,-nx+1)
kinmasterSH(-nx, nx) = kinmasterSH(-nx+1, nx-1)
kinmasterSH( nx, nx) = kinmasterSH( nx-1, nx-1)

allocate(stereoNH(-nx:nx,-nx:nx),stat=istat)
allocate(stereoSH(-nx:nx,-nx:nx),stat=istat)
! get stereographic projections
Radius = 1.0
tav = sum(kinmasterNH)/float(2*nx+1)**2
do i=-nx,nx 
  do j=-nx,nx 
    Lambert = Lambert_T( xy = (/ float(i), float(j) /) / float(knl%nx) )
    ierr = Lambert%StereoGraphicInverse( xyz, Radius )
    xyz = xyz/vecnorm(xyz)
    if (ierr.ne.0) then 
      stereoNH(i,j) = tav
      stereoSH(i,j) = tav
    else
      stereoNH(i,j) = InterpolateLambert(xyz, kinmasterNH, nx)
      stereoSH(i,j) = InterpolateLambert(xyz, kinmasterSH, nx)
    end if
  end do
end do

! then save this in an HDF file 
call timer%makeTimeStamp()
dstr = timer%getDateString() 
tstre = timer%getTimeString()

!------------------------------
! write the output to an HDF5 file
!------------------------------
call Message%printMessage('opening '//trim(knl%datafile), frm = "(A)" )
call HDFnames%set_ProgramData(SC_EMkinematical) 
call HDFnames%set_NMLlist(SC_EMkinematicalNameList) 
call HDFnames%set_NMLfilename(SC_EMkinematicalNML) 

! open the output file
datafile = EMsoft%generateFilePath('EMdatapathname',trim(knl%datafile))

! open the file using the default properties.
hdferr =  HDF%createFile(datafile)

! write the EMheader to the file
groupname = trim(HDFnames%get_ProgramData())
call HDF%writeEMheader(dstr, tstrb, tstre, progname, groupname)

! add the CrystalData group at the top level of the file
call cell%addXtalDataGroup(SG, EMsoft, HDF)

! open or create a namelist group to write all the namelist files into
hdferr = HDF%createGroup(HDFnames%get_NMLfiles())

! read the text file and write the array to the file
dataset = trim(HDFnames%get_NMLfilename())
hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%nmldeffile)

! leave this group
call HDF%pop()

! create a namelist group to write all the namelist files into
hdferr = HDF%createGroup(HDFnames%get_NMLparameters())

call self%writeHDFNameList(HDF, HDFnames)

! leave this group
call HDF%pop()

! then the remainder of the data in a EMData group
hdferr = HDF%createGroup(HDFnames%get_EMData())
hdferr = HDF%createGroup(HDFnames%get_ProgramData())

dataset = SC_numg
hdferr = HDF%writeDatasetInteger(dataset, icnt)

dataset = SC_hklmax
hdferr = HDF%writeDatasetIntegerArray(dataset, (/ imh, imk, iml /), 3)

! dataset = SC_hklarray
! hdferr = HDF%writeDatasetIntegerArray(dataset, gvec(1:icnt,1:3), icnt, 3)

dataset = SC_modFsquared
hdferr = HDF%writeDatasetFloatArray(dataset, Vgg(1:icnt), icnt)

! dataset = SC_unitGvectors
! hdferr = HDF%writeDatasetFloatArray(dataset, unitvec(1:icnt, 1:3), icnt, 3)

dataset = SC_BraggAngle
hdferr = HDF%writeDatasetFloatArray(dataset, th(1:icnt), icnt)

dataset = SC_masterNH
hdferr = HDF%writeDatasetFloatArray(dataset, kinmasterNH, 2*nx+1, 2*nx+1)

dataset = SC_masterSH
hdferr = HDF%writeDatasetFloatArray(dataset, kinmasterSH, 2*nx+1, 2*nx+1)

dataset = SC_stereoNH
hdferr = HDF%writeDatasetFloatArray(dataset, stereoNH, 2*nx+1, 2*nx+1)

dataset = SC_stereoSH
hdferr = HDF%writeDatasetFloatArray(dataset, stereoSH, 2*nx+1, 2*nx+1)

! and close everything
call HDF%pop(.TRUE.)

call Message%printMessage(' Output data stored in file '//trim(datafile))

end associate

end subroutine KinematicalBands

!--------------------------------------------------------------------------
subroutine AntiAlias(master,ixy,nix,niy,nx,inten)
!! author: MDG 
!! version: 1.0 
!! date: 03/26/20
!!
!! put intensities in array using simple anti-aliasing solution to reduce jaggies
!!
!! we may want to replace this by a better algorithm, e.g.,  Xiaolin Wu's anti-aliasing 
!! line drawing algorithm

IMPLICIT NONE

integer(kind=irg),INTENT(IN)    :: nx
real(kind=sgl),INTENT(INOUT)    :: master(-nx:nx,-nx:nx)
real(kind=sgl),INTENT(IN)       :: ixy(2)
integer(kind=irg),INTENT(IN)    :: nix
integer(kind=irg),INTENT(IN)    :: niy
real(kind=sgl),INTENT(IN)       :: inten

real(kind=sgl)                  :: dx, dy, d

! look for the nearest pixel and share the intensity with that pixel
! only for the interior pixels of the pattern; to make things faster,
! we'll use the simple Manhattan distance w.r.t. nearest neighbors and
! apply the lever rule to set the intensities

if ((abs(nix).lt.nx).and.(abs(niy).lt.nx)) then 
  dx = ixy(1)-nix
  dy = ixy(2)-niy
  if (abs(dx).gt.abs(dy)) then
    d = abs(dx)
    if (dx.lt.0.0) then  
      master(nix-1,niy) = master(nix-1,niy) - d * inten
      master(nix,niy) = master(nix,niy) - (1.0-d) * inten
    else
      master(nix,niy) = master(nix,niy) - (1.0-d) * inten
      master(nix+1,niy) = master(nix+1,niy) - d * inten
    end if
  else
    d = abs(dy)
    if (dy.lt.0.0) then  
      master(nix,niy-1) = master(nix,niy-1) - d * inten
      master(nix,niy) = master(nix,niy) - (1.0-d) * inten
    else
      master(nix,niy) = master(nix,niy) - (1.0-d) * inten
      master(nix,niy+1) = master(nix,niy+1) - d * inten
    end if
  end if
end if

end subroutine AntiAlias





end module mod_kinematical