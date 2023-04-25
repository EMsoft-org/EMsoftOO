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

module mod_ISEmaster
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/17/22
  !!
  !! class definition for the EMISEmaster program

use mod_kinds
use mod_global

IMPLICIT NONE 
private 

type, private :: poslist 
  real(kind=dbl)          :: xyz(3)
  real(kind=dbl)          :: radius
  type(poslist), pointer  :: next
end type poslist 

! namelist for the EMISEmaster program
type, public :: ISEmasterNameListType
  integer(kind=irg)       :: npx
  integer(kind=irg)       :: nthreads
  real(kind=sgl)          :: iscale(3)
  real(kind=sgl)          :: multiplier
  character(fnlen)        :: outname
  character(fnlen)        :: tiffname
  character(fnlen)        :: xtalname
end type ISEmasterNameListType

! class definition
type, public :: ISEmaster_T
private 
  character(fnlen)              :: nmldeffile = 'EMISEmaster.nml'
  type(ISEmasterNameListType)   :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: ISEmaster_
  procedure, pass(self) :: get_npx_
  procedure, pass(self) :: get_multiplier_
  procedure, pass(self) :: get_nthreads_
  procedure, pass(self) :: get_iscale_
  procedure, pass(self) :: get_outname_
  procedure, pass(self) :: get_tiffname_
  procedure, pass(self) :: get_xtalname_
  procedure, pass(self) :: set_npx_
  procedure, pass(self) :: set_multiplier_
  procedure, pass(self) :: set_nthreads_
  procedure, pass(self) :: set_iscale_
  procedure, pass(self) :: set_outname_
  procedure, pass(self) :: set_tiffname_
  procedure, pass(self) :: set_xtalname_
  procedure, pass(self) :: getISEintensity_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: ISEmaster => ISEmaster_
  generic, public :: get_npx => get_npx_
  generic, public :: get_multiplier => get_multiplier_
  generic, public :: get_nthreads => get_nthreads_
  generic, public :: get_iscale => get_iscale_
  generic, public :: get_outname => get_outname_
  generic, public :: get_tiffname => get_tiffname_
  generic, public :: get_xtalname => get_xtalname_
  generic, public :: set_npx => set_npx_
  generic, public :: set_multiplier => set_multiplier_
  generic, public :: set_nthreads => set_nthreads_
  generic, public :: set_iscale => set_iscale_
  generic, public :: set_outname => set_outname_
  generic, public :: set_tiffname => set_tiffname_
  generic, public :: set_xtalname => set_xtalname_
  generic, public :: getISEintensity => getISEintensity_

end type ISEmaster_T

! the constructor routine for this class 
interface ISEmaster_T
  module procedure ISEmaster_constructor
end interface ISEmaster_T

contains

!--------------------------------------------------------------------------
type(ISEmaster_T) function ISEmaster_constructor( nmlfile ) result(ISEmaster)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! constructor for the ISEmaster_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call ISEmaster%readNameList(nmlfile)

end function ISEmaster_constructor

!--------------------------------------------------------------------------
subroutine ISEmaster_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! destructor for the ISEmaster_T Class
 
IMPLICIT NONE

type(ISEmaster_T), INTENT(INOUT)  :: self 

call reportDestructor('ISEmaster_T')

end subroutine ISEmaster_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! read the namelist from an nml file for the ISEmaster_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(ISEmaster_T), INTENT(INOUT)    :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg)       :: npx
integer(kind=irg)       :: nthreads
real(kind=sgl)          :: iscale(3)
real(kind=sgl)          :: multiplier
character(fnlen)        :: outname
character(fnlen)        :: tiffname
character(fnlen)        :: xtalname

! define the IO namelist to facilitate passing variables to the program.
namelist /ISEmastervars/ npx,nthreads,xtalname,outname,iscale,tiffname,multiplier 

! set the input parameters to default values (except for xtalname, which must be present)
npx = 500
nthreads = 1
iscale = (/ 2.0, 0.15, 2.0 /)
multiplier = 3.0
tiffname = 'undefined'
outname = 'undefined'
xtalname = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=ISEmastervars)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(xtalname).eq.'undefined') then
  call Message%printError('readNameList_:',' crystal structure file name is undefined in '//nmlfile)
 end if

 if (trim(outname).eq.'undefined') then
  call Message%printError('readNameList_:',' output file name is undefined in '//nmlfile)
 end if

  if (trim(tiffname).eq.'undefined') then
  call Message%printError('readNameList_:',' tiff file name is undefined in '//nmlfile)
 end if
end if

! if we get here, then all appears to be ok, and we need to fill in the nml fields
self%nml%npx = npx
self%nml%nthreads = nthreads
self%nml%iscale = iscale
self%nml%multiplier = multiplier
self%nml%outname= outname
self%nml%tiffname= tiffname
self%nml%xtalname = xtalname

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! pass the namelist for the ISEmaster_T Class to the calling program

IMPLICIT NONE 

class(ISEmaster_T), INTENT(INOUT)          :: self
type(ISEmasterNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(ISEmaster_T), INTENT(INOUT)       :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 2, n_real = 4
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( nml => self%nml )

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! write all the single integers
io_int = (/ nml%npx, nml%nthreads /)
intlist(1) = 'npx'
intlist(2) = 'nthreads'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ nml%iscale(1), nml%iscale(2), nml%iscale(3), nml%multiplier  /)
reallist(1) = 'a'
reallist(2) = 'b'
reallist(3) = 'c'
reallist(4) = 'multiplier'
call HDF%writeNMLreals(io_real, reallist, n_real)

! write all the strings
dataset = SC_outname
line2(1) = trim(nml%outname)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create outname dataset', hdferr)

dataset = SC_tiffname
line2(1) = trim(nml%tiffname)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create tiffname dataset', hdferr)

dataset = SC_xtalname
line2(1) = trim(nml%xtalname)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create xtalname dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
function get_npx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_npx_
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! get npx from the ISEmaster_T class

IMPLICIT NONE 

class(ISEmaster_T), INTENT(INOUT)     :: self
integer(kind=irg)                     :: out

out = self%nml%npx

end function get_npx_

!--------------------------------------------------------------------------
subroutine set_npx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_npx_
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! set npx in the ISEmaster_T class

IMPLICIT NONE 

class(ISEmaster_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%npx = inp

end subroutine set_npx_

!--------------------------------------------------------------------------
function get_multiplier_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_multiplier_
!! author: MDG 
!! version: 1.0 
!! date: 02/18/22
!!
!! get multiplier from the ISEmaster_T class

IMPLICIT NONE 

class(ISEmaster_T), INTENT(INOUT)     :: self
real(kind=sgl)                        :: out

out = self%nml%multiplier

end function get_multiplier_

!--------------------------------------------------------------------------
subroutine set_multiplier_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_multiplier_
!! author: MDG 
!! version: 1.0 
!! date: 02/18/22
!!
!! set multiplier in the ISEmaster_T class

IMPLICIT NONE 

class(ISEmaster_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)            :: inp

self%nml%multiplier = inp

end subroutine set_multiplier_

!--------------------------------------------------------------------------
function get_nthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nthreads_
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! get nthreads from the ISEmaster_T class

IMPLICIT NONE 

class(ISEmaster_T), INTENT(INOUT)     :: self
integer(kind=irg)                     :: out

out = self%nml%nthreads

end function get_nthreads_

!--------------------------------------------------------------------------
subroutine set_nthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nthreads_
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! set nthreads in the ISEmaster_T class

IMPLICIT NONE 

class(ISEmaster_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%nthreads = inp

end subroutine set_nthreads_

!--------------------------------------------------------------------------
function get_iscale_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_iscale_
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! get iscale from the ISEmaster_T class

IMPLICIT NONE 

class(ISEmaster_T), INTENT(INOUT)     :: self
real(kind=sgl)                        :: out(3)

out = self%nml%iscale

end function get_iscale_

!--------------------------------------------------------------------------
subroutine set_iscale_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_iscale_
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! set iscale in the ISEmaster_T class

IMPLICIT NONE 

class(ISEmaster_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)            :: inp(3)

self%nml%iscale = inp

end subroutine set_iscale_

!--------------------------------------------------------------------------
function get_outname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_outname_
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! get outname from the ISEmaster_T class

IMPLICIT NONE 

class(ISEmaster_T), INTENT(INOUT)     :: self
character(fnlen)                      :: out

out = self%nml%outname

end function get_outname_

!--------------------------------------------------------------------------
subroutine set_outname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_outname_
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! set outname in the ISEmaster_T class

IMPLICIT NONE 

class(ISEmaster_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%outname = inp

end subroutine set_outname_

!--------------------------------------------------------------------------
function get_tiffname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_tiffname_
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! get tiffname from the ISEmaster_T class

IMPLICIT NONE 

class(ISEmaster_T), INTENT(INOUT)     :: self
character(fnlen)                      :: out

out = self%nml%tiffname

end function get_tiffname_

!--------------------------------------------------------------------------
subroutine set_tiffname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_tiffname_
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! set tiffname in the ISEmaster_T class

IMPLICIT NONE 

class(ISEmaster_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%tiffname = inp

end subroutine set_tiffname_

!--------------------------------------------------------------------------
function get_xtalname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_xtalname_
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! get xtalname from the ISEmaster_T class

IMPLICIT NONE 

class(ISEmaster_T), INTENT(INOUT)     :: self
character(fnlen)                      :: out

out = self%nml%xtalname

end function get_xtalname_

!--------------------------------------------------------------------------
subroutine set_xtalname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_xtalname_
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! set xtalname in the ISEmaster_T class

IMPLICIT NONE 

class(ISEmaster_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%xtalname = inp

end subroutine set_xtalname_

!--------------------------------------------------------------------------
recursive function getISEintensity_(self, kloc, atomcnt, atomlist, atomrad, rsphere, a, b, usehex) result(inten)
!DEC$ ATTRIBUTES DLLEXPORT :: getISEintensity_
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! compute the ISE intensity (based on iCHORD forward model)

use mod_quaternions
use mod_rotations
use mod_others

IMPLICIT NONE

class(ISEmaster_T), INTENT(INOUT)   :: self
real(kind=dbl), INTENT(IN)          :: kloc(3)             ! this is already a unit vector
integer(kind=irg), INTENT(IN)       :: atomcnt 
real(kind=dbl), INTENT(IN)          :: atomlist(3,atomcnt) ! in Angstrom units
real(kind=dbl), INTENT(IN)          :: atomrad(atomcnt)    ! in Angstrom units
real(kind=dbl), INTENT(IN)          :: rsphere             ! in Angstrom units
real(kind=dbl), INTENT(IN)          :: a
real(kind=dbl), INTENT(IN)          :: b
logical,INTENT(IN),OPTIONAL         :: usehex
real(kind=dbl)                      :: inten 

type(a_T)                           :: ax
type(q_T)                           :: qu
type(Quaternion_T)                  :: quat

real(kind=dbl)                      :: apos(3,atomcnt), dp, xyz(3), Dsphere, Appix, shft, px, py, t, k(3), &
                                       dis(atomcnt), ndis(atomcnt), Gd(atomcnt), adisk, bdisk, arad(atomcnt)
real(kind=dbl),allocatable          :: pplane(:,:) 
integer(kind=irg)                   :: i, iatom, ss(atomcnt), spsize, ipx, ipy, ix, iy, maxrad

call setRotationPrecision('d')

! take care of the hexagonal reference frame
k = kloc
if (present(usehex)) then 
  if (usehex.eqv..TRUE.) then 
    k(2) = (k(1)+2.D0*k(2))/sqrt(3.D0)
  end if 
end if 
k = k/sqrt(sum(k*k))

! use the k vector to determine a rotation quaternion
dp = DOT_PRODUCT(k, (/ 0.D0, 0.D0, 1.D0 /) )
if (abs(dp).eq.1.D0) then ! no rotation needed
! just copy the atom coordinates
    apos = atomlist 
else 
    t = sqrt( k(1)*k(1) + k(2)*k(2))
    if (t.ne.0.D0) then 
      ax = a_T( adinp = (/ -k(2)/t, k(1)/t, 0.D0, acos( dp ) /) )
    else
      ax = a_T( adinp = (/ 0.D0, 0.D0, 1.D0, acos( dp ) /) )
    end if 
    qu = ax%aq()
    quat = Quaternion_T( qd = qu%q_copyd() )
    do i=1,atomcnt 
        apos(1:3,i) = quat%quat_Lp( atomlist(1:3,i) )
    end do
end if 

! then compute the distances to the projection plane (normal to the z-axis)
do i=1,atomcnt 
    dis(i) = rsphere - apos(3,i)
end do 

! sort in ascending order
call qsortd( dis, ss, atomcnt)
do i=1,atomcnt 
    ndis(i) = dis(ss(i))
end do 
dis = ndis
! and reverse the ndis array
do i=1,atomcnt
    ndis(i) = dis(atomcnt-i+1)
end do 

! get the brightness values 
Gd = a * ndis**b

! get the projection plane size
Dsphere = 2.D0 * rsphere
Appix = 0.25D0
spsize = int(Dsphere * 8.D0/Appix )
arad = (atomrad/Appix)**2

! initialize the projection plane
allocate(pplane(spsize,spsize)) 
pplane = 0.D0 

! get the largest disk radius in units of pixels ... 
maxrad = int(sqrt(maxval(arad))) + 1

! move all the atom x and y coordinates so that the center of the sphere projects on the center of the plane
shft = spsize*0.5D0
apos = 6.D0 * apos / Appix + shft 

! for all the atoms in the list, replicate the intensity in the projection plane in a small disk
! but only if the intensity is equal to zero to guarantee that the brightest disk is kept
do iatom=1,atomcnt 
  ! get the projected position in pixel units (nearest point) 
  px = apos(1,ss(iatom))
  py = apos(2,ss(iatom))
  ipx = nint(px)
  ipy = nint(py)
  do ix = ipx-maxrad,ipx+maxrad
    do iy = ipy-maxrad,ipy+maxrad
      if ( ((ix-px)**2+(iy-py)**2).le.arad(iatom) ) then
       ! if (pplane(ix,iy).eq.0.0) then
       if (Gd(iatom).gt.pplane(ix,iy)) then
         pplane(ix,iy)=Gd(iatom)
       end if
      end if
    end do
  end do
end do

! and sum all the values to get the ISE intensity 
inten = sum(pplane)

deallocate(pplane)

end function getISEintensity_

!--------------------------------------------------------------------------
subroutine ISEmaster_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: ISEmaster_
!! author: MDG 
!! version: 1.0 
!! date: 02/17/22
!!
!! perform the computations

use mod_EMsoft
use mod_initializers
use mod_memory
use mod_symmetry
use mod_crystallography
use mod_kvectors
use mod_diffraction
use mod_io
use mod_math
use mod_timing
use mod_Lambert
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_image
use mod_postscript
use stringconstants
use omp_lib
use ISO_C_BINDING
use, intrinsic :: iso_fortran_env

IMPLICIT NONE 

class(ISEmaster_T), INTENT(INOUT)       :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 
type(HDFnames_T), INTENT(INOUT)         :: HDFnames                

type(Cell_T)            :: cell
type(HDF_T)             :: HDF
type(Diffraction_T)     :: Diff
type(IO_T)              :: Message
type(SpaceGroup_T)      :: SG
type(timing_T)          :: timer
type(DynType)           :: Dyn
type(memory_T)          :: mem
type(Lambert_T)         :: L
type(kvectors_T)        :: kvec
type(kvectorlist), pointer :: ktmp
type(poslist),pointer   :: plist, ptmp 

real(kind=dbl)          :: arg, Radius, xyz(3)
integer(HSIZE_T)        :: dims2(2), cnt2(2), offset2(2)
integer(kind=irg)       :: isym,i,j,ik,npy,ipx,ipy,ipz,debug,iE,izz, izzmax, iequiv(3,48), nequiv, num_el, MCnthreads, & ! counters
                           numk, timestart, timestop, numsites, nthreads, k2, icnt, k, mm, mp, & ! number of independent incident beam directions
                           ir,nat(maxpasym),kk(3), skip, ijmax, one, NUMTHREADS, TID, SamplingType, npx, &
                           numset,n,ix,iy,iz, io_int(6), maxlat(1), numcells(3), ncells, atomcnt, ii, TIFF_nx, TIFF_ny, &
                           istat,gzero,ic,ip,ikk, ll, totstrong, totweak, jh, ierr, nix, niy, nixp, niyp     ! counters
real(kind=dbl)          :: tpi,Znsq, kkl, DBWF, kin, delta, h, lambda, omtl, srt, dc(3), xy(2), edge, scl, tmp, dx, dxm, &
                           dy, dym, ca, cb, cc, rsphere, kv(3), g(3), rr(3), ff(3) !
real(kind=sgl)          :: io_real(5), kkk(3), tstop, r2, ISEinten, sxyz(3), ma, mi, dmin
real(kind=sgl),allocatable  :: auxNH(:,:), auxSH(:,:) 
real(kind=dbl),allocatable  :: atomlist(:,:), atomrad(:), ATOM_pos(:,:), ctmp(:,:)  ! results
real(kind=sgl),allocatable  :: mLPNH(:,:), mLPSH(:,:), masterSPNH(:,:), masterSPSH(:,:)
integer(kind=irg), allocatable  :: atp(:)
logical                 :: usehex, switchmirror, verbose
character(fnlen)        :: xtalname

integer(kind=irg)       :: hdferr, nlines
character(fnlen)        :: oldprogname, groupname, energyfile, outname, datagroupname, attributename, HDF_FileVersion, &
                           fname, TIFF_filename 
character(11)           :: dstr
character(15)           :: tstrb
character(15)           :: tstre
logical                 :: f_exists, readonly, overwrite=.TRUE., insert=.TRUE., stereog, g_exists, xtaldataread, FL, doLegendre
character(fnlen,kind=c_char),allocatable,TARGET   :: stringarray(:)
character(fnlen,kind=c_char)                      :: line2(1)

! type(gnode),save                :: rlp
! type(reflisttype),pointer       :: reflist,firstw, rltmp

real(kind=dbl),allocatable      :: karray(:,:)
integer(kind=irg),allocatable   :: kij(:,:)
character(fnlen)                :: dataset, instring

! declare variables for use in object oriented image module
integer                         :: iostat
character(len=128)              :: iomsg
logical                         :: isInteger
type(image_t)                   :: im
integer(int8)                   :: i8 (3,4)
integer(int8), allocatable      :: TIFF_image(:,:)

call openFortranHDFInterface()

associate( nml => self%nml )

stereog = .TRUE.

timer = timing_T()
tstrb = timer%getTimeString()

tpi = 2.D0*cPi

!=============================================
outname = trim(EMsoft%generateFilePath('EMdatapathname',nml%outname))

!=============================================
!=============================================
! crystallography section; 
verbose = .TRUE.
call cell%setFileName(nml%xtalname)
call Diff%setrlpmethod('WK')
dmin = 0.05
call Diff%setV(30.D0)
call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, dmin, verbose)

! check the crystal system and setting; abort the program for trigonal with rhombohedral setting with
! an explanation for the user

if ((SG%getSpaceGroupXtalSystem().eq.5).and.(cell%getLatParm('b').eq.cell%getLatParm('c'))) then
    call Message%printMessage( (/ &
    '                                                                         ', &
    ' ========Program Aborted========                                         ', &
    ' The ISE master pattern simulation for rhombohedral/trigonal structures  ', &
    ' requires that the structure be described using the hexagonal reference  ', &
    ' frame.  Please re-enter the crystal structure in this setting.          '/) )
    stop
end if

! determine the point group number
 j=0
 do i=1,32
  if (SGPG(i).le.SG%getSpaceGroupNumber()) j=i
 end do
 isym = j

! Here, we encode isym into a new number that describes the sampling scheme; the new schemes are 
! described in detail in the EBSD manual pdf file.  We use them here to generate a set of unique
! beam directions which we then convert to rotation matrices to rotate the set of atoms into the 
! correct orientation for the ballistic channeling projection algorithm.

SamplingType = PGSamplingType(isym)

! next, intercept the special cases (hexagonal vs. rhombohedral cases that require special treatment)
if ((SamplingType.eq.-1).or.(isym.eq.14).or.(isym.eq.26)) then 
  SamplingType = SG%getHexvsRho(isym)
end if 

! if the point group is trigonal or hexagonal, we need to switch usehex to .TRUE. so that
! the program will use the hexagonal sampling method
usehex = .FALSE.
if ((SG%getSpaceGroupXtalSystem().eq.4).or.(SG%getSpaceGroupXtalSystem().eq.5)) usehex = .TRUE.

! write (*,*) '========================'
! write (*,*) 'isym = ',isym
! write (*,*) 'SamplingType = ', SamplingType
! write (*,*) 'usehex = ', usehex
! write (*,*) '========================'

! ---------- end of symmetry and crystallography section
!=============================================
!=============================================

!=============================================
!=============================================
npx = nml%npx
npy = npx
gzero = 1  ! index of incident beam
!=============================================
!=============================================

! ---------- allocate memory for the master patterns
mem = memory_T()
call mem%alloc(mLPNH, (/ npx, npy /), 'mLPNH', initval=0.0, startdims = (/ -npx, -npy /))
call mem%alloc(mLPSH, (/ npx, npy /), 'mLPSH', initval=0.0, startdims = (/ -npx, -npy /))
call mem%alloc(masterSPNH, (/ npx, npy /), 'masterSPNH', initval=0.0, startdims = (/ -npx, -npy /))
call mem%alloc(masterSPSH, (/ npx, npy /), 'masterSPSH', initval=0.0, startdims = (/ -npx, -npy /))
! ---------- end allocate memory for the master patterns
!=============================================
!=============================================

!=============================================
! create the HDF5 output file
!=============================================

HDF = HDF_T()

! Create a new file using the default properties.
hdferr =  HDF%createFile(outname)

! write the EMheader to the file
datagroupname = trim(HDFnames%get_ProgramData())
call HDF%writeEMheader(EMsoft,dstr, tstrb, tstre, progname, datagroupname)

! open or create a namelist group to write all the namelist files into
groupname = SC_NMLfiles
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

! create the EBSDmaster group and add a HDF_FileVersion attribbute to it 
  HDF_FileVersion = '4.0'
  HDF_FileVersion = cstringify(HDF_FileVersion)
  attributename = SC_HDFFileVersion
  hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)

! =====================================================
! The following write commands constitute HDF_FileVersion = 4.0
! =====================================================
dataset = SC_xtalname
  allocate(stringarray(1))
  stringarray(1)= trim(nml%xtalname)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1, overwrite)
  else
    hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)
  end if

! create the hyperslabs and write zeroes to them for now
dataset = SC_mLPNH
  dims2 = (/  2*nml%npx+1, 2*nml%npx+1 /)
  cnt2 = (/ 2*nml%npx+1, 2*nml%npx+1 /)
  offset2 = (/ 0, 0 /)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims2, offset2, cnt2, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims2, offset2, cnt2)
  end if

dataset = SC_mLPSH
  dims2 = (/  2*nml%npx+1, 2*nml%npx+1 /)
  cnt2 = (/ 2*nml%npx+1, 2*nml%npx+1 /)
  offset2 = (/ 0, 0 /)
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims2, offset2, cnt2, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims2, offset2, cnt2)
  end if

dataset = SC_masterSPNH
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then 
    hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPNH, dims2, offset2, cnt2, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPNH, dims2, offset2, cnt2)
  end if

dataset = SC_masterSPSH
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then 
    hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPSH, dims2, offset2, cnt2, insert)
  else
    hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPSH, dims2, offset2, cnt2)
  end if

! =====================================================
! end of HDF_FileVersion = 4.0 write statements
! =====================================================

  call HDF%popall()


!=============================================
!=============================================
! ---------- from here on, we need to repeat the entire computation for each energy value
! so this is where we could in principle implement an OpenMP approach; alternatively, 
! we could do the inner loop over the incident beam directions in OpenMP (probably simpler)

call timer%Time_tick(2)

! print a message to indicate where we are in the computation
call Message%printMessage(' -> Initializing incident beam direction list')

!=============================================
! ---------- create the incident beam directions list
! determine all independent incident beam directions (use a linked list starting at khead)
! numk is the total number of k-vectors to be included in this computation;  each beam direction
! will then be converted into the rotation matrix that takes the z-axis into the beam direction.
kvec = kvectors_T()   ! initialize the wave vector list
call kvec%set_kinp( (/ 0.D0, 0.D0, 1.D0 /) )
call kvec%set_ktmax( 0.D0 )
call kvec%set_SamplingType( SamplingType )

call kvec%set_mapmode('RoscaLambert')
if (usehex) then
  call kvec%Calckvectors(cell, SG, Diff, (/ 0.D0, 0.D0, 0.D0 /),npx,npy, ijmax,usehex)
else
  call kvec%Calckvectors(cell, SG, Diff, (/ 0.D0, 0.D0, 0.D0 /),npx,npy, ijmax,usehex)
end if
numk = kvec%get_numk()
io_int(1)=numk
call Message%WriteValue(' # independent beam directions to be considered = ', io_int, 1, "(I8)")

! convert part of the kvector linked list into arrays for OpenMP
call mem%alloc(karray, (/4,numk/), 'karray', 0.D0)
call mem%alloc(kij, (/3,numk/), 'kij', 0)

! point to the first beam direction
ktmp => kvec%get_ListHead()
! and loop through the list, keeping a normalized k and i,j
karray(1:3,1) = ktmp%k(1:3)
kij(1:3,1) = (/ ktmp%i, ktmp%j, ktmp%hs /)
do ik=2,numk
  ktmp => ktmp%next
  karray(1:3,ik) = ktmp%k(1:3)
  kij(1:3,ik) = (/ ktmp%i, ktmp%j, ktmp%hs /)
end do
! and remove the linked list
call kvec%Delete_kvectorlist()

verbose = .FALSE.
totstrong = 0
totweak = 0

! ---------- end of "create the incident beam directions list"
!=============================================

!=============================================
!=============================================
! generate the crystal lattice for the spherical ISE projection 
! we need to generate enough atoms to fill a sphere with radius
! three times the largest lattice parameter, and then keep only
! those atoms (so count them first and then allocate the coordinate array)
ca = cell%getLatParm('a')
cb = cell%getLatParm('b')
cc = cell%getLatParm('c')
mp = nint(self%get_multiplier())
rsphere = self%get_multiplier() * maxval( (/ ca, cb, cc /) )
r2 = rsphere*rsphere
maxlat = maxloc( (/ ca, cb, cc /) )
if (maxlat(1).eq.1) numcells = (/ mp , nint( ca/cb ) * mp , nint( ca/cc ) * mp /) + 1
if (maxlat(1).eq.2) numcells = (/ nint( cb/ca ) * mp, mp, nint( cb/cc ) * mp /) + 1
if (maxlat(1).eq.3) numcells = (/ nint( cc/ca ) * mp, nint( cc/cb ) * mp, mp /) + 1

! we could use the calcPositions routine for this, but instead we make a linked list 
! with all the atom positions inside the sphere... and then, once we know how many atoms
! we have, we convert this to the appropriately sized cell%apos array
allocate(plist)
nullify(plist%next)
ptmp => plist
atomcnt = 0
atp = cell%getAtomtype()
ATOM_pos = cell%getAsymPosData()

do i=1,cell%getNatomtype()
! for each atom in the asymmetric unit
  call SG%calcOrbit(ATOM_pos(i,1:3),n,ctmp)
! replicate in all cells in crystal coordinates
  do j=-numcells(1),numcells(1)
   ff(1)=dble(j)
   do k=-numcells(2),numcells(2)
    ff(2)=dble(k)
    do ll=-numcells(3),numcells(3)
     ff(3)=dble(ll)
     do k2=1,n
      do mm=1,3
       rr(mm)=ctmp(k2,mm)+ff(mm)
      end do
      call cell%TransSpace(rr,g,'d','c')
      if (sum(g*g).le.r2) then   ! it is inside the sphere
        ptmp%xyz = g * 10.0
        ptmp%radius = ATOM_MTradii( atp(i) ) * 10.0 * nml%iscale(3)
        allocate(ptmp%next)
        ptmp => ptmp%next
        nullify(ptmp%next)
        atomcnt = atomcnt+1
      end if 
     end do ! kk
    end do ! l
   end do ! k
  end do ! j
  deallocate(ctmp)
end do 
! that completes the atom list inside the truncation sphere

call mem%alloc(atomlist, (/ 3,atomcnt /), 'atomlist', initval=0.D0) 
call mem%alloc(atomrad, (/ atomcnt /), 'atomrad', initval=0.D0)

ptmp => plist 
icnt = 1
do
  atomlist(1:3,icnt) = ptmp%xyz(1:3)
  atomrad(icnt) = ptmp%radius
  if ((.not.associated(ptmp%next)).or.(icnt.eq.atomcnt)) EXIT
  icnt = icnt + 1
  plist => ptmp%next
  deallocate(ptmp)
  ptmp => plist
end do

rsphere = rsphere * 10.0   ! the ISE algorithm expects Angstrom units ...

io_int(1) = atomcnt
call Message%WriteValue(' Number of atoms generated : ', io_int, 1)

!=============================================
!=============================================

! here's where we introduce the OpenMP calls, to speed up the overall calculations...

! set the number of OpenMP threads 
  if (nml%nthreads.eq.0) then 
    nthreads = OMP_GET_MAX_THREADS()
  else
    nthreads = nml%nthreads
  end if
  call OMP_SET_NUM_THREADS(nthreads)
  io_int(1) = nthreads
  call Message%WriteValue(' Attempting to set number of threads to ',io_int, 1, frm = "(I4)")

! use OpenMP to run on multiple cores ... 
!$OMP PARALLEL PRIVATE(ik,TID,ipx,ipy,ipz,ix,iequiv,nequiv,io_int,ISEinten)

  NUMTHREADS = OMP_GET_NUM_THREADS()
  TID = OMP_GET_THREAD_NUM()

!$OMP DO SCHEDULE(DYNAMIC,100)    
! ---------- and here we start the beam direction loop
   beamloop:do ik = 1,numk

     ISEinten = self%getISEintensity(karray(1:3,ik), atomcnt, atomlist, atomrad, rsphere, &
                                     dble(nml%iscale(1)), dble(nml%iscale(2)),usehex)

! and store the resulting svals values, applying point group symmetry where needed.
     ipx = kij(1,ik)
     ipy = kij(2,ik)
     ipz = kij(3,ik)
!
     if (usehex.eqv..TRUE.) then
       call L%Apply3DPGSymmetry(cell,SG,ipx,ipy,ipz,npx,iequiv,nequiv,usehex)
     else
       if ((SG%getSpaceGroupNumber().ge.195).and.(SG%getSpaceGroupNumber().le.230)) then
         call L%Apply3DPGSymmetry(cell,SG,ipx,ipy,ipz,npx,iequiv,nequiv,cubictype=SamplingType)
       else
         call L%Apply3DPGSymmetry(cell,SG,ipx,ipy,ipz,npx,iequiv,nequiv)
       end if
     end if

!$OMP CRITICAL
     do ix=1,nequiv
       if (iequiv(3,ix).eq.-1) mLPSH(iequiv(1,ix),iequiv(2,ix)) = ISEinten
       if (iequiv(3,ix).eq.1) mLPNH(iequiv(1,ix),iequiv(2,ix)) = ISEinten
     end do
!$OMP END CRITICAL
  
     if (mod(ik,5000).eq.0) then
       io_int(1) = ik
       io_int(2) = numk
       call Message%WriteValue('  completed beam direction ',io_int, 2, "(I8,' of ',I8)")
     end if

    end do beamloop
  
! end of OpenMP portion
!$OMP END PARALLEL

  call mem%dealloc(karray, 'karray')
  call mem%dealloc(kij, 'kij')

! we need to normalize the master pattern arrays to a max value of 1
  r2 = maxval( (/ maxval(mLPSH), maxval(mLPNH) /) )
  mLPNH = mLPNH / r2
  mLPSH = mLPSH / r2

  if (usehex) then
! and finally, we convert the hexagonally sampled array to a square Lambert projection which will be used 
! for all pattern interpolations;  we need to do this for both the Northern and Southern hemispheres

! we begin by allocating auxiliary arrays to hold copies of the hexagonal data; the original arrays will
! then be overwritten with the newly interpolated data.
    call mem%alloc(auxNH, (/npx,npy/), 'auxNH', startdims=(/-npx,-npy/))
    call mem%alloc(auxSH, (/npx,npy/), 'auxSH', startdims=(/-npx,-npy/))
    auxNH = mLPNH
    auxSH = mLPSH

! 
    edge = 1.D0 / dble(nml%npx)
    scl = float(nml%npx) 
    do i=-npx,npx
      do j=-npy,npy
! determine the spherical direction for this point
        L = Lambert_T( xyd = (/ dble(i), dble(j) /) * edge )
        ierr = L%LambertSquareToSphere(dc)
! convert direction cosines to hexagonal Lambert projections
        L = Lambert_T( xyzd = dc )
        ierr = L%LambertSphereToHex(xy)
        xy = xy * scl
! interpolate intensity from the neighboring points
        if (ierr.eq.0) then 
          nix = floor(xy(1))
          niy = floor(xy(2))
          nixp = nix+1
          niyp = niy+1
          if (nixp.gt.npx) nixp = nix
          if (niyp.gt.npx) niyp = niy
          dx = xy(1) - nix
          dy = xy(2) - niy
          dxm = 1.D0 - dx
          dym = 1.D0 - dy
          mLPNH(i,j) = auxNH(nix,niy)*dxm*dym + auxNH(nixp,niy)*dx*dym + auxNH(nix,niyp)*dxm*dy + auxNH(nixp,niyp)*dx*dy
          mLPSH(i,j) = auxSH(nix,niy)*dxm*dym + auxSH(nixp,niy)*dx*dym + auxSH(nix,niyp)*dxm*dy + auxSH(nixp,niyp)*dx*dy
        end if
      end do
    end do
    call mem%dealloc(auxNH, 'auxNH')
    call mem%dealloc(auxSH, 'auxSH')
  end if

! make sure that the outer pixel rim of the mLPSH patterns is identical to
! that of the mLPNH array.
  mLPSH(-npx,-npx:npx) = mLPNH(-npx,-npx:npx)
  mLPSH( npx,-npx:npx) = mLPNH( npx,-npx:npx)
  mLPSH(-npx:npx,-npx) = mLPNH(-npx:npx,-npx)
  mLPSH(-npx:npx, npx) = mLPNH(-npx:npx, npx)

! get stereographic projections (summed over the atomic positions)
  Radius = 1.0
  do i=-npx,npx 
    do j=-npx,npx 
      L = Lambert_T( xyd = (/ dble(i), dble(j) /) / dble(npx) )
      ierr = L%StereoGraphicInverse( xyz, Radius )
      xyz = xyz/vecnorm(xyz)
      if (ierr.ne.0) then 
        masterSPNH(i,j) = 0.0
        masterSPSH(i,j) = 0.0
      else
        sxyz = sngl(xyz)
        masterSPNH(i,j) = InterpolateLambert(sxyz, mLPNH, npx)
        masterSPSH(i,j) = InterpolateLambert(sxyz, mLPSH, npx)
      end if
    end do
  end do

! and here is where the major changes are for version 5.0: all output now in HDF5 format
  call timer%makeTimeStamp()
  dstr = timer%getDateString()
  tstre = timer%getTimeString()

  datagroupname = SC_ISEmaster

 ! open the existing file using the default properties.
  hdferr =  HDF%openFile(outname)

! update the time string
  hdferr = HDF%openGroup(HDFnames%get_EMheader())
  hdferr = HDF%openGroup(datagroupname)

dataset = SC_StopTime
  call timer%Time_tock(2)
  tstop = timer%getInterval(2)
  call timer%Time_reset(2)
  line2(1) = dstr//', '//tstre
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)

  io_int(1) = tstop
  call Message%WriteValue(' Execution time [s]: ',io_int,1)

dataset = SC_Duration
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetFloat(dataset, tstop, overwrite)
  else
    hdferr = HDF%writeDatasetFloat(dataset, tstop)
  end if

  call HDF%pop()
  call HDF%pop()
  
  hdferr = HDF%openGroup(HDFnames%get_EMData())
  hdferr = HDF%openGroup(datagroupname)

! add data to the hyperslab
dataset = SC_mLPNH
  dims2 = (/  2*npx+1, 2*npx+1 /)
  cnt2 = (/ 2*npx+1, 2*npx+1 /)
  offset2 = (/ 0, 0 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, mLPNH, dims2, offset2, cnt2, insert)

dataset = SC_mLPSH
  hdferr = HDF%writeHyperslabFloatArray(dataset, mLPSH, dims2, offset2, cnt2, insert)

dataset = SC_masterSPNH
  hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPNH, dims2, offset2, cnt2, insert)

dataset = SC_masterSPSH
  hdferr = HDF%writeHyperslabFloatArray(dataset, masterSPSH, dims2, offset2, cnt2, insert)

  call HDF%popall()

! and close the fortran hdf interface
 call closeFortranHDFInterface()

  call Message%printMessage('Final data stored in file '//trim(nml%outname), frm = "(A/)")

! output the stereographic projection as a tiff file 
fname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(nml%tiffname)//'_SPNH.tiff'
TIFF_filename = trim(fname)
TIFF_nx = 2*npx+1 
TIFF_ny = 2*npx+1 

! allocate memory for image
allocate(TIFF_image( TIFF_nx,TIFF_ny ))

! fill the image with whatever data you have (between 0 and 255)
ma = maxval(masterSPNH)
mi = minval(masterSPNH)

TIFF_image = int(255 * (masterSPNH-mi)/(ma-mi))

! set up the image_t structure
im = image_t(TIFF_image)
if(im%empty()) call Message%printMessage("ComputeMasterPattern","failed to convert array to image")

! create the file
call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage("failed to write image to file : "//iomsg)
else  
  call Message%printMessage('ISE masterSPNH map written to '//trim(TIFF_filename))
end if 
deallocate(TIFF_image)

call mem%dealloc(mLPNH, 'mLPNH')
call mem%dealloc(mLPSH, 'mLPSH')
call mem%dealloc(masterSPNH, 'masterSPNH')
call mem%dealloc(masterSPSH, 'masterSPSH')

end associate

end subroutine ISEmaster_



end module mod_ISEmaster