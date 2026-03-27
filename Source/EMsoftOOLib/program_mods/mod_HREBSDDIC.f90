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

module mod_HREBSDDIC
  !! author: MDG (OO version)
  !! version: 1.0 
  !! date: 11/07/24
  !!
  !! class definition for the EMHREBSDDIC program
  !! uses the inverse compositional Gauss–Newton algorithm from ATEX
  !!
  !! this program does all the pattern pre-processing; this was originally
  !! a separate program EMppEBSD, but there is really no need to do that 
  !! and store another possibly very large file...
  !! Instead, we have the user fill out a normal namelist file with a few 
  !! parameters that are also part of a DI file; there should not be any
  !! need to carry out an indexing run and refinement before running the 
  !! DIC algorithm (although in a future version we will need to have the 
  !! orientations so that we can use simulated reference patterns...)

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMHREBSDDIC program
type, public :: HREBSDDICNameListType
   integer(kind=irg)    :: patx
   integer(kind=irg)    :: paty
   integer(kind=irg)    :: nbx
   integer(kind=irg)    :: nby
   integer(kind=irg)    :: maxnumit
   integer(kind=irg)    :: ipf_ht
   integer(kind=irg)    :: ipf_wd
   integer(kind=irg)    :: maskradius
   integer(kind=irg)    :: exptnumsx
   integer(kind=irg)    :: exptnumsy
   integer(kind=irg)    :: binning
   integer(kind=irg)    :: nthreads
   integer(kind=irg)    :: nregions
   integer(kind=irg)    :: logparam
   real(kind=dbl)       :: hipassw
   real(kind=sgl)       :: L
   real(kind=sgl)       :: thetac
   real(kind=sgl)       :: sigma
   real(kind=sgl)       :: delta
   real(kind=sgl)       :: omega
   real(kind=sgl)       :: xpc
   real(kind=sgl)       :: ypc
   real(kind=sgl)       :: energymin        ! not used for now
   real(kind=sgl)       :: energymax        ! not used for now
   real(kind=sgl)       :: stepX
   real(kind=sgl)       :: stepY
   real(kind=sgl)       :: C11
   real(kind=sgl)       :: C12
   real(kind=sgl)       :: C44
   real(kind=sgl)       :: C13
   real(kind=sgl)       :: C33
   real(kind=dbl)       :: mindeltap
   real(kind=dbl)       :: scalingfactor
   character(fnlen)     :: pixelornormalized
   character(fnlen)     :: exptfile
   character(fnlen)     :: datafile
   character(fnlen)     :: tmpfile
   character(fnlen)     :: inputtype
   character(3)         :: crystal
   character(3)         :: filtertype
   character(1)         :: keeptmpfile
   character(1)         :: maskpattern
   logical              :: verbose
   integer(kind=irg)    :: cross(4)
   integer(kind=irg)    :: ROI(4)
   integer(kind=irg)    :: refpatpos(2)
   character(fnlen)     :: HDFstrings(10)
end type HREBSDDICNameListType

! class definition
type, public :: HREBSD_DIC_T
private 
  character(fnlen)                        :: nmldeffile = 'EMHREBSDDIC.nml'
  type(HREBSDDICNameListType)             :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: HREBSD_DIC_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: HREBSD_DIC => HREBSD_DIC_

end type HREBSD_DIC_T

! the constructor routine for this class 
interface HREBSD_DIC_T
  module procedure HREBSDDIC_constructor
end interface HREBSD_DIC_T

contains

!--------------------------------------------------------------------------
type(HREBSD_DIC_T) function HREBSDDIC_constructor( nmlfile ) result(HREBSDDIC)
!DEC$ ATTRIBUTES DLLEXPORT :: HREBSDDIC_constructor
!! author: MDG 
!! version: 1.0 
!! date: 11/07/24
!!
!! constructor for the HREBSD_DIC_T Class; 

IMPLICIT NONE

character(fnlen), INTENT(IN)  :: nmlfile 

call HREBSDDIC%readNameList(nmlfile)

end function HREBSDDIC_constructor

!--------------------------------------------------------------------------
subroutine HREBSDDIC_destructor(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: HREBSDDIC_destructor
!! author: MDG 
!! version: 1.0 
!! date: 11/07/24
!!
!! destructor for the HREBSD_DIC_T Class
 
IMPLICIT NONE

type(HREBSD_DIC_T), INTENT(INOUT)  :: self 

call reportDestructor('HREBSD_DIC_T')

end subroutine HREBSDDIC_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 11/07/24
!!
!! read the namelist from an nml file for the HREBSD_DIC_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(HREBSD_DIC_T), INTENT(INOUT) :: self
character(fnlen),INTENT(IN)        :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)        :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                     :: EMsoft 
type(IO_T)                         :: Message       
logical                            :: skipread = .FALSE.

integer(kind=irg)                  :: patx
integer(kind=irg)                  :: paty
integer(kind=irg)                  :: nbx
integer(kind=irg)                  :: nby
integer(kind=irg)                  :: maxnumit
integer(kind=irg)                  :: nthreads
integer(kind=irg)                  :: ipf_ht
integer(kind=irg)                  :: ipf_wd
integer(kind=irg)                  :: maskradius
integer(kind=irg)                  :: exptnumsx
integer(kind=irg)                  :: exptnumsy
integer(kind=irg)                  :: binning
integer(kind=irg)                  :: nregions
integer(kind=irg)                  :: logparam
real(kind=dbl)                     :: hipassw
real(kind=sgl)                     :: L
real(kind=sgl)                     :: thetac
real(kind=sgl)                     :: sigma
real(kind=sgl)                     :: delta
real(kind=sgl)                     :: omega
real(kind=sgl)                     :: xpc
real(kind=sgl)                     :: ypc
real(kind=sgl)                     :: energymin        ! not used for now
real(kind=sgl)                     :: energymax        ! not used for now
real(kind=sgl)                     :: stepX
real(kind=sgl)                     :: stepY
real(kind=sgl)                     :: C11
real(kind=sgl)                     :: C12
real(kind=sgl)                     :: C44
real(kind=sgl)                     :: C13
real(kind=sgl)                     :: C33
real(kind=dbl)                     :: mindeltap
real(kind=dbl)                     :: scalingfactor
character(fnlen)                   :: pixelornormalized
character(fnlen)                   :: exptfile
character(fnlen)                   :: datafile
character(fnlen)                   :: tmpfile
character(fnlen)                   :: inputtype
character(3)                       :: crystal
character(3)                       :: filtertype
character(1)                       :: keeptmpfile
character(1)                       :: maskpattern
logical                            :: verbose
integer(kind=irg)                  :: cross(4)
integer(kind=irg)                  :: ROI(4)
integer(kind=irg)                  :: refpatpos(2)
character(fnlen)                   :: HDFstrings(10)


namelist / HREBSDDICdata / patx, paty, nthreads, C11, C12, C44, C13, C33, mindeltap, verbose, sigma, ROI, binning, & 
                           datafile, exptfile, crystal, nbx, nby, maxnumit, scalingfactor, maskpattern, &
                           pixelornormalized, cross, refpatpos, ipf_ht, ipf_wd, maskradius, exptnumsx, &
                           exptnumsy, nregions, logparam, hipassw, L, thetac, delta, omega, xpc, ypc, filtertype, &
                           energymin, energymax, StepX, StepY, tmpfile, inputtype, keeptmpfile, HDFstrings

patx = 0
paty = 0
nbx = 0
nby = 0
maxnumit = 50
nthreads = 1
ipf_ht = 100
ipf_wd = 100
maskradius = 240
exptnumsx = 640
exptnumsy = 480
binning = 1
nregions = 10
logparam = 10
hipassw = 0.05
L = 15000.0
thetac = 10.0
sigma = 70.0
delta = 50.0
omega = 0.0
xpc = 0.0
ypc = 0.0
energymin = 10.0      ! not used for now = 
energymax = 20.0      ! not used for now = 
stepX = 1.0
stepY = 1.0
C11 = 276.0
C12 = 159.0
C44 = 132.0
C13 = 0.0
C33 = 0.0
mindeltap = 0.001D0
scalingfactor = 1.5D0
pixelornormalized = 'normalized'
exptfile = 'undefined'
datafile = 'undefined'
tmpfile = 'undefined'
inputtype = 'Binary'
crystal = 'cub'
filtertype = 'fft'
keeptmpfile = 'n'
maskpattern = 'n'
verbose = .FALSE.
cross = (/ 0, 0, 0, 0 /)
ROI = (/ 0, 0, 0, 0 /)
refpatpos = (/ -1, -1 /)
HDFstrings = (/ '', '', '', '', '', '', '', '', '', '' /)

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=HREBSDDICdata)
    close(UNIT=dataunit,STATUS='keep')

! check for required entries
    if (minval(refpatpos).lt.0) then 
      call Message%printError('readNameList:',' refpatpos must be positive coordinates in '//nmlfile)
    end if

    if (trim(exptfile).eq.'undefined') then
      call Message%printError('readNameList:',' exptfile name is undefined in '//nmlfile)
    end if

    if (trim(datafile).eq.'undefined') then
      call Message%printError('readNameList:',' datafile file name is undefined in '//nmlfile)
    end if

end if

self%nml%patx = patx
self%nml%paty = paty
self%nml%nbx = nbx
self%nml%nby = nby
self%nml%maxnumit = maxnumit
self%nml%nthreads = nthreads
self%nml%ipf_ht = ipf_ht
self%nml%ipf_wd = ipf_wd
self%nml%maskradius = maskradius
self%nml%exptnumsx = exptnumsx
self%nml%exptnumsy = exptnumsy
self%nml%binning = binning
self%nml%nregions = nregions
self%nml%logparam = logparam
self%nml%hipassw = hipassw
self%nml%L = L
self%nml%thetac = thetac
self%nml%sigma = sigma
self%nml%delta = delta
self%nml%omega = omega
self%nml%xpc = xpc
self%nml%ypc = ypc
self%nml%energymin  = energymin
self%nml%energymax  = energymax
self%nml%stepX = stepX
self%nml%stepY = stepY
self%nml%C11 = C11
self%nml%C12 = C12
self%nml%C44 = C44
self%nml%C13 = C13
self%nml%C33 = C33
self%nml%mindeltap = mindeltap
self%nml%scalingfactor = scalingfactor
self%nml%pixelornormalized = pixelornormalized
self%nml%exptfile = exptfile
self%nml%datafile = datafile
self%nml%tmpfile = tmpfile
self%nml%inputtype = inputtype
self%nml%crystal = crystal
self%nml%filtertype = filtertype
self%nml%keeptmpfile = keeptmpfile
self%nml%maskpattern = maskpattern
self%nml%verbose = verbose
self%nml%cross = cross
self%nml%ROI = ROI
self%nml%refpatpos = refpatpos
self%nml%HDFstrings = HDFstrings

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 11/07/24
!!
!! pass the namelist for the HREBSD_DIC_T Class to the calling program

IMPLICIT NONE 

class(HREBSD_DIC_T), INTENT(INOUT)           :: self
type(HREBSDDICNameListType)                  :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 11/07/24
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(HREBSD_DIC_T), INTENT(INOUT)      :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 15, n_real = 19
integer(kind=irg)                       :: hdferr,  io_int(n_int), vb
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1),line10(10)

associate( enl => self%nml )

vb = 0
if (self%nml%verbose.eqv..TRUE.) vb = 1 

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! write all the single integers
io_int = (/ enl%nthreads, enl%patx, enl%paty, enl%nbx, enl%nby, enl%maxnumit, vb, enl%ipf_ht, enl%ipf_wd, &
            enl%maskradius, enl%exptnumsx, enl%exptnumsy, enl%nregions, enl%logparam, enl%binning /)
intlist(1) = 'nthreads'
intlist(2) = 'patx'
intlist(3) = 'paty'
intlist(4) = 'nbx'
intlist(5) = 'nby'
intlist(6) = 'maxnumit'
intlist(7) = 'verbose'
intlist(8) = 'ipf_ht'
intlist(9) = 'ipf_wd'
intlist(10) = 'maskradius'
intlist(11) = 'exptnumsx'
intlist(12) = 'exptnumsy'
intlist(13) = 'nregions'
intlist(14) = 'logparam'
intlist(15) = 'binning'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ enl%C11, enl%C12, enl%C44, enl%C13, enl%C33, real(enl%mindeltap),real(enl%scalingfactor), &
             real(enl%hipassw), enl%L, enl%thetac, enl%delta, enl%omega, enl%xpc, enl%ypc, enl%energymin, enl%energymax, &
             enl%stepX, enl%stepY, enl%sigma /)
reallist(1) = 'C11'
reallist(2) = 'C12'
reallist(3) = 'C44'
reallist(4) = 'C13'
reallist(5) = 'C33'
reallist(6) = 'mindeltap'
reallist(7) = 'scalingfactor'
reallist(8) = 'hipassw'
reallist(9) = 'L'
reallist(10) = 'thetac'
reallist(11) = 'delta'
reallist(12) = 'omega'
reallist(13) = 'xpc'
reallist(14) = 'ypc'
reallist(15) = 'energymin'
reallist(16) = 'energymax'
reallist(17) = 'stepX'
reallist(18) = 'stepY'
reallist(19) = 'sigma'
call HDF%writeNMLreals(io_real, reallist, n_real)

dataset = 'cross'
hdferr = HDF%writeDatasetIntegerArray(dataset, enl%cross, 4)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create cross dataset', hdferr)

dataset = 'ROI'
hdferr = HDF%writeDatasetIntegerArray(dataset, enl%ROI, 4)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create ROI dataset', hdferr)

dataset = 'refpatpos'
hdferr = HDF%writeDatasetIntegerArray(dataset, enl%refpatpos, 2)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create cross dataset', hdferr)

! write all the strings
dataset = SC_datafile
line2(1) = trim(enl%datafile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create datafile dataset', hdferr)

dataset = 'exptfile' 
line2(1) = trim(enl%exptfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create exptfile dataset', hdferr)

dataset = 'tmpfile'
line2(1) = trim(enl%tmpfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create tmpfile dataset', hdferr)

dataset = 'inputtype'
line2(1) = trim(enl%inputtype)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create inputtype dataset', hdferr)

dataset = 'crystal'
line2(1) = trim(enl%crystal)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create crystal dataset', hdferr)

dataset = 'pixelornormalized'
line2(1) = trim(enl%pixelornormalized)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create pixelornormalized dataset', hdferr)

dataset = 'filtertype'
line2(1) = trim(enl%filtertype)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create filtertype dataset', hdferr)

dataset = 'keeptmpfile'
line2(1) = trim(enl%keeptmpfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create keeptmpfile dataset', hdferr)

dataset = 'maskpattern'
line2(1) = trim(enl%maskpattern)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create maskpattern dataset', hdferr)

dataset = SC_HDFstrings
line10 = enl%HDFstrings
hdferr = HDF%writeDatasetStringArray(dataset, line10, 10)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create HDFstrings dataset', hdferr)

! pop this group off the stack
call HDF%pop()
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine HREBSD_DIC_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: HREBSD_DIC_
!! author: MDG 
!! version: 1.0 
!! date: 11/21/24
!!
!! perform the computations

use mod_EMsoft
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_OMPsupport
use mod_quaternions
use mod_rotations 
use mod_io
use mod_HDFsupport
use mod_math
use HDF5
use mod_ppEBSD
use mod_DIC
use mod_DIfiles
use mod_patterns
use mod_timing
use mod_memory
use mod_fftw3
use stringconstants
use mod_platformsupport
use bspline_module
use bspline_kinds_module, only: wp, ip
use ISO_C_BINDING

IMPLICIT NONE 

class(HREBSD_DIC_T), INTENT(INOUT)      :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

type(DIfile_T)                          :: DIFT
type(QuaternionArray_T)                 :: qAR
type(IO_T)                              :: Message
type(Timing_T)                          :: timer
type(memory_T)                          :: mem, memth
type(HDF_T)                             :: HDF
type(e_T)                               :: eu
type(o_T)                               :: om
type(q_T)                               :: qu
type(a_T)                               :: ax
type(Quaternion_T)                      :: quat
type(DIC_T)                             :: DIC
type(HDFnames_T)                        :: HDFnames2
type(EBSDDINameListType)                :: dinl

character(fnlen)                        :: inpfile, HDFstring
real(kind=dbl)                          :: stepsizes(3), sig, totaltilt, ave, cosang, sinang, patcentx, patcenty, &
                                           stepx, stepy, Dref, deltaD, alphainv, gamma1, gamma2
real(kind=sgl)                          :: io_real(6), ss, tstop
real(kind=dbl),allocatable              :: delta(:,:), homographies(:,:), normdp(:), residuals(:), Fehat(:,:,:), &
                                           correctedhomographies(:,:), deltaDD(:), alpha(:)
real(kind=sgl),allocatable              :: exptpattern(:,:), targetpattern(:,:)
integer(kind=irg),allocatable           :: nit(:)
integer(kind=irg)                       :: hdferr, binx, biny, L, recordsize, patsz, ROI_size, sz(2), i, j, kk, numangles, &
                                           istat, interp_grid, interp_size, info, offset, nbx, nby, ii, jj, ierr, io_int(2), &
                                           numpats, TID, cnt, correctsize, recordsize_correct, totnumexpt

real(kind=dbl)                          :: interp_step, std, sf
real(kind=dbl)                          :: Wnew(3,3), Winv(3,3), oldW(3,3)
real(wp)                                :: hg(8), W(3,3), ndp, oldnorm, SOL(8), Hessian(8,8), CIC, PCx, PCy, hpartial(8)
character(11)                           :: dstr
character(15)                           :: tstrb
character(15)                           :: tstre
! integer(HSIZE_T)                        :: dims2(2), dims3(3), offset3(3)
real(kind=dbl)                          :: C(6,6)
character(fnlen)                        :: DIfile, groupname, dataset, datagroupname, outname, fname, gname, hostname
logical                                 :: f_exists, g_exists, overwrite = .TRUE.
integer(kind=4)                         :: hnstat
real(kind=sgl),allocatable              :: mask(:,:),masklin(:)

! is this program executed from home or work ?  [code will be removed in Release version]
hnstat = system_hostnm(hostname) 
write (*,*) trim(hostname)

! this program reads a dot product EMDI file, including all the experimental parameters
! as well as a pre-processed pattern file generated by EMppEBSD.

associate( enl => self%nml )

call setRotationPrecision('d')
call openFortranHDFInterface()

! use the memory class for most array allocations
mem = Memory_T()

! make sure that the exptfile actually exists 
f_exists = .FALSE.
fname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(enl%exptfile)
inquire(file=trim(fname), exist=f_exists)

if (.not.f_exists) then 
  call Message%printError('HREBSD_DIC_',' pattern file not found ')
end if

! define some pattern dimension parameters
binx = enl%exptnumsx
biny = enl%exptnumsy
L = binx * biny
recordsize = 4 * L
patsz = L
if (sum(enl%ROI).eq.0) then 
  totnumexpt = enl%ipf_wd*enl%ipf_ht
else
  totnumexpt = enl%ROI(3)*enl%ROI(4)
end if

! make sure that correctsize is a multiple of 16; if not, make it so
! (probably not necessary in this program but kep in place for consistency)
if (mod(L,16) .ne. 0) then
    correctsize = 16*ceiling(float(L)/16.0)
else
    correctsize = L
end if

recordsize_correct = correctsize*4

if (trim(enl%pixelornormalized).eq.'normalized') then
  PCx = 0.5_wp
  PCy = 0.5_wp
else
  PCx = real(binx,wp)/2.0_wp
  PCy = real(biny,wp)/2.0_wp
end if

! initialize the timing routines
timer = Timing_T()
tstrb = timer%getTimeString()

!=====================================================
! Preprocess all the experimental patterns and store
! them in a temporary file as vectors
!=====================================================
! first, make sure that this file does not already exist
g_exists = .FALSE.
gname = trim(EMsoft%generateFilePath('EMtmppathname'))//trim(enl%tmpfile)
inquire(file=trim(fname), exist=f_exists)

! if the file exists, delete and recreate it
if (g_exists) then
    open(unit=dataunit,file=trim(gname),&
      status='unknown',form='unformatted',access='direct',recl=recordsize_correct,iostat=ierr)
    close(unit=dataunit,status='delete')
end if

call mem%alloc(mask, (/ binx,biny /), 'mask', 1.0)
call mem%alloc(masklin, (/ L /), 'masklin', 0.0)
if (enl%maskpattern.eq.'y') then
  do ii = 1,biny
      do jj = 1,binx
          if((ii-biny/2)**2 + (jj-binx/2)**2 .ge. enl%maskradius**2) then
              mask(jj,ii) = 0.0
          end if
      end do
  end do
end if

! convert the mask to a linear (1D) array
do ii = 1,biny
    do jj = 1,binx
        masklin((ii-1)*binx+jj) = mask(jj,ii)
    end do
end do

! copy the relevant enl parameters into the dinl structure
dinl%nthreads = enl%nthreads
dinl%hipassw = enl%hipassw
dinl%ROI = enl%ROI
dinl%ipf_wd = enl%ipf_wd
dinl%ipf_ht = enl%ipf_ht
dinl%tmpfile = trim(gname)
dinl%exptfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(enl%exptfile)
dinl%inputtype = enl%inputtype
dinl%HDFstrings = enl%HDFstrings
dinl%nregions = enl%nregions
dinl%DIModality = 'EBSD'

HDF = HDF_T() 

if (enl%filtertype.eq.'fft') then 
! standard hipass and adaptive histogram equalization pre-processing step
  call PreProcessPatterns(EMsoft, HDF, .FALSE., dinl, binx, biny, masklin, correctsize, totnumexpt)
else
  ! logarithmic high-pass filter only
  call PreProcessPatterns(EMsoft, HDF, .FALSE., dinl, binx, biny, masklin, correctsize, totnumexpt, &
                          log=.TRUE., logparam=enl%logparam)
end if 
!=====================================================
! end of pattern preprocessing
!=====================================================


!=====================================================
! get the reference pattern 
!=====================================================
call Message%printMessage(' loading reference pattern')
! file exists to the first thing we do is extract the reference pattern; 
call mem%alloc(exptpattern, (/ binx,biny /), 'exptpattern', 0.0 )
open(unit=dataunit,file=trim(gname),&
     status='unknown',form='unformatted',access='direct',recl=recordsize,iostat=ierr)
if (sum(enl%ROI).ne.0) then 
! these are now offsets into the reduced pre-processed file !!!
  offset = enl%paty * enl%ROI(3) + enl%patx + 1
  numpats = enl%ROI(3) * enl%ROI(4)
else
  offset = enl%paty * enl%ipf_wd + enl%patx + 1
  numpats = enl%ipf_wd * enl%ipf_ht
end if
read(dataunit,rec=offset) exptpattern
close(dataunit,status='keep')

exptpattern = exptpattern - minval(exptpattern)
exptpattern = exptpattern/maxval(exptpattern)

! allocate global arrays to store the output of the DIC routine.
call mem%alloc(normdp, (/ numpats /), 'normdp', initval=0.D0 ) 
call mem%alloc(homographies, (/ 8, numpats /), 'homographies', initval=0.D0 ) 
call mem%alloc(residuals, (/ numpats /), 'residuals', initval=0.D0 ) 
call mem%alloc(nit, (/ numpats /), 'nit', initval=0 ) 

! open the pattern file so that each thread can read from it
open(unit=dataunit,file=trim(gname),&
     status='unknown',form='unformatted',access='direct',recl=recordsize,iostat=ierr)

call OMP_setNThreads(enl%nthreads)

memth = memory_T( nt = enl%nthreads )

sf = enl%scalingfactor

call Message%printMessage(' starting DIC analysis')

! here we go parallel with OpenMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID, DIC, nbx, nby, targetpattern, ii, jj, hg, W)&
!$OMP& PRIVATE(oldnorm, oldW, SOL, ndp, Wnew, Winv, io_int, i, hpartial) 

! set the thread ID
TID = OMP_GET_THREAD_NUM()

!copy the exptpattern into the thread's DIC class
if (sum(enl%cross).gt.0) then 
  DIC = DIC_T( binx, biny, normalize = .TRUE., cross=enl%cross )
else
  DIC = DIC_T( binx, biny, normalize = .TRUE. )
end if
call DIC%setverbose( enl%verbose )
call DIC%setpattern( 'r', dble(exptpattern) )

! define the border widths nbx and nby for the subregion;
! if requested, a cross-shaped area (vertical & horizontal) will
! be eliminated from the subregion (the DIC algorithm does not depend
! on the actual shape of the subregion...)
nbx = enl%nbx
nby = enl%nby
call DIC%defineSR( nbx, nby, PCx, PCy) 

! generate the b-splines for the reference pattern and verify the accuracy
! also get the derivatives of the reference pattern to determine the Hessian
call DIC%getbsplines(refp=.TRUE., verify=.TRUE., grads=.TRUE.)

! zero-mean and normalize the referenceSR array
call DIC%applyZMN(doreference=.TRUE.)

! compute the Hessian matrix
call DIC%getHessian(Hessian)

! we'll do a 'fake' run with a null homography to make sure all parameters and 
! auxiliary arrays are properly defined.
hpartial = (/ (0.0_wp, i=1,8) /)
call DIC%applyHomography(hpartial, PCx, PCy)

! allocate the targetpattern array
call memth%alloc(targetpattern, (/ binx,biny /), 'targetpattern', TID=TID )

! in the main loop we iterate over all the experimental patterns in the input file
!$OMP DO SCHEDULE(DYNAMIC)
do ii=1, numpats

!$OMP CRITICAL
  read(dataunit,rec=ii) targetpattern
!$OMP END CRITICAL

  targetpattern = targetpattern - minval(targetpattern)
  targetpattern = targetpattern/maxval(targetpattern)

  call DIC%setpattern( 'd', dble(targetpattern) )
  call DIC%getbsplines(defp=.TRUE.)
! ! zero-mean and normalize the referenceSR and targetSR arrays
  call DIC%applyZMN(dotarget=.TRUE.)
  call DIC%getresiduals()

! loop until max iterations or convergence
  do jj = 1, enl%maxnumit   ! initialize with the null homography
    if (jj.eq.1) then 
      hpartial = (/ (0.0_wp, i=1,8) /)
      W = DIC%getShapeFunction(hpartial)
    end if
    call DIC%applyHomography(hpartial, PCx, PCy, dotarget=.TRUE.)
    call DIC%applyZMN(dotarget=.TRUE.)
    call DIC%getresiduals(CIC)
    call DIC%solveHessian(SOL, ndp)

! convert to a shape function and take the inverse
! we use a multiplicative factor to increase the value of the increment; testing 
! shows that in ideal conditions (simulated patterns), a value of sf=1.5 can
! speed up convergence by a factor of 2
    Wnew = DIC%getShapeFunction(reshape(SOL, (/8/))*sf)
    Winv = matrixInvert_wp( Wnew )
    W = matmul( W, Winv )
    W = W / W(3,3)
    hpartial = reshape(SOL, (/8/)) * sf
    if (ndp.lt.enl%mindeltap) then
        EXIT
    end if
  end do

  if (mod(ii,250).eq.0) then 
    io_int(1) = ii
    io_int(2) = numpats
    call Message%WriteValue(' completed # patterns/total ',io_int,2)
  end if  

! store the results
  hg = DIC%getHomography(W)
  if (jj.eq.enl%maxnumit+1) then  ! zero solution if no convergence is reached
    homographies(1:8,ii) = (/ (0.0_wp, i=1,8) /)
  else
    homographies(1:8,ii) = dble(hg)
  end if
  normdp(ii) = dble(ndp)
  residuals(ii) = CIC
  nit(ii) = jj

  ! do some cleanup before the next one...
  call DIC%cleanup()
end do
!$OMP END DO
call memth%dealloc(targetpattern, 'targetpattern', TID=TID)
!$OMP BARRIER
!$OMP END PARALLEL

if (enl%keeptmpfile.eq.'n') then 
  close(dataunit,status='delete')
else
  close(dataunit,status='keep')
end if

call Message%printMessage(' DIC analysis complete')

!=========================================
! geometry correction of the homographies
!=========================================
call Message%printMessage(' starting detector geometry correction')

! homographies(2:4,:) = -homographies(2:4,:)
! homographies(6,:) = -homographies(6,:)
! homographies(8,:) = -homographies(8,:)

homographies = -homographies 

! first we need the pattern geometry for each pattern 
! since we're using normalized coordinates, these are 
! all in units of the pattern dimensions
! this then allows us to correct for the sampling location
! and remove isotropic scaling due  to changes in DD
! as well as changes in pattern center (which are uniform translations)
allocate( delta(2,numpats), deltaDD(numpats), correctedhomographies(8, numpats), alpha(numpats) )

! we need to make sure we get the units correctly converted; the EMsoft (pcx,pcy)
! values are in units of pixels, and the y-direction is opposite to the one used
! in the DIC implementation; converting these to normalized coordinates [0,1]
! we have the following relations:
patcentx = 0.5D0 + enl%xpc ! / dble(binx)
patcenty = 0.5D0 - enl%ypc ! / dble(biny)

! step sizes in normalized units; they are in microns in the DI file
stepx = enl%StepX / dble( enl%delta ) ! * binx )
stepy = enl%StepY / dble( enl%delta ) ! * biny )

! L and step size in normalized coordinates [ angle sigma needs to be imported !!! ]
Dref = enl%L / dble( enl%delta ) ! * binx )
deltaD = -stepy * sin( (enl%thetac+90.D0-70.D0) * dtor )

! DD and PC arrays
cnt = 1
if (sum(enl%ROI).eq.0) then
  do j=0,enl%ipf_ht-1
    do i=0,enl%ipf_wd-1
      deltaDD(cnt) =   (dble(j)-enl%refpatpos(2)) * deltaD !  + Dref for the full value
      delta(1:2,cnt) = (/ (dble(i)-enl%refpatpos(1)) * stepx, & ! + patcentx
                          (dble(j)-enl%refpatpos(2)) * stepy /) ! + patcenty
      alpha(cnt) = (Dref-deltaDD(cnt)) / Dref
      cnt = cnt+1
    end do 
  end do 
else
  do j=0,enl%ROI(4)-1
    do i=0,enl%ROI(3)-1
      deltaDD(cnt) =   (dble(j)-enl%refpatpos(2)) * deltaD !  + Dref for the full value
      delta(1:2,cnt) = (/ (dble(i)-enl%refpatpos(1)) * stepx, & ! + patcentx
                          (dble(j)-enl%refpatpos(2)) * stepy /) ! + patcenty
      alpha(cnt) = (Dref-deltaDD(cnt)) / Dref
      cnt = cnt+1
    end do 
  end do 
end if

! next we do the actual homography corrections ... 
! in terms of the parameters used in Ernould Chapter 2, section 3.3.2, PC[1:2,*] is equivalent 
! to delta_1 and delta_2 for all sampling points.  x_{0i} are the negative of the normalized
! coordinates of the pattern center 
do ii=1,numpats
  alphainv = 1.0D0 / alpha(ii)
  ! gamma1 = (delta(1,ii) - patcentx*(alpha(ii) - 1.D0))
  ! gamma2 = (delta(2,ii) - patcenty*(alpha(ii) - 1.D0))
  gamma1 = (delta(1,ii) + (delta(1,ii)-patcentx)*(alpha(ii) - 1.D0))
  gamma2 = (delta(2,ii) + (delta(2,ii)-patcenty)*(alpha(ii) - 1.D0))
  hg = homographies(1:8,ii) ! we actually find the opposite homography due to the interpolation routines
  correctedhomographies(1:8,ii) = (/ alphainv*(hg(1)+1.D0-gamma1*hg(7))-1.D0, &
                                     alphainv*(hg(2)-gamma1*hg(8)), &
                                     alphainv*(hg(3)-gamma1), &
                                     alphainv*(hg(4)-gamma2*hg(7)), &
                                     alphainv*(hg(5)+1.D0-gamma2*hg(8))-1.D0, &
                                     alphainv*(hg(6)-gamma2), &
                                     hg(7), &
                                     hg(8) /)
end do

call Message%printMessage(' computing Fehat deformation tensor')

! next we get the reduced deformation tensor in the detector reference frame
allocate( Fehat(3, 3, numpats) )
DIC = DIC_T( binx, biny, normalize = .TRUE. )
do ii=1, numpats
  ! Fehat(1:3,1:3,ii) = DIC%homography2Fe(correctedhomographies(1:8,ii), &
  !                                       (/patcentx,patcenty/)+delta(1:2,ii), Dref + deltaDD(ii))
  hg = homographies(1:8,ii)
  Fehat(1:3,1:3,ii) = DIC%homography2Fe(hg, &
                                        (/patcentx,patcenty/)+delta(1:2,ii), Dref + deltaDD(ii))
end do 



! ! stiffness matrix in the crystal frame
! C = 0.D0  
! C(1,1)=dble(enl%C11)
! C(2,2)=dble(enl%C11)
! C(1,2)=dble(enl%C12)
! C(2,1)=dble(enl%C12)
! C(4,4)=dble(enl%C44)
! C(5,5)=dble(enl%C44)

! ! add equivalent entries for cubic crystal
! if (enl%crystal.eq.'cub') then
!   C(3,3)=dble(enl%C11)    
!   C(1,3)=dble(enl%C12) 
!   C(2,3)=dble(enl%C12)
!   C(3,1)=dble(enl%C12)
!   C(3,2)=dble(enl%C12) 
!   C(6,6)=dble(enl%C44)
! ! or for hexagonal close packed crystal
! else if (enl%crystal.eq.'hex') then
!   C(3,3)=dble(enl%C33)    
!   C(1,3)=dble(enl%C13) 
!   C(2,3)=dble(enl%C13)
!   C(3,1)=dble(enl%C13)
!   C(3,2)=dble(enl%C13)
!   C(6,6)=dble(enl%C11-enl%C12)/2.D0  
! else
!   call Message%printError('HREBSD_','Undefined crystal type '//trim(enl%crystal))
! end if

! call Message%printMessage(' Stiffness tensor (unit:GPa) = ')

! do i = 1, 6
!   io_real(1:6) = sngl(C(i,1:6))
!   call Message%WriteValue('',io_real,6)
! end do







! prepare the HDF5 output file
call timer%makeTimeStamp()
dstr = timer%getDateString()
tstre = timer%getTimeString()

! Create a new file using the default properties.
fname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(enl%datafile)

hdferr =  HDF%createFile(fname)

! write the EMheader to the file
  datagroupname = trim(HDFnames%get_ProgramData())
  call HDF%writeEMheader(EMsoft,dstr, tstrb, tstre, progname, datagroupname)

! add the Duration field to the EMheader group
  hdferr = HDF%openGroup(HDFnames%get_EMheader())
  hdferr = HDF%openGroup(HDFnames%get_ProgramData())

dataset = SC_Duration
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  tstop = 0
  if (g_exists) then
    hdferr = HDF%writeDatasetFloat(dataset, tstop, overwrite)
  else
    hdferr = HDF%writeDatasetFloat(dataset, tstop)
  end if
  call HDF%pop()
  call HDF%pop()

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

! then the remainder of the data in a EMData group
  hdferr = HDF%createGroup(HDFnames%get_EMData())
  hdferr = HDF%createGroup(HDFnames%get_ProgramData())

dataset = 'homographies'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetDoubleArray(dataset, homographies, 8, numpats, overwrite)
  else
    hdferr = HDF%writeDatasetDoubleArray(dataset, homographies, 8, numpats)
  end if

dataset = 'correctedhomographies'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetDoubleArray(dataset, correctedhomographies, 8, numpats, overwrite)
  else
    hdferr = HDF%writeDatasetDoubleArray(dataset, correctedhomographies, 8, numpats)
  end if

dataset = 'Fehat'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetDoubleArray(dataset, Fehat, 3, 3, numpats, overwrite)
  else
    hdferr = HDF%writeDatasetDoubleArray(dataset, Fehat, 3, 3, numpats)
  end if

dataset = 'delta'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetDoubleArray(dataset, delta, 2, numpats, overwrite)
  else
    hdferr = HDF%writeDatasetDoubleArray(dataset, delta, 2, numpats)
  end if

dataset = 'deltaDD'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetDoubleArray(dataset, deltaDD, numpats, overwrite)
  else
    hdferr = HDF%writeDatasetDoubleArray(dataset, deltaDD, numpats)
  end if

dataset = 'alpha'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetDoubleArray(dataset, alpha, numpats, overwrite)
  else
    hdferr = HDF%writeDatasetDoubleArray(dataset, alpha, numpats)
  end if

dataset = 'residuals'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetDoubleArray(dataset, residuals, numpats, overwrite)
  else
    hdferr = HDF%writeDatasetDoubleArray(dataset, residuals, numpats)
  end if

dataset = 'ndp'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetDoubleArray(dataset, normdp, numpats, overwrite)
  else
    hdferr = HDF%writeDatasetDoubleArray(dataset, normdp, numpats)
  end if

dataset = 'nit'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetIntegerArray(dataset, nit, numpats, overwrite)
  else
    hdferr = HDF%writeDatasetIntegerArray(dataset, nit, numpats)
  end if

call HDF%popall()

call Message%printMessage(' results stored in '//trim(fname))

call closeFortranHDFInterface()

end associate 

end subroutine HREBSD_DIC_


end module mod_HREBSDDIC
