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

module mod_MCOpenCL
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/04/20
  !!
  !! class definition for the EMMCOpenCL program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMMCOpenCL program
type, public :: MCOpenCLNameListType
  integer(kind=irg) :: stdout
  integer(kind=irg) :: numsx
  integer(kind=irg) :: ivolx 
  integer(kind=irg) :: ivoly 
  integer(kind=irg) :: ivolz 
  integer(kind=irg) :: globalworkgrpsz
  integer(kind=irg) :: num_el
  integer(kind=irg) :: totnum_el
  integer(kind=irg) :: multiplier
  integer(kind=irg) :: devid
  integer(kind=irg) :: platid
  real(kind=sgl)    :: ivolstepx 
  real(kind=sgl)    :: ivolstepy 
  real(kind=sgl)    :: ivolstepz 
  real(kind=dbl)    :: sig
  real(kind=dbl)    :: sigstart
  real(kind=dbl)    :: sigend
  real(kind=dbl)    :: sigstep
  real(kind=dbl)    :: omega
  real(kind=dbl)    :: EkeV
  real(kind=dbl)    :: Ehistmin
  real(kind=dbl)    :: Ebinsize
  real(kind=dbl)    :: depthmax
  real(kind=dbl)    :: depthstep
  real(kind=dbl)    :: thickness
  real(kind=dbl)    :: radius
  real(kind=dbl)    :: incloc
  character(3)      :: Notify
  character(4)      :: MCmode
  character(fnlen)  :: xtalname
  character(fnlen)  :: dataname
  character(fnlen)  :: mode
end type MCOpenCLNameListType

! class definition
type, public :: MCOpenCL_T
private 
  character(fnlen)            :: nmldeffile = 'EMMCOpenCL.nml'
  type(MCOpenCLNameListType)  :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: MCOpenCL_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: MCOpenCL => MCOpenCL_

end type MCOpenCL_T

! the constructor routine for this class 
interface MCOpenCL_T
  module procedure MCOpenCL_constructor
end interface MCOpenCL_T

contains

!--------------------------------------------------------------------------
type(MCOpenCL_T) function MCOpenCL_constructor( nmlfile ) result(MCOpenCL)
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! constructor for the MCOpenCL_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call MCOpenCL%readNameList(nmlfile)

end function MCOpenCL_constructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly, writetofile)
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! read the namelist from an nml file for the MCOpenCL_T Class 

use mod_io 

IMPLICIT NONE 

class(MCOpenCL_T), INTENT(INOUT)            :: self
character(fnlen),INTENT(IN)                 :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)                 :: initonly
 !! fill in the default values only; do not read the file
character(fnlen),INTENT(IN),optional        :: writetofile
 !! file name to which to write out the entire namelist

type(IO_T)                                  :: Message       
logical                                     :: skipread = .FALSE.

integer(kind=irg) :: stdout
integer(kind=irg) :: numsx
integer(kind=irg) :: ivolx 
integer(kind=irg) :: ivoly 
integer(kind=irg) :: ivolz 
integer(kind=irg) :: globalworkgrpsz
integer(kind=irg) :: num_el
integer(kind=irg) :: totnum_el
integer(kind=irg) :: multiplier
integer(kind=irg) :: devid
integer(kind=irg) :: platid
real(kind=sgl)    :: ivolstepx 
real(kind=sgl)    :: ivolstepy 
real(kind=sgl)    :: ivolstepz 
real(kind=dbl)    :: sig
real(kind=dbl)    :: sigstart
real(kind=dbl)    :: sigend
real(kind=dbl)    :: sigstep
real(kind=dbl)    :: omega
real(kind=dbl)    :: EkeV
real(kind=dbl)    :: Ehistmin
real(kind=dbl)    :: Ebinsize
real(kind=dbl)    :: depthmax
real(kind=dbl)    :: depthstep
character(3)      :: Notify
character(4)      :: MCmode
character(fnlen)  :: xtalname
character(fnlen)  :: dataname
character(fnlen)  :: mode

! define the IO namelist to facilitate passing variables to the program.
namelist  / MCCLdata / stdout, xtalname, sigstart, numsx, num_el, globalworkgrpsz, EkeV, multiplier, &
dataname, totnum_el, Ehistmin, Ebinsize, depthmax, depthstep, omega, MCmode, mode, devid, platid, &
sigend, sigstep, sig, Notify, ivolx, ivoly, ivolz, ivolstepx, ivolstepy, ivolstepz

if (present(writetofile)) then
  if (trim(writetofile).ne.'') then 
    xtalname = trim(self%nml%xtalname)
    mode = self%nml%mode
    ivolx = self%nml%ivolx
    ivoly = self%nml%ivoly
    ivolz = self%nml%ivolz
    ivolstepx = self%nml%ivolstepx
    ivolstepy = self%nml%ivolstepy
    ivolstepz = self%nml%ivolstepz
    globalworkgrpsz = self%nml%globalworkgrpsz
    num_el = self%nml%num_el
    totnum_el = self%nml%totnum_el 
    multiplier = self%nml%multiplier
    devid = self%nml%devid 
    platid = self%nml%platid
    sig = self%nml%sig  
    omega = self%nml%omega
    EkeV = self%nml%EkeV 
    Ehistmin = self%nml%Ehistmin
    Ebinsize = self%nml%Ebinsize
    dataname = self%nml%dataname

    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='unknown')
    write(UNIT=dataunit,NML=MCCLdata)
    close(UNIT=dataunit,STATUS='keep')
    return 
  end if 
end if 

! set the input parameters to default values (except for xtalname, which must be present)
stdout = 6
ivolx = 1001
ivoly = 1001
ivolz = 101
numsx = 1501
globalworkgrpsz = 100
num_el = 10
totnum_el = 2000000000
multiplier = 1
devid = 1
platid = 1
ivolstepx = 1.0
ivolstepy = 1.0
ivolstepz = 1.0
sig = 70.D0
sigstart = 70.D0
sigend = 70.D0
sigstep = 1.D0
omega = 0.D0
EkeV = 30.D0
Ehistmin = 5.D0
Ebinsize = 0.5D0
depthmax = 100.D0
depthstep = 1.0D0
Notify = 'Off'
MCmode = 'CSDA'
xtalname = 'undefined'
dataname = 'MCoutput.data'
mode = 'full'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
  open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
  read(UNIT=dataunit,NML=MCCLdata)
  close(UNIT=dataunit,STATUS='keep')

! check for required entries
  if (trim(xtalname).eq.'undefined') then
  call Message%printError('EMMC:',' structure file name is undefined in '//nmlfile)
  end if
end if

! if we get here, then all appears to be ok, and we need to fill in the mcnl fields
self%nml%stdout = stdout
self%nml%numsx = numsx
self%nml%ivolx = ivolx
self%nml%ivoly = ivoly 
self%nml%ivolz = ivolz
self%nml%globalworkgrpsz = globalworkgrpsz
self%nml%num_el = num_el
self%nml%totnum_el = totnum_el
self%nml%multiplier = multiplier
self%nml%devid = devid
self%nml%platid = platid
self%nml%ivolstepx = ivolstepx 
self%nml%ivolstepy = ivolstepy
self%nml%ivolstepz = ivolstepz
self%nml%sigstart = sigstart
self%nml%sigend = sigend
self%nml%sigstep = sigstep
self%nml%sig = sig
self%nml%omega = omega
self%nml%EkeV = EkeV
self%nml%Ehistmin = Ehistmin
self%nml%Ebinsize = Ebinsize
self%nml%depthmax = depthmax
self%nml%depthstep = depthstep
self%nml%Notify= Notify
self%nml%MCmode = MCmode
self%nml%xtalname = xtalname
self%nml%dataname = dataname
self%nml%mode = mode

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! pass the namelist for the MCOpenCL_T Class to the calling program

IMPLICIT NONE 

class(MCOpenCL_T), INTENT(INOUT)          :: self
type(MCOpenCLNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF)
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! write namelist to HDF file

use mod_HDFsupport
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(MCOpenCL_T), INTENT(INOUT)        :: self 
type(HDF_T), INTENT(INOUT)              :: HDF

integer(kind=irg),parameter             :: n_int = 11, n_real_bse1 = 9, n_real_full = 7, n_real_ivol= 6
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=dbl)                          :: io_real_bse1(n_real_bse1), io_real_full(n_real_full), &
                                           io_real_ivol(n_real_ivol)
character(20)                           :: reallist_bse1(n_real_bse1), reallist_full(n_real_full), &
                                           reallist_ivol(n_real_ivol)
character(20)                           :: intlist(n_int)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( mcnl => self%nml )

! create the group for this namelist
groupname = SC_MCCLNameList
hdferr = HDF%createGroup(groupname)

! write all the single integers
io_int = (/ mcnl%stdout, mcnl%numsx, mcnl%globalworkgrpsz, mcnl%num_el, mcnl%totnum_el, mcnl%multiplier, mcnl%devid, &
            mcnl%platid, mcnl%ivolx, mcnl%ivoly, mcnl%ivolz /)
intlist(1) = 'stdout'
intlist(2) = 'numsx'
intlist(3) = 'globalworkgrpsz'
intlist(4) = 'num_el'
intlist(5) = 'totnum_el'
intlist(6) = 'multiplier'
intlist(7) = 'devid'
intlist(8) = 'platid'
intlist(9) = 'ivolx'
intlist(10) = 'ivoly'
intlist(11) = 'ivolz'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single doubles
if (mcnl%mode .eq. 'bse1') then
   io_real_bse1 = (/ mcnl%sigstart, mcnl%sigend, mcnl%sigstep, mcnl%omega, mcnl%EkeV, mcnl%Ehistmin, &
             mcnl%Ebinsize, mcnl%depthmax, mcnl%depthstep /)
   reallist_bse1(1) = 'sigstart'
   reallist_bse1(2) = 'sigend'
   reallist_bse1(3) = 'sigstep'
   reallist_bse1(4) = 'omega'
   reallist_bse1(5) = 'EkeV'
   reallist_bse1(6) = 'Ehistmin'
   reallist_bse1(7) = 'Ebinsize'
   reallist_bse1(8) = 'depthmax'
   reallist_bse1(9) = 'depthstep'
   call HDF%writeNMLdbles(io_real_bse1, reallist_bse1, n_real_bse1)
else if (mcnl%mode .eq. 'full') then
   io_real_full = (/ mcnl%sig, mcnl%omega, mcnl%EkeV, mcnl%Ehistmin, &
             mcnl%Ebinsize, mcnl%depthmax, mcnl%depthstep /)
   reallist_full(1) = 'sig'
   reallist_full(2) = 'omega'
   reallist_full(3) = 'EkeV'
   reallist_full(4) = 'Ehistmin'
   reallist_full(5) = 'Ebinsize'
   reallist_full(6) = 'depthmax'
   reallist_full(7) = 'depthstep'
   call HDF%writeNMLdbles(io_real_full, reallist_full, n_real_full)
else if (mcnl%mode .eq. 'Ivol') then
   io_real_ivol = (/ mcnl%sig, mcnl%omega, mcnl%EkeV, dble(mcnl%ivolstepx), dble(mcnl%ivolstepy), dble(mcnl%ivolstepz) /)
   reallist_ivol(1) = 'sig'
   reallist_ivol(2) = 'omega'
   reallist_ivol(3) = 'EkeV'
   reallist_ivol(4) = 'ivolstepx'
   reallist_ivol(5) = 'ivolstepy'
   reallist_ivol(6) = 'ivolstepz'
   call HDF%writeNMLdbles(io_real_ivol, reallist_ivol, n_real_ivol)
end if

! write all the strings
dataset = SC_MCmode
sval(1) = mcnl%MCmode
hdferr = HDF%writeDatasetStringArray(dataset, sval, 1)
if (hdferr.ne.0) call error_check_(self,'writeHDFNameList: unable to create MCmode dataset',.TRUE., hdferr)

dataset = SC_xtalname
line2(1) = mcnl%xtalname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call error_check_(self,'writeHDFNameList: unable to create xtalname dataset',.TRUE., hdferr)

dataset = SC_dataname
line2(1) = mcnl%dataname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call error_check_(self,'writeHDFNameList: unable to create dataname dataset',.TRUE., hdferr)

dataset = SC_mode
sval(1) = mcnl%mode
hdferr = HDF%writeDatasetStringArray(dataset, sval, 1)
if (hdferr.ne.0) call error_check_(self,'writeHDFNameList: unable to create mode dataset',.TRUE., hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine MCOpenCL_(self, EMsoft, progname)
!! author: MDG 
!! version: 1.0 
!! date: 01/24/20
!!
!! perform the Monte Carlo computations

use mod_EMsoft
use mod_crystallography
use mod_symmetry
use mod_io
use stringconstants
! use initializersHDF
use mod_initializers
use mod_timing
use mod_diffraction
use mod_Lambert
use clfortran
use mod_CLsupport
use HDF5
use mod_HDFsupport
use mod_notifications
use mod_math
use ISO_C_BINDING

IMPLICIT NONE 

class(MCOpenCL_T), INTENT(INOUT)   :: self
type(EMsoft_T), INTENT(INOUT)      :: EMsoft
character(fnlen),INTENT(IN)        :: progname

type(Cell_T)            :: cell
type(DynType)           :: Dyn
type(Timing_T)          :: timer
type(IO_T)              :: Message
type(OpenCL_T)          :: CL
type(Lambert_T)         :: Lambert 
type(HDF_T)             :: HDF
type(SpaceGroup_T)      :: SG
type(Diffraction_T)     :: Diff

integer(kind=irg)       :: numsy        ! number of Lambert map points along y
integer(kind=irg)       :: numEbins     ! number of energy bins
integer(kind=irg)       :: numzbins     ! number of depth bins
integer(kind=irg)       :: nx           ! no. of pixels
integer(kind=irg)       :: j,k,l,ip,istat, ivx, ivy, ivz
integer(kind=ill)       :: i, io_int(1), num_max, totnum_el_nml, multiplier
real(kind=4),target     :: Ze           ! average atomic number
real(kind=4),target     :: density      ! density in g/cm^3
real(kind=4),target     :: at_wt        ! average atomic weight in g/mole
logical                 :: verbose
real(kind=4)            :: dens, avA, avZ, io_real(3), dmin, Radius  ! used with CalcDensity routine
real(kind=4),target     :: EkeV, sig, omega ! input values to the kernel. Can only be real kind=4 otherwise values are not properly passed
integer(kind=ill)       :: totnum_el, bse     ! total number of electrons to simulate and no. of backscattered electrons
integer(kind=4)         :: prime ! input values to the kernel
integer(kind=4),target  :: globalworkgrpsz, num_el, steps ! input values to the kernel
integer(kind=8)         :: size_in_bytes,size_in_bytes_seeds ! size of arrays passed to kernel. Only accepts kind=8 integers by clCreateBuffer etc., so donot change
integer(kind=8),target  :: globalsize(2), localsize(2) ! size of global and local work groups. Again only kind=8 is accepted by clEnqueueNDRangeKernel
character(4)            :: mode
! results from kernel stored here
real(kind=4),allocatable, target :: Lamresx(:), Lamresy(:), Lamresz(:), depthres(:), energyres(:)

! final results stored here
integer(kind=4),allocatable :: accum_e(:,:,:), accum_z(:,:,:,:), accum_xyz(:,:,:), rnseeds(:)
real(kind=sgl),allocatable  :: accumSP(:,:,:)
integer(kind=ill),allocatable :: accum_e_ill(:,:,:)
integer(kind=4),allocatable,target  :: init_seeds(:)
integer(kind=4)         :: idxy(2), iE, px, py, iz, nseeds, hdferr, tstart, tstop ! auxiliary variables
real(kind=4)            :: cxyz(3), edis, xy(2) ! auxiliary variables
integer(kind=irg)       :: xs, ys, zs
real(kind=8)            :: delta,rand, xyz(3)
character(11)           :: dstr
character(15)           :: tstrb
character(15)           :: tstre
logical                 :: f_exists

integer(c_size_t),target       :: slocal(2), localout

! OpenCL variables
integer(c_intptr_t),allocatable, target  :: platform(:)
integer(c_intptr_t),allocatable, target  :: device(:)
integer(c_intptr_t),target     :: context
integer(c_intptr_t),target     :: command_queue
integer(c_intptr_t),target     :: prog
integer(c_intptr_t),target     :: kernel
integer(c_intptr_t),target     :: LamX, LamY, LamZ, depth, energy, seeds
type(c_ptr)                    :: event
integer(c_int32_t)             :: ierr, pcnt, ierr2
integer(c_size_t),target       :: slength
integer(c_intptr_t),target     :: ctx_props(3)
character(3),target            :: kernelname 
character(5),target            :: kernelname2
character(19),target           :: progoptions
character(fnlen),target        :: info ! info about the GPU
integer(c_int64_t)             :: cmd_queue_props

integer, parameter      :: iunit = 10
integer, parameter      :: source_length = 50000
character(len=source_length),target  :: source
character(len=source_length, KIND=c_char),TARGET :: csource
type(c_ptr), target :: psource
integer(c_int)         :: nump, numd, irec, val,val1 ! auxiliary variables
integer(c_size_t)      :: cnum, cnuminfo
character(fnlen)        :: groupname, dataset, instring, dataname, fname, sourcefile, datagroupname, attributename, HDF_FileVersion
integer(kind=irg)       :: numangle, iang

character(fnlen),ALLOCATABLE      :: MessageLines(:)
integer(kind=irg)                 :: NumLines
character(fnlen)                  :: SlackUsername, exectime
character(100)                    :: c

associate (mcnl => self%nml )

numsy = mcnl%numsx

timer = Timing_T(showDateTime=.TRUE.)
tstre = timer%getTimeString()

! get the crystal structure from the *.xtal file
verbose = .TRUE.
dmin = 0.05
val = 0
val1 = 0
call cell%setFileName(mcnl%xtalname)
Diff = Diffraction_T( mcnl%EkeV, cell)
call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, dmin, noLUT=.TRUE., verbose=verbose)

! then calculate density, average atomic number and average atomic weight
call cell%calcDensity()
io_real(1:3) = cell%getDensity()
density = io_real(1)
at_wt = io_real(2)
Ze = io_real(3)
call WriteValue('Density, avA, avZ = ',io_real,3,"(2f10.5,',',f10.5)")
mode = mcnl%mode

if (mode .eq. 'full') then
    steps = 300
else if (mode .eq. 'bse1') then
    steps = 1
else if (mode .eq. 'Ivol') then
    steps = 300
else
    stop 'Unknown mode specified in namelist file'
end if

! various parameters
EkeV = mcnl%EkeV
!sig = mcnl%sig*dtor
omega = mcnl%omega*dtor
globalworkgrpsz = mcnl%globalworkgrpsz
num_el = mcnl%num_el ! no. of electron simulation by one work item
num_max = globalworkgrpsz*globalworkgrpsz*num_el ! total simulation in one loop
totnum_el_nml = mcnl%totnum_el
multiplier =  mcnl%multiplier
totnum_el = totnum_el_nml * multiplier ! total number of electrons to simulate
globalsize = (/ mcnl%globalworkgrpsz, mcnl%globalworkgrpsz /)
numEbins =  int((mcnl%EkeV-mcnl%Ehistmin)/mcnl%Ebinsize)+1
numzbins =  int(mcnl%depthmax/mcnl%depthstep)+1
nx = (mcnl%numsx-1)/2

! allocate result arrays for GPU part
if (mode.eq.'Ivol') then 
  allocate(Lamresx(num_max), Lamresy(num_max), Lamresz(num_max), stat=istat)
  Lamresz = 0.0
else
  allocate(Lamresx(num_max), Lamresy(num_max), depthres(num_max), energyres(num_max), stat=istat)
  depthres = 0.0
  energyres = 0.0
end if
Lamresx = 0.0
Lamresy = 0.0
size_in_bytes = num_max*sizeof(EkeV)
size_in_bytes_seeds = 4*globalworkgrpsz*globalworkgrpsz*sizeof(EkeV)

if (mode .eq. 'bse1') then
    if (mcnl%sigstep .ne. 0.D0) then
       numangle = nint((mcnl%sigend - mcnl%sigstart)/mcnl%sigstep)+1
    else
       call Message%printError('MCOpenCL:','zero step size for sigma values')
    end if
end if

if (mode .eq. 'full') then
   numangle = 1
   allocate(accum_e(numEbins,-nx:nx,-nx:nx),accum_z(numEbins,numzbins,-nx/10:nx/10,-nx/10:nx/10),stat=istat)
   accum_e = 0
   accum_z = 0
else if (mode .eq. 'bse1') then
   allocate(accum_e(numangle,-nx:nx,-nx:nx),accum_z(numangle,numzbins,-nx/10:nx/10,-nx/10:nx/10),stat=istat)
   accum_e = 0
   accum_z = 0
else if (mode .eq. 'Ivol') then
   ivx = (mcnl%ivolx-1)/2
   ivy = (mcnl%ivoly-1)/2
   ivz = mcnl%ivolz
   numangle = 1
   allocate(accum_xyz(-ivx:ivx,-ivy:ivy,ivz),stat=istat)
   accum_xyz = 0
else
   call Message%printError('MCOpenCL:','Unknown mode specified in namelist file')
end if

! changed by MDG [09/01/15] after extensive modifications to Lambert routines
! old code delta = dble(nx)/LPs%sPio2
delta = dble(nx)

!=====================
! INITIALIZATION
!=====================
CL = OpenCL_T()
call CL%init_PDCCQ(platform, nump, mcnl%platid, device, numd, mcnl%devid, info, context, command_queue)

!=====================
! BUILD THE KERNEL
!=====================

! read the source file
if (mode .eq. 'Ivol') then 
  sourcefile = 'EMMCxyz.cl'
else
  sourcefile = 'EMMC.cl'
end if
call Message%printMessage('OpenCL source file set to : '//trim(sourcefile))
call CL%read_source_file(EMsoft, sourcefile, csource, slength)

! create the program
io_int(1) = slength
call WriteValue('Kernel source length (characters) : ',io_int,1)
pcnt = 1
psource = C_LOC(csource)
prog = clCreateProgramWithSource(context, pcnt, C_LOC(psource), C_LOC(slength), ierr)
call CL%error_check('DoMCsimulation:clCreateProgramWithSource', ierr)

! build the program
! progoptions = '-cl-no-signed-zeros'
! ierr = clBuildProgram(prog, numd, C_LOC(device), C_LOC(progoptions), C_NULL_FUNPTR, C_NULL_PTR)
ierr = clBuildProgram(prog, numd, C_LOC(device), C_NULL_PTR, C_NULL_FUNPTR, C_NULL_PTR)

! get the compilation log
ierr2 = clGetProgramBuildInfo(prog, device(mcnl%devid), CL_PROGRAM_BUILD_LOG, sizeof(source), C_LOC(source), cnum)
if(len(trim(source)) > 0) call Message%printMessage(trim(source(1:cnum)),frm='(A)')
call CL%error_check('DoMCsimulation:clBuildProgram', ierr)
call CL%error_check('DoMCsimulation:clGetProgramBuildInfo', ierr2)

! if we get here, then the program build was successful and we can proceed with the creation of the kernel
call Message%printMessage('Program Build Successful... Creating kernel')

! finally get the kernel and release the program
if (mode.eq.'Ivol') then
  kernelname2 = 'MCxyz'//CHAR(0)
  kernel = clCreateKernel(prog, C_LOC(kernelname2), ierr)
  call CL%error_check('DoMCsimulation:clCreateKernel:MCxyz', ierr)
else
  kernelname = 'MC'//CHAR(0)
  kernel = clCreateKernel(prog, C_LOC(kernelname), ierr)
  call CL%error_check('DoMCsimulation:clCreateKernel:MC', ierr)
end if

ierr = clReleaseProgram(prog)
call CL%error_check('DoMCsimulation:clReleaseProgram', ierr)

! get the random number generator seeds for use on the GPU
fname = EMsoft%generateFilePath('Randomseedfilename')
open(unit = iunit, file = trim(fname), form='unformatted', status='old')
read(iunit) nseeds
allocate(rnseeds(nseeds))
read(iunit) rnseeds
close(unit=iunit,status='keep')

if (4*globalworkgrpsz**2 .gt. nseeds) then
  call Message%printMessage('------------------------------')
  io_int(1) = nseeds 
  call Message%WriteValue('Total number of prime number seeds available = ',io_int,1)
  io_int(1) = 4*globalworkgrpsz**2
  call Message%WriteValue('Total number of prime number seeds needed    = ',io_int,1)
  call Message%printMessage(' ')
  call Message%printMessage('Please reduce the globalworkgrpsz parameter or increase the number of seeds')
  call Message%printMessage('in the '//trim(fname)//' file. The total')
  call Message%printMessage('number of prime seeds needed is equal to 4*globalworkgrpsz*globalworkgrpsz.')
  call Message%printError('MCOpenCL:','insufficient prime number seeds')
end if

allocate(init_seeds(4*globalworkgrpsz*globalworkgrpsz),stat=istat)
init_seeds = 0
do i = 1,globalworkgrpsz
    do j = 1,globalworkgrpsz
        do k = 1,4
            init_seeds(4*((i-1)*globalworkgrpsz+(j-1))+k) = rnseeds(4*((i-1)*globalworkgrpsz+j)+k)
        end do
    end do
end do

! create device memory buffers
LamX = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, C_NULL_PTR, ierr)
call CL%error_check('DoMCsimulation:clCreateBuffer:LamX', ierr)

LamY = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, C_NULL_PTR, ierr)
call CL%error_check('DoMCsimulation:clCreateBuffer:LamY', ierr)

if (mode.eq.'Ivol') then 
  LamZ = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, C_NULL_PTR, ierr)
  call CL%error_check('DoMCsimulation:clCreateBuffer:LamZ', ierr)
else
  depth = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, C_NULL_PTR, ierr)
  call CL%error_check('DoMCsimulation:clCreateBuffer:depth', ierr)

  energy = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, C_NULL_PTR, ierr)
  call CL%error_check('DoMCsimulation:clCreateBuffer:energy', ierr)
end if

seeds = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, C_NULL_PTR, ierr)
call CL%error_check('DoMCsimulation:clCreateBuffer:seeds', ierr)

ierr = clEnqueueWriteBuffer(command_queue, seeds, CL_TRUE, 0_8, size_in_bytes_seeds, C_LOC(init_seeds(1)), &
                            0, C_NULL_PTR, C_NULL_PTR)
call CL%error_check('DoMCsimulation:clEnqueueWriteBuffer', ierr)

if (mode .eq. 'bse1') then
   call Message%printMessage('Monte Carlo mode set to bse1. Calculating statistics for tilt series...',frm='(A/)')
else if (mode .eq. 'full') then
   call Message%printMessage('Monte Carlo mode set to full. Performing full calculation...',frm='(A/)')
else if (mode .eq. 'Ivol') then 
   call Message%printMessage('Monte Carlo mode set to Ivol. Performing full calculation...',frm='(A/)')
else
   call Message%printError('DoMCSimulation','Unknown mode specified in namelist/json file')
end if

call timer%Time_tick(tstart)

angleloop: do iang = 1,numangle

    if (mode .eq. 'bse1') then
        io_int(1) = iang
        call Message%Writevalue('Angle loop #',io_int,1,'(I3)')
        sig = (mcnl%sigstart + (iang-1)*mcnl%sigstep)*dtor
    else if (mode .eq. 'full') then
        sig = mcnl%sig*dtor
    else if (mode .eq. 'Ivol') then
        sig = mcnl%sig*dtor
    end if

    mainloop: do i = 1,(totnum_el/num_max+1)

! set the kernel arguments
if (mode.ne.'Ivol') then 
        ierr = clSetKernelArg(kernel, 0, sizeof(LamX), C_LOC(LamX))
        call CL%error_check('DoMCsimulation:clSetKernelArg:LamX', ierr)

        ierr = clSetKernelArg(kernel, 1, sizeof(LamY), C_LOC(LamY))
        call CL%error_check('DoMCsimulation:clSetKernelArg:LamY', ierr)

        ierr = clSetKernelArg(kernel, 2, sizeof(EkeV), C_LOC(EkeV))
        call CL%error_check('DoMCsimulation:clSetKernelArg:EkeV', ierr)

        ierr = clSetKernelArg(kernel, 3, sizeof(globalworkgrpsz), C_LOC(globalworkgrpsz))
        call CL%error_check('DoMCsimulation:clSetKernelArg:globalworkgrpsz', ierr)

        ierr = clSetKernelArg(kernel, 4, sizeof(Ze), C_LOC(Ze))
        call CL%error_check('DoMCsimulation:clSetKernelArg:Ze', ierr)

        ierr = clSetKernelArg(kernel, 5, sizeof(density), C_LOC(density))
        call CL%error_check('DoMCsimulation:clSetKernelArg:density', ierr)

        ierr = clSetKernelArg(kernel, 6, sizeof(at_wt), C_LOC(at_wt))
        call CL%error_check('DoMCsimulation:clSetKernelArg:at_wt', ierr)

        ierr = clSetKernelArg(kernel, 7, sizeof(num_el), C_LOC(num_el))
        call CL%error_check('DoMCsimulation:clSetKernelArg:num_el', ierr)

        ierr = clSetKernelArg(kernel, 8, sizeof(seeds), C_LOC(seeds))
        call CL%error_check('DoMCsimulation:clSetKernelArg:seeds', ierr)

        ierr = clSetKernelArg(kernel, 9, sizeof(sig), C_LOC(sig))
        call CL%error_check('DoMCsimulation:clSetKernelArg:sig', ierr)

        ierr = clSetKernelArg(kernel, 10, sizeof(omega), C_LOC(omega))
        call CL%error_check('DoMCsimulation:clSetKernelArg:omega', ierr)

        ierr = clSetKernelArg(kernel, 11, sizeof(depth), C_LOC(depth))
        call CL%error_check('DoMCsimulation:clSetKernelArg:depth', ierr)

        ierr = clSetKernelArg(kernel, 12, sizeof(energy), C_LOC(energy))
        call CL%error_check('DoMCsimulation:clSetKernelArg:energy', ierr)

        ierr = clSetKernelArg(kernel, 13, sizeof(steps), C_LOC(steps))
        call CL%error_check('DoMCsimulation:clSetKernelArg:steps', ierr)
else
        ierr = clSetKernelArg(kernel, 0, sizeof(LamX), C_LOC(LamX))
        call CL%error_check('DoMCsimulation:clSetKernelArg:LamX', ierr)

        ierr = clSetKernelArg(kernel, 1, sizeof(LamY), C_LOC(LamY))
        call CL%error_check('DoMCsimulation:clSetKernelArg:LamY', ierr)

        ierr = clSetKernelArg(kernel, 2, sizeof(LamZ), C_LOC(LamZ))
        call CL%error_check('DoMCsimulation:clSetKernelArg:LamZ', ierr)

        ierr = clSetKernelArg(kernel, 3, sizeof(EkeV), C_LOC(EkeV))
        call CL%error_check('DoMCsimulation:clSetKernelArg:EkeV', ierr)

        ierr = clSetKernelArg(kernel, 4, sizeof(globalworkgrpsz), C_LOC(globalworkgrpsz))
        call CL%error_check('DoMCsimulation:clSetKernelArg:globalworkgrpsz', ierr)

        ierr = clSetKernelArg(kernel, 5, sizeof(Ze), C_LOC(Ze))
        call CL%error_check('DoMCsimulation:clSetKernelArg:Ze', ierr)

        ierr = clSetKernelArg(kernel, 6, sizeof(density), C_LOC(density))
        call CL%error_check('DoMCsimulation:clSetKernelArg:density', ierr)

        ierr = clSetKernelArg(kernel, 7, sizeof(at_wt), C_LOC(at_wt))
        call CL%error_check('DoMCsimulation:clSetKernelArg:at_wt', ierr)

        ierr = clSetKernelArg(kernel, 8, sizeof(num_el), C_LOC(num_el))
        call CL%error_check('DoMCsimulation:clSetKernelArg:num_el', ierr)

        ierr = clSetKernelArg(kernel, 9, sizeof(seeds), C_LOC(seeds))
        call CL%error_check('DoMCsimulation:clSetKernelArg:seeds', ierr)

        ierr = clSetKernelArg(kernel, 10, sizeof(sig), C_LOC(sig))
        call CL%error_check('DoMCsimulation:clSetKernelArg:sig', ierr)

        ierr = clSetKernelArg(kernel, 11, sizeof(omega), C_LOC(omega))
        call CL%error_check('DoMCsimulation:clSetKernelArg:omega', ierr)

        ierr = clSetKernelArg(kernel, 12, sizeof(steps), C_LOC(steps))
        call CL%error_check('DoMCsimulation:clSetKernelArg:steps', ierr)
end if

! execute the kernel
        ierr = clEnqueueNDRangeKernel(command_queue, kernel, 2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
                                      0, C_NULL_PTR, C_NULL_PTR)
        call CL%error_check('DoMCsimulation:clEnqueueNDRangeKernel', ierr)

! wait for the commands to finish
        ierr = clFinish(command_queue)
        call CL%error_check('DoMCsimulation:clFinish', ierr)

! read the resulting vector from device memory
if (mode.ne.'Ivol') then
        ierr = clEnqueueReadBuffer(command_queue,LamX,CL_TRUE,0_8,size_in_bytes,C_LOC(Lamresx(1)),0,C_NULL_PTR,C_NULL_PTR)
        call CL%error_check('DoMCsimulation:clEnqueueReadBuffer:Lamresx', ierr)
        ierr = clEnqueueReadBuffer(command_queue,LamY,CL_TRUE,0_8,size_in_bytes,C_LOC(Lamresy(1)),0,C_NULL_PTR,C_NULL_PTR)
        call CL%error_check('DoMCsimulation:clEnqueueReadBuffer:Lamresy', ierr)
        ierr = clEnqueueReadBuffer(command_queue,depth,CL_TRUE,0_8,size_in_bytes,C_LOC(depthres(1)),0,C_NULL_PTR,C_NULL_PTR)
        call CL%error_check('DoMCsimulation:clEnqueueReadBuffer:depthres', ierr)
        ierr = clEnqueueReadBuffer(command_queue,energy,CL_TRUE,0_8,size_in_bytes,C_LOC(energyres(1)),0,C_NULL_PTR,C_NULL_PTR)
        call CL%error_check('DoMCsimulation:clEnqueueReadBuffer:energyres', ierr)
else
        ierr = clEnqueueReadBuffer(command_queue,LamX,CL_TRUE,0_8,size_in_bytes,C_LOC(Lamresx(1)),0,C_NULL_PTR,C_NULL_PTR)
        call CL%error_check('DoMCsimulation:clEnqueueReadBuffer:Lamresx', ierr)
        ierr = clEnqueueReadBuffer(command_queue,LamY,CL_TRUE,0_8,size_in_bytes,C_LOC(Lamresy(1)),0,C_NULL_PTR,C_NULL_PTR)
        call CL%error_check('DoMCsimulation:clEnqueueReadBuffer:Lamresy', ierr)
        ierr = clEnqueueReadBuffer(command_queue,LamZ,CL_TRUE,0_8,size_in_bytes,C_LOC(Lamresz(1)),0,C_NULL_PTR,C_NULL_PTR)
        call CL%error_check('DoMCsimulation:clEnqueueReadBuffer:Lamresz', ierr)
end if

!    call clEnqueueReadBuffer(command_queue, seeds, cl_bool(.true.), 0_8, size_in_bytes_seeds, init_seeds(1), ierr)
        if (mode .eq. 'full') then
           subloopfull: do j = 1, num_max

               if ((Lamresx(j) .ne. -10.0) .and. (Lamresy(j) .ne. -10.0) &
               .and. (depthres(j) .ne. 10.0) .and. (energyres(j) .ne. 0.0) &
               .and. .not.isnan(Lamresx(j)) .and. .not.isnan(Lamresy(j))) then
! and get the nearest pixel [ take into account reversal of coordinate frame (x,y) -> (y,-x) ]
                   if ((nint(delta*Lamresy(j)) .eq. 0.0) .and. (nint(-delta*Lamresx(j)) .eq. 0.0)) then
                       val1 = val1 + 1
                   end if

                   val = val + 1
                   idxy = (/ nint(delta*Lamresy(j)), nint(-delta*Lamresx(j)) /)

                   if (maxval(abs(idxy)).le.nx) then
! If Ec larger than Emin, then we should count this electron
                       if (energyres(j).gt.mcnl%Ehistmin) then

                           iE = nint((energyres(j)-mcnl%Ehistmin)/mcnl%Ebinsize)+1
! first add this electron to the correct exit distance vs. energy bin (coarser than the angular plot)
                           edis = abs(depthres(j))  ! distance from last scattering point to surface along trajectory
                           iz = nint(edis/mcnl%depthstep) +1
                           if ( (iz.gt.0).and.(iz.le.numzbins) ) then

                               px = nint(idxy(1)/10.0)
                               py = nint(idxy(2)/10.0)
                               accum_z(iE,iz,px,py) = accum_z(iE,iz,px,py) + 1

                           end if
! then add it to the modified Lambert accumulator array.
                           accum_e(iE,idxy(1),idxy(2)) = accum_e(iE,idxy(1),idxy(2)) + 1
                       end if
                   end if
               end if
           end do subloopfull
        end if
 
        if (mode .eq. 'bse1') then
           subloopbse1: do j = 1, num_max

               if ((Lamresx(j) .ne. -10.0) .and. (Lamresy(j) .ne. -10.0) &
               .and. (depthres(j) .ne. 10.0) .and. (energyres(j) .ne. 0.0) &
               .and. .not.isnan(Lamresx(j)) .and. .not.isnan(Lamresy(j))) then
! and get the nearest pixel [ take into account reversal of coordinate frame (x,y) -> (y,-x) ]
                   if ((nint(delta*Lamresy(j)) .eq. 0.0) .and. (nint(-delta*Lamresx(j)) .eq. 0.0)) then
                       val1 = val1 + 1
                   end if

                   val = val + 1
                   idxy = (/ nint(delta*Lamresy(j)), nint(-delta*Lamresx(j)) /)

                   if (maxval(abs(idxy)).le.nx) then
! first add this electron to the correct exit distance vs. sigma (coarser than the angular plot)
                       edis = abs(depthres(j))  ! distance from last scattering point to surface along trajectory
                       iz = nint(edis/mcnl%depthstep) +1
                       if ( (iz.gt.0).and.(iz.le.numzbins) ) then
                           px = nint(idxy(1)/10.0)
                           py = nint(idxy(2)/10.0)
                           accum_z(iang,iz,px,py) = accum_z(iang,iz,px,py) + 1

                       end if
! then add it to the modified Lambert accumulator array.
                       accum_e(iang,idxy(1),idxy(2)) = accum_e(iang,idxy(1),idxy(2)) + 1
                   end if
               end if
           end do subloopbse1
        end if

! this simulation mode produces a 3D histogram of the interaction volume with potentially different step
! sizes in the plane as opposed to the depth direction.  The scaling parameters are new name list parameters
! (new as of version 5.0.2).  
        if (mode .eq. 'Ivol') then
           subloopIvol: do j = 1, num_max
               if ((Lamresx(j) .ne. -100000.0) .and. (Lamresy(j) .ne. -100000.0) &
               .and. (Lamresz(j) .ne. -100000.0) &
               .and. .not.isnan(Lamresx(j)) .and. .not.isnan(Lamresy(j)) .and. .not.isnan(Lamresz(j))) then
                  xs = nint( Lamresx(j) / mcnl%ivolstepx )
                  ys = nint( Lamresy(j) / mcnl%ivolstepy )
                  zs = nint( Lamresz(j) / mcnl%ivolstepz )
                  if ((abs(xs).lt.ivx).and.(abs(ys).lt.ivy).and.(zs.lt.ivz) ) accum_xyz(xs,ys,zs+1) = accum_xyz(xs,ys,zs+1) + 1
               end if
           end do subloopIvol
        end if

        if (mod(i,50).eq.0) then
            io_int(1) = i*num_max
            call Message%WriteValue(' Total number of electrons incident = ',io_int, 1, "(I15)")
            if (mode .eq. 'bse1') then
                io_int(1) = sum(accum_e(iang,:,:))
                call Message%WriteValue(' Number of BSE1 electrons = ',io_int, 1, "(I15)")
            else if(mode .eq. 'full') then
                allocate(accum_e_ill(numEbins,-nx:nx,-nx:nx),stat=istat)
                accum_e_ill = accum_e
                io_int(1) = sum(accum_e_ill)
                deallocate(accum_e_ill)
                call Message%WriteValue(' Number of BSE electrons = ',io_int, 1, "(I15)")
            else if(mode .eq. 'Ivol') then
                io_int(1) = sum(accum_xyz)
                call Message%WriteValue(' Number of electrons in interaction volume = ',io_int, 1, "(I15)")
            else
                call Message%printError('DoMCSimulations','Unknown mode specified in namelist/json file')
            end if

            if (i.eq.50) then

            end if

        end if


    end do mainloop
! and write some infgormation to the console

    io_int(1) = totnum_el
    call Message%WriteValue('Total number of incident electrons = ',io_int,1,'(I15)')
    if (mode .eq. 'bse1') then
        io_int(1) = sum(accum_e(iang,:,:))
        call Message%WriteValue('Total number of BSE1 electrons = ',io_int,1,'(I15)')
        bse = sum(accum_e(iang,:,:))
        io_real(1) = dble(bse)/dble(totnum_el)
        call Message%WriteValue('Backscatter yield = ',io_real,1,'(F15.6)')
    else if (mode .eq. 'full') then
! note that we need to prevent integer overflows !
        allocate(accum_e_ill(numEbins,-nx:nx,-nx:nx),stat=istat)
        accum_e_ill = accum_e
        io_int(1) = sum(accum_e_ill)
        deallocate(accum_e_ill)
        call Message%WriteValue('Total number of BSE electrons = ',io_int,1,'(I15)')
        io_real(1) = dble(io_int(1))/dble(totnum_el)
        call Message%WriteValue('Backscatter yield = ',io_real,1,'(F15.6)')
    else if (mode .eq. 'Ivol') then
        io_int(1) = sum(accum_xyz)
        call Message%WriteValue('Total number of electrons in interaction volume = ',io_int,1,'(I15)')
    else 
        call Message%printError('DoMCSimulations','Unknown mode specified in namelist/json file')
    end if
 

end do angleloop

if (mode .eq. 'full') then
! get stereographic projections from the accum_e array
  allocate(accumSP(numEbins,-nx:nx,-nx:nx))
  Radius = 1.0
  do i=-nx,nx 
    do j=-nx,nx 
      Lambert = Lambert_T( xy = (/ float(i), float(j) /) / float(nx) )
      ierr = Lambert%StereoGraphicInverse(xyz, dble(Radius))
      xyz = xyz/vecnorm(xyz)
      if (ierr.ne.0) then 
        accumSP(1:numEbins,i,j) = 0.0
      else
        accumSP(1:numEbins,i,j) = InterpolateLambert(xyz, accum_e, nx, numEbins)
      end if
    end do
  end do
end if

call timer%Time_tock() 
io_int(1) = timer%getInterval()
call WriteValue('Total execution time [s] = ',io_int,1)

io_int(1) = totnum_el/num_max
totnum_el = (io_int(1)+1)*num_max

timer = Timing_T()
dstr = timer%getDateString()
tstre = timer%getTimeString()

! initialize the HDF class
HDF = HDF_T()

! get the filename; if it already exists, then delete it and create a new one
dataname = EMsoft%generateFilePath('EMdatapathname', mcnl%dataname)
inquire(file=trim(dataname), exist=f_exists)

if (f_exists) then
  open(unit=dataunit, file=trim(dataname), status='old',form='unformatted')
  close(unit=dataunit, status='delete')
end if

! Create a new file using the default properties.
hdferr = HDF%createFile(dataname)

! write the EMheader to the file
datagroupname = 'MCOpenCL'
call HDF%writeEMheader(dstr, tstrb, tstre, progname, datagroupname)

! add the CrystalData group at the top level of the file
call cell%SaveDataHDF(SG, EMsoft)

! create a namelist group to write all the namelist files into
groupname = SC_NMLfiles
hdferr = HDF%createGroup(groupname)

! read the text file and write the array to the file
dataset = SC_MCOpenCLNML
hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%nmldeffile)

! leave this group
call HDF%pop()

! create a namelist group to write all the namelist files into
groupname = SC_NMLparameters
hdferr = HDF%createGroup(groupname)
call self%writeHDFNameList(HDF)

! leave this group
call HDF%pop()

! then the remainder of the data in a EMData group
groupname = SC_EMData
hdferr = HDF%createGroup(groupname)

! here we add the data groupname MCOpenCL and we attach to it a HDF_FileVersion attribute 
hdferr = HDF%createGroup(datagroupname)
HDF_FileVersion = '4.0'
attributename = SC_HDFFileVersion
hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)

! =====================================================
! The following write commands constitute HDF_FileVersion = 4.0
! =====================================================

dataset = SC_numzbins
hdferr = HDF%writeDatasetInteger(dataset, numzbins)

! modified using multiplier
dataset = SC_totnumel
hdferr = HDF%writeDatasetInteger(dataset, mcnl%totnum_el)

dataset = SC_multiplier
hdferr = HDF%writeDatasetInteger(dataset, mcnl%multiplier)

if (mode .eq. 'full') then

dataset = SC_numEbins
    hdferr = HDF%writeDatasetInteger(dataset, numEbins)

dataset = SC_accume
    hdferr = HDF%writeDatasetIntegerArray(dataset, accum_e, numEbins, 2*nx+1, 2*nx+1)

dataset = SC_accumz
    hdferr = HDF%writeDatasetIntegerArray(dataset, accum_z, numEbins, numzbins, 2*(nx/10)+1, 2*(nx/10)+1)

dataset = SC_accumSP
    hdferr = HDF%writeDatasetFloatArray(dataset, accumSP, numEbins, 2*nx+1, 2*nx+1)

else if (mode .eq. 'bse1') then

dataset = SC_numangle
    hdferr = HDF%writeDatasetInteger(dataset, numangle)

dataset = SC_accume
    hdferr = HDF%writeDatasetIntegerArray(dataset, accum_e, numangle, 2*nx+1, 2*nx+1)

dataset = SC_accumz
    hdferr = HDF%writeDatasetIntegerArray(dataset, accum_z, numangle, numzbins, 2*(nx/10)+1, 2*(nx/10)+1)

else if (mode .eq. 'Ivol') then

dataset = SC_accumxyz
    hdferr = HDF%writeDatasetIntegerArray(dataset, accum_xyz, 2*ivx+1, 2*ivy+1, ivz)

end if

! =====================================================
! end of HDF_FileVersion = 4.0 write statements
! =====================================================

call HDF_pop(.TRUE.)
!
!=====================
! RELEASE EVERYTHING
!=====================

ierr = clReleaseKernel(kernel)
call CL%error_check('DoMCsimulation:clReleaseKernel', ierr)
ierr = clReleaseCommandQueue(command_queue)
call CL%error_check('DoMCsimulation:clReleaseCommandQueue', ierr)
ierr = clReleaseContext(context)
call CL%error_check('DoMCsimulation:clReleaseContext', ierr)
ierr = clReleaseMemObject(LamX)
call CL%error_check('DoMCsimulation:clReleaseMemObject:LamX', ierr)
ierr = clReleaseMemObject(LamY)
call CL%error_check('DoMCsimulation:clReleaseMemObject:LamY', ierr)
if (mode.eq.'Ivol') then 
  ierr = clReleaseMemObject(LamZ)
  call CL%error_check('DoMCsimulation:clReleaseMemObject:LamZ', ierr)
else
  ierr = clReleaseMemObject(depth)
  call CL%error_check('DoMCsimulation:clReleaseMemObject:depth', ierr)
  ierr = clReleaseMemObject(energy)
  call CL%error_check('DoMCsimulation:clReleaseMemObject:energy', ierr)
end if
ierr = clReleaseMemObject(seeds)
call CL%error_check('DoMCsimulation:clReleaseMemObject:seeds', ierr)


! if requested, we notify the user that this program has completed its run
if (trim(EMsoft%getConfigParameter('Notify')).ne.'Off') then
  if (trim(mcnl%Notify).eq.'On') then 
    NumLines = 3
    allocate(MessageLines(NumLines))

    call hostnm(c)

    MessageLines(1) = 'EMMCOpenCL program has ended successfully'
    MessageLines(2) = 'Monte Carlo data stored in '//trim(dataname)
    write (exectime,"(I10)") tstop  
    MessageLines(3) = 'Total execution time [s]: '//trim(exectime)
    SlackUsername = 'EMsoft on '//trim(c)
    i = PostMessage(EMsoft, MessageLines, NumLines, SlackUsername)
  end if
end if

end associate

end subroutine MCOpenCL_



end module mod_MCOpenCL