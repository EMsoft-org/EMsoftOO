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

module mod_MCOpenCL
  !! author: MDG
  !! version: 1.0
  !! date: 02/04/20
  !!
  !! class definition for the EMMCOpenCL program

use mod_kinds
use mod_global
use mod_MCfiles
use mod_platformsupport

IMPLICIT NONE

! class definition
type, public :: MCOpenCL_T
private
  character(fnlen)            :: nmldeffile = 'EMMCOpenCL.nml'
  type(MCOpenCLNameListType)  :: nml

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: MCOpenCL_
  procedure, pass(self) :: testMCOpenCLWrapper_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: MCOpenCL => MCOpenCL_
  generic, public :: testMCOpenCLWrapper => testMCOpenCLWrapper_

end type MCOpenCL_T

! the constructor routine for this class
interface MCOpenCL_T
  module procedure MCOpenCL_constructor
end interface MCOpenCL_T

contains

!--------------------------------------------------------------------------
type(MCOpenCL_T) function MCOpenCL_constructor( nmlfile ) result(MCOpenCL)
!DEC$ ATTRIBUTES DLLEXPORT :: MCOpenCL_constructor
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
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
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
real(kind=dbl)    :: thickness
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
                       sigend, sigstep, sig, Notify, ivolx, ivoly, ivolz, ivolstepx, ivolstepy, ivolstepz, thickness

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
thickness = 100.D0
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
self%nml%thickness = thickness
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
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
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
subroutine MCOpenCL_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: MCOpenCL_
!! author: MDG
!! version: 1.0
!! date: 01/24/20
!!
!! perform the Monte Carlo computations for EBSD, ECP, TKD, or Ivol modes

use mod_EMsoft
use mod_crystallography
use mod_symmetry
use mod_io
use stringconstants
use mod_initializers
use mod_timing
use mod_diffraction
use mod_Lambert
use clfortran
use mod_CLsupport
use HDF5
use mod_HDFsupport
use mod_HDFnames
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
type(HDFnames_T)        :: HDFnames
type(SpaceGroup_T)      :: SG
type(Diffraction_T)     :: Diff
type(MCfile_T)          :: MCFT

integer(kind=irg)       :: numsy        ! number of Lambert map points along y
integer(kind=irg)       :: nx           ! no. of pixels
integer(kind=irg)       :: j,k,l,ip,istat, ivx, ivy, ivz
integer(kind=ill)       :: i, io_int(1), num_max, totnum_el_nml, multiplier
real(kind=4),target     :: Ze           ! average atomic number
real(kind=4),target     :: density      ! density in g/cm^3
real(kind=4),target     :: at_wt        ! average atomic weight in g/mole
logical                 :: verbose
real(kind=4)            :: dens, avA, avZ, io_real(3), dmin, Radius  ! used with CalcDensity routine
real(kind=8)            :: io_dble(3)
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
! integer(kind=4),allocatable :: accum_e(:,:,:), accum_z(:,:,:,:), accum_xyz(:,:,:), rnseeds(:)
! real(kind=sgl),allocatable  :: accumSP(:,:,:)
integer(kind=4),allocatable :: rnseeds(:)
integer(kind=ill),allocatable :: accum_e_ill(:,:,:)
integer(kind=4),allocatable,target  :: init_seeds(:)
integer(kind=4)         :: idxy(2), iE, px, py, iz, nseeds, hdferr, tstart, tstop ! auxiliary variables
real(kind=4)            :: cxyz(3), edis, xy(2) ! auxiliary variables
integer(kind=irg)       :: xs, ys, zs
real(kind=8)            :: delta,rand, xyz(3)
real(kind=4), target    :: thickness
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
integer(kind=irg)       :: iang

character(fnlen),ALLOCATABLE      :: MessageLines(:)
integer(kind=irg)                 :: NumLines
character(fnlen)                  :: SlackUsername, exectime
character(100)                    :: c
integer(kind=4)                   :: hnStat


associate (mcnl => self%nml, MCDT => MCFT%MCDT )

numsy = mcnl%numsx

timer = Timing_T()
tstrb = timer%getTimeString()

call openFortranHDFInterface()
HDF = HDF_T()
HDFnames = HDFnames_T()

! set the HDF group names for reading the MC input file
call HDFnames%set_ProgramData(SC_MCOpenCL)
call HDFnames%set_NMLlist(SC_MCCLNameList)
call HDFnames%set_NMLfilename(SC_MCOpenCLNML)

! get the crystal structure from the *.xtal file
verbose = .TRUE.
dmin = 0.05
val = 0
val1 = 0
call cell%setFileName(mcnl%xtalname)
call Diff%setV(mcnl%EkeV)
call Diff%setrlpmethod('WK')
call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, dmin, noLUT=.TRUE., verbose=verbose, useHDF = HDF)

! then calculate density, average atomic number and average atomic weight
call cell%calcDensity()
io_dble(1:3) = cell%getDensity()
density = io_dble(1)
at_wt = io_dble(2)
Ze = io_dble(3)
call Message%WriteValue(' Density, avA, avZ = ',io_dble,3,"(/2f10.5,',',f10.5)")
mode = mcnl%mode

if ( (mode .eq. 'full') .or. (mode .eq. 'foil') ) then
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
omega = mcnl%omega*dtor
globalworkgrpsz = mcnl%globalworkgrpsz
num_el = mcnl%num_el ! no. of electron simulation by one work item
num_max = globalworkgrpsz*globalworkgrpsz*num_el ! total simulation in one loop
totnum_el_nml = mcnl%totnum_el
multiplier =  mcnl%multiplier
totnum_el = totnum_el_nml * multiplier ! total number of electrons to simulate
globalsize = (/ mcnl%globalworkgrpsz, mcnl%globalworkgrpsz /)
MCDT%numEbins =  int((mcnl%EkeV-mcnl%Ehistmin)/mcnl%Ebinsize)+1
MCDT%numzbins =  int(mcnl%depthmax/mcnl%depthstep)+1
nx = (mcnl%numsx-1)/2
if (mode .eq. 'foil') then
  thickness = sngl(mcnl%thickness)
end if

! allocate result arrays for GPU part
if (mode.eq.'Ivol') then
  allocate(Lamresx(num_max), Lamresy(num_max), Lamresz(num_max), stat=istat)
  Lamresx = 0.0
  Lamresy = 0.0
  Lamresz = 0.0
else
  allocate(Lamresx(num_max), Lamresy(num_max), depthres(num_max), energyres(num_max), stat=istat)
  Lamresx = 0.0
  Lamresy = 0.0
  depthres = 0.0
  energyres = 0.0
end if
size_in_bytes = num_max*sizeof(EkeV)
size_in_bytes_seeds = 4*globalworkgrpsz*globalworkgrpsz*sizeof(EkeV)

if (mode .eq. 'bse1') then
    if (mcnl%sigstep .ne. 0.D0) then
       MCDT%numangle = nint((mcnl%sigend - mcnl%sigstart)/mcnl%sigstep)+1
    else
       call Message%printError('MCOpenCL:','zero step size for sigma values')
    end if
end if

if ( (mode .eq. 'full') .or. (mode .eq. 'foil') ) then
   MCDT%numangle = 1
   allocate(MCDT%accum_e(MCDT%numEbins,-nx:nx,-nx:nx), &
            MCDT%accum_z(MCDT%numEbins,MCDT%numzbins,-nx/10:nx/10,-nx/10:nx/10),stat=istat)
   MCDT%accum_e = 0
   MCDT%accum_z = 0
else if (mode .eq. 'bse1') then
   allocate(MCDT%accum_e(MCDT%numangle,-nx:nx,-nx:nx), &
            MCDT%accum_z(MCDT%numangle,MCDT%numzbins,-nx/10:nx/10,-nx/10:nx/10),stat=istat)
   MCDT%accum_e = 0
   MCDT%accum_z = 0
else if (mode .eq. 'Ivol') then
   ivx = (mcnl%ivolx-1)/2
   ivy = (mcnl%ivoly-1)/2
   ivz = mcnl%ivolz
   MCDT%numangle = 1
   allocate(MCDT%accum_xyz(-ivx:ivx,-ivy:ivy,ivz),stat=istat)
   MCDT%accum_xyz = 0
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
else if (mode .eq. 'foil') then
  sourcefile = 'EMMCfoil.cl'
else
  sourcefile = 'EMMC.cl'
end if
call Message%printMessage(' OpenCL source file set to : '//trim(sourcefile))
call CL%read_source_file(EMsoft, sourcefile, csource, slength)

! create the program
io_int(1) = slength
call Message%WriteValue(' Kernel source length (characters) : ',io_int,1)
pcnt = 1
psource = C_LOC(csource)
prog = clCreateProgramWithSource(context, pcnt, C_LOC(psource), C_LOC(slength), ierr)
call CL%error_check('DoMCsimulation:clCreateProgramWithSource', ierr)

! build the program
ierr = clBuildProgram(prog, numd, C_LOC(device), C_NULL_PTR, C_NULL_FUNPTR, C_NULL_PTR)

! get the compilation log
ierr2 = clGetProgramBuildInfo(prog, device(mcnl%devid), CL_PROGRAM_BUILD_LOG, sizeof(source), C_LOC(source), cnum)
if(len(trim(source)) > 0) call Message%printMessage(trim(source(1:cnum)),frm='(A)')
call CL%error_check('DoMCsimulation:clBuildProgram', ierr)
call CL%error_check('DoMCsimulation:clGetProgramBuildInfo', ierr2)

! if we get here, then the program build was successful and we can proceed with the creation of the kernel
call Message%printMessage(' Program Build Successful... Creating kernel')

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
   call Message%printMessage(' Monte Carlo mode set to bse1. Calculating statistics for tilt series...',frm='(/A/)')
else if (mode .eq. 'full') then
   call Message%printMessage(' Monte Carlo mode set to full. Performing full calculation...',frm='(/A/)')
else if (mode .eq. 'Ivol') then
   call Message%printMessage(' Monte Carlo mode set to Ivol. Performing full calculation...',frm='(/A/)')
else if (mode .eq. 'foil') then
   call Message%printMessage(' Monte Carlo mode set to foil. Performing full calculation...',frm='(/A/)')
else
   call Message%printError('DoMCSimulation','Unknown mode specified in namelist/json file')
end if

call timer%Time_tick()

angleloop: do iang = 1,MCDT%numangle

    if (mode .eq. 'bse1') then
        io_int(1) = iang
        call Message%Writevalue(' Angle loop #',io_int,1,'(I3)')
        sig = (mcnl%sigstart + (iang-1)*mcnl%sigstep)*dtor
    else
        sig = mcnl%sig*dtor
    end if

    mainloop: do i = 1,(totnum_el/num_max+1)

! set the kernel arguments
if (mode.ne.'Ivol') then
  if (mode .ne. 'foil') then
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
        ierr = clSetKernelArg(kernel, 0, sizeof(EkeV), C_LOC(EkeV))
        call CL%error_check('DoMCsimulation:clSetKernelArg:EkeV', ierr)

        ierr = clSetKernelArg(kernel, 1, sizeof(globalworkgrpsz), C_LOC(globalworkgrpsz))
        call CL%error_check('DoMCsimulation:clSetKernelArg:globalworkgrpsz', ierr)

        ierr = clSetKernelArg(kernel, 2, sizeof(Ze), C_LOC(Ze))
        call CL%error_check('DoMCsimulation:clSetKernelArg:Ze', ierr)

        ierr = clSetKernelArg(kernel, 3, sizeof(density), C_LOC(density))
        call CL%error_check('DoMCsimulation:clSetKernelArg:density', ierr)

        ierr = clSetKernelArg(kernel, 4, sizeof(at_wt), C_LOC(at_wt))
        call CL%error_check('DoMCsimulation:clSetKernelArg:at_wt', ierr)

        ierr = clSetKernelArg(kernel, 5, sizeof(num_el), C_LOC(num_el))
        call CL%error_check('DoMCsimulation:clSetKernelArg:num_el', ierr)

        ierr = clSetKernelArg(kernel, 6, sizeof(seeds), C_LOC(seeds))
        call CL%error_check('DoMCsimulation:clSetKernelArg:seeds', ierr)

        ierr = clSetKernelArg(kernel, 7, sizeof(sig), C_LOC(sig))
        call CL%error_check('DoMCsimulation:clSetKernelArg:sig', ierr)

        ierr = clSetKernelArg(kernel, 8, sizeof(omega), C_LOC(omega))
        call CL%error_check('DoMCsimulation:clSetKernelArg:omega', ierr)

        ierr = clSetKernelArg(kernel, 9, sizeof(depth), C_LOC(depth))
        call CL%error_check('DoMCsimulation:clSetKernelArg:depth', ierr)

        ierr = clSetKernelArg(kernel, 10, sizeof(energy), C_LOC(energy))
        call CL%error_check('DoMCsimulation:clSetKernelArg:energy', ierr)

        ierr = clSetKernelArg(kernel, 11, sizeof(steps), C_LOC(steps))
        call CL%error_check('DoMCsimulation:clSetKernelArg:steps', ierr)

        ierr = clSetKernelArg(kernel, 12, sizeof(thickness), C_LOC(thickness))
        call CL%error_check('DoMCsimulation:clSetKernelArg:thickness', ierr)

        ierr = clSetKernelArg(kernel, 13, sizeof(LamX), C_LOC(LamX))
        call CL%error_check('DoMCsimulation:clSetKernelArg:LamXSH', ierr)

        ierr = clSetKernelArg(kernel, 14, sizeof(LamY), C_LOC(LamY))
        call CL%error_check('DoMCsimulation:clSetKernelArg:LamYSH', ierr)
      end if
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
        if ( (mode .eq. 'full') .or. (mode .eq. 'foil') ) then
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
! occasionally, the GPU returns an energy value that is larger than the maximum energy 
! likely du to an error in the OpenCL code; we intercept these rare events here
                       if ( (energyres(j).gt.mcnl%Ehistmin).and.(energyres(j).le.mcnl%EkeV) ) then

                           iE = nint((energyres(j)-mcnl%Ehistmin)/mcnl%Ebinsize)+1
! first add this electron to the correct exit distance vs. energy bin (coarser than the angular plot)
                           edis = abs(depthres(j))  ! distance from last scattering point to surface along trajectory
                           iz = nint(edis/mcnl%depthstep) +1
                           if ( (iz.gt.0).and.(iz.le.MCDT%numzbins) ) then

                               px = nint(idxy(1)/10.0)
                               py = nint(idxy(2)/10.0)
                               MCDT%accum_z(iE,iz,px,py) = MCDT%accum_z(iE,iz,px,py) + 1

                           end if
! then add it to the modified Lambert accumulator array.
                           MCDT%accum_e(iE,idxy(1),idxy(2)) = MCDT%accum_e(iE,idxy(1),idxy(2)) + 1
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
                       if ( (iz.gt.0).and.(iz.le.MCDT%numzbins) ) then
                           px = nint(idxy(1)/10.0)
                           py = nint(idxy(2)/10.0)
                           MCDT%accum_z(iang,iz,px,py) = MCDT%accum_z(iang,iz,px,py) + 1

                       end if
! then add it to the modified Lambert accumulator array.
                       MCDT%accum_e(iang,idxy(1),idxy(2)) = MCDT%accum_e(iang,idxy(1),idxy(2)) + 1
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
                  if ((abs(xs).lt.ivx).and.(abs(ys).lt.ivy).and.(zs.lt.ivz) ) then
                    MCDT%accum_xyz(xs,ys,zs+1) = MCDT%accum_xyz(xs,ys,zs+1) + 1
                  end if
               end if
           end do subloopIvol
        end if

        if (mod(i,50).eq.0) then
            io_int(1) = i*num_max
            call Message%WriteValue(' Total number of electrons incident = ',io_int, 1, "(I15)")
            if (mode .eq. 'bse1') then
                io_int(1) = sum(MCDT%accum_e(iang,:,:))
                call Message%WriteValue(' Number of BSE1 electrons = ',io_int, 1, "(I15)")
            else if ( (mode .eq. 'full') .or. (mode .eq. 'foil') ) then
                allocate(accum_e_ill(MCDT%numEbins,-nx:nx,-nx:nx),stat=istat)
                accum_e_ill = MCDT%accum_e
                io_int(1) = sum(accum_e_ill)
                deallocate(accum_e_ill)
                call Message%WriteValue(' Number of BSE electrons = ',io_int, 1, "(I15)")
            else if(mode .eq. 'Ivol') then
                io_int(1) = sum(MCDT%accum_xyz)
                call Message%WriteValue(' Number of electrons in interaction volume = ',io_int, 1, "(I15)")
            else
                call Message%printError('DoMCSimulations','Unknown mode specified in namelist/json file')
            end if
        end if

    end do mainloop
! and write some infgormation to the console

    io_int(1) = totnum_el
    call Message%printMessage(' ')
    call Message%WriteValue(' Total number of incident electrons = ',io_int,1,'(I15)')
    if (mode .eq. 'bse1') then
        io_int(1) = sum(MCDT%accum_e(iang,:,:))
        call Message%WriteValue(' Total number of BSE1 electrons = ',io_int,1,'(I15)')
        bse = sum(MCDT%accum_e(iang,:,:))
        io_real(1) = dble(bse)/dble(totnum_el)
        call Message%WriteValue(' Backscatter yield = ',io_real,1,'(F15.6)')
    else if ( (mode .eq. 'full') .or. (mode .eq. 'foil') ) then
! note that we need to prevent integer overflows !
        allocate(accum_e_ill(MCDT%numEbins,-nx:nx,-nx:nx),stat=istat)
        accum_e_ill = MCDT%accum_e
        io_int(1) = sum(accum_e_ill)
        deallocate(accum_e_ill)
        call Message%WriteValue(' Total number of BSE electrons = ',io_int,1,'(I15)')
        io_real(1) = dble(io_int(1))/dble(totnum_el)
        call Message%WriteValue(' Backscatter yield = ',io_real,1,'(F15.6)')
    else if (mode .eq. 'Ivol') then
        io_int(1) = sum(MCDT%accum_xyz)
        call Message%WriteValue(' Total number of electrons in interaction volume = ',io_int,1,'(I15)')
    else
        call Message%printError('DoMCSimulations','Unknown mode specified in namelist/json file')
    end if


end do angleloop

if ( (mode .eq. 'full') .or. (mode .eq. 'foil') ) then
! get stereographic projections from the accum_e array
  allocate(MCDT%accumSP(MCDT%numEbins,-nx:nx,-nx:nx))
  Radius = 1.0
  do i=-nx,nx
    do j=-nx,nx
      Lambert = Lambert_T( xy = (/ float(i), float(j) /) / float(nx) )
      ierr = Lambert%StereoGraphicInverse(xyz, dble(Radius))
      xyz = xyz/vecnorm(xyz)
      if (ierr.ne.0) then
        MCDT%accumSP(1:MCDT%numEbins,i,j) = 0.0
      else
        MCDT%accumSP(1:MCDT%numEbins,i,j) = InterpolateLambert(xyz, MCDT%accum_e, nx, MCDT%numEbins)
      end if
    end do
  end do
end if

call timer%Time_tock()
io_real(1) = timer%getInterval()
call Message%printMessage(' ')
call Message%WriteValue(' Total execution time [s] = ',io_real,1)

io_int(1) = totnum_el/num_max
totnum_el = (io_int(1)+1)*num_max

timer = Timing_T()
dstr = timer%getDateString()
tstre = timer%getTimeString()

! set a few variables, copy the namelist to the MCFT class
! and save the data to an HDF5 file
HDF = HDF_T()
call MCFT%copynml(mcnl)
call MCFT%writeMCfile(EMsoft, cell, SG, HDF, HDFnames, progname, dstr, tstrb, tstre)

call closeFortranHDFInterface()

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

    hnStat = system_hostnm(c)

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

!--------------------------------------------------------------------------
subroutine testMCOpenCLWrapper_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: testMCOpenCLWrapper_
!! author: MDG
!! version: 1.0
!! date: 04/16/21
!!
!! driver routine to test the EMsoftCgetMCOpenCL routine in mod_SEMCLwrappers.f90

use mod_EMsoft
use mod_crystallography
use mod_symmetry
use mod_io
use stringconstants
use mod_initializers
use mod_diffraction
use mod_Lambert
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_notifications
use mod_math
use ISO_C_BINDING
use mod_SEMCLwrappers

IMPLICIT NONE

class(MCOpenCL_T), INTENT(INOUT)   :: self
type(EMsoft_T), INTENT(INOUT)      :: EMsoft
character(fnlen),INTENT(IN)        :: progname

type(Cell_T)            :: cell
type(DynType)           :: Dyn
type(IO_T)              :: Message
type(HDF_T)             :: HDF
type(HDFnames_T)        :: HDFnames
type(SpaceGroup_T)      :: SG
type(Diffraction_T)     :: Diff
type(MCfile_T)          :: MCFT

integer(kind=irg)       :: numsy        ! number of Lambert map points along y
integer(kind=irg)       :: nx           ! no. of pixels
integer(kind=irg)       :: j,k,l,ip,istat, ivx, ivy, ivz, numsites
integer(kind=ill)       :: i, io_int(1), num_max, totnum_el_nml, multiplier, s4(4), s3(3)
real(kind=4),target     :: Ze           ! average atomic number
real(kind=4),target     :: density      ! density in g/cm^3
real(kind=4),target     :: at_wt        ! average atomic weight in g/mole
logical                 :: verbose
real(kind=4)            :: dens, avA, avZ, io_real(3), dmin, Radius  ! used with CalcDensity routine
real(kind=8)            :: io_dble(3)
real(kind=4),target     :: EkeV, sig, omega ! input values to the kernel. Can only be real kind=4 otherwise values are not properly passed
integer(kind=ill)       :: totnum_el, bse     ! total number of electrons to simulate and no. of backscattered electrons
integer(kind=4)         :: prime ! input values to the kernel
integer(kind=4),target  :: globalworkgrpsz, num_el, steps ! input values to the kernel
integer(kind=8)         :: size_in_bytes,size_in_bytes_seeds ! size of arrays passed to kernel. Only accepts kind=8 integers by clCreateBuffer etc., so donot change
integer(kind=8),target  :: globalsize(2), localsize(2) ! size of global and local work groups. Again only kind=8 is accepted by clEnqueueNDRangeKernel
character(4)            :: mode

integer(kind=4)         :: idxy(2), iE, px, py, iz, nseeds, hdferr, tstart, tstop ! auxiliary variables
real(kind=4)            :: cxyz(3), edis, xy(2) ! auxiliary variables
integer(kind=irg)       :: xs, ys, zs
real(kind=8)            :: delta,rand, xyz(3)
real(kind=4), target    :: thickness
character(11)           :: dstr
character(15)           :: tstrb
character(15)           :: tstre
logical                 :: f_exists

integer(c_size_t),target       :: slocal(2), localout

integer(c_int)         :: nump, numd, irec, val,val1 ! auxiliary variables
integer(c_size_t)      :: cnum, cnuminfo
character(fnlen)       :: s, pdesc, pname, outname, dataset
integer(kind=irg)      :: iang

real(kind=sgl), allocatable     :: atompos(:,:)
integer(kind=irg),allocatable   :: atomtypes(:)
real(kind=sgl),allocatable      :: atpos(:,:)
integer(kind=irg),allocatable   :: attp(:)
real(kind=sgl)                  :: latparm(6)
integer(c_int32_t)              :: ipar(wraparraysize)
real(kind=sgl)                  :: fpar(wraparraysize)
character(kind=c_char, len=1)   :: spar(wraparraysize*fnlen)
integer(c_size_t)               :: objAddress
character(len=1)                :: cancel


associate (mcnl => self%nml, MCDT => MCFT%MCDT )

numsy = mcnl%numsx

call openFortranHDFInterface()
HDF = HDF_T()
HDFnames = HDFnames_T()

! set the HDF group names for reading the MC input file
call HDFnames%set_ProgramData(SC_MCOpenCL)
call HDFnames%set_NMLlist(SC_MCCLNameList)
call HDFnames%set_NMLfilename(SC_MCOpenCLNML)

! get the crystal structure from the *.xtal file
verbose = .TRUE.
dmin = 0.05
val = 0
val1 = 0
call cell%setFileName(mcnl%xtalname)
call Diff%setV(mcnl%EkeV)
call Diff%setrlpmethod('WK')
call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, dmin, noLUT=.TRUE., verbose=verbose, useHDF = HDF)
numsites = cell%getNatomtype()

MCDT%numEbins =  int((mcnl%EkeV-mcnl%Ehistmin)/mcnl%Ebinsize)+1
MCDT%numzbins =  int(mcnl%depthmax/mcnl%depthstep)+1

! next we set up the input parameters for the EMsoftCgetMCOpenCL routine 
ipar = 0
ipar(1)  = (numsy-1)/2
ipar(2)  = mcnl%globalworkgrpsz
ipar(3)  = mcnl%num_el
ipar(4)  = mcnl%totnum_el
ipar(5)  = mcnl%multiplier
ipar(6)  = mcnl%devid
ipar(7)  = mcnl%platid
ipar(8)  = SG%getSpaceGroupXtalSystem()
ipar(9)  = cell%getNatomtype()
ipar(10) = SG%getSpaceGroupNumber()
ipar(11) = SG%getSpaceGroupSetting()
ipar(12) = MCFT%getnumEbins()
ipar(13) = MCFT%getnumzbins()
ipar(14) = 1
ipar(15) = 1
ipar(16) = ipar(1)/10

fpar = 0
fpar(1)  = mcnl%sig
fpar(2)  = mcnl%omega
fpar(3)  = mcnl%EkeV
fpar(4)  = mcnl%Ehistmin
fpar(5)  = mcnl%Ebinsize
fpar(6)  = mcnl%depthmax
fpar(7)  = mcnl%depthstep
fpar(8)  = mcnl%sigstart
fpar(9)  = mcnl%sigend
fpar(10) = mcnl%sigstep

! the calling program passes a c-string array spar that we need to convert to the 
! standard EMsoft config structure for use inside this routine
pname = ''
pdesc = ''
EMsoft = EMsoft_T( pname, pdesc, silent=.TRUE. )
! copy the necessary strings into the spar array
spar = ' '
s = trim(EMsoft%getConfigParameter('OpenCLpathname'))
j = 22*fnlen+1
k = 1
do i=j,j+len(trim(s)) 
  spar(i) = s(k:k)
  k=k+1
end do 
spar(j+len(trim(s))+1) = C_NULL_CHAR

s = trim(EMsoft%getConfigParameter('Randomseedfilename'))
j = 25*fnlen+1
k = 1
do i=j,j+len(trim(s)) 
  spar(i) = s(k:k)
  k=k+1
end do 
spar(j+len(trim(s))+1) = C_NULL_CHAR

atomtypes = cell%getAtomtype()
atompos = cell%getAsymPosData()
latparm = cell%getLatParm()
allocate(atpos(numsites,5), attp(numsites))
atpos = atompos(1:numsites,1:5)
attp = atomtypes(1:numsites) 

cancel = char(0)
objAddress = 0_c_size_t

MCDT%numangle = 1
nx = ipar(1)
allocate(MCDT%accum_e(MCDT%numEbins,-nx:nx,-nx:nx), &
         MCDT%accum_z(MCDT%numEbins,MCDT%numzbins,-nx/10:nx/10,-nx/10:nx/10),stat=istat)
MCDT%accum_e = 0
MCDT%accum_z = 0

! and call the routine 
write (*,*) 'calling C wrapper'
call EMsoftCgetMCOpenCL(ipar, fpar, spar, atpos, attp, latparm, MCDT%accum_e, MCDT%accum_z, C_NULL_FUNPTR, objAddress, cancel)
write (*,*)  ' done.'

!=============================================
! create a simple HDF5 output file
!=============================================
HDF = HDF_T()

! Open an existing file or create a new file using the default properties.
outname = EMsoft%generateFilePath('EMdatapathname', mcnl%dataname)
hdferr =  HDF%createFile(outname)

! create the hyperslabs and write zeroes to them for now
s3 = shape(MCDT%accum_e)
nx = s3(2)

dataset = SC_accume
    hdferr = HDF%writeDatasetIntegerArray(dataset, MCDT%accum_e, ipar(12), nx, nx)

s4 = shape(MCDT%accum_z)
nx = s4(3)

dataset = SC_accumz
    hdferr = HDF%writeDatasetIntegerArray(dataset, MCDT%accum_z, ipar(12), ipar(13), nx, nx)

call HDF%popall()

call closeFortranHDFInterface()

write (*,*) ' output stored in '//trim(outname)

end associate

end subroutine testMCOpenCLWrapper_



end module mod_MCOpenCL
