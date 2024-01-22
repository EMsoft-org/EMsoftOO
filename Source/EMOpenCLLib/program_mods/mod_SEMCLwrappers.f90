! ###################################################################
! Copyright (c) 2013-2024, Marc De Graef Research Group/Carnegie Mellon University
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
!--------------------------------------------------------------------------
! EMsoft:mod_SEMCLwrappers.f90
!--------------------------------------------------------------------------
!
! MODULE: mod_SEMCLwrappers
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief routines that can be called by external code; all routines requiring HDF are in EMdymodHDF.f90
!
!> @date  01/21/18 MDG 1.0 separated C/C++ callable SEM routines out from original EMdymod module
!--------------------------------------------------------------------------
!
! general information: the ipar and fpar arrays for all the routines that are C-callable
! are identical, so we document here their component definitions; to allow for future expansion, each
! array has 40 entries, of which about half are currently (April 2016) used.
!
! integer(kind=irg) :: ipar(40)  components 
! ipar(1) : nx  = (numsx-1)/2
! ipar(2) : globalworkgrpsz
! ipar(3) : num_el
! ipar(4) : totnum_el
! ipar(5) : multiplier
! ipar(6) : devid
! ipar(7) : platid
! ipar(8) : CrystalSystem
! ipar(9) : Natomtypes
! ipar(10): SpaceGroupNumber
! ipar(11): SpaceGroupSetting
! ipar(12): numEbins
! ipar(13): numzbins
! ipar(14): mcmode  ( 1 = 'full', 2 = 'bse1' )
! ipar(15): numangle
! ipar(16): nxten = nx/10
! the following are only used in the master routine
! ipar(17): npx
! ipar(18): nthreads
! the following are only used in the EBSD pattern routine
! ipar(19): numx of detector pixels
! ipar(20): numy of detector pixels
! ipar(21): number of orientation in quaternion set
! ipar(22): binning factor (0-3)
! ipar(23): binned x-dimension
! ipar(24): binned y-dimension
! ipar(25): anglemode  (0 for quaternions, 1 for Euler angles)
! ipar(26): already initialized 
! ipar(27:40) : 0 (unused for now)

! real(kind=dbl) :: fpar(40)  components
! fpar(1) : sig
! fpar(2) : omega
! fpar(3) : EkeV
! fpar(4) : Ehistmin
! fpar(5) : Ebinsize
! fpar(6) : depthmax
! fpar(7) : depthstep
! fpar(8) : sigstart
! fpar(9) : sigend
! fpar(10): sigstep
! parameters only used in the master pattern routine
! fpar(11) : dmin
! fpar(12) : Bethe  c1
! fpar(13) : Bethe  c2
! fpar(14) : Bethe  c3
! parameters only used in the EBSD pattern routine
! fpar(15): pattern center x
! fpar(16): pattern center y
! fpar(17): scintillator pixel size
! fpar(18): detector tilt angle
! fpar(19): sample-scintillator distance
! fpar(20): beam current
! fpar(21): dwelltime
! fpar(22): gamma value
! fpar(23:40): 0 (unused for now)

! newly added in version 3.2, to facilitate passing EMsoft configuration
! strings back and forth to C/C++ programs that call EMdymod routines...
! character(fnlen)  :: spar(40)   configuration string components
! spar(1): EMsoftpathname
! spar(2): EMXtalFolderpathname
! spar(3): EMdatapathname
! spar(4): EMtmppathname
! spar(5): EMsoftLibraryLocation
! spar(6): EMSlackWebHookURL
! spar(7): EMSlackChannel
! spar(8): UserName
! spar(9): UserLocation
! spar(10): UserEmail
! spar(11): EMNotify
! spar(12): Develop
! spar(13): Release
! spar(14): h5copypath
! spar(15): EMsoftplatform
! spar(16): EMsofttestpath
! spar(17): EMsoftTestingPath
! spar(18): EMsoftversion
! spar(19): Configpath
! spar(20): Templatepathname
! spar(21): Resourcepathname
! spar(22): Homepathname
! spar(23): OpenCLpathname
! spar(24): Templatecodefilename
! spar(25): WyckoffPositionsfilename
! spar(26): Randomseedfilename
! spar(27): EMsoftnativedelimiter
! spar(28:40): '' (unused for now)


!
module mod_SEMCLwrappers

    !--------------------------------------------------------------------------
    ! Callback routine(s) to communicate progress with DREAM.3D package
    
    ! Define interface of call-back routine
    ! arguments are:
    !  objAddress: unique 8-byte integer to identify the calling class in DREAM.3D
    !  patternCompleted: integer indicating the current pattern ID number
    !
    ABSTRACT INTERFACE
       SUBROUTINE ProgressCallBack(objAddress, patternCompleted) bind(C)
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(c_size_t),INTENT(IN), VALUE          :: objAddress
        INTEGER(KIND=4), INTENT(IN), VALUE           :: patternCompleted
       END SUBROUTINE ProgressCallBack
    END INTERFACE
    
    
    ! similar callback routine, with two integer arguments
    ABSTRACT INTERFACE
       SUBROUTINE ProgressCallBack2(objAddress, loopCompleted, totalLoops, bseYield) bind(C)
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(c_size_t),INTENT(IN), VALUE          :: objAddress
        INTEGER(KIND=4), INTENT(IN), VALUE           :: loopCompleted
        INTEGER(KIND=4), INTENT(IN), VALUE           :: totalLoops
        REAL(KIND=4),INTENT(IN), VALUE              :: bseYield
       END SUBROUTINE ProgressCallBack2
    END INTERFACE
    
    ! similar callback routine, with two integer arguments
    ABSTRACT INTERFACE
       SUBROUTINE ProgressCallBack3(objAddress, loopCompleted, totalLoops, EloopCompleted, totalEloops) bind(C)
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(c_size_t),INTENT(IN), VALUE          :: objAddress
        INTEGER(KIND=4), INTENT(IN), VALUE           :: loopCompleted
        INTEGER(KIND=4), INTENT(IN), VALUE           :: totalLoops
        INTEGER(KIND=4), INTENT(IN), VALUE           :: EloopCompleted
        INTEGER(KIND=4), INTENT(IN), VALUE           :: totalELoops
       END SUBROUTINE ProgressCallBack3
    END INTERFACE
    
    !--------------------------------------------------------------------------

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! the routines starting with EMsoftC are callable from C/C++
! programs and can handle progress callback and a cancel request.
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! SUBROUTINE:EMsoftCgetMCOpenCL
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief This subroutine can be called by a C/C++ program as a standalone routine to compute Monte Carlo data
!
!> @details This subroutine provides a method to compute a Monte Carlo data set, normally computed
!> with the EMMCOpenCL.f90 program.  The routine can be called from an external C/C++ program; 
!> the routine provides a callback mechanism to update the calling program about computational 
!> progress, as well as a cancel option.
!>
!> The routine is intended to be called from a C/C++ program, e.g., DREAM.3D.  This routine is a 
!> simplified version of the core of the EMMCOpenCL program. 
!>
!> Since the HDF5 library with fortran90 support can only be a static library on Mac OS X, we must
!> have the calling program read the .xtal HDF5 file and pass the necessary information on to
!> this routine.  This is a workaround until the HDF group fixes the static library issue; DREAM.3D
!> requires a dynamical HDF5 library, so for DREAM.3D and EMsoft to properly work together, the 
!> callable routines in this file may not depend on any HDF code at all, either directly or indirectly.
!>
!> @param ipar array with integer input parameters
!> @param fpar array with float input parameters
!> @param atdata atom coordinate array
!> @param attypes atom type array
!> @param latparm lattice parameter array
!> @param accum_e output array with Monte Carlo energy histogram
!> @param accum_z output array with Monte Carlo depth histogram
!
!> @date 03/08/16 MDG 1.0 original
!> @date 03/19/16 MDG 1.1 corrections to a few variable types
!> @date 04/13/16 MDG 1.2 correction to accum_z array size due to changes in calling DREAM.3D filter
!> @date 04/18/16 MDG 1.3 increased number of entries in ipar, fpar for compatibility with EMsoftCgetEBSDmaster routine
!> @date 04/28/16 MDG 1.4 corrected error in indexing of init_seeds array; caused DREAM.3D to crash randomly
!> @date 11/09/17 MDG 2.0 added spar string array to pass EMsoft configuration strings into the routine
!> @date 04/14/21 MDG 3.0 rewrite for Object Oriented library 
!--------------------------------------------------------------------------
recursive subroutine EMsoftCgetMCOpenCL(ipar, fpar, spar, atompos, atomtypes, latparm, accum_e, accum_z, cproc, &
objAddress, cancel) bind(c, name='EMsoftCgetMCOpenCL')    ! this routine is callable from a C/C++ program
!DEC$ ATTRIBUTES DLLEXPORT :: EMsoftCgetMCOpenCL

! ipar components
! ipar(1) : integer(kind=irg)       :: nx  = (numsx-1)/2
! ipar(2) : integer(kind=irg)       :: globalworkgrpsz
! ipar(3) : integer(kind=irg)       :: num_el
! ipar(4) : integer(kind=irg)       :: totnum_el
! ipar(5) : integer(kind=irg)       :: multiplier
! ipar(6) : integer(kind=irg)       :: devid
! ipar(7) : integer(kind=irg)       :: platid
! ipar(8) : integer(kind=irg)       :: CrystalSystem
! ipar(9) : integer(kind=irg)       :: Natomtypes
! ipar(10): integer(kind=irg)       :: SpaceGroupNumber
! ipar(11): integer(kind=irg)       :: SpaceGroupSetting
! ipar(12): integer(kind=irg)       :: numEbins
! ipar(13): integer(kind=irg)       :: numzbins
! ipar(14): integer(kind=irg)       :: mcmode  ( 1 = 'full', 2 = 'bse1' )
! ipar(15): integer(kind=irg)       :: numangle
! ipar(16): integer(kind=irg)       :: nxten = nx/10
! other entries are not used

! fpar components
! fpar(1) : real(kind=dbl)          :: sig
! fpar(2) : real(kind=dbl)          :: omega
! fpar(3) : real(kind=dbl)          :: EkeV
! fpar(4) : real(kind=dbl)          :: Ehistmin
! fpar(5) : real(kind=dbl)          :: Ebinsize
! fpar(6) : real(kind=dbl)          :: depthmax
! fpar(7) : real(kind=dbl)          :: depthstep
! fpar(8) : real(kind=dbl)          :: sigstart
! fpar(9) : real(kind=dbl)          :: sigend
! fpar(10): real(kind=dbl)          :: sigstep
! other entries are not used

! spar components
! this routine needs the following parameters to be set:
! spar(23): OpenCLpathname
! spar(26): Randomseedfilename
! 

use mod_EMsoft
use mod_crystallography
use mod_symmetry
use mod_io
use stringconstants
! use mod_initializers
use mod_diffraction
use mod_Lambert
use clfortran
use mod_CLsupport
use mod_notifications
use mod_math
use mod_memory
use,INTRINSIC :: ISO_C_BINDING


IMPLICIT NONE

integer(c_int32_t),INTENT(IN)           :: ipar(wraparraysize)
real(kind=sgl),INTENT(IN)               :: fpar(wraparraysize)
character(kind=c_char, len=1), target, INTENT(IN) :: spar(wraparraysize*fnlen)
real(kind=sgl),INTENT(IN)               :: atompos(ipar(9),5)
integer(kind=irg),INTENT(IN)            :: atomtypes(ipar(9))
real(kind=sgl),INTENT(IN)               :: latparm(6)
integer(kind=irg),INTENT(OUT)           :: accum_e(ipar(12),-ipar(1):ipar(1),-ipar(1):ipar(1))
integer(kind=irg),INTENT(OUT)           :: accum_z(ipar(12),ipar(13),-ipar(16):ipar(16),-ipar(16):ipar(16))
TYPE(C_FUNPTR), INTENT(IN), VALUE       :: cproc
integer(c_size_t),INTENT(IN), VALUE     :: objAddress
character(fnlen)                        :: clPath=''
character(len=1),INTENT(IN)             :: cancel

type(EMsoft_T)          :: EMsoft
type(Cell_T)            :: cell
type(DynType)           :: Dyn
type(IO_T)              :: Message
type(OpenCL_T)          :: CL
type(Lambert_T)         :: Lambert
type(SpaceGroup_T)      :: SG
type(Diffraction_T)     :: Diff
type(memory_T)          :: mem
! type(MCfile_T)          :: MCFT

! local variables and parameters
character(4)                            :: mode
character(fnlen)                        :: pname, pdesc
integer(kind=ill)                       :: i=0, j=0, k=0, io_int(1)=0, num_max=0, totnum_el=0, ipg=0, isave=0, istat=0
integer(kind=irg)                       :: nx=0, numEbins=0, numzbins=0, numangle=0, iang=0, cn=0, dn=0, totn=0
integer(kind=irg),target                :: globalworkgrpsz=0, num_el=0, steps=0
integer(kind=8),target                  :: globalsize(2)=0, localsize(2)=0
integer(kind=8)                         :: size_in_bytes=0,size_in_bytes_seeds=0

real(kind=sgl),target                   :: dens=0, avA=0, avZ=0, omega=0, EkeV=0, sig=0, bseyield=0, io_real(3)
real(kind=4),target                     :: density=0, Ze=0, at_wt=0, delta=0
real(kind=4),allocatable, target        :: Lamresx(:), Lamresy(:), depthres(:), energyres(:)

integer(kind=4),allocatable             :: rnseeds(:)
integer(kind=4),allocatable,target      :: init_seeds(:)
integer(kind=4)                         :: idxy(2), iE=0, px=0, py=0, iz=0, nseeds=0, hdferr=0, tstart=0, trimSpace=0 ! auxiliary variables
real(kind=4)                            :: cxyz(3), edis=0, bse=0, xy(2), xs=0, ys=0, zs=0, sclf=0 ! auxiliary variables
logical                                 :: f_exists=.FALSE.
real(kind=dbl)                          :: apositions(maxpasym,5)

! OpenCL variables
integer(c_intptr_t),allocatable, target :: platform(:)
integer(c_intptr_t),allocatable, target :: device(:)
integer(c_intptr_t),target              :: context=0
integer(c_intptr_t),target              :: command_queue=0
integer(c_intptr_t),target              :: prog=0
integer(c_intptr_t),target              :: kernel=0
integer(c_intptr_t),target              :: LamX=0, LamY=0, LamZ=0, depth=0, energy=0, seeds=0
type(c_ptr)                             :: event
integer(c_int32_t)                      :: ierr=0, pcnt=0
integer(c_size_t),target                :: slength=0
integer(c_intptr_t),target              :: ctx_props(3)
character(2),target                     :: kernelname=''
character(19),target                    :: progoptions=''
character(fnlen),target                 :: info='' ! info about the GPU
integer(c_int64_t)                      :: cmd_queue_props=0

integer, parameter                      :: iunit = 10
integer, parameter                      :: source_length = 50000
character(len=source_length),target     :: source=''
character(len=source_length, KIND=c_char),TARGET :: csource=''
type(c_ptr), target                     :: psource
integer(c_int)                          :: nump=0, numd=0, irec=0, val=0,val1=0 ! auxiliary variables
integer(c_size_t)                       :: cnum=0, cnuminfo=0
character(fnlen)                        :: instring='', dataname='', fname='', sourcefile=''
PROCEDURE(ProgressCallBack2), POINTER   :: proc
character(250),target                   :: currentDir=''
character(fnlen)                        :: emmcPath='', outname=''
character(fnlen)                        :: randomSeedPath=''
integer(c_int32_t),target               :: filestat=0
INTEGER(kind=irg)                       :: status
CHARACTER(LEN=30)                       :: Format=''
real(kind=dbl)                          :: dza(3) 

nullify(proc)

! link the proc procedure to the cproc argument
CALL C_F_PROCPOINTER (cproc, proc)

! the calling program passes a c-string array spar that we need to convert to the 
! standard EMsoft config structure for use inside this routine
pname = ''
pdesc = ''
EMsoft = EMsoft_T( pname, pdesc, silent=.TRUE. )
! then we override all values with the ones in the spar array
call EMsoft%C2F_configuration_strings(C_LOC(spar))

! the following is necessitated by the fact that none of this code may 
! depend on HDF5 routines, so we need to cut-and-paste from various 
! other library routines to set things up so that we can compute the 
! density, and the average atomic number and atomic mass...

! copy all the unit cell parameters into the proper fields and compute the 
! density parameters needed by the Monte Carlo routine; 
! initialize cell and SG classes 
cell = Cell_T( dble(latparm) )
SG = SpaceGroup_T( SGnumber = ipar(10), xtalSystem = ipar(8), setting = ipar(11), &
                   dmt = cell%getdmt(), rmt = cell%getrmt() )
! fill in additional symmetry and cell parameters
if ((ipar(10).ge.143).and.(ipar(10).le.167)) then
  call SG%setSpaceGrouptrigonal(.TRUE.)
else
  call SG%setSpaceGrouptrigonal(.FALSE.)
end if 
! atom type and coordinate parameters
call cell%setNatomtype( ipar(9) )
call cell%setAtomtype( atomtypes(1:ipar(9)) )
apositions(1:ipar(9), 1:5) = atompos
call cell%setAtomPos( apositions )
! generate the symmetry operations
call SG%setSpaceGrouphexset( .FALSE. )
if (ipar(8).eq.4) call SG%setSpaceGrouphexset( .TRUE. )
if ((ipar(8).eq.5).AND.(ipar(11).ne.2)) call SG%setSpaceGrouphexset( .TRUE. )
! Get the symmorphic space group corresponding to the point group
! of the actual space group
 ipg = SG%getPGnumber()
! if the actual group is also the symmorphic group, then both 
! steps can be done simultaneously, otherwise two calls to 
! GenerateSymmetry are needed.
 if (SGPG(ipg).eq.ipar(10)) then
  call SG%GenerateSymmetry( .TRUE., cell%getdmt(), cell%getrmt() )
 else
  isave = SG%getSpaceGroupNumber()
  call SG%setSpaceGroupNumber( SGPG(ipg) )
  call SG%GenerateSymmetry( .TRUE., cell%getdmt(), cell%getrmt() )
  call SG%setSpaceGroupNumber( int(isave) )
  call SG%GenerateSymmetry( .FALSE., cell%getdmt(), cell%getrmt() )
 end if
! next we get all the atom positions
call cell%CalcPositions( SG, 'v' )

! and now we have all we need to compute the density, average A and average Z
call cell%CalcDensity()
dza = cell%getDensity()

! and copy these values into the desired variables
density = dza(1) 
Ze = dza(2)
at_wt = dza(3)

! define a number of parameters
steps = 300
mode = 'full'
if (ipar(14).ne.1) mode = 'bse1'   
EkeV = sngl(fpar(3))
omega = sngl(fpar(2))*dtor
globalworkgrpsz = ipar(2)
num_el = int(ipar(3))  ! no. of electron simulation by one work item
num_max = globalworkgrpsz*globalworkgrpsz*num_el ! total simulation in one loop
totnum_el = ipar(4) * ipar(5) ! total number of electrons to simulate
globalsize = (/ globalworkgrpsz, globalworkgrpsz /)
numEbins =  int(ipar(12))
numzbins =  int(ipar(13))
nx = int(ipar(1))
delta = dble(nx)
size_in_bytes = num_max*sizeof(EkeV)
size_in_bytes_seeds = 4*globalworkgrpsz*globalworkgrpsz*sizeof(EkeV)
numangle = int(ipar(15))

! next allocate and initialize a couple of arrays
mem = memory_T()
call mem%alloc(Lamresx, (/ int(num_max) /), 'Lamresx', 0.0)
call mem%alloc(Lamresy, (/ int(num_max) /), 'Lamresy', 0.0)
call mem%alloc(depthres, (/ int(num_max) /), 'depthres', 0.0)
call mem%alloc(energyres, (/ int(num_max) /), 'energyres', 0.0)
accum_e = 0
accum_z = 0

!======================
! OpenCL INITIALIZATION
!======================
CL = OpenCL_T( verb = .FALSE. )
call CL%init_PDCCQ(platform, nump, int(ipar(7)), device, numd, int(ipar(6)), info, context, command_queue)

!=====================
! BUILD THE KERNEL
!=====================
! read the source file
sourcefile='/EMMC.cl'
emmcPath=trim(EMsoft%getConfigParameter('OpenCLpathname'))//trim(sourcefile)
emmcPath=EMsoft%toNativePath(emmcPath)

call CL%read_source_file(EMsoft, emmcPath, csource, slength)

! create the program
pcnt = 1
psource = C_LOC(csource)
prog = clCreateProgramWithSource(context, pcnt, C_LOC(psource), C_LOC(slength), ierr)

! build the program
ierr = clBuildProgram(prog, numd, C_LOC(device), C_NULL_PTR, C_NULL_FUNPTR, C_NULL_PTR)
if (ierr.le.0) then
  ierr = clGetProgramBuildInfo(prog, device(ipar(6)), CL_PROGRAM_BUILD_LOG, sizeof(source), C_LOC(source), cnum)
endif

! get the compilation log
ierr = clGetProgramBuildInfo(prog, device(ipar(6)), CL_PROGRAM_BUILD_LOG, sizeof(source), C_LOC(source), cnum)

! if we get here, then the program build was successful and we can proceed with the creation of the kernel
! finally get the kernel and release the program
kernelname = 'MC'//CHAR(0)
kernel = clCreateKernel(prog, C_LOC(kernelname), ierr)
ierr = clReleaseProgram(prog)

open(unit = iunit, file = trim(EMsoft%toNativePath(EMsoft%getConfigParameter('Randomseedfilename'))), &
     form='unformatted', status='old')
read(iunit) nseeds
call mem%alloc(rnseeds, (/ nseeds /), 'rnseeds')
read(iunit) rnseeds
close(unit=iunit,status='keep')

! the next error needs to be checked in the calling program
!if (globalworkgrpsz**2 .gt. nseeds) call FatalError('EMMCOpenCL:','insufficient prime numbers')

call mem%alloc(init_seeds, (/ 4*globalworkgrpsz*globalworkgrpsz /),'init_seeds', 0)
do i = 1,globalworkgrpsz
    do j = 1,globalworkgrpsz
        do k = 1,4
            init_seeds(4*((i-1)*globalworkgrpsz+(j-1))+k) = rnseeds(4*((i-1)*globalworkgrpsz+j)+k)
        end do
    end do
end do

! create device memory buffers
LamX = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, C_NULL_PTR, ierr)

LamY = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, C_NULL_PTR, ierr)

depth = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, C_NULL_PTR, ierr)

energy = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, C_NULL_PTR, ierr)

seeds = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, C_NULL_PTR, ierr)

ierr = clEnqueueWriteBuffer(command_queue, seeds, CL_TRUE, 0_8, size_in_bytes_seeds, C_LOC(init_seeds(1)), &
                            0, C_NULL_PTR, C_NULL_PTR)

! set the callback parameters
dn = 1
cn = dn
totn = numangle * (totnum_el/num_max+1)

! loop over angles (used for BSE1, single run for full)
angleloop: do iang = 1,numangle

  if (mode .eq. 'bse1') then
    sig = (fpar(8) + (iang-1)*fpar(10))*dtor
  else 
    sig = fpar(1)*dtor
  end if

  mainloop: do i = 1,(totnum_el/num_max+1)

! set the kernel arguments
    ierr = clSetKernelArg(kernel, 0, sizeof(LamX), C_LOC(LamX))

    ierr = clSetKernelArg(kernel, 1, sizeof(LamY), C_LOC(LamY))

    ierr = clSetKernelArg(kernel, 2, sizeof(EkeV), C_LOC(EkeV))

    ierr = clSetKernelArg(kernel, 3, sizeof(globalworkgrpsz), C_LOC(globalworkgrpsz))

    ierr = clSetKernelArg(kernel, 4, sizeof(Ze), C_LOC(Ze))

    ierr = clSetKernelArg(kernel, 5, sizeof(density), C_LOC(density))

    ierr = clSetKernelArg(kernel, 6, sizeof(at_wt), C_LOC(at_wt))

    ierr = clSetKernelArg(kernel, 7, sizeof(num_el), C_LOC(num_el))

    ierr = clSetKernelArg(kernel, 8, sizeof(seeds), C_LOC(seeds))

    ierr = clSetKernelArg(kernel, 9, sizeof(sig), C_LOC(sig))

    ierr = clSetKernelArg(kernel, 10, sizeof(omega), C_LOC(omega))

    ierr = clSetKernelArg(kernel, 11, sizeof(depth), C_LOC(depth))

    ierr = clSetKernelArg(kernel, 12, sizeof(energy), C_LOC(energy))

    ierr = clSetKernelArg(kernel, 13, sizeof(steps), C_LOC(steps))

! execute the kernel
    ierr = clEnqueueNDRangeKernel(command_queue, kernel, 2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
                                  0, C_NULL_PTR, C_NULL_PTR)

! wait for the commands to finish
    ierr = clFinish(command_queue)

! read the resulting vector from device memory
    ierr = clEnqueueReadBuffer(command_queue,LamX,CL_TRUE,0_8,size_in_bytes,C_LOC(Lamresx(1)),0,C_NULL_PTR,C_NULL_PTR)
    ierr = clEnqueueReadBuffer(command_queue,LamY,CL_TRUE,0_8,size_in_bytes,C_LOC(Lamresy(1)),0,C_NULL_PTR,C_NULL_PTR)
    ierr = clEnqueueReadBuffer(command_queue,depth,CL_TRUE,0_8,size_in_bytes,C_LOC(depthres(1)),0,C_NULL_PTR,C_NULL_PTR)
    ierr = clEnqueueReadBuffer(command_queue,energy,CL_TRUE,0_8,size_in_bytes,C_LOC(energyres(1)),0,C_NULL_PTR,C_NULL_PTR)

    if (mode .eq. 'full') then
      val = 0
      subloopfull: do j = 1, num_max
        if ((Lamresx(j) .ne. -10.0) .and. (Lamresy(j) .ne. -10.0) &
          .and. (depthres(j) .ne. 10.0) .and. (energyres(j) .ne. 0.0) &
          .and. .not.isnan(Lamresx(j)) .and. .not.isnan(Lamresy(j))) then
! and get the nearest pixel [ take into account reversal of coordinate frame (x,y) -> (y,-x) ]

             idxy = (/ nint(delta*Lamresy(j)), nint(-delta*Lamresx(j)) /)

             if (maxval(abs(idxy)).le.nx) then
! If Ec larger than Emin, then we should count this electron
               if (energyres(j).gt.fpar(4)) then

                 val = val + 1
                 iE = nint((energyres(j)-fpar(4))/fpar(5))+1
! first add this electron to the correct exit distance vs. energy bin (coarser than the angular plot)
                 edis = abs(depthres(j))  ! distance from last scattering point to surface along trajectory
                 iz = nint(edis/fpar(7)) +1
                 if ( (iz.gt.0).and.(iz.le.ipar(13)) ) then

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
            iz = nint(edis/fpar(7)) +1
            if ( (iz.gt.0).and.(iz.le.ipar(13)) ) then
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

! has the cancel flag been set by the calling program ?
    if(cancel.ne.char(0)) then 
      EXIT angleloop
    end if

! update the progress counter and report it to the calling program via the proc callback routine
    if(objAddress.ne.0) then
      bseyield = 100.0*float(sum(accum_e))/float(i*num_max)
      call proc(objAddress, cn, totn, bseyield)
      cn = cn+dn
    end if

! write (*,*) ' completed ',i,' out of ', (totnum_el/num_max+1)
  end do mainloop
end do angleloop 

!=====================
! RELEASE EVERYTHING
!=====================

ierr = clReleaseKernel(kernel)
ierr = clReleaseCommandQueue(command_queue)
ierr = clReleaseContext(context)
ierr = clReleaseMemObject(LamX)
ierr = clReleaseMemObject(LamY)
ierr = clReleaseMemObject(depth)
ierr = clReleaseMemObject(energy)
ierr = clReleaseMemObject(seeds)

! and deallocate all arrays 
call mem%dealloc(Lamresx, 'Lamresx')
call mem%dealloc(Lamresy, 'Lamresy')
call mem%dealloc(depthres, 'depthres')
call mem%dealloc(energyres, 'energyres')
call mem%dealloc(rnseeds, 'rnseeds')
call mem%dealloc(init_seeds, 'init_seeds')

end subroutine EMsoftCgetMCOpenCL



    
end module mod_SEMCLwrappers
    
