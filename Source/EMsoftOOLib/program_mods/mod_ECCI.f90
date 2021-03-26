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

module mod_ECCI
  !! author: MDG
  !! version: 1.0
  !! date: 03/14/20
  !!
  !! class definition for the EMECCI program

use mod_kinds
use mod_global

IMPLICIT NONE


! namelist for the EMECP program
type, public :: ECCINameListType
  integer(kind=irg)       :: stdout
  integer(kind=irg)       :: nthreads
  integer(kind=irg)       :: k(3)
  integer(kind=irg)       :: nktstep
  integer(kind=irg)       :: DF_npix
  integer(kind=irg)       :: DF_npiy
  real(kind=sgl)          :: voltage
  real(kind=sgl)          :: dkt
  real(kind=sgl)          :: ktmax
  real(kind=sgl)          :: euler(3)
  real(kind=sgl)          :: lauec(2)
  real(kind=sgl)          :: lauec2(2)
  real(kind=sgl)          :: dmin
  real(kind=sgl)          :: DF_L
  real(kind=sgl)          :: DF_slice
  character(4)            :: dispmode
  character(4)            :: summode
  character(5)            :: progmode
  character(4)            :: mode
  character(fnlen)        :: xtalname
  character(fnlen)        :: montagename
  character(fnlen)        :: defectfilename
  character(fnlen)        :: dispfile
  character(fnlen)        :: energyfile
  character(fnlen)        :: dataname
  character(fnlen)        :: ECPname
  character(fnlen)        :: sgname
  character(fnlen)        :: BetheParametersFile
end type ECCINameListType

! class definition
type, public :: ECCI_T
private
  character(fnlen)                :: nmldeffile = 'EMECCI.nml'
  type(ECCINameListType),public    :: nml

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: ECCI_
  procedure, pass(self) :: Calckvectorcircle_
  procedure, pass(self) :: Calckvectorcone_
  procedure, pass(self) :: Calckvectortrace_
  procedure, pass(self) :: getDisplacementField_
  procedure, pass(self) :: CalcSgarray_
  procedure, pass(self) :: CalcExpval_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: readDispfileHDF_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: ECCI => ECCI_
  generic, public :: Calckvectorcircle => Calckvectorcircle_
  generic, public :: Calckvectorcone => Calckvectorcone_
  generic, public :: Calckvectortrace => Calckvectortrace_
  generic, public :: getDisplacementField => getDisplacementField_
  generic, public :: CalcSgarray => CalcSgarray_
  generic, public :: CalcExpval => CalcExpval_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readDispfileHDF => readDispfileHDF_

end type ECCI_T

! the constructor routine for this class
interface ECCI_T
  module procedure ECCI_constructor
end interface ECCI_T

contains

!--------------------------------------------------------------------------
type(ECCI_T) function ECCI_constructor( nmlfile ) result(ECCI)
!DEC$ ATTRIBUTES DLLEXPORT :: ECCI_constructor
!! author: MDG
!! version: 1.0
!! date: 03/14/20
!!
!! constructor for the ECCI_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

if (present(nmlfile)) call ECCI%readNameList(nmlfile)

end function ECCI_constructor

!--------------------------------------------------------------------------
subroutine ECCI_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: ECCI_destructor
!! author: MDG
!! version: 1.0
!! date: 03/14/20
!!
!! destructor for the ECCI_T Class

IMPLICIT NONE

type(ECCI_T), INTENT(INOUT)  :: self

call reportDestructor('ECCI_T')

end subroutine ECCI_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 03/14/20
!!
!! read the namelist from an nml file for the ECCI_T Class

use mod_io
use mod_EMsoft

IMPLICIT NONE

class(ECCI_T), INTENT(INOUT)      :: self
character(fnlen),INTENT(IN)      :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)      :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                   :: EMsoft
type(IO_T)                       :: Message
logical                          :: skipread = .FALSE.

! parameters for the standard ECCI program
integer(kind=irg)       :: stdout
integer(kind=irg)       :: nthreads
integer(kind=irg)       :: k(3)
integer(kind=irg)       :: nktstep
integer(kind=irg)       :: DF_npix
integer(kind=irg)       :: DF_npiy
real(kind=sgl)          :: voltage
real(kind=sgl)          :: dkt
real(kind=sgl)          :: ktmax
real(kind=sgl)          :: euler(3)
real(kind=sgl)          :: lauec(2)
real(kind=sgl)          :: lauec2(2)
real(kind=sgl)          :: dmin
real(kind=sgl)          :: DF_L
real(kind=sgl)          :: DF_slice
character(4)            :: dispmode
character(4)            :: summode
character(5)            :: progmode
character(4)            :: mode
character(fnlen)        :: xtalname
character(fnlen)        :: montagename
character(fnlen)        :: defectfilename
character(fnlen)        :: dispfile
character(fnlen)        :: energyfile
character(fnlen)        :: dataname
character(fnlen)        :: ECPname
character(fnlen)        :: sgname
character(fnlen)        :: BetheParametersFile

namelist / ECCIlist / DF_L, DF_npix, DF_npiy, DF_slice, dmin, sgname, stdout, &
                      progmode, dispfile, ktmax, dkt, ECPname, summode, lauec, lauec2, &
                      dispmode, mode, nthreads, xtalname, voltage, k, nktstep, &
                      dataname, BetheParametersFile, defectfilename, montagename, energyfile, euler

! set the input parameters to default values (except for xtalname, which must be present)
stdout = 6
nthreads = 1
k = (/ 0,0,1 /)
nktstep = 20
DF_npix = 256
DF_npiy = 256
voltage = 30.
dkt = 0.1
ktmax = 5.0
euler = (/ 0.0, 0.0, 0.0 /)
lauec = (/ 0.0, 0.0 /)
lauec2 = (/ 0.0, 0.0 /)
dmin = 0.1
DF_L = 1.0
DF_slice = 1.0
dispmode = 'not'
summode = 'diag'
progmode = 'array'
mode = 'full'
xtalname = 'undefined'
montagename = 'undefined'
defectfilename = 'undefined'
dispfile = 'displacements.data'
energyfile = 'undefined'
dataname = 'ECCIout.data'
ECPname = 'undefined'
sgname = 'nofile'
BetheParametersFile = 'undefined'


! read the name list, depending on the class type
if (.not.skipread) then
! read the namelist file
  open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
  read(UNIT=dataunit,NML=ECCIlist)
  close(UNIT=dataunit,STATUS='keep')

! check for required entries
  if (trim(xtalname).eq.'undefined') then
    call Message%printError('EMECCI:',' crystal file name is undefined in '//nmlfile)
  end if

  ! make sure the ECPname variable has been properly defined
  if (trim(ECPname).eq.'undefined') then
    call Message%printError('EMECCI:',' ECP pattern file name is undefined in '//nmlfile)
  end if
end if


! if we get here, then all appears to be ok, and we need to fill in the emnl fields
self%nml%stdout = stdout
self%nml%nthreads = nthreads
self%nml%k = k
self%nml%nktstep = nktstep
self%nml%DF_npix = DF_npix
self%nml%DF_npiy = DF_npiy
self%nml%voltage = voltage
self%nml%dkt = dkt
self%nml%ktmax = ktmax
self%nml%euler = euler
self%nml%lauec = lauec
self%nml%lauec2 = lauec2
self%nml%dmin = dmin
self%nml%DF_L = DF_L
self%nml%DF_slice = DF_slice
self%nml%dispmode = dispmode
self%nml%summode = summode
self%nml%mode = mode
self%nml%progmode = progmode
self%nml%xtalname = xtalname
self%nml%montagename = montagename
self%nml%defectfilename = defectfilename
self%nml%dispfile = dispfile
self%nml%energyfile = energyfile
self%nml%dataname = dataname
self%nml%ECPname = ECPname
self%nml%sgname = sgname
self%nml%BetheParametersFile = BetheParametersFile

end subroutine readNameList_


!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 03/14/20
!!
!! pass the namelist for the ECCI_T Class to the calling program

IMPLICIT NONE

class(ECCI_T), INTENT(INOUT)   :: self
type(ECCINameListType)         :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG
!! version: 1.0
!! date: 03/14/20
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants
use ISO_C_BINDING

IMPLICIT NONE

class(ECCI_T), INTENT(INOUT)            :: self
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 5, n_real = 6
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( emnl => self%nml )

! create the group for this namelist
groupname = trim(HDFnames%get_NMLlist())
hdferr = HDF%createGroup(groupname)

! write all the single integers
io_int = (/ emnl%stdout, emnl%nthreads, emnl%nktstep, emnl%DF_npix, emnl%DF_npiy /)
intlist(1) = 'stdout'
intlist(2) = 'nthreads'
intlist(3) = 'nktstep'
intlist(4) = 'DF_npix'
intlist(5) = 'DF_npiy'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! vectors
dataset = SC_k
hdferr = HDF%writeDatasetIntegerArray(dataset, emnl%k, 3)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECCINameList: unable to create k dataset', hdferr)


! write all the single reals
io_real = (/ emnl%voltage, emnl%dkt, emnl%ktmax, emnl%dmin, emnl%DF_L, emnl%DF_slice /)
reallist(1) = 'voltage'
reallist(2) = 'dkt'
reallist(3) = 'ktmax'
reallist(4) = 'dmin'
reallist(5) = 'DF_L'
reallist(6) = 'DF_slice'
call HDF%writeNMLreals(io_real, reallist, n_real)

! 2-vectors
dataset = SC_lauec
hdferr = HDF%writeDatasetFloatArray(dataset, emnl%lauec, 2)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECCINameList: unable to create lauec dataset', hdferr)

dataset = SC_lauec2
hdferr = HDF%writeDatasetFloatArray(dataset, emnl%lauec2, 2)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECCINameList: unable to create lauec2 dataset', hdferr)

! write all the strings
dataset = SC_dispmode
line2(1) = emnl%dispmode
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECCINameList: unable to create dispmode dataset', hdferr)

dataset = SC_summode
line2(1) = emnl%summode
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECCINameList: unable to create summode dataset', hdferr)

dataset = SC_progmode
line2(1) = emnl%progmode
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECCINameList: unable to create progmode dataset', hdferr)

dataset = SC_xtalname
line2(1) = emnl%xtalname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECCINameList: unable to create xtalname dataset', hdferr)

dataset = SC_montagename
line2(1) = emnl%montagename
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECCINameList: unable to create montagename dataset', hdferr)

dataset = SC_defectfilename
line2(1) = emnl%defectfilename
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECCINameList: unable to create defectfilename dataset', hdferr)

dataset = SC_dispfile
line2(1) = emnl%dispfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECCINameList: unable to create dispfile dataset', hdferr)

dataset = SC_dataname
line2(1) = emnl%dataname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECCINameList: unable to create dataname dataset', hdferr)

dataset = SC_ECPname
line2(1) = emnl%ECPname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECCINameList: unable to create ECPname dataset', hdferr)

dataset = SC_sgname
line2(1) = emnl%sgname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteECCINameList: unable to create sgname dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
recursive subroutine readDispfileHDF_(self, ipar, dispfield, fname)
!DEC$ ATTRIBUTES DLLEXPORT :: readDispfileHDF_
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use HDF5
use stringconstants 
use mod_io
use mod_rotations

use ISO_C_BINDING

IMPLICIT NONE

class(ECCI_T), INTENT(INOUT)           :: self 

type(HDF_T)                                 :: HDF
type(HDFnames_T)                            :: HDFnames
type(io_T)                                  :: Message
character(fnlen),INTENT(IN)                 :: fname
integer(kind=irg),INTENT(INOUT)             :: ipar(3)
real(kind=dbl),allocatable,INTENT(INOUT)    :: dispfield(:,:,:,:)
integer(kind=irg)                           :: io_int(1), i, hdferr

real(kind=sgl),parameter                    :: dtor = 0.0174533  ! convert from degrees to radians
integer(kind=irg)                           :: istat
character(fnlen)                            :: deformationfile, groupname, dataset
logical                                     :: g_exists 

real(kind=dbl),allocatable                  :: pcxy(:,:), x(:,:)
integer(HSIZE_T)                            :: dims1(1), dims2(2), dims3(3), dims4(4)

call openFortranHDFInterface()
HDF = HDF_T()

deformationfile = fname

hdferr =  HDF%openFile(deformationfile)

! write the EMheader to the file
groupname = 'EMdata'
hdferr = HDF%openGroup(groupname)

! read the single integer parameters
ipar = 0
dataset = 'npx'  ! ipar(1)

call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if(g_exists) call HDF%readDatasetInteger(dataset, hdferr, ipar(1))

dataset = 'npy'  ! ipar(2)

call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if(g_exists) call HDF%readDatasetInteger(dataset, hdferr, ipar(2))

dataset = 'npz'  ! ipar(3)

call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if(g_exists) call HDF%readDatasetInteger(dataset, hdferr, ipar(3)) 

allocate(dispfield(ipar(3),ipar(1),ipar(2),3))

  ! and finally, read the deformation field dataset
dataset = 'disparray'
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if(g_exists) call HDF%readDatasetDoubleArray(dataset, dims4, hdferr, dispfield) 

  ! close the group and file
call HDF%pop(.TRUE.)

call Message%printMessage('')
call Message%printMessage(' -> completed reading displacement field info from file '//trim(deformationfile))
call Message%printMessage('')

end subroutine readDispfileHDF_

!--------------------------------------------------------------------------
subroutine ECCI_(self, EMsoft, progname, nmldeffile)
!DEC$ ATTRIBUTES DLLEXPORT :: ECCI_
!! author: MDG
!! version: 1.0
!! date: 03/14/20
!!
!! perform the computations

use mod_EMsoft
use mod_initializers
use mod_symmetry
use mod_crystallography
use mod_gvectors
use mod_kvectors
use mod_io
use mod_math
use mod_diffraction
use mod_timing
use mod_Lambert
use HDF5
use mod_HDFsupport
use mod_HDFnames
use ISO_C_BINDING
use omp_lib
use mod_OMPsupport
use mod_notifications
use stringconstants
use mod_MCfiles
use mod_defect
use mod_image
use mod_rotations
use mod_memory

use, intrinsic :: iso_fortran_env

IMPLICIT NONE

class(ECCI_T), INTENT(INOUT)             :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname
character(fnlen), INTENT(INOUT)         :: nmldeffile

type(HDF_T)                             :: HDF
type(HDFnames_T)                        :: HDFnames
type(Cell_T)                            :: cell
type(DynType)                           :: Dyn
type(SpaceGroup_T)                      :: SG
type(Diffraction_T)                     :: Diff
type(IO_T)                              :: Message
type(Defect_T)                          :: Defects
type(gvectors_T)                        :: reflist
type(kvectors_T)                        :: kvec
type(reflisttype),pointer               :: firstw, rltmp, rltmpa, rltmpb
type(image_t)                           :: im
type(Timing_T)                          :: timer
type(MCfile_T)                          :: MCFT
type(e_T)                               :: eu
type(q_t)                               :: qu
type(Quaternion_T)                      :: quat
type(memory_T)                          :: mem 

type(MCOpenCLNameListType) :: mcnl

integer(kind=irg)       :: numangle, numzbins, nx, ny, npx, npy, totnum_el, numsites ! reading from MC file
real(kind=dbl)          :: EkeV, Ehistmin, Ebinsize, depthmax, depthstep, sig, omega  ! reading from MC file
integer(kind=irg), allocatable :: acc_z(:,:,:,:),accum_z(:,:,:,:) ! reading from MC file
integer(kind=irg)       ::num_el, etotal, nsx, nsy, izz, iz

integer(kind=irg)       :: pgnum, SamplingType, nns, nnw, ipar(3)
real(kind=dbl),allocatable      :: dispfield(:,:,:,:)

integer(kind=irg)                       :: kkk, nn,i,j,npix,npiy,ii,jj, numset, t_interval,nat(maxpasym), montage_nx, montage_ny, &
                                           DF_nums_new,DF_npix_new,DF_npiy_new, numstart,numstop, isg, TID, tickstart, &
                                           NTHR, SETNTHR, isym, ir, ga(3), gb(3),ic,g,numd,ix,iy,nkt,nbeams, ik, ig, &
                                           k, numk,ixp,iyp, io_int(6), skip, gg(3), error_cnt, dinfo, nref, maxXY, hdferr, ijmax

real(kind=dbl)                          :: io_real(5), lambda ! real output variable

integer(kind=irg),parameter             :: numdd=360 ! 180

real(kind=sgl)                          :: thetacr,  nabsl, thick, X(2), bragg, thetac, kstar(3), gperp(3), av, mi, ma, &
                                           gdotR,DF_gf(3), tpi, DM(2,2), DD, c(3), gx(3), gy(3), &
                                           gac(3), gbc(3),zmax, tstop, kk(3), FN(3), gc(3) 
real(kind=dbl)                          ::  kt(3), delta , ktmax, arg, glen, DynFN(3), xx
complex(kind=dbl),allocatable           :: amp(:),amp2(:),Azz(:,:),DF_R(:,:)
real(kind=sgl),allocatable              :: lambdaZ(:), disparray(:,:,:,:),imatvals(:,:), ECCIimages(:,:,:), XYarray(:,:), &
                                           ECCIstore(:,:,:)
real(kind=sgl),allocatable              :: svals(:), sgarray(:), klist(:,:), knlist(:)
integer(kind=sgl),allocatable           :: expval(:,:,:),  hklarray(:,:), nab(:,:), XYint(:,:)
complex(kind=dbl)                       :: para(0:numdd),dx,dy,dxm,dym, xgp
complex(kind=dbl),allocatable           :: DHWM(:,:),DHWMvoid(:,:),DDD(:,:),Sarray(:,:,:,:), Sarrayk(:,:,:,:,:)
complex(kind=dbl),allocatable           :: Lgh(:), Lgh2(:,:), Sgh(:), Sghtmp(:), DHWMz(:,:), Sgh2(:,:), Sghtmp2(:,:)
complex(kind=dbl)                       :: czero=cmplx(0.D0,0.D0),cone=dcmplx(1.D0,0.D0)
type(kvectorlist),pointer               :: khead, ktmp
logical                                 :: verbose

integer(kind=irg),allocatable           :: kij(:,:)

character(fnlen)                        :: groupname, dataset, datagroupname,  outname, Image_filename, fname, energyfile

type(gnode)                             :: rlp
type(BetheParameterType)                :: BetheParameters

character(11)                           :: dstr
character(15)                           :: tstrb
character(15)                           :: tstre

integer(HSIZE_T)                        :: cnt3(3), offset3(3),cnt4(4), offset4(4)

integer(int8), allocatable              :: montage(:,:)

integer                                 :: iostat
character(len=128)                      :: iomsg

logical                                 :: overwrite = .TRUE., insert = .TRUE., g_exists, usehex

integer(HSIZE_T), dimension(1:3)        :: hdims, offset
integer(HSIZE_T), dimension(1:5)        :: hdims5, offset5
integer(HSIZE_T)                        :: dims3(3), dims4(4), dims5(5)
character(fnlen,kind=c_char)            :: line2(1)

character (len=10) :: file_name
character(len=10) :: file_id

call openFortranHDFInterface()
HDF = HDF_T()

! set the HDF group names for this program
HDFnames = HDFnames_T()

! simplify the notation a little
associate( emnl => self%nml )

! initialize the timing routines
timer = Timing_T()
tstrb = timer%getTimeString()

!=============================================
!=============================================
! ---------- read Monte Carlo .h5 output file and extract necessary parameters
call HDFnames%set_ProgramData(SC_MCOpenCL)
call HDFnames%set_NMLlist(SC_MCCLNameList)
call HDFnames%set_NMLfilename(SC_MCOpenCLNML)
energyfile = ''
energyfile = EMsoft%generateFilePath('EMdatapathname',trim(emnl%energyfile))
outname = trim(energyfile)
call MCFT%setFileName(energyfile)
call MCFT%readMCfile(HDF, HDFnames, getAccumz=.TRUE.)
mcnl = MCFT%getnml()
call MCFT%copyaccumz(accum_z)

numzbins = MCFT%getnumzbins()
numangle = MCFT%getnumangles()

nsx = (mcnl%numsx - 1)/2
nsy = nsx
etotal = sum(accum_z(numangle,:,:,:))

io_int(1) = mcnl%totnum_el
call Message%WriteValue(' --> Total number of BSE electrons in MC data set ', io_int, 1)

!=============================================
!=============================================
! crystallography section
verbose = .TRUE.

call cell%setFileName(emnl%xtalname)

call Diff%setrlpmethod('WK')
call Diff%setV(dble(emnl%voltage))

call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, emnl%dmin, verbose, useHDF=HDF)
call Diff%Printrlp()
lambda = Diff%getWaveLength()

! check the crystal system and setting; abort the program for trigonal with rhombohedral setting with
! an explanation for the user
if ((SG%getSpaceGroupXtalSystem().eq.5).and.(cell%getLatParm('b').eq.cell%getLatParm('c'))) then
    call Message%printMessage( (/ &
    '                                                                         ', &
    ' ========Program Aborted========                                         ', &
    ' The ECP master pattern simulation for rhombohedral/trigonal structures  ', &
    ' requires that the structure be described using the hexagonal reference  ', &
    ' frame.  Please re-enter the crystal structure in this setting and re-run', &
    ' the Monte Carlo calculation and this master pattern program.            '/) )
    stop
end if

! allocate and compute the Sgh loop-up table
numset = cell%getNatomtype()
call Diff%Initialize_SghLUT(cell, SG, emnl%dmin, numset, nat, verbose)

write(*,*) sum(nat(1:numset))
 
! determine the point group number
j=0
do i=1,32
    if (SGPG(i).le.SG%getSpaceGroupNumber()) j=i
end do
isym = j

! here is new code dealing with all the special cases (quite a few more compared to the
! Laue group case)...  isym is the point group number. Once the symmetry case has been
! fully determined (taking into account things like 31m and 3m1 an such), then the only places
! that symmetry is handled are the modified Calckvectors routine, and the filling of the modified
! Lambert projections after the dynamical simulation step.  We are also changing the name of the
! sr array (or srhex) to mLPNH and mLPSH (modified Lambert Projection Northern/Southern Hemisphere),
! and we change the output HDF5 file a little as well. We need to make sure that the EMEBSD program
! issues a warning when an old format HDF5 file is read.

! Here, we encode isym into a new number that describes the sampling scheme; the new schemes are
! described in detail in the EBSD manual pdf file.

SamplingType = PGSamplingType(isym)

! next, intercept the special cases (hexagonal vs. rhombohedral cases that require special treatment)
if ((SamplingType.eq.-1).or.(isym.eq.14).or.(isym.eq.26)) then
  SamplingType = SG%getHexvsRho(isym)
end if

! if the point group is trigonal or hexagonal, we need to switch usehex to .TRUE. so that
! the program will use the hexagonal sampling method
usehex = .FALSE.
if ((SG%getSpaceGroupXtalSystem().eq.4).or.(SG%getSpaceGroupXtalSystem().eq.5)) usehex = .TRUE.

! ---------- end of symmetry and crystallography section
!=============================================
!=============================================

! set the Bethe parameters
! force dynamical matrix routine to read new Bethe parameters from file
call Diff%SetBetheParameters(EMsoft, .TRUE., emnl%BetheParametersFile)

!=============================================
! ---------- create the incident beam directions list
! determine all independent incident beam directions (use a linked list starting at khead)
! numk is the total number of k-vectors to be included in this computation;
! note that this needs to be redone for each energy, since the wave vector changes with energy
kvec = kvectors_T()   ! initialize the wave vector list

!call eu%e_setd( emnl%euler*dtor )
!qu = eu%eq()
!call qu%q_print('Quaternion            : ')
!quat = Quaternion_T( qd = qu%q_copyd() )

call SG%BFsymmetry(emnl%k,j,isym,ir)

! determine and display the shortest reciprocal lattice vectors for this zone
call cell%ShortestG(SG,emnl%k,ga,gb,isym)

! here we figure out how many beams there are
if (trim(emnl%progmode).eq.'array') then
  call kvec%set_kinp( dble(emnl%k) )
  call kvec%set_ktmax( dble(emnl%ktmax) )
  call kvec%set_SamplingType( SamplingType )

  npx = nint(emnl%ktmax/emnl%dkt)
  npy = npx
  ijmax = float(npx)**2  

  call kvec%set_mapmode('ECCI')
  if (usehex) then
  call kvec%Calckvectors(cell, SG, Diff, dble(ga),npx,npy, ijmax,usehex)
  else
  call kvec%Calckvectors(cell, SG, Diff, dble(ga),npx,npy, ijmax,usehex)
  end if
  numk = kvec%get_numk()
end if 

if (trim(emnl%progmode).eq.'circl') then
  ! get the maximum beam tilt angle
  nkt = nint( emnl%ktmax / emnl%dkt )
  ! Circular trajectory mode
  call self%Calckvectorcircle(cell,lambda,khead,emnl%k,ga,emnl%lauec,emnl%ktmax,nkt,numk)
  io_int(1)=numk
end if

io_int(1)=numk
call Message%WriteValue('# independent beam directions to be considered = ', io_int, 1, "(I8)")

! Convert to array for OpenMP
mem = memory_T()

call mem%alloc( kij, (/3, numk/), 'kij')
call mem%alloc( klist, (/3, numk/), 'klist')
call mem%alloc( knlist, (/numk/), 'knlist')
call mem%alloc( XYarray, (/2, numk/), 'XYarray')

! ! point to the first beam direction
if (trim(emnl%progmode).eq.'array') ktmp => kvec%get_ListHead()
if (trim(emnl%progmode).eq.'circl') ktmp => khead

! ! and loop through the list, keeping k, kn, and i,j
kij(1:3,1) = (/ ktmp%i, ktmp%j, ktmp%hs /)
klist(1:3,1) = ktmp%k
knlist(1) = ktmp%kn

do i = 2,numk
  ktmp => ktmp%next
  kij(1:3,i) = (/ ktmp%i, ktmp%j, ktmp%hs /)
  klist(1:3,i) = ktmp%k
  knlist(i) = ktmp%kn
end do

! number of reflections, and associated information (hkl, ...)
call cell%TransSpace(float(emnl%k),c,'d','c')
call cell%NormVec(c,'c')
! then make ga the x-axis
call cell%TransSpace(float(ga),gx,'r','c')
call cell%NormVec(gx,'c')
! compute the cross product between k and gx; this is the y-axis
call cell%CalcCross(c,gx,gy,'c','c',0)

! this needs to be fixed !!!!!
call cell%TransSpace(float(emnl%k),kstar,'d','r')        ! transform incident direction to reciprocal space
call cell%CalcCross(float(ga),kstar,gperp,'r','r',0)! compute g_perp = ga x k
call cell%NormVec(gperp,'r')                        ! normalize g_perp

glen = cell%CalcLength(float(ga),'r')

if (trim(emnl%progmode).eq.'array') ktmp => kvec%get_ListHead()
if (trim(emnl%progmode).eq.'circl') ktmp => khead 

do ic=1,numk
    XYarray(1,ic) = -cell%CalcDot(sngl(ktmp%kt),float(ga),'c') ! / glen
    XYarray(2,ic) = -cell%CalcDot(sngl(ktmp%kt),gperp,'c') * glen
    ktmp => ktmp%next
end do

call kvec%Delete_kvectorlist()
nullify(ktmp)

!=============================================
! ---------- Compute lambdaZ

call Diff%setrlpmethod('WK')
call Diff%CalcUcg(cell,(/0,0,0/))
rlp = Diff%getrlp()
nabsl = rlp%xgp    
io_real(1) = nabsl

call Message%WriteValue('Normal absorption length : ', io_real, 1, "(F10.5/)")

call mem%alloc(lambdaZ, (/numzbins/), 'lambdaZ', 0.0_sgl)

do iz=1,numzbins
  lambdaZ(iz) = float(sum(accum_z(numangle,iz,:,:)))/float(etotal)
  lambdaZ(iz) = lambdaZ(iz) * exp(2.0*sngl(cPi)*(iz-1)*mcnl%depthstep/nabsl)
end do

!=============================================
!=============================================
! Defects section

! copy some of the namelist parameters into the defects structure
Defects%DF_npix = emnl%DF_npix
Defects%DF_npiy = emnl%DF_npiy
Defects%DF_L = emnl%DF_L
Defects%DF_slice = emnl%DF_slice

! If there is no displacement field file we compute displacement field
if (emnl%dispfile.eq.'undefined') then

  call Message%printMessage(' --> Compute displacement field')
  
  ! determine the cartesian components of ga
  Defects%DF_gf = float(ga)
  Defects%DF_gstar = Defects%DF_gf/cell%CalcLength(defects%DF_gf,'r')**2    ! define G* such that G.G* = 1
  call cell%TransSpace(Defects%DF_gf,Defects%DF_gc,'r','c')         ! convert to Cartesian reference frame

  ! next, we read all the foil and defect data using the new InitializeDefects routine in defectmodule.f90
  verbose = .FALSE.
  call Defects%InitializeDefects(EMsoft,cell,emnl%defectfilename,emnl%DF_npix,emnl%DF_npiy,emnl%DF_L,DF_gf,error_cnt,verbose)

  DynFN = Defects%foil%F  

  ! define the foil thickness, attenuation, and number slices per column
  thick = Defects%foil%zb    ! this is the same everywhere for this version; needs to be updated in the next version
  Defects%DF_nums = nint(thick/Defects%DF_slice)  ! this is the number of slices for this particular column

  allocate(Defects%DF_R(Defects%DF_nums,3))     ! each thread has its own DF_R array

  allocate(disparray(2,Defects%DF_nums,Defects%DF_npix,Defects%DF_npiy))
  disparray = 0.0

  call self%getDisplacementField(cell, Defects, ga, gb, disparray)

else

  ! read the displacement field arrays from the displacement HDF file
  !=================================================================
  fname = EMsoft%generateFilePath('EMdatapathname',trim(self%nml%dispfile))
  call self%readDispfileHDF(ipar, dispfield, fname)

  ! define the foil thickness, attenuation, and number slices per column
  DynFN = (/0,0,1/)
  thick = ipar(3)
  Defects%DF_nums = nint(thick/Defects%DF_slice)  ! this is the number of slices for this particular column
  
  allocate(imatvals(2,Defects%DF_nums))
  imatvals = 0

  call cell%TransSpace(float(ga),gac,'r','c')
  call cell%TransSpace(float(gb),gbc,'r','c')

  numd = 360

  allocate(disparray(2,Defects%DF_nums,ipar(1),ipar(2)))
  disparray = 0.0

  do i=1,ipar(1)
    do j=1,ipar(2)
   ! loop over the fixed thickness slices
      do ik=1,Defects%DF_nums
        gdotR = sum(gac*dispfield(ik,i,j,1:3))
        if (NANCHK(gdotR)) gdotR = 0.0
        xx = numd*amod(gdotR+1000.0,1.0)
        imatvals(1,ik) = xx
        gdotR = sum(gbc*dispfield(ik,i,j,1:3))
        if (NANCHK(gdotR)) gdotR = 0.0
        xx = numd*amod(gdotR+1000.0,1.0)
        imatvals(2,ik) = xx
      end do
      disparray(1:2,1:Defects%DF_nums,i,j) = imatvals(1:2,1:Defects%DF_nums)
    end do
  end do

  io_real(1) = minval(disparray)
  io_real(2) = maxval(disparray)
  call Message%WriteValue(' --> Disparray bounds: ', io_real, 2, "(2(F10.5,' '))")

end if

! ok, all the set up is now complete;
npix = Defects%DF_npix
npiy = Defects%DF_npiy
allocate(svals(defects%DF_nums))
allocate(ECCIimages(npix,npiy,1))
ECCIimages = 0.0

!if (trim(emnl%montagename).ne.'undefined') then
  allocate(ECCIstore(npix,npiy,numk))
  ECCIstore = 0.0
!end if

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! store necessary data in data file
  call Message%printMessage(' --> Storing output ECCIs in '//trim(emnl%dataname),"(A/)")

  outname = EMsoft%generateFilePath('EMdatapathname',trim(emnl%dataname))
!=============================================
! create or update the HDF5 output file
!=============================================
  call HDFnames%set_ProgramData(SC_ECCI)
  call HDFnames%set_NMLlist(SC_ECCINameList)
  call HDFnames%set_NMLfilename(SC_ECCImasterNML)

! Open an existing file or create a new file using the default properties.
  hdferr =  HDF%createFile(outname)

! write the EMheader to the file
  datagroupname = trim(HDFnames%get_ProgramData())
  call HDF%writeEMheader(EMsoft,dstr, tstrb, tstre, progname, datagroupname)

  ! create a namelist group to write all the namelist files into
  groupname = SC_NMLfiles
    hdferr = HDF%createGroup(groupname)

  ! read the text file and write the array to the file
  dataset = SC_ECCImasterNML
    hdferr = HDF%writeDatasetTextFile(dataset, nmldeffile)

  ! in this case, we also put the defect json file inside the HDF file
  dataset = SC_ECCIdefectJSON
    outname = EMsoft%generateFilePath('EMdatapathname',trim(emnl%defectfilename))
    hdferr = HDF%writeDatasetTextFile(dataset, outname)
! If there is no displacement field file we compute displacement field
if (emnl%dispfile.eq.'undefined') then
  ! and also the foil descriptor file
  dataset = SC_ECCIfoilJSON
    outname = EMsoft%generateFilePath('EMdatapathname',trim(Defects%foilname))
    hdferr = HDF%writeDatasetTextFile(dataset, outname)
end if
 ! leave this group
  call HDF%pop()

  ! create a namelist group to write all the namelist files into
  hdferr = HDF%createGroup(HDFnames%get_NMLparameters())
  call self%writeHDFNameList(HDF, HDFnames)

  ! we do not need to parse the defect information, since the defect and foil JSON files
  ! are both included in the NMLfiles section of the HDF5 output file., and we have routines in
  ! JSONsupport to parse those files.

  ! leave this group
  call HDF%pop()

  ! then the remainder of the data in a EMData group
  groupname = SC_EMData
    hdferr = HDF%createGroup(groupname)

  dataset = SC_Braggga
    hdferr = HDF%writeDatasetFloat(dataset, Diff%CalcDiffAngle(cell,(/ga(1),ga(2),ga(3)/))*0.5)

  dataset = SC_thetac
    hdferr = HDF%writeDatasetFloat(dataset, thetac)

  dataset = SC_nref
    hdferr = HDF%writeDatasetInteger(dataset, nn)

  dataset = SC_numk
    hdferr = HDF%writeDatasetInteger(dataset, numk)

  dataset = SC_npix
    hdferr = HDF%writeDatasetInteger(dataset, emnl%DF_npix)

  dataset = SC_npiy
    hdferr = HDF%writeDatasetInteger(dataset, emnl%DF_npiy)

! create the ECCIimage hyperslab and write zeroes to them for now

    dataset ='Disparray'
    dims4 = (/  2, Defects%DF_nums, Defects%DF_npix, Defects%DF_npiy /)
    cnt4 = (/ 2,Defects%DF_nums,Defects%DF_npix,Defects%DF_npiy /)
    offset4 = (/ 0, 0, 0,0 /)
    hdferr = HDF%writeHyperslabFloatArray(dataset, disparray, dims4, offset4, cnt4)

  ! create the ECCIimage hyperslab and write zeroes to them for now
  dataset = SC_ECCIimages
    dims3 = (/  npix, npiy, numk /)
    cnt3 = (/ npix, npiy, 1 /)
    offset3 = (/ 0, 0, 0 /)
    hdferr = HDF%writeHyperslabFloatArray(dataset, ECCIimages, dims3, offset3, cnt3)

    call HDF%pop(.TRUE.)
    ! close the Fortran interface
    call closeFortranHDFInterface()

! leave the HDF5 file open for final output at the end of the program
!--------------------------------------------------------------
!--------------------------------------------------------------
!--------------------------------------------------------------


!=============================================
! ---------- Main Loop Initialization

numstart = 1
numstop = numk
io_int(1) = numk
call Message%WriteValue(' --> ECCI: number of beam directions =  ', io_int, 1, "(I5)")

!DynFN = (/0.0,0.0,1.0/)

! define the numd complex defect parameters (i.e., precompute the sin and cos arrays
numd = 360
para = czero
do i=0,numd
  arg = 2.D0*cPi*dble(i)/dble(numd)
  para(i) = cmplx(dcos(arg),-dsin(arg))
end do
!rite(*,*) para
cone = cmplx(1.D0,0.D0)

mainloop: do isg = numstart,numstop 
  ECCIimages = 0.0

  !=============================================
! ---------- create the master reflection list for this beam direction
  reflist = gvectors_T()

  kk = klist(1:3,isg)
  call cell%NormVec(kk,'c')
  kk = kk/lambda
  call cell%TransSpace(kk,kstar,'c','r')
  kk = kstar
  FN = kstar
  
  call reflist%Initialize_ReflectionList(cell, SG, Diff, sngl(DynFN), kk, sngl(emnl%dmin), verbose)
  nn = reflist%get_nref()

  ! ---------- end of "create the master reflection list"
  !=============================================
  
  nullify(firstw)
  
  ! There is now two modes available: 'full' which is the standard full dynmical mode (slow)
  ! and 'fast' that use Bethe Parameters

  if (self%nml%mode.eq.'fast') then
    nns = 0
    nnw = 0
    call reflist%Apply_BethePotentials(Diff, firstw, nns, nnw)
    nn = nns
  end if

  ! allocate the various DHW Matrices
  !call mem%alloc(DHWMz, (/nn,nn/), 'DHWMz')
  if (allocated(DHWMz)) deallocate(DHWMz)
  allocate(DHWMz(nn,nn))
  DHWMz = czero
  
  !call mem%alloc(DHWMvoid, (/nn,nn/), 'DHWMvoid')
  if (allocated(DHWMvoid)) deallocate(DHWMvoid)
  allocate(DHWMvoid(nn,nn))

  DM = 0.0
  DD = 0.0
  DHWMvoid = czero; DHWMz=czero
  DM(1,1) = cell%CalcDot(float(gb),float(gb),'c')
  DM(1,2) = -cell%CalcDot(float(ga),float(gb),'c')
  DM(2,1) = DM(1,2)
  DM(2,2) = cell%CalcDot(float(ga),float(ga),'c')
  DD = DM(1,1)*DM(2,2) - DM(1,2)*DM(2,1)


  if (self%nml%mode.eq.'Full') then
      nullify(rltmpa)
      nullify(rltmpb)
      call reflist%GetDynMatDHW(cell, Diff, firstw, rltmpa, rltmpb, DHWMz, nn, DM, ga, gb)

  else if (self%nml%mode.eq.'fast') then
      call reflist%GetDynMat(cell, Diff, firstw, DHWMz, nns, nnw)
      ! Conversion from Bloch wave matrix to scattering matrix formalism
      DHWMz = DHWMz * dcmplx(0.D0,cPi * lambda)
  end if 
  
  ! loop over all reflections to get the appropriate powers
  !call mem%alloc(expval, (/2,nn,nn/), 'expval')
  if (allocated(expval)) deallocate(expval)
  allocate(expval(2,nn,nn))
  expval = 0.0
  
  call reflist%GetExpval_ECCI(cell, expval, Diff, firstw, nn, DM, ga, gb )

  ! Compute Sgh
  ! Only diagonals terms are computed gives inverted contrast...
  !call mem%alloc(Sgh, (/nn/), 'Sgh')
  !if (allocated(Sgh)) deallocate(Sgh)
  !allocate(Sgh(nn))

  !nat = 0
  !call Diff%preCalcSghECCI(cell, SG, nn, nat, Sgh)

  ! Computation with non diagonals terms of Sgh
  ! then we need to initialize the Sgh arrays
  if (allocated(Sghtmp2)) deallocate(Sghtmp2)
  allocate(Sghtmp2(nn,nn))
  Sghtmp2 = czero
  call reflist%getSghfromLUTsum(Diff,nn,numset,nat,Sghtmp2)

  call Message%printMessage(' --> Done',"(A)")

  if (self%nml%mode.eq.'Full') then

  ! compute the excitation error for the incident beam directions
  !call mem%alloc(sgarray, (/nn/), 'sgarray')
  if (allocated(sgarray)) deallocate(sgarray)
  allocate(sgarray(nn))

  call reflist%GetSgArray_ECCI(cell, sgarray, dble(klist), numk, isg, DynFN, Diff, firstw, nn)
  forall (i=1:nn)
    DHWMz(i,i)= cmplx(0.D0,2.D0*cPi*sgarray(i)) + xgp ! xgp already has i Pi in it.
    DHWMvoid(i,i) = DHWMz(i,i)
  end forall

  else if (self%nml%mode.eq.'fast') then
    forall (i=1:nn)
      DHWMvoid(i,i) = DHWMz(i,i)
    end forall
  end if

  allocate(Sarray(nn,nn,0:numd,0:numd))
  Sarray = czero
  !call mem%alloc(Sarray, (/nn,nn,numd,numd/), 'Sarray')
  NTHR = emnl%nthreads
  !$OMP  PARALLEL NUM_THREADS(NTHR) DEFAULT(SHARED) PRIVATE(TID,i,j,k,ii,jj,ic,ir,g,Azz,DDD,zmax,Sarrayk)
  TID = OMP_GET_THREAD_NUM()

  allocate(Azz(nn,nn), DDD(nn,nn))   ! these are private variables, so each thread must allocate them !

  if (TID.eq.0) then
    io_int(1) = isg
    call Message%WriteValue(' -> ',io_int,1,"(I4,' ')",advance="no")
    call Message%printMessage('starting Sarray computation',"(A)",advance="no")
  end if

  !call mem%alloc(Azz, (/nn,nn/), 'Azz')
  !call mem%alloc(DDD, (/nn,nn/), 'DDD')
  !$OMP DO SCHEDULE(STATIC) 
  do j=0,numd
    do i=0,numd
    ! loop over all reflections in the array DD using the information in expval
    ! ic is the column index
      do ic=1,nn
    ! ir is the row index
        do ir=1,nn
          if (ic.ne.ir) then  ! exclude the diagonal
            DDD(ir,ic) = DHWMz(ir,ic) 
          else
            DDD(ir,ic) = DHWMz(ir,ic) * para(i)**expval(1,ir,ic) * para(j)**expval(2,ir,ic)
          end if
        end do
      end do
      
      call MatrixExponential(DDD, Azz, dble(defects%DF_slice), 'Pade', nn)

      Sarray(1:nn,1:nn,i,j) = Azz(1:nn,1:nn)
    end do
  end do
  !$OMP END DO
  !call mem%dealloc2(Azz, 'Azz')
  !call mem%dealloc2(DDD, 'DDD')
  deallocate(Azz,DDD)
  !$OMP END PARALLEL
  
  call Message%printMessage(' --> Done; scattering matrix computation ',"(A)",advance="no")

  !----------------------------------------------------!
  ! Finally, here it is: the actual image computation  !
  ! This was modified from the regular (S)TEM mode to  !
  ! reflect the depth integration, which requires a    !
  ! summation of the product of Sgh and Lgh.           !
  !----------------------------------------------------!

  NTHR = 12
  !$OMP  PARALLEL NUM_THREADS(NTHR) DEFAULT(SHARED) &
  !$OMP& PRIVATE(TID,i,j,k,ii,Azz,amp,amp2,ix,iy,dx,dy,dxm,dym,ixp,iyp,Lgh,Lgh2,ir,ic,svals)
  TID = OMP_GET_THREAD_NUM()
  !call mem%alloc(Azz, (/nn,nn/), 'Azz')
  !call mem%alloc(amp, (/nn/), 'amp')
  !call mem%alloc(amp2, (/nn/), 'amp2')
  !call mem%alloc(Lgh, (/nn/), 'Lgh')
  allocate(Azz(nn,nn),amp(nn),amp2(nn),Lgh(nn),Lgh2(nn,nn))

  !$OMP DO SCHEDULE (STATIC)
  donpix: do i=1,npix
    donpiy:   do j=1,npiy
      ! initialize the wave function for this pixel with (1.0,0.0) for the incident beam
      Lgh = czero
      Lgh2 = czero
      amp = czero
      amp(1) = cone
      amp2 = czero
      doslices: do k=1,defects%DF_nums    ! loop over the fixed thickness slices
    ! compute the appropriate scattering matrix to propagate with (see section 8.3.3 in the book)
          if (disparray(1,k,i,j).eq.-10000) then  ! this is point inside a void
            Azz = DHWMvoid    ! so we use the void propagator matrix
          else  ! it is not a void
    ! in this version, we use complex bi-variate interpolation to select the appropriate Azz
    ! matrix from the pre-computed Sarray.
            ix = int(disparray(1,k,i,j))
            ixp = ix+1
            if (ix.eq.numd) ixp=1
            iy = int(disparray(2,k,i,j))
            iyp = iy+1
            if (iy.eq.numd) iyp=1
            dx = cmplx(amod(disparray(1,k,i,j),1.0),0.0)
            dy = cmplx(amod(disparray(2,k,i,j),1.0),0.0)
            dxm = cone-dx
            dym = cone-dy
            Azz = dxm*dym*Sarray(1:nn,1:nn,ix,iy)+dx*dym*Sarray(1:nn,1:nn,ixp,iy)+ &
                  dxm*dy*Sarray(1:nn,1:nn,ix,iyp)+dx*dy*Sarray(1:nn,1:nn,ixp,iyp)
          end if
          amp2 = matmul(Azz,amp)

          if (k.eq.1) then
            Lgh2(1:nn,1:nn) = (mcnl%depthstep/mcnl%depthmax)*lambdaZ(k)*spread(amp2(1:nn),dim=2,ncopies=nn)*&
                                spread(conjg(amp2(1:nn)),dim=1,ncopies=nn)
          else
            Lgh2(1:nn,1:nn) = Lgh2(1:nn,1:nn)+(mcnl%depthstep/mcnl%depthmax)*lambdaZ(k)*spread(amp2(1:nn),dim=2,ncopies=nn)*&
                              spread(conjg(amp2(1:nn)),dim=1,ncopies=nn)
          end if
          
          amp = amp2
      end do doslices ! loop over slices

      ! Commented section used to test Sgh non diagonals terms use...

      svals = 0.0
      svals = real(sum(Lgh2(1:nn,1:nn)*Sghtmp2(1:nn,1:nn)))

      svals = svals/float(sum(nat(1:numset)))
      
    ! then we need to multiply Sgh and Lgh, sum, and take the real part which will
    ! produce the desired BSE intensity
      ECCIimages(i,j,1) =  sngl(sum(svals))

    end do donpiy
  end do donpix
  !$OMP END DO
  !call mem%dealloc2(Azz, 'Azz')
  !call mem%dealloc1(amp, 'amp')
  !call mem%dealloc1(amp2, 'amp2')
  !call mem%dealloc1(Lgh, 'Lgh') 
  deallocate(Azz,amp,amp2,Lgh)
  !$OMP END PARALLEL
  deallocate(Sarray)
  !call mem%dealloc4(Sarray, 'Sarray')

  ECCIimages = ECCIimages / float(Defects%DF_nums) !/ float(sum(nat))

  if ((trim(emnl%montagename).ne.'undefined')) then
    do i=1,emnl%DF_npix
      do j=1,emnl%DF_npiy
        ECCIstore(i,j,isg) = ECCIimages(i,j,1)
      end do
    end do
  end if
  
  call reflist%Delete_gvectorlist()

  ! open the HDF interface
  call openFortranHDFInterface()

  outname = EMsoft%generateFilePath('EMdatapathname',trim(emnl%dataname))

  ! Open an existing file
  hdferr = HDF%openFile(outname)

  ! and update the end time
  call timer%makeTimeStamp()
  tstre = timer%getTimeString()

  groupname = trim(HDFnames%get_EMheader())
  hdferr = HDF%openGroup(groupname)

  groupname = trim(HDFnames%get_ProgramData())
  hdferr = HDF%openGroup(groupname)

  ! and update the end time
  call timer%makeTimeStamp()
  tstre = timer%getTimeString()

  ! stop time /EMheader/StopTime 'character'
  dataset = SC_StopTime
  call timer%Time_tock(1)
  tstop = timer%getInterval(1)
  line2(1) = dstr//', '//tstre
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)

  dataset = SC_Duration
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetFloat(dataset, tstop, overwrite)
  else
    hdferr = HDF%writeDatasetFloat(dataset, tstop)
  end if

  ! close all groups
  call HDF%pop()
  call HDF%pop()

  groupname = SC_EMData
    hdferr = HDF%openGroup(groupname)

  ! add data to the hyperslab
  dataset = SC_ECCIimages

  offset = (/ 0, 0, isg-1 /)
  hdims = (/  npix, npiy, numk /)
  dims3 = (/  npix, npiy, 1 /)
  hdferr = HDF%writeHyperslabFloatArray(dataset, ECCIimages, hdims, offset, dims3, insert)

  call HDF%pop(.TRUE.)

  ! close the Fortran interface
  call closeFortranHDFInterface()
200 end do mainloop

if (trim(emnl%montagename).ne.'undefined') then

  call Message%printMessage('Generating image montage')

! output the montage as an image file (tiff, jpeg, or png)
  fname = EMsoft%generateFilePath('EMdatapathname',trim(emnl%montagename))
  Image_filename = trim(fname)


! allocate the montage array
  if ((emnl%progmode.eq.'array').or.(emnl%progmode.eq.'circl')) then
! divide the beam offsets by the dkt step size so that we can turn them into integers
    XYarray = XYarray / emnl%dkt
    allocate(XYint(2,numk))
    XYint = nint(XYarray)
    maxXY = maxval(XYint)
    montage_nx = emnl%DF_npix*(2*maxXY+1)
    montage_ny = emnl%DF_npiy*(2*maxXY+1)
    allocate(montage(montage_nx,montage_ny))
  else ! progmode = 'trace'
    montage_nx = emnl%DF_npix*numk
    montage_ny = emnl%DF_npiy
    allocate(montage(montage_nx,montage_ny))
  end if

! assign the average value of the ECCIimages array to the montage
  av = sum(ECCIstore)/float(emnl%DF_npix)/float(emnl%DF_npiy)/float(numk)
  montage = av

! get the intensity scaling parameters
  ma = maxval(ECCIstore)
  mi = minval(ECCIstore)

! fill and scale the montage array
  do kkk=1,numk
    do j=1,emnl%DF_npiy
      do i=1,emnl%DF_npix
        if ((emnl%progmode.eq.'array').or.(emnl%progmode.eq.'circl'))  then
          ii = emnl%DF_npix * (maxXY + XYint(1,kkk)) + i
          jj = emnl%DF_npiy * (maxXY + XYint(2,kkk)) + j
          montage(ii,jj) = int(255 * (ECCIstore(i,j,kkk)-mi)/(ma-mi))
        else
          ii = emnl%DF_npix * (kkk-1) + i
          jj = j
          montage(ii,jj) = int(255 * (ECCIstore(i,j,kkk)-mi)/(ma-mi))
        end if
      end do
    end do
  end do
  deallocate(ECCIstore)


! set up the image_t structure
  im = image_t(montage)
  if (im%empty()) call Message%printMessage("EMECCI","failed to convert array to image")

! create the file
  call im%write(trim(Image_filename), iostat, iomsg) ! format automatically detected from extension
  if (0.ne.iostat) then
    call Message%printMessage("Failed to write image to file : "//iomsg)
  else
    call Message%printMessage('ECCI image montage written to '//trim(Image_filename))
  end if
  deallocate(montage)

else
  ! Circle mode section in test to export output in an images stack
  ! get the intensity scaling parameters
  ma = maxval(ECCIstore)
  mi = minval(ECCIstore)

  do kkk=1,numk

    ECCIstore(1:emnl%DF_npix,1:emnl%DF_npiy,kkk) = int(255 * (ECCIstore(1:emnl%DF_npix,1:emnl%DF_npiy,kkk)-mi)/(ma-mi))

    write(file_id, '(i0)') kkk
    ! output the montage as an image file (tiff, jpeg, or png)
    file_name = trim(adjustl(file_id)) // '.tiff'
    fname = EMsoft%generateFilePath('EMdatapathname',trim(file_name))
    Image_filename = trim(fname)
    ! set up the image_t structure
    im = image_t(ECCIstore(1:emnl%DF_npix,1:emnl%DF_npiy,:))
    if (im%empty()) call Message%printMessage("EMECCI","failed to convert array to image")

    ! create the file
      call im%write(trim(Image_filename), iostat, iomsg) ! format automatically detected from extension
      if (0.ne.iostat) then
        call Message%printMessage("Failed to write image to file : "//iomsg)
      else
        call Message%printMessage('ECCI image montage written to '//trim(Image_filename))
      end if
  end do
end if

end associate

end subroutine ECCI_

!--------------------------------------------------------------------------
subroutine  CalcExpval_(self, expval,  nab, nn)
  !DEC$ ATTRIBUTES DLLEXPORT ::  CalcExpval_
  !! author: MDG
  !! version: 1.0
  !! date: 03/14/20
  !!
  !! Compute the exponential value
  
  use mod_crystallography
  use mod_diffraction
  use mod_defect
  use mod_io

  class(ECCI_T), INTENT(INOUT)        :: self
  integer(kind=irg),INTENT(INOUT)     :: expval(2,nn,nn)
  integer(kind=irg),INTENT(IN)        :: nn, nab(2,nn)

  type(io_T)                          :: Message
  integer(kind=irg)                   :: ir, ic

  do ir=1,nn
! ic is the column index
    do ic=1,nn
      if (ic.ne.ir) then  ! exclude the diagonal
        expval(1,ir,ic) = nab(1,ir)-nab(1,ic)
        expval(2,ir,ic) = nab(2,ir)-nab(2,ic)
      end if
   end do
 end do
 

  end subroutine CalcExpval_

!--------------------------------------------------------------------------
subroutine  CalcSgarray_(self, cell, Diff, sgarray, hklarray, klist, DynFN, numk, nn)
  !DEC$ ATTRIBUTES DLLEXPORT ::  CalcSgarray_
  !! author: MDG
  !! version: 1.0
  !! date: 03/14/20
  !!
  !! Compute the sg array
  
  use mod_crystallography
  use mod_diffraction
  use mod_defect
  use mod_io

  class(ECCI_T), INTENT(INOUT)        :: self
  type(Cell_T),INTENT(INOUT)          :: cell
  type(Diffraction_T),INTENT(INOUT)   :: Diff
  real(kind=sgl),INTENT(INOUT)        :: sgarray(nn, numk)
  real(kind=sgl),INTENT(IN)           :: klist(3,numk)
  real(kind=dbl),INTENT(IN)           :: DynFN(3)
  integer(kind=irg),INTENT(IN)        :: numk, nn, hklarray(3,nn)

  type(io_T)                          :: Message
  integer(kind=irg)                   :: ik, ig, gg(3)

  beamloopCL: do ik=1,numk
    reflectionloopCL: do ig=1,nn
      gg = float(hklarray(1:3,ig))
      sgarray(ig,ik) = Diff%Calcsg(cell,float(gg),sngl(klist(1:3,ik)),sngl(DynFN))
    end do reflectionloopCL
  end do beamloopCL

  end subroutine CalcSgarray_


!--------------------------------------------------------------------------
subroutine getDisplacementField_(self, cell, Defects, ga, gb, disparray)
  !DEC$ ATTRIBUTES DLLEXPORT :: getDisplacementField_
  !! author: MDG
  !! version: 1.0
  !! date: 03/14/20
  !!
  !! Compute the displacement field
  
  use mod_crystallography
  use mod_defect
  use mod_io

  class(ECCI_T), INTENT(INOUT)        :: self
  type(Cell_T),INTENT(IN)             :: cell
  type(Defect_T),INTENT(INOUT)        :: Defects
  real(kind=sgl),INTENT(INOUT)        :: disparray(2,Defects%DF_nums,Defects%DF_npix,Defects%DF_npiy)
  integer(kind=irg),INTENT(IN)        :: ga(3), gb(3)

  type(io_T)                          :: Message
  real(kind=sgl),allocatable          :: imatvals(:,:)
  real(kind=sgl)                      :: gdotR, gac(3), gbc(3)
  real(kind=dbl)                      :: xx
  integer(kind=irg)                   :: i, j, ik, numd
  real(kind=dbl)                      :: io_real(2)

  ! precompute ALL the defect columns and, if needed, store them in dispfile
 ! this portion should be carried out in multi-threaded mode as much as possible
  allocate(imatvals(2,Defects%DF_nums))
  imatvals = 0

  call cell%TransSpace(float(ga),gac,'r','c')
  call cell%TransSpace(float(gb),gbc,'r','c')

  numd = 360

    do i=1,Defects%DF_npix
      do j=1,Defects%DF_npiy
        Defects%DF_R = 0.0
  ! compute the displacement vectors DF_R for all points in the column
        call Defects%CalcR(cell,i,j)
  ! loop over the fixed thickness slices
        do ik=1,Defects%DF_nums
  ! then convert to the dot-product
        if (Defects%DF_R(ik,1).eq.-10000.0) then  ! this is point inside a void
          imatvals(1:2,ik) = -10000
        else  ! it is not a void, so use the full dot product g.R (all vectors must be Cartesian !)
  ! use gac and gbc to get the two dot products and store both of them as integers mapped onto the 0..numd range
  !
  ! due to the logarithmic singularity at a dislocation core, it is possible for gdotR to be NaN; we need
  ! to intercept these cases, and replace the value of gdotR by 0.0
          !write(*,*) Defects%DF_R(ik,1:3)
          gdotR = sum(gac*Defects%DF_R(ik,1:3))
          if (NANCHK(gdotR)) gdotR = 0.0
          xx = numd*amod(gdotR+1000.0,1.0)
          imatvals(1,ik) = xx
          gdotR = sum(gbc*Defects%DF_R(ik,1:3))
          if (NANCHK(gdotR)) gdotR = 0.0
          xx = numd*amod(gdotR+1000.0,1.0)
          imatvals(2,ik) = xx
        end if
      end do ! ik loop
      disparray(1:2,1:Defects%DF_nums,i,j) = imatvals(1:2,1:Defects%DF_nums)
    end do
  end do

  io_real(1) = minval(disparray)
  io_real(2) = maxval(disparray)
  call Message%WriteValue(' --> Disparray bounds: ', io_real, 2, "(2(F10.5,' '))")

  end subroutine getDisplacementField_

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
recursive subroutine Calckvectorcircle_(self, cell, mLambda, khead, k, ga, ktxy, ktrad, ktstep, numk)
!DEC$ ATTRIBUTES DLLEXPORT :: Calckvectorcircle_

use mod_io
use mod_diffraction
use mod_crystallography
use mod_kvectors

IMPLICIT NONE

class(ECCI_T), INTENT(INOUT)        :: self
type(Cell_T),INTENT(IN)             :: cell
type(IO_T)                          :: Message
real(kind=dbl), INTENT(IN)          :: mlambda
type(kvectorlist),pointer           :: khead
integer(kind=irg),INTENT(IN)        :: k(3)
integer(kind=irg),INTENT(IN)        :: ga(3)
real(kind=sgl),INTENT(IN)           :: ktxy(2)
real(kind=sgl),INTENT(IN)           :: ktrad
integer(kind=irg),INTENT(IN)        :: ktstep
integer(kind=irg),INTENT(OUT)       :: numk

type(kvectorlist),pointer           :: ktmp,ktail
integer                             :: istat,imin,imax,jmin,jmax,ijmax,i,j,ic,jc, iang
real                                :: kr(3),glen,delta, dtang, kstar(3),kt(3),gan(3),gperp(3),ktlen, &
                                       dkt, ii, jj, ii1, ii2, jj1, jj2
real(kind=dbl)                      :: ki, iangrad
! compute geometrical factors
 glen = cell%CalcLength(float(ga),'r')         ! length of ga
 gan = ga/glen                                 ! normalized ga
 delta = ktrad*glen/float(ktstep)              ! grid step size in nm-1
 dkt = ktrad/float(ktstep)
 call cell%TransSpace(float(k),kstar,'d','r')       ! transform incident direction to reciprocal space
 call cell%CalcCross(float(ga),kstar,gperp,'r','r',0)! compute g_perp = ga x k
 call cell%NormVec(gperp,'r')                       ! normalize g_perp
 call cell%NormVec(kstar,'r')                       ! normalize reciprocal beam vector

! deal only with the incident beam (parallel illumination)
if (ktstep.eq.0) then
 if (.not.associated(khead)) then     ! allocate the head and ktail of the linked list
   allocate(khead,stat=istat)         ! allocate new value
   if (istat.ne.0) call Message%printError('Calckvectorcone: unable to allocate head pointer',' ')
   ktail => khead                      ! ktail points to new value
   nullify(ktail%next)                ! nullify next in new value
   numk = 1                          ! keep track of number of k-vectors so far
 ! this should be the center vector of the illumination cone !!!
   kt = - glen * (ktxy(1)*gan + ktxy(2) * gperp)
   ktail%kt = kt                           ! store tangential component of k
   ktlen = glen**2*(ktxy(1)**2+ktxy(2)**2)         ! squared length of tangential component

   kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar ! complete wave vector
   ktail%k = kr                            ! store in pointer list
   ktail%kn = cell%CalcDot(ktail%k,dble(kstar),'r')    ! normal component of k
 end if
else
! next, put the center of the cone in units of (i,j) (original ECP "screen" coordinates)
  ic = int(ktxy(1)*glen/delta)
  jc = int(ktxy(2)*glen/delta)

  imin =  -ki; imax = ki; jmin = -ki; jmax = ki;
  ijmax = (ki)**2

  if (.not.associated(khead)) then 
    allocate(khead,stat=istat)
    ktail => khead
  end if

  imax = 360
  ki = dble(ktrad)*dtor
  numk = 0

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
      
      numk = numk + 1                 ! keep track of number of k-vectors so far
      kt = - (ii)*gan*delta- (jj)*gperp*delta  ! tangential component of k
      ktail%kt = kt                    ! store tangential component of k
      ktlen = delta**2*(ii**2+jj**2)         ! squared length of tangential component

      kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar ! complete wave vector
      ktail%k = kr                     ! store in pointer list
      ktail%kn = cell%CalcDot(ktail%k,dble(kstar),'r')
      ktail => ktail%next               ! ktail points to new value
      nullify(ktail%next)              ! nullify next in new value
    end if
  end do
 end if

end subroutine Calckvectorcircle_

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
recursive subroutine Calckvectorcone_(self, cell, mLambda, khead, k, ga, ktxy, ktrad, ktstep, numk)
!DEC$ ATTRIBUTES DLLEXPORT :: Calckvectorcone_

use mod_io
use mod_diffraction
use mod_crystallography
use mod_kvectors

IMPLICIT NONE

class(ECCI_T), INTENT(INOUT)        :: self
type(Cell_T),INTENT(IN)             :: cell
type(IO_T)                          :: Message
real(kind=dbl), INTENT(IN)          :: mlambda
type(kvectorlist),pointer           :: khead
integer(kind=irg),INTENT(IN)        :: k(3)
integer(kind=irg),INTENT(IN)        :: ga(3)
real(kind=sgl),INTENT(IN)           :: ktxy(2)
real(kind=sgl),INTENT(IN)           :: ktrad
integer(kind=irg),INTENT(IN)        :: ktstep
integer(kind=irg),INTENT(OUT)       :: numk

type(kvectorlist),pointer           :: ktmp,ktail
integer                             :: istat,imin,imax,jmin,jmax,ijmax,i,j,ic,jc,ki
real                                :: kr(3),glen,delta,kstar(3),kt(3),gan(3),gperp(3),ktlen, dkt

! compute geometrical factors
 glen = cell%CalcLength(float(ga),'r')         ! length of ga
 gan = ga/glen                                 ! normalized ga
 delta = ktrad*glen/float(ktstep)              ! grid step size in nm-1
 dkt = ktrad/float(ktstep)
 call cell%TransSpace(float(k),kstar,'d','r')       ! transform incident direction to reciprocal space
 call cell%CalcCross(float(ga),kstar,gperp,'r','r',0)! compute g_perp = ga x k
 call cell%NormVec(gperp,'r')                       ! normalize g_perp
 call cell%NormVec(kstar,'r')                       ! normalize reciprocal beam vector

! deal only with the incident beam (parallel illumination)
if (ktstep.eq.0) then
 if (.not.associated(khead)) then     ! allocate the head and ktail of the linked list
   allocate(khead,stat=istat)         ! allocate new value
   if (istat.ne.0) call Message%printError('Calckvectorcone: unable to allocate head pointer',' ')
   ktail => khead                      ! ktail points to new value
   nullify(ktail%next)                ! nullify next in new value
   numk = 1                          ! keep track of number of k-vectors so far
 ! this should be the center vector of the illumination cone !!!
   kt = - glen * (ktxy(1)*gan + ktxy(2) * gperp)
   ktail%kt = kt                           ! store tangential component of k
   ktlen = glen**2*(ktxy(1)**2+ktxy(2)**2)         ! squared length of tangential component

   kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar ! complete wave vector
   ktail%k = kr                            ! store in pointer list
   ktail%kn = cell%CalcDot(ktail%k,dble(kstar),'r')    ! normal component of k
 end if
else
! next, put the center of the cone in units of (i,j) (original ECP "screen" coordinates)
  ic = int(ktxy(1)*glen/delta)
  jc = int(ktxy(2)*glen/delta)
  ki = ktstep

 if (.not.associated(khead)) then     ! allocate the head and ktail of the linked list
   allocate(khead,stat=istat)         ! allocate new value
   if (istat.ne.0) call Message%printError('Calckvectorcone: unable to allocate head pointer',' ')
   ktail => khead                      ! ktail points to new value
   nullify(ktail%next)                ! nullify next in new value
   numk = 1                          ! keep track of number of k-vectors so far
 ! this should be the center vector of the illumination cone !!!
   ktail%i = ic                            ! i-index of beam
   ktail%j = jc                            ! j-index of beam
   kt = -float(ktail%i)*delta*gan - float(ktail%j)*delta*gperp  ! tangential component of k
   ktail%kt = kt                           ! store tangential component of k
   ktlen = delta**2*(ktail%i**2+ktail%j**2)         ! squared length of tangential component

   kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar ! complete wave vector
   ktail%k = kr                            ! store in pointer list
   ktail%kn = cell%CalcDot(ktail%k,dble(kstar),'r')    ! normal component of k
 else
   call Message%printError('Calckvectorcone: pointer head already allocated',' ')
 end if

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
      numk = numk + 1                 ! keep track of number of k-vectors so far
      ktail%i = ic+i                   ! i-index of beam
      ktail%j = jc+j                   ! j-index of beam
      kt = - float(ktail%i)*delta*gan - float(ktail%j)*delta*gperp  ! tangential component of k
      ktail%kt = kt                    ! store tangential component of k
      ktlen = delta**2*(ktail%i**2+ktail%j**2)         ! squared length of tangential component

      kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar ! complete wave vector
      ktail%k = kr                     ! store in pointer list
      ktail%kn = cell%CalcDot(ktail%k,dble(kstar),'r')    ! normal component of k
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
!> @param khead head of k-vector list
!> @param k incident wave vector (zone axis)
!> @param ga principal g vector
!> @param ktxy tangential components of trace start point
!> @param ktxy2 tangential components of trace end point
!> @param ktrad cone opening angle
!> @param ktstep number of steps along cone radius
!> @param numk resulting number of incident beam directions
!
!> @date 12/08/13 MDG 1.0 original
!> @date 12/01/15 MDG 1.1 simplifcation of input variables
!--------------------------------------------------------------------------
recursive subroutine Calckvectortrace_(self,cell,mlambda,khead,k,ga,ktxy,ktxy2,ktrad,ktstep,numk)
!DEC$ ATTRIBUTES DLLEXPORT :: Calckvectortrace_

use mod_io
use mod_diffraction
use mod_crystallography
use mod_kvectors

IMPLICIT NONE

class(ECCI_T), INTENT(INOUT)        :: self
type(Cell_T),INTENT(IN)             :: cell
type(IO_T)                          :: Message
real(kind=dbl), INTENT(IN)          :: mlambda
type(kvectorlist),pointer           :: khead
integer(kind=irg),INTENT(IN)        :: k(3)
integer(kind=irg),INTENT(IN)        :: ga(3)
real(kind=sgl),INTENT(IN)           :: ktxy(2)
real(kind=sgl),INTENT(IN)           :: ktxy2(2)
real(kind=sgl),INTENT(IN)           :: ktrad
integer(kind=irg),INTENT(IN)        :: ktstep
integer(kind=irg),INTENT(OUT)       :: numk

type(kvectorlist),pointer           :: ktail
integer                             :: istat,j
real                                :: kr(3),glen,delta,kstar(3),kt(3),gan(3),gperp(3),ktlen, dktx, dkty

! compute geometrical factors
 glen = cell%CalcLength(float(ga),'r')              ! length of ga
 gan = ga/glen                                 ! normalized ga
 delta = 2.0*ktrad*glen/float(2*ktstep+1)      ! grid step size in nm-1
 call cell%TransSpace(float(k),kstar,'d','r')       ! transform incident direction to reciprocal space
 call cell%CalcCross(float(ga),kstar,gperp,'r','r',0)! compute g_perp = ga x k
 call cell%NormVec(gperp,'r')                       ! normalize g_perp
 call cell%NormVec(kstar,'r')                       ! normalize reciprocal beam vector

 j = 0
 if (.not.associated(khead)) then     ! allocate the head and ktail of the linked list
   allocate(khead,stat=istat)         ! allocate new value
   if (istat.ne.0) call Message%printError('Calckvectortrace: unable to allocate head pointer',' ')
   ktail => khead                     ! ktail points to new value
   nullify(ktail%next)                ! nullify next in new value
   numk = 1                           ! keep track of number of k-vectors so far
! this should be the starting point of the line trace
!   kt = - glen * ( ktxy(1)*gan + ktxy(2) * gperp)
   kt = - glen * ( ktxy(1)*gan - ktxy(2) * gperp)
   ktail%kt = kt                           ! store tangential component of k
   ktlen = glen**2*(ktxy(1)**2+ktxy(2)**2)         ! squared length of tangential component
   kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar ! complete wave vector
   ktail%k = kr                            ! store in pointer list
   ktail%kn = cell%CalcDot(ktail%k,dble(kstar),'r')    ! normal component of k
 end if

 dktx = (ktxy2(1) - ktxy(1))/float(ktstep-1)
 dkty = (ktxy2(2) - ktxy(2))/float(ktstep-1)

 do j=1,ktstep-1
      allocate(ktail%next,stat=istat)  ! allocate new value
      if (istat.ne.0) call Message%printError('Calckvectortrace: unable to allocate pointer',' ')
      ktail => ktail%next              ! ktail points to new value
      nullify(ktail%next)              ! nullify next in new value
      numk = numk + 1                  ! keep track of number of k-vectors so far
!     kt = - glen * ( (ktxy(1)+float(j)*dktx)*gan + (ktxy(2)+float(j)*dkty) * gperp) ! tangential component of k
      kt = - glen * ( (ktxy(1)+float(j)*dktx)*gan - (ktxy(2)+float(j)*dkty) * gperp) ! tangential component of k
      ktail%kt = kt                    ! store tangential component of k
      ktlen = delta**2*(ktail%i**2+ktail%j**2)         ! squared length of tangential component
      kr = kt + sqrt(1.0/mLambda**2 - ktlen)*kstar ! complete wave vector
      ktail%k = kr                     ! store in pointer list
      ktail%kn = cell%CalcDot(ktail%k,dble(kstar),'r')    ! normal component of k
 end do

end subroutine Calckvectortrace_


!C***********************************************************************
!C
!C                        naninfchk.f
!C
!C      *****************************************************************
!C      *                                                               *
!C	* 	Absoft Corporation 					*
!C 	*	2781 Bond Street					*
!C	*	Rochester Hills, MI  48309				*
!C	*								*
!C	*	This file contains example code for demonstration	*
!C	*	purposes only.  Absoft makes no warranty of the	*
!C	*	suitability of this code for any purpose.		*
!C	*								*
!C	*	In no event shall Absoft be liable for any incidental,*
!C	*	indirect, special, or consequential damages arising	*
!C	*	out of the use of this code.				*
!C	*								*
!C	*****************************************************************
!C
!C Routines to test real and double values against NaN and INF
!C
!C            NANCHK(X) - tests REAL*4 value X against NaN
!!C            DNANCHK(X) - tests REAL*8 value X against NaN
!!C            INFCHK(X) - tests REAL*4 value X against INF
!!C            DINFCHK(X) - test REAL*8 value X against INF
!C
!C For little endian machines (Intel x86), compile with
!C
!C      f77 -c -DBYTE_SWAPPED=1 naninfchk.f
!C	or
!C      f90 -c -DBYTE_SWAPPED=1 naninfchk.f -YBOZTYPE=INT
!C
!C For big endian machines (PowerPC), compile with
!C
!C      f77 -c naninfchk.f
!C	or
!C      f90 -c naninfchk.f -YBOZTYPE=INT
!C
!C***********************************************************************
RECURSIVE LOGICAL FUNCTION NANCHK(X)
!DEC$ ATTRIBUTES DLLEXPORT :: NANCHK

IMPLICIT NONE
REAL,INTENT(IN)      :: X
REAL                 :: Y
INTEGER              :: I
EQUIVALENCE(Y,I)

Y = X
NANCHK = isnan(Y) !((I .AND. z'7f800000') .EQ. z'7f800000') .AND.((I .AND. z'007fffff') .NE. z'00000000')

RETURN
END

end module mod_ECCI
