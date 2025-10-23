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

module mod_sampler
  !! author: MDG 
  !! version: 1.0 
  !! date: 09/02/25
  !!
  !! class definition for the EMsampler program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMsampler program
type, public :: samplerNameListType
  real(kind=dbl)    :: kappa(2)           ! concentration parameters
  real(kind=dbl)    :: dir1(3)            ! Rodrigues mean direction
  real(kind=dbl)    :: dir2(3)            ! Rodrigues mean direction
  real(kind=dbl)    :: dir3(3)            ! Rodrigues mean direction
  real(kind=dbl)    :: dir4(3)            ! Rodrigues mean direction
  real(kind=dbl)    :: viewangle          ! view angle of the observer in PoVray renderings
  integer(kind=irg) :: norientations      ! number of orientations per dataset
  integer(kind=irg) :: pgnum              ! point group number (for Laue point group)
  integer(kind=irg) :: seed1              ! seed 1 for pseudo-random number generator [sampling]
  integer(kind=irg) :: seed2              ! seed 2 for pseudo-random number generator [averaging]
  logical           :: reduce             ! reduce dir# parameters to RFZ before sampling ?
  character(fnlen)  :: hdfname            ! name of output HDF5 file
  character(fnlen)  :: imagefolder        ! folder path w.r.t EMdatapathname for image files
  character(fnlen)  :: prefix             ! prefix for image files (will be png files)
  character(fnlen)  :: PVexec             ! path to PoVray executable
  character(fnlen)  :: PVincludepath      ! path to PoVray include files
end type samplerNameListType

! class definition
type, public :: sampler_T
private 
  character(fnlen)            :: nmldeffile = 'EMsampler.nml'
  type(samplerNameListType)   :: nml 
  character(fnlen)            :: PVexec
  character(fnlen)            :: PVincludepath

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: sampler_
  procedure, pass(self) :: renderdata_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: sampler => sampler_
  generic, public :: renderdata => renderdata_

end type sampler_T

! the constructor routine for this class 
interface sampler_T
  module procedure sampler_constructor
end interface sampler_T

contains

!--------------------------------------------------------------------------
type(sampler_T) function sampler_constructor( nmlfile ) result(sampler)
!! author: MDG 
!! version: 1.0 
!! date: 09/02/25
!!
!! constructor for the sampler_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call sampler%readNameList(nmlfile)

end function sampler_constructor

!--------------------------------------------------------------------------
subroutine sampler_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 09/02/25
!!
!! destructor for the sampler_T Class
 
IMPLICIT NONE

type(sampler_T), INTENT(INOUT)  :: self 


end subroutine sampler_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 09/02/25
!!
!! read the namelist from an nml file for the sampler_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(sampler_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

real(kind=dbl)    :: kappa(2)           ! concentration parameters
real(kind=dbl)    :: dir1(3)            ! Rodrigues mean direction
real(kind=dbl)    :: dir2(3)            ! Rodrigues mean direction
real(kind=dbl)    :: dir3(3)            ! Rodrigues mean direction
real(kind=dbl)    :: dir4(3)            ! Rodrigues mean direction
real(kind=dbl)    :: viewangle          ! view angle of the observer in PoVray renderings
integer(kind=irg) :: norientations      ! number of orientations per dataset
integer(kind=irg) :: pgnum              ! point group number (for Laue point group)
integer(kind=irg) :: seed1              ! seed 1 for pseudo-random number generator [sampling]
integer(kind=irg) :: seed2              ! seed 2 for pseudo-random number generator [averaging]
logical           :: reduce             ! reduce dir# parameters to RFZ before sampling ?
character(fnlen)  :: hdfname            ! name of output HDF5 file
character(fnlen)  :: imagefolder        ! folder path w.r.t EMdatapathname for image files
character(fnlen)  :: prefix             ! prefix for image files (will be png files)
character(fnlen)  :: PVexec             ! path to PoVray executable
character(fnlen)  :: PVincludepath      ! path to PoVray include files

namelist / EMsampler / kappa, dir1, dir2, dir3, dir4, norientations, pgnum, hdfname, reduce, &
                       imagefolder, prefix, PVexec, PVincludepath, seed1, seed2, viewangle 

kappa = (/ 512.D0, 3456.D0 /)
dir1 = (/ 0.D0, 0.D0, 0.D0 /)
dir2 = (/ 0.D0, 0.D0, 0.D0 /) 
dir3 = (/ 0.D0, 0.D0, 0.D0 /) 
dir4 = (/ 0.D0, 0.D0, 0.D0 /) 
norientations = 5000
pgnum = 1
seed1 = 54321
seed2 = 3514
reduce = .TRUE.
hdfname = 'undefined'
imagefolder = 'undefined'
prefix = 'undefined'
PVexec = 'undefined'
PVincludepath = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
  open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
  read(UNIT=dataunit,NML=EMsampler)
  close(UNIT=dataunit,STATUS='keep')

  if (trim(hdfname).eq.'undefined') then
      call Message%printError('readNameList:',' hdfname is undefined in '//nmlfile)
  end if

  if (PVexec.ne.'undefined') then
    if (trim(imagefolder).eq.'undefined') then
        call Message%printError('readNameList:',' imagefolder is undefined in '//nmlfile)
    end if
    if (trim(prefix).eq.'undefined') then
        call Message%printError('readNameList:',' prefix is undefined in '//nmlfile)
    end if
    if (trim(PVincludepath).eq.'undefined') then
        call Message%printError('readNameList:',' PVincludepath is undefined in '//nmlfile)
    end if
    self%PVexec = PVexec
    self%PVincludepath = PVincludepath
  end if 
end if

self%nml%kappa = kappa
self%nml%dir1 = dir1
self%nml%dir2 = dir2
self%nml%dir3 = dir3
self%nml%dir4 = dir4
self%nml%viewangle = viewangle
self%nml%norientations = norientations
self%nml%pgnum = pgnum  
self%nml%seed1 = seed1  
self%nml%seed2 = seed2  
self%nml%reduce = reduce
self%nml%hdfname = hdfname
self%nml%imagefolder = imagefolder
self%nml%prefix = prefix
self%nml%PVexec = PVexec
self%nml%PVincludepath = PVincludepath

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 09/02/25
!!
!! pass the namelist for the sampler_T Class to the calling program

IMPLICIT NONE 

class(sampler_T), INTENT(INOUT)          :: self
type(samplerNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 09/02/25
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(sampler_T), INTENT(INOUT)         :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 4, n_real = 9
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( enl => self%nml )

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! write all the single integers
io_int = (/ enl%norientations, enl%pgnum, enl%seed1, enl%seed2 /)
intlist(1) = 'norientations'
intlist(2) = 'pgnum'
intlist(3) = 'seed1'
intlist(4) = 'seed2'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! a 2-vector
dataset = 'kappa'
hdferr = HDF%writeDatasetDoubleArray(dataset, enl%kappa, 2)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create kappa dataset', hdferr)

! 3-vectors
dataset = 'dir1'
hdferr = HDF%writeDatasetDoubleArray(dataset, enl%dir1, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create dir1 dataset', hdferr)

dataset = 'dir2'
hdferr = HDF%writeDatasetDoubleArray(dataset, enl%dir2, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create dir2 dataset', hdferr)

dataset = 'dir3'
hdferr = HDF%writeDatasetDoubleArray(dataset, enl%dir3, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create dir3 dataset', hdferr)

dataset = 'dir4'
hdferr = HDF%writeDatasetDoubleArray(dataset, enl%dir4, 3)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create dir4 dataset', hdferr)

dataset = 'hdfname'
line2(1) = trim(enl%hdfname)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create hdfname dataset', hdferr)

dataset = 'imagefolder'
line2(1) = trim(enl%imagefolder)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create imagefolder dataset', hdferr)

dataset = 'prefix'
line2(1) = trim(enl%prefix)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create prefix dataset', hdferr)

dataset = 'PVexec'
line2(1) = trim(enl%PVexec)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create PVexec dataset', hdferr)

dataset = 'PVincludepath'
line2(1) = trim(enl%PVincludepath)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create PVincludepath dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine sampler_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: sampler_
!! author: MDG 
!! version: 1.0 
!! date: 09/02/25
!!
!! perform the computations
!!
!! this routine computes four orientation data sets each for two different 
!! concentration parameters kappa and for both Watson and von Mises-Fisher 
!! distributions (16 orientations sets in all).
!! Then for each orientation data set, both WAT and vMF averaging are performed.
!! All results are stored in an HDF file and if PVexec is set, then rendered
!! 3D stereographic projections are generated for every dataset using PoVray

use mod_EMsoft
use mod_HDFnames
use mod_crystallography
use mod_quaternions
use mod_dirstats
use mod_io
use mod_symmetry
use mod_rotations
use mod_timing
use mod_so3
use mod_povray
use mod_OrientationViz
use HDF5
use mod_HDFsupport
use stringconstants
use ISO_C_BINDING

IMPLICIT NONE 

class(sampler_T), INTENT(INOUT)         :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

type(DirStat_T)                         :: dictVMF, dictWAT
type(Quaternion_T)                      :: mu, muhat, meanquat(4)
type(QuaternionArray_T)                 :: qAR, dummy, qsym
type(IO_T)                              :: Message
type(SpaceGroup_T)                      :: SG
type(HDF_T)                             :: HDF
type(Cell_T)                            :: cell
type(q_T)                               :: qFZ, q
type(a_T)                               :: a
type(r_T)                               :: r
type(o_T)                               :: o
type(so3_T)                             :: SO
type(Timing_T)                          :: timer
type(PoVRay_T)                          :: PoVst

real(kind=dbl),allocatable              :: vMFquatarray(:,:,:), WATquatarray(:,:,:), vMFkappahat(:), inputquat(:,:), &
                                           WATkappahat(:), vMFmuhat(:,:), WATmuhat(:,:), misor(:,:)
integer(kind=irg)                       :: setcnt, seed1, seed2, hdferr, i, j, k, io_int(2), FZtype, FZorder, seedarray(2,2,8)
real(kind=dbl)                          :: rod(4), rodL, kappahat
real(kind=sgl)                          :: tstop
character(11)                           :: dstr
character(15)                           :: tstrb
character(15)                           :: tstre
character(1)                            :: ch
character(fnlen)                        :: dataset, datagroupname, fname, attributename, HDF_FileVersion, nmldeffile, outname
logical                                 :: overwrite = .TRUE.
character(fnlen,kind=c_char)            :: line2(1)


call openFortranHDFInterface()
nmldeffile = trim(EMsoft%nmldeffile)

timer = Timing_T()
tstrb = timer%getTimeString()
dstr = timer%getDateString()

! first get the name list
associate( enl => self%nml )

allocate(vMFquatarray( 4, enl%norientations, 8), WATquatarray( 4, enl%norientations, 8) )
allocate( vMFkappahat(8), WATkappahat(8), vMFmuhat(4,8), WATmuhat(4,8), misor(2,8), inputquat(4,4) )

! put the 4 mean directions in quaternion format 
! dir1 
rodL = sqrt(sum(enl%dir1**2))
rod(1:3) = enl%dir1(1:3)/rodL
rod(4) = rodL
r = r_T( rdinp = rod )
call r%r_print(' rod 1: ')
q = r%rq()
mu = Quaternion_T( qd = q%q_copyd() )
call mu%quat_print(' dir 1: ')
meanquat(1) = mu 
! dir2 
rodL = sqrt(sum(enl%dir2**2))
rod(1:3) = enl%dir2(1:3)/rodL
rod(4) = rodL
r = r_T( rdinp = rod )
call r%r_print(' rod 2: ')
q = r%rq()
mu = Quaternion_T( qd = q%q_copyd() )
call mu%quat_print(' dir 2: ')
meanquat(2) = mu 
! dir3 
rodL = sqrt(sum(enl%dir3**2))
rod(1:3) = enl%dir3(1:3)/rodL
rod(4) = rodL
r = r_T( rdinp = rod )
call r%r_print(' rod 3: ')
q = r%rq()
mu = Quaternion_T( qd = q%q_copyd() )
call mu%quat_print(' dir 3: ')
meanquat(3) = mu 
! dir4 
rodL = sqrt(sum(enl%dir4**2))
rod(1:3) = enl%dir4(1:3)/rodL
rod(4) = rodL
r = r_T( rdinp = rod )
call r%r_print(' rod 4: ')
q = r%rq()
mu = Quaternion_T( qd = q%q_copyd() )
call mu%quat_print(' dir 4: ')
meanquat(4) = mu 

!====================================
!====================================
!====================================
! first we'll create the HDF file and leave it open 
call Message%printMessage(' Initializing HDF output file')
fname = EMsoft%generateFilePath('EMdatapathname',enl%hdfname)
HDF = HDF_T()

hdferr =  HDF%createFile(fname)
if (hdferr.ne.0) call HDF%error_check('HDF_createFile ', hdferr)

!====================================
! new in Release 4.3: add a Manufacturer string (null terminated)
dataset = SC_Manufacturer
line2(1) = 'EMsoftOO'
line2(1) = cstringify(line2(1))
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
!====================================

! write the EMheader to the file
datagroupname = trim(HDFnames%get_ProgramData())
call HDF%writeEMheader(EMsoft, dstr, tstrb, tstre, progname, datagroupname)

! create a namelist group to write all the namelist files into
hdferr = HDF%createGroup(HDFnames%get_NMLfiles())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup NMLfiles', hdferr)

! read the text file and write the array to the file
dataset = SC_SamplerNML
hdferr = HDF%writeDatasetTextFile(dataset, nmldeffile)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetTextFile ', hdferr)

call HDF%pop()

! create a NMLparameters group to write all the namelist entries into
hdferr = HDF%createGroup(HDFnames%get_NMLparameters())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup NMLparameters', hdferr)

call self%writeHDFNameList_(HDF, HDFnames)

! and leave this group
call HDF%pop()

! then the remainder of the data in a EMData group
hdferr = HDF%createGroup(HDFnames%get_EMData())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup EMData', hdferr)

! create the EMSampler group and add a HDF_FileVersion attribute to it
hdferr = HDF%createGroup(datagroupname)
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup EBSD/TKD', hdferr)
HDF_FileVersion = '4.1'
attributename = SC_HDFFileVersion
hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)

! and we leave this file open for further output
!====================================
!====================================
!====================================

! start the computations
call Message%printMessage(' Starting computation for Laue point group '//PGTHD(enl%pgnum))

! first put the mean directions inside the RFZ
SO = so3_T( pgnum = enl%pgnum )
! get the symmetry operator quaternions for the point group
call dummy%QSym_Init(enl%pgnum, qsym)

if (enl%reduce.eqv..TRUE.) then
  do i=1,4
    q = q_T( qdinp = meanquat(i)%get_quatd() )
    call SO%ReduceOrientationtoRFZ( q, qsym, r )
    q = r%rq()
    mu = Quaternion_T( qd = q%q_copyd() )
    meanquat(i) = mu 
  end do 
end if

do i=1,4
  inputquat(1:4,i) = meanquat(i)%get_quatd()
end do  

! we'll do vMF sampling first and generate 8 orientation data sets 
dictVMF = DirStat_T( DStype='VMF', PGnum=enl%pgnum )
call dictVMF%setNumEM(25)
call dictVMF%setNumIter(30)

setcnt = 1 
do i=1,2      ! loop over the kappa concentration parameter values 
  do j=1,4    ! loop over the mean quaternion directions
    seedarray(1:2,1,setcnt) = (/ seed1, seed2 /)
    qAR = dictVMF%SampleDS( enl%norientations, seed1, meanquat(j), enl%kappa(i) )
! copy the orientations into the vMFquatarray for storage in the output HDF file
    do k=1,enl%norientations 
      mu = qAR%getQuatfromArray(k)
      vMFquatarray(1:4,k,setcnt) = mu%get_quatd()
    end do
    io_int = (/ i, j /)
    call Message%WriteValue(' Averaging vMF data set ', io_int, 2)
! and we might as well do the averaging at this point ...
    call dictVMF%setQuatArray( qAR )
    muhat = Quaternion_T( qd=(/ 1.D0, 0.D0, 0.D0,0.D0 /) )
    call dictVMF%EMforDS( seed2, muhat, kappahat )
    vMFkappahat(setcnt) = kappahat
    vMFmuhat(1:4,setcnt) = muhat%get_quatd()
    misor(1,setcnt) = 2.0 * acos( sum( muhat%get_quatd() * meanquat(j)%get_quatd() ) )/dtor
    write (*,*) ' input    : ', meanquat(j)%get_quatd(), enl%kappa(i)
    write (*,*) ' averaged : ', muhat%get_quatd(), kappahat
    write (*,*) ' misor (°): ', misor(1,setcnt)
! generate the PoVray file and run the program
    if (trim(enl%prefix).ne.'undefined') then 
      write (ch,"(I1)") setcnt
      outname = trim(enl%imagefolder)//'/'//trim(enl%prefix)//'-vMF-'//ch
      outname = EMsoft%generateFilePath('EMdatapathname', outname)
      call self%renderdata_(EMsoft,outname,enl%imagefolder,qAR,enl%pgnum)
    end if
    setcnt = setcnt + 1
 end do 
end do

! and write these results to the HDF file
dataset = 'vMFquatarray'
hdferr = HDF%writeDatasetDoubleArray(dataset, vMFquatarray, 4, enl%norientations, 8)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetDoubleArray3D vMFquatarray', hdferr)

dataset = 'vMFkappahat'
hdferr = HDF%writeDatasetDoubleArray(dataset, vMFkappahat, 8)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetDoubleArray1D vMFkappahat', hdferr)

dataset = 'vMFmuhat'
hdferr = HDF%writeDatasetDoubleArray(dataset, vMFmuhat, 4, 8)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetDoubleArray2D vMFmuhat', hdferr)


! next do the Watson sampling and averaging
dictWAT = DirStat_T( DStype='WAT', PGnum=enl%pgnum )
call dictWAT%setNumEM(25)
call dictWAT%setNumIter(30)

setcnt = 1 
do i=1,2      ! loop over the kappa concentration parameter values 
  do j=1,4    ! loop over the mean quaternion directions
    seedarray(1:2,2,setcnt) = (/ seed1, seed2 /)
    qAR = dictWAT%SampleDS( enl%norientations, seed1, meanquat(j), enl%kappa(i) )
! copy the orientations into the WATquatarray for storage in the output HDF file
    do k=1,enl%norientations 
      mu = qAR%getQuatfromArray(k)
      WATquatarray(1:4,k,setcnt) = mu%get_quatd()
    end do
    io_int = (/ i, j /)
    call Message%WriteValue(' Averaging WAT data set ', io_int, 2)
! and we might as well do the averaging at this point ...
    call dictWAT%setQuatArray( qAR )
    muhat = Quaternion_T( qd=(/ 1.D0, 0.D0, 0.D0,0.D0 /) )
    call dictWAT%EMforDS( seed2, muhat, kappahat )
    WATkappahat(setcnt) = kappahat
    WATmuhat(1:4,setcnt) = muhat%get_quatd()
    misor(2,setcnt) = 2.0 * acos( sum( muhat%get_quatd() * meanquat(j)%get_quatd() ) )/dtor
    write (*,*) ' input    : ', meanquat(j)%get_quatd(), enl%kappa(i)
    write (*,*) ' averaged : ', muhat%get_quatd(), kappahat
    write (*,*) ' misor (°): ', misor(2,setcnt)
! generate the PoVray file and run the program
    if (trim(enl%prefix).ne.'undefined') then 
      write (ch,"(I1)") setcnt
      outname = trim(enl%imagefolder)//'/'//trim(enl%prefix)//'-WAT-'//ch
      outname = EMsoft%generateFilePath('EMdatapathname', outname)
      call self%renderdata_(EMsoft,outname,enl%imagefolder,qAR,enl%pgnum)
    end if
    setcnt = setcnt + 1
 end do 
end do

! and write these results to the HDF file
dataset = 'WATquatarray'
hdferr = HDF%writeDatasetDoubleArray(dataset, WATquatarray, 4, enl%norientations, 8)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetDoubleArray3D WATquatarray', hdferr)

dataset = 'WATkappahat'
hdferr = HDF%writeDatasetDoubleArray(dataset, WATkappahat, 8)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetDoubleArray1D WATkappahat', hdferr)

dataset = 'WATmuhat'
hdferr = HDF%writeDatasetDoubleArray(dataset, WATmuhat, 4, 8)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetDoubleArray2D WATmuhat', hdferr)

dataset = 'seedarray'
hdferr = HDF%writeDatasetIntegerArray(dataset, seedarray, 2, 2, 8)
if (hdferr.ne.0) call HDF%error_check('writeDatasetIntegerArray seedarray', hdferr)

dataset = 'misor'
hdferr = HDF%writeDatasetDoubleArray(dataset, misor, 2, 8)
if (hdferr.ne.0) call HDF%error_check('writeDatasetDoubleArray misor', hdferr)

dataset = 'inputquat'
hdferr = HDF%writeDatasetDoubleArray(dataset, inputquat, 4, 4)
if (hdferr.ne.0) call HDF%error_check('writeDatasetDoubleArray inputquat', hdferr)

call HDF%pop() 
call HDF%pop() 

! and update the end time
call timer%Time_tock()
tstop = timer%getInterval()

timer = timing_T()
tstre = timer%getTimeString()

hdferr = HDF%openGroup(HDFnames%get_EMheader())
if (hdferr.ne.0) call HDF%error_check('HDF_openGroup EMheader', hdferr)

hdferr = HDF%openGroup(HDFnames%get_ProgramData())
if (hdferr.ne.0) call HDF%error_check('HDF_openGroup ProgramData', hdferr)

! stop time /EMheader/StopTime 'character'
dataset = SC_StopTime
line2(1) = dstr//', '//tstre
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetStringArray StopTime', hdferr)

dataset = SC_Duration
hdferr = HDF%writeDatasetFloat(dataset, tstop)
if (hdferr.ne.0) call HDF%error_check('HDF_writeDatasetFloat Duration', hdferr)

! close the datafile
call HDF%popall()

call closeFortranHDFInterface()

end associate

end subroutine sampler_

!--------------------------------------------------------------------------
subroutine renderdata_(self, EMsoft, povname, subfolder, qAR, pgnum)
!DEC$ ATTRIBUTES DLLEXPORT :: renderdata_
!! author: MDG 
!! version: 1.0 
!! date: 09/04/25
!!
!! render the orientation data into a stereographic projection using PoVray

use mod_EMsoft
use mod_crystallography
use mod_quaternions
use mod_io
use mod_symmetry
use mod_rotations
use mod_so3
use mod_povray
use mod_OrientationViz
use stringconstants
use ISO_C_BINDING

IMPLICIT NONE 

class(sampler_T), INTENT(INOUT)         :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(IN)            :: povname
character(fnlen), INTENT(IN)            :: subfolder
type(QuaternionArray_T), INTENT(INOUT)  :: qAR
integer(kind=irg), INTENT(IN)           :: pgnum

type(so3_T)                             :: SO
type(IO_T)                              :: Message
type(QuaternionArray_T)                 :: dummy, qsym 
type(PoVRay_T)                          :: PoV
type(q_T)                               :: q
type(s_T)                               :: st

integer(kind=irg)                       :: FZtype, FZorder, num, norientations, dFZ, ix 
real(kind=sgl)                          :: eyepos(3), sphrad, xyz(3), dd
real(kind=dbl)                          :: cylr
character(fnlen)                        :: locationline, lightline, skyline, colorstring, rgbstring, str, pvcmd, povfile
character(9)                            :: px, py, pz, pd
character(21)                           :: p1, p2
type(FZpointd),pointer                  :: FZtmp
logical                                 :: fexists 

! set up the SO3 class
SO = so3_T( pgnum, zerolist='FZ')
call SO%getFZtypeandorder(FZtype, FZorder)
call SO%QuaternionArraytonewlist(qAR, 'FZ')

! get the symmetry operator quaternions for the point group
call dummy%QSym_Init(pgnum, qsym)
num = qsym%getQnumber()
norientations = qAR%getQnumber()

! initialize PoVray parameters
dd = 2.5
write (pd,"(F9.3)") dd
locationline = "location < "
eyepos = (/ 0.911259, 0.0, 0.412 /)
eyepos = eyepos/sqrt( sum( eyepos*eyepos))
write (px,"(F9.3)") eyepos(1)
write (py,"(F9.3)") eyepos(2)
write (pz,"(F9.3)") eyepos(3)

p1 = "*cos(clck*0.0174533)"
p2 = "*sin(clck*0.0174533)"

locationline = trim(locationline)//px//p1//"-"//py//p2//","//px//p2//"+"//py//p1//","//pz//">*"//pd

povfile = trim(povname)//'-st.pov'
PoV = PoVRay_T( EMsoft, povfile, locationline=locationline, viewangle = self%nml%viewangle )
call PoV%toggleVerbose()
if (FZorder.lt.0) call PoV%set_roto(abs(FZorder))

! reduce to RFZ
call SO%ReducelisttoRFZ(qsym)
FZtmp => SO%getListHead('FZ')          ! point to the top of the list

! draw RFZ outline
cylr = 0.0015D0
sphrad = 0.005
dFZ = 3  ! for stereographic projections
call PoV%drawFZ(SO, dFZ, cylr, outline=1)
! open the union of spheres...
write (90,"('union { ')")

pointloop: do ix = 1,norientations
  q = FZtmp%qu
  st = q%qs()
  xyz = sngl(st%s_copyd())
  write (90,"('sphere { <',2(F14.6,','),F14.6,'>,',F6.4,' }')") xyz(1:3), sphrad
  FZtmp => FZtmp%next
end do pointloop

write(rgbstring,"(F8.6,',',F8.6,',',F8.6)") (/ 0.0, 0.0, 1.0 /)
colorstring = 'material { texture { pigment { rgb <'//trim(rgbstring)//'> filter 0.95 }'

write (90,"(A)") trim(colorstring)
! this next line also closes the union of spheres...
write (90,"(' finish { diffuse 0.6, 0.6 brilliance 1.0 }  } } }')")
write (90,"(A)") 'background {   color rgb <0.9, 0.9, 0.9> }'
call PoV%closeFile()
if (POV%verbose) call Message%printMessage('PoVray rendering script stored in '//trim(povname)//'-st.pov')

!---------------------------------------------------------------------
! next we generate the PoVRay.ini file with the rendering instructions
open(unit=dataunit,file='povray.ini',status='unknown',form='formatted')
write(dataunit,"('Input_File_Name=',A)") trim(povfile)
write(dataunit,"('Output_File_Name=',A)") trim(povname)//'-st.png'

str = '+L'//trim(self%PVincludepath)
write(dataunit,"(A)") trim(str)
write(dataunit,"('+W',I4,' +H',I4)") 1024, 1024
write(dataunit,"('Initial_Clock=1')")
write(dataunit,"('Initial_Frame=1')")
write(dataunit,"('Final_Clock=1')")
write(dataunit,"('Final_Frame=1')")
write(dataunit,"('Work_Threads=6')")
close(unit=dataunit,status='keep')
if (POV%verbose) call Message%printMessage(' --> povray.ini file created ')

pvcmd = trim(self%PVexec)
inquire(file=trim(pvcmd),exist=fexists)
if (fexists.eqv..TRUE.) then
  if (POV%verbose) call Message%printMessage('Found PovRay command line executable; rendering frame')
  if (POV%verbose) call Message%printMessage('Executing '//trim(pvcmd)//' povray.ini >/dev/null 2>/dev/null')
  call system(trim(pvcmd)//' povray.ini >/dev/null 2>/dev/null')
end if


end subroutine renderdata_


end module mod_sampler