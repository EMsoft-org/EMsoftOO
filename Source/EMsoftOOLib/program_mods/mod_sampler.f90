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
  integer(kind=irg) :: norientations      ! number of orientations per dataset
  integer(kind=irg) :: pgnum              ! point group number (for Laue point group)
  character(fnlen)  :: hdfname            ! name of output HDF5 file
end type samplerNameListType

! class definition
type, public :: sampler_T
private 
  character(fnlen)            :: nmldeffile = 'EMsampler.nml'
  type(samplerNameListType)   :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: sampler_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: sampler => sampler_

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
integer(kind=irg) :: norientations      ! number of orientations per dataset
integer(kind=irg) :: pgnum              ! point group number (for Laue point group)
character(fnlen)  :: hdfname            ! name of output HDF5 file

namelist / EMsampler / kappa, dir1, dir2, dir3, dir4, norientations, pgnum, hdfname 

kappa = (/ 512.D0, 3456.D0 /)
dir1 = (/ 0.D0, 0.D0, 0.D0 /)
dir2 = (/ 0.D0, 0.D0, 0.D0 /) 
dir3 = (/ 0.D0, 0.D0, 0.D0 /) 
dir4 = (/ 0.D0, 0.D0, 0.D0 /) 
norientations = 5000
pgnum = 1
hdfname = 'undefined'

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
end if

self%nml%kappa = kappa
self%nml%dir1 = dir1
self%nml%dir2 = dir2
self%nml%dir3 = dir3
self%nml%dir4 = dir4
self%nml%norientations = norientations
self%nml%pgnum = pgnum  
self%nml%hdfname = hdfname

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

integer(kind=irg),parameter             :: n_int = 2, n_real = 9
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( enl => self%nml )

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! write all the single integers
io_int = (/ enl%norientations, enl%pgnum /)
intlist(1) = 'norientations'
intlist(2) = 'pgnum'
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
!! All results are stored in an HDF file.

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
type(QuaternionArray_T)                 :: qAR
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

real(kind=dbl),allocatable              :: vMFquatarray(:,:,:), WATquatarray(:,:,:), vMFkappahat(:), &
                                           WATkappahat(:), vMFmuhat(:,:), WATmuhat(:,:)
integer(kind=irg)                       :: setcnt, seed1, seed2, hdferr, i, j, k, io_int(2)
real(kind=dbl)                          :: rod(4), rodL, kappahat
real(kind=sgl)                          :: tstop
character(11)                           :: dstr
character(15)                           :: tstrb
character(15)                           :: tstre
character(fnlen)                        :: dataset, datagroupname, fname, attributename, HDF_FileVersion, nmldeffile
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
allocate( vMFkappahat(8), WATkappahat(8), vMFmuhat(4,8), WATmuhat(4,8) )

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

! we'll do vMF sampling first and generate 8 orientation data sets 
dictVMF = DirStat_T( DStype='VMF', PGnum=enl%pgnum )
call dictVMF%setNumEM(25)
call dictVMF%setNumIter(30)
seed1 = 54321
seed2 = 43514

setcnt = 1 
do i=1,2      ! loop over the kappa concentration parameter values 
  do j=1,4    ! loop over the mean quaternion directions
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
    write (*,*) muhat%get_quatd()
    write (*,*) meanquat(j)%get_quatd()
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
seed1 = 54321
seed2 = 43514

setcnt = 1 
do i=1,2      ! loop over the kappa concentration parameter values 
  do j=1,4    ! loop over the mean quaternion directions
    qAR = dictWAT%SampleDS( enl%norientations, seed1, meanquat(j), enl%kappa(i) )
! copy the orientations into the vMFquatarray for storage in the output HDF file
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
    write (*,*) muhat%get_quatd()
    write (*,*) meanquat(j)%get_quatd()
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



end module mod_sampler