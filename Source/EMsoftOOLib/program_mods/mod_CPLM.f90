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

module mod_CPLM
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! class definition for the EMCPLM program

use mod_kinds
use mod_global
use mod_MuellerCalculus

IMPLICIT NONE 

! namelist for the EMCPLM program
type, public :: CPLMNameListType
  integer(kind=irg)    :: phinum
  integer(kind=irg)    :: numpx
  integer(kind=irg)    :: numpy
  character(fnlen)     :: tiffprefix
  character(fnlen)     :: masterfile
  character(fnlen)     :: anglefile
  character(fnlen)     :: outputfile
end type CPLMNameListType

! class definition
type, public :: CPLM_T
private 
  character(fnlen)        :: nmldeffile = 'EMCPLM.nml'
  type(CPLMNameListType)  :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: CPLM_
  procedure, pass(self) :: setphinum_
  procedure, pass(self) :: getphinum_
  procedure, pass(self) :: setnumpx_
  procedure, pass(self) :: getnumpx_
  procedure, pass(self) :: setnumpy_
  procedure, pass(self) :: getnumpy_
  procedure, pass(self) :: settiffprefix_
  procedure, pass(self) :: gettiffprefix_
  procedure, pass(self) :: setmasterfile_
  procedure, pass(self) :: getmasterfile_
  procedure, pass(self) :: setanglefile_
  procedure, pass(self) :: getanglefile_
  procedure, pass(self) :: setoutputfile_
  procedure, pass(self) :: getoutputfile_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: CPLM => CPLM_
  generic, public :: setphinum => setphinum_
  generic, public :: getphinum => getphinum_
  generic, public :: setnumpx => setnumpx_
  generic, public :: getnumpx => getnumpx_
  generic, public :: setnumpy => setnumpy_
  generic, public :: getnumpy => getnumpy_
  generic, public :: settiffprefix => settiffprefix_
  generic, public :: gettiffprefix => gettiffprefix_
  generic, public :: setmasterfile => setmasterfile_
  generic, public :: getmasterfile => getmasterfile_
  generic, public :: setanglefile => setanglefile_
  generic, public :: getanglefile => getanglefile_
  generic, public :: setoutputfile => setoutputfile_
  generic, public :: getoutputfile => getoutputfile_

end type CPLM_T

! the constructor routine for this class 
interface CPLM_T
  module procedure CPLM_constructor
end interface CPLM_T

contains

!--------------------------------------------------------------------------
type(CPLM_T) function CPLM_constructor( nmlfile ) result(CPLM)
!DEC$ ATTRIBUTES DLLEXPORT :: CPLM_constructor
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! constructor for the CPLM_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call CPLM%readNameList(nmlfile)

end function CPLM_constructor

!--------------------------------------------------------------------------
subroutine CPLM_destructor(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: CPLM_destructor
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! destructor for the CPLM_T Class
 
IMPLICIT NONE

type(CPLM_T), INTENT(INOUT)  :: self 

call reportDestructor('CPLM_T')

end subroutine CPLM_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! read the namelist from an nml file for the CPLM_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(CPLM_T), INTENT(INOUT)         :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg)                    :: phinum
integer(kind=irg)                    :: numpx
integer(kind=irg)                    :: numpy
character(fnlen)                     :: tiffprefix
character(fnlen)                     :: masterfile
character(fnlen)                     :: anglefile
character(fnlen)                     :: outputfile

namelist / CPLMData / phinum, numpx, numpy, tiffprefix, masterfile, anglefile, outputfile

phinum = 36
numpx = 100
numpy = 100
tiffprefix = 'undefined'
masterfile = 'undefined'
anglefile = 'undefined'
outputfile = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=CPLMData)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(tiffprefix).eq.'undefined') then
  call Message%printError('readNameList:',' tiffprefix is undefined in '//nmlfile)
 end if
 if (trim(masterfile).eq.'undefined') then
  call Message%printError('readNameList:',' master output file name is undefined in '//nmlfile)
 end if
 if (trim(anglefile).eq.'undefined') then
  call Message%printError('readNameList:',' anglefile file name is undefined in '//nmlfile)
 end if
 if (trim(outputfile).eq.'undefined') then
  call Message%printError('readNameList:',' outputfile file name is undefined in '//nmlfile)
 end if
end if

self%nml%phinum = phinum
self%nml%numpx = numpx
self%nml%numpy = numpy
self%nml%tiffprefix = tiffprefix
self%nml%masterfile = masterfile
self%nml%anglefile = anglefile
self%nml%outputfile = outputfile

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! pass the namelist for the CPLM_T Class to the calling program

IMPLICIT NONE 

class(CPLM_T), INTENT(INOUT)          :: self
type(CPLMNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(CPLM_T), INTENT(INOUT)            :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 3, n_real = 9
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( nml => self%nml )

! create the group for this namelist
groupname = trim(HDFnames%get_NMLlist())
hdferr = HDF%createGroup(groupname)

! write all the single integers
io_int = (/ nml%phinum, nml%numpx, nml%numpy /)
intlist(1) = 'phinum'
intlist(2) = 'numpx' 
intlist(3) = 'numpy'
call HDF%writeNMLintegers(io_int, intlist, n_int)

dataset = 'tiffprefix'
line2(1) = nml%tiffprefix 
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList_: unable to create tiffprefix dataset', hdferr)

dataset = SC_masterfile
line2(1) = nml%masterfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList_: unable to create masterfile dataset', hdferr)

dataset = SC_anglefile
line2(1) = nml%anglefile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList_: unable to create anglefile dataset', hdferr)

dataset = SC_outputfile
line2(1) = nml%outputfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList_: unable to create outputfile dataset', hdferr)

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine setphinum_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setphinum_
!! author: MDG
!! version: 1.0
!! date: 08/10/23
!!
!! set phinum in the CPLM class

IMPLICIT NONE

class(CPLM_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%phinum = inp

end subroutine setphinum_

!--------------------------------------------------------------------------
function getphinum_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getphinum_
!! author: MDG
!! version: 1.0
!! date: 08/10/23
!!
!! get phinum from the CPLM class

IMPLICIT NONE

class(CPLM_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%phinum

end function getphinum_

!--------------------------------------------------------------------------
subroutine setnumpx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnumpx_
!! author: MDG
!! version: 1.0
!! date: 08/10/23
!!
!! set numpx in the CPLM class

IMPLICIT NONE

class(CPLM_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%numpx = inp

end subroutine setnumpx_

!--------------------------------------------------------------------------
function getnumpx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnumpx_
!! author: MDG
!! version: 1.0
!! date: 08/10/23
!!
!! get numpx from the CPLM class

IMPLICIT NONE

class(CPLM_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%numpx

end function getnumpx_

!--------------------------------------------------------------------------
subroutine setnumpy_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnumpy_
!! author: MDG
!! version: 1.0
!! date: 08/10/23
!!
!! set numpy in the CPLM class

IMPLICIT NONE

class(CPLM_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%numpy = inp

end subroutine setnumpy_

!--------------------------------------------------------------------------
function getnumpy_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnumpy_
!! author: MDG
!! version: 1.0
!! date: 08/10/23
!!
!! get numpy from the CPLM class

IMPLICIT NONE

class(CPLM_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%numpy

end function getnumpy_

!--------------------------------------------------------------------------
subroutine settiffprefix_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: settiffprefix_
!! author: MDG
!! version: 1.0
!! date: 08/10/23
!!
!! set tiffprefix in the CPLM class

IMPLICIT NONE

class(CPLM_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%tiffprefix = trim(inp)

end subroutine settiffprefix_

!--------------------------------------------------------------------------
function gettiffprefix_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: gettiffprefix_
!! author: MDG
!! version: 1.0
!! date: 08/10/23
!!
!! get tiffprefix from the CPLM class

IMPLICIT NONE

class(CPLM_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%tiffprefix)

end function gettiffprefix_

!--------------------------------------------------------------------------
subroutine setmasterfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setmasterfile_
!! author: MDG
!! version: 1.0
!! date: 08/10/23
!!
!! set masterfile in the CPLM class

IMPLICIT NONE

class(CPLM_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%masterfile = trim(inp)

end subroutine setmasterfile_

!--------------------------------------------------------------------------
function getmasterfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getmasterfile_
!! author: MDG
!! version: 1.0
!! date: 08/10/23
!!
!! get masterfile from the CPLM class

IMPLICIT NONE

class(CPLM_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%masterfile)

end function getmasterfile_

!--------------------------------------------------------------------------
subroutine setanglefile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setanglefile_
!! author: MDG
!! version: 1.0
!! date: 08/10/23
!!
!! set anglefile in the CPLM class

IMPLICIT NONE

class(CPLM_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%anglefile = trim(inp)

end subroutine setanglefile_

!--------------------------------------------------------------------------
function getanglefile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getanglefile_
!! author: MDG
!! version: 1.0
!! date: 08/10/23
!!
!! get anglefile from the CPLM class

IMPLICIT NONE

class(CPLM_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%anglefile)

end function getanglefile_

!--------------------------------------------------------------------------
subroutine setoutputfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setoutputfile_
!! author: MDG
!! version: 1.0
!! date: 08/10/23
!!
!! set outputfile in the CPLM class

IMPLICIT NONE

class(CPLM_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%outputfile = trim(inp)

end subroutine setoutputfile_

!--------------------------------------------------------------------------
function getoutputfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getoutputfile_
!! author: MDG
!! version: 1.0
!! date: 08/10/23
!!
!! get outputfile from the CPLM class

IMPLICIT NONE

class(CPLM_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%outputfile)

end function getoutputfile_

!--------------------------------------------------------------------------
subroutine CPLM_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: CPLM_
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! perform the computations

use mod_EMsoft
use mod_HDFnames
use mod_HDFsupport
use HDF5
use mod_MuellerCalculus
use mod_CPLMmaster
use mod_so3
use mod_crystallography
use mod_symmetry
use stringconstants
use mod_rotations
use mod_quaternions
use mod_memory
use mod_Lambert
use mod_timing

IMPLICIT NONE 

class(CPLM_T), INTENT(INOUT)            :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

type(IO_T)                              :: Message 
type(CPLMmaster_T)                      :: CPLMmaster
type(HDFnames_T)                        :: saveHDFnames
type(so3_T)                             :: SO
type(Cell_T)                            :: cell
type(SpaceGroup_T)                      :: SG
type(HDF_T)                             :: HDF
type(QuaternionArray_T)                 :: qAR
type(Quaternion_T)                      :: quat
type(a_T)                               :: ax
type(q_T)                               :: qu
type(memory_T)                          :: mem
type(StokesVectorType)                  :: SVin, SV, SVout
type(MuellerCalculus_T)                 :: MC
type(MuellerMatrixType)                 :: MMsample, MMchain
type(Timing_T)                          :: timer
  
character(fnlen)                        :: fname, oname, descriptor, datafile, dataset, groupname, attributename, &
                                           datagroupname, HDF_FileVersion 
logical                                 :: f_exists, g_exists, overwrite = .TRUE.
integer(kind=irg)                       :: pgnum, hdferr, npx, numpoints, i, j, k, nix, niy, nixp, niyp, io_int(1)
real(kind=dbl),allocatable              :: intensities(:,:), qrot(:,:), images(:,:,:)
real(kind=dbl)                          :: vc(3), vr(3), dc(3), phistepsize, scl, dx, dy, dxm, dym
real(kind=sgl)                          :: tstop
character(11)                           :: dstr
character(15)                           :: tstrb
character(15)                           :: tstre
character(fnlen,kind=c_char)            :: line2(1)

call openFortranHDFInterface()

HDF = HDF_T() 
mem = memory_T() 

! initialize the timing routines
timer = Timing_T()
tstrb = timer%getTimeString()  
call timer%Time_tick(1)

associate(nml => self%nml)

!-------------------------------------
! read the CPLM master pattern from the masterfile
fname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(nml%masterfile)
inquire(file=trim(fname), exist=f_exists)
if (f_exists.eqv..FALSE.) then
  call Message%printError('CPLM',' masterfile does not exist')
end if

saveHDFnames = HDFnames 
call HDFnames%set_ProgramData(SC_CPLMmaster)
call HDFnames%set_NMLlist(SC_CPLMmasterNameList)
call HDFnames%set_NMLfilename(SC_CPLMmasterNML)
! call HDFnames%get_AllNames()
call CPLMmaster%readCPLMmasterfile(EMsoft, HDF, HDFnames, fname, hdferr, getMaster=.TRUE.)
HDFnames = saveHDFnames
! call HDFnames%get_AllNames()
npx = CPLMmaster%getnpx()

!-------------------------------------
! get the point group number from the xtal file 
call cell%getCrystalData(CPLMmaster%getxtalname(), SG, EMsoft)

! initialize the SO3 class
SO = so3_T( SG%getPGnumber(), zerolist='FZ' )

! read all the orientations from the anglefile
oname = EMsoft%generateFilePath('EMdatapathname',nml%anglefile)
call Message%printMessage(' Reading orientations from file '//trim(oname))
call SO%getOrientationsfromFile( oname, listN=10 )
numpoints = SO%getListCount( 'FZ' )

! and convert the list of orientations to a quaternion array
call SO%listtoQuaternionArray( qAR )

! we no longer need the linked list 
call SO%delete_FZlist()
call Message%printMessage('    ---> completed reading orientations ')

!-------------------------------------
! allocate array for intensity results
call mem%alloc(intensities, (/ nml%phinum, numpoints /), 'intensities', initval=0.D0 )
call mem%alloc(images, (/ nml%phinum, nml%numpx, nml%numpy /), 'images', initval=0.D0 )

!===============================================
! for each orientation in the list, generate phinum sets of direction cosines for the rotating [001] axis
! and for each of those, compute the Muller matrix by interpolation from the master pattern.  
! Then multiply the Mueller matrix with the analyzer matrix and the linearly polarized Stokes vector 
! and keep the first element of the resulting vector as the intensity.
vc = (/ 0.D0, 0.D0, 1.D0 /)  ! the original c-axis

! set up the sample rotations around the ND axis
call mem%alloc(qrot, (/ 4, nml%phinum /), 'qrot')
phistepsize = 2.D0*cPi/dble(nml%phinum)
do i=1,nml%phinum
  if (i.le.nml%phinum/2) then 
    ax = a_T( adinp = (/0.D0, 0.D0, 1.D0, dble(i-1) * phistepsize /) )
  else ! flip the rotation axis and complement the angle
    ax = a_T( adinp = (/0.D0, 0.D0,-1.D0, 2.D0*cPi - dble(i-1) * phistepsize /) )
  end if 
  qu = ax%aq()
  qrot(1:4,i) = qu%q_copyd()
end do

descriptor = 'incident on sample'

! define the input Stokes vector
SVin%S = (/ 1.D0, 0.D0, 0.D0, 0.D0 /)
SVin%descriptor = 'input Stokes vector'

call MC%print_StokesVector( SVin ) 

SVin = MC%propagateStokesVector(MC%concatenateMuellerMatrices(MC%get_basicMuellerMatrix(1), &
        MC%get_basicMuellerMatrix(9)), SVin, descriptor)
!SVin = MC%propagateStokesVector(MC%get_basicMuellerMatrix(3), SVin, descriptor)

call MC%print_StokesVector( SVin ) 

! define the analyzer optics
MMchain = MC%concatenateMuellerMatrices(MC%get_basicMuellerMatrix(9), &
        MC%rotate_MuellerMatrix(MC%get_basicMuellerMatrix(2), 1.D0*cPi/180.D0) )
!MMchain = MC%rotate_MuellerMatrix(MC%get_basicMuellerMatrix(4), 5.D0*cPi/180.D0)

call MC%print_MuellerMatrix( MMchain )

! loop over all orientations
scl = dble(npx)
do i = 1,numpoints
  quat = qAR%getQuatfromArray( i )
  vr = quat%quat_Lp( vc )
!  write(53,"(2(F12.8,','),F12.8)") vr
  if (vr(3).ne.1.D0) then 
! next we use these  to get the Mueller matrix by interpolation for each of the sample rotation steps
    do j = 1, nml%phinum
      quat = Quaternion_T( qd = qrot(1:4,j) ) 
      dc = quat%quat_Lp( vr )
      if (dc(3).le.0.D0) dc = -dc
! convert these direction cosines to coordinates in the Rosca-Lambert projection
      call LambertgetInterpolation(dc, scl, npx, npx, nix, niy, nixp, niyp, dx, dy, dxm, dym)

      MMsample%M = CPLMmaster%MPNH(1:4,1:4,nix,niy) * dxm * dym + &
                   CPLMmaster%MPNH(1:4,1:4,nixp,niy) * dx * dym + &
                   CPLMmaster%MPNH(1:4,1:4,nix,niyp) * dxm * dy + &
                   CPLMmaster%MPNH(1:4,1:4,nixp,niyp) * dx * dy 
! next we apply this to the incident Stokes vector 
      SV = MC%propagateStokesVector(MMsample, SVin, descriptor)
! and we propagate to the detector
      SVout = MC%propagateStokesVector(MMchain, SV, descriptor)
      intensities(j,i) = SVout%S(0)
    end do
  end if
end do

! convert intensities to images array 
do j = 1,nml%numpy 
  do i = 1,nml%numpx 
    k = (j-1)*nml%numpx + i 
    images(1:nml%phinum, i, j) = intensities(1:nml%phinum,k)
  end do 
end do


call timer%Time_tock(1)
tstop = timer%getInterval(1)
call timer%Time_reset(1)

call timer%makeTimeStamp()
dstr = timer%getDateString()
tstre = timer%getTimeString()

io_int(1) = tstop
call Message%WriteValue('Total execution time [s] = ',io_int,1)

! Create a new file using the default properties.
datafile = trim(EMsoft%generateFilePath('EMdatapathname',nml%outputfile))
hdferr =  HDF%createFile(datafile)

! write the EMheader to the file
  datagroupname = trim(HDFnames%get_ProgramData())
  call HDF%writeEMheader(EMsoft, dstr, tstrb, tstre, progname, datagroupname)

! open or create a namelist group to write all the namelist files into
  hdferr = HDF%createGroup(HDFnames%get_NMLfiles())

! read the text file and write the array to the file
  dataset = HDFnames%get_NMLfilename()
  hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%nmldeffile)

! leave this group
  call HDF%pop()

! create a namelist group to write all the namelist files into
  hdferr = HDF%createGroup(HDFnames%get_NMLparameters())
  call self%writeHDFNameList(HDF, HDFnames)

  ! leave this group
  call HDF%pop()
  call HDF%pop()

! then the remainder of the data in a EMData group
  hdferr = HDF%createGroup(HDFnames%get_EMData())

! create the CPLMmaster group and add a HDF_FileVersion attribbute to it
  hdferr = HDF%createGroup(HDFnames%get_ProgramData())
  HDF_FileVersion = '4.0'
  HDF_FileVersion = cstringify(HDF_FileVersion)
  attributename = SC_HDFFileVersion
  hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)
  
  dataset = SC_Manufacturer
  line2(1) = 'EMsoft'
  line2(1) = cstringify(line2(1))
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)

  ! and start writing the data arrays
  dataset = 'intensities'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetDoubleArray(dataset, intensities, nml%phinum, numpoints, overwrite)
  else
    hdferr = HDF%writeDatasetDoubleArray(dataset, intensities, nml%phinum, numpoints)
  end if

  dataset = 'images'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetDoubleArray(dataset, images, nml%phinum, nml%numpx, nml%numpy, overwrite)
  else
    hdferr = HDF%writeDatasetDoubleArray(dataset, images, nml%phinum, nml%numpx, nml%numpy)
  end if


  call HDF%popall()


! close the fortran HDF interface
call closeFortranHDFInterface()

end associate

end subroutine CPLM_



end module mod_CPLM