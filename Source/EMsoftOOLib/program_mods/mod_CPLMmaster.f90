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

module mod_CPLMmaster
    !! author: MDG 
    !! version: 1.0 
    !! date: 08/10/23
    !!
    !! class definition for the EMCPLMmaster program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the CPLMmaster program
type, public :: CPLMmasterNameListType
  integer(kind=irg)       :: npx
  integer(kind=irg)       :: nthreads
  real(kind=sgl)          :: eps1Re
  real(kind=sgl)          :: eps1Im
  real(kind=sgl)          :: eps2Re
  real(kind=sgl)          :: eps2Im
  real(kind=sgl)          :: wl
  real(kind=sgl)          :: theta
  logical                 :: normalize
  character(3)            :: Notify
  character(fnlen)        :: xtalname
  character(fnlen)        :: masterfile
end type CPLMmasterNameListType

! class definition
type, public :: CPLMmaster_T
private 
  character(fnlen)              :: nmldeffile = 'EMCPLMmaster.nml'
  type(CPLMmasterNameListType)  :: nml 
  real(kind=dbl),allocatable, public    :: MPNH(:,:,:,:)

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: readCPLMmasterfile_
  procedure, pass(self) :: CPLMmaster_
  procedure, pass(self) :: setnpx_
  procedure, pass(self) :: getnpx_
  procedure, pass(self) :: setnthreads_
  procedure, pass(self) :: getnthreads_
  procedure, pass(self) :: seteps1Re_
  procedure, pass(self) :: geteps1Re_
  procedure, pass(self) :: seteps1Im_
  procedure, pass(self) :: geteps1Im_
  procedure, pass(self) :: seteps2Re_
  procedure, pass(self) :: geteps2Re_
  procedure, pass(self) :: seteps2Im_
  procedure, pass(self) :: geteps2Im_
  procedure, pass(self) :: setwl_
  procedure, pass(self) :: getwl_
  procedure, pass(self) :: settheta_
  procedure, pass(self) :: gettheta_
  procedure, pass(self) :: setnormalize_
  procedure, pass(self) :: getnormalize_
  procedure, pass(self) :: setNotify_
  procedure, pass(self) :: getNotify_
  procedure, pass(self) :: setxtalname_
  procedure, pass(self) :: getxtalname_
  procedure, pass(self) :: setmasterfile_
  procedure, pass(self) :: getmasterfile_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: readCPLMmasterfile => readCPLMmasterfile_
  generic, public :: CPLMmaster => CPLMmaster_
  generic, public :: setnpx => setnpx_
  generic, public :: getnpx => getnpx_
  generic, public :: setnthreads => setnthreads_
  generic, public :: getnthreads => getnthreads_
  generic, public :: seteps1Re => seteps1Re_
  generic, public :: geteps1Re => geteps1Re_
  generic, public :: seteps1Im => seteps1Im_
  generic, public :: geteps1Im => geteps1Im_
  generic, public :: seteps2Re => seteps2Re_
  generic, public :: geteps2Re => geteps2Re_
  generic, public :: seteps2Im => seteps2Im_
  generic, public :: geteps2Im => geteps2Im_
  generic, public :: setwl => setwl_
  generic, public :: getwl => getwl_
  generic, public :: settheta => settheta_
  generic, public :: gettheta => gettheta_
  generic, public :: setnormalize => setnormalize_
  generic, public :: getnormalize => getnormalize_
  generic, public :: setNotify => setNotify_
  generic, public :: getNotify => getNotify_
  generic, public :: setxtalname => setxtalname_
  generic, public :: getxtalname => getxtalname_
  generic, public :: setmasterfile => setmasterfile_
  generic, public :: getmasterfile => getmasterfile_

end type CPLMmaster_T

! the constructor routine for this class 
interface CPLMmaster_T
  module procedure CPLMmaster_constructor
end interface CPLMmaster_T

contains

!--------------------------------------------------------------------------
type(CPLMmaster_T) function CPLMmaster_constructor( nmlfile ) result(CPLMmaster)
!DEC$ ATTRIBUTES DLLEXPORT :: CPLMmaster_constructor
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! constructor for the CPLMmaster_T Class; optionally reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

if (present(nmlfile)) then 
  call CPLMmaster%readNameList(nmlfile)
end if

end function CPLMmaster_constructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! read the namelist from an nml file for the CPLMmaster_T Class 

use mod_io 

IMPLICIT NONE 

class(CPLMmaster_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg)                    :: npx
integer(kind=irg)                    :: nthreads
real(kind=sgl)                       :: eps1Re
real(kind=sgl)                       :: eps1Im
real(kind=sgl)                       :: eps2Re
real(kind=sgl)                       :: eps2Im
real(kind=sgl)                       :: wl
real(kind=sgl)                       :: theta 
logical                              :: normalize
character(3)                         :: Notify
character(fnlen)                     :: xtalname
character(fnlen)                     :: masterfile

! define the IO namelist to facilitate passing variables to the program.
namelist  / CPLMMasterData / npx, nthreads, eps1Re, eps1Im, eps2Re, eps2Im, wl, theta, &
          Notify, xtalname, masterfile, normalize  

xtalname = 'undefined'
theta = 0.0 
wl = 750.0
eps1Re =  1.0
eps1Im =  0.0
eps2Re =  1.0
eps2Im =  0.0
npx = 360
normalize = .FALSE.
nthreads = 1
masterfile = 'undefined'
Notify = 'Off' 

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=CPLMMasterData)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(xtalname).eq.'undefined') then
  call Message%printError('readNameList:',' structure file name is undefined in '//nmlfile)
 end if
 if (trim(masterfile).eq.'undefined') then
  call Message%printError('readNameList:',' master output file name is undefined in '//nmlfile)
 end if
end if

self%nml%npx = npx
self%nml%nthreads = nthreads
self%nml%eps1Re = eps1Re
self%nml%eps1Im = eps1Im
self%nml%eps2Re = eps2Re
self%nml%eps2Im = eps2Im
self%nml%wl = wl
self%nml%theta = theta
self%nml%normalize = normalize
self%nml%Notify = Notify 
self%nml%xtalname = xtalname
self%nml%masterfile = masterfile

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! pass the namelist for the CPLMmaster_T Class to the calling program

IMPLICIT NONE 

class(CPLMmaster_T), INTENT(INOUT)          :: self
type(CPLMmasterNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)        :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 3, n_real = 6
integer(kind=irg)                       :: hdferr,  io_int(n_int), nm
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( omnl => self%nml )

  ! create the group for this namelist
groupname = trim(HDFnames%get_NMLlist())
hdferr = HDF%createGroup(groupname)

! convert normalize parameter to integer
nm = 0
if (omnl%normalize.eqv..TRUE.) nm = 1

! write all the single integers
io_int = (/ omnl%npx, omnl%nthreads, nm /)
intlist(1) = 'npx'
intlist(2) = 'nthreads' 
intlist(3) = 'normalize'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ omnl%eps1Re, omnl%eps1Im, omnl%eps2Re, omnl%eps2Im, omnl%wl, omnl%theta /)
reallist(1) = 'eps1Re'
reallist(2) = 'eps1Im'
reallist(3) = 'eps2Re'
reallist(4) = 'eps2Im'
reallist(5) = 'wl' 
reallist(6) = 'theta' 
call HDF%writeNMLreals(io_real, reallist, n_real)

dataset = SC_xtalname
line2(1) = omnl%xtalname 
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteCPLMNameList: unable to create maskpattern dataset', hdferr)

dataset = SC_masterfile
line2(1) = omnl%masterfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteCPLMNameList: unable to create maskpattern dataset', hdferr)


end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
recursive subroutine readCPLMmasterfile_(self, EMsoft, HDF, HDFnames, mpfile, hdferr, getMaster)
!DEC$ ATTRIBUTES DLLEXPORT :: readCPLMmasterfile_
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! read data from CPLM master HDF file

use mod_EMsoft
use HDF5
use mod_HDFsupport
use mod_HDFnames
use stringconstants 
use mod_io

use ISO_C_BINDING

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)        :: self 
type(EMsoft_T), INTENT(INOUT)             :: EMsoft
type(HDF_T), INTENT(INOUT)                :: HDF
type(HDFnames_T), INTENT(INOUT)           :: HDFnames
character(fnlen),INTENT(IN)               :: mpfile
integer(kind=irg),INTENT(OUT)             :: hdferr
logical,INTENT(IN),OPTIONAL               :: getMaster

type(IO_T)                                :: Message 

character(fnlen)                          :: groupname, dataset, tmpnmlname 
logical                                   :: g_exists, stat 
integer(kind=irg)                         :: nlines, i, j, npx
integer(HSIZE_T)                          :: sz(1), dims4(4)
real(kind=dbl),allocatable                :: master(:,:,:,:)

character(fnlen, KIND=c_char),allocatable,TARGET    :: stringarray(:)

! is this a proper HDF5 file ?
call h5fis_hdf5_f(trim(mpfile), stat, hdferr)

! open the file using readonly access
hdferr =  HDF%openFile(mpfile, readonly = .TRUE.)

groupname = trim(HDFnames%get_NMLfiles())
write (*,*) 'attempting to open group '//trim(groupname)
hdferr = HDF%openGroup(groupname)

dataset = trim(HDFnames%get_NMLfilename())
write (*,*) 'attempting to open data set '//trim(dataset)
  call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)

if (g_exists.eqv..TRUE.) then
!====================================
! read all NMLparameters group datasets by writing the NMLfiles string array to a
! temporary file and calling the regular nml read routine. This should replace a
! large number of individual read routines, so it simplifies the code ...
!====================================
    call HDF%readdatasetstringarray(dataset, nlines, hdferr, stringarray)
    sz = shape(stringarray)
    tmpnmlname = trim(EMsoft%generateFilePath('EMtmppathname'))//'tmp.nml'
    open(unit=65,file=trim(tmpnmlname),status='unknown',form='formatted')
    write (65,"(A)") '&CPLMMasterData'
    do i=2,sz(1)
      write (65,"(A)") trim(stringarray(i))
    end do
    close(unit=65,status='keep')
    call self%readNameList(tmpnmlname)
! delete the tmp file
    open(unit=65,file=trim(tmpnmlname),status='unknown',form='formatted')
    close(unit=65,status='delete')
    call HDF%pop()
else
    call Message%printError('readCPLMmasterfile_','could not read namelist')
end if

npx = self%nml%npx

groupname = SC_EMData
hdferr = HDF%openGroup(groupname)

groupname = trim(HDFnames%get_ProgramData())
hdferr = HDF%openGroup(groupname)

if (present(getMaster)) then 
  if (getMaster.eqv..TRUE.) then
    dataset = SC_CPLMmasterLPNH
    call HDF%readDatasetDoubleArray(dataset, dims4, hdferr, master)
    allocate(self%MPNH(4, 4, -npx:npx, -npx:npx))
    do i=-npx,npx
      do j=-npx,npx
        self%MPNH(1:4,1:4,i,j) = master(1:4,1:4,i+1+npx,j+1+npx)
      end do
    end do
  end if 
end if 

call Message%printMessage(' read data from'//trim(mpfile))

! and close everything
call HDF%popall()

end subroutine readCPLMmasterfile_


!--------------------------------------------------------------------------
subroutine setnpx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnpx_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! set npx in the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%npx = inp

end subroutine setnpx_

!--------------------------------------------------------------------------
function getnpx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnpx_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! get npx from the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%npx

end function getnpx_

!--------------------------------------------------------------------------
subroutine setnthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnthreads_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! set nthreads in the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%nthreads = inp

end subroutine setnthreads_

!--------------------------------------------------------------------------
function getnthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnthreads_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! get nthreads from the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%nthreads

end function getnthreads_

!--------------------------------------------------------------------------
subroutine seteps1Re_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: seteps1Re_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! set eps1Re in the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%eps1Re = inp

end subroutine seteps1Re_

!--------------------------------------------------------------------------
function geteps1Re_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: geteps1Re_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! get eps1Re from the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%eps1Re

end function geteps1Re_

!--------------------------------------------------------------------------
subroutine seteps1Im_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: seteps1Im_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! set eps1Im in the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%eps1Im = inp

end subroutine seteps1Im_

!--------------------------------------------------------------------------
function geteps1Im_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: geteps1Im_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! get eps1Im from the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%eps1Im

end function geteps1Im_

!--------------------------------------------------------------------------
subroutine seteps2Re_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: seteps2Re_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! set eps2Re in the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%eps2Re = inp

end subroutine seteps2Re_

!--------------------------------------------------------------------------
function geteps2Re_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: geteps2Re_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! get eps2Re from the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%eps2Re

end function geteps2Re_

!--------------------------------------------------------------------------
subroutine seteps2Im_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: seteps2Im_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! set eps2Im in the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%eps2Im = inp

end subroutine seteps2Im_

!--------------------------------------------------------------------------
function geteps2Im_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: geteps2Im_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! get eps2Im from the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%eps2Im

end function geteps2Im_

!--------------------------------------------------------------------------
subroutine setwl_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setwl_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! set wl in the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%wl = inp

end subroutine setwl_

!--------------------------------------------------------------------------
function getwl_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getwl_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! get wl from the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%wl

end function getwl_

!--------------------------------------------------------------------------
subroutine settheta_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: settheta_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! set theta in the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%theta = inp

end subroutine settheta_

!--------------------------------------------------------------------------
function gettheta_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: gettheta_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! get theta from the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%theta

end function gettheta_

!--------------------------------------------------------------------------
subroutine setnormalize_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnormalize_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! set normalize in the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
logical, INTENT(IN)       :: inp

self%nml%normalize = inp

end subroutine setnormalize_

!--------------------------------------------------------------------------
function getnormalize_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnormalize_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! get normalize from the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
logical                   :: out

out = self%nml%normalize

end function getnormalize_

!--------------------------------------------------------------------------
subroutine setNotify_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setNotify_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! set Notify in the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
character(3), INTENT(IN)       :: inp

self%nml%Notify = trim(inp)

end subroutine setNotify_

!--------------------------------------------------------------------------
function getNotify_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getNotify_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! get Notify from the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
character(3)                   :: out

out = trim(self%nml%Notify)

end function getNotify_

!--------------------------------------------------------------------------
subroutine setxtalname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setxtalname_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! set xtalname in the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%xtalname = trim(inp)

end subroutine setxtalname_

!--------------------------------------------------------------------------
function getxtalname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getxtalname_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! get xtalname from the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%xtalname)

end function getxtalname_

!--------------------------------------------------------------------------
subroutine setmasterfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setmasterfile_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! set masterfile in the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%masterfile = trim(inp)

end subroutine setmasterfile_

!--------------------------------------------------------------------------
function getmasterfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getmasterfile_
!! author: MDG
!! version: 1.0
!! date: 08/13/23
!!
!! get masterfile from the CPLMmaster_T class

IMPLICIT NONE

class(CPLMmaster_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%masterfile)

end function getmasterfile_

!--------------------------------------------------------------------------
subroutine CPLMmaster_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: CPLMmaster_  
!! author: MDG 
!! version: 1.0 
!! date: 08/10/23
!!
!! perform the computations

use mod_EMsoft
use mod_initializers
use mod_symmetry
use mod_crystallography
use mod_Lambert
use mod_io
use mod_diffraction
use HDF5
use mod_HDFsupport
use mod_HDFnames
use iso_c_binding
use mod_MuellerCalculus
use mod_timing
use mod_math
use stringconstants
use mod_memory

IMPLICIT NONE 

class(CPLMmaster_T), INTENT(INOUT)      :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 

type(Cell_T)                            :: cell
type(DynType)                           :: Dyn
type(Diffraction_T)                     :: Diff
type(HDF_T)                             :: HDF
type(HDFnames_T)                        :: HDFnames
type(IO_T)                              :: Message
type(Lambert_T)                         :: L
type(MuellerCalculus_T)                 :: MC
type(SpaceGroup_T)                      :: SG
type(Timing_T)                          :: timer
type(memory_T)                          :: mem

complex(kind=dbl)                       :: eps(2)    ! principal dielectric constants
complex(kind=dbl)                       :: rvals(4)  ! reflectivity coefficients
real(kind=dbl)                          :: lambda, xy(2), dc(3), edge, MM(4,4), xyz(3), Radius, theta, phii, tmp
real(kind=sgl)                          :: dmin, EkeV, tstop, mi, ma

real(kind=dbl),allocatable              :: LPNH(:,:,:,:)
real(kind=sgl),allocatable              :: SPNH(:,:,:,:)
real(kind=sgl),allocatable              :: combo(:,:), combo_unscaled(:,:)

integer(kind=irg)                       :: hdferr, i, ii, jj, io_int(1), nx, ny, ierr

integer(HSIZE_T), dimension(1:3)        :: hdims, offset
integer(HSIZE_T)                        :: dims3(3)
character(fnlen,kind=c_char)            :: line2(1)
character(fnlen)                        :: groupname, dataset, attributename, datagroupname, HDF_FileVersion, fname
character(11)                           :: dstr
character(15)                           :: tstrb
character(15)                           :: tstre
character(fnlen)                        :: datafile
logical                                 :: overwrite = .TRUE., insert = .TRUE., g_exists
logical                                 :: verbose

call openFortranHDFInterface()
HDF = HDF_T()

! set the HDF group names for this program
HDFnames = HDFnames_T()

! simplify the notation a little
associate( omnl => self%nml )

! initialize the timing routines
timer = Timing_T()
tstrb = timer%getTimeString()  

! initialize the memory class
mem = memory_T()

!=============================================
!=============================================
! crystallography section
dmin = 0.1
EkeV = 1.0

call cell%setFileName(omnl%xtalname)
call Diff%setrlpmethod('WK')

call Diff%setV(dble(EkeV))
call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, dmin, verbose, useHDF=HDF)

! illumination wave length [m]
lambda = omnl%wl * 1.0D-9 

! principal dielectric constants as complex numbers
eps = (/ cmplx(dble(omnl%eps1Re),dble(omnl%eps1Im)), cmplx(dble(omnl%eps2Re),dble(omnl%eps2Im)) /)

!=============================================
!=============================================
! allocate the output arrays
nx = 2*omnl%npx+1
ny = 2*omnl%npx+1

call mem%alloc(LPNH, (/ 4,4,omnl%npx,omnl%npx /), 'LPNH', initval = 0.D0, startdims = (/ 1, 1, -omnl%npx,-omnl%npx /) )
call mem%alloc(SPNH, (/ 4,4,omnl%npx,omnl%npx /), 'SPNH', initval = 0.0, startdims = (/ 1, 1, -omnl%npx,-omnl%npx /) )

!===============================================
!  perform the actual computation of the master Mueller Matrix patterns
!===============================================
call timer%Time_tick(1)
edge = 1.D0 / dble(omnl%npx)
call Message%printMessage('Starting master Mueller Matrix computation')

do ii=-omnl%npx,omnl%npx
  do jj=-omnl%npx,omnl%npx
    ! determine the spherical direction for this point
    L = Lambert_T( xyd = (/ dble(ii), dble(jj) /) * edge )
    ierr = L%LambertSquareToSphere(dc)
    rvals = MC%get_UniaxialReflectivities(lambda, eps, 1.D0, dc, dble(omnl%theta))
    MM = MC%get_SampleMuellerMatrix(rvals)
    ! normalize against M[1,1] ; make this optional !!!
    if ((MM(1,1).ne.0.D0).and.(omnl%normalize.eqv..TRUE.)) MM = MM / MM(1,1)
    LPNH(1:4,1:4,ii,jj) = MM(1:4,1:4)
  end do
end do

! convert the arrays to stereographic projections
call Message%printMessage('Starting stereographic projection conversion')
Radius = 1.D0
do ii=-omnl%npx,omnl%npx 
  do jj=-omnl%npx,omnl%npx 
    L = Lambert_T( xyd = (/ dble(ii), dble(jj) /) * edge )
    ierr = L%StereoGraphicInverse( xyz, Radius )
    xyz = xyz/vecnorm(xyz)
    if (ierr.ne.0) then 
      SPNH(1:4,1:4,ii,jj) = 0.0
    else
      SPNH(1:4,1:4,ii,jj) = sngl(InterpolateLambert(xyz, LPNH, omnl%npx))
    end if
  end do
end do

! finally, place the SPs all together in a single array on a common gray scale 
call mem%alloc( combo, (/ 4*nx, 4*ny /), 'combo', initval = 0.0 )
do ii=1,4
  do jj=1,4
    combo((jj-1)*nx+1:jj*nx,(ii-1)*ny+1:ii*ny) = transpose(SPNH(ii,jj,:,:))
  end do 
end do 

! then, place them all together in a single array scaled individually between 0.0 and 1.0
call mem%alloc( combo_unscaled, (/ 4*nx, 4*ny /), 'combo_unscaled', initval = 0.0 )
do ii=1,4
  do jj=1,4
    mi = minval(SPNH(ii,jj,:,:))
    ma = maxval(SPNH(ii,jj,:,:))
    if ((ma-mi).ne.0.0) then 
      combo_unscaled((jj-1)*nx+1:jj*nx,(ii-1)*ny+1:ii*ny) = (transpose(SPNH(ii,jj,:,:))-mi)/(ma-mi)
    end if 
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

!=============================================
! create or update the HDF5 output file
!=============================================
call HDFnames%set_ProgramData(SC_CPLMmaster)
call HDFnames%set_NMLlist(SC_CPLMmasterNameList)
call HDFnames%set_NMLfilename(SC_CPLMmasterNML)

! Create a new file using the default properties.
datafile = trim(EMsoft%generateFilePath('EMdatapathname',omnl%masterfile))
hdferr =  HDF%createFile(datafile)

! write the EMheader to the file
  datagroupname = trim(HDFnames%get_ProgramData())
  call HDF%writeEMheader(EMsoft,dstr, tstrb, tstre, progname, datagroupname)

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
  dataset = SC_CPLMmasterLPNH
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetDoubleArray(dataset, LPNH, 4, 4, nx, ny, overwrite)
  else
    hdferr = HDF%writeDatasetDoubleArray(dataset, LPNH, 4, 4, nx, ny)
  end if

  ! and also the stereographic projection version of these arrays

  dataset = SC_CPLMmasterSPNH
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetFloatArray(dataset, SPNH, 4, 4, nx, ny, overwrite)
  else
    hdferr = HDF%writeDatasetFloatArray(dataset, SPNH, 4, 4, nx, ny)
  end if

  ! add the large output arrays
  dataset = 'combo'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetFloatArray(dataset, combo, 4*nx, 4*ny, overwrite)
  else
    hdferr = HDF%writeDatasetFloatArray(dataset, combo, 4*nx, 4*ny)
  end if

  dataset = 'combo_unscaled'
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetFloatArray(dataset, combo_unscaled, 4*nx, 4*ny, overwrite)
  else
    hdferr = HDF%writeDatasetFloatArray(dataset, combo_unscaled, 4*nx, 4*ny)
  end if

  call HDF%popall()

  call Message%printMessage(' Final data stored in file '//trim(omnl%masterfile), frm = "(A/)")  
end associate

call mem%dealloc(LPNH, 'LPNH')
call mem%dealloc(SPNH, 'SPNH')
call mem%dealloc(combo, 'combo')

end subroutine CPLMmaster_

end module mod_CPLMmaster
