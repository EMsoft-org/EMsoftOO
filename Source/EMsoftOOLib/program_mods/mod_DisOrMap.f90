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

module mod_DisOrMap
  !! author: MDG 
  !! version: 1.0 
  !! date: 04/22/22
  !!
  !! class definition for the EMDisOrMap program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMDisOrMap program
type, public :: DisOrMapNameListType
  integer(kind=irg)   :: px(10)
  integer(kind=irg)   :: py(10)
  character(fnlen)    :: dotproductfile
  character(fnlen)    :: DisOrMapfile
end type DisOrMapNameListType

! class definition
type, public :: DisOrMap_T
private 
  character(fnlen)            :: nmldeffile = 'EMDisOrMap.nml'
  type(DisOrMapNameListType)  :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: DisOrMap_
  procedure, pass(self) :: get_px_
  procedure, pass(self) :: get_py_
  procedure, pass(self) :: get_dotproductfile_
  procedure, pass(self) :: get_DisOrMapfile_
  procedure, pass(self) :: set_px_
  procedure, pass(self) :: set_py_
  procedure, pass(self) :: set_dotproductfile_
  procedure, pass(self) :: set_DisOrMapfile_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: DisOrMap => DisOrMap_
  generic, public :: get_px => get_px_
  generic, public :: get_py => get_py_
  generic, public :: get_dotproductfile => get_dotproductfile_
  generic, public :: get_DisOrMapfile => get_DisOrMapfile_
  generic, public :: set_px => set_px_
  generic, public :: set_py => set_py_
  generic, public :: set_dotproductfile => set_dotproductfile_
  generic, public :: set_DisOrMapfile => set_DisOrMapfile_

end type DisOrMap_T

! the constructor routine for this class 
interface DisOrMap_T
  module procedure DisOrMap_constructor
end interface DisOrMap_T

contains

!--------------------------------------------------------------------------
type(DisOrMap_T) function DisOrMap_constructor( nmlfile ) result(DisOrMap)
!! author: MDG 
!! version: 1.0 
!! date: 04/22/22
!!
!! constructor for the DisOrMap_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call DisOrMap%readNameList(nmlfile)

end function DisOrMap_constructor

!--------------------------------------------------------------------------
subroutine DisOrMap_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 04/22/22
!!
!! destructor for the DisOrMap_T Class
 
IMPLICIT NONE

type(DisOrMap_T), INTENT(INOUT)  :: self 

call reportDestructor('DisOrMap_T')

end subroutine DisOrMap_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 04/22/22
!!
!! read the namelist from an nml file for the DisOrMap_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(DisOrMap_T), INTENT(INOUT)     :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg)   :: px(10)
integer(kind=irg)   :: py(10)
character(fnlen)    :: dotproductfile
character(fnlen)    :: DisOrMapfile

namelist / DisOrMap / px, py, dotproductfile, DisOrMapfile

! set the input parameters to default values
px = (/ 0,0,0,0,0,0,0,0,0,0 /)
py = (/ 0,0,0,0,0,0,0,0,0,0 /)
dotproductfile = 'undefined'
DisOrMapfile = 'undefined'

! read the name list, depending on the class type
if (.not.skipread) then
! read the namelist file
  open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
  read(UNIT=dataunit,NML=DisOrMap)
  close(UNIT=dataunit,STATUS='keep')

! check for required entries
  if (trim(dotproductfile).eq.'undefined') then
    call Message%printError('EMgetDisOrMap:',' dotproductfile file name is undefined in '//nmlfile)
  end if

  ! make sure the energyfile variable has been properly defined
  if (trim(DisOrMapfile).eq.'undefined') then
    call Message%printError('EMgetDisOrMap:',' DisOrMapfile parameter is undefined in '//nmlfile)
  end if
end if

self%nml%px = px
self%nml%py = py
self%nml%dotproductfile = dotproductfile
self%nml%DisOrMapfile = DisOrMapfile

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 04/22/22
!!
!! pass the namelist for the DisOrMap_T Class to the calling program

IMPLICIT NONE 

class(DisOrMap_T), INTENT(INOUT)          :: self
type(DisOrMapNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
function get_px_(self, i) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_px_
!! author: MDG 
!! version: 1.0 
!! date: 04/22/22
!!
!! get px from the DisOrMap_T class

IMPLICIT NONE 

class(DisOrMap_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)        :: i
integer(kind=irg)                    :: out

out = self%nml%px(i)

end function get_px_

!--------------------------------------------------------------------------
subroutine set_px_(self,inp, i)
!DEC$ ATTRIBUTES DLLEXPORT :: set_px_
!! author: MDG 
!! version: 1.0 
!! date: 04/22/22
!!
!! set px in the DisOrMap_T class

IMPLICIT NONE 

class(DisOrMap_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)        :: inp
integer(kind=irg), INTENT(IN)        :: i

self%nml%px(i) = inp

end subroutine set_px_

!--------------------------------------------------------------------------
function get_py_(self, i) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_py_
!! author: MDG 
!! version: 1.0 
!! date: 04/22/22
!!
!! get py from the DisOrMap_T class

IMPLICIT NONE 

class(DisOrMap_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)        :: i
integer(kind=irg)                    :: out

out = self%nml%py(i)

end function get_py_

!--------------------------------------------------------------------------
subroutine set_py_(self, inp, i)
!DEC$ ATTRIBUTES DLLEXPORT :: set_py_
!! author: MDG 
!! version: 1.0 
!! date: 04/22/22
!!
!! set py in the DisOrMap_T class

IMPLICIT NONE 

class(DisOrMap_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)        :: inp
integer(kind=irg), INTENT(IN)        :: i

self%nml%py(i) = inp

end subroutine set_py_

!--------------------------------------------------------------------------
function get_dotproductfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dotproductfile_
!! author: MDG 
!! version: 1.0 
!! date: 04/22/22
!!
!! get dotproductfile from the DisOrMap_T class

IMPLICIT NONE 

class(DisOrMap_T), INTENT(INOUT)     :: self
character(fnlen)                     :: out

out = self%nml%dotproductfile

end function get_dotproductfile_

!--------------------------------------------------------------------------
subroutine set_dotproductfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_dotproductfile_
!! author: MDG 
!! version: 1.0 
!! date: 04/22/22
!!
!! set dotproductfile in the DisOrMap_T class

IMPLICIT NONE 

class(DisOrMap_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)         :: inp

self%nml%dotproductfile = inp

end subroutine set_dotproductfile_

!--------------------------------------------------------------------------
function get_DisOrMapfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_DisOrMapfile_
!! author: MDG 
!! version: 1.0 
!! date: 04/22/22
!!
!! get DisOrMapfile from the DisOrMap_T class

IMPLICIT NONE 

class(DisOrMap_T), INTENT(INOUT)     :: self
character(fnlen)                     :: out

out = self%nml%DisOrMapfile

end function get_DisOrMapfile_

!--------------------------------------------------------------------------
subroutine set_DisOrMapfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_DisOrMapfile_
!! author: MDG 
!! version: 1.0 
!! date: 04/22/22
!!
!! set DisOrMapfile in the DisOrMap_T class

IMPLICIT NONE 

class(DisOrMap_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)         :: inp

self%nml%DisOrMapfile = inp

end subroutine set_DisOrMapfile_

!--------------------------------------------------------------------------
subroutine DisOrMap_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: DisOrMap_
!! author: MDG 
!! version: 1.0 
!! date: 04/22/22
!!
!! perform the computations

use mod_EMsoft
use mod_HDFnames
use HDF5
use mod_HDFsupport
use mod_DIfiles
use mod_memory
use mod_so3
use mod_quaternions
use mod_rotations
use ISO_C_BINDING
use mod_image
use, intrinsic :: iso_fortran_env

IMPLICIT NONE 

class(DisOrMap_T), INTENT(INOUT)      :: self
type(EMsoft_T), INTENT(INOUT)         :: EMsoft
character(fnlen), INTENT(INOUT)       :: progname 

type(HDF_T)                           :: HDF
type(HDFnames_T)                      :: HDFnames
type(IO_T)                            :: Message
type(DIfile_T)                        :: DIFT
type(memory_T)                        :: mem
type(DictionaryIndexingNameListType)  :: dinl
type(so3_T)                           :: SO
type(QuaternionArray_T)               :: qAR, qdummy, qInp
type(Quaternion_T)                    :: quat
type(e_T)                             :: eu 
type(q_T)                             :: qu, quref
type(a_T)                             :: disax

character(fnlen)                      :: DIfile, TIFF_filename
integer(kind=irg)                     :: hdferr, ipf_x, ipf_y, ipf_wd, ipf_ht, i, j, nump, &
                                         io_int(2), numdis, Pmdims, FZt, FZo, ix, iy, mp(1), ipos
real(kind=dbl),allocatable            :: disor(:,:), qdinp(:,:), dmap(:,:)
integer(kind=irg),allocatable         :: indexmap(:,:), line(:)
real(kind=dbl)                        :: ax(4), mi, ma, io_real(2), qq(4)

! declare variables for use in object oriented image module
integer                               :: iostat
character(len=128)                    :: iomsg
logical                               :: isInteger
type(image_t)                         :: im
integer(int8)                         :: i8 (3,4), ival
integer(int8), allocatable            :: TIFF_image(:,:)

associate(domnl=>self%nml, DIDT=>DIFT%DIDT)

mem = memory_T() 

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()
HDFnames = HDFnames_T()

call HDFnames%set_NMLfiles(SC_NMLfiles)
call HDFnames%set_NMLfilename(SC_DictionaryIndexingNML)
call HDFnames%set_NMLparameters(SC_NMLparameters)
call HDFnames%set_NMLlist(SC_DictionaryIndexingNameListType)

! read data from the dot product file
DIfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(domnl%dotproductfile)
call DIFT%readDotProductFile(EMsoft, HDF, HDFnames, DIfile, hdferr, &
                             getRefinedEulerAngles = .TRUE.)
dinl = DIFT%getNameList()

! get the ROI if there is one
if (sum(dinl%ROI).ne.0) then ! the dp file has an ROI that is smaller than the complete field of view
  ipf_x = dinl%ROI(1)
  ipf_y = dinl%ROI(2)
  ipf_wd = dinl%ROI(3)
  ipf_ht = dinl%ROI(4)
else ! we use the full range
  ipf_x = 1
  ipf_y = 1
  ipf_wd = dinl%ipf_wd
  ipf_ht = dinl%ipf_ht
end if 
numdis = ipf_wd * ipf_ht
io_int(1:2) = (/ ipf_wd, ipf_ht /)
call Message%WriteValue('Size of the IPF map : ', io_int, 2, frm="(I5,' by ',I5)")

! check that all the requested (px,py) pairs are inside the ROI and count how many
! reference points there are
nump = 0
do i = 1,10
  if (maxval((/ domnl%px(i), domnl%py(i) /)).ne.0) then
    nump = nump + 1
    if (.not.((domnl%px(i).le.ipf_wd).and.(domnl%py(i).le.ipf_ht))) then
      io_int(1:2) = (/ domnl%px(i), domnl%py(i) /)
      call Message%WriteValue('selected point : ', io_int, 2)
      call Message%printError('EMgetDisOrMap','this input point is outside the valid range')
    end if
  end if
end do

io_int(1) = nump 
call Message%WriteValue(' Number of valid reference points found : ', io_int, 1)

! perform all rotation operations in double precision
call setRotationPrecision( 'd' )
call mem%alloc(qdinp, (/ 4, numdis /), 'qdinp', initval=0.D0)
do i=1,numdis
  eu = e_T( edinp = dble(DIDT%RefinedEulerAngles(1:3,i)) )
  qu = eu%eq()
  qdinp(1:4,i) = qu%q_copyd()
end do

! put the refined orientations in a quaternion array
qInp = QuaternionArray_T( n = numdis, qd = qdinp )

! allocate array for the disorientation values
call mem%alloc(disor, (/ numdis, nump /), 'disor', initval=0.D0)

! set up the correct SO3 parameters for fundamentalzone reduction
SO = so3_T( DIDT%pgnum )
call qdummy%QSym_Init( DIDT%pgnum, qAR )
Pmdims = qAR%getQnumber()
call SO%getFZtypeandOrder(FZt, FZo)

! here is the main loop over all reference points
do i=1,nump 
! get the reference point Euler angles and convert to quaternion quref
  ipos = (domnl%py(i)-1) * ipf_wd + domnl%px(i)
  eu = e_T( edinp = dble(DIDT%RefinedEulerAngles(1:3, ipos)) )
  quref = eu%eq()
! loop over all the input orientations to get the disorientation w.r.t. the 
! reference orientation quref
  do j=1,numdis
    quat = qInp%getQuatfromArray(j)
    qu = q_T( qdinp = quat%get_quatd() )
    call SO%getDisorientation(qAR, quref, qu, disax, fix1=.TRUE.)
    ax = disax%a_copyd()
    disor(j, i) = ax(4)/dtor
  end do 
end do 

! for each pixel, we determine for which reference point the disorientation 
! angle is the smallest one (if nump gt 1); we then create a map of the indices
! of the reference points
if (nump.eq.1) then 
  call mem%alloc(indexmap, (/ ipf_wd, ipf_ht /), 'indexmap', initval=1)
else
  call mem%alloc(indexmap, (/ ipf_wd, ipf_ht /), 'indexmap', initval=0)
  call mem%alloc(line, (/ nump /), 'line', initval=0)
  do iy=1, ipf_ht 
    do ix=1, ipf_wd 
      i = (iy-1) * ipf_wd + ix
      line = disor(i,1:nump)
      mp = minloc(line)
      indexmap(ix, iy) = mp(1)
    end do 
  end do 
  call mem%dealloc(line, 'line')
end if 

! next we generate the disorientation maps for each reference point
call mem%alloc(dmap, (/ ipf_wd, ipf_ht /), 'dmap', initval = 0.D0)
allocate(TIFF_image(ipf_wd, ipf_ht))

do iy=1, ipf_ht 
  do ix=1, ipf_wd 
    dmap(ix,iy) = disor( (iy-1) * ipf_wd + ix, indexmap(ix,iy) )
  end do 
end do 

TIFF_filename = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(domnl%DisOrMapfile)

! scale and output the disorientation image
mi = 0.D0
ma = maxval(dmap)
io_real(1:2) = (/ mi, ma /)
call Message%WriteValue(' Range of disorientations in montage : ',io_real, 2)

TIFF_image = int( 255 * (dmap-mi)/(ma-mi) ) 

! draw a small cross at each of the reference point positions
ival = 127_int8
do i=1,nump
  TIFF_image(maxval((/domnl%px(i)-3,1/)):minval((/domnl%px(i)+3,ipf_wd/)),domnl%py(i)) = ival
  TIFF_image(domnl%px(i),maxval((/domnl%py(i)-3,1/)):minval((/domnl%py(i)+3,ipf_ht/))) = ival
end do

im = image_t(TIFF_image)
if(im%empty()) call Message%printMessage("EMgetDisOrMap","failed to convert array to image")

! create the file
call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
    call Message%printMessage("failed to write image to file : "//iomsg)
else
    call Message%printMessage('Disorientation map written to '//trim(TIFF_filename))
end if

call im%clear()

deallocate(TIFF_image)
call mem%dealloc(dmap, 'dmap')
call mem%dealloc(disor, 'disor')
call mem%dealloc(qdinp, 'qdinp')

end associate

end subroutine DisOrMap_

end module mod_DisOrMap