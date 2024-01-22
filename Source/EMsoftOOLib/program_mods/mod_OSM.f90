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

module mod_OSM
  !! author: MDG
  !! version: 1.0
  !! date: 04/06/20
  !!
  !! class definition for the EMgetOSM program

use mod_kinds
use mod_global

IMPLICIT NONE

! namelist for the EMgetOSM program
type, public :: OSMNameListType
  integer(kind=irg)  :: nmatch(5)
  character(fnlen)   :: dotproductfile
  character(fnlen)   :: tiffname
end type OSMNameListType

! class definition
type, public :: OSM_T
private
  character(fnlen)       :: nmldeffile = 'EMgetOSM.nml'
  type(OSMNameListType)  :: nml

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: OSM_
  procedure, pass(self) :: get_nmatch_
  procedure, pass(self) :: get_dotproductfile_
  procedure, pass(self) :: get_tiffname_
  procedure, pass(self) :: set_nmatch_
  procedure, pass(self) :: set_dotproductfile_
  procedure, pass(self) :: set_tiffname_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: OSM => OSM_
  generic, public :: get_nmatch => get_nmatch_
  generic, public :: get_dotproductfile => get_dotproductfile_
  generic, public :: get_tiffname => get_tiffname_
  generic, public :: set_nmatch => set_nmatch_
  generic, public :: set_dotproductfile => set_dotproductfile_
  generic, public :: set_tiffname => set_tiffname_

end type OSM_T

! the constructor routine for this class
interface OSM_T
  module procedure OSM_constructor
end interface OSM_T

contains

!--------------------------------------------------------------------------
type(OSM_T) function OSM_constructor( nmlfile ) result(OSM)
!DEC$ ATTRIBUTES DLLEXPORT :: OSM_constructor
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! constructor for the OSM_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call OSM%readNameList(nmlfile)

end function OSM_constructor

!--------------------------------------------------------------------------
subroutine OSM_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: OSM_destructor
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! destructor for the OSM_T Class

IMPLICIT NONE

type(OSM_T), INTENT(INOUT)  :: self

call reportDestructor('OSM_T')

end subroutine OSM_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! read the namelist from an nml file for the OSM_T Class

use mod_io
use mod_EMsoft

IMPLICIT NONE

class(OSM_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft
type(IO_T)                           :: Message
logical                              :: skipread = .FALSE.

integer(kind=irg)       :: nmatch(5)
character(fnlen)        :: dotproductfile
character(fnlen)        :: tiffname

! define the IO namelist to facilitate passing variables to the program.
namelist  / getOSM / nmatch, dotproductfile, tiffname

! set the input parameters to default values
nmatch = (/ 20, 0, 0, 0, 0 /)
dotproductfile = 'undefined'
tiffname = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=getOSM)
    close(UNIT=dataunit,STATUS='keep')

! check for required entries
    if (trim(dotproductfile).eq.'undefined') then
        call Message%printError('readNameList:',' dot product file name is undefined in '//nmlfile)
    end if

    if (trim(tiffname).eq.'undefined') then
        call Message%printError('readNameList:',' output tiff file name is undefined in '//nmlfile)
    end if
 end if

! if we get here, then all appears to be ok, and we need to fill in the enl fields
self%nml%nmatch = nmatch
self%nml%dotproductfile = dotproductfile
self%nml%tiffname = tiffname

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! pass the namelist for the OSM_T Class to the calling program

IMPLICIT NONE

class(OSM_T), INTENT(INOUT)     :: self
type(OSMNameListType)           :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
function get_nmatch_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nmatch_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! get nmatch(5) from the OSM_T class

IMPLICIT NONE

class(OSM_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out

out = self%nml%nmatch(5)

end function get_nmatch_

!--------------------------------------------------------------------------
subroutine set_nmatch_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nmatch_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! set nmatch(5) in the OSM_T class

IMPLICIT NONE

class(OSM_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp(5)

self%nml%nmatch = inp

end subroutine set_nmatch_

!--------------------------------------------------------------------------
function get_dotproductfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dotproductfile_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! get dotproductfile from the OSM_T class

IMPLICIT NONE

class(OSM_T), INTENT(INOUT)     :: self
character(fnlen)                :: out

out = self%nml%dotproductfile

end function get_dotproductfile_

!--------------------------------------------------------------------------
subroutine set_dotproductfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_dotproductfile_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! set dotproductfile in the OSM_T class

IMPLICIT NONE

class(OSM_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)    :: inp

self%nml%dotproductfile = inp

end subroutine set_dotproductfile_

!--------------------------------------------------------------------------
function get_tiffname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_tiffname_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! get tiffname from the OSM_T class

IMPLICIT NONE

class(OSM_T), INTENT(INOUT)     :: self
character(fnlen)                :: out

out = self%nml%tiffname

end function get_tiffname_

!--------------------------------------------------------------------------
subroutine set_tiffname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_tiffname_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! set tiffname in the OSM_T class

IMPLICIT NONE

class(OSM_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)    :: inp

self%nml%tiffname = inp

end subroutine set_tiffname_

!--------------------------------------------------------------------------
subroutine OSM_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: OSM_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! perform the computations

use mod_EMsoft
use mod_io
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_DIsupport
use mod_DIfiles
use ISO_C_BINDING
use mod_image
use, intrinsic :: iso_fortran_env

IMPLICIT NONE

class(OSM_T), INTENT(INOUT)             :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname

type(HDF_T)                             :: HDF
type(HDFnames_T)                        :: HDFnames
type(IO_T)                              :: Message
type(DIfile_T)                          :: DIFT
type(DictionaryIndexingNameListType)    :: dinl

real(kind=sgl),allocatable              :: OSMmap(:,:)
integer(kind=irg)                       :: dims(2), dimsOSM(2), hdferr, io_int(2), osmnum, i
character(fnlen)                        :: fname, TIFF_filename, dpfile, groupname, dataset, DIfile
character(2)                            :: fnum
real(kind=sgl)                          :: ma, mi

! declare variables for use in object oriented image module
integer                                 :: iostat
character(len=128)                      :: iomsg
logical                                 :: isInteger
type(image_t)                           :: im
integer(int8)                           :: i8 (3,4)
integer(int8), allocatable              :: TIFF_image(:,:)

associate(osmnl=>self%nml, DIDT=>DIFT%DIDT)

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()
HDFnames = HDFnames_T()

call HDFnames%set_NMLfiles(SC_NMLfiles)
call HDFnames%set_NMLfilename(SC_DictionaryIndexingNML)
call HDFnames%set_NMLparameters(SC_NMLparameters)
call HDFnames%set_NMLlist(SC_DictionaryIndexingNameListType)

DIfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(osmnl%dotproductfile)
call DIFT%readDotProductFile(EMsoft, HDF, HDFnames, DIfile, hdferr, &
                             getTopMatchIndices = .TRUE.)
dinl = DIFT%getNameList()

! first get the number of different OSM values in the list (non-zero entries)
osmnum = 0
do i=1,5
  if (osmnl%nmatch(i).ne.0) osmnum = osmnum+1
end do

! check to ake sure that the requested osmnl%nmatch values is <= the available number
dims = shape( DIDT%TopMatchIndices )
do i=1,osmnum
  if (osmnl%nmatch(i).gt.dinl%nnk) then
   io_int(1) = osmnl%nmatch(i)
   io_int(2) = dims(1)
   call Message%WriteValue(' Number of requested OSM levels = ',io_int, 2, "(I3,'; available number = ',I3)")
   call Message%printMessage('   --> Resetting requested number to maximum available')
   osmnl%nmatch(i) = dims(1)
  end if
end do

! allocate memory for OSM map
if (sum(dinl%ROI).ne.0) then
  allocate(OSMmap( dinl%ROI(3), dinl%ROI(4) ) )
else
  allocate(OSMmap( dinl%ipf_wd, dinl%ipf_ht ) )
end if
OSMmap = 0.0

! allocate memory for image
dimsOSM = shape(OSMmap)
allocate(TIFF_image( dimsOSM(1), dimsOSM(2) ))

do i=1,osmnum
! compute the Orientation Similarity Map
  if (sum(dinl%ROI).ne.0) then
    call getOrientationSimilarityMap( dims, DIDT%TopMatchIndices, osmnl%nmatch(i), dinl%ROI(3), dinl%ROI(4), OSMmap)
  else
    call getOrientationSimilarityMap( dims, DIDT%TopMatchIndices, osmnl%nmatch(i), dinl%ipf_wd, dinl%ipf_ht, OSMmap)
  end if

  ! output the ADP map as a tiff file
  write(fnum,"(I2.2)") osmnl%nmatch(i)

! we need to add this as a dataset to the dot product file so that it becomes available
! to other programs...

  hdferr =  HDF%openFile(DIfile)

! open the Scan 1/EBSD/Data group; dictionary indexing files only have one "scan" in them...
  groupname = 'Scan 1'
    hdferr = HDF%openGroup(groupname)
  groupname = SC_EBSD
    hdferr = HDF%openGroup(groupname)
  groupname = SC_Data
    hdferr = HDF%openGroup(groupname)

  dataset = 'OSM_'//fnum
  if (sum(dinl%ROI).ne.0) then
    hdferr = HDF%writeDatasetFloatArray(dataset, OSMmap, dinl%ROI(3), dinl%ROI(4))
  else
    hdferr = HDF%writeDatasetFloatArray(dataset, OSMmap, dinl%ipf_wd, dinl%ipf_ht)
  end if

! and close the HDF5 dot product file
  call HDF%popall()

  fname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(osmnl%tiffname)//fnum//'.tiff'
  TIFF_filename = trim(fname)

  ! fill the image with whatever data you have (between 0 and 255)
  ma = maxval(OSMmap)
  mi = minval(OSMmap)

  TIFF_image = int(255 * (OSMmap-mi)/(ma-mi))
  im = image_t(TIFF_image)
  if(im%empty()) call Message%printMessage("EMgetOSM","failed to convert array to image")

  ! create the file
  call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
  if(0.ne.iostat) then
    call Message%printMessage("failed to write image to file : "//iomsg)
  else
    call Message%printMessage('OSM map written to '//trim(TIFF_filename))
  end if

  call im%clear()
  OSMmap = 0.0
end do

call closeFortranHDFInterface()

end associate

end subroutine OSM_



end module mod_OSM
