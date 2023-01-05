! ###################################################################
! Copyright (c) 2018-2023, Marc De Graef Research Group/Carnegie Mellon University
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
module mod_vendors
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! a class to read and write various vendor file formats (used to be called "patternmod.f90")
!! with contributions from Saransh Singh, Michael Atkinson (Oxford Instruments binary pattern files)
!! and Håkon Ånes (NORDIF binary pattern file).
!!
!! This class also contains all the routines for writing .ang and .ctf files.

use mod_kinds
use mod_global
use mod_io
use HDF5
use mod_HDFsupport
use mod_math

IMPLICIT NONE

! the following arrays are only used for the Bruker HDF5 format, since the order of the
! EBSD patterns in the RawPatterns array is not necessarily the correct order !  We use these
! two index arrays to obtain the correct order ...
integer(kind=irg),allocatable,save    :: semix(:)
integer(kind=irg),allocatable,save    :: semiy(:)
integer(HSIZE_T),save                 :: semixydims(1)

! these are used to keep track of the even/odd patterns start locations in the .up1 and .up2 input formats
logical,save                          :: up1wdLeven, up1halfshift
logical,save                          :: up2wdLeven, up2halfshift
integer(kind=ill),save                :: offset

private :: get_num_HDFgroups_

type, public :: Vendor_T
  private
    character(fnlen)                  :: inputtype
    integer(kind=irg)                 :: itype
    character(fnlen)                  :: filename
    character(fnlen)                  :: Modality
    integer(kind=irg)                 :: funit = 55

contains
  private
    procedure, pass(self) :: determine_input_type_
    procedure, pass(self) :: invert_ordering_arrays_
    procedure, pass(self) :: openExpPatternFile_
    procedure, pass(self) :: getExpPatternRow_
    procedure, pass(self) :: getSingleExpPattern_
    procedure, pass(self) :: closeExpPatternFile_
    procedure, pass(self) :: get_Modality_
    procedure, pass(self) :: set_Modality_
    procedure, pass(self) :: get_inputtype_
    procedure, pass(self) :: get_itype_
    procedure, pass(self) :: get_filename_
    procedure, pass(self) :: get_funit_
    procedure, pass(self) :: set_inputtype_
    procedure, pass(self) :: set_itype_
    procedure, pass(self) :: set_filename_
    procedure, pass(self) :: set_funit_
    procedure, pass(self) :: ctf_writeFile_
    procedure, pass(self) :: ang_writeFile_
    procedure, pass(self) :: ctfmerge_writeFile_
    procedure, pass(self) :: angmerge_writeFile_
    final :: Vendor_destructor

    generic, public :: openExpPatternFile => openExpPatternFile_
    generic, public :: getExpPatternRow => getExpPatternRow_
    generic, public :: getSingleExpPattern => getSingleExpPattern_
    generic, public :: closeExpPatternFile => closeExpPatternFile_
    generic, public :: get_Modality => get_Modality_
    generic, public :: set_Modality => set_Modality_
    generic, public :: get_inputtype => get_inputtype_
    generic, public :: get_itype => get_itype_
    generic, public :: get_filename => get_filename_
    generic, public :: get_funit => get_funit_
    generic, public :: set_inputtype => set_inputtype_
    generic, public :: set_itype => set_itype_
    generic, public :: set_filename => set_filename_
    generic, public :: set_funit => set_funit_
    generic, public :: ctf_writeFile => ctf_writeFile_
    generic, public :: ang_writeFile => ang_writeFile_
    generic, public :: ctfmerge_writeFile => ctfmerge_writeFile_
    generic, public :: angmerge_writeFile => angmerge_writeFile_

end type Vendor_T

! the constructor routine for this class
interface Vendor_T
  module procedure Vendor_constructor
end interface Vendor_T

contains

!--------------------------------------------------------------------------
type(Vendor_T) function Vendor_constructor( inputtype ) result(VT)
!DEC$ ATTRIBUTES DLLEXPORT :: Vendor_constructor
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! constructor for the Vendor_T Class

IMPLICIT NONE

character(fnlen), INTENT(IN), OPTIONAL  :: inputtype

if (present(inputtype)) then
  call VT%set_inputtype( inputtype )
end if

end function Vendor_constructor

!--------------------------------------------------------------------------
subroutine Vendor_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: Vendor_destructor
!! author: MDG
!! version: 1.0
!! date: 02/02/20
!!
!! destructor for the MRC_T Class

IMPLICIT NONE

type(Vendor_T), INTENT(INOUT)   :: self

call reportDestructor('Vendor_T')

if (allocated(semix)) deallocate(semix)
if (allocated(semiy)) deallocate(semiy)

end subroutine Vendor_destructor

!--------------------------------------------------------------------------
!
! FUNCTION: determine_input_type_
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert inputtype string to integer value
!
!> @param inputtype
!
!> @date 02/13/18 MDG 1.0 original
!--------------------------------------------------------------------------
recursive function determine_input_type_(self) result(itype)
!DEC$ ATTRIBUTES DLLEXPORT :: determine_input_type_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! convert inputtype string to integer value

IMPLICIT NONE

class(Vendor_T), INTENT(INOUT)          :: self
integer(kind=irg)                       :: itype

type(IO_T)                              :: Message

self%itype = -1

if (trim(self%inputtype).eq."Binary") call self%set_itype(1)
if (trim(self%inputtype).eq."TSLup1") call self%set_itype(2)
if (trim(self%inputtype).eq."TSLup2") call self%set_itype(3)
if (trim(self%inputtype).eq."TSLHDF") call self%set_itype(4)
if (trim(self%inputtype).eq."OxfordBinary") call self%set_itype(5)
if (trim(self%inputtype).eq."OxfordHDF") call self%set_itype(6)
if (trim(self%inputtype).eq."EMEBSD") call self%set_itype(7)
if (trim(self%inputtype).eq."BrukerHDF") call self%set_itype(8)
if (trim(self%inputtype).eq."NORDIF") call self%set_itype(9)
if (trim(self%inputtype).eq."EMEBSD32i") call self%set_itype(10) ! for 32-bit integer pattern files
if (trim(self%inputtype).eq."EMEBSD32f") call self%set_itype(11) ! for 32-bit float pattern files

if (self%itype.eq.-1) call Message%printError('get_input_type','invalid file input type')
itype = self%get_itype()

end function determine_input_type_

!--------------------------------------------------------------------------
recursive function get_num_HDFgroups_(HDFstrings) result(numg)
!DEC$ ATTRIBUTES DLLEXPORT :: get_num_HDFgroups_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! extract the number of HDF groups from the HDFstrings array

IMPLICIT NONE

character(fnlen),INTENT(IN)             :: HDFstrings(10)
integer(kind=irg)                       :: numg

integer(kind=irg)                       :: i

numg = 0
do i=1,10
  if (len(trim(HDFstrings(i))).gt.0) numg = numg+1
end do
numg = numg-1   ! the last one should be a data set name

end function get_num_HDFgroups_

!--------------------------------------------------------------------------
recursive subroutine invert_ordering_arrays_(self, npat)
!DEC$ ATTRIBUTES DLLEXPORT :: invert_ordering_arrays_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! invert the pattern reordering arrays for Bruker HDF5 format

IMPLICIT NONE

class(Vendor_T), INTENT(INOUT)     :: self
integer(kind=irg),INTENT(IN)       :: npat

integer(kind=irg),allocatable      :: semixnew(:), semiynew(:)
integer(kind=irg)                  :: i, ix, iy, ipos

! allocate the new reordering arrays
allocate(semixnew(semixydims(1)), semiynew(semixydims(1)))

! invert the coordinate arrays  [tested on 2/16/18, MDG]
do i=1,semixydims(1)
  ix = mod(i, npat)-1
  iy = i/npat
  if (ix.lt.0) then
    ix = npat-1
    iy = iy-1
  end if
  ipos = semiy(i) * npat + semix(i) + 1
  semixnew(ipos) = ix
  semiynew(ipos) = iy
end do

! copy the new arrays over the old ones
semix = semixnew
semiy = semiynew
deallocate(semixnew, semiynew)

end subroutine invert_ordering_arrays_

!--------------------------------------------------------------------------
function get_Modality_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Modality_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! get Modality from the Vendor_T class

IMPLICIT NONE

class(Vendor_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = self%Modality

end function get_Modality_

!--------------------------------------------------------------------------
subroutine set_Modality_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Modality_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! set Modality in the Vendor_T class

IMPLICIT NONE

class(Vendor_T), INTENT(INOUT)     :: self
character(*), INTENT(IN)           :: inp

self%Modality = inp

end subroutine set_Modality_

!--------------------------------------------------------------------------
function get_inputtype_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_inputtype_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! get inputtype from the Vendor_T class

IMPLICIT NONE

class(Vendor_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = self%inputtype

end function get_inputtype_

!--------------------------------------------------------------------------
subroutine set_inputtype_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_inputtype_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! set inputtype in the Vendor_T class

IMPLICIT NONE

class(Vendor_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%inputtype = inp
self%itype = self%determine_input_type_()

end subroutine set_inputtype_

!--------------------------------------------------------------------------
function get_itype_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_itype_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! get itype from the Vendor_T class

IMPLICIT NONE

class(Vendor_T), INTENT(INOUT)     :: self
integer(kind=irg)                  :: out

out = self%itype

end function get_itype_

!--------------------------------------------------------------------------
subroutine set_itype_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_itype_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! set itype in the Vendor_T class

IMPLICIT NONE

class(Vendor_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)      :: inp

self%itype = inp

end subroutine set_itype_

!--------------------------------------------------------------------------
function get_filename_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_filename_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! get filename from the Vendor_T class

IMPLICIT NONE

class(Vendor_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = self%filename

end function get_filename_

!--------------------------------------------------------------------------
subroutine set_filename_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_filename_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! set filename in the Vendor_T class

IMPLICIT NONE

class(Vendor_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%filename = inp

end subroutine set_filename_

!--------------------------------------------------------------------------
function get_funit_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_funit_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! get funit from the Vendor_T class

IMPLICIT NONE

class(Vendor_T), INTENT(INOUT)     :: self
integer(kind=irg)                  :: out

out = self%funit

end function get_funit_

!--------------------------------------------------------------------------
subroutine set_funit_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_funit_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! set funit in the Vendor_T class (defaults to 55)

IMPLICIT NONE

class(Vendor_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)      :: inp

self%funit = inp

end subroutine set_funit_

!--------------------------------------------------------------------------
recursive function openExpPatternFile_(self, EMsoft, npat, L, recsize, HDFstrings, HDF) result(istat)
!DEC$ ATTRIBUTES DLLEXPORT :: openExpPatternFile_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! open a file with experimental patterns for a given input file type

use mod_EMsoft

IMPLICIT NONE

class(Vendor_T), INTENT(INOUT)        :: self
type(EMsoft_T), INTENT(INOUT)         :: EMsoft
integer(kind=irg),INTENT(IN)          :: npat
integer(kind=irg),INTENT(IN)          :: L
integer(kind=irg),INTENT(IN)          :: recsize
character(fnlen),INTENT(IN),OPTIONAL  :: HDFstrings(10)
type(HDF_T),INTENT(INOUT),OPTIONAL    :: HDF
integer(kind=irg)                     :: istat

type(IO_T)                            :: Message
character(fnlen)                      :: ename
integer(kind=irg)                     :: i, ierr, io_int(1), itype, hdferr, hdfnumg, recordsize, up2header(4), &
                                         ios, up1header(4), version, patx, paty, myoffset
character(fnlen)                      :: groupname, dataset, platform
logical                               :: f_exists

istat = 0

! for HDF mode, determine how many HDFgroups there are; the last entry in HDFstrings should be the data set name
if (present(HDFstrings)) then
  hdfnumg = get_num_HDFgroups_(HDFstrings)
end if

if (self%filename(1:1).eq.trim(EMsoft%getConfigParameter('EMsoftnativedelimiter'))) then
  ename = trim(self%filename)
else
  ename = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(self%filename)
end if

f_exists = .FALSE.
inquire(file=trim(ename), exist=f_exists)

if (.not.f_exists) then
   call Message%printMessage(' Input file '//trim(ename)//' does not exist in this location ... ')
   call Message%printMessage(' Please check the input parameters in the namelist file.')
   call Message%printMessage(' ')
   call Message%printError('openExpPatternFile','Unrecoverable error; file not found')
end if

call Message%printMessage(' Pattern input file '//trim(ename))
call Message%printMessage('   input file type '//trim(self%inputtype))

platform = trim(EMsoft%getConfigParameter('EMsoftplatform'))

! depending on the inputtype, we open the input file in the appropriate way
select case (self%itype)
    case(1)  ! "Binary"
        recordsize = L*4      ! all platforms use record length in units of bytes
        open(unit=self%funit,file=trim(ename),&
            status='old',form='unformatted',access='direct',recl=recordsize,iostat=ierr)
        if (ierr.ne.0) then
            io_int(1) = ierr
            call Message%WriteValue("File open error; error type ",io_int,1)
            call Message%printError("openExpPatternFile","Cannot continue program")
        end if

    case(2,3)  ! "TSLup1", TSLup2"
        ! open the file in STREAM access mode to allow for byte-level access
        open(unit=self%funit,file=trim(ename), status='old',access='stream',iostat=ios)
        if (ios.ne.0) then
            io_int(1) = ios
            call Message%WriteValue("File open error; error type ",io_int,1)
            call Message%printError("openExpPatternFile","Cannot continue program")
        end if
        ! the first four 4-byte entries form a header with a version number (unimportant), then
        ! the two dimensions of patterns, and finally an offset parameter indicating at which byte
        ! the first pattern starts.  We don't really need the other parameters, but we'll read them anyway.
        read(unit=self%funit, iostat=ios) version, patx, paty, myoffset
        offset = myoffset + 1_ill
        if (ios.ne.0) then
            io_int(1) = ios
            call Message%WriteValue("Read error in .up1/2 file header",io_int,1)
            call Message%printError("openExpPatternFile","Cannot continue program")
        end if

    case(5)  ! "OxfordBinary"
        ! open the file in STREAM access mode to allow for byte-level access
        open(unit=self%funit,file=trim(ename), status='old',access='stream',iostat=ios)
        if (ios.ne.0) then
            io_int(1) = ios
            call Message%WriteValue("File open error; error type ",io_int,1)
            call Message%printError("openExpPatternFile","Cannot continue program")
        end if

    case(6)  ! "OxfordHDF"
! Fall 2022, Oxford has an HDF5 file with extension .h5oina; some of the patterns in the 
! "Processed Patterns" dataset may have been compressed using the lzf compression scheme.
! Reading them does not require anything special when the most recent version of the 
! EMsoftOO_SDK is used.  It does require that the HDF5_PLUGIN_PATH environment variable is 
! set to the location of the plugin dynamical libraries; in a development environment, this
! would be done with the following shell command (csh): 
!    setenv HDF5_PLUGIN_PATH /full/path/to/EMsoftOO_SDK/hdf5-1.12.2-Release/lib/plugin
! We remind the user here that this parameter needs to be set in order for the uncompression
! dynamically loaded library to be found...  The actual reading of the data set is identical
! to that of the TSLHDF format, but some or all patterns may be in compressed form.
! [MDG, 11/09/22]

        call Message%printMessage((/ " ================================================================ ", &
                                     "The Oxford HDF5 format uses the lzf compression algorithm to      ", &
                                     "compress patterns; to properly read those patterns and uncompress ", &
                                     "them, the HDF5_PLUGIN_PATH environmental parameter must be set.   ", &
                                     "Please make sure that this parameter points to the correct folder;", &
                                     "for instance, EMsoftOO_SDK/hdf5-1.12.2-Release/lib/plugin or      ", &
                                     "similar; the folder should contain a liblzf.so or similar library.", &
                                     " ================================================================ " /))

        ! HDF = HDF_T()   this needs to be done in the calling program !
        ! open the file
        hdferr =  HDF%openFile(ename, readonly=.TRUE.)
        if (hdferr.ne.0) call HDF%error_check('openExpPatternFile:HDF%openFile', hdferr)
        ! open all the groups to the correct level of the data set
        do i=1,hdfnumg
            groupname = trim(HDFstrings(i))
            hdferr = HDF%openGroup(groupname)
            if (hdferr.ne.0) call HDF%error_check('openExpPatternFile:HDF%openGroup: groupname issue, check for typos.', hdferr)
        end do
        ! and here we leave this file open so that we can read data blocks using the hyperslab mechanism;

    case(4, 7, 10, 11)  ! "TSLHDF", "EMEBSD", "EMEBSD32i", "EMEBSD32f"
        ! HDF = HDF_T()   this needs to be done in the calling program !
        ! open the file
        hdferr =  HDF%openFile(ename, readonly=.TRUE.)
        if (hdferr.ne.0) call HDF%error_check('openExpPatternFile:HDF%openFile', hdferr)
        ! open all the groups to the correct level of the data set
        do i=1,hdfnumg
            groupname = trim(HDFstrings(i))
            hdferr = HDF%openGroup(groupname)
            if (hdferr.ne.0) call HDF%error_check('openExpPatternFile:HDF%openGroup: groupname issue, check for typos.', hdferr)
        end do
        ! and here we leave this file open so that we can read data blocks using the hyperslab mechanism;

    case(8)  !  "BrukerHDF"
        ! HDF = HDF_T()   this needs to be done in the calling program !
        ! open the file
        hdferr =  HDF%openFile(ename, readonly=.TRUE.)
        if (hdferr.ne.0) call HDF%error_check('openExpPatternFile:HDF%openFile', hdferr)

        ! open all the groups to the correct level of the data set
        do i=1,hdfnumg
            groupname = trim(HDFstrings(i))
            hdferr = HDF%openGroup(groupname)
             if (hdferr.ne.0) call HDF%error_check('openExpPatternFile:HDF%openGroup: groupname issue, check for typos.', hdferr)
            !  this part is different from the other vendors: the patterns are not necessarily in the correct order
            !  so we need to read the reordering arrays here...  The reordering arrays are always in the SEM group,
            !  which is one level down from the top (i.e., where we are right now).  Both arrays have the SAVE attribute.
            if (i.eq.1) then
               groupname = 'SEM'
               hdferr = HDF%openGroup(groupname)
               dataset = 'SEM IX'
               call HDF%readDatasetIntegerArray(dataset, semixydims, hdferr, semix)
               if (hdferr.ne.0) &
                 call HDF%error_check('openExpPatternFile:HDF%readDatasetIntegerArray: problem reading SEM IX array', hdferr)
               dataset = 'SEM IY'
               call HDF%readDatasetIntegerArray(dataset, semixydims, hdferr, semiy)
               if (hdferr.ne.0) &
                 call HDF%error_check('openExpPatternFile:HDF%readDatasetIntegerArray: problem reading SEM IY array', hdferr)
               call self%invert_ordering_arrays_(npat)
               call Message%printMessage('  found pattern reordering arrays')
               ! and leave this group
               call HDF%pop()
            end if
        end do
        ! and here we leave this file open so that we can read data blocks using the hyperslab mechanism;

    case(9)  !  "NORDIF"
        open(unit=self%funit, file=trim(ename), status='old', access='stream', iostat=ios)
        if (ios.ne.0) then
            io_int(1) = ios
            call Message%WriteValue("File open error; error type ", io_int, 1)
            call Message%printError("openExpPatternFile", "Cannot continue program")
        end if

    case default
        istat = -1
end select

end function openExpPatternFile_

!--------------------------------------------------------------------------
recursive subroutine getExpPatternRow_(self, iii, wd, patsz, L, dims3, offset3, exppatarray, ROI, flipy, HDFstrings, HDF)
!DEC$ ATTRIBUTES DLLEXPORT :: getExpPatternRow_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! read a single row of patterns from the input file(s)

IMPLICIT NONE

class(Vendor_T),INTENT(INOUT)           :: self
integer(kind=irg),INTENT(IN)            :: iii
integer(kind=irg),INTENT(IN)            :: wd
integer(kind=irg),INTENT(IN)            :: patsz
integer(kind=irg),INTENT(IN)            :: L
integer(HSIZE_T),INTENT(IN)             :: dims3(3)
integer(HSIZE_T),INTENT(IN)             :: offset3(3)
real(kind=sgl),INTENT(INOUT)            :: exppatarray(patsz * wd)
!f2py intent(in,out) ::  exppatarray
integer(kind=irg),OPTIONAL,INTENT(IN)   :: ROI(4)
logical,OPTIONAL,INTENT(IN)             :: flipy
character(fnlen),INTENT(IN),OPTIONAL    :: HDFstrings(10)
type(HDF_T),INTENT(INOUT),OPTIONAL      :: HDF

integer(kind=irg)                       :: hdfnumg, ierr, ios
real(kind=sgl)                          :: imageexpt(L), z
character(fnlen)                        :: dataset
character(kind=c_char),allocatable      :: EBSDpat(:,:,:)
integer(kind=irg),allocatable           :: EBSDpat32i(:,:,:)
real(kind=sgl),allocatable              :: EBSDpat32f(:,:,:)
integer(kind=C_INT16_T),allocatable     :: EBSDpatint(:,:,:)
character(1),allocatable                :: buffer(:)
integer(kind=ish),allocatable           :: pairs(:)
integer(kind=irg)                       :: sng, pixcnt
integer(kind=ish)                       :: pair(2)
integer(HSIZE_T)                        :: dims3new(3), offset3new(3), newspot
integer(kind=ill)                       :: recpos, ii, jj, kk, ispot, liii, lpatsz, lwd, lL, buffersize, kspot, jspot, &
                                           kkstart, kkend, multfactor
integer(kind=8)                         :: patoffsets(wd)
logical                                 :: flip

flip = .FALSE.
if (present(flipy)) then
  if (flipy.eqv..TRUE.) flip = .TRUE.
end if

if (present(HDFstrings)) then
  hdfnumg = get_num_HDFgroups_(HDFstrings)
  if (hdfnumg.gt.0) dataset = trim(HDFstrings(hdfnumg+1))
end if

! are we dealing with a smaller ROI or the full field view ?
! we need to use ill-type integers since the numbers can get pretty large...
liii = iii
lwd = wd
lpatsz = patsz
lL = L
if (self%itype.eq.2) multfactor = 1_ill
if (self%itype.eq.3) multfactor = 2_ill

if (present(ROI)) then
 kkstart = ROI(1)
 kkend = kkstart + ROI(3) - 1_ill
! for the TSL up1 and up2 formats we need to skip the first ROI(2)-1
! rows and set the correct offset value (in bytes)
 if (((self%itype.eq.2).or.(self%itype.eq.3)).and.(iii.eq.ROI(2))) then   ! make sure we do this only once ...
   do ii=1,ROI(2)-1
     offset = offset + (lwd * lL) * multfactor   ! this is in units of bytes
   end do
 end if
else
 kkstart = 1_ill
 kkend = dims3(3)
end if

select case (self%itype)
    case(1)  ! "Binary"
! This is the original EMsoft binary format that we used initially for indexing runs
! when the experimental patterns were only available in individual image file format.
! This file would have been created using a Matlab or IDL routine.  We anticipate that
! this format will not be used for much longer.
! In view of the pattern flip resolution, the user must ensure that the Matlab script
!  DOES NOT flip the pattern upside down !
      do jj=kkstart,kkend
        read(self%funit,rec=(liii-1)*lwd + jj) imageexpt
        exppatarray((jj-kkstart)*patsz+1:(jj-1)*patsz+L) = imageexpt(1:L)
      end do

    case(2,3)  ! "TSLup1", TSLup2"
! up1 file has single bytes as entries, up2 has 2-byte unsigned integers
      ! generate a buffer of the correct size ...
      buffersize = (lwd * lL) * multfactor
      allocate(buffer(buffersize))
! first we read the entire buffer as bytes
      read(unit=self%funit, pos=offset, iostat=ios) buffer

! then we convert the byte values into single byte or 2-byte integers
      if (multfactor.eq.2_ill) then ! .up2 format
        allocate(pairs(buffersize/2_ill))
        pairs = transfer(buffer,pairs)
      else ! .up1 format
        allocate(pairs(buffersize))
        do jj=1_ill,buffersize
         pairs(jj) = ichar(buffer(jj))
        end do
      end if
      deallocate(buffer)

! then we need to place them in the exppatarray array
      exppatarray = 0.0
      pixcnt = (kkstart-1)*dims3(1)*dims3(2)+1
      do kk=kkstart,kkend   ! loop over all the patterns in this row/ROI
        kspot = (kk-kkstart)*patsz
        do jj=1,dims3(2)
          if (flip.eqv..TRUE.) then
            jspot = (dims3(2)-jj)*dims3(1)
          else
            jspot = (jj-1)*dims3(1)
          end if
          do ii=1,dims3(1)
            exppatarray(kspot+jspot+ii) = float(pairs(pixcnt))
            pixcnt = pixcnt + 1
          end do
        end do
      end do

! increment the row offset parameter (in bytes)
      offset = offset + (lwd * lL) * multfactor
      deallocate(pairs)

! finally, correct for the fact that the original values were unsigned integers
      if (self%itype.eq.3) then
        where(exppatarray.lt.0.0) exppatarray = exppatarray + 65536.0
      else
        where(exppatarray.lt.0.0) exppatarray = exppatarray + 256.0
      end if

    case(5)  ! "OxfordBinary"

! read position of patterns in file for a single row from the header
      read(unit=self%funit, pos=(liii-1)*lwd*8+9, iostat=ios) patoffsets

! generate a buffer to load individual patterns into
      buffersize = lL
      allocate(buffer(buffersize))

! allocate pairs to store all patterns in a row
      buffersize = lwd * lL
      allocate(pairs(buffersize))

      do ii=1,lwd
! read each pattern into buffer with the 16 bytes of metadata skipped
        read(unit=self%funit, pos=patoffsets(ii)+17_8, iostat=ios) buffer

! loop over pixels and convert the byte values into single byte integers
        do jj=1_ill,lL
          pairs((ii - 1)*lL + jj) = ichar(buffer(jj))
        enddo
      end do

      deallocate(buffer)

 ! then we need to place them in the exppatarray array
      exppatarray = 0.0
      pixcnt = (kkstart-1)*dims3(1)*dims3(2)+1
      do kk=kkstart,kkend   ! loop over all the patterns in this row/ROI
        kspot = (kk-kkstart)*patsz
        do jj=1,dims3(2)
          jspot = (jj-1)*dims3(1)
          do ii=1,dims3(1)
            exppatarray(kspot+jspot+ii) = float(pairs(pixcnt))
            pixcnt = pixcnt + 1
          end do
        end do
      end do

      deallocate(pairs)

! finally, correct for the fact that the original values were unsigned integers
      where(exppatarray.lt.0.0) exppatarray = exppatarray + 256.0

    case(6)  ! "OxfordHDF"
! read a hyperslab section from the HDF5 input file; note that these patterns may be 
! in lzf compressed form; the HDF read routine will transparently take care of the 
! uncompressing, but only if the proper library can be found; make sure the HDF5_PLUGIN_PATH
! environmental parameter has been set correctly.  [MDG, 11/09/22]
        EBSDpatint = HDF%readHyperslabIntegerArray3D(dataset, offset3, dims3)
        exppatarray = 0.0
        do kk=kkstart,kkend
            do jj=1,dims3(2)
                do ii=1,dims3(1)
                   z = float(EBSDpatint(ii,jj,kk))
                   if (z.lt.0.0) z = z+2.0**16
                   exppatarray((kk-kkstart)*patsz+(jj-1)*dims3(1)+ii) = z
                end do
            end do
        end do

    case(4)  ! "TSLHDF" passed tests on 2/14/18 by MDG
! read a hyperslab section from the HDF5 input file
        EBSDpatint = HDF%readHyperslabIntegerArray3D(dataset, offset3, dims3)
        exppatarray = 0.0
        do kk=kkstart,kkend
            do jj=1,dims3(2)
                do ii=1,dims3(1)
                   z = float(EBSDpatint(ii,jj,kk))
                   if (z.lt.0.0) z = z+2.0**16
                   exppatarray((kk-kkstart)*patsz+(jj-1)*dims3(1)+ii) = z
                end do
            end do
        end do


    case(7)  ! "EMEBSD" passed tests on 2/14/18 by MDG
! read a hyperslab section from the HDF5 input file
        EBSDpat = HDF%readHyperslabCharArray3D(dataset, offset3, dims3)
        exppatarray = 0.0
        do kk=kkstart,kkend
            do jj=1,dims3(2)
                do ii=1,dims3(1)
                      exppatarray((kk-kkstart)*patsz+(jj-1)*dims3(1)+ii) = float(ichar(EBSDpat(ii,jj,kk)))
                end do
            end do
        end do

    case(10)  ! "EMEBSD32i"
! read a hyperslab section from the HDF5 input file
        EBSDpat32i = HDF%readHyperslabIntegerArray3D(dataset, offset3, dims3)
        exppatarray = 0.0
        do kk=kkstart,kkend
            do jj=1,dims3(2)
                do ii=1,dims3(1)
                      exppatarray((kk-kkstart)*patsz+(jj-1)*dims3(1)+ii) = float(EBSDpat32i(ii,jj,kk))
                end do
            end do
        end do

    case(11)  ! "EMEBSD32f" passed tests on 2/14/18 by MDG
! read a hyperslab section from the HDF5 input file
        EBSDpat32f = HDF%readHyperslabFloatArray3D(dataset, offset3, dims3)
        exppatarray = 0.0
        do kk=kkstart,kkend
            do jj=1,dims3(2)
                do ii=1,dims3(1)
                      exppatarray((kk-kkstart)*patsz+(jj-1)*dims3(1)+ii) = EBSDpat32f(ii,jj,kk)
                end do
            end do
        end do

    case(8)  ! "BrukerHDF"  passed tests on 2/16/18 by MDG
! since the pattern order in the Bruker data file is not necessarily the order in which the patterns
! were acquired, we need to read each patttern separately from the file using the appropriate offset, which
! is calculated using the semix and semiy arrays.  That means that we have to redefine both dims3 and offset3
! and loop over an entire row using the original pattern coordinate (ispot) as an index into the reordering arrays.
        exppatarray = 0.0
        dims3new = (/ dims3(1), dims3(2), 1_HSIZE_T /)
        do kk=kkstart,kkend  ! loop over all spots in the row/ROI
            ispot = (iii-1)*wd + kk
            newspot = semiy(ispot) * wd + semix(ispot)
            offset3new = (/ offset3(1), offset3(2),  newspot /)
            EBSDpat = HDF%readHyperslabCharArray3D(dataset, offset3new, dims3new)
            do jj=1,dims3(2)
                do ii=1,dims3(1)
                    exppatarray((kk-kkstart)*patsz+(jj-1)*dims3(1)+ii) = float(ichar(EBSDpat(ii,jj,1)))
                end do
            end do
        end do

    case(9)  !  "NORDIF"
        ! Buffer for single patterns
        buffersize = lL
        allocate(buffer(buffersize))

        ! Pairs to store all patterns in a row
        buffersize = lwd * lL
        allocate(pairs(buffersize))

        ! Loop over pixels and convert byte values into single byte integers
        do ii = 1, lwd
            ! pos = [(row-1)*scan_width + column - 1]*pattern_size + 1
            read(unit=self%funit, pos=((liii-1)*lwd + ii - 1)*lL + 1, iostat=ios) buffer
            do jj = 1_ill, lL
                pairs((ii-1)*lL + jj) = ichar(buffer(jj))
            end do
        end do
        deallocate(buffer)

        ! Place patterns in experimental pattern array
        exppatarray = 0.0
        ! Pattern pixels to read (might not be full pattern, depending on ROI)
        pixcnt = (kkstart-1)*dims3(1)*dims3(2) + 1
        ! Loop over row (might not be full row, depending on ROI)
        do kk = kkstart, kkend
            kspot = (kk-kkstart) * patsz
            ! Loop over rows of pattern pixels
            do jj = 1, dims3(2)
                jspot = (jj-1) * dims3(1)
                ! Loop over columns of pattern pixels, converting into float32
                do ii = 1, dims3(1)
                    exppatarray(kspot + jspot + ii) = float(pairs(pixcnt))
                    pixcnt = pixcnt + 1
                end do
            end do
        end do
        deallocate(pairs)

    case default
end select

end subroutine getExpPatternRow_

!--------------------------------------------------------------------------
recursive subroutine getSingleExpPattern_(self, iii, wd, patsz, L, dims3, offset3, exppat, HDFstrings, HDF)
!DEC$ ATTRIBUTES DLLEXPORT :: getSingleExpPattern_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! read a single experimental pattern from the input file

IMPLICIT NONE

class(Vendor_T),INTENT(INOUT)           :: self
integer(kind=irg),INTENT(IN)            :: iii
integer(kind=irg),INTENT(IN)            :: wd
integer(kind=irg),INTENT(IN)            :: patsz
integer(kind=irg),INTENT(IN)            :: L
integer(HSIZE_T),INTENT(IN)             :: dims3(3)
integer(HSIZE_T),INTENT(IN)             :: offset3(3)
real(kind=sgl),INTENT(INOUT)            :: exppat(patsz)
!f2py intent(in,out) ::  exppat
character(fnlen),INTENT(IN),OPTIONAL    :: HDFstrings(10)
type(HDF_T),INTENT(INOUT),OPTIONAL      :: HDF

integer(kind=irg)                       :: hdfnumg, ierr, ios, itype
real(kind=sgl)                          :: imageexpt(L), z
character(fnlen)                        :: dataset
character(kind=c_char),allocatable      :: EBSDpat(:,:,:)
integer(kind=irg),allocatable           :: EBSDpat32i(:,:,:)
real(kind=sgl),allocatable              :: EBSDpat32f(:,:,:)
integer(kind=C_INT16_T),allocatable     :: EBSDpatint(:,:,:)
character(1),allocatable                :: buffer(:)
integer(kind=ish),allocatable           :: pairs(:)
integer(kind=irg)                       :: sng, pixcnt
integer(kind=ish)                       :: pair(2)
integer(HSIZE_T)                        :: dims3new(3), offset3new(3), newspot
integer(kind=ill)                       :: recpos, ii, jj, kk, ispot, liii, lpatsz, lwd, lL, buffersize, kspot, jspot, &
                                           l1, l2, multfactor
integer(kind=8)                         :: patoffsets(wd)

if (present(HDFstrings)) then
  hdfnumg = get_num_HDFgroups_(HDFstrings)
  if (hdfnumg.gt.0) dataset = trim(HDFstrings(hdfnumg+1))
end if

if (itype.eq.2) multfactor = 1_ill
if (itype.eq.3) multfactor = 2_ill

select case (self%itype)
    case(1)  ! "Binary"
! This is the original EMsoft binary format that we used initially for indexing runs
! when the experimental patterns were only available in individual image file format.
! This file would have been created using a Matlab or IDL routine.  We anticipate that
! this format will not be used for much longer.  To call the routine for a single pattern,
! simply place y*wd+x in the third entry of the offset3 array.
        read(self%funit,rec=offset3(3)) imageexpt
        exppat(1:L) = imageexpt(1:L)


    case(2,3)  ! "TSL_up1", TSLup2"
! tested and compared to IDL version of read_up2 routine on 2/20/18, MDG.
! the requested pattern should be in the encoded in the third entry of the offset3 array.
! we need to use ill-type integers since the numbers can get pretty large...
        lwd = wd
        lpatsz = patsz
        lL = L
        l1 = offset3(3)
        if (itype.eq.2) then
         multfactor = 1_ill
        else
         multfactor = 2_ill
        end if

        offset = offset + (l1-1_ill) * lpatsz * multfactor
        buffersize = lpatsz * multfactor
        allocate(buffer(buffersize))
        read(unit=self%funit, pos=offset, iostat=ios) buffer

! then we convert the byte values into single byte or 2-byte integers
        if (multfactor.eq.2_ill) then ! .up2 format
          allocate(pairs(buffersize/2_ill))
          pairs = transfer(buffer,pairs)
        else ! .up1 format
          allocate(pairs(buffersize))
          do jj=1_ill,buffersize
           pairs(jj) = ichar(buffer(jj))
          end do
        end if
        deallocate(buffer)

! ! then we need to place them in the exppatarray array with the proper offsets if patsz ne L
        exppat = 0.0
        pixcnt = 1
        do jj=1,dims3(2)
          jspot = (jj-1)*dims3(1)
          do ii=1,dims3(1)
            exppat(jspot+ii) = float(pairs(pixcnt))
            pixcnt = pixcnt + 1
          end do
        end do
        deallocate(pairs)

! finally, correct for the fact that the original values were unsigned integers
        if (itype.eq.3) then
          where(exppat.lt.0.0) exppat = exppat + 65536.0
        else
          where(exppat.lt.0.0) exppat = exppat + 256.0
        end if

    case(5)  ! "OxfordBinary" ! [added/tested MDG 07/13/19]
! read position of patterns in file for a single row from the header
      liii = iii
      l1 = mod(offset3(3),wd)
      lL = L
      lwd = wd
      read(unit=self%funit, pos=(liii-1)*lwd*8+9, iostat=ios) patoffsets

! generate buffers to load individual pattern into
      buffersize = lL
      allocate(buffer(buffersize), pairs(buffersize))

! read single pattern into buffer with the 16 bytes of metadata skipped
      read(unit=self%funit, pos=patoffsets(l1)+17_8, iostat=ios) buffer

! convert the byte values into single byte integers
      pairs = ichar(buffer)
      deallocate(buffer)

 ! then we need to place it in the exppat array
      exppat = 0.0
      pixcnt = 1
      do jj=1,dims3(2)
        jspot = (jj-1)*dims3(1)
        do ii=1,dims3(1)
          exppat(jspot+ii) = float(pairs(pixcnt))
          pixcnt = pixcnt + 1
        end do
      end do
      deallocate(pairs)

! finally, correct for the fact that the original values were unsigned integers
      where(exppat.lt.0.0) exppat = exppat + 256.0

    case(6)  ! "OxfordHDF"
! read a hyperslab section from the HDF5 input file; note that these patterns may be 
! in lzf compressed form; the HDF read routine will transparently take care of the 
! uncompressing, but only if the proper library can be found; make sure the HDF5_PLUGIN_PATH
! environmental parameter has been set correctly.  [MDG, 11/09/22]

! read a hyperslab single pattern section from the HDF5 input file
! dims3 should have the pattern dimensions and then 1_HSIZE_T for the third dimension
! offset3 should have (0,0) and then the offset of the pattern (0-based)
        EBSDpatint = HDF%readHyperslabIntegerArray3D(dataset, offset3, dims3)
        exppat = 0.0
        do jj=1,dims3(2)
            do ii=1,dims3(1)
                  z = float(EBSDpatint(ii,jj,1))
                  if (z.lt.0.0) z = z+2.0**16
                  exppat((jj-1)*dims3(1)+ii) = z
            end do
        end do

    case(4)  ! "TSLHDF" passed tests on 2/20/18 by MDG
! read a hyperslab single pattern section from the HDF5 input file
! dims3 should have the pattern dimensions and then 1_HSIZE_T for the third dimension
! offset3 should have (0,0) and then the offset of the pattern (0-based)
        EBSDpatint = HDF%readHyperslabIntegerArray3D(dataset, offset3, dims3)
        exppat = 0.0
        do jj=1,dims3(2)
            do ii=1,dims3(1)
                  z = float(EBSDpatint(ii,jj,1))
                  if (z.lt.0.0) z = z+2.0**16
                  exppat((jj-1)*dims3(1)+ii) = z
            end do
        end do

    case(7)  ! "EMEBSD" passed tests on 2/20/18 by MDG
! read a hyperslab single pattern section from the HDF5 input file
! dims3 should have the pattern dimensions and then 1_HSIZE_T for the third dimension
! offset3 should have (0,0) and then the offset of the pattern (0-based)
        EBSDpat = HDF%readHyperslabCharArray3D(dataset, offset3, dims3)
        exppat = 0.0
        do jj=1,dims3(2)
            do ii=1,dims3(1)
                  exppat((jj-1)*dims3(1)+ii) = float(ichar(EBSDpat(ii,jj,1)))
            end do
        end do


    case(10)  ! "EMEBSD32i"
! read a hyperslab single pattern section from the HDF5 input file
! dims3 should have the pattern dimensions and then 1_HSIZE_T for the third dimension
! offset3 should have (0,0) and then the offset of the pattern (0-based)
        EBSDpat32i = HDF%readHyperslabIntegerArray3D(dataset, offset3, dims3)
        exppat = 0.0
        do jj=1,dims3(2)
            do ii=1,dims3(1)
                  exppat((jj-1)*dims3(1)+ii) = float(EBSDpat32i(ii,jj,1))
            end do
        end do

    case(11)  ! "EMEBSD32f"
! read a hyperslab single pattern section from the HDF5 input file
! dims3 should have the pattern dimensions and then 1_HSIZE_T for the third dimension
! offset3 should have (0,0) and then the offset of the pattern (0-based)
        EBSDpat32f = HDF%readHyperslabFloatArray3D(dataset, offset3, dims3)
        exppat = 0.0
        do jj=1,dims3(2)
            do ii=1,dims3(1)
                  exppat((jj-1)*dims3(1)+ii) = EBSDpat32f(ii,jj,1)
            end do
        end do

    case(8)  ! "BrukerHDF"  to be tested
! since the pattern order in the Bruker data file is not necessarily the order in which the patterns
! were acquired, we need to read each patttern separately from the file using the appropriate offset, which
! is calculated using the semix and semiy arrays.  That means that we have to redefine both dims3 and offset3
! and use the original pattern coordinate (ispot) as an index into the reordering arrays.
        exppat = 0.0
        dims3new = (/ dims3(1), dims3(2), 1_HSIZE_T /)
        ispot = offset3(3)
        newspot = semiy(ispot) * wd + semix(ispot)
        offset3new = (/ offset3(1), offset3(2),  newspot /)
        EBSDpat = HDF%readHyperslabCharArray3D(dataset, offset3new, dims3new)
        do jj=1,dims3(2)
            do ii=1,dims3(1)
                exppat((jj-1)*dims3(1)+ii) = float(ichar(EBSDpat(ii,jj,1)))
            end do
        end do

    case(9)  !  "NORDIF"
        ! Use ill-type integers
        lL = L
        lwd = wd

        ! Buffers for single patterns
        buffersize = lL
        allocate(buffer(buffersize), pairs(buffersize))

        ! Read single pattern into buffer
        ! offset3(3) = row * scan width + column
        offset = offset3(3)*lL + 1
        read(unit=self%funit, pos=offset, iostat=ios) buffer

        ! Convert byte values into single byte integers
        pairs = ichar(buffer)
        deallocate(buffer)

        ! Place pattern in experimental pattern array
        exppat = 0.0
        pixcnt = 1
        ! Loop over rows of pattern pixels
        do jj = 1, dims3(2)
            jspot = (jj-1)*dims3(1)
            ! Loop over columns of pattern pixels, converting into float32
            do ii = 1, dims3(1)
                exppat(jspot + ii) = float(pairs(pixcnt))
                pixcnt = pixcnt + 1
            end do
        end do
        deallocate(pairs)

    case default
end select

end subroutine getSingleExpPattern_

!--------------------------------------------------------------------------
recursive subroutine closeExpPatternFile_(self, HDF)
!DEC$ ATTRIBUTES DLLEXPORT :: closeExpPatternFile_
!! author: MDG
!! version: 1.0
!! date: 03/31/20
!!
!! close a file with experimental patterns for a given input file type

use mod_HDFsupport
use HDF5

IMPLICIT NONE

class(Vendor_T),INTENT(INOUT)           :: self
type(HDF_T),INTENT(INOUT),OPTIONAL      :: HDF

type(IO_T)                              :: Message
integer(kind=irg)                       :: hdferr

select case (self%itype)
    case(1, 2, 3, 5, 9)  ! "Binary" "TSLup1" "TSLup2" "OxfordBinary" "NORDIF"
        close(unit=self%funit,status='keep')

    case(4, 6, 7, 10, 11)  ! "TSLHDF" "EMEBSD"
        if (present(HDF)) call HDF%pop(.TRUE.)

    case(8)  !  "BrukerHDF"
        if (present(HDF)) call HDF%pop(.TRUE.)
        deallocate(semix, semiy)

    case default
        call Message%printError("closeExpPatternFile","unknown input format")
end select

end subroutine closeExpPatternFile_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! the following 4 routines replace the old EBSDiomod module
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine ctf_writeFile_(self,EMsoft,cell,SG,nml,ipar,fpar,indexmain,eulerarray,resultmain,OSMmap,IQmap,noindex, & 
                                    orthoset,orthoSG)
!DEC$ ATTRIBUTES DLLEXPORT :: ctf_writeFile_
!! author: MDG
!! version: 1.0
!! date: 04/01/20
!!
!! Write a *.ctf output file with EBSD/TKD data (HKL format)

use mod_EMsoft
use mod_crystallography
use mod_symmetry
use mod_io
use mod_DIfiles

IMPLICIT NONE

class(Vendor_T),INTENT(INOUT)                       :: self
type(EMsoft_T),INTENT(INOUT)                        :: EMsoft
type(cell_T),INTENT(INOUT)                          :: cell
type(SpaceGroup_T),INTENT(INOUT)                    :: SG
class(DictionaryIndexingNameListType),INTENT(INOUT) :: nml
!f2py intent(in,out) ::  nml
integer(kind=irg),INTENT(IN)                        :: ipar(10)
real(kind=sgl),INTENT(IN)                           :: fpar(2)
integer(kind=irg),INTENT(IN)                        :: indexmain(ipar(1),ipar(2))
real(kind=sgl),INTENT(IN)                           :: eulerarray(3,ipar(4))
real(kind=sgl),INTENT(IN)                           :: resultmain(ipar(1),ipar(2))
real(kind=sgl),INTENT(IN)                           :: OSMmap(ipar(7),ipar(8))
real(kind=sgl),INTENT(IN)                           :: IQmap(ipar(3))
logical,INTENT(IN),OPTIONAL                         :: noindex
integer(kind=irg),INTENT(IN),OPTIONAL               :: orthoset
integer(kind=irg),INTENT(IN),OPTIONAL               :: orthoSG

type(IO_T)                                          :: Message
integer(kind=irg)                                   :: ierr, i, ii, indx, hdferr, SGnum, LaueGroup, BCval, BSval
character(fnlen)                                    :: ctfname, xtalname, modality
character                                           :: TAB = CHAR(9)
character(fnlen)                                    :: str1,str2,str3,str4,str5,str6,str7,str8,str9,str10
real(kind=sgl)                                      :: euler(3), eu, mi, ma
logical                                             :: donotuseindexarray, isEBSD=.FALSE., isTKD=.FALSE., isECP=.FALSE.
real(kind=dbl)                                      :: cellparams(6), a, b, c
integer(kind=irg),allocatable                       :: osm(:), iq(:)

donotuseindexarray = .FALSE.
if (present(noindex)) then
  if (noindex.eqv..TRUE.) then
    donotuseindexarray = .TRUE.
  end if
end if

modality = trim(self%get_Modality())

select case(modality)
  case('EBSD')
    isEBSD = .TRUE.
  case('TKD')
    isTKD = .TRUE.
  case('ECP')
    isECP = .TRUE.
  case default
    call Message%printError('ctf_writeFile_', 'unknown name list type requested')
end select


! get the OSMmap into 1D format and scale to the range [0..255]
allocate(osm(ipar(3)))
mi = minval(OSMmap)
ma = maxval(OSMmap)
if (mi.eq.ma) then
  osm = 0
else
  indx = 1
  do i=1,ipar(8)
    do ii=1,ipar(7)
      osm(indx) = nint(255.0 * (OSMmap(ii,i)-mi)/(ma-mi))
      indx = indx+1
    end do
  end do
end if

! scale the IQmap to the range [0..255]
allocate(iq(ipar(3)))
iq = nint(255.0 * IQmap)

! open the file (overwrite old one if it exists)
ctfname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(nml%ctffile)
open(unit=dataunit2,file=trim(ctfname),status='unknown',action='write',iostat=ierr)

write(dataunit2,'(A)') 'Channel Text File'
str1 = 'EMsoft v. '//trim(EMsoft%getConfigParameter('EMsoftversion'))// &
       '; BANDS=pattern index, MAD=CI, BC=OSM, BS=IQ; Diffraction Modality: '//trim(modality) 
if (present(orthoset)) then
       str1 = trim(str1)//'; orthorhombic space group setting '//extendedHMOrthsymbols(orthoset,orthoSG)
end if 
write(dataunit2,'(A)') trim(str1)
write(dataunit2,'(A)') 'Author  '//trim(EMsoft%getConfigParameter('UserName'))
write(dataunit2,'(A)') 'JobMode Grid'
write(dataunit2,'(2A,I5)') 'XCells',TAB, ipar(7)
write(dataunit2,'(2A,I5)') 'YCells',TAB, ipar(8)
write(dataunit2,'(2A,F6.2)') 'XStep',TAB, nml%StepX
write(dataunit2,'(2A,F6.2)') 'YStep',TAB, nml%StepY
write(dataunit2,'(A)') 'AcqE1'//TAB//'0'
write(dataunit2,'(A)') 'AcqE2'//TAB//'0'
write(dataunit2,'(A)') 'AcqE3'//TAB//'0'
write(dataunit2,'(A,A,$)') 'Euler angles refer to Sample Coordinate system (CS0)!',TAB
str1 = 'Mag'//TAB//'30'//TAB//'Coverage'//TAB//'100'//TAB//'Device'//TAB//'0'//TAB//'KV'
write(str2,'(F4.1)') fpar(1)  ! EkeV
str1 = trim(str1)//TAB//trim(str2)//TAB//'TiltAngle'
write(str2,'(F5.2)') fpar(2)  ! MCsig
str2 = adjustl(str2)
str1 = trim(str1)//TAB//trim(str2)//TAB//'TiltAxis'//TAB//'0'
write(dataunit2,'(A)') trim(str1)
write(dataunit2,'(A)') 'Phases'//TAB//'1'

cellparams = cell%getLatParm()
SGnum = SG%getSpaceGroupNumber()
xtalname = trim(cell%getFileName())

! unit cell size
cellparams(1:3) = cellparams(1:3)*10.0  ! convert to Angstrom
if (present(orthoset)) then   ! permute the lattice parameters to the correct orthorhombic setting 
! settings  " a  b  c", " b  a -c", " c  a  b", "-c  b  a", " b  c  a", " a -c  b" 
  select case(orthoset)
    case(1) 
      a = cellparams(1)
      b = cellparams(2)
      c = cellparams(3)
    case(2) 
      a = cellparams(2)
      b = cellparams(1)
      c = cellparams(3)
    case(3) 
      a = cellparams(3)
      b = cellparams(1)
      c = cellparams(2)
    case(4) 
      a = cellparams(3)
      b = cellparams(2)
      c = cellparams(1)
    case(5) 
      a = cellparams(2)
      b = cellparams(3)
      c = cellparams(1)
    case(6) 
      a = cellparams(1)
      b = cellparams(3)
      c = cellparams(2)
  end select 
  write(str1,'(F8.3)') a
  write(str2,'(F8.3)') b
  write(str3,'(F8.3)') c
else
  write(str1,'(F8.3)') cellparams(1)
  write(str2,'(F8.3)') cellparams(2)
  write(str3,'(F8.3)') cellparams(3)
end if 
str1 = adjustl(str1)
str2 = adjustl(str2)
str3 = adjustl(str3)
str1 = trim(str1)//';'//trim(str2)//';'//trim(str3)

! unit cell angles
write(str4,'(F8.3)') cellparams(4)
write(str5,'(F8.3)') cellparams(5)
write(str6,'(F8.3)') cellparams(6)
str4 = adjustl(str4)
str5 = adjustl(str5)
str6 = adjustl(str6)
str1 = trim(str1)//TAB//trim(str4)//';'//trim(str5)//';'//trim(str6)

! structure name
str3 = ''
ii = len(trim(xtalname))-5
do i=1,ii
  str3(i:i) = xtalname(i:i)
end do
str1 = trim(str1)//TAB//trim(str3)

! rotational symmetry group
str4 = ''
LaueGroup = SG%getLaueGroupNumber()
write(str4,'(I2)') LaueGroup
str1 = trim(str1)//TAB//trim(adjustl(str4))

! space group
str2 = ''
write(str2,'(I3)') SGnum
str1 = trim(str1)//TAB//trim(adjustl(str2))

! and now collect them all into a single string
write(dataunit2,'(A)') trim(str1)

! this is the table header
write(dataunit2,'(A)') 'Phase'//TAB//'X'//TAB//'Y'//TAB//'Bands'//TAB//'Error'//TAB//'Euler1'//TAB//'Euler2'//TAB//'Euler3' &
                      //TAB//'MAD'//TAB//'BC'//TAB//'BS'

! go through the entire array and write one line per sampling point
do ii = 1,ipar(3)
    BCval = osm(ii)
    BSval = iq(ii)
    if (donotuseindexarray.eqv..TRUE.) then
      indx = 0
      euler = eulerarray(1:3,ii)
    else
      indx = indexmain(1,ii)
      euler = eulerarray(1:3,indx)
    end if
! changed order of coordinates to conform with ctf standard
    if (sum(nml%ROI).ne.0) then
      write(str2,'(F12.3)') float(floor(float(ii-1)/float(nml%ROI(3))))*nml%StepY
      write(str1,'(F12.3)') float(MODULO(ii-1,nml%ROI(3)))*nml%StepX
    else
      write(str2,'(F12.3)') float(floor(float(ii-1)/float(nml%ipf_wd)))*nml%StepY
      write(str1,'(F12.3)') float(MODULO(ii-1,nml%ipf_wd))*nml%StepX
    end if

    write(str3,'(I8)') indx  ! pattern index into dictionary list of discrete orientations
    write(str8,'(I8)') 0 ! integer zero error; was indx, which is now moved to BANDS
    eu = euler(1) - 90.0 ! conversion from TSL to Oxford convention
    if (eu.lt.0) eu = eu + 360.0
    write(str5,'(F12.3)') eu
    eu = euler(2)
    if (eu.lt.0) eu = eu + 360.0
    write(str6,'(F12.3)') eu
! intercept the hexagonal case, for which we need to subtract 30° from the third Euler angle
! Note: after working with Lionel Germain, we concluded that we do not need to subtract 30°
! in the ctf file, because the fundamental zone is already oriented according to the Oxford
! convention... That means that we need to subtract the angle for the .ang file (to be implemented)
! [modified by MDG on 3/5/18]
    if ((LaueGroup.eq.8).or.(LaueGroup.eq.9)) euler(3) = euler(3) - 30.0
    eu = euler(3)
    if (eu.lt.0) eu = eu + 360.0
    write(str7,'(F12.3)') eu
    write(str4,'(F12.6)') resultmain(1,ii)   ! this replaces MAD
! the following two parameters need to be modified to contain more meaningful information
    write(str9,'(I8)') BCval   ! OSM value in range [0..255]
    write(str10,'(I8)') BSval  !  IQ value in range [0..255]
! some Oxford 3D files have four additional integer columns;
! GrainIndex
! GrainRandomColourR
! GrainRandomColourG
! GrainRandomColourB
!
    write(dataunit2,'(A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A)')'1',TAB,trim(adjustl(str1)),TAB,&
    trim(adjustl(str2)),TAB,trim(adjustl(str3)),TAB,trim(adjustl(str8)),TAB,trim(adjustl(str5)),&
    TAB,trim(adjustl(str6)),TAB,trim(adjustl(str7)),TAB,trim(adjustl(str4)),TAB,trim(adjustl(str9)),&
    TAB,trim(adjustl(str10))
end do

close(dataunit2,status='keep')

end subroutine ctf_writeFile_

!--------------------------------------------------------------------------
!
! SUBROUTINE:ang_writeFile_
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief
!
!> @param nml namelist
!> @param ipar  series of integer dimensions
!> @param indexmain indices into the main Euler array
!> @param eulerarray array of Euler angle triplets
!> @param resultmain dot product array
!
!> @date 02/07/15  SS 1.0 original
!> @date 03/10/16 MDG 1.1 moved from program to module and updated [TO BE COMPLETED]
!> @date 11/08/18 MDG 2.0 rewrite and testing
!--------------------------------------------------------------------------
recursive subroutine ang_writeFile_(self,EMsoft,cell,SG,nml,ipar,fpar,indexmain,eulerarray,resultmain,IQmap,noindex,&
                                    orthoset,orthoSG)
!DEC$ ATTRIBUTES DLLEXPORT :: ang_writeFile_
!! author: MDG
!! version: 1.0
!! date: 04/01/20
!!
!! Write a *.ang output file with EBSD data (TSL format)

use mod_EMsoft
use mod_crystallography
use mod_symmetry
use mod_io
use mod_DIfiles

IMPLICIT NONE

class(Vendor_T),INTENT(INOUT)                       :: self
type(EMsoft_T),INTENT(INOUT)                        :: EMsoft
type(cell_T),INTENT(INOUT)                          :: cell
type(SpaceGroup_T),INTENT(INOUT)                    :: SG
class(DictionaryIndexingNameListType),INTENT(INOUT) :: nml
integer(kind=irg),INTENT(IN)                        :: ipar(10)
real(kind=sgl),INTENT(IN)                           :: fpar(1)
integer(kind=irg),INTENT(IN)                        :: indexmain(ipar(1),ipar(2))
real(kind=sgl),INTENT(IN)                           :: eulerarray(3,ipar(4))
real(kind=sgl),INTENT(IN)                           :: resultmain(ipar(1),ipar(2))
real(kind=sgl),INTENT(IN)                           :: IQmap(ipar(3))
logical,INTENT(IN),OPTIONAL                         :: noindex
integer(kind=irg),INTENT(IN),OPTIONAL               :: orthoset 
integer(kind=irg),INTENT(IN),OPTIONAL               :: orthoSG

type(IO_T)                                          :: Message
integer(kind=irg)                                   :: ierr, ii, indx, SGnum
character(fnlen)                                    :: angname, xtalname, modality
character(fnlen)                                    :: str1,str2,str3,str4,str5,str6,str7,str8,str9,str10
character                                           :: TAB = CHAR(9)
character(2)                                        :: TSLsymmetry
real(kind=sgl)                                      :: euler(3), s, BSval
real(kind=dbl)                                      :: cellparams(6), a, b, c
logical                                             :: donotuseindexarray, isEBSD=.FALSE., isTKD=.FALSE., isECP=.FALSE.

donotuseindexarray = .FALSE.
if (present(noindex)) then
  if (noindex.eqv..TRUE.) then
    donotuseindexarray = .TRUE.
  end if
end if

modality = trim(self%get_Modality())

select case(modality)
  case('EBSD')
    isEBSD = .TRUE.
  case('TKD')
    isTKD = .TRUE.
  case('ECP')
    isECP = .TRUE.
  case default
    call Message%printError('ang_writeFile', 'unknown name list type requested')
end select

! open the file (overwrite old one if it exists)
angname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(nml%angfile)
open(unit=dataunit2,file=trim(angname),status='unknown',action='write',iostat=ierr)

! this requires a lot of information...
write(dataunit2,'(A)') '# TEM_PIXperUM          1.000000'
s = ( float(nml%numsx)*0.5 + nml%xpc ) / float(nml%numsx)      ! x-star
write(dataunit2,'(A,F9.6)') '# x-star                ', s
s = ( float(nml%numsy)*0.5 + nml%ypc ) / float(nml%numsy)      ! y-star
write(dataunit2,'(A,F9.6)') '# y-star                ', s
s = nml%L / ( nml%delta * float(nml%numsx) )                   ! z-star
write(dataunit2,'(A,F9.6)') '# z-star                ', s
write(dataunit2,'(A,F9.6)') '# WorkingDistance       ', fpar(1)  ! WD ! this quantity is not used in EMsoft
write(dataunit2,'(A)') '#'
write(dataunit2,'(A)') '# Phase 1'

cellparams = cell%getLatParm()
SGnum = SG%getSpaceGroupNumber()
xtalname = trim(cell%getFileName())

ii = scan(trim(xtalname),'.')
angname = xtalname(1:ii-1)
write(dataunit2,'(A)') '# MaterialName    '//trim(angname)
write(dataunit2,'(A)') '# Formula       '//trim(angname)
str1 = '# Info          patterns indexed using EMsoft::EMDI; Diffraction Modality: '//trim(modality)
if (present(orthoset)) str1 = trim(str1)//'; orthorhombic space group setting '//extendedHMOrthsymbols(orthoset,orthoSG)
write(dataunit2,'(A)') trim(str1)

! and get the TSL symmetry string from the TSLsymtype array
TSLsymmetry = TSLsymtype(SG%getPGnumber())

! symmetry string
write(dataunit2,'(A)') '# Symmetry              '//TSLsymmetry

! lattice parameters
cellparams(1:3) = cellparams(1:3)*10.0  ! convert to Angstrom
if (present(orthoset)) then   ! permute the lattice parameters to the correct orthorhombic setting 
! settings  " a  b  c", " b  a -c", " c  a  b", "-c  b  a", " b  c  a", " a -c  b" 
  select case(orthoset)
    case(1) 
      a = cellparams(1)
      b = cellparams(2)
      c = cellparams(3)
    case(2) 
      a = cellparams(2)
      b = cellparams(1)
      c = cellparams(3)
    case(3) 
      a = cellparams(3)
      b = cellparams(1)
      c = cellparams(2)
    case(4) 
      a = cellparams(3)
      b = cellparams(2)
      c = cellparams(1)
    case(5) 
      a = cellparams(2)
      b = cellparams(3)
      c = cellparams(1)
    case(6) 
      a = cellparams(1)
      b = cellparams(3)
      c = cellparams(2)
  end select 
  write(str1,'(F8.3)') a
  write(str2,'(F8.3)') b
  write(str3,'(F8.3)') c
else
  write(str1,'(F8.3)') cellparams(1)
  write(str2,'(F8.3)') cellparams(2)
  write(str3,'(F8.3)') cellparams(3)
end if 
str1 = adjustl(str1)
str2 = adjustl(str2)
str3 = adjustl(str3)
str1 = trim(str1)//' '//trim(str2)//' '//trim(str3)

! unit cell angles
write(str4,'(F8.3)') cellparams(4)
write(str5,'(F8.3)') cellparams(5)
write(str6,'(F8.3)') cellparams(6)
str4 = adjustl(str4)
str5 = adjustl(str5)
str6 = adjustl(str6)
str1 = trim(str1)//TAB//trim(str4)//' '//trim(str5)//' '//trim(str6)

write(dataunit2,'(A)') '# LatticeConstants      '//trim(str1)
!==========================

! next we need to get the hklFamilies ranked by kinematical intensity, going out to some value
! this is probably not necessary [based on Stuart's feedback]
write(dataunit2,'(A)') '# NumberFamilies        0'

!==========================
write(dataunit2,'(A)') '# Categories 0 0 0 0 0'
write(dataunit2,'(A)') '#'
write(dataunit2,'(A)') '# GRID: SqrGrid'
write(dataunit2,'(A,F12.6)') '# XSTEP: ', nml%StepX
write(dataunit2,'(A,F12.6)') '# YSTEP: ', nml%StepY
write(dataunit2,'(A,I5)') '# NCOLS_ODD: ',ipar(7)
write(dataunit2,'(A,I5)') '# NCOLS_EVEN: ',ipar(7)
write(dataunit2,'(A,I5)') '# NROWS: ', ipar(8)
write(dataunit2,'(A)') '#'
write(dataunit2,'(A,A)') '# OPERATOR:   ', trim(EMsoft%getConfigParameter('UserName'))
write(dataunit2,'(A)') '#'
write(dataunit2,'(A)') '# SAMPLEID:'
write(dataunit2,'(A)') '#'
write(dataunit2,'(A)') '# SCANID:'
write(dataunit2,'(A)') '#'

! ok, next we have the actual data, which is in the following order
! * phi1                      -> Phi1
! * phi                       -> Phi
! * phi2                      -> Phi2
! * x pos                     -> pixel position
! * y pos                     -> pixel position
! * image quality             -> iq
! * confidence index          -> resultmain
! * phase                     -> 1 (starts at 1 for multi-phase file)
! the second entry after the arrow is the EMsoft parameter that we write into that location
! these 8 entries must be present; they represent Version 3 of the EDAX/TSL .ang format 

! go through the entire array and write one line per sampling point
do ii = 1,ipar(3)
    BSval = 255.0 * IQmap(ii)
! should we use the index array or not?
    if (donotuseindexarray.eqv..TRUE.) then
      indx = 0
      euler = eulerarray(1:3,ii)
    else
      indx = indexmain(1,ii)
      euler = eulerarray(1:3,indx)
    end if
    write(str1,'(A,F8.5)') ' ',euler(1)*dtor
    write(str2,'(A,F8.5)') ' ',euler(2)*dtor
    write(str3,'(A,F8.5)') ' ',euler(3)*dtor
! sampling coordinates [interchanged x and y on 05/28/19, MDG]
    if (sum(nml%ROI).ne.0) then
      write(str4,'(A,F12.5)') ' ',float(MODULO(ii-1,nml%ROI(3)))*nml%StepX
      write(str5,'(A,F12.5)') ' ',float(floor(float(ii-1)/float(nml%ROI(3))))*nml%StepY
    else
      write(str4,'(A,F12.5)') ' ',float(MODULO(ii-1,nml%ipf_wd))*nml%StepX
      write(str5,'(A,F12.5)') ' ',float(floor(float(ii-1)/float(nml%ipf_wd)))*nml%StepY
    end if
! Image Quality (using the Krieger Lassen pattern sharpness parameter iq)
    if (resultmain(1,ii).eq.0.0) BSval = 0.0  ! to prevent NaN values in certain cases
    write(str6,'(A,F6.1)') ' ',BSval  !  IQ value in range [0.0 .. 255.0]
    write(str7,'(A,F6.3)') ' ',resultmain(1,ii)   ! this replaces MAD
    write(str8,'(A,I1)') '  ',1
!
    write(dataunit2,"(A,' ',A,' ',A,' ',A,' ',A,' ',A,' ',A,' ',A)") trim(adjustl(str1)),trim(adjustl(str2)),&
                                            trim(adjustl(str3)),trim(adjustl(str4)),trim(adjustl(str5)),&
                                            trim(adjustl(str6)),trim(adjustl(str7)),trim(adjustl(str8))
end do

close(dataunit2,status='keep')

end subroutine ang_writeFile_


!--------------------------------------------------------------------------
recursive subroutine ctfmerge_writeFile_(self,EMsoft,cells,SGs,nml,ipar,fpar,eangles,phaseID,dplist,OSMlist,IQmap)
!DEC$ ATTRIBUTES DLLEXPORT :: ctfmerge_writeFile_
!! author: MDG
!! version: 1.0
!! date: 04/01/20
!!
!! Write a *.ctf output file with multiphase EBSD/TKD data (HKL format)

use mod_EMsoft
use mod_crystallography
use mod_symmetry
use mod_io
use mod_DIfiles

IMPLICIT NONE

class(Vendor_T),INTENT(INOUT)                       :: self
type(EMsoft_T),INTENT(INOUT)                        :: EMsoft
integer(kind=irg),INTENT(IN)                        :: ipar(4)
type(cell_T),INTENT(INOUT)                          :: cells(ipar(2))
type(SpaceGroup_T),INTENT(INOUT)                    :: SGs(ipar(2))
class(DictionaryIndexingNameListType),INTENT(INOUT) :: nml
!f2py intent(in,out) ::  nml
real(kind=sgl),INTENT(IN)                           :: fpar(2)
real(kind=sgl),INTENT(INOUT)                        :: eangles(3,ipar(1),ipar(2))
!f2py intent(in,out) ::  eangles
integer(kind=irg),INTENT(IN)                        :: phaseID(ipar(1))
real(kind=sgl),INTENT(IN)                           :: dplist(ipar(1),ipar(2))
real(kind=sgl),INTENT(IN)                           :: OSMlist(ipar(1),ipar(2))
real(kind=sgl),INTENT(IN)                           :: IQmap(ipar(1))

type(IO_T)                                          :: Message
integer(kind=irg)                                   :: ierr, i, ii, indx, hdferr, SGnum(5), LaueGroup(5), BCval, BSval, &
                                                       iph, numph
character(fnlen)                                    :: ctfname, xtalname, modality
character                                           :: TAB = CHAR(9)
character(fnlen)                                    :: str1,str2,str3,str4,str5,str6,str7,str8,str9,str10
character(1)                                        :: np
real(kind=sgl)                                      :: euler(3), eu, mi, ma
logical                                             :: donotuseindexarray, isEBSD=.FALSE., isTKD=.FALSE., isECP=.FALSE.
real(kind=dbl)                                      :: cellparams(6)
integer(kind=irg),allocatable                       :: osm(:), iq(:)
real(kind=sgl),allocatable                          :: osmr(:)

modality = trim(self%get_Modality())

select case(modality)
  case('EBSD')
    isEBSD = .TRUE.
  case('TKD')
    isTKD = .TRUE.
  case('ECP')
    isECP = .TRUE.
  case default
    call Message%printError('ctfmerge_writeFile_', 'unknown name list type requested')
end select

numph = ipar(2)

! get the OSMmap into 1D format and scale to the range [0..255]
allocate(osm(ipar(1)), osmr(ipar(1)))
do i=1,ipar(1)
  osmr(i) = OSMlist(i,phaseID(i))
end do
mi = minval(osmr)
ma = maxval(osmr)
osm = nint(255.0 * (osmr-mi)/(ma-mi))
deallocate(osmr)

! scale the IQmap to the range [0..255]
allocate(iq(ipar(1)))
iq = nint(255.0 * IQmap)

! open the file (overwrite old one if it exists)
ctfname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(nml%ctffile)
open(unit=dataunit2,file=trim(ctfname),status='unknown',action='write',iostat=ierr)

write(dataunit2,'(A)') 'Channel Text File'
write(dataunit2,'(A)') 'EMsoft v. '//trim(EMsoft%getConfigParameter('EMsoftversion'))// &
                       '; BANDS=pattern index, MAD=CI, BC=OSM, BS=IQ; Diffraction Modality: '//trim(modality)
write(dataunit2,'(A)') 'Author  '//trim(EMsoft%getConfigParameter('Username'))
write(dataunit2,'(A)') 'JobMode Grid'
write(dataunit2,'(2A,I5)') 'XCells',TAB, ipar(3)
write(dataunit2,'(2A,I5)') 'YCells',TAB, ipar(4)
write(dataunit2,'(2A,F6.2)') 'XStep',TAB, nml%StepX
write(dataunit2,'(2A,F6.2)') 'YStep',TAB, nml%StepY
write(dataunit2,'(A)') 'AcqE1'//TAB//'0'
write(dataunit2,'(A)') 'AcqE2'//TAB//'0'
write(dataunit2,'(A)') 'AcqE3'//TAB//'0'
write(dataunit2,'(A,A,$)') 'Euler angles refer to Sample Coordinate system (CS0)!',TAB
str1 = 'Mag'//TAB//'30'//TAB//'Coverage'//TAB//'100'//TAB//'Device'//TAB//'0'//TAB//'KV'
write(str2,'(F4.1)') fpar(1)  ! EkeV
str1 = trim(str1)//TAB//trim(str2)//TAB//'TiltAngle'
write(str2,'(F5.2)') fpar(2)  ! MCsig
str2 = adjustl(str2)
str1 = trim(str1)//TAB//trim(str2)//TAB//'TiltAxis'//TAB//'0'
write(dataunit2,'(A)') trim(str1)
write(np,"(I1)") ipar(2)
write(dataunit2,'(A)') 'Phases'//TAB//np

do iph=1,numph
  cellparams = cells(iph)%getLatParm()
  SGnum(iph) = SGs(iph)%getSpaceGroupNumber()
  xtalname = trim(cells(iph)%getFileName())

  ! unit cell size
  cellparams(1:3) = cellparams(1:3)*10.0  ! convert to Angstrom
  write(str1,'(F8.3)') cellparams(1)
  write(str2,'(F8.3)') cellparams(2)
  write(str3,'(F8.3)') cellparams(3)
  str1 = adjustl(str1)
  str2 = adjustl(str2)
  str3 = adjustl(str3)
  str1 = trim(str1)//';'//trim(str2)//';'//trim(str3)

  ! unit cell angles
  write(str4,'(F8.3)') cellparams(4)
  write(str5,'(F8.3)') cellparams(5)
  write(str6,'(F8.3)') cellparams(6)
  str4 = adjustl(str4)
  str5 = adjustl(str5)
  str6 = adjustl(str6)
  str1 = trim(str1)//TAB//trim(str4)//';'//trim(str5)//';'//trim(str6)

  ! structure name
  str3 = ''
  ii = len(trim(xtalname))-5
  do i=1,ii
    str3(i:i) = xtalname(i:i)
  end do
  str1 = trim(str1)//TAB//trim(str3)

  ! rotational symmetry group
  str4 = ''
  LaueGroup(iph) = SGs(iph)%getLaueGroupNumber()
  write(str4,'(I2)') LaueGroup(iph)
  str1 = trim(str1)//TAB//trim(adjustl(str4))

  ! space group
  str2 = ''
  write(str2,'(I3)') SGnum(iph)
  str1 = trim(str1)//TAB//trim(adjustl(str2))

  ! and now collect them all into a single string
  write(dataunit2,'(A)') trim(str1)
end do

! this is the table header
write(dataunit2,'(A)') 'Phase'//TAB//'X'//TAB//'Y'//TAB//'Bands'//TAB//'Error'//TAB//'Euler1'//TAB//'Euler2'//TAB//'Euler3' &
                      //TAB//'MAD'//TAB//'BC'//TAB//'BS'

! Euler angles are always in degrees
eangles = eangles * rtod

! go through the entire array and write one line per sampling point
do ii = 1,ipar(1)
    BCval = osm(ii)
    BSval = iq(ii)
    euler = eangles(1:3,ii,phaseID(ii))

! changed order of coordinates to conform with ctf standard
    if (sum(nml%ROI).ne.0) then
      write(str2,'(F12.3)') float(floor(float(ii-1)/float(nml%ROI(3))))*nml%StepY
      write(str1,'(F12.3)') float(MODULO(ii-1,nml%ROI(3)))*nml%StepX
    else
      write(str2,'(F12.3)') float(floor(float(ii-1)/float(nml%ipf_wd)))*nml%StepY
      write(str1,'(F12.3)') float(MODULO(ii-1,nml%ipf_wd))*nml%StepX
    end if

    write(str3,'(I8)') 0
    write(str8,'(I8)') 0 ! integer zero error; was indx, which is now moved to BANDS
    eu = euler(1) - 90.0 ! conversion from TSL to Oxford convention
    if (eu.lt.0) eu = eu + 360.0
    write(str5,'(F12.3)') eu
    eu = euler(2)
    if (eu.lt.0) eu = eu + 360.0
    write(str6,'(F12.3)') eu
! intercept the hexagonal case, for which we need to subtract 30° from the third Euler angle
! Note: after working with Lionel Germain, we concluded that we do not need to subtract 30°
! in the ctf file, because the fundamental zone is already oriented according to the Oxford
! convention... That means that we need to subtract the angle for the .ang file (to be implemented)
! [modified by MDG on 3/5/18]
    if (( LaueGroup(phaseID(ii)).eq.8 ).or.( LaueGroup(phaseID(ii)).eq.9 ) ) euler(3) = euler(3) - 30.0
    eu = euler(3)
    if (eu.lt.0) eu = eu + 360.0
    write(str7,'(F12.3)') eu
    write(str4,'(F12.6)') dplist(ii,phaseID(ii))   ! this replaces MAD
! the following two parameters need to be modified to contain more meaningful information
    write(str9,'(I8)') BCval   ! OSM value in range [0..255]
    write(str10,'(I8)') BSval  !  IQ value in range [0..255]
! some Oxford 3D files have four additional integer columns;
! GrainIndex
! GrainRandomColourR
! GrainRandomColourG
! GrainRandomColourB
!
    write(np,"(I1)") phaseID(ii)
    write(dataunit2,'(A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A)') np,TAB,trim(adjustl(str1)),TAB,&
    trim(adjustl(str2)),TAB,trim(adjustl(str3)),TAB,trim(adjustl(str8)),TAB,trim(adjustl(str5)),&
    TAB,trim(adjustl(str6)),TAB,trim(adjustl(str7)),TAB,trim(adjustl(str4)),TAB,trim(adjustl(str9)),&
    TAB,trim(adjustl(str10))
end do

close(dataunit2,status='keep')

! reset the Euler angles to radians
eangles = eangles * dtor

end subroutine ctfmerge_writeFile_

!--------------------------------------------------------------------------
recursive subroutine angmerge_writeFile_(self,EMsoft,cells,SGs,nml,ipar,fpar,eangles,phaseID,dplist,IQmap)
!DEC$ ATTRIBUTES DLLEXPORT :: angmerge_writeFile_
!! author: MDG
!! version: 1.0
!! date: 04/01/20
!!
!! Write a *.ang output file with multiphase EBSD/TKD data (TSL format)

use mod_EMsoft
use mod_crystallography
use mod_symmetry
use mod_io
use mod_DIfiles

IMPLICIT NONE

class(Vendor_T),INTENT(INOUT)                       :: self
type(EMsoft_T),INTENT(INOUT)                        :: EMsoft
integer(kind=irg),INTENT(IN)                        :: ipar(4)
type(cell_T),INTENT(INOUT)                          :: cells(ipar(2))
type(SpaceGroup_T),INTENT(INOUT)                    :: SGs(ipar(2))
class(DictionaryIndexingNameListType),INTENT(INOUT) :: nml
real(kind=sgl),INTENT(INOUT)                        :: fpar(1)
real(kind=sgl),INTENT(IN)                           :: eangles(3,ipar(1),ipar(2))
integer(kind=irg),INTENT(IN)                        :: phaseID(ipar(1))
real(kind=sgl),INTENT(IN)                           :: dplist(ipar(1),ipar(2))
real(kind=sgl),INTENT(IN)                           :: IQmap(ipar(1))

type(IO_T)                                          :: Message
integer(kind=irg)                                   :: ierr, ii, indx, SGnum, iph, numph
character(fnlen)                                    :: angname, xtalname, modality
character(fnlen)                                    :: str1,str2,str3,str4,str5,str6,str7,str8,str9,str10
character                                           :: TAB = CHAR(9)
character(2)                                        :: TSLsymmetry
character(1)                                        :: np
real(kind=sgl)                                      :: euler(3), s, BSval
real(kind=dbl)                                      :: cellparams(6)
logical                                             :: donotuseindexarray, isEBSD=.FALSE., isTKD=.FALSE., isECP=.FALSE.

modality = trim(self%get_Modality())

select case(modality)
  case('EBSD')
    isEBSD = .TRUE.
  case('TKD')
    isTKD = .TRUE.
  case('ECP')
    isECP = .TRUE.
  case default
    call Message%printError('angmerge_writeFile_', 'unknown name list type requested')
end select

numph = ipar(2)

! open the file (overwrite old one if it exists)
angname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(nml%angfile)
open(unit=dataunit2,file=trim(angname),status='unknown',action='write',iostat=ierr)

! this requires a lot of information...
write(dataunit2,'(A)') '# TEM_PIXperUM          1.000000'
s = ( float(nml%numsx)*0.5 + nml%xpc ) / float(nml%numsx)      ! x-star
write(dataunit2,'(A,F9.6)') '# x-star                ', s
s = ( float(nml%numsy)*0.5 + nml%ypc ) / float(nml%numsy)      ! y-star
write(dataunit2,'(A,F9.6)') '# y-star                ', s
s = nml%L / ( nml%delta * float(nml%numsx) )                   ! z-star
write(dataunit2,'(A,F9.6)') '# z-star                ', s
write(dataunit2,'(A,F9.6)') '# WorkingDistance       ', fpar(1) ! WD ! this quantity is not used in EMsoft
write(dataunit2,'(A)') '#'

do iph=1,ipar(2)
  write (np,"(I1)") iph
  write(dataunit2,'(A)') '# Phase '//np

  cellparams = cells(iph)%getLatParm()
  SGnum = SGs(iph)%getSpaceGroupNumber()
  xtalname = trim(cells(iph)%getFileName())

  ii = scan(trim(xtalname),'.')
  angname = xtalname(1:ii-1)
  write(dataunit2,'(A)') '# MaterialName    '//trim(angname)
  write(dataunit2,'(A)') '# Formula       '//trim(angname)
  write(dataunit2,'(A)') '# Info          patterns indexed using EMsoft::EMDI; Diffraction Modality: '//trim(modality)

  ! and get the TSL symmetry string from the TSLsymtype array
  TSLsymmetry = TSLsymtype(SGs(iph)%getPGnumber())

  ! symmetry string
  write(dataunit2,'(A)') '# Symmetry              '//TSLsymmetry

  ! lattice parameters
  cellparams(1:3) = cellparams(1:3)*10.0  ! convert to Angstrom
  write(str1,'(F8.3)') cellparams(1)
  write(str2,'(F8.3)') cellparams(2)
  write(str3,'(F8.3)') cellparams(3)
  str1 = adjustl(str1)
  str2 = adjustl(str2)
  str3 = adjustl(str3)
  str1 = trim(str1)//' '//trim(str2)//' '//trim(str3)

  ! unit cell angles
  write(str4,'(F8.3)') cellparams(4)
  write(str5,'(F8.3)') cellparams(5)
  write(str6,'(F8.3)') cellparams(6)
  str4 = adjustl(str4)
  str5 = adjustl(str5)
  str6 = adjustl(str6)
  str1 = trim(str1)//TAB//trim(str4)//' '//trim(str5)//' '//trim(str6)

  write(dataunit2,'(A)') '# LatticeConstants      '//trim(str1)
  !==========================

  ! next we need to get the hklFamilies ranked by kinematical intensity, going out to some value
  ! this is probably not necessary [based on Stuart's feedback]
  write(dataunit2,'(A)') '# NumberFamilies        0'
end do
!==========================
! write(dataunit2,'(A)') '# Categories 0 0 0 0 0'
! write(dataunit2,'(A)') '#'
write(dataunit2,'(A)') '# GRID: SqrGrid'
write(dataunit2,'(A,F12.6)') '# XSTEP: ', nml%StepX
write(dataunit2,'(A,F12.6)') '# YSTEP: ', nml%StepY
write(dataunit2,'(A,I5)') '# NCOLS_ODD: ',ipar(3)
write(dataunit2,'(A,I5)') '# NCOLS_EVEN: ',ipar(3)
write(dataunit2,'(A,I5)') '# NROWS: ', ipar(4)
write(dataunit2,'(A)') '#'
write(dataunit2,'(A,A)') '# OPERATOR:   ', trim(EMsoft%getConfigParameter('Username'))
write(dataunit2,'(A)') '#'
write(dataunit2,'(A)') '# SAMPLEID:'
write(dataunit2,'(A)') '#'
write(dataunit2,'(A)') '# SCANID:'
write(dataunit2,'(A)') '#'

! ok, next we have the actual data, which is in the following order
! * phi1                      -> Phi1
! * phi                       -> Phi
! * phi2                      -> Phi2
! * x pos                     -> pixel position
! * y pos                     -> pixel position
! * image quality             -> iq
! * confidence index          -> resultmain
! * phase                     -> 0 (since there is only one phase in each indexing run)
! the second entry after the arrow is the EMsoft parameter that we write into that location
! these 8 entries must be present; they represent Version 3 of the EDAX/TSL .ang format 

! go through the entire array and write one line per sampling point
do ii = 1,ipar(3)*ipar(4)
    write (np,"(I1)") phaseID(ii)
    BSval = 255.0 * IQmap(ii)
    euler = eangles(1:3,ii,phaseID(ii))
    write(str1,'(A,F8.5)') ' ',euler(1)
    write(str2,'(A,F8.5)') ' ',euler(2)
    write(str3,'(A,F8.5)') ' ',euler(3)
! sampling coordinates [interchanged x and y on 05/28/19, MDG]
    if (sum(nml%ROI).ne.0) then
      write(str4,'(A,F12.5)') ' ',float(MODULO(ii-1,nml%ROI(3)))*nml%StepX
      write(str5,'(A,F12.5)') ' ',float(floor(float(ii-1)/float(nml%ROI(3))))*nml%StepY
    else
      write(str4,'(A,F12.5)') ' ',float(MODULO(ii-1,nml%ipf_wd))*nml%StepX
      write(str5,'(A,F12.5)') ' ',float(floor(float(ii-1)/float(nml%ipf_wd)))*nml%StepY
    end if
! Image Quality (using the Krieger Lassen pattern sharpness parameter iq)
    write(str6,'(A,F6.1)') ' ',BSval  !  IQ value in range [0.0 .. 255.0]
    write(str7,'(A,F6.3)') ' ',dplist(ii,phaseID(ii))   ! this replaces MAD
    write(str8,'(A)') '  '//np
!
    write(dataunit2,"(A,' ',A,' ',A,' ',A,' ',A,' ',A,' ',A,' ',A)") trim(adjustl(str1)),trim(adjustl(str2)),&
                                            trim(adjustl(str3)),trim(adjustl(str4)),trim(adjustl(str5)),&
                                            trim(adjustl(str6)),trim(adjustl(str7)),trim(adjustl(str8))
end do

close(dataunit2,status='keep')

end subroutine angmerge_writeFile_

end module mod_vendors
