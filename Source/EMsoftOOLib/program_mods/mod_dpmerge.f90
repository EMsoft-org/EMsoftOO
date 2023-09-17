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

module mod_dpmerge
  !! author: MDG
  !! version: 1.0
  !! date: 04/07/20
  !!
  !! class definition for the EMdpmerge program

use mod_kinds
use mod_global

IMPLICIT NONE

! namelist for the EMdpmerge program
type, public :: dpmergeNameListType
  character(fnlen)  :: dotproductfile(5)
  character(fnlen)  :: ctfname
  character(fnlen)  :: angname
  character(fnlen)  :: phasemapname
  integer(kind=irg) :: phasecolors(5)
  integer(kind=irg) :: scaling
  real(kind=sgl)    :: scalefactors(5)
  character(8)      :: usedp
  character(2)      :: indexingmode
end type dpmergeNameListType

! class definition
type, public :: dpmerge_T
private
  character(fnlen)            :: nmldeffile = 'EMdpmerge.nml'
  type(dpmergeNameListType)   :: nml

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: dpmerge_
  procedure, pass(self) :: get_dotproductfile_
  procedure, pass(self) :: get_ctfname_
  procedure, pass(self) :: get_angname_
  procedure, pass(self) :: get_phasemapname_
  procedure, pass(self) :: get_phasecolors_
  procedure, pass(self) :: get_scalefactors_
  procedure, pass(self) :: get_scaling_
  procedure, pass(self) :: get_usedp_
  procedure, pass(self) :: get_indexingmode_
  procedure, pass(self) :: set_dotproductfile_
  procedure, pass(self) :: set_ctfname_
  procedure, pass(self) :: set_angname_
  procedure, pass(self) :: set_phasemapname_
  procedure, pass(self) :: set_phasecolors_
  procedure, pass(self) :: set_scaling_
  procedure, pass(self) :: set_scalefactors_
  procedure, pass(self) :: set_usedp_
  procedure, pass(self) :: set_indexingmode_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: dpmerge => dpmerge_
  generic, public :: get_dotproductfile => get_dotproductfile_
  generic, public :: get_ctfname => get_ctfname_
  generic, public :: get_angname => get_angname_
  generic, public :: get_phasemapname => get_phasemapname_
  generic, public :: get_phasecolors => get_phasecolors_
  generic, public :: get_scalefactors => get_scalefactors_
  generic, public :: get_scaling => get_scaling_
  generic, public :: get_usedp => get_usedp_
  generic, public :: get_indexingmode => get_indexingmode_
  generic, public :: set_dotproductfile => set_dotproductfile_
  generic, public :: set_ctfname => set_ctfname_
  generic, public :: set_angname => set_angname_
  generic, public :: set_phasemapname => set_phasemapname_
  generic, public :: set_phasecolors => set_phasecolors_
  generic, public :: set_scaling => set_scaling_
  generic, public :: set_scalefactors => set_scalefactors_
  generic, public :: set_usedp => set_usedp_
  generic, public :: set_indexingmode => set_indexingmode_

end type dpmerge_T

! the constructor routine for this class
interface dpmerge_T
  module procedure dpmerge_constructor
end interface dpmerge_T

contains

!--------------------------------------------------------------------------
type(dpmerge_T) function dpmerge_constructor( nmlfile ) result(dpmerge)
!DEC$ ATTRIBUTES DLLEXPORT :: dpmerge_constructor
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! constructor for the dpmerge_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call dpmerge%readNameList(nmlfile)

end function dpmerge_constructor

!--------------------------------------------------------------------------
subroutine dpmerge_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: dpmerge_destructor
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! destructor for the dpmerge_T Class

IMPLICIT NONE

type(dpmerge_T), INTENT(INOUT)  :: self

call reportDestructor('dpmerge_T')

end subroutine dpmerge_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! read the namelist from an nml file for the dpmerge_T Class

use mod_io
use mod_EMsoft

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)      :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft
type(IO_T)                           :: Message
logical                              :: skipread = .FALSE.

character(fnlen)        :: dotproductfile(5)
character(fnlen)        :: ctfname
character(fnlen)        :: angname
character(fnlen)        :: phasemapname
integer(kind=irg)       :: phasecolors(5)
integer(kind=irg)       :: scaling
real(kind=sgl)          :: scalefactors(5)  
character(8)            :: usedp
character(2)            :: indexingmode

! define the IO namelist to facilitate passing variables to the program.
namelist  / dpmerge / dotproductfile, ctfname, angname, usedp, indexingmode, phasemapname, phasecolors, scalefactors, scaling

! set the input parameters to default values
dotproductfile = (/ 'undefined','undefined','undefined','undefined','undefined' /)
ctfname = 'undefined'
angname = 'undefined'
phasemapname = 'undefined'
phasecolors = (/ 1, 2, 0, 0, 0 /)
scalefactors = (/ 1.0, 1.0, 1.0, 1.0, 1.0 /)
scaling = 1
usedp = 'original'
indexingmode = 'DI'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=dpmerge)
    close(UNIT=dataunit,STATUS='keep')

! check for required entries
    if ((trim(dotproductfile(1)).eq.'undefined').or.(trim(dotproductfile(2)).eq.'undefined')) then
        call Message%printError('readNameList:',' at least two dot product file names must be defined in '//nmlfile)
    end if

    if ((trim(ctfname).eq.'undefined').and.(trim(angname).eq.'undefined')) then
        call Message%printError('readNameList:',' either ctfname or angname must be defined in '//nmlfile)
    end if
 end if

! if we get here, then all appears to be ok, and we need to fill in the dpmnl fields
self%nml%dotproductfile = dotproductfile
self%nml%ctfname = ctfname
self%nml%angname = angname
self%nml%phasemapname = phasemapname
self%nml%phasecolors = phasecolors
self%nml%scalefactors = scalefactors
self%nml%scaling = scaling
self%nml%indexingmode = indexingmode
self%nml%usedp = usedp

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! pass the namelist for the dpmerge_T Class to the calling program

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)          :: self
type(dpmergeNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
function get_dotproductfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dotproductfile_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! get dotproductfile from the dpmerge_T class

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)     :: self
character(fnlen)                    :: out(5)

out = self%nml%dotproductfile

end function get_dotproductfile_

!--------------------------------------------------------------------------
subroutine set_dotproductfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_dotproductfile_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! set dotproductfile in the dpmerge_T class

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)        :: inp(5)

self%nml%dotproductfile = inp

end subroutine set_dotproductfile_

!--------------------------------------------------------------------------
function get_ctfname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ctfname_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! get ctfname from the dpmerge_T class

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)     :: self
character(fnlen)                    :: out

out = self%nml%ctfname

end function get_ctfname_

!--------------------------------------------------------------------------
subroutine set_ctfname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ctfname_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! set ctfname in the dpmerge_T class

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)        :: inp

self%nml%ctfname = inp

end subroutine set_ctfname_

!--------------------------------------------------------------------------
function get_angname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_angname_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! get angname from the dpmerge_T class

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)     :: self
character(fnlen)                    :: out

out = self%nml%angname

end function get_angname_

!--------------------------------------------------------------------------
subroutine set_angname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_angname_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! set angname in the dpmerge_T class

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)        :: inp

self%nml%angname = inp

end subroutine set_angname_

!--------------------------------------------------------------------------
function get_phasemapname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_phasemapname_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! get phasemapname from the dpmerge_T class

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)     :: self
character(fnlen)                    :: out

out = self%nml%phasemapname

end function get_phasemapname_

!--------------------------------------------------------------------------
subroutine set_phasemapname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_phasemapname_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! set phasemapname in the dpmerge_T class

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)        :: inp

self%nml%phasemapname = inp

end subroutine set_phasemapname_

!--------------------------------------------------------------------------
function get_phasecolors_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_phasecolors_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! get phasecolors from the dpmerge_T class

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out(5)

out = self%nml%phasecolors

end function get_phasecolors_

!--------------------------------------------------------------------------
subroutine set_phasecolors_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_phasecolors_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! set phasecolors in the dpmerge_T class

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp(5)

self%nml%phasecolors = inp

end subroutine set_phasecolors_

!--------------------------------------------------------------------------
function get_scalefactors_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_scalefactors_
!! author: MDG
!! version: 1.0
!! date: 11/17/22
!!
!! get scalefactors from the dpmerge_T class

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)     :: self
real(kind=sgl)                      :: out(5)

out = self%nml%scalefactors

end function get_scalefactors_

!--------------------------------------------------------------------------
subroutine set_scalefactors_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_scalefactors_
!! author: MDG
!! version: 1.0
!! date: 11/17/22
!!
!! set scalefactors in the dpmerge_T class

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)          :: inp(5)

self%nml%scalefactors = inp

end subroutine set_scalefactors_

!--------------------------------------------------------------------------
function get_scaling_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_scaling_
!! author: MDG
!! version: 1.0
!! date: 11/17/22
!!
!! get scaling from the dpmerge_T class

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%scaling

end function get_scaling_

!--------------------------------------------------------------------------
subroutine set_scaling_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_scaling_
!! author: MDG
!! version: 1.0
!! date: 11/17/22
!!
!! set scaling in the dpmerge_T class

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: inp

self%nml%scaling = inp

end subroutine set_scaling_

!--------------------------------------------------------------------------
function get_usedp_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_usedp_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! get usedp from the dpmerge_T class

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)     :: self
character(8)                        :: out

out = self%nml%usedp

end function get_usedp_

!--------------------------------------------------------------------------
subroutine set_usedp_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_usedp_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! set usedp in the dpmerge_T class

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)     :: self
character(8), INTENT(IN)            :: inp

self%nml%usedp = inp

end subroutine set_usedp_

!--------------------------------------------------------------------------
function get_indexingmode_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_indexingmode_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! get indexingmode from the dpmerge_T class

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)     :: self
character(2)                        :: out

out = self%nml%indexingmode

end function get_indexingmode_

!--------------------------------------------------------------------------
subroutine set_indexingmode_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_indexingmode_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! set indexingmode in the dpmerge_T class

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)     :: self
character(2), INTENT(IN)            :: inp

self%nml%indexingmode = inp

end subroutine set_indexingmode_

!--------------------------------------------------------------------------
subroutine dpmerge_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: dpmerge_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! perform the computations

use mod_EMsoft
use mod_io
use mod_crystallography
use mod_symmetry
use mod_HDFsupport
use mod_HDFnames
use mod_DIfiles
use mod_MPfiles
use mod_vendors
use ISO_C_BINDING
use mod_image

use, intrinsic :: iso_fortran_env

IMPLICIT NONE

class(dpmerge_T), INTENT(INOUT)         :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname

type(HDF_T)                             :: HDF
type(HDFnames_T)                        :: HDFnames, MPHDFnames
type(IO_T)                              :: Message
type(DIfile_T)                          :: DIFT
type(MPfile_T)                          :: MPFT
type(DictionaryIndexingNameListType)    :: dinl
type(Vendor_T)                          :: VT
type(cell_T),allocatable                :: cells(:)
type(SpaceGroup_T),allocatable          :: SGs(:)
type(EBSDmasterNameListType)            :: mpnl

real(kind=sgl),allocatable              :: dplist(:,:), OSMlist(:,:), exptIQ(:), eangles(:,:,:), pfrac(:), pID(:), dpmap(:,:)
integer(kind=irg),allocatable           :: phaseID(:), pnum(:)
integer(kind=irg)                       :: ipf_wd, ipf_ht, irow, numpat, ml(1), ipar(4)
integer(kind=irg)                       :: dims(1), hdferr, io_int(2), i, j, ii, numdp
real(kind=sgl)                          :: io_real(1), mi, ma, fpar1(1), fpar2(2)
character(fnlen)                        :: fname, xtalname(5), infile, rdxtalname, TIFF_filename, DIfile, modality
logical                                 :: f_exists
integer(kind=ish)                       :: imax = 255
integer(kind=ish),allocatable           :: dpmapi(:,:)

! declare variables for use in object oriented image module
integer                                 :: iostat
character(len=128)                      :: iomsg
logical                                 :: isInteger
type(image_t)                           :: im
integer(int8), allocatable              :: TIFF_image(:,:)
integer                                 :: dim2(2)
integer(c_int32_t)                      :: result


associate(dpmnl=>self%nml, DIDT=>DIFT%DIDT, MPDT=>MPFT%MPDT)

! how many input files are there ?
numdp = 0
do i=1,5
  if (trim(dpmnl%dotproductfile(i)).ne.'') numdp = numdp + 1
end do
io_int(1) = numdp
call Message%WriteValue('',io_int,1,"(' Found ',I2,' dot product file names'/)")

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()

HDFnames = HDFnames_T()
MPHDFnames = HDFnames_T()

call HDFnames%set_NMLfiles(SC_NMLfiles)
call HDFnames%set_NMLfilename(SC_DictionaryIndexingNML)
call HDFnames%set_NMLparameters(SC_NMLparameters)
call HDFnames%set_NMLlist(SC_DictionaryIndexingNameListType)

call MPHDFnames%set_NMLfiles(SC_NMLfiles)
call MPHDFnames%set_ProgramData(SC_EBSDmaster)
call MPHDFnames%set_NMLfilename(SC_EBSDmasterNML)
call MPHDFnames%set_NMLparameters(SC_NMLparameters)
call MPHDFnames%set_NMLlist(SC_EBSDmasterNameList)
modality = 'EBSD'

! if (dpmnl%indexingmode.eq.'DI') then
! loop over the input files and extract all the necessary data; at the same time, check to make
! sure that they cover the same data (after the first one has been read, allocate the data arrays)
  do i=1,numdp
    call Message%printMessage(' Reading data from '//trim(dpmnl%dotproductfile(i)) )
    DIfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(dpmnl%dotproductfile(i))
    if (trim(dpmnl%usedp).eq.'original') then ! read the original DI results
      call DIFT%readDotProductFile(EMsoft, HDF, HDFnames, DIfile, hdferr, &
                                   getEulerAngles = .TRUE., &
                                   getIQ = .TRUE., &
                                   getOSM = .TRUE., &
                                   getCI = .TRUE.)
    else   ! read the results from the refinement run
      call DIFT%readDotProductFile(EMsoft, HDF, HDFnames, DIfile, hdferr, &
                                   getRefinedEulerAngles = .TRUE., &
                                   getIQ = .TRUE., &
                                   getOSM = .TRUE., &
                                   getRefinedDotProducts = .TRUE., &
                                   getCI = .TRUE.)
    end if
    dinl = DIFT%getNameList()

    if (i.eq.1) then
  ! get the ROI dimensions and allocate the arrays
      if (sum(dinl%ROI).ne.0) then
        ipf_wd = dinl%ROI(3)
        ipf_ht = dinl%ROI(4)
      else
        ipf_wd = dinl%ipf_wd
        ipf_ht = dinl%ipf_ht
      end if
      numpat = ipf_wd * ipf_ht
      allocate( dplist(numpat,numdp), OSMlist(numpat, numdp), exptIQ(numpat), eangles(3,numpat,numdp) )
    else
  ! check dimensions of the ROI; they must be the same.  if they are, then add the data to the various arrays
      dims = shape(DIDT%CI)
      if (dims(1).ne.numpat) then
        call Message%printError('EMdpmerge','inconsistent ROI dimensions in dot product input files ' )
      end if
    end if
! copy the data
    do irow = 1, ipf_ht
      OSMlist( (irow-1)*ipf_wd+1:irow*ipf_wd,i) = DIDT%OSM(1:ipf_wd,irow)
    end do
    exptIQ(:) = DIDT%IQ(:)
    deallocate( DIDT%OSM, DIDT%IQ )
    if (trim(dpmnl%usedp).eq.'original') then
      eangles(1:3,1:numpat,i) = DIDT%EulerAngles(1:3,1:numpat)
      dplist(1:numpat,i) = DIDT%CI(1:numpat)
      deallocate( DIDT%EulerAngles, DIDT%CI )
    else
      eangles(1:3,1:numpat,i) = DIDT%RefinedEulerAngles(1:3,1:numpat)
      dplist(1:numpat,i) = DIDT%RefinedDotProducts(1:numpat)
      deallocate( DIDT%RefinedEulerAngles, DIDT%RefinedDotProducts )
    end if

! finally, get the name of the xtal file from the master pattern file
! if that file can not be found, ask the user interactively to enter the xtalname parameter
    infile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(dinl%masterfile)
    inquire(file=trim(infile), exist=f_exists)

    if (f_exists.eqv..TRUE.) then
      call MPFT%setFileName(infile)
      call MPFT%setModality(modality)
      call MPFT%readMPfile(HDF, MPHDFnames, mpnl)
      xtalname(i) = trim(MPDT%xtalname)
    else
      call Message%printMessage('***************************')
      call Message%printMessage('Master pattern file '//trim(infile)//' can not be found')
      call Message%printMessage('***************************')
      call Message%WriteValue(' Current dot product file :', trim(dpmnl%dotproductfile(i)))
      call Message%ReadValue('  Enter crystal structure file name (with extension) ', rdxtalname)
      xtalname(i) = trim(rdxtalname)
    end if

  end do

  ! apply the scale factors to each of the dot product arrays
  ! this can be used to slightly adjust one phase with respect to another
  do i=1,numdp
    dplist(1:numpat,i) = dplist(1:numpat,i) * dpmnl%scalefactors(i)
  end do

  ! determine which phase has the largest confidence index for each ROI sampling point
  allocate(phaseID(numpat), pID(numdp))
  do i=1,numpat
    pID(1:numdp) = dplist(i,1:numdp)
    ml = maxloc(pID)
    phaseID(i) = ml(1)
  end do
  deallocate(pID)

  ! determine the phase fractions and print that information
  allocate(pnum(numdp),pfrac(numdp))
  pnum = 0
  do i=1,numpat
    pnum(phaseID(i)) = pnum(phaseID(i)) + 1
  end do
  pfrac = float(pnum)/float(numpat) * 100.0
  call Message%printMessage(' Phase fractions :',"(/A)")
  call Message%printMessage(' -----------------')
  do i=1,numdp
    io_real(1) = pfrac(i)
    call Message%WriteValue('  Phase '//trim(xtalname(i)), io_real, 1, "(F6.2)")
  end do

! get the crystallographic and symmetry information for all phases.
  allocate(cells(numdp), SGs(numdp))
  do i=1,numdp
    call cells(i)%setFileName(xtalname(i))
    call cells(i)%readDataHDF(SGs(i), EMsoft, HDF)
  end do

! write a new .ctf or .ang file, if requested
  ipar(1) = numpat
  ipar(2) = numdp
  ipar(3) = ipf_wd
  ipar(4) = ipf_ht
  fpar1(1) = 10.0
  fpar2(1) = dinl%energymax
  fpar2(2) = DIDT%MCsig
  VT = Vendor_T()
  call VT%set_Modality(modality)

  if (trim(dpmnl%ctfname).ne.'undefined') then
    dinl%ctffile = trim(dpmnl%ctfname)
    call VT%ctfmerge_writeFile(EMsoft,cells,SGs,dinl,ipar,fpar2,eangles, phaseID, dplist, OSMlist, exptIQ)
    call Message%printMessage(' Merged orientation data stored in ctf file : '//trim(dpmnl%ctfname),"(/A)")
  end if

! write a new .ang file, if requested
  if (trim(dpmnl%angname).ne.'undefined') then
    dinl%angfile = trim(dpmnl%angname)
    call VT%angmerge_writeFile(EMsoft,cells,SGs,dinl,ipar,fpar1,eangles, phaseID, dplist, exptIQ)
    call Message%printMessage(' Merged orientation data stored in ang file : '//trim(dpmnl%angname),"(/A)")
  end if

! else ! indexing mode must be SI
! ! the files are Spherical Indexing files, so they do not have an OSM map in them, and some other
! ! things are different, so we need a somewhat different approach.


! end if

if (trim(dpmnl%phasemapname).ne.'undefined') then
  ! output the phase map as a tiff/bmp/ file
  fname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(dpmnl%phasemapname)
  TIFF_filename = trim(fname)

  ! allocate memory for a color image; each pixel has 3 bytes [RGB]
  allocate(TIFF_image(3*ipf_wd,ipf_ht))
  TIFF_image = 0_int8

  ! fill the image with whatever data you have (between 0 and 255)
  allocate( dpmap(ipf_wd, ipf_ht), dpmapi(ipf_wd, ipf_ht) )
  do j=1,ipf_ht
   do i=1,ipf_wd
    ii = (j-1) * ipf_wd + i
    dpmap(i,j) = dplist(ii,phaseID(ii))
   end do
  end do

  ma = maxval(dpmap)
  mi = minval(dpmap)
  dpmap = int(255*(dpmap-mi)/(ma-mi))

! the pre-defined colors are red, green, blue, yellow, cyan, fushia, and white
! if dpmnl$scaling=1 then each color is weighted by the maximum dot product 
! value to make the image a bit more realistic
  do j=1,ipf_ht
   do i=1,ipf_wd
    ii = (j-1) * ipf_wd + i
    if (dpmnl%scaling.eq.1) then 
      select case(dpmnl%phasecolors(phaseID(ii)))
        case(1)
          TIFF_image(1+3*(i-1),j) = dpmap(i,j)

        case(2)
          TIFF_image(2+3*(i-1),j) = dpmap(i,j)

        case(3)
          TIFF_image(3+3*(i-1),j) = dpmap(i,j)

        case(4)
          TIFF_image(1+3*(i-1),j) = dpmap(i,j)
          TIFF_image(2+3*(i-1),j) = dpmap(i,j)

        case(5)
          TIFF_image(2+3*(i-1),j) = dpmap(i,j)
          TIFF_image(3+3*(i-1),j) = dpmap(i,j)

        case(6)
          TIFF_image(1+3*(i-1),j) = dpmap(i,j)
          TIFF_image(3+3*(i-1),j) = dpmap(i,j)

        case(7)
          TIFF_image(1+3*(i-1),j) = dpmap(i,j)
          TIFF_image(2+3*(i-1),j) = dpmap(i,j)
          TIFF_image(3+3*(i-1),j) = dpmap(i,j)
      end select
    else
      select case(dpmnl%phasecolors(phaseID(ii)))
        case(1)
          TIFF_image(1+3*(i-1),j) = imax

        case(2)
          TIFF_image(2+3*(i-1),j) = imax

        case(3)
          TIFF_image(3+3*(i-1),j) = imax

        case(4)
          TIFF_image(1+3*(i-1),j) = imax
          TIFF_image(2+3*(i-1),j) = imax

        case(5)
          TIFF_image(2+3*(i-1),j) = imax
          TIFF_image(3+3*(i-1),j) = imax

        case(6)
          TIFF_image(1+3*(i-1),j) = imax
          TIFF_image(3+3*(i-1),j) = imax

        case(7)
          TIFF_image(1+3*(i-1),j) = imax
          TIFF_image(2+3*(i-1),j) = imax
          TIFF_image(3+3*(i-1),j) = imax
       end select

     end if
   end do
  end do

  ! set up the image_t structure
  im = image_t(TIFF_image)
  im%dims = (/ ipf_wd, ipf_ht /)
  im%samplesPerPixel = 3
  if(im%empty()) call Message%printMessage("EMdpmerge: failed to convert array to rgb image")

  ! create the file
  call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
  if(0.ne.iostat) then
    call Message%printMessage(" Failed to write image to file : "//iomsg)
  else
    call Message%printMessage(' Color phase map written to '//trim(TIFF_filename),"(/A)")
  end if
deallocate(TIFF_image)
end if

! close the fortran HDF interface
call closeFortranHDFInterface()

end associate

end subroutine dpmerge_



end module mod_dpmerge
