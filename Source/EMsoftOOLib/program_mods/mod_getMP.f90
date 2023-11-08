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

module mod_getMP
  !! author: MDG
  !! version: 1.0
  !! date: 04/12/20
  !!
  !! class definition for the EMgetMP program

use mod_kinds
use mod_global

IMPLICIT NONE

! namelist for the EMgetMP program
type, public :: getMPNameListType
  character(fnlen)  :: masterfile
  character(fnlen)  :: projectionmode
  character(fnlen)  :: outputfile
  logical           :: ratioimage
end type getMPNameListType

! class definition
type, public :: getMP_T
private
  character(fnlen)          :: nmldeffile = 'EMgetMP.nml'
  type(getMPNameListType)   :: nml

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: getMP_
  procedure, pass(self) :: get_masterfile_
  procedure, pass(self) :: get_projectionmode_
  procedure, pass(self) :: get_outputfile_
  procedure, pass(self) :: set_masterfile_
  procedure, pass(self) :: set_projectionmode_
  procedure, pass(self) :: set_outputfile_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: getMP => getMP_
  generic, public :: get_masterfile => get_masterfile_
  generic, public :: get_projectionmode => get_projectionmode_
  generic, public :: get_outputfile => get_outputfile_
  generic, public :: set_masterfile => set_masterfile_
  generic, public :: set_projectionmode => set_projectionmode_
  generic, public :: set_outputfile => set_outputfile_
end type getMP_T

! the constructor routine for this class
interface getMP_T
  module procedure getMP_constructor
end interface getMP_T

contains

!--------------------------------------------------------------------------
type(getMP_T) function getMP_constructor( nmlfile ) result(getMP)
!DEC$ ATTRIBUTES DLLEXPORT :: getMP_constructor
!! author: MDG
!! version: 1.0
!! date: 04/12/20
!!
!! constructor for the getMP_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call getMP%readNameList(nmlfile)

end function getMP_constructor

!--------------------------------------------------------------------------
subroutine getMP_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: getMP_destructor
!! author: MDG
!! version: 1.0
!! date: 04/12/20
!!
!! destructor for the getMP_T Class

IMPLICIT NONE

type(getMP_T), INTENT(INOUT)  :: self

call reportDestructor('getMP_T')

end subroutine getMP_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 04/12/20
!!
!! read the namelist from an nml file for the getMP_T Class

use mod_io
use mod_EMsoft

IMPLICIT NONE

class(getMP_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft
type(IO_T)                           :: Message
logical                              :: skipread = .FALSE.

character(fnlen)  :: masterfile
character(fnlen)  :: projectionmode
character(fnlen)  :: outputfile
logical           :: ratioimage

namelist / getMPlist / masterfile, projectionmode, outputfile, ratioimage

masterfile = 'undefined'
projectionmode = 'stereographic'
outputfile = 'undefined'
ratioimage = .FALSE.

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=getMPlist)
    close(UNIT=dataunit,STATUS='keep')

    if (trim(masterfile).eq.'undefined') then
        call Message%printError('readNameList:',' master pattern file name is undefined in '//nmlfile)
    end if

    if (trim(outputfile).eq.'undefined') then
        call Message%printError('readNameList:',' output file name is undefined in '//nmlfile)
    end if
end if

self%nml%masterfile = trim(masterfile)
self%nml%projectionmode = trim(projectionmode)
self%nml%outputfile = trim(outputfile)
self%nml%ratioimage = ratioimage

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 04/12/20
!!
!! pass the namelist for the getMP_T Class to the calling program

IMPLICIT NONE

class(getMP_T), INTENT(INOUT)          :: self
type(getMPNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
function get_masterfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_masterfile_
!! author: MDG
!! version: 1.0
!! date: 04/12/20
!!
!! get masterfile from the getMP_T class

IMPLICIT NONE

class(getMP_T), INTENT(INOUT)     :: self
character(fnlen)                  :: out

out = self%nml%masterfile

end function get_masterfile_

!--------------------------------------------------------------------------
subroutine set_masterfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_masterfile_
!! author: MDG
!! version: 1.0
!! date: 04/12/20
!!
!! set masterfile in the getMP_T class

IMPLICIT NONE

class(getMP_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)      :: inp

self%nml%masterfile = inp

end subroutine set_masterfile_

!--------------------------------------------------------------------------
function get_projectionmode_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_projectionmode_
!! author: MDG
!! version: 1.0
!! date: 04/12/20
!!
!! get projectionmode from the getMP_T class

IMPLICIT NONE

class(getMP_T), INTENT(INOUT)     :: self
character(fnlen)                  :: out

out = self%nml%projectionmode

end function get_projectionmode_

!--------------------------------------------------------------------------
subroutine set_projectionmode_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_projectionmode_
!! author: MDG
!! version: 1.0
!! date: 04/12/20
!!
!! set projectionmode in the getMP_T class

IMPLICIT NONE

class(getMP_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)      :: inp

self%nml%projectionmode = inp

end subroutine set_projectionmode_

!--------------------------------------------------------------------------
function get_outputfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_outputfile_
!! author: MDG
!! version: 1.0
!! date: 04/12/20
!!
!! get outputfile from the getMP_T class

IMPLICIT NONE

class(getMP_T), INTENT(INOUT)     :: self
character(fnlen)                  :: out

out = self%nml%outputfile

end function get_outputfile_

!--------------------------------------------------------------------------
subroutine set_outputfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_outputfile_
!! author: MDG
!! version: 1.0
!! date: 04/12/20
!!
!! set outputfile in the getMP_T class

IMPLICIT NONE

class(getMP_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)      :: inp

self%nml%outputfile = inp

end subroutine set_outputfile_

!--------------------------------------------------------------------------
subroutine getMP_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: getMP_
!! author: MDG
!! version: 1.0
!! date: 04/12/20
!!
!! perform the computations

use mod_EMsoft
use mod_MCfiles
use mod_MPfiles
use mod_HDFsupport
use HDF5
use mod_HDFnames
use mod_math
use mod_io
use mod_image
use mod_Lambert
use stringconstants
use ISO_C_BINDING
use, intrinsic :: iso_fortran_env

IMPLICIT NONE

class(getMP_T), INTENT(INOUT)           :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname

type(IO_T)                              :: Message
type(MCfile_T)                          :: MCFT
type(MPfile_T)                          :: MPFT
type(HDF_T)                             :: HDF
type(HDFnames_T)                        :: HDFnames
type(Lambert_T)                         :: L
type(SEMmasterNameListType)             :: mpnl

character(fnlen)                        :: fname, modality, TIFF_filename1, TIFF_filename2, TIFF_filename3
integer(kind=irg)                       :: i, j, ierr, n, d, TIFF_nx, TIFF_ny
real(kind=sgl),allocatable              :: mLPNH(:,:), mLPSH(:,:), weights(:), masterSPNH(:,:), masterSPSH(:,:), ratio(:,:)
real(kind=sgl)                          :: avNH, avSH, mi, ma, xyzs(3)
real(kind=dbl)                          :: xyz(3), Radius

! declare variables for use in object oriented image module
integer                                 :: iostat
character(len=128)                      :: iomsg
logical                                 :: isInteger
type(image_t)                           :: im
integer(int8)                           :: i8 (3,4)
integer(int8), allocatable              :: TIFF_image1(:,:), TIFF_image2(:,:), TIFF_ratio(:,:)


associate(enl=>self%nml, MCDT=>MCFT%MCDT, MPDT=>MPFT%MPDT)

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()

! set the HDF group names
HDFnames = HDFnames_T()

! get the data file modality
fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfile))
call MPFT%determineModality(HDF, fname)
modality = trim(MPFT%getModality())
call Message%printMessage(' Master Pattern modality : '//trim(modality))

! 1. read Monte Carlo data so that we can compute the energy-weighted master pattern
call HDFnames%set_ProgramData(SC_MCOpenCL)
call HDFnames%set_NMLlist(SC_MCCLNameList)
call HDFnames%set_NMLfilename(SC_MCOpenCLNML)
fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfile))
call MCFT%setFileName(fname)
call MCFT%readMCfile(HDF, HDFnames, getAccume=.TRUE.)

! 2. read the master pattern file
if (trim(modality).eq.'TKD') then
  call HDFnames%set_ProgramData(SC_TKDmaster)
  call HDFnames%set_NMLlist(SC_TKDmasterNameList)
  call HDFnames%set_NMLfilename(SC_TKDmasterNML)
end if
if (trim(modality).eq.'EBSD') then
  call HDFnames%set_ProgramData(SC_EBSDmaster)
  call HDFnames%set_NMLlist(SC_EBSDmasterNameList)
  call HDFnames%set_NMLfilename(SC_EBSDmasterNML)
end if
if (trim(modality).eq.'ECP') then
  call HDFnames%set_ProgramData(SC_ECPmaster)
  call HDFnames%set_NMLlist(SC_ECPmasterNameList)
  call HDFnames%set_NMLfilename(SC_ECPmasterNML)
end if
call HDFnames%set_Variable(SC_MCOpenCL)

fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfile))
call MPFT%setFileName(fname)
call MPFT%readMPfile(HDF, HDFnames, mpnl, &
                     getmLPNH=.TRUE., &
                     getmLPSH=.TRUE.)

! 3. build the energy averaged master pattern
! ! first make sure that we do some appropriate energy weighting to make the master patterns 2D instead of 3D
! ! sum MC counts over rectangle [-1/2,1/2] and [-1,-1/3] in square Lambert, then normalize and use as
! ! weight factors, then store in 2D MPs
n = size(MCDT%accum_e, 1) ! get number of energy bins
allocate(weights(n)) ! allocate space for energy histogram
do i = 1, n
  weights(i) = sum(MCDT%accum_e(i,:,:)) ! this could be modified to sum over partial rectangle
enddo
weights = weights / sum(weights) ! this is currently wieghted over the full square Lambert

d = mpnl%npx
allocate(mLPNH(-d:d,-d:d))
allocate(mLPSH(-d:d,-d:d))
mLPNH = 0.D0
mLPSH = 0.D0
do i = 1, n
  mLPNH = mLPNH + MPDT%mLPNH(:,:,i) * weights(i)
  mLPSH = mLPSH + MPDT%mLPSH(:,:,i) * weights(i)
enddo
deallocate(MCDT%accum_e, weights, MPDT%mLPNH, MPDT%mLPSH)

! and finally we generate the output arrays, depending on what the user asked for
fname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(enl%outputfile)
TIFF_nx = 2*d+1
TIFF_ny = 2*d+1
! allocate memory for image
allocate(TIFF_image1(TIFF_nx,TIFF_ny),TIFF_image2(TIFF_nx,TIFF_ny),TIFF_ratio(TIFF_nx,TIFF_ny))

if (trim(enl%projectionmode).eq.'modifiedLambert') then
! this is the standard mode for master patterns, so no conversion necessary
  TIFF_filename1 = trim(fname) //'_mLPNH.tiff'
  ma = maxval(mLPNH)
  mi = minval(mLPNH)
  TIFF_image1 = int(255 * (mLPNH - mi)/ (ma-mi))

  TIFF_filename2 = trim(fname) //'_mLPSH.tiff'
  ma = maxval(mLPSH)
  mi = minval(mLPSH)
  TIFF_image2 = int(255 * (mLPSH - mi)/ (ma-mi))

! the ratio image can be useful if the structure is non-centrosymmetric 
  if (enl%ratioimage.eqv..TRUE.) then 
    TIFF_filename3 = trim(fname) //'_ratio.tiff'
    mLPNH = mLPNH / mLPSH 
    ma = maxval(mLPNH)
    mi = minval(mLPNH)
    TIFF_ratio = int(255 * (mLPNH - mi)/ (ma-mi))
  end if 
else
  if (trim(enl%projectionmode).eq.'stereographic') then
    TIFF_filename1 = trim(fname) //'_SPNH.tiff'
    TIFF_filename2 = trim(fname) //'_SPSH.tiff'
    Radius = 1.D0
  else  ! regular Lambert
    TIFF_filename1 = trim(fname) //'_LPNH.tiff'
    TIFF_filename2 = trim(fname) //'_LPSH.tiff'
    Radius = 1.D0/sqrt(2.D0)
  end if

! we need to set the area outside the projection circles to the average
! value of the master pattern array to avoid having contrast issues
  avNH = sum(mLPNH) / float(2*d+1)**2
  avSH = sum(mLPSH) / float(2*d+1)**2

! transform the mLPNH and mLPSH arrays to stereographic/Lambert projections
  allocate(masterSPNH(-d:d,-d:d),masterSPSH(-d:d,-d:d))
  do i=-d,d
    do j=-d,d
      if ((i*i+j*j).gt.d*d) then
        masterSPNH(i,j) = avNH
        masterSPSH(i,j) = avSH
      else
        L = Lambert_T( xyd = (/ dble(i), dble(j) /) / dble(d) )
        if (trim(enl%projectionmode).eq.'stereographic') then
          ierr = L%StereoGraphicInverse( xyz, Radius )
        else
          ierr = L%LambertInverse( xyz, Radius )
        end if
        if (ierr.eq.0) then
          xyzs = sngl(xyz/vecnorm(xyz))
          masterSPNH(i,j) = InterpolateLambert(xyzs, mLPNH, d)
          masterSPSH(i,j) = InterpolateLambert(xyzs, mLPSH, d)
        end if
      end if
    end do
  end do

! and convert to the output arrays
  ma = maxval(masterSPNH)
  mi = minval(masterSPNH)
  TIFF_image1 = int(255 * (masterSPNH - mi)/ (ma-mi))

  ma = maxval(masterSPSH)
  mi = minval(masterSPSH)
  TIFF_image2 = int(255 * (masterSPSH - mi)/ (ma-mi))

! the ratio image can be useful if the structure is non-centrosymmetric 
  if (enl%ratioimage.eqv..TRUE.) then 
    TIFF_filename3 = trim(fname) //'_ratio.tiff'
    masterSPNH = masterSPNH / masterSPSH 
    ma = maxval(masterSPNH)
    mi = minval(masterSPNH)
    TIFF_ratio = int(255 * (masterSPNH - mi)/ (ma-mi))
  end if 

  deallocate(masterSPSH, masterSPNH)
end if
! set up the image_t structure
im = image_t(TIFF_image1)
if(im%empty()) call Message%printMessage("EMgetMP","failed to convert array to image")

! create the file
call im%write(trim(TIFF_filename1), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage("failed to write image to file : "//iomsg)
else
  call Message%printMessage(' projection written to '//trim(TIFF_filename1))
end if

im = image_t(TIFF_image2)
if(im%empty()) call Message%printMessage("EMgetMP","failed to convert array to image")

! create the file
call im%write(trim(TIFF_filename2), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage("failed to write image to file : "//iomsg)
else
  call Message%printMessage(' projection written to '//trim(TIFF_filename2))
end if

if (enl%ratioimage.eqv..TRUE.) then 
  im = image_t(TIFF_ratio)
  if(im%empty()) call Message%printMessage("EMgetMP","failed to convert array to image")

! create the file
  call im%write(trim(TIFF_filename3), iostat, iomsg) ! format automatically detected from extension
  if(0.ne.iostat) then
    call Message%printMessage("failed to write image to file : "//iomsg)
  else
    call Message%printMessage(' projection written to '//trim(TIFF_filename3))
  end if
end if

deallocate(TIFF_image1, TIFF_image2, TIFF_ratio)

end associate

end subroutine getMP_



end module mod_getMP
