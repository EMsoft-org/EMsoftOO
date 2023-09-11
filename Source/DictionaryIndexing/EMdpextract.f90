! ###################################################################
! Copyright (c) 2019-2023, Marc De Graef Research Group/Carnegie Mellon University
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
program EMdpextract
!! author: MDG 
!! version: 1.0 
!! date: 04/05/20
!!
!! extract a number of images from a dot product file and store them in a folder

use mod_kinds 
use mod_global
use mod_EMsoft
use HDF5
use mod_HDFsupport
use mod_io
use mod_DIfiles

IMPLICIT NONE

character(fnlen)      :: progname = 'EMdpextract.f90' 
character(fnlen)      :: progdesc = 'Extract a number of images from a dictionary indexing dot product file'

type(EMsoft_T)        :: EMsoft 
type(IO_T)            :: Message 

character(fnlen)      :: dpfilebase, cwd
integer(kind=irg)     :: another, io_int(1), numarg, i, hdferr

! print the EMsoft header and handle any command line arguments  
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 907 /) )

! open the fortran HDF interface
call openFortranHDFInterface()

! get the current folder 
call getcwd(cwd)

! number of command line arguments
numarg = command_argument_count()

if (numarg.gt.0) then
  do i=1,numarg
    call get_command_argument(i,dpfilebase)
    call extractImages(EMsoft, dpfilebase, cwd)
  end do
else ! no command line arguments, so we work interactively 
  another = 1

  do while (another.eq.1) 
    ! ask the user for the dot product file name without the .h5 extension
    dpfilebase = ''
    call Message%ReadValue('Enter the name of the dot product file without the .h5 extension:', dpfilebase)

    call extractImages(EMsoft, dpfilebase, cwd)

    call Message%printMessage('----')
    call Message%ReadValue(' Another one ? (1/0) :', io_int)
    another = io_int(1)

  end do 

end if

! close the fortran HDF interface
call closeFortranHDFInterface()

end program EMdpextract

!--------------------------------------------------------------------------
subroutine extractImages(EMsoft, dpfilebase, cwd)
!! author: MDG 
!! version: 1.0 
!! date: 04/05/20
!!
!! Subroutine to extract images from a dot product file

use mod_EMsoft
use HDF5
use h5im
use h5lt
use mod_HDFsupport
use mod_HDFnames
use mod_io
use mod_image
use, intrinsic :: iso_fortran_env
use mod_DIfiles
use mod_platformsupport

IMPLICIT NONE 

type(EMsoft_T),INTENT(INOUT)  :: EMsoft
character(fnlen),INTENT(IN)   :: dpfilebase 
character(fnlen),INTENT(IN)   :: cwd

type(IO_T)                    :: Message
type(DIfile_T)                :: DIFT 
type(HDF_T)                   :: HDF 
type(HDFnames_T)              :: HDFnames

character(fnlen)              :: dirname, image_filename, dpfile
integer(kind=irg)             :: hdferr, shp(2), jj, nx, ny, status
real(kind=sgl)                :: mi, ma
logical                       :: fexists

! declare variables for use in object oriented image module
integer                       :: iostat
character(len=128)            :: iomsg
logical                       :: isInteger
type(image_t)                 :: im, im2
integer(int8)                 :: i8 (3,4), int8val
integer(int8), allocatable    :: output_image(:,:)

! read all the image-type arrays from the file
dpfile = trim(cwd)//'/'//trim(dpfilebase)//'.h5'
call Message%printMessage(' looking for file '//trim(dpfile))

HDF = HDF_T() 
HDFnames = HDFnames_T() 

call HDFnames%set_NMLfiles(SC_NMLfiles)
call HDFnames%set_NMLfilename(SC_DictionaryIndexingNML)
call HDFnames%set_NMLparameters(SC_NMLparameters)
call HDFnames%set_NMLlist(SC_DictionaryIndexingNameListType)

call DIFT%readDotProductFile(EMsoft, HDF, HDFnames, dpfile, hdferr, &
                             getRefinedDotProducts=.TRUE., &
                             getADP=.TRUE., &
                             getKAM=.TRUE., &
                             getCI=.TRUE., &
                             getIQ=.TRUE., & 
                             getOSM=.TRUE. )

associate(DIDT=>DIFT%DIDT)

! take these arrays and generate image files for each of them; place them in a folder with the dpfilebase name 
dirname = trim(dpfilebase)
inquire(file=trim(dirname),exist=fexists)
if (.not.(fexists)) then
  status = system_system('mkdir '//trim(dirname))
  call Message%printMessage(' '//trim(dirname)//' folder has been created')
end if
status = system_chdir(trim(dirname))

! ==============================
! ==============================
! ==============================
! ADP map
image_filename = trim(dpfilebase)//'_ADP.tiff'
shp = shape(DIDT%ADP)
nx = shp(1)
ny = shp(2)
allocate(output_image(nx,ny))
output_image = DIDT%ADP
im2 = image_t(output_image)
if(im2%empty()) call Message%printMessage("EMEBSDDIpreview","failed to convert array to image")

! this data set is already scaled to byte range

! create the file
call im2%write(trim(image_filename), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage(" failed to write image to file : "//iomsg)
else  
  call Message%printMessage('  ADP array written to '//trim(image_filename))
end if 
deallocate(output_image,DIDT%ADP)

! ==============================
! ==============================
! ==============================
! CI map
image_filename = trim(dpfilebase)//'_CI.tiff'
allocate(output_image(nx,ny))

mi = minval(DIDT%CI)
DIDT%CI = DIDT%CI - mi
ma = maxval(DIDT%CI)
output_image = 0

! CI is a 1-D array; map onto 2-D array
do jj = 1,ny
    output_image(1:nx,jj) = int(255.0*DIDT%CI((jj-1)*nx+1:jj*nx)/ma)
end do

im2 = image_t(output_image)
if(im2%empty()) call Message%printMessage("EMEBSDDIpreview","failed to convert array to image")

! create the file
call im2%write(trim(image_filename), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage(" failed to write image to file : "//iomsg)
else  
  call Message%printMessage('  CI array written to '//trim(image_filename))
end if 
deallocate(output_image,DIDT%CI)

! ==============================
! ==============================
! ==============================
! refined CI map
if (allocated(DIDT%RefinedDotProducts)) then 
  image_filename = trim(dpfilebase)//'_CIrefined.tiff'
  allocate(output_image(nx,ny))

  mi = minval(DIDT%RefinedDotProducts)
  DIDT%RefinedDotProducts = DIDT%RefinedDotProducts - mi
  ma = maxval(DIDT%RefinedDotProducts)
  output_image = 0

  ! CI is a 1-D array; map onto 2-D array
  do jj = 1,ny
      output_image(1:nx,jj) = int(255.0*DIDT%RefinedDotProducts((jj-1)*nx+1:jj*nx)/ma)
  end do

  im2 = image_t(output_image)
  if(im2%empty()) call Message%printMessage("EMEBSDDIpreview","failed to convert array to image")

  ! create the file
  call im2%write(trim(image_filename), iostat, iomsg) ! format automatically detected from extension
  if(0.ne.iostat) then
    call Message%printMessage(" failed to write image to file : "//iomsg)
  else  
    call Message%printMessage('  CI array written to '//trim(image_filename))
  end if 
  deallocate(output_image,DIDT%RefinedDotProducts)
end if


! ==============================
! ==============================
! ==============================
! IQ map
image_filename = trim(dpfilebase)//'_IQ.tiff'
allocate(output_image(nx,ny))

mi = minval(DIDT%IQ)
DIDT%IQ = DIDT%IQ - mi
ma = maxval(DIDT%IQ)
output_image = 0

! IQ is a 1-D array; map onto 2-D array
do jj = 1,ny
    output_image(1:nx,jj) = int(255.0*DIDT%IQ((jj-1)*nx+1:jj*nx)/ma)
end do

im2 = image_t(output_image)
if(im2%empty()) call Message%printMessage("EMEBSDDIpreview","failed to convert array to image")

! create the file
call im2%write(trim(image_filename), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage(" failed to write image to file : "//iomsg)
else  
  call Message%printMessage('  IQ array written to '//trim(image_filename))
end if 
deallocate(output_image,DIDT%IQ)

! ==============================
! ==============================
! ==============================
! KAM map
image_filename = trim(dpfilebase)//'_KAM.tiff'
allocate(output_image(nx,ny))

mi = minval(DIDT%KAM)
DIDT%KAM = DIDT%KAM - mi
ma = maxval(DIDT%KAM)
output_image = 0

output_image = int(255.0*DIDT%KAM/ma)

im2 = image_t(output_image)
if(im2%empty()) call Message%printMessage("EMEBSDDIpreview","failed to convert array to image")

! create the file
call im2%write(trim(image_filename), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage(" failed to write image to file : "//iomsg)
else  
  call Message%printMessage('  KAM array written to '//trim(image_filename))
end if 
deallocate(output_image,DIDT%KAM)

! ==============================
! ==============================
! ==============================
! OSM map
image_filename = trim(dpfilebase)//'_OSM.tiff'
allocate(output_image(nx,ny))

mi = minval(DIDT%OSM)
DIDT%OSM = DIDT%OSM - mi
ma = maxval(DIDT%OSM)
output_image = 0

output_image = int(255.0*DIDT%OSM/ma)

im2 = image_t(output_image)
if(im2%empty()) call Message%printMessage("EMEBSDDIpreview","failed to convert array to image")

! create the file
call im2%write(trim(image_filename), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage(" failed to write image to file : "//iomsg)
else  
  call Message%printMessage('  OSM array written to '//trim(image_filename))
end if 
deallocate(output_image,DIDT%OSM)

status = system_chdir(trim(cwd))

end associate 

end subroutine extractImages
