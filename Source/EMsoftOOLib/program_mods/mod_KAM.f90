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

module mod_KAM
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/27/24
  !!
  !! class definition for the EMKAM program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMKAM program
type, public :: KAMNameListType
  real(kind=sgl)    :: kamrange(2)
  integer(kind=irg) :: orav
  character(fnlen)  :: dotproductfile
  character(fnlen)  :: kamtiff
end type KAMNameListType

! class definition
type, public :: KAM_T
private 
  character(fnlen)       :: nmldeffile = 'EMKAM.nml'
  type(KAMNameListType)  :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: KAM_
  procedure, pass(self) :: setkamrange_
  procedure, pass(self) :: getkamrange_
  procedure, pass(self) :: setorav_
  procedure, pass(self) :: getorav_
  procedure, pass(self) :: setdotproductfile_
  procedure, pass(self) :: getdotproductfile_
  procedure, pass(self) :: setkamtiff_
  procedure, pass(self) :: getkamtiff_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: KAM => KAM_
  generic, public :: setkamcutoff => setkamrange_
  generic, public :: getkamcutoff => getkamrange_
  generic, public :: setorav => setorav_
  generic, public :: getorav => getorav_
  generic, public :: setdotproductfile => setdotproductfile_
  generic, public :: getdotproductfile => getdotproductfile_
  generic, public :: setkamtiff => setkamtiff_
  generic, public :: getkamtiff => getkamtiff_

end type KAM_T

! the constructor routine for this class 
interface KAM_T
  module procedure KAM_constructor
end interface KAM_T

contains

!--------------------------------------------------------------------------
type(KAM_T) function KAM_constructor( nmlfile ) result(KAM)
!DEC$ ATTRIBUTES DLLEXPORT :: KAM_constructor
!! author: MDG 
!! version: 1.0 
!! date: 12/27/24
!!
!! constructor for the KAM_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call KAM%readNameList(nmlfile)

end function KAM_constructor

!--------------------------------------------------------------------------
subroutine KAM_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 12/27/24
!!
!! destructor for the KAM_T Class
 
IMPLICIT NONE

type(KAM_T), INTENT(INOUT)  :: self 

call reportDestructor('KAM_T')

end subroutine KAM_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 12/27/24
!!
!! read the namelist from an nml file for the KAM_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(KAM_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

real(kind=sgl)                       :: kamrange(2)
integer(kind=irg)                    :: orav
character(fnlen)                     :: dotproductfile
character(fnlen)                     :: kamtiff

! define the IO namelist to facilitate passing variables to the program.
namelist  / getKAM / kamrange, orav, dotproductfile, kamtiff

! set the input parameters to default values
kamrange = (/ 0.0, 5.0 /)
orav = 20
dotproductfile = 'undefined'
KAMtiff = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=getKAM)
    close(UNIT=dataunit,STATUS='keep')

! check for required entries
    if (trim(dotproductfile).eq.'undefined') then
        call Message%printError('readNameList:',' dot product file name is undefined in '//nmlfile)
    end if

    if (trim(KAMtiff).eq.'undefined') then
        call Message%printError('readNameList:',' output tiff file name is undefined in '//nmlfile)
    end if
 end if

! if we get here, then all appears to be ok, and we need to fill in the enl fields
self%nml%kamrange = kamrange
self%nml%orav = orav
self%nml%dotproductfile = dotproductfile
self%nml%KAMtiff = KAMtiff

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 12/27/24
!!
!! pass the namelist for the KAM_T Class to the calling program

IMPLICIT NONE 

class(KAM_T), INTENT(INOUT)          :: self
type(KAMNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine setkamrange_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setkamrange_
!! author: MDG
!! version: 1.0
!! date: 12/27/24
!!
!! set kamrange in the KAM_T class

IMPLICIT NONE

class(KAM_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)      :: inp(2)

self%nml%kamrange = inp

end subroutine setkamrange_

!--------------------------------------------------------------------------
function getkamrange_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getkamrange_
!! author: MDG
!! version: 1.0
!! date: 12/27/24
!!
!! get kamrange from the KAM_T class

IMPLICIT NONE

class(KAM_T), INTENT(INOUT)     :: self
real(kind=sgl)                  :: out(2)

out = self%nml%kamrange

end function getkamrange_

!--------------------------------------------------------------------------
subroutine setorav_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setorav_
!! author: MDG
!! version: 1.0
!! date: 12/27/24
!!
!! set orav in the KAM_T class

IMPLICIT NONE

class(KAM_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%orav = inp

end subroutine setorav_

!--------------------------------------------------------------------------
function getorav_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getorav_
!! author: MDG
!! version: 1.0
!! date: 12/27/24
!!
!! get orav from the KAM_T class

IMPLICIT NONE

class(KAM_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%orav

end function getorav_

!--------------------------------------------------------------------------
subroutine setdotproductfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdotproductfile_
!! author: MDG
!! version: 1.0
!! date: 12/27/24
!!
!! set dotproductfile in the KAM_T class

IMPLICIT NONE

class(KAM_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%dotproductfile = trim(inp)

end subroutine setdotproductfile_

!--------------------------------------------------------------------------
function getdotproductfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdotproductfile_
!! author: MDG
!! version: 1.0
!! date: 12/27/24
!!
!! get dotproductfile from the KAM_T class

IMPLICIT NONE

class(KAM_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%dotproductfile)

end function getdotproductfile_

!--------------------------------------------------------------------------
subroutine setkamtiff_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setkamtiff_
!! author: MDG
!! version: 1.0
!! date: 12/27/24
!!
!! set kamtiff in the KAM_T class

IMPLICIT NONE

class(KAM_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%kamtiff = trim(inp)

end subroutine setkamtiff_

!--------------------------------------------------------------------------
function getkamtiff_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getkamtiff_
!! author: MDG
!! version: 1.0
!! date: 12/27/24
!!
!! get kamtiff from the KAM_T class

IMPLICIT NONE

class(KAM_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%kamtiff)

end function getkamtiff_

!--------------------------------------------------------------------------
subroutine KAM_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: KAM_
!! author: MDG 
!! version: 1.0 
!! date: 12/27/24
!!
!! perform the computations

use mod_EMsoft
use mod_io
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_DIsupport
use mod_DIfiles
use mod_rotations 
use mod_quaternions
use ISO_C_BINDING
use mod_image
use, intrinsic :: iso_fortran_env

IMPLICIT NONE 

class(KAM_T), INTENT(INOUT)             :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 

type(HDF_T)                             :: HDF
type(HDFnames_T)                        :: HDFnames
type(IO_T)                              :: Message
type(DIfile_T)                          :: DIFT
type(DictionaryIndexingNameListType)    :: dinl
! type(QuaternionArray_T)                 :: qAR, sym
! type(e_T)                               :: e 
! type(q_T)                               :: q 
! type(Quaternion_T)                      :: qu

logical                                 :: stat, readonly, noindex
integer(kind=irg)                       :: hdferr, nlines, FZcnt, Nexp, nnm, nnk, Pmdims, i, j, k, olabel, Nd, Ne, ipar(10), &
                                           ipar2(6), pgnum, ipat, ipf_wd, ipf_ht, idims2(2), io_int(2), TIFF_nx, TIFF_ny
character(fnlen)                        :: groupname, dataset, fname, TIFF_filename, DIfile
integer(HSIZE_T)                        :: dims2(2)
real(kind=sgl)                          :: testeu(3), eu1(3), eu2(3), ro1(4), ro2(4), s, da

character(fnlen),allocatable            :: stringarray(:)
integer(kind=irg),allocatable           :: tmi(:,:), tmitmp(:,:), indexmain(:,:)
real(kind=sgl),allocatable              :: kam(:,:), eulers(:,:), Eulerstmp(:,:), Eulervals(:,:), avEuler(:,:), &
                                           dplist(:,:), dplisttmp(:,:)

! declare variables for use in object oriented image module
integer                                 :: iostat
character(len=128)                      :: iomsg
logical                                 :: isInteger
type(image_t)                           :: im
integer(int8)                           :: i8 (3,4)
integer(int8), allocatable              :: TIFF_image(:,:)

associate(kamnl=>self%nml, dinl=>DIFT%nml, DIDT=>DIFT%DIDT)

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()
HDFnames = HDFnames_T()

call HDFnames%set_NMLfiles(SC_NMLfiles)
call HDFnames%set_NMLfilename(SC_DictionaryIndexingNML)
call HDFnames%set_NMLparameters(SC_NMLparameters)
call HDFnames%set_NMLlist(SC_DictionaryIndexingNameListType)

DIfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(kamnl%dotproductfile)
call DIFT%readDotProductFile(EMsoft, HDF, HDFnames, DIfile, hdferr, &
                             getEulerAngles=.TRUE., &
                             getRefinedEulerAngles=.TRUE., &
                             getRefinedDotProducts = .TRUE., &
                             getTopDotProductList=.TRUE., &
                             getTopMatchIndices=.TRUE.) 

dinl = DIFT%getNameList()

if (sum(dinl%ROI).ne.0) then 
  ipf_wd = dinl%ROI(3)
  ipf_ht = dinl%ROI(4)
else
  ipf_ht = dinl%ipf_ht
  ipf_wd = dinl%ipf_wd
end if 
DIDT%Nexp = ipf_wd*ipf_ht
FZcnt = DIDT%FZcnt 

allocate(eulers(3,DIDT%Nexp), aveuler(3,DIDT%Nexp))
if (allocated(DIDT%RefinedEulerAngles)) then 
  eulers = DIDT%RefinedEulerAngles
else
  eulers = DIDT%EulerAngles
end if

if (maxval(eulers).gt.2.0*cPi) eulers = eulers*dtor

call setRotationPrecision('d')
call setCheckValidity(.FALSE.)

! qAR = QuaternionArray_T( DIDT%Nexp, s = 'd')

! if (allocated(DIDT%RefinedEulerAngles)) then 
!   do i=1,DIDT%Nexp 
!    e = e_T( edinp = dble(DIDT%RefinedEulerAngles(1:3,i)) )
!    q = e%eq()
!    qu = quaternion_T( qd = q%q_copyd() )
!    call qAR%insertQuatinArray( i, qu )
!   end do 
!   deallocate(DIDT%RefinedEulerAngles)
! else 
!   do i=1,DIDT%Nexp 
!    e = e_T( edinp = dble(DIDT%EulerAngles(1:3,i)) )
!    q = e%eq()
!    qu = quaternion_T( qd = q%q_copyd() )
!    call qAR%insertQuatinArray( i, qu )
!   end do 
!   deallocate(DIDT%EulerAngles)
! end if 

! FZcnt = qAR%getListCount()

allocate(dplist(dinl%nnk,DIDT%Nexp))
do i=1,DIDT%Nexp
  dplist(1:dinl%nnk,i) = DIDT%TopDotProductList(1:dinl%nnk,i)
end do
deallocate(DIDT%TopDotProductList)
nnm = dinl%nnk

allocate(tmi(nnm,DIDT%Nexp))
do i=1,DIDT%Nexp
  tmi(1:nnm,i) = DIDT%TopMatchIndices(1:nnm,i)
end do
deallocate(DIDT%TopMatchIndices)

! and next we compute the KAM map (kam)
allocate(kam(ipf_wd,ipf_ht))

! do we need to do an orientation average first ?
if (kamnl%orav.ne.0) then
  ipar2(1) = DIDT%pgnum
  ipar2(2) = FZcnt
  ipar2(3) = DIDT%Nexp
  ipar2(4) = dinl%nnk
  ipar2(5) = DIDT%Nexp * ceiling(float(ipf_wd*ipf_ht)/float(DIDT%Nexp))
  ! to average we need at least two values so check the value of orav
  if (kamnl%orav.eq.1) then
    ipar2(6) = 2
  else
    ipar2(6) = kamnl%orav
  end if 
  call Message%printMessage('Computing orientation averages ... ')
  call DIgetAverageOrientations(ipar2, eulers, tmi, dplist, avEuler)
  ! eulers = eulers*sngl(cPi)/180.0
end if

call getKAMMap(DIDT%Nexp, eulers, ipf_wd, ipf_ht, DIDT%pgnum, kam)

! put the KAM map in degrees and set any values outside the KAMrange to zero
kam = kam*180.0/sngl(cPi)

write (*,*) ' KAM range : ', minval(kam), maxval(kam)

where (kam.lt.kamnl%kamrange(1)) kam = 0.0
where (kam.gt.kamnl%kamrange(2)) kam = 0.0
kam = kam-minval(kam)
kam = kam/maxval(kam)
kam = kam*255.0

! and finally we generate the output arrays, depending on what the user asked for
TIFF_filename = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(kamnl%kamtiff)
TIFF_nx = ipf_wd
TIFF_ny = ipf_ht

! allocate memory for image
allocate(TIFF_image(TIFF_nx,TIFF_ny))

TIFF_image = int(kam)

! set up the image_t structure
im = image_t(TIFF_image)
if(im%empty()) call Message%printMessage("EMgetKAM","failed to convert array to image")

! create the file
call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage(" failed to write image to file : "//iomsg)
else
  call Message%printMessage(' KAM map written to '//trim(TIFF_filename))
end if

end associate

end subroutine KAM_

end module mod_KAM
