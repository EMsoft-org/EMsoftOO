! ###################################################################
! Copyright (c) 2013-2026, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_ANG
  !! author: MDG
  !! version: 1.0
  !! date: 04/06/20
  !!
  !! class definition for the EMgetANG program

use mod_kinds
use mod_global

IMPLICIT NONE

! namelist for the EMgetANG program
type, public :: ANGNameListType
  character(4)      :: modality
  character(8)      :: angledataset   ! 'original' or 'refined'
  logical           :: isHROSM
  character(fnlen)  :: xtalname
  character(fnlen)  :: newangfile
  character(fnlen)  :: dotproductfile
end type ANGNameListType

! class definition
type, public :: ANG_T
private
  character(fnlen)       :: nmldeffile = 'EMgetANG.nml'
  type(ANGNameListType)  :: nml

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: ANG_
  procedure, pass(self) :: get_mod_
  procedure, pass(self) :: get_isHROSM_
  procedure, pass(self) :: get_angledataset_
  procedure, pass(self) :: get_xtalname_
  procedure, pass(self) :: get_newangfile_
  procedure, pass(self) :: get_dotproductfile_
  procedure, pass(self) :: set_mod_
  procedure, pass(self) :: set_isHROSM_
  procedure, pass(self) :: set_angledataset_
  procedure, pass(self) :: set_xtalname_
  procedure, pass(self) :: set_newangfile_
  procedure, pass(self) :: set_dotproductfile_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: ANG => ANG_
  generic, public :: get_mod => get_mod_
  generic, public :: get_isHROSM => get_isHROSM_
  generic, public :: get_angledataset => get_angledataset_
  generic, public :: get_xtalname => get_xtalname_
  generic, public :: get_newangfile => get_newangfile_
  generic, public :: get_dotproductfile => get_dotproductfile_
  generic, public :: set_mod => set_mod_
  generic, public :: set_isHROSM => set_isHROSM_
  generic, public :: set_angledataset => set_angledataset_
  generic, public :: set_xtalname => set_xtalname_
  generic, public :: set_newangfile => set_newangfile_
  generic, public :: set_dotproductfile => set_dotproductfile_
end type ANG_T

! the constructor routine for this class
interface ANG_T
  module procedure ANG_constructor
end interface ANG_T

contains

!--------------------------------------------------------------------------
type(ANG_T) function ANG_constructor( nmlfile ) result(ANG)
!DEC$ ATTRIBUTES DLLEXPORT :: ANG_constructor
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! constructor for the ANG_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call ANG%readNameList(nmlfile)

end function ANG_constructor

!--------------------------------------------------------------------------
subroutine ANG_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: ANG_destructor
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! destructor for the ANG_T Class

IMPLICIT NONE

type(ANG_T), INTENT(INOUT)  :: self

call reportDestructor('ANG_T')

end subroutine ANG_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! read the namelist from an nml file for the ANG_T Class

use mod_io
use mod_EMsoft

IMPLICIT NONE

class(ANG_T), INTENT(INOUT) :: self
character(fnlen),INTENT(IN) :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN) :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)              :: EMsoft
type(IO_T)                  :: Message
logical                     :: skipread = .FALSE.

character(4)                :: modality
character(8)                :: angledataset   ! 'original' or 'refined'
logical                     :: isHROSM
character(fnlen)            :: xtalname
character(fnlen)            :: newangfile
character(fnlen)            :: dotproductfile

namelist /ANGlist/ modality, xtalname, angledataset, dotproductfile, newangfile, isHROSM

dotproductfile = 'undefined'
newangfile = 'undefined'
isHROSM = .FALSE.
xtalname = 'undefined'
modality = 'EBSD'
angledataset = 'original'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=ANGlist)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(dotproductfile).eq.'undefined') then
  call Message%printError('readNameList:',' dotproductfile name is undefined in '//nmlfile)
 end if

 if (trim(xtalname).eq.'undefined') then
  call Message%printError('readNameList:',' xtalname name is undefined in '//nmlfile)
 end if

 if (trim(newangfile).eq.'undefined') then
  call Message%printError('readNameList:',' newangfile name is undefined in '//nmlfile)
 end if
end if

self%nml%dotproductfile = dotproductfile
self%nml%newangfile = newangfile
self%nml%isHROSM = isHROSM
self%nml%xtalname = xtalname
self%nml%modality = trim(modality)
self%nml%angledataset = trim(angledataset)

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! pass the namelist for the ANG_T Class to the calling program

IMPLICIT NONE

class(ANG_T), INTENT(INOUT)          :: self
type(ANGNameListType)                :: nml

nml = self%nml

end function getNameList_


!--------------------------------------------------------------------------
function get_mod_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_mod_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! get modality from the ANG_T class

IMPLICIT NONE

class(ANG_T), INTENT(INOUT)     :: self
character(4)                    :: out

out = self%nml%modality

end function get_mod_

!--------------------------------------------------------------------------
subroutine set_mod_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_mod_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! set modality in the ANG_T class

IMPLICIT NONE

class(ANG_T), INTENT(INOUT)     :: self
character(4), INTENT(IN)        :: inp

self%nml%modality = inp

end subroutine set_mod_

!--------------------------------------------------------------------------
function get_isHROSM_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_isHROSM_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! get isHROSM from the ANG_T class

IMPLICIT NONE

class(ANG_T), INTENT(INOUT)     :: self
logical                         :: out

out = self%nml%isHROSM

end function get_isHROSM_

!--------------------------------------------------------------------------
subroutine set_isHROSM_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_isHROSM_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! set isHROSM in the ANG_T class

IMPLICIT NONE

class(ANG_T), INTENT(INOUT)     :: self
logical, INTENT(IN)             :: inp

self%nml%isHROSM = inp

end subroutine set_isHROSM_

!--------------------------------------------------------------------------
function get_angledataset_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_angledataset_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! get angledataset from the ANG_T class

IMPLICIT NONE

class(ANG_T), INTENT(INOUT)     :: self
character(8)                    :: out

out = self%nml%angledataset

end function get_angledataset_

!--------------------------------------------------------------------------
subroutine set_angledataset_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_angledataset_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! set angledataset in the ANG_T class

IMPLICIT NONE

class(ANG_T), INTENT(INOUT)     :: self
character(8), INTENT(IN)        :: inp

self%nml%angledataset = inp

end subroutine set_angledataset_

!--------------------------------------------------------------------------
function get_xtalname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_xtalname_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! get xtalname from the ANG_T class

IMPLICIT NONE

class(ANG_T), INTENT(INOUT)     :: self
character(fnlen)                :: out

out = self%nml%xtalname

end function get_xtalname_

!--------------------------------------------------------------------------
subroutine set_xtalname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_xtalname_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! set xtalname in the ANG_T class

IMPLICIT NONE

class(ANG_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)    :: inp

self%nml%xtalname = inp

end subroutine set_xtalname_

!--------------------------------------------------------------------------
function get_newangfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_newangfile_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! get newangfile from the ANG_T class

IMPLICIT NONE

class(ANG_T), INTENT(INOUT)     :: self
character(fnlen)                :: out

out = self%nml%newangfile

end function get_newangfile_

!--------------------------------------------------------------------------
subroutine set_newangfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_newangfile_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! set newangfile in the ANG_T class

IMPLICIT NONE

class(ANG_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)    :: inp

self%nml%newangfile = inp

end subroutine set_newangfile_

!--------------------------------------------------------------------------
function get_dotproductfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dotproductfile_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! get dotproductfile from the ANG_T class

IMPLICIT NONE

class(ANG_T), INTENT(INOUT)     :: self
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
!! set dotproductfile in the ANG_T class

IMPLICIT NONE

class(ANG_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)    :: inp

self%nml%dotproductfile = inp

end subroutine set_dotproductfile_

!--------------------------------------------------------------------------
subroutine ANG_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: ANG_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! perform the computations

use mod_EMsoft
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_vendors
use mod_DIfiles
use mod_crystallography
use mod_symmetry
use mod_io

IMPLICIT NONE

class(ANG_T), INTENT(INOUT)                      :: self
type(EMsoft_T), INTENT(INOUT)                    :: EMsoft
character(fnlen), INTENT(INOUT)                  :: progname

type(HDF_T)                                      :: HDF
type(HDFnames_T)                                 :: HDFnames
type(IO_T)                                       :: Message
type(DIfile_T)                                   :: DIFT
type(cell_T)                                     :: cell
type(SpaceGroup_T)                               :: SG
type(Vendor_T)                                   :: VT
type(DictionaryIndexingNameListType)             :: dinl

logical                                          :: stat, noindex, g_exists
character(fnlen)                                 :: dpfile, masterfile, energyfile
integer(kind=irg)                                :: hdferr, i, ii, jj, kk, iii, istat

real(kind=sgl),allocatable                       :: euler_best(:,:), CIlist(:), neweulers(:,:,:), newCI(:,:)
integer(kind=irg),allocatable                    :: indexmain(:,:)
real(kind=sgl),allocatable                       :: resultmain(:,:)
integer(HSIZE_T)                                 :: dims(1), sz(1), dims3(3), ddims2(2)

integer(kind=irg)                                :: numk, numdictsingle, numexptsingle

character(fnlen)                                 :: dataset, groupname
character(fnlen)                                 :: ename, fname

logical                                          :: verbose, readonly = .TRUE.
logical                                          :: f_exists, init=.TRUE., overwrite =.TRUE., refined
real(kind=sgl)                                   :: quat(4), ma, mi, dp, tstart, tstop, io_real(1), tmp, totnum_el, genfloat, &
                                                    vlen, fpar(1)
integer(kind=irg)                                :: ipar(10), Emin, Emax, nthreads, TID, io_int(2), tick, tock, ierr, L, pos
integer(kind=irg)                                :: ll, mm, jpar(7), Nexp, pgnum, FZcnt, dims2(2), ss(1), nlines, stlen 
real(kind=dbl)                                   :: prefactor, F
character(fnlen)                                 :: modality, DIfile, tmpnmlname, dotproductfile, targ
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)

VT = Vendor_T()
call VT%set_Modality(self%nml%modality)

associate(enl=>self%nml, DIDT=>DIFT%DIDT)

!====================================
! read the relevant fields from the dot product HDF5 file
dotproductfile = trim(enl%dotproductfile)

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()

HDFnames = HDFnames_T()

! if this is an EMHROSM output file, then we need to get the original dpfile
! parameter from its NMLparameters; then we make sure that file exists
! 
if (enl%isHROSM.eqv..TRUE.) then 
    call HDFnames%set_NMLfiles(SC_NMLfiles)
    call HDFnames%set_NMLfilename(SC_HROSMNML)
    call HDFnames%set_NMLparameters(SC_NMLparameters)
    call HDFnames%set_NMLlist(SC_HROSMNameListType)
    DIfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(enl%dotproductfile)
    DIFT = DIfile_T( fname = trim(DIfile) )
! is this a proper HDF5 file ?
    call h5fis_hdf5_f(trim(DIfile), stat, hdferr)
    if (stat.eqv..FALSE.) then ! the file exists, so let's open it an first make sure it is an EBSD dot product file
       call Message%printError('ANG_','This is not a proper HDF5 file: '//trim(DIfile))
    end if
    hdferr =  HDF%openFile(DIfile, readonly)
! make sure this is an HROSM dot product file
    hdferr = HDF%openGroup(HDFnames%get_NMLfiles())

    dataset = trim(HDFnames%get_NMLfilename())
      call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..FALSE.) then
       call Message%printError('ANG_','No HROSMNML data set found in '//trim(DIfile))
    end if
    call HDF%pop()
!====================================
! find the line with dpfile in the namelist field
!====================================
    hdferr = HDF%openGroup(HDFnames%get_NMLfiles())
    dataset = trim(HDFnames%get_NMLfilename())
    call HDF%readdatasetstringarray(dataset, nlines, hdferr, stringarray)
! close the HDF file
    call HDF%popall()
    sz = shape(stringarray)
    targ = 'dpfile'
! find the string that has dpfile in it...
    do i=1,sz(1)
      pos = index(stringarray(i), trim(targ))
      if (pos.ne.0) then
        pos = index(stringarray(i), " = ")
        stlen = len_trim(stringarray(i))
        targ = stringarray(i)(pos+4:stlen-2)
        call self%set_dotproductfile( targ )        
      end if 
    end do
end if 

call HDFnames%set_NMLfiles(SC_NMLfiles)
call HDFnames%set_NMLfilename(SC_DictionaryIndexingNML)
call HDFnames%set_NMLparameters(SC_NMLparameters)
call HDFnames%set_NMLlist(SC_DictionaryIndexingNameListType)

! if (trim(modalityname) .eq. 'EBSD') then
refined = .FALSE.
DIfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(enl%dotproductfile)
if (trim(enl%angledataset).eq.'refined') then
    call DIFT%readDotProductFile(EMsoft, HDF, HDFnames, DIfile, hdferr, &
                                 getCI=.TRUE., &
                                 getRefinedDotProducts=.TRUE., &
                                 getIQ=.TRUE., &
                                 getOSM=.TRUE., &
                                 getRefinedEulerAngles=.TRUE., &
                                 getPhi1=.TRUE., &
                                 getPhi=.TRUE., &
                                 getPhi2=.TRUE.)
    refined = .TRUE.
else
    call DIFT%readDotProductFile(EMsoft, HDF, HDFnames, DIfile, hdferr, &
                                 getCI=.TRUE., &
                                 getIQ=.TRUE., &
                                 getOSM=.TRUE., &
                                 getPhi1=.TRUE., &
                                 getPhi=.TRUE., &
                                 getPhi2=.TRUE.)
end if

dinl = DIFT%getNameList()

if (DIDT%Nexp.eq.-1) then
    ss = shape(DIDT%Phi1)
    Nexp = ss(1)
else
    Nexp = DIDT%Nexp
end if
allocate(euler_best(3,Nexp),CIlist(Nexp),stat=istat)
if (istat .ne. 0) then
    call Message%printError('ANG:',' Failed to allocate CIlist_new and/or euler_bestmatch array')
end if
euler_best = 0.0
CIlist = 0.0
! is this an EMHROSM output file ?
! if it is, then we need to replace the Euler angles by a new set of Euler angles...
if (enl%isHROSM.eqv..FALSE.) then 
    if (refined.eqv..FALSE.) then
        euler_best(1,1:Nexp) = DIDT%Phi1(1:Nexp)*rtod
        euler_best(2,1:Nexp) = DIDT%Phi(1:Nexp)*rtod
        euler_best(3,1:Nexp) = DIDT%Phi2(1:Nexp)*rtod
        if (allocated(DIDT%Phi1)) deallocate(DIDT%Phi1)
        if (allocated(DIDT%Phi)) deallocate(DIDT%Phi)
        if (allocated(DIDT%Phi2)) deallocate(DIDT%Phi2)
        call Message%printMessage(' Using original Euler angles from dot product/SI file')
        if (allocated(DIDT%CI)) then 
          CIlist(1:Nexp) = DIDT%CI(1:Nexp)
          deallocate(DIDT%CI)
      end if
    else
        euler_best(1,1:Nexp) = DIDT%RefinedEulerAngles(1,1:Nexp)*rtod
        euler_best(2,1:Nexp) = DIDT%RefinedEulerAngles(2,1:Nexp)*rtod
        euler_best(3,1:Nexp) = DIDT%RefinedEulerAngles(3,1:Nexp)*rtod
        if (allocated(DIDT%RefinedEulerAngles)) deallocate(DIDT%RefinedEulerAngles)
        if (allocated(DIDT%RefinedDotProducts)) then 
          CIlist(1:Nexp) = DIDT%RefinedDotProducts(1:Nexp)
          deallocate(DIDT%RefinedDotProducts)
        end if 
        call Message%printMessage(' Using refined Euler angles from dot product/SI file')
    end if
else
! directly read the newEuler array from the EMHROSM output file and generate the proper arrays
    call HDFnames%set_NMLfiles(SC_NMLfiles)
    call HDFnames%set_NMLfilename(SC_HROSMNML)
    call HDFnames%set_NMLparameters(SC_NMLparameters)
    call HDFnames%set_NMLlist(SC_HROSMNameListType)
    DIfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(dotproductfile)
    hdferr = HDF%openFile(DIfile, readonly)
    groupname = 'EMData'
    hdferr = HDF%openGroup(groupname)
    groupname = 'HROSM'
    hdferr = HDF%openGroup(groupname)
    dataset = 'newEuler'
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetFloatArray(dataset, dims3, hdferr, neweulers)
    else
      call Message%printError('ANG_ ', '  --> no newEuler data set found')
    end if
    euler_best = reshape( neweulers, (/ 3, Nexp /) ) * rtod
    deallocate(neweulers)
!
    dataset = 'newCI'
    call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetFloatArray(dataset, ddims2, hdferr, newCI)
    else
      call Message%printError('ANG_ ', '  --> no newCI data set found')
    end if
    CIlist = reshape( newCI, (/ Nexp /) )
    deallocate(newCI)
    call HDF%popall()
end if 

call cell%setFileName(enl%xtalname)
call cell%readDataHDF(SG, EMsoft, HDF)
pgnum = SG%getPGnumber()

! and prepare the .ang output file
dinl%angfile = trim(enl%newangfile)

ipar = 0
ipar(1) = 1
ipar(2) = Nexp
ipar(3) = Nexp
ipar(4) = Nexp
ipar(5) = FZcnt
ipar(6) = pgnum
if (sum(dinl%ROI).ne.0) then
    ipar(7) = dinl%ROI(3)
    ipar(8) = dinl%ROI(4)
else
    ipar(7) = dinl%ipf_wd
    ipar(8) = dinl%ipf_ht
end if

fpar(1) = 10.0

allocate(indexmain(ipar(1),1:ipar(2)),resultmain(ipar(1),1:ipar(2)))
indexmain = 0
resultmain(1,1:ipar(2)) = CIlist(1:Nexp)

if (dinl%angfile.ne.'undefined') then
  call VT%ang_writeFile(EMsoft,cell,SG,dinl,ipar,fpar,indexmain,euler_best,resultmain,DIDT%IQ,noindex=.TRUE.)
  call Message%printMessage('Data stored in ang file : '//trim(enl%newangfile))
end if

! close the fortran HDF interface
call closeFortranHDFInterface()

end associate

end subroutine ANG_

end module mod_ANG
