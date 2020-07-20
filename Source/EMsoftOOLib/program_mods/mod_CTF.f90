! ###################################################################
! Copyright (c) 2013-2020, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_CTF
  !! author: MDG
  !! version: 1.0
  !! date: 04/06/20
  !!
  !! class definition for the EMgetCTF program

use mod_kinds
use mod_global

IMPLICIT NONE

! namelist for the EMgetCTF program
type, public :: CTFNameListType
  character(4)      :: modality
  character(8)      :: angledataset   ! 'original' or 'refined'
  character(fnlen)  :: xtalname
  character(fnlen)  :: newctffile
  character(fnlen)  :: dotproductfile
end type CTFNameListType

! class definition
type, public :: CTF_T
private
  character(fnlen)       :: nmldeffile = 'EMgetCTF.nml'
  type(CTFNameListType)  :: nml

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: CTF_
  procedure, pass(self) :: get_mod_
  procedure, pass(self) :: get_angledataset_
  procedure, pass(self) :: get_xtalname_
  procedure, pass(self) :: get_newctffile_
  procedure, pass(self) :: get_dotproductfile_
  procedure, pass(self) :: set_mod_
  procedure, pass(self) :: set_angledataset_
  procedure, pass(self) :: set_xtalname_
  procedure, pass(self) :: set_newctffile_
  procedure, pass(self) :: set_dotproductfile_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: CTF => CTF_
  generic, public :: get_mod => get_mod_
  generic, public :: get_angledataset => get_angledataset_
  generic, public :: get_xtalname => get_xtalname_
  generic, public :: get_newctffile => get_newctffile_
  generic, public :: get_dotproductfile => get_dotproductfile_
  generic, public :: set_mod => set_mod_
  generic, public :: set_angledataset => set_angledataset_
  generic, public :: set_xtalname => set_xtalname_
  generic, public :: set_newctffile => set_newctffile_
  generic, public :: set_dotproductfile => set_dotproductfile_
end type CTF_T

! the constructor routine for this class
interface CTF_T
  module procedure CTF_constructor
end interface CTF_T

contains

!--------------------------------------------------------------------------
type(CTF_T) function CTF_constructor( nmlfile ) result(CTF)
!DEC$ ATTRIBUTES DLLEXPORT :: CTF_constructor
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! constructor for the CTF_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call CTF%readNameList(nmlfile)

end function CTF_constructor

!--------------------------------------------------------------------------
subroutine CTF_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: CTF_destructor
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! destructor for the CTF_T Class

IMPLICIT NONE

type(CTF_T), INTENT(INOUT)  :: self

call reportDestructor('CTF_T')

end subroutine CTF_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! read the namelist from an nml file for the CTF_T Class

use mod_io
use mod_EMsoft

IMPLICIT NONE

class(CTF_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft
type(IO_T)                           :: Message
logical                              :: skipread = .FALSE.

character(4)            :: modality
character(8)            :: angledataset   ! 'original' or 'refined'
character(fnlen)        :: xtalname
character(fnlen)        :: newctffile
character(fnlen)        :: dotproductfile

namelist /CTFlist/ modality, xtalname, angledataset, dotproductfile, newctffile

dotproductfile = 'undefined'
newctffile = 'undefined'
xtalname = 'undefined'
modality = 'EBSD'
angledataset = 'original'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=CTFlist)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(dotproductfile).eq.'undefined') then
  call Message%printError('readNameList:',' dotproductfile name is undefined in '//nmlfile)
 end if

 if (trim(xtalname).eq.'undefined') then
  call Message%printError('readNameList:',' xtalname name is undefined in '//nmlfile)
 end if

 if (trim(newctffile).eq.'undefined') then
  call Message%printError('readNameList:',' newctffile name is undefined in '//nmlfile)
 end if
end if

self%nml%dotproductfile = dotproductfile
self%nml%newctffile = newctffile
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
!! pass the namelist for the CTF_T Class to the calling program

IMPLICIT NONE

class(CTF_T), INTENT(INOUT)          :: self
type(CTFNameListType)                :: nml

nml = self%nml

end function getNameList_


!--------------------------------------------------------------------------
function get_mod_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_mod_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! get modality from the CTF_T class

IMPLICIT NONE

class(CTF_T), INTENT(INOUT)     :: self
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
!! set modality in the CTF_T class

IMPLICIT NONE

class(CTF_T), INTENT(INOUT)     :: self
character(4), INTENT(IN)        :: inp

self%nml%modality = inp

end subroutine set_mod_

!--------------------------------------------------------------------------
function get_angledataset_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_angledataset_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! get angledataset from the CTF_T class

IMPLICIT NONE

class(CTF_T), INTENT(INOUT)     :: self
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
!! set angledataset in the CTF_T class

IMPLICIT NONE

class(CTF_T), INTENT(INOUT)     :: self
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
!! get xtalname from the CTF_T class

IMPLICIT NONE

class(CTF_T), INTENT(INOUT)     :: self
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
!! set xtalname in the CTF_T class

IMPLICIT NONE

class(CTF_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)    :: inp

self%nml%xtalname = inp

end subroutine set_xtalname_

!--------------------------------------------------------------------------
function get_newctffile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_newctffile_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! get newctffile from the CTF_T class

IMPLICIT NONE

class(CTF_T), INTENT(INOUT)     :: self
character(fnlen)                :: out

out = self%nml%newctffile

end function get_newctffile_

!--------------------------------------------------------------------------
subroutine set_newctffile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_newctffile_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! set newctffile in the CTF_T class

IMPLICIT NONE

class(CTF_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)    :: inp

self%nml%newctffile = inp

end subroutine set_newctffile_

!--------------------------------------------------------------------------
function get_dotproductfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dotproductfile_
!! author: MDG
!! version: 1.0
!! date: 04/06/20
!!
!! get dotproductfile from the CTF_T class

IMPLICIT NONE

class(CTF_T), INTENT(INOUT)     :: self
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
!! set dotproductfile in the CTF_T class

IMPLICIT NONE

class(CTF_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)    :: inp

self%nml%dotproductfile = inp

end subroutine set_dotproductfile_

!--------------------------------------------------------------------------
subroutine CTF_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: CTF_
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

class(CTF_T), INTENT(INOUT)             :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname

type(HDF_T)                             :: HDF
type(HDFnames_T)                        :: HDFnames
type(IO_T)                              :: Message
type(DIfile_T)                          :: DIFT
type(cell_T)                            :: cell
type(SpaceGroup_T)                      :: SG
type(Vendor_T)                          :: VT
type(DictionaryIndexingNameListType)    :: dinl

logical                                 :: stat, readonly, noindex, g_exists
character(fnlen)                        :: dpfile, masterfile, energyfile
integer(kind=irg)                       :: hdferr, ii, jj, kk, iii, istat

real(kind=sgl),allocatable              :: euler_best(:,:), CIlist(:)
integer(kind=irg),allocatable           :: indexmain(:,:)
real(kind=sgl),allocatable              :: resultmain(:,:)
integer(HSIZE_T)                        :: dims(1)

integer(kind=irg)                       :: numk, numdictsingle, numexptsingle

character(fnlen, KIND=c_char),allocatable,TARGET    :: stringarray(:)
character(fnlen)                        :: dataset, groupname
character(fnlen)                        :: ename, fname

logical                                 :: verbose
logical                                 :: f_exists, init=.TRUE., overwrite =.TRUE., refined
real(kind=sgl)                          :: quat(4), ma, mi, dp, tstart, tstop, io_real(1), tmp, totnum_el, genfloat, &
                                           vlen, fpar(2)
integer(kind=irg)                       :: ipar(10), Emin, Emax, nthreads, TID, io_int(2), tick, tock, ierr, L
integer(kind=irg)                       :: ll, mm, jpar(7), Nexp, pgnum, FZcnt, nlines, dims2(2), ss(1)
real(kind=dbl)                          :: prefactor, F
character(fnlen)                        :: modality, DIfile

VT = Vendor_T()
call VT%set_Modality(self%nml%modality)

associate(enl=>self%nml, DIDT=>DIFT%DIDT)

!====================================
! read the relevant fields from the dot product HDF5 file

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()

HDFnames = HDFnames_T()

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
    call Message%printError('CTF:',' Failed to allocate CIlist_new and/or euler_bestmatch array')
end if
euler_best = 0.0
CIlist = 0.0
if (refined.eqv..FALSE.) then
    euler_best(1,1:Nexp) = DIDT%Phi1(1:Nexp)*rtod
    euler_best(2,1:Nexp) = DIDT%Phi(1:Nexp)*rtod
    euler_best(3,1:Nexp) = DIDT%Phi2(1:Nexp)*rtod
    if (allocated(DIDT%Phi1)) deallocate(DIDT%Phi1)
    if (allocated(DIDT%Phi)) deallocate(DIDT%Phi)
    if (allocated(DIDT%Phi2)) deallocate(DIDT%Phi2)
    call Message%printMessage(' Using original Euler angles from dot product/SI file')
else
    euler_best(1,1:Nexp) = DIDT%RefinedEulerAngles(1,1:Nexp)*rtod
    euler_best(2,1:Nexp) = DIDT%RefinedEulerAngles(2,1:Nexp)*rtod
    euler_best(3,1:Nexp) = DIDT%RefinedEulerAngles(3,1:Nexp)*rtod
    if (allocated(DIDT%RefinedEulerAngles)) deallocate(DIDT%RefinedEulerAngles)
    call Message%printMessage(' Using refined Euler angles from dot product/SI file')
end if
if (allocated(DIDT%CI)) then
    CIlist(1:Nexp) = DIDT%CI(1:Nexp)
    deallocate(DIDT%CI)
else
    CIlist = 0.0
end if

call cell%setFileName(enl%xtalname)
call cell%readDataHDF(SG, EMsoft, HDF)
pgnum = SG%getPGnumber()

! and prepare the .ctf output file
dinl%ctffile = trim(enl%newctffile)

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

fpar(1) = dinl%energymax
fpar(2) = DIDT%MCsig

allocate(indexmain(ipar(1),1:ipar(2)),resultmain(ipar(1),1:ipar(2)))
indexmain = 0
resultmain(1,1:ipar(2)) = CIlist(1:Nexp)

if (dinl%ctffile.ne.'undefined') then
  call VT%ctf_writeFile(EMsoft,cell,SG,dinl,ipar,fpar,indexmain,euler_best,resultmain,DIDT%OSM,DIDT%IQ,noindex=.TRUE.)
  call Message%printMessage('Data stored in ctf file : '//trim(enl%newctffile))
end if

! close the fortran HDF interface
call closeFortranHDFInterface()

end associate

end subroutine CTF_

end module mod_CTF
