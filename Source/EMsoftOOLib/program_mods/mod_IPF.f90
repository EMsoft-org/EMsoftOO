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

module mod_IPF
  !! author: MDG 
  !! version: 1.0 
  !! date: 09/03/21
  !!
  !! class definition for the EMgetIPF program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMIPF program
type, public :: IPFNameListType
  character(fnlen)        :: dotproductfile
  character(fnlen)        :: IPFfilename
  character(fnlen)        :: IPFmode
  integer(kind=irg)       :: sampleDir(3)
  integer(kind=irg)       :: nthreads
  integer(kind=irg)       :: pgnum
end type IPFNameListType

! class definition
type, public :: IPF_T
private 
  character(fnlen)       :: nmldeffile = 'EMgetIPF.nml'
  type(IPFNameListType)  :: nml 

contains
private 
  procedure, pass(self) :: get_dotproductfile_
  procedure, pass(self) :: get_IPFfilename_
  procedure, pass(self) :: get_sampleDir_
  procedure, pass(self) :: get_nthreads_
  procedure, pass(self) :: get_IPFmode_
  procedure, pass(self) :: get_pgnum_
  procedure, pass(self) :: set_dotproductfile_
  procedure, pass(self) :: set_IPFfilename_
  procedure, pass(self) :: set_sampleDir_
  procedure, pass(self) :: set_nthreads_
  procedure, pass(self) :: set_IPFmode_
  procedure, pass(self) :: set_pgnum_
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: getIPF_
  procedure, pass(self) :: updateIPFmap_

  generic, public :: get_dotproductfile => get_dotproductfile_
  generic, public :: get_IPFfilename => get_IPFfilename_
  generic, public :: get_sampleDir => get_sampleDir_
  generic, public :: get_nthreads => get_nthreads_
  generic, public :: get_pgnum => get_pgnum_
  generic, public :: get_IPFmode => get_IPFmode_
  generic, public :: set_dotproductfile => set_dotproductfile_
  generic, public :: set_IPFfilename => set_IPFfilename_
  generic, public :: set_sampleDir => set_sampleDir_
  generic, public :: set_pgnum => set_pgnum_
  generic, public :: set_nthreads => set_nthreads_
  generic, public :: set_IPFmode => set_IPFmode_
  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: getIPF => getIPF_
  generic, public :: updateIPFmap => updateIPFmap_

end type IPF_T

! the constructor routine for this class 
interface IPF_T
  module procedure IPF_constructor
end interface IPF_T

contains

!--------------------------------------------------------------------------
type(IPF_T) function IPF_constructor( nmlfile ) result(IPF)
!DEC$ ATTRIBUTES DLLEXPORT :: IPF_constructor
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! constructor for the IPF_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

! the calling program must assign values to all entries if no nmlfile parameter is present
if (present(nmlfile)) then ! read the namelist file
  call IPF%readNameList(nmlfile)
end if 

end function IPF_constructor

!--------------------------------------------------------------------------
subroutine IPF_destructor(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: IPF_destructor
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! destructor for the IPF_T Class
 
IMPLICIT NONE

type(IPF_T), INTENT(INOUT)  :: self 

call reportDestructor('IPF_T')

end subroutine IPF_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! read the namelist from an nml file for the IPF_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

character(fnlen)        :: dotproductfile
character(fnlen)        :: IPFfilename
character(fnlen)        :: IPFmode
integer(kind=irg)       :: sampleDir(3)
integer(kind=irg)       :: nthreads
integer(kind=irg)       :: pgnum

! define the IO namelist to facilitate passing variables to the program.
namelist  / getIPF / dotproductfile, IPFfilename, sampleDir, nthreads, IPFmode, pgnum

dotproductfile = 'undefined'
sampleDir = (/ 0, 0, 1 /)
IPFfilename = 'undefined'
IPFmode = 'TSL' 
nthreads = 1
pgnum = 32

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=getIPF)
    close(UNIT=dataunit,STATUS='keep')

    if (trim(dotproductfile).eq.'undefined') then
        call Message%printError('readNameList:',' dot product file name is undefined in '//nmlfile)
    end if

    if (trim(IPFfilename).eq.'undefined') then
        call Message%printError('readNameList:',' IPF map file name is undefined in '//nmlfile)
    end if
end if 

self%nml%dotproductfile = dotproductfile
self%nml%IPFfilename = IPFfilename
self%nml%IPFmode = IPFmode
self%nml%sampleDir = sampleDir
self%nml%nthreads = nthreads
self%nml%pgnum = pgnum

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! pass the namelist for the IPF_T Class to the calling program

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)          :: self
type(IPFNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
function get_dotproductfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dotproductfile_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get dotproductfile from the IPF_T class

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)     :: self
character(fnlen)                :: out

out = self%nml%dotproductfile

end function get_dotproductfile_

!--------------------------------------------------------------------------
subroutine set_dotproductfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_dotproductfile_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set dotproductfile in the IPF_T class

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)    :: inp

self%nml%dotproductfile = inp

end subroutine set_dotproductfile_

!--------------------------------------------------------------------------
function get_IPFfilename_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_IPFfilename_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get IPFfilename from the IPF_T class

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)     :: self
character(fnlen)                :: out

out = self%nml%IPFfilename

end function get_IPFfilename_

!--------------------------------------------------------------------------
subroutine set_IPFfilename_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_IPFfilename_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set IPFfilename in the IPF_T class

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)    :: inp

self%nml%IPFfilename = inp

end subroutine set_IPFfilename_

!--------------------------------------------------------------------------
function get_IPFmode_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_IPFmode_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get IPFmode from the IPF_T class

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)     :: self
character(fnlen)                :: out

out = self%nml%IPFmode

end function get_IPFmode_

!--------------------------------------------------------------------------
subroutine set_IPFmode_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_IPFmode_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set IPFmode in the IPF_T class

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)    :: inp

self%nml%IPFmode = inp

end subroutine set_IPFmode_

!--------------------------------------------------------------------------
function get_sampleDir_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_sampleDir_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get sampleDir from the IPF_T class

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out(3)

out = self%nml%sampleDir

end function get_sampleDir_

!--------------------------------------------------------------------------
subroutine set_sampleDir_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_sampleDir_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set sampleDir in the IPF_T class

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp(3)

self%nml%sampleDir = inp

end subroutine set_sampleDir_

!--------------------------------------------------------------------------
function get_nthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nthreads_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get nthreads from the IPF_T class

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out

out = self%nml%nthreads

end function get_nthreads_

!--------------------------------------------------------------------------
subroutine set_nthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nthreads_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set nthreads in the IPF_T class

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp

self%nml%nthreads = inp

end subroutine set_nthreads_

!--------------------------------------------------------------------------
function get_pgnum_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_pgnum_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get pgnum from the IPF_T class

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out

out = self%nml%pgnum

end function get_pgnum_

!--------------------------------------------------------------------------
subroutine set_pgnum_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_pgnum_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set pgnum in the IPF_T class

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp

self%nml%pgnum = inp

end subroutine set_pgnum_

!--------------------------------------------------------------------------
subroutine getIPF_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: getIPF_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! perform the computations

use mod_EMsoft
use mod_IPFsupport
use mod_rotations
use mod_quaternions
use mod_DIfiles
use mod_MPfiles
use mod_symmetry
use HDF5
use mod_HDFsupport
use mod_rotations
use mod_HDFnames
use ISO_C_BINDING
use mod_io

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)             :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 

type(IPFmap_T)                          :: IPFmap 
type(HDF_T)                             :: HDF
type(HDFnames_T)                        :: HDFnames
type(DIfile_T)                          :: DIFT
type(QuaternionArray_T)                 :: qAR, sym
type(e_T)                               :: e 
type(q_T)                               :: q 
type(Quaternion_T)                      :: qu
type(IO_T)                              :: Message

character(fnlen)                        :: nmldeffile, progdesc, DIfile, fname
integer(kind=irg)                       :: hdferr, i, LaueClass, Pm, io_int(2)

associate(csnl=>self%nml, dinl=>DIFT%nml, DIDT=>DIFT%DIDT)

call setRotationPrecision('d')

!====================================
! read the relevant fields from the dot product HDF5 file

! open the fortran HDF interface
call openFortranHDFInterface()
HDF = HDF_T()

HDFnames = HDFnames_T()

call HDFnames%set_NMLfiles(SC_NMLfiles)
call HDFnames%set_NMLfilename(SC_DictionaryIndexingNML)
call HDFnames%set_NMLparameters(SC_NMLparameters)
call HDFnames%set_NMLlist(SC_DictionaryIndexingNameListType)

DIfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(csnl%dotproductfile)

DIDT%Nexp = dinl%ipf_wd*dinl%ipf_ht

call DIFT%readDotProductFile(EMsoft, HDF, HDFnames, DIfile, hdferr, &
                             getRefinedEulerAngles=.TRUE., &
                             getEulerAngles=.TRUE. )

qAR = QuaternionArray_T( DIDT%Nexp, s = 'd')

if (allocated(DIDT%RefinedEulerAngles)) then 
  do i=1,DIDT%Nexp 
   e = e_T( edinp = dble(DIDT%RefinedEulerAngles(1:3,i)) )
   q = e%eq()
   qu = quaternion_T( qd = q%q_copyd() )
   call qAR%insertQuatinArray( i, qu )
  end do 
  deallocate(DIDT%RefinedEulerAngles)
else 
  do i=1,DIDT%Nexp 
   e = e_T( edinp = dble(DIDT%EulerAngles(1:3,i)) )
   q = e%eq()
   qu = quaternion_T( qd = q%q_copyd() )
   call qAR%insertQuatinArray( i, qu )
  end do 
  deallocate(DIDT%EulerAngles)
end if 

! get the symmetry operator quaternions for the point group
io_int = (/ csnl%pgnum, PGLaueinv(csnl%pgnum) /)
call Message%WriteValue(' initializing symmetry operators for point/Laue group ', io_int, 2)

call qAR%QSym_Init(csnl%pgnum, sym)

! create the IPFmap class, set the parameters, and generate the IPF map 
IPFmap = IPFmap_T()

call IPFmap%set_ipf_LaueClass(PGLaueinv(csnl%pgnum))
call IPFmap%set_ipf_wd(dinl%ipf_wd)
call IPFmap%set_ipf_ht(dinl%ipf_ht)
call IPFmap%set_ipf_mode(csnl%IPFmode)
call IPFmap%set_ipf_filename(csnl%IPFfilename)
call IPFmap%set_ipf_nthreads(csnl%nthreads)

call IPFmap%get_IPFMap(EMsoft, csnl%sampleDir, qAR, sym)

end associate

end subroutine getIPF_

!--------------------------------------------------------------------------
subroutine updateIPFmap_(self, EMsoft, progname, ipf_wd, ipf_ht, pgnum, IPFname, qAR, sym)
!DEC$ ATTRIBUTES DLLEXPORT :: updateIPFmap_
!! author: MDG 
!! version: 1.0 
!! date: 01/23/22
!!
!! updates an IPF map file during dictionary indexing ...

use mod_EMsoft
use mod_IPFsupport
use mod_quaternions
use mod_symmetry
use ISO_C_BINDING

IMPLICIT NONE 

class(IPF_T), INTENT(INOUT)             :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 
integer(kind=irg), INTENT(IN)           :: ipf_wd
integer(kind=irg), INTENT(IN)           :: ipf_ht
integer(kind=irg), INTENT(IN)           :: pgnum 
character(fnlen),INTENT(INOUT)          :: IPFname 
type(QuaternionArray_T),INTENT(INOUT)   :: qAR
type(QuaternionArray_T),INTENT(INOUT)   :: sym

type(IPFmap_T)                          :: IPFmap 
character(fnlen)                        :: IPFmode

!====================================
! this is called from inside the dictionary indexing routine, so 
! we do not read any HDF5 file; instead, we use the arrays that 
! are passed into this routine

! create the IPFmap class, set the parameters, and generate the IPF-Z map 
IPFmap = IPFmap_T()

call IPFmap%set_ipf_LaueClass(PGLaueinv(pgnum))
call IPFmap%set_ipf_wd(ipf_wd)
call IPFmap%set_ipf_ht(ipf_ht)
IPFmode = 'TSL'
call IPFmap%set_ipf_mode(IPFmode)
call IPFmap%set_ipf_filename(IPFname)
call IPFmap%set_ipf_nthreads(1)

call IPFmap%get_IPFMap(EMsoft, (/ 0, 0, 1/), qAR, sym, cDir=.TRUE.)

end subroutine updateIPFmap_

end module mod_IPF
