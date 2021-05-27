! ###################################################################
! Copyright (c) 2013-2021, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_DIsetting
  !! author: MDG
  !! version: 1.0
  !! date: 04/07/20
  !!
  !! class definition for the EMDIsetting program

use mod_kinds
use mod_global

IMPLICIT NONE

! namelist for the EMDIsetting program
type, public :: DIsettingNameListType
  integer(kind=irg) :: nthreads
  integer(kind=irg) :: orthorhombicSetting
  character(fnlen)  :: dotproductfile
  character(fnlen)  :: newctffile
end type DIsettingNameListType

! class definition
type, public :: DIsetting_T
private
  character(fnlen)       :: nmldeffile = 'EMDIsetting.nml'
  type(DIsettingNameListType)  :: nml

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: DIsetting_
  procedure, pass(self) :: get_nthreads_
  procedure, pass(self) :: get_orthorhombicSetting_
  procedure, pass(self) :: get_dotproductfile_
  procedure, pass(self) :: get_newctffile_
  procedure, pass(self) :: set_nthreads_
  procedure, pass(self) :: set_orthorhombicSetting_
  procedure, pass(self) :: set_dotproductfile_
  procedure, pass(self) :: set_newctffile_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: DIsetting => DIsetting_
  generic, public :: get_nthreads => get_nthreads_
  generic, public :: get_orthorhombicSetting => get_orthorhombicSetting_
  generic, public :: get_dotproductfile => get_dotproductfile_
  generic, public :: get_newctffile => get_newctffile_
  generic, public :: set_nthreads => set_nthreads_
  generic, public :: set_orthorhombicSetting => set_orthorhombicSetting_
  generic, public :: set_dotproductfile => set_dotproductfile_
  generic, public :: set_newctffile => set_newctffile_

end type DIsetting_T

! the constructor routine for this class
interface DIsetting_T
  module procedure DIsetting_constructor
end interface DIsetting_T

contains

!--------------------------------------------------------------------------
type(DIsetting_T) function DIsetting_constructor( nmlfile ) result(DIsetting)
!DEC$ ATTRIBUTES DLLEXPORT :: DIsetting_constructor
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! constructor for the DIsetting_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call DIsetting%readNameList(nmlfile)

end function DIsetting_constructor

!--------------------------------------------------------------------------
subroutine DIsetting_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: DIsetting_destructor
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! destructor for the DIsetting_T Class

IMPLICIT NONE

type(DIsetting_T), INTENT(INOUT)  :: self

call reportDestructor('DIsetting_T')

end subroutine DIsetting_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! read the namelist from an nml file for the DIsetting_T Class

use mod_io
use mod_EMsoft

IMPLICIT NONE

class(DIsetting_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft
type(IO_T)                           :: Message
logical                              :: skipread = .FALSE.

integer(kind=irg)  :: orthorhombicSetting
integer(kind=irg)  :: nthreads
character(fnlen)   :: dotproductfile
character(fnlen)   :: newctffile

namelist /ChangeSettingslist/ nthreads, orthorhombicSetting, dotproductfile, newctffile

dotproductfile = 'undefined'
newctffile = 'undefined'
orthorhombicSetting = 1
nthreads = 1

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=ChangeSettingslist)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(dotproductfile).eq.'undefined') then
  call Message%printError('readNameList:',' dotproductfile name is undefined in '//nmlfile)
 end if

 if (trim(newctffile).eq.'undefined') then
  call Message%printError('readNameList:',' newctffile name is undefined in '//nmlfile)
 end if
end if

self%nml%dotproductfile = dotproductfile
self%nml%newctffile = newctffile
self%nml%nthreads = nthreads
self%nml%orthorhombicSetting = orthorhombicSetting

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! pass the namelist for the DIsetting_T Class to the calling program

IMPLICIT NONE

class(DIsetting_T), INTENT(INOUT)          :: self
type(DIsettingNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants

use ISO_C_BINDING

IMPLICIT NONE

class(DIsetting_T), INTENT(INOUT)        :: self
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 11, n_real = 9
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( mcnl => self%nml )

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
function get_nthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nthreads_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! get nthreads from the DIsetting_T class

IMPLICIT NONE

class(DIsetting_T), INTENT(INOUT)     :: self
integer(kind=irg)                     :: out

out = self%nml%nthreads

end function get_nthreads_

!--------------------------------------------------------------------------
subroutine set_nthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nthreads_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! set nthreads in the DIsetting_T class

IMPLICIT NONE

class(DIsetting_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%nthreads = inp

end subroutine set_nthreads_

!--------------------------------------------------------------------------
function get_orthorhombicSetting_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_orthorhombicSetting_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! get orthorhombicSetting from the DIsetting_T class

IMPLICIT NONE

class(DIsetting_T), INTENT(INOUT)     :: self
integer(kind=irg)                     :: out

out = self%nml%orthorhombicSetting

end function get_orthorhombicSetting_

!--------------------------------------------------------------------------
subroutine set_orthorhombicSetting_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_orthorhombicSetting_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! set orthorhombicSetting in the DIsetting_T class

IMPLICIT NONE

class(DIsetting_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)         :: inp

self%nml%orthorhombicSetting = inp

end subroutine set_orthorhombicSetting_

!--------------------------------------------------------------------------
function get_dotproductfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dotproductfile_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! get dotproductfile from the DIsetting_T class

IMPLICIT NONE

class(DIsetting_T), INTENT(INOUT)     :: self
character(fnlen)                      :: out

out = self%nml%dotproductfile

end function get_dotproductfile_

!--------------------------------------------------------------------------
subroutine set_dotproductfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_dotproductfile_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! set dotproductfile in the DIsetting_T class

IMPLICIT NONE

class(DIsetting_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%dotproductfile = inp

end subroutine set_dotproductfile_

!--------------------------------------------------------------------------
function get_newctffile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_newctffile_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! get newctffile from the DIsetting_T class

IMPLICIT NONE

class(DIsetting_T), INTENT(INOUT)     :: self
character(fnlen)                      :: out

out = self%nml%newctffile

end function get_newctffile_

!--------------------------------------------------------------------------
subroutine set_newctffile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_newctffile_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! set newctffile in the DIsetting_T class

IMPLICIT NONE

class(DIsetting_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)          :: inp

self%nml%newctffile = inp

end subroutine set_newctffile_

!--------------------------------------------------------------------------
subroutine DIsetting_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: DIsetting_
!! author: MDG
!! version: 1.0
!! date: 04/07/20
!!
!! perform the computations

use mod_EMsoft
use mod_io
use mod_crystallography
use mod_symmetry
use mod_rotations
use mod_quaternions
use mod_DIfiles
use mod_MPfiles
use mod_vendors
use HDF5
use mod_HDFsupport
use mod_rotations
use omp_lib
use mod_OMPsupport
use mod_HDFnames
use ISO_C_BINDING

IMPLICIT NONE

class(DIsetting_T), INTENT(INOUT)           :: self
type(EMsoft_T), INTENT(INOUT)               :: EMsoft
character(fnlen), INTENT(INOUT)             :: progname

type(IO_T)                                  :: Message
type(HDF_T)                                 :: HDF
type(HDFnames_T)                            :: HDFnames
type(Vendor_T)                              :: VT
type(DIfile_T)                              :: DIFT
type(MPfile_T)                              :: MPFT
type(cell_T)                                :: cell
type(SpaceGroup_T)                          :: SG
type(e_T)                                   :: eu
type(q_T)                                   :: qrot, qin
type(Quaternion_T)                          :: qr, qrin, qnew

type(EBSDmasterNameListType)                :: mpnl

character(fnlen)                            :: nmldeffile, progdesc, DIfile, fname
integer(kind=irg)                           :: hdferr, pgnum, sgnum, orthonum, i, TID, ipar(10)
character(fnlen)                            :: outstring, dataset, infile, groupname, comment, modality
real(kind=dbl),allocatable                  :: newEulers(:,:), newAvOr(:,:), oldEulers(:,:), newRefined(:,:)
real(kind=sgl),allocatable                  :: eulers(:,:), ang(:), resultmain(:,:)
real(kind=sgl)                              :: fpar2(2)
logical                                     :: g_exists, readonly, verbose, overwrite = .TRUE., transformRefined
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)

associate(csnl=>self%nml, dinl=>DIFT%nml, DIDT=>DIFT%DIDT, MPDT=>MPFT%MPDT)

call setRotationPrecision('d')

!====================================
! 1. read the relevant fields from the dot product HDF5 file
! open the fortran HDF interface
! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()

HDFnames = HDFnames_T()

call HDFnames%set_NMLfiles(SC_NMLfiles)
call HDFnames%set_NMLfilename(SC_DictionaryIndexingNML)
call HDFnames%set_NMLparameters(SC_NMLparameters)
call HDFnames%set_NMLlist(SC_DictionaryIndexingNameListType)

DIfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(csnl%dotproductfile)
call DIFT%readDotProductFile(EMsoft, HDF, HDFnames, DIfile,  hdferr, &
                             getCI=.TRUE., &
                             getIQ=.TRUE., &
                             getOSM=.TRUE., &
                             getPhi1=.TRUE., &
                             getPhi=.TRUE., &
                             getPhi2=.TRUE., &
                             getRefinedEulerAngles=.TRUE., &
                             getAverageOrientations=.TRUE., &
                             getEulerAngles=.TRUE., &
                             getTopMatchIndices=.TRUE.)

allocate(oldEulers(3,DIDT%FZcnt))
allocate(newEulers(3,DIDT%FZcnt))
allocate(newAvOr(3,DIDT%Nexp))

transformRefined = .FALSE.
if (allocated(DIDT%RefinedEulerAngles)) then
  transformRefined = .TRUE.
  allocate(newRefined(3,DIDT%Nexp))
end if

do i=1,DIDT%FZcnt
  oldEulers(1:3,i) = DIDT%EulerAngles(1:3,i)
end do
newEulers = oldEulers * dtor
newAvOr = DIDT%AverageOrientations * dtor

call Message%printMessage(' ')
call Message%printMessage(' -> completed reading of dot product file')

!====================================
! 2. read EBSD master pattern file (including HDF format)

call HDFnames%set_ProgramData(SC_EBSDmaster)
call HDFnames%set_NMLlist(SC_EBSDmasterNameList)
call HDFnames%set_NMLfilename(SC_EBSDmasterNML)

fname = EMsoft%generateFilePath('EMdatapathname',trim(dinl%masterfile))
call MPFT%setFileName(fname)
call MPFT%readMPfile(HDF, HDFnames, mpnl)

!====================================
! 3. check to make sure that this structure is actually orthorhombic...
!    (monoclinic settings will be added in a later version)
call cell%setFileName(MPDT%xtalname)
call cell%readDataHDF(SG, EMsoft, HDF)
pgnum = SG%getPGnumber()

if ((pgnum.lt.6).or.(pgnum.gt.8)) then
    call Message%printError('DIsetting','Crystal structure point group # must be 6, 7, or 8 (orthorhombic')
end if

! get the sequential space group number within the orthorhombic system
orthonum = sgnum - SGPG(6) + 1

! OK, this is an orthorhombic structure, so we will first convert the EulerAngle representation to
! rotation matrices, then apply the appropriate permutation matrix, copy the EulerAngles array to a new
! EulerAnglesOriginal array, overwrite the Phi1, Phi, and Phi2 arrays with the new values, write the
! new EulerAngles along with a data set attribute to indicate that this is a derived array, and
! finally generate the appropriate .ctf files

call Message%printMessage(' The original orthorhombic setting is '//extendedOrthsettings(1)//' for '//SYM_SGname(sgnum))
outstring = ' The requested setting is '//extendedOrthsettings(csnl%orthorhombicSetting)// &
            ' for '//extendedHMOrthsymbols(csnl%orthorhombicSetting,orthonum)
call Message%printMessage(outstring)
call Message%printMessage(' ')
call Message%printMessage(' Starting conversion of orientation data to new setting ')
call Message%printMessage('  (original EulerAngles array will be renamed to EulerAnglesOriginal)')
call Message%printMessage(' ')

! select the quaternion that will carry out the transformation to the new setting
select case(csnl%orthorhombicSetting)
  case(1)
    eu = e_T( edinp = (/ 0.D0, 0.D0, 0.D0 /) )
  case(2)
    eu = e_T( edinp = (/ cPi*0.5, cPi, 0.D0/) )
  case(3)
    eu = e_T( edinp = (/ cPi*0.5D0, cPi*0.5D0, 0.D0/) )
  case(4)
    eu = e_T( edinp = (/ cPi*1.5D0, cPi*0.5D0, cPi*0.5D0/) )
  case(5)
    eu = e_T( edinp = (/ cPi, cPi*0.5D0, cPi*0.5D0/) )
  case(6)
    eu = e_T( edinp = (/ 0.0D0, cPi*0.5D0, 0.D0/) )
  case default
    call Message%printError('EMEBSDDIchangesetting','Unknown orthorhombic setting; please check for typos')
end select
qrot = eu%eq()
qr = Quaternion_T( qd = qrot%q_copyd() )
qr = conjg(qr)

! parallel section starts here
call OMP_setNThreads(csnl%nthreads)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID, i, qin, eu, qrin, qnew)
TID = OMP_GET_THREAD_NUM()

!$OMP DO SCHEDULE(DYNAMIC)
do i=1,DIDT%FZcnt
  eu = e_T( edinp = newEulers(1:3,i) )
  qin = eu%eq()
  qrin = Quaternion_T( qd = qin%q_copyd() )
  qnew = qr*qrin
  qin = q_T( qdinp = qnew%get_quatd() )
  eu = qin%qe()
  newEulers(1:3,i) = eu%e_copyd()
end do
!$OMP END DO

if (TID.eq.0) call Message%printMessage('  -> completed transformation of EulerAngles array')

! rotate the average orientations to the new setting

!$OMP DO SCHEDULE(DYNAMIC)
do i=1,DIDT%Nexp
  eu = e_T( edinp = newAvOr(1:3,i) )
  qin = eu%eq()
  qrin = Quaternion_T( qd = qin%q_copyd() )
  qnew = qr*qrin
  qin = q_T( qdinp = qnew%get_quatd() )
  eu = qin%qe()
  newAvOr(1:3,i) = eu%e_copyd()
end do
!$OMP END DO
if (TID.eq.0) call Message%printMessage('  -> completed transformation of AverageOrientations array')

if (transformRefined.eqv..TRUE.) then
!$OMP DO SCHEDULE(DYNAMIC)
  do i=1,DIDT%Nexp
    eu = e_T( edinp = dble(DIDT%RefinedEulerAngles(1:3,i)) )
    qin = eu%eq()
    qrin = Quaternion_T( qd = qin%q_copyd() )
    qnew = qr*qrin
    qin = q_T( qdinp = qnew%get_quatd() )
    eu = qin%qe()
    newRefined(1:3,i) = eu%e_copyd()
  end do
!$OMP END DO
  if (TID.eq.0) call Message%printMessage('  -> completed transformation of RefinedEulerAngles array')
end if

!$OMP END PARALLEL

newEulers = newEulers * rtod
newAvOr = newAvOr * rtod
! and ends here...

!===================================================================================
! open the dot product file
hdferr =  HDF%openFile(DIfile)

! open the Scan 1/EBSD/Data group
groupname = 'Scan 1'
hdferr = HDF%openGroup(groupname)
groupname = SC_EBSD
hdferr = HDF%openGroup(groupname)
groupname = SC_Data
hdferr = HDF%openGroup(groupname)

! next, create a new EulerAnglesOriginal array and write the original angles to its
! same thing for the AverageOrientations array
dataset = 'EulerAnglesOriginal'
call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetFloatArray(dataset, DIDT%EulerAngles, 3, DIDT%FZcnt, overwrite)
else
  hdferr = HDF%writeDatasetFloatArray(dataset, DIDT%EulerAngles, 3, DIDT%FZcnt)
end if

! then, overwrite the existing EulerAngles data set with the rotated orientations
dataset = SC_EulerAngles
hdferr = HDF%writeDatasetFloatArray(dataset, sngl(newEulers), 3, DIDT%FZcnt, overwrite)
! add an explanatory attribute to the data set

! do the same with the AverageOrientations data set
dataset = 'AverageOrientationsOriginal'
call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetFloatArray(dataset, DIDT%AverageOrientations, 3, DIDT%Nexp, overwrite)
else
  hdferr = HDF%writeDatasetFloatArray(dataset, DIDT%AverageOrientations, 3, DIDT%Nexp)
end if

! then, overwrite the existing EulerAngles data set with the rotated orientations
dataset = SC_AverageOrientations
hdferr = HDF%writeDatasetFloatArray(dataset, sngl(newAvOr), 3, DIDT%Nexp, overwrite)

if (transformRefined.eqv..TRUE.) then
  ! next, create a new RefinedEulerAnglesOriginal array and write the original angles to its
  dataset = 'RefinedEulerAnglesOriginal'
  call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists) then
    hdferr = HDF%writeDatasetFloatArray(dataset, DIDT%RefinedEulerAngles, 3, DIDT%Nexp, overwrite)
  else
    hdferr = HDF%writeDatasetFloatArray(dataset, DIDT%RefinedEulerAngles, 3, DIDT%Nexp)
  end if

  ! then, overwrite the existing EulerAngles data set with the rotated orientations
  dataset = SC_RefinedEulerAngles
  hdferr = HDF%writeDatasetFloatArray(dataset, sngl(newRefined), 3, DIDT%Nexp, overwrite)
end if

! add an explanatory data set
comment = 'Original orthorhombic setting changed to '//extendedOrthsettings(csnl%orthorhombicSetting)// &
          ', corresponding to space group symbol '//extendedHMOrthsymbols(csnl%orthorhombicSetting,orthonum)
allocate(stringarray(1))
stringarray(1)= trim(comment)
dataset = 'Comment'
call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
if (g_exists) then
  hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1, overwrite)
else
  hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)
end if
deallocate(stringarray)

! also, update the Phi1, Phi, and Phi2 data sets
newEulers = newEulers * sngl(dtor)
allocate(eulers(3,DIDT%Nexp), ang(DIDT%Nexp))
do i=1,DIDT%Nexp
  eulers(1:3,i) = newEulers(1:3,DIDT%TopMatchIndices(1,i))
end do

dataset = SC_Phi1
ang(:) = eulers(1,:)
hdferr = HDF%writeDatasetFloatArray(dataset, ang, DIDT%Nexp, overwrite)

dataset = SC_Phi
ang(:) = eulers(2,:)
hdferr = HDF%writeDatasetFloatArray(dataset, ang, DIDT%Nexp, overwrite)

dataset = SC_Phi2
ang(:) = eulers(3,:)
hdferr = HDF%writeDatasetFloatArray(dataset, ang, DIDT%Nexp, overwrite)

! leave this group and file
call HDF%pop(.TRUE.)

! finally, we need to write a new .ctf file as well...
dinl%ctffile = trim(csnl%newctffile)

ipar = 0
ipar(1) = 1
ipar(2) = DIDT%Nexp
ipar(3) = DIDT%Nexp
ipar(4) = DIDT%Nexp
ipar(5) = DIDT%FZcnt
ipar(6) = pgnum
if (sum(dinl%ROI).ne.0) then
    ipar(7) = dinl%ROI(3)
    ipar(8) = dinl%ROI(4)
else
    ipar(7) = dinl%ipf_wd
    ipar(8) = dinl%ipf_ht
end if
fpar2(1) = dinl%energymax
fpar2(2) = DIDT%MCsig

allocate(resultmain(1,ipar(2)))
resultmain(1,:) = DIDT%CI(:)

eulers = eulers / sngl(dtor)

VT = Vendor_T()
modality = 'EBSD'
call VT%set_Modality(modality)

call VT%ctf_writeFile(EMsoft,cell,SG,dinl,ipar,fpar2,DIDT%TopMatchIndices, &
                      eulers,resultmain,DIDT%OSM,DIDT%IQ,noindex=.TRUE.)
call Message%printMessage('Data stored in ctf file : '//trim(dinl%ctffile))

call closeFortranHDFInterface()

end associate

end subroutine DIsetting_



end module mod_DIsetting
