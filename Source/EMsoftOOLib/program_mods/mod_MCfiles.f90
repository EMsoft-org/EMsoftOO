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
! AND ANY EXPRESS OR IMCLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
! IMCLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMCLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
! USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! ###################################################################

module mod_MCfiles
  !! author: MDG
  !! version: 1.0 
  !! date: 02/12/20
  !!
  !! class for Monte Carlo file handling 

use mod_kinds
use mod_global
use stringconstants

IMPLICIT NONE 
private

! namelist for the EMMCOpenCL program
type, public :: MCOpenCLNameListType
  integer(kind=irg) :: stdout
  integer(kind=irg) :: numsx
  integer(kind=irg) :: ivolx 
  integer(kind=irg) :: ivoly 
  integer(kind=irg) :: ivolz 
  integer(kind=irg) :: globalworkgrpsz
  integer(kind=irg) :: num_el
  integer(kind=irg) :: totnum_el
  integer(kind=irg) :: multiplier
  integer(kind=irg) :: devid
  integer(kind=irg) :: platid
  real(kind=sgl)    :: ivolstepx 
  real(kind=sgl)    :: ivolstepy 
  real(kind=sgl)    :: ivolstepz 
  real(kind=dbl)    :: sig
  real(kind=dbl)    :: sigstart
  real(kind=dbl)    :: sigend
  real(kind=dbl)    :: sigstep
  real(kind=dbl)    :: omega
  real(kind=dbl)    :: EkeV
  real(kind=dbl)    :: Ehistmin
  real(kind=dbl)    :: Ebinsize
  real(kind=dbl)    :: depthmax
  real(kind=dbl)    :: depthstep
  real(kind=dbl)    :: thickness
  real(kind=dbl)    :: radius
  real(kind=dbl)    :: incloc
  character(3)      :: Notify
  character(4)      :: MCmode
  character(fnlen)  :: xtalname
  character(fnlen)  :: dataname
  character(fnlen)  :: mode
end type MCOpenCLNameListType

! this structure covers EBSD, ECP, and interaction volume data
type, public :: MCdataType
  integer(kind=irg)               :: numEbins
  integer(kind=irg)               :: numzbins
  integer(kind=irg)               :: numangle
  integer(kind=irg),allocatable   :: accum_e(:,:,:)
  integer(kind=irg),allocatable   :: accum_z(:,:,:,:)
  real(kind=sgl),allocatable      :: accumSP(:,:,:)
  integer(kind=irg),allocatable   :: accum_xyz(:,:,:)
end type MCdataType

type, public :: MCfile_T 
  private 
    type(MCdataType),public             :: MCDT
    type(MCOpenCLNameListType),public   :: nml 
    character(fnlen)                    :: MCfile

  contains
  private 

    procedure, pass(self) :: readMCfile_
    procedure, pass(self) :: writeMCfile_
    procedure, pass(self) :: writeHDFNameList_
    procedure, pass(self) :: copynml_
    procedure, pass(self) :: copyaccume_
    procedure, pass(self) :: copyaccumz_
    procedure, pass(self) :: copyaccumSP_
    procedure, pass(self) :: copyaccumxyz_
    procedure, pass(self) :: getnml_
    procedure, pass(self) :: getnumEbins_
    procedure, pass(self) :: getnumzbins_
    procedure, pass(self) :: getnumangles_
    procedure, pass(self) :: setFileName_
    procedure, pass(self) :: copyMCdata_

    generic, public :: readMCfile => readMCfile_
    generic, public :: writeHDFNameList => writeHDFNameList_
    generic, public :: writeMCfile => writeMCfile_
    generic, public :: copynml => copynml_
    generic, public :: copyaccume => copyaccume_
    generic, public :: copyaccumz => copyaccumz_
    generic, public :: copyaccumSP => copyaccumSP_
    generic, public :: copyaccumxyz => copyaccumxyz_
    generic, public :: getnml => getnml_
    generic, public :: getnumEbins => getnumEbins_
    generic, public :: getnumzbins => getnumzbins_
    generic, public :: getnumangles => getnumangles_
    generic, public :: setFileName => setFileName_
    generic, public :: copyMCdata => copyMCdata_

end type MCfile_T

! the constructor routine for this class 
interface MCfile_T
  module procedure MCfile_constructor
end interface MCfile_T

contains 

!--------------------------------------------------------------------------
type(MCfile_T) function MCfile_constructor( ) result(MCfile)
!! author: MDG 
!! version: 1.0 
!! date: 02/05/20
!!
!! constructor for the MCfile_T Class
 
IMPLICIT NONE


end function MCfile_constructor

!--------------------------------------------------------------------------
subroutine copynml_(self, nml)
!! author: MDG 
!! version: 1.0 
!! date: 02/06/20
!!
!! copy the namelist into the MCfile_T class for writing to file

IMPLICIT NONE 

class(MCfile_T), INTENT(INOUT)          :: self
type(MCOpenCLNameListType),INTENT(IN)   :: nml

self%nml = nml

end subroutine copynml_

!--------------------------------------------------------------------------
subroutine copyaccume_(self, acc, keep)
!! author: MDG 
!! version: 1.0 
!! date: 02/12/20
!!
!! copy the accumulator array

IMPLICIT NONE 

class(MCfile_T), INTENT(INOUT)                :: self
integer(kind=irg), allocatable, INTENT(OUT)   :: acc(:,:,:)
logical, INTENT(IN), OPTIONAL                 :: keep 

integer(kind=irg)                             :: s(3), nx

s = shape(self%MCDT%accum_e)
nx = (s(2)-1)/2
allocate(acc(s(1),-nx:nx,-nx:nx))

acc = self%MCDT%accum_e 

if (.not.present(keep)) deallocate(self%MCDT%accum_e)

end subroutine copyaccume_

!--------------------------------------------------------------------------
subroutine copyaccumz_(self, acc, keep)
!! author: MDG 
!! version: 1.0 
!! date: 02/12/20
!!
!! copy the accumulator array

IMPLICIT NONE 

class(MCfile_T), INTENT(INOUT)                :: self
integer(kind=irg), allocatable, INTENT(OUT)   :: acc(:,:,:,:)
logical, INTENT(IN), OPTIONAL                 :: keep 

integer(kind=irg)                             :: s(4), nx, ny

s = shape(self%MCDT%accum_z)
nx = (s(3)-1)/2
ny = (s(4)-1)/2
allocate( acc(s(1),s(2),-nx:nx,-ny:ny) )

acc = self%MCDT%accum_z 

if (.not.present(keep)) deallocate(self%MCDT%accum_z)

end subroutine copyaccumz_

!--------------------------------------------------------------------------
subroutine copyaccumSP_(self, acc, keep)
!! author: MDG 
!! version: 1.0 
!! date: 02/12/20
!!
!! copy the accumulator array

IMPLICIT NONE 

class(MCfile_T), INTENT(INOUT)                :: self
integer(kind=irg), allocatable, INTENT(OUT)   :: acc(:,:,:)
logical, INTENT(IN), OPTIONAL                 :: keep 

integer(kind=irg)                             :: s(3)

s = shape(self%MCDT%accumSP)
allocate(acc(s(1),s(2),s(3)))

acc = self%MCDT%accumSP 

if (.not.present(keep)) deallocate(self%MCDT%accumSP)

end subroutine copyaccumSP_

!--------------------------------------------------------------------------
subroutine copyaccumxyz_(self, acc, keep)
!! author: MDG 
!! version: 1.0 
!! date: 02/12/20
!!
!! copy the accumulator array

IMPLICIT NONE 

class(MCfile_T), INTENT(INOUT)                :: self
integer(kind=irg), allocatable, INTENT(OUT)   :: acc(:,:,:)
logical, INTENT(IN), OPTIONAL                 :: keep 

integer(kind=irg)                             :: s(3), nx, ny

s = shape(self%MCDT%accum_xyz)
nx = (s(1)-1)/2
ny = (s(2)-1)/2
allocate(acc(-nx:nx,-ny:ny,s(3)))

acc = self%MCDT%accum_xyz 

if (.not.present(keep)) deallocate(self%MCDT%accum_xyz)

end subroutine copyaccumxyz_

!--------------------------------------------------------------------------
function getnml_(self) result(nml)
!! author: MDG 
!! version: 1.0 
!! date: 02/06/20
!!
!! get the namelist from the MCfile_T class

IMPLICIT NONE 

class(MCfile_T), INTENT(INOUT)          :: self
type(MCOpenCLNameListType)              :: nml

nml = self%nml

end function getnml_

!--------------------------------------------------------------------------
function getnumEbins_(self) result(n)
!! author: MDG 
!! version: 1.0 
!! date: 02/06/20
!!
!! get numEbins from the MCfile_T class

IMPLICIT NONE 

class(MCfile_T), INTENT(INOUT)    :: self
integer(kind=irg)                 :: n

n = self%MCDT%numEbins

end function getnumEbins_

!--------------------------------------------------------------------------
function getnumzbins_(self) result(n)
!! author: MDG 
!! version: 1.0 
!! date: 02/06/20
!!
!! get numzbins from the MCfile_T class

IMPLICIT NONE 

class(MCfile_T), INTENT(INOUT)    :: self
integer(kind=irg)                 :: n

n = self%MCDT%numzbins

end function getnumzbins_

!--------------------------------------------------------------------------
function getnumangles_(self) result(n)
!! author: MDG 
!! version: 1.0 
!! date: 02/06/20
!!
!! get numangles from the MCfile_T class

IMPLICIT NONE 

class(MCfile_T), INTENT(INOUT)    :: self
integer(kind=irg)                 :: n

n = self%MCDT%numangle

end function getnumangles_

!--------------------------------------------------------------------------
subroutine setFileName_(self, MCfile)
!! author: MDG 
!! version: 1.0 
!! date: 02/06/20
!!
!! set the Monte Carlo file name 

IMPLICIT NONE 

class(MCfile_T), INTENT(INOUT)          :: self
character(fnlen), INTENT(IN)            :: MCfile

self%MCfile = trim(MCfile)

end subroutine setFileName_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF)
!! author: MDG 
!! version: 1.0 
!! date: 01/22/20
!!
!! write namelist to HDF file

use mod_HDFsupport
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(MCfile_T), INTENT(INOUT)          :: self 
type(HDF_T), INTENT(INOUT)              :: HDF

integer(kind=irg),parameter             :: n_int = 11, n_real_bse1 = 9, n_real_full = 7, n_real_ivol= 6
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=dbl)                          :: io_real_bse1(n_real_bse1), io_real_full(n_real_full), &
                                           io_real_ivol(n_real_ivol)
character(20)                           :: reallist_bse1(n_real_bse1), reallist_full(n_real_full), &
                                           reallist_ivol(n_real_ivol)
character(20)                           :: intlist(n_int)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( mcnl => self%nml )

! create the group for this namelist
groupname = SC_MCCLNameList
hdferr = HDF%createGroup(groupname)

! write all the single integers
io_int = (/ mcnl%stdout, mcnl%numsx, mcnl%globalworkgrpsz, mcnl%num_el, mcnl%totnum_el, mcnl%multiplier, mcnl%devid, &
            mcnl%platid, mcnl%ivolx, mcnl%ivoly, mcnl%ivolz /)
intlist(1) = 'stdout'
intlist(2) = 'numsx'
intlist(3) = 'globalworkgrpsz'
intlist(4) = 'num_el'
intlist(5) = 'totnum_el'
intlist(6) = 'multiplier'
intlist(7) = 'devid'
intlist(8) = 'platid'
intlist(9) = 'ivolx'
intlist(10) = 'ivoly'
intlist(11) = 'ivolz'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single doubles
if (mcnl%mode .eq. 'bse1') then
   io_real_bse1 = (/ mcnl%sigstart, mcnl%sigend, mcnl%sigstep, mcnl%omega, mcnl%EkeV, mcnl%Ehistmin, &
             mcnl%Ebinsize, mcnl%depthmax, mcnl%depthstep /)
   reallist_bse1(1) = 'sigstart'
   reallist_bse1(2) = 'sigend'
   reallist_bse1(3) = 'sigstep'
   reallist_bse1(4) = 'omega'
   reallist_bse1(5) = 'EkeV'
   reallist_bse1(6) = 'Ehistmin'
   reallist_bse1(7) = 'Ebinsize'
   reallist_bse1(8) = 'depthmax'
   reallist_bse1(9) = 'depthstep'
   call HDF%writeNMLdbles(io_real_bse1, reallist_bse1, n_real_bse1)
else if (mcnl%mode .eq. 'full') then
   io_real_full = (/ mcnl%sig, mcnl%omega, mcnl%EkeV, mcnl%Ehistmin, &
             mcnl%Ebinsize, mcnl%depthmax, mcnl%depthstep /)
   reallist_full(1) = 'sig'
   reallist_full(2) = 'omega'
   reallist_full(3) = 'EkeV'
   reallist_full(4) = 'Ehistmin'
   reallist_full(5) = 'Ebinsize'
   reallist_full(6) = 'depthmax'
   reallist_full(7) = 'depthstep'
   call HDF%writeNMLdbles(io_real_full, reallist_full, n_real_full)
else if (mcnl%mode .eq. 'Ivol') then
   io_real_ivol = (/ mcnl%sig, mcnl%omega, mcnl%EkeV, dble(mcnl%ivolstepx), dble(mcnl%ivolstepy), dble(mcnl%ivolstepz) /)
   reallist_ivol(1) = 'sig'
   reallist_ivol(2) = 'omega'
   reallist_ivol(3) = 'EkeV'
   reallist_ivol(4) = 'ivolstepx'
   reallist_ivol(5) = 'ivolstepy'
   reallist_ivol(6) = 'ivolstepz'
   call HDF%writeNMLdbles(io_real_ivol, reallist_ivol, n_real_ivol)
end if

! write all the strings
dataset = SC_MCmode
sval(1) = mcnl%MCmode
hdferr = HDF%writeDatasetStringArray(dataset, sval, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create MCmode dataset', hdferr)

dataset = SC_xtalname
line2(1) = mcnl%xtalname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create xtalname dataset', hdferr)

dataset = SC_dataname
line2(1) = mcnl%dataname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create dataname dataset', hdferr)

dataset = SC_mode
sval(1) = mcnl%mode
hdferr = HDF%writeDatasetStringArray(dataset, sval, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create mode dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
recursive subroutine readMCfile_(self, HDF, getAccume, getAccumz, getAccumSP, getAccumxyz) 
!! author: MDG 
!! version: 1.0 
!! date: 02/05/20
!!
!! read namelist and selected data from Monte Carlo file 

use HDF5 
use mod_HDFsupport
use mod_io
use stringconstants
use ISO_C_BINDING 

IMPLICIT NONE 

class(MCfile_T), INTENT(INOUT)  :: self 
type(HDF_T),INTENT(INOUT)       :: HDF
logical,INTENT(IN),OPTIONAL     :: getAccume     
 !! energy accumulator array switch
logical,INTENT(IN),OPTIONAL     :: getAccumz     
 !! depth accumulator array switch
logical,INTENT(IN),OPTIONAL     :: getAccumSP    
 !! stereographic accumulator array switch
logical,INTENT(IN),OPTIONAL     :: getAccumxyz   
 !! interaction volume array switch

type(IO_T)                                          :: Message
character(fnlen)                                    :: groupname, datagroupname, dataset
logical                                             :: stat, readonly, g_exists, f_exists, FL
integer(kind=irg)                                   :: ii, nlines, nx, ny, nz, hdferr
real(kind=dbl)                                      :: x
integer(kind=irg),allocatable                       :: iarray(:)
real(kind=sgl),allocatable                          :: farray(:)
integer(kind=irg),allocatable                       :: accum_e(:,:,:)
integer(kind=irg),allocatable                       :: accum_xyz(:,:,:)
integer(kind=irg),allocatable                       :: accum_z(:,:,:,:)
integer(HSIZE_T)                                    :: dims(1), dims2(2), dims3(3), offset3(3), dims4(4) 
character(fnlen, KIND=c_char),allocatable,TARGET    :: stringarray(:)

associate( nml => self%nml, MCDT => self%MCDT )

! we assume that calling program has set the MCfile variable
! using the setFileName(MCfile) method (full path name)

inquire(file=trim(self%MCfile), exist=f_exists)

if (.not.f_exists) then
  call Message%printError('readMCfile','Monte Carlo input file does not exist')
end if

! is this a proper HDF5 file ?
call h5fis_hdf5_f(trim(self%MCfile), stat, hdferr)

if (stat.eqv..FALSE.) then ! the file exists, so let's open it an first make sure it is an EBSD dot product file
   call Message%printError('readMCfile','This is not a proper HDF5 file')
end if 
   
! open the Monte Carlo file 
readonly = .TRUE.
hdferr =  HDF%openFile(self%MCfile, readonly)

! check whether or not the MC file was generated using DREAM.3D
! this is necessary so that the proper reading of fixed length vs. variable length strings will occur.
! this test sets a flag in side the HDFsupport module so that the proper reading routines will be employed
datagroupname = '/EMheader/MCOpenCL'
call H5Lexists_f(HDF%getobjectID(),trim(datagroupname),g_exists, hdferr)
if (.not.g_exists) then
  call Message%printError('readMCfile','This HDF file does not contain Monte Carlo header data')
end if

groupname = SC_EMheader
hdferr = HDF%openGroup(groupname)
groupname = SC_MCOpenCL
hdferr = HDF%openGroup(groupname)
FL = .FALSE.
datagroupname = 'FixedLength'
FL = HDF%CheckFixedLengthflag(datagroupname)
if (FL.eqv..TRUE.) then 
  call Message%printMessage('Input file was generated by a program using fixed length strings')
end if
call HDF%pop()
call HDF%pop()

!====================================
! make sure this is a Monte Carlo file
!====================================
groupname = SC_NMLfiles
    hdferr = HDF%openGroup(groupname )
dataset = 'MCOpenCLNML'
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists.eqv..FALSE.) then
    call HDF%pop(.TRUE.)
    call Message%printError('readMCfile','this is not an EBSD Monte Carlo file')
end if
call HDF%pop()

!====================================
! read all NMLparameters group datasets
!====================================
groupname = SC_NMLparameters
    hdferr = HDF%openGroup(groupname)
groupname = SC_MCCLNameList
    hdferr = HDF%openGroup(groupname)

! we'll read these roughly in the order that the HDFView program displays them...

dataset = SC_mode
    call HDF%readDatasetStringArray(dataset, nlines, hdferr, stringarray)
    nml%mode = trim(stringarray(1))
    deallocate(stringarray)

dataset = 'ivolx'
  call HDF%readDatasetInteger(dataset, hdferr, nml%ivolx) 

dataset = 'ivoly'
  call HDF%readDatasetInteger(dataset, hdferr, nml%ivoly) 

dataset = 'ivolz'
  call HDF%readDatasetInteger(dataset, hdferr, nml%ivolz) 

dataset = 'ivolstepx'
  call HDF%readDatasetDouble(dataset, hdferr, x)
  nml%ivolstepx = sngl(x)

dataset = 'ivolstepy'
  call HDF%readDatasetDouble(dataset, hdferr, x)
  nml%ivolstepy = sngl(x)

dataset = 'ivolstepz'
  call HDF%readDatasetDouble(dataset, hdferr, x)
  nml%ivolstepz = sngl(x)

dataset = SC_Ebinsize
  call HDF%readDatasetDouble(dataset, hdferr, nml%Ebinsize)

dataset = SC_Ehistmin
  call HDF%readDatasetDouble(dataset, hdferr, nml%Ehistmin)

dataset = SC_depthmax
  call HDF%readDatasetDouble(dataset, hdferr, nml%depthmax)

dataset = SC_depthstep
  call HDF%readDatasetDouble(dataset, hdferr, nml%depthstep)

dataset = SC_EkeV
    call HDF%readDatasetDouble(dataset, hdferr, nml%EkeV)

dataset = SC_dataname
    call HDF%readDatasetStringArray(dataset, nlines, hdferr, stringarray)
    nml%dataname = trim(stringarray(1))
    deallocate(stringarray)

dataset = SC_MCmode
    call HDF%readDatasetStringArray(dataset, nlines, hdferr, stringarray)
    nml%MCmode = trim(stringarray(1))
    deallocate(stringarray)

dataset = SC_devid
    call HDF%readDatasetInteger(dataset, hdferr, nml%devid)

dataset = SC_globalworkgrpsz
    call HDF%readDatasetInteger(dataset, hdferr, nml%globalworkgrpsz)

dataset = SC_multiplier
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetInteger(dataset, hdferr, nml%multiplier)
else
    nml%multiplier  = 1
end if

dataset = 'num_el'
    call HDF%readDatasetInteger(dataset, hdferr, nml%num_el)

dataset = SC_numsx
    call HDF%readDatasetInteger(dataset, hdferr, nml%numsx)

dataset = SC_omega
    call HDF%readDatasetDouble(dataset, hdferr, nml%omega)

dataset = SC_platid
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetInteger(dataset, hdferr, nml%platid)
else
    nml%platid  = 1
end if

dataset = SC_sig
    call HDF%readDatasetDouble(dataset, hdferr, nml%sig)

dataset = SC_stdout
    call HDF%readDatasetInteger(dataset, hdferr, nml%stdout)

dataset = SC_totnumel
    call HDF%readDatasetInteger(dataset, hdferr, nml%totnum_el)

dataset = SC_xtalname
    call HDF%readDatasetStringArray(dataset, nlines, hdferr, stringarray)
    nml%xtalname = trim(stringarray(1))
    deallocate(stringarray)

! and close the NMLparameters group
    call HDF%pop()
    call HDF%pop()

!====================================
!====================================

! open the Monte Carlo data group
groupname = SC_EMData
    hdferr = HDF%openGroup(groupname)
groupname = SC_MCOpenCL
    hdferr = HDF%openGroup(groupname)

! integers
dataset = SC_multiplier
if (nml%multiplier.eq.1) then 
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  if (g_exists.eqv..TRUE.) then
      call HDF%readDatasetInteger(dataset, hdferr, nml%multiplier)
  else
      nml%multiplier = 1
  end if
end if 

if (trim(nml%mode).eq.'full') then 
  dataset = SC_numEbins
    call HDF%readDatasetInteger(dataset, hdferr, MCDT%numEbins)
  dataset = SC_numzbins
    call HDF%readDatasetInteger(dataset, hdferr, MCDT%numzbins)
else if (trim(nml%mode).eq.'bse1') then 
  dataset = SC_numangle
    call HDF%readDatasetInteger(dataset, hdferr, MCDT%numangle)
  dataset = SC_numzbins
    call HDF%readDatasetInteger(dataset, hdferr, MCDT%numzbins)
end if


! various optional arrays
if (present(getAccume)) then 
  if (getAccume.eqv..TRUE.) then
    dataset = SC_accume
    call HDF%readDatasetIntegerArray(dataset, dims3, hdferr, accum_e)
    nx = (dims3(2)-1)/2
    allocate(MCDT%accum_e(1:dims3(1),-nx:nx,-nx:nx))
    MCDT%accum_e = accum_e
    deallocate(accum_e)
  end if 
end if

if (present(getAccumz)) then 
  if (getAccumz.eqv..TRUE.) then
    dataset = SC_accumz
    call HDF%readDatasetIntegerArray(dataset, dims4, hdferr, accum_z)
    nx = (dims4(3)-1)/2
    allocate(MCDT%accum_z(1:dims4(1),1:dims4(2),-nx:nx, -nx:nx))
    MCDT%accum_z = accum_z
    deallocate(accum_z)  
  end if 
end if

if (present(getAccumSP)) then 
  if (getAccumSP.eqv..TRUE.) then
    dataset = SC_accumSP
    call HDF%readDatasetFloatArray(dataset, dims3, hdferr, MCDT%accumSP)
  end if 
end if

if (present(getAccumxyz)) then 
  if (getAccumxyz.eqv..TRUE.) then
    dataset = SC_accumxyz
    call HDF%readDatasetIntegerArray(dataset, dims3, hdferr, accum_xyz)
    nx = (dims3(1)-1)/2
    ny = (dims3(2)-1)/2
    nz = dims3(3)
    allocate(MCDT%accum_xyz(-nx:nx, -ny:ny, nz))
    MCDT%accum_xyz = accum_xyz
    deallocate(accum_xyz)  
  end if 
end if

! and close the HDF5 Monte Carloe file
call HDF%pop(.TRUE.)

call Message%printMessage(' -> completed reading Monte Carlo data from '//trim(self%MCfile), frm = "(A/)")

end associate 

end subroutine readMCfile_

!--------------------------------------------------------------------------
recursive subroutine writeMCfile_(self, EMsoft, cell, SG, HDF, progname, dstr, tstrb, tstre)
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/06/20
  !!
  !! write Monte Carlo data to an HDF5 file

use mod_EMsoft
use mod_crystallography
use mod_symmetry
use HDF5
use mod_HDFsupport
use stringconstants

IMPLICIT NONE 

class(MCfile_T), INTENT(INOUT)      :: self 
type(EMsoft_T), INTENT(INOUT)       :: EMsoft 
type(Cell_T), INTENT(INOUT)         :: cell 
type(SpaceGroup_T), INTENT(INOUT)   :: SG
type(HDF_T),INTENT(INOUT)           :: HDF 
character(fnlen), INTENT(IN)        :: progname
character(11), INTENT(INOUT)        :: dstr
character(15), INTENT(IN)           :: tstrb
character(15), INTENT(IN)           :: tstre

character(fnlen)                    :: dataname, datagroupname, groupname, attributename, dataset
logical                             :: f_exists
integer(kind=irg)                   :: hdferr, nx, s(3), s4(4)  
character(fnlen,kind=c_char)        :: HDF_FileVersion

associate ( nml => self%nml, MCDT => self%MCDT )

! get the filename; if it already exists, then delete it and create a new one
dataname = EMsoft%generateFilePath('EMdatapathname', nml%dataname)
inquire(file=trim(dataname), exist=f_exists)

if (f_exists) then
  open(unit=dataunit, file=trim(dataname), status='old',form='unformatted')
  close(unit=dataunit, status='delete')
end if

! Create a new file using the default properties.
hdferr = HDF%createFile(dataname)

! write the EMheader to the file
datagroupname = 'MCOpenCL'
call HDF%writeEMheader(dstr, tstrb, tstre, progname, datagroupname)

! add the CrystalData group at the top level of the file
call cell%addXtalDataGroup(SG, EMsoft, HDF)

! create a namelist group to write all the namelist files into
groupname = SC_NMLfiles
hdferr = HDF%createGroup(groupname)

! read the text file and write the array to the file
dataset = SC_MCOpenCLNML
hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%nmldeffile)

! leave this group
call HDF%pop()

! create a namelist group to write all the namelist files into
groupname = SC_NMLparameters
hdferr = HDF%createGroup(groupname)
call self%writeHDFNameList(HDF)

! leave this group
call HDF%pop()

! then the remainder of the data in a EMData group
groupname = SC_EMData
hdferr = HDF%createGroup(groupname)

! here we add the data groupname MCOpenCL and we attach to it a HDF_FileVersion attribute 
hdferr = HDF%createGroup(datagroupname)
HDF_FileVersion = '4.0'
HDF_FileVersion = cstringify(HDF_FileVersion)
attributename = SC_HDFFileVersion
hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)

! =====================================================
! The following write commands constitute HDF_FileVersion = 4.0
! =====================================================

dataset = SC_numzbins
hdferr = HDF%writeDatasetInteger(dataset, MCDT%numzbins)

dataset = SC_numEbins
hdferr = HDF%writeDatasetInteger(dataset, MCDT%numEbins)

! modified using multiplier
dataset = SC_totnumel
hdferr = HDF%writeDatasetInteger(dataset, nml%totnum_el)

dataset = SC_multiplier
hdferr = HDF%writeDatasetInteger(dataset, nml%multiplier)

if (nml%mode .eq. 'full') then
  s = shape(MCDT%accum_e)
  nx = s(2)
  
  dataset = SC_numEbins
      hdferr = HDF%writeDatasetInteger(dataset, MCDT%numEbins)
  
  dataset = SC_accume
      hdferr = HDF%writeDatasetIntegerArray(dataset, MCDT%accum_e, MCDT%numEbins, nx, nx)
  
  dataset = SC_accumSP
      hdferr = HDF%writeDatasetFloatArray(dataset, MCDT%accumSP, MCDT%numEbins, nx, nx)
  
  s4 = shape(MCDT%accum_z)
  nx = s4(3)
  
  dataset = SC_accumz
      hdferr = HDF%writeDatasetIntegerArray(dataset, MCDT%accum_z, MCDT%numEbins, MCDT%numzbins, nx, nx)
  

else if (nml%mode .eq. 'bse1') then
  s = shape(MCDT%accum_e)
  nx = s(2)
  
  dataset = SC_numangle
      hdferr = HDF%writeDatasetInteger(dataset, MCDT%numangle)
  
  dataset = SC_accume
      hdferr = HDF%writeDatasetIntegerArray(dataset, MCDT%accum_e, MCDT%numangle, nx, nx)
  
  s4 = shape(MCDT%accum_z)
  nx = s4(3)

  dataset = SC_accumz
      hdferr = HDF%writeDatasetIntegerArray(dataset, MCDT%accum_z, MCDT%numangle, MCDT%numzbins, nx, nx)
  
else if (nml%mode .eq. 'Ivol') then

  s = shape(MCDT%accum_xyz)

  dataset = SC_accumxyz
      hdferr = HDF%writeDatasetIntegerArray(dataset, MCDT%accum_xyz, s(1), s(2), s(3))

end if

! =====================================================
! end of HDF_FileVersion = 4.0 write statements
! =====================================================

call HDF%pop(.TRUE.)

end associate 

end subroutine writeMCfile_

!--------------------------------------------------------------------------
recursive subroutine copyMCdata_(self, EMsoft, HDF, inputfile, outputfile, h5)
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/06/20
  !!
  !! copy Monte Carlo data from one file to a new file using h5copy

use mod_EMsoft
use HDF5
use mod_HDFsupport
use mod_io
use stringconstants

IMPLICIT NONE

class(MCfile_T),INTENT(INOUT)     :: self
type(EMsoft_T),INTENT(INOUT)      :: EMsoft
type(HDF_T),INTENT(INOUT)         :: HDF 
character(fnlen),INTENT(IN)       :: inputfile
character(fnlen),INTENT(IN)       :: outputfile
character(fnlen),INTENT(IN)       :: h5

type(IO_T)                        :: Message
character(fnlen)                  :: infile, outfile, h5copypath, groupname
character(512)                    :: cmd, cmd2
logical                           :: f_exists, readonly, developer
integer(kind=irg)                 :: hdferr
character(fnlen)                  :: dev

! first we make sure that we actually have the h5copy program available
! check for EMDevelop parameter 
developer = .FALSE.
dev = EMsoft%getConfigParameter('Develop')
if (trim(dev).eq.'Yes') developer = .TRUE.

if (developer.eqv..TRUE.) then 
! if TRUE, use EMsoft_geth5copypath which is defined at configure time 
  h5copypath = trim(EMsoft%getConfigParameter('h5copypath'))//' -p -v '
  h5copypath = EMsoft%toNativePath(h5copypath)
else 
! if FALSE, check name list h5copypath parameter 
  if (trim(h5).ne.'undefined') then 
    h5copypath = trim(h5)//' -p -v '
    h5copypath = EMsoft%toNativePath(h5copypath)
  else 
! if undefined, then fail
    call Message%printError('copyMCdata','h5copypath must be set in the name list file ')
  end if
end if

call Message%printMessage(' Using '//trim(h5copypath)//' to copy Monte Carlo data to new file')

! first make sure that the input file exists and has MC data in it
infile = trim(EMsoft%generateFilePath('EMdatapathname',inputfile))
inquire(file=trim(infile), exist=f_exists)

outfile = trim(EMsoft%generateFilePath('EMdatapathname',outputfile))

! if the file does not exist, abort the program with an error message
if (f_exists.eqv..FALSE.) then 
  call Message%printError('copyMCdata','Monte Carlo copyfromenergyfile does not exist: '//trim(infile))
end if

! make sure it has MCopenCL data in it; hdf open is done in the calling program
readonly = .TRUE.
hdferr =  HDF%openFile(infile, readonly)

groupname = SC_EMData
hdferr = HDF%openGroup(groupname)
if (hdferr.eq.-1) then 
  call Message%printError('copyMCdata','EMData group does not exist in '//trim(infile))
end if

groupname = SC_MCOpenCL
hdferr = HDF%openGroup(groupname)
if (hdferr.eq.-1) then 
  call Message%printError('copyMCdata','MCOpenCL group does not exist in '//trim(infile))
end if

call HDF%pop(.TRUE.)

! OK, if we get here, then the file does exist and it contains Monte Carlo data, so we let
! the user know
call Message%printMessage('--> Input file contains Monte Carlo data')

! next, we copy the necessary groups into the new Monte Carlo file
cmd = trim(h5copypath)//' -i "'//trim(infile)
cmd = trim(cmd)//'" -o "'//trim(outfile)

cmd2 = trim(cmd)//'" -s "/CrystalData" -d "/CrystalData"'
call system(trim(cmd2))

cmd2 = trim(cmd)//'" -s "/EMData/MCOpenCL" -d "/EMData/MCOpenCL"'
call system(trim(cmd2))

cmd2 = trim(cmd)//'" -s "/EMheader/MCOpenCL" -d "/EMheader/MCOpenCL"'
call system(trim(cmd2))

cmd2 = trim(cmd)//'" -s "/NMLfiles/MCOpenCLNML" -d "/NMLfiles/MCOpenCLNML"'
call system(trim(cmd2))

cmd2 = trim(cmd)//'" -s "/NMLparameters/MCCLNameList" -d "/NMLparameters/MCCLNameList"'
call system(trim(cmd2))

call Message%printMessage('--> Output file generated with Monte Carlo data copied from '//trim(infile))

end subroutine copyMCdata_




end module mod_MCfiles