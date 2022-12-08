!! ###################################################################
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

module mod_ISE
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/08/22
  !!
  !! class definition for the EMISE program

use mod_kinds
use mod_global
use mod_ISEmaster

IMPLICIT NONE 

type, private :: ISEAccumType
  integer(kind=irg)               :: numdet
  integer(kind=irg)               :: ipf_wd
  integer(kind=irg)               :: ipf_ht
  real(kind=sgl),allocatable      :: mLPNH(:,:), mLPSH(:,:)
  real(kind=sgl),allocatable      :: masterSPNH(:,:), masterSPSH(:,:)
end type ISEAccumType

! namelist for the EMISE program
type, public :: ISENameListType
  real(kind=sgl)    :: gammavalue
  integer(kind=irg) :: nthreads
  character(3)      :: scalingmode
  character(fnlen)  :: useangles
  character(fnlen)  :: masterfile
  character(fnlen)  :: datafile
  character(fnlen)  :: imagefile
end type ISENameListType

! class definition
type, public :: ISE_T
private 
  character(fnlen)            :: nmldeffile = 'EMISE.nml'
  type(ISENameListType)       :: nml 
  type(ISEmasterNameListType) :: mpnml
  character(fnlen)            :: ISEMPfile
  type(ISEAccumType)          :: det
  character(fnlen, KIND=c_char),allocatable   :: nmlstrings(:)

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: ISE_
  procedure, pass(self) :: setgammavalue_
  procedure, pass(self) :: getgammavalue_
  procedure, pass(self) :: setnthreads_
  procedure, pass(self) :: getnthreads_
  procedure, pass(self) :: setscalingmode_
  procedure, pass(self) :: getscalingmode_
  procedure, pass(self) :: setuseangles_
  procedure, pass(self) :: getuseangles_
  procedure, pass(self) :: setISEMPfile_
  procedure, pass(self) :: getISEMPfile_
  procedure, pass(self) :: setmasterfile_
  procedure, pass(self) :: getmasterfile_
  procedure, pass(self) :: setdatafile_
  procedure, pass(self) :: getdatafile_
  procedure, pass(self) :: setimagefile_
  procedure, pass(self) :: getimagefile_
  procedure, pass(self) :: readISEMPfile_
  procedure, pass(self) :: ComputeISEimage_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: ISE => ISE_
  generic, public :: setgammavalue => setgammavalue_
  generic, public :: getgammavalue => getgammavalue_
  generic, public :: setnthreads => setnthreads_
  generic, public :: getnthreads => getnthreads_
  generic, public :: setscalingmode => setscalingmode_
  generic, public :: getscalingmode => getscalingmode_
  generic, public :: setuseangles => setuseangles_
  generic, public :: getuseangles => getuseangles_
  generic, public :: setISEMPfile => setISEMPfile_
  generic, public :: getISEMPfile => getISEMPfile_
  generic, public :: setmasterfile => setmasterfile_
  generic, public :: getmasterfile => getmasterfile_
  generic, public :: setdatafile => setdatafile_
  generic, public :: getdatafile => getdatafile_
  generic, public :: setimagefile => setimagefile_
  generic, public :: getimagefile => getimagefile_
  generic, public :: readISEMPfile => readISEMPfile_
  generic, public :: ComputeISEimage => ComputeISEimage_

end type ISE_T

! the constructor routine for this class 
interface ISE_T
  module procedure ISE_constructor
end interface ISE_T

contains

!--------------------------------------------------------------------------
type(ISE_T) function ISE_constructor( nmlfile ) result(ISE)
!! author: MDG 
!! version: 1.0 
!! date: 12/04/22
!!
!! constructor for the ISE_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call ISE%readNameList(nmlfile)

end function ISE_constructor

!--------------------------------------------------------------------------
subroutine ISE_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 12/04/22
!!
!! destructor for the ISE_T Class
 
IMPLICIT NONE

type(ISE_T), INTENT(INOUT)  :: self 

call reportDestructor('ISE_T')

end subroutine ISE_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 12/04/22
!!
!! read the namelist from an nml file for the ISE_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(ISE_T), INTENT(INOUT)          :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

  real(kind=sgl)    :: gammavalue
  integer(kind=irg) :: nthreads
  character(3)      :: scalingmode
  character(fnlen)  :: useangles
  character(fnlen)  :: masterfile
  character(fnlen)  :: datafile
  character(fnlen)  :: imagefile

! define the IO namelist to facilitate passing variables to the program.
namelist  / ISEdata / gammavalue, nthreads, useangles, scalingmode, masterfile, datafile, imagefile

! set the input parameters to default values
masterfile = 'undefined'
datafile = 'undefined'
useangles = 'original'
scalingmode = 'not'
gammavalue = 1.0
imagefile = 'undefined'
nthreads = 1

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=ISEdata)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(masterfile).eq.'undefined') then
  call Message%printError('readNameList:',' master pattern file name is undefined in '//nmlfile)
 end if

 if (trim(datafile).eq.'undefined') then
  call Message%printError('readNameList:',' output file name is undefined in '//nmlfile)
 end if

 if (trim(imagefile).eq.'undefined') then
  call Message%printError('readNameList:',' image file name is undefined in '//nmlfile)
 end if
end if 

self%nml%gammavalue = gammavalue
self%nml%nthreads = nthreads
self%nml%useangles = useangles
self%nml%scalingmode = scalingmode
self%nml%masterfile = masterfile
self%nml%datafile = datafile
self%nml%imagefile = imagefile

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 12/04/22
!!
!! pass the namelist for the ISE_T Class to the calling program

IMPLICIT NONE 

class(ISE_T), INTENT(INOUT)          :: self
type(ISENameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine setgammavalue_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setgammavalue_
!! author: MDG
!! version: 1.0
!! date: 12/08/22
!!
!! set gammavalue in the ISE_T class

IMPLICIT NONE

class(ISE_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%gammavalue = inp

end subroutine setgammavalue_

!--------------------------------------------------------------------------
function getgammavalue_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getgammavalue_
!! author: MDG
!! version: 1.0
!! date: 12/08/22
!!
!! get gammavalue from the ISE_T class

IMPLICIT NONE

class(ISE_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%gammavalue

end function getgammavalue_

!--------------------------------------------------------------------------
subroutine setnthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnthreads_
!! author: MDG
!! version: 1.0
!! date: 12/08/22
!!
!! set nthreads in the ISE_T class

IMPLICIT NONE

class(ISE_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%nthreads = inp

end subroutine setnthreads_

!--------------------------------------------------------------------------
function getnthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnthreads_
!! author: MDG
!! version: 1.0
!! date: 12/08/22
!!
!! get nthreads from the ISE_T class

IMPLICIT NONE

class(ISE_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%nthreads

end function getnthreads_

!--------------------------------------------------------------------------
subroutine setscalingmode_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setscalingmode_
!! author: MDG
!! version: 1.0
!! date: 12/08/22
!!
!! set scalingmode in the ISE_T class

IMPLICIT NONE

class(ISE_T), INTENT(INOUT)     :: self
character(3), INTENT(IN)       :: inp

self%nml%scalingmode = trim(inp)

end subroutine setscalingmode_

!--------------------------------------------------------------------------
function getscalingmode_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getscalingmode_
!! author: MDG
!! version: 1.0
!! date: 12/08/22
!!
!! get scalingmode from the ISE_T class

IMPLICIT NONE

class(ISE_T), INTENT(INOUT)     :: self
character(3)                   :: out

out = trim(self%nml%scalingmode)

end function getscalingmode_

!--------------------------------------------------------------------------
subroutine setuseangles_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setuseangles_
!! author: MDG
!! version: 1.0
!! date: 12/08/22
!!
!! set useangles in the ISE_T class

IMPLICIT NONE

class(ISE_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%useangles = trim(inp)

end subroutine setuseangles_

!--------------------------------------------------------------------------
function getuseangles_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getuseangles_
!! author: MDG
!! version: 1.0
!! date: 12/08/22
!!
!! get useangles from the ISE_T class

IMPLICIT NONE

class(ISE_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%useangles)

end function getuseangles_

!--------------------------------------------------------------------------
subroutine setmasterfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setmasterfile_
!! author: MDG
!! version: 1.0
!! date: 12/08/22
!!
!! set masterfile in the ISE_T class

IMPLICIT NONE

class(ISE_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%masterfile = trim(inp)

end subroutine setmasterfile_

!--------------------------------------------------------------------------
function getmasterfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getmasterfile_
!! author: MDG
!! version: 1.0
!! date: 12/08/22
!!
!! get masterfile from the ISE_T class

IMPLICIT NONE

class(ISE_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%masterfile)

end function getmasterfile_

!--------------------------------------------------------------------------
subroutine setISEMPfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setISEMPfile_
!! author: MDG
!! version: 1.0
!! date: 12/08/22
!!
!! set ISEMPfile in the ISE_T class

IMPLICIT NONE

class(ISE_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%ISEMPfile = trim(inp)

end subroutine setISEMPfile_

!--------------------------------------------------------------------------
function getISEMPfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getISEMPfile_
!! author: MDG
!! version: 1.0
!! date: 12/08/22
!!
!! get ISEMPfile from the ISE_T class

IMPLICIT NONE

class(ISE_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%ISEMPfile)

end function getISEMPfile_

!--------------------------------------------------------------------------
subroutine setdatafile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdatafile_
!! author: MDG
!! version: 1.0
!! date: 12/08/22
!!
!! set datafile in the ISE_T class

IMPLICIT NONE

class(ISE_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%datafile = trim(inp)

end subroutine setdatafile_

!--------------------------------------------------------------------------
function getdatafile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdatafile_
!! author: MDG
!! version: 1.0
!! date: 12/08/22
!!
!! get datafile from the ISE_T class

IMPLICIT NONE

class(ISE_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%datafile)

end function getdatafile_

!--------------------------------------------------------------------------
subroutine setimagefile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setimagefile_
!! author: MDG
!! version: 1.0
!! date: 12/08/22
!!
!! set imagefile in the ISE_T class

IMPLICIT NONE

class(ISE_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)       :: inp

self%nml%imagefile = trim(inp)

end subroutine setimagefile_

!--------------------------------------------------------------------------
function getimagefile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getimagefile_
!! author: MDG
!! version: 1.0
!! date: 12/08/22
!!
!! get imagefile from the ISE_T class

IMPLICIT NONE

class(ISE_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = trim(self%nml%imagefile)

end function getimagefile_

!--------------------------------------------------------------------------
recursive subroutine readISEMPfile_(self, HDF, HDFnames, nml, getmLPNH, getmLPSH, getmasterSPNH, &
                                    getmasterSPSH, getstrings, silent)
!DEC$ ATTRIBUTES DLLEXPORT :: readISEMPfile_
!! author: MDG 
!! version: 1.0 
!! date: 12/08/22
!!
!! read ISE master pattern file 

use mod_HDFsupport
use HDF5
use mod_HDFnames
use stringconstants 
use mod_io

use ISO_C_BINDING

IMPLICIT NONE

class(ISE_T), INTENT(INOUT)             :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames
type(ISEmasterNameListType)             :: nml
logical,INTENT(IN),OPTIONAL             :: getmLPNH
logical,INTENT(IN),OPTIONAL             :: getmLPSH
logical,INTENT(IN),OPTIONAL             :: getmasterSPNH
logical,INTENT(IN),OPTIONAL             :: getmasterSPSH
logical,INTENT(IN),OPTIONAL             :: getstrings
logical,INTENT(IN),OPTIONAL             :: silent

type(IO_T)                                       :: Message

character(fnlen)                                 :: infile, groupname, datagroupname, dataset
logical                                          :: stat, readonly, g_exists, f_exists
real(kind=sgl),allocatable                       :: mLPNH(:,:)
integer(HSIZE_T)                                 :: dims2(2)
integer(kind=irg)                                :: hdferr, istat, nlines
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)

associate( nml => self%nml )

! we assume that the calling program has opened the HDF interface
inquire(file=trim(self%getISEMPfile()), exist=f_exists)

if (.not.f_exists) then
  call Message%printError('readISEMPfile','ISE Master Pattern input file does not exist '//trim(self%ISEMPfile))
end if

! is this a proper HDF5 file ?
call h5fis_hdf5_f(trim(self%ISEMPfile), stat, hdferr)

if (stat.eqv..FALSE.) then ! the file exists, so let's open it an first make sure it is an EBSD dot product file
   call Message%printError('readISEMPfile','This is not a proper HDF5 file')
end if

! open the Master Pattern file
readonly = .TRUE.
hdferr =  HDF%openFile(self%ISEMPfile, readonly)

hdferr = HDF%openGroup(HDFnames%get_EMheader())
datagroupname = trim(HDFnames%get_ProgramData())
call H5Lexists_f(HDF%getobjectID(),trim(datagroupname),g_exists, hdferr)
if (.not.g_exists) then
  call Message%printError('readISEMPfile','This HDF file does not contain Master Pattern header data')
end if
call HDF%pop()

!====================================
! make sure this is a Master Pattern file
!====================================
hdferr = HDF%openGroup(HDFnames%get_NMLfiles())
dataset = trim(HDFnames%get_NMLfilename())
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists.eqv..FALSE.) then
    call HDF%pop(.TRUE.)
    call Message%printError('readISEMPfile','this is not a valid ISE Master Pattern file')
end if
call HDF%pop()

! get the name list strings if requested
if (present(getstrings)) then
  if (getstrings.eqv..TRUE.) then
    hdferr = HDF%openGroup(HDFnames%get_NMLfiles())
    dataset = trim(HDFnames%get_NMLfilename())
    call HDF%readdatasetstringarray(dataset, nlines, hdferr, self%nmlstrings)
    call HDF%pop()
  end if
end if

!====================================
! read all NMLparameters group datasets
!====================================
hdferr = HDF%openGroup(HDFnames%get_NMLparameters())
hdferr = HDF%openGroup(HDFnames%get_NMLlist())

! we need to first read all the parameters that are common to EBSD, ECP, and TKD namelists
dataset = 'a'
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetFloat(dataset, hdferr, self%mpnml%iscale(1))
end if

dataset = 'b'
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetFloat(dataset, hdferr, self%mpnml%iscale(2))
end if

dataset = 'c'
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetFloat(dataset, hdferr, self%mpnml%iscale(3))
end if

dataset = 'multiplier'
call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
if (g_exists.eqv..TRUE.) then
    call HDF%readDatasetFloat(dataset, hdferr, self%mpnml%multiplier)
end if

dataset = SC_npx
    call HDF%readDatasetInteger(dataset, hdferr, self%mpnml%npx)

dataset = SC_nthreads
    call HDF%readDatasetInteger(dataset, hdferr, self%mpnml%nthreads)

dataset = SC_outname
    call HDF%readDatasetStringArray(dataset, nlines, hdferr, stringarray)
    self%mpnml%outname = trim(stringarray(1))
    deallocate(stringarray)

dataset = SC_tiffname
    call HDF%readDatasetStringArray(dataset, nlines, hdferr, stringarray)
    self%mpnml%tiffname = trim(stringarray(1))
    deallocate(stringarray)

dataset = SC_xtalname
    call HDF%readDatasetStringArray(dataset, nlines, hdferr, stringarray)
    self%mpnml%xtalname = trim(stringarray(1))
    deallocate(stringarray)

! and close the NMLparameters group
call HDF%pop()
call HDF%pop()

! open the Master Pattern data group
hdferr = HDF%openGroup(HDFnames%get_EMData())
hdferr = HDF%openGroup(HDFnames%get_ProgramData())

if (present(getmLPNH)) then
  if (getmLPNH.eqv..TRUE.) then
    dataset = SC_mLPNH
    call HDF%readDatasetFloatArray(dataset, dims2, hdferr, mLPNH)
    allocate(self%det%mLPNH(-self%mpnml%npx:self%mpnml%npx,-self%mpnml%npx:self%mpnml%npx),stat=istat)
    self%det%mLPNH = mLPNH
    deallocate(mLPNH)
  end if
end if

if (present(getmLPSH)) then
  if (getmLPSH.eqv..TRUE.) then
    dataset = SC_mLPSH
    call HDF%readDatasetFloatArray(dataset, dims2, hdferr, mLPNH)
    allocate(self%det%mLPSH(-self%mpnml%npx:self%mpnml%npx,-self%mpnml%npx:self%mpnml%npx),stat=istat)
    self%det%mLPSH = mLPNH
    deallocate(mLPNH)
  end if
end if

if (present(getmasterSPNH)) then
  if (getmasterSPNH.eqv..TRUE.) then
    dataset = SC_masterSPNH
    call HDF%readDatasetFloatArray(dataset, dims2, hdferr, mLPNH)
    allocate(self%det%masterSPNH(-self%mpnml%npx:self%mpnml%npx,-self%mpnml%npx:self%mpnml%npx),stat=istat)
    self%det%masterSPNH = mLPNH
    deallocate(mLPNH)
  end if
end if

if (present(getmasterSPSH)) then
  if (getmasterSPSH.eqv..TRUE.) then
    dataset = SC_masterSPSH
    call HDF%readDatasetFloatArray(dataset, dims2, hdferr, mLPNH)
    allocate(self%det%masterSPSH(-self%mpnml%npx:self%mpnml%npx,-self%mpnml%npx:self%mpnml%npx),stat=istat)
    self%det%masterSPSH = mLPNH
    deallocate(mLPNH)
  end if
end if

! and close the HDF5 Master Pattern file
call HDF%pop(.TRUE.)

if (.not.present(silent)) then
  call Message%printMessage(' --> Completed reading master pattern data from '//trim(self%ISEMPfile), frm = "(A/)")
end if

end associate

end subroutine readISEMPfile_

!--------------------------------------------------------------------------
subroutine ComputeISEimage_(self, EMsoft, numang, Eangles)
!DEC$ ATTRIBUTES DLLEXPORT :: ComputeISEimage_
!! author: MDG 
!! version: 1.0 
!! date: 12/06/22
!!
!! compute an ISE image for a given ROI

use mod_EMsoft
use mod_io
use mod_quaternions 
use mod_rotations
use mod_Lambert
use mod_image 
use omp_lib
use ISO_C_BINDING
use, intrinsic :: iso_fortran_env

IMPLICIT NONE

class(ISE_T), INTENT(INOUT)                 :: self
type(EMsoft_T),INTENT(INOUT)                :: EMsoft
integer(kind=irg),INTENT(IN)                :: numang 
real(kind=sgl),INTENT(IN)                   :: Eangles(3,numang)

type(IO_T)                                  :: Message
type(Quaternion_T)                          :: quat, dquat 
type(e_T)                                   :: eu
type(q_T)                                   :: qu
type(Lambert_T)                             :: L

real(kind=sgl),allocatable                  :: ISEimage(:,:)
real(kind=sgl)                              :: s, mi, ma, ECPfactor, q(4)
integer(kind=irg)                           :: ix, iy, icnt, jd, sz(3), nxmc
real(kind=sgl)                              :: dc(3), avdc(3), newavdc(3), ixy(2), scl, sclmc, io_real(2)
real(kind=dbl)                              :: ddc(3)
real(kind=sgl)                              :: dx, dy, dxm, dym, x, y, z
integer(kind=irg)                           :: ii, jj, kk, istat
integer(kind=irg)                           :: nix, niy, nixp, niyp, nixmc, niymc, TID

! declare variables for use in object oriented image module
character(fnlen)                            :: TIFF_filename, fname 
integer                                     :: iostat
character(len=128)                          :: iomsg
logical                                     :: isInteger
type(image_t)                               :: im
integer(int8), allocatable                  :: TIFF_image(:,:)
integer                                     :: dim2(2)
integer(c_int32_t)                          :: result

call setRotationPrecision('d')

associate( enl => self%nml, mpnml => self%mpnml, ISE => self%det )

allocate(ISEimage(ISE%ipf_wd,ISE%ipf_ht))
ISEimage = 0.0
scl = float(self%mpnml%npx) 

! THIS ALL NEEDS TO BE REVIEWED AND IMPROVED !!!

! loop over all the image pixels 
icnt = 0
call OMP_SET_NUM_THREADS(enl%nthreads)

!$OMP PARALLEL default(shared) private(ix,iy,s,icnt,qu,jd,kk,dc,nix,niy,nixp,niyp,dx,dy,dxm,dym,newavdc,ECPfactor)&
!$OMP& private(nixmc, niymc, eu, quat, dquat, ddc) 

TID = OMP_GET_THREAD_NUM()

!$OMP DO SCHEDULE(DYNAMIC)
do iy = 1, ISE%ipf_ht
    do ix = 1, ISE%ipf_wd
        s = 0.0
        icnt = ISE%ipf_wd*(iy-1) + ix
! get the orientation and determine the quaternion to be applied to all the 
! detector vectors
        eu = e_T( edinp = dble( Eangles(1:3, icnt) ) )
        qu = eu%eq()
        quat = Quaternion_T( qd = qu%q_copyd() )

! get the pixel direction cosines 
        ddc = (/ 0.D0, 0.D0, 1.D0 /)
! apply the grain rotation to the detector direction cosines
        ddc = quat%quat_Lp( ddc )
        dc = sngl(ddc/sqrt(sum(ddc*ddc)))
! convert these direction cosines to interpolation coordinates in the Rosca-Lambert projection
        dc = dc / sqrt(sum(dc*dc))
        call LambertgetInterpolation(dc, scl, mpnml%npx, mpnml%npx, nix, niy, nixp, niyp, dx, dy, dxm, dym)

        if (dc(3) .ge. 0.0) then
                s = ISE%mLPNH(nix,niy) * dxm * dym + &
                    ISE%mLPNH(nixp,niy) * dx * dym + &
                    ISE%mLPNH(nix,niyp) * dxm * dy + &
                    ISE%mLPNH(nixp,niyp) * dx * dy 
          else
                s = ISE%mLPSH(nix,niy) * dxm * dym + &
                    ISE%mLPSH(nixp,niy) * dx * dym + &
                    ISE%mLPSH(nix,niyp) * dxm * dy + &
                    ISE%mLPSH(nixp,niyp) * dx * dy 
        end if
        ISEimage(ix,iy) = s 
    end do 
    if (mod(iy,10).eq.0) write (*,*) ' working on line ', iy
end do 
!$OMP END DO
!$OMP END PARALLEL

! and save the resulting ISE image to a tiff file
! output the ADP map as a tiff file 
fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%imagefile))
TIFF_filename = trim(fname)

! allocate memory for image
allocate(TIFF_image(ISE%ipf_wd,ISE%ipf_ht))

! fill the image with whatever data you have (between 0 and 255)
ma = maxval(ISEimage)
mi = minval(ISEimage)
io_real(1:2) = (/ mi, ma /)
call Message%WriteValue('Intensity range : ', io_real, 2) 

TIFF_image = int(255 * (ISEimage-mi)/(ma-mi))

! set up the image_t structure
im = image_t(TIFF_image)
if(im%empty()) call Message%printMessage("EMISE","failed to convert array to image")

! create the file
call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage("failed to write image to file : "//iomsg)
else  
  call Message%printMessage('ISE image written to '//trim(TIFF_filename))
end if 
deallocate(TIFF_image)

end associate 

end subroutine ComputeISEimage_

!--------------------------------------------------------------------------
subroutine ISE_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: ISE_
!! author: MDG 
!! version: 1.0 
!! date: 12/08/22
!!
!! perform the computations

use mod_EMsoft
use mod_so3
use mod_quaternions
use mod_DIfiles
use mod_ISEmaster
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_io
use mod_rotations
use stringconstants
use mod_vendors
use mod_memory

IMPLICIT NONE 

class(ISE_T), INTENT(INOUT)         :: self
type(EMsoft_T), INTENT(INOUT)       :: EMsoft
character(fnlen), INTENT(INOUT)     :: progname 

type(HDF_T)                         :: HDF
type(HDFnames_T)                    :: HDFnames
type(so3_T)                         :: SO
type(IO_T)                          :: Message
type(Quaternion_T)                  :: quat
type(QuaternionArray_T)             :: qAR
type(memory_T)                      :: mem
type(DIfile_T)                      :: DIFT
type(DictionaryIndexingNameListType):: dinl
type(Vendor_T)                      :: VT

type(ISEmasterNameListType)         :: mpnml

integer(kind=irg)                   :: i, sz(3), nx, hdferr, resang, resctf
character(fnlen)                    :: fname, DIfile
logical                             :: refined
real(kind=sgl),allocatable          :: Eangles(:,:)

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()
HDFnames = HDFnames_T()

associate( enl => self%nml, ISEdetector => self%det, DIdata => DIFT%DIDT )

! 1. read ISE master pattern file (HDF format)
call HDFnames%set_ProgramData(SC_ISEmaster)
call HDFnames%set_NMLlist(SC_ISEmasterNameList)
call HDFnames%set_NMLfilename(SC_ISEmasterNML)
! call HDFnames%set_Variable(SC_MCOpenCL)
fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfile))
call self%setISEMPfile(fname)
call self%readISEMPfile(HDF, HDFnames, mpnml, getmLPNH=.TRUE., getmLPSH=.TRUE.)

! 2. read the angular arrays from the HDF5 file (DI only for now)
resang = index(enl%datafile, '.ang')
resctf = index(enl%datafile, '.ctf')
DIfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(enl%datafile)
if ( (resang.eq.0).and.(resctf.eq.0) ) then 
  refined = .FALSE.
  call HDFnames%set_NMLfiles(SC_NMLfiles)
  call HDFnames%set_NMLfilename(SC_DictionaryIndexingNML)
  call HDFnames%set_NMLparameters(SC_NMLparameters)
  call HDFnames%set_NMLlist(SC_DictionaryIndexingNameListType)
  if (trim(enl%useangles).eq.'refined') then
      call DIFT%readDotProductFile(EMsoft, HDF, HDFnames, DIfile, hdferr, &
                                   getRefinedEulerAngles=.TRUE.)
      refined = .TRUE.
  else
      call DIFT%readDotProductFile(EMsoft, HDF, HDFnames, DIfile, hdferr, &
                                   getPhi1=.TRUE., &
                                   getPhi=.TRUE., &
                                   getPhi2=.TRUE.)
  end if

  dinl = DIFT%getNameList()

  if (sum(dinl%ROI).eq.0) then 
      ISEdetector%ipf_wd = dinl%ipf_wd
      ISEdetector%ipf_ht = dinl%ipf_ht
  else 
      ISEdetector%ipf_wd = dinl%ROI(3)
      ISEdetector%ipf_ht = dinl%ROI(4)
  end if

  nx = ISEdetector%ipf_wd * ISEdetector%ipf_ht
  allocate(Eangles(3, nx))
  if (trim(enl%useangles).eq.'original') then
      do i=1,nx 
        Eangles(1:3,i) = (/ DIdata%Phi1(i), DIdata%Phi(i), DIdata%Phi2(i) /)
      end do
      deallocate(DIdata%Phi1, DIdata%Phi, DIdata%Phi2)
  else 
      Eangles = DIdata%RefinedEulerAngles
      deallocate(DIdata%RefinedEulerAngles)
  end if

  if (maxval(Eangles).gt.(2.D0*cPi)) Eangles = Eangles * dtor
else
  VT = Vendor_T()
  if (resang.ne.0) then ! we have an .ang file
    call VT%getAnglesfromANGfile(DIfile, dinl%ipf_wd, dinl%ipf_ht, dinl%StepX, dinl%StepY, Eangles)
  else  ! we must have a .ctf file
    call VT%getAnglesfromCTFfile(DIfile, dinl%ipf_wd, dinl%ipf_ht, dinl%StepX, dinl%StepY, Eangles)
  end if 
  ISEdetector%ipf_wd = dinl%ipf_wd
  ISEdetector%ipf_ht = dinl%ipf_ht
  dinl%ROI = (/ 0, 0, 0, 0 /)
end if 

! 3. and finally perform the image computations
call self%ComputeISEimage(EMsoft, nx, Eangles)

end associate

! open the HDF interface
call closeFortranHDFInterface()

end subroutine ISE_



end module mod_ISE