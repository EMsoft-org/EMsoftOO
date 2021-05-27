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

module mod_reflectors
  !! author: MDG
  !! version: 1.0
  !! date: 03/24/20
  !!
  !! class definition for the EMreflectors program

use mod_kinds
use mod_global

IMPLICIT NONE

! namelist for the EMreflectors program
type, public :: reflectorsNameListType
  real(kind=sgl)          :: increment
  real(kind=sgl)          :: dmin
  integer(kind=irg)       :: numlist
  integer(kind=irg)       :: nthreads
  character(fnlen)        :: outputformat
  character(fnlen)        :: masterfile
  character(fnlen)        :: listfile
  logical                 :: kinematical
end type reflectorsNameListType

! class definition
type, public :: reflectors_T
private
  character(fnlen)              :: nmldeffile = 'EMreflectors.nml'
  type(reflectorsNameListType)  :: nml

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: reflectors_
  procedure, pass(self) :: get_increment_
  procedure, pass(self) :: get_dmin_
  procedure, pass(self) :: get_numlist_
  procedure, pass(self) :: get_nthreads_
  procedure, pass(self) :: get_outputformat_
  procedure, pass(self) :: get_masterfile_
  procedure, pass(self) :: get_listfile_
  procedure, pass(self) :: get_kinematical_
  procedure, pass(self) :: set_increment_
  procedure, pass(self) :: set_dmin_
  procedure, pass(self) :: set_numlist_
  procedure, pass(self) :: set_nthreads_
  procedure, pass(self) :: set_outputformat_
  procedure, pass(self) :: set_masterfile_
  procedure, pass(self) :: set_listfile_
  procedure, pass(self) :: set_kinematical_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: reflectors => reflectors_
  generic, public :: get_increment => get_increment_
  generic, public :: get_dmin => get_dmin_
  generic, public :: get_numlist => get_numlist_
  generic, public :: get_nthreads => get_nthreads_
  generic, public :: get_outputformat => get_outputformat_
  generic, public :: get_masterfile => get_masterfile_
  generic, public :: get_listfile => get_listfile_
  generic, public :: get_kinematical => get_kinematical_
  generic, public :: set_increment => set_increment_
  generic, public :: set_dmin => set_dmin_
  generic, public :: set_numlist => set_numlist_
  generic, public :: set_nthreads => set_nthreads_
  generic, public :: set_outputformat => set_outputformat_
  generic, public :: set_masterfile => set_masterfile_
  generic, public :: set_listfile => set_listfile_
  generic, public :: set_kinematical => set_kinematical_
end type reflectors_T

! the constructor routine for this class
interface reflectors_T
  module procedure reflectors_constructor
end interface reflectors_T

contains

!--------------------------------------------------------------------------
type(reflectors_T) function reflectors_constructor( nmlfile ) result(reflectors)
!DEC$ ATTRIBUTES DLLEXPORT :: reflectors_constructor
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! constructor for the reflectors_T Class; reads the name list

IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile

call reflectors%readNameList(nmlfile)

end function reflectors_constructor

!--------------------------------------------------------------------------
subroutine reflectors_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: reflectors_destructor
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! destructor for the reflectors_T Class

IMPLICIT NONE

type(reflectors_T), INTENT(INOUT)  :: self

call reportDestructor('reflectors_T')

end subroutine reflectors_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! read the namelist from an nml file for the reflectors_T Class

use mod_io
use mod_EMsoft

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)   :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft
type(IO_T)                           :: Message
logical                              :: skipread = .FALSE.

real(kind=sgl)                       :: increment
real(kind=sgl)                       :: dmin
integer(kind=irg)                    :: numlist
integer(kind=irg)                    :: nthreads
character(fnlen)                     :: outputformat
character(fnlen)                     :: masterfile
character(fnlen)                     :: listfile
logical                              :: kinematical

! define the IO namelist to facilitate passing variables to the program.
namelist /reflectors/ increment, dmin, masterfile, listfile, numlist, nthreads, outputformat, kinematical

! set the input parameters to default values (except for xtalname, which must be present)
increment = 0.025               ! angular increment [°]
dmin = 0.05                     ! smallest d-spacing to include in dynamical matrix [nm]
numlist = 20
nthreads = 1
outputformat = 'csv'            ! options: 'latex', 'csv', and 'markdown'
masterfile = 'undefined'        ! master pattern filename (EBSD/ECP/TKD)
listfile = 'undefined'          ! filename for output (no extension)
kinematical = .FALSE.           ! if .TRUE., a kinematical master pattern will be generated

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
  open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
  read(UNIT=dataunit,NML=reflectors)
  close(UNIT=dataunit,STATUS='keep')

! check for required entries
  if (trim(listfile).eq.'undefined') then
    call Message%printError('readNameList:',' output file name is undefined in '//nmlfile)
  end if
  if (trim(masterfile).eq.'undefined') then
    call Message%printError('readNameList:',' master file name is undefined in '//nmlfile)
  end if
end if

! if we get here, then all appears to be ok, and we need to fill in the emnl fields
self%nml%increment = increment
self%nml%dmin = dmin
self%nml%outputformat = trim(outputformat)
self%nml%numlist = numlist
self%nml%nthreads = nthreads
self%nml%masterfile = trim(masterfile)
self%nml%listfile = trim(listfile)
self%nml%kinematical = kinematical

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! pass the namelist for the reflectors_T Class to the calling program

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)   :: self
type(reflectorsNameListType)         :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
function get_increment_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_increment_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! get increment from the reflectors_T class

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)     :: self
real(kind=sgl)                         :: out

out = self%nml%increment

end function get_increment_

!--------------------------------------------------------------------------
subroutine set_increment_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_increment_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! set increment in the reflectors_T class

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)             :: inp

self%nml%increment = inp

end subroutine set_increment_

!--------------------------------------------------------------------------
function get_dmin_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dmin_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! get dmin from the reflectors_T class

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)     :: self
real(kind=sgl)                         :: out

out = self%nml%dmin

end function get_dmin_

!--------------------------------------------------------------------------
subroutine set_dmin_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_dmin_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! set dmin in the reflectors_T class

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)             :: inp

self%nml%dmin = inp

end subroutine set_dmin_

!--------------------------------------------------------------------------
function get_numlist_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_numlist_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! get numlist from the reflectors_T class

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)     :: self
integer(kind=irg)                      :: out

out = self%nml%numlist

end function get_numlist_

!--------------------------------------------------------------------------
subroutine set_numlist_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_numlist_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! set numlist in the reflectors_T class

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)          :: inp

self%nml%numlist = inp

end subroutine set_numlist_

!--------------------------------------------------------------------------
function get_nthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nthreads_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! get nthreads from the reflectors_T class

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)     :: self
integer(kind=irg)                      :: out

out = self%nml%nthreads

end function get_nthreads_

!--------------------------------------------------------------------------
subroutine set_nthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nthreads_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! set nthreads in the reflectors_T class

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)          :: inp

self%nml%nthreads = inp

end subroutine set_nthreads_

!--------------------------------------------------------------------------
function get_outputformat_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_outputformat_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! get outputformat from the reflectors_T class

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)     :: self
character(fnlen)                       :: out

out = self%nml%outputformat

end function get_outputformat_

!--------------------------------------------------------------------------
subroutine set_outputformat_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_outputformat_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! set outputformat in the reflectors_T class

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)           :: inp

self%nml%outputformat = inp

end subroutine set_outputformat_

!--------------------------------------------------------------------------
function get_masterfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_masterfile_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! get masterfile from the reflectors_T class

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)     :: self
character(fnlen)                       :: out

out = self%nml%masterfile

end function get_masterfile_

!--------------------------------------------------------------------------
subroutine set_masterfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_masterfile_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! set masterfile in the reflectors_T class

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)           :: inp

self%nml%masterfile = inp

end subroutine set_masterfile_

!--------------------------------------------------------------------------
function get_listfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_listfile_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! get listfile from the reflectors_T class

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)     :: self
character(fnlen)                       :: out

out = self%nml%listfile

end function get_listfile_

!--------------------------------------------------------------------------
subroutine set_listfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_listfile_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! set listfile in the reflectors_T class

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)           :: inp

self%nml%listfile = inp

end subroutine set_listfile_

!--------------------------------------------------------------------------
function get_kinematical_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_kinematical_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! get kinematical from the reflectors_T class

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)     :: self
logical                                :: out

out = self%nml%kinematical

end function get_kinematical_

!--------------------------------------------------------------------------
subroutine set_kinematical_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_kinematical_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! set kinematical in the reflectors_T class

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)     :: self
logical, INTENT(IN)                    :: inp

self%nml%kinematical = inp

end subroutine set_kinematical_


!--------------------------------------------------------------------------
subroutine reflectors_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: reflectors_
!! author: MDG
!! version: 1.0
!! date: 03/24/20
!!
!! perform the computations

use mod_EMsoft
use mod_initializers
use mod_crystallography
use mod_diffraction
use mod_symmetry
use mod_io
use mod_rotations
use mod_quaternions
use mod_Lambert
use mod_MCfiles
use mod_MPfiles
use mod_timing
use mod_Lambert
use clfortran
! use mod_CLsupport
use HDF5
use mod_HDFsupport
use mod_HDFnames
use ISO_C_BINDING
use stringconstants
use omp_lib
use mod_OMPsupport

IMPLICIT NONE

class(reflectors_T), INTENT(INOUT)      :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname

type(cell_T)                            :: cell
type(HDF_T)                             :: HDF
type(HDFnames_T)                        :: HDFnames
type(MCfile_T)                          :: MCFT
type(MPfile_T)                          :: MPFT
type(Lambert_T)                         :: Lambert
type(IO_T)                              :: Message
type(Diffraction_T)                     :: Diff
type(SpaceGroup_T)                      :: SG
type(Quaternion_T)                      :: qu
type(q_T)                               :: qq
type(a_T)                               :: ax

type(MCOpenCLNameListType)              :: mcnl
type(EBSDmasterNameListType)            :: mpnlEBSD
type(ECPmasterNameListType)             :: mpnlECP
type(TKDmasterNameListType)             :: mpnlTKD

character(fnlen)                        :: listfile, masterfile, groupname, dataset, xtalname, outputfile, infile, fname
logical                                 :: f_exists, readonly, verbose
integer(kind=irg)                       :: hdferr, nlines, i, istat, ix, iy, nx, io_int(1), nkeep, nl2, k2
integer(HSIZE_T)                        :: dims3(3), dims4(3)
real(kind=dbl)                          :: EkeV, mLambda
real(kind=sgl)                          :: m

integer(kind=irg)                       :: imh, imk, iml, ii, j, num, nums, mhkl, valpos, numphi,numtheta,iequiv(3,48),nequiv
integer(kind=irg),allocatable           :: family(:,:,:),numfam(:),idx(:), idx2(:), sfi(:), itmp(:,:)
integer(kind=irg)                       :: h,k,l,totfam,ind(3),icnt, oi_int(1), g1(3), g2(3), NUMTHREADS, TID
logical                                 :: first
logical,allocatable                     :: z(:,:,:)
real(kind=sgl)                          :: g(3), thr, dphi, gc(3), gax(3), gz(3), v(3), x, val1, val2, valmax
real(kind=sgl),allocatable              :: Vgg(:),ddg(:),gg(:),th(:), gcart(:,:), gcrys(:,:), cp(:), sp(:), ca(:), sa(:), &
                                           phi(:), theta(:), dc(:,:), Vg(:), VgX(:), VggX(:), cosnorm(:)
character(1)                            :: space
real(kind=sgl)                          :: dhkl, incrad, glen, sd
real(kind=sgl)                          :: ixy(2),scl, dmin
real(kind=sgl)                          :: dx,dy,dxm,dym
integer(kind=irg)                       :: jj,kk
integer(kind=irg)                       :: nix,niy,nixp,niyp
logical,allocatable                     :: keep(:)

integer(kind=irg),allocatable           :: acc_e(:,:,:)
real(kind=sgl),allocatable              :: Eweights(:)
real(kind=sgl),allocatable              :: srtmp(:,:,:,:), mLPNH(:,:,:), mLPSH(:,:,:), masterNH(:,:), masterSH(:,:), &
                                           KBI(:), kinmasterNH(:,:), kinmasterSH(:,:), kinNH(:,:), kinSH(:,:)
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)

type(DynType),save                      :: Dyn
type(gnode),save                        :: rlp

thr = 1.E-5
space = 'r'
call setRotationPrecision('s')

associate( enl=>self%nml )

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T()

! set the HDF group names for this program
HDFnames = HDFnames_T()

! 1. determine the master pattern file modality
fname = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfile))
call MPFT%determineModality(HDF, fname)
call Message%printMessage(' Input file modality is '//trim(MPFT%getModality()))

! 2. read the Monte Carlo data file (HDF format)
call HDFnames%set_ProgramData(SC_MCOpenCL)
call HDFnames%set_NMLlist(SC_MCCLNameList)
call HDFnames%set_NMLfilename(SC_MCOpenCLNML)
call MCFT%setFileName(fname)
call MCFT%readMCfile(HDF, HDFnames, getAccume=.TRUE.)
mcnl = MCFT%getnml()

!------------------------------
! compute the energy weight factors by integrating the lower rectangular portion
! of the Lambert projection; we'll take the lower quarter in vertical dimension
! and a similar distance to left and right in the horizontal direction for EBSD/TKD.
! For ECP we'll need a different integration area for the weights; to be implemented.
!------------------------------
dims3 = shape(MCFT%MCDT%accum_e)
allocate(Eweights(dims3(1)))
Eweights = 0.0
do i=1,dims3(1)
  Eweights(i) = sum(MCFT%MCDT%accum_e(i,-dims3(2)/4:dims3(2)/4,-dims3(3)/2:-dims3(3)/4))
end do
Eweights = Eweights/maxval(Eweights)
deallocate(MCFT%MCDT%accum_e)

! 3. read EBSD/ECP/TKD/Kossel master pattern file (HDF format)
call HDFnames%set_Variable(SC_MCOpenCL)
call MPFT%setFileName(fname)
if (trim(MPFT%getModality()).eq.'EBSD') then
  call HDFnames%set_ProgramData(SC_EBSDmaster)
  call HDFnames%set_NMLlist(SC_EBSDmasterNameList)
  call HDFnames%set_NMLfilename(SC_EBSDmasterNML)
  call MPFT%readMPfile(HDF, HDFnames, mpnlEBSD, getmLPNH=.TRUE., getmLPSH=.TRUE.)
  dmin = mpnlEBSD%dmin
  else if (trim(MPFT%getModality()).eq.'ECP') then
    call HDFnames%set_ProgramData(SC_ECPmaster)
    call HDFnames%set_NMLlist(SC_ECPmasterNameList)
    call HDFnames%set_NMLfilename(SC_ECPmasterNML)
    call MPFT%readMPfile(HDF, HDFnames, mpnlECP, getmLPNH=.TRUE., getmLPSH=.TRUE.)
    dmin = mpnlECP%dmin
    else if (trim(MPFT%getModality()).eq.'TKD') then
      call HDFnames%set_ProgramData(SC_TKDmaster)
      call HDFnames%set_NMLlist(SC_TKDmasterNameList)
      call HDFnames%set_NMLfilename(SC_TKDmasterNML)
      call MPFT%readMPfile(HDF, HDFnames, mpnlTKD, getmLPNH=.TRUE., getmLPSH=.TRUE.)
      dmin = mpnlTKD%dmin
      else if (trim(MPFT%getModality()).eq.'Kossel') then
        ! call HDFnames%set_ProgramData(SC_Kosselmaster)
        ! call HDFnames%set_NMLlist(SC_KosselmasterNameList)
        ! call HDFnames%set_NMLfilename(SC_KosselmasterNML)
        ! call MPFT%readMPfile(HDF, HDFnames, mpnl, getmLPNH=.TRUE., getmLPSH=.TRUE.)
      end if

dims4 = shape(MPFT%MPDT%mLPNH)
nx = (dims4(1)-1)/2

! perform E-weighted averaging to get a single NH+SH master pattern
allocate(masterNH(-nx:nx,-nx:nx),stat=istat)
allocate(masterSH(-nx:nx,-nx:nx),stat=istat)
masterNH = 0.0
masterSH = 0.0
do ix=-nx,nx
  do iy=-nx,nx
    masterNH(ix,iy) = sum(MPFT%MPDT%mLPNH(ix,iy,1:dims4(3))*Eweights(1:dims3(1)))
    masterSH(ix,iy) = sum(MPFT%MPDT%mLPSH(ix,iy,1:dims4(3))*Eweights(1:dims3(1)))
  end do
end do
deallocate(MPFT%MPDT%mLPNH, MPFT%MPDT%mLPSH)

! do we need to generate a kinematical master pattern ?
if (enl%kinematical.eqv..TRUE.) then
  allocate(kinmasterNH(-nx:nx, -nx:nx), kinmasterSH(-nx:nx,-nx:nx) )
  kinmasterNH = 0.0
  kinmasterSH = 0.0
end if

! subtract the average value from the master pattern arrays and divide by the standard deviation
m = sum(masterNH)/float((2*nx+1)**2)
masterNH = masterNH - m
sd = sqrt( sum(masterNH**2) / (float((2*nx+1)**2 - 1)))
masterNH = masterNH / sd

m = sum(masterSH)/float((2*nx+1)**2)
masterSH = masterSH - m
sd = sqrt( sum(masterSH**2) / (float((2*nx+1)**2 - 1)))
masterSH = masterSH / sd

! ok, now we have the averaged master pattern; next we need to init the crystal
! structure, get a list of unique reflectors, and for each one integrate the
! Kikuchi band...

! get the crystal structure from the *.xtal file
verbose = .TRUE.
mcnl = MCFT%getnml()
xtalname = trim(mcnl%xtalname)
call cell%setFileName(xtalname)
call Diff%setV(mcnl%EkeV)
call Diff%setrlpmethod('WK')
call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, dmin, verbose, useHDF=HDF)

! generate a list of hkl indices for which we need to compute the integral over the Kikuchi band
! since the commercial packages use the kinematical structure factor to generate this list, we
! will do the same...

! compute the range of reflections for the lookup table and allocate the table
! The master list is easily created by brute force
 imh = 1
 do
   dhkl = 1.0/cell%CalcLength((/float(imh) ,0.0_sgl,0.0_sgl/), 'r')
   if (dhkl.lt.dmin) EXIT
   imh = imh + 1
 end do
 imk = 1
 do
   dhkl = 1.0/cell%CalcLength((/0.0_sgl,float(imk),0.0_sgl/), 'r')
   if (dhkl.lt.dmin) EXIT
   imk = imk + 1
 end do
 iml = 1
 do
   dhkl = 1.0/cell%CalcLength((/0.0_sgl,0.0_sgl,float(iml)/), 'r')
   if (dhkl.lt.dmin) EXIT
   iml = iml + 1
 end do

! allocate all arrays
 allocate(z(-2*imh:2*imh,-2*imk:2*imk,-2*iml:2*iml))
 ii = (2*imh+1)*(2*imk+1)*(2*iml+1)
 allocate(family(ii,48,3))
 allocate(numfam(ii))
 allocate(Vgg(ii))
 allocate(VggX(ii))
 allocate(ddg(ii))
 allocate(gg(ii))
 allocate(th(ii))
! determine the families of reflections with (hkl)<=(imh,imk,iml)
! first initialize the boolean array z
 z = .FALSE.
! then loop through all (hkl) values
 first = .TRUE.
 icnt = 1
 totfam=0
 call Diff%setrlpmethod('DT') ! we're computing simple Doyle-Turner or Smith-Burge scattering factors to get the list of reflectors
 do h=-imh,imh
  ind(1)=-h
  do k=-imk,imk
   ind(2)=-k
   do l=-iml,iml
    ind(3)=-l

! make sure we have not already done this one in another family
    if (.not.z(-h,-k,-l)) then

! if it is a new one, then determine the entire family
     call Diff%setrlpmethod('DT')
     call Diff%CalcUcg(cell,ind)
     rlp = Diff%getrlp()

! but ignore the reciprocal lattice point if Vgg is small
     if (abs(rlp%Ucg).ge.thr) then

! copy family in array and label all its members and multiples in z-array
      call SG%CalcFamily(ind,num,space,itmp)
      do i=1,num
       do j=1,3
        family(icnt,i,j)=itmp(i,j)
       end do
       z(itmp(i,1),itmp(i,2),itmp(i,3))=.TRUE.
      end do

! compute the X-ray kinematical structure factor
      call Diff%setrlpmethod('XR')
      call Diff%CalcUcg(cell,ind)
      rlp = Diff%getrlp()
      VggX(icnt) = abs(rlp%Ucg)

! also get the structure factor with the WK parameters and absorption
      call Diff%setrlpmethod('WK')
      call Diff%CalcUcg(cell,ind)
      rlp = Diff%getrlp()
      Vgg(icnt) = rlp%Vmod

! increment family counter
      numfam(icnt)=num
      totfam=totfam+num-1
      icnt=icnt+1
     end if
    end if
   end do
  end do
 end do

 icnt=icnt-1
 oi_int(1)=icnt
 call Message%WriteValue(' Total number of families        = ', oi_int, 1, "(I6)")
 oi_int(1)=totfam
 call Message%WriteValue(' Total number of family members  = ', oi_int, 1, "(I6)")

! compute d-spacings, g-spacings, theta
 call Message%printMessage(' Computing d-spacings, g-spacings, and scattering angles')
 allocate(gcart(3,icnt),gcrys(3,icnt))
 mLambda = Diff%getWaveLength()
 do k=1,icnt
  g(1:3)=float(family(k,1,1:3))
  gg(k)=cell%CalcLength(g,'r')
  gcrys(1:3,k) = g(1:3)
  call cell%TransSpace(g,gc,'r','c')
  call cell%NormVec(gc,'c')
  gcart(1:3,k) = gc(1:3)
  th(k)=asin(0.5*mLambda*gg(k))
 end do

! here we need to eliminate those entries that are multiples of a smaller hkl
! and only keep the one that has the largest structure factor.

allocate(idx(icnt))
call SPSORT(gg,icnt,idx,1,istat)

allocate(keep(icnt))
keep = .TRUE.
keep(idx(1)) = .FALSE.   ! eliminate (000) from the list

mhkl = int(maxval(gcrys))

call Message%printMessage(' Selecting lowest hkl values with largest structure factor ')
do k=2,icnt-1
 if (keep(idx(k)).eqv..TRUE.) then
  valpos = idx(k)
  g1 = int(gcrys(:,valpos))
  val1 = VggX(valpos)
  valmax = val1
  keep(valpos) = .TRUE.
! scan through the multiples
  do j=2,mhkl
    g2 = j * g1
    do i=k+1,icnt
      if (sum(abs(g2-int(gcrys(:,idx(i))))).eq.0) then
        if (VggX(idx(i)).gt.valmax) then
          keep(valpos) = .FALSE.
          valpos = idx(i)
          valmax = VggX(valpos)
          keep(valpos) = .TRUE.
        else
          keep(idx(i)) = .FALSE.
        end if
      end if
    end do
  end do
 end if
end do

nkeep = 0
do i=1,icnt
  if (keep(i).eqv..TRUE.) nkeep = nkeep+1
end do

!=======================================
!=======================================
!=======================================
! and here is the main part of this program: Kikuchi band integration for
! each unique family (one member per family).

! allocate the output arrays (KikuchiBandIntegral = KBI)
allocate(KBI(icnt),Vg(icnt),VgX(icnt))
KBI = 0.0

! azimuthal integration angle
numphi = 360.0/enl%increment
allocate(phi(numphi), cp(numphi),sp(numphi))
dphi = enl%increment * sngl(cPi)/180.0
phi = (/ (float(i-1)*dphi,i=1,numphi) /)
cp = cos(phi)
sp = sin(phi)

incrad = enl%increment * sngl(cPi)/180.0
scl = float(nx)

! set the number of OpenMP threads
call OMP_setNThreads(enl%nthreads)
io_int(1) = enl%nthreads
call Message%WriteValue(' Setting # threads to ',io_int,1,"(I3)")
io_int(1) = nkeep
call Message%WriteValue(' Total number of integrations to be carried out ',io_int,1,"(I6)")

if (enl%kinematical.eqv..TRUE.) then
  call Message%printMessage(' Computation of symmetrized kinematical pattern will slow things down a bit ... ')
end if

call Message%printMessage(' Starting parallel integrations... (.=100, |=1000) ')

! use OpenMP to run on multiple cores ...
!$OMP PARALLEL DEFAULT(PRIVATE) &
!$OMP& SHARED(k, nx, cp, sp, icnt, keep, th, incrad, numphi, gcart, cell, scl, masterNH, masterSH) &
!$OMP& SHARED(Vg, VgX, Vgg, VggX, KBI, nkeep, kinmasterNH, kinmasterSH)

NUMTHREADS = OMP_GET_NUM_THREADS()
TID = OMP_GET_THREAD_NUM()

allocate(kinNH(-nx:nx,-nx:nx), kinSH(-nx:nx,-nx:nx))
kinNH = 0.0
kinSH = 0.0

!$OMP DO SCHEDULE(STATIC,enl%nthreads)
do k=1,icnt-1   ! ignore the last point
 if (keep(k)) then
  ii = nint(th(k)/incrad)
  numtheta = 2*ii+1
  allocate(theta(numtheta),ca(numtheta),sa(numtheta))
  nums = numphi * numtheta
  allocate( dc(3,nums), cosnorm(nums) )


  theta = (/ (float(i),i=-ii,ii) /) * incrad
  gz = (/ 0.0, 0.0, 1.0 /)

! get the unrotated direction cosines of the sampling points on the sphere
! also initialize the cosine term in the surface integration, and the segment normalization
    ca = cos( theta )
    sa = sin( theta )
    ii = 1
    do i=1,numphi
      do j=1,numtheta
        dc(1,ii) = ca(j) * cp(i)
        dc(2,ii) = ca(j) * sp(i)
        dc(3,ii) = sa(j)
        cosnorm(ii) = ca(j) / ( 4.0 * sngl(cPi) * sin(th(k)) )
        ii = ii+1
      end do
    end do


! then determine the rotation quaternion to bring the z axis onto the g direction (cartesian)
    v = gcart(1:3,k)
    x = cell%CalcDot(gz,v,'c')
    if (x.ne.1.0) then   ! the cross product exists
      call cell%CalcCross(v,gz,gax,'c','c',0) ! gax is the rotation axis
      gax = gax / sqrt(sum(gax*gax))
      x = acos(x)
      if (x.lt.0.0) then
        ax = a_T( ainp = (/-gax(1),-gax(2),-gax(3), -x /) )
      else
        ax = a_T( ainp = (/ gax(1), gax(2), gax(3), x /) )
      end if
      qq = ax%aq()
      call qu%set_quats( qq%q_copy() )
      qu = conjg(qu)
      do i=1,nums
        v = qu%quat_Lp( dc(:,i) )
        v = v / sqrt(sum(v*v))
        dc(:,i) = v(:)
      end do
    end if

! now that all the points have been rotated, we simply transform the direction cosines
! into square Lambert coordinates and interpolate from the master patterns in the usual way...
    do i=1,nums
  ! convert these direction cosines to coordinates in the Rosca-Lambert projection
      v(:) = dc(:,i)
      call LambertgetInterpolation(v, scl, nx, nx, nix, niy, nixp, niyp, dx, dy, dxm, dym)

  ! interpolate the intensity
      if (dc(3,i) .ge. 0.0) then
         KBI(k) = KBI(k)+ ( masterNH(nix,niy) * dxm * dym +  masterNH(nixp,niy) * dx * dym + &
                            masterNH(nix,niyp) * dxm * dy +  masterNH(nixp,niyp) * dx * dy ) * 0.25 * cosnorm(i)
      else
         KBI(k) = KBI(k)+ ( masterSH(nix,niy) * dxm * dym +  masterSH(nixp,niy) * dx * dym + &
                            masterSH(nix,niyp) * dxm * dy +  masterSH(nixp,niyp) * dx * dy ) * 0.25 * cosnorm(i)
      end if
      Vg(k) = Vgg(k)
      VgX(k) = VggX(k)
      if (enl%kinematical.eqv..TRUE.) then ! add the kinematical intensity and symmetrize it
        call Lambert%Apply3DPGSymmetry(cell,SG,nix,niy,1,nx,iequiv,nequiv)
        do ix=1,nequiv
          if (iequiv(3,ix).eq.-1) then
            kinSH(iequiv(1,ix),iequiv(2,ix)) = kinSH(iequiv(1,ix),iequiv(2,ix)) + Vg(k)
          else
            kinNH(iequiv(1,ix),iequiv(2,ix)) = kinNH(iequiv(1,ix),iequiv(2,ix)) + Vg(k)
          end if
        end do
      end if
    end do
    deallocate(theta,ca,sa,cosnorm,dc)
 else
    Vg(k) = 0.0
    VgX(k) = 0.0
 end if
 if (mod(k,100).eq.0) then
  if (mod(k,1000).eq.0) then
     write (*,"('|')",advance="no")
   else
     write (*,"('.')",advance="no")
   end if
 end if
end do
!$OMP END DO

if (enl%kinematical.eqv..TRUE.) then
!$OMP CRITICAL
  kinmasterNH = kinmasterNH + kinNH
  kinmasterSH = kinmasterSH + kinSH
!$OMP END CRITICAL
end if

!$OMP END PARALLEL

call Message%printMessage(' done ')


x = maxval(KBI)
KBI = KBI * 100.0/x
Vg(icnt) = 0.0
x = maxval(Vg)
Vg = Vg * 100.0/x
VgX(icnt) = 0.0
x = maxval(VgX)
VgX = VgX * 100.0/x

allocate(idx2(icnt))
call SPSORT(Vg,icnt,idx2,-1,istat)
call SPSORT(KBI,icnt,idx,-1,istat)

allocate(sfi(icnt))
do i=1,icnt
  k = idx2(i)
  sfi(k) = i
end do

listfile = EMsoft%generateFilePath('EMdatapathname',trim(enl%listfile))

if ((trim(enl%outputformat).eq.'latex').or.(trim(enl%outputformat).eq.'all')) then
  outputfile = trim(listfile)//'.tex'
  open(unit=80,file=trim(outputfile),status='unknown',form='formatted')

! format everything as a LaTeX table, with rank, hkl, |g|, KBI, Vg (sfi)
  write (80,"('\begin{table}[th]\caption{reflector ranking}\centering\leavevmode\begin{tabular}{llrrrcllrrr}')")
  write (80,"('\hline $\#$ & $(hkl)$ & $\beta_{hkl}$ & $I^{\text{abs}}_{hkl}$ & $I^{\text{X}}_{hkl}$ & $\quad$ &')")
  write (80,"('$\#$ & $(hkl)$ & $\beta_{hkl}$ & $I^{\text{abs}}_{hkl}$ & $I^{\text{X}}_{hkl}$\\')")
  write (80,"('\hline')")
  if (mod(enl%numlist,2).eq.0) then 
    nl2 = enl%numlist/2
  else 
    nl2 = (enl%numlist+1)/2
  end if 
  do i=1,enl%numlist
    k = idx(i)
    if ((i+nl2).le.enl%numlist) then 
      k2 = idx(i+nl2)
    else
      k2 = -1
    end if
! the first reflection for this line in the table
    if ((sum(abs(gcrys(:,k))).ne.0.0).and.(KBI(k).ne.0.0)) then
      glen = cell%CalcLength(gcrys(:,k),'r')
        write (80,"(I2,'& $(')",advance="no") i
        do jj=1,3
          if (int(gcrys(jj,k)).lt.0) then
            if (jj.lt.3) then
              write (80,"('\bar{',I3,'}\,')",advance="no") abs(int(gcrys(jj,k)))
            else
              write (80,"('\bar{',I3,'}')",advance="no") abs(int(gcrys(jj,k)))
            end if
          else
            if (jj.lt.3) then
              write (80,"(I3,'\,')",advance="no") int(gcrys(jj,k))
            else
              write (80,"(I3)",advance="no") int(gcrys(jj,k))
            end if
          end if
        end do
        write (80,"(')$ & ',F6.2,' & ',F6.2,' & ',F6.2,' & &')") KBI(k), Vg(k), VgX(k)
    end if
! the second reflection for this line in the table
    if ((sum(abs(gcrys(:,k2))).ne.0.0).and.(KBI(k2).ne.0.0).and.(k2.ne.-1)) then
      glen = cell%CalcLength(gcrys(:,k2),'r')
        write (80,"(I2,'& $(')",advance="no") i+nl2
        do jj=1,3
          if (int(gcrys(jj,k2)).lt.0) then
            if (jj.lt.3) then
              write (80,"('\bar{',I3,'}\,')",advance="no") abs(int(gcrys(jj,k2)))
            else
              write (80,"('\bar{',I3,'}')",advance="no") abs(int(gcrys(jj,k2)))
            end if
          else
            if (jj.lt.3) then
              write (80,"(I3,'\,')",advance="no") int(gcrys(jj,k2))
            else
              write (80,"(I3)",advance="no") int(gcrys(jj,k2))
            end if
          end if
        end do
        write (80,"(')$ & ',F6.2,' & ',F6.2,' & ',F6.2,'\\ ')") KBI(k2), Vg(k2), VgX(k2)
      end if
      if (k2.eq.-1) write (80,"('\\')") 
  end do
  write(80,"('\hline\end{tabular}\end{table}')")
  close(unit=80,status='keep')
  call Message%printMessage('Data stored in LaTeX file '//trim(outputfile))
end if


if ((trim(enl%outputformat).eq.'csv').or.(trim(enl%outputformat).eq.'all')) then
  outputfile = trim(listfile)//'.csv'
  open(unit=80,file=trim(outputfile),status='unknown',form='formatted')

  write (80,"(A)") '#,h,k,l,KBI,Ikin+abs,IX'

  do i=1,enl%numlist
    k = idx(i)
    if ((sum(abs(gcrys(:,k))).ne.0.0).and.(KBI(k).ne.0.0)) then
      write (80,"(I4,',',3(I3,','),F6.2,',',F6.2,',',F6.2)") i,int(gcrys(:,k)),KBI(k),Vg(k),VgX(k)
    end if
  end do
  close(unit=80,status='keep')
  call Message%printMessage('Data stored in .csv file '//trim(outputfile))
end if

if ((trim(enl%outputformat).eq.'markdown').or.(trim(enl%outputformat).eq.'all')) then
  outputfile = trim(listfile)//'.md'
  open(unit=80,file=trim(outputfile),status='unknown',form='formatted')

  write (80,"(A)") '##EBSD Dynamical Reflector Ranking for '//trim(mcnl%xtalname)
  write (80,"(A)") ' '
  write (80,"(A)") '|\# | (hkl) | beta_hkl | Ikin+abs |    IX   |'
  write (80,"(A)") '|---|-------|----------|----------|---------|'

  do i=1,enl%numlist
    k = idx(i)
    if ((sum(abs(gcrys(:,k))).ne.0.0).and.(KBI(k).ne.0.0)) then
        write (80,"('|',I2,'| (')",advance="no") i
        do jj=1,3
          write (80,"(I3)",advance="no") int(gcrys(jj,k))
        end do
        write (80,"(')|',F6.2,'|',F6.2,'|',F6.2,'|')") KBI(k), Vg(k), VgX(k)
    end if
  end do
  close(unit=80,status='keep')
  call Message%printMessage('Data stored in .md file '//trim(outputfile))
end if

! do we need to store the kinematical patterns in the MP file ?
if (enl%kinematical.eqv..TRUE.) then

  infile = EMsoft%generateFilePath('EMdatapathname',trim(enl%masterfile))

  hdferr =  HDF%openFile(infile)

  groupname = SC_EMData
  hdferr = HDF%openGroup(groupname)
  groupname = SC_EBSDmaster
  hdferr = HDF%openGroup(groupname)

  dataset = 'kinmasterNH'
  hdferr = HDF%writeDatasetFloatArray(dataset, kinmasterNH, 2*nx+1, 2*nx+1)

  dataset = 'kinmasterSH'
  hdferr = HDF%writeDatasetFloatArray(dataset, kinmasterSH, 2*nx+1, 2*nx+1)

  call HDF%pop(.TRUE.)
end if

! close the HDF interface
call closeFortranHDFInterface()

end associate

end subroutine reflectors_



end module mod_reflectors
