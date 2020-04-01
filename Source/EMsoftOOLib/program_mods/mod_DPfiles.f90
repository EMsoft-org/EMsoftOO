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

module mod_DPfiles
  !! author: MDG 
  !! version: 1.0 
  !! date: 03/31/20
  !!
  !! class definition for Dot Product file routines; also contains
  !! everything dealing with dictionary indexing namelists.

use mod_kinds
use mod_global

IMPLICIT NONE 

! parent namelist for dictionary indexing programs; this is the minimum needed
! for EBSD indexing, but other modalities may need additional parameters that 
! are defined via inherited classes
type, public :: DictionaryIndexingNameListType
  integer(kind=irg)  :: ncubochoric
  integer(kind=irg)  :: numexptsingle
  integer(kind=irg)  :: numdictsingle
  integer(kind=irg)  :: ipf_ht
  integer(kind=irg)  :: ipf_wd
  integer(kind=irg)  :: ROI(4)
  integer(kind=irg)  :: nnk
  integer(kind=irg)  :: nnav
  integer(kind=irg)  :: nosm
  integer(kind=irg)  :: nism
  integer(kind=irg)  :: maskradius
  integer(kind=irg)  :: numsx
  integer(kind=irg)  :: numsy
  integer(kind=irg)  :: binning
  integer(kind=irg)  :: nthreads
  integer(kind=irg)  :: devid
  integer(kind=irg)  :: usenumd
  integer(kind=irg)  :: multidevid(8)
  integer(kind=irg)  :: platid
  integer(kind=irg)  :: nregions
  integer(kind=irg)  :: nlines
  real(kind=sgl)     :: L
  real(kind=sgl)     :: thetac
  real(kind=sgl)     :: delta
  real(kind=sgl)     :: omega
  real(kind=sgl)     :: xpc
  real(kind=sgl)     :: ypc
  real(kind=sgl)     :: isangle
  real(kind=sgl)     :: energymin
  real(kind=sgl)     :: energymax
  real(kind=sgl)     :: gammavalue
  real(kind=sgl)     :: axisangle(4)
  real(kind=sgl)     :: beamcurrent
  real(kind=sgl)     :: dwelltime
  real(kind=sgl)     :: hipassw
  real(kind=sgl)     :: stepX
  real(kind=sgl)     :: stepY
  character(1)       :: maskpattern
  character(3)       :: scalingmode
  character(3)       :: Notify
  character(1)       :: keeptmpfile
  character(fnlen)   :: anglefile
  character(fnlen)   :: exptfile
  character(fnlen)   :: masterfile
  character(fnlen)   :: energyfile
  character(fnlen)   :: datafile
  character(fnlen)   :: tmpfile
  character(fnlen)   :: ctffile
  character(fnlen)   :: avctffile
  character(fnlen)   :: angfile
  character(fnlen)   :: eulerfile
  character(fnlen)   :: dictfile
  character(fnlen)   :: maskfile
  character(fnlen)   :: indexingmode
  character(fnlen)   :: refinementNMLfile
  character(fnlen)   :: inputtype
  character(fnlen)   :: HDFstrings(10)
end type DictionaryIndexingNameListType

type, public, extends(DictionaryIndexingNameListType) :: EBSDDINameListType
end type EBSDDINameListType

type, public, extends(DictionaryIndexingNameListType) :: ECPDINameListType
end type ECPDINameListType

type, public, extends(DictionaryIndexingNameListType) :: TKDDINameListType
end type TKDDINameListType


type, public :: MPdataType

end type MPdataType

! class definition
type, public :: DPfile_T
  private 
    ! type(DPdataType),public       :: DPDT
    character(fnlen)              :: DPfile
    character(fnlen)              :: modality = 'unknown'

contains
  private 
    procedure, pass(self) :: setFileName_
    procedure, pass(self) :: writeHDFNameList_
    procedure, pass(self) :: set_Modality_
    procedure, pass(self) :: get_Modality_
    ! procedure, pass(self) :: readDPfile_
    final :: DPfile_destructor

    generic, public :: writeHDFNameList => writeHDFNameList_
    ! generic, public :: readDPfile => readDPfile_
    generic, public :: setFileName => setFileName_
    generic, public :: setModality => set_Modality_
    generic, public :: getModality => get_Modality_

end type DPfile_T

! the constructor routine for this class 
interface DPfile_T
  module procedure DPfile_constructor
end interface DPfile_T

contains

!--------------------------------------------------------------------------
type(DPfile_T) function DPfile_constructor( ) result(DPfile)
!! author: MDG 
!! version: 1.0 
!! date: 03/31/20
!!
!! constructor for the DPfile_T Class; reads the name list 
 
IMPLICIT NONE

end function DPfile_constructor

!--------------------------------------------------------------------------
subroutine DPfile_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 03/31/20
!!
!! destructor for the DPfile_T Class
 
IMPLICIT NONE

type(DPfile_T), INTENT(INOUT)  :: self 

end subroutine DPfile_destructor

!--------------------------------------------------------------------------
function get_Modality_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/31/20
!!
!! get Modality from the DPfile_T class

IMPLICIT NONE 

class(DPfile_T), INTENT(INOUT)     :: self
character(fnlen)                   :: out

out = self%Modality

end function get_Modality_

!--------------------------------------------------------------------------
subroutine set_Modality_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/31/20
!!
!! set Modality in the DPfile_T class

IMPLICIT NONE 

class(DPfile_T), INTENT(INOUT)     :: self
character(*), INTENT(IN)           :: inp

self%Modality = inp

end subroutine set_Modality_

!--------------------------------------------------------------------------
subroutine setFileName_(self, DPfile)
!! author: MDG 
!! version: 1.0 
!! date: 03/31/20
!!
!! set the Master Pattern file name 

IMPLICIT NONE 

class(DPfile_T), INTENT(INOUT)   :: self
character(fnlen), INTENT(IN)     :: DPfile

self%DPfile = trim(DPfile)

end subroutine setFileName_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames, emnl)
!! author: MDG 
!! version: 1.0 
!! date: 03/31/20
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 
use mod_IO
use ISO_C_BINDING

IMPLICIT NONE

class(DPfile_T), INTENT(INOUT)                      :: self 
type(HDF_T), INTENT(INOUT)                          :: HDF
type(HDFnames_T), INTENT(INOUT)                     :: HDFnames
class(DictionaryIndexingNameListType), INTENT(INOUT):: emnl 

type(IO_T)                                          :: Message
integer(kind=irg)                                   :: n_int, n_real
integer(kind=irg)                                   :: hdferr
integer(kind=irg),allocatable                       :: io_int(:)
real(kind=sgl),allocatable                          :: io_real(:)
character(20),allocatable                           :: intlist(:), reallist(:)
character(fnlen)                                    :: dataset, sval(1),groupname
character(fnlen,kind=c_char)                        :: line2(1), line10(10)
logical                                             :: g_exists, overwrite=.TRUE., isEBSD=.FALSE., &
                                                       isECP=.FALSE., isTKD=.FALSE., isEBSDSHT=.FALSE.

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

select type(emnl) 
  type is (EBSDDINameListType)
    isEBSD = .TRUE.
    n_int = 19
    n_real = 15
    allocate( io_int(n_int), intlist(n_int), io_real(n_real), reallist(n_real) )
  type is (ECPDINameListType)
    isECP = .TRUE.
    n_int = 19
    n_real = 15
    allocate( io_int(n_int), intlist(n_int), io_real(n_real), reallist(n_real) )
  class default 
    call Message%printError('writeHDFNameList', 'unknown name list type requested')
end select

! write all the single integers
io_int = (/ emnl%ncubochoric, emnl%numexptsingle, emnl%numdictsingle, emnl%ipf_ht, &
            emnl%ipf_wd, emnl%nnk, emnl%maskradius, emnl%numsx, emnl%numsy, emnl%binning, &
            emnl%nthreads, emnl%devid, emnl%platid, emnl%nregions, emnl%nnav, &
            emnl%nosm, emnl%nlines, emnl%usenumd, emnl%nism /)
intlist(1) = 'Ncubochoric'
intlist(2) = 'numexptsingle'
intlist(3) = 'numdictsingle'
intlist(4) = 'ipf_ht'
intlist(5) = 'ipf_wd '
intlist(6) = 'nnk'
intlist(7) = 'maskradius'
intlist(8) = 'numsx'
intlist(9) = 'numsy'
intlist(10) = 'binning'
intlist(11) = 'nthreads'
intlist(12) = 'devid'
intlist(13) = 'platid'
intlist(14) = 'nregions'
intlist(15) = 'nnav'
intlist(16) = 'nosm'
intlist(17) = 'nlines'
intlist(18) = 'usenumd'
intlist(19) = 'nism'
call HDF%writeNMLintegers(io_int, intlist, n_int)

io_real = (/ emnl%L, emnl%thetac, emnl%delta, emnl%omega, emnl%xpc, &
             emnl%ypc, emnl%energymin, emnl%energymax, emnl%gammavalue, emnl%StepX, &
             emnl%stepY, emnl%isangle, emnl%beamcurrent, emnl%dwelltime, emnl%hipassw /)
reallist(1) = 'L'
reallist(2) = 'thetac'
reallist(3) = 'delta'
reallist(4) = 'omega'
reallist(5) = 'xpc'
reallist(6) = 'ypc'
reallist(7) = 'energymin'
reallist(8) = 'energymax'
reallist(9) = 'gammavalue'
reallist(10) = 'stepX'
reallist(11) = 'stepY'
reallist(12) = 'isangle'
reallist(13) = 'beamcurrent'
reallist(14) = 'dwelltime'
reallist(15) = 'hipassw'
call HDF%writeNMLreals(io_real, reallist, n_real)

! a 4-vector
dataset = SC_axisangle
hdferr = HDF%writeDatasetFloatArray(dataset, emnl%axisangle, 4)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create axisangle dataset', hdferr)

! an integer 4-vector
dataset = SC_ROI
hdferr = HDF%writeDatasetIntegerArray(dataset, emnl%ROI, 4)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create ROI dataset', hdferr)

! an integer 8-vector
dataset = 'multidevid'
hdferr = HDF%writeDatasetIntegerArray(dataset, emnl%multidevid, 8)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create multidevid dataset', hdferr)

! strings
dataset = SC_maskpattern
line2(1) = emnl%maskpattern
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create maskpattern dataset', hdferr)

dataset = SC_scalingmode
line2(1) = emnl%scalingmode
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create scalingmode dataset', hdferr)

dataset = SC_exptfile
line2(1) = emnl%exptfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create exptfile dataset', hdferr)

dataset = SC_masterfile
line2(1) = emnl%masterfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create masterfile dataset', hdferr)

dataset = SC_energyfile
line2(1) = emnl%energyfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create energyfile dataset', hdferr)

dataset = SC_datafile
line2(1) = emnl%datafile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create datafile dataset', hdferr)

dataset = SC_tmpfile
line2(1) = emnl%tmpfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create tmpfile dataset', hdferr)

dataset = SC_ctffile
line2(1) = emnl%ctffile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create ctffile dataset', hdferr)

dataset = SC_angfile
line2(1) = emnl%angfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create angfile dataset', hdferr)

dataset = SC_anglefile
line2(1) = emnl%anglefile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create anglefile dataset', hdferr)

dataset = SC_eulerfile
line2(1) = emnl%eulerfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create eulerfile dataset', hdferr)

dataset = SC_maskfile
line2(1) = emnl%maskfile
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create maskfile dataset', hdferr)

dataset = SC_inputtype
line2(1) = emnl%inputtype
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create inputtype dataset', hdferr)

dataset = SC_HDFstrings
line10 = emnl%HDFstrings
hdferr = HDF%writeDatasetStringArray(dataset, line10, 10)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create HDFstrings dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end subroutine writeHDFNameList_




!--------------------------------------------------------------------------
subroutine DPfiles_(self, EMsoft, progname)
!! author: MDG 
!! version: 1.0 
!! date: 03/31/20
!!
!! perform the computations

use mod_EMsoft

IMPLICIT NONE 

class(DPfile_T), INTENT(INOUT)       :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 


end subroutine DPfiles_



end module mod_DPfiles