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

module mod_hh
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/25/22
  !!
  !! class definition for the EMhh4 program

use mod_kinds
use mod_global

IMPLICIT NONE 

private :: NEWTON, ANCALC, PANCALC, RKM, DERIV

!=======================================
!=======================================
!=======================================
! below are type definitions for the f90 EMsoft version of the original
! Head&Humble hh4.f90 program, now called EMhh4.f90. The comment line
! before each type def is the original COMMON block definition
!
! COMMON/MAPN/NEW,ZR,ZI,QR(9),QI(9),KRASH 
type MAPN_block
  integer(kind=irg)        :: KRASH, NEW
  real(kind=sgl)           :: ZR, ZI, QR(9), QI(9)
end type MAPN_block

! COMMON/MA/PR(4),PI(4),AR(4,4),AI(4,4),EMR(4,4),EMI(4,4),H(4,4) 
type MA_block
  real(kind=sgl)           :: PR(4), PI(4), AR(4,4), AI(4,4), EMR(4,4), EMI(4,4), H(4,4) 
end type MA_block

! COMMON/MKAP/D1(6,6),EP(3,6),EA(3,3) 
type MKAP_block
 real(kind=sgl)            :: D1(6,6), EP(3,6), EA(3,3) 
end type MKAP_block

! COMMON/MRD/CN(61),X,X1,Y(8),ERROR,Q,KOUNT,D(8),YT(8),DT(8,4),ANO,SKIP 
type MRD_block
 real(kind=sgl)            :: CN(61), X, X1, Y(8), ERROR, Q, D(8), YT(8), DT(8,4), ANO, SKIP 
 integer(kind=irg)         :: KOUNT
end type MRD_block

! COMMON/MT/LU(3),LG(3),LBM(3),LFN(3),LB(3),LB2(3),LB3(3),LB4(3), 
!           LFP(3),LFP1(3),LFP3(3),LS1(3),LS2(3),LS3(3),TLU(3),TLG(3),
!           TLBM(3),TLFN(3),TLB(3),TLB2(3),TLB3(3),TLB4(3),TLFP(3),
!           TLFP1(3),TLFP3(3),TLS1(3),TLS2(3),TLS3(3),LF1(3),
!           LF2(3),LF3(3),LF4(3),TLF1(3),TLF2(3),TLF3(3),TLF4(3) 
type MT_block
 real(kind=sgl)            :: TLU(3), TLG(3), TLBM(3), TLFN(3), TLB(3), TLB2(3), TLB3(3), TLB4(3), TLFP(3), &
                              TLFP1(3), TLFP3(3), TLS1(3), TLS2(3), TLS3(3), TLF1(3), TLF2(3), TLF3(3), TLF4(3) 
 integer(kind=irg)         :: LU(3), LG(3), LBM(3), LFN(3), LB(3), LB2(3), LB3(3), LB4(3), LD, LD2, LD3, LD4, &
                              LFP(3), LFP1(3), LFP3(3), LS1(3), LS2(3), LS3(3), LF1(3), LF2(3), LF3(3), LF4(3), &
                              LQ1, LQ2, LQ3
end type MT_block

! COMMON/MKT/AT(3,3),ATR(3,3)
type MKT_block
 real(kind=sgl)            :: AT(3,3), ATR(3,3)
end type MKT_block

! COMMON/SCALE30/LTEST
type SCALE30_block
 integer(kind=irg)         :: LTEST
end type SCALE30_block

! COMMON/MP/PC(4),AS(4,4),EL(4,4) 
type MP_block
 complex(kind=sgl)         :: PC(4), AS(4,4), EL(4,4) 
end type MP_block

! COMMON/MAP/DC(3,3)
type MAP_block
 real(kind=sgl)            :: DC(3,3)
end type MAP_block

! namelist for the EMhh program
type, public :: hhNameListType
  integer(kind=irg)   :: IROW
  integer(kind=irg)   :: ICOL
  integer(kind=irg)   :: wnum
  integer(kind=sgl)   :: LTEST
  integer(kind=irg)   :: LB(3), LD 
  integer(kind=irg)   :: LB2(3), LD2
  integer(kind=irg)   :: LB3(3), LD3
  integer(kind=irg)   :: LB4(3), LD4
  integer(kind=irg)   :: LU(3)
  integer(kind=irg)   :: LG(3)
  integer(kind=irg)   :: LBM(3)
  integer(kind=irg)   :: LFN(3)
  integer(kind=irg)   :: LFP1(3), LFP(3), LFP3(3)
  integer(kind=irg)   :: LS1(3), LQ1 
  integer(kind=irg)   :: LS2(3), LQ2 
  integer(kind=irg)   :: LS3(3), LQ3 
  real(kind=sgl)      :: kV
  real(kind=sgl)      :: THICK, START, FINISH
  real(kind=sgl)      :: wmin, wmax
  real(kind=sgl)      :: SEP, SEP2
  real(kind=sgl)      :: FAP1, FAP3
  real(kind=sgl)      :: D1row1(6)
  real(kind=sgl)      :: D1row2(6)
  real(kind=sgl)      :: D1row3(6)
  real(kind=sgl)      :: D1row4(6)
  real(kind=sgl)      :: D1row5(6)
  real(kind=sgl)      :: D1row6(6)
  character(fnlen)    :: xtalname
  character(fnlen)    :: outname
  character(fnlen)    :: imageprefix
  character(fnlen)    :: imagetype 
end type hhNameListType

! class definition
type, public :: hh_T
private 
  character(fnlen)        :: nmldeffile = 'EMhh.nml'
  type(hhNameListType)    :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: writeHH4_HDFfile_
  procedure, pass(self) :: hh_writeInfo_
  procedure, pass(self) :: dohh_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: writeHH4_HDFfile => writeHH4_HDFfile_
  generic, public :: hh_writeInfo => hh_writeInfo_
  generic, public :: dohh => dohh_

end type hh_T

! the constructor routine for this class 
interface hh_T
  module procedure hh_constructor
end interface hh_T

contains

!--------------------------------------------------------------------------
type(hh_T) function hh_constructor( nmlfile ) result(hh)
!DEC$ ATTRIBUTES DLLEXPORT :: hh_constructor
!! author: MDG 
!! version: 1.0 
!! date: 01/25/22
!!
!! constructor for the hh_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

if (present(nmlfile)) then 
  call hh%readNameList(nmlfile)
end if 

end function hh_constructor

!--------------------------------------------------------------------------
subroutine hh_destructor(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: hh_destructor
!! author: MDG 
!! version: 1.0 
!! date: 01/25/22
!!
!! destructor for the hh_T Class
 
IMPLICIT NONE

type(hh_T), INTENT(INOUT)  :: self 

call reportDestructor('hh_T')

end subroutine hh_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 01/25/22
!!
!! read the namelist from an nml file for the hh_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(hh_T), INTENT(INOUT)           :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg)       :: IROW
integer(kind=irg)       :: ICOL
integer(kind=irg)       :: LB(3), LD 
integer(kind=irg)       :: LB2(3), LD2
integer(kind=irg)       :: LB3(3), LD3
integer(kind=irg)       :: LB4(3), LD4
integer(kind=irg)       :: LU(3)
integer(kind=irg)       :: LG(3)
integer(kind=irg)       :: LBM(3)
integer(kind=irg)       :: LFN(3)
integer(kind=irg)       :: wnum
integer(kind=irg)       :: LFP1(3), LFP(3), LFP3(3)
integer(kind=irg)       :: LS1(3), LQ1 
integer(kind=irg)       :: LS2(3), LQ2 
integer(kind=irg)       :: LS3(3), LQ3 
integer(kind=sgl)       :: LTEST
real(kind=sgl)          :: kV
real(kind=sgl)          :: THICK, START, FINISH
real(kind=sgl)          :: wmin, wmax
real(kind=sgl)          :: SEP, SEP2
real(kind=sgl)          :: FAP1, FAP3
real(kind=sgl)          :: D1row1(6)
real(kind=sgl)          :: D1row2(6)
real(kind=sgl)          :: D1row3(6)
real(kind=sgl)          :: D1row4(6)
real(kind=sgl)          :: D1row5(6)
real(kind=sgl)          :: D1row6(6)
character(fnlen)        :: xtalname
character(fnlen)        :: outname
character(fnlen)        :: imageprefix
character(fnlen)        :: imagetype 

namelist /hhlist/ IROW, ICOL, LB, LD , LB2, LD2, LB3, LD3, LB4, LD4, LU, LG, LBM, LFN, &
                  wnum, LFP1, LFP, LFP3, LS1, LQ1 , LS2, LQ2 , LS3, LQ3 , LTEST, kV, THICK, START, FINISH, &
                  wmin, wmax, SEP, SEP2, FAP1, FAP3, D1row1, D1row2, D1row3, D1row4, D1row5, D1row6,&
                  xtalname, outname, imageprefix, imagetype

 xtalname = 'undefined'
 outname = 'undefined'
 imageprefix = 'undefined'
 imagetype = 'tiff'
 IROW = 160
 ICOL = 256
 kV = 200.0
 LB = (/1, 0, 1/)
 LD = 2
 LB2 = (/0, 0, 0/) 
 LD2 = 1
 LB3 = (/0, 0, 0/)
 LD3 = 1
 LB4 = (/0, 0, 0/)
 LD4 = 1
 LU = (/1, 1, 1/)
 LG = (/2, 0, 0/)
 LBM = (/0, 0, 1/)
 LFN = (/0, 0, 1/)
 THICK = 5.0
 START = 0.0
 FINISH = 6.0
 wmin = -1.0
 wmax =  1.0
 wnum =  5
 LFP1 = (/0, 0, 0/)
 LFP = (/0, 0, 0/)
 LFP3 = (/0, 0, 0/) 
 LS1 = (/0, 0, 0/) 
 LQ1 = 2
 LS2 = (/0, 0, 0/) 
 LQ2 = 2
 LS3 = (/0, 0, 0/) 
 LQ3 = 2      
 SEP = 2.0
 SEP2 = 2.0
 FAP1 = 0.0
 FAP3 = 0.0
 D1row1 = (/100.0, 80.0, 80.0,  0.0,  0.0,  0.0/)
 D1row2 = (/ 80.0,100.0, 80.0,  0.0,  0.0,  0.0/)
 D1row3 = (/ 80.0, 80.0,100.0,  0.0,  0.0,  0.0/)
 D1row4 = (/  0.0,  0.0,  0.0, 50.0,  0.0,  0.0/)
 D1row5 = (/  0.0,  0.0,  0.0,  0.0, 50.0,  0.0/)
 D1row6 = (/  0.0,  0.0,  0.0,  0.0,  0.0, 50.0/)
 LTEST = 0

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=hhlist)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(outname).eq.'undefined') then
  call Message%printError('EMhh4:',' output HDF file name is undefined in '//nmlfile)
 end if

 if (trim(xtalname).eq.'undefined') then
  call Message%printError('EMhh4:',' xtalname name is undefined in '//nmlfile)
 end if

end if

self%nml%IROW = IROW
self%nml%ICOL = ICOL
self%nml%LB = LB
self%nml%LD  = LD
self%nml%LB2 = LB2
self%nml%LD2 = LD2
self%nml%LB3 = LB3
self%nml%LD3 = LD3
self%nml%LB4 = LB4
self%nml%LD4 = LD4
self%nml%LU = LU
self%nml%LG = LG
self%nml%LBM = LBM
self%nml%LFN = LFN
self%nml%wnum = wnum
self%nml%LFP1 = LFP1
self%nml%LFP = LFP
self%nml%LFP3 = LFP3
self%nml%LS1 = LS1
self%nml%LQ1  = LQ1
self%nml%LS2 = LS2
self%nml%LQ2  = LQ2
self%nml%LS3 = LS3
self%nml%LQ3  = LQ3
self%nml%LTEST = LTEST
self%nml%kV = kV
self%nml%THICK = THICK
self%nml%START = START
self%nml%FINISH = FINISH
self%nml%wmin = wmin
self%nml%wmax = wmax
self%nml%SEP = SEP
self%nml%SEP2 = SEP2
self%nml%FAP1 = FAP1
self%nml%FAP3 = FAP3
self%nml%D1row1 = D1row1 
self%nml%D1row2 = D1row2 
self%nml%D1row3 = D1row3 
self%nml%D1row4 = D1row4 
self%nml%D1row5 = D1row5 
self%nml%D1row6 = D1row6 
self%nml%xtalname = xtalname
self%nml%outname = outname
self%nml%imageprefix = imageprefix
self%nml%imagetype = imagetype 

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 01/25/22
!!
!! pass the namelist for the hh_T Class to the calling program

IMPLICIT NONE 

class(hh_T), INTENT(INOUT)          :: self
type(hhNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 01/25/22
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(hh_T), INTENT(INOUT)        :: self 
type(HDF_T), INTENT(INOUT)        :: HDF
type(HDFnames_T), INTENT(INOUT)   :: HDFnames

integer(kind=irg),parameter       :: n_int = 11, n_real = 10 
integer(kind=irg)                 :: hdferr,  io_int(n_int)
real(kind=sgl)                    :: io_real(n_real)
character(20)                     :: intlist(n_int), reallist(n_real)
character(fnlen)                  :: dataset, sval(1),groupname
character(fnlen,kind=c_char)      :: line2(1)

associate( hhnl => self%nml )

! create the group for this namelist
groupname = trim(HDFnames%get_NMLlist())
hdferr = HDF%createGroup(groupname)

! write all the single integers
io_int      = (/ hhnl%IROW, hhnl%ICOL, hhnl%wnum, hhnl%LTEST, hhnl%LD, hhnl%LD2, &
                 hhnl%LD3, hhnl%LD4, hhnl%LQ1, hhnl%LQ2, hhnl%LQ3 /)
intlist(1)  = 'IROW'
intlist(2)  = 'ICOL'
intlist(3)  = 'wnum'
intlist(4)  = 'LTEST'
intlist(5)  = 'LD'
intlist(6)  = 'LD2'
intlist(7)  = 'LD3'
intlist(8)  = 'LD4'
intlist(9)  = 'LQ1'
intlist(10)  = 'LQ2'
intlist(11)  = 'LQ3'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the integer arrays 
dataset = 'LB'
hdferr = HDF%writeDatasetIntegerArray(dataset, hhnl%LB, 3)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create LB dataset',hdferr)

dataset = 'LB2'
hdferr = HDF%writeDatasetIntegerArray(dataset, hhnl%LB2, 3)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create LB2 dataset',hdferr)

dataset = 'LB3'
hdferr = HDF%writeDatasetIntegerArray(dataset, hhnl%LB3, 3)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create LB3 dataset',hdferr)

dataset = 'LB4'
hdferr = HDF%writeDatasetIntegerArray(dataset, hhnl%LB4, 3)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create LB4 dataset',hdferr)

dataset = 'LU'
hdferr = HDF%writeDatasetIntegerArray(dataset, hhnl%LU, 3)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create LU dataset',hdferr)

dataset = 'LG'
hdferr = HDF%writeDatasetIntegerArray(dataset, hhnl%LG, 3)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create LG dataset',hdferr)

dataset = 'LBM'
hdferr = HDF%writeDatasetIntegerArray(dataset, hhnl%LBM, 3)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create LBM dataset',hdferr)

dataset = 'LFN'
hdferr = HDF%writeDatasetIntegerArray(dataset, hhnl%LFN, 3)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create LFN dataset',hdferr)

dataset = 'LFP1'
hdferr = HDF%writeDatasetIntegerArray(dataset, hhnl%LFP1, 3)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create LFP1 dataset',hdferr)

dataset = 'LFP'
hdferr = HDF%writeDatasetIntegerArray(dataset, hhnl%LFP, 3)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create LFP dataset',hdferr)

dataset = 'LFP3'
hdferr = HDF%writeDatasetIntegerArray(dataset, hhnl%LFP3, 3)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create LFP3 dataset',hdferr)

dataset = 'LS1'
hdferr = HDF%writeDatasetIntegerArray(dataset, hhnl%LS1, 3)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create LS1 dataset',hdferr)

dataset = 'LS2'
hdferr = HDF%writeDatasetIntegerArray(dataset, hhnl%LS2, 3)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create LS2 dataset',hdferr)

dataset = 'LS3'
hdferr = HDF%writeDatasetIntegerArray(dataset, hhnl%LS3, 3)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create LS3 dataset',hdferr)


! write all the single reals
io_real     = (/ hhnl%kV, hhnl%THICK, hhnl%START, hhnl%FINISH, hhnl%wmin, hhnl%wmax, hhnl%SEP, hhnl%SEP2, hhnl%FAP1, hhnl%FAP3 /)
reallist(1) = 'kV'
reallist(2) = 'THICK'
reallist(3) = 'START'
reallist(4) = 'FINISH'
reallist(5) = 'wmin'
reallist(6) = 'wmax'
reallist(7) = 'SEP'
reallist(8) = 'SEP2'
reallist(9) = 'FAP1'
reallist(10) = 'FAP3'
call HDF%writeNMLreals(io_real, reallist, n_real)

! real arrays
dataset = 'D1row1'
hdferr = HDF%writeDatasetFloatArray(dataset, hhnl%D1row1, 6)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create D1row1 dataset',hdferr)

dataset = 'D1row2'
hdferr = HDF%writeDatasetFloatArray(dataset, hhnl%D1row2, 6)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create D1row2 dataset',hdferr)

dataset = 'D1row3'
hdferr = HDF%writeDatasetFloatArray(dataset, hhnl%D1row3, 6)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create D1row3 dataset',hdferr)

dataset = 'D1row4'
hdferr = HDF%writeDatasetFloatArray(dataset, hhnl%D1row4, 6)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create D1row4 dataset',hdferr)

dataset = 'D1row5'
hdferr = HDF%writeDatasetFloatArray(dataset, hhnl%D1row5, 6)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create D1row5 dataset',hdferr)

dataset = 'D1row6'
hdferr = HDF%writeDatasetFloatArray(dataset, hhnl%D1row6, 6)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create D1row6 dataset',hdferr)

! write all the strings
dataset = 'outname'
line2(1) = hhnl%outname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create outname dataset',hdferr)

dataset = SC_xtalname
line2(1) = hhnl%xtalname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create xtalname dataset',hdferr)

dataset = 'imageprefix'
line2(1) = hhnl%imageprefix
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create imageprefix dataset',hdferr)

dataset = 'imagetype'
line2(1) = hhnl%imagetype
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('HDFwriteHH4NameList: unable to create imagetype dataset',hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHH4_HDFfile_(self, EMsoft, hhnl, BF, DF, dstr, tstrb, tstre, legendfiles, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHH4_HDFfile_
!! author: MDG 
!! version: 1.0 
!! date: 08/22/19 adapted from similar code in dictionary indexing program
!!
!! generate HDF5 output file for the EMhh4 program

use mod_EMsoft
use HDF5
use mod_HDFsupport 
use stringconstants

IMPLICIT NONE 

class(hh_T),INTENT(INOUT)                 :: self
type(EMsoft_T),INTENT(INOUT)              :: EMsoft
type(hhNameListType),INTENT(INOUT)        :: hhnl
!f2py intent(in,out) ::  hhnl
real(kind=sgl),INTENT(IN)                 :: BF(hhnl%ICOL, hhnl%IROW, hhnl%wnum)
real(kind=sgl),INTENT(IN)                 :: DF(hhnl%ICOL, hhnl%IROW, hhnl%wnum)
character(11),INTENT(INOUT)               :: dstr
!f2py intent(in,out) ::  dstr
character(15),INTENT(IN)                  :: tstrb
character(15),INTENT(IN)                  :: tstre
character(fnlen),INTENT(IN)               :: legendfiles(hhnl%wnum)
character(fnlen),INTENT(IN)               :: progname

type(HDF_T)                               :: HDF

integer(kind=irg)                         :: irow, icol, imnum, hdferr, i
character(fnlen)                          :: groupname, dataset, hhfile, nmlname
character(3)                              :: lnum

irow = hhnl%IROW
icol = hhnl%ICOL
imnum = hhnl%wnum

HDF = HDF_T()

! Create a new file using the default properties.
hhfile = trim(EMsoft%generateFilePath('EMdatapathname',hhnl%outname))
hdferr =  HDF%createFile(hhfile)
if (hdferr.ne.0) call HDF%error_check('HDF_createFile ', hdferr)

call self%hh_writeInfo(EMsoft, HDF, dstr, tstrb, tstre, progname, hhnl)

!===============================
! create a namelist group to write all the legend files into
groupname = 'LegendFiles'
  hdferr = HDF%createGroup(groupname)

! and write the nml file for this program to the HDF5 file
! read the text file and write the array to the file
  do i=1,hhnl%wnum
    write (lnum,"(I3.3)") i
    dataset = 'Legend_'//lnum
    hdferr = HDF%writeDatasetTextFile(dataset, legendfiles(i))
  end do

! leave this group
  call HDF%pop()
  
!===============================
! and finally write all the actual data sets 
groupname = SC_EMdata
  hdferr = HDF%createGroup(groupname)

dataset = 'BF'
hdferr = HDF%writeDatasetFloatArray(dataset, BF, hhnl%ICOL, hhnl%IROW, hhnl%wnum)
if (hdferr.ne.0) call HDF%error_check('writeHH4_HDFfile: unable to create BF dataset',hdferr)

dataset = 'DF'
hdferr = HDF%writeDatasetFloatArray(dataset, DF, hhnl%ICOL, hhnl%IROW, hhnl%wnum)
if (hdferr.ne.0) call HDF%error_check('writeHH4_HDFfile: unable to create DF dataset',hdferr)

! close the file 
call HDF%popall()
  
end subroutine writeHH4_HDFfile_

!--------------------------------------------------------------------------
subroutine hh_writeInfo_(self, EMsoft, HDF, dstr, tstrb, tstre, progname, hhnl)
!DEC$ ATTRIBUTES DLLEXPORT :: hh_writeInfo_
!! author: MDG 
!! version: 1.0 
!! date: 08/22/19 adapted from similar code in dictionary indexing program
!!
!! write general information fields to the h5ebsd file, including EMsoft specific fields

use mod_EMsoft
use HDF5
use mod_HDFsupport
use mod_HDFnames
use stringconstants

IMPLICIT NONE

class(hh_T),INTENT(INOUT)                           :: self
type(EMsoft_T),INTENT(INOUT)                        :: EMsoft
type(HDF_T),INTENT(INOUT)                           :: HDF
character(11),INTENT(INOUT)                         :: dstr
!f2py intent(in,out) ::  dstr
character(15),INTENT(IN)                            :: tstrb
character(15),INTENT(IN)                            :: tstre
character(fnlen),INTENT(IN)                         :: progname
type(hhNameListType),INTENT(INOUT)                  :: hhnl
!f2py intent(in,out) ::  hhnl

type(HDFnames_T)                                    :: HDFnames

character(fnlen, KIND=c_char),allocatable,TARGET    :: stringarray(:)
character(fnlen)                                    :: dataset, groupname 
integer(kind=irg)                                   :: hdferr

allocate(stringarray(1))

! add the EMsoft header 
  call HDF%writeEMheader(EMsoft, dstr, tstrb, tstre, progname)

! create a namelist group to write all the namelist files into
groupname = SC_NMLfiles
  hdferr = HDF%createGroup(groupname)

! and write the nml file for this program to the HDF5 file
! read the text file and write the array to the file
dataset = 'EMHH4NML'
  hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%nmldeffile)

! leave this group
  call HDF%pop()
  
! create a namelist group to write all the namelist files into
HDFnames = HDFnames_T()
  call HDFnames%set_NMLlist(SC_NMLparameters)
  call self%writeHDFNameList_(HDF, HDFnames)

end subroutine hh_writeInfo_

!--------------------------------------------------------------------------
recursive subroutine NEWTON(MAPN) 
!DEC$ ATTRIBUTES DLLEXPORT :: NEWTON
!! author: Head, A.K. and Humble, P. and Clarebrough, L.M. and Morton, A.J. and Forwood, C.T. 
!! version: 1.0 
!! date: 08/13/19 adapted from .f77 original
!!
!! Newton routine to find root of 8th order polynomial

IMPLICIT NONE

type(MAPN_block),INTENT(INOUT)    :: MAPN
!f2py intent(in,out) ::  MAPN

integer(kind=irg)                 :: KOUNT, KONVRG, J, M
real(kind=sgl)                    :: XR, XI, YR, YI, TR, TI, F
  
! FIND A ROOT OF THE POLYNOMIAL OF THE EIGHTH ORDER
 KONVRG=0
loop: do KOUNT=1,70 
  XR=0.0
  XI=0.0
  YR=0.0
  YI=0.0
  do J=1,MAPN%NEW
   TR=MAPN%ZR*YR-MAPN%ZI*YI+XR 
   YI=MAPN%ZR*YI+MAPN%ZI*YR+XI 
   YR=TR 
   M=MAPN%NEW+1-J 
   TR=MAPN%ZR*XR-MAPN%ZI*XI+MAPN%QR(M)
   TI=MAPN%ZR*XI+MAPN%ZI*XR+MAPN%QI(M)
   if (KONVRG.ne.0) then 
    MAPN%QR(M)=XR
    MAPN%QI(M)=XI
   end if
   XR=TR 
   XI=TI 
  end do
  if (KONVRG.eq.0) then 
   F=1.0/(YR**2+YI**2) 
   TR=F*(XR*YR+XI*YI)
   TI=F*(XI*YR-XR*YI)
   MAPN%ZR=MAPN%ZR-TR
   MAPN%ZI=MAPN%ZI-TI
   if (((TR**2+TI**2)/(MAPN%ZR**2+MAPN%ZI**2)-1.0E-12).le.0.0) then
    KONVRG=1
   end if
  else
   EXIT loop 
  end if
 end do loop
!
 if (KOUNT.eq.70) then 
   MAPN%KRASH=-70 
 else
  if ((ABS(MAPN%ZR)-1.0E5*ABS(MAPN%ZI)).le.0.0) then 
   MAPN%NEW = MAPN%NEW-1
  else
   MAPN%KRASH=10-MAPN%NEW
  end if
 end if
end subroutine NEWTON

!--------------------------------------------------------------------------
recursive subroutine ANCALC(MAP, MKAP, MAPN, MA, SCALE30) 
!DEC$ ATTRIBUTES DLLEXPORT :: ANCALC
!! author: Head, A.K. and Humble, P. and Clarebrough, L.M. and Morton, A.J. and Forwood, C.T. 
!! version: 1.0 
!! date: 08/13/19 adapted from .f77 original
!!
!! Displacement calculations for anisotropic elasticity

IMPLICIT NONE

type(MAP_block),INTENT(INOUT)     :: MAP
!f2py intent(in,out) ::  MAP
type(MKAP_block),INTENT(INOUT)    :: MKAP
!f2py intent(in,out) ::  MKAP
type(MAPN_block),INTENT(INOUT)    :: MAPN
!f2py intent(in,out) ::  MAPN
type(MA_block),INTENT(INOUT)      :: MA
!f2py intent(in,out) ::  MA
type(SCALE30_block),INTENT(INOUT) :: SCALE30
!f2py intent(in,out) ::  SCALE30

integer(kind=irg)              :: I, J, K, L, M, N, LT, KQ, KR, KS, KT, NJ, I1, I2, J1, J2, K1, K2, LP, LQ, KP, NL, ML
integer(kind=irg),parameter    :: L1(6)=(/1,2,3,2,3,1/), L2(6)=(/1,2,3,3,1,2/), &
                                  L3(3,3)= reshape( (/1,6,5,6,2,4,5,4,3/), (/3,3/) ), & 
                                  N1(4)=(/2,4,2,1/), N2(4)=(/3,1,4,2/), N3(4)=(/4,3,1,3/), &
                                  NN(3)=(/6,2,4/), MM(3)=(/1,6,5/), NP(3) = (/2,3,1/), NQ(3) = (/3,1,2/)

real(kind=sgl)                 :: C(6,6), EE(3,6), EN(3,3), QM(7,4), G(9), E(9), F(9), HI(12), &
                                  DR(4,4), DI(4,4), B(4,4), ELR(4,4), ELI(4,4), X, Y, Z, PRK, PIK, SQR, SQI, XR, XI, YR, YI, &
                                  DELR, DELI, AUMR, AUMI, DEL
  
! CALCULATE DISPLACEMENTS 
 do M=1,6
  I=L1(M)
  J=L2(M)
  do N=1,M
   K=L1(N)
   L=L2(N)
   X=0.0
   DO LP=1,3 
    Y=0.0
    DO LQ=1,3
     LT=L3(LP,LQ)
     Y=Y+MAP%DC(J,LQ)*(MAP%DC(K,1)*(MAP%DC(L,1)*MKAP%D1(LT,1)+MAP%DC(L,2)*MKAP%D1(LT,6)+MAP%DC(L,3)*MKAP%D1(LT,5)) &
        +MAP%DC(K,2)*(MAP%DC(L,1)*MKAP%D1(LT,6)+MAP%DC(L,2)*MKAP%D1(LT,2)+MAP%DC(L,3)*MKAP%D1(LT,4))  &
        +MAP%DC(K,3)*(MAP%DC(L,1)*MKAP%D1(LT,5)+MAP%DC(L,2)*MKAP%D1(LT,4)+MAP%DC(L,3)*MKAP%D1(LT,3)))
    end do
    X=X+MAP%DC(I,LP)*Y
   end do
   C(M,N)=X
   C(N,M)=X
  end do
 end do
 if (SCALE30%LTEST.eq.1) then
  write(dataunit,"('  C-DC ',6F10.6)") ((C(I,J),J=1,6),I=1,6)
 end if 
 G(1)=C(5,5)
 G(2)=2.0*C(4,5)
 G(3)=C(4,4)
 G(4)=C(6,6)
 G(5)=2.0*C(2,6)
 G(6)=C(2,2)
 G(7)=C(1,1)
 G(8)=2.0*C(1,6)
 G(9)=C(6,6)
 E(1)=C(5,6)
 E(2)=C(2,5)+C(4,6)
 E(3)=C(2,4)
 E(4)=C(1,5)
 E(5)=C(5,6)+C(1,4)
 E(6)=C(4,6)
 E(7)=C(1,6)
 E(8)=C(6,6)+C(1,2)
 E(9)=C(2,6)
 MAPN%QR = 0.0
 MAPN%QI = 0.0
 do KQ=1,3
  do KR=1,3
   do KS=1,3
    KT=KQ+KR+KS-2
    MAPN%QR(KT)=MAPN%QR(KT)+G(KQ)*G(KR+3)*G(KS+6)+2.0*E(KQ)*E(KR+3)*E(KS+6) &
                -E(KQ)*E(KR)*G(KS+6)-E(KQ+3)*E(KR+3)*G(KS+3)-E(KQ+6)*E(KR+6)*G(KS)
   end do
  end do
 end do
 if (SCALE30%LTEST.eq.1) then
  write(dataunit,"('  QR- HH4 ',F30.15)") (MAPN%QR(KP),KP=1,7)
 end if 
 MAPN%QR=MAPN%QR/MAPN%QR(7)
 MAPN%KRASH=0
 MAPN%NEW=7
 MAPN%ZR=0.1
 MAPN%ZI=1.0
 CALL NEWTON(MAPN)
 if (MAPN%KRASH.lt.0) then
  STOP ' Solution of sextic does not converge '
 end if
 if (MAPN%KRASH.gt.0) then
  STOP ' Sextic has real root '
 end if
 MA%PR(1)=MAPN%ZR
 MA%PI(1)=ABS(MAPN%ZI)
 MAPN%ZI=-MAPN%ZI
 CALL NEWTON(MAPN)
 if (MAPN%KRASH.lt.0) then
  STOP ' Solution of sextic does not converge '
 end if
 if (MAPN%KRASH.gt.0) then
  STOP ' Sextic has real root '
 end if
 MAPN%ZR=0.5
 MAPN%ZI=0.9
 CALL NEWTON(MAPN)
 if (MAPN%KRASH.lt.0) then
  STOP ' Solution of sextic does not converge '
 end if
 if (MAPN%KRASH.gt.0) then
  STOP ' Sextic has real root '
 end if
 MA%PR(2)=MAPN%ZR
 MA%PI(2)=ABS(MAPN%ZI)
 MAPN%ZI=-MAPN%ZI
 CALL NEWTON(MAPN)
 if (MAPN%KRASH.lt.0) then
  STOP ' Solution of sextic does not converge '
 end if
 if (MAPN%KRASH.gt.0) then
  STOP ' Sextic has real root '
 end if
 MAPN%ZR=-MAPN%ZR
 CALL NEWTON(MAPN)
 if (MAPN%KRASH.lt.0) then
  STOP ' Solution of sextic does not converge '
 end if
 if (MAPN%KRASH.gt.0) then
  STOP ' Sextic has real root '
 end if
 MA%PR(3)=MAPN%ZR
 MA%PI(3)=ABS(MAPN%ZI) 
 MAPN%ZR=-C(4,5)/C(4,4) 
 MAPN%ZI=SQRT(ABS(C(4,4)*C(5,5)-C(4,5)**2))/C(4,4)
 do N=1,2 
   if (((MAPN%ZR-MA%PR(N))**2+(MAPN%ZI-MA%PI(N))**2-(MAPN%ZR-MA%PR(N+1))**2-(MAPN%ZI-MA%PI(N+1))**2).lt.0.0) then
    Z=MA%PR(N) 
    MA%PR(N)=MA%PR(N+1) 
    MA%PR(N+1)=Z 
    Z=MA%PI(N) 
    MA%PI(N)=MA%PI(N+1) 
    MA%PI(N+1)=Z 
   end if
 end do
!
 if (SCALE30%LTEST.eq.1) then
  do I=1,3
   WRITE(dataunit,"(' Roots HH4 ',2F20.15)") MA%PR(I),MA%PI(I)
  end do 
 end if 
 DO K=1,3 
  I=NP(K)
  L=NQ(K)
  PRK=MA%PR(K)
  PIK=MA%PI(K)
  SQR=PRK**2-PIK**2
  SQI=2.0*PRK*PIK
  DR(1,1)=C(1,1)+PRK*2.0*C(1,6)+SQR*C(6,6)
  DR(2,2)=C(6,6)+PRK*2.0*C(2,6)+SQR*C(2,2)
  DR(3,3)=C(5,5)+PRK*2.0*C(4,5)+SQR*C(4,4)
  DR(1,2)=C(1,6)+PRK*(C(1,2)+C(6,6))+SQR*C(2,6)
  DR(2,1)=DR(1,2)
  DR(1,3)=C(1,5)+PRK*(C(1,4)+C(5,6))+SQR*C(4,6)
  DR(3,1)=DR(1,3)
  DR(2,3)=C(5,6)+PRK*(C(4,6)+C(2,5))+SQR*C(2,4)
  DR(3,2)=DR(2,3)
  DI(1,1)=PIK*2.0*C(1,6)+SQI*C(6,6)
  DI(2,2)=PIK*2.0*C(2,6)+SQI*C(2,2)
  DI(3,3)=PIK*2.0*C(4,5)+SQI*C(4,4)
  DI(1,2)=PIK*(C(1,2)+C(6,6))+SQI*C(2,6)
  DI(2,1)=DI(1,2)
  DI(1,3)=PIK*(C(1,4)+C(5,6))+SQI*C(4,6)
  DI(3,1)=DI(1,3)
  DI(2,3)=PIK*(C(4,6)+C(2,5))+SQI*C(2,4)
  DI(3,2)=DI(2,3)
  do J=1,3
   M=NP(J)
   N=NQ(J)
   MA%AR(J,K)=DR(I,M)*DR(L,N)-DI(I,M)*DI(L,N)-DR(I,N)*DR(L,M)+DI(I,N)*DI(L,M)
   MA%AI(J,K)=DR(I,M)*DI(L,N)+DI(I,M)*DR(L,N)-DR(I,N)*DI(L,M)-DI(I,N)*DR(L,M)
  end do
 end do
 if (SCALE30%LTEST.eq.1) then
  do J=1,3
   do K=1,3
    WRITE(dataunit,"('   Vektoren as HH4 ',2F20.12)") MA%AR(J,K),MA%AI(J,K)
   end do
  end do 
 end if 
 do J=1,3
  NJ=NN(J)
  do K=1,3
   XR=0.0
   XI=0.0
   do L=1,3
    NL=NN(L)
    ML=MM(L)
    YR=C(NJ,ML)+C(NJ,NL)*MA%PR(K)
    YI=C(NJ,NL)*MA%PI(K) 
    XR=XR+YR*MA%AR(L,K)-YI*MA%AI(L,K)
    XI=XI+YI*MA%AR(L,K)+YR*MA%AI(L,K)
   end do
   ELR(J,K)=XR
   ELI(J,K)=XI
  end do
 end do
 do J=1,3
  J1=NP(J)
  J2=NQ(J)
  DO K=1,3
   K1=NP(K)
   K2=NQ(K)
   MA%EMR(K,J)=ELR(J1,K1)*ELR(J2,K2)-ELI(J1,K1)*ELI(J2,K2)-ELR(J1,K2)*ELR(J2,K1)+ELI(J1,K2)*ELI(J2,K1)
   MA%EMI(K,J)=ELR(J1,K1)*ELI(J2,K2)+ELI(J1,K1)*ELR(J2,K2)-ELR(J1,K2)*ELI(J2,K1)-ELI(J1,K2)*ELR(J2,K1)
  end do
 end do
 DELR=0.0
 DELI=0.0
 do J=1,3
  DELR=DELR+ELR(3,J)*MA%EMR(J,3)-ELI(3,J)*MA%EMI(J,3)
  DELI=DELI+ELR(3,J)*MA%EMI(J,3)+ELI(3,J)*MA%EMR(J,3)
 end do
 AUMR=DELR/(DELR**2+DELI**2)
 AUMI=-DELI/(DELR**2+DELI**2)
 do J=1,3
  do K=1,3
   X=MA%EMR(J,K)*AUMR-MA%EMI(J,K)*AUMI
   MA%EMI(J,K)=MA%EMR(J,K)*AUMI+MA%EMI(J,K)*AUMR
   MA%EMR(J,K)=X
  end do
 end do
 do I=1,3
  do J=1,3
   B(I,J)=-sum(MA%AR(I,1:3)*MA%EMI(1:3,J))-sum(MA%AI(I,1:3)*MA%EMR(1:3,J))
  end do
 end do
 do I=1,3
  I1=NP(I)
  I2=NQ(I)
  do J=1,3
   J1=NP(J)
   J2=NQ(J)
   MA%H(I,J)=B(I1,J1)*B(I2,J2)-B(I1,J2)*B(I2,J1)
  end do
 end do
 DEL=B(3,1)*MA%H(3,1)+B(3,2)*MA%H(3,2)+B(3,3)*MA%H(3,3)
 MA%H = MA%H/DEL
 MAPN%QR(8:9)=0.0 
 MAPN%QI(8:9)=0.0 
 MA%PR(4)=0.0 
 MA%PI(4)=0.0 
 MA%AR(1:4,4)=0.0 
 MA%AR(4,1:4)=0.0 
 MA%AI(1:4,4)=0.0 
 MA%AI(4,1:4)=0.0 
 MA%EMR(1:4,4)=0.0
 MA%EMR(4,1:4)=0.0
 MA%EMI(1:4,4)=0.0
 MA%EMI(4,1:4)=0.0
 MA%H(1:4,4)=0.0
 MA%H(4,1:4)=0.0
end subroutine ANCALC

!--------------------------------------------------------------------------
recursive subroutine PANCALC(MAP, MKAP, MAPN, MA, MP, SCALE30) 
!DEC$ ATTRIBUTES DLLEXPORT :: PANCALC
!! author: Head, A.K. and Humble, P. and Clarebrough, L.M. and Morton, A.J. and Forwood, C.T. 
!! version: 1.0 
!! date: 08/13/19 adapted from .f77 original
!!
!! Displacement calculations for anisotropic elasticity with piezoelectric effect

IMPLICIT NONE

type(MAP_block),INTENT(INOUT)     :: MAP
!f2py intent(in,out) ::  MAP
type(MKAP_block),INTENT(INOUT)    :: MKAP
!f2py intent(in,out) ::  MKAP
type(MAPN_block),INTENT(INOUT)    :: MAPN
!f2py intent(in,out) ::  MAPN
type(MA_block),INTENT(INOUT)      :: MA
!f2py intent(in,out) ::  MA
type(MP_block),INTENT(INOUT)      :: MP
!f2py intent(in,out) ::  MP
type(SCALE30_block),INTENT(INOUT) :: SCALE30
!f2py intent(in,out) ::  SCALE30


!**************************************************** 
!*     SUBROUTINE PANCALC has been extended for a   * 
!*     piezoelectric crystal with a core-charge "Q" * 
!**************************************************** 
complex(kind=sgl)           :: PC, AS, EL
complex(kind=sgl)           :: A(4,4), AY(4), AXA(4,4), AXL(4,4), LXL(4,4), MXX(4), MXXX(4,4)
real(kind=sgl)              :: C(6,6), EE(3,6), EN(3,3), QM(7,4), G(9), E(9), F(9), HI(12),  &
                               PR(4), PI(4), X, Y, Z, PR1, PI1, PR2, PI2, PR3, PI3, PR4, PI4, AZ, &
                               BY, CY, D, AZZ, BZ, CZ, DZ
integer(kind=irg)           :: I, J, K, L, M, N, LP, LQ, LT, LV, KP, KQ, KR, KS, KT, &
                               MX, MY, MZ, NX, NY, NZ, NJ, NL, ML, K1, K2, K3, J1, I1, I2
integer(kind=irg),parameter :: L1(6)=(/1,2,3,2,3,1/), L2(6)=(/1,2,3,3,1,2/), &
                               L3(3,3)= reshape( (/1,6,5,6,2,4,5,4,3/), (/3,3/)), & 
                               N1(4)=(/2,1,1,1/), N2(4)=(/3,3,2,2/), N3(4)=(/4,4,4,3/), &
                               NN(3)=(/6,2,4/), MM(3)=(/1,6,5/)
logical                     :: goto1000

 if (SCALE30%LTEST.eq.1) then
  write(6,"('  D1--',6F8.4)") ((MKAP%D1(I,J),J=1,6),I=1,6) 
  write(6,"('  EP--',6F8.4)") ((MKAP%EP(K,L),L=1,6),K=1,3) 
  write(6,"('  EA--',3F8.4)") ((MKAP%EA(M,N),N=1,3),M=1,3) 
 end if 
!*******************************************************************
!*      Elasticity tensor in DC-reference frame                    *
!*******************************************************************
 do M=1,6 
  I=L1(M) 
  J=L2(M) 
  do N=1,M 
   K=L1(N) 
   L=L2(N) 
   X=0.0 
   do LP=1,3
    Y=0.0 
    do LQ=1,3
     LT=L3(LP,LQ)
     Y=Y+MAP%DC(J,LQ)*(MAP%DC(K,1)*(MAP%DC(L,1)*MKAP%D1(LT,1)+MAP%DC(L,2)*MKAP%D1(LT,6)+MAP%DC(L,3)*MKAP%D1(LT,5)) &
        +MAP%DC(K,2)*(MAP%DC(L,1)*MKAP%D1(LT,6)+MAP%DC(L,2)*MKAP%D1(LT,2)+MAP%DC(L,3)*MKAP%D1(LT,4)) &
        +MAP%DC(K,3)*(MAP%DC(L,1)*MKAP%D1(LT,5)+MAP%DC(L,2)*MKAP%D1(LT,4)+MAP%DC(L,3)*MKAP%D1(LT,3)))
    end do
    X=X+MAP%DC(I,LP)*Y
   end do
   C(M,N)=X
   C(N,M)=X
  end do
 end do
 where (abs(C).lt.1.0e-8) C=0.0
 if (SCALE30%LTEST.eq.1) then 
  write(6,"(' C-DC',6F8.4)") ((C(M,N),N=1,6),M=1,6)
 end if 
!*******************************************************************
!*       Piezo-tensor --E-- in DC-system                           *
!*******************************************************************
 do L=1,3 
  do M=1,6 
   I=L1(M) 
   J=L2(M) 
   X=0.0 
   do LV=1,3
    X=X+MAP%DC(L,LV)*(MAP%DC(I,1)*(MAP%DC(J,1)*MKAP%EP(LV,1)+MAP%DC(J,2)*MKAP%EP(LV,6)+MAP%DC(J,3)*MKAP%EP(LV,5)) &
       +MAP%DC(I,2)*(MAP%DC(J,1)*MKAP%EP(LV,6)+MAP%DC(J,2)*MKAP%EP(LV,2)+MAP%DC(J,3)*MKAP%EP(LV,4)) &
       +MAP%DC(I,3)*(MAP%DC(J,1)*MKAP%EP(LV,5)+MAP%DC(J,2)*MKAP%EP(LV,4)+MAP%DC(J,3)*MKAP%EP(LV,3)))
   end do
   EE(L,M)=X 
  end do
 end do
 where (ABS(EE).lt.1E-8) EE=0.0
 if (SCALE30%LTEST.eq.1) then
  write(6,"(' E-DC',6F8.4)")  ((EE(M,N),N=1,6),M=1,3)
 end if 
!*******************************************************************
!*        Dielectric tensor  in DC-system                          *
!*******************************************************************
 do I=1,3 
  do J=1,3 
   EN(I,J)=0.0
   do K=1,3 
    EN(I,J)=EN(I,J)+MAP%DC(I,K)*(MAP%DC(J,1)*MKAP%EA(K,1)+MAP%DC(J,2)*MKAP%EA(K,2)+MAP%DC(J,3)*MKAP%EA(K,3)) 
   end do
   EN(J,I)=EN(I,J) 
  end do
 end do
 where (ABS(EN).lt.1E-8) EN=0.0
 if (SCALE30%LTEST.eq.1) then
  write(6,"(' DI-DC',3F8.4)") ((EN(M,N),N=1,3),M=1,3) 
 end if 
!*******************************************************************
!        Computation of eighth-order polynomial.  The 
!        4 by 4 matrix is decomposed into four sixth-order
!        polynomials.  (Original German text below)
!*******************************************************************
!*       BERECHNUNG DES POLYNOMS ACHTEN GRADES. DIE 4 X 4          *
!*       MATRIX WIRD ENTWICKELT, D.H. ES WERDEN ZUERST VIER        *
!*       POLYNOME SECHSTEN GRADES BERECHNET, DANN JEWEILS          *
!*       MIT DEM ENTSPRECHENDEN ELEMENT DER VIERTEN ZEILE          *
!*       MULTIPLIZIERT. ADDITION DER VIER POLYNOME ACHTEN          *
!*       GRADES ERGEBEN DANN DAS GESUCHTE POLYNOM.                 *
!*******************************************************************
  MAPN%QR=0.0
  MAPN%QI=0.0
!*******************************************************************
!    BERUECKSICHTIGUNG DER ENTKOPPELTEN  4 X 4   MATRIX            *
! 
!                         X  X  0  0                               *
!                         X  X  0  0                               *
!                         0  0  X  X                               *
!                         0  0  X  X                                * 
!*******************************************************************
 X=C(1,4)+C(1,5)+C(2,4)+C(4,6)+C(5,6)+EE(1,1)+EE(1,2)+EE(1,6)+EE(2,1)+EE(2,2)+EE(2,6)
 goto1000 = .FALSE.
 if (X.eq.0.0) then
  MAPN%QR(1)=C(1,1)*C(6,6)-C(1,6)**2
  MAPN%QR(2)=2*(C(1,1)*C(2,6)-C(1,6)*C(1,2))
  MAPN%QR(3)=C(1,1)*C(2,2)+2*(C(1,6)*C(2,6)-C(1,2)*C(6,6))-C(1,2)**2
  MAPN%QR(4)=2*(C(1,6)*C(2,2)-C(1,2)*C(2,6))
  MAPN%QR(5)=C(2,2)*C(6,6)-C(2,6)**2
  do KP=1,5
   MAPN%QR(KP)=MAPN%QR(KP)/MAPN%QR(5)
   if (SCALE30%LTEST.eq.1) then 
    write(6,"('  1 QR=',F15.5)") MAPN%QR(KP) 
   end if
  end do
  MAPN%KRASH=0
  MAPN%NEW=5
  MAPN%ZR=0.1 
  MAPN%ZI=1.0 
  CALL NEWTON(MAPN)
  if (MAPN%KRASH.gt.0) then 
   STOP 'POLYNOM(4) HAS A REAL ROOT '
  end if
  if (MAPN%KRASH.lt.0) then 
   STOP 'SOLUTION OF POLYNOM(4) DOES NOT CONVERGE '
  end if
  PR(1)=MAPN%ZR 
  PI(1)=ABS(MAPN%ZI)
  MAPN%ZI=-MAPN%ZI 
  CALL NEWTON(MAPN)
  if (MAPN%KRASH.gt.0) then 
   STOP 'POLYNOM(4) HAS A REAL ROOT '
  end if
  if (MAPN%KRASH.lt.0) then 
   STOP  'SOLUTION OF POLYNOM(4) DOES NOT CONVERGE '
  end if
  MAPN%ZR=0.5 
  MAPN%ZI=0.9 
  CALL NEWTON(MAPN)
  if (MAPN%KRASH.gt.0) then 
   STOP 'POLYNOM(4) HAS A REAL ROOT '
  end if
  if (MAPN%KRASH.lt.0) then 
   STOP 'SOLUTION OF POLYNOM(4) DOES NOT CONVERGE '
  end if
  PR(2)=MAPN%ZR 
  PI(2)=ABS(MAPN%ZI)
  IF (ABS(PR(1)).LT.1E-5)  PR(1)=0.0 
  IF (ABS(1.0-PI(1)).LT.1E-5)  PI(1)=1.0 
  IF (ABS(PR(2)).LT.1E-5)  PR(2)=0.0 
  IF (ABS(1.0-PI(2)).LT.1E-5)  PI(2)=1.0 
  MP%PC(1)=CMPLX(PR(1),PI(1)) 
  MP%PC(2)=CMPLX(PR(2),PI(2)) 
  MP%AS(1,1)=C(6,6)+2*MP%PC(1)*C(2,6)+C(2,2)*MP%PC(1)**2
  MP%AS(2,1)=-C(1,6)-MP%PC(1)*(C(1,2)+C(6,6))-C(2,6)*MP%PC(1)**2
  MP%AS(3,1)=CMPLX(0.0,0.0) 
  MP%AS(4,1)=CMPLX(0.0,0.0) 
  MP%AS(1,2)=-C(1,6)-MP%PC(2)*(C(1,2)+C(6,6))-C(2,6)*MP%PC(2)**2
  MP%AS(2,2)=C(1,1)+2*MP%PC(2)*C(1,6)+C(6,6)*MP%PC(2)**2
  MP%AS(3,2)=CMPLX(0.0,0.0) 
  MP%AS(4,2)=CMPLX(0.0,0.0) 
  do KT=1,5
   MAPN%QR(KT)=0.0 
   MAPN%QI(KT)=0.0 
  end do
  MAPN%QR(1)=C(5,5)*EN(1,1)+EE(1,5)**2
  MAPN%QR(2)=2*(C(5,5)*EN(1,2)+C(4,5)*EN(1,1)+EE(1,5)*EE(1,4)+EE(1,5)*EE(2,5))
  MAPN%QR(3)=C(5,5)*EN(2,2)+4*C(4,5)*EN(1,2)+C(4,4)*EN(1,1)+EE(1,4)**2+2*EE(1,4)*EE(2,5)+EE(2,5)**2+2*EE(1,5)*EE(2,4)
  MAPN%QR(4)=2*(C(4,5)*EN(2,2)+C(4,4)*EN(1,2)+EE(2,4)*EE(1,4)+EE(2,4)*EE(2,5)) 
  MAPN%QR(5)=C(4,4)*EN(2,2)+EE(2,4)**2
  do KP=1,5
   MAPN%QR(KP)=MAPN%QR(KP)/MAPN%QR(5)
   if (SCALE30%LTEST.eq.1) then 
    write(6,"('  2 QR=',F15.5)") MAPN%QR(KP) 
   end if
  end do
  MAPN%KRASH=0
  MAPN%NEW=5
  MAPN%ZR=0.1 
  MAPN%ZI=1.0 
  CALL NEWTON(MAPN)
  if (MAPN%KRASH.gt.0) then 
   STOP 'POLYNOM(4) HAS A REAL ROOT '
  end if
  if (MAPN%KRASH.lt.0) then 
   STOP 'SOLUTION OF POLYNOM(4) DOES NOT CONVERGE '
  end if
  PR(3)=MAPN%ZR 
  PI(3)=ABS(MAPN%ZI)
  MAPN%ZI=-MAPN%ZI 
  CALL NEWTON(MAPN)
  if (MAPN%KRASH.gt.0) then 
   STOP 'POLYNOM(4) HAS A REAL ROOT '
  end if
  if (MAPN%KRASH.lt.0) then 
   STOP 'SOLUTION OF POLYNOM(4) DOES NOT CONVERGE '
  end if
  MAPN%ZR=0.5 
  MAPN%ZI=0.9 
  CALL NEWTON(MAPN)
  if (MAPN%KRASH.gt.0) then 
   STOP 'POLYNOM(4) HAS A REAL ROOT '
  end if
  if (MAPN%KRASH.lt.0) then 
   STOP 'SOLUTION OF POLYNOM(4) DOES NOT CONVERGE '
  end if
  PR(4)=MAPN%ZR 
  PI(4)=ABS(MAPN%ZI)
  if (ABS(PR(3)).LT.1E-5)  PR(3)=0.0 
  if (ABS(PR(4)).LT.1E-5)  PR(4)=0.0 
  if (ABS(1.0-PI(3)).LT.1E-5) PI(3)=1.0
  if (ABS(1.0-PI(4)).LT.1E-5) PI(4)=1.0
  MP%PC(3)=CMPLX(PR(3),PI(3)) 
  MP%PC(4)=CMPLX(PR(4),PI(4)) 
  MP%AS(1,3)=CMPLX(0.0,0.0) 
  MP%AS(2,3)=CMPLX(0.0,0.0) 
  MP%AS(3,3)=-EN(1,1)-2*MP%PC(3)*EN(1,2)-EN(2,2)*MP%PC(3)**2
  MP%AS(4,3)=-EE(1,5)-MP%PC(3)*(EE(1,4)+EE(2,5))-EE(2,4)*MP%PC(3)**2
  if (ABS(CMPLX(0.0,0.0)-MP%AS(3,3)).lt.1E-5) then
   if (ABS(CMPLX(0.0,0.0)-MP%AS(4,3)).lt.1E-5) then 
    MP%AS(3,3)=CMPLX(1.0,1.0)*(1./(2.*SQRT(C(4,4))))
    MP%AS(4,3)=CMPLX(0.0,0.0) 
   end if 
  end if
  MP%AS(1,4)=CMPLX(0.0,0.0) 
  MP%AS(2,4)=CMPLX(0.0,0.0) 
  MP%AS(3,4)=-EE(1,5)-MP%PC(4)*(EE(1,4)+EE(2,5))-EE(2,4)*MP%PC(4)**2
  MP%AS(4,4)=C(5,5)+2*MP%PC(4)*C(4,5)+C(4,4)*MP%PC(4)**2
  if (ABS(CMPLX(0.0,0.0)-MP%AS(4,4)).lt.1E-5) then
   if (ABS(CMPLX(0.0,0.0)-MP%AS(3,4)).lt.1E-5) then 
    MP%AS(3,4)=CMPLX(0.0,0.0) 
    MP%AS(4,4)=CMPLX(0.0,0.0) 
   end if 
  end if
  if (SCALE30%LTEST.eq.1) then 
   do I=1,4
    WRITE(6,"('  QUARTIC  PC= ',2F15.5)")  MP%PC(I)
   end do
   do I=1,4
    do J=1,4
     write(6,"('  QUARTIC  AS ',I1,' = ',2F20.12)") I, MP%AS(J,I)
    end do
   end do
  end if
! the original code has a goto statement here, skipping many lines below.
! we will replace this with a logical test.
!  GOTO 1000
  goto1000 = .TRUE.
 end if 
! 
!***********************************************************************
! 
!     AB  HIER WERDEN DIE "NORMALEN" FAELLE BEHANDELT                  *
!                                                                      *
!           X  X  X  X                   X  X  0  X                    *
!           X  X  X  X                   X  X  0  X                    *
!           X  X  X  X      ODER         0  0  X  0                    *
!           X  X  X  X                   X  X  0  X                    *
! 
!*********************************************************************
if (goto1000.eqv..FALSE.) then 
 QM=0.0
 G(1)=C(5,5) 
 G(2)=2.0*C(4,5) 
 G(3)=C(4,4) 
 G(4)=C(6,6) 
 G(5)=2.0*C(2,6) 
 G(6)=C(2,2) 
 G(7)=C(1,1) 
 G(8)=2.0*C(1,6) 
 G(9)=C(6,6) 
 E(1)=C(5,6) 
 E(2)=C(4,6)+C(2,5)
 E(3)=C(2,4) 
 E(4)=C(1,5) 
 E(5)=C(1,4)+C(5,6)
 E(6)=C(4,6) 
 E(7)=C(1,6) 
 E(8)=C(1,2)+C(6,6)
 E(9)=C(2,6) 
 do KQ=1,3
  do KR=1,3
   do KS=1,3
    KT=KQ+KR+KS-2 
    QM(KT,1)=QM(KT,1)+G(KQ)*G(KR+3)*G(KS+6)+2.0*E(KQ)*E(KR+3)*E(KS+6)- &
             E(KQ)*E(KR)*G(KS+6)-E(KQ+3)*E(KR+3)*G(KS+3)-E(KQ+6)*E(KR+6)*G(KS)
   end do 
  end do 
 end do 
 G(1)=EE(1,5)
 G(2)=EE(2,5)+EE(1,4)
 G(3)=EE(2,4)
 E(1)=EE(1,6)
 E(2)=EE(2,6)+EE(1,2)
 E(3)=EE(2,2)
 E(4)=EE(1,1)
 E(5)=EE(2,1)+EE(1,6)
 E(6)=EE(2,6)
 F(1)=C(5,6) 
 F(2)=C(4,6)+C(2,5)
 F(3)=C(2,4) 
 F(4)=C(1,5) 
 F(5)=C(1,4)+C(5,6)
 F(6)=C(4,6) 
 F(7)=C(1,6) 
 F(8)=C(1,2)+C(6,6)
 F(9)=C(2,6) 
 do KQ=1,3
  do KR=1,3
   do KS=1,3
    KT=KQ+KR+KS-2 
    QM(KT,2)=QM(KT,2)-(G(KQ)*G(KR+3)*G(KS+6)+E(KQ)*E(KR+6)*F(KS+3)+  &
                       E(KQ+3)*F(KR+6)*F(KS)-E(KQ+3)*F(KR+3)*G(KS+3)- &
                       E(KQ)*F(KR)*G(KS+6)-E(KQ+6)*F(KR+6)*G(KS))
   end do
  end do
 end do
 G(4)=C(5,6) 
 G(5)=C(4,6)+C(2,5)
 G(6)=C(2,4) 
 E(7)=C(1,5) 
 E(8)=C(1,4)+C(5,6)
 E(9)=C(4,6) 
 F(1)=C(5,5) 
 F(2)=2.0*C(4,5) 
 F(3)=C(4,4) 
 do KQ=1,3
  do KR=1,3
   do KS=1,3
    KT=KQ+KR+KS-2 
    QM(KT,3)=QM(KT,3)+(G(KQ)*G(KR+3)*G(KS+6)+E(KQ)*E(KR+6)*F(KS+3)+  &
                    E(KQ+3)*F(KR+6)*F(KS)-E(KQ+3)*F(KR+3)*G(KS+3)- &
                    E(KQ)*F(KR)*G(KS+6)-E(KQ+6)*F(KR+6)*G(KS))
   end do
  end do
 end do
 G(7)=C(1,6) 
 G(8)=C(1,2)+C(6,6)
 G(9)=C(2,6) 
 F(7)=C(6,6) 
 F(8)=2.0*C(2,6) 
 F(9)=C(2,2) 
 F(4)=C(5,6) 
 F(5)=C(4,6)+C(2,5)
 F(6)=C(2,4) 
 do KQ=1,3
  do KR=1,3
   do KS=1,3
    KT=KQ+KR+KS-2 
    QM(KT,4)=QM(KT,4)-(G(KQ)*G(KR+3)*G(KS+6)+E(KQ)*E(KR+6)*F(KS+3)+  &
                    E(KQ+3)*F(KR+6)*F(KS)-E(KQ+3)*F(KR+3)*G(KS+3)- &
                    E(KQ)*F(KR)*G(KS+6)-E(KQ+6)*F(KR+6)*G(KS))
   end do
  end do
 end do
 HI(1)=-EN(1,1)
 HI(5)=-2.0*EN(1,2)
 HI(9)=-EN(2,2)
 HI(2)=G(1)
 HI(6)=G(2)
 HI(10)=G(3) 
 HI(3)=E(1)
 HI(7)=E(2)
 HI(11)=E(3) 
 HI(4)=E(4)
 HI(8)=E(5)
 HI(12)=E(6) 
 do J=1,4 
  MAPN%QR(1)=MAPN%QR(1)+QM(1,J)*HI(J) 
  MAPN%QR(2)=MAPN%QR(2)+QM(2,J)*HI(J)+QM(1,J)*HI(J+4) 
  MAPN%QR(3)=MAPN%QR(3)+QM(3,J)*HI(J)+QM(2,J)*HI(J+4)+QM(1,J)*HI(J+8) 
  MAPN%QR(4)=MAPN%QR(4)+QM(4,J)*HI(J)+QM(3,J)*HI(J+4)+QM(2,J)*HI(J+8) 
  MAPN%QR(5)=MAPN%QR(5)+QM(5,J)*HI(J)+QM(4,J)*HI(J+4)+QM(3,J)*HI(J+8) 
  MAPN%QR(6)=MAPN%QR(6)+QM(6,J)*HI(J)+QM(5,J)*HI(J+4)+QM(4,J)*HI(J+8) 
  MAPN%QR(7)=MAPN%QR(7)+QM(7,J)*HI(J)+QM(6,J)*HI(J+4)+QM(5,J)*HI(J+8) 
  MAPN%QR(8)=MAPN%QR(8)+QM(7,J)*HI(J+4)+QM(6,J)*HI(J+8) 
  MAPN%QR(9)=MAPN%QR(9)+QM(7,J)*HI(J+8) 
 end do
 do KP=1,9
  MAPN%QR(KP)=MAPN%QR(KP)/MAPN%QR(9) 
  if (SCALE30%LTEST.eq.1) then
   write(6,"('  QR=',F15.5)") MAPN%QR(KP)
  end if 
 end do
!*******************************************************************
!*       Find zeroes of polynomials                                *
!*******************************************************************
 MAPN%KRASH=0 
 MAPN%NEW=9 
 MAPN%ZR=0.1
 MAPN%ZI=1.0
 CALL NEWTON(MAPN) 
 if (MAPN%KRASH.gt.0) then 
  STOP 'POLYNOM(4) HAS A REAL ROOT '
 end if
 if (MAPN%KRASH.lt.0) then 
  STOP 'SOLUTION OF POLYNOMIAL DOES NOT CONVERGE '
 end if
 PR(1)=MAPN%ZR
 PI(1)=ABS(MAPN%ZI) 
 MAPN%ZI=-MAPN%ZI
 CALL NEWTON(MAPN) 
 if (MAPN%KRASH.gt.0) then 
  STOP 'POLYNOM(4) HAS A REAL ROOT '
 end if
 if (MAPN%KRASH.lt.0) then 
  STOP 'SOLUTION OF POLYNOMIAL DOES NOT CONVERGE '
 end if
 MAPN%ZR=0.5
 MAPN%ZI=0.9
 CALL NEWTON(MAPN) 
 if (MAPN%KRASH.gt.0) then 
  STOP 'POLYNOM(4) HAS A REAL ROOT '
 end if
 if (MAPN%KRASH.lt.0) then 
  STOP 'SOLUTION OF POLYNOMIAL DOES NOT CONVERGE '
 end if
 PR(2)=MAPN%ZR
 PI(2)=ABS (MAPN%ZI)
 MAPN%ZI=-MAPN%ZI
 CALL NEWTON(MAPN) 
 if (MAPN%KRASH.gt.0) then 
  STOP 'POLYNOM(4) HAS A REAL ROOT '
 end if
 if (MAPN%KRASH.lt.0) then 
  STOP 'SOLUTION OF POLYNOMIAL DOES NOT CONVERGE '
 end if
 MAPN%ZR=0.1
 MAPN%ZI=1.0
 CALL NEWTON(MAPN) 
 if (MAPN%KRASH.gt.0) then 
  STOP 'POLYNOM(4) HAS A REAL ROOT '
 end if
 if (MAPN%KRASH.lt.0) then 
  STOP 'SOLUTION OF POLYNOMIAL DOES NOT CONVERGE '
 end if
 PR(3)=MAPN%ZR
 PI(3)=ABS(MAPN%ZI) 
 MAPN%ZI=-MAPN%ZI
 CALL NEWTON(MAPN) 
 if (MAPN%KRASH.gt.0) then 
  STOP 'POLYNOM(4) HAS A REAL ROOT '
 end if
 if (MAPN%KRASH.lt.0) then 
  STOP 'SOLUTION OF POLYNOMIAL DOES NOT CONVERGE '
 end if
 MAPN%ZR=-MAPN%ZR
 CALL NEWTON(MAPN) 
 if (MAPN%KRASH.gt.0) then 
  STOP 'POLYNOM(4) HAS A REAL ROOT '
 end if
 if (MAPN%KRASH.lt.0) then 
  STOP 'SOLUTION OF POLYNOMIAL DOES NOT CONVERGE '
 end if
 PR(4)=MAPN%ZR
 PI(4)=ABS(MAPN%ZI) 
 MAPN%ZR=-C(4,5)/C(4,4) 
 MAPN%ZI=SQRT(ABS(C(4,4)*C(5,5)-C(4,5)**2))/C(4,4)
 do N=1,2
  if (((MAPN%ZR-PR(N))**2+(MAPN%ZI-PI(N))**2-(MAPN%ZR-PR(N+1))**2-(MAPN%ZI-PI(N+1))**2).lt.0.0) then 
   Z=PR(N) 
   PR(N)=PR(N+1) 
   PR(N+1)=Z 
   Z=PI(N) 
   PI(N)=PI(N+1) 
   PI(N+1)=Z 
  end if
 end do
 MP%PC=CMPLX(PR,PI)
 if (SCALE30%LTEST.eq.1) then
  write(6,"('  SOLUTIONS OF THE POLYNOM  ')") 
  write(6,"('  PC= ',2F15.8)") (MP%PC(J),J=1,4) 
 end if 
!*******************************************************************
!*        BESTIMMUNG DER EIGENVEKTOREN. VERWENDET WIRD DER         *
!*        " KOFAKTOR "ZUR BESTIMMUNG DER ELEMENTE DER EIGENVEKTOREN*
!*******************************************************************
 do K=1,4 
  A(1,1)=C(1,1)+2.0*MP%PC(K)*C(1,6)+C(6,6)*MP%PC(K)**2
  A(1,2)=C(1,6)+MP%PC(K)*(C(1,2)+C(6,6))+C(2,6)*MP%PC(K)**2 
  A(1,3)=C(1,5)+MP%PC(K)*(C(1,4)+C(5,6))+C(4,6)*MP%PC(K)**2 
  A(1,4)=EE(1,1)+MP%PC(K)*(EE(2,1)+EE(1,6))+EE(2,6)*MP%PC(K)**2 
  A(2,2)=C(6,6)+2.0*MP%PC(K)*C(2,6)+C(2,2)*MP%PC(K)**2
  A(2,3)=C(5,6)+MP%PC(K)*(C(4,6)+C(2,5))+C(2,4)*MP%PC(K)**2 
  A(2,4)=EE(1,6)+MP%PC(K)*(EE(2,6)+EE(1,2))+EE(2,2)*MP%PC(K)**2 
  A(3,3)=C(5,5)+2.0*MP%PC(K)*C(4,5)+C(4,4)*MP%PC(K)**2
  A(3,4)=EE(1,5)+MP%PC(K)*(EE(2,5)+EE(1,4))+EE(2,4)*MP%PC(K)**2 
  A(4,4)=-(EN(1,1)+2.0*MP%PC(K)*EN(1,2)+EN(2,2)*MP%PC(K)**2)
  A(2,1)=A(1,2) 
  A(3,1)=A(1,3) 
  A(3,2)=A(2,3) 
  A(4,1)=A(1,4) 
  A(4,2)=A(2,4) 
  A(4,3)=A(3,4) 
  MX=N1(K)
  MY=N2(K)
  MZ=N3(K)
  do J=1,4 
   NX=N1(J)
   NY=N2(J)
   NZ=N3(J)
   MP%AS(J,K)=((-1)**(J+K))*(A(NX,MX)*A(NY,MY)*A(NZ,MZ)+A(NX,MY)*A(NY,MZ)*A(NZ,MX)+A(NX,MZ)*A(NY,MX)*A(NZ,MY)- &
               A(NX,MZ)*A(NY,MY)*A(NZ,MX)-A(NX,MX)*A(NY,MZ)*A(NZ,MY)-A(NX,MY)*A(NY,MX)*A(NZ,MZ))
  end do
  if (SCALE30%LTEST.eq.1) then
   write(6,"('  AS--EIGENVEKTOR',8F8.4)") MP%AS(1,K),MP%AS(2,K),MP%AS(3,K),MP%AS(4,K)
   do K1=1,4 
    AY(K1)=CMPLX(0.0,0.0) 
    do K2=1,4 
     AY(K1)=AY(K1)+A(K1,K2)*MP%AS(K2,K) 
    end do
   end do
   write(6,"('  PROBE--',F15.8)") (AY(K3),K3=1,4)
  endif 
 end do 
 do I=1,4 
  X=0.0 
  do J=1,4 
   X=X+MP%AS(J,I)*CONJG(MP%AS(J,I))
  end do
  X=SQRT(X) 
  do J=1,4
   MP%AS(J,I)=MP%AS(J,I)/X 
  end do
 end do
end if  ! (if (goto1000.eqv..FALSE.))
! original code has a continue statement here ...
!1000 CONTINUE
!*******************************************************************
!*        BERECHNUNG DER VEKTOREN  L(I,K)                          *
!*        NACH BARNETT UND LOTHE                                   *
!*******************************************************************
 do J=1,3 
  NJ=NN(J)
  do K=1,4 
   MP%EL(J,K)=CMPLX(0.0,0.0)
   do L=1,3 
    NL=NN(L)
    ML=MM(L)
    MP%EL(J,K)=MP%EL(J,K)+(C(NJ,ML)+MP%PC(K)*C(NJ,NL))*MP%AS(L,K) 
   end do
   MP%EL(J,K)=-(MP%EL(J,K)+(EE(1,NJ)+MP%PC(K)*EE(2,NJ))*MP%AS(4,K))
  end do
 end do
 do K=1,4 
  MP%EL(4,K)=CMPLX(0.0,0.0)
  do L=1,3 
   NL=NN(L)
   ML=MM(L)
   MP%EL(4,K)=MP%EL(4,K)+(EE(2,ML)+MP%PC(K)*EE(2,NL))*MP%AS(L,K) 
  end do
  MP%EL(4,K)=-(MP%EL(4,K)-(EN(2,1)+MP%PC(K)*EN(2,2))*MP%AS(4,K))
 end do
!*******************************************************************
!*      NORMIERUNG     2*A(I,K)*L(J,K)=1                           *
!*******************************************************************
 do K=1,4 
  MXX(K)=CMPLX(0.0,0.0) 
  do L=1,4 
   MXX(K)=MXX(K)+2.0*MP%AS(L,K)*MP%EL(L,K) 
  end do
  IF (MXX(K).EQ.CMPLX(0.0,0.0)) MXX(K)=CMPLX(1.0,0.0) 
  do J1=1,4
   MP%EL(J1,K)=MP%EL(J1,K)/CSQRT(MXX(K)) 
   MP%AS(J1,K)=MP%AS(J1,K)/CSQRT(MXX(K)) 
  end do
 end do
 if (SCALE30%LTEST.eq.1) then
  write(6,"('  AS-LOTHE',4F15.8)") ((MP%AS(I,J),J=1,4),I=1,4)
  write(6,"('  EL-LOTHE',4F15.8)")  ((MP%EL(I,J),J=1,4),I=1,4) 
  do K=1,4
   do L=1,4
    AXA(K,L)=CMPLX(0.0,0.0) 
    AXL(K,L)=CMPLX(0.0,0.0) 
    LXL(K,L)=CMPLX(0.0,0.0) 
    do I=1,4
     AXA(K,L)=AXA(K,L)+MP%AS(K,I)*MP%AS(L,I) 
     LXL(K,L)=LXL(K,L)+MP%EL(K,I)*MP%EL(L,I) 
     AXL(K,L)=AXL(K,L)+MP%AS(K,I)*MP%EL(L,I) 
    end do
    AXA(K,L)=2.0*REAL(AXA(K,L)) 
    LXL(K,L)=2.0*REAL(LXL(K,L)) 
    AXL(K,L)=2.0*REAL(AXL(K,L)) 
   end do
  end do
  write(6,"('  AXA',4F20.8)") ((AXA(I1,I2),I2=1,4),I1=1,4) 
  write(6,"('  LXL',4F20.8)") ((LXL(I1,I2),I2=1,4),I1=1,4) 
  write(6,"('  AXL',4F20.8)") ((AXL(I1,I2),I2=1,4),I1=1,4) 
  PR1=PR(1)**2
  PI1=PI(1)**2
  PR2=PR(2)**2
  PI2=PI(2)**2
  PR3=PR(3)**2
  PI3=PI(3)**2
  PR4=PR(4)**2
  PI4=PI(4)**2
  AZ=-2.0*(PR(1)+PR(2)) 
  BY=PR1+PI1+PR2+PI2+4.0*PR(1)*PR(2)
  CY=-2.0*((PR1+PI1)*PR(2)+(PR2+PI2)*PR(1)) 
  D=(PR1+PI1)*(PR2+PI2) 
  AZZ=-2.0*(PR(3)+PR(4))
  BZ=PR3+PI3+PR4+PI4+4.0*PR(3)*PR(4)
  CZ=-2.0*((PR3+PI3)*PR(4)+(PR4+PI4)*PR(3)) 
  DZ=(PR3+PI3)*(PR4+PI4)
  MAPN%QR(9)=1.0 
  MAPN%QR(8)=AZ+AZZ
  MAPN%QR(7)=BY+AZ*AZZ+BZ
  MAPN%QR(6)=CY+AZZ*BY+AZ*BZ+CZ
  MAPN%QR(5)=D+AZZ*CY+BZ*BY+AZ*CZ+DZ 
  MAPN%QR(4)=AZZ*D+BZ*CY+BY*CZ+AZ*DZ 
  MAPN%QR(3)=BZ*D+CY*CZ+BY*DZ
  MAPN%QR(2)=D*CZ+CY*DZ
  MAPN%QR(1)=D*DZ
  write(6,"('  KONTROLLE DER WURZELN DES POLYNOMS DURCH BERECHNUNG DER POLYNOMSKOEFFIZIENTEN ') ")
  write(6,"(3F10.6)") (MAPN%QR(K),K=1,9)
 end if 
end subroutine PANCALC

!--------------------------------------------------------------------------
recursive subroutine RKM(MRD)
!DEC$ ATTRIBUTES DLLEXPORT :: RKM
!! author: Head, A.K. and Humble, P. and Clarebrough, L.M. and Morton, A.J. and Forwood, C.T. 
!! version: 1.0 
!! date: 08/13/19 adapted from .f77 original
!!
!! Runge-Kutta integration of Howie-Whelan equations
!!
!! Note:
!! The original RKM routine was rather complicated due to the extensive
!! use of computed goto statements.  The new routine is still ugly, but 
!! no longer uses any goto statements, computed or otherwise...
!! 

IMPLICIT NONE

type(MRD_block),INTENT(INOUT)     :: MRD
!f2py intent(in,out) ::  MRD

integer(kind=irg)                 :: M1, M2, J, M, leave
real(kind=sgl)                    :: ERHIGH, ERLOW, H1, H2, H3, XT, TEST
logical                           :: verbose = .FALSE.

M1=4 
MRD%KOUNT=0  
ERHIGH=5.0*MRD%ERROR 
ERLOW=0.03125*ERHIGH  
J = 1
do while (J.ne.100) 
 select case(J)
  case(1)
    if (verbose) write (*,*) 'case 1'
    if (MRD%Q.ne.0.0) then
      J = 2
    else
      J = 20
    end if

  case(2)
    if (verbose) write (*,*) 'case 2'
    IF ( ( (MRD%X1-MRD%X)/MRD%Q - 1.0000001 ).le.0.0 ) then ! 29,29,7
      J=29
    else 
      J=7
    end if  

  case(3)
    if (verbose) write (*,*) 'case 3'
    MRD%Q=H1; J=11

  case(4)
    if (verbose) write (*,*) 'case 4'
    MRD%Q=2.0*MRD%Q; H2=0.5*MRD%Q; H3=MRD%Q/3.0; J=2

  case(5)
    if (verbose) write (*,*) 'case 5'
    IF ((H1-2.0*MRD%Q).lt.0.0) then ! 30,3,3
      J = 30
    else 
      J = 3
    end if

  case(6)
    if (verbose) write (*,*) 'case 6'
    M1=5; H2=0.5*MRD%Q; H3=MRD%Q/3.0; J=7

  case(7)
    if (verbose) write (*,*) 'case 7'
    XT=MRD%X; MRD%YT(:)=MRD%Y(:); J=8

  case(8)
    if (verbose) write (*,*) 'case 8 '
    CALL DERIV(MRD)
    MRD%DT(:,1)=H3*MRD%D(:); MRD%Y(:)=MRD%Y(:)+MRD%DT(:,1); MRD%X=MRD%X+H3
    CALL DERIV(MRD)
    MRD%Y(:)=MRD%YT(:)+0.5*(MRD%DT(:,1)+H3*MRD%D(:)); MRD%SKIP=1.0
    CALL DERIV(MRD)
    MRD%SKIP=0.0; MRD%DT(:,2)=MRD%Q*MRD%D(:); MRD%Y(:)=MRD%YT(:)+0.375*(MRD%DT(:,1)+MRD%DT(:,2)); MRD%X=XT+H2
    CALL DERIV(MRD)
    MRD%DT(:,3)=4.0*H3*MRD%D(:); MRD%Y(:)=MRD%YT(:)+1.5*(MRD%DT(:,1)+MRD%DT(:,3)-MRD%DT(:,2)); MRD%X=XT+MRD%Q
    CALL DERIV(MRD)
    M2=0
    J = 10

  case(9)
    if (verbose) write (*,*) 'case 9'
    M1=4; MRD%Q=0.5*MRD%Q; H2=0.5*MRD%Q; H3=MRD%Q/3.0; MRD%Y(:)=MRD%YT(:); MRD%X=XT; J=8

  case(10)
    if (verbose) write (*,*) 'case 10'
    leave = 0
    doten: do M=1,8 
      MRD%DT(M,4)=H3*MRD%D(M)
      TEST=ABS (MRD%DT(M,1)+MRD%DT(M,3)-0.5*(MRD%DT(M,4)+3.0*MRD%DT(M,2)))
      IF ((TEST-ERHIGH).lt.0.0) then ! 26,9,9
        IF ((TEST-ERLOW).lt.0.0) then ! 10,27,27
         CYCLE doten
        END IF 
        M2=-2     ! old line 27
      ELSE 
        J = 9
        leave = 1
        EXIT doten     ! jump out of case select
      END IF
    end do doten  ! old line 10
    if (leave.eq.0) then 
      MRD%Y(:)=MRD%YT(:)+0.5*(MRD%DT(:,1)+MRD%DT(:,3)+MRD%DT(:,4))
      MRD%KOUNT=MRD%KOUNT+1
      J=M1+M2
    end if

  case(11)
    if (verbose) write (*,*) 'case 11'
    H2=0.5*MRD%Q; H3=MRD%Q/3.0; J=100

  case(20)
    if (verbose) write (*,*) 'case 20'
    MRD%Q=MRD%X1-MRD%X; H1=MRD%Q; J=6 

  case(29)
    if (verbose) write (*,*) 'case 29'
    H1=MRD%Q; MRD%Q=MRD%X1-MRD%X; J=6

  case(30)
    if (verbose) write (*,*) 'case 30'
    MRD%Q=2.0*MRD%Q; J=11

  case default 
    write (*,*) 'RKM: case value does not exist'
    stop

 end select
end do 
! if J=100, then we return from this routine     

end subroutine RKM

!--------------------------------------------------------------------------
recursive subroutine DERIV(MRD)
!DEC$ ATTRIBUTES DLLEXPORT :: DERIV
!! author: Head, A.K. and Humble, P. and Clarebrough, L.M. and Morton, A.J. and Forwood, C.T. 
!! version: 1.0 
!! date: 08/13/19 adapted from .f77 original
!!
!! Derivative calculation
!!
! 
!************************************************************** 
!*     SUBROUTINE DERIV   BERECHNUNG DER ABLEITUNG VON        * 
!*     G IN R NACH DER FORMEL VON STROH, ERWEITERT FUER       * 
!*     PIEZOELEKTRISCHE KRISTALLE. BERUECKSICHTIGUNG DER      * 
!*     ANOMALEN ABSORPTION                                    * 
!************************************************************** 

IMPLICIT NONE

type(MRD_block),INTENT(INOUT)     :: MRD
!f2py intent(in,out) ::  MRD

real       :: X11, X22, X33, X44, R1, R2, R3, R4, BETA1, BETA2, BETA3, BETA4, Z
real,save  :: BETA
! 
! 
 if (MRD%SKIP.eq.0.0) then
!*************************************
!*     ERSTE VERSETZUNG              *
!*************************************
  if ((MRD%X-MRD%CN(21)).eq.0.0) then 
   X11=1.0E-10 
  else  
   X11=MRD%X-MRD%CN(21)
  end if
  R1=(MRD%CN(30)-MRD%CN(31))/X11
!*************************************
!*     ZWEITE VERSETZUNG             *
!*************************************
  if ((MRD%X+MRD%CN(21)).eq.0.0) then 
   X22=1.0E-10 
  else  
   X22=MRD%X+MRD%CN(21)
  end if
  R2=(MRD%CN(30)+MRD%CN(31))/X22
!*************************************
!*     DRITTE VERSETZUNG             *
!*************************************
  if ((MRD%X+MRD%CN(43)).eq.0.0) then
   X33=1.0E-10 
  else 
   X33= MRD%X+MRD%CN(43) 
  end if
  R3=(MRD%CN(30)+MRD%CN(44))/X33
!*************************************
!*     VIERTE VERSETZUNG             *
!*************************************
  if ((MRD%X-MRD%CN(43)).eq.0.0) then  
   X44=1.0E-10 
  else 
   X44=MRD%X-MRD%CN(43)
  end if
  R4=(MRD%CN(30)-MRD%CN(44))/X44
  BETA1=(((R1*MRD%CN(1)+MRD%CN(5))/((R1+MRD%CN(9))**2+MRD%CN(13)))+&
         ((R1*MRD%CN(2)+MRD%CN(6))/((R1+MRD%CN(10))**2+MRD%CN(14)))+&
         ((R1*MRD%CN(3)+MRD%CN(7))/((R1+MRD%CN(11))**2+MRD%CN(15)))+&
         ((R1*MRD%CN(4)+MRD%CN(8))/((R1+MRD%CN(12))**2+MRD%CN(16))))/X11 
  BETA2=(((R2*MRD%CN(22)+MRD%CN(26))/((R2+MRD%CN(9))**2+MRD%CN(13)))+&
         ((R2*MRD%CN(23)+MRD%CN(27))/((R2+MRD%CN(10))**2+MRD%CN(14)))+&
         ((R2*MRD%CN(24)+MRD%CN(28))/((R2+MRD%CN(11))**2+MRD%CN(15)))+&
         ((R2*MRD%CN(25)+MRD%CN(29))/((R2+MRD%CN(12))**2+MRD%CN(16))))/X22 
  BETA3=(((R3*MRD%CN(32)+MRD%CN(36))/((R3+MRD%CN(9))**2+MRD%CN(13)))+&
         ((R3*MRD%CN(33)+MRD%CN(37))/((R3+MRD%CN(10))**2+MRD%CN(14)))+&
         ((R3*MRD%CN(34)+MRD%CN(38))/((R3+MRD%CN(11))**2+MRD%CN(15)))+&
         ((R3*MRD%CN(35)+MRD%CN(39))/((R3+MRD%CN(12))**2+MRD%CN(16))))/X33 
  BETA4=(((R4*MRD%CN(51)+MRD%CN(55))/((R4+MRD%CN(9))**2+MRD%CN(13)))+&
         ((R4*MRD%CN(52)+MRD%CN(56))/((R4+MRD%CN(10))**2+MRD%CN(14)))+&
         ((R4*MRD%CN(53)+MRD%CN(57))/((R4+MRD%CN(11))**2+MRD%CN(15)))+&
         ((R4*MRD%CN(54)+MRD%CN(58))/((R4+MRD%CN(12))**2+MRD%CN(16))))/X44 
  BETA=MRD%CN(18)+BETA1+BETA2+BETA3+BETA4 
end if
!*******************************************************************
!*     BERUECKSICHTIGUNG DER ANOMALEN ABSORPTION                   *
!*******************************************************************
  Z=MRD%ANO*(MRD%Y(1)+MRD%Y(3)) 
  MRD%D(1)=Z-MRD%Y(4) 
  MRD%D(3)=-BETA*MRD%Y(4)+Z-MRD%Y(2)
  Z=MRD%ANO*(MRD%Y(2)+MRD%Y(4)) 
  MRD%D(2)=Z+MRD%Y(3) 
  MRD%D(4)=BETA*MRD%Y(3)+Z+MRD%Y(1) 
  Z=MRD%ANO*(MRD%Y(5)+MRD%Y(7)) 
  MRD%D(5)=Z-MRD%Y(8) 
  MRD%D(7)=-BETA*MRD%Y(8)+Z-MRD%Y(6)
  Z=MRD%ANO*(MRD%Y(6)+MRD%Y(8)) 
  MRD%D(6)=Z+MRD%Y(7) 
  MRD%D(8)=BETA*MRD%Y(7)+Z+MRD%Y(5) 
end subroutine DERIV

!--------------------------------------------------------------------------
subroutine dohh_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: dohh_
!! author: MDG 
!! version: 1.0 
!! date: 01/25/22
!!
!! perform the computations

!
! PROGRAM: EMhh4
!
!> @author Original hh.f77 code history, see
!> Author = {Head, A.K. and Humble, P. and Clarebrough, L.M. and Morton, A.J. and Forwood, C.T.},
!> Publisher = {North Holland Publishing Company},
!> Series = {Defects in {C}rystalline {S}olids, edited by Amelinckx, S. and Gevers, R. and Nihoul, J.},
!> Title = {Computed Electron Micrographs and Defect Identification},
!> Volume = {7},
!> Year = 1973}
!> Translated to fortran-90 by MDG. All original variable names are kept in upper case
!> notation, all EMsoft add-ons and new code is in lower case.
!
!> @note improvements to the code:
!> - Integrated with EMsoft modules.
!> - Added HDF5 and image output + support for tilt series.
!> - Added standard EMsoft namelist handling to replace the older input format
!> - Added class structure for EMsoftOO
!
! The new f90 code is a bit more readable than the original f77 code,
! mostly because the modern Fortran language has better control structures.
! All ordinary and computed goto statements have been replaced with
! case and other control statements.  No line labels are used at all.
! Variable names have been left unchanged (mostly) to facilitate comparison
! with original hh4 code as well as debugging.  The original f77 code (a version
! generated in the '80s by the group of Prof. Skalicky at the Technical University
! of Vienna, Austria.) can be downloaded from:
!
! http://ctem.web.cmu.edu/software/HeadHumble.tar.gz
! 
!> @date 01/15/03 MDG 1.0 original conversion from f77 version
!> @date 08/13/19 MDG 2.0 new version, integrated with EMsoft 4.3
!> @date 08/23/19 MDG 2.1 debugged and tested; added HDF5 and image output formats
!> @date 01/25/22 MDG 3.0 converted to object oriented fortran 2018
!--------------------------------------------------------------------------
! original program header 
!** (INPUT,OUTPUT,TAPE5=INPUT,TAPE6=OUTPUT,TAPE10) 
!***********************************************************************
!*           HEADHUMBLE 4 - PIEZOELEKTRISCHE KRISTALLE                 *
!*           PROGRAMM FUER 4 SYMMETRISCHE VERSETZUNGEN MIT             *
!*           3 STAPELFEHLERN WOBEI FUER DIE FAULTPLANES                *
!*           GILT: FP1 = FP3. KONFIGURATION SYMMETRISCH ZUM            *
!*           URSPRUNG. KEIN D - UND M - DRUCK                          *
!*            30. AUGUST 1983                                          *
!***********************************************************************

use mod_EMsoft
use mod_crystallography
use mod_symmetry
use mod_diffraction
use mod_io
use mod_timing
use HDF5
use mod_HDFsupport
use mod_HDFnames
use stringconstants
use mod_image
use mod_memory
use, intrinsic :: iso_fortran_env

IMPLICIT NONE 

class(hh_T), INTENT(INOUT)              :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 

type(Cell_T)                            :: cell
type(SpaceGroup_T)                      :: SG
type(HDF_T)                             :: HDF
type(HDFnames_T)                        :: HDFnames
type(Timing_T)                          :: timer
type(IO_T)                              :: Message
type(memory_T)                          :: mem
type(Diffraction_T)                     :: Diff

! all original COMMON blocks are replaced by user-defined structures in this module
type(MAPN_block)              :: MAPN
type(MA_block)                :: MA
type(MKAP_block)              :: MKAP
type(MRD_block)               :: MRD
type(MT_block)                :: MT
type(MKT_block)               :: MKT
type(SCALE30_block)           :: SCALE30
type(MP_block)                :: MP
type(MAP_block)               :: MAP

!=======================
! original hh4.f variables 
!=======================
! regular integers
integer(kind=irg)             :: LLQ, NNN, NNNN, I, J, JB, JC, K, L, KMIN, KMAX, KTOT, MOVE, LUCK, ISTORE, LSWITC, &
                                 IFLAG, JT, JM, KOUNTF, INDL, KK, JZ, IND, LD, LQ, LPIEZO
! regular integer arrays
integer(kind=irg)             :: ITYPE(4)
! integer constants
integer(kind=irg),parameter   :: NP(3) = (/2,3,1/), NQ(3) = (/3,1,2/), addon(4) = (/ 1, 0, 2, 3 /)
integer(kind=irg)             :: ICOL, IROW, ICOLP 
! integer allocatable
integer(kind=irg),allocatable :: IX(:), IXX(:)

! regular reals
real(kind=sgl)                :: GINB, GINBXU, Z, SBM1, SBM2, SBM3, SFN1, SFN2, SFN3, FNBM, PT, SL, PT2, SL2, &
                                 THBM, EXT1, EXT2, EXT3, EXT4, EXTRA, FRACTI, DIVISO, DELT, WL, DELW, DELL, BACK, BACKD, &
                                 STORE, DELTA, DEL, DEL2, DISTR, DISTRA, DISTL, DISTLA, XXX, YYY, ZZZ, VVV, FAULT1, &
                                 ALPHA, FAULT2, FAULT3, STARTA, SURFAC, XX1, ANO, FINISH, START, THICK, XIGEE, &
                                 TRAMP3, TRAMP7,DNR, DNI, DNN, TTB, TTD, WW, LC1, LC2, LC3, LC4
! regular real arrays
real(kind=sgl)                :: GD(3), BD(4), B2D(4), B3D(4), B4D(4), BM(3), FN(3), FP1X(3),  &
                                 FPX(3), FP3X(3), FP(3), FP3(3), FNX(3), DCX(3,3), DR(4), DI(4), &
                                 UR(4,4), UI(4,4), VR(4,4), VI(4,4), DD(3), SUR(4), SUI(4), COSA(3), SINA(3), &
                                 UX(3), AB(3), AB1(3), POSA(4), POSB(4), COORD(4), HANDL(4), &
                                 HANDR(4), TEMPY(8), QL1(4), QL2(4), QL3(4), QL4(4), S(4,4), QS(4,4), RR(3)
! real allocatable 
real(kind=sgl),allocatable    :: FX(:,:), TBD(:,:), TQB(:),TQD(:), BFINTENS(:,:,:), DFINTENS(:,:,:)

! regular complex
complex(kind=sgl)             :: SU(4), CNX(8), MXXX(4,4), MYYY(4,4), MZZZ(4,4)

! character variables
character(15)                 :: IY
character(fnlen)              :: diagfile = 'HHdiagnostics.txt'

!=======================
! additional EMsoft variables (not originally in hh4.f code)
!=======================
type(gnode)                   :: rlp
character(fnlen)              :: mess, fname
character(fnlen),allocatable  :: legendfiles(:)
character(3)                  :: filenum
character(11)                 :: dstr
character(15)                 :: tstrb
character(15)                 :: tstre
integer(kind=irg)             :: numim, imnum, hdferr
real(kind=sgl),allocatable    :: wvalues(:)
logical                       :: isallowed
real(kind=dbl)                :: eps = 1.0D-6
real(kind=sgl)                :: wstep, io_real(1), mimi, mama

! declare variables for use in object oriented image module
integer                       :: iostat
character(len=128)            :: iomsg
logical                       :: isInteger
type(image_t)                 :: im
integer(int8), allocatable    :: output_image(:,:)

! open the HDF interface
 call openFortranHDFInterface()

! construct the HDF and HDFnames classes
 HDF = HDF_T()
 HDFnames = HDFnames_T()

 associate( hhnl => self%nml )

! initialize the timing routines
 timer = Timing_T()
 tstrb = timer%getTimeString()
 dstr = timer%getDateString()
 tstre = ''
 call timer%Time_tick(1)

! set the crystal structure filename
 cell = cell_T()
 call cell%setFileName(hhnl%xtalname)

! get some parameters from the namelist
 ICOL = hhnl%ICOL 
 IROW = hhnl%IROW
 ICOLP = hhnl%ICOL+1 
 numim = hhnl%wnum

! allocate the string array for the legend files 
 mem = memory_T()
 call mem%alloc( legendfiles, (/ numim /), 'legendfiles' )
 do i=1,numim
   write (filenum,"(I3.3)") i
   legendfiles(i) = 'legend_'//filenum//'.txt'
 end do

! allocate all allocatable arrays here 
 call mem%alloc( IX, (/ ICOLP /), 'IX') 
 call mem%alloc( IXX, (/ ICOLP /), 'IXX')
 call mem%alloc( FX, (/ ICOL,4 /), 'FX') 
 call mem%alloc( TBD, (/ IROW,ICOLP /), 'TBD') 
 call mem%alloc( TQB, (/ ICOLP /), 'TQB') 
 call mem%alloc( TQD, (/ ICOLP /), 'TQD')
 call mem%alloc( BFINTENS, (/ ICOL, IROW, numim /), 'BFINTENS')
 call mem%alloc( DFINTENS, (/ ICOL, IROW, numim /), 'DFINTENS')
 call mem%alloc( wvalues, (/ numim /), 'wvalues' )

! set original hh4.f parameters that are not used in this f90 version of the code
! (we're setting these out of sheer nostalgia for the old f77 code ...)
 IND = 0
 LPIEZO = 0
 LQ=2
 LLQ=1 
 NNN=ICOL
 NNNN=30+10*LQ 

! copy namelist entries to appropriate variables

! elastic moduli
 MKAP%D1(1,1:6) = hhnl%D1row1 
 MKAP%D1(2,1:6) = hhnl%D1row2 
 MKAP%D1(3,1:6) = hhnl%D1row3 
 MKAP%D1(4,1:6) = hhnl%D1row4 
 MKAP%D1(5,1:6) = hhnl%D1row5 
 MKAP%D1(6,1:6) = hhnl%D1row6 
! other parameters
 MT%LU = hhnl%LU
 MT%LG = hhnl%LG
 MT%LBM = hhnl%LBM
 MT%LFN = hhnl%LFN
 MT%LB = hhnl%LB
 MT%LB2 = hhnl%LB2
 MT%LB3 = hhnl%LB3
 MT%LB4 = hhnl%LB4
 MT%LD = hhnl%LD
 MT%LD2 = hhnl%LD2
 MT%LD3 = hhnl%LD3
 MT%LD4 = hhnl%LD4
 MT%LFP1 = hhnl%LFP1
 MT%LFP = hhnl%LFP
 MT%LFP3 = hhnl%LFP3 
 MT%LS1 = hhnl%LS1
 MT%LS2 = hhnl%LS2
 MT%LS3 = hhnl%LS3
 MT%LQ1 = hhnl%LQ1
 MT%LQ2 = hhnl%LQ2
 MT%LQ3 = hhnl%LQ3
 ! MT%LF1 = hhnl%LF1
 ! MT%LF2 = hhnl%LF2
 ! MT%LF3 = hhnl%LF3
 ! MT%LF4 = hhnl%LF4
 SCALE30%LTEST = hhnl%LTEST
 LD = hhnl%LD 

! we don't do piezoelectric contributions in this version of the program
 MKAP%EP = 0.0
 MKAP%EA = 0.0

! open diagnostic file for output if hhnl%LTEST.eq.1
 if (hhnl%LTEST.eq.1) then 
   OPEN(dataunit,FILE=trim(diagfile), status='unknown')
 end if

!*********************************************************************
!*      Select the crystal system                                    *
!*      for trigonal we can use both trigonal and hexagonal indices  *
!*      (0=trigonal) (1=hexagonal)                                   *
!
! the original hh.f program called an individual routine for each crystal system
! and essentially computed the direct and reciprocal structure matrices.
! In terms of the crystallographic variables of EMsoft, the following
! relations hold:
!   hh.f        EMsoft
!   AT          cell%dsm/cell%a
!   ATR         cell%rsm*cell%a
! We will mostly stick to the hh.f variable names for easy translation.
! The routine CrystalData does everything that was originally done
! by the following hh.f subroutines:  TRICLIN, MONOCLI, RHOMBIS,
! TRIGONA, TETRAGONA, HEXAGON, and CUBIC.
!*********************************************************************
 call cell%getCrystalData(hhnl%xtalname, SG, EMsoft, verbose=.TRUE.)
 ! scale the direct and reciprocal structure matrices
 MKT%AT  = cell%getdsm()/cell%getLatParm('a')
 MKT%ATR = cell%getrsm()*cell%getLatParm('a')

!*********************************************************************
! transform all crystal variables to Cartesian frame (former TRAFO routine)
!*********************************************************************
 call cell%TransSpace(float(MT%LU),  MT%TLU,  'd', 'c')
 call cell%TransSpace(float(MT%LF1), MT%TLF1, 'd', 'c')
 call cell%TransSpace(float(MT%LF2), MT%TLF2, 'd', 'c')
 call cell%TransSpace(float(MT%LF3), MT%TLF3, 'd', 'c')
 call cell%TransSpace(float(MT%LF4), MT%TLF4, 'd', 'c')
 call cell%TransSpace(float(MT%LBM), MT%TLBM, 'd', 'c')
 call cell%TransSpace(float(MT%LFN), MT%TLFN, 'd', 'c')
 call cell%TransSpace(float(MT%LB),  MT%TLB,  'd', 'c')
 call cell%TransSpace(float(MT%LB2), MT%TLB2, 'd', 'c')
 call cell%TransSpace(float(MT%LB3), MT%TLB3, 'd', 'c')
 call cell%TransSpace(float(MT%LB4), MT%TLB4, 'd', 'c')
 call cell%TransSpace(float(MT%LS1), MT%TLS1, 'd', 'c')
 call cell%TransSpace(float(MT%LS2), MT%TLS2, 'd', 'c')
 call cell%TransSpace(float(MT%LS3), MT%TLS3, 'd', 'c')

 call cell%TransSpace(float(MT%LG),   MT%TLG,   'r', 'c')
 call cell%TransSpace(float(MT%LFP),  MT%TLFP,  'r', 'c')
 call cell%TransSpace(float(MT%LFP1), MT%TLFP1, 'r', 'c')
 call cell%TransSpace(float(MT%LFP3), MT%TLFP3, 'r', 'c')

!***********************************************************
!*      Computation of  G.B   and  G.(B X U)               *
!***********************************************************
! this uses EMsoft routines rather than the original code
 GINB=cell%CalcDot(MT%TLG,MT%TLB,'c')/FLOAT(LD) 
 call cell%CalcCross(MT%TLB,MT%TLU,RR,'c','c',0)
 GINBXU=cell%CalcDot(MT%TLG,RR,'c')/FLOAT(LD)

!***********************************************************************
! new code: compute the extinction distance and the ano ratio
! for the selected g-vector
!***********************************************************************
call Diff%setrlpmethod('WK')
call Diff%setV(dble(hhnl%kV))
call Diff%CalcWaveLength(cell)
call cell%CalcPositions(SG,'v')
call Diff%CalcUcg(cell,MT%LG,applyqgshift=.TRUE.)

isallowed = SG%IsGAllowed(MT%LG)
rlp = Diff%getrlp()
! check whether this is a lattice extinction, a symmetry extinction, or an allowed reflection
if ((rlp%Umod.lt.eps).and.(isallowed.eqv..FALSE.)) then
  call Message%printError('HHComputeImages','This reflection is absent due to lattice centering.')
end if
if ((rlp%Umod.lt.eps).and.(isallowed.eqv..TRUE.)) then
  call Message%printError('HHComputeImages','This reflection is absent due to glide or screw symmetry elements.')
end if
MRD%CN = 0.0

! get the extinction distance and the absorption ratio for this reflection
XIGEE = 10.0*rlp%xg
MRD%ANO = -1.0/rlp%ar
io_real(1) = XIGEE
call Message%WriteValue(' Extinction distance [nm] ',io_real, 1)
io_real(1) = MRD%ANO
call Message%WriteValue(' Absorption ratio         ',io_real, 1)

!***********************************************************************
!*       Make sure that the input data make geometric sense            *
!*       and initialize some default values; a lot of this was         *
!*       typical spaghetti-code....                                    *
!***********************************************************************
! make sure the START and FINISH values make sense
 FINISH = hhnl%FINISH 
 START = hhnl%START 
 THICK = hhnl%THICK 
 if ((FINISH.eq.0.0).AND.(START.eq.0.0)) then 
  START=0.0 
  FINISH=THICK
 end if
 if (FINISH.le.START) then
   call Message%printError('HHComputeImages','START after FINISH')
 end if

! intercept zero denominators by setting them to 1
 if (MT%LD.eq.0)  MT%LD=1
 if (MT%LD2.eq.0)  MT%LD2=1
 if (MT%LD3.eq.0)  MT%LD3=1
 if (MT%LD4.eq.0)  MT%LD4=1
 if (MT%LQ1.eq.0)  MT%LQ1=1
 if (MT%LQ2.eq.0)  MT%LQ2=1
 if (MT%LQ3.eq.0)  MT%LQ3=1

!***********************************************************************
!***     Default foil normal LFN = LBM                               ***
!***********************************************************************
 if (sum(MT%TLFN**2).eq.0.0) then 
  MT%TLFN=MT%TLBM
 end if

!***********************************************************************
!***     Make sure that FP1 and FP3 are identical                    ***
!***********************************************************************
 if (sum(MT%TLFP3**2).eq.0.0) then 
  MT%TLFP3=MT%TLFP1 
 else
  if (sum((MT%TLFP3-MT%TLFP1)**2).ne.0.0) then 
    call Message%printMessage('Fault planes 1 and 3 not identical in input; identity 3=1 imposed')
  end if
  MT%TLFP3=MT%TLFP1 
 end if

!***********************************************************************
!*       Create and normalize the dislocation reference frame (DC)     *
!***********************************************************************
 MAP%DC(3,1:3)=MT%TLU(1:3)
 do J=1,3 
  K=NP(J) 
  L=NQ(J) 
  MAP%DC(1,J)=MT%TLBM(K)*MT%TLU(L)-MT%TLBM(L)*MT%TLU(K) 
 end do
 do J=1,3 
  K=NP(J) 
  L=NQ(J) 
  MAP%DC(2,J)=MT%TLU(K)*MAP%DC(1,L)-MT%TLU(L)*MAP%DC(1,K) 
 end do
 do J=1,3 
  DD(J)=0.0 
  do K=1,3 
   DD(J)=DD(J)+MAP%DC(J,K)**2
  end do
  if ((DD(J)-0.0001).le.0.0) then
   call Message%printError('HHComputeImages', ' Beam parallel to line direction ')
  end if
 end do
 if ((sum(MT%LU*MT%LFP1)**2+sum(MT%LU*MT%LFP)**2+sum(MT%LU*MT%LFP3)**2).ne.0.0) then
  call Message%printError('HHComputeImages',' Line direction not in fault planes ')
 end if
 do J=1,3 
  MAP%DC(J,1:3)=MAP%DC(J,1:3)/SQRT(DD(J)) 
 end do

!********************************************************************** 
!***   Create and normalize the reference frame attached to the    **** 
!***   incident beam direction (DCX)                               **** 
!********************************************************************** 
 do J=1,3 
  DCX(1,J)=-MAP%DC(1,J) 
  DCX(2,J)=-MT%TLBM(J) 
 end do
 DCX(3,1)=DCX(1,2)*DCX(2,3)-DCX(2,2)*DCX(1,3)
 DCX(3,2)=DCX(2,1)*DCX(1,3)-DCX(1,1)*DCX(2,3)
 DCX(3,3)=DCX(1,1)*DCX(2,2)-DCX(2,1)*DCX(1,2)
 do J=1,3 
  DD(J)=sum(DCX(J,1:3)**2) 
 end do
 do J=1,3 
  DCX(J,1:3)=DCX(J,1:3)/SQRT(DD(J)) 
 end do

! 
!********************************************************************** 
!***   Transformation of input data from crystal reference frame   **** 
!***   to DC and DCX reference frames                              **** 
!********************************************************************** 
! 
 BD=0.0; QL1=0.0; QL2=0.0; QL3=0.0; QL4=0.0
 GD=0.0 ; BM=0.0 ; FN=0.0 ; FNX=0.0; B2D=0.0
 B3D=0.0; B4D=0.0; FPX=0.0; UX=0.0 ; FP1X=0.0 
 FP3X=0.0 ; FP=0.0 ; FP3=0.0
 do J=1,3 
  BD(J)=BD(J)+sum((MT%TLB(1:3))*MAP%DC(J,1:3))/FLOAT(MT%LD)
  QL1(J)=QL1(J)+sum(MT%TLF1(1:3)*MAP%DC(J,1:3))*LC1 
  QL2(J)=QL2(J)+sum(MT%TLF2(1:3)*MAP%DC(J,1:3))*LC2 
  QL3(J)=QL3(J)+sum(MT%TLF3(1:3)*MAP%DC(J,1:3))*LC3 
  QL4(J)=QL4(J)+sum(MT%TLF4(1:3)*MAP%DC(J,1:3))*LC4 
  BM(J)=BM(J)+sum(MT%TLBM(1:3)*MAP%DC(J,1:3))
  FN(J)=FN(J)+sum(MT%TLFN(1:3)*MAP%DC(J,1:3))
  FNX(J)=FNX(J)+sum(MT%TLFN(1:3)*DCX(J,1:3))
  FPX(J)=FPX(J)+sum(MT%TLFP(1:3)*DCX(J,1:3))
  FP(J)=FP(J)+sum(MT%TLFP(1:3)*MAP%DC(J,1:3))
  FP3(J)=FP3(J)+sum(MT%TLFP3(1:3)*MAP%DC(J,1:3))
  FP1X(J)=FP1X(J)+sum(MT%TLFP1(1:3)*DCX(J,1:3)) 
  FP3X(J)=FP3X(J)+sum(MT%TLFP3(1:3)*DCX(J,1:3)) 
  B2D(J)=B2D(J)+sum(MT%TLB2(1:3)*MAP%DC(J,1:3))/FLOAT(MT%LD2)
  B3D(J)=B3D(J)+sum(MT%TLB3(1:3)*MAP%DC(J,1:3))/FLOAT(MT%LD3)
  B4D(J)=B4D(J)+sum(MT%TLB4(1:3)*MAP%DC(J,1:3))/FLOAT(MT%LD4)
  UX(J)=UX(J)+sum(MT%TLU(1:3)*DCX(J,1:3)) 
  GD(J)=GD(J)+sum(MT%TLG(1:3)*MAP%DC(J,1:3))
 end do
!
 Z=SQRT(sum(FN**2))
 FN = FN/Z
 Z=SQRT(sum(BM(2:3)**2))
 BM(2:3)=BM(2:3)/Z 
 Z=sum(MT%TLBM**2)
 Z=SQRT(Z) 
 SBM1=MT%TLBM(1)/Z
 SBM2=MT%TLBM(2)/Z
 SBM3=MT%TLBM(3)/Z
 Z=sum(MT%TLFN**2)
 Z=SQRT(Z) 
 SFN1=MT%TLFN(1)/Z
 SFN2=MT%TLFN(2)/Z
 SFN3=MT%TLFN(3)/Z
!
 if (FN(3).le.0.0) then
  call Message%printError('HHComputeImages',' Line direction parallel to surface ')
 end if
 FNBM=SFN1*SBM1+SFN2*SBM2+SFN3*SBM3
 if (FNBM.le.0.0) then 
  call Message%printError('HHComputeImages',' Foil normal and beam not acute ')
 end if 

!***********************************************************************
!*       Computation of image size and positions of dislocations       *
!***********************************************************************
 Z=SQRT(FP(1)**2+FP(2)**2) 
 if(Z.eq.0.0) Z=1. 
 ! if(hhnl%SEP.ne.0.0) then
 !  call Message(' FAULT PLANE 2 ZERO WITH SEP NONZERO')
 ! end if
 AB(1)=FP(2)/Z 
 AB(2)=-FP(1)/Z
 AB(3)=0.0
 PT=hhnl%SEP*AB(1)
 SL=hhnl%SEP*AB(2)/BM(2)
 Z=SQRT(FP3(1)**2+FP3(2)**2) 
 if(Z.eq.0.0) Z=1. 
 AB1(1)=FP3(2)/Z 
 AB1(2)=-FP3(1)/Z
 AB1(3)=0.0
 PT2=hhnl%SEP2*AB1(1) 
 SL2=-hhnl%SEP2*AB1(2)/BM(2)
! we do not consider piezoelectric effects in this version of the program
!if (LPIEZO.eq.0) then 
  call ANCALC(MAP, MKAP, MAPN, MA, SCALE30) 
  if (MAPN%KRASH.eq.0) then
   if (SCALE30%LTEST.eq.1) then
     do I=1,3
      write(dataunit,*) (MA%AR(I,J),MA%AI(I,J),J=1,3)
     end do
     write(dataunit,*) (GD(I),I=1,3)
     do I=1,3
      write(dataunit,*) MA%PR(I),MA%PI(I)
     end do
     do I=1,3
      write(dataunit,*) (MA%EMR(I,J),MA%EMI(I,J),J=1,3)
     end do
     write(dataunit,"('  H  ',3F20.10)") ((MA%H(I,J),J=1,3),I=1,3)
   end if
   do JB=1,3 
    SUR(JB)=0.0
    SUI(JB)=0.0
    do K=1,3
     SUR(JB)=SUR(JB)+GD(K)*MA%AR(K,JB) 
     SUI(JB)=SUI(JB)+GD(K)*MA%AI(K,JB) 
    end do
   end do
   do JB=1,3 
    DR(JB)=SUR(JB)*MA%PR(JB)-SUI(JB)*MA%PI(JB) 
    DI(JB)=SUR(JB)*MA%PI(JB)+SUI(JB)*MA%PR(JB) 
   end do
   do JB=1,3 
    do L=1,3
     UR(JB,L)=0.0 
     UI(JB,L)=0.0 
     do J=1,3
      UR(JB,L)=UR(JB,L)+MA%EMR(JB,J)*MA%H(J,L) 
      UI(JB,L)=UI(JB,L)+MA%EMI(JB,J)*MA%H(J,L) 
     end do
    end do 
   end do 
   do JB=1,3 
    do L=1,3
     VR(JB,L)=DR(JB)*UR(JB,L)-DI(JB)*UI(JB,L) 
     VI(JB,L)=DR(JB)*UI(JB,L)+DI(JB)*UR(JB,L) 
    end do
   end do
   do JB=1,3 
    do L=1,3
     UI(JB,L)=VR(JB,L)*MA%PR(JB)+VI(JB,L)*MA%PI(JB) 
    end do
   end do
!***********************************************************
!*            First Dislocation                            *
!***********************************************************
   do J=1,3
    MRD%CN(J+8)=MA%PR(J)
    MRD%CN(J+12)=MA%PI(J)**2
    MRD%CN(J)=0.0
    MRD%CN(J+4)=0.0
    do L=1,3
     MRD%CN(J)=MRD%CN(J)+VR(J,L)*BD(L)
     MRD%CN(J+4)=MRD%CN(J+4)+UI(J,L)*BD(L)
    end do
   end do
   if (SCALE30%LTEST.eq.1) then
     write(dataunit,*) (MRD%CN(I),I=1,16)
   end if
!***********************************************************
!*        Second Dislocation                               *
!***********************************************************
   do J=1,3
    MRD%CN(J+21)=0.0 
    MRD%CN(J+25)=0.0 
    do L=1,3
     MRD%CN(J+21)=MRD%CN(J+21)+VR(J,L)*B2D(L) 
     MRD%CN(J+25)=MRD%CN(J+25)+UI(J,L)*B2D(L) 
    end do
   end do
!***********************************************************
!*       Third Dislocation                                 *
!***********************************************************
   do J=1,3
    MRD%CN(J+31)=0.0 
    MRD%CN(J+35)=0.0 
    do L=1,3
     MRD%CN(J+31)=MRD%CN(J+31)+VR(J,L)*B3D(L) 
     MRD%CN(J+35)=MRD%CN(J+35)+UI(J,L)*B3D(L) 
    end do
   end do
!***********************************************************
!*       Fourth Dislocation                                *
!***********************************************************
   do J=1,3
    MRD%CN(J+50)=0.0 
    MRD%CN(J+54)=0.0 
    do L=1,3
     MRD%CN(J+50)=MRD%CN(J+50)+VR(J,L)*B4D(L) 
     MRD%CN(J+54)=MRD%CN(J+54)+UI(J,L)*B4D(L) 
    end do
   end do 
   MRD%CN(4)=0.0
   MRD%CN(8)=0.0
   MRD%CN(12)=0.0 
   MRD%CN(16)=0.0 
   MRD%CN(25)=0.0 
   MRD%CN(29)=0.0 
   MRD%CN(35)=0.0 
   MRD%CN(39)=0.0 
   MRD%CN(54)=0.0 
   MRD%CN(58)=0.0 
   if (SCALE30%LTEST.eq.1) then
     write(dataunit,"(' H-ANCALC ',4F10.4)") ((MA%H(I,J),J=1,4),I=1,4) 
   end if
  else
   STOP 
  end if 
!else  
! the following lines are all commented out because we do not 
! want to consider piezoelectric effects at this point in time (1/25/2022)
!   call PANCALC(MAP, MKAP, MAPN, MA, MP, SCALE30)
!   if (MAPN%KRASH.eq.0) then
!    do K=1,4
!     SU(K)=CMPLX(0.0,0.0)
!     do I=1,3
!      SU(K)=SU(K)+GD(I)*MP%AS(I,K) 
!     end do
!    end do
!    do K=1,4
!     SU(K)=SU(K)*MP%PC(K) 
!    end do
! !*******************************************************************
! !   First Dislocation                                              *
! !*******************************************************************
!    do J=1,4 
!     MRD%CN(J+8)=REAL(MP%PC(J)) 
!     MRD%CN(J+12)=AIMAG(MP%PC(J))**2
!     CNX(J)=CMPLX(0.0,0.0) 
!     do L=1,4
!      CNX(J)=CNX(J)+(MP%EL(L,J)*BD(L)-MP%AS(L,J)*QL1(L))
!     end do
!     MRD%CN(J)=AIMAG(CNX(J)*SU(J))*2.0 
!     MRD%CN(J+4)=AIMAG(CNX(J)*SU(J)*CONJG(MP%PC(J)))*2.0
!    end do
! !*******************************************************************
! !   Second Dislocation                                             *
! !*******************************************************************
!    do J=1,4
!     CNX(J)=CMPLX(0.0,0.0) 
!     do L=1,4
!      CNX(J)=CNX(J)+(MP%EL(L,J)*B2D(L)-MP%AS(L,J)*QL2(L)) 
!     end do
!     MRD%CN(J+21)=AIMAG(CNX(J)*SU(J))*2.0
!     MRD%CN(J+25)=AIMAG(CNX(J)*SU(J)*CONJG(MP%PC(J)))*2.0 
!    end do
! !*******************************************************************
! !   Third Dislocation                                                *
! !****************************************************************** 
!    do J=1,4
!     CNX(J)=CMPLX(0.0,0.0) 
!     do L=1,4
!      CNX(J)=CNX(J)+(MP%EL(L,J)*B3D(L)-MP%AS(L,J)*QL3(L)) 
!     end do
!     MRD%CN(J+31)=AIMAG(CNX(J)*SU(J))*2.0
!     MRD%CN(J+35)=AIMAG(CNX(J)*SU(J)*CONJG(MP%PC(J)))*2.0 
!    end do
! !*******************************************************************
! !   Fourth Dislocation                                               *
! !*******************************************************************
!    do J=1,4
!     CNX(J)=CMPLX(0.0,0.0) 
!     do L=1,4
!      CNX(J)=CNX(J)+(MP%EL(L,J)*B4D(L)-MP%AS(L,J)*QL4(L)) 
!     end do
!     MRD%CN(J+50)=AIMAG(CNX(J)*SU(J))*2.0
!     MRD%CN(J+54)=AIMAG(CNX(J)*SU(J)*CONJG(MP%PC(J)))*2.0 
!    end do
! !*********************************************************************
! !*    BERECHNEN DER MATRIZEN   H, S, Q   UND AUSDRUCKEN              *
! !*       SIEHE DEFINITION NACH                                       *
! !*                - BARNETT UND LOTHE -                              *
! !*            PHYS.STAT.SOL(B) 67,105(1975)                          *
! !*********************************************************************
! !
!    do I=1,4
!     do J=1,4
!      MXXX(I,J)=CMPLX(0.0,0.0)
!      MYYY(I,J)=CMPLX(0.0,0.0)
!      MZZZ(I,J)=CMPLX(0.0,0.0)
!      do K=1,4
!       MXXX(I,J)=MXXX(I,J)+MP%EL(I,K)*MP%EL(J,K) 
!       MYYY(I,J)=MYYY(I,J)+MP%AS(I,K)*MP%EL(J,K) 
!       MZZZ(I,J)=MZZZ(I,J)+MP%AS(I,K)*MP%AS(J,K) 
!      end do
!      MA%H(I,J)=-(1.0/(2.0*cPi))*AIMAG(MXXX(I,J)) 
!      S(I,J)=-2.0*AIMAG(MYYY(I,J))
!      QS(I,J)=-2.0*AIMAG(MZZZ(I,J)) 
!     end do
!    end do
!    write(6,"(' ',4F10.4,' H- MATRIX')") ((MA%H(K,L),L=1,4),K=1,4) 
!    write(6,"(' ',4F10.4,' S- MATRIX')") ((S(K,L),L=1,4),K=1,4)
!    write(6,"(' ',4F10.4,' Q- MATRIX')") ((QS(K,L),L=1,4),K=1,4) 
!   else 
!    STOP
!   end if
!  end if 
 THBM=THICK/FNBM 
 KMIN=999
 KMAX=0
 KTOT=0

! 
!********************************************************************** 
!  BESTIMMUNG DER BILDLAENGE                                          * 
!********************************************************************** 
! 
 EXT1=ABS(SL+PT*FNX(1)/FNX(2)) 
 EXT2=ABS(SL2-PT2*FNX(1)/FNX(2)) 
 EXT3=ABS((SL2-PT2*FNX(1)/FNX(2))-(SL+PT*FNX(1)/FNX(2))) 
 if (SL.gt.0.0) then
  EXT4=((SL+PT*FNX(1)/FNX(2))+2*(SL2-PT2*FNX(1)/FNX(2))-(SL+PT*FNX(1)/FNX(2))) 
 else
  EXT4=EXT1+2*EXT2
 end if
 EXTRA=AMAX1(EXT1,EXT2,EXT3,EXT4)
 FRACTI=(FINISH-START)/THICK 
 DIVISO=BM(3)/BM(2)-FNX(3)/FNX(2)
 DELT=cPi*FRACTI*(THBM+EXTRA)/FLOAT(NNN)
 WL=((THICK*BM(2)/FN(3))+EXTRA/DIVISO)*FRACTI
 DELW=0.7650*cPi*WL/FLOAT(IROW-1)
 DELL=DELW/2.0+0.00000001
 MRD%CN(20)=cPi*PT/2.0
 MRD%CN(21)=cPi*SL/2.0
 MRD%CN(31)=MRD%CN(20)/BM(2) 
 MRD%CN(40)=cPi*PT2 
 MRD%CN(41)=cPi*SL2 
 MRD%CN(42)=MRD%CN(20)+MRD%CN(40)
 MRD%CN(43)=MRD%CN(21)-MRD%CN(41)
 MRD%CN(44)=MRD%CN(42)/BM(2) 

! 
!********************************************************************** 
!*       AUSDRUCKEN DER KONSTANTEN, KOORDINATEN UND DER GRAU-          *
!*       SKALA  (PARAMETER "LTEST"=1)                         * 
!********************************************************************** 
! 
 if (SCALE30%LTEST.eq.1) then
  write (dataunit,"(1H1,16H DC-KOORDINATEN )") 
  write (dataunit,"(1X,3F12.8)") ((MAP%DC(I,J),J=1,3),I=1,3) 
  write (dataunit,"(1H ,17H DCX-KOORDINATEN )") 
  write (dataunit,"(1X,3F12.8)") ((DCX(I,J),J=1,3),I=1,3)
  write (dataunit,"(1H ,48H KOORDINATEN DER BURGERSVEKTOREN UND FAULTPLANES)") 
  write (dataunit,"(1H0,3F12.8,'  BD',3F12.8,' B2D',3F12.8,' B3D'/1H ,3F12.8,'  BM',3F12.8,'  GD',3F12.8,'  FN'/1H ,3F12.8,' FNX', &
             3F12.8,'FP1X',3F12.8,' FPX'/1H ,3F12.8,'FP3X',3F12.8,'  FP',3F12.8,' FP3')") &
             BD,B2D,B3D,BM,GD,FN,FNX,FP1X,FPX,FP3X,FP,FP3
  write (dataunit,"(1H0,26H ERSTER SEPARATIONSVEKTOR )") 
  write (dataunit,"(1H ,3F12.8,3H AB,F12.8,3H SL,F12.8,3H PT,F12.8,5H DELT)") AB,SL,PT,DELT 
  write (dataunit,"(1H0,27H ZWEITER SEPARATIONSVEKTOR )") 
  write (dataunit,"(1H ,3F12.8,3HAB1,F12.8,3HPT2/1H ,F12.8,4HEXT1,1F12.8,4HEXT2,F12.8,4HEXT3,F12.8,5HEXTRA)") &
             AB1,SL2,PT2,EXT1,EXT2,EXT3,EXTRA
  write (dataunit,"(1H ,15H CN-KONSTANTEN )") 
  write (dataunit,"(1H ,4F12.8)") (MRD%CN(J),J=1,16)
  write (dataunit,"(' ',4F12.6)") (MRD%CN(J+21),J=1,8)
  write (dataunit,"(' ',4F12.6)") (MRD%CN(J+31),J=1,8)
  write (dataunit,"(' ',4F12.6)") (MRD%CN(J+50),J=1,8)
 end if
 MRD%CN(30)=1000.0 

 if (hhnl%wnum.eq.1) then 
   wvalues(1) = hhnl%wmin
 else  
   wstep = (hhnl%wmax-hhnl%wmin) / float(hhnl%wnum-1) 
   wvalues = hhnl%wmin + (/ (i*wstep, i=0,hhnl%wnum-1) /)
 end if

! loop over all the image pairs to be computed 
do imnum=1,hhnl%wnum
! set the excitation error parameter and related quantities
 MRD%CN(17) = wvalues(imnum)
 MRD%CN(18) = 2.0*MRD%CN(17) 
 MRD%X=0.0 
 MRD%Q=0.0 
 MRD%D=0.0
 MRD%DT=0.0
 MRD%YT=0.0
 MRD%ERROR=0.0001
 MRD%SKIP=0.0

 io_real(1) = MRD%CN(17)
 MRD%KOUNT = 0
 call Message%WriteValue(' starting image computation for w ', io_real,1)

! 
!********************************************************************** 
!*       BERECHNUNG DER UNTERGRUNDINTENSITAET DURCH INTEGRATION        *
!*       UEBER DIE GESAMMTE DICKE DES UNGESTOERTEN KRISTALLS           *
!********************************************************************** 
! 
 MRD%Y(1:8)=0.0 
 MRD%Y(1)=1.0
 MRD%Y(7)=1.0
 MRD%X1=THBM*cPi/FLOAT(ICOL/2) 
 CALL RKM(MRD)
 MRD%X1=cPi*THBM
 CALL RKM(MRD)
 BACK=1.0
 BACKD=1.0

!  This is a very long do-loop;  could be rewritten with function
!  and subroutine calls... and really should be parallelized using OpenMP
 do JC=1,IROW

  MRD%CN(19)=(FLOAT(JC)-FLOAT(IROW/2)-0.5)*DELW
  MOVE=0
  COORD(1)=-MRD%CN(42)
  COORD(2)=-MRD%CN(20)
  COORD(3)=+MRD%CN(20)
  COORD(4)=+MRD%CN(42) 
! 
!***************************************************************************
!*          ETWAIGE KORREKTUR DER VERSETZUNGSREIHENFOLGE 1-4 VON LINKS     * 
!*          NACH RECHTS. UMSPEICHERN DER KOORDINATEN                       *
!***************************************************************************
! 
  if ((COORD(2)-COORD(3)).ne.0.0) then
   STORE=COORD(3)
   COORD(3)=COORD(2) 
   COORD(2)=STORE
   if ((COORD(3)-COORD(4)).gt.0.0) then  
    STORE=COORD(4)
    COORD(4)=COORD(3) 
    COORD(3)=STORE
    STORE=COORD(2)
    COORD(2)=COORD(1) 
    COORD(1)=STORE
    if ((COORD(2)-COORD(3)).gt.0.0) then
     STORE=COORD(3)
     COORD(3)=COORD(2) 
     COORD(2)=STORE
    end if
   end if
  end if
!***************************************************************************
!*        DEFINITION DER SCHUTZZONE = 5 A                                  *
!***************************************************************************
  DELTA=5.0 
  DEL=DELTA*cPi/XIGEE
  DEL2=2.0*DEL
  HANDR=COORD-DEL 
  HANDL=COORD+DEL 
  if ((HANDL(1)-HANDR(2)).ge.0.0) then 
   if ((HANDL(2)-HANDR(3)).ge.0.0) then 
!*Case 1
    HANDL(1)=HANDL(4) 
    K=1 
    L=1 
!*Case 2
   else
    HANDL(1)=HANDL(2) 
    HANDL(3)=HANDL(4) 
    K=3 
    L=2 
   end if 
!*Case 3 
  else
   if ((HANDL(2)-HANDR(3)).ge.0.0) then 
    HANDL(2)=HANDL(3) 
    HANDL(3)=HANDL(4) 
    HANDR(3)=HANDR(4) 
    K=3 
    L=1 
!*Case 4
   else 
    K=4 
    L=1 
   end if
  end if
  do KK=1,K,L 
   DISTR=MRD%CN(19)+HANDL(KK)
   DISTRA=ABS(DISTR) 
   DISTL=-(HANDR(KK)+MRD%CN(19)) 
   DISTLA=ABS(DISTL) 
   if (DISTR.gt.0.0) then 
    if (DISTL.gt.0.0) then
     if ((DISTRA-DEL).gt.0.0) then
      if ((DISTLA-DEL).gt.0.0) then  
       if ((DISTRA-DEL2).gt.0.0) then 
        if ((DISTLA-DEL2).gt.0.0) then 
         call Message%printError('dohh_','DISLOCATIONS TOO CLOSE TOGETHER FOR COLUMN TO BE MOVED MEANINGFULLY ')
        else
         MRD%CN(15)=-HANDR(KK) 
         MOVE=1
         EXIT  ! the do KK loop
        end if
       else
         MRD%CN(15)=-HANDL(KK) 
         MOVE=1
         EXIT  ! the do KK loop
       end if
      else
       MRD%CN(15)=-HANDR(KK) 
       MOVE=1
       EXIT  ! the do KK loop
      end if
     else
       MRD%CN(15)=-HANDL(KK) 
       MOVE=1
       EXIT  ! the do KK loop
     end if
    end if ! else cycle KK loop
   end if ! else cycle KK loop
  end do 
! 
!***********************************************************************
!*      BEGRENZUNGEN DER STAPELFEHLER DURCH VERSETZUNGEN               *
!***********************************************************************
  XXX=MRD%CN(19)+MRD%CN(20) 
  YYY=MRD%CN(19)-MRD%CN(20) 
  ZZZ=MRD%CN(19)+MRD%CN(42) 
  VVV=MRD%CN(19)-MRD%CN(42) 
  MRD%CN(30)=MRD%CN(19)/BM(2) 
  COSA = 0.0
  SINA = 0.0
!***********************************
!*   First Stacking Fault          *
!***********************************
  FAULT1=10000.0
  if (YYY*VVV.le.0.0) then 
   if (sum(MT%TLS1**2).ne.0.0) then 
    ALPHA=cPi*sum(MT%TLG*MT%TLS1)*2.0/FLOAT(MT%LQ1) 
    COSA(1)=COS(ALPHA)
    SINA(1)=SIN(ALPHA)
    if (FP1X(2).ne.0.0) then 
     FAULT1=MRD%CN(21)-(MRD%CN(19)-MRD%CN(20))*FP1X(1)/FP1X(2)+hhnl%FAP1
    end if
   end if
  end if
!*************************************
!*  Second stacking fault          ***
!*************************************
  FAULT2=10001.0
  if (XXX*YYY.lt.0.0) then
   if (sum(MT%TLS2**2).ne.0.0) then
    ALPHA=cPi*sum(MT%TLG*MT%TLS2)*2.0/FLOAT(MT%LQ2) 
    COSA(2)=COS(ALPHA)
    SINA(2)=SIN(ALPHA)
    if (FPX(2).ne.0.0) then
     FAULT2=-MRD%CN(19)*FPX(1)/FPX(2)
    end if
   end if
  end if

!*************************************
!*  Third stacking fault           ***
!*************************************
  FAULT3=10002.0
  if (XXX.lt.0.0) then
   if (ZZZ.ge.0.0) then
    if (sum(MT%TLS3**2).ne.0.0) then 
     ALPHA=cPi*sum(MT%TLG*MT%TLS3)*2.0/FLOAT(MT%LQ3)
     COSA(3)=COS(ALPHA)
     SINA(3)=SIN(ALPHA)
     if (FP3X(2).ne.0.0) then
      FAULT3=-MRD%CN(21)-(MRD%CN(19)+MRD%CN(20))*FP3X(1)/FP3X(2)+hhnl%FAP3 
     end if
    end if 
   end if 
  end if
!
  STARTA=cPi*(EXTRA/2.0-(THBM+EXTRA)*FINISH/THICK)-(MRD%CN(19)*FNX(1)/FNX(2))
  SURFAC=STARTA+cPi*THBM 
  POSA = (/ FAULT1, FAULT2, FAULT3, SURFAC /)
  ITYPE = (/ 1, 2, 3, 4 /)
  do J=1,3
   LUCK=0
   do K=1,3
    if ((POSA(K)-POSA(K+1)).gt.0.0) then 
     STORE=POSA(K+1) 
     POSA(K+1)=POSA(K) 
     POSA(K)=STORE 
     ISTORE=ITYPE(K+1) 
     ITYPE(K+1)=ITYPE(K) 
     ITYPE(K)=ISTORE 
     LUCK=1
    end if
   end do
   if (LUCK.eq.0) EXIT 
  end do
  LSWITC=0
  do J=1,4
   if ((ITYPE(J)-4).ne.0) then 
    if (LSWITC.ne.0) then
     POSB(J)=POSA(J) 
     CYCLE
    else
     POSB(J)=-10050.0+FLOAT(J) 
    end if 
   else
    LSWITC=1
    POSB(J)=-10050.0+FLOAT(J) 
   end if
  end do
  IFLAG=0 
  MRD%X=STARTA
  MRD%X1=MRD%X+DELT 
  MRD%Y(1:8)=0.0 
  MRD%Y(1)=1.0
  MRD%Y(7)=1.0
  do JT=1,NNN 
   KOUNTF=1
   if ((JT-1).ne.0) then 
    MRD%X1=MRD%X1+DELT
   end if
 ifkount: do
    if ((KOUNTF-5).eq.0) EXIT ifkount
    if ((POSA(KOUNTF)-MRD%X).lt.0.0) then 
     KOUNTF=KOUNTF+1 
     CYCLE ifkount
    else
     if ((MRD%X1-POSA(KOUNTF)).ge.0.0) then 
      XX1=MRD%X1
      I=ITYPE(KOUNTF) 
      MRD%X1=POSA(KOUNTF) 
      call RKM(MRD)
      KTOT=KTOT+MRD%KOUNT 
      if (I.lt.4) then 
        TRAMP3=MRD%Y(3) 
        TRAMP7=MRD%Y(7) 
        MRD%Y(3)=MRD%Y(3)*COSA(I)-MRD%Y(4)*SINA(I)
        MRD%Y(7)=MRD%Y(7)*COSA(I)-MRD%Y(8)*SINA(I)
        MRD%Y(4)=MRD%Y(4)*COSA(I)+TRAMP3*SINA(I)
        MRD%Y(8)=MRD%Y(8)*COSA(I)+TRAMP7*SINA(I)
        MRD%X1=XX1
        POSA(KOUNTF)=-9000 - addon(I) 
        KOUNTF=KOUNTF+1 
      else 
        TEMPY(1:8)=MRD%Y(1:8) 
        MRD%X1=XX1
        POSA(KOUNTF)=-9000 - addon(4) 
        KOUNTF=KOUNTF+1 
        IFLAG=1
      end if 
     else 
      EXIT ifkount
     end if
    end if
   end do ifkount
!
   call RKM(MRD)
   KTOT=KTOT+MRD%KOUNT 
   DNR=MRD%Y(1)*MRD%Y(7)-MRD%Y(2)*MRD%Y(8)-MRD%Y(3)*MRD%Y(5)+MRD%Y(4)*MRD%Y(6) 
   DNI=MRD%Y(1)*MRD%Y(8)+MRD%Y(2)*MRD%Y(7)-MRD%Y(3)*MRD%Y(6)-MRD%Y(4)*MRD%Y(5) 
   DNN=1.0/(DNR**2+DNI**2) 
   FX(JT,1)=DNN*(MRD%Y(7)*DNR+MRD%Y(8)*DNI)
   FX(JT,2)=DNN*(MRD%Y(8)*DNR-MRD%Y(7)*DNI)
   FX(JT,3)=-DNN*(MRD%Y(3)*DNR+MRD%Y(4)*DNI) 
   FX(JT,4)=DNN*(MRD%Y(3)*DNI-MRD%Y(4)*DNR)
  end do
  if (IFLAG.eq.0) then 
    MRD%X1=SURFAC 
    KOUNTF=1
 ifkount2:  do
     if ((KOUNTF-5).eq.0) EXIT ifkount2
     if ((POSA(KOUNTF)-MRD%X).lt.0.0) then 
      KOUNTF=KOUNTF+1 
      CYCLE ifkount2
     else
      if ((MRD%X1-POSA(KOUNTF)).ge.0.0) then 
       XX1=MRD%X1
       I=ITYPE(KOUNTF) 
       MRD%X1=POSA(KOUNTF) 
       CALL RKM(MRD)
       KTOT=KTOT+MRD%KOUNT 
       if (I.eq.4) EXIT ifkount2
       TRAMP3=MRD%Y(3) 
       TRAMP7=MRD%Y(7) 
       MRD%Y(3)=MRD%Y(3)*COSA(I)-MRD%Y(4)*SINA(I)
       MRD%Y(7)=MRD%Y(7)*COSA(I)-MRD%Y(8)*SINA(I)
       MRD%Y(4)=MRD%Y(4)*COSA(I)+TRAMP3*SINA(I)
       MRD%Y(8)=MRD%Y(8)*COSA(I)+TRAMP7*SINA(I)
       MRD%X1=XX1
       POSA(KOUNTF)=-9050 - (I-1) 
       KOUNTF=KOUNTF+1 
       CYCLE ifkount2
      end if
     end if
    end do ifkount2
    if (I.ne.4) then
     CALL RKM(MRD)
     KTOT=KTOT+MRD%KOUNT 
    end if
    TEMPY(1:8)=MRD%Y(1:8) 
  end if
! 
  MRD%X=SURFAC
  MRD%X1=MRD%X+DELT 
  MRD%Y(1:8)=TEMPY(1:8) 
  do JM=1,NNN 
    KOUNTF=1
    if ((JM-1).ne.0) then 
     MRD%X1=MRD%X1+DELT
    end if
 ifkount3: do
     if ((KOUNTF-5).ne.0) then 
      if ((POSB(KOUNTF)-MRD%X).lt.0.0) then 
       KOUNTF=KOUNTF+1 
       CYCLE ifkount3
      else 
       if ((MRD%X1-POSB(KOUNTF)).lt.0.0) EXIT ifkount3
       XX1=MRD%X1
       I=ITYPE(KOUNTF) 
       MRD%X1=POSB(KOUNTF) 
       CALL RKM(MRD)
       KTOT=KTOT+MRD%KOUNT 
       if (I.eq.4) then
        write(6,"(1HG/26H1ITYPE(4) IN B INTEGRATION)")
        STOP
       end if
       TRAMP3=MRD%Y(3) 
       TRAMP7=MRD%Y(7) 
       MRD%Y(3)=MRD%Y(3)*COSA(I)-MRD%Y(4)*SINA(I)
       MRD%Y(7)=MRD%Y(7)*COSA(I)-MRD%Y(8)*SINA(I)
       MRD%Y(4)=MRD%Y(4)*COSA(I)+TRAMP3*SINA(I)
       MRD%Y(8)=MRD%Y(8)*COSA(I)+TRAMP7*SINA(I)
       MRD%X1=XX1
       POSA(KOUNTF)=-8050 - (I-1) 
       KOUNTF=KOUNTF+1 
       CYCLE ifkount3
      end if
     end if
    end do ifkount3
    call RKM(MRD)
    KTOT=KTOT+MRD%KOUNT 
    INDL=LLQ*JM+1 
! 
!********************************************************************** 
!***   BERECHNUNG DER  INTENSITAETEN UND GLEICHZEITIGES AB-         *** 
!***   SPEICHERN DER UNTERGRUNDINTENSITAETEB (TTD)                  *** 
!********************************************************************** 
! 
    TTB=(FX(JM,1)*MRD%Y(1)-FX(JM,2)*MRD%Y(2)+FX(JM,3)*MRD%Y(5)-FX(JM,4)*MRD%Y(6))**2+ &
        (FX(JM,1)*MRD%Y(2)+FX(JM,2)*MRD%Y(1)+FX(JM,3)*MRD%Y(6)+FX(JM,4)*MRD%Y(5))**2 
    TTD=(FX(JM,1)*MRD%Y(3)-FX(JM,2)*MRD%Y(4)+FX(JM,3)*MRD%Y(7)-FX(JM,4)*MRD%Y(8))**2+ &
        (FX(JM,1)*MRD%Y(4)+FX(JM,2)*MRD%Y(3)+FX(JM,3)*MRD%Y(8)+FX(JM,4)*MRD%Y(7))**2 
    TQB(INDL)=TTB
    TQD(INDL)=TTD
  end do  ! 
  TQB(1)=(TEMPY(1)**2+TEMPY(2)**2)
  TQD(1)=(TEMPY(3)**2+TEMPY(4)**2)

  if (LQ.eq.0) then
   do JZ=2,ICOL,2 
      TQB(JZ)=0.5*(TQB(JZ-1)+TQB(JZ+1))
      TQD(JZ)=0.5*(TQD(JZ-1)+TQD(JZ+1))
   end do
  end if
  ! do J=1,3
  !  if ((ABS(MRD%CN(19)+COORD(J))-DELL).le.0.0) then 
  !  end if
  ! end do
! 
!  AUSDRUCK EINER BILDZEILE 
! 
  do j=1,ICOL
   BFINTENS(j,JC,imnum)=TQB(j)
   DFINTENS(j,JC,imnum)=TQD(j)
  end do
 end do  ! from several pages back !!!


 WW=79.0*DELW/cPi 

 if (SCALE30%LTEST.eq.1) then  
! close the diagnostics file if this is the last image 
   if (imnum.eq.hhnl%wnum) CLOSE (UNIT=dataunit)
 end if
! 
!*******************************************************
!*       AUSDRUCK DER BILDLEGENDE                *******
!*******************************************************
!

! open a temporary file with the image legend for addition to the output HDF5 file 
   call Message%printMessage('  --> writing image legend to '//trim(legendfiles(imnum)) )
   OPEN(dataunit,FILE=trim(legendfiles(imnum)), status='UNKNOWN')

   write(dataunit,"('    HH.f90  ',' TWO-BEAM ',F6.2,' WL',F6.2,' WW', &
           F5.2,' STR ',F5.2,' FIN ',F7.3,'TH',F7.3,'THBM')") WL,WW,START,FINISH,THICK,THBM 
  
   write(dataunit,"(' ',3I2,'/',I1,'B   ',F8.4,' Q1 ',3I2,'U    ',3I2,  &
           'G    ',3I2,'BM   ',3I2,'FN',F7.3,'W',F9.3,'BACK')") &
           MT%LB,LD,QL1(4),MT%LU,MT%LG,MT%LBM,MT%LFN,MRD%CN(17),BACK
  
   write(dataunit,"(' ',3I2,'/',I1,'B2  ',F8.4,' Q2 ',F5.2,'SEP  ',3I2, &
           'FP1  ',3I2,'FP2  ',3I2,'FP3  ',3I2,'/',I1,'SH1  ',3I2,'/',I1, &
           'SH2  ',3I2,'/',I1,'SH3  ',' 4VS-SYM ',I2,' NNNN')") &
           MT%LB2,MT%LD2,QL2(4),hhnl%SEP,MT%LFP1,MT%LFP,MT%LFP3,MT%LS1,MT%LQ1,MT%LS2,MT%LQ2,MT%LS3,MT%LQ3,NNNN
  
   write(dataunit,"(' ',3I2,'/',I1,'B3  ',F8.4,' Q3 ',F5.2,'SEP2  ',F5.2, &
           ' FAP1 ',F5.2,' FAP3 ',F8.1,' XIGEE ',F4.1,' DELTA ',F8.5, &
            ' ANO ')") MT%LB3,MT%LD3,QL3(4),hhnl%SEP2,hhnl%FAP1,hhnl%FAP3,XIGEE,DELTA,MRD%ANO
  
   write(dataunit,"(' ',3I2,'/',I1,'B4  ',F8.4,' Q4 ',' PIEZO=',I1,' IND=',I1)") &
             MT%LB4,MT%LD4,QL4(4),LPIEZO,IND 
  
   write(dataunit,"('  G.B=',F10.4,'  G.(B X U)=',F10.4)") GINB,GINBXU

   CLOSE(dataunit, status='keep')
!
end do ! loop over all images 

call timer%makeTimeStamp()
dstr = timer%getDateString()
tstre = timer%getTimeString()

call Message%printMessage(' Writing HDF5 outfile file ',"(/A)")

call openFortranHDFInterface()
call self%writeHH4_HDFfile_(EMsoft, hhnl, BFINTENS, DFINTENS, dstr, tstrb, tstre, legendfiles, progname)
call closeFortranHDFInterface()

! cleanup: remove the legendfiles 
call Message%printMessage(' Cleaning up legend files (they have been added to the HDF5 output file)',"(/A/)")
do i=1,numim
   OPEN(dataunit,FILE=trim(legendfiles(i)), status='OLD')
   CLOSE(dataunit,status='delete')
end do

! image output requested ?  If so, then BF and DF pairs are put together in single images but with 
! a common gray scale for the entire series.
if (trim(hhnl%imageprefix).ne.'undefined') then 
  allocate(output_image(2*ICOL, IROW))
  mimi = minval ( (/ minval(BFINTENS), minval(DFINTENS) /) )
  mama = maxval ( (/ maxval(BFINTENS), maxval(DFINTENS) /) )
  BFINTENS = 255.0*(BFINTENS-mimi)/(mama-mimi)
  DFINTENS = 255.0*(DFINTENS-mimi)/(mama-mimi)

  do i=1,numim
   write (filenum,"(I3.3)") i
   fname = EMsoft%generateFilePath('EMdatapathname',trim(hhnl%imageprefix)//filenum//'.'//trim(hhnl%imagetype) )
   
   output_image(1:ICOL,1:IROW) = int(BFINTENS(1:ICOL,1:IROW,i))
   output_image(ICOL+1:2*ICOL,1:IROW) = int(DFINTENS(1:ICOL,1:IROW,i))

   im = image_t(output_image)
   if(im%empty()) call Message%printMessage(" HHComputeImages: failed to convert array to image")

   call im%write(trim(fname), iostat, iomsg) ! format automatically detected from extension
   if(0.ne.iostat) then
      call Message%printMessage(" HHComputeImages: failed to write image to file : "//iomsg)
   else  
      call Message%printMessage(' BF/DF images written to '//trim(fname))
   end if 
  end do 

end if

end associate 

end subroutine dohh_



end module mod_hh